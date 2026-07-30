//! This will convert data to a writable form and output the results. Needs some updates to account for the new data structures.
extern crate log;

use log::*;
use noodles::bgzf;
use std::collections::HashMap;
use std::io::Write;
use std::io::{self, BufReader, Lines, Read};
use std::num::ParseIntError;
use std::path::PathBuf;
use thiserror::Error;

use crate::file_tools::file_io::{
    create_output_file, is_gzipped_file, read_gzip_lines, read_lines,
};
use crate::structs::mutated_map::{AdCounter, MutatedMap};
use crate::structs::nucleotides::sequence_array_to_string;
use crate::structs::variants::{AlternateType, Provenance, Variant, VariantError};

#[derive(Debug, Error)]
pub enum VcfToolsError {
    #[error("Error reading the vcf file. Check file format: {0}")]
    MalformedVcf(String),
    #[error("IO error from the bed reader: {0}")]
    IoError(#[from] io::Error),
    #[error("Error parsing ints from the bed file: {0}")]
    ParserError(#[from] ParseIntError),
    #[error("File has an unknown or missing file extension: {0}")]
    FileExtensionUnknown(String),
    #[error("Error creating a variant from file: {0}")]
    VariantError(#[from] VariantError),
}

pub fn write_vcf(
    mutated_maps: &HashMap<String, Vec<MutatedMap>>,
    contig_order: &Vec<String>,
    fasta_lengths: &HashMap<String, usize>,
    reference_path: &PathBuf,
    overwrite_output: bool,
    output_vcf: &PathBuf,
    ad_counters: &HashMap<String, AdCounter>,
) -> io::Result<()> {
    // Takes:
    // mutated_maps: A map of contig names keyed to lists of mutated maps holding variants
    //     in that contig consisting of a tuple of (position, alt base, ref base).
    // fasta_order: A vector of contig names in the order of the reference fasta.
    // faste_lengths: lengths of contigs, required for vcf.
    // reference_path: The location of the reference file this vcf is showing variants from.
    // overwrite output: Whether to overwrite output for the vcf
    // output_vcf: The PathBuf object with the path to the output file
    // ad_counters: Per-contig allelic-depth counters from the gen-reads fragment loop.
    //     Drives FORMAT/AD, FORMAT/DP, FORMAT/AF. Empty (or missing per-contig) is fine —
    //     variants with no observed coverage emit `.` placeholders.
    // Result:
    // Throws and error if there's a problem, or else returns nothing.
    //
    // Output is wrapped in BGZF (block-gzip) — fully gzip-compatible for any
    // gzip-aware reader, AND tabix-indexable (`bcftools index -t`, which
    // plain-gzip output silently breaks). #206. Uses the noodles::bgzf writer
    // already pulled in by bam_writer.
    let outfile = create_output_file(output_vcf, overwrite_output)?;
    let mut buffer = bgzf::io::Writer::new(outfile);
    // add the vcf header
    writeln!(&mut buffer, "##fileformat=VCFv4.1")?;
    writeln!(&mut buffer, "##reference={}", reference_path.display())?;
    for contig in contig_order {
        let length = fasta_lengths[contig];
        writeln!(&mut buffer, "##contig=<ID={},length={}>", &contig, length)?;
    }
    // VCF spec §1.4.1: meta-information lines must be KEY=VALUE. #207.
    writeln!(
        &mut buffer,
        "##source=eidolon-{}",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(&mut buffer, "##ALT=<ID=DEL,Description=\"Deletion\">")?;
    writeln!(
        &mut buffer,
        "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">"
    )?;
    writeln!(&mut buffer, "##ALT=<ID=DUP,Description=\"Duplication\">")?;
    writeln!(&mut buffer, "##ALT=<ID=INV,Description=\"Inversion\">")?;
    writeln!(
        &mut buffer,
        "##ALT=<ID=CNV,Description=\"Copy number variant\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=EIDOLON_PROVENANCE,Number=1,Type=String,\
        Description=\"Origin of variant in this gen-reads run: \
        'denovo' = sampled by the simulator, 'input' = supplied via input_vcf:\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=EIDOLON_CCF,Number=1,Type=Float,\
        Description=\"Intended cellular fraction of this de-novo variant \
        (cancer-cell fraction of its assigned subclone; #405). The realized alt \
        fraction is EIDOLON_CCF x allele dosage; absent when no subclonal model is set.\">"
    )?;
    writeln!(
        &mut buffer,
        "##INFO=<ID=EIDOLON_VAF,Number=1,Type=Float,\
        Description=\"Intended observed variant allele fraction of this somatic \
        variant in the tumor/normal-mixed output (purity x dosage x CCF; #405). \
        Directly comparable to a caller's VAF on the merged reads — unlike FORMAT/AF, \
        which is measured per-pass (tumor-only). Absent outside cancer somatic runs.\">"
    )?;
    writeln!(
        &mut buffer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        &mut buffer,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,\
        Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
    )?;
    writeln!(
        &mut buffer,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
    )?;
    writeln!(
        &mut buffer,
        "##FORMAT=<ID=AF,Number=A,Type=Float,\
        Description=\"Allele frequency for each alt allele (AD[1] / DP)\">"
    )?;
    // Add a neat sample column
    writeln!(
        &mut buffer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEIDOLON_simulated_sample"
    )?;
    // write out mutations
    //
    // Per-contig: collect every record (literal + symbolic) from every block
    // into one Vec and sort by position before writing. The underlying
    // storage is a HashMap (variant_map), so iteration order is non-
    // deterministic and unsorted output breaks `bcftools index -t` (tabix
    // requires position-sorted input) — which broke the cancer_simulate.sh
    // merge step (see #185). Sorting in the writer is the simplest fix and
    // matches what every other VCF-producing tool does.
    //
    // Literals and SVs share a single sorted output stream because that's
    // how VCF spec expects them. Ties at the same position are stable-sorted
    // in collection order; rare in practice and tolerable.
    for contig in contig_order {
        let block_maps = &mutated_maps[contig];
        let contig_ad = ad_counters.get(contig);
        let mut records: Vec<&Variant> = Vec::new();
        for block_map in block_maps {
            records.extend(block_map.variant_map.values());
            records.extend(block_map.sv_records.iter());
        }
        records.sort_by_key(|v| v.location);
        // Defense in depth: if a variant somehow ends up in both
        // variant_map AND sv_records (bypassing MutatedMap::from_interval's
        // bin-routing), dedup by (location, ref, alt) before writing so we
        // don't emit it twice. See test_write_vcf_avoid_duplicating_symbolic_variants.
        records.dedup_by(|a, b| {
            a.location == b.location && a.reference == b.reference && a.alternate == b.alternate
        });

        for variant in records {
            let alt_str = match &variant.alternate {
                AlternateType::Literal(bases) => sequence_array_to_string(bases),
                AlternateType::Symbolic(sv) => sv.raw_alt.clone(),
            };
            let prov = variant.provenance.as_str();
            let is_symbolic = variant.alternate.is_symbolic();
            // INFO: symbolic records preserve any pre-existing INFO (SVLEN / END /
            // CN / SVTYPE) and append EIDOLON_PROVENANCE. Literal records have no
            // upstream INFO, so just emit EIDOLON_PROVENANCE alone.
            let info_str = if is_symbolic {
                match variant.info.as_deref() {
                    None | Some(".") | Some("") => format!("EIDOLON_PROVENANCE={prov}"),
                    Some(existing) => match strip_own_output_tokens(existing) {
                        Some(kept) => format!("{kept};EIDOLON_PROVENANCE={prov}"),
                        None => format!("EIDOLON_PROVENANCE={prov}"),
                    },
                }
            } else {
                // Literal records: EIDOLON_PROVENANCE is the base INFO. Simulator-emitted
                // tags (EIDOLON_CCF, EIDOLON_VAF) live in `info` on de-novo and reproductive
                // somatic (SomaticVcf) variants — pass those through. Germline input-VCF
                // literals deliberately drop their source INFO (only provenance emitted),
                // so existing golden output for input_vcf: runs is unchanged.
                match (variant.provenance, variant.info.as_deref()) {
                    (Provenance::Denovo | Provenance::SomaticVcf, Some(extra))
                        if !extra.is_empty() && extra != "." =>
                    {
                        // Denovo `info` is generated fresh (EIDOLON_CCF / EIDOLON_VAF),
                        // so in practice there is nothing to strip on this path — the
                        // SomaticVcf arm never carries source INFO either, since
                        // gen_reads::runner clears it on ingest. Kept as a cheap
                        // invariant so a future path that does populate `info` from a
                        // file can't reintroduce a duplicate provenance key.
                        match strip_own_output_tokens(extra) {
                            Some(kept) => format!("{kept};EIDOLON_PROVENANCE={prov}"),
                            None => format!("EIDOLON_PROVENANCE={prov}"),
                        }
                    }
                    _ => format!("EIDOLON_PROVENANCE={prov}"),
                }
            };
            // SAMPLE: for literal variants, AD/DP/AF come from the per-contig
            // counter populated by the gen-reads fragment loop; positions with
            // no observed coverage emit DP=0, AD=0,0, AF=`.`. Symbolic SVs use
            // span-based depth semantics that don't fit a point-based counter,
            // so they emit `.` placeholders to keep FORMAT consistent across
            // all records.
            let sample_str = if is_symbolic {
                format!("{}:.,.:.:.", variant.genotype_str)
            } else {
                let (ref_count, alt_count) = contig_ad
                    .and_then(|c| c.get(&variant.location).copied())
                    .unwrap_or((0, 0));
                let depth = ref_count + alt_count;
                let af_str = if depth > 0 {
                    format!("{:.4}", alt_count as f64 / depth as f64)
                } else {
                    ".".to_string()
                };
                format!(
                    "{}:{},{}:{}:{}",
                    variant.genotype_str, ref_count, alt_count, depth, af_str
                )
            };
            let line = format!(
                "{}\t{}\t.\t{}\t{}\t37\tPASS\t{}\tGT:AD:DP:AF\t{}",
                contig,
                variant.location + 1,
                sequence_array_to_string(&variant.reference),
                alt_str,
                info_str,
                sample_str,
            );
            writeln!(&mut buffer, "{}", line)?;
        }
    }
    Ok(())
}

/// Drop the eidolon-emitted routing tokens that this writer owns, from an INFO field
/// about to be passed through into new output.
///
/// `write_vcf` preserves a source record's INFO verbatim on symbolic (and de-novo
/// literal) records, so feeding a golden VCF back in via `input_vcf:` reproduces
/// whatever tokens it carried. Two classes must not survive that round trip:
///
/// * **`*_PROVENANCE`** — this writer re-emits it from `variant.provenance`. A stale
///   copy yields a duplicate key with *contradictory* values (`=denovo` from the source
///   next to the fresh `=input` describing how this pass actually used the record).
///   Worse, `vcf_merge::tumor_origin` tests for the exact field
///   `EIDOLON_PROVENANCE=denovo`, so a stale copy makes it label a variant present in
///   **both** passes as `somatic` — putting a germline carry-through into the somatic
///   ground truth that every bundled benchmark scores against.
/// * **`*_ORIGIN`** — added *after* this writer by `vcf_merge`, which is also the only
///   place that declares it. Passing one through emits an undeclared INFO tag: fine when
///   streaming VCF→VCF, but a hard failure on BCF translation, which
///   `bcftools concat | sort` does internally.
///
/// Both the current `EIDOLON_*` spellings and the pre-v3.0.0 `NEAT_*` ones are stripped.
/// Scoping this to only the legacy names would fix the shrinking half of the problem
/// (v2.0.0 goldens) and leave the growing half (v3.0.0 goldens, which is what users
/// produce from now on).
///
/// **`EIDOLON_CCF` / `EIDOLON_VAF` are deliberately kept.** Unlike provenance, they are
/// #405's ground truth and live *in* `variant.info` — `gen_reads::runner` appends them
/// there before the write — so stripping them would destroy the tags on every de-novo
/// somatic record. They are declared by this writer, so passing one through is harmless.
/// `NEAT_CCF`/`NEAT_VAF`/`NEAT_REASON` are still dropped: those spellings are never
/// re-emitted, so a survivor would be an undeclared tag.
///
/// Returns `None` when nothing is left to keep, so callers can emit provenance alone
/// instead of a leading `;`. Only these exact keys are removed — third-party INFO
/// (`DP`, `AF`, `SVTYPE`, …) is untouched.
fn strip_own_output_tokens(info: &str) -> Option<String> {
    const OWNED_KEYS: [&str; 7] = [
        // Current spellings: re-emitted here, or declared only downstream.
        "EIDOLON_PROVENANCE",
        "EIDOLON_ORIGIN",
        // Pre-v3.0.0 spellings: never re-emitted, so any survivor is undeclared.
        "NEAT_PROVENANCE",
        "NEAT_ORIGIN",
        "NEAT_CCF",
        "NEAT_VAF",
        "NEAT_REASON",
    ];
    let kept: Vec<&str> = info
        .split(';')
        .filter(|field| {
            let key = field.split('=').next().unwrap_or(field);
            !OWNED_KEYS.contains(&key)
        })
        .filter(|field| !field.is_empty())
        .collect();
    if kept.is_empty() {
        None
    } else {
        Some(kept.join(";"))
    }
}

pub fn read_vcf(vcf_file: PathBuf) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    if is_gzipped_file(&vcf_file)? {
        process_gzip_vcf(&vcf_file)
    } else {
        process_vcf(&vcf_file)
    }
}

/// Lean counterpart to [`read_vcf`] for callers that only need the structural
/// fields of each [`Variant`] (variant_type, location, reference, alternate,
/// genotype). Skips populating `info`, `id`, `filter`, `genotype_str`,
/// `quality_score`, and the FORMAT/SAMPLE vectors — none of which are read
/// by gen_mut_model. On gnomAD-scale inputs this cuts peak RSS by several
/// GB (the INFO column alone is hundreds of bytes to several KB per record).
///
/// Use [`read_vcf`] when the caller needs to round-trip records back to VCF
/// (e.g. gen_reads) or read `v.filter` (compare_vcfs).
pub fn read_vcf_lean(vcf_file: PathBuf) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    if is_gzipped_file(&vcf_file)? {
        let reader = read_gzip_lines(&vcf_file)?;
        read_open_vcf_lean(reader)
    } else {
        let reader = read_lines(&vcf_file)?;
        read_open_vcf_lean(reader)
    }
}

fn process_gzip_vcf(filename: &PathBuf) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    let reader = read_gzip_lines(filename)?;
    read_open_vcf(reader)
}

fn process_vcf(filename: &PathBuf) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    let reader = read_lines(filename)?;
    read_open_vcf(reader)
}

fn read_open_vcf<P: Read>(
    reader: Lines<BufReader<P>>,
) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    let mut file_records: HashMap<String, Vec<Variant>> = HashMap::new();
    for result in reader {
        match result {
            Ok(line) => {
                if line.starts_with("#") {
                    continue;
                } else {
                    let split_line: Vec<&str> = line.split("\t").collect();
                    let chrom = split_line[0];
                    let pos = split_line[1].parse();
                    let location: usize = {
                        match pos {
                            Ok(num) => num,
                            Err(error) => return Err(VcfToolsError::ParserError(error)),
                        }
                    };
                    let id = split_line[2];
                    let vcf_reference = split_line[3];
                    let vcf_alternate = split_line[4];
                    if vcf_alternate.contains(',') {
                        debug!("Skipping multi-allelic site at {}:{}", chrom, location);
                        continue;
                    }
                    // QUAL may be "." in real-world VCFs
                    let quality_score: usize = split_line[5].parse().unwrap_or(0);
                    let filter = split_line[6];
                    let info_field = split_line[7];
                    let mut format: Vec<String> = Vec::new();
                    let mut sample: Vec<String> = Vec::new();
                    if split_line.len() > 9 {
                        let fmt = split_line[8];
                        let smp = split_line[9];
                        // "." means no FORMAT/SAMPLE data — handled by the
                        // sites-only fallback below.
                        if fmt != "." && smp != "." {
                            let format_split: Vec<&str> = fmt.split(":").collect();
                            let sample_split: Vec<&str> = smp.split(":").collect();
                            let format_len = format_split.len();
                            // The VCF spec allows trailing per-sample FORMAT fields to be
                            // dropped (all but GT), so a sample may carry FEWER values than
                            // FORMAT declares — pad the missing trailing values as "." rather
                            // than rejecting spec-compliant VCFs (e.g. the GIAB truth sets:
                            // FORMAT=GT:PS with an unphased sample of just "1/1"). Only MORE
                            // sample values than FORMAT keys is genuinely malformed.
                            if sample_split.len() > format_len {
                                return Err(VcfToolsError::MalformedVcf(
                                    "sample has more fields than FORMAT declares, invalid VCF"
                                        .to_string(),
                                ));
                            }
                            for i in 0..format_len {
                                format.push(format_split[i].to_string());
                                sample
                                    .push(sample_split.get(i).copied().unwrap_or(".").to_string());
                            }
                        }
                    }
                    if split_line.len() > 10 {
                        warn!("Currently eidolon can only read one sample per record")
                    }
                    // Sites-only VCFs (gnomAD-SV, dbSNP, etc.) have no
                    // FORMAT / SAMPLE columns at all. Without a genotype
                    // we can't tell hom from het, but rejecting the record
                    // outright would drop all of gnomAD-SV. Default to
                    // heterozygous ("0/1") — the trainer will then report
                    // `homozygous_frequency = 0` for the corpus, which is
                    // an honest "we don't know" rather than a synthesized
                    // 50/50 split.
                    if !format.iter().any(|f| f == "GT") {
                        format.push("GT".to_string());
                        sample.push("0/1".to_string());
                    }
                    let variant = match Variant::from_file(
                        location,
                        id,
                        filter,
                        info_field,
                        vcf_reference,
                        vcf_alternate,
                        quality_score,
                        format,
                        sample,
                    ) {
                        Ok(v) => v,
                        Err(e) => {
                            warn!("Skipping variant at {}:{} — {}", chrom, location, e);
                            continue;
                        }
                    };
                    // Standard hashmap safety check
                    if !file_records.contains_key(chrom) {
                        file_records.insert(chrom.to_string(), Vec::new());
                    }
                    file_records.get_mut(chrom).unwrap().push(variant);
                }
            }
            Err(error) => return Err(VcfToolsError::IoError(error)),
        }
    }
    Ok(file_records)
}

/// Lean parser used by [`read_vcf_lean`]. Mirrors [`read_open_vcf`] but
/// extracts just the GT token (without allocating a `Vec<String>` for
/// FORMAT/SAMPLE) and hands it to [`Variant::from_file_lean`], which leaves
/// the unused string fields empty. The cumulative effect is several GB
/// less peak RSS on a gnomAD-scale ingest.
fn read_open_vcf_lean<P: Read>(
    reader: Lines<BufReader<P>>,
) -> Result<HashMap<String, Vec<Variant>>, VcfToolsError> {
    let mut file_records: HashMap<String, Vec<Variant>> = HashMap::new();
    for result in reader {
        match result {
            Ok(line) => {
                if line.starts_with("#") {
                    continue;
                }
                let split_line: Vec<&str> = line.split("\t").collect();
                let chrom = split_line[0];
                let pos = split_line[1].parse();
                let location: usize = match pos {
                    Ok(num) => num,
                    Err(error) => return Err(VcfToolsError::ParserError(error)),
                };
                let vcf_reference = split_line[3];
                let vcf_alternate = split_line[4];
                if vcf_alternate.contains(',') {
                    debug!("Skipping multi-allelic site at {}:{}", chrom, location);
                    continue;
                }
                let info_field = split_line[7];
                // Pull just the GT field if present. Sites-only VCFs (no
                // FORMAT/SAMPLE) and records whose FORMAT lacks GT default
                // to "0/1" — matching the full-fat reader's behavior.
                let gt_str: &str = if split_line.len() > 9 {
                    match extract_gt_str(split_line[8], split_line[9])? {
                        Some(gt) => gt,
                        None => "0/1",
                    }
                } else {
                    "0/1"
                };
                if split_line.len() > 10 {
                    warn!("Currently eidolon can only read one sample per record")
                }
                let variant = match Variant::from_file_lean(
                    location,
                    info_field,
                    vcf_reference,
                    vcf_alternate,
                    gt_str,
                ) {
                    Ok(v) => v,
                    Err(e) => {
                        warn!("Skipping variant at {}:{} — {}", chrom, location, e);
                        continue;
                    }
                };
                file_records
                    .entry(chrom.to_string())
                    .or_default()
                    .push(variant);
            }
            Err(error) => return Err(VcfToolsError::IoError(error)),
        }
    }
    Ok(file_records)
}

/// Find the GT column in a `FORMAT`/`SAMPLE` pair without allocating a
/// `Vec<String>` for either. Returns the GT token (e.g. `"0/1"`) by
/// reference into the input. Mirrors the validation in `read_open_vcf`:
///   - `fmt == "."` or `smp == "."`  → `Ok(None)` (sites-only fallback)
///   - sample longer than FORMAT     → `MalformedVcf` error (trailing-dropped is allowed)
///   - no GT key                     → `Ok(None)` (defaults applied by caller)
fn extract_gt_str<'a>(fmt: &'a str, smp: &'a str) -> Result<Option<&'a str>, VcfToolsError> {
    if fmt == "." || smp == "." {
        return Ok(None);
    }
    let format_parts: Vec<&str> = fmt.split(':').collect();
    let sample_parts: Vec<&str> = smp.split(':').collect();
    // Trailing FORMAT fields may be dropped per-sample (all but GT), so the sample
    // can be shorter than FORMAT; only more sample values than keys is malformed.
    if sample_parts.len() > format_parts.len() {
        return Err(VcfToolsError::MalformedVcf(
            "sample has more fields than FORMAT declares, invalid VCF".to_string(),
        ));
    }
    for (i, key) in format_parts.iter().enumerate() {
        if *key == "GT" {
            // GT is not droppable, but guard against a truncated sample anyway.
            return Ok(sample_parts.get(i).copied());
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::{
        nucleotides::Nucleotide,
        variants::{Genotype, Provenance, SvData, SvType, Variant, VariantType},
    };
    use std::io::Write;

    #[test]
    pub fn test_write_vcf() {
        let temp_dir = tempfile::tempdir().unwrap();
        let variant1 = Variant::new(
            VariantType::SNP,
            3,
            &vec![Nucleotide::A],
            &vec![Nucleotide::G],
            &mut vec![0, 1],
        )
        .unwrap();
        let variant2 = Variant::new(
            VariantType::SNP,
            7,
            &vec![Nucleotide::T],
            &vec![Nucleotide::G],
            &mut vec![1, 1],
        )
        .unwrap();
        let mutated_map = MutatedMap::from_interval(100, 200, vec![variant1, variant2]).unwrap();
        let mut mutated_maps = HashMap::new();
        mutated_maps.insert("chr1".to_string(), vec![mutated_map]);
        let fasta_order = Vec::from(["chr1".to_string()]);
        let fasta_length = HashMap::from([("chr1".to_string(), 20)]);
        let reference_path = PathBuf::from("test_data/references/H1N1.fa");
        let output_filename = temp_dir.path().join("good_test.vcf");
        let result = write_vcf(
            &mutated_maps,
            &fasta_order,
            &fasta_length,
            &reference_path,
            false,
            &output_filename,
            &HashMap::new(),
        );
        result.unwrap();
        assert!(output_filename.exists());

        // Read back and verify content
        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output_filename)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        // Header lines must be present
        assert!(lines.iter().any(|l| l == "##fileformat=VCFv4.1"));
        assert!(
            lines
                .iter()
                .any(|l| l.contains("##contig=<ID=chr1,length=20>"))
        );
        // ##source= is the spec-blessed origin tag (was the malformed
        // "##Generated by eidolon" line in v1.10.x — see #207).
        assert!(lines.iter().any(|l| l.starts_with("##source=eidolon-")));
        assert!(
            lines
                .iter()
                .any(|l| l.contains("#CHROM\tPOS\tID\tREF\tALT"))
        );

        // Both variant lines must be present (HashMap iteration order is non-deterministic).
        // Variants from `Variant::new` carry `Provenance::Denovo` (→ INFO/EIDOLON_PROVENANCE=denovo)
        // and an empty ad_counter means no reads were simulated (→ AD=0,0 DP=0 AF=`.`).
        assert!(
            lines.iter().any(|l| l
                == "chr1\t4\t.\tA\tG\t37\tPASS\tEIDOLON_PROVENANCE=denovo\tGT:AD:DP:AF\t0/1:0,0:0:."),
            "missing het SNP line; got: {:?}",
            lines
        );
        assert!(
            lines.iter().any(|l| l
                == "chr1\t8\t.\tT\tG\t37\tPASS\tEIDOLON_PROVENANCE=denovo\tGT:AD:DP:AF\t1/1:0,0:0:."),
            "missing hom SNP line; got: {:?}",
            lines
        );
    }

    #[test]
    pub fn test_write_vcf_multi_contig() {
        let temp_dir = tempfile::tempdir().unwrap();

        let map1 = MutatedMap::from_interval(
            100,
            200,
            vec![
                Variant::new(
                    VariantType::SNP,
                    5,
                    &vec![Nucleotide::A],
                    &vec![Nucleotide::C],
                    &mut vec![0, 1],
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let map2 = MutatedMap::from_interval(
            100,
            200,
            vec![
                Variant::new(
                    VariantType::SNP,
                    10,
                    &vec![Nucleotide::G],
                    &vec![Nucleotide::T],
                    &mut vec![1, 1],
                )
                .unwrap(),
            ],
        )
        .unwrap();

        let mut mutated_maps = HashMap::new();
        mutated_maps.insert("chr1".to_string(), vec![map1]);
        mutated_maps.insert("chr2".to_string(), vec![map2]);

        let contig_order = vec!["chr1".to_string(), "chr2".to_string()];
        let fasta_lengths = HashMap::from([
            ("chr1".to_string(), 300usize),
            ("chr2".to_string(), 500usize),
        ]);

        let output_path = temp_dir.path().join("multi_contig.vcf");
        write_vcf(
            &mutated_maps,
            &contig_order,
            &fasta_lengths,
            &PathBuf::from("ref.fa"),
            false,
            &output_path,
            &HashMap::new(),
        )
        .unwrap();

        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        // Both contigs must appear in the header
        assert!(
            lines
                .iter()
                .any(|l| l.contains("##contig=<ID=chr1,length=300>"))
        );
        assert!(
            lines
                .iter()
                .any(|l| l.contains("##contig=<ID=chr2,length=500>"))
        );

        // Variant lines for both contigs must be present
        assert!(
            lines.iter().any(|l| l.starts_with("chr1\t6\t")),
            "missing chr1 variant; got: {:?}",
            lines
        );
        assert!(
            lines.iter().any(|l| l.starts_with("chr2\t11\t")),
            "missing chr2 variant; got: {:?}",
            lines
        );
    }

    #[test]
    fn test_read_vcf_lean_drops_unused_string_fields() {
        // Lean reader produces structurally-correct Variants but leaves
        // the per-record string fields empty/None. Locks in the memory-
        // saving contract: a future refactor that re-populates info /
        // filter / etc. here would silently re-introduce gnomAD-scale
        // bloat in gen_mut_model.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\trsBIG\tA\tT\t60\tPASS\tAF=0.01;AC=2;AN=100;long_info_string_that_would_bloat_memory\tGT:DP:GQ\t0/1:30:99\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf_lean(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").expect("chr1 not parsed");
        assert_eq!(variants.len(), 1);
        let v = &variants[0];
        // Structural fields still populated.
        assert_eq!(v.location, 100);
        assert_eq!(v.variant_type, crate::structs::variants::VariantType::SNP);
        assert!(matches!(
            v.genotype,
            crate::structs::variants::Genotype::Heterozygous
        ));
        // String fields all empty / None.
        assert!(v.id.is_none(), "id should be None");
        assert!(v.filter.is_none(), "filter should be None");
        assert!(v.info.is_none(), "info should be None (the big one)");
        assert!(v.quality_score.is_none(), "quality_score should be None");
        assert!(v.genotype_str.is_empty(), "genotype_str should be empty");
        assert!(v.format.is_empty(), "format should be empty");
        assert!(v.sample.is_empty(), "sample should be empty");
    }

    #[test]
    fn test_read_vcf_lean_sites_only_defaults_to_heterozygous() {
        // Sites-only VCFs (8 columns, no FORMAT/SAMPLE) must still produce
        // a Variant via the lean path — matching read_vcf's contract.
        let vcf_content = "\
##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=500\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf_lean(tmp.path().to_path_buf()).unwrap();
        let v = &records.get("chr1").expect("chr1 missing")[0];
        // SV metadata still survives via SvData (parsed at construction
        // time from the INFO column, then info: Option<String> dropped).
        let sv = v.alternate.as_symbolic().expect("symbolic ALT");
        assert_eq!(sv.end, Some(500));
        assert_eq!(sv.sv_type, SvType::Del);
        assert!(matches!(
            v.genotype,
            crate::structs::variants::Genotype::Heterozygous
        ));
    }

    #[test]
    fn test_read_vcf_qual_dot() {
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let result = read_vcf(tmp.path().to_path_buf());
        assert!(result.is_ok(), "Expected Ok, got {:?}", result);
        let records = result.unwrap();
        let variants = records.get("chr1").expect("chr1 not found");
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].quality_score, Some(0));
    }

    #[test]
    fn test_read_vcf_trailing_format_fields_dropped() {
        // Per the VCF spec, trailing per-sample FORMAT fields may be dropped
        // (all but GT). GIAB truth VCFs do this — e.g. FORMAT=GT:PS with an
        // unphased sample carrying only the GT. The parser must accept it and
        // still recover the genotype, not reject it as malformed (job 19899126:
        // gen-mut-model on the HG002 v4.2.1 truth VCF died here).
        let vcf_content = "\
##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:PS\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf())
            .expect("trailing-dropped FORMAT must parse, not error");
        let variants = records.get("chr1").expect("chr1 not found");
        assert_eq!(variants.len(), 1);
        assert!(matches!(variants[0].genotype, Genotype::Heterozygous));
    }

    #[test]
    fn test_read_vcf_lean_trailing_format_fields_dropped() {
        // Same spec allowance via the lean path (extract_gt_str).
        let vcf_content = "\
##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:PS:DP\t1/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf_lean(tmp.path().to_path_buf())
            .expect("trailing-dropped FORMAT must parse via the lean path");
        let v = &records.get("chr1").expect("chr1 missing")[0];
        assert!(matches!(v.genotype, Genotype::Homozygous));
    }

    #[test]
    fn test_read_vcf_sample_longer_than_format_rejected() {
        // The inverse is still malformed: MORE sample values than FORMAT keys.
        let vcf_content = "\
##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1:99\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let result = read_vcf(tmp.path().to_path_buf());
        assert!(
            matches!(result, Err(VcfToolsError::MalformedVcf(_))),
            "sample longer than FORMAT should be rejected, got {:?}",
            result
        );
    }

    #[test]
    fn test_read_vcf_multi_allelic_alt_skipped() {
        // Multi-allelic records (ALT field containing ',') are explicitly skipped
        // by the parser. This test locks in that contract.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT,C\t60\tPASS\t.\tGT\t1|2\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").expect("chr1 not found");
        assert_eq!(
            variants.len(),
            1,
            "multi-allelic record should be skipped; only the simple SNP at 200 should remain",
        );
        assert_eq!(variants[0].location, 200);
    }

    #[test]
    fn test_read_vcf_symbolic_del_parses_to_symbolic_with_info() {
        // A `<DEL>` ALT becomes AlternateType::Symbolic carrying the parsed
        // INFO fields (SVTYPE / END), with the raw ALT preserved verbatim so
        // the writer can round-trip it.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=200\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").expect("chr1 not parsed");
        assert_eq!(variants.len(), 1);
        let sv = variants[0]
            .alternate
            .as_symbolic()
            .expect("symbolic ALT should classify as Symbolic");
        assert_eq!(sv.raw_alt, "<DEL>");
        assert_eq!(sv.sv_type, SvType::Del);
        assert_eq!(sv.end, Some(200));
    }

    #[test]
    fn test_read_vcf_symbolic_dup_with_cn_and_svlen() {
        // <DUP> with both SVLEN and CN should populate the corresponding
        // SvData fields. svlen is signed and stored as-is.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t1000\t.\tT\t<DUP:TANDEM>\t60\tPASS\tSVTYPE=DUP;SVLEN=500;CN=3\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let sv = records.get("chr1").unwrap()[0]
            .alternate
            .as_symbolic()
            .unwrap();
        assert_eq!(sv.sv_type, SvType::Dup);
        assert_eq!(sv.raw_alt, "<DUP:TANDEM>");
        assert_eq!(sv.svlen, Some(500));
        assert_eq!(sv.copy_number, Some(3));
        // END wasn't supplied, so it stays None.
        assert!(sv.end.is_none());
    }

    #[test]
    fn test_read_vcf_symbolic_cnv() {
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t500\t.\tA\t<CNV>\t60\tPASS\tSVTYPE=CNV;END=2500;CN=4\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let v = &records.get("chr1").unwrap()[0];
        assert_eq!(v.variant_type, VariantType::Complex);
        let sv = v.alternate.as_symbolic().unwrap();
        assert_eq!(sv.sv_type, SvType::Cnv);
        assert_eq!(sv.end, Some(2500));
        assert_eq!(sv.copy_number, Some(4));
        // span() = END - POS + 1 = 2500 - 500 + 1 = 2001
        assert_eq!(sv.span(500), Some(2001));
    }

    #[test]
    fn test_read_vcf_symbolic_round_trip_through_writer() {
        // End-to-end: a symbolic ALT read from input survives through the
        // writer back to a verbatim ALT string — read + write together,
        // driving the real reader path (the companion writer-only test below
        // stubs in a synthetic Variant).
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t150\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=300\tGT\t1/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let v = records.get("chr1").unwrap()[0].clone();

        let temp_dir = tempfile::tempdir().unwrap();
        let map = MutatedMap::from_interval(0, 1000, vec![v]).unwrap();
        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);
        let output = temp_dir.path().join("rt.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 1000usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();
        let body: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        // The ALT column survives verbatim through the read+write cycle.
        assert!(
            body.iter().any(|l| l.contains("\tG\t<DEL>\t")),
            "expected <DEL> ALT to round-trip; got: {:?}",
            body
        );
    }

    #[test]
    fn test_read_vcf_symbolic_inv() {
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t500\t.\tA\t<INV>\t60\tPASS\tSVTYPE=INV;END=2500\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let v = &records.get("chr1").unwrap()[0];
        assert_eq!(v.variant_type, VariantType::Complex);
        let sv = v.alternate.as_symbolic().unwrap();
        assert_eq!(sv.sv_type, SvType::Inv);
        assert_eq!(sv.end, Some(2500));
        // span() = END - POS + 1 = 2500 - 500 + 1 = 2001
        assert_eq!(sv.span(500), Some(2001));
    }

    #[test]
    fn test_read_vcf_symbolic_svtype_disagreement_still_parses() {
        // If SVTYPE conflicts with the ALT tag, we trust the ALT and warn.
        // The variant is still produced (we don't drop the record).
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\t<DEL>\t60\tPASS\tSVTYPE=DUP;END=200\tGT\t0/1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let sv = records.get("chr1").unwrap()[0]
            .alternate
            .as_symbolic()
            .unwrap();
        // ALT wins.
        assert_eq!(sv.sv_type, SvType::Del);
        assert_eq!(sv.raw_alt, "<DEL>");
    }

    #[test]
    fn test_read_vcf_info_with_embedded_semicolons_preserved() {
        // INFO is split only on \t, so embedded ';' in the INFO column should not corrupt
        // downstream fields (FORMAT and SAMPLE). Lock in the round-trip.
        let info = "DESC=\"contains;semis\";END=300;NS=1";
        let vcf_content = format!(
            "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t{}\tGT\t0/1\n",
            info,
        );
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").unwrap();
        assert_eq!(variants.len(), 1);
        // The INFO column survives intact; GT (in SAMPLE column) is unaffected.
        assert_eq!(variants[0].info.as_deref(), Some(info));
    }

    #[test]
    fn test_read_vcf_gzipped_input() {
        // A plain-gzip-encoded VCF (.vcf.gz) — not bgzf — should still be
        // parsed identically. This tests reader compatibility with external
        // VCFs that were gzipped (not bgzipped) before being handed to eidolon;
        // since v1.11.1 our own writer emits bgzf, but the reader must keep
        // accepting plain gzip too.
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n\
chr2\t250\t.\tG\tC\t60\tPASS\t.\tGT\t1|1\n";

        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("input.vcf.gz");
        let f = std::fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(f, Compression::default());
        encoder.write_all(vcf_content.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let records = read_vcf(path).unwrap();
        assert_eq!(records.get("chr1").unwrap().len(), 1);
        assert_eq!(records.get("chr2").unwrap().len(), 1);
    }

    /// `MutatedMap::from_interval` routes symbolic ALTs to `sv_records`
    /// and literal ALTs to `variant_map` — a variant is never in both bins
    /// via the public API. This test manually constructs the malformed
    /// state (variant in both bins) to confirm the writer is robust to it
    /// even when someone bypasses `from_interval`. write_vcf dedups by
    /// `variant.location` before emission as defense-in-depth.
    #[test]
    fn test_write_vcf_avoid_duplicating_symbolic_variants() {
        use crate::structs::variants::{AlternateType, Genotype, Provenance, SvData, SvType};
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 100,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };

        // Manually construct a MutatedMap where the same SV is in both variant_map AND sv_records.
        // This simulates the problematic state described in the issue.
        let mut variant_map = HashMap::new();
        variant_map.insert(100usize, sv.clone());
        let map = MutatedMap {
            flagged_positions: vec![100],
            block_interval: (0, 1000),
            variant_map,
            sv_records: vec![sv],
        };

        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);

        let temp_dir = tempfile::tempdir().unwrap();
        let output = temp_dir.path().join("dupe.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 1000usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();

        let body: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        let variant_lines: Vec<_> = body.iter().filter(|l| l.contains("\t101\t")).collect();
        assert_eq!(
            variant_lines.len(),
            1,
            "Expected only 1 entry for the variant at position 101, but found: {:?}",
            variant_lines
        );
    }

    #[test]
    fn test_read_vcf_all_breakend_forms_round_trip() {
        // Test all 4 VCF 4.2 breakend forms
        let alts = vec![
            "T[chr2:300[",
            "G]chr2:400]",
            "]chr1:100]C",
            "[chr1:200[A",
            "G.",
            ".A",
            "GTC[2:1000[",
        ];

        let vcf_content = format!(
            "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100\tbnd_U\tT\t{}\t60\tPASS\tSVTYPE=BND;MATEID=bnd_V\tGT\t0/1
chr1\t200\tbnd_W\tG\t{}\t60\tPASS\tSVTYPE=BND;MATEID=bnd_X\tGT\t0/1
chr2\t300\tbnd_V\tC\t{}\t60\tPASS\tSVTYPE=BND;MATEID=bnd_U\tGT\t0/1
chr2\t400\tbnd_X\tA\t{}\t60\tPASS\tSVTYPE=BND;MATEID=bnd_W\tGT\t0/1
chr1\t500\t.\tG\t{}\t60\tPASS\tSVTYPE=BND\tGT\t0/1
chr1\t600\t.\tA\t{}\t60\tPASS\tSVTYPE=BND\tGT\t0/1
chr2\t700\t.\tGTC\t{}\t60\tPASS\tSVTYPE=BND\tGT\t0/1
",
            alts[0], alts[1], alts[2], alts[3], alts[4], alts[5], alts[6]
        );

        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();

        let chr1_variants = records.get("chr1").unwrap();
        let chr2_variants = records.get("chr2").unwrap();

        assert_eq!(chr1_variants.len(), 4);
        assert_eq!(chr2_variants.len(), 3);

        // Round trip check
        let temp_dir = tempfile::tempdir().unwrap();
        let mut maps = HashMap::new();

        let mut chr1_records = Vec::new();
        for v in chr1_variants {
            chr1_records.push(v.clone());
        }
        let map1 = MutatedMap::from_interval(0, 1000, chr1_records).unwrap();
        maps.insert("chr1".to_string(), vec![map1]);

        let mut chr2_records = Vec::new();
        for v in chr2_variants {
            chr2_records.push(v.clone());
        }
        let map2 = MutatedMap::from_interval(0, 1000, chr2_records).unwrap();
        maps.insert("chr2".to_string(), vec![map2]);

        let output = temp_dir.path().join("bnd_rt.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string(), "chr2".to_string()],
            &HashMap::from([
                ("chr1".to_string(), 1000usize),
                ("chr2".to_string(), 1000usize),
            ]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();

        let body: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        for alt in alts {
            assert!(
                body.iter().any(|l| l.contains(alt)),
                "expected ALT {} to round-trip; got: {:?}",
                alt,
                body
            );
        }
    }

    #[test]
    fn test_write_vcf_bnd_variant_round_trips() {
        use crate::structs::variants::{AlternateType, Genotype, Provenance, SvData, SvType};
        let bnd_alt = "G]17:198982]";
        let sv = Variant {
            variant_type: VariantType::BND,
            location: 100,
            reference: vec![Nucleotide::G],
            alternate: AlternateType::Symbolic(SvData::new(bnd_alt, SvType::Bnd)),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: Some("MATEID=test".to_string()),
            format: vec!["GT".to_string()],
            sample: vec!["0/1".to_string()],
            provenance: Provenance::InputVcf,
        };

        let map = MutatedMap::from_interval(0, 1000, vec![sv]).unwrap();
        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);

        let temp_dir = tempfile::tempdir().unwrap();
        let output = temp_dir.path().join("bnd.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 1000usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();

        let body: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        let bnd_line = body
            .iter()
            .find(|l| l.contains(bnd_alt))
            .expect("BND ALT not found in output VCF");
        assert!(
            bnd_line.contains("MATEID=test"),
            "INFO field (MATEID) not preserved in output VCF: {}",
            bnd_line
        );
    }

    #[test]
    fn test_read_vcf_all_dot_genotype_treated_as_homozygous() {
        // Today `gt_from_str` returns Homozygous when no allele parses as zero — `./.` falls
        // into that bucket. This is arguably a bug (missing GT is not homozygous), but the
        // test pins the current behavior so any future fix is a deliberate, reviewed change.
        use crate::structs::variants::Genotype;
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t./.\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").unwrap();
        assert_eq!(variants.len(), 1);
        assert!(matches!(variants[0].genotype, Genotype::Homozygous));
    }

    #[test]
    fn test_read_vcf_mixed_phased_unphased_gt_both_supported() {
        // Both `0|1` (phased) and `0/1` (unphased) must parse and yield the same genotype.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0|1\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        use crate::structs::variants::Genotype;
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").unwrap();
        assert_eq!(variants.len(), 2);
        assert!(matches!(variants[0].genotype, Genotype::Heterozygous));
        assert!(matches!(variants[1].genotype, Genotype::Heterozygous));
        // Both genotype strings are stored verbatim (not normalized).
        assert_eq!(variants[0].genotype_str, "0/1");
        assert_eq!(variants[1].genotype_str, "0|1");
    }

    #[test]
    fn test_read_vcf_missing_gt_defaults_to_heterozygous() {
        // Records whose FORMAT lacks GT (e.g. only DP:GQ) used to be silently
        // dropped. v1.10 instead defaults them to heterozygous so sites-only
        // and partial-FORMAT VCFs train without losing every record. The
        // trainer's homozygous_frequency will then be 0 for the corpus —
        // an honest "we don't know" rather than a synthesized 50/50 split.
        let vcf_content = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tDP:GQ\t30:99\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").expect("chr1 should have one variant");
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].genotype_str, "0/1");
        assert!(matches!(
            variants[0].genotype,
            crate::structs::variants::Genotype::Heterozygous
        ));
    }

    #[test]
    fn test_read_vcf_sites_only_record_defaults_to_heterozygous() {
        // gnomAD-SV (and other sites-only callsets) carry only 8 columns —
        // no FORMAT, no SAMPLE. The reader must still produce a Variant per
        // record so the SV trainer can fit a model from a sites-only VCF.
        let vcf_content = "\
##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=500\n";
        let mut tmp = tempfile::Builder::new().suffix(".vcf").tempfile().unwrap();
        write!(tmp, "{}", vcf_content).unwrap();
        let records = read_vcf(tmp.path().to_path_buf()).unwrap();
        let variants = records.get("chr1").expect("chr1 should have one variant");
        assert_eq!(variants.len(), 1);
        assert!(variants[0].alternate.is_symbolic());
        assert_eq!(variants[0].genotype_str, "0/1");
    }

    #[test]
    fn test_write_vcf_symbolic_alt_round_trips_raw_string() {
        // Writer contract: an AlternateType::Symbolic payload serializes
        // verbatim via SvData.raw_alt. Driven from a synthetic Variant
        // (not the reader) so the writer is exercised in isolation —
        // the reader+writer round-trip is covered separately.
        let temp_dir = tempfile::tempdir().unwrap();
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 99, // 0-based; will appear as POS=100 in output
            reference: vec![Nucleotide::A],
            alternate: super::AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        let mutated_map = MutatedMap::from_interval(0, 200, vec![sv]).unwrap();
        let mut mutated_maps = HashMap::new();
        mutated_maps.insert("chr1".to_string(), vec![mutated_map]);
        let output_path = temp_dir.path().join("sv_round_trip.vcf");
        write_vcf(
            &mutated_maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 200usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output_path,
            &HashMap::new(),
        )
        .unwrap();
        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        // Symbolic SVs use span-based depth semantics that don't fit the
        // point-based AD/DP/AF counter — they emit `.` placeholders. The
        // `Provenance::Denovo` constructor on the test fixture surfaces as
        // INFO/EIDOLON_PROVENANCE=denovo.
        assert!(
            lines.iter().any(|l| l
                == "chr1\t100\t.\tA\t<DEL>\t37\tPASS\tEIDOLON_PROVENANCE=denovo\tGT:AD:DP:AF\t0/1:.,.:.:."),
            "expected symbolic ALT to round-trip as `<DEL>` with denovo provenance and `.` AD/DP/AF; got: {:?}",
            lines
        );
    }

    /// When the per-contig AdCounter has real (ref_count, alt_count) values
    /// at a variant's position, the writer must surface them as AD = a,b,
    /// DP = a+b, AF = b/(a+b) rounded to 4 decimal places. After the #185
    /// merge, every literal record also carries `EIDOLON_PROVENANCE=denovo`
    /// (since `Variant::new` is the de-novo constructor).
    #[test]
    fn test_write_vcf_emits_ad_dp_af_from_counter() {
        let temp_dir = tempfile::tempdir().unwrap();
        let v_het = Variant::new(
            VariantType::SNP,
            10,
            &vec![Nucleotide::A],
            &vec![Nucleotide::G],
            &mut vec![0, 1],
        )
        .unwrap();
        let v_hom = Variant::new(
            VariantType::SNP,
            20,
            &vec![Nucleotide::C],
            &vec![Nucleotide::T],
            &mut vec![1, 1],
        )
        .unwrap();
        let map = MutatedMap::from_interval(0, 100, vec![v_het, v_hom]).unwrap();
        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);

        // Stub AdCounter: het pos=10 → 18 ref, 12 alt (AF=0.4); hom pos=20 → 0,30 (AF=1.0)
        let mut ad_chr1: AdCounter = HashMap::new();
        ad_chr1.insert(10, (18, 12));
        ad_chr1.insert(20, (0, 30));
        let mut ad_counters = HashMap::new();
        ad_counters.insert("chr1".to_string(), ad_chr1);

        let output = temp_dir.path().join("ad_dp_af.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 100usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &ad_counters,
        )
        .unwrap();
        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();

        // Het: position 11 (1-based), AD=18,12, DP=30, AF=0.4000
        assert!(
            lines.iter().any(|l| l
                == "chr1\t11\t.\tA\tG\t37\tPASS\tEIDOLON_PROVENANCE=denovo\tGT:AD:DP:AF\t0/1:18,12:30:0.4000"),
            "missing het record with AD=18,12 DP=30 AF=0.4000; got: {:?}",
            lines
        );
        // Hom: position 21 (1-based), AD=0,30, DP=30, AF=1.0000
        assert!(
            lines.iter().any(|l| l
                == "chr1\t21\t.\tC\tT\t37\tPASS\tEIDOLON_PROVENANCE=denovo\tGT:AD:DP:AF\t1/1:0,30:30:1.0000"),
            "missing hom record with AD=0,30 DP=30 AF=1.0000; got: {:?}",
            lines
        );

        // Header lines must declare the three new FORMAT fields.
        assert!(
            lines.iter().any(|l| l.starts_with("##FORMAT=<ID=AD,")),
            "missing AD header"
        );
        assert!(
            lines.iter().any(|l| l.starts_with("##FORMAT=<ID=DP,")),
            "missing DP header"
        );
        assert!(
            lines.iter().any(|l| l.starts_with("##FORMAT=<ID=AF,")),
            "missing AF header"
        );
    }

    /// `Provenance::InputVcf` literal variants must surface as
    /// `EIDOLON_PROVENANCE=input` in the golden VCF — this is the signal the
    /// cancer-simulator merge step relies on to detect germline carry-through
    /// (a variant that was supplied via `input_vcf:` and appears in both
    /// passes is the "shared" case).
    #[test]
    fn test_write_vcf_input_provenance_surfaces_as_input_tag() {
        let temp_dir = tempfile::tempdir().unwrap();
        // Hand-roll the variant directly so we exercise the InputVcf path
        // without round-tripping through from_file (which has its own tests).
        let v = Variant {
            variant_type: VariantType::SNP,
            location: 49, // 0-based; will appear as POS=50 in output
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::G]),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::InputVcf,
        };
        let mutated_map = MutatedMap::from_interval(0, 200, vec![v]).unwrap();
        let mut mutated_maps = HashMap::new();
        mutated_maps.insert("chr1".to_string(), vec![mutated_map]);
        let output_path = temp_dir.path().join("input_prov.vcf");
        write_vcf(
            &mutated_maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 200usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output_path,
            &HashMap::new(),
        )
        .unwrap();
        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        // Empty ad_counter → AD=0,0 DP=0 AF=`.`.
        assert!(
            lines.iter().any(|l| l
                == "chr1\t50\t.\tA\tG\t37\tPASS\tEIDOLON_PROVENANCE=input\tGT:AD:DP:AF\t0/1:0,0:0:."),
            "expected EIDOLON_PROVENANCE=input on InputVcf literal; got: {:?}",
            lines
        );
        // Header line must declare the new INFO field so downstream
        // parsers (bcftools, hap.py, …) don't drop it as unknown.
        assert!(
            lines
                .iter()
                .any(|l| l.starts_with("##INFO=<ID=EIDOLON_PROVENANCE,")),
            "missing EIDOLON_PROVENANCE INFO header line; got: {:?}",
            lines
        );
        // The sample-column name is an emitted output token too (renamed from
        // NEAT_simulated_sample in v3.0.0). Without this assertion it is the one
        // renamed token no test would catch being reverted.
        assert!(
            lines
                .iter()
                .any(|l| l.starts_with("#CHROM") && l.ends_with("\tEIDOLON_simulated_sample")),
            "missing/incorrect EIDOLON_simulated_sample column header; got: {:?}",
            lines
        );
    }

    #[test]
    fn strip_own_output_tokens_drops_routing_keys_both_spellings() {
        // Third-party / SV INFO must survive untouched.
        assert_eq!(
            strip_own_output_tokens("SVTYPE=DEL;END=200;DP=30").as_deref(),
            Some("SVTYPE=DEL;END=200;DP=30")
        );
        // Pre-v3.0.0 keys removed, order of the rest preserved.
        assert_eq!(
            strip_own_output_tokens("SVTYPE=DEL;NEAT_PROVENANCE=denovo;END=200").as_deref(),
            Some("SVTYPE=DEL;END=200")
        );
        assert_eq!(
            strip_own_output_tokens(
                "NEAT_ORIGIN=somatic;NEAT_CCF=1.0;NEAT_VAF=0.3;NEAT_REASON=unknown;AF=0.5"
            )
            .as_deref(),
            Some("AF=0.5")
        );
        // Nothing left to keep -> None, so the caller emits provenance alone
        // rather than a record starting with ';'.
        assert_eq!(strip_own_output_tokens("NEAT_PROVENANCE=input"), None);

        // The CURRENT spellings must be stripped too. Round-tripping a v3.0.0 golden
        // through input_vcf: is the common case, not the legacy one.
        assert_eq!(
            strip_own_output_tokens("SVTYPE=DEL;EIDOLON_PROVENANCE=denovo;END=200").as_deref(),
            Some("SVTYPE=DEL;END=200")
        );
        assert_eq!(
            strip_own_output_tokens("SVTYPE=DEL;EIDOLON_ORIGIN=somatic").as_deref(),
            Some("SVTYPE=DEL")
        );
        // A stale EIDOLON_PROVENANCE=denovo surviving here is what made
        // vcf_merge::tumor_origin label a germline carry-through as `somatic`.
        assert_eq!(
            strip_own_output_tokens("EIDOLON_PROVENANCE=denovo;EIDOLON_ORIGIN=shared"),
            None
        );

        // #405 ground truth is KEPT: gen_reads::runner puts EIDOLON_CCF/EIDOLON_VAF in
        // `variant.info` before the write, so stripping them would destroy the tags on
        // every de-novo somatic record. This writer declares both.
        assert_eq!(
            strip_own_output_tokens("EIDOLON_CCF=1.0000;EIDOLON_VAF=0.3000").as_deref(),
            Some("EIDOLON_CCF=1.0000;EIDOLON_VAF=0.3000")
        );
        assert_eq!(
            strip_own_output_tokens("EIDOLON_CCF=1.0000;EIDOLON_PROVENANCE=denovo").as_deref(),
            Some("EIDOLON_CCF=1.0000")
        );

        // Exact-key matching: lookalikes are not owned tokens.
        assert_eq!(
            strip_own_output_tokens("XNEAT_ORIGIN=x;NEAT_ORIGINAL=y;EIDOLON_ORIGINS=z").as_deref(),
            Some("XNEAT_ORIGIN=x;NEAT_ORIGINAL=y;EIDOLON_ORIGINS=z")
        );
        // Bare flags (no '=') are preserved.
        assert_eq!(
            strip_own_output_tokens("IMPRECISE;NEAT_VAF=0.1").as_deref(),
            Some("IMPRECISE")
        );
    }

    /// Symbolic SVs typically arrive with a populated INFO column
    /// (SVLEN, END, SVTYPE, CN). The writer must *append* EIDOLON_PROVENANCE
    /// to that existing field rather than overwrite it — otherwise span
    /// metadata is lost on the symbolic record's round trip.
    #[test]
    fn test_write_vcf_symbolic_appends_provenance_to_existing_info() {
        let temp_dir = tempfile::tempdir().unwrap();
        let sv = Variant {
            variant_type: VariantType::Complex,
            location: 199, // 0-based; will appear as POS=200 in output
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: Some("SVTYPE=DEL;END=300;SVLEN=-100".to_string()),
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::InputVcf,
        };
        let mutated_map = MutatedMap::from_interval(0, 400, vec![sv]).unwrap();
        let mut mutated_maps = HashMap::new();
        mutated_maps.insert("chr1".to_string(), vec![mutated_map]);
        let output_path = temp_dir.path().join("sv_appended.vcf");
        write_vcf(
            &mutated_maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 400usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output_path,
            &HashMap::new(),
        )
        .unwrap();
        let lines: Vec<String> = crate::file_tools::file_io::read_gzip_lines(&output_path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        // Symbolic SVs always emit `.,.:.:.` for AD/DP/AF (span-based depth,
        // not point-based). Existing INFO must survive with EIDOLON_PROVENANCE
        // appended.
        assert!(
            lines.iter().any(|l| l
                == "chr1\t200\t.\tA\t<DEL>\t37\tPASS\t\
                 SVTYPE=DEL;END=300;SVLEN=-100;EIDOLON_PROVENANCE=input\tGT:AD:DP:AF\t1/1:.,.:.:."),
            "expected SVTYPE/END/SVLEN to survive with EIDOLON_PROVENANCE appended; got: {:?}",
            lines
        );
    }

    /// `write_vcf` must emit records in position-sorted order — `bcftools
    /// index -t` (tabix) requires it, and the cancer_simulate.sh merge
    /// step calls `bcftools index -t` on each per-pass output. Pre-fix the
    /// writer iterated the variant_map HashMap and emitted records in
    /// nondeterministic order, which silently broke the merge step.
    ///
    /// Mixes literals (HashMap-backed) and SVs (Vec-backed) to verify the
    /// sort runs across both bins, not just within each.
    #[test]
    fn test_write_vcf_records_emitted_in_position_sorted_order() {
        let temp_dir = tempfile::tempdir().unwrap();
        let v_500 = Variant::new(
            VariantType::SNP,
            500,
            &vec![Nucleotide::A],
            &vec![Nucleotide::G],
            &mut vec![0, 1],
        )
        .unwrap();
        let v_100 = Variant::new(
            VariantType::SNP,
            100,
            &vec![Nucleotide::C],
            &vec![Nucleotide::T],
            &mut vec![0, 1],
        )
        .unwrap();
        let v_300 = Variant::new(
            VariantType::SNP,
            300,
            &vec![Nucleotide::G],
            &vec![Nucleotide::A],
            &mut vec![0, 1],
        )
        .unwrap();
        // Symbolic SV at position 200 — must interleave with the literals,
        // not get bucketed to the end of the file.
        let sv_200 = Variant {
            variant_type: VariantType::Complex,
            location: 200,
            reference: vec![Nucleotide::T],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            genotype_str: "0/1".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: Some("SVTYPE=DEL;END=250".to_string()),
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };

        let map = MutatedMap::from_interval(0, 1000, vec![v_500, v_100, v_300, sv_200]).unwrap();
        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);

        let output = temp_dir.path().join("sorted.vcf");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 1000usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();

        let positions: Vec<usize> = crate::file_tools::file_io::read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap().parse::<usize>().unwrap())
            .collect();
        assert_eq!(
            positions,
            vec![101, 201, 301, 501],
            "records must be in position-sorted order regardless of insertion order"
        );
    }

    /// Pin BGZF framing on the output (issue #206). Plain gzip and BGZF are
    /// both gzip-compatible — `read_gzip_lines` (MultiGzDecoder under the
    /// hood) accepts either — so a regression to plain-gzip output would be
    /// invisible to most eidolon tests but would silently break the cancer-
    /// simulator merge step's `bcftools index -t` call. Check the two
    /// BGZF-defining markers directly:
    ///   - The first 4 bytes have FLG.FEXTRA=1 (gzip magic + 0x04 flag byte)
    ///   - The "BC" subfield identifier appears in the gzip extra field
    /// Plain gzip has FLG=0x00 and no extra field, so this test fails loudly
    /// if anyone swaps the writer back.
    #[test]
    fn test_write_vcf_output_is_bgzf_framed() {
        let temp_dir = tempfile::tempdir().unwrap();
        let map = MutatedMap::from_interval(0, 100, vec![]).unwrap();
        let mut maps = HashMap::new();
        maps.insert("chr1".to_string(), vec![map]);
        let output = temp_dir.path().join("bgzf_probe.vcf.gz");
        write_vcf(
            &maps,
            &vec!["chr1".to_string()],
            &HashMap::from([("chr1".to_string(), 100usize)]),
            &PathBuf::from("ref.fa"),
            false,
            &output,
            &HashMap::new(),
        )
        .unwrap();
        let mut bytes = [0u8; 16];
        use std::io::Read;
        std::fs::File::open(&output)
            .unwrap()
            .read_exact(&mut bytes)
            .unwrap();
        assert_eq!(
            &bytes[0..4],
            &[0x1f, 0x8b, 0x08, 0x04],
            "output should be BGZF (gzip magic + FLG.FEXTRA=4); got: {:02x?}",
            &bytes[0..4]
        );
        // bytes[10..12] is XLEN (LE u16), bytes[12..14] is the subfield ID
        // ("BC" for BGZF per the htslib spec).
        assert_eq!(
            &bytes[12..14],
            b"BC",
            "expected BGZF 'BC' subfield identifier at bytes 12-13; got: {:?}",
            std::str::from_utf8(&bytes[12..14])
        );
    }
}
