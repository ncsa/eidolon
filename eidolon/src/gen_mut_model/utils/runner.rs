use log::*;
use std::{collections::HashMap, path::PathBuf};

use eidolon_core::{
    file_tools::fasta_stream::{FastaStream, non_n_regions},
    models::{mutation_model::MutationModel, snp_trinuc_model::TrinucFrame},
    structs::{
        bed_record::BedRecord,
        nucleotides::Nucleotide,
        sv_model::SvModel,
        variants::{Genotype, Variant, VariantType},
    },
};

use crate::gen_mut_model::errors::GenMutationModelError;

pub fn runner(
    reference: &PathBuf,
    filtered_mutations: HashMap<String, Vec<Variant>>,
    bed_table: HashMap<String, Vec<BedRecord>>,
    output_file: &PathBuf,
) -> Result<(), GenMutationModelError> {
    let mut trinuc_count: HashMap<TrinucFrame, usize> = HashMap::new();
    let mut trinuc_transition_count: HashMap<(TrinucFrame, TrinucFrame), usize> = HashMap::new();
    let mut snp_count = 0;
    let mut insertion_count: HashMap<usize, usize> = HashMap::new();
    let mut deletion_count: HashMap<usize, usize> = HashMap::new();
    let mut homozygous_count = 0;
    let mut bed_track_len: usize = 0;
    let mut total_reflen: usize = 0;
    // Symbolic / structural ALTs (`<DEL>`, `<DUP>`, `<CNV>`, ...) bypass
    // the per-base SNP/indel accounting entirely and feed
    // `SvModel::fit_from_observations` after the main loop.
    let mut sv_observations: Vec<Variant> = Vec::new();

    let use_bed = !bed_table.is_empty();

    for result in FastaStream::open(reference)? {
        let (contig_name, raw) = result?;
        // IUPAC codes map to N here intentionally — model-building from real VCF data
        // doesn't need stochastic resolution because variant callers skip ambiguous positions.
        let sequence: Vec<Nucleotide> = raw.chars().map(Nucleotide::from).collect();

        let non_n = non_n_regions(&sequence);

        // Accumulate total_reflen over non-N regions regardless of BED mode.
        for &(start, end) in &non_n {
            total_reflen += end - start;
        }

        // Count trinucleotides.
        if use_bed {
            if let Some(bed_regions) = bed_table.get(&contig_name) {
                for region in bed_regions {
                    let r_start = region.start;
                    let r_end = region.end.min(sequence.len());
                    if r_end.saturating_sub(r_start) < 3 {
                        continue;
                    }
                    bed_track_len += r_end - r_start;
                    for i in (r_start + 1)..(r_end - 1) {
                        let frame =
                            TrinucFrame::from((sequence[i - 1], sequence[i], sequence[i + 1]));
                        *trinuc_count.entry(frame).or_default() += 1;
                    }
                }
            }
        } else {
            for &(start, end) in &non_n {
                for i in (start + 1)..(end - 1) {
                    let frame = TrinucFrame::from((sequence[i - 1], sequence[i], sequence[i + 1]));
                    *trinuc_count.entry(frame).or_default() += 1;
                }
            }
        }

        // Process variants for this contig.
        let matching_variants = match filtered_mutations.get(&contig_name) {
            Some(v) if !v.is_empty() => v,
            _ => {
                debug!("No variants found for {}", contig_name);
                continue;
            }
        };

        for variant in matching_variants {
            // Symbolic / structural ALTs (`<DEL>`, `<DUP>`, `<CNV>`, ...)
            // have no literal base content for the trinucleotide / indel
            // stats to consume; instead they go through the SV model fit
            // after the main loop. Routing them here keeps the
            // homozygous_count and SNP/indel totals clean.
            if variant.alternate.is_symbolic() {
                sv_observations.push(variant.clone());
                continue;
            }
            match variant.variant_type {
                VariantType::SNP => {
                    snp_count += 1;
                    // VCF POS is 1-based; skip variants too close to contig edges.
                    if variant.location < 2 {
                        debug!("Skipping edge variant at position {}", variant.location);
                        continue;
                    }
                    let loc = variant.location - 1; // 0-based
                    if loc + 1 >= sequence.len() {
                        debug!(
                            "Skipping edge variant at position {} (out of contig bounds)",
                            variant.location
                        );
                        continue;
                    }
                    let canon = |n: Nucleotide| match n {
                        Nucleotide::Maskeda => Nucleotide::A,
                        Nucleotide::Maskedc => Nucleotide::C,
                        Nucleotide::Maskedg => Nucleotide::G,
                        Nucleotide::Maskedt => Nucleotide::T,
                        other => other,
                    };
                    let n0 = canon(sequence[loc - 1]);
                    let n1 = canon(sequence[loc]);
                    let n2 = canon(sequence[loc + 1]);
                    if n1 != variant.reference[0] {
                        warn!(
                            "Reference mismatch at position {}: VCF ref {:?}, FASTA base {:?}; skipping",
                            variant.location, variant.reference[0], n1
                        );
                        continue;
                    }
                    let ref_frame = TrinucFrame::from((n0, n1, n2));
                    debug_assert!(
                        variant.alternate.is_literal(),
                        "symbolic ALT reached SNP arm in gen_mut_model"
                    );
                    let alt = variant.alternate.as_literal().unwrap();
                    let alt_frame = TrinucFrame::from((n0, alt[0], n2));
                    *trinuc_transition_count
                        .entry((ref_frame, alt_frame))
                        .or_default() += 1;
                }
                VariantType::Insertion => {
                    debug_assert!(
                        variant.alternate.is_literal(),
                        "symbolic ALT reached Insertion arm in gen_mut_model"
                    );
                    let variant_len =
                        variant.alternate.as_literal().unwrap().len() - variant.reference.len();
                    *insertion_count.entry(variant_len).or_default() += 1;
                }
                VariantType::Deletion => {
                    debug_assert!(
                        variant.alternate.is_literal(),
                        "symbolic ALT reached Deletion arm in gen_mut_model"
                    );
                    let variant_len =
                        variant.reference.len() - variant.alternate.as_literal().unwrap().len();
                    *deletion_count.entry(variant_len).or_default() += 1;
                }
                _ => debug!("Unknown variant type, skipping for this analysis."),
            }
            match variant.genotype {
                Genotype::Homozygous => homozygous_count += 1,
                Genotype::Heterozygous => {}
            }
        }
    }

    if trinuc_count.is_empty() {
        error!("Trinucleotide counts were empty");
        return Err(GenMutationModelError::TrinucCountError(
            "Trinuc counts are empty. Unknown error".to_string(),
        ));
    }

    // Pre-group transition counts by reference frame to avoid an O(n²) scan.
    let mut trans_by_ref: HashMap<TrinucFrame, HashMap<TrinucFrame, usize>> = HashMap::new();
    for (&(ref_f, alt_f), &count) in &trinuc_transition_count {
        trans_by_ref.entry(ref_f).or_default().insert(alt_f, count);
    }

    // Compute probabilities.
    let mut trinuc_mut_prob: HashMap<TrinucFrame, f64> = HashMap::new();
    let mut trinuc_trans_prob: HashMap<(TrinucFrame, TrinucFrame), f64> = HashMap::new();

    for (frame, count) in &trinuc_count {
        if *count == 0 {
            trinuc_mut_prob.insert(*frame, 0.0);
            continue;
        }
        let alts = trans_by_ref.get(frame);
        let frame_count: usize = alts.map_or(0, |m| m.values().sum());
        trinuc_mut_prob.insert(*frame, (frame_count as f64) / (*count as f64));
        if frame_count > 0
            && let Some(alts) = alts
        {
            for (&alt_f, &tc) in alts {
                trinuc_trans_prob.insert((*frame, alt_f), (tc as f64) / (frame_count as f64));
            }
        }
    }

    let total_insertions: usize = insertion_count.values().sum();
    let total_deletions: usize = deletion_count.values().sum();
    let allowed_variant_count = (snp_count + total_insertions + total_deletions) as f64;

    // Guard the SNP/indel frequency math: without literal variants every
    // ratio is 0/0 = NaN, which then serializes as JSON `null` and fails
    // to deserialize back as f64. Two paths through:
    //   - corpus had literal variants → compute the real frequencies
    //   - corpus is SV-only (e.g. gnomAD-SV) → emit a model with
    //     mutation_rate = 0 and uniform-fallback variant_probs so the
    //     SvModel fit downstream can still proceed. gen-reads loaded
    //     with this model will generate zero SNPs/indels and rely
    //     entirely on de novo SVs (via sv_rate_scale > 0) or
    //     input_vcf records.
    let sv_only_corpus = allowed_variant_count == 0.0;
    if sv_only_corpus {
        info!(
            "Training VCF has no SNP / insertion / deletion observations \
             (saw {} symbolic SV(s) only). Producing an SV-only model with \
             mutation_rate=0 — gen-reads loaded with this model will generate \
             zero SNPs/indels and rely entirely on SVs or input_vcf records.",
            filtered_mutations.values().map(|v| v.len()).sum::<usize>()
        );
    }

    let (variant_probs, homozygous_frequency, average_mutation_rate) = if sv_only_corpus {
        // Uniform fallback for variant_probs — `DiscreteDistribution::new`
        // rejects all-zero weights, so we feed it 1/3 each. Doesn't
        // matter at sample time because mutation_rate is 0.
        (vec![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], 0.0, 0.0)
    } else {
        let snp_freq = (snp_count as f64) / allowed_variant_count;
        let del_freq = (total_deletions as f64) / allowed_variant_count;
        let ins_freq = (total_insertions as f64) / allowed_variant_count;
        let probs = vec![snp_freq, ins_freq, del_freq];
        let hom_freq = if homozygous_count > 0 {
            (homozygous_count as f64) / allowed_variant_count
        } else {
            0.001 / allowed_variant_count
        };
        let mut_rate = if use_bed {
            allowed_variant_count / (bed_track_len as f64)
        } else {
            allowed_variant_count / (total_reflen as f64)
        };
        (probs, hom_freq, mut_rate)
    };

    let ins_lengths: Vec<usize> = insertion_count.keys().cloned().collect();
    let ins_weights: Vec<f64> = insertion_count.values().map(|&x| x as f64).collect();
    let del_lengths: Vec<usize> = deletion_count.keys().cloned().collect();
    let del_weights: Vec<f64> = deletion_count.values().map(|&x| x as f64).collect();

    let result = MutationModel::from_raw_data(
        average_mutation_rate,
        homozygous_frequency,
        variant_probs,
        trinuc_mut_prob,
        trinuc_trans_prob,
        ins_lengths,
        ins_weights,
        del_lengths,
        del_weights,
    );

    match result {
        Ok(mut mut_model) => {
            // Fit the optional SV component from the parallel accumulator.
            // The fitter returns `None` when there isn't enough signal
            // (no observations, every type too thin, etc.) — leaving
            // `sv_model = None` matches the v1.9 model shape so older
            // gen-reads builds keep loading the file unchanged.
            if !sv_observations.is_empty() {
                let sv_denom = if use_bed { bed_track_len } else { total_reflen };
                mut_model.sv_model = SvModel::fit_from_observations(&sv_observations, sv_denom);
            }
            mut_model.write_to_file(output_file)?;
            info!("Mutation model success! Wrote model to {:?}", output_file);
            Ok(())
        }
        Err(error) => Err(GenMutationModelError::MutModelError(error)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eidolon_core::{
        file_tools::vcf_tools::read_vcf, models::mutation_model::MutationModel,
        structs::bed_record::BedRecord, structs::variants::SvType,
    };
    use tempfile::tempdir;

    #[test]
    fn test_runner_with_snps() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let vcf_path = PathBuf::from(format!("{}/test_data/vcfs/small_snps.vcf", manifest_dir));
        let mutations = read_vcf(vcf_path).unwrap();
        let out_dir = tempdir().unwrap();
        let output_file = out_dir.path().join("test_model.json.gz");
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());
        let model = MutationModel::from_file(&output_file).unwrap();
        assert!(
            model.mutation_rate > 0.0,
            "mutation_rate should be positive"
        );
        assert!(
            model.homozygous_frequency > 0.0,
            "homozygous_frequency should be positive"
        );
    }

    #[test]
    fn test_runner_with_indels() {
        // VCF containing an insertion and a deletion alongside SNPs. The runner must accept
        // non-SNP variant types without erroring and produce a writable model. Catches any
        // regression that drops the insertion/deletion match arms.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let out_dir = tempdir().unwrap();

        // Position 22 of H1N1_HA is 'C' (after 19 leading Ns: T T C ...). Build a VCF with
        // a SNP, an insertion, and a deletion at sensible coordinates.
        let vcf_path = out_dir.path().join("indels.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\n\
H1N1_HA\t50\t.\tT\tTAG\t60\tPASS\t.\tGT\t0/1\n\
H1N1_HA\t80\t.\tACG\tA\t60\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();

        let mutations = read_vcf(vcf_path).unwrap();
        let output_file = out_dir.path().join("indel_model.json.gz");
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());

        let model = MutationModel::from_file(&output_file).unwrap();
        assert!(model.mutation_rate > 0.0);
    }

    #[test]
    fn test_runner_skips_reference_mismatch_variant() {
        // When VCF REF doesn't agree with the FASTA base at the same position, the runner
        // logs a warning and skips that variant. As long as at least one good variant
        // remains, the model still builds.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let out_dir = tempdir().unwrap();

        // First record matches (pos 22 is C). Second record's REF=A disagrees with FASTA.
        let vcf_path = out_dir.path().join("mismatch.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\n\
H1N1_HA\t25\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();

        let mutations = read_vcf(vcf_path).unwrap();
        let output_file = out_dir.path().join("mismatch_model.json.gz");
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());
    }

    #[test]
    fn test_runner_bed_unknown_contig_succeeds_or_errors_cleanly() {
        // BED entry references a contig that does not exist in the FASTA. The runner must
        // not panic; it either succeeds with bed_track_len=0 (which makes mutation_rate
        // diverge — usually surfaces as a model-construction error) or errors cleanly.
        // Whatever the exact behavior is today, lock it in.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let vcf_path = PathBuf::from(format!("{}/test_data/vcfs/small_snps.vcf", manifest_dir));
        let mutations = read_vcf(vcf_path).unwrap();
        let out_dir = tempdir().unwrap();
        let output_file = out_dir.path().join("bed_unknown.json.gz");

        let bed_record = BedRecord::new_bed_record("chrZ_nonexistent".to_string(), 1, 100).unwrap();
        let bed_table = HashMap::from([("chrZ_nonexistent".to_string(), vec![bed_record])]);

        // We don't enforce one outcome here — just that it terminates without panicking.
        let result = runner(&reference, mutations, bed_table, &output_file);
        let _ = result; // either Ok or Err is acceptable today
    }

    #[test]
    fn test_runner_with_bed_table() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let vcf_path = PathBuf::from(format!("{}/test_data/vcfs/small_snps.vcf", manifest_dir));
        let mutations = read_vcf(vcf_path).unwrap();
        let out_dir = tempdir().unwrap();
        let output_file = out_dir.path().join("bed_model.json.gz");

        // BED region covering SNP positions 22, 25, 28 on H1N1_HA (1-based VCF coords)
        let bed_record = BedRecord::new_bed_record("H1N1_HA".to_string(), 1, 100).unwrap();
        let bed_table = HashMap::from([("H1N1_HA".to_string(), vec![bed_record])]);

        runner(&reference, mutations, bed_table, &output_file).unwrap();

        assert!(output_file.exists());
        let model = MutationModel::from_file(&output_file).unwrap();
        assert!(
            model.mutation_rate > 0.0,
            "mutation_rate should be positive"
        );
        // 3 variants / 99 bp ≈ 0.030; whole-reference denominator is ~14 000 bp total
        assert!(
            model.mutation_rate > 0.01,
            "BED-constrained rate should exceed 0.01, got {}",
            model.mutation_rate
        );
    }

    #[test]
    fn test_runner_skips_symbolic_variants_without_panicking() {
        // A mixed input VCF (literal SNP + symbolic <DEL>) used to risk a panic
        // at as_literal().unwrap() if a symbolic record ever reached the SNP /
        // Insertion / Deletion arms. They're tagged VariantType::Complex today
        // so the match arm wouldn't fire, but the homozygous_count tally would
        // still count them and bias the model. The explicit symbolic-skip at
        // the top of the loop must drop them before any per-base accounting.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let out_dir = tempdir().unwrap();

        let vcf_path = out_dir.path().join("mixed_sv.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\n\
H1N1_HA\t100\t.\tA\t<DEL>\t60\tPASS\tEND=500\tGT\t1/1\n\
H1N1_HA\t600\t.\tT\t<DUP>\t60\tPASS\tEND=700\tGT\t1/1\n",
        )
        .unwrap();

        let mutations = read_vcf(vcf_path).unwrap();
        let output_file = out_dir.path().join("mixed_sv_model.json.gz");
        // Must not panic — symbolic records skip past the as_literal().unwrap()
        // sites and never reach the homozygous_count tally.
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());
        let model = MutationModel::from_file(&output_file).unwrap();
        // Only the SNP contributes — and it's heterozygous, so homozygous_count
        // stays at 0 and the model falls back to the tiny default frequency
        // (0.001 / allowed_variant_count) rather than picking up the symbolic homs.
        assert!(model.mutation_rate > 0.0);
        assert!(
            model.homozygous_frequency < 0.5,
            "symbolic homozygous SVs must not be counted as homozygous SNPs; got {}",
            model.homozygous_frequency
        );
        // With only one DEL and one DUP observed, both types fall below
        // the 2-observation per-type fit bar and the trainer writes
        // sv_model = None rather than a half-populated stub.
        assert!(
            model.sv_model.is_none(),
            "sv_model should be None when no type has enough observations; got {:?}",
            model.sv_model
        );
    }

    #[test]
    fn test_runner_fits_sv_model_from_sv_rich_vcf() {
        // Train on a VCF with enough <DEL> / <DUP> / <CNV> / <BND> observations
        // to clear the per-type fit bar. The produced MutationModel must
        // carry a populated sv_model whose distributions match the input
        // counts (within tolerance for the log-normal fit).
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let out_dir = tempdir().unwrap();

        let vcf_path = out_dir.path().join("sv_rich.vcf");
        // Includes a SNP so the SNP/indel side of the model still
        // builds — feeding the trainer a VCF with no literal variants
        // would NaN out the SNP-frequency math. Real training corpora
        // never look like that.
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\n\
H1N1_HA\t100\t.\tA\t<DEL>\t60\tPASS\tEND=500\tGT\t1/1\n\
H1N1_HA\t1000\t.\tA\t<DEL>\t60\tPASS\tEND=1200\tGT\t0/1\n\
H1N1_HA\t1500\t.\tA\t<DEL>\t60\tPASS\tEND=2000\tGT\t0/1\n\
H1N1_HA\t2500\t.\tT\t<DUP>\t60\tPASS\tEND=2700\tGT\t1/1\n\
H1N1_HA\t3000\t.\tT\t<DUP>\t60\tPASS\tEND=3300\tGT\t0/1\n\
H1N1_HA\t4000\t.\tA\t<CNV>\t60\tPASS\tEND=4500;CN=3\tGT\t1/1\n\
H1N1_HA\t5000\t.\tA\t<CNV>\t60\tPASS\tEND=5500;CN=0\tGT\t1/1\n\
H1N1_HA\t6000\t.\tG\t<BND>\t60\tPASS\t.\tGT\t1/1\n\
H1N1_HA\t7000\t.\tG\t<BND>\t60\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();

        let mutations = read_vcf(vcf_path).unwrap();
        let output_file = out_dir.path().join("sv_rich_model.json.gz");
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());

        let model = MutationModel::from_file(&output_file).unwrap();
        let sv = model
            .sv_model
            .as_ref()
            .expect("sv_model must be populated when all types clear the fit bar");
        assert!(sv.is_usable());
        // 3 DELs + 2 DUPs + 2 CNVs + 2 BNDs, all clearing the bar — all four survive.
        assert_eq!(sv.type_probabilities.len(), 4);
        assert!(sv.type_probabilities.contains_key(&SvType::Del));
        assert!(sv.type_probabilities.contains_key(&SvType::Dup));
        assert!(sv.type_probabilities.contains_key(&SvType::Cnv));
        assert!(sv.type_probabilities.contains_key(&SvType::Bnd));
        // Each type also has a length distribution.
        assert!(sv.length_log_normal.contains_key(&SvType::Del));
        assert!(sv.length_log_normal.contains_key(&SvType::Dup));
        assert!(sv.length_log_normal.contains_key(&SvType::Cnv));
        assert!(sv.length_log_normal.contains_key(&SvType::Bnd));
        // BND length distribution should be fixed at mu=0, sigma=0.
        let (mu, sigma) = sv.length_log_normal[&SvType::Bnd];
        assert_eq!(mu, 0.0);
        assert_eq!(sigma, 0.0);
        // CN distribution: CN=0 and CN=3 each at 1/2.
        assert_eq!(sv.cnv_copy_number_distribution.len(), 2);
        assert!((sv.cnv_copy_number_distribution[&0] - 0.5).abs() < 1e-12);
        assert!((sv.cnv_copy_number_distribution[&3] - 0.5).abs() < 1e-12);
        // 5 homozygous of 9 records.
        assert!((sv.homozygous_frequency - 5.0 / 9.0).abs() < 1e-6);
        // Per-base rate is positive and tiny (7 events over ~14kb of
        // H1N1 reference).
        assert!(sv.per_base_rate > 0.0);
        assert!(sv.per_base_rate < 1e-2);
    }

    #[test]
    fn test_runner_builds_sv_only_model_from_vcf_with_no_literals() {
        // gnomAD-SV and similar SV-dedicated callsets carry symbolic SVs
        // but zero literal SNPs/indels. v1.9.1 errored out on this path
        // (NoLiteralVariants) defensively against NaN-poisoned models;
        // v1.10 instead produces an SV-only model with mutation_rate=0
        // and a populated sv_model, so the trainer is usable against
        // real-world SV callsets.
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let reference = PathBuf::from(format!("{}/test_data/references/H1N1.fa", manifest_dir));
        let out_dir = tempdir().unwrap();

        // Two DELs (enough to clear the SvModel per-type fit bar of 2
        // observations) and a DUP+CNV that fall below the bar and get
        // dropped from the sv_model itself but still count as
        // "observations seen". Fixture mimics gnomAD-SV shape (zero
        // literal SNP/indel records).
        let vcf_path = out_dir.path().join("sv_only.vcf");
        std::fs::write(
            &vcf_path,
            "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
H1N1_HA\t100\t.\tA\t<DEL>\t60\tPASS\tEND=200\tGT\t1/1\n\
H1N1_HA\t300\t.\tA\t<DEL>\t60\tPASS\tEND=450\tGT\t0/1\n\
H1N1_HA\t600\t.\tT\t<DUP>\t60\tPASS\tEND=700\tGT\t1/1\n",
        )
        .unwrap();

        let mutations = read_vcf(vcf_path).unwrap();
        let output_file = out_dir.path().join("sv_only_model.json.gz");
        runner(&reference, mutations, HashMap::new(), &output_file).unwrap();
        assert!(output_file.exists());

        // Re-read the produced model — must deserialize cleanly (no NaN
        // poisoning the JSON), have mutation_rate=0, homozygous_frequency=0,
        // and a populated sv_model with the DEL fit.
        let model = MutationModel::from_file(&output_file).unwrap();
        assert_eq!(model.mutation_rate, 0.0);
        assert_eq!(model.homozygous_frequency, 0.0);
        let sv = model
            .sv_model
            .as_ref()
            .expect("SV-only training must produce a populated sv_model");
        assert!(sv.is_usable());
        assert!(sv.type_probabilities.contains_key(&SvType::Del));
        // DUP had only 1 observation → dropped from the model.
        assert!(!sv.type_probabilities.contains_key(&SvType::Dup));
    }
}
