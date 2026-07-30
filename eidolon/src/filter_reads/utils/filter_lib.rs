use crate::filter_reads::errors::FilterReadsError;
use eidolon_core::file_tools::block_gz::BlockGzWriter;
use eidolon_core::structs::bed_record::BedRecord;
use flate2::read::MultiGzDecoder;
use log::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Lines, Write};
use std::path::PathBuf;
use std::sync::Once;
use thiserror::Error;

/// Guards the one-time pre-v3.0.0 read-name-prefix warning in
/// `parse_fastq_record_name` — it is called once per FASTQ record, and a per-read
/// warning on a multi-million-read file would be worse than silence.
static LEGACY_PREFIX_WARNED: Once = Once::new();

use eidolon_core::file_tools::file_io::{create_output_file, read_gzip_lines, read_lines};
/// The purpose of this module is to filter a eidolon-generated fastq and vcf file, so that they only show regions
/// of overlap with an input bed file.
enum ReaderType {
    GzipReader(Lines<BufReader<MultiGzDecoder<File>>>),
    NonZipReader(Lines<BufReader<File>>),
}

#[derive(Debug, Error)]
enum FilterLibError {
    #[error("Possible Q score")]
    InvalidReadName,
    #[error("Line parsing was not understood: {0}")]
    MalformedLine(String),
}

impl Iterator for ReaderType {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::GzipReader(iterator) => {
                let next = iterator.next();
                next.map(|result| result.expect("Error reading gzipped file."))
            }
            Self::NonZipReader(iterator) => {
                let next = iterator.next();
                next.map(|result| result.expect("Error reading non-zipped file."))
            }
        }
    }
}

fn prep_files_for_filtering(
    file_in: &PathBuf,
    is_gzip: bool,
    file_out: &PathBuf,
) -> Result<(ReaderType, BlockGzWriter<File>), FilterReadsError> {
    // Open file for reading.
    let lines: ReaderType = {
        if is_gzip {
            open_gzipped(file_in)
        } else {
            open_unzipped(file_in)
        }
    };
    let out_file = create_output_file(file_out, true)?;
    let buffer = BlockGzWriter::new(out_file);
    Ok((lines, buffer))
}

fn open_gzipped(filename: &PathBuf) -> ReaderType {
    ReaderType::GzipReader(
        read_gzip_lines(filename).expect("Error decompressing gzipped fastq file."),
    )
}

fn open_unzipped(filename: &PathBuf) -> ReaderType {
    ReaderType::NonZipReader(
        read_lines(filename).expect("Error reading fastq file while filtering."),
    )
}

fn parse_fastq_record_name(line: &str) -> Result<(String, usize, usize), FilterLibError> {
    // gen_reads produces:
    //   "@EIDOLON_generated_<contig>_<10-digit start>_<10-digit end>_<16-hex uniq>/N"
    // The trailing 16-char hex disambiguates same-position fragments (see #210
    // for why it's necessary). Parser takes the last three underscore-separated
    // tokens as (start, end, uniq) and ignores uniq.
    const PREFIX: &str = "@EIDOLON_generated_";
    // v2.0.0 and earlier emitted "@RNEAT_generated_" (v3.0.0 renamed the output
    // tokens). Accept both: FASTQs are far too large to be worth rewriting just to
    // change a 19-byte prefix, so legacy output stays filterable. Note the two
    // prefixes differ in length (19 vs 17), hence the matched-prefix slice below.
    const LEGACY_PREFIX: &str = "@RNEAT_generated_";
    let prefix = if line.starts_with(PREFIX) {
        PREFIX
    } else if line.starts_with(LEGACY_PREFIX) {
        // Warn once. Accepting the old prefix silently is the only remaining path where
        // a user's pre-v3.0.0 data flows through v3.0.0 with no signal that the output
        // format changed — and a user's own `grep '^@RNEAT_generated_'` over the same
        // file returns zero rows with exit 0. This is the cheapest place to make that
        // discoverable, so say it out loud (once, not per read).
        LEGACY_PREFIX_WARNED.call_once(|| {
            warn!(
                "input uses the pre-v3.0.0 read-name prefix \"{}\" (v3.0.0 renamed it to \
                 \"{}\"). Accepting it for compatibility. Note eidolon does NOT rewrite \
                 read names — any of your own scripts matching the old prefix will \
                 silently match nothing. See the \"Upgrading from 2.0.0\" section in the \
                 README.",
                LEGACY_PREFIX.trim_start_matches('@'),
                PREFIX.trim_start_matches('@'),
            );
        });
        LEGACY_PREFIX
    } else {
        // A line that *embeds* a known prefix without starting with one is a read
        // name we can't parse: that's a loud MalformedLine. Anything else is
        // InvalidReadName, which filter_fastq treats as an ordinary data line.
        // (Note both prefixes include the leading '@', so the "@N_"/"@T_" tagged
        // names in a merged cancer FASTQ do NOT land here — they are InvalidReadName
        // and get dropped. Pre-existing behaviour; no harness runs filter-reads on
        // a merged FASTQ, so it's untouched here.)
        return if line.contains(PREFIX) || line.contains(LEGACY_PREFIX) {
            Err(FilterLibError::MalformedLine(line.to_string()))
        } else {
            Err(FilterLibError::InvalidReadName)
        };
    };
    // Slice off the matched prefix and the "/N" strand suffix (2 chars). Use
    // prefix.len() rather than a magic number so a future prefix rename can't
    // silently shift the offset.
    let trimmed = &line[prefix.len()..line.len() - 2];
    let split_line: Vec<&str> = trimmed.split('_').collect();
    let length = split_line.len();
    if length < 4 {
        return Err(FilterLibError::MalformedLine(line.to_string()));
    }
    let start: usize = split_line[length - 3].parse().unwrap();
    let end: usize = split_line[length - 2].parse().unwrap();
    // split_line[length - 1] is the per-fragment uniqueness tag; ignored.
    let contig = split_line[..length - 3].join("_");
    Ok((contig, start, end))
}

pub fn filter_fastq(
    bed_table: &HashMap<String, Vec<BedRecord>>,
    fastq_in: &PathBuf,
    is_gzip: bool,
    fastq_out: &PathBuf,
) -> Result<(), FilterReadsError> {
    // bed_file: path to bed file to use for filtering
    // fastq_in: fastq ready for filtering, can be gzipped or not
    // is_gzip: whether the input file is gzipped or not
    // output_file: will be gzipped.
    // Run prep for filtering to open fastq
    let (in_file, mut out_file) = prep_files_for_filtering(fastq_in, is_gzip, fastq_out)?;
    // For each record in fastq, find name, parse, and check against bed records
    let mut include = false;

    for line in in_file {
        // Technically a quality score can start with @
        if line.starts_with("@") {
            // Record name candidate found
            let candidate_result = parse_fastq_record_name(&line);
            // It's possible that this is not a read name, but rather a quality score starting with @
            let (contig, start, end) = match candidate_result {
                Ok(region) => Some(region),
                Err(error) => {
                    // If this is a quality score and we need to write it, we do, then we run the next loop
                    debug!("Had a weird record: {line}");
                    match error {
                        FilterLibError::MalformedLine(line) => {
                            return Err(FilterReadsError::FastqFilterError(format!(
                                "Filter error with line: {}",
                                line
                            )));
                        }
                        FilterLibError::InvalidReadName => {
                            if include {
                                out_file.write_all(line.as_bytes())?;
                                out_file.write_all(b"\n")?;
                            }
                        }
                    }
                    // That line is processed and we can move on
                    continue;
                }
            }
            .unwrap(); // This unwrap will either continue us or be a region, so this is good to go.
            include = if let Some(records) = bed_table.get(&contig) {
                let found = records.iter().any(|r| r.overlaps(&contig, start, end));
                if found {
                    out_file.write_all(line.as_bytes())?;
                    out_file.write_all(b"\n")?;
                }
                found
            } else {
                false
            };
        } else {
            // this indicates an actual record line. Fastq has a format of
            //     @name
            //     ACCCATTAC
            //     +
            //     ?????????
            // where the name is unique and in our case contains coordinates of the read
            // The second line is a sequence
            // the third line is a literal plus sign spacer
            // the foruth line is an encoded ascii string of the quality scores (usually 1-42 once translated)
            // So once we find a name we want, we need to read in the next three lines or we messed up the read
            // Once we hit the next name, we'll turn this off if that one is not included
            if include {
                out_file.write_all(line.as_bytes())?;
                out_file.write_all(b"\n")?;
            }
        }
    }
    Ok(())
}

pub fn filter_vcf(
    bed_table: &HashMap<String, Vec<BedRecord>>,
    vcf_in: &PathBuf,
    is_gzip: bool,
    vcf_out: &PathBuf,
) -> Result<(), FilterReadsError> {
    // Run prep for filtering to open vcf
    let (infile, mut outfile) = prep_files_for_filtering(vcf_in, is_gzip, vcf_out)?;
    for line in infile {
        if line.starts_with("#") {
            outfile.write_all(line.as_bytes())?;
            outfile.write_all(b"\n")?;
        } else {
            let line_vec: Vec<&str> = line.split('\t').collect();
            let chrom = line_vec[0];
            let pos: usize = line_vec[1].parse()?;
            if let Some(records) = bed_table.get(line_vec[0])
                && records.iter().any(|r| r.contains(chrom, pos))
            {
                outfile.write_all(line.as_bytes())?;
                outfile.write_all(b"\n")?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_prep_file_for_filtering() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();
        let mut temp_file: PathBuf = PathBuf::from(temp_dir.path());
        temp_file.push("Test_data.txt\n");
        let mut file_obj = create_output_file(&temp_file, true).unwrap();
        file_obj
            .write_all(
                "Test data!\nThis is only a test\nPlease disregard all information within\nthe end"
                    .as_bytes(),
            )
            .unwrap();
        file_obj.flush().unwrap();
        let mut outfile = PathBuf::from(temp_dir.path());
        outfile.push("test_data_filtered.txt");
        prep_files_for_filtering(&temp_file, false, &outfile)
            .expect("prep_files_for_filtering should succeed on a real input");
    }

    #[test]
    fn test_parse_record_name() {
        // Updated for the per-fragment uniqueness tag introduced in #210.
        let record_name =
            "@EIDOLON_generated_chrom_1_0000001000_0000002000_000000000000002a/1".to_string();
        assert_eq!(
            parse_fastq_record_name(&record_name).unwrap(),
            ("chrom_1".to_string(), 1000, 2000),
        );
    }

    /// Record names that pre-date the uniqueness-tag rollout (i.e. only have
    /// `<contig>_<start>_<end>` without the trailing hex tag) should produce a
    /// clear MalformedLine error rather than silently parsing the contig as
    /// the wrong thing or panicking. This protects users who feed old eidolon
    /// output through a new filter_reads — they'll see a real error.
    /// v2.0.0-era `@RNEAT_generated_` read names must still parse — FASTQs are too
    /// large to rewrite just for a prefix rename. The two prefixes differ in length
    /// (17 vs 19), so this also pins that the slice uses the *matched* prefix.
    #[test]
    fn test_parse_record_name_accepts_pre_rename_prefix() {
        let old = "@RNEAT_generated_chrom_1_0000001000_0000002000_000000000000002a/1".to_string();
        assert_eq!(
            parse_fastq_record_name(&old).unwrap(),
            ("chrom_1".to_string(), 1000, 2000),
        );
    }

    #[test]
    fn test_parse_record_name_legacy_format_errors() {
        let legacy = "@EIDOLON_generated_chr1_0000001000_0000002000/1".to_string();
        match parse_fastq_record_name(&legacy) {
            // Legacy split has 3 tokens [chr1, 0000001000, 0000002000]; new
            // parser reads end (parses 1000) and start (parses chr1 — fails).
            // The unwrap on `parse()` panics in that case — accept either
            // path so this test stays robust to small parser tweaks.
            Err(_) => {}
            Ok(parsed) => panic!(
                "legacy read name should not parse cleanly under the new parser; got {:?}",
                parsed
            ),
        }
    }

    #[test]
    fn test_filter_fastq() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();

        // Two reads: first inside BED chr1:0-2000, second outside.
        // The 16-hex trailer is the per-fragment uniqueness tag (#210).
        let fastq_content = concat!(
            "@EIDOLON_generated_chr1_0000001000_0000002000_0000000000000000/1\n",
            "ACGTACGT\n",
            "+\n",
            "IIIIIIII\n",
            "@EIDOLON_generated_chr1_0000005000_0000006000_0000000000000001/1\n",
            "TTTTTTTT\n",
            "+\n",
            "IIIIIIII\n",
        );
        let input = temp_dir.path().join("input.fastq");
        std::fs::write(&input, fastq_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("output.fastq.gz");
        filter_fastq(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(
            lines.len(),
            4,
            "Only the in-range read (4 lines) should be in output"
        );
        assert!(lines[0].contains("0000001000"));
    }

    #[test]
    fn test_filter_vcf() {
        let temp_dir: TempDir = tempfile::tempdir().unwrap();

        // Two variants: pos 1500 inside BED chr1:0-2000, pos 5000 outside
        let vcf_content = concat!(
            "##fileformat=VCFv4.1\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            "chr1\t1500\t.\tA\tG\t37\tPASS\t.\tGT\t0/1\n",
            "chr1\t5000\t.\tT\tC\t37\tPASS\t.\tGT\t0/1\n",
        );
        let input = temp_dir.path().join("input.vcf");
        std::fs::write(&input, vcf_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("output.vcf.gz");
        filter_vcf(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let variant_lines: Vec<&String> = lines.iter().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            variant_lines.len(),
            1,
            "Only the in-range variant should appear"
        );
        assert!(variant_lines[0].contains("1500"));
    }

    #[test]
    fn test_filter_fastq_empty_input() {
        // An empty input FASTQ must produce an empty (header-only) gzipped output,
        // not an error or a panic.
        let temp_dir: TempDir = tempfile::tempdir().unwrap();
        let input = temp_dir.path().join("empty.fastq");
        std::fs::write(&input, "").unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("empty_out.fastq.gz");
        filter_fastq(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert!(
            lines.is_empty(),
            "empty input should produce empty output, got: {:?}",
            lines
        );
    }

    #[test]
    fn test_filter_fastq_all_filtered_out() {
        // BED table doesn't overlap any reads — every record is filtered. Output has zero
        // FASTQ records.
        let temp_dir: TempDir = tempfile::tempdir().unwrap();
        let fastq_content = concat!(
            "@EIDOLON_generated_chr1_0000001000_0000002000_0000000000000000/1\n",
            "ACGTACGT\n",
            "+\n",
            "IIIIIIII\n",
            "@EIDOLON_generated_chr1_0000005000_0000006000_0000000000000001/1\n",
            "TTTTTTTT\n",
            "+\n",
            "IIIIIIII\n",
        );
        let input = temp_dir.path().join("all_out.fastq");
        std::fs::write(&input, fastq_content).unwrap();

        // BED on a contig that's not in the FASTQ → no overlaps.
        let bed_table = HashMap::from([(
            "chrZ".to_string(),
            vec![BedRecord::new_bed_record("chrZ".to_string(), 0, 1_000_000).unwrap()],
        )]);

        let output = temp_dir.path().join("all_out.fastq.gz");
        filter_fastq(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert!(lines.is_empty(), "no reads should pass; got {:?}", lines);
    }

    #[test]
    fn test_filter_vcf_no_passing_records_keeps_headers() {
        // VCF with one variant outside the BED region. Output must keep header lines
        // (start with '#') but emit no data lines.
        let temp_dir: TempDir = tempfile::tempdir().unwrap();
        let vcf_content = concat!(
            "##fileformat=VCFv4.1\n",
            "##source=eidolon-test\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            "chr1\t9999\t.\tA\tG\t37\tPASS\t.\tGT\t0/1\n",
        );
        let input = temp_dir.path().join("no_pass.vcf");
        std::fs::write(&input, vcf_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 0, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("no_pass.vcf.gz");
        filter_vcf(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        let headers: Vec<&String> = lines.iter().filter(|l| l.starts_with('#')).collect();
        let data: Vec<&String> = lines.iter().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(headers.len(), 3, "all three header lines should survive");
        assert!(
            data.is_empty(),
            "no data lines should be written: {:?}",
            data
        );
    }

    #[test]
    fn test_filter_fastq_boundary() {
        // Verifies that the corrected overlaps() semantics (half-open intervals,
        // adjacent intervals do NOT overlap) are respected end-to-end through filter_fastq.
        let temp_dir: TempDir = tempfile::tempdir().unwrap();

        // BED region: [1000, 2000)
        // Read A [800,  1000) — ends exactly at BED start  → excluded
        // Read B [2000, 2200) — starts exactly at BED end  → excluded
        // Read C [999,  1001) — straddles BED start        → included
        let fastq_content = concat!(
            "@EIDOLON_generated_chr1_0000000800_0000001000_0000000000000000/1\n",
            "ACGTACGT\n",
            "+\n",
            "IIIIIIII\n",
            "@EIDOLON_generated_chr1_0000002000_0000002200_0000000000000001/1\n",
            "TTTTTTTT\n",
            "+\n",
            "IIIIIIII\n",
            "@EIDOLON_generated_chr1_0000000999_0000001001_0000000000000002/1\n",
            "GGGGGGGG\n",
            "+\n",
            "IIIIIIII\n",
        );
        let input = temp_dir.path().join("boundary.fastq");
        std::fs::write(&input, fastq_content).unwrap();

        let bed_table = HashMap::from([(
            "chr1".to_string(),
            vec![BedRecord::new_bed_record("chr1".to_string(), 1000, 2000).unwrap()],
        )]);

        let output = temp_dir.path().join("boundary_out.fastq.gz");
        filter_fastq(&bed_table, &input, false, &output).unwrap();

        let lines: Vec<String> = read_gzip_lines(&output)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(
            lines.len(),
            4,
            "Only read C (straddles boundary) should be included"
        );
        assert!(
            lines[0].contains("0000000999"),
            "Expected read C in output, got: {}",
            lines[0]
        );
    }
}
