//! `filter-reads` on uncompressed input in a directory whose name contains a dot (#824).
//!
//! Output names were built by splitting the whole path on `.` and indexing from the end, which
//! assumed a `.fastq.gz`/`.vcf.gz` file in a dot-free path. A plain `reads.fastq` panicked
//! (`attempt to subtract with overflow`), and a dotted directory rebuilt the wrong input path.
//! This runs the binary on exactly that shape and reads the filtered records back.

mod common;
use common::{eidolon, fresh_workdir};
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};
use std::path::Path;

fn gz_lines(path: &Path) -> Vec<String> {
    let f = std::fs::File::open(path)
        .unwrap_or_else(|e| panic!("expected output {}: {e}", path.display()));
    BufReader::new(MultiGzDecoder::new(f))
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

#[test]
fn plain_input_in_a_dotted_directory_is_filtered_to_a_gz_output() {
    let (_tmp, work) = fresh_workdir();
    let dir = work.join("run.1");
    std::fs::create_dir(&dir).unwrap();

    // One read and one variant inside the BED region chr1:0-2000, one of each outside it.
    std::fs::write(
        dir.join("reads.fastq"),
        concat!(
            "@EIDOLON_generated_chr1_0000001000_0000002000_0000000000000000/1\n",
            "ACGTACGT\n+\nIIIIIIII\n",
            "@EIDOLON_generated_chr1_0000005000_0000006000_0000000000000001/1\n",
            "TTTTTTTT\n+\nIIIIIIII\n",
        ),
    )
    .unwrap();
    std::fs::write(
        dir.join("calls.vcf"),
        concat!(
            "##fileformat=VCFv4.1\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            "chr1\t1500\t.\tA\tG\t37\tPASS\t.\tGT\t0/1\n",
            "chr1\t5000\t.\tT\tC\t37\tPASS\t.\tGT\t0/1\n",
        ),
    )
    .unwrap();
    std::fs::write(dir.join("regions.bed"), "chr1\t0\t2000\n").unwrap();
    let config = dir.join("filter.yml");
    std::fs::write(
        &config,
        format!(
            "bed_file: {d}/regions.bed\nfiles_to_filter:\n  - {d}/reads.fastq\n  - {d}/calls.vcf\n\
             filter_key: .\noverwrite_output: true\n",
            d = dir.display()
        ),
    )
    .unwrap();

    eidolon()
        .arg("--log-dest")
        .arg(dir.join("run.log"))
        .args(["filter-reads", "-c"])
        .arg(&config)
        .assert()
        .success();

    // The writer compresses, so plain input gets a `.gz` output beside it.
    let reads = gz_lines(&dir.join("reads_filter.fastq.gz"));
    assert_eq!(reads.len(), 4, "exactly the in-region read: {reads:?}");
    assert!(reads[0].contains("chr1_0000001000_0000002000"), "{reads:?}");

    let calls = gz_lines(&dir.join("calls_filter.vcf.gz"));
    let records: Vec<&String> = calls.iter().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records.len(), 1, "exactly the in-region variant: {calls:?}");
    assert!(records[0].starts_with("chr1\t1500\t"), "{calls:?}");
}
