mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_with_bnd_variant_produces_chimeric_reads_in_fastq() {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_bnd.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )
        .unwrap();
        writeln!(
            f,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        )
        .unwrap();
        // BND variant: H1N1_HA:500 joined to H1N1_HA:1000
        writeln!(
            f,
            "H1N1_HA\t500\t.\tG\tG[H1N1_HA:1000[\t60\tPASS\t.\tGT\t1/1"
        )
        .unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "bnd_fastq");
    config.coverage = 100; // High coverage to ensure we get some chimeric reads
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq_1 = work.join("bnd_fastq_r1.fastq.gz");
    assert!(
        out_fastq_1.exists(),
        "output FASTQ 1 not produced at {:?}",
        out_fastq_1
    );

    // Check if any read name contains "chimeric"
    let fastq_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_fastq_1).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };

    let has_chimeric = fastq_lines.iter().any(|l| l.contains("EIDOLON_chimeric"));
    assert!(
        has_chimeric,
        "Expected to find chimeric reads in FASTQ output. Lines: {:?}",
        &fastq_lines[..fastq_lines.len().min(10)]
    );
}
