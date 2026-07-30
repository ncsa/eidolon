mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

#[test]
fn gen_reads_with_inv_variant_produces_chimeric_reads_in_fastq() {
    let (_dir, work) = fresh_workdir();

    let input_vcf = work.join("input_inv.vcf");
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
        // INV variant: H1N1_HA:500 to 1000
        writeln!(f, "H1N1_HA\t500\t.\tG\t<INV>\t60\tPASS\tEND=1000\tGT\t1/1").unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "inv_fastq");
    config.coverage = 100; // High coverage to ensure we get some chimeric reads
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_fastq_1 = work.join("inv_fastq_r1.fastq.gz");
    assert!(
        out_fastq_1.exists(),
        "output FASTQ 1 not produced at {:?}",
        out_fastq_1
    );

    // Check if any read name contains "chimeric_INV"
    let fastq_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_fastq_1).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };

    let has_chimeric = fastq_lines
        .iter()
        .any(|l| l.contains("EIDOLON_chimeric_INV"));
    assert!(
        has_chimeric,
        "Expected to find chimeric INV reads in FASTQ output. Lines: {:?}",
        &fastq_lines[..fastq_lines.len().min(10)]
    );

    // Check both junctions are represented. QNAMEs now end with
    //   EIDOLON_chimeric_INV_<contig>_<pos>_<end>_<junction>_<16-hex>/<mate>
    // (the 16-hex frag_idx tag was added for QNAME-uniqueness per #210).
    // The substring `_<junction>_0` is unambiguous for these synthetic
    // H1N1 fixtures — the hex frag_idx values all start with `0` since
    // num_frags is well under 2^60.
    let has_junction_1 = fastq_lines.iter().any(|l| l.contains("_1_0"));
    let has_junction_2 = fastq_lines.iter().any(|l| l.contains("_2_0"));

    assert!(has_junction_1, "Expected to find Junction 1 chimeric reads");
    assert!(has_junction_2, "Expected to find Junction 2 chimeric reads");
}
