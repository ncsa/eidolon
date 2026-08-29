//! A breakend ALT may carry novel sequence inserted at the junction (VCF 4.2 §5.4). The
//! reads must contain it (#498).
//!
//! WHY THIS FILE EXISTS: the truth VCF round-trips the full ALT verbatim, so every
//! VCF-level check passes whether or not the insert reached a single read. #498 measured 20
//! chimeric reads at a junction and not one carrying the declared insert — a benchmark built
//! from that data scores a caller wrong for correctly reporting no inserted sequence, and
//! right for hallucinating it. Only a test that reads the SEQUENCE can see this, which is
//! the same lesson as #451 and #516.

mod common;
use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use std::io::Write as _;

/// Inserted at the junction. Chosen to be absent from the H1N1 fixture so a hit cannot come
/// from the reference: 24 bases is ~4^-24, and it is asserted below rather than assumed.
const INSERT: &str = "GATTACAGATTACAGGCCTTAAGC";

fn write_bnd_vcf(path: &std::path::Path, alt_left: &str, alt_right: &str) {
    let mut f = std::fs::File::create(path).unwrap();
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
    // REF bases are the ACTUAL reference at those positions. The reads are stitched from
    // the reference, not from the REF column, so a wrong REF here would make every probe
    // below miss for a reason that has nothing to do with the insert.
    writeln!(
        f,
        "H1N1_HA\t600\t.\tT\t{alt_left}\t60\tPASS\tSVTYPE=BND\tGT\t1/1"
    )
    .unwrap();
    writeln!(
        f,
        "H1N1_PB2\t900\t.\tG\t{alt_right}\t60\tPASS\tSVTYPE=BND\tGT\t1/1"
    )
    .unwrap();
}

fn run(work: &std::path::Path, name: &str, input_vcf: std::path::PathBuf) -> Vec<String> {
    let mut config = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), name);
    config.coverage = 100;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    config.sv_rate_scale = Some(0.0); // only the supplied junction, no de novo SVs
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    use flate2::read::MultiGzDecoder;
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(work.join(format!("{name}_r1.fastq.gz"))).unwrap();
    BufReader::new(MultiGzDecoder::new(f))
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

/// Read one contig from the bundled H1N1 FASTA. Probes are derived from the reference
/// rather than hardcoded so they cannot silently rot if the fixture changes.
fn contig(name: &str) -> String {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut out = String::new();
    let mut in_contig = false;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            in_contig = h.split_whitespace().next() == Some(name);
        } else if in_contig {
            out.push_str(line.trim());
        }
    }
    assert!(!out.is_empty(), "contig {name} not found in the fixture");
    out
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            o => o,
        })
        .collect()
}

fn count_hits(lines: &[String], probe: &str) -> usize {
    let rc = revcomp(probe);
    lines
        .iter()
        .filter(|l| !l.starts_with('@') && !l.starts_with('+'))
        .filter(|l| l.contains(probe) || l.contains(&rc))
        .count()
}

#[test]
fn bnd_inserted_sequence_reaches_the_reads() {
    let (_dir, work) = fresh_workdir();
    let vcf = work.join("bnd_insert.vcf");
    // t[p[ form: `t` starts with the REF base, so the insert follows it.
    // ]p]t form: `t` ends with the REF base, so the insert precedes it.
    write_bnd_vcf(
        &vcf,
        &format!("A{INSERT}[H1N1_PB2:900["),
        &format!("]H1N1_HA:600]{INSERT}C"),
    );
    let lines = run(&work, "bnd_ins", vcf);

    let chimeric = lines
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric"))
        .count();
    assert!(
        chimeric > 0,
        "no chimeric reads generated — the junction itself did not fire, so this test \
         cannot say anything about the insert"
    );

    let hits = count_hits(&lines, INSERT);
    assert!(
        hits > 0,
        "{chimeric} chimeric read(s) at the junction and NOT ONE carries the {} bp \
         inserted sequence the ALT declares. The truth VCF asserts an insertion the reads \
         do not contain (#498).",
        INSERT.len()
    );
}

/// The insert must be spliced BETWEEN the reference pieces, not appended or prepended.
///
/// A probe spanning [tail of the local reference piece | head of the insert] can only match
/// if the two are adjacent in the read. Counting the insert alone would pass even if it were
/// dumped at the end of the fragment.
#[test]
fn inserted_sequence_is_adjacent_to_the_reference_piece() {
    let (_dir, work) = fresh_workdir();
    let vcf = work.join("bnd_adj.vcf");
    write_bnd_vcf(
        &vcf,
        &format!("T{INSERT}[H1N1_PB2:900["),
        &format!("]H1N1_HA:600]{INSERT}G"),
    );
    let lines = run(&work, "bnd_adj", vcf);

    // Last 12 reference bases before the junction, joined to the first 12 inserted bases.
    // This can only match if the two are adjacent in the read — counting the insert alone
    // would pass even if it were dumped at the end of the fragment.
    let ha = contig("H1N1_HA");
    let probe = format!("{}{}", &ha[588..600], &INSERT[..12]);
    assert!(
        count_hits(&lines, &probe) > 0,
        "the inserted sequence is present but NOT adjacent to the reference piece — it is \
         being spliced in the wrong position"
    );

    // And the bare junction those 24 bases would form WITHOUT an insert must now be gone:
    // the insert splits them. This is what distinguishes "spliced between" from "appended".
    let pb2 = contig("H1N1_PB2");
    let bare = format!("{}{}", &ha[588..600], &pb2[899..911]);
    assert_eq!(
        count_hits(&lines, &bare),
        0,
        "reads still contain the two reference pieces flush against each other, so the \
         insert is being added somewhere other than the junction"
    );
}

/// MUST NOT FIRE: a breakend with no inserted sequence must leave the two reference pieces
/// flush. If the anchor base were ever mistaken for an insert, it would be duplicated at
/// every junction in every run — a one-base corruption invisible to any count.
#[test]
fn a_bare_breakend_inserts_nothing() {
    let (_dir, work) = fresh_workdir();
    let vcf = work.join("bnd_bare.vcf");
    write_bnd_vcf(&vcf, "T[H1N1_PB2:900[", "]H1N1_HA:600]G");
    let lines = run(&work, "bnd_bare", vcf);

    let chimeric = lines
        .iter()
        .filter(|l| l.contains("EIDOLON_chimeric"))
        .count();
    assert!(chimeric > 0, "no chimeric reads to inspect");

    assert_eq!(
        count_hits(&lines, INSERT),
        0,
        "a bare breakend produced reads containing sequence no record declares"
    );

    // The two reference pieces must be FLUSH. If the anchor were ever mistaken for an
    // insert it would be duplicated here, breaking this exact 24-mer — a one-base
    // corruption at every junction in every run that no count could see.
    let ha = contig("H1N1_HA");
    let pb2 = contig("H1N1_PB2");
    let bare = format!("{}{}", &ha[588..600], &pb2[899..911]);
    assert!(
        count_hits(&lines, &bare) > 0,
        "a bare breakend did not put the two reference pieces flush at the junction"
    );
}
