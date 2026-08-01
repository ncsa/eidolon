//! Known-answer tests on the SEQUENCE of chimeric junction reads.
//!
//! Why this file exists: an audit found that **no test anywhere in the repo read a
//! chimeric read's bases**. Every chimeric test asserts on QNAMEs and counts, so 11 of 11
//! sequence-affecting mutations passed the full suite. The two that matter most:
//!
//!   * `get_stitched_sequence` could be changed to ignore BOTH reverse-complement flags
//!     (`if false && rev1`) with nothing failing — including the unit tests that pin
//!     `get_bnd_pieces`' output. Those pin the flag's VALUE; nothing checked it was
//!     HONOURED. The chain is: sv_model sets the flags -> get_bnd_pieces returns the
//!     layout -> get_stitched_sequence applies the revcomp, and only the last link
//!     actually produces bases.
//!   * `del_chimeric.rs` claims its QNAME check "pins the coordinate encoding so a future
//!     refactor that shifts position or end by ±1 fails the test". It cannot: the QNAME is
//!     formatted from the same `location`/`end` INPUTS, never from the emitted geometry.
//!     Both ±1 mutations passed. `dup_chimeric.rs` carries the identical false claim.
//!
//! The expectations below are derived from VCF 4.2 semantics, not read off the
//! implementation, so this is a known-answer test rather than a blessed baseline.
//!
//! Identity is compared with a tolerance because the default sequencing-error model is
//! active. That costs nothing in discrimination: a correctly stitched piece matches its
//! expected reference window at ~98%, while a piece in the wrong orientation or at the
//! wrong locus scores ~25% (random over 4 bases). The threshold sits far from both.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use std::collections::HashMap;
use std::io::Write as _;

const CONTIG: &str = "H1N1_HA";
/// A correctly placed piece scores ~0.98 against its window; a wrong one ~0.25.
const MIN_IDENTITY: f64 = 0.90;

fn load_reference() -> HashMap<String, Vec<u8>> {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut out: HashMap<String, Vec<u8>> = HashMap::new();
    let mut name = String::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            name = h.split_whitespace().next().unwrap().to_string();
            out.insert(name.clone(), Vec::new());
        } else {
            out.get_mut(&name).unwrap().extend(line.trim().bytes());
        }
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn identity(a: &[u8], b: &[u8]) -> f64 {
    if a.is_empty() || a.len() != b.len() {
        return 0.0;
    }
    let same = a
        .iter()
        .zip(b)
        .filter(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
        .count();
    same as f64 / a.len() as f64
}

/// Build the DERIVED haplotype the SV implies, from VCF 4.2 semantics.
///
/// This is the criterion the whole file rests on: every read emitted for a chimeric
/// fragment — whether it straddles the breakpoint (a split read) or sits wholly inside
/// one piece (the discordant-pair mate) — must be a contiguous slice of this sequence.
/// One rule, both read classes, and it is derived from what the SV MEANS rather than
/// from what the code does.
fn derived_haplotype(reference: &[u8], sv: &Sv) -> Vec<u8> {
    let mut out = Vec::new();
    match sv {
        // POS=p (1-based) keeps the anchor; bases p+1..=end are removed.
        Sv::Del { pos, end } => {
            out.extend_from_slice(&reference[..*pos]);
            out.extend_from_slice(&reference[*end..]);
        }
        // VCF 4.2 case 2, `t]p]`: REF[..=pos] + revcomp(MATE[..=mate_pos]).
        Sv::BndCase2 { pos, mate_pos } => {
            out.extend_from_slice(&reference[..*pos]);
            out.extend(revcomp(&reference[..*mate_pos]));
        }
        // VCF 4.2 case 3, `[p[t`: revcomp(MATE[mate_pos..]) + REF[pos..]. Here it is the
        // FIRST piece that is reverse-complemented, which case 2 never exercises.
        Sv::BndCase3 { pos, mate_pos } => {
            out.extend(revcomp(&reference[*mate_pos - 1..]));
            out.extend_from_slice(&reference[*pos - 1..]);
        }
    }
    out
}

enum Sv {
    Del { pos: usize, end: usize },
    BndCase2 { pos: usize, mate_pos: usize },
    BndCase3 { pos: usize, mate_pos: usize },
}

/// Best contiguous match of `read` anywhere in `hay`, returning `(offset, identity)`.
fn best_contiguous(read: &[u8], hay: &[u8]) -> (usize, f64) {
    let mut best = (0usize, 0.0f64);
    if hay.len() < read.len() {
        return best;
    }
    for offset in 0..=hay.len() - read.len() {
        let id = identity(read, &hay[offset..offset + read.len()]);
        if id > best.1 {
            best = (offset, id);
        }
    }
    best
}

fn write_sv_vcf(path: &std::path::Path, record: &str) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        f,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )
    .unwrap();
    writeln!(
        f,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"
    )
    .unwrap();
    writeln!(
        f,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">"
    )
    .unwrap();
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
    )
    .unwrap();
    writeln!(f, "{record}").unwrap();
}

/// Every chimeric read paired by QNAME, as `(r1, r2)`.
fn chimeric_pairs(tag: &str, record: &str) -> Vec<(Vec<u8>, Vec<u8>)> {
    let (_dir, work) = fresh_workdir();
    let input_vcf = work.join(format!("input_{tag}.vcf"));
    write_sv_vcf(&input_vcf, record);
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), tag);
    config.read_len = 50;
    config.coverage = 30;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    config.mutation_rate = Some(0.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let grab = |suffix: &str| -> HashMap<String, Vec<u8>> {
        let lines = read_gzip_fastq_lines(&work.join(format!("{tag}_{suffix}.fastq.gz")));
        let mut m = HashMap::new();
        for i in (0..lines.len()).step_by(4) {
            if lines[i].starts_with('@') && lines[i].contains("EIDOLON_chimeric") {
                let qname = lines[i].split('/').next().unwrap().to_string();
                m.insert(qname, lines[i + 1].as_bytes().to_vec());
            }
        }
        m
    };
    let (r1, r2) = (grab("r1"), grab("r2"));
    let mut out: Vec<(Vec<u8>, Vec<u8>)> = r1
        .iter()
        .filter_map(|(q, a)| r2.get(q).map(|b| (a.clone(), b.clone())))
        .collect();
    out.sort();
    out
}

/// Run gen-reads over one homozygous SV and return every chimeric read's sequence.
fn chimeric_reads(tag: &str, record: &str) -> Vec<Vec<u8>> {
    let (_dir, work) = fresh_workdir();
    let input_vcf = work.join(format!("input_{tag}.vcf"));
    write_sv_vcf(&input_vcf, record);

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), tag);
    config.read_len = 50;
    config.coverage = 30;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    // No de novo mutations: any mismatch against the derived haplotype must come from
    // the junction, not from a planted SNP.
    config.mutation_rate = Some(0.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let lines = read_gzip_fastq_lines(&work.join(format!("{tag}_r1.fastq.gz")));
    let mut out = Vec::new();
    for i in (0..lines.len()).step_by(4) {
        if lines[i].starts_with('@') && lines[i].contains("EIDOLON_chimeric") {
            out.push(lines[i + 1].as_bytes().to_vec());
        }
    }
    out
}

/// Every chimeric read must be a contiguous slice of the derived haplotype, and enough
/// of them must actually straddle the breakpoint to prove a junction was built at all.
///
/// `junction` is the offset in the derived haplotype where the two pieces meet. A read
/// covering it is a split read; one that does not is the discordant-pair mate, and both
/// are legitimate output of a chimeric fragment.
fn assert_derived(tag: &str, record: &str, sv: Sv, junction: usize) {
    let reference = &load_reference()[CONTIG];
    let derived = derived_haplotype(reference, &sv);
    let reads = chimeric_reads(tag, record);
    assert!(
        reads.len() >= 10,
        "{tag}: only {} chimeric read(s) — too few to conclude anything",
        reads.len()
    );

    let (mut scored, mut spanning) = (0usize, 0usize);
    for read in &reads {
        let (off, id) = best_contiguous(read, &derived);
        assert!(
            id >= MIN_IDENTITY,
            "{tag}: chimeric read is not a slice of the derived haplotype \
             (best identity {id:.2} at {off}). The reads do not encode the \
             rearrangement the VCF declares.\n  read: {}",
            String::from_utf8_lossy(read)
        );
        scored += 1;

        // The negative control only means anything for reads that straddle the
        // breakpoint SUBSTANTIALLY. A read clipping the junction by 3 bases is still
        // ~94% identical to the reference — correctly so — and asserting otherwise
        // would fail on healthy output. Requiring a quarter of the read either side
        // is also the length BWA needs to split-align it (#224).
        let anchor = read.len() / 4;
        if off + anchor <= junction && off + read.len() >= junction + anchor {
            spanning += 1;
            // THE assertion this file exists for: a read genuinely spanning the
            // junction cannot also be ordinary contiguous reference. Dropping a
            // reverse complement, inverting a junction dispatch, or collapsing a
            // junction back onto reference all surface here.
            let (roff, rid) = best_contiguous(read, reference);
            assert!(
                rid < MIN_IDENTITY,
                "{tag}: a read spanning the breakpoint by >={anchor}bp either side also \
                 matches the UNBROKEN reference at {roff} (identity {rid:.2}) — it \
                 carries no junction signal.\n  read: {}",
                String::from_utf8_lossy(read)
            );
        }
    }
    assert_eq!(
        scored,
        reads.len(),
        "{tag}: scored {scored} of {} chimeric reads",
        reads.len()
    );
    assert!(
        spanning >= 3,
        "{tag}: only {spanning} of {} chimeric reads straddle the breakpoint — \
         without split reads there is no junction signal to find",
        reads.len()
    );
}

/// `<DEL>` POS=500 END=800: bases 501..=800 removed, anchor at 500 kept, so the derived
/// haplotype is `REF[..500] + REF[800..]` and the junction sits at offset 500.
#[test]
fn del_chimeric_reads_match_the_derived_deletion_haplotype() {
    assert_derived(
        "chim_del",
        "H1N1_HA\t500\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=800\tGT\t1/1",
        Sv::Del { pos: 500, end: 800 },
        500,
    );
}

/// BND `G]H1N1_HA:1500]` at POS=500 is VCF 4.2 case 2: `REF[..=pos] +
/// revcomp(MATE[..=mate_pos])`. The mate piece is REVERSE-COMPLEMENTED, and until now
/// nothing verified that flag was applied — only that it was set.
#[test]
fn bnd_chimeric_reads_reverse_complement_the_mate_piece() {
    assert_derived(
        "chim_bnd",
        "H1N1_HA\t500\t.\tG\tG]H1N1_HA:1500]\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
        Sv::BndCase2 {
            pos: 500,
            mate_pos: 1500,
        },
        500,
    );
}

/// BND `[H1N1_HA:1500[G` at POS=500 is VCF 4.2 case 3: `revcomp(MATE[mate_pos..]) +
/// REF[pos..]`. Included specifically because the reverse complement falls on the FIRST
/// piece — case 2 leaves that code path untested, and mutating it was silently safe.
#[test]
fn bnd_case3_chimeric_reads_reverse_complement_the_leading_piece() {
    assert_derived(
        "chim_bnd3",
        "H1N1_HA\t500\t.\tG\t[H1N1_HA:1500[G\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
        Sv::BndCase3 {
            pos: 500,
            mate_pos: 1500,
        },
        // revcomp(REF[1499..]) is 1701-1499 = 202 bases, so the junction sits at 202.
        202,
    );
}

/// R1 and R2 must be sequenced from OPPOSITE strands of the same fragment (FR). Removing
/// `reverse_complement` from the four `generate_*_pair` sites passed every existing test,
/// including `fastq_validation.rs::paired_ended_fastq_pair_invariants` — both mates would
/// be same-strand, so there would be no proper pairs and no discordant-pair signal for a
/// caller to find, which is half of what an SV simulator exists to produce.
///
/// The discriminating assertion is the negative one: R2 must NOT match the derived
/// haplotype in the forward orientation. Observed on healthy output: 0/30 forward,
/// 27/30 reverse-complement.
#[test]
fn chimeric_mates_are_sequenced_from_opposite_strands() {
    let reference = &load_reference()[CONTIG];
    let derived = derived_haplotype(reference, &Sv::Del { pos: 500, end: 800 });
    let pairs = chimeric_pairs(
        "chim_pair",
        "H1N1_HA\t500\t.\tG\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=800\tGT\t1/1",
    );
    assert!(
        pairs.len() >= 10,
        "only {} chimeric pair(s) — too few to conclude anything",
        pairs.len()
    );

    let (mut r2_forward, mut r2_reversed) = (0usize, 0usize);
    for (r1, r2) in &pairs {
        assert!(
            best_contiguous(r1, &derived).1 >= MIN_IDENTITY,
            "R1 is not a slice of the derived haplotype: {}",
            String::from_utf8_lossy(r1)
        );
        if best_contiguous(r2, &derived).1 >= MIN_IDENTITY {
            r2_forward += 1;
        }
        if best_contiguous(&revcomp(r2), &derived).1 >= MIN_IDENTITY {
            r2_reversed += 1;
        }
    }
    assert_eq!(
        r2_forward,
        0,
        "{r2_forward} of {} R2 reads match the derived haplotype FORWARD — the mates are \
         same-strand, so the pairs are not FR and carry no discordant-pair signal",
        pairs.len()
    );
    assert!(
        r2_reversed * 10 >= pairs.len() * 8,
        "only {r2_reversed} of {} R2 reads reverse-complement-match the derived haplotype",
        pairs.len()
    );
}
