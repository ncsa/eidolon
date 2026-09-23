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
//! These tests use an explicit high-quality sequencing-error model so their known-answer
//! assertions remain about junction geometry rather than the current global defaults.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};
use eidolon_core::models::{
    quality_scores::QualityScoreModel, sequencing_error_model::SequencingErrorModel,
};
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;

const CONTIG: &str = "H1N1_HA";
/// A correctly placed piece scores ~0.98 against its window; a wrong one ~0.25.
const MIN_IDENTITY: f64 = 0.90;

/// Configure a Phred-60 error model for geometry tests.  Its per-base error probability
/// is 1e-6, so the fixture does not depend on the production default error rates.
fn configure_geometry_error_model(config: &mut GenReadsConfig, work: &Path) {
    let quality_model =
        QualityScoreModel::from_counts(vec![60], 50, vec![1.0], vec![vec![vec![1.0]]; 50], false)
            .unwrap();
    let model = SequencingErrorModel::from_raw_data(0.0, quality_model, None).unwrap();
    let path = work.join("geometry-sequencing-error-model.json.gz");
    model.write_model(&path).unwrap();
    config.sequence_error_model = Some(path);
}

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
        // <INV>: POS is the base BEFORE the event, so the inverted block is 1-based
        // [POS+1, END] = 0-based [pos, end), carried reverse-complemented. TWO novel
        // junctions, at offsets `pos` and `end`.
        Sv::Inv { pos, end } => {
            out.extend_from_slice(&reference[..*pos]);
            out.extend(revcomp(&reference[*pos..*end]));
            out.extend_from_slice(&reference[*end..]);
        }
        // Tandem <DUP>: POS is the base BEFORE the event, so the duplicated bases
        // are 1-based [POS+1, END] = 0-based [pos, end), and the derived haplotype
        // carries that block twice. The novel junction is at offset `end`.
        Sv::Dup { pos, end } => {
            out.extend_from_slice(&reference[..*end]);
            out.extend_from_slice(&reference[*pos..*end]);
            out.extend_from_slice(&reference[*end..]);
        }
        // VCF 4.2 case 2, `t]p]`: REF[..=pos] + revcomp(MATE[..=mate_pos]).
        Sv::BndCase2 { pos, mate_pos } => {
            out.extend_from_slice(&reference[..*pos]);
            out.extend(revcomp(&reference[..*mate_pos]));
        }
        // VCF 4.2 case 1, `t[p[`: REF[..=pos] + MATE[mate_pos..]. A DIRECT join —
        // nothing reverse-complemented. This is what a deletion-like junction looks
        // like, and until geometry sampling existed eidolon could not emit one.
        Sv::BndCase1 { pos, mate_pos } => {
            out.extend_from_slice(&reference[..*pos]);
            out.extend_from_slice(&reference[*mate_pos - 1..]);
        }
        // VCF 4.2 case 4, `]p]t`: MATE[..=mate_pos] + REF[pos..]. Also direct, with the
        // mate piece LEADING. The reciprocal of case 1.
        Sv::BndCase4 { pos, mate_pos } => {
            out.extend_from_slice(&reference[..*mate_pos]);
            out.extend_from_slice(&reference[*pos - 1..]);
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
    Inv { pos: usize, end: usize },
    Dup { pos: usize, end: usize },
    BndCase2 { pos: usize, mate_pos: usize },
    BndCase3 { pos: usize, mate_pos: usize },
    BndCase1 { pos: usize, mate_pos: usize },
    BndCase4 { pos: usize, mate_pos: usize },
}

/// Contiguous match, or a match interrupted by ONE short indel.
///
/// The default sequencing-error model emits indels, and a purely positional
/// comparison frameshifts on them: 7 of 80 INV chimeric reads scored below 0.90 for
/// that reason alone, and all 7 were explained by a single 1-2bp indel against the
/// SAME derived haplotype. Tolerating one small gap keeps these tests measuring
/// GEOMETRY rather than the error model — without it, a correct implementation
/// looks ~9% broken and the noise floor hides real regressions.
fn matches_with_one_small_indel(read: &[u8], hay: &[u8]) -> bool {
    let find = |needle: &[u8]| -> Option<usize> {
        if needle.is_empty() || hay.len() < needle.len() {
            return None;
        }
        hay.windows(needle.len()).position(|w| w == needle)
    };
    for k in 8..read.len().saturating_sub(7) {
        let Some(p) = find(&read[..k]) else { continue };
        let rest = &read[k..];
        for g in -3i64..=3 {
            let q = (p + k) as i64 + g;
            if q < 0 || q as usize + rest.len() > hay.len() {
                continue;
            }
            if identity(rest, &hay[q as usize..q as usize + rest.len()]) >= 0.95 {
                return true;
            }
        }
    }
    false
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

/// Chimeric reads with their QNAMEs, for assertions about which junction a read came from.
fn chimeric_reads_with_names(tag: &str, record: &str) -> Vec<(String, Vec<u8>)> {
    let (_dir, work) = fresh_workdir();
    let input_vcf = work.join(format!("input_{tag}.vcf"));
    write_sv_vcf(&input_vcf, record);
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), tag);
    config.read_len = 50;
    // Label-routing assertions require independent spanning reads at both inversion
    // junctions. Use a denser deterministic fixture rather than weakening that proof.
    config.coverage = 60;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.input_vcf = Some(input_vcf);
    config.mutation_rate = Some(0.0);
    configure_geometry_error_model(&mut config, &work);
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
            out.push((lines[i].clone(), lines[i + 1].as_bytes().to_vec()));
        }
    }
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
    configure_geometry_error_model(&mut config, &work);
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
    let mut unmatched: Vec<String> = Vec::new();
    let mut exact_spanning = 0usize;
    for read in &reads {
        let (off, id) = best_contiguous(read, &derived);
        if id >= MIN_IDENTITY || matches_with_one_small_indel(read, &derived) {
            scored += 1;
        } else {
            // Not a geometry failure on its own: the sequencing-error model is active,
            // and a read carrying BOTH a substitution and an indel cannot be anchored
            // by an exact-prefix search. Measured on INV, 4 of 5 such reads were a
            // single indel and the fifth was a substitution plus a 1bp deletion — all
            // against this same haplotype. Collected and bounded below rather than
            // failed individually; a real geometry error puts EVERY read at ~0.62, far
            // under the bound.
            unmatched.push(format!(
                "identity {id:.2} at {off}: {}",
                String::from_utf8_lossy(read)
            ));
            continue;
        }

        // The negative control only means anything for reads that straddle the
        // breakpoint SUBSTANTIALLY. A read clipping the junction by 3 bases is still
        // ~94% identical to the reference — correctly so — and asserting otherwise
        // would fail on healthy output. Requiring a quarter of the read either side
        // is also the length BWA needs to split-align it (#224).
        let anchor = read.len() / 4;
        if off + anchor <= junction && off + read.len() >= junction + anchor {
            spanning += 1;
            // Indel tolerance is needed for sequencing errors, but it also absorbs a
            // 1bp geometry shift — reverting a junction offset by one passed this file
            // until this check existed. An error-free read spanning the junction must
            // match EXACTLY, and a shifted junction frameshifts every such read.
            if id >= 0.98 {
                exact_spanning += 1;
            }
            // THE assertion this file exists for: a read genuinely spanning the
            // junction cannot also be ordinary contiguous reference. Dropping a
            // reverse complement, inverting a junction dispatch, or collapsing a
            // junction back onto reference all surface here.
            let (roff, rid) = best_contiguous(read, reference);
            assert!(
                rid < MIN_IDENTITY && !matches_with_one_small_indel(read, reference),
                "{tag}: a read spanning the breakpoint by >={anchor}bp either side also \
                 matches the UNBROKEN reference at {roff} (identity {rid:.2}) — it \
                 carries no junction signal.\n  read: {}",
                String::from_utf8_lossy(read)
            );
        }
    }
    // Coverage of the harness's own input, reported rather than assumed. A geometry
    // defect fails essentially every read, so this bound is not a fudge factor: the
    // mutations this file guards against land at ~0.62 identity across the board.
    let matched_pct = 100 * scored / reads.len();
    assert!(
        matched_pct >= 90,
        "{tag}: only {scored} of {} chimeric reads match the derived haplotype \
         ({matched_pct}%). Sequencing-error noise accounts for a few percent; this is \
         a geometry failure.\n  unmatched:\n    {}",
        reads.len(),
        unmatched.join("\n    ")
    );
    assert!(
        exact_spanning >= 3,
        "{tag}: only {exact_spanning} read(s) match the junction EXACTLY (of {spanning} \
         spanning it). Error-free reads crossing a correctly-built junction match at \
         1.00; a junction shifted by even one base leaves none of them exact."
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

/// A tandem `<DUP>` must duplicate 1-based [POS+1, END], not [POS, END].
///
/// `get_dup_pieces` treated `location` as the FIRST AFFECTED base while
/// `get_del_pieces` treated it as the last unaffected one, so an input <DUP> and an
/// input <DEL> in the same file meant different things by POS and the duplicated
/// block started one base early. Nothing tested it: `dup_chimeric.rs` asserts QNAMEs,
/// and `get_dup_pieces` had no unit tests at all.
#[test]
fn dup_chimeric_reads_duplicate_the_vcf_conventional_block() {
    assert_derived(
        "chim_dup",
        "H1N1_HA\t600\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=1400\tGT\t1/1",
        Sv::Dup {
            pos: 600,
            end: 1400,
        },
        // The novel junction sits where the duplicated block is re-entered.
        //
        // The block is 800bp against a ~200bp fragment. That was once believed to matter:
        // a fragment spanning the WHOLE duplication was thought to need a three-piece
        // stitch get_dup_pieces cannot express, costing "4 of 30" reads. MEASURED FALSE —
        // see short_dup_spanned_by_a_fragment_still_matches_the_derived_haplotype below,
        // which sweeps 60..800bp and passes throughout. The 13% was almost certainly the
        // sequencing-error noise that matches_with_one_small_indel was later added to
        // absorb. The old note also cited #474, which is closed and about anchor-base
        // conventions, not stitching.
        1400,
    );
}

/// A DUP block SHORTER than a fragment — the "three-piece stitch" case, which turns out NOT
/// to be broken.
///
/// The test above uses an 800bp block against a ~200bp fragment deliberately, on the stated
/// grounds that a fragment spanning the WHOLE duplication "needs a three-piece stitch (left +
/// block + right) that get_dup_pieces cannot express — it returns two — and 4 of 30 reads then
/// match no haplotype at all". `docs/sv_support_matrix.md` carried that as an open defect.
///
/// MEASURED, and it does not reproduce. Blocks of 60/100/150/200/400/800bp against a ~200bp
/// fragment — every one of which a fragment can span — all match the derived haplotype at the
/// >=90% threshold the rest of this file uses.
///
/// The likely explanation is in this file's own header: `matches_with_one_small_indel` exists
/// because a purely positional comparison frameshifts on sequencing-error indels, and "without
/// it, a correct implementation looks ~9% broken". The claimed defect was 4 of 30 = 13%, inside
/// that band — so the original measurement probably predates that tolerance and was the error
/// model rather than a geometry failure.
///
/// WHAT THIS DOES NOT SHOW: that no fragment is SKIPPED. It asserts the reads that exist are
/// correct, not that every read that should exist does. If fragments needing three pieces were
/// silently dropped rather than mis-stitched, coverage across a short DUP would dip and this
/// test would stay green. That needs a depth comparison against a no-variant control — the
/// method `sv_support_matrix.rs` uses — before the cell is declared fully clean.
#[test]
fn short_dup_spanned_by_a_fragment_still_matches_the_derived_haplotype() {
    assert_derived(
        "chim_dup_short",
        "H1N1_HA\t600\t.\tG\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=750\tGT\t1/1",
        Sv::Dup { pos: 600, end: 750 },
        750,
    );
}

/// `<INV>` inverts 1-based [POS+1, END] and must reverse-complement it.
///
/// `get_inv_pieces` had NO tests: four separate mutations passed the whole suite,
/// including inverting the junction dispatch (`junction == 1` -> `== 2`) and turning
/// off either reverse-complement flag. With the revcomp off, the "inversion" junction
/// carries no inversion signal at all, and `inv_fastq.rs` stayed green because it only
/// greps QNAMEs for `_1_0`.
#[test]
fn inv_chimeric_reads_reverse_complement_the_inverted_block() {
    assert_derived(
        "chim_inv",
        "H1N1_HA\t600\t.\tG\t<INV>\t60\tPASS\tSVTYPE=INV;END=1400\tGT\t1/1",
        Sv::Inv {
            pos: 600,
            end: 1400,
        },
        600,
    );
}

/// The junction a read is LABELLED with must be the junction it actually spans.
///
/// `assert_derived` only asks whether a read is a slice of the derived haplotype, and
/// both INV junctions produce valid slices — so swapping the junction dispatch
/// (`junction == 1` -> `== 2`) passes it. `inv_fastq.rs` greps QNAMEs for `_1_0`/`_2_0`
/// without checking where those reads came from, so nothing caught it anywhere.
#[test]
fn inv_junction_labels_match_the_junction_the_read_spans() {
    let reference = &load_reference()[CONTIG];
    let (pos, end) = (600usize, 1400usize);
    let derived = derived_haplotype(reference, &Sv::Inv { pos, end });
    let reads = chimeric_reads_with_names(
        "chim_inv_lbl",
        "H1N1_HA\t600\t.\tG\t<INV>\t60\tPASS\tSVTYPE=INV;END=1400\tGT\t1/1",
    );
    let (mut j1, mut j2) = (0usize, 0usize);
    for (name, read) in &reads {
        let (off, id) = best_contiguous(read, &derived);
        if id < 0.98 {
            continue; // sequencing-error noise; geometry is asserted elsewhere
        }
        let anchor = read.len() / 4;
        let spans = |j: usize| off + anchor <= j && off + read.len() >= j + anchor;
        if name.contains("_1_0") && spans(pos) {
            j1 += 1;
        }
        if name.contains("_2_0") && spans(end) {
            j2 += 1;
        }
        // A read labelled for one junction must never be found spanning the other.
        assert!(
            !(name.contains("_1_0") && spans(end)),
            "read labelled junction 1 spans junction 2 at {end}: {name}"
        );
        assert!(
            !(name.contains("_2_0") && spans(pos)),
            "read labelled junction 2 spans junction 1 at {pos}: {name}"
        );
    }
    assert!(
        j1 >= 3 && j2 >= 3,
        "expected reads spanning both labelled junctions; got j1={j1}, j2={j2}"
    );
}

/// The two DIRECT forms, which geometry sampling (#458 item 2) newly makes reachable.
///
/// Until then eidolon emitted only `t]p]`, so cases 1 and 4 could not appear in output
/// at all and nothing tested them. They matter because a direct join carries no
/// inversion signature: Manta classifies these as DEL/DUP:TANDEM rather than putting
/// them in its BND bucket, which is precisely the separation this change exists to
/// restore.
#[test]
fn bnd_case1_direct_join_carries_no_reverse_complement() {
    assert_derived(
        "chim_bnd1",
        "H1N1_HA\t500\t.\tG\tG[H1N1_HA:1500[\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
        Sv::BndCase1 {
            pos: 500,
            mate_pos: 1500,
        },
        500,
    );
}

#[test]
fn bnd_case4_direct_join_leads_with_the_mate_piece() {
    assert_derived(
        "chim_bnd4",
        "H1N1_HA\t500\t.\tG\t]H1N1_HA:1500]G\t60\tPASS\tSVTYPE=BND\tGT\t1/1",
        Sv::BndCase4 {
            pos: 500,
            mate_pos: 1500,
        },
        1500,
    );
}
