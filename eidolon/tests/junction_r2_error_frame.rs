//! #777 — SV junction R2 reads must draw substitution errors in the same frame as every other R2.
//!
//! CORRECTNESS CRITERION. The substitution matrix is defined in REFERENCE orientation: the
//! ordinary paired writer generates R2 forward over the fragment, draws its errors, and only then
//! reverse-complements the record, and the BAM fit (`transition_matrix_round_trip.rs`) counts
//! mismatches the same way. So on ANY R2 read, a mismatch against the molecule's top strand at a
//! top-strand base `b` must be `matrix[b]`, whether the read came from the ordinary writer or
//! from an SV junction writer. It is FALSIFIED by a junction R2 whose mismatches follow the
//! matrix's reverse complement.
//!
//! THE PLANTED MATRIX makes the two frames disjoint: A->C, C->G, G->T, T->A, each with
//! probability 1. Its reverse complement is A->T, C->A, G->C, T->G, which shares no cell with it,
//! so the in-frame fraction is 1.0 when the frame is right and 0.0 when it is wrong. Nothing in
//! between is possible, which is what lets a few hundred errors decide it.
//!
//! THE KNOWN ANSWER is independent of the code: each junction read is placed on the derived
//! haplotype the SV implies under VCF 4.2 semantics (built here from the reference, not from
//! eidolon's pieces), and ordinary R2 reads are read against the reference from the golden BAM.
//! Indel errors are switched off in the model, so an ungapped placement is exact.
//!
//! QUALITY ORIENTATION rides along (#734): the quality model is Q40 for cycles 1-50 and Q7 for
//! cycles 51-100, so a junction R2 must carry `I` then `(` in the FASTQ and its mismatches must
//! sit in its LAST 50 cycles. Errors are drawn from each base's own quality, so a fix that
//! flipped the read without laying R2's quality down backwards first would fail both.
//!
//! MUST NOT FIRE. Junction R1 and ordinary R2 are measured in the same run and must be in frame
//! too. Every stratum's denominator is asserted.
//!
//! Fixture: ecoli (4.6 Mb, one contig), per CLAUDE.md. BND is planted intra-contig from the
//! input VCF, which reaches the BND writer without needing a second contig; case 1 (`t[p[`,
//! forward second piece) and case 2 (`t]p]`, reverse-complemented second piece) are both
//! planted because they end R2's window on opposite strands of the reference.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, read_gzip_fastq_lines};
use eidolon_core::models::{
    quality_scores::QualityScoreModel, sequencing_error_model::SequencingErrorModel,
};
use eidolon_core::structs::transition_matrix::TransitionMatrix;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

const READ_LEN: usize = 100;
/// Cycles 1..=HALF are Q40, the rest Q7.
const HALF: usize = 50;
const Q_HIGH: usize = 40;
const Q_LOW: usize = 7;
/// Events per SV class.
const PER_TYPE: usize = 10;
/// Flank of the derived haplotype kept either side of an event.
const FLANK: usize = 1_500;
const SEED: usize = 16;

fn ecoli() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/ecoli.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
}

fn load_reference() -> (String, Vec<u8>) {
    let text = std::fs::read_to_string(ecoli()).unwrap();
    let mut name = String::new();
    let mut seq = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            name = h.split_whitespace().next().unwrap().to_string();
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    (name, seq)
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

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The planted substitution: A->C, C->G, G->T, T->A.
fn planted(from: usize) -> usize {
    (from + 1) % 4
}

/// The planted matrix's reverse complement, which is what a read-orientation draw produces.
fn reverse_complement_frame(from: usize) -> usize {
    3 - planted(3 - from)
}

/// The planted matrix, Q40 then Q7 by cycle. `indel_probability` is the chance an error is an
/// indel, and `insertion_fraction` the share of those that are insertions.
fn write_model(work: &Path, indel_probability: f64, insertion_fraction: f64) -> PathBuf {
    let options = vec![Q_LOW, Q_HIGH];
    let trans: Vec<Vec<Vec<f64>>> = (0..READ_LEN)
        .map(|p| {
            // distros[p] draws the score at index p + 1.
            let row = if p + 1 < HALF {
                vec![0.0, 1.0]
            } else {
                vec![1.0, 0.0]
            };
            vec![row.clone(), row]
        })
        .collect();
    let quality =
        QualityScoreModel::from_counts(options, READ_LEN, vec![0.0, 1.0], trans, false).unwrap();
    let mut rows = [[0.0f64; 4]; 4];
    for (from, row) in rows.iter_mut().enumerate() {
        row[planted(from)] = 1.0;
    }
    let matrix = TransitionMatrix::from(rows[0], rows[1], rows[2], rows[3]).unwrap();
    let model = SequencingErrorModel::from_raw_data(0.0, quality, Some(matrix)).unwrap();
    let raw = work.join("raw_model.json.gz");
    model.write_model(&raw).unwrap();

    // `from_raw_data` fixes indel_probability at 0.01 and insertion_fraction at 0.4.
    let mut v: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(
        std::fs::File::open(&raw).unwrap(),
    ))
    .unwrap();
    for key in ["indel_probability", "insertion_fraction"] {
        assert!(v.get(key).is_some(), "model schema changed: no {key} field");
    }
    v["indel_probability"] = serde_json::json!(indel_probability);
    v["insertion_fraction"] = serde_json::json!(insertion_fraction);
    let path = work.join("planted_model.json.gz");
    let mut enc = flate2::write::GzEncoder::new(
        std::fs::File::create(&path).unwrap(),
        flate2::Compression::default(),
    );
    enc.write_all(serde_json::to_string(&v).unwrap().as_bytes())
        .unwrap();
    enc.finish().unwrap();
    // The patched file must load and carry the patch.
    let back = SequencingErrorModel::from_file(&path).unwrap();
    assert_eq!(
        back.quality_score_model().quality_score_options,
        [Q_LOW, Q_HIGH]
    );
    path
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Class {
    Del,
    Dup,
    Inv,
    Bnd,
}

struct Event {
    class: Class,
    /// QNAME prefix the writer gives this event's reads.
    qname_prefix: String,
    vcf: String,
    /// The derived haplotype around the event, top strand.
    haplotype: Vec<u8>,
}

/// Plant PER_TYPE of each class, 80 kb apart. Derived haplotypes follow VCF 4.2 exactly as
/// `chimeric_sequence_content.rs::derived_haplotype` does, windowed to FLANK.
fn events(contig: &str, reference: &[u8]) -> Vec<Event> {
    let base_at = |p1: usize| reference[p1 - 1] as char;
    let mut out = Vec::new();
    for t in 0..5 {
        for i in 0..PER_TYPE {
            let pos = 200_000 + (t * PER_TYPE + i) * 80_000; // 1-based POS
            let r = base_at(pos);
            let ev = match t {
                0 => {
                    let end = pos + 2_000;
                    Event {
                        class: Class::Del,
                        qname_prefix: format!("EIDOLON_chimeric_DEL_{contig}_{pos}_{end}_"),
                        vcf: format!(
                            "{contig}\t{pos}\t.\t{r}\t<DEL>\t60\tPASS\tSVTYPE=DEL;END={end}\tGT\t1/1"
                        ),
                        haplotype: [&reference[pos - FLANK..pos], &reference[end..end + FLANK]]
                            .concat(),
                    }
                }
                1 => {
                    let end = pos + 3_000;
                    Event {
                        class: Class::Dup,
                        qname_prefix: format!("EIDOLON_chimeric_DUP_{contig}_{pos}_{end}_"),
                        vcf: format!(
                            "{contig}\t{pos}\t.\t{r}\t<DUP>\t60\tPASS\tSVTYPE=DUP;END={end}\tGT\t1/1"
                        ),
                        // ref[..end] + ref[pos..end] + ref[end..]
                        haplotype: [
                            &reference[end - FLANK..end],
                            &reference[pos..end],
                            &reference[end..end + FLANK],
                        ]
                        .concat(),
                    }
                }
                2 => {
                    let end = pos + 3_000;
                    Event {
                        class: Class::Inv,
                        qname_prefix: format!("EIDOLON_chimeric_INV_{contig}_{pos}_{end}_"),
                        vcf: format!(
                            "{contig}\t{pos}\t.\t{r}\t<INV>\t60\tPASS\tSVTYPE=INV;END={end}\tGT\t1/1"
                        ),
                        // ref[..pos] + rc(ref[pos..end]) + ref[end..]
                        haplotype: [
                            reference[pos - FLANK..pos].to_vec(),
                            revcomp(&reference[pos..end]),
                            reference[end..end + FLANK].to_vec(),
                        ]
                        .concat(),
                    }
                }
                3 => {
                    // Case 1, t[p[: REF[..=pos] + MATE[mate_pos..]
                    let mate = pos + 30_000;
                    Event {
                        class: Class::Bnd,
                        qname_prefix: format!(
                            "EIDOLON_chimeric_{contig}_{}_{contig}_{mate}_",
                            pos - 1
                        ),
                        vcf: format!(
                            "{contig}\t{pos}\t.\t{r}\t{r}[{contig}:{mate}[\t60\tPASS\tSVTYPE=BND\tGT\t1/1"
                        ),
                        haplotype: [
                            &reference[pos - FLANK..pos],
                            &reference[mate - 1..mate - 1 + FLANK],
                        ]
                        .concat(),
                    }
                }
                _ => {
                    // Case 2, t]p]: REF[..=pos] + revcomp(MATE[..=mate_pos])
                    let mate = pos + 40_000;
                    Event {
                        class: Class::Bnd,
                        qname_prefix: format!(
                            "EIDOLON_chimeric_{contig}_{}_{contig}_{mate}_",
                            pos - 1
                        ),
                        vcf: format!(
                            "{contig}\t{pos}\t.\t{r}\t{r}]{contig}:{mate}]\t60\tPASS\tSVTYPE=BND\tGT\t1/1"
                        ),
                        haplotype: [
                            reference[pos - FLANK..pos].to_vec(),
                            revcomp(&reference[mate - FLANK..mate]),
                        ]
                        .concat(),
                    }
                }
            };
            out.push(ev);
        }
    }
    out
}

fn write_vcf(path: &Path, events: &[Event]) {
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
    for e in events {
        writeln!(f, "{}", e.vcf).unwrap();
    }
}

/// Best ungapped placement of `read` on `hay`, seeded by exact k-mers. Returns (offset,
/// identity). Seeds are tried along the whole read, so the clean half always contributes some.
fn place(read: &[u8], hay: &[u8], index: &HashMap<&[u8], Vec<usize>>) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    let mut tried = std::collections::HashSet::new();
    for off in (0..=read.len() - SEED).step_by(4) {
        let Some(hits) = index.get(&read[off..off + SEED]) else {
            continue;
        };
        for &p in hits {
            if p < off || p - off + read.len() > hay.len() || !tried.insert(p - off) {
                continue;
            }
            let start = p - off;
            let same = read
                .iter()
                .zip(&hay[start..start + read.len()])
                .filter(|(a, b)| a == b)
                .count();
            let id = same as f64 / read.len() as f64;
            if best.is_none_or(|(_, b)| id > b) {
                best = Some((start, id));
            }
        }
    }
    best
}

fn kmer_index(hay: &[u8]) -> HashMap<&[u8], Vec<usize>> {
    let mut m: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for p in 0..=hay.len() - SEED {
        m.entry(&hay[p..p + SEED]).or_default().push(p);
    }
    m
}

/// Mismatch tallies for one stratum: `cells[from][to]` against the top strand, plus mismatches
/// by CYCLE of the emitted read.
#[derive(Default)]
struct Tally {
    reads: usize,
    cells: [[usize; 4]; 4],
    by_cycle: Vec<usize>,
}

impl Tally {
    fn new() -> Self {
        Self {
            by_cycle: vec![0; READ_LEN],
            ..Default::default()
        }
    }

    /// `top` is the read on the molecule's top strand; `truth` the matching top-strand bases.
    /// `reverse` says the emitted read is the reverse complement of `top` (an R2).
    fn add(&mut self, top: &[u8], truth: &[u8], reverse: bool) {
        self.reads += 1;
        for i in 0..top.len() {
            let (Some(h), Some(b)) = (base_index(truth[i]), base_index(top[i])) else {
                continue;
            };
            if h != b {
                self.cells[h][b] += 1;
                let cycle = if reverse { top.len() - 1 - i } else { i };
                self.by_cycle[cycle] += 1;
            }
        }
    }

    fn total(&self) -> usize {
        self.cells.iter().flatten().sum()
    }

    fn in_frame(&self) -> usize {
        (0..4).map(|h| self.cells[h][planted(h)]).sum()
    }

    fn rc_frame(&self) -> usize {
        (0..4)
            .map(|h| self.cells[h][reverse_complement_frame(h)])
            .sum()
    }
}

struct Measured {
    junction_r1: HashMap<Class, Tally>,
    junction_r2: HashMap<Class, Tally>,
    /// Junction R2 quality strings, for the orientation check.
    junction_r2_quals: Vec<String>,
    ordinary_r2: Tally,
    unplaced: Vec<String>,
    /// Chimeric reads the test could not attribute to a planted event.
    unattributed: usize,
}

/// One gen-reads run shared by the first two tests: it is the expensive part.
fn measured() -> &'static Measured {
    static RUN: std::sync::OnceLock<Measured> = std::sync::OnceLock::new();
    // Substitution-only, so every junction read places ungapped.
    RUN.get_or_init(|| {
        let (_dir, work) = fresh_workdir();
        let events = simulate(&work, 0.0, 0.0);
        measure(&work, &events)
    })
}

/// Run gen-reads over the planted events; returns them for attribution.
fn simulate(work: &Path, indel_probability: f64, insertion_fraction: f64) -> Vec<Event> {
    let (contig, reference) = load_reference();
    let events = events(&contig, &reference);
    let vcf = work.join("events.vcf");
    write_vcf(&vcf, &events);

    let mut config = GenReadsConfig::new(ecoli(), work.to_path_buf(), "frame");
    config.read_len = READ_LEN;
    config.coverage = 4;
    config.paired_ended = true;
    config.fragment_mean = Some(350.0);
    config.fragment_st_dev = Some(30.0);
    config.produce_fastq = true;
    config.produce_bam = true;
    config.input_vcf = Some(vcf);
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);
    config.sequence_error_model = Some(write_model(work, indel_probability, insertion_fraction));
    config.rng_seed = "777".to_string();
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    events
}

fn measure(work: &Path, events: &[Event]) -> Measured {
    let (_, reference) = load_reference();

    let indexes: Vec<HashMap<&[u8], Vec<usize>>> =
        events.iter().map(|e| kmer_index(&e.haplotype)).collect();

    let mut m = Measured {
        junction_r1: HashMap::new(),
        junction_r2: HashMap::new(),
        junction_r2_quals: Vec::new(),
        ordinary_r2: Tally::new(),
        unplaced: Vec::new(),
        unattributed: 0,
    };
    for (mate, file) in [(1, "frame_r1.fastq.gz"), (2, "frame_r2.fastq.gz")] {
        let lines = read_gzip_fastq_lines(&work.join(file));
        for rec in lines.chunks(4) {
            let name = rec[0].trim_start_matches('@');
            if !name.starts_with("EIDOLON_chimeric") {
                continue;
            }
            let Some(ev) = events
                .iter()
                .position(|e| name.starts_with(&e.qname_prefix))
            else {
                m.unattributed += 1;
                continue;
            };
            let emitted = rec[1].as_bytes();
            assert_eq!(emitted.len(), READ_LEN, "{name}: read length");
            let top = if mate == 2 {
                revcomp(emitted)
            } else {
                emitted.to_vec()
            };
            let hay = &events[ev].haplotype;
            match place(&top, hay, &indexes[ev]) {
                Some((start, id)) if id >= 0.75 => {
                    let tally = if mate == 2 {
                        m.junction_r2_quals.push(rec[3].clone());
                        &mut m.junction_r2
                    } else {
                        &mut m.junction_r1
                    };
                    tally
                        .entry(events[ev].class)
                        .or_insert_with(Tally::new)
                        .add(&top, &hay[start..start + READ_LEN], mate == 2);
                }
                other => m.unplaced.push(format!("{name}: {other:?}")),
            }
        }
    }

    // Ordinary R2 off the golden BAM, which stores reverse records top-strand oriented.
    let mut reader = bam::io::reader::Builder
        .build_from_path(work.join("frame.bam"))
        .unwrap();
    let _ = reader.read_header().unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let flags = record.flags();
        if !flags.is_last_segment() || !flags.is_reverse_complemented() {
            continue;
        }
        let name = record.name().map(|n| n.to_string()).unwrap_or_default();
        if name.starts_with("EIDOLON_chimeric") {
            continue;
        }
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        let ops: Vec<_> = record.cigar().iter().map(|o| o.unwrap()).collect();
        if ops.len() != 1 || ops[0].kind() != Kind::Match || ops[0].len() != READ_LEN {
            continue;
        }
        let begin = usize::from(start) - 1;
        if begin + READ_LEN > reference.len() {
            continue;
        }
        let seq: Vec<u8> = record.sequence().iter().collect();
        m.ordinary_r2
            .add(&seq, &reference[begin..begin + READ_LEN], true);
    }
    m
}

fn describe(label: &str, t: &Tally) -> String {
    format!(
        "{label}: reads={} mismatches={} in-frame={} reverse-complement-frame={}",
        t.reads,
        t.total(),
        t.in_frame(),
        t.rc_frame()
    )
}

fn assert_in_frame(label: &str, t: &Tally, min_errors: usize) {
    assert!(
        t.total() >= min_errors,
        "{label}: only {} mismatches, too few to decide the frame",
        t.total()
    );
    let frac = t.in_frame() as f64 / t.total() as f64;
    assert!(
        frac >= 0.99,
        "{label}: {:.3} of mismatches follow the planted matrix in reference orientation and \
         {} follow its reverse complement. R2 must draw its errors on the top strand and flip \
         afterwards, like the ordinary paired writer (#777).\n{}",
        frac,
        t.rc_frame(),
        describe(label, t)
    );
}

#[test]
fn junction_r2_reads_draw_substitutions_in_reference_orientation() {
    let m = measured();

    // Denominators first: every chimeric read attributed and placed.
    eprintln!(
        "junction reads unattributed: {}, unplaced: {}",
        m.unattributed,
        m.unplaced.len()
    );
    assert_eq!(m.unattributed, 0, "chimeric reads matched no planted event");
    let placed: usize = m
        .junction_r1
        .values()
        .chain(m.junction_r2.values())
        .map(|t| t.reads)
        .sum();
    assert!(
        m.unplaced.len() * 100 <= placed,
        "{} of {} junction reads could not be placed on their derived haplotype:\n{}",
        m.unplaced.len(),
        placed + m.unplaced.len(),
        m.unplaced.join("\n")
    );

    // Every stratum printed before any assertion, so a failure shows the whole picture.
    for (label, map) in [("R1", &m.junction_r1), ("R2", &m.junction_r2)] {
        for (class, t) in map {
            eprintln!("{}", describe(&format!("{class:?} junction {label}"), t));
        }
    }
    eprintln!("{}", describe("ordinary R2", &m.ordinary_r2));

    for class in [Class::Del, Class::Dup, Class::Inv, Class::Bnd] {
        let r2 = m
            .junction_r2
            .get(&class)
            .unwrap_or_else(|| panic!("{class:?}: no junction R2 reads at all"));
        let r1 = m
            .junction_r1
            .get(&class)
            .unwrap_or_else(|| panic!("{class:?}: no junction R1 reads at all"));
        // PER_TYPE events x coverage 4 frags each (INV: two junctions), less nothing: the
        // model has no indels, so no pair can truncate.
        assert!(
            r2.reads >= PER_TYPE * 4 && r1.reads == r2.reads,
            "{class:?}: {} R1 / {} R2 junction reads placed, expected >= {} of each",
            r1.reads,
            r2.reads,
            PER_TYPE * 4
        );
        assert_in_frame(&format!("{class:?} junction R1"), r1, 200);
        assert_in_frame(&format!("{class:?} junction R2"), r2, 200);
    }
    assert_in_frame("ordinary R2", &m.ordinary_r2, 10_000);
}

/// #734 on the junction path: R2 cycle 1 is its first emitted base, so the Q40 half comes first
/// in the FASTQ and the errors sit in the last HALF cycles.
#[test]
fn junction_r2_quality_and_errors_run_in_cycle_order() {
    let m = measured();
    assert!(
        m.junction_r2_quals.len() >= 4 * PER_TYPE * 4,
        "only {} junction R2 reads",
        m.junction_r2_quals.len()
    );
    let hi = (Q_HIGH as u8 + 33) as char;
    let lo = (Q_LOW as u8 + 33) as char;
    let expected: String = std::iter::repeat_n(hi, HALF)
        .chain(std::iter::repeat_n(lo, READ_LEN - HALF))
        .collect();
    for q in &m.junction_r2_quals {
        assert_eq!(
            q, &expected,
            "junction R2 quality must be in cycle order: Q{Q_HIGH} for cycles 1-{HALF}"
        );
    }
    for (class, t) in &m.junction_r2 {
        let head: usize = t.by_cycle[..HALF].iter().sum();
        let tail: usize = t.by_cycle[HALF..].iter().sum();
        let bases = (t.reads * HALF) as f64;
        eprintln!(
            "{class:?} junction R2: error rate cycles 1-{HALF} {:.5}, cycles {}-{READ_LEN} {:.5}",
            head as f64 / bases,
            HALF + 1,
            tail as f64 / bases
        );
        // Q40 is 1e-4 per base and Q7 is 0.20. The planted matrix has a zero diagonal, so
        // every error drawn is a visible mismatch.
        assert!(
            (head as f64 / bases) < 0.01 && (tail as f64 / bases) > 0.1,
            "{class:?}: junction R2 errors must follow cycle quality: head {head}, tail {tail} \
             over {} reads",
            t.reads
        );
    }
}

/// MUST NOT DROP. R2's forward window ends exactly at the fragment's right end, so a
/// sequencing deletion there needs bases past the end or the read truncates and the pair is
/// silently skipped (`TruncatedRead` is caught as a debug message). `fragment_tail_pad` supplies
/// them. Here a quarter of errors are deletions, ~2.5 per R2 in its Q7 half, and every planted
/// fragment must still come out as a full-length pair: `num_frags` is the coverage for a
/// homozygous event, so the count is known in advance (4 per event, per INV junction).
///
/// The load is deliberately moderate. The buffer is `read_len` bases, the same one the
/// ordinary writer pads its fragments with, so R2 truncates once its deletions total more than
/// a read length. Measured with EVERY error a deletion (~10 per R2): one pair in 240 still
/// truncated that way, as ordinary R2s do under that load. That is parity with the ordinary
/// path, not a defect of the buffer, and this test does not pin it.
#[test]
fn junction_pairs_survive_r2_deletion_errors() {
    let (_dir, work) = fresh_workdir();
    let events = simulate(&work, 0.25, 0.0);
    let mut expected: HashMap<Class, usize> = HashMap::new();
    for e in &events {
        *expected.entry(e.class).or_default() += if e.class == Class::Inv { 8 } else { 4 };
    }
    for (mate, file) in [(1, "frame_r1.fastq.gz"), (2, "frame_r2.fastq.gz")] {
        let lines = read_gzip_fastq_lines(&work.join(file));
        let mut got: HashMap<Class, usize> = HashMap::new();
        for rec in lines.chunks(4) {
            let name = rec[0].trim_start_matches('@');
            let Some(e) = events.iter().find(|e| name.starts_with(&e.qname_prefix)) else {
                continue;
            };
            assert_eq!(rec[1].len(), READ_LEN, "{name}: truncated read emitted");
            *got.entry(e.class).or_default() += 1;
        }
        eprintln!("R{mate}: junction reads per class {got:?}, expected {expected:?}");
        assert_eq!(
            got, expected,
            "R{mate}: junction pairs were dropped under deletion errors. R2 needs a buffer past \
             the fragment end to consume (#777)."
        );
    }

    // Non-vacuity: the deletions were really drawn on junction R2. The golden BAM records each
    // as a `D` op in the CIGAR.
    let mut reader = bam::io::reader::Builder
        .build_from_path(work.join("frame.bam"))
        .unwrap();
    let _ = reader.read_header().unwrap();
    let (mut r2, mut with_del) = (0usize, 0usize);
    for result in reader.records() {
        let record = result.unwrap();
        let name = record.name().map(|n| n.to_string()).unwrap_or_default();
        if !record.flags().is_last_segment() || !name.starts_with("EIDOLON_chimeric") {
            continue;
        }
        r2 += 1;
        if record
            .cigar()
            .iter()
            .any(|op| op.unwrap().kind() == Kind::Deletion)
        {
            with_del += 1;
        }
    }
    eprintln!("junction R2 records carrying a deletion: {with_del} of {r2}");
    assert!(
        // 1 - exp(-2.5) = 0.92 of R2s draw at least one deletion.
        r2 > 0 && with_del * 4 >= r2 * 3,
        "only {with_del} of {r2} junction R2 records carry a deletion; the fixture is not \
         exercising the buffer"
    );
}
