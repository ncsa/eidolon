//! Model-parity regression tests.
//!
//! Each model builder (`gen-mut-model`, `gen-seq-error-model`,
//! `gen-frag-length-model`, `gen-gc-bias-model`) is run against a fixed input
//! set and the resulting model's canonical-JSON form is compared against a
//! checked-in baseline under `test_data/baseline_models/`.
//!
//! Why canonical JSON, not raw bytes: model structs contain `HashMap` fields
//! (e.g. `MutationModel::statistical_models`). Rust's default HashMap uses a
//! randomized hasher seed per process, so two runs of the model builder
//! serialize the same data with different key/array order. The canonical
//! form parses the gzipped JSON, sorts object keys, and sorts arrays-of-
//! pairs (the JSON shape that serde uses for tuple-keyed HashMaps) by their
//! serialized key. The output is stable across processes.
//!
//! Regenerating baselines: see `scripts/regenerate_model_baselines.sh`. Only
//! do this when an intentional algorithmic change is being made; treat any
//! unexpected baseline drift as a regression.

mod common;

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use common::eidolon;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde_json::Value;

// ── Canonicalization ─────────────────────────────────────────────────────────

/// Recursively canonicalize a JSON `Value`:
///   - Object keys: already lexically sorted by `serde_json::Map`'s default
///     BTreeMap backing — recurse into values.
///   - Arrays: recurse into items, then if every item is a 2-element array
///     whose first element is itself a composite (Array or Object), treat
///     the parent as a HashMap-with-tuple-key serialization and sort by the
///     serialized first element. This catches `HashMap<(K1,K2), V>` and
///     `HashMap<EnumKey, V>` without disturbing legitimate ordered lists
///     of scalar pairs.
fn canonicalize(v: Value) -> Value {
    match v {
        Value::Object(map) => {
            let entries: Vec<(String, Value)> =
                map.into_iter().map(|(k, v)| (k, canonicalize(v))).collect();
            Value::Object(entries.into_iter().collect())
        }
        Value::Array(items) => {
            let mut canonical: Vec<Value> = items.into_iter().map(canonicalize).collect();
            let looks_like_map_entries = !canonical.is_empty()
                && canonical.iter().all(|e| {
                    matches!(
                        e,
                        Value::Array(a) if a.len() == 2 && matches!(a[0], Value::Array(_) | Value::Object(_))
                    )
                });
            if looks_like_map_entries {
                canonical.sort_by_key(|item| {
                    let arr = item.as_array().unwrap();
                    serde_json::to_string(&arr[0]).unwrap()
                });
            }
            Value::Array(canonical)
        }
        other => other,
    }
}

/// Read a gzipped JSON model, parse, canonicalize, and pretty-print.
fn canonical_model_json(path: &Path) -> String {
    let mut decoder = GzDecoder::new(fs::File::open(path).unwrap_or_else(|e| {
        panic!("could not open {}: {e}", path.display());
    }));
    let mut raw = String::new();
    decoder
        .read_to_string(&mut raw)
        .expect("model file must decompress to valid UTF-8 JSON");
    let value: Value = serde_json::from_str(&raw).expect("model file must be valid JSON");
    let canonical = canonicalize(value);
    serde_json::to_string_pretty(&canonical).unwrap()
}

// ── Baselines ────────────────────────────────────────────────────────────────

/// `BLESS_BASELINES=1 cargo test --test model_parity` writes the freshly-built
/// canonical JSON to the baseline path instead of asserting equality. Used to
/// (re)generate baselines after an intentional algorithmic change. Treat any
/// drift not preceded by a BLESS commit as a regression.
fn bless_enabled() -> bool {
    std::env::var("BLESS_BASELINES").is_ok_and(|v| v == "1")
}

fn baseline_path(name: &str) -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/baseline_models/{name}.canonical.json.gz",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

fn read_gzipped_baseline(path: &Path) -> std::io::Result<String> {
    let f = fs::File::open(path)?;
    let mut decoder = GzDecoder::new(f);
    let mut s = String::new();
    decoder.read_to_string(&mut s)?;
    Ok(s)
}

fn write_gzipped_baseline(path: &Path, canonical: &str) -> std::io::Result<()> {
    fs::create_dir_all(path.parent().unwrap())?;
    let f = fs::File::create(path)?;
    let mut encoder = GzEncoder::new(f, Compression::default());
    encoder.write_all(canonical.as_bytes())?;
    encoder.finish()?;
    Ok(())
}

fn assert_or_bless(canonical: &str, baseline_name: &str) {
    let path = baseline_path(baseline_name);
    if bless_enabled() {
        write_gzipped_baseline(&path, canonical).unwrap();
        eprintln!("[BLESS] wrote {}", path.display());
        return;
    }
    let expected = read_gzipped_baseline(&path).unwrap_or_else(|e| {
        panic!(
            "baseline {} missing or unreadable: {e}\n\
             Run with BLESS_BASELINES=1 to generate.",
            path.display()
        );
    });
    if expected != canonical {
        let diff_path = path.with_extension("").with_extension("actual.json");
        fs::write(&diff_path, canonical).ok();
        panic!(
            "Model parity drift detected for `{baseline_name}`.\n\
             Baseline: {}\n\
             Actual (written for inspection): {}\n\
             If this is intentional, rebuild baselines with BLESS_BASELINES=1.",
            path.display(),
            diff_path.display(),
        );
    }
}

// ── Input fixtures ───────────────────────────────────────────────────────────

fn h1n1_reference() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/references/H1N1.fa",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

fn small_snps_vcf() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/vcfs/small_snps.vcf",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

fn h1n1_read1_fastq_gz() -> PathBuf {
    PathBuf::from(format!(
        "{}/test_data/H1N1_read1.fq.gz",
        env!("CARGO_MANIFEST_DIR"),
    ))
}

/// Write a deterministic paired-end BAM for the H1N1_HA contig: `n_pairs`
/// proper pairs all at the same start/mate positions, plus `n_extra` unpaired
/// reads at position 1 to inflate coverage. Records are written in a fixed
/// order with fixed fields, so the resulting model fit is reproducible.
///
/// This intentionally mirrors the helper used in
/// `gen_bam_models::utils::runner` tests, so the test inputs match the
/// in-tree parity test fixtures.
fn write_test_bam(
    path: &Path,
    contig: &[u8],
    contig_len: usize,
    n_pairs: usize,
    tlen: usize,
    read_len: usize,
    n_extra_coverage_reads: usize,
) {
    use noodles::bam;
    use noodles::core::Position;
    use noodles::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags, MappingQuality,
                cigar::{Op, op::Kind},
            },
            record_buf::{Cigar, Sequence},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    let header = sam::Header::builder()
        .add_reference_sequence(
            contig.to_vec(),
            Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(contig_len).unwrap()),
        )
        .build();
    let file = fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
    let seq = vec![b'A'; read_len];

    let paired_first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    let paired_second = Flags::SEGMENTED;
    for _ in 0..n_pairs {
        let mate_start = tlen - read_len + 1;
        let mut r1 = RecordBuf::default();
        *r1.flags_mut() = paired_first;
        *r1.cigar_mut() = cigar.clone();
        *r1.reference_sequence_id_mut() = Some(0);
        *r1.alignment_start_mut() = Position::new(1);
        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = Position::new(mate_start);
        *r1.template_length_mut() = tlen as i32;
        *r1.sequence_mut() = Sequence::from(seq.as_slice());
        *r1.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r1).unwrap();

        let mut r2 = RecordBuf::default();
        *r2.flags_mut() = paired_second;
        *r2.cigar_mut() = cigar.clone();
        *r2.reference_sequence_id_mut() = Some(0);
        *r2.alignment_start_mut() = Position::new(mate_start);
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = Position::new(1);
        *r2.template_length_mut() = -(tlen as i32);
        *r2.sequence_mut() = Sequence::from(seq.as_slice());
        *r2.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r2).unwrap();
    }

    for _ in 0..n_extra_coverage_reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::empty();
        *r.cigar_mut() = cigar.clone();
        *r.reference_sequence_id_mut() = Some(0);
        *r.alignment_start_mut() = Position::new(1);
        *r.sequence_mut() = Sequence::from(seq.as_slice());
        *r.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
        writer.write_alignment_record(&header, &r).unwrap();
    }
}

fn write_temp_yaml(content: &str) -> tempfile::NamedTempFile {
    use std::io::Write as _;
    let mut f = tempfile::Builder::new().suffix(".yml").tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

// ── Tests ────────────────────────────────────────────────────────────────────

/// `gen-frag-length-model --distribution normal` fits a Normal to the TLEN
/// distribution observed in a paired BAM. The synthetic BAM has 50 pairs all
/// at TLEN=200, so the fit collapses to mean=200, st_dev=0. The baseline file
/// pins this exact fit; any change to the histogram→Normal pipeline surfaces here.
///
/// EXPLICITLY `normal`, not the default. Since v3.3.0 the default is `discrete`,
/// and one distinct TLEN is refused there -- a distribution over a single value
/// has no shape. This test therefore pins the ESCAPE HATCH; the default is pinned
/// by `frag_length_discrete_model_matches_baseline` below, so neither family can
/// drift unnoticed.
#[test]
fn frag_length_model_matches_baseline() {
    let tmp = tempfile::tempdir().unwrap();
    let bam = tmp.path().join("input.bam");
    write_test_bam(&bam, b"H1N1_HA", 1700, 50, 200, 100, 0);
    let out = tmp.path().join("frag.json.gz");
    let yaml = write_temp_yaml(&format!(
        "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 2\n\
         distribution: normal\n",
        bam.display(),
        out.display()
    ));
    eidolon()
        .args(["gen-frag-length-model", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    let canonical = canonical_model_json(&out);
    assert_or_bless(&canonical, "frag_length");
}

/// The DEFAULT path since v3.3.0: `gen-frag-length-model` keeps the observed shape
/// instead of collapsing it to two floats. Pins the whole weight vector, so a silent
/// return to a Normal fit -- which is what discarded the right tail of every real
/// library -- cannot pass unnoticed.
#[test]
fn frag_length_discrete_model_matches_baseline() {
    let tmp = tempfile::tempdir().unwrap();
    let bam = tmp.path().join("input_spread.bam");
    // A spread of TLENs, deliberately with holes in it: lengths 150..250 in steps of 4.
    // The gaps are the point -- a discrete model built from this must come out
    // contiguous, and the baseline pins the bridged result.
    let tlens: Vec<usize> = (150..=250)
        .step_by(4)
        .flat_map(|l| std::iter::repeat_n(l, 6))
        .collect();
    write_test_bam_with_tlens(&bam, b"H1N1_HA", 1700, &tlens, 100);
    let out = tmp.path().join("frag_discrete.json.gz");
    let yaml = write_temp_yaml(&format!(
        "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 2\n",
        bam.display(),
        out.display()
    ));
    eidolon()
        .args(["gen-frag-length-model", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    // Assert the SHAPE of the artifact before pinning its bytes. A baseline that was
    // blessed from a Normal would otherwise keep passing forever.
    let mut raw = String::new();
    flate2::read::GzDecoder::new(fs::File::open(&out).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    let v: serde_json::Value = serde_json::from_str(&raw).unwrap();
    let dist = v
        .get("Discrete")
        .expect("the default must build a Discrete model")
        .get("distribution")
        .unwrap();
    let values: Vec<i64> = dist["values"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_i64().unwrap())
        .collect();
    assert_eq!(
        values.len() as i64,
        values[values.len() - 1] - values[0] + 1,
        "the built model must have no gaps, got {} bins spanning {}..{}",
        values.len(),
        values[0],
        values[values.len() - 1]
    );

    let canonical = canonical_model_json(&out);
    assert_or_bless(&canonical, "frag_length_discrete");
}

/// Like `write_test_bam` but lays one pair per supplied TLEN, so a test can control the
/// SHAPE of the observed distribution rather than only its centre.
fn write_test_bam_with_tlens(
    path: &Path,
    contig: &[u8],
    contig_len: usize,
    tlens: &[usize],
    read_len: usize,
) {
    use noodles::bam;
    use noodles::core::Position;
    use noodles::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags, MappingQuality,
                cigar::{Op, op::Kind},
            },
            record_buf::{Cigar, Sequence},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    let header = sam::Header::builder()
        .add_reference_sequence(
            contig.to_vec(),
            Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(contig_len).unwrap()),
        )
        .build();
    let file = fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    let cigar: Cigar = [Op::new(Kind::Match, read_len)].into_iter().collect();
    let seq = vec![b'A'; read_len];
    let paired_first = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
    let paired_second = Flags::SEGMENTED;

    for &tlen in tlens {
        let mate_start = tlen - read_len + 1;
        let mut r1 = RecordBuf::default();
        *r1.flags_mut() = paired_first;
        *r1.cigar_mut() = cigar.clone();
        *r1.reference_sequence_id_mut() = Some(0);
        *r1.alignment_start_mut() = Position::new(1);
        *r1.mate_reference_sequence_id_mut() = Some(0);
        *r1.mate_alignment_start_mut() = Position::new(mate_start);
        *r1.template_length_mut() = tlen as i32;
        *r1.mapping_quality_mut() = MappingQuality::new(60);
        *r1.sequence_mut() = Sequence::from(seq.clone());
        writer.write_alignment_record(&header, &r1).unwrap();

        let mut r2 = RecordBuf::default();
        *r2.flags_mut() = paired_second;
        *r2.cigar_mut() = cigar.clone();
        *r2.reference_sequence_id_mut() = Some(0);
        *r2.alignment_start_mut() = Position::new(mate_start);
        *r2.mate_reference_sequence_id_mut() = Some(0);
        *r2.mate_alignment_start_mut() = Position::new(1);
        *r2.template_length_mut() = -(tlen as i32);
        *r2.mapping_quality_mut() = MappingQuality::new(60);
        *r2.sequence_mut() = Sequence::from(seq.clone());
        writer.write_alignment_record(&header, &r2).unwrap();
    }
}

/// `gen-gc-bias-model` walks the reference in windows, accumulates per-GC%
/// mean coverage, and writes a 101-bin weight vector. With uniform coverage
/// (one unpaired read per base across the full contig), all bins should
/// collapse to neutral weight 1.0. Pins the entire weight vector via the
/// canonical baseline.
#[test]
fn gc_bias_model_matches_baseline() {
    let tmp = tempfile::tempdir().unwrap();
    let bam = tmp.path().join("input.bam");
    // H1N1_HA is 1701 bp. Write 5 unpaired reads stacked at position 1 — this
    // produces flat coverage across the first 100 bp of the contig only, but
    // that's enough for `min_windows_per_bin=1` to populate the model. Real
    // BAMs would obviously cover more, but model fitting is deterministic on
    // input regardless.
    write_test_bam(&bam, b"H1N1_HA", 1701, 0, 0, 100, 5);
    let out = tmp.path().join("gc.json.gz");
    let yaml = write_temp_yaml(&format!(
        "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n\
         min_mapq: 0\nbed_file: .\n\
         window_size: 100\nwindow_stride: 100\nmin_windows_per_bin: 1\n",
        h1n1_reference().display(),
        bam.display(),
        out.display()
    ));
    eidolon()
        .args(["gen-gc-bias-model", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    let canonical = canonical_model_json(&out);
    assert_or_bless(&canonical, "gc_bias");
}

/// `gen-mut-model` fits transition matrices and mutation-rate scalars from a
/// VCF against a reference. Inputs are tiny but the algorithmic surface is
/// large (trinuc counting, HashMap-keyed transitions). Any unintended
/// rewrite of the counting pipeline will surface here.
#[test]
fn mut_model_matches_baseline() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("mut.json.gz");
    let yaml = write_temp_yaml(&format!(
        "reference: {}\nvcf_file: {}\noutput_file: {}\noverwrite_output: true\nbed_file: .\n",
        h1n1_reference().display(),
        small_snps_vcf().display(),
        out.display()
    ));
    eidolon()
        .args(["gen-mut-model", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    let canonical = canonical_model_json(&out);
    assert_or_bless(&canonical, "mut_model");
}

/// `gen-seq-error-model` aggregates per-position quality-score distributions
/// from a FASTQ. The bundled H1N1 read1 FASTQ is the input; any change to
/// quality binning or the default transition matrix surfaces here.
#[test]
fn seq_error_model_matches_baseline() {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("seq_err.json.gz");
    let yaml = write_temp_yaml(&format!(
        "fastq_file: {}\noutput_file: {}\noverwrite_output: true\n",
        h1n1_read1_fastq_gz().display(),
        out.display()
    ));
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    let canonical = canonical_model_json(&out);
    assert_or_bless(&canonical, "seq_error");
}

/// The unified `gen-bam-models` GC-bias path must produce the same model as
/// the standalone `gen-gc-bias-model` against the same BAM. Complements the
/// in-tree fragment-length parity test in `gen_bam_models/utils/runner.rs`,
/// which only covers `frag_length`.
#[test]
fn gen_bam_models_gc_bias_byte_equal_to_standalone() {
    let tmp = tempfile::tempdir().unwrap();
    let bam = tmp.path().join("input.bam");
    write_test_bam(&bam, b"H1N1_HA", 1701, 0, 0, 100, 5);

    let unified_out = tmp.path().join("gc_unified.json.gz");
    let unified_yaml = write_temp_yaml(&format!(
        "bam_file: {}\nmin_mapq: 0\ngc_bias:\n  reference: {}\n  output_file: {}\n  \
         overwrite_output: true\n  window_size: 100\n  window_stride: 100\n  \
         min_windows_per_bin: 1\n",
        bam.display(),
        h1n1_reference().display(),
        unified_out.display()
    ));
    eidolon()
        .args(["gen-bam-models", "-c"])
        .arg(unified_yaml.path())
        .assert()
        .success();

    let standalone_out = tmp.path().join("gc_standalone.json.gz");
    let standalone_yaml = write_temp_yaml(&format!(
        "reference: {}\nbam_file: {}\noutput_file: {}\noverwrite_output: true\n\
         min_mapq: 0\nbed_file: .\n\
         window_size: 100\nwindow_stride: 100\nmin_windows_per_bin: 1\n",
        h1n1_reference().display(),
        bam.display(),
        standalone_out.display()
    ));
    eidolon()
        .args(["gen-gc-bias-model", "-c"])
        .arg(standalone_yaml.path())
        .assert()
        .success();

    let unified_canonical = canonical_model_json(&unified_out);
    let standalone_canonical = canonical_model_json(&standalone_out);
    assert_eq!(
        unified_canonical, standalone_canonical,
        "unified gen-bam-models GC-bias output diverged from standalone gen-gc-bias-model"
    );
}
