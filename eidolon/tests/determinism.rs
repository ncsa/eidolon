//! Determinism guarantees for the gen-reads pipeline.
//!
//! Reproducibility is a load-bearing promise of any read simulator — without it,
//! variant-calling pipelines can't be debugged or regression-tested. These tests run the
//! real binary twice with the same seed and assert the outputs match.
//!
//! The invariant we test is **same seed → byte-identical output**, and separately the
//! weaker **same record multiset**.
//!
//! This file used to claim that line order could vary even at `num_threads=1`, because
//! eidolon iterated HashMap keys when assembling per-contig outputs and Rust randomizes
//! that per process. That was true, and #599 fixed it: `ContigContext` now carries an
//! explicit `contig_order`. The rationale outlived the defect, and left the README's
//! headline claim — "the same seed always produces the same reads **in the same order**"
//! — pinned by nothing, since every test here sorted before comparing.
//!
//! Measured on H1N1 at `sv_rate_scale=50` (38 SVs actually placed, so not vacuous):
//! decompressed FASTQ is byte-identical across repeated runs and across
//! `num_threads` 1/4/8, and so is the no-SV path. `byte_identical_output_*` below pins
//! that. The multiset tests are kept as the weaker floor.
//!
//! We exercise the invariant under both num_threads=1 and the rayon default so a
//! single-thread-only determinism leak doesn't slip past us, and also assert that the
//! seed argument actually changes output (otherwise "determinism" is a coincidence).
//!
//! We compare decompressed FASTQ contents rather than raw gzip bytes: gzip headers can
//! embed mtime/OS metadata that varies across runs even when the payload is identical.

mod common;

use common::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines};

fn run_gen_reads(seed: &str, output_dir: &std::path::Path, name: &str, threads: Option<usize>) {
    run_gen_reads_cfg(seed, output_dir, name, threads, None)
}

/// As above, but able to switch on de novo SVs — which is the path #599 lived on.
fn run_gen_reads_cfg(
    seed: &str,
    output_dir: &std::path::Path,
    name: &str,
    threads: Option<usize>,
    sv_rate_scale: Option<f64>,
) {
    let mut config = GenReadsConfig::new(h1n1_reference(), output_dir.to_path_buf(), name);
    config.rng_seed = seed.to_string();
    config.num_threads = threads;
    config.sv_rate_scale = sv_rate_scale;
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
}

/// Group a flat FASTQ line list into 4-line records, then sort lexicographically. This
/// turns line order into a multiset comparison.
fn sorted_records(lines: Vec<String>) -> Vec<[String; 4]> {
    assert!(
        lines.len().is_multiple_of(4),
        "FASTQ line count must be multiple of 4"
    );
    let mut records: Vec<[String; 4]> = lines
        .chunks(4)
        .map(|c| [c[0].clone(), c[1].clone(), c[2].clone(), c[3].clone()])
        .collect();
    records.sort();
    records
}

#[test]
fn same_seed_single_threaded_produces_same_record_multiset() {
    // num_threads=1 removes inter-contig parallelism. If determinism leaks even here,
    // some shared HashMap or atomic counter is consuming RNG non-deterministically and
    // the issue is *not* just a benign per-thread reordering — it's a real bug.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    let seed = "eidolon single threaded";
    run_gen_reads(seed, &dir_a, "st_a", Some(1));
    run_gen_reads(seed, &dir_b, "st_b", Some(1));

    let a = sorted_records(read_gzip_fastq_lines(&dir_a.join("st_a_r1.fastq.gz")));
    let b = sorted_records(read_gzip_fastq_lines(&dir_b.join("st_b_r1.fastq.gz")));

    assert!(!a.is_empty(), "first run produced no records");
    assert_eq!(
        a.len(),
        b.len(),
        "same seed (single-threaded) produced different record counts: a={}, b={}",
        a.len(),
        b.len(),
    );
    assert_eq!(
        a, b,
        "same seed (single-threaded) produced different record multisets. \
         A record was dropped, duplicated, or its content differs.",
    );
}

#[test]
fn same_seed_multi_threaded_produces_same_record_multiset() {
    // With parallel processing the per-contig output order can vary, but the set of
    // records (including names and quality scores) must be identical. If this fails,
    // the RNG is being shared across threads in a way that races, or reads are being
    // dropped/duplicated under contention.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    let seed = "eidolon multi threaded";
    run_gen_reads(seed, &dir_a, "mt_a", None);
    run_gen_reads(seed, &dir_b, "mt_b", None);

    let a = sorted_records(read_gzip_fastq_lines(&dir_a.join("mt_a_r1.fastq.gz")));
    let b = sorted_records(read_gzip_fastq_lines(&dir_b.join("mt_b_r1.fastq.gz")));

    assert!(!a.is_empty(), "first run produced no records");
    assert_eq!(
        a.len(),
        b.len(),
        "same seed produced different record counts (multi-threaded): a={}, b={}",
        a.len(),
        b.len(),
    );
    assert_eq!(
        a, b,
        "same seed produced different record multisets (multi-threaded). \
         Some record was dropped, duplicated, or differs between runs.",
    );
}

#[test]
fn different_seeds_produce_different_output() {
    // Sanity check that the seed argument is actually load-bearing. If this fails, the
    // RNG isn't actually being seeded from the config string and "determinism" is a
    // coincidence rather than a guarantee.
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();

    run_gen_reads("alpha seed text", &dir_a, "diff_a", Some(1));
    run_gen_reads("bravo seed text", &dir_b, "diff_b", Some(1));

    let a = read_gzip_fastq_lines(&dir_a.join("diff_a_r1.fastq.gz"));
    let b = read_gzip_fastq_lines(&dir_b.join("diff_b_r1.fastq.gz"));

    assert!(!a.is_empty() && !b.is_empty());
    assert_ne!(
        a, b,
        "different seeds produced identical FASTQ — seed argument is not load-bearing",
    );
}

/// **#599.** With SVs switched on, the same seed produced a DIFFERENT SET of reads on
/// every run — not merely a different order.
///
/// Measured before the fix, four consecutive single-threaded runs at one seed:
/// 2259, 2265, 2271, 2268 R1 records. The truth VCF was byte-identical all four times
/// (213 SVs), so SV sampling was fine; the divergence was in chimeric junction reads.
/// `process_chimeric_variants` walked `ctx.mutated_maps` — a HashMap, iteration order
/// randomized per process — while threading ONE `rng` through the whole loop, so contig
/// order decided which draws each junction received. The fragment-length retry loop could
/// then settle after a different number of draws, changing the read COUNT too.
///
/// Every other test in this file leaves SVs off, which is exactly why none of them saw it.
/// The multiset comparison used here is the same invariant the rest of the file asserts —
/// line order may still vary, the record set may not.
/// The README's actual claim: same seed, same reads, SAME ORDER — and unchanged by thread
/// count. Every other test here sorts first, so an ordering regression (a reintroduced
/// HashMap iteration, a rayon collect that no longer preserves order) would pass all of
/// them. #599 was exactly that class of defect.
///
/// Compares decompressed bytes: gzip headers can carry mtime/OS metadata that varies even
/// when the payload does not.
#[test]
fn byte_identical_output_with_svs_across_thread_counts() {
    let seed = "sv byte identity seed";
    let mut baseline: Option<Vec<String>> = None;

    // 1 and 8 threads, each run twice: catches a per-process ordering leak (same thread
    // count, different run) and a thread-count-dependent one.
    for threads in [1usize, 8] {
        for rep in 0..2 {
            let (_g, dir) = fresh_workdir();
            let name = format!("bytes_t{threads}_{rep}");
            run_gen_reads_cfg(seed, &dir, &name, Some(threads), Some(50.0));
            let lines = read_gzip_fastq_lines(&dir.join(format!("{name}_r1.fastq.gz")));
            assert!(
                !lines.is_empty(),
                "no reads generated at threads={threads} — nothing to compare"
            );
            match &baseline {
                None => baseline = Some(lines),
                Some(first) => assert_eq!(
                    *first, lines,
                    "same seed produced different output at threads={threads} (rep {rep}). \
                     The README claims byte-identical output in the same order regardless \
                     of thread count; this is #599's class of defect."
                ),
            }
        }
    }
}

/// The README says output is in "the same order", and `runner.rs` documents that order as
/// **reference-file order, NOT HashMap order** (#599). Nothing pinned it: reversing a contig
/// iteration is still perfectly deterministic, so every run-to-run comparison in this file —
/// the byte-identity one above included — passes a mutant that emits contigs backwards.
///
/// Scope, measured rather than assumed: the FASTQ has two sections. Ordinary reads come
/// first, grouped by contig in reference-file order. Chimeric reads are appended after, in
/// per-event order — a breakend joins two contigs, so a chimeric read's name carries one of
/// them and that section is legitimately not sorted by contig. Only the first section is
/// asserted here; chimeric block ordering is left unpinned rather than pinned to something
/// this test would have to guess at.
#[test]
fn ordinary_reads_are_emitted_in_reference_file_order() {
    let (_g, dir) = fresh_workdir();
    run_gen_reads_cfg("contig order seed", &dir, "order", Some(1), Some(50.0));
    let lines = read_gzip_fastq_lines(&dir.join("order_r1.fastq.gz"));

    let reference_order = [
        "H1N1_HA", "H1N1_MP", "H1N1_NA", "H1N1_NP", "H1N1_NS", "H1N1_PA", "H1N1_PB1",
        "H1N1_PB2",
    ];

    let mut emitted: Vec<String> = Vec::new();
    for line in lines
        .iter()
        .filter(|l| l.starts_with("@EIDOLON_generated_"))
    {
        let contig = reference_order
            .iter()
            .find(|c| line.contains(*c))
            .unwrap_or_else(|| panic!("read name names no known contig: {line}"));
        if emitted.last().map(String::as_str) != Some(*contig) {
            assert!(
                !emitted.iter().any(|e| e == contig),
                "ordinary reads for {contig} appear in more than one block — they are \
                 interleaved, and reference-file order is not being preserved"
            );
            emitted.push((*contig).to_string());
        }
    }

    assert!(
        emitted.len() > 1,
        "only {} contig(s) produced ordinary reads — the ordering claim would be vacuous",
        emitted.len()
    );
    let expected: Vec<&str> = reference_order
        .iter()
        .copied()
        .filter(|c| emitted.iter().any(|e| e == c))
        .collect();
    assert_eq!(
        emitted, expected,
        "contigs were emitted in a different order than they appear in the reference file. \
         runner.rs documents reference-file order (#599); the README promises it."
    );
}

#[test]
fn same_seed_with_svs_produces_the_same_record_multiset() {
    let (_a, dir_a) = fresh_workdir();
    let (_b, dir_b) = fresh_workdir();
    let seed = "sv determinism seed";

    // High enough that the H1N1 fixture actually carries junctions to generate reads for;
    // with none planted the test would pass vacuously.
    run_gen_reads_cfg(seed, &dir_a, "sv_a", Some(1), Some(50.0));
    run_gen_reads_cfg(seed, &dir_b, "sv_b", Some(1), Some(50.0));

    // The shared default config is single-end, so there is one mate file, matching the
    // other tests in this file.
    let a = sorted_records(read_gzip_fastq_lines(&dir_a.join("sv_a_r1.fastq.gz")));
    let b = sorted_records(read_gzip_fastq_lines(&dir_b.join("sv_b_r1.fastq.gz")));

    assert!(
        !a.is_empty(),
        "no reads generated — the SV run produced nothing to compare"
    );
    assert_eq!(
        a.len(),
        b.len(),
        "same seed produced {} records in one run and {} in the other — the read SET is \
         unstable, not just its order (#599)",
        a.len(),
        b.len()
    );
    assert_eq!(
        a, b,
        "same seed produced a different record multiset across two runs (#599)"
    );
}
