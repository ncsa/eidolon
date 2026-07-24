//! End-to-end test for the native `eidolon gen-cancer-reads` subcommand (#239).
//!
//! Drives a tumor/normal simulation on H1N1 and checks the orchestration +
//! merges: per-pass golden VCFs + tagged/concatenated FASTQs + an origin-tagged
//! truth VCF (germline | somatic | shared).

mod common;

use common::{eidolon, fresh_workdir, h1n1_reference};
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};

fn read_gz_lines(path: &std::path::Path) -> Vec<String> {
    let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(path).unwrap()));
    r.lines().map(|l| l.unwrap()).collect()
}

#[test]
fn gen_cancer_reads_produces_tagged_fastqs_and_origin_truth() {
    let (_dir, work) = fresh_workdir();
    let yaml = work.join("cancer.yml");
    // High per-pass mutation rates so H1N1's tiny (~14 kb) reference yields plenty
    // of germline + somatic SNVs to classify. SVs left off (finicky on H1N1's
    // short segments); the merge logic is type-agnostic.
    std::fs::write(
        &yaml,
        format!(
            "reference: {ref}\n\
             output_dir: {out}\n\
             output_prefix: ctest\n\
             total_coverage: 30\n\
             purity: 0.5\n\
             read_len: 70\n\
             paired_ended: true\n\
             fragment_mean: 250\n\
             fragment_st_dev: 30\n\
             normal_mutation_rate: 0.01\n\
             tumor_mutation_rate: 0.01\n\
             overwrite_output: true\n\
             rng_seed: cancer-e2e\n",
            ref = h1n1_reference().display(),
            out = work.display(),
        ),
    )
    .unwrap();

    eidolon()
        .args(["gen-cancer-reads", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    // Per-pass + merged outputs exist.
    for f in [
        "ctest_normal.vcf.gz",
        "ctest_tumor.vcf.gz",
        "ctest_merged_r1.fastq.gz",
        "ctest_merged_r2.fastq.gz",
        "ctest_merged_truth.vcf.gz",
    ] {
        assert!(work.join(f).is_file(), "expected output {f}");
    }

    // Merged FASTQ read names carry both N_ and T_ tags.
    let r1 = read_gz_lines(&work.join("ctest_merged_r1.fastq.gz"));
    let headers: Vec<&String> = r1.iter().step_by(4).collect();
    assert!(
        headers.iter().any(|h| h.starts_with("@N_")),
        "no N_-tagged reads"
    );
    assert!(
        headers.iter().any(|h| h.starts_with("@T_")),
        "no T_-tagged reads"
    );

    // Truth VCF: NEAT_ORIGIN present, with somatic AND shared (germline carried
    // through the tumor pass) classes represented.
    let truth = read_gz_lines(&work.join("ctest_merged_truth.vcf.gz"));
    assert!(
        truth.iter().any(|l| l.contains("##INFO=<ID=NEAT_ORIGIN")),
        "truth VCF missing NEAT_ORIGIN header declaration"
    );
    let body: Vec<&String> = truth.iter().filter(|l| !l.starts_with('#')).collect();
    let has = |tag: &str| {
        body.iter()
            .any(|l| l.contains(&format!("NEAT_ORIGIN={tag}")))
    };
    assert!(has("somatic"), "no somatic records in truth");
    assert!(
        has("shared"),
        "no shared (germline-carried) records in truth"
    );
    // every body record must be origin-tagged
    assert!(
        body.iter().all(|l| l.contains("NEAT_ORIGIN=")),
        "some truth records lack NEAT_ORIGIN"
    );
}

/// End-to-end contract for #405 (generative subclonal VAF): with a `subclones:`
/// architecture, de-novo somatic variants must spread across the configured
/// cancer-cell fractions instead of collapsing to the Genotype default
/// (het 0.5 / hom 1.0). A minor subclone at CCF 0.3 is the discriminating signal:
/// its measured AF (~0.3, in the tumor-pass golden VCF the truth VCF carries
/// through) is a value the pre-#405 single-fraction model cannot produce.
#[test]
fn subclonal_somatic_variants_span_a_vaf_spectrum() {
    let (_dir, work) = fresh_workdir();
    let yaml = work.join("cancer_subclonal.yml");
    // High somatic rate → many somatic SNVs; high coverage → tight per-site AF.
    // Two equal-weight subclones: clonal (1.0) and a minor subclone (0.3).
    std::fs::write(
        &yaml,
        format!(
            "reference: {ref}\n\
             output_dir: {out}\n\
             output_prefix: sctest\n\
             total_coverage: 160\n\
             purity: 0.5\n\
             read_len: 70\n\
             paired_ended: true\n\
             fragment_mean: 250\n\
             fragment_st_dev: 30\n\
             normal_mutation_rate: 0.001\n\
             tumor_mutation_rate: 0.02\n\
             subclones:\n\
             \x20 - {{ccf: 1.0, weight: 1.0}}\n\
             \x20 - {{ccf: 0.3, weight: 1.0}}\n\
             overwrite_output: true\n\
             rng_seed: cancer-subclonal-e2e\n",
            ref = h1n1_reference().display(),
            out = work.display(),
        ),
    )
    .unwrap();

    eidolon()
        .args(["gen-cancer-reads", "-c"])
        .arg(&yaml)
        .assert()
        .success();

    // Measured AF (final subfield of the GT:AD:DP:AF sample column) for every
    // somatic-tagged truth record.
    let truth = read_gz_lines(&work.join("sctest_merged_truth.vcf.gz"));
    let somatic_afs: Vec<f64> = truth
        .iter()
        .filter(|l| !l.starts_with('#') && l.contains("NEAT_ORIGIN=somatic"))
        .filter_map(|l| {
            let sample = l.split('\t').next_back()?;
            sample.split(':').next_back()?.parse::<f64>().ok()
        })
        .collect();

    assert!(
        somatic_afs.len() >= 10,
        "need a meaningful somatic sample to judge the spectrum; got {}: {somatic_afs:?}",
        somatic_afs.len()
    );

    // The minor-subclone cluster: AFs near 0.3 that the het(0.5)/hom(1.0) default
    // cannot produce. Window [0.15, 0.45] is clear of both defaults.
    let minor = somatic_afs
        .iter()
        .filter(|&&a| (0.15..=0.45).contains(&a))
        .count();
    // The clonal cluster: high AFs.
    let clonal = somatic_afs.iter().filter(|&&a| a >= 0.7).count();

    assert!(
        minor > 0,
        "no somatic AF near the 0.3 subclone — subclonal CCF not applied: {somatic_afs:?}"
    );
    assert!(
        clonal > 0,
        "no somatic AF near the clonal (1.0) subclone: {somatic_afs:?}"
    );
    // A genuine spectrum, not a single fraction: both clusters populated.
    assert!(
        minor >= somatic_afs.len() / 5 && clonal >= somatic_afs.len() / 5,
        "expected both subclones well-represented (minor={minor}, clonal={clonal}, \
         n={}): {somatic_afs:?}",
        somatic_afs.len()
    );
}
