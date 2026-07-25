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

/// End-to-end contract for #405 (generative subclonal VAF): a `subclones:`
/// architecture spreads de-novo somatic variants across cancer-cell fractions,
/// and each CCF *composes with the variant's dosage* (alt = dosage × CCF) rather
/// than replacing it. With two subclones (CCF 1.0, 0.3) over het/hom de-novo
/// variants the measured somatic AFs land at {0.5·1, 1·1, 0.5·0.3, 1·0.3} =
/// {0.5, 1.0, 0.15, 0.3}. The discriminating signal is any AF below ~0.4 — a
/// value the pre-#405 dosage-only model (het 0.5 / hom 1.0) cannot produce.
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

    // The minor-subclone cluster: dosage × 0.3 ∈ {0.15, 0.3}, i.e. AFs at/below
    // ~0.4 that the dosage-only default (het 0.5 / hom 1.0) cannot produce.
    let minor = somatic_afs
        .iter()
        .filter(|&&a| (0.05..=0.4).contains(&a))
        .count();
    // The clonal cluster: dosage × 1.0 ∈ {0.5, 1.0}.
    let clonal = somatic_afs.iter().filter(|&&a| a >= 0.45).count();

    assert!(
        minor > 0,
        "no somatic AF below ~0.4 — subclonal CCF not composed in: {somatic_afs:?}"
    );
    assert!(
        clonal > 0,
        "no somatic AF at the clonal (CCF 1.0) fractions: {somatic_afs:?}"
    );
    // A genuine spectrum, not a single fraction: both clusters populated. With
    // equal subclone weights ~half the somatic burden falls in each.
    assert!(
        minor >= somatic_afs.len() / 5 && clonal >= somatic_afs.len() / 5,
        "expected both subclones well-represented (minor={minor}, clonal={clonal}, \
         n={}): {somatic_afs:?}",
        somatic_afs.len()
    );

    // ── INFO/NEAT_CCF ground-truth tag (#405) ────────────────────────────────
    // The header must declare it, every somatic record must carry the *intended*
    // CCF (one of the two configured values), and no germline/shared record may.
    assert!(
        truth.iter().any(|l| l.contains("##INFO=<ID=NEAT_CCF")),
        "merged truth missing NEAT_CCF header declaration"
    );

    let ccf_of = |line: &str| -> Option<f64> {
        let info = line.split('\t').nth(7)?;
        info.split(';')
            .find_map(|kv| kv.strip_prefix("NEAT_CCF="))?
            .parse::<f64>()
            .ok()
    };
    // Dosage from the GT string (alt copies / total called), mirroring
    // Variant::dosage_fraction — the sampler's per-copy alt probability.
    let dosage_of = |line: &str| -> Option<f64> {
        let gt = line.split('\t').next_back()?.split(':').next()?;
        let (mut alt, mut total) = (0u32, 0u32);
        for a in gt.split(['/', '|']) {
            match a {
                "." | "" => {}
                "0" => total += 1,
                _ => {
                    alt += 1;
                    total += 1;
                }
            }
        }
        (total > 0).then(|| alt as f64 / total as f64)
    };
    let dp_of = |line: &str| -> Option<u32> {
        line.split('\t')
            .next_back()?
            .split(':')
            .nth(2)?
            .parse()
            .ok()
    };

    let body: Vec<&String> = truth.iter().filter(|l| !l.starts_with('#')).collect();

    // Germline / shared records never carry a somatic CCF.
    for l in body.iter().filter(|l| !l.contains("NEAT_ORIGIN=somatic")) {
        assert!(
            ccf_of(l).is_none(),
            "non-somatic record carries NEAT_CCF: {l}"
        );
    }

    // Every somatic record carries a NEAT_CCF from the configured architecture.
    let somatic: Vec<&&String> = body
        .iter()
        .filter(|l| l.contains("NEAT_ORIGIN=somatic"))
        .collect();
    let mut checked = 0;
    let mut err_sum = 0.0;
    for l in &somatic {
        let ccf = ccf_of(l).unwrap_or_else(|| panic!("somatic record lacks NEAT_CCF: {l}"));
        assert!(
            (ccf - 1.0).abs() < 1e-6 || (ccf - 0.3).abs() < 1e-6,
            "unexpected NEAT_CCF {ccf} (configured 1.0 / 0.3): {l}"
        );
        // Ground-truth relationship: measured AF ≈ dosage × NEAT_CCF. Average the
        // absolute error over adequately-covered sites to stay robust to per-site
        // binomial noise (the assertion is on the aggregate, not each record).
        if let (Some(d), Some(dp)) = (dosage_of(l), dp_of(l))
            && dp >= 25
        {
            let sample = l.split('\t').next_back().unwrap();
            let af: f64 = sample.split(':').next_back().unwrap().parse().unwrap();
            err_sum += (af - d * ccf).abs();
            checked += 1;
        }
    }
    assert!(
        checked >= 10,
        "too few well-covered somatic sites to validate ({checked})"
    );
    let mean_err = err_sum / checked as f64;
    assert!(
        mean_err < 0.05,
        "measured AF should track dosage × NEAT_CCF; mean |err| = {mean_err:.4} over {checked} sites"
    );
}
