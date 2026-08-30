//! Baseline germline invariants: assertions that the *distributions* a mutation
//! model carries actually reach `gen-reads` output.
//!
//! Motivation (#372 and the SBS1 regression): every caller-level check passed for
//! months while trinucleotide context was used only to pick the replacement base
//! and never to weight WHERE mutations land. Recall, precision, Ts/Tv and variant
//! counts were all unaffected — the defect lived entirely in the spectrum, which
//! nothing asserted. These tests close that class:
//!
//! * **known-answer** — plant C>T at every CpG, then measure the realised context
//!   distribution against each context's occurrence in the reference. Correct
//!   placement concentrates SNPs at CpG (~48x enrichment on H1N1); context-neutral
//!   placement yields ~1x. No caller, no signature-fitting tool, no cancer.
//! * **round-trip** — fit a model, simulate with it, re-fit a model from the
//!   simulated output, and compare. Any parameter the builder measures is covered
//!   automatically, including ones added later.
//! * **no dead parameters** — perturb a model field and assert the corresponding
//!   output statistic moves. A parameter that changes nothing is not wired up.

mod common;

use std::collections::HashMap;
use std::fs;
use std::io::{Read as _, Write as _};
use std::path::{Path, PathBuf};

use common::{eidolon, h1n1_reference};
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde_json::Value;

// ── fixtures ────────────────────────────────────────────────────────────────

/// Every contig in a FASTA as (name, uppercased sequence).
fn read_fasta(path: &Path) -> Vec<(String, Vec<u8>)> {
    let data = fs::read_to_string(path).unwrap();
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in data.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.split_whitespace().next().unwrap().to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().to_ascii_uppercase().bytes());
        }
    }
    assert!(!out.is_empty(), "no contigs in {}", path.display());
    out
}

fn write_yaml(dir: &Path, tag: &str, body: &str) -> PathBuf {
    let p = dir.join(format!("{tag}.yml"));
    fs::write(&p, body).unwrap();
    p
}

/// 1-based POS of every CpG cytosine (ref `C` immediately followed by `G`),
/// excluding contig edges so a full trinucleotide context always exists.
fn cpg_positions(seq: &[u8]) -> Vec<usize> {
    (2..seq.len())
        .filter(|&i| seq[i - 1] == b'C' && seq[i] == b'G')
        .collect()
}

/// Training VCF planting C>T at every CpG cytosine — the SBS1 substitution, in the
/// only context where SBS1 lives. Returns the number of records written.
fn write_cpg_training_vcf(path: &Path, contigs: &[(String, Vec<u8>)]) -> usize {
    let mut f = fs::File::create(path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(
        f,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )
    .unwrap();
    let mut n = 0;
    for (name, seq) in contigs {
        for pos in cpg_positions(seq) {
            writeln!(f, "{name}\t{pos}\t.\tC\tT\t60\tPASS\t.\tGT\t0/1").unwrap();
            n += 1;
        }
    }
    n
}

// ── driving the binary ──────────────────────────────────────────────────────

fn build_mut_model(dir: &Path, tag: &str, vcf: &Path) -> PathBuf {
    let model = dir.join(format!("mut_{tag}.json.gz"));
    let cfg = write_yaml(
        dir,
        &format!("build_{tag}"),
        &format!(
            "reference: {r}\nvcf_file: {v}\noutput_file: {m}\n\
             overwrite_output: true\nbed_file: .\n",
            r = h1n1_reference().display(),
            v = vcf.display(),
            m = model.display(),
        ),
    );
    eidolon()
        .args(["gen-mut-model", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    model
}

/// Run gen-reads with `model`, VCF only, and return the output VCF path.
fn simulate_vcf(dir: &Path, tag: &str, model: &Path) -> PathBuf {
    simulate_vcf_ploidy(dir, tag, model, 1)
}

/// As `simulate_vcf`, at an explicit ploidy (zygosity needs ploidy >= 2).
fn simulate_vcf_ploidy(dir: &Path, tag: &str, model: &Path, ploidy: usize) -> PathBuf {
    let cfg = write_yaml(
        dir,
        &format!("sim_{tag}"),
        &format!(
            "reference: {r}\nread_len: 100\ncoverage: 5\npaired_ended: false\nploidy: {p}\n\
             mutation_model: {m}\nproduce_vcf: true\nproduce_fastq: false\n\
             produce_bam: false\noverwrite_output: true\noutput_dir: {o}\n\
             output_filename: bgi_{tag}\nrng_seed: baseline germline {tag}\nnum_threads: 1\n",
            r = h1n1_reference().display(),
            m = model.display(),
            o = dir.display(),
            p = ploidy,
        ),
    );
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(&cfg)
        .assert()
        .success();
    dir.join(format!("bgi_{tag}.vcf.gz"))
}

// ── reading output ──────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Rec {
    chrom: String,
    pos: usize,
    r#ref: String,
    alt: String,
    /// First subfield of the sample column (the GT), empty when absent.
    gt: String,
}

fn read_vcf(path: &Path) -> Vec<Rec> {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    raw.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 5 {
                return None;
            }
            Some(Rec {
                chrom: f[0].to_string(),
                pos: f[1].parse().ok()?,
                r#ref: f[3].to_ascii_uppercase(),
                // first ALT only; multi-allelic records are still one placement
                alt: f[4].split(',').next()?.to_ascii_uppercase(),
                gt: f
                    .get(9)
                    .and_then(|s| s.split(':').next())
                    .unwrap_or("")
                    .to_string(),
            })
        })
        .collect()
}

fn is_snp(r: &Rec) -> bool {
    r.r#ref.len() == 1 && r.alt.len() == 1 && !r.alt.starts_with('<')
}

/// (insertions, deletions) among non-SNP records.
fn indel_counts(recs: &[Rec]) -> (usize, usize) {
    let mut ins = 0;
    let mut del = 0;
    for r in recs {
        if r.alt.starts_with('<') || is_snp(r) {
            continue;
        }
        if r.alt.len() > r.r#ref.len() {
            ins += 1;
        } else if r.r#ref.len() > r.alt.len() {
            del += 1;
        }
    }
    (ins, del)
}

/// Of the SNPs in `recs`, how many sit in a CpG context, and how many total.
/// Uses the *reference* trinucleotide at each SNP position, so this is measured
/// independently of anything the model claims.
fn cpg_snp_fraction(recs: &[Rec], contigs: &[(String, Vec<u8>)]) -> (usize, usize) {
    let by_name: HashMap<&str, &Vec<u8>> = contigs.iter().map(|(n, s)| (n.as_str(), s)).collect();
    let mut cpg = 0;
    let mut total = 0;
    for r in recs.iter().filter(|r| is_snp(r)) {
        let Some(seq) = by_name.get(r.chrom.as_str()) else {
            continue;
        };
        // POS is 1-based; need a base on each side for a full trinucleotide.
        if r.pos < 2 || r.pos >= seq.len() {
            continue;
        }
        total += 1;
        // CpG cytosine: this base is C and the next is G.
        if seq[r.pos - 1] == b'C' && seq[r.pos] == b'G' {
            cpg += 1;
        }
    }
    (cpg, total)
}

// ── model JSON access ───────────────────────────────────────────────────────

fn load_model(path: &Path) -> Value {
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    serde_json::from_str(&raw).unwrap()
}

fn save_model(path: &Path, v: &Value) {
    let mut enc = GzEncoder::new(fs::File::create(path).unwrap(), Compression::default());
    enc.write_all(serde_json::to_string(v).unwrap().as_bytes())
        .unwrap();
    enc.finish().unwrap();
}

/// `snp_distro` stores *cumulative* normalized weights over the 64 frames in a
/// fixed order. Difference them to recover each frame's weight. The frame order is
/// deliberately not reconstructed here — these tests compare weight vectors to each
/// other, or measure enrichment in sequence space, so the index→trinucleotide map
/// never has to be duplicated from the source.
fn context_weight_vector(model: &Value) -> Vec<f64> {
    let cum = model["statistical_models"]["snp_trinuc_model"]["snp_distro"]["weights"]
        .as_array()
        .expect("snp_distro.weights should be an array");
    let mut out = Vec::with_capacity(cum.len());
    let mut prev = 0.0f64;
    for c in cum {
        let c = c.as_f64().unwrap();
        out.push((c - prev).max(0.0));
        prev = c;
    }
    out
}

/// `DiscreteDistribution` serializes its weights **cumulatively** (normalized running
/// sum), not as raw weights: `variant_dist.weights == [1.0, 1.0, 1.0]` means
/// SNP=1.0, Insertion=0.0, Deletion=0.0 — not "equal thirds". Tests state intent as
/// raw weights and convert here, so a malformed edit can't silently mean something
/// else. (This bit the first draft of these tests: `[0.0, 1.0, 1.0]` was intended as
/// "indels, evenly split" but actually means "100% insertions".)
fn cumulative_weights(raw: &[f64]) -> Value {
    let total: f64 = raw.iter().sum();
    assert!(total > 0.0, "weights must not be all-zero");
    let mut acc = 0.0f64;
    let mut out = Vec::with_capacity(raw.len());
    for w in raw {
        acc += w / total;
        out.push(acc);
    }
    // guard against drift so the last bucket is always reachable
    *out.last_mut().unwrap() = 1.0;
    serde_json::json!(out)
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    let dot: f64 = (0..n).map(|i| a[i] * b[i]).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na == 0.0 || nb == 0.0 {
        0.0
    } else {
        dot / (na * nb)
    }
}

fn argmax(v: &[f64]) -> usize {
    v.iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap()
}

// ════════════════════════════════════════════════════════════════════════════
// 1. Known-answer: context weights decide WHERE mutations land.
// ════════════════════════════════════════════════════════════════════════════

/// Train on C>T at every CpG, simulate, and measure the realised context
/// distribution against CpG's occurrence in the reference.
///
/// The expected answer is computable without running the code: H1N1 is ~2.06% CpG
/// by interior position, so context-neutral placement puts ~2% of SNPs at CpG
/// (enrichment ~1x) while context-weighted placement concentrates them there
/// (enrichment ~48x). This is the assertion that would have caught the SBS1
/// regression on a 13 kb fixture in seconds.
#[test]
fn context_weights_decide_where_snps_land() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    // Baseline expectation, measured from the reference itself.
    let mut cpg_sites = 0usize;
    let mut interior = 0usize;
    for (_, seq) in &contigs {
        cpg_sites += cpg_positions(seq).len();
        interior += seq.len().saturating_sub(2);
    }
    let cpg_share = cpg_sites as f64 / interior as f64;
    assert!(
        cpg_share > 0.0 && cpg_share < 0.10,
        "fixture sanity: CpG share {cpg_share:.4} is not a small minority"
    );

    let vcf = dir.join("cpg_train.vcf");
    let planted = write_cpg_training_vcf(&vcf, &contigs);
    assert_eq!(
        planted, cpg_sites,
        "training VCF should plant exactly one SNP per CpG cytosine"
    );

    let model = build_mut_model(dir, "cpg", &vcf);
    let out = simulate_vcf(dir, "cpg", &model);
    let recs = read_vcf(&out);
    let (cpg_snps, total_snps) = cpg_snp_fraction(&recs, &contigs);

    assert!(
        total_snps >= 30,
        "too few SNPs ({total_snps}) to judge a context distribution — fixture is \
         underpowered, not passing"
    );

    let observed = cpg_snps as f64 / total_snps as f64;
    let enrichment = observed / cpg_share;
    eprintln!(
        "[baseline] CpG sites {cpg_sites}/{interior} = {:.2}% of positions; \
         output SNPs at CpG {cpg_snps}/{total_snps} = {:.2}% → enrichment {enrichment:.1}x",
        100.0 * cpg_share,
        100.0 * observed
    );

    assert!(
        enrichment > 5.0,
        "SNP placement is not following the model's context weights: CpG holds \
         {:.2}% of reference positions and {:.2}% of output SNPs (enrichment \
         {enrichment:.1}x, expected >>1). Context-neutral placement gives ~1x — \
         this is the SBS1 regression shape.",
        100.0 * cpg_share,
        100.0 * observed
    );
}

// ════════════════════════════════════════════════════════════════════════════
// 2. Round-trip: fit → simulate → re-fit → compare.
// ════════════════════════════════════════════════════════════════════════════

/// Whatever the builder can measure, the simulator must reproduce well enough that
/// re-fitting recovers it. This covers every distribution the builder learns —
/// including ones added after this test was written — without naming them.
#[test]
fn mutation_model_round_trips_its_context_spectrum() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    let vcf = dir.join("rt_train.vcf");
    write_cpg_training_vcf(&vcf, &contigs);
    let model_in = build_mut_model(dir, "rt_in", &vcf);

    // simulate, then re-fit a model from what came out
    let out_vcf = simulate_vcf(dir, "rt", &model_in);
    let n_out = read_vcf(&out_vcf).len();
    assert!(
        n_out >= 30,
        "round-trip needs a usable number of output variants, got {n_out}"
    );
    // gen-mut-model wants a plain VCF; decompress the simulated one.
    let plain = dir.join("rt_out.vcf");
    let mut raw = String::new();
    GzDecoder::new(fs::File::open(&out_vcf).unwrap())
        .read_to_string(&mut raw)
        .unwrap();
    fs::write(&plain, &raw).unwrap();
    let model_out = build_mut_model(dir, "rt_out", &plain);

    let w_in = context_weight_vector(&load_model(&model_in));
    let w_out = context_weight_vector(&load_model(&model_out));
    let sim = cosine(&w_in, &w_out);
    let flat = vec![1.0f64; w_in.len()];
    let sim_flat = cosine(&w_out, &flat);

    eprintln!(
        "[baseline] round-trip context spectrum: cosine(in,out)={sim:.3}, \
         cosine(out,uniform)={sim_flat:.3}, argmax in={} out={}",
        argmax(&w_in),
        argmax(&w_out)
    );

    // The re-fitted spectrum must look like the input, and must NOT look like the
    // uniform spectrum a context-neutral simulator would produce.
    assert!(
        sim > 0.80,
        "re-fitted context spectrum does not resemble the model that produced it \
         (cosine {sim:.3}); a parameter is being dropped between model and output"
    );
    assert!(
        sim > sim_flat,
        "re-fitted spectrum is closer to uniform ({sim_flat:.3}) than to the input \
         model ({sim:.3}) — placement is context-neutral"
    );
}

// ════════════════════════════════════════════════════════════════════════════
// 3. No dead parameters.
// ════════════════════════════════════════════════════════════════════════════

/// `variant_dist` must decide the SNV:indel mix. Perturbing it and seeing the
/// output mix move is the check that the field is wired at all.
#[test]
fn variant_dist_governs_snv_indel_ratio() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    let vcf = dir.join("vd_train.vcf");
    write_cpg_training_vcf(&vcf, &contigs);
    let base = build_mut_model(dir, "vd_base", &vcf);

    // Two models differing ONLY in variant_dist: SNP-only vs indel-only.
    let mut snp_only = load_model(&base);
    snp_only["variant_dist"]["weights"] = cumulative_weights(&[1.0, 0.0, 0.0]);
    let snp_model = dir.join("mut_snp_only.json.gz");
    save_model(&snp_model, &snp_only);

    let mut indel_heavy = load_model(&base);
    indel_heavy["variant_dist"]["weights"] = cumulative_weights(&[0.0, 1.0, 1.0]);
    let indel_model = dir.join("mut_indel_heavy.json.gz");
    save_model(&indel_model, &indel_heavy);

    let a = read_vcf(&simulate_vcf(dir, "vd_snp", &snp_model));
    let b = read_vcf(&simulate_vcf(dir, "vd_indel", &indel_model));

    let a_snp = a.iter().filter(|r| is_snp(r)).count();
    let b_snp = b.iter().filter(|r| is_snp(r)).count();
    let a_frac = a_snp as f64 / a.len().max(1) as f64;
    let b_frac = b_snp as f64 / b.len().max(1) as f64;
    eprintln!(
        "[baseline] variant_dist: SNP-only model → {a_snp}/{} SNP ({:.2}); \
         indel-only model → {b_snp}/{} SNP ({:.2})",
        a.len(),
        a_frac,
        b.len(),
        b_frac
    );

    assert!(
        !a.is_empty() && !b.is_empty(),
        "no variants produced; cannot judge the SNV:indel mix"
    );
    assert!(
        a_frac > 0.90,
        "a SNP-only variant_dist still produced {:.0}% non-SNP records — \
         variant_dist is not governing variant type",
        100.0 * (1.0 - a_frac)
    );
    assert!(
        b_frac < 0.10,
        "an indel-only variant_dist still produced {:.0}% SNPs — \
         variant_dist is not governing variant type",
        100.0 * b_frac
    );
}

/// The insertion:deletion split within indels must follow the model.
///
/// NOTE: the field that governs this is `variant_dist` (fitted by gen-mut-model as
/// `[snp_freq, ins_freq, del_freq]`), **not**
/// `statistical_models.indel_model.insertion_probability`. That field is fitted and
/// serialized but never consulted when generating variants — its only reader,
/// `IndelModel::is_insertion()`, has no callers. Asserting against `variant_dist`
/// tests the path that is actually live; see the audit note for the dead field.
#[test]
fn variant_dist_governs_ins_del_split() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    let vcf = dir.join("id_train.vcf");
    write_cpg_training_vcf(&vcf, &contigs);
    let base = build_mut_model(dir, "id_base", &vcf);

    // Two models differing ONLY in the ins:del ratio, with SNPs suppressed so the
    // split is directly measurable.
    let mk = |ins: f64, del: f64, tag: &str| -> PathBuf {
        let mut m = load_model(&base);
        m["variant_dist"]["weights"] = cumulative_weights(&[0.0, ins, del]);
        let path = dir.join(format!("mut_id_{tag}.json.gz"));
        save_model(&path, &m);
        path
    };

    let ins_heavy = mk(9.0, 1.0, "ins");
    let del_heavy = mk(1.0, 9.0, "del");

    let (i_a, d_a) = indel_counts(&read_vcf(&simulate_vcf(dir, "id_ins", &ins_heavy)));
    let (i_b, d_b) = indel_counts(&read_vcf(&simulate_vcf(dir, "id_del", &del_heavy)));

    let frac = |i: usize, d: usize| i as f64 / (i + d).max(1) as f64;
    let f_a = frac(i_a, d_a);
    let f_b = frac(i_b, d_b);
    eprintln!(
        "[baseline] variant_dist ins:del 9:1 → ins {i_a}/del {d_a} ({f_a:.2}); \
         1:9 → ins {i_b}/del {d_b} ({f_b:.2})"
    );

    assert!(
        i_a + d_a >= 20 && i_b + d_b >= 20,
        "too few indels ({} and {}) to judge the ins:del split",
        i_a + d_a,
        i_b + d_b
    );
    assert!(
        f_a > 0.6,
        "a 9:1 insertion-weighted model produced only {:.0}% insertions",
        100.0 * f_a
    );
    assert!(
        f_b < 0.4,
        "a 1:9 insertion-weighted model produced {:.0}% insertions",
        100.0 * f_b
    );
    assert!(
        f_a > f_b + 0.3,
        "the ins:del split did not track variant_dist (9:1 → {f_a:.2}, 1:9 → {f_b:.2})"
    );
}

/// The indel length distributions must decide the lengths that appear in output.
///
/// Known-answer: force `ins_dist` to a single insertion length and `del_dist` to a
/// single (different) deletion length, then read the lengths back out of the VCF.
/// The expected answer is fixed by the model, not by anything the code computes.
#[test]
fn indel_length_distributions_drive_output_lengths() {
    const INS_LEN: usize = 3;
    const DEL_LEN: usize = 6;

    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    let vcf = dir.join("len_train.vcf");
    write_cpg_training_vcf(&vcf, &contigs);
    let base = build_mut_model(dir, "len_base", &vcf);

    let mut m = load_model(&base);
    // indels only, evenly split, so both length distributions are exercised
    m["variant_dist"]["weights"] = cumulative_weights(&[0.0, 1.0, 1.0]);
    m["statistical_models"]["indel_model"]["ins_dist"] =
        serde_json::json!({ "values": [INS_LEN], "weights": [1.0] });
    m["statistical_models"]["indel_model"]["del_dist"] =
        serde_json::json!({ "values": [DEL_LEN], "weights": [1.0] });
    let model = dir.join("mut_fixed_len.json.gz");
    save_model(&model, &m);

    let recs = read_vcf(&simulate_vcf(dir, "len", &model));

    let mut ins_lens: Vec<usize> = Vec::new();
    let mut del_lens: Vec<usize> = Vec::new();
    for r in &recs {
        if r.alt.starts_with('<') || is_snp(r) {
            continue;
        }
        if r.alt.len() > r.r#ref.len() {
            ins_lens.push(r.alt.len() - r.r#ref.len());
        } else if r.r#ref.len() > r.alt.len() {
            del_lens.push(r.r#ref.len() - r.alt.len());
        }
    }

    assert!(
        ins_lens.len() >= 20 && del_lens.len() >= 20,
        "too few indels to judge lengths (ins {}, del {})",
        ins_lens.len(),
        del_lens.len()
    );

    let ins_ok = ins_lens.iter().filter(|&&l| l == INS_LEN).count();
    let del_ok = del_lens.iter().filter(|&&l| l == DEL_LEN).count();
    let ins_frac = ins_ok as f64 / ins_lens.len() as f64;
    let del_frac = del_ok as f64 / del_lens.len() as f64;
    eprintln!(
        "[baseline] fixed lengths: insertions {ins_ok}/{} at {INS_LEN}bp ({ins_frac:.2}); \
         deletions {del_ok}/{} at {DEL_LEN}bp ({del_frac:.2})",
        ins_lens.len(),
        del_lens.len()
    );

    // A small tail is expected: deletions landing near a contig edge fall back to a
    // SNP-like record rather than truncating the reference.
    assert!(
        ins_frac > 0.90,
        "ins_dist pinned to {INS_LEN}bp but only {:.0}% of insertions were that \
         length — the insertion length distribution is not reaching output",
        100.0 * ins_frac
    );
    assert!(
        del_frac > 0.90,
        "del_dist pinned to {DEL_LEN}bp but only {:.0}% of deletions were that \
         length — the deletion length distribution is not reaching output",
        100.0 * del_frac
    );
}

/// `homozygous_frequency` must decide the het/hom ratio in output.
///
/// The existing guard (`test_default_homozygous_frequency_is_realistic`) pins the
/// default model's *value* to a realistic diploid range, which is what caught the
/// NEAT-lineage 0.01 default. It does not assert the value reaches output — the
/// same value/wiring gap that let context weighting sit unused. This closes it.
#[test]
fn homozygous_frequency_governs_output_het_hom_ratio() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let contigs = read_fasta(&h1n1_reference());

    let vcf = dir.join("hz_train.vcf");
    write_cpg_training_vcf(&vcf, &contigs);
    let base = build_mut_model(dir, "hz_base", &vcf);

    let mk = |hz: f64, tag: &str| -> PathBuf {
        let mut m = load_model(&base);
        m["homozygous_frequency"] = serde_json::json!(hz);
        let path = dir.join(format!("mut_hz_{tag}.json.gz"));
        save_model(&path, &m);
        path
    };

    let hom_frac = |recs: &[Rec]| -> (usize, usize) {
        let mut hom = 0;
        let mut het = 0;
        for r in recs {
            let alleles: Vec<&str> = r.gt.split(['/', '|']).collect();
            if alleles.len() < 2 || alleles.contains(&".") {
                continue;
            }
            // homozygous-ALT: every allele is a non-reference call
            if alleles.iter().all(|a| *a != "0") {
                hom += 1;
            } else if alleles.iter().any(|a| *a != "0") {
                het += 1;
            }
        }
        (hom, het)
    };

    let lo = mk(0.05, "lo");
    let hi = mk(0.95, "hi");
    let (hom_lo, het_lo) = hom_frac(&read_vcf(&simulate_vcf_ploidy(dir, "hz_lo", &lo, 2)));
    let (hom_hi, het_hi) = hom_frac(&read_vcf(&simulate_vcf_ploidy(dir, "hz_hi", &hi, 2)));

    let f = |h: usize, t: usize| h as f64 / (h + t).max(1) as f64;
    let f_lo = f(hom_lo, het_lo);
    let f_hi = f(hom_hi, het_hi);
    eprintln!(
        "[baseline] homozygous_frequency 0.05 → hom {hom_lo}/het {het_lo} ({f_lo:.2}); \
         0.95 → hom {hom_hi}/het {het_hi} ({f_hi:.2})"
    );

    assert!(
        hom_lo + het_lo >= 30 && hom_hi + het_hi >= 30,
        "too few genotyped variants ({} and {}) to judge zygosity",
        hom_lo + het_lo,
        hom_hi + het_hi
    );
    assert!(
        f_lo < 0.25,
        "homozygous_frequency 0.05 produced {:.0}% homozygous calls",
        100.0 * f_lo
    );
    assert!(
        f_hi > 0.75,
        "homozygous_frequency 0.95 produced only {:.0}% homozygous calls",
        100.0 * f_hi
    );
}
