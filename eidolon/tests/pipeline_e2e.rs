//! End-to-end pipeline tests: train a model with one eidolon subcommand, then consume it
//! from another. These exercise the real binary boundary (CLI → config parser →
//! model file → consumer) — the kind of wiring breakage that unit tests will miss
//! because they call internal functions directly.

mod common;

use common::{
    GenReadsConfig, eidolon, fresh_workdir, h1n1_reference, read_gzip_fastq_lines,
    write_gen_seq_error_model_config, write_iupac_fasta, write_n_tract_fasta,
    write_synthetic_fastq_gz,
};
use std::collections::HashSet;

/// Convert Phred-encoded quality bytes back to raw Q-scores under Phred+33.
fn qual_bytes_to_scores(line: &str) -> Vec<u8> {
    line.bytes().map(|b| b.saturating_sub(33)).collect()
}

/// Extract every quality line from a multi-stream gzipped FASTQ. Quality lines are the
/// 4th of every 4-line record group.
fn quality_lines(fastq_path: &std::path::Path) -> Vec<String> {
    let lines = read_gzip_fastq_lines(fastq_path);
    assert!(
        lines.len().is_multiple_of(4),
        "FASTQ line count must be a multiple of 4, got {}",
        lines.len(),
    );
    lines
        .into_iter()
        .enumerate()
        .filter_map(|(i, l)| if i % 4 == 3 { Some(l) } else { None })
        .collect()
}

#[test]
fn gen_reads_output_is_byte_identical_across_thread_counts() {
    // Regression guard for the sub-contig chunking path: with a small chunk_size
    // the multi-segment H1N1 reference splits into several chunks per contig.
    // Output (R1 AND R2) must be byte-identical regardless of num_threads, and
    // free of duplicate read names. This catches the class of bug where a chunk
    // wrote to a shared temp file / a path was merged once per chunk, which
    // duplicated reads and raced under parallelism.
    let run = |threads: usize, prefix: &str| -> (std::path::PathBuf, std::path::PathBuf) {
        let (dir, out) = fresh_workdir();
        std::mem::forget(dir); // keep the temp dir alive for the duration of the test
        let mut config = GenReadsConfig::new(h1n1_reference(), out.clone(), prefix);
        config.paired_ended = true;
        config.coverage = 20;
        config.read_len = 50;
        config.rng_seed = "chunk determinism seed".to_string();
        config.num_threads = Some(threads);
        config.chunk_size = Some(500); // force several chunks per ~1.7 kb segment
        let yaml = config.write_yaml();
        eidolon()
            .args(["gen-reads", "-c"])
            .arg(yaml.path())
            .assert()
            .success();
        (
            out.join(format!("{prefix}_r1.fastq.gz")),
            out.join(format!("{prefix}_r2.fastq.gz")),
        )
    };

    let (r1_t1, r2_t1) = run(1, "t1");
    let (r1_t4, r2_t4) = run(4, "t4");

    for (a, b, which) in [(&r1_t1, &r1_t4, "R1"), (&r2_t1, &r2_t4, "R2")] {
        let la = read_gzip_fastq_lines(a);
        let lb = read_gzip_fastq_lines(b);
        assert!(!la.is_empty(), "{which}: no reads produced");
        assert_eq!(
            la, lb,
            "{which} differs between 1-thread and 4-thread runs (same seed) — chunking is not thread-count-invariant",
        );
        // No duplicate read names (the duplication symptom of the temp-file bug).
        let names: Vec<&String> = la.iter().step_by(4).collect();
        let unique: HashSet<&&String> = names.iter().collect();
        assert_eq!(
            names.len(),
            unique.len(),
            "{which}: duplicate read names present",
        );
    }
}

#[test]
fn gen_reads_with_default_model_produces_well_formed_fastq() {
    // Smoke test: run gen-reads end-to-end through the real binary with all defaults
    // (no custom error model, no input VCF) and verify the output is well-formed
    // (multiple of 4 lines, seq length == qual length, names start with '@').
    let (_dir, out) = fresh_workdir();
    let config = GenReadsConfig::new(h1n1_reference(), out.clone(), "default_model");
    let yaml = config.write_yaml();

    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let fastq = out.join("default_model_r1.fastq.gz");
    assert!(fastq.exists(), "FASTQ not produced at {fastq:?}");

    let lines = read_gzip_fastq_lines(&fastq);
    assert!(!lines.is_empty(), "FASTQ has no records");
    assert_eq!(lines.len() % 4, 0, "incomplete record at FASTQ tail");

    for chunk in lines.chunks(4) {
        assert!(
            chunk[0].starts_with('@'),
            "name line must start with @: {:?}",
            chunk[0]
        );
        assert_eq!(chunk[2], "+", "third line must be '+': {:?}", chunk[2]);
        assert_eq!(
            chunk[1].len(),
            chunk[3].len(),
            "seq/qual length mismatch in {chunk:?}",
        );
    }
}

#[test]
fn train_binned_model_then_gen_reads_emits_only_bin_qualities() {
    // Cross-binary integration: train a binned quality model with gen-seq-error-model,
    // hand it to gen-reads, then prove gen-reads honors the binning by inspecting the
    // FASTQ output. This is the user-visible promise of binned-quality support.
    let (_dir, work) = fresh_workdir();

    // Step 1: build a tiny FASTQ for training. Quality scores cycle 33..=42 — all of
    // which snap to bin 37 under the bin set we'll request.
    let train_fastq = work.join("train.fastq.gz");
    write_synthetic_fastq_gz(&train_fastq, 100, 50);

    // Step 2: train a binned sequencing error model via the real binary.
    let model_path = work.join("binned_model.json.gz");
    let bins = [2usize, 12, 23, 37];
    let model_cfg = write_gen_seq_error_model_config(&train_fastq, &model_path, Some(&bins));
    eidolon()
        .args(["gen-seq-error-model", "-c"])
        .arg(model_cfg.path())
        .assert()
        .success();
    assert!(model_path.exists(), "model not produced at {model_path:?}");

    // Step 3: generate reads using the binned model.
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "binned_pipeline");
    config.sequence_error_model = Some(model_path);
    let reads_cfg = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(reads_cfg.path())
        .assert()
        .success();

    // Step 4: verify every emitted quality score is in the bin set.
    let fastq = work.join("binned_pipeline_r1.fastq.gz");
    let bin_set: HashSet<u8> = bins.iter().map(|&b| b as u8).collect();
    let mut total_qualities = 0usize;
    for ql in quality_lines(&fastq) {
        for score in qual_bytes_to_scores(&ql) {
            assert!(
                bin_set.contains(&score),
                "gen-reads emitted non-bin quality score {score} (bins = {bins:?})",
            );
            total_qualities += 1;
        }
    }
    assert!(total_qualities > 0, "no quality scores observed");
}

#[test]
fn gen_reads_with_symbolic_del_modulates_depth_and_round_trips_to_vcf() {
    // End-to-end contract: a symbolic <DEL> in the input VCF must (a) be
    // preserved verbatim in the output golden VCF, and (b) drive per-region
    // depth modulation in gen_reads. We park a homozygous <DEL> at
    // H1N1_HA:500-1500 (1-based, inclusive), which corresponds to 0-based
    // half-open [499, 1500). At a hom DEL the multiplier is 0, so zero reads
    // should start in that span; reads outside it should still be generated
    // at full coverage.
    use std::io::Write as _;

    let (_dir, work) = fresh_workdir();

    // Coverage modulation needs the SV's span — supply it via INFO/END.
    let input_vcf = work.join("input_sv.vcf");
    {
        let mut f = std::fs::File::create(&input_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "##ALT=<ID=DEL,Description=\"Deletion\">").unwrap();
        writeln!(
            f,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
        )
        .unwrap();
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
        // Homozygous <DEL> at H1N1_HA, 1-based POS=500, END=1500.
        writeln!(f, "H1N1_HA\t500\t.\tA\t<DEL>\t60\tPASS\tEND=1500\tGT\t1/1").unwrap();
    }

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "sv_run");
    config.coverage = 50;
    config.produce_vcf = true;
    config.input_vcf = Some(input_vcf);
    // Drive enough de-novo mutations across the contig that, pre-fix, several
    // would have been sampled inside the DEL span. Post-fix: zero — the rate
    // is zeroed over hom-DEL spans before generate_variants runs.
    config.mutation_rate = Some(0.05);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    // (a) Output VCF round-trip — the <DEL> should appear verbatim, with the
    // original INFO field preserved (END=1500).
    let out_vcf_gz = work.join("sv_run.vcf.gz");
    assert!(
        out_vcf_gz.exists(),
        "output VCF not produced: {out_vcf_gz:?}"
    );
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_vcf_gz).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let del_line = vcf_lines
        .iter()
        .find(|l| l.starts_with("H1N1_HA\t500\t") && l.contains("<DEL>"))
        .unwrap_or_else(|| panic!("expected <DEL> record in output VCF; got: {vcf_lines:?}"));
    assert!(
        del_line.contains("END=1500"),
        "expected INFO/END to round-trip; got: {del_line}"
    );
    // No de-novo SNPs should be sampled inside the deleted span: the mutation
    // rate over [POS+1, END] in 1-based (= [501, 1500]) is zeroed for hom DEL.
    let snps_inside_del: Vec<&String> = vcf_lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.starts_with("H1N1_HA\t"))
        .filter(|l| !l.contains("<DEL>"))
        .filter(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            cols.get(1)
                .and_then(|p| p.parse::<usize>().ok())
                .map(|p| (501..=1500).contains(&p))
                .unwrap_or(false)
        })
        .collect();
    assert!(
        snps_inside_del.is_empty(),
        "expected no de-novo SNPs inside the hom <DEL> span [501, 1500]; got: {snps_inside_del:?}"
    );

    // (b) Depth modulation — count read start positions on H1N1_HA. Names look
    // like `@EIDOLON_generated_{contig}_{abs_start:010}_{abs_end:010}`.
    let fastq = work.join("sv_run_r1.fastq.gz");
    let lines = read_gzip_fastq_lines(&fastq);
    let mut inside_del = 0usize;
    let mut outside_del = 0usize;
    for chunk in lines.chunks(4) {
        let name = &chunk[0];
        if !name.starts_with("@EIDOLON_generated_H1N1_HA_") {
            continue;
        }
        // Strip leading '@' and prefix; the remainder is "{start}_{end}".
        let rest = name.trim_start_matches("@EIDOLON_generated_H1N1_HA_");
        let abs_start: usize = rest
            .split('_')
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or_else(|| panic!("unparseable read name: {name}"));
        // DEL span (0-based half-open): [499, 1500). Use inset margins so a
        // single-base off-by-one in coverage modulation wouldn't pass.
        if (500..1450).contains(&abs_start) {
            inside_del += 1;
        } else if abs_start < 449 || abs_start >= 1550 {
            outside_del += 1;
        }
    }
    assert_eq!(
        inside_del, 0,
        "expected ZERO reads starting inside the hom <DEL> span; got {inside_del}"
    );
    assert!(
        outside_del > 200,
        "expected plenty of reads outside the DEL span at coverage=50 over ~675 unaffected \
         bases of H1N1_HA; got only {outside_del}"
    );
}

#[test]
fn gen_reads_with_iupac_reference_produces_no_iupac_in_output() {
    // Regression test: a reference containing IUPAC ambiguity codes (R/Y/M/K/S/W/H/B/V/D)
    // must be processed without crashing and must not leak any IUPAC letter into output
    // FASTQ sequence lines — every base emitted must be in [ACGTNacgtn].
    let (_dir, work) = fresh_workdir();
    let ref_path = work.join("iupac_ref.fa");
    write_iupac_fasta(&ref_path);

    let config = GenReadsConfig::new(ref_path, work.clone(), "iupac_run");
    let yaml = config.write_yaml();

    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let fastq = work.join("iupac_run_r1.fastq.gz");
    assert!(fastq.exists(), "FASTQ not produced at {fastq:?}");

    let lines = read_gzip_fastq_lines(&fastq);
    assert!(!lines.is_empty(), "FASTQ has no records");
    assert_eq!(lines.len() % 4, 0);

    let iupac: HashSet<char> = "RYMKSWHBVDrymkswhbvd".chars().collect();
    for chunk in lines.chunks(4) {
        for c in chunk[1].chars() {
            assert!(
                !iupac.contains(&c),
                "IUPAC character {:?} leaked into read sequence: {:?}",
                c,
                chunk[1]
            );
        }
    }
}

#[test]
fn gen_reads_treats_n_tract_as_gap_no_variants_inside() {
    // Regression for the N-handling fix: reference N bases are assembly gaps, not
    // sequence to fill. A central N tract must be excluded from variant placement
    // by the non-N region machinery. Pre-fix, `apply_n_substitution` overwrote every
    // N with a random base at load time, so the entire contig became one NonNRegion
    // and de-novo variants were sampled across the gap. The reference has flanks at
    // 1-based [1, 240] and [481, 720] with a hard N tract at [241, 480].
    let (_dir, work) = fresh_workdir();
    let ref_path = work.join("n_tract_ref.fa");
    write_n_tract_fasta(&ref_path);

    let mut config = GenReadsConfig::new(ref_path, work.clone(), "n_tract_run");
    config.coverage = 50;
    config.produce_vcf = true;
    // High rate so, pre-fix, many variants would have landed in the 240 bp gap.
    config.mutation_rate = Some(0.05);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_vcf_gz = work.join("n_tract_run.vcf.gz");
    assert!(
        out_vcf_gz.exists(),
        "output VCF not produced: {out_vcf_gz:?}"
    );
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(
            std::fs::File::open(&out_vcf_gz).unwrap(),
        ));
        r.lines().map(|l| l.unwrap()).collect()
    };

    let variants_in_gap: Vec<&String> = vcf_lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| {
            l.split('\t')
                .nth(1)
                .and_then(|p| p.parse::<usize>().ok())
                // 1-based N tract span [241, 480].
                .map(|p| (241..=480).contains(&p))
                .unwrap_or(false)
        })
        .collect();
    assert!(
        variants_in_gap.is_empty(),
        "variants were placed inside the N tract [241, 480]; N regions must be excluded: {variants_in_gap:?}"
    );
}

#[test]
fn gen_reads_with_trained_sv_model_emits_de_novo_symbolic_svs() {
    // Phase 4 e2e: train a custom SV-rich mutation model with
    // `gen-mut-model`, then run `gen-reads` against the H1N1 reference
    // with `sv_rate_scale: 10.0` and the trained model. The output VCF
    // must contain at least one de novo symbolic record carrying both
    // INFO/SVTYPE and INFO/END, and that record's POS must NOT match
    // any of the input training records (which aren't part of the
    // generation run at all — the trained model is its own input).
    use std::io::Write as _;

    let (_dir, work) = fresh_workdir();

    // Synthetic training VCF: multiple DELs and DUPs spread across the
    // H1N1 contigs. Lengths are small (60-180bp) so the resulting
    // log-normal fit produces SVs that comfortably fit inside the
    // largest H1N1 contig (HA at ~1700bp; sampler upper-bound is
    // contig_len / 4 ≈ 425bp). Also includes a literal SNP so the
    // SNP/indel side of the trainer math doesn't divide by zero
    // (caught by v1.9.1's NoLiteralVariants guard).
    let train_vcf = work.join("sv_train.vcf");
    {
        let mut f = std::fs::File::create(&train_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
        )
        .unwrap();
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
        // Literal SNP — keeps the SNP/indel denominator non-zero.
        writeln!(f, "H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1").unwrap();
        // Small DELs: 60, 80, 90, 120bp. Mean ~88bp → log-normal μ ≈ ln(88) ≈ 4.5.
        writeln!(f, "H1N1_HA\t100\t.\tA\t<DEL>\t60\tPASS\tEND=160\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t300\t.\tA\t<DEL>\t60\tPASS\tEND=380\tGT\t1/1").unwrap();
        writeln!(f, "H1N1_HA\t500\t.\tA\t<DEL>\t60\tPASS\tEND=590\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t700\t.\tA\t<DEL>\t60\tPASS\tEND=820\tGT\t0/1").unwrap();
        // Small DUPs: 80, 100, 130, 150bp.
        writeln!(f, "H1N1_HA\t900\t.\tT\t<DUP>\t60\tPASS\tEND=980\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t1100\t.\tT\t<DUP>\t60\tPASS\tEND=1200\tGT\t1/1").unwrap();
        writeln!(f, "H1N1_HA\t1300\t.\tT\t<DUP>\t60\tPASS\tEND=1430\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t1500\t.\tT\t<DUP>\t60\tPASS\tEND=1650\tGT\t0/1").unwrap();
    }

    // gen-mut-model needs a tiny YAML config pointing at the reference,
    // VCF, and output path.
    let trained_model = work.join("trained_sv.json.gz");
    let train_yaml = work.join("train.yml");
    std::fs::write(
        &train_yaml,
        format!(
            "reference: {ref_}\n\
             vcf_file: {vcf}\n\
             output_file: {out}\n\
             overwrite_output: true\n",
            ref_ = h1n1_reference().display(),
            vcf = train_vcf.display(),
            out = trained_model.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-mut-model", "-c"])
        .arg(&train_yaml)
        .assert()
        .success();
    assert!(trained_model.exists(), "gen-mut-model didn't write output");

    // Now drive gen-reads with the trained model and a high
    // sv_rate_scale so de novo sampling fires reliably.
    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "sv_sampled");
    config.coverage = 25;
    config.produce_vcf = true;
    config.mutation_model = Some(trained_model);
    config.sv_rate_scale = Some(10.0);
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    // Inspect the output VCF for symbolic records. At least one DEL or
    // DUP must appear, and each must carry both SVTYPE and END so the
    // round-trip downstream (compare-vcfs, depth modulation, etc.)
    // recognizes it.
    let out_vcf = work.join("sv_sampled.vcf.gz");
    assert!(out_vcf.exists(), "gen-reads did not produce output VCF");
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out_vcf).unwrap()));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let symbolic_records: Vec<&String> = vcf_lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.contains("<DEL>") || l.contains("<DUP>") || l.contains("<CNV>"))
        .collect();
    assert!(
        !symbolic_records.is_empty(),
        "expected at least one de novo symbolic SV in output; got VCF lines: {vcf_lines:#?}"
    );
    // Presence is not enough: `contains("END=")` is satisfied by ANY value, and shifting
    // every emitted END by +1000 passed the entire workspace. END is what downstream
    // consumers size the event from (depth modulation, compare-vcfs, truvari), so assert
    // it is INTERNALLY CONSISTENT with the record's own SVLEN and POS.
    let info_int = |line: &str, key: &str| -> Option<i64> {
        line.split('\t')
            .nth(7)?
            .split(';')
            .find_map(|f| f.strip_prefix(key).and_then(|v| v.parse::<i64>().ok()))
    };
    for line in &symbolic_records {
        let pos: i64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        let end = info_int(line, "END=")
            .unwrap_or_else(|| panic!("symbolic record has no parseable INFO/END: {line}"));
        let svlen = info_int(line, "SVLEN=")
            .unwrap_or_else(|| panic!("symbolic record has no parseable INFO/SVLEN: {line}"));
        assert!(end > pos, "INFO/END={end} is not after POS={pos}: {line}");
        // The span END implies must equal |SVLEN|, using the convention `place_sv`
        // documents: for <DEL> POS is the base BEFORE the event, so the affected range is
        // [POS+1, END] and END-POS == |SVLEN|; for <DUP>/<CNV>/<INV> "POS is the first
        // affected base", so the range is [POS, END] inclusive and END-POS+1 == |SVLEN|.
        // Whether that second convention should match VCF's is issue #474; this asserts
        // the record is self-consistent under whichever one applies, which is what a
        // downstream consumer needs and what shifting every END by +1000 broke.
        let svtype = line
            .split('\t')
            .nth(7)
            .and_then(|i| i.split(';').find_map(|f| f.strip_prefix("SVTYPE=")))
            .unwrap_or_else(|| panic!("symbolic record has no INFO/SVTYPE: {line}"));
        let implied_span = if svtype == "DEL" {
            end - pos
        } else {
            end - pos + 1
        };
        assert_eq!(
            implied_span,
            svlen.abs(),
            "{svtype} record: the span INFO/END implies ({implied_span}) disagrees with \
             |SVLEN| ({}) — a consumer sizing this event off END sees a different \
             variant than one sizing it off SVLEN: {line}",
            svlen.abs()
        );
    }
    for line in &symbolic_records {
        assert!(
            line.contains("SVTYPE="),
            "symbolic record missing INFO/SVTYPE: {line}"
        );
        assert!(
            line.contains("END="),
            "symbolic record missing INFO/END: {line}"
        );
    }
}

#[test]
fn gen_reads_with_sv_rate_scale_zero_emits_no_de_novo_svs() {
    // Inverse of the sampling test: the same trained model + the same
    // contig, but `sv_rate_scale: 0.0` (the default) must produce zero
    // symbolic records. Locks in the opt-in semantics — a v1.9 upgrader
    // who doesn't touch their YAML can't suddenly see SVs.
    use std::io::Write as _;

    let (_dir, work) = fresh_workdir();
    let train_vcf = work.join("svs.vcf");
    {
        let mut f = std::fs::File::create(&train_vcf).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
        )
        .unwrap();
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
        writeln!(f, "H1N1_HA\t22\t.\tC\tT\t60\tPASS\t.\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t100\t.\tA\t<DEL>\t60\tPASS\tEND=160\tGT\t0/1").unwrap();
        writeln!(f, "H1N1_HA\t300\t.\tA\t<DEL>\t60\tPASS\tEND=380\tGT\t1/1").unwrap();
    }
    let trained_model = work.join("m.json.gz");
    let train_yaml = work.join("train.yml");
    std::fs::write(
        &train_yaml,
        format!(
            "reference: {r}\nvcf_file: {v}\noutput_file: {o}\noverwrite_output: true\n",
            r = h1n1_reference().display(),
            v = train_vcf.display(),
            o = trained_model.display(),
        ),
    )
    .unwrap();
    eidolon()
        .args(["gen-mut-model", "-c"])
        .arg(&train_yaml)
        .assert()
        .success();

    let mut config = GenReadsConfig::new(h1n1_reference(), work.clone(), "no_svs");
    config.coverage = 25;
    config.produce_vcf = true;
    config.mutation_model = Some(trained_model);
    // Default behavior — sv_rate_scale left at None (= absent from YAML
    // = 0.0 in RunConfiguration::default()).
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let out_vcf = work.join("no_svs.vcf.gz");
    let vcf_lines: Vec<String> = {
        use flate2::read::MultiGzDecoder;
        use std::io::{BufRead, BufReader};
        let r = BufReader::new(MultiGzDecoder::new(std::fs::File::open(&out_vcf).unwrap()));
        r.lines().map(|l| l.unwrap()).collect()
    };
    let symbolic: Vec<&String> = vcf_lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.contains("<DEL>") || l.contains("<DUP>") || l.contains("<CNV>"))
        .collect();
    assert!(
        symbolic.is_empty(),
        "sv_rate_scale=0 must produce no symbolic SVs; got: {symbolic:?}"
    );
}

// ── Stage 2 of #189: input-VCF-driven aneuploidy depth-modulation ────────
//
// These tests pin that `build_coverage_multipliers` and `coverage_multiplier_for`
// produce the correct depth profile for whole-contig SVs (the proxy for whole-
// chromosome aneuploidy events). The math is contig-length-independent — if it
// works at H1N1_HA's 13,325 bp it works at chr1's 250 Mb because there are no
// contig-size-dependent branches in the path.
//
// All four tests use the same shape:
//   1. Run gen-reads with NO input VCF to establish a per-test baseline read
//      count (RNG is seeded, so this is reproducible).
//   2. Run gen-reads with an input VCF carrying a whole-contig <DEL> / <DUP> /
//      <CNV> at the genotype / CN we want to exercise.
//   3. Assert that the modulated count is in the expected ratio of baseline,
//      within ±15% tolerance to absorb statistical variation (σ/N ~ 1/√N ≈
//      0.9% at baseline ~13k reads; 15% leaves ample margin while still
//      catching any sign / off-by-one error in the multiplier arithmetic).

/// Count reads whose names start with `@EIDOLON_generated_H1N1_HA_` (i.e., were
/// emitted from the H1N1_HA contig). Other test fixtures might add reads from
/// other contigs to the same FASTQ; this filter keeps the count contig-local.
fn count_h1n1_reads(fastq: &std::path::Path) -> usize {
    let lines = read_gzip_fastq_lines(fastq);
    lines
        .chunks(4)
        .filter(|c| c[0].starts_with("@EIDOLON_generated_H1N1_HA_"))
        .count()
}

/// Write a one-record VCF carrying a single symbolic SV spanning the whole
/// H1N1_HA contig. `pos_1based` / `end_1based` use VCF convention (1-based
/// inclusive); for `<DEL>` the deleted span is `[pos+1, end]`, for `<DUP>` /
/// `<CNV>` the affected span is `[pos, end]`.
fn write_whole_contig_sv_vcf(
    work: &std::path::Path,
    pos_1based: usize,
    end_1based: usize,
    alt: &str,
    info_extras: &str,
    gt: &str,
) -> std::path::PathBuf {
    use std::io::Write as _;
    let path = work.join("input_sv.vcf");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##ALT=<ID=DEL,Description=\"Deletion\">").unwrap();
    writeln!(f, "##ALT=<ID=DUP,Description=\"Duplication\">").unwrap();
    writeln!(f, "##ALT=<ID=CNV,Description=\"Copy number variable\">").unwrap();
    writeln!(
        f,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
    )
    .unwrap();
    writeln!(
        f,
        "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">"
    )
    .unwrap();
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
    let info = format!("END={end_1based}{info_extras}");
    writeln!(
        f,
        "H1N1_HA\t{pos_1based}\t.\tA\t{alt}\t60\tPASS\t{info}\tGT\t{gt}"
    )
    .unwrap();
    path
}

/// Run gen-reads at the harness defaults (coverage=50, single-ended, mutation
/// rate disabled) with the optional `input_vcf`, returning the count of reads
/// emitted on H1N1_HA. The mutation rate is zeroed so de-novo SNPs don't
/// confound the read-count comparison — we're measuring the SV multiplier's
/// effect on coverage, not the variant generator's.
fn run_h1n1_reads(
    work: &std::path::Path,
    prefix: &str,
    input_vcf: Option<std::path::PathBuf>,
) -> usize {
    let mut config = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), prefix);
    config.coverage = 50;
    config.mutation_rate = Some(0.0);
    config.input_vcf = input_vcf;
    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();
    count_h1n1_reads(&work.join(format!("{prefix}_r1.fastq.gz")))
}

#[test]
fn whole_contig_hom_del_zeroes_all_reads() {
    // Aneuploidy-scale check #1: a homozygous <DEL> spanning the entire
    // H1N1_HA contig must produce zero reads. The depth multiplier for hom
    // <DEL> is 0.0 regardless of ploidy; the mut==0 override at runner.rs:486
    // then suppresses de-novo SNPs in the span too, so the read count must
    // be exactly zero.
    let (_dir, work) = fresh_workdir();
    let input_vcf = write_whole_contig_sv_vcf(&work, 1, 13325, "<DEL>", "", "1/1");
    let modulated = run_h1n1_reads(&work, "hom_del", Some(input_vcf));
    assert_eq!(
        modulated, 0,
        "homozygous <DEL> spanning whole H1N1_HA contig must produce zero reads; got {modulated}"
    );
}

#[test]
fn whole_contig_hom_dup_doubles_reads() {
    // Aneuploidy-scale check #2: a homozygous <DUP> spanning the entire
    // contig should produce ~2× the baseline read count. coverage_multiplier_for
    // returns (ploidy_f + ploidy_f) / ploidy_f = 2.0 at any ploidy.
    let (_dir, work) = fresh_workdir();
    let baseline = run_h1n1_reads(&work, "dup_baseline", None);
    assert!(
        baseline > 1000,
        "baseline H1N1 read count suspiciously low ({baseline})"
    );
    let input_vcf = write_whole_contig_sv_vcf(&work, 1, 13325, "<DUP>", "", "1/1");
    let modulated = run_h1n1_reads(&work, "dup_hom", Some(input_vcf));
    let ratio = modulated as f64 / baseline as f64;
    assert!(
        (1.7..=2.3).contains(&ratio),
        "hom <DUP> ratio {ratio:.3} outside [1.7, 2.3] — expected ~2.0× baseline. \
         baseline={baseline}, modulated={modulated}"
    );
}

#[test]
fn whole_contig_het_dup_yields_1_5x_reads() {
    // Aneuploidy-scale check #3: heterozygous <DUP> at ploidy 2 should yield
    // 1.5× baseline coverage_multiplier_for returns (ploidy + 1) / ploidy =
    // 1.5 for het DUP at ploidy 2.
    let (_dir, work) = fresh_workdir();
    let baseline = run_h1n1_reads(&work, "het_dup_baseline", None);
    let input_vcf = write_whole_contig_sv_vcf(&work, 1, 13325, "<DUP>", "", "0/1");
    let modulated = run_h1n1_reads(&work, "het_dup", Some(input_vcf));
    let ratio = modulated as f64 / baseline as f64;
    assert!(
        (1.3..=1.7).contains(&ratio),
        "het <DUP> ratio {ratio:.3} outside [1.3, 1.7] — expected ~1.5× baseline. \
         baseline={baseline}, modulated={modulated}"
    );
}

#[test]
fn whole_contig_cnv_cn6_triples_reads() {
    // Aneuploidy-scale check #4: a <CNV> with INFO/CN=6 at ploidy 2 should
    // yield 3× baseline. coverage_multiplier_for has an early return at
    // sv_model_defaults.rs:929 — when CN is present it overrides the
    // genotype-based default, returning CN / ploidy = 6 / 2 = 3.0.
    let (_dir, work) = fresh_workdir();
    let baseline = run_h1n1_reads(&work, "cnv_baseline", None);
    let input_vcf = write_whole_contig_sv_vcf(&work, 1, 13325, "<CNV>", ";CN=6", "0/1");
    let modulated = run_h1n1_reads(&work, "cnv_cn6", Some(input_vcf));
    let ratio = modulated as f64 / baseline as f64;
    assert!(
        (2.7..=3.3).contains(&ratio),
        "<CNV> CN=6 ratio {ratio:.3} outside [2.7, 3.3] — expected ~3.0× baseline. \
         baseline={baseline}, modulated={modulated}"
    );
}
