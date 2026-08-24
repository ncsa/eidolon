//! Shared machinery for **Gate 2** — realigning eidolon's FASTQ with a real aligner and asking
//! whether the result carries the evidence a structural-variant caller keys on.
//!
//! See `docs/sv_polish_roadmap.md` for what the gates are. The short version: every other SV
//! test inspects eidolon describing its own work (its FASTQ, its golden BAM, its truth VCF). A
//! caller sees none of those — it sees reads somebody else aligned. This module is the seam.
//!
//! Analysis is deliberately done in Rust over the SAM with `noodles` rather than by piping
//! `samtools` through `awk`: the arithmetic *is* the assertion, so it should be readable and
//! debuggable. That choice caught a bug on its first run — depth was being accumulated across
//! all eight H1N1 contigs, which backfilled a deleted window and turned a homozygous deletion
//! into an apparent 1.2x *enrichment*.

#![allow(dead_code)]

use super::{GenReadsConfig, eidolon, fresh_workdir, h1n1_reference};
use noodles::sam;
use noodles::sam::alignment::record::Sequence as _;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::io::{BufRead, Write as _};
use std::path::{Path, PathBuf};
use std::process::Command;

/// A symbolic SV to plant via `input_vcf`, with the windows its depth is measured over.
///
/// The windows are per-event on purpose. A 1.2 kb inversion cannot share a 300 bp deletion's
/// interior window, and it cannot share its flank either — `H1N1_HA:1000-1400` sits *inside* an
/// `H1N1_PB2:500-1700` inversion. Getting this wrong measures the wrong thing quietly, which is
/// how the first version of this harness turned a homozygous deletion into a 1.2x enrichment.
pub struct SvSpec {
    pub svtype: &'static str,
    pub contig: &'static str,
    pub pos: usize,
    pub end: usize,
    /// `"1/1"` for homozygous. Gates use hom so a signature is present or absent rather than
    /// halved — a het event turns every assertion below into a ratio judgement.
    pub gt: &'static str,
    /// Depth window strictly inside the event. Must clear each breakpoint by more than a
    /// fragment length where the event is large enough to allow it: coverage within ~one
    /// fragment of a junction is depressed by junction effects, not by the event's dosage
    /// (`docs/sv_support_matrix.md` measures 0.74–0.82 there for an inversion).
    pub interior: (usize, usize),
    /// Unaffected baseline on the same contig, so it shares that contig's coverage.
    pub flank: (usize, usize),
    /// Optional third window, for events whose expectation differs on either side of a
    /// breakpoint. A breakend needs it: downstream of the junction must collapse while
    /// immediately upstream must stay intact, and one ratio cannot express both.
    pub probe: Option<(usize, usize)>,
    /// For a breakend: the mate locus `(contig, pos)`. When set, `generate_reads` emits a
    /// **mated pair** of bracket-ALT records rather than one symbolic record, because a lone
    /// breakend is a different (and broken) thing — an unpaired `A.` destroys local coverage
    /// (#500). `svtype` is ignored for these; the geometry lives in the ALT.
    pub mate: Option<(&'static str, usize)>,
}

/// The reference base at `contig:pos` (1-based).
///
/// A breakend's REF must be the real anchor base: #451 was a truth VCF carrying a literal `N`
/// there, and a mismatched REF is a malformed record rather than a harmless approximation. The
/// H1N1 fixture has CRLF line endings, so `\r` is stripped — a detail that has already cost one
/// debugging session via `samtools faidx`.
pub fn ref_base(contig: &str, pos: usize) -> char {
    let text = std::fs::read_to_string(h1n1_reference()).unwrap();
    let mut seq = String::new();
    let mut in_contig = false;
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(name) = line.strip_prefix('>') {
            if in_contig {
                break;
            }
            in_contig = name.split_whitespace().next() == Some(contig);
            continue;
        }
        if in_contig {
            seq.push_str(line);
        }
    }
    assert!(
        !seq.is_empty(),
        "contig {contig} not found in the H1N1 fixture"
    );
    seq.chars()
        .nth(pos - 1)
        .unwrap_or_else(|| panic!("{contig}:{pos} is past the end of a {}bp contig", seq.len()))
        .to_ascii_uppercase()
}

/// Locate bwa-mem2, or fail with something actionable. **Never returns a "skip"** — a gate that
/// silently passes when its aligner is missing is worth less than no gate at all.
pub fn bwa_mem2() -> String {
    if let Ok(p) = std::env::var("BWA_MEM2") {
        assert!(
            Path::new(&p).is_file(),
            "BWA_MEM2={p} does not exist. Gate 2 cannot run without an aligner."
        );
        return p;
    }
    let found = Command::new("bwa-mem2")
        .arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    assert!(
        found,
        "bwa-mem2 not found on PATH.\n\
         Gate 2 realigns eidolon's FASTQ with a real aligner; there is nothing to assert \
         without one.\n\
         On this workstation it lives in the `aln` conda environment:\n\
         \n    conda activate aln\n\
         \nor point at it directly:\n\
         \n    BWA_MEM2=/path/to/bwa-mem2 cargo test --test <gate> -- --ignored\n"
    );
    "bwa-mem2".to_string()
}

/// Generate paired reads over H1N1, optionally planting one symbolic SV. The control run uses
/// the same reference, seed and coverage, so the two differ **only** by the variant.
pub fn generate_reads(
    work: &Path,
    tag: &str,
    sv: Option<&SvSpec>,
    novel: Option<&str>,
) -> (PathBuf, PathBuf) {
    let mut config = GenReadsConfig::new(h1n1_reference(), work.to_path_buf(), tag);
    config.coverage = 60;
    config.read_len = 100;
    config.paired_ended = true;
    config.produce_fastq = true;
    config.produce_bam = false;
    config.produce_vcf = true;
    // No de novo variants: the only difference between run and control must be the planted SV.
    config.mutation_rate = Some(0.0);
    config.sv_rate_scale = Some(0.0);

    if let Some(sv) = sv {
        let input_vcf = work.join(format!("{tag}.vcf"));
        let mut f = std::fs::File::create(&input_vcf).unwrap();
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
        match sv.mate {
            // A literal insertion: REF is the anchor base, ALT is that base plus the novel
            // sequence. Symbolic <INS> is deliberately NOT used — it carries no bases, so the
            // test would have no known answer to probe the reads for, which is the whole point
            // of an insertion gate.
            None if novel.is_some() => {
                let novel = novel.unwrap();
                let base = ref_base(sv.contig, sv.pos);
                writeln!(
                    f,
                    "{}\t{}\t.\t{base}\t{base}{novel}\t60\tPASS\tSVTYPE=INS;SVLEN={};END={}\tGT\t{}",
                    sv.contig,
                    sv.pos,
                    novel.len(),
                    sv.pos,
                    sv.gt
                )
                .unwrap();
            }
            None => {
                writeln!(
                    f,
                    "{}\t{}\t.\tG\t<{}>\t60\tPASS\tSVTYPE={};END={}\tGT\t{}",
                    sv.contig, sv.pos, sv.svtype, sv.svtype, sv.end, sv.gt
                )
                .unwrap();
            }
            Some((mc, mp)) => {
                // VCF 4.2 §5.4: a `t[p[` breakend at c1:p1 joining c2:p2 has as its mate the
                // record `]c1:p1]t'` at c2:p2. Case 1, a direct join with no reverse
                // complement, which is the geometry a caller reports as DEL/DUP-like rather
                // than converting to <INV> — chosen so this gate measures the breakend itself
                // and not the inversion-conversion path.
                let a = ref_base(sv.contig, sv.pos);
                let b = ref_base(mc, mp);
                writeln!(
                    f,
                    "{}\t{}\tbnd_a\t{a}\t{a}[{mc}:{mp}[\t60\tPASS\tSVTYPE=BND;MATEID=bnd_b\tGT\t{}",
                    sv.contig, sv.pos, sv.gt
                )
                .unwrap();
                writeln!(
                    f,
                    "{mc}\t{mp}\tbnd_b\t{b}\t]{}:{}]{b}\t60\tPASS\tSVTYPE=BND;MATEID=bnd_a\tGT\t{}",
                    sv.contig, sv.pos, sv.gt
                )
                .unwrap();
            }
        }
        config.input_vcf = Some(input_vcf);
    }

    let yaml = config.write_yaml();
    eidolon()
        .args(["gen-reads", "-c"])
        .arg(yaml.path())
        .assert()
        .success();

    let r1 = work.join(format!("{tag}_r1.fastq.gz"));
    let r2 = work.join(format!("{tag}_r2.fastq.gz"));
    assert!(r1.is_file() && r2.is_file(), "expected {r1:?} and {r2:?}");
    (r1, r2)
}

/// Align a FASTQ pair with bwa-mem2 and return the path to the SAM.
pub fn align(bwa: &str, work: &Path, tag: &str, r1: &Path, r2: &Path) -> PathBuf {
    // Index into the work dir so the repo's test_data is never written to.
    let local_ref = work.join("ref.fa");
    if !local_ref.exists() {
        std::fs::copy(h1n1_reference(), &local_ref).unwrap();
        let out = Command::new(bwa)
            .arg("index")
            .arg(&local_ref)
            .output()
            .expect("bwa-mem2 index failed to spawn");
        assert!(
            out.status.success(),
            "bwa-mem2 index failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }

    let sam_path = work.join(format!("{tag}.sam"));
    let out = Command::new(bwa)
        .args(["mem", "-t", "2"])
        .arg(&local_ref)
        .arg(r1)
        .arg(r2)
        .output()
        .expect("bwa-mem2 mem failed to spawn");
    assert!(
        out.status.success(),
        "bwa-mem2 mem failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::write(&sam_path, &out.stdout).unwrap();

    // A SAM with a header and no alignments would satisfy every "no signal" branch below, so
    // establish that reads were actually placed before anything is measured.
    let mapped = std::io::BufReader::new(std::fs::File::open(&sam_path).unwrap())
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.starts_with('@'))
        .count();
    assert!(
        mapped > 100,
        "{tag}: bwa-mem2 emitted only {mapped} alignment record(s) — the alignment step failed, \
         so nothing measured from it would mean anything"
    );
    sam_path
}

/// What a caller can see, extracted from a realigned SAM.
#[derive(Debug, Default)]
pub struct Signatures {
    /// Mean depth strictly inside the event (20 bp in from each breakpoint, so a read clipped
    /// at a junction cannot contribute and the window measures the event, not its edges).
    pub depth_inside: f64,
    /// Mean depth over `FLANK_WINDOW`, well clear of the event.
    pub depth_outside: f64,
    /// Soft clips whose clip point is within ±25 bp of either breakpoint.
    pub clips_at_breakpoints: usize,
    /// Soft clips anywhere else — the background the above must stand out from.
    pub clips_elsewhere: usize,
    /// Leftmost-read-of-pair with |TLEN| inflated well past the fragment mean. Necessary for a
    /// deletion but **not sufficient to identify one**: measured, a 1.2 kb inversion produces
    /// ~99 of these too, because a mate landing inside the inverted block maps to a mirrored
    /// position. Orientation is what discriminates — DEL shows long pairs with none of the
    /// others, DUP shows everted, INV shows same-orientation. Which is exactly why callers use
    /// insert size and orientation together rather than either alone.
    pub long_pairs: usize,
    /// Leftmost read reverse, mate forward — "everted" / RF orientation. The
    /// **tandem-duplication** signature: a pair spanning the duplication junction reads out of
    /// the second copy into the first, so the mates appear swapped.
    pub everted_pairs: usize,
    /// Both mates on the SAME strand (FF or RR). The **inversion** signature: one mate falls
    /// inside the inverted block and is therefore read from the opposite strand to its partner.
    /// Distinct from everted — an RF pair has mates on opposite strands, an inversion pair does
    /// not, and conflating them would let a DUP satisfy an INV assertion.
    pub same_orientation_pairs: usize,
    /// Read whose mate aligned to a DIFFERENT contig. The **translocation** signature for
    /// paired-end callers: an inter-chromosomal junction is the only thing that legitimately
    /// produces these.
    pub cross_contig_pairs: usize,
    /// Read carrying an `SA` tag whose supplementary alignment is on the mate contig — a split
    /// read spanning the junction. This is what a split-read caller assembles a breakend from,
    /// and it is strictly stronger evidence than a clip, which says only "something ends here".
    pub sa_to_mate_contig: usize,
    /// Mean depth over `SvSpec::probe`, if set.
    pub depth_probe: f64,
    pub reads: usize,
    /// Read that failed to align while its mate aligned fine. The **insertion** signature: a
    /// fragment landing wholly inside inserted sequence has no reference home, so the aligner
    /// can place its mate and not it. Callers (Manta, Delly) collect exactly these to assemble
    /// the inserted allele. Counted across all contigs, since an unmapped read has none.
    pub unmapped_with_mapped_mate: usize,
    /// Interior 30-mers of the novel sequence that appear in at least one read, either strand.
    pub novel_probe_hits: usize,
    /// How many interior probes were tried — the denominator `novel_probe_hits` is over. A
    /// hit count without it is a metric over an unknown denominator (CLAUDE.md rule 4).
    pub novel_probes: usize,
}

/// Accumulate the signatures over `contig` for an event spanning `pos..=end`.
pub fn analyse(sam_path: &Path, sv: &SvSpec, novel: Option<&str>) -> Signatures {
    let (contig, pos, end) = (sv.contig, sv.pos, sv.end);
    let mut reader = std::fs::File::open(sam_path)
        .map(std::io::BufReader::new)
        .map(sam::io::Reader::new)
        .unwrap();
    let _header = reader.read_header().unwrap();

    // Longest H1N1 contig is 2280 bp (PB2). Sized past it so no contig is silently truncated.
    let contig_len = 3_000usize;
    let mut depth = vec![0usize; contig_len];
    let mut sig = Signatures::default();

    // Only populated when there is a novel sequence to look for; the H1N1 fixture at 60x is
    // ~8 k reads, so holding them is cheap, and searching once per probe beats re-reading.
    let mut read_seqs: Vec<String> = Vec::new();

    for result in reader.records() {
        let record = result.unwrap();

        // Counted before the on-target filter below: an unmapped read has no reference name,
        // so filtering on one would discard every record this signature is about.
        if let Ok(flags) = record.flags() {
            if flags.is_unmapped() && !flags.is_mate_unmapped() {
                sig.unmapped_with_mapped_mate += 1;
            }
        }
        if novel.is_some() {
            read_seqs.push(
                record
                    .sequence()
                    .iter()
                    .map(|b| char::from(b).to_ascii_uppercase())
                    .collect(),
            );
        }

        // Breakend evidence lives on BOTH contigs, so these three counters are accumulated
        // before the on-target filter below (which exists for depth, a single-contig measure).
        if let Some((mate_contig, mate_pos)) = sv.mate {
            let this: Option<Vec<u8>> = record.reference_sequence_name().map(|n| n.to_vec());
            let mate: Option<Vec<u8>> = record.mate_reference_sequence_name().map(|n| n.to_vec());
            let ours = |n: &Option<Vec<u8>>| {
                n.as_deref()
                    .is_some_and(|n| n == contig.as_bytes() || n == mate_contig.as_bytes())
            };
            if ours(&this) && ours(&mate) && this != mate {
                sig.cross_contig_pairs += 1;
            }
            if ours(&this) {
                let sa = record
                    .data()
                    .get(&noodles::sam::alignment::record::data::field::Tag::OTHER_ALIGNMENTS);
                if let Some(Ok(value)) = sa {
                    let text = format!("{value:?}");
                    let other = if this.as_deref() == Some(contig.as_bytes()) {
                        mate_contig
                    } else {
                        contig
                    };
                    if text.contains(other) {
                        sig.sa_to_mate_contig += 1;
                    }
                }
                let _ = mate_pos;
            }
        }

        // H1N1 has EIGHT contigs. Without this filter every contig's reads land in one depth
        // array, the event window gets backfilled by unrelated contigs, and the signal vanishes.
        let on_target = matches!(
            record.reference_sequence_name(),
            Some(n) if n == contig.as_bytes()
        );
        if !on_target {
            continue;
        }
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        sig.reads += 1;
        let start = usize::from(start); // 1-based

        let mut ref_pos = start;
        let ops: Vec<_> = record.cigar().iter().map(|o| o.unwrap()).collect();
        for (i, op) in ops.iter().enumerate() {
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for p in ref_pos..(ref_pos + len).min(contig_len) {
                        depth[p] += 1;
                    }
                    ref_pos += len;
                }
                Kind::Deletion | Kind::Skip => ref_pos += len,
                Kind::SoftClip => {
                    // A leading clip points at the read's start; a trailing clip at ref_pos.
                    let clip_at = if i == 0 { start } else { ref_pos };
                    if clip_at.abs_diff(pos) <= 25 || clip_at.abs_diff(end) <= 25 {
                        sig.clips_at_breakpoints += 1;
                    } else {
                        sig.clips_elsewhere += 1;
                    }
                }
                _ => {}
            }
        }

        // Pair geometry, counted once per pair from the leftmost mate (TLEN > 0).
        let tlen = record.template_length().map(|t| t as i64).unwrap_or(0);
        if tlen > 0 {
            let event_len = end.saturating_sub(pos);
            if tlen as usize > 250 + event_len / 2 {
                sig.long_pairs += 1;
            }
            if let Ok(flags) = record.flags() {
                let rev = flags.is_reverse_complemented();
                let mate_rev = flags.is_mate_reverse_complemented();
                if rev && !mate_rev {
                    sig.everted_pairs += 1;
                }
                if rev == mate_rev {
                    sig.same_orientation_pairs += 1;
                }
            }
        }
    }

    let mean = |lo: usize, hi: usize| -> f64 {
        let slice = &depth[lo.min(contig_len)..hi.min(contig_len)];
        if slice.is_empty() {
            return 0.0;
        }
        slice.iter().sum::<usize>() as f64 / slice.len() as f64
    };
    sig.depth_inside = mean(sv.interior.0, sv.interior.1);
    sig.depth_outside = mean(sv.flank.0, sv.flank.1);
    if let Some((lo, hi)) = sv.probe {
        sig.depth_probe = mean(lo, hi);
    }

    // Interior probes only. The HEAD of an insertion is present even when the insertion is
    // only partially realized — that is exactly how #516 hid for eight releases — so probing
    // it would report success on a broken run. Five interior points, matching the Delta
    // harness's `verify_planted_ins`, so a local failure and a cluster failure mean the same
    // thing.
    if let Some(novel) = novel {
        let n = novel.len();
        if n >= 30 {
            let mut offsets: Vec<usize> = (1..=5)
                .map(|k| (n * k / 6).saturating_sub(15).min(n - 30))
                .collect();
            offsets.dedup();
            sig.novel_probes = offsets.len();
            for off in offsets {
                let probe = &novel[off..off + 30];
                let rc = super::revcomp(probe);
                if read_seqs
                    .iter()
                    .any(|r| r.contains(probe) || r.contains(&rc))
                {
                    sig.novel_probe_hits += 1;
                }
            }
        }
    }
    sig
}

/// Run a gate: generate with and without the SV, align both, and return `(with_sv, control)`.
pub fn run_gate(sv: &SvSpec) -> (Signatures, Signatures, tempfile::TempDir) {
    run_gate_inner(sv, None)
}

/// Run an **insertion** gate: the same run-versus-control shape, but the planted record is a
/// literal insertion of `novel` rather than a symbolic SV, so the reads can be probed for the
/// inserted bases themselves.
pub fn run_ins_gate(sv: &SvSpec, novel: &str) -> (Signatures, Signatures, tempfile::TempDir) {
    run_gate_inner(sv, Some(novel))
}

fn run_gate_inner(sv: &SvSpec, novel: Option<&str>) -> (Signatures, Signatures, tempfile::TempDir) {
    let bwa = bwa_mem2();
    let (dir, work) = fresh_workdir();

    let (vr1, vr2) = generate_reads(&work, "withsv", Some(sv), novel);
    let (cr1, cr2) = generate_reads(&work, "control", None, None);

    let with_sv = analyse(&align(&bwa, &work, "withsv", &vr1, &vr2), sv, novel);
    // The control is probed for the SAME novel sequence: it has no insertion, so every hit
    // there is a false positive and the probes are worthless. Without this the probe count
    // could be satisfied by sequence that was in the fixture all along.
    let control = analyse(&align(&bwa, &work, "control", &cr1, &cr2), sv, novel);
    (with_sv, control, dir)
}
