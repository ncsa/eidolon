//! **Gate 2 for INS** — does a real aligner, given eidolon's FASTQ, produce the evidence a
//! structural-variant caller actually looks for at a long insertion?
//!
//! This is the gate that did not exist. DEL, DUP, INV and BND each had one; INS — the type
//! #516 was about — did not, and that is not a coincidence. #516 was a *read-level* defect:
//! eidolon emitted a truth VCF declaring insertions whose bases never reached the reads. Every
//! test that existed compared eidolon's output to eidolon's own truth, so both sides agreed and
//! the defect survived from v1.13.1 to v3.1.0, through three multi-hour Delta campaigns that
//! scored VCF against VCF and deleted the BAMs without ever asking whether the reads contained
//! what the truth declared.
//!
//! An insertion's signatures are not a deletion's. There is no depth collapse — the inserted
//! bases are absent from the reference, so nothing about reference coverage marks them. What a
//! caller keys on instead:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | the inserted bases are actually in the reads | without this there is nothing to assemble |
//! | soft clips pile up at the single breakpoint | split-read callers (Manta, Delly, LUMPY) |
//! | reads unmapped with a mapped mate | the insertion interior has no reference home |
//! | reference depth is NOT depleted | what separates an insertion from a deletion |
//!
//! The first is the #516 assertion and the reason this file exists. It probes **five interior
//! 30-mers**, never the head: a partially realized insertion has a perfectly good head, so a
//! head probe reports success on exactly the broken output this gate is meant to catch. The
//! same five-point rule is what `verify_planted_ins` applies on Delta, so a failure here and a
//! failure there mean the same thing.
//!
//! Every assertion is made against a **no-variant control** built from the same reference, seed
//! and coverage, so no signature can be satisfied by background. The control is probed for the
//! same novel sequence and must score zero — otherwise the probes are matching the fixture
//! rather than the insertion, and the headline assertion would prove nothing.
//!
//! ## Why this test is `#[ignore]`d
//!
//! CI has no aligner. Rather than skip silently when `bwa-mem2` is missing — a pass that means
//! nothing — the test **fails** with instructions if the binary is absent, and is `#[ignore]`d
//! so CI never reaches it. Run it deliberately:
//!
//! ```text
//! conda activate aln          # or: export BWA_MEM2=/path/to/bwa-mem2
//! cargo test --test gate2_realigned_ins -- --ignored --nocapture
//! ```

mod common;

use common::gate2::{SvSpec, run_ins_gate};
use common::synthetic_insert;

/// Homozygous 600 bp insertion. Hom so a signature is present or absent rather than halved.
///
/// 600 bp is chosen against the 250 bp fragment mean: comfortably longer, so fragments land
/// wholly inside the inserted sequence and have no reference home. At 250 bp or below the
/// unmapped-mate signature would be weak for a reason that has nothing to do with correctness.
///
/// `interior` and `flank` are both far from the anchor at 800. For an insertion they are not
/// measuring the event — they establish that reference depth is *undisturbed*, which is the
/// must-not-fire half of the gate. Coverage within a fragment length of the anchor is depressed
/// by junction effects rather than by the event, so a window straddling it would measure the
/// wrong thing.
const INS: SvSpec = SvSpec {
    svtype: "INS",
    contig: "H1N1_PB2",
    pos: 800,
    end: 800,
    gt: "1/1",
    interior: (400, 700),
    mate: None,
    probe: None,
    flank: (1_400, 1_900),
};

const INS_LEN: usize = 600;

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_ins_produces_the_evidence_a_caller_needs() {
    let novel = synthetic_insert(INS_LEN);
    let (ins, ctl, _dir) = run_ins_gate(&INS, &novel);
    println!("  INS run: {ins:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: the inserted bases reached the reads (this is #516) ──────────────
    assert!(
        ins.novel_probes >= 5,
        "the gate must try at least 5 interior probes; it tried {} — with no denominator the \
         hit count below is a metric over an unknown base",
        ins.novel_probes
    );
    assert_eq!(
        ins.novel_probe_hits, ins.novel_probes,
        "a HOMOZYGOUS {INS_LEN} bp insertion at 60x must put every interior probe in the reads; \
         only {}/{} were found. This is the #516 signature: the truth VCF declares an insertion \
         whose bases are not in the output.",
        ins.novel_probe_hits, ins.novel_probes
    );
    // MUST-NOT-FIRE: the control has no insertion, so the same probes must find nothing. If
    // they hit here the probes are matching the fixture and signature 1 proves nothing.
    assert_eq!(
        ctl.novel_probe_hits, 0,
        "the control has no insertion, yet {}/{} interior probes matched its reads — the novel \
         sequence is not novel with respect to the reference, so every hit above is suspect",
        ctl.novel_probe_hits, ctl.novel_probes
    );

    // ── Signature 2: split reads pile up at the breakpoint ────────────────────────────
    assert!(
        ins.clips_at_breakpoints >= 5,
        "a split-read caller needs clipped reads at the insertion point; found {}",
        ins.clips_at_breakpoints
    );
    assert!(
        ins.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "breakpoint clipping ({}) is not clearly above the control background ({}) — a caller \
         could not distinguish it either",
        ins.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 3: fragments inside the insertion have no reference home ────────────
    assert!(
        ins.unmapped_with_mapped_mate >= 5,
        "a {INS_LEN} bp insertion over a 250 bp fragment mean must strand whole fragments \
         inside itself, leaving reads unmapped with a mapped mate; found {}",
        ins.unmapped_with_mapped_mate
    );
    assert!(
        ins.unmapped_with_mapped_mate > ctl.unmapped_with_mapped_mate * 3,
        "unmapped-with-mapped-mate ({}) is not clearly above the control background ({})",
        ins.unmapped_with_mapped_mate,
        ctl.unmapped_with_mapped_mate
    );

    // ── Signature 4: reference depth is NOT depleted ──────────────────────────────────
    // The one that separates an insertion from a deletion. Without it a run that simply
    // dropped reads over the locus would satisfy signatures 2 and 3 and look like a pass.
    let ins_ratio = ins.depth_inside / ins.depth_outside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    assert!(
        ins_ratio > 0.70,
        "an insertion adds sequence and must not deplete reference coverage, but interior/flank \
         was {ins_ratio:.3} (inside {:.1}x, outside {:.1}x) — that is deletion-shaped",
        ins.depth_inside,
        ins.depth_outside
    );
    assert!(
        (ins_ratio - ctl_ratio).abs() < 0.30,
        "reference coverage away from the insertion should look like the control's: \
         {ins_ratio:.3} vs {ctl_ratio:.3}"
    );
}
