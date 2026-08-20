//! **Gate 2 for DEL** — does a real aligner, given eidolon's FASTQ, produce the evidence a
//! structural-variant caller actually looks for?
//!
//! Every other SV test in this repo inspects eidolon's own output: the FASTQ bases, the golden
//! BAM, the truth VCF. All three are eidolon describing its own work. A caller never sees any of
//! them — it sees reads that someone else aligned. **Nothing has ever checked that step.** The
//! whole-genome campaigns jump straight from FASTQ to "did Manta find it", which conflates
//! "the evidence is correct" with "the evidence was sufficient", and cannot tell you which
//! failed when the answer is no.
//!
//! So this aligns eidolon's reads with **bwa-mem2** and asserts, on the resulting alignment, the
//! three signatures every read-based DEL caller keys on:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | depth collapses over the deleted interval | depth-based callers (CNVnator, GATK gCNV) |
//! | soft clips pile up at both breakpoints | split-read callers (Manta, Delly, LUMPY) |
//! | spanning pairs have TLEN inflated by ~SVLEN | discordant-pair callers (all of the above) |
//!
//! Each is asserted against a **no-variant control** built from the same reference, seed and
//! coverage, so "the signature is present" cannot be satisfied by background noise — which is
//! how a 13% artifact once read as a real DUP defect (see `docs/sv_support_matrix.md`).
//!
//! ## Why this test is `#[ignore]`d
//!
//! CI has no aligner. Rather than skip silently when `bwa-mem2` is missing — a pass that means
//! nothing, the exact shape this project keeps re-earning — the test **fails** with instructions
//! if the binary is absent, and is `#[ignore]`d so CI never reaches it. Run it deliberately:
//!
//! ```text
//! conda activate aln          # or: export BWA_MEM2=/path/to/bwa-mem2
//! cargo test --test gate2_realigned_del -- --ignored --nocapture
//! ```
//!
//! The analysis is done in Rust over the SAM with `noodles`, not by piping samtools through awk,
//! so it is debuggable and its arithmetic is visible.

mod common;

use common::gate2::{SvSpec, run_gate};

/// Homozygous 299 bp deletion. Hom so every fragment over the locus carries it — a het event
/// halves every signal below and turns each assertion into a ratio judgement.
const DEL: SvSpec = SvSpec {
    svtype: "DEL",
    contig: "H1N1_HA",
    pos: 500,
    end: 799,
    gt: "1/1",
    // 20 bp in from each breakpoint: enough that a read clipped at a junction cannot contribute,
    // and this event is too short (300 bp vs 200 bp fragments) for a fragment-length margin.
    interior: (520, 779),
    mate: None,
    probe: None,
    flank: (1_000, 1_400),
};

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_del_produces_the_three_signatures_a_caller_needs() {
    let (del, ctl, _dir) = run_gate(&DEL);
    println!("  DEL run: {del:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: depth collapses over the deletion ────────────────────────────────
    let del_ratio = del.depth_inside / del.depth_outside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    assert!(
        del_ratio < 0.10,
        "depth over a HOMOZYGOUS deletion should collapse; interior/flank was {del_ratio:.3} \
         (inside {:.1}x, outside {:.1}x)",
        del.depth_inside,
        del.depth_outside
    );
    // MUST-NOT-FIRE: the control has no deletion, so the same window must be normally covered.
    // Without this, a run producing no reads there for an unrelated reason would satisfy the
    // assertion above — which is exactly how the first version of this harness fooled itself.
    assert!(
        ctl_ratio > 0.80,
        "control has no deletion, yet the same window is depleted (ratio {ctl_ratio:.3}) — \
         the window itself is unreliable, so signature 1 proves nothing"
    );

    // ── Signature 2: split reads pile up at the breakpoints ───────────────────────────
    assert!(
        del.clips_at_breakpoints >= 5,
        "a split-read caller needs clipped reads at the junction; found {}",
        del.clips_at_breakpoints
    );
    assert!(
        del.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "breakpoint clipping ({}) is not clearly above the control background ({}) — a caller \
         could not distinguish it either",
        del.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 3: spanning pairs are discordant by ~SVLEN ──────────────────────────
    assert!(
        del.long_pairs >= 5,
        "a paired-end caller needs pairs whose insert size is inflated by the deletion; found {}",
        del.long_pairs
    );
    assert!(
        del.long_pairs > ctl.long_pairs * 3,
        "long pairs ({}) are not clearly above the control background ({})",
        del.long_pairs,
        ctl.long_pairs
    );
}
