//! **Gate 2 for DUP** — does a real aligner produce the evidence a duplication caller needs?
//!
//! Motivated by a concrete failure. In SV validation array 21072620, Manta and Delly missed the
//! **identical** set of duplications (pooled TP=91, FN=14 for both). That is not caller
//! behaviour — two independent tools agreeing to the record points at a property of the events.
//! Investigation ruled out a size floor (one miss was 1.7 Mb) and assembly gaps (three of four
//! were 0.0% N), and one of the four turned out to be a `--passonly` artifact (#541). The rest
//! are unexplained, and the question Gate 2 answers is the one nobody had asked: *is the
//! evidence a duplication caller looks for actually present in the reads?*
//!
//! A tandem duplication has a signature a deletion does not:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | depth **rises** over the duplicated span | depth/CNV callers |
//! | soft clips at the junction | split-read callers |
//! | **everted (RF) pairs** — mates appear swapped | the discriminating tandem-DUP signal |
//!
//! The everted pair is the discriminating one: a fragment spanning the junction reads out of the
//! *second* copy back into the *first*, so its mates align in reverse order. A duplication that
//! raised depth without producing everted pairs would look like a copy-number gain with no
//! breakpoint — detectable by a depth caller and invisible to Manta or Delly, which is
//! consistent with what array 21072620 reported.
//!
//! ## Measured excess over the declared dosage — evidence for #499
//!
//! A homozygous DUP without `CN` gets `coverage_multiplier_for` = **2.0**, but the realigned
//! depth over this 300 bp event reads **2.55x** its flank. Mutating the multiplier gives a clean
//! dose-response:
//!
//! | multiplier | depth inside | flank | ratio |
//! |---|---|---|---|
//! | 2.0 (shipped) | 168.60 | 66.05 | 2.553 |
//! | 1.0 | 99.14 | 65.25 | 1.519 |
//! | 0.5 | 69.08 | 65.55 | 1.054 |
//!
//! That is `ratio ~= multiplier + 0.53` — a constant offset independent of the multiplier, and
//! ~44 everted pairs plus 64 clipped reads is about the right number of extra reads to account
//! for it. The junction reads appear to be emitted *in addition to* the coverage-multiplied
//! regular reads, rather than being part of them. Being a fixed number of reads, the excess
//! scales inversely with event length, which predicts the ~8% seen on #499's 1200 bp event
//! against the 27% seen here on 300 bp.
//!
//! **This gate deliberately does not assert the excess away.** It asserts the signature is
//! present; the dosage question belongs to #499 and to the Phase 1 size sweep, which this
//! measurement now gives a mechanism to test.
//!
//! Requires bwa-mem2; `#[ignore]`d so CI never reaches it. See `common/gate2.rs`.
//!
//! ```text
//! conda activate aln
//! cargo test --test gate2_realigned_dup -- --ignored --nocapture
//! ```

mod common;

use common::gate2::{SvSpec, run_gate};

/// Homozygous 300 bp tandem duplication. Hom so the depth change is a doubling rather than a
/// 1.5x, which keeps every assertion a presence/absence question.
const DUP: SvSpec = SvSpec {
    svtype: "DUP",
    contig: "H1N1_HA",
    pos: 500,
    end: 799,
    gt: "1/1",
};

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_dup_produces_the_three_signatures_a_caller_needs() {
    let (dup, ctl, _dir) = run_gate(&DUP);
    println!("  DUP run: {dup:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: depth rises over the duplicated span ─────────────────────────────
    // Homozygous DUP without CN means one extra copy per haplotype, so the multiplier is 2.0.
    // Asserted as a ratio against the same run's own flank, which cancels any coverage
    // difference between the two runs.
    let dup_ratio = dup.depth_inside / dup.depth_outside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    // Threshold is 1.9, just under the physically expected 2.0, NOT under the observed 2.55.
    // The first version used 1.5 and a mutation setting the multiplier to 1.0 landed at 1.519 --
    // surviving by 0.019. Anchor a threshold to what the value SHOULD be, never to a margin
    // below what it happens to be.
    assert!(
        dup_ratio > 1.9,
        "depth over a homozygous duplication should roughly double; interior/flank was \
         {dup_ratio:.3} (inside {:.1}x, outside {:.1}x)",
        dup.depth_inside,
        dup.depth_outside
    );
    // MUST-NOT-FIRE: the control has no duplication, so the same window must be flat. Without
    // this, a window that happens to be over-covered would satisfy the assertion above.
    assert!(
        (0.8..1.2).contains(&ctl_ratio),
        "control has no duplication, yet the same window reads {ctl_ratio:.3} — the window is \
         unreliable, so signature 1 proves nothing"
    );

    // ── Signature 2: split reads at the junction ──────────────────────────────────────
    assert!(
        dup.clips_at_breakpoints >= 5,
        "a split-read caller needs clipped reads at the duplication junction; found {}",
        dup.clips_at_breakpoints
    );
    assert!(
        dup.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "breakpoint clipping ({}) is not clearly above the control background ({})",
        dup.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 3: everted pairs — the tandem-duplication discriminator ─────────────
    // This is the one that separates "a duplication" from "a region with more coverage".
    assert!(
        dup.everted_pairs >= 5,
        "a tandem duplication must produce everted (RF) read pairs at its junction; found {}. \
         Without them the event is a copy-number gain with no breakpoint — visible to a depth \
         caller and invisible to Manta/Delly",
        dup.everted_pairs
    );
    assert!(
        dup.everted_pairs > ctl.everted_pairs * 3,
        "everted pairs ({}) are not clearly above the control background ({})",
        dup.everted_pairs,
        ctl.everted_pairs
    );
}
