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
    // 20 bp in from each breakpoint: enough that a read clipped at a junction cannot contribute,
    // and this event is too short (300 bp vs 200 bp fragments) for a fragment-length margin.
    interior: (520, 779),
    mate: None,
    probe: None,
    flank: (1_000, 1_400),
};

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_dup_produces_the_three_signatures_a_caller_needs() {
    let (dup, ctl, _dir) = run_gate(&DUP);
    println!("  DUP run: {dup:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: depth rises over the duplicated span ─────────────────────────────
    // Homozygous DUP without CN means one extra copy per haplotype, so the multiplier is 2.0.
    //
    // The denominator is THE SAME WINDOW IN A SEPARATE NO-VARIANT CONTROL RUN, not this run's
    // own flank. #499's closing note prescribes exactly that, and #582 is what ignoring it
    // costs: dividing by the in-run flank, this gate read 2.553 before the fragment-placement
    // rewrite and 1.739 after, and failed. The duplication had not regressed — measured
    // against the control run the same two builds give 2.671 and 1.970, so the dosage
    // actually moved from 34% over-delivered to within 1.5% of physical expectation. What
    // changed was the denominator: placement went from an artificially uniform tiler (depth
    // VMR 0.226) to genuine Poisson (~1.0), and a 400 bp flank window is small enough to
    // scatter. An in-run flank is not a control; it is another sample of the same run.
    //
    // The `ratio ~= multiplier + 0.53` offset described above is NOT simply gone. Measured
    // at multiplier 2.0 the ratio is 1.970, so no excess is visible there — but a mutation
    // run forcing the homozygous multiplier to 1.0 lands at 1.520, i.e. +0.52, exactly the
    // documented offset. So the junction reads are still emitted in addition to the
    // multiplied ones; whatever cancels them at 2.0 is not understood, and #499 should not
    // be considered closed on the strength of the 1.970 alone.
    let dup_ratio = dup.depth_inside / ctl.depth_inside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    println!(
        "  interior/control-interior = {dup_ratio:.3}  (in-run flank ratio would be {:.3})",
        dup.depth_inside / dup.depth_outside.max(1e-9)
    );
    // A BAND around the physically expected 2.0, not a one-sided floor. Anchored to what the
    // value should be, never to a margin below what it happens to be — and two-sided so an
    // over-delivering multiplier fails too, which a `> 1.9` floor would have waved through at
    // the 2.671 this gate used to see.
    assert!(
        (1.8..2.2).contains(&dup_ratio),
        "depth over a homozygous duplication should double; interior vs the control run's \
         same window was {dup_ratio:.3} (dup {:.1}x, control {:.1}x), outside 1.8-2.2",
        dup.depth_inside,
        ctl.depth_inside
    );
    // MUST-NOT-FIRE: the control has no duplication, so its own interior and flank must
    // agree. If they do not, the window is unreliable and the ratio above proves nothing.
    assert!(
        (0.8..1.2).contains(&ctl_ratio),
        "control has no duplication, yet its interior/flank reads {ctl_ratio:.3} — the window \
         is unreliable, so signature 1 proves nothing"
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
