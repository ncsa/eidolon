//! **Gate 2 for INV** — does a real aligner produce the evidence an inversion caller needs?
//!
//! Two questions here, one of them load-bearing for how we report results.
//!
//! **1. Are the signatures present?** An inversion is *balanced* — no sequence is gained or lost —
//! so depth is not the signal. What a caller sees is:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | depth **unchanged** over the interior | distinguishes INV from DEL/DUP |
//! | soft clips at both breakpoints | split-read callers |
//! | **same-orientation (FF/RR) pairs** | the discriminating inversion signal |
//!
//! The same-orientation pair is the discriminator: one mate falls inside the inverted block and
//! is therefore read from the opposite strand to its partner, so both align to the same strand.
//! That is distinct from a duplication's everted (RF) pair, whose mates remain on *opposite*
//! strands — and the two are counted separately here so a DUP cannot satisfy an INV assertion.
//!
//! **2. Is the INV precision artifact really ours?** Array 21072620 pooled `manta_INV` at
//! precision 0.513 and `delly_INV` at 0.527 — both about half, with FP=56 and FP=52. We have
//! been attributing that to our own representation: inversion-oriented junctions convert to
//! `<INV>` via Manta's `convertInversion.py` and land as INV false positives. Two independent
//! callers agreeing that closely is good evidence, but it is still inference. If eidolon's
//! inversions produce clean FF/RR pairs and clean breakpoint clips, then the reads are not the
//! problem and the representational explanation survives a real test rather than a plausible one.
//!
//! ## Fixture choice matters here
//!
//! `docs/sv_support_matrix.md` records that coverage within ~one fragment of an inversion
//! junction sits at 0.74–0.82, and warns explicitly: *an inversion shorter than the fragment
//! length is not a clean fixture* — a 300 bp inversion with 200 bp fragments is entirely inside
//! its own junction-dip zone, which is what once made a correct implementation read as 0.63. So
//! this uses the 1.2 kb inversion on `H1N1_PB2` and measures 300 bp in from each breakpoint.
//!
//! Requires bwa-mem2; `#[ignore]`d so CI never reaches it. See `common/gate2.rs`.
//!
//! ```text
//! conda activate aln
//! cargo test --test gate2_realigned_inv -- --ignored --nocapture
//! ```

mod common;

use common::gate2::{SvSpec, run_gate};

/// Homozygous 1200 bp inversion — the fixture `docs/sv_support_matrix.md` used to establish that
/// the inverted sequence is covered at full depth and the loss is junction-proximal only.
const INV: SvSpec = SvSpec {
    svtype: "INV",
    contig: "H1N1_PB2",
    pos: 500,
    end: 1_700,
    gt: "1/1",
    // 300 bp in from each breakpoint: more than one fragment (mean 200), so the window measures
    // the inverted sequence rather than the junction dip.
    interior: (800, 1_400),
    mate: None,
    probe: None,
    flank: (1_750, 2_050),
};

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_inv_produces_the_signatures_a_caller_needs() {
    let (inv, ctl, _dir) = run_gate(&INV);
    println!("  INV run: {inv:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: depth is UNCHANGED — an inversion is balanced ─────────────────────
    // This is the assertion that separates INV from DEL and DUP. Note it is a must-NOT-fire in
    // the depth channel: a depth change here would mean the inversion is losing or duplicating
    // sequence, which is a different (and worse) defect than failing to be detectable.
    let inv_ratio = inv.depth_inside / inv.depth_outside.max(1e-9);
    let ctl_ratio = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    assert!(
        (0.85..1.15).contains(&inv_ratio),
        "an inversion is balanced, so interior depth must match the flank; got {inv_ratio:.3} \
         (inside {:.1}x, outside {:.1}x). Below 1.0 means sequence is being lost, above means \
         it is being duplicated",
        inv.depth_inside,
        inv.depth_outside
    );
    // The control establishes the windows themselves are comparable. Without it, an interior and
    // flank that happen to differ would make the assertion above either vacuous or unpassable.
    assert!(
        (0.85..1.15).contains(&ctl_ratio),
        "control interior/flank is {ctl_ratio:.3} — the two windows are not comparable, so \
         signature 1 measures the windows rather than the inversion"
    );

    // ── Signature 2: split reads at BOTH breakpoints ───────────────────────────────────
    assert!(
        inv.clips_at_breakpoints >= 5,
        "a split-read caller needs clipped reads at the inversion junctions; found {}",
        inv.clips_at_breakpoints
    );
    assert!(
        inv.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "breakpoint clipping ({}) is not clearly above the control background ({})",
        inv.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 3: same-orientation pairs — the inversion discriminator ──────────────
    assert!(
        inv.same_orientation_pairs >= 5,
        "an inversion must produce same-orientation (FF/RR) read pairs at its junctions; found \
         {}. Without them a caller sees clipped reads with no orientation evidence and cannot \
         classify the event as an inversion",
        inv.same_orientation_pairs
    );
    assert!(
        inv.same_orientation_pairs > ctl.same_orientation_pairs * 3,
        "same-orientation pairs ({}) are not clearly above the control background ({})",
        inv.same_orientation_pairs,
        ctl.same_orientation_pairs
    );

    // ── Must-not-fire across signatures: an inversion is not a duplication ─────────────
    // Everted (RF) pairs are the tandem-DUP signature. If an inversion produced them, a caller
    // would have grounds to call DUP here, and our INV false-positive story would be about the
    // reads rather than about representation.
    assert!(
        inv.everted_pairs <= ctl.everted_pairs + 2,
        "inversion produced {} everted (RF) pairs against the control's {} — that is a \
         DUPLICATION signature, and its presence would explain caller confusion by the reads \
         rather than by our representation",
        inv.everted_pairs,
        ctl.everted_pairs
    );
}
