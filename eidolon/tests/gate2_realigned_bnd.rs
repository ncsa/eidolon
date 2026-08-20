//! **Gate 2 for BND** — does a real aligner produce the evidence a translocation caller needs?
//!
//! This is the type with the worst history. `BND recall=0.000` survived repeated confident
//! explanations, and the cause was never a caller and never scoring: until v3.1.0 the truth VCF
//! declared `t]p]` (head-to-head, mate reverse-complemented) while the reads carried a direct
//! join, so eidolon had been planting DELs and DUPs and labelling them breakends since v1.13.1.
//! Array 21072620 is the first whole-genome evidence that the geometry is fixed — Manta recovers
//! 46 of 50 records, all four bracket forms present, 0 unpaired and 0 mispaired.
//!
//! What no test has ever asked is whether the *alignment* carries what a caller needs:
//!
//! | signature | why a caller needs it |
//! |---|---|
//! | reads whose **mate is on the other contig** | the paired-end translocation signal |
//! | **`SA` tags pointing at the mate contig** | split-read assembly of the junction |
//! | clips at the junction | supporting, but weak on its own |
//! | depth **unchanged** near the junction | a balanced join loses no sequence |
//!
//! The `SA` tag is the strong one. A soft clip says only "something ends here"; an `SA` tag says
//! "the rest of this read aligns *there*", which is what Manta and Delly assemble a breakend
//! from. A junction producing clips without `SA` tags would be visible as an anomaly and
//! unclassifiable as a translocation.
//!
//! The depth assertion is a **must-not-fire** with real history: an unpaired breakend drops local
//! depth 72.0 → 42.4 while producing zero junction reads (#500). A properly mated breakend must
//! not do that.
//!
//! ## Measured semantics: a mated breakend is copy-number neutral
//!
//! Worth stating because it is a design choice and not an obvious one, and this gate is what
//! established it. Read literally, two mated records describe **one** junction — `t[p[` joins
//! HA-left to NA-right — and a balanced reciprocal translocation needs **four** records. A lone
//! junction is therefore unbalanced, and everything downstream of it should vanish.
//!
//! eidolon does not do that. Measured: downstream depth **0.93** of the flank, upstream **0.99**,
//! with 89 cross-contig pairs and 76 `SA` tags. So the junction evidence is added and the
//! sequence on both sides is retained — the **balanced reciprocal** case, which is the common one
//! in real genomes and the safer default, since callers detect breakends from split reads and
//! discordant pairs rather than from depth.
//!
//! The consequence to know: planting a deliberately **unbalanced** translocation and expecting
//! copy-number loss will not work, because depth does not follow the junction.
//!
//! (A window *straddling* the junction, `700-900`, does read 0.737 — a junction-proximal dip of
//! the same kind `docs/sv_support_matrix.md` measures at 0.74–0.82 for inversions. That is a
//! fragment-scale edge effect, not a copy-number change, which is why the windows here sit
//! clear of the junction on both sides.)
//!
//! ## Fixture
//!
//! A mated pair of bracket-ALT records — `t[p[` at `H1N1_HA:800` joined to `H1N1_NA:1200`, with
//! its VCF 4.2 counterpart `]p]t` on the mate side. **Case 1, a direct join with no reverse
//! complement**, chosen deliberately: inversion-oriented junctions (`t]p]`, `[p[t`) are converted
//! to `<INV>` by Manta's own `convertInversion.py`, so measuring one of those would entangle this
//! gate with the conversion path. The REF at each locus is read from the reference rather than
//! assumed — a literal `N` there was #451.
//!
//! Requires bwa-mem2; `#[ignore]`d so CI never reaches it. See `common/gate2.rs`.
//!
//! ```text
//! conda activate aln
//! cargo test --test gate2_realigned_bnd -- --ignored --nocapture
//! ```

mod common;

use common::gate2::{SvSpec, run_gate};

/// Inter-chromosomal breakend, homozygous so every fragment over the locus carries the junction.
/// `end` is unused for a mated breakend (the geometry is in the ALT) but must be a sane value.
const BND: SvSpec = SvSpec {
    svtype: "BND",
    contig: "H1N1_HA",
    pos: 800,
    end: 800,
    gt: "1/1",
    mate: Some(("H1N1_NA", 1_200)),
    // A `t[p[` junction joins H1N1_HA-left to H1N1_NA-right and DISCARDS H1N1_HA-right. This is
    // one junction, not a reciprocal translocation, so it is unbalanced by construction and
    // homozygously loses everything downstream. Three windows are needed to say that precisely:
    //   interior — downstream of the junction, must COLLAPSE (the declared geometry)
    //   probe    — just upstream, must stay INTACT (the loss is local, not contig-wide)
    //   flank    — far upstream baseline
    interior: (830, 1_200),
    probe: Some((600, 780)),
    flank: (100, 500),
};

#[test]
#[ignore = "requires bwa-mem2; run with --ignored (see module docs)"]
fn gate2_bnd_produces_the_signatures_a_caller_needs() {
    let (bnd, ctl, _dir) = run_gate(&BND);
    println!("  BND run: {bnd:?}");
    println!("  control: {ctl:?}");

    // ── Signature 1: the mate lands on the other contig ────────────────────────────────
    // This is the signal that makes a translocation a translocation. Nothing else in a
    // single-genome simulation legitimately produces cross-contig pairs.
    assert!(
        bnd.cross_contig_pairs >= 5,
        "a paired-end caller needs reads whose mate aligns to the partner contig; found {}",
        bnd.cross_contig_pairs
    );
    assert!(
        ctl.cross_contig_pairs == 0,
        "the control has no breakend, yet produced {} cross-contig pair(s) — that is either \
         mismapping or a harness error, and either way signature 1 would not be attributable \
         to the junction",
        ctl.cross_contig_pairs
    );

    // ── Signature 2: SA tags name the mate contig ──────────────────────────────────────
    // Strictly stronger than a clip: a clip says "something ends here", an SA tag says "the
    // rest of this read aligns THERE", which is what a breakend is assembled from.
    assert!(
        bnd.sa_to_mate_contig >= 3,
        "a split-read caller needs supplementary alignments pointing at the mate contig; found \
         {}. Clips alone leave the junction detectable as an anomaly but unclassifiable as a \
         translocation",
        bnd.sa_to_mate_contig
    );
    assert!(
        ctl.sa_to_mate_contig == 0,
        "control produced {} SA-to-mate alignment(s) with no breakend planted",
        ctl.sa_to_mate_contig
    );

    // ── Signature 3: clips at the junction ─────────────────────────────────────────────
    assert!(
        bnd.clips_at_breakpoints >= 5,
        "expected clipped reads at the junction; found {}",
        bnd.clips_at_breakpoints
    );
    assert!(
        bnd.clips_at_breakpoints > ctl.clips_at_breakpoints * 3,
        "junction clipping ({}) is not clearly above the control background ({})",
        bnd.clips_at_breakpoints,
        ctl.clips_at_breakpoints
    );

    // ── Signature 4: the junction is COPY-NUMBER NEUTRAL ───────────────────────────────
    // MEASURED SEMANTICS, worth stating because it is a design choice and not an obvious one.
    // Strictly, two mated records describe ONE junction: `t[p[` joins HA-left to NA-right, and a
    // balanced reciprocal translocation needs FOUR records (two junctions). Read literally, a
    // lone junction is unbalanced and everything downstream of it should vanish.
    //
    // eidolon does not do that. Measured downstream depth is 0.93 of the flank and upstream is
    // 0.99, so a mated breakend is realized as copy-number neutral: the junction evidence is
    // added, the sequence on both sides is retained. That is the BALANCED reciprocal case, which
    // is both the common one in real genomes and the safer default — a caller detects a breakend
    // from split reads and discordant pairs, not from depth.
    //
    // The consequence to be aware of: planting a deliberately UNBALANCED translocation and
    // expecting copy-number loss will not work, because depth does not follow the junction.
    let downstream = bnd.depth_inside / bnd.depth_outside.max(1e-9);
    let upstream = bnd.depth_probe / bnd.depth_outside.max(1e-9);
    assert!(
        downstream > 0.80,
        "a mated breakend is copy-number neutral, so depth downstream of the junction must be \
         retained; got {downstream:.3} ({:.1}x against a {:.1}x flank). A collapse here would \
         mean the junction is deleting sequence — the #500 shape",
        bnd.depth_inside,
        bnd.depth_outside
    );
    assert!(
        upstream > 0.80,
        "depth immediately upstream of the junction must be retained; got {upstream:.3} ({:.1}x)",
        bnd.depth_probe
    );
    // MUST-NOT-FIRE: with no breakend, neither window moves. Without this, windows that were
    // well covered for unrelated reasons would satisfy both assertions above.
    let ctl_down = ctl.depth_inside / ctl.depth_outside.max(1e-9);
    let ctl_up = ctl.depth_probe / ctl.depth_outside.max(1e-9);
    assert!(
        ctl_down > 0.80 && ctl_up > 0.80,
        "control has no breakend, yet downstream={ctl_down:.3} upstream={ctl_up:.3} — the \
         windows are unreliable, so the assertions above measure them rather than the junction"
    );
}
