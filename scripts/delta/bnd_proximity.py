#!/usr/bin/env python3
"""Score breakend (BND) recovery by breakpoint proximity.

`truvari bench` matches sequence-resolved SVs by position *and size overlap*. A
breakend has no SVLEN and `END == POS`, so truvari cannot match one at all and
reports recall 0.000 for BND no matter how correct the truth is (#457). This scores
BNDs the way the CNV direction check in sv_pipeline.sbatch already scores CNVs:
against the thing the caller can actually be judged on.

A truth *junction* is recovered if the caller reports a junction whose two
breakpoints both land within `--window` bp of the truth's two breakpoints. Two
properties that matter:

* **Junctions, not records.** A junction is one event described by two mate records,
  so counting records would weight BND 2x against one-record types like DEL. Both
  sides are canonicalized to a single junction key before counting.

* **Representation-agnostic.** Callers legitimately re-represent a junction as
  something other than a BND pair — Manta reports tandem-duplication junctions as
  `<DUP:TANDEM>` at the correct coordinates. Scoring only caller BND records would
  charge that as a miss when the caller in fact found the event. So a caller
  junction is derived from *either* a BND ALT's mate coordinate *or* a symbolic
  record's (POS, END) span.

Output is deliberately labelled: recall over truth junctions is the meaningful
number. "Unmatched caller junctions" is reported for context but is NOT precision —
a somatic SV caller emits junctions for events outside the truth set by design.
"""

import argparse
import gzip
import re
import sys

# t[p[ / t]p] / [p[t / ]p]t — capture the mate locus regardless of orientation.
BND_ALT = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _info(field):
    out = {}
    if field in (".", ""):
        return out
    for kv in field.split(";"):
        if "=" in kv:
            k, _, v = kv.partition("=")
            out[k] = v
        else:
            out[kv] = True
    return out


def junction_key(c1, p1, c2, p2):
    """Canonical, orientation-independent key for one junction."""
    a, b = (c1, p1), (c2, p2)
    return (a, b) if a <= b else (b, a)


def read_junctions(path, bnd_only):
    """Collect junctions from a VCF.

    `bnd_only=True` for the truth (we know it is a BND set). For caller VCFs we also
    accept symbolic records with an END, so a re-represented junction still counts.
    Returns a dict of canonical key -> a representative description string.
    """
    out = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, alt, info = f[0], int(f[1]), f[4], _info(f[7])
            m = BND_ALT.search(alt)
            if m:
                mate_c, mate_p = m.group(1), int(m.group(2))
                out.setdefault(junction_key(chrom, pos, mate_c, mate_p), line.rstrip("\n"))
                continue
            if bnd_only:
                continue
            # Symbolic record with a span: its two edges are a junction.
            end = info.get("END")
            if alt.startswith("<") and end is not None:
                try:
                    end = int(end)
                except ValueError:
                    continue
                if end > pos:
                    out.setdefault(junction_key(chrom, pos, chrom, end), line.rstrip("\n"))
    return out


def matches(truth_key, call_key, window):
    """Both breakpoints within `window`, in the same pairing."""
    (tc1, tp1), (tc2, tp2) = truth_key
    (cc1, cp1), (cc2, cp2) = call_key
    return (
        tc1 == cc1
        and tc2 == cc2
        and abs(tp1 - cp1) <= window
        and abs(tp2 - cp2) <= window
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, help="truth BND VCF (SVTYPE=BND subset)")
    ap.add_argument("--calls", required=True, help="caller VCF (all types)")
    ap.add_argument("--window", type=int, default=200,
                    help="max bp offset per breakpoint (default 200)")
    ap.add_argument("--label", default="bnd", help="label for the summary line")
    args = ap.parse_args()

    truth = read_junctions(args.truth, bnd_only=True)
    calls = read_junctions(args.calls, bnd_only=False)

    matched, unmatched_truth = [], []
    used = set()
    for tk in truth:
        hit = next((ck for ck in calls if ck not in used and matches(tk, ck, args.window)), None)
        if hit is None:
            unmatched_truth.append(tk)
        else:
            used.add(hit)
            matched.append((tk, hit))

    n_truth, n_match = len(truth), len(matched)
    recall = (n_match / n_truth) if n_truth else 0.0

    print(f"truth junctions={n_truth}  matched={n_match}  "
          f"missed={len(unmatched_truth)}  window={args.window}bp")
    print(f"  {args.label}_BND_proximity recall={recall:.3f}  "
          f"(junctions, not breakend records)")
    print(f"  caller junctions considered={len(calls)}  "
          f"unmatched={len(calls) - len(used)}  "
          "(NOT precision — a somatic caller reports junctions outside the truth "
          "set by design)")
    if unmatched_truth:
        print("  missed truth junctions:")
        for (c1, p1), (c2, p2) in unmatched_truth[:10]:
            print(f"    {c1}:{p1} <-> {c2}:{p2}")
        if len(unmatched_truth) > 10:
            print(f"    ... and {len(unmatched_truth) - 10} more")

    # A truth set with no junctions means the caller was never actually tested;
    # say so rather than printing recall=0.000 as though it were a measurement.
    if n_truth == 0:
        print("  NOTE: no truth junctions — this is not a measurement of the caller.",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
