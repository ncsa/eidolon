#!/usr/bin/env python3
"""Compare per-site alt-allele fractions between two VCFs.

Purpose: validate SCN Phase 2 (issue #398 / realism epic #311). Feed the REAL
pool's observed per-site AFs (A, the input to gen-reads) and the SIMULATED golden
VCF gen-reads emitted (B). A high correlation + low per-site error means the
per-variant `allele_fraction` feature replayed the observed spectrum faithfully —
i.e. simulated AF ~= real AF, which is exactly Joao's question.

Sites are matched by (CHROM, POS, REF, ALT). For each side the AF is taken from,
in order: FORMAT/AF (the golden VCF's measured field), INFO/AF, or FORMAT/AD
(alt_depth / total_depth). Use AD-derived AF for pooled truth sets, where GT is
not diploid-meaningful.

Dependency-free (Python 3 stdlib only): no numpy/scipy/pysam. Pearson and
Spearman are computed directly.

Usage:
  python3 scripts/delta/scn_af_compare.py \
      --truth pool.af.vcf.gz --sim af_run.vcf.gz [--min-depth 20]
"""
import argparse
import collections
import gzip
import math
import sys
from collections import OrderedDict


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _field_index(fmt, key):
    parts = fmt.split(":")
    return parts.index(key) if key in parts else None


def _pick_per_allele(raw, alt_idx, n_alts, kind):
    """Select one allele's value from a per-allele VCF field.

    `kind` is "A" (one value per ALT, e.g. AF) or "R" (ref first, then one per ALT,
    e.g. AD). Returns None unless the field's arity actually matches the ALT count:
    a field of the wrong length cannot be indexed safely, and guessing would
    silently attribute one allele's number to another. The previous code took
    `[0]` unconditionally, which reported the FIRST allele's fraction for every
    allele of a multi-allelic record.
    """
    parts = raw.split(",")
    want = n_alts if kind == "A" else n_alts + 1
    if len(parts) != want:
        return None
    try:
        return float(parts[alt_idx if kind == "A" else alt_idx + 1])
    except ValueError:
        return None


def _af_for_allele(info, fmt, sample, alt_idx, n_alts):
    """Alt fraction for the ALT at 0-based `alt_idx`, and the site's total depth.

    Order: FORMAT/AF, INFO/AF, then FORMAT/AD. Returns (None, depth) when no
    fraction can be derived for this allele.
    """
    depth = _depth_from_ad(fmt, sample)
    if fmt and sample:
        af_i = _field_index(fmt, "AF")
        if af_i is not None:
            vals = sample.split(":")
            if af_i < len(vals):
                v = _pick_per_allele(vals[af_i], alt_idx, n_alts, "A")
                if v is not None:
                    return v, depth
    for kv in info.split(";"):
        if kv.startswith("AF="):
            v = _pick_per_allele(kv[3:], alt_idx, n_alts, "A")
            if v is not None:
                return v, depth
    # FORMAT/AD (Number=R). This is the path #450 depends on: reading the observed
    # side straight from `bcftools mpileup` means AD carries one count per allele,
    # so this allele's own count divided by the site total is its observed fraction.
    if fmt and sample:
        ad_i = _field_index(fmt, "AD")
        if ad_i is not None:
            vals = sample.split(":")
            if ad_i < len(vals):
                parts = vals[ad_i].split(",")
                if len(parts) == n_alts + 1:
                    try:
                        counts = [float(x) for x in parts]
                    except ValueError:
                        return None, depth
                    total = sum(counts)
                    if total > 0:
                        return counts[alt_idx + 1] / total, total
    return None, depth


def _depth_from_ad(fmt, sample):
    if not fmt or not sample:
        return None
    ad_i = _field_index(fmt, "AD")
    if ad_i is None:
        return None
    vals = sample.split(":")
    if ad_i >= len(vals):
        return None
    try:
        return sum(float(x) for x in vals[ad_i].split(","))
    except ValueError:
        return None


def load_af(path):
    """Return (alleles, sites) from a VCF.

    `alleles` maps (chrom, pos, ref, alt) -> (af, depth), with **multi-allelic
    records expanded one entry per literal ALT**. That expansion is what makes the
    mpileup-derived observed side joinable (#450): `bcftools mpileup` reports
    `ALT=T,<*>` with `AD=93,7,0`, so the truth's single ALT `T` has to be matched
    against the ALT *list* and its own AD element selected. Keying on the whole ALT
    string could never match, and skipping any record containing a symbolic
    alternative discarded the real base sitting beside `<*>`.

    `sites` maps (chrom, pos, ref) -> total depth, for every record with coverage.
    The caller uses it to score a truth ALT the observed side never listed as
    **0.0 observed** rather than as a missing site — see main().
    """
    alleles = OrderedDict()
    sites = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt_field = f[0], f[1], f[3], f[4]
            info = f[7] if len(f) > 7 else "."
            fmt = f[8] if len(f) > 8 else ""
            sample = f[9] if len(f) > 9 else ""
            alts = alt_field.split(",")
            depth = _depth_from_ad(fmt, sample)
            if depth is not None:
                sites[(chrom, pos, ref)] = depth
            for i, alt in enumerate(alts):
                # Skip symbolic / breakend ALTERNATIVES but keep literal ones from
                # the same record: mpileup's `<*>` non-ref placeholder sits beside a
                # real base, and dropping the whole record over it was the bug.
                if alt.startswith("<") or "[" in alt or "]" in alt or alt == ".":
                    continue
                af, dp = _af_for_allele(info, fmt, sample, i, len(alts))
                if af is not None:
                    alleles[(chrom, pos, ref, alt)] = (af, dp)
    return alleles, sites


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def _ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    if len(xs) < 2:
        return float("nan")
    return pearson(_ranks(xs), _ranks(ys))


def main():
    ap = argparse.ArgumentParser(description="Per-site AF correlation of two VCFs (truth vs sim)")
    ap.add_argument("--truth", required=True, help="VCF A — real/input per-site AFs")
    ap.add_argument("--sim", required=True, help="VCF B — simulated golden VCF")
    ap.add_argument("--min-depth", type=float, default=0.0,
                    help="skip sites whose total AD is below this on either side (default 0)")
    ap.add_argument("--max-uncovered-frac", type=float, default=0.10,
                    help="fail if more than this fraction of planted truth alleles go "
                         "unscored (default 0.10). A metric over an unknown denominator "
                         "is not a result, so this is enforced rather than warned about.")
    args = ap.parse_args()

    a, _a_sites = load_af(args.truth)
    b, b_sites = load_af(args.sim)

    # A truth ALT the observed side never listed, at a position it DID cover, means
    # zero reads carried that allele — an observed fraction of 0.0, which is a
    # measurement. Dropping those would recreate exactly the VAF-dependent exclusion
    # this fixes (#450), just at a lower threshold: the sites most likely to have no
    # alt reads at all are the lowest-VAF ones, so excluding them biases the result
    # optimistically. Fill them in and report how many.
    shared, zero_filled = [], 0
    for k in a:
        if k in b:
            shared.append(k)
            continue
        chrom, pos, ref, _alt = k
        dp = b_sites.get((chrom, pos, ref))
        if dp is not None and dp > 0:
            b[k] = (0.0, dp)
            shared.append(k)
            zero_filled += 1
    only_a = len(a) - len(shared)
    only_b = len(b) - len(shared)

    xs, ys, gated = [], [], 0
    for k in shared:
        (af_a, dp_a), (af_b, dp_b) = a[k], b[k]
        if args.min_depth > 0:
            if (dp_a is not None and dp_a < args.min_depth) or (
                dp_b is not None and dp_b < args.min_depth
            ):
                gated += 1
                continue
        xs.append(af_a)
        ys.append(af_b)

    print(f"truth sites={len(a)}  sim sites={len(b)}  shared={len(shared)}"
          f"  (only-truth={only_a}, only-sim={only_b})")
    if zero_filled:
        print(f"  {zero_filled} truth allele(s) had coverage but zero observed reads "
              "-> scored as observed AF 0.0 (not dropped)")
    if only_a:
        print(f"  {only_a} truth allele(s) had NO coverage on the observed side and are "
              f"excluded ({only_a / len(a):.1%} of the planted set).", file=sys.stderr)
    if args.min_depth > 0:
        print(f"min-depth {args.min_depth:g}: {gated} shared sites gated, {len(xs)} compared")
    if len(xs) < 2:
        sys.exit("fewer than 2 comparable sites — check inputs / --min-depth")

    diffs = [y - x for x, y in zip(xs, ys)]
    mae = sum(abs(d) for d in diffs) / len(diffs)
    rmse = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    print(f"n={len(xs)}  Pearson r={pearson(xs, ys):.4f}  Spearman rho={spearman(xs, ys):.4f}")
    print(f"MAE={mae:.4f}  RMSE={rmse:.4f}  mean(sim-truth)={sum(diffs)/len(diffs):+.4f}")
    print("  target: r>=0.95 and per-bin MAE within AF-estimation noise at that coverage")

    def _bin(af):
        return min(int(af * 10), 9)

    planted = collections.Counter(_bin(af) for af, _dp in a.values())

    print("per-AF-decile coverage and MAE (truth bin):")
    print(f"  {'bin':<12} {'planted':>8} {'scored':>8} {'unscored':>9}   MAE")
    shortfall = 0
    for i in range(10):
        if not planted.get(i):
            continue
        lo, hi = i / 10, i / 10 + 0.1
        bd = [abs(y - x) for x, y in zip(xs, ys) if _bin(x) == i]
        n_p, n_s = planted[i], len(bd)
        shortfall += n_p - n_s
        mae_s = f"{sum(bd) / len(bd):.4f}" if bd else "n/a — NOTHING SCORED"
        flag = "" if n_p == n_s else "  <-- incomplete"
        print(f"  [{lo:.1f},{hi:.1f})   {n_p:>8} {n_s:>8} {n_p - n_s:>9}   {mae_s}{flag}")
    print(f"  total planted={sum(planted.values())}  scored={len(xs)}  "
          f"unscored={shortfall} ({shortfall / sum(planted.values()):.1%})")

    # Enforced, not advised: a bias/MAE computed over a subset that silently omits the
    # lowest-VAF stratum is exactly the #450 failure, and it read as a clean PASS.
    uncovered = (len(a) - len(xs)) / len(a)
    if uncovered > args.max_uncovered_frac:
        sys.exit(f"FAIL: {uncovered:.1%} of planted truth alleles went unscored "
                 f"(limit {args.max_uncovered_frac:.1%}). The reported bias/MAE cover a "
                 f"subset of the planted set and must not be quoted. Attribute the "
                 f"shortfall per stratum above before believing any of it.")


if __name__ == "__main__":
    main()
