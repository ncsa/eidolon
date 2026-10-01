#!/usr/bin/env python3
"""Assemble the data-derived cancer SvModel for #218 and splice it into the
bundled tumor model.

This is the normalization half of the refit, and the data-derived successor to
tools/inject_cancer_sv_model.py (which hand-injected literature values). It
takes:

  * the gen-mut-model FIT (tools/build_pcawg_sv_vcf.py -> gen-mut-model) for the
    distributions that ARE directly usable from the fit: per-type
    length_log_normal and the CNV copy-number distribution; and
  * the adapter SIDECAR (.counts.json) for per-donor event counts, which the
    fit cannot turn into a per-tumor rate/type-mix because the source corpora
    span different denominators (PCAWG SV ~2748 donors, CNA ~2778) and gnomAD
    INS is a germline site catalogue, not per-donor somatic counts.

Corrections applied (all documented, all overridable):

  type_probabilities: built from per-tumor event counts, NOT the fit's raw
    record ratios. DEL/DUP/INV/BND use PCAWG per-donor means. CNV and INS use
    literature per-tumor counts (--cnv-per-tumor / --ins-per-tumor) because the
    raw focal-CNA count (~135/donor) is segmentation, not discrete events, and
    INS rate was deliberately not taken from gnomAD's germline fraction.

  per_base_rate: (sum of per-tumor event counts) / haploid reference length.
    Divides the corpus aggregate down to a single-tumor rate (~3.5e-8), which is
    what the cosmic_bundle.rs <1e-6 gate expects.

  cnv_copy_number_distribution: from the fit, but with the diploid-neutral CN
    (--exclude-cn, default 2) dropped and the high tail capped (--cnv-cap),
    then renormalized. (The raw CNA tail runs to CN>500 — ecDNA/amplicon
    territory, #192, not the focal-CNV model here.)

  homozygous_frequency: preserved from --base-model. BEDPE/CNA carry no
    per-event zygosity, so it is not data-derivable here.

  length_log_normal[Bnd] = [0,0] (BNDs are point events; no length).

Usage:
  tools/normalize_pcawg_sv_model.py \
      --fitted-model  ~/code/data/pcawg/pcawg_sv_fit.json.gz \
      --sidecar       ~/code/data/pcawg/pcawg_sv_corpus.counts.json \
      --base-model    tools/cosmic_v104_pancancer_model.json.gz \
      --out           tools/cosmic_v104_pancancer_model.json.gz
"""
import argparse
import gzip
import json
import sys

# SvType JSON keys (TitleCase, matching serde output).
SV_TYPES = ["Del", "Dup", "Inv", "Bnd", "Cnv", "Ins"]
SIDECAR_KEY = {"Del": "DEL", "Dup": "DUP", "Inv": "INV",
               "Bnd": "BND", "Cnv": "CNV", "Ins": "INS"}


def load_gz_json(path):
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def write_gz_json(path, obj):
    with gzip.open(path, "wt") as fh:
        json.dump(obj, fh)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fitted-model", required=True,
                    help="gen-mut-model output JSON.gz (for length/CN dists)")
    ap.add_argument("--sidecar", required=True,
                    help="adapter .counts.json (per-donor counts)")
    ap.add_argument("--base-model", required=True,
                    help="model to splice into (for non-SV fields + homozygous_frequency)")
    ap.add_argument("--out", required=True, help="output model JSON.gz")
    ap.add_argument("--haploid-reflen", type=float, default=3.1e9,
                    help="haploid reference length for per_base_rate denominator "
                         "(default 3.1e9 = GRCh38 primary chromosomes)")
    ap.add_argument("--cnv-per-tumor", type=float, default=6.0,
                    help="literature focal-CNV events per tumor (default 6). "
                         "Raw focal-CNA count (~135/donor) is segmentation, not "
                         "discrete events; see #218.")
    ap.add_argument("--ins-per-tumor", type=float, default=3.0,
                    help="literature somatic INS (MEI) events per tumor "
                         "(default 3). gnomAD supplies INS length, not rate.")
    ap.add_argument("--cnv-cap", type=int, default=10,
                    help="collapse CN values above this into the cap when "
                         "building cnv_copy_number_distribution (default 10)")
    ap.add_argument("--exclude-cn", type=int, nargs="*", default=[2],
                    help="CN values to drop as diploid-neutral (default: 2)")
    args = ap.parse_args()

    fit = load_gz_json(args.fitted_model)
    base = load_gz_json(args.base_model)
    with open(args.sidecar) as fh:
        side = json.load(fh)

    fit_sv = fit.get("sv_model")
    if not fit_sv:
        sys.exit("fitted model has no sv_model — did gen-mut-model fit SVs?")
    base_sv = base.get("sv_model")
    if not base_sv:
        # INS length, the CNV fallback and homozygous_frequency come from the base's
        # sv_model. Without one they would default silently: INS length to a degenerate
        # [0, 0], homozygous_frequency to 0.1 (#761). A model from
        # `fetch_cosmic_corpus.sh --train` has no sv_model; graft the shipped pan-cancer
        # one onto it first (graft_sv_model.py).
        sys.exit(f"base model {args.base_model} has no sv_model, so INS length and "
                 "homozygous_frequency have no source. Use a base that carries one, such "
                 "as tools/cosmic_v104_pancancer_model.json.gz.")

    means = side["notes"]["per_donor_means"]  # DEL/DUP/INV/BND (+ CNV raw)

    # ── type_probabilities + per_base_rate from per-tumor event counts ──────
    per_tumor = {
        "Del": means["DEL"], "Dup": means["DUP"],
        "Inv": means["INV"], "Bnd": means["BND"],
        "Cnv": args.cnv_per_tumor, "Ins": args.ins_per_tumor,
    }
    total_per_tumor = sum(per_tumor.values())
    type_probabilities = {t: per_tumor[t] / total_per_tumor for t in SV_TYPES}
    per_base_rate = total_per_tumor / args.haploid_reflen

    # ── length_log_normal: from the fit; BND fixed at [0,0] ────────────────
    length_log_normal = {}
    fit_len = fit_sv.get("length_log_normal", {})
    for t in SV_TYPES:
        if t == "Bnd":
            length_log_normal[t] = [0.0, 0.0]
        elif t in fit_len:
            length_log_normal[t] = fit_len[t]
        else:
            # Fall back to the base model if the fit dropped a thin type.
            length_log_normal[t] = base_sv.get("length_log_normal", {}).get(
                t, [0.0, 0.0])

    # ── cnv_copy_number_distribution: drop neutral, cap tail, renormalize ──
    raw_cn = side.get("cnv_copy_number_hist", {})
    capped = {}
    for cn_str, n in raw_cn.items():
        cn = int(cn_str)
        if cn in args.exclude_cn:
            continue
        cn = min(cn, args.cnv_cap)
        capped[cn] = capped.get(cn, 0) + n
    cn_total = sum(capped.values())
    if cn_total > 0:
        cnv_cn_dist = {str(cn): capped[cn] / cn_total
                       for cn in sorted(capped)}
    else:
        cnv_cn_dist = base_sv.get("cnv_copy_number_distribution", {})

    homozygous_frequency = base_sv.get("homozygous_frequency", 0.1)

    new_sv = {
        "per_base_rate": per_base_rate,
        "type_probabilities": type_probabilities,
        "length_log_normal": length_log_normal,
        "cnv_copy_number_distribution": cnv_cn_dist,
        "homozygous_frequency": homozygous_frequency,
    }

    # ── report old vs new + gate checks ────────────────────────────────────
    old_sv = base_sv

    def fmt_probs(d):
        return " ".join(f"{t}={d.get(t, 0):.3f}" for t in SV_TYPES)

    print("── type_probabilities ──", file=sys.stderr)
    print(f"  old: {fmt_probs(old_sv.get('type_probabilities', {}))}", file=sys.stderr)
    print(f"  new: {fmt_probs(type_probabilities)}", file=sys.stderr)
    print(f"── per_base_rate ──  old={old_sv.get('per_base_rate'):.3e}  "
          f"new={per_base_rate:.3e}", file=sys.stderr)
    print("── length_log_normal (mu, sigma) ──", file=sys.stderr)
    for t in SV_TYPES:
        o = old_sv.get("length_log_normal", {}).get(t)
        print(f"  {t}: old={o}  new={[round(x,3) for x in length_log_normal[t]]}",
              file=sys.stderr)
    print(f"── cnv CN dist (new) ── {cnv_cn_dist}", file=sys.stderr)
    print(f"── homozygous_frequency ── {homozygous_frequency} (preserved)",
          file=sys.stderr)

    prob_sum = sum(type_probabilities.values())
    bnd = type_probabilities["Bnd"]
    cn_sum = sum(cnv_cn_dist.values())
    print("\n── GATE CHECKS (cosmic_bundle.rs) ──", file=sys.stderr)
    print(f"  all 6 types present: {all(t in type_probabilities for t in SV_TYPES)}",
          file=sys.stderr)
    print(f"  type_probabilities sum = {prob_sum:.6f}  (need ~1.0): "
          f"{'OK' if abs(prob_sum-1) < 1e-6 else 'FAIL'}", file=sys.stderr)
    print(f"  cnv CN dist sum = {cn_sum:.6f}  (need ~1.0): "
          f"{'OK' if abs(cn_sum-1) < 1e-6 else 'FAIL'}", file=sys.stderr)
    print(f"  per_base_rate < 1e-6: {'OK' if per_base_rate < 1e-6 else 'FAIL'} "
          f"({per_base_rate:.3e})", file=sys.stderr)
    print(f"  BND fraction >= 0.20: {'OK' if bnd >= 0.20 else 'FAIL'} "
          f"({bnd:.3f})", file=sys.stderr)

    base["sv_model"] = new_sv
    write_gz_json(args.out, base)
    print(f"\n>> wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
