# `compare-af` fixtures

`truth.vcf.gz` / `sim.vcf.gz` / `sim_450.vcf.gz` reproduce Delta job 20675479's
architecture: three VAF clusters of 104 / 183 / 280 sites at 0.05 / 0.15 / 0.45.
`sim.vcf.gz` drops 7 / 9 / 8 sites per cluster (the 4.2% shortfall actually observed);
`sim_450.vcf.gz` drops the ENTIRE lowest cluster, reproducing #450.

The `py_reference_*` files are the output of `scripts/delta/scn_af_compare.py` — the
Python implementation this subcommand replaced — captured at the moment of the port.
They exist so the migration's evidence survives deleting the Python: the Rust must
still produce these bytes exactly, on both the passing and the failing case.

That the #450 reproduction has BETTER headline numbers than the complete set
(bias +0.0011 vs +0.0008, MAE 0.0227 vs 0.0240) is the whole reason the coverage gate
is enforced rather than advisory.
