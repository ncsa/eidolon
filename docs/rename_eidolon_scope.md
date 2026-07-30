# Rename scope — `rusty-neat` / `rneat` → `eidolon`

Planning doc for the project rename. **Status: Phase A DONE (PR #425, merged to `develop`;
ships in the v2.0.0 release PR #426). Phases B–C pending. Decision-1 output-token migration
DONE (v3.0.0) — see Decisions §1.**
Target release for the rename: **v2.0.0**, cut *after* v1.21.1 is fully released. Do not
entangle the rename with any behavior change.

## Why `eidolon`

Greek εἴδωλον — "a phantom; an exact, insubstantial double of a real thing" (the *eidolon
of Helen* was a perfect phantom sent to Troy in her place). That is precisely what the tool
produces: a realistic phantom double of real sequencing data. It survived the deepest
vetting of the candidates:

- **Namespaces:** bioconda **free** (our distribution channel); crates.io taken → the published
  crate is named `eidolon-core` (see below), tool/binary/conda package stay `eidolon`.
- **No in-domain clash.** Every real-world use of "Eidolon" is in gaming/entertainment
  (Steam titles, a games studio, tabletop RPGs) — a different domain that will neither confuse
  a genomics audience nor compete for bioinformatics search results.
- Runner-up **`verism`** (free in *every* registry; only weakness is phonetic proximity to the
  simulator *VarSim*) is retained as the fallback if `eidolon` ever hits a snag.
- `galatea` was rejected: it collides in-domain with **Galatea Bio**, a funded genomics/NGS
  company, plus the Inkitt "Galatea" reading app.

The NEAT pedigree stays documented in the README and `CITATION.cff` — the rename changes the
project's identity, not its heritage.

## The one rule that governs the rename

**No blanket `s/neat/eidolon/`.** Three different kinds of "neat" live in the tree and must be
treated differently:

| token class | examples | action |
|---|---|---|
| **Identity** | `rneat`, `rusty-neat`, `RNEAT_BIN`, `RNEAT_REPO`, `rusty-neat-tests.yml` | → rename to `eidolon` / `EIDOLON_` |
| **Pedigree** (keep) | `NEAT`, `NEAT2`, `NEAT4`, `neat-genreads`, comparison text | → **leave untouched** |
| **Output-format tokens** | `NEAT_PROVENANCE`, `NEAT_simulated_sample`, `NEAT_ORIGIN`, `RNEAT_chimeric_*` | → kept for v2.0.0; **migrated to `EIDOLON_*` in v3.0.0** (decision 1) |

Codebase inventory at scoping time: `rneat` 866× / 129 files, `rusty-neat` 93× / 32 files,
`RNEAT` 204× / 46 files (mostly `RNEAT_BIN`/`RNEAT_REPO` env vars + `RNEAT_chimeric_` prefixes).

## Decisions (locked)

1. **Output-format tokens — keep in v2.0.0, migrate later. → DONE (v3.0.0).** The `NEAT_*` INFO tags
   (`NEAT_PROVENANCE`, `NEAT_ORIGIN`, plus the later `NEAT_REASON`/`NEAT_CCF`/`NEAT_VAF`), the
   `NEAT_simulated_sample` VCF sample column, and the `RNEAT_generated_*` / `RNEAT_chimeric_*`
   FASTQ read-name prefixes appear in generated output, are asserted byte-for-byte in tests
   (`vcf_tools.rs`, `dup_chimeric.rs`), and are parsed downstream (`cancer_benchmark.sh` filters
   `INFO/NEAT_ORIGIN`). v2.0.0 kept them as-is for output stability; **v3.0.0 migrates all of
   them to `EIDOLON_*`** as a deliberate, separately-versioned breaking format change — bundled
   scripts and the `filter_reads` read-name parser moved in lockstep and the parity/baseline
   assertions re-blessed. Pedigree `NEAT`/`NEAT2`/`NEAT4` and internal env vars (`NEAT_DATA`,
   `NEAT_tumor_corpus`) are intentionally out of scope and unchanged.
2. **Ship a deprecation shim.** A `rneat` alias binary that prints a "renamed to `eidolon`"
   warning and forwards all args to `eidolon`, kept for at least one release as a safety net.
3. **Version: v2.0.0.** A CLI/crate identity change is breaking, even though there is no
   behavior change.
4. **Rename the library crate `common` → `eidolon-core`.** (`common` is unpublishable and
   generic; this is also the crates.io-safe published name, with the binary/conda package as
   `eidolon`.)

## Work breakdown

### Phase A — code / docs / CI (one PR → `develop`) — DONE (PR #425)
- `Cargo.toml`: bin `rneat` → `eidolon`; lib `common` → `eidolon-core`; `version` → `2.0.0`;
  update `Cargo.lock`. Add the `rneat` deprecation-shim bin target (decision 2).
- Tests: `assert_cmd::cargo_bin("rneat")` → `"eidolon"` (`cancer_parity.rs`,
  `tests/common/mod.rs`); fixture path references.
- Source: log/help/`about` strings; `RNEAT_*` → `EIDOLON_*` env-var names in `scripts/` + `tools/`.
- Leave all `NEAT_*` / `RNEAT_chimeric_*` output tokens untouched **in this phase** (decision 1);
  they were migrated later, in v3.0.0.
- CI: `rust_binaries.yml` (artifact names/paths), `rusty-neat-tests.yml` → `eidolon-tests.yml`
  (job name + filename), `thirdparty-licenses.yml`.
- Docs: `README.md` (keep the NEAT-pedigree section; add a "formerly rusty-neat / rneat" line),
  `CLAUDE.md` / `AGENTS.md`, `docs/*.md`, `conda-recipe/meta.yaml`.
- `CHANGELOG.md`: v2.0.0 rename entry.

### Phase B — GitHub repo rename
- Rename `ncsa/rusty-neat` → `ncsa/eidolon`. GitHub auto-redirects old URLs; CI and branch
  protection carry over. Update local remotes and README badge paths (the Zenodo badge is
  unchanged).

### Phase C — distribution / external
- **bioconda:** new `eidolon` feedstock (PR to bioconda-recipes); mark `rneat` deprecated.
  Update in-repo `conda-recipe/meta.yaml` (`name`, source URL, `about.home`, test commands).
  Longest external step (bioconda review latency).
- **Zenodo:** publish a new version under the *stable* concept DOI `10.5281/zenodo.20100558`,
  title "eidolon (formerly rusty-neat) …"; update `CITATION.cff` (title/name, keep authors +
  NEAT pedigree). Update the README DOI text.
- **Delta (operator env):** rename `/projects/bhrd/jallen17/rusty-neat` and
  `$SCRATCH/cargo-target/rusty-neat`; update `scripts/delta/*` and the `reference_delta_paths`
  memory.
- **Courtesy:** heads-up to the `fversaci/crusty-neat` fork owner (the one known downstream).

## Sequencing
v1.21.1 fully released (tag correct on `main`, conda recipe bumped) **→** Phase A PR **→** merge
+ v2.0.0 release **→** Phase B repo rename **→** Phase C bioconda + Zenodo **→** post-2.0.0
output-token migration (decision 1 follow-up).

## Risk notes
- The only way this goes wrong is a careless global replace nuking pedigree refs or the
  `NEAT_*` format tags. Phase A must use explicit token allow/deny lists, never `sed -i s/neat/`.
- Output parity tests stay green in v2.0.0 *because* the output tokens are unchanged; they will
  need a deliberate re-bless in the post-2.0.0 token-migration release.
