---
name: release
description: "Use when cutting an eidolon release or fixing one that shipped wrong — version bumps, tagging, binary assets, the conda sha256 follow-up, and merging back. Examples: \"cut v3.3.0\", \"the macOS binary is missing from the release\", \"do the post-release steps\", \"tag a patch\"."
---

# Releasing eidolon

Run these without being asked to enumerate them. The user should say "cut a release" or
"fix the release" and get the whole sequence, including the follow-ups.

**Nothing here is done until the assets are verified on the release page.** A pushed tag is
not evidence, a green workflow run is not evidence — `gh release view` listing all five
binaries is. This whole procedure exists because v3.2.0 was declared shipped while it was
missing its macOS asset.

## Which branch

| change | branch from | PR targets |
|---|---|---|
| normal feature / fix | `develop` | `develop` |
| release cut | `develop` | `main` |
| **post-release fix** | **`main`** | **`main`** |

A post-release fix is cut from `main` and carries **only** that fix. Do not merge `develop`
into `main` to deliver one — that drags in everything else that landed since, and the user
has said so explicitly. Confirm the scope with `git diff --stat origin/main..HEAD` before
opening the PR and state the file count in the PR body.

## Sequence

1. **Establish what actually changed** since the last tag, and say it plainly:
   ```bash
   git diff --stat <lasttag>..HEAD -- eidolon/src eidolon-core/src
   ```
   If that is empty, the binary is unchanged and the release notes must say so — otherwise a
   patch release reads as a behaviour fix.

2. **Bump the version in both places**, then refresh the lockfile:
   ```bash
   sed -i "s/^version = '<old>'/version = '<new>'/" Cargo.toml        # [workspace.package]
   sed -i 's/{% set version = "<old>" %}/{% set version = "<new>" %}/' conda-recipe/meta.yaml
   cargo check --workspace -q     # rewrites Cargo.lock; do not hand-edit it
   ```

3. **CHANGELOG entry at the top**, matching the existing format exactly — a bare `M/D/YYYY`
   line, `=========`, then `## eidolon vX.Y.Z — <short title>`. Lead with what changed for a
   user. If the binary is unchanged, that is the first sentence.

4. **PR, wait for green, and let the user merge.** Merging is theirs: `--admin` bypasses
   branch protection and is blocked in this environment anyway. `release-blockers` runs only
   on `main`-targeted PRs — on a develop PR it reports `SKIPPED`, which is not a pass.

5. **Tag from `main` after the merge**, never from the branch:
   ```bash
   git checkout main && git pull --ff-only
   git tag -a vX.Y.Z -m "vX.Y.Z — <title>" && git push origin vX.Y.Z
   ```

6. **Watch the build and verify every asset.** Five are expected: linux-gnu, linux-gnu-rhel8,
   aarch64-rhel8, apple-darwin, windows-msvc.exe.
   ```bash
   gh run list --workflow rust_binaries.yml --limit 1
   gh release view vX.Y.Z --json assets -q '[.assets[].name] | join("\n")'
   ```
   Fewer than five is a failed release, not a partial success. `fail-fast: false` is
   deliberate so one broken target does not take out the others — which also means the run
   can go green-ish while an asset is missing. Count them.

7. **Conda sha256 follow-up** — the tarball only exists once the tag is pushed, so this is
   always a second PR:
   ```bash
   curl -sL https://github.com/ncsa/eidolon/archive/refs/tags/vX.Y.Z.tar.gz | sha256sum
   ```
   Put it in `conda-recipe/meta.yaml`, PR to `main`.

8. **Merge back so `develop` has the bump**, or its next release cut starts from a stale
   version.

9. **Report** the tag, the asset count, and anything not verified.

## Gotchas that have actually bitten

- **A tag cannot be rebuilt into correctness.** Actions runs the workflow file *from the ref
  being built*, so re-pushing a tag re-runs that tag's workflow, bugs included.
  `workflow_dispatch` does not rescue it either: the trigger must exist both on the default
  branch (`main`) and at the dispatched ref. A broken release needs a new tag.
- **`dtolnay/rust-toolchain@stable` sets `RUSTUP_TOOLCHAIN` and silently overrides
  `rust-toolchain.toml`.** It installs targets into stable while cargo follows the pin —
  which is exactly how v3.2.0 lost its macOS binary (`can't find crate for core`).
- **`gh pr edit` is broken on this repo** (deprecated Projects-classic GraphQL). Patch a body
  with `gh api -X PATCH repos/ncsa/eidolon/pulls/<N> -f body=…`. `gh issue edit` works.
- **`gh pr view --json commits` has served stale commit lists.** Use
  `gh api repos/ncsa/eidolon/pulls/<N>/commits` to check whether a push landed.
- **Adding a workspace member changes `THIRDPARTY.yml`'s `root_name`** and turns
  `check-thirdparty` red. Regenerate rather than hand-edit:
  ```bash
  cargo-bundle-licenses --format yaml --output THIRDPARTY.yml --previous THIRDPARTY.yml
  ```
- **Confirm the merge landed** before tagging: `git merge-base --is-ancestor <sha> origin/main`.
