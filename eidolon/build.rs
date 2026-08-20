//! Stamp the build's git commit into the binary.
//!
//! WHY THIS EXISTS: `eidolon --version` reported only `CARGO_PKG_VERSION`, and the Delta
//! validation pipeline's stale-binary guard compared that against `Cargo.toml`. Almost all
//! work in this repo lands without a version bump — 3.1.0 has carried a dozen Rust commits
//! and was never even tagged — so the common case, *a binary built from an older commit of
//! the same version*, passed the guard silently. #513. The guard's own comment names the
//! incident: job 20682989 reproduced job 20675480's numbers exactly and the log could not
//! answer "was the fix even in this build?".
//!
//! Deliberately shells out to `git` rather than adding a crate (`vergen` et al). CLAUDE.md
//! requires the shipped artifact stay pure Rust with no runtime requirements and a
//! `Cargo.lock` free of extra dependencies; a build script invoking a subprocess adds
//! neither.
//!
//! Degrades to "unknown" when git is unavailable or `.git` is absent — a conda build from a
//! source tarball, for instance. It is then the CONSUMER's job to decide whether an
//! unverifiable build is acceptable; the validation pipeline refuses one, because "cannot
//! verify" must not read as "verified".

use std::process::Command;

fn main() {
    // Re-run when HEAD moves or the index changes, so the stamp cannot go stale within a
    // working tree. Both may be absent (tarball build); cargo tolerates missing paths here.
    println!("cargo:rerun-if-changed=../.git/HEAD");
    println!("cargo:rerun-if-changed=../.git/index");
    println!("cargo:rerun-if-env-changed=EIDOLON_GIT_SHA");

    // An explicit override wins, so a packager can stamp a known provenance without git.
    if let Ok(forced) = std::env::var("EIDOLON_GIT_SHA") {
        if !forced.trim().is_empty() {
            println!("cargo:rustc-env=EIDOLON_GIT_SHA={}", forced.trim());
            return;
        }
    }

    println!("cargo:rustc-env=EIDOLON_GIT_SHA={}", git_describe());
}

fn git_describe() -> String {
    let sha = match run(&["rev-parse", "--short=7", "HEAD"]) {
        Some(s) if !s.is_empty() => s,
        // No git, no repo, or a tarball. Say so rather than inventing provenance.
        _ => return "unknown".to_string(),
    };
    // A dirty tree means the SHA does NOT describe what was compiled. Saying so is the
    // whole point: a guard that accepts a clean-looking SHA for uncommitted code is the
    // same false-pass shape this stamp exists to prevent.
    match run(&["status", "--porcelain", "--untracked-files=no"]) {
        Some(s) if !s.is_empty() => format!("{sha}-dirty"),
        _ => sha,
    }
}

fn run(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}
