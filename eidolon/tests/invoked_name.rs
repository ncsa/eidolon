//! The binary must run under whatever file name it is invoked by (#794).
//!
//! `main.rs` parses in clap's multicall mode, which takes the binary's file name as the
//! command. Release assets are named `eidolon-<target>`, so as downloaded every one of them
//! stopped with "unrecognized subcommand". These tests copy the REAL built binary to each
//! name, because the defect exists only at the file-name level: invoking `eidolon` through
//! `assert_cmd` can never see it.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The built binary copied into `dir` as `name`, executable.
fn copy_as(dir: &Path, name: &str) -> PathBuf {
    let src = assert_cmd::cargo::cargo_bin("eidolon");
    let dst = dir.join(name);
    fs::copy(&src, &dst).unwrap_or_else(|e| panic!("copy {src:?} -> {dst:?}: {e}"));
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(&dst, fs::Permissions::from_mode(0o755)).unwrap();
    }
    dst
}

fn version_under(dir: &Path, name: &str) -> String {
    let out = Command::new(copy_as(dir, name))
        .arg("--version")
        .output()
        .unwrap_or_else(|e| panic!("run as {name}: {e}"));
    assert!(
        out.status.success(),
        "run as `{name}`, `--version` failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

/// Known answer: every release asset's name, and an arbitrary one, print exactly what the
/// binary prints as `eidolon`.
#[test]
fn any_file_name_runs_as_eidolon() {
    let tmp = tempfile::tempdir().unwrap();
    let expected = version_under(tmp.path(), "eidolon");
    assert!(
        expected.starts_with("eidolon "),
        "baseline `eidolon --version` printed {expected:?}"
    );
    for name in [
        "eidolon-x86_64-unknown-linux-gnu",
        "eidolon-x86_64-unknown-linux-gnu-rhel8",
        "eidolon-aarch64-unknown-linux-gnu-rhel8",
        "eidolon-x86_64-apple-darwin",
        "foo",
    ] {
        assert_eq!(version_under(tmp.path(), name), expected, "run as `{name}`");
    }
}

/// Must-not-fire: a binary named after a subcommand still runs THAT subcommand, as multicall
/// mode already did. Falling back to `eidolon` for every name would turn `gen-reads --help`
/// into the top-level help.
#[test]
fn a_subcommand_name_still_runs_that_subcommand() {
    let tmp = tempfile::tempdir().unwrap();
    let out = Command::new(copy_as(tmp.path(), "gen-reads"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let help = String::from_utf8(out.stdout).unwrap();
    assert!(
        help.contains("Usage: gen-reads"),
        "a binary named gen-reads did not run gen-reads:\n{help}"
    );
}
