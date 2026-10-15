//! Various file tools needed throughout the code.
use log::*;
use std::path::PathBuf;
use std::{fs, io};

/// Checks that `filename`'s directory exists, creating it when `create` is set. A missing
/// directory with `create` unset is a `NotFound` error naming it, not a panic: callers pass
/// user-supplied paths.
pub fn check_parent(filename: &PathBuf, create: bool) -> io::Result<&PathBuf> {
    // A bare file name has an empty parent, meaning the current directory.
    if let Some(parent) = filename.parent().filter(|p| !p.as_os_str().is_empty())
        && !parent.exists()
    {
        if create {
            check_create_dir(&parent.to_path_buf())?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("directory {} does not exist", parent.display()),
            ));
        }
    }
    Ok(filename)
}

/// Creates `dir_to_check`, and any missing parents, if it does not exist.
pub fn check_create_dir(dir_to_check: &PathBuf) -> io::Result<()> {
    if !dir_to_check.is_dir() {
        warn!("Directory not found, creating: {:?}", dir_to_check);
        fs::create_dir_all(dir_to_check)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_parent() {
        let filename = PathBuf::from("test_data/H1N1.fa");
        assert!(PathBuf::from("test_data").is_dir());
        // An existing parent hands back the same path it was given.
        let returned = check_parent(&filename, false).unwrap();
        assert_eq!(returned, &filename);
    }

    #[test]
    fn test_check_parent_missing_without_create_is_an_error() {
        let temp_dir = tempfile::tempdir().unwrap();
        let missing = temp_dir.path().join("absent");
        let filename = missing.join("test.fa");
        assert!(!missing.exists());
        let err = check_parent(&filename, false).expect_err("a missing directory is an error");
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
        assert!(err.to_string().contains("absent"), "{err}");
        assert!(!missing.exists(), "nothing is created without `create`");
    }

    #[test]
    fn test_check_parent_creates_a_missing_directory() {
        let temp_dir = tempfile::tempdir().unwrap();
        let missing = temp_dir.path().join("nested").join("absent");
        let filename = missing.join("test.fa");
        check_parent(&filename, true).unwrap();
        assert!(missing.is_dir(), "missing parents are created too");
    }

    #[test]
    fn test_check_parent_bare_file_name_is_the_current_directory() {
        let filename = PathBuf::from("test.fa");
        assert_eq!(check_parent(&filename, false).unwrap(), &filename);
    }
}
