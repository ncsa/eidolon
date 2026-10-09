//! This contains only one function at the moment, a general opener and reader
use std::fs::File;
use std::io::{self, BufRead, BufReader, Lines, Read, Result, Write};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

/// Errors with `ErrorKind::AlreadyExists` if `path` is an existing regular
/// file and `allow_overwrite` is false. `label` is the human-facing name of
/// the field that produced `path` (e.g. `"frag_length.output_file"`) and is
/// embedded in the error message so the user can locate the offending
/// config entry.
///
/// Each subcommand's config parser typically wraps the returned `io::Error`
/// into its own module error via `.map_err(...)`.
pub fn check_overwrite(label: &str, path: &Path, allow_overwrite: bool) -> Result<()> {
    if !allow_overwrite && path.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "{label} already exists and overwrite_output is false: {}",
                path.display()
            ),
        ));
    }
    Ok(())
}

pub fn read_lines(filename: &PathBuf) -> Result<Lines<BufReader<File>>> {
    // This opens file and creates a buffer to read lines
    let file = File::open(filename)?;
    Ok(BufReader::new(file).lines())
}

pub fn read_gzip_lines(filename: &PathBuf) -> Result<Lines<BufReader<MultiGzDecoder<File>>>> {
    // `MultiGzDecoder` reads concatenated gzip streams, which is what `bgzip`
    // produces (one stream per block, ~64 KB each). Plain `GzDecoder` stops
    // after the first stream and silently returns no further data — that
    // path silently dropped 100% of records for every bgzipped VCF, the
    // standard format the bioinformatics ecosystem uses for indexed VCFs
    // (gnomAD, dbSNP, anything via `bcftools view -Oz` or `tabix`).
    // `MultiGzDecoder` is a strict superset: it also handles single-stream
    // plain gzip correctly, so there's no downside to using it everywhere.
    let file = File::open(filename)?;
    Ok(BufReader::new(MultiGzDecoder::new(file)).lines())
}

pub fn open_safe(filename: &PathBuf) -> Result<BufReader<File>> {
    // This opens a file for reading without compression
    let file = File::open(filename)?;
    Ok(BufReader::new(file))
}

pub fn create_output_file(filename: &PathBuf, overwrite_file: bool) -> Result<File> {
    if filename.is_file() && !overwrite_file {
        // The file already exists and we're in non overwrite mode
        panic!(
            "Attempting to overwrite an existing file: {}",
            filename.display()
        )
    } else {
        File::create(filename)
    }
}

pub fn append_to_file(filename: &PathBuf) -> Result<File> {
    // if the file doesn't exist, we'll create it. If it does, we'll append to it
    if filename.is_file() {
        File::options().append(true).open(filename)
    } else {
        File::create(filename)
    }
}

pub struct VectorBuffer {
    // The purpose of this struct is to mock a buffer so we can run single-ended read generation
    buffer: Vec<u8>,
}

impl Default for VectorBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl VectorBuffer {
    pub fn new() -> VectorBuffer {
        VectorBuffer { buffer: Vec::new() }
    }
}

impl Write for VectorBuffer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buffer.extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // In this case, flushing does nothing because
        // VectorBuffer only simulates buffering to stdout or other medium
        Ok(())
    }
}

impl Read for VectorBuffer {
    fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
        let val = self.buffer[0];
        self.buffer = Vec::from(&self.buffer[1..]);
        Ok(val as usize)
    }
}

pub fn is_gzipped_file(filename: &PathBuf) -> Result<bool> {
    let mut file = File::open(filename)?;
    let mut magic = [0u8; 2];

    let bytes_read = file.read(&mut magic)?;

    Ok(bytes_read == 2 && magic == [0x1f, 0x8b])
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;

    #[test]
    fn test_read_lines() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.txt");
        std::fs::write(&path, "line1\nline2\nline3\n").unwrap();
        let lines: Vec<String> = read_lines(&path).unwrap().map(|l| l.unwrap()).collect();
        assert_eq!(lines, vec!["line1", "line2", "line3"]);
    }

    #[test]
    fn test_read_gzip_lines() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.txt.gz");
        let file = File::create(&path).unwrap();
        let mut enc = GzEncoder::new(file, Compression::default());
        enc.write_all(b"alpha\nbeta\n").unwrap();
        enc.finish().unwrap();
        let lines: Vec<String> = read_gzip_lines(&path)
            .unwrap()
            .map(|l| l.unwrap())
            .collect();
        assert_eq!(lines, vec!["alpha", "beta"]);
    }

    #[test]
    fn test_open_safe() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.txt");
        std::fs::write(&path, "hello").unwrap();
        let mut reader = open_safe(&path).unwrap();
        let mut contents = String::new();
        reader.read_to_string(&mut contents).unwrap();
        assert_eq!(contents, "hello");
    }

    #[test]
    fn test_create_output_file_new() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("out.txt");
        assert!(!path.exists());
        let mut f = create_output_file(&path, false).unwrap();
        f.write_all(b"data").unwrap();
        drop(f);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "data");
    }

    #[test]
    fn test_create_output_file_overwrite() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("out.txt");
        std::fs::write(&path, "original").unwrap();
        let mut f = create_output_file(&path, true).unwrap();
        f.write_all(b"new").unwrap();
        drop(f);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "new");
    }

    #[test]
    #[should_panic(expected = "Attempting to overwrite an existing file")]
    fn test_create_output_file_no_overwrite_panics() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("out.txt");
        std::fs::write(&path, "original").unwrap();
        create_output_file(&path, false).unwrap();
    }

    #[test]
    fn test_append_to_file_creates_new() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("new.txt");
        assert!(!path.exists());
        let mut f = append_to_file(&path).unwrap();
        f.write_all(b"first").unwrap();
        drop(f);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "first");
    }

    #[test]
    fn test_append_to_file_appends() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("existing.txt");
        std::fs::write(&path, "line1\n").unwrap();
        let mut f = append_to_file(&path).unwrap();
        f.write_all(b"line2\n").unwrap();
        drop(f);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "line1\nline2\n");
    }

    #[test]
    fn test_vector_buffer_accumulates_writes() {
        let mut buf = VectorBuffer::new();
        buf.write_all(b"hello").unwrap();
        buf.write_all(b" world").unwrap();
        assert_eq!(buf.buffer, b"hello world");
    }
}
