//! Combine per-chunk temp FASTQs into the final output in random order (#870).
//!
//! **Why.** gen-reads generates reads chunk by chunk, so concatenating the temp files gives a
//! FASTQ in genome order, which no sequencer produces. Aligners read FASTQs in batches, and a
//! batch filled with one region's reads behaves differently from a mixed one: bwa-mem2 needed
//! 166 GB on a 21x T2T-CHM13 simulation as written and 90 GB on the same reads shuffled (#870).
//! The original in-memory shuffle was removed because it held every read at once.
//!
//! **How: a two-pass external shuffle with bounded memory.**
//! 1. Scatter. Read the temp files in order, and send each record (each pair, when paired) to one
//!    of `K` bucket files chosen at random. Both mates of a pair go to the same bucket, in the
//!    same position. Each input file is deleted once read, so temp disk stays near one copy.
//! 2. Shuffle. Load one bucket at a time, permute it in memory, and write it as gzip; then
//!    concatenate the buckets (gzip members concatenate into one valid stream).
//!
//! Bucket assignment and then a uniform permutation within each bucket give a uniform
//! permutation of the whole file. `K` is chosen so a bucket fits the memory budget.
//!
//! **Determinism.** The scatter draws from one RNG, sequentially. Each bucket's permutation draws
//! from its own child RNG. Threads only change which bucket is processed when, so output is
//! byte-identical across thread counts.

use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::read::MultiGzDecoder;
use libdeflater::{CompressionLvl, Compressor};

use crate::file_tools::block_gz::BlockGzWriter;
use crate::file_tools::fastq_tools::FastqToolsError;
use crate::rng::NeatRng;

/// Uncompressed bytes per bucket, both mates together. Pass 2 holds this much per thread,
/// roughly doubled by per-record allocation.
const BUCKET_TARGET_BYTES: u64 = 256 * 1024 * 1024;
/// Upper bound on buckets. The scatter holds one buffer per bucket per mate, so this bounds that
/// memory (512 * 2 * `BUCKET_BUFFER_BYTES` = 256 MiB). Past it, buckets grow instead.
const MAX_BUCKETS: usize = 512;
/// Scatter buffer per bucket per mate. A full buffer is compressed as one gzip member and
/// appended to its bucket file, which is then closed again, so the scatter holds at most one
/// bucket file open at any bucket count. One writer per bucket per mate would need
/// 2 * MAX_BUCKETS descriptors and failed under `ulimit -n` 256 or 512.
const BUCKET_BUFFER_BYTES: usize = 256 * 1024;
/// Gzipped FASTQ expands about 3.5-4.5x; overestimating only adds buckets.
const EXPANSION_ESTIMATE: u64 = 5;
/// Pass-2 threads, capped because each holds a whole bucket in memory.
const MAX_SHUFFLE_THREADS: usize = 8;
/// Pass 2 runs at most as many buckets at once as fit in this much memory, judged by the
/// largest bucket. Small runs get full parallelism; at genome scale, where `MAX_BUCKETS` makes
/// buckets larger, concurrency drops so memory stays near this figure.
const SHUFFLE_MEMORY_BUDGET: u64 = 4 * 1024 * 1024 * 1024;

fn read_err(path: &Path, e: impl std::fmt::Display) -> FastqToolsError {
    FastqToolsError::FastqReadError(format!("{}: {e}", path.display()))
}

fn write_err(path: &Path, e: impl std::fmt::Display) -> FastqToolsError {
    FastqToolsError::FastqWriteError(format!("{}: {e}", path.display()))
}

/// Reads one 4-line FASTQ record into `buf` (appending). `Ok(false)` at a clean end of input;
/// an error if the input ends partway through a record.
fn read_record(
    reader: &mut impl BufRead,
    buf: &mut Vec<u8>,
    path: &Path,
) -> Result<bool, FastqToolsError> {
    let start = buf.len();
    for line in 0..4 {
        let n = reader
            .read_until(b'\n', buf)
            .map_err(|e| read_err(path, e))?;
        if n == 0 {
            if line == 0 {
                return Ok(false);
            }
            return Err(read_err(
                path,
                format!("input ends {line} lines into a FASTQ record"),
            ));
        }
    }
    if buf.last() != Some(&b'\n') {
        buf.push(b'\n');
    }
    if buf[start] != b'@' {
        return Err(read_err(path, "a FASTQ record does not start with '@'"));
    }
    Ok(true)
}

fn open_gz(path: &Path) -> Result<BufReader<MultiGzDecoder<File>>, FastqToolsError> {
    let file = File::open(path).map_err(|e| read_err(path, e))?;
    Ok(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
}

/// One bucket file being scattered into: records accumulate in `buf` and leave as a gzip member
/// appended to `path`, so the file is open only while a member is written.
struct BucketSink {
    path: PathBuf,
    buf: Vec<u8>,
    wrote: bool,
}

impl BucketSink {
    fn new(path: PathBuf) -> Self {
        Self {
            path,
            buf: Vec::new(),
            wrote: false,
        }
    }

    fn push(
        &mut self,
        rec: &[u8],
        c: &mut Compressor,
        scratch: &mut Vec<u8>,
    ) -> Result<(), FastqToolsError> {
        self.buf.extend_from_slice(rec);
        if self.buf.len() >= BUCKET_BUFFER_BYTES {
            self.flush(c, scratch)?;
        }
        Ok(())
    }

    /// Appends the buffer as one gzip member. An empty member is written only for a bucket that
    /// received nothing, so every bucket file is valid gzip.
    fn flush(&mut self, c: &mut Compressor, scratch: &mut Vec<u8>) -> Result<(), FastqToolsError> {
        if self.buf.is_empty() && self.wrote {
            return Ok(());
        }
        let bound = c.gzip_compress_bound(self.buf.len());
        if scratch.len() < bound {
            scratch.resize(bound, 0);
        }
        let n = c
            .gzip_compress(&self.buf, scratch)
            .map_err(|e| write_err(&self.path, format!("gzip compression failed: {e:?}")))?;
        let mut f = fs::OpenOptions::new()
            .append(true)
            .create(true)
            .open(&self.path)
            .map_err(|e| write_err(&self.path, e))?;
        f.write_all(&scratch[..n])
            .map_err(|e| write_err(&self.path, e))?;
        self.buf.clear();
        self.wrote = true;
        Ok(())
    }

    fn finish(&mut self, c: &mut Compressor, scratch: &mut Vec<u8>) -> Result<(), FastqToolsError> {
        self.flush(c, scratch)?;
        self.buf = Vec::new();
        Ok(())
    }
}

/// Pass-2 concurrency: the requested threads, capped at `MAX_SHUFFLE_THREADS` and at however
/// many of the largest bucket fit in `SHUFFLE_MEMORY_BUDGET`. Always at least one.
fn shuffle_threads(requested: usize, largest_bucket_bytes: u64) -> usize {
    let fit = (SHUFFLE_MEMORY_BUDGET / largest_bucket_bytes.max(1)) as usize;
    requested.min(MAX_SHUFFLE_THREADS).min(fit).max(1)
}

/// Number of buckets for inputs totalling `compressed_bytes` on disk.
fn bucket_count(compressed_bytes: u64) -> usize {
    let estimate = compressed_bytes.saturating_mul(EXPANSION_ESTIMATE);
    (estimate.div_ceil(BUCKET_TARGET_BYTES) as usize).clamp(1, MAX_BUCKETS)
}

/// Writes `files_r1` (and `files_r2`, record for record, when paired) to `out_r1` (and
/// `out_r2`) in a uniformly random order, keeping mates in register. `rng` is consumed for the
/// scatter; each bucket's permutation uses `rng.derive_child(bucket)`. Temp files go in
/// `work_dir`, and every input file is deleted once read.
pub fn combine_temp_fastqs_shuffled(
    files_r1: Vec<PathBuf>,
    files_r2: Vec<PathBuf>,
    out_r1: &Path,
    out_r2: Option<&Path>,
    rng: &mut NeatRng,
    work_dir: &Path,
    threads: usize,
) -> Result<(), FastqToolsError> {
    let compressed: u64 = files_r1
        .iter()
        .chain(files_r2.iter())
        .map(|p| fs::metadata(p).map(|m| m.len()).unwrap_or(0))
        .sum();
    shuffle_into_buckets(
        files_r1,
        files_r2,
        out_r1,
        out_r2,
        rng,
        work_dir,
        threads,
        bucket_count(compressed),
    )
}

/// The shuffle with an explicit bucket count, so tests can exercise the multi-bucket path on
/// inputs small enough that `bucket_count` would choose one.
#[allow(clippy::too_many_arguments)]
fn shuffle_into_buckets(
    files_r1: Vec<PathBuf>,
    files_r2: Vec<PathBuf>,
    out_r1: &Path,
    out_r2: Option<&Path>,
    rng: &mut NeatRng,
    work_dir: &Path,
    threads: usize,
    k: usize,
) -> Result<(), FastqToolsError> {
    let paired = out_r2.is_some();
    if paired && files_r1.len() != files_r2.len() {
        return Err(FastqToolsError::FastqWriteError(format!(
            "shuffle: {} read-1 temp files but {} read-2",
            files_r1.len(),
            files_r2.len()
        )));
    }
    let bucket_dir = work_dir.join("shuffle_buckets");
    fs::create_dir_all(&bucket_dir).map_err(|e| write_err(&bucket_dir, e))?;
    let bucket_path = |b: usize, mate: u8, stage: &str| {
        bucket_dir.join(format!("{stage}_{b:04}_r{mate}.fastq.gz"))
    };

    // ── Pass 1: scatter ───────────────────────────────────────────────────────────────────────
    let mates: Vec<u8> = if paired { vec![1, 2] } else { vec![1] };
    let mut sinks: Vec<Vec<BucketSink>> = (0..k)
        .map(|b| {
            mates
                .iter()
                .map(|&m| BucketSink::new(bucket_path(b, m, "bucket")))
                .collect()
        })
        .collect();
    let mut compressor = Compressor::new(CompressionLvl::fastest());
    let mut scratch = Vec::new();
    let mut counts = vec![0u64; k];
    let mut bytes = vec![0u64; k];
    let (mut rec1, mut rec2) = (Vec::with_capacity(1024), Vec::with_capacity(1024));
    for (i, p1) in files_r1.iter().enumerate() {
        let mut in1 = open_gz(p1)?;
        let mut in2 = if paired {
            Some((open_gz(&files_r2[i])?, &files_r2[i]))
        } else {
            None
        };
        loop {
            rec1.clear();
            let more1 = read_record(&mut in1, &mut rec1, p1)?;
            let more2 = match in2.as_mut() {
                Some((r, p)) => {
                    rec2.clear();
                    read_record(r, &mut rec2, p)?
                }
                None => more1,
            };
            if more1 != more2 {
                return Err(read_err(
                    p1,
                    "read-1 and read-2 temp files hold different numbers of records",
                ));
            }
            if !more1 {
                break;
            }
            let b = ((rng.random()? * k as f64) as usize).min(k - 1);
            sinks[b][0].push(&rec1, &mut compressor, &mut scratch)?;
            if paired {
                sinks[b][1].push(&rec2, &mut compressor, &mut scratch)?;
            }
            counts[b] += 1;
            bytes[b] += (rec1.len() + if paired { rec2.len() } else { 0 }) as u64;
        }
        drop(in1);
        drop(in2);
        let _ = fs::remove_file(p1);
        if paired {
            let _ = fs::remove_file(&files_r2[i]);
        }
    }
    for sink in sinks.iter_mut().flatten() {
        sink.finish(&mut compressor, &mut scratch)?;
    }
    drop(sinks);

    // ── Pass 2: shuffle each bucket ───────────────────────────────────────────────────────────
    let next = AtomicUsize::new(0);
    let first_error: Mutex<Option<FastqToolsError>> = Mutex::new(None);
    let base = *rng;
    let shuffle_bucket = |b: usize| -> Result<(), FastqToolsError> {
        // One contiguous buffer per mate, with record start offsets: a Vec per record would
        // cost an allocation and its slack for every read.
        let mut records: Vec<(Vec<u8>, Vec<usize>)> = Vec::with_capacity(mates.len());
        for &m in &mates {
            let p = bucket_path(b, m, "bucket");
            let mut r = open_gz(&p)?;
            let mut data = Vec::new();
            let mut starts = Vec::with_capacity(counts[b] as usize + 1);
            loop {
                starts.push(data.len());
                if !read_record(&mut r, &mut data, &p)? {
                    break;
                }
            }
            if starts.len() as u64 - 1 != counts[b] {
                return Err(read_err(
                    &p,
                    format!(
                        "bucket holds {} records, {} were written",
                        starts.len() - 1,
                        counts[b]
                    ),
                ));
            }
            drop(r);
            let _ = fs::remove_file(&p);
            records.push((data, starts));
        }
        let mut order: Vec<usize> = (0..counts[b] as usize).collect();
        base.derive_child(b as u64).shuffle_in_place(&mut order)?;
        for (j, &m) in mates.iter().enumerate() {
            let p = bucket_path(b, m, "shuffled");
            let f = File::create(&p).map_err(|e| write_err(&p, e))?;
            let mut w = BlockGzWriter::new(BufWriter::new(f));
            let (data, starts) = &records[j];
            for &i in &order {
                w.write_all(&data[starts[i]..starts[i + 1]])
                    .map_err(|e| write_err(&p, e))?;
            }
            w.finish().map_err(|e| write_err(&p, e))?;
        }
        Ok(())
    };
    std::thread::scope(|s| {
        for _ in 0..shuffle_threads(threads, bytes.iter().copied().max().unwrap_or(0)) {
            s.spawn(|| {
                loop {
                    let b = next.fetch_add(1, Ordering::SeqCst);
                    if b >= k || first_error.lock().unwrap().is_some() {
                        break;
                    }
                    if let Err(e) = shuffle_bucket(b) {
                        first_error.lock().unwrap().get_or_insert(e);
                        break;
                    }
                }
            });
        }
    });
    if let Some(e) = first_error.into_inner().unwrap() {
        return Err(e);
    }

    // ── Concatenate, in bucket order ──────────────────────────────────────────────────────────
    for (j, &m) in mates.iter().enumerate() {
        let out = if j == 0 {
            out_r1
        } else {
            out_r2.expect("paired")
        };
        let mut w = BufWriter::new(
            fs::OpenOptions::new()
                .append(true)
                .create(true)
                .open(out)
                .map_err(|e| write_err(out, e))?,
        );
        for b in 0..k {
            let p = bucket_path(b, m, "shuffled");
            let mut f = BufReader::new(File::open(&p).map_err(|e| read_err(&p, e))?);
            std::io::copy(&mut f, &mut w).map_err(|e| write_err(out, e))?;
            drop(f);
            let _ = fs::remove_file(&p);
        }
        w.flush().map_err(|e| write_err(out, e))?;
    }
    let _ = fs::remove_dir(&bucket_dir);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    /// Writes `n` records named by index (and their mates) as `chunks` temp-file pairs, in order.
    fn temp_inputs(
        dir: &Path,
        n: usize,
        chunks: usize,
        paired: bool,
    ) -> (Vec<PathBuf>, Vec<PathBuf>) {
        let (mut f1, mut f2) = (Vec::new(), Vec::new());
        let per = n.div_ceil(chunks);
        for c in 0..chunks {
            let p1 = dir.join(format!("t{c}_r1.fastq.gz"));
            let p2 = dir.join(format!("t{c}_r2.fastq.gz"));
            let mut w1 = BlockGzWriter::new(File::create(&p1).unwrap());
            let mut w2 = BlockGzWriter::new(File::create(&p2).unwrap());
            for i in c * per..((c + 1) * per).min(n) {
                write!(w1, "@read_{i:08}/1\nACGT{i}\n+\nIIII{i}\n").unwrap();
                write!(w2, "@read_{i:08}/2\nTGCA{i}\n+\nJJJJ{i}\n").unwrap();
            }
            w1.finish().unwrap();
            w2.finish().unwrap();
            f1.push(p1);
            if paired {
                f2.push(p2);
            }
        }
        (f1, f2)
    }

    fn names(path: &Path) -> Vec<String> {
        let mut s = String::new();
        MultiGzDecoder::new(File::open(path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        s.lines().step_by(4).map(str::to_string).collect()
    }

    fn index_of(name: &str) -> usize {
        name[6..14].parse().unwrap()
    }

    /// Shuffles `n` records from `chunks` temp-file pairs through `k` buckets. Tests pass k > 1
    /// explicitly: at these sizes `bucket_count` would always choose 1, leaving the scatter untested.
    fn run_k(
        n: usize,
        chunks: usize,
        seed: &str,
        threads: usize,
        k: usize,
    ) -> (Vec<String>, Vec<String>) {
        let dir = tempfile::tempdir().unwrap();
        let (f1, f2) = temp_inputs(dir.path(), n, chunks, true);
        let (o1, o2) = (
            dir.path().join("o1.fastq.gz"),
            dir.path().join("o2.fastq.gz"),
        );
        let mut rng = NeatRng::new_from_seed(&vec![seed.to_string()]).unwrap();
        shuffle_into_buckets(f1, f2, &o1, Some(&o2), &mut rng, dir.path(), threads, k).unwrap();
        (names(&o1), names(&o2))
    }

    /// The output files' raw bytes, for comparing gzip streams rather than their contents.
    fn raw_k(n: usize, chunks: usize, seed: &str, threads: usize, k: usize) -> (Vec<u8>, Vec<u8>) {
        let dir = tempfile::tempdir().unwrap();
        let (f1, f2) = temp_inputs(dir.path(), n, chunks, true);
        let (o1, o2) = (
            dir.path().join("o1.fastq.gz"),
            dir.path().join("o2.fastq.gz"),
        );
        let mut rng = NeatRng::new_from_seed(&vec![seed.to_string()]).unwrap();
        shuffle_into_buckets(f1, f2, &o1, Some(&o2), &mut rng, dir.path(), threads, k).unwrap();
        (fs::read(&o1).unwrap(), fs::read(&o2).unwrap())
    }

    fn run(n: usize, chunks: usize, seed: &str, threads: usize) -> (Vec<String>, Vec<String>) {
        run_k(n, chunks, seed, threads, 16)
    }

    #[test]
    fn output_is_a_permutation_with_mates_in_register() {
        let n = 20_000;
        let (r1, r2) = run(n, 7, "perm", 4);
        assert_eq!(r1.len(), n);
        let mut idx: Vec<usize> = r1.iter().map(|s| index_of(s)).collect();
        for (a, b) in r1.iter().zip(&r2) {
            assert_eq!(
                a.trim_end_matches("/1"),
                b.trim_end_matches("/2"),
                "mates out of register"
            );
        }
        idx.sort_unstable();
        assert_eq!(idx, (0..n).collect::<Vec<_>>(), "every record exactly once");
    }

    #[test]
    fn output_order_is_random_not_the_input_order() {
        // Adjacent records: in input order every neighbor is the next index. Shuffled, a neighbor
        // is within 100 indexes with probability ~ 2 * 100 / n = 0.01.
        let n = 20_000;
        let (r1, _) = run(n, 7, "order", 4);
        let idx: Vec<i64> = r1.iter().map(|s| index_of(s) as i64).collect();
        let near = idx
            .windows(2)
            .filter(|w| (w[1] - w[0]).abs() <= 100)
            .count() as f64
            / (n - 1) as f64;
        assert!(
            near < 0.03,
            "adjacent records within 100 indexes: {near:.4}; the output is not shuffled"
        );
        assert_ne!(
            idx[..100],
            (0..100).collect::<Vec<i64>>()[..],
            "starts in input order"
        );
    }

    #[test]
    fn same_seed_same_output_across_thread_counts_and_different_seed_differs() {
        let n = 5_000;
        let a = run(n, 3, "det", 1);
        assert_eq!(a, run(n, 3, "det", 8), "thread count changed the output");
        // The gzip bytes too: member boundaries or headers that varied with scheduling would
        // pass the record comparison above.
        assert!(
            raw_k(n, 3, "det", 1, 16) == raw_k(n, 3, "det", 8, 16),
            "thread count changed the gzip bytes"
        );
        assert_ne!(
            a.0,
            run(n, 3, "det2", 1).0,
            "a different seed gave the same order"
        );
    }

    #[test]
    fn single_ended_and_empty_inputs() {
        let dir = tempfile::tempdir().unwrap();
        let (f1, _) = temp_inputs(dir.path(), 1_000, 2, false);
        let o1 = dir.path().join("o1.fastq.gz");
        let mut rng = NeatRng::new_from_seed(&vec!["se".to_string()]).unwrap();
        combine_temp_fastqs_shuffled(f1, vec![], &o1, None, &mut rng, dir.path(), 2).unwrap();
        assert_eq!(names(&o1).len(), 1_000);

        let (e1, e2) = temp_inputs(dir.path(), 0, 1, true);
        let (o1, o2) = (
            dir.path().join("e1.fastq.gz"),
            dir.path().join("e2.fastq.gz"),
        );
        combine_temp_fastqs_shuffled(e1, e2, &o1, Some(&o2), &mut rng, dir.path(), 2).unwrap();
        assert!(names(&o1).is_empty() && names(&o2).is_empty());
    }

    #[test]
    fn unequal_mate_files_are_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let (f1, f2) = temp_inputs(dir.path(), 100, 1, true);
        // Drop the last read-2 record.
        let mut s = String::new();
        MultiGzDecoder::new(File::open(&f2[0]).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        let kept: Vec<&str> = s.lines().collect();
        let mut w = BlockGzWriter::new(File::create(&f2[0]).unwrap());
        for l in &kept[..kept.len() - 4] {
            writeln!(w, "{l}").unwrap();
        }
        w.finish().unwrap();
        let (o1, o2) = (
            dir.path().join("o1.fastq.gz"),
            dir.path().join("o2.fastq.gz"),
        );
        let mut rng = NeatRng::new_from_seed(&vec!["x".to_string()]).unwrap();
        let err = combine_temp_fastqs_shuffled(f1, f2, &o1, Some(&o2), &mut rng, dir.path(), 1)
            .unwrap_err();
        assert!(
            err.to_string().contains("different numbers of records"),
            "{err}"
        );
    }

    #[test]
    fn one_bucket_and_many_buckets_both_give_a_shuffled_permutation() {
        // k = 1 is all in-bucket shuffle; k = 64 is mostly scatter. Both must permute and mix.
        for k in [1, 64] {
            let n = 10_000;
            let (r1, r2) = run_k(n, 5, "k", 3, k);
            let mut idx: Vec<i64> = r1.iter().map(|s| index_of(s) as i64).collect();
            let near = idx
                .windows(2)
                .filter(|w| (w[1] - w[0]).abs() <= 100)
                .count() as f64
                / (n - 1) as f64;
            assert!(near < 0.05, "k={k}: adjacent within 100: {near:.4}");
            assert!(
                r1.iter()
                    .zip(&r2)
                    .all(|(a, b)| a.trim_end_matches("/1") == b.trim_end_matches("/2")),
                "k={k}"
            );
            idx.sort_unstable();
            assert_eq!(idx, (0..n as i64).collect::<Vec<_>>(), "k={k}");
        }
    }

    #[test]
    fn shuffle_threads_fit_the_memory_budget() {
        let gib = 1024 * 1024 * 1024;
        assert_eq!(
            shuffle_threads(64, 256 * 1024 * 1024),
            MAX_SHUFFLE_THREADS,
            "small buckets: full cap"
        );
        assert_eq!(shuffle_threads(64, gib), 4, "1 GiB buckets: 4 fit in 4 GiB");
        assert_eq!(
            shuffle_threads(64, 10 * gib),
            1,
            "a bucket over budget still runs, alone"
        );
        assert_eq!(shuffle_threads(2, 1), 2, "never more than requested");
        assert_eq!(shuffle_threads(0, 0), 1, "at least one");
    }

    const FD_CHILD_ENV: &str = "EIDOLON_SHUFFLE_FD_CHILD";
    /// A file-descriptor limit well under what one writer per bucket per mate would need at the
    /// cap (2 * MAX_BUCKETS), and well over what the shuffle needs.
    const FD_LIMIT: usize = 48;

    /// Runs only as the child of `cap_path_runs_under_a_low_file_descriptor_limit`, which sets
    /// the limit; a no-op otherwise.
    #[test]
    fn cap_path_child() {
        if std::env::var(FD_CHILD_ENV).is_err() {
            return;
        }
        let n = 3_000;
        let (r1, r2) = run_k(n, 4, "fd", MAX_SHUFFLE_THREADS, MAX_BUCKETS);
        assert_eq!(r1.len(), n);
        assert!(
            r1.iter()
                .zip(&r2)
                .all(|(a, b)| a.trim_end_matches("/1") == b.trim_end_matches("/2"))
        );
        println!("CAP-PATH-OK {n} pairs through {MAX_BUCKETS} buckets");
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn cap_path_runs_under_a_low_file_descriptor_limit() {
        // Paired, at MAX_BUCKETS, under `ulimit -n 48`. The child prints a marker so a filter
        // that matched nothing cannot pass for a run that succeeded.
        let exe = std::env::current_exe().unwrap();
        let out = std::process::Command::new("bash")
            .arg("-c")
            .arg(format!(
                "ulimit -n {FD_LIMIT} && exec \"$0\" --exact file_tools::fastq_shuffle::tests::cap_path_child --nocapture --test-threads=1"
            ))
            .arg(&exe)
            .env(FD_CHILD_ENV, "1")
            .output()
            .unwrap();
        let (stdout, stderr) = (
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr),
        );
        assert!(
            out.status.success() && stdout.contains("CAP-PATH-OK"),
            "shuffle at {MAX_BUCKETS} buckets failed under ulimit -n {FD_LIMIT}:\n{stdout}\n{stderr}"
        );
    }

    #[test]
    fn bucket_count_scales_with_input_and_is_capped() {
        assert_eq!(bucket_count(0), 1);
        assert_eq!(bucket_count(1_000), 1);
        assert_eq!(
            bucket_count(BUCKET_TARGET_BYTES),
            EXPANSION_ESTIMATE as usize
        );
        assert_eq!(bucket_count(u64::MAX), MAX_BUCKETS);
    }
}
