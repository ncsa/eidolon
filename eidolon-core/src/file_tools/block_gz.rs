//! Block-wise gzip writer backed by libdeflate.
//!
//! libdeflate compresses substantially faster than streaming zlib/zlib-ng for
//! the same gzip output, but its API is one-shot (compress a whole buffer),
//! not streaming. `BlockGzWriter` bridges that: it buffers incoming bytes and,
//! whenever it has accumulated ~`BLOCK_SIZE`, compresses that block into an
//! independent gzip *member* and writes it out. Concatenated gzip members form
//! a standard multi-member gzip stream — `gunzip`, `zcat`, `MultiGzDecoder`,
//! and every common bioinformatics tool read them transparently.
//!
//! Memory is bounded to one block regardless of total output size, preserving
//! eidolon's streaming, low-memory behaviour while cutting compression CPU.
//!
//! ## Why this and not BGZF / `bgzip`?
//!
//! This is intentionally *plain* multi-member gzip (à la `pigz`), not BGZF —
//! even though eidolon already uses BGZF (`noodles::bgzf`) for BAM and `.vcf.gz`.
//! BGZF adds a `BC` extra-field block-size header, caps blocks at 64 KiB, and
//! writes an EOF marker, all to support **random-access indexing** (`.gzi`,
//! tabix). FASTQ output needs none of that — it is streamed start-to-finish by
//! aligners, and FASTQ in the wild is conventionally plain gzip. Choosing plain
//! gzip here also keeps two things simple that BGZF would complicate:
//!
//! 1. **The combine step can byte-concatenate temp files.** Plain gzip members
//!    concatenate into a valid stream with no ceremony; BGZF temps each carry an
//!    EOF marker that would land mid-stream on concatenation (readable, but not
//!    clean BGZF — you'd have to strip them, as `samtools cat` does for BAM).
//! 2. **No competing thread pool.** `noodles`'s multithreaded BGZF writer spawns
//!    its own compression threads, which would oversubscribe against eidolon's
//!    per-contig rayon parallelism (and read generation is memory-bandwidth
//!    bound, so extra compression threads don't help). This writer compresses
//!    inline within the calling thread.
//!
//! So this is not a re-implementation of `bgzip`; it's the smaller, format-
//! appropriate tool for FASTQ. Use `noodles::bgzf` when you need an indexable
//! BGZF stream (BAM / tabix VCF).

use libdeflater::{CompressionLvl, Compressor};
use log::error;
use std::io::{self, Write};

/// Uncompressed bytes buffered before a block is emitted. 256 KiB balances
/// compression ratio (bigger = slightly better) against memory and per-member
/// header overhead (smaller = more gzip headers).
const BLOCK_SIZE: usize = 256 * 1024;

/// A `Write` adapter that emits multi-member gzip via libdeflate.
///
/// Call [`BlockGzWriter::finish`] to flush the final partial block and observe
/// any error. As a safety net `Drop` also flushes (matching `flate2::GzEncoder`,
/// whose `Drop` likewise finalizes and swallows errors), so existing code that
/// relies on drop-finalization keeps working.
pub struct BlockGzWriter<W: Write> {
    inner: W,
    buf: Vec<u8>,
    /// Reused scratch buffer for compressed output, so we don't allocate a fresh
    /// buffer per block (which adds allocator pressure that hurts multi-thread).
    out: Vec<u8>,
    compressor: Compressor,
    finished: bool,
    /// Whether any gzip member has been emitted. If none has by finalize time,
    /// we emit an empty member so the output is still a valid (empty) gzip
    /// stream — matching `GzEncoder`, and avoiding 0-byte files that a decoder
    /// would reject with "unexpected end of file" (e.g. a fully deleted contig).
    wrote_any: bool,
}

impl<W: Write> BlockGzWriter<W> {
    pub fn new(inner: W) -> Self {
        // CompressionLvl::default() is libdeflate level 6 — the same ratio
        // target as flate2's Compression::default(), so output size is
        // comparable and only the speed changes.
        Self {
            inner,
            buf: Vec::with_capacity(BLOCK_SIZE),
            out: Vec::new(),
            compressor: Compressor::new(CompressionLvl::default()),
            finished: false,
            wrote_any: false,
        }
    }

    /// Compress `block` (which may be empty) as one gzip member and write it.
    fn emit_member(&mut self, block: &[u8]) -> io::Result<()> {
        let bound = self.compressor.gzip_compress_bound(block.len());
        if self.out.len() < bound {
            self.out.resize(bound, 0);
        }
        let n = self
            .compressor
            .gzip_compress(block, &mut self.out)
            .map_err(|e| io::Error::other(format!("libdeflate gzip compression failed: {e:?}")))?;
        self.inner.write_all(&self.out[..n])?;
        self.wrote_any = true;
        Ok(())
    }

    /// Emit the buffered remainder, and — if nothing was ever written — an empty
    /// member, so the stream is always valid gzip.
    fn flush_buffered(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            let block = std::mem::take(&mut self.buf);
            self.emit_member(&block)?;
        }
        if !self.wrote_any {
            self.emit_member(&[])?;
        }
        Ok(())
    }

    /// Flush the final partial block and the inner writer. After this, `Drop`
    /// is a no-op.
    pub fn finish(mut self) -> io::Result<()> {
        self.flush_buffered()?;
        self.finished = true;
        self.inner.flush()
        // self drops here; Drop sees `finished` and does nothing.
    }
}

impl<W: Write> Write for BlockGzWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.buf.extend_from_slice(data);
        // Emit full blocks, keeping any remainder buffered.
        while self.buf.len() >= BLOCK_SIZE {
            let rest = self.buf.split_off(BLOCK_SIZE);
            let block = std::mem::replace(&mut self.buf, rest);
            self.emit_member(&block)?;
        }
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // Do not force a block boundary on every flush() — that would create
        // many tiny gzip members. Blocks are emitted at BLOCK_SIZE or finish().
        self.inner.flush()
    }
}

impl<W: Write> Drop for BlockGzWriter<W> {
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        if let Err(e) = self.flush_buffered() {
            error!("BlockGzWriter: failed to flush final gzip block on drop: {e}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;

    /// Decompress a (possibly multi-member) gzip buffer back to raw bytes,
    /// failing if it is not valid gzip.
    fn gunzip(bytes: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        MultiGzDecoder::new(bytes)
            .read_to_end(&mut out)
            .expect("output must be valid (multi-member) gzip");
        out
    }

    /// Compress `input` through a BlockGzWriter, finalizing with `finish()`.
    fn compress_finish(input: &[u8]) -> Vec<u8> {
        let mut buf = Vec::new();
        {
            let mut w = BlockGzWriter::new(&mut buf);
            w.write_all(input).unwrap();
            w.finish().unwrap();
        }
        buf
    }

    fn pattern(n: usize) -> Vec<u8> {
        // Deterministic, mildly compressible byte pattern.
        (0..n).map(|i| (i % 251) as u8).collect()
    }

    #[test]
    fn empty_input_yields_valid_empty_gzip() {
        let buf = compress_finish(&[]);
        // Must be a non-empty, valid gzip stream (not a 0-byte file, which a
        // decoder rejects) that decompresses to nothing.
        assert!(!buf.is_empty(), "empty input must still emit a gzip member");
        assert!(buf.starts_with(&[0x1f, 0x8b]), "must be gzip");
        assert_eq!(gunzip(&buf), Vec::<u8>::new());
    }

    #[test]
    fn small_input_round_trips() {
        let input = b"@read1\nACGTACGT\n+\nIIIIIIII\n";
        assert_eq!(gunzip(&compress_finish(input)), input);
    }

    #[test]
    fn multi_block_input_round_trips() {
        // > 2 blocks exercises the block-splitting path.
        let input = pattern(BLOCK_SIZE * 2 + 1234);
        assert_eq!(gunzip(&compress_finish(&input)), input);
    }

    #[test]
    fn exact_block_boundaries_round_trip() {
        for n in [BLOCK_SIZE - 1, BLOCK_SIZE, BLOCK_SIZE + 1, BLOCK_SIZE * 2] {
            let input = pattern(n);
            assert_eq!(gunzip(&compress_finish(&input)), input, "failed at len {n}");
        }
    }

    #[test]
    fn many_small_writes_reassemble() {
        // Byte-at-a-time writes must reassemble identically (covers the
        // cross-write buffering and block boundaries landing mid-write).
        let input = pattern(BLOCK_SIZE + 777);
        let mut buf = Vec::new();
        {
            let mut w = BlockGzWriter::new(&mut buf);
            for byte in &input {
                w.write_all(std::slice::from_ref(byte)).unwrap();
            }
            w.finish().unwrap();
        }
        assert_eq!(gunzip(&buf), input);
    }

    #[test]
    fn drop_without_finish_still_finalizes() {
        // runner.rs relies on Drop to finalize the per-contig temp writers.
        let input = pattern(BLOCK_SIZE + 50);
        let mut buf = Vec::new();
        {
            let mut w = BlockGzWriter::new(&mut buf);
            w.write_all(&input).unwrap();
            // no finish() — Drop must flush the final block.
        }
        assert_eq!(gunzip(&buf), input);
    }

    #[test]
    fn concatenated_streams_are_one_valid_stream() {
        // The combine step byte-concatenates per-contig temp files; the result
        // must decode as the concatenation of the inputs.
        let a = compress_finish(b"first contig reads\n");
        let b = compress_finish(b"second contig reads\n");
        let mut joined = a.clone();
        joined.extend_from_slice(&b);
        assert_eq!(
            gunzip(&joined),
            b"first contig reads\nsecond contig reads\n"
        );
    }
}
