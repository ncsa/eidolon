use crate::{
    file_tools::file_io::{is_gzipped_file, read_gzip_lines, read_lines},
    rng::{NeatRng, NeatRngError},
    structs::{
        nucleotides::{
            Nucleotide,
            Nucleotide::{N, X},
        },
        sequence_block::{RegionType, SequenceMap},
    },
};
use std::io;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FastaStreamError {
    #[error("I/O error reading FASTA: {0}")]
    IoError(#[from] io::Error),
    /// The file has content but no `>` header line, or has data before the first one.
    ///
    /// Silently skipping pre-header content let a gzipped FASTQ be read as a reference:
    /// `>` is ASCII 62, an ordinary Phred+33 quality character (Q29), so the reader
    /// walked past the reads until a QUALITY line happened to start with `>` and then
    /// used that quality string as a contig name. The run produced contigs named
    /// `=D2DFGGGCGGG...` and completed without complaint.
    #[error(
        "{path}: not a FASTA — {detail}. The first non-blank line of a FASTA must be a \
         '>' header. (If this is a FASTQ, pass the reference FASTA instead.)"
    )]
    NotFasta { path: String, detail: String },
}

/// Streaming FASTA reader that yields `(contig_name, raw_sequence)` one contig at a time.
///
/// The sequence is returned as a raw `String` (the concatenated FASTA lines) so that
/// callers with access to an RNG can call [`resolve_iupac_bases`] to stochastically
/// resolve IUPAC ambiguity codes before converting to [`Nucleotide`]. Callers that do
/// not need IUPAC resolution can convert directly with
/// `raw.chars().map(Nucleotide::from).collect()`, which maps IUPAC codes to `N`.
///
/// Only one contig's sequence is held in memory at a time, making this suitable for
/// large reference genomes without the temp-file overhead of `read_fasta`.
pub struct FastaStream {
    lines: Box<dyn Iterator<Item = io::Result<String>> + Send>,
    pending_name: Option<String>,
}

impl FastaStream {
    pub fn open(path: &PathBuf) -> Result<Self, FastaStreamError> {
        let lines: Box<dyn Iterator<Item = io::Result<String>> + Send> = if is_gzipped_file(path)? {
            Box::new(read_gzip_lines(path)?)
        } else {
            Box::new(read_lines(path)?)
        };

        let mut stream = FastaStream {
            lines,
            pending_name: None,
        };

        // Advance to the first header line.
        //
        // Anything other than blank lines before it is a hard error. This used to skip
        // silently, which meant a file that is not a FASTA at all could still "open":
        // a gzipped FASTQ got walked past its reads until a QUALITY line beginning with
        // '>' (ASCII 62 = Q29 in Phred+33) was mistaken for a header, and the run
        // proceeded with quality strings as contig names.
        let display = path.display().to_string();
        let mut skipped_nonblank: Option<String> = None;
        loop {
            match stream.lines.next() {
                None => {
                    return match skipped_nonblank {
                        Some(first) => Err(FastaStreamError::NotFasta {
                            path: display,
                            detail: format!(
                                "no '>' header line anywhere in the file; first content \
                                 line was {first:?}"
                            ),
                        }),
                        // A genuinely empty (or blank-only) file yields no contigs, which
                        // downstream code already handles; that is not a format error.
                        None => Ok(stream),
                    };
                }
                Some(Err(e)) => return Err(e.into()),
                Some(Ok(line)) if line.starts_with('>') => {
                    if let Some(first) = skipped_nonblank {
                        return Err(FastaStreamError::NotFasta {
                            path: display,
                            detail: format!("content appears before the first header: {first:?}"),
                        });
                    }
                    stream.pending_name = Some(parse_contig_name(&line));
                    break;
                }
                Some(Ok(line)) => {
                    if !line.trim().is_empty() && skipped_nonblank.is_none() {
                        skipped_nonblank = Some(line.chars().take(60).collect::<String>());
                    }
                }
            }
        }

        Ok(stream)
    }
}

impl Iterator for FastaStream {
    type Item = Result<(String, String), FastaStreamError>;

    fn next(&mut self) -> Option<Self::Item> {
        let name = self.pending_name.take()?;
        let mut sequence = String::new();

        loop {
            match self.lines.next() {
                None => break,
                Some(Err(e)) => return Some(Err(e.into())),
                Some(Ok(line)) => {
                    if line.starts_with('>') {
                        self.pending_name = Some(parse_contig_name(&line));
                        break;
                    }
                    sequence.push_str(&line);
                }
            }
        }

        Some(Ok((name, sequence)))
    }
}

fn parse_contig_name(header: &str) -> String {
    header[1..]
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

/// Builds a `Vec<SequenceMap>` that partitions `sequence` into contiguous N/non-N runs.
///
/// The returned map covers the full sequence without gaps or overlaps.
/// Never returns an error; the `Result` wrapper is kept for call-site compatibility.
pub fn map_buffer(sequence: &[Nucleotide]) -> Vec<SequenceMap> {
    if sequence.is_empty() {
        return Vec::new();
    }
    let mut map: Vec<SequenceMap> = Vec::new();
    let mut region_start = 0;
    let mut region_end = 1;
    // Soft-masked bases (Maskeda/c/g/t) represent repeat-annotated but valid sequence;
    // treat them as regular bases for region mapping. Only N and X mark true gaps.
    let mut inside_n_region = matches!(sequence[0], N | X);

    if inside_n_region {
        for base in &sequence[1..] {
            match base {
                N | X => region_end += 1,
                _ => {
                    map.push(SequenceMap::from(RegionType::NRegion, 0, region_end));
                    region_start = region_end;
                    region_end = region_start + 1;
                    inside_n_region = false;
                    break;
                }
            }
        }
        if inside_n_region {
            map.push(SequenceMap::from(RegionType::NRegion, 0, region_end));
            return map;
        }
    }

    for base in &sequence[region_end..] {
        match base {
            N | X => {
                if inside_n_region {
                    region_end += 1;
                } else {
                    inside_n_region = true;
                    map.push(SequenceMap::from(
                        RegionType::NonNRegion,
                        region_start,
                        region_end,
                    ));
                    region_start = region_end;
                    region_end = region_start + 1;
                }
            }
            _ => {
                if inside_n_region {
                    inside_n_region = false;
                    map.push(SequenceMap::from(
                        RegionType::NRegion,
                        region_start,
                        region_end,
                    ));
                    region_start = region_end;
                    region_end = region_start + 1;
                } else {
                    region_end += 1;
                }
            }
        }
    }

    let region_type = if inside_n_region {
        RegionType::NRegion
    } else {
        RegionType::NonNRegion
    };
    map.push(SequenceMap::from(region_type, region_start, region_end));
    map
}

/// Converts a raw FASTA sequence string to nucleotides, stochastically resolving
/// IUPAC ambiguity codes (R/Y/M/K/S/W/H/B/V/D) to their constituent bases via `rng`.
/// N is preserved as-is (assembly gap, not ambiguity). Soft-masked bases (a/c/g/t) pass
/// through unchanged. Returns the resolved sequence and the count of IUPAC bases resolved.
pub fn resolve_iupac_bases(
    raw: &str,
    rng: &mut NeatRng,
) -> Result<(Vec<Nucleotide>, usize), NeatRngError> {
    let mut sequence = Vec::with_capacity(raw.len());
    let mut iupac_count = 0usize;

    for c in raw.chars() {
        let resolved = match c.to_ascii_uppercase() {
            'R' => {
                iupac_count += 1;
                iupac_pick(&['G', 'A'], rng)?
            }
            'Y' => {
                iupac_count += 1;
                iupac_pick(&['C', 'T'], rng)?
            }
            'M' => {
                iupac_count += 1;
                iupac_pick(&['A', 'C'], rng)?
            }
            'K' => {
                iupac_count += 1;
                iupac_pick(&['G', 'T'], rng)?
            }
            'S' => {
                iupac_count += 1;
                iupac_pick(&['C', 'G'], rng)?
            }
            'W' => {
                iupac_count += 1;
                iupac_pick(&['A', 'T'], rng)?
            }
            'H' => {
                iupac_count += 1;
                iupac_pick(&['A', 'C', 'T'], rng)?
            }
            'B' => {
                iupac_count += 1;
                iupac_pick(&['C', 'G', 'T'], rng)?
            }
            'V' => {
                iupac_count += 1;
                iupac_pick(&['A', 'C', 'G'], rng)?
            }
            'D' => {
                iupac_count += 1;
                iupac_pick(&['A', 'G', 'T'], rng)?
            }
            _ => c,
        };
        sequence.push(Nucleotide::from(resolved));
    }

    Ok((sequence, iupac_count))
}

fn iupac_pick(bases: &[char], rng: &mut NeatRng) -> Result<char, NeatRngError> {
    let rand = rng.random()?;
    let idx = (rand * bases.len() as f64) as usize;
    Ok(bases[idx.min(bases.len() - 1)])
}

/// Returns the contiguous non-N/non-X regions of a sequence as `(start, end)` pairs.
pub fn non_n_regions(sequence: &[Nucleotide]) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut region_start: Option<usize> = None;

    for (i, &nuc) in sequence.iter().enumerate() {
        match (nuc, region_start) {
            (Nucleotide::N | Nucleotide::X, Some(s)) => {
                regions.push((s, i));
                region_start = None;
            }
            (Nucleotide::N | Nucleotide::X, None) => {}
            (_, None) => region_start = Some(i),
            (_, Some(_)) => {}
        }
    }
    if let Some(s) = region_start {
        regions.push((s, sequence.len()));
    }
    regions
}

/// Scans a FASTA file and returns `(contig_name, length)` for each contig
/// without storing sequences. Used to build BAM headers before streaming reads.
pub fn scan_fasta_lengths(path: &PathBuf) -> Result<Vec<(String, usize)>, FastaStreamError> {
    let lines: Box<dyn Iterator<Item = io::Result<String>>> = if is_gzipped_file(path)? {
        Box::new(read_gzip_lines(path)?)
    } else {
        Box::new(read_lines(path)?)
    };

    let mut result: Vec<(String, usize)> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_len: usize = 0;

    for line in lines {
        let line = line?;
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                result.push((name, current_len));
            }
            current_name = Some(parse_contig_name(&line));
            current_len = 0;
        } else {
            current_len += line.chars().filter(|c| !c.is_whitespace()).count();
        }
    }
    if let Some(name) = current_name {
        result.push((name, current_len));
    }

    Ok(result)
}

#[cfg(test)]
mod tests {

    /// A FASTQ must be rejected, not silently parsed as a reference.
    ///
    /// This is the real defect, reproduced: `>` is ASCII 62, an ordinary Phred+33
    /// quality character (Q29), so a quality line CAN begin with it. The old reader
    /// skipped everything before the first `>` without comment, walked past the reads,
    /// and adopted a quality string as a contig name. A whole gen-reads run then
    /// completed with contigs called `=D2DFGGGCGGGBFBGG...` and nothing flagged it.
    #[test]
    fn a_fastq_is_rejected_rather_than_parsed_as_a_reference() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.fq");
        // Quality line deliberately begins with '>' — the character that made this work.
        std::fs::write(
            &path,
            "@read1\nACGTACGTAC\n+\n>GGFGGFGGFG\n@read2\nTTTTGGGGCC\n+\nIIIIIIIIII\n",
        )
        .unwrap();
        let err = match FastaStream::open(&path) {
            Ok(_) => panic!("a FASTQ must not open as a FASTA"),
            Err(e) => e,
        };
        let msg = err.to_string();
        assert!(
            msg.contains("not a FASTA"),
            "error should say the file is not a FASTA; got: {msg}"
        );
        assert!(
            msg.contains("@read1"),
            "error should quote the offending first content line; got: {msg}"
        );
    }

    /// No header anywhere is also a hard error, not an empty contig set.
    #[test]
    fn a_headerless_file_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("plain.txt");
        std::fs::write(&path, "ACGTACGT\nACGTACGT\n").unwrap();
        let err = match FastaStream::open(&path) {
            Ok(_) => panic!("a headerless file must not open"),
            Err(e) => e,
        };
        assert!(err.to_string().contains("no '>' header"), "{err}");
    }

    /// MUST NOT FIRE: a well-formed FASTA still opens, including one that starts with
    /// blank lines. Rejecting those would break real files.
    #[test]
    fn a_wellformed_fasta_still_opens_even_with_leading_blank_lines() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ref.fa");
        std::fs::write(&path, "\n\n>chr1 desc here\nACGT\nACGT\n>chr2\nTTTT\n").unwrap();
        let got: Vec<(String, String)> = FastaStream::open(&path)
            .expect("a valid FASTA must open")
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(
            got,
            vec![
                ("chr1".to_string(), "ACGTACGT".to_string()),
                ("chr2".to_string(), "TTTT".to_string()),
            ]
        );
    }
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{}", content).unwrap();
        f
    }

    fn write_gzipped_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        let mut enc = GzEncoder::new(&mut f, Compression::default());
        enc.write_all(content.as_bytes()).unwrap();
        enc.finish().unwrap();
        f
    }

    fn make_sequence(bases: &str) -> Vec<Nucleotide> {
        bases.chars().map(Nucleotide::from).collect()
    }

    // FastaStream tests

    #[test]
    fn test_reads_single_contig() {
        let f = write_fasta(">chr1\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[0].1, "ACGT");
    }

    #[test]
    fn test_reads_multiple_contigs() {
        let f = write_fasta(">chr1\nACGT\n>chr2\nTTTT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[1].0, "chr2");
    }

    #[test]
    fn test_parses_multiline_sequence() {
        let f = write_fasta(">chr1\nACGT\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].1.len(), 8);
    }

    #[test]
    fn test_strips_description_from_header() {
        let f = write_fasta(">chr1 some description here\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].0, "chr1");
    }

    #[test]
    fn test_reads_gzipped_fasta() {
        let f = write_gzipped_fasta(">chr1\nACGT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].0, "chr1");
        assert_eq!(contigs[0].1.len(), 4);
    }

    #[test]
    fn test_n_bases_preserved_as_raw_char() {
        let f = write_fasta(">chr1\nACNT\n");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(contigs[0].1.chars().nth(2).unwrap(), 'N');
    }

    #[test]
    fn test_empty_file_yields_no_contigs() {
        let f = write_fasta("");
        let contigs: Vec<_> = FastaStream::open(&f.path().to_path_buf())
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(contigs.is_empty());
    }

    // non_n_regions tests

    #[test]
    fn test_non_n_regions_no_ns() {
        let seq = make_sequence("ACGTACGT");
        assert_eq!(non_n_regions(&seq), vec![(0, 8)]);
    }

    #[test]
    fn test_non_n_regions_leading_n() {
        let seq = make_sequence("NNACGT");
        assert_eq!(non_n_regions(&seq), vec![(2, 6)]);
    }

    #[test]
    fn test_non_n_regions_interior_n() {
        let seq = make_sequence("ACGTNACGT");
        assert_eq!(non_n_regions(&seq), vec![(0, 4), (5, 9)]);
    }

    #[test]
    fn test_non_n_regions_all_n() {
        let seq = make_sequence("NNNN");
        assert!(non_n_regions(&seq).is_empty());
    }

    #[test]
    fn test_non_n_regions_x_closes_region() {
        let mut seq = make_sequence("ACGT");
        seq.push(Nucleotide::X);
        seq.extend(make_sequence("ACGT"));
        assert_eq!(non_n_regions(&seq), vec![(0, 4), (5, 9)]);
    }

    #[test]
    fn test_non_n_regions_x_at_start_is_skipped() {
        let mut seq = vec![Nucleotide::X, Nucleotide::X];
        seq.extend(make_sequence("ACGT"));
        assert_eq!(non_n_regions(&seq), vec![(2, 6)]);
    }

    // scan_fasta_lengths tests

    #[test]
    fn test_scan_fasta_lengths_single_contig() {
        let f = write_fasta(">chr1\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths, vec![("chr1".to_string(), 4)]);
    }

    #[test]
    fn test_scan_fasta_lengths_multiple_contigs() {
        let f = write_fasta(">chr1\nACGT\n>chr2\nGGGGGGGGGG\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(
            lengths,
            vec![("chr1".to_string(), 4), ("chr2".to_string(), 10),]
        );
    }

    #[test]
    fn test_scan_fasta_lengths_multiline() {
        let f = write_fasta(">chr1\nACGT\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths, vec![("chr1".to_string(), 8)]);
    }

    #[test]
    fn test_scan_fasta_lengths_strips_description() {
        let f = write_fasta(">chr1 some description\nACGT\n");
        let lengths = scan_fasta_lengths(&f.path().to_path_buf()).unwrap();
        assert_eq!(lengths[0].0, "chr1");
    }

    // resolve_iupac_bases tests

    #[test]
    fn test_resolve_iupac_all_codes_yield_only_acgtn() {
        let mut rng = NeatRng::new_from_seed(&vec!["iupac_test".to_string()]).unwrap();
        let (seq, count) = resolve_iupac_bases("RYMKSWBVDHACGTN", &mut rng).unwrap();
        assert_eq!(count, 10);
        assert_eq!(seq.len(), 15);
        for nuc in &seq {
            assert!(
                matches!(
                    nuc,
                    Nucleotide::A | Nucleotide::C | Nucleotide::G | Nucleotide::T | Nucleotide::N
                ),
                "unexpected nucleotide {:?}",
                nuc
            );
        }
    }

    /// Every IUPAC code must resolve to EXACTLY its own constituents.
    ///
    /// `test_resolve_iupac_all_codes_yield_only_acgtn` asserts membership in
    /// {A,C,G,T,N}, which ANY wrong mapping satisfies — resolving Y (pyrimidine, C/T)
    /// to purines passed it and every other test in the workspace. Only R and H had
    /// distribution tests, leaving Y, M, K, S, W, B, V and D unconstrained.
    #[test]
    fn resolve_iupac_yields_exactly_each_code_constituents() {
        // IUPAC nucleotide codes, per the IUBMB/IUPAC-IUB definition.
        let codes: [(char, &[Nucleotide]); 10] = [
            ('R', &[Nucleotide::A, Nucleotide::G]), // puRine
            ('Y', &[Nucleotide::C, Nucleotide::T]), // pYrimidine
            ('M', &[Nucleotide::A, Nucleotide::C]), // aMino
            ('K', &[Nucleotide::G, Nucleotide::T]), // Keto
            ('S', &[Nucleotide::C, Nucleotide::G]), // Strong (3 H-bonds)
            ('W', &[Nucleotide::A, Nucleotide::T]), // Weak (2 H-bonds)
            ('B', &[Nucleotide::C, Nucleotide::G, Nucleotide::T]), // not A
            ('D', &[Nucleotide::A, Nucleotide::G, Nucleotide::T]), // not C
            ('H', &[Nucleotide::A, Nucleotide::C, Nucleotide::T]), // not G
            ('V', &[Nucleotide::A, Nucleotide::C, Nucleotide::G]), // not T
        ];
        let mut rng = NeatRng::new_from_seed(&vec!["iupac constituents".to_string()]).unwrap();
        for (code, expected) in codes {
            let input: String = std::iter::repeat(code).take(600).collect();
            let (seq, count) = resolve_iupac_bases(&input, &mut rng).unwrap();
            assert_eq!(count, 600, "{code}: wrong ambiguity count");

            let mut seen: Vec<Nucleotide> = Vec::new();
            for n in &seq {
                if !seen.contains(n) {
                    seen.push(*n);
                }
            }
            // Nothing outside the code's constituents...
            for n in &seen {
                assert!(
                    expected.contains(n),
                    "{code} resolved to {n:?}, which is not one of its constituents \
                     {expected:?}"
                );
            }
            // ...and every constituent actually appears, so a mapping that collapses a
            // code onto a strict subset (e.g. always picking the first base) fails too.
            for n in expected {
                assert!(
                    seen.contains(n),
                    "{code} never resolved to {n:?} in 600 draws; observed only {seen:?}"
                );
            }
        }
    }

    #[test]
    fn test_resolve_iupac_preserves_acgtn_unchanged() {
        let mut rng = NeatRng::new_from_seed(&vec!["iupac_test".to_string()]).unwrap();
        let (seq, count) = resolve_iupac_bases("ACGTN", &mut rng).unwrap();
        assert_eq!(count, 0);
        assert_eq!(
            seq,
            vec![
                Nucleotide::A,
                Nucleotide::C,
                Nucleotide::G,
                Nucleotide::T,
                Nucleotide::N
            ]
        );
    }

    #[test]
    fn test_resolve_iupac_r_uniform_distribution() {
        // R resolves to G or A; over many draws expect roughly 50/50
        let mut rng = NeatRng::new_from_seed(&vec!["iupac_stat_test".to_string()]).unwrap();
        let raw: String = "R".repeat(1000);
        let (seq, count) = resolve_iupac_bases(&raw, &mut rng).unwrap();
        assert_eq!(count, 1000);
        let g_count = seq.iter().filter(|&&n| n == Nucleotide::G).count();
        let a_count = seq.iter().filter(|&&n| n == Nucleotide::A).count();
        assert_eq!(g_count + a_count, 1000);
        assert!(
            (400..=600).contains(&g_count),
            "G count {} outside [400, 600]",
            g_count
        );
    }

    #[test]
    fn test_resolve_iupac_h_uniform_distribution() {
        // H resolves to A, C, or T; over many draws expect roughly equal thirds
        let mut rng = NeatRng::new_from_seed(&vec!["iupac_h_test".to_string()]).unwrap();
        let raw: String = "H".repeat(3000);
        let (seq, count) = resolve_iupac_bases(&raw, &mut rng).unwrap();
        assert_eq!(count, 3000);
        let a_count = seq.iter().filter(|&&n| n == Nucleotide::A).count();
        let c_count = seq.iter().filter(|&&n| n == Nucleotide::C).count();
        let t_count = seq.iter().filter(|&&n| n == Nucleotide::T).count();
        assert_eq!(a_count + c_count + t_count, 3000);
        // Each constituent should be ~1000; allow ±20% (200)
        assert!(
            (800..=1200).contains(&a_count),
            "A count {} outside [800, 1200]",
            a_count
        );
        assert!(
            (800..=1200).contains(&c_count),
            "C count {} outside [800, 1200]",
            c_count
        );
        assert!(
            (800..=1200).contains(&t_count),
            "T count {} outside [800, 1200]",
            t_count
        );
    }

    #[test]
    fn test_resolve_iupac_masked_bases_pass_through() {
        let mut rng = NeatRng::new_from_seed(&vec!["iupac_mask_test".to_string()]).unwrap();
        let (seq, count) = resolve_iupac_bases("acgt", &mut rng).unwrap();
        assert_eq!(count, 0);
        assert_eq!(
            seq,
            vec![
                Nucleotide::Maskeda,
                Nucleotide::Maskedc,
                Nucleotide::Maskedg,
                Nucleotide::Maskedt,
            ]
        );
    }
}
