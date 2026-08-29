use eidolon_core::file_tools::block_gz::BlockGzWriter;
use eidolon_core::file_tools::file_io::create_output_file;
use eidolon_core::rng::NeatRng;
use eidolon_core::structs::variants::{Genotype, Provenance, SvType, VariantType};
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::Arc;
use tempfile;

use crate::{
    eidolon_core::{
        file_tools::{
            bam_writer::{BamBodyWriter, BamContext, BamRecordStager, concat_temp_bams},
            bed_reader::read_bed,
            fasta_stream::{FastaStream, map_buffer, resolve_iupac_bases},
            fastq_tools::{
                HaplotypeContext, PlacedFragment, Strand, combine_temp_fastqs, generate_read,
                reverse_complement, write_block_fastq, write_read_to_fastq,
            },
            file_io::{VectorBuffer, append_to_file},
            vcf_tools::{read_vcf, write_vcf},
        },
        models::{
            fragment_length::FragmentLengthModel, mutation_model::MutationModel,
            quality_scores::QualityScoreModel, sequencing_error_model::SequencingErrorModel,
        },
        structs::{
            bed_record::BedRecord,
            haplotype_map::InsertionCoordinateMap,
            mutated_map::{AdCounter, MutatedMap},
            nucleotides::Nucleotide,
            read_record::ReadRecord,
            sequence_block::{RegionType, SequenceBlock, SequenceMap},
            sv_model::MIN_SV_LENGTH_BP,
            variants::{SvData, Variant},
        },
    },
    gen_reads::{
        errors::GenerateReadsError,
        utils::{
            config::RunConfiguration,
            generate_fragments::{generate_fragments, generate_weighted_fragments},
            generate_variants::generate_variants,
            subclone::SubcloneModel,
        },
    },
};
use eidolon_core::models::gc_bias_model::GcBiasModel;
use flate2::{Compression, write::GzEncoder};

/// Dedicated RNG sub-stream for the genome-wide translocation pass, so adding it leaves
/// every per-contig sampling decision bit-for-bit unchanged.
const TRANSLOCATION_STREAM: u64 = 9_100_000;
const TRANSLOCATION_CCF_STREAM: u64 = TRANSLOCATION_STREAM + 1;

struct ContigContext<'a> {
    config: &'a RunConfiguration,
    target_bed: &'a Option<HashMap<String, Vec<BedRecord>>>,
    mutation_regions: &'a Option<HashMap<String, Vec<BedRecord>>>,
    default_run_mutation_rate: f64,
    fragment_length_model: &'a FragmentLengthModel,
    gc_bias_model: &'a GcBiasModel,
    quality_score_model: &'a QualityScoreModel,
    seq_error_model: &'a SequencingErrorModel,
    working_dir: &'a std::path::Path,
    base_rng: NeatRng,
    bam_context: Option<Arc<BamContext>>,
    reference: Arc<HashMap<String, Vec<Nucleotide>>>,
    mutated_maps: Arc<HashMap<String, MutatedMap>>,
    // Contigs in REFERENCE FILE order. `mutated_maps` and `reference` are HashMaps, whose
    // iteration order Rust randomizes per process; anything that walks them while drawing
    // from a shared RNG becomes nondeterministic run to run. That is #599: chimeric read
    // generation iterated `mutated_maps` directly and produced a different read set every
    // run at a fixed seed, single-threaded. Iterate this instead and look the map up.
    contig_order: Arc<Vec<String>>,
    max_del_lens: HashMap<String, usize>,
}

pub fn run_neat(
    config: &RunConfiguration,
    rng: &mut NeatRng,
) -> Result<Vec<PathBuf>, GenerateReadsError> {
    let working_dir = tempfile::tempdir().unwrap();
    info!("Created temp dir at {:?}", working_dir);

    info!("Generate mutation model");
    let mutation_model = {
        match &config.mutation_model {
            Some(filename) => MutationModel::from_file(filename)?,
            None => MutationModel::default()?,
        }
    };
    let mutation_regions = match &config.mutation_regions {
        Some(path) => {
            info!("Loading mutation regions BED: {:?}", path);
            Some(read_bed(path, true)?)
        }
        None => None,
    };
    let default_run_mutation_rate = match config.mutation_rate {
        Some(rate) => rate,
        None => mutation_model.mutation_rate,
    };

    info!("Generate fragment length model");
    let fragment_length_model: FragmentLengthModel = {
        match &config.fragment_model {
            Some(filename) => FragmentLengthModel::discrete_from_file(filename)?,
            None => match config.fragment_mean {
                Some(mean) => {
                    FragmentLengthModel::new_normal(mean, config.fragment_st_dev.unwrap())?
                }
                None => FragmentLengthModel::default()?,
            },
        }
    };

    info!("Generate sequencing error model");
    let seq_error_model: SequencingErrorModel = {
        match &config.sequence_error_model {
            Some(filename) => SequencingErrorModel::from_file(filename)?,
            None => SequencingErrorModel::default()?,
        }
    };

    info!("Generate quality score model");
    let quality_score_model: QualityScoreModel = {
        match &config.quality_score_model {
            Some(filename) => QualityScoreModel::from_file(filename)?,
            // Fall through to the quality model embedded in the sequencing-error model
            // when the user hasn't supplied an explicit override. This matches the
            // documented contract of `sequence_error_model:` and ensures binned-quality
            // training (gen-seq-error-model with binned_quality_bins) actually drives
            // gen-reads sampling — otherwise the binned QSM sits unused inside the
            // SeqErrorModel while the default continuous model is used here.
            None => seq_error_model.quality_score_model().clone(),
        }
    };

    let gc_bias_model = match &config.gc_bias_model {
        Some(path) => {
            info!("Loading GC Bias model: {}", path.display());
            GcBiasModel::from_file(path)?
        }
        None => GcBiasModel::default(),
    };

    let target_bed = match &config.target_bed {
        Some(path) => {
            info!("Loading target BED: {:?}", path);
            Some(read_bed(path, false)?)
        }
        None => None,
    };

    let input_variants: Option<HashMap<String, Vec<Variant>>> = match &config.input_vcf {
        Some(path) => {
            info!("Loading input VCF: {}", path.display());
            let raw = read_vcf(path.to_path_buf())?;
            Some(filter_input_vcf(raw))
        }
        None => None,
    };

    // #405 reproductive somatic: replay a supplied somatic VCF in this pass. Tag its
    // variants `SomaticVcf` (so the tumor/normal merge resolves origin `somatic`, not
    // `shared`) and scale their allele_fraction by `somatic_af_scale` — the cancer
    // tumor pass sets 1/purity, so a supplied *observed* VAF reproduces after mixing.
    let input_variants = if let Some(path) = &config.somatic_vcf {
        info!("Loading reproductive somatic VCF: {}", path.display());
        let mut som = filter_input_vcf(read_vcf(path.to_path_buf())?);
        let mut clamped = 0usize;
        for variants in som.values_mut() {
            for v in variants.iter_mut() {
                v.provenance = Provenance::SomaticVcf;
                // Drop the source VCF's INFO (like germline input_vcf) so it doesn't
                // leak into the golden; re-populate only with our own ground-truth tag.
                v.info = None;
                if let Some(af) = v.allele_fraction {
                    let raw = af * config.somatic_af_scale;
                    let scaled = if raw > 1.0 {
                        clamped += 1;
                        1.0
                    } else {
                        raw
                    };
                    v.allele_fraction = Some(scaled);
                    // EIDOLON_VAF = intended observed VAF after mixing = purity × scaled =
                    // the input observed VAF (or purity, if it was clamped).
                    if let Some(p) = config.merged_vaf_purity {
                        append_info_tag(&mut v.info, format!("EIDOLON_VAF={:.4}", p * scaled));
                    }
                }
            }
        }
        if clamped > 0 {
            warn!(
                "somatic_vcf: clamped {clamped} scaled allele fraction(s) > 1.0 \
                 (observed VAF exceeded purity)"
            );
        }
        // Merge the somatic variants into the germline input map (or stand alone).
        Some(match input_variants {
            Some(mut germline) => {
                for (contig, mut vs) in som {
                    germline.entry(contig).or_default().append(&mut vs);
                }
                germline
            }
            None => som,
        })
    } else {
        input_variants
    };

    info!("Reading fasta file: {}", config.reference.display());
    let mut reference_map = HashMap::new();
    let mut contig_order_in_file = Vec::new();
    let fasta = FastaStream::open(&config.reference)?;
    for (idx, result) in fasta.enumerate() {
        let (name, raw) = result?;
        // Always record the contig in file order so downstream per-contig RNG
        // derivations (IUPAC at `idx`, mutated maps at `idx`, chunks at
        // contig_idx) key off the contig's ORIGINAL position — preserving
        // byte-identical output regardless of which contigs are loaded.
        contig_order_in_file.push(name.clone());
        // target_bed-aware loading: only materialize the sequence for contigs in
        // the target BED. target_bed already restricts read GENERATION to those
        // contigs (process_chunk skips the rest), so non-target contigs were
        // loaded into RAM but never used — wasteful at genome scale (a chr1
        // shard otherwise holds all ~3 GB of GRCh38 to simulate 25 Mb). Skipping
        // them here changes no output (those contigs produce no reads either
        // way); it only drops the dead weight.
        if let Some(bed) = &target_bed
            && !bed.contains_key(&name)
        {
            continue;
        }
        let mut child_rng = rng.derive_child(idx as u64);
        let (seq, iupac_count) = resolve_iupac_bases(&raw, &mut child_rng)?;
        if iupac_count > 0 {
            warn!(
                "Contig {}: resolved {} IUPAC ambiguity base(s) to ACGT",
                name, iupac_count
            );
        }
        // N bases are left as-is. The non-N region machinery (map_buffer /
        // get_non_n_regions) excludes them from read anchoring, mutation
        // placement, and SV anchoring, so assembly gaps stay as coverage
        // dropouts rather than being filled with fabricated sequence.
        reference_map.insert(name, seq);
    }
    let reference = Arc::new(reference_map);

    let bam_context: Option<Arc<BamContext>> = if config.produce_bam {
        // Only loaded (target) contigs have sequence; with target_bed set the
        // BAM header covers the targeted contigs (the only ones that get reads).
        let contig_lengths: Vec<(String, usize)> = contig_order_in_file
            .iter()
            .filter_map(|name| reference.get(name).map(|s| (name.clone(), s.len())))
            .collect();
        Some(Arc::new(BamContext::new(&contig_lengths)))
    } else {
        None
    };

    // Phase 0: place inter-chromosomal translocations across the WHOLE genome.
    //
    // A breakend's two ends live on different contigs and each generates reads from its
    // own contig, so both records must exist before either contig is processed. The
    // per-contig sampler cannot do this and used to hardcode the mate to the anchor's own
    // contig, which made every "translocation" same-contig (466 of 466, job 20719077).
    // BND's share of the per-base rate is subtracted from the per-contig budget inside
    // `sample_variants`, so total SV yield is unchanged — it is split, not duplicated.
    let mut translocations: HashMap<String, Vec<Variant>> = if config.sv_rate_scale > 0.0
        && let Some(sv_model) = mutation_model.sv_model.as_ref()
        && sv_model.is_usable()
    {
        let mut tra_rng = rng.derive_child(TRANSLOCATION_STREAM);
        let t = sv_model.sample_translocations(
            &contig_order_in_file,
            reference.as_ref(),
            config.ploidy,
            config.sv_rate_scale,
            &mut tra_rng,
        );
        let n: usize = t.values().map(|v| v.len()).sum();
        if n > 0 {
            info!(
                "Placed {} inter-chromosomal translocation(s) ({} breakend records) across {} contig(s)",
                n / 2,
                n,
                t.len()
            );
        }
        t
    } else {
        HashMap::new()
    };
    if config.subclone_model.is_some() && !translocations.is_empty() {
        let mut ccf_rng = rng.derive_child(TRANSLOCATION_CCF_STREAM);
        apply_translocation_subclone_model(
            &mut translocations,
            config.subclone_model.as_ref(),
            config.ploidy,
            config.merged_vaf_purity,
            &mut ccf_rng,
        )?;
    }

    // Phase 1: Generate MutatedMaps for all contigs
    info!("Generating mutations for all contigs");
    let mut all_mutated_maps = HashMap::new();
    let mut max_del_lens = HashMap::new();
    for (m_idx, name) in contig_order_in_file.iter().enumerate() {
        // Skip contigs whose sequence wasn't loaded (non-target under
        // target_bed). m_idx stays the original file position, so the loaded
        // contigs' mutation RNG (derive_child(m_idx + 1000000)) is unchanged.
        let Some(seq) = reference.get(name) else {
            continue;
        };
        let m_rng = rng.derive_child((m_idx + 1000000) as u64);
        let (m_map, max_del) = generate_mutated_map(
            name,
            seq,
            config,
            &target_bed,
            &mutation_regions,
            &input_variants,
            translocations.get(name),
            &mutation_model,
            default_run_mutation_rate,
            m_rng,
        )?;
        all_mutated_maps.insert(name.clone(), m_map);
        max_del_lens.insert(name.clone(), max_del);
    }
    let shared_mutated_maps = Arc::new(all_mutated_maps);

    let ctx = ContigContext {
        config,
        target_bed: &target_bed,
        mutation_regions: &mutation_regions,
        default_run_mutation_rate,
        fragment_length_model: &fragment_length_model,
        gc_bias_model: &gc_bias_model,
        quality_score_model: &quality_score_model,
        seq_error_model: &seq_error_model,
        working_dir: working_dir.path(),
        base_rng: *rng,
        bam_context,
        reference,
        mutated_maps: shared_mutated_maps,
        contig_order: Arc::new(contig_order_in_file.clone()),
        max_del_lens,
    };

    let mut mutated_maps: HashMap<String, Vec<MutatedMap>> = HashMap::new();
    let mut all_fastq_files: HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)> = HashMap::new();
    let mut contig_order: Vec<String> = Vec::new();
    let mut fasta_lengths: HashMap<String, usize> = HashMap::new();
    // Multiple BAM body files per contig (one per chunk); value is (chunk_start, path).
    let mut bam_body_files: HashMap<String, Vec<(usize, PathBuf)>> = HashMap::new();
    let mut ad_counters: HashMap<String, AdCounter> = HashMap::new();

    info!("Generating simulated dataset");

    // Flatten contigs into fixed-size sub-contig chunks so work parallelizes
    // *within* a large chromosome, not just across contigs. Chunk size is
    // independent of num_threads, so output is identical regardless of thread
    // count. Each contig's sequence is shared across its chunks via the cache.
    let chunk_size = resolve_chunk_size(config);
    match config.chunk_size {
        None | Some(0) => info!("\t>Sub-contig chunking: disabled (one chunk per contig)"),
        Some(n) => info!("\t>Sub-contig chunking: enabled, chunk size {} bp", n),
    }
    let seq_cache = SeqBlockCache::new(&ctx.reference, &contig_order_in_file);
    let chunk_work: Vec<ChunkWork> = contig_order_in_file
        .iter()
        .enumerate()
        .flat_map(|(contig_idx, name)| {
            // Non-target contigs aren't loaded (target_bed-aware loading), so
            // they produce no chunks — skip them rather than emitting a chunk
            // that would panic in process_chunk's seq_cache.get. contig_idx
            // (the enumerate position) is preserved for loaded contigs, so their
            // per-chunk RNG (derive_child(contig_idx)) is unchanged.
            let contig_len = match ctx.reference.get(name) {
                Some(seq) => seq.len(),
                None => 0,
            };
            let chunks = if contig_len == 0 {
                Vec::new()
            } else {
                split_contig_into_chunks(contig_len, chunk_size)
            };
            chunks
                .into_iter()
                .enumerate()
                .map(move |(chunk_idx, (chunk_start, chunk_end))| ChunkWork {
                    contig_idx,
                    chunk_idx,
                    name: name.clone(),
                    chunk_start,
                    chunk_end,
                })
        })
        .collect();

    let parallel_iter =
        chunk_work
            .into_par_iter()
            .map(|w| -> Result<ChunkResult, GenerateReadsError> {
                // Deterministic per-chunk seed. Chunk 0 (the only chunk when
                // chunking is disabled) uses the plain per-contig derivation so its
                // RNG stream — and therefore its reads — matches the pre-chunking
                // behaviour exactly. Later chunks derive a sub-seed from it.
                let child_rng = if w.chunk_idx == 0 {
                    ctx.base_rng.derive_child(w.contig_idx as u64)
                } else {
                    ctx.base_rng
                        .derive_child(w.contig_idx as u64)
                        .derive_child(w.chunk_idx as u64)
                };
                let block = seq_cache.get(&w.name);
                process_chunk(
                    w.contig_idx,
                    w.chunk_idx,
                    w.name,
                    w.chunk_start,
                    w.chunk_end,
                    block,
                    &ctx,
                    child_rng,
                )
            });
    let collected: Result<Vec<ChunkResult>, _> = match config.num_threads {
        Some(n) => rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .map_err(|e| {
                GenerateReadsError::CliError(format!("Failed to build thread pool: {}", e))
            })?
            .install(|| parallel_iter.collect()),
        None => parallel_iter.collect(),
    };
    let mut results = collected?;
    // Order by (contig, chunk_start) so BAM bodies concatenate coordinate-sorted.
    results.sort_unstable_by_key(|r| (r.contig_idx, r.chunk_start));

    // Phase 3: Generate chimeric reads for BND junctions
    if config.produce_fastq || config.produce_bam {
        let chimeric_rng = ctx.base_rng.derive_child(999999);
        let chimeric_res = process_chimeric_variants(&ctx, chimeric_rng)?;
        results.push(chimeric_res);
    }

    for cr in results {
        collect_chunk_result(
            cr,
            &mut contig_order,
            &mut fasta_lengths,
            &mut all_fastq_files,
            &mut bam_body_files,
            &mut ad_counters,
        );
    }

    // The golden VCF reads per-contig MutatedMaps from the shared Arc (chunks no
    // longer carry their own clone).
    for name in &contig_order {
        if let Some(m) = ctx.mutated_maps.get(name) {
            mutated_maps.insert(name.clone(), vec![m.clone()]);
        }
    }

    info!("Read generation complete, producing output files");

    if config.produce_fastq {
        info!("Producing final fastq(s) file(s)");

        // Concatenate in contig order (each contig's chunk files are already in
        // chunk_start order), then the chimeric reads, so the FASTQ is
        // byte-identical regardless of thread count — not just the same read set.
        let mut all_r1: Vec<PathBuf> = Vec::new();
        let mut all_r2: Vec<PathBuf> = Vec::new();
        for name in &contig_order {
            if let Some((r1_files, r2_files)) = all_fastq_files.remove(name) {
                all_r1.extend(r1_files);
                all_r2.extend(r2_files);
            }
        }
        if let Some((r1_files, r2_files)) = all_fastq_files.remove("chimeric") {
            all_r1.extend(r1_files);
            all_r2.extend(r2_files);
        }
        // Safety net for any unexpected leftover keys.
        for (r1_files, r2_files) in all_fastq_files.into_values() {
            all_r1.extend(r1_files);
            all_r2.extend(r2_files);
        }

        match &config.output_fastq_1 {
            Some(filename1) => {
                create_output_file(filename1, config.overwrite_output)?;
                if config.paired_ended {
                    match &config.output_fastq_2 {
                        Some(filename2) => {
                            create_output_file(filename2, config.overwrite_output)?;
                            combine_temp_fastqs(all_r1, all_r2, filename1, Some(filename2))?;
                        }
                        None => {
                            error!(
                                "Produce fastq true and paired-ended true, but output_fastq_2 was missing."
                            );
                            return Err(GenerateReadsError::ConfigError);
                        }
                    }
                } else {
                    combine_temp_fastqs(all_r1, vec![], filename1, None)?;
                }
            }
            None => {
                error!("Produce fastq true but output_fastq_1 was missing.");
                return Err(GenerateReadsError::ConfigError);
            }
        }
    }

    let mut files_written = Vec::new();
    if config.paired_ended {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", filename1);
            files_written.push(filename1.clone());
            if let Some(filename2) = &config.output_fastq_2 {
                info!("Successfully wrote fastq file: {:?}", filename2);
                files_written.push(filename2.clone());
            }
        }
    } else {
        if let Some(filename1) = &config.output_fastq_1 {
            info!("Successfully wrote fastq file: {:?}", filename1);
            files_written.push(filename1.clone());
        }
    }

    if config.produce_bam
        && let (Some(bam_ctx), Some(bam_path)) = (ctx.bam_context.as_ref(), &config.output_bam)
    {
        info!(
            "Assembling BAM from {} temp body file(s)",
            bam_body_files.len()
        );
        // Concatenate bodies in (contig order, chunk_start) order so the
        // assembled BAM stays coordinate-sorted across sub-contig chunks.
        let mut ordered_bodies: Vec<PathBuf> = Vec::new();
        for name in &contig_order {
            if let Some(mut bodies) = bam_body_files.remove(name) {
                bodies.sort_by_key(|(start, _)| *start);
                ordered_bodies.extend(bodies.into_iter().map(|(_, path)| path));
            }
        }
        // Chimeric reads are staged under a pseudo-contig so they do not appear in the BAM
        // header's reference dictionary. Their records still carry real contig coordinates and
        // must be appended after the regular per-contig bodies, just as their FASTQs are appended
        // in the chimeric pass above.
        if let Some(mut bodies) = bam_body_files.remove("chimeric") {
            bodies.sort_by_key(|(start, _)| *start);
            ordered_bodies.extend(bodies.into_iter().map(|(_, path)| path));
        }
        concat_temp_bams(bam_ctx, &ordered_bodies, bam_path)?;
        info!("Successfully wrote BAM file: {:?}", bam_path);
        files_written.push(bam_path.clone());
    }

    if let Some(filename) = &config.output_vcf {
        info!("Writing output vcf file");
        let result = write_vcf(
            &mutated_maps,
            &contig_order,
            &fasta_lengths,
            &config.reference,
            config.overwrite_output,
            filename,
            &ad_counters,
        );
        match result {
            Ok(()) => {
                info!("Successfully wrote vcf file: {:?}", filename);
                files_written.push(filename.clone());
            }
            Err(error) => {
                error!("Error writing vcf file!");
                return Err(GenerateReadsError::IoError(error));
            }
        }
    }
    Ok(files_written.clone())
}

/// Resolve the effective chunk size (bp) from config.
///
/// Sub-contig chunking is **disabled by default**: benchmarking showed eidolon's
/// read-generation loop is memory-bandwidth bound, so splitting a contig across
/// cores does not improve wall time (and can mildly regress it). The machinery
/// is kept behind an opt-in for CPU-bound or very-many-core scenarios.
///
/// - `None` (omitted)  → disabled: one chunk spans the whole contig (default).
/// - `Some(0)`         → disabled (explicit).
/// - `Some(n)` (n > 0) → fixed chunk size of `n` bp (opt-in).
fn resolve_chunk_size(config: &RunConfiguration) -> usize {
    match config.chunk_size {
        None | Some(0) => usize::MAX, // one chunk spans the whole contig
        Some(n) => n,
    }
}

/// One unit of parallel read-generation work: a sub-range of a contig.
struct ChunkWork {
    contig_idx: usize,
    chunk_idx: usize,
    name: String,
    chunk_start: usize,
    chunk_end: usize,
}

/// Result of generating reads for one chunk.
struct ChunkResult {
    contig_idx: usize,
    chunk_start: usize,
    name: String,
    len: usize, // full contig length (not chunk length)
    data: Option<ChunkData>,
}

struct ChunkData {
    r1_files: Vec<PathBuf>,
    r2_files: Vec<PathBuf>,
    bam_body_file: Option<PathBuf>,
    ad_counter: AdCounter,
}

/// Split a contig of `len` bp into evenly-sized chunks of ~`chunk_size` bp.
/// An empty contig yields a single `[0, 0)` chunk so it still produces a result.
fn split_contig_into_chunks(len: usize, chunk_size: usize) -> Vec<(usize, usize)> {
    if len == 0 {
        return vec![(0, 0)];
    }
    let chunk_size = chunk_size.max(1);
    let n = len.div_ceil(chunk_size).max(1);
    let base = len / n;
    let rem = len % n;
    let mut out = Vec::with_capacity(n);
    let mut start = 0;
    for i in 0..n {
        let this = base + usize::from(i < rem);
        out.push((start, start + this));
        start += this;
    }
    out
}

/// Lazily-built cache of full-contig `SequenceBlock`s shared across a contig's
/// chunks. Each contig's sequence is cloned at most once (the `OnceLock`
/// serializes the first build); all chunks of that contig share the `Arc`, so
/// memory is one sequence copy per *built* contig — not per chunk.
struct SeqBlockCache<'a> {
    reference: &'a HashMap<String, Vec<Nucleotide>>,
    cells: HashMap<String, std::sync::OnceLock<Arc<SequenceBlock>>>,
}

impl<'a> SeqBlockCache<'a> {
    fn new(reference: &'a HashMap<String, Vec<Nucleotide>>, names: &[String]) -> Self {
        let cells = names
            .iter()
            .map(|n| (n.clone(), std::sync::OnceLock::new()))
            .collect();
        Self { reference, cells }
    }

    /// Contigs are pre-validated (names come from the reference), so the build
    /// closure is infallible.
    fn get(&self, name: &str) -> Arc<SequenceBlock> {
        let cell = self
            .cells
            .get(name)
            .expect("chunk contig must be present in the sequence-block cache");
        Arc::clone(cell.get_or_init(|| {
            let sequence = self
                .reference
                .get(name)
                .expect("chunk contig must exist in reference")
                .clone();
            let sequence_map = map_buffer(&sequence);
            let ref_end = sequence.len();
            Arc::new(SequenceBlock {
                contig: name.to_string(),
                ref_start: 0,
                ref_end,
                sequence,
                sequence_map,
            })
        }))
    }
}

fn process_chunk(
    contig_idx: usize,
    chunk_idx: usize,
    contig_name: String,
    chunk_start: usize,
    chunk_end: usize,
    block: Arc<SequenceBlock>,
    ctx: &ContigContext,
    mut rng: NeatRng,
) -> Result<ChunkResult, GenerateReadsError> {
    let _ = chunk_idx; // reserved for diagnostics
    let contig_len = block.ref_end;
    debug!(
        "Processing {} chunk [{}, {})",
        contig_name, chunk_start, chunk_end
    );

    if let Some(bed) = ctx.target_bed
        && !bed.contains_key(&contig_name)
    {
        debug!("Skipping {} — not in target BED", contig_name);
        return Ok(ChunkResult {
            contig_idx,
            chunk_start,
            name: contig_name,
            len: contig_len,
            data: None,
        });
    }

    if block.sequence.is_empty() {
        warn!("Contig {} has empty sequence, skipping", contig_name);
        return Ok(ChunkResult {
            contig_idx,
            chunk_start,
            name: contig_name,
            len: contig_len,
            data: None,
        });
    }

    // The full-contig SequenceBlock is shared across this contig's chunks.
    let current_block = &*block;

    debug!("    > Generating bias map.");
    let raw_regions = current_block.get_non_n_regions();
    // Unclipped non-N interval bounds, captured BEFORE the BED intersection and
    // chunk clipping below, because those two narrow the intervals for
    // *ownership* while a fragment's END is deliberately allowed to run past
    // them (see the extension_budget doc on generate_fragments). An assembly
    // N-gap is the one boundary it must NOT run past: those bases are excluded
    // from read generation on purpose, and a fragment spilling into one emits
    // reads containing fabricated gap sequence. Bounding extension by the
    // enclosing non-N interval keeps chunk- and multiplier-boundary extension
    // (both wanted) while making a gap a true terminus (like the contig end).
    let non_n_bounds: Vec<(usize, usize)> = raw_regions.iter().map(|r| (r.start, r.end)).collect();
    let bed_regions: Vec<SequenceMap> = if let Some(bed) = ctx.target_bed {
        let contig_beds = bed.get(&contig_name).map(|v| v.as_slice()).unwrap_or(&[]);
        intersect_with_bed(&raw_regions, contig_beds, 0)
    } else {
        raw_regions.into_iter().cloned().collect()
    };
    // Restrict this chunk to fragments ANCHORED in [chunk_start, chunk_end).
    // Reads may still extend past chunk_end into the full shared sequence, so
    // no read-content stitching is needed and coverage stays uniform: every
    // fragment is owned by exactly one chunk (the one containing its anchor).
    let regions_of_interest: Vec<SequenceMap> = bed_regions
        .into_iter()
        .filter_map(|r| {
            let s = r.start.max(chunk_start);
            let e = r.end.min(chunk_end);
            (s < e).then(|| SequenceMap::from(r.region_type, s, e))
        })
        .collect();
    if regions_of_interest.is_empty() {
        return Ok(ChunkResult {
            contig_idx,
            chunk_start,
            name: contig_name,
            len: contig_len,
            data: None,
        });
    }

    // Build a compact segment list instead of a per-position Vec<f64>.
    // Each segment is (start, end, rate); N-regions and gaps are simply absent.
    // This replaces an O(chromosome_length) allocation with O(regions + BED_records).
    let mut rate_segments: Vec<(usize, usize, f64)> = regions_of_interest
        .iter()
        .map(|r| (r.start, r.end, ctx.default_run_mutation_rate))
        .collect();

    if let Some(mut_beds) = ctx.mutation_regions
        && let Some(records) = mut_beds.get(&contig_name)
    {
        for rec in records {
            if let Some(custom_rate) = rec.mut_rate {
                rate_segments = apply_rate_override(rate_segments, rec.start, rec.end, custom_rate);
            }
        }
    }

    // Borrow the shared per-contig MutatedMap (no per-chunk clone). The golden
    // VCF reads it back from the shared Arc after the parallel phase.
    let mutated_map = ctx.mutated_maps.get(&contig_name).ok_or_else(|| {
        GenerateReadsError::IoError(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("MutatedMap for {} not found", contig_name),
        ))
    })?;
    let max_del_len = *ctx.max_del_lens.get(&contig_name).unwrap_or(&0);

    // Keep short inserts when adapters are on (readthrough pads them) OR when the
    // adapter-free short-insert control is requested. Drives both fragment retention
    // (below) and the insert-length read cap in write_block_fastq.
    let keep_short = ctx.config.adapters.enabled || ctx.config.keep_short_fragments;

    // Altered haplotypes for long literal insertions in this chunk. Handed to
    // write_block_fastq alongside the fragments so ONE writer serves both
    // coordinate spaces (#516).
    let mut haplotypes: Vec<HaplotypeContext> = Vec::new();
    let block_fragments: Vec<PlacedFragment> = {
        let mut block_frags: Vec<PlacedFragment> = Vec::new();
        // SV coverage multipliers are needed here to scale fragment counts.
        // Even though they are also in MutatedMap, we need them as intervals.
        let sv_variants: Vec<Variant> = mutated_map.sv_records.iter().cloned().collect();
        let coverage_multipliers = build_coverage_multipliers(
            &sv_variants,
            ctx.config.ploidy,
            contig_len,
            ctx.config.subclone_model.is_some(),
        );

        for (region_start, region_end) in regions_of_interest.into_iter().map(|r| (r.start, r.end))
        {
            for (sub_start, sub_end, mult) in
                split_region_by_multipliers(region_start, region_end, &coverage_multipliers)
            {
                let scaled = scale_coverage(ctx.config.coverage, mult);
                if scaled == 0 {
                    continue;
                }
                // How far real, materializable sequence continues past this
                // sub-region's own right edge -- bounded by the enclosing non-N
                // interval, NOT by contig_len: an assembly gap's N bases are
                // excluded from read generation deliberately, so a fragment
                // extending across one would emit reads carrying fabricated gap
                // sequence (measured: 103 of 6000 reads contained N on a
                // reference with a 2kb gap, against 0 both before this branch
                // and at step 1, before extension was wired in). Chunk and
                // coverage-multiplier boundaries are still crossed freely --
                // only a gap (or the contig end) is a true terminus.
                // That measurably redistributes some of a narrow segment's
                // own declared depth to its neighbor (confirmed a real
                // redistribution, not a bug: total read count and R1/R2
                // totals are exactly conserved against a Terminal-everywhere
                // placement; ownership counts and cover_dataset's own
                // placement are both exact) -- and confirmed to vanish at
                // realistic scale: a 1200bp event on H1N1's ~400bp flanks
                // needed a 1.13x correction, the same event on a real chr22
                // window (megabase flanks) needed 1.07x, a realistic 100kb
                // event needed 1.003x. No compensation belongs here; a
                // multiplier-boundary guard was tried and reverted, since it
                // fixes small-scale depth accuracy at the cost of
                // reintroducing the zero/near-zero-output cliff for SVs
                // narrower than the fragment-length scale. See
                // docs/claude_engineering_audit.md §5.6's 2026-08-22 addendum
                // and docs/sv_polish_roadmap.md's Phase 1 item 1.
                let materializable_end = non_n_bounds
                    .iter()
                    .find(|&&(ns, ne)| sub_start >= ns && sub_start < ne)
                    .map(|&(_, ne)| ne)
                    .unwrap_or(contig_len);
                let extension_budget = materializable_end.saturating_sub(sub_end);

                // A long literal insertion gives this sub-region a SECOND, longer
                // molecule to sample from. Model it that way rather than patching
                // reference-sampled reads afterwards: draw the alt haplotype's share
                // of coverage in haplotype coordinates (where the inserted sequence
                // has width, so fragments can begin inside it -- the whole of #516)
                // and the rest from the reference. Zygosity then falls out of the
                // sampling instead of needing a per-read coin, which is what the
                // reverted attempt got wrong: it emitted every insertion at full
                // depth regardless of genotype (measured het/hom ratio 1.02).
                //
                // Only insertions at least a read long need this. Shorter ones are
                // already fully realized by inline variant application, and routing
                // them here would change output for no benefit.
                let long_ins: Vec<&Variant> = mutated_map
                    .variant_map
                    .values()
                    .filter(|v| {
                        v.variant_type == VariantType::Insertion
                            && v.location >= sub_start
                            && v.location < sub_end
                            && v.alternate
                                .as_literal()
                                .and_then(|a| a.len().checked_sub(v.reference.len()))
                                .is_some_and(|novel| novel >= ctx.config.read_len)
                    })
                    .collect();

                // Build ONE alt haplotype carrying every long insertion in this
                // sub-region. An earlier version handled a single insertion and fell
                // back to the pre-#516 head-only behaviour for the rest; because
                // sub-regions are large (a whole contig, split only by coverage
                // multipliers), that caught the ordinary case of planting a size range
                // at once via input_vcf -- measured on 200/600/1200bp sharing a contig:
                // middle and tail both zero for all three.
                //
                // Insertions are grouped by their alt fraction. Two insertions at the
                // same fraction belong on the same molecule; two at different fractions
                // (a het and a hom, say) do not, and giving each group its own haplotype
                // keeps each event's dosage independent rather than averaging them.
                // SV-scale literal deletions belong on a haplotype for the mirror-image
                // reason (#590). The per-base path can only SKIP deleted bases while
                // rendering a read; it cannot stop fragments being PLACED across a span
                // that, on this molecule, does not exist. Measured on a homozygous
                // deletion against a no-variant control, depth inside the deleted span
                // as a fraction of flank: 0.42 at 120 bp and 0.70 at 500 bp, where the
                // identical event written as symbolic <DEL> gives 0.00.
                //
                // MIN_SV_LENGTH_BP is the codebase's own SV threshold, and the line is
                // drawn there rather than at read_len (as insertions are) because the
                // per-base path fails well before a read length. Below it the residual
                // is small (measured 0.09 of flank at 30 bp) and routing every
                // indel-model deletion through a materialized haplotype would change
                // every run for no measurable gain.
                let sv_dels: Vec<&Variant> = mutated_map
                    .variant_map
                    .values()
                    .filter(|v| {
                        v.location >= sub_start
                            && v.location < sub_end
                            && v.alternate.as_literal().is_some_and(|alt| {
                                v.reference.len().saturating_sub(alt.len()) >= MIN_SV_LENGTH_BP
                            })
                    })
                    .collect();

                let mut ref_cov = scaled;
                let mut alt_plans: Vec<(usize, usize)> = Vec::new(); // (hap_idx, alt_cov)
                if !long_ins.is_empty() || !sv_dels.is_empty() {
                    // Mirror the allele semantics generate_read uses for a literal
                    // variant, so haplotype sampling and inline application agree: an
                    // explicit allele_fraction wins, otherwise homozygous is always alt
                    // and heterozygous is a half-and-half split.
                    let fraction_of = |v: &Variant| -> f64 {
                        match v.allele_fraction {
                            Some(f) => f.clamp(0.0, 1.0),
                            None => match v.genotype {
                                Genotype::Homozygous => 1.0,
                                Genotype::Heterozygous => 0.5,
                            },
                        }
                    };
                    // Group by fraction, keyed to the nearest thousandth so ordinary
                    // float noise does not split a group that is conceptually one.
                    let mut groups: std::collections::BTreeMap<u64, Vec<&Variant>> =
                        std::collections::BTreeMap::new();
                    // Insertions and deletions share the grouping: two events at the same
                    // alt fraction belong on the SAME molecule, whichever kind they are.
                    for v in long_ins.iter().chain(sv_dels.iter()) {
                        let key = (fraction_of(v) * 1000.0).round() as u64;
                        groups.entry(key).or_default().push(v);
                    }
                    let mut spent = 0usize;
                    for (key, members) in groups {
                        let f = key as f64 / 1000.0;
                        let mut ins_entries: Vec<(usize, Vec<Nucleotide>)> = Vec::new();
                        let mut del_entries: Vec<(usize, usize)> = Vec::new();
                        for v in &members {
                            let Some(alt) = v.alternate.as_literal() else {
                                continue;
                            };
                            if alt.len() > v.reference.len() {
                                if let Some(novel) = alt.get(v.reference.len()..) {
                                    ins_entries.push((v.location, novel.to_vec()));
                                }
                            } else if v.reference.len() > alt.len() {
                                del_entries.push((v.location, v.reference.len() - alt.len()));
                            }
                        }
                        let Some(map) = InsertionCoordinateMap::with_deletions(
                            contig_len,
                            ins_entries,
                            del_entries,
                        ) else {
                            // Refused only for genuinely ambiguous input (two insertions
                            // anchored at the same base, an anchor off the contig). Say
                            // so: these events keep the pre-#516 head-only behaviour and
                            // must not be mistaken for working ones.
                            warn!(
                                "{contig_name}: could not build a haplotype for {} long \
                                 insertion(s) in [{sub_start},{sub_end}) — they keep the \
                                 pre-#516 head-only behaviour (#516)",
                                members.len()
                            );
                            continue;
                        };
                        let alt_cov = scale_coverage(scaled, f);
                        if alt_cov == 0 {
                            continue;
                        }
                        spent += alt_cov;
                        haplotypes.push(HaplotypeContext { map });
                        alt_plans.push((haplotypes.len() - 1, alt_cov));
                    }
                    ref_cov = scaled.saturating_sub(spent);
                }

                for (hap_idx, alt_cov) in alt_plans {
                    // The same sub-region expressed on the alt haplotype. Ask the map
                    // to project both boundaries rather than computing the shift by
                    // hand: with several insertions the span grows by however many of
                    // them fall inside, and that arithmetic is exactly what the map
                    // exists to get right. The extra width is what lets a fragment
                    // BEGIN inside inserted sequence, which reference coordinates
                    // cannot express at all.
                    let h = &haplotypes[hap_idx];
                    let hap_start = h
                        .map
                        .reference_base_to_haplotype(sub_start)
                        .unwrap_or(sub_start);
                    // sub_end is exclusive, so project the last included base and add
                    // one; projecting sub_end itself would fall off the contig end.
                    let hap_end = h
                        .map
                        .reference_base_to_haplotype(sub_end.saturating_sub(1))
                        .map(|p| p + 1)
                        .unwrap_or_else(|| h.map.haplotype_len());
                    let hap_span = hap_end.saturating_sub(hap_start);
                    let alt_frags = generate_fragments(
                        extension_budget,
                        hap_span,
                        ctx.config.read_len,
                        max_del_len,
                        hap_start,
                        alt_cov,
                        ctx.config.paired_ended,
                        ctx.config.long_reads,
                        keep_short,
                        ctx.fragment_length_model,
                        &mut rng,
                    )?;
                    block_frags.extend(alt_frags.into_iter().map(|(s, e)| PlacedFragment {
                        start: s,
                        end: e,
                        haplotype: Some(hap_idx),
                    }));
                }

                let scaled = ref_cov;
                if scaled == 0 {
                    continue;
                }
                let frags = if ctx.gc_bias_model.is_uniform() {
                    generate_fragments(
                        extension_budget,
                        sub_end - sub_start,
                        ctx.config.read_len,
                        max_del_len,
                        sub_start,
                        scaled,
                        ctx.config.paired_ended,
                        ctx.config.long_reads,
                        keep_short,
                        ctx.fragment_length_model,
                        &mut rng,
                    )?
                } else {
                    generate_weighted_fragments(
                        extension_budget,
                        current_block,
                        sub_start,
                        sub_end,
                        ctx.config.read_len,
                        max_del_len,
                        scaled,
                        ctx.gc_bias_model,
                        ctx.fragment_length_model,
                        ctx.config.gc_bias_normalize_coverage,
                        ctx.config.paired_ended,
                        ctx.config.long_reads,
                        keep_short,
                        &mut rng,
                    )?
                };
                if !frags.is_empty() {
                    block_frags.extend(frags.into_iter().map(PlacedFragment::from));
                }
            }
        }
        block_frags
    };

    // Breakpoint double-counting fix: the chimeric pass emits junction-spanning
    // read-pairs for every chimeric SV junction (BND/INV, plus DEL and CNV-loss),
    // but the regular pass above also covers those breakpoints from the unbroken
    // reference (BND/INV are coverage-neutral; DEL/CNV-loss only zero the deleted
    // interior — flank reads crossing the breakpoint still leak). Left alone, a
    // homozygous junction lands at ~2x coverage (regular + junction reads), a
    // het at ~1.5x. Drop the broken-allele fraction of regular pairs that cross a
    // junction so total junction depth ≈ coverage. DUP / CNV-gain create a novel
    // tandem adjacency that linear reads never reproduce, so they are not
    // suppressed. No-op (no RNG drawn) when the contig has no suppressible SVs.
    let block_fragments = suppress_junction_double_count(
        block_fragments,
        &mutated_map.sv_records,
        ctx.config.ploidy,
        ctx.config.read_len,
        ctx.config.paired_ended,
        ctx.config.subclone_model.is_some(),
        &mut rng,
    )?;

    let mut contig_files_r1: Vec<PathBuf> = Vec::new();
    let mut contig_files_r2: Vec<PathBuf> = Vec::new();

    // Per-contig allelic-depth counter. Threaded into every write_block_fastq
    // call below; each variant-overlapping read increments the (ref|alt) slot
    // for that variant. The fully-populated counter is handed off via
    // ProcessedContigData → run_neat → write_vcf for FORMAT/AD/DP/AF.
    let mut ad_counter: AdCounter = AdCounter::new();

    // Create a per-contig BAM body writer if BAM output is requested.
    let mut bam_body_writer: Option<BamBodyWriter> = if let Some(bam_ctx) = &ctx.bam_context {
        let bam_temp_path = PathBuf::from(ctx.working_dir).join(format!(
            "temp_bam_{:06}_{:010}_{}.bam",
            contig_idx, chunk_start, contig_name
        ));
        Some(BamBodyWriter::new(bam_temp_path, Arc::clone(bam_ctx))?)
    } else {
        None
    };

    let read_name_prefix = format!("EIDOLON_generated_{}", current_block.contig);

    // Resolve 3' adapter sequences once (#125). Empty vecs = disabled, which makes
    // write_block_fastq take its unchanged code path (output byte-identical when off).
    let (r1_adapter, r2_adapter): (Vec<Nucleotide>, Vec<Nucleotide>) =
        if ctx.config.adapters.enabled {
            (
                ctx.config
                    .adapters
                    .r1
                    .chars()
                    .map(Nucleotide::from)
                    .collect(),
                ctx.config
                    .adapters
                    .r2
                    .chars()
                    .map(Nucleotide::from)
                    .collect(),
            )
        } else {
            (Vec::new(), Vec::new())
        };

    if ctx.config.produce_fastq {
        let mut file_to_write_1 = PathBuf::from(ctx.working_dir);
        file_to_write_1.push(format!(
            "temp_{}_{:010}_{:010}_r1_tmp.fastq.gz",
            contig_name, chunk_start, chunk_end,
        ));
        let file1 = append_to_file(&file_to_write_1)?;
        let writer1 = BufWriter::new(&file1);
        let mut buffer1 = BlockGzWriter::new(writer1);
        let bam_stager: Option<&mut dyn BamRecordStager> = bam_body_writer
            .as_mut()
            .map(|w| w as &mut dyn BamRecordStager);
        if ctx.config.paired_ended {
            let mut file_to_write_2 = PathBuf::from(ctx.working_dir);
            file_to_write_2.push(format!(
                "temp_{}_{:010}_{:010}_r2_tmp.fastq.gz",
                contig_name, chunk_start, chunk_end,
            ));
            let file2 = append_to_file(&file_to_write_2)?;
            let writer2 = BufWriter::new(&file2);
            let mut buffer2 = BlockGzWriter::new(writer2);
            debug!("Writing paired-ended contig fastq files");
            write_block_fastq(
                block_fragments.into_iter().map(Into::into).collect(),
                &haplotypes,
                mutated_map,
                current_block,
                true,
                &mut buffer1,
                &mut buffer2,
                ctx.config.read_len,
                ctx.config.long_reads,
                keep_short,
                &read_name_prefix,
                ctx.quality_score_model,
                ctx.seq_error_model,
                &mut rng,
                bam_stager,
                &mut ad_counter,
                &r1_adapter,
                &r2_adapter,
            )?;
            contig_files_r1.push(file_to_write_1);
            contig_files_r2.push(file_to_write_2);
        } else {
            debug!("Writing single-ended contig fastq file");
            let dummy_data: VectorBuffer = VectorBuffer::new();
            let mut buffer2 = GzEncoder::new(dummy_data, Compression::default());
            write_block_fastq(
                block_fragments.into_iter().map(Into::into).collect(),
                &haplotypes,
                mutated_map,
                current_block,
                false,
                &mut buffer1,
                &mut buffer2,
                ctx.config.read_len,
                ctx.config.long_reads,
                keep_short,
                &read_name_prefix,
                ctx.quality_score_model,
                ctx.seq_error_model,
                &mut rng,
                bam_stager,
                &mut ad_counter,
                &r1_adapter,
                &r2_adapter,
            )?;
            contig_files_r1.push(file_to_write_1);
        }
    } else if ctx.config.produce_bam {
        // BAM-only: generate reads and stage them into the BAM body writer.
        // The FASTQ bytes are discarded, so write them to a null sink rather
        // than compressing them into a throwaway buffer (the records still flow
        // to the BAM writer via bam_stager).
        let bam_stager: Option<&mut dyn BamRecordStager> = bam_body_writer
            .as_mut()
            .map(|w| w as &mut dyn BamRecordStager);
        let mut buf1 = std::io::sink();
        let mut buf2 = std::io::sink();
        debug!("BAM-only: generating reads for {}", contig_name);
        write_block_fastq(
            block_fragments.into_iter().map(Into::into).collect(),
            &haplotypes,
            mutated_map,
            current_block,
            ctx.config.paired_ended,
            &mut buf1,
            &mut buf2,
            ctx.config.read_len,
            ctx.config.long_reads,
            keep_short,
            &read_name_prefix,
            ctx.quality_score_model,
            ctx.seq_error_model,
            &mut rng,
            bam_stager,
            &mut ad_counter,
            &r1_adapter,
            &r2_adapter,
        )?;
    }

    // Flush and finalize the BAM body file; the bgzf EOF is written on drop.
    let bam_body_file = if let Some(mut bw) = bam_body_writer {
        bw.flush_all()?;
        Some(bw.path.clone())
    } else {
        None
    };

    Ok(ChunkResult {
        contig_idx,
        chunk_start,
        name: contig_name,
        len: contig_len,
        data: Some(ChunkData {
            r1_files: contig_files_r1,
            r2_files: contig_files_r2,
            bam_body_file,
            ad_counter,
        }),
    })
}

/// Return the fraction of cellular copies carrying an SV. For CNVs this is the
/// magnitude of the copy-number deviation; for other SVs it is the genotype
/// dosage. This is the quantity that CCF scales for the tumor pass.
fn sv_dosage_fraction(v: &Variant, ploidy: usize) -> f64 {
    if let Some(cn) = v.alternate.as_symbolic().and_then(|sv| sv.copy_number) {
        let p = ploidy.max(1) as f64;
        return (cn as f64 - p).abs() / p;
    }
    v.dosage_fraction()
}

/// Apply the existing subclone model to de-novo SVs only. Input/germline SVs
/// remain unchanged, and the no-model path consumes no additional RNG draws.
fn apply_sv_subclone_model(
    variants: &mut [Variant],
    model: Option<&SubcloneModel>,
    ploidy: usize,
    merged_vaf_purity: Option<f64>,
    rng: &mut NeatRng,
) -> Result<(), GenerateReadsError> {
    let Some(model) = model else {
        return Ok(());
    };
    for variant in variants {
        let ccf = model.sample_ccf(rng).map_err(GenerateReadsError::from)?;
        stamp_sv_subclone(variant, ccf, ploidy, merged_vaf_purity);
    }
    Ok(())
}

fn stamp_sv_subclone(
    variant: &mut Variant,
    ccf: f64,
    ploidy: usize,
    merged_vaf_purity: Option<f64>,
) {
    let af = sv_dosage_fraction(variant, ploidy) * ccf;
    variant.allele_fraction = Some(af);
    append_info_tag(&mut variant.info, format!("EIDOLON_CCF={ccf:.4}"));
    if let Some(purity) = merged_vaf_purity {
        append_info_tag(&mut variant.info, format!("EIDOLON_VAF={:.4}", purity * af));
    }
}

/// Stamp both ends of each genome-wide translocation with one shared CCF. A
/// translocation is one biological event even though it has two VCF records;
/// sampling the ends independently would create an impossible allele balance.
fn apply_translocation_subclone_model(
    translocations: &mut HashMap<String, Vec<Variant>>,
    model: Option<&SubcloneModel>,
    ploidy: usize,
    merged_vaf_purity: Option<f64>,
    rng: &mut NeatRng,
) -> Result<(), GenerateReadsError> {
    let Some(model) = model else {
        return Ok(());
    };
    // Keep pair traversal ordered so CCF assignment remains reproducible across
    // processes (HashSet iteration order is intentionally randomized).
    let mut pairs = BTreeSet::new();
    for variants in translocations.values() {
        for variant in variants {
            let Some(id) = variant.id.as_ref() else {
                continue;
            };
            let mate = variant.info.as_deref().and_then(|info| {
                info.split(';')
                    .find_map(|field| field.strip_prefix("MATEID="))
            });
            let Some(mate) = mate else { continue };
            let key = if id.as_str() <= mate {
                (id.clone(), mate.to_string())
            } else {
                (mate.to_string(), id.clone())
            };
            pairs.insert(key);
        }
    }
    for (id, mate) in pairs {
        let ccf = model.sample_ccf(rng).map_err(GenerateReadsError::from)?;
        for variants in translocations.values_mut() {
            for variant in variants {
                if variant.id.as_deref() == Some(id.as_str())
                    || variant.id.as_deref() == Some(mate.as_str())
                {
                    stamp_sv_subclone(variant, ccf, ploidy, merged_vaf_purity);
                }
            }
        }
    }
    Ok(())
}

/// Effective alternate fraction for SV read evidence. With no subclone model,
/// this deliberately returns the historical genotype-based value.
fn sv_effective_fraction(v: &Variant, ploidy: usize, use_subclone: bool) -> f64 {
    if use_subclone {
        v.allele_fraction
            .unwrap_or_else(|| sv_dosage_fraction(v, ploidy))
    } else {
        match v.genotype {
            Genotype::Homozygous => 1.0,
            Genotype::Heterozygous => 1.0 / (ploidy.max(1) as f64),
        }
    }
}

fn generate_mutated_map(
    contig_name: &str,
    sequence: &[Nucleotide],
    config: &RunConfiguration,
    target_bed: &Option<HashMap<String, Vec<BedRecord>>>,
    mutation_regions: &Option<HashMap<String, Vec<BedRecord>>>,
    input_variants: &Option<HashMap<String, Vec<Variant>>>,
    // Breakend records already placed for this contig by the genome-wide translocation
    // pass. Seeded into `sv_variants` before de novo sampling so overlap rejection sees
    // them, exactly like input-VCF SVs.
    preplaced_svs: Option<&Vec<Variant>>,
    mutation_model: &MutationModel,
    default_run_mutation_rate: f64,
    mut rng: NeatRng,
) -> Result<(MutatedMap, usize), GenerateReadsError> {
    let contig_len = sequence.len();
    if contig_len == 0 {
        return Ok((
            MutatedMap::from_interval(0, 0, vec![]).map_err(GenerateReadsError::from)?,
            0,
        ));
    }

    let sequence_map = map_buffer(sequence);
    let current_block = SequenceBlock {
        contig: contig_name.to_string(),
        ref_start: 0,
        ref_end: contig_len,
        sequence: sequence.to_vec(),
        sequence_map,
    };

    let raw_regions = current_block.get_non_n_regions();
    let regions_of_interest: Vec<SequenceMap> = if let Some(bed) = target_bed {
        let contig_beds = bed.get(contig_name).map(|v| v.as_slice()).unwrap_or(&[]);
        intersect_with_bed(&raw_regions, contig_beds, 0)
    } else {
        raw_regions.into_iter().cloned().collect()
    };

    if regions_of_interest.is_empty() {
        return Ok((
            MutatedMap::from_interval(0, contig_len, vec![]).map_err(GenerateReadsError::from)?,
            0,
        ));
    }

    let mut rate_segments: Vec<(usize, usize, f64)> = regions_of_interest
        .iter()
        .map(|r| (r.start, r.end, default_run_mutation_rate))
        .collect();

    if let Some(mut_beds) = mutation_regions
        && let Some(records) = mut_beds.get(contig_name)
    {
        for rec in records {
            if let Some(custom_rate) = rec.mut_rate {
                rate_segments = apply_rate_override(rate_segments, rec.start, rec.end, custom_rate);
            }
        }
    }

    let mut num_mutations_sum: f64 = rate_segments
        .iter()
        .map(|&(s, e, r)| (e - s) as f64 * r)
        .sum();

    let mut block_variants: Vec<Variant> = Vec::new();
    let mut sv_variants: Vec<Variant> = Vec::new();
    // Genome-wide translocations land here first: they are already placed, and seeding
    // them before de novo sampling makes overlap rejection treat them as occupied.
    if let Some(pre) = preplaced_svs {
        sv_variants.extend(pre.iter().cloned());
    }
    if let Some(iv) = input_variants
        && let Some(vs) = iv.get(contig_name)
    {
        let mut excluded: Vec<usize> = Vec::new();
        let mut seen: HashSet<usize> = HashSet::new();
        for v in vs {
            let pos0 = v.location.saturating_sub(1);
            if pos0 >= contig_len {
                continue;
            }
            let local_pos = v.location - 1;
            let mut v2 = v.clone();
            v2.location = local_pos;
            if v.alternate.is_symbolic() {
                sv_variants.push(v2);
                continue;
            }
            if seen.insert(local_pos) {
                let rate = rate_at(&rate_segments, local_pos);
                if rate > 0.0 {
                    num_mutations_sum -= rate;
                    excluded.push(local_pos);
                }
                block_variants.push(v2);
            }
        }
        if !excluded.is_empty() {
            excluded.sort_unstable();
            rate_segments = exclude_positions(rate_segments, &excluded);
        }
    }

    if config.sv_rate_scale > 0.0
        && let Some(sv_model) = mutation_model.sv_model.as_ref()
        && sv_model.is_usable()
    {
        let de_novo = sv_model.sample_variants(
            contig_name,
            contig_len,
            &sv_variants,
            sequence,
            config.ploidy,
            config.sv_rate_scale,
            config.sv_max_length_fraction,
            &mut rng,
        );
        // The #516 interim cap that used to sit here is GONE. It dropped every de novo
        // insertion longer than `read_len - 1` on the grounds that "reads carry at most
        // read_len - 1 of an insertion's novel bases (fragments are placed in reference
        // offsets, so none can begin inside a zero-reference-width event)", and it said
        // it would stand "until the fragment sampler can place reads in haplotype
        // coordinates". That is now exactly what happens: a long insertion is sampled on
        // its own altered haplotype, where the inserted sequence has width and fragments
        // begin inside it. Keeping the cap would leave the de novo path unable to reach
        // the fix at all — it was the reason a de novo campaign could never plant a long
        // insertion no matter what the model asked for.
        //
        // The cap's own argument is what retires it, so it is removed rather than
        // raised: there is no longer a length beyond which the reads cannot support the
        // declared SVLEN. `drop_unrealizable_insertions` and its tests are deleted with
        // it — an unused guard kept behind an allow(dead_code) is the "silently
        // disabled" shape this repo's deny-warnings policy exists to catch. The distinct
        // question of an insertion too large to MATERIALIZE is already handled generally
        // by `sv_max_length_fraction`, which bounds de novo SV length to a fraction of
        // the contig (default 0.25).
        let mut de_novo = de_novo;
        apply_sv_subclone_model(
            &mut de_novo,
            config.subclone_model.as_ref(),
            config.ploidy,
            config.merged_vaf_purity,
            &mut rng,
        )?;
        sv_variants.extend(de_novo);
    }

    let coverage_multipliers = build_coverage_multipliers(
        &sv_variants,
        config.ploidy,
        contig_len,
        config.subclone_model.is_some(),
    );
    let mut zeroed = false;
    for &(s, e, mult) in &coverage_multipliers {
        if mult == 0.0 && s < e {
            rate_segments = apply_rate_override(rate_segments, s, e, 0.0);
            zeroed = true;
        }
    }
    if zeroed {
        num_mutations_sum = rate_segments
            .iter()
            .map(|&(s, e, r)| (e - s) as f64 * r)
            .sum();
    }

    let mut max_del_len = 0;
    for v in &block_variants {
        if v.variant_type == VariantType::Deletion && v.reference.len() > 1 {
            max_del_len = max_del_len.max(v.reference.len() - 1);
        }
    }

    let num_mutations = num_mutations_sum.trunc() as usize;
    if num_mutations > 0 {
        let result = generate_variants(
            &current_block,
            &rate_segments,
            mutation_model,
            num_mutations,
            config.ploidy,
            &mut rng,
        )?;
        if let Some(vec) = result {
            for mut variant in vec {
                // #405: in the cancer tumor pass, distribute de-novo somatic variants
                // across subclones. A subclone's CCF is a *cellular-fraction factor*
                // that composes with the variant's allele dosage — it does not replace
                // it. So observed alt fraction = dosage × CCF (× purity via the
                // tumor/normal coverage split at merge time): a heterozygous somatic
                // SNV at CCF f lands at f/2, which is what subclonal-deconvolution
                // tools invert. Multiplying (not overwriting) also lets polyploid
                // dosage (#266/#267) flow in for free once genotype_str carries a real
                // per-copy spread. `None` elsewhere leaves output byte-identical.
                if let Some(model) = &config.subclone_model {
                    let ccf = model.sample_ccf(&mut rng)?;
                    let base = variant
                        .allele_fraction
                        .unwrap_or_else(|| variant.dosage_fraction());
                    let af = base * ccf;
                    variant.allele_fraction = Some(af);
                    // Ground truth in the golden VCF: EIDOLON_CCF = intended cellular
                    // fraction (observed AD/AF tracks dosage × CCF within the tumor
                    // pass); EIDOLON_VAF = intended observed fraction after tumor/normal
                    // mixing (purity × af), directly comparable to a caller's VAF.
                    append_info_tag(&mut variant.info, format!("EIDOLON_CCF={ccf:.4}"));
                    if let Some(p) = config.merged_vaf_purity {
                        append_info_tag(&mut variant.info, format!("EIDOLON_VAF={:.4}", p * af));
                    }
                }
                if variant.variant_type == VariantType::Deletion
                    && variant.reference.len() - 1 > max_del_len
                {
                    max_del_len = variant.reference.len() - 1;
                }
                block_variants.push(variant);
            }
        }
    }

    block_variants.extend(sv_variants);
    let mutated_map = MutatedMap::from_interval(0, contig_len, block_variants)
        .map_err(GenerateReadsError::from)?;
    Ok((mutated_map, max_del_len))
}

fn process_chimeric_variants(
    ctx: &ContigContext,
    mut rng: NeatRng,
) -> Result<ChunkResult, GenerateReadsError> {
    let mut all_reads = Vec::new();
    let mut processed_ids = HashSet::new();

    // Reference-file order, NOT HashMap order — see `ContigContext::contig_order` (#599).
    // A single `rng` is threaded through this whole loop, so the order in which contigs
    // are visited decides which draws each junction gets; a randomized order changed the
    // emitted read set, and even its size, on every run.
    for contig_name in ctx.contig_order.iter() {
        let Some(m_map) = ctx.mutated_maps.get(contig_name) else {
            continue;
        };
        for sv_rec in &m_map.sv_records {
            let sv = match sv_rec.alternate.as_symbolic() {
                Some(s) => s,
                _ => continue,
            };

            if sv.sv_type == SvType::Bnd {
                let mate_contig = match &sv.mate_contig {
                    Some(c) => c,
                    None => continue,
                };
                let mate_pos = match sv.mate_pos {
                    Some(p) => p,
                    None => continue,
                };

                // BND records come in pairs — each side describes the same
                // junction from its own contig+position. Canonicalize the
                // (contig, pos) tuples so the "smaller" side comes first,
                // then use that as the dedup key. Tuple ordering handles
                // both the cross-contig case (compare contig names
                // lexicographically) and the same-contig case (compare
                // positions) uniformly. Stored in `processed_ids` keyed by
                // type-prefixed string so BND and INV share one HashSet.
                //
                // BOTH sides must be expressed in the SAME coordinate base or the
                // two records canonicalize to different keys and the junction is
                // processed twice (2x chimeric reads). `sv_rec.location` is
                // 0-based; `sv.mate_pos` is 1-based (it comes from the VCF ALT
                // string, and line ~2333 does saturating_sub(1) before indexing
                // the sequence). Normalize the mate to 0-based here.
                let here = (contig_name.as_str(), sv_rec.location);
                let mate_pos_0based = mate_pos.saturating_sub(1);
                let mate = (mate_contig.as_str(), mate_pos_0based);
                // The key body must also use the normalized coordinate, not the
                // 1-based mate_pos, or the two sides still produce different keys.
                let bnd_id = if here <= mate {
                    (
                        contig_name.clone(),
                        sv_rec.location,
                        mate_contig.clone(),
                        mate_pos_0based,
                    )
                } else {
                    (
                        mate_contig.clone(),
                        mate_pos_0based,
                        contig_name.clone(),
                        sv_rec.location,
                    )
                };
                if !processed_ids.insert(format!("BND_{:?}", bnd_id)) {
                    continue;
                }

                // Coverage model for chimeric BND reads:
                //   - Homozygous BND: every allele carries the junction, so
                //     every read covering the breakpoint should be a
                //     junction read. mult = 1.0 → num_frags = full coverage.
                //   - Heterozygous: half the alleles carry the junction; the
                //     other half are unbroken reference. mult = 1/ploidy.
                //
                // Important caveat: the *regular* per-contig pass also
                // generates reads covering the breakpoint position (it just
                // reads from the unbroken reference, oblivious to the BND).
                // For a homozygous BND this means the breakpoint locus ends
                // up with regular reads PLUS junction reads — roughly double
                // the true biological coverage there. For a heterozygous
                // BND it lands closer to correct (regular pass covers both
                // alleles, chimeric pass adds 1/ploidy worth of junction
                // reads). A proper fix would teach the regular pass to skip
                // the broken-allele fraction of reads at BND positions;
                // tracked as a v2 follow-up.
                let mult = sv_effective_fraction(
                    sv_rec,
                    ctx.config.ploidy,
                    ctx.config.subclone_model.is_some(),
                );

                let num_frags = scale_coverage(ctx.config.coverage, mult);
                if num_frags == 0 {
                    continue;
                }

                for frag_idx in 0..num_frags {
                    // Fragment length picked the same way the main read-gen
                    // path does it (paired-end: model-sampled with a
                    // length-floor retry; single-end: read_len + a 32bp pad
                    // so sequencing-error deletions don't truncate the read).
                    let se_pad = if ctx.config.paired_ended { 0 } else { 32 };
                    let frag_len = if ctx.config.paired_ended {
                        let mut attempts = 0;
                        let mut f = 0;
                        while attempts < 100 {
                            let rand_val = rng.random().map_err(GenerateReadsError::from)?;
                            f = ctx
                                .fragment_length_model
                                .generate_fragment(rand_val)
                                .map_err(GenerateReadsError::from)?
                                as usize;
                            if ctx.config.long_reads || f >= ctx.config.read_len + 10 {
                                break;
                            }
                            attempts += 1;
                        }
                        if f < ctx.config.read_len && !ctx.config.long_reads {
                            f = ctx.config.read_len + 10;
                        }
                        f
                    } else {
                        ctx.config.read_len + se_pad
                    };

                    let offset = balanced_chimeric_offset(frag_len, ctx.config.read_len, &mut rng)?;

                    let result = generate_chimeric_pair(
                        ctx,
                        contig_name,
                        sv_rec.location,
                        sv,
                        frag_len,
                        offset,
                        frag_idx,
                        &mut rng,
                    );

                    match result {
                        Ok((read1, read2)) => {
                            all_reads.push(read1);
                            if let Some(r2) = read2 {
                                all_reads.push(r2);
                            }
                        }
                        Err(GenerateReadsError::FqToolsError(
                            eidolon_core::file_tools::fastq_tools::FastqToolsError::TruncatedRead(
                                msg,
                            ),
                        )) => {
                            debug!("Skipping truncated chimeric read: {}", msg);
                        }
                        Err(e) => return Err(e),
                    }
                }
            } else if sv.sv_type == SvType::Inv {
                let end = match sv.end {
                    Some(e) => e,
                    None => {
                        if let Some(span) = sv.span(sv_rec.location) {
                            sv_rec.location + span - 1
                        } else {
                            continue;
                        }
                    }
                };

                let inv_id = (contig_name.clone(), sv_rec.location, end);
                if !processed_ids.insert(format!("INV_{:?}", inv_id)) {
                    continue;
                }

                // Same coverage model as BND — full coverage for homozygous,
                // 1/ploidy for heterozygous. The double-counting caveat from
                // the BND branch applies here too: the regular per-contig
                // pass still generates reads spanning the inversion's two
                // breakpoints (it reads from the unbroken forward reference),
                // so a homozygous inversion ends up with regular + junction
                // coverage at each breakpoint.
                let mult = sv_effective_fraction(
                    sv_rec,
                    ctx.config.ploidy,
                    ctx.config.subclone_model.is_some(),
                );

                let num_frags = scale_coverage(ctx.config.coverage, mult);
                if num_frags == 0 {
                    continue;
                }

                for frag_idx in 0..num_frags {
                    let se_pad = if ctx.config.paired_ended { 0 } else { 32 };
                    let frag_len = if ctx.config.paired_ended {
                        let mut attempts = 0;
                        let mut f = 0;
                        while attempts < 100 {
                            let rand_val = rng.random().map_err(GenerateReadsError::from)?;
                            f = ctx
                                .fragment_length_model
                                .generate_fragment(rand_val)
                                .map_err(GenerateReadsError::from)?
                                as usize;
                            if ctx.config.long_reads || f >= ctx.config.read_len + 10 {
                                break;
                            }
                            attempts += 1;
                        }
                        if f < ctx.config.read_len && !ctx.config.long_reads {
                            f = ctx.config.read_len + 10;
                        }
                        f
                    } else {
                        ctx.config.read_len + se_pad
                    };

                    // An inversion has two breakpoints (junction=1 at the
                    // start, junction=2 at the end), and each one needs
                    // its own junction-spanning reads. Both junctions share
                    // the frag_idx (it's a per-INV-record counter) but the
                    // junction number disambiguates the QNAMEs that
                    // generate_inv_pair emits.
                    for junction in 1..=2 {
                        let offset =
                            balanced_chimeric_offset(frag_len, ctx.config.read_len, &mut rng)?;
                        let result = generate_inv_pair(
                            ctx,
                            contig_name,
                            sv_rec.location,
                            end,
                            junction,
                            frag_len,
                            offset,
                            frag_idx,
                            &mut rng,
                        );

                        match result {
                            Ok((read1, read2)) => {
                                all_reads.push(read1);
                                if let Some(r2) = read2 {
                                    all_reads.push(r2);
                                }
                            }
                            Err(GenerateReadsError::FqToolsError(
                                eidolon_core::file_tools::fastq_tools::FastqToolsError::TruncatedRead(
                                    msg,
                                ),
                            )) => {
                                debug!("Skipping truncated chimeric read: {}", msg);
                            }
                            Err(e) => return Err(e),
                        }
                    }
                }
            } else if sv.sv_type == SvType::Cnv {
                // CNV with INFO/CN: dispatch to the DEL or DUP chimeric path
                // based on whether the total copy number is below or above
                // the diploid baseline. CN < ploidy → loss signature
                // (DEL-like); CN > ploidy → gain signature (DUP-like).
                // CN == ploidy contributes no SV signal and is skipped.
                // A CNV without INFO/CN can't be classified — skip and log.
                let cn = match sv.copy_number {
                    Some(c) => c as usize,
                    None => {
                        debug!(
                            "Skipping <CNV> at {}:{} with no INFO/CN — chimeric path needs CN to choose DEL vs DUP signature",
                            contig_name,
                            sv_rec.location + 1
                        );
                        continue;
                    }
                };
                let ploidy = ctx.config.ploidy;
                if cn == ploidy {
                    continue;
                }

                let end = match sv.end {
                    Some(e) => e,
                    None => {
                        if let Some(span) = sv.span(sv_rec.location) {
                            sv_rec.location + span - 1
                        } else {
                            continue;
                        }
                    }
                };

                let cnv_id = (contig_name.clone(), sv_rec.location, end);
                if !processed_ids.insert(format!("CNV_{:?}", cnv_id)) {
                    continue;
                }

                // For CNVs we coverage-scale by the magnitude of the CN
                // deviation, not by genotype. A CN=0 (full loss) gets
                // 1.0× junction reads (every haplotype carries the
                // junction); CN=1 gets 0.5× at diploid (half the
                // haplotypes); CN=4 gets 2× (one extra tandem copy per
                // ref haplotype, so two tandem junctions). This mirrors
                // how coverage_multiplier_for treats CNVs as cn/ploidy
                // for depth.
                let mult = if ctx.config.subclone_model.is_some() {
                    sv_effective_fraction(sv_rec, ctx.config.ploidy, true)
                } else if cn < ploidy {
                    (ploidy - cn) as f64 / ploidy as f64
                } else {
                    (cn - ploidy) as f64 / ploidy as f64
                };

                let num_frags = scale_coverage(ctx.config.coverage, mult);
                if num_frags == 0 {
                    continue;
                }

                for frag_idx in 0..num_frags {
                    let se_pad = if ctx.config.paired_ended { 0 } else { 32 };
                    let frag_len = if ctx.config.paired_ended {
                        let mut attempts = 0;
                        let mut f = 0;
                        while attempts < 100 {
                            let rand_val = rng.random().map_err(GenerateReadsError::from)?;
                            f = ctx
                                .fragment_length_model
                                .generate_fragment(rand_val)
                                .map_err(GenerateReadsError::from)?
                                as usize;
                            if ctx.config.long_reads || f >= ctx.config.read_len + 10 {
                                break;
                            }
                            attempts += 1;
                        }
                        if f < ctx.config.read_len && !ctx.config.long_reads {
                            f = ctx.config.read_len + 10;
                        }
                        f
                    } else {
                        ctx.config.read_len + se_pad
                    };

                    let offset = balanced_chimeric_offset(frag_len, ctx.config.read_len, &mut rng)?;

                    // CN < ploidy → emit DEL-like junction reads (loss).
                    // CN > ploidy → emit DUP-like junction reads (gain).
                    let result = if cn < ploidy {
                        generate_del_pair(
                            ctx,
                            contig_name,
                            sv_rec.location,
                            end,
                            frag_len,
                            offset,
                            frag_idx,
                            &mut rng,
                        )
                    } else {
                        generate_dup_pair(
                            ctx,
                            contig_name,
                            sv_rec.location,
                            end,
                            frag_len,
                            offset,
                            frag_idx,
                            &mut rng,
                        )
                    };

                    match result {
                        Ok((read1, read2)) => {
                            all_reads.push(read1);
                            if let Some(r2) = read2 {
                                all_reads.push(r2);
                            }
                        }
                        Err(GenerateReadsError::FqToolsError(
                            eidolon_core::file_tools::fastq_tools::FastqToolsError::TruncatedRead(
                                msg,
                            ),
                        )) => {
                            debug!("Skipping truncated chimeric read: {}", msg);
                        }
                        Err(e) => return Err(e),
                    }
                }
            } else if sv.sv_type == SvType::Dup {
                // Symbolic <DUP> generates a tandem duplication: one extra
                // copy of REF[POS..END] inserted immediately after END.
                // The new tandem boundary creates a single junction where
                // the LAST bases of the original dup region butt up
                // against the FIRST bases of the duplicated copy. A read
                // spanning this junction carries left context REF[end-k..end]
                // (end of dup) stitched to right context REF[POS..POS+k]
                // (start of next copy). When BWA aligns these against the
                // unbroken reference the two halves map to "end of dup" and
                // "start of dup" respectively, in the WRONG order — that's
                // the inverted-insert / split-read signature Manta uses to
                // call somatic duplications.
                let end = match sv.end {
                    Some(e) => e,
                    None => {
                        if let Some(span) = sv.span(sv_rec.location) {
                            sv_rec.location + span - 1
                        } else {
                            continue;
                        }
                    }
                };

                let dup_id = (contig_name.clone(), sv_rec.location, end);
                if !processed_ids.insert(format!("DUP_{:?}", dup_id)) {
                    continue;
                }

                let mult = sv_effective_fraction(
                    sv_rec,
                    ctx.config.ploidy,
                    ctx.config.subclone_model.is_some(),
                );

                let num_frags = scale_coverage(ctx.config.coverage, mult);
                if num_frags == 0 {
                    continue;
                }

                for frag_idx in 0..num_frags {
                    let se_pad = if ctx.config.paired_ended { 0 } else { 32 };
                    let frag_len = if ctx.config.paired_ended {
                        let mut attempts = 0;
                        let mut f = 0;
                        while attempts < 100 {
                            let rand_val = rng.random().map_err(GenerateReadsError::from)?;
                            f = ctx
                                .fragment_length_model
                                .generate_fragment(rand_val)
                                .map_err(GenerateReadsError::from)?
                                as usize;
                            if ctx.config.long_reads || f >= ctx.config.read_len + 10 {
                                break;
                            }
                            attempts += 1;
                        }
                        if f < ctx.config.read_len && !ctx.config.long_reads {
                            f = ctx.config.read_len + 10;
                        }
                        f
                    } else {
                        ctx.config.read_len + se_pad
                    };

                    let offset = balanced_chimeric_offset(frag_len, ctx.config.read_len, &mut rng)?;

                    let result = generate_dup_pair(
                        ctx,
                        contig_name,
                        sv_rec.location,
                        end,
                        frag_len,
                        offset,
                        frag_idx,
                        &mut rng,
                    );

                    match result {
                        Ok((read1, read2)) => {
                            all_reads.push(read1);
                            if let Some(r2) = read2 {
                                all_reads.push(r2);
                            }
                        }
                        Err(GenerateReadsError::FqToolsError(
                            eidolon_core::file_tools::fastq_tools::FastqToolsError::TruncatedRead(
                                msg,
                            ),
                        )) => {
                            debug!("Skipping truncated chimeric read: {}", msg);
                        }
                        Err(e) => return Err(e),
                    }
                }
            } else if sv.sv_type == SvType::Del {
                // Symbolic <DEL> generates a single junction at the anchor
                // (POS). Junction reads carry left context REF[..=POS] plus
                // right context REF[END..] (post-deletion). When aligned by
                // BWA against the unbroken reference, these reads produce
                // discordant-PE pairs (mate distance overshoots by END-POS)
                // and split-read alignments (soft-clip at POS realigns to
                // REF[END..]) — the exact signals Manta uses to call somatic
                // DELs. Before #220 the symbolic-DEL path emitted only depth
                // modulation, leaving Manta with 0% recall on the simulated
                // deletions.
                // Same convention as INV: `end` is 1-based VCF END, which
                // numerically equals the 0-based index of the first post-
                // deletion base. The sv.span() / fallback math matches INV's
                // line 961-969 (the off-by-one in passing a 0-based location
                // to sv.span exactly cancels with the `- 1` in the formula).
                let end = match sv.end {
                    Some(e) => e,
                    None => {
                        if let Some(span) = sv.span(sv_rec.location) {
                            sv_rec.location + span - 1
                        } else {
                            continue;
                        }
                    }
                };

                // Use the (contig, location, end) tuple as the dedup key.
                // BND uses pair-canonicalization because both sides describe
                // the same junction; DEL has a single record per deletion,
                // so the same-position dedup is enough.
                let del_id = (contig_name.clone(), sv_rec.location, end);
                if !processed_ids.insert(format!("DEL_{:?}", del_id)) {
                    continue;
                }

                // Same coverage model as BND/INV — full coverage for
                // homozygous (every allele carries the junction), 1/ploidy
                // for heterozygous. The double-counting caveat noted on
                // the BND branch applies here too: the regular per-contig
                // pass still generates reads spanning POS (it reads from
                // the unbroken reference). For homozygous DELs the
                // breakpoint locus ends up with regular reads PLUS
                // junction reads. Tracked at #220's "reconcile depth-
                // modulation" follow-up.
                let mult = match sv_rec.genotype {
                    Genotype::Homozygous => 1.0,
                    Genotype::Heterozygous => 1.0 / (ctx.config.ploidy as f64),
                };

                let num_frags = scale_coverage(ctx.config.coverage, mult);
                if num_frags == 0 {
                    continue;
                }

                for frag_idx in 0..num_frags {
                    let se_pad = if ctx.config.paired_ended { 0 } else { 32 };
                    let frag_len = if ctx.config.paired_ended {
                        let mut attempts = 0;
                        let mut f = 0;
                        while attempts < 100 {
                            let rand_val = rng.random().map_err(GenerateReadsError::from)?;
                            f = ctx
                                .fragment_length_model
                                .generate_fragment(rand_val)
                                .map_err(GenerateReadsError::from)?
                                as usize;
                            if ctx.config.long_reads || f >= ctx.config.read_len + 10 {
                                break;
                            }
                            attempts += 1;
                        }
                        if f < ctx.config.read_len && !ctx.config.long_reads {
                            f = ctx.config.read_len + 10;
                        }
                        f
                    } else {
                        ctx.config.read_len + se_pad
                    };

                    let offset = balanced_chimeric_offset(frag_len, ctx.config.read_len, &mut rng)?;

                    let result = generate_del_pair(
                        ctx,
                        contig_name,
                        sv_rec.location,
                        end,
                        frag_len,
                        offset,
                        frag_idx,
                        &mut rng,
                    );

                    match result {
                        Ok((read1, read2)) => {
                            all_reads.push(read1);
                            if let Some(r2) = read2 {
                                all_reads.push(r2);
                            }
                        }
                        Err(GenerateReadsError::FqToolsError(
                            eidolon_core::file_tools::fastq_tools::FastqToolsError::TruncatedRead(
                                msg,
                            ),
                        )) => {
                            debug!("Skipping truncated chimeric read: {}", msg);
                        }
                        Err(e) => return Err(e),
                    }
                }
            }
        }
    }

    // Write all chimeric reads to a temp fastq and BAM
    let contig_name = "chimeric".to_string();
    let idx = 999999;

    let mut contig_files_r1 = Vec::new();
    let mut contig_files_r2 = Vec::new();
    let mut bam_body_file = None;

    if !all_reads.is_empty() {
        if ctx.config.produce_fastq {
            let mut file_to_write_1 = PathBuf::from(ctx.working_dir);
            file_to_write_1.push("temp_chimeric_r1.fastq.gz");
            let file1 = append_to_file(&file_to_write_1)?;
            let writer1 = BufWriter::new(&file1);
            let mut buffer1 = GzEncoder::new(writer1, Compression::default());

            if ctx.config.paired_ended {
                let mut file_to_write_2 = PathBuf::from(ctx.working_dir);
                file_to_write_2.push("temp_chimeric_r2.fastq.gz");
                let file2 = append_to_file(&file_to_write_2)?;
                let writer2 = BufWriter::new(&file2);
                let mut buffer2 = GzEncoder::new(writer2, Compression::default());

                for i in (0..all_reads.len()).step_by(2) {
                    write_read_to_fastq(&all_reads[i], &mut buffer1)
                        .map_err(GenerateReadsError::from)?;
                    if i + 1 < all_reads.len() {
                        write_read_to_fastq(&all_reads[i + 1], &mut buffer2)
                            .map_err(GenerateReadsError::from)?;
                    }
                }
                contig_files_r1.push(file_to_write_1);
                contig_files_r2.push(file_to_write_2);
            } else {
                for read in &all_reads {
                    write_read_to_fastq(read, &mut buffer1).map_err(GenerateReadsError::from)?;
                }
                contig_files_r1.push(file_to_write_1);
            }
        }

        if let Some(bam_ctx) = &ctx.bam_context {
            let bam_temp_path = PathBuf::from(ctx.working_dir).join("temp_bam_chimeric.bam");
            let mut writer = BamBodyWriter::new(bam_temp_path.clone(), Arc::clone(bam_ctx))?;
            for read in &all_reads {
                writer.stage_read_record(read).map_err(|e| {
                    GenerateReadsError::IoError(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        e.to_string(),
                    ))
                })?;
            }
            // BamBodyWriter::stage_read_record buffers for deferred
            // coordinate-sorted output. flush_all drains the buffer; without
            // it, the staged chimeric reads never reach disk. The regular
            // per-contig path calls flush_all at process_contig's tail; the
            // chimeric path needs the same.
            writer.flush_all()?;
            bam_body_file = Some(bam_temp_path);
        }
    }

    Ok(ChunkResult {
        contig_idx: idx,
        chunk_start: 0,
        name: contig_name,
        len: 0,
        data: Some(ChunkData {
            r1_files: contig_files_r1,
            r2_files: contig_files_r2,
            bam_body_file,
            ad_counter: AdCounter::new(),
        }),
    })
}

fn generate_chimeric_pair(
    ctx: &ContigContext,
    contig: &str,
    pos: usize,
    sv: &SvData,
    frag_len: usize,
    offset: usize,
    frag_idx: usize,
    rng: &mut NeatRng,
) -> Result<(ReadRecord, Option<ReadRecord>), GenerateReadsError> {
    // Offset is where the junction is relative to the start of the fragment.
    // frag_len = L1 + L2
    // L1 = offset
    // L2 = frag_len - offset

    let ((c1, s1, e1, rev1), (c2, s2, e2, rev2)) =
        get_bnd_pieces(contig, pos, sv, offset, frag_len - offset, &ctx.reference)?;

    let seq1 = get_stitched_sequence(ctx, &c1, s1, e1, rev1, &c2, s2, e2, rev2, rng)?;

    let read_len = ctx.config.read_len;
    // The trailing 16-hex `frag_idx` matches the uniqueness-tag pattern that
    // write_block_fastq uses for regular reads (#210). Without it, two
    // chimeric reads spawned from the same BND (num_frags > 1) would share
    // a QNAME and Picard MarkDuplicates would drop one as a "PCR duplicate".
    let base_name = format!(
        "EIDOLON_chimeric_{}_{}_{}_{}_{:016x}",
        c1,
        pos,
        c2,
        sv.mate_pos.unwrap_or(0),
        frag_idx,
    );

    let quality_scores_1 = ctx
        .quality_score_model
        .generate_quality_scores(read_len, rng)
        .map_err(GenerateReadsError::from)?;

    // Chimeric pairs don't drive AD/DP/AF (BND junction reads are span/
    // discordant-pair signal, not point-coverage signal). Pass a throwaway
    // local AdCounter so generate_read still increments somewhere but the
    // values don't leak into the per-contig counter used by write_vcf.
    let mut throwaway_ad = AdCounter::new();
    let r1 = generate_read(
        &seq1,
        // BND junction reads are stitched from reference pieces: every base is 'M', and
        // no haplotype deletion applies.
        None,
        None,
        &[], // Mutations already applied in get_stitched_sequence
        &HashMap::new(),
        read_len,
        format!("{}/1", base_name),
        Strand::Forward,
        quality_scores_1,
        ctx.seq_error_model,
        rng,
        c1.clone(),
        s1,
        c2.clone(),
        s2,
        frag_len as i32,
        ctx.config.paired_ended,
        &mut throwaway_ad,
    )
    .map_err(GenerateReadsError::from)?;

    let mut r2 = None;
    if ctx.config.paired_ended {
        let quality_scores_2 = ctx
            .quality_score_model
            .generate_quality_scores(read_len, rng)
            .map_err(GenerateReadsError::from)?;
        let r2_record = generate_read(
            &reverse_complement(seq1),
            // Reference-derived bases only: no haplotype mask, no haplotype deletion.
            None,
            None,
            &[],
            &HashMap::new(),
            read_len,
            format!("{}/2", base_name),
            Strand::Reverse,
            quality_scores_2,
            ctx.seq_error_model,
            rng,
            c2.clone(),
            s2,
            c1.clone(),
            s1,
            -(frag_len as i32),
            true,
            &mut throwaway_ad,
        )
        .map_err(GenerateReadsError::from)?;
        r2 = Some(r2_record);
    }

    Ok((r1, r2))
}

fn generate_inv_pair(
    ctx: &ContigContext,
    contig: &str,
    location: usize, // 0-based
    end: usize,      // 1-based
    junction: usize,
    frag_len: usize,
    offset: usize,
    frag_idx: usize,
    rng: &mut NeatRng,
) -> Result<(ReadRecord, Option<ReadRecord>), GenerateReadsError> {
    let ((c1, s1, e1, rev1), (c2, s2, e2, rev2)) = get_inv_pieces(
        contig,
        location,
        end,
        junction,
        offset,
        frag_len - offset,
        ctx,
    )?;

    let seq1 = get_stitched_sequence(ctx, &c1, s1, e1, rev1, &c2, s2, e2, rev2, rng)?;

    let read_len = ctx.config.read_len;
    // The trailing 16-hex `frag_idx` matches the uniqueness-tag pattern that
    // write_block_fastq + generate_chimeric_pair use (#210). The `junction`
    // tag already disambiguates the two breakpoints of a single inversion;
    // frag_idx disambiguates fragments at the same breakpoint when
    // num_frags > 1.
    let base_name = format!(
        "EIDOLON_chimeric_INV_{}_{}_{}_{}_{:016x}",
        contig,
        location + 1,
        end,
        junction,
        frag_idx,
    );

    let quality_scores_1 = ctx
        .quality_score_model
        .generate_quality_scores(read_len, rng)
        .map_err(GenerateReadsError::from)?;

    // Inversion-junction reads contribute to junction signal, not point
    // coverage — same rationale as generate_chimeric_pair. Use a throwaway
    // AdCounter so generate_read's increment site has somewhere to write
    // without leaking values into the per-contig counter used by write_vcf.
    let mut throwaway_ad = AdCounter::new();

    let r1 = generate_read(
        &seq1,
        // Reference-derived bases only: no haplotype insertion mask, no haplotype
        // deletion.
        None,
        None,
        &[],
        &HashMap::new(),
        read_len,
        format!("{}/1", base_name),
        Strand::Forward,
        quality_scores_1,
        ctx.seq_error_model,
        rng,
        c1.clone(),
        s1,
        c2.clone(),
        s2,
        frag_len as i32,
        ctx.config.paired_ended,
        &mut throwaway_ad,
    )
    .map_err(GenerateReadsError::from)?;

    let mut r2 = None;
    if ctx.config.paired_ended {
        let quality_scores_2 = ctx
            .quality_score_model
            .generate_quality_scores(read_len, rng)
            .map_err(GenerateReadsError::from)?;
        let r2_record = generate_read(
            &reverse_complement(seq1),
            // Reference-derived bases only: no haplotype mask, no haplotype deletion.
            None,
            None,
            &[],
            &HashMap::new(),
            read_len,
            format!("{}/2", base_name),
            Strand::Reverse,
            quality_scores_2,
            ctx.seq_error_model,
            rng,
            c2.clone(),
            s2,
            c1.clone(),
            s1,
            -(frag_len as i32),
            true,
            &mut throwaway_ad,
        )
        .map_err(GenerateReadsError::from)?;
        r2 = Some(r2_record);
    }

    Ok((r1, r2))
}

/// Generate a chimeric junction read-pair for a symbolic <DEL>. The
/// deletion creates a single junction at the anchor: left context
/// REF[..=location] stitched to right context REF[end..]. When BWA
/// aligns these reads back to the unbroken reference, the contiguous
/// stitched fragment surfaces as either a discordant PE pair (mate
/// distance overshoots by end-location) or a split-read alignment
/// (soft-clip at the junction realigns to REF[end..]). Both are the
/// signals Manta consumes to call somatic deletions.
fn generate_del_pair(
    ctx: &ContigContext,
    contig: &str,
    location: usize, // 0-based location (POS-1 = anchor index)
    end: usize,      // 1-based VCF END (= 0-based start of post-DEL right piece)
    frag_len: usize,
    offset: usize,
    frag_idx: usize,
    rng: &mut NeatRng,
) -> Result<(ReadRecord, Option<ReadRecord>), GenerateReadsError> {
    let ((c1, s1, e1, rev1), (c2, s2, e2, rev2)) =
        get_del_pieces(contig, location, end, offset, frag_len - offset, ctx)?;

    let seq1 = get_stitched_sequence(ctx, &c1, s1, e1, rev1, &c2, s2, e2, rev2, rng)?;

    let read_len = ctx.config.read_len;
    // QNAME format mirrors BND/INV's `EIDOLON_chimeric_*` scheme so the
    // FASTQ-validation tests in fastq_validation.rs and the read-name
    // parser in filter_lib.rs handle them uniformly. The `DEL` tag
    // disambiguates from BND (`EIDOLON_chimeric_<c1>_<pos>_<c2>_...`)
    // and INV (`EIDOLON_chimeric_INV_<contig>_<pos>_<end>_<junction>_...`).
    let base_name = format!(
        "EIDOLON_chimeric_DEL_{}_{}_{}_{:016x}",
        contig,
        location + 1,
        end,
        frag_idx,
    );

    let quality_scores_1 = ctx
        .quality_score_model
        .generate_quality_scores(read_len, rng)
        .map_err(GenerateReadsError::from)?;

    // DEL junction reads are span/discordant-pair signal, not point-
    // coverage signal — same rationale as generate_chimeric_pair. Use a
    // throwaway AdCounter so the increment site has somewhere to write
    // without leaking into the per-contig counter used by write_vcf.
    let mut throwaway_ad = AdCounter::new();

    let r1 = generate_read(
        &seq1,
        // Reference-derived bases only: no haplotype insertion mask, no haplotype
        // deletion.
        None,
        None,
        &[],
        &HashMap::new(),
        read_len,
        format!("{}/1", base_name),
        Strand::Forward,
        quality_scores_1,
        ctx.seq_error_model,
        rng,
        c1.clone(),
        s1,
        c2.clone(),
        s2,
        frag_len as i32,
        ctx.config.paired_ended,
        &mut throwaway_ad,
    )
    .map_err(GenerateReadsError::from)?;

    let mut r2 = None;
    if ctx.config.paired_ended {
        let quality_scores_2 = ctx
            .quality_score_model
            .generate_quality_scores(read_len, rng)
            .map_err(GenerateReadsError::from)?;
        let r2_record = generate_read(
            &reverse_complement(seq1),
            // Reference-derived bases only: no haplotype mask, no haplotype deletion.
            None,
            None,
            &[],
            &HashMap::new(),
            read_len,
            format!("{}/2", base_name),
            Strand::Reverse,
            quality_scores_2,
            ctx.seq_error_model,
            rng,
            c2.clone(),
            s2,
            c1.clone(),
            s1,
            -(frag_len as i32),
            true,
            &mut throwaway_ad,
        )
        .map_err(GenerateReadsError::from)?;
        r2 = Some(r2_record);
    }

    Ok((r1, r2))
}

/// The two pieces that make a DEL junction. Left piece is REF[..=location]
/// (anchor included, mirroring BND case 1's `t[p[` anchor-preserving
/// shape). Right piece is REF[end..] — the first post-deletion base.
/// Both pieces are forward-strand on the same contig.
fn get_del_pieces(
    contig: &str,
    location: usize, // 0-based location (POS-1 = anchor index)
    end: usize, // 1-based VCF END (= 0-based exclusive end of deletion = 0-based start of post-DEL region)
    len1: usize,
    len2: usize,
    ctx: &ContigContext,
) -> Result<((String, usize, usize, bool), (String, usize, usize, bool)), GenerateReadsError> {
    let c_len = ctx.reference.get(contig).map(|s| s.len()).ok_or_else(|| {
        GenerateReadsError::CliError(format!(
            "DEL at {contig}:{location} references contig {contig} but that contig is not in the reference"
        ))
    })?;
    // Left piece ends at index location+1 (exclusive), so the last base is
    // REF[location] — the anchor base. s1 = e1 - len1 (or 0 if it'd go
    // negative).
    let e1 = (location + 1).min(c_len);
    let s1 = e1.saturating_sub(len1);
    // Right piece starts at 0-based index `end` (the first post-deletion
    // base) and extends len2 bases forward (capped at contig length).
    let s2 = end.min(c_len);
    let e2 = (s2 + len2).min(c_len);
    Ok((
        (contig.to_string(), s1, e1, false),
        (contig.to_string(), s2, e2, false),
    ))
}

/// Generate a chimeric junction read-pair for a symbolic <DUP> (tandem).
/// The duplication creates one extra copy of REF[POS..END]; the tandem
/// boundary is between the last base of the original copy and the first
/// base of the new copy. Junction reads carry left context REF[end-k..end]
/// stitched to right context REF[location..location+k]. When BWA aligns
/// them against the unbroken reference the two halves map to "end of
/// dup" and "start of dup" — in the inverted order — surfacing as
/// split-read or wrong-orientation PE signal.
fn generate_dup_pair(
    ctx: &ContigContext,
    contig: &str,
    location: usize, // 0-based POS (first base of dup region)
    end: usize,      // 1-based END (= 0-based exclusive end of dup region)
    frag_len: usize,
    offset: usize,
    frag_idx: usize,
    rng: &mut NeatRng,
) -> Result<(ReadRecord, Option<ReadRecord>), GenerateReadsError> {
    let ((c1, s1, e1, rev1), (c2, s2, e2, rev2)) =
        get_dup_pieces(contig, location, end, offset, frag_len - offset, ctx)?;

    let seq1 = get_stitched_sequence(ctx, &c1, s1, e1, rev1, &c2, s2, e2, rev2, rng)?;

    let read_len = ctx.config.read_len;
    let base_name = format!(
        "EIDOLON_chimeric_DUP_{}_{}_{}_{:016x}",
        contig,
        location + 1,
        end,
        frag_idx,
    );

    let quality_scores_1 = ctx
        .quality_score_model
        .generate_quality_scores(read_len, rng)
        .map_err(GenerateReadsError::from)?;
    let mut throwaway_ad = AdCounter::new();

    let r1 = generate_read(
        &seq1,
        // Reference-derived bases only: no haplotype insertion mask, no haplotype
        // deletion.
        None,
        None,
        &[],
        &HashMap::new(),
        read_len,
        format!("{}/1", base_name),
        Strand::Forward,
        quality_scores_1,
        ctx.seq_error_model,
        rng,
        c1.clone(),
        s1,
        c2.clone(),
        s2,
        frag_len as i32,
        ctx.config.paired_ended,
        &mut throwaway_ad,
    )
    .map_err(GenerateReadsError::from)?;

    let mut r2 = None;
    if ctx.config.paired_ended {
        let quality_scores_2 = ctx
            .quality_score_model
            .generate_quality_scores(read_len, rng)
            .map_err(GenerateReadsError::from)?;
        let r2_record = generate_read(
            &reverse_complement(seq1),
            // Reference-derived bases only: no haplotype mask, no haplotype deletion.
            None,
            None,
            &[],
            &HashMap::new(),
            read_len,
            format!("{}/2", base_name),
            Strand::Reverse,
            quality_scores_2,
            ctx.seq_error_model,
            rng,
            c2.clone(),
            s2,
            c1.clone(),
            s1,
            -(frag_len as i32),
            true,
            &mut throwaway_ad,
        )
        .map_err(GenerateReadsError::from)?;
        r2 = Some(r2_record);
    }

    Ok((r1, r2))
}

/// The two pieces that make a tandem-DUP junction. Left piece is the
/// END of the duplicated region (REF[end-len1..end]); right piece is the
/// START of the duplicated region (REF[location..location+len2]) —
/// representing the first bases of the new tandem copy.
///
/// For small DUPs (span smaller than len1 or len2) the pieces are
/// capped by the contig boundaries only; the resulting stitched
/// fragment may then look like the unbroken reference (if the windows
/// extend into surrounding context on both sides), in which case BWA
/// aligns it normally and the read contributes no SV signal. That's
/// acceptable — tiny DUPs are below most callers' resolution anyway.
fn get_dup_pieces(
    contig: &str,
    location: usize, // 0-based POS (first base of dup region)
    end: usize,      // 1-based END (= 0-based exclusive end of dup region)
    len1: usize,
    len2: usize,
    ctx: &ContigContext,
) -> Result<((String, usize, usize, bool), (String, usize, usize, bool)), GenerateReadsError> {
    let c_len = ctx.reference.get(contig).map(|s| s.len()).ok_or_else(|| {
        GenerateReadsError::CliError(format!(
            "DUP at {contig}:{location} references contig {contig} but that contig is not in the reference"
        ))
    })?;
    // Left piece: last len1 bases of the duplicated region (or earlier
    // if the dup is shorter than len1, in which case we extend into the
    // pre-dup context — see fn doc for the implication).
    let e1 = end.min(c_len);
    let s1 = e1.saturating_sub(len1);
    // Right piece: first len2 bases of the duplicated region. POS is the base
    // BEFORE the event (VCF 4.2), so the region begins at location+1.
    let s2 = location + 1;
    let e2 = (s2 + len2).min(c_len);
    Ok((
        (contig.to_string(), s1, e1, false),
        (contig.to_string(), s2, e2, false),
    ))
}

fn get_inv_pieces(
    contig: &str,
    location: usize, // 0-based location (POS-1)
    end: usize,      // 1-based END
    junction: usize, // 1 or 2
    len1: usize,
    len2: usize,
    ctx: &ContigContext,
) -> Result<((String, usize, usize, bool), (String, usize, usize, bool)), GenerateReadsError> {
    // Error rather than silently defaulting to zero-length sequences if the
    // inversion's contig is missing from the reference — same defense-in-
    // depth as get_bnd_pieces.
    let c_len = ctx.reference.get(contig).map(|s| s.len()).ok_or_else(|| {
        GenerateReadsError::CliError(format!(
            "INV at {contig}:{location} references contig {contig} but that contig is not in the reference"
        ))
    })?;
    Ok(if junction == 1 {
        // Junction 1: REF[..POS] | RC(REF[POS+1..END])
        // POS is the base BEFORE the inverted block (VCF 4.2), so the block begins
        // at location+1 and the left piece ends there.
        let e1 = location + 1;
        let s1 = e1.saturating_sub(len1);

        let e2 = end.min(c_len);
        let s2 = e2.saturating_sub(len2).max(location + 1);
        (
            (contig.to_string(), s1, e1, false),
            (contig.to_string(), s2, e2, true),
        )
    } else {
        // Junction 2: RC(REF[POS+1..END]) | REF[END+1..]
        // Same convention as junction 1: the inverted block begins at location+1.
        let s1 = location + 1;
        let e1 = (s1 + len1).min(end).min(c_len);

        let s2 = end.min(c_len);
        let e2 = (s2 + len2).min(c_len);
        (
            (contig.to_string(), s1, e1, true),
            (contig.to_string(), s2, e2, false),
        )
    })
}

/// Split a BND into the two reference pieces a chimeric fragment is stitched from.
///
/// Takes the reference map rather than the whole `ContigContext` because contig lengths
/// are all it needs — and because the geometry it selects is the thing most worth
/// testing directly. The `bool` in each tuple is "reverse-complement this piece".
fn get_bnd_pieces(
    contig: &str,
    pos: usize, // 0-based
    sv: &SvData,
    len1: usize,
    len2: usize,
    reference: &HashMap<String, Vec<Nucleotide>>,
) -> Result<((String, usize, usize, bool), (String, usize, usize, bool)), GenerateReadsError> {
    let mate_contig = sv.mate_contig.as_ref().unwrap().clone();
    let mate_pos = sv.mate_pos.unwrap().saturating_sub(1);

    // BNDs can legitimately point at a contig outside the reference (a real
    // VCF data quality issue). Surface that as an error rather than silently
    // producing zero-length sequences via `unwrap_or(0)`.
    let c1_len = reference.get(contig).map(|s| s.len()).ok_or_else(|| {
        GenerateReadsError::CliError(format!(
            "BND at {contig}:{pos} references its own contig {contig} but that contig is not in the reference"
        ))
    })?;
    let c2_len = reference.get(&mate_contig).map(|s| s.len()).ok_or_else(|| {
        GenerateReadsError::CliError(format!(
            "BND at {contig}:{pos} has mate on contig {mate_contig} but that contig is not in the reference"
        ))
    })?;

    Ok(if sv.bnd_join_after {
        if sv.bnd_mate_extends_right {
            // Case 1: t[p[ -> REF[..=pos] + MATE[mate_pos..]
            let s1 = pos.saturating_sub(len1.saturating_sub(1));
            let e1 = pos + 1;
            let s2 = mate_pos;
            let e2 = (mate_pos + len2).min(c2_len);
            (
                (contig.to_string(), s1, e1, false),
                (mate_contig, s2, e2, false),
            )
        } else {
            // Case 2: t]p] -> REF[..=pos] + revcomp(MATE[..=mate_pos])
            let s1 = pos.saturating_sub(len1.saturating_sub(1));
            let e1 = pos + 1;
            let e2 = mate_pos + 1;
            let s2 = e2.saturating_sub(len2);
            (
                (contig.to_string(), s1, e1, false),
                (mate_contig, s2, e2, true),
            )
        }
    } else if sv.bnd_mate_extends_right {
        // Case 3: [p[t -> revcomp(MATE[mate_pos..]) + REF[pos..]
        let s1 = mate_pos;
        let e1 = (mate_pos + len1).min(c2_len);
        let s2 = pos;
        let e2 = (pos + len2).min(c1_len);
        (
            (mate_contig, s1, e1, true),
            (contig.to_string(), s2, e2, false),
        )
    } else {
        // Case 4: ]p]t -> MATE[..=mate_pos] + REF[pos..]
        let e1 = mate_pos + 1;
        let s1 = e1.saturating_sub(len1);
        let s2 = pos;
        let e2 = (pos + len2).min(c1_len);
        (
            (mate_contig, s1, e1, false),
            (contig.to_string(), s2, e2, false),
        )
    })
}

fn get_stitched_sequence(
    ctx: &ContigContext,
    c1: &str,
    s1: usize,
    e1: usize,
    rev1: bool,
    c2: &str,
    s2: usize,
    e2: usize,
    rev2: bool,
    rng: &mut NeatRng,
) -> Result<Vec<Nucleotide>, GenerateReadsError> {
    let mut seq1 = get_mutated_subseq(ctx, c1, s1, e1, rng)?;
    if rev1 {
        seq1 = reverse_complement(seq1);
    }
    let mut seq2 = get_mutated_subseq(ctx, c2, s2, e2, rng)?;
    if rev2 {
        seq2 = reverse_complement(seq2);
    }
    seq1.extend(seq2);
    Ok(seq1)
}

fn get_mutated_subseq(
    ctx: &ContigContext,
    contig: &str,
    start: usize,
    end: usize,
    rng: &mut NeatRng,
) -> Result<Vec<Nucleotide>, GenerateReadsError> {
    let ref_seq = ctx.reference.get(contig).ok_or_else(|| {
        GenerateReadsError::IoError(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("Contig {} not found", contig),
        ))
    })?;
    let m_map = ctx.mutated_maps.get(contig).ok_or_else(|| {
        GenerateReadsError::IoError(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("MutatedMap for {} not found", contig),
        ))
    })?;

    let mut seq = Vec::new();
    for i in start..end {
        if i >= ref_seq.len() {
            seq.push(Nucleotide::N);
            continue;
        }
        if m_map.is_flagged(&i) {
            let variants = m_map
                .mutate_position(i, rng)
                .map_err(GenerateReadsError::from)?;
            seq.extend(variants);
        } else {
            seq.push(ref_seq[i]);
        }
    }
    Ok(seq)
}

fn collect_chunk_result(
    cr: ChunkResult,
    contig_order: &mut Vec<String>,
    fasta_lengths: &mut HashMap<String, usize>,
    all_fastq_files: &mut HashMap<String, (Vec<PathBuf>, Vec<PathBuf>)>,
    bam_body_files: &mut HashMap<String, Vec<(usize, PathBuf)>>,
    ad_counters: &mut HashMap<String, AdCounter>,
) {
    // `process_chimeric_variants` returns a pseudo ChunkResult named
    // "chimeric" with len=0 — a control-flow tag for organizing the
    // junction-spanning FASTQ/BAM outputs, NOT a real reference contig.
    // Its reads' positions are already keyed to real contigs (the BND
    // mate contig + position), so we only want this pseudo-contig's
    // FASTQ/BAM outputs to flow into the final merge — never its name
    // into VCF-relevant accumulators. Leaving it in `contig_order` /
    // `fasta_lengths` produces a malformed `##contig=<ID=chimeric,length=0>`
    // header line that breaks strict downstream parsers like truvari.
    //
    // Results arrive sorted by (contig_idx, chunk_start), so the first chunk of
    // each contig records its order/length; later chunks of the same contig
    // accumulate their reads, AD counts, and BAM bodies.
    let is_chimeric_pseudo = cr.name == "chimeric" && cr.len == 0;
    if !is_chimeric_pseudo && !fasta_lengths.contains_key(&cr.name) {
        contig_order.push(cr.name.clone());
        fasta_lengths.insert(cr.name.clone(), cr.len);
    }
    if let Some(data) = cr.data {
        if !is_chimeric_pseudo {
            // Merge per-variant allelic depth across this contig's chunks.
            let acc = ad_counters.entry(cr.name.clone()).or_default();
            for (key, (ref_n, alt_n)) in data.ad_counter {
                let slot = acc.entry(key).or_insert((0, 0));
                slot.0 += ref_n;
                slot.1 += alt_n;
            }
        }
        // FASTQ outputs (order-independent) accumulate per contig; chimeric
        // reads carry real per-real-contig records and still flow into the merge.
        let entry = all_fastq_files.entry(cr.name.clone()).or_default();
        entry.0.extend(data.r1_files);
        entry.1.extend(data.r2_files);
        if let Some(bam_path) = data.bam_body_file {
            bam_body_files
                .entry(cr.name)
                .or_default()
                .push((cr.chunk_start, bam_path));
        }
    }
}

/// Removes literal multi-base Complex variants from the input map, emitting a
/// warning for each one. SNPs, insertions, and deletions are kept as-is.
/// Symbolic / structural variants (`<DEL>`, `<DUP>`, `<CNV>`, ...) are also
/// kept — gen_reads uses them downstream to modulate coverage and to round-
/// trip into the output VCF; they never go through per-base mutation.
/// Append a `KEY=value` tag to a variant's optional INFO string, merging with any
/// existing content (`;`-joined). Used to attach simulator ground-truth tags
/// (EIDOLON_CCF, EIDOLON_VAF) to somatic variants for the golden VCF (#405).
fn append_info_tag(info: &mut Option<String>, tag: String) {
    *info = Some(match info.take() {
        Some(e) if !e.is_empty() && e != "." => format!("{e};{tag}"),
        _ => tag,
    });
}

/// Split a GT into its allele tokens, handling both `/` and `|` separators.
fn gt_alleles(gt: &str) -> Vec<&str> {
    if gt.contains('/') {
        gt.split('/').collect()
    } else {
        gt.split('|').collect()
    }
}

/// True when a GT calls the reference on every allele (`0/0`, `0|0`, `0`).
///
/// Such a record says this sample does NOT carry the alt — it is common in cohort VCFs,
/// where a site is present because some other sample has it. `gt_from_str` maps it to
/// `Heterozygous` (it sees a zero allele and stops), so the alt was generated at ~0.5
/// while the truth VCF carried `GT=0/0` beside `AF=0.4857` (#591).
fn is_hom_ref(gt: &str) -> bool {
    let alleles = gt_alleles(gt);
    !alleles.is_empty() && alleles.iter().all(|a| a.trim() == "0")
}

/// Rewrite each no-call allele in a PARTIAL no-call as a reference allele: `./1` -> `0/1`,
/// `././1` -> `0/0/1`. Returns `None` when nothing changed.
///
/// Two reasons, and the downstream one decides it. **Semantically**, `./1` tells us about
/// exactly one alternate allele; generating `1/1` invents a second we were never given,
/// while `0/1` asserts precisely what the input said and nothing more.
///
/// **Downstream**, `bcftools` classifies a partial no-call as neither het nor hom —
/// measured, `GT="het"` selects only `0/1` and `GT="hom"` only `1/1`, while `./1` appears
/// solely under `GT="mis"`. A `./1` truth record therefore drops out of *any*
/// zygosity-stratified analysis without a word. `0/1` is understood everywhere.
///
/// This also fixes the allele fraction for free: `dosage_fraction` skips `.` alleles, so
/// `./1` scored alt=1 of total=1 = 1.0 and the variant was generated as fully homozygous
/// (#591). `0/1` scores 1 of 2 = 0.5 with no special case, and generalises to any ploidy.
fn normalize_partial_no_call(gt: &str) -> Option<String> {
    if !gt.contains('.') {
        return None;
    }
    let sep = if gt.contains('/') { '/' } else { '|' };
    let out: Vec<String> = gt
        .split(sep)
        .map(|a| {
            if a.trim() == "." {
                "0".to_string()
            } else {
                a.trim().to_string()
            }
        })
        .collect();
    let joined = out.join(&sep.to_string());
    if joined == gt { None } else { Some(joined) }
}

/// True when every allele of a GT is `.` — a no-call.
///
/// `gt_from_str` cannot express this: it skips `.` alleles, never sets `found_zero`, and
/// falls through to `Homozygous`. So a `./.` record silently becomes homozygous ALT and
/// the truth VCF then carries `GT=./.` next to `AF=1.0000` — a record contradicting
/// itself (#591). Detect it on the raw string, before that lossy conversion.
fn is_no_call(gt: &str) -> bool {
    let alleles = gt_alleles(gt);
    !alleles.is_empty() && alleles.iter().all(|a| a.trim() == ".")
}

fn filter_input_vcf(raw: HashMap<String, Vec<Variant>>) -> HashMap<String, Vec<Variant>> {
    let mut out: HashMap<String, Vec<Variant>> = HashMap::new();
    let (mut n_complex, mut n_null, mut n_nocall, mut n_homref) = (0usize, 0usize, 0usize, 0usize);
    let mut n_partial = 0usize;
    for (contig, variants) in raw {
        let mut kept = Vec::new();
        for v in variants {
            // Order matters only for which warning a doubly-bad record gets; each arm drops.
            if v.variant_type == VariantType::Complex && v.alternate.is_literal() {
                n_complex += 1;
                warn!(
                    "Skipping complex variant at {}:{} (multi-base REF and ALT that is not \
                     a simple indel) — not yet supported",
                    contig, v.location
                );
            } else if v.alternate.as_literal() == Some(v.reference.as_slice()) {
                // REF == ALT is not a variant: there is no alternate allele to carry. Keeping
                // it would put an unachievable record in the truth VCF, where it becomes a
                // false negative for every caller and inflates the recall denominator (#591).
                n_null += 1;
                warn!(
                    "Skipping non-variant at {}:{} (ALT is identical to REF, so there is no \
                     alternate allele to generate)",
                    contig, v.location
                );
            } else if is_hom_ref(&v.genotype_str) {
                // 0/0 says this sample carries the reference on every allele. Generating the
                // ALT anyway produced AF~0.49 beside GT=0/0 (#591). Common in cohort VCFs,
                // where the site exists because a DIFFERENT sample carries it.
                n_homref += 1;
                warn!(
                    "Skipping homozygous-reference record at {}:{} (GT={}) — this sample \
                     does not carry the alternate allele",
                    contig, v.location, v.genotype_str
                );
            } else if is_no_call(&v.genotype_str) {
                // A no-call says nothing about what to generate. Applying it as homozygous
                // ALT — which is what the genotype parser does by default — asserts a
                // variant the GT explicitly declines to call (#591).
                n_nocall += 1;
                warn!(
                    "Skipping no-call genotype at {}:{} (GT={}) — a no-call carries no \
                     instruction about what to generate",
                    contig, v.location, v.genotype_str
                );
            } else {
                // A PARTIAL no-call is a real call on at least one allele, so it is kept —
                // but the `.` is rewritten to `0` first. Done after the full-no-call and
                // hom-ref arms above so those keep their own, clearer warnings.
                let mut v = v;
                if let Some(fixed) = normalize_partial_no_call(&v.genotype_str) {
                    n_partial += 1;
                    debug!(
                        "Normalizing partial no-call at {}:{}: GT={} -> {} (an uncalled \
                         allele is generated as reference, not as a second alt)",
                        contig, v.location, v.genotype_str, fixed
                    );
                    v.genotype_str = fixed;
                    // Normalization replaced at least one `.` with `0`, so the result
                    // carries a reference allele by construction and the coarse label is
                    // Heterozygous. The precise per-copy dosage comes from
                    // `dosage_fraction()` reading the string, which now counts that `0`.
                    v.genotype = Genotype::Heterozygous;
                }
                kept.push(v);
            }
        }
        if !kept.is_empty() {
            out.insert(contig, kept);
        }
    }
    // A step that drops input deliberately has to say how much, or a shrunken truth set
    // looks like a small one.
    if n_partial > 0 {
        info!(
            "input_vcf: normalized {n_partial} partial no-call genotype(s) — an uncalled \
             allele is generated as reference (./1 -> 0/1)"
        );
    }
    let dropped = n_complex + n_null + n_nocall + n_homref;
    if dropped > 0 {
        warn!(
            "input_vcf: dropped {dropped} record(s) — {n_complex} complex, {n_null} \
             ALT-equals-REF, {n_nocall} no-call genotype, {n_homref} homozygous-reference"
        );
    }
    out
}

/// Overrides the mutation rate in [ovr_start, ovr_end) within existing segments,
/// splitting segment boundaries where needed. Positions outside existing segments
/// (N-regions, gaps) are not affected.
fn apply_rate_override(
    segments: Vec<(usize, usize, f64)>,
    ovr_start: usize,
    ovr_end: usize,
    ovr_rate: f64,
) -> Vec<(usize, usize, f64)> {
    let mut result = Vec::with_capacity(segments.len() + 2);
    for (s, e, rate) in segments {
        if ovr_end <= s || ovr_start >= e {
            result.push((s, e, rate));
            continue;
        }
        let isect_s = s.max(ovr_start);
        let isect_e = e.min(ovr_end);
        if s < isect_s {
            result.push((s, isect_s, rate));
        }
        result.push((isect_s, isect_e, ovr_rate));
        if isect_e < e {
            result.push((isect_e, e, rate));
        }
    }
    result
}

/// Returns the mutation rate at `pos`, or 0.0 if the position falls in an N-region
/// or gap. Segments must be sorted by start and non-overlapping.
fn rate_at(segments: &[(usize, usize, f64)], pos: usize) -> f64 {
    let idx = segments.partition_point(|&(s, _, _)| s <= pos);
    if idx == 0 {
        return 0.0;
    }
    let (_, e, rate) = segments[idx - 1];
    if pos < e { rate } else { 0.0 }
}

/// Splits segments to remove individual excluded positions (e.g. positions already
/// occupied by input variants). `excluded` must be sorted.
fn exclude_positions(
    segments: Vec<(usize, usize, f64)>,
    excluded: &[usize],
) -> Vec<(usize, usize, f64)> {
    if excluded.is_empty() {
        return segments;
    }
    let mut result = Vec::with_capacity(segments.len() + excluded.len());
    let mut ei = 0;
    for (s, e, rate) in segments {
        let mut cur = s;
        while ei < excluded.len() && excluded[ei] < e {
            let pos = excluded[ei];
            ei += 1;
            if pos < cur {
                continue;
            }
            if cur < pos {
                result.push((cur, pos, rate));
            }
            cur = pos + 1;
        }
        if cur < e {
            result.push((cur, e, rate));
        }
    }
    result
}

/// Collects `(position_0based, broken_fraction)` for every BND and INV junction
/// in `sv_records`, sorted by position. A BND contributes one junction (its
/// POS); an INV contributes two (its POS and its END) — mirroring the two
/// junctions `process_chimeric_variants` emits for an inversion. `broken_fraction`
/// is the chimeric pass's `mult`: 1.0 homozygous, 1/ploidy heterozygous — the
/// fraction of alleles that carry the junction.
///
/// END for an INV is computed exactly as the chimeric INV branch does
/// (`sv.end`, else POS + span − 1 via `SvData::span`) so the suppression window
/// references the same base the chimeric reads were placed against.
fn collect_suppressible_junctions(
    sv_records: &[Variant],
    ploidy: usize,
    use_subclone: bool,
) -> Vec<(usize, f64)> {
    let ploidy_f = (ploidy.max(1)) as f64;
    let fraction = |v: &Variant| sv_effective_fraction(v, ploidy, use_subclone);
    let mut junctions: Vec<(usize, f64)> = Vec::new();
    for v in sv_records {
        let sv = match v.alternate.as_symbolic() {
            Some(s) => s,
            None => continue,
        };
        match sv.sv_type {
            // Single point junction at POS. The chimeric pass emits genotype-mult
            // junction reads here; the regular pass also covers POS from the
            // unbroken reference (BND is coverage-neutral, DEL's coverage
            // multiplier only zeros the *interior* — flank reads crossing POS
            // still leak). broken_fraction = the chimeric pass's genotype mult.
            SvType::Bnd | SvType::Del => junctions.push((v.location, fraction(v))),
            // Two junctions (POS and END) — both breakpoints get junction reads.
            SvType::Inv => {
                let bf = fraction(v);
                junctions.push((v.location, bf));
                let end = match sv.end {
                    Some(e) => e,
                    None => match sv.span(v.location) {
                        Some(span) if span > 0 => v.location + span - 1,
                        _ => continue,
                    },
                };
                junctions.push((end, bf));
            }
            // CNV loss (cn < ploidy) is DEL-like: one POS junction at the
            // CN-derived loss fraction (matches the chimeric CNV branch's
            // (ploidy − cn)/ploidy). CNV gain (cn ≥ ploidy) is DUP-like and, like
            // <DUP>, only creates a NOVEL tandem adjacency that the regular
            // linear pass never reproduces — so it is NOT double-counted and is
            // intentionally left unsuppressed.
            SvType::Cnv => {
                if let Some(cn) = sv.copy_number {
                    let cn = cn as usize;
                    if cn < ploidy {
                        let bf = if use_subclone {
                            fraction(v)
                        } else {
                            (ploidy - cn) as f64 / ploidy_f
                        };
                        junctions.push((v.location, bf));
                    }
                }
            }
            // <DUP> tandem boundary is a novel adjacency (not in linear reads);
            // <INS> is a literal insertion. Neither double-counts the regular pass.
            SvType::Dup | SvType::Ins | SvType::Unknown => continue,
        }
    }
    junctions.sort_unstable_by_key(|&(p, _)| p);
    junctions
}

/// Largest `broken_fraction` among junctions whose position falls in the
/// half-open read window `[lo, hi)`. `positions` must be `junctions` projected
/// to positions (kept sorted) so the window can be located by binary search.
fn max_broken_fraction_in_window(
    junctions: &[(usize, f64)],
    positions: &[usize],
    lo: usize,
    hi: usize,
) -> f64 {
    let i = positions.partition_point(|&p| p < lo);
    let j = positions.partition_point(|&p| p < hi);
    junctions[i..j]
        .iter()
        .map(|&(_, bf)| bf)
        .fold(0.0_f64, f64::max)
}

/// Removes the broken-allele fraction of regular read-pairs that cross a chimeric
/// SV junction (BND, INV, DEL, and CNV-loss — see `collect_suppressible_junctions`),
/// so junction depth isn't double-counted against the chimeric pass.
///
/// A pair `(start, end)` crosses junction `j` when its R1 window
/// `[start, start + read_len)` — or, paired-end, its R2 window
/// `[end − read_len, end)` — contains `j`. Junctions in the unsequenced insert
/// gap don't cross a sequenced read and are left alone. Junctions at
/// broken_fraction = 1.0 (homozygous, or CN=0 loss) drop the pair outright and
/// consume no RNG; fractional junctions drop with probability `broken_fraction`.
/// When the contig has no suppressible junctions the input is returned untouched
/// (no RNG drawn), so non-SV and DUP/CNV-gain-only runs are byte-for-byte
/// unaffected.
fn suppress_junction_double_count(
    fragments: Vec<PlacedFragment>,
    sv_records: &[Variant],
    ploidy: usize,
    read_len: usize,
    paired_ended: bool,
    use_subclone: bool,
    rng: &mut NeatRng,
) -> Result<Vec<PlacedFragment>, GenerateReadsError> {
    let junctions = collect_suppressible_junctions(sv_records, ploidy, use_subclone);
    if junctions.is_empty() {
        return Ok(fragments);
    }
    let positions: Vec<usize> = junctions.iter().map(|&(p, _)| p).collect();
    let mut kept: Vec<PlacedFragment> = Vec::with_capacity(fragments.len());
    for placed in fragments {
        let (start, end) = (placed.start, placed.end);
        let mut bf = max_broken_fraction_in_window(&junctions, &positions, start, start + read_len);
        if paired_ended && end > read_len {
            bf = bf.max(max_broken_fraction_in_window(
                &junctions,
                &positions,
                end - read_len,
                end,
            ));
        }
        if bf <= 0.0 {
            kept.push(placed); // crosses no junction
        } else if bf >= 1.0 {
            continue; // homozygous: every allele broken — drop, no RNG draw
        } else if rng.random()? >= bf {
            kept.push(placed); // heterozygous: survived the drop coin
        }
        // else: heterozygous and dropped
    }
    Ok(kept)
}

/// Builds a sorted, contiguous list of `(start, end, multiplier)` coverage
/// segments spanning `[0, block_end)`. Default multiplier is `1.0`; each
/// symbolic SV multiplies the multiplier in its span (overlapping SVs compose
/// multiplicatively). SVs without a usable span (no END / SVLEN) or with a
/// multiplier of `1.0` are skipped silently.
fn build_coverage_multipliers(
    sv_variants: &[Variant],
    ploidy: usize,
    block_end: usize,
    use_subclone: bool,
) -> Vec<(usize, usize, f64)> {
    let mut segments: Vec<(usize, usize, f64)> = if block_end > 0 {
        vec![(0, block_end, 1.0)]
    } else {
        Vec::new()
    };
    for v in sv_variants {
        let sv = match v.alternate.as_symbolic() {
            Some(s) => s,
            None => continue,
        };
        // Convert the 0-based-stored location back to the VCF's 1-based POS
        // so SvData::span() (which expects 1-based) returns the right count.
        let pos_1based = v.location.saturating_add(1);
        let span_bases = match sv.span(pos_1based) {
            Some(n) if n > 0 => n,
            _ => {
                warn!(
                    "Symbolic SV at 1-based POS {} has no END/SVLEN — skipping coverage modulation",
                    pos_1based
                );
                continue;
            }
        };
        let base_mult =
            match coverage_multiplier_for(sv.sv_type, sv.copy_number, &v.genotype, ploidy) {
                Some(m) => m,
                None => {
                    warn!(
                        "CNV at 1-based POS {} has no INFO/CN — cannot determine copy number; \
                     skipping coverage modulation",
                        pos_1based
                    );
                    continue;
                }
            };
        let mult = if use_subclone {
            let dosage = sv_dosage_fraction(v, ploidy);
            let ccf = if dosage > 0.0 {
                (v.allele_fraction.unwrap_or(dosage) / dosage).clamp(0.0, 1.0)
            } else {
                1.0
            };
            1.0 + (base_mult - 1.0) * ccf
        } else {
            base_mult
        };
        if (mult - 1.0).abs() < f64::EPSILON {
            continue;
        }
        let (mod_start, mod_end) =
            sv_modulation_range(v.location, sv.sv_type, span_bases, block_end);
        if mod_start >= mod_end {
            continue;
        }
        segments = apply_coverage_factor(segments, mod_start, mod_end, mult);
    }
    segments
}

/// Returns the 0-based half-open coordinate range over which a symbolic SV
/// modulates coverage, given the variant's 0-based stored `location` (= VCF
/// POS − 1) and the `span_bases` reported by [`SvData::span`].
///
/// One convention for every symbolic type (VCF 4.2): POS is the anchor base
/// immediately *before* the event and is not itself affected, so the affected
/// bases run from POS+1 to END (1-based, inclusive) and the modulated range is
/// `[location + 1, location + span)` in 0-based half-open coords.
///
/// This used to differ by type — DUP/CNV/INV modulated `[location, location + span)`
/// on the basis that "POS is conventionally *inside* the affected region". That is
/// not the VCF convention, and it meant an input `<DUP>` was modulated one base
/// early while a `<DEL>` in the same file was not.
fn sv_modulation_range(
    location_0based: usize,
    sv_type: SvType,
    span_bases: usize,
    block_end: usize,
) -> (usize, usize) {
    let raw_end = location_0based.saturating_add(span_bases);
    let end = raw_end.min(block_end);
    let start = match sv_type {
        // Point events never modulate coverage today (their multiplier is 1.0, so
        // build_coverage_multipliers skips them before reaching here), but they have
        // no POS+1..END range, so don't let a future multiplier change silently shift
        // them by one.
        SvType::Ins | SvType::Bnd => location_0based.min(block_end),
        _ => location_0based.saturating_add(1).min(block_end),
    };
    (start, end)
}

/// Returns the coverage multiplier for a single symbolic SV given its type,
/// optional `INFO/CN`, genotype, and ploidy. Returns `None` only for `<CNV>`
/// records that have no copy number — those need an explicit CN, so we skip
/// rather than guess. Non-depth-modulating SV types (`<INS>`, `<INV>`,
/// breakends, `<...>` unknown tags) return `Some(1.0)`.
fn coverage_multiplier_for(
    sv_type: SvType,
    copy_number: Option<u32>,
    genotype: &Genotype,
    ploidy: usize,
) -> Option<f64> {
    let ploidy_f = (ploidy.max(1)) as f64;
    if let Some(cn) = copy_number {
        return Some(cn as f64 / ploidy_f);
    }
    match sv_type {
        SvType::Del => Some(match genotype {
            Genotype::Homozygous => 0.0,
            Genotype::Heterozygous => ((ploidy.saturating_sub(1)) as f64) / ploidy_f,
        }),
        SvType::Dup => Some(match genotype {
            // Homozygous DUP without CN: assume one extra copy per haplotype
            // (total ploidy * 2 copies) → multiplier = 2.0.
            Genotype::Homozygous => (ploidy_f + ploidy_f) / ploidy_f,
            // Heterozygous DUP without CN: one extra copy on a single haplotype
            // → multiplier = (ploidy + 1) / ploidy.
            Genotype::Heterozygous => (ploidy_f + 1.0) / ploidy_f,
        }),
        SvType::Cnv => None,
        SvType::Ins | SvType::Inv | SvType::Bnd | SvType::Unknown => Some(1.0),
    }
}

/// Multiplies the coverage factor in `[ovr_start, ovr_end)` by `factor`,
/// splitting segment boundaries as needed. Behaves like `apply_rate_override`
/// except composition is multiplicative instead of replacement.
fn apply_coverage_factor(
    segments: Vec<(usize, usize, f64)>,
    ovr_start: usize,
    ovr_end: usize,
    factor: f64,
) -> Vec<(usize, usize, f64)> {
    let mut result = Vec::with_capacity(segments.len() + 2);
    for (s, e, mult) in segments {
        if ovr_end <= s || ovr_start >= e {
            result.push((s, e, mult));
            continue;
        }
        let isect_s = s.max(ovr_start);
        let isect_e = e.min(ovr_end);
        if s < isect_s {
            result.push((s, isect_s, mult));
        }
        result.push((isect_s, isect_e, mult * factor));
        if isect_e < e {
            result.push((isect_e, e, mult));
        }
    }
    result
}

/// Intersects `[region_start, region_end)` with the multiplier segments,
/// returning the sub-regions clipped to the region in coordinate order.
/// Sub-regions outside any segment (which shouldn't happen because the
/// segments span `[0, block_end)`) implicitly get multiplier 1.0 only if
/// `coverage_multipliers` is empty — in which case the whole region is
/// returned as one piece.
fn split_region_by_multipliers(
    region_start: usize,
    region_end: usize,
    segments: &[(usize, usize, f64)],
) -> Vec<(usize, usize, f64)> {
    if segments.is_empty() {
        return vec![(region_start, region_end, 1.0)];
    }
    let mut out = Vec::new();
    for &(s, e, m) in segments {
        let lo = s.max(region_start);
        let hi = e.min(region_end);
        if lo < hi {
            out.push((lo, hi, m));
        }
    }
    out
}

/// Pick a chimeric-junction offset that splits `frag_len` so both
/// pieces are at least `read_len / 4` bases long. The offset controls
/// how the stitched fragment is divided into left and right reference
/// pieces (`len1 = offset`, `len2 = frag_len - offset`). When offset is
/// too small or too large, BWA aligns the read as a regular alignment
/// with the dominant piece and emits a short / low-MAPQ split alignment
/// for the minority piece — which Manta filters out of its candidate
/// pool (MAPQ < 20). The `read_len / 4` floor matches BWA-MEM's typical
/// anchor threshold of ~30-40 bp for confident split alignment.
///
/// Added for #224. Before this constraint the offset was sampled
/// uniformly from `[1, min(frag_len-1, read_len-1)]`, which produced
/// chimeric reads with offset=1 or 2 about 1-2% of the time —
/// individually small but enough to drag the BND candidate-pool below
/// Manta's somatic threshold across many junctions.
fn balanced_chimeric_offset(
    frag_len: usize,
    read_len: usize,
    rng: &mut NeatRng,
) -> Result<usize, GenerateReadsError> {
    let min_offset = (read_len / 4).max(1);
    // Don't constrain past where the fragment can support both pieces.
    let max_offset = frag_len.saturating_sub(min_offset);
    let lo = min_offset as i64;
    let hi = max_offset.max(min_offset + 1) as i64;
    rng.range_i64(lo, hi)
        .map(|v| v as usize)
        .map_err(GenerateReadsError::CliError)
}

/// Scales a base coverage value by a non-negative multiplier, rounding to the
/// nearest integer and clamping at 0. Returns 0 for non-finite or negative
/// inputs (defensive; build_coverage_multipliers shouldn't produce either).
fn scale_coverage(base: usize, multiplier: f64) -> usize {
    if !multiplier.is_finite() || multiplier <= 0.0 {
        return 0;
    }
    (base as f64 * multiplier).round() as usize
}

/// Intersects a set of non-N regions (block-local coordinates) with BED records
/// (global contig coordinates) and returns only the overlapping sub-intervals,
/// still expressed in block-local coordinates.
///
/// `block_offset` is `SequenceBlock::ref_start` — the global contig position at
/// which this block begins.
fn intersect_with_bed(
    regions: &[&SequenceMap],
    bed_records: &[BedRecord],
    block_offset: usize,
) -> Vec<SequenceMap> {
    let mut out = Vec::new();
    for region in regions {
        let global_start = region.start + block_offset;
        let global_end = region.end + block_offset;
        for bed in bed_records {
            let isect_start = global_start.max(bed.start);
            let isect_end = global_end.min(bed.end);
            if isect_start < isect_end {
                out.push(SequenceMap::from(
                    RegionType::NonNRegion,
                    isect_start - block_offset,
                    isect_end - block_offset,
                ));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    // Junction-suppression tests predate haplotype-aware placement and are written
    // in plain (start, end) spans, which is still the right level for them:
    // suppression is about where a fragment sits, not which haplotype it came from.
    // Convert at the boundary rather than rewriting the cases.
    fn placed(v: Vec<(usize, usize)>) -> Vec<PlacedFragment> {
        v.into_iter().map(PlacedFragment::from).collect()
    }
    fn spans(v: &[PlacedFragment]) -> Vec<(usize, usize)> {
        v.iter().map(|f| (f.start, f.end)).collect()
    }

    use super::*;
    use eidolon_core::structs::bed_record::BedRecord;
    use eidolon_core::structs::sequence_block::{RegionType, SequenceMap};
    use eidolon_core::structs::variants::AlternateType;

    /// The four VCF 4.2 §5.4 breakend forms, and specifically WHICH PIECE gets reverse-
    /// complemented. This is the property nothing tested: `bnd_fastq.rs` drives the input
    /// path (where the parser sets the flags) and only asserts that some read is named
    /// "EIDOLON_chimeric", so a direct join and a head-to-head join were indistinguishable
    /// to the whole suite.
    ///
    /// That gap let de novo BNDs ship a truth VCF saying `t]p]` (reverse-complemented)
    /// while the reads were built as case 4 (a direct join) — which is simply a deletion
    /// or duplication, and is what Manta correctly called them on Delta.
    fn bnd_sv(join_after: bool, mate_extends_right: bool) -> SvData {
        let mut sv = SvData::new("N]chr1:5001]", SvType::Bnd);
        sv.mate_contig = Some("chr1".to_string());
        sv.mate_pos = Some(5001); // 1-based, so 5000 0-based
        sv.bnd_join_after = join_after;
        sv.bnd_mate_extends_right = mate_extends_right;
        sv
    }

    fn bnd_reference() -> HashMap<String, Vec<Nucleotide>> {
        let mut r = HashMap::new();
        r.insert("chr1".to_string(), vec![Nucleotide::A; 10_000]);
        r
    }

    /// #224: a chimeric fragment whose junction sits too close to either end produces a
    /// read with an anchor too short for BWA to split-align — the regression that shipped
    /// in v1.13.0 and was fixed by the `read_len / 4` floor. Nothing tested that floor,
    /// so replacing it with `1` was silently safe.
    #[test]
    fn balanced_chimeric_offset_keeps_a_split_alignable_anchor_either_side() {
        let mut rng = NeatRng::new_from_seed(&vec!["chimeric offset".to_string()]).unwrap();
        let (read_len, frag_len) = (151usize, 400usize);
        let floor = read_len / 4; // 37
        for _ in 0..1000 {
            let off = balanced_chimeric_offset(frag_len, read_len, &mut rng).unwrap();
            assert!(
                off >= floor,
                "offset {off} leaves a {off}bp anchor before the junction, below the \
                 read_len/4 floor of {floor} (#224)"
            );
            assert!(
                frag_len - off >= floor,
                "offset {off} leaves a {}bp anchor after the junction, below {floor} (#224)",
                frag_len - off
            );
        }
    }

    /// Degenerate case: a fragment no longer than the read cannot honour the floor on both
    /// sides. It must still return a usable offset rather than panic or return 0.
    #[test]
    fn balanced_chimeric_offset_survives_a_fragment_no_longer_than_the_read() {
        let mut rng = NeatRng::new_from_seed(&vec!["degenerate".to_string()]).unwrap();
        for _ in 0..100 {
            let off = balanced_chimeric_offset(151, 151, &mut rng).unwrap();
            assert!(off >= 1, "offset must be at least 1, got {off}");
        }
    }

    #[test]
    fn bnd_pieces_reverse_complement_exactly_the_spec_cases() {
        let reference = bnd_reference();
        let (pos, len1, len2) = (1000usize, 100usize, 50usize);

        // Case 1, t[p[ : REF[..=pos] + MATE[mate_pos..]. Direct, nothing reversed.
        let (a, b) =
            get_bnd_pieces("chr1", pos, &bnd_sv(true, true), len1, len2, &reference).unwrap();
        assert_eq!((a.1, a.2, a.3), (901, 1001, false), "case 1 anchor piece");
        assert_eq!(
            (b.1, b.2, b.3),
            (5000, 5050, false),
            "case 1 mate piece must NOT be reversed"
        );

        // Case 2, t]p] : REF[..=pos] + revcomp(MATE[..=mate_pos]). This is the form de
        // novo BNDs declare, so the mate piece MUST be reverse-complemented.
        let (a, b) =
            get_bnd_pieces("chr1", pos, &bnd_sv(true, false), len1, len2, &reference).unwrap();
        assert_eq!((a.1, a.2, a.3), (901, 1001, false), "case 2 anchor piece");
        assert_eq!(
            (b.1, b.2, b.3),
            (4951, 5001, true),
            "case 2 mate piece MUST be reverse-complemented"
        );

        // Case 3, [p[t : revcomp(MATE[mate_pos..]) + REF[pos..]. Mate piece comes first.
        let (a, b) =
            get_bnd_pieces("chr1", pos, &bnd_sv(false, true), len1, len2, &reference).unwrap();
        assert_eq!(
            (a.1, a.2, a.3),
            (5000, 5100, true),
            "case 3 mate piece leads, reversed"
        );
        assert_eq!((b.1, b.2, b.3), (1000, 1050, false), "case 3 anchor piece");

        // Case 4, ]p]t : MATE[..=mate_pos] + REF[pos..]. Direct, nothing reversed — the
        // layout a de novo BND was silently getting from the false/false defaults.
        let (a, b) =
            get_bnd_pieces("chr1", pos, &bnd_sv(false, false), len1, len2, &reference).unwrap();
        assert_eq!(
            (a.1, a.2, a.3),
            (4901, 5001, false),
            "case 4 mate piece leads, NOT reversed"
        );
        assert_eq!((b.1, b.2, b.3), (1000, 1050, false), "case 4 anchor piece");
    }

    /// The regression proper: the geometry a de novo BND's ALT declares must be the
    /// geometry its reads are built from. Distinct from the case table above, which pins
    /// the generator; this pins the two ends AGREEING.
    #[test]
    fn denovo_bnd_alt_form_produces_a_reverse_complemented_junction() {
        let reference = bnd_reference();
        // `t]p]` is what sv_model.rs emits for every de novo BND.
        let (_mc, _mp, join_after, mate_right) =
            eidolon_core::structs::variants::parse_bnd_alt_for_test("N]chr1:5001]");
        let (_a, b) = get_bnd_pieces(
            "chr1",
            1000,
            &bnd_sv(join_after, mate_right),
            100,
            50,
            &reference,
        )
        .unwrap();
        assert!(
            b.3,
            "a `t]p]` breakend joins a REVERSE-COMPLEMENTED piece (VCF 4.2 §5.4); a \
             direct join here is a deletion or duplication, not a breakend"
        );
    }

    #[test]
    fn test_split_contig_into_chunks_covers_contig_without_gaps_or_overlap() {
        // 1 Mbp chunks over a 2.5 Mbp contig → 3 even chunks, contiguous, covering [0, len).
        let chunks = split_contig_into_chunks(2_500_000, 1_000_000);
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks.first().unwrap().0, 0);
        assert_eq!(chunks.last().unwrap().1, 2_500_000);
        for w in chunks.windows(2) {
            assert_eq!(w[0].1, w[1].0, "chunks must be contiguous (no gap/overlap)");
        }
        // Even split: sizes differ by at most 1.
        let sizes: Vec<usize> = chunks.iter().map(|(s, e)| e - s).collect();
        assert!(sizes.iter().max().unwrap() - sizes.iter().min().unwrap() <= 1);
    }

    #[test]
    fn test_split_contig_into_chunks_edge_cases() {
        // Contig smaller than the chunk size → a single whole-contig chunk.
        assert_eq!(split_contig_into_chunks(500, 1_000_000), vec![(0, 500)]);
        // Empty contig → one [0,0) chunk so it still yields a result.
        assert_eq!(split_contig_into_chunks(0, 1_000_000), vec![(0, 0)]);
        // chunk_size 0 is treated as 1 (guarded) and must not divide-by-zero.
        assert!(!split_contig_into_chunks(10, 0).is_empty());
    }

    #[test]
    fn test_resolve_chunk_size_config_modes() {
        let mut cfg = RunConfiguration::default();
        // None (default) → disabled: one chunk spans the whole contig.
        cfg.chunk_size = None;
        assert_eq!(resolve_chunk_size(&cfg), usize::MAX);
        // Some(0) → disabled (explicit).
        cfg.chunk_size = Some(0);
        assert_eq!(resolve_chunk_size(&cfg), usize::MAX);
        // Some(n) → fixed opt-in chunk size.
        cfg.chunk_size = Some(5_000_000);
        assert_eq!(resolve_chunk_size(&cfg), 5_000_000);
    }

    #[test]
    fn test_intersect_with_bed() {
        let r1 = SequenceMap::from(RegionType::NonNRegion, 100, 200);
        let r2 = SequenceMap::from(RegionType::NonNRegion, 300, 400);
        let regions = vec![&r1, &r2];

        let b1 = BedRecord::new_bed_record("chr1".to_string(), 150, 350).unwrap();
        let bed_records = vec![b1];

        // block_offset = 0
        let result = intersect_with_bed(&regions, &bed_records, 0);
        assert_eq!(result.len(), 2);
        // Intersection with r1 [100, 200] and b1 [150, 350] -> [150, 200]
        assert_eq!(result[0].start, 150);
        assert_eq!(result[0].end, 200);
        // Intersection with r2 [300, 400] and b1 [150, 350] -> [300, 350]
        assert_eq!(result[1].start, 300);
        assert_eq!(result[1].end, 350);

        // block_offset = 1000
        // r1 global [1100, 1200], r2 global [1300, 1400]
        let b2 = BedRecord::new_bed_record("chr1".to_string(), 1150, 1350).unwrap();
        let result2 = intersect_with_bed(&regions, &[b2], 1000);
        assert_eq!(result2.len(), 2);
        assert_eq!(result2[0].start, 150); // global 1150 - 1000
        assert_eq!(result2[0].end, 200); // global 1200 - 1000
        assert_eq!(result2[1].start, 300); // global 1300 - 1000
        assert_eq!(result2[1].end, 350); // global 1350 - 1000
    }

    /// #591: an input record that asserts no alternate allele must not become one.
    ///
    /// Three shapes, all of which used to be silently GENERATED and then written into the
    /// truth VCF contradicting themselves. Measured on H1N1 before the fix:
    ///
    ///   REF==ALT   ->  GT=1/1  AD=0,42  AF=1.0000   (there is no alt allele at all)
    ///   GT=./.     ->  GT=./.  AD=0,44  AF=1.0000   (the GT declines to call it)
    ///   GT=0/0     ->  GT=0/0  AD=18,17 AF=0.4857   (the GT calls reference twice)
    ///
    /// Each is a false negative for every caller scored against that truth, and inflates
    /// the recall denominator by a record no caller could ever produce.
    #[test]
    fn input_records_that_assert_no_alternate_allele_are_dropped() {
        let mk = |location: usize,
                  reference: Vec<Nucleotide>,
                  alt: Vec<Nucleotide>,
                  gt: &str,
                  vt: VariantType| Variant {
            location,
            reference,
            alternate: AlternateType::Literal(alt),
            variant_type: vt,
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            genotype_str: gt.to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
            provenance: Provenance::InputVcf,
        };
        use Nucleotide::{A, C, T};

        let raw = HashMap::from([(
            "chr1".to_string(),
            vec![
                // dropped: no alternate allele exists
                mk(100, vec![A], vec![A], "1/1", VariantType::SNP),
                // dropped: the genotype declines to call
                mk(200, vec![A], vec![T], "./.", VariantType::SNP),
                mk(250, vec![A], vec![T], ".|.", VariantType::SNP),
                // dropped: the genotype calls reference on every allele
                mk(300, vec![A], vec![T], "0/0", VariantType::SNP),
                mk(350, vec![A], vec![T], "0|0", VariantType::SNP),
                // MUST NOT FIRE — every one of these is a real variant call
                mk(400, vec![A], vec![T], "1/1", VariantType::SNP),
                mk(500, vec![A], vec![T], "0/1", VariantType::SNP),
                mk(600, vec![A], vec![T], "1|0", VariantType::SNP),
                // a PARTIAL no-call still calls one alt allele, so it is a variant
                mk(700, vec![A], vec![T], "./1", VariantType::SNP),
                // REF and ALT differ in length: an indel, not a null variant
                mk(800, vec![A, C], vec![A], "1/1", VariantType::Deletion),
            ],
        )]);

        let kept = filter_input_vcf(raw);
        let mut locs: Vec<usize> = kept
            .get("chr1")
            .map(|v| v.iter().map(|x| x.location).collect())
            .unwrap_or_default();
        locs.sort_unstable();
        assert_eq!(
            locs,
            vec![400, 500, 600, 700, 800],
            "expected only the records that actually assert an alternate allele to survive; \
             dropping a real call here would silently shrink the truth set"
        );
    }

    /// A partial no-call is kept, but its `.` becomes `0` — #591 follow-up.
    ///
    /// Before: `./1` was written to the truth VCF verbatim and generated at AF 1.0000,
    /// because `dosage_fraction` skips `.` alleles and scored it alt=1 of total=1.
    /// Two things were wrong with that. It invents a second alt allele the input never
    /// mentioned, and — measured with bcftools — `./1` is classified as neither het nor
    /// hom: `GT="het"` selects only `0/1`, `GT="hom"` only `1/1`, and `./1` appears solely
    /// under `GT="mis"`. A truth record like that drops out of any zygosity-stratified
    /// analysis silently.
    #[test]
    fn a_partial_no_call_is_normalized_to_a_reference_allele() {
        assert_eq!(normalize_partial_no_call("./1").as_deref(), Some("0/1"));
        assert_eq!(normalize_partial_no_call("1/.").as_deref(), Some("1/0"));
        assert_eq!(normalize_partial_no_call(".|1").as_deref(), Some("0|1"));
        // generalises past diploid: one called alt of three copies
        assert_eq!(normalize_partial_no_call("././1").as_deref(), Some("0/0/1"));

        // MUST NOT FIRE: nothing to normalize, so nothing is rewritten.
        for gt in ["0/1", "1/1", "1|0", "0/0/1", "0"] {
            assert_eq!(
                normalize_partial_no_call(gt),
                None,
                "{gt} has no uncalled allele and must be left exactly as supplied"
            );
        }
    }

    /// The normalized record must reach the kept set carrying the rewritten GT — a unit
    /// test on the helper alone would pass even if it were never wired in.
    #[test]
    fn the_normalized_partial_no_call_is_what_gets_kept() {
        let v = Variant {
            location: 700,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::T]),
            variant_type: VariantType::SNP,
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            genotype_str: "./1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
            provenance: Provenance::InputVcf,
        };
        let kept = filter_input_vcf(HashMap::from([("chr1".to_string(), vec![v])]));
        let got = &kept.get("chr1").expect("record kept")[0];
        assert_eq!(got.genotype_str, "0/1", "the `.` must be rewritten to `0`");
        assert_eq!(
            got.genotype,
            Genotype::Heterozygous,
            "one alt of two copies is heterozygous, not homozygous"
        );
        // 1 of 2 copies, not 1 of 1 — this is the AF that reaches read generation.
        assert!(
            (got.dosage_fraction() - 0.5).abs() < 1e-9,
            "expected dosage 0.5, got {}",
            got.dosage_fraction()
        );
    }

    #[test]
    fn no_call_and_hom_ref_genotypes_are_recognised_on_the_raw_string() {
        // gt_from_str cannot express either: it skips '.' alleles (so "./." falls through to
        // Homozygous) and stops at the first '0' (so "0/0" reads as Heterozygous). Both have
        // to be caught before that conversion.
        for gt in ["./.", ".|.", ".", "./././."] {
            assert!(is_no_call(gt), "{gt} is a no-call");
            assert!(!is_hom_ref(gt), "{gt} is not homozygous reference");
        }
        for gt in ["0/0", "0|0", "0", "0/0/0"] {
            assert!(is_hom_ref(gt), "{gt} is homozygous reference");
            assert!(!is_no_call(gt), "{gt} is not a no-call");
        }
        // MUST NOT FIRE: each of these calls at least one alternate allele.
        for gt in ["0/1", "1/1", "1|0", "./1", "1/.", "0/2"] {
            assert!(!is_no_call(gt), "{gt} calls an allele and is not a no-call");
            assert!(!is_hom_ref(gt), "{gt} carries an alt and is not hom-ref");
        }
    }

    #[test]
    fn test_filter_input_vcf() {
        use crate::eidolon_core::structs::variants::{Genotype, Provenance, VariantType};
        use eidolon_core::structs::nucleotides::Nucleotide;
        let mut raw = HashMap::new();
        let v1 = Variant {
            location: 100,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::T]),
            variant_type: VariantType::SNP,
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            genotype_str: "1/1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
            provenance: Provenance::InputVcf,
        };
        let v2 = Variant {
            location: 200,
            reference: vec![Nucleotide::A, Nucleotide::T],
            alternate: AlternateType::Literal(vec![Nucleotide::C, Nucleotide::G]),
            variant_type: VariantType::Complex,
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            genotype_str: "1/1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
            provenance: Provenance::InputVcf,
        };
        // Symbolic SV — tagged Complex but must NOT be dropped by
        // filter_input_vcf: gen_reads uses it downstream for coverage
        // modulation and round-trips the raw ALT to the output VCF.
        use eidolon_core::structs::variants::{SvData, SvType};
        let v3 = Variant {
            location: 500,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(SvData::new("<DEL>", SvType::Del)),
            variant_type: VariantType::Complex,
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            genotype_str: "1/1".to_string(),
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            sample: vec![],
            provenance: Provenance::InputVcf,
        };
        raw.insert("chr1".to_string(), vec![v1.clone(), v2, v3]);

        let filtered = filter_input_vcf(raw);
        assert_eq!(filtered.len(), 1);
        // Literal Complex (v2) is dropped; SNP (v1) and symbolic <DEL> (v3) are kept.
        assert_eq!(filtered["chr1"].len(), 2);
        let locs: Vec<usize> = filtered["chr1"].iter().map(|v| v.location).collect();
        assert!(locs.contains(&100));
        assert!(locs.contains(&500));
    }

    fn sv_variant_with_span(
        location_0based: usize,
        end_1based: usize,
        sv_type: SvType,
        genotype: Genotype,
        copy_number: Option<u32>,
    ) -> Variant {
        use eidolon_core::structs::nucleotides::Nucleotide;
        use eidolon_core::structs::variants::SvData;
        let mut sv = SvData::new(
            match sv_type {
                SvType::Del => "<DEL>",
                SvType::Dup => "<DUP>",
                SvType::Cnv => "<CNV>",
                SvType::Ins => "<INS>",
                SvType::Inv => "<INV>",
                _ => "<UNKNOWN>",
            },
            sv_type,
        );
        sv.end = Some(end_1based);
        sv.copy_number = copy_number;
        let genotype_str = match genotype {
            Genotype::Homozygous => "1/1".to_string(),
            Genotype::Heterozygous => "0/1".to_string(),
        };
        Variant {
            location: location_0based,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Symbolic(sv),
            variant_type: VariantType::Complex,
            genotype,
            genotype_str,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: vec![],
            provenance: eidolon_core::structs::variants::Provenance::Denovo,
            sample: vec![],
        }
    }

    #[test]
    fn sv_subclone_ccf_scales_depth_and_junction_fraction() {
        let mut del = sv_variant_with_span(100, 200, SvType::Del, Genotype::Heterozygous, None);
        del.allele_fraction = Some(0.25); // 0.5 dosage × 0.5 CCF

        let historical = build_coverage_multipliers(&[del.clone()], 2, 300, false);
        assert!(historical.iter().any(|&(_, _, m)| (m - 0.5).abs() < 1e-9));

        let subclonal = build_coverage_multipliers(&[del.clone()], 2, 300, true);
        assert!(subclonal.iter().any(|&(_, _, m)| (m - 0.75).abs() < 1e-9));

        let junctions = collect_suppressible_junctions(&[del], 2, true);
        assert_eq!(junctions, vec![(100, 0.25)]);
    }

    #[test]
    fn applying_sv_subclone_model_stamps_ccf_and_vaf() {
        use crate::gen_reads::utils::subclone::Subclone;

        let model = SubcloneModel::new(vec![Subclone {
            ccf: 0.5,
            weight: 1.0,
        }])
        .unwrap();
        let mut variants = vec![sv_variant_with_span(
            100,
            200,
            SvType::Bnd,
            Genotype::Heterozygous,
            None,
        )];
        let mut rng = NeatRng::new_from_seed(&vec!["sv-ccf".to_string()]).unwrap();

        apply_sv_subclone_model(&mut variants, Some(&model), 2, Some(0.6), &mut rng).unwrap();

        assert!((variants[0].allele_fraction.unwrap() - 0.25).abs() < 1e-9);
        let info = variants[0].info.as_deref().unwrap();
        assert!(info.contains("EIDOLON_CCF=0.5000"));
        assert!(info.contains("EIDOLON_VAF=0.1500"));
    }

    #[test]
    fn no_sv_subclone_model_preserves_variant_and_rng_path() {
        let original = sv_variant_with_span(100, 200, SvType::Bnd, Genotype::Heterozygous, None);
        let mut variants = vec![original.clone()];
        let mut rng = NeatRng::new_from_seed(&vec!["sv-no-ccf".to_string()]).unwrap();

        apply_sv_subclone_model(&mut variants, None, 2, Some(0.6), &mut rng).unwrap();

        assert_eq!(variants[0], original);
    }

    #[test]
    fn test_collect_suppressible_junctions() {
        // BND (homozygous) at 100 → one junction (100, 1.0).
        // INV (heterozygous, ploidy 2) over [200, 300] → (200, 0.5) and (300, 0.5).
        let bnd = sv_variant_with_span(100, 0, SvType::Bnd, Genotype::Homozygous, None);
        let inv = sv_variant_with_span(200, 300, SvType::Inv, Genotype::Heterozygous, None);
        let j = collect_suppressible_junctions(&[inv, bnd], 2, false);
        assert_eq!(j.len(), 3, "BND→1 junction, INV→2");
        // sorted by position
        assert_eq!(j[0].0, 100);
        assert!(
            (j[0].1 - 1.0).abs() < 1e-9,
            "homozygous broken_fraction = 1.0"
        );
        assert_eq!(j[1].0, 200);
        assert!(
            (j[1].1 - 0.5).abs() < 1e-9,
            "het broken_fraction = 1/ploidy"
        );
        assert_eq!(j[2].0, 300);
        assert!((j[2].1 - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_collect_del_loss_suppressed_dup_gain_not() {
        // DEL → one POS junction at the genotype mult.
        // CNV-loss (cn < ploidy) → one POS junction at (ploidy−cn)/ploidy.
        // DUP and CNV-gain (cn ≥ ploidy) create novel tandem adjacencies → NO junction.
        let del = sv_variant_with_span(100, 200, SvType::Del, Genotype::Homozygous, None);
        let dup = sv_variant_with_span(300, 400, SvType::Dup, Genotype::Homozygous, None);
        let cnv_loss = sv_variant_with_span(500, 600, SvType::Cnv, Genotype::Homozygous, Some(1));
        let cnv_gain = sv_variant_with_span(700, 800, SvType::Cnv, Genotype::Homozygous, Some(4));
        let j = collect_suppressible_junctions(&[del, dup, cnv_loss, cnv_gain], 2, false);
        // Only DEL (100) and CNV-loss (500) contribute.
        assert_eq!(
            j.len(),
            2,
            "expected DEL + CNV-loss junctions only, got {j:?}"
        );
        assert_eq!(j[0].0, 100);
        assert!((j[0].1 - 1.0).abs() < 1e-9, "homozygous DEL → 1.0");
        assert_eq!(j[1].0, 500);
        // CNV cn=1, ploidy 2 → (2−1)/2 = 0.5.
        assert!(
            (j[1].1 - 0.5).abs() < 1e-9,
            "CNV-loss cn=1 → (ploidy−cn)/ploidy = 0.5"
        );
    }

    #[test]
    fn test_collect_cnv_full_loss_is_one() {
        // CN=0 (full loss) → broken_fraction (ploidy−0)/ploidy = 1.0.
        let cnv0 = sv_variant_with_span(500, 600, SvType::Cnv, Genotype::Homozygous, Some(0));
        let j = collect_suppressible_junctions(&[cnv0], 2, false);
        assert_eq!(j.len(), 1);
        assert!((j[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_suppress_homozygous_drops_crossing_keeps_flank() {
        // Single-end, read_len 100, homozygous BND junction at 500.
        let bnd = sv_variant_with_span(500, 0, SvType::Bnd, Genotype::Homozygous, None);
        let frags = vec![
            (450, 550), // R1 [450,550) contains 500 → crosses → drop
            (300, 400), // flank-left → keep
            (600, 700), // flank-right → keep
        ];
        let mut rng = NeatRng::new_from_seed(&vec!["bp".to_string()]).unwrap();
        let kept =
            suppress_junction_double_count(placed(frags), &[bnd], 2, 100, false, false, &mut rng)
                .unwrap();
        assert_eq!(spans(&kept), vec![(300, 400), (600, 700)]);
    }

    #[test]
    fn test_suppress_no_sv_returns_unchanged() {
        // No BND/INV → input returned verbatim, no RNG consumed.
        let frags = vec![(0, 100), (200, 300)];
        let mut rng = NeatRng::new_from_seed(&vec!["bp".to_string()]).unwrap();
        let kept = suppress_junction_double_count(
            placed(frags.clone()),
            &[],
            2,
            100,
            false,
            false,
            &mut rng,
        )
        .unwrap();
        assert_eq!(spans(&kept), frags);
    }

    #[test]
    fn test_suppress_paired_gap_not_crossed() {
        // Paired, read_len 100. Fragment [400,700]: R1 [400,500), R2 [600,700).
        // A junction at 550 sits in the unsequenced gap [500,600) → NOT crossed.
        let in_gap = sv_variant_with_span(550, 0, SvType::Bnd, Genotype::Homozygous, None);
        let mut rng = NeatRng::new_from_seed(&vec!["bp".to_string()]).unwrap();
        let kept = suppress_junction_double_count(
            placed(vec![(400, 700)]),
            &[in_gap],
            2,
            100,
            true,
            false,
            &mut rng,
        )
        .unwrap();
        assert_eq!(
            spans(&kept),
            vec![(400, 700)],
            "gap junction must not suppress"
        );
        // A junction at 450 sits in R1 → crossed → dropped.
        let in_r1 = sv_variant_with_span(450, 0, SvType::Bnd, Genotype::Homozygous, None);
        let kept2 = suppress_junction_double_count(
            placed(vec![(400, 700)]),
            &[in_r1],
            2,
            100,
            true,
            false,
            &mut rng,
        )
        .unwrap();
        assert!(kept2.is_empty(), "R1-crossing pair must be suppressed");
    }

    #[test]
    fn test_suppress_heterozygous_partial() {
        // Het junction (broken_fraction 0.5): over many identical crossing pairs
        // ~half survive. Deterministic seed; assert neither ~0 nor ~all.
        let het = sv_variant_with_span(500, 0, SvType::Bnd, Genotype::Heterozygous, None);
        let frags: Vec<(usize, usize)> = vec![(450usize, 550usize); 1000];
        let mut rng = NeatRng::new_from_seed(&vec!["bp het".to_string()]).unwrap();
        let kept =
            suppress_junction_double_count(placed(frags), &[het], 2, 100, false, false, &mut rng)
                .unwrap();
        assert!(
            (300..700).contains(&kept.len()),
            "expected ~50% of 1000 het-crossing pairs kept, got {}",
            kept.len()
        );
    }

    #[test]
    fn test_sv_modulation_range_del_skips_anchor() {
        // POS=101 (1-based) → location_0based=100. END=200 (1-based incl) so
        // span = 100. DEL anchor (POS, 0-based 100) is NOT deleted: range
        // starts at 101. End is 200 (= POS + span - 1 + 1 = 0-based half-open).
        assert_eq!(sv_modulation_range(100, SvType::Del, 100, 1000), (101, 200));
        // Single-base DEL where span_bases==1 collapses to empty (just the
        // anchor) — caller must skip via mod_start >= mod_end.
        let (s, e) = sv_modulation_range(100, SvType::Del, 1, 1000);
        assert!(
            s >= e,
            "expected empty range for span==1 DEL, got [{s}, {e})"
        );
    }

    #[test]
    fn test_sv_modulation_range_excludes_the_anchor_for_every_type() {
        // location=100 -> POS=101 (1-based); span=100 -> END = 101+100-1 = 200.
        // VCF 4.2: the anchor at POS is NOT affected, so the affected bases are
        // 1-based 102..200, i.e. 0-based [101, 200) — identical to DEL.
        //
        // These used to be (100, 200): DUP/CNV/INV modulated one base early, so an
        // input <DUP> and an input <DEL> in the same file used different conventions.
        assert_eq!(sv_modulation_range(100, SvType::Dup, 100, 1000), (101, 200));
        assert_eq!(sv_modulation_range(100, SvType::Cnv, 100, 1000), (101, 200));
        assert_eq!(sv_modulation_range(100, SvType::Inv, 100, 1000), (101, 200));
        // Same input, same answer as DEL — that equality is the property that was
        // missing, not an incidental detail.
        assert_eq!(
            sv_modulation_range(100, SvType::Dup, 100, 1000),
            sv_modulation_range(100, SvType::Del, 100, 1000)
        );
    }

    #[test]
    fn test_sv_modulation_range_clipped_to_block_end() {
        // SV running past block_end gets clipped on both ends.
        assert_eq!(sv_modulation_range(95, SvType::Del, 100, 110), (96, 110));
        // Was (95, 110): DUP now excludes its anchor like every other type.
        assert_eq!(sv_modulation_range(95, SvType::Dup, 100, 110), (96, 110));
        // Start clipped above block_end → empty range.
        let (s, e) = sv_modulation_range(150, SvType::Del, 50, 100);
        assert!(s >= e);
    }

    #[test]
    fn test_coverage_multiplier_for() {
        // DEL without CN: hom = 0, het = (ploidy-1)/ploidy
        assert_eq!(
            coverage_multiplier_for(SvType::Del, None, &Genotype::Homozygous, 2),
            Some(0.0)
        );
        assert_eq!(
            coverage_multiplier_for(SvType::Del, None, &Genotype::Heterozygous, 2),
            Some(0.5)
        );

        // DUP without CN: hom = 2.0, het = 1.5 on diploid
        assert_eq!(
            coverage_multiplier_for(SvType::Dup, None, &Genotype::Homozygous, 2),
            Some(2.0)
        );
        assert_eq!(
            coverage_multiplier_for(SvType::Dup, None, &Genotype::Heterozygous, 2),
            Some(1.5)
        );

        // CN-driven (any type, but most useful for CNV): multiplier = CN / ploidy
        assert_eq!(
            coverage_multiplier_for(SvType::Cnv, Some(4), &Genotype::Homozygous, 2),
            Some(2.0)
        );
        assert_eq!(
            coverage_multiplier_for(SvType::Cnv, Some(1), &Genotype::Heterozygous, 2),
            Some(0.5)
        );
        assert_eq!(
            coverage_multiplier_for(SvType::Cnv, Some(0), &Genotype::Homozygous, 2),
            Some(0.0)
        );

        // CNV without CN: None — caller must skip and warn
        assert_eq!(
            coverage_multiplier_for(SvType::Cnv, None, &Genotype::Homozygous, 2),
            None
        );

        // Non-depth-modulating SVs: 1.0
        for t in [SvType::Ins, SvType::Inv, SvType::Bnd, SvType::Unknown] {
            assert_eq!(
                coverage_multiplier_for(t, None, &Genotype::Heterozygous, 2),
                Some(1.0)
            );
        }
    }

    #[test]
    fn test_apply_coverage_factor_composes_multiplicatively() {
        // Whole-block default segment, halve [20, 50) -> still gives 1.0 outside.
        let segs = vec![(0usize, 100usize, 1.0f64)];
        let halved = apply_coverage_factor(segs, 20, 50, 0.5);
        assert_eq!(halved, vec![(0, 20, 1.0), (20, 50, 0.5), (50, 100, 1.0)]);

        // Apply a 2× factor on a sub-range that already has 0.5: composes to 1.0.
        let doubled = apply_coverage_factor(halved, 30, 40, 2.0);
        assert_eq!(
            doubled,
            vec![
                (0, 20, 1.0),
                (20, 30, 0.5),
                (30, 40, 1.0),
                (40, 50, 0.5),
                (50, 100, 1.0),
            ]
        );
    }

    #[test]
    fn test_build_coverage_multipliers_hom_del_zeros_span() {
        // <DEL> at 1-based POS=101 (0-based location=100), END=200 (1-based inclusive).
        // VCF semantics: POS (base at index 100) is the anchor and is NOT deleted;
        // bases POS+1..=END are. So modulation runs over 0-based [101, 200).
        let svs = vec![sv_variant_with_span(
            100,
            200,
            SvType::Del,
            Genotype::Homozygous,
            None,
        )];
        let segs = build_coverage_multipliers(&svs, 2, 500, false);
        assert_eq!(segs, vec![(0, 101, 1.0), (101, 200, 0.0), (200, 500, 1.0)]);
    }

    #[test]
    fn test_build_coverage_multipliers_het_dup_inflates_span() {
        // <DUP> heterozygous on diploid → multiplier 1.5.
        let svs = vec![sv_variant_with_span(
            50,
            149,
            SvType::Dup,
            Genotype::Heterozygous,
            None,
        )];
        // POS = 51 (1-based); span = 149 - 51 + 1 = 99. The anchor at POS is not
        // duplicated, so the affected range is [51, 50 + 99) = [51, 149).
        // Was [50, 149) — one base early.
        let segs = build_coverage_multipliers(&svs, 2, 300, false);
        assert_eq!(segs, vec![(0, 51, 1.0), (51, 149, 1.5), (149, 300, 1.0)]);
    }

    #[test]
    fn test_build_coverage_multipliers_cnv_uses_cn() {
        // <CNV> with INFO/CN=4 on diploid → multiplier 2.0.
        let svs = vec![sv_variant_with_span(
            0,
            99,
            SvType::Cnv,
            Genotype::Homozygous,
            Some(4),
        )];
        // POS = 1 (1-based); span = 99 - 1 + 1 = 99. Anchor excluded, so the
        // affected range is [1, 0 + 99) = [1, 99). Was [0, 99).
        let segs = build_coverage_multipliers(&svs, 2, 200, false);
        assert_eq!(segs, vec![(0, 1, 1.0), (1, 99, 2.0), (99, 200, 1.0)]);
    }

    #[test]
    fn test_build_coverage_multipliers_cnv_without_cn_is_skipped() {
        let svs = vec![sv_variant_with_span(
            0,
            99,
            SvType::Cnv,
            Genotype::Homozygous,
            None,
        )];
        let segs = build_coverage_multipliers(&svs, 2, 200, false);
        assert_eq!(segs, vec![(0, 200, 1.0)]);
    }

    #[test]
    fn test_build_coverage_multipliers_ins_inv_unchanged() {
        let svs = vec![
            sv_variant_with_span(10, 19, SvType::Ins, Genotype::Heterozygous, None),
            sv_variant_with_span(50, 99, SvType::Inv, Genotype::Homozygous, None),
        ];
        let segs = build_coverage_multipliers(&svs, 2, 200, false);
        // Both should be skipped (multiplier == 1.0), leaving the default segment.
        assert_eq!(segs, vec![(0, 200, 1.0)]);
    }

    #[test]
    fn test_split_region_by_multipliers_intersects_correctly() {
        let segs = vec![(0usize, 100usize, 1.0f64), (100, 200, 0.5), (200, 300, 1.0)];
        // Region fully inside one segment.
        let r1 = split_region_by_multipliers(120, 150, &segs);
        assert_eq!(r1, vec![(120, 150, 0.5)]);
        // Region spanning two segments — split at the boundary.
        let r2 = split_region_by_multipliers(80, 180, &segs);
        assert_eq!(r2, vec![(80, 100, 1.0), (100, 180, 0.5)]);
        // Region outside any segment yields nothing (when segs are non-empty).
        let r3 = split_region_by_multipliers(400, 500, &segs);
        assert!(r3.is_empty());
        // Empty segments → return the whole region with multiplier 1.0.
        let r4 = split_region_by_multipliers(10, 20, &[]);
        assert_eq!(r4, vec![(10, 20, 1.0)]);
    }

    #[test]
    fn test_scale_coverage_rounds_and_clamps() {
        assert_eq!(scale_coverage(10, 0.0), 0);
        assert_eq!(scale_coverage(10, 0.5), 5);
        assert_eq!(scale_coverage(10, 1.5), 15);
        assert_eq!(scale_coverage(10, 2.0), 20);
        // Negative / non-finite → 0.
        assert_eq!(scale_coverage(10, -1.0), 0);
        assert_eq!(scale_coverage(10, f64::NAN), 0);
        assert_eq!(scale_coverage(10, f64::INFINITY), 0);
        // 0.49 * 10 = 4.9 → rounds to 5.
        assert_eq!(scale_coverage(10, 0.49), 5);
        // 0.04 * 10 = 0.4 → rounds to 0.
        assert_eq!(scale_coverage(10, 0.04), 0);
    }

    #[test]
    fn test_apply_rate_override() {
        let segs = vec![(0usize, 100usize, 0.001f64), (200, 400, 0.001)];

        // No overlap: override entirely before first segment
        let result = apply_rate_override(segs.clone(), 0, 0, 0.01);
        assert_eq!(result, segs);

        // No overlap: override entirely after last segment
        let result = apply_rate_override(segs.clone(), 500, 600, 0.01);
        assert_eq!(result, segs);

        // Partial overlap at start of first segment: [0,50) gets new rate, [50,100) keeps old
        let result = apply_rate_override(segs.clone(), 0, 50, 0.01);
        assert_eq!(
            result,
            vec![(0, 50, 0.01), (50, 100, 0.001), (200, 400, 0.001)]
        );

        // Partial overlap at end of first segment: [0,80) keeps old, [80,100) gets new rate
        let result = apply_rate_override(segs.clone(), 80, 150, 0.01);
        assert_eq!(
            result,
            vec![(0, 80, 0.001), (80, 100, 0.01), (200, 400, 0.001)]
        );

        // Full containment of first segment: entire [0,100) replaced
        let result = apply_rate_override(segs.clone(), 0, 100, 0.02);
        assert_eq!(result, vec![(0, 100, 0.02), (200, 400, 0.001)]);

        // Override spanning both segments (gap between them is unaffected)
        let result = apply_rate_override(segs.clone(), 50, 300, 0.05);
        assert_eq!(
            result,
            vec![
                (0, 50, 0.001),
                (50, 100, 0.05),
                (200, 300, 0.05),
                (300, 400, 0.001)
            ]
        );
    }

    #[test]
    fn test_rate_at() {
        let segs = vec![(10usize, 50usize, 0.001f64), (100, 200, 0.005)];

        // Inside first segment
        assert_eq!(rate_at(&segs, 25), 0.001);
        // At start of first segment (inclusive)
        assert_eq!(rate_at(&segs, 10), 0.001);
        // At end of first segment (exclusive — gap)
        assert_eq!(rate_at(&segs, 50), 0.0);
        // In gap between segments
        assert_eq!(rate_at(&segs, 75), 0.0);
        // Before all segments
        assert_eq!(rate_at(&segs, 0), 0.0);
        // Inside second segment
        assert_eq!(rate_at(&segs, 150), 0.005);
        // At end of second segment (exclusive)
        assert_eq!(rate_at(&segs, 200), 0.0);
    }

    #[test]
    fn test_exclude_positions() {
        let segs = vec![(0usize, 100usize, 0.001f64), (200, 400, 0.001)];

        // Empty excluded list — segments unchanged
        assert_eq!(exclude_positions(segs.clone(), &[]), segs);

        // Exclude middle of first segment — splits into two
        let result = exclude_positions(segs.clone(), &[50]);
        assert_eq!(
            result,
            vec![(0, 50, 0.001), (51, 100, 0.001), (200, 400, 0.001)]
        );

        // Exclude start of segment — trims left boundary
        let result = exclude_positions(segs.clone(), &[0]);
        assert_eq!(result, vec![(1, 100, 0.001), (200, 400, 0.001)]);

        // Exclude last position of segment — trims right boundary
        let result = exclude_positions(segs.clone(), &[99]);
        assert_eq!(result, vec![(0, 99, 0.001), (200, 400, 0.001)]);

        // Exclude position in gap — no change to segments
        let result = exclude_positions(segs.clone(), &[150]);
        assert_eq!(result, segs);

        // Exclude multiple positions across both segments
        let result = exclude_positions(segs.clone(), &[50, 250]);
        assert_eq!(
            result,
            vec![
                (0, 50, 0.001),
                (51, 100, 0.001),
                (200, 250, 0.001),
                (251, 400, 0.001),
            ]
        );
    }
}
