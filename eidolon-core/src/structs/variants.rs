//! These enums and structs help us keep track of variants added to the reference
//! It stores all the necessary data to write a variant out.
use crate::structs::nucleotides::{Nucleotide, is_acgtn};
use log::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum VariantError {
    #[error("ERROR: Malformed indel variant!")]
    MalformedIndel,
    #[error("ERROR: Malformed SNP variant!")]
    MalformedSnp,
    #[error("ERROR: Malformed variant ref!")]
    MalformedRef,
    #[error("ERROR: Malformed variant alt!")]
    MalformedAlt,
    #[error("Error Generating the genotype string!")]
    GenoStringError,
    #[error("Invalid Vcf Input: {0}")]
    InvalidVcf(String),
    #[error("Complex variants are not yet supported.")]
    ComplexNotSupported,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Genotype {
    Heterozygous,
    Homozygous,
}

/// Where a variant came from in the current `gen-reads` run.
///
/// Emitted as `INFO/NEAT_PROVENANCE` in the golden VCF so downstream
/// merges (e.g. `tools/cancer_simulate.sh`'s normal/tumor truth merge)
/// can resolve cancer-specific origin labels (`germline` / `somatic` /
/// `shared`) without `gen-reads` itself needing any cancer vocabulary.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub enum Provenance {
    /// Variant was sampled de novo by the simulator.
    Denovo,
    /// Variant was carried through from the user's `input_vcf:` file.
    InputVcf,
    /// Variant was supplied via a reproductive somatic VCF (#405). Distinct from
    /// `InputVcf` (germline) so a tumor/normal merge resolves it to `somatic`, not
    /// `shared`, even though — like `InputVcf` — it came from a file.
    SomaticVcf,
}

impl Provenance {
    /// Lowercase string used in `INFO/NEAT_PROVENANCE=` values.
    pub fn as_str(&self) -> &'static str {
        match self {
            Provenance::Denovo => "denovo",
            Provenance::InputVcf => "input",
            Provenance::SomaticVcf => "somatic_input",
        }
    }
}

fn gt_from_str(input: &str) -> Genotype {
    let mut result: Vec<&str> = input.split("/").collect();
    if result.len() == 1 {
        // Split may have failed, trying a different delimiter
        result = input.split("|").collect();
    }
    let mut found_zero = false;
    for item in result {
        if item == "." {
            continue;
        }
        if item.parse::<u64>().unwrap_or(1) == 0 {
            found_zero = true;
            break;
        }
    }
    if found_zero {
        Genotype::Heterozygous
    } else {
        Genotype::Homozygous
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub enum VariantType {
    SNP,
    Insertion,
    Deletion,
    Complex,
    BND,
}

impl From<usize> for VariantType {
    fn from(i: usize) -> Self {
        match i {
            0 => Self::SNP,
            1 => Self::Insertion,
            2 => Self::Deletion,
            3 => Self::Complex,
            4 => Self::BND,
            _ => panic!("Index out of range!"),
        }
    }
}

/// Classifies a symbolic / structural ALT (VCF 4.2 §1.4.1 — §1.4.5).
///
/// `Unknown` is the catch-all for symbolic ALTs we recognize as symbolic but
/// can't classify — e.g. a tag like `<NOVEL_THING>` that isn't in the spec. We
/// still round-trip the raw ALT string verbatim from `SvData::raw_alt`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum SvType {
    Del,
    Dup,
    Cnv,
    Ins,
    Inv,
    Bnd,
    Unknown,
}

impl SvType {
    /// Parse a VCF symbolic / breakend ALT string into a coarse type tag.
    /// Returns `None` for literal-base ALTs (callers should not invoke this
    /// on `AlternateType::Literal` payloads).
    ///
    /// Recognizes:
    ///   - `<DEL>`, `<DEL:...>` → Del
    ///   - `<DUP>`, `<DUP:TANDEM>`, `<DUP:...>` → Dup
    ///   - `<CNV>`, `<CNV:...>` → Cnv
    ///   - `<INS>`, `<INS:ME:ALU>`, `<INS:...>` → Ins
    ///   - `<INV>` → Inv
    ///   - Breakend notation (`G]17:198982]`, `]13:123456]T`, etc.) → Bnd
    ///   - Other `<TAG>` forms → Unknown
    pub fn from_alt_string(alt: &str) -> Option<Self> {
        if alt.is_empty() {
            return None;
        }
        if alt.starts_with('<') && alt.ends_with('>') {
            let body = &alt[1..alt.len() - 1];
            let head = body.split(':').next().unwrap_or(body);
            return Some(match head {
                "DEL" => SvType::Del,
                "DUP" => SvType::Dup,
                "CNV" => SvType::Cnv,
                "INS" => SvType::Ins,
                "INV" => SvType::Inv,
                "BND" => SvType::Bnd,
                _ => SvType::Unknown,
            });
        }
        // Breakend notation contains `[` or `]` and a `:`-separated mate locus.
        if (alt.contains('[') || alt.contains(']')) && alt.contains(':') {
            return Some(SvType::Bnd);
        }
        // Single breakends (VCF 4.2 §1.4.2)
        if alt.contains('.') && !alt.starts_with('<') {
            // Check if it matches t. or .t where t is a sequence of bases
            if alt.ends_with('.') {
                let t = &alt[..alt.len() - 1];
                if !t.is_empty() && is_acgtn(t) {
                    return Some(SvType::Bnd);
                }
            } else if alt.starts_with('.') {
                let t = &alt[1..];
                if !t.is_empty() && is_acgtn(t) {
                    return Some(SvType::Bnd);
                }
            }
        }
        None
    }
}

/// Parsed payload for a symbolic / structural ALT.
///
/// Stores the original ALT string for verbatim round-trip, plus optional
/// INFO-derived fields (`END=`, `SVLEN=`, `CN=`). `Variant::from_file`
/// populates the optional fields from the VCF INFO column via `parse_sv_info`;
/// other constructors (e.g. `SvData::new`) start them as `None` and let the
/// caller fill in whatever it has.
#[derive(Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub struct SvData {
    /// The original ALT field, verbatim — used to round-trip back to the
    /// output VCF without re-synthesizing the angle-bracket notation.
    pub raw_alt: String,
    pub sv_type: SvType,
    /// 1-based inclusive end coordinate from `INFO/END`, if present.
    pub end: Option<usize>,
    /// `INFO/SVLEN` if present. Signed per VCF 4.2: negative for deletions,
    /// positive for insertions/duplications.
    pub svlen: Option<i64>,
    /// `INFO/CN` if present. Total copy number for `<CNV>` / `<DUP>` records.
    pub copy_number: Option<u32>,
    /// Mate contig for breakends (BND).
    pub mate_contig: Option<String>,
    /// 1-based mate position for breakends (BND).
    pub mate_pos: Option<usize>,
    /// For BND: true if the breakend is after the REF base (e.g. t[p[ or t]p]).
    pub bnd_join_after: bool,
    /// For BND: true if the mate piece extends to the right ([p[).
    pub bnd_mate_extends_right: bool,
}

impl SvData {
    /// Construct an `SvData` with only the raw ALT and type tag populated.
    /// Optional INFO fields default to `None` — set them on the returned
    /// value if available.
    pub fn new(raw_alt: impl Into<String>, sv_type: SvType) -> Self {
        SvData {
            raw_alt: raw_alt.into(),
            sv_type,
            end: None,
            svlen: None,
            copy_number: None,
            mate_contig: None,
            mate_pos: None,
            bnd_join_after: false,
            bnd_mate_extends_right: false,
        }
    }

    /// Compute the reference span (in bases) of this SV, given the variant's
    /// 1-based start `location`. Prefers `END` (span = end - location + 1);
    /// falls back to `|SVLEN|`. Returns `None` if neither is set, or if `END`
    /// is present but `< location` — defense-in-depth against an upstream
    /// caller that constructed `SvData` directly with a malformed `END`.
    /// (`Variant::from_file` already drops invalid `END` values at parse
    /// time with a warning, so this branch normally never fires.)
    pub fn span(&self, location: usize) -> Option<usize> {
        if let Some(end) = self.end {
            if end < location {
                return None;
            }
            return Some(end - location + 1);
        }
        if let Some(svlen) = self.svlen {
            return Some(svlen.unsigned_abs() as usize);
        }
        None
    }

    /// Compute the total length of the variant event. For non-insertions, this
    /// is identical to the reference span. For insertions, it is the length
    /// of the inserted sequence (derived from `SVLEN`).
    pub fn event_length(&self, location: usize) -> Option<usize> {
        match self.sv_type {
            SvType::Ins => self.svlen.map(|l| l.unsigned_abs() as usize),
            _ => self.span(location),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub enum AlternateType {
    Literal(Vec<Nucleotide>),
    Symbolic(SvData),
}

impl AlternateType {
    /// Returns the literal base sequence if this is a `Literal` variant,
    /// or `None` for symbolic / structural ALTs.
    pub fn as_literal(&self) -> Option<&[Nucleotide]> {
        match self {
            AlternateType::Literal(vector) => Some(vector),
            AlternateType::Symbolic(_) => None,
        }
    }

    /// Returns the parsed symbolic payload if this is a `Symbolic` variant,
    /// or `None` for literal-base ALTs.
    pub fn as_symbolic(&self) -> Option<&SvData> {
        match self {
            AlternateType::Symbolic(sv) => Some(sv),
            AlternateType::Literal(_) => None,
        }
    }

    pub fn is_literal(&self) -> bool {
        matches!(self, AlternateType::Literal(_))
    }

    pub fn is_symbolic(&self) -> bool {
        matches!(self, AlternateType::Symbolic(_))
    }
}

// `allele_fraction: Option<f64>` means Variant can't derive `Eq`/`Hash` (f64 is
// only `PartialEq`). Neither was used — Variant is always a value in a Vec/HashMap,
// never a key or set member; compare_vcfs hashes a separate `VariantKey` projection.
#[derive(Debug, PartialEq, Clone)]
pub struct Variant {
    // This is a basic holder for a variant. There are several aspects of the variant that we
    // don't need to store, such as which ploid the variant appears on and the context sequence,
    // which we only need for the output files, we can generate at the time we need it, such as the
    // insertion sequence.
    //
    // This is the type of variant. Any new types will need to be added to VariantType
    // to be implemented here.
    pub variant_type: VariantType,
    // the location, in reference coordinates, of the start of this variant
    pub location: usize,
    // The reference allele of interest. This is either one base or several bases for a deletion.
    pub reference: Vec<Nucleotide>,
    // The alternate allele of interest. This is either one base or several bases for an indel, or
    // it is a key string from the vcf spec indicating a complex variant.
    pub alternate: AlternateType,
    // the genotype string is for writing in the vcf.
    pub genotype_str: String,
    // The genotype is either heterozygous (not on all alleles) or homozygous (on all alleles)
    pub genotype: Genotype,
    // Optional per-variant alt-allele fraction in (0, 1], parsed from an input VCF
    // (INFO/AF, else FORMAT/AD). When present, read generation emits the alt allele on
    // this fraction of overlapping reads instead of the {homozygous=1.0, het=0.5} default —
    // letting an input VCF drive a continuous AF spectrum (e.g. pooled/somatic data, #398).
    // `None` for model-generated (de novo) variants, which keep the Genotype-based behavior.
    pub allele_fraction: Option<f64>,

    // VCF files specific. These are variables specific to VCF files read as input,
    // either for variant inputs or for generating models.

    // Included id, if present
    pub id: Option<String>,
    // Included quality score
    pub quality_score: Option<usize>,
    // Included filter, if present
    pub filter: Option<String>,
    // Included Info, if present
    pub info: Option<String>,
    // Other included information, if present
    pub format: Vec<String>,
    // Sample, if present, must be same length as format
    pub sample: Vec<String>,
    // Where this variant came from in the current gen-reads run.
    // Emitted as INFO/NEAT_PROVENANCE in the golden VCF.
    pub provenance: Provenance,
}

/// Clamps a parsed allele fraction to the valid `(0, 1]` domain used by read
/// generation. Non-finite or non-positive values become `None` (fall back to the
/// Genotype-based fraction); values above 1.0 are clamped to 1.0.
fn sanitize_fraction(f: f64) -> Option<f64> {
    if !f.is_finite() || f <= 0.0 {
        None
    } else if f > 1.0 {
        Some(1.0)
    } else {
        Some(f)
    }
}

/// Parses an optional alt-allele fraction from an input-VCF record. Prefers
/// `INFO/AF` (first value if comma-separated); falls back to `FORMAT/AD`
/// (sum of alt depths / total depth). Returns `None` when neither is present or
/// usable — the variant then keeps its Genotype-based fraction. AD is the correct
/// source for pooled data, where GT is not diploid-meaningful (see #398).
fn parse_allele_fraction(info_field: &str, format: &[String], sample: &[String]) -> Option<f64> {
    if let Some(af) = info_field
        .split(';')
        .find_map(|kv| kv.strip_prefix("AF="))
        .and_then(|v| v.split(',').next())
        .and_then(|v| v.trim().parse::<f64>().ok())
    {
        return sanitize_fraction(af);
    }
    if let Some(pos) = format.iter().position(|f| f == "AD")
        && let Some(ad) = sample.get(pos)
    {
        let counts: Vec<f64> = ad
            .split(',')
            .filter_map(|x| x.trim().parse::<f64>().ok())
            .collect();
        if counts.len() >= 2 {
            let total: f64 = counts.iter().sum();
            let alt: f64 = counts[1..].iter().sum();
            if total > 0.0 {
                return sanitize_fraction(alt / total);
            }
        }
    }
    None
}

impl Variant {
    pub fn from_file(
        location: usize,
        id: &str,
        filter: &str,
        info_field: &str,
        vcf_reference: &str,
        vcf_alternate: &str,
        quality_score: usize,
        format: Vec<String>,
        sample: Vec<String>,
    ) -> Result<Self, VariantError> {
        let genotype_str = {
            if format.contains(&String::from("GT")) {
                trace!("Found genotype in vcf file");
                let gt_pos = format.clone().into_iter().position(|x| x == "GT").unwrap();
                sample[gt_pos].to_string()
            } else {
                return Err(VariantError::InvalidVcf(
                    "Missing GT field. Cannot process input VCF without GT in FORMAT field"
                        .to_string(),
                ));
            }
        };

        let (variant_type, alternate, reference) =
            parse_alt_payload(location, info_field, vcf_reference, vcf_alternate)?;
        let genotype = gt_from_str(&genotype_str);
        let allele_fraction = parse_allele_fraction(info_field, &format, &sample);
        Ok(Variant {
            variant_type,
            location,
            reference,
            alternate,
            genotype_str,
            genotype,
            allele_fraction,
            id: Some(id.to_string()),
            quality_score: Some(quality_score),
            filter: Some(filter.to_string()),
            info: Some(info_field.to_string()),
            format,
            sample,
            provenance: Provenance::InputVcf,
        })
    }

    /// Minimal constructor for callers that only need structural fields
    /// (variant_type, location, reference, alternate, genotype) — e.g.
    /// gen_mut_model, which trains from millions of records and never
    /// round-trips them back to VCF.
    ///
    /// Skips allocating the per-record String fields (`info`, `id`,
    /// `filter`, `genotype_str`, plus the FORMAT/SAMPLE vectors) that
    /// gen_mut_model never reads. On gnomAD-scale inputs the `info`
    /// column alone is hundreds of bytes to several KB per record, so
    /// dropping it cuts memory by ~2–5 GB on a full-chromosome fit.
    pub fn from_file_lean(
        location: usize,
        info_field: &str,
        vcf_reference: &str,
        vcf_alternate: &str,
        gt_str: &str,
    ) -> Result<Self, VariantError> {
        let (variant_type, alternate, reference) =
            parse_alt_payload(location, info_field, vcf_reference, vcf_alternate)?;
        let genotype = gt_from_str(gt_str);
        Ok(Variant {
            variant_type,
            location,
            reference,
            alternate,
            genotype_str: String::new(),
            genotype,
            // gen_mut_model (the only from_file_lean caller) trains on positions/types
            // and never reads allele_fraction, so skip the parse.
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::InputVcf,
        })
    }

    pub fn new(
        variant_type: VariantType,
        location: usize,
        reference: &Vec<Nucleotide>,
        alternate: &Vec<Nucleotide>,
        genotype: &mut Vec<usize>,
    ) -> Result<Self, VariantError> {
        // Check for empty vecs before any index access
        if reference.is_empty() {
            return Err(VariantError::MalformedRef);
        }
        if alternate.is_empty() {
            return Err(VariantError::MalformedAlt);
        }

        match variant_type {
            VariantType::Insertion | VariantType::Deletion => {
                if (reference.len() == alternate.len()) || reference[0] != alternate[0] {
                    return Err(VariantError::MalformedIndel);
                }
            }
            VariantType::SNP => {
                if reference.len() != alternate.len() {
                    return Err(VariantError::MalformedSnp);
                }
            }
            _ => {
                return Err(VariantError::ComplexNotSupported);
            }
        }

        let mut genotype_label = Genotype::Homozygous;
        // Homozygous would mean all alleles have the variant, so all are 1's
        // the presence of a single 0 indicates heterozygous.
        // This is a little shakey biologically, but we're going with this.
        if genotype.contains(&0) {
            genotype_label = Genotype::Heterozygous;
        }

        Ok(Variant {
            variant_type,
            location,
            reference: reference.clone(),
            alternate: AlternateType::Literal(alternate.clone()),
            genotype_str: genotype_to_string(genotype)?,
            genotype: genotype_label,
            // De novo variants use the Genotype-based fraction (Het/Hom), not an input AF.
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        })
    }

    pub fn get_loc(&self) -> Result<usize, VariantError> {
        Ok(self.location)
    }

    /// Allele **dosage** as a fraction: alt-copy count over total called copies,
    /// derived from the N-wide `genotype_str` (`0/1` → 0.5, `1/1` → 1.0,
    /// `0/0/0/1` → 0.25). Any non-zero allele index counts as alt; missing (`.`)
    /// alleles are excluded from the total. This is the biological "what fraction
    /// of this cell's homologous copies carry the variant" and, for a diploid,
    /// reproduces the read sampler's het-0.5 / hom-1.0 default.
    ///
    /// It is the composable **base** for `allele_fraction`: a subclonal CCF (#405)
    /// or other cellular-fraction factor multiplies onto it, and polyploid dosage
    /// (#266/#267) flows in automatically once `genotype_str` carries a real
    /// per-copy spread. Falls back to the coarse `Genotype` label (Hom → 1.0,
    /// Het → 0.5) if `genotype_str` has no parseable alleles.
    pub fn dosage_fraction(&self) -> f64 {
        let (mut alt, mut total) = (0usize, 0usize);
        for allele in self.genotype_str.split(['/', '|']) {
            match allele.trim() {
                "" | "." => {}
                "0" => total += 1,
                _ => {
                    alt += 1;
                    total += 1;
                }
            }
        }
        if total == 0 {
            return match self.genotype {
                Genotype::Homozygous => 1.0,
                Genotype::Heterozygous => 0.5,
            };
        }
        alt as f64 / total as f64
    }
}

/// Shared core of `Variant::from_file` and `Variant::from_file_lean`: parse
/// REF into a `Vec<Nucleotide>`, classify the ALT, and (for symbolic ALTs)
/// lift `END` / `SVLEN` / `CN` / `SVTYPE` from the INFO column into `SvData`.
/// Both constructors handle the trailing string-field population themselves.
fn parse_alt_payload(
    location: usize,
    info_field: &str,
    vcf_reference: &str,
    vcf_alternate: &str,
) -> Result<(VariantType, AlternateType, Vec<Nucleotide>), VariantError> {
    let mut reference = Vec::new();
    for char in vcf_reference.chars() {
        reference.push(Nucleotide::from(char))
    }
    let (variant_type, alternate) = match parse_alternate(vcf_reference, vcf_alternate)? {
        ParsedAlt::Literal(vt) => {
            let alt_bases: Vec<Nucleotide> = vcf_alternate.chars().map(Nucleotide::from).collect();
            (vt, AlternateType::Literal(alt_bases))
        }
        ParsedAlt::Symbolic { sv_type, raw_alt } => {
            let mut sv_data = SvData::new(raw_alt, sv_type);
            let parsed_info = parse_sv_info(info_field);
            // If SVTYPE is present and disagrees with the ALT tag, trust
            // the ALT (it's the more canonical signal) but warn so the
            // user can clean up their VCF.
            if let Some(svtype_str) = &parsed_info.svtype
                && !sv_type_matches_svtype(sv_type, svtype_str)
            {
                warn!(
                    "INFO/SVTYPE={svtype_str} disagrees with ALT {} at position {location} — trusting ALT",
                    sv_data.raw_alt
                );
            }
            // Drop an INFO/END that's before POS — keeping it would let
            // span() silently produce a bogus 1-base modulation range for
            // <DUP>/<CNV>/<INV> (the <DEL> path masks the bug via an
            // empty-range check downstream). Fall back to SVLEN if the
            // user supplied that too, otherwise the record passes through
            // with no depth modulation and gen_reads logs a "no span" warn.
            sv_data.end = match parsed_info.end {
                Some(end) if end < location => {
                    warn!(
                        "INFO/END={end} is before POS={location} for ALT {} — ignoring END",
                        sv_data.raw_alt
                    );
                    None
                }
                other => other,
            };
            sv_data.svlen = parsed_info.svlen;
            sv_data.copy_number = parsed_info.copy_number;
            if matches!(sv_type, SvType::Bnd) {
                let (m_contig, m_pos, j_after, m_extends_right) = parse_bnd_alt(&sv_data.raw_alt);
                sv_data.mate_contig = m_contig;
                sv_data.mate_pos = m_pos;
                sv_data.bnd_join_after = j_after;
                sv_data.bnd_mate_extends_right = m_extends_right;
            }
            // BND (breakend) variants are classified with VariantType::BND
            // to allow downstream code to specifically handle them.
            // Other symbolic SVs share the `Complex` variant_type with literal
            // multi-base substitutions; downstream code distinguishes them
            // by the `AlternateType` enum (literal vs symbolic), not by
            // variant_type.
            let vt = if matches!(sv_type, SvType::Bnd) {
                VariantType::BND
            } else {
                VariantType::Complex
            };
            (vt, AlternateType::Symbolic(sv_data))
        }
    };
    Ok((variant_type, alternate, reference))
}

/// Classification of a parsed VCF ALT field.
///
/// Returned by [`parse_alternate`] so callers can branch on whether the
/// record should be built as an `AlternateType::Literal` (with base content)
/// or an `AlternateType::Symbolic` (with parsed SV metadata).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParsedAlt {
    /// Literal-base ALT (SNP / Insertion / Deletion / Complex multi-base).
    /// The caller still constructs the `Vec<Nucleotide>` from the raw string.
    Literal(VariantType),
    /// Symbolic / structural / breakend ALT. The `raw_alt` field is the
    /// original VCF ALT string, preserved verbatim for round-trip output.
    Symbolic { sv_type: SvType, raw_alt: String },
}

fn parse_alternate(reference: &str, alt: &str) -> Result<ParsedAlt, VariantError> {
    if reference.is_empty() {
        return Err(VariantError::InvalidVcf(
            "REF field in vcf with empty ref. Invalid VCF.".to_string(),
        ));
    }
    if alt.is_empty() {
        return Err(VariantError::InvalidVcf(
            "ALT field in vcf with empty alt. Invalid VCF.".to_string(),
        ));
    }
    if !is_acgtn(reference) {
        return Err(VariantError::InvalidVcf(format!(
            "REF must contain only A/C/G/T/N; got {reference}"
        )));
    }
    // Symbolic / structural / breakend ALTs (e.g. <DUP:TANDEM>, <DEL>,
    // <INS:ME:ALU>, `G]17:198982]`) classify into ParsedAlt::Symbolic with
    // the raw ALT preserved verbatim. gen_reads uses the resulting SvData
    // for per-region depth modulation and round-trips the raw ALT to the
    // output VCF; eidolon does not yet generate symbolic SVs de novo.
    if !is_acgtn(alt) {
        return match SvType::from_alt_string(alt) {
            Some(sv_type) => Ok(ParsedAlt::Symbolic {
                sv_type,
                raw_alt: alt.to_string(),
            }),
            None => Err(VariantError::InvalidVcf(format!(
                "ALT is neither literal bases nor a recognized symbolic / breakend form: {alt}"
            ))),
        };
    }
    Ok(ParsedAlt::Literal(match (reference.len(), alt.len()) {
        (1, 1) => VariantType::SNP,
        (1, _) => VariantType::Insertion,
        (_, 1) => VariantType::Deletion,
        (_, _) => VariantType::Complex,
    }))
}

/// Parsed fields lifted from the VCF INFO column for a structural variant
/// record. All fields are optional — callers should treat absence as
/// "not specified" rather than as a hard error.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct ParsedSvInfo {
    pub end: Option<usize>,
    pub svlen: Option<i64>,
    pub copy_number: Option<u32>,
    /// The `SVTYPE=` value verbatim (e.g. `DEL`, `DUP`, `CNV`).
    /// `Variant::from_file` uses this only to log a warning if it disagrees
    /// with the symbolic ALT tag — the ALT is the source of truth.
    pub svtype: Option<String>,
}

/// Scan an INFO column for the structural-variant-related fields
/// (`END`, `SVLEN`, `CN`, `SVTYPE`). The INFO column is a `;`-separated
/// list of `KEY=VALUE` (or bare-flag) pairs per VCF 4.2 §1.6.1.
///
/// Malformed values for a known key are logged at warn level and dropped
/// (the field stays `None`) — a bad `SVLEN=` shouldn't lose the entire
/// variant. Unknown keys are silently ignored.
pub fn parse_sv_info(info_field: &str) -> ParsedSvInfo {
    let mut out = ParsedSvInfo::default();
    if info_field == "." || info_field.is_empty() {
        return out;
    }
    for entry in info_field.split(';') {
        let (key, value) = match entry.split_once('=') {
            Some(kv) => kv,
            None => continue, // bare flag like IMPRECISE
        };
        match key {
            "END" => match value.parse::<usize>() {
                Ok(n) => out.end = Some(n),
                Err(e) => warn!("Skipping malformed INFO/END={value}: {e}"),
            },
            "SVLEN" => match value.parse::<i64>() {
                Ok(n) => out.svlen = Some(n),
                Err(e) => warn!("Skipping malformed INFO/SVLEN={value}: {e}"),
            },
            "CN" => match value.parse::<u32>() {
                Ok(n) => out.copy_number = Some(n),
                Err(e) => warn!("Skipping malformed INFO/CN={value}: {e}"),
            },
            "SVTYPE" => out.svtype = Some(value.to_string()),
            _ => {}
        }
    }
    out
}

/// Cross-check a [`SvType`] (derived from the ALT tag) against the textual
/// `SVTYPE=` value (lifted from INFO). Returns `true` if they describe the
/// same kind of SV, `false` if they disagree. An `Unknown` ALT tag is
/// treated as compatible with anything (we can't contradict a tag we
/// didn't recognize).
fn sv_type_matches_svtype(sv_type: SvType, svtype_str: &str) -> bool {
    if matches!(sv_type, SvType::Unknown) {
        return true;
    }
    let expected = match sv_type {
        SvType::Del => "DEL",
        SvType::Dup => "DUP",
        SvType::Cnv => "CNV",
        SvType::Ins => "INS",
        SvType::Inv => "INV",
        SvType::Bnd => "BND",
        SvType::Unknown => return true,
    };
    svtype_str.eq_ignore_ascii_case(expected)
}

fn parse_bnd_alt(alt: &str) -> (Option<String>, Option<usize>, bool, bool) {
    // VCF 4.2 breakend notation (Section 1.4.2)
    // 4 possible forms involving brackets:
    // 1. t[p[  2. t]p]  3. ]p]t  4. [p[t
    // where p is "contig:pos"
    let bracket_start = alt.find('[').or_else(|| alt.find(']'));
    let bracket_end = alt.rfind('[').or_else(|| alt.rfind(']'));

    if let (Some(s), Some(e)) = (bracket_start, bracket_end) {
        if s < e {
            let p = &alt[s + 1..e];
            let join_after = s > 0;
            let bracket = alt.chars().nth(s).unwrap();
            let mate_extends_right = bracket == '[';

            if let Some((contig, pos_str)) = p.split_once(':') {
                if let Ok(pos) = pos_str.parse::<usize>() {
                    return (
                        Some(contig.to_string()),
                        Some(pos),
                        join_after,
                        mate_extends_right,
                    );
                }
            }
        }
    }
    (None, None, false, false)
}

fn genotype_to_string(genotype: &mut Vec<usize>) -> Result<String, VariantError> {
    // Converts a vector of 0s and 1s representing genotype to a standard
    // vcf genotype string.
    //
    // In order to make the vcf look a little more realistic, we shuffle the genotypes, that way
    // it doesn't look like one ploid got all the mutations. I actually don't know if this is
    // realistic or not, but that's what we're doing.
    let geno_str = genotype
        .iter()
        .map(|x| x.to_string() + "/")
        .collect::<String>();
    match geno_str.strip_suffix("/") {
        Some(str) => Ok(str.to_string()),
        None => Err(VariantError::GenoStringError),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::variants::VariantType::{Complex, Deletion, Insertion, SNP};

    #[test]
    fn test_variant_creation() {
        let variant = Variant {
            // Not actually a valid variant, but just testing the constructor
            variant_type: SNP,
            location: 222,
            reference: vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            alternate: AlternateType::Literal(vec![Nucleotide::A]),
            genotype_str: "1/1/1".to_string(),
            genotype: Genotype::Homozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        assert_eq!(variant.variant_type, SNP);
        assert_eq!(
            variant.reference,
            vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G]
        );
        assert_eq!(variant.genotype, Genotype::Homozygous);
    }

    #[test]
    fn test_variant_new() {
        let variant = Variant::new(
            Insertion,
            10,
            &vec![Nucleotide::A],
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &mut vec![0, 1],
        )
        .unwrap();
        assert_eq!(variant.genotype, Genotype::Heterozygous);
        // Variant::new is the de-novo constructor — gen-reads samples
        // germline and somatic variants through this path, so the writer
        // can rely on Denovo as the default provenance for unannotated
        // outputs.
        assert_eq!(variant.provenance, Provenance::Denovo);
    }

    #[test]
    fn dosage_fraction_from_genotype_string() {
        let snp = |gt: &mut Vec<usize>| {
            Variant::new(SNP, 10, &vec![Nucleotide::A], &vec![Nucleotide::T], gt).unwrap()
        };
        // Diploid: het → 0.5, hom → 1.0 (reproduces the read-sampler default).
        assert_eq!(snp(&mut vec![0, 1]).dosage_fraction(), 0.5);
        assert_eq!(snp(&mut vec![1, 1]).dosage_fraction(), 1.0);
        // Polyploid dosage classes fall straight out of the N-wide GT string.
        assert_eq!(snp(&mut vec![0, 0, 0, 1]).dosage_fraction(), 0.25);
        assert_eq!(snp(&mut vec![0, 1, 1, 1]).dosage_fraction(), 0.75);
    }

    #[test]
    fn dosage_fraction_handles_missing_and_phased_alleles() {
        let v = Variant {
            variant_type: SNP,
            location: 1,
            reference: vec![Nucleotide::A],
            alternate: AlternateType::Literal(vec![Nucleotide::T]),
            // Phased with a missing allele: alt=1 over called=2 → 0.5.
            genotype_str: "0|1|.".to_string(),
            genotype: Genotype::Heterozygous,
            allele_fraction: None,
            id: None,
            quality_score: None,
            filter: None,
            info: None,
            format: Vec::new(),
            sample: Vec::new(),
            provenance: Provenance::Denovo,
        };
        assert_eq!(v.dosage_fraction(), 0.5);
    }

    #[test]
    fn test_variant_from_file_carries_input_provenance() {
        // Variants read from a user-supplied input_vcf must be tagged
        // InputVcf so the cancer-simulator merge step can distinguish
        // them from de-novo somatic draws on the tumor pass.
        let v = Variant::from_file(
            99,
            "rs1",
            "PASS",
            ".",
            "A",
            "G",
            60,
            vec!["GT".to_string()],
            vec!["0/1".to_string()],
        )
        .unwrap();
        assert_eq!(v.provenance, Provenance::InputVcf);
    }

    #[test]
    fn test_deletion_new() {
        let variant = Variant::new(
            Deletion,
            22,
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &vec![Nucleotide::A],
            &mut vec![1, 1],
        )
        .unwrap();
        assert_eq!(variant.genotype, Genotype::Homozygous)
    }

    #[test]
    #[should_panic(expected = "MalformedSnp")]
    fn test_bad_variant_creation() {
        Variant::new(
            SNP,
            22,
            &vec![Nucleotide::A, Nucleotide::C, Nucleotide::T, Nucleotide::G],
            &vec![Nucleotide::A],
            &mut vec![1, 1, 1],
        )
        .unwrap();
    }

    #[test]
    fn test_genotype_to_string() {
        let mut genotype = vec![0, 1, 0];
        assert_eq!(
            String::from("0/1/0"),
            genotype_to_string(&mut genotype).unwrap()
        );
    }

    #[test]
    fn test_empty_reference_returns_malformed_ref() {
        let result = Variant::new(SNP, 0, &vec![], &vec![Nucleotide::G], &mut vec![1, 1]);
        assert!(matches!(result, Err(VariantError::MalformedRef)));
    }

    #[test]
    fn test_empty_alternate_returns_malformed_alt() {
        let result = Variant::new(SNP, 0, &vec![Nucleotide::A], &vec![], &mut vec![1, 1]);
        assert!(matches!(result, Err(VariantError::MalformedAlt)));
    }

    fn assert_literal(parsed: ParsedAlt, expected: VariantType) {
        match parsed {
            ParsedAlt::Literal(vt) => assert_eq!(vt, expected),
            ParsedAlt::Symbolic { sv_type, raw_alt } => {
                panic!("expected Literal({expected:?}), got Symbolic({sv_type:?}, {raw_alt})")
            }
        }
    }

    #[test]
    fn parse_alternate_snp() {
        assert_literal(parse_alternate("A", "C").unwrap(), SNP);
    }

    #[test]
    fn parse_alternate_insertion() {
        assert_literal(parse_alternate("A", "ACGT").unwrap(), Insertion);
    }

    #[test]
    fn parse_alternate_deletion() {
        assert_literal(parse_alternate("ACGT", "A").unwrap(), Deletion);
    }

    #[test]
    fn parse_alternate_complex_multibase_substitution() {
        assert_literal(parse_alternate("AC", "GT").unwrap(), Complex);
    }

    #[test]
    fn parse_alternate_complex_unequal_length() {
        // Both > 1 but different lengths is still Complex (e.g., 2bp ref → 3bp alt).
        assert_literal(parse_alternate("AC", "GTA").unwrap(), Complex);
    }

    #[test]
    fn parse_alternate_accepts_lowercase_and_n() {
        // Soft-masked bases and N should still classify as literal-base variants.
        assert_literal(parse_alternate("a", "n").unwrap(), SNP);
        assert_literal(parse_alternate("A", "acgN").unwrap(), Insertion);
    }

    #[test]
    fn parse_alternate_empty_reference_is_invalid_vcf() {
        let err = parse_alternate("", "A").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("REF")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }

    #[test]
    fn parse_alternate_empty_alt_is_invalid_vcf() {
        let err = parse_alternate("A", "").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("ALT")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }

    #[test]
    fn parse_alternate_symbolic_alts_round_trip_with_correct_type() {
        // VCF 4.2 symbolic ALTs classify into ParsedAlt::Symbolic with the
        // correct SvType derived from the tag. The raw_alt string is preserved
        // verbatim so the writer can round-trip the original.
        let cases = [
            ("<DUP:TANDEM>", SvType::Dup),
            ("<DEL>", SvType::Del),
            ("<INS:ME:ALU>", SvType::Ins),
            ("<CNV>", SvType::Cnv),
            ("<INV>", SvType::Inv),
        ];
        for (alt, expected) in cases {
            match parse_alternate("A", alt).unwrap() {
                ParsedAlt::Symbolic { sv_type, raw_alt } => {
                    assert_eq!(sv_type, expected, "wrong SvType for {alt}");
                    assert_eq!(raw_alt, alt, "raw_alt should be verbatim");
                }
                ParsedAlt::Literal(vt) => {
                    panic!("expected Symbolic for {alt}, got Literal({vt:?})")
                }
            }
        }
    }

    #[test]
    fn parse_alternate_breakend_alt_classifies_as_bnd() {
        // VCF 4.2 breakend notation (Section 1.4.2)
        // 4 possible forms:
        // 1. t[p[ (piece extending to the right of p is joined after t)
        // 2. t]p] (piece extending to the left of p is joined after t)
        // 3. ]p]t (piece extending to the left of p is joined before t)
        // 4. [p[t (piece extending to the right of p is joined before t)

        let cases = vec![
            ("G", "G]17:198982]", "t]p]"),
            ("G", "G[17:198982[", "t[p["),
            ("A", "]13:123456]A", "]p]t"),
            ("A", "[13:123456[A", "[p[t"),
            ("GTC", "GTC[2:1000[", "multi-base t [p["),
            ("G", "G.", "single breakend t."),
            ("A", ".A", "single breakend .t"),
            ("A", "]1:100]ATG", "BND with insertion ]p]t"),
            ("A", "ATG[1:100[", "BND with insertion t[p["),
        ];

        for (ref_base, alt, desc) in cases {
            match parse_alternate(ref_base, alt).expect(desc) {
                ParsedAlt::Symbolic { sv_type, raw_alt } => {
                    assert_eq!(sv_type, SvType::Bnd, "form {} should be Bnd", desc);
                    assert_eq!(raw_alt, alt, "form {} raw_alt mismatch", desc);
                }
                ParsedAlt::Literal(vt) => {
                    panic!(
                        "expected Symbolic(Bnd) for form {}, got Literal({vt:?})",
                        desc
                    )
                }
            }

            // Check if VariantType is also BND when fully parsed
            let (vt, _, _) = parse_alt_payload(100, ".", ref_base, alt).expect(desc);
            assert_eq!(vt, VariantType::BND, "form {} VariantType mismatch", desc);
        }
    }

    #[test]
    fn parse_alternate_unrecognized_nonliteral_alt_is_invalid_vcf() {
        // A non-IUPAC ALT that's neither symbolic nor breakend should still
        // be rejected — e.g. a malformed bracket-only string.
        let err = parse_alternate("A", "<broken").unwrap_err();
        assert!(matches!(err, VariantError::InvalidVcf(_)));
    }

    #[test]
    fn parse_bnd_alt_extracts_mate_info() {
        assert_eq!(
            parse_bnd_alt("G]17:198982]"),
            (Some("17".to_string()), Some(198982), true, false)
        );
        assert_eq!(
            parse_bnd_alt("[chrX:123[A"),
            (Some("chrX".to_string()), Some(123), false, true)
        );
        assert_eq!(
            parse_bnd_alt("]2:500]T"),
            (Some("2".to_string()), Some(500), false, false)
        );
        assert_eq!(
            parse_bnd_alt("A[1:100["),
            (Some("1".to_string()), Some(100), true, true)
        );
        // Single breakends have no mate info in brackets
        assert_eq!(parse_bnd_alt("G."), (None, None, false, false));
        assert_eq!(parse_bnd_alt(".A"), (None, None, false, false));
    }

    #[test]
    fn parse_sv_info_extracts_all_known_fields() {
        let info = "SVTYPE=DEL;END=200;SVLEN=-50;CN=1;IMPRECISE;NS=1";
        let parsed = parse_sv_info(info);
        assert_eq!(parsed.svtype.as_deref(), Some("DEL"));
        assert_eq!(parsed.end, Some(200));
        assert_eq!(parsed.svlen, Some(-50));
        assert_eq!(parsed.copy_number, Some(1));
    }

    #[test]
    fn parse_sv_info_dot_or_empty_returns_default() {
        // Bare "." (the "no info" sentinel) and "" both produce all-None.
        assert_eq!(parse_sv_info("."), ParsedSvInfo::default());
        assert_eq!(parse_sv_info(""), ParsedSvInfo::default());
    }

    #[test]
    fn parse_sv_info_ignores_unknown_keys_and_bare_flags() {
        let info = "DP=30;IMPRECISE;END=500;NS=42";
        let parsed = parse_sv_info(info);
        assert_eq!(parsed.end, Some(500));
        assert!(parsed.svtype.is_none());
        assert!(parsed.svlen.is_none());
        assert!(parsed.copy_number.is_none());
    }

    #[test]
    fn parse_sv_info_malformed_values_warn_and_drop() {
        // Malformed numeric values for known keys should produce None (not
        // a panic, not a returned error). Other fields still parse.
        let info = "END=not_a_number;SVLEN=10;CN=oops";
        let parsed = parse_sv_info(info);
        assert!(parsed.end.is_none());
        assert_eq!(parsed.svlen, Some(10));
        assert!(parsed.copy_number.is_none());
    }

    #[test]
    fn parse_sv_info_tolerates_empty_entries_and_trailing_semicolons() {
        // Real-world VCFs sometimes carry trailing semicolons or doubled
        // separators (`END=200;;SVLEN=-50`). The parser should treat each
        // empty entry as a no-op and still pick up the surrounding fields
        // rather than dropping the variant or panicking on `split_once`.
        let parsed = parse_sv_info("END=200;;SVLEN=-50;");
        assert_eq!(parsed.end, Some(200));
        assert_eq!(parsed.svlen, Some(-50));
    }

    #[test]
    fn parse_sv_info_last_value_wins_on_duplicate_keys() {
        // Spec-violating but real: some upstream tools emit the same key
        // twice. Lock in last-value-wins so a future refactor doesn't
        // silently flip to first-value-wins and break round-tripping.
        let parsed = parse_sv_info("END=100;END=999;CN=2");
        assert_eq!(parsed.end, Some(999));
        assert_eq!(parsed.copy_number, Some(2));
    }

    #[test]
    fn from_file_literal_complex_with_end_info_stays_literal() {
        // A literal-base Complex variant (multi-base REF *and* literal multi-
        // base ALT) is unusual but legal. parse_sv_info is gated on the ALT
        // being symbolic, so an `END=` field on a literal Complex record must
        // NOT promote it to Symbolic — the AlternateType stays Literal and the
        // raw INFO string is preserved verbatim for the writer. Locks in that
        // the two paths don't bleed.
        let v = Variant::from_file(
            42,
            ".",
            "PASS",
            "END=999;SVLEN=-3",
            "ACG",
            "TTT",
            60,
            vec!["GT".to_string()],
            vec!["0/1".to_string()],
        )
        .unwrap();
        assert_eq!(v.variant_type, VariantType::Complex);
        assert!(v.alternate.is_literal());
        assert!(v.alternate.as_symbolic().is_none());
        assert_eq!(
            v.alternate.as_literal().unwrap(),
            &[Nucleotide::T, Nucleotide::T, Nucleotide::T][..]
        );
        // INFO survives verbatim for the writer to round-trip; nothing has
        // re-interpreted END/SVLEN on a literal record.
        assert_eq!(v.info.as_deref(), Some("END=999;SVLEN=-3"));
    }

    #[test]
    fn sv_type_matches_svtype_strict_and_unknown_passthrough() {
        // Exact match (case-insensitive)
        assert!(sv_type_matches_svtype(SvType::Del, "DEL"));
        assert!(sv_type_matches_svtype(SvType::Del, "del"));
        assert!(sv_type_matches_svtype(SvType::Cnv, "CNV"));
        // Disagreement
        assert!(!sv_type_matches_svtype(SvType::Del, "DUP"));
        // Unknown tag is permissive (we can't contradict what we didn't recognize)
        assert!(sv_type_matches_svtype(SvType::Unknown, "WHATEVER"));
    }

    #[test]
    fn alternate_type_literal_as_literal_returns_slice() {
        let alt = AlternateType::Literal(vec![Nucleotide::A, Nucleotide::C]);
        assert_eq!(alt.as_literal(), Some(&[Nucleotide::A, Nucleotide::C][..]));
    }

    #[test]
    fn alternate_type_symbolic_as_literal_returns_none() {
        let alt = AlternateType::Symbolic(SvData::new("<DUP:TANDEM>", SvType::Dup));
        assert!(alt.as_literal().is_none());
        // The dispatcher round-trip: symbolic ALTs go through as_symbolic().
        assert!(alt.is_symbolic());
        assert_eq!(alt.as_symbolic().unwrap().raw_alt, "<DUP:TANDEM>");
        assert_eq!(alt.as_symbolic().unwrap().sv_type, SvType::Dup);
    }

    #[test]
    fn alternate_type_literal_helpers() {
        let alt = AlternateType::Literal(vec![Nucleotide::A]);
        assert!(alt.is_literal());
        assert!(!alt.is_symbolic());
        assert!(alt.as_symbolic().is_none());
    }

    #[test]
    fn sv_type_from_alt_string_recognizes_standard_tags() {
        assert_eq!(SvType::from_alt_string("<DEL>"), Some(SvType::Del));
        assert_eq!(SvType::from_alt_string("<DEL:ME>"), Some(SvType::Del));
        assert_eq!(SvType::from_alt_string("<DUP>"), Some(SvType::Dup));
        assert_eq!(SvType::from_alt_string("<DUP:TANDEM>"), Some(SvType::Dup));
        assert_eq!(SvType::from_alt_string("<CNV>"), Some(SvType::Cnv));
        assert_eq!(SvType::from_alt_string("<INS>"), Some(SvType::Ins));
        assert_eq!(SvType::from_alt_string("<INS:ME:ALU>"), Some(SvType::Ins));
        assert_eq!(SvType::from_alt_string("<INV>"), Some(SvType::Inv));
    }

    #[test]
    fn sv_type_from_alt_string_breakend() {
        // VCF 4.2 breakend forms — any of these should be tagged Bnd.
        assert_eq!(SvType::from_alt_string("G]17:198982]"), Some(SvType::Bnd));
        assert_eq!(SvType::from_alt_string("]13:123456]T"), Some(SvType::Bnd));
        assert_eq!(SvType::from_alt_string("[2:321682[A"), Some(SvType::Bnd));
    }

    #[test]
    fn sv_type_from_alt_string_unknown_tag() {
        // Angle-bracket form we don't have a dedicated variant for —
        // still classified as symbolic (Unknown), not literal.
        assert_eq!(
            SvType::from_alt_string("<NOVEL_THING>"),
            Some(SvType::Unknown)
        );
    }

    #[test]
    fn sv_type_from_alt_string_literal_returns_none() {
        // Literal-base ALTs shouldn't pass through this parser at all;
        // callers gate it on AlternateType, but defend the contract.
        assert_eq!(SvType::from_alt_string("A"), None);
        assert_eq!(SvType::from_alt_string("ACGT"), None);
        assert_eq!(SvType::from_alt_string(""), None);
    }

    #[test]
    fn sv_data_span_from_end_field() {
        // END is 1-based inclusive — for a deletion at POS=100, END=109
        // covers 10 bases (100..=109).
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(109);
        assert_eq!(sv.span(100), Some(10));
    }

    #[test]
    fn sv_data_span_falls_back_to_svlen() {
        // SVLEN is signed per VCF 4.2: negative for deletions, positive
        // for insertions/duplications. span() takes the absolute value.
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.svlen = Some(-50);
        assert_eq!(sv.span(100), Some(50));

        let mut sv = SvData::new("<DUP>", SvType::Dup);
        sv.svlen = Some(200);
        assert_eq!(sv.span(100), Some(200));
    }

    #[test]
    fn sv_data_span_prefers_end_over_svlen() {
        // If both are populated and disagree, END wins (it's the more
        // authoritative coordinate in VCF 4.2).
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(109);
        sv.svlen = Some(-999);
        assert_eq!(sv.span(100), Some(10));
    }

    #[test]
    fn sv_data_span_none_when_neither_set() {
        let sv = SvData::new("<DEL>", SvType::Del);
        assert!(sv.span(100).is_none());
    }

    #[test]
    fn sv_data_span_returns_none_when_end_before_location() {
        // INFO/END < POS is malformed (typo, swapped fields, etc.). The old
        // implementation used saturating_sub and silently produced Some(1),
        // which then gave <DUP>/<CNV>/<INV> records a bogus 1-base modulation
        // range. span() must now return None so downstream callers warn and
        // skip coverage modulation instead of acting on garbage.
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(50);
        assert!(sv.span(100).is_none());

        // Equality is allowed — END == POS means a degenerate 1-base span,
        // which is still semantically meaningful for DUP/INV (and gets
        // collapsed to empty for DEL by sv_modulation_range downstream).
        sv.end = Some(100);
        assert_eq!(sv.span(100), Some(1));
    }

    #[test]
    fn sv_data_span_invalid_end_does_not_fall_back_to_svlen() {
        // If END is present but malformed (< location), treat the whole record
        // as un-spannable rather than silently substituting SVLEN. The user
        // should fix their VCF; Variant::from_file already drops invalid END
        // at parse time with a warn, so this branch is defense-in-depth for
        // anything constructing SvData directly.
        let mut sv = SvData::new("<DEL>", SvType::Del);
        sv.end = Some(50);
        sv.svlen = Some(-200);
        assert!(sv.span(100).is_none());
    }

    #[test]
    fn from_file_invalid_end_is_dropped_and_falls_back_to_svlen() {
        // POS=100, END=50 (malformed) but SVLEN=-30 — Variant::from_file
        // clears the bad END, span() then falls back to |SVLEN|.
        let v = Variant::from_file(
            100,
            ".",
            "PASS",
            "SVTYPE=DEL;END=50;SVLEN=-30",
            "A",
            "<DEL>",
            60,
            vec!["GT".to_string()],
            vec!["1/1".to_string()],
        )
        .unwrap();
        let sv = v.alternate.as_symbolic().unwrap();
        assert!(sv.end.is_none(), "invalid END should be cleared");
        assert_eq!(sv.svlen, Some(-30));
        assert_eq!(sv.span(100), Some(30));
    }

    #[test]
    fn parse_alternate_non_iupac_ref_is_invalid_vcf() {
        // REF must be literal bases — a stray bracket here means a malformed VCF row.
        let err = parse_alternate("A<", "C").unwrap_err();
        match err {
            VariantError::InvalidVcf(msg) => assert!(msg.contains("REF")),
            other => panic!("expected InvalidVcf, got {other:?}"),
        }
    }

    #[test]
    fn allele_fraction_prefers_info_af() {
        assert_eq!(
            parse_allele_fraction("DP=100;AF=0.25", &[], &[]),
            Some(0.25)
        );
    }

    #[test]
    fn allele_fraction_info_af_takes_first_of_multiallelic() {
        assert_eq!(parse_allele_fraction("AF=0.3,0.5", &[], &[]), Some(0.3));
    }

    #[test]
    fn allele_fraction_falls_back_to_ad() {
        let format = vec!["GT".to_string(), "AD".to_string(), "DP".to_string()];
        let sample = vec!["0/1".to_string(), "18,12".to_string(), "30".to_string()];
        // No INFO/AF → 12 alt / 30 total = 0.4.
        let af = parse_allele_fraction("DP=30", &format, &sample).unwrap();
        assert!(
            (af - 0.4).abs() < 1e-9,
            "AD-derived AF should be 0.4, got {af}"
        );
    }

    #[test]
    fn allele_fraction_none_when_absent() {
        assert_eq!(parse_allele_fraction("DP=30", &[], &[]), None);
    }

    #[test]
    fn allele_fraction_sanitizes_out_of_range() {
        assert_eq!(parse_allele_fraction("AF=0", &[], &[]), None);
        assert_eq!(parse_allele_fraction("AF=-0.1", &[], &[]), None);
        assert_eq!(parse_allele_fraction("AF=1.5", &[], &[]), Some(1.0));
    }

    #[test]
    fn from_file_populates_allele_fraction_from_info() {
        let v = Variant::from_file(
            100,
            ".",
            "PASS",
            "AF=0.15",
            "A",
            "C",
            37,
            vec!["GT".to_string()],
            vec!["0/1".to_string()],
        )
        .unwrap();
        assert_eq!(v.allele_fraction, Some(0.15));
    }

    #[test]
    fn from_file_ad_derived_allele_fraction() {
        let v = Variant::from_file(
            100,
            ".",
            "PASS",
            "DP=40",
            "A",
            "C",
            37,
            vec!["GT".to_string(), "AD".to_string()],
            vec!["0/1".to_string(), "30,10".to_string()],
        )
        .unwrap();
        // 10 alt / 40 total = 0.25.
        assert_eq!(v.allele_fraction, Some(0.25));
    }

    #[test]
    fn from_file_no_af_yields_none() {
        let v = Variant::from_file(
            100,
            ".",
            "PASS",
            "DP=40",
            "A",
            "C",
            37,
            vec!["GT".to_string()],
            vec!["0/1".to_string()],
        )
        .unwrap();
        assert_eq!(v.allele_fraction, None);
    }
}
