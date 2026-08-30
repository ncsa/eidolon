pub struct ReadRecord {
    pub name: String,
    pub sequence: String,
    pub quality_scores: Vec<usize>,
    pub cigar_ops: Vec<char>,
    pub is_paired: bool,
    pub is_reverse: bool,
    pub contig: String,
    pub position: usize,
    pub mate_contig: String,
    pub mate_position: usize,
    pub template_length: i32,
    /// True when this read has no honest reference alignment — it lies wholly
    /// inside novel inserted sequence, which has no reference coordinate at all.
    /// Such a record is emitted with the SAM UNMAPPED flag and an empty CIGAR
    /// rather than being placed at the insertion's anchor with an all-insertion
    /// CIGAR, which consumes no reference and is not a valid alignment.
    ///
    /// By SAM convention an unmapped read with a mapped mate still carries the
    /// mate's RNAME/POS so the pair sorts together; that keeps the golden BAM
    /// usable as an answer key (#449) — the read is still located at the event
    /// it came from — without asserting an alignment that does not exist.
    pub is_unmapped: bool,
}
