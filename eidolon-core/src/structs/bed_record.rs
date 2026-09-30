use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum BedErrors {
    #[error("Error creating bed record: {0}")]
    BedRecordError(String),
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct BedRecord {
    contig: String,
    pub start: usize,
    pub end: usize,
    pub mut_rate: Option<f64>,
}

impl BedRecord {
    pub fn new_bed_record(contig: String, start: usize, end: usize) -> Result<Self, BedErrors> {
        if start >= end {
            return Err(BedErrors::BedRecordError(
                "Region with start pos >= end_pos".to_string(),
            ));
        }
        // skip string parsing if we know this is a regular target bed
        Ok(BedRecord {
            contig,
            start,
            end,
            mut_rate: None,
        })
    }

    pub fn new_mut_region_record(
        contig: String,
        start: usize,
        end: usize,
        other: &str,
    ) -> Result<Self, BedErrors> {
        if start >= end {
            return Err(BedErrors::BedRecordError(
                "Region with start pos >= end_pos".to_string(),
            ));
        }
        let mut_rate = Some(Self::parse_other_for_mut(other)?);
        Ok(BedRecord {
            contig,
            start,
            end,
            mut_rate,
        })
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    pub fn get_contig(&self) -> String {
        self.contig.to_owned()
    }

    pub fn contains(&self, contig: &str, position: usize) -> bool {
        if self.contig == contig {
            (position >= self.start) && (position < self.end)
        } else {
            false
        }
    }

    /// Whether a VCF record at 1-based `POS` falls inside this 0-based,
    /// half-open BED interval. `[s, e)` covers `POS` values `s + 1 ..= e`.
    ///
    /// Use this, not [`contains`](Self::contains), for any position read from a
    /// VCF's `POS` column (or a `Variant::location`, which stores it
    /// unconverted). `POS = 0` (the VCF telomere convention) is outside every
    /// interval.
    pub fn contains_vcf_pos(&self, contig: &str, pos: usize) -> bool {
        pos >= 1 && self.contains(contig, pos - 1)
    }

    pub fn overlaps(&self, contig: &str, start: usize, end: usize) -> bool {
        if self.contig != contig {
            return false;
        }
        // Intervals overlap if they are not strictly one after another.
        // self: [self.start, self.end), other: [start, end)
        self.start < end && start < self.end
    }

    pub fn parse_other_for_mut(other: &str) -> Result<f64, BedErrors> {
        if let Some(index) = other.to_lowercase().find("mut_rate=") {
            let start_index = index + "mut_rate=".len();
            let remainder = &other[start_index..];
            // The value might be followed by a space, tab, or semicolon if there are other fields
            let end_index = remainder
                .find(|c: char| c.is_whitespace() || c == ';' || c == ',' || c == '|')
                .unwrap_or(remainder.len());
            let value_str = &remainder[..end_index];
            value_str.parse::<f64>().map_err(|_| {
                BedErrors::BedRecordError(format!(
                    "Error parsing text after 'mut_rate': {}",
                    value_str
                ))
            })
        } else {
            Err(BedErrors::BedRecordError(format!(
                "Failed to locate 'mut_rate' text after column 3: {}",
                other
            )))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_other_finds_mut_rate() {
        let other = "mut_rate=0.003";
        assert_eq!(BedRecord::parse_other_for_mut(other).unwrap(), 0.003);

        let other = "name=gene1;mut_rate=0.005;score=100";
        assert_eq!(BedRecord::parse_other_for_mut(other).unwrap(), 0.005);

        let other = "some_other_field mut_rate=0.001\tmore_fields";
        assert_eq!(BedRecord::parse_other_for_mut(other).unwrap(), 0.001);
    }

    #[test]
    fn test_mut_rate_parser_fails() {
        let other = "no_rate_here";
        assert!(BedRecord::parse_other_for_mut(other).is_err());
    }

    #[test]
    fn test_not_a_number_error() {
        let other = "mut_rate=not_a_number";
        assert!(BedRecord::parse_other_for_mut(other).is_err());
    }

    #[test]
    fn test_new_mut_region_record() {
        let contig = "chr1".to_string();
        let start = 100;
        let end = 200;
        let other = "name=gene1;mut_rate=0.005;score=100";
        let result = BedRecord::new_mut_region_record(contig, start, end, other).unwrap();
        assert_eq!(result.start, 100);
        assert_eq!(result.end, 200);
        assert_eq!(result.mut_rate, Some(0.005));
    }

    #[test]
    fn test_new_mut_region_record_bad_rate() {
        let contig = "chr1".to_string();
        let start = 100;
        let end = 200;
        let other = "mut_rate=not_a_number";
        assert!(BedRecord::new_mut_region_record(contig, start, end, other).is_err());
    }

    #[test]
    fn test_overlaps() {
        let record = BedRecord::new_bed_record("chr1".to_string(), 100, 200).unwrap();
        assert!(record.overlaps("chr1", 150, 250)); // end inside
        assert!(record.overlaps("chr1", 50, 150)); // start inside
        assert!(record.overlaps("chr1", 150, 160)); // fully inside
        assert!(record.overlaps("chr1", 50, 250)); // fully contains record
        assert!(!record.overlaps("chr1", 200, 300)); // fully after (end is exclusive)
        assert!(!record.overlaps("chr1", 0, 100)); // fully before
        assert!(!record.overlaps("chr2", 150, 250)); // wrong contig
    }

    #[test]
    fn test_contains() {
        let record = BedRecord::new_bed_record("chr1".to_string(), 100, 200).unwrap();
        assert!(record.contains("chr1", 100)); // start is inclusive
        assert!(record.contains("chr1", 199)); // last valid position
        assert!(!record.contains("chr1", 200)); // end is exclusive
        assert!(!record.contains("chr1", 99)); // before start
        assert!(!record.contains("chr2", 150)); // wrong contig
    }

    /// #770 known answer: `chr1 100 200` covers 1-based POS 101..=200.
    #[test]
    fn contains_vcf_pos_converts_one_based_pos() {
        let record = BedRecord::new_bed_record("chr1".to_string(), 100, 200).unwrap();
        assert!(record.contains_vcf_pos("chr1", 101), "first base of region");
        assert!(record.contains_vcf_pos("chr1", 200), "last base of region");
        assert!(!record.contains_vcf_pos("chr1", 100), "base before region");
        assert!(!record.contains_vcf_pos("chr1", 201), "base after region");
        assert!(!record.contains_vcf_pos("chr2", 150), "wrong contig");
    }

    /// POS 0 is the VCF telomere convention; it must not underflow or match a
    /// region starting at 0.
    #[test]
    fn contains_vcf_pos_zero_is_outside() {
        let record = BedRecord::new_bed_record("chr1".to_string(), 0, 10).unwrap();
        assert!(!record.contains_vcf_pos("chr1", 0));
        assert!(record.contains_vcf_pos("chr1", 1));
        assert!(record.contains_vcf_pos("chr1", 10));
        assert!(!record.contains_vcf_pos("chr1", 11));
    }

    #[test]
    fn new_bedrecord() {
        let contig = "chr1".to_string();
        let start = 0;
        let end = 900;
        let result = BedRecord::new_bed_record(contig, start, end).unwrap();
        assert_eq!(result.start, 0);
        assert_eq!(result.end, 900)
    }

    #[test]
    fn test_new_bed_record_edge_cases() {
        // Zero length region should fail
        assert!(BedRecord::new_bed_record("chr1".to_string(), 100, 100).is_err());
        // Negative length (start > end) should fail
        assert!(BedRecord::new_bed_record("chr1".to_string(), 200, 100).is_err());

        // Large coordinates
        let rec = BedRecord::new_bed_record("chr1".to_string(), 0, usize::MAX).unwrap();
        assert_eq!(rec.len(), usize::MAX);
    }

    #[test]
    fn test_get_contig() {
        let rec = BedRecord::new_bed_record("chrX".to_string(), 10, 20).unwrap();
        assert_eq!(rec.get_contig(), "chrX");
    }

    #[test]
    fn test_parse_other_case_insensitivity() {
        let other = "MUT_RATE=0.01";
        assert_eq!(BedRecord::parse_other_for_mut(other).unwrap(), 0.01);
    }

    #[test]
    fn test_parse_other_delimiters() {
        assert_eq!(
            BedRecord::parse_other_for_mut("mut_rate=0.1,other=1").unwrap(),
            0.1
        );
        assert_eq!(
            BedRecord::parse_other_for_mut("mut_rate=0.2|other=1").unwrap(),
            0.2
        );
    }

    #[test]
    fn test_parse_other_missing_rate() {
        let result = BedRecord::parse_other_for_mut("no_rate_here");
        assert!(result.is_err());
    }

    #[test]
    fn test_overlaps_more() {
        let record = BedRecord::new_bed_record("chr1".to_string(), 100, 200).unwrap();
        // Fully contains the record
        assert!(record.overlaps("chr1", 50, 250));
        // Exactly matches the record
        assert!(record.overlaps("chr1", 100, 200));
        // Touches end (non-overlapping because end is exclusive)
        assert!(!record.overlaps("chr1", 200, 250));
        // Touches start (non-overlapping because start is inclusive and we check [50, 100) and [100, 200))
        assert!(!record.overlaps("chr1", 50, 100));
    }
}
