pub mod errors;
pub mod utils;

use crate::gen_frag_length_model::{
    errors::GenFragLengthModelError, utils::config::RunConfiguration,
};
use std::path::PathBuf;

pub fn main(config_file: &PathBuf) -> Result<(), GenFragLengthModelError> {
    let config = RunConfiguration::from(config_file)?;
    utils::runner::runner(&config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_main_end_to_end() {
        // Smoke test: config YAML → RunConfiguration → runner → serialized model.
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("frags.bam");
        // A spread of lengths, not three values repeated: a discrete model over three
        // distinct fragment lengths is refused now, and rightly -- it has no shape.
        let tlens: Vec<usize> = (300..=360)
            .flat_map(|l| std::iter::repeat_n(l, 4))
            .collect();
        write_frag_bam(&bam_path, &tlens);
        let output = temp.path().join("model.json.gz");
        let yaml = format!(
            "input_file: {}\noutput_file: {}\noverwrite_output: true\nmin_reads: 2\n",
            bam_path.display(),
            output.display()
        );
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "{}", yaml).unwrap();
        main(&tmp.path().to_path_buf()).unwrap();
        assert!(
            output.exists(),
            "output model file should have been written"
        );
    }

    fn write_frag_bam(path: &PathBuf, tlens: &[usize]) {
        use noodles::bam;
        use noodles::sam::{
            self as sam,
            alignment::{
                RecordBuf,
                io::Write as _,
                record::{
                    Flags, MappingQuality,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
        };
        let header = sam::Header::builder()
            .add_reference_sequence(
                b"chr1".to_vec(),
                Map::<ReferenceSequence>::new(std::num::NonZero::<usize>::new(1_000_000).unwrap()),
            )
            .build();
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).unwrap();
        for &tlen in tlens {
            let cigar: Cigar = [Op::new(Kind::Match, 4)].into_iter().collect();
            let mut record = RecordBuf::default();
            *record.flags_mut() = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
            *record.cigar_mut() = cigar;
            *record.sequence_mut() = Sequence::from(b"ACGT".as_ref());
            *record.mapping_quality_mut() = Some(MappingQuality::try_from(30u8).unwrap());
            *record.reference_sequence_id_mut() = Some(0);
            *record.mate_reference_sequence_id_mut() = Some(0);
            *record.template_length_mut() = tlen as i32;
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
}
