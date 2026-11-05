//! These tools give some general interactions at the base level with file_io and folder_tools and some more specific functionality
//! with reading bioinformatics filetypes. As we expand the capabilities, we'll add in more filetype tools, such as bam and vcf.

pub mod bam_reader;
pub mod bam_writer;
pub mod bed_reader;
pub mod block_gz;
pub mod fasta_stream;
pub mod fastq_shuffle;
pub mod fastq_tools;
pub mod file_io;
pub mod folder_tools;
pub mod vcf_tools;
