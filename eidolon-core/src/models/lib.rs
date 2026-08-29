//! This library contains some data and functions that are useful in the other models
use crate::file_tools::file_io::create_output_file;
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use serde::{Deserialize, Serialize};
use std;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

pub fn model_writer<T: Serialize>(model: T, filename: &PathBuf) -> std::io::Result<()> {
    // This will take any serializable model and write it to file.
    let data = serde_json::to_vec(&model).unwrap();
    // Overwrite check is the caller's responsibility (e.g., config builder)
    let fileout = create_output_file(filename, true)
        .unwrap_or_else(|_| panic!("Error creating output {:?}", filename));
    let writer = BufWriter::new(fileout);
    let mut encoder = GzEncoder::new(writer, Compression::default());
    encoder.write_all(&data)?;
    Ok(())
}

pub fn model_reader<T>(filename: &PathBuf) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein =
        File::open(filename).unwrap_or_else(|_| panic!("Error opening file {:?}", filename));
    let reader = GzDecoder::new(BufReader::new(filein));
    let value: T = serde_json::from_reader(reader)?;
    Ok(value)
}

#[allow(unused)]
pub fn model_unzipped_reader<T>(filename: &PathBuf) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein =
        File::open(filename).unwrap_or_else(|_| panic!("Error opening file {:?}", filename));
    let reader = BufReader::new(filein);
    let value: T = serde_json::from_reader(reader)?;
    Ok(value)
}
