//! This library contains some data and functions that are useful in the other models
use crate::file_tools::file_io::create_output_file;
use flate2::{Compression, read::GzDecoder, write::GzEncoder};
use log::info;
use serde::{Deserialize, Serialize};
use std;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

/// Key under which the provenance stamp is stored in a model file (#708).
///
/// Chosen to sort and read as metadata rather than model content, and to be a name no model
/// struct will ever define.
pub const PROVENANCE_KEY: &str = "_eidolon";

/// Bumped ONLY when the model format changes such that a reader at the previous version would
/// misread a file rather than merely miss a field.
///
/// Adding a field with `#[serde(default)]` does NOT require a bump: an older reader ignores it
/// and a newer reader supplies the default, which is the compatibility the README documents.
/// A bump is for changes that silently alter meaning -- a field whose units change, a
/// distribution whose parameterization changes, a value that stops meaning what it did.
pub const MODEL_FORMAT_VERSION: u32 = 1;

/// Provenance stamp embedded in every model file written by `model_writer`.
///
/// WHY (#708): model-file compatibility is part of eidolon's public API, but a model carried
/// nothing identifying the version that wrote it, so neither direction could be checked.
/// Backward compatibility worked by convention -- each added field annotated
/// `#[serde(default)]` -- while forward compatibility failed SILENTLY: model structs do not set
/// `deny_unknown_fields`, so an older binary reading a newer model discarded what it did not
/// recognize, fell back to its own defaults, and reported nothing.
///
/// This stamp makes the second case detectable. Note what it cannot do: a binary released
/// BEFORE the stamp existed will never check for it, so it only protects mismatches from this
/// version forward.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ModelProvenance {
    /// Format version of the file, compared against [`MODEL_FORMAT_VERSION`] on read.
    pub format_version: u32,
    /// eidolon version that wrote the file. Informational -- a run's records should be able to
    /// say which build produced the model it used.
    pub written_by: String,
}

impl ModelProvenance {
    /// The stamp this build writes.
    pub fn current() -> Self {
        Self {
            format_version: MODEL_FORMAT_VERSION,
            written_by: env!("CARGO_PKG_VERSION").to_string(),
        }
    }
}

pub fn model_writer<T: Serialize>(model: T, filename: &PathBuf) -> std::io::Result<()> {
    // This will take any serializable model and write it to file.
    //
    // The provenance stamp is injected as a sibling key of the model's own fields rather than
    // wrapping them, so the file shape is unchanged for every existing reader: a reader that
    // does not know about `_eidolon` ignores it as an unknown field, exactly as it always did.
    let data = serde_json::to_vec(&stamped(&model)).unwrap();
    // Overwrite check is the caller's responsibility (e.g., config builder)
    let fileout = create_output_file(filename, true)
        .unwrap_or_else(|_| panic!("Error creating output {:?}", filename));
    let writer = BufWriter::new(fileout);
    let mut encoder = GzEncoder::new(writer, Compression::default());
    encoder.write_all(&data)?;
    Ok(())
}

/// Attach the provenance stamp to a model's serialized form.
///
/// A model that does not serialize to a JSON object is returned untouched: there is nowhere to
/// put the key, and silently changing such a file's shape would be worse than not stamping it.
fn stamped<T: Serialize>(model: &T) -> serde_json::Value {
    let mut value = serde_json::to_value(model).unwrap();
    if let Some(map) = value.as_object_mut() {
        map.insert(
            PROVENANCE_KEY.to_string(),
            serde_json::to_value(ModelProvenance::current()).unwrap(),
        );
    }
    value
}

/// Read a model, checking its provenance stamp first.
///
/// A file with no stamp was written before stamps existed. That is the documented backward-
/// compatible case: it loads, and the absence is logged rather than treated as an error.
pub fn model_reader<T>(filename: &PathBuf) -> Result<T, std::io::Error>
where
    T: for<'de> Deserialize<'de>,
{
    // This will take any serializable model and write it to file.
    let filein =
        File::open(filename).unwrap_or_else(|_| panic!("Error opening file {:?}", filename));
    let reader = GzDecoder::new(BufReader::new(filein));
    let mut value: serde_json::Value = serde_json::from_reader(reader)?;
    check_provenance(&mut value, filename)?;
    let value: T = serde_json::from_value(value)?;
    Ok(value)
}

/// Validate and remove the provenance stamp, leaving the model's own fields behind.
///
/// A file whose `format_version` is AHEAD of this build is refused. That is the whole point of
/// the stamp: previously such a file loaded, lost every field this build does not know, and
/// simulated with local defaults without saying so.
fn check_provenance(
    value: &mut serde_json::Value,
    filename: &PathBuf,
) -> Result<Option<ModelProvenance>, std::io::Error> {
    let Some(map) = value.as_object_mut() else {
        return Ok(None);
    };
    let Some(raw) = map.remove(PROVENANCE_KEY) else {
        info!(
            "Model {:?} carries no provenance stamp, so it was written before eidolon recorded \
             one. Fields added since deserialize to their shipped defaults.",
            filename
        );
        return Ok(None);
    };
    let provenance: ModelProvenance = serde_json::from_value(raw)?;
    if provenance.format_version > MODEL_FORMAT_VERSION {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "model {:?} is in format version {}, but this build of eidolon ({}) understands \
                 up to version {}. It was written by eidolon {}. Reading it here would silently \
                 discard everything this build does not recognize and simulate with local \
                 defaults instead, so it is refused. Use eidolon {} or newer, or rebuild the \
                 model with this version.",
                filename,
                provenance.format_version,
                env!("CARGO_PKG_VERSION"),
                MODEL_FORMAT_VERSION,
                provenance.written_by,
                provenance.written_by,
            ),
        ));
    }
    info!(
        "Model {:?} was written by eidolon {} (format version {})",
        filename, provenance.written_by, provenance.format_version
    );
    Ok(Some(provenance))
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::io::Read;

    /// A real pre-stamp model file. These shipped defaults are loaded by every eidolon run and
    /// were written long before provenance existed, which makes them the honest backward-
    /// compatibility fixture rather than a synthetic one built to pass.
    static PRE_STAMP_MODEL: &[u8] =
        include_bytes!("model_data/default_fragment_length_model.json.gz");

    #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
    struct Toy {
        name: String,
        value: u32,
    }

    fn write_raw(bytes: &[u8]) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("model.json.gz");
        std::fs::write(&path, bytes).unwrap();
        (dir, path)
    }

    /// Gzip a JSON value the way `model_writer` does, WITHOUT adding a stamp, so a file can be
    /// constructed exactly as an older eidolon would have written it.
    fn write_unstamped(value: &serde_json::Value) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("model.json.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder
            .write_all(&serde_json::to_vec(value).unwrap())
            .unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();
        (dir, path)
    }

    /// THE backward-compatibility case the README promises: a model written before stamps
    /// existed still loads. Uses a real shipped model file, not a fixture.
    #[test]
    fn a_model_written_before_stamps_existed_still_loads() {
        let (_dir, path) = write_raw(PRE_STAMP_MODEL);
        let loaded: serde_json::Value = model_reader(&path).unwrap();

        // Loading is necessary but not sufficient -- the content must be the file's own, with
        // nothing dropped. Compare against a parse that bypasses the provenance path entirely.
        let mut direct = String::new();
        GzDecoder::new(PRE_STAMP_MODEL)
            .read_to_string(&mut direct)
            .unwrap();
        let expected: serde_json::Value = serde_json::from_str(&direct).unwrap();
        assert_eq!(loaded, expected, "a stampless model must load unchanged");
        assert!(
            expected.get(PROVENANCE_KEY).is_none(),
            "fixture precondition: this file must genuinely predate the stamp"
        );
    }

    /// The stamp round-trips, and says which build wrote the file.
    #[test]
    fn a_stamped_model_round_trips_and_reports_this_build() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("toy.json.gz");
        let toy = Toy {
            name: "t".to_string(),
            value: 7,
        };
        model_writer(&toy, &path).unwrap();

        // The model itself survives the stamp.
        let back: Toy = model_reader(&path).unwrap();
        assert_eq!(back, toy);

        // And the stamp is present, with this build's values.
        let raw: serde_json::Value = {
            let mut s = String::new();
            GzDecoder::new(std::fs::File::open(&path).unwrap())
                .read_to_string(&mut s)
                .unwrap();
            serde_json::from_str(&s).unwrap()
        };
        let stamp: ModelProvenance = serde_json::from_value(
            raw.get(PROVENANCE_KEY)
                .expect("stamp must be written")
                .clone(),
        )
        .unwrap();
        assert_eq!(stamp.format_version, MODEL_FORMAT_VERSION);
        assert_eq!(stamp.written_by, env!("CARGO_PKG_VERSION"));
    }

    /// The forward case, which is the whole reason the stamp exists: a model from a FUTURE
    /// format is refused rather than silently stripped of everything this build cannot read.
    #[test]
    fn a_model_from_a_future_format_version_is_refused() {
        let (_dir, path) = write_unstamped(&json!({
            "name": "t",
            "value": 7,
            PROVENANCE_KEY: {
                "format_version": MODEL_FORMAT_VERSION + 1,
                "written_by": "99.0.0",
            },
        }));
        let err = model_reader::<Toy>(&path).unwrap_err();
        let msg = err.to_string();
        // Assert CONTENT: a reader must learn what the file is, what this build handles, and
        // what wrote it.
        assert!(
            msg.contains(&(MODEL_FORMAT_VERSION + 1).to_string()),
            "must name the file's format version: {msg}"
        );
        assert!(
            msg.contains(&MODEL_FORMAT_VERSION.to_string()),
            "must name the version this build understands: {msg}"
        );
        assert!(
            msg.contains("99.0.0"),
            "must name the writing eidolon: {msg}"
        );
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    /// Must not fire: a model AT the current format version loads normally. An off-by-one in
    /// the comparison would reject every file this build itself writes.
    #[test]
    fn a_model_at_the_current_format_version_loads() {
        let (_dir, path) = write_unstamped(&json!({
            "name": "t",
            "value": 7,
            PROVENANCE_KEY: {
                "format_version": MODEL_FORMAT_VERSION,
                "written_by": "3.0.0",
            },
        }));
        let back: Toy = model_reader(&path).unwrap();
        assert_eq!(
            back,
            Toy {
                name: "t".to_string(),
                value: 7
            }
        );
    }

    /// The stamp must be additive: strip it and the file is what it would have been without it.
    /// This is what lets `model_parity` keep its baselines unchanged.
    #[test]
    fn the_stamp_does_not_alter_model_content() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("toy.json.gz");
        let toy = Toy {
            name: "t".to_string(),
            value: 7,
        };
        model_writer(&toy, &path).unwrap();

        let mut s = String::new();
        GzDecoder::new(std::fs::File::open(&path).unwrap())
            .read_to_string(&mut s)
            .unwrap();
        let mut written: serde_json::Value = serde_json::from_str(&s).unwrap();
        written.as_object_mut().unwrap().remove(PROVENANCE_KEY);

        assert_eq!(written, serde_json::to_value(&toy).unwrap());
    }
}
