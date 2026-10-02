//! Reject config keys a subcommand does not read (#496).
//!
//! Every subcommand loads its YAML into a map and looks up the keys it knows. Any other
//! key was silently ignored, so a typo (`tumor_mutation_model:` for `tumor_model:`) left
//! that setting at its default and the run still reported success. Each parser now
//! passes every key it was given through `check_keys` before using any of them.

/// Return an error message naming each key in `keys` that is not in `known`, or `Ok`
/// when there are none. `context` names the subcommand or section, for the message.
pub fn check_keys<'a>(
    keys: impl IntoIterator<Item = &'a str>,
    known: &[&str],
    context: &str,
) -> Result<(), String> {
    check_keys_with_hints(keys, known, &[], context)
}

/// `check_keys`, plus `hints`: (key, advice) pairs for keys known to come from elsewhere,
/// such as a NEAT config. A hinted key is still rejected; the advice replaces the
/// edit-distance suggestion in the message.
pub fn check_keys_with_hints<'a>(
    keys: impl IntoIterator<Item = &'a str>,
    known: &[&str],
    hints: &[(&str, &str)],
    context: &str,
) -> Result<(), String> {
    let mut unknown: Vec<&str> = keys.into_iter().filter(|k| !known.contains(k)).collect();
    if unknown.is_empty() {
        return Ok(());
    }
    unknown.sort_unstable();
    let named: Vec<String> = unknown
        .iter()
        .map(|k| {
            if let Some((_, advice)) = hints.iter().find(|(h, _)| h == k) {
                return format!("`{k}` ({advice})");
            }
            match closest(k, known) {
                Some(s) => format!("`{k}` (did you mean `{s}`?)"),
                None => format!("`{k}`"),
            }
        })
        .collect();
    let mut accepted: Vec<&str> = known.to_vec();
    accepted.sort_unstable();
    Err(format!(
        "unknown config key(s) in {context}: {}. An unread key would leave its setting at \
         the default without warning, so the run stops here. Accepted keys: {}",
        named.join(", "),
        accepted.join(", ")
    ))
}

/// The accepted key nearest to `key` by edit distance, if it is close enough to be a typo.
fn closest<'k>(key: &str, known: &[&'k str]) -> Option<&'k str> {
    let limit = (key.len() / 3).clamp(1, 3);
    known
        .iter()
        .map(|k| (levenshtein(key, k), *k))
        .filter(|(d, _)| *d <= limit)
        .min()
        .map(|(_, k)| k)
}

fn levenshtein(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, cb) in b.iter().enumerate() {
            let sub = prev[j] + usize::from(ca != *cb);
            cur[j + 1] = sub.min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    const KNOWN: &[&str] = &["tumor_model", "normal_model", "purity"];

    #[test]
    fn known_keys_pass() {
        assert_eq!(check_keys(["purity", "tumor_model"], KNOWN, "x"), Ok(()));
        assert_eq!(check_keys([], KNOWN, "x"), Ok(()));
    }

    #[test]
    fn an_unknown_key_is_named_with_its_context() {
        let err =
            check_keys(["purity", "totally_bogus_key"], KNOWN, "gen-cancer-reads").unwrap_err();
        assert!(err.contains("in gen-cancer-reads:"), "{err}");
        assert!(err.contains("`totally_bogus_key`"), "{err}");
        assert!(!err.contains("did you mean"), "no near match exists: {err}");
        assert!(
            !err.contains("`purity`"),
            "a known key is not reported: {err}"
        );
    }

    #[test]
    fn a_typo_gets_a_suggestion() {
        let err = check_keys(["tumour_model"], KNOWN, "x").unwrap_err();
        assert!(
            err.contains("`tumour_model` (did you mean `tumor_model`?)"),
            "{err}"
        );
    }

    #[test]
    fn every_unknown_key_is_reported_in_sorted_order() {
        let err = check_keys(["zzz", "aaa", "purity"], KNOWN, "x").unwrap_err();
        let (a, z) = (err.find("`aaa`").unwrap(), err.find("`zzz`").unwrap());
        assert!(a < z, "{err}");
    }

    #[test]
    fn a_hint_replaces_the_suggestion() {
        let hints = [("tumour_model", "use `tumor_model`, the American spelling")];
        let err = check_keys_with_hints(["tumour_model"], KNOWN, &hints, "x").unwrap_err();
        assert!(
            err.contains("`tumour_model` (use `tumor_model`, the American spelling)"),
            "{err}"
        );
        assert!(!err.contains("did you mean"), "{err}");
    }

    #[test]
    fn levenshtein_known_answers() {
        assert_eq!(levenshtein("", ""), 0);
        assert_eq!(levenshtein("abc", ""), 3);
        assert_eq!(levenshtein("kitten", "sitting"), 3);
        assert_eq!(levenshtein("tumour_model", "tumor_model"), 1);
    }
}
