//! The 6 mm il Cannone card with source hair compliance,
//! `data/instruments/descriptions/violin/cannone-h6-bow-hair.material-strings`, is bound to the two committed
//! cards it is made from: the 6 mm card `cannone-dominant-g3-h6.material-strings` with exactly two rows added
//! after its 22nd line, `bow_ribbon` and `bow_compliance`, each copied verbatim from
//! `cannone-hair-compliance.material-strings` (Pitteroff and Woodhouse 1998, Part I).
//!
//! A data-form check, never a physical claim. An edit to one of the three files that is not carried to the
//! others fails here, so the open-string renders made with this card stay reproducible from the committed cards.

use std::fs;
use std::path::{Path, PathBuf};

const H6: &str =
    "crates/fs-couple/examples/material_strings/cannone-dominant-g3-h6.material-strings";
const HAIR: &str = "data/instruments/descriptions/violin/cannone-hair-compliance.material-strings";
const CARD: &str = "data/instruments/descriptions/violin/cannone-h6-bow-hair.material-strings";
/// The 6 mm card's line after which the two bow rows are added.
const INSERT_AFTER: usize = 22;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn read(file: &str) -> String {
    let path = root().join(file);
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn keyword(line: &str) -> Option<&str> {
    line.split_whitespace().next()
}

/// The one line of `card` whose first word is `key`.
fn row<'a>(card: &'a str, key: &str) -> &'a str {
    let rows: Vec<&str> = card.lines().filter(|l| keyword(l) == Some(key)).collect();
    assert_eq!(
        rows.len(),
        1,
        "{HAIR} must hold exactly one `{key}` row, found {}",
        rows.len()
    );
    rows[0]
}

#[test]
fn h6_bow_hair_card_is_the_h6_card_plus_the_two_committed_bow_rows() {
    let h6 = read(H6);
    let hair = read(HAIR);
    let card = read(CARD);
    let added = [row(&hair, "bow_ribbon"), row(&hair, "bow_compliance")];
    assert!(
        h6.lines().count() > INSERT_AFTER,
        "{H6} has no line {INSERT_AFTER}"
    );
    for key in ["bow_ribbon", "bow_compliance", "bow-residual-flexibility"] {
        assert!(
            !h6.lines().any(|l| keyword(l) == Some(key)),
            "{H6} already names `{key}`"
        );
    }
    let mut expected: Vec<&str> = Vec::new();
    for (i, line) in h6.lines().enumerate() {
        expected.push(line);
        if i + 1 == INSERT_AFTER {
            expected.extend(added);
        }
    }
    let got: Vec<&str> = card.lines().collect();
    assert_eq!(
        got, expected,
        "{CARD} must be {H6} with the `bow_ribbon` and `bow_compliance` rows of {HAIR} after line {INSERT_AFTER}"
    );
    assert!(card.ends_with('\n'), "{CARD} must end with a newline");
}
