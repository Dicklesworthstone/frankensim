//! Enumeration ratchet (bead frankensim-rc-1109-root-3gc12.8.1).
//!
//! AGENTS.md operating rule 1: coupled capability lands as a Component, a
//! connection or a description, never as a new pairing entry point. This check
//! counts the marks of enumeration in the two crates where it happened and
//! compares them with a tracked file. A count above the file is a violation.
//! A count below the file is also a violation: the commit that removed code
//! lowers the file, so the mark can only move down.
//!
//! Counted by a text scan of non-test source, one per declaration:
//! - `fn` items whose name contains `_with_`, at any visibility, so that trait
//!   methods and a second method of the same name on another type count;
//! - `struct`, `enum` and `type` items whose name ends in `Candidate`;
//! - `struct`, `enum` and `type` items whose name ends in `Error`;
//! - trait definitions.
//!
//! A name that is exactly `Error` or `Candidate` is not counted (it is an
//! associated type). Scanning of a file stops at a `#[cfg(test)]` line that is
//! followed by a `mod`. Source lines are recorded and not gated.
//!
//! No claim: this is a line scan, not a parser. A declaration inside a block
//! comment or a string literal is counted, and one produced by a macro, or
//! with an attribute on the same line, is not.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use crate::Violation;

pub(crate) const CHECK: &str = "enumeration-ratchet";
const RATCHET_PATH: &str = "enumeration-ratchet.json";
const SCOPES: [(&str, &str); 2] = [
    ("fs-couple", "crates/fs-couple/src"),
    ("fs-material", "crates/fs-material/src"),
];
const GATED: [&str; 4] = ["with_fns", "candidates", "errors", "traits"];

/// The marks found in one crate.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct Marks {
    with_fns: usize,
    candidates: usize,
    errors: usize,
    traits: usize,
    lines: usize,
}

impl Marks {
    fn count(&self, key: &str) -> usize {
        match key {
            "with_fns" => self.with_fns,
            "candidates" => self.candidates,
            "errors" => self.errors,
            "traits" => self.traits,
            _ => self.lines,
        }
    }
}

/// Name of the item a line declares with `keyword`, at any visibility.
fn declared<'a>(line: &'a str, keyword: &str) -> Option<&'a str> {
    let mut rest = line.trim_start();
    if let Some(after) = rest.strip_prefix("pub") {
        let after = match after.strip_prefix('(') {
            Some(scoped) => &scoped[scoped.find(')')? + 1..],
            None => after,
        };
        // `public_thing` is not `pub`: keep the line as it was.
        if after.starts_with(' ') {
            rest = after.trim_start();
        }
    }
    for qualifier in [
        "default ",
        "const ",
        "async ",
        "unsafe ",
        "extern \"C\" ",
        "extern \"system\" ",
    ] {
        rest = rest.strip_prefix(qualifier).unwrap_or(rest);
    }
    let name = rest.strip_prefix(keyword)?.strip_prefix(' ')?;
    let end = name
        .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
        .unwrap_or(name.len());
    (end > 0).then(|| &name[..end])
}

/// Whether `name` is longer than `suffix` and ends in it.
fn names_a(name: &str, suffix: &str) -> bool {
    name.len() > suffix.len() && name.ends_with(suffix)
}

/// Scan one source text into `marks`.
pub(crate) fn scan(text: &str, marks: &mut Marks) {
    let mut lines = text.lines().peekable();
    let mut in_tests = false;
    while let Some(line) = lines.next() {
        marks.lines += 1;
        if line.trim() == "#[cfg(test)]"
            && lines
                .peek()
                .is_some_and(|next| declared(next, "mod").is_some())
        {
            in_tests = true;
        }
        if in_tests {
            continue;
        }
        if declared(line, "fn").is_some_and(|name| name.contains("_with_")) {
            marks.with_fns += 1;
        }
        for keyword in ["struct", "enum", "type"] {
            if let Some(name) = declared(line, keyword) {
                marks.candidates += usize::from(names_a(name, "Candidate"));
                marks.errors += usize::from(names_a(name, "Error"));
            }
        }
        marks.traits += usize::from(declared(line, "trait").is_some());
    }
}

/// Scan every non-test `.rs` file under `dir`.
fn scan_tree(dir: &Path) -> Result<Marks, String> {
    let mut marks = Marks::default();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(current) = stack.pop() {
        let entries = std::fs::read_dir(&current)
            .map_err(|error| format!("cannot read {}: {error}", current.display()))?;
        let mut paths: Vec<PathBuf> = entries.flatten().map(|entry| entry.path()).collect();
        paths.sort();
        for path in paths {
            let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
            if path.is_dir() {
                if name != "tests" {
                    stack.push(path);
                }
            } else if path.extension().is_some_and(|extension| extension == "rs")
                && name != "tests.rs"
                && !name.ends_with("_tests.rs")
            {
                let text = std::fs::read_to_string(&path)
                    .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
                scan(&text, &mut marks);
            }
        }
    }
    Ok(marks)
}

/// The integer recorded under `"<scope>.<key>"`.
fn recorded(text: &str, scope: &str, key: &str) -> Option<usize> {
    let needle = format!("\"{scope}.{key}\":");
    let rest = text[text.find(&needle)? + needle.len()..].trim_start();
    let end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

fn render(all: &[(&str, Marks)]) -> String {
    let mut out = String::from(
        "{\n  \"schema\": \"frankensim-enumeration-ratchet-v2\",\n  \"bead\": \"frankensim-rc-1109-root-3gc12.8.1\",\n",
    );
    for (scope, marks) in all {
        for key in GATED.iter().chain(&["lines"]) {
            let _ = writeln!(out, "  \"{scope}.{key}\": {},", marks.count(key));
        }
    }
    out.push_str(
        "  \"no_claim\": \"line-scan counts of declarations in non-test source: fn names containing _with_ at any visibility; struct, enum and type names ending in Candidate; struct, enum and type names ending in Error; trait definitions. Lines are recorded and not gated. A count says nothing about whether the code is correct or used.\"\n}\n",
    );
    out
}

fn scan_all(root: &Path) -> Result<Vec<(&'static str, Marks)>, String> {
    SCOPES
        .iter()
        .map(|(scope, dir)| Ok((*scope, scan_tree(&root.join(dir))?)))
        .collect()
}

/// Rewrite the ratchet file from the tree.
pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let all = scan_all(root)?;
    std::fs::write(root.join(RATCHET_PATH), render(&all))
        .map_err(|error| format!("cannot write {RATCHET_PATH}: {error}"))
}

/// Compare scanned marks with the text of the ratchet file.
pub(crate) fn compare(all: &[(&str, Marks)], tracked: &str) -> Vec<String> {
    let mut problems = Vec::new();
    for (scope, marks) in all {
        for key in GATED {
            let now = marks.count(key);
            match recorded(tracked, scope, key) {
                None => problems.push(format!("{RATCHET_PATH} has no \"{scope}.{key}\"")),
                Some(mark) if now > mark => problems.push(format!(
                    "{scope} {key} rose from {mark} to {now}: coupled work lands as a \
                     Component, a connection or a description (AGENTS.md operating rule 1); \
                     the added lines are in `git diff` of the crate's src"
                )),
                Some(mark) if now < mark => problems.push(format!(
                    "{scope} {key} fell from {mark} to {now}: lower the mark in the same commit \
                     with `cargo run -p xtask -- generate-enumeration-ratchet`"
                )),
                Some(_) => {}
            }
        }
    }
    problems
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let violation = |detail: String| Violation {
        check: CHECK,
        crate_name: RATCHET_PATH.to_string(),
        detail,
    };
    let all = match scan_all(root) {
        Ok(all) => all,
        Err(error) => return vec![violation(error)],
    };
    match std::fs::read_to_string(root.join(RATCHET_PATH)) {
        Ok(tracked) => compare(&all, &tracked).into_iter().map(violation).collect(),
        Err(error) => vec![violation(format!("cannot read {RATCHET_PATH}: {error}"))],
    }
}

#[cfg(test)]
mod tests {
    use super::{Marks, compare, render, scan};

    const SOURCE: &str = "\
pub fn advance_x_with_y() {}
    pub(crate) fn step_with_heat(a: f64) {}
fn private_with_nothing() {}
pub struct StepCandidate;
pub enum StepError {}
pub struct OtherError;
trait Hidden {}
pub trait Shown {
    fn drive_with_load(&self);
}
// pub fn note_with_text
impl Shown for StepCandidate {
    type Error = StepError;
    fn drive_with_load(&self) {}
}
pub type AliasError = StepError;
pub enum TrialCandidate {}
pub const unsafe fn raw_with_parts() {}
public_with_macro!(fn);
";

    const TESTS: &str = "\
#[cfg(test)]
mod tests {
    pub fn helper_with_fixture() {}
    pub struct FixtureError;
    trait Probe {}
}
";

    fn marks(text: &str) -> Marks {
        let mut marks = Marks::default();
        scan(text, &mut marks);
        marks
    }

    #[test]
    fn scan_counts_one_per_declaration_at_any_visibility() {
        let found = marks(SOURCE);
        // advance_x_with_y, step_with_heat, private_with_nothing, the trait
        // method, its implementation, raw_with_parts. Not the comment and
        // not the macro call.
        assert_eq!(found.with_fns, 6);
        // StepCandidate, TrialCandidate.
        assert_eq!(found.candidates, 2);
        // StepError, OtherError, AliasError. Not the associated `type Error`.
        assert_eq!(found.errors, 3);
        assert_eq!(found.traits, 2);
        assert_eq!(found.lines, 19);
    }

    #[test]
    fn a_trailing_test_module_is_not_counted() {
        let with_tests = format!("{SOURCE}{TESTS}");
        let found = marks(&with_tests);
        let plain = marks(SOURCE);
        assert_eq!(
            (found.with_fns, found.candidates, found.errors, found.traits),
            (plain.with_fns, plain.candidates, plain.errors, plain.traits)
        );
        assert_eq!(found.lines, plain.lines + 6);
        // The attribute alone, without a module after it, hides nothing.
        let gated = format!("#[cfg(test)]\npub fn only_with_tests() {{}}\n{SOURCE}");
        assert_eq!(marks(&gated).with_fns, plain.with_fns + 1);
    }

    #[test]
    fn a_new_pairing_fails_and_lowering_the_file_passes() {
        let tracked = render(&[("fs-couple", marks(SOURCE))]);
        assert!(compare(&[("fs-couple", marks(SOURCE))], &tracked).is_empty());

        // A second method with a name that already exists is a new pairing.
        let grown = format!("{SOURCE}impl Other {{\n    pub fn advance_x_with_y() {{}}\n}}\n");
        let problems = compare(&[("fs-couple", marks(&grown))], &tracked);
        assert_eq!(problems.len(), 1);
        assert!(
            problems[0].contains("with_fns rose from 6 to 7"),
            "{}",
            problems[0]
        );

        // Removing a pairing without lowering the file is refused too...
        let shrunk = SOURCE.replace("pub fn advance_x_with_y() {}\n", "");
        let problems = compare(&[("fs-couple", marks(&shrunk))], &tracked);
        assert!(
            problems[0].contains("with_fns fell from 6 to 5"),
            "{}",
            problems[0]
        );
        // ...and passes once the file is regenerated from the smaller tree.
        let lowered = render(&[("fs-couple", marks(&shrunk))]);
        assert!(compare(&[("fs-couple", marks(&shrunk))], &lowered).is_empty());
        // Lines alone never fail the gate.
        let longer = format!("{SOURCE}// one more line\n");
        assert!(compare(&[("fs-couple", marks(&longer))], &tracked).is_empty());
    }

    #[test]
    fn a_file_without_a_mark_is_refused() {
        let problems = compare(&[("fs-couple", marks(SOURCE))], "{}");
        assert_eq!(problems.len(), 4);
    }
}
