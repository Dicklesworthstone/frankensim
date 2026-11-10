//! Enumeration ratchet (bead frankensim-rc-1109-root-3gc12.8.1).
//!
//! AGENTS.md operating rule 1: coupled capability lands as a Component, a
//! connection or a description, never as a new pairing entry point. This check
//! counts the marks of enumeration in the two crates where it happened and
//! compares them with a tracked file. A count above the file is a violation.
//! A count below the file is also a violation: the commit that removed code
//! lowers the file, so the mark can only move down.
//!
//! Counted by a text scan of non-test source: distinct `pub fn` names
//! containing `_with_`, distinct `pub struct` names ending in `Candidate`,
//! distinct `pub struct` or `pub enum` names ending in `Error`, and trait
//! definitions. Source lines are recorded for the record and are not gated.

use std::collections::BTreeSet;
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
    with_fns: BTreeSet<String>,
    candidates: BTreeSet<String>,
    errors: BTreeSet<String>,
    traits: usize,
    lines: usize,
}

impl Marks {
    fn count(&self, key: &str) -> usize {
        match key {
            "with_fns" => self.with_fns.len(),
            "candidates" => self.candidates.len(),
            "errors" => self.errors.len(),
            "traits" => self.traits,
            _ => self.lines,
        }
    }
}

/// Name of the item a line declares with `keyword`, with or without `pub`.
fn declared<'a>(line: &'a str, keyword: &str, public_only: bool) -> Option<&'a str> {
    let mut rest = line.trim_start();
    if let Some(after) = rest.strip_prefix("pub") {
        rest = match after.strip_prefix('(') {
            Some(scoped) => &scoped[scoped.find(')')? + 1..],
            None => after,
        };
        if !rest.starts_with(' ') {
            return None;
        }
    } else if public_only {
        return None;
    }
    rest = rest.trim_start();
    for qualifier in ["const ", "async ", "unsafe "] {
        rest = rest.strip_prefix(qualifier).unwrap_or(rest);
    }
    let name = rest.strip_prefix(keyword)?.strip_prefix(' ')?;
    let end = name
        .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
        .unwrap_or(name.len());
    (end > 0).then(|| &name[..end])
}

/// Scan one source text into `marks`.
pub(crate) fn scan(text: &str, marks: &mut Marks) {
    for line in text.lines() {
        marks.lines += 1;
        if let Some(name) = declared(line, "fn", true)
            && name.contains("_with_")
        {
            marks.with_fns.insert(name.to_owned());
        }
        for keyword in ["struct", "enum"] {
            if let Some(name) = declared(line, keyword, true) {
                if keyword == "struct" && name.ends_with("Candidate") {
                    marks.candidates.insert(name.to_owned());
                }
                if name.ends_with("Error") {
                    marks.errors.insert(name.to_owned());
                }
            }
        }
        if declared(line, "trait", false).is_some() {
            marks.traits += 1;
        }
    }
}

/// Scan every non-test `.rs` file under `dir`, in sorted order.
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
        "{\n  \"schema\": \"frankensim-enumeration-ratchet-v1\",\n  \"bead\": \"frankensim-rc-1109-root-3gc12.8.1\",\n",
    );
    for (scope, marks) in all {
        for key in GATED.iter().chain(&["lines"]) {
            let _ = writeln!(out, "  \"{scope}.{key}\": {},", marks.count(key));
        }
    }
    out.push_str(
        "  \"no_claim\": \"text-scan counts of non-test source: pub fn names containing _with_, pub struct names ending in Candidate, pub struct or enum names ending in Error, trait definitions; lines are recorded and not gated; a count says nothing about whether the code is correct or used\"\n}\n",
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

    const SOURCE: &str = "pub fn advance_x_with_y() {}\n    pub(crate) fn step_with_heat(a: f64) {}\n\
        fn private_with_nothing() {}\npub struct StepCandidate;\npub enum StepError {}\n\
        pub struct OtherError;\ntrait Hidden {}\npub trait Shown {}\n// pub fn note_with_text\n";

    fn marks(text: &str) -> Marks {
        let mut marks = Marks::default();
        scan(text, &mut marks);
        marks
    }

    #[test]
    fn scan_counts_public_marks_only() {
        let found = marks(SOURCE);
        let names: Vec<&str> = found.with_fns.iter().map(String::as_str).collect();
        assert_eq!(names, ["advance_x_with_y", "step_with_heat"]);
        assert_eq!(
            (found.candidates.len(), found.errors.len(), found.traits),
            (1, 2, 2)
        );
        assert_eq!(found.lines, 9);
    }

    #[test]
    fn a_new_pairing_fails_and_lowering_the_file_passes() {
        let tracked = render(&[("fs-couple", marks(SOURCE))]);
        assert!(compare(&[("fs-couple", marks(SOURCE))], &tracked).is_empty());

        let grown = format!("{SOURCE}pub fn advance_x_with_y_and_z() {{}}\n");
        let problems = compare(&[("fs-couple", marks(&grown))], &tracked);
        assert_eq!(problems.len(), 1);
        assert!(
            problems[0].contains("with_fns rose from 2 to 3"),
            "{}",
            problems[0]
        );

        // Removing a pairing without lowering the file is refused too...
        let shrunk = SOURCE.replace("pub fn advance_x_with_y() {}\n", "");
        let problems = compare(&[("fs-couple", marks(&shrunk))], &tracked);
        assert!(
            problems[0].contains("with_fns fell from 2 to 1"),
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
