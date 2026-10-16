//! Closure-evidence lint (bead hx4p, the failure-compounding admission
//! rule): a TECHNICAL-FAILURE bead (type `bug`) may not close without
//! naming its regression evidence or an explicit non-repro disposition.
//! A fixed bug that leaves no permanent test surface is a failure the
//! Gauntlet never compounds on — exactly what `fs_bisect::compound`
//! (bead 6nb.9) exists to prevent; this lint is its enforcement arm at
//! the tracker boundary.
//!
//! The accepted close-reason vocabulary (documented convention, checked
//! case-insensitively): evidence words — regression, test, battery,
//! suite, conformance, lane, gate, golden, famil(y), oracle, proof,
//! pass, covered, verified — or disposition words — no-repro / not
//! reproducible, duplicate, superseded, invalid, wontfix, false
//! positive, already fixed, upstream (tracked in another repo). The
//! historical base (37 closed bugs at introduction) satisfies this with
//! zero violations, so enforcement is universal, not grandfathered.

use crate::Violation;
use std::path::Path;

/// Words that count as regression evidence or an explicit disposition.
const ACCEPTED: &[&str] = &[
    "regression",
    "test",
    "battery",
    "suite",
    "conformance",
    "lane",
    "gate",
    "golden",
    "famil",
    "oracle",
    "proof",
    "pass",
    "covered",
    "verified",
    "no-repro",
    "no repro",
    "not reproducib",
    "duplicate",
    "supersed",
    "invalid",
    "wontfix",
    "false positive",
    "already fixed",
    "upstream",
];

/// Extract a JSON string field from one JSONL line, honoring backslash
/// escapes (close reasons contain quotes and newlines). Returns the raw
/// escaped content — good enough for case-insensitive word search.
fn json_str_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let tag = format!("\"{key}\":\"");
    let start = line.find(&tag)? + tag.len();
    let bytes = line.as_bytes();
    let mut i = start;
    while i < bytes.len() {
        match bytes[i] {
            b'\\' => i += 2,
            b'"' => return Some(&line[start..i]),
            _ => i += 1,
        }
    }
    None
}

/// The closure-evidence lint over `.beads/issues.jsonl`.
pub fn check_closures(root: &Path) -> Vec<Violation> {
    let mut violations = Vec::new();
    let path = root.join(".beads/issues.jsonl");
    let Ok(text) = std::fs::read_to_string(&path) else {
        // No tracker in this tree (e.g. a test fixture): nothing to lint.
        return violations;
    };
    for line in text.lines() {
        // Cheap pre-filters before any string extraction.
        if !(line.contains("\"issue_type\":\"bug\"") && line.contains("\"status\":\"closed\"")) {
            continue;
        }
        let id = json_str_field(line, "id").unwrap_or("<unparsed>");
        let reason = json_str_field(line, "close_reason").unwrap_or("");
        let lower = reason.to_ascii_lowercase();
        if reason.trim().is_empty() || !ACCEPTED.iter().any(|w| lower.contains(w)) {
            violations.push(Violation {
                check: "closure-evidence",
                crate_name: id.to_string(),
                detail: format!(
                    "closed bug {id} names no regression evidence or explicit disposition in \
                     its close reason — cite the regression test/battery/golden/family it \
                     leaves behind, or say why none exists (no-repro, duplicate, superseded, \
                     upstream); bead hx4p / fs_bisect::compound is the workflow"
                ),
            });
        }
    }
    violations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escaped_field_extraction_survives_quotes() {
        let line =
            r#"{"id":"x-1","close_reason":"fixed \"the\" bug; battery green","status":"closed"}"#;
        assert_eq!(
            json_str_field(line, "close_reason"),
            Some(r#"fixed \"the\" bug; battery green"#)
        );
        assert_eq!(json_str_field(line, "id"), Some("x-1"));
        assert_eq!(json_str_field(line, "missing"), None);
    }

    #[test]
    fn closure_lint_end_to_end_on_fixture_tracker() {
        let base = std::env::temp_dir().join(format!("fsim-closure-test-{}", std::process::id()));
        std::fs::create_dir_all(base.join(".beads")).unwrap();
        let jsonl = concat!(
            // Good: closed bug with regression evidence.
            r#"{"id":"t-good","issue_type":"bug","status":"closed","close_reason":"root-caused; regression pinned in the conformance battery"}"#,
            "\n",
            // Good: explicit disposition.
            r#"{"id":"t-dup","issue_type":"bug","status":"closed","close_reason":"duplicate of t-good"}"#,
            "\n",
            // BAD: closed bug, no evidence language.
            r#"{"id":"t-bad","issue_type":"bug","status":"closed","close_reason":"done"}"#,
            "\n",
            // Exempt: closed feature.
            r#"{"id":"t-feat","issue_type":"feature","status":"closed","close_reason":"shipped"}"#,
            "\n",
            // Exempt: open bug.
            r#"{"id":"t-open","issue_type":"bug","status":"open","close_reason":""}"#,
            "\n",
        );
        std::fs::write(base.join(".beads/issues.jsonl"), jsonl).unwrap();
        let v = check_closures(&base);
        assert_eq!(v.len(), 1, "exactly the seeded violation: {v:?}");
        assert!(v[0].crate_name == "t-bad", "{v:?}");
        assert!(v[0].detail.contains("regression"), "teaching hint: {v:?}");
        let _ = std::fs::remove_dir_all(&base);
    }
}
