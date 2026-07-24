//! Abstraction-consolidation governance (bead
//! `frankensim-extreal-program-f85xj.16.8`).
//!
//! `consolidation-review.json` records which crates no supported workflow
//! exercises, and what was decided about each. A hand-maintained list of that
//! kind rots silently: crates gain and lose consumers constantly, so a record
//! written once describes the tree it was written against and nothing later.
//!
//! This checker re-derives the usage sweep from `crates/*/Cargo.toml` on every
//! run and refuses three drift classes:
//!
//! 1. **Undispositioned accretion** — a crate is exercised by no supported
//!    workflow and no review has decided anything about it.
//! 2. **Stale disposition** — a crate carries a disposition but is now reached
//!    by a workflow. A `FREEZE` in that state is the important one: the crate
//!    gained a consumer, so the parked label is now false.
//! 3. **A broken sweep** — the declared known-exercised control crate is itself
//!    unreached, which means the closure is not working and an empty candidate
//!    list would otherwise look like a clean inventory.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub const CHECK: &str = "consolidation-review";
pub const RECORD_FILE: &str = "consolidation-review.json";
const RECORD_SCHEMA: &str = "frankensim-consolidation-review-v1";
const PROTOCOL_FILE: &str = "docs/CONSOLIDATION_REVIEW.md";
const POLICY_BEAD: &str = "frankensim-extreal-program-f85xj.16.8";
const MAX_RECORD_BYTES: usize = 1024 * 1024;

const ROOT_FIELDS: [&str; 4] = ["policy_bead", "protocol", "reviews", "schema"];

/// Dispositions. `RETIRE` is recordable but is a proposal only: agents never
/// execute a removal, per the repository's no-deletion rule.
const DISPOSITIONS: [&str; 5] = [
    "CONSOLIDATE",
    "FREEZE",
    "KEEP",
    "REPAIR-OR-QUARANTINE",
    "RETIRE",
];

/// Dispositions that assert the crate is parked *because* nothing consumes it.
/// If such a crate becomes reachable, the record is making a false statement.
const PARKED: [&str; 3] = ["CONSOLIDATE", "FREEZE", "REPAIR-OR-QUARANTINE"];

const REQUIRED_PROTOCOL: [&str; 6] = [
    RECORD_FILE,
    "check-consolidation",
    "usage sweep",
    "green precondition",
    "L3-promotion falsifier",
    "never execute",
];

pub struct ConsolidationReport {
    pub violations: Vec<Violation>,
    pub decisions: Vec<PolicyNote>,
}

fn violation(entity: &str, detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: entity.to_string(),
        detail: detail.into(),
    }
}

fn note(entity: &str, verdict: &'static str, detail: impl Into<String>) -> PolicyNote {
    PolicyNote {
        check: CHECK,
        crate_name: entity.to_string(),
        verdict,
        detail: detail.into(),
    }
}

fn obj(value: &JsonValue) -> Option<&BTreeMap<String, JsonValue>> {
    match value {
        JsonValue::Object(map) => Some(map),
        _ => None,
    }
}

fn arr(value: &JsonValue) -> Option<&[JsonValue]> {
    match value {
        JsonValue::Array(items) => Some(items),
        _ => None,
    }
}

fn text(value: &JsonValue) -> Option<&str> {
    match value {
        JsonValue::String(value) => Some(value),
        _ => None,
    }
}

/// One crate manifest reduced to the only thing the sweep needs: which `fs-*`
/// crates it depends on, over every dependency section.
///
/// Deliberately section-agnostic. A crate used only as a dev-dependency oracle
/// by a workflow *is* exercised by that workflow, so excluding dev edges would
/// overstate the candidate set.
pub fn manifest_edges(manifest: &str) -> BTreeSet<String> {
    let mut edges = BTreeSet::new();
    let mut in_deps = false;
    for line in manifest.lines() {
        let trimmed = line.trim();
        if let Some(section) = trimmed.strip_prefix('[').and_then(|s| s.strip_suffix(']')) {
            // `dependencies`, `dev-dependencies`, `build-dependencies`, and any
            // `target."cfg(...)".dependencies` form.
            in_deps = section.ends_with("dependencies");
            continue;
        }
        if !in_deps {
            continue;
        }
        let Some((name, _)) = trimmed.split_once('=') else {
            continue;
        };
        let name = name.trim().trim_matches('"');
        if name.starts_with("fs-")
            && name
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
        {
            edges.insert(name.to_string());
        }
    }
    edges
}

/// Crates reachable from the workflow roots. A crate outside this set is
/// exercised by no supported workflow.
fn reachable(roots: &[String], graph: &BTreeMap<String, BTreeSet<String>>) -> BTreeSet<String> {
    let mut seen = BTreeSet::new();
    let mut stack: Vec<String> = roots.to_vec();
    while let Some(name) = stack.pop() {
        if !graph.contains_key(&name) || !seen.insert(name.clone()) {
            continue;
        }
        if let Some(edges) = graph.get(&name) {
            stack.extend(edges.iter().cloned());
        }
    }
    seen
}

pub fn check_sources(
    record: &str,
    graph: &BTreeMap<String, BTreeSet<String>>,
    protocol: Option<&str>,
) -> ConsolidationReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();

    let parsed = match JsonParser::with_string_limit(record, MAX_RECORD_BYTES).finish() {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(RECORD_FILE, format!("invalid JSON: {error}")));
            return ConsolidationReport {
                violations,
                decisions,
            };
        }
    };
    let Some(root) = obj(&parsed) else {
        violations.push(violation(
            RECORD_FILE,
            "document root must be a JSON object",
        ));
        return ConsolidationReport {
            violations,
            decisions,
        };
    };
    let found: BTreeSet<&str> = root.keys().map(String::as_str).collect();
    let want: BTreeSet<&str> = ROOT_FIELDS.iter().copied().collect();
    if found != want {
        let missing: Vec<&str> = want.difference(&found).copied().collect();
        let extra: Vec<&str> = found.difference(&want).copied().collect();
        violations.push(violation(
            RECORD_FILE,
            format!("field set must be exact; missing={missing:?} unexpected={extra:?}"),
        ));
        return ConsolidationReport {
            violations,
            decisions,
        };
    }
    if root.get("schema").and_then(text) != Some(RECORD_SCHEMA) {
        violations.push(violation(
            RECORD_FILE,
            format!("`schema` must be {RECORD_SCHEMA:?}"),
        ));
    }
    if root.get("policy_bead").and_then(text) != Some(POLICY_BEAD) {
        violations.push(violation(
            RECORD_FILE,
            format!("`policy_bead` must be {POLICY_BEAD:?}"),
        ));
    }
    if root.get("protocol").and_then(text) != Some(PROTOCOL_FILE) {
        violations.push(violation(
            RECORD_FILE,
            format!("`protocol` must be {PROTOCOL_FILE:?}"),
        ));
    }

    let reviews = root.get("reviews").and_then(arr).unwrap_or(&[]);
    let Some(latest) = reviews.last().and_then(obj) else {
        violations.push(violation(
            RECORD_FILE,
            "`reviews` must contain at least one review object; the record is the review history",
        ));
        return ConsolidationReport {
            violations,
            decisions,
        };
    };
    let review_id = latest
        .get("review_id")
        .and_then(text)
        .unwrap_or("<unnamed>");

    let roots: Vec<String> = latest
        .get("roots")
        .and_then(arr)
        .unwrap_or(&[])
        .iter()
        .filter_map(text)
        .map(str::to_string)
        .collect();
    if roots.is_empty() {
        violations.push(violation(
            review_id,
            "`roots` must name the supported-workflow crates; an empty root list makes every crate look unreached",
        ));
        return ConsolidationReport {
            violations,
            decisions,
        };
    }
    for root_crate in &roots {
        if !graph.contains_key(root_crate) {
            violations.push(violation(
                review_id,
                format!("workflow root {root_crate:?} is not a crate in the tree"),
            ));
        }
    }

    let reached = reachable(&roots, graph);
    let unreached: BTreeSet<String> = graph
        .keys()
        .filter(|c| !reached.contains(*c))
        .cloned()
        .collect();

    // A sweep that silently stops working would report an empty candidate set,
    // which reads as a clean inventory. The control makes that fail loudly.
    let control = latest
        .get("validation")
        .and_then(obj)
        .and_then(|v| v.get("known_exercised_control"))
        .and_then(text)
        .unwrap_or("");
    let control_crate = control.split_whitespace().next().unwrap_or("");
    if control_crate.is_empty() || !graph.contains_key(control_crate) {
        violations.push(violation(
            review_id,
            "`validation.known_exercised_control` must start with a crate name that exists",
        ));
    } else if !reached.contains(control_crate) {
        violations.push(violation(
            review_id,
            format!(
                "control crate {control_crate:?} is NOT reached by the sweep; the closure is broken \
                 and an empty candidate list cannot be trusted"
            ),
        ));
    }

    let mut dispositioned: BTreeMap<String, String> = BTreeMap::new();
    for (index, entry) in latest
        .get("candidates")
        .and_then(arr)
        .unwrap_or(&[])
        .iter()
        .enumerate()
    {
        let entity = format!("{review_id}#candidates[{index}]");
        let Some(entry) = obj(entry) else {
            violations.push(violation(&entity, "candidate entries must be objects"));
            continue;
        };
        let Some(name) = entry.get("crate").and_then(text) else {
            violations.push(violation(&entity, "candidate needs a `crate` name"));
            continue;
        };
        let disposition = entry.get("disposition").and_then(text).unwrap_or("");
        if !DISPOSITIONS.contains(&disposition) {
            violations.push(violation(
                name,
                format!(
                    "disposition {disposition:?} is not one of {}",
                    DISPOSITIONS.join(", ")
                ),
            ));
        }
        if entry
            .get("rationale")
            .and_then(text)
            .is_none_or(|r| r.trim().is_empty())
        {
            violations.push(violation(
                name,
                "every disposition needs a non-empty `rationale`; an unexplained decision cannot be re-reviewed",
            ));
        }
        if !graph.contains_key(name) {
            violations.push(violation(
                name,
                format!(
                    "candidate is not a crate in the tree; remove the stale row from {review_id}"
                ),
            ));
        }
        if dispositioned
            .insert(name.to_string(), disposition.to_string())
            .is_some()
        {
            violations.push(violation(name, "candidate appears twice in one review"));
        }
    }

    // Drift 1: unreached and undecided.
    for name in &unreached {
        if !dispositioned.contains_key(name) {
            violations.push(violation(
                name,
                format!(
                    "no supported workflow exercises this crate and {review_id} records no disposition; \
                     add one (KEEP with a named consumer or trust rationale, CONSOLIDATE, FREEZE, or a \
                     RETIRE proposal) to {RECORD_FILE}"
                ),
            ));
        }
    }

    // Drift 2: a parked disposition that the tree has outgrown.
    for (name, disposition) in &dispositioned {
        if PARKED.contains(&disposition.as_str()) && reached.contains(name) {
            violations.push(violation(
                name,
                format!(
                    "disposition is {disposition} but a supported workflow now reaches this crate; \
                     the parked label is false -- unfreeze it and record the change"
                ),
            ));
        }
    }

    if let Some(protocol) = protocol {
        for required in REQUIRED_PROTOCOL {
            if !protocol.contains(required) {
                violations.push(violation(
                    PROTOCOL_FILE,
                    format!("consolidation protocol is missing {required:?}"),
                ));
            }
        }
    } else {
        violations.push(violation(PROTOCOL_FILE, "document is unreadable"));
    }

    decisions.push(note(
        "repository",
        "inventory",
        format!(
            "{review_id}: {} crates, {} workflow roots, {} reached, {} exercised by no supported workflow, {} dispositioned",
            graph.len(),
            roots.len(),
            reached.len(),
            unreached.len(),
            dispositioned.len()
        ),
    ));
    for disposition in DISPOSITIONS {
        let count = dispositioned
            .values()
            .filter(|value| value.as_str() == disposition)
            .count();
        if count > 0 {
            decisions.push(note(
                disposition,
                "disposition",
                format!("{count} crate(s)"),
            ));
        }
    }

    ConsolidationReport {
        violations,
        decisions,
    }
}

fn read_graph(root: &Path) -> BTreeMap<String, BTreeSet<String>> {
    let mut graph = BTreeMap::new();
    let Ok(entries) = std::fs::read_dir(root.join("crates")) else {
        return graph;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let manifest = path.join("Cargo.toml");
        if !manifest.is_file() {
            continue;
        }
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if let Ok(text) = std::fs::read_to_string(&manifest) {
            graph.insert(name.to_string(), manifest_edges(&text));
        }
    }
    graph
}

pub fn check_consolidation(root: &Path) -> ConsolidationReport {
    let record = match std::fs::read_to_string(root.join(RECORD_FILE)) {
        Ok(record) => record,
        Err(error) => {
            return ConsolidationReport {
                violations: vec![violation(
                    RECORD_FILE,
                    format!("file is unreadable: {error}"),
                )],
                decisions: Vec::new(),
            };
        }
    };
    let protocol = std::fs::read_to_string(root.join(PROTOCOL_FILE)).ok();
    check_sources(&record, &read_graph(root), protocol.as_deref())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(pairs: &[(&str, &[&str])]) -> BTreeMap<String, BTreeSet<String>> {
        pairs
            .iter()
            .map(|(name, deps)| {
                (
                    name.to_string(),
                    deps.iter().map(|d| d.to_string()).collect(),
                )
            })
            .collect()
    }

    fn protocol() -> String {
        REQUIRED_PROTOCOL
            .iter()
            .map(|term| format!("{term}\n"))
            .collect()
    }

    fn record(candidates: &str) -> String {
        format!(
            r#"{{
              "schema": "{RECORD_SCHEMA}",
              "policy_bead": "{POLICY_BEAD}",
              "protocol": "{PROTOCOL_FILE}",
              "reviews": [{{
                "review_id": "TEST-A",
                "roots": ["fs-demo-e2e"],
                "validation": {{"known_exercised_control": "fs-used (1 dependent)"}},
                "candidates": [{candidates}]
              }}]
            }}"#
        )
    }

    const ORPHAN: &str =
        r#"{"crate": "fs-orphan", "disposition": "FREEZE", "rationale": "no consumer"}"#;

    /// The tree the bead's TESTS section describes: one exercised crate, one
    /// synthetic orphan.
    fn seeded() -> BTreeMap<String, BTreeSet<String>> {
        graph(&[
            ("fs-demo-e2e", &["fs-used"]),
            ("fs-used", &[]),
            ("fs-orphan", &[]),
        ])
    }

    #[test]
    fn a_known_exercised_crate_is_never_flagged() {
        let report = check_sources(&record(ORPHAN), &seeded(), Some(&protocol()));
        assert!(
            !report
                .violations
                .iter()
                .any(|item| item.crate_name == "fs-used"),
            "the control crate must not be flagged: {:?}",
            report.violations
        );
        assert!(
            report.violations.is_empty(),
            "seeded tree with a dispositioned orphan is clean: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_synthetic_orphan_without_a_disposition_is_flagged() {
        let report = check_sources(&record(""), &seeded(), Some(&protocol()));
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.crate_name == "fs-orphan"
                    && item.detail.contains("records no disposition")),
            "expected undispositioned-accretion violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_frozen_crate_that_gained_a_consumer_is_caught() {
        // fs-orphan is now reached through the workflow root.
        let grown = graph(&[
            ("fs-demo-e2e", &["fs-used", "fs-orphan"]),
            ("fs-used", &[]),
            ("fs-orphan", &[]),
        ]);
        let report = check_sources(&record(ORPHAN), &grown, Some(&protocol()));
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("the parked label is false")),
            "a FREEZE must not survive gaining a consumer: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_broken_sweep_cannot_masquerade_as_a_clean_inventory() {
        // The control crate is unreachable, so the closure itself is suspect.
        let broken = graph(&[("fs-demo-e2e", &[]), ("fs-used", &[]), ("fs-orphan", &[])]);
        let report = check_sources(&record(ORPHAN), &broken, Some(&protocol()));
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("the closure is broken")),
            "expected control-crate violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_disposition_without_a_rationale_is_refused() {
        let bare = r#"{"crate": "fs-orphan", "disposition": "KEEP", "rationale": "  "}"#;
        let report = check_sources(&record(bare), &seeded(), Some(&protocol()));
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("needs a non-empty `rationale`")),
            "expected rationale violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn dev_dependency_edges_count_as_exercise() {
        let manifest = "[package]\nname = \"x\"\n\n[dependencies]\nfs-a = { path = \"../fs-a\" }\n\n[dev-dependencies]\nfs-b = { path = \"../fs-b\" }\n";
        let edges = manifest_edges(manifest);
        assert!(
            edges.contains("fs-a") && edges.contains("fs-b"),
            "{edges:?}"
        );
    }

    #[test]
    fn non_dependency_sections_do_not_contribute_edges() {
        let manifest = "[package]\nname = \"fs-thing\"\n\n[lints.rust]\nfs-not-a-dep = \"warn\"\n";
        assert!(manifest_edges(manifest).is_empty());
    }

    #[test]
    fn the_live_consolidation_review_is_clean() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let report = check_consolidation(root);
        assert!(
            report.violations.is_empty(),
            "live consolidation-review.json must be clean: {:?}",
            report.violations
        );
        assert!(
            report
                .decisions
                .iter()
                .any(|note| note.verdict == "inventory")
        );
    }
}
