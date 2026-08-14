//! Portfolio-candidate integrity gate (bead f85xj.16.11, slice 2).
//!
//! `candidate-portfolio-snapshot.json` is the machine-readable candidate
//! projection of PORTFOLIO_DISPOSITION_CONTRACT_V1. This check keeps the
//! committed candidate structurally sound against the LIVE bead graph:
//! every referenced issue and owner must still exist, the role vocabulary
//! must stay closed, no row may claim the primary role, the inertness
//! declaration must never be dropped, and the embedded producer validation
//! must record zero failures.
//!
//! Deliberate non-check, learned from the governance drill-rot incidents:
//! the snapshot's COUNTS are not compared against the live graph. They are
//! bound to their recorded beads data root and lag honestly; demanding
//! count freshness would make this gate red on every bead mutation and
//! train agents to regenerate mechanically (or ignore it). Referential
//! integrity cannot rot that way — a vanished issue is a real defect.

use std::collections::BTreeSet;
use std::path::Path;

use crate::depgraph::{JsonParser, JsonValue};

use super::Violation;

pub(crate) const CHECK: &str = "portfolio-candidate";
const CANDIDATE_PATH: &str = "candidate-portfolio-snapshot.json";
const BEADS_FILE: &str = ".beads/issues.jsonl";
const MAX_BYTES: usize = 1 << 20;
const SCHEMA: &str = "extreal.portfolio-disposition.v1";
const KIND: &str = "CandidatePortfolioAuthoritySnapshot";

const ROLE_VOCABULARY: [&str; 9] = [
    "primary_product_spine",
    "primary_enhancement",
    "prerequisite_consumed",
    "shared_infrastructure",
    "capped_research_or_second_vertical",
    "secondary_candidate",
    "deferred_by_gate",
    "superseded_or_absorbed",
    "historical_evidence",
];

fn violation(detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: CANDIDATE_PATH.to_string(),
        detail: detail.into(),
    }
}

fn live_issue_ids(beads: &str) -> BTreeSet<String> {
    let mut ids = BTreeSet::new();
    for line in beads.lines() {
        if line.trim().is_empty() || line.len() > MAX_BYTES {
            continue;
        }
        // Malformed rows are moonshot_policy's territory; this gate only
        // needs the id set of rows that DO parse.
        if let Ok(parsed) = JsonParser::with_string_limit(line, MAX_BYTES).finish() {
            if let JsonValue::Object(map) = parsed {
                if let Some(JsonValue::String(id)) = map.get("id") {
                    ids.insert(id.clone());
                }
            }
        }
    }
    ids
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    check_sources(
        &match std::fs::read_to_string(root.join(CANDIDATE_PATH)) {
            Ok(text) => text,
            Err(error) => {
                return vec![violation(format!("cannot read {CANDIDATE_PATH}: {error}"))];
            }
        },
        &match std::fs::read_to_string(root.join(BEADS_FILE)) {
            Ok(text) => text,
            Err(error) => return vec![violation(format!("cannot read {BEADS_FILE}: {error}"))],
        },
    )
}

fn check_sources(candidate: &str, beads: &str) -> Vec<Violation> {
    let mut v = Vec::new();
    let parsed = match JsonParser::with_string_limit(candidate, MAX_BYTES).finish() {
        Ok(value) => value,
        Err(error) => return vec![violation(format!("candidate is not valid JSON: {error}"))],
    };
    let JsonValue::Object(rootmap) = &parsed else {
        return vec![violation("candidate must be a JSON object")];
    };
    let text_field = |key: &str| -> Option<&str> {
        match rootmap.get(key) {
            Some(JsonValue::String(value)) => Some(value.as_str()),
            _ => None,
        }
    };

    match text_field("schema") {
        Some(SCHEMA) => {}
        other => v.push(violation(format!(
            "schema is {other:?}, expected {SCHEMA:?}"
        ))),
    }
    match text_field("kind") {
        Some(KIND) => {}
        other => v.push(violation(format!("kind is {other:?}, expected {KIND:?}"))),
    }
    // The inertness declaration is load-bearing: a candidate that stops
    // saying it cannot confer authority is on its way to being treated as
    // if it does.
    if !text_field("authority_statement").is_some_and(|s| s.contains("INERT CANDIDATE DATA")) {
        v.push(violation(
            "authority_statement must retain the INERT CANDIDATE DATA declaration",
        ));
    }

    let live = live_issue_ids(beads);
    let mut require_live = |label: &str, id: Option<&str>| match id {
        Some(id) if live.contains(id) => {}
        Some(id) => v.push(violation(format!(
            "{label} {id:?} does not resolve in the live bead graph"
        ))),
        None => v.push(violation(format!("{label} is missing"))),
    };
    require_live("primary_spine_id", text_field("primary_spine_id"));
    require_live(
        "primary_acceptance_owner",
        text_field("primary_acceptance_owner"),
    );
    require_live("activation_owner", text_field("activation_owner"));

    let Some(JsonValue::Array(dispositions)) = rootmap.get("dispositions") else {
        v.push(violation("dispositions must be a non-empty array"));
        return v;
    };
    if dispositions.is_empty() {
        v.push(violation("dispositions must be a non-empty array"));
    }
    let mut seen = BTreeSet::new();
    for (index, row) in dispositions.iter().enumerate() {
        let entity = format!("dispositions[{index}]");
        let JsonValue::Object(row) = row else {
            v.push(violation(format!("{entity} must be an object")));
            continue;
        };
        let row_text = |key: &str| -> Option<&str> {
            match row.get(key) {
                Some(JsonValue::String(value)) => Some(value.as_str()),
                _ => None,
            }
        };
        let Some(issue_id) = row_text("issue_id") else {
            v.push(violation(format!("{entity} has no issue_id")));
            continue;
        };
        if !seen.insert(issue_id.to_string()) {
            v.push(violation(format!("duplicate disposition for {issue_id:?}")));
        }
        for key in ["issue_id", "integration_owner", "acceptance_owner"] {
            match row_text(key) {
                Some(id) if live.contains(id) => {}
                Some(id) => v.push(violation(format!(
                    "{entity}.{key} {id:?} does not resolve in the live bead graph"
                ))),
                None => v.push(violation(format!("{entity} has no {key}"))),
            }
        }
        match row_text("spine_role") {
            Some("primary_product_spine") => v.push(violation(format!(
                "{entity} claims primary_product_spine; primary identity lives only in primary_spine_id"
            ))),
            Some(role) if ROLE_VOCABULARY.contains(&role) => {}
            Some(role) => v.push(violation(format!(
                "{entity} role {role:?} is outside the closed vocabulary"
            ))),
            None => v.push(violation(format!("{entity} has no spine_role"))),
        }
    }

    // The embedded producer validation must exist and record zero failures;
    // a candidate whose own producer saw failures should never have been
    // emitted, let alone committed.
    match rootmap.get("producer_validation") {
        Some(JsonValue::Object(pv)) => match pv.get("checks_failed") {
            Some(JsonValue::Number(n)) if n == "0" => {}
            other => v.push(violation(format!(
                "producer_validation.checks_failed must be 0, found {other:?}"
            ))),
        },
        _ => v.push(violation("producer_validation section is missing")),
    }
    // The live-graph binding must carry its identity roots (counts are NOT
    // freshness-checked — see the module doc for why).
    match rootmap.get("live_graph_binding") {
        Some(JsonValue::Object(binding)) => {
            for key in ["beads_data_root_sha256", "head_sha", "captured_at"] {
                match binding.get(key) {
                    Some(JsonValue::String(value)) if !value.is_empty() => {}
                    _ => v.push(violation(format!(
                        "live_graph_binding.{key} must be a non-empty string"
                    ))),
                }
            }
        }
        _ => v.push(violation("live_graph_binding section is missing")),
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn live_sources() -> (String, String) {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root")
            .to_path_buf();
        (
            std::fs::read_to_string(root.join(CANDIDATE_PATH)).expect("candidate"),
            std::fs::read_to_string(root.join(BEADS_FILE)).expect("beads"),
        )
    }

    #[test]
    fn the_live_candidate_is_clean() {
        let (candidate, beads) = live_sources();
        let violations = check_sources(&candidate, &beads);
        assert!(violations.is_empty(), "{violations:?}");
    }

    #[test]
    fn seeded_faults_fail_closed() {
        // Baseline-independent mutations with must-apply assertions — the
        // drill-rot lesson from the moonshot/program-metrics gates, applied
        // from birth here.
        let (candidate, beads) = live_sources();

        let dropped_inert = candidate.replacen("INERT CANDIDATE DATA", "advisory data", 1);
        assert_ne!(dropped_inert, candidate, "mutation must apply");
        assert!(
            check_sources(&dropped_inert, &beads)
                .iter()
                .any(|item| item.detail.contains("INERT")),
            "dropping the inertness declaration must fail"
        );

        let phantom = candidate.replacen(
            "frankensim-ext-flagship-geneva-b4hj",
            "frankensim-ext-flagship-phantom-zzzz",
            1,
        );
        assert_ne!(phantom, candidate, "mutation must apply");
        assert!(
            check_sources(&phantom, &beads)
                .iter()
                .any(|item| item.detail.contains("does not resolve")),
            "a phantom disposition target must fail"
        );

        let primary_grab = candidate.replacen(
            "\"spine_role\": \"primary_enhancement\"",
            "\"spine_role\": \"primary_product_spine\"",
            1,
        );
        assert_ne!(primary_grab, candidate, "mutation must apply");
        assert!(
            check_sources(&primary_grab, &beads)
                .iter()
                .any(|item| item.detail.contains("primary identity lives only")),
            "a row claiming the primary role must fail"
        );

        let unvalidated = candidate.replacen("\"checks_failed\": 0", "\"checks_failed\": 3", 1);
        assert_ne!(unvalidated, candidate, "mutation must apply");
        assert!(
            check_sources(&unvalidated, &beads)
                .iter()
                .any(|item| item.detail.contains("checks_failed")),
            "recorded producer failures must fail the gate"
        );
    }
}
