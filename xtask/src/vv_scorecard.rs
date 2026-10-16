//! Deterministic public V&V scorecard artifact lane.
//!
//! `generate` renders the committed scorecard artifacts from the seeded
//! validation corpus and the built-in adversarial registry. No retained
//! ledgered run-result or executed adversarial-assessment store exists yet
//! (e05/e07 scope), so the honest current projection supplies empty run and
//! assessment sets and the committed artifact renders loud NO-DATA cells
//! instead of fabricated outcomes. `check` regenerates both artifacts and
//! requires byte identity with the tracked files; it composes into
//! `check-all` and therefore into the DSR quality gate.

use std::path::Path;

use fs_vvreg::adversarial::{AdversarialAssessment, adversarial_registry};
use fs_vvreg::corpus::corpus;
use fs_vvreg::scorecard::{ScorecardRunRecord, build_scorecard};

use super::Violation;

pub(crate) const CHECK: &str = "vv-scorecard";
const MARKDOWN_PATH: &str = "vv-scorecard.md";
const JSON_PATH: &str = "vv-scorecard.json";

/// NO RUN-RECORD PRODUCER EXISTS (bead frankensim-uqbkz): this is a WIRING
/// gap, not a data gap. Nothing in the workspace constructs a
/// `ScorecardRunRecord` from a real solve — the solve spine's QoI stage
/// (frankensim-s2l9v) is the intended producer, and until it executes and
/// an admission rule binds a retained run to a corpus row, every per-QoI
/// cell reads NO-DATA by construction. This named empty source exists so
/// the next reader cannot mistake the scorecard for "an instrument
/// awaiting data": it is an instrument whose input port has no producer.
/// Replace this constant with the authenticated ingestion path when the
/// spine can supply a real run; do not populate it any other way.
const NO_RUN_RECORD_PRODUCER_EXISTS: &[ScorecardRunRecord] = &[];

/// Same wiring gap for executed adversarial assessments (e07 scope):
/// `0/8 executed` in the rendered artifact is this constant, not a result.
const NO_ASSESSMENT_PRODUCER_EXISTS: &[AdversarialAssessment] = &[];

fn render() -> Result<(String, String), String> {
    let scorecard = build_scorecard(
        corpus(),
        adversarial_registry(),
        NO_RUN_RECORD_PRODUCER_EXISTS,
        NO_ASSESSMENT_PRODUCER_EXISTS,
    )
    .map_err(|error| format!("cannot build the V&V scorecard: {error}"))?;
    Ok((scorecard.render_markdown(), scorecard.render_json()))
}

pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let (markdown, json) = render()?;
    std::fs::write(root.join(MARKDOWN_PATH), markdown)
        .map_err(|error| format!("cannot write {MARKDOWN_PATH}: {error}"))?;
    std::fs::write(root.join(JSON_PATH), json)
        .map_err(|error| format!("cannot write {JSON_PATH}: {error}"))?;
    Ok(())
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let (markdown, json) = match render() {
        Ok(artifacts) => artifacts,
        Err(detail) => {
            return vec![Violation {
                check: CHECK,
                crate_name: "<repo>".to_string(),
                detail,
            }];
        }
    };
    [(MARKDOWN_PATH, markdown), (JSON_PATH, json)]
        .into_iter()
        .filter_map(
            |(path, expected)| match std::fs::read_to_string(root.join(path)) {
                Ok(actual) if actual == expected => None,
                Ok(_) => Some(Violation {
                    check: CHECK,
                    crate_name: path.to_string(),
                    detail: "tracked scorecard is stale; run cargo run -p xtask -- generate-vv-scorecard"
                        .to_string(),
                }),
                Err(error) => Some(Violation {
                    check: CHECK,
                    crate_name: path.to_string(),
                    detail: format!("cannot read retained scorecard artifact: {error}"),
                }),
            },
        )
        .collect()
}
