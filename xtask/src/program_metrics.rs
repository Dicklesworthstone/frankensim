//! Deterministic program metrics dashboard artifact lane.
//!
//! `generate` renders the committed dashboard from the seeded validation
//! corpus, the built-in adversarial registry, and the capability maturity
//! registry. `check` regenerates both artifacts and requires byte identity
//! with the tracked files; it composes into `check-all` and therefore into the
//! DSR quality gate.
//!
//! `record` is deliberately a SEPARATE verb. The trend basis is the retained
//! history file, so generation never mutates its own basis and `check` stays
//! idempotent: rendering the dashboard twice cannot move the trend cells.
//! Recording a generation is an explicit act with a caller-supplied label.

use std::collections::BTreeMap;
use std::path::Path;

use fs_govern::program_metrics::{
    HistoryGeneration, MetricCell, MetricHistory, MetricObservation, ProgramDashboard,
    ProgramSources, build_dashboard, frankensim_rows,
};
use fs_vvreg::ContentHash;
use fs_vvreg::adversarial::adversarial_registry;
use fs_vvreg::corpus::corpus;
use fs_vvreg::scorecard::build_scorecard;

use super::Violation;
use crate::depgraph::{JsonParser, JsonValue};
use crate::maturity;

pub(crate) const CHECK: &str = "program-metrics";
const MARKDOWN_PATH: &str = "program-metrics.md";
const JSON_PATH: &str = "program-metrics.json";
const HISTORY_PATH: &str = "program-metrics-history.jsonl";
const MAX_HISTORY_BYTES: u64 = 8 * 1024 * 1024;

fn obj(value: &JsonValue) -> Option<&BTreeMap<String, JsonValue>> {
    match value {
        JsonValue::Object(map) => Some(map),
        _ => None,
    }
}

fn text(value: &JsonValue) -> Option<&str> {
    match value {
        JsonValue::String(value) => Some(value.as_str()),
        _ => None,
    }
}

fn number(value: &JsonValue) -> Option<u64> {
    match value {
        JsonValue::Number(digits) => digits.parse::<u64>().ok(),
        _ => None,
    }
}

fn content_hash(hex: &str) -> Option<ContentHash> {
    if hex.len() != 64 {
        return None;
    }
    let mut bytes = [0_u8; 32];
    for (index, chunk) in hex.as_bytes().chunks_exact(2).enumerate() {
        let pair = std::str::from_utf8(chunk).ok()?;
        bytes[index] = u8::from_str_radix(pair, 16).ok()?;
    }
    Some(ContentHash(bytes))
}

/// One recorded metric value. An unrecognized shape is a refusal, never a
/// silent default: a misparsed history would fabricate a trend.
fn parse_cell(value: &JsonValue) -> Result<MetricCell, String> {
    let map = obj(value).ok_or_else(|| "history metric value is not an object".to_string())?;
    let status = map
        .get("status")
        .and_then(text)
        .ok_or_else(|| "history metric value has no status".to_string())?;
    match status {
        "no-data" => Ok(MetricCell::NoData {
            reason: map
                .get("reason")
                .and_then(text)
                .unwrap_or("recorded without a reason")
                .to_string(),
            unblocked_by: map.get("unblocked_by").and_then(text).map(str::to_string),
        }),
        "count" => {
            let value = map
                .get("value")
                .and_then(number)
                .ok_or_else(|| "recorded count has no value".to_string())?;
            Ok(MetricCell::Measured(MetricObservation::count(value)))
        }
        "ratio" => {
            let numerator = map
                .get("numerator")
                .and_then(number)
                .ok_or_else(|| "recorded ratio has no numerator".to_string())?;
            let denominator = map
                .get("denominator")
                .and_then(number)
                .and_then(std::num::NonZeroU64::new)
                .ok_or_else(|| {
                    "recorded ratio has no nonzero denominator; an empty population is NO-DATA, \
                     never a zero rate"
                        .to_string()
                })?;
            Ok(MetricCell::Measured(MetricObservation::ratio(
                numerator,
                denominator,
            )))
        }
        other => Err(format!("unrecognized recorded metric status `{other}`")),
    }
}

/// Read the retained generation history. A missing file is an empty history,
/// not an error: the first dashboard legitimately has no prior generation.
fn read_history(root: &Path) -> Result<MetricHistory, String> {
    let path = root.join(HISTORY_PATH);
    let source = match std::fs::metadata(&path) {
        Err(_) => return Ok(MetricHistory::empty()),
        Ok(metadata) if metadata.len() > MAX_HISTORY_BYTES => {
            return Err(format!(
                "{HISTORY_PATH} exceeds the admitted {MAX_HISTORY_BYTES}-byte bound"
            ));
        }
        Ok(_) => std::fs::read_to_string(&path)
            .map_err(|error| format!("cannot read {HISTORY_PATH}: {error}"))?,
    };

    let mut generations = Vec::new();
    for (index, line) in source.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let line_number = index + 1;
        let parsed = JsonParser::new(line)
            .finish()
            .map_err(|error| format!("{HISTORY_PATH}:{line_number} is not valid JSON: {error}"))?;
        let map = obj(&parsed)
            .ok_or_else(|| format!("{HISTORY_PATH}:{line_number} is not a JSON object"))?;
        let generation = map
            .get("generation")
            .and_then(number)
            .ok_or_else(|| format!("{HISTORY_PATH}:{line_number} has no generation number"))?;
        let label = map
            .get("label")
            .and_then(text)
            .ok_or_else(|| format!("{HISTORY_PATH}:{line_number} has no label"))?;
        let identity = map
            .get("source_identity")
            .and_then(text)
            .and_then(content_hash)
            .ok_or_else(|| {
                format!("{HISTORY_PATH}:{line_number} has no 32-byte hex source_identity")
            })?;
        let recorded = map
            .get("metrics")
            .and_then(obj)
            .ok_or_else(|| format!("{HISTORY_PATH}:{line_number} has no metrics object"))?;
        let mut metrics = BTreeMap::new();
        for (id, value) in recorded {
            let cell = parse_cell(value)
                .map_err(|error| format!("{HISTORY_PATH}:{line_number} metric `{id}`: {error}"))?;
            metrics.insert(id.clone(), cell);
        }
        generations.push(
            HistoryGeneration::try_new(generation, label, identity, metrics).map_err(|error| {
                format!("{HISTORY_PATH}:{line_number} is not admissible: {error}")
            })?,
        );
    }
    MetricHistory::try_new(generations)
        .map_err(|error| format!("{HISTORY_PATH} is not a valid trend basis: {error}"))
}

/// Project the live registries into the dashboard.
fn dashboard(root: &Path) -> Result<ProgramDashboard, String> {
    let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[])
        .map_err(|error| format!("cannot build the V&V scorecard: {error}"))?;
    let levels = maturity::capability_levels(root)?.current;
    let sources =
        ProgramSources::from_registries(corpus(), adversarial_registry(), &scorecard, &levels)
            .map_err(|error| format!("cannot read the program sources: {error}"))?;
    let rows = frankensim_rows(sources)
        .map_err(|error| format!("cannot project the program metrics: {error}"))?;
    let history = read_history(root)?;
    build_dashboard(&rows, &history)
        .map_err(|error| format!("cannot build the program metrics dashboard: {error}"))
}

fn render(root: &Path) -> Result<(String, String), String> {
    let dashboard = dashboard(root)?;
    Ok((dashboard.render_markdown(), dashboard.render_json()))
}

pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let (markdown, json) = render(root)?;
    std::fs::write(root.join(MARKDOWN_PATH), markdown)
        .map_err(|error| format!("cannot write {MARKDOWN_PATH}: {error}"))?;
    std::fs::write(root.join(JSON_PATH), json)
        .map_err(|error| format!("cannot write {JSON_PATH}: {error}"))?;
    Ok(())
}

/// Append the current values to the retained history as a future trend basis.
///
/// This never rewrites an existing generation: the history is append-only, so
/// a past dashboard state cannot be retconned to make a later one look better.
pub(crate) fn record(root: &Path, label: &str) -> Result<(), String> {
    let dashboard = dashboard(root)?;
    let history = read_history(root)?;
    let generation = history.next_generation();
    let line = dashboard.render_history_line(generation, label);
    let path = root.join(HISTORY_PATH);
    let mut existing = std::fs::read_to_string(&path).unwrap_or_default();
    if !existing.is_empty() && !existing.ends_with('\n') {
        existing.push('\n');
    }
    existing.push_str(&line);
    std::fs::write(&path, existing)
        .map_err(|error| format!("cannot write {HISTORY_PATH}: {error}"))?;
    eprintln!("recorded program metrics generation {generation} ({label})");
    Ok(())
}

pub(crate) fn check(root: &Path) -> Vec<Violation> {
    let (markdown, json) = match render(root) {
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
                    detail: "tracked program metrics dashboard is stale; run cargo run -p xtask \
                             -- generate-program-metrics"
                        .to_string(),
                }),
                Err(error) => Some(Violation {
                    check: CHECK,
                    crate_name: path.to_string(),
                    detail: format!("cannot read retained dashboard artifact: {error}"),
                }),
            },
        )
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("fsim-program-metrics-{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir is creatable");
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root")
            .to_path_buf();
        std::fs::copy(
            root.join("capability-maturity.json"),
            dir.join("capability-maturity.json"),
        )
        .expect("the real maturity registry is readable");
        dir
    }

    fn repo_root() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root")
            .to_path_buf()
    }

    /// The e2e lane: project the REAL registries and log every row with its
    /// sources, so a reviewer can see what the dashboard actually read.
    #[test]
    fn real_registry_e2e_logs_every_row_and_its_sources() {
        let dashboard = dashboard(&repo_root()).expect("the real dashboard builds");
        for entry in dashboard.rows() {
            let metric = entry.metric();
            println!(
                "family={} id={} value={} trend={} sources={:?}",
                metric.family().slug(),
                metric.id(),
                match metric.cell() {
                    MetricCell::NoData { reason, .. } => format!("NO-DATA ({reason})"),
                    MetricCell::Measured(observation) => observation.render(),
                },
                entry.trend().slug(),
                metric.sources()
            );
        }

        // The bead's DONE-WHEN floor: at least five live rows plus honest gaps.
        assert!(
            dashboard.measured_count() >= 5,
            "expected >=5 measured rows, got {}",
            dashboard.measured_count()
        );
        assert!(
            dashboard.no_data_count() > 0,
            "a dashboard with no NO-DATA rows would mean the program has no gaps"
        );

        // Every measured row cites a source; every NO-DATA row explains itself.
        for entry in dashboard.rows() {
            let metric = entry.metric();
            assert!(
                !metric.honesty().trim().is_empty(),
                "row {} has no honesty statement",
                metric.id()
            );
            match metric.cell() {
                MetricCell::Measured(_) => assert!(
                    !metric.sources().is_empty(),
                    "measured row {} cites no source",
                    metric.id()
                ),
                MetricCell::NoData { reason, .. } => assert!(
                    !reason.trim().is_empty(),
                    "NO-DATA row {} gives no reason",
                    metric.id()
                ),
            }
        }
    }

    /// A stale artifact is caught, and a regenerated one is accepted. This is
    /// the gate that keeps the committed dashboard honest.
    #[test]
    fn check_catches_a_stale_artifact() {
        let root = repo_root();
        assert!(
            check(&root).is_empty(),
            "the committed dashboard artifacts are stale; run generate-program-metrics"
        );

        let dir = scratch("stale");
        std::fs::write(dir.join(MARKDOWN_PATH), "stale markdown").expect("write");
        std::fs::write(dir.join(JSON_PATH), "stale json").expect("write");
        let violations = check(&dir);
        assert_eq!(violations.len(), 2, "both artifacts should be reported");
        for violation in &violations {
            assert_eq!(violation.check, CHECK);
            assert!(violation.detail.contains("generate-program-metrics"));
        }
    }

    /// Recording is append-only and drives the trend across two generations.
    #[test]
    fn recorded_generations_drive_the_trend() {
        let dir = scratch("record");
        // Generation 1: no prior basis at all.
        let first = dashboard(&dir).expect("dashboard builds");
        assert_eq!(first.prior_generation(), None);

        record(&dir, "first baseline").expect("record succeeds");
        let history = read_history(&dir).expect("history parses");
        assert_eq!(history.len(), 1);
        assert_eq!(history.next_generation(), 2);

        // Generation 2 over unchanged sources is PROVABLY unchanged.
        let second = dashboard(&dir).expect("dashboard builds");
        assert_eq!(second.prior_generation(), Some(1));
        assert!(second.sources_identical_to_prior());
        assert_eq!(second.source_identity(), first.source_identity());

        // Appending never rewrites the earlier generation.
        record(&dir, "second baseline").expect("record succeeds");
        let grown = read_history(&dir).expect("history parses");
        assert_eq!(grown.len(), 2);
        assert_eq!(grown.next_generation(), 3);
        let raw = std::fs::read_to_string(dir.join(HISTORY_PATH)).expect("history readable");
        assert_eq!(raw.lines().count(), 2);
        assert!(raw.contains("first baseline"), "generation 1 was rewritten");
    }

    /// A malformed history is refused rather than silently treated as empty:
    /// a dropped basis would render every trend as `no prior generation` and
    /// quietly erase a regression.
    #[test]
    fn malformed_history_is_refused() {
        let dir = scratch("malformed");
        let cases = [
            (
                "{\"generation\":1,\"label\":\"x\",\"metrics\":{}}",
                "source_identity",
            ),
            (
                "{\"generation\":1,\"source_identity\":\"00\",\"label\":\"x\",\"metrics\":{}}",
                "source_identity",
            ),
            (
                "{\"label\":\"x\",\"source_identity\":\"{hash}\",\"metrics\":{}}",
                "generation",
            ),
            (
                "{\"generation\":1,\"label\":\"x\",\"source_identity\":\"{hash}\"}",
                "metrics",
            ),
        ];
        let hash = "0".repeat(64);
        for (line, expected) in cases {
            std::fs::write(dir.join(HISTORY_PATH), line.replace("{hash}", &hash))
                .expect("write history");
            let error = read_history(&dir).expect_err("malformed history must be refused");
            assert!(
                error.contains(expected),
                "expected `{expected}` in refusal, got `{error}`"
            );
        }

        // An empty population recorded as a zero-denominator ratio is refused:
        // that is the NO-DATA-versus-zero invariant, enforced on the way IN.
        std::fs::write(
            dir.join(HISTORY_PATH),
            format!(
                "{{\"generation\":1,\"label\":\"x\",\"source_identity\":\"{hash}\",\"metrics\":\
                 {{\"m\":{{\"status\":\"ratio\",\"numerator\":0,\"denominator\":0}}}}}}"
            ),
        )
        .expect("write history");
        let error = read_history(&dir).expect_err("zero denominator must be refused");
        assert!(error.contains("nonzero denominator"), "{error}");

        // A missing history file is an empty basis, not an error.
        let fresh = scratch("missing-history");
        assert!(
            read_history(&fresh)
                .expect("missing history is empty")
                .is_empty()
        );
    }
}
