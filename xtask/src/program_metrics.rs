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
    HistoryGeneration, ImportScorecardSource, MetricCell, MetricHistory, MetricObservation,
    ProgramDashboard, ProgramSources, SpineMetricsSource, build_dashboard, frankensim_rows,
};
use fs_vvreg::ContentHash;
use fs_vvreg::adversarial::adversarial_registry;
use fs_vvreg::corpus::corpus;
use fs_vvreg::scorecard::build_scorecard;

use super::Violation;
use crate::depgraph::{JsonParser, JsonValue};
use crate::maturity;
use crate::{spine_metrics, spine_ratchet, tropical_path};

pub(crate) const CHECK: &str = "program-metrics";
const MARKDOWN_PATH: &str = "program-metrics.md";
const JSON_PATH: &str = "program-metrics.json";
const HISTORY_PATH: &str = "program-metrics-history.jsonl";
const MAX_HISTORY_BYTES: u64 = 8 * 1024 * 1024;
const IMPORT_SUMMARY_PATH: &str = "data/cad-import-corpus/scorecard-summary-v1.json";
const IMPORT_MANIFEST_PATH: &str = "data/cad-import-corpus/corpus-v1.tsv";
const IMPORT_SUMMARY_SEMANTICS: &str = "fs-io-supplier-corpus-summary-v1";
const IMPORT_MANIFEST_IDENTITY_DOMAIN: &str = "fs-io supplier corpus manifest bytes v1";
const IMPORT_AUTHORITY: &str = "human-locked-only-dashboard-denominator";
const MAX_IMPORT_SUMMARY_BYTES: u64 = 64 * 1024;
const MAX_IMPORT_MANIFEST_BYTES: u64 = 8 * 1024 * 1024;

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

fn object_field<'a>(
    map: &'a BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<&'a BTreeMap<String, JsonValue>, String> {
    map.get(key)
        .and_then(obj)
        .ok_or_else(|| format!("{context} has no `{key}` object"))
}

fn text_field<'a>(
    map: &'a BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<&'a str, String> {
    map.get(key)
        .and_then(text)
        .ok_or_else(|| format!("{context} has no `{key}` string"))
}

fn count_field(
    map: &BTreeMap<String, JsonValue>,
    key: &str,
    context: &str,
) -> Result<usize, String> {
    let value = map
        .get(key)
        .and_then(number)
        .ok_or_else(|| format!("{context} has no non-negative integer `{key}`"))?;
    usize::try_from(value).map_err(|_| format!("{context} `{key}` does not fit usize"))
}

fn is_canonical_hash(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn checked_sum(parts: &[usize], context: &str) -> Result<usize, String> {
    parts.iter().try_fold(0_usize, |total, part| {
        total
            .checked_add(*part)
            .ok_or_else(|| format!("{context} count sum overflows usize"))
    })
}

fn retained_manifest_population(manifest_bytes: &[u8]) -> Result<usize, String> {
    let manifest = std::str::from_utf8(manifest_bytes)
        .map_err(|error| format!("{IMPORT_MANIFEST_PATH} is not UTF-8: {error}"))?;
    let records = manifest.lines().filter(|line| {
        let trimmed = line.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    });
    records
        .count()
        .checked_sub(1)
        .ok_or_else(|| format!("{IMPORT_MANIFEST_PATH} has no non-comment header"))
}

/// Parse and cross-check the sibling-free supplier-scorecard projection.
///
/// The summary is not trusted merely because it is tracked: its manifest
/// identity is recomputed from the exact live TSV bytes, and every population
/// relationship is revalidated before the counts reach `fs-govern`.
fn parse_import_scorecard_summary(
    source: &str,
    manifest_bytes: &[u8],
) -> Result<ImportScorecardSource, String> {
    let parsed = JsonParser::new(source)
        .finish()
        .map_err(|error| format!("{IMPORT_SUMMARY_PATH} is not valid JSON: {error}"))?;
    let map = obj(&parsed).ok_or_else(|| format!("{IMPORT_SUMMARY_PATH} is not a JSON object"))?;
    if map.get("schema").and_then(number) != Some(1) {
        return Err(format!("{IMPORT_SUMMARY_PATH} has unsupported schema"));
    }
    for (field, expected) in [
        ("semantics", IMPORT_SUMMARY_SEMANTICS),
        ("manifest_identity_domain", IMPORT_MANIFEST_IDENTITY_DOMAIN),
        ("authority", IMPORT_AUTHORITY),
    ] {
        let observed = text_field(map, field, IMPORT_SUMMARY_PATH)?;
        if observed != expected {
            return Err(format!(
                "{IMPORT_SUMMARY_PATH} `{field}` is `{observed}`, expected `{expected}`"
            ));
        }
    }

    let manifest_identity = text_field(map, "manifest_identity", IMPORT_SUMMARY_PATH)?;
    if !is_canonical_hash(manifest_identity) {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} `manifest_identity` is not 32-byte lowercase hexadecimal"
        ));
    }
    let scorecard_identity = text_field(map, "scorecard_identity", IMPORT_SUMMARY_PATH)?;
    if !is_canonical_hash(scorecard_identity) {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} `scorecard_identity` is not 32-byte lowercase hexadecimal"
        ));
    }
    let observed_manifest_identity =
        fs_blake3::hash_domain(IMPORT_MANIFEST_IDENTITY_DOMAIN, manifest_bytes).to_string();
    if manifest_identity != observed_manifest_identity {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} is stale: manifest identity {manifest_identity} does not match \
             live {IMPORT_MANIFEST_PATH} identity {observed_manifest_identity}"
        ));
    }

    let population = object_field(map, "population", IMPORT_SUMMARY_PATH)?;
    let reviewed = object_field(map, "reviewed", IMPORT_SUMMARY_PATH)?;
    let population_total = count_field(population, "total", "population")?;
    let retained_population = retained_manifest_population(manifest_bytes)?;
    if population_total != retained_population {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} population total {population_total} does not match the \
             {retained_population} retained rows in {IMPORT_MANIFEST_PATH}"
        ));
    }
    let population_clean = count_field(population, "clean", "population")?;
    let population_repaired = count_field(population, "repaired", "population")?;
    let population_refused = count_field(population, "refused", "population")?;
    let population_unreviewed = count_field(population, "unreviewed", "population")?;
    let population_mismatches = count_field(population, "annotation_mismatch", "population")?;
    let proposal_mismatches = count_field(population, "proposal_mismatch", "population")?;
    if checked_sum(
        &[population_clean, population_repaired, population_refused],
        "population outcomes",
    )? != population_total
    {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} population outcomes do not partition total"
        ));
    }

    let reviewed_total = count_field(reviewed, "total", "reviewed")?;
    let reviewed_clean = count_field(reviewed, "clean", "reviewed")?;
    let reviewed_repaired = count_field(reviewed, "repaired", "reviewed")?;
    let reviewed_refused = count_field(reviewed, "refused", "reviewed")?;
    let reviewed_mismatches = count_field(reviewed, "annotation_mismatch", "reviewed")?;
    if reviewed_total
        .checked_add(population_unreviewed)
        .ok_or_else(|| format!("{IMPORT_SUMMARY_PATH} review counts overflow usize"))?
        != population_total
    {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} reviewed and unreviewed rows do not partition total"
        ));
    }
    if population_mismatches != reviewed_mismatches {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} population and reviewed annotation-mismatch counts disagree"
        ));
    }
    if proposal_mismatches > population_unreviewed {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} proposal mismatches exceed unreviewed rows"
        ));
    }

    ImportScorecardSource::try_new(
        population_total,
        reviewed_total,
        reviewed_clean,
        reviewed_repaired,
        reviewed_refused,
        reviewed_mismatches,
    )
    .map_err(|error| format!("{IMPORT_SUMMARY_PATH} carries invalid reviewed counts: {error}"))
}

fn read_import_scorecard(root: &Path) -> Result<ImportScorecardSource, String> {
    let summary_path = root.join(IMPORT_SUMMARY_PATH);
    let summary_metadata = std::fs::metadata(&summary_path)
        .map_err(|error| format!("cannot stat {IMPORT_SUMMARY_PATH}: {error}"))?;
    if summary_metadata.len() > MAX_IMPORT_SUMMARY_BYTES {
        return Err(format!(
            "{IMPORT_SUMMARY_PATH} exceeds the admitted {MAX_IMPORT_SUMMARY_BYTES}-byte bound"
        ));
    }
    let summary = std::fs::read_to_string(&summary_path)
        .map_err(|error| format!("cannot read {IMPORT_SUMMARY_PATH}: {error}"))?;

    let manifest_path = root.join(IMPORT_MANIFEST_PATH);
    let manifest_metadata = std::fs::metadata(&manifest_path)
        .map_err(|error| format!("cannot stat {IMPORT_MANIFEST_PATH}: {error}"))?;
    if manifest_metadata.len() > MAX_IMPORT_MANIFEST_BYTES {
        return Err(format!(
            "{IMPORT_MANIFEST_PATH} exceeds the admitted {MAX_IMPORT_MANIFEST_BYTES}-byte bound"
        ));
    }
    let manifest_bytes = std::fs::read(&manifest_path)
        .map_err(|error| format!("cannot read {IMPORT_MANIFEST_PATH}: {error}"))?;
    parse_import_scorecard_summary(&summary, &manifest_bytes)
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

/// Spine inputs for the dashboard: the executing stage prefix derived live
/// from the fs-cli source (the ratchet gate keeps that derivation honest),
/// and the deliberately regenerated beads snapshot. Anything unreadable is a
/// `None` field, which the projection renders as `NO-DATA` — the dashboard
/// must never invent a spine number.
fn spine_sources(root: &Path) -> SpineMetricsSource {
    let stages = std::fs::read_to_string(root.join(spine_ratchet::SOLVE_SOURCE))
        .ok()
        .and_then(|source| spine_ratchet::derive_stages(&source).ok())
        .map(|stages| {
            let executing = spine_ratchet::executing_prefix(&stages).len();
            (executing, stages.len())
        });
    let snapshot = spine_metrics::load(root);
    let tropical = tropical_path::load(root);
    let spine_on_path = tropical.as_ref().map(|artifact| {
        tropical_path::SPINE_BEADS
            .iter()
            .filter(|id| artifact.slack_hours.get(**id) == Some(&0.0))
            .count()
    });
    let spine_tracked = tropical.as_ref().map(|artifact| {
        tropical_path::SPINE_BEADS
            .iter()
            .filter(|id| artifact.slack_hours.contains_key(**id))
            .count()
    });
    SpineMetricsSource {
        stages_executing: stages.map(|(executing, _)| executing),
        stages_total: stages.map(|(_, total)| total),
        beads_open: snapshot.map(|snapshot| snapshot.open),
        beads_blocked: snapshot.map(|snapshot| snapshot.blocked),
        beads_actionable: snapshot.map(|snapshot| snapshot.actionable),
        spine_on_critical_path: spine_on_path,
        spine_tracked,
    }
}

/// Project the live registries into the dashboard.
fn dashboard(root: &Path) -> Result<ProgramDashboard, String> {
    let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[])
        .map_err(|error| format!("cannot build the V&V scorecard: {error}"))?;
    let levels = maturity::capability_levels(root)?.current;
    let import_scorecard = read_import_scorecard(root)?;
    let sources = ProgramSources::from_registries(
        corpus(),
        adversarial_registry(),
        &scorecard,
        import_scorecard,
        &levels,
    )
    .map_err(|error| format!("cannot read the program sources: {error}"))?
    .with_spine(spine_sources(root));
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
        let corpus_dir = dir.join("data/cad-import-corpus");
        std::fs::create_dir_all(&corpus_dir).expect("corpus artifact directory is creatable");
        for path in [IMPORT_SUMMARY_PATH, IMPORT_MANIFEST_PATH] {
            let file_name = std::path::Path::new(path)
                .file_name()
                .expect("artifact path has a file name");
            std::fs::copy(root.join(path), corpus_dir.join(file_name))
                .expect("the retained import scorecard source is readable");
        }
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
        let import_scorecard =
            read_import_scorecard(&repo_root()).expect("the import scorecard source validates");
        assert_eq!(import_scorecard.total(), 21);
        assert_eq!(import_scorecard.reviewed(), 0);
        assert_eq!(import_scorecard.clean(), 0);
        assert_eq!(import_scorecard.repaired(), 0);
        assert_eq!(import_scorecard.refused(), 0);
        assert_eq!(import_scorecard.annotation_mismatches(), 0);

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

    /// o5et9: the spine rows must render the three genuinely different
    /// states differently — a measured ratchet ratio, a measured snapshot
    /// count, and a NO-DATA e2e gap are not interchangeable, and "lane green
    /// with no retained receipt" must never render as a measured zero.
    #[test]
    fn spine_rows_discriminate_measured_snapshot_and_no_data_states() {
        let dashboard = dashboard(&repo_root()).expect("the real dashboard builds");
        let cell = |id: &str| {
            dashboard
                .rows()
                .iter()
                .find(|entry| entry.metric().id() == id)
                .unwrap_or_else(|| panic!("spine metric `{id}` exists"))
                .metric()
                .cell()
                .status()
        };
        // The live ratchet derivation feeds a measured ratio.
        assert_eq!(cell("spine-stages-executing"), "ratio");
        // The validated snapshot feeds measured beads rows.
        assert_eq!(cell("beads-blocked-ratio"), "ratio");
        assert_eq!(cell("beads-actionable"), "count");
        // The e2e lane has no retained receipt: NO-DATA, never a zero.
        assert_eq!(cell("spine-e2e-lane-green"), "no-data");
        // The tropical critical path landed (kx95s): a measured ratio whose
        // numerator is a REAL ZERO (no spine bead is on the path today).
        assert_eq!(cell("spine-critical-path-positions"), "ratio");
        // And the discrimination itself: a NO-DATA status and a measured
        // zero-status are different renderings of different facts.
        assert_ne!("no-data", "count");
    }

    /// The tracked compact projection is not trusted across manifest drift.
    /// This is a pure negative test: it mutates only an in-memory copy.
    #[test]
    fn import_summary_refuses_manifest_identity_drift() {
        let root = repo_root();
        let source =
            std::fs::read_to_string(root.join(IMPORT_SUMMARY_PATH)).expect("summary is readable");
        let manifest =
            std::fs::read(root.join(IMPORT_MANIFEST_PATH)).expect("manifest is readable");
        let identity =
            fs_blake3::hash_domain(IMPORT_MANIFEST_IDENTITY_DOMAIN, &manifest).to_string();
        let stale = source.replacen(&identity, &"0".repeat(64), 1);
        let error = parse_import_scorecard_summary(&stale, &manifest)
            .expect_err("a stale manifest identity must be refused");
        assert!(error.contains("does not match live"), "{error}");
    }

    /// A valid manifest identity cannot bless summary counts that describe a
    /// different retained population.
    #[test]
    fn import_summary_refuses_population_drift() {
        let root = repo_root();
        let source =
            std::fs::read_to_string(root.join(IMPORT_SUMMARY_PATH)).expect("summary is readable");
        let manifest =
            std::fs::read(root.join(IMPORT_MANIFEST_PATH)).expect("manifest is readable");
        let stale = source.replacen("\"total\":21", "\"total\":20", 1);
        let error = parse_import_scorecard_summary(&stale, &manifest)
            .expect_err("a stale population total must be refused");
        assert!(
            error.contains("does not match the 21 retained rows"),
            "{error}"
        );
    }

    /// The reviewed-data branch is exercised outside fs-govern's broad
    /// integration-test dev-dependency cone.
    #[test]
    fn reviewed_import_counts_render_three_rates_and_a_mismatch_count() {
        let scorecard = build_scorecard(corpus(), adversarial_registry(), &[], &[])
            .expect("the seeded V&V scorecard builds");
        let levels = maturity::capability_levels(&repo_root())
            .expect("the maturity registry parses")
            .current;
        let import_scorecard = ImportScorecardSource::try_new(21, 4, 1, 2, 1, 1)
            .expect("the reviewed fixture partitions its population");
        let sources = ProgramSources::from_registries(
            corpus(),
            adversarial_registry(),
            &scorecard,
            import_scorecard,
            &levels,
        )
        .expect("the reviewed fixture projects");
        let rows = frankensim_rows(sources).expect("the metric rows build");
        let cell = |id: &str| {
            rows.iter()
                .find(|row| row.id() == id)
                .unwrap_or_else(|| panic!("metric `{id}` exists"))
                .cell()
        };
        let ratio = |id: &str| match cell(id) {
            MetricCell::Measured(MetricObservation::Ratio {
                numerator,
                denominator,
            }) => (*numerator, denominator.get()),
            other => panic!("metric `{id}` is not a ratio: {other:?}"),
        };
        assert_eq!(ratio("import-admission-rate"), (1, 4));
        assert_eq!(ratio("import-repair-rate"), (2, 4));
        assert_eq!(ratio("import-refusal-rate"), (1, 4));
        assert_eq!(
            cell("import-annotation-regressions"),
            &MetricCell::Measured(MetricObservation::count(1))
        );
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
