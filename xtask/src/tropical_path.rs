//! Certified tropical critical path over the open-bead DAG (bead frankensim-kx95s).
//!
//! Why this exists. The 2026-07-28 reality check measured the swarm aiming
//! away from the product spine: bv ranks by topological centrality (PageRank,
//! betweenness, HITS), which answers "what is well-connected" — not the
//! delivery question "given estimates and the dependency DAG, which beads set
//! the makespan, and which buy nothing?" `fs-tropical` is the in-house
//! max-plus instrument built for exactly that question; this module points it
//! at the project's own tracker.
//!
//! Projection (recorded in the artifact's provenance, never silent):
//!
//! - nodes: every issue whose status is not `closed`;
//! - edges: `blocks` and `parent-child` dependencies where both ends are
//!   open, oriented completion-order (blocker → blocked, child → parent);
//! - weights: `estimated_minutes` as hours, falling back to
//!   [`DEFAULT_ESTIMATE_MINUTES`] when absent — defaulted weights are COUNTED
//!   in the artifact, because an unestimated bead silently weighted zero
//!   would corrupt the makespan;
//! - cycles: `fs_tropical` refuses them (`TropicalError::Cyclic`) and so do
//!   we — the projection never silently linearises a cyclic tracker.
//!
//! The artifact is deliberately regenerated, like the spine ratchet and the
//! beads snapshot: the check gate validates it and renders tracker movement
//! as a visible `stale-snapshot`-style note, not a wedge. A live CYCLE,
//! however, is a violation: the critical path is undefined until it is
//! broken, and pretending otherwise is how the swarm steers by fiction.

use std::collections::BTreeMap;
use std::path::Path;

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation, fnv1a64};

pub(crate) const CHECK: &str = "tropical-critical-path";
const ARTIFACT_PATH: &str = "tropical-critical-path.json";
const ISSUES_PATH: &str = ".beads/issues.jsonl";
const SCHEMA: &str = "frankensim-tropical-critical-path-v1";
const MAX_ISSUES_BYTES: u64 = 128 * 1024 * 1024;
const MAX_ARTIFACT_BYTES: u64 = 4 * 1024 * 1024;
const MAX_TRACKER_STRING_BYTES: usize = 4 * 1024 * 1024;
/// Weight given to a bead with no `estimated_minutes`. One workday: large
/// enough that an unestimated bead cannot vanish off the makespan, small
/// enough that a defaulted chain cannot dominate a genuinely estimated one.
/// Recorded in the artifact; change it deliberately.
const DEFAULT_ESTIMATE_MINUTES: u64 = 480;

/// The five beads the 2026-07-28 reality check named as standing between the
/// project and its first end-to-end answer. Their slack is the steering
/// number this artifact exists to produce.
pub(crate) const SPINE_BEADS: [&str; 5] = [
    "frankensim-frn2i",
    "frankensim-s93ej",
    "frankensim-s2l9v",
    "frankensim-extreal-program-f85xj.6.9",
    "frankensim-extreal-program-f85xj.6.10",
];

/// One open bead as a DAG node.
#[derive(Debug, Clone)]
struct BeadNode {
    id: String,
    estimate_hours: f64,
    weight_defaulted: bool,
    /// Completion-order successors: beads this one blocks, plus its parent.
    successors: Vec<String>,
}

/// The computed projection, before rendering.
#[derive(Debug, Clone)]
pub(crate) struct Projection {
    issues_fnv: u64,
    open: usize,
    edges_blocks: usize,
    edges_parent_child: usize,
    defaulted_weights: usize,
    makespan_hours: f64,
    path_is_unique: bool,
    /// Critical path as (bead id, latency hours), source → sink.
    critical_path: Vec<(String, f64)>,
    /// Per-bead slack hours over every open bead.
    slack_hours: BTreeMap<String, f64>,
    /// Slack for each of [`SPINE_BEADS`], when the bead is open.
    spine_positions: BTreeMap<String, f64>,
    /// bv's topological top-10 at generation time, with our slack alongside:
    /// (id, bv score, slack hours). Empty when bv was unavailable.
    bv_top: Vec<(String, f64, Option<f64>)>,
    bv_available: bool,
}

fn violation(detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: ARTIFACT_PATH.to_string(),
        detail: detail.into(),
    }
}

/// Parse tracker bytes into open-bead nodes. Pure so tests can fixture it.
fn parse_nodes(text: &str) -> Result<Vec<BeadNode>, String> {
    let mut nodes = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parsed = JsonParser::with_string_limit(line, MAX_TRACKER_STRING_BYTES)
            .finish()
            .map_err(|error| {
                format!(
                    "{ISSUES_PATH} line {} is not valid JSON: {error}",
                    index + 1
                )
            })?;
        let JsonValue::Object(map) = &parsed else {
            return Err(format!(
                "{ISSUES_PATH} line {} is not a JSON object",
                index + 1
            ));
        };
        let string_field = |key: &str| match map.get(key) {
            Some(JsonValue::String(value)) => Some(value.clone()),
            _ => None,
        };
        let status = string_field("status")
            .ok_or_else(|| format!("{ISSUES_PATH} line {} has no string `status`", index + 1))?;
        if status == "closed" {
            continue;
        }
        let id = string_field("id")
            .ok_or_else(|| format!("{ISSUES_PATH} line {} has no string `id`", index + 1))?;
        let (estimate_hours, weight_defaulted) = match map.get("estimated_minutes") {
            Some(JsonValue::Number(raw)) => {
                let minutes = raw.parse::<f64>().map_err(|error| {
                    format!(
                        "{ISSUES_PATH} line {index}: bad estimated_minutes: {error}",
                        index = index + 1
                    )
                })?;
                if !(minutes.is_finite() && minutes >= 0.0) {
                    return Err(format!(
                        "{ISSUES_PATH} line {} has non-finite or negative estimated_minutes",
                        index + 1
                    ));
                }
                (minutes / 60.0, false)
            }
            _ => (DEFAULT_ESTIMATE_MINUTES as f64 / 60.0, true),
        };
        let mut successors = Vec::new();
        if let Some(JsonValue::Array(deps)) = map.get("dependencies") {
            for dep in deps {
                let JsonValue::Object(dep) = dep else {
                    return Err(format!(
                        "{ISSUES_PATH} line {} has a non-object dependency",
                        index + 1
                    ));
                };
                let dep_type = match dep.get("type") {
                    Some(JsonValue::String(t)) => t.as_str(),
                    _ => continue,
                };
                let target = match dep.get("depends_on_id") {
                    Some(JsonValue::String(target)) => target.clone(),
                    _ => continue,
                };
                match dep_type {
                    // A blocks dependency `self depends_on target` completes
                    // target first: edge target → self.
                    "blocks" | "parent-child" => successors.push(format!("{dep_type}:{target}")),
                    _ => {}
                }
            }
        }
        nodes.push(BeadNode {
            id,
            estimate_hours,
            weight_defaulted,
            successors,
        });
    }
    if nodes.is_empty() {
        return Err(format!(
            "{ISSUES_PATH} parsed to zero open issues; a gate with no inputs cannot report a path"
        ));
    }
    Ok(nodes)
}

/// Project the tracker bytes into a certified critical path.
///
/// Errors on cycles (never silently linearised) and on malformed input.
fn project(text: &str) -> Result<Projection, String> {
    let nodes = parse_nodes(text)?;
    let index_of: BTreeMap<&str, usize> = nodes
        .iter()
        .enumerate()
        .map(|(index, node)| (node.id.as_str(), index))
        .collect();
    let mut edges_blocks = 0usize;
    let mut edges_parent_child = 0usize;
    let mut dag = fs_tropical::TaskDag::new(nodes.iter().map(|node| node.estimate_hours).collect());
    for (to_index, node) in nodes.iter().enumerate() {
        for successor in &node.successors {
            let (kind, target) = successor
                .split_once(':')
                .ok_or_else(|| format!("malformed successor encoding {successor}"))?;
            let Some(&other_index) = index_of.get(target) else {
                continue; // target closed or absent: no completion constraint
            };
            if other_index == to_index {
                continue; // a self-loop carries no schedule constraint
            }
            // Completion order differs by dependency kind:
            // - `blocks`: the blocker finishes first (blocker -> blocked);
            // - `parent-child`: the CHILD finishes first, because a parent
            //   closes when its children close (child -> parent). Reversing
            //   the second manufactures cycles out of ordinary epic trees.
            let (from_index, to) = match kind {
                "blocks" => {
                    edges_blocks += 1;
                    (other_index, to_index)
                }
                "parent-child" => {
                    edges_parent_child += 1;
                    (to_index, other_index)
                }
                _ => continue,
            };
            dag = dag.with_edge(from_index, to);
        }
    }
    let path = dag.critical_path().map_err(|error| match error {
        fs_tropical::TropicalError::Cyclic => format!(
            "the open-bead dependency graph contains a CYCLE; the critical path is undefined \
             and this projection refuses to silently linearise it — break the cycle in the \
             tracker deliberately ({error})"
        ),
        other => format!("cannot compute the bead critical path: {other}"),
    })?;
    let slack_hours: BTreeMap<String, f64> = nodes
        .iter()
        .enumerate()
        .map(|(index, node)| (node.id.clone(), path.slack[index]))
        .collect();
    let spine_positions = SPINE_BEADS
        .iter()
        .filter_map(|id| {
            slack_hours
                .get(*id)
                .map(|slack| ((*id).to_string(), *slack))
        })
        .collect();
    Ok(Projection {
        issues_fnv: fnv1a64(text.as_bytes()),
        open: nodes.len(),
        edges_blocks,
        edges_parent_child,
        defaulted_weights: nodes.iter().filter(|node| node.weight_defaulted).count(),
        makespan_hours: path.makespan,
        path_is_unique: path.path_is_unique,
        critical_path: path
            .path
            .iter()
            .map(|&index| (nodes[index].id.clone(), nodes[index].estimate_hours))
            .collect(),
        slack_hours,
        spine_positions,
        bv_top: Vec::new(),
        bv_available: false,
    })
}

/// Attach bv's topological top-10 when bv is runnable, so the disagreement
/// between centrality and schedule slack is retained on the record. bv is a
/// generate-time input only; the check gate never shells out.
fn attach_bv_comparison(root: &Path, projection: &mut Projection) {
    let output = std::process::Command::new("bv")
        .arg("--robot-triage")
        .current_dir(root)
        .env("BV_NO_COLOR", "1")
        .output();
    let Ok(output) = output else {
        return; // bv not installed here: comparison is NO-DATA, path stands alone
    };
    if !output.status.success() {
        return;
    }
    let Ok(text) = String::from_utf8(output.stdout) else {
        return;
    };
    let Ok(parsed) = JsonParser::with_string_limit(&text, MAX_TRACKER_STRING_BYTES).finish() else {
        return;
    };
    let JsonValue::Object(map) = &parsed else {
        return;
    };
    let Some(JsonValue::Object(triage)) = map.get("triage") else {
        return;
    };
    let Some(JsonValue::Array(recommendations)) = triage.get("recommendations") else {
        return;
    };
    for recommendation in recommendations.iter().take(10) {
        let JsonValue::Object(recommendation) = recommendation else {
            continue;
        };
        let (Some(JsonValue::String(id)), Some(JsonValue::Number(score))) =
            (recommendation.get("id"), recommendation.get("score"))
        else {
            continue;
        };
        let Ok(score) = score.parse::<f64>() else {
            continue;
        };
        projection
            .bv_top
            .push((id.clone(), score, projection.slack_hours.get(id).copied()));
    }
    projection.bv_available = !projection.bv_top.is_empty();
}

fn render_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value}")
    } else {
        "null".to_string()
    }
}

/// Render the canonical artifact bytes.
fn render(projection: &Projection) -> String {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"schema\": \"{SCHEMA}\",\n"));
    out.push_str("  \"bead\": \"frankensim-kx95s\",\n");
    out.push_str(&format!(
        "  \"source\": {{\n    \"issues_fnv1a64\": \"{:016x}\",\n    \"open\": {},\n    \
         \"edges_blocks\": {},\n    \"edges_parent_child\": {},\n    \"defaulted_weights\": {},\n    \
         \"default_estimate_minutes\": {}\n  }},\n",
        projection.issues_fnv,
        projection.open,
        projection.edges_blocks,
        projection.edges_parent_child,
        projection.defaulted_weights,
        DEFAULT_ESTIMATE_MINUTES
    ));
    out.push_str(&format!(
        "  \"makespan_hours\": {},\n  \"path_is_unique\": {},\n",
        render_f64(projection.makespan_hours),
        projection.path_is_unique
    ));
    out.push_str("  \"critical_path\": [\n");
    for (index, (id, latency)) in projection.critical_path.iter().enumerate() {
        let comma = if index + 1 == projection.critical_path.len() {
            ""
        } else {
            ","
        };
        out.push_str(&format!(
            "    {{\"id\": \"{id}\", \"latency_hours\": {}}}{comma}\n",
            render_f64(*latency)
        ));
    }
    out.push_str("  ],\n");
    out.push_str("  \"spine_positions\": {\n");
    let spine_len = projection.spine_positions.len();
    for (index, (id, slack)) in projection.spine_positions.iter().enumerate() {
        let comma = if index + 1 == spine_len { "" } else { "," };
        out.push_str(&format!(
            "    \"{id}\": {{\"slack_hours\": {}, \"on_critical_path\": {}}}{comma}\n",
            render_f64(*slack),
            *slack == 0.0
        ));
    }
    out.push_str("  },\n");
    out.push_str("  \"slack_hours\": {\n");
    let slack_len = projection.slack_hours.len();
    for (index, (id, slack)) in projection.slack_hours.iter().enumerate() {
        let comma = if index + 1 == slack_len { "" } else { "," };
        out.push_str(&format!("    \"{id}\": {}{comma}\n", render_f64(*slack)));
    }
    out.push_str("  },\n");
    if projection.bv_available {
        out.push_str(
            "  \"bv_comparison\": {\n    \"status\": \"measured\",\n    \"bv_top10\": [\n",
        );
        let top_len = projection.bv_top.len();
        for (index, (id, score, slack)) in projection.bv_top.iter().enumerate() {
            let comma = if index + 1 == top_len { "" } else { "," };
            let slack = slack.map_or_else(|| "null".to_string(), render_f64);
            out.push_str(&format!(
                "      {{\"id\": \"{id}\", \"bv_score\": {}, \"slack_hours\": {slack}}}{comma}\n",
                render_f64(*score)
            ));
        }
        out.push_str("    ]\n  },\n");
    } else {
        out.push_str(
            "  \"bv_comparison\": {\n    \"status\": \"no-bv-available\",\n    \"bv_top10\": []\n  },\n",
        );
    }
    out.push_str(
        "  \"no_claim\": \"weights are bead estimates with a recorded 480-minute default for \
         unestimated beads; the makespan is a schedule bound over ESTIMATES, not a promise. \
         slack 0 means on the critical path under those weights. bv's centrality ranking is \
         retained beside slack where bv was runnable; disagreement between the two is the \
         finding, not an error\"\n",
    );
    out.push_str("}\n");
    out
}

fn json_number(map: &BTreeMap<String, JsonValue>, key: &str) -> Result<f64, String> {
    match map.get(key) {
        Some(JsonValue::Number(raw)) => raw
            .parse::<f64>()
            .map_err(|error| format!("artifact `{key}` is not a number: {error}")),
        _ => Err(format!("artifact has no number `{key}`")),
    }
}

fn json_usize(map: &BTreeMap<String, JsonValue>, key: &str) -> Result<usize, String> {
    let value = json_number(map, key)?;
    if value < 0.0 || value.fract() != 0.0 || value > (usize::MAX as f64) {
        return Err(format!("artifact `{key}` is not a count: {value}"));
    }
    Ok(value as usize)
}

/// The validated artifact contents.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TropicalArtifact {
    pub issues_fnv: u64,
    pub open: usize,
    pub defaulted_weights: usize,
    pub makespan_hours: f64,
    pub path_is_unique: bool,
    pub critical_path: Vec<(String, f64)>,
    pub spine_positions: BTreeMap<String, f64>,
    pub slack_hours: BTreeMap<String, f64>,
}

/// Parse and fully validate artifact text.
fn parse_artifact(text: &str) -> Result<TropicalArtifact, String> {
    let parsed = JsonParser::with_string_limit(text, MAX_TRACKER_STRING_BYTES)
        .finish()
        .map_err(|error| format!("{ARTIFACT_PATH} is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err(format!("{ARTIFACT_PATH} is not a JSON object"));
    };
    match map.get("schema") {
        Some(JsonValue::String(schema)) if schema == SCHEMA => {}
        Some(JsonValue::String(schema)) => {
            return Err(format!(
                "{ARTIFACT_PATH} schema is `{schema}`, expected `{SCHEMA}`; refusing to read a \
                 foreign artifact"
            ));
        }
        _ => return Err(format!("{ARTIFACT_PATH} has no schema string")),
    }
    let source = match map.get("source") {
        Some(JsonValue::Object(source)) => source,
        _ => return Err(format!("{ARTIFACT_PATH} has no `source` object")),
    };
    let issues_fnv = match source.get("issues_fnv1a64") {
        Some(JsonValue::String(hex)) => u64::from_str_radix(hex, 16)
            .map_err(|error| format!("artifact `issues_fnv1a64` is not hex: {error}"))?,
        _ => return Err("artifact has no `issues_fnv1a64` string".to_string()),
    };
    let open = json_usize(source, "open")?;
    let defaulted_weights = json_usize(source, "defaulted_weights")?;
    if defaulted_weights > open {
        return Err(format!(
            "artifact is internally inconsistent: defaulted_weights ({defaulted_weights}) \
             exceeds open ({open})"
        ));
    }
    let makespan_hours = json_number(map, "makespan_hours")?;
    if !(makespan_hours.is_finite() && makespan_hours >= 0.0) {
        return Err("artifact makespan is not a finite non-negative number".to_string());
    }
    let path_is_unique = matches!(map.get("path_is_unique"), Some(JsonValue::Bool));
    let mut critical_path = Vec::new();
    match map.get("critical_path") {
        Some(JsonValue::Array(entries)) => {
            for entry in entries {
                let JsonValue::Object(entry) = entry else {
                    return Err("artifact critical_path has a non-object entry".to_string());
                };
                let id = match entry.get("id") {
                    Some(JsonValue::String(id)) => id.clone(),
                    _ => return Err("artifact critical_path entry has no id".to_string()),
                };
                let latency = json_number(entry, "latency_hours")?;
                critical_path.push((id, latency));
            }
        }
        _ => return Err("artifact has no `critical_path` array".to_string()),
    }
    let mut slack_hours = BTreeMap::new();
    match map.get("slack_hours") {
        Some(JsonValue::Object(slack)) => {
            for (id, value) in slack {
                match value {
                    JsonValue::Number(raw) => {
                        let parsed = raw.parse::<f64>().map_err(|error| {
                            format!("artifact slack for {id} is not a number: {error}")
                        })?;
                        if !(parsed.is_finite() && parsed >= 0.0) {
                            return Err(format!(
                                "artifact slack for {id} is not finite and non-negative"
                            ));
                        }
                        slack_hours.insert(id.clone(), parsed);
                    }
                    _ => return Err(format!("artifact slack for {id} is not numeric")),
                }
            }
        }
        _ => return Err("artifact has no `slack_hours` object".to_string()),
    }
    if slack_hours.len() != open {
        return Err(format!(
            "artifact is internally inconsistent: slack covers {} beads but source records \
             {open} open",
            slack_hours.len()
        ));
    }
    let mut spine_positions = BTreeMap::new();
    match map.get("spine_positions") {
        Some(JsonValue::Object(spine)) => {
            for (id, value) in spine {
                let JsonValue::Object(value) = value else {
                    return Err(format!("artifact spine position for {id} is not an object"));
                };
                let slack = json_number(value, "slack_hours")?;
                match slack_hours.get(id) {
                    Some(recorded) if recorded == &slack => {}
                    Some(recorded) => {
                        return Err(format!(
                            "artifact is internally inconsistent: spine slack for {id} is \
                             {slack} but the slack map says {recorded}"
                        ));
                    }
                    None => {
                        return Err(format!(
                            "artifact spine position {id} is absent from the slack map"
                        ));
                    }
                }
                spine_positions.insert(id.clone(), slack);
            }
        }
        _ => return Err("artifact has no `spine_positions` object".to_string()),
    }
    // The makespan must equal the critical path's latency sum, in path order.
    let path_sum: f64 = critical_path.iter().map(|(_, latency)| latency).sum();
    let tolerance = makespan_hours.max(1.0) * 1e-9;
    if (path_sum - makespan_hours).abs() > tolerance {
        return Err(format!(
            "artifact is internally inconsistent: critical path latencies sum to {path_sum} \
             but makespan is {makespan_hours}"
        ));
    }
    // Every path bead must carry zero slack in the slack map.
    for (id, _) in &critical_path {
        match slack_hours.get(id) {
            Some(slack) if *slack == 0.0 => {}
            Some(slack) => {
                return Err(format!(
                    "artifact is internally inconsistent: critical-path bead {id} has slack \
                     {slack}"
                ));
            }
            None => {
                return Err(format!(
                    "artifact critical-path bead {id} is absent from the slack map"
                ));
            }
        }
    }
    Ok(TropicalArtifact {
        issues_fnv,
        open,
        defaulted_weights,
        makespan_hours,
        path_is_unique,
        critical_path,
        spine_positions,
        slack_hours,
    })
}

fn read_bounded(root: &Path, relative: &str, max: u64) -> Result<String, String> {
    let bytes = std::fs::read(root.join(relative))
        .map_err(|error| format!("{relative} unreadable: {error}"))?;
    if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > max {
        return Err(format!("{relative} exceeds the {max}-byte gate bound"));
    }
    String::from_utf8(bytes).map_err(|error| format!("{relative} is not UTF-8: {error}"))
}

/// Load the tracked artifact for the dashboard; absence or invalidity is
/// NO-DATA there and a violation here.
pub(crate) fn load(root: &Path) -> Option<TropicalArtifact> {
    let text = read_bounded(root, ARTIFACT_PATH, MAX_ARTIFACT_BYTES).ok()?;
    parse_artifact(&text).ok()
}

/// Regenerate the artifact from the live tracker (and bv, when runnable).
pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let issues = read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES)?;
    let mut projection = project(&issues)?;
    attach_bv_comparison(root, &mut projection);
    std::fs::write(root.join(ARTIFACT_PATH), render(&projection))
        .map_err(|error| format!("cannot write {ARTIFACT_PATH}: {error}"))
}

/// The standing gate: validate the artifact, refuse a live cycle, and render
/// tracker movement as a visible note rather than a wedge.
pub(crate) fn check(root: &Path) -> (Vec<Violation>, Vec<PolicyNote>) {
    let text = match read_bounded(root, ARTIFACT_PATH, MAX_ARTIFACT_BYTES) {
        Ok(text) => text,
        Err(error) => {
            return (
                vec![violation(format!(
                    "{error}; a missing or unreadable artifact is not a critical path — run \
                     `cargo run -p xtask -- generate-tropical-path`"
                ))],
                Vec::new(),
            );
        }
    };
    let artifact = match parse_artifact(&text) {
        Ok(artifact) => artifact,
        Err(error) => return (vec![violation(error)], Vec::new()),
    };
    let mut notes = Vec::new();
    match read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES) {
        Ok(issues) => {
            // Re-project the LIVE tracker: a cycle is a violation even when
            // the retained artifact predates it, because the certified path
            // is undefined until the cycle breaks.
            if let Err(error) = project(&issues) {
                if error.contains("CYCLE") {
                    return (vec![violation(format!("live tracker: {error}"))], notes);
                }
                // A non-cycle live parse failure is a tracker-format problem
                // the snapshot gate already owns; stay a note here.
                notes.push(PolicyNote {
                    check: CHECK,
                    crate_name: ARTIFACT_PATH.to_string(),
                    verdict: "live-unreadable",
                    detail: format!("the live tracker could not be re-projected: {error}"),
                });
            }
            let live_fnv = fnv1a64(issues.as_bytes());
            if live_fnv == artifact.issues_fnv {
                notes.push(PolicyNote {
                    check: CHECK,
                    crate_name: ARTIFACT_PATH.to_string(),
                    verdict: "current",
                    detail: format!(
                        "critical path is current (makespan {}h over {} open beads, {} \
                         defaulted weights)",
                        artifact.makespan_hours, artifact.open, artifact.defaulted_weights
                    ),
                });
            } else {
                notes.push(PolicyNote {
                    check: CHECK,
                    crate_name: ARTIFACT_PATH.to_string(),
                    verdict: "stale-snapshot",
                    detail: format!(
                        "the tracker moved since the critical path was computed (recorded \
                         {:016x}, live {:016x}); regenerate with `cargo run -p xtask -- \
                         generate-tropical-path` when the move should be on the record",
                        artifact.issues_fnv, live_fnv
                    ),
                });
            }
        }
        Err(error) => {
            return (
                vec![violation(format!(
                    "{error}; cannot compare the artifact against the tracker, refusing to \
                     certify currency"
                ))],
                notes,
            );
        }
    }
    (Vec::new(), notes)
}

#[cfg(test)]
mod tests {
    //! G0/G3: projection laws, slack semantics, cycle refusal, artifact
    //! anti-silent-disable, and live validation.

    use super::*;

    fn issue(id: &str, status: &str, estimate: Option<u64>, deps: &[(&str, &str)]) -> String {
        let deps: Vec<String> = deps
            .iter()
            .map(|(kind, target)| {
                format!(
                    "{{\"issue_id\":\"{id}\",\"depends_on_id\":\"{target}\",\"type\":\"{kind}\"}}"
                )
            })
            .collect();
        let estimate = estimate.map_or_else(|| "null".to_string(), |e| e.to_string());
        format!(
            "{{\"id\":\"{id}\",\"status\":\"{status}\",\"estimated_minutes\":{estimate},\"dependencies\":[{}]}}",
            deps.join(",")
        )
    }

    #[test]
    fn g0_the_critical_path_is_the_longest_chain_not_the_loudest_bead() {
        // a(2h) -> b(10h) -> c(1h) is the makespan chain at 13h; d(8h) alone
        // has 5h of slack despite being the second-largest single task.
        let text = [
            issue("a", "open", Some(120), &[]),
            issue("b", "open", Some(600), &[("blocks", "a")]),
            issue("c", "open", Some(60), &[("blocks", "b")]),
            issue("d", "open", Some(480), &[]),
        ]
        .join("\n");
        let projection = project(&text).expect("projection");
        assert_eq!(projection.makespan_hours, 13.0);
        let path: Vec<&str> = projection
            .critical_path
            .iter()
            .map(|(id, _)| id.as_str())
            .collect();
        assert_eq!(path, vec!["a", "b", "c"]);
        assert_eq!(projection.slack_hours["d"], 5.0);
        assert_eq!(projection.slack_hours["b"], 0.0);
        assert_eq!(projection.defaulted_weights, 0);
    }

    #[test]
    fn g0_unestimated_beads_get_the_recorded_default_not_zero() {
        let text = [
            issue("a", "open", None, &[]),
            issue("b", "open", Some(60), &[("blocks", "a")]),
        ]
        .join("\n");
        let projection = project(&text).expect("projection");
        assert_eq!(projection.defaulted_weights, 1);
        // 480-minute default + 60 minutes = 9 hours of makespan; a silent
        // zero would have reported 1.
        assert_eq!(projection.makespan_hours, 9.0);
    }

    #[test]
    fn g0_closed_blockers_and_closed_nodes_drop_out() {
        let text = [
            issue("a", "open", Some(60), &[]),
            issue("b", "closed", Some(6000), &[]),
            issue("c", "open", Some(60), &[("blocks", "b")]),
        ]
        .join("\n");
        let projection = project(&text).expect("projection");
        assert_eq!(projection.open, 2);
        assert_eq!(projection.makespan_hours, 1.0);
        assert!(!projection.slack_hours.contains_key("b"));
    }

    #[test]
    fn g0_a_cycle_refuses_rather_than_linearising() {
        let text = [
            issue("a", "open", Some(60), &[("blocks", "c")]),
            issue("b", "open", Some(60), &[("blocks", "a")]),
            issue("c", "open", Some(60), &[("blocks", "b")]),
        ]
        .join("\n");
        let error = project(&text).expect_err("cycle must refuse");
        assert!(error.contains("CYCLE"), "{error}");
    }

    #[test]
    fn g0_parent_child_edges_chain_completion() {
        let text = [
            issue("parent", "open", Some(0), &[]),
            issue("child", "open", Some(300), &[("parent-child", "parent")]),
        ]
        .join("\n");
        let projection = project(&text).expect("projection");
        assert_eq!(projection.edges_parent_child, 1);
        assert_eq!(projection.makespan_hours, 5.0);
        let path: Vec<&str> = projection
            .critical_path
            .iter()
            .map(|(id, _)| id.as_str())
            .collect();
        assert_eq!(path, vec!["child", "parent"]);
    }

    #[test]
    fn g0_artifact_round_trips_and_stays_canonical() {
        let text = [
            issue("a", "open", Some(120), &[]),
            issue("b", "open", Some(240), &[("blocks", "a")]),
        ]
        .join("\n");
        let projection = project(&text).expect("projection");
        let rendered = render(&projection);
        let parsed = parse_artifact(&rendered).expect("round trip");
        assert_eq!(parsed.open, 2);
        assert_eq!(parsed.makespan_hours, 6.0);
        assert_eq!(parsed.critical_path.len(), 2);
        assert_eq!(parsed.slack_hours["a"], 0.0);
    }

    #[test]
    fn g0_foreign_or_inconsistent_artifacts_refuse() {
        assert!(parse_artifact("").is_err());
        assert!(parse_artifact("{\"schema\":\"something-else\"}").is_err());
        // Slack map / open count mismatch.
        let text = [issue("a", "open", Some(60), &[])].join("\n");
        let projection = project(&text).expect("projection");
        let mut rendered = render(&projection);
        rendered = rendered.replace("\"open\": 1", "\"open\": 2");
        let error = parse_artifact(&rendered).expect_err("count mismatch must refuse");
        assert!(error.contains("inconsistent"), "{error}");
        // A critical-path bead with nonzero slack is impossible.
        let mut rendered = render(&projection);
        rendered = rendered.replace("\"frankensim-a\"", "\"frankensim-a\"");
        let tampered = rendered.replace("\"a\": 0", "\"a\": 5");
        if tampered != rendered {
            let error = parse_artifact(&tampered).expect_err("slack tamper must refuse");
            assert!(
                error.contains("inconsistent") || error.contains("slack"),
                "{error}"
            );
        }
    }

    #[test]
    fn g0_the_live_tracker_projects_and_the_artifact_validates() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let artifact = load(root).expect("the tracked artifact must exist and validate");
        let issues = read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES).expect("tracker");
        let live = project(&issues).expect("the live tracker must project");
        assert_eq!(artifact.slack_hours.len(), artifact.open);
        assert!(live.open > 0);
        // The reality check's question, answered as a number: each spine
        // bead's slack is either recorded, or the bead is closed/absent.
        for id in SPINE_BEADS {
            if artifact.slack_hours.contains_key(id) {
                assert!(artifact.spine_positions.contains_key(id), "{id}");
            }
        }
        let (violations, notes) = check(root);
        assert!(violations.is_empty(), "{violations:?}");
        assert!(
            notes
                .iter()
                .any(|note| matches!(note.verdict, "current" | "stale-snapshot")),
            "currency must be rendered explicitly: {notes:?}"
        );
    }
}
