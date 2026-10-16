//! Checked EXTREAL critical-path projection (bead
//! `frankensim-extreal-program-f85xj.16.2`).
//!
//! Beads remains authoritative for work status and dependency edges. The root
//! `vertical-capability-graph.json` is deliberately a projection: it binds the
//! maturity registry's product capabilities to real issue ids, names owners
//! for the four load-bearing integration seams, and retains one robot-triage
//! receipt. This check rejects drift between those authorities without
//! pretending that membership or closure proves scientific maturity.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub const GRAPH_FILE: &str = "vertical-capability-graph.json";
const REGISTRY_FILE: &str = "capability-maturity.json";
const BEADS_FILE: &str = ".beads/issues.jsonl";
const GRAPH_SCHEMA: &str = "frankensim-vertical-capability-graph-v1";
const PROGRAM_PREFIX: &str = "frankensim-extreal-program-f85xj";
const PROGRAM_LABEL: &str = "extreal";
pub const CHECK: &str = "vertical-critical-path";
const MAX_BEAD_ROW_BYTES: usize = 1024 * 1024;
const EXPECTED_SEAMS: [&str; 4] = [
    "cli-session",
    "corpus-scorecard",
    "physics-evidence",
    "schema-scenario",
];
const EXPECTED_SHAPE: [&str; 4] = ["data_hash", "status", "summary", "tracks"];

pub struct CriticalPathReport {
    pub violations: Vec<Violation>,
    pub decisions: Vec<PolicyNote>,
}

#[derive(Clone)]
struct Issue {
    status: String,
    labels: BTreeSet<String>,
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

fn number(value: &JsonValue) -> Option<&str> {
    match value {
        JsonValue::Number(value) => Some(value),
        _ => None,
    }
}

fn field<'a>(map: &'a BTreeMap<String, JsonValue>, key: &str) -> Option<&'a JsonValue> {
    map.get(key)
}

fn string_set(value: Option<&JsonValue>) -> Option<BTreeSet<String>> {
    let items = arr(value?)?;
    let strings: Option<BTreeSet<_>> = items
        .iter()
        .map(|item| text(item).map(str::to_string))
        .collect();
    strings
}

fn parse_root(source: &str, entity: &str, violations: &mut Vec<Violation>) -> Option<JsonValue> {
    match JsonParser::new(source).finish() {
        Ok(value) if obj(&value).is_some() => Some(value),
        Ok(_) => {
            violations.push(violation(entity, "document root must be a JSON object"));
            None
        }
        Err(error) => {
            violations.push(violation(entity, format!("invalid JSON: {error}")));
            None
        }
    }
}

fn parse_registry(source: &str, violations: &mut Vec<Violation>) -> BTreeMap<String, String> {
    let Some(parsed) = parse_root(source, REGISTRY_FILE, violations) else {
        return BTreeMap::new();
    };
    let root = obj(&parsed).expect("parse_root returns an object");
    let Some(capabilities) = field(root, "capabilities").and_then(arr) else {
        violations.push(violation(
            REGISTRY_FILE,
            "missing array field `capabilities`",
        ));
        return BTreeMap::new();
    };

    let mut result = BTreeMap::new();
    for (index, capability) in capabilities.iter().enumerate() {
        let Some(capability) = obj(capability) else {
            violations.push(violation(
                REGISTRY_FILE,
                format!("capabilities[{index}] must be an object"),
            ));
            continue;
        };
        let (Some(id), Some(level)) = (
            field(capability, "id").and_then(text),
            field(capability, "level").and_then(text),
        ) else {
            violations.push(violation(
                REGISTRY_FILE,
                format!("capabilities[{index}] needs string `id` and `level`"),
            ));
            continue;
        };
        if result.insert(id.to_string(), level.to_string()).is_some() {
            violations.push(violation(
                REGISTRY_FILE,
                format!("capability {id:?} occurs more than once"),
            ));
        }
    }
    result
}

fn parse_issues(source: &str, violations: &mut Vec<Violation>) -> BTreeMap<String, Issue> {
    let mut issues = BTreeMap::new();
    for (line_index, line) in source.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let entity = format!("{BEADS_FILE}:{}", line_index + 1);
        if line.len() > MAX_BEAD_ROW_BYTES {
            violations.push(violation(
                &entity,
                format!(
                    "issue row is {} bytes; bounded checker limit is {MAX_BEAD_ROW_BYTES}",
                    line.len()
                ),
            ));
            continue;
        }
        let parsed = match JsonParser::with_string_limit(line, MAX_BEAD_ROW_BYTES).finish() {
            Ok(value) => value,
            Err(error) => {
                violations.push(violation(&entity, format!("invalid JSONL row: {error}")));
                continue;
            }
        };
        let Some(row) = obj(&parsed) else {
            violations.push(violation(&entity, "issue row must be an object"));
            continue;
        };
        let (Some(id), Some(status)) = (
            field(row, "id").and_then(text),
            field(row, "status").and_then(text),
        ) else {
            violations.push(violation(
                &entity,
                "issue row needs string `id` and `status`",
            ));
            continue;
        };
        let labels = string_set(field(row, "labels")).unwrap_or_default();
        if issues
            .insert(
                id.to_string(),
                Issue {
                    status: status.to_string(),
                    labels,
                },
            )
            .is_some()
        {
            violations.push(violation(BEADS_FILE, format!("duplicate issue id {id:?}")));
        }
    }
    issues
}

fn check_issue_ref(
    id: &str,
    context: &str,
    issues: &BTreeMap<String, Issue>,
    violations: &mut Vec<Violation>,
) -> bool {
    if issues.contains_key(id) {
        true
    } else {
        violations.push(violation(
            context,
            format!("references missing Bead {id:?}"),
        ));
        false
    }
}

fn check_bindings(
    root: &BTreeMap<String, JsonValue>,
    capabilities: &BTreeMap<String, String>,
    issues: &BTreeMap<String, Issue>,
    violations: &mut Vec<Violation>,
) -> usize {
    let Some(bindings) = field(root, "capability_bindings").and_then(arr) else {
        violations.push(violation(GRAPH_FILE, "missing array `capability_bindings`"));
        return 0;
    };
    let mut bound = BTreeSet::new();
    for (index, binding) in bindings.iter().enumerate() {
        let entity = format!("capability_bindings[{index}]");
        let Some(binding) = obj(binding) else {
            violations.push(violation(&entity, "binding must be an object"));
            continue;
        };
        let Some(id) = field(binding, "capability_id").and_then(text) else {
            violations.push(violation(&entity, "missing string `capability_id`"));
            continue;
        };
        if !bound.insert(id.to_string()) {
            violations.push(violation(&entity, format!("duplicate capability {id:?}")));
        }
        let Some(level) = capabilities.get(id) else {
            violations.push(violation(
                &entity,
                format!("capability {id:?} is absent from {REGISTRY_FILE}"),
            ));
            continue;
        };
        let Some(beads) = string_set(field(binding, "implementing_beads")) else {
            violations.push(violation(
                &entity,
                "`implementing_beads` must be a string array",
            ));
            continue;
        };
        if beads.is_empty() {
            violations.push(violation(&entity, "`implementing_beads` must not be empty"));
        }
        let mut has_closed = false;
        for bead in beads {
            if check_issue_ref(&bead, &entity, issues, violations)
                && issues
                    .get(&bead)
                    .is_some_and(|issue| issue.status == "closed")
            {
                has_closed = true;
            }
        }
        if level != "L1" && !has_closed {
            violations.push(violation(
                &entity,
                format!(
                    "registry level {level} requires at least one closed implementing Bead; open work alone cannot justify delivered maturity"
                ),
            ));
        }
    }

    let registered: BTreeSet<_> = capabilities.keys().cloned().collect();
    for missing in registered.difference(&bound) {
        violations.push(violation(
            GRAPH_FILE,
            format!("registry capability {missing:?} has no Bead binding"),
        ));
    }
    for extra in bound.difference(&registered) {
        violations.push(violation(
            GRAPH_FILE,
            format!("binding {extra:?} has no registry capability"),
        ));
    }
    bindings.len()
}

fn check_seams(
    root: &BTreeMap<String, JsonValue>,
    issues: &BTreeMap<String, Issue>,
    violations: &mut Vec<Violation>,
) -> usize {
    let Some(seams) = field(root, "seam_owners").and_then(arr) else {
        violations.push(violation(GRAPH_FILE, "missing array `seam_owners`"));
        return 0;
    };
    let mut actual = BTreeSet::new();
    for (index, seam) in seams.iter().enumerate() {
        let entity = format!("seam_owners[{index}]");
        let Some(seam) = obj(seam) else {
            violations.push(violation(&entity, "seam must be an object"));
            continue;
        };
        let Some(name) = field(seam, "seam").and_then(text) else {
            violations.push(violation(&entity, "missing string `seam`"));
            continue;
        };
        if !actual.insert(name.to_string()) {
            violations.push(violation(&entity, format!("duplicate seam {name:?}")));
        }
        if field(seam, "owner")
            .and_then(text)
            .is_none_or(str::is_empty)
        {
            violations.push(violation(&entity, "seam must have a named `owner`"));
        }
        if field(seam, "responsibility")
            .and_then(text)
            .is_none_or(str::is_empty)
        {
            violations.push(violation(
                &entity,
                "seam must state its integration `responsibility`",
            ));
        }
        let Some(beads) = string_set(field(seam, "beads")) else {
            violations.push(violation(&entity, "`beads` must be a string array"));
            continue;
        };
        if beads.is_empty() {
            violations.push(violation(&entity, "seam `beads` must not be empty"));
        }
        for bead in beads {
            check_issue_ref(&bead, &entity, issues, violations);
        }
    }

    let expected: BTreeSet<_> = EXPECTED_SEAMS
        .iter()
        .map(|value| (*value).to_string())
        .collect();
    if actual != expected {
        violations.push(violation(
            GRAPH_FILE,
            format!("seam set must be exactly {expected:?}, found {actual:?}"),
        ));
    }
    seams.len()
}

fn check_receipt(root: &BTreeMap<String, JsonValue>, violations: &mut Vec<Violation>) {
    let Some(receipt) = field(root, "retained_triage_receipt").and_then(obj) else {
        violations.push(violation(
            GRAPH_FILE,
            "missing object `retained_triage_receipt`",
        ));
        return;
    };
    if field(receipt, "command").and_then(text) != Some("bv --robot-plan --label extreal") {
        violations.push(violation(
            "retained_triage_receipt",
            "command must be exactly `bv --robot-plan --label extreal`",
        ));
    }
    if field(receipt, "scope_label").and_then(text) != Some(PROGRAM_LABEL) {
        violations.push(violation(
            "retained_triage_receipt",
            "scope_label must be `extreal`",
        ));
    }
    let generated = field(receipt, "generated_at").and_then(text).unwrap_or("");
    if generated.len() != 20
        || !generated.ends_with('Z')
        || generated.as_bytes().get(10) != Some(&b'T')
    {
        violations.push(violation(
            "retained_triage_receipt",
            "generated_at must have UTC `YYYY-MM-DDTHH:MM:SSZ` shape",
        ));
    }
    let hash = field(receipt, "data_hash").and_then(text).unwrap_or("");
    if hash.len() != 16
        || !hash
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        violations.push(violation(
            "retained_triage_receipt",
            "data_hash must be 16 lowercase hexadecimal digits",
        ));
    }
    for key in [
        "issue_count",
        "open",
        "closed",
        "blocked_count",
        "health",
        "total_actionable",
        "total_blocked",
    ] {
        let valid = field(receipt, key)
            .and_then(number)
            .is_some_and(|value| value.bytes().all(|byte| byte.is_ascii_digit()));
        if !valid {
            violations.push(violation(
                "retained_triage_receipt",
                format!("{key} must be a non-negative integer"),
            ));
        }
    }
    let shape = string_set(field(receipt, "documented_shape")).unwrap_or_default();
    let expected: BTreeSet<_> = EXPECTED_SHAPE
        .iter()
        .map(|value| (*value).to_string())
        .collect();
    if shape != expected {
        violations.push(violation(
            "retained_triage_receipt",
            format!("documented_shape must be exactly {expected:?}"),
        ));
    }
}

fn check_sources(graph: &str, registry: &str, beads: &str) -> CriticalPathReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();
    let capabilities = parse_registry(registry, &mut violations);
    let issues = parse_issues(beads, &mut violations);
    let Some(parsed) = parse_root(graph, GRAPH_FILE, &mut violations) else {
        return CriticalPathReport {
            violations,
            decisions,
        };
    };
    let root = obj(&parsed).expect("parse_root returns an object");

    for (key, expected) in [
        ("schema", GRAPH_SCHEMA),
        ("graph_authority", BEADS_FILE),
        ("program_label", PROGRAM_LABEL),
        ("program_id_prefix", PROGRAM_PREFIX),
    ] {
        if field(root, key).and_then(text) != Some(expected) {
            violations.push(violation(GRAPH_FILE, format!("{key} must be {expected:?}")));
        }
    }
    for key in ["purpose", "no_claim"] {
        if field(root, key).and_then(text).is_none_or(str::is_empty) {
            violations.push(violation(GRAPH_FILE, format!("{key} must be stated")));
        }
    }

    let linked_gates = string_set(field(root, "linked_gates")).unwrap_or_default();
    for gate in ["frankensim-ty23", "frankensim-v6dn"] {
        if !linked_gates.contains(gate) {
            violations.push(violation(
                GRAPH_FILE,
                format!("missing linked gate {gate:?}"),
            ));
        }
        check_issue_ref(gate, "linked_gates", &issues, &mut violations);
    }
    let policy = field(root, "off_path_policy_bead")
        .and_then(text)
        .unwrap_or("");
    if policy != "frankensim-extreal-program-f85xj.16.3" {
        violations.push(violation(
            GRAPH_FILE,
            "off_path_policy_bead must be f85xj.16.3",
        ));
    }
    check_issue_ref(policy, "off_path_policy_bead", &issues, &mut violations);

    let native: Vec<_> = issues
        .iter()
        .filter(|(id, _)| id.starts_with(PROGRAM_PREFIX))
        .collect();
    for (id, issue) in &native {
        if !issue.labels.contains(PROGRAM_LABEL) {
            violations.push(violation(
                id,
                format!("native EXTREAL Bead lacks required label {PROGRAM_LABEL:?}"),
            ));
        }
    }
    let binding_count = check_bindings(root, &capabilities, &issues, &mut violations);
    let seam_count = check_seams(root, &issues, &mut violations);
    check_receipt(root, &mut violations);

    decisions.push(note(
        "<repo>",
        "inventory",
        format!(
            "{} native EXTREAL Beads, {} capability bindings, and {} named seams checked; dependency/status authority remains {BEADS_FILE}",
            native.len(), binding_count, seam_count
        ),
    ));
    CriticalPathReport {
        violations,
        decisions,
    }
}

fn check_documentation(root: &Path, violations: &mut Vec<Violation>) {
    let path = root.join("docs/CONVENTIONS.md");
    let Ok(source) = std::fs::read_to_string(&path) else {
        violations.push(violation("docs/CONVENTIONS.md", "document is unreadable"));
        return;
    };
    for required in [
        "bv --robot-plan --label extreal",
        "bv --robot-triage --label extreal",
        "vertical-capability-graph.json",
        "schema-scenario",
        "physics-evidence",
        "corpus-scorecard",
        "cli-session",
        "f85xj.16.3",
    ] {
        if !source.contains(required) {
            violations.push(violation(
                "docs/CONVENTIONS.md",
                format!("critical-path doctrine is missing {required:?}"),
            ));
        }
    }
}

pub fn check_critical_path(root: &Path) -> CriticalPathReport {
    let mut read_violations = Vec::new();
    let read = |relative: &str, violations: &mut Vec<Violation>| {
        std::fs::read_to_string(root.join(relative)).map_err(|error| {
            violations.push(violation(relative, format!("file is unreadable: {error}")));
        })
    };
    let graph = read(GRAPH_FILE, &mut read_violations);
    let registry = read(REGISTRY_FILE, &mut read_violations);
    let beads = read(BEADS_FILE, &mut read_violations);
    let (Ok(graph), Ok(registry), Ok(beads)) = (graph, registry, beads) else {
        return CriticalPathReport {
            violations: read_violations,
            decisions: Vec::new(),
        };
    };
    let mut report = check_sources(&graph, &registry, &beads);
    report.violations.extend(read_violations);
    check_documentation(root, &mut report.violations);
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeded_cross_reference_and_label_faults_fail_closed() {
        let registry = r#"{"capabilities":[{"id":"cap.one","level":"L2"}]}"#;
        let beads = [
            r#"{"id":"frankensim-extreal-program-f85xj.1","status":"closed","labels":["extreal"]}"#,
            r#"{"id":"frankensim-ty23","status":"open","labels":["extreal"]}"#,
            r#"{"id":"frankensim-v6dn","status":"open","labels":["extreal"]}"#,
            r#"{"id":"frankensim-extreal-program-f85xj.16.3","status":"open","labels":["extreal"]}"#,
        ]
        .join("\n");
        let seams = EXPECTED_SEAMS
            .iter()
            .map(|name| format!(r#"{{"seam":"{name}","owner":"o","responsibility":"r","beads":["frankensim-extreal-program-f85xj.1"]}}"#))
            .collect::<Vec<_>>()
            .join(",");
        let graph = format!(
            r#"{{"schema":"{GRAPH_SCHEMA}","graph_authority":"{BEADS_FILE}","program_label":"{PROGRAM_LABEL}","program_id_prefix":"{PROGRAM_PREFIX}","purpose":"p","no_claim":"n","linked_gates":["frankensim-ty23","frankensim-v6dn"],"off_path_policy_bead":"frankensim-extreal-program-f85xj.16.3","capability_bindings":[{{"capability_id":"cap.one","implementing_beads":["frankensim-extreal-program-f85xj.1"]}}],"seam_owners":[{seams}],"retained_triage_receipt":{{"command":"bv --robot-plan --label extreal","generated_at":"2026-07-22T21:19:43Z","data_hash":"0123456789abcdef","scope_label":"extreal","issue_count":4,"open":3,"closed":1,"blocked_count":0,"health":100,"total_actionable":3,"total_blocked":0,"documented_shape":["summary","tracks","status","data_hash"]}}}}"#
        );
        assert!(
            check_sources(&graph, registry, &beads)
                .violations
                .is_empty()
        );

        let missing_label = beads.replacen(r#"["extreal"]"#, "[]", 1);
        let report = check_sources(&graph, registry, &missing_label);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("lacks required label")),
            "seeded label fault must be detected: {:?}",
            report.violations
        );

        let missing_bead = graph.replacen(
            r#""implementing_beads":["frankensim-extreal-program-f85xj.1"]"#,
            r#""implementing_beads":["ghost"]"#,
            1,
        );
        let report = check_sources(&missing_bead, registry, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("missing Bead")),
            "seeded bead fault must be detected: {:?}",
            report.violations
        );

        let unknown_capability = graph.replacen("cap.one", "cap.ghost", 1);
        let report = check_sources(&unknown_capability, registry, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("absent from capability-maturity.json")),
            "seeded capability fault must be detected: {:?}",
            report.violations
        );
    }

    #[test]
    fn the_live_vertical_projection_is_clean() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let report = check_critical_path(root);
        assert!(
            report.violations.is_empty(),
            "live vertical-capability-graph.json must be clean: {:?}",
            report.violations
        );
        assert!(
            report
                .decisions
                .iter()
                .any(|item| item.verdict == "inventory")
        );
    }
}
