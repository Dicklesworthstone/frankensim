//! Capability maturity registry check (bead
//! `frankensim-extreal-program-f85xj.16.1`).
//!
//! "Is this capability experimental, verified, integrated, validated, or
//! supported?" had no queryable answer: status lived in README prose, CONTRACT
//! no-claim sections, and folklore. `capability-maturity.json` records the
//! answer and this check keeps the record honest:
//!
//! 1. the registry parses, matches its schema, and every entry is well formed;
//! 2. every evidence ref RESOLVES — the registry may not cite a test, contract,
//!    lane, or document that does not exist;
//! 3. each level's own evidence bar is met (see `docs/MATURITY_LEVELS.md`);
//! 4. PROMOTIONS since the last committed registry are surfaced as policy
//!    notes, because bead `.2.3`'s claim-integrity gate consumes exactly that
//!    signal;
//! 5. DEMOTIONS are always allowed and merely logged.
//!
//! The asymmetry in (4)/(5) is deliberate. Lowering a claim is how the registry
//! stays truthful, so it must never be procedurally harder than raising one; a
//! system that resists demotion accumulates false claims by construction.
//!
//! What this check does NOT do: judge whether a cited test actually exercises
//! the capability, or whether a level is deserved. It proves the paperwork is
//! present and internally consistent. Claiming more would make this check
//! itself a claim-integrity defect.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

pub const REGISTRY_FILE: &str = "capability-maturity.json";
const REGISTRY_SCHEMA: &str = "frankensim-capability-maturity-v1";
const SPINE_E2E_RECEIPT_SCHEMA: &str = "frankensim-spine-e2e-receipt-v1";
const CHECK: &str = "capability-maturity";
const LEVELS: [&str; 5] = ["L1", "L2", "L3", "L4", "L5"];
const README_MATRIX_BEGIN: &str = "<!-- BEGIN GENERATED FRANKENSIM CAPABILITY MATRIX -->";
const README_MATRIX_END: &str = "<!-- END GENERATED FRANKENSIM CAPABILITY MATRIX -->";

/// Evidence kinds and whether the check can resolve them against the tree.
/// `corpus` is recorded but unresolvable until the V&V corpus registry (e04)
/// exists — an honest gap, and the reason nothing is L4 today.
const RESOLVABLE_KINDS: [&str; 5] = ["test", "contract", "lane", "doc", "receipt"];
const RECORDED_ONLY_KINDS: [&str; 1] = ["corpus"];

pub struct MaturityReport {
    pub violations: Vec<Violation>,
    pub decisions: Vec<PolicyNote>,
}

fn violation(entity: &str, detail: String) -> Violation {
    Violation {
        check: CHECK,
        crate_name: entity.to_string(),
        detail,
    }
}

fn note(entity: &str, verdict: &'static str, detail: String) -> PolicyNote {
    PolicyNote {
        check: CHECK,
        crate_name: entity.to_string(),
        verdict,
        detail,
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
        JsonValue::String(s) => Some(s.as_str()),
        _ => None,
    }
}

fn field<'a>(map: &'a BTreeMap<String, JsonValue>, key: &str) -> Option<&'a JsonValue> {
    map.get(key)
}

/// `YYYY-MM-DD`, validated structurally (no calendar arithmetic: the check
/// stays deterministic and wall-clock-free so `check-all` output does not
/// change with the date).
fn is_iso_date(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() == 10
        && bytes[4] == b'-'
        && bytes[7] == b'-'
        && bytes
            .iter()
            .enumerate()
            .all(|(index, byte)| matches!(index, 4 | 7) || byte.is_ascii_digit())
}

fn level_index(level: &str) -> Option<usize> {
    LEVELS.iter().position(|candidate| *candidate == level)
}

/// Ordinal of a level name, for callers comparing two levels. `None` for an
/// unrecognized name — the caller must treat that as "cannot compare", never
/// as "equal".
pub fn level_rank(level: &str) -> Option<usize> {
    level_index(level)
}

/// The registry's levels now and as last committed, plus each capability's
/// crate scope. The claim-integrity promotion gate (bead `.2.3`) consumes this
/// to decide which capabilities are being promoted and what a defect must
/// overlap to block one.
pub struct CapabilityLevels {
    pub current: BTreeMap<String, String>,
    pub committed: BTreeMap<String, String>,
    pub crates: BTreeMap<String, BTreeSet<String>>,
}

/// Pull `(id -> level)` and `(id -> crate scopes)` out of a registry document.
/// Structural defects are the business of `check_maturity`; this extraction
/// simply skips what it cannot read, because the gate must not double-report
/// the registry's own validity problems.
fn levels_and_scopes(
    source: &str,
) -> (BTreeMap<String, String>, BTreeMap<String, BTreeSet<String>>) {
    let mut levels = BTreeMap::new();
    let mut scopes = BTreeMap::new();
    let Ok(parsed) = JsonParser::new(source).finish() else {
        return (levels, scopes);
    };
    let Some(items) = obj(&parsed)
        .and_then(|map| field(map, "capabilities"))
        .and_then(arr)
    else {
        return (levels, scopes);
    };
    for item in items {
        let Some(map) = obj(item) else { continue };
        let (Some(id), Some(level)) = (
            field(map, "id").and_then(text),
            field(map, "level").and_then(text),
        ) else {
            continue;
        };
        levels.insert(id.to_string(), level.to_string());
        let crates = field(map, "crates")
            .and_then(arr)
            .map(|items| items.iter().filter_map(text).map(str::to_string).collect())
            .unwrap_or_default();
        scopes.insert(id.to_string(), crates);
    }
    (levels, scopes)
}

/// Read the working registry and its last committed state.
///
/// A missing committed predecessor is not an error: the registry is new, so
/// nothing in it is a promotion. An unreadable working registry IS an error,
/// because a gate that cannot see the levels must refuse rather than conclude
/// that nothing is being promoted.
pub fn capability_levels(root: &Path) -> Result<CapabilityLevels, String> {
    let source = std::fs::read_to_string(root.join(REGISTRY_FILE)).map_err(|error| {
        format!(
            "{REGISTRY_FILE} is unreadable ({error}); the promotion gate cannot conclude that \
             nothing is being promoted from a registry it could not read"
        )
    })?;
    let (current, crates) = levels_and_scopes(&source);

    let output = std::process::Command::new("git")
        .args(["show", &format!("HEAD:{REGISTRY_FILE}")])
        .current_dir(root)
        .output();
    let committed = match output {
        Ok(output) if output.status.success() => String::from_utf8(output.stdout)
            .map(|text| levels_and_scopes(&text).0)
            .unwrap_or_default(),
        _ => BTreeMap::new(),
    };

    Ok(CapabilityLevels {
        current,
        committed,
        crates,
    })
}

/// One registry entry, reduced to what the check reasons about.
struct Entry {
    id: String,
    title: String,
    level: String,
    crates: Vec<String>,
    notes: String,
    kinds: BTreeSet<String>,
    lanes: Vec<String>,
    receipts: Vec<ReceiptEvidence>,
}

/// The receipt evidence an L3+ entry must carry. `stage` is deliberate: a
/// successful adjacent stage cannot silently stand in for this capability.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ReceiptEvidence {
    path: String,
    stage: String,
}

/// Parse the registry into entries, pushing a violation for every structural
/// defect. Returns entries for whatever parsed cleanly so one bad row does not
/// hide the rest.
fn parse_registry(source: &str, entity: &str, violations: &mut Vec<Violation>) -> Vec<Entry> {
    let parsed = match JsonParser::new(source).finish() {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(
                entity,
                format!("{REGISTRY_FILE} is not valid JSON: {error}"),
            ));
            return Vec::new();
        }
    };
    let Some(root) = obj(&parsed) else {
        violations.push(violation(
            entity,
            format!("{REGISTRY_FILE} is not a JSON object"),
        ));
        return Vec::new();
    };
    match field(root, "schema").and_then(text) {
        Some(REGISTRY_SCHEMA) => {}
        Some(other) => violations.push(violation(
            entity,
            format!("{REGISTRY_FILE} declares schema {other:?}, expected {REGISTRY_SCHEMA:?}"),
        )),
        None => violations.push(violation(
            entity,
            format!("{REGISTRY_FILE} has no string \"schema\" field"),
        )),
    }
    let Some(items) = field(root, "capabilities").and_then(arr) else {
        violations.push(violation(
            entity,
            format!("{REGISTRY_FILE} has no \"capabilities\" array"),
        ));
        return Vec::new();
    };

    let mut entries = Vec::new();
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for (index, item) in items.iter().enumerate() {
        let Some(map) = obj(item) else {
            violations.push(violation(
                entity,
                format!("capability #{index} is not a JSON object"),
            ));
            continue;
        };
        let Some(id) = field(map, "id").and_then(text).filter(|id| !id.is_empty()) else {
            violations.push(violation(
                entity,
                format!("capability #{index} has no non-empty string \"id\""),
            ));
            continue;
        };
        if !seen.insert(id.to_string()) {
            violations.push(violation(
                id,
                format!("capability id {id:?} appears more than once; ids are the registry's key"),
            ));
            continue;
        }
        for required in ["title", "owner", "level", "last_review"] {
            if field(map, required)
                .and_then(text)
                .is_none_or(str::is_empty)
            {
                violations.push(violation(
                    id,
                    format!("capability {id:?} has no non-empty string {required:?}"),
                ));
            }
        }
        let level = field(map, "level").and_then(text).unwrap_or_default();
        if !level.is_empty() && level_index(level).is_none() {
            violations.push(violation(
                id,
                format!("capability {id:?} declares level {level:?}; expected one of {LEVELS:?}"),
            ));
        }
        if let Some(review) = field(map, "last_review").and_then(text)
            && !is_iso_date(review)
        {
            violations.push(violation(
                id,
                format!("capability {id:?} last_review {review:?} is not YYYY-MM-DD"),
            ));
        }
        let mut crates = Vec::new();
        match field(map, "crates").and_then(arr) {
            Some(items) if !items.is_empty() => {
                for (crate_index, item) in items.iter().enumerate() {
                    match text(item).filter(|name| !name.is_empty()) {
                        Some(name) => crates.push(name.to_string()),
                        None => violations.push(violation(
                            id,
                            format!(
                                "capability {id:?} crate scope #{crate_index} is not a non-empty string"
                            ),
                        )),
                    }
                }
            }
            _ => violations.push(violation(
                id,
                format!(
                    "capability {id:?} has no non-empty \"crates\" scope array; scope is what the \
                     claim-integrity promotion gate matches on, and an unscoped capability would \
                     have to be treated as global"
                ),
            )),
        }
        let (kinds, lanes, receipts) = collect_evidence(map, id, violations);
        entries.push(Entry {
            id: id.to_string(),
            title: field(map, "title")
                .and_then(text)
                .unwrap_or_default()
                .to_string(),
            level: level.to_string(),
            crates,
            notes: field(map, "notes")
                .and_then(text)
                .unwrap_or_default()
                .to_string(),
            kinds,
            lanes,
            receipts,
        });
    }
    entries
}

fn markdown_cell(value: &str) -> String {
    value
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .replace('|', "\\|")
}

fn render_readme_matrix(entries: &[Entry]) -> String {
    let mut ordered: Vec<&Entry> = entries.iter().collect();
    ordered.sort_by(|left, right| left.id.cmp(&right.id));
    let mut output = format!(
        "{README_MATRIX_BEGIN}\n\
| Capability | Registry level | Crate scope | Registry boundary |\n\
|------------|----------------|-------------|-------------------|\n"
    );
    for entry in ordered {
        let crates = entry
            .crates
            .iter()
            .map(|name| format!("`{}`", markdown_cell(name)))
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(
            output,
            "| `{}` — {} | {} | {crates} | {} |",
            markdown_cell(&entry.id),
            markdown_cell(&entry.title),
            entry.level,
            markdown_cell(&entry.notes),
        )
        .ok();
    }
    output.push_str(README_MATRIX_END);
    output
}

fn readme_matrix_block(source: &str) -> Result<&str, String> {
    let starts: Vec<usize> = source
        .match_indices(README_MATRIX_BEGIN)
        .map(|(index, _)| index)
        .collect();
    let ends: Vec<usize> = source
        .match_indices(README_MATRIX_END)
        .map(|(index, _)| index)
        .collect();
    if starts.len() != 1 || ends.len() != 1 {
        return Err(format!(
            "expected exactly one generated capability-matrix marker pair, found {} starts and {} ends",
            starts.len(),
            ends.len()
        ));
    }
    let start = starts[0];
    let finish = ends[0]
        .checked_add(README_MATRIX_END.len())
        .ok_or_else(|| "capability-matrix end offset overflow".to_string())?;
    if ends[0] <= start {
        return Err("capability-matrix end marker precedes its start marker".to_string());
    }
    Ok(&source[start..finish])
}

fn check_readme_matrix_text(readme: &str, entries: &[Entry]) -> Vec<Violation> {
    let expected = render_readme_matrix(entries);
    match readme_matrix_block(readme) {
        Ok(actual) if actual == expected => Vec::new(),
        Ok(_) => vec![violation(
            "README.md",
            format!(
                "README generated capability matrix is stale; replace it with the exact registry projection:\n{expected}"
            ),
        )],
        Err(error) => vec![violation(
            "README.md",
            format!("README generated capability matrix is malformed: {error}"),
        )],
    }
}

fn check_readme_summary_counts(readme: &str, entries: &[Entry]) -> Vec<Violation> {
    let mut violations = Vec::new();
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for entry in entries {
        *counts.entry(entry.level.as_str()).or_default() += 1;
    }
    for level in LEVELS {
        let rows: Vec<&str> = readme
            .lines()
            .filter(|line| line.starts_with(&format!("| {level} |")))
            .collect();
        if rows.len() != 1 {
            violations.push(violation(
                "README.md",
                format!(
                    "README maturity summary must contain exactly one {level} row, found {}",
                    rows.len()
                ),
            ));
            continue;
        }
        let cells: Vec<&str> = rows[0].split('|').map(str::trim).collect();
        let claimed = cells.get(3).and_then(|cell| cell.parse::<usize>().ok());
        let actual = counts.get(level).copied().unwrap_or(0);
        if claimed != Some(actual) {
            violations.push(violation(
                "README.md",
                format!(
                    "README maturity summary claims {level}={claimed:?}, but {REGISTRY_FILE} has {actual}"
                ),
            ));
        }
    }
    let total_lines: Vec<&str> = readme
        .lines()
        .filter(|line| line.contains(" product-meaningful capabilities"))
        .collect();
    if total_lines.len() != 1 {
        violations.push(violation(
            "README.md",
            format!(
                "README capability summary must contain exactly one product-meaningful-capabilities total, found {}",
                total_lines.len()
            ),
        ));
    } else {
        let line = total_lines[0];
        let marker = " product-meaningful capabilities";
        let Some(position) = line.find(marker) else {
            return violations;
        };
        let digits: String = line[..position]
            .chars()
            .rev()
            .take_while(char::is_ascii_digit)
            .collect::<String>()
            .chars()
            .rev()
            .collect();
        if digits.parse::<usize>().ok() != Some(entries.len()) {
            violations.push(violation(
                "README.md",
                format!(
                    "README capability total {digits:?} does not match the {} registry entries",
                    entries.len()
                ),
            ));
        }
    }
    violations
}

fn check_readme_projection_entries(root: &Path, entries: &[Entry]) -> MaturityReport {
    let readme = match std::fs::read_to_string(root.join("README.md")) {
        Ok(readme) => readme,
        Err(error) => {
            return MaturityReport {
                violations: vec![violation(
                    "README.md",
                    format!("cannot read README.md for capability-matrix drift check: {error}"),
                )],
                decisions: Vec::new(),
            };
        }
    };
    let mut violations = check_readme_matrix_text(&readme, entries);
    violations.extend(check_readme_summary_counts(&readme, entries));
    let decisions = if violations.is_empty() {
        entries
            .iter()
            .map(|entry| {
                note(
                    &entry.id,
                    "verified",
                    format!(
                        "README capability projection matches {REGISTRY_FILE}: level={} crates={}",
                        entry.level,
                        entry.crates.join(",")
                    ),
                )
            })
            .collect()
    } else {
        Vec::new()
    };
    MaturityReport {
        violations,
        decisions,
    }
}

pub fn check_readme_projection(root: &Path) -> MaturityReport {
    let mut violations = Vec::new();
    let path = root.join(REGISTRY_FILE);
    let source = match std::fs::read_to_string(&path) {
        Ok(source) => source,
        Err(error) => {
            return MaturityReport {
                violations: vec![violation(
                    REGISTRY_FILE,
                    format!("cannot read {REGISTRY_FILE}: {error}"),
                )],
                decisions: Vec::new(),
            };
        }
    };
    let entries = parse_registry(&source, REGISTRY_FILE, &mut violations);
    if !violations.is_empty() {
        return MaturityReport {
            violations,
            decisions: Vec::new(),
        };
    }
    check_readme_projection_entries(root, &entries)
}

/// Validate the evidence array and return the set of kinds present. Refs are
/// resolved against the tree by `resolve_refs`, which needs the repo root.
fn collect_evidence(
    map: &BTreeMap<String, JsonValue>,
    id: &str,
    violations: &mut Vec<Violation>,
) -> (BTreeSet<String>, Vec<String>, Vec<ReceiptEvidence>) {
    let mut kinds = BTreeSet::new();
    let mut lanes = Vec::new();
    let mut receipts = Vec::new();
    let Some(items) = field(map, "evidence").and_then(arr) else {
        violations.push(violation(
            id,
            format!("capability {id:?} has no \"evidence\" array"),
        ));
        return (kinds, lanes, receipts);
    };
    for (index, item) in items.iter().enumerate() {
        let Some(entry) = obj(item) else {
            violations.push(violation(
                id,
                format!("capability {id:?} evidence #{index} is not an object"),
            ));
            continue;
        };
        let kind = field(entry, "kind").and_then(text).unwrap_or_default();
        let reference = field(entry, "ref").and_then(text).unwrap_or_default();
        if kind.is_empty() || reference.is_empty() {
            violations.push(violation(
                id,
                format!("capability {id:?} evidence #{index} needs non-empty \"kind\" and \"ref\""),
            ));
            continue;
        }
        if !RESOLVABLE_KINDS.contains(&kind) && !RECORDED_ONLY_KINDS.contains(&kind) {
            violations.push(violation(
                id,
                format!(
                    "capability {id:?} evidence #{index} has unknown kind {kind:?}; expected one \
                     of {RESOLVABLE_KINDS:?} or {RECORDED_ONLY_KINDS:?}"
                ),
            ));
            continue;
        }
        kinds.insert(kind.to_string());
        if kind == "lane" {
            lanes.push(reference.to_string());
        }
        if kind == "receipt" {
            let Some(stage) = field(entry, "stage")
                .and_then(text)
                .filter(|stage| !stage.is_empty())
            else {
                violations.push(violation(
                    id,
                    format!(
                        "capability {id:?} receipt evidence {reference:?} needs a non-empty \"stage\"; \
                         L3 proof is capability-stage-specific"
                    ),
                ));
                continue;
            };
            receipts.push(ReceiptEvidence {
                path: reference.to_string(),
                stage: stage.to_string(),
            });
        }
    }
    (kinds, lanes, receipts)
}

/// Resolve every evidence ref against the tree. A registry that cites evidence
/// which is not there is exactly the defect class this program exists to stop.
fn resolve_refs(root: &Path, source: &str, violations: &mut Vec<Violation>) {
    let Ok(parsed) = JsonParser::new(source).finish() else {
        return;
    };
    let Some(items) = obj(&parsed)
        .and_then(|map| field(map, "capabilities"))
        .and_then(arr)
    else {
        return;
    };
    for item in items {
        let Some(map) = obj(item) else { continue };
        let id = field(map, "id").and_then(text).unwrap_or("<no-id>");
        let Some(evidence) = field(map, "evidence").and_then(arr) else {
            continue;
        };
        for entry in evidence {
            let Some(entry) = obj(entry) else { continue };
            let (Some(kind), Some(reference)) = (
                field(entry, "kind").and_then(text),
                field(entry, "ref").and_then(text),
            ) else {
                continue;
            };
            if !RESOLVABLE_KINDS.contains(&kind) {
                continue;
            }
            let (path, symbol) = match reference.split_once("::") {
                Some((path, symbol)) => (path, Some(symbol)),
                None => (reference, None),
            };
            let full = root.join(path);
            if !full.is_file() {
                violations.push(violation(
                    id,
                    format!(
                        "capability {id:?} cites {kind} evidence {reference:?} but {path} does not \
                         exist — the registry may not cite evidence that is not there"
                    ),
                ));
                continue;
            }
            let Some(symbol) = symbol else { continue };
            let Ok(body) = std::fs::read_to_string(&full) else {
                violations.push(violation(
                    id,
                    format!(
                        "capability {id:?} cites {reference:?} but {path} is not readable UTF-8"
                    ),
                ));
                continue;
            };
            if !body.contains(&format!("fn {symbol}")) {
                violations.push(violation(
                    id,
                    format!(
                        "capability {id:?} cites {kind} evidence {reference:?} but {path} contains \
                         no `fn {symbol}` — a renamed or deleted test silently voids the level it \
                         justifies"
                    ),
                ));
            }
        }
    }
}

/// Level bars that are mechanically checkable from the evidence kinds present.
/// The qualitative bars in `docs/MATURITY_LEVELS.md` (independent oracle,
/// stated coverage, written support policy) are reviewer obligations; this
/// check only enforces the parts a machine can see.
struct ReceiptRun {
    head: String,
    script: String,
}

fn receipt_stage_run(source: &str, capability: &str, stage: &str) -> Result<ReceiptRun, String> {
    let parsed = JsonParser::new(source)
        .finish()
        .map_err(|error| format!("receipt is not valid JSON: {error}"))?;
    let root = obj(&parsed).ok_or_else(|| "receipt is not a JSON object".to_string())?;
    match field(root, "schema").and_then(text) {
        Some(SPINE_E2E_RECEIPT_SCHEMA) => {}
        Some(schema) => {
            return Err(format!(
                "receipt schema {schema:?} is not the admitted {SPINE_E2E_RECEIPT_SCHEMA:?}"
            ));
        }
        None => return Err("receipt has no schema string".to_string()),
    }
    let run = field(root, "run")
        .and_then(obj)
        .ok_or_else(|| "receipt has no run object".to_string())?;
    let head = field(run, "head_sha")
        .and_then(text)
        .filter(|head| !head.is_empty())
        .ok_or_else(|| "receipt has no non-empty run.head_sha".to_string())?;
    let script = field(run, "script")
        .and_then(text)
        .filter(|script| !script.is_empty())
        .ok_or_else(|| "receipt has no non-empty run.script".to_string())?;
    let stages = field(root, "stages")
        .and_then(arr)
        .ok_or_else(|| "receipt has no stages array".to_string())?;
    let matching: Vec<&BTreeMap<String, JsonValue>> = stages
        .iter()
        .filter_map(obj)
        .filter(|stage_entry| {
            matches!(
                (
                    field(stage_entry, "capability").and_then(text),
                    field(stage_entry, "stage").and_then(text),
                ),
                (Some(receipt_capability), Some(receipt_stage))
                    if receipt_capability == capability && receipt_stage == stage
            )
        })
        .collect();
    match matching.as_slice() {
        [stage_entry] if field(*stage_entry, "status").and_then(text) == Some("executed") => {
            Ok(ReceiptRun {
                head: head.to_string(),
                script: script.to_string(),
            })
        }
        [] => Err(format!(
            "receipt does not list capability {capability:?} stage {stage:?} as executed"
        )),
        [_] => Err(format!(
            "receipt does not mark capability {capability:?} stage {stage:?} as executed"
        )),
        _ => Err(format!(
            "receipt lists capability {capability:?} stage {stage:?} more than once"
        )),
    }
}

fn receipt_head_covers_lanes(
    root: &Path,
    receipt_head: &str,
    receipt_script: &str,
    lanes: &[String],
) -> Result<(), String> {
    if !lanes.is_empty() && !lanes.iter().any(|lane| lane == receipt_script) {
        return Err(format!(
            "receipt run.script {receipt_script:?} does not match any cited lane {lanes:?}"
        ));
    }
    for lane in lanes {
        let output = std::process::Command::new("git")
            .args(["log", "-1", "--format=%H", "--", lane])
            .current_dir(root)
            .output()
            .map_err(|error| format!("cannot read last commit for cited lane {lane:?}: {error}"))?;
        let lane_head = String::from_utf8(output.stdout)
            .map_err(|_| format!("git returned non-UTF-8 for cited lane {lane:?}"))?;
        let lane_head = lane_head.trim();
        if lane_head.is_empty() {
            return Err(format!("cited lane {lane:?} has no committed history"));
        }
        let covered = std::process::Command::new("git")
            .args(["merge-base", "--is-ancestor", lane_head, receipt_head])
            .current_dir(root)
            .status()
            .map_err(|error| {
                format!("cannot compare receipt HEAD to cited lane {lane:?}: {error}")
            })?;
        if !covered.success() {
            return Err(format!(
                "receipt HEAD {receipt_head} is not at or after cited lane {lane:?} commit {lane_head}"
            ));
        }
    }
    let current_head = std::process::Command::new("git")
        .args(["rev-parse", "--verify", "HEAD^{commit}"])
        .current_dir(root)
        .output()
        .map_err(|error| format!("cannot read current HEAD for receipt admission: {error}"))?;
    let current_head = String::from_utf8(current_head.stdout)
        .map_err(|_| "git returned non-UTF-8 for current HEAD".to_string())?;
    let current_head = current_head.trim();
    if current_head.is_empty() {
        return Err("current HEAD has no commit for receipt admission".to_string());
    }
    let reachable = std::process::Command::new("git")
        .args(["merge-base", "--is-ancestor", receipt_head, current_head])
        .current_dir(root)
        .status()
        .map_err(|error| format!("cannot compare receipt HEAD to current HEAD: {error}"))?;
    if !reachable.success() {
        return Err(format!(
            "receipt HEAD {receipt_head} is not reachable from current HEAD {current_head}"
        ));
    }
    Ok(())
}

fn check_level_bars(
    root: &Path,
    entries: &[Entry],
    violations: &mut Vec<Violation>,
    decisions: &mut Vec<PolicyNote>,
) {
    for entry in entries {
        let Some(index) = level_index(&entry.level) else {
            continue;
        };
        // L2+ : must cite at least one resolvable test.
        if index >= 1 && !entry.kinds.contains("test") {
            violations.push(violation(
                &entry.id,
                format!(
                    "capability {:?} claims {} but cites no `test` evidence; L2 and above require \
                     a named, resolvable test (docs/MATURITY_LEVELS.md)",
                    entry.id, entry.level
                ),
            ));
        }
        // L3+ : must cite an e2e lane and a retained receipt for this exact
        // capability stage. A script exists only as an intention until its
        // recorded run proves the stage executed at a current-enough HEAD.
        if index >= 2 && !entry.kinds.contains("lane") {
            violations.push(violation(
                &entry.id,
                format!(
                    "capability {:?} claims {} but cites no `lane` evidence; L3 and above require \
                     an end-to-end lane, not only unit tests",
                    entry.id, entry.level
                ),
            ));
        }
        if index >= 2 && entry.receipts.is_empty() {
            violations.push(violation(
                &entry.id,
                format!(
                    "capability {:?} claims {} but cites no `receipt` evidence; L3 and above \
                     require a retained executed-stage receipt",
                    entry.id, entry.level
                ),
            ));
        }
        if index >= 2 {
            for receipt in &entry.receipts {
                let verdict = std::fs::read_to_string(root.join(&receipt.path))
                    .map_err(|error| format!("cannot read receipt {:?}: {error}", receipt.path))
                    .and_then(|source| receipt_stage_run(&source, &entry.id, &receipt.stage))
                    .and_then(|run| {
                        receipt_head_covers_lanes(root, &run.head, &run.script, &entry.lanes)?;
                        Ok(run)
                    });
                match verdict {
                    Ok(run) => decisions.push(note(
                        &entry.id,
                        "executed-receipt",
                        format!(
                            "L3 receipt={} HEAD={} stage={} verdict=executed",
                            receipt.path, run.head, receipt.stage
                        ),
                    )),
                    Err(error) => violations.push(violation(
                        &entry.id,
                        format!(
                            "capability {:?} L3 receipt {:?} stage {:?} is not admitted: {error}",
                            entry.id, receipt.path, receipt.stage
                        ),
                    )),
                }
            }
        }
        // L4+ : must cite corpus validation.
        if index >= 3 && !entry.kinds.contains("corpus") {
            violations.push(violation(
                &entry.id,
                format!(
                    "capability {:?} claims {} but cites no `corpus` evidence; L4 is validation \
                     against an external corpus over a stated domain",
                    entry.id, entry.level
                ),
            ));
        }
        // L5 : must cite a written support policy document.
        if index >= 4 && !entry.kinds.contains("doc") {
            violations.push(violation(
                &entry.id,
                format!(
                    "capability {:?} claims L5 but cites no `doc` evidence; L5 requires a written \
                     support policy, not an intention",
                    entry.id
                ),
            ));
        }
    }
}

/// Compare against the last committed registry and classify each level change.
/// Promotions are the signal bead `.2.3` gates on; demotions always pass.
fn check_transitions(root: &Path, current: &[Entry], decisions: &mut Vec<PolicyNote>) {
    let output = std::process::Command::new("git")
        .args(["show", &format!("HEAD:{REGISTRY_FILE}")])
        .current_dir(root)
        .output();
    let Ok(output) = output else { return };
    if !output.status.success() {
        // No committed registry yet: the whole file is new, which is not a
        // promotion of anything.
        decisions.push(note(
            "<repo>",
            "baseline",
            format!("{REGISTRY_FILE} has no committed predecessor; recording the initial baseline"),
        ));
        return;
    }
    let Ok(previous_text) = String::from_utf8(output.stdout) else {
        return;
    };
    let mut ignored = Vec::new();
    let previous = parse_registry(&previous_text, "<committed>", &mut ignored);
    let baseline: BTreeMap<&str, &str> = previous
        .iter()
        .map(|entry| (entry.id.as_str(), entry.level.as_str()))
        .collect();

    for entry in current {
        let Some(index) = level_index(&entry.level) else {
            continue;
        };
        match baseline.get(entry.id.as_str()) {
            None => decisions.push(note(
                &entry.id,
                "introduced",
                format!(
                    "capability {:?} is new to the registry at {}",
                    entry.id, entry.level
                ),
            )),
            Some(before) => {
                let Some(before_index) = level_index(before) else {
                    continue;
                };
                if index > before_index {
                    decisions.push(note(
                        &entry.id,
                        "promotion",
                        format!(
                            "capability {:?} is being PROMOTED {before} -> {}; the claim-integrity \
                             gate (bead f85xj.2.3) must clear this before it lands",
                            entry.id, entry.level
                        ),
                    ));
                } else if index < before_index {
                    decisions.push(note(
                        &entry.id,
                        "demotion",
                        format!(
                            "capability {:?} is being DEMOTED {before} -> {}; demotions are always \
                             allowed and are logged, never blocked",
                            entry.id, entry.level
                        ),
                    ));
                }
            }
        }
    }
    for (id, before) in &baseline {
        if !current.iter().any(|entry| entry.id == *id) {
            decisions.push(note(
                id,
                "withdrawn",
                format!("capability {id:?} was removed from the registry (was {before})"),
            ));
        }
    }
}

/// Registry check: see the module docs for the five rules.
pub fn check_maturity(root: &Path) -> MaturityReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();

    let path = root.join(REGISTRY_FILE);
    let Ok(source) = std::fs::read_to_string(&path) else {
        violations.push(violation(
            "<repo>",
            format!(
                "{REGISTRY_FILE} is missing — capability maturity is the spine of program \
                 governance and the claim-integrity promotion gate reads it (bead f85xj.16.1)"
            ),
        ));
        return MaturityReport {
            violations,
            decisions,
        };
    };

    let entries = parse_registry(&source, REGISTRY_FILE, &mut violations);
    if violations.is_empty() {
        let projection = check_readme_projection_entries(root, &entries);
        violations.extend(projection.violations);
        decisions.extend(projection.decisions);
    }
    resolve_refs(root, &source, &mut violations);
    check_level_bars(root, &entries, &mut violations, &mut decisions);
    check_transitions(root, &entries, &mut decisions);

    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for entry in &entries {
        *counts.entry(entry.level.as_str()).or_default() += 1;
    }
    decisions.push(note(
        "<repo>",
        "inventory",
        format!(
            "{} capabilities recorded: {}",
            entries.len(),
            LEVELS
                .iter()
                .map(|level| format!("{level}={}", counts.get(level).copied().unwrap_or(0)))
                .collect::<Vec<_>>()
                .join(" ")
        ),
    ));

    MaturityReport {
        violations,
        decisions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn registry(level: &str, evidence: &str) -> String {
        format!(
            r#"{{"schema":"{REGISTRY_SCHEMA}","capabilities":[
                {{"id":"a.b","title":"T","owner":"o","level":"{level}",
                  "last_review":"2026-07-22","crates":["fs-x"],
                  "evidence":[{evidence}]}}]}}"#
        )
    }

    #[test]
    fn a_well_formed_entry_passes_structural_checks() {
        let mut v = Vec::new();
        let entries = parse_registry(
            &registry(
                "L1",
                r#"{"kind":"contract","ref":"crates/fs-x/CONTRACT.md"}"#,
            ),
            "t",
            &mut v,
        );
        assert!(v.is_empty(), "{v:?}");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].level, "L1");
    }

    #[test]
    fn structural_defects_are_each_caught() {
        // Bad schema.
        let mut v = Vec::new();
        parse_registry(r#"{"schema":"wrong","capabilities":[]}"#, "t", &mut v);
        assert!(
            v.iter().any(|x| x.detail.contains("declares schema")),
            "{v:?}"
        );

        // Not JSON.
        let mut v = Vec::new();
        parse_registry("{not json", "t", &mut v);
        assert!(
            v.iter().any(|x| x.detail.contains("not valid JSON")),
            "{v:?}"
        );

        // Unknown level.
        let mut v = Vec::new();
        parse_registry(
            &registry("L9", r#"{"kind":"doc","ref":"README.md"}"#),
            "t",
            &mut v,
        );
        assert!(
            v.iter().any(|x| x.detail.contains("expected one of")),
            "{v:?}"
        );

        // Bad date.
        let mut v = Vec::new();
        let bad =
            registry("L1", r#"{"kind":"doc","ref":"README.md"}"#).replace("2026-07-22", "22/07/26");
        parse_registry(&bad, "t", &mut v);
        assert!(v.iter().any(|x| x.detail.contains("YYYY-MM-DD")), "{v:?}");

        // Missing crate scope.
        let mut v = Vec::new();
        let unscoped =
            registry("L1", r#"{"kind":"doc","ref":"README.md"}"#).replace(r#"["fs-x"]"#, "[]");
        parse_registry(&unscoped, "t", &mut v);
        assert!(v.iter().any(|x| x.detail.contains("crates")), "{v:?}");

        // Unknown evidence kind.
        let mut v = Vec::new();
        parse_registry(
            &registry("L1", r#"{"kind":"vibes","ref":"x"}"#),
            "t",
            &mut v,
        );
        assert!(v.iter().any(|x| x.detail.contains("unknown kind")), "{v:?}");

        // Duplicate ids.
        let mut v = Vec::new();
        let dup = format!(
            r#"{{"schema":"{REGISTRY_SCHEMA}","capabilities":[
              {{"id":"d","title":"T","owner":"o","level":"L1","last_review":"2026-07-22",
                "crates":["c"],"evidence":[]}},
              {{"id":"d","title":"T","owner":"o","level":"L1","last_review":"2026-07-22",
                "crates":["c"],"evidence":[]}}]}}"#
        );
        parse_registry(&dup, "t", &mut v);
        assert!(
            v.iter().any(|x| x.detail.contains("more than once")),
            "{v:?}"
        );
    }

    #[test]
    fn level_bars_require_their_evidence_kinds() {
        let bar = |level: &str, kinds: &[&str]| {
            let entries = vec![Entry {
                id: "cap".to_string(),
                title: "Capability".to_string(),
                level: level.to_string(),
                crates: vec!["fs-cap".to_string()],
                notes: "Boundary".to_string(),
                kinds: kinds.iter().map(|k| (*k).to_string()).collect(),
                lanes: Vec::new(),
                receipts: Vec::new(),
            }];
            let mut v = Vec::new();
            let mut decisions = Vec::new();
            check_level_bars(Path::new("."), &entries, &mut v, &mut decisions);
            v
        };
        assert!(bar("L1", &[]).is_empty(), "L1 needs no test evidence");
        assert!(bar("L2", &[]).iter().any(|x| x.detail.contains("`test`")));
        assert!(bar("L2", &["test"]).is_empty());
        assert!(
            bar("L3", &["test"])
                .iter()
                .any(|x| x.detail.contains("`lane`"))
        );
        assert!(
            bar("L4", &["test", "lane"])
                .iter()
                .any(|x| x.detail.contains("`corpus`"))
        );
        assert!(
            bar("L5", &["test", "lane", "corpus"])
                .iter()
                .any(|x| x.detail.contains("support policy"))
        );
    }

    #[test]
    fn l3_receipt_requires_its_capability_stage_to_be_executed() {
        let root =
            std::env::temp_dir().join(format!("fsim-maturity-receipt-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("scripts/e2e")).unwrap();
        std::fs::write(root.join("scripts/e2e/cooling_01.sh"), "#!/bin/sh\n").unwrap();
        let git = |args: &[&str]| {
            let output = std::process::Command::new("git")
                .args(args)
                .current_dir(&root)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "git {args:?}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8(output.stdout).unwrap().trim().to_string()
        };
        git(&["init", "-q", "-b", "main"]);
        git(&["config", "user.email", "maturity-test@example.invalid"]);
        git(&["config", "user.name", "maturity-test"]);
        git(&["add", "scripts/e2e/cooling_01.sh"]);
        git(&["commit", "-qm", "lane"]);
        let lane_head = git(&["rev-parse", "HEAD"]);

        let mut parse_violations = Vec::new();
        let entries = parse_registry(
            &registry(
                "L3",
                r#"{"kind":"test","ref":"xtask/src/maturity.rs::l3_receipt_requires_its_capability_stage_to_be_executed"},{"kind":"lane","ref":"scripts/e2e/cooling_01.sh"},{"kind":"receipt","ref":"spine-e2e-summary.json","stage":"conduction"}"#,
            ),
            "synthetic",
            &mut parse_violations,
        );
        assert!(
            parse_violations.is_empty(),
            "synthetic registry: {parse_violations:?}"
        );
        let mut entries = entries;
        entries[0].receipts[0].path = "receipt.json".to_string();

        let missing_stage = format!(
            r#"{{
            "schema":"frankensim-spine-e2e-receipt-v1",
            "run":{{"head_sha":"{lane_head}","script":"scripts/e2e/cooling_01.sh"}},
            "stages":[{{"capability":"a.b","stage":"flow","status":"executed"}}]
        }}"#
        );
        std::fs::write(root.join("receipt.json"), missing_stage).unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations.iter().any(|violation| violation
                .detail
                .contains("does not list capability \"a.b\" stage \"conduction\" as executed")),
            "missing stage must be a gate violation: {violations:?}"
        );

        let executed = format!(
            r#"{{
            "schema":"frankensim-spine-e2e-receipt-v1",
            "run":{{"head_sha":"{lane_head}","script":"scripts/e2e/cooling_01.sh"}},
            "stages":[{{"capability":"a.b","stage":"conduction","status":"executed"}}]
        }}"#
        );
        std::fs::write(root.join("receipt.json"), &executed).unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations.is_empty(),
            "executed stage must pass the maturity gate: {violations:?}"
        );
        assert_eq!(
            decisions
                .iter()
                .filter(|decision| decision.verdict == "executed-receipt")
                .count(),
            1,
            "executed stage must emit one receipt decision: {decisions:?}"
        );

        let wrong_script =
            executed.replace("scripts/e2e/cooling_01.sh", "scripts/e2e/unrelated.sh");
        std::fs::write(root.join("receipt.json"), wrong_script).unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("run.script")),
            "a receipt from another script must be a gate violation: {violations:?}"
        );

        git(&["checkout", "-qb", "foreign"]);
        std::fs::write(root.join("foreign"), "divergent receipt\n").unwrap();
        git(&["add", "foreign"]);
        git(&["commit", "-qm", "foreign receipt"]);
        let foreign_head = git(&["rev-parse", "HEAD"]);
        git(&["checkout", "-q", "main"]);
        std::fs::write(
            root.join("receipt.json"),
            executed.replace(&lane_head, &foreign_head),
        )
        .unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("not reachable from current HEAD")),
            "a divergent receipt HEAD must be a gate violation: {violations:?}"
        );

        let foreign_schema =
            executed.replace("frankensim-spine-e2e-receipt-v1", "unadmitted-receipt-v1");
        std::fs::write(root.join("receipt.json"), foreign_schema).unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("is not the admitted")),
            "foreign receipt schema must be a gate violation: {violations:?}"
        );

        let duplicate_stage = r#"{
            "schema":"frankensim-spine-e2e-receipt-v1",
            "run":{"head_sha":"0123456789abcdef0123456789abcdef01234567","script":"scripts/e2e/cooling_01.sh"},
            "stages":[
                {"capability":"a.b","stage":"conduction","status":"executed"},
                {"capability":"a.b","stage":"conduction","status":"failed"}
            ]
        }"#;
        std::fs::write(root.join("receipt.json"), duplicate_stage).unwrap();
        let mut violations = Vec::new();
        let mut decisions = Vec::new();
        check_level_bars(&root, &entries, &mut violations, &mut decisions);
        assert!(
            violations
                .iter()
                .any(|violation| violation.detail.contains("more than once")),
            "duplicate stage rows must be a gate violation: {violations:?}"
        );
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn unresolvable_evidence_refs_are_violations() {
        let base = std::env::temp_dir().join(format!("fsim-maturity-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(base.join("crates/fs-x/tests")).unwrap();
        std::fs::write(base.join("crates/fs-x/tests/t.rs"), "fn real_case() {}\n").unwrap();

        // Present file + present symbol resolves.
        let mut v = Vec::new();
        resolve_refs(
            &base,
            &registry(
                "L2",
                r#"{"kind":"test","ref":"crates/fs-x/tests/t.rs::real_case"}"#,
            ),
            &mut v,
        );
        assert!(v.is_empty(), "{v:?}");

        // Missing symbol is caught — a renamed test voids the level.
        let mut v = Vec::new();
        resolve_refs(
            &base,
            &registry(
                "L2",
                r#"{"kind":"test","ref":"crates/fs-x/tests/t.rs::ghost_case"}"#,
            ),
            &mut v,
        );
        assert!(
            v.iter().any(|x| x.detail.contains("no `fn ghost_case`")),
            "{v:?}"
        );

        // Missing file is caught.
        let mut v = Vec::new();
        resolve_refs(
            &base,
            &registry(
                "L2",
                r#"{"kind":"test","ref":"crates/fs-x/tests/gone.rs::real_case"}"#,
            ),
            &mut v,
        );
        assert!(
            v.iter().any(|x| x.detail.contains("does not exist")),
            "{v:?}"
        );

        // corpus refs are recorded, never resolved.
        let mut v = Vec::new();
        resolve_refs(
            &base,
            &registry("L4", r#"{"kind":"corpus","ref":"no-such-dataset"}"#),
            &mut v,
        );
        assert!(v.is_empty(), "corpus refs must not be resolved yet: {v:?}");

        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn generated_capability_matrix_is_sorted_and_fails_on_one_stale_level() {
        let entries = vec![
            Entry {
                id: "z.last".to_string(),
                title: "Last".to_string(),
                level: "L1".to_string(),
                crates: vec!["fs-z".to_string()],
                notes: "Experimental boundary".to_string(),
                kinds: BTreeSet::new(),
                lanes: Vec::new(),
                receipts: Vec::new(),
            },
            Entry {
                id: "a.first".to_string(),
                title: "First".to_string(),
                level: "L2".to_string(),
                crates: vec!["fs-a".to_string(), "fs-b".to_string()],
                notes: "Verified against A | B".to_string(),
                kinds: BTreeSet::from(["test".to_string()]),
                lanes: Vec::new(),
                receipts: Vec::new(),
            },
        ];
        let generated = render_readme_matrix(&entries);
        assert!(
            generated.find("`a.first`").unwrap() < generated.find("`z.last`").unwrap(),
            "the projection is canonical by capability id"
        );
        assert!(generated.contains("A \\| B"), "Markdown pipes are escaped");
        assert!(check_readme_matrix_text(&generated, &entries).is_empty());

        let stale = generated.replacen("| L2 |", "| L3 |", 1);
        let violations = check_readme_matrix_text(&stale, &entries);
        assert_eq!(
            violations.len(),
            1,
            "one seeded maturity drift: {violations:?}"
        );
        assert!(violations[0].detail.contains("matrix is stale"));
    }

    #[test]
    fn readme_maturity_summary_is_exact_checked_against_registry_counts() {
        let entries = vec![
            Entry {
                id: "a".to_string(),
                title: "A".to_string(),
                level: "L1".to_string(),
                crates: vec!["fs-a".to_string()],
                notes: String::new(),
                kinds: BTreeSet::new(),
                lanes: Vec::new(),
                receipts: Vec::new(),
            },
            Entry {
                id: "b".to_string(),
                title: "B".to_string(),
                level: "L2".to_string(),
                crates: vec!["fs-b".to_string()],
                notes: String::new(),
                kinds: BTreeSet::new(),
                lanes: Vec::new(),
                receipts: Vec::new(),
            },
        ];
        let summary = concat!(
            "it registers 2 product-meaningful capabilities:\n",
            "| L1 | Experimental | 1 | boundary |\n",
            "| L2 | Verified | 1 | boundary |\n",
            "| L3 | Integrated | 0 | boundary |\n",
            "| L4 | Validated | 0 | boundary |\n",
            "| L5 | Supported | 0 | boundary |\n",
        );
        assert!(check_readme_summary_counts(summary, &entries).is_empty());
        let stale = summary.replacen("| L2 | Verified | 1 |", "| L2 | Verified | 9 |", 1);
        let violations = check_readme_summary_counts(&stale, &entries);
        assert_eq!(
            violations.len(),
            1,
            "one stale summary count: {violations:?}"
        );
        assert!(violations[0].detail.contains("L2"));

        let duplicated_total =
            format!("{summary}it also registers 2 product-meaningful capabilities:\n");
        let violations = check_readme_summary_counts(&duplicated_total, &entries);
        assert_eq!(
            violations.len(),
            1,
            "one duplicate-total defect: {violations:?}"
        );
        assert!(violations[0].detail.contains("exactly one"));
    }

    #[test]
    fn the_live_registry_is_clean() {
        // The repo's own registry must satisfy every rule; this is the check
        // that keeps the shipped file honest as capabilities move.
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let report = check_maturity(root);
        assert!(
            report.violations.is_empty(),
            "live capability-maturity.json must be clean: {:?}",
            report.violations
        );
        // The inventory the registry reports today: one L3 admitted on the
        // retained lane receipt (bead frankensim-rc-root-q61wp.13). A change
        // here must come with the registry change that justifies it.
        assert!(
            report
                .decisions
                .iter()
                .any(|note| note.verdict == "inventory" && note.detail.contains("L1=3 L2=11 L3=1")),
            "the live registry must report the L1=3 L2=11 L3=1 inventory: {:?}",
            report.decisions
        );
        assert!(
            report.decisions.iter().any(|note| {
                note.verdict == "executed-receipt"
                    && note.crate_name == "thermal.conduction-solve"
                    && note.detail.contains("stage=conduction verdict=executed")
            }),
            "the L3 must be admitted from an executed-stage receipt, never by assertion: {:?}",
            report.decisions
        );
    }
}
