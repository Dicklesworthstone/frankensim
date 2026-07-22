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
use std::path::Path;

pub const REGISTRY_FILE: &str = "capability-maturity.json";
const REGISTRY_SCHEMA: &str = "frankensim-capability-maturity-v1";
const CHECK: &str = "capability-maturity";
const LEVELS: [&str; 5] = ["L1", "L2", "L3", "L4", "L5"];
const README_MATRIX_BEGIN: &str = "<!-- BEGIN GENERATED FRANKENSIM CAPABILITY MATRIX -->";
const README_MATRIX_END: &str = "<!-- END GENERATED FRANKENSIM CAPABILITY MATRIX -->";

/// Evidence kinds and whether the check can resolve them against the tree.
/// `corpus` is recorded but unresolvable until the V&V corpus registry (e04)
/// exists — an honest gap, and the reason nothing is L4 today.
const RESOLVABLE_KINDS: [&str; 4] = ["test", "contract", "lane", "doc"];
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
        let kinds = collect_evidence(map, id, violations);
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
        output.push_str(&format!(
            "| `{}` — {} | {} | {crates} | {} |\n",
            markdown_cell(&entry.id),
            markdown_cell(&entry.title),
            entry.level,
            markdown_cell(&entry.notes),
        ));
    }
    output.push_str(README_MATRIX_END);
    output
}

fn readme_matrix_block<'a>(source: &'a str) -> Result<&'a str, String> {
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
) -> BTreeSet<String> {
    let mut kinds = BTreeSet::new();
    let Some(items) = field(map, "evidence").and_then(arr) else {
        violations.push(violation(
            id,
            format!("capability {id:?} has no \"evidence\" array"),
        ));
        return kinds;
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
    }
    kinds
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
fn check_level_bars(entries: &[Entry], violations: &mut Vec<Violation>) {
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
        // L3+ : must cite an e2e lane.
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
    check_level_bars(&entries, &mut violations);
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
            }];
            let mut v = Vec::new();
            check_level_bars(&entries, &mut v);
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
            },
            Entry {
                id: "a.first".to_string(),
                title: "First".to_string(),
                level: "L2".to_string(),
                crates: vec!["fs-a".to_string(), "fs-b".to_string()],
                notes: "Verified against A | B".to_string(),
                kinds: BTreeSet::from(["test".to_string()]),
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
            },
            Entry {
                id: "b".to_string(),
                title: "B".to_string(),
                level: "L2".to_string(),
                crates: vec!["fs-b".to_string()],
                notes: String::new(),
                kinds: BTreeSet::new(),
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
        assert!(
            report
                .decisions
                .iter()
                .any(|note| note.verdict == "inventory"),
            "an inventory note is always emitted"
        );
    }
}
