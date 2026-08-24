//! Instrument-claims registry check (bead `frankensim-music-claims-registry-mc31g`).
//!
//! The music program's doctrine as code: a musical claim is a
//! `(filling, image, qoi)` row in `instrument-claims.json`, never prose.
//!
//! - D19: claims are rows; every row names its owner crates and exactness
//!   class, so "clarinet TMM hits Ernoult cents" and "modal ZOH renders a
//!   pluck" are queryable objects rather than folklore.
//! - D21 (menus, not winners): a row may TRANSITION gate status
//!   (`ungated -> green -> refused`, and back when new evidence appears);
//!   it may never be deleted relative to the committed predecessor.
//!   Deleting a passing orthogonal image is how registries lie, so removal
//!   is a violation, not an edit.
//! - D25: a `live_default` image without BOTH a green gate and a measured
//!   budget row is a policy violation this check refuses. Real-time claims
//!   are measured headroom, never vibes.
//!
//! What this check does NOT do: judge whether cited evidence actually earns
//! a gate (that is the per-track gates beads' review job), or resolve
//! receipt/bake-off/listening/corpus references whose stores are later
//! beads (recorded-only, exactly like maturity's `corpus` kind). It proves
//! the paperwork is present, internally consistent, and never silently
//! shrunk. Claiming more would make the gate itself a claim-integrity
//! defect.
//!
//! Parser note: the in-house JSON parser collapses `true`/`false` into a
//! payload-less `Bool`, so flags in this registry are string enums
//! (`"live_default": "yes" | "no"`) — deterministic, greppable, lossless.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub const REGISTRY_FILE: &str = "instrument-claims.json";
pub const CHECK: &str = "instrument-claims";
const REGISTRY_SCHEMA: &str = "frankensim-instrument-claims-v1";
const GATE_VALUES: [&str; 3] = ["ungated", "green", "refused"];
const EXACTNESS_VALUES: [&str; 4] = ["X-Exact", "X-Consist", "X-Struct", "X-Est"];
const DETERMINISM_VALUES: [&str; 5] = [
    "one-host",
    "same-isa",
    "cross-isa",
    "statistical",
    "fast-mode",
];
const LIVE_DEFAULT_VALUES: [&str; 2] = ["yes", "no"];
/// Evidence kinds whose `ref` is a tree path (optionally `path::item`); the
/// path part must exist in the tree.
const RESOLVABLE_EVIDENCE_KINDS: [&str; 3] = ["test", "contract", "doc"];
/// Evidence kinds recorded but not yet resolvable — their stores are later
/// beads (bake-off receipts, listening receipts, corpus rows, budget rows).
/// An honest gap, mirrored on maturity's `corpus` kind.
const RECORDED_ONLY_EVIDENCE_KINDS: [&str; 4] = ["receipt", "bakeoff", "listening", "corpus"];

/// Corpus-reference lint (bead `frankensim-music-v8-root-3ez8g.1.1`).
///
/// `corpus_refs` entries with the `vvreg:` prefix must name a case id
/// registered in the tracked acoustic corpus manifest — a gate citing an
/// unregistered corpus is citing nothing. Refs without the prefix are
/// refused outright: every music corpus registers through fs-vvreg (the
/// reuse law — no parallel music-validation registry), so there is no other
/// legal namespace. Absence rows (`absent-hunt` / `refused-retention`) ARE
/// resolvable targets: citing one records "this gate is data-blocked on a
/// named hunt", which is exactly what the population signal is for.
pub const ACOUSTIC_MANIFEST_FILE: &str = "data/vv-corpus/acoustic/acoustic-v1.tsv";
const CORPUS_REF_PREFIX: &str = "vvreg:";

/// Parse the case-id column out of the tracked acoustic manifest.
fn manifest_case_ids(manifest: &str) -> BTreeSet<String> {
    manifest
        .lines()
        .skip(1) // header
        .filter_map(|line| line.split('\t').nth(1))
        .filter(|id| !id.trim().is_empty())
        .map(str::to_string)
        .collect()
}

/// Determinism composition lint (bead `frankensim-music-v8-root-3ez8g.1.3`).
///
/// Replay strength composes like evidence colour: a composed image inherits
/// the WEAKEST operand in its owner chain, and a row may never claim more.
/// The ladder is `one-host < same-isa < cross-isa`; `statistical` and
/// `fast-mode` are orthogonal declarations (replay-in-distribution and
/// declared-deviation respectively), exempt from the ladder but surfaced as
/// notes so reviewers see them.
///
/// Ceilings are sourced from crate CONTRACTs and recorded golden evidence,
/// never from optimism. Current truth (2026-08-14): ZERO cross-ISA goldens
/// exist anywhere in the music stack, so the DEFAULT ceiling is `one-host`.
/// The promotion path is golden evidence: when the cross-ISA audit (bead
/// 3ez8g.13.4) records matching digests on both reference ISA families for
/// a crate, its ceiling row is raised HERE, citing the goldens, in the same
/// commit (the golden-bump protocol).
const DETERMINISM_DEFAULT_CEILING: &str = "one-host";
/// Per-crate ceilings that differ from the default, with the reason.
/// fs-tribo: its CONTRACT declares NO cross-ISA bit-stability (norms use the
/// platform `hypot` sequence for overflow safety) — a deliberate, permanent
/// cap until bead 3ez8g.7.3 routes or re-declares it. Listing it explicitly
/// (even while equal to the default) makes the cap survive any future
/// default raise.
const DETERMINISM_CEILINGS: [(&str, &str); 6] = [
    ("fs-tribo", "one-host"),
    // Cross-ISA audit bead 3ez8g.13.4 (2026-08-23): matching digests on
    // aarch64-apple and x86_64-linux in BOTH build modes from the committed
    // tests/cross_isa_golden.rs per crate —
    // fs-phs 0x798c84cbeb3c39b9 (step ledger: duffing + 3-mode modal bank),
    // fs-vfit 0xd00d69e2b740e56b (vector fit -> bilinear -> filter steps),
    // fs-modal 0x11dced94c7f67115 (certified spectral slice windows),
    // fs-couple 0x432320b3c2bf06d9 (exact-ZOH modal render trajectory),
    // fs-duct 0xd6e9724c5414cf8d (TMM impedance sweep + peak finder).
    ("fs-duct", "cross-isa"),
    ("fs-vfit", "cross-isa"),
    ("fs-phs", "cross-isa"),
    ("fs-modal", "cross-isa"),
    ("fs-couple", "cross-isa"),
];

/// Ladder position; `None` for the orthogonal declarations.
fn determinism_strength(class: &str) -> Option<u8> {
    match class {
        "one-host" => Some(0),
        "same-isa" => Some(1),
        "cross-isa" => Some(2),
        _ => None,
    }
}

fn determinism_ceiling(owner: &str) -> &'static str {
    DETERMINISM_CEILINGS
        .iter()
        .find(|(name, _)| *name == owner)
        .map_or(DETERMINISM_DEFAULT_CEILING, |(_, ceiling)| ceiling)
}

pub struct InstrumentClaimsReport {
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

/// One registry row, reduced to what the check reasons about. (Exactness
/// classes are validated during parse; no downstream rule reads them yet —
/// the determinism-rows lint, bead 3ez8g.1.3, will.)
struct Row {
    key: String,
    owner_crates: Vec<String>,
    gate: String,
    live_default: String,
    determinism: String,
    evidence: Vec<(String, String)>,
    budget_row: Option<String>,
    corpus_refs: Vec<String>,
    notes: String,
}

/// Parse the registry into rows, pushing a violation for every structural
/// defect. Returns whatever parsed cleanly so one bad row does not hide the
/// rest (the maturity-check convention).
fn parse_registry(source: &str, violations: &mut Vec<Violation>) -> Vec<Row> {
    let parsed = match JsonParser::new(source).finish() {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(
                REGISTRY_FILE,
                format!("{REGISTRY_FILE} is not valid JSON: {error}"),
            ));
            return Vec::new();
        }
    };
    let Some(top) = obj(&parsed) else {
        violations.push(violation(
            REGISTRY_FILE,
            format!("{REGISTRY_FILE} is not a JSON object"),
        ));
        return Vec::new();
    };
    match top.get("schema").and_then(text) {
        Some(REGISTRY_SCHEMA) => {}
        Some(other) => violations.push(violation(
            REGISTRY_FILE,
            format!("{REGISTRY_FILE} declares schema {other:?}, expected {REGISTRY_SCHEMA:?}"),
        )),
        None => violations.push(violation(
            REGISTRY_FILE,
            format!("{REGISTRY_FILE} has no string \"schema\" field"),
        )),
    }
    let Some(rows) = top.get("rows").and_then(arr) else {
        violations.push(violation(
            REGISTRY_FILE,
            format!("{REGISTRY_FILE} has no \"rows\" array"),
        ));
        return Vec::new();
    };

    let mut out = Vec::new();
    for (index, value) in rows.iter().enumerate() {
        let entity = format!("{REGISTRY_FILE}#rows[{index}]");
        let Some(map) = obj(value) else {
            violations.push(violation(&entity, "row is not a JSON object".to_string()));
            continue;
        };
        let mut broken = false;
        let mut required = |field: &str| -> String {
            match map.get(field).and_then(text) {
                Some(s) if !s.trim().is_empty() => s.to_string(),
                Some(_) => {
                    violations.push(violation(
                        &entity,
                        format!("field \"{field}\" is empty; every row field is load-bearing"),
                    ));
                    broken = true;
                    String::new()
                }
                None => {
                    violations.push(violation(
                        &entity,
                        format!("field \"{field}\" is missing or not a string"),
                    ));
                    broken = true;
                    String::new()
                }
            }
        };
        let filling = required("filling");
        let image = required("image");
        let qoi = required("qoi");
        let gate = required("gate");
        let live_default = required("live_default");
        let determinism = required("determinism");
        let key = format!("{filling}/{image}/{qoi}");
        let entity = if broken { entity } else { key.clone() };

        if !GATE_VALUES.contains(&gate.as_str()) && !gate.is_empty() {
            violations.push(violation(
                &entity,
                format!("gate {gate:?} is not one of {GATE_VALUES:?}"),
            ));
        }
        if !LIVE_DEFAULT_VALUES.contains(&live_default.as_str()) && !live_default.is_empty() {
            violations.push(violation(
                &entity,
                format!(
                    "live_default {live_default:?} is not one of {LIVE_DEFAULT_VALUES:?} \
                     (string enum; the in-house parser collapses JSON booleans)"
                ),
            ));
        }
        if !DETERMINISM_VALUES.contains(&determinism.as_str()) && !determinism.is_empty() {
            violations.push(violation(
                &entity,
                format!("determinism {determinism:?} is not one of {DETERMINISM_VALUES:?}"),
            ));
        }

        let mut owner_crates = Vec::new();
        match map.get("owner_crates").and_then(arr) {
            Some(items) if !items.is_empty() => {
                for item in items {
                    match text(item) {
                        Some(name) if !name.trim().is_empty() => {
                            owner_crates.push(name.to_string());
                        }
                        _ => violations.push(violation(
                            &entity,
                            "owner_crates entries must be non-empty strings".to_string(),
                        )),
                    }
                }
            }
            _ => violations.push(violation(
                &entity,
                "owner_crates must be a non-empty array (D23: every image names its owner)"
                    .to_string(),
            )),
        }

        match map.get("exactness").and_then(arr) {
            Some(items) if (1..=2).contains(&items.len()) => {
                let mut seen = BTreeSet::new();
                for item in items {
                    match text(item) {
                        Some(class) if EXACTNESS_VALUES.contains(&class) => {
                            if !seen.insert(class.to_string()) {
                                violations.push(violation(
                                    &entity,
                                    format!("exactness class {class:?} listed twice"),
                                ));
                            }
                        }
                        other => violations.push(violation(
                            &entity,
                            format!("exactness entry {other:?} is not one of {EXACTNESS_VALUES:?}"),
                        )),
                    }
                }
            }
            Some(items) => violations.push(violation(
                &entity,
                format!(
                    "exactness has {} entries; a row states 1 or 2 classes (2 when an image \
                     splits its limit from its loss model, e.g. exact delay + ZK losses)",
                    items.len()
                ),
            )),
            None => violations.push(violation(
                &entity,
                "exactness must be an array of 1..=2 classes".to_string(),
            )),
        }

        let mut evidence = Vec::new();
        match map.get("evidence").and_then(arr) {
            Some(items) => {
                for item in items {
                    let Some(entry) = obj(item) else {
                        violations.push(violation(
                            &entity,
                            "evidence entries must be {kind, ref} objects".to_string(),
                        ));
                        continue;
                    };
                    let kind = entry.get("kind").and_then(text).unwrap_or_default();
                    let reference = entry.get("ref").and_then(text).unwrap_or_default();
                    let known = RESOLVABLE_EVIDENCE_KINDS.contains(&kind)
                        || RECORDED_ONLY_EVIDENCE_KINDS.contains(&kind);
                    if !known {
                        violations.push(violation(
                            &entity,
                            format!(
                                "evidence kind {kind:?} is not one of {RESOLVABLE_EVIDENCE_KINDS:?} \
                                 (resolved) or {RECORDED_ONLY_EVIDENCE_KINDS:?} (recorded-only)"
                            ),
                        ));
                    }
                    if reference.trim().is_empty() {
                        violations.push(violation(
                            &entity,
                            "evidence \"ref\" is missing or empty".to_string(),
                        ));
                    }
                    evidence.push((kind.to_string(), reference.to_string()));
                }
            }
            None => violations.push(violation(
                &entity,
                "evidence must be an array (empty is fine for ungated rows)".to_string(),
            )),
        }

        let budget_row = match map.get("budget_row") {
            Some(JsonValue::Null) | None => None,
            Some(JsonValue::String(s)) if !s.trim().is_empty() => Some(s.clone()),
            Some(_) => {
                violations.push(violation(
                    &entity,
                    "budget_row must be null or a non-empty content-address string".to_string(),
                ));
                None
            }
        };

        let mut corpus_refs = Vec::new();
        match map.get("corpus_refs").and_then(arr) {
            Some(items) => {
                for item in items {
                    match text(item) {
                        Some(reference) if !reference.trim().is_empty() => {
                            corpus_refs.push(reference.to_string());
                        }
                        _ => violations.push(violation(
                            &entity,
                            "corpus_refs entries must be non-empty strings".to_string(),
                        )),
                    }
                }
            }
            None => violations.push(violation(
                &entity,
                "corpus_refs must be an array (empty until corpora register)".to_string(),
            )),
        }

        let notes = map
            .get("notes")
            .and_then(text)
            .unwrap_or_default()
            .to_string();

        if broken {
            continue;
        }
        out.push(Row {
            key,
            owner_crates,
            gate,
            live_default,
            determinism,
            evidence,
            budget_row,
            corpus_refs,
            notes,
        });
    }
    out
}

/// Reduce a registry source to `key -> (gate, live_default)` for the
/// predecessor diff. Parse failures collapse to an empty map: a predecessor
/// that never existed or never parsed constrains nothing (the maturity
/// convention — a new registry is not a deletion).
fn keyed_gates(source: &str) -> BTreeMap<String, (String, String)> {
    let mut ignored = Vec::new();
    parse_registry(source, &mut ignored)
        .into_iter()
        .map(|row| (row.key, (row.gate, row.live_default)))
        .collect()
}

/// Pure core: validate `current` against the rules and its `committed`
/// predecessor, resolving tree references through the injected closures so
/// unit tests need neither a filesystem nor git.
fn check_registry_sources(
    current: &str,
    committed: Option<&str>,
    path_exists: &dyn Fn(&str) -> bool,
    owner_exists: &dyn Fn(&str) -> bool,
    corpus_ids: Option<&BTreeSet<String>>,
) -> InstrumentClaimsReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();
    let rows = parse_registry(current, &mut violations);

    // Key uniqueness: a duplicated (filling, image, qoi) makes gate status
    // ambiguous, which is worse than absent.
    let mut seen = BTreeSet::new();
    for row in &rows {
        if !seen.insert(row.key.clone()) {
            violations.push(violation(
                &row.key,
                "duplicate (filling, image, qoi) key; gate status must be unambiguous".to_string(),
            ));
        }
    }

    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut live_defaults = 0usize;
    for row in &rows {
        *counts.entry(row.gate.as_str()).or_insert(0) += 1;
        if row.live_default == "yes" {
            live_defaults += 1;
        }

        // D25: live defaults need a green gate AND a measured budget row.
        if row.live_default == "yes" {
            if row.gate != "green" {
                violations.push(violation(
                    &row.key,
                    format!(
                        "live_default row has gate {:?}; D25 requires gate=green before an \
                         image may be the runtime default",
                        row.gate
                    ),
                ));
            }
            if row.budget_row.is_none() {
                violations.push(violation(
                    &row.key,
                    "live_default row has no budget_row; D25: real-time is a measured claim \
                     (samples/sec + headroom on a named machine), never prose"
                        .to_string(),
                ));
            }
        }
        // Lifecycle statics: green needs evidence; refused needs its reason.
        if row.gate == "green" && row.evidence.is_empty() {
            violations.push(violation(
                &row.key,
                "gate=green with empty evidence; a green gate cites the tests/receipts that \
                 earned it"
                    .to_string(),
            ));
        }
        if row.gate == "refused" && row.notes.trim().is_empty() {
            violations.push(violation(
                &row.key,
                "gate=refused without a reason in notes; a refusal records why (and the \
                 counter-argument), per audit discipline"
                    .to_string(),
            ));
        }
        // D23: owners must exist as workspace crates.
        for owner in &row.owner_crates {
            if !owner_exists(owner) {
                violations.push(violation(
                    &row.key,
                    format!("owner crate {owner:?} does not exist under crates/"),
                ));
            }
        }
        // Resolvable evidence must resolve; recorded-only kinds are honest gaps.
        for (kind, reference) in &row.evidence {
            if RESOLVABLE_EVIDENCE_KINDS.contains(&kind.as_str()) {
                let path_part = reference.split("::").next().unwrap_or_default();
                if !path_exists(path_part) {
                    violations.push(violation(
                        &row.key,
                        format!(
                            "evidence {kind}:{reference:?} does not resolve (path {path_part:?} \
                             is not in the tree); stale evidence pointers are how gates rot"
                        ),
                    ));
                }
            }
        }
        // Corpus-reference lint (bead 3ez8g.1.1): vvreg: refs resolve
        // against the tracked acoustic manifest or refuse by name.
        for reference in &row.corpus_refs {
            let Some(case_id) = reference.strip_prefix(CORPUS_REF_PREFIX) else {
                violations.push(violation(
                    &row.key,
                    format!(
                        "corpus ref {reference:?} lacks the {CORPUS_REF_PREFIX} prefix; every \
                         music corpus registers through fs-vvreg (no parallel registry), so \
                         there is no other legal namespace"
                    ),
                ));
                continue;
            };
            match corpus_ids {
                Some(ids) if ids.contains(case_id) => {}
                Some(_) => violations.push(violation(
                    &row.key,
                    format!(
                        "corpus ref {reference:?} names no case in {ACOUSTIC_MANIFEST_FILE}; \
                         a gate citing an unregistered corpus is citing nothing"
                    ),
                )),
                None => violations.push(violation(
                    &row.key,
                    format!(
                        "corpus ref {reference:?} cannot be resolved: {ACOUSTIC_MANIFEST_FILE} \
                         is unreadable; a gate that cannot see the corpus registry refuses"
                    ),
                )),
            }
        }
        // Weakest-operand determinism lint (bead 3ez8g.1.3): a row's ladder
        // class may not exceed the weakest ceiling in its owner chain.
        match determinism_strength(&row.determinism) {
            Some(claimed) => {
                for owner in &row.owner_crates {
                    let ceiling = determinism_ceiling(owner);
                    let allowed = determinism_strength(ceiling).unwrap_or(0);
                    if claimed > allowed {
                        violations.push(violation(
                            &row.key,
                            format!(
                                "determinism {:?} exceeds owner crate {owner:?}'s ceiling \
                                 {ceiling:?}; replay strength composes like evidence colour \
                                 (weakest operand wins) and ceilings rise only on recorded \
                                 golden evidence (bead 3ez8g.13.4)",
                                row.determinism
                            ),
                        ));
                    }
                }
            }
            None => {
                // statistical / fast-mode: orthogonal declarations, exempt
                // from the ladder but surfaced so reviewers see them.
                decisions.push(note(
                    &row.key,
                    "review",
                    format!(
                        "row declares determinism {:?} (off the replay ladder); verify the \
                         declaration matches the owning contract's stated class",
                        row.determinism
                    ),
                ));
            }
        }
    }

    // D21: the committed predecessor's keys may never disappear.
    if let Some(committed) = committed {
        let previous = keyed_gates(committed);
        let current_keys: BTreeSet<&str> = rows.iter().map(|row| row.key.as_str()).collect();
        for (key, (old_gate, old_live)) in &previous {
            if !current_keys.contains(key.as_str()) {
                violations.push(violation(
                    key,
                    format!(
                        "row deleted (was gate={old_gate:?}); D21: menus, not winners — an \
                         image row may move to gate=refused but may never vanish"
                    ),
                ));
                continue;
            }
            if let Some(row) = rows.iter().find(|row| &row.key == key) {
                if &row.gate != old_gate {
                    decisions.push(note(
                        key,
                        "transition",
                        format!("gate {old_gate} -> {}", row.gate),
                    ));
                }
                if &row.live_default != old_live {
                    decisions.push(note(
                        key,
                        "transition",
                        format!("live_default {old_live} -> {}", row.live_default),
                    ));
                }
            }
        }
    }

    decisions.push(note(
        REGISTRY_FILE,
        "summary",
        format!(
            "rows={} ungated={} green={} refused={} live_defaults={}",
            rows.len(),
            counts.get("ungated").copied().unwrap_or(0),
            counts.get("green").copied().unwrap_or(0),
            counts.get("refused").copied().unwrap_or(0),
            live_defaults
        ),
    ));

    InstrumentClaimsReport {
        violations,
        decisions,
    }
}

/// Filesystem/git entry point used by `check-instrument-claims` and
/// `check-all`.
pub fn check(root: &Path) -> InstrumentClaimsReport {
    let current = match std::fs::read_to_string(root.join(REGISTRY_FILE)) {
        Ok(source) => source,
        Err(error) => {
            return InstrumentClaimsReport {
                violations: vec![violation(
                    REGISTRY_FILE,
                    format!(
                        "{REGISTRY_FILE} is unreadable ({error}); a gate that cannot read the \
                         registry refuses rather than concluding nothing is claimed"
                    ),
                )],
                decisions: Vec::new(),
            };
        }
    };
    let output = std::process::Command::new("git")
        .args(["show", &format!("HEAD:{REGISTRY_FILE}")])
        .current_dir(root)
        .output();
    let committed = match output {
        Ok(output) if output.status.success() => String::from_utf8(output.stdout).ok(),
        _ => None,
    };
    let path_exists = |path: &str| {
        let joined = root.join(path);
        joined.is_file() || joined.is_dir()
    };
    let owner_exists = |name: &str| root.join("crates").join(name).is_dir();
    let corpus_ids = std::fs::read_to_string(root.join(ACOUSTIC_MANIFEST_FILE))
        .ok()
        .map(|manifest| manifest_case_ids(&manifest));
    check_registry_sources(
        &current,
        committed.as_deref(),
        &path_exists,
        &owner_exists,
        corpus_ids.as_ref(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const NO_PATHS: fn(&str) -> bool = |_| false;
    const ALL_PATHS: fn(&str) -> bool = |_| true;
    const ALL_OWNERS: fn(&str) -> bool = |_| true;

    fn run(current: &str, committed: Option<&str>) -> InstrumentClaimsReport {
        check_registry_sources(current, committed, &ALL_PATHS, &ALL_OWNERS, None)
    }

    fn run_with_corpus(current: &str, ids: &[&str]) -> InstrumentClaimsReport {
        let ids: BTreeSet<String> = ids.iter().map(|id| (*id).to_string()).collect();
        check_registry_sources(current, None, &ALL_PATHS, &ALL_OWNERS, Some(&ids))
    }

    /// A minimal valid row with per-field overrides (the parser refuses
    /// duplicate keys, so overrides substitute rather than append). Each
    /// override is `(key, json_literal)`.
    fn row_with(overrides: &[(&str, &str)]) -> String {
        let mut fields: Vec<(&str, String)> = vec![
            ("filling", r#""wind-reed""#.to_string()),
            ("image", r#""tmm""#.to_string()),
            ("qoi", r#""peaks""#.to_string()),
            ("owner_crates", r#"["fs-duct"]"#.to_string()),
            ("exactness", r#"["X-Consist"]"#.to_string()),
            ("gate", r#""ungated""#.to_string()),
            ("live_default", r#""no""#.to_string()),
            ("determinism", r#""one-host""#.to_string()),
            ("evidence", "[]".to_string()),
            ("budget_row", "null".to_string()),
            ("corpus_refs", "[]".to_string()),
            ("notes", r#""seed""#.to_string()),
        ];
        for (key, value) in overrides {
            if let Some(slot) = fields.iter_mut().find(|(name, _)| name == key) {
                slot.1 = (*value).to_string();
            } else {
                fields.push((key, (*value).to_string()));
            }
        }
        let body = fields
            .iter()
            .map(|(key, value)| format!("\"{key}\":{value}"))
            .collect::<Vec<_>>()
            .join(",");
        format!("{{{body}}}")
    }

    fn row(overrides: &str) -> String {
        assert!(overrides.is_empty(), "use row_with for overrides");
        row_with(&[])
    }

    fn registry(rows: &[String]) -> String {
        format!(
            "{{\"schema\":\"{REGISTRY_SCHEMA}\",\"rows\":[{}]}}",
            rows.join(",")
        )
    }

    fn details(report: &InstrumentClaimsReport) -> String {
        report
            .violations
            .iter()
            .map(|violation| violation.detail.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn seed_row_round_trips_clean() {
        let report = run(&registry(&[row("")]), None);
        assert!(
            report.violations.is_empty(),
            "unexpected: {}",
            details(&report)
        );
        let summary = report
            .decisions
            .iter()
            .find(|note| note.verdict == "summary")
            .expect("summary note");
        assert!(
            summary
                .detail
                .contains("rows=1 ungated=1 green=0 refused=0")
        );
    }

    #[test]
    fn invalid_json_refuses() {
        let report = run("{not json", None);
        assert!(details(&report).contains("not valid JSON"));
    }

    #[test]
    fn wrong_schema_refuses() {
        let report = run("{\"schema\":\"other\",\"rows\":[]}", None);
        assert!(details(&report).contains("declares schema"));
    }

    #[test]
    fn missing_rows_refuses() {
        let report = run(&format!("{{\"schema\":\"{REGISTRY_SCHEMA}\"}}"), None);
        assert!(details(&report).contains("no \"rows\" array"));
    }

    #[test]
    fn duplicate_key_refuses() {
        let report = run(&registry(&[row(""), row("")]), None);
        assert!(details(&report).contains("duplicate (filling, image, qoi)"));
    }

    #[test]
    fn bad_enums_refuse() {
        let report = run(
            &registry(&[row_with(&[
                ("gate", r#""open""#),
                ("determinism", r#""vibes""#),
                ("live_default", r#""true""#),
            ])]),
            None,
        );
        let text = details(&report);
        assert!(text.contains("gate \"open\""), "{text}");
        assert!(text.contains("determinism \"vibes\""), "{text}");
        assert!(text.contains("live_default \"true\""), "{text}");
    }

    #[test]
    fn exactness_bounds_refuse() {
        let empty = run(&registry(&[row_with(&[("exactness", "[]")])]), None);
        assert!(details(&empty).contains("exactness has 0 entries"));
        let three = run(
            &registry(&[row_with(&[(
                "exactness",
                r#"["X-Exact","X-Struct","X-Est"]"#,
            )])]),
            None,
        );
        assert!(details(&three).contains("exactness has 3 entries"));
        let dup = run(
            &registry(&[row_with(&[("exactness", r#"["X-Exact","X-Exact"]"#)])]),
            None,
        );
        assert!(details(&dup).contains("listed twice"));
        let unknown = run(
            &registry(&[row_with(&[("exactness", r#"["exactish"]"#)])]),
            None,
        );
        assert!(details(&unknown).contains("is not one of"));
    }

    #[test]
    fn green_requires_evidence() {
        let report = run(&registry(&[row_with(&[("gate", r#""green""#)])]), None);
        assert!(details(&report).contains("gate=green with empty evidence"));
    }

    #[test]
    fn refused_requires_reason() {
        let report = run(
            &registry(&[row_with(&[("gate", r#""refused""#), ("notes", r#""  ""#)])]),
            None,
        );
        assert!(details(&report).contains("gate=refused without a reason"));
    }

    #[test]
    fn live_default_requires_green_and_budget() {
        let ungated = run(
            &registry(&[row_with(&[("live_default", r#""yes""#)])]),
            None,
        );
        let text = details(&ungated);
        assert!(text.contains("D25 requires gate=green"), "{text}");
        assert!(text.contains("no budget_row"), "{text}");

        let evidence = r#"[{"kind":"test","ref":"crates/fs-duct/src/lib.rs"}]"#;
        let green_no_budget = run(
            &registry(&[row_with(&[
                ("live_default", r#""yes""#),
                ("gate", r#""green""#),
                ("evidence", evidence),
            ])]),
            None,
        );
        let text = details(&green_no_budget);
        assert!(!text.contains("D25 requires gate=green"), "{text}");
        assert!(text.contains("no budget_row"), "{text}");

        let complete = run(
            &registry(&[row_with(&[
                ("live_default", r#""yes""#),
                ("gate", r#""green""#),
                ("budget_row", r#""blake3:abc""#),
                ("evidence", evidence),
            ])]),
            None,
        );
        assert!(complete.violations.is_empty(), "{}", details(&complete));
    }

    #[test]
    fn dangling_resolvable_evidence_refuses() {
        let report = check_registry_sources(
            &registry(&[row_with(&[
                ("gate", r#""green""#),
                (
                    "evidence",
                    r#"[{"kind":"test","ref":"crates/nope/tests/x.rs::t"}]"#,
                ),
            ])]),
            None,
            &NO_PATHS,
            &ALL_OWNERS,
            None,
        );
        assert!(details(&report).contains("does not resolve"));
    }

    #[test]
    fn recorded_only_evidence_is_not_resolved() {
        let report = check_registry_sources(
            &registry(&[row_with(&[
                ("gate", r#""green""#),
                (
                    "evidence",
                    r#"[{"kind":"listening","ref":"blake3:receipt"}]"#,
                ),
            ])]),
            None,
            &NO_PATHS,
            &ALL_OWNERS,
            None,
        );
        assert!(report.violations.is_empty(), "{}", details(&report));
    }

    #[test]
    fn unknown_evidence_kind_refuses() {
        let report = run(
            &registry(&[row_with(&[("evidence", r#"[{"kind":"vibe","ref":"x"}]"#)])]),
            None,
        );
        assert!(details(&report).contains("evidence kind \"vibe\""));
    }

    #[test]
    fn missing_owner_crate_refuses() {
        let report = check_registry_sources(
            &registry(&[row("")]),
            None,
            &ALL_PATHS,
            &|name| name != "fs-duct",
            None,
        );
        assert!(details(&report).contains("owner crate \"fs-duct\""));
    }

    #[test]
    fn deleting_a_row_refuses_but_transitions_note() {
        let before = registry(&[row(""), row_with(&[("image", r#""char-line""#)])]);
        let after_deleted = registry(&[row("")]);
        let report = run(&after_deleted, Some(&before));
        assert!(details(&report).contains("row deleted"));
        assert!(details(&report).contains("D21"));

        let after_refused = registry(&[
            row(""),
            row_with(&[
                ("image", r#""char-line""#),
                ("gate", r#""refused""#),
                ("notes", r#""lost the bake-off; kept per D21""#),
            ]),
        ]);
        let report = run(&after_refused, Some(&before));
        assert!(report.violations.is_empty(), "{}", details(&report));
        assert!(
            report
                .decisions
                .iter()
                .any(|note| note.verdict == "transition"
                    && note.detail.contains("ungated -> refused")),
            "expected a transition note"
        );
    }

    #[test]
    fn unparsable_predecessor_constrains_nothing() {
        let report = run(&registry(&[row("")]), Some("{broken"));
        assert!(report.violations.is_empty(), "{}", details(&report));
    }

    #[test]
    fn budget_row_type_refuses() {
        let report = run(&registry(&[row_with(&[("budget_row", "42")])]), None);
        assert!(details(&report).contains("budget_row must be null or"));
    }

    #[test]
    fn determinism_above_the_default_ceiling_refuses() {
        // Zero cross-ISA goldens exist in the music stack, so the default
        // ceiling is one-host: both higher ladder classes refuse, naming
        // the capping crate and the promotion path.
        for claimed in ["same-isa", "cross-isa"] {
            let report = run(
                &registry(&[row_with(&[("determinism", &format!("\"{claimed}\""))])]),
                None,
            );
            let text = details(&report);
            assert!(text.contains("exceeds owner crate \"fs-duct\""), "{text}");
            assert!(text.contains("3ez8g.13.4"), "{text}");
        }
    }

    #[test]
    fn tribo_cap_names_the_capping_crate() {
        let report = run(
            &registry(&[row_with(&[
                ("owner_crates", r#"["fs-tribo","fs-couple"]"#),
                ("determinism", r#""cross-isa""#),
            ])]),
            None,
        );
        let text = details(&report);
        assert!(text.contains("owner crate \"fs-tribo\""), "{text}");
    }

    #[test]
    fn off_ladder_declarations_are_exempt_but_surfaced() {
        for declared in ["statistical", "fast-mode"] {
            let report = run(
                &registry(&[row_with(&[("determinism", &format!("\"{declared}\""))])]),
                None,
            );
            assert!(report.violations.is_empty(), "{}", details(&report));
            assert!(
                report.decisions.iter().any(|note| note.verdict == "review"
                    && note.detail.contains("off the replay ladder")),
                "missing review note for {declared}"
            );
        }
    }

    #[test]
    fn one_host_rows_pass_the_ceiling() {
        let report = run(&registry(&[row("")]), None);
        assert!(report.violations.is_empty(), "{}", details(&report));
    }

    #[test]
    fn corpus_refs_resolve_or_refuse() {
        // Registered id: passes. Unregistered: refuses by name. Missing
        // prefix: refuses (no parallel namespace). Unreadable manifest with
        // any ref: refuses (a gate that cannot see the corpus registry
        // refuses rather than concluding the citation is fine).
        let cited = registry(&[row_with(&[(
            "corpus_refs",
            r#"["vvreg:acoustic-ernoult-2021-xxxx"]"#,
        )])]);
        let good = run_with_corpus(&cited, &["acoustic-ernoult-2021-xxxx"]);
        assert!(good.violations.is_empty(), "{}", details(&good));

        let dangling = run_with_corpus(&cited, &["some-other-case"]);
        assert!(
            details(&dangling).contains("names no case"),
            "{}",
            details(&dangling)
        );

        let unprefixed = run_with_corpus(
            &registry(&[row_with(&[("corpus_refs", r#"["ernoult-2021"]"#)])]),
            &["acoustic-ernoult-2021-xxxx"],
        );
        assert!(
            details(&unprefixed).contains("lacks the vvreg: prefix"),
            "{}",
            details(&unprefixed)
        );

        let unreadable = run(&cited, None);
        assert!(
            details(&unreadable).contains("is unreadable"),
            "{}",
            details(&unreadable)
        );
    }

    #[test]
    fn manifest_case_ids_parse_the_tsv_shape() {
        let manifest = "schema_version\tcase_id\tfamily\n1\tacoustic-a\tx\n1\tacoustic-b\ty\n";
        let ids = manifest_case_ids(manifest);
        assert!(ids.contains("acoustic-a") && ids.contains("acoustic-b"));
        assert_eq!(ids.len(), 2);
    }
}
