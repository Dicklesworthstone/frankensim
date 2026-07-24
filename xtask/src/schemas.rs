//! Public schema freeze governance (bead
//! `frankensim-extreal-program-f85xj.16.5`).
//!
//! `schema-policy.json` names the small set of serialized schemas FrankenSim
//! promises not to break, and classifies every other public version constant on
//! the product boundary as internal. This checker refuses three failure modes
//! that prose alone cannot catch:
//!
//! 1. **Drift** — a frozen schema's version constant moves without the policy
//!    record moving with it, so the published promise describes a format that no
//!    longer exists.
//! 2. **Unbacked obligation** — a record declares a migration obligation whose
//!    evidence test is missing, so "we migrate old documents" is unproven.
//! 3. **Accretion** — a new public version constant appears in a product-boundary
//!    crate and is neither promised nor explicitly disclaimed, drifting toward
//!    accidental-public status.
//!
//! The version values, lockstep values, and evidence tests are read out of the
//! actual sources, so this check is a statement about the tree rather than about
//! the registry's self-description.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

pub const CHECK: &str = "schema-policy";
pub const POLICY_FILE: &str = "schema-policy.json";
const POLICY_SCHEMA: &str = "frankensim-schema-policy-v1";
const DOCTRINE_FILE: &str = "docs/SCHEMA_POLICY.md";
const POLICY_BEAD: &str = "frankensim-extreal-program-f85xj.16.5";
const MAX_POLICY_BYTES: usize = 1024 * 1024;

const ROOT_FIELDS: [&str; 7] = [
    "accretion_scope",
    "doctrine",
    "frozen",
    "not_promised",
    "policy_bead",
    "ratified",
    "schema",
];
const FROZEN_FIELDS: [&str; 14] = [
    "current_version",
    "deprecation_horizon",
    "id",
    "lockstep",
    "major_bump",
    "migration_evidence",
    "migration_obligation",
    "minor_bump",
    "owner",
    "surface",
    "title",
    "version_constant",
    "version_kind",
    "version_location",
];
const LOCKSTEP_FIELDS: [&str; 3] = ["constant", "location", "value"];
const EVIDENCE_FIELDS: [&str; 2] = ["location", "test"];
const NOT_PROMISED_FIELDS: [&str; 4] = ["constant", "location", "reason", "reason_class"];

const INTEGER_KIND: &str = "integer-const";
const STRING_KIND: &str = "string-const";
const VERSION_KINDS: [&str; 2] = [INTEGER_KIND, STRING_KIND];

const AUTO_MIGRATION: &str = "auto-migration-receipt";
const REFUSE_UNMIGRATABLE: &str = "refuse-unmigratable";
const NO_PREDECESSOR: &str = "no-predecessor";
const OBLIGATIONS: [&str; 3] = [AUTO_MIGRATION, REFUSE_UNMIGRATABLE, NO_PREDECESSOR];

const REASON_CLASSES: [&str; 5] = [
    "derived-lockstep",
    "identity-domain",
    "internal-receipt",
    "internal-row",
    "internal-wire",
];

/// Doctrine terms `docs/SCHEMA_POLICY.md` must carry, so the published promise
/// cannot quietly lose the machinery that makes it enforceable.
const REQUIRED_DOCTRINE: [&str; 8] = [
    "schema-policy.json",
    "check-schemas",
    AUTO_MIGRATION,
    REFUSE_UNMIGRATABLE,
    NO_PREDECESSOR,
    "accretion",
    "deprecation horizon",
    "internal and breakable",
];

pub struct SchemaPolicyReport {
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

/// Require an exact field set: a missing field leaves a promise undefined, and
/// an unexpected one is an unreviewed extension of the policy vocabulary.
fn exact_fields(
    map: &BTreeMap<String, JsonValue>,
    expected: &[&str],
    entity: &str,
    violations: &mut Vec<Violation>,
) -> bool {
    let found: BTreeSet<&str> = map.keys().map(String::as_str).collect();
    let want: BTreeSet<&str> = expected.iter().copied().collect();
    let missing: Vec<&str> = want.difference(&found).copied().collect();
    let extra: Vec<&str> = found.difference(&want).copied().collect();
    if !missing.is_empty() {
        violations.push(violation(
            entity,
            format!("missing required field(s): {}", missing.join(", ")),
        ));
    }
    if !extra.is_empty() {
        violations.push(violation(
            entity,
            format!("unexpected field(s): {}", extra.join(", ")),
        ));
    }
    missing.is_empty() && extra.is_empty()
}

fn nonempty<'a>(
    map: &'a BTreeMap<String, JsonValue>,
    key: &str,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> Option<&'a str> {
    let value = map.get(key).and_then(text).unwrap_or("");
    if value.trim().is_empty() {
        violations.push(violation(
            entity,
            format!("missing non-empty string `{key}`"),
        ));
        None
    } else {
        Some(value)
    }
}

fn one_of(
    value: &str,
    allowed: &[&str],
    key: &str,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> bool {
    if allowed.contains(&value) {
        return true;
    }
    violations.push(violation(
        entity,
        format!(
            "`{key}` is {value:?}; expected one of {}",
            allowed.join(", ")
        ),
    ));
    false
}

fn is_date(value: &str) -> bool {
    value.len() == 10
        && value.as_bytes()[4] == b'-'
        && value.as_bytes()[7] == b'-'
        && value
            .bytes()
            .enumerate()
            .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit())
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ConstKind {
    Integer,
    Str,
}

impl ConstKind {
    fn declared(self) -> &'static str {
        match self {
            ConstKind::Integer => INTEGER_KIND,
            ConstKind::Str => STRING_KIND,
        }
    }
}

/// Parse a single `const NAME: TY = LITERAL;` declaration.
///
/// Deliberately line-oriented and literal-only: a version constant that is
/// computed, aliased, or spread across lines is not a frozen-schema constant
/// this policy can bind, and returning `None` makes that refuse loudly at the
/// call site rather than resolve to a guess.
fn parse_const_decl(line: &str) -> Option<(&str, ConstKind, String)> {
    let trimmed = line.trim();
    let rest = trimmed
        .strip_prefix("pub const ")
        .or_else(|| trimmed.strip_prefix("pub(crate) const "))
        .or_else(|| trimmed.strip_prefix("pub(super) const "))
        .or_else(|| trimmed.strip_prefix("const "))?;
    let (name, rest) = rest.split_once(':')?;
    let name = name.trim();
    if name.is_empty()
        || !name
            .bytes()
            .all(|byte| byte.is_ascii_uppercase() || byte.is_ascii_digit() || byte == b'_')
    {
        return None;
    }
    let (ty, literal) = rest.split_once('=')?;
    let literal = literal.trim().strip_suffix(';')?.trim();
    match ty.trim() {
        "u8" | "u16" | "u32" | "u64" | "i64" | "usize" => {
            let digits = literal.replace('_', "");
            (!digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())).then_some((
                name,
                ConstKind::Integer,
                digits,
            ))
        }
        "&str" | "&'static str" => {
            let value = literal.strip_prefix('"')?.strip_suffix('"')?;
            (!value.contains('\\')).then(|| (name, ConstKind::Str, value.to_string()))
        }
        _ => None,
    }
}

fn resolve_const(source: &str, name: &str) -> Option<(ConstKind, String)> {
    source.lines().find_map(|line| {
        let (found, kind, value) = parse_const_decl(line)?;
        (found == name).then_some((kind, value))
    })
}

/// Public integer version constants declared in one source file. This is the
/// accretion surface: `pub` because a private constant is not reachable by a
/// downstream consumer, integer-typed because that is the shape a serialized
/// format version takes in this workspace.
fn scan_public_version_constants(source: &str) -> Vec<(String, String)> {
    source
        .lines()
        .filter_map(|line| {
            if !line.trim_start().starts_with("pub const ") {
                return None;
            }
            let (name, kind, value) = parse_const_decl(line)?;
            (kind == ConstKind::Integer && name.contains("VERSION"))
                .then(|| (name.to_string(), value))
        })
        .collect()
}

/// Where a constant is classified. A constant must be classified exactly once:
/// two homes means two different promises for one symbol.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Classification {
    Frozen(String),
    Lockstep(String),
    NotPromised,
}

impl Classification {
    fn describe(&self) -> String {
        match self {
            Classification::Frozen(id) => format!("frozen schema {id}"),
            Classification::Lockstep(id) => format!("lockstep of frozen schema {id}"),
            Classification::NotPromised => "not_promised".to_string(),
        }
    }
}

/// Resolve one declared constant against the tree and compare it with the
/// registry's recorded value.
fn verify_constant(
    sources: &BTreeMap<String, String>,
    location: &str,
    constant: &str,
    expected: &str,
    expected_kind: Option<ConstKind>,
    entity: &str,
    violations: &mut Vec<Violation>,
) {
    let Some(source) = sources.get(location) else {
        violations.push(violation(
            entity,
            format!("declared location {location:?} is not a readable tracked source"),
        ));
        return;
    };
    let Some((kind, value)) = resolve_const(source, constant) else {
        violations.push(violation(
            entity,
            format!("constant `{constant}` is not declared as a literal const in {location}"),
        ));
        return;
    };
    if let Some(expected_kind) = expected_kind
        && kind != expected_kind
    {
        violations.push(violation(
            entity,
            format!(
                "`{constant}` is declared {} but the record says {}",
                kind.declared(),
                expected_kind.declared()
            ),
        ));
        return;
    }
    if value != expected {
        violations.push(violation(
            entity,
            format!(
                "{location} declares `{constant}` = {value:?} but the policy record says {expected:?}; \
                 bump the record in the same commit as the schema"
            ),
        ));
    }
}

fn audit_frozen(
    record: &BTreeMap<String, JsonValue>,
    index: usize,
    sources: &BTreeMap<String, String>,
    classified: &mut BTreeMap<String, Classification>,
    ids: &mut BTreeSet<String>,
    violations: &mut Vec<Violation>,
    decisions: &mut Vec<PolicyNote>,
) {
    let fallback = format!("frozen[{index}]");
    let entity = record
        .get("id")
        .and_then(text)
        .filter(|id| !id.trim().is_empty())
        .unwrap_or(&fallback)
        .to_string();
    if !exact_fields(record, &FROZEN_FIELDS, &entity, violations) {
        return;
    }
    let Some(id) = nonempty(record, "id", &entity, violations) else {
        return;
    };
    if !ids.insert(id.to_string()) {
        violations.push(violation(&entity, format!("duplicate frozen id {id:?}")));
    }
    for key in [
        "title",
        "owner",
        "minor_bump",
        "major_bump",
        "deprecation_horizon",
    ] {
        nonempty(record, key, &entity, violations);
    }

    let surface = record.get("surface").and_then(arr).unwrap_or(&[]);
    if surface.is_empty() {
        violations.push(violation(
            &entity,
            "`surface` must state at least one promised element; an unstated surface cannot be honoured",
        ));
    }
    for item in surface {
        if text(item).is_none_or(|value| value.trim().is_empty()) {
            violations.push(violation(
                &entity,
                "`surface` entries must be non-empty strings",
            ));
        }
    }

    let kind = nonempty(record, "version_kind", &entity, violations).unwrap_or("");
    let kind_ok = one_of(kind, &VERSION_KINDS, "version_kind", &entity, violations);
    let expected_kind = kind_ok.then(|| {
        if kind == INTEGER_KIND {
            ConstKind::Integer
        } else {
            ConstKind::Str
        }
    });

    let constant = nonempty(record, "version_constant", &entity, violations);
    let location = nonempty(record, "version_location", &entity, violations);
    let current = nonempty(record, "current_version", &entity, violations);
    if let (Some(constant), Some(location), Some(current)) = (constant, location, current) {
        verify_constant(
            sources,
            location,
            constant,
            current,
            expected_kind,
            &entity,
            violations,
        );
        match classified.insert(constant.to_string(), Classification::Frozen(id.to_string())) {
            None => {}
            Some(previous) => violations.push(violation(
                &entity,
                format!(
                    "`{constant}` is already classified as {}; a constant must be classified exactly once",
                    previous.describe()
                ),
            )),
        }
    }

    for (position, entry) in record
        .get("lockstep")
        .and_then(arr)
        .unwrap_or(&[])
        .iter()
        .enumerate()
    {
        let lock_entity = format!("{entity}#lockstep[{position}]");
        let Some(entry) = obj(entry) else {
            violations.push(violation(&lock_entity, "lockstep entries must be objects"));
            continue;
        };
        if !exact_fields(entry, &LOCKSTEP_FIELDS, &lock_entity, violations) {
            continue;
        }
        let constant = nonempty(entry, "constant", &lock_entity, violations);
        let location = nonempty(entry, "location", &lock_entity, violations);
        let value = nonempty(entry, "value", &lock_entity, violations);
        if let (Some(constant), Some(location), Some(value)) = (constant, location, value) {
            verify_constant(
                sources,
                location,
                constant,
                value,
                None,
                &lock_entity,
                violations,
            );
            if let Some(previous) = classified.insert(
                constant.to_string(),
                Classification::Lockstep(id.to_string()),
            ) {
                violations.push(violation(
                    &lock_entity,
                    format!(
                        "`{constant}` is already classified as {}; a constant must be classified exactly once",
                        previous.describe()
                    ),
                ));
            }
        }
    }

    let obligation = nonempty(record, "migration_obligation", &entity, violations).unwrap_or("");
    one_of(
        obligation,
        &OBLIGATIONS,
        "migration_obligation",
        &entity,
        violations,
    );
    let evidence = record
        .get("migration_evidence")
        .and_then(arr)
        .unwrap_or(&[]);
    for (position, entry) in evidence.iter().enumerate() {
        let evidence_entity = format!("{entity}#migration_evidence[{position}]");
        let Some(entry) = obj(entry) else {
            violations.push(violation(
                &evidence_entity,
                "migration_evidence entries must be objects",
            ));
            continue;
        };
        if !exact_fields(entry, &EVIDENCE_FIELDS, &evidence_entity, violations) {
            continue;
        }
        let test = nonempty(entry, "test", &evidence_entity, violations);
        let location = nonempty(entry, "location", &evidence_entity, violations);
        if let (Some(test), Some(location)) = (test, location) {
            match sources.get(location) {
                None => violations.push(violation(
                    &evidence_entity,
                    format!("evidence location {location:?} is not a readable tracked source"),
                )),
                Some(source) => {
                    if !source.contains(&format!("fn {test}(")) {
                        violations.push(violation(
                            &evidence_entity,
                            format!("{location} does not declare evidence test `{test}`"),
                        ));
                    }
                }
            }
        }
    }

    // The obligation is what makes the promise enforceable, so each kind has a
    // distinct, checkable consequence.
    match obligation {
        AUTO_MIGRATION | REFUSE_UNMIGRATABLE => {
            if evidence.is_empty() {
                violations.push(violation(
                    &entity,
                    format!(
                        "migration_obligation {obligation:?} requires at least one named evidence test"
                    ),
                ));
            }
        }
        NO_PREDECESSOR => {
            if !evidence.is_empty() {
                violations.push(violation(
                    &entity,
                    "migration_obligation \"no-predecessor\" must not cite evidence; there is no predecessor to migrate",
                ));
            }
            // A no-predecessor schema that has moved past its first version has
            // silently skipped the obligation to declare how the predecessor is
            // handled. Refusing here is the whole point of the classification.
            let current = record.get("current_version").and_then(text).unwrap_or("");
            let first_version = current == "1" || current.ends_with("-v1");
            if !first_version {
                violations.push(violation(
                    &entity,
                    format!(
                        "migration_obligation \"no-predecessor\" but current_version is {current:?}; \
                         a bumped schema must declare auto-migration-receipt or refuse-unmigratable \
                         and cite its evidence"
                    ),
                ));
            }
        }
        _ => {}
    }

    decisions.push(note(
        &entity,
        "frozen",
        format!(
            "{} at {} = {} ({}; {} evidence test(s))",
            record.get("owner").and_then(text).unwrap_or("?"),
            record.get("version_constant").and_then(text).unwrap_or("?"),
            record.get("current_version").and_then(text).unwrap_or("?"),
            obligation,
            evidence.len()
        ),
    ));
}

fn audit_not_promised(
    entries: &[JsonValue],
    sources: &BTreeMap<String, String>,
    classified: &mut BTreeMap<String, Classification>,
    violations: &mut Vec<Violation>,
) {
    let mut order: Vec<(String, String)> = Vec::new();
    for (index, entry) in entries.iter().enumerate() {
        let entity = format!("not_promised[{index}]");
        let Some(entry) = obj(entry) else {
            violations.push(violation(&entity, "not_promised entries must be objects"));
            continue;
        };
        if !exact_fields(entry, &NOT_PROMISED_FIELDS, &entity, violations) {
            continue;
        }
        let constant = nonempty(entry, "constant", &entity, violations);
        let location = nonempty(entry, "location", &entity, violations);
        nonempty(entry, "reason", &entity, violations);
        let class = nonempty(entry, "reason_class", &entity, violations).unwrap_or("");
        one_of(class, &REASON_CLASSES, "reason_class", &entity, violations);
        let (Some(constant), Some(location)) = (constant, location) else {
            continue;
        };
        order.push((location.to_string(), constant.to_string()));
        // A disclaimer for a constant that no longer exists is registry rot: it
        // makes the list look complete while covering nothing.
        match sources.get(location) {
            None => violations.push(violation(
                &entity,
                format!("declared location {location:?} is not a readable tracked source"),
            )),
            Some(source) => {
                if resolve_const(source, constant).is_none() {
                    violations.push(violation(
                        &entity,
                        format!(
                            "{location} no longer declares `{constant}`; remove the stale not_promised row"
                        ),
                    ));
                }
            }
        }
        if let Some(previous) = classified.insert(constant.to_string(), Classification::NotPromised)
        {
            violations.push(violation(
                &entity,
                format!(
                    "`{constant}` is already classified as {}; a constant must be classified exactly once",
                    previous.describe()
                ),
            ));
        }
    }
    let mut sorted = order.clone();
    sorted.sort();
    if order != sorted {
        violations.push(violation(
            POLICY_FILE,
            "`not_promised` must be sorted by (location, constant) so additions are reviewable as a diff",
        ));
    }
}

/// Every public integer version constant in a product-boundary crate must be
/// either promised or explicitly disclaimed.
fn audit_accretion(
    scope: &[String],
    sources: &BTreeMap<String, String>,
    classified: &BTreeMap<String, Classification>,
    violations: &mut Vec<Violation>,
) -> usize {
    let mut scanned = 0usize;
    for crate_name in scope {
        let prefix = format!("crates/{crate_name}/src/");
        let mut seen_any = false;
        for (path, source) in sources.range(prefix.clone()..) {
            if !path.starts_with(&prefix) {
                break;
            }
            seen_any = true;
            for (constant, value) in scan_public_version_constants(source) {
                scanned += 1;
                if !classified.contains_key(&constant) {
                    violations.push(violation(
                        path,
                        format!(
                            "public version constant `{constant}` = {value} is not classified in {POLICY_FILE}; \
                             add it to `frozen` (a promise) or `not_promised` (internal and breakable)"
                        ),
                    ));
                }
            }
        }
        if !seen_any {
            violations.push(violation(
                POLICY_FILE,
                format!(
                    "accretion_scope names {crate_name:?} but no source was read under {prefix}"
                ),
            ));
        }
    }
    scanned
}

fn check_doctrine(sources: &BTreeMap<String, String>, violations: &mut Vec<Violation>) {
    let Some(source) = sources.get(DOCTRINE_FILE) else {
        violations.push(violation(DOCTRINE_FILE, "document is unreadable"));
        return;
    };
    for required in REQUIRED_DOCTRINE {
        if !source.contains(required) {
            violations.push(violation(
                DOCTRINE_FILE,
                format!("schema-freeze doctrine is missing {required:?}"),
            ));
        }
    }
}

pub fn check_sources(policy: &str, sources: &BTreeMap<String, String>) -> SchemaPolicyReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();

    let parsed = match JsonParser::with_string_limit(policy, MAX_POLICY_BYTES).finish() {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(POLICY_FILE, format!("invalid JSON: {error}")));
            return SchemaPolicyReport {
                violations,
                decisions,
            };
        }
    };
    let Some(root) = obj(&parsed) else {
        violations.push(violation(
            POLICY_FILE,
            "document root must be a JSON object",
        ));
        return SchemaPolicyReport {
            violations,
            decisions,
        };
    };
    if !exact_fields(root, &ROOT_FIELDS, POLICY_FILE, &mut violations) {
        return SchemaPolicyReport {
            violations,
            decisions,
        };
    }
    if root.get("schema").and_then(text) != Some(POLICY_SCHEMA) {
        violations.push(violation(
            POLICY_FILE,
            format!("`schema` must be {POLICY_SCHEMA:?}"),
        ));
    }
    if root.get("policy_bead").and_then(text) != Some(POLICY_BEAD) {
        violations.push(violation(
            POLICY_FILE,
            format!("`policy_bead` must be {POLICY_BEAD:?}"),
        ));
    }
    if root.get("doctrine").and_then(text) != Some(DOCTRINE_FILE) {
        violations.push(violation(
            POLICY_FILE,
            format!("`doctrine` must be {DOCTRINE_FILE:?}"),
        ));
    }
    let ratified = root.get("ratified").and_then(text).unwrap_or("");
    if !is_date(ratified) {
        violations.push(violation(
            POLICY_FILE,
            format!("`ratified` must be an ISO date, found {ratified:?}"),
        ));
    }

    let scope: Vec<String> = root
        .get("accretion_scope")
        .and_then(arr)
        .unwrap_or(&[])
        .iter()
        .filter_map(text)
        .map(str::to_string)
        .collect();
    if scope.is_empty() {
        violations.push(violation(
            POLICY_FILE,
            "`accretion_scope` must name at least one product-boundary crate",
        ));
    }
    let mut sorted_scope = scope.clone();
    sorted_scope.sort();
    sorted_scope.dedup();
    if sorted_scope != scope {
        violations.push(violation(
            POLICY_FILE,
            "`accretion_scope` must be sorted and free of duplicates",
        ));
    }

    let mut classified: BTreeMap<String, Classification> = BTreeMap::new();
    let mut ids = BTreeSet::new();
    let frozen = root.get("frozen").and_then(arr).unwrap_or(&[]);
    if frozen.is_empty() {
        violations.push(violation(POLICY_FILE, "`frozen` must not be empty"));
    }
    for (index, entry) in frozen.iter().enumerate() {
        match obj(entry) {
            Some(record) => audit_frozen(
                record,
                index,
                sources,
                &mut classified,
                &mut ids,
                &mut violations,
                &mut decisions,
            ),
            None => violations.push(violation(
                &format!("frozen[{index}]"),
                "frozen entries must be objects",
            )),
        }
    }

    let not_promised = root.get("not_promised").and_then(arr).unwrap_or(&[]);
    audit_not_promised(not_promised, sources, &mut classified, &mut violations);

    let scanned = audit_accretion(&scope, sources, &classified, &mut violations);
    check_doctrine(sources, &mut violations);

    decisions.push(note(
        "repository",
        "inventory",
        format!(
            "{} frozen schema(s) promised, {} constant(s) explicitly not promised, \
             {} public version constant(s) scanned across {} product-boundary crate(s)",
            frozen.len(),
            not_promised.len(),
            scanned,
            scope.len()
        ),
    ));

    SchemaPolicyReport {
        violations,
        decisions,
    }
}

/// Collect every source the policy binds: the declared locations plus the whole
/// accretion surface.
fn gather_sources(root: &Path, policy: &str) -> BTreeMap<String, String> {
    let mut wanted: BTreeSet<String> = BTreeSet::new();
    wanted.insert(DOCTRINE_FILE.to_string());
    // Declared locations are discovered by a light scan of the registry text so
    // that an unparsable or malformed registry still produces field-level
    // violations rather than a cascade of "unreadable source" noise.
    let mut remainder = policy;
    while let Some(start) = remainder.find("\"location\"") {
        remainder = &remainder[start + "\"location\"".len()..];
        let Some(open) = remainder.find('"') else {
            break;
        };
        let after = &remainder[open + 1..];
        let Some(close) = after.find('"') else { break };
        wanted.insert(after[..close].to_string());
        remainder = &after[close..];
    }
    let mut remainder = policy;
    while let Some(start) = remainder.find("\"version_location\"") {
        remainder = &remainder[start + "\"version_location\"".len()..];
        let Some(open) = remainder.find('"') else {
            break;
        };
        let after = &remainder[open + 1..];
        let Some(close) = after.find('"') else { break };
        wanted.insert(after[..close].to_string());
        remainder = &after[close..];
    }

    let mut sources = BTreeMap::new();
    for relative in wanted {
        if let Ok(contents) = std::fs::read_to_string(root.join(&relative)) {
            sources.insert(relative, contents);
        }
    }

    for crate_dir in scope_crates(policy) {
        let base = root.join("crates").join(&crate_dir).join("src");
        let mut stack = vec![base];
        while let Some(dir) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&dir) else {
                continue;
            };
            let mut paths: Vec<PathBuf> = entries.flatten().map(|entry| entry.path()).collect();
            paths.sort();
            for path in paths {
                if path.is_dir() {
                    stack.push(path);
                } else if path.extension().is_some_and(|ext| ext == "rs")
                    && let Ok(contents) = std::fs::read_to_string(&path)
                    && let Ok(relative) = path.strip_prefix(root)
                {
                    sources.insert(relative.to_string_lossy().replace('\\', "/"), contents);
                }
            }
        }
    }
    sources
}

fn scope_crates(policy: &str) -> Vec<String> {
    let Ok(parsed) = JsonParser::with_string_limit(policy, MAX_POLICY_BYTES).finish() else {
        return Vec::new();
    };
    obj(&parsed)
        .and_then(|root| root.get("accretion_scope"))
        .and_then(arr)
        .unwrap_or(&[])
        .iter()
        .filter_map(text)
        .map(str::to_string)
        .collect()
}

pub fn check_schema_policy(root: &Path) -> SchemaPolicyReport {
    let policy = match std::fs::read_to_string(root.join(POLICY_FILE)) {
        Ok(policy) => policy,
        Err(error) => {
            return SchemaPolicyReport {
                violations: vec![violation(
                    POLICY_FILE,
                    format!("file is unreadable: {error}"),
                )],
                decisions: Vec::new(),
            };
        }
    };
    let sources = gather_sources(root, &policy);
    check_sources(&policy, &sources)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sources(pairs: &[(&str, &str)]) -> BTreeMap<String, String> {
        let mut map: BTreeMap<String, String> = pairs
            .iter()
            .map(|(path, body)| (path.to_string(), body.to_string()))
            .collect();
        map.entry(DOCTRINE_FILE.to_string()).or_insert_with(|| {
            REQUIRED_DOCTRINE
                .iter()
                .map(|term| format!("{term}\n"))
                .collect()
        });
        map
    }

    fn policy(frozen: &str, not_promised: &str) -> String {
        format!(
            r#"{{
              "schema": "{POLICY_SCHEMA}",
              "policy_bead": "{POLICY_BEAD}",
              "ratified": "2026-07-24",
              "doctrine": "{DOCTRINE_FILE}",
              "accretion_scope": ["fs-demo"],
              "frozen": [{frozen}],
              "not_promised": [{not_promised}]
            }}"#
        )
    }

    const DEMO_FROZEN: &str = r#"{
        "id": "demo.format",
        "title": "Demo",
        "owner": "fs-demo",
        "version_kind": "integer-const",
        "version_constant": "DEMO_FORMAT_VERSION",
        "version_location": "crates/fs-demo/src/lib.rs",
        "current_version": "3",
        "lockstep": [],
        "surface": ["the demo bytes"],
        "minor_bump": "additive",
        "major_bump": "anything else",
        "migration_obligation": "auto-migration-receipt",
        "migration_evidence": [
          {"test": "demo_v2_migrates", "location": "crates/fs-demo/tests/demo.rs"}
        ],
        "deprecation_horizon": "one major"
      }"#;

    const DEMO_SRC: &str = "pub const DEMO_FORMAT_VERSION: u32 = 3;\n";
    const DEMO_TEST: &str = "#[test]\nfn demo_v2_migrates() { }\n";

    fn clean_sources() -> BTreeMap<String, String> {
        sources(&[
            ("crates/fs-demo/src/lib.rs", DEMO_SRC),
            ("crates/fs-demo/tests/demo.rs", DEMO_TEST),
        ])
    }

    #[test]
    fn a_clean_policy_passes_and_reports_an_inventory() {
        let report = check_sources(&policy(DEMO_FROZEN, ""), &clean_sources());
        assert!(
            report.violations.is_empty(),
            "expected no violations: {:?}",
            report.violations
        );
        assert!(report.decisions.iter().any(|note| note.verdict == "frozen"));
        assert!(
            report
                .decisions
                .iter()
                .any(|note| note.verdict == "inventory")
        );
    }

    #[test]
    fn a_version_bump_without_a_policy_bump_is_caught() {
        let bumped = sources(&[
            (
                "crates/fs-demo/src/lib.rs",
                "pub const DEMO_FORMAT_VERSION: u32 = 4;\n",
            ),
            ("crates/fs-demo/tests/demo.rs", DEMO_TEST),
        ]);
        let report = check_sources(&policy(DEMO_FROZEN, ""), &bumped);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("but the policy record says")),
            "expected drift violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_missing_evidence_test_is_caught() {
        let missing = sources(&[
            ("crates/fs-demo/src/lib.rs", DEMO_SRC),
            ("crates/fs-demo/tests/demo.rs", "fn unrelated() {}\n"),
        ]);
        let report = check_sources(&policy(DEMO_FROZEN, ""), &missing);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("does not declare evidence test")),
            "expected evidence violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn an_unclassified_public_version_constant_is_flagged_as_accretion() {
        let accreted = sources(&[
            (
                "crates/fs-demo/src/lib.rs",
                "pub const DEMO_FORMAT_VERSION: u32 = 3;\npub const SNEAKY_WIRE_VERSION: u32 = 1;\n",
            ),
            ("crates/fs-demo/tests/demo.rs", DEMO_TEST),
        ]);
        let report = check_sources(&policy(DEMO_FROZEN, ""), &accreted);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("SNEAKY_WIRE_VERSION")
                    && item.detail.contains("is not classified")),
            "expected accretion violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_private_or_non_version_constant_is_not_accretion() {
        let quiet = sources(&[
            (
                "crates/fs-demo/src/lib.rs",
                "pub const DEMO_FORMAT_VERSION: u32 = 3;\nconst PRIVATE_WIRE_VERSION: u32 = 1;\npub const DEMO_LIMIT: u32 = 5;\n",
            ),
            ("crates/fs-demo/tests/demo.rs", DEMO_TEST),
        ]);
        let report = check_sources(&policy(DEMO_FROZEN, ""), &quiet);
        assert!(
            report.violations.is_empty(),
            "private and non-version constants are not a public promise: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_stale_not_promised_row_is_caught() {
        let stale = r#"{
            "constant": "GONE_VERSION",
            "location": "crates/fs-demo/src/lib.rs",
            "reason_class": "internal-wire",
            "reason": "removed last week"
          }"#;
        let report = check_sources(&policy(DEMO_FROZEN, stale), &clean_sources());
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("remove the stale not_promised row")),
            "expected stale-row violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_constant_cannot_be_classified_twice() {
        let duplicate = r#"{
            "constant": "DEMO_FORMAT_VERSION",
            "location": "crates/fs-demo/src/lib.rs",
            "reason_class": "internal-wire",
            "reason": "also claimed as internal"
          }"#;
        let report = check_sources(&policy(DEMO_FROZEN, duplicate), &clean_sources());
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("classified exactly once")),
            "expected double-classification violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn a_bumped_no_predecessor_schema_must_declare_a_migration_path() {
        let frozen = DEMO_FROZEN
            .replace("\"auto-migration-receipt\"", "\"no-predecessor\"")
            .replace(
                r#"[
          {"test": "demo_v2_migrates", "location": "crates/fs-demo/tests/demo.rs"}
        ]"#,
                "[]",
            );
        let report = check_sources(&policy(&frozen, ""), &clean_sources());
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("must declare auto-migration-receipt")),
            "a no-predecessor schema at version 3 has skipped its obligation: {:?}",
            report.violations
        );
    }

    #[test]
    fn an_unknown_policy_field_is_refused() {
        let extended = policy(DEMO_FROZEN, "").replace(
            "\"accretion_scope\"",
            "\"surprise\": 1, \"accretion_scope\"",
        );
        let report = check_sources(&extended, &clean_sources());
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("unexpected field")),
            "expected unexpected-field violation: {:?}",
            report.violations
        );
    }

    #[test]
    fn const_parsing_handles_the_shapes_this_repo_actually_uses() {
        assert_eq!(
            parse_const_decl("pub const FSIM_VERSION: u32 = 1;").map(|(n, k, v)| (n, k, v)),
            Some(("FSIM_VERSION", ConstKind::Integer, "1".to_string()))
        );
        assert_eq!(
            parse_const_decl("pub const SCHEMA_VERSION: i64 = 20;").map(|(_, _, v)| v),
            Some("20".to_string())
        );
        assert_eq!(
            parse_const_decl(r#"const SCHEMA: &str = "frankensim-source-manifest-v1";"#)
                .map(|(_, k, v)| (k, v)),
            Some((ConstKind::Str, "frankensim-source-manifest-v1".to_string()))
        );
        // A computed or aliased constant is deliberately unresolvable.
        assert!(parse_const_decl("pub const DERIVED: u32 = OTHER + 1;").is_none());
        assert!(parse_const_decl("// pub const COMMENTED: u32 = 1;").is_none());
    }

    /// Drive the real entry point — registry read, source gathering, and audit —
    /// over a fabricated tree, so the wiring between `gather_sources` and
    /// `check_sources` is covered rather than only the pure core.
    #[test]
    fn the_audit_runs_end_to_end_over_a_fabricated_tree() {
        let base = std::env::temp_dir().join(format!("fsim-schema-policy-{}", std::process::id()));
        let write = |relative: &str, body: &str| {
            let path = base.join(relative);
            std::fs::create_dir_all(path.parent().expect("parent")).expect("create dirs");
            std::fs::write(&path, body).expect("write fixture");
        };
        write("crates/fs-demo/src/lib.rs", DEMO_SRC);
        write("crates/fs-demo/tests/demo.rs", DEMO_TEST);
        write(
            DOCTRINE_FILE,
            &REQUIRED_DOCTRINE
                .iter()
                .map(|term| format!("{term}\n"))
                .collect::<String>(),
        );
        write(POLICY_FILE, &policy(DEMO_FROZEN, ""));

        let clean = check_schema_policy(&base);
        assert!(
            clean.violations.is_empty(),
            "fabricated clean tree must pass: {:?}",
            clean.violations
        );
        assert!(
            clean
                .decisions
                .iter()
                .any(|note| note.verdict == "inventory"),
            "an inventory note is always emitted"
        );

        // Accrete a new public serialized-format version in the same crate; the
        // audit must refuse it without any registry edit.
        write(
            "crates/fs-demo/src/wire.rs",
            "pub const LATE_ARRIVAL_VERSION: u32 = 1;\n",
        );
        let accreted = check_schema_policy(&base);
        assert!(
            accreted.violations.iter().any(|item| {
                item.crate_name == "crates/fs-demo/src/wire.rs"
                    && item.detail.contains("LATE_ARRIVAL_VERSION")
            }),
            "a newly added public version constant must be flagged: {:?}",
            accreted.violations
        );

        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn the_live_schema_policy_is_clean() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let report = check_schema_policy(root);
        assert!(
            report.violations.is_empty(),
            "live schema-policy.json must be clean: {:?}",
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
