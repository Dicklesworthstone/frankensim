//! Live moonshot-portfolio governance (bead
//! `frankensim-extreal-program-f85xj.16.3`).
//!
//! `moonshot-portfolio.json` records the status-quo cap and complete
//! declarations for active `[M]` work. Beads remains authoritative for labels,
//! status, and dependency edges; `vertical-capability-graph.json` names the
//! product critical path. This checker joins those sources, calls the pure
//! `fs-govern` policy algebra, and independently rejects any active moonshot
//! that transitively blocks a critical-path Bead.

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation};
use fs_govern::{
    DisplacementRecord, MOONSHOT_V1_INITIAL_CAP, MoonshotBudget, MoonshotDeclaration,
    MoonshotDisposition, MoonshotPortfolio, NamedFalsifier, QuarterlyReview, ReplacementAdmission,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::Path;

pub const CHECK: &str = "moonshot-portfolio";
pub const PORTFOLIO_FILE: &str = "moonshot-portfolio.json";
const GRAPH_FILE: &str = "vertical-capability-graph.json";
const BEADS_FILE: &str = ".beads/issues.jsonl";
const PORTFOLIO_SCHEMA: &str = "frankensim-moonshot-portfolio-v1";
const POLICY_BEAD: &str = "frankensim-extreal-program-f85xj.16.3";
const ACTIVE_STATUS: &str = "in_progress";
const MAX_BEAD_ROW_BYTES: usize = 1024 * 1024;
const EXACT_MOONSHOT_LABELS: [&str; 2] = ["ambition-m", "ambition:M"];
const MOONSHOT_LABEL_PREFIX: &str = "moonshot";
const DECLARED_MOONSHOT_LABELS: [&str; 3] = ["ambition-m", "ambition:M", "moonshot*"];
const INITIAL_MOONSHOT_BEADS: [&str; 6] = [
    "frankensim-epic-ascent-7tv.22",
    "frankensim-epic-bedrock-6ys.20",
    "frankensim-leapfrog-2026-program-i94v.2.3.2.2",
    "frankensim-leapfrog-2026-program-i94v.5.1.2.1",
    "frankensim-leapfrog-2026-program-i94v.5.2.11.1",
    "frankensim-leapfrog-2026-program-i94v.5.4.1.2",
];
const _: [(); MOONSHOT_V1_INITIAL_CAP as usize] = [(); INITIAL_MOONSHOT_BEADS.len()];

pub struct MoonshotPolicyReport {
    pub violations: Vec<Violation>,
    pub decisions: Vec<PolicyNote>,
}

fn date_key(value: &str) -> Option<(u16, u8, u8)> {
    if value.len() != 10
        || value.as_bytes().get(4) != Some(&b'-')
        || value.as_bytes().get(7) != Some(&b'-')
        || !value
            .bytes()
            .enumerate()
            .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit())
    {
        return None;
    }
    let year = value.get(..4)?.parse().ok()?;
    let month: u8 = value.get(5..7)?.parse().ok()?;
    let day: u8 = value.get(8..)?.parse().ok()?;
    let max_day = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if year % 400 == 0 || (year % 4 == 0 && year % 100 != 0) => 29,
        2 => 28,
        _ => return None,
    };
    (day != 0 && day <= max_day).then_some((year, month, day))
}

fn is_date(value: &str) -> bool {
    date_key(value).is_some()
}

#[derive(Clone, Debug)]
struct Issue {
    status: String,
    assignee: Option<String>,
    labels: BTreeSet<String>,
    blocking_dependencies: Vec<String>,
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

fn number(value: &JsonValue) -> Option<u64> {
    match value {
        JsonValue::Number(value) => value.parse().ok(),
        _ => None,
    }
}

fn field<'a>(map: &'a BTreeMap<String, JsonValue>, key: &str) -> Option<&'a JsonValue> {
    map.get(key)
}

fn nonempty<'a>(
    map: &'a BTreeMap<String, JsonValue>,
    key: &str,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> Option<&'a str> {
    let value = field(map, key).and_then(text).unwrap_or("");
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

fn parse_root(source: &str, entity: &str, violations: &mut Vec<Violation>) -> Option<JsonValue> {
    match JsonParser::with_string_limit(source, MAX_BEAD_ROW_BYTES).finish() {
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

fn string_set(value: Option<&JsonValue>) -> Option<BTreeSet<String>> {
    arr(value?)?
        .iter()
        .map(|item| text(item).map(str::to_string))
        .collect()
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
        let blocking_dependencies = field(row, "dependencies")
            .and_then(arr)
            .unwrap_or_default()
            .iter()
            .filter_map(obj)
            .filter(|dependency| field(dependency, "type").and_then(text) == Some("blocks"))
            .filter_map(|dependency| field(dependency, "depends_on_id").and_then(text))
            .map(str::to_string)
            .collect();
        if issues
            .insert(
                id.to_string(),
                Issue {
                    status: status.to_string(),
                    assignee: field(row, "assignee").and_then(text).map(str::to_string),
                    labels,
                    blocking_dependencies,
                },
            )
            .is_some()
        {
            violations.push(violation(BEADS_FILE, format!("duplicate issue id {id:?}")));
        }
    }
    issues
}

fn is_moonshot(issue: &Issue) -> bool {
    EXACT_MOONSHOT_LABELS
        .iter()
        .any(|label| issue.labels.contains(*label))
        || issue
            .labels
            .iter()
            .any(|label| label.starts_with(MOONSHOT_LABEL_PREFIX))
}

fn parse_falsifier(
    root: &BTreeMap<String, JsonValue>,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> Option<NamedFalsifier> {
    let nested = field(root, "falsifier").and_then(obj);
    let Some(nested) = nested else {
        violations.push(violation(entity, "missing object `falsifier`"));
        return None;
    };
    let (Some(id), Some(observation), Some(decision_rule)) = (
        nonempty(nested, "id", entity, violations),
        nonempty(nested, "observation", entity, violations),
        nonempty(nested, "decision_rule", entity, violations),
    ) else {
        return None;
    };
    match NamedFalsifier::new(id, observation, decision_rule) {
        Ok(value) => Some(value),
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            None
        }
    }
}

fn parse_budget(
    root: &BTreeMap<String, JsonValue>,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> Option<MoonshotBudget> {
    let Some(nested) = field(root, "budget").and_then(obj) else {
        violations.push(violation(entity, "missing object `budget`"));
        return None;
    };
    let effort = field(nested, "effort_minutes")
        .and_then(number)
        .unwrap_or(0);
    let Some(deadline) = nonempty(nested, "calendar_deadline", entity, violations) else {
        return None;
    };
    match MoonshotBudget::new(effort, deadline) {
        Ok(value) => Some(value),
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            None
        }
    }
}

fn parse_review(
    root: &BTreeMap<String, JsonValue>,
    entity: &str,
    violations: &mut Vec<Violation>,
) -> Option<QuarterlyReview> {
    let Some(nested) = field(root, "quarterly_review").and_then(obj) else {
        violations.push(violation(entity, "missing object `quarterly_review`"));
        return None;
    };
    let (Some(next_review), Some(status), Some(required_evidence)) = (
        nonempty(nested, "next_review", entity, violations),
        nonempty(nested, "status", entity, violations),
        nonempty(nested, "required_evidence", entity, violations),
    ) else {
        return None;
    };
    match QuarterlyReview::new(next_review, status, required_evidence) {
        Ok(value) => Some(value),
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            None
        }
    }
}

fn parse_declaration(
    root: &BTreeMap<String, JsonValue>,
    entity: &str,
    expected_legacy_baseline: Option<bool>,
    violations: &mut Vec<Violation>,
) -> Option<MoonshotDeclaration> {
    let legacy_baseline = match field(root, "admission_kind").and_then(text) {
        Some("legacy-baseline") => true,
        Some("replacement") => false,
        _ => {
            violations.push(violation(
                entity,
                "admission_kind must be `legacy-baseline` or `replacement`",
            ));
            return None;
        }
    };
    if let Some(expected_legacy_baseline) = expected_legacy_baseline
        && legacy_baseline != expected_legacy_baseline
    {
        let expected_kind = if expected_legacy_baseline {
            "legacy-baseline"
        } else {
            "replacement"
        };
        violations.push(violation(
            entity,
            format!("admission_kind must be {expected_kind:?}"),
        ));
    }
    let (Some(bead_id), Some(lane), Some(owner), Some(disjointness)) = (
        nonempty(root, "bead_id", entity, violations),
        nonempty(root, "lane", entity, violations),
        nonempty(root, "owner", entity, violations),
        nonempty(root, "critical_path_disjointness", entity, violations),
    ) else {
        return None;
    };
    let falsifier = parse_falsifier(root, entity, violations)?;
    let budget = parse_budget(root, entity, violations)?;
    let review = parse_review(root, entity, violations)?;
    match MoonshotDeclaration::new(
        bead_id,
        lane,
        owner,
        legacy_baseline,
        falsifier,
        budget,
        disjointness,
        review,
    ) {
        Ok(value) => Some(value),
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            None
        }
    }
}

fn replay_portfolio_history(
    root: &BTreeMap<String, JsonValue>,
    issues: &BTreeMap<String, Issue>,
    violations: &mut Vec<Violation>,
) -> BTreeSet<String> {
    let expected_initial: BTreeSet<String> = INITIAL_MOONSHOT_BEADS
        .iter()
        .map(|id| (*id).to_string())
        .collect();
    let declared_initial = string_set(field(root, "initial_active_beads")).unwrap_or_default();
    if declared_initial != expected_initial {
        violations.push(violation(
            PORTFOLIO_FILE,
            "initial_active_beads must retain the immutable six-Bead v1 baseline",
        ));
    }
    let mut active = expected_initial;
    let mut terminal = BTreeSet::new();
    let history = field(root, "disposition_history")
        .and_then(arr)
        .unwrap_or_default();
    for (index, row) in history.iter().enumerate() {
        let entity = format!("disposition_history[{index}]");
        let Some(row) = obj(row) else {
            violations.push(violation(&entity, "history row must be an object"));
            continue;
        };
        let (
            Some(action),
            Some(displaced),
            Some(disposition),
            Some(_state_artifact),
            Some(_reason),
            Some(recorded_at),
        ) = (
            nonempty(row, "action", &entity, violations),
            nonempty(row, "displaced_bead_id", &entity, violations),
            nonempty(row, "disposition", &entity, violations),
            nonempty(row, "state_artifact", &entity, violations),
            nonempty(row, "reason", &entity, violations),
            nonempty(row, "recorded_at", &entity, violations),
        )
        else {
            continue;
        };
        if parse_disposition(disposition).is_none() {
            violations.push(violation(&entity, "unknown terminal disposition"));
        }
        if !recorded_at.get(..10).is_some_and(is_date) {
            violations.push(violation(
                &entity,
                "recorded_at must begin with a valid YYYY-MM-DD date",
            ));
        }
        if !active.remove(displaced) {
            violations.push(violation(
                &entity,
                format!("displaced Bead {displaced:?} is not active at this history step"),
            ));
        }
        terminal.insert(displaced.to_string());
        match action {
            "liquidate" => {
                if field(row, "admitted_bead_id").is_some() {
                    violations.push(violation(
                        &entity,
                        "liquidation must not carry admitted_bead_id",
                    ));
                }
            }
            "replace" => {
                let Some(admitted) = nonempty(row, "admitted_bead_id", &entity, violations) else {
                    continue;
                };
                if terminal.contains(admitted) {
                    violations.push(violation(
                        &entity,
                        format!("terminal Bead {admitted:?} cannot be revived as a replacement"),
                    ));
                }
                if !active.insert(admitted.to_string()) {
                    violations.push(violation(
                        &entity,
                        format!("replacement Bead {admitted:?} is already active"),
                    ));
                }
                match issues.get(admitted) {
                    None => violations.push(violation(
                        &entity,
                        format!("replacement Bead {admitted:?} does not exist"),
                    )),
                    Some(issue) if !is_moonshot(issue) => violations.push(violation(
                        &entity,
                        format!("replacement Bead {admitted:?} lacks a moonshot label"),
                    )),
                    Some(_) => {}
                }
            }
            _ => violations.push(violation(
                &entity,
                "action must be `replace` or `liquidate`",
            )),
        }
    }
    active
}

fn critical_path_ids(graph: &str, violations: &mut Vec<Violation>) -> BTreeSet<String> {
    let Some(parsed) = parse_root(graph, GRAPH_FILE, violations) else {
        return BTreeSet::new();
    };
    let root = obj(&parsed).expect("parse_root returns object");
    let mut ids = string_set(field(root, "linked_gates")).unwrap_or_default();
    for (array_name, bead_field) in [
        ("capability_bindings", "implementing_beads"),
        ("seam_owners", "beads"),
    ] {
        let Some(rows) = field(root, array_name).and_then(arr) else {
            violations.push(violation(
                GRAPH_FILE,
                format!("missing array `{array_name}` for moonshot reachability"),
            ));
            continue;
        };
        for (index, row) in rows.iter().enumerate() {
            let Some(row) = obj(row) else {
                violations.push(violation(
                    GRAPH_FILE,
                    format!("{array_name}[{index}] must be an object"),
                ));
                continue;
            };
            let Some(beads) = string_set(field(row, bead_field)) else {
                violations.push(violation(
                    GRAPH_FILE,
                    format!("{array_name}[{index}].{bead_field} must be a string array"),
                ));
                continue;
            };
            ids.extend(beads);
        }
    }
    ids
}

fn blocking_path(
    start: &str,
    critical: &BTreeSet<String>,
    issues: &BTreeMap<String, Issue>,
) -> Option<Vec<String>> {
    if critical.contains(start) {
        return Some(vec![start.to_string()]);
    }
    let mut dependents: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for (issue_id, issue) in issues {
        for dependency in &issue.blocking_dependencies {
            dependents.entry(dependency).or_default().push(issue_id);
        }
    }
    let mut queue = VecDeque::from([(start.to_string(), vec![start.to_string()])]);
    let mut seen = BTreeSet::from([start.to_string()]);
    while let Some((current, path)) = queue.pop_front() {
        for dependent in dependents.get(current.as_str()).into_iter().flatten() {
            if !seen.insert((*dependent).to_string()) {
                continue;
            }
            let mut next_path = path.clone();
            next_path.push((*dependent).to_string());
            if critical.contains(*dependent) {
                return Some(next_path);
            }
            queue.push_back(((*dependent).to_string(), next_path));
        }
    }
    None
}

fn parse_disposition(value: &str) -> Option<MoonshotDisposition> {
    match value {
        "completed" => Some(MoonshotDisposition::Completed),
        "falsified" => Some(MoonshotDisposition::Falsified),
        "shelved-with-state" => Some(MoonshotDisposition::ShelvedWithState),
        _ => None,
    }
}

fn check_rehearsal(
    root: &BTreeMap<String, JsonValue>,
    portfolio: MoonshotPortfolio,
    issues: &BTreeMap<String, Issue>,
    critical: &BTreeSet<String>,
    violations: &mut Vec<Violation>,
) {
    let entity = "admission_rehearsal";
    let Some(rehearsal) = field(root, entity).and_then(obj) else {
        violations.push(violation(PORTFOLIO_FILE, "missing admission rehearsal"));
        return;
    };
    if field(rehearsal, "kind").and_then(text) != Some("hypothetical-only") {
        violations.push(violation(entity, "kind must be `hypothetical-only`"));
    }
    if field(rehearsal, "expected_verdict").and_then(text) != Some("admissible-replacement") {
        violations.push(violation(
            entity,
            "expected_verdict must be `admissible-replacement`",
        ));
    }
    let transcript = field(rehearsal, "transcript")
        .and_then(arr)
        .unwrap_or_default();
    if transcript.len() < 7 || transcript.iter().any(|line| text(line).is_none()) {
        violations.push(violation(
            entity,
            "transcript must retain all seven string decision steps",
        ));
    }
    let Some(candidate) = field(rehearsal, "candidate").and_then(obj) else {
        violations.push(violation(entity, "missing object `candidate`"));
        return;
    };
    let Some(candidate) = parse_declaration(
        candidate,
        "admission_rehearsal.candidate",
        Some(false),
        violations,
    ) else {
        return;
    };
    match issues.get(candidate.bead_id()) {
        None => violations.push(violation(
            entity,
            "candidate must name a real open moonshot Bead so graph disjointness is checkable",
        )),
        Some(issue) if !is_moonshot(issue) => violations.push(violation(
            entity,
            "candidate Bead does not carry a recognized moonshot label",
        )),
        Some(issue) if issue.status != "open" => violations.push(violation(
            entity,
            format!(
                "candidate Bead must be unadmitted open backlog, found {}",
                issue.status
            ),
        )),
        Some(issue) => {
            if let Some(assignee) = &issue.assignee
                && assignee != candidate.owner()
            {
                violations.push(violation(
                    entity,
                    format!(
                        "candidate owner {:?} differs from live Bead assignee {assignee:?}",
                        candidate.owner()
                    ),
                ));
            }
            if let Some(path) = blocking_path(candidate.bead_id(), critical, issues) {
                violations.push(violation(
                    entity,
                    format!(
                        "replacement candidate blocks the product critical path through {}",
                        path.join(" -> ")
                    ),
                ));
            }
        }
    }
    let Some(displaces) = field(rehearsal, "displaces").and_then(obj) else {
        violations.push(violation(entity, "missing object `displaces`"));
        return;
    };
    let (Some(bead_id), Some(disposition), Some(state_artifact), Some(reason)) = (
        nonempty(displaces, "bead_id", entity, violations),
        nonempty(displaces, "disposition", entity, violations),
        nonempty(displaces, "state_artifact", entity, violations),
        nonempty(displaces, "reason", entity, violations),
    ) else {
        return;
    };
    let Some(disposition) = parse_disposition(disposition) else {
        violations.push(violation(entity, "unknown displacement disposition"));
        return;
    };
    let displacement = match DisplacementRecord::new(bead_id, disposition, state_artifact, reason) {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            return;
        }
    };
    let admission = match ReplacementAdmission::new(candidate, displacement) {
        Ok(value) => value,
        Err(error) => {
            violations.push(violation(entity, error.to_string()));
            return;
        }
    };
    if let Err(error) = portfolio.assess_replacement(&admission) {
        violations.push(violation(entity, error.to_string()));
    }
}

fn check_sources(portfolio_source: &str, graph: &str, beads: &str) -> MoonshotPolicyReport {
    let mut violations = Vec::new();
    let mut decisions = Vec::new();
    let issues = parse_issues(beads, &mut violations);
    let critical = critical_path_ids(graph, &mut violations);
    let Some(parsed) = parse_root(portfolio_source, PORTFOLIO_FILE, &mut violations) else {
        return MoonshotPolicyReport {
            violations,
            decisions,
        };
    };
    let root = obj(&parsed).expect("parse_root returns object");
    for (key, expected) in [
        ("schema", PORTFOLIO_SCHEMA),
        ("policy_bead", POLICY_BEAD),
        ("issue_authority", BEADS_FILE),
        ("critical_path_graph", GRAPH_FILE),
        ("active_status", ACTIVE_STATUS),
    ] {
        if field(root, key).and_then(text) != Some(expected) {
            violations.push(violation(
                PORTFOLIO_FILE,
                format!("{key} must be {expected:?}"),
            ));
        }
    }
    if field(root, "policy_version").and_then(number) != Some(1) {
        violations.push(violation(PORTFOLIO_FILE, "policy_version must be 1"));
    }
    for key in ["captured_at", "cap_rule", "no_claim"] {
        let _ = nonempty(root, key, PORTFOLIO_FILE, &mut violations);
    }
    let labels = string_set(field(root, "moonshot_labels")).unwrap_or_default();
    let expected_labels = DECLARED_MOONSHOT_LABELS
        .iter()
        .map(|label| (*label).to_string())
        .collect();
    if labels != expected_labels {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!("moonshot_labels must be exactly {expected_labels:?}"),
        ));
    }

    let cap = field(root, "status_quo_cap")
        .and_then(number)
        .and_then(|value| u32::try_from(value).ok())
        .unwrap_or(0);
    if cap > MOONSHOT_V1_INITIAL_CAP {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!(
                "status_quo_cap {cap} exceeds the immutable v1 initial ceiling {MOONSHOT_V1_INITIAL_CAP}"
            ),
        ));
    }
    let tagged: Vec<_> = issues
        .iter()
        .filter(|(_, issue)| is_moonshot(issue))
        .collect();
    let active_ids: BTreeSet<_> = tagged
        .iter()
        .filter(|(_, issue)| issue.status == ACTIVE_STATUS)
        .map(|(id, _)| (*id).clone())
        .collect();
    let replayed_active = replay_portfolio_history(root, &issues, &mut violations);
    if replayed_active != active_ids {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!(
                "disposition_history replays to {replayed_active:?}, but live active moonshots are {active_ids:?}"
            ),
        ));
    }
    if usize::try_from(cap).ok() != Some(active_ids.len()) {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!(
                "status_quo_cap must equal the live active count so completed or shelved work cannot leave a refillable slot: cap {cap}, active {}",
                active_ids.len()
            ),
        ));
    }
    let open_count = tagged
        .iter()
        .filter(|(_, issue)| issue.status == "open")
        .count();
    let closed_count = tagged
        .iter()
        .filter(|(_, issue)| issue.status == "closed")
        .count();
    let known_count = active_ids.len() + open_count + closed_count;
    if known_count != tagged.len() {
        violations.push(violation(
            BEADS_FILE,
            "moonshot inventory contains a status outside open/in_progress/closed",
        ));
    }
    if let Some(snapshot) = field(root, "inventory_snapshot").and_then(obj) {
        for (key, actual) in [
            ("tagged_total", tagged.len()),
            ("active", active_ids.len()),
            ("open_backlog", open_count),
            ("closed_terminal", closed_count),
        ] {
            let expected = field(snapshot, key).and_then(number);
            if expected != u64::try_from(actual).ok() {
                violations.push(violation(
                    "inventory_snapshot",
                    format!("{key} is stale: recorded {expected:?}, live {actual}"),
                ));
            }
        }
        let _ = nonempty(
            snapshot,
            "inactive_rule",
            "inventory_snapshot",
            &mut violations,
        );
    } else {
        violations.push(violation(
            PORTFOLIO_FILE,
            "missing object `inventory_snapshot`",
        ));
    }

    let mut declarations = Vec::new();
    let mut declared_ids = BTreeSet::new();
    let rows = field(root, "active_lanes")
        .and_then(arr)
        .unwrap_or_default();
    for (index, row) in rows.iter().enumerate() {
        let entity = format!("active_lanes[{index}]");
        let Some(row) = obj(row) else {
            violations.push(violation(&entity, "declaration must be an object"));
            continue;
        };
        let Some(declaration) = parse_declaration(row, &entity, None, &mut violations) else {
            continue;
        };
        let bead_id = declaration.bead_id().to_string();
        let should_be_legacy = INITIAL_MOONSHOT_BEADS.contains(&bead_id.as_str());
        if declaration.is_legacy_baseline() != should_be_legacy {
            violations.push(violation(
                &entity,
                if should_be_legacy {
                    "an undisplaced v1 baseline Bead must retain admission_kind `legacy-baseline`"
                } else {
                    "a post-baseline active Bead must carry admission_kind `replacement` and a history row"
                },
            ));
        }
        if !declared_ids.insert(bead_id.clone()) {
            violations.push(violation(
                &entity,
                format!("duplicate declaration {bead_id:?}"),
            ));
        }
        match issues.get(&bead_id) {
            None => violations.push(violation(&entity, format!("unknown Bead {bead_id:?}"))),
            Some(issue) if !is_moonshot(issue) => violations.push(violation(
                &entity,
                format!("Bead {bead_id:?} lacks a moonshot label"),
            )),
            Some(issue) if issue.status != ACTIVE_STATUS => violations.push(violation(
                &entity,
                format!("Bead {bead_id:?} is {}, not active", issue.status),
            )),
            Some(issue) => {
                if let Some(assignee) = &issue.assignee
                    && assignee != declaration.owner()
                {
                    violations.push(violation(
                        &entity,
                        format!(
                            "declared owner {:?} differs from live Bead assignee {assignee:?}",
                            declaration.owner()
                        ),
                    ));
                }
            }
        }
        if let Some(path) = blocking_path(&bead_id, &critical, &issues) {
            violations.push(violation(
                &entity,
                format!(
                    "active moonshot blocks the product critical path through {}",
                    path.join(" -> ")
                ),
            ));
        }
        declarations.push(declaration);
    }
    for missing in active_ids.difference(&declared_ids) {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!("active moonshot Bead {missing:?} lacks a declaration"),
        ));
    }
    for extra in declared_ids.difference(&active_ids) {
        violations.push(violation(
            PORTFOLIO_FILE,
            format!("declared moonshot Bead {extra:?} is not active"),
        ));
    }

    let portfolio = match MoonshotPortfolio::new(cap, declarations) {
        Ok(portfolio) => Some(portfolio),
        Err(error) => {
            violations.push(violation(PORTFOLIO_FILE, error.to_string()));
            None
        }
    };
    if let Some(portfolio) = portfolio {
        check_rehearsal(root, portfolio, &issues, &critical, &mut violations);
    }
    decisions.push(note(
        "<repo>",
        "inventory",
        format!(
            "{} tagged moonshot Beads: {} active of cap {}, {} open backlog, {} closed; {} critical-path Bead ids checked transitively",
            tagged.len(),
            active_ids.len(),
            cap,
            open_count,
            closed_count,
            critical.len()
        ),
    ));
    decisions.push(note(
        "admission_rehearsal",
        "hypothetical-only",
        "replacement assessment retains a named displacement, terminal disposition, state artifact, falsifier, effort/calendar budget, and path-disjointness declaration; no Bead state changed",
    ));
    MoonshotPolicyReport {
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
        "moonshot-portfolio.json",
        "check-moonshots",
        "disposition history",
        "shelved-with-state",
        "quarterly falsifier review",
        "vertical-capability-graph.json",
    ] {
        if !source.contains(required) {
            violations.push(violation(
                "docs/CONVENTIONS.md",
                format!("moonshot doctrine is missing {required:?}"),
            ));
        }
    }
}

pub fn check_moonshot_policy(root: &Path) -> MoonshotPolicyReport {
    let mut read_violations = Vec::new();
    let read = |relative: &str, violations: &mut Vec<Violation>| {
        std::fs::read_to_string(root.join(relative)).map_err(|error| {
            violations.push(violation(relative, format!("file is unreadable: {error}")));
        })
    };
    let portfolio = read(PORTFOLIO_FILE, &mut read_violations);
    let graph = read(GRAPH_FILE, &mut read_violations);
    let beads = read(BEADS_FILE, &mut read_violations);
    let (Ok(portfolio), Ok(graph), Ok(beads)) = (portfolio, graph, beads) else {
        return MoonshotPolicyReport {
            violations: read_violations,
            decisions: Vec::new(),
        };
    };
    let mut report = check_sources(&portfolio, &graph, &beads);
    report.violations.extend(read_violations);
    check_documentation(root, &mut report.violations);
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transitive_moonshot_blocking_path_is_caught() {
        let issues = BTreeMap::from([
            (
                "moonshot".to_string(),
                Issue {
                    status: ACTIVE_STATUS.to_string(),
                    assignee: None,
                    labels: BTreeSet::from(["ambition:M".to_string()]),
                    blocking_dependencies: Vec::new(),
                },
            ),
            (
                "bridge".to_string(),
                Issue {
                    status: "open".to_string(),
                    assignee: None,
                    labels: BTreeSet::new(),
                    blocking_dependencies: vec!["moonshot".to_string()],
                },
            ),
            (
                "critical".to_string(),
                Issue {
                    status: "open".to_string(),
                    assignee: None,
                    labels: BTreeSet::new(),
                    blocking_dependencies: vec!["bridge".to_string()],
                },
            ),
        ]);
        let path = blocking_path(
            "moonshot",
            &BTreeSet::from(["critical".to_string()]),
            &issues,
        )
        .expect("seeded blocking path");
        assert_eq!(path, ["moonshot", "bridge", "critical"]);
    }

    #[test]
    fn live_moonshot_portfolio_is_clean() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let report = check_moonshot_policy(root);
        assert!(
            report.violations.is_empty(),
            "live moonshot portfolio must be clean: {:?}",
            report.violations
        );
        assert!(
            report
                .decisions
                .iter()
                .any(|item| item.verdict == "hypothetical-only")
        );
    }

    #[test]
    fn seeded_inventory_and_rehearsal_faults_fail_closed() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let portfolio = std::fs::read_to_string(root.join(PORTFOLIO_FILE)).expect("portfolio");
        let graph = std::fs::read_to_string(root.join(GRAPH_FILE)).expect("graph");
        let beads = std::fs::read_to_string(root.join(BEADS_FILE)).expect("beads");

        let stale_count = portfolio.replacen("\"active\": 6", "\"active\": 5", 1);
        let report = check_sources(&stale_count, &graph, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("active is stale")),
            "stale inventory count must fail: {:?}",
            report.violations
        );

        let unrecorded_live_mismatch = portfolio.replacen(
            "\"disposition_history\": []",
            "\"disposition_history\": [{\"action\":\"liquidate\",\"displaced_bead_id\":\"frankensim-epic-ascent-7tv.22\",\"disposition\":\"completed\",\"state_artifact\":\"artifact:test\",\"reason\":\"seeded history mutation\",\"recorded_at\":\"2026-07-22T00:00:00Z\"}]",
            1,
        );
        let report = check_sources(&unrecorded_live_mismatch, &graph, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("replays to")),
            "history and live Beads state must agree: {:?}",
            report.violations
        );

        let inflated_cap = portfolio.replacen("\"status_quo_cap\": 6", "\"status_quo_cap\": 7", 1);
        let report = check_sources(&inflated_cap, &graph, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("initial ceiling")),
            "the initial cap must never grow: {:?}",
            report.violations
        );

        let untracked_candidate = portfolio.replacen(
            "frankensim-epic-bedrock-6ys.5.1",
            "rehearsal:untracked-candidate",
            1,
        );
        let report = check_sources(&untracked_candidate, &graph, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("real open moonshot Bead")),
            "a rehearsal without live graph authority must fail: {:?}",
            report.violations
        );

        let missing_state = portfolio.replacen(
            "bead:frankensim-leapfrog-2026-program-i94v.5.2.11.1@2026-07-22T21:56:15Z",
            "",
            1,
        );
        let report = check_sources(&missing_state, &graph, &beads);
        assert!(
            report
                .violations
                .iter()
                .any(|item| item.detail.contains("state_artifact")),
            "state-less shelving must fail: {:?}",
            report.violations
        );
    }
}
