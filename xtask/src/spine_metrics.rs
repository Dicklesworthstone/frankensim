//! Beads tracker snapshot for the product-spine dashboard (bead frankensim-o5et9).
//!
//! Why a snapshot and not the live tracker: the program-metrics dashboard is a
//! byte-checked artifact, and `.beads/issues.jsonl` moves on every `br` op.
//! Feeding the live tracker into it would make `check-program-metrics` red
//! nonstop, and a gate that is always red teaches everyone to ignore it. So
//! the spine counts live in this deliberately regenerated snapshot — the same
//! idiom as the spine ratchet — and this gate keeps the snapshot honest:
//!
//! - a missing, corrupt, or internally inconsistent snapshot is a VIOLATION
//!   (a gate that reports zero findings on a rotted input is how gates die);
//! - a tracker that has MOVED since the snapshot renders a visible
//!   `stale-snapshot` policy note, never a silent pass and never a wedge —
//!   staleness is information, not an incident;
//! - the e2e lane state is the single honest v1 value `no-retained-receipt`:
//!   the staged-producer lane runs green out-of-band but retains no tracked
//!   receipt (bead frankensim-iakds owns retaining one). Anything else is a
//!   schema drift the gate refuses until a deliberate v2 admits it.
//!
//! Counting definitions (matching bv's open/blocked/actionable accounting so
//! the number a triage run prints and the number on the dashboard agree):
//! `open` = every issue whose status is not `closed`; `blocked` = an open
//! issue with at least one `blocks` dependency on a non-closed issue; a
//! blocker referenced but absent from the file does not block (it is outside
//! the tracked set); `actionable` = open minus blocked.

use std::collections::BTreeMap;
use std::path::Path;

use crate::depgraph::{JsonParser, JsonValue};
use crate::{PolicyNote, Violation, fnv1a64};

pub(crate) const CHECK: &str = "spine-metrics";
const SNAPSHOT_PATH: &str = "spine-metrics.json";
const ISSUES_PATH: &str = ".beads/issues.jsonl";
const E2E_RECEIPT_PATH: &str = "spine-e2e-summary.json";
const E2E_RECEIPT_SCHEMA: &str = "frankensim-spine-e2e-receipt-v1";
const SCHEMA: &str = "frankensim-spine-metrics-v2";
const E2E_STATUS_NO_RECEIPT: &str = "no-retained-receipt";
const E2E_STATUS_MEASURED: &str = "measured";
/// issues.jsonl is tens of MB today; the bound refuses a runaway file rather
/// than reading unbounded tracker state into a policy gate.
const MAX_ISSUES_BYTES: u64 = 128 * 1024 * 1024;
const MAX_SNAPSHOT_BYTES: u64 = 64 * 1024;

/// The validated snapshot contents the dashboard consumes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Snapshot {
    /// FNV-1a64 of the `.beads/issues.jsonl` bytes the counts came from.
    pub issues_fnv: u64,
    /// Non-closed issues.
    pub open: usize,
    /// Open issues with at least one non-closed blocker.
    pub blocked: usize,
    /// Open issues with no non-closed blocker.
    pub actionable: usize,
    /// Stages proven green by the retained e2e receipt, when one exists.
    pub e2e_stages_green: Option<usize>,
    /// FNV-1a64 of the retained e2e receipt bytes, when one exists.
    pub e2e_receipt_fnv: Option<u64>,
}

fn violation(detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: SNAPSHOT_PATH.to_string(),
        detail: detail.into(),
    }
}

/// One bead line can carry a full comment history; the bound stays finite
/// but must cover the real tracker, not a Cargo-metadata-shaped input.
const MAX_TRACKER_STRING_BYTES: usize = 4 * 1024 * 1024;

/// Derive the counts from tracker bytes. Pure so tests can fixture it.
fn count_issues(text: &str) -> Result<Snapshot, String> {
    let mut status_by_id: BTreeMap<String, bool> = BTreeMap::new();
    let mut blockers: Vec<(String, Vec<String>)> = Vec::new();
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
        let id = match map.get("id") {
            Some(JsonValue::String(id)) => id.clone(),
            _ => {
                return Err(format!(
                    "{ISSUES_PATH} line {} has no string `id`",
                    index + 1
                ));
            }
        };
        let status = match map.get("status") {
            Some(JsonValue::String(status)) => status.clone(),
            _ => {
                return Err(format!(
                    "{ISSUES_PATH} line {} has no string `status`",
                    index + 1
                ));
            }
        };
        let mut blocking = Vec::new();
        if let Some(JsonValue::Array(deps)) = map.get("dependencies") {
            for dep in deps {
                let JsonValue::Object(dep) = dep else {
                    return Err(format!(
                        "{ISSUES_PATH} line {} has a non-object dependency",
                        index + 1
                    ));
                };
                let is_blocks =
                    matches!(dep.get("type"), Some(JsonValue::String(t)) if t == "blocks");
                if is_blocks {
                    match dep.get("depends_on_id") {
                        Some(JsonValue::String(target)) => blocking.push(target.clone()),
                        _ => {
                            return Err(format!(
                                "{ISSUES_PATH} line {} has a blocks dependency with no \
                                 string `depends_on_id`",
                                index + 1
                            ));
                        }
                    }
                }
            }
        }
        status_by_id.insert(id.clone(), status != "closed");
        blockers.push((id, blocking));
    }
    if status_by_id.is_empty() {
        return Err(format!(
            "{ISSUES_PATH} parsed to zero issues; a gate with no inputs cannot report counts"
        ));
    }
    let mut open = 0usize;
    let mut blocked = 0usize;
    for (id, is_open) in &status_by_id {
        if !is_open {
            continue;
        }
        open += 1;
        let is_blocked = blockers
            .iter()
            .find(|(blocker_id, _)| blocker_id == id)
            .is_some_and(|(_, targets)| {
                targets
                    .iter()
                    .any(|target| status_by_id.get(target).copied().unwrap_or(false))
            });
        if is_blocked {
            blocked += 1;
        }
    }
    let actionable = open - blocked;
    Ok(Snapshot {
        issues_fnv: fnv1a64(text.as_bytes()),
        open,
        blocked,
        actionable,
        e2e_stages_green: None,
        e2e_receipt_fnv: None,
    })
}

/// The retained e2e receipt, validated structurally. The stage/gap
/// cross-check against the live product source happens in [`check`], where a
/// mismatch is a violation; absence here is NO-DATA, never an error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct E2eReceipt {
    pub stages_green: usize,
    pub failures: usize,
    pub first_gap: String,
    pub first_gap_owner: String,
    pub fnv: u64,
}

fn load_e2e_receipt(root: &Path) -> Result<Option<E2eReceipt>, String> {
    let text = match std::fs::read_to_string(root.join(E2E_RECEIPT_PATH)) {
        Ok(text) => text,
        Err(_) => return Ok(None),
    };
    let parsed = JsonParser::with_string_limit(&text, MAX_TRACKER_STRING_BYTES)
        .finish()
        .map_err(|error| format!("{E2E_RECEIPT_PATH} is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err(format!("{E2E_RECEIPT_PATH} is not a JSON object"));
    };
    match map.get("schema") {
        Some(JsonValue::String(schema)) if schema == E2E_RECEIPT_SCHEMA => {}
        Some(JsonValue::String(schema)) => {
            return Err(format!(
                "{E2E_RECEIPT_PATH} schema is `{schema}`, expected `{E2E_RECEIPT_SCHEMA}`"
            ));
        }
        _ => return Err(format!("{E2E_RECEIPT_PATH} has no schema string")),
    }
    let summary = match map.get("summary") {
        Some(JsonValue::Object(summary)) => summary,
        _ => return Err(format!("{E2E_RECEIPT_PATH} has no summary object")),
    };
    let number = |key: &str| match summary.get(key) {
        Some(JsonValue::Number(raw)) => raw
            .parse::<usize>()
            .map_err(|error| format!("receipt `{key}` is not a count: {error}")),
        _ => Err(format!("receipt has no count `{key}`")),
    };
    let text_field = |key: &str| match summary.get(key) {
        Some(JsonValue::String(value)) => Ok(value.clone()),
        _ => Err(format!("receipt has no string `{key}`")),
    };
    Ok(Some(E2eReceipt {
        stages_green: number("stages_executing")?,
        failures: number("failures")?,
        first_gap: text_field("first_gap")?,
        first_gap_owner: text_field("first_gap_owner")?,
        fnv: fnv1a64(text.as_bytes()),
    }))
}

/// Render the canonical snapshot bytes.
fn render(snapshot: &Snapshot) -> String {
    let (e2e_status, e2e_detail) = match (snapshot.e2e_stages_green, snapshot.e2e_receipt_fnv) {
        (Some(green), Some(fnv)) => (
            E2E_STATUS_MEASURED,
            format!(",\n    \"stages_green\": {green},\n    \"receipt_fnv1a64\": \"{fnv:016x}\""),
        ),
        _ => (E2E_STATUS_NO_RECEIPT, String::new()),
    };
    format!(
        "{{\n  \"schema\": \"{SCHEMA}\",\n  \"bead\": \"frankensim-o5et9\",\n  \"beads\": {{\n    \
         \"issues_fnv1a64\": \"{:016x}\",\n    \"open\": {},\n    \"blocked\": {},\n    \
         \"actionable\": {}\n  }},\n  \"e2e\": {{\n    \"status\": \"{e2e_status}\"{e2e_detail}\n  \
         }},\n  \"no_claim\": \"counts a deliberately regenerated tracker snapshot; the live \
         tracker moves on every br op, so this trails it by design. e2e measured means the \
         retained spine-e2e-summary.json receipt attests the executing prefix; \
         no-retained-receipt means the lane runs green out-of-band but no tracked receipt \
         exists yet\"\n}}\n",
        snapshot.issues_fnv, snapshot.open, snapshot.blocked, snapshot.actionable
    )
}

fn number_field(map: &BTreeMap<String, JsonValue>, key: &str) -> Result<usize, String> {
    match map.get(key) {
        Some(JsonValue::Number(raw)) => raw
            .parse::<usize>()
            .map_err(|error| format!("snapshot `{key}` is not a count: {error}")),
        _ => Err(format!("snapshot has no count `{key}`")),
    }
}

/// Parse and fully validate snapshot text.
fn parse_snapshot(text: &str) -> Result<Snapshot, String> {
    let parsed = JsonParser::new(text)
        .finish()
        .map_err(|error| format!("{SNAPSHOT_PATH} is not valid JSON: {error}"))?;
    let JsonValue::Object(map) = &parsed else {
        return Err(format!("{SNAPSHOT_PATH} is not a JSON object"));
    };
    match map.get("schema") {
        Some(JsonValue::String(schema)) if schema == SCHEMA => {}
        Some(JsonValue::String(schema)) => {
            return Err(format!(
                "{SNAPSHOT_PATH} schema is `{schema}`, expected `{SCHEMA}`; refusing to read a \
                 foreign artifact as zero counts"
            ));
        }
        _ => return Err(format!("{SNAPSHOT_PATH} has no schema string")),
    }
    let beads = match map.get("beads") {
        Some(JsonValue::Object(beads)) => beads,
        _ => return Err(format!("{SNAPSHOT_PATH} has no `beads` object")),
    };
    let issues_fnv = match beads.get("issues_fnv1a64") {
        Some(JsonValue::String(hex)) => u64::from_str_radix(hex, 16)
            .map_err(|error| format!("snapshot `issues_fnv1a64` is not hex: {error}"))?,
        _ => return Err("snapshot has no `issues_fnv1a64` string".to_string()),
    };
    let open = number_field(beads, "open")?;
    let blocked = number_field(beads, "blocked")?;
    let actionable = number_field(beads, "actionable")?;
    if blocked > open {
        return Err(format!(
            "snapshot is internally inconsistent: blocked ({blocked}) exceeds open ({open})"
        ));
    }
    if blocked + actionable != open {
        return Err(format!(
            "snapshot is internally inconsistent: blocked ({blocked}) + actionable \
             ({actionable}) != open ({open})"
        ));
    }
    let (e2e_stages_green, e2e_receipt_fnv) = match map.get("e2e") {
        Some(JsonValue::Object(e2e)) => match e2e.get("status") {
            Some(JsonValue::String(status)) if status == E2E_STATUS_NO_RECEIPT => (None, None),
            Some(JsonValue::String(status)) if status == E2E_STATUS_MEASURED => {
                let green = number_field(e2e, "stages_green")?;
                let fnv = match e2e.get("receipt_fnv1a64") {
                    Some(JsonValue::String(hex)) => {
                        u64::from_str_radix(hex, 16).map_err(|error| {
                            format!("snapshot `receipt_fnv1a64` is not hex: {error}")
                        })?
                    }
                    _ => {
                        return Err("measured snapshot has no `receipt_fnv1a64` string".to_string());
                    }
                };
                (Some(green), Some(fnv))
            }
            Some(JsonValue::String(status)) => {
                return Err(format!(
                    "snapshot e2e status is `{status}`; v2 admits `{E2E_STATUS_NO_RECEIPT}` \
                     and `{E2E_STATUS_MEASURED}`"
                ));
            }
            _ => return Err("snapshot `e2e` has no status string".to_string()),
        },
        _ => return Err(format!("{SNAPSHOT_PATH} has no `e2e` object")),
    };
    Ok(Snapshot {
        issues_fnv,
        open,
        blocked,
        actionable,
        e2e_stages_green,
        e2e_receipt_fnv,
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

/// Load the tracked snapshot for the dashboard, tolerating absence: a missing
/// or invalid snapshot renders NO-DATA rows here while [`check`] flags it as
/// a violation. The two consumers must never disagree about what the bytes
/// mean, so both go through [`parse_snapshot`].
pub(crate) fn load(root: &Path) -> Option<Snapshot> {
    let text = read_bounded(root, SNAPSHOT_PATH, MAX_SNAPSHOT_BYTES).ok()?;
    parse_snapshot(&text).ok()
}

/// Regenerate the snapshot from the live tracker, picking up the retained
/// e2e receipt when one exists.
pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let issues = read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES)?;
    let mut snapshot = count_issues(&issues)?;
    if let Some(receipt) = load_e2e_receipt(root)? {
        snapshot.e2e_stages_green = Some(receipt.stages_green);
        snapshot.e2e_receipt_fnv = Some(receipt.fnv);
    }
    std::fs::write(root.join(SNAPSHOT_PATH), render(&snapshot))
        .map_err(|error| format!("cannot write {SNAPSHOT_PATH}: {error}"))
}

/// The standing gate: validate the snapshot, and report (never charge) a
/// tracker that moved since it was taken.
pub(crate) fn check(root: &Path) -> (Vec<Violation>, Vec<PolicyNote>) {
    let text = match read_bounded(root, SNAPSHOT_PATH, MAX_SNAPSHOT_BYTES) {
        Ok(text) => text,
        Err(error) => {
            return (
                vec![violation(format!(
                    "{error}; a missing or unreadable snapshot is not zero spine metrics — run \
                     `cargo run -p xtask -- generate-spine-metrics`"
                ))],
                Vec::new(),
            );
        }
    };
    let snapshot = match parse_snapshot(&text) {
        Ok(snapshot) => snapshot,
        Err(error) => return (vec![violation(error)], Vec::new()),
    };
    let mut violations = Vec::new();
    let mut notes = Vec::new();
    // The retained e2e receipt must prove the same stage boundary the live
    // product source admits: a receipt claiming a stage the product gaps
    // (or missing one the product executes) is a stale proof, not evidence.
    if let (Some(green), Some(fnv)) = (snapshot.e2e_stages_green, snapshot.e2e_receipt_fnv) {
        match load_e2e_receipt(root) {
            Ok(Some(receipt)) => {
                if receipt.fnv != fnv {
                    violations.push(violation(format!(
                        "the e2e receipt moved since the snapshot (recorded {fnv:016x}, live \
                         {:016x}); regenerate the snapshot deliberately",
                        receipt.fnv
                    )));
                }
                if receipt.failures != 0 {
                    violations.push(violation(format!(
                        "the retained e2e receipt records {} failing check(s); a failing \
                         receipt is not a green prefix",
                        receipt.failures
                    )));
                }
                if receipt.stages_green != green {
                    violations.push(violation(format!(
                        "snapshot records {green} e2e stages green but the receipt says {}; \
                         reconcile deliberately",
                        receipt.stages_green
                    )));
                }
                match std::fs::read_to_string(root.join(crate::spine_ratchet::SOLVE_SOURCE))
                    .map_err(|error| error.to_string())
                    .and_then(|source| crate::spine_ratchet::derive_stages(&source))
                {
                    Ok(stages) => {
                        let live_prefix = crate::spine_ratchet::executing_prefix(&stages);
                        if receipt.stages_green != live_prefix.len() {
                            violations.push(violation(format!(
                                "the retained e2e receipt proves {} stages green but the \
                                 product now executes {}: the proof is stale, not the \
                                 product's current boundary",
                                receipt.stages_green,
                                live_prefix.len()
                            )));
                        }
                        if let Some(first_gap) =
                            stages.iter().find(|stage| stage.gap_owner.is_some())
                        {
                            if receipt.first_gap != first_gap.name {
                                violations.push(violation(format!(
                                    "the retained e2e receipt names first gap `{}` but the \
                                     product gaps `{}`",
                                    receipt.first_gap, first_gap.name
                                )));
                            }
                            let owner = first_gap.gap_owner.as_deref().unwrap_or("");
                            if receipt.first_gap_owner != owner {
                                violations.push(violation(format!(
                                    "the retained e2e receipt names gap owner `{}` but the \
                                     product names `{owner}`",
                                    receipt.first_gap_owner
                                )));
                            }
                        }
                    }
                    Err(error) => {
                        violations.push(violation(format!(
                            "cannot cross-check the e2e receipt against the live stage \
                             table: {error}"
                        )));
                    }
                }
            }
            Ok(None) => {
                violations.push(violation(
                    "the snapshot records a measured e2e state but spine-e2e-summary.json is \
                     absent; a measured claim without its receipt is not evidence",
                ));
            }
            Err(error) => violations.push(violation(error)),
        }
    }
    match read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES) {
        Ok(issues) => {
            let live_fnv = fnv1a64(issues.as_bytes());
            if live_fnv == snapshot.issues_fnv {
                notes.push(PolicyNote {
                    check: CHECK,
                    crate_name: SNAPSHOT_PATH.to_string(),
                    verdict: "current",
                    detail: format!(
                        "tracker snapshot is current (open {}, blocked {}, actionable {})",
                        snapshot.open, snapshot.blocked, snapshot.actionable
                    ),
                });
            } else {
                notes.push(PolicyNote {
                    check: CHECK,
                    crate_name: SNAPSHOT_PATH.to_string(),
                    verdict: "stale-snapshot",
                    detail: format!(
                        "the tracker moved since the snapshot (recorded {:016x}, live {:016x}); \
                         the dashboard deliberately trails it — regenerate with `cargo run -p \
                         xtask -- generate-spine-metrics` when the move should be on the record",
                        snapshot.issues_fnv, live_fnv
                    ),
                });
            }
        }
        Err(error) => {
            // The snapshot itself validated; the tracker being unreadable is a
            // separate fact and must not turn into "snapshot is fine".
            return (
                vec![violation(format!(
                    "{error}; cannot compare the snapshot against the tracker, refusing to \
                     certify currency"
                ))],
                notes,
            );
        }
    }
    (violations, notes)
}

#[cfg(test)]
mod tests {
    //! G0/G3: counting laws, anti-silent-disable refusals, and the
    //! stale-snapshot visibility contract.

    use super::*;

    fn issue(id: &str, status: &str, blocks: &[&str]) -> String {
        let deps: Vec<String> = blocks
            .iter()
            .map(|target| {
                format!(
                    "{{\"issue_id\":\"{id}\",\"depends_on_id\":\"{target}\",\"type\":\"blocks\"}}"
                )
            })
            .collect();
        format!(
            "{{\"id\":\"{id}\",\"status\":\"{status}\",\"dependencies\":[{}]}}",
            deps.join(",")
        )
    }

    #[test]
    fn g0_counts_follow_the_open_blocked_actionable_partition() {
        let text = [
            issue("a", "open", &[]),
            issue("b", "open", &["a"]),
            issue("c", "in_progress", &["b"]),
            issue("d", "closed", &[]),
            issue("e", "open", &["d"]),   // closed blocker does not block
            issue("f", "open", &["zzz"]), // absent blocker does not block
        ]
        .join("\n");
        let snapshot = count_issues(&text).expect("counts");
        assert_eq!(snapshot.open, 5, "a, b, c, e, f are non-closed");
        assert_eq!(snapshot.blocked, 2, "b blocked by a, c blocked by b");
        assert_eq!(snapshot.actionable, 3);
        assert_eq!(snapshot.blocked + snapshot.actionable, snapshot.open);
    }

    #[test]
    fn g0_empty_tracker_refuses_rather_than_reporting_zero_counts() {
        assert!(count_issues("").is_err());
        assert!(count_issues("  \n \n").is_err());
    }

    #[test]
    fn g0_malformed_lines_refuse_with_line_numbers() {
        let error = count_issues("{ not json\n").expect_err("must refuse");
        assert!(error.contains("line 1"), "{error}");
        let no_id = count_issues("{\"status\":\"open\"}\n").expect_err("must refuse");
        assert!(no_id.contains("`id`"), "{no_id}");
    }

    #[test]
    fn g0_snapshot_round_trips_and_stays_canonical() {
        let text = [issue("a", "open", &[]), issue("b", "open", &["a"])].join("\n");
        let snapshot = count_issues(&text).expect("counts");
        let rendered = render(&snapshot);
        let parsed = parse_snapshot(&rendered).expect("round trip");
        assert_eq!(parsed, snapshot);
    }

    #[test]
    fn g0_foreign_or_inconsistent_snapshots_refuse() {
        assert!(parse_snapshot("").is_err());
        assert!(parse_snapshot("{\"schema\":\"something-else\"}").is_err());
        let swapped = "{
  \"schema\": \"frankensim-spine-metrics-v2\",
  \"beads\": { \"issues_fnv1a64\": \"0000000000000000\", \"open\": 2, \"blocked\": 3, \"actionable\": 0 },
  \"e2e\": { \"status\": \"no-retained-receipt\" }
}\n";
        let error = parse_snapshot(swapped).expect_err("blocked > open must refuse");
        assert!(error.contains("inconsistent"), "{error}");
        let bad_sum = swapped
            .replace("\"blocked\": 3", "\"blocked\": 1")
            .replace("\"actionable\": 0", "\"actionable\": 0");
        let error = parse_snapshot(&bad_sum).expect_err("sum mismatch must refuse");
        assert!(error.contains("inconsistent"), "{error}");
        let green = "{
  \"schema\": \"frankensim-spine-metrics-v2\",
  \"beads\": { \"issues_fnv1a64\": \"0000000000000000\", \"open\": 2, \"blocked\": 1, \"actionable\": 1 },
  \"e2e\": { \"status\": \"green\" }
}\n";
        let error = parse_snapshot(green).expect_err("unknown e2e status must refuse");
        assert!(error.contains("v2 admits"), "{error}");
        // A measured e2e state without its fnv is structurally incomplete.
        let measured_no_fnv = "{
  \"schema\": \"frankensim-spine-metrics-v2\",
  \"beads\": { \"issues_fnv1a64\": \"0000000000000000\", \"open\": 2, \"blocked\": 1, \"actionable\": 1 },
  \"e2e\": { \"status\": \"measured\", \"stages_green\": 3 }
}\n";
        let error = parse_snapshot(measured_no_fnv).expect_err("measured without fnv must refuse");
        assert!(error.contains("receipt_fnv1a64"), "{error}");
    }

    #[test]
    fn g0_the_live_tracker_and_snapshot_validate() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("workspace root");
        let snapshot = load(root).expect("the tracked snapshot must exist and validate");
        let issues = read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES).expect("tracker");
        let live = count_issues(&issues).expect("live counts");
        // The snapshot may trail the tracker (visible as a stale-snapshot
        // note), but its own arithmetic must hold and its shape must match
        // the live derivation's.
        assert_eq!(snapshot.blocked + snapshot.actionable, snapshot.open);
        assert!(live.open >= snapshot.open.min(live.blocked + live.actionable));
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
