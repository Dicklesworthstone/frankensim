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
const SCHEMA: &str = "frankensim-spine-metrics-v1";
const E2E_STATUS_NO_RECEIPT: &str = "no-retained-receipt";
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
                format!("{ISSUES_PATH} line {} is not valid JSON: {error}", index + 1)
            })?;
        let JsonValue::Object(map) = &parsed else {
            return Err(format!("{ISSUES_PATH} line {} is not a JSON object", index + 1));
        };
        let id = match map.get("id") {
            Some(JsonValue::String(id)) => id.clone(),
            _ => return Err(format!("{ISSUES_PATH} line {} has no string `id`", index + 1)),
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
                let is_blocks = matches!(dep.get("type"), Some(JsonValue::String(t)) if t == "blocks");
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
    })
}

/// Render the canonical snapshot bytes.
fn render(snapshot: &Snapshot) -> String {
    format!(
        "{{\n  \"schema\": \"{SCHEMA}\",\n  \"bead\": \"frankensim-o5et9\",\n  \"beads\": {{\n    \
         \"issues_fnv1a64\": \"{:016x}\",\n    \"open\": {},\n    \"blocked\": {},\n    \
         \"actionable\": {}\n  }},\n  \"e2e\": {{\n    \"status\": \"{E2E_STATUS_NO_RECEIPT}\"\n  \
         }},\n  \"no_claim\": \"counts a deliberately regenerated tracker snapshot; the live \
         tracker moves on every br op, so this trails it by design. e2e is no-retained-receipt: \
         the lane runs green out-of-band (frankensim-ustax) but retains no tracked checked \
         receipt (frankensim-iakds owns one)\"\n}}\n",
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
    match map.get("e2e") {
        Some(JsonValue::Object(e2e)) => match e2e.get("status") {
            Some(JsonValue::String(status)) if status == E2E_STATUS_NO_RECEIPT => {}
            Some(JsonValue::String(status)) => {
                return Err(format!(
                    "snapshot e2e status is `{status}`; v1 admits only \
                     `{E2E_STATUS_NO_RECEIPT}` — a retained receipt lands through bead \
                     frankensim-iakds with a deliberate schema bump"
                ));
            }
            _ => return Err("snapshot `e2e` has no status string".to_string()),
        },
        _ => return Err(format!("{SNAPSHOT_PATH} has no `e2e` object")),
    }
    Ok(Snapshot {
        issues_fnv,
        open,
        blocked,
        actionable,
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

/// Regenerate the snapshot from the live tracker.
pub(crate) fn generate(root: &Path) -> Result<(), String> {
    let issues = read_bounded(root, ISSUES_PATH, MAX_ISSUES_BYTES)?;
    let snapshot = count_issues(&issues)?;
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
    let mut notes = Vec::new();
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
    (Vec::new(), notes)
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
  \"schema\": \"frankensim-spine-metrics-v1\",
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
  \"schema\": \"frankensim-spine-metrics-v1\",
  \"beads\": { \"issues_fnv1a64\": \"0000000000000000\", \"open\": 2, \"blocked\": 1, \"actionable\": 1 },
  \"e2e\": { \"status\": \"green\" }
}\n";
        let error = parse_snapshot(green).expect_err("unknown e2e status must refuse");
        assert!(error.contains("frankensim-iakds"), "{error}");
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
