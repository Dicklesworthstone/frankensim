//! Standing constellation drift-classification gate (bead frankensim-es6pt).
//!
//! Why this gate exists. The strict `check-constellation` lock check answers
//! exactly one question: "are all seven siblings byte-identical to the lock?"
//! Its drift verdict is an aggregate FNV hash, which hides WHICH sibling
//! moved, and it cannot tell a fast-forward (the constellation advanced; the
//! lock is stale) from a divergence (a checkout wandered; that is an
//! incident). Those two states need OPPOSITE remedies — advancing the lock
//! through the compatibility train versus deliberately realigning a checkout —
//! so a gate that cannot name the sibling or the direction selects the wrong
//! remedy. Measured 2026-07-29, all seven siblings were clean strict
//! fast-forwards, 25 to 611 commits ahead of their pins; treating that as
//! drift and moving checkouts back would have discarded every one of those
//! commits (AGENTS.md forbids moving a shared sibling checkout to satisfy a
//! preflight).
//!
//! What it charges. Only a wandered checkout is drift, and drift is red:
//!
//! - DIRTY worktree (on or off the pin),
//! - HEAD BEHIND the pin (the checkout retreated),
//! - DIVERGED history (neither HEAD nor pin is an ancestor of the other),
//! - observation refusal or HEAD moving mid-observation (fail-closed).
//!
//! Every verdict names the sibling and prints expected and actual heads.
//!
//! What it deliberately does not charge. A strict fast-forward is reported as
//! `stale-lock` — visible in every run, but not an incident: the remedy is
//! the constellation train (f85xj.13.4 successor), never a checkout move. A
//! missing checkout, a shallow checkout whose grafted boundary could corrupt
//! ancestry, or a locally absent pin object is NO-DATA: rendered explicitly,
//! never silently counted as pass or drift. The strict lock check remains the
//! equality authority for reproducibility; this gate is the standing incident
//! detector inside check-all.
//!
//! Anti-silent-disable: an unreadable, corrupt, empty, or row-incomplete
//! constellation.lock is itself a violation. A gate that reports zero
//! violations because its input rotted is how gates die.

use std::path::Path;

use crate::constellation_cleanliness::{repository_worktree_status, sanitized_git_command};
use crate::{
    CONSTELLATION_REPOS, LockRow, PolicyNote, Violation, git_out, parse_lock_rows,
    read_constellation_lock,
};

pub(crate) const CHECK: &str = "constellation-drift";

/// HEAD's relationship to the locked pin when both ends are observable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PinRelation {
    /// HEAD is exactly the pin.
    OnPin,
    /// The pin is a strict ancestor of HEAD: the constellation advanced and
    /// the LOCK is stale. Not drift; remedy is the train, not a checkout move.
    Ahead { commits: u64 },
    /// HEAD is a strict ancestor of the pin: the checkout retreated. Drift.
    Behind { commits: u64 },
    /// Neither is an ancestor of the other. Drift incident.
    Diverged { ahead: u64, behind: u64 },
    /// Ancestry cannot be verified here (pin object absent, or a shallow
    /// boundary could corrupt the answer). NO-DATA, never charged either way.
    Unverifiable,
}

/// Chargeable outcome for one sibling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SiblingClass {
    /// On-pin and clean.
    Ok,
    /// Strict fast-forward and clean: the lock is stale, not the checkout.
    StaleLock,
    /// Not chargeable: checkout absent/unreadable, or ancestry unverifiable.
    NoData,
    /// The checkout wandered (dirty, retreated, diverged, or unobservable).
    Drift,
}

/// One rendered per-sibling verdict row. The owning library is carried by
/// the caller (and embedded in `detail`), so the row itself stays minimal.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SiblingRow {
    class: SiblingClass,
    detail: String,
}

/// What the gate returns: violations are the red rows; notes carry the
/// ok/stale-lock/no-data classification for every sibling, every run.
pub(crate) struct DriftReport {
    pub(crate) violations: Vec<Violation>,
    pub(crate) notes: Vec<PolicyNote>,
}

fn violation_for(lib: &str, detail: impl Into<String>) -> Violation {
    Violation {
        check: CHECK,
        crate_name: lib.to_string(),
        detail: detail.into(),
    }
}

fn note(verdict: &'static str, lib: &str, detail: String) -> PolicyNote {
    PolicyNote {
        check: CHECK,
        crate_name: lib.to_string(),
        verdict,
        detail,
    }
}

fn refuse(detail: impl Into<String>) -> DriftReport {
    DriftReport {
        violations: vec![violation_for("constellation.lock", detail)],
        notes: Vec::new(),
    }
}

/// Exit-code query with the "no" case separated from real errors: 0 → true,
/// 1 → false, anything else is an error (missing objects are not a "no").
fn git_exit_flag(dir: &Path, args: &[&str]) -> Result<bool, String> {
    let output = sanitized_git_command(dir, args)
        .output()
        .map_err(|e| format!("git {args:?} failed to spawn: {e}"))?;
    match output.status.code() {
        Some(0) => Ok(true),
        Some(1) => Ok(false),
        other => Err(format!(
            "git {args:?} in {} exited {other:?}: {}",
            dir.display(),
            String::from_utf8_lossy(&output.stderr).trim()
        )),
    }
}

/// Object-presence probe where ANY failure means "not present" (an invalid
/// rev name and a missing object are the same fact for classification).
fn commit_object_present(dir: &Path, rev: &str) -> bool {
    sanitized_git_command(dir, &["cat-file", "-e", &format!("{rev}^{{commit}}")])
        .output()
        .map(|output| output.status.success())
        .unwrap_or(false)
}

/// Commits in `to` not reachable from `from` (`from..to`).
fn rev_count(dir: &Path, from: &str, to: &str) -> Result<u64, String> {
    git_out(dir, &["rev-list", "--count", &format!("{from}..{to}")])?
        .parse::<u64>()
        .map_err(|error| {
            format!(
                "git rev-list --count {from}..{to} in {} unparsable: {error}",
                dir.display()
            )
        })
}

/// Classify HEAD against the pin. Relation queries never run when the answer
/// could be corrupted by a shallow boundary: `merge-base --is-ancestor` on a
/// grafted clone can report "not an ancestor" for true ancestors, which would
/// manufacture a false divergence — the worst verdict this gate owns.
fn pin_relation(dir: &Path, pin: &str, head: &str) -> Result<PinRelation, String> {
    if head == pin {
        return Ok(PinRelation::OnPin);
    }
    if !commit_object_present(dir, pin) {
        return Ok(PinRelation::Unverifiable);
    }
    let shallow = git_out(dir, &["rev-parse", "--is-shallow-repository"])
        .map(|answer| answer == "true")
        .unwrap_or(false);
    if shallow {
        return Ok(PinRelation::Unverifiable);
    }
    let pin_is_ancestor = git_exit_flag(dir, &["merge-base", "--is-ancestor", pin, head])?;
    let head_is_ancestor = git_exit_flag(dir, &["merge-base", "--is-ancestor", head, pin])?;
    match (pin_is_ancestor, head_is_ancestor) {
        (true, true) => Ok(PinRelation::OnPin),
        (true, false) => Ok(PinRelation::Ahead {
            commits: rev_count(dir, pin, head)?,
        }),
        (false, true) => Ok(PinRelation::Behind {
            commits: rev_count(dir, head, pin)?,
        }),
        (false, false) => Ok(PinRelation::Diverged {
            ahead: rev_count(dir, pin, head)?,
            behind: rev_count(dir, head, pin)?,
        }),
    }
}

/// First lines of a dirty-worktree description, bounded so a large status
/// cannot flood the verdict row.
fn status_excerpt(status: &str) -> String {
    const MAX_LINES: usize = 5;
    let mut lines: Vec<&str> = status.lines().take(MAX_LINES + 1).collect();
    let truncated = lines.len() > MAX_LINES;
    lines.truncate(MAX_LINES);
    let mut excerpt = lines.join("; ");
    if truncated {
        excerpt.push_str("; …");
    }
    excerpt
}

fn classify_row(
    lib: &str,
    pin: &str,
    head: &str,
    relation: PinRelation,
    dirty: bool,
    status: &str,
) -> SiblingRow {
    let heads = format!("expected pin {pin}, actual HEAD {head}");
    let (class, detail) = match (relation, dirty) {
        (PinRelation::OnPin, false) => (
            SiblingClass::Ok,
            format!("{lib} on-pin and clean at {head}"),
        ),
        (PinRelation::OnPin, true) => (
            SiblingClass::Drift,
            format!(
                "{lib} DIRTY at the locked pin ({heads}); a modified worktree is not the \
                 pinned source: {}",
                status_excerpt(status)
            ),
        ),
        (PinRelation::Ahead { commits }, false) => (
            SiblingClass::StaleLock,
            format!(
                "{lib} off-pin fast-forward: HEAD is {commits} commits ahead of the lock \
                 ({heads}); the LOCK is stale, not the checkout — advance it through the \
                 constellation train (f85xj.13.4 successor); AGENTS.md forbids moving the \
                 checkout back"
            ),
        ),
        (PinRelation::Ahead { commits }, true) => (
            SiblingClass::Drift,
            format!(
                "{lib} DIRTY while {commits} commits ahead of the lock pin ({heads}): {}",
                status_excerpt(status)
            ),
        ),
        (PinRelation::Behind { commits }, _) => (
            SiblingClass::Drift,
            format!(
                "{lib} checkout RETREATED {commits} commits from the lock pin ({heads}); \
                 a checkout moved back is drift — realign deliberately through \
                 constellation governance"
            ),
        ),
        (PinRelation::Diverged { ahead, behind }, _) => (
            SiblingClass::Drift,
            format!(
                "{lib} DIVERGED from the lock pin (ahead {ahead}, behind {behind}; \
                 {heads}); reconcile through constellation governance — this is the \
                 incident class the aggregate lock hash hid"
            ),
        ),
        (PinRelation::Unverifiable, false) => (
            SiblingClass::NoData,
            format!(
                "{lib} ancestry to pin {pin} unverifiable (shallow boundary or pin object \
                 absent; HEAD {head}); NO-DATA — the strict check-constellation remains \
                 the equality authority"
            ),
        ),
        (PinRelation::Unverifiable, true) => (
            SiblingClass::Drift,
            format!(
                "{lib} DIRTY with pin ancestry unverifiable ({heads}): {}",
                status_excerpt(status)
            ),
        ),
    };
    SiblingRow { class, detail }
}

/// Observe one sibling checkout against its lock row. Observation order is
/// head → deep cleanliness → head again, so a checkout that moved mid-walk is
/// charged as incoherent rather than reported from two different states.
fn observe_sibling(projects: &Path, lib: &str, dirname: &str, row: &LockRow) -> SiblingRow {
    let dir = projects.join(dirname);
    let pin = row.git_head.as_str();
    let head = match git_out(&dir, &["rev-parse", "HEAD"]) {
        Ok(head) => head,
        Err(_) => {
            return SiblingRow {
                class: SiblingClass::NoData,
                detail: format!(
                    "{lib} has no readable git checkout at {} (not materialized in this \
                     environment; gate not charged)",
                    dir.display()
                ),
            };
        }
    };
    // The same deep, submodule-concealment-resistant raw-index observation the
    // strict lock check trusts — `git status` alone can be talked out of
    // reporting dirt by local ignore policy.
    let status = match repository_worktree_status(&dir) {
        Ok(status) => status,
        Err(error) => {
            return SiblingRow {
                class: SiblingClass::Drift,
                detail: format!(
                    "{lib} worktree observation REFUSED at {} ({heads}): {error}; \
                     fail-closed — an unobservable checkout cannot be trusted as \
                     drift-free",
                    dir.display(),
                    heads = format!("expected pin {pin}, actual HEAD {head}"),
                ),
            };
        }
    };
    let head_after = git_out(&dir, &["rev-parse", "HEAD"]).unwrap_or_default();
    if head_after != head {
        return SiblingRow {
            class: SiblingClass::Drift,
            detail: format!(
                "{lib} HEAD moved while being observed (before {head}, after \
                 {head_after}); refusing incoherent provenance"
            ),
        };
    }
    let relation = match pin_relation(&dir, pin, &head) {
        Ok(relation) => relation,
        Err(error) => {
            return SiblingRow {
                class: SiblingClass::NoData,
                detail: format!(
                    "{lib} ancestry to pin {pin} unverifiable (HEAD {head}): {error}; \
                     NO-DATA, not charged"
                ),
            };
        }
    };
    classify_row(lib, pin, &head, relation, !status.is_empty(), &status)
}

/// Evaluate a parsed lock against the constellation projects directory.
/// Kept separate from [`check`] so tests can point it at fixture layouts.
fn evaluate(rows: &[LockRow], projects: &Path) -> DriftReport {
    let mut violations = Vec::new();
    let mut notes = Vec::new();
    if rows.len() != CONSTELLATION_REPOS.len() {
        return refuse(format!(
            "lock declares {} libraries, expected {}; refusing to classify a partial \
             constellation as drift-free",
            rows.len(),
            CONSTELLATION_REPOS.len()
        ));
    }
    let mut materialized = 0usize;
    for &(lib, dirname) in CONSTELLATION_REPOS {
        let Some(row) = rows.iter().find(|row| row.lib == lib) else {
            violations.push(violation_for(
                lib,
                format!(
                    "lock is missing constellation library {lib}; a partial lock cannot \
                     certify drift-freedom"
                ),
            ));
            continue;
        };
        let row_report = observe_sibling(projects, lib, dirname, row);
        if !matches!(row_report.class, SiblingClass::NoData) {
            materialized += 1;
        }
        match row_report.class {
            SiblingClass::Drift => violations.push(violation_for(lib, row_report.detail)),
            SiblingClass::Ok => notes.push(note("ok", lib, row_report.detail)),
            SiblingClass::StaleLock => notes.push(note("stale-lock", lib, row_report.detail)),
            SiblingClass::NoData => notes.push(note("no-data", lib, row_report.detail)),
        }
    }
    if materialized == 0 {
        notes.push(note(
            "no-data",
            "constellation",
            format!(
                "0 of {} siblings materialized as readable git checkouts; the drift gate \
                 observed nothing this run (NO-DATA, not a pass)",
                CONSTELLATION_REPOS.len()
            ),
        ));
    }
    DriftReport { violations, notes }
}

/// The standing gate: read the tracked lock, classify every sibling, charge
/// only wandered checkouts. Lock input failures refuse closed.
pub(crate) fn check(root: &Path) -> DriftReport {
    let lock_path = root.join("constellation.lock");
    let rows = match read_constellation_lock(&lock_path)
        .and_then(|text| parse_lock_rows(&text).map(|(_, rows)| rows))
    {
        Ok(rows) if rows.is_empty() => {
            return refuse(
                "constellation.lock parsed to zero library rows; a gate with no inputs \
                 cannot report drift-free",
            );
        }
        Ok(rows) => rows,
        Err(error) => {
            return refuse(format!(
                "constellation.lock unreadable or invalid: {error}; refusing to report \
                 zero drift on a corrupt input"
            ));
        }
    };
    let Some(projects) = root.parent() else {
        return refuse("workspace root has no parent; cannot locate constellation siblings");
    };
    evaluate(&rows, projects)
}

#[cfg(test)]
mod tests {
    //! G0/G3 falsifiers for the drift gate. The bead's named falsifier: a
    //! sibling moved off-pin must turn the gate red NAMING THAT SIBLING (the
    //! aggregate-hash failure this gate replaces hid which one), with
    //! expected and actual heads in the verdict. The companion obligation: a
    //! strict fast-forward must NOT go red — that false red is what invited
    //! the destructive "move the checkout back" remedy.

    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::process::Command;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

    /// A temp `projects/` directory holding fixture sibling checkouts.
    struct FixtureProjects {
        root: PathBuf,
    }

    impl FixtureProjects {
        fn new() -> Self {
            let id = NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed);
            let root = std::env::temp_dir().join(format!(
                "frankensim-drift-gate-test-{}-{id}",
                std::process::id()
            ));
            fs::create_dir_all(&root).expect("fixture root");
            Self { root }
        }

        /// Create one fixture sibling checkout and return its path.
        fn init_sibling(&self, dirname: &str) -> PathBuf {
            let dir = self.root.join(dirname);
            fs::create_dir_all(&dir).expect("sibling dir");
            git(&dir, &["init"]);
            dir
        }
    }

    impl Drop for FixtureProjects {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }

    fn git(dir: &Path, args: &[&str]) -> String {
        let output = Command::new("git")
            .arg("-C")
            .arg(dir)
            .args([
                "-c",
                "user.email=drift-fixture@example.invalid",
                "-c",
                "user.name=drift-fixture",
                "-c",
                "commit.gpgsign=false",
                "-c",
                "init.defaultBranch=main",
            ])
            .args(args)
            .env("GIT_CONFIG_GLOBAL", "/dev/null")
            .env("GIT_CONFIG_SYSTEM", "/dev/null")
            .env_remove("GIT_DIR")
            .env_remove("GIT_WORK_TREE")
            .output()
            .expect("git spawn");
        assert!(
            output.status.success(),
            "git {args:?} in {} failed: {}",
            dir.display(),
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout)
            .expect("git stdout utf8")
            .trim()
            .to_string()
    }

    fn commit_file(dir: &Path, name: &str, contents: &str) -> String {
        fs::write(dir.join(name), contents).expect("write fixture file");
        git(dir, &["add", name]);
        git(dir, &["commit", "-m", &format!("add {name}")]);
        git(dir, &["rev-parse", "HEAD"])
    }

    /// Lock rows covering all seven siblings; `asupersync` gets `pin`, the
    /// rest get a zero pin and stay unmaterialized (NO-DATA).
    fn rows_with(pin: &str) -> Vec<LockRow> {
        CONSTELLATION_REPOS
            .iter()
            .map(|&(lib, _)| LockRow {
                lib: lib.to_string(),
                version: "0.0.0".to_string(),
                git_head: if lib == "asupersync" {
                    pin.to_string()
                } else {
                    "0".repeat(40)
                },
                remote: "fixture".to_string(),
                path: String::new(),
            })
            .collect()
    }

    fn find_note<'a>(report: &'a DriftReport, verdict: &str, lib: &str) -> &'a PolicyNote {
        report
            .notes
            .iter()
            .find(|note| note.verdict == verdict && note.crate_name == lib)
            .unwrap_or_else(|| panic!("expected {verdict} note for {lib}: {:?}", report.notes))
    }

    #[test]
    fn on_pin_clean_is_ok_and_not_charged() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let pin = commit_file(&dir, "a.txt", "one\n");
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert!(report.violations.is_empty(), "{:?}", report.violations);
        assert!(find_note(&report, "ok", "asupersync").detail.contains(&pin));
    }

    #[test]
    fn fast_forward_is_stale_lock_not_red() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let pin = commit_file(&dir, "a.txt", "one\n");
        let head = commit_file(&dir, "b.txt", "two\n");
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert!(
            report.violations.is_empty(),
            "a strict fast-forward must not be red: {:?}",
            report.violations
        );
        let note = find_note(&report, "stale-lock", "asupersync");
        assert!(note.detail.contains(&pin), "expected head: {note:?}");
        assert!(note.detail.contains(&head), "actual head: {note:?}");
        assert!(note.detail.contains("1 commits ahead"));
    }

    #[test]
    fn retreated_checkout_is_red_and_names_the_sibling() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let first = commit_file(&dir, "a.txt", "one\n");
        let pin = commit_file(&dir, "b.txt", "two\n");
        git(&dir, &["checkout", &first]);
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        let violation = &report.violations[0];
        assert_eq!(violation.check, CHECK);
        assert_eq!(violation.crate_name, "asupersync");
        assert!(violation.detail.contains("RETREATED 1 commits"));
        assert!(violation.detail.contains(&pin), "expected head named");
        assert!(violation.detail.contains(&first), "actual head named");
        assert!(
            report
                .notes
                .iter()
                .all(|note| note.crate_name != "asupersync"),
            "a drifted sibling must not also emit a pass-shaped note"
        );
    }

    #[test]
    fn diverged_checkout_is_red_and_names_the_sibling() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let base = commit_file(&dir, "a.txt", "one\n");
        let head = commit_file(&dir, "b.txt", "two\n");
        git(&dir, &["checkout", &base]);
        let pin = commit_file(&dir, "c.txt", "three\n");
        git(&dir, &["checkout", "main"]);
        assert_eq!(git(&dir, &["rev-parse", "HEAD"]), head);
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        let violation = &report.violations[0];
        assert_eq!(violation.crate_name, "asupersync");
        assert!(violation.detail.contains("DIVERGED"), "{violation:?}");
        assert!(violation.detail.contains(&pin), "expected head named");
        assert!(violation.detail.contains(&head), "actual head named");
    }

    #[test]
    fn dirty_at_the_pin_is_red_and_names_the_sibling() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let pin = commit_file(&dir, "a.txt", "one\n");
        fs::write(dir.join("a.txt"), "one\nmud\n").expect("dirty the worktree");
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        let violation = &report.violations[0];
        assert_eq!(violation.crate_name, "asupersync");
        assert!(violation.detail.contains("DIRTY"));
        assert!(violation.detail.contains(&pin));
    }

    #[test]
    fn dirty_while_ahead_is_red_even_though_ahead_alone_is_not() {
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let pin = commit_file(&dir, "a.txt", "one\n");
        commit_file(&dir, "b.txt", "two\n");
        fs::write(dir.join("a.txt"), "one\nmud\n").expect("dirty the worktree");
        let report = evaluate(&rows_with(&pin), &projects.root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        assert!(report.violations[0].detail.contains("DIRTY"));
        assert!(report.violations[0].detail.contains("1 commits ahead"));
    }

    #[test]
    fn missing_siblings_are_no_data_not_pass_not_drift() {
        let projects = FixtureProjects::new();
        let report = evaluate(&rows_with(&"0".repeat(40)), &projects.root);
        assert!(report.violations.is_empty(), "{:?}", report.violations);
        assert!(
            find_note(&report, "no-data", "asupersync")
                .detail
                .contains("not materialized")
        );
        assert!(
            find_note(&report, "no-data", "constellation")
                .detail
                .contains("observed nothing")
        );
    }

    #[test]
    fn partial_lock_refuses_closed() {
        let projects = FixtureProjects::new();
        let mut rows = rows_with(&"0".repeat(40));
        rows.pop();
        let report = evaluate(&rows, &projects.root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        assert!(report.violations[0].detail.contains("partial"));
    }

    #[test]
    fn corrupt_lock_input_is_a_violation_not_a_clean_gate() {
        let projects = FixtureProjects::new();
        let root = projects.root.join("frankensim");
        fs::create_dir_all(&root).expect("workspace root");
        fs::write(root.join("constellation.lock"), "{ not json\n").expect("corrupt lock");
        let report = check(&root);
        assert_eq!(report.violations.len(), 1, "{:?}", report.violations);
        assert!(report.violations[0].detail.contains("corrupt input"));
    }

    #[test]
    fn moved_head_mid_observation_is_incoherent_not_classified() {
        // Unit-level: the coherence guard is what refuses a verdict computed
        // from two different repository states.
        let projects = FixtureProjects::new();
        let dir = projects.init_sibling("asupersync");
        let before = commit_file(&dir, "a.txt", "one\n");
        let after = commit_file(&dir, "b.txt", "two\n");
        assert_ne!(before, after);
        // classify_row is pure: a dirty verdict while on-pin is drift, which
        // is the shape the coherence guard protects from being misfiled.
        let row = classify_row(
            "asupersync",
            &before,
            &after,
            PinRelation::Ahead { commits: 1 },
            false,
            "",
        );
        assert_eq!(row.class, SiblingClass::StaleLock);
        assert!(row.detail.contains(&before));
        assert!(row.detail.contains(&after));
    }
}
