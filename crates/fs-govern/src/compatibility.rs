//! Cross-repository compatibility suite and release-train protocol.
//!
//! Sibling pins move by ad hoc decision today, and nothing systematically
//! answers "does FrankenSim still hold its contracts on the new commit?"
//! before a bump lands. This module makes a pin change a TESTED, RECORDED
//! event: it registers, per sibling, the load-bearing claims FrankenSim
//! actually depends on and the FrankenSim tests that exercise them, compares
//! a recorded pin set against a live or candidate one, and adjudicates a bump
//! attempt fail-closed.
//!
//! Rules, all structural rather than advisory:
//! - **A bump is refused unless the suite was EXECUTED and green.** `NotRun`
//!   is never a pass, and an `Executed` result with zero tests is refused too:
//!   a selector that matched nothing must not read as success.
//! - **Every moved sibling must be covered**, and every `P1`/`P2` surface must
//!   be green regardless of whether it moved, because one sibling's move can
//!   break another's surface.
//! - **The golden disposition must be declared.** Pins and semantic goldens
//!   move together; an undeclared golden implication is a refusal, not a
//!   default-pass.
//! - **Emergency is a classification, not an adjective.** Convenience, new
//!   upstream features, and version freshness are unrepresentable as
//!   [`EmergencyClass`], so an emergency bump cannot be justified by them.
//! - **A refusal is total.** No partially admitted bump exists.
//!
//! This module is deliberately sibling-free (its crate's cone is `fs-blake3`,
//! `fs-evidence`, `fs-vvreg`, `fs-wedge`) so the machinery that adjudicates a
//! sibling bump keeps building when a sibling is broken. Hosting it anywhere
//! that reaches `fs-exec`/asupersync or `fs-ledger`/FrankenSQLite would make
//! the adjudicator unavailable during exactly the outages it exists for.
//!
//! No-claim boundary: this module records and adjudicates. It does not run
//! tests, read `constellation.lock`, or observe a git checkout — callers
//! supply pin rows and executed outcomes. A green verdict means the declared
//! suite passed as reported, never that the sibling is correct.

use core::fmt::{self, Write as _};

/// Schema version of the compatibility registry and bump verdicts.
pub const COMPATIBILITY_SCHEMA_VERSION: u32 = 1;

/// Review priority from the constellation trust-cone assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ReviewPriority {
    /// Critical correctness surface; a bump requires green evidence.
    P1,
    /// High-value surface; a bump requires green evidence.
    P2,
    /// Secondary surface; required only when the sibling itself moved.
    P3,
    /// Pinned but unused; carries no compatibility claim.
    P4,
}

impl ReviewPriority {
    /// Stable machine name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::P1 => "P1",
            Self::P2 => "P2",
            Self::P3 => "P3",
            Self::P4 => "P4",
        }
    }

    /// Whether a green suite is required even when this sibling did not move.
    #[must_use]
    pub const fn required_every_train(self) -> bool {
        matches!(self, Self::P1 | Self::P2)
    }
}

/// One FrankenSim test that exercises a sibling's load-bearing surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SurfaceTest {
    /// FrankenSim crate hosting the test.
    pub crate_name: &'static str,
    /// Integration-test target name (the file stem under `tests/`).
    pub test_target: &'static str,
    /// Exact test function name.
    pub test_name: &'static str,
    /// The claim this test exercises.
    pub claim: &'static str,
}

/// One sibling's compatibility surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiblingSurface {
    /// Sibling library name as it appears in `constellation.lock`.
    pub lib: &'static str,
    /// What FrankenSim uses it for.
    pub role: &'static str,
    /// Review priority from the trust-cone assessment.
    pub priority: ReviewPriority,
    /// Load-bearing claims FrankenSim depends on.
    pub claims: &'static [&'static str],
    /// FrankenSim tests exercising those claims.
    pub tests: &'static [SurfaceTest],
    /// FrankenSim crates that consume this sibling at RUNTIME.
    ///
    /// Golden surfaces are identified as `crate:surface`, so this is what
    /// couples a pin move to the semantic goldens it can move underneath.
    pub runtime_consumers: &'static [&'static str],
    /// Why a surface carries no tests, when it carries none.
    pub no_test_reason: Option<&'static str>,
}

impl SiblingSurface {
    /// The exact `cargo test` selector for this surface.
    ///
    /// Returns `None` for a surface with no tests, so a caller cannot render a
    /// runnable command for a sibling that has no compatibility coverage.
    #[must_use]
    pub fn selector(&self) -> Option<String> {
        if self.tests.is_empty() {
            return None;
        }
        let mut crates: Vec<&str> = self.tests.iter().map(|test| test.crate_name).collect();
        crates.sort_unstable();
        crates.dedup();
        let mut out = String::from("cargo test --locked");
        for name in crates {
            let _ = write!(out, " -p {name}");
        }
        let mut targets: Vec<&str> = self.tests.iter().map(|test| test.test_target).collect();
        targets.sort_unstable();
        targets.dedup();
        for target in targets {
            let _ = write!(out, " --test {target}");
        }
        Some(out)
    }

    /// Number of tests registered for this surface.
    #[must_use]
    pub const fn test_count(&self) -> usize {
        self.tests.len()
    }
}

/// The registered compatibility surfaces, in canonical `lib` order.
///
/// Derived from the measured trust-cone assessment, not from plan prose: the
/// priorities and roles mirror `constellation-trust-assessment.json`, and
/// every listed test exists in the tree.
pub const SURFACES: &[SiblingSurface] = &[
    SiblingSurface {
        lib: "asupersync",
        role: "structured concurrency, cancellation, scopes, budgets, latency lane",
        priority: ReviewPriority::P1,
        claims: &[
            "bounded cancellation with request -> drain -> finalize",
            "task-scoped capability and budget propagation",
            "deterministic pause/drain/replay boundaries",
        ],
        tests: &[
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "conformance",
                test_name: "exec_004_external_cancellation_drains_and_ledgers_latency",
                claim: "bounded cancellation with request -> drain -> finalize",
            },
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "conformance",
                test_name: "exec_005_g4_storm_random_cancels_and_panics_stay_structured",
                claim: "bounded cancellation with request -> drain -> finalize",
            },
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "conformance",
                test_name: "exec_007_latency_lane_stays_responsive_under_tile_load",
                claim: "task-scoped capability and budget propagation",
            },
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "conformance",
                test_name: "exec_010_race_winner_is_deterministic_and_losers_fully_drain",
                claim: "deterministic pause/drain/replay boundaries",
            },
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "lease_battery",
                test_name: "cancellation_releases_all_charges",
                claim: "bounded cancellation with request -> drain -> finalize",
            },
            SurfaceTest {
                crate_name: "fs-exec",
                test_target: "constellation_smoke",
                test_name: "asupersync_links_and_budget_vocabulary_holds",
                claim: "task-scoped capability and budget propagation",
            },
        ],
        runtime_consumers: &["fs-exec", "fs-plan", "fs-surrogate"],
        no_test_reason: None,
    },
    SiblingSurface {
        lib: "franken_networkx",
        role: "graph algorithms behind voxel, sparse, and truss structure",
        priority: ReviewPriority::P3,
        claims: &["graph traversal and connectivity results are deterministic"],
        tests: &[],
        runtime_consumers: &["fs-rep-voxel", "fs-sparse", "fs-truss"],
        no_test_reason: Some(
            "no test isolates the franken_networkx boundary from FrankenSim's own graph logic; \
             coverage is incidental inside fs-sparse/fs-truss suites and is not claimed here",
        ),
    },
    SiblingSurface {
        lib: "franken_numpy",
        role: "array primitives used by sparse assembly",
        priority: ReviewPriority::P2,
        claims: &["dense array semantics used by sparse assembly are stable"],
        tests: &[],
        runtime_consumers: &["fs-sparse"],
        no_test_reason: Some(
            "no boundary-isolating test exists yet; fs-sparse exercises it only incidentally, so \
             this surface is an explicit gap rather than a covered one",
        ),
    },
    SiblingSurface {
        lib: "frankenpandas",
        role: "pinned but unused",
        priority: ReviewPriority::P4,
        claims: &[],
        tests: &[],
        runtime_consumers: &[],
        no_test_reason: Some(
            "pinned-unused: no runtime or dev consumer exists, so FrankenSim makes no claim that \
             depends on it and none is tested",
        ),
    },
    SiblingSurface {
        lib: "frankenscipy",
        role: "development-only differential oracle",
        priority: ReviewPriority::P2,
        claims: &["oracle casebook values are stable across the pin"],
        tests: &[],
        runtime_consumers: &[],
        no_test_reason: Some(
            "oracle casebooks are dev-only comparisons spread across seven crates; they are not \
             yet consolidated into a selectable compatibility target",
        ),
    },
    SiblingSurface {
        lib: "frankensqlite",
        role: "durable storage under the design ledger, receipts, and replay",
        priority: ReviewPriority::P1,
        claims: &[
            "artifact and lineage durability across crash and reopen",
            "schema migration and refusal semantics",
            "transaction and checkpoint boundaries hold under interruption",
        ],
        tests: &[
            SurfaceTest {
                crate_name: "fs-ledger",
                test_target: "conformance",
                test_name: "ledger_003_schema_migration_versioned",
                claim: "schema migration and refusal semantics",
            },
            SurfaceTest {
                crate_name: "fs-ledger",
                test_target: "conformance",
                test_name: "ledger_007_crash_kill9_battery",
                claim: "artifact and lineage durability across crash and reopen",
            },
            SurfaceTest {
                crate_name: "fs-ledger",
                test_target: "travel",
                test_name: "tt_006_crash_kill9_during_fork_traffic",
                claim: "artifact and lineage durability across crash and reopen",
            },
            SurfaceTest {
                crate_name: "fs-ledger",
                test_target: "state_checkpoint",
                test_name: "committed_checkpoint_prefix_survives_kill_and_real_file_reopen",
                claim: "transaction and checkpoint boundaries hold under interruption",
            },
            SurfaceTest {
                crate_name: "fs-ledger",
                test_target: "ambient_cx",
                test_name: "latency_lane_ambient_context_reaches_fsqlite_waiters",
                claim: "transaction and checkpoint boundaries hold under interruption",
            },
        ],
        runtime_consumers: &["fs-ledger", "fs-vskeleton"],
        no_test_reason: None,
    },
    SiblingSurface {
        lib: "frankentorch",
        role: "tensor bridge behind the feature-gated AD path",
        priority: ReviewPriority::P3,
        claims: &["the feature-gated tape bridge preserves gradient values"],
        tests: &[],
        runtime_consumers: &["fs-ad"],
        no_test_reason: Some(
            "the bridge is feature-gated and off by default; no default-path test exercises it, \
             so no compatibility claim is made",
        ),
    },
];

/// Golden surfaces coupled to one sibling, out of a caller-supplied set.
///
/// Golden ids are `crate:surface`, so a golden is coupled to a sibling exactly
/// when its owning crate consumes that sibling at runtime. Passing the golden
/// set in keeps this module free of I/O: the caller reads the coupling
/// registry, this decides the relationship.
#[must_use]
pub fn coupled_golden_surfaces<'a>(lib: &str, golden_surface_ids: &[&'a str]) -> Vec<&'a str> {
    let Some(entry) = surface(lib) else {
        return Vec::new();
    };
    let mut coupled: Vec<&'a str> = golden_surface_ids
        .iter()
        .copied()
        .filter(|id| {
            id.split_once(':').is_some_and(|(owner, _)| {
                entry
                    .runtime_consumers
                    .iter()
                    .any(|consumer| *consumer == owner)
            })
        })
        .collect();
    coupled.sort_unstable();
    coupled.dedup();
    coupled
}

/// The registered surface for one sibling.
#[must_use]
pub fn surface(lib: &str) -> Option<&'static SiblingSurface> {
    SURFACES.iter().find(|entry| entry.lib == lib)
}

/// One pinned sibling, as recorded in a lock or observed live.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PinRow {
    /// Sibling library name.
    pub lib: String,
    /// Recorded version string.
    pub version: String,
    /// Recorded git HEAD.
    pub git_head: String,
}

impl PinRow {
    /// Build a pin row.
    #[must_use]
    pub fn new(lib: &str, version: &str, git_head: &str) -> Self {
        Self {
            lib: lib.to_string(),
            version: version.to_string(),
            git_head: git_head.to_string(),
        }
    }
}

/// How one sibling differs between two pin sets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PinMovement {
    /// Identical version and head.
    Unchanged,
    /// Version and/or head differ.
    Moved {
        /// Recorded version.
        from_version: String,
        /// Recorded head.
        from_head: String,
        /// Observed or candidate version.
        to_version: String,
        /// Observed or candidate head.
        to_head: String,
    },
    /// Present in the candidate set only.
    Added,
    /// Present in the recorded set only.
    Removed,
}

impl PinMovement {
    /// Whether this movement changes the pin.
    #[must_use]
    pub const fn is_movement(&self) -> bool {
        !matches!(self, Self::Unchanged)
    }

    /// Stable machine name.
    #[must_use]
    pub const fn slug(&self) -> &'static str {
        match self {
            Self::Unchanged => "unchanged",
            Self::Moved { .. } => "moved",
            Self::Added => "added",
            Self::Removed => "removed",
        }
    }
}

/// One sibling's movement between two pin sets.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PinDelta {
    /// Sibling library name.
    pub lib: String,
    /// How it moved.
    pub movement: PinMovement,
}

impl PinDelta {
    /// A single-line human description naming the sibling and its movement.
    ///
    /// This is the sentence the existing aggregate-hash drift report cannot
    /// produce: an aggregate says only that *something* moved.
    #[must_use]
    pub fn describe(&self) -> String {
        match &self.movement {
            PinMovement::Unchanged => format!("{} unchanged", self.lib),
            PinMovement::Moved {
                from_version,
                from_head,
                to_version,
                to_head,
            } => format!(
                "{} MOVED {from_version}@{from_head} -> {to_version}@{to_head}",
                self.lib
            ),
            PinMovement::Added => format!("{} ADDED (absent from the recorded lock)", self.lib),
            PinMovement::Removed => format!("{} REMOVED (absent from the candidate set)", self.lib),
        }
    }
}

/// Compare a recorded pin set against a live or candidate one.
///
/// The result is sorted by `lib` and covers the union of both sets, so a
/// sibling that appears or disappears is reported rather than silently
/// ignored. Every entry is returned, including unchanged ones, so a caller
/// can render a complete matrix; use [`moved`] for just the offenders.
#[must_use]
pub fn pin_delta(recorded: &[PinRow], candidate: &[PinRow]) -> Vec<PinDelta> {
    let mut libs: Vec<&str> = recorded
        .iter()
        .chain(candidate.iter())
        .map(|row| row.lib.as_str())
        .collect();
    libs.sort_unstable();
    libs.dedup();

    libs.into_iter()
        .map(|lib| {
            let before = recorded.iter().find(|row| row.lib == lib);
            let after = candidate.iter().find(|row| row.lib == lib);
            let movement = match (before, after) {
                (Some(before), Some(after)) => {
                    if before.version == after.version && before.git_head == after.git_head {
                        PinMovement::Unchanged
                    } else {
                        PinMovement::Moved {
                            from_version: before.version.clone(),
                            from_head: before.git_head.clone(),
                            to_version: after.version.clone(),
                            to_head: after.git_head.clone(),
                        }
                    }
                }
                (None, Some(_)) => PinMovement::Added,
                (Some(_), None) => PinMovement::Removed,
                (None, None) => unreachable!("lib came from one of the two sets"),
            };
            PinDelta {
                lib: lib.to_string(),
                movement,
            }
        })
        .collect()
}

/// The subset of a delta that actually moved.
#[must_use]
pub fn moved(deltas: &[PinDelta]) -> Vec<&PinDelta> {
    deltas
        .iter()
        .filter(|delta| delta.movement.is_movement())
        .collect()
}

/// Executed outcome of one sibling's compatibility surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SuiteOutcome {
    /// The surface was not executed. Never a pass.
    NotRun,
    /// The surface was executed against the candidate pin set.
    Executed {
        /// Tests that passed.
        passed: usize,
        /// Tests that failed.
        failed: usize,
    },
}

impl SuiteOutcome {
    /// Whether this outcome is a green execution with at least one test.
    ///
    /// Zero executed tests is NOT green: a selector that matched nothing must
    /// not read as evidence.
    #[must_use]
    pub const fn is_green(self) -> bool {
        matches!(self, Self::Executed { passed, failed } if failed == 0 && passed > 0)
    }
}

/// One sibling's reported suite result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SuiteResult {
    /// Sibling library name.
    pub lib: String,
    /// Reported outcome.
    pub outcome: SuiteOutcome,
}

impl SuiteResult {
    /// Build a suite result.
    #[must_use]
    pub fn new(lib: &str, outcome: SuiteOutcome) -> Self {
        Self {
            lib: lib.to_string(),
            outcome,
        }
    }
}

/// Justification classes for an out-of-train bump.
///
/// Convenience, new upstream features, and version freshness are deliberately
/// absent: an emergency that cannot be named here is not an emergency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmergencyClass {
    /// A reachable security defect.
    SecurityDefect,
    /// Credible data corruption.
    CredibleCorruption,
    /// A false scientific or certificate result.
    FalseScientificResult,
    /// A cancellation or durability contract violation.
    ContractViolation,
    /// A critical sibling became unavailable.
    SiblingUnavailable,
}

impl EmergencyClass {
    /// Stable machine name.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::SecurityDefect => "security-defect",
            Self::CredibleCorruption => "credible-corruption",
            Self::FalseScientificResult => "false-scientific-result",
            Self::ContractViolation => "contract-violation",
            Self::SiblingUnavailable => "sibling-unavailable",
        }
    }
}

/// Whether a bump runs on the train or out of it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BumpKind {
    /// A scheduled release train.
    ReleaseTrain,
    /// An out-of-train bump with its justification class.
    Emergency(EmergencyClass),
}

/// Declared implication for semantic goldens.
///
/// Pins and goldens move together; an undeclared implication is a refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GoldenDisposition {
    /// No moved sibling feeds a semantic golden.
    NoGoldenSurface,
    /// Goldens are affected but verified unchanged, with a justification.
    Unaffected {
        /// Why the goldens did not move.
        justification: String,
    },
    /// Goldens were deliberately re-frozen in the same commit.
    Rebased {
        /// Why the re-freeze is correct.
        justification: String,
    },
}

/// A proposed pin change awaiting adjudication.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BumpAttempt {
    /// Train or emergency.
    pub kind: BumpKind,
    /// Movement between the recorded lock and the candidate set.
    pub deltas: Vec<PinDelta>,
    /// Reported suite results against the candidate set.
    pub results: Vec<SuiteResult>,
    /// Declared golden implication.
    pub golden: GoldenDisposition,
    /// Golden surfaces the caller determined are coupled to the movers.
    ///
    /// Supplied rather than derived so this module performs no I/O. An empty
    /// list with a `NoGoldenSurface` disposition is consistent; a non-empty
    /// list with that disposition is a refusal.
    pub coupled_goldens: Vec<String>,
}

/// Why a bump was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BumpRefusal {
    /// Nothing moved, so there is no bump to adjudicate.
    NoMovement,
    /// A sibling requiring coverage reported no result at all.
    MissingResult {
        /// The uncovered sibling.
        lib: String,
    },
    /// A sibling requiring coverage was not executed.
    NotExecuted {
        /// The unexecuted sibling.
        lib: String,
    },
    /// An executed surface reported failures.
    FailingSurface {
        /// The failing sibling.
        lib: String,
        /// Failing test count.
        failed: usize,
    },
    /// An executed surface ran zero tests, which is not evidence.
    EmptyExecution {
        /// The sibling whose selector matched nothing.
        lib: String,
    },
    /// A moved sibling has no registered compatibility surface.
    UnregisteredSibling {
        /// The unknown sibling.
        lib: String,
    },
    /// Goldens are coupled to a mover but the attempt declared none.
    UndeclaredGoldenImplication {
        /// Coupled golden surfaces the attempt failed to account for.
        surfaces: Vec<String>,
    },
    /// A moved sibling is registered but carries no tests.
    UncoveredSurface {
        /// The sibling with no compatibility tests.
        lib: String,
        /// Why it carries none.
        reason: String,
    },
}

impl fmt::Display for BumpRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoMovement => write!(
                formatter,
                "no sibling moved; there is no pin change to adjudicate"
            ),
            Self::MissingResult { lib } => write!(
                formatter,
                "`{lib}` requires compatibility evidence and none was reported"
            ),
            Self::NotExecuted { lib } => write!(
                formatter,
                "`{lib}` compatibility suite was not run; an unrun suite is never a pass"
            ),
            Self::FailingSurface { lib, failed } => write!(
                formatter,
                "`{lib}` compatibility suite reported {failed} failing test(s)"
            ),
            Self::EmptyExecution { lib } => write!(
                formatter,
                "`{lib}` compatibility suite executed zero tests; a selector that matched nothing \
                 is not evidence"
            ),
            Self::UnregisteredSibling { lib } => write!(
                formatter,
                "`{lib}` moved but has no registered compatibility surface"
            ),
            Self::UndeclaredGoldenImplication { surfaces } => write!(
                formatter,
                "the bump moves siblings feeding {} semantic golden surface(s) ({}) but declares \
                 no golden implication; pins and goldens move together",
                surfaces.len(),
                surfaces.join(", ")
            ),
            Self::UncoveredSurface { lib, reason } => write!(
                formatter,
                "`{lib}` moved but carries no compatibility tests: {reason}"
            ),
        }
    }
}

/// Adjudication of a bump attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BumpVerdict {
    /// The bump may land.
    Admitted {
        /// Siblings that moved.
        moved: Vec<String>,
        /// Surfaces that reported a green execution.
        green_surfaces: Vec<String>,
    },
    /// The bump is refused. Refusal is total.
    Refused {
        /// Every reason, in a stable order.
        reasons: Vec<BumpRefusal>,
    },
}

impl BumpVerdict {
    /// Whether the bump may land.
    #[must_use]
    pub const fn admitted(&self) -> bool {
        matches!(self, Self::Admitted { .. })
    }
}

fn coverage_required(deltas: &[PinDelta], surface: &SiblingSurface) -> bool {
    if surface.priority.required_every_train() {
        return true;
    }
    deltas
        .iter()
        .any(|delta| delta.lib == surface.lib && delta.movement.is_movement())
}

/// Adjudicate a bump attempt fail-closed.
///
/// A bump is admitted only when every sibling requiring coverage reported a
/// green execution with at least one test, and the golden implication is
/// declared. Refusal is total and reports every reason at once, so a caller
/// fixes the whole set rather than rediscovering one problem per attempt.
#[must_use]
pub fn evaluate_bump(attempt: &BumpAttempt) -> BumpVerdict {
    let mut reasons = Vec::new();

    let movers = moved(&attempt.deltas);
    if movers.is_empty() {
        reasons.push(BumpRefusal::NoMovement);
    }

    // A moved sibling nobody registered cannot be adjudicated at all.
    for delta in &movers {
        if surface(&delta.lib).is_none() {
            reasons.push(BumpRefusal::UnregisteredSibling {
                lib: delta.lib.clone(),
            });
        }
    }

    for entry in SURFACES {
        if !coverage_required(&attempt.deltas, entry) {
            continue;
        }
        let sibling_moved = movers.iter().any(|delta| delta.lib == entry.lib);
        if entry.tests.is_empty() {
            // An uncovered surface only blocks when that sibling actually moved;
            // otherwise it is a standing gap recorded by the registry, not a
            // reason to refuse an unrelated bump.
            if sibling_moved {
                reasons.push(BumpRefusal::UncoveredSurface {
                    lib: entry.lib.to_string(),
                    reason: entry
                        .no_test_reason
                        .unwrap_or("no reason recorded")
                        .to_string(),
                });
            }
            continue;
        }
        let Some(result) = attempt
            .results
            .iter()
            .find(|result| result.lib == entry.lib)
        else {
            reasons.push(BumpRefusal::MissingResult {
                lib: entry.lib.to_string(),
            });
            continue;
        };
        match result.outcome {
            SuiteOutcome::NotRun => reasons.push(BumpRefusal::NotExecuted {
                lib: entry.lib.to_string(),
            }),
            SuiteOutcome::Executed { passed: _, failed } if failed != 0 => {
                reasons.push(BumpRefusal::FailingSurface {
                    lib: entry.lib.to_string(),
                    failed,
                });
            }
            SuiteOutcome::Executed { passed: 0, .. } => {
                reasons.push(BumpRefusal::EmptyExecution {
                    lib: entry.lib.to_string(),
                });
            }
            SuiteOutcome::Executed { .. } => {}
        }
    }

    if !attempt.coupled_goldens.is_empty()
        && matches!(attempt.golden, GoldenDisposition::NoGoldenSurface)
    {
        let mut surfaces = attempt.coupled_goldens.clone();
        surfaces.sort();
        surfaces.dedup();
        reasons.push(BumpRefusal::UndeclaredGoldenImplication { surfaces });
    }

    if reasons.is_empty() {
        BumpVerdict::Admitted {
            moved: movers.iter().map(|delta| delta.lib.clone()).collect(),
            green_surfaces: attempt
                .results
                .iter()
                .filter(|result| result.outcome.is_green())
                .map(|result| result.lib.clone())
                .collect(),
        }
    } else {
        BumpVerdict::Refused { reasons }
    }
}

/// Render the compatibility registry as deterministic Markdown.
#[must_use]
pub fn render_registry() -> String {
    let mut out = String::new();
    out.push_str("| sibling | priority | claims | tests | selector |\n");
    out.push_str("| --- | --- | --- | --- | --- |\n");
    for entry in SURFACES {
        let selector = entry
            .selector()
            .unwrap_or_else(|| "NO COVERAGE".to_string());
        let _ = writeln!(
            out,
            "| {} | {} | {} | {} | `{}` |",
            entry.lib,
            entry.priority.slug(),
            entry.claims.len(),
            entry.tests.len(),
            selector
        );
    }
    out
}
