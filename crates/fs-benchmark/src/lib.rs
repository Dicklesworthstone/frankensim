//! fs-benchmark — the wedge-vertical benchmark & trace corpus (plan addendum,
//! Proposal 7). Layer: UTIL (versioned data + measurement helpers).
//!
//! Governance Rule 2 (doctrine): "a proposal whose kill measurement was never
//! instrumented counts as killed — unmeasured survival is not survival." Many
//! addendum kill criteria measure against "the wedge vertical's benchmark set /
//! recorded traces / merge trials", yet no other bead owns that artifact. THIS
//! is that single, shared, versioned, DETERMINISTIC corpus — so at least six
//! proposals are instrumented instead of killed-by-default.
//!
//! It bundles, for the conjugate-heat-transfer (electronics cooling) wedge:
//! a [`QueryCase`] set (QoIs, tolerances, reference answers with their COLOR
//! and reference cost), an optimizer [`DesignTask`] set, recorded [`EditTrace`]s
//! with known-correct skip sets, an [`MmsCase`] elliptic battery, and swarm
//! [`MergeTrial`]s — plus the measurement helpers ([`speedup`], [`win_rate`],
//! [`conflict_rate`], [`accept_rate`]) each kill-measurement bead uses. Every
//! datum is `const`, so [`corpus_digest`] is bit-stable across runs — the
//! replayability the acceptance criteria demand.

pub use fs_evidence::{Color, ColorRank};

/// The corpus version (measurements are only comparable within a version).
pub const BENCHMARK_VERSION: u32 = 1;

/// One conjugate-heat-transfer design query with a certified reference answer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QueryCase {
    /// Stable id.
    pub id: &'static str,
    /// The quantity of interest.
    pub qoi: &'static str,
    /// The requested tolerance.
    pub tolerance: f64,
    /// The reference (certified) answer.
    pub reference_answer: f64,
    /// The reference compute cost (arbitrary units) to reach it.
    pub reference_cost: f64,
    /// The color CLASS of the reference answer (the full `Color` payload is
    /// materialized by the reference solver; the corpus records its class).
    pub reference_color: ColorRank,
}

/// One optimizer benchmark design task with a known optimum (for adjoint vs
/// derivative-free comparison).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DesignTask {
    /// Stable id.
    pub id: &'static str,
    /// Design dimension.
    pub dimension: usize,
    /// The known optimal objective value.
    pub optimum: f64,
}

/// A recorded agent design-iteration edit trace with a known-correct skip set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EditTrace {
    /// Stable id.
    pub id: &'static str,
    /// Total ops in the DAG for this edit.
    pub total_ops: usize,
    /// Ops that are CERTIFIABLY skippable (the known-correct skip set size).
    pub correct_skips: usize,
}

/// An elliptic manufactured-solution case with a high-precision reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MmsCase {
    /// Stable id.
    pub id: &'static str,
    /// The exact solution sampled at the domain center (reference value).
    pub exact_center: f64,
}

/// A swarm-concurrency merge trial with its harmonic-conflict count.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeTrial {
    /// Stable id.
    pub id: &'static str,
    /// Total merges attempted.
    pub total_merges: usize,
    /// Merges that surfaced a harmonic (structural) conflict.
    pub harmonic_conflicts: usize,
}

const QUERY_SET: [QueryCase; 3] = [
    QueryCase {
        id: "cht-q1",
        qoi: "max junction temperature (K)",
        tolerance: 0.5,
        reference_answer: 358.2,
        reference_cost: 1000.0,
        reference_color: ColorRank::Verified,
    },
    QueryCase {
        id: "cht-q2",
        qoi: "board-to-ambient pressure drop (Pa)",
        tolerance: 5.0,
        reference_answer: 42.0,
        reference_cost: 600.0,
        reference_color: ColorRank::Validated,
    },
    QueryCase {
        id: "cht-q3",
        qoi: "hotspot thermal margin (K)",
        tolerance: 1.0,
        reference_answer: 7.4,
        reference_cost: 250.0,
        reference_color: ColorRank::Estimated,
    },
];

const DESIGN_TASKS: [DesignTask; 3] = [
    DesignTask {
        id: "heatsink-fin-pitch",
        dimension: 4,
        optimum: 351.0,
    },
    DesignTask {
        id: "duct-routing",
        dimension: 8,
        optimum: 38.5,
    },
    DesignTask {
        id: "component-placement",
        dimension: 16,
        optimum: 9.1,
    },
];

const EDIT_TRACES: [EditTrace; 2] = [
    EditTrace {
        id: "raise-fin-count",
        total_ops: 120,
        correct_skips: 96,
    },
    EditTrace {
        id: "move-inlet",
        total_ops: 80,
        correct_skips: 40,
    },
];

const MMS_BATTERY: [MmsCase; 2] = [
    MmsCase {
        id: "poisson-sin",
        exact_center: 1.0,
    },
    MmsCase {
        id: "aniso-diffusion",
        exact_center: 0.25,
    },
];

const MERGE_TRIALS: [MergeTrial; 2] = [
    MergeTrial {
        id: "two-agent-fin-vs-duct",
        total_merges: 40,
        harmonic_conflicts: 6,
    },
    MergeTrial {
        id: "three-agent-layout",
        total_merges: 60,
        harmonic_conflicts: 13,
    },
];

/// The conjugate-heat-transfer query set.
#[must_use]
pub fn query_set() -> &'static [QueryCase] {
    &QUERY_SET
}

/// The optimizer benchmark design tasks.
#[must_use]
pub fn design_tasks() -> &'static [DesignTask] {
    &DESIGN_TASKS
}

/// The recorded edit traces.
#[must_use]
pub fn edit_traces() -> &'static [EditTrace] {
    &EDIT_TRACES
}

/// The elliptic MMS battery.
#[must_use]
pub fn mms_battery() -> &'static [MmsCase] {
    &MMS_BATTERY
}

/// The swarm merge trials.
#[must_use]
pub fn merge_trials() -> &'static [MergeTrial] {
    &MERGE_TRIALS
}

// -- Measurement helpers (each addendum kill measurement) ------------------

/// Speedup `baseline / candidate` (Proposal 8 planner, Proposal 2 skip-yield;
/// kill thresholds are `>= 2×`). Zero if the candidate cost is non-positive.
#[must_use]
pub fn speedup(baseline_cost: f64, candidate_cost: f64) -> f64 {
    if candidate_cost > 0.0 {
        baseline_cost / candidate_cost
    } else {
        0.0
    }
}

/// Win rate: the fraction of `true` outcomes (Proposal 1 adjoint beats
/// derivative-free on `>= 70%` of tasks). Zero for an empty slice.
#[must_use]
pub fn win_rate(outcomes: &[bool]) -> f64 {
    if outcomes.is_empty() {
        return 0.0;
    }
    outcomes.iter().filter(|&&b| b).count() as f64 / outcomes.len() as f64
}

/// A generic rate `count / total` (Proposal 10 harmonic-conflict rate `< 25%`,
/// Proposal 9 accept-rate, Proposal D catch-rate). Zero if `total == 0`.
#[must_use]
pub fn rate(count: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}

/// The harmonic (structural) conflict rate of a merge trial (Proposal 10).
#[must_use]
pub fn conflict_rate(trial: &MergeTrial) -> f64 {
    rate(trial.harmonic_conflicts, trial.total_merges)
}

/// The accept-rate `accepts / attempts` (Proposal 9 certified speculation).
#[must_use]
pub fn accept_rate(accepts: usize, attempts: usize) -> f64 {
    rate(accepts, attempts)
}

// -- Governance Rule 2 discharge -------------------------------------------

/// One proposal whose kill measurement this corpus instruments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InstrumentedProposal {
    /// The proposal id.
    pub proposal: &'static str,
    /// The corpus dataset its kill number is computed against.
    pub dataset: &'static str,
    /// The kill metric.
    pub kill_metric: &'static str,
}

/// The proposals this corpus keeps from being killed-by-default.
const INSTRUMENTED: [InstrumentedProposal; 8] = [
    InstrumentedProposal {
        proposal: "8",
        dataset: "query_set",
        kill_metric: "planner >=2x vs baseline at equal certified accuracy",
    },
    InstrumentedProposal {
        proposal: "1",
        dataset: "design_tasks",
        kill_metric: "adjoint beats derivative-free on >=70% of tasks",
    },
    InstrumentedProposal {
        proposal: "2",
        dataset: "edit_traces",
        kill_metric: "certified skip-yield >=2x vs hash memoization",
    },
    InstrumentedProposal {
        proposal: "D",
        dataset: "design_tasks",
        kill_metric: "guard endpoint catch-rate > random-design catch-rate",
    },
    InstrumentedProposal {
        proposal: "F",
        dataset: "design_tasks",
        kill_metric: "robust optima not dominated by nominal+safety on realized cost",
    },
    InstrumentedProposal {
        proposal: "9",
        dataset: "mms_battery",
        kill_metric: "accept-rate >30% AND warm-start >=1.5x on the elliptic/FEEC battery",
    },
    InstrumentedProposal {
        proposal: "10",
        dataset: "merge_trials",
        kill_metric: "<25% of realistic merges surface harmonic conflicts",
    },
    InstrumentedProposal {
        proposal: "A",
        dataset: "query_set",
        kill_metric: "RB-certified regions cover >=20% of wedge query volume",
    },
];

/// The instrumented proposals (Governance Rule 2 discharge).
#[must_use]
pub fn instrumented_proposals() -> &'static [InstrumentedProposal] {
    &INSTRUMENTED
}

// -- Determinism -----------------------------------------------------------

/// FNV-1a over a canonical serialization of the entire corpus — a bit-stable
/// content digest proving the measurements are replayable.
#[must_use]
pub fn corpus_digest() -> u64 {
    use core::fmt::Write as _;
    let mut buf = format!("v{BENCHMARK_VERSION};");
    for q in &QUERY_SET {
        write!(
            buf,
            "Q{}:{}:{}:{}:{}",
            q.id,
            q.tolerance.to_bits(),
            q.reference_answer.to_bits(),
            q.reference_cost.to_bits(),
            q.reference_color as u8
        )
        .expect("write");
    }
    for d in &DESIGN_TASKS {
        write!(buf, "D{}:{}:{}", d.id, d.dimension, d.optimum.to_bits()).expect("write");
    }
    for e in &EDIT_TRACES {
        write!(buf, "E{}:{}:{}", e.id, e.total_ops, e.correct_skips).expect("write");
    }
    for m in &MMS_BATTERY {
        write!(buf, "M{}:{}", m.id, m.exact_center.to_bits()).expect("write");
    }
    for t in &MERGE_TRIALS {
        write!(buf, "T{}:{}:{}", t.id, t.total_merges, t.harmonic_conflicts).expect("write");
    }
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for b in buf.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

// -- Audit -----------------------------------------------------------------

/// The corpus completeness audit.
#[derive(Debug, Clone, PartialEq)]
pub struct CorpusAudit {
    /// The corpus version.
    pub version: u32,
    /// Number of instrumented proposals.
    pub instrumented: usize,
    /// Gaps (human-readable) — empty when the corpus is complete.
    pub gaps: Vec<String>,
}

impl CorpusAudit {
    /// Is the corpus complete (no gaps)?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }
}

/// Audit the corpus: every dataset non-empty, every query answer colored,
/// every instrumented proposal referencing a real, non-empty dataset.
#[must_use]
pub fn audit() -> CorpusAudit {
    let mut gaps = Vec::new();
    if QUERY_SET.is_empty() {
        gaps.push("query set is empty".to_string());
    }
    if DESIGN_TASKS.is_empty() {
        gaps.push("design tasks are empty".to_string());
    }
    if EDIT_TRACES.is_empty() {
        gaps.push("edit traces are empty".to_string());
    }
    if MMS_BATTERY.is_empty() {
        gaps.push("MMS battery is empty".to_string());
    }
    if MERGE_TRIALS.is_empty() {
        gaps.push("merge trials are empty".to_string());
    }
    // every instrumented proposal must name a real, non-empty dataset.
    for ip in &INSTRUMENTED {
        let n = dataset_len(ip.dataset);
        if n == 0 {
            gaps.push(format!(
                "proposal {} references empty/unknown dataset '{}'",
                ip.proposal, ip.dataset
            ));
        }
    }
    CorpusAudit {
        version: BENCHMARK_VERSION,
        instrumented: INSTRUMENTED.len(),
        gaps,
    }
}

fn dataset_len(name: &str) -> usize {
    match name {
        "query_set" => QUERY_SET.len(),
        "design_tasks" => DESIGN_TASKS.len(),
        "edit_traces" => EDIT_TRACES.len(),
        "mms_battery" => MMS_BATTERY.len(),
        "merge_trials" => MERGE_TRIALS.len(),
        _ => 0,
    }
}
