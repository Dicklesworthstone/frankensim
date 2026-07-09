//! fs-schedule-e2e — CampaignSchedule: certified scheduling of a design
//! campaign, driven by value of information. Layer: L6 (HELM).
//!
//! # The campaign
//!
//! A design campaign asks two orthogonal questions, and this answers BOTH with
//! certificates, composing crates never designed to meet:
//!
//! - **WHEN does it finish?** ([`fs_tropical`]) — the refinement studies form a
//!   precedence DAG with per-study latencies. The completion time is the
//!   longest weighted path — a MAX-PLUS (tropical) critical-path computation,
//!   exact by construction. It also names the bottleneck study (tuning anything
//!   with positive slack cannot move the finish).
//! - **WHETHER to keep going?** ([`fs_voi`]) — the candidate designs have
//!   uncertain performance from numerical / statistical / model sources. The
//!   Expected Value of Perfect Information (EVPI) measures how much the current
//!   ranking ambiguity is worth resolving; `recommend` picks the highest
//!   value-per-cost study, OR says STOP when the decision is already robust
//!   (EVPI below threshold) — the anytime "don't gather data you don't need".
//! - **Honest colors** ([`fs_evidence`]) — the makespan is `Verified` (an exact
//!   tropical computation); the EVPI-driven recommendation is `Estimated` (a
//!   decision-theoretic model).
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::Color;
use fs_tropical::TaskDag;
use fs_voi::{Action, DesignEstimate, Recommendation, evpi, ranking_flip_probability, recommend};

/// One refinement study: a node in the schedule DAG.
#[derive(Debug, Clone)]
pub struct Study {
    /// Study name.
    pub name: String,
    /// Latency (compute cost / wall-clock units).
    pub latency: f64,
    /// Prerequisite study indices (must finish first).
    pub deps: Vec<usize>,
}

impl Study {
    /// A study.
    #[must_use]
    pub fn new(name: impl Into<String>, latency: f64, deps: Vec<usize>) -> Study {
        Study {
            name: name.into(),
            latency,
            deps,
        }
    }
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct ScheduleReport {
    /// The campaign makespan (tropical critical-path length).
    pub makespan: f64,
    /// The makespan's color (`Verified` — exact tropical algebra).
    pub makespan_color: Color,
    /// The critical path (study indices, in order).
    pub critical_path: Vec<usize>,
    /// The bottleneck study name (highest-latency critical study).
    pub bottleneck: Option<String>,
    /// Studies with positive slack (safe to defer) — names.
    pub slack_studies: Vec<String>,
    /// The current Expected Value of Perfect Information (decision ambiguity).
    pub evpi: f64,
    /// The leading design by mean performance.
    pub leading_design: String,
    /// Probability the ranking flips between the top two designs.
    pub flip_risk: f64,
    /// The VoI recommendation: `Act: <study>` or `Stop: <reason>`.
    pub recommendation: String,
    /// Whether the campaign should STOP (decision already robust).
    pub should_stop: bool,
}

/// Run the CampaignSchedule campaign.
///
/// # Panics
/// If `studies` is empty.
#[must_use]
pub fn run_campaign(
    studies: &[Study],
    designs: &[DesignEstimate],
    actions: &[Action],
    stop_threshold: f64,
) -> ScheduleReport {
    assert!(!studies.is_empty(), "need at least one study");

    // --- WHEN: the tropical critical path over the study precedence DAG. ---
    let latencies: Vec<f64> = studies.iter().map(|s| s.latency).collect();
    let mut dag = TaskDag::new(latencies);
    for (v, s) in studies.iter().enumerate() {
        for &u in &s.deps {
            dag = dag.with_edge(u, v);
        }
    }
    let cp = dag.critical_path().expect("the study DAG is acyclic");
    let bottleneck = dag.bottleneck(&cp).map(|i| studies[i].name.clone());
    let slack_studies: Vec<String> = cp
        .slack
        .iter()
        .enumerate()
        .filter(|(_, slack)| **slack > 1e-9)
        .map(|(i, _)| studies[i].name.clone())
        .collect();
    // The makespan is an EXACT max-plus computation → Verified.
    let makespan_color = Color::Verified {
        lo: cp.makespan,
        hi: cp.makespan,
    };

    // --- WHETHER: value of information over the candidate designs. ---
    let current_evpi = evpi(designs);
    // Leading design + the top-two flip risk (the decision's fragility).
    let mut ranked: Vec<&DesignEstimate> = designs.iter().collect();
    ranked.sort_by(|a, b| b.mean.total_cmp(&a.mean));
    let leading_design = ranked.first().map(|d| d.name.clone()).unwrap_or_default();
    let flip_risk = if ranked.len() >= 2 {
        ranking_flip_probability(ranked[0], ranked[1])
    } else {
        0.0
    };
    let (recommendation, should_stop) = match recommend(designs, actions, stop_threshold) {
        Recommendation::Act {
            action,
            value_per_cost,
        } => (
            format!("Act: {action} (value/cost {value_per_cost:.3})"),
            false,
        ),
        Recommendation::Stop { reason } => (format!("Stop: {reason}"), true),
    };

    ScheduleReport {
        makespan: cp.makespan,
        makespan_color,
        critical_path: cp.path,
        bottleneck,
        slack_studies,
        evpi: current_evpi,
        leading_design,
        flip_risk,
        recommendation,
        should_stop,
    }
}
