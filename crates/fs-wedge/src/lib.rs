//! fs-wedge — go-to-market wedge selection as data (plan addendum,
//! Proposal 7). Layer: UTIL (pure data + audit; no dependencies).
//!
//! The wedge is the beachhead. The load-bearing DOCTRINE is a NEGATIVE one:
//!
//! > DO NOT SELL AGAINST PEAK SINGLE-PHYSICS FIDELITY ANYWHERE.
//!
//! Unification at the solver level loses to specialized codes on every
//! individual physics; nobody buys a beautifully glued assembly of second-rate
//! solvers, and nobody needs FrankenSim where one mature code owns the whole
//! problem. The wedge must be a WORKFLOW that is today three tools, lossy
//! handoffs, and week-long iteration — where certified seams (the sheaf),
//! incremental re-solve of variants (Proposal 2), and autonomous gradient
//! exploration (Proposal 1) dominate EVEN WITH merely-decent kernels.
//!
//! This crate encodes the selection: the chosen V1 vertical (conjugate heat
//! transfer for electronics cooling) scored on FOUR criteria, the named
//! second/third verticals with their proposal-exercise mapping, and the
//! [`CycleTimeBaseline`] that makes the `>=3×` cycle-time kill criterion
//! MEASURABLE. This changes nothing about what the system can do and everything
//! about whether anyone finds out — hence its own commercial kill criterion,
//! separate from the platform roadmap.

/// The load-bearing negative doctrine of wedge selection.
pub const WEDGE_DOCTRINE: &str = "Do not sell against peak single-physics fidelity anywhere; the wedge is a \
     multi-tool workflow with lossy handoffs where certified seams + incremental \
     re-solve + autonomous gradients win even with merely-decent kernels.";

/// The four criteria a wedge vertical is scored on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WedgeCriterion {
    /// Kernels are individually MATURE AND MODEST (no peak-fidelity arms race);
    /// correlation-based bottom rungs make the fidelity ladder immediately real.
    KernelMaturity,
    /// The cross-team iteration loop is the ACKNOWLEDGED, quantified pain today.
    IterationPain,
    /// ROI is QUANTIFIABLE per design cycle.
    QuantifiableRoi,
    /// Regulatory friction is LOW (the evidence-package story matures on
    /// friendly ground before facing the FAA).
    LowRegulatoryFriction,
}

impl WedgeCriterion {
    /// All four criteria, in order.
    pub const ALL: [WedgeCriterion; 4] = [
        WedgeCriterion::KernelMaturity,
        WedgeCriterion::IterationPain,
        WedgeCriterion::QuantifiableRoi,
        WedgeCriterion::LowRegulatoryFriction,
    ];

    /// A stable slug.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            WedgeCriterion::KernelMaturity => "kernel-maturity",
            WedgeCriterion::IterationPain => "iteration-pain",
            WedgeCriterion::QuantifiableRoi => "quantifiable-roi",
            WedgeCriterion::LowRegulatoryFriction => "low-regulatory-friction",
        }
    }
}

/// A criterion score (`0..=10`) with its rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriterionScore {
    /// Which criterion.
    pub criterion: WedgeCriterion,
    /// The score, `0..=10`.
    pub score: u8,
    /// Why.
    pub rationale: &'static str,
}

const fn s(criterion: WedgeCriterion, score: u8, rationale: &'static str) -> CriterionScore {
    CriterionScore {
        criterion,
        score,
        rationale,
    }
}

/// A candidate vertical with its rank, four-criteria scores, the proposals it
/// exercises, and a rationale.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Vertical {
    /// A stable slug.
    pub name: &'static str,
    /// A human name.
    pub display: &'static str,
    /// Rank (1 = chosen beachhead, then 2, 3).
    pub rank: u8,
    /// The four-criteria scores (in `WedgeCriterion::ALL` order).
    pub scores: [CriterionScore; 4],
    /// The proposals this vertical progressively exercises.
    pub exercises: &'static [&'static str],
    /// Why this vertical, at this rank.
    pub rationale: &'static str,
}

impl Vertical {
    /// This vertical's score for a criterion.
    #[must_use]
    pub fn score(&self, criterion: WedgeCriterion) -> u8 {
        self.scores
            .iter()
            .find(|s| s.criterion == criterion)
            .map_or(0, |s| s.score)
    }

    /// The minimum score across all four criteria (a wedge is only as good as
    /// its weakest criterion).
    #[must_use]
    pub fn weakest_criterion_score(&self) -> u8 {
        self.scores.iter().map(|s| s.score).min().unwrap_or(0)
    }
}

use WedgeCriterion::{IterationPain, KernelMaturity, LowRegulatoryFriction, QuantifiableRoi};

/// The ranked verticals: V1 conjugate heat transfer, then aeroelastic
/// screening, then additive-manufacturing distortion.
const VERTICALS: [Vertical; 3] = [
    Vertical {
        name: "conjugate-heat-transfer",
        display: "Conjugate heat transfer for electronics cooling",
        rank: 1,
        scores: [
            s(
                KernelMaturity,
                8,
                "conduction FEM + forced-convection CFD with correlation-based Nusselt rungs (the fs-ladder cht() bottom rung — makes Proposal 3 real)",
            ),
            s(
                IterationPain,
                9,
                "the thermal<->mechanical/layout iteration loop is the acknowledged pain: today 3 tools, lossy handoffs, week-long cycles",
            ),
            s(
                QuantifiableRoi,
                9,
                "ROI is quantifiable per design cycle (cycle-time reduction directly measurable)",
            ),
            s(
                LowRegulatoryFriction,
                9,
                "low regulatory friction — the evidence-package story matures on friendly ground before the FAA",
            ),
        ],
        exercises: &["2", "1", "3", "12"],
        rationale: "the beachhead: modest mature kernels, acknowledged cross-team pain, quantifiable ROI, friendly regulatory ground",
    },
    Vertical {
        name: "aeroelastic-screening",
        display: "Aeroelastic screening",
        rank: 2,
        scores: [
            s(
                KernelMaturity,
                6,
                "structural + aerodynamic kernels are mature but the coupling is where handoffs hurt",
            ),
            s(
                IterationPain,
                8,
                "flutter/divergence screening iterates across structures and aero teams",
            ),
            s(
                QuantifiableRoi,
                7,
                "ROI via faster screening of the design envelope",
            ),
            s(
                LowRegulatoryFriction,
                5,
                "moderate friction — closer to certification-sensitive aerospace",
            ),
        ],
        exercises: &["1"],
        rationale: "second vertical: progressively exercises Proposal 1 (autonomous gradient exploration across the coupled loop)",
    },
    Vertical {
        name: "additive-manufacturing-distortion",
        display: "Additive-manufacturing distortion",
        rank: 3,
        scores: [
            s(
                KernelMaturity,
                6,
                "thermo-mechanical distortion kernels exist but validation against builds is the pain",
            ),
            s(
                IterationPain,
                8,
                "print-measure-recompensate loops are slow and physical",
            ),
            s(
                QuantifiableRoi,
                7,
                "ROI via fewer scrapped builds / compensation iterations",
            ),
            s(
                LowRegulatoryFriction,
                6,
                "moderate friction depending on the end-use part",
            ),
        ],
        exercises: &["11", "4"],
        rationale: "third vertical: exercises Proposal 11 (reality as another chart — registration against scans) and Proposal 4 (extend the complex into time)",
    },
];

/// The ranked verticals.
#[must_use]
pub fn verticals() -> &'static [Vertical] {
    &VERTICALS
}

/// The four wedge-selection criteria.
#[must_use]
pub fn four_criteria() -> [WedgeCriterion; 4] {
    WedgeCriterion::ALL
}

/// The chosen beachhead (rank 1): conjugate heat transfer.
#[must_use]
pub fn chosen_wedge() -> &'static Vertical {
    VERTICALS
        .iter()
        .find(|v| v.rank == 1)
        .expect("a rank-1 wedge")
}

/// The baseline that makes the cycle-time kill criterion measurable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleTimeBaseline {
    /// Which vertical.
    pub vertical: &'static str,
    /// Today's baseline design-cycle time (days) for the acknowledged loop.
    pub baseline_days: f64,
    /// The cycle-time reduction factor the kill criterion demands (`3.0`).
    pub target_reduction: f64,
    /// The window (quarters after GA) to hit it or re-select the wedge.
    pub kill_within_quarters: u8,
}

impl CycleTimeBaseline {
    /// Does a measured cycle time meet the `>=target_reduction×` kill
    /// criterion? (`baseline / measured >= target_reduction`.)
    #[must_use]
    pub fn meets_kill_criterion(&self, measured_days: f64) -> bool {
        measured_days > 0.0 && self.baseline_days / measured_days >= self.target_reduction
    }
}

/// The conjugate-heat-transfer cycle-time baseline (a week-long loop today).
pub const CHT_BASELINE: CycleTimeBaseline = CycleTimeBaseline {
    vertical: "conjugate-heat-transfer",
    baseline_days: 5.0,
    target_reduction: 3.0,
    kill_within_quarters: 2,
};

/// One named go-to-market audit check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuditCheck {
    /// The check name.
    pub name: &'static str,
    /// Did it pass?
    pub passed: bool,
}

/// The go-to-market audit result.
#[derive(Debug, Clone, PartialEq)]
pub struct WedgeAudit {
    /// The named checks (chosen-strong-on-all, ranks-complete,
    /// all-exercise-proposals, kill-criterion-measurable).
    pub checks: Vec<AuditCheck>,
    /// Any gaps (human-readable).
    pub gaps: Vec<String>,
}

impl WedgeAudit {
    /// Is the go-to-market story complete and self-consistent?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty()
    }

    /// Did a named check pass?
    #[must_use]
    pub fn passed(&self, name: &str) -> bool {
        self.checks.iter().any(|c| c.name == name && c.passed)
    }
}

/// The threshold a chosen wedge must clear on EVERY criterion.
pub const STRONG_THRESHOLD: u8 = 8;

/// Audit the wedge selection: the chosen wedge must be strong on all four
/// criteria, the verticals ranked 1/2/3, every vertical must exercise a
/// proposal, and the kill criterion must be measurable (`>= 3×`).
#[must_use]
pub fn audit() -> WedgeAudit {
    let mut gaps = Vec::new();

    let chosen = chosen_wedge();
    let chosen_strong_on_all = chosen.weakest_criterion_score() >= STRONG_THRESHOLD;
    if !chosen_strong_on_all {
        gaps.push(format!(
            "chosen wedge '{}' is weak on its weakest criterion ({} < {STRONG_THRESHOLD})",
            chosen.name,
            chosen.weakest_criterion_score()
        ));
    }

    let mut ranks: Vec<u8> = VERTICALS.iter().map(|v| v.rank).collect();
    ranks.sort_unstable();
    let ranks_complete = ranks == vec![1, 2, 3];
    if !ranks_complete {
        gaps.push("verticals are not ranked exactly 1, 2, 3".to_string());
    }

    let all_exercise_proposals = VERTICALS.iter().all(|v| !v.exercises.is_empty());
    if !all_exercise_proposals {
        gaps.push("a vertical names no exercised proposal".to_string());
    }

    let kill_criterion_measurable = (CHT_BASELINE.target_reduction - 3.0).abs() < f64::EPSILON;
    if !kill_criterion_measurable {
        gaps.push("cycle-time kill criterion is not the required 3x".to_string());
    }

    WedgeAudit {
        checks: vec![
            AuditCheck {
                name: "chosen-strong-on-all",
                passed: chosen_strong_on_all,
            },
            AuditCheck {
                name: "ranks-complete",
                passed: ranks_complete,
            },
            AuditCheck {
                name: "all-exercise-proposals",
                passed: all_exercise_proposals,
            },
            AuditCheck {
                name: "kill-criterion-measurable",
                passed: kill_criterion_measurable,
            },
        ],
        gaps,
    }
}

/// Emit the wedge selection as deterministic machine-readable JSON.
#[must_use]
pub fn to_json() -> String {
    use core::fmt::Write as _;
    let mut out = String::from("{\"verticals\":[");
    for (i, v) in VERTICALS.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"name\":\"{}\",\"rank\":{},\"weakest_score\":{},\"exercises\":[",
            v.name,
            v.rank,
            v.weakest_criterion_score()
        )
        .expect("write to String");
        for (j, p) in v.exercises.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            write!(out, "\"{p}\"").expect("write");
        }
        out.push_str("]}");
    }
    write!(
        out,
        "],\"baseline_days\":{},\"target_reduction\":{}}}",
        CHT_BASELINE.baseline_days, CHT_BASELINE.target_reduction
    )
    .expect("write");
    out
}
