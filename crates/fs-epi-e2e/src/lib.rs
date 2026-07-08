//! fs-epi-e2e — the epistemic type-system end-to-end suite (plan addendum,
//! Proposal 3's Layer-2 conformance harness). Layer: L6.
//!
//! A runnable battery that exercises Layer 2 AS A WHOLE — not each proposal's
//! unit laws in isolation — and is the artifact of record that the type system
//! "fails safe, not just correct". Five stages, each emitting structured log
//! events (returned as data, not printed):
//!
//! 1. **Laundering** — composition cannot upgrade a color (verified∘estimated →
//!    estimated), and a validated claim out of its regime auto-demotes to
//!    estimated ([`fs_evidence`]).
//! 2. **Falsifier economy** — a certificate class with no falsifier is blocked
//!    by the no-falsifier-no-ship gate, and the consequence×doubt allocator
//!    spends monotonically (cold-start = max doubt).
//! 3. **Goodhart guard** — a discretization-exploit endpoint is REFUSED
//!    (`Failed`), a genuine smooth optimum is honored (`Cleared`), and an
//!    unavailable step leaves the endpoint `Provisional`, never false-cleared
//!    ([`fs_opt`]).
//! 4. **Objective epistemics** — no optimization against an un-colored
//!    objective, the weakest-input rule on headlines, a colored fragility curve
//!    ([`fs_robust`]).
//! 5. **Evidence round-trip** — a package re-verifies through the solver-free
//!    checker, renders its budget pie, and a tampered package fails with a
//!    localized finding ([`fs_package`] → [`fs_checker`]).
//!
//! [`run_battery`] runs all five and returns a structured [`EpiE2eReport`];
//! each fail-closed assertion is the load-bearing check for its stage.

use std::collections::BTreeMap;

use fs_checker::{check, check_against_root};
use fs_evidence::{
    ClaimContext, Color, ColorRank, FalsifierHistory, FalsifierRegistry, IntervalOp,
    ValidityDomain, allocate_budget, check_regime, compose,
};
use fs_opt::{
    DeltaPerturbationStep, Endpoint, EscalationKind, EscalationStep, GoodhartGuard, GuardStatus,
    StepOutcome, converged_and_guard_cleared,
};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_robust::{ColoredObjective, RobustError, fragility_curve, robust_optimum, weakest_color};

/// One stage's structured result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageLog {
    /// The stage name.
    pub stage: &'static str,
    /// Did every fail-closed assertion in the stage hold?
    pub passed: bool,
    /// The structured log events (what each check observed).
    pub events: Vec<String>,
}

/// The full end-to-end report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EpiE2eReport {
    /// The five stage logs, in order.
    pub stages: Vec<StageLog>,
}

impl EpiE2eReport {
    /// Did the whole battery pass?
    #[must_use]
    pub fn passed(&self) -> bool {
        self.stages.iter().all(|s| s.passed)
    }

    /// A named stage, if present.
    #[must_use]
    pub fn stage(&self, name: &str) -> Option<&StageLog> {
        self.stages.iter().find(|s| s.stage == name)
    }
}

/// Run the full Layer-2 battery.
#[must_use]
pub fn run_battery() -> EpiE2eReport {
    EpiE2eReport {
        stages: vec![
            stage_laundering(),
            stage_falsifier(),
            stage_goodhart_guard(),
            stage_objective_epistemics(),
            stage_evidence_roundtrip(),
        ],
    }
}

fn verified() -> Color {
    Color::Verified { lo: -1.0, hi: 1.0 }
}
fn estimated() -> Color {
    Color::Estimated {
        estimator: "surrogate".to_string(),
        dispersion: 2.0,
    }
}
fn cht_regime() -> ValidityDomain {
    ValidityDomain::unconstrained().with("Re", 1e5, 3e5)
}

/// Stage 1: laundering refusal + out-of-regime demotion.
#[must_use]
pub fn stage_laundering() -> StageLog {
    let mut events = Vec::new();
    let mut passed = true;

    // composition takes the MINIMUM rank — an estimate cannot be laundered up.
    let composed = compose(&verified(), &estimated(), IntervalOp::Add);
    let no_upgrade = composed.rank() == ColorRank::Estimated;
    events.push(format!(
        "compose(verified, estimated) -> {} (no laundering: min rank)",
        rank_name(composed.rank())
    ));
    passed &= no_upgrade;

    // a validated claim OUT of its regime auto-demotes to estimated.
    let validated = Color::Validated {
        regime: cht_regime(),
        dataset: "wt-2026".to_string(),
    };
    let mut outside = BTreeMap::new();
    outside.insert("Re".to_string(), 5e5);
    let (demoted, demotion) = check_regime(&validated, &outside);
    let demoted_ok = demoted.rank() == ColorRank::Estimated && demotion.is_some();
    events.push(format!(
        "validated @ Re=5e5 (regime [1e5,3e5]) -> {}, demotion={}",
        rank_name(demoted.rank()),
        demotion.is_some()
    ));
    passed &= demoted_ok;

    // a validated claim INSIDE its regime is preserved (no spurious demotion).
    let validated2 = Color::Validated {
        regime: cht_regime(),
        dataset: "wt-2026".to_string(),
    };
    let mut inside = BTreeMap::new();
    inside.insert("Re".to_string(), 2e5);
    let (kept, no_dem) = check_regime(&validated2, &inside);
    let kept_ok = kept.rank() == ColorRank::Validated && no_dem.is_none();
    events.push(format!(
        "validated @ Re=2e5 (in regime) -> {} (kept)",
        rank_name(kept.rank())
    ));
    passed &= kept_ok;

    StageLog {
        stage: "laundering",
        passed,
        events,
    }
}

/// Stage 2: no-falsifier-no-ship gate + consequence×doubt allocation.
#[must_use]
pub fn stage_falsifier() -> StageLog {
    let mut events = Vec::new();
    let mut passed = true;

    // a class with no registered falsifier is blocked from shipping.
    let registry = FalsifierRegistry::standard();
    let blocked = registry.ship_gate(&["totally-unregistered-class"]);
    let gate_ok = !blocked.is_empty();
    events.push(format!(
        "ship_gate([unregistered]) blocked {blocked:?} (no-falsifier-no-ship)"
    ));
    passed &= gate_ok;

    // the allocator spends monotonically in consequence (cold-start = max doubt).
    let history = FalsifierHistory::new();
    let claims = vec![
        ClaimContext {
            class: "elliptic".to_string(),
            regime: "Re-2e5".to_string(),
            consequence: 10.0,
        },
        ClaimContext {
            class: "elliptic".to_string(),
            regime: "Re-2e5".to_string(),
            consequence: 1.0,
        },
    ];
    let budget = allocate_budget(100.0, &claims, &history);
    let monotone = budget.len() == 2 && budget[0] > budget[1];
    events.push(format!(
        "allocate_budget: high-consequence {:.2} > low-consequence {:.2}",
        budget.first().copied().unwrap_or(0.0),
        budget.get(1).copied().unwrap_or(0.0)
    ));
    passed &= monotone;

    // zero claims -> zero spend (boundary).
    let empty = allocate_budget(100.0, &[], &history);
    let boundary_ok = empty.is_empty();
    events.push(format!(
        "allocate_budget([]) = {} entries (zero spend)",
        empty.len()
    ));
    passed &= boundary_ok;

    StageLog {
        stage: "falsifier",
        passed,
        events,
    }
}

/// A trivially-passing escalation step of a given kind (so the full escalation
/// set can be present — the guard only CLEARS when a step ran for every kind).
struct PassStep(EscalationKind);
impl EscalationStep for PassStep {
    fn kind(&self) -> EscalationKind {
        self.0
    }
    fn evaluate(&self, _endpoint: &Endpoint) -> StepOutcome {
        StepOutcome::Passed
    }
}

/// A full escalation set: the three non-perturbation kinds trivially pass, plus
/// a real δ-perturbation step over `objective`.
fn full_guard<F: Fn(&[f64]) -> f64 + 'static>(
    better_tol: f64,
    sharpness_tol: f64,
    objective: F,
) -> GoodhartGuard {
    GoodhartGuard::new()
        .with_step(Box::new(PassStep(EscalationKind::RungKPlus1)))
        .with_step(Box::new(PassStep(EscalationKind::CrossRepresentation)))
        .with_step(Box::new(PassStep(EscalationKind::EstimatorIndependence)))
        .with_step(Box::new(DeltaPerturbationStep::new(
            0.1,
            better_tol,
            sharpness_tol,
            objective,
        )))
}

/// Stage 3: the Goodhart guard refuses exploits, honors genuine optima, stays
/// provisional when it cannot check.
#[must_use]
pub fn stage_goodhart_guard() -> StageLog {
    let mut events = Vec::new();
    let mut passed = true;

    // EXPLOIT: the endpoint claims objective 0.0, but a δ-perturbation of the
    // honest objective (-x) finds a strictly LOWER value -> the "optimum" is a
    // discretization artifact -> the guard vetoes (Failed, not honored) even
    // though the other escalation steps pass.
    let exploit = full_guard(0.05, 1e9, |x: &[f64]| -x[0]).evaluate(&Endpoint::new(
        "exploit",
        vec![0.0],
        0.0,
    ));
    let exploit_refused = exploit.status == GuardStatus::Failed && !exploit.is_honored();
    events.push(format!(
        "discretization-exploit endpoint -> honored={} (refused)",
        exploit.is_honored()
    ));
    passed &= exploit_refused;

    // SMOOTH: x^2 has a genuine minimum at 0; the full escalation set passes
    // -> Cleared -> honored.
    let smooth = full_guard(0.05, 0.1, |x: &[f64]| x[0] * x[0]).evaluate(&Endpoint::new(
        "smooth",
        vec![0.0],
        0.0,
    ));
    let smooth_honored =
        smooth.status == GuardStatus::Cleared && converged_and_guard_cleared(true, &smooth);
    events.push(format!(
        "genuine smooth optimum -> honored={} (no false veto)",
        smooth.is_honored()
    ));
    passed &= smooth_honored;

    // UNAVAILABLE: a guard with no steps cannot check anything -> Provisional,
    // NEVER false-cleared.
    let bare = GoodhartGuard::new().evaluate(&Endpoint::new("unchecked", vec![0.0], 0.0));
    let provisional = bare.status == GuardStatus::Provisional && !bare.is_honored();
    events.push(format!(
        "unavailable step -> provisional (honored={})",
        bare.is_honored()
    ));
    passed &= provisional;

    StageLog {
        stage: "goodhart-guard",
        passed,
        events,
    }
}

/// Stage 4: objective epistemics — the contract, the weakest-input rule, a
/// colored fragility curve.
#[must_use]
pub fn stage_objective_epistemics() -> StageLog {
    let mut events = Vec::new();
    let mut passed = true;

    // no optimization against an un-colored objective.
    let uncolored = ColoredObjective::new("fiction", vec![1.0, 2.0, 3.0], vec![]);
    let refused = matches!(
        robust_optimum(std::slice::from_ref(&uncolored), 0.9),
        Err(RobustError::UncoloredObjective { .. })
    );
    events.push(format!("robust_optimum(un-colored) refused = {refused}"));
    passed &= refused;

    // the weakest input colors the headline.
    let headline = weakest_color(&[verified(), estimated()]).map(|c| c.rank());
    let weakest_ok = headline == Some(ColorRank::Estimated);
    events.push(format!(
        "weakest_color(verified, estimated) -> {:?}",
        headline.map(rank_name)
    ));
    passed &= weakest_ok;

    // the seismic deliverable is a colored, monotone fragility curve.
    let frag = fragility_curve(&[3.0, 4.0, 5.0, 6.0, 7.0], &[1.0, 4.0, 9.0], estimated());
    let frag_ok = frag.as_ref().is_ok_and(|f| {
        f.color.rank() == ColorRank::Estimated
            && f.curve.first().map(|p| p.prob_failure) <= f.curve.last().map(|p| p.prob_failure)
    });
    events.push(format!("colored fragility curve ok = {frag_ok}"));
    passed &= frag_ok;

    StageLog {
        stage: "objective-epistemics",
        passed,
        events,
    }
}

/// Stage 5: evidence package round-trip through the solver-free checker.
#[must_use]
pub fn stage_evidence_roundtrip() -> StageLog {
    let mut events = Vec::new();
    let mut passed = true;

    let pkg = EvidencePackage::new(Provenance::new("commit-abc", "lock-def"))
        .with_claim(Claim::new("c1", "stress <= sigma*", verified()))
        .with_claim(Claim::new("c2", "surrogate says ok", estimated()));
    let root = pkg.merkle_root();

    // solver-free re-verification passes.
    let good = check(&pkg);
    events.push(format!(
        "checker re-verify (no solver) passed = {}",
        good.passed()
    ));
    passed &= good.passed();

    // a tampered package fails with a LOCALIZED finding.
    let tampered = check_against_root(&pkg, root ^ 0xdead_beef);
    let tamper_caught = !tampered.passed()
        && tampered
            .findings
            .iter()
            .any(|f| f.kind == "content-address-mismatch");
    events.push(format!(
        "tampered package caught = {tamper_caught} (localized)"
    ));
    passed &= tamper_caught;

    // the budget pie renders both color classes present.
    let pie = good.render_pie();
    let pie_ok = pie.contains("verified") && pie.contains("estimated");
    events.push(format!("budget pie rendered = {pie_ok}"));
    passed &= pie_ok;

    StageLog {
        stage: "evidence-roundtrip",
        passed,
        events,
    }
}

fn rank_name(r: ColorRank) -> &'static str {
    match r {
        ColorRank::Verified => "verified",
        ColorRank::Validated => "validated",
        ColorRank::Estimated => "estimated",
    }
}
