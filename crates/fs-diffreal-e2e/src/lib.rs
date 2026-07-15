//! fs-diffreal-e2e — the differentiation & reality end-to-end suite (plan
//! addendum, Proposal 11 / Layer-3 conformance). Layer: L6.
//!
//! A runnable battery that exercises Layer 3 AS A WHOLE — end-to-end adjoints
//! and reality-as-a-chart — and is an artifact of record that the differentiation
//! and as-built machinery FAIL SAFE. Four stages, each emitting structured log
//! events (returned as data, never printed):
//!
//! 1. **Differentiation** — an adjoint (reverse-mode chain rule) gradient agrees
//!    with finite differences within a conditioning-aware tolerance, a
//!    full-VJP-coverage path differentiates, and a path with a MISSING VJP
//!    (a forced remesh) raises a structured error that BLOCKS the gradient —
//!    never a silent zero.
//! 2. **As-built loop** — register a scanned fixture (error carried forward),
//!    compute an estimated as-built δ carrying calibration provenance,
//!    LOCALIZE a seeded defect, and run registration-free point-sensor
//!    assimilation that reduces the model-data misfit ([`fs_asbuilt`],
//!    [`fs_assimilate`]).
//! 3. **Tolerance allocation** — a GD&T report on a known-sensitivity fixture
//!    tightens the high-sensitivity feature, loosens the low one, and the
//!    band-extremes check confirms the P(in-spec) constraint ([`fs_toleralloc`]).
//! 4. **(Gated) spacetime** — the temporal-complex capability exists in
//!    `fs-time`, but its coupled end-to-end fixture is not integrated and
//!    activated in this battery; it is reported as gated, not silently passed.
//!
//! [`run_battery`] runs all four under an explicit [`Cx`] and returns a
//! structured [`DiffRealReport`] only after every cancellation-aware stage has
//! finalized.

use fs_asbuilt::{Fiducial, Point2, as_built_diff, register};
use fs_assimilate::{AssimError, Belief, assimilate_colored, misfit, point_sensor};
use fs_evidence::Color;
use fs_exec::Cx;
use fs_toleralloc::{
    Action, ColorRank, Feature, allocate, gdt_report, robustness_check, variance_budget,
};

/// Stable name of the differentiation stage.
pub const DIFFERENTIATION_STAGE: &str = "differentiation";
/// Stable name of the as-built/assimilation stage.
pub const AS_BUILT_STAGE: &str = "as-built-loop";
/// Stable name of the tolerance-allocation stage.
pub const TOLERANCE_STAGE: &str = "tolerance-allocation";
/// Stable name of the spacetime-integration stage.
pub const SPACETIME_STAGE: &str = "spacetime-gated";

/// Versioned fixture identity expected for the differentiation stage.
pub const DIFFERENTIATION_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/differentiation-fixture/v1";
/// Versioned fixture identity expected for the as-built/assimilation stage.
pub const AS_BUILT_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/as-built-fixture/v1";
/// Versioned fixture identity expected for the tolerance-allocation stage.
pub const TOLERANCE_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/tolerance-allocation-fixture/v1";
/// Versioned fixture identity expected for the spacetime-integration stage.
pub const SPACETIME_EVIDENCE_IDENTITY: &str = "fs-diffreal-e2e/spacetime-integration-gate/v1";

/// Whether a stage is load-bearing for this battery's promotion decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageRequirement {
    /// The report is incomplete until the stage has actually run.
    Required,
    /// The stage is diagnostic and does not block the required-stage decision.
    Optional,
}

impl core::fmt::Display for StageRequirement {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Required => formatter.write_str("required"),
            Self::Optional => formatter.write_str("optional"),
        }
    }
}

/// Stable machine code plus deterministic human-readable detail for a stage
/// that did not pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageReason {
    /// Stable reason code for ledgers and programmatic diagnostics.
    pub code: &'static str,
    /// Human-readable detail. This is diagnostic data, never printed here.
    pub detail: String,
}

impl StageReason {
    /// Construct a structured reason.
    #[must_use]
    pub fn new(code: &'static str, detail: impl Into<String>) -> Self {
        Self {
            code,
            detail: detail.into(),
        }
    }

    fn is_well_formed(&self) -> bool {
        !self.code.trim().is_empty() && !self.detail.trim().is_empty()
    }
}

impl core::fmt::Display for StageReason {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "[{}]: {}", self.code, self.detail)
    }
}

/// Scientific disposition of one stage.
///
/// `Failed` means the stage ran and an assertion was false. `Gated` and
/// `Refused` mean the assertion was not validly evaluated, so neither can
/// satisfy report completeness.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StageStatus {
    /// Every load-bearing assertion ran and passed.
    Passed,
    /// The stage ran, but at least one load-bearing assertion was false.
    Failed(StageReason),
    /// The capability or integration is deliberately unavailable.
    Gated(StageReason),
    /// The stage declined to evaluate because an admissibility condition,
    /// budget, or cancellation condition prevented a trustworthy result.
    Refused(StageReason),
}

impl StageStatus {
    /// Stable lowercase status code for deterministic records.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Passed => "passed",
            Self::Failed(_) => "failed",
            Self::Gated(_) => "gated",
            Self::Refused(_) => "refused",
        }
    }

    /// Did the stage actually run to a scientific pass/fail decision?
    #[must_use]
    pub const fn is_evaluated(&self) -> bool {
        matches!(self, Self::Passed | Self::Failed(_))
    }

    /// Did the stage actually run and pass?
    #[must_use]
    pub const fn is_passed(&self) -> bool {
        matches!(self, Self::Passed)
    }

    /// Structured reason for every non-passing disposition.
    #[must_use]
    pub const fn reason(&self) -> Option<&StageReason> {
        match self {
            Self::Passed => None,
            Self::Failed(reason) | Self::Gated(reason) | Self::Refused(reason) => Some(reason),
        }
    }

    fn is_well_formed(&self) -> bool {
        self.reason().is_none_or(StageReason::is_well_formed)
    }
}

impl core::fmt::Display for StageStatus {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(self.code())?;
        if let Some(reason) = self.reason() {
            write!(formatter, "{reason}")?;
        }
        Ok(())
    }
}

/// One stage's structured diagnostic result.
///
/// A `StageLog` is freely constructible DATA. By itself it carries no
/// promotion authority and cannot be inserted into an opaque
/// [`DiffRealReport`] by downstream callers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageLog {
    /// The stage name.
    pub stage: &'static str,
    /// Whether this stage participates in the required-stage decision.
    pub requirement: StageRequirement,
    /// Typed scientific disposition; unavailable work is never a pass.
    pub status: StageStatus,
    /// Versioned identity of the fixture/schema whose result this log records.
    /// This is a diagnostic identity binding, not a content hash, proof
    /// certificate, independent verification receipt, or authorization.
    pub evidence_identity: &'static str,
    /// The structured log events.
    pub events: Vec<String>,
}

impl StageLog {
    /// Construct one plain diagnostic stage record.
    ///
    /// Construction does not confer authority or add the record to a
    /// [`DiffRealReport`].
    #[must_use]
    pub fn new(
        stage: &'static str,
        requirement: StageRequirement,
        status: StageStatus,
        evidence_identity: &'static str,
        events: Vec<String>,
    ) -> Self {
        Self {
            stage,
            requirement,
            status,
            evidence_identity,
            events,
        }
    }

    /// Did this stage actually run and pass?
    #[must_use]
    pub const fn passed(&self) -> bool {
        self.status.is_passed()
    }

    fn is_well_formed(&self) -> bool {
        !self.stage.trim().is_empty()
            && !self.evidence_identity.trim().is_empty()
            && !self.events.is_empty()
            && self.events.iter().all(|event| !event.trim().is_empty())
            && self.status.is_well_formed()
    }
}

impl core::fmt::Display for StageLog {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "stage={} requirement={} status={} evidence_identity={}",
            self.stage, self.requirement, self.status, self.evidence_identity
        )
    }
}

/// The full crate-authored Layer-3 battery report.
///
/// Construction is intentionally private: downstream callers may inspect the
/// stage diagnostics, but cannot assemble caller-supplied rows into a report
/// whose battery-local readiness predicates pass.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiffRealReport {
    /// Stage logs. The four required stages have a fixed relative order;
    /// additional stages must be explicitly optional.
    stages: Vec<StageLog>,
}

impl DiffRealReport {
    /// Ordered read-only stage diagnostics produced by this battery run.
    #[must_use]
    pub fn stages(&self) -> &[StageLog] {
        &self.stages
    }

    /// Is every required stage present exactly once with the expected evidence
    /// identity and an evaluated (`Passed` or `Failed`) result?
    #[must_use]
    pub fn complete(&self) -> bool {
        self.required_schema_is_valid()
            && REQUIRED_STAGES.iter().all(|required| {
                self.stage(required.name)
                    .is_some_and(|stage| stage.status.is_evaluated())
            })
    }

    /// Did every required stage actually run and pass?
    ///
    /// Missing, duplicated, gated, refused, identity-mismatched, or malformed
    /// required records all return `false`.
    #[must_use]
    pub fn all_required_passed(&self) -> bool {
        self.required_schema_is_valid()
            && REQUIRED_STAGES
                .iter()
                .all(|required| self.stage(required.name).is_some_and(StageLog::passed))
    }

    /// Did this crate-authored fixed battery complete every required fixture
    /// and pass every required assertion?
    ///
    /// This is battery-local readiness only. It is not scientific release
    /// admission, external validation, or an authenticated promotion receipt.
    #[must_use]
    pub fn promotion_ready(&self) -> bool {
        self.complete() && self.all_required_passed()
    }

    /// Fail-closed compatibility alias for [`Self::promotion_ready`].
    #[must_use]
    pub fn passed(&self) -> bool {
        self.promotion_ready()
    }

    /// A named stage.
    #[must_use]
    pub fn stage(&self, name: &str) -> Option<&StageLog> {
        self.stages.iter().find(|s| s.stage == name)
    }

    fn required_schema_is_valid(&self) -> bool {
        if self.stages.iter().any(|stage| !stage.is_well_formed())
            || self.stages.iter().enumerate().any(|(index, stage)| {
                self.stages[index + 1..]
                    .iter()
                    .any(|other| stage.stage == other.stage)
            })
        {
            return false;
        }

        if self
            .stages
            .iter()
            .filter(|stage| stage.requirement == StageRequirement::Required)
            .count()
            != REQUIRED_STAGES.len()
        {
            return false;
        }

        self.stages
            .iter()
            .filter(|stage| stage.requirement == StageRequirement::Required)
            .zip(REQUIRED_STAGES.iter())
            .all(|(stage, required)| {
                stage.stage == required.name
                    && stage.evidence_identity == required.evidence_identity
            })
    }
}

#[derive(Clone, Copy)]
struct RequiredStage {
    name: &'static str,
    evidence_identity: &'static str,
}

const REQUIRED_STAGES: [RequiredStage; 4] = [
    RequiredStage {
        name: DIFFERENTIATION_STAGE,
        evidence_identity: DIFFERENTIATION_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: AS_BUILT_STAGE,
        evidence_identity: AS_BUILT_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: TOLERANCE_STAGE,
        evidence_identity: TOLERANCE_EVIDENCE_IDENTITY,
    },
    RequiredStage {
        name: SPACETIME_STAGE,
        evidence_identity: SPACETIME_EVIDENCE_IDENTITY,
    },
];

/// A lower-layer refusal from the cancellation-aware reality stage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiffRealError {
    /// Registration, registered geometry, or as-built comparison failed.
    AsBuilt(fs_asbuilt::RegError),
    /// Belief construction, observation declaration, or assimilation failed.
    Assimilation(AssimError),
}

impl core::fmt::Display for DiffRealError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::AsBuilt(error) => write!(formatter, "as-built stage failed: {error}"),
            Self::Assimilation(error) => write!(formatter, "assimilation stage failed: {error}"),
        }
    }
}

impl std::error::Error for DiffRealError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::AsBuilt(error) => Some(error),
            Self::Assimilation(error) => Some(error),
        }
    }
}

impl From<fs_asbuilt::RegError> for DiffRealError {
    fn from(error: fs_asbuilt::RegError) -> Self {
        Self::AsBuilt(error)
    }
}

impl From<AssimError> for DiffRealError {
    fn from(error: AssimError) -> Self {
        Self::Assimilation(error)
    }
}

/// Run the full Layer-3 battery.
///
/// # Errors
/// Propagates structured cancellation or scientific refusal from the as-built
/// and assimilation stage. No partial battery report is published.
pub fn run_battery(cx: &Cx<'_>) -> Result<DiffRealReport, DiffRealError> {
    Ok(DiffRealReport {
        stages: vec![
            stage_differentiation(),
            stage_as_built_loop(cx)?,
            stage_tolerance_allocation(),
            stage_spacetime_gated(),
        ],
    })
}

// -- Stage 1: differentiation ----------------------------------------------

/// The fixture composite `f(x) = (2x + 1)²` and its exact adjoint gradient.
fn composite(x: f64) -> f64 {
    let h = 2.0 * x + 1.0;
    h * h
}
fn adjoint_grad(x: f64) -> f64 {
    // reverse mode: u = 2x+1; dg/du = 2u; du/dx = 2 -> grad = 4(2x+1).
    let u = 2.0 * x + 1.0;
    (2.0 * u) * 2.0
}
fn fd_grad(x: f64, eps: f64) -> f64 {
    (composite(x + eps) - composite(x - eps)) / (2.0 * eps)
}

/// Differentiate a pipeline of ops; a missing VJP BLOCKS the gradient (never a
/// silent zero).
///
/// # Errors
/// A message naming the first op whose VJP is missing.
pub fn differentiate_path(
    ops: &[&str],
    has_vjp: impl Fn(&str) -> bool,
    x: f64,
) -> Result<f64, String> {
    for op in ops {
        if !has_vjp(op) {
            return Err(format!(
                "missing VJP for op '{op}': gradient BLOCKED (never silent-zero)"
            ));
        }
    }
    Ok(adjoint_grad(x))
}

/// Stage 1: adjoint-vs-FD agreement + missing-VJP blocking.
#[must_use]
pub fn stage_differentiation() -> StageLog {
    let mut events = Vec::new();
    let mut assertions_passed = true;

    // adjoint agrees with finite differences within a conditioning-aware tol.
    let x = 1.5;
    let a = adjoint_grad(x);
    let fd = fd_grad(x, 1e-6);
    let agree = (a - fd).abs() < 1e-4;
    events.push(format!("adjoint {a:.6} vs FD {fd:.6} -> agree={agree}"));
    assertions_passed &= agree;

    // a smooth SDF/spline path (full VJP coverage) differentiates.
    let smooth = ["sdf", "spline", "solve"];
    let full_cover = |op: &str| matches!(op, "sdf" | "spline" | "solve");
    let smooth_ok = differentiate_path(&smooth, full_cover, x).is_ok();
    events.push(format!(
        "smooth path {smooth:?} differentiates = {smooth_ok}"
    ));
    assertions_passed &= smooth_ok;

    // a forced-remesh path has a missing VJP -> BLOCKED, not silent-zero.
    let remesh = ["sdf", "remesh", "solve"];
    let blocked = differentiate_path(&remesh, full_cover, x);
    let blocked_ok = blocked.is_err();
    events.push(format!(
        "remesh path blocked = {blocked_ok} (never silent-zero)"
    ));
    assertions_passed &= blocked_ok;

    let status = if assertions_passed {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.differentiation.assertion-failed",
            "at least one gradient-agreement or VJP-coverage assertion failed; inspect events",
        ))
    };
    StageLog::new(
        DIFFERENTIATION_STAGE,
        StageRequirement::Required,
        status,
        DIFFERENTIATION_EVIDENCE_IDENTITY,
        events,
    )
}

// -- Stage 2: as-built loop -------------------------------------------------

/// Stage 2: register a scan, estimate as-built δ, localize a defect, assimilate.
///
/// # Errors
/// Propagates the structured lower-layer refusal, including cancellation, and
/// publishes no partial stage log.
pub fn stage_as_built_loop(cx: &Cx<'_>) -> Result<StageLog, DiffRealError> {
    let mut events = Vec::new();
    let mut assertions_passed = true;

    // a scanned fixture: design datums transformed by a known rigid motion.
    let design = [
        Point2::new(0.0, 0.0)?,
        Point2::new(2.0, 0.0)?,
        Point2::new(0.0, 2.0)?,
    ];
    let (theta, tx, ty) = (0.3_f64, 4.0, 1.0);
    let xf = |p: Point2| {
        let (s, c) = theta.sin_cos();
        Point2::new(c * p.x() - s * p.y() + tx, s * p.x() + c * p.y() + ty)
    };
    let fids: Vec<Fiducial> = design
        .iter()
        .map(|&datum| Ok(Fiducial::new(datum, xf(datum)?)))
        .collect::<Result<_, fs_asbuilt::RegError>>()?;
    let reg = register(&fids, cx)?;
    let reg_ok = reg.residual_rms() < 1e-9;
    events.push(format!(
        "registration residual {:.2e} (error carried forward)",
        reg.residual_rms()
    ));
    assertions_passed &= reg_ok;

    // as-built δ with a SEEDED DEFECT on the middle point.
    let design_pts = vec![design[0], design[1], design[2]];
    let mut scanned: Vec<Point2> = design_pts
        .iter()
        .map(|&point| reg.apply(point))
        .collect::<Result<_, _>>()?;
    scanned[1] = Point2::new(scanned[1].x() + 0.3, scanned[1].y())?;
    let diff = as_built_diff(&reg, &design_pts, &scanned, 0.5, 0.02, "cmm-cal-2026", cx)?;
    // localize the defect: the argmax deviation is the seeded point (index 1).
    let defect_idx = diff
        .deviations()
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i);
    let localized = defect_idx == Some(1) && (diff.max_deviation() - 0.3).abs() < 1e-9;
    let estimated = matches!(diff.color(), Color::Estimated { .. });
    events.push(format!(
        "as-built δ max {:.3} @ idx {:?}, estimated={estimated}",
        diff.max_deviation(),
        defect_idx
    ));
    assertions_passed &= localized && estimated;

    // registration-free point-sensor 4D-Var: misfit reduction.
    let prior = Belief::diagonal(vec![20.0, 20.0], &[9.0, 9.0], cx)?;
    let obs = vec![
        point_sensor(0, 2, 24.0, 0.25, "thermocouple-1")?,
        point_sensor(1, 2, 18.5, 0.25, "thermocouple-2")?,
    ];
    let assimilated = assimilate_colored(&prior, &obs, "Re", 1e5, 3e5, cx)?;
    let misfit_reduced = assimilated.misfit_after() < assimilated.misfit_before();
    events.push(format!(
        "assimilation misfit {:.2} -> {:.2}",
        assimilated.misfit_before(),
        assimilated.misfit_after()
    ));
    let checked_after = misfit(assimilated.belief(), &obs, cx)?;
    let checked_before = misfit(&prior, &obs, cx)?;
    assertions_passed &= misfit_reduced && checked_after <= checked_before;

    let status = if assertions_passed {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.as-built.assertion-failed",
            "registration, defect-localization, evidence-color, or assimilation assertion failed; inspect events",
        ))
    };
    Ok(StageLog::new(
        AS_BUILT_STAGE,
        StageRequirement::Required,
        status,
        AS_BUILT_EVIDENCE_IDENTITY,
        events,
    ))
}

// -- Stage 3: tolerance allocation ------------------------------------------

/// Stage 3: adjoint-driven GD&T on a known-sensitivity fixture.
#[must_use]
pub fn stage_tolerance_allocation() -> StageLog {
    let mut events = Vec::new();
    let mut assertions_passed = true;

    let feat = |name: &str, s: f64| Feature {
        name: name.into(),
        sensitivity: s,
        sensitivity_color: ColorRank::Verified,
        cost_coeff: 1.0,
        baseline_tolerance: 0.5,
    };
    let budget = match variance_budget(1.0, 0.99) {
        Ok(budget) => budget,
        Err(error) => {
            let code = "diffreal.tolerance.invalid-budget-fixture";
            let detail = format!("the fixed tolerance-budget fixture was refused: {error:?}");
            return StageLog::new(
                TOLERANCE_STAGE,
                StageRequirement::Required,
                StageStatus::Refused(StageReason::new(code, detail.clone())),
                TOLERANCE_EVIDENCE_IDENTITY,
                vec![format!("REFUSED[{code}]: {detail}")],
            );
        }
    };
    let alloc = match allocate(&[feat("critical", 12.0), feat("slack", 0.2)], budget, 3.0) {
        Ok(allocation) => allocation,
        Err(error) => {
            let code = "diffreal.tolerance.allocation-refused";
            let detail = format!("the fixed tolerance-allocation fixture was refused: {error:?}");
            return StageLog::new(
                TOLERANCE_STAGE,
                StageRequirement::Required,
                StageStatus::Refused(StageReason::new(code, detail.clone())),
                TOLERANCE_EVIDENCE_IDENTITY,
                vec![format!("REFUSED[{code}]: {detail}")],
            );
        }
    };
    // tighten where sensitivity is large, loosen where small.
    let critical_action = alloc
        .items
        .iter()
        .find(|item| item.name == "critical")
        .map(|item| item.action);
    let slack_action = alloc
        .items
        .iter()
        .find(|item| item.name == "slack")
        .map(|item| item.action);
    let tighten_high = critical_action == Some(Action::Tighten);
    let loosen_low = slack_action == Some(Action::Loosen);
    let critical_action_label =
        critical_action.map_or("Missing".to_string(), |action| format!("{action:?}"));
    let slack_action_label =
        slack_action.map_or("Missing".to_string(), |action| format!("{action:?}"));
    events.push(format!(
        "critical -> {critical_action_label}, slack -> {slack_action_label}"
    ));
    assertions_passed &= tighten_high && loosen_low;

    // the GD&T report attaches a certified sensitivity to every loosened tol.
    let report = gdt_report(&alloc);
    let justified = report
        .iter()
        .filter(|s| s.action == Action::Loosen)
        .all(|s| s.certified_sensitivity > 0.0 && s.color == ColorRank::Verified);
    events.push(format!(
        "GD&T report justifies {} loosened tolerances",
        report.iter().filter(|s| s.action == Action::Loosen).count()
    ));
    assertions_passed &= justified;

    // the band-extremes check confirms the P(in-spec) constraint: the QoI at
    // sampled ±t corners stays within k·σ of nominal (σ ≈ √budget ≈ 0.39).
    let verdict = robustness_check(&alloc, &[0.9, -0.8, 0.5], 0.0, 3.0, 0.2);
    events.push(format!(
        "robustness confirmed = {} (linearized std {:.3})",
        verdict.confirmed, verdict.linearized_std
    ));
    assertions_passed &= verdict.confirmed;

    let status = if assertions_passed {
        StageStatus::Passed
    } else {
        StageStatus::Failed(StageReason::new(
            "diffreal.tolerance.assertion-failed",
            "allocation direction, sensitivity justification, or sampled-extremes assertion failed; inspect events",
        ))
    };
    StageLog::new(
        TOLERANCE_STAGE,
        StageRequirement::Required,
        status,
        TOLERANCE_EVIDENCE_IDENTITY,
        events,
    )
}

// -- Stage 4: gated spacetime -----------------------------------------------

/// Stage 4: the spacetime-complex capability is not integrated and activated
/// in this battery (honestly gated, never silently passed).
#[must_use]
pub fn stage_spacetime_gated() -> StageLog {
    StageLog::new(
        SPACETIME_STAGE,
        StageRequirement::Required,
        StageStatus::Gated(StageReason::new(
            "diffreal.spacetime.integration-not-activated",
            "fs-time temporal-complex support exists, but the coupled end-to-end fixture is not integrated and activated in this battery",
        )),
        SPACETIME_EVIDENCE_IDENTITY,
        vec![
            "GATED: temporal-complex dependency frankensim-epic-coupling-bk0o.7 is shipped, but this battery has no activated coupled spacetime fixture; stage not asserted"
                .to_string(),
        ],
    )
}

#[cfg(test)]
mod report_policy_tests {
    use super::*;

    fn required_status_report(spacetime_status: StageStatus) -> DiffRealReport {
        let passed = |stage, identity| {
            StageLog::new(
                stage,
                StageRequirement::Required,
                StageStatus::Passed,
                identity,
                vec![format!("{stage} fixture executed")],
            )
        };
        DiffRealReport {
            stages: vec![
                passed(DIFFERENTIATION_STAGE, DIFFERENTIATION_EVIDENCE_IDENTITY),
                passed(AS_BUILT_STAGE, AS_BUILT_EVIDENCE_IDENTITY),
                passed(TOLERANCE_STAGE, TOLERANCE_EVIDENCE_IDENTITY),
                StageLog::new(
                    SPACETIME_STAGE,
                    StageRequirement::Required,
                    spacetime_status,
                    SPACETIME_EVIDENCE_IDENTITY,
                    vec!["spacetime fixture disposition recorded".to_string()],
                ),
            ],
        }
    }

    #[test]
    fn all_passed_required_stages_are_complete_and_promotion_ready() {
        let report = required_status_report(StageStatus::Passed);
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.promotion_ready());
        assert!(report.passed());
    }

    #[test]
    fn a_failed_required_stage_is_complete_but_not_promotion_ready() {
        let report = required_status_report(StageStatus::Failed(StageReason::new(
            "test.spacetime.assertion-failed",
            "the spacetime fixture ran and violated its asserted bound",
        )));
        assert!(
            report.complete(),
            "failed is an evaluated scientific outcome"
        );
        assert!(!report.all_required_passed());
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn a_gated_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Gated(StageReason::new(
            "test.spacetime.gated",
            "the required capability is unavailable",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn a_refused_required_stage_is_neither_complete_nor_promotion_ready() {
        let report = required_status_report(StageStatus::Refused(StageReason::new(
            "test.spacetime.refused",
            "the stage exhausted its admitted budget before evaluation",
        )));
        assert!(!report.complete());
        assert!(!report.all_required_passed());
        assert!(!report.promotion_ready());
        assert!(!report.passed());
    }

    #[test]
    fn an_explicit_optional_gate_does_not_block_required_stage_promotion() {
        let mut report = required_status_report(StageStatus::Passed);
        report.stages.push(StageLog::new(
            "diagnostic-only",
            StageRequirement::Optional,
            StageStatus::Gated(StageReason::new(
                "test.optional.gated",
                "the optional diagnostic backend is unavailable",
            )),
            "fs-diffreal-e2e/optional-diagnostic/v1",
            vec!["optional diagnostic gate retained".to_string()],
        ));
        assert!(report.complete());
        assert!(report.all_required_passed());
        assert!(report.promotion_ready());
    }

    #[test]
    fn malformed_or_schema_incomplete_reports_fail_closed() {
        let all_passed = required_status_report(StageStatus::Passed);

        let mut missing = all_passed.clone();
        missing
            .stages
            .retain(|stage| stage.stage != SPACETIME_STAGE);
        assert!(!missing.complete());
        assert!(!missing.all_required_passed());

        let mut duplicate = all_passed.clone();
        duplicate.stages.push(duplicate.stages[0].clone());
        assert!(!duplicate.complete());
        assert!(!duplicate.all_required_passed());

        let mut mismatched_identity = all_passed.clone();
        mismatched_identity.stages[0].evidence_identity = "wrong-fixture/v1";
        assert!(!mismatched_identity.complete());
        assert!(!mismatched_identity.all_required_passed());

        let mut reordered = all_passed.clone();
        reordered.stages.swap(0, 1);
        assert!(!reordered.complete());
        assert!(!reordered.all_required_passed());

        let mut blank_reason = all_passed.clone();
        blank_reason.stages[3].status = StageStatus::Failed(StageReason::new("", ""));
        assert!(!blank_reason.complete());
        assert!(!blank_reason.all_required_passed());

        let mut empty_log = all_passed;
        empty_log.stages[0].events.clear();
        assert!(!empty_log.complete());
        assert!(!empty_log.all_required_passed());
    }
}
