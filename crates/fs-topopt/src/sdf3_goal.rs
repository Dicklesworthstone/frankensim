//! Goal-driven refinement of an accepted adaptive density design.
//! Enrichment inherits physical stiffness scales, not a refiltered raw design.
//! The accepted study is only borrowed; failed fine solves cannot replace it.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::{ElasticityError3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octant3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{AdaptiveJacobi3, AdaptivePreconditionError3};
use fs_solver::op::two_level::{AdditiveTwoLevel, TwoLevelBudget, TwoLevelError};
use fs_sparse::precond::Precond;
use fs_dwr::elasticity3::{GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3, dorfler3, estimate_goal3};
use crate::{EvaluationStop, SolveControl, SolveWork};
use crate::sdf3::CutDensityStudy3;

/// One independent dead body-force density with its objective weight.
/// The callback must be pure and identical to the accepted coarse solve's law.
#[derive(Clone, Copy)]
pub struct GoalBodyLoad3<'a> {
    /// Force per volume in the geometry's declared units.
    pub density: &'a dyn Fn([f64; 3]) -> [f64; 3],
    /// Nonnegative weight; weights are not automatically normalized.
    pub weight: f64,
}
/// Explicit enriched-solve policy. The default preserves the original solver;
/// the runnable adaptive example opts into bounded two-level preparation.
#[derive(Debug, Clone, Copy)]
pub enum GoalPreconditioner3 {
    /// Original identity action, without setup applications.
    Identity,
    /// Exact density/constraint-aware Jacobi, with bounded accumulation.
    Jacobi { max_contributions: usize },
    /// Fixed SPD Galerkin correction, prepared once for the entire load family.
    TwoLevel { budget: TwoLevelBudget, max_diagonal_contributions: usize },
}
enum PreparedGoal3<'a> {
    Identity,
    Jacobi(AdaptiveJacobi3<'a>),
    TwoLevel(AdditiveTwoLevel<'a, AdaptiveElasticity3>),
}
impl Precond for PreparedGoal3<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self {
            Self::Identity => z.copy_from_slice(r),
            Self::Jacobi(p) => p.apply(r, z),
            Self::TwoLevel(p) => p.apply(r, z),
        }
    }
}

/// Limits beyond the existing shared linear-work control and geometry budget.
#[derive(Debug, Clone, Copy)]
pub struct GoalRefinementOptions3 {
    /// Total scalar coefficients in the sparse inter-grid transfer.
    pub max_transfer_terms: usize,
    /// Maximum independently solved cases in this request.
    pub max_load_cases: usize,
    /// Actual-field and two-grid numerical identity gates.
    pub numerical: GoalOptions3,
    /// Reused across all independent enriched loads; never silently downgraded.
    pub preconditioner: GoalPreconditioner3,
}
impl Default for GoalRefinementOptions3 {
    fn default() -> Self {
        Self { max_transfer_terms: 2_000_000, max_load_cases: 64, numerical: GoalOptions3::default(),
            preconditioner: GoalPreconditioner3::Identity }
    }
}
/// No partial family of estimates is usable on a returned error.
#[derive(Debug, Clone, PartialEq)]
pub enum GoalRefinementError3 {
    /// Linear-work exhaustion, cancellation, or numerical breakdown.
    Evaluation(EvaluationStop),
    /// Incompatible geometry/fields or a failed physical callback.
    Physics(ElasticityError3),
    /// Refused setup size, diagonal, or coarse factorization.
    Preconditioner(AdaptivePreconditionError3),
    /// Numerical estimator refused its assumptions or identity.
    Estimate(GoalError3),
    /// Invalid load family or explicit resource limits.
    Invalid(&'static str),
}
impl From<EvaluationStop> for GoalRefinementError3 { fn from(e: EvaluationStop) -> Self { Self::Evaluation(e) } }
impl From<ElasticityError3> for GoalRefinementError3 {
    fn from(e: ElasticityError3) -> Self {
        if matches!(e, ElasticityError3::Cancelled) { Self::Evaluation(EvaluationStop::Cancelled) } else { Self::Physics(e) }
    }
}
impl From<AdaptivePreconditionError3> for GoalRefinementError3 {
    fn from(e: AdaptivePreconditionError3) -> Self {
        match e {
            AdaptivePreconditionError3::Physics(e) => e.into(),
            AdaptivePreconditionError3::Coarse(TwoLevelError::Cancelled) => Self::Evaluation(EvaluationStop::Cancelled),
            other => Self::Preconditioner(other),
        }
    }
}
impl From<GoalError3> for GoalRefinementError3 {
    fn from(e: GoalError3) -> Self {
        match e { GoalError3::Physics(e) => e.into(), other => Self::Estimate(other) }
    }
}
impl std::fmt::Display for GoalRefinementError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "3-D goal refinement stopped: {self:?}") }
}
impl std::error::Error for GoalRefinementError3 {}
fn poll(control: &mut SolveControl<'_>) -> ControlFlow<()> {
    if control.checkpoint("sdf3-goal").is_ok() { ControlFlow::Continue(()) } else { ControlFlow::Break(()) }
}

/// Complete independent-load estimates at one accepted physical design.
#[derive(Debug, Clone)]
pub struct ComplianceRefinement3 {
    /// Case-local signed correction decompositions, in original load order.
    pub cases: Vec<GoalEstimate3>,
    /// Corresponding original weights, without implicit normalization.
    pub weights: Vec<f64>,
    /// Weighted compliance from revalidated coarse fields.
    pub coarse_value: f64,
    /// Weighted enriched compliance at the same inherited stiffness field.
    pub fine_value: f64,
    /// Weighted two-level difference including consistency and algebraic terms.
    pub correction: f64,
    /// Sum of weighted absolute local masses. Opposite independent loads and
    /// signed residuals cannot cancel each other's refinement signal.
    pub marking_mass: BTreeMap<Octant3, f64>,
    /// Cumulative solve work, including earlier optimization and this enrichment.
    pub work: SolveWork,
    /// Actual explicitly selected solve policy. Setup applications are included
    /// separately in `work.preconditioner_operator_applications`.
    pub preconditioner: GoalPreconditioner3,
}
impl ComplianceRefinement3 {
    /// Choose actual coarse cells to refine under a caller-specified mark cap.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, theta, max_marks, checkpoint)
    }
}
impl CutDensityStudy3<AdaptiveElasticity3> {
    /// Estimate the current accepted design using one enriched solve per case.
    ///
    /// `coarse_displacements` must be the complete accepted fields, in body-load
    /// order. Recomputed residuals reject stale states or different load laws.
    /// `enriched` must come from a refinement of the SAME implicit design domain,
    /// reference material and box supports. It is consumed so preparation and
    /// failure cannot mutate a caller's other accepted operator. Its existing
    /// stiffness scales are replaced by exact parent-inherited physical scales.
    /// No density filtering, optimization step or coarse re-solve occurs here.
    ///
    /// All enriched Krylov work uses `control`, including failed solves. No
    /// partial estimate family or mark set escapes on cancellation/exhaustion.
    /// All coarse fields are checked before preconditioner setup. Preparation
    /// happens once per load family and is invalidated with the consumed fine
    /// operator. Setup fine-applications remain visible in `control.work()` even
    /// on an interrupted setup or a later failed solve; they do not consume the
    /// outer-iteration allowance and have their own explicit setup budget.
    /// The returned differences are not certified continuum-error bounds.
    pub fn estimate_compliance_enrichment(&self, mut enriched: AdaptiveElasticity3,
        loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<ComplianceRefinement3, GoalRefinementError3> {
        control.checkpoint("sdf3-goal-start")?;
        if loads.is_empty() || loads.len() > options.max_load_cases || loads.len() != coarse_displacements.len()
            || !loads.iter().any(|l| l.weight > 0.0) || loads.iter().any(|l| !l.weight.is_finite() || l.weight < 0.0)
            || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance].iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0) {
            return Err(GoalRefinementError3::Invalid("invalid load family or numerical/resource policy"));
        }
        let scales = AdaptiveTransfer3::new(self.operator(), &enriched, options.max_transfer_terms, || poll(control))?.inherited_scales();
        enriched.set_scales(&scales)?;
        let transfer = AdaptiveTransfer3::new(self.operator(), &enriched, options.max_transfer_terms, || poll(control))?;
        for (load, coarse) in loads.iter().zip(coarse_displacements) {
            // Reject the entire stale family BEFORE spending setup or a fine solve.
            let rhs_coarse = self.operator().body_load(load.density, || poll(control))?;
            let residual = self.operator().field_residual(coarse, &rhs_coarse, || poll(control))?;
            if residual > options.numerical.residual_tolerance {
                return Err(GoalError3::FieldResidual { field: "coarse-primal", value: residual }.into());
            }
        }
        let prepared = match options.preconditioner {
            GoalPreconditioner3::Identity => PreparedGoal3::Identity,
            GoalPreconditioner3::Jacobi { max_contributions } => {
                PreparedGoal3::Jacobi(enriched.prepare_jacobi(max_contributions, || poll(control))?)
            }
            GoalPreconditioner3::TwoLevel { budget, max_diagonal_contributions } => {
                let mut recorded = 0usize;
                let mut setup_stop = None;
                let result = transfer.prepare_two_level(budget, max_diagonal_contributions, |work| {
                    let Some(additional) = work.operator_applications.checked_sub(recorded) else {
                        setup_stop = Some(EvaluationStop::Breakdown { stage: "preconditioner-accounting" });
                        return ControlFlow::Break(());
                    };
                    recorded = work.operator_applications;
                    match control.record_preconditioner_applications(additional) {
                        Ok(()) => ControlFlow::Continue(()),
                        Err(stop) => { setup_stop = Some(stop); ControlFlow::Break(()) }
                    }
                });
                if let Some(stop) = setup_stop { return Err(stop.into()); }
                PreparedGoal3::TwoLevel(result?)
            }
        };
        let mut report = ComplianceRefinement3 { cases: Vec::with_capacity(loads.len()), weights: loads.iter().map(|l| l.weight).collect(),
            coarse_value: 0.0, fine_value: 0.0, correction: 0.0, marking_mass: BTreeMap::new(), work: control.work(), preconditioner: options.preconditioner };
        for (load, coarse) in loads.iter().zip(coarse_displacements) {
            let rhs = enriched.body_load(load.density, || poll(control))?;
            let fine = control.solve_preconditioned(&enriched, &prepared, &rhs, 1e-12, 50_000, "sdf3-goal-elasticity")?;
            // The DWR owner rechecks true residuals before using this CG output.
            let estimate = estimate_goal3(&transfer, load.density, load.density,
                GoalFields3::compliance(coarse, &fine), options.numerical, || poll(control))?;
            report.coarse_value += load.weight * estimate.coarse_value;
            report.fine_value += load.weight * estimate.fine_value;
            report.correction += load.weight * estimate.correction();
            for (&cell, value) in &estimate.cells {
                *report.marking_mass.entry(cell).or_insert(0.0) += load.weight * value.marking_mass;
            }
            report.cases.push(estimate);
        }
        if ![report.coarse_value, report.fine_value, report.correction].iter().all(|v| v.is_finite())
            || !report.marking_mass.values().all(|v| v.is_finite()) { return Err(GoalRefinementError3::Invalid("weighted goal overflow")); }
        control.checkpoint("sdf3-goal-publish")?;
        report.work = control.work(); Ok(report)
    }
}
