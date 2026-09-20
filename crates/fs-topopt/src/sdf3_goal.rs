//! Goal-driven refinement of an accepted adaptive density design.
//! Enrichment inherits physical stiffness scales, not a refiltered raw design.
//! The accepted study is only borrowed; failed fine solves cannot replace it.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::{ElasticityError3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octant3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptiveJacobi3, AdaptivePreconditionError3, AdaptiveMultilevelOptions3, AdaptiveSolveSpace3,
};
use fs_solver::op::two_level::{AdditiveTwoLevel, TwoLevelBudget, TwoLevelError};
use fs_solver::op::multilevel::MultilevelError;
use fs_sparse::precond::Precond;
use fs_dwr::elasticity3::{GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3, dorfler3, estimate_goal3};
use crate::{EvaluationStop, SolveControl, SolveWork};
use crate::sdf3::{AdaptiveSdf3Elasticity, Sdf3Elasticity, CutDensityStudy3};

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
    /// Recursive sparse correction. The accepted grid is its first coarse
    /// space; supply still-coarser geometries through the explicit ladder API
    /// when that grid is larger than the bottom-factor allowance.
    Multilevel { options: AdaptiveMultilevelOptions3 },
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
            AdaptivePreconditionError3::Coarse(TwoLevelError::Cancelled)
                | AdaptivePreconditionError3::Hierarchy(MultilevelError::Cancelled) => Self::Evaluation(EvaluationStop::Cancelled),
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
    /// Selected solve policy. Work separates sparse Galerkin products, fine
    /// setup applications and outer iterations; they are different cost units.
    pub preconditioner: GoalPreconditioner3,
}
impl ComplianceRefinement3 {
    /// Choose actual coarse cells to refine under a caller-specified mark cap.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, theta, max_marks, checkpoint)
    }
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Estimate an accepted design using one independently solved enriched RHS
    /// per case. All coarse fields are revalidated before numerical preparation.
    /// The enriched model inherits physical parent scales, not a refiltered raw
    /// design. Failed setup or solves cannot mutate the accepted study and do not
    /// publish a partial estimate family. Differences are not continuum bounds.
    ///
    /// This entry point supplies no extra coarse spaces. Recursive preparation
    /// can still use the accepted grid as its bottom if that grid fits the cap;
    /// otherwise use `estimate_compliance_enrichment_with_coarse_levels`.
    pub fn estimate_compliance_enrichment(&self, enriched: AdaptiveElasticity3,
        loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_compliance_enrichment_with_coarse_levels(enriched, &[], loads, coarse_displacements, options, control)
    }

    /// Recursive goal solve using an explicit sequence of geometries BELOW the
    /// accepted grid, nearest first. The accepted grid itself is always the
    /// first correction space for the enriched operator. Every deeper operator
    /// is Galerkin-coarsened from that same fine stiffness, not rediscretized.
    /// Extra geometries are rejected for other policies rather than ignored.
    /// Preparation is shared across loads; setup limits apply to the whole
    /// hierarchy per call and consumed work remains visible after failure.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_compliance_enrichment_with_coarse_levels(&self, mut enriched: AdaptiveElasticity3,
        coarser: &[&AdaptiveElasticity3], loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>],
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ComplianceRefinement3, GoalRefinementError3> {
        control.checkpoint("sdf3-goal-start")?;
        if loads.is_empty() || loads.len() > options.max_load_cases || loads.len() != coarse_displacements.len()
            || !loads.iter().any(|l| l.weight > 0.0) || loads.iter().any(|l| !l.weight.is_finite() || l.weight < 0.0)
            || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance].iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0)
            || coarser.len() > 18 || (!coarser.is_empty() && !matches!(options.preconditioner, GoalPreconditioner3::Multilevel { .. })) {
            return Err(GoalRefinementError3::Invalid("invalid load family or numerical/resource policy"));
        }
        let coarse_operator = self.operator().adaptive();
        let scales = AdaptiveTransfer3::new(coarse_operator, &enriched, options.max_transfer_terms, || poll(control))?.inherited_scales();
        enriched.set_scales(&scales)?;
        for (load, coarse) in loads.iter().zip(coarse_displacements) {
            // Reject the entire stale family BEFORE spending setup or a fine solve.
            let rhs_coarse = coarse_operator.body_load(load.density, || poll(control))?;
            let residual = coarse_operator.field_residual(coarse, &rhs_coarse, || poll(control))?;
            if residual > options.numerical.residual_tolerance {
                return Err(GoalError3::FieldResidual { field: "coarse-primal", value: residual }.into());
            }
        }
        if let GoalPreconditioner3::Multilevel { options: hierarchy_options } = options.preconditioner {
            let mut levels = Vec::with_capacity(coarser.len()+1);
            levels.push(coarse_operator); levels.extend_from_slice(coarser);
            let space = AdaptiveSolveSpace3::multilevel(enriched, &levels, hierarchy_options, || poll(control))?;
            let transfer = AdaptiveTransfer3::new(coarse_operator, space.elasticity(), options.max_transfer_terms, || poll(control))?;
            let prepared = space.prepare_elasticity(control)?;
            return complete_family(&transfer, loads, coarse_displacements, options, &prepared, control);
        }
        let transfer = AdaptiveTransfer3::new(coarse_operator, &enriched, options.max_transfer_terms, || poll(control))?;
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
            GoalPreconditioner3::Multilevel { .. } => unreachable!("recursive branch returns above"),
        };
        complete_family(&transfer, loads, coarse_displacements, options, &prepared, control)
    }
}

// All policies share the exact independent-load solve and DWR acceptance path.
fn complete_family(transfer: &AdaptiveTransfer3<'_>, loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>],
    options: GoalRefinementOptions3, prepared: &impl Precond, control: &mut SolveControl<'_>)
    -> Result<ComplianceRefinement3, GoalRefinementError3> {
    let enriched = transfer.fine();
    let mut report = ComplianceRefinement3 { cases: Vec::with_capacity(loads.len()), weights: loads.iter().map(|l| l.weight).collect(),
        coarse_value: 0.0, fine_value: 0.0, correction: 0.0, marking_mass: BTreeMap::new(), work: control.work(), preconditioner: options.preconditioner };
    for (load, coarse) in loads.iter().zip(coarse_displacements) {
        let rhs = enriched.body_load(load.density, || poll(control))?;
        let fine = control.solve_preconditioned(enriched, prepared, &rhs, 1e-12, 50_000, "sdf3-goal-elasticity")?;
        // The DWR owner rechecks true residuals before using this CG output.
        let estimate = estimate_goal3(transfer, load.density, load.density,
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
