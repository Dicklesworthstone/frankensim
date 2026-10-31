//! Goal-driven refinement of an accepted adaptive density design.
//! Enrichment inherits physical stiffness and occupancy without refiltering.
//! The accepted study is only borrowed; failed fine solves cannot replace it.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::{ElasticityError3, MaterialBodyLoad3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octant3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptiveJacobi3, AdaptivePreconditionError3, AdaptiveMultilevelOptions3, AdaptiveSolveSpace3,
};
use fs_solver::op::two_level::{AdditiveTwoLevel, TwoLevelBudget, TwoLevelError};
use fs_solver::op::multilevel::MultilevelError;
use fs_sparse::precond::Precond;
use fs_dwr::elasticity3::{GoalError3, GoalEstimate3, GoalFields3, GoalMarking3, GoalOptions3, dorfler3, estimate_reference_goal3};
use fs_dwr::elasticity3::reaction::{AffineGoalForm3, MaterialBodyLoads3, estimate_affine_assembled_goal3};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
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
/// One independent reference load, potentially combining a body force with
/// surface traction or pressure. Its pure laws must match the accepted solve.
/// Reintegrating the law on the enriched grid is NOT prolongating nodal forces.
#[derive(Clone, Copy)]
pub struct GoalReferenceLoad3<'a> {
    /// Physical reference-configuration load law, not a density-dependent load.
    pub load: ReferenceLoad3<'a>,
    /// Nonnegative objective weight, without implicit normalization.
    pub weight: f64,
}
// Keep the body-only public API and mixed-load API on one implementation,
// without allocating/copying a second load family or exposing a callback solver.
trait GoalLaw3 {
    fn load(&self) -> ReferenceLoad3<'_>;
    fn weight(&self) -> f64;
}
impl GoalLaw3 for GoalBodyLoad3<'_> {
    fn load(&self) -> ReferenceLoad3<'_> { ReferenceLoad3::body(self.density) }
    fn weight(&self) -> f64 { self.weight }
}
impl GoalLaw3 for GoalReferenceLoad3<'_> {
    fn load(&self) -> ReferenceLoad3<'_> { self.load }
    fn weight(&self) -> f64 { self.weight }
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

/// Prepared once per two-grid design and shared by every independent case.
pub(crate) struct MaterialGoal3<'a> {
    coarse: &'a MaterialBodyLoad3,
    fine: MaterialBodyLoad3,
    coarse_occupancy: &'a [f64],
    fine_occupancy: Vec<f64>,
    pub coarse_force: Vec<f64>,
    pub fine_force: Vec<f64>,
}
impl MaterialGoal3<'_> {
    pub(crate) fn loads(&self) -> MaterialBodyLoads3<'_> {
        MaterialBodyLoads3 { coarse: self.coarse, fine: &self.fine,
            coarse_occupancy: self.coarse_occupancy, fine_occupancy: &self.fine_occupancy }
    }
}

pub(crate) fn prepare_material_goal3<'a, O: Sdf3Elasticity, F: Sdf3Elasticity>(
    study: &'a CutDensityStudy3<O>, fine: &F, parents: &[usize],
    coarse_projected: Option<&'a [f64]>, control: &mut SolveControl<'_>,
) -> Result<Option<MaterialGoal3<'a>>, GoalRefinementError3> {
    let Some(coarse) = study.material_body_load() else { return Ok(None); };
    let projected = coarse_projected.ok_or(GoalRefinementError3::Invalid(
        "material body-load refinement requires explicit coarse physical occupancy"))?;
    if projected.len() != study.cells()
        || projected.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
        return Err(GoalRefinementError3::Invalid("invalid coarse material body-load occupancy"));
    }
    let fine_load = study.material_body_load_for(fine, control)?.ok_or(
        GoalRefinementError3::Invalid("missing prepared fine material body load"))?;
    let mut fine_occupancy = Vec::with_capacity(parents.len());
    for (cell, &parent) in parents.iter().enumerate() {
        if cell % 64 == 0 { control.checkpoint("sdf3-goal-material-occupancy")?; }
        fine_occupancy.push(*projected.get(parent).ok_or(
            GoalRefinementError3::Invalid("invalid material body-load parent"))?);
    }
    let coarse_force = study.material_body_force(projected, control)?.ok_or(
        GoalRefinementError3::Invalid("missing prepared coarse material body load"))?;
    let fine_force = crate::sdf3::body::controlled(control,
        |poll| fine_load.load(&fine_occupancy, poll))?;
    Ok(Some(MaterialGoal3 { coarse, fine: fine_load, coarse_occupancy: projected,
        fine_occupancy, coarse_force, fine_force }))
}

pub(crate) fn controlled_goal3<T>(control: &mut SolveControl<'_>, stage: &'static str,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, GoalError3>)
    -> Result<T, GoalRefinementError3> {
    let mut stopped = None;
    let result = action(&mut || match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => { stopped = Some(error); ControlFlow::Break(()) },
    });
    if let Some(error) = stopped { return Err(error.into()); }
    result.map_err(Into::into)
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
/// Internal retention for full load-basis assessments. Old callers allocate no
/// retained family. Only a complete successful estimator call returns these.
#[derive(Debug, Default)]
pub(crate) struct EnrichedLoadFields3 {
    pub forces: Vec<Vec<f64>>,
    pub displacements: Vec<Vec<f64>>,
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Body-only compliance enrichment using the shared reference-load path.
    /// Each coarse field is revalidated; no partially solved family is returned.
    pub fn estimate_compliance_enrichment(&self, enriched: AdaptiveElasticity3,
        loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, &[], loads, coarse_displacements, None, options, None, control)
    }
    /// Body-only enrichment with extra correction geometries BELOW the accepted
    /// grid. Extra levels are admitted only for the multilevel solve policy.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_compliance_enrichment_with_coarse_levels(&self, enriched: AdaptiveElasticity3,
        coarser: &[&AdaptiveElasticity3], loads: &[GoalBodyLoad3<'_>], coarse_displacements: &[Vec<f64>],
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, coarser, loads, coarse_displacements, None, options, None, control)
    }

    /// Refine pressure-, traction-, body-, or mixed-loaded compliance studies.
    /// Both operators must retain surface rules when any law uses the interface.
    /// Each law is reintegrated on its actual grid; coarse nodal loads are never
    /// prolonged and loads are never rescaled with density. Fine geometry and
    /// pure load laws must describe the same reference design as the coarse one.
    ///
    /// Validate the COMPLETE accepted load family before numeric preparation;
    /// share one fixed preconditioner and one cumulative solve budget. Stale
    /// fields, unavailable surface rules, cancellation or an exhausted solve
    /// return no partial estimates and cannot mutate the accepted study.
    /// Surface quadrature/transfer differences remain in the two-grid identity;
    /// the reported correction is not a continuum-error or shape certificate.
    /// A configured material acceleration requires the explicit-occupancy
    /// entry point so physical density is never inferred from stiffness.
    pub fn estimate_reference_compliance_enrichment(&self, enriched: AdaptiveElasticity3,
        loads: &[GoalReferenceLoad3<'_>], coarse_displacements: &[Vec<f64>], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, &[], loads, coarse_displacements, None, options, None, control)
    }

    /// Reference-load compliance with the accepted physical occupancy supplied
    /// explicitly for a configured material acceleration. Fine occupancy is
    /// parent-inherited directly; it is never refiltered or inferred from SIMP.
    /// The retained acceleration recipe prepares the fine operator's own map.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_reference_compliance_enrichment_with_occupancy(&self, enriched: AdaptiveElasticity3,
        loads: &[GoalReferenceLoad3<'_>], coarse_displacements: &[Vec<f64>], coarse_projected: &[f64],
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, &[], loads, coarse_displacements,
            Some(coarse_projected), options, None, control)
    }

    /// The same occupancy-aware experiment with explicit coarse correction
    /// geometries. Additional levels remain a multilevel-solver policy only.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_reference_compliance_enrichment_with_occupancy_and_coarse_levels(&self,
        enriched: AdaptiveElasticity3, coarser: &[&AdaptiveElasticity3], loads: &[GoalReferenceLoad3<'_>],
        coarse_displacements: &[Vec<f64>], coarse_projected: &[f64], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, coarser, loads, coarse_displacements,
            Some(coarse_projected), options, None, control)
    }
    /// Mixed reference-load enrichment with additional ordered coarse geometries.
    /// The accepted grid is always the first correction space. Additional levels
    /// are geometry only: stiffness is inherited/Galerkin-coarsened from the fine
    /// physical design. No current optimizer model or accepted field is replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_reference_compliance_enrichment_with_coarse_levels(&self, enriched: AdaptiveElasticity3,
        coarser: &[&AdaptiveElasticity3], loads: &[GoalReferenceLoad3<'_>], coarse_displacements: &[Vec<f64>],
        options: GoalRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ComplianceRefinement3, GoalRefinementError3> {
        self.estimate_load_family(enriched, coarser, loads, coarse_displacements, None, options, None, control)
    }

    pub(crate) fn estimate_reference_compliance_with_fields(&self, enriched: AdaptiveElasticity3,
        loads: &[GoalReferenceLoad3<'_>], coarse: &[Vec<f64>], options: GoalRefinementOptions3,
        control: &mut SolveControl<'_>) -> Result<(ComplianceRefinement3, EnrichedLoadFields3), GoalRefinementError3> {
        let mut fields = EnrichedLoadFields3::default();
        let estimate = self.estimate_load_family(enriched, &[], loads, coarse, None, options, Some(&mut fields), control)?;
        Ok((estimate, fields))
    }

    #[allow(clippy::too_many_arguments)]
    fn estimate_load_family<L: GoalLaw3>(&self, mut enriched: AdaptiveElasticity3,
        coarser: &[&AdaptiveElasticity3], loads: &[L], coarse_displacements: &[Vec<f64>],
        coarse_projected: Option<&[f64]>, options: GoalRefinementOptions3,
        retained: Option<&mut EnrichedLoadFields3>, control: &mut SolveControl<'_>)
        -> Result<ComplianceRefinement3, GoalRefinementError3> {
        control.checkpoint("sdf3-goal-start")?;
        if self.has_material_body_load() && coarse_projected.is_none() {
            return Err(GoalRefinementError3::Invalid(
                "material body-load refinement requires explicit coarse physical occupancy"));
        }
        self.require_material_body_recipe("sdf3-goal-material-recipe")?;
        if loads.is_empty() || loads.len() > options.max_load_cases || loads.len() != coarse_displacements.len()
            || !loads.iter().any(|l| l.weight() > 0.0) || loads.iter().any(|l| !l.weight().is_finite() || l.weight() < 0.0)
            || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance].iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0)
            || coarser.len() > 18 || (!coarser.is_empty() && !matches!(options.preconditioner, GoalPreconditioner3::Multilevel { .. })) {
            return Err(GoalRefinementError3::Invalid("invalid load family or numerical/resource policy"));
        }
        let coarse_operator = self.operator().adaptive();
        let initial_transfer = AdaptiveTransfer3::new(coarse_operator, &enriched, options.max_transfer_terms, || poll(control))?;
        let scales = initial_transfer.inherited_scales();
        let material = prepare_material_goal3(self, &enriched, initial_transfer.parents(), coarse_projected, control)?;
        drop(initial_transfer);
        enriched.set_scales(&scales)?;
        for (load, coarse) in loads.iter().zip(coarse_displacements) {
            // Reject the entire stale family BEFORE spending setup or a fine solve.
            let mut rhs_coarse = coarse_operator.reference_load(load.load(), || poll(control))?;
            crate::sdf3::body::add_force(&mut rhs_coarse, material.as_ref().map(|m| m.coarse_force.as_slice()), control)?;
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
            return complete_family(&transfer, loads, coarse_displacements, options, &prepared, material.as_ref(), retained, control);
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
        complete_family(&transfer, loads, coarse_displacements, options, &prepared, material.as_ref(), retained, control)
    }
}

// All policies and load kinds share the actual independent-load/DWR path.
#[allow(clippy::too_many_arguments)]
fn complete_family<L: GoalLaw3>(transfer: &AdaptiveTransfer3<'_>, loads: &[L], coarse_displacements: &[Vec<f64>],
    options: GoalRefinementOptions3, prepared: &impl Precond, material: Option<&MaterialGoal3<'_>>,
    mut retained: Option<&mut EnrichedLoadFields3>,
    control: &mut SolveControl<'_>)
    -> Result<ComplianceRefinement3, GoalRefinementError3> {
    let enriched = transfer.fine();
    let mut report = ComplianceRefinement3 { cases: Vec::with_capacity(loads.len()), weights: loads.iter().map(|l| l.weight()).collect(),
        coarse_value: 0.0, fine_value: 0.0, correction: 0.0, marking_mass: BTreeMap::new(), work: control.work(), preconditioner: options.preconditioner };
    for (load, coarse) in loads.iter().zip(coarse_displacements) {
        let law = load.load();
        let mut rhs = enriched.reference_load(law, || poll(control))?;
        crate::sdf3::body::add_force(&mut rhs, material.map(|m| m.fine_force.as_slice()), control)?;
        let fine = control.solve_preconditioned(enriched, prepared, &rhs, 1e-12, 50_000, "sdf3-goal-elasticity")?;
        // The DWR owner rechecks true residuals before using this CG output.
        let estimate = if let Some(body) = material {
            let mut coarse_force = transfer.coarse().reference_load(law, || poll(control))?;
            crate::sdf3::body::add_force(&mut coarse_force, Some(&body.coarse_force), control)?;
            let coarse_goal = AffineGoalForm3 { gradient: coarse_force, offset: 0.0 };
            let fine_goal = AffineGoalForm3 { gradient: rhs.clone(), offset: 0.0 };
            controlled_goal3(control, "sdf3-goal-material-estimate", |poll|
                estimate_affine_assembled_goal3(transfer, EquilibriumLoad3 { external: law, prescribed: None },
                    Some(body.loads()), &coarse_goal, &fine_goal,
                    GoalFields3::compliance(coarse, &fine), options.numerical, poll))?.linearized
        } else {
            estimate_reference_goal3(transfer, law, law,
                GoalFields3::compliance(coarse, &fine), options.numerical, || poll(control))?
        };
        report.coarse_value += load.weight() * estimate.coarse_value;
        report.fine_value += load.weight() * estimate.fine_value;
        report.correction += load.weight() * estimate.correction();
        for (&cell, value) in &estimate.cells {
            *report.marking_mass.entry(cell).or_insert(0.0) += load.weight() * value.marking_mass;
        }
        report.cases.push(estimate);
        if let Some(fields) = retained.as_deref_mut() {
            fields.forces.push(rhs);
            fields.displacements.push(fine);
        }
    }
    if ![report.coarse_value, report.fine_value, report.correction].iter().all(|v| v.is_finite())
        || !report.marking_mass.values().all(|v| v.is_finite()) { return Err(GoalRefinementError3::Invalid("weighted goal overflow")); }
    control.checkpoint("sdf3-goal-publish")?;
    report.work = control.work(); Ok(report)
}
