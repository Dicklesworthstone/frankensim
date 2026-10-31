//! Per-limit two-grid checks for minimum-material design. A summed loss or the
//! coarse active aggregate cannot stand in for an individual physical bound.
//! Every original limit gets its own actual adjoint pair; reaction offsets,
//! load/goal transfer, and algebraic defects use the existing affine DWR owner.
//! A selected mass objective adds direct native inventory corrections on that
//! same parent map, without adding a displacement or inventory adjoint solve.
use super::*;
use std::collections::BTreeMap;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::Octant3;
use fs_dwr::elasticity3::{GoalError3, GoalFields3, GoalMarking3, GoalOptions3, dorfler3};
use fs_dwr::elasticity3::reaction::{AffineGoalEstimate3, AffineMotionGoal3,
    MaterialBodyLoads3, assemble_affine_motion_goal3, estimate_affine_assembled_goal3};
use fs_dwr::elasticity3::response::EquilibriumLoad3;
use crate::sdf3::body::{add_force, controlled};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3};
use super::design::MinimumMaterialObjective3;

pub mod adaptive;
mod inventory;
mod refit;
pub use refit::{MinimumMaterialRefit3, ReferenceMinimumMaterialFit3};

/// One pure observation law, independently integrated on each grid. No nodal
/// observation is interpolated; a reaction retains its material dependence.
#[derive(Clone, Copy)]
pub enum ReferenceBoundedResponse3<'a> {
    Displacement(ReferenceLoad3<'a>),
    Reaction(&'a PrescribedMotion3<'a>),
}
impl<'a> ReferenceBoundedResponse3<'a> {
    fn goal(self) -> AffineMotionGoal3<'a> {
        match self {
            Self::Displacement(load) => AffineMotionGoal3 { displacement: load, reaction_mode: None },
            Self::Reaction(mode) => AffineMotionGoal3 { displacement: ReferenceLoad3::default(), reaction_mode: Some(mode) },
        }
    }
}
#[derive(Clone, Copy)]
pub struct ReferenceResponseLimit3<'a> {
    pub response: ReferenceBoundedResponse3<'a>,
    pub sense: ResponseLimitSense3,
    pub bound: f64,
    pub scale: f64,
}
#[derive(Clone, Copy)]
pub struct ReferenceResponseLimitCase3<'a> {
    pub load: ReferenceLoad3<'a>,
    pub prescribed: Option<&'a PrescribedMotion3<'a>>,
    pub limits: &'a [ReferenceResponseLimit3<'a>],
}
#[derive(Debug, Clone, Copy)]
pub struct ResponseLimitRefinementOptions3 {
    /// Same smoothing and family caps as the retained limit evaluation.
    pub limits: ResponseLimitOptions3,
    /// Actual field residual and affine two-grid identity gates.
    pub numerical: GoalOptions3,
    pub max_transfer_terms: usize,
    /// Same inventory as the retained design. Physical mass requires its
    /// native material densities; it is never inferred from geometric volume.
    pub objective: MinimumMaterialObjective3,
    /// Numerical normalization for the volume objective, as in the optimizer.
    pub volume_scale: f64,
}
impl Default for ResponseLimitRefinementOptions3 {
    fn default() -> Self {
        Self { limits: ResponseLimitOptions3::default(), numerical: GoalOptions3::default(), max_transfer_terms: 2_000_000,
            objective: MinimumMaterialObjective3::VolumeFraction, volume_scale: 1.0 }
    }
}

/// One bound assessed on both grids, in the original case/limit order.
#[derive(Debug, Clone)]
pub struct ResponseLimitChange3 {
    pub sense: ResponseLimitSense3,
    pub bound: f64,
    pub scale: f64,
    /// Numerical signed normalized inequalities (nonpositive means satisfied).
    pub coarse_residual: f64,
    pub fine_residual: f64,
    /// The physical affine goal; reaction offset transfer is kept separately
    /// inside this report, never fabricated into a localized residual indicator.
    pub goal: AffineGoalEstimate3,
}
impl ResponseLimitChange3 {
    /// Reconstructed fine-minus-coarse SIGNED residual, not a continuum bound.
    pub fn correction(&self) -> f64 { self.sense.sign()*self.goal.correction()/self.scale }
}
#[derive(Debug, Clone)]
pub struct ResponseLimitCaseRefinement3 {
    pub limits: Vec<ResponseLimitChange3>,
    pub fine_displacement: Vec<f64>,
}
/// Complete family for inherited PHYSICAL stiffness. Fine feasibility is a
/// numerical diagnostic of this enriched discretization, not a claim about a
/// refiltered fine design, all finer grids, shape perturbations or the continuum.
#[derive(Debug, Clone)]
pub struct ResponseLimitRefinement3 {
    pub cases: Vec<ResponseLimitCaseRefinement3>,
    pub coarse_constraint: f64,
    pub fine_constraint: f64,
    pub coarse_worst_residual: f64,
    pub fine_worst_residual: f64,
    /// Parent-inherited projected density integrated against each grid's measures.
    pub coarse_volume_fraction: f64,
    pub fine_volume_fraction: f64,
    pub objective: MinimumMaterialObjective3,
    pub volume_scale: f64,
    /// Actual inventory at the inherited physical design. None in volume mode.
    pub coarse_mass_kg: Option<f64>,
    pub fine_mass_kg: Option<f64>,
    /// Values normalized by the explicitly selected volume or kg scale.
    pub coarse_objective: f64,
    pub fine_objective: f64,
    /// Signed fine-minus-coarse kg correction for each ORIGINAL-grid parent.
    /// This is a direct inventory quadrature change, not an elasticity DWR
    /// term or a continuum bound. Empty for the legacy volume objective.
    pub mass_cell_corrections_kg: BTreeMap<Octant3, f64>,
    /// Parent-inherited projected density in the enriched operator's exact cell order.
    pub fine_projected_rho: Vec<f64>,
    /// Sum of normalized absolute response indicators over ALL limits, including those
    /// inactive in the coarse smooth maximum. Mass mode also adds each native
    /// inventory correction's absolute value divided by its explicit kg scale.
    /// Bound offsets do not localize.
    pub marking_mass: BTreeMap<Octant3, f64>,
    pub work: SolveWork,
}
impl ResponseLimitRefinement3 {
    pub fn fine_feasible(&self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.fine_constraint.is_finite()
            && self.fine_worst_residual.is_finite() && self.fine_constraint <= tolerance
            && self.fine_worst_residual <= tolerance
            && !self.cases.is_empty() && self.cases.iter().all(|c| !c.limits.is_empty()
                && c.limits.iter().all(|l| l.fine_residual.is_finite() && l.fine_residual <= tolerance))
    }
    /// Refuse cancellation between upper/lower bounds or independent cases.
    pub fn max_normalized_change(&self) -> f64 {
        self.cases.iter().flat_map(|c| &c.limits)
            .map(|l| (l.fine_residual-l.coarse_residual).abs()).fold(0.0, f64::max)
    }
    /// Marks ORIGINAL-grid cells, not the globally enriched probe. Empty marks
    /// with a failed bound are possible (e.g. offset transfer); never a pass.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut()->ControlFlow<()>)
        -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.marking_mass, theta, max_marks, checkpoint)
    }
}
fn invalid(s: &'static str) -> GoalRefinementError3 { GoalRefinementError3::Invalid(s) }
fn number(v: f64) -> Result<f64, GoalRefinementError3> {
    if v.is_finite() { Ok(v) } else { Err(invalid("nonfinite response-limit refinement arithmetic")) }
}
fn matches(a: f64, b: f64, scale: f64, tolerance: f64) -> bool {
    a.is_finite() && b.is_finite() && ((a-b)/scale).abs() <= tolerance
}
fn envelope(values: &[f64], sharpness: f64) -> Result<(f64, f64), GoalRefinementError3> {
    aggregate(values, sharpness).map(|(g, m, _)| (g, m))
        .map_err(|_| invalid("invalid enriched response-limit aggregate"))
}
pub(super) fn admit_reference_limits<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>,
    cases: &[ReferenceResponseLimitCase3<'_>], options: ResponseLimitOptions3) -> Result<(), ResponseError3> {
    let invalid = ResponseError3::Invalid;
    if cases.is_empty() || cases.len() > options.max_cases || !options.sharpness.is_finite()
        || options.sharpness <= 0.0 || !(1.0/options.sharpness).is_finite() {
        return Err(invalid("invalid reference limit family or sharpness"));
    }
    let mut count = 0usize;
    for case in cases {
        count = count.checked_add(case.limits.len()).filter(|n| *n <= options.max_limits)
            .ok_or(invalid("reference response-limit allowance exhausted"))?;
        if case.limits.is_empty() { return Err(invalid("empty reference response-limit case")); }
        if (case.prescribed.is_some() || case.limits.iter().any(|l| matches!(l.response, ReferenceBoundedResponse3::Reaction(_))))
            && study.operator.adaptive().embedded_dirichlet_penalty().is_none() {
            return Err(invalid("reference limit motion/reaction requires embedded supports"));
        }
        for limit in case.limits {
            if !limit.bound.is_finite() || !limit.scale.is_finite() || limit.scale <= 0.0 || !(1.0/limit.scale).is_finite() {
                return Err(invalid("invalid reference response bound or scale"));
            }
        }
    }
    Ok(())
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Check every physical bound on an enriched grid at the retained material.
    /// Revalidate the COMPLETE source family before preparing either backend.
    /// Each backend prepares once; one new primal per case and an adjoint pair
    /// per limit are charged to the original control. The design's aggregate
    /// adjoint cannot replace these individual observation adjoints.
    ///
    /// On success AND returned failure incoming source scales are restored.
    /// The probe is consumed; no partial result family escapes. Pure laws and
    /// unchanged geometry/supports/material are caller obligations. This does
    /// not assume nonlinear aggregate differences are a linear DWR sum.
    pub fn estimate_response_limit_enrichment<F: AdaptiveSdf3Elasticity>(&mut self, enriched: F,
        accepted: &ResponseLimitEvaluation3, cases: &[ReferenceResponseLimitCase3<'_>],
        options: ResponseLimitRefinementOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseLimitRefinement3, GoalRefinementError3> {
        self.estimate_response_limit_enrichment_retaining(enriched, accepted, cases, options, control)
            .map(|(report, _)| report)
    }

    /// Internal joint-design seam: retain the already-solved enriched operator
    /// so another diagnostic can reuse its fresh primal fields without another
    /// equilibrium family. The public response-only API still consumes it.
    pub(in crate::sdf3) fn estimate_response_limit_enrichment_retaining<F: AdaptiveSdf3Elasticity>(
        &mut self, mut enriched: F, accepted: &ResponseLimitEvaluation3,
        cases: &[ReferenceResponseLimitCase3<'_>], options: ResponseLimitRefinementOptions3,
        control: &mut SolveControl<'_>)
        -> Result<(ResponseLimitRefinement3, F), GoalRefinementError3> {
        control.checkpoint("response-limit-goal-start")?;
        self.require_material_body_recipe("response-limit-goal-requires-material-acceleration")?;
        admit_reference_limits(self, cases, options.limits).map_err(|_| invalid("invalid reference response-limit family"))?;
        inventory::admit(&self.operator, &enriched, accepted, options)?;
        if options.max_transfer_terms == 0 || cases.len() != accepted.displacements.len()
            || cases.len() != accepted.measurements.len() || cases.len() != accepted.residuals.len()
            || accepted.rho.len() != self.cells() || accepted.rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || ![options.numerical.residual_tolerance, options.numerical.identity_tolerance]
                .iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0) {
            return Err(invalid("invalid retained limit design or refinement controls"));
        }
        for (i, case) in cases.iter().enumerate() {
            if case.limits.len() != accepted.measurements[i].len() || case.limits.len() != accepted.residuals[i].len() {
                return Err(invalid("retained response-limit row count mismatch"));
            }
        }
        let design = self.design(&accepted.rho, control)?;
        let tolerance = options.numerical.identity_tolerance;
        if design.scales != accepted.scales || design.projected != accepted.projected_rho
            || !matches(design.volume, accepted.volume_fraction, 1.0, tolerance) {
            return Err(invalid("retained response-limit design belongs to another pipeline"));
        }
        let (coarse_mass_kg, coarse_objective) = inventory::measure(&self.operator,
            &design.projected, design.volume, options, control)?;
        if coarse_mass_kg != accepted.mass_kg {
            return Err(invalid("retained response-limit mass belongs to another material field"));
        }
        let coarse_body_force = self.material_body_force(&design.projected, control)?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let result = (|| {
            let coarse = self.operator.adaptive();
            let mut all_coarse = Vec::new();
            for (i, case) in cases.iter().enumerate() {
                let u = &accepted.displacements[i];
                let mut rhs = coarse.reference_load_with_motion(case.load, case.prescribed, || poll(control))?;
                add_force(&mut rhs, coarse_body_force.as_deref(), control)?;
                let residual = coarse.field_residual(u, &rhs, || poll(control))?;
                if !residual.is_finite() || residual > options.numerical.residual_tolerance {
                    return Err(GoalError3::FieldResidual { field: "coarse-limit-primal", value: residual }.into());
                }
                for (j, limit) in case.limits.iter().enumerate() {
                    let form = assemble_affine_motion_goal3(coarse, limit.response.goal(), case.prescribed, || poll(control))?;
                    let value = number(dot(&form.gradient, u)+form.offset)?;
                    let r = number(limit.sense.sign()*number(value-limit.bound)?/limit.scale)?;
                    if !matches(value, accepted.measurements[i][j], limit.scale, tolerance)
                        || !matches(r, accepted.residuals[i][j], 1.0, tolerance) {
                        return Err(invalid("retained limit law, bound, sense or scale mismatch"));
                    }
                    all_coarse.push(r);
                }
            }
            let (cg, cm) = envelope(&all_coarse, options.limits.sharpness)?;
            if !matches(cg, accepted.constraint, 1.0, tolerance) || !matches(cm, accepted.worst_residual, 1.0, tolerance) {
                return Err(invalid("retained response-limit aggregate mismatch"));
            }
            let inherited = AdaptiveTransfer3::new(coarse, enriched.adaptive(), options.max_transfer_terms,
                || poll(control))?.inherited_scales();
            enriched.set_scales(&inherited)?;
            let transfer = AdaptiveTransfer3::new(coarse, enriched.adaptive(), options.max_transfer_terms, || poll(control))?;
            if !enriched.adaptive().leaves().iter().zip(transfer.parents())
                .any(|(leaf, &parent)| leaf.level() > coarse.leaves()[parent].level()) {
                return Err(invalid("limit estimate requires an enriched active space"));
            }
            let volumes = enriched.volumes(); let total = number(volumes.iter().sum())?;
            if total <= 0.0 { return Err(invalid("nonpositive enriched cut volume")); }
            let fine_projected_rho: Vec<f64> = transfer.parents().iter()
                .map(|&p| design.projected[p]).collect();
            let fine_body = self.material_body_load_for(&enriched, control)?;
            let fine_body_force = fine_body.as_ref().map(|load|
                controlled(control, |poll| load.load(&fine_projected_rho, poll))).transpose()?;
            let fine_volume = number(volumes.iter().zip(&fine_projected_rho)
                .map(|(v, r)| (v/total)*r).sum())?;
            let (fine_mass_kg, fine_objective) = inventory::measure(&enriched,
                &fine_projected_rho, fine_volume, options, control)?;
            let mass_cell_corrections_kg = inventory::cell_corrections(&self.operator, &enriched,
                &design.projected, transfer.parents(), options, control)?;
            let mut marking_mass = BTreeMap::new();
            if let MinimumMaterialObjective3::MassKg { scale_kg } = options.objective {
                for (&leaf, &change) in &mass_cell_corrections_kg {
                    marking_mass.insert(leaf, inventory::normalized(change.abs(), scale_kg)?);
                }
            }
            let pc = self.operator.prepare_elasticity(control)?;
            let pf = enriched.prepare_elasticity(control)?;
            let mut report = ResponseLimitRefinement3 { cases: Vec::with_capacity(cases.len()),
                coarse_constraint: cg, coarse_worst_residual: cm, fine_constraint: 0.0, fine_worst_residual: 0.0,
                coarse_volume_fraction: design.volume, fine_volume_fraction: fine_volume,
                objective: options.objective, volume_scale: options.volume_scale,
                coarse_mass_kg, fine_mass_kg, coarse_objective, fine_objective, mass_cell_corrections_kg,
                fine_projected_rho, marking_mass, work: control.work() };
            let mut all_fine = Vec::new();
            for (i, case) in cases.iter().enumerate() {
                let mut rhs = enriched.adaptive().reference_load_with_motion(case.load, case.prescribed, || poll(control))?;
                add_force(&mut rhs, fine_body_force.as_deref(), control)?;
                let uf = checked_solve_preconditioned(&enriched, &pf, &rhs, 1e-11, "response-limit-goal-primal", control)?;
                let uc = &accepted.displacements[i];
                let mut changes = Vec::with_capacity(case.limits.len());
                for limit in case.limits {
                    control.checkpoint("response-limit-goal-observation")?;
                    let goal = limit.response.goal();
                    let qc = assemble_affine_motion_goal3(coarse, goal, case.prescribed, || poll(control))?;
                    let qf = assemble_affine_motion_goal3(enriched.adaptive(), goal, case.prescribed, || poll(control))?;
                    let zc = checked_solve_preconditioned(&self.operator, &pc, &qc.gradient, 1e-11,
                        "response-limit-goal-coarse-adjoint", control)?;
                    let zf = checked_solve_preconditioned(&enriched, &pf, &qf.gradient, 1e-11,
                        "response-limit-goal-fine-adjoint", control)?;
                    let body = self.material_body_load().zip(fine_body.as_ref()).map(|(coarse, fine)|
                        MaterialBodyLoads3 { coarse, fine, coarse_occupancy: &design.projected,
                            fine_occupancy: &report.fine_projected_rho });
                    let estimate = controlled_goal3(control, "response-limit-goal-estimate", |poll|
                        estimate_affine_assembled_goal3(&transfer,
                            EquilibriumLoad3 { external: case.load, prescribed: case.prescribed }, body, &qc, &qf,
                            GoalFields3 { coarse_primal: uc, fine_primal: &uf, coarse_adjoint: &zc, fine_adjoint: &zf },
                            options.numerical, poll))?;
                    let cr = number(limit.sense.sign()*number(estimate.coarse_value-limit.bound)?/limit.scale)?;
                    let fr = number(limit.sense.sign()*number(estimate.fine_value-limit.bound)?/limit.scale)?;
                    number(fr-cr)?; number(estimate.correction()/limit.scale)?;
                    for (&leaf, cell) in &estimate.linearized.cells {
                        let mass = report.marking_mass.entry(leaf).or_insert(0.0);
                        *mass = number(*mass+cell.marking_mass/limit.scale)?;
                    }
                    all_fine.push(fr);
                    changes.push(ResponseLimitChange3 { sense: limit.sense, bound: limit.bound, scale: limit.scale,
                        coarse_residual: cr, fine_residual: fr, goal: estimate });
                }
                report.cases.push(ResponseLimitCaseRefinement3 { limits: changes, fine_displacement: uf });
            }
            (report.fine_constraint, report.fine_worst_residual) = envelope(&all_fine, options.limits.sharpness)?;
            control.checkpoint("response-limit-goal-publish")?;
            report.work = control.work(); Ok(report)
        })();
        self.operator.set_scales(&previous).expect("previous admitted scales remain valid");
        result.map(|report| (report, enriched))
    }
}
