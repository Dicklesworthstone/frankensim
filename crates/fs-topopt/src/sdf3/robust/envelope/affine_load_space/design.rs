//! Constrained design for a common load plus a continuous uncertain load set.
//! One existing projected-AL state owns both compliance and inventory modes.
use super::constraints::{AffineLoadSpaceConstraintPolicy3, AffineLoadSpaceConstraints3};
use super::response::AffineResponseLimits3;
use super::{AffineLoadSpace3, AffineLoadSpaceEvaluation3};
use crate::gradient_check::{GradientSample, audit_design_map};
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3::robust::RobustOptions3;
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{EvaluationStop, GradientCheckOptions, SimpParams, SolveControl, SolveWork};
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlOptions, ProjectedAlReport, ProjectedAlSample, ProjectedAlState,
    ProjectedAlStop, ProjectedAlWork,
};
use std::{cell::RefCell, ops::ControlFlow};

/// The same physical load set serves either constrained design formulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AffineLoadSpaceDesignGoal3 {
    /// Minimize normalized worst compliance with a DILATED volume cap.
    /// The policy reference is a numerical compliance scale, not a requirement.
    WorstCompliance { volume_limit: f64 },
    /// Minimize DILATED inventory subject to worst compliance <= policy.reference().
    /// `volume_scale` normalizes volume mode; MassKg uses its explicit kg scale.
    /// There is no additional volume constraint in this formulation.
    MinimumMaterial {
        objective: MinimumMaterialObjective3,
        volume_scale: f64,
    },
}
impl AffineLoadSpaceDesignGoal3 {
    fn needs_mass(self) -> bool {
        matches!(
            self,
            Self::MinimumMaterial {
                objective: MinimumMaterialObjective3::MassKg { .. },
                ..
            }
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AffineLoadSpaceDesignOptions3 {
    pub robustness: RobustOptions3,
    pub density_floor: f64,
    /// Actual full-chain objective AND complete-constraint check.
    pub gradient_check: GradientCheckOptions,
    pub optimizer: ProjectedAlOptions,
}
impl Default for AffineLoadSpaceDesignOptions3 {
    fn default() -> Self {
        Self {
            robustness: RobustOptions3::default(),
            density_floor: 1e-3,
            gradient_check: GradientCheckOptions::default(),
            optimizer: ProjectedAlOptions::default(),
        }
    }
}

/// Numerical directional evidence for the ACTUAL optimizer objective and
/// constraint, including any joint response limits and their aggregation.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceGradientCheck3 {
    pub params: SimpParams,
    pub goal: AffineLoadSpaceDesignGoal3,
    pub options: GradientCheckOptions,
    pub constraints: AffineLoadSpaceConstraintPolicy3,
    pub baseline_objective: f64,
    pub baseline_constraint: f64,
    pub objective_relative_error: f64,
    pub constraint_relative_error: f64,
    pub directions: usize,
    pub work: SolveWork,
}
impl AffineLoadSpaceGradientCheck3 {
    #[must_use]
    pub fn passed(&self) -> bool {
        self.directions > 0
            && [
                self.objective_relative_error,
                self.constraint_relative_error,
            ]
            .iter()
            .all(|e| e.is_finite() && *e >= 0.0 && *e <= self.options.relative_tolerance)
    }
}

#[derive(Debug, Clone)]
pub struct AffineLoadSpaceDesignPoint3 {
    pub rho: Vec<f64>,
    pub evaluation: AffineLoadSpaceEvaluation3,
}
#[derive(Debug, Clone, PartialEq)]
pub struct AffineLoadSpaceDesignIteration3 {
    pub iteration: usize,
    pub objective: f64,
    pub constraint: f64,
    pub maximum: f64,
    pub worst_compliance: f64,
    pub dilated_volume_fraction: f64,
    pub dilated_mass_kg: Option<f64>,
    /// The physical cap passes directly, independently of AL's KKT tolerance.
    pub feasible: bool,
}

fn normalize(value: f64, scale: f64) -> Result<f64, EvaluationStop> {
    let normalized = value / scale;
    if normalized.is_finite() && (value == 0.0 || normalized != 0.0) {
        Ok(normalized)
    } else {
        Err(failure("sdf3-affine-design-inventory-range"))
    }
}
fn inventory(
    e: &AffineLoadSpaceEvaluation3,
    goal: AffineLoadSpaceDesignGoal3,
) -> Result<(f64, Vec<f64>), EvaluationStop> {
    let p = &e.physical;
    let (value, gradient, scale) = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => (
            p.dilated_volume_fraction,
            &p.dilated_volume_gradient,
            volume_limit,
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction,
            volume_scale,
        } => (
            p.dilated_volume_fraction,
            &p.dilated_volume_gradient,
            volume_scale,
        ),
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg },
            ..
        } => (
            p.dilated_mass_kg
                .ok_or(failure("sdf3-affine-design-mass"))?,
            p.dilated_mass_gradient
                .as_ref()
                .ok_or(failure("sdf3-affine-design-mass"))?,
            scale_kg,
        ),
    };
    if !value.is_finite() || value < 0.0 {
        return Err(failure("sdf3-affine-design-inventory"));
    }
    Ok((
        normalize(value, scale)?,
        gradient
            .iter()
            .map(|&g| normalize(g, scale))
            .collect::<Result<_, _>>()?,
    ))
}
pub(super) fn physical_feasible(
    e: &AffineLoadSpaceEvaluation3,
    goal: AffineLoadSpaceDesignGoal3,
    policy: &AffineLoadSpace3,
) -> bool {
    let original = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => {
            e.physical.dilated_volume_fraction.is_finite()
                && e.physical.dilated_volume_fraction <= volume_limit
        }
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => {
            e.maximum.is_finite()
                && (0.0..=1.0).contains(&e.maximum)
                && e.worst_compliance.is_finite()
                && e.worst_compliance >= 0.0
                && e.worst_compliance <= policy.reference()
        }
    };
    original
        && e.constraints.matches(e)
        && e.responses.as_ref().is_none_or(|r| r.physically_feasible())
        && e.stress.as_ref().is_none_or(|s| s.physically_feasible())
        && constraint_value(e, goal).is_ok_and(|constraint| constraint <= 0.0)
}

fn aggregate_constraints(
    base: f64,
    e: &AffineLoadSpaceEvaluation3,
) -> Result<(f64, [f64; 3]), EvaluationStop> {
    if !e.constraints.matches(e) {
        return Err(failure("sdf3-affine-design-constraint-binding"));
    }
    super::constraints::aggregate(
        base,
        e.responses.as_ref().map(|r| r.constraint),
        e.stress.as_ref().map(|s| s.constraint),
        e.constraints.sharpness,
    )
}

fn constraint_value(
    e: &AffineLoadSpaceEvaluation3,
    goal: AffineLoadSpaceDesignGoal3,
) -> Result<f64, EvaluationStop> {
    let base = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => {
            normalize(e.physical.dilated_volume_fraction, volume_limit)? - 1.0
        }
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => e.maximum - 1.0,
    };
    aggregate_constraints(base, e).map(|(value, _)| value)
}

// Both constructors, gradient probes, line-search trials and fresh rechecks
// use the same physics dispatch. A bound family cannot disappear on a trial.
#[allow(clippy::too_many_arguments)]
fn evaluate<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    rho: &[f64],
    dead: &[f64],
    basis: &[&[f64]],
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
    constraints: AffineLoadSpaceConstraints3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
    let result = study.affine_load_space_constraints_with_mass(
        rho,
        dead,
        basis,
        options.robustness,
        policy,
        constraints,
        goal.needs_mass(),
        control,
    )?;
    if let Some(limit) = constraints.stress {
        for (i, &r) in result.physical.eroded_rho.iter().enumerate() {
            let prescribed_void = r == 0.0
                && study
                    .physical_regions()
                    .is_some_and(|regions| regions[i] == crate::sdf3::PhysicalRegion3::Void);
            if !prescribed_void
                && !crate::sdf3::design::monotone_stress_branch(
                    r,
                    study.params(),
                    limit.options.relaxation_power,
                )
            {
                return Err(failure("sdf3-affine-design-stress-turnover"));
            }
        }
    }
    Ok(result)
}

fn optimization_sample(
    e: &AffineLoadSpaceEvaluation3,
    goal: AffineLoadSpaceDesignGoal3,
    gradients: bool,
) -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let gradient = match &e.gradient {
        Some(gradient) => gradient.clone(),
        None if !gradients => Vec::new(),
        None => return Ok(None),
    };
    let (material, material_gradient) = inventory(e, goal)?;
    let mut sample = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => ProjectedAlSample {
            objective: e.maximum,
            gradient,
            constraint: material - 1.0,
            constraint_gradient: material_gradient,
        },
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => ProjectedAlSample {
            objective: material,
            gradient: material_gradient,
            constraint: e.maximum - 1.0,
            constraint_gradient: gradient,
        },
    };
    if e.constraints.has_limits() {
        let (constraint, weights) = aggregate_constraints(sample.constraint, e)?;
        if gradients {
            for gradient in &mut sample.constraint_gradient {
                *gradient *= weights[0];
            }
            for (weight, component) in [
                (weights[1], e.responses.as_ref().map(|r| &r.gradient)),
                (weights[2], e.stress.as_ref().map(|s| &s.gradient)),
            ] {
                let Some(component) = component else {
                    continue;
                };
                let Some(component) = component else {
                    return Ok(None);
                };
                if component.len() != sample.constraint_gradient.len() {
                    return Err(failure("sdf3-affine-design-constraint-gradient-shape"));
                }
                for (gradient, value) in sample.constraint_gradient.iter_mut().zip(component) {
                    *gradient = super::grid::finite(
                        *gradient + weight * value,
                        "sdf3-affine-design-constraint-gradient",
                    )?;
                }
            }
        } else {
            sample.constraint_gradient.clear();
        }
        sample.constraint = constraint;
    }
    Ok(Some(sample))
}
// Rank physical inventory directly: normalizing two distinct feasible masses
// can round them to the same quotient. AL/history keep their declared scales.
fn archive_objective(e: &AffineLoadSpaceEvaluation3, goal: AffineLoadSpaceDesignGoal3) -> f64 {
    match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => e.maximum,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction,
            ..
        } => e.physical.dilated_volume_fraction,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { .. },
            ..
        } => e
            .physical
            .dilated_mass_kg
            .expect("admitted native mass objective"),
    }
}
fn row(
    point: &AffineLoadSpaceDesignPoint3,
    iteration: usize,
    sample: &ProjectedAlSample,
    goal: AffineLoadSpaceDesignGoal3,
    policy: &AffineLoadSpace3,
) -> AffineLoadSpaceDesignIteration3 {
    let e = &point.evaluation;
    AffineLoadSpaceDesignIteration3 {
        iteration,
        objective: sample.objective,
        constraint: sample.constraint,
        maximum: e.maximum,
        worst_compliance: e.worst_compliance,
        dilated_volume_fraction: e.physical.dilated_volume_fraction,
        dilated_mass_kg: e.physical.dilated_mass_kg,
        feasible: physical_feasible(e, goal, policy),
    }
}

#[allow(clippy::too_many_arguments)]
fn sample<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    rho: &[f64],
    dead: &[f64],
    basis: &[&[f64]],
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
    constraints: AffineLoadSpaceConstraints3<'_>,
    last: &mut Option<AffineLoadSpaceDesignPoint3>,
    control: &mut SolveControl<'_>,
) -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let e = evaluate(
        study,
        rho,
        dead,
        basis,
        policy,
        goal,
        options,
        constraints,
        control,
    )?;
    // A valid worst-load witness can exist without a unique differentiable
    // envelope. Let the existing bounded line search reject that trial.
    let Some(sample) = optimization_sample(&e, goal, true)? else {
        return Ok(None);
    };
    *last = Some(AffineLoadSpaceDesignPoint3 {
        rho: rho.to_vec(),
        evaluation: e,
    });
    Ok(Some(sample))
}

#[allow(clippy::too_many_arguments)]
fn admit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    dead: &[f64],
    basis: &[&[f64]],
    rho: &[f64],
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
) -> Result<(), ProjectedAlError<EvaluationStop>> {
    admit_model(study, rho, basis.len(), policy, goal, options)?;
    if dead.len() != study.operator().n()
        || dead.iter().any(|f| !f.is_finite())
        || basis
            .iter()
            .any(|f| f.len() != dead.len() || f.iter().any(|v| !v.is_finite()))
    {
        return Err(ProjectedAlError::Invalid(
            "invalid affine load-space design inputs or gradient policy",
        ));
    }
    Ok(())
}

// Shared with reference-law fitting so raw bounds, numerical controls, native
// inventory and protected-solid constraints are admitted before callbacks.
pub(super) fn admit_model<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    rho: &[f64],
    basis_count: usize,
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
) -> Result<(), ProjectedAlError<EvaluationStop>> {
    let invalid = ProjectedAlError::Invalid;
    options
        .optimizer
        .validate::<EvaluationStop>(study.cells())?;
    if !options.robustness.valid(study.params())
        || basis_count != policy.dimension()
        || !options.density_floor.is_finite()
        || options.density_floor <= 0.0
        || options.density_floor >= 1.0
        || rho.len() != study.cells()
        || rho
            .iter()
            .any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0)
        || !options.gradient_check.step.is_finite()
        || options.gradient_check.step <= 0.0
        || options.gradient_check.step > 0.125
        || !options.gradient_check.relative_tolerance.is_finite()
        || options.gradient_check.relative_tolerance <= 0.0
        || options.gradient_check.relative_tolerance >= 1.0
    {
        return Err(invalid(
            "invalid affine load-space design inputs or gradient policy",
        ));
    }
    let scale = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit } => {
            if volume_limit > 1.0 || study.prescribed_solid_fraction() > volume_limit {
                return Err(invalid("invalid or unattainable dilated volume cap"));
            }
            volume_limit
        }
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction,
            volume_scale,
        } => volume_scale,
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg },
            ..
        } => {
            let masses = study.operator().reference_cell_masses_kg().ok_or(invalid(
                "affine mass design requires native material density",
            ))?;
            if masses.len() != study.cells() || masses.iter().any(|m| !m.is_finite() || *m <= 0.0) {
                return Err(invalid("invalid native material inventory"));
            }
            scale_kg
        }
    };
    if !scale.is_finite() || scale <= 0.0 || !(1.0 / scale).is_finite() {
        return Err(invalid(
            "affine design requires a positive invertible inventory scale",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    rho: &[f64],
    dead: &[f64],
    basis: &[&[f64]],
    policy: &AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
    constraints: AffineLoadSpaceConstraints3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<AffineLoadSpaceGradientCheck3, EvaluationStop> {
    let params = study.params();
    let numerical = audit_design_map(
        params,
        rho,
        options.gradient_check,
        control,
        |r, baseline, control| {
            let e = evaluate(
                study,
                r,
                dead,
                basis,
                policy,
                goal,
                options,
                constraints,
                control,
            )?;
            let sample = optimization_sample(&e, goal, baseline)?
                .ok_or(failure("sdf3-affine-design-nondifferentiable-baseline"))?;
            // With a single normalized cap, audit before subtracting its
            // constant one. Otherwise a loose cap can round away a perfectly
            // resolvable tiny compliance variation during finite differences.
            // This additive constant does not change the audited derivative.
            let constraint = if !constraints.policy().has_limits() {
                match goal {
                    AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => inventory(&e, goal)?.0,
                    AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => e.maximum,
                }
            } else {
                sample.constraint
            };
            Ok(GradientSample {
                compliance: sample.objective,
                compliance_gradient: sample.gradient,
                volume: constraint,
                volume_gradient: sample.constraint_gradient,
            })
        },
    )?;
    Ok(AffineLoadSpaceGradientCheck3 {
        params,
        goal,
        options: options.gradient_check,
        constraints: constraints.policy(),
        baseline_objective: numerical.baseline_compliance,
        baseline_constraint: numerical.baseline_volume_fraction
            - if !constraints.policy().has_limits() {
                1.0
            } else {
                0.0
            },
        objective_relative_error: numerical
            .probes
            .iter()
            .map(|p| p.compliance_relative_error)
            .fold(0.0, f64::max),
        constraint_relative_error: numerical
            .probes
            .iter()
            .map(|p| p.volume_relative_error)
            .fold(0.0, f64::max),
        directions: numerical.probes.len(),
        work: numerical.work,
    })
}
fn poll(control: &mut SolveControl<'_>, stopped: &mut Option<EvaluationStop>) -> ControlFlow<()> {
    match control.checkpoint("sdf3-affine-design-control") {
        Ok(()) => ControlFlow::Continue(()),
        Err(e) => {
            *stopped = Some(e);
            ControlFlow::Break(())
        }
    }
}

/// Fixed-grid design for f = f_dead + B*rho_eroded + F*a, ||a||<=1.
/// All uncertainty amplitudes/correlations, the common acceleration, material,
/// filter and projection policy stay fixed throughout a run.
///
/// The native global load oracle supplies the exact-discrete formulation's
/// numerical maximum. A gradient is used only on its separated regular branch.
/// A hard/near-hard baseline refuses optimization; a nonsmooth trial is rejected
/// by the existing bounded AL line search. No chosen-witness derivative is
/// presented as a smooth envelope gradient. All actual physics, gradient probes
/// and rejected trials share the original SolveControl.
///
/// In-memory run segments preserve the original AL state. KKT reports belong
/// to the latest accepted point; best_feasible may be earlier. Inventory is
/// DILATED, while stiffness and self-weight are ERODED. Incoming operator scales
/// are restored by each evaluation. There is no monotone feasibility-restoration
/// assumption: adding material also changes self-weight. A finite AL allowance
/// does not establish global design optimality or continuum feasibility.
pub struct AffineLoadSpaceDesign3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    dead: &'a [f64],
    basis: &'a [&'a [f64]],
    constraints: AffineLoadSpaceConstraints3<'a>,
    control: &'a mut SolveControl<'callback>,
    policy: AffineLoadSpace3,
    goal: AffineLoadSpaceDesignGoal3,
    options: AffineLoadSpaceDesignOptions3,
    check: AffineLoadSpaceGradientCheck3,
    state: ProjectedAlState,
    accepted: AffineLoadSpaceDesignPoint3,
    best: Option<AffineLoadSpaceDesignPoint3>,
    history: Vec<AffineLoadSpaceDesignIteration3>,
}
impl<'a, 'callback, O: Sdf3Elasticity> AffineLoadSpaceDesign3<'a, 'callback, O> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        rho: &[f64],
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::new_with_constraints(
            study,
            dead_offset,
            basis,
            rho,
            policy,
            goal,
            options,
            AffineLoadSpaceConstraints3::default(),
            control,
        )
    }

    /// Keep every signed displacement/support-force bound over the complete
    /// affine load ball in the SAME optimizer. The original compliance/volume
    /// constraint and all response limits form one conservative upper aggregate;
    /// raw physical checks additionally protect the feasible archive. Observation
    /// laws, bounds, numerical scales and smoothing stay fixed across all trials.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_response_limits(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        rho: &[f64],
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        responses: AffineResponseLimits3<'a>,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::new_with_constraints(
            study,
            dead_offset,
            basis,
            rho,
            policy,
            goal,
            options,
            AffineLoadSpaceConstraints3::from_responses(responses),
            control,
        )
    }

    /// Bind compliance/volume, signed response limits and native sampled
    /// strength to one existing AL owner. Every trial and fresh recheck carries
    /// the entire family; no component's adverse direction stands in for another.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_constraints(
        study: &'a mut CutDensityStudy3<O>,
        dead_offset: &'a [f64],
        basis: &'a [&'a [f64]],
        rho: &[f64],
        policy: AffineLoadSpace3,
        goal: AffineLoadSpaceDesignGoal3,
        options: AffineLoadSpaceDesignOptions3,
        constraints: AffineLoadSpaceConstraints3<'a>,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        admit(study, dead_offset, basis, rho, &policy, goal, options)?;
        constraints
            .admit(study)
            .map_err(ProjectedAlError::Evaluation)?;
        if let Some(limit) = constraints.stress {
            let params = study.params();
            let floor = crate::filter::heaviside(
                options.density_floor,
                params.beta,
                params.eta + options.robustness.threshold_offset,
            );
            if !crate::sdf3::design::monotone_stress_branch(
                floor,
                params,
                limit.options.relaxation_power,
            ) {
                return Err(ProjectedAlError::Invalid(
                    "eroded density floor permits ersatz stress foldback",
                ));
            }
        }
        let check = gradient_check(
            study,
            rho,
            dead_offset,
            basis,
            &policy,
            goal,
            options,
            constraints,
            control,
        )
        .map_err(ProjectedAlError::Evaluation)?;
        if !check.passed() {
            return Err(ProjectedAlError::Invalid(
                "affine objective or complete constraint gradient gate failed",
            ));
        }
        let mut last = None;
        let mut stopped = None;
        let n = study.cells();
        let result = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(
                rho,
                &vec![options.density_floor; n],
                &vec![1.0; n],
                options.optimizer,
                &mut |r| {
                    sample(
                        study,
                        r,
                        dead_offset,
                        basis,
                        &policy,
                        goal,
                        options,
                        constraints,
                        &mut last,
                        &mut **ledger.borrow_mut(),
                    )
                },
                |_| poll(&mut **ledger.borrow_mut(), &mut stopped),
            )
        };
        if let Some(e) = stopped {
            return Err(ProjectedAlError::Evaluation(e));
        }
        let state = result?;
        let accepted = last.expect("admitted affine sample retains complete physical fields");
        let first = row(&accepted, 0, state.sample(), goal, &policy);
        let best = first.feasible.then(|| accepted.clone());
        Ok(Self {
            study,
            dead: dead_offset,
            basis,
            constraints,
            control,
            policy,
            goal,
            options,
            check,
            state,
            accepted,
            best,
            history: vec![first],
        })
    }
    #[must_use]
    pub fn study(&self) -> &CutDensityStudy3<O> {
        self.study
    }
    #[must_use]
    pub fn policy(&self) -> &AffineLoadSpace3 {
        &self.policy
    }
    #[must_use]
    pub fn goal(&self) -> AffineLoadSpaceDesignGoal3 {
        self.goal
    }
    #[must_use]
    pub fn options(&self) -> AffineLoadSpaceDesignOptions3 {
        self.options
    }
    #[must_use]
    pub fn response_limits(&self) -> Option<AffineResponseLimits3<'a>> {
        self.constraints.responses
    }
    #[must_use]
    pub fn constraints(&self) -> AffineLoadSpaceConstraints3<'a> {
        self.constraints
    }
    #[must_use]
    pub fn accepted(&self) -> &AffineLoadSpaceDesignPoint3 {
        &self.accepted
    }
    #[must_use]
    pub fn best_feasible(&self) -> Option<&AffineLoadSpaceDesignPoint3> {
        self.best.as_ref()
    }
    #[must_use]
    pub fn history(&self) -> &[AffineLoadSpaceDesignIteration3] {
        &self.history
    }
    #[must_use]
    pub fn gradient_check(&self) -> &AffineLoadSpaceGradientCheck3 {
        &self.check
    }
    #[must_use]
    pub fn work(&self) -> SolveWork {
        self.control.work()
    }
    #[must_use]
    pub fn optimizer_work(&self) -> ProjectedAlWork {
        self.state.work()
    }

    pub fn run(
        &mut self,
        additional_steps: usize,
    ) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let mut stopped = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let dead = self.dead;
                let basis = self.basis;
                let policy = &self.policy;
                let goal = self.goal;
                let options = self.options;
                let constraints = self.constraints;
                self.state.try_run(
                    remaining.min(1),
                    &mut |r| {
                        sample(
                            study,
                            r,
                            dead,
                            basis,
                            policy,
                            goal,
                            options,
                            constraints,
                            &mut last,
                            &mut **ledger.borrow_mut(),
                        )
                    },
                    |_| poll(&mut **ledger.borrow_mut(), &mut stopped),
                )
            };
            if self.state.work().iterations > before {
                let accepted =
                    last.expect("accepted affine update retains complete physical fields");
                assert_eq!(
                    accepted.rho,
                    self.state.point(),
                    "affine physical fields detached from AL state"
                );
                self.accepted = accepted;
                let next = row(
                    &self.accepted,
                    self.state.work().iterations,
                    self.state.sample(),
                    self.goal,
                    &self.policy,
                );
                if next.feasible
                    && self.best.as_ref().is_none_or(|b| {
                        archive_objective(&self.accepted.evaluation, self.goal)
                            < archive_objective(&b.evaluation, self.goal)
                    })
                {
                    self.best = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control
                    .checkpoint("sdf3-affine-design-accepted")
                    .map_err(ProjectedAlError::Evaluation)?;
            }
            if let Some(e) = stopped {
                return Err(ProjectedAlError::Evaluation(e));
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 {
                return Ok(report);
            }
            remaining -= 1;
        }
    }

    /// Recompute the full affine load envelope at the retained feasible point,
    /// not only the previously maximizing direction. State/history stay intact.
    pub fn recheck_best(&mut self) -> Result<Option<AffineLoadSpaceDesignPoint3>, EvaluationStop> {
        let Some(best) = &self.best else {
            return Ok(None);
        };
        self.control.checkpoint("sdf3-affine-design-recheck")?;
        let evaluation = evaluate(
            self.study,
            &best.rho,
            self.dead,
            self.basis,
            &self.policy,
            self.goal,
            self.options,
            self.constraints,
            self.control,
        )?;
        if !physical_feasible(&evaluation, self.goal, &self.policy) {
            return Err(failure("sdf3-affine-design-recheck"));
        }
        Ok(Some(AffineLoadSpaceDesignPoint3 {
            rho: best.rho.clone(),
            evaluation,
        }))
    }
}
