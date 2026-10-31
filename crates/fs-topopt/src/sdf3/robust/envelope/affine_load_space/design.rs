//! Constrained design for a common load plus a continuous uncertain load set.
//! One existing projected-AL state owns both compliance and inventory modes.
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
    /// Actual full-chain worst-compliance AND selected-inventory check.
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

/// Numerical directional evidence at the initial design. The shared audit's
/// second scalar is named inventory here because it may be physical kg.
#[derive(Debug, Clone)]
pub struct AffineLoadSpaceGradientCheck3 {
    pub params: SimpParams,
    pub goal: AffineLoadSpaceDesignGoal3,
    pub options: GradientCheckOptions,
    pub baseline_maximum: f64,
    pub baseline_inventory: f64,
    pub maximum_relative_error: f64,
    pub inventory_relative_error: f64,
    pub directions: usize,
    pub work: SolveWork,
}
impl AffineLoadSpaceGradientCheck3 {
    #[must_use]
    pub fn passed(&self) -> bool {
        self.directions > 0
            && [self.maximum_relative_error, self.inventory_relative_error]
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
fn physical_feasible(
    e: &AffineLoadSpaceEvaluation3,
    goal: AffineLoadSpaceDesignGoal3,
    policy: &AffineLoadSpace3,
) -> bool {
    match goal {
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
    }
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
    last: &mut Option<AffineLoadSpaceDesignPoint3>,
    control: &mut SolveControl<'_>,
) -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let e = study.affine_load_space_with_mass(
        rho,
        dead,
        basis,
        options.robustness,
        policy,
        goal.needs_mass(),
        control,
    )?;
    // A valid worst-load witness can exist without a unique differentiable
    // envelope. Let the existing bounded line search reject that trial.
    let Some(gradient) = e.gradient.as_ref() else {
        return Ok(None);
    };
    let (material, material_gradient) = inventory(&e, goal)?;
    let sample = match goal {
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => ProjectedAlSample {
            objective: e.maximum,
            gradient: gradient.clone(),
            constraint: material - 1.0,
            constraint_gradient: material_gradient,
        },
        AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => ProjectedAlSample {
            objective: material,
            gradient: material_gradient,
            constraint: e.maximum - 1.0,
            constraint_gradient: gradient.clone(),
        },
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
    let invalid = ProjectedAlError::Invalid;
    options
        .optimizer
        .validate::<EvaluationStop>(study.cells())?;
    if !options.robustness.valid(study.params())
        || dead.len() != study.operator().n()
        || dead.iter().any(|f| !f.is_finite())
        || basis.len() != policy.dimension()
        || basis
            .iter()
            .any(|f| f.len() != dead.len() || f.iter().any(|v| !v.is_finite()))
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
    control: &mut SolveControl<'_>,
) -> Result<AffineLoadSpaceGradientCheck3, EvaluationStop> {
    let params = study.params();
    let numerical = audit_design_map(
        params,
        rho,
        options.gradient_check,
        control,
        |r, baseline, control| {
            let e = study.affine_load_space_with_mass(
                r,
                dead,
                basis,
                options.robustness,
                policy,
                goal.needs_mass(),
                control,
            )?;
            let gradient = match &e.gradient {
                Some(g) => g.clone(),
                None if !baseline => Vec::new(),
                None => return Err(failure("sdf3-affine-design-nondifferentiable-baseline")),
            };
            let (material, material_gradient) = inventory(&e, goal)?;
            Ok(GradientSample {
                compliance: e.maximum,
                compliance_gradient: gradient,
                volume: material,
                volume_gradient: material_gradient,
            })
        },
    )?;
    Ok(AffineLoadSpaceGradientCheck3 {
        params,
        goal,
        options: options.gradient_check,
        baseline_maximum: numerical.baseline_compliance,
        baseline_inventory: numerical.baseline_volume_fraction,
        maximum_relative_error: numerical
            .probes
            .iter()
            .map(|p| p.compliance_relative_error)
            .fold(0.0, f64::max),
        inventory_relative_error: numerical
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
        admit(study, dead_offset, basis, rho, &policy, goal, options)?;
        let check = gradient_check(
            study,
            rho,
            dead_offset,
            basis,
            &policy,
            goal,
            options,
            control,
        )
        .map_err(ProjectedAlError::Evaluation)?;
        if !check.passed() {
            return Err(ProjectedAlError::Invalid(
                "affine worst-load or inventory gradient gate failed",
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
        let evaluation = self.study.affine_load_space_with_mass(
            &best.rho,
            self.dead,
            self.basis,
            self.options.robustness,
            &self.policy,
            self.goal.needs_mass(),
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
