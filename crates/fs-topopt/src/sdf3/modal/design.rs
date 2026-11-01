//! Frequency-constrained design on the same material-bound CutFEM pencil.
//! The existing projected augmented-Lagrangian owner supplies bounded line
//! searches, KKT diagnostics and resumable iteration state; this module supplies
//! physical objectives and their full filter/projection/material derivatives.

use super::{Sdf3EigenfrequencyEvaluation3, Sdf3EigenfrequencyOptions3};
use crate::gradient_check::{GradientSample, audit_design_map};
use crate::modal::MatrixFreeEigenError;
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity, failure};
use crate::{EvaluationStop, GradientCheckOptions, SolveControl, SolveWork};
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlOptions, ProjectedAlReport, ProjectedAlSample,
    ProjectedAlState, ProjectedAlStop, ProjectedAlWork,
};
use fs_cutfem::elastic3::MaterialMass3;
use fs_spectral::generalized::GeneralizedError;
use std::{cell::RefCell, ops::ControlFlow};

/// Two formulations of the same free-vibration design problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Sdf3EigenfrequencyGoal3 {
    /// Maximize the smooth minimum of squared angular frequencies, subject to
    /// the original physical volume fraction cap. The eigenvalue reference is
    /// a numerical scale here, not a minimum-frequency requirement.
    MaximumFrequency { volume_limit: f64 },
    /// Minimize physical volume or kg while the smooth minimum AND the lowest
    /// returned eigenvalue remain at least the eigenvalue reference.
    /// Physical inventory uses linear occupancy, independently of the modal
    /// low-density mass interpolation.
    MinimumMaterial {
        objective: MinimumMaterialObjective3,
        volume_scale: f64,
    },
}

#[derive(Debug, Clone, Copy)]
pub struct Sdf3EigenfrequencyDesignOptions3 {
    pub frequency: Sdf3EigenfrequencyOptions3,
    /// Positive squared angular frequency [s^-2] for SI geometry/materials.
    /// A minimum frequency f in Hz corresponds to (2*pi*f)^2.
    pub eigenvalue_reference: f64,
    pub density_floor: f64,
    /// Numerical full-chain check of BOTH the actual objective and constraint.
    pub gradient_check: GradientCheckOptions,
    pub optimizer: ProjectedAlOptions,
}
impl Default for Sdf3EigenfrequencyDesignOptions3 {
    fn default() -> Self {
        Self {
            frequency: Sdf3EigenfrequencyOptions3::default(),
            eigenvalue_reference: 1.0,
            density_floor: 1e-3,
            gradient_check: GradientCheckOptions::default(),
            optimizer: ProjectedAlOptions::default(),
        }
    }
}

/// Directional numerical evidence for the bound optimization formulation.
#[derive(Debug, Clone)]
pub struct Sdf3EigenfrequencyGradientCheck3 {
    pub goal: Sdf3EigenfrequencyGoal3,
    pub eigenvalue_reference: f64,
    pub options: GradientCheckOptions,
    pub baseline_objective: f64,
    pub baseline_constraint: f64,
    pub objective_relative_error: f64,
    pub constraint_relative_error: f64,
    pub directions: usize,
    pub work: SolveWork,
}
impl Sdf3EigenfrequencyGradientCheck3 {
    #[must_use]
    pub fn passed(&self) -> bool {
        self.directions > 0
            && [self.objective_relative_error, self.constraint_relative_error]
                .iter()
                .all(|e| e.is_finite() && *e >= 0.0 && *e <= self.options.relative_tolerance)
    }
}

#[derive(Debug, Clone)]
pub struct Sdf3EigenfrequencyDesignPoint3 {
    pub rho: Vec<f64>,
    pub evaluation: Sdf3EigenfrequencyEvaluation3,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Sdf3EigenfrequencyDesignIteration3 {
    pub iteration: usize,
    pub objective: f64,
    pub constraint: f64,
    pub aggregate: f64,
    pub lowest_eigenvalue: f64,
    pub volume_fraction: f64,
    pub physical_mass_kg: f64,
    /// Checks the original physical requirement independently of AL tolerance.
    pub feasible: bool,
}

fn normalized(value: f64, scale: f64) -> Result<f64, MatrixFreeEigenError> {
    let result = value / scale;
    if result.is_finite() && (value == 0.0 || result != 0.0) {
        Ok(result)
    } else {
        Err(GeneralizedError::Breakdown("unrepresentable modal design normalization"))
    }
}
fn inventory(
    e: &Sdf3EigenfrequencyEvaluation3,
    goal: Sdf3EigenfrequencyGoal3,
) -> (f64, &[f64], f64) {
    match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } => {
            (e.volume_fraction, &e.volume_gradient, volume_limit)
        }
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction, volume_scale,
        } => (e.volume_fraction, &e.volume_gradient, volume_scale),
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg }, ..
        } => (e.physical_mass_kg, &e.physical_mass_gradient, scale_kg),
    }
}
fn optimization_sample(
    e: &Sdf3EigenfrequencyEvaluation3,
    goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyDesignOptions3,
) -> Result<ProjectedAlSample, MatrixFreeEigenError> {
    let (material, derivative, scale) = inventory(e, goal);
    let material = normalized(material, scale)?;
    let material_gradient = derivative.iter()
        .map(|&g| normalized(g, scale)).collect::<Result<Vec<_>, _>>()?;
    let frequency = normalized(e.objective.aggregate, options.eigenvalue_reference)?;
    let negative_frequency_gradient = e.objective.gradient.iter()
        .map(|&g| normalized(-g, options.eigenvalue_reference))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } => ProjectedAlSample {
            objective: -frequency,
            gradient: negative_frequency_gradient,
            constraint: material - 1.0,
            constraint_gradient: material_gradient,
        },
        Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } => ProjectedAlSample {
            objective: material,
            gradient: material_gradient,
            constraint: 1.0 - frequency,
            constraint_gradient: negative_frequency_gradient,
        },
    })
}

fn physical_feasible(
    e: &Sdf3EigenfrequencyEvaluation3,
    goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyDesignOptions3,
) -> bool {
    let original = match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } =>
            e.volume_fraction.is_finite() && e.volume_fraction <= volume_limit,
        Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } =>
            e.objective.aggregate.is_finite()
                && e.objective.aggregate >= options.eigenvalue_reference
                && e.objective.modal.values.first().is_some_and(|v|
                    v.is_finite() && *v >= options.eigenvalue_reference),
    };
    original && optimization_sample(e, goal, options)
        .is_ok_and(|sample| sample.constraint <= 0.0)
}
fn archive_objective(e: &Sdf3EigenfrequencyEvaluation3, goal: Sdf3EigenfrequencyGoal3) -> f64 {
    match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } => -e.objective.aggregate,
        Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } => inventory(e, goal).0,
    }
}
fn row(
    point: &Sdf3EigenfrequencyDesignPoint3, iteration: usize,
    sample: &ProjectedAlSample, goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyDesignOptions3,
) -> Sdf3EigenfrequencyDesignIteration3 {
    let e = &point.evaluation;
    Sdf3EigenfrequencyDesignIteration3 {
        iteration, objective: sample.objective, constraint: sample.constraint,
        aggregate: e.objective.aggregate, lowest_eigenvalue: e.objective.modal.values[0],
        volume_fraction: e.volume_fraction, physical_mass_kg: e.physical_mass_kg,
        feasible: physical_feasible(e, goal, options),
    }
}

fn admit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>, rho: &[f64], goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyDesignOptions3,
) -> Result<(), ProjectedAlError<MatrixFreeEigenError>> {
    options.optimizer.validate(study.cells())?;
    let scale = inventory_scale(goal);
    if rho.len() != study.cells()
        || !options.density_floor.is_finite()
        || options.density_floor <= 0.0 || options.density_floor >= 1.0
        || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0)
        || [options.eigenvalue_reference, scale].iter().any(|s|
            !s.is_finite() || *s <= 0.0 || !(1.0 / s).is_finite())
        || !options.gradient_check.step.is_finite()
        || options.gradient_check.step <= 0.0 || options.gradient_check.step > 0.125
        || !options.gradient_check.relative_tolerance.is_finite()
        || options.gradient_check.relative_tolerance <= 0.0
        || options.gradient_check.relative_tolerance >= 1.0
    {
        return Err(ProjectedAlError::Invalid("invalid modal design inputs or numerical scales"));
    }
    if let Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } = goal {
        if volume_limit > 1.0 || study.prescribed_solid_fraction() > volume_limit {
            return Err(ProjectedAlError::Invalid("invalid or unattainable modal design volume cap"));
        }
    }
    Ok(())
}
fn inventory_scale(goal: Sdf3EigenfrequencyGoal3) -> f64 {
    match goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } => volume_limit,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction, volume_scale,
        } => volume_scale,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg }, ..
        } => scale_kg,
    }
}

fn gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, mass: &MaterialMass3, rho: &[f64],
    goal: Sdf3EigenfrequencyGoal3, options: Sdf3EigenfrequencyDesignOptions3,
    control: &mut SolveControl<'_>,
) -> Result<Sdf3EigenfrequencyGradientCheck3, MatrixFreeEigenError> {
    let mut spectral_failure = None;
    let numerical = audit_design_map(
        study.params(), rho, options.gradient_check, control,
        |r, _, control| {
            let sample = study.evaluate_eigenfrequency(mass, r, options.frequency, control)
                .and_then(|e| {
                    let sample = optimization_sample(&e, goal, options)?;
                    // Audit the varying ratio before adding/subtracting one,
                    // while it still retains small resolvable variations.
                    let constraint = match goal {
                        Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } => {
                            let (value, _, scale) = inventory(&e, goal);
                            normalized(value, scale)?
                        }
                        Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } =>
                            normalized(-e.objective.aggregate, options.eigenvalue_reference)?,
                    };
                    Ok((sample, constraint))
                });
            match sample {
                Ok((s, constraint)) => {
                    Ok(GradientSample {
                        compliance: s.objective, compliance_gradient: s.gradient,
                        volume: constraint, volume_gradient: s.constraint_gradient,
                    })
                }
                Err(e) => {
                    spectral_failure = Some(e);
                    Err(failure("sdf3-modal-design-gradient"))
                }
            }
        },
    );
    if let Some(error) = spectral_failure { return Err(error); }
    let numerical = numerical.map_err(GeneralizedError::Operator)?;
    Ok(Sdf3EigenfrequencyGradientCheck3 {
        goal, eigenvalue_reference: options.eigenvalue_reference,
        options: options.gradient_check,
        baseline_objective: numerical.baseline_compliance,
        baseline_constraint: numerical.baseline_volume_fraction
            + match goal { Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } => -1.0, _ => 1.0 },
        objective_relative_error: numerical.probes.iter()
            .map(|p| p.compliance_relative_error).fold(0.0, f64::max),
        constraint_relative_error: numerical.probes.iter()
            .map(|p| p.volume_relative_error).fold(0.0, f64::max),
        directions: numerical.probes.len(), work: numerical.work,
    })
}

fn sample<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, mass: &MaterialMass3, rho: &[f64],
    goal: Sdf3EigenfrequencyGoal3, options: Sdf3EigenfrequencyDesignOptions3,
    last: &mut Option<Sdf3EigenfrequencyDesignPoint3>, control: &mut SolveControl<'_>,
) -> Result<Option<ProjectedAlSample>, MatrixFreeEigenError> {
    let evaluation = match study.evaluate_eigenfrequency(mass, rho, options.frequency, control) {
        Ok(evaluation) => evaluation,
        // A separated baseline is required, but a trial can cross a spectral
        // cluster or the piecewise mass-law switch. Let the bounded line search
        // reject that trial without inventing a derivative.
        Err(error) if super::nondifferentiable(&error) => return Ok(None),
        Err(error) => return Err(error),
    };
    let sample = optimization_sample(&evaluation, goal, options)?;
    *last = Some(Sdf3EigenfrequencyDesignPoint3 { rho: rho.to_vec(), evaluation });
    Ok(Some(sample))
}
fn poll(control: &mut SolveControl<'_>, stopped: &mut Option<EvaluationStop>) -> ControlFlow<()> {
    match control.checkpoint("sdf3-modal-design-control") {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => { *stopped = Some(error); ControlFlow::Break(()) }
    }
}

/// A fixed-grid frequency design run. Every accepted point, finite-difference
/// probe and rejected line-search trial uses the same native material field,
/// consistent mass cache, support map, filter and frequency aggregation.
///
/// Modes are infinitesimal free vibrations about the undeformed state; external
/// loads do not create geometric stiffness or prestress in this formulation.
/// The selected low cluster must be separated from an omitted mode. Numerical
/// residuals and a checked discrete gradient do not certify spectral completeness,
/// global design optimality, or continuum frequency limits.
///
/// Evaluation restores incoming operator scales on success and returned errors.
/// The retained best feasible point is checked against original physical bounds,
/// independently of the optimizer's KKT tolerance.
pub struct Sdf3EigenfrequencyDesign3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    mass: &'a MaterialMass3,
    control: &'a mut SolveControl<'callback>,
    goal: Sdf3EigenfrequencyGoal3,
    options: Sdf3EigenfrequencyDesignOptions3,
    check: Sdf3EigenfrequencyGradientCheck3,
    state: ProjectedAlState,
    accepted: Sdf3EigenfrequencyDesignPoint3,
    best: Option<Sdf3EigenfrequencyDesignPoint3>,
    history: Vec<Sdf3EigenfrequencyDesignIteration3>,
}
impl<'a, 'callback, O: Sdf3Elasticity> Sdf3EigenfrequencyDesign3<'a, 'callback, O> {
    pub fn new(
        study: &'a mut CutDensityStudy3<O>, mass: &'a MaterialMass3, rho: &[f64],
        goal: Sdf3EigenfrequencyGoal3, options: Sdf3EigenfrequencyDesignOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<MatrixFreeEigenError>> {
        admit(study, rho, goal, options)?;
        let check = gradient_check(study, mass, rho, goal, options, control)
            .map_err(ProjectedAlError::Evaluation)?;
        if !check.passed() {
            return Err(ProjectedAlError::Invalid("modal objective or constraint gradient gate failed"));
        }
        let mut last = None;
        let mut stopped = None;
        let n = study.cells();
        let result = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(
                rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |r| sample(study, mass, r, goal, options, &mut last, &mut **ledger.borrow_mut()),
                |_| poll(&mut **ledger.borrow_mut(), &mut stopped),
            )
        };
        if let Some(error) = stopped {
            return Err(ProjectedAlError::Evaluation(GeneralizedError::Operator(error)));
        }
        let state = result?;
        let accepted = last.expect("admitted modal sample retains its physical evaluation");
        let first = row(&accepted, 0, state.sample(), goal, options);
        let best = first.feasible.then(|| accepted.clone());
        Ok(Self { study, mass, control, goal, options, check, state, accepted, best,
            history: vec![first] })
    }

    #[must_use]
    pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use]
    pub fn goal(&self) -> Sdf3EigenfrequencyGoal3 { self.goal }
    #[must_use]
    pub fn options(&self) -> Sdf3EigenfrequencyDesignOptions3 { self.options }
    #[must_use]
    pub fn accepted(&self) -> &Sdf3EigenfrequencyDesignPoint3 { &self.accepted }
    #[must_use]
    pub fn best_feasible(&self) -> Option<&Sdf3EigenfrequencyDesignPoint3> { self.best.as_ref() }
    #[must_use]
    pub fn history(&self) -> &[Sdf3EigenfrequencyDesignIteration3] { &self.history }
    #[must_use]
    pub fn gradient_check(&self) -> &Sdf3EigenfrequencyGradientCheck3 { &self.check }
    #[must_use]
    pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use]
    pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }

    pub fn run(&mut self, additional_steps: usize)
        -> Result<ProjectedAlReport, ProjectedAlError<MatrixFreeEigenError>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let mut stopped = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let mass = self.mass;
                let goal = self.goal;
                let options = self.options;
                self.state.try_run(
                    remaining.min(1),
                    &mut |r| sample(study, mass, r, goal, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| poll(&mut **ledger.borrow_mut(), &mut stopped),
                )
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted modal update retains complete physical fields");
                assert_eq!(accepted.rho, self.state.point(), "modal fields detached from AL state");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations,
                    self.state.sample(), self.goal, self.options);
                if next.feasible && self.best.as_ref().is_none_or(|best|
                    archive_objective(&self.accepted.evaluation, self.goal)
                        < archive_objective(&best.evaluation, self.goal)) {
                    self.best = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("sdf3-modal-design-accepted")
                    .map_err(|e| ProjectedAlError::Evaluation(GeneralizedError::Operator(e)))?;
            }
            if let Some(error) = stopped {
                return Err(ProjectedAlError::Evaluation(GeneralizedError::Operator(error)));
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 {
                return Ok(report);
            }
            remaining -= 1;
        }
    }

    /// Freshly solve the complete selected modal cluster at the retained best
    /// feasible design. The existing optimizer state and history remain intact.
    pub fn recheck_best(&mut self) -> Result<Option<Sdf3EigenfrequencyDesignPoint3>, MatrixFreeEigenError> {
        let Some(best) = &self.best else { return Ok(None); };
        let evaluation = self.study.evaluate_eigenfrequency(
            self.mass, &best.rho, self.options.frequency, self.control)?;
        if !physical_feasible(&evaluation, self.goal, self.options) {
            return Err(GeneralizedError::Breakdown("modal design failed fresh physical recheck"));
        }
        Ok(Some(Sdf3EigenfrequencyDesignPoint3 { rho: best.rho.clone(), evaluation }))
    }
}
