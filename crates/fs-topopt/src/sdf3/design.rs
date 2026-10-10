//! Minimum projected material volume under an explicit relaxed stress cap.
//!
//! Reuses fs-ascent's linear-storage projected augmented Lagrangian. Every
//! sample is the real retained raw-SDF elasticity/stress/adjoint evaluation;
//! there is no dense KKT system, surrogate, numerical-gradient optimizer or
//! separate optimizer implementation. Geometry and material model remain fixed.

use super::stress::{self, StressError3, StressEvaluation3, StressMeasure3, StressOptions3};
use super::{CutDensityStudy3, PhysicalRegion3, Sdf3Elasticity};
use crate::pipeline::LoadCase;
use crate::{SimpParams, SolveControl, SolveWork};
use fs_ascent::projected_al::{
    ProjectedAlCheckpoint, ProjectedAlError, ProjectedAlOptions, ProjectedAlReport, ProjectedAlSample, ProjectedAlState,
    ProjectedAlStop, ProjectedAlWork,
};
use std::{cell::RefCell, ops::ControlFlow};

#[derive(Debug, Clone, Copy)]
pub struct StressDesignOptions3 {
    pub stress: StressOptions3,
    /// Cap on the declared quadrature aggregate, in reference-material stress
    /// units. The default is a normalized average; an explicitly selected
    /// sampled-peak norm bounds the retained samples, never the continuum.
    pub stress_limit: f64,
    /// Exact raw-density box [density_floor, 1]. Prescribed physical regions
    /// are applied after the filter/projection map, so an authored void is zero.
    /// The projected lower bound of free material must remain above the
    /// ersatz-modulus qp turnover; a small floor or sharp projection is refused.
    pub density_floor: f64,
    /// The single inequality is dimensionless: aggregate/stress_limit - 1.
    /// Feasibility and numerical KKT tolerances use these coordinates.
    pub optimizer: ProjectedAlOptions,
}
impl Default for StressDesignOptions3 {
    fn default() -> Self {
        Self {
            stress: StressOptions3::default(),
            stress_limit: 1.0,
            density_floor: 0.05,
            optimizer: ProjectedAlOptions::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StressDesignIteration3 {
    pub iteration: usize,
    pub volume_fraction: f64,
    pub stress_aggregate: f64,
    pub sampled_relaxed_max: f64,
    pub sampled_physical_max: f64,
    /// Positive part of aggregate/stress_limit - 1.
    pub constraint_violation: f64,
    /// Feasibility of this aggregate at the declared numerical tolerance.
    pub feasible: bool,
}

/// Accepted optimizer state and the best feasible design, without serializing
/// an operator or trusting cached displacement/gradient fields after restart.
/// The consumer must bind this to the original geometry, loads and policy.
#[derive(Debug, Clone)]
pub struct StressDesignCheckpoint3 {
    pub optimizer: ProjectedAlCheckpoint,
    pub best_feasible_density: Option<Vec<f64>>,
    pub history: Vec<StressDesignIteration3>,
    /// Physical endpoint re-evaluations already included in optimizer work.
    pub restoration_evaluations: usize,
}
impl StressDesignCheckpoint3 {
    /// One accepted-design solve, and a second only for a distinct incumbent.
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(self.best_feasible_density.as_ref()
            .is_some_and(|rho| rho != &self.optimizer.point))
    }
}
fn row(
    e: &StressEvaluation3,
    iteration: usize,
    options: StressDesignOptions3,
) -> StressDesignIteration3 {
    let violation = (e.aggregate / options.stress_limit - 1.0).max(0.0);
    StressDesignIteration3 {
        iteration,
        volume_fraction: e.volume_fraction,
        stress_aggregate: e.aggregate,
        sampled_relaxed_max: e.sampled_relaxed_max,
        sampled_physical_max: e.sampled_physical_max,
        constraint_violation: violation,
        feasible: violation <= options.optimizer.tolerance,
    }
}

// Under a homogeneous force-controlled density change, relaxed stress is
// proportional to r^q/(e_min+(1-e_min)r^p). Below this turning point it decreases
// when material is removed, because the ersatz floor carries the load. Such a
// branch cannot be admitted as a strength-constrained minimum-volume design.
fn monotone_stress_branch(r: f64, params: SimpParams, q: f64) -> bool {
    let turnover_power = q * params.e_min / ((params.penal - q) * (1.0 - params.e_min));
    r.is_finite() && r > 0.0 && fs_math::det::pow(r, params.penal) > turnover_power
}
fn sample<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[LoadCase<'_>],
    rho: &[f64],
    options: StressDesignOptions3,
    measure: StressMeasure3,
    last: &mut Option<StressEvaluation3>,
    control: &mut SolveControl<'_>,
) -> Result<Option<ProjectedAlSample>, StressError3> {
    control.checkpoint("sdf3-stress-design-evaluate")?;
    let evaluated = study.evaluate_stress_with_measure(rho, loads, options.stress, measure, control)?;
    // The graph filter's exact maximum principle preserves the raw box. Check
    // the actually solved/projection values too, rather than assume numerical
    // filter error cannot cross the strength-model boundary. An authored
    // physical void is an immutable zero in the material map, not an optimizer
    // step along the ersatz foldback branch. No free density gets this exception.
    if evaluated.projected_rho.iter().enumerate().any(|(i, r)| {
        let prescribed_void = *r == 0.0
            && study.physical_regions()
                .is_some_and(|regions| regions[i] == PhysicalRegion3::Void);
        !prescribed_void
            && !monotone_stress_branch(*r, study.params(), options.stress.relaxation_power)
    })
    {
        return Err(StressError3::Invalid(
            "projected density entered the ersatz stress turnover",
        ));
    }
    let sample = ProjectedAlSample {
        objective: evaluated.volume_fraction,
        gradient: evaluated.volume_gradient.clone(),
        constraint: evaluated.aggregate / options.stress_limit - 1.0,
        constraint_gradient: evaluated
            .gradient
            .iter()
            .map(|g| g / options.stress_limit)
            .collect(),
    };
    *last = Some(evaluated);
    Ok(Some(sample))
}

/// In-memory resumable constrained design, bound to one geometry, force family,
/// SIMP/filter/projection model, stress measure and cumulative work ledger.
/// Accepted AL steps need not be feasible or decrease original material volume.
/// `best_feasible()` retains the least-volume accepted feasible design separately
/// so a budget/cancellation stop does not discard a previously useful result.
///
/// Trial evaluations always restore scales. Only accepted optimizer points are
/// installed; a later failed poll preserves their matching complete fields and
/// counters. A returned error leaves this session available for inspection and
/// continuation. This numerical KKT system does not prove global optimality,
/// continuum stress safety, manufacturability, or a certified geometric boundary.
pub struct StressDesignStudy3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    loads: &'a [LoadCase<'a>],
    control: &'a mut SolveControl<'callback>,
    options: StressDesignOptions3,
    measure: StressMeasure3,
    state: ProjectedAlState,
    accepted: StressEvaluation3,
    best_feasible: Option<StressEvaluation3>,
    history: Vec<StressDesignIteration3>,
    restoration_evaluations: usize,
}
impl<'a, 'callback, O: Sdf3Elasticity> StressDesignStudy3<'a, 'callback, O> {
    pub fn new(
        study: &'a mut CutDensityStudy3<O>,
        loads: &'a [LoadCase<'a>],
        rho: &[f64],
        options: StressDesignOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<StressError3>> {
        Self::new_with_measure(study, loads, rho, options, StressMeasure3::NormalizedAverage, control)
    }

    /// Start with an explicit stress functional. This changes only the scalar
    /// constraint and its exact adjoint, not the AL optimizer or density bounds.
    /// Keep this measure bound to the source when persisting `checkpoint()`.
    pub fn new_with_measure(
        study: &'a mut CutDensityStudy3<O>,
        loads: &'a [LoadCase<'a>],
        rho: &[f64],
        options: StressDesignOptions3,
        measure: StressMeasure3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<StressError3>> {
        let n = study.cells();
        options.optimizer.validate(n)?;
        if !options.stress_limit.is_finite()
            || options.stress_limit <= 0.0
            || !options.density_floor.is_finite()
            || options.density_floor <= 0.0
            || options.density_floor >= 1.0
            || rho.len() != n
            || rho
                .iter()
                .any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0)
        {
            return Err(ProjectedAlError::Invalid(
                "invalid stress design cap, bounds or initial density",
            ));
        }
        stress::admit(study, rho, loads, options.stress).map_err(ProjectedAlError::Evaluation)?;
        let params = study.params();
        let projected_floor =
            crate::filter::heaviside(options.density_floor, params.beta, params.eta);
        if !monotone_stress_branch(projected_floor, params, options.stress.relaxation_power) {
            return Err(ProjectedAlError::Invalid(
                "density floor and projection permit ersatz stress foldback",
            ));
        }
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(
                rho,
                &vec![options.density_floor; n],
                &vec![1.0; n],
                options.optimizer,
                &mut |rho| {
                    sample(
                        study,
                        loads,
                        rho,
                        options,
                        measure,
                        &mut last,
                        &mut **ledger.borrow_mut(),
                    )
                },
                |_| {
                    if ledger
                        .borrow_mut()
                        .checkpoint("sdf3-stress-design-control")
                        .is_ok()
                    {
                        ControlFlow::Continue(())
                    } else {
                        ControlFlow::Break(())
                    }
                },
            )?
        };
        control
            .checkpoint("sdf3-stress-design-initialize")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        let accepted = last.expect("initial optimizer sample contains complete stress evidence");
        study
            .operator
            .set_scales(&accepted.scales)
            .expect("evaluated scales are valid");
        let first = row(&accepted, 0, options);
        let best_feasible = first.feasible.then(|| accepted.clone());
        Ok(Self {
            study,
            loads,
            control,
            options,
            measure,
            state,
            accepted,
            best_feasible,
            history: vec![first],
            restoration_evaluations: 0,
        })
    }

    /// Rebuild accepted and incumbent fields from the SAME admitted physical
    /// problem, then restore spectral steps, multipliers and spent work exactly.
    /// Endpoint checks cost one or two real stress evaluations, included in the
    /// original evaluation allowance; all Krylov work uses the supplied control.
    /// Failed admission restores incoming scales and changes no checkpoint.
    pub fn restore(
        study: &'a mut CutDensityStudy3<O>,
        loads: &'a [LoadCase<'a>],
        checkpoint: StressDesignCheckpoint3,
        options: StressDesignOptions3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<StressError3>> {
        Self::restore_with_measure(study, loads, checkpoint, options, StressMeasure3::NormalizedAverage, control)
    }

    /// Restore using the original explicitly bound measure. As with geometry,
    /// loads and the other stress options, the checkpoint consumer must retain
    /// this policy. Both accepted and distinct incumbent endpoints are freshly
    /// evaluated with it; a changed physical response rejects restoration.
    pub fn restore_with_measure(
        study: &'a mut CutDensityStudy3<O>,
        loads: &'a [LoadCase<'a>],
        checkpoint: StressDesignCheckpoint3,
        options: StressDesignOptions3,
        measure: StressMeasure3,
        control: &'a mut SolveControl<'callback>,
    ) -> Result<Self, ProjectedAlError<StressError3>> {
        let n = study.cells();
        let cost = checkpoint.restoration_cost();
        let restoration_evaluations = checkpoint.restoration_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("stress restoration counter overflow"))?;
        let mut optimizer = checkpoint.optimizer;
        let previous_evaluations = optimizer.work.evaluations;
        optimizer.work.evaluations = previous_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("stress restoration evaluation overflow"))?;
        let restored = ProjectedAlState::try_restore(
            optimizer, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
        )?;
        let history = checkpoint.history;
        if checkpoint.restoration_evaluations >= previous_evaluations
            || history.len() != restored.work().iterations + 1
            || history.iter().enumerate().any(|(i, r)| {
                r.iteration != i || !r.volume_fraction.is_finite()
                    || !(0.0..=1.0).contains(&r.volume_fraction)
                    || [r.stress_aggregate, r.sampled_relaxed_max, r.sampled_physical_max,
                        r.constraint_violation].iter().any(|v| !v.is_finite() || *v < 0.0)
                    || r.constraint_violation != (r.stress_aggregate / options.stress_limit - 1.0).max(0.0)
                    || r.feasible != (r.constraint_violation <= options.optimizer.tolerance)
            })
        {
            return Err(ProjectedAlError::Invalid("invalid stress checkpoint history or work"));
        }
        let best_row = history.iter().filter(|r| r.feasible).fold(None,
            |best: Option<&StressDesignIteration3>, r| {
                if best.is_none_or(|b| r.volume_fraction < b.volume_fraction) { Some(r) } else { best }
            });
        if best_row.is_some() != checkpoint.best_feasible_density.is_some() {
            return Err(ProjectedAlError::Invalid("stress checkpoint lost its feasible incumbent"));
        }
        let previous_scales = study.operator.scales().to_vec();
        let mut result = Self::new_with_measure(study, loads, restored.point(), options, measure, control)?;
        let rebuilt = (|| {
            if result.state.sample() != restored.sample()
                || history.last() != Some(&row(&result.accepted, restored.work().iterations, options))
            {
                return Err(ProjectedAlError::Invalid("restored stress point changed its physical response"));
            }
            let best = match checkpoint.best_feasible_density {
                None => None,
                Some(rho) if rho == result.accepted.rho => Some(result.accepted.clone()),
                Some(rho) => {
                    if rho.len() != n || rho.iter().any(|r| !r.is_finite()
                        || !(options.density_floor..=1.0).contains(r)) {
                        return Err(ProjectedAlError::Invalid("invalid retained feasible density"));
                    }
                    let mut evaluated = None;
                    sample(result.study, loads, &rho, options, measure, &mut evaluated, result.control)
                        .map_err(ProjectedAlError::Evaluation)?;
                    evaluated
                }
            };
            if let (Some(evaluated), Some(expected)) = (&best, best_row) {
                if row(evaluated, expected.iteration, options) != *expected {
                    return Err(ProjectedAlError::Invalid("retained feasible design changed its physical response"));
                }
            }
            Ok(best)
        })();
        let best = match rebuilt {
            Ok(best) => best,
            Err(error) => {
                result.study.operator.set_scales(&previous_scales).expect("incoming scales remain valid");
                return Err(error);
            }
        };
        result.state = restored;
        result.best_feasible = best;
        result.history = history;
        result.restoration_evaluations = restoration_evaluations;
        Ok(result)
    }

    #[must_use]
    pub fn checkpoint(&self) -> StressDesignCheckpoint3 {
        StressDesignCheckpoint3 {
            optimizer: self.state.checkpoint(),
            best_feasible_density: self.best_feasible.as_ref().map(|e| e.rho.clone()),
            history: self.history.clone(),
            restoration_evaluations: self.restoration_evaluations,
        }
    }
    /// Functional used by initialization, trials, adjoints and restored fields.
    #[must_use]
    pub fn stress_measure(&self) -> StressMeasure3 { self.measure }
    #[must_use]
    pub fn accepted(&self) -> &StressEvaluation3 {
        &self.accepted
    }
    #[must_use]
    pub fn best_feasible(&self) -> Option<&StressEvaluation3> {
        self.best_feasible.as_ref()
    }
    #[must_use]
    pub fn study(&self) -> &CutDensityStudy3<O> {
        self.study
    }
    #[must_use]
    pub fn point(&self) -> &[f64] {
        self.state.point()
    }
    #[must_use]
    pub fn history(&self) -> &[StressDesignIteration3] {
        &self.history
    }
    #[must_use]
    pub fn work(&self) -> SolveWork {
        self.control.work()
    }
    #[must_use]
    pub fn optimizer_work(&self) -> ProjectedAlWork {
        self.state.work()
    }
    #[must_use]
    pub fn constraint_violation(&self) -> f64 {
        self.state.sample().constraint.max(0.0)
    }
    #[must_use]
    pub fn feasible(&self) -> bool {
        self.constraint_violation() <= self.options.optimizer.tolerance
    }

    /// Advance the existing AL state. Zero steps use cached derivatives without
    /// repeating physics or changing multipliers. Segments preserve spectral and
    /// dual state; rejected trials remain charged to both cumulative budgets.
    pub fn run(
        &mut self,
        additional_steps: usize,
    ) -> Result<ProjectedAlReport, ProjectedAlError<StressError3>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let loads = self.loads;
                let options = self.options;
                let measure = self.measure;
                self.state.try_run(
                    remaining.min(1),
                    &mut |rho| {
                        sample(
                            study,
                            loads,
                            rho,
                            options,
                            measure,
                            &mut last,
                            &mut **ledger.borrow_mut(),
                        )
                    },
                    |_| {
                        if ledger
                            .borrow_mut()
                            .checkpoint("sdf3-stress-design-control")
                            .is_ok()
                        {
                            ControlFlow::Continue(())
                        } else {
                            ControlFlow::Break(())
                        }
                    },
                )
            };
            if self.state.work().iterations > before {
                let accepted =
                    last.expect("accepted optimizer point contains complete stress evidence");
                assert_eq!(
                    accepted.rho,
                    self.state.point(),
                    "accepted stress field and optimizer point must match"
                );
                self.study
                    .operator
                    .set_scales(&accepted.scales)
                    .expect("accepted scales remain valid");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations, self.options);
                if next.feasible
                    && self
                        .best_feasible
                        .as_ref()
                        .is_none_or(|best| self.accepted.volume_fraction < best.volume_fraction)
                {
                    self.best_feasible = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control
                    .checkpoint("sdf3-stress-design-accepted")
                    .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 {
                return Ok(report);
            }
            remaining -= 1;
        }
    }
}
