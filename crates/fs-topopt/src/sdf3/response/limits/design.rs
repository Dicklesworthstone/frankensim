//! Minimum material subject to all declared force/displacement inequalities.
//! Reuse the existing linear-storage projected AL, not a second optimizer.
use super::*;
use std::cell::RefCell;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions, ProjectedAlReport,
    ProjectedAlSample, ProjectedAlState, ProjectedAlStop, ProjectedAlWork};

mod checkpoint;
pub use checkpoint::{MinimumMaterialCheckpoint3, MinimumMaterialEndpoint3, MinimumMaterialLimit3};
mod wire;
pub use wire::{MinimumMaterialCodecBudget3, MinimumMaterialCodecError3};

/// Material inventory minimized while retaining every original response limit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MinimumMaterialObjective3 {
    /// Projected geometric volume fraction, normalized by `volume_scale`.
    VolumeFraction,
    /// Physical inventory from each native source material's reference density
    /// and cut-cell volume. `scale_kg` is a positive numerical objective scale;
    /// reported mass remains in kg. Legacy E/nu-only operators have no density
    /// authority and cannot use this objective.
    MassKg { scale_kg: f64 },
}

#[derive(Debug, Clone, Copy)]
pub struct MinimumMaterialOptions3 {
    pub limits: ResponseLimitOptions3,
    /// Raw controls remain exactly in [density_floor,1]. Physical protected
    /// solids/voids and the filter/projection are governed by the original study.
    pub density_floor: f64,
    /// Positive fixed scale for the numerical volume-fraction objective.
    /// It never rescales reported physical volume or individual response limits.
    pub volume_scale: f64,
    pub objective: MinimumMaterialObjective3,
    pub optimizer: ProjectedAlOptions,
}
impl Default for MinimumMaterialOptions3 {
    fn default() -> Self {
        Self { limits: ResponseLimitOptions3::default(), density_floor: 1e-3,
            volume_scale: 1.0, objective: MinimumMaterialObjective3::VolumeFraction,
            optimizer: ProjectedAlOptions::default() }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct MinimumMaterialIteration3 {
    pub iteration: usize,
    pub volume_fraction: f64,
    /// Present for a physical-mass objective, in kg before normalization.
    pub mass_kg: Option<f64>,
    /// The normalized objective actually supplied to the optimizer.
    pub objective: f64,
    pub constraint: f64,
    pub worst_residual: f64,
}
fn row(e: &ResponseLimitEvaluation3, iteration: usize, objective: f64) -> MinimumMaterialIteration3 {
    MinimumMaterialIteration3 { iteration, volume_fraction: e.volume_fraction,
        mass_kg: e.mass_kg, objective,
        constraint: e.constraint, worst_residual: e.worst_residual }
}
fn objective_inventory(volume: f64, mass_kg: Option<f64>, objective: MinimumMaterialObjective3) -> Option<f64> {
    let value = match objective {
        MinimumMaterialObjective3::VolumeFraction if mass_kg.is_none() => volume,
        MinimumMaterialObjective3::MassKg { .. } => mass_kg?,
        _ => return None,
    };
    (value.is_finite() && value >= 0.0).then_some(value)
}
fn objective_value(volume: f64, mass_kg: Option<f64>, options: MinimumMaterialOptions3) -> Option<f64> {
    let inventory = objective_inventory(volume, mass_kg, options.objective)?;
    let scale = match options.objective {
        MinimumMaterialObjective3::VolumeFraction => options.volume_scale,
        MinimumMaterialObjective3::MassKg { scale_kg } => scale_kg,
    };
    let value = inventory/scale;
    (value.is_finite() && (inventory == 0.0 || value != 0.0)).then_some(value)
}
// Admission also serves reference-load refits, before invoking any callback.
pub(super) fn admit_objective<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>,
    options: MinimumMaterialOptions3) -> Result<(), ProjectedAlError<ResponseError3>> {
    if let MinimumMaterialObjective3::MassKg { scale_kg } = options.objective {
        if !scale_kg.is_finite() || scale_kg <= 0.0 || !(1.0/scale_kg).is_finite() {
            return Err(ProjectedAlError::Invalid("mass objective requires a positive invertible kg scale"));
        }
        let masses = study.operator.reference_cell_masses_kg()
            .ok_or(ProjectedAlError::Invalid("mass objective requires native source material densities"))?;
        if masses.len() != study.cells() || masses.iter().any(|m| !m.is_finite() || *m <= 0.0) {
            return Err(ProjectedAlError::Invalid("mass objective has invalid reference cell inventory"));
        }
    }
    Ok(())
}
fn optimization_sample(e: &ResponseLimitEvaluation3, options: MinimumMaterialOptions3)
    -> Result<ProjectedAlSample, ResponseError3> {
    let objective = objective_value(e.volume_fraction, e.mass_kg, options)
        .ok_or(ResponseError3::Invalid("missing or unrepresentable minimum-material objective"))?;
    let (gradient, scale) = match options.objective {
        MinimumMaterialObjective3::VolumeFraction => (&e.volume_gradient, options.volume_scale),
        MinimumMaterialObjective3::MassKg { scale_kg } =>
            (e.mass_gradient.as_ref().ok_or(ResponseError3::Invalid("missing material mass gradient"))?, scale_kg),
    };
    Ok(ProjectedAlSample { objective,
        gradient: gradient.iter().map(|g| {
            let normalized = finite(g/scale)?;
            if *g != 0.0 && normalized == 0.0 {
                return Err(ResponseError3::Invalid("minimum-material objective gradient underflow"));
            }
            Ok(normalized)
        }).collect::<Result<_, _>>()?,
        constraint: e.constraint, constraint_gradient: e.constraint_gradient.clone() })
}
fn evaluate<O: AdaptiveSdf3Elasticity>(study: &mut CutDensityStudy3<O>, cases: &[ResponseLimitCase3<'_>],
    rho: &[f64], options: MinimumMaterialOptions3, control: &mut SolveControl<'_>)
    -> Result<ResponseLimitEvaluation3, ResponseError3> {
    study.response_limits_with_mass(rho, cases, options.limits,
        matches!(options.objective, MinimumMaterialObjective3::MassKg { .. }), control)
}
fn sample<O: AdaptiveSdf3Elasticity>(study: &mut CutDensityStudy3<O>, cases: &[ResponseLimitCase3<'_>],
    rho: &[f64], options: MinimumMaterialOptions3, last: &mut Option<ResponseLimitEvaluation3>,
    control: &mut SolveControl<'_>) -> Result<Option<ProjectedAlSample>, ResponseError3> {
    let e = evaluate(study, cases, rho, options, control)?;
    let s = optimization_sample(&e, options)?;
    *last = Some(e); Ok(Some(s))
}

/// Resumable fixed-grid design minimizing projected material volume or mass while
/// bounding independently solved displacement and embedded reaction responses.
/// The one smooth inequality is a conservative aggregate of all original limits;
/// individual feasibility is checked again when retaining the best endpoint.
///
/// Ordinary AL accepted steps may be infeasible. `best_feasible()` stores the
/// least-objective ACCEPTED feasible point separately from the last optimizer
/// iterate. A stop/KKT result belongs to the latter; never to an earlier best.
/// Only fully accepted steps install physical scales. Late cancellation after
/// an optimizer transition still synchronizes its complete fields. Probes and
/// rejected trials retain work but never replace accepted material.
///
/// Smoothing stays fixed across resumed runs. A tight window can be infeasible
/// for the upper aggregate even when the original inequalities are feasible.
/// This is not an exact equality-constrained solver or continuum certificate.
/// Mass uses the retained heterogeneous native material densities and the same
/// projected design as stiffness. The geometric filter and reported volume stay
/// unchanged. Reference material properties remain fixed during a session;
/// material-catalog switching and geometric shape changes are not modeled.
pub struct MinimumMaterialStudy3<'a, 'callback, O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    cases: &'a [ResponseLimitCase3<'a>],
    control: &'a mut SolveControl<'callback>,
    options: MinimumMaterialOptions3,
    state: ProjectedAlState,
    accepted: ResponseLimitEvaluation3,
    best_feasible: Option<ResponseLimitEvaluation3>,
    history: Vec<MinimumMaterialIteration3>,
    restoration_evaluations: usize,
}
impl<'a, 'callback, O: AdaptiveSdf3Elasticity> MinimumMaterialStudy3<'a, 'callback, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseLimitCase3<'a>], rho: &[f64],
        options: MinimumMaterialOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<ResponseError3>> {
        let n = study.cells(); options.optimizer.validate(n)?;
        if rho.len() != n || !options.density_floor.is_finite() || options.density_floor <= 0.0
            || options.density_floor >= 1.0 || !options.volume_scale.is_finite() || options.volume_scale <= 0.0
            || !(1.0/options.volume_scale).is_finite()
            || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("invalid minimum-material options or initial design"));
        }
        admit_objective(study, options)?;
        admit_limits(study, cases, options.limits).map_err(ProjectedAlError::Evaluation)?;
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |x| sample(study, cases, x, options, &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("response-limit-design-control").is_ok() {
                    ControlFlow::Continue(()) } else { ControlFlow::Break(()) })?
        };
        control.checkpoint("response-limit-design-initialize")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        let accepted = last.expect("accepted initial sample owns complete physical fields");
        study.operator.set_scales(&accepted.scales).expect("evaluated scales are admitted");
        let best_feasible = accepted.feasible(options.optimizer.tolerance).then(|| accepted.clone());
        let history = vec![row(&accepted, 0, state.sample().objective)];
        Ok(Self { study, cases, control, options, state, accepted, best_feasible, history, restoration_evaluations: 0 })
    }
    #[must_use] pub fn accepted(&self) -> &ResponseLimitEvaluation3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&ResponseLimitEvaluation3> { self.best_feasible.as_ref() }
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use] pub fn history(&self) -> &[MinimumMaterialIteration3] { &self.history }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }

    /// Continue exactly the same multiplier, penalty and spectral-step state.
    /// Zero updates use retained derivatives. Outer limits and all physical
    /// solve/preparation work remain cumulative, including unsuccessful probes.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<ResponseError3>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations; let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study; let cases = self.cases; let options = self.options;
                self.state.try_run(remaining.min(1),
                    &mut |x| sample(study, cases, x, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("response-limit-design-control").is_ok() {
                        ControlFlow::Continue(()) } else { ControlFlow::Break(()) })
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted material step owns complete physical fields");
                assert_eq!(accepted.rho, self.state.point(), "accepted limit fields match the optimizer point");
                self.study.operator.set_scales(&accepted.scales).expect("accepted scales are admitted");
                self.accepted = accepted;
                if self.accepted.feasible(self.options.optimizer.tolerance)
                    && self.best_feasible.as_ref().is_none_or(|best|
                        objective_inventory(self.accepted.volume_fraction, self.accepted.mass_kg, self.options.objective)
                        < objective_inventory(best.volume_fraction, best.mass_kg, self.options.objective)) {
                    self.best_feasible = Some(self.accepted.clone());
                }
                self.history.push(row(&self.accepted, self.state.work().iterations, self.state.sample().objective));
                self.control.checkpoint("response-limit-design-accepted")
                    .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
    /// End the optimizer and install the least-objective accepted feasible design.
    /// Consuming the session prevents further iterations from confusing this
    /// earlier endpoint with the latest AL state. No extra solve is performed.
    /// A missing incumbent is an error, never a feasible-result placeholder.
    /// Returned work includes all spent trials, not only that endpoint's solves.
    pub fn into_best_feasible(self) -> Result<ResponseLimitEvaluation3, ResponseError3> {
        self.control.checkpoint("response-limit-design-select")?;
        let mut best = self.best_feasible.ok_or(ResponseError3::Invalid("no accepted response-limit-feasible design"))?;
        self.study.operator.set_scales(&best.scales)?;
        best.work = self.control.work();
        Ok(best)
    }
}
