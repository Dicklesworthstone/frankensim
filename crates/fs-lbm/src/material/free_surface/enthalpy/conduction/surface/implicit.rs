//! Opt-in stiff surface exchange through the original thermal residual owner.
//!
//! Geometry and actual masses remain fixed during each thermal/flow split.
//! Newton unknowns are specific-enthalpy increments, so a numerical Jacobian
//! probe does not impose the same joules on cells with very different masses.
//! A frozen positive row scaling conditions the solve; acceptance always uses
//! the original unscaled per-cell and total energy residuals.

use super::*;
use fs_conduction::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::Cell;

/// Numerical work and physical acceptance controls for prescribed reservoirs.
/// This does not change interface geometry, mass, heat capacity or source law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceImplicitSolve {
    pub maximum_unknowns: usize,
    pub maximum_krylov_restart: usize,
    pub maximum_newton_iterations: usize,
    /// Shared by all intervals, Jacobian probes and final residual evaluations.
    pub maximum_residual_evaluations: usize,
    /// Maximum-component numerical Jacobian perturbation in J/kg, not joules.
    pub difference_step_j_kg: f64,
    /// Independent maximum actual cell energy-balance defect in joules.
    /// The total gate remains `FreeSurfaceConductionControls::energy_tolerance_j`.
    pub cell_energy_tolerance_j: f64,
    /// Residual coordinates are J/kg after frozen positive row normalization.
    pub newton: NewtonKrylovConfig,
}

impl Default for FreeSurfaceImplicitSolve {
    fn default() -> Self {
        Self {
            maximum_unknowns: 100_000,
            maximum_krylov_restart: 64,
            maximum_newton_iterations: 16,
            maximum_residual_evaluations: 4_096,
            difference_step_j_kg: 0.01,
            cell_energy_tolerance_j: 1e-10,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-6,
                // Material/interface roots and directional differences impose
                // a numerical tangent floor. Do not demand an inner Krylov
                // solve below it; outer Newton and raw joule gates still apply.
                forcing_minimum: 1e-5,
                ..NewtonKrylovConfig::default()
            },
        }
    }
}

/// Accepted implicit work. Numerical convergence alone is not an accuracy test.
#[derive(Clone, Debug, Default)]
pub struct FreeSurfaceImplicitReport {
    pub intervals: usize,
    pub residual_evaluations: usize,
    pub newton_iterations: usize,
    pub maximum_cell_residual_j: f64,
    pub last_nonlinear: Option<NewtonReport>,
}

/// Numerical refusal; original material and thermal errors retain their types.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceImplicitFailure {
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    EvaluationBudget,
}

impl core::fmt::Display for FreeSurfaceImplicitFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Newton(e) => write!(f, "surface thermal Newton: {e}"),
            Self::NotConverged(r) => write!(
                f,
                "surface thermal residual {} after {} attempts: {:?}",
                r.residual_norm, r.iterations, r.diagnosis
            ),
            Self::Interrupted(p) => write!(f, "surface thermal solve interrupted: {p:?}"),
            Self::EvaluationBudget => {
                write!(f, "surface thermal residual evaluation budget exhausted")
            }
        }
    }
}
impl std::error::Error for FreeSurfaceImplicitFailure {}

fn numerical(error: FiniteDifferenceError<Error>) -> Error {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(reason) => Error::Invalid(reason),
    }
}

fn validate(solve: FreeSurfaceImplicitSolve, count: usize) -> Result<(), Error> {
    solve
        .newton
        .validate()
        .map_err(|e| Error::ImplicitSurface(FreeSurfaceImplicitFailure::Newton(e)))?;
    if count == 0
        || count > solve.maximum_unknowns
        || solve.maximum_krylov_restart == 0
        || solve.newton.linear_restart > solve.maximum_krylov_restart
        || solve.maximum_newton_iterations == 0
        || solve.maximum_residual_evaluations == 0
        || !solve.difference_step_j_kg.is_finite()
        || solve.difference_step_j_kg <= 0.0
        || !solve.cell_energy_tolerance_j.is_finite()
        || solve.cell_energy_tolerance_j < 0.0
    {
        return Err(Error::Invalid(
            "invalid implicit surface dimension, work, difference or energy controls",
        ));
    }
    Ok(())
}

// This adapter changes only nonlinear coordinates and conditioning. Every
// residual, material-backed storage increment and final joule gate belongs to
// BackwardEulerStep, including in the independent tiny-capacity tests below.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn solve_interval<'t, T>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    scales_kg: &[f64],
    solve: FreeSurfaceImplicitSolve,
    total_tolerance_j: f64,
    evaluations: &Cell<usize>,
    resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<T>, Error>,
    cancelled: impl Fn() -> bool,
) -> Result<(DepositionSimulation<'t>, T, EnergyAudit, NewtonReport, f64), Error> {
    let check = || {
        if cancelled() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    validate(solve, thermal.cells().len())?;
    let step = BackwardEulerStep::new(thermal, target_s, check)?;
    if scales_kg.len() != step.cell_ids().len()
        || scales_kg.iter().any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(Error::Invalid(
            "implicit thermal row scales must be complete and positive",
        ));
    }
    let masses: Vec<_> = step
        .cell_ids()
        .iter()
        .map(|id| thermal.cells()[id].mass_kg)
        .collect();
    let evaluate = |increments: &[f64]| {
        check()?;
        if evaluations.get() >= solve.maximum_residual_evaluations {
            return Err(Error::ImplicitSurface(
                FreeSurfaceImplicitFailure::EvaluationBudget,
            ));
        }
        evaluations.set(evaluations.get() + 1);
        let energy = increments
            .iter()
            .zip(&masses)
            .map(|(&h, &mass)| {
                let joules = environment::finite(h * mass)?;
                if h != 0.0 && joules == 0.0 {
                    return Err(Error::Invalid(
                        "implicit specific-energy increment underflows actual cell energy",
                    ));
                }
                Ok(joules)
            })
            .collect::<Result<Vec<_>, Error>>()?;
        step.evaluate(&energy, &resolve, check)
    };
    let problem = FallibleFiniteDifference::new(
        step.cell_ids().len(),
        solve.difference_step_j_kg,
        |x: &[f64]| {
            evaluate(x).map(|trial| {
                trial
                    .residual_j()
                    .iter()
                    .zip(scales_kg)
                    .map(|(&r, &s)| r / s)
                    .collect()
            })
        },
    )
    .map_err(numerical)?;
    let construction =
        NewtonKrylovState::new(&problem, vec![0.0; step.cell_ids().len()], solve.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver =
        construction.map_err(|e| Error::ImplicitSurface(FreeSurfaceImplicitFailure::Newton(e)))?;
    let outcome = solver.run_controlled(&problem, solve.maximum_newton_iterations, &mut || {
        cancelled() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    check()?;
    if outcome.progress != SolverRunProgress::Complete {
        return Err(Error::ImplicitSurface(
            FreeSurfaceImplicitFailure::Interrupted(outcome.progress),
        ));
    }
    if !outcome.report.converged {
        return Err(Error::ImplicitSurface(
            FreeSurfaceImplicitFailure::NotConverged(outcome.report),
        ));
    }
    let trial = evaluate(&solver.x)?;
    let maximum = trial.maximum_cell_residual_j();
    let (next, payload, audit) =
        trial.accept(solve.cell_energy_tolerance_j, total_tolerance_j, check)?;
    Ok((next, payload, audit, outcome.report, maximum))
}

// Compensated sums preserve small exchange increments over many intervals.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let next = environment::finite(self.value + value)?;
        self.correction = environment::finite(
            self.correction
                + if self.value.abs() >= value.abs() {
                    (self.value - next) + value
                } else {
                    (value - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn finish(&self) -> Result<f64, Error> {
        environment::finite(self.value + self.correction)
    }
}

/// Original timeline and observer, advanced at implicit endpoint quadrature.
/// Callers have already processed the initial birth/contact events. Other
/// histories require a separate coupled residual and therefore refuse here.
pub(in super::super) fn advance<'t, 'a, 'm>(
    mut thermal: DepositionSimulation<'t>,
    mut observer: Conduction<'a, 'm>,
    target_s: f64,
    controls: FreeSurfaceConductionControls,
    solve: FreeSurfaceImplicitSolve,
) -> Result<(DepositionSimulation<'t>, Conduction<'a, 'm>, EnergyAudit), Error> {
    validate(solve, thermal.cells().len())?;
    if observer.surfaces.is_empty()
        || observer.finite_atmosphere.is_some()
        || observer.electrical.is_some()
        || !observer.bodies.is_empty()
        || !observer.contacts.is_empty()
        || !observer.solid_links.is_empty()
        || !observer.environment.is_empty()
        || !observer.radiation.is_empty()
    {
        return Err(Error::Invalid(
            "implicit surface exchange requires prescribed surface reservoirs without finite gas, bodies, current or separate environment/radiation sources",
        ));
    }
    if !target_s.is_finite()
        || target_s <= thermal.time_s()
        || !controls.maximum_dt_s.is_finite()
        || controls.maximum_dt_s <= 0.0
        || !controls.energy_tolerance_j.is_finite()
        || controls.energy_tolerance_j < 0.0
    {
        return Err(Error::Invalid(
            "implicit surface exchange needs a later finite target and valid thermal controls",
        ));
    }
    let gate = observer.gate;
    let initial_h: Vec<_> = thermal
        .cells()
        .values()
        .map(|c| (c.id, c.specific_h_j_kg))
        .collect();
    let evaluations = Cell::new(0);
    let mut external = Sum::default();
    let mut report = FreeSurfaceImplicitReport::default();
    while thermal.time_s() < target_s {
        checkpoint(gate)?;
        if report.intervals >= controls.maximum_substeps {
            return Err(DepositionError::Budget("implicit surface intervals").into());
        }
        let start = thermal.time_s();
        let dt = controls.maximum_dt_s.min(target_s - start);
        let end = if dt == target_s - start {
            target_s
        } else {
            environment::finite(start + dt)?
        };
        if end <= start {
            return Err(DepositionError::Arithmetic.into());
        }
        // Pin the original gas source claims before any numerical trial. The
        // same samples give conditioning, never substitute for endpoint rates.
        observer.sources(&thermal, &mut || checkpoint(gate))?;
        let conductances = observer.conductances(&thermal, &mut || checkpoint(gate))?;
        let slope = observer.caloric.maximum_temperature_slope();
        let mut scales: BTreeMap<_, _> = thermal
            .cells()
            .iter()
            .map(|(&id, c)| (id, c.mass_kg))
            .collect();
        for (&id, link) in thermal.links() {
            let g = conductances
                .get(&id)
                .copied()
                .unwrap_or(link.conductance_w_k);
            for cell in [link.a, link.b] {
                *scales.get_mut(&cell).expect("active link endpoint") += dt * g * slope;
            }
        }
        for (&id, sample) in observer
            .surface_endpoints
            .iter()
            .zip(&observer.surface_last)
        {
            *scales.get_mut(&id).expect("active surface endpoint") +=
                dt * sample.liquid_decay_conductance_w_k * slope;
        }
        let scales: Vec<_> = scales.into_values().collect();
        let resolve = |endpoint: &DepositionSimulation<'t>| {
            let mut payload = observer.clone();
            let sources_w = payload.sources(endpoint, &mut || checkpoint(gate))?;
            let conductances_w_k = payload.conductances(endpoint, &mut || checkpoint(gate))?;
            Ok(StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            })
        };
        let (next, mut accepted, audit, nonlinear, maximum) = solve_interval(
            &thermal,
            end,
            &scales,
            solve,
            controls.energy_tolerance_j,
            &evaluations,
            resolve,
            || gate.is_requested(),
        )?;
        // Only the accepted payload carries rates at this BE endpoint. The
        // original observer now integrates those rates ONCE and checks sources;
        // no nonlinear trial calls step or advances any accepted history.
        accepted.step(start, end - start, &BTreeMap::new(), &next, &mut || {
            checkpoint(gate)
        })?;
        external.add(audit.external_heat_j)?;
        report.intervals += 1;
        report.newton_iterations += nonlinear.iterations;
        report.maximum_cell_residual_j = report.maximum_cell_residual_j.max(maximum);
        report.last_nonlinear = Some(nonlinear);
        thermal = next;
        observer = accepted;
    }
    let mut stored = Sum::default();
    for (id, h) in initial_h {
        checkpoint(gate)?;
        let cell = &thermal.cells()[&id];
        stored.add(environment::finite(
            cell.mass_kg * (cell.specific_h_j_kg - h),
        )?)?;
    }
    let (stored_energy_change_j, external_heat_j) = (stored.finish()?, external.finish()?);
    let residual_j = environment::finite(stored_energy_change_j - external_heat_j)?;
    if residual_j.abs() > controls.energy_tolerance_j {
        return Err(DepositionError::EnergyBalance(residual_j).into());
    }
    report.residual_evaluations = evaluations.get();
    let audit = EnergyAudit {
        incoming_mass_kg: 0.0,
        incoming_energy_j: 0.0,
        external_heat_j,
        stored_energy_change_j,
        residual_j,
        substeps: report.intervals,
        events: 0,
    };
    observer.implicit_surface_report = Some(report);
    checkpoint(gate)?;
    Ok((thermal, observer, audit))
}

#[cfg(test)]
#[path = "implicit/tests.rs"]
mod tests;
