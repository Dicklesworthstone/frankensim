//! Simultaneous material/thermal feedback for the existing six-axis bushings.
//!
//! Trial endpoint temperatures resolve the original interface claims. Native
//! implicit mechanics then supplies actual midpoint dissipation to the original
//! backward-Euler cell-energy residual. No trial advances accepted history.
//! Endpoint-frozen coefficients and backward Euler remain first-order in time;
//! nonlinear convergence is not a trajectory-accuracy certificate.

use super::*;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use std::cell::{Cell, RefCell};

/// Physical, nonlinear and work limits for one event-free coupled interval.
#[derive(Clone, Copy, Debug)]
pub struct CoupledBushingHeatControls {
    /// Original simultaneous mechanical solve and energy/parameter-work gates.
    pub mechanics: BushingHeatControls,
    /// Explicit accuracy limit for endpoint-frozen material coefficients [s].
    pub maximum_step_s: f64,
    pub maximum_cells: usize,
    /// Includes every source/mechanical evaluation and numerical stencil.
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Maximum-component directional-difference perturbation in cell joules.
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub thermal_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

/// A concrete source/physics failure is retained before any solver outcome.
#[derive(Debug)]
pub enum CoupledBushingHeatError {
    Physical(BushingHeatError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    Invalid(&'static str),
    Budget(&'static str),
}
impl From<BushingHeatError> for CoupledBushingHeatError {
    fn from(e: BushingHeatError) -> Self {
        Self::Physical(e)
    }
}
impl From<DepositionError> for CoupledBushingHeatError {
    fn from(e: DepositionError) -> Self {
        Self::Physical(e.into())
    }
}
impl From<NewtonError> for CoupledBushingHeatError {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl core::fmt::Display for CoupledBushingHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled heated bushings: {self:?}")
    }
}
impl std::error::Error for CoupledBushingHeatError {}

/// The original physical accounts and the actual converged residual evidence.
pub struct CoupledBushingHeatReport {
    pub physical: BushingHeatReport,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    pub trials: usize,
}

impl<'a> HeatedBushings<'a> {
    /// Close damping heat and source-temperature feedback in the SAME interval.
    ///
    /// Unknowns are cell-energy increments [J] in the original sorted CellId
    /// order. Every trial starts from this accepted state, refreshes all original
    /// material claims at proposed endpoint temperatures, and advances all
    /// bushings in one existing native implicit mechanical solve. Its midpoint
    /// loss supplies each original heat destination exactly once.
    ///
    /// Source stiffness changes are applied at the initial mechanical strain;
    /// U_new(q_initial)-U_old(q_initial) is explicit external parameter work.
    /// It is never deposited as heat. This convention does not supply reciprocal
    /// thermoelastic free energy or entropy. Carrier geometry/inertia stay fixed.
    /// Refine time and native mechanical tolerances independently of this solve.
    ///
    /// Start/interior events refuse; target events remain pending. The returned
    /// state can continue with this method or the existing lagged scheduler.
    /// Arbitrary additional observers are not inferred into this closed block.
    /// Cancellation, source refusal, any exceeded budget, unconverged iteration
    /// or failed physical balance returns neither proposed history.
    pub fn step_to_coupled(
        &self,
        target_s: f64,
        loads: SpatialMotionLoads<'_>,
        controls: CoupledBushingHeatControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, CoupledBushingHeatReport), CoupledBushingHeatError> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let count = self.thermal.cells().len();
        let dt = target_s - self.thermal.time_s();
        if !loads.connections.is_empty()
            || !controls.maximum_step_s.is_finite()
            || controls.maximum_step_s <= 0.0
            || dt > controls.maximum_step_s
            || controls.maximum_cells == 0
            || count > controls.maximum_cells
            || initial_energy_increments_j.len() != count
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || [
                controls.mechanics.maximum_parameter_work_j,
                controls.mechanics.energy_tolerance_j,
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
            ]
            .into_iter()
            .any(|x| !nonnegative(x))
        {
            return Err(CoupledBushingHeatError::Invalid(
                "unrouted loads or invalid coupled step/energy/work budgets",
            ));
        }
        let thermal_step = BackwardEulerStep::new(&self.thermal, target_s, checkpoint)?;
        let trials = Cell::new(0usize);
        let evaluate = |increments: &[f64]| {
            checkpoint()?;
            if trials.get() >= controls.maximum_trials {
                return Err(CoupledBushingHeatError::Budget("source/mechanical trials"));
            }
            trials.set(trials.get() + 1);
            let increments = thermal_step.resolvable_increments_j(increments);
            thermal_step.evaluate(
                &increments,
                |thermal| {
                    self.coupled_coefficients(thermal, loads, dt, controls.mechanics, &mut || {
                        checkpoint()
                    })
                    .map_err(CoupledBushingHeatError::from)
                },
                checkpoint,
            )
        };
        let problem =
            FallibleFiniteDifference::new(count, controls.difference_step_j, |x: &[f64]| {
                Ok::<_, CoupledBushingHeatError>(evaluate(x)?.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(
            &problem,
            initial_energy_increments_j.to_vec(),
            controls.newton,
        );
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let outcome =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if outcome.progress != SolverRunProgress::Complete {
            return Err(CoupledBushingHeatError::Interrupted(outcome.progress));
        }
        if !outcome.report.converged {
            return Err(CoupledBushingHeatError::NotConverged(outcome.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (thermal, payload, audit) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        // Ask the original thermal owner for external supply independently of
        // damping heat. The trial source is not installed as a persistent load.
        let external = finite(
            dt * thermal.external_step_power_w(&BTreeMap::new(), &BTreeMap::new(), checkpoint)?,
        )?;
        let native = &payload.step.spatial;
        let supplied = finite(native.motion.external_work + native.motion.actuator_work)?;
        let mechanical_change = finite(
            native.motion.momentum_after.kinetic_energy
                - native.motion.momentum_before.kinetic_energy
                + native.stored_energy_after
                - native.stored_energy_before
                + payload.parameter,
        )?;
        let residual = finite(
            mechanical_change + audit.stored_energy_change_j
                - supplied
                - payload.parameter
                - external,
        )?;
        if (audit.external_heat_j - external - payload.heat).abs()
            > controls.mechanics.energy_tolerance_j
            || residual.abs() > controls.mechanics.energy_tolerance_j
        {
            return Err(CoupledBushingHeatError::Invalid(
                "coupled native loss, thermal supply or total energy balance failed",
            ));
        }
        checkpoint()?;
        let physical = BushingHeatReport {
            thermal: audit,
            dissipated_heat_j: payload.heat,
            external_heat_j: external,
            elastic_parameter_work_j: payload.parameter,
            absolute_parameter_work_j: payload.absolute_parameter,
            mechanical_supplied_work_j: supplied,
            mechanical_energy_change_j: mechanical_change,
            mechanical_work_residual_j: native.mechanical_work_residual,
            total_energy_residual_j: residual,
        };
        let next = Self {
            reference: self.reference,
            thermal,
            regions: self.regions.clone(),
            bindings: payload.bindings,
            state: payload.step.spatial.motion.state_after,
        };
        Ok((
            next,
            CoupledBushingHeatReport {
                physical,
                nonlinear: outcome.report,
                maximum_cell_residual_j,
                trials: trials.get(),
            },
        ))
    }

    fn coupled_coefficients(
        &self,
        thermal: &DepositionSimulation<'a>,
        loads: SpatialMotionLoads<'_>,
        dt: f64,
        controls: BushingHeatControls,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<StepCoefficients<CoupledTrial<'a>>, BushingHeatError> {
        let mut bindings = self.bindings.clone();
        let mut parameter = 0.0;
        let mut absolute_parameter = 0.0;
        for binding in &mut bindings {
            checkpoint()?;
            let point = actual_point(binding, &self.regions, thermal)?;
            let rebound = checked(checkpoint, |cancel| {
                binding.bushing.rebind_at_state(
                    binding.point.interface,
                    &point,
                    &self.state,
                    controls.maximum_parameter_work_j - absolute_parameter,
                    cancel,
                )
            })?;
            parameter = finite(parameter + rebound.parameter_work_j)?;
            absolute_parameter = finite(absolute_parameter + rebound.parameter_work_j.abs())?;
            if absolute_parameter > controls.maximum_parameter_work_j {
                return Err(BushingHeatError::Invalid(
                    "elastic parameter work exceeds coupled step budget",
                ));
            }
            binding.bushing = rebound.binding;
        }
        let refs = bindings.iter().map(|b| &b.bushing).collect::<Vec<_>>();
        let step = checked(checkpoint, |cancel| {
            step_bushings_implicit(&refs, &self.state, loads, dt, controls.mechanics, cancel)
        })?;
        if step.spatial.mechanical_work_residual.abs() > controls.energy_tolerance_j {
            return Err(BushingHeatError::Invalid(
                "native mechanical work defect exceeds budget",
            ));
        }
        let mut sources = BTreeMap::new();
        let mut heat = 0.0;
        for binding in &bindings {
            let power = checked(checkpoint, |cancel| {
                binding.bushing.evaluate(&step.midpoint_state, cancel)
            })?
            .dissipation_power;
            let loss = finite(dt * power)?;
            if loss < 0.0 {
                return Err(BushingHeatError::Invalid("negative native damping heat"));
            }
            let first = finite(loss * binding.fraction_to_a)?;
            for (cell, energy) in binding.cells.into_iter().zip([first, loss - first]) {
                let watts = finite(energy / dt)?;
                let value = sources.entry(cell).or_insert(0.0);
                *value = finite(*value + watts)?;
            }
            heat = finite(heat + loss)?;
        }
        let routed = finite(dt * sources.values().sum::<f64>())?;
        if (heat - step.spatial.dissipated_energy).abs() > controls.energy_tolerance_j
            || (routed - heat).abs() > controls.energy_tolerance_j
        {
            return Err(BushingHeatError::Invalid(
                "native loss and explicit heat routes differ",
            ));
        }
        checkpoint()?;
        Ok(StepCoefficients {
            sources_w: sources,
            conductances_w_k: BTreeMap::new(),
            payload: CoupledTrial {
                bindings,
                step,
                heat,
                parameter,
                absolute_parameter,
            },
        })
    }
}

struct CoupledTrial<'a> {
    bindings: Vec<BushingHeatBinding<'a>>,
    step: ImplicitSpatialStep,
    heat: f64,
    parameter: f64,
    absolute_parameter: f64,
}

fn numerical(error: FiniteDifferenceError<CoupledBushingHeatError>) -> CoupledBushingHeatError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => CoupledBushingHeatError::Invalid(reason),
    }
}
