//! Prescribed-flow nonlinear heat transport on the original thermal residual.
//!
//! Limited donor reconstruction is nonlinear in cell enthalpy. The existing
//! native Newton/Krylov owner solves the original backward-Euler cell-energy
//! equations; native calorics, reciprocal links and directed fluxes retain
//! their owners. No second time integrator or replacement heat capacity enters.

use std::cell::Cell;

use fs_exec::Cx;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport};
use fs_thermal::enthalpy::backward_euler::{
    BackwardEulerStep, BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, EnergyAudit};

use super::distributed::MAX_DISTRIBUTED_FLUID_CELLS;
use super::{GraphEnthalpyTransport, TransientTransportError, thermal_checkpoint};

/// Numerical controls for one original implicit interval. A small algebraic
/// residual does not establish time/space accuracy or physical source support.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitTransportControls {
    pub solver: NewtonKrylovConfig,
    /// Includes native wall/solid cells as well as the bound fluid inventory.
    pub maximum_cells: usize,
    pub maximum_newton_iterations: usize,
    /// Counts all residual, numerical derivative and final acceptance queries.
    pub maximum_endpoint_evaluations: usize,
    pub energy_scale_j: f64,
    pub finite_difference_step: f64,
    pub cell_tolerance_j: f64,
    pub total_tolerance_j: f64,
}

#[derive(Debug)]
pub struct ImplicitTransportReport<T> {
    /// Source payload at the same accepted unrelaxed thermal endpoint.
    pub heat: T,
    pub thermal: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    pub endpoint_evaluations: usize,
    pub nonlinear: NewtonReport,
}

#[derive(Debug)]
pub enum ImplicitTransportError<E> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Transport(TransientTransportError),
    Solver(NewtonError),
    Resolve(E),
    WorkBudget { performed: usize, maximum: usize },
    NotConverged(NewtonReport),
}
impl<E> From<DepositionError> for ImplicitTransportError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<TransientTransportError> for ImplicitTransportError<E> {
    fn from(error: TransientTransportError) -> Self {
        Self::Transport(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ImplicitTransportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit prescribed-flow heat transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ImplicitTransportError<E> {}

impl GraphEnthalpyTransport {
    /// Advance fixed native signed flow with endpoint-resolved thermal sources.
    ///
    /// This entry accepts donor-cell and explicitly prepared limited pipe
    /// transport. `resolve` supplies ONLY nonadvective heat coefficients; the
    /// bound native mass graph supplies internal transport and external stream
    /// enthalpy once. Original fixed cell sources and boundary/link laws retain
    /// their normal meaning. Other active thermal cells may store wall heat.
    /// The callback must be pure, source-supported and admit its own physical
    /// applicability at each queried state. Flow itself is prescribed over
    /// this interval; use the coupled drivers when momentum must respond.
    ///
    /// Every residual starts from the same accepted thermal state. Original
    /// local/global energy gates decide acceptance. Source refusal, exhausted
    /// work, cancellation and nonlinear failure publish no thermal state or
    /// source payload. Interior scheduled events retain the native refusal;
    /// split there before calling. Time and spatial refinement remain required.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_implicit_with_heat<'m, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitTransportControls,
        energy_guess_j: &[f64],
        resolve: impl Fn(&Cx<'_>, &DepositionSimulation<'m>) -> Result<BoundaryStepCoefficients<T>, E>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitTransportReport<T>), ImplicitTransportError<E>>
    {
        thermal_checkpoint(cx)?;
        controls
            .solver
            .validate()
            .map_err(ImplicitTransportError::Solver)?;
        self.validate_storage(cx, thermal)?;
        if !(1..=MAX_DISTRIBUTED_FLUID_CELLS).contains(&controls.maximum_cells)
            || thermal.cells().is_empty()
            || thermal.cells().len() > controls.maximum_cells
            || !(1..=1024).contains(&controls.maximum_newton_iterations)
            || !(1..=1_000_000).contains(&controls.maximum_endpoint_evaluations)
            || !positive(controls.energy_scale_j)
            || !positive(controls.finite_difference_step)
            || !nonnegative(controls.cell_tolerance_j)
            || !nonnegative(controls.total_tolerance_j)
        {
            return Err(ImplicitTransportError::Invalid(
                "bounded original thermal inventory, scales and physical residual limits required",
            ));
        }
        // Admit the fixed mass graph and any prepared reconstruction before
        // entering the source callback, including initially stationary flow.
        self.endpoint_transport(cx, thermal)?;
        let step = BackwardEulerStep::new(thermal, target_s, || thermal_checkpoint(cx))?;
        let count = step.cell_ids().len();
        if energy_guess_j.len() != count || energy_guess_j.iter().any(|v| !v.is_finite()) {
            return Err(ImplicitTransportError::Invalid(
                "complete finite cell-energy guess required",
            ));
        }
        let evaluations = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            thermal_checkpoint(cx)?;
            if evaluations.get() == controls.maximum_endpoint_evaluations {
                return Err(ImplicitTransportError::WorkBudget {
                    performed: evaluations.get(),
                    maximum: controls.maximum_endpoint_evaluations,
                });
            }
            evaluations.set(evaluations.get() + 1);
            let energies: Vec<_> = x.iter().map(|v| v * controls.energy_scale_j).collect();
            if energies.iter().any(|value| !value.is_finite()) {
                return Err(ImplicitTransportError::Invalid(
                    "unrepresentable cell-energy coordinates",
                ));
            }
            let energies = step.resolvable_increments_j(&energies);
            step.evaluate_with_transport(
                &energies,
                |candidate| {
                    let advected = self.endpoint_transport(cx, candidate)?;
                    let BoundaryStepCoefficients { bulk, boundaries } =
                        resolve(cx, candidate).map_err(ImplicitTransportError::Resolve)?;
                    thermal_checkpoint(cx)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = bulk;
                    for (id, watts) in advected.heat.bulk.sources_w {
                        thermal_checkpoint(cx)?;
                        let value = sources_w.entry(id).or_insert(0.0);
                        *value += watts;
                        if !value.is_finite() {
                            return Err(ImplicitTransportError::Invalid(
                                "combined endpoint source is nonfinite",
                            ));
                        }
                    }
                    Ok(TransportStepCoefficients {
                        heat: BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w,
                                conductances_w_k,
                                payload,
                            },
                            boundaries,
                        },
                        transfers: advected.transfers,
                    })
                },
                || thermal_checkpoint(cx),
            )
        };
        let difference =
            FallibleFiniteDifference::new(count, controls.finite_difference_step, |x: &[f64]| {
                let trial = evaluate(x)?;
                Ok(trial
                    .residual_j()
                    .iter()
                    .map(|v| v / controls.energy_scale_j)
                    .collect())
            })
            .map_err(difference_error)?;
        let guess = energy_guess_j
            .iter()
            .map(|v| v / controls.energy_scale_j)
            .collect();
        let initialized = NewtonKrylovState::new(&difference, guess, controls.solver);
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialized.map_err(ImplicitTransportError::Solver)?;
        let run =
            solver.run_controlled(&difference, controls.maximum_newton_iterations, &mut || {
                difference.has_failed() || cx.checkpoint().is_err()
            });
        if let Some(error) = difference.take_failure() {
            return Err(difference_error(error));
        }
        thermal_checkpoint(cx)?;
        if !run.report.converged {
            return Err(ImplicitTransportError::NotConverged(run.report));
        }
        let final_trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = final_trial.maximum_cell_residual_j();
        let (next, heat, thermal) = final_trial.accept(
            controls.cell_tolerance_j,
            controls.total_tolerance_j,
            || thermal_checkpoint(cx),
        )?;
        thermal_checkpoint(cx)?;
        Ok((
            next,
            ImplicitTransportReport {
                heat,
                thermal,
                maximum_cell_residual_j,
                endpoint_evaluations: evaluations.get(),
                nonlinear: run.report,
            },
        ))
    }
}

fn difference_error<E>(
    error: FiniteDifferenceError<ImplicitTransportError<E>>,
) -> ImplicitTransportError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => ImplicitTransportError::Invalid(message),
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
