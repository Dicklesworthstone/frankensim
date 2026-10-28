//! Full two-axis modal recoil, retained shear and original finite-body heat.
//!
//! Every thermal trial repeats the same matrix-coupled mechanical proposal
//! from accepted history. Its irreversible contact work enters the original
//! caloric cells once; both components of reversible shear remain mechanical.

use super::{
    PlanarTangentialModalCandidate, PlanarTangentialModalDrive, PlanarTangentialModalError,
    PlanarTangentialModalFrame, PlanarTangentialModalSystem, TangentialLaw,
};
use crate::friction_heat::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
};
use fs_exec::CancelGate;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};

pub use super::super::heat::{ModalFrictionHeatSolve, TangentialHeatSource};

/// A refused source, full matrix contact, thermal or nonlinear proposal.
#[derive(Debug)]
pub enum PlanarTangentialHeatError<E> {
    Invalid(&'static str),
    Mechanical(PlanarTangentialModalError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}

impl<E> From<PlanarTangentialModalError> for PlanarTangentialHeatError<E> {
    fn from(error: PlanarTangentialModalError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for PlanarTangentialHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for PlanarTangentialHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for PlanarTangentialHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar modal friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for PlanarTangentialHeatError<E> {}

/// One accepted two-axis mechanical and finite thermal interval.
pub struct HeatedPlanarTangentialModalFrame<P> {
    pub mechanical: PlanarTangentialModalFrame,
    pub heat: FrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Mechanical plus caloric storage/loss minus outside supplies. Physical
    /// friction heat cancels as an internal transfer exactly once.
    pub combined_energy_residual_j: f64,
}

struct Staged<S, P> {
    mechanical: PlanarTangentialModalCandidate<S>,
    source: P,
}

fn numerical<E>(
    error: FiniteDifferenceError<PlanarTangentialHeatError<E>>,
) -> PlanarTangentialHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => PlanarTangentialHeatError::Invalid(message),
    }
}

impl<L: TangentialLaw> PlanarTangentialModalSystem<L> {
    /// Solve both contact directions, body motion and finite-body heat together.
    ///
    /// The route must represent this contact's actual bodies. Its geometry and
    /// the original endpoint conductivities determine their signed heat flows,
    /// including ordinary conduction when friction loss is zero. The callback
    /// may return source receipts, but must not publish external state during
    /// probes or supply the mechanical contact loss a second time.
    ///
    /// `drive.contact_temperature_k` must be `None`: the common temperature is
    /// solved here. The temperature guess is only a numerical starting point;
    /// normal load, pressure and both driver velocities remain explicit. The
    /// original full 2-by-2 modal response is retained in every thermal trial.
    ///
    /// The same law and accepted history own all constitutive evaluations.
    /// Temperature-dependent dissipation is supported when that law admits it.
    /// A frozen-temperature elastic source must refuse a changed temperature;
    /// this method never rebuilds its parameters or resets stored shear to make
    /// a trial admissible. Laws with variable elastic free energy require an
    /// explicit reversible heat/entropy coupling beyond this dissipative port.
    /// They must refuse unsupported parameter changes in their own `advance`.
    ///
    /// Unknown enthalpy increments follow sorted CellId order. Original
    /// backward-Euler conduction and held-effort mechanics require refinement;
    /// tangent-frame transport and normal evolution are outside this owner.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        route: &FrictionHeatRoute,
        external: &[f64],
        drive: PlanarTangentialModalDrive,
        controls: ModalFrictionHeatSolve,
        energy_guess_j: &[f64],
        temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<TangentialHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedPlanarTangentialModalFrame<P>, PlanarTangentialHeatError<E>> {
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let sample = self.network.samples_rendered();
        let next_sample = sample
            .checked_add(1)
            .ok_or(PlanarTangentialHeatError::Invalid(
                "modal sample clock overflow",
            ))?;
        let dt = self.network.sample_period_s();
        let start = sample as f64 * dt;
        let target = next_sample as f64 * dt;
        if thermal.time_s() != start
            || !target.is_finite()
            || target <= start
            || drive.contact_temperature_k.is_some()
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !controls.temperature_scale_j_per_k.is_finite()
            || controls.temperature_scale_j_per_k <= 0.0
            || !controls.combined_energy_tolerance_j.is_finite()
            || controls.combined_energy_tolerance_j < 0.0
            || !temperature_guess_k.is_finite()
            || temperature_guess_k <= 0.0
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
        {
            return Err(PlanarTangentialHeatError::Invalid(
                "synchronized clocks, solved temperature and bounded heat coordinates required",
            ));
        }
        let step = FrictionHeatStep::new(thermal, target, route, checkpoint)?;
        let cells = step.cell_ids().len();
        let trials = Cell::new(0usize);
        let owner = RefCell::new(&mut *self);
        let evaluate = |x: &[f64]| {
            checkpoint()?;
            if trials.get() >= controls.maximum_trials {
                return Err(PlanarTangentialHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != cells + 1 {
                return Err(PlanarTangentialHeatError::Invalid(
                    "incomplete planar tangential friction/heat iterate",
                ));
            }
            let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
            step.evaluate(
                &x[..cells],
                temperature,
                |endpoint, interface_temperature| {
                    let TangentialHeatSource {
                        conductivity_w_m_k,
                        transport,
                    } = resolve(endpoint, interface_temperature)
                        .map_err(PlanarTangentialHeatError::Source)?;
                    let StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    let mechanical = owner.borrow_mut().stage_candidate(
                        external,
                        PlanarTangentialModalDrive {
                            contact_temperature_k: Some(interface_temperature),
                            ..drive
                        },
                        Some(gate),
                    )?;
                    Ok(FrictionHeatInput {
                        dissipated_j: mechanical.frame.friction_dissipation_j,
                        conductivity_w_m_k,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: Staged {
                                mechanical,
                                source: payload,
                            },
                        },
                    })
                },
                checkpoint,
            )
        };
        let problem =
            FallibleFiniteDifference::new(cells + 1, controls.difference_step_j, |x: &[f64]| {
                evaluate(x).map(|trial| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.push(0.0);
        let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                gate.is_requested() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(PlanarTangentialHeatError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(PlanarTangentialHeatError::NotConverged(result.report));
        }
        // A finite-difference probe may leave scratch at a different root.
        // Restage both efforts and the complete contact state at accepted x.
        let trial = evaluate(&solver.x)?;
        let combined = trial.payload().mechanical.frame.energy_residual_j
            + trial.report().total_energy_residual_j;
        if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
            return Err(PlanarTangentialHeatError::Invalid(
                "joint planar mechanical/thermal work does not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        let mechanical = staged.mechanical.frame;
        // All fallible work finishes before this joint publication.
        owner.publish_candidate(staged.mechanical);
        *thermal = next_thermal;
        Ok(HeatedPlanarTangentialModalFrame {
            mechanical,
            heat,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }
}
