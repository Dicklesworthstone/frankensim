//! Retained tangential memory and finite-body conduction in one transaction.
//!
//! The immutable law sees the solved common interface temperature on every
//! trial. Only irreversible contact loss enters the original caloric cells;
//! reversible shear storage stays in the mechanical account. No heat fraction,
//! extra thermal store, replacement law or implicit history transfer is added.

use super::{
    TangentialLaw, TangentialModalCandidate, TangentialModalDrive, TangentialModalFrame,
    TangentialModalSystem,
};
use crate::friction_heat::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
};
use fs_exec::CancelGate;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};

pub use super::super::heat::{ModalFrictionHeatError, ModalFrictionHeatSolve};

/// Original material transport at a candidate endpoint; no constitutive reset.
pub struct TangentialHeatSource<P> {
    /// Conductivities of the two original bodies at their own cell temperatures.
    pub conductivity_w_m_k: [f64; 2],
    /// Other heat/conduction and source receipts. Do not reinsert friction heat.
    pub transport: StepCoefficients<P>,
}

/// One accepted contact/body/conduction interval on the original sample clock.
pub struct HeatedTangentialModalFrame<P> {
    pub mechanical: TangentialModalFrame,
    pub heat: FrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Complete mechanical plus caloric storage/loss minus external supplies.
    /// Friction heat cancels as an internal transfer, exactly once.
    pub combined_energy_residual_j: f64,
}

struct Staged<S, P> {
    mechanical: TangentialModalCandidate<S>,
    source: P,
}

fn numerical<E>(
    error: FiniteDifferenceError<ModalFrictionHeatError<E>>,
) -> ModalFrictionHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => ModalFrictionHeatError::Invalid(message),
    }
}

impl<L: TangentialLaw> TangentialModalSystem<L> {
    /// Solve retained shear, body motion and original finite-body heat together.
    ///
    /// The route must identify this mechanical port's actual contacting bodies.
    /// Its geometry and the endpoint conductivities determine the heat split.
    /// The callback may return original source receipts but must not publish
    /// external state during a trial or supply the contact loss a second time.
    ///
    /// `drive.contact_temperature_k` must be `None`: the common temperature is
    /// solved here, not prescribed. `temperature_guess_k` is only a numerical
    /// initial guess. Pressure, normal load and drive speed remain explicit.
    ///
    /// The law and accepted history are unchanged throughout Newton probes.
    /// Temperature-dependent dissipation is supported when the law admits it.
    /// A frozen-temperature source law can refuse a changed temperature; this
    /// operation never rebuilds its source or erases retained shear to proceed.
    /// Temperature-dependent elastic coefficients need an additional physical
    /// free-energy/entropy coupling and are not manufactured by this adapter.
    /// Backward-Euler conduction and held-effort mechanics require refinement.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        route: &FrictionHeatRoute,
        external: &[f64],
        drive: TangentialModalDrive,
        controls: ModalFrictionHeatSolve,
        energy_guess_j: &[f64],
        temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<TangentialHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedTangentialModalFrame<P>, ModalFrictionHeatError<E>> {
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
            .ok_or(ModalFrictionHeatError::Invalid(
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
            return Err(ModalFrictionHeatError::Invalid(
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
                return Err(ModalFrictionHeatError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != cells + 1 {
                return Err(ModalFrictionHeatError::Invalid(
                    "incomplete tangential friction/heat iterate",
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
                        .map_err(ModalFrictionHeatError::Source)?;
                    let StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    let mechanical = owner.borrow_mut().stage_candidate(
                        external,
                        TangentialModalDrive {
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
            return Err(ModalFrictionHeatError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(ModalFrictionHeatError::NotConverged(result.report));
        }
        // Restage once at the actual accepted root: numerical probes may have
        // left mechanical scratch at another temperature, but never published it.
        let trial = evaluate(&solver.x)?;
        let combined = trial.payload().mechanical.frame.energy_residual_j
            + trial.report().total_energy_residual_j;
        if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
            return Err(ModalFrictionHeatError::Invalid(
                "joint tangential mechanical/thermal work does not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        owner.network.validate_staged_material_points(Some(gate))
            .map_err(super::super::ModalFrictionError::from)?;
        let mechanical = staged.mechanical.frame;
        // All fallible work and callbacks precede this joint publication.
        owner.publish_candidate(staged.mechanical);
        *thermal = next_thermal;
        Ok(HeatedTangentialModalFrame {
            mechanical,
            heat,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }
}
