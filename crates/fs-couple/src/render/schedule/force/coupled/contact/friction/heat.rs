//! The original modal friction step and finite thermal network accept together.
//!
//! Every thermal Newton trial resolves the law at the solved common interface
//! temperature and repeats the ORIGINAL mechanical proposal from accepted
//! history. Only physical friction loss enters the existing caloric cells.
//! Their heat split is a result of conduction, never a prescribed fraction.
//! Nonlinear mechanical storage, source applicability and branch history retain
//! their original owners. No state publishes until every owner accepts.

use super::*;
use crate::friction_heat::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};

/// Same-step modal friction with an existing finite thermal layer.
pub mod third_body;

/// Common-slip Coulomb composition with local moving-contact temperatures.
pub mod moving;

/// Endpoint constitutive and transport data, all evaluated at the same trial.
pub struct ModalHeatSource<P> {
    /// Original source's dry-friction law at the common interface temperature.
    /// Caller-supplied laws retain caller authority; no material source is inferred.
    pub law: FrictionLaw,
    /// The two original material conductivities at their own cell temperatures.
    pub conductivity_w_m_k: [f64; 2],
    /// Other heat/conduction only. Friction work is supplied by mechanics.
    pub transport: StepCoefficients<P>,
}

/// Bounded nonlinear work and independent thermal/whole-system energy gates.
#[derive(Clone, Copy, Debug)]
pub struct ModalFrictionHeatSolve {
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scale [J/K], not an added heat capacity or store.
    /// The last Newton coordinate is this scale times interface-temperature change.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    /// The inner forcing floor must respect finite-difference accuracy; an
    /// analytic-Jacobian forcing floor can stall below the directional-probe
    /// precision. Outer convergence and every independent energy gate still
    /// apply. The resolved moving-contact fixture uses a 1e-6 forcing floor.
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum ModalFrictionHeatError<E> {
    Invalid(&'static str),
    Mechanical(ModalFrictionError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for ModalFrictionHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<ModalFrictionError> for ModalFrictionHeatError<E> {
    fn from(error: ModalFrictionError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<NewtonError> for ModalFrictionHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ModalFrictionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "modal friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ModalFrictionHeatError<E> {}

fn numerical<E>(
    error: FiniteDifferenceError<ModalFrictionHeatError<E>>,
) -> ModalFrictionHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => ModalFrictionHeatError::Invalid(message),
    }
}

/// One accepted mechanical sample and its original finite-network heat account.
pub struct HeatedModalFrictionFrame<P> {
    pub mechanical: FrictionModalFrame,
    pub heat: FrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Mechanical storage/loss plus caloric storage, minus outside supplies.
    /// Physical friction transfer cancels internally exactly once.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: FrictionModalFrame,
    nonlinear: Option<nonlinear::NonlinearFrictionFrame>,
    law: FrictionLaw,
    binding: Option<interface::DryFrictionBinding>,
    branch: ContactBranch,
    load_line: Option<StribeckLoadLine>,
    source: P,
}

impl FrictionModalSystem {
    /// Advance mechanics and the existing thermal network in one transaction.
    ///
    /// The route's A/B cells must represent this port's actual contacting
    /// bodies. Resolve conductivity from their original material sources and
    /// friction from the COMMON interface temperature supplied to `resolve`.
    /// The callback may stage source receipts but must not publish history.
    ///
    /// A retained dry-interface card keeps its exact law and applicability;
    /// its temperature is updated through that same card. Replacing its source
    /// or coefficients is refused. Unbound ports can use explicit temperature-
    /// dependent source callbacks, which retain their original caller authority.
    /// The source must not insert the friction heat again in its cell powers.
    ///
    /// Thermal unknowns use sorted CellId order; the temperature guess is a
    /// numerical starting point, never an imposed interface temperature.
    /// Backward-Euler heat accuracy requires timestep and mesh refinement.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_friction_heat<'m, P, E>(
        &mut self,
        thermal: &mut DepositionSimulation<'m>,
        route: &FrictionHeatRoute,
        external: &[f64],
        drive: FrictionDrive,
        controls: ModalFrictionHeatSolve,
        energy_guess_j: &[f64],
        temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<ModalHeatSource<P>, E>,
        gate: &CancelGate,
    ) -> Result<HeatedModalFrictionFrame<P>, ModalFrictionHeatError<E>> {
        if self.adhesive.is_some() {
            return Err(ModalFrictionHeatError::Invalid(
                "dry-law thermal refresh cannot replace a frozen adhesive normal/shear source",
            ));
        }
        let checkpoint = || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let sample = self.samples_rendered();
        let next_sample = sample
            .checked_add(1)
            .ok_or(ModalFrictionHeatError::Invalid(
                "modal sample clock overflow",
            ))?;
        let start = sample as f64 * self.sample_period_s();
        let target = next_sample as f64 * self.sample_period_s();
        if thermal.time_s() != start
            || !target.is_finite()
            || target <= start
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
                "synchronized clocks and finite bounded friction/heat coordinates required",
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
                    "incomplete friction/heat iterate",
                ));
            }
            let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
            step.evaluate(
                &x[..cells],
                temperature,
                |endpoint, interface_temperature| {
                    let ModalHeatSource {
                        law,
                        conductivity_w_m_k,
                        transport,
                    } = resolve(endpoint, interface_temperature)
                        .map_err(ModalFrictionHeatError::Source)?;
                    let StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    let staged = owner.borrow_mut().stage_at_temperature(
                        law,
                        interface_temperature,
                        external,
                        drive,
                        payload,
                        gate,
                    )?;
                    Ok(FrictionHeatInput {
                        dissipated_j: staged.mechanical.friction_dissipation_j,
                        conductivity_w_m_k,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: staged,
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
        // Last evaluation leaves the ORIGINAL mechanical owner's scratch state
        // at exactly the same endpoint accepted by the thermal owner below.
        let trial = evaluate(&solver.x)?;
        let combined =
            trial.payload().mechanical.energy_residual_j + trial.report().total_energy_residual_j;
        if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
            return Err(ModalFrictionHeatError::Invalid(
                "joint mechanical/thermal work does not balance",
            ));
        }
        let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
        checkpoint()?;
        drop(problem);
        let owner = owner.into_inner();
        // No fallible work or user callback follows this publication boundary.
        owner.port.law = staged.law;
        owner.interface = staged.binding;
        owner.load_line = staged.load_line;
        owner.staged_branch = staged.branch;
        owner.publish_candidate(staged.mechanical, staged.nonlinear);
        *thermal = next_thermal;
        Ok(HeatedModalFrictionFrame {
            mechanical: staged.mechanical,
            heat,
            source: staged.source,
            nonlinear: result.report,
            trials: trials.get(),
            combined_energy_residual_j: combined,
        })
    }

    // Temporarily install only the endpoint source. Always restore accepted
    // law/binding/branch cache, including when the physical proposal refuses.
    // Each thermal trial gets a fresh load line: its parameters can change with
    // temperature even when normal load and mechanical mobility are unchanged.
    fn stage_at_temperature<P>(
        &mut self,
        law: FrictionLaw,
        temperature: f64,
        external: &[f64],
        drive: FrictionDrive,
        source: P,
        gate: &CancelGate,
    ) -> Result<Staged<P>, ModalFrictionError> {
        law.static_coefficient()?;
        let binding = if let Some(old) = &self.interface {
            if law != *old.card().friction_law() {
                return Err(invalid(
                    "thermal callback cannot replace a retained dry-interface source",
                )
                .into());
            }
            Some(interface::DryFrictionBinding::new(
                old.card().clone(),
                temperature,
                old.nominal_area_m2(),
            )?)
        } else {
            None
        };
        let saved_law = std::mem::replace(&mut self.port.law, law.clone());
        let saved_binding = std::mem::replace(&mut self.interface, binding.clone());
        let saved_line = self.load_line.take();
        let saved_branch = self.staged_branch;
        let result = self.stage_candidate(external, drive, Some(gate));
        let line = self.load_line.take();
        let branch = self.staged_branch;
        self.port.law = saved_law;
        self.interface = saved_binding;
        self.load_line = saved_line;
        self.staged_branch = saved_branch;
        let (mechanical, nonlinear) = result?;
        Ok(Staged {
            mechanical,
            nonlinear,
            law,
            binding,
            branch,
            load_line: line,
            source,
        })
    }
}
