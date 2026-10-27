//! Same-endpoint chemistry and original material-derived elastic wall motion.
//!
//! Eliminate volume with the native held-pressure modal response. At each
//! pressure trial, its anchored wall sweep supplies V1 and W=p_trial*(V1-V0).
//! The original implicit mechanism then solves at U1=U0+Q-W and V1. The outer
//! native Newton residual is p_EOS(U1,V1,n1)-p_trial. No independent piston,
//! reaction law, caloric law or nonlinear iteration owner is introduced.
//!
//! Rates and gas pressure use backward-Euler endpoints; the structural held
//! pressure response remains native ZOH. The coupled method is first order,
//! not the exact exponential of a reacting fluid/structure system.

use super::super::{
    GasValveMotion, ModalGasError, ModalGasSolve, admit_clock, modal_energy, poll, positive,
};
use super::ReactiveModalGasError;
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve,
    boundary::{BoundaryGasReactionCandidate, BoundaryGasReactionStep, GasReactionBoundary},
};
use fs_exec::Cx;
use fs_solid::reduce::thermal::motion::{ThermalValveModalMotion, ThermalValveMotionFrame};
use fs_solver::{NewtonKrylovState, NewtonReport, NonlinearProblem, SolverRunProgress};
use std::cell::{Cell, RefCell};

type Error = ReactiveModalGasError;
fn invalid(message: &'static str) -> Error {
    Error::Invalid(message)
}

/// Native motion/caloric gates plus explicit pressure-root controls.
#[derive(Clone, Copy, Debug)]
pub struct CoupledReactiveModalGasSolve {
    pub mechanical: ModalGasSolve,
    /// Positive endpoint pressure interval containing the initial gas pressure.
    /// Every nonlinear/derivative probe must remain in the admitted source domain.
    pub pressure_interval_pa: [f64; 2],
    pub pressure_tolerance_pa: f64,
    /// Relative centered/one-sided pressure stencil, in [32 epsilon, 0.01].
    /// The caller must resolve source/inversion noise and refine independently.
    pub pressure_difference_fraction: f64,
}

/// One gas/species/pressure/work and native structural accounting window.
#[derive(Clone, Debug)]
pub struct CoupledReactiveModalGasStep {
    pub reaction: BoundaryGasReactionStep,
    pub motion: ThermalValveMotionFrame,
    pub held_pressure_pa: f64,
    pub pressure_residual_pa: f64,
    pub modal_energy_change_j: f64,
    pub exterior_pressure_work_j: f64,
    /// Native modal input plus exterior work minus the gas work debit [J].
    pub pressure_work_residual_j: f64,
    /// EOS/held-pressure mismatch times actual boundary sweep [J].
    pub pressure_eos_work_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub volume_residual_m3: f64,
    pub nonlinear: NewtonReport,
    /// All outer chemistry/motion probes, including Jacobian/final probes.
    pub trials: usize,
}

/// Publication requires all original source, species, motion, work and entropy gates.
#[derive(Clone, Debug)]
pub struct CoupledReactiveModalGasCandidate {
    motion: GasValveMotion,
    step: CoupledReactiveModalGasStep,
}
impl CoupledReactiveModalGasCandidate {
    pub const fn motion(&self) -> &GasValveMotion {
        &self.motion
    }
    pub const fn step(&self) -> &CoupledReactiveModalGasStep {
        &self.step
    }
    pub fn into_parts(self) -> (GasValveMotion, CoupledReactiveModalGasStep) {
        (self.motion, self.step)
    }
}

struct Trial {
    chemistry: BoundaryGasReactionCandidate,
    motion: ThermalValveModalMotion,
    frame: ThermalValveMotionFrame,
    pressure_pa: f64,
    residual_pa: f64,
}

// This is a residual/Jacobian adapter for the existing Newton-Krylov owner.
// A failed source or work probe is never converted into numerical convergence.
struct PressureProblem<F> {
    evaluate: F,
    scale_pa: f64,
    interval_pa: [f64; 2],
    difference_fraction: f64,
    failed: Cell<bool>,
    failure: RefCell<Option<Error>>,
    derivative: RefCell<Option<(u64, f64)>>,
}
impl<F: Fn(f64) -> Result<Trial, Error>> PressureProblem<F> {
    fn fail(&self, error: Error) {
        if !self.failed.replace(true) {
            *self.failure.borrow_mut() = Some(error);
        }
    }
    fn sample(&self, coordinate: f64) -> Option<Trial> {
        if self.failed.get() {
            return None;
        }
        match (self.evaluate)(coordinate) {
            Ok(trial) => Some(trial),
            Err(error) => {
                self.fail(error);
                None
            }
        }
    }
    fn slope(&self, coordinate: f64) -> Option<f64> {
        if self.failed.get() {
            return None;
        }
        if let Some((key, value)) = *self.derivative.borrow() {
            if key == coordinate.to_bits() {
                return Some(value);
            }
        }
        let pressure = self.scale_pa.mul_add(coordinate, self.scale_pa);
        let h = self.difference_fraction * pressure.abs().max(self.scale_pa);
        let lo = ((pressure - h).max(self.interval_pa[0]) - self.scale_pa) / self.scale_pa;
        let hi = ((pressure + h).min(self.interval_pa[1]) - self.scale_pa) / self.scale_pa;
        if !lo.is_finite() || !hi.is_finite() || lo >= hi {
            self.fail(invalid(
                "coupled pressure Jacobian stencil is unrepresentable",
            ));
            return None;
        }
        let left = self.sample(lo)?.residual_pa;
        let right = self.sample(hi)?.residual_pa;
        let derivative = ((right - left) / self.scale_pa) / (hi - lo);
        if !derivative.is_finite() || derivative == 0.0 {
            self.fail(invalid(
                "coupled pressure Jacobian is singular or unrepresentable",
            ));
            return None;
        }
        *self.derivative.borrow_mut() = Some((coordinate.to_bits(), derivative));
        Some(derivative)
    }
}
impl<F: Fn(f64) -> Result<Trial, Error>> NonlinearProblem for PressureProblem<F> {
    fn dimension(&self) -> usize {
        1
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if x.len() == 1 && output.len() == 1 {
            if let Some(trial) = self.sample(x[0]) {
                output[0] = trial.residual_pa / self.scale_pa;
            }
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if x.len() == 1 && direction.len() == 1 && output.len() == 1 {
            if let Some(slope) = self.slope(x[0]) {
                output[0] = slope * direction[0];
            }
        }
    }
}

impl GasValveMotion {
    /// Resolve chemistry, heat, pressure work and native wall motion together.
    ///
    /// Volume is determined by the original mesh-derived pressure footprint and
    /// anchored modal coordinates. Every trial restarts both physical owners at
    /// their accepted initial state. Its original implicit mechanism evaluates
    /// all channel rates at the proposed final composition, volume and energy.
    /// Heat is external only; formation energy is already inside gas storage.
    /// Native structural damping leaves this enclosure and is not gas heating.
    ///
    /// Requires the original coupled implicit mechanism and one shared gas
    /// inversion policy. Chemistry-only subcycling is refused; refine the whole
    /// coupled sample. Outer work is bounded by mechanical.maximum_trials, and
    /// each inner mechanism retains its own reaction/caloric limits. A finite
    /// difference is used only for the root; physical endpoint residuals decide
    /// acceptance independently. No uniqueness or unconditional convergence is
    /// claimed. The source-fixed structural basis, equilibrium and observations
    /// remain intact; no flame, shock, wall heat transfer, moving mesh rebuild
    /// or composition-dependent acoustic radiation transfer is added.
    /// The native frozen-composition gas-loaded clock gate is retained; it is
    /// not an error estimate for the complete reacting system.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_reactive_step_candidate(
        &self,
        reaction: &GasReaction,
        heat_j: f64,
        exterior_pressure_pa: f64,
        chemistry: GasReactionSolve,
        controls: CoupledReactiveModalGasSolve,
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveModalGasCandidate, Error> {
        poll(cx)?;
        let native = controls.mechanical;
        let p0 = self.gas.state().state().pressure;
        let v0 = self.gas.volume_m3();
        let [pmin, pmax] = controls.pressure_interval_pa;
        let [vmin, vmax] = native.volume_interval_m3;
        if !heat_j.is_finite()
            || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
            || !positive(pmin)
            || !positive(pmax)
            || pmin >= pmax
            || p0 < pmin
            || p0 > pmax
            || !positive(vmin)
            || !positive(vmax)
            || vmin >= vmax
            || v0 < vmin
            || v0 > vmax
            || !positive(controls.pressure_tolerance_pa)
            || !controls.pressure_difference_fraction.is_finite()
            || !(32.0 * f64::EPSILON..=0.01).contains(&controls.pressure_difference_fraction)
            || !positive(native.volume_tolerance_m3)
            || !positive(native.energy_tolerance_j)
            || native.maximum_trials == 0
            || !(1..=128).contains(&native.maximum_newton_iterations)
            || chemistry.gas != native.gas
        {
            return Err(invalid(
                "bounded pressure/volume controls and one shared reaction/gas policy required",
            ));
        }
        native.gas.validate().map_err(ModalGasError::from)?;
        native.newton.validate().map_err(ModalGasError::from)?;
        admit_clock(&self.motion, self.gas.state(), v0)?;
        let dt = self.motion.model().sample_period_s();
        let area = self.motion.valve().valve().mass_normalized_pressure_load();
        let initial_energy = modal_energy(&self.motion);
        if !initial_energy.is_finite() {
            return Err(invalid("initial native modal energy is unrepresentable"));
        }
        let trials = Cell::new(0usize);
        let evaluate = |coordinate: f64| -> Result<Trial, Error> {
            poll(cx)?;
            if trials.get() >= native.maximum_trials {
                return Err(ModalGasError::TrialBudget.into());
            }
            trials.set(trials.get() + 1);
            let pressure_pa = p0.mul_add(coordinate, p0);
            if !pressure_pa.is_finite() || pressure_pa < pmin || pressure_pa > pmax {
                return Err(invalid("coupled pressure trial left its declared interval"));
            }
            let mut motion = self.motion.clone();
            let frame = motion
                .step_pressure(pressure_pa - exterior_pressure_pa, cx)
                .map_err(ModalGasError::from)?;
            let swept: f64 = area
                .iter()
                .zip(&self.reference_coordinates)
                .zip(motion.model().states())
                .map(|((&a, &reference), state)| a * (state.displacement_m_sqrt_kg - reference))
                .sum();
            let volume = self.reference_volume_m3 + swept;
            if !volume.is_finite() || volume < vmin || volume > vmax {
                return Err(invalid(
                    "native reacting wall motion left the declared gas volume interval",
                ));
            }
            let work_by_gas_j = pressure_pa * (volume - v0);
            if !work_by_gas_j.is_finite() || (volume != v0 && work_by_gas_j == 0.0) {
                return Err(invalid("reacting gas boundary work is unrepresentable"));
            }
            let chemistry = reaction.mechanical_implicit_trial(
                &self.gas,
                dt,
                GasReactionBoundary {
                    final_volume_m3: volume,
                    heat_j,
                    work_by_gas_j,
                },
                chemistry,
                || cx.checkpoint().is_err(),
            )?;
            admit_clock(&self.motion, chemistry.gas().state(), volume)?;
            let residual_pa = chemistry.gas().state().state().pressure - pressure_pa;
            if !residual_pa.is_finite() {
                return Err(invalid(
                    "coupled reacting pressure residual is unrepresentable",
                ));
            }
            Ok(Trial {
                chemistry,
                motion,
                frame,
                pressure_pa,
                residual_pa,
            })
        };
        let problem = PressureProblem {
            evaluate: &evaluate,
            scale_pa: p0,
            interval_pa: controls.pressure_interval_pa,
            difference_fraction: controls.pressure_difference_fraction,
            failed: Cell::new(false),
            failure: RefCell::new(None),
            derivative: RefCell::new(None),
        };
        let construction = NewtonKrylovState::new(&problem, vec![0.0], native.newton);
        if let Some(error) = problem.failure.borrow_mut().take() {
            return Err(error);
        }
        let mut solver = construction.map_err(ModalGasError::from)?;
        let result = solver.run_controlled(&problem, native.maximum_newton_iterations, &mut || {
            cx.checkpoint().is_err() || problem.failed.get()
        });
        if let Some(error) = problem.failure.borrow_mut().take() {
            return Err(error);
        }
        poll(cx)?;
        if result.progress != SolverRunProgress::Complete {
            return Err(ModalGasError::Interrupted(result.progress).into());
        }
        if !result.report.converged {
            return Err(ModalGasError::NotConverged(result.report).into());
        }
        let Trial {
            chemistry: endpoint,
            motion,
            frame,
            pressure_pa,
            residual_pa,
        } = evaluate(solver.x[0])?;
        let (gas, reaction) = endpoint.into_parts();
        let actual_volume: f64 = self.reference_volume_m3
            + area
                .iter()
                .zip(&self.reference_coordinates)
                .zip(motion.model().states())
                .map(|((&a, &reference), state)| a * (state.displacement_m_sqrt_kg - reference))
                .sum::<f64>();
        let volume_residual_m3 = gas.volume_m3() - actual_volume;
        let modal_energy_change_j = frame.modal.total_modal_energy_j - initial_energy;
        let exterior_pressure_work_j = exterior_pressure_pa * frame.swept_volume_change_m3;
        let pressure_work_residual_j =
            (frame.modal.input_work_j + exterior_pressure_work_j) - reaction.boundary.work_by_gas_j;
        let pressure_eos_work_residual_j = residual_pa * frame.swept_volume_change_m3;
        let total_energy_residual_j = reaction.reaction.stored_energy_change_j
            + modal_energy_change_j
            + frame.modal.viscous_dissipation_j
            + exterior_pressure_work_j
            - heat_j;
        if residual_pa.abs() > controls.pressure_tolerance_pa
            || volume_residual_m3.abs() > native.volume_tolerance_m3
            || [
                modal_energy_change_j,
                exterior_pressure_work_j,
                pressure_work_residual_j,
                pressure_eos_work_residual_j,
                total_energy_residual_j,
                volume_residual_m3,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || pressure_work_residual_j.abs() > native.energy_tolerance_j
            || pressure_eos_work_residual_j.abs() > native.energy_tolerance_j
            || total_energy_residual_j.abs() > native.energy_tolerance_j
        {
            return Err(invalid(
                "reacting pressure, anchored volume or gas/modal work did not close",
            ));
        }
        // Outer work probes deliberately defer this physical inequality. It is
        // mandatory now that the same endpoint pressure/work residuals close.
        if reaction.entropy_residual_j_k < -chemistry.entropy_tolerance_j_k {
            return Err(GasReactionError::EntropyDecrease {
                change_j_k: reaction.entropy_residual_j_k,
                tolerance_j_k: chemistry.entropy_tolerance_j_k,
            }
            .into());
        }
        poll(cx)?;
        let mut candidate = self.clone();
        candidate.gas = gas;
        candidate.motion = motion;
        Ok(CoupledReactiveModalGasCandidate {
            motion: candidate,
            step: CoupledReactiveModalGasStep {
                reaction,
                motion: frame,
                held_pressure_pa: pressure_pa,
                pressure_residual_pa: residual_pa,
                modal_energy_change_j,
                exterior_pressure_work_j,
                pressure_work_residual_j,
                pressure_eos_work_residual_j,
                total_energy_residual_j,
                volume_residual_m3,
                nonlinear: result.report,
                trials: trials.get(),
            },
        })
    }

    /// Publish the native gas and wall histories only after the complete coupled
    /// chemical/mechanical candidate passes. Refusal preserves deterministic retry.
    #[allow(clippy::too_many_arguments)]
    pub fn step_coupled_reactive(
        &mut self,
        reaction: &GasReaction,
        heat_j: f64,
        exterior_pressure_pa: f64,
        chemistry: GasReactionSolve,
        controls: CoupledReactiveModalGasSolve,
        cx: &Cx<'_>,
    ) -> Result<CoupledReactiveModalGasStep, Error> {
        let (candidate, step) = self
            .coupled_reactive_step_candidate(
                reaction,
                heat_j,
                exterior_pressure_pa,
                chemistry,
                controls,
                cx,
            )?
            .into_parts();
        poll(cx)?;
        *self = candidate;
        Ok(step)
    }
}
