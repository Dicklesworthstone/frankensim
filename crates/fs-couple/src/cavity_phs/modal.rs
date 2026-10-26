//! Finite gas pressure/work coupled to source-derived elastic wall motion.
//!
//! The original gas owner advances heat and volume with its discrete midpoint
//! pressure. The original solid modal owner advances the same held pressure
//! through its exact ZOH map. A scalar volume residual joins their endpoints;
//! the existing Newton-Krylov owner solves it using the exact discrete tangent.

use std::cell::{Cell, RefCell};

use fs_exec::Cx;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::midpoint::{
    GasHeatVolumeCandidate, GasHeatVolumeStep,
};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_solid::reduce::thermal::motion::{
    ThermalValveModalMotion, ThermalValveMotionError, ThermalValveMotionFrame,
};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
    SolverRunProgress,
};

/// Source and numerical controls; the native modal owner supplies the timestep.
#[derive(Clone, Copy, Debug)]
pub struct ModalGasSolve {
    pub gas: GasEnergySolve,
    /// Every trial volume must lie in this positive interval [m3]. This is a
    /// source-domain bound, not a promise that a physical root is bracketed.
    pub volume_interval_m3: [f64; 2],
    pub volume_tolerance_m3: f64,
    /// Independent pressure-transfer and enclosing energy allowances [J].
    pub energy_tolerance_j: f64,
    /// Gas plus native-motion evaluations, including Jacobian and final checks.
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Residual and Newton coordinate are scaled by the current gas volume.
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum ModalGasError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Motion(ThermalValveMotionError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TimeResolution(fs_phs::linear_exact::ResolvedBandError),
    TrialBudget,
    Cancelled,
}
impl core::fmt::Display for ModalGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite gas/modal wall coupling: {self:?}")
    }
}
impl std::error::Error for ModalGasError {}
impl From<GasEnergyError> for ModalGasError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(error)
        }
    }
}
impl From<ThermalValveMotionError> for ModalGasError {
    fn from(error: ThermalValveMotionError) -> Self {
        match error {
            ThermalValveMotionError::Cancelled => Self::Cancelled,
            other => Self::Motion(other),
        }
    }
}
impl From<NewtonError> for ModalGasError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
fn invalid(what: &'static str) -> ModalGasError {
    ModalGasError::Invalid(what)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cx: &Cx<'_>) -> Result<(), ModalGasError> {
    cx.checkpoint().map_err(|_| ModalGasError::Cancelled)
}

/// One closed gas and its material/mesh-derived pressure boundary.
///
/// The retained pressure footprint defines an affine swept-volume chart. The
/// supplied gas volume is the actual cavity volume at the supplied modal state,
/// including its equilibrium displacement. Positive pressure opens this volume.
/// Geometry and the thermal equilibrium/basis stay fixed throughout this run;
/// no hidden modal rebuild, equilibrium reset or gas-volume reset is performed.
#[derive(Clone, Debug)]
pub struct GasValveMotion {
    gas: ClosedGasParcel,
    motion: ThermalValveModalMotion,
    reference_volume_m3: f64,
    reference_coordinates: Vec<f64>,
}

/// Original stage accounts and independently checked interconnection balances.
#[derive(Clone, Debug)]
pub struct ModalGasStep {
    pub gas: GasHeatVolumeStep,
    /// Shared source-bearing pressure effort, also available to acoustic owners.
    pub midpoint: ResolvedThermochemicalGasState,
    pub motion: ThermalValveMotionFrame,
    pub modal_energy_change_j: f64,
    /// Work against constant exterior pressure on the actual wall sweep [J].
    pub exterior_pressure_work_j: f64,
    /// Native modal input plus exterior work minus native gas pressure work [J].
    pub pressure_work_residual_j: f64,
    /// Gas storage + modal storage + original modal loss + exterior work - Q.
    /// Structural damping leaves this enclosure; it is not credited to the gas.
    pub total_energy_residual_j: f64,
    /// Gas volume minus the original anchored modal volume chart [m3].
    pub volume_residual_m3: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

fn modal_energy(motion: &ThermalValveModalMotion) -> f64 {
    motion
        .model()
        .states()
        .iter()
        .zip(motion.model().modes())
        .map(|(x, mode)| {
            let elastic = mode.angular_frequency_rad_s * x.displacement_m_sqrt_kg;
            0.5 * (x.velocity_m_sqrt_kg_per_s * x.velocity_m_sqrt_kg_per_s + elastic * elastic)
        })
        .sum()
}

// Conservative rate bound for diagonal structural stiffness plus the gas
// rank-one stiffness (rho*c^2/V) a a^T. Exact structural ZOH does not make the
// nonlinear held-pressure interconnection the full coupled exponential.
fn admit_clock(
    motion: &ThermalValveModalMotion,
    gas: &ResolvedThermochemicalGasState,
    volume: f64,
) -> Result<(), ModalGasError> {
    let modes = motion.model().modes();
    let stiffness = modes
        .iter()
        .map(|m| m.angular_frequency_rad_s.powi(2))
        .fold(0.0, f64::max);
    let damping = modes
        .iter()
        .map(|m| 2.0 * m.damping_ratio * m.angular_frequency_rad_s)
        .fold(0.0, f64::max);
    let area2: f64 = motion
        .valve()
        .valve()
        .mass_normalized_pressure_load()
        .iter()
        .map(|a| a * a)
        .sum();
    let gas = gas.state();
    let rate = (stiffness + gas.density * gas.sound_speed * gas.sound_speed * area2 / volume)
        .sqrt()
        + damping;
    let dt = motion.model().sample_period_s();
    if !rate.is_finite()
        || rate * dt >= core::f64::consts::PI * motion.model().budget().nyquist_guard_fraction
    {
        return Err(invalid(
            "gas-loaded modes exceed the original sample-frequency guard",
        ));
    }
    fs_phs::linear_exact::admit_midpoint_band(rate, dt).map_err(ModalGasError::TimeResolution)
}

struct Trial {
    gas: GasHeatVolumeCandidate,
    motion: ThermalValveModalMotion,
    frame: ThermalValveMotionFrame,
    residual_m3: f64,
    tangent: f64,
}

// Private physical residual adapter, not another nonlinear iteration owner.
// Preserve the first source refusal; a failed trial cannot become convergence.
struct VolumeProblem<F> {
    evaluate: F,
    scale: f64,
    failed: Cell<bool>,
    failure: RefCell<Option<ModalGasError>>,
}
impl<F: Fn(f64) -> Result<Trial, ModalGasError>> VolumeProblem<F> {
    fn sample(&self, coordinate: f64) -> Option<Trial> {
        if self.failed.get() {
            return None;
        }
        match (self.evaluate)(coordinate) {
            Ok(trial) => Some(trial),
            Err(error) => {
                self.failed.set(true);
                *self.failure.borrow_mut() = Some(error);
                None
            }
        }
    }
}
impl<F: Fn(f64) -> Result<Trial, ModalGasError>> NonlinearProblem for VolumeProblem<F> {
    fn dimension(&self) -> usize {
        1
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if x.len() == 1 && output.len() == 1 {
            if let Some(trial) = self.sample(x[0]) {
                output[0] = trial.residual_m3 / self.scale;
            }
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if x.len() == 1 && direction.len() == 1 && output.len() == 1 {
            if let Some(trial) = self.sample(x[0]) {
                output[0] = trial.tangent * direction[0];
            }
        }
    }
}

impl GasValveMotion {
    /// Bind a finite source gas to an accepted native elastic motion.
    /// The gas parcel supplies the physical cavity volume; no area or effective
    /// mass is copied into an independent piston model.
    pub fn new(
        gas: ClosedGasParcel,
        motion: ThermalValveModalMotion,
        cx: &Cx<'_>,
    ) -> Result<Self, ModalGasError> {
        poll(cx)?;
        admit_clock(&motion, gas.state(), gas.volume_m3())?;
        let reference_volume_m3 = gas.volume_m3();
        let reference_coordinates = motion
            .model()
            .states()
            .iter()
            .map(|state| state.displacement_m_sqrt_kg)
            .collect();
        Ok(Self {
            gas,
            motion,
            reference_volume_m3,
            reference_coordinates,
        })
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn motion(&self) -> &ThermalValveModalMotion {
        &self.motion
    }
    /// Preserve original histories for a separately admitted source update.
    /// Changed basis/equilibrium needs native state transfer and explicit
    /// parameter/geometry work accounting before binding another volume chart.
    pub fn into_parts(self) -> (ClosedGasParcel, ThermalValveModalMotion) {
        (self.gas, self.motion)
    }

    /// Solve finite gas heat/work and actual elastic boundary motion together.
    ///
    /// Heat is external energy into the gas over one native modal sample.
    /// Exterior pressure is constant and absolute; zero denotes vacuum.
    /// Every trial starts from the same accepted gas and motion. Newton uses
    /// (V1-V0)/V0, with exact tangent
    /// 1 - (sum_i a_i^2 bq_i) d(p_mid)/d(V1), from native ZOH coefficients and
    /// the gas owner's discrete derivative, not a fixed-gamma approximation.
    ///
    /// Gas is spatially uniform, frozen-composition and closed. This closes
    /// pressure/work feedback, not heat transfer to the wall. The fixed source
    /// epoch retains its original structural loss and observation transfers;
    /// exterior radiation transfers are not refreshed from the new gas here.
    /// Reciprocal thermoelastic heat, reaction, nozzles, contact and finite
    /// geometry are absent. The held-pressure interconnection needs temporal
    /// refinement and is not an exact coupled exponential or entropy theorem.
    /// Failure publishes neither participant nor its clock.
    pub fn step(
        &mut self,
        heat_j: f64,
        exterior_pressure_pa: f64,
        controls: ModalGasSolve,
        cx: &Cx<'_>,
    ) -> Result<ModalGasStep, ModalGasError> {
        poll(cx)?;
        let v0 = self.gas.volume_m3();
        let [lower, upper] = controls.volume_interval_m3;
        if !heat_j.is_finite()
            || !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
            || !positive(lower)
            || !positive(upper)
            || lower >= upper
            || v0 < lower
            || v0 > upper
            || !positive(controls.volume_tolerance_m3)
            || !positive(controls.energy_tolerance_j)
            || controls.maximum_trials == 0
            || !(1..=128).contains(&controls.maximum_newton_iterations)
        {
            return Err(invalid(
                "finite heat, nonnegative exterior pressure and bounded volume/work controls required",
            ));
        }
        controls.gas.validate()?;
        controls.newton.validate()?;
        admit_clock(&self.motion, self.gas.state(), v0)?;
        let area = self.motion.valve().valve().mass_normalized_pressure_load();
        let compliance: f64 = area
            .iter()
            .enumerate()
            .map(|(i, a)| a * a * self.motion.model().probe_linear_map(i).4)
            .sum();
        let initial_energy = modal_energy(&self.motion);
        if !compliance.is_finite() || compliance < 0.0 || !initial_energy.is_finite() {
            return Err(invalid(
                "native held-pressure compliance or modal storage is unrepresentable",
            ));
        }
        let trials = Cell::new(0usize);
        let evaluate = |fraction: f64| {
            poll(cx)?;
            if trials.get() >= controls.maximum_trials {
                return Err(ModalGasError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let volume = v0 + v0 * fraction;
            if !volume.is_finite() || volume < lower || volume > upper {
                return Err(invalid(
                    "nonlinear trial left the declared gas volume interval",
                ));
            }
            let gas = self
                .gas
                .midpoint_heat_volume(volume, heat_j, controls.gas, || cx.checkpoint().is_err())?;
            admit_clock(&self.motion, gas.midpoint(), 0.5 * v0 + 0.5 * volume)?;
            admit_clock(&self.motion, gas.endpoint().state(), volume)?;
            let mut motion = self.motion.clone();
            let frame =
                motion.step_pressure(gas.midpoint().state().pressure - exterior_pressure_pa, cx)?;
            let swept: f64 = area
                .iter()
                .zip(&self.reference_coordinates)
                .zip(motion.model().states())
                .map(|((&a, &reference), state)| a * (state.displacement_m_sqrt_kg - reference))
                .sum();
            // Anchor geometry at construction; repeated tolerances must not
            // silently accumulate a separate gas-volume trajectory.
            let residual_m3 = (volume - self.reference_volume_m3) - swept;
            let tangent = 1.0 - compliance * gas.step().pressure_final_volume_derivative_pa_m3;
            if !residual_m3.is_finite() || !positive(tangent) {
                return Err(invalid(
                    "coupled volume residual or exact discrete tangent is unrepresentable",
                ));
            }
            Ok(Trial {
                gas,
                motion,
                frame,
                residual_m3,
                tangent,
            })
        };
        let problem = VolumeProblem {
            evaluate: &evaluate,
            scale: v0,
            failed: Cell::new(false),
            failure: RefCell::new(None),
        };
        let construction = NewtonKrylovState::new(&problem, vec![0.0], controls.newton);
        if let Some(error) = problem.failure.borrow_mut().take() {
            return Err(error);
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                cx.checkpoint().is_err() || problem.failed.get()
            });
        if let Some(error) = problem.failure.borrow_mut().take() {
            return Err(error);
        }
        poll(cx)?;
        if result.progress != SolverRunProgress::Complete {
            return Err(ModalGasError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(ModalGasError::NotConverged(result.report));
        }
        let Trial {
            gas,
            motion,
            frame,
            residual_m3,
            ..
        } = evaluate(solver.x[0])?;
        let midpoint = gas.midpoint().clone();
        let modal_energy_change_j = frame.modal.total_modal_energy_j - initial_energy;
        // The exterior acts on the actual wall sweep, not a potentially
        // different approximate gas endpoint. Independently gate that transfer.
        let exterior_pressure_work_j = exterior_pressure_pa * frame.swept_volume_change_m3;
        let pressure_work_residual_j =
            (frame.modal.input_work_j + exterior_pressure_work_j) - gas.step().work_by_gas_j;
        let total_energy_residual_j = gas.step().stored_energy_change_j
            + modal_energy_change_j
            + frame.modal.viscous_dissipation_j
            + exterior_pressure_work_j
            - heat_j;
        if residual_m3.abs() > controls.volume_tolerance_m3
            || [
                modal_energy_change_j,
                exterior_pressure_work_j,
                pressure_work_residual_j,
                total_energy_residual_j,
            ]
            .iter()
            .any(|x| !x.is_finite())
            || pressure_work_residual_j.abs() > controls.energy_tolerance_j
            || total_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "native gas volume, pressure work or enclosing energy did not close",
            ));
        }
        poll(cx)?;
        let (gas, gas_report) = gas.into_parts();
        self.gas = gas;
        self.motion = motion;
        Ok(ModalGasStep {
            gas: gas_report,
            midpoint,
            motion: frame,
            modal_energy_change_j,
            exterior_pressure_work_j,
            pressure_work_residual_j,
            total_energy_residual_j,
            volume_residual_m3: residual_m3,
            nonlinear: result.report,
            trials: trials.get(),
        })
    }
}
