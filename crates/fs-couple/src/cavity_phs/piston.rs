//! A moving inertial wall coupled to a closed, source-backed gas cavity.
//!
//! The gas owns caloric state, EOS, source support, transport and acoustic
//! properties. This module owns ONE translational boundary coordinate and its
//! momentum. Both use the same midpoint pressure: there is no prescribed-motion
//! animation, nominal-air load, or independent pressure/work approximation.

/// Shared source-resolved rigid-body inertia coupled to the same gas dynamics.
pub mod material;
/// Reciprocal finite-body enthalpy/phase exchange on the same gas/momentum step.
pub mod thermal;

use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::midpoint::{
    GasHeatVolumeCandidate, GasHeatVolumeStep,
};

/// Closed constant-area cavity with a translating inertial end wall.
/// Volume is the single position owner: the effective wall coordinate is V/A.
/// Other rigid-body coordinates, contact stops and crank linkages are not modeled.
#[derive(Clone, Debug, PartialEq)]
pub struct GasPiston {
    gas: ClosedGasParcel,
    mass_kg: f64,
    area_m2: f64,
    velocity_m_s: f64,
}

/// Explicit external loading during one physical timestep.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonInput {
    /// Positive elapsed physical time [s].
    pub duration_s: f64,
    /// Integrated heat INTO the gas [J]; not a combustion model or watts.
    pub heat_j: f64,
    /// Constant absolute pressure outside the piston [Pa], zero for vacuum.
    pub exterior_pressure_pa: f64,
    /// Additional constant force [N], positive towards increasing gas volume.
    /// Does not include the gas or exterior-pressure load, which are added here.
    pub external_force_n: f64,
}

/// Nonlinear work limits and independent physical residual admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonSolve {
    /// Final-volume bracket [m3], containing current volume. Both endpoints and
    /// all queried states must be supported by the original gas material.
    pub volume_bracket_m3: [f64; 2],
    /// Existing gas temperature, source-energy and iteration controls.
    pub gas: GasEnergySolve,
    /// Maximum absolute impulse residual [N s].
    pub momentum_tolerance_n_s: f64,
    /// Maximum absolute combined gas/kinetic/external-work residual [J].
    pub energy_tolerance_j: f64,
    /// Interior safeguarded Newton/bisection trials in 1..=128. In addition,
    /// at most three material proposals check current volume and the bracket.
    pub maximum_iterations: usize,
}

/// No error changes accepted gas, geometry or velocity.
#[derive(Clone, Debug, PartialEq)]
pub enum PistonError {
    /// Invalid physical data, controls or unrepresentable balances.
    Invalid(&'static str),
    /// The shared thermochemical/material owner refused a candidate.
    Gas(GasEnergyError),
    /// Cancellation observed before atomic publication.
    Cancelled,
    /// The physical impulse residual has no sign change within caller bounds.
    OutsideBracket { lower_impulse_n_s: f64, upper_impulse_n_s: f64 },
    /// Work budget or floating-point progress exhausted.
    NotConverged { iterations: usize, impulse_residual_n_s: f64 },
}
impl core::fmt::Display for PistonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed gas/piston coupling: {self:?}")
    }
}
impl std::error::Error for PistonError {}
impl From<GasEnergyError> for PistonError {
    fn from(error: GasEnergyError) -> Self {
        if error == GasEnergyError::Cancelled { Self::Cancelled } else { Self::Gas(error) }
    }
}

/// Actual accepted material and mechanical work, not a topology-only certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct PistonStep {
    /// Gas owner's actual heat/work/entropy and stage accounting.
    pub gas: GasHeatVolumeStep,
    /// Source-bearing midpoint used for force, work and acoustic observations.
    pub midpoint: ResolvedThermochemicalGasState,
    /// Signed wall displacement [m], positive on gas expansion.
    pub displacement_m: f64,
    /// Accepted initial and final wall velocities [m/s].
    pub velocity_m_s: [f64; 2],
    /// Actual change of 1/2 M v^2 [J].
    pub kinetic_energy_change_j: f64,
    /// Additional prescribed force times actual displacement [J].
    pub external_force_work_j: f64,
    /// Work delivered against the external pressure, pe*(V1-V0) [J].
    pub exterior_pressure_work_j: f64,
    /// M*(v1-v0) - dt*(A*(pm-pe)+F) [N s].
    pub momentum_residual_n_s: f64,
    /// Kinetic change minus gas work plus exterior work minus force work [J].
    pub mechanical_energy_residual_j: f64,
    /// Gas internal plus kinetic change plus exterior work minus force work
    /// minus prescribed heat [J], independently measured from accepted states.
    pub total_energy_residual_j: f64,
    /// Interior mechanical trials; current-volume/bracket checks are separate.
    pub iterations: usize,
    /// Actual material proposals, each bounded by its own gas controls.
    pub material_trials: usize,
}

struct Trial {
    gas: GasHeatVolumeCandidate,
    displacement: f64,
    velocity: f64,
    residual: f64,
    tangent: f64,
}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn invalid(s: &'static str) -> PistonError { PistonError::Invalid(s) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PistonError> {
    if cancelled() { Err(PistonError::Cancelled) } else { Ok(()) }
}

impl GasPiston {
    /// Couple a resolved gas parcel to explicit inertial mass and constant area.
    /// Positive velocity expands the cavity. No material or atmosphere preset.
    pub fn new(gas: ClosedGasParcel, mass_kg: f64, area_m2: f64, velocity_m_s: f64)
        -> Result<Self, PistonError>
    {
        if !positive(mass_kg) || !positive(area_m2) || !velocity_m_s.is_finite()
            || !positive(gas.volume_m3()/area_m2)
            || !(0.5*mass_kg*velocity_m_s*velocity_m_s).is_finite()
        { return Err(invalid("positive represented mass/area/position and finite velocity/kinetic energy required")); }
        Ok(Self { gas, mass_kg, area_m2, velocity_m_s })
    }
    /// Accepted common gas state, geometry and sources; no independent pressure.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    /// Conserved inertial mass [kg], not the gas mass.
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Constant work-conjugate piston area [m2].
    #[must_use]
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Current accepted wall velocity [m/s].
    #[must_use]
    pub const fn velocity_m_s(&self) -> f64 { self.velocity_m_s }
    /// Effective wall position V/A [m], with no second stored geometry state.
    #[must_use]
    pub fn position_m(&self) -> f64 { self.gas.volume_m3()/self.area_m2 }

    fn trial(&self, volume: f64, input: PistonInput, controls: PistonSolve,
        cancelled: &mut impl FnMut() -> bool) -> Result<Trial, PistonError>
    {
        poll(cancelled)?;
        let gas = self.gas.midpoint_heat_volume(volume, input.heat_j, controls.gas, &mut *cancelled)?;
        let displacement = (volume-self.gas.volume_m3())/self.area_m2;
        let velocity = 2.0*(displacement/input.duration_s)-self.velocity_m_s;
        let force = self.area_m2*(gas.midpoint().state().pressure-input.exterior_pressure_pa)
            + input.external_force_n;
        let residual = self.mass_kg*(velocity-self.velocity_m_s)-input.duration_s*force;
        let tangent = (2.0*(self.mass_kg/self.area_m2))/input.duration_s
            - (input.duration_s*self.area_m2)*gas.step().pressure_final_volume_derivative_pa_m3;
        if !displacement.is_finite() || !velocity.is_finite() || !force.is_finite()
            || !residual.is_finite() || !positive(tangent) || !positive(volume/self.area_m2)
            || !(0.5*self.mass_kg*velocity*velocity).is_finite()
        { return Err(invalid("piston candidate position, velocity, impulse or tangent is unrepresentable")); }
        Ok(Trial { gas, displacement, velocity, residual, tangent })
    }

    fn finish(&mut self, trial: Trial, input: PistonInput, controls: PistonSolve,
        iterations: usize, material_trials: usize, cancelled: &mut impl FnMut() -> bool)
        -> Result<PistonStep, PistonError>
    {
        let before = self.velocity_m_s;
        let kinetic = (0.5*self.mass_kg*(trial.velocity-before))*(trial.velocity+before);
        let exterior = input.exterior_pressure_pa*(trial.gas.endpoint().volume_m3()-self.gas.volume_m3());
        let external = input.external_force_n*trial.displacement;
        let gas_step = trial.gas.step();
        let mechanical = (kinetic-gas_step.work_by_gas_j)+exterior-external;
        let total = (gas_step.stored_energy_change_j+kinetic)+exterior-external-input.heat_j;
        if ![kinetic, exterior, external, mechanical, total].iter().all(|v| v.is_finite())
            || trial.residual.abs() > controls.momentum_tolerance_n_s
            || total.abs() > controls.energy_tolerance_j
        { return Err(invalid("accepted piston/gas impulse or independently measured total energy did not close")); }
        let midpoint = trial.gas.midpoint().clone();
        let (gas, gas_step) = trial.gas.into_parts();
        let report = PistonStep { gas: gas_step, midpoint, displacement_m: trial.displacement,
            velocity_m_s: [before, trial.velocity], kinetic_energy_change_j: kinetic,
            external_force_work_j: external, exterior_pressure_work_j: exterior,
            momentum_residual_n_s: trial.residual, mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total, iterations, material_trials };
        poll(cancelled)?;
        self.gas = gas;
        self.velocity_m_s = trial.velocity;
        Ok(report)
    }

    /// Solve gas energy and piston momentum together; publish only after both
    /// the physical impulse and independent total-energy residual pass.
    ///
    /// V1-V0 = A*dt*(v0+v1)/2 and M*(v1-v0) = dt*(A*(pm-pe)+F).
    /// The gas computes pm and its exact DISCRETE volume tangent. Every trial
    /// starts from the same accepted gas state; rejected trials never advance it.
    /// Exterior pressure and additional force are constant during this call.
    ///
    /// Caller cancellation can be bound to the existing runtime CancelGate.
    /// At most maximum_iterations+3 gas proposals are made, each with the gas
    /// owner's bounded source work. Error/cancellation leaves all state intact.
    /// This is uniform closed-gas dynamics, not piston contact, valves, friction,
    /// combustion, shocks or an exact-isentropy/nonnegative-entropy guarantee.
    pub fn step(&mut self, input: PistonInput, controls: PistonSolve,
        mut cancelled: impl FnMut() -> bool) -> Result<PistonStep, PistonError>
    {
        poll(&mut cancelled)?;
        let [mut lo, mut hi] = controls.volume_bracket_m3;
        if !positive(input.duration_s) || !input.heat_j.is_finite()
            || !input.exterior_pressure_pa.is_finite() || input.exterior_pressure_pa < 0.0
            || !input.external_force_n.is_finite() || !positive(lo) || !positive(hi) || lo >= hi
            || self.gas.volume_m3() < lo || self.gas.volume_m3() > hi
            || !positive(controls.momentum_tolerance_n_s) || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
        { return Err(invalid("invalid piston drive, volume bracket, residual tolerances or iteration cap")); }
        let first = self.trial(self.gas.volume_m3(), input, controls, &mut cancelled)?;
        if first.residual.abs() <= controls.momentum_tolerance_n_s {
            return self.finish(first, input, controls, 0, 1, &mut cancelled);
        }
        let lower = self.trial(lo, input, controls, &mut cancelled)?;
        let upper = self.trial(hi, input, controls, &mut cancelled)?;
        if lower.residual.abs() <= controls.momentum_tolerance_n_s {
            return self.finish(lower, input, controls, 0, 3, &mut cancelled);
        }
        if upper.residual.abs() <= controls.momentum_tolerance_n_s {
            return self.finish(upper, input, controls, 0, 3, &mut cancelled);
        }
        if lower.residual > 0.0 || upper.residual < 0.0 || lower.residual > upper.residual {
            return Err(PistonError::OutsideBracket { lower_impulse_n_s: lower.residual,
                upper_impulse_n_s: upper.residual });
        }
        let mut volume = self.gas.volume_m3();
        let proposed = volume-first.residual/first.tangent;
        if proposed > lo && proposed < hi { volume = proposed; }
        else { volume = lo+0.5*(hi-lo); }
        let mut residual = first.residual;
        for iteration in 1..=controls.maximum_iterations {
            let trial = self.trial(volume, input, controls, &mut cancelled)?;
            residual = trial.residual;
            if residual.abs() <= controls.momentum_tolerance_n_s {
                return self.finish(trial, input, controls, iteration, iteration+3, &mut cancelled);
            }
            if residual < 0.0 { lo = volume; } else { hi = volume; }
            let proposed = volume-residual/trial.tangent;
            let margin = 0.1*(hi-lo);
            let next = if proposed > lo+margin && proposed < hi-margin { proposed }
                else { lo+0.5*(hi-lo) };
            if next <= lo || next >= hi || next == volume {
                return Err(PistonError::NotConverged { iterations: iteration, impulse_residual_n_s: residual });
            }
            volume = next;
        }
        Err(PistonError::NotConverged { iterations: controls.maximum_iterations, impulse_residual_n_s: residual })
    }
}

#[cfg(test)]
mod tests;
