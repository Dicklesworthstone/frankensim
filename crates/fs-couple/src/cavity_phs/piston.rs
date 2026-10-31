//! A free or fixed inertial wall coupled to a source-backed gas cavity.
//!
//! The gas owns caloric state, EOS, source support, transport and acoustic
//! properties. This module owns ONE translational boundary coordinate and its
//! momentum or fixed-support reaction. Both use the same midpoint pressure: there is no prescribed-motion
//! animation, nominal-air load, or independent pressure/work approximation.

/// Two source-resolved gases acting on one inertial wall and one momentum state.
pub mod double_acting;
mod impulse;
mod support;
pub use support::PistonSupport;
mod constitutive;
pub use constitutive::{PistonGas, PistonGasState};
/// Conservative nozzle exchange with a finite reservoir and the same wall step.
pub mod exchange;

/// Shared source-resolved rigid-body inertia coupled to the same gas dynamics.
pub mod material;
/// Source-bound reaction progress and finite gas energy driving the same wall.
pub mod reactive;
/// Reciprocal finite-body enthalpy/phase exchange on the same gas/momentum step.
pub mod thermal;

use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::midpoint::{
    GasHeatVolumeCandidate, GasHeatVolumeStep,
};

/// Constant-area cavity with a translating inertial end wall.
/// `step` conserves chamber mass; `step_with_nozzle` admits conservative exchange
/// with a finite reservoir before advancing the same gas and wall equations.
/// Volume is the single position owner: the effective wall coordinate is V/A.
/// Other rigid-body coordinates, contact stops and crank linkages are not modeled.
/// The default is the existing ideal thermochemical parcel. A source-backed
/// `IsochoricRealGasParcel` uses the SAME mechanics and acceptance with its
/// nonideal midpoint pressure, calorics and acoustic compressibility. Existing
/// nozzle and reactive adapters retain their ideal-parcel domain; finite-body
/// thermal exchange also accepts this same real-gas owner.
#[derive(Clone, Debug)]
pub struct GasPiston<G: PistonGas = ClosedGasParcel> {
    gas: G,
    mass_kg: f64,
    area_m2: f64,
    velocity_m_s: f64,
    support: PistonSupport,
    ports: crate::port_adoption::ModulePorts,
}

// Equality remains equality of accepted physical state. Observing a stationary
// step may add an account without changing the gas, geometry or momentum.
impl<G: PistonGas> PartialEq for GasPiston<G> {
    fn eq(&self, other: &Self) -> bool {
        self.gas == other.gas
            && self.mass_kg == other.mass_kg
            && self.area_m2 == other.area_m2
            && self.velocity_m_s == other.velocity_m_s
            && self.support == other.support
    }
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
    /// Fixed support accepts [V, V] and queries only the unchanged volume.
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
    /// Source-backed gas compressibility is unresolved at this wall clock.
    TimeResolution(fs_phs::linear_exact::ResolvedBandError),
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
        write!(f, "gas/piston coupling: {self:?}")
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
pub struct PistonStep<S = ResolvedThermochemicalGasState> {
    /// Gas owner's actual heat/work/entropy and stage accounting.
    pub gas: GasHeatVolumeStep,
    /// Source-bearing midpoint used for force, work and acoustic observations.
    pub midpoint: S,
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
    /// Signed impulse exerted ON the wall by its fixed support [N s].
    /// Zero for a freely translating wall; positive towards gas expansion.
    pub support_impulse_n_s: f64,
    /// Support impulse times actual midpoint wall velocity [J]. Both available
    /// modes give exactly zero: free has no support, fixed has no motion.
    pub support_work_j: f64,
    /// M*(v1-v0) - dt*(A*(pm-pe)+F) - support impulse [N s].
    pub momentum_residual_n_s: f64,
    /// Kinetic change minus gas work plus exterior work minus force/support work [J].
    pub mechanical_energy_residual_j: f64,
    /// Gas internal plus kinetic change plus exterior work minus force work
    /// minus support work and prescribed heat [J], measured from accepted states.
    pub total_energy_residual_j: f64,
    /// Interior mechanical trials; current-volume/bracket checks are separate.
    pub iterations: usize,
    /// Actual material proposals, each bounded by its own gas controls.
    pub material_trials: usize,
}

struct Trial<G: PistonGas = ClosedGasParcel> {
    gas: GasHeatVolumeCandidate<G, G::State>,
    displacement: f64,
    velocity: f64,
    residual: f64,
    tangent: f64,
    support_impulse: f64,
}

// Shared native kinematics and impulse account. The gas constitutive owner
// chooses the pressure sample; neither chemistry nor heat owns a second wall.
struct MechanicalTrial {
    displacement: f64,
    velocity: f64,
    residual: f64,
    support_impulse: f64,
}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn invalid(s: &'static str) -> PistonError { PistonError::Invalid(s) }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PistonError> {
    if cancelled() { Err(PistonError::Cancelled) } else { Ok(()) }
}

impl<G: PistonGas> GasPiston<G> {
    /// Couple a resolved gas parcel to explicit inertial mass and constant area.
    /// Positive velocity expands the cavity. No material or atmosphere preset.
    pub fn new(
        gas: G,
        mass_kg: f64,
        area_m2: f64,
        velocity_m_s: f64,
    ) -> Result<Self, PistonError> {
        if !positive(mass_kg) || !positive(area_m2) || !velocity_m_s.is_finite()
            || !positive(gas.volume_m3()/area_m2)
            || !(0.5*mass_kg*velocity_m_s*velocity_m_s).is_finite()
        { return Err(invalid("positive represented mass/area/position and finite velocity/kinetic energy required")); }
        let ports = crate::port_adoption::ModulePorts::new(
            "gas-piston",
            vec![
                ("wall", fs_ports::PortKind::MechanicalForceVelocity),
                (
                    "heat/primary",
                    fs_ports::PortKind::ThermalTemperatureEntropy,
                ),
            ],
        );
        Ok(Self {
            gas,
            mass_kg,
            area_m2,
            velocity_m_s,
            support: PistonSupport::Free,
            ports,
        })
    }
    /// Accepted common gas state, geometry and sources; no independent pressure.
    #[must_use]
    pub const fn gas(&self) -> &G { &self.gas }
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

    fn admit_linear_clock(&self, state: &G::State, volume_m3: f64,
        duration_s: f64) -> Result<(),PistonError>
    {
        // A fixed constraint removes this wall-compressibility oscillator.
        // Gas caloric, reaction and transport owners retain their own controls.
        if self.support == PistonSupport::Fixed { return Ok(()); }
        let density=state.density_kg_m3();
        let sound_speed=state.sound_speed_m_s();
        let omega=(density*sound_speed*sound_speed*self.area_m2*self.area_m2
            /(self.mass_kg*volume_m3)).sqrt();
        fs_phs::linear_exact::admit_midpoint_band(omega,duration_s)
            .map_err(PistonError::TimeResolution)
    }

    fn trial(&self, volume: f64, input: PistonInput, controls: PistonSolve,
        cancelled: &mut impl FnMut() -> bool) -> Result<Trial<G>, PistonError>
    {
        poll(cancelled)?;
        let gas = self.gas.midpoint_heat_volume(volume, input.heat_j, controls.gas, &mut *cancelled)?;
        let motion = self.mechanical_trial(volume, gas.midpoint().pressure_pa(), input)?;
        let tangent = (2.0*(self.mass_kg/self.area_m2))/input.duration_s
            - (input.duration_s*self.area_m2)*gas.step().pressure_final_volume_derivative_pa_m3;
        if !positive(tangent) {
            return Err(invalid("piston candidate tangent is unrepresentable"));
        }
        Ok(Trial { gas, displacement: motion.displacement, velocity: motion.velocity,
            residual: motion.residual, tangent, support_impulse: motion.support_impulse })
    }

    fn mechanical_trial(&self, volume: f64, pressure: f64, input: PistonInput)
        -> Result<MechanicalTrial, PistonError>
    {
        let displacement = (volume-self.gas.volume_m3())/self.area_m2;
        let velocity = 2.0*(displacement/input.duration_s)-self.velocity_m_s;
        let force = self.area_m2*(pressure-input.exterior_pressure_pa)
            + input.external_force_n;
        let support_impulse = if self.support == PistonSupport::Fixed {
            if self.velocity_m_s != 0.0 || volume != self.gas.volume_m3() {
                return Err(invalid("fixed support requires its original volume and zero velocity"));
            }
            -input.duration_s*force
        } else { 0.0 };
        let residual = self.mass_kg*(velocity-self.velocity_m_s)-input.duration_s*force-support_impulse;
        if !positive(pressure) || !displacement.is_finite() || !velocity.is_finite() || !force.is_finite()
            || !support_impulse.is_finite() || !residual.is_finite() || !positive(volume/self.area_m2)
            || !(0.5*self.mass_kg*velocity*velocity).is_finite()
        { return Err(invalid("piston candidate position, velocity or impulse is unrepresentable")); }
        Ok(MechanicalTrial { displacement, velocity, residual, support_impulse })
    }

    fn finish(
        &mut self,
        trial: Trial<G>,
        input: PistonInput,
        controls: PistonSolve,
        iterations: usize,
        material_trials: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PistonStep<G::State>, PistonError> {
        self.admit_linear_clock(trial.gas.midpoint(),
            f64::midpoint(self.gas.volume_m3(),trial.gas.endpoint().volume_m3()),input.duration_s)?;
        self.admit_linear_clock(trial.gas.endpoint().state(),trial.gas.endpoint().volume_m3(),input.duration_s)?;
        let before = self.velocity_m_s;
        let kinetic = (0.5*self.mass_kg*(trial.velocity-before))*(trial.velocity+before);
        let exterior = input.exterior_pressure_pa*(trial.gas.endpoint().volume_m3()-self.gas.volume_m3());
        let external = input.external_force_n*trial.displacement;
        let support_work = trial.support_impulse * f64::midpoint(before, trial.velocity);
        let force = self.area_m2*(trial.gas.midpoint().pressure_pa()-input.exterior_pressure_pa)
            + input.external_force_n;
        let momentum = self.mass_kg*(trial.velocity-before)-input.duration_s*force-trial.support_impulse;
        let gas_step = trial.gas.step();
        let mechanical = (kinetic-gas_step.work_by_gas_j)+exterior-external-support_work;
        let total = (gas_step.stored_energy_change_j+kinetic)+exterior-external-support_work-input.heat_j;
        if ![kinetic, exterior, external, support_work, momentum, mechanical, total].iter().all(|v| v.is_finite())
            || momentum.abs() > controls.momentum_tolerance_n_s
            || total.abs() > controls.energy_tolerance_j
        { return Err(invalid("accepted piston/gas impulse or independently measured total energy did not close")); }
        let midpoint = trial.gas.midpoint().clone();
        let (gas, gas_step) = trial.gas.into_parts();
        let report = PistonStep { gas: gas_step, midpoint, displacement_m: trial.displacement,
            velocity_m_s: [before, trial.velocity], kinetic_energy_change_j: kinetic,
            external_force_work_j: external, exterior_pressure_work_j: exterior,
            support_impulse_n_s: trial.support_impulse, support_work_j: support_work,
            momentum_residual_n_s: momentum, mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total, iterations, material_trials };
        poll(cancelled)?;
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: self.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(
                    report.gas.stored_energy_change_j + report.kinetic_energy_change_j,
                ),
                loss_j: Some(0.0),
                included_work_j: 0.0,
                entropy_change_j_k: report.gas.entropy_change_j_k,
                tolerance_j: Some(controls.energy_tolerance_j),
                unavailable: None,
            },
            [
                crate::port_adoption::PortTransfer::work(external - exterior + support_work, true),
                crate::port_adoption::PortTransfer::heat(
                    input.heat_j,
                    report.midpoint.temperature_k(),
                ),
            ],
        );
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
    /// With fixed support, V and velocity remain exactly unchanged/zero; the
    /// support supplies the balancing impulse with zero mechanical work. The
    /// same gas owner still resolves all heat, material and pressure changes.
    ///
    /// Caller cancellation can be bound to the existing runtime CancelGate.
    /// At most maximum_iterations+3 gas proposals are made, each with the gas
    /// owner's bounded source work. Error/cancellation leaves all state intact.
    /// This is uniform closed-gas dynamics, not piston contact, valves, friction,
    /// combustion, shocks or an exact-isentropy/nonnegative-entropy guarantee.
    pub fn step(&mut self, input: PistonInput, controls: PistonSolve,
        mut cancelled: impl FnMut() -> bool) -> Result<PistonStep<G::State>, PistonError>
    {
        poll(&mut cancelled)?;
        let [lo, hi] = controls.volume_bracket_m3;
        if !positive(input.duration_s) || !input.heat_j.is_finite()
            || !input.exterior_pressure_pa.is_finite() || input.exterior_pressure_pa < 0.0
            || !input.external_force_n.is_finite() || !positive(lo) || !positive(hi) || lo > hi
            || (lo == hi && self.support == PistonSupport::Free)
            || self.gas.volume_m3() < lo || self.gas.volume_m3() > hi
            || !positive(controls.momentum_tolerance_n_s) || !positive(controls.energy_tolerance_j)
            || !(1..=128).contains(&controls.maximum_iterations)
        { return Err(invalid("invalid piston drive, volume bracket, residual tolerances or iteration cap")); }
        if self.support == PistonSupport::Fixed {
            let trial = self.fixed_trial(input, controls, &mut cancelled)?;
            return self.finish(trial, input, controls, 0, 1, &mut cancelled);
        }
        self.admit_linear_clock(self.gas.state(),self.gas.volume_m3(),input.duration_s)?;
        let solved = impulse::solve(self.gas.volume_m3(), controls, |volume| {
            let trial = self.trial(volume, input, controls, &mut cancelled)?;
            Ok(impulse::Trial { residual: trial.residual, tangent: trial.tangent, value: trial })
        })?;
        self.finish(solved.value, input, controls, solved.iterations, solved.proposals, &mut cancelled)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod warping;

// Both source-material families expose the same original accepted port account.
impl<G: PistonGas> fs_ports::PortedModule for GasPiston<G> {
    fn port_schemas(&self) -> &[fs_ports::PortSchema] { &self.ports.schemas }
    fn port_coordinate_charts(&self) -> &[fs_ports::PortCoordinateChart] { &self.ports.charts }
    fn junction_schemas(&self) -> &[fs_ports::JunctionSchema] { &self.ports.junctions }
    fn last_window_duration_s(&self) -> Option<f64> { self.ports.duration_s() }
    fn last_window_transfers(&self) -> Option<&[fs_ports::IntegratedWindowTransfer]> {
        self.ports.transfers()
    }
    fn last_window_transfer_error(&self) -> Option<&fs_ports::CoupleError> {
        self.ports.transfer_error()
    }
    fn audit_last_window(&self) -> Result<fs_ports::ModuleWindowAudit, fs_ports::CoupleError> {
        self.ports.audit()
    }
}

#[cfg(test)]
mod nonideal_tests;

