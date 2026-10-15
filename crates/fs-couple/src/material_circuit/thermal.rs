//! Finite circuit storage coupled to the existing body/gas/piston transaction.
//!
//! The electrical source is NOT drained during trial evaluations. Its accepted
//! charge/flux state is held until the same midpoint solves electrical loss,
//! material heating, phase storage, heat exchange and mechanical pressure work.

use std::cell::RefCell;
use super::{MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitError,
    MaterialCircuitSolve, MaterialCircuitStep, MaterialResistor};
use crate::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve,
    ThermalPistonCandidate, ThermalPistonError, ThermalPistonStep};
use crate::cavity_phs::piston::thermal::power::{PowerSample, PowerFeedbackStep, PoweredThermalCandidate};
use crate::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

/// One uniform thermal body receiving ALL resistive heat from a circuit. The
/// circuit can contain multiple C/L/source/transformer branches, but exactly one
/// material resistor. Multi-body circuits use the general per-branch candidate;
/// their other resistor losses must never be silently assigned to this body.
pub struct CircuitThermalSource<'a> {
    /// Accepted circuit storage; borrowed immutably through all trials.
    pub circuit: &'a MaterialCircuit,
    /// Original graph branch index of this whole-body resistor.
    pub resistor_branch: usize,
    /// Same enthalpy-card binding used by ideal and circuit-driven Ohmic heating.
    pub body: MaterialOhmicBody<'a>,
    /// Existing circuit external ports held during this step. Empty for a closed
    /// finite supply. Source work is retained separately from stored energy.
    pub inputs: &'a [f64],
}

/// Electrical work and combined physical energy admission, independent of the
/// existing source, heat-exchange and piston iteration/tolerance controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircuitThermalSolve {
    /// Original branch-power and electrical-storage checks.
    pub circuit: MaterialCircuitSolve,
    /// Maximum absolute residual over electrical + thermal + gas + kinetic [J].
    pub total_energy_tolerance_j: f64,
    /// Maximum actual circuit source evaluations, in 1..=130. Each has the
    /// bounded graph/state work of `MaterialCircuit::step_candidate`.
    pub maximum_circuit_trials: usize,
}

/// Original domain-specific failures remain inspectable.
#[derive(Debug)]
pub enum CircuitHeatingError<E> {
    /// Original thermal interface callback refusal.
    Interface(E),
    /// Original same-card conductor/phase/geometry refusal.
    Conductor(JouleHeatingError),
    /// Original circuit or material-heat-partition refusal.
    Circuit(MaterialCircuitError),
    /// Requested global circuit evaluation budget exhausted.
    Budget,
    /// Incomplete physical binding or unrepresentable integrated power.
    Invalid(&'static str),
}
impl<E: core::fmt::Debug> core::fmt::Display for CircuitHeatingError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"circuit thermal source: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CircuitHeatingError<E> {}

/// Three physical owners accepted together; energy transfer is already included.
pub struct CircuitThermalCandidate<I> {
    powered: PoweredThermalCandidate<I, MaterialCircuitCandidate>,
    total_energy_residual_j: f64,
    circuit_trials: usize,
}
/// Original electrical, thermomechanical and nonlinear-source reports.
pub struct CircuitThermalStep {
    /// Circuit storage, external supply work and source-bound resistor heat.
    pub electrical: MaterialCircuitStep,
    /// Gas, body and mechanical energy/momentum accounting.
    pub thermal: ThermalPistonStep,
    /// Applied versus evaluated electrical heat and original work counts.
    pub power: PowerFeedbackStep,
    /// Independently evaluated combined storage/supply balance [J].
    pub total_energy_residual_j: f64,
    /// All circuit proposals, not only the final accepted one.
    pub circuit_trials: usize,
}
impl<I> CircuitThermalCandidate<I> {
    /// Same gas, body and momentum candidate produced by the existing owner.
    pub fn coupled(&self) -> &ThermalPistonCandidate<I> { self.powered.coupled() }
    /// Corresponding tentative charge/flux state and electrical source receipts.
    pub fn electrical(&self) -> &MaterialCircuitCandidate { &self.powered.power().source }
    /// Combined energy balance, including finite supply depletion [J].
    pub const fn total_energy_residual_j(&self) -> f64 { self.total_energy_residual_j }
    /// Move ALL accepted endpoints together. Do not apply either Joule heat or
    /// interface heat again; those transfers are already in the returned states.
    pub fn into_parts(self) -> (GasPiston, UniformEnthalpyStep<I>, MaterialCircuit, CircuitThermalStep) {
        let (coupled,supply,power)=self.powered.into_parts();
        let (piston,body,thermal)=coupled.into_parts();
        let (circuit,electrical)=supply.source.into_parts();
        (piston,body,circuit,CircuitThermalStep { electrical,thermal,power,
            total_energy_residual_j:self.total_energy_residual_j,circuit_trials:self.circuit_trials })
    }
}

impl GasPiston {
    /// Solve finite electrical storage, state-dependent Joule loss, body phase,
    /// heat transfer and gas momentum using the EXISTING owners. The resistor is
    /// queried at the same body energy midpoint used by heat exchange. Every
    /// circuit trial starts with identical accepted charge/flux; no callback
    /// mutates or debits the supply. Capacitor/inductor energy is not an extra
    /// prescribed heat source, and external electrical port work is counted once.
    ///
    /// C/L, topology, conductor geometry and the thermal body's geometry are
    /// fixed. This supports dynamic circuit storage, not a changing magnetic or
    /// dielectric material law, switches/arcs, mass transport or resolved thermal
    /// gradients. Source-law support and timestep refinement remain required.
    /// One whole-body resistor is explicit; unmodeled resistive sinks refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_thermal_candidate<I,E>(
        &self, input:PistonInput, body:UniformEnthalpyStepInput<'_>,
        piston:PistonSolve, exchange:HeatExchangeSolve, source:HeatExchangeSolve,
        supply:CircuitThermalSource<'_>, electrical:CircuitThermalSolve,
        mut conductance:impl FnMut(&ResolvedThermochemicalGasState,EquilibriumPhaseState)
            -> Result<ConductanceSample<I>,E>,
        cancelled:impl FnMut()->bool,
    ) -> Result<CircuitThermalCandidate<I>,ThermalPistonError<CircuitHeatingError<E>>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        if !electrical.total_energy_tolerance_j.is_finite() || electrical.total_energy_tolerance_j<=0.0
            || !(1..=130).contains(&electrical.maximum_circuit_trials)
            || supply.circuit.resistors().len()!=1
            || supply.circuit.resistors()[0].branch!=supply.resistor_branch {
            return Err(ThermalPistonError::Invalid("one complete thermal resistor binding and bounded circuit work required"));
        }
        let mut trials=0usize;
        let proposed=self.powered_thermal_candidate(input,body,piston,exchange,source,
            |gas,phase| conductance(gas,phase).map_err(CircuitHeatingError::Interface),
            |_,phase| {
                if trials==electrical.maximum_circuit_trials { return Err(CircuitHeatingError::Budget); }
                trials+=1;
                let conductor=supply.body.resolve_conductor(phase,body.volume_m3)
                    .map_err(CircuitHeatingError::Conductor)?;
                let candidate=supply.circuit.step_candidate(
                    vec![MaterialResistor {branch:supply.resistor_branch,conductor}],
                    supply.inputs,input.duration_s,electrical.circuit,|| (*cancel.borrow_mut())())
                    .map_err(CircuitHeatingError::Circuit)?;
                let heat=candidate.step().heat_j;
                let power_w=heat/input.duration_s;
                if !power_w.is_finite() || (heat!=0.0 && power_w==0.0) {
                    return Err(CircuitHeatingError::Invalid("integrated circuit heat rate is unrepresentable"));
                }
                Ok(PowerSample {power_w,source:candidate})
            },|| (*cancel.borrow_mut())());
        let powered=match proposed {
            Err(ThermalPistonError::Interface(CircuitHeatingError::Circuit(MaterialCircuitError::Cancelled))) =>
                return Err(ThermalPistonError::Cancelled),
            other => other?,
        };
        let thermal=powered.coupled().step();
        let mechanics=&thermal.piston;
        let circuit=powered.power().source.step();
        let total=((circuit.stored_after_j-circuit.stored_before_j)
            +mechanics.gas.stored_energy_change_j+mechanics.kinetic_energy_change_j
            +thermal.thermal_energy_change_j)+mechanics.exterior_pressure_work_j
            -mechanics.external_force_work_j-input.heat_j-body.internal_heat_j-circuit.native.supplied;
        if !total.is_finite() || total.abs()>electrical.total_energy_tolerance_j {
            return Err(ThermalPistonError::Invalid("finite circuit and all thermal/mechanical storage do not close"));
        }
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        Ok(CircuitThermalCandidate {powered,total_energy_residual_j:total,circuit_trials:trials})
    }
}
