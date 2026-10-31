//! Finite circuit storage coupled to the existing body/gas/piston transaction.
//!
//! The electrical source is NOT drained during trial evaluations. Its accepted
//! charge/flux state is held until the same midpoint solves electrical loss,
//! material heating, phase storage, heat exchange and mechanical pressure work.

use std::cell::RefCell;
use super::{MaterialCircuit, MaterialCircuitCandidate, MaterialCircuitError,
    MaterialCircuitSolve, MaterialCircuitStep, MaterialResistor};
use super::conductivity::{ConductionAxis, ConductorResponse, ConductivitySolve,
    MaterialConductivityCandidate, MaterialConductivityStep, MAX_CONDUCTIVITY_ITERATIONS};
use crate::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve,
    ThermalPistonCandidate, ThermalPistonError, ThermalPistonStep};
use crate::cavity_phs::piston::thermal::power::{PowerSample, PowerFeedbackStep, PoweredThermalCandidate};
use crate::cavity_phs::piston::thermal::power::joule::{JouleHeatingError, MaterialOhmicBody};
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::QuantitySpec;

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
    /// Maximum actual native circuit evaluations, including discarded electrical
    /// iterates. The static route admits 1..=130; the nonlinear conductivity
    /// route admits up to 130 * MAX_CONDUCTIVITY_ITERATIONS. Each evaluation has
    /// the bounded graph/state work of `MaterialCircuit::step_candidate`.
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
pub struct CircuitThermalCandidate<I, C = MaterialCircuitCandidate> {
    powered: PoweredThermalCandidate<I, C>,
    total_energy_residual_j: f64,
    circuit_trials: usize,
}
/// Original electrical, thermomechanical and nonlinear-source reports.
pub struct CircuitThermalStep<E = MaterialCircuitStep> {
    /// Circuit storage, external supply work and source-bound resistor heat.
    pub electrical: E,
    /// Gas, body and mechanical energy/momentum accounting.
    pub thermal: ThermalPistonStep,
    /// Applied versus evaluated electrical heat and original work counts.
    pub power: PowerFeedbackStep,
    /// Independently evaluated combined storage/supply balance [J].
    pub total_energy_residual_j: f64,
    /// All circuit proposals, not only the final accepted one.
    pub circuit_trials: usize,
}
impl<I, C> CircuitThermalCandidate<I, C> {
    /// Same gas, body and momentum candidate produced by the existing owner.
    pub fn coupled(&self) -> &ThermalPistonCandidate<I> { self.powered.coupled() }
    /// Corresponding tentative charge/flux state and electrical source receipts.
    pub fn electrical(&self) -> &C { &self.powered.power().source }
    /// Combined energy balance, including finite supply depletion [J].
    pub const fn total_energy_residual_j(&self) -> f64 { self.total_energy_residual_j }
}
impl<I> CircuitThermalCandidate<I> {
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

impl<I> CircuitThermalCandidate<I, MaterialConductivityCandidate> {
    /// Publish the same complete thermal/mechanical transaction and original
    /// nonlinear electrical candidate. Constitutive residuals remain separate
    /// from the native Joule heat; they are never a second thermal source.
    pub fn into_parts(self) -> (
        GasPiston, UniformEnthalpyStep<I>, MaterialCircuit,
        CircuitThermalStep<MaterialConductivityStep>,
    ) {
        let (coupled, supply, power) = self.powered.into_parts();
        let (piston, body, thermal) = coupled.into_parts();
        let (circuit, electrical) = supply.source.into_parts();
        (piston, body, circuit, CircuitThermalStep { electrical, thermal, power,
            total_energy_residual_j: self.total_energy_residual_j,
            circuit_trials: self.circuit_trials })
    }
}

fn admit<E>(supply: &CircuitThermalSource<'_>, electrical: CircuitThermalSolve,
    maximum_trials: usize) -> Result<(), ThermalPistonError<CircuitHeatingError<E>>> {
    if !electrical.total_energy_tolerance_j.is_finite() || electrical.total_energy_tolerance_j <= 0.0
        || !(1..=maximum_trials).contains(&electrical.maximum_circuit_trials)
        || supply.circuit.resistors().len() != 1
        || supply.circuit.resistors()[0].branch != supply.resistor_branch {
        return Err(ThermalPistonError::Invalid(
            "one complete thermal resistor binding and bounded circuit work required"));
    }
    Ok(())
}

fn total_energy<E>(thermal: &ThermalPistonStep, circuit: &MaterialCircuitStep,
    input: PistonInput, body: UniformEnthalpyStepInput<'_>, tolerance: f64)
    -> Result<f64, ThermalPistonError<CircuitHeatingError<E>>> {
    let mechanics = &thermal.piston;
    let total = (circuit.stored_energy_change_j()
        + mechanics.gas.stored_energy_change_j + mechanics.kinetic_energy_change_j
        + thermal.thermal_energy_change_j) + mechanics.exterior_pressure_work_j
        - mechanics.external_force_work_j - input.heat_j - body.internal_heat_j - circuit.native.supplied;
    if !total.is_finite() || total.abs() > tolerance {
        return Err(ThermalPistonError::Invalid(
            "finite circuit and all thermal/mechanical storage do not close"));
    }
    Ok(total)
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
        conductance:impl FnMut(&ResolvedThermochemicalGasState,EquilibriumPhaseState)
            -> Result<ConductanceSample<I>,E>,
        cancelled:impl FnMut()->bool,
    ) -> Result<CircuitThermalCandidate<I>,ThermalPistonError<CircuitHeatingError<E>>> {
        self.circuit_thermal_candidate_with_basis(input,body,piston,exchange,source,supply,electrical,
            ElectricalPropertyBasis::Resistivity,conductance,cancelled)
    }

    /// The same finite-circuit/body/gas/momentum transaction with an explicitly
    /// chosen conductivity or resistivity source. The body's original selection
    /// and pins remain authoritative at each shared thermal midpoint; no failed
    /// source is replaced by the other electrical representation.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_thermal_candidate_with_basis<I,E>(
        &self, input:PistonInput, body:UniformEnthalpyStepInput<'_>,
        piston:PistonSolve, exchange:HeatExchangeSolve, source:HeatExchangeSolve,
        supply:CircuitThermalSource<'_>, electrical:CircuitThermalSolve, basis:ElectricalPropertyBasis,
        mut conductance:impl FnMut(&ResolvedThermochemicalGasState,EquilibriumPhaseState)
            -> Result<ConductanceSample<I>,E>,
        cancelled:impl FnMut()->bool,
    ) -> Result<CircuitThermalCandidate<I>,ThermalPistonError<CircuitHeatingError<E>>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        admit(&supply, electrical, 130)?;
        let mut trials=0usize;
        let proposed=self.powered_thermal_candidate(input,body,piston,exchange,source,
            |gas,phase| conductance(gas,phase).map_err(CircuitHeatingError::Interface),
            |_,phase| {
                if trials==electrical.maximum_circuit_trials { return Err(CircuitHeatingError::Budget); }
                trials+=1;
                let conductor=supply.body.resolve_conductor_with_basis(phase,body.volume_m3,basis)
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
        let total=total_energy(powered.coupled().step(),powered.power().source.step(),
            input,body,electrical.total_energy_tolerance_j)?;
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        Ok(CircuitThermalCandidate {powered,total_energy_residual_j:total,circuit_trials:trials})
    }

    /// Close the original E/J-dependent conductor law inside the SAME finite
    /// circuit/body/gas/momentum feedback. Every electrical iterate uses the
    /// actual body enthalpy midpoint supplied by the thermal owner; the existing
    /// conductivity owner then resolves its declared electric-field/current-
    /// density axes from the native circuit's solved algebraic midpoint.
    ///
    /// The accepted circuit's complete source bundle, selected claims, geometry
    /// and conductivity/resistivity basis remain authoritative. `axes` must bind
    /// every electrical operating coordinate; none may stay at the seed value.
    /// This body supplies temperature and any declared phase coordinate. Other
    /// material conditions stay prescribed by its original context. The source
    /// selection policy cannot replace the accepted circuit's pinned claims.
    ///
    /// Both electrical constitutive residuals AND the existing thermal/source
    /// residuals must converge. Only native circuit Joule heat enters the body;
    /// all rejected electrical trials start from identical accepted charge/flux.
    /// The global circuit-work cap includes every discarded conductivity iterate.
    /// Cancellation, source refusal, work exhaustion or nonconvergence publishes
    /// no candidate and drains no electrical or thermal state.
    ///
    /// This retains the static route's fixed geometry/C/L, one-whole-resistor
    /// and uniform-body restrictions. Contacts require their own jointly solved
    /// interface law and are refused by the original conductivity owner. Nested
    /// convergence does not establish timestep accuracy or uniqueness.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_conductivity_thermal_candidate<I, E>(
        &self, input: PistonInput, body: UniformEnthalpyStepInput<'_>,
        piston: PistonSolve, exchange: HeatExchangeSolve, source: HeatExchangeSolve,
        supply: CircuitThermalSource<'_>, electrical: CircuitThermalSolve,
        axes: &[ConductionAxis<'_>], conductivity: ConductivitySolve,
        mut conductance: impl FnMut(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<I>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CircuitThermalCandidate<I, MaterialConductivityCandidate>,
        ThermalPistonError<CircuitHeatingError<E>>> {
        let cancel = RefCell::new(cancelled);
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        admit(&supply, electrical, 130 * MAX_CONDUCTIVITY_ITERATIONS)?;
        conductivity.validate().map_err(|error| {
            ThermalPistonError::Interface(CircuitHeatingError::Circuit(error))
        })?;
        if axes.len() > 2 {
            return Err(ThermalPistonError::Invalid(
                "at most one electric-field and one current-density axis required"));
        }
        let original = &supply.circuit.resistors()[0].conductor;
        let operating_axes: Vec<_> = axes.iter().map(|axis| {
            (axis.name, QuantitySpec::dimensional(axis.coordinate.dimensions()))
        }).collect();
        let responses = [ConductorResponse { branch: supply.resistor_branch,
            card: supply.body.material.material, axes }];
        let mut trials = 0usize;
        let proposed = self.powered_thermal_candidate(input, body, piston, exchange, source,
            |gas, phase| conductance(gas, phase).map_err(CircuitHeatingError::Interface),
            |_, phase| {
                let remaining = electrical.maximum_circuit_trials - trials;
                if remaining == 0 { return Err(CircuitHeatingError::Budget); }
                let point = supply.body.material.at(phase)
                    .map_err(JouleHeatingError::from).map_err(CircuitHeatingError::Conductor)?;
                let material = original.material()
                    .at_query_point(supply.body.material.material, &point)
                    .map_err(fs_material::conductor::ConductorError::MaterialState)
                    .map_err(JouleHeatingError::from).map_err(CircuitHeatingError::Conductor)?;
                let conductor = supply.body.conductor_from_material_with_axes(
                    phase, body.volume_m3, material, original.source_basis(), &operating_axes,
                ).map_err(CircuitHeatingError::Conductor)?;
                let candidate = supply.circuit.conductivity_with_materials_candidate(
                    vec![MaterialResistor { branch: supply.resistor_branch, conductor }],
                    &responses, supply.inputs, input.duration_s, electrical.circuit,
                    ConductivitySolve { maximum_iterations: conductivity.maximum_iterations.min(remaining),
                        ..conductivity }, || (*cancel.borrow_mut())(),
                ).map_err(|error| {
                    if remaining < conductivity.maximum_iterations
                        && matches!(error, MaterialCircuitError::ConductivityConvergence { .. }) {
                        CircuitHeatingError::Budget
                    } else { CircuitHeatingError::Circuit(error) }
                })?;
                trials += candidate.step().iterations;
                let heat = candidate.step().electrical.heat_j;
                let power_w = heat / input.duration_s;
                if !power_w.is_finite() || (heat != 0.0 && power_w == 0.0) {
                    return Err(CircuitHeatingError::Invalid("integrated circuit heat rate is unrepresentable"));
                }
                Ok(PowerSample { power_w, source: candidate })
            }, || (*cancel.borrow_mut())());
        let powered = match proposed {
            Err(ThermalPistonError::Interface(CircuitHeatingError::Circuit(MaterialCircuitError::Cancelled))) =>
                return Err(ThermalPistonError::Cancelled),
            other => other?,
        };
        let total = total_energy(powered.coupled().step(), &powered.power().source.step().electrical,
            input, body, electrical.total_energy_tolerance_j)?;
        if (*cancel.borrow_mut())() { return Err(ThermalPistonError::Cancelled); }
        Ok(CircuitThermalCandidate { powered, total_energy_residual_j: total, circuit_trials: trials })
    }
}
