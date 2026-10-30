//! A thermal gradient can power a material field actuator, and prescribed
//! mechanical work can return through the same circuit and thermal stores.
//! Seebeck/Peltier/Thomson transfers use the original thermoelectric owner;
//! force and changing C/L use the original reciprocal field-motion descriptor.
//! No electrical, material or thermal history is advanced outside that shared
//! candidate. This is the original endpoint thermal discretization, not a new
//! time integrator or a temperature-dependent magnetic constitutive law.

use super::*;
use crate::circuit_enthalpy::motion::FieldMotionEnthalpyReport;
use crate::material_circuit::motion::{FieldMotionCandidate, FieldMotionInput,
    MaterialFieldMotion};

impl ThermoelectricDriver for FieldMotionCandidate {
    fn native(&self) -> &StepRecord { &self.step().native }
}

/// All native transfers plus the complete external-supply balance. Internal
/// Seebeck work and its opposite thermal transfer appear in `coupled` exactly
/// once and are removed together from the external supply accounts below.
pub struct ThermoelectricMotionReport<S = MaterialNetworkSample> {
    /// The original field/kinetic/circuit and finite enthalpy accounts.
    pub coupled: FieldMotionEnthalpyReport<ThermoelectricSample<S>>,
    /// Original electrical source work, excluding the internal Seebeck ports.
    pub external_electrical_work_j: f64,
    /// Held external force times the ACTUAL accepted gap displacement [J].
    pub external_mechanical_work_j: f64,
    /// Thermal boundary/source work, excluding Joule and Kelvin transfers [J].
    pub external_heat_j: f64,
    /// Signed numerical net of closed enclosure exchange [J]. Internal
    /// radiation remains in the joint balance and is never external supply.
    pub closed_radiation_net_j: f64,
    /// Field + mechanical + thermal storage minus all external supply [J].
    pub total_energy_residual_j: f64,
}

/// One publishable field/circuit, moving-body and thermal endpoint. Original
/// material sources, charge/flux, momentum and thermal event history survive.
pub struct ThermoelectricMotionCandidate<'m, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'m>,
    motion: MaterialFieldMotion,
    report: ThermoelectricMotionReport<S>,
}
impl<'m, S> ThermoelectricMotionCandidate<'m, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    pub fn circuit(&self) -> &MaterialCircuit { self.motion.circuit() }
    pub const fn report(&self) -> &ThermoelectricMotionReport<S> { &self.report }
    /// Move all accepted histories together, without replaying any physical law.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialFieldMotion,
        ThermoelectricMotionReport<S>) {
        (self.thermal, self.motion, self.report)
    }
}

/// Advance a sourced thermocouple/Peltier network and material field actuator.
///
/// Finite junction and bulk temperatures drive the original integrated S(T)
/// voltage. The original descriptor evolves that circuit's charge/flux, field
/// geometry and momentum together. Its ACTUAL source current determines the
/// reciprocal junction/interior heat; actual material Joule loss enters the
/// same finite enthalpy residual once. Every nonlinear trial starts from the
/// same accepted circuit and body, with freshly evaluated pinned transport.
///
/// `ports` follows the original circuit order and has zero placeholders for
/// internal Seebeck ports. External force is positive toward increasing gap.
/// Full balance includes field storage, body/support kinetic and elastic energy,
/// thermal storage, held
/// force work, remaining electrical ports and external heat before publication.
/// Refusal or cancellation returns no endpoint and preserves every input.
///
/// The existing thermoelectric solid-leg and moving-gap scope still applies:
/// no electrical contacts, field entropy, moving-conductor Lorentz heat,
/// temperature-dependent permeability/permittivity, expansion, saturation or
/// hysteresis is inferred. Thermal cells describe the stationary transport
/// network; the guided body and field medium retain their original admitted
/// source snapshots and inertia. No second heat credit for field work is made.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
    motion: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
    legs: &[ThermoelectricLeg], electrical_bases: &[(usize, ElectricalPropertyBasis)],
    controls: CircuitEnthalpySolve, guess_j: &[f64], mut cancelled: impl FnMut() -> bool,
) -> Result<ThermoelectricMotionCandidate<'m>, Error> {
    let (thermal, motion, coupled) = coupled_candidate(material, thermal, target_s,
        motion.circuit(), ports, legs, electrical_bases, controls, guess_j,
        |sources, inputs, dt, cancelled| motion.step_candidate(sources,
            FieldMotionInput { electrical_ports: inputs, external_force_n, duration_s: dt },
            controls.circuit, cancelled), &mut cancelled)?;
    finish_candidate(thermal, motion, coupled, 0.0, controls, cancelled)
}

// Fixed transport and complete contact/layer/radiation assemblies publish through
// this SAME field/Kelvin/enthalpy account. Only their retained source payloads
// differ; source resolvers never advance mechanical or electrical history.
pub(crate) fn finish_candidate<'m, S>(
    thermal: DepositionSimulation<'m>, motion: MaterialFieldMotion,
    coupled: FieldMotionEnthalpyReport<ThermoelectricSample<S>>,
    closed_radiation_net_j: f64, controls: CircuitEnthalpySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ThermoelectricMotionCandidate<'m, S>, Error> {
    let closed_radiation_net_j = finite(closed_radiation_net_j)?;
    let field = &coupled.electrical;
    let external_electrical_work_j = finite(field.electrical_supply_j
        - coupled.source.reversible_electrical_work_j)?;
    let external_mechanical_work_j = finite(field.mechanical_supply_j)?;
    let external_heat_j = finite(coupled.non_electrical_heat_j
        - coupled.source.reversible_heat_j - closed_radiation_net_j)?;
    let total_energy_residual_j = finite(field.electrical_energy_j[1]-field.electrical_energy_j[0]
        + field.kinetic_energy_j[1]-field.kinetic_energy_j[0]
        + field.elastic_energy_j[1]-field.elastic_energy_j[0]
        + coupled.thermal.stored_energy_change_j
        - external_electrical_work_j - external_mechanical_work_j - external_heat_j)?;
    if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
        return Err(invalid("thermoelectric, field-motion and thermal storage do not close external supply"));
    }
    if cancelled() { return Err(CircuitEnthalpyError::Cancelled); }
    Ok(ThermoelectricMotionCandidate { thermal, motion, report: ThermoelectricMotionReport {
        coupled, external_electrical_work_j, external_mechanical_work_j, external_heat_j,
        closed_radiation_net_j, total_energy_residual_j } })
}
