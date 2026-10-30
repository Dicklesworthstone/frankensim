//! A material actuator's electrical, mechanical and thermal histories advance
//! together. The existing field-motion descriptor supplies reciprocal gap work;
//! the existing enthalpy residual supplies finite-body heating and transport.
//! Endpoint conductivity/resistivity and contact resistance feed back into the
//! SAME electrical trajectory that produces force and changes C/L geometry.

/// Joint temperature/field/current-dependent conduction and actuator motion.
pub mod conductivity;
/// Live source-derived expansion and elastic response of the actuator's member.
pub mod support;

use super::*;
use super::contacts::{ContactHeatCells, ContactThermalCoefficients, ContactThermalSample};
use crate::material_circuit::motion::{FieldMotionCandidate, FieldMotionInput,
    FieldMotionStep, MaterialFieldMotion};

/// One accepted actuator, held external drives and complete bulk heat routes.
#[derive(Clone, Copy)]
pub struct FieldMotionEnthalpyInput<'a> {
    /// Original field/circuit/body history, never advanced by a material callback.
    pub motion: &'a MaterialFieldMotion,
    /// Original circuit's electrical port order.
    pub ports: &'a [f64],
    /// Held external force [N], positive toward increasing gap.
    pub external_force_n: f64,
    /// Exactly one receiving thermal cell for each original resistor branch.
    pub destinations: &'a [ResistorCell],
}

/// Original field-motion work plus the shared local/global enthalpy account.
pub type FieldMotionEnthalpyReport<S> = CircuitEnthalpyReport<S, FieldMotionStep>;

/// One publishable thermal, electrical and moving-body endpoint.
pub struct FieldMotionEnthalpyCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    motion: MaterialFieldMotion,
    report: FieldMotionEnthalpyReport<S>,
}
impl<'m, S> FieldMotionEnthalpyCandidate<'m, S> {
    /// All receiving bodies and original thermal event history.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Joint charge/flux, gap, momentum and retained material snapshots.
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    /// Circuit view of the SAME accepted moving-gap history.
    pub fn circuit(&self) -> &MaterialCircuit { self.motion.circuit() }
    /// Original motion step, deposited losses and exact endpoint material sample.
    pub const fn report(&self) -> &FieldMotionEnthalpyReport<S> { &self.report }
    /// Publish all three histories and their work account together.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialFieldMotion, FieldMotionEnthalpyReport<S>) {
        (self.thermal, self.motion, self.report)
    }
}

impl ThermalDriver for FieldMotionCandidate {
    type Model = MaterialFieldMotion;
    type Step = FieldMotionStep;
    fn account(&self) -> DriverAccount<'_> {
        field_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) { FieldMotionCandidate::into_parts(self) }
}

fn field_account(step: &FieldMotionStep) -> DriverAccount<'_> {
    DriverAccount { resistors: &step.resistors, contacts: &step.contacts,
        heat_j: step.heat_j,
        numerical_dissipation_j: 0.0,
        stored_change_j: (step.electrical_energy_j[1]-step.electrical_energy_j[0])
            +(step.kinetic_energy_j[1]-step.kinetic_energy_j[0])
            +(step.elastic_energy_j[1]-step.elastic_energy_j[0]),
        supplied_j: step.electrical_supply_j+step.mechanical_supply_j }
}

/// Close material resistance, Joule heating and reciprocal gap motion in one
/// interval. Optional support elastic energy and consistent inertia stay in that
/// same driver balance; its original frozen material state is not heated here.
/// Use a prepared material network's `resolve` to retain source pins.
/// Every thermal trial restarts from the SAME accepted flux/charge, gap and
/// momentum, evaluates transport at its proposed enthalpy state, and advances
/// the original field-motion descriptor for the common physical duration.
///
/// Bulk/contact heat is deposited once. Electrical field energy, moving-body
/// kinetic energy, finite thermal storage, original electrical work, external
/// force times actual displacement and authored heat must close before return.
/// Material/convergence/domain/budget refusal and cancellation return no partial
/// state. Thermal events retain the existing left-limit and interval-cut rules.
///
/// Field coefficients and moving-body inertia keep the original motion model's
/// fixed, admitted material snapshots. This heats stationary winding/electrode
/// resistors and electrical interfaces; it does not add thermal expansion,
/// moving-conductor Joule laws, magnetic hysteresis or field entropy. Do not bind
/// a heated winding cell as an independently changing inertial mass. A thermal
/// coefficient law cannot silently substitute for missing reversible field heat.
/// Coupled temporal accuracy remains first order in general: thermal rates and
/// material transport use the endpoint, while the native motion step supplies
/// its own discrete-gradient work. No lagged-coefficient fallback is used.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: FieldMotionEnthalpyInput<'_>, controls: CircuitEnthalpySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldMotionEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    coupled_candidate(thermal, target_s, input, &[], controls, guess_j, |state| {
        Ok(ContactThermalCoefficients { bulk: resolve(state)?, contacts: Vec::new() })
    }, cancelled)
}

/// Original reciprocal field motion and finite heating with endpoint total
/// ambient/reservoir conductances. The thermal owner receives replacements,
/// not a second heat source; all original electrical and mechanical work enters
/// the same joint acceptance. Contact-bound circuits need their contact owner.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_boundary_conductances<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: FieldMotionEnthalpyInput<'_>, controls: CircuitEnthalpySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<BoundaryCircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldMotionEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let original = input.motion.circuit();
    let circuit_input = CircuitEnthalpyInput { circuit: original, ports: input.ports, destinations: input.destinations };
    let (thermal, motion, report) = coupled_driver_candidate(thermal, target_s, circuit_input,
        &[], controls, guess_j, |state, _| {
            let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let bulk = coefficients.bulk;
            contacts::admit_sources(original, &[], || (*cancel.borrow_mut())())?;
            let circuit = input.motion.step_with_contacts_candidate(bulk.resistors, Vec::new(),
                FieldMotionInput { electrical_ports: input.ports, external_force_n: input.external_force_n,
                    duration_s: target_s-thermal.time_s() }, controls.circuit, || (*cancel.borrow_mut())())?;
            Ok(DrivenCoefficients { circuit: BoundaryDriver { inner: circuit, boundaries: coefficients.boundaries },
                sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
                source: bulk.source, circuit_trials: 1 })
        }, || (*cancel.borrow_mut())())?;
    Ok(FieldMotionEnthalpyCandidate { thermal, motion, report })
}

/// The same interval with refreshed material electrical contacts and explicitly
/// declared A/B heat destinations. Source identity and selected contact claims
/// use the SAME continuity check as the stationary electrothermal circuit.
/// All contact temperatures, resistances and thermal links come from the actual
/// proposed endpoint. Contact forces, opening/arcs and interface storage remain
/// outside the original zero-storage electrical-interface model.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: FieldMotionEnthalpyInput<'_>, routes: &[ContactHeatCells],
    controls: CircuitEnthalpySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldMotionEnthalpyCandidate<'m, ContactThermalSample<S>>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let check = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
    let routes = contacts::admit_routes(thermal, input.motion.circuit(), routes, check)?
        .into_values().collect::<Vec<_>>();
    coupled_candidate(thermal, target_s, input, &routes, controls, guess_j, |state| {
        let coefficients = resolve(state)?;
        let bulk = coefficients.bulk;
        Ok(ContactThermalCoefficients { contacts: coefficients.contacts,
            bulk: CircuitThermalCoefficients { resistors: bulk.resistors, sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k,
                source: ContactThermalSample { material: bulk.source, destinations: routes.clone() } } })
    }, || (*cancel.borrow_mut())())
}

#[allow(clippy::too_many_arguments)]
fn coupled_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64,
    input: FieldMotionEnthalpyInput<'_>, routes: &[ContactHeatCells],
    controls: CircuitEnthalpySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldMotionEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    let cancel = RefCell::new(cancelled);
    let original = input.motion.circuit();
    let circuit_input = CircuitEnthalpyInput { circuit: original, ports: input.ports, destinations: input.destinations };
    let (thermal, motion, report) = coupled_driver_candidate(thermal, target_s, circuit_input,
        routes, controls, guess_j, |state, _remaining| {
            let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let bulk = coefficients.bulk;
            contacts::admit_sources(original, &coefficients.contacts, || (*cancel.borrow_mut())())?;
            let circuit = input.motion.step_with_contacts_candidate(bulk.resistors, coefficients.contacts,
                FieldMotionInput { electrical_ports: input.ports, external_force_n: input.external_force_n,
                    duration_s: target_s-thermal.time_s() }, controls.circuit, || (*cancel.borrow_mut())())?;
            Ok(DrivenCoefficients { circuit, sources_w: bulk.sources_w,
                conductances_w_k: bulk.conductances_w_k, source: bulk.source, circuit_trials: 1 })
        }, || (*cancel.borrow_mut())())?;
    Ok(FieldMotionEnthalpyCandidate { thermal, motion, report })
}
