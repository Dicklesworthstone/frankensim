//! Actual conductor temperature, electrical operating point and field motion
//! in the original finite-enthalpy residual. All three histories commit together.

use super::*;
use crate::material_circuit::conductivity::{ConductivitySolve, ConductorResponse};
use crate::material_circuit::conductivity::contact::ContactResponse;
use crate::material_circuit::motion::conductivity::{FieldConductivityCandidate,
    FieldConductivityStep};

/// Native field/body work, actual E/J material data and finite-enthalpy balance.
pub type FieldConductivityEnthalpyReport<S> = CircuitEnthalpyReport<S, FieldConductivityStep>;

/// One admitted thermal, electrical and moving-body endpoint.
pub struct FieldConductivityEnthalpyCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    motion: MaterialFieldMotion,
    report: FieldConductivityEnthalpyReport<S>,
}
impl<'m, S> FieldConductivityEnthalpyCandidate<'m, S> {
    /// Original finite-body temperatures, enthalpies and event cursor.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Joint charge/flux, gap and momentum with actual applied material sources.
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    /// Original field-motion work and joint thermal/material admission.
    pub const fn report(&self) -> &FieldConductivityEnthalpyReport<S> { &self.report }
    /// Publish all histories together; native Joule heat is already deposited.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialFieldMotion,
        FieldConductivityEnthalpyReport<S>) {
        (self.thermal, self.motion, self.report)
    }
}

impl ThermalDriver for FieldConductivityCandidate {
    type Model = MaterialFieldMotion;
    type Step = FieldConductivityStep;
    fn account(&self) -> DriverAccount<'_> { field_account(&self.step().field) }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        FieldConductivityCandidate::into_parts(self)
    }
}

/// Resolve nonlinear conduction and reciprocal actuator motion at the SAME
/// finite thermal endpoint. Original resistor bundles supplied by `resolve`
/// retain their complete temperature/material context; `responses` binds every
/// electrical operating coordinate to the actual moving descriptor's E/J.
/// Only that original descriptor's Joule loss enters the receiving cells.
/// All native descriptor trials, including discarded thermal stencils, count
/// against `maximum_circuit_trials`; no iteration debits accepted history.
///
/// Conductor geometry, body inertia and field coefficients retain the original
/// motion restrictions. This heats stationary electrode/winding resistors, not
/// moving conductors or the lossless field medium. Series contacts refuse this
/// bulk-conduction path. Optional support retains its original fixed thermal
/// state. Endpoint thermal/midpoint electrical coupling is generally first order;
/// source support, numerical convergence and work gates do not certify accuracy.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64, input: FieldMotionEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>], controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    step_candidate_with_transport(thermal, target_s, input, responses, controls,
        conductivity, guess_j, resolve, |_, _, transport| Ok(transport), cancelled)
}

/// Refresh any E/J-dependent thermal links and source observations at the same
/// admitted native field trajectory. `transport` receives the original thermal
/// proposal, actual electrical material sample, and seed transport coefficients;
/// it supplies only independent external heat and original link replacements.
/// It may not add Joule heat or advance physical history. The ordinary entry
/// point keeps thermal transport independent of electrical operating coordinates.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_transport<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64, input: FieldMotionEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>], controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    transport: impl Fn(&DepositionSimulation<'m>, &FieldConductivityStep,
        StepCoefficients<S>) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    coupled_candidate(thermal, target_s, input, responses, None, &[], controls,
        conductivity, guess_j, |state| resolve(state).map(|bulk| ContactThermalCoefficients {
            bulk, contacts: Vec::new(),
        }), transport, cancelled)
}

/// Close source-dependent bulk/contact conduction, reciprocal field motion and
/// finite heating in the SAME thermal residual. The original motion descriptor
/// supplies actual bulk E/J and each contact's own V/I/J, together with gap work,
/// charge/flux and momentum. Complete source bundles refresh at that electrical
/// operating point and the proposed thermal endpoint; no separability is assumed.
///
/// `resolve` supplies every original bulk and contact bundle at the thermal trial.
/// `contact_responses` declares the original interface cards, ordered-surface
/// orientation and electrical coordinates; `None` accepts only electrically
/// independent contacts. Every contact has an explicit A/B heat partition in
/// `routes`. Those native losses reach the finite cells once, separately from
/// bulk losses, reciprocal field work and physical storage.
///
/// `transport` has [`step_candidate_with_transport`]'s semantics and can consume
/// the same admitted contact operating point. Every inner descriptor solve counts
/// against the existing outer work budget, including stencils and final admission.
/// Source, domain, route, budget or cancellation refusal publishes none of the
/// electrical, mechanical or thermal histories.
///
/// Conductors and contact geometry remain stationary; the field gap and original
/// guided body evolve. Constant field coefficients, fixed inertia, zero interface
/// storage, frozen support state and the original endpoint-thermal timing remain
/// the model limits. This does not add motional conductor heat, opening/arcing,
/// field entropy, a live thermoelastic support or a convergence guarantee.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contact_response_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64, input: FieldMotionEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>], contact_responses: Option<&[ContactResponse<'_>]>,
    routes: &[ContactHeatCells], controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(&DepositionSimulation<'m>, &FieldConductivityStep,
        StepCoefficients<S>) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldConductivityEnthalpyCandidate<'m, ContactThermalSample<S>>, CircuitEnthalpyError<E>> {
    coupled_candidate(thermal, target_s, input, responses, contact_responses, routes,
        controls, conductivity, guess_j, |state| resolve(state).map(|resolved| {
            let bulk = resolved.bulk;
            ContactThermalCoefficients { contacts: resolved.contacts,
                bulk: CircuitThermalCoefficients { resistors: bulk.resistors,
                    sources_w: bulk.sources_w, conductances_w_k: bulk.conductances_w_k,
                    source: ContactThermalSample { material: bulk.source, destinations: routes.to_vec() } } }
        }), |state, electrical, coefficients| {
            let destinations = coefficients.payload.destinations;
            transport(state, electrical, StepCoefficients { sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                payload: coefficients.payload.material }).map(|transport| StepCoefficients {
                    sources_w: transport.sources_w, conductances_w_k: transport.conductances_w_k,
                    payload: ContactThermalSample { material: transport.payload, destinations } })
        }, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn coupled_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>, target_s: f64, input: FieldMotionEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>], contact_responses: Option<&[ContactResponse<'_>]>,
    routes: &[ContactHeatCells], controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(&DepositionSimulation<'m>, &FieldConductivityStep,
        StepCoefficients<S>) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<FieldConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    conductivity.validate()?;
    let cancel = RefCell::new(cancelled);
    let circuit_input = CircuitEnthalpyInput { circuit: input.motion.circuit(),
        ports: input.ports, destinations: input.destinations };
    let (thermal, motion, report) = coupled_driver_candidate(thermal, target_s,
        circuit_input, routes, controls, guess_j, |state, remaining| {
            let resolved = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let coefficients = resolved.bulk;
            let bounded = ConductivitySolve {
                maximum_iterations: conductivity.maximum_iterations.min(remaining), ..conductivity };
            let circuit = input.motion.conductivity_with_contact_response_candidate(
                coefficients.resistors, resolved.contacts, responses, contact_responses,
                FieldMotionInput { electrical_ports: input.ports,
                    external_force_n: input.external_force_n, duration_s: target_s-thermal.time_s() },
                controls.circuit, bounded, || (*cancel.borrow_mut())()).map_err(|error| match error {
                    MaterialCircuitError::ConductivityConvergence { .. }
                        if bounded.maximum_iterations < conductivity.maximum_iterations =>
                            CircuitEnthalpyError::CircuitBudget,
                    other => CircuitEnthalpyError::from(other),
                })?;
            let circuit_trials = circuit.step().iterations;
            let transport = transport(state, circuit.step(), StepCoefficients {
                sources_w: coefficients.sources_w, conductances_w_k: coefficients.conductances_w_k,
                payload: coefficients.source }).map_err(CircuitEnthalpyError::Source)?;
            Ok(DrivenCoefficients { circuit, sources_w: transport.sources_w,
                conductances_w_k: transport.conductances_w_k, source: transport.payload, circuit_trials })
        }, || (*cancel.borrow_mut())())?;
    Ok(FieldConductivityEnthalpyCandidate { thermal, motion, report })
}
