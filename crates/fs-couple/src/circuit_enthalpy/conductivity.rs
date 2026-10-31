//! Close electrical field/current response and finite-body heating together.
//!
//! The existing conductivity iterator resolves original material cards at the
//! circuit's actual midpoint fields. The existing enthalpy residual supplies
//! each trial's temperature and deposits the resulting native Joule heat once.
//! Neither trial publishes electrical or thermal history before joint admission.

use super::*;
use crate::material_circuit::conductivity::{
    ConductivitySolve, ConductorResponse, MaterialConductivityCandidate, MaterialConductivityStep,
};
use crate::material_circuit::conductivity::contact::ContactResponse;
use super::contacts::{ContactHeatCells, ContactThermalCoefficients, ContactThermalSample};

/// Original nonlinear electrical response and shared local/global heat account.
pub type ConductivityEnthalpyReport<S> = CircuitEnthalpyReport<S, MaterialConductivityStep>;

/// Electrical and thermal endpoints belonging to the same physical interval.
pub struct ConductivityEnthalpyCandidate<'m, S> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    report: ConductivityEnthalpyReport<S>,
}

impl<'m, S> ConductivityEnthalpyCandidate<'m, S> {
    /// Complete proposed thermal history, preserving the original event cursor.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }

    /// Proposed charge and flux with the original circuit topology and sources.
    pub const fn circuit(&self) -> &MaterialCircuit {
        &self.circuit
    }

    /// Native work, actual-field material queries and finite-body energy balance.
    pub const fn report(&self) -> &ConductivityEnthalpyReport<S> {
        &self.report
    }

    /// Publish both histories together; branch losses are already deposited.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'m>,
        MaterialCircuit,
        ConductivityEnthalpyReport<S>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl ThermalDriver for MaterialConductivityCandidate {
    type Model = MaterialCircuit;
    type Step = MaterialConductivityStep;

    fn account(&self) -> DriverAccount<'_> {
        let electrical = &self.step().electrical;
        DriverAccount {
            resistors: &electrical.resistors,
            contacts: &electrical.contacts,
            heat_j: electrical.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: electrical.stored_after_j - electrical.stored_before_j,
            supplied_j: electrical.native.supplied,
        }
    }

    fn into_parts(self) -> (Self::Model, Self::Step) {
        MaterialConductivityCandidate::into_parts(self)
    }
}

/// Solve heat transport, material response and nonlinear conduction together.
///
/// `resolve` evaluates the original complete resistor bundles at the proposed
/// thermal endpoint through the shared material/thermal state owners.
/// `responses` then binds the actual solved electric field/current density to
/// explicitly named source coordinates. All other coordinates come from those
/// thermal seeds. No separable temperature/field law, metal identification or
/// extrapolation is inferred when a source does not supply such a response.
///
/// Every inner trial starts from the same accepted charge/flux. All inner native
/// circuit solves count against the original outer circuit-work budget, including
/// finite-difference stencils and fresh final admission. Native branch loss is
/// the heat transfer; constitutive iteration error is never deposited as heat.
/// The returned report retains both the applied and actual-field material data.
///
/// This stationary bulk-conductor path holds topology, geometry, C/L, caloric
/// laws and source ports fixed during the interval. Series electrical contacts
/// require a separately coupled contact refresh and are refused by the electrical
/// owner. Thermal timing remains backward Euler with midpoint circuit work:
/// generally first-order coupled accuracy, without a uniqueness or unconditional
/// stability claim. Source limits, budgets and cancellation preserve both inputs.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    step_candidate_with_transport(
        thermal,
        target_s,
        input,
        responses,
        controls,
        conductivity,
        initial_energy_increments_j,
        resolve,
        |_, _, transport| Ok(transport),
        cancelled,
    )
}

/// Resolve thermal transport at the SAME solved electrical operating point.
///
/// The seed callback has the same contract as [`step_candidate`]. After its
/// complete bundles have been refreshed by the original electrical owner,
/// `transport` receives that admitted candidate and the seed thermal data.
/// It must replace any field/current-dependent heat conductances and source
/// observations before the shared backward-Euler residual is evaluated.
/// `StepCoefficients::payload` retains the updated source observations.
/// Its powers are independent external sources, never additional Joule heat.
///
/// This callback advances no state and introduces no additional electrical or
/// thermal solve. Every nonlinear stencil and final admission performs it again;
/// refusal leaves both accepted histories unchanged. The ordinary entry point
/// uses the identity transport update for coefficients independent of E/J.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_transport<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
    transport: impl Fn(
        &DepositionSimulation<'m>,
        &MaterialConductivityStep,
        StepCoefficients<S>,
    ) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    coupled_with_contacts_candidate(
        thermal, target_s, input, responses, None, &[], controls, conductivity,
        initial_energy_increments_j,
        |state| resolve(state).map(|bulk| ContactThermalCoefficients {
            bulk, contacts: Vec::new(),
        }),
        transport, cancelled,
    )
}

/// Close nonlinear bulk conduction, source-refreshed contact drops and finite
/// heat together. Contact Joule heat enters both original bodies through the
/// explicitly authored partitions; bulk losses keep their original destinations.
///
/// `resolve` supplies complete bulk and contact bundles at the SAME proposed
/// thermal state. The inner native circuit/source closure uses the actual bulk
/// voltage, excluding series contact drops. Contact sources may respond to the
/// supplied temperature/phase, but their electrical operating coordinates require
/// a separate interface closure and are refused. Every residual evaluation and
/// final admission restarts from the same accepted physical histories.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    contact_cells: &[ContactHeatCells],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, ContactThermalSample<S>>, CircuitEnthalpyError<E>> {
    step_with_contacts_candidate_with_transport(
        thermal, target_s, input, responses, contact_cells, controls, conductivity,
        initial_energy_increments_j, resolve, |_, _, transport| Ok(transport), cancelled,
    )
}

/// The same contact/bulk closure with thermal transport refreshed at the solved
/// electrical operating point. `transport` has [`step_candidate_with_transport`]'s
/// semantics: it replaces original thermal coefficients and independent powers,
/// never deposits Joule heat or advances state. Explicit contact partitions are
/// retained with the accepted material observations.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts_candidate_with_transport<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    contact_cells: &[ContactHeatCells],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(
        &DepositionSimulation<'m>,
        &MaterialConductivityStep,
        StepCoefficients<S>,
    ) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, ContactThermalSample<S>>, CircuitEnthalpyError<E>> {
    step_with_contact_response_candidate(
        thermal, target_s, input, responses, None, contact_cells, controls,
        conductivity, initial_energy_increments_j, resolve, transport, cancelled,
    )
}

/// Close contact voltage/current response, nonlinear bulk conduction and finite
/// thermal storage in one transaction. `contact_responses` binds every original
/// interface card and its explicit ordered-surface electrical axes. The same
/// native midpoint supplies bulk E/J and individual contact V/I/J; thermal seed
/// temperatures stay fixed throughout that inner constitutive iteration.
///
/// Source temperature and operating-point effects need not be separable: each
/// complete interface bundle is re-evaluated at both sets of coordinates.
/// `transport` may refresh thermal coefficients from the admitted operating
/// point, with the same semantics as [`step_candidate_with_transport`]. Every
/// native bulk/contact loss reaches its authored heat destination exactly once.
/// The original outer circuit budget counts every inner native solve, including
/// failed stencils and final admission. No candidate advances accepted history.
///
/// This retains fixed geometry, instantaneous passive scalar interface response,
/// zero contact storage and the endpoint-thermal/midpoint-electrical time model.
/// It adds no contact thickness/field, switching, breakdown or convergence claim.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contact_response_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    contact_responses: Option<&[ContactResponse<'_>]>,
    contact_cells: &[ContactHeatCells],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(
        &DepositionSimulation<'m>,
        &MaterialConductivityStep,
        StepCoefficients<S>,
    ) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, ContactThermalSample<S>>, CircuitEnthalpyError<E>> {
    coupled_with_contacts_candidate(
        thermal, target_s, input, responses, contact_responses, contact_cells, controls, conductivity,
        initial_energy_increments_j,
        |state| resolve(state).map(|resolved| {
            let bulk = resolved.bulk;
            ContactThermalCoefficients {
                contacts: resolved.contacts,
                bulk: CircuitThermalCoefficients {
                    resistors: bulk.resistors,
                    sources_w: bulk.sources_w,
                    conductances_w_k: bulk.conductances_w_k,
                    source: ContactThermalSample {
                        material: bulk.source,
                        destinations: contact_cells.to_vec(),
                    },
                },
            }
        }),
        |state, electrical, coefficients| {
            let destinations = coefficients.payload.destinations;
            transport(state, electrical, StepCoefficients {
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                payload: coefficients.payload.material,
            }).map(|transport| StepCoefficients {
                sources_w: transport.sources_w,
                conductances_w_k: transport.conductances_w_k,
                payload: ContactThermalSample { material: transport.payload, destinations },
            })
        },
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn coupled_with_contacts_candidate<'m, S, E>(
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    input: CircuitEnthalpyInput<'_>,
    responses: &[ConductorResponse<'_>],
    contact_responses: Option<&[ContactResponse<'_>]>,
    contact_cells: &[ContactHeatCells],
    controls: CircuitEnthalpySolve,
    conductivity: ConductivitySolve,
    initial_energy_increments_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'m>) -> Result<ContactThermalCoefficients<S>, E>,
    transport: impl Fn(
        &DepositionSimulation<'m>,
        &MaterialConductivityStep,
        StepCoefficients<S>,
    ) -> Result<StepCoefficients<S>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<ConductivityEnthalpyCandidate<'m, S>, CircuitEnthalpyError<E>> {
    // Validate the authored policy before the outer work budget narrows it.
    conductivity.validate()?;
    let cancel = RefCell::new(cancelled);
    let (thermal, circuit, report) = coupled_driver_candidate(
        thermal,
        target_s,
        input,
        contact_cells,
        controls,
        initial_energy_increments_j,
        |state, remaining| {
            let resolved = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let coefficients = resolved.bulk;
            let bounded = ConductivitySolve {
                maximum_iterations: conductivity.maximum_iterations.min(remaining),
                ..conductivity
            };
            let circuit = input
                .circuit
                .conductivity_with_contact_response_candidate(
                    coefficients.resistors,
                    resolved.contacts,
                    responses,
                    contact_responses,
                    input.ports,
                    target_s - thermal.time_s(),
                    controls.circuit,
                    bounded,
                    || (*cancel.borrow_mut())(),
                )
                .map_err(|error| match error {
                    MaterialCircuitError::ConductivityConvergence { .. }
                        if bounded.maximum_iterations < conductivity.maximum_iterations =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => CircuitEnthalpyError::from(other),
                })?;
            let circuit_trials = circuit.step().iterations;
            let transport = transport(
                state,
                circuit.step(),
                StepCoefficients {
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    payload: coefficients.source,
                },
            )
            .map_err(CircuitEnthalpyError::Source)?;
            Ok(DrivenCoefficients {
                circuit,
                sources_w: transport.sources_w,
                conductances_w_k: transport.conductances_w_k,
                source: transport.payload,
                circuit_trials,
            })
        },
        || (*cancel.borrow_mut())(),
    )?;
    Ok(ConductivityEnthalpyCandidate {
        thermal,
        circuit,
        report,
    })
}
