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
    // Validate the authored policy before the outer work budget narrows it.
    conductivity.validate()?;
    let cancel = RefCell::new(cancelled);
    let (thermal, circuit, report) = coupled_driver_candidate(
        thermal,
        target_s,
        input,
        &[],
        controls,
        initial_energy_increments_j,
        |state, remaining| {
            let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
            let bounded = ConductivitySolve {
                maximum_iterations: conductivity.maximum_iterations.min(remaining),
                ..conductivity
            };
            let circuit = input
                .circuit
                .conductivity_with_materials_candidate(
                    coefficients.resistors,
                    responses,
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
            Ok(DrivenCoefficients {
                circuit,
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                source: coefficients.source,
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
