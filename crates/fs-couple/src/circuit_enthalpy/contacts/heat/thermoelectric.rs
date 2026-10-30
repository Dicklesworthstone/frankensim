//! Kelvin conversion, sourced interfaces and radiation on one accepted endpoint.

use super::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use crate::circuit_enthalpy::contacts::advance_circuit;
use crate::circuit_enthalpy::thermoelectric::{
    Error, ThermoelectricCandidate, ThermoelectricLeg, ThermoelectricNetworkError,
    coupled_resolved_candidate, finish_candidate_with_internal_heat,
};
use crate::material_circuit::MaterialCircuit;
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ThermoelectricNetworkError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Solve Seebeck/Peltier/Thomson conversion with the complete heat assembly.
    ///
    /// Every nonlinear trial uses the SAME material temperatures for the original
    /// Kelvin legs, bulk/contact resistances, thermal interfaces, integrated layer
    /// fluxes and enclosure radiosity. Native circuit currents supply bulk and
    /// contact Joule heat exactly once; reversible leg heat retains its original
    /// junction/bulk destinations. The accepted source payload retains every heat
    /// mechanism. Closed radiation is internal exchange and is excluded from the
    /// reported external heat under the original total-energy tolerance.
    ///
    /// Kelvin legs retain their original ideal-junction temperature trace.
    /// Thermal contacts and integrated layer replacements therefore belong on
    /// SEPARATE heat-sink/junction-loading links, and electrical contacts on
    /// ordinary load branches. A contact directly on a leg requires independently
    /// resolved interface temperatures and thermoelectric jump physics; neither
    /// is inferred from an Ohmic or Fourier resistance. Represent distinct face
    /// temperatures with explicit native thermal cells. Surface radiation may
    /// exchange heat directly with the original leg or junction cells.
    ///
    /// Source pins, geometry, phase domains, finite thermal storage, native
    /// circuit dynamics and first-order endpoint-coupling limits are unchanged.
    /// No coefficient is selected by component name. Failure or cancellation
    /// returns no partially accepted electrical or thermal history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_thermoelectric_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
        controls: CircuitEnthalpySolve, initial_energy_increments_j: &[f64],
    ) -> Result<ThermoelectricCandidate<'m, HeatNetworkSample>, Error> {
        super::check(cx).map_err(heat_error)?;
        for leg in legs {
            super::check(cx).map_err(heat_error)?;
            if self.network.electrical.contacts.iter()
                .any(|contact| contact.branch == leg.resistor_branch) {
                return Err(CircuitEnthalpyError::Invalid(
                    "electrical contact on a Kelvin leg requires a resolved thermoelectric interface",
                ));
            }
        }
        let routes = self.network.electrical.heat_destinations(thermal)
            .map_err(HeatNetworkError::from).map_err(heat_error)?;
        let (next, electrical, report, _) = coupled_resolved_candidate(
            self.network.electrical.bulk, thermal, target_s, circuit, ports, legs,
            &routes, controls, initial_energy_increments_j,
            |state, _| {
                let coefficients = self.resolve(cx, state).map_err(heat_error)?;
                let sample = &coefficients.bulk.source;
                for leg in legs {
                    super::check(cx).map_err(heat_error)?;
                    if sample.thermal_contacts.iter().any(|x| leg.thermal_links.contains(&x.link))
                        || sample.integrated_conduction.iter().any(|x| leg.thermal_links.contains(&x.link)) {
                        return Err(CircuitEnthalpyError::Invalid(
                            "a Kelvin edge requires its original junction-to-bulk temperature trace",
                        ));
                    }
                }
                Ok(coefficients)
            }, |source| &source.material,
            |resistors, contacts, inputs, dt, cancelled| {
                advance_circuit(circuit, resistors, contacts, inputs, dt,
                    controls.circuit, cancelled)
            }, false, || cx.checkpoint().is_err(),
        )?;
        let internal_heat_j = (target_s - thermal.time_s())
            * report.source.network.radiation_net_power_w();
        finish_candidate_with_internal_heat(next, electrical, report,
            internal_heat_j, controls, || cx.checkpoint().is_err())
    }
}
