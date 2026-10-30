//! Kelvin conversion, sourced interfaces and radiation on one accepted endpoint.

use super::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve};
use crate::circuit_enthalpy::contacts::{admit_sources, advance_circuit};
use crate::circuit_enthalpy::thermoelectric::{
    Error, ThermoelectricCandidate, ThermoelectricDriver, ThermoelectricLeg,
    ThermoelectricNetworkError, ThermoelectricSample,
    coupled_resolved_candidate, finish_candidate_with_internal_heat,
};
use crate::circuit_enthalpy::thermoelectric::motion::{self, ThermoelectricMotionCandidate};
use crate::material_circuit::{MaterialCircuit, MaterialResistor};
use crate::material_circuit::contacts::SeriesElectricalContact;
use crate::material_circuit::motion::{FieldMotionInput, MaterialFieldMotion};
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
        let (next, electrical, report) = self.thermoelectric_resolved_candidate(
            cx, thermal, target_s, circuit, ports, legs, controls, initial_energy_increments_j,
            |resistors, contacts, inputs, dt, cancelled| {
                advance_circuit(circuit, resistors, contacts, inputs, dt,
                    controls.circuit, cancelled)
            },
        )?;
        let internal_heat_j = (target_s - thermal.time_s())
            * report.source.network.radiation_net_power_w();
        finish_candidate_with_internal_heat(next, electrical, report,
            internal_heat_j, controls, || cx.checkpoint().is_err())
    }

    /// Convert a finite thermal gradient into reciprocal field-actuator motion
    /// while contacts, conductive layers and radiation load the SAME junctions.
    ///
    /// The original Kelvin residual computes Seebeck voltage and deposits signed
    /// Peltier/Thomson heat from the actual native circuit current. The original
    /// field descriptor advances charge/flux, gap and momentum together; its
    /// electrical back-reaction is part of that same current. Bulk and electrical
    /// contact Joule heat reaches the original caloric cells exactly once.
    /// Material thermal interfaces, integrated k(T) layers and enclosure radiation
    /// are refreshed inside every proposed enthalpy solve, so changing their data
    /// changes temperatures, generated current and physical actuator motion.
    ///
    /// All restrictions of `step_thermoelectric_candidate` apply, including
    /// refusing unresolved contact temperature jumps on Kelvin legs. The moving
    /// field medium, guided body and optional axial support retain their original
    /// admitted material snapshots and inertia. This does not infer field entropy,
    /// live support expansion, hysteresis, Lorentz transport or deforming thermal
    /// geometry. The original first-order thermal split and motion limits apply.
    ///
    /// Closed radiation is internal, not external heat. Joint publication requires
    /// field, kinetic, support and finite thermal storage to close against actual
    /// external electrical work, force-displacement work and authored heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_thermoelectric_motion_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        model: &MaterialFieldMotion, ports: &[f64], external_force_n: f64,
        legs: &[ThermoelectricLeg], controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ThermoelectricMotionCandidate<'m, HeatNetworkSample>, Error> {
        let circuit = model.circuit();
        let (next, model, report) = self.thermoelectric_resolved_candidate(
            cx, thermal, target_s, circuit, ports, legs, controls, initial_energy_increments_j,
            |resistors, contacts, inputs, dt, cancelled| {
                admit_sources(circuit, &contacts, &mut *cancelled)?;
                Ok(model.step_with_contacts_candidate(resistors, contacts,
                    FieldMotionInput { electrical_ports: inputs, external_force_n,
                        duration_s: dt }, controls.circuit, cancelled)?)
            },
        )?;
        let internal_heat_j = (target_s - thermal.time_s())
            * report.source.network.radiation_net_power_w();
        motion::finish_candidate(next, model, report, internal_heat_j,
            controls, || cx.checkpoint().is_err())
    }

    // One admission/resolution path for both physical drivers. This prevents
    // future thermal mechanisms or Kelvin-interface restrictions from silently
    // disappearing in one composition while remaining present in the other.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn thermoelectric_resolved_candidate<'m, D: ThermoelectricDriver>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], legs: &[ThermoelectricLeg],
        controls: CircuitEnthalpySolve, initial_energy_increments_j: &[f64],
        advance: impl Fn(Vec<MaterialResistor>, Vec<SeriesElectricalContact>, &[f64], f64,
            &mut dyn FnMut() -> bool) -> Result<D, Error>,
    ) -> Result<(DepositionSimulation<'m>, D::Model,
        CircuitEnthalpyReport<ThermoelectricSample<HeatNetworkSample>, D::Step>), Error> {
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
            }, |source| &source.material, advance, false, || cx.checkpoint().is_err(),
        )?;
        Ok((next, electrical, report))
    }
}
