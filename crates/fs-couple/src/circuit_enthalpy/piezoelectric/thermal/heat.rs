//! Thermal contacts, layers and closed radiation on the original ceramic solve.
use super::{
    Error, PiezoelectricCircuit, PiezoelectricThermalBinding, PiezoelectricThermalCandidate,
    PiezoelectricThermalError, resolved_candidate,
};
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork,
};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(PiezoelectricThermalError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Advance the original ceramic, circuit and thermal bodies with their
    /// complete sourced heat network at one common nonlinear endpoint.
    ///
    /// Temperature-dependent compliance, signed piezoelectric coupling and
    /// permittivity determine the same rod motion, charge and reversible heat
    /// as the bulk-only thermal entry. Thermal contacts, source-integrated
    /// conduction layers and gray-diffuse radiation resolve once at every trial,
    /// retaining their original source selections and full accepted sample in
    /// `report().coupled.source.network`.
    ///
    /// Closed radiation is internal exchange. Its computed net remains visible
    /// in `closed_radiation_net_j`, is excluded from `external_heat_j`, and
    /// passes the same total-energy gate with physical ceramic internal energy,
    /// kinetic energy, other circuit storage and caloric storage. Neither
    /// reversible ceramic heat nor resistor Joule heat is deposited twice.
    ///
    /// The original small-signal, fully solid, uniform ceramic and positive
    /// capacity envelope remain mandatory. Heat-network geometry, contact area
    /// and view factors stay fixed over this step; this does not infer changing
    /// contact mechanics, thermal expansion, dielectric loss or depoling.
    /// Electrical interfaces require a native work/heat owner and refuse here.
    /// Endpoint thermal and midpoint mechanical/electrical composition is first
    /// order in general. Refusal/cancellation returns no partial history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_piezoelectric_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        model: &PiezoelectricCircuit,
        electrical_ports: &[f64],
        force_n: f64,
        ceramic: &PiezoelectricThermalBinding,
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
    ) -> Result<PiezoelectricThermalCandidate<'t, HeatNetworkSample>, Error> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let material = self.field_material_network().map_err(heat_error)?;
        resolved_candidate(
            material,
            thermal,
            target_s,
            model,
            electrical_ports,
            force_n,
            ceramic,
            controls,
            guess_j,
            |trial| {
                self.resolve(cx, trial)
                    .map(|resolved| resolved.bulk)
                    .map_err(heat_error)
            },
            HeatNetworkSample::radiation_net_power_w,
            || cx.checkpoint().is_err(),
        )
    }
}
