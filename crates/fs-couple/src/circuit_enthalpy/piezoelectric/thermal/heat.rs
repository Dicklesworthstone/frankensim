//! Complete sourced heat transport on the original temperature-dependent ceramic solve.
use super::{
    Error, PiezoelectricCircuit, PiezoelectricThermalBinding, PiezoelectricThermalCandidate,
    PiezoelectricThermalError, resolved_candidate,
};
use crate::circuit_enthalpy::contacts::heat::{
    HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork,
};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, admit_driver_unknowns,
    field_storage::admit_circuit_sources,
};
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(PiezoelectricThermalError::Heat(other.into())),
    }
}

/// Close coolant, ceramic response, rod motion and circuit history on the same
/// original enthalpy residual using any complete material heat assembly.
///
/// The assembly resolves material transport at each proposed thermal endpoint.
/// Its finite coolant links remain internal to the original thermal inventory;
/// prescribed reservoirs remain external heat. Compliance, signed coupling and
/// permittivity are re-evaluated from the ceramic's original thermal source at
/// that same endpoint. The existing driver alone computes field entropy,
/// reversible ceramic heat, native Joule heat and charge/strain evolution.
///
/// The exact prepared owner supplies material bindings and heat destinations.
/// Accepted conductor identities are checked before any trial; electrical
/// contacts refuse because this rod/circuit owner has no contact-work law.
/// The complete accepted assembly sample remains in
/// `report().coupled.source.network`, with typed coolant/source failures
/// preserved under `PiezoelectricThermalError::Heat`. Closed-radiation net heat
/// is excluded from external supply by the same physical internal-energy gate.
///
/// The existing uniform, fully solid, small-signal ceramic, source domains,
/// whole-volume binding and positive combined heat-capacity envelope apply.
/// Neither thermal strain, depoling nor dielectric loss is inferred. Heat
/// geometry and contact/view-factor topology stay fixed; midpoint motion with
/// endpoint thermal coefficients is first order in general. Every trial starts
/// from the same accepted history, and refusal publishes no partial endpoint.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t, A: HeatAssembly + ?Sized>(
    assembly: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    model: &PiezoelectricCircuit,
    electrical_ports: &[f64],
    force_n: f64,
    ceramic: &PiezoelectricThermalBinding,
    controls: CircuitEnthalpySolve,
    guess_j: &[f64],
) -> Result<PiezoelectricThermalCandidate<'t, A::Sample>, Error<A::Error>> {
    cx.checkpoint()
        .map_err(|_| CircuitEnthalpyError::Cancelled)?;
    admit_driver_unknowns(thermal, controls, guess_j)?;
    let prepared = assembly.prepared_heat_network();
    let material = prepared
        .field_material_network()
        .map_err(heat_error::<A::Error>)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error::<A::Error>)?;
    let mut initial = assembly
        .resolve_heat(cx, thermal)
        .map_err(|error| error.map_source(PiezoelectricThermalError::Heat))?;
    if !initial.contacts.is_empty() {
        return Err(CircuitEnthalpyError::Invalid(
            "thermal piezoelectric circuit cannot discard electrical contact work or heat",
        ));
    }
    admit_circuit_sources(
        model.circuit(),
        &mut initial.bulk.resistors,
        &mut initial.contacts,
        || cx.checkpoint().is_err(),
    )?;
    drop(initial);
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
            let resolved = assembly
                .resolve_heat(cx, trial)
                .map_err(|error| error.map_source(PiezoelectricThermalError::Heat))?;
            if !resolved.contacts.is_empty() {
                return Err(CircuitEnthalpyError::Invalid(
                    "thermal piezoelectric circuit cannot discard electrical contact work or heat",
                ));
            }
            Ok(resolved.bulk)
        },
        |sample: &A::Sample| sample.as_ref().radiation_net_power_w(),
        || cx.checkpoint().is_err(),
    )
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
        step_candidate(
            self,
            cx,
            thermal,
            target_s,
            model,
            electrical_ports,
            force_n,
            ceramic,
            controls,
            guess_j,
        )
    }
}
