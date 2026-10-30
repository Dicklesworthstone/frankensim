//! Complete material heat transport with the original field/circuit energy owner.
use super::{FieldThermalBinding, FieldThermalCandidate, FieldThermalError, driven_candidate};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use crate::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitStep};
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

/// The original field candidate carrying the COMPLETE accepted heat-network
/// sample. No thermal state, circuit history or material data is copied to a
/// parallel owner; consume `into_parts()` and publish both original histories.
pub type HeatFieldCandidate<'m, S = HeatNetworkSample> =
    FieldThermalCandidate<'m, MaterialCircuit, MaterialCircuitStep, S>;

/// Retain the original typed heat-transfer and field-source failures.
#[derive(Debug)]
pub enum HeatFieldError<E = HeatNetworkError> {
    Heat(E),
    Field(FieldThermalError),
}
impl<E> From<FieldThermalError> for HeatFieldError<E> {
    fn from(error: FieldThermalError) -> Self { Self::Field(error) }
}
impl<E: core::fmt::Debug> core::fmt::Display for HeatFieldError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material field and heat transport: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for HeatFieldError<E> {}
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> CircuitEnthalpyError<HeatFieldError<E>> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HeatFieldError::Heat(other.into())),
    }
}

/// Close the original field/circuit and caloric histories through any complete
/// material heat assembly, including source-resolved finite coolant.
///
/// The assembly's original prepared owner supplies physical source identity and
/// heat routes. Its one complete evaluation at each trial supplies material
/// transport, extra reservoirs and finite links; the accepted payload remains in
/// `report.coupled.source.network`. The original field driver alone owns C(T),
/// L(T), conserved-coordinate rebinding, reversible heat and entropy.
///
/// Finite-link exchange remains internal native conduction. Prescribed sources
/// remain external heat. The retained enclosure's signed numerical net is
/// excluded by the SAME physical field-plus-caloric energy gate. No heat or
/// field state is advanced outside that existing coupled residual.
///
/// Original solid field-source domains, fixed geometry, positive combined
/// capacities, event cuts and temporal-refinement requirements still apply.
#[allow(clippy::too_many_arguments)]
pub fn field_storage_candidate<'m, A: HeatAssembly + ?Sized>(
    assembly: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    circuit: &MaterialCircuit,
    ports: &[f64],
    fields: &[FieldThermalBinding],
    controls: CircuitEnthalpySolve,
    guess_j: &[f64],
) -> Result<HeatFieldCandidate<'m, A::Sample>, CircuitEnthalpyError<HeatFieldError<A::Error>>> {
    cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
    let prepared = assembly.prepared_heat_network();
    let destinations = prepared.resistor_destinations();
    let routes = prepared.heat_destinations(thermal).map_err(heat_error::<A::Error>)?;
    driven_candidate(
        prepared.field_storage_material_network(), thermal, target_s, circuit,
        ports, fields, &destinations, &routes, controls, guess_j,
        |trial| assembly.resolve_heat(cx, trial)
            .map_err(|error| error.map_source(HeatFieldError::Heat)),
        |sample: &A::Sample| sample.as_ref().radiation_net_power_w(),
        || cx.checkpoint().is_err(),
    )
}

impl PreparedHeatNetwork<'_, '_> {
    /// Close material C(T)/L(T), electrical interfaces, thermal contacts,
    /// source-integrated layer conduction and gray-diffuse radiation at ONE
    /// finite-thermal endpoint. The same accepted charge/flux starts every
    /// trial; field parameter work and entropy supply reversible caloric heat,
    /// separately from native bulk/interface Joule loss.
    ///
    /// This uses the same field driver as bulk-only and contact-only circuits.
    /// The existing heat network resolves transport once per trial, preserving
    /// its source pins, area/contact identities and original enclosure. Report
    /// `coupled.source.network` retains every accepted transport calculation.
    ///
    /// Closed radiation is internal exchange, including its numerical net:
    /// `external_heat_j` excludes that net and reversible field heat. The SAME
    /// total-energy tolerance gates physical field U plus caloric storage. No
    /// source correction or second thermal solve hides enclosure roundoff.
    ///
    /// Geometry, topology and view factors are fixed over the step. Existing
    /// whole-volume field bindings, solid affine field laws, positive capacity
    /// envelopes and coordinate limits remain mandatory. Endpoint thermal and
    /// midpoint circuit coupling is first order and requires time refinement;
    /// no contact mechanics, arc, magnetic force or dispersive loss is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn field_storage_candidate<'m>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        circuit: &MaterialCircuit,
        ports: &[f64],
        fields: &[FieldThermalBinding],
        controls: CircuitEnthalpySolve,
        guess_j: &[f64],
    ) -> Result<HeatFieldCandidate<'m>, CircuitEnthalpyError<HeatFieldError>> {
        field_storage_candidate(
            self, cx, thermal, target_s, circuit, ports, fields, controls, guess_j,
        )
    }
}
