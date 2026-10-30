//! Complete material heat transport with the original field/circuit energy owner.
use super::{FieldThermalBinding, FieldThermalCandidate, FieldThermalError, driven_candidate};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpySolve};
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork};
use crate::material_circuit::{MaterialCircuit, MaterialCircuitStep};
use fs_exec::Cx;
use fs_thermal::enthalpy::DepositionSimulation;

/// The original field candidate carrying the COMPLETE accepted heat-network
/// sample. No thermal state, circuit history or material data is copied to a
/// parallel owner; consume `into_parts()` and publish both original histories.
pub type HeatFieldCandidate<'m> =
    FieldThermalCandidate<'m, MaterialCircuit, MaterialCircuitStep, HeatNetworkSample>;

/// Retain the original typed heat-transfer and field-source failures.
#[derive(Debug)]
pub enum HeatFieldError {
    Heat(HeatNetworkError),
    Field(FieldThermalError),
}
impl From<FieldThermalError> for HeatFieldError {
    fn from(error: FieldThermalError) -> Self { Self::Field(error) }
}
impl core::fmt::Display for HeatFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material field and heat transport: {self:?}")
    }
}
impl std::error::Error for HeatFieldError {}
fn heat_error(error: HeatNetworkError) -> CircuitEnthalpyError<HeatFieldError> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(HeatFieldError::Heat(other)),
    }
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
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let destinations = self.resistor_destinations();
        let routes = self.heat_destinations(thermal).map_err(heat_error)?;
        driven_candidate(
            self.field_storage_material_network(), thermal, target_s, circuit,
            ports, fields, &destinations, &routes, controls, guess_j,
            |trial| self.resolve(cx, trial).map_err(heat_error),
            HeatNetworkSample::radiation_net_power_w,
            || cx.checkpoint().is_err(),
        )
    }
}
