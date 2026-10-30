//! Source consistency between the native structural and thermal owners.
use super::*;

// The original residual consumes coefficients and retains one source payload.
// A complete heat sample remains intact; its material snapshot is borrowed by
// the gas patches rather than cloned into a competing report field.
pub(super) trait GasTransportSample {
    fn material(&self) -> &MaterialNetworkSample;
    fn radiation_net_power_w(&self) -> f64 {
        0.0
    }
}
impl GasTransportSample for MaterialNetworkSample {
    fn material(&self) -> &MaterialNetworkSample {
        self
    }
}
impl GasTransportSample for HeatNetworkSample {
    fn material(&self) -> &MaterialNetworkSample {
        &self.material
    }
    fn radiation_net_power_w(&self) -> f64 {
        HeatNetworkSample::radiation_net_power_w(self)
    }
}

pub(super) struct GasTransportState<S> {
    pub sources_w: BTreeMap<CellId, f64>,
    pub conductances_w_k: BTreeMap<u64, f64>,
    pub sample: Option<S>,
}
pub(super) trait GasHeatTransport {
    type Sample: GasTransportSample;
    fn resolve(
        &self,
        solid: &SolidEnthalpy<'_>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<GasTransportState<Self::Sample>, SolidGasError>;
}

impl GasHeatTransport for Option<&PreparedMaterialNetwork<'_, '_>> {
    type Sample = MaterialNetworkSample;
    fn resolve(
        &self,
        solid: &SolidEnthalpy<'_>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<GasTransportState<Self::Sample>, SolidGasError> {
        let (conductances_w_k, sample) = solid.resolve_gas_transport(*self, thermal, cx)?;
        Ok(GasTransportState {
            sources_w: BTreeMap::new(),
            conductances_w_k,
            sample,
        })
    }
}

impl GasHeatTransport for &PreparedHeatNetwork<'_, '_> {
    type Sample = HeatNetworkSample;
    fn resolve(
        &self,
        solid: &SolidEnthalpy<'_>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<GasTransportState<Self::Sample>, SolidGasError> {
        // No graph circuit is being advanced. Even a zero-current electrical
        // binding must not silently disappear from this composition.
        let network = self.thermal_material_network()?;
        for bound in &solid.cells {
            cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
            let cell = network
                .cells
                .iter()
                .find(|cell| cell.cell == bound.source.cell)
                .ok_or(SolidGasError::Invalid(
                    "structural cell lacks its original heat-network binding",
                ))?;
            let accepted = thermal
                .cells()
                .get(&cell.cell)
                .ok_or(SolidGasError::Invalid(
                    "heat-network structural cell is inactive",
                ))?;
            if cell.material_id != accepted.material_id
                || !std::ptr::eq(cell.curve, bound.source.curve)
                || cell.point.material.content_hash() != bound.source.point.material.content_hash()
                || cell.point.temperature_axis != bound.source.point.temperature_axis
                || cell.point.liquid_fraction_axis != bound.source.point.liquid_fraction_axis
            {
                return Err(SolidGasError::Invalid(
                    "heat and solid must retain the same cell material and caloric owner",
                ));
            }
            let phase = bound
                .source
                .curve
                .state_at_specific_enthalpy(accepted.specific_h_j_kg)
                .map_err(SolidEnthalpyError::from)?;
            if cell.point.at(phase)? != bound.source.point.at(phase)? {
                return Err(SolidGasError::Invalid(
                    "heat and solid must retain the same complete physical conditions",
                ));
            }
        }
        let resolved = (*self).resolve(cx, thermal)?;
        if !resolved.contacts.is_empty() || !resolved.bulk.resistors.is_empty() {
            return Err(SolidGasError::Invalid(
                "finite gas/solid heat cannot discard electrical work",
            ));
        }
        solid.admit_gas_transport(&resolved.bulk.source.material, cx)?;
        Ok(GasTransportState {
            sources_w: resolved.bulk.sources_w,
            conductances_w_k: resolved.bulk.conductances_w_k,
            sample: Some(resolved.bulk.source),
        })
    }
}

impl SolidEnthalpy<'_> {
    pub(super) fn resolve_gas_transport(
        &self,
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        thermal: &DepositionSimulation<'_>,
        cx: &Cx<'_>,
    ) -> Result<(BTreeMap<u64, f64>, Option<MaterialNetworkSample>), SolidGasError> {
        let Some(transport) = transport else {
            return Ok((BTreeMap::new(), None));
        };
        let resolved = transport.resolve(thermal, || cx.checkpoint().is_err())?;
        if !resolved.resistors.is_empty() || !resolved.sources_w.is_empty() {
            return Err(SolidGasError::Invalid(
                "gas/solid transport requires a resistor-free network with no additional work source",
            ));
        }
        self.admit_gas_transport(&resolved.source, cx)?;
        Ok((resolved.conductances_w_k, Some(resolved.source)))
    }

    fn admit_gas_transport(
        &self,
        source: &MaterialNetworkSample,
        cx: &Cx<'_>,
    ) -> Result<(), SolidGasError> {
        for bound in &self.cells {
            cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
            let sample = source
                .cells
                .binary_search_by_key(&bound.source.cell, |sample| sample.cell)
                .ok()
                .map(|index| &source.cells[index])
                .ok_or(SolidGasError::Invalid(
                    "structural cell lacks its original transport source",
                ))?;
            if let Some(material) = &sample.conductivity {
                let expected = bound.source.point.at(sample.phase)?;
                if material.card_identity() != bound.source.point.material.content_hash()
                    || material.query() != expected
                {
                    return Err(SolidGasError::Invalid(
                        "structural and transport sources must share the original card and physical coordinates",
                    ));
                }
            }
        }
        Ok(())
    }
}

pub(super) fn cell_material(
    transport: Option<&MaterialNetworkSample>,
    cell: CellId,
) -> Result<Option<&ResolvedMaterialStatePoint>, SolidGasError> {
    let Some(transport) = transport else {
        return Ok(None);
    };
    let sample = transport
        .cells
        .binary_search_by_key(&cell, |sample| sample.cell)
        .ok()
        .map(|index| &transport.cells[index])
        .ok_or(SolidGasError::Invalid(
            "gas patch lacks an original material cell",
        ))?;
    // Isolated cells may omit bulk conductivity under the original network's
    // rules. Gas patches then retain their original source-resolving path.
    // Otherwise the layer owner consumes this exact snapshot and admits its
    // complete physical point, property key and explicit reference claim pin.
    Ok(sample.conductivity.as_ref())
}
