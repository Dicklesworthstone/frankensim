//! Source consistency between the native structural and thermal owners.
use super::*;

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
        for bound in &self.cells {
            cx.checkpoint().map_err(|_| SolidGasError::Cancelled)?;
            let sample = resolved
                .source
                .cells
                .binary_search_by_key(&bound.source.cell, |sample| sample.cell)
                .ok()
                .map(|index| &resolved.source.cells[index])
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
        Ok((resolved.conductances_w_k, Some(resolved.source)))
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
