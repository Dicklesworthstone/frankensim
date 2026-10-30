//! Re-evaluate a source-bearing gas at a different property temperature.
use super::*;

impl ThermochemicalGasModel {
    /// Resolve a property temperature at the snapshot's ACTUAL pressure,
    /// composition, fixed source conditions and selected transport claims.
    ///
    /// `self` is a source bank: its initial fractions need not match the
    /// snapshot, but it must contain every active species with the exact
    /// caloric/transport sources. Inactive bank species remain inactive.
    /// This supports films around transported or reacting gas without replacing
    /// the finite gas's energy, density, acoustic state or composition.
    ///
    /// The original bulk evaluation is checked against the retained laws before
    /// changing temperature. No claim is reselected and source-domain failures
    /// are not extrapolated. The result is a property evaluation, not a physical
    /// update of the original parcel or a finite layer with independent storage.
    pub fn evaluate_temperature_from_snapshot(
        &self,
        bulk: &ResolvedThermochemicalGasState,
        temperature: Temperature,
    ) -> Result<ResolvedThermochemicalGasState, MaterialError> {
        if !positive(temperature.value()) {
            return Err(invalid("property temperature must be finite and positive"));
        }
        let receipt = bulk.thermodynamics().receipt();
        let bank = self.thermodynamics.components();
        if receipt
            .species()
            .iter()
            .any(|id| !bank.iter().any(|x| x.species() == id))
        {
            return Err(invalid("property source bank lacks an active gas species"));
        }
        let fractions = bank
            .iter()
            .map(|component| {
                receipt
                    .species()
                    .iter()
                    .position(|id| id == component.species())
                    .map_or(0.0, |i| f64::from_bits(receipt.declared_fraction_bits()[i]))
            })
            .collect();
        let composition = fs_thermochem::Composition::new(receipt.declared_basis(), fractions)
            .map_err(|e| failure("snapshot composition", e))?;
        let mut model = self.with_composition(composition)?;
        if model.transport.len() != bulk.components.len() || model.source_packs != bulk.source_packs
        {
            return Err(invalid(
                "property source bank differs from the original active gas sources",
            ));
        }
        for (data, resolved) in model.transport.iter_mut().zip(&bulk.components) {
            if data.material.content_hash() != resolved.card_identity() {
                return Err(invalid(
                    "property source bank changes a gas transport material",
                ));
            }
            data.selection = MaterialPropertySelection::PinnedByProperty(
                resolved
                    .properties()
                    .iter()
                    .map(|property| {
                        (
                            property.requirement().name().to_owned(),
                            property.answer().receipt.selected,
                        )
                    })
                    .collect(),
            );
        }
        let context = bulk
            .components
            .first()
            .ok_or_else(|| invalid("gas snapshot lacks its original material context"))?
            .query();
        let admitted = model.evaluate(&context)?;
        // Pinning changes query-policy receipts, not selected physics. Compare
        // original caloric receipts as well as the complete numerical gas view.
        if admitted.thermodynamics != bulk.thermodynamics
            || admitted.state != bulk.state
            || admitted.specific_internal_energy_j_kg != bulk.specific_internal_energy_j_kg
            || admitted.specific_enthalpy_j_kg != bulk.specific_enthalpy_j_kg
        {
            return Err(invalid(
                "property source bank cannot reproduce the original gas snapshot",
            ));
        }
        let point = context
            .with("temperature", temperature.value())
            .map_err(|e| failure("property temperature coordinate", e))?;
        model.evaluate(&point)
    }
}
