//! Preserve bulk and ordered-interface sources across thermal topology growth.
use super::*;

impl PreparedElectricalContactNetwork<'_, '_> {
    /// Extend the thermal network around the SAME circuit and electrical contacts.
    ///
    /// Existing bulk source selections are transferred by the original bulk
    /// extension. Interface claims, A/B material order, fixed conditions, branch
    /// assignments, contact area and explicit heat fractions remain unchanged.
    /// New thermal cells/links are explicit; no electrical branch or interface is
    /// inferred. Values are evaluated at the supplied state, never cached.
    ///
    /// Use after accepted deposition events, retaining the accepted circuit's
    /// charges/fluxes. A newly eligible competing interface fit cannot replace a
    /// retained claim; an expired claim refuses. Intentional changes to electrical
    /// contact topology, geometry or source conditions require explicit preparation
    /// and circuit reconfiguration, not this source-continuity operation.
    pub fn extend_fixed_circuit<'n, 'a>(
        &self,
        network: &'n ElectricalContactNetwork<'a>,
        thermal: &DepositionSimulation<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedElectricalContactNetwork<'n, 'a>, ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        self.admit_geometry(thermal, &mut cancelled)?;
        if network.contacts.len() != self.network.contacts.len() {
            return Err(Error::Invalid("extension cannot add or remove electrical contacts"));
        }
        let old_bases: BTreeMap<_, _> = self.network.electrical_bases.iter().copied().collect();
        let new_bases: BTreeMap<_, _> = network.electrical_bases.iter().copied().collect();
        if network.electrical_bases.len() != old_bases.len()
            || new_bases.len() != network.electrical_bases.len()
            || old_bases != new_bases
        {
            return Err(Error::Invalid("extension must retain each electrical property basis"));
        }
        let mut remaining: BTreeMap<_, _> = self.network.contacts.iter()
            .map(|contact| (contact.link, contact)).collect();
        for contact in network.contacts {
            if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
            let old = remaining.remove(&contact.link)
                .ok_or(Error::Invalid("duplicate or changed electrical contact ID"))?;
            let a = old.point;
            let b = contact.point;
            let evolving = |axis: &str| a.temperature_axes.contains(&axis)
                || a.liquid_fraction_axes.contains(&Some(axis))
                || a.mean_temperature_axis == Some(axis);
            if contact.branch != old.branch
                || contact.heat_fraction_to_first.to_bits() != old.heat_fraction_to_first.to_bits()
                || a.interface.content_hash() != b.interface.content_hash()
                || a.materials.iter().zip(b.materials).any(|(x, y)| x.content_hash() != y.content_hash())
                || a.temperature_axes != b.temperature_axes
                || a.liquid_fraction_axes != b.liquid_fraction_axes
                || a.mean_temperature_axis != b.mean_temperature_axis
                || a.context.axis_quantities() != b.context.axis_quantities()
                || a.context.axes().len() != b.context.axes().len()
                || a.context.axes().iter().any(|(key, value)| b.context.axes().get(key)
                    .is_none_or(|other| !evolving(key.as_str()) && value.to_bits() != other.to_bits()))
            {
                return Err(Error::Invalid("extension must retain ordered interface sources, conditions and heat routes"));
            }
        }
        if !remaining.is_empty() {
            return Err(Error::Invalid("existing electrical contact is absent from successor"));
        }
        network.routes(thermal)?;
        let bulk = self.bulk.extend_fixed_circuit(network.bulk, thermal, &mut cancelled)
            .map_err(Error::Bulk)?;
        let candidate = PreparedElectricalContactNetwork {
            network, bulk, pins: self.pins.clone(), geometry: self.geometry.clone(),
        };
        // The original evaluator admits the complete successor with retained
        // claims. Do not first reselect its interfaces using successor policies.
        candidate.resolve(thermal, &mut cancelled)?;
        if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
        Ok(candidate)
    }

    // Source identity alone cannot detect a retargeted link between equal-material
    // cells or a changed area. Keep immutable geometry, not physical history.
    pub(super) fn admit_geometry(
        &self, thermal: &DepositionSimulation<'_>, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
        for (&id, &(cells, area)) in &self.geometry {
            if cancelled() { return Err(Error::Thermal(DepositionError::Cancelled)); }
            let link = thermal.links().get(&id)
                .ok_or(Error::Invalid("prepared electrical contact link is no longer active"))?;
            if [link.a, link.b] != cells || link.area_m2.to_bits() != area.to_bits() {
                return Err(Error::Invalid("prepared electrical contact endpoints and area must remain fixed"));
            }
        }
        Ok(())
    }
}
