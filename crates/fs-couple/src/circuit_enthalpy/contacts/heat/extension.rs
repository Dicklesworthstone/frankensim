//! Continue complete material heat transfer after an explicitly authored topology change.
use super::*;
use std::collections::{BTreeMap, BTreeSet};
use fs_conduction::transient::deposition::contact::MAX_MATERIAL_NETWORK_CONTACTS;
use fs_matdb::QueryPoint;

impl PreparedHeatNetwork<'_, '_> {
    // This consumer owns both original electrical contact work and two-sided
    // heat. Unlike field_material_network(), it must not discard those ports.
    pub(crate) fn field_storage_material_network(
        &self,
    ) -> &crate::circuit_enthalpy::material::MaterialNetwork<'_> {
        self.network.electrical.bulk
    }

    /// Retain bulk, electrical-interface, thermal-interface and emissivity sources
    /// while binding a complete successor heat network around the SAME circuit.
    ///
    /// `surface_successors[i]` identifies old enclosure surface `i` in the new
    /// enclosure. `None` explicitly retires that surface; every old surface needs
    /// an entry, and surviving surfaces map one-to-one. Several faces may belong
    /// to the same cell, so matching by cell alone would be ambiguous. Unmapped
    /// new surfaces use their authored policies. Retained surfaces keep their
    /// exact source and material conditions even when their index/area changes.
    ///
    /// Active old thermal interfaces cannot silently become bulk-only links.
    /// Removed thermal links retire their interface; new links use explicit new
    /// policies. Retained interfaces keep ordered bodies and fixed conditions.
    /// The electrical owner preserves its stricter fixed branch/contact geometry.
    ///
    /// The caller supplies new trace geometry, exposed areas and reciprocal view
    /// factors; this does not infer occlusion, contact or a moving enclosure.
    /// No time, enthalpy, charge, flux or constitutive history advances here.
    /// All original material/geometry evaluators admit the result before return;
    /// expired retained claims refuse instead of selecting a fallback.
    pub fn extend_fixed_circuit<'n, 'a>(
        &self, cx: &Cx<'_>, network: &'n MaterialHeatNetwork<'a>,
        thermal: &DepositionSimulation<'_>, surface_successors: &[Option<usize>],
    ) -> Result<PreparedHeatNetwork<'n, 'a>, HeatNetworkError> {
        check(cx)?;
        let invalid = HeatNetworkError::Invalid;
        if self.integrated_conduction.is_some() {
            return Err(invalid("integrated conduction geometry must be explicitly rebound after topology changes"));
        }
        let old_count = self.network.radiation.map_or(0, |r| r.enclosure.surfaces().len());
        let new_count = network.radiation.map_or(0, |r| r.enclosure.surfaces().len());
        if surface_successors.len() != old_count
            || network.radiation.is_some_and(|r| r.cells.len() != new_count)
            || network.thermal_contacts.len() > MAX_MATERIAL_NETWORK_CONTACTS {
            return Err(invalid("complete bounded surface transition and thermal contacts required"));
        }
        // Original bulk/contact extension owns circuit and material continuity.
        let electrical = self.electrical.extend_fixed_circuit(
            network.electrical, thermal, || cx.checkpoint().is_err())?;
        let mut contacts = network.thermal_contacts.to_vec();
        let mut contact_indices = BTreeMap::new();
        for (i, contact) in contacts.iter().enumerate() {
            check(cx)?;
            if contact_indices.insert(contact.link, i).is_some() {
                return Err(invalid("duplicate successor thermal interface"));
            }
        }
        for old in &self.thermal_contacts {
            check(cx)?;
            if !thermal.links().contains_key(&old.link) { continue; }
            let i = *contact_indices.get(&old.link)
                .ok_or(invalid("an active thermal interface cannot lose its source binding"))?;
            let next = &mut contacts[i];
            let a = old.point;
            let b = next.point;
            if old.sides.iter().zip(next.sides).any(|(x, y)|
                    x.cell != y.cell || x.material_id != y.material_id || !std::ptr::eq(x.curve, y.curve))
                || a.interface.content_hash() != b.interface.content_hash()
                || a.materials.iter().zip(b.materials).any(|(x, y)| x.content_hash() != y.content_hash())
                || a.temperature_axes != b.temperature_axes
                || a.liquid_fraction_axes != b.liquid_fraction_axes
                || a.mean_temperature_axis != b.mean_temperature_axis
                || !same_context(a.context, b.context, |axis| a.temperature_axes.contains(&axis)
                    || a.liquid_fraction_axes.contains(&Some(axis)) || a.mean_temperature_axis == Some(axis)) {
                return Err(invalid("retained thermal interfaces require the same ordered materials and conditions"));
            }
            next.selection = old.selection;
        }
        let mut surfaces = Vec::with_capacity(new_count);
        let mut surface_indices = BTreeMap::new();
        if let Some(radiation) = network.radiation {
            for binding in radiation.cells {
                check(cx)?;
                if binding.surface >= new_count
                    || surface_indices.insert(binding.surface, surfaces.len()).is_some() {
                    return Err(invalid("duplicate or foreign successor radiative surface"));
                }
                surfaces.push(RadiatingCell { surface: binding.surface, cell: binding.cell,
                    material_id: binding.material_id, curve: binding.curve, material: binding.material,
                    selection: binding.selection.clone() });
            }
        }
        let mut targets = BTreeSet::new();
        for old in &self.radiating_cells {
            check(cx)?;
            let Some(target) = surface_successors[old.surface] else { continue; };
            if !targets.insert(target) {
                return Err(invalid("radiative surface continuation must be one-to-one"));
            }
            let i = *surface_indices.get(&target)
                .ok_or(invalid("continued radiative surface is absent from successor enclosure"))?;
            let next = &mut surfaces[i];
            let a = old.material;
            let b = next.material;
            if old.cell != next.cell || old.material_id != next.material_id
                || !std::ptr::eq(old.curve, next.curve)
                || a.material.content_hash() != b.material.content_hash()
                || a.temperature_axis != b.temperature_axis
                || a.liquid_fraction_axis != b.liquid_fraction_axis
                || !same_context(a.context, b.context, |axis|
                    axis == a.temperature_axis || Some(axis) == a.liquid_fraction_axis) {
                return Err(invalid("continued radiative surface must retain its material cell and fixed conditions"));
            }
            next.selection = old.selection.clone();
        }
        // Query original laws with retained pins BEFORE considering new policies.
        // No temporary stage borrow escapes; the returned object borrows network.
        network.prepare_bound(cx, thermal, electrical, &contacts, &surfaces)
    }
}

// Compare only fixed coordinates; actual thermal coordinates are rebuilt by the
// original state owners. Missing axes and changed quantity conventions still fail.
fn same_context(a: &QueryPoint, b: &QueryPoint, evolving: impl Fn(&str) -> bool) -> bool {
    a.axis_quantities() == b.axis_quantities() && a.axes().len() == b.axes().len()
        && a.axes().iter().all(|(key, value)| b.axes().get(key)
            .is_some_and(|other| evolving(key.as_str()) || value.to_bits() == other.to_bits()))
}
