//! Actual material inventory and enthalpy leaving the original thermal cells.
//!
//! The physical owner supplies removed volume, original material density and
//! complete replacement geometry. This owner checks their consistency, removes
//! the represented mass and volume, and exports that material's CURRENT specific
//! enthalpy. The remaining specific enthalpy is unchanged. There is no latent
//! heat, additional thermal sink, inferred density, or debris heating law.
//!
//! Every old link incident on a changed cell needs an explicit patch update,
//! including an explicitly retained patch when its conductance remains valid.
//! The geometry owner supplies poses, surface measures and new conductances;
//! volume alone determines none of them. Existing coverage and caloric admission
//! apply to the complete endpoint before publication. Phase, strain, momentum,
//! absolute entropy and subsequent debris transport remain caller-owned.

use std::collections::{BTreeMap, BTreeSet};

use super::{
    CellGeometry, CellId, ContactPatchUpdate, DepositedCell, DepositionError, DepositionSimulation,
    MaterialId, ThermalStorageBasis, finite, nonnegative, positive,
};

/// One original material inventory and its physically supplied removal endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialRemoval {
    /// Complete accepted cell; any stale identity, mass, h, geometry or boundary refuses.
    pub expected: DepositedCell,
    /// Requested volume leaving this material [m³], not a heat or damage proxy.
    pub removed_volume_m3: f64,
    /// Positive density [kg/m³] supplied by the original physical material owner.
    /// Both stored inventories must agree with it within the numerical tolerance.
    pub density_kg_m3: f64,
    /// Caller-owned pose, volume and surface measure after this removal.
    pub replacement_geometry: CellGeometry,
    /// Original registered enthalpy law, including its reference-energy datum.
    pub caloric_identity: [u8; 32],
}

/// Absolute arithmetic discrepancy limits; these do not cover physical model error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialRemovalTolerances {
    /// Density/inventory consistency and represented mass decrement [kg].
    pub mass_kg: f64,
    /// Requested versus represented volume decrement [m³].
    pub volume_m3: f64,
}

/// Actual material transport from one cell, in its original caloric convention.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialRemovalCellReceipt {
    pub cell_id: CellId,
    pub material_id: MaterialId,
    pub caloric_identity: [u8; 32],
    pub density_kg_m3: f64,
    pub previous_mass_kg: f64,
    pub mass_kg: f64,
    pub previous_volume_m3: f64,
    pub volume_m3: f64,
    pub requested_volume_m3: f64,
    /// Difference of the actual stored volumes, not an unrepresented request.
    pub removed_volume_m3: f64,
    /// Difference of the actual stored masses; this is the exported inventory.
    pub removed_mass_kg: f64,
    /// Unchanged specific enthalpy of both remaining and removed material [J/kg].
    pub carried_h_j_kg: f64,
    /// Signed outward enthalpy [J]; negative values are valid for negative datums.
    pub outgoing_enthalpy_j: f64,
    /// Stable change in retained m*h [J], evaluated without subtracting gross energies.
    pub retained_energy_change_j: f64,
    /// Actual mass decrement minus density times requested volume [kg].
    pub mass_roundoff_kg: f64,
    /// Actual volume decrement minus requested volume [m³].
    pub volume_roundoff_m3: f64,
    /// Largest |m-rho*V| over original, remaining and removed inventory [kg].
    pub maximum_density_residual_kg: f64,
}

/// Canonically ordered outflow and stable retained-plus-exported balances.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialRemovalAudit {
    pub cells: Vec<MaterialRemovalCellReceipt>,
    pub removed_mass_kg: f64,
    pub removed_volume_m3: f64,
    pub outgoing_enthalpy_j: f64,
    pub retained_energy_change_j: f64,
    pub mass_balance_residual_kg: f64,
    pub energy_residual_j: f64,
}

impl DepositionSimulation<'_> {
    /// Remove well-mixed material and replace its geometry at the current time.
    ///
    /// Process all due schedule events first. Supplied material phase, density,
    /// geometry and conductances remain the external physical owner's choices;
    /// only enthalpy storage can supply this open-material transport. Retained
    /// temperature is unchanged because the SAME caloric law and specific h
    /// survive. This operation advances neither heat nor time.
    ///
    /// Every incident old contact on a changed cell must appear in `patches`.
    /// New coverage is validated after all cells and patches change together.
    /// Cell exhaustion and unrepresentable positive removal refuse. One graph
    /// clone and O(selected cells + patches) staging retain transactional failure.
    pub fn remove_material_candidate(
        &self,
        removals: &[MaterialRemoval],
        patches: &[ContactPatchUpdate],
        tolerances: MaterialRemovalTolerances,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, MaterialRemovalAudit), DepositionError> {
        self.remove_material_with_event_order(removals, patches, tolerances, false, checkpoint)
    }

    /// Remove material at an explicitly chosen LEFT limit of pending endpoint events.
    ///
    /// A coupled heat owner may finish an interval with events exactly at the
    /// current time still pending. This ingress permits those events, rejects
    /// past-due events, and leaves the original clock, schedule and cursor intact.
    /// It does not process or skip a boundary/contact event. The caller must
    /// process the same pending events before starting the next heat interval.
    pub fn remove_material_before_pending_events(
        &self,
        removals: &[MaterialRemoval],
        patches: &[ContactPatchUpdate],
        tolerances: MaterialRemovalTolerances,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, MaterialRemovalAudit), DepositionError> {
        self.remove_material_with_event_order(removals, patches, tolerances, true, checkpoint)
    }

    fn remove_material_with_event_order(
        &self,
        removals: &[MaterialRemoval],
        patches: &[ContactPatchUpdate],
        tolerances: MaterialRemovalTolerances,
        before_pending_events: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, MaterialRemovalAudit), DepositionError> {
        checkpoint()?;
        if !nonnegative(tolerances.mass_kg) || !nonnegative(tolerances.volume_m3) {
            return Err(invalid(
                "finite nonnegative material-removal tolerances required",
            ));
        }
        if removals.len() > self.budget.max_cells || patches.len() > self.budget.max_events {
            return Err(DepositionError::Budget(
                "material-removal cells or contact updates",
            ));
        }
        if self.schedule.get(self.cursor).is_some_and(|event| {
            event.time_s < self.time_s || (!before_pending_events && event.time_s == self.time_s)
        }) {
            return Err(invalid(
                "process due schedule events before material removal",
            ));
        }
        let mut ordered = BTreeMap::new();
        for (index, removal) in removals.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if ordered.insert(removal.expected.id, removal).is_some() {
                return Err(invalid("duplicate material-removal cell"));
            }
        }
        let mut next = self.clone();
        let mut changed = BTreeSet::new();
        let mut receipts = Vec::with_capacity(ordered.len());
        let mut sums = [[0.0; 2]; 6];
        for removal in ordered.values() {
            checkpoint()?;
            let id = removal.expected.id;
            let cell = self
                .cells
                .get(&id)
                .ok_or(invalid("inactive material-removal cell"))?;
            if cell != &removal.expected {
                return Err(invalid("stale material-removal cell"));
            }
            let law = self.materials[&cell.material_id];
            if law.storage_basis() != ThermalStorageBasis::Enthalpy
                || law.caloric_identity() != Some(removal.caloric_identity)
            {
                return Err(invalid(
                    "material removal requires the exact original enthalpy law and datum",
                ));
            }
            let volume = removal.removed_volume_m3;
            let density = removal.density_kg_m3;
            let geometry = removal.replacement_geometry;
            if !nonnegative(volume)
                || volume >= cell.geometry.volume_m3
                || !positive(density)
                || !positive(geometry.volume_m3)
                || !positive(geometry.surface_area_m2)
                || !geometry
                    .start_m
                    .iter()
                    .chain(&geometry.end_m)
                    .all(|x| x.is_finite())
                || geometry.start_m == geometry.end_m
            {
                return Err(invalid(
                    "positive remaining material and finite replacement geometry required",
                ));
            }
            let density_mass = finite(density * cell.geometry.volume_m3)?;
            let requested_mass = finite(density * volume)?;
            let mass = finite(cell.mass_kg - requested_mass)?;
            if !positive(density_mass) || !positive(mass) || requested_mass >= cell.mass_kg {
                return Err(invalid(
                    "material removal leaves no positive original inventory",
                ));
            }
            let removed_mass = finite(cell.mass_kg - mass)?;
            let removed_volume = finite(cell.geometry.volume_m3 - geometry.volume_m3)?;
            if volume == 0.0 {
                if removed_mass != 0.0 || removed_volume != 0.0 {
                    return Err(invalid(
                        "zero requested removal must retain mass and occupied volume",
                    ));
                }
            } else if requested_mass == 0.0 || removed_mass <= 0.0 || removed_volume <= 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            let mass_roundoff = finite(removed_mass - requested_mass)?;
            let volume_roundoff = finite(removed_volume - volume)?;
            let remaining_density_mass = finite(density * geometry.volume_m3)?;
            let removed_density_mass = finite(density * removed_volume)?;
            if !positive(remaining_density_mass)
                || (removed_volume > 0.0 && removed_density_mass == 0.0)
            {
                return Err(DepositionError::Arithmetic);
            }
            let density_residual = finite(cell.mass_kg - density_mass)?
                .abs()
                .max(finite(mass - remaining_density_mass)?.abs())
                .max(finite(removed_mass - removed_density_mass)?.abs());
            if mass_roundoff.abs() > tolerances.mass_kg
                || volume_roundoff.abs() > tolerances.volume_m3
                || density_residual > tolerances.mass_kg
            {
                return Err(invalid(
                    "material removal exceeds its density, mass or volume discrepancy limit",
                ));
            }
            if !positive(law.temperature_k(cell.specific_h_j_kg)?) {
                return Err(DepositionError::MaterialDomain);
            }
            let outgoing = finite(removed_mass * cell.specific_h_j_kg)?;
            if removed_mass > 0.0 && cell.specific_h_j_kg != 0.0 && outgoing == 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            let retained = -outgoing;
            let mut candidate = cell.clone();
            candidate.mass_kg = mass;
            candidate.geometry = geometry;
            law.validate_cell(&candidate)?;
            if &candidate != cell {
                changed.insert(id);
            }
            next.cells.insert(id, candidate);
            let values = [
                removed_mass,
                removed_volume,
                outgoing,
                retained,
                finite((mass - cell.mass_kg) + removed_mass)?,
                finite(retained + outgoing)?,
            ];
            for (sum, value) in sums.iter_mut().zip(values) {
                add(sum, value)?;
            }
            receipts.push(MaterialRemovalCellReceipt {
                cell_id: id,
                material_id: cell.material_id,
                caloric_identity: removal.caloric_identity,
                density_kg_m3: density,
                previous_mass_kg: cell.mass_kg,
                mass_kg: mass,
                previous_volume_m3: cell.geometry.volume_m3,
                volume_m3: geometry.volume_m3,
                requested_volume_m3: volume,
                removed_volume_m3: removed_volume,
                removed_mass_kg: removed_mass,
                carried_h_j_kg: cell.specific_h_j_kg,
                outgoing_enthalpy_j: outgoing,
                retained_energy_change_j: retained,
                mass_roundoff_kg: mass_roundoff,
                volume_roundoff_m3: volume_roundoff,
                maximum_density_residual_kg: density_residual,
            });
        }
        let mut updated = BTreeSet::new();
        for (index, patch) in patches.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if ![patch.expected, patch.replacement]
                .into_iter()
                .flatten()
                .any(|p| ordered.contains_key(&p.link.a) || ordered.contains_key(&p.link.b))
            {
                return Err(invalid(
                    "material removal cannot replace unrelated contact patches",
                ));
            }
            updated.insert(patch.id);
        }
        for (index, link) in self.links.values().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            if (changed.contains(&link.a) || changed.contains(&link.b))
                && !updated.contains(&link.id)
            {
                return Err(invalid(
                    "changed material geometry requires every incident contact to be explicitly updated",
                ));
            }
        }
        // This shared in-place operation belongs only to the private candidate.
        // It owns stale-patch, topology, conductance, budget and coverage checks.
        next.apply_contact_patch_updates(patches, &mut checkpoint)?;
        // Explicitly retained patches may not be in that operation's changed
        // set; their coverage must also fit every newly reduced cell surface.
        let selected: Vec<_> = changed.into_iter().collect();
        next.exposed_areas_m2(&selected, &mut checkpoint)?;
        let [
            mass,
            volume,
            outgoing,
            retained,
            mass_residual,
            energy_residual,
        ] = sums.map(|s| s[0]);
        checkpoint()?;
        Ok((
            next,
            MaterialRemovalAudit {
                cells: receipts,
                removed_mass_kg: mass,
                removed_volume_m3: volume,
                outgoing_enthalpy_j: outgoing,
                retained_energy_change_j: retained,
                mass_balance_residual_kg: mass_residual,
                energy_residual_j: energy_residual,
            },
        ))
    }
}

fn invalid(message: &'static str) -> DepositionError {
    DepositionError::InvalidInput(message)
}

fn add(sum: &mut [f64; 2], value: f64) -> Result<(), DepositionError> {
    let corrected = finite(value - sum[1])?;
    let next = finite(sum[0] + corrected)?;
    sum[1] = finite((next - sum[0]) - corrected)?;
    sum[0] = next;
    Ok(())
}

#[cfg(test)]
#[path = "removal/tests.rs"]
mod tests;
