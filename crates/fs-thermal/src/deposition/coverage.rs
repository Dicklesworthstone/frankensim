//! Separate material-surface measures with one conservative heat exchange.
use super::*;

/// A thermal exchange link and the area it covers on each incident cell.
/// Areas use the SAME measure as each cell's `surface_area_m2` and boundary
/// coefficient. They can be different material-reference areas. The link's
/// nominal area can, for example, be the common projected footprint; it does
/// not replace either covered area and never implicitly scales conductance.
///
/// One supplied total conductance still exchanges exactly opposite powers.
/// This does not infer conductance, pressure, adhesion or an area Jacobian.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactPatch {
    pub link: ThermalLink,
    /// Covered areas [on link.a, on link.b], m2; positive even if insulating.
    pub covered_area_m2: [f64; 2],
}
impl ContactPatch {
    /// The original symmetric link interpretation, without extra state.
    pub const fn symmetric(link: ThermalLink) -> Self {
        Self { covered_area_m2: [link.area_m2; 2], link }
    }
}

/// Complete expected and replacement interface, including both surface measures.
/// None denotes separation. Repartition all affected patches in one batch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactPatchUpdate {
    pub id: u64,
    pub expected: Option<ContactPatch>,
    pub replacement: Option<ContactPatch>,
}

impl DepositionSimulation<'_> {
    /// Authoritative current link AND per-side coverage; None means absent.
    pub fn contact_patch(&self, id: u64) -> Option<ContactPatch> {
        self.links.get(&id).map(|link| ContactPatch {
            link: *link, covered_area_m2: self.contact_areas(link),
        })
    }

    pub(super) fn contact_areas(&self, link: &ThermalLink) -> [f64; 2] {
        self.covered_areas.get(&link.id).copied().unwrap_or([link.area_m2; 2])
    }

    /// Current free areas for a bounded, unique selection of active cells.
    ///
    /// One ordered contact sweep replaces a full link scan per selected cell.
    /// Each cell retains the exact scalar `exposed_area_m2` subtraction order:
    /// reservoir first, then incident links in canonical link-ID order. Distinct
    /// per-side coverage is authoritative even on zero-conductance contacts.
    /// Input order cannot change output bits; returned keys are canonical.
    ///
    /// Uses O(S) additional storage for S selected cells and one contact sweep
    /// with bounded tree-map lookups, rather than S full scans. Empty selection
    /// skips the link sweep. Counts use the original cell budget. Cancellation
    /// is checked every 256 selected cells/links and before publication.
    pub fn exposed_areas_m2(
        &self,
        selected: &[CellId],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        checkpoint()?;
        if selected.len() > self.budget.max_cells {
            return Err(DepositionError::Budget("exposed surface selection"));
        }
        let mut areas = BTreeMap::new();
        for (index, &id) in selected.iter().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            let cell = self.cells.get(&id)
                .ok_or(DepositionError::InvalidInput("inactive cell"))?;
            let area = cell.geometry.surface_area_m2
                - cell.boundary.reservoir.map_or(0.0, |r| r.area_m2);
            if areas.insert(id, area).is_some() {
                return Err(DepositionError::InvalidInput("duplicate exposed surface cell"));
            }
        }
        if !areas.is_empty() {
            for (index, link) in self.links.values().enumerate() {
                if index % 256 == 0 { checkpoint()?; }
                let covered = self.contact_areas(link);
                for (id, cover) in [(link.a, covered[0]), (link.b, covered[1])] {
                    if let Some(area) = areas.get_mut(&id) { *area -= cover; }
                }
            }
        }
        for (index, &area) in areas.values().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if !nonnegative(area) {
                return Err(DepositionError::InvalidInput("contacts exceed surface area"));
            }
        }
        checkpoint()?;
        Ok(areas)
    }

    // Shared admission for scheduled symmetric links and explicit pair patches.
    // The caller owns a scratch candidate; a refusal never escapes as a state.
    pub(super) fn insert_contact_patch(&mut self, patch: ContactPatch) -> Result<(), DepositionError> {
        self.insert_contact_patch_structure(patch)?;
        self.exposed_area_m2(patch.link.a)?;
        self.exposed_area_m2(patch.link.b)?;
        Ok(())
    }

    // Batch repartition uses the SAME structural checks, then admits final
    // coverage in one sweep before its private candidate can be published.
    fn insert_contact_patch_structure(&mut self, patch: ContactPatch) -> Result<(), DepositionError> {
        let link = patch.link;
        if self.links.len() >= self.budget.max_links {
            return Err(DepositionError::Budget("active links"));
        }
        if self.links.contains_key(&link.id) || link.a == link.b
            || !positive(link.area_m2) || !nonnegative(link.conductance_w_k)
            || !self.cells.contains_key(&link.a) || !self.cells.contains_key(&link.b) {
            return Err(DepositionError::InvalidInput(
                "contact requires distinct active cells and a unique ID"));
        }
        if patch.covered_area_m2.iter().any(|&area| !positive(area)) {
            return Err(DepositionError::InvalidInput("positive finite covered area required on each side"));
        }
        self.links.insert(link.id, link);
        if patch.covered_area_m2 != [link.area_m2; 2] {
            self.covered_areas.insert(link.id, patch.covered_area_m2);
        }
        Ok(())
    }

    /// Replace a full batch at the current time; process due schedule events first.
    /// Expected values include coverage, so changed surface measures cannot pass
    /// a stale proposal simply because the exchange conductance is unchanged.
    /// Removing every changed patch before installing any replacement permits
    /// nonmatching-mesh repartition without transient double coverage.
    ///
    /// State, masses, enthalpies, boundary coefficients and time are unchanged on
    /// refusal. Successful changes affect the NEXT thermal interval. Both exposed
    /// area and its CFL contribution use the stored side-specific measure. A
    /// scheduled separation also clears coverage; a later ordinary Contact event
    /// with the same ID creates a new symmetric patch, not stale retained area.
    ///
    /// Coverage is checked once on the complete candidate over all touched cells,
    /// not once per inserted patch. O(updates + active links + active cells)
    /// storage including the original candidate clone; area admission uses one
    /// ordered contact sweep. Loops poll cancellation and obey existing budgets.
    /// This is the explicit deposition lane, not the separate implicit parcel lane.
    pub fn replace_contact_patches(&self, updates: &[ContactPatchUpdate],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>) -> Result<Self, DepositionError> {
        checkpoint()?;
        if updates.len() > self.budget.max_events {
            return Err(DepositionError::Budget("contact updates"));
        }
        if self.schedule.get(self.cursor).is_some_and(|e| e.time_s <= self.time_s) {
            return Err(DepositionError::InvalidInput("process due scheduled events before contact updates"));
        }
        let mut ids = std::collections::BTreeSet::new();
        let mut touched = std::collections::BTreeSet::new();
        for update in updates {
            checkpoint()?;
            if !ids.insert(update.id)
                || update.expected.is_some_and(|patch| patch.link.id != update.id)
                || update.replacement.is_some_and(|patch| patch.link.id != update.id) {
                return Err(DepositionError::InvalidInput("duplicate or inconsistent contact update identity"));
            }
            if self.contact_patch(update.id) != update.expected {
                return Err(DepositionError::InvalidInput("stale thermal contact update"));
            }
            if update.expected != update.replacement {
                for patch in [update.expected, update.replacement].into_iter().flatten() {
                    touched.insert(patch.link.a);
                    touched.insert(patch.link.b);
                }
            }
        }
        let mut next = self.clone();
        let mut audit = EnergyAudit { incoming_mass_kg: 0.0, incoming_energy_j: 0.0,
            external_heat_j: 0.0, stored_energy_change_j: 0.0, residual_j: 0.0,
            substeps: 0, events: 0 };
        for update in updates {
            checkpoint()?;
            if update.expected != update.replacement && update.expected.is_some() {
                next.apply(Event::Separation(update.id), &mut audit)?;
            }
        }
        for update in updates {
            checkpoint()?;
            if update.expected != update.replacement {
                if let Some(patch) = update.replacement { next.insert_contact_patch_structure(patch)?; }
            }
        }
        let selected: Vec<_> = touched.into_iter().collect();
        next.exposed_areas_m2(&selected, &mut checkpoint)?;
        checkpoint()?;
        Ok(next)
    }
}

#[cfg(test)]
#[path = "coverage/tests.rs"]
mod tests;
