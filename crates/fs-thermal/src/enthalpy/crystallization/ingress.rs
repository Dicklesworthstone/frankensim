//! Explicit original phase history at the thermal network's mass boundary.

use std::collections::BTreeSet;

use super::*;

impl<'a> CrystallizationSimulation<'a> {
    /// Adopt active cells and scheduled feed with their original phase history.
    ///
    /// Each supplied state belongs to exactly one active cell or one unconsumed
    /// birth in the original timeline. The material's complete kinetic/caloric
    /// identity and reconstructed state must validate, and stored specific
    /// enthalpy must match the cell exactly. This transfers an accepted history;
    /// it neither infers age from temperature nor replaces incoming energy with
    /// a caloric round trip. Omitted kinetic cells retain the explicitly
    /// all-amorphous convention of [`Self::new`]. Ordinary materials cannot carry
    /// a crystal state.
    ///
    /// Future histories remain immutable until their own birth, including births
    /// exactly at the requested endpoint. Imported crystals enter the original
    /// mass/enthalpy account once and do not count as newly crystallized mass.
    /// A clone retains the complete pending feed; implicit candidates and event
    /// advances use the same ingress. No residence time before the declared
    /// birth is simulated or inferred.
    ///
    /// Supplied keys must be unique and all must match. Ambiguous repeated
    /// births, stale enthalpy, another law, missing cells and exhausted original
    /// cell budgets refuse the whole attachment. Birth geometry, original
    /// material validation, boundary coverage and all remaining event rules
    /// still belong to the original thermal owner. The caloric controls only
    /// validate supplied state reconstruction; they do not bound later time
    /// integration or grant new source-domain authority.
    pub fn new_with_incoming_states(
        thermal: DepositionSimulation<'a>,
        materials: &[CrystallizationMaterial<'a>],
        incoming_states: &[(CellId, CrystallizationState)],
        caloric: CrystallizationSolve,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, CrystallizationTimelineError> {
        checkpoint()?;
        caloric.validate()?;
        if materials.len() > thermal.budget.max_cells {
            return Err(DepositionError::Budget("crystallization materials").into());
        }
        if incoming_states.len() > thermal.budget.max_cells {
            return Err(DepositionError::Budget("incoming crystallization states").into());
        }
        let mut map = BTreeMap::new();
        for &material in materials {
            checkpoint()?;
            if material.amorphous.storage_basis() != ThermalStorageBasis::Enthalpy
                || !positive(material.maximum_temperature_slope)
                || material.maximum_temperature_slope
                    < material.amorphous.maximum_temperature_slope()
                || !thermal.uses_material_law(material.material_id, material.amorphous)
                || map.insert(material.material_id, material).is_some()
            {
                return Err(DepositionError::InvalidInput(
                    "crystallization needs unique original enthalpy laws and a valid mixture slope bound",
                ).into());
            }
            let original_identity = material.amorphous.caloric_identity();
            if original_identity.is_none()
                || original_identity != material.law.parent_caloric_identity()
            {
                return Err(DepositionError::InvalidInput(
                    "crystallization requires the same identified parent calorics as the original thermal law",
                ).into());
            }
        }

        let mut incoming = BTreeMap::new();
        let mut unmatched = BTreeSet::new();
        for &(id, state) in incoming_states {
            checkpoint()?;
            if incoming.insert(id, state).is_some() {
                return Err(DepositionError::InvalidInput(
                    "incoming crystallization cell IDs must be unique",
                )
                .into());
            }
            unmatched.insert(id);
        }
        let mut states = BTreeMap::new();
        for (&id, cell) in &thermal.cells {
            checkpoint()?;
            if let Some(&state) = incoming.get(&id) {
                admit_incoming(cell, state, &map, caloric, &mut unmatched)?;
                states.insert(id, state);
            } else if let Some(model) = map.get(&cell.material_id) {
                states.insert(id, amorphous_state(&thermal, id, model, &mut checkpoint)?);
            }
        }
        // Scan the immutable future schedule once, rather than searching it
        // independently for every feed cell. The default empty ingress adds no
        // schedule walk or per-cell history allocation.
        if !incoming.is_empty() {
            for timed in thermal.pending_events() {
                checkpoint()?;
                if let Event::Birth(cell) = &timed.event {
                    if let Some(&state) = incoming.get(&cell.id) {
                        admit_incoming(cell, state, &map, caloric, &mut unmatched)?;
                    }
                }
            }
        }
        if !unmatched.is_empty() {
            return Err(DepositionError::InvalidInput(
                "incoming crystallization history has no active cell or pending birth",
            )
            .into());
        }
        checkpoint()?;
        Ok(Self {
            thermal,
            materials: map,
            states,
            incoming_states: Arc::new(incoming),
            remelting: BTreeMap::new(),
            adaptivity: BTreeMap::new(),
        })
    }
}

fn admit_incoming(
    cell: &DepositedCell,
    state: CrystallizationState,
    materials: &BTreeMap<MaterialId, CrystallizationMaterial<'_>>,
    caloric: CrystallizationSolve,
    unmatched: &mut BTreeSet<CellId>,
) -> Result<(), CrystallizationTimelineError> {
    if !unmatched.remove(&cell.id) {
        return Err(DepositionError::InvalidInput(
            "incoming crystallization history has ambiguous repeated births",
        )
        .into());
    }
    let material = materials
        .get(&cell.material_id)
        .ok_or(DepositionError::InvalidInput(
            "incoming crystal history requires a bound crystallization material",
        ))?;
    material.law.validate_state(state, caloric)?;
    if state.specific_enthalpy_j_kg().to_bits() != cell.specific_h_j_kg.to_bits() {
        return Err(CrystallizationError::StateMismatch.into());
    }
    Ok(())
}
