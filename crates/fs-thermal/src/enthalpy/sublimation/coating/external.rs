//! Thermal admission for a joint substrate/coating/externally owned gas solve.
//!
//! The original substrate, deposited material and native heat/event history stay
//! together. The gas owner supplies complete supported endpoints and owns their
//! phase-transfer conservation; this boundary never adds a second latent heater.
//! Prescribed reservoir values are pressure/temperature observations only. No
//! finite gas inventory is copied into a coating or advanced by the heat path.

use super::{AdvanceMode, CoatingSimulation, CoatingState, CoatingSurface, initial_state};
use crate::enthalpy::sublimation::ThermalSublimationError;
use crate::enthalpy::{
    CellId, DepositedCell, DepositionBudget, DepositionError, EnergyAudit, StepConfig, nonnegative,
};
use crate::sublimation::coating::SupportedSurfaceState;
use crate::sublimation::{SolverConfig, VaporReservoir};
use std::collections::{BTreeMap, BTreeSet};

/// One prepared supported endpoint bound to both original stored materials.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatingCellExchange {
    /// Complete post-heat substrate, including identity, geometry and boundary.
    pub expected_substrate: DepositedCell,
    /// Complete post-heat coating; a stale inventory, h or gas observation refuses.
    pub expected_coating: CoatingState,
    /// Exposure owned by the preceding physical interval, before endpoint events.
    pub interval_exposed_area_m2: f64,
    /// The substrate mass must remain unchanged; deposit mass may be exactly zero.
    /// Deposit h is derived from the original registered chart at this endpoint's
    /// substrate temperature, not supplied as an independent coordinate.
    pub replacement: SupportedSurfaceState,
    /// Prescribed pressure/temperature observation of the accepted original gas.
    /// Finite reservoir copies are refused, even for a single coated substrate.
    pub replacement_reservoir: VaporReservoir,
}

impl<'a> CoatingSimulation<'a> {
    /// Immutable original coating laws, including selected future substrates.
    /// Initial masses and reservoir observations describe registration, not the
    /// current inventories in `coatings()`.
    #[must_use]
    pub fn coating_surfaces(&self) -> &BTreeMap<CellId, CoatingSurface<'a>> {
        &self.surfaces
    }

    /// Original work ceilings. An external owner must carry one remaining
    /// substep allowance across every split of the requested physical advance.
    #[must_use]
    pub const fn work_budget(&self) -> DepositionBudget {
        self.thermal.work_budget()
    }

    /// Evaluate an active substrate's ORIGINAL caloric law without replacing
    /// its heat capacity, exposing mutable bare storage, or changing its state.
    pub fn temperature_at_specific_enthalpy_k(
        &self,
        id: CellId,
        specific_h_j_kg: f64,
    ) -> Result<f64, DepositionError> {
        self.thermal
            .temperature_at_specific_enthalpy_k(id, specific_h_j_kg)
    }

    /// Bounded active-cell selection with one original, ordered contact sweep.
    /// Capture these areas BEFORE heat and endpoint events for a split phase step.
    pub fn exposed_areas_m2(
        &self,
        selected: &[CellId],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        self.thermal.exposed_areas_m2(selected, checkpoint)
    }

    /// Advance native heat and events without any vapor/phase exchange.
    ///
    /// Sources, boundaries and links warm BOTH the original substrate and its
    /// actual coating at their common temperature. Births include the original
    /// substrate AND any supplied initial deposit in incoming mass and energy.
    /// `coating_config` controls this shared-temperature caloric solve only.
    ///
    /// Every configured and active reservoir must be prescribed; shared or
    /// private finite reservoirs would introduce an ambiguous second gas owner.
    /// Their observations are unchanged by this heat-only operation. An external
    /// owner solves phase exchange against its one gas and refreshes observations
    /// through `apply_coating_exchanges` at the accepted endpoint.
    ///
    /// The allowance cannot exceed the original per-advance budget. A zero
    /// allowance still permits due events and empty pre-birth intervals. Failure
    /// or cancellation returns no candidate and leaves this owner unchanged.
    pub fn advance_heat_to_limited(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        coating_config: SolverConfig,
        maximum_substeps: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), ThermalSublimationError> {
        self.require_external_observations(&mut checkpoint)?;
        let (next, audit) = self.advance_mode(
            target_s,
            thermal_config,
            coating_config,
            AdvanceMode::External { maximum_substeps },
            checkpoint,
        )?;
        Ok((next, audit.thermal))
    }

    /// Atomically admit externally solved supported endpoints on this history.
    ///
    /// The external owner admits gas/condensed mass and energy conservation and
    /// interval timing. Capture expected states after the heat advance and its
    /// endpoint events, while retaining START exposure and excluding endpoint
    /// births from their preceding molecular interval.
    ///
    /// Both expected states, immutable substrate mass, reference area, original
    /// caloric laws, pure-solid deposit domain and reservoir observations are
    /// validated before any write. Accepted endpoints preserve the clock,
    /// schedule, links, geometry, material identity and all unselected cells.
    pub fn apply_coating_exchanges(
        &mut self,
        exchanges: &[CoatingCellExchange],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        self.require_external_observations(&mut checkpoint)?;
        self.require_processed_events()?;
        if exchanges.len() > self.thermal.budget.max_cells {
            return Err(DepositionError::Budget("external coating exchange cells").into());
        }
        let mut ids = BTreeSet::new();
        let mut prepared = Vec::with_capacity(exchanges.len());
        for exchange in exchanges {
            checkpoint()?;
            let id = exchange.expected_substrate.id;
            let cell = self
                .thermal
                .cells
                .get(&id)
                .ok_or(DepositionError::InvalidInput(
                    "external coating exchange requires an active substrate",
                ))?;
            let coating = self.coatings.get(&id).ok_or(DepositionError::InvalidInput(
                "external coating exchange requires a registered active coating",
            ))?;
            if !ids.insert(id)
                || cell != &exchange.expected_substrate
                || coating != &exchange.expected_coating
            {
                return Err(DepositionError::InvalidInput(
                    "duplicate or stale external coating exchange",
                )
                .into());
            }
            let replacement = exchange.replacement;
            let area = exchange.interval_exposed_area_m2;
            if replacement.substrate_mass_kg != cell.mass_kg
                || !replacement.substrate_specific_enthalpy_j_kg.is_finite()
                || !nonnegative(replacement.deposit_mass_kg)
                || !nonnegative(area)
                || area > cell.geometry.surface_area_m2
                || replacement.area_m2 != area
            {
                return Err(DepositionError::InvalidInput(
                    "external coating endpoint changes substrate mass or leaves inventory/reference area",
                ).into());
            }
            require_prescribed(exchange.replacement_reservoir)?;
            let mut candidate = cell.clone();
            candidate.specific_h_j_kg = replacement.substrate_specific_enthalpy_j_kg;
            let mut model = self.surfaces[&id];
            model.initial_mass_kg = replacement.deposit_mass_kg;
            model.reservoir = exchange.replacement_reservoir;
            // Reuse the zero-time coating admission to derive h_d(T_s) and
            // validate the original pure-solid law, kinetics and observation.
            // No phase or heat transfer occurs in this constitutive check.
            let accepted = initial_state(&self.thermal, &candidate, &model, &mut checkpoint)?;
            self.thermal.materials[&cell.material_id].validate_cell(&candidate)?;
            prepared.push((id, candidate.specific_h_j_kg, accepted));
        }
        checkpoint()?;
        // No fallible work, constitutive callbacks or cancellation after this gate.
        for (id, substrate_h, coating) in prepared {
            self.thermal.cells.get_mut(&id).unwrap().specific_h_j_kg = substrate_h;
            *self.coatings.get_mut(&id).unwrap() = coating;
        }
        Ok(())
    }

    fn require_external_observations(
        &self,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        checkpoint()?;
        if self.chamber.is_some() {
            return Err(DepositionError::InvalidInput(
                "external coating exchange cannot own a shared finite reservoir",
            )
            .into());
        }
        // Include future births so a later event cannot silently install a
        // private finite stock into an otherwise externally coupled history.
        for (index, reservoir) in self
            .surfaces
            .values()
            .map(|surface| surface.reservoir)
            .chain(self.coatings.values().map(|coating| coating.reservoir))
            .enumerate()
        {
            if index % 256 == 0 {
                checkpoint()?;
            }
            require_prescribed(reservoir)?;
        }
        Ok(())
    }
}

fn require_prescribed(reservoir: VaporReservoir) -> Result<(), DepositionError> {
    if let VaporReservoir::Prescribed {
        partial_pressure_pa,
        temperature_k,
    } = reservoir
    {
        if nonnegative(partial_pressure_pa) && crate::enthalpy::positive(temperature_k) {
            return Ok(());
        }
    }
    Err(DepositionError::InvalidInput(
        "external coating exchange requires prescribed gas observations",
    ))
}

#[cfg(test)]
#[path = "external/tests.rs"]
mod tests;
