//! Several partially frozen bodies sharing one original finite atmosphere.
//!
//! Every explicitly exposed LIQUID interface uses its own original sharp
//! solid/liquid caloric source. All films see the same initial gas, and the
//! parent chamber accumulates their signed species and carried energy before
//! reconstructing that gas once. Concurrent evaporation cannot finance an
//! excessive withdrawal of old vapor at another surface.
//!
//! The original numeric constituent evaluator owns source compatibility,
//! liquid density and actual liquid donor inventory. Bulk enthalpy owns both
//! sensible and fusion storage. Per-body fusion reports are diagnostics of
//! those same endpoints, never additional heat sources or inferred interface
//! areas. Each body retains a liquid seed. Complete freezing, nucleation,
//! moving geometry and reciprocal gas heat require their respective owners.
//! This is the existing frozen-rate reduction and needs timestep refinement.

use super::{
    ChamberCaloric, ChamberSurface, ClosedGasParcel, EquilibriumEvaporation, EvaporationLaw,
    GasEvaporationError, GasLiquidBinding, GasLiquidChamberCandidate, GasLiquidChamberSolve,
    GasLiquidChamberStep, GasLiquidChamberSurfaceStep, GasSpeciesBoundaryCandidate,
    LiquidSurfaceState, PhaseInventoryStep, ThermochemicalGasModel, advance_chamber,
};

/// One original condensed inventory and its explicitly exposed liquid area.
/// Stable IDs identify surfaces, not material sources; multiple bodies may
/// share a source while retaining separate mass and enthalpy histories.
#[derive(Clone, Copy)]
pub struct EquilibriumChamberSurface<'a> {
    pub id: u64,
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    /// TOTAL condensed mass and original bulk h; area is liquid exposure.
    pub bulk: LiquidSurfaceState,
    pub inventory: EquilibriumEvaporation<'a>,
    pub mass_transfer_coefficient_m_s: f64,
}

impl<'a> From<EquilibriumChamberSurface<'a>> for ChamberSurface<'a> {
    fn from(surface: EquilibriumChamberSurface<'a>) -> Self {
        Self {
            id: surface.id,
            law: surface.law,
            binding: surface.binding,
            liquid: surface.bulk,
            caloric: ChamberCaloric::Equilibrium(surface.inventory),
            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
        }
    }
}

/// Original-source fusion inventory at the same accepted chamber endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumChamberPhaseStep {
    pub id: u64,
    pub phase: PhaseInventoryStep,
}

/// One original chamber endpoint and its separately stored phase partitions.
/// Pure-liquid users keep their original report layout and allocations. This
/// lane allocates one contiguous partition vector, without per-body boxes.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumChamberCandidate<'catalog> {
    exchange: GasLiquidChamberCandidate<'catalog>,
    phases: Vec<EquilibriumChamberPhaseStep>,
}

impl<'catalog> EquilibriumChamberCandidate<'catalog> {
    pub const fn exchange(&self) -> &GasLiquidChamberCandidate<'catalog> {
        &self.exchange
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.exchange.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.exchange.catalog()
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        self.exchange.gas_exchange()
    }
    /// Canonical exchange rows. Legacy `liquid` fields describe TOTAL condensed
    /// mass and bulk enthalpy, not just the liquid constituent.
    pub fn surfaces(&self) -> &[GasLiquidChamberSurfaceStep] {
        self.exchange.surfaces()
    }
    pub const fn step(&self) -> &GasLiquidChamberStep {
        self.exchange.step()
    }
    pub fn phases(&self) -> &[EquilibriumChamberPhaseStep] {
        &self.phases
    }
    /// Allocation-free lookup in the same canonical surface-ID order.
    pub fn phase_inventory(&self, id: u64) -> Option<&PhaseInventoryStep> {
        self.phases
            .binary_search_by_key(&id, |row| row.id)
            .ok()
            .map(|index| &self.phases[index].phase)
    }
    pub fn into_parts(
        self,
    ) -> (
        GasLiquidChamberCandidate<'catalog>,
        Vec<EquilibriumChamberPhaseStep>,
    ) {
        (self.exchange, self.phases)
    }
}

/// Advance all original fusion inventories and their one finite gas together.
///
/// Reports are in canonical surface-ID order. The phase rows contain each
/// complete original-source partition; the exchange rows' legacy `liquid`
/// fields contain TOTAL condensed storage. The original one-chamber mass and
/// sum-of-absolute-energy budgets do not grow with the number of surfaces.
/// Empty rosters admit/preserve the gas, and cancellation returns no endpoint.
pub fn advance_chamber_with_gas<'catalog>(
    surfaces: &[EquilibriumChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumChamberCandidate<'catalog>, GasEvaporationError> {
    advance_chamber_with_gas_mass_selection(
        surfaces,
        gas,
        catalog,
        duration_s,
        controls,
        |_, mass| Ok(mass),
        cancelled,
    )
}

/// The same shared phase owner with one native mass selector per stable ID.
///
/// Requested and selected debits retain the original actual-liquid and old-gas
/// donor limits. Selection runs in canonical order and may only choose a nearby
/// representable signed increment under the existing absolute mass budget.
/// External callback effects are caller-owned and cannot be rolled back.
#[allow(clippy::too_many_arguments)]
pub fn advance_chamber_with_gas_mass_selection<'catalog>(
    surfaces: &[EquilibriumChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasLiquidChamberSolve,
    select_mass: impl FnMut(u64, f64) -> Result<f64, GasEvaporationError>,
    cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumChamberCandidate<'catalog>, GasEvaporationError> {
    advance_chamber(
        surfaces.iter().copied().map(ChamberSurface::from),
        gas,
        catalog,
        duration_s,
        controls,
        select_mass,
        cancelled,
    )
    .map(|(exchange, phases)| EquilibriumChamberCandidate { exchange, phases })
}

#[cfg(test)]
mod tests;
