//! Adaptive sharp-phase inventories on the ORIGINAL chamber history owner.
//!
//! Temperature alone cannot detect error on a fusion plateau. Compare actual
//! SOLID and LIQUID masses from the original equilibrium chart in addition to
//! the existing temperature, total condensed mass and gas-species error norm.
//! Each phase receives its own absolute/relative mass scale, so a trace solid
//! population cannot hide behind the complete condensed inventory.
//!
//! All rates, old-inventory limits, gas reconstruction, accepted-stream sums,
//! cancellation and whole-window energy/species gates remain with the existing
//! owners. Source charts and accommodation maps are borrowed across trials.
//! This supplies a local first-order error estimate, not a rigorous global
//! bound, phase-diagram validation or automatic melting/extinction topology.

use super::{
    ClosedGasParcel, GasChamberAdvanceCandidate, GasChamberAdvanceError,
    GasChamberAdvanceSolve, GasChamberCandidate, GasChamberSurface,
    SharpMeltingCaloric, ThermalAccommodation, ThermochemicalGasModel,
    advance_history_with_calorics, error_ratio, phase_inventory, poll,
};
use phase_inventory::{CondensedInventory, CondensedPhaseStep};
use std::collections::BTreeMap;

// Canonical surface ordering is supplied by the original adaptive driver.
// No new caloric interpolation or constant heat capacity is evaluated here.
pub(super) fn compare(
    before: &[GasChamberSurface<'_>],
    calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    coarse: &GasChamberCandidate<'_>,
    fine: &GasChamberCandidate<'_>,
    controls: GasChamberAdvanceSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, GasChamberAdvanceError> {
    let mut maximum = 0.0_f64;
    for ((original, a), b) in before.iter().zip(coarse.surfaces()).zip(fine.surfaces()) {
        poll(cancelled)?;
        let source = calorics[&original.id];
        let initial = CondensedInventory::new(source, original.solid)?;
        let coarse = CondensedInventory::new(source, a.solid)?;
        let fine = CondensedInventory::new(source, b.solid)?;
        for (old, a, b) in [
            (initial.solid_mass_kg, coarse.solid_mass_kg, fine.solid_mass_kg),
            (initial.liquid_mass_kg, coarse.liquid_mass_kg, fine.liquid_mass_kg),
        ] {
            maximum = maximum.max(error_ratio(a, b, old,
                controls.absolute_mass_tolerance_kg, controls.relative_tolerance)?);
        }
    }
    Ok(maximum)
}

/// Original accepted history plus NET melting/freezing across its full duration.
/// Phase reports compare initial and final inventories, never the last half step.
#[derive(Clone, Debug, PartialEq)]
pub struct PartitionedChamberAdvanceCandidate<'catalog> {
    history: GasChamberAdvanceCandidate<'catalog>,
    phases: BTreeMap<u64, CondensedPhaseStep>,
}
impl<'catalog> PartitionedChamberAdvanceCandidate<'catalog> {
    /// Complete original stream, conservation, error and work accounting.
    /// Condensed-body storage includes fusion; it must not be added again.
    #[must_use]
    pub const fn history(&self) -> &GasChamberAdvanceCandidate<'catalog> { &self.history }
    /// Original gas with its final composition, calorics, EOS and transport.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { self.history.gas() }
    /// Actual constituent changes in canonical surface-id order.
    /// melted_mass_kg is positive for melting and negative for freezing.
    #[must_use]
    pub fn phases(&self) -> &BTreeMap<u64, CondensedPhaseStep> { &self.phases }
    /// Consume all accepted endpoints and their original-source phase accounts.
    #[must_use]
    pub fn into_parts(self) -> (GasChamberAdvanceCandidate<'catalog>, BTreeMap<u64, CondensedPhaseStep>) {
        (self.history, self.phases)
    }
}

/// Advance mixed solid/liquid bodies and a shared gas for one complete duration.
///
/// Every surface retains its explicit solid-vapor law, species binding and
/// original sharp melting chart. Callbacks must evaluate those same charts.
/// The actual exposed SOLID area is held fixed; no fraction-to-area conversion,
/// moving pure-body geometry or liquid evaporation is inferred. Optional
/// accommodation uses the same area and old gas as the sticking streams.
///
/// Both condensed phase masses participate in adaptive error control. Gross
/// sublimation cannot spend old liquid mass or material frozen later in the
/// interval. Accepted fine steps alone contribute to stream/reflection history.
/// Failed or cancelled trials publish no state; inventory errors alone retry.
/// Exact all-liquid completion may be a final endpoint, but a subsequent solid
/// interface interval refuses until an external topology/evaporation owner acts.
#[allow(clippy::too_many_arguments)]
pub fn advance_partitioned_chamber_adaptive<'catalog>(
    surfaces: &[GasChamberSurface<'_>],
    calorics: &BTreeMap<u64, SharpMeltingCaloric<'_>>,
    accommodation: Option<&BTreeMap<u64, ThermalAccommodation>>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration_s: f64,
    controls: GasChamberAdvanceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PartitionedChamberAdvanceCandidate<'catalog>, GasChamberAdvanceError> {
    let history = advance_history_with_calorics(surfaces, accommodation, None, Some(calorics),
        gas, catalog, duration_s, controls, &mut cancelled)?;
    let phases = phase_inventory::phase_reports(surfaces, calorics, history.endpoint(),
        controls.step.surface.gas.mass_tolerance_kg, &mut cancelled)?;
    poll(&mut cancelled)?;
    Ok(PartitionedChamberAdvanceCandidate { history, phases })
}

#[cfg(test)]
mod tests;
