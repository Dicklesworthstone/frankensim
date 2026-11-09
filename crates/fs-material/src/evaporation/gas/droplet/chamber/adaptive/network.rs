//! Whole-network accuracy and physical accounts for evolving droplet chambers.
//!
//! Gas transport and evaporation retain their original constitutive and caloric
//! owners. These helpers compare complete physical states and accumulate only
//! explicitly accepted original intervals. Every gas parcel participates, even
//! if no droplet occupies it. Fixed volume, stable IDs and positive droplets
//! remain required; no gas entropy statement is promoted to the coupled liquid
//! system, and no absolute formation-energy subtraction defines a defect.

mod account;
pub use account::{
    GasDropletNetworkAudit, GasDropletNetworkParcelAudit,
    GasDropletNetworkParcelSpeciesAudit, GasDropletNetworkWindow,
};

use std::collections::{BTreeMap, BTreeSet};

use super::{
    AdaptiveGasDropletControls, ClosedGasParcel, GasDropletSeed, GasDropletSurfaceStep,
    GasDropletTimeError, GasEvaporationError, SphericalDropletGeometry, component_error, poll,
};

/// Explicit storage/work ceilings for the complete network, not one chamber.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GasDropletNetworkLimits {
    pub maximum_parcels: usize,
    /// Active species entries across one complete gas map, and the retained
    /// per-parcel account. New arrivals consume this same account allowance.
    pub maximum_species_entries: usize,
}
impl Default for GasDropletNetworkLimits {
    fn default() -> Self {
        Self {
            maximum_parcels: 1024,
            maximum_species_entries: 65536,
        }
    }
}

/// Bound the original roster before cloning maps or allocating inventory work.
/// This inspects existing model metadata without querying caloric/source laws.
pub fn admit_work(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    limits: GasDropletNetworkLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), GasEvaporationError> {
    poll(cancelled)?;
    if limits.maximum_parcels == 0
        || limits.maximum_species_entries == 0
        || parcels.is_empty()
        || parcels.len() > limits.maximum_parcels
    {
        return Err(invalid("droplet network exceeds its nonempty parcel work ceiling"));
    }
    let mut entries = 0usize;
    for gas in parcels.values() {
        poll(cancelled)?;
        entries = entries
            .checked_add(gas.model().thermodynamics().components().len())
            .ok_or_else(|| invalid("droplet network species work count overflow"))?;
        if entries > limits.maximum_species_entries {
            return Err(invalid("droplet network exceeds its species work ceiling"));
        }
    }
    Ok(())
}

/// First-order coarse/fine discrepancy on all liquid and gas coordinates.
///
/// Liquid rows and seeds are in canonical stable-ID order. The original coarse
/// interval supplies the initial liquid temperature and geometry, so this helper
/// never queries a source before the actual transport-then-phase owner admits it.
/// Gas comparisons include total mass and every species in the union of all
/// three states, including a species that was initially absent at a receiving
/// node. The gas-temperature field reports the maximum across all nodes.
///
/// Relative scales are changes from the original comparison state, not absolute
/// enthalpy or Kelvin temperature. The caller supplies the original method-order
/// correction for its actual fine spans. This is a local numerical estimate.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn time_error(
    initial: &[GasDropletSeed<'_>],
    before: &BTreeMap<u64, ClosedGasParcel>,
    coarse: (&[GasDropletSurfaceStep], &BTreeMap<u64, ClosedGasParcel>),
    fine: (&[GasDropletSurfaceStep], &BTreeMap<u64, ClosedGasParcel>),
    controls: AdaptiveGasDropletControls,
    correction: f64,
    limits: GasDropletNetworkLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasDropletTimeError, GasEvaporationError> {
    poll(cancelled)?;
    controls.validate(0.0, initial.len())?;
    admit_work(before, limits, cancelled)?;
    admit_work(coarse.1, limits, cancelled)?;
    admit_work(fine.1, limits, cancelled)?;
    if !correction.is_finite()
        || correction <= 0.0
        || initial.len() != coarse.0.len()
        || initial.len() != fine.0.len()
        || initial.windows(2).any(|rows| rows[0].id >= rows[1].id)
        || !before.keys().eq(coarse.1.keys())
        || !before.keys().eq(fine.1.keys())
    {
        return Err(invalid("droplet network time comparison changed its physical roster"));
    }
    let mut error = GasDropletTimeError::default();
    for ((old, coarse), fine) in initial.iter().zip(coarse.0).zip(fine.0) {
        poll(cancelled)?;
        let geometry = SphericalDropletGeometry::from_mass(
            old.droplet.mass_kg,
            coarse.geometry[0].liquid_density_kg_m3(),
        )?;
        if old.id != coarse.id || old.id != fine.id || geometry != coarse.geometry[0] {
            return Err(invalid("droplet network comparison changed an input ID or mass"));
        }
        let (difference, ratio) = component_error(
            old.droplet.specific_enthalpy_j_kg,
            coarse.droplet.specific_enthalpy_j_kg,
            fine.droplet.specific_enthalpy_j_kg,
            controls.enthalpy_absolute_tolerance_j_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_enthalpy_difference_j_kg =
            error.maximum_enthalpy_difference_j_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            coarse.exchange.surface_temperature_k[0],
            coarse.exchange.surface_temperature_k[1],
            fine.exchange.surface_temperature_k[1],
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_temperature_difference_k =
            error.maximum_temperature_difference_k.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            old.droplet.mass_kg,
            coarse.droplet.mass_kg,
            fine.droplet.mass_kg,
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_liquid_mass_difference_kg =
            error.maximum_liquid_mass_difference_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (difference, ratio) = component_error(
            geometry.radius_m(),
            coarse.geometry[1].radius_m(),
            fine.geometry[1].radius_m(),
            controls.radius_absolute_tolerance_m,
            controls.relative_tolerance,
            correction,
        )?;
        error.maximum_radius_difference_m = error.maximum_radius_difference_m.max(difference);
        error.normalized = error.normalized.max(ratio);
    }
    for (&id, old) in before {
        poll(cancelled)?;
        let coarse = &coarse.1[&id];
        let fine = &fine.1[&id];
        if coarse.volume_m3().to_bits() != old.volume_m3().to_bits()
            || fine.volume_m3().to_bits() != old.volume_m3().to_bits()
        {
            return Err(invalid("droplet network time comparison requires fixed gas volumes"));
        }
        let (difference, ratio) = component_error(
            old.state().state().temperature,
            coarse.state().state().temperature,
            fine.state().state().temperature,
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
            correction,
        )?;
        error.gas_temperature_difference_k = error.gas_temperature_difference_k.max(difference);
        error.normalized = error.normalized.max(ratio);
        let (_, ratio) = component_error(
            old.mass_kg(),
            coarse.mass_kg(),
            fine.mass_kg(),
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
            correction,
        )?;
        error.normalized = error.normalized.max(ratio);
        let old_inventory = old.species_inventory_kg()?;
        let coarse_inventory = coarse.species_inventory_kg()?;
        let fine_inventory = fine.species_inventory_kg()?;
        let species: BTreeSet<_> = old_inventory
            .keys()
            .chain(coarse_inventory.keys())
            .chain(fine_inventory.keys())
            .collect();
        for id in species {
            poll(cancelled)?;
            let (difference, ratio) = component_error(
                old_inventory.get(id).copied().unwrap_or(0.0),
                coarse_inventory.get(id).copied().unwrap_or(0.0),
                fine_inventory.get(id).copied().unwrap_or(0.0),
                controls.mass_absolute_tolerance_kg,
                controls.relative_tolerance,
                correction,
            )?;
            error.maximum_species_mass_difference_kg =
                error.maximum_species_mass_difference_kg.max(difference);
            error.normalized = error.normalized.max(ratio);
        }
    }
    poll(cancelled)?;
    Ok(error)
}

fn invalid(message: &'static str) -> GasEvaporationError {
    GasEvaporationError::Invalid(message)
}

#[cfg(test)]
mod tests;
