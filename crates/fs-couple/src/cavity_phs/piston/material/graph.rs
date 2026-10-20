//! Multiport flow and independent motion retaining each original solid specimen.
//!
//! No numeric piston-mass override, new material resolver or reconstructed port
//! history. The same moving-network implementation advances the existing inner
//! dynamics; this boundary retains the immutable regional body/source binding.

use std::collections::BTreeMap;
use super::{ClosedGasParcel, MaterialGasPiston};
use super::super::exchange::graph;
pub use graph::{MovingGasCandidate, MovingGasDrive, MovingGasError, MovingGasSolve, MovingGasStep};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;

/// Run the original shared-state flow/momentum composition with material inertia.
///
/// Each `MaterialGasPiston` retains its original resolved regional body, canonical
/// mass, source conditions and accepted momentum/port history. The returned
/// bodies are not rebuilt from density names, heated by gas temperature or given
/// new inertias after the step. Material replacement is an initialization choice,
/// not an unaccounted live mass/energy exchange. All gas transport and mechanics
/// gates, first-order splitting limits and open-port limitations are unchanged.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    pistons: &BTreeMap<u64, MaterialGasPiston>, reservoirs: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, drives: &BTreeMap<u64, MovingGasDrive>,
    controls: MovingGasSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<MovingGasCandidate<MaterialGasPiston>, MovingGasError> {
    if cancelled() { return Err(MovingGasError::Cancelled); }
    if pistons.len().checked_add(reservoirs.len()).is_none_or(|n| n > MAX_SNAPSHOT_PARCELS)
        || links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(MovingGasError::Invalid("material moving graph exceeds the shared work bounds"));
    }
    let mut numeric = BTreeMap::new();
    let mut bodies = BTreeMap::new();
    for (&id, source) in pistons {
        if cancelled() { return Err(MovingGasError::Cancelled); }
        numeric.insert(id, source.piston.clone());
        bodies.insert(id, source.body.clone());
    }
    let candidate = graph::step_candidate(&numeric, reservoirs, links, duration_s,
        endpoints, drives, controls, &mut cancelled)?;
    // The shared owner preserves precisely the moving-id set. Move each saved
    // body rather than cloning its potentially large regional source records again.
    let candidate = candidate.map_pistons(|id, piston| MaterialGasPiston {
        body: bodies.remove(&id).expect("moving network preserves validated material identities"), piston,
    });
    if cancelled() { return Err(MovingGasError::Cancelled); }
    Ok(candidate)
}

#[cfg(test)]
mod tests;
