//! Multiport flow and independent motion retaining each original solid specimen.
//!
//! No numeric piston-mass override, new material resolver or reconstructed port
//! history. The same moving-network implementation advances the existing inner
//! dynamics; this boundary retains the immutable regional body/source binding.

use std::collections::BTreeMap;
use super::{ClosedGasParcel, GasPiston, MaterialGasPiston, ResolvedBody};
use super::super::exchange::graph;
pub use graph::{MovingGasCandidate, MovingGasDrive, MovingGasError, MovingGasSolve, MovingGasStep};
pub use graph::reactive::{ChamberReaction, ReactiveMovingGasCandidate};
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
    let (numeric, mut bodies) = prepare(pistons, reservoirs, links, &mut cancelled)?;
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

// Both entry points preserve the same accepted numeric owner and solid binding.
#[allow(clippy::type_complexity)]
fn prepare(
    pistons: &BTreeMap<u64, MaterialGasPiston>, reservoirs: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink], cancelled: &mut impl FnMut() -> bool,
) -> Result<(BTreeMap<u64, GasPiston>, BTreeMap<u64, ResolvedBody>), MovingGasError> {
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
    Ok((numeric, bodies))
}

/// Flow, chemistry and motion with original material-derived wall inertia.
/// Incoming mixtures react once before pressure/momentum trials; all chemical
/// catalogs, consumed-species limits and full-network conservation gates belong
/// to the shared owners. The existing regional body/source records are retained
/// exactly. There is no mass override or second chemical/thermal/mechanical law.
///
/// `reactions` selects existing moving or fixed gas ids. Gas chemistry does not
/// implicitly heat, expand, melt or re-resolve the rigid wall material. Changing
/// density data is a model-initialization choice, not a live inertia exchange.
/// All source records, gas states, momenta and one-interval accounts publish
/// together. Empty bindings retain the original nonreacting material path.
#[allow(clippy::too_many_arguments)]
pub fn reactive_step_candidate(
    pistons: &BTreeMap<u64, MaterialGasPiston>, reservoirs: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, drives: &BTreeMap<u64, MovingGasDrive>,
    reactions: &[ChamberReaction<'_>], controls: MovingGasSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReactiveMovingGasCandidate<MaterialGasPiston>, MovingGasError> {
    let (numeric, mut bodies) = prepare(pistons, reservoirs, links, &mut cancelled)?;
    let candidate = graph::reactive::step_candidate(&numeric, reservoirs, links, duration_s,
        endpoints, drives, reactions, controls, &mut cancelled)?;
    let candidate = candidate.map_pistons(|id, piston| MaterialGasPiston {
        body: bodies.remove(&id).expect("moving network preserves validated material identities"), piston,
    });
    if cancelled() { return Err(MovingGasError::Cancelled); }
    Ok(candidate)
}

#[cfg(test)]
mod tests;
