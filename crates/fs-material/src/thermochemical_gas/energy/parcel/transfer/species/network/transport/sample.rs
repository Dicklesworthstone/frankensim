//! Sample physical gas paths without reconstructing any inventory. The closed
//! gas driver and joint external-phase driver consume these SAME rates.
use super::*;
use super::super::SpeciesNetworkHeatTransfer;

/// Source-bearing mechanism rates. No closed-gas entropy claim is attached:
/// an enclosing phase coupling may legitimately remove gas entropy.
#[derive(Clone, Debug, PartialEq)]
pub struct GasTransportRates {
    original_parcels: Vec<u64>,
    pub duration_s: f64,
    pub nozzles: Vec<SnapshotNozzleLinkStep>,
    pub diffusion: Vec<DiffusionNetworkLinkStep>,
    pub porous: Vec<PorousNetworkLinkStep>,
    pub membranes: Vec<MembraneNetworkLinkStep>,
    pub dusty_gas: Vec<DustyGasNetworkLinkStep>,
    pub effusion: Vec<EffusionNetworkLinkStep>,
}
impl GasTransportRates {
    /// Check the final coupled endpoint against the original path directions.
    /// This neither resamples a rate nor spends any inventory a second time.
    pub fn admit_endpoint_driving_forces(&self, parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool) -> Result<(), GasTransportError> {
        check(&mut cancelled)?;
        if !parcels.keys().copied().eq(self.original_parcels.iter().copied())
            || self.nozzles.iter().any(|row| !parcels.contains_key(&row.donor)
                || !parcels.contains_key(&row.receiver))
            || self.porous.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
            || self.membranes.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
            || self.dusty_gas.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
            || self.effusion.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id))) {
            return Err(Error::Invalid("final coupled transport endpoint needs every original gas identity"));
        }
        pressure::admit_pressures(parcels, &self.nozzles, &mut cancelled)?;
        porous::admit_pressures(parcels, &self.porous, &mut cancelled)?;
        membrane::admit_partial_pressures(parcels, &self.membranes, &mut cancelled)?;
        check(&mut cancelled)
    }
    pub(super) fn into_closed_step(self, exchange: SpeciesNetworkStep, entropy_change_j_k: f64) -> GasTransportStep {
        GasTransportStep { duration_s:self.duration_s, nozzles:self.nozzles, diffusion:self.diffusion,
            porous:self.porous, membranes:self.membranes, dusty_gas:self.dusty_gas,
            effusion:self.effusion, exchange, entropy_change_j_k }
    }
}

pub(crate) struct GasTransportSample {
    pub transfers: Vec<SpeciesNetworkTransfer>,
    pub heat: Vec<SpeciesNetworkHeatTransfer>,
    pub rates: GasTransportRates,
}

pub(crate) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>, paths: GasTransportPaths<'_>,
    duration_s: f64, endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<GasTransportSample, Error> {
    check(&mut cancelled)?;
    let count = [paths.nozzles.len(), paths.diffusion.len(), paths.porous.len(),
        paths.membranes.len(), paths.dusty_gas.len(), paths.effusion.len()]
        .into_iter().try_fold(0usize, |total, n| total.checked_add(n));
    if count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS) {
        return Err(Error::Invalid("combined connection count exceeds the shared transport bound"));
    }
    let (bulk, nozzle_reports) = pressure::sample_transfers(
        parcels, paths.nozzles, duration_s, endpoints, &mut cancelled)?;
    let (mut transfers, diffusion_reports) = molecular::sample_transfers(
        parcels, paths.diffusion, duration_s, endpoints,
        DiffusionNetworkSolve { exchange: controls.exchange,
            entropy_tolerance_j_k: controls.entropy_tolerance_j_k }, &mut cancelled)?;
    let (porous_bulk, porous_reports) = porous::sample_transfers(
        parcels, paths.porous, duration_s, &mut cancelled)?;
    let (selective, membrane_reports) = membrane::sample_transfers(
        parcels, paths.membranes, duration_s, &mut cancelled)?;
    let (pore_species, dusty_gas_reports) = dusty_gas::sample_transfers(
        parcels, paths.dusty_gas, duration_s, &mut cancelled)?;
    let (molecular_beams, molecular_energy, effusion_reports) = effusion::sample_transfers(
        parcels, paths.effusion, duration_s, &mut cancelled)?;
    if [transfers.len(), selective.len(), pore_species.len(), molecular_beams.len(), molecular_energy.len()]
        .into_iter().try_fold(0usize, |total, n| total.checked_add(n))
        .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS) {
        return Err(Error::Invalid("combined selective flux count exceeds the shared work bound"));
    }
    transfers.extend(selective);
    transfers.extend(pore_species);
    transfers.extend(molecular_beams);
    for row in bulk.into_iter().chain(porous_bulk) {
        check(&mut cancelled)?;
        if row.mass_kg == 0.0 { continue; }
        let parts = inventory(parcels[&row.donor].model(), row.mass_kg)
            .map_err(|error| super::super::material(row.donor, error))?;
        for (species, mass_kg) in parts {
            check(&mut cancelled)?;
            if transfers.len() + molecular_energy.len() == MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid("combined species flux count exceeds the shared work bound"));
            }
            transfers.push(SpeciesNetworkTransfer { donor:row.donor, receiver:row.receiver, species, mass_kg });
        }
    }
    check(&mut cancelled)?;
    Ok(GasTransportSample { transfers, heat:molecular_energy, rates:GasTransportRates {
        original_parcels:parcels.keys().copied().collect(), duration_s, nozzles:nozzle_reports,
        diffusion:diffusion_reports, porous:porous_reports, membranes:membrane_reports,
        dusty_gas:dusty_gas_reports, effusion:effusion_reports } })
}
