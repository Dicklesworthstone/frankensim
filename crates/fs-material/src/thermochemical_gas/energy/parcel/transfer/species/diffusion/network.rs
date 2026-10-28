//! One source snapshot for all diffusion faces, one update per finite volume.
//! Reuse the original Maxwell-Stefan face solver and selective inventory owner.
use std::collections::BTreeMap;
use super::{BinaryDiffusionMaterial, ClosedGasParcel, GasDiffusionFlux, GasDiffusionGeometry,
    GasDiffusionSolve, GasEnergyError, GasSpeciesExchangeSolve, positive};
use super::super::network::{self, SpeciesNetworkError, SpeciesNetworkSolve,
    SpeciesNetworkStep, SpeciesNetworkTransfer};
use super::super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use crate::thermochemical_gas::energy::GasEnergySolve;

/// Face-law admission only; endpoint/withdrawal limits belong to the whole graph.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionFaceSolve {
    pub maximum_temperature_difference_k: f64,
    pub maximum_pressure_difference_pa: f64,
    pub gradient_tolerance_per_m: f64,
    pub molar_flux_tolerance_mol_m2_s: f64,
}
/// One oriented face, with actual unordered-pair material bindings.
#[derive(Clone, Debug)]
pub struct DiffusionNetworkLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub geometry: GasDiffusionGeometry,
    pub pairs: &'a [BinaryDiffusionMaterial<'a>],
    pub solve: DiffusionFaceSolve,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffusionNetworkSolve {
    pub exchange: SpeciesNetworkSolve,
    /// Maximum rounding decrease in the ENTIRE closed network's entropy [J/K].
    pub entropy_tolerance_j_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum DiffusionNetworkError {
    Invalid(&'static str),
    Face { id: u64, error: GasEnergyError },
    Exchange(SpeciesNetworkError),
    Parcel { id: u64, error: GasEnergyError },
    Cancelled,
}
impl core::fmt::Display for DiffusionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"material diffusion network: {self:?}")
    }
}
impl std::error::Error for DiffusionNetworkError {}
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkLinkStep {
    pub id: u64,
    /// Left/right orientation used by this original face law.
    pub parcels: [u64; 2],
    pub flux: GasDiffusionFlux,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkStep {
    pub duration_s: f64,
    pub links: Vec<DiffusionNetworkLinkStep>,
    pub exchange: SpeciesNetworkStep,
    /// Conserved-species entropy change; numerical species-mass residuals
    /// remain in the separate exchange account.
    pub entropy_change_j_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DiffusionNetworkCandidate {
    parcels: BTreeMap<u64,ClosedGasParcel>,
    step: DiffusionNetworkStep,
}
impl DiffusionNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64,ClosedGasParcel> { &self.parcels }
    pub const fn step(&self) -> &DiffusionNetworkStep { &self.step }
    pub fn into_parts(self) -> (BTreeMap<u64,ClosedGasParcel>,DiffusionNetworkStep) {
        (self.parcels,self.step)
    }
}
use DiffusionNetworkError as Error;
fn check(cancelled:&mut impl FnMut()->bool) -> Result<(),Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

/// Resolve all molecular diffusion faces, then update the common inventories.
///
/// Changing link IDs or input order cannot change the endpoint arithmetic.
/// Report order follows IDs, but the inventory owner sums canonical physical
/// species transfers. Incoming species cannot replenish an outgoing allowance
/// or travel across a second face within the same explicit physical interval.
/// Pairwise-safe flows may jointly exceed a donor's limit: the whole graph then
/// refuses. No clipping, per-edge remixture, or additional enthalpy deposition.
///
/// All original isothermal/isobaric, ideal-gas, molar-frame and source-domain
/// restrictions remain. The entire network's final entropy is checked; this
/// is not a proof of unconditional stability or spatial/time accuracy. Closed
/// faces do not resolve fictitious pair coefficients. Refine time and space.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels:&BTreeMap<u64,ClosedGasParcel>, links:&[DiffusionNetworkLink<'_>], duration_s:f64,
    endpoints:&BTreeMap<u64,GasEnergySolve>, controls:DiffusionNetworkSolve,
    mut cancelled:impl FnMut()->bool,
) -> Result<DiffusionNetworkCandidate,Error> {
    let (transfers, reports) = sample_transfers(parcels, links, duration_s, endpoints, controls, &mut cancelled)?;
    let candidate=network::exchange_candidate(parcels,&transfers,endpoints,controls.exchange,&mut cancelled)
        .map_err(|error|if error==SpeciesNetworkError::Cancelled {Error::Cancelled} else {Error::Exchange(error)})?;
    let entropy_change_j_k = admit_entropy(parcels, candidate.parcels(),
        controls.entropy_tolerance_j_k, controls.exchange.mass_tolerance_kg, &mut cancelled)?;
    check(&mut cancelled)?;
    let (parcels,exchange)=candidate.into_parts();
    Ok(DiffusionNetworkCandidate {parcels,step:DiffusionNetworkStep {duration_s,links:reports,
        exchange,entropy_change_j_k}})
}

// One face sampler serves both diffusion-only and mixed-path transport. It
// resolves original sources without spending any accepted species inventory.
pub(crate) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>, links: &[DiffusionNetworkLink<'_>], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, controls: DiffusionNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<SpeciesNetworkTransfer>, Vec<DiffusionNetworkLinkStep>), Error> {
    check(&mut cancelled)?;
    if !positive(duration_s) || !positive(controls.entropy_tolerance_j_k) || parcels.is_empty()
        || parcels.len()>MAX_SNAPSHOT_PARCELS || links.len()>MAX_SNAPSHOT_TRANSFERS
        || !parcels.keys().eq(endpoints.keys()) {
        return Err(Error::Invalid("bounded nonempty graph, exact endpoint map and positive time/entropy allowance required"));
    }
    let mut ordered=BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0]==link.parcels[1] || link.parcels.iter().any(|id|!parcels.contains_key(id))
            || ordered.insert(link.id,link).is_some() {
            return Err(Error::Invalid("unique face IDs and two distinct existing gas owners required"));
        }
    }
    let mut transfers=Vec::new(); let mut reports=Vec::with_capacity(links.len());
    for (&id,link) in &ordered {
        check(&mut cancelled)?;
        let [left,right]=link.parcels;
        // Per-node inversion controls and aggregate allowances have ONE owner.
        let solve=GasDiffusionSolve {
            maximum_temperature_difference_k:link.solve.maximum_temperature_difference_k,
            maximum_pressure_difference_pa:link.solve.maximum_pressure_difference_pa,
            gradient_tolerance_per_m:link.solve.gradient_tolerance_per_m,
            molar_flux_tolerance_mol_m2_s:link.solve.molar_flux_tolerance_mol_m2_s,
            entropy_tolerance_j_k:controls.entropy_tolerance_j_k,
            exchange:GasSpeciesExchangeSolve { gas:[endpoints[&left],endpoints[&right]],
                maximum_species_fraction:controls.exchange.maximum_species_fraction,
                mass_tolerance_kg:controls.exchange.mass_tolerance_kg,
                energy_tolerance_j:controls.exchange.energy_tolerance_j },
        };
        let flux=parcels[&left].diffusion_flux_to(&parcels[&right],link.geometry,link.pairs,solve,&mut cancelled)
            .map_err(|error|if error==GasEnergyError::Cancelled {Error::Cancelled} else {Error::Face {id,error}})?;
        for species in &flux.species {
            check(&mut cancelled)?;
            let mass=species.mass_rate_left_to_right_kg_s*duration_s;
            if !mass.is_finite() || (species.mass_rate_left_to_right_kg_s!=0.0 && mass==0.0) {
                return Err(Error::Invalid("integrated diffusion species mass is unrepresentable"));
            }
            if mass==0.0 {continue;}
            if transfers.len()==MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid("diffusion graph exceeds the selective-transfer work bound"));
            }
            let (donor,receiver)=if mass>0.0 {(left,right)} else {(right,left)};
            transfers.push(SpeciesNetworkTransfer {donor,receiver,species:species.species.clone(),mass_kg:mass.abs()});
        }
        reports.push(DiffusionNetworkLinkStep {id,parcels:link.parcels,flux});
    }
    Ok((transfers, reports))
}

// Cancel conserved species reference datums before rounding. Ideal mixing,
// source caloric entropy and finite-volume mass redistribution all remain.
pub(crate) fn admit_entropy(
    parcels: &BTreeMap<u64, ClosedGasParcel>, next_parcels: &BTreeMap<u64, ClosedGasParcel>,
    tolerance_j_k: f64, mass_tolerance_kg: f64, mut cancelled: impl FnMut() -> bool,
) -> Result<f64, Error> {
    use super::super::entropy;
    if !parcels.keys().eq(next_parcels.keys()) {
        return Err(Error::Invalid("closed entropy audit requires the original parcel identities"));
    }
    let change = entropy::closed_entropy_change(
        parcels.iter().map(|(&id, old)| (id, old, &next_parcels[&id])),
        mass_tolerance_kg, &mut cancelled,
    ).map_err(|error| match error {
        entropy::Error::Invalid(reason) => Error::Invalid(reason),
        entropy::Error::Parcel { id, error } => Error::Parcel { id, error },
        entropy::Error::Cancelled => Error::Cancelled,
    })?;
    if change < -tolerance_j_k {
        return Err(Error::Invalid("diffusion network entropy decreased; refine the physical interval"));
    }
    Ok(change)
}
