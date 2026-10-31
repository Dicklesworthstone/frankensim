//! Unsplit gas transport and competitive adsorption on one original snapshot.
//!
//! Every physical gas path and finite adsorbent samples the same accepted gas.
//! All gross gas withdrawals are admitted together; one original gas inversion
//! per affected volume follows. Desorbed/incoming gas cannot be forwarded in the
//! interval. Independent finite solids retain their own sources and temperatures.
//!
//! Joint entropy includes ALL gas volumes, bound species, dry skeletons and one
//! vacancy population per solid. A gas-only entropy gate is invalid for this
//! open subsystem and is deliberately not called. No absorption heat is added
//! twice. The physical path owners still admit their final driving-force order.
//!
//! This is first-order frozen-forcing coupling, not an implicit bath solve or
//! unconditional time-error guarantee. Refine the interval for feedback. Fixed
//! gas volumes and rigid skeletons exclude swelling, pore-volume work, external
//! heating, and simultaneous surface hopping. Callers must not also apply a
//! standalone adsorption update for these same surfaces in this interval.

use super::*;
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::network::{self as gas_network, SpeciesNetworkError, SpeciesNetworkStep,
        boundary::SpeciesNetworkBoundary,
        transport::{self, GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve}},
}};
use std::collections::BTreeSet;
mod entropy;
/// Source-integrated heat exchange between the same finite dry skeletons.
pub mod conduction;
/// Source-bound sensible heat between each finite adsorbent and its bound gas.
pub mod gas_heat;

/// A distinct finite site population exposed to one original gas volume.
/// Multiple populations may share a volume, but may not alias the same solid.
#[derive(Clone, Copy, Debug)]
pub struct AdsorptionNetworkBinding<'a> {
    pub id: u64,
    pub parcel: u64,
    pub solid: &'a CompetitiveSolid,
    /// The gas caloric controls must equal this volume's endpoint controls.
    /// Local solid/kinetic limits remain active; the strictest entropy allowance
    /// is applied to the complete closed system, not to an open pair.
    pub solve: CompetitiveSolve,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdsorptionNetworkError {
    Invalid(&'static str),
    Transport(GasTransportError),
    Gas(SpeciesNetworkError),
    Adsorption { id: u64, error: AdsorptionError },
    Thermodynamic(AdsorptionError),
    /// Original source/geometry/flux refusal, attributed to the physical link.
    Conduction { id: u64, reason: String },
    Cancelled,
}
impl core::fmt::Display for AdsorptionNetworkError {
    fn fmt(&self, f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"joint adsorbing gas network: {self:?}")
    }
}
impl std::error::Error for AdsorptionNetworkError {}
impl From<AdsorptionError> for AdsorptionNetworkError {
    fn from(error:AdsorptionError)->Self {
        if error==AdsorptionError::Cancelled { Self::Cancelled } else { Self::Thermodynamic(error) }
    }
}
impl From<SpeciesNetworkError> for AdsorptionNetworkError {
    fn from(error:SpeciesNetworkError)->Self {
        if error==SpeciesNetworkError::Cancelled { Self::Cancelled } else { Self::Gas(error) }
    }
}
impl From<GasTransportError> for AdsorptionNetworkError {
    fn from(error:GasTransportError)->Self {
        if error==GasTransportError::Cancelled { Self::Cancelled } else { Self::Transport(error) }
    }
}
type NetworkResult<T> = core::result::Result<T, AdsorptionNetworkError>;
use AdsorptionNetworkError as Error;
fn site_error(id:u64,error:AdsorptionError)->Error {
    if error==AdsorptionError::Cancelled { Error::Cancelled } else { Error::Adsorption { id,error } }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionNetworkSpeciesStep {
    pub species: SpeciesId,
    pub source: [LangmuirProperties;2],
    pub equilibrium_coverage: [f64;2],
}
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionNetworkSolidStep {
    pub id: u64,
    pub parcel: u64,
    /// Canonical species order matches the original solid and species below.
    pub sites: CompetitiveSiteStep,
    pub species: Vec<AdsorptionNetworkSpeciesStep>,
    pub gas_temperature_k: [f64;2],
    pub solid_temperature_k: [f64;2],
    pub dry_sensible_energy_change_j: f64,
    pub bound_sensible_energy_change_j: f64,
    pub solid_energy_change_j: f64,
    /// Signed physical heat received from the other retained finite solids.
    pub conductive_heat_j: f64,
    /// Signed non-advective heat received from this solid's finite gas volume.
    pub gas_sensible_heat_j: f64,
    pub energy_residual_j: f64,
    pub thermal_iterations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionNetworkStep {
    pub duration_s: f64,
    pub transport: GasTransportRates,
    pub exchange: SpeciesNetworkStep,
    pub solids: Vec<AdsorptionNetworkSolidStep>,
    /// Original layer/interface sources and the actually applied internal heat.
    pub conduction: Vec<conduction::AdsorbentConductionSample>,
    /// Explicit sensible heat and pinned sources, separate from species streams.
    pub gas_heat: Vec<gas_heat::AdsorbentGasHeatSample>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Reference-cancelled sum of original gas and finite solid caloric defects.
    pub total_energy_residual_j: f64,
    /// Entropy change of the complete CLOSED gas+solid system [J/K].
    pub entropy_change_j_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct AdsorptionNetworkCandidate {
    parcels: BTreeMap<u64,ClosedGasParcel>,
    solids: BTreeMap<u64,CompetitiveSolid>,
    step: AdsorptionNetworkStep,
}
impl AdsorptionNetworkCandidate {
    pub fn parcels(&self)->&BTreeMap<u64,ClosedGasParcel> { &self.parcels }
    pub fn solids(&self)->&BTreeMap<u64,CompetitiveSolid> { &self.solids }
    pub const fn step(&self)->&AdsorptionNetworkStep { &self.step }
    pub fn into_parts(self)->(BTreeMap<u64,ClosedGasParcel>,BTreeMap<u64,CompetitiveSolid>,AdsorptionNetworkStep) {
        (self.parcels,self.solids,self.step)
    }
}

/// Advance every named physical mechanism and adsorbent exactly once.
/// The common retained catalog must cover all active gas and adsorbed species
/// with their original thermochemistry. It is not a replacement gas composition.
/// No gas or solid is published if any late source, balance or cancellation gate
/// fails. Reports separate mechanism sampling from joint thermodynamic admission.
pub fn step_candidate(
    parcels:&BTreeMap<u64,ClosedGasParcel>, paths:GasTransportPaths<'_>,
    bindings:&[AdsorptionNetworkBinding<'_>], catalog:&ThermochemicalGasModel,
    duration_s:f64, endpoints:&BTreeMap<u64,GasEnergySolve>, controls:GasTransportSolve,
    cancelled:impl FnMut()->bool,
)->NetworkResult<AdsorptionNetworkCandidate> {
    step_candidate_with_conduction(parcels,paths,bindings,catalog,duration_s,endpoints,
        controls,&[],cancelled)
}

/// Advance original gas paths, finite adsorption heat and internal skeleton
/// conduction from one accepted snapshot. Source layers and optional ordered
/// interfaces are sampled at the actual solid temperatures. Their equal and
/// opposite heat increments enter the SAME original finite caloric inversion.
/// Endpoint source and complete-system entropy admission remain transactional.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_conduction(
    parcels:&BTreeMap<u64,ClosedGasParcel>, paths:GasTransportPaths<'_>,
    bindings:&[AdsorptionNetworkBinding<'_>], catalog:&ThermochemicalGasModel,
    duration_s:f64, endpoints:&BTreeMap<u64,GasEnergySolve>, controls:GasTransportSolve,
    links:&[conduction::AdsorbentConductionLink<'_, '_>],
    cancelled:impl FnMut()->bool,
)->NetworkResult<AdsorptionNetworkCandidate> {
    step_candidate_with_heat(parcels,paths,bindings,catalog,duration_s,endpoints,
        controls,links,&[],cancelled)
}

/// Advance skeleton conduction and gas/solid sensible heat together with the
/// original transport and adsorption snapshot. Each gas heat path debits the
/// SAME gas inversion and credits the SAME dry+bound inversion once. Source
/// coefficients see the actual gas/solid temperatures and pressure; endpoint
/// support, original absolute energy budgets and complete entropy stay active.
/// Heat is frozen over this interval; its accepted value is not a second
/// adsorption enthalpy or an unconditional stability/time-accuracy guarantee.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate_with_heat(
    parcels:&BTreeMap<u64,ClosedGasParcel>, paths:GasTransportPaths<'_>,
    bindings:&[AdsorptionNetworkBinding<'_>], catalog:&ThermochemicalGasModel,
    duration_s:f64, endpoints:&BTreeMap<u64,GasEnergySolve>, controls:GasTransportSolve,
    links:&[conduction::AdsorbentConductionLink<'_, '_>],
    gas_links:&[gas_heat::AdsorbentGasHeatLink<'_>],
    mut cancelled:impl FnMut()->bool,
)->NetworkResult<AdsorptionNetworkCandidate> {
    poll(&mut cancelled)?;
    let work=bindings.iter().try_fold(bindings.len(),|n,b|n.checked_add(b.solid.moles.len()));
    if parcels.is_empty() || parcels.len()>MAX_SNAPSHOT_PARCELS || bindings.len()>MAX_SNAPSHOT_PARCELS
        || work.is_none_or(|n|n>MAX_SNAPSHOT_TRANSFERS)
        || catalog.thermodynamics().components().len()>MAX_COMPETITIVE_SPECIES
        || !parcels.keys().eq(endpoints.keys()) || !positive(duration_s)
        || !positive(controls.entropy_tolerance_j_k)
        || !positive(controls.exchange.mass_tolerance_kg) || !positive(controls.exchange.energy_tolerance_j)
        || !positive(controls.exchange.maximum_species_fraction) || controls.exchange.maximum_species_fraction>=1.0 {
        return Err(Error::Invalid("bounded shared gas/site graph, complete catalog/endpoints and positive physical budgets required"));
    }
    entropy::admit_catalog(parcels,catalog,&mut cancelled)?;
    let mut ordered=BTreeMap::new(); let mut owners=BTreeSet::new();
    let mut entropy_limit=controls.entropy_tolerance_j_k;
    for binding in bindings {
        poll(&mut cancelled)?;
        if !parcels.contains_key(&binding.parcel) || ordered.insert(binding.id,binding).is_some()
            || !owners.insert(binding.solid as *const CompetitiveSolid)
            || endpoints[&binding.parcel]!=binding.solve.gas.gas {
            return Err(Error::Invalid("adsorption needs unique nonaliased solids, existing volumes and identical endpoint caloric controls"));
        }
        binding.solve.validate(duration_s).map_err(|e|site_error(binding.id,e))?;
        entropy_limit=entropy_limit.min(binding.solve.entropy_tolerance_j_k);
    }
    let heat=conduction::sample(&ordered,links,duration_s,&mut cancelled)?;
    let mut gas_heat=gas_heat::sample(parcels,&ordered,gas_links,duration_s,&mut cancelled)?;
    let sampled=transport::sample_transfers(parcels,paths,duration_s,endpoints,controls,&mut cancelled)?;
    let mut prepared=Vec::with_capacity(bindings.len());
    let mut heat_rounding=Vec::with_capacity(bindings.len());
    for (&id,&binding) in &ordered {
        let skeleton=heat.into_solid_j.get(&id).copied().unwrap_or(0.0);
        let sensible=gas_heat.into_solid_j.get(&id).copied().unwrap_or(0.0);
        let total=sum([skeleton,sensible])?;
        heat_rounding.push(sum([total,-skeleton,-sensible])?);
        let proposal=prepare_adsorption(binding.solid,&parcels[&binding.parcel],catalog,
            duration_s,binding.solve,total,
            &mut cancelled).map_err(|e|site_error(id,e))?;
        prepared.push((binding,proposal));
    }
    let boundaries:Vec<_>=prepared.iter().map(|(binding,proposal)|SpeciesNetworkBoundary {
        id:binding.id,parcel:binding.parcel,catalog,species:&proposal.boundary }).collect();
    let gas=gas_network::exchange_candidate_with_external_heat(parcels,&sampled.transfers,&sampled.heat,
        &boundaries,&gas_heat.into_gas,endpoints,controls.exchange,&mut cancelled)?;
    sampled.rates.admit_endpoint_driving_forces(gas.parcels(),&mut cancelled)?;
    let mut solid_reports=Vec::with_capacity(prepared.len());
    // Include any representational residual in the antisymmetric link lowering
    // in the one closed-system gate; it is never an external heat supply.
    let mut defects=vec![gas.step().total_energy_residual_j,heat.balance_residual_j,
        gas_heat.balance_residual_j,sum(heat_rounding)?];
    for (binding,proposal) in &prepared {
        poll(&mut cancelled)?;
        let old=binding.solid; let next=&proposal.next;
        let original_gas=&parcels[&binding.parcel]; let next_gas=&gas.parcels()[&binding.parcel];
        let tg0=original_gas.state().state().temperature; let tg1=next_gas.state().state().temperature;
        if (tg1-tg0).abs()>binding.solve.maximum_temperature_change_k
            || (tg1-next.temperature_k).abs()>binding.solve.maximum_gas_surface_temperature_difference_k {
            return Err(site_error(binding.id,invalid("joint gas endpoint exceeds the adsorption thermal regime")));
        }
        let gas_account=gas.step().parcels.iter().find(|r|r.id==binding.parcel)
            .ok_or(Error::Invalid("joint exchange lost a bound gas account"))?;
        if gas_account.energy_residual_j.abs()>0.5*binding.solve.gas.energy_tolerance_j
            || gas_account.mass_residual_kg.abs()>binding.solve.gas.mass_tolerance_kg
            || gas_account.maximum_species_mass_residual_kg>binding.solve.gas.mass_tolerance_kg {
            return Err(site_error(binding.id,invalid("joint gas endpoint exceeds a bound surface's local balance budget")));
        }
        let mut after=Vec::with_capacity(old.moles.len());
        for law in &old.model.channels {
            poll(&mut cancelled)?;
            after.push(law.properties(next.temperature_k,next_gas.state().state().pressure)
                .map_err(|e|site_error(binding.id,e))?);
        }
        let target0=equilibrium_coverages(&proposal.before,&old.model,original_gas)?;
        let target1=equilibrium_coverages(&after,&old.model,next_gas)?;
        let species=proposal.before.iter().cloned().zip(after).enumerate().map(|(i,(a,b))|
            AdsorptionNetworkSpeciesStep { species:old.model.channels[i].species.clone(),
                source:[a,b],equilibrium_coverage:[target0[i],target1[i]] }).collect();
        let defect=finite(sum(proposal.sensible)?-proposal.heat)?;
        defects.push(defect);
        solid_reports.push(AdsorptionNetworkSolidStep { id:binding.id,parcel:binding.parcel,
            sites:proposal.sites.clone(),species,gas_temperature_k:[tg0,tg1],
            solid_temperature_k:[old.temperature_k,next.temperature_k],
            dry_sensible_energy_change_j:proposal.sensible[0],bound_sensible_energy_change_j:proposal.sensible[1],
            solid_energy_change_j:finite(sum(proposal.sensible)?+sum(proposal.binding.iter().copied())?)?,
            conductive_heat_j:heat.into_solid_j.get(&binding.id).copied().unwrap_or(0.0),
            gas_sensible_heat_j:gas_heat.into_solid_j.get(&binding.id).copied().unwrap_or(0.0),
            energy_residual_j:defect,thermal_iterations:proposal.iterations });
    }
    conduction::admit_endpoint(&prepared,links,&mut cancelled)?;
    gas_heat::admit_endpoint(gas.parcels(),&prepared,gas_links,&mut gas_heat.samples,&mut cancelled)?;
    let energy_residual=sum(defects)?;
    if energy_residual.abs()>controls.exchange.energy_tolerance_j {
        return Err(Error::Invalid("complete gas/solid energy residual exceeds the single interval budget"));
    }
    let (mass_residual,maximum_species)=entropy::mass_balance(parcels,gas.parcels(),&prepared,catalog,&mut cancelled)?;
    if mass_residual.abs()>controls.exchange.mass_tolerance_kg || maximum_species>controls.exchange.mass_tolerance_kg {
        return Err(Error::Invalid("complete gas/solid species inventory balance failed"));
    }
    let entropy=entropy::change(parcels,gas.parcels(),&prepared,catalog,gas.step(),&mut cancelled)?;
    if entropy < -entropy_limit { return Err(AdsorptionError::EntropyDecreased { change_j_k:entropy }.into()); }
    poll(&mut cancelled)?;
    let (parcels,exchange)=gas.into_parts();
    let solids=prepared.into_iter().map(|(binding,proposal)|(binding.id,proposal.next)).collect();
    Ok(AdsorptionNetworkCandidate { parcels,solids,step:AdsorptionNetworkStep { duration_s,
        transport:sampled.rates,exchange,solids:solid_reports,conduction:heat.samples,gas_heat:gas_heat.samples,
        total_mass_residual_kg:mass_residual,
        maximum_species_mass_residual_kg:maximum_species,total_energy_residual_j:energy_residual,
        entropy_change_j_k:entropy } })
}

#[cfg(test)]
mod tests;
