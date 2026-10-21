//! Simultaneous selective transport on a graph of original material inventories.
//! All edges use original donor enthalpy; combine every flux before rebuilding
//! an endpoint. No newly received species can leave again within this interval.
/// Pressure restrictions and molecular mixing share one inventory update.
pub mod transport;

use std::collections::{BTreeMap, BTreeSet};
use super::{ClosedGasParcel, GasEnergyError, GasEnergySolve, SpeciesId, Sum,
    endpoint, inventory, mixed_model_from_portions, model_for_inventory, positive,
    source_enthalpy, stored_change};
use super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};

/// Nonnegative integrated mass of ONE species, not the whole donor mixture.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkTransfer {
    pub donor: u64,
    pub receiver: u64,
    pub species: SpeciesId,
    pub mass_kg: f64,
}
/// Absolute local and global limits, independent of formation-energy magnitude.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesNetworkSolve {
    /// Gross withdrawal divided by the ORIGINAL individual species mass.
    pub maximum_species_fraction: f64,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum SpeciesNetworkError {
    Invalid(&'static str),
    Parcel { id: u64, error: GasEnergyError },
    Withdrawal { id: u64, species: SpeciesId, fraction: f64, maximum: f64 },
    Cancelled,
}
impl core::fmt::Display for SpeciesNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "selective gas network: {self:?}")
    }
}
impl std::error::Error for SpeciesNetworkError {}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkFlux {
    pub transfer: SpeciesNetworkTransfer,
    pub source_enthalpy_j_kg: f64,
    pub carried_enthalpy_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkParcelStep {
    pub id: u64,
    pub incoming_mass_kg: f64,
    pub outgoing_mass_kg: f64,
    pub net_enthalpy_j: f64,
    pub energy_change_j: f64,
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub iterations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkStep {
    pub fluxes: Vec<SpeciesNetworkFlux>,
    pub parcels: Vec<SpeciesNetworkParcelStep>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub total_energy_residual_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: SpeciesNetworkStep,
}
impl SpeciesNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    pub const fn step(&self) -> &SpeciesNetworkStep { &self.step }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, SpeciesNetworkStep) {
        (self.parcels, self.step)
    }
}
use SpeciesNetworkError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn material(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled { Error::Cancelled } else { Error::Parcel { id, error } }
}
#[derive(Default)]
struct Budget {
    incoming: BTreeSet<u64>,
    mass_in: Sum,
    mass_out: Sum,
    enthalpy: Sum,
    relative_energy: Sum,
    species: BTreeMap<SpeciesId, Sum>,
    withdrawn: BTreeMap<SpeciesId, Sum>,
}

/// Update all species-selective connections from ONE unchanged network snapshot.
///
/// Parallel records are combined in canonical (donor,receiver,species,mass)
/// order. Total and species-specific GROSS outflow is admitted before any
/// endpoint is evaluated; incoming mass cannot conceal depletion. A full source
/// catalog join for each receiving node uses the original material owner, then
/// only the actual target species inventories set its composition. No endpoint
/// is repeatedly mixed as edges are traversed. Idle nodes stay bitwise unchanged.
///
/// Fixed-volume upstream-state transport only. Prescribed fluxes alone do not
/// establish a diffusion/pressure law, pumping work or nonnegative entropy.
/// Call the physical transport owner for rates and refine its timestep.
pub fn exchange_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>, transfers: &[SpeciesNetworkTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>, controls: SpeciesNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    check(&mut cancelled)?;
    if parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS
        || transfers.len() > MAX_SNAPSHOT_TRANSFERS || !parcels.keys().eq(endpoints.keys())
        || !positive(controls.maximum_species_fraction) || controls.maximum_species_fraction >= 1.0
        || !positive(controls.mass_tolerance_kg) || !positive(controls.energy_tolerance_j) {
        return Err(Error::Invalid("bounded graph, exact endpoint map and positive physical limits required"));
    }
    for (&id, solve) in endpoints {
        check(&mut cancelled)?;
        solve.validate().map_err(|e| material(id, e))?;
    }
    let mut rows = Vec::with_capacity(transfers.len());
    for row in transfers {
        check(&mut cancelled)?;
        if row.donor == row.receiver || !parcels.contains_key(&row.donor)
            || !parcels.contains_key(&row.receiver) || !row.mass_kg.is_finite() || row.mass_kg < 0.0 {
            return Err(Error::Invalid("finite nonnegative species mass and distinct existing owners required"));
        }
        if !parcels[&row.donor].model().thermodynamics().components().iter()
            .any(|s| s.species() == &row.species) {
            return Err(Error::Invalid("selective edge names a species absent from its original donor"));
        }
        if row.mass_kg > 0.0 { rows.push(row); }
    }
    rows.sort_by(|a,b| (a.donor,a.receiver,&a.species).cmp(&(b.donor,b.receiver,&b.species))
        .then_with(|| a.mass_kg.total_cmp(&b.mass_kg)));
    let mut combined: BTreeMap<(u64,u64,SpeciesId), Sum> = BTreeMap::new();
    for row in rows {
        check(&mut cancelled)?;
        combined.entry((row.donor,row.receiver,row.species.clone())).or_default()
            .add(row.mass_kg).map_err(|e| material(row.donor,e))?;
    }
    let mut budgets: BTreeMap<u64, Budget> = BTreeMap::new();
    let mut fluxes = Vec::with_capacity(combined.len());
    for ((donor,receiver,species), amount) in combined {
        check(&mut cancelled)?;
        let dm = amount.value;
        let h = source_enthalpy(&parcels[&donor], &species).map_err(|e| material(donor,e))?;
        let carried = dm*h;
        for (id, sign) in [(donor,-1.0),(receiver,1.0)] {
            let map = |e| material(id,e);
            let budget = budgets.entry(id).or_default();
            if id == donor {
                budget.mass_out.add(dm).map_err(map)?;
                budget.withdrawn.entry(species.clone()).or_default().add(dm).map_err(map)?;
            } else {
                budget.incoming.insert(donor);
                budget.mass_in.add(dm).map_err(map)?;
            }
            budget.enthalpy.add(sign*carried).map_err(map)?;
            budget.relative_energy.add(sign*dm*(h-parcels[&id].state().specific_internal_energy_j_kg())).map_err(map)?;
            budget.species.entry(species.clone()).or_default().add(sign*dm).map_err(map)?;
        }
        fluxes.push(SpeciesNetworkFlux { transfer: SpeciesNetworkTransfer { donor,receiver,species,mass_kg:dm },
            source_enthalpy_j_kg:h, carried_enthalpy_j:carried });
    }
    // Check aggregate donor allowances before the first nonlinear material solve.
    let mut originals = BTreeMap::new();
    for (&id,budget) in &budgets {
        check(&mut cancelled)?;
        let before = inventory(parcels[&id].model(),parcels[&id].mass_kg()).map_err(|e| material(id,e))?;
        for (species,withdrawn) in &budget.withdrawn {
            let fraction = withdrawn.value/before.get(species).copied().unwrap_or(0.0);
            if !fraction.is_finite() || fraction > controls.maximum_species_fraction {
                return Err(Error::Withdrawal { id,species:species.clone(),fraction,maximum:controls.maximum_species_fraction });
            }
        }
        originals.insert(id,before);
    }
    let mut next = BTreeMap::new(); let mut accounts = Vec::with_capacity(parcels.len());
    let (mut total_mass,mut total_energy) = (Sum::default(),Sum::default());
    let mut global_species: BTreeMap<SpeciesId,Sum> = BTreeMap::new();
    for (&id,old) in parcels {
        check(&mut cancelled)?;
        let Some(budget) = budgets.get(&id) else {
            next.insert(id,old.clone());
            accounts.push(SpeciesNetworkParcelStep { id,incoming_mass_kg:0.0,outgoing_mass_kg:0.0,
                net_enthalpy_j:0.0,energy_change_j:0.0,energy_residual_j:0.0,mass_residual_kg:0.0,
                maximum_species_mass_residual_kg:0.0,iterations:0 });
            continue;
        };
        let map = |e| material(id,e);
        let before = &originals[&id];
        let mut target = before.clone();
        for (species,change) in &budget.species {
            let value = before.get(species).copied().unwrap_or(0.0);
            let updated = value+change.value;
            if !positive(updated) || (change.value != 0.0 && updated == value) {
                return Err(Error::Invalid("selective network exhausted or lost a species increment"));
            }
            target.insert(species.clone(),updated);
        }
        // Join source catalogs, not physical temperatures or transport values.
        let mut portions = vec![(old.model(),old.mass_kg())];
        let mut catalog_mass = Sum::default(); catalog_mass.add(old.mass_kg()).map_err(map)?;
        for donor in &budget.incoming {
            check(&mut cancelled)?;
            let original = &parcels[donor];
            portions.push((original.model(),original.mass_kg()));
            catalog_mass.add(original.mass_kg()).map_err(map)?;
        }
        let catalog = mixed_model_from_portions(&portions,catalog_mass.value).map_err(map)?;
        let delta_mass = budget.mass_in.value-budget.mass_out.value;
        let mass = old.mass_kg()+delta_mass;
        if !positive(mass) { return Err(Error::Invalid("selective network final mass is unrepresentable")); }
        let (new,iterations) = if target == *before && delta_mass == 0.0 && budget.relative_energy.value == 0.0 {
            (old.clone(),0)
        } else {
            endpoint(old,model_for_inventory(&catalog,&target).map_err(map)?,mass,
                budget.relative_energy.value/mass,endpoints[&id],&mut cancelled).map_err(map)?
        };
        let du = stored_change(old,&new);
        let dm = new.mass_kg()-old.mass_kg();
        let er = du-budget.enthalpy.value; let mr = dm-delta_mass;
        if ![du,dm,er,mr].iter().all(|v| v.is_finite())
            || er.abs() > controls.energy_tolerance_j || mr.abs() > controls.mass_tolerance_kg {
            return Err(Error::Invalid("selective network local mass or energy balance failed"));
        }
        let actual = inventory(new.model(),new.mass_kg()).map_err(map)?;
        let mut maximum: f64 = 0.0;
        for (species,expected) in &target {
            let value = actual.get(species).copied().unwrap_or(0.0);
            maximum = maximum.max((value-expected).abs());
            global_species.entry(species.clone()).or_default()
                .add(value-before.get(species).copied().unwrap_or(0.0)).map_err(map)?;
        }
        if maximum > controls.mass_tolerance_kg {
            return Err(Error::Invalid("selective network local species balance failed"));
        }
        total_mass.add(dm).map_err(map)?; total_energy.add(du).map_err(map)?;
        accounts.push(SpeciesNetworkParcelStep { id,incoming_mass_kg:budget.mass_in.value,
            outgoing_mass_kg:budget.mass_out.value,net_enthalpy_j:budget.enthalpy.value,
            energy_change_j:du,energy_residual_j:er,mass_residual_kg:mr,
            maximum_species_mass_residual_kg:maximum,iterations });
        next.insert(id,new);
    }
    let maximum = global_species.values().map(|s|s.value.abs()).fold(0.0_f64,f64::max);
    if maximum > controls.mass_tolerance_kg || total_mass.value.abs() > controls.mass_tolerance_kg
        || total_energy.value.abs() > controls.energy_tolerance_j {
        return Err(Error::Invalid("selective network global species, mass or energy balance failed"));
    }
    check(&mut cancelled)?;
    Ok(SpeciesNetworkCandidate { parcels:next,step:SpeciesNetworkStep { fluxes,parcels:accounts,
        total_mass_residual_kg:total_mass.value,maximum_species_mass_residual_kg:maximum,
        total_energy_residual_j:total_energy.value } })
}
