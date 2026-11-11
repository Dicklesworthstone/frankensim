//! Canonical original-snapshot streams, heat and gross inventory admission.
use super::boundary::{self, BoundaryAccounts, RealGasSpeciesNetworkBoundary};
use super::{
    Error, GasEnergySolve, Parcels, SpeciesNetworkExternalHeat, SpeciesNetworkHeatTransfer,
    SpeciesNetworkSolve, SpeciesNetworkTransfer, Sum, check, material,
};
use crate::thermochemical_gas::energy::parcel::transfer::{
    inventory,
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
};
use crate::thermochemical_gas::positive;
use fs_thermochem::SpeciesId;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Default)]
pub(super) struct Budget {
    pub incoming: BTreeSet<u64>,
    pub mass_in: Sum,
    pub mass_out: Sum,
    pub species: BTreeMap<SpeciesId, Sum>,
    pub withdrawn: BTreeMap<SpeciesId, Sum>,
    pub heat: Sum,
    pub relative_stream: Sum,
}
pub(super) struct Snapshot<'a> {
    pub flows: Vec<SpeciesNetworkTransfer>,
    pub heat: Vec<SpeciesNetworkHeatTransfer>,
    pub external_heat: Vec<SpeciesNetworkExternalHeat>,
    pub budgets: BTreeMap<u64, Budget>,
    pub before: BTreeMap<u64, BTreeMap<SpeciesId, f64>>,
    pub boundary: BoundaryAccounts<'a>,
    /// Per-parcel bound on the rounding of a law-balanced net mass [kg].
    pub mass_allowance_kg: BTreeMap<u64, f64>,
}

// Public entry points call this before building even the borrowed input map.
pub(super) fn admit_work(
    parcels: usize,
    transfers: usize,
    heat: usize,
    external: usize,
    boundary: usize,
) -> Result<(), Error> {
    if parcels == 0
        || parcels > MAX_SNAPSHOT_PARCELS
        || transfers
            .checked_add(heat)
            .and_then(|n| n.checked_add(external))
            .and_then(|n| n.checked_add(boundary))
            .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
    {
        return Err(Error::Invalid(
            "real-gas network exceeds the original snapshot work bound",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare<'a>(
    parcels: &Parcels<'_>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    external: &[SpeciesNetworkExternalHeat],
    boundaries: &[RealGasSpeciesNetworkBoundary<'a>],
    mass_allowance_kg: &BTreeMap<u64, f64>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Snapshot<'a>, Error> {
    check(cancelled)?;
    admit_work(
        parcels.len(),
        transfers.len(),
        heat.len(),
        external.len(),
        boundary::work_count(boundaries)?,
    )?;
    if !parcels.keys().eq(endpoints.keys())
        || !positive(controls.maximum_species_fraction)
        || controls.maximum_species_fraction >= 1.0
        || !positive(controls.mass_tolerance_kg)
        || !positive(controls.energy_tolerance_j)
        || mass_allowance_kg
            .iter()
            .any(|(id, kg)| !parcels.contains_key(id) || !kg.is_finite() || *kg < 0.0)
    {
        return Err(Error::Invalid(
            "bounded real-gas graph, exact endpoint map and positive physical limits required",
        ));
    }
    for (&id, solve) in endpoints {
        check(cancelled)?;
        solve.validate().map_err(|e| material(id, e))?;
    }
    let mut ordered = Vec::with_capacity(transfers.len());
    for row in transfers {
        check(cancelled)?;
        if row.donor == row.receiver
            || !parcels.contains_key(&row.donor)
            || !parcels.contains_key(&row.receiver)
            || !row.mass_kg.is_finite()
            || row.mass_kg < 0.0
        {
            return Err(Error::Invalid(
                "real-gas stream requires nonnegative mass and distinct original owners",
            ));
        }
        if !parcels[&row.donor]
            .model()
            .source()
            .thermodynamics()
            .components()
            .iter()
            .any(|source| source.species() == &row.species)
        {
            return Err(Error::Invalid(
                "species absent from original real-gas donor",
            ));
        }
        if row.mass_kg > 0.0 {
            ordered.push(row);
        }
    }
    ordered.sort_by(|a, b| {
        (a.donor, a.receiver, &a.species)
            .cmp(&(b.donor, b.receiver, &b.species))
            .then_with(|| a.mass_kg.total_cmp(&b.mass_kg))
    });
    let mut combined: BTreeMap<(u64, u64, SpeciesId), Sum> = BTreeMap::new();
    for row in ordered {
        check(cancelled)?;
        combined
            .entry((row.donor, row.receiver, row.species.clone()))
            .or_default()
            .add(row.mass_kg)
            .map_err(|e| material(row.donor, e))?;
    }
    let mut flows = Vec::with_capacity(combined.len());
    let mut budgets: BTreeMap<u64, Budget> = BTreeMap::new();
    for ((donor, receiver, species), value) in combined {
        check(cancelled)?;
        for (id, sign) in [(donor, -1.0), (receiver, 1.0)] {
            let map = |e| material(id, e);
            let budget = budgets.entry(id).or_default();
            if id == donor {
                budget.mass_out.add(value.value).map_err(map)?;
                budget
                    .withdrawn
                    .entry(species.clone())
                    .or_default()
                    .add(value.value)
                    .map_err(map)?;
            } else {
                budget.mass_in.add(value.value).map_err(map)?;
                budget.incoming.insert(donor);
            }
            budget
                .species
                .entry(species.clone())
                .or_default()
                .add(sign * value.value)
                .map_err(map)?;
        }
        flows.push(SpeciesNetworkTransfer {
            donor,
            receiver,
            species,
            mass_kg: value.value,
        });
    }
    let mut ordered_heat = Vec::with_capacity(heat.len());
    for row in heat {
        check(cancelled)?;
        let [left, right] = row.parcels;
        if left == right
            || !parcels.contains_key(&left)
            || !parcels.contains_key(&right)
            || !row.heat_left_to_right_j.is_finite()
        {
            return Err(Error::Invalid(
                "real-gas heat requires finite energy and distinct original owners",
            ));
        }
        if row.heat_left_to_right_j != 0.0 {
            ordered_heat.push(if left < right {
                ([left, right], row.heat_left_to_right_j)
            } else {
                ([right, left], -row.heat_left_to_right_j)
            });
        }
    }
    ordered_heat.sort_by(|(a, aq), (b, bq)| a.cmp(b).then_with(|| aq.total_cmp(bq)));
    let mut combined_heat: BTreeMap<[u64; 2], Sum> = BTreeMap::new();
    for (pair, value) in ordered_heat {
        check(cancelled)?;
        combined_heat
            .entry(pair)
            .or_default()
            .add(value)
            .map_err(|e| material(pair[0], e))?;
    }
    let mut heat = Vec::with_capacity(combined_heat.len());
    for (parcels, value) in combined_heat {
        if value.value == 0.0 {
            continue;
        }
        for (id, sign) in [(parcels[0], -1.0), (parcels[1], 1.0)] {
            budgets
                .entry(id)
                .or_default()
                .heat
                .add(sign * value.value)
                .map_err(|e| material(id, e))?;
        }
        heat.push(SpeciesNetworkHeatTransfer {
            parcels,
            heat_left_to_right_j: value.value,
        });
    }
    let mut ordered_external = BTreeMap::new();
    for row in external {
        check(cancelled)?;
        if !parcels.contains_key(&row.parcel)
            || !row.heat_into_gas_j.is_finite()
            || ordered_external.insert(row.id, row).is_some()
        {
            return Err(Error::Invalid(
                "external real-gas heat needs unique identities, original owners and finite energy",
            ));
        }
    }
    let mut external_heat = Vec::with_capacity(ordered_external.len());
    for row in ordered_external.into_values() {
        check(cancelled)?;
        if row.heat_into_gas_j == 0.0 {
            continue;
        }
        budgets
            .entry(row.parcel)
            .or_default()
            .heat
            .add(row.heat_into_gas_j)
            .map_err(|e| material(row.parcel, e))?;
        external_heat.push(row.clone());
    }
    let boundary = boundary::accumulate(parcels, boundaries, &mut budgets, controls, cancelled)?;
    // Admit aggregate ORIGINAL donor inventories, including every external
    // withdrawal, before any partial-caloric
    // query, source mixing or endpoint solve. Inflow cannot fund an outflow.
    let mut before = BTreeMap::new();
    for (&id, budget) in &budgets {
        check(cancelled)?;
        let gas = parcels[&id];
        let original =
            inventory(gas.model().source(), gas.mass_kg()).map_err(|e| material(id, e))?;
        for (species, amount) in &budget.withdrawn {
            let fraction = amount.value / original.get(species).copied().unwrap_or(0.0);
            if !fraction.is_finite() || fraction > controls.maximum_species_fraction {
                return Err(Error::Withdrawal {
                    id,
                    species: species.clone(),
                    fraction,
                    maximum: controls.maximum_species_fraction,
                });
            }
        }
        before.insert(id, original);
    }
    Ok(Snapshot {
        flows,
        heat,
        external_heat,
        budgets,
        before,
        boundary,
        mass_allowance_kg: mass_allowance_kg.clone(),
    })
}
