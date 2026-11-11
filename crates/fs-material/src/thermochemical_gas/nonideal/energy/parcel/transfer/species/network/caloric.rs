//! Original partial-enthalpy streams and one finite PR caloric solve per node.
use super::super::super::mixing;
use super::super::{endpoint, endpoint_density};
use super::snapshot::Snapshot;
use super::{
    Error, GasEnergySolve, Parcels, RealGasSpeciesNetworkCandidate, RealGasSpeciesNetworkStep,
    SpeciesNetworkFlux, SpeciesNetworkParcelStep, SpeciesNetworkSolve, Sum, check, material,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::model_for_inventory;
use crate::thermochemical_gas::energy::parcel::transfer::{inventory, mixed_model_from_portions};
use crate::thermochemical_gas::{failure, positive};
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

pub(super) fn streams(
    parcels: &Parcels<'_>,
    snapshot: &mut Snapshot<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<SpeciesNetworkFlux>, Error> {
    // Exactly one bounded partial-caloric evaluation per ORIGINAL active donor.
    let mut partial = BTreeMap::new();
    for row in &snapshot.flows {
        check(cancelled)?;
        if let std::collections::btree_map::Entry::Vacant(entry) = partial.entry(row.donor) {
            let gas = parcels[&row.donor];
            entry.insert(
                gas.model()
                    .partial_calorics_at_density(gas.density_kg_m3(), gas.state().temperature_k())
                    .map_err(|e| material(row.donor, e))?,
            );
        }
    }
    let mut fluxes = Vec::with_capacity(snapshot.flows.len());
    for row in &snapshot.flows {
        check(cancelled)?;
        let upstream = parcels[&row.donor];
        let (index, source) = upstream
            .model()
            .source()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
            .find(|(_, species)| species.species() == &row.species)
            .ok_or(Error::Invalid(
                "original real-gas donor lost a source species",
            ))?;
        let h = partial[&row.donor].partial_specific_enthalpies_j_kg()[index];
        let carried_enthalpy_j = row.mass_kg * h;
        if !carried_enthalpy_j.is_finite() {
            return Err(Error::Invalid(
                "real-gas gross stream enthalpy is unrepresentable",
            ));
        }
        let departure =
            partial[&row.donor].departure().components()[index].residual_enthalpy_j_mol();
        for (id, sign) in [(row.donor, -1.0), (row.receiver, 1.0)] {
            let delta_h = source
                .molar_enthalpy_change(
                    Temperature::new(parcels[&id].state().temperature_k()),
                    Temperature::new(upstream.state().temperature_k()),
                )
                .map_err(|e| {
                    material(id, failure("real-gas network reference enthalpy", e).into())
                })?
                .value();
            let amount = sign * row.mass_kg / source.molar_mass().value();
            if !amount.is_finite() || amount == 0.0 {
                return Err(Error::Invalid(
                    "real-gas network molar stream is unrepresentable",
                ));
            }
            let budget = snapshot
                .budgets
                .get_mut(&id)
                .expect("admitted stream endpoint");
            budget
                .relative_stream
                .add(amount * delta_h)
                .map_err(|e| material(id, e))?;
            budget
                .relative_stream
                .add(amount * departure)
                .map_err(|e| material(id, e))?;
        }
        fluxes.push(SpeciesNetworkFlux {
            transfer: row.clone(),
            source_enthalpy_j_kg: h,
            carried_enthalpy_j,
        });
    }
    Ok(fluxes)
}

pub(super) fn finish(
    parcels: &Parcels<'_>,
    snapshot: Snapshot<'_>,
    fluxes: Vec<SpeciesNetworkFlux>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    let boundaries = snapshot.boundary.fluxes(parcels, cancelled)?;
    let mut next = BTreeMap::new();
    let mut accounts = Vec::with_capacity(parcels.len());
    let (mut total_mass, mut total_energy) = (Sum::default(), Sum::default());
    let mut global_species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for (&id, &old) in parcels {
        check(cancelled)?;
        let Some(budget) = snapshot.budgets.get(&id) else {
            next.insert(id, old.clone());
            accounts.push(SpeciesNetworkParcelStep {
                id,
                incoming_mass_kg: 0.0,
                outgoing_mass_kg: 0.0,
                net_enthalpy_j: 0.0,
                net_additional_heat_j: 0.0,
                energy_change_j: 0.0,
                energy_residual_j: 0.0,
                mass_residual_kg: 0.0,
                maximum_species_mass_residual_kg: 0.0,
                iterations: 0,
            });
            continue;
        };
        let map = |e| material(id, e);
        let before = &snapshot.before[&id];
        let mut target = before.clone();
        let mut delta_mass = Sum::default();
        for (species, change) in &budget.species {
            let previous = before.get(species).copied().unwrap_or(0.0);
            let updated = previous + change.value;
            if !positive(updated)
                || (change.value != 0.0 && updated == previous)
                || ((updated - previous) - change.value).abs() > controls.mass_tolerance_kg
            {
                return Err(Error::Invalid(
                    "real-gas network species inventory is unrepresentable or exhausted",
                ));
            }
            target.insert(species.clone(), updated);
            delta_mass.add(change.value).map_err(map)?;
        }
        // A net within the parcel's allowance is the rounding of a mass its
        // streams' own law balances (a mass-frame face; the lowering derives
        // the bound), so the parcel keeps its mass exactly. With no allowance
        // this is the exact test it always was. delta_mass itself stays the
        // computed value for the reference energy and the mass residual below.
        let allowance = snapshot.mass_allowance_kg.get(&id).copied().unwrap_or(0.0);
        let net = if delta_mass.value.abs() <= allowance {
            0.0
        } else {
            delta_mass.value
        };
        let mass = old.mass_kg() + net;
        let density = endpoint_density(old, mass);
        if !positive(mass)
            || !positive(density)
            || (net != 0.0 && mass == old.mass_kg())
            || (mass != old.mass_kg() && density == old.density_kg_m3())
            || ((mass - old.mass_kg()) - net).abs() > controls.mass_tolerance_kg
        {
            return Err(Error::Invalid(
                "real-gas network mass/density increment is unrepresentable",
            ));
        }
        // This is a catalog join, never an upstream state update. Unrelated
        // nodes and outgoing receivers cannot introduce extra source demands.
        let mut originals = vec![old.model()];
        let external_catalogs = snapshot.boundary.catalogs.get(&id);
        let catalog = if budget.incoming.is_empty() && external_catalogs.is_none() {
            old.model().source().clone()
        } else {
            let mut portions = vec![(old.model().source(), old.mass_kg())];
            let mut catalog_mass = Sum::default();
            catalog_mass.add(old.mass_kg()).map_err(map)?;
            for donor in &budget.incoming {
                check(cancelled)?;
                let upstream = parcels[donor];
                originals.push(upstream.model());
                portions.push((upstream.model().source(), upstream.mass_kg()));
                catalog_mass.add(upstream.mass_kg()).map_err(map)?;
            }
            if let Some(external_catalogs) = external_catalogs {
                for catalog in external_catalogs {
                    check(cancelled)?;
                    originals.push(catalog);
                    // Positive catalog weights gather source ownership only;
                    // the authoritative target inventory sets composition.
                    portions.push((catalog.source(), old.mass_kg()));
                    catalog_mass.add(old.mass_kg()).map_err(map)?;
                }
            }
            mixed_model_from_portions(&portions, catalog_mass.value).map_err(map)?
        };
        let same_inventory = target == *before;
        let source = if same_inventory {
            old.model().source().clone()
        } else {
            model_for_inventory(&catalog, &target).map_err(map)?
        };
        let model = if same_inventory && budget.incoming.is_empty() && external_catalogs.is_none() {
            old.model().clone()
        } else {
            mixing::model_for_source(&originals, source).map_err(map)?
        };
        let temperature = old.state().temperature_k();
        let anchor = model
            .caloric_at_density(density, temperature)
            .map_err(map)?
            .residual_energy_j_kg;
        let departure = old.state().root().residual_internal_energy_j_mol()
            / old.model().mechanics().molar_mass().value();
        let mut net_moles = Sum::default();
        let mut net_enthalpy = Sum::default();
        // Group signed species mass BEFORE multiplying its formation reference.
        // Gross incoming/outgoing enthalpy reports are never subtracted here.
        for (species, change) in &budget.species {
            if change.value == 0.0 {
                continue;
            }
            let original = catalog
                .thermodynamics()
                .components()
                .iter()
                .find(|s| s.species() == species)
                .ok_or(Error::Invalid(
                    "real-gas network catalog lost a target species",
                ))?;
            let amount = change.value / original.molar_mass().value();
            if !amount.is_finite() || amount == 0.0 {
                return Err(Error::Invalid(
                    "real-gas network net species amount is unrepresentable",
                ));
            }
            let properties = original
                .evaluate_properties(Temperature::new(temperature))
                .map_err(|e| {
                    material(id, failure("real-gas network reference calorics", e).into())
                })?;
            net_moles.add(amount).map_err(map)?;
            net_enthalpy
                .add(amount * properties.h().value())
                .map_err(map)?;
        }
        net_enthalpy
            .add(budget.relative_stream.value)
            .map_err(map)?;
        let mut caloric = Sum::default();
        caloric
            .add(old.mass_kg() * (departure - anchor))
            .map_err(map)?;
        caloric.add(-delta_mass.value * anchor).map_err(map)?;
        caloric.add(budget.relative_stream.value).map_err(map)?;
        caloric
            .add(R * temperature * net_moles.value)
            .map_err(map)?;
        caloric.add(budget.heat.value).map_err(map)?;
        let (candidate, residual, iterations) =
            if mass == old.mass_kg() && model == *old.model() && caloric.value == 0.0 {
                (old.clone(), 0.0, 0)
            } else {
                endpoint(old, model, mass, caloric.value, endpoints[&id], cancelled).map_err(map)?
            };
        if residual.abs() > controls.energy_tolerance_j {
            return Err(Error::Invalid(
                "real-gas network local energy budget exceeded",
            ));
        }
        let after = inventory(candidate.model().source(), candidate.mass_kg()).map_err(map)?;
        let mass_change = candidate.mass_kg() - old.mass_kg();
        let mass_residual_kg = mass_change - delta_mass.value;
        let mut local_species = 0.0_f64;
        for original in catalog.thermodynamics().components() {
            let species = original.species();
            let actual = after.get(species).copied().unwrap_or(0.0)
                - before.get(species).copied().unwrap_or(0.0);
            let intended = budget.species.get(species).map_or(0.0, |sum| sum.value);
            let defect = actual - intended;
            if !defect.is_finite() {
                return Err(Error::Invalid(
                    "real-gas network species defect is unrepresentable",
                ));
            }
            local_species = local_species.max(defect.abs());
            global_species
                .entry(species.clone())
                .or_default()
                .add(actual)
                .map_err(map)?;
        }
        if mass_residual_kg.abs() > controls.mass_tolerance_kg
            || local_species > controls.mass_tolerance_kg
        {
            return Err(Error::Invalid(
                "real-gas network local mass/species budget exceeded",
            ));
        }
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(local_species);
        total_mass.add(mass_change).map_err(map)?;
        total_energy.add(residual).map_err(map)?;
        let mut stored = Sum::default();
        stored.add(net_enthalpy.value).map_err(map)?;
        stored.add(budget.heat.value).map_err(map)?;
        stored.add(residual).map_err(map)?;
        accounts.push(SpeciesNetworkParcelStep {
            id,
            incoming_mass_kg: budget.mass_in.value,
            outgoing_mass_kg: budget.mass_out.value,
            net_enthalpy_j: net_enthalpy.value,
            net_additional_heat_j: budget.heat.value,
            energy_change_j: stored.value,
            energy_residual_j: residual,
            mass_residual_kg,
            maximum_species_mass_residual_kg: local_species,
            iterations,
        });
        next.insert(id, candidate);
    }
    total_mass
        .add(-snapshot.boundary.mass.value)
        .map_err(|e| material(*parcels.keys().next().expect("admitted nonempty graph"), e))?;
    for (species, external) in &snapshot.boundary.species {
        global_species
            .entry(species.clone())
            .or_default()
            .add(-external.value)
            .map_err(|e| material(*parcels.keys().next().expect("admitted nonempty graph"), e))?;
    }
    for sum in global_species.into_values() {
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(sum.value.abs());
    }
    if total_mass.value.abs() > controls.mass_tolerance_kg
        || maximum_species_mass_residual_kg > controls.mass_tolerance_kg
        || total_energy.value.abs() > controls.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "real-gas network did not close global mass/species/energy",
        ));
    }
    check(cancelled)?;
    Ok(RealGasSpeciesNetworkCandidate {
        parcels: next,
        step: RealGasSpeciesNetworkStep {
            fluxes,
            heat: snapshot.heat,
            external_heat: snapshot.external_heat,
            boundaries,
            parcels: accounts,
            total_mass_residual_kg: total_mass.value,
            maximum_species_mass_residual_kg,
            total_energy_residual_j: total_energy.value,
        },
    })
}
