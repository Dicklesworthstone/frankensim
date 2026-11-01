//! Open ideal-gas storage in the held-bath partial-enthalpy/entropy chart.
//!
//! Remove u_ref and (partial_s_ref - R) per source-qualified species before
//! arithmetic. Their bounded inventory residuals have separate mass budgets.
//! A bath contributes h_b-u_ref = Delta u(T_ref,T_b)+R*T_b and
//! partial_s_b-partial_s_ref+R. No absolute formation/entropy datum is summed.

use super::super::super::snapshot::MAX_SNAPSHOT_TRANSFERS;
use super::super::network::transport::{GasReservoirExchange, GasReservoirWindowAudit};
use super::{
    Accumulator, ClosedGasParcel, Error, GasEnergyError, MAX_SNAPSHOT_PARCELS,
    Nasa9StandardStateModelV1, R, SpeciesId, Temperature, amount, check, log_ratio,
    mixing_remainder, same_source, source_error, thermal_change,
};
use crate::thermochemical_gas::ThermochemicalGasModel;
use std::collections::{BTreeMap, BTreeSet};

struct Reference<'a> {
    model: &'a Nasa9StandardStateModelV1,
    catalog: &'a ThermochemicalGasModel,
    index: usize,
    temperature_k: f64,
    pressure_pa: f64,
    fraction: f64,
    mass_residual: Accumulator,
}

fn same_binding(
    a: &ThermochemicalGasModel,
    i: usize,
    b: &ThermochemicalGasModel,
    j: usize,
) -> bool {
    let source = &a.thermodynamics().components()[i];
    same_source(source, &b.thermodynamics().components()[j])
        && a.transport_materials()[i] == b.transport_materials()[j]
        && a.source_packs()
            .iter()
            .find(|p| &p.species == source.species())
            == b.source_packs()
                .iter()
                .find(|p| &p.species == source.species())
        && (a.water_species() == Some(source.species()))
            == (b.water_species() == Some(source.species()))
}

fn register<'a>(
    references: &mut BTreeMap<SpeciesId, Vec<Reference<'a>>>,
    catalog: &'a ThermochemicalGasModel,
    index: usize,
    temperature_k: f64,
    pressure_pa: f64,
    fraction: f64,
) -> Result<(), Error> {
    if [temperature_k, pressure_pa, fraction]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
    {
        return Err(Error::Invalid(
            "reservoir reference needs occupied supported species",
        ));
    }
    let model = &catalog.thermodynamics().components()[index];
    let group = references.entry(model.species().clone()).or_default();
    if !group
        .iter()
        .any(|r| same_binding(r.catalog, r.index, catalog, index))
    {
        group.push(Reference {
            model,
            catalog,
            index,
            temperature_k,
            pressure_pa,
            fraction,
            mass_residual: Accumulator::default(),
        });
    }
    Ok(())
}

fn potential(
    reference: &Reference<'_>,
    temperature_k: f64,
    pressure_pa: f64,
    fraction: f64,
    id: u64,
) -> Result<f64, Error> {
    // Separate ratios avoid underflow in the product x*p and preserve small
    // relative changes. No unoccupied species is sent to this logarithm.
    Ok(
        thermal_change(reference.model, reference.temperature_k, temperature_k, id)?
            - R * (log_ratio(pressure_pa, reference.pressure_pa)
                + log_ratio(fraction, reference.fraction)),
    )
}

fn internal_change(
    model: &Nasa9StandardStateModelV1,
    from: f64,
    to: f64,
    id: u64,
) -> Result<f64, Error> {
    model
        .molar_internal_energy_change(Temperature::new(from), Temperature::new(to))
        .map(|value| value.value())
        .map_err(|error| Error::Parcel {
            id,
            error: source_error(error),
        })
}

pub(in super::super) fn open_inventory_change<'a>(
    before: &'a BTreeMap<u64, ClosedGasParcel>,
    after: &'a BTreeMap<u64, ClosedGasParcel>,
    reservoirs: &[GasReservoirExchange<'a>],
    mass_tolerance_kg: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasReservoirWindowAudit, Error> {
    check(&mut cancelled)?;
    if before.len() > MAX_SNAPSHOT_PARCELS || reservoirs.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(Error::Invalid(
            "reservoir endpoint audit exceeds native work bound",
        ));
    }
    let mut work = reservoirs.len();
    for bath in reservoirs {
        work = work
            .checked_add(bath.species.len())
            .filter(|&n| n <= MAX_SNAPSHOT_TRANSFERS)
            .ok_or(Error::Invalid(
                "reservoir species audit exceeds native work bound",
            ))?;
    }
    let mut references: BTreeMap<SpeciesId, Vec<Reference<'a>>> = BTreeMap::new();
    let mut mass = Accumulator::default();
    let mut energy = Accumulator::default();
    let mut entropy = Accumulator::default();
    for (&id, old) in before {
        check(&mut cancelled)?;
        let next = &after[&id];
        if old.volume_m3() != next.volume_m3() {
            return Err(Error::Invalid(
                "reservoir audit requires unchanged gas volumes",
            ));
        }
        mass.add(next.mass_kg() - old.mass_kg())?;
        // Only original inventories and explicitly declared bath streams may
        // introduce a source. A final endpoint cannot authorize its own data,
        // even when its species mass is smaller than the numerical allowance.
        {
            let gas = old;
            let state = gas.state().state();
            for (index, &fraction) in gas
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions()
                .iter()
                .enumerate()
            {
                check(&mut cancelled)?;
                register(
                    &mut references,
                    gas.model(),
                    index,
                    state.temperature,
                    state.pressure,
                    fraction,
                )?;
            }
        }
    }

    // Immutable identity caches a completed source check for a bath shared by
    // many faces. Addresses never enter reports or numerical reduction order.
    let mut qualified_baths = BTreeSet::new();
    for bath in reservoirs {
        check(&mut cancelled)?;
        let old = before.get(&bath.parcel).ok_or(Error::Invalid(
            "reservoir exchange names an absent gas parcel",
        ))?;
        let next = &after[&bath.parcel];
        if !bath.reduced_heat_into_gas_j.is_finite() {
            return Err(Error::Invalid("reservoir heat must be finite"));
        }
        if qualified_baths.insert((
            core::ptr::from_ref(bath.catalog),
            core::ptr::from_ref(bath.state),
        )) {
            let point = bath
                .state
                .components()
                .first()
                .ok_or(Error::Invalid(
                    "reservoir state has no original source query",
                ))?
                .query();
            let resolved = bath
                .catalog
                .evaluate(&point)
                .map_err(|error| Error::Parcel {
                    id: bath.parcel,
                    error: GasEnergyError::Material(error),
                })?;
            if &resolved != bath.state {
                return Err(Error::Invalid(
                    "reservoir catalog substitutes its resolved state sources",
                ));
            }
        }
        let state = bath.state.state();
        energy.add(-bath.reduced_heat_into_gas_j)?;
        entropy.add(-bath.reduced_heat_into_gas_j / state.temperature)?;
        let mut seen = BTreeSet::new();
        for row in bath.species {
            check(&mut cancelled)?;
            if !row.mass_into_gas_kg.is_finite() || !seen.insert(&row.species) {
                return Err(Error::Invalid(
                    "reservoir streams need unique species and finite masses",
                ));
            }
            if row.mass_into_gas_kg == 0.0 {
                continue;
            }
            let index = bath
                .catalog
                .thermodynamics()
                .components()
                .iter()
                .position(|model| model.species() == &row.species)
                .ok_or(Error::Invalid(
                    "nonzero reservoir stream has zero bath species activity",
                ))?;
            let model = &bath.catalog.thermodynamics().components()[index];
            let fraction = bath.catalog.thermodynamics().mole_composition().fractions()[index];
            if fraction <= 0.0 {
                return Err(Error::Invalid(
                    "nonzero reservoir stream has zero bath species activity",
                ));
            }
            for gas in [old, next] {
                if let Some(local) = gas
                    .model()
                    .thermodynamics()
                    .components()
                    .iter()
                    .position(|candidate| candidate.species() == &row.species)
                {
                    if !same_binding(bath.catalog, index, gas.model(), local) {
                        return Err(Error::Invalid(
                            "reservoir stream substitutes an original gas species source",
                        ));
                    }
                }
            }
            register(
                &mut references,
                bath.catalog,
                index,
                state.temperature,
                state.pressure,
                fraction,
            )?;
            let reference = references
                .get_mut(&row.species)
                .unwrap()
                .iter_mut()
                .find(|r| same_binding(r.catalog, r.index, bath.catalog, index))
                .unwrap();
            let dn = row.mass_into_gas_kg / model.molar_mass().value();
            if !dn.is_finite() || dn == 0.0 {
                return Err(Error::Invalid(
                    "reservoir species amount is unrepresentable",
                ));
            }
            energy.add(
                -dn * internal_change(
                    model,
                    reference.temperature_k,
                    state.temperature,
                    bath.parcel,
                )?,
            )?;
            energy.add(-(dn * R) * state.temperature)?;
            entropy.add(
                -dn * potential(
                    reference,
                    state.temperature,
                    state.pressure,
                    fraction,
                    bath.parcel,
                )?,
            )?;
            entropy.add(-dn * R)?;
            reference.mass_residual.add(-row.mass_into_gas_kg)?;
            mass.add(-row.mass_into_gas_kg)?;
        }
    }

    for (&id, old) in before {
        let next = &after[&id];
        let initial = old.model().thermodynamics();
        let final_model = next.model().thermodynamics();
        let (mut i, mut j) = (0, 0);
        while i < initial.components().len() || j < final_model.components().len() {
            check(&mut cancelled)?;
            let ordering = match (initial.components().get(i), final_model.components().get(j)) {
                (Some(a), Some(b)) => a.species().cmp(b.species()),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => unreachable!(),
            };
            let (catalog, index, n0, n1, x0, x1) = match ordering {
                std::cmp::Ordering::Less => {
                    let row = (
                        old.model(),
                        i,
                        amount(old, i)?,
                        0.0,
                        initial.mole_composition().fractions()[i],
                        0.0,
                    );
                    i += 1;
                    row
                }
                std::cmp::Ordering::Greater => {
                    let row = (
                        next.model(),
                        j,
                        0.0,
                        amount(next, j)?,
                        0.0,
                        final_model.mole_composition().fractions()[j],
                    );
                    j += 1;
                    row
                }
                std::cmp::Ordering::Equal => {
                    if !same_binding(old.model(), i, next.model(), j) {
                        return Err(Error::Invalid(
                            "reservoir endpoint changed an original species source",
                        ));
                    }
                    let row = (
                        old.model(),
                        i,
                        amount(old, i)?,
                        amount(next, j)?,
                        initial.mole_composition().fractions()[i],
                        final_model.mole_composition().fractions()[j],
                    );
                    i += 1;
                    j += 1;
                    row
                }
            };
            let model = &catalog.thermodynamics().components()[index];
            let reference = references
                .get_mut(model.species())
                .and_then(|group| {
                    group
                        .iter_mut()
                        .find(|r| same_binding(r.catalog, r.index, catalog, index))
                })
                .ok_or(Error::Invalid(
                    "reservoir endpoint introduced an undeclared species source",
                ))?;
            let t0 = old.state().state().temperature;
            let t1 = next.state().state().temperature;
            if n0 == 0.0 {
                energy.add(n1 * internal_change(model, reference.temperature_k, t1, id)?)?;
                entropy.add(
                    n1 * (potential(reference, t1, next.state().state().pressure, x1, id)? + R),
                )?;
            } else if n1 == 0.0 {
                energy.add(-n0 * internal_change(model, reference.temperature_k, t0, id)?)?;
                entropy.add(
                    -n0 * (potential(reference, t0, old.state().state().pressure, x0, id)? + R),
                )?;
            } else {
                energy.add(n1 * internal_change(model, t0, t1, id)?)?;
                if n1 != n0 {
                    energy.add(
                        (n1 - n0) * internal_change(model, reference.temperature_k, t0, id)?,
                    )?;
                    entropy.add(
                        (n1 - n0) * potential(reference, t0, old.state().state().pressure, x0, id)?,
                    )?;
                }
                entropy.add(n1 * (thermal_change(model, t0, t1, id)? - R * log_ratio(t1, t0)))?;
                entropy.add(-R * mixing_remainder(n0, n1))?;
            }
            reference
                .mass_residual
                .add((n1 - n0) * model.molar_mass().value())?;
        }
    }

    let mut maximum_species_mass_residual_kg = 0.0_f64;
    for group in references.into_values() {
        let mut species_mass = Accumulator::default();
        for reference in group {
            check(&mut cancelled)?;
            let residual = reference.mass_residual.total()?;
            maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(residual.abs());
            species_mass.add(residual)?;
            if residual.abs() > mass_tolerance_kg {
                return Err(Error::Invalid(
                    "reservoir window violated a source-qualified species budget",
                ));
            }
        }
        let residual = species_mass.total()?.abs();
        maximum_species_mass_residual_kg = maximum_species_mass_residual_kg.max(residual);
        if residual > mass_tolerance_kg {
            return Err(Error::Invalid(
                "reservoir window violated a combined species budget",
            ));
        }
    }
    check(&mut cancelled)?;
    Ok(GasReservoirWindowAudit {
        total_mass_residual_kg: mass.total()?,
        maximum_species_mass_residual_kg,
        energy_residual_j: energy.total()?,
        entropy_production_j_k: entropy.total()?,
    })
}
