//! Source-preserving species union and reference-free PR receiving energy.

use super::{GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, positive, updated_point};
use crate::thermochemical_gas::energy::parcel::transfer::{
    GasSpeciesTransfer, inventory, mixed_model_from_portions,
};
use crate::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::PengRobinsonComponentV1;
use std::collections::BTreeMap;

fn invalid(message: &'static str) -> GasEnergyError {
    GasEnergyError::Invalid(message)
}

/// The existing source join owns NASA, transport, water and pack compatibility.
/// This adds only the original PR critical constants and KNOWN pair entries.
/// A zero interaction inside one catalog supplies no missing cross-catalog data.
pub(super) fn receiver_model(
    donor: &IsochoricRealGasParcel,
    receiver: &IsochoricRealGasParcel,
    mass_kg: f64,
    receiver_mass_kg: f64,
) -> Result<ThermochemicalPengRobinsonModel, GasEnergyError> {
    let source = mixed_model_from_portions(
        &[
            (receiver.model().source(), receiver.mass_kg()),
            (donor.model().source(), mass_kg),
        ],
        receiver_mass_kg,
    )?;
    let mut critical: BTreeMap<SpeciesId, PengRobinsonComponentV1> = BTreeMap::new();
    let mut pairs = BTreeMap::new();
    for model in [receiver.model(), donor.model()] {
        let species = model.source().thermodynamics().components();
        let mechanics = model.mechanics();
        for (i, (species_i, component)) in species.iter().zip(mechanics.components()).enumerate() {
            let id = species_i.species();
            if let Some(previous) = critical.insert(id.clone(), *component) {
                if previous != *component {
                    return Err(invalid(
                        "shared gas species has conflicting Peng-Robinson critical properties",
                    ));
                }
            }
            for (j, species_j) in species.iter().enumerate() {
                let pair = (id.clone(), species_j.species().clone());
                let interaction = mechanics.binary_interactions()[i * species.len() + j];
                if let Some(previous) = pairs.insert(pair, interaction) {
                    if previous != interaction {
                        return Err(invalid(
                            "shared gas species pair has conflicting Peng-Robinson interaction data",
                        ));
                    }
                }
            }
        }
    }
    let species = source.thermodynamics().components();
    let mut bindings = Vec::with_capacity(species.len());
    let mut interactions = Vec::with_capacity(species.len() * species.len());
    for row in species {
        let component = critical.get(row.species()).ok_or(invalid(
            "mixed gas lost an original Peng-Robinson species binding",
        ))?;
        bindings.push(PengRobinsonSpeciesBinding {
            species: row.species().clone(),
            critical_temperature: component.critical_temperature,
            critical_pressure: component.critical_pressure,
            acentric_factor: component.acentric_factor,
        });
        for column in species {
            interactions.push(*pairs.get(&(row.species().clone(), column.species().clone()))
                .ok_or(invalid("mixed gas needs explicit Peng-Robinson interaction data for every new species pair"))?);
        }
    }
    Ok(ThermochemicalPengRobinsonModel::new(
        source,
        bindings,
        interactions,
    )?)
}

/// Solve the receiving state using the ORIGINAL increment/Cv/root owner.
/// At receiver T0, conserved ideal species energies cancel before arithmetic.
/// Only ideal temperature increments, PR departures and upstream p/rho remain.
pub(super) fn receiver_endpoint(
    donor: &IsochoricRealGasParcel,
    old: &IsochoricRealGasParcel,
    model: ThermochemicalPengRobinsonModel,
    mass_kg: f64,
    received_mass_kg: f64,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(IsochoricRealGasParcel, f64, usize), GasEnergyError> {
    let density = mass_kg / old.volume_m3();
    if !positive(density) || density == old.density_kg_m3() {
        return Err(invalid(
            "mixed real-gas receiving density is unrepresentable",
        ));
    }
    let temperature = old.state().temperature_k();
    // The reference is the NEW composition at its NEW density and the old
    // receiver temperature. Use the original lightweight caloric owner here;
    // reconstruction below admits the complete anchor and its material cards.
    let anchor_departure = model
        .caloric_at_density(density, temperature)?
        .residual_energy_j_kg;
    let departure = |gas: &IsochoricRealGasParcel| {
        gas.state().root().residual_internal_energy_j_mol()
            / gas.model().mechanics().molar_mass().value()
    };
    let donor_sensible = donor
        .model()
        .source()
        .thermodynamics()
        .specific_internal_energy_change(
            Temperature::new(temperature),
            Temperature::new(donor.state().temperature_k()),
        )
        .map_err(|error| {
            crate::thermochemical_gas::failure("real-gas mixing caloric increment", error)
        })?
        .value();
    let target = (old.mass_kg() / mass_kg) * (departure(old) - anchor_departure)
        + (received_mass_kg / mass_kg)
            * (donor_sensible + departure(donor) - anchor_departure
                + donor.state().pressure_pa() / donor.density_kg_m3());
    if !target.is_finite() {
        return Err(invalid(
            "mixed real-gas reference-free energy target is unrepresentable",
        ));
    }
    let solved = model.reconstruct_density_internal_energy_increment(
        old.point(),
        density,
        target,
        controls,
        &mut *cancelled,
    )?;
    let residual = solved.residual_j_kg();
    let iterations = solved.iterations();
    let state = solved.into_state();
    let point = updated_point(old.point(), state.temperature_k(), state.pressure_pa())?;
    let mut next = old.clone();
    next.model = model;
    next.point = point;
    next.state = state;
    next.mass_kg = mass_kg;
    next.density_kg_m3 = density;
    Ok((next, residual, iterations))
}

/// Audit actual endpoint inventories against EACH upstream species flux.
/// Formation datums cancel in the intended ledger; rounding of conserved
/// species amounts has this separate absolute physical mass budget.
pub(super) fn species_balance(
    old_donor: &IsochoricRealGasParcel,
    old_receiver: &IsochoricRealGasParcel,
    donor: &IsochoricRealGasParcel,
    receiver: &IsochoricRealGasParcel,
    transferred_mass_kg: f64,
) -> Result<(Vec<GasSpeciesTransfer>, f64), GasEnergyError> {
    let before = [
        inventory(old_donor.model().source(), old_donor.mass_kg())?,
        inventory(old_receiver.model().source(), old_receiver.mass_kg())?,
    ];
    let after = [
        inventory(donor.model().source(), donor.mass_kg())?,
        inventory(receiver.model().source(), receiver.mass_kg())?,
    ];
    let moved = inventory(old_donor.model().source(), transferred_mass_kg)?;
    let mut maximum = 0.0_f64;
    for id in before[0].keys().chain(before[1].keys()) {
        let transferred = moved.get(id).copied().unwrap_or(0.0);
        let donor_change =
            after[0].get(id).copied().unwrap_or(0.0) - before[0].get(id).copied().unwrap_or(0.0);
        let receiver_change =
            after[1].get(id).copied().unwrap_or(0.0) - before[1].get(id).copied().unwrap_or(0.0);
        for defect in [
            donor_change + transferred,
            receiver_change - transferred,
            donor_change + receiver_change,
        ] {
            if !defect.is_finite() {
                return Err(invalid(
                    "real-gas species inventory defect is unrepresentable",
                ));
            }
            maximum = maximum.max(defect.abs());
        }
    }
    Ok((
        moved
            .into_iter()
            .map(|(species, mass_kg)| GasSpeciesTransfer { species, mass_kg })
            .collect(),
        maximum,
    ))
}
