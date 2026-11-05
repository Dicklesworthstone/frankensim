//! Explicit external species/phase streams on the existing NASA reference chart.
//!
//! `E_into = net_mass * h_ideal_species(T_original) + relative_energy`.
//! This is deliberately the SAME coordinate as the ideal network boundary:
//! the external phase authors its complete energy, including any declared
//! uptake law. A negative row does not silently acquire a gas PR departure.
//! `RealGasBoundaryStream` supplies actual upstream gas partial enthalpy when
//! that is the intended transport law, without subtracting formation energies.

use super::snapshot::Budget;
use super::{
    Error, GasEnergyError, IsochoricRealGasParcel, Parcels, SpeciesNetworkSolve, Sum, check,
    material,
};
use crate::thermochemical_gas::energy::parcel::transfer::snapshot::MAX_SNAPSHOT_TRANSFERS;
use crate::thermochemical_gas::failure;
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_qty::Temperature;
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

pub use crate::thermochemical_gas::energy::parcel::transfer::species::network::boundary::{
    BoundarySpeciesTransfer, SpeciesNetworkBoundaryFlux,
};

/// One external physical owner. Its complete retained source binds incoming
/// NASA/transport data, critical constants and known PR interactions. This
/// catalog supplies no external phase EOS, latent heat or finite inventory.
#[derive(Clone, Copy, Debug)]
pub struct RealGasSpeciesNetworkBoundary<'a> {
    pub id: u64,
    pub parcel: u64,
    pub catalog: &'a ThermochemicalPengRobinsonModel,
    /// Complete authored energy in the original IDEAL species enthalpy chart.
    pub species: &'a [BoundarySpeciesTransfer],
}

/// One actual gas stream with its borrowed original donor source. Its external
/// owner must debit the opposite mass/energy and admit the donor's aggregate
/// withdrawal. The receiver passed to `inflow` must be the same original
/// snapshot subsequently supplied to the network; it sets the energy chart.
#[derive(Clone, Debug)]
pub struct RealGasBoundaryStream<'a> {
    catalog: &'a ThermochemicalPengRobinsonModel,
    transfer: BoundarySpeciesTransfer,
}
impl<'a> RealGasBoundaryStream<'a> {
    /// An authored positive incoming mass at the ORIGINAL donor gas's actual
    /// partial enthalpy, including its original PR departure and flow work.
    pub fn inflow(
        donor: &'a IsochoricRealGasParcel,
        receiver: &IsochoricRealGasParcel,
        species: &SpeciesId,
        mass_kg: f64,
    ) -> Result<Self, GasEnergyError> {
        Self::gas_stream(
            donor,
            receiver.state().temperature_k(),
            species,
            mass_kg,
            1.0,
        )
    }

    /// Actual outflow at this original gas's partial enthalpy. Use a raw
    /// `BoundarySpeciesTransfer` instead for a separately declared phase
    /// uptake energy; the network never replaces that law with gas enthalpy.
    pub fn outflow(
        gas: &'a IsochoricRealGasParcel,
        species: &SpeciesId,
        mass_kg: f64,
    ) -> Result<Self, GasEnergyError> {
        Self::gas_stream(gas, gas.state().temperature_k(), species, mass_kg, -1.0)
    }

    fn gas_stream(
        donor: &'a IsochoricRealGasParcel,
        reference_temperature_k: f64,
        species: &SpeciesId,
        mass_kg: f64,
        sign: f64,
    ) -> Result<Self, GasEnergyError> {
        let invalid = super::super::invalid;
        if !mass_kg.is_finite() || mass_kg < 0.0 {
            return Err(invalid(
                "gas boundary stream requires finite nonnegative mass",
            ));
        }
        let (index, source) = donor
            .model()
            .source()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
            .find(|(_, source)| source.species() == species)
            .ok_or(invalid(
                "gas boundary stream needs its original donor species",
            ))?;
        let mut relative = Sum::default();
        if mass_kg > 0.0 {
            if donor.model().mechanics().mole_fractions()[index] <= 0.0 {
                return Err(invalid(
                    "absent original donor species cannot supply a gas boundary stream",
                ));
            }
            let partial = donor.model().partial_calorics_at_density(
                donor.density_kg_m3(),
                donor.state().temperature_k(),
            )?;
            let change = source
                .molar_enthalpy_change(
                    Temperature::new(reference_temperature_k),
                    Temperature::new(donor.state().temperature_k()),
                )
                .map_err(|e| {
                    GasEnergyError::from(failure("real-gas boundary reference enthalpy", e))
                })?
                .value();
            let amount = sign * mass_kg / source.molar_mass().value();
            if !amount.is_finite() || amount == 0.0 {
                return Err(invalid(
                    "gas boundary stream molar amount is unrepresentable",
                ));
            }
            relative.add(amount * change)?;
            relative
                .add(amount * partial.departure().components()[index].residual_enthalpy_j_mol())?;
        }
        Ok(Self {
            catalog: donor.model(),
            transfer: BoundarySpeciesTransfer {
                species: species.clone(),
                incoming_mass_kg: if sign > 0.0 { mass_kg } else { 0.0 },
                outgoing_mass_kg: if sign < 0.0 { mass_kg } else { 0.0 },
                net_mass_into_gas_kg: sign * mass_kg,
                relative_energy_into_gas_j: relative.value,
            },
        })
    }

    pub const fn transfer(&self) -> &BoundarySpeciesTransfer {
        &self.transfer
    }
    pub const fn catalog(&self) -> &'a ThermochemicalPengRobinsonModel {
        self.catalog
    }
    pub fn boundary(&self, id: u64, parcel: u64) -> RealGasSpeciesNetworkBoundary<'_> {
        RealGasSpeciesNetworkBoundary {
            id,
            parcel,
            catalog: self.catalog,
            species: std::slice::from_ref(&self.transfer),
        }
    }
}

pub(super) struct BoundaryAccounts<'a> {
    pub catalogs: BTreeMap<u64, Vec<&'a ThermochemicalPengRobinsonModel>>,
    rows: Vec<(
        u64,
        u64,
        &'a ThermochemicalPengRobinsonModel,
        BoundarySpeciesTransfer,
    )>,
    pub species: BTreeMap<SpeciesId, Sum>,
    pub mass: Sum,
}

pub(super) fn work_count(boundaries: &[RealGasSpeciesNetworkBoundary<'_>]) -> Result<usize, Error> {
    boundaries.iter().try_fold(boundaries.len(), |total, row| {
        total
            .checked_add(row.species.len())
            .filter(|&n| n <= MAX_SNAPSHOT_TRANSFERS)
            .ok_or(Error::Invalid(
                "real-gas boundary owner/species count exceeds shared work bound",
            ))
    })
}

// Only schema, inventory and source-presence checks here. No caloric query or
// source join precedes the combined original-inventory withdrawal admission.
pub(super) fn accumulate<'a>(
    parcels: &Parcels<'_>,
    boundaries: &[RealGasSpeciesNetworkBoundary<'a>],
    budgets: &mut BTreeMap<u64, Budget>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BoundaryAccounts<'a>, Error> {
    let mut ordered = BTreeMap::new();
    for owner in boundaries {
        check(cancelled)?;
        if !parcels.contains_key(&owner.parcel)
            || owner.species.is_empty()
            || ordered.insert(owner.id, owner).is_some()
        {
            return Err(Error::Invalid(
                "real-gas boundaries need unique IDs, existing nodes and explicit species",
            ));
        }
    }
    let mut result = BoundaryAccounts {
        catalogs: BTreeMap::new(),
        rows: Vec::new(),
        species: BTreeMap::new(),
        mass: Sum::default(),
    };
    let (mut gross_total, mut gross_species) = (Sum::default(), BTreeMap::<SpeciesId, Sum>::new());
    for (id, owner) in ordered {
        check(cancelled)?;
        let map = |e| material(owner.parcel, e);
        let mut rows = BTreeMap::new();
        for row in owner.species {
            check(cancelled)?;
            if rows.insert(&row.species, row).is_some()
                || !row.incoming_mass_kg.is_finite()
                || row.incoming_mass_kg < 0.0
                || !row.outgoing_mass_kg.is_finite()
                || row.outgoing_mass_kg < 0.0
                || !row.net_mass_into_gas_kg.is_finite()
                || !row.relative_energy_into_gas_j.is_finite()
            {
                return Err(Error::Invalid(
                    "real-gas boundary needs unique species, finite energy and nonnegative gross mass",
                ));
            }
        }
        let mut active = false;
        for row in rows.into_values() {
            check(cancelled)?;
            let residual = (row.incoming_mass_kg - row.outgoing_mass_kg) - row.net_mass_into_gas_kg;
            if !residual.is_finite()
                || residual.abs() > controls.mass_tolerance_kg
                || (row.net_mass_into_gas_kg > 0.0 && row.incoming_mass_kg == 0.0)
                || (row.net_mass_into_gas_kg < 0.0 && row.outgoing_mass_kg == 0.0)
            {
                return Err(Error::Invalid(
                    "real-gas boundary net inventory disagrees with gross streams",
                ));
            }
            if !owner
                .catalog
                .source()
                .thermodynamics()
                .components()
                .iter()
                .any(|s| s.species() == &row.species)
            {
                return Err(Error::Invalid(
                    "real-gas boundary species is absent from its retained catalog",
                ));
            }
            gross_total.add(residual).map_err(map)?;
            gross_species
                .entry(row.species.clone())
                .or_default()
                .add(residual)
                .map_err(map)?;
            if row.incoming_mass_kg == 0.0
                && row.outgoing_mass_kg == 0.0
                && row.relative_energy_into_gas_j == 0.0
            {
                continue;
            }
            active = true;
            let budget = budgets.entry(owner.parcel).or_default();
            budget.mass_in.add(row.incoming_mass_kg).map_err(map)?;
            budget.mass_out.add(row.outgoing_mass_kg).map_err(map)?;
            if row.outgoing_mass_kg > 0.0 {
                budget
                    .withdrawn
                    .entry(row.species.clone())
                    .or_default()
                    .add(row.outgoing_mass_kg)
                    .map_err(map)?;
            }
            if row.net_mass_into_gas_kg != 0.0 {
                budget
                    .species
                    .entry(row.species.clone())
                    .or_default()
                    .add(row.net_mass_into_gas_kg)
                    .map_err(map)?;
                result
                    .species
                    .entry(row.species.clone())
                    .or_default()
                    .add(row.net_mass_into_gas_kg)
                    .map_err(map)?;
            }
            budget
                .relative_stream
                .add(row.relative_energy_into_gas_j)
                .map_err(map)?;
            result.mass.add(row.net_mass_into_gas_kg).map_err(map)?;
            result
                .rows
                .push((id, owner.parcel, owner.catalog, row.clone()));
        }
        if active {
            result
                .catalogs
                .entry(owner.parcel)
                .or_default()
                .push(owner.catalog);
        }
    }
    if gross_total.value.abs() > controls.mass_tolerance_kg
        || gross_species
            .values()
            .any(|sum| sum.value.abs() > controls.mass_tolerance_kg)
    {
        return Err(Error::Invalid(
            "aggregate real-gas boundary gross/net species balance failed",
        ));
    }
    Ok(result)
}

impl BoundaryAccounts<'_> {
    // Reports use the explicit ideal-reference baseline. PR energy supplied by
    // actual-gas convenience streams is already in `relative_energy`, once.
    pub(super) fn fluxes(
        &self,
        parcels: &Parcels<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<SpeciesNetworkBoundaryFlux>, Error> {
        self.rows
            .iter()
            .map(|(id, parcel, catalog, transfer)| {
                check(cancelled)?;
                let source = catalog
                    .source()
                    .thermodynamics()
                    .components()
                    .iter()
                    .find(|source| source.species() == &transfer.species)
                    .expect("admitted original boundary species");
                let h = source
                    .evaluate_properties(Temperature::new(parcels[parcel].state().temperature_k()))
                    .map_err(|e| {
                        material(
                            *parcel,
                            failure("real-gas boundary reference report", e).into(),
                        )
                    })?
                    .h()
                    .value()
                    / source.molar_mass().value();
                let energy =
                    transfer.net_mass_into_gas_kg * h + transfer.relative_energy_into_gas_j;
                if !h.is_finite() || !energy.is_finite() {
                    return Err(Error::Invalid(
                        "real-gas boundary complete energy report is unrepresentable",
                    ));
                }
                Ok(SpeciesNetworkBoundaryFlux {
                    id: *id,
                    parcel: *parcel,
                    transfer: transfer.clone(),
                    original_partial_enthalpy_j_kg: h,
                    energy_into_gas_j: energy,
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests;
