//! External phase streams enter the SAME gas snapshot budget as internal edges.
//!
//! The energy coordinate is relative to the original gas's partial enthalpy:
//! E_into_gas = net_mass*h_species(T_original) + relative_energy. This avoids
//! subtracting large formation energies to recover a small thermal increment.
//! A boundary supplies rates and a retained source catalog, not another gas EOS.
//! External phases must accept their opposite mass/energy change atomically.

use super::*;

/// One species' integrated boundary streams. Gross in/out and the authoritative
/// net inventory increment are separate so opposing streams cannot hide uptake.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundarySpeciesTransfer {
    pub species: SpeciesId,
    pub incoming_mass_kg: f64,
    pub outgoing_mass_kg: f64,
    pub net_mass_into_gas_kg: f64,
    /// Energy beyond net_mass*h_species at this gas's ORIGINAL temperature [J].
    pub relative_energy_into_gas_j: f64,
}

/// One external physical owner. IDs are unique within the boundary list.
/// Several owners may exchange with the same gas on the same snapshot.
#[derive(Clone, Copy, Debug)]
pub struct SpeciesNetworkBoundary<'a> {
    pub id: u64,
    pub parcel: u64,
    pub catalog: &'a ThermochemicalGasModel,
    pub species: &'a [BoundarySpeciesTransfer],
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkBoundaryFlux {
    pub id: u64,
    pub parcel: u64,
    pub transfer: BoundarySpeciesTransfer,
    pub original_partial_enthalpy_j_kg: f64,
    pub energy_into_gas_j: f64,
}

pub(super) struct BoundaryAccounts<'a> {
    pub catalogs: BTreeMap<u64, Vec<&'a ThermochemicalGasModel>>,
    pub fluxes: Vec<SpeciesNetworkBoundaryFlux>,
    pub species: BTreeMap<SpeciesId, Sum>,
    pub mass: Sum,
}

pub(super) fn work_count(boundaries: &[SpeciesNetworkBoundary<'_>]) -> Result<usize, Error> {
    boundaries.iter().try_fold(boundaries.len(), |total, row| {
        total.checked_add(row.species.len()).filter(|&n| n <= MAX_SNAPSHOT_TRANSFERS)
            .ok_or(Error::Invalid("boundary owner/species count exceeds shared transport work bound"))
    })
}

/// Insert every external exchange before the first caloric inversion. Original
/// catalogs are joined by the parent owner, which also rejects source conflicts.
pub(super) fn accumulate<'a>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    boundaries: &[SpeciesNetworkBoundary<'a>],
    budgets: &mut BTreeMap<u64, Budget>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BoundaryAccounts<'a>, Error> {
    let mut ordered = BTreeMap::new();
    for row in boundaries {
        check(cancelled)?;
        if !parcels.contains_key(&row.parcel) || row.species.is_empty()
            || ordered.insert(row.id, row).is_some() {
            return Err(Error::Invalid("boundary owners need unique IDs, an existing gas and explicit species rows"));
        }
    }
    let mut result = BoundaryAccounts { catalogs: BTreeMap::new(), fluxes: Vec::new(),
        species: BTreeMap::new(), mass: Sum::default() };
    let mut gross_residuals: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    let mut gross_total = Sum::default();
    for (id, owner) in ordered {
        check(cancelled)?;
        let gas = &parcels[&owner.parcel];
        let map = |e| material(owner.parcel, e);
        result.catalogs.entry(owner.parcel).or_default().push(owner.catalog);
        let mut species = BTreeMap::new();
        for row in owner.species {
            check(cancelled)?;
            if species.insert(&row.species, row).is_some()
                || !row.incoming_mass_kg.is_finite() || row.incoming_mass_kg < 0.0
                || !row.outgoing_mass_kg.is_finite() || row.outgoing_mass_kg < 0.0
                || !row.net_mass_into_gas_kg.is_finite()
                || !row.relative_energy_into_gas_j.is_finite() {
                return Err(Error::Invalid("boundary streams need unique species and finite nonnegative gross masses"));
            }
        }
        for row in species.into_values() {
            check(cancelled)?;
            let residual = (row.incoming_mass_kg - row.outgoing_mass_kg) - row.net_mass_into_gas_kg;
            if !residual.is_finite() || residual.abs() > controls.mass_tolerance_kg
                || (row.net_mass_into_gas_kg > 0.0 && row.incoming_mass_kg == 0.0)
                || (row.net_mass_into_gas_kg < 0.0 && row.outgoing_mass_kg == 0.0) {
                return Err(Error::Invalid("boundary net inventory does not agree with its gross directional streams"));
            }
            gross_residuals.entry(row.species.clone()).or_default().add(residual).map_err(map)?;
            gross_total.add(residual).map_err(map)?;
            let source = owner.catalog.thermodynamics().components().iter()
                .find(|s| s.species() == &row.species)
                .ok_or(Error::Invalid("boundary species is absent from its retained source catalog"))?;
            let h = source.evaluate_properties(Temperature::new(gas.state().state().temperature))
                .map_err(|e| map(source_error(e)))?.h().value() / source.molar_mass().value();
            let energy = row.net_mass_into_gas_kg * h + row.relative_energy_into_gas_j;
            if !h.is_finite() || !energy.is_finite() {
                return Err(Error::Invalid("boundary original partial enthalpy or energy report is unrepresentable"));
            }
            let budget = budgets.entry(owner.parcel).or_default();
            budget.boundary_mass_in.add(row.incoming_mass_kg).map_err(map)?;
            budget.boundary_mass_out.add(row.outgoing_mass_kg).map_err(map)?;
            budget.boundary_delta_mass.add(row.net_mass_into_gas_kg).map_err(map)?;
            if row.outgoing_mass_kg > 0.0 {
                budget.withdrawn.entry(row.species.clone()).or_default()
                    .add(row.outgoing_mass_kg).map_err(map)?;
            }
            if row.net_mass_into_gas_kg != 0.0 {
                budget.species.entry(row.species.clone()).or_default()
                    .add(row.net_mass_into_gas_kg).map_err(map)?;
                result.species.entry(row.species.clone()).or_default()
                    .add(row.net_mass_into_gas_kg).map_err(map)?;
            }
            budget.sensible_enthalpy.add(row.relative_energy_into_gas_j).map_err(map)?;
            budget.relative_energy.add(row.relative_energy_into_gas_j).map_err(map)?;
            result.mass.add(row.net_mass_into_gas_kg).map_err(map)?;
            result.fluxes.push(SpeciesNetworkBoundaryFlux { id, parcel: owner.parcel,
                transfer: row.clone(), original_partial_enthalpy_j_kg: h, energy_into_gas_j: energy });
        }
    }
    if gross_total.value.abs() > controls.mass_tolerance_kg
        || gross_residuals.values().any(|s| s.value.abs() > controls.mass_tolerance_kg) {
        return Err(Error::Invalid("aggregate boundary gross/net species balance failed"));
    }
    check(cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
