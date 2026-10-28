//! Simultaneous selective transport on a graph of original material inventories.
//! All edges use original donor enthalpy; combine every flux before rebuilding
//! an endpoint. No newly received species can leave again within this interval.
/// Pressure restrictions and molecular mixing share one inventory update.
pub mod transport;
/// External phases share the original gas budgets and caloric reconstruction.
pub mod boundary;

use super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use super::{
    ClosedGasParcel, GasEnergyError, GasEnergySolve, SpeciesId, Sum, ThermochemicalGasModel,
    inventory, mixed_model_from_portions, model_for_inventory, positive, source_enthalpy,
    source_error,
};
use crate::thermochemical_gas::energy::updated_point;
use fs_qty::Temperature;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::collections::{BTreeMap, BTreeSet};
use boundary::{SpeciesNetworkBoundary, SpeciesNetworkBoundaryFlux};

/// Nonnegative integrated mass of ONE species, not the whole donor mixture.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkTransfer {
    pub donor: u64,
    pub receiver: u64,
    pub species: SpeciesId,
    pub mass_kg: f64,
}
/// Antisymmetric integrated energy in addition to the original species enthalpy.
/// Positive heat leaves parcels[0] and enters parcels[1]. Prescribing this
/// correction does not itself establish a physical heat/work law or entropy.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkHeatTransfer {
    pub parcels: [u64; 2],
    pub heat_left_to_right_j: f64,
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
    Parcel {
        id: u64,
        error: GasEnergyError,
    },
    Withdrawal {
        id: u64,
        species: SpeciesId,
        fraction: f64,
        maximum: f64,
    },
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
    /// Net original species enthalpy, using grouped species mass changes and
    /// source-relative incoming enthalpy; rounded gross reports are not subtracted.
    pub net_enthalpy_j: f64,
    /// Additional signed energy received, separate from species enthalpy.
    pub net_additional_heat_j: f64,
    /// Intended species-reference energy change plus the actual original-law
    /// sensible storage. Actual inventory rounding has a separate mass budget.
    pub energy_change_j: f64,
    /// Original-caloric closure after reference cancellation, equivalent to
    /// Delta U minus net species enthalpy minus additional heat.
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub iterations: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkStep {
    pub fluxes: Vec<SpeciesNetworkFlux>,
    /// Nonzero combined edges in ascending parcel-pair order; endpoints are
    /// canonical ascending IDs and the retained heat may have either sign.
    pub heat: Vec<SpeciesNetworkHeatTransfer>,
    /// External source streams in owner-ID/species order. Opposite phase
    /// storage must be accounted by the enclosing coupled owner.
    pub boundaries: Vec<SpeciesNetworkBoundaryFlux>,
    pub parcels: Vec<SpeciesNetworkParcelStep>,
    /// Actual gas mass change minus the declared external net mass.
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Sum of local original-caloric residuals after conserved species datums
    /// cancel in the intended paired ledger. Actual inventory rounding is
    /// admitted separately by the total/individual mass residual budgets.
    pub total_energy_residual_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: SpeciesNetworkStep,
}
impl SpeciesNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &SpeciesNetworkStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, SpeciesNetworkStep) {
        (self.parcels, self.step)
    }
}
use SpeciesNetworkError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn material(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Parcel { id, error }
    }
}
#[derive(Default)]
struct Budget {
    incoming: BTreeSet<u64>,
    mass_in: Sum,
    mass_out: Sum,
    boundary_mass_in: Sum,
    boundary_mass_out: Sum,
    boundary_delta_mass: Sum,
    sensible_enthalpy: Sum,
    heat: Sum,
    relative_energy: Sum,
    species: BTreeMap<SpeciesId, Sum>,
    withdrawn: BTreeMap<SpeciesId, Sum>,
}

// The original incremental caloric owner cancels NASA formation constants
// BEFORE evaluating changes. Rebuild the same gas parcel exactly once.
fn relative_endpoint(
    old: &ClosedGasParcel,
    model: ThermochemicalGasModel,
    mass: f64,
    thermal_energy: f64,
    controls: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, usize), GasEnergyError> {
    let change = thermal_energy / mass;
    if !change.is_finite() || (thermal_energy != 0.0 && change == 0.0) {
        return Err(GasEnergyError::Invalid(
            "network sensible energy increment is unrepresentable",
        ));
    }
    let solved = model.reconstruct_density_internal_energy_increment(
        old.point(),
        mass / old.volume_m3(),
        change,
        controls,
        cancelled,
    )?;
    let mut parcel = old.parcel.clone();
    parcel.point = updated_point(
        old.point(),
        solved.state().state().temperature,
        solved.state().state().pressure,
    )?;
    parcel.state = solved.state().clone();
    parcel.model = model;
    parcel.mass_kg = mass;
    Ok((ClosedGasParcel::from_isochoric(parcel), solved.iterations()))
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
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    exchange_candidate_with_heat(parcels, transfers, &[], endpoints, controls, cancelled)
}

/// Reconstruct each affected gas once from the same original species streams
/// and antisymmetric additional energy. Heat-only nodes use this same owner;
/// no already-mixed candidate receives a second temperature/pressure update.
///
/// Signed heat edges are normalized to ascending endpoint IDs, sorted by pair
/// and signed value, then combined with compensated summation. Their joint row
/// count with species streams is bounded before allocation. Gross original
/// species withdrawal admission is unchanged; heat cannot conceal depletion.
/// Original species enthalpy and additional energy are reported independently.
/// Per-species net mass increments multiply reference energies only after all
/// gross streams are combined; temperature differences use original NASA
/// relative integrals. Incoming species must support the receiver's original
/// temperature and the declared common caloric bracket. No source extrapolation
/// is added to make a newly mixed state admissible.
pub fn exchange_candidate_with_heat(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    exchange_candidate_with_boundaries(parcels, transfers, heat, &[], endpoints, controls, cancelled)
}

/// Joint gas/gas and external-phase exchange on one original snapshot.
/// Every boundary's gross uptake joins every internal donor allowance BEFORE
/// one original-caloric inversion per gas. Incoming/desorbed mass cannot fund
/// outgoing flow. Gas-only entropy is not asserted for this open subsystem;
/// an enclosing owner must audit its actual finite external phase endpoints.
pub fn exchange_candidate_with_boundaries(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    boundaries: &[SpeciesNetworkBoundary<'_>],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    exchange_bounded(parcels, transfers, heat, boundaries, endpoints, controls,
        MAX_SNAPSHOT_TRANSFERS, cancelled)
}

// A bounded bulk edge expands into at most the existing active-species bound.
// Keep its 4096-edge capacity without relaxing the public selective API's
// separately declared work limit. Snapshot lowering admits bulk rows first.
pub(in crate::thermochemical_gas::energy::parcel::transfer) fn exchange_bulk_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    exchange_bounded(parcels, transfers, &[], &[], endpoints, controls,
        MAX_SNAPSHOT_TRANSFERS * super::MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn exchange_bounded(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    boundaries: &[SpeciesNetworkBoundary<'_>],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    maximum_transfers: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpeciesNetworkCandidate, Error> {
    check(&mut cancelled)?;
    let boundary_work = boundary::work_count(boundaries)?;
    if parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || transfers
            .len()
            .checked_add(heat.len())
            .and_then(|n| n.checked_add(boundary_work))
            .is_none_or(|n| n > maximum_transfers)
        || !parcels.keys().eq(endpoints.keys())
        || !positive(controls.maximum_species_fraction)
        || controls.maximum_species_fraction >= 1.0
        || !positive(controls.mass_tolerance_kg)
        || !positive(controls.energy_tolerance_j)
    {
        return Err(Error::Invalid(
            "bounded graph, exact endpoint map and positive physical limits required",
        ));
    }
    for (&id, solve) in endpoints {
        check(&mut cancelled)?;
        solve.validate().map_err(|e| material(id, e))?;
    }
    let mut rows = Vec::with_capacity(transfers.len());
    for row in transfers {
        check(&mut cancelled)?;
        if row.donor == row.receiver
            || !parcels.contains_key(&row.donor)
            || !parcels.contains_key(&row.receiver)
            || !row.mass_kg.is_finite()
            || row.mass_kg < 0.0
        {
            return Err(Error::Invalid(
                "finite nonnegative species mass and distinct existing owners required",
            ));
        }
        if !parcels[&row.donor]
            .model()
            .thermodynamics()
            .components()
            .iter()
            .any(|s| s.species() == &row.species)
        {
            return Err(Error::Invalid(
                "selective edge names a species absent from its original donor",
            ));
        }
        if row.mass_kg > 0.0 {
            rows.push(row);
        }
    }
    let mut heat_rows = Vec::with_capacity(heat.len());
    for row in heat {
        check(&mut cancelled)?;
        let [left, right] = row.parcels;
        if left == right
            || !parcels.contains_key(&left)
            || !parcels.contains_key(&right)
            || !row.heat_left_to_right_j.is_finite()
        {
            return Err(Error::Invalid(
                "additional heat needs finite energy and distinct existing owners",
            ));
        }
        if row.heat_left_to_right_j != 0.0 {
            heat_rows.push(if left < right {
                ([left, right], row.heat_left_to_right_j)
            } else {
                ([right, left], -row.heat_left_to_right_j)
            });
        }
    }
    heat_rows.sort_by(|(a, aq), (b, bq)| a.cmp(b).then_with(|| aq.total_cmp(bq)));
    let mut combined_heat: BTreeMap<[u64; 2], Sum> = BTreeMap::new();
    for (pair, value) in heat_rows {
        check(&mut cancelled)?;
        combined_heat
            .entry(pair)
            .or_default()
            .add(value)
            .map_err(|e| material(pair[0], e))?;
    }
    rows.sort_by(|a, b| {
        (a.donor, a.receiver, &a.species)
            .cmp(&(b.donor, b.receiver, &b.species))
            .then_with(|| a.mass_kg.total_cmp(&b.mass_kg))
    });
    let mut combined: BTreeMap<(u64, u64, SpeciesId), Sum> = BTreeMap::new();
    for row in rows {
        check(&mut cancelled)?;
        combined
            .entry((row.donor, row.receiver, row.species.clone()))
            .or_default()
            .add(row.mass_kg)
            .map_err(|e| material(row.donor, e))?;
    }
    let mut budgets: BTreeMap<u64, Budget> = BTreeMap::new();
    let mut fluxes = Vec::with_capacity(combined.len());
    for ((donor, receiver, species), amount) in combined {
        check(&mut cancelled)?;
        let dm = amount.value;
        let h = source_enthalpy(&parcels[&donor], &species).map_err(|e| material(donor, e))?;
        let carried = dm * h;
        if !carried.is_finite() {
            return Err(Error::Invalid(
                "original gross enthalpy report is unrepresentable",
            ));
        }
        let original_source = parcels[&donor]
            .model()
            .thermodynamics()
            .components()
            .iter()
            .find(|source| source.species() == &species)
            .ok_or_else(|| {
                material(
                    donor,
                    GasEnergyError::Invalid(
                        "original donor species is absent from its retained catalog",
                    ),
                )
            })?;
        let sensible_h = original_source
            .molar_enthalpy_change(
                Temperature::new(parcels[&receiver].state().state().temperature),
                Temperature::new(parcels[&donor].state().state().temperature),
            )
            .map_err(|e| material(receiver, source_error(e)))?
            .value()
            / original_source.molar_mass().value();
        for (id, sign) in [(donor, -1.0), (receiver, 1.0)] {
            let map = |e| material(id, e);
            let budget = budgets.entry(id).or_default();
            if id == donor {
                budget.mass_out.add(dm).map_err(map)?;
                budget
                    .withdrawn
                    .entry(species.clone())
                    .or_default()
                    .add(dm)
                    .map_err(map)?;
            } else {
                budget.incoming.insert(donor);
                budget.mass_in.add(dm).map_err(map)?;
                budget.sensible_enthalpy.add(dm * sensible_h).map_err(map)?;
                budget.relative_energy.add(dm * sensible_h).map_err(map)?;
            }
            budget
                .species
                .entry(species.clone())
                .or_default()
                .add(sign * dm)
                .map_err(map)?;
        }
        fluxes.push(SpeciesNetworkFlux {
            transfer: SpeciesNetworkTransfer {
                donor,
                receiver,
                species,
                mass_kg: dm,
            },
            source_enthalpy_j_kg: h,
            carried_enthalpy_j: carried,
        });
    }
    let mut heat_fluxes = Vec::with_capacity(combined_heat.len());
    for (pair, amount) in combined_heat {
        check(&mut cancelled)?;
        let value = amount.value;
        if value == 0.0 {
            continue;
        }
        for (id, sign) in [(pair[0], -1.0), (pair[1], 1.0)] {
            let budget = budgets.entry(id).or_default();
            budget.heat.add(sign * value).map_err(|e| material(id, e))?;
            budget
                .relative_energy
                .add(sign * value)
                .map_err(|e| material(id, e))?;
        }
        heat_fluxes.push(SpeciesNetworkHeatTransfer {
            parcels: pair,
            heat_left_to_right_j: value,
        });
    }
    let external = boundary::accumulate(parcels, boundaries, &mut budgets, controls, &mut cancelled)?;
    // Check aggregate donor allowances before the first nonlinear material solve.
    let mut originals = BTreeMap::new();
    for (&id, budget) in &budgets {
        check(&mut cancelled)?;
        let before =
            inventory(parcels[&id].model(), parcels[&id].mass_kg()).map_err(|e| material(id, e))?;
        for (species, withdrawn) in &budget.withdrawn {
            let fraction = withdrawn.value / before.get(species).copied().unwrap_or(0.0);
            if !fraction.is_finite() || fraction > controls.maximum_species_fraction {
                return Err(Error::Withdrawal {
                    id,
                    species: species.clone(),
                    fraction,
                    maximum: controls.maximum_species_fraction,
                });
            }
        }
        originals.insert(id, before);
    }
    let mut next = BTreeMap::new();
    let mut accounts = Vec::with_capacity(parcels.len());
    let (mut total_mass, mut total_energy) = (Sum::default(), Sum::default());
    let mut global_species: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    for (&id, old) in parcels {
        check(&mut cancelled)?;
        let Some(budget) = budgets.get(&id) else {
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
        let before = &originals[&id];
        let mut target = before.clone();
        for (species, change) in &budget.species {
            let value = before.get(species).copied().unwrap_or(0.0);
            let updated = value + change.value;
            // A net increment below half an ulp of the inventory rounds away. That is
            // refused only when it exceeds the declared mass budget; smaller residues
            // stay bounded by the global species balance below.
            if !positive(updated)
                || (updated == value && change.value.abs() > controls.mass_tolerance_kg)
            {
                return Err(Error::Invalid(
                    "selective network exhausted or lost a species increment",
                ));
            }
            target.insert(species.clone(), updated);
        }
        // Join source catalogs, not physical temperatures or transport values.
        let mut portions = vec![(old.model(), old.mass_kg())];
        let mut catalog_mass = Sum::default();
        catalog_mass.add(old.mass_kg()).map_err(map)?;
        for donor in &budget.incoming {
            check(&mut cancelled)?;
            let original = &parcels[donor];
            portions.push((original.model(), original.mass_kg()));
            catalog_mass.add(original.mass_kg()).map_err(map)?;
        }
        if let Some(catalogs) = external.catalogs.get(&id) {
            for &catalog in catalogs {
                check(&mut cancelled)?;
                portions.push((catalog, old.mass_kg()));
                catalog_mass.add(old.mass_kg()).map_err(map)?;
            }
        }
        let catalog = mixed_model_from_portions(&portions, catalog_mass.value).map_err(map)?;
        let delta_mass = (budget.mass_in.value - budget.mass_out.value) + budget.boundary_delta_mass.value;
        let mass = old.mass_kg() + delta_mass;
        if !positive(mass) {
            return Err(Error::Invalid(
                "selective network final mass is unrepresentable",
            ));
        }
        let t0 = old.state().state().temperature;
        let mut thermal_target = Sum::default();
        thermal_target
            .add(budget.relative_energy.value)
            .map_err(map)?;
        let mut net_enthalpy = Sum::default();
        net_enthalpy
            .add(budget.sensible_enthalpy.value)
            .map_err(map)?;
        let mut reference_change = Sum::default();
        for (species, change) in &budget.species {
            check(&mut cancelled)?;
            if change.value == 0.0 {
                continue;
            }
            let source = catalog
                .thermodynamics()
                .components()
                .iter()
                .find(|source| source.species() == species)
                .ok_or_else(|| {
                    material(
                        id,
                        GasEnergyError::Invalid(
                            "target catalog lost an original transported species",
                        ),
                    )
                })?;
            let molecular_mass = source.molar_mass().value();
            let evaluated = source
                .evaluate(Temperature::new(t0))
                .map_err(|e| material(id, source_error(e)))?;
            let properties = evaluated.properties();
            // Σ δm_i R_i T0 is evaluated only after gross counterstreams cancel.
            // The remaining stream term is Σ incoming dm_i Δh_i(T0,Tdonor).
            thermal_target
                .add(change.value * (R * t0 / molecular_mass))
                .map_err(map)?;
            net_enthalpy
                .add(change.value * (properties.h().value() / molecular_mass))
                .map_err(map)?;
            reference_change
                .add(change.value * (properties.u().value() / molecular_mass))
                .map_err(map)?;
        }
        let (new, iterations) =
            if target == *before && delta_mass == 0.0 && thermal_target.value == 0.0 {
                (old.clone(), 0)
            } else {
                let model = if target == *before && delta_mass == 0.0 {
                    old.model().clone()
                } else {
                    model_for_inventory(&catalog, &target).map_err(map)?
                };
                relative_endpoint(
                    old,
                    model,
                    mass,
                    thermal_target.value,
                    endpoints[&id],
                    &mut cancelled,
                )
                .map_err(map)?
            };
        let thermal_change = new.mass_kg()
            * new
                .model()
                .thermodynamics()
                .specific_internal_energy_change(
                    Temperature::new(t0),
                    Temperature::new(new.state().state().temperature),
                )
                .map_err(|e| material(id, source_error(e)))?
                .value();
        let mut stored = reference_change;
        stored.add(thermal_change).map_err(map)?;
        let du = stored.value;
        let dm = new.mass_kg() - old.mass_kg();
        // This is ΔU−Σ dm*h−Q with all conserved formation constants already
        // canceled, not subtraction of the potentially enormous displayed totals.
        let er = thermal_change - thermal_target.value;
        let mr = dm - delta_mass;
        if ![du, dm, er, mr].iter().all(|v| v.is_finite())
            || er.abs() > controls.energy_tolerance_j
            || mr.abs() > controls.mass_tolerance_kg
        {
            return Err(Error::Invalid(
                "selective network local mass or energy balance failed",
            ));
        }
        let actual = inventory(new.model(), new.mass_kg()).map_err(map)?;
        let mut maximum: f64 = 0.0;
        for (species, expected) in &target {
            let value = actual.get(species).copied().unwrap_or(0.0);
            let intended = budget.species.get(species).map_or(0.0, |s| s.value);
            let representation = (value - before.get(species).copied().unwrap_or(0.0)) - intended;
            maximum = maximum
                .max((value - expected).abs())
                .max(representation.abs());
            global_species
                .entry(species.clone())
                .or_default()
                .add(value - before.get(species).copied().unwrap_or(0.0))
                .map_err(map)?;
        }
        if maximum > controls.mass_tolerance_kg {
            return Err(Error::Invalid(
                "selective network local species balance failed",
            ));
        }
        total_mass.add(dm).map_err(map)?;
        total_energy.add(er).map_err(map)?;
        let incoming_mass_kg = budget.mass_in.value + budget.boundary_mass_in.value;
        let outgoing_mass_kg = budget.mass_out.value + budget.boundary_mass_out.value;
        if !incoming_mass_kg.is_finite() || !outgoing_mass_kg.is_finite() {
            return Err(Error::Invalid("combined gross gas mass reports overflowed"));
        }
        accounts.push(SpeciesNetworkParcelStep {
            id,
            incoming_mass_kg,
            outgoing_mass_kg,
            net_enthalpy_j: net_enthalpy.value,
            net_additional_heat_j: budget.heat.value,
            energy_change_j: du,
            energy_residual_j: er,
            mass_residual_kg: mr,
            maximum_species_mass_residual_kg: maximum,
            iterations,
        });
        next.insert(id, new);
    }
    total_mass.add(-external.mass.value).map_err(|e| material(0, e))?;
    for (species, amount) in &external.species {
        global_species.entry(species.clone()).or_default()
            .add(-amount.value).map_err(|e| material(0, e))?;
    }
    let maximum = global_species
        .values()
        .map(|s| s.value.abs())
        .fold(0.0_f64, f64::max);
    if maximum > controls.mass_tolerance_kg
        || total_mass.value.abs() > controls.mass_tolerance_kg
        || total_energy.value.abs() > controls.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "selective network global species, mass or energy balance failed",
        ));
    }
    check(&mut cancelled)?;
    Ok(SpeciesNetworkCandidate {
        parcels: next,
        step: SpeciesNetworkStep {
            fluxes,
            heat: heat_fluxes,
            boundaries: external.fluxes,
            parcels: accounts,
            total_mass_residual_kg: total_mass.value,
            maximum_species_mass_residual_kg: maximum,
            total_energy_residual_j: total_energy.value,
        },
    })
}

#[cfg(test)]
mod heat_tests;
