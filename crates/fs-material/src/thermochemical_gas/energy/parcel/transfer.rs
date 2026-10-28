//! Conservative, species-resolved transfer between two quiescent gas volumes.
//!
//! A prescribed mass carries the donor's upstream specific ENTHALPY, including
//! flow work, not its internal energy. Both endpoints use the existing NASA-9
//! energy inversion, EOS and transport owners. The receiver combines species
//! inventories, never Kelvin temperatures or unrelated material coefficients.
//! Shared species must retain identical caloric and transport source models.
//!
//! This is a first-order, upstream-state finite-volume exchange at fixed volumes.
//! It supplies no flow-rate law, pressure-direction rule, mixing time, reaction,
//! moving-boundary step or entropy-production guarantee. Jet kinetic energy is
//! assumed thermalized in the receiving well-mixed volume; no resolved momentum
//! is stored here. A rate/port owner must bound and refine its physical timestep.

/// Source-backed converging-nozzle flow and pressure-driven paired exchange.
pub mod nozzle;
/// Shared-snapshot multi-inlet/outlet transfer with one reconstruction per volume.
pub mod snapshot;
/// Selective counterflow with individual species enthalpy and inventory conservation.
pub mod species;

use std::collections::BTreeMap;
use super::ClosedGasParcel;
use super::super::{GasEnergyError, GasEnergySolve, poll, positive, updated_point};
use crate::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_thermochem::{Composition, CompositionBasis, Nasa9StandardStateModelV1, SpeciesId};
use fs_thermochem::mixture::{FrozenIdealGasMixtureModelV1, MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1};

/// Explicit endpoint solve budgets and absolute extensive-balance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasMassTransferSolve {
    /// Donor temperature bracket and internal-energy residual tolerance.
    pub donor: GasEnergySolve,
    /// Receiver temperature bracket and internal-energy residual tolerance.
    pub receiver: GasEnergySolve,
    /// Maximum absolute total AND per-species transfer discrepancy [kg].
    pub mass_tolerance_kg: f64,
    /// Maximum absolute per-volume AND pair internal-energy discrepancy [J].
    pub energy_tolerance_j: f64,
}
impl GasMassTransferSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        self.donor.validate()?;
        self.receiver.validate()?;
        if !positive(self.mass_tolerance_kg) || !positive(self.energy_tolerance_j) {
            return Err(invalid("positive finite transfer mass and energy tolerances required"));
        }
        Ok(())
    }
}

/// Original canonical species and its upstream transferred mass.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSpeciesTransfer {
    /// Existing species identity; no material-name inference.
    pub species: SpeciesId,
    /// Upstream species mass delivered from donor to receiver [kg].
    pub mass_kg: f64,
}

/// Actual extensive balances of a complete two-volume transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct GasMassTransferStep {
    /// Prescribed mass out of donor and into receiver [kg].
    pub transferred_mass_kg: f64,
    /// Upstream donor specific enthalpy [J/kg], with original reference retained.
    pub source_enthalpy_j_kg: f64,
    /// Signed energy carried into the receiver, dm*h_upstream [J].
    /// May be negative under a valid formation-energy convention.
    pub transferred_enthalpy_j: f64,
    /// Donor species-reference change plus original-law sensible storage [J].
    /// Actual inventory rounding is separately bounded by the mass budgets.
    pub donor_energy_change_j: f64,
    /// Receiver species-reference change plus original-law sensible storage [J].
    pub receiver_energy_change_j: f64,
    /// Sum of local energy residuals after conserved species datums cancel [J].
    /// Do not reconstruct this small residual by subtracting rounded reports.
    pub energy_residual_j: f64,
    /// Actual net mass change across the two volumes [kg].
    pub mass_residual_kg: f64,
    /// Largest donor/receiver discrepancy from an individual species flux [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Actual transferred species inventory in canonical order.
    pub species: Vec<GasSpeciesTransfer>,
    /// Sum of interior iterations in the two existing energy inversions.
    pub iterations: usize,
}

/// Tentative endpoints. Publish BOTH only after the surrounding coupling accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasMassTransferCandidate {
    donor: ClosedGasParcel,
    receiver: ClosedGasParcel,
    step: GasMassTransferStep,
}
impl GasMassTransferCandidate {
    /// Donor after withdrawal, with its original source laws and volume.
    #[must_use]
    pub const fn donor(&self) -> &ClosedGasParcel { &self.donor }
    /// Receiver with its new composition, density, energy and complete medium.
    #[must_use]
    pub const fn receiver(&self) -> &ClosedGasParcel { &self.receiver }
    /// Actual transfer and conservation checks.
    #[must_use]
    pub const fn step(&self) -> &GasMassTransferStep { &self.step }
    /// Move both endpoints and the accounting without another material query.
    #[must_use]
    pub fn into_parts(self) -> (ClosedGasParcel, ClosedGasParcel, GasMassTransferStep) {
        (self.donor, self.receiver, self.step)
    }
}
fn invalid(message: &'static str) -> GasEnergyError { GasEnergyError::Invalid(message) }
fn source_error(error: impl core::fmt::Display) -> GasEnergyError {
    crate::thermochemical_gas::failure("mass transfer", error).into()
}

fn mass_fractions(model: &ThermochemicalGasModel) -> Result<Composition, GasEnergyError> {
    let masses: Vec<_> = model.thermodynamics().components().iter().map(|s| s.molar_mass()).collect();
    model.thermodynamics().mole_composition().to_mass_fractions(&masses).map_err(source_error)
}
fn inventory(model: &ThermochemicalGasModel, mass: f64) -> Result<BTreeMap<SpeciesId, f64>, GasEnergyError> {
    model.thermodynamics().components().iter().zip(mass_fractions(model)?.fractions())
        .map(|(species, fraction)| {
            let value = mass * fraction;
            if !positive(value) { return Err(invalid("positive species mass was lost or overflowed")); }
            Ok((species.species().clone(), value))
        }).collect()
}
fn water(model: &ThermochemicalGasModel) -> Option<SpeciesId> {
    model.water_component.map(|i| model.thermodynamics().components()[i].species().clone())
}

// Each portion is an extensive mass with its ORIGINAL source model. Both the
// pair and multiport paths use this single species/transport/pack join.
fn mixed_model_from_portions(
    portions: &[(&ThermochemicalGasModel, f64)], total_mass: f64,
) -> Result<ThermochemicalGasModel, GasEnergyError> {
    if portions.is_empty() || !positive(total_mass)
        || portions.iter().any(|(_, amount)| !positive(*amount)) {
        return Err(invalid("positive finite material portions and total mass required"));
    }
    let mut components: BTreeMap<SpeciesId, (Nasa9StandardStateModelV1, SpeciesTransportMaterial, f64)> = BTreeMap::new();
    for &(model, amount) in portions {
        let fractions = mass_fractions(model)?;
        for ((species, transport), fraction) in model.thermodynamics().components().iter()
            .zip(model.transport_materials()).zip(fractions.fractions()) {
            let added = amount * fraction;
            if !positive(added) { return Err(invalid("unrepresentable species transfer")); }
            if let Some((existing, data, weight)) = components.get_mut(species.species()) {
                if &*existing != species || &*data != transport {
                    return Err(invalid("shared gas species has conflicting thermochemical or transport source laws"));
                }
                *weight += added;
                if !positive(*weight) { return Err(invalid("mixed species mass overflow")); }
            } else {
                if components.len() == MAX_FROZEN_IDEAL_GAS_COMPONENTS_V1 {
                    return Err(invalid("mixed gas exceeds the existing species work bound"));
                }
                components.insert(species.species().clone(), (species.clone(), transport.clone(), added));
            }
        }
    }
    let mut explicit_water = None;
    for &(source, _) in portions {
        if let Some(id) = water(source) {
            if explicit_water.as_ref().is_some_and(|previous| previous != &id) {
                return Err(invalid("gas sources disagree about the explicit water species"));
            }
            explicit_water = Some(id);
        }
    }
    let mut species = Vec::with_capacity(components.len());
    let mut transport = Vec::with_capacity(components.len());
    let mut fractions = Vec::with_capacity(components.len());
    for (_, (model, data, amount)) in components {
        species.push(model); transport.push(data); fractions.push(amount / total_mass);
    }
    let masses: Vec<_> = species.iter().map(|s| s.molar_mass()).collect();
    // The shared quantity owner performs mass/mole conversion. Its roundoff-
    // tolerant carrier and the thermochemical owner's EXACT canonical sum have
    // different contracts. Reconcile only the final rounding residue; never
    // drop a positive trace species. Actual per-species mass is checked below.
    let composition = Composition::new(CompositionBasis::MassFraction, fractions)
        .and_then(|c| c.to_mole_fractions(&masses)).map_err(source_error)?;
    let mut mole = composition.fractions().to_vec();
    let sum: f64 = mole.iter().sum();
    if (sum - 1.0).abs() > 256.0 * f64::EPSILON {
        return Err(invalid("mixed composition conversion exceeds roundoff"));
    }
    if sum != 1.0 {
        let last = mole.len() - 1;
        let prefix: f64 = mole[..last].iter().sum();
        mole[last] = 1.0 - prefix;
    }
    if mole.iter().any(|x| !positive(*x)) {
        return Err(invalid("a mixed trace species is not representable"));
    }
    let composition = Composition::new(CompositionBasis::MoleFraction, mole).map_err(source_error)?;
    let thermo = FrozenIdealGasMixtureModelV1::new(species, composition).map_err(source_error)?;
    let mut model = ThermochemicalGasModel::new(thermo, transport, explicit_water.as_ref())?;
    let mut receipts = BTreeMap::new();
    for &(source, _) in portions {
        for receipt in source.source_packs() {
            if let Some(previous) = receipts.insert(receipt.species.clone(), receipt.clone()) {
                if previous != *receipt {
                    return Err(invalid("shared gas species has conflicting original pack bindings"));
                }
            }
        }
    }
    model.source_packs = receipts.into_values().collect();
    Ok(model)
}

fn endpoint(
    old: &ClosedGasParcel, model: ThermochemicalGasModel, mass: f64, specific_change: f64,
    controls: GasEnergySolve, cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, usize), GasEnergyError> {
    let solved = model.solve_energy_change(&old.parcel.point, mass / old.volume_m3(),
        old.state().specific_internal_energy_j_kg(), specific_change, controls, cancelled)?;
    let mut parcel = old.parcel.clone();
    parcel.point = updated_point(&old.parcel.point, solved.state().state().temperature, solved.state().state().pressure)?;
    parcel.state = solved.state().clone();
    parcel.model = model;
    parcel.mass_kg = mass;
    Ok((ClosedGasParcel::from_isochoric(parcel), solved.iterations()))
}
fn stored_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    // Avoid subtracting two large total formation energies when possible.
    next.mass_kg() * (next.state().specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg())
        + (next.mass_kg() - old.mass_kg()) * old.state().specific_internal_energy_j_kg()
}

impl ClosedGasParcel {
    /// Transfer a prescribed upstream mass to another fixed gas volume.
    /// Existing accepted parcels are never mutated; the two proposed endpoints
    /// must be published together. The returned receiver remains compatible
    /// with existing heat, piston and acoustic consumers.
    ///
    /// dm*h_d leaves the donor and enters the receiver exactly once. Each
    /// species crosses in the upstream mass fraction. Caloric reference and
    /// transport conflicts refuse; unsupported new mixtures cannot silently
    /// fall back to air. Other source axes stay local to each volume.
    ///
    /// This explicit upstream step can cool the donor below its admitted
    /// temperature for a large withdrawal. Reduce the physical step instead
    /// of clipping energy or inventing a gas phase. No vacuum endpoint is
    /// admitted. Zero transfer preserves both endpoints bit-for-bit.
    pub fn transfer_mass_to(
        &self, receiver: &Self, mass_kg: f64, controls: GasMassTransferSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMassTransferCandidate, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !mass_kg.is_finite() || mass_kg < 0.0 || mass_kg >= self.mass_kg() {
            return Err(invalid("finite nonnegative transfer smaller than donor inventory required"));
        }
        let source_h = self.state().specific_enthalpy_j_kg();
        if mass_kg == 0.0 {
            return Ok(GasMassTransferCandidate { donor: self.clone(), receiver: receiver.clone(),
                step: GasMassTransferStep { transferred_mass_kg: 0.0, source_enthalpy_j_kg: source_h,
                    transferred_enthalpy_j: 0.0, donor_energy_change_j: 0.0, receiver_energy_change_j: 0.0,
                    energy_residual_j: 0.0, mass_residual_kg: 0.0, maximum_species_mass_residual_kg: 0.0,
                    species: vec![], iterations: 0 } });
        }
        let donor_mass = self.mass_kg() - mass_kg;
        let receiver_mass = receiver.mass_kg() + mass_kg;
        if mass_kg <= controls.mass_tolerance_kg || !positive(donor_mass) || !positive(receiver_mass)
            || donor_mass == self.mass_kg() || receiver_mass == receiver.mass_kg()
            || ((self.mass_kg()-donor_mass)-mass_kg).abs() > controls.mass_tolerance_kg
            || ((receiver_mass-receiver.mass_kg())-mass_kg).abs() > controls.mass_tolerance_kg {
            return Err(invalid("transfer is not representable within the absolute mass tolerance"));
        }
        if std::ptr::eq(self, receiver) {
            return Err(invalid("nonzero mass transfer needs two distinct physical gas owners"));
        }
        // Bulk pairs and multiport flows share the species-relative caloric
        // owner. In particular h-u is R*T, even when an arbitrary formation
        // datum makes subtracting the absolute values lose all flow work.
        let candidate = snapshot::transfer_candidate(
            &BTreeMap::from([(0, self.clone()), (1, receiver.clone())]),
            &[snapshot::GasSnapshotTransfer { donor: 0, receiver: 1, mass_kg }],
            &BTreeMap::from([(0, controls.donor), (1, controls.receiver)]),
            snapshot::GasSnapshotSolve {
                maximum_withdrawn_fraction: 1.0_f64.next_down(),
                mass_tolerance_kg: controls.mass_tolerance_kg,
                energy_tolerance_j: controls.energy_tolerance_j,
            },
            &mut cancelled,
        ).map_err(|error| match error {
            snapshot::GasSnapshotError::Cancelled => GasEnergyError::Cancelled,
            snapshot::GasSnapshotError::Parcel { error, .. } => error,
            snapshot::GasSnapshotError::Invalid(message) => invalid(message),
            snapshot::GasSnapshotError::Withdrawal { .. } => invalid("transfer exhausts its original donor"),
        })?;
        let (mut endpoints, report) = candidate.into_parts();
        let donor_account = &report.parcels[0];
        let receiver_account = &report.parcels[1];
        let maximum = report.maximum_species_mass_residual_kg
            .max(donor_account.maximum_species_mass_residual_kg)
            .max(receiver_account.maximum_species_mass_residual_kg);
        let flux = report.fluxes.into_iter().next().expect("one positive bulk transfer");
        let step = GasMassTransferStep { transferred_mass_kg: mass_kg, source_enthalpy_j_kg: source_h,
            transferred_enthalpy_j: flux.enthalpy_j,
            donor_energy_change_j: donor_account.energy_change_j,
            receiver_energy_change_j: receiver_account.energy_change_j,
            energy_residual_j: report.total_energy_residual_j,
            mass_residual_kg: report.total_mass_residual_kg,
            maximum_species_mass_residual_kg: maximum,
            species: flux.species,
            iterations: donor_account.iterations + receiver_account.iterations };
        poll(&mut cancelled)?;
        Ok(GasMassTransferCandidate {
            donor: endpoints.remove(&0).expect("original donor endpoint"),
            receiver: endpoints.remove(&1).expect("original receiver endpoint"), step })
    }
}

#[cfg(test)]
mod tests;
