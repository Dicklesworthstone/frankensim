//! Steady source-backed solution-diffusion through a dense gas membrane.
//!
//! For dilute Henry sorption and permeability constant at fixed temperature,
//! J_i = (P_i / thickness) (x_i,left p_left - x_i,right p_right).
//! P [mol/(m s Pa)], permeance [mol/(m2 s Pa)], and intrinsic porous k [m2]
//! are distinct quantities. No mass or STP-volume conversion is inferred.
//! SI basis: https://polymer.nims.go.jp/PoLyInfo/guide/en/property.html .
//!
//! This is a steady negligible-storage film, with explicitly authored species
//! sources and mixture/humidity context. No competitive sorption, swelling,
//! plasticization, interface resistance, nonideal fugacity, thermal
//! transpiration or transient breakthrough is modeled. Pure-gas measurements
//! do not automatically qualify mixed-gas operation.

use std::collections::BTreeMap;
use super::super::SpeciesNetworkTransfer;
use super::super::super::{catalog, positive, source_error};
use crate::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use crate::thermochemical_gas::energy::{GasEnergyError, updated_point};
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use crate::thermochemical_gas::energy::parcel::transfer::snapshot::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Dims;
use fs_thermochem::SpeciesId;

pub const MOLAR_PERMEABILITY_DIMS: Dims = Dims([0, -1, 1, 0, 0, 1]);
pub const MOLAR_PERMEANCE_DIMS: Dims = Dims([-1, -1, 1, 0, 0, 1]);
pub const MOLAR_PERMEABILITY_PROPERTY: &str = "gas_molar_permeability";
pub const MAX_MEMBRANE_SPECIES: usize = 32;

/// Dense-film diffusivity times molar Henry solubility, in mol/(m s Pa).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolarPermeability(f64);
impl MolarPermeability {
    /// Zero is an explicitly impermeable species.
    pub fn try_new(value_mol_m_s_pa: f64) -> Result<Self, GasEnergyError> {
        if !value_mol_m_s_pa.is_finite() || value_mol_m_s_pa < 0.0 {
            return Err(GasEnergyError::Invalid("molar permeability must be finite and nonnegative"));
        }
        Ok(Self(value_mol_m_s_pa))
    }
    pub const fn value_si(self) -> f64 { self.0 }
    /// Thickness conversion is explicit; intrinsic porous k cannot enter here.
    pub fn permeance(self, thickness_m: f64) -> Result<MolarPermeance, GasEnergyError> {
        if !positive(thickness_m) {
            return Err(GasEnergyError::Invalid("membrane thickness must be positive and finite"));
        }
        let value = self.0 / thickness_m;
        if !value.is_finite() || (self.0 > 0.0 && value == 0.0) {
            return Err(GasEnergyError::Invalid("membrane permeance is unrepresentable"));
        }
        Ok(MolarPermeance(value))
    }
}
/// Pressure-normalized molar flux, in mol/(m2 s Pa).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolarPermeance(f64);
impl MolarPermeance {
    pub const fn value_si(self) -> f64 { self.0 }
}

#[derive(Clone, Debug)]
pub struct MembranePermeabilityMaterial<'a> {
    pub species: SpeciesId,
    pub material: &'a MaterialCard,
    /// Other authored process, mixture and humidity conditions stay fixed.
    pub context: &'a QueryPoint,
    pub selection: MaterialPropertySelection,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneGeometry {
    /// Zero closes the path without querying inactive source coefficients.
    pub area_m2: f64,
    pub thickness_m: f64,
}
#[derive(Clone, Debug)]
pub struct MembraneNetworkLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub geometry: MembraneGeometry,
    /// One explicit permeability, including zeros, per catalog species.
    pub materials: &'a [MembranePermeabilityMaterial<'a>],
    /// Physical approximation limit, not an error estimate; zero requires equal T.
    pub maximum_temperature_difference_k: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneSpeciesFlux {
    pub species: SpeciesId,
    pub partial_pressures_pa: [f64; 2],
    pub permeability: MolarPermeability,
    pub permeance: MolarPermeance,
    pub molar_flux_mol_m2_s: f64,
    pub mass_rate_left_to_right_kg_s: f64,
    pub material: ResolvedMaterialStatePoint,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneFlux {
    pub face_temperature_k: f64,
    pub face_pressure_pa: f64,
    pub species: Vec<MembraneSpeciesFlux>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneNetworkLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub flux: MembraneFlux,
}
#[derive(Clone, Debug, PartialEq)]
pub enum MembraneNetworkError {
    Invalid(&'static str),
    Face { id: u64, error: GasEnergyError },
    EquilibriumCrossed { id: u64, species: SpeciesId },
    Cancelled,
}
impl core::fmt::Display for MembraneNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material membrane transport: {self:?}")
    }
}
impl std::error::Error for MembraneNetworkError {}

fn partial_pressure(parcel: &ClosedGasParcel, species: &SpeciesId) -> Result<f64, GasEnergyError> {
    let Some((_, &fraction)) = parcel.model().thermodynamics().components().iter()
        .zip(parcel.model().thermodynamics().mole_composition().fractions())
        .find(|(s, _)| s.species() == species) else { return Ok(0.0); };
    let pressure = fraction * parcel.state().state().pressure;
    if !pressure.is_finite() || pressure < 0.0 || (fraction > 0.0 && pressure == 0.0) {
        return Err(GasEnergyError::Invalid("membrane species partial pressure is unrepresentable"));
    }
    Ok(pressure)
}

impl ClosedGasParcel {
    /// Sample a steady film without changing gas inventories.
    /// Opposite species fluxes coexist; no zero-sum molar/mass condition applies.
    /// The shared exchange owner later carries original partial enthalpies once.
    pub fn membrane_flux_to(
        &self, right: &Self, geometry: MembraneGeometry,
        materials: &[MembranePermeabilityMaterial<'_>],
        maximum_temperature_difference_k: f64, mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneFlux, GasEnergyError> {
        let invalid = GasEnergyError::Invalid;
        let mut check = || if cancelled() { Err(GasEnergyError::Cancelled) } else { Ok(()) };
        check()?;
        if !geometry.area_m2.is_finite() || geometry.area_m2 < 0.0
            || !positive(geometry.thickness_m) || materials.len() > MAX_MEMBRANE_SPECIES
            || !maximum_temperature_difference_k.is_finite() || maximum_temperature_difference_k < 0.0 {
            return Err(invalid("bounded membrane materials, nonnegative area/temperature difference and positive thickness required"));
        }
        let (a, b) = (self.state().state(), right.state().state());
        let mut flux = MembraneFlux {
            face_temperature_k: 0.5 * a.temperature + 0.5 * b.temperature,
            face_pressure_pa: 0.5 * a.pressure + 0.5 * b.pressure, species: Vec::new(),
        };
        if geometry.area_m2 == 0.0 { check()?; return Ok(flux); }
        if (a.temperature - b.temperature).abs() > maximum_temperature_difference_k {
            return Err(invalid("gas states exceed the declared isothermal membrane support"));
        }
        let joined = catalog(self, right)?;
        let components = joined.thermodynamics().components();
        if components.len() > MAX_MEMBRANE_SPECIES || materials.len() != components.len() {
            return Err(invalid("membrane needs one explicit permeability per bounded catalog species"));
        }
        let mut ordered = BTreeMap::new();
        for binding in materials {
            check()?;
            if ordered.insert(&binding.species, binding).is_some() {
                return Err(invalid("duplicate membrane species permeability"));
            }
        }
        let requirement = ScalarPropertyRequirement::try_new(
            MOLAR_PERMEABILITY_PROPERTY, MOLAR_PERMEABILITY_DIMS, ScalarAdmissibility::NonNegative,
        ).map_err(source_error)?;
        for species in components {
            check()?;
            let binding = ordered.get(species.species())
                .ok_or(invalid("missing or foreign membrane species permeability"))?;
            let point = updated_point(binding.context, flux.face_temperature_k, flux.face_pressure_pa)?;
            let material = resolve_material_state_point(binding.material, &point,
                std::slice::from_ref(&requirement), binding.selection.clone()).map_err(source_error)?;
            let property = material.property(MOLAR_PERMEABILITY_PROPERTY)
                .ok_or(invalid("resolved membrane permeability is missing"))?;
            // A sampled coefficient cannot freeze an arbitrary sorption law.
            // Only scalar/canonical-temperature laws admit this Henry rung.
            if !property.is_constant_at_fixed_temperature()
                || property.source_axes().iter().any(|(name, _)| name != "temperature")
                || !property.is_static_at_fixed_external_conditions()
                || property.answer().evidence.model.validity.bounds().keys().any(|name| {
                    name == "T" || property.answer().receipt.axis_quantities.get(name).is_some_and(|q|
                        (q.dims() == fs_qty::Temperature::DIMS && name != "temperature")
                        || (q.dims() == fs_qty::Pressure::DIMS && name != "pressure"))
                }) {
                return Err(invalid("linear solution-diffusion needs a constant-at-fixed-temperature source with canonical T/p axes"));
            }
            let permeability = MolarPermeability::try_new(property.value_si())?;
            let permeance = permeability.permeance(geometry.thickness_m)?;
            // Refresh the exact selected claim, never reselect under pressure.
            // Mean T defines the caller-admitted near-isothermal film law.
            for gas in [a, b] {
                check()?;
                let boundary = updated_point(binding.context, gas.temperature, gas.pressure)?;
                material.at_query_point(binding.material, &boundary).map_err(source_error)?;
            }
            let partial_pressures_pa = [partial_pressure(self, species.species())?,
                partial_pressure(right, species.species())?];
            let difference = partial_pressures_pa[0] - partial_pressures_pa[1];
            let molar_flux_mol_m2_s = permeance.value_si() * difference;
            let rate = (molar_flux_mol_m2_s * species.molar_mass().value()) * geometry.area_m2;
            if !molar_flux_mol_m2_s.is_finite() || !rate.is_finite()
                || (permeance.value_si() != 0.0 && difference != 0.0 && rate == 0.0) {
                return Err(invalid("membrane species flux is unrepresentable"));
            }
            flux.species.push(MembraneSpeciesFlux { species: species.species().clone(),
                partial_pressures_pa, permeability, permeance, molar_flux_mol_m2_s,
                mass_rate_left_to_right_kg_s: rate, material });
        }
        flux.species.sort_by(|a, b| a.species.cmp(&b.species));
        check()?;
        Ok(flux)
    }
}

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), MembraneNetworkError> {
    if cancelled() { Err(MembraneNetworkError::Cancelled) } else { Ok(()) }
}
/// Sample all links from the same unchanged inventory.
pub(super) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>, links: &[MembraneNetworkLink<'_>],
    duration_s: f64, mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<SpeciesNetworkTransfer>, Vec<MembraneNetworkLinkStep>), MembraneNetworkError> {
    use MembraneNetworkError as Error;
    check(&mut cancelled)?;
    if !positive(duration_s) || parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS {
        return Err(Error::Invalid("bounded nonempty membrane graph and positive duration required"));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1] || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || ordered.insert(link.id, link).is_some() {
            return Err(Error::Invalid("unique membrane IDs and distinct existing gas owners required"));
        }
    }
    let mut transfers = Vec::new();
    let mut reports = Vec::with_capacity(links.len());
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let [left, right] = link.parcels;
        let flux = parcels[&left].membrane_flux_to(&parcels[&right], link.geometry,
            link.materials, link.maximum_temperature_difference_k, &mut cancelled)
            .map_err(|error| if error == GasEnergyError::Cancelled { Error::Cancelled }
                else { Error::Face { id, error } })?;
        for row in &flux.species {
            check(&mut cancelled)?;
            let mass = row.mass_rate_left_to_right_kg_s * duration_s;
            if !mass.is_finite() || (row.mass_rate_left_to_right_kg_s != 0.0 && mass == 0.0) {
                return Err(Error::Invalid("integrated membrane species mass is unrepresentable"));
            }
            if mass == 0.0 { continue; }
            if transfers.len() == MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid("membrane graph exceeds the selective-transfer bound"));
            }
            let (donor, receiver) = if mass > 0.0 { (left, right) } else { (right, left) };
            transfers.push(SpeciesNetworkTransfer { donor, receiver,
                species: row.species.clone(), mass_kg: mass.abs() });
        }
        reports.push(MembraneNetworkLinkStep { id, parcels: link.parcels, flux });
    }
    Ok((transfers, reports))
}
/// Check original species directions on the final combined state.
/// No rate is resampled or inventory spent a second time.
pub(super) fn admit_partial_pressures(
    parcels: &BTreeMap<u64, ClosedGasParcel>, reports: &[MembraneNetworkLinkStep],
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), MembraneNetworkError> {
    use MembraneNetworkError as Error;
    check(&mut cancelled)?;
    for report in reports {
        check(&mut cancelled)?;
        let [left, right] = report.parcels;
        let (Some(left), Some(right)) = (parcels.get(&left), parcels.get(&right)) else {
            return Err(Error::Invalid("membrane admission needs original parcel identities"));
        };
        for row in &report.flux.species {
            check(&mut cancelled)?;
            let map = |error| Error::Face { id: report.id, error };
            let difference = partial_pressure(left, &row.species).map_err(map)?
                - partial_pressure(right, &row.species).map_err(map)?;
            if (row.mass_rate_left_to_right_kg_s > 0.0 && difference < 0.0)
                || (row.mass_rate_left_to_right_kg_s < 0.0 && difference > 0.0) {
                return Err(Error::EquilibriumCrossed { id: report.id, species: row.species.clone() });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
