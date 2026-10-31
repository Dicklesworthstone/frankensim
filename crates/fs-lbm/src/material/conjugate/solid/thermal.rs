//! Fixed-mass source calorimetry for the original transient conjugate owner.

use super::{SolidConductivityLaw, SolidConductivitySample};
use crate::conjugate::SolidMaterial;
use fs_matdb::ClaimSelection;
use fs_material::state_point::{MaterialStatePointError, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_material::thermal::transport::{
    ThermalTransportError,
    storage::{HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage},
};
use fs_qty::Density;

/// Preparation refusals retain the original material and caloric diagnoses.
#[derive(Clone, Debug, PartialEq)]
pub enum SolidThermalError {
    /// Missing source role, incompatible context or invalid reference state.
    Invalid(&'static str),
    /// Original complete material-bundle failure.
    Material(MaterialStatePointError),
    /// Original source-capacity path or integration failure.
    Storage(ThermalTransportError),
}

impl core::fmt::Display for SolidThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "solid thermal material: {self:?}")
    }
}
impl std::error::Error for SolidThermalError {}
impl From<MaterialStatePointError> for SolidThermalError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
impl From<ThermalTransportError> for SolidThermalError {
    fn from(error: ThermalTransportError) -> Self {
        Self::Storage(error)
    }
}

/// One solid's original conductivity and exact single-phase caloric integral.
///
/// The prepared storage uses the fixed reference mass in ONE cubic metre.
/// A native voxel scales this integral by its reference volume, so changing
/// source density with temperature cannot create or remove specimen mass.
/// The caller selects a specific-heat law compatible with the held mechanical
/// and environmental conditions. No Cp/Cv conversion, latent plateau, phase
/// change, thermal expansion work or deformation is inferred.
#[derive(Clone, Debug)]
pub struct SolidThermalLaw<'a> {
    conductivity: SolidConductivityLaw<'a>,
    storage: PreparedHeatStorage,
    reference_density_kg_m3: f64,
    reference_temperature_k: f64,
}

impl<'a> SolidConductivityLaw<'a> {
    /// Bind the SAME original density and specific-heat claims to exact storage.
    ///
    /// `storage_temperature_range_k` is an explicit supported single-phase
    /// interval. The existing material capacity owner retains every original
    /// source breakpoint and integrates it directly; no sampled Cp averages or
    /// resampled enthalpy table are introduced. Both property names must occur
    /// in the retained complete source bundle. Only temperature changes, and
    /// the Cp claim is pinned to the original selected claim throughout.
    pub fn with_heat_storage(
        self,
        density_property: &str,
        specific_heat_property: &str,
        storage_temperature_range_k: [f64; 2],
    ) -> Result<SolidThermalLaw<'a>, SolidThermalError> {
        let invalid = SolidThermalError::Invalid;
        let density = self
            .reference
            .property(density_property)
            .ok_or(invalid("solid storage needs its original density property"))?;
        let specific_heat = self
            .reference
            .property(specific_heat_property)
            .ok_or(invalid(
                "solid storage needs its original specific-heat property",
            ))?;
        if density.requirement().dims() != Density::DIMS
            || specific_heat.requirement().dims() != SPECIFIC_HEAT_CAPACITY_DIMS
            || !density.is_static_at_fixed_external_conditions()
            || !specific_heat.is_static_at_fixed_external_conditions()
            || !density.value_si().is_finite()
            || density.value_si() <= 0.0
            || !specific_heat.value_si().is_finite()
            || specific_heat.value_si() <= 0.0
        {
            return Err(invalid(
                "solid storage requires positive static density and specific heat at fixed external conditions",
            ));
        }
        let reference_density_kg_m3 = density.value_si();
        let reference_temperature_k = self.reference.query().axes()[&self.temperature_axis];
        let [low, high] = storage_temperature_range_k;
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || reference_temperature_k < low
            || reference_temperature_k > high
        {
            return Err(invalid(
                "solid storage needs an ordered positive temperature span containing its reference state",
            ));
        }
        let lower = self
            .reference
            .query()
            .with(&self.temperature_axis, low)
            .map_err(|_| {
                invalid("solid storage lower temperature violates the original source context")
            })?;
        let upper = self
            .reference
            .query()
            .with(&self.temperature_axis, high)
            .map_err(|_| {
                invalid("solid storage upper temperature violates the original source context")
            })?;
        let selected = specific_heat.answer().receipt.selected;
        let claim = self.card.claims().claim(selected).ok_or(invalid(
            "solid storage Cp claim is absent from its original card",
        ))?;
        let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
            id: 0,
            card: self.card,
            capacity: &claim.key,
            amount: HeatStorageAmount::MassKg(reference_density_kg_m3),
            temperature_axis: &self.temperature_axis,
            lower: &lower,
            upper: &upper,
            selection: ClaimSelection::Pinned(selected),
        }])?;
        storage.heat_capacity_j_k(reference_temperature_k)?;
        Ok(SolidThermalLaw {
            conductivity: self,
            storage,
            reference_density_kg_m3,
            reference_temperature_k,
        })
    }
}

impl<'a> SolidThermalLaw<'a> {
    /// Original pinned conductivity and complete source-bundle owner.
    pub const fn conductivity(&self) -> &SolidConductivityLaw<'a> {
        &self.conductivity
    }
    /// Exact source capacity/integral for one cubic metre of reference material.
    /// Multiply by the native reference voxel volume, never by current density.
    pub const fn storage(&self) -> &PreparedHeatStorage {
        &self.storage
    }
    /// Original density that fixes each voxel's invariant mass.
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density_kg_m3
    }

    /// Existing topology/radiation metadata; the nonlinear heat solve replaces
    /// these reference coefficients with the actual source state each iteration.
    pub(crate) fn reference_material(&self) -> Result<SolidMaterial, SolidThermalError> {
        let k = self.conductivity.reference_conductivity_w_m_k();
        Ok(SolidMaterial::new("source-bound thermal solid", k[0])
            .with_orthotropic(k)
            .with_heat_capacity(
                self.storage
                    .heat_capacity_j_k(self.reference_temperature_k)?,
            ))
    }

    pub(crate) fn at_temperature(
        &self,
        cell: usize,
        material_index: u16,
        temperature_k: f64,
    ) -> Result<SolidConductivitySample, SolidThermalError> {
        self.storage.heat_capacity_j_k(temperature_k)?;
        Ok(self
            .conductivity
            .at_temperature(cell, material_index, temperature_k)?)
    }
}
