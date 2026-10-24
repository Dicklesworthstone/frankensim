//! Spatial linear magnetic thermodynamics at fixed flux density and geometry.
//! Each cell owns psi=V|B|²/(2mu(T)); it uses exactly the same affine-source
//! admission and Helmholtz derivative functions as native lumped field storage.

use super::{
    FieldStorageError, FieldThermalState, affine_path, capacity_bound, energy, finite, positive,
    thermal_state,
};
use crate::{
    field_storage::{MAGNETIC_PERMEABILITY_PROPERTY, linear_magnetic_permeability},
    state_point::ResolvedMaterialStatePoint,
};
use fs_matdb::MaterialCard;

#[derive(Clone, Debug, PartialEq)]
pub struct ThermalMagneticMaterial {
    card: MaterialCard,
    material: ResolvedMaterialStatePoint,
    temperature_axis: String,
    temperature_range_k: [f64; 2],
    reference_temperature_k: f64,
    slope_h_m_k: f64,
    minimum_permeability_h_m: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalMagneticMaterial {
    material: ResolvedMaterialStatePoint,
    temperature_k: f64,
    permeability_h_m: f64,
    slope_h_m_k: f64,
}

impl ThermalMagneticMaterial {
    /// Pin a complete original source bundle and positive globally affine mu(T)
    /// over the supplied absolute-temperature range. Other context stays fixed.
    pub fn new(
        card: &MaterialCard,
        material: &ResolvedMaterialStatePoint,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        linear_magnetic_permeability(material)?;
        let (reference_temperature_k, slope_h_m_k, minimum_permeability_h_m) = affine_path(
            card,
            material,
            MAGNETIC_PERMEABILITY_PROPERTY,
            temperature_axis,
            temperature_range_k,
        )?;
        Ok(Self {
            card: card.clone(),
            material: material.clone(),
            temperature_axis: temperature_axis.to_owned(),
            temperature_range_k,
            reference_temperature_k,
            slope_h_m_k,
            minimum_permeability_h_m,
        })
    }
    pub fn card(&self) -> &MaterialCard {
        &self.card
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }
    pub fn permeability_slope_h_m_k(&self) -> f64 {
        self.slope_h_m_k
    }

    /// Refresh the complete source bundle, preserving source selections and all
    /// nonthermal coordinates. Missing data or a changed card never substitutes.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalMagneticMaterial, FieldStorageError> {
        if !temperature_k.is_finite()
            || !(self.temperature_range_k[0]..=self.temperature_range_k[1]).contains(&temperature_k)
        {
            return Err(FieldStorageError::Invalid(
                "magnetic temperature is outside its supported source interval",
            ));
        }
        let point = self
            .material
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| FieldStorageError::Invalid("magnetic thermal query is unrepresentable"))?;
        let material = self.material.at_query_point(&self.card, &point)?;
        let permeability_h_m = linear_magnetic_permeability(&material)?;
        let expected = self.slope_h_m_k.mul_add(
            temperature_k - self.reference_temperature_k,
            linear_magnetic_permeability(&self.material)?,
        );
        if !expected.is_finite()
            || (expected - permeability_h_m).abs()
                > 512.0 * f64::EPSILON * expected.abs().max(permeability_h_m.abs())
        {
            return Err(FieldStorageError::Invalid(
                "magnetic source left its admitted affine law",
            ));
        }
        Ok(ResolvedThermalMagneticMaterial {
            material,
            temperature_k,
            permeability_h_m,
            slope_h_m_k: self.slope_h_m_k,
        })
    }

    /// Full-interval lower bound on fixed-B field heat capacity for the actual
    /// cell volume and declared |B| envelope. The thermal owner must sum all cell
    /// contributions and add its caloric capacity before accepting the model.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_flux_density_t: f64,
        volume_m3: f64,
    ) -> Result<f64, FieldStorageError> {
        positive(volume_m3, "magnetic volume must be finite positive")?;
        if !maximum_flux_density_t.is_finite() || maximum_flux_density_t < 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic flux density envelope must be finite nonnegative",
            ));
        }
        let worst = finite(
            energy(maximum_flux_density_t, self.minimum_permeability_h_m)? * volume_m3,
            "magnetic capacity energy bound is unrepresentable",
        )?;
        if maximum_flux_density_t != 0.0 && worst == 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic capacity energy bound underflows",
            ));
        }
        let ratio = finite(
            self.slope_h_m_k / self.minimum_permeability_h_m,
            "magnetic capacity derivative is unrepresentable",
        )?;
        if self.slope_h_m_k != 0.0 && ratio == 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic capacity derivative underflows",
            ));
        }
        capacity_bound(worst, self.temperature_range_k[1], ratio)
    }
}
impl ResolvedThermalMagneticMaterial {
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn permeability_h_m(&self) -> f64 {
        self.permeability_h_m
    }

    /// Integrated field potentials at the actual local gauge-invariant |B|.
    /// Signed entropy/internal-energy contributions use the zero-field caloric
    /// chart as reference; neither contribution is a standalone positive store.
    pub fn state(
        &self,
        flux_density_t: f64,
        volume_m3: f64,
    ) -> Result<FieldThermalState, FieldStorageError> {
        positive(volume_m3, "magnetic volume must be finite positive")?;
        if !flux_density_t.is_finite() || flux_density_t < 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic magnitude must be finite nonnegative",
            ));
        }
        let psi = finite(
            energy(flux_density_t, self.permeability_h_m)? * volume_m3,
            "magnetic cell free energy is unrepresentable",
        )?;
        if flux_density_t != 0.0 && psi == 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic cell free energy underflows",
            ));
        }
        let ratio = finite(
            self.slope_h_m_k / self.permeability_h_m,
            "magnetic thermal derivative is unrepresentable",
        )?;
        if self.slope_h_m_k != 0.0 && ratio == 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic thermal derivative underflows",
            ));
        }
        thermal_state(psi, self.temperature_k, ratio)
    }
}
