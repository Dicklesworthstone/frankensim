//! Explicit separable reversible magnetism: H(B,T)=s(T) H_ref(B).
//!
//! The caller chooses this constitutive reduction and supplies BOTH original
//! source claims. The reference curve is temperature independent; a separate
//! positive dimensionless scale is affine on the admitted temperature interval
//! and equals one at the reference temperature. No general H(B,T) surface is
//! inferred or fitted. For psi=s psi_ref, S=-s' psi_ref, U=psi+T S and C_B=0.

use super::*;
use crate::field_storage::thermal::{FieldThermalState, affine_source_path};
use fs_matdb::PropertyKey;
use fs_qty::Dims;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct TemperatureScale {
    pub property: String,
    pub value: f64,
    pub slope_per_k: f64,
    pub reference_temperature_k: f64,
    pub temperature_range_k: [f64; 2],
}

/// Source-pinned separable model, not a general two-coordinate magnetic law.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalReversibleMagneticMaterial {
    reference: ReversibleMagneticMaterial,
    material: ResolvedMaterialStatePoint,
    temperature_axis: String,
    scale: TemperatureScale,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermalReversibleMagneticMaterial {
    magnetic: ReversibleMagneticMaterial,
    temperature_k: f64,
}

fn finite(value: f64) -> Result<f64, FieldStorageError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FieldStorageError::Invalid(
            "separable magnetic thermal arithmetic is unrepresentable",
        ))
    }
}

impl ThermalReversibleMagneticMaterial {
    /// Declare H=s(T)H_ref(B) using the SAME complete source bundle, including
    /// any original conductor/thermal roles. No requirement or pin is dropped.
    /// The separate scale must be positive dimensionless data with s(Tref)=1.
    /// The original H_ref retains H(0)=0, positive differential reluctivity and
    /// complete source support over the B/T strip. A positive affine scale
    /// preserves monotonicity at EVERY interior temperature, not just samples.
    pub fn new(
        reference: &ReversibleMagneticMaterial,
        material: &ResolvedMaterialStatePoint,
        scale_property: &PropertyKey,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
    ) -> Result<Self, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        let source = reference
            .material
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .ok_or(invalid("reference magnetic source is absent"))?;
        if reference.temperature_scale.is_some()
            || source
                .source_axes()
                .iter()
                .any(|(axis, dims)| axis == temperature_axis || *dims == fs_qty::Temperature::DIMS)
            || scale_property.quantity() != QuantitySpec::dimensional(Dims([0; 6]))
            || scale_property.name() == REVERSIBLE_FIELD_STRENGTH_PROPERTY
        {
            return Err(invalid(
                "separable magnetic scaling needs an unscaled temperature-independent reference and distinct dimensionless scale",
            ));
        }
        reference.admit_fixed_material_refresh(material)?;
        let checked = material.at_query_point(&reference.card, &material.query())?;
        super::super::admit_source_temperature_refresh(material, &checked)?;
        if material.query() != reference.material.query()
            || material
                .properties()
                .iter()
                .zip(checked.properties())
                .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "reference magnetic law and complete thermal bundle must start at the same original point",
            ));
        }
        let property = material.property(scale_property.name()).ok_or(invalid(
            "complete magnetic bundle lacks its explicitly named scale",
        ))?;
        if property.requirement().quantity() != scale_property.quantity()
            || property.value_si() != 1.0
        {
            return Err(invalid(
                "reference magnetic scale must be dimensionless and exactly one at the reference temperature",
            ));
        }
        let (reference_temperature_k, slope_per_k, _) = affine_source_path(
            &reference.card,
            material,
            scale_property.name(),
            temperature_axis,
            temperature_range_k,
            true,
        )?;
        let law = Self {
            reference: reference.clone(),
            material: material.clone(),
            temperature_axis: temperature_axis.to_owned(),
            scale: TemperatureScale {
                property: scale_property.name().to_owned(),
                value: 1.0,
                slope_per_k,
                reference_temperature_k,
                temperature_range_k,
            },
        };
        // Canonical full-B admission at both endpoints, together with H_ref's
        // absent temperature value axis and rectangular validity, covers the
        // complete strip. No endpoint-only claim about a triangular H(B,T) fit.
        for temperature in temperature_range_k {
            law.at_temperature(temperature)?;
        }
        Ok(law)
    }

    pub fn card(&self) -> &MaterialCard {
        &self.reference.card
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    pub fn temperature_range_k(&self) -> [f64; 2] {
        self.scale.temperature_range_k
    }
    pub fn reference(&self) -> &ReversibleMagneticMaterial {
        &self.reference
    }
    pub fn scale_slope_per_k(&self) -> f64 {
        self.scale.slope_per_k
    }

    /// Refresh every original source role, then scale the original field law.
    /// H_ref receipts remain reference-curve evidence: no synthesized H(B,T)
    /// samples or replacement material claims are constructed.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ResolvedThermalReversibleMagneticMaterial, FieldStorageError> {
        let invalid = FieldStorageError::Invalid;
        if !temperature_k.is_finite()
            || !(self.scale.temperature_range_k[0]..=self.scale.temperature_range_k[1])
                .contains(&temperature_k)
        {
            return Err(invalid(
                "separable magnetic temperature left its complete source interval",
            ));
        }
        let query = self
            .material
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("separable magnetic thermal query is unrepresentable"))?;
        let material = self.material.at_query_point(&self.reference.card, &query)?;
        self.reference.admit_fixed_material_refresh(&material)?;
        let value = material
            .property(&self.scale.property)
            .expect("pinned magnetic scale")
            .value_si();
        let expected = finite(
            self.scale
                .slope_per_k
                .mul_add(temperature_k - self.scale.reference_temperature_k, 1.0),
        )?;
        if !value.is_finite()
            || value <= 0.0
            || expected <= 0.0
            || (value - expected).abs() > 512.0 * f64::EPSILON * value.abs().max(expected.abs())
        {
            return Err(invalid(
                "magnetic scale left its admitted positive affine source path",
            ));
        }
        let mut scale = self.scale.clone();
        scale.value = value;
        let mut magnetic = self.reference.clone();
        magnetic.material = material;
        magnetic.temperature_scale = Some(scale);
        for &slope in &magnetic.slopes_a_m_t {
            if finite(slope * value)? <= 0.0 {
                return Err(invalid("scaled differential reluctivity underflows"));
            }
        }
        magnetic.evaluate(magnetic.maximum_flux_density_t())?;
        Ok(ResolvedThermalReversibleMagneticMaterial {
            magnetic,
            temperature_k,
        })
    }

    /// Exactly zero fixed-B field capacity for this declared affine-scale law.
    /// The shared thermal owner still admits its positive caloric capacity and
    /// original field envelope; zero is not a missing-data fallback.
    pub fn minimum_heat_capacity_j_k(
        &self,
        maximum_flux_density_t: f64,
        volume_m3: f64,
    ) -> Result<f64, FieldStorageError> {
        if !volume_m3.is_finite()
            || volume_m3 <= 0.0
            || !maximum_flux_density_t.is_finite()
            || maximum_flux_density_t < 0.0
            || maximum_flux_density_t > self.reference.maximum_flux_density_t()
        {
            return Err(FieldStorageError::Invalid(
                "finite positive magnetic volume and source-supported B envelope required",
            ));
        }
        for temperature in self.temperature_range_k() {
            self.at_temperature(temperature)?
                .state(maximum_flux_density_t, volume_m3)?;
        }
        Ok(0.0)
    }
}

impl ResolvedThermalReversibleMagneticMaterial {
    pub fn magnetic(&self) -> &ReversibleMagneticMaterial {
        &self.magnetic
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.magnetic.material()
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Zero-field caloric data own the B-independent thermal reference. Field
    /// entropy/internal energy may be signed relative to that datum; neither
    /// is a separate nonnegative store or an additional Joule source.
    pub fn state(
        &self,
        flux_density_t: f64,
        volume_m3: f64,
    ) -> Result<FieldThermalState, FieldStorageError> {
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(FieldStorageError::Invalid(
                "magnetic volume must be finite positive",
            ));
        }
        let scale = self
            .magnetic
            .temperature_scale
            .as_ref()
            .expect("prepared separable model");
        let free_energy_j =
            finite(self.magnetic.evaluate(flux_density_t)?.energy_density_j_m3 * volume_m3)?;
        let reference_energy = finite(
            self.magnetic
                .reference_response(flux_density_t)?
                .energy_density_j_m3
                * volume_m3,
        )?;
        let entropy_j_k = finite(-scale.slope_per_k * reference_energy)?;
        if flux_density_t > 0.0 && (free_energy_j == 0.0 || reference_energy == 0.0)
            || scale.slope_per_k != 0.0 && reference_energy != 0.0 && entropy_j_k == 0.0
        {
            return Err(FieldStorageError::Invalid(
                "separable magnetic thermal energy or entropy underflows",
            ));
        }
        let internal_energy_j = finite(self.temperature_k.mul_add(entropy_j_k, free_energy_j))?;
        Ok(FieldThermalState {
            free_energy_j,
            entropy_j_k,
            internal_energy_j,
            heat_capacity_j_k: 0.0,
        })
    }
}
