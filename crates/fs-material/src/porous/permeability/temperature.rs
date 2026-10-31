//! Original intrinsic-permeability claims on a fixed skeleton and material frame.
//!
//! This is a constitutive source owner, not a Darcy solver. It changes only the
//! declared absolute temperature and reuses the complete original source bundle.
//! Consumers retain their own pressure, velocity, storage and geometry limits.

use super::{INTRINSIC_PERMEABILITY_DIMS, IntrinsicPermeability, MaterialError, invalid};
use crate::liquid_temperature::validate_temperature;
use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::{MaterialCard, PropertyKey};
use fs_qty::{QuantitySpec, Temperature};

/// A complete source result and its original principal tensor convention.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPermeability {
    material: ResolvedMaterialStatePoint,
    permeability: IntrinsicPermeability,
}
impl ResolvedPermeability {
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn permeability(&self) -> &IntrinsicPermeability {
        &self.permeability
    }
}

/// Temperature-dependent intrinsic permeability, without inferred gas slip,
/// saturation, pressure evolution, changing axes or deformation.
#[derive(Clone, Debug)]
pub struct PermeabilityTemperatureLaw<'a> {
    card: &'a MaterialCard,
    reference: ResolvedPermeability,
    keys: [PropertyKey; 3],
    temperature_axis: String,
}
impl<'a> PermeabilityTemperatureLaw<'a> {
    /// Retain three exact principal property keys and orthonormal world axes.
    /// Repeated keys explicitly declare equal principal values. Zero intrinsic
    /// permeability remains an exact seal; a consumer may refuse sealed axes
    /// when its native momentum model has no matching velocity constraint.
    pub fn try_new(
        card: &'a MaterialCard,
        material: ResolvedMaterialStatePoint,
        keys: [PropertyKey; 3],
        axes: [[f64; 3]; 3],
        temperature_axis: &str,
    ) -> Result<Self, MaterialError> {
        if material.card_identity() != card.content_hash() {
            return Err(invalid(
                "permeability refresh requires its original material card",
            ));
        }
        validate_temperature(&material.query(), temperature_axis)
            .map_err(|e| invalid(&e.to_string()))?;
        let reference = resolve(material, &keys, axes)?;
        Ok(Self {
            card,
            reference,
            keys,
            temperature_axis: temperature_axis.into(),
        })
    }

    /// Pin every original selected source while updating only absolute T.
    pub fn at_temperature(
        &self,
        temperature: Temperature,
    ) -> Result<ResolvedPermeability, MaterialError> {
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(invalid(
                "permeability needs a positive finite absolute temperature",
            ));
        }
        let point = self
            .reference
            .material
            .query()
            .with(&self.temperature_axis, temperature.value())
            .map_err(|e| invalid(&e.to_string()))?;
        let material = self
            .reference
            .material
            .at_query_point(self.card, &point)
            .map_err(|e| invalid(&e.to_string()))?;
        resolve(
            material,
            &self.keys,
            *self.reference.permeability.material_axes_in_world(),
        )
    }

    pub const fn reference(&self) -> &ResolvedPermeability {
        &self.reference
    }
    pub const fn source_keys(&self) -> &[PropertyKey; 3] {
        &self.keys
    }
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
}

fn resolve(
    material: ResolvedMaterialStatePoint,
    keys: &[PropertyKey; 3],
    axes: [[f64; 3]; 3],
) -> Result<ResolvedPermeability, MaterialError> {
    let mut values = [0.0; 3];
    for (value, key) in values.iter_mut().zip(keys) {
        let property = material
            .property_by_key(key)
            .ok_or_else(|| invalid("missing exact principal intrinsic-permeability source"))?;
        if property.requirement().quantity()
            != QuantitySpec::dimensional(INTRINSIC_PERMEABILITY_DIMS)
            || !property.is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "intrinsic permeability requires a static source with exact m2 quantity",
            ));
        }
        *value = property.value_si();
    }
    let permeability = IntrinsicPermeability::orthotropic(values, axes)?;
    Ok(ResolvedPermeability {
        material,
        permeability,
    })
}
