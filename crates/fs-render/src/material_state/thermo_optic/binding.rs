//! Thermo-optic observations of an original material's accepted state.
//!
//! One immutable reference binding remains at the SCHOTT law's reference
//! temperature. Each successor reads the named absolute-temperature coordinate
//! of the accepted material and recomputes the native glass from that original
//! reference. The resulting ordinary material binding feeds existing closed
//! tracer primitives and slab emission with the exact accepted-state identity.
//! Reference absorption and surface finish are explicitly held fixed. This
//! supplies no absorption-temperature law, expansion, phase transition, or
//! permission to treat an open regional surface as a closed dielectric volume.

use super::{SchottThermoOpticLaw, ThermoOpticFieldError, schott_glass_at_temperature};
use crate::material_state::{MaterialOpticalBinding, OpticalAppearance};
use fs_blake3::DomainHasher;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq)]
struct Reference {
    binding: MaterialOpticalBinding,
    state: ResolvedMaterialStatePoint,
    law: SchottThermoOpticLaw,
    temperature_axis: String,
}

/// Original material optics and thermal law, with one accepted observation.
///
/// Cloning or refreshing shares the immutable original source. Only endpoint
/// admission evaluates state coordinates; native wavelength queries retain the
/// allocation-free phase/group model and perform no material resolution.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialThermoOpticBinding {
    reference: Arc<Reference>,
    binding: MaterialOpticalBinding,
}

impl MaterialThermoOpticBinding {
    /// Bind an existing dielectric observation at the law's exact reference
    /// temperature. Its mechanical source must be `reference_state` itself.
    ///
    /// The named coordinate must already be typed static absolute temperature;
    /// an untyped `T` or a reservoir temperature is not substituted. The source
    /// law and reference dispersion must describe the same physical specimen;
    /// their association is explicit caller input, not a measured-data claim.
    /// The original absorption and all surface properties remain fixed during
    /// subsequent thermal observations. An already corrected reference refuses.
    pub fn try_new(
        reference: MaterialOpticalBinding,
        reference_state: &ResolvedMaterialStatePoint,
        law: SchottThermoOpticLaw,
        temperature_axis: &str,
    ) -> Result<Self, MaterialThermoOpticBindingError> {
        if reference.accepted_state_identity() != reference_state.identity() {
            return Err(MaterialThermoOpticBindingError::ReferenceStateMismatch);
        }
        if !matches!(reference.appearance(), OpticalAppearance::Dielectric { .. }) {
            return Err(MaterialThermoOpticBindingError::ReferenceAppearance);
        }
        if temperature(reference_state, temperature_axis)? != law.reference_temperature_k() {
            return Err(MaterialThermoOpticBindingError::ReferenceTemperature);
        }
        let original = Reference {
            binding: reference,
            state: reference_state.clone(),
            law,
            temperature_axis: temperature_axis.to_owned(),
        };
        let binding = observe(&original, reference_state)?;
        Ok(Self {
            reference: Arc::new(original),
            binding,
        })
    }

    /// Observe an accepted endpoint on the same card and physical condition,
    /// allowing only the declared material temperature to change.
    ///
    /// Every other coordinate, coordinate name and quantity meaning must agree
    /// with the original source. The accepted bundle may contain the physical
    /// model's own property requirements; optical requirements are not imposed
    /// on mechanics. Refusal leaves this observation and the source unchanged.
    /// No correction accumulates over heating/cooling cycles.
    pub fn at_material_state(
        &self,
        accepted: &ResolvedMaterialStatePoint,
    ) -> Result<Self, MaterialThermoOpticBindingError> {
        let binding = observe(&self.reference, accepted)?;
        Ok(Self {
            reference: Arc::clone(&self.reference),
            binding,
        })
    }

    /// Native observation carrying the actual accepted material-state identity.
    pub const fn binding(&self) -> &MaterialOpticalBinding {
        &self.binding
    }

    pub fn reference_binding(&self) -> &MaterialOpticalBinding {
        &self.reference.binding
    }

    pub fn reference_state(&self) -> &ResolvedMaterialStatePoint {
        &self.reference.state
    }

    pub fn law(&self) -> &SchottThermoOpticLaw {
        &self.reference.law
    }

    pub fn temperature_axis(&self) -> &str {
        &self.reference.temperature_axis
    }
}

fn observe(
    original: &Reference,
    accepted: &ResolvedMaterialStatePoint,
) -> Result<MaterialOpticalBinding, MaterialThermoOpticBindingError> {
    let reference = &original.state;
    if accepted.material() != reference.material()
        || accepted.card_identity() != reference.card_identity()
        || accepted.axis_quantities() != reference.axis_quantities()
        || accepted.query_point().len() != reference.query_point().len()
        || accepted
            .query_point()
            .iter()
            .zip(reference.query_point())
            .any(|((axis, value), (reference_axis, reference_value))| {
                axis != reference_axis
                    || (axis != &original.temperature_axis && value != reference_value)
            })
    {
        return Err(MaterialThermoOpticBindingError::StateMismatch);
    }
    let temperature_k = temperature(accepted, &original.temperature_axis)?;
    let OpticalAppearance::Dielectric { glass, surface } = original.binding.appearance() else {
        return Err(MaterialThermoOpticBindingError::ReferenceAppearance);
    };
    let glass = schott_glass_at_temperature(&glass, &original.law, temperature_k)
        .map_err(MaterialThermoOpticBindingError::Optical)?;
    let mut source = DomainHasher::new("org.frankensim.render.material-thermo-optic-binding.v1");
    source.update(original.binding.identity().as_bytes());
    source.update(original.law.identity().as_bytes());
    source.update(&(original.temperature_axis.len() as u64).to_le_bytes());
    source.update(original.temperature_axis.as_bytes());
    let mut binding = original.binding;
    binding.mechanical_state_identity = accepted.identity();
    binding.optical_state_identity = source.finalize();
    binding.appearance = OpticalAppearance::Dielectric { glass, surface };
    Ok(binding)
}

fn temperature(
    state: &ResolvedMaterialStatePoint,
    axis: &str,
) -> Result<f64, MaterialThermoOpticBindingError> {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if state.axis_quantities().get(axis) != Some(&absolute) {
        return Err(MaterialThermoOpticBindingError::TemperatureAxis);
    }
    state
        .state_coordinate(axis)
        .filter(|value| value.is_finite() && *value > 0.0)
        .ok_or(MaterialThermoOpticBindingError::TemperatureAxis)
}

/// Refusal before a thermo-optic material observation is published.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaterialThermoOpticBindingError {
    ReferenceStateMismatch,
    StateMismatch,
    TemperatureAxis,
    ReferenceTemperature,
    ReferenceAppearance,
    Optical(ThermoOpticFieldError),
}

impl core::fmt::Display for MaterialThermoOpticBindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material thermo-optic binding refused: {self:?}")
    }
}

impl core::error::Error for MaterialThermoOpticBindingError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Optical(error) => Some(error),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests;
