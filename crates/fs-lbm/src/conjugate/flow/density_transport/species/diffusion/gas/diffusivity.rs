//! Explicit common mixing data or the original binary molecular source.
//!
//! Binary Fick diffusion uses its own m²/s material property. Viscosity and
//! thermal conductivity do not determine it without another declared model.
//! Sampling follows the actual reconstructed T/p at each prepared stage.
//! Soret, barodiffusion and arbitrary multicomponent friction remain separate.

use super::super::{CancelGate, ChtError, poll, refuse};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
};
use fs_qty::QuantitySpec;
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Original authority for the scalar coefficient in the native mass frame.
///
/// Molecular binary data are admitted only for their exact two species.
/// A larger mixture requires an explicitly supplied common mixing field;
/// it cannot silently use one pair coefficient for every constituent.
#[derive(Clone, Copy, Debug)]
pub enum GasSpeciesDiffusivity<'a> {
    /// Authored m²/s values in original voxel order, held for the interval.
    Frozen(&'a [f64]),
    /// Original two-species material law, re-evaluated at each stage's T/p.
    Binary(&'a BinaryDiffusionMaterial<'a>),
}

impl GasSpeciesDiffusivity<'_> {
    pub(super) fn validate_species<const N: usize>(
        &self,
        species: &[SpeciesId; N],
    ) -> Result<(), ChtError> {
        if N == 0 {
            return Err(refuse("gas diffusion requires a nonempty species bank"));
        }
        if let Self::Binary(binding) = self {
            if N != 2 || binding.species[0] == binding.species[1] {
                return Err(refuse(
                    "binary molecular diffusion requires exactly two distinct species",
                ));
            }
            let direct = binding.species[0] == species[0] && binding.species[1] == species[1];
            let reversed = binding.species[0] == species[1] && binding.species[1] == species[0];
            if !direct && !reversed {
                return Err(refuse(
                    "binary diffusivity species do not match the original gas bank",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn sample(
        &self,
        cell: usize,
        temperature_k: f64,
        pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        poll(gate)?;
        if !(temperature_k.is_finite()
            && temperature_k > 0.0
            && pressure_pa.is_finite()
            && pressure_pa > 0.0)
        {
            return Err(refuse(
                "diffusivity sampling requires positive physical temperature and pressure",
            ));
        }
        let value = match self {
            Self::Frozen(values) => *values
                .get(cell)
                .ok_or_else(|| refuse("diffusivity field does not cover the original gas cell"))?,
            Self::Binary(binding) => {
                let typed =
                    |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
                let point = binding
                    .context
                    .clone()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        temperature_k,
                    )
                    .and_then(|point| {
                        point.with_quantity("pressure", typed(QuantityKind::Pressure), pressure_pa)
                    })
                    .map_err(|error| refuse(&format!("binary diffusivity coordinates: {error}")))?;
                let requirement = ScalarPropertyRequirement::try_new(
                    BINARY_DIFFUSIVITY_PROPERTY,
                    BINARY_DIFFUSIVITY_DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .map_err(|error| refuse(&format!("binary diffusivity requirement: {error}")))?;
                let resolved = resolve_material_state_point(
                    binding.material,
                    &point,
                    std::slice::from_ref(&requirement),
                    binding.selection.clone(),
                )
                .map_err(|error| refuse(&format!("binary diffusivity source: {error}")))?;
                resolved
                    .property(BINARY_DIFFUSIVITY_PROPERTY)
                    .ok_or_else(|| refuse("original binary diffusivity is missing"))?
                    .value_si()
            }
        };
        if !value.is_finite() || value < 0.0 {
            return Err(refuse("gas diffusivity must be finite and nonnegative"));
        }
        poll(gate)?;
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
