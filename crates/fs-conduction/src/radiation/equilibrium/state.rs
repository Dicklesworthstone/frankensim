//! Material-state ingress for a fresh enclosure equilibrium evaluation.
//! Geometry and the radiosity equations remain with the original enclosure.

use super::super::{ConductionError, GrayDiffuseEnclosure, SurfaceEmissivity,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, radiation_error};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_material::thermal::source::ThermalMaterialPoint;

impl SurfaceEmissivity {
    /// Resolve total hemispherical emissivity at the SAME material/phase state
    /// as caloric storage and transport. Only the explicitly declared T/fraction
    /// coordinates evolve; pressure, finish/lot conditions and quantity schemas
    /// are retained. No wavelength sample becomes a total emissivity.
    ///
    /// Returns the original complete material bundle alongside the radiation
    /// coefficient. Retain it for source/condition inspection. The coefficient
    /// is a frozen snapshot, NOT a law for a later temperature: re-resolve at
    /// every trial state when the material response changes.
    ///
    /// # Errors
    /// Missing, ambiguous, wrong-unit, dynamic-only, foreign-card or unsupported
    /// material states refuse. Total emissivity must lie in (0,1]; zero is not
    /// silently converted to a black/reflecting material or a dropped surface.
    pub fn from_thermal_state(
        surface_name: &str, binding: ThermalMaterialPoint<'_>,
        phase: EquilibriumPhaseState, selection: MaterialPropertySelection,
    ) -> Result<(Self, ResolvedMaterialStatePoint), ConductionError> {
        let failure = |error: String| ConductionError::MaterialQuery {
            property: SURFACE_EMISSIVITY_PROPERTY.to_owned(),
            temperature: phase.temperature_k(), upstream: error,
        };
        if surface_name.trim().is_empty() {
            return Err(radiation_error("<unnamed>", "surface name is blank",
                "name the radiating surface"));
        }
        let point = binding.at(phase).map_err(|e| failure(e.to_string()))?;
        let requirement = ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, ScalarAdmissibility::StrictlyPositive,
        ).map_err(|e| failure(e.to_string()))?;
        let source = resolve_material_state_point(binding.material, &point, &[requirement], selection)
            .map_err(|e| failure(e.to_string()))?;
        let property = source.property(SURFACE_EMISSIVITY_PROPERTY)
            .ok_or_else(|| failure("resolved total emissivity is absent".into()))?;
        binding.admit_static(property).map_err(|e| failure(e.to_string()))?;
        let value = property.value_si();
        if !value.is_finite() || value <= 0.0 || value > 1.0 {
            return Err(radiation_error(surface_name, "total hemispherical emissivity must lie in (0,1]",
                "provide an admitted static total-emissivity material claim"));
        }
        let answer = property.answer();
        let emissivity = Self {
            value, uncertainty: answer.evidence.value.uncertainty.clone(),
            temperature_k: phase.temperature_k(),
            temperature_validity_k: answer.evidence.model.validity.bound(binding.temperature_axis),
            card_identity: binding.material.content_hash(),
            material_state: binding.material.id().to_string(), receipt: answer.receipt.clone(),
        };
        Ok((emissivity, source))
    }
}

impl GrayDiffuseEnclosure {
    /// Rebind a complete material snapshot while retaining the exact surface
    /// geometry, names, trace slots, view factors and their evidence. Input order
    /// is the existing surface/matrix order. A material change needs a NEW
    /// factorization; no previously prepared operator is reused or mutated.
    ///
    /// # Errors
    /// Incomplete emissivity replacement or original enclosure admission failure.
    pub fn with_emissivities(&self, emissivities: Vec<SurfaceEmissivity>)
        -> Result<Self, ConductionError>
    {
        if emissivities.len() != self.surfaces.len() {
            return Err(radiation_error("enclosure", "one emissivity per original surface is required",
                "resolve the complete material snapshot in surface order"));
        }
        let mut surfaces = self.surfaces.clone();
        for (surface, emissivity) in surfaces.iter_mut().zip(emissivities) {
            surface.emissivity = emissivity;
        }
        Self::new(surfaces, self.view_factors.clone())
    }
}
