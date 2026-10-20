//! Material-state ingress for a fresh enclosure equilibrium evaluation.
//! Geometry and the radiosity equations remain with the original enclosure.

use super::super::{ConductionError, GrayDiffuseEnclosure, SurfaceEmissivity,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, radiation_error};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

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
        let point = binding.at(phase).map_err(|e| ConductionError::MaterialQuery {
            property: SURFACE_EMISSIVITY_PROPERTY.to_owned(),
            temperature: phase.temperature_k(), upstream: e.to_string(),
        })?;
        Self::from_condition(surface_name, ThermalMaterialPoint { context: &point, ..binding }, selection)
    }

    /// Resolve the same static surface law from an explicit complete condition
    /// point, without fabricating a caloric curve or phase state for a boundary.
    /// The binding's existing absolute-temperature coordinate is the physical
    /// SURFACE temperature. All pressure, processing and other coordinates keep
    /// their exact descriptors. The caller owns any phase-fraction coordinate;
    /// enthalpy consumers should use `from_thermal_state` for its identity check.
    ///
    /// This shares property resolution, static-law admission and receipt
    /// construction with the phase-state path. The returned coefficient is a
    /// frozen sample; re-resolve after a surface condition changes.
    ///
    /// # Errors
    /// Absent/nonpositive/nonabsolute temperature, blank surface, missing or
    /// inadmissible emissivity and unsupported source state refuse atomically.
    pub fn from_condition(
        surface_name: &str, binding: ThermalMaterialPoint<'_>,
        selection: MaterialPropertySelection,
    ) -> Result<(Self, ResolvedMaterialStatePoint), ConductionError> {
        let temperature = binding.context.axes().get(binding.temperature_axis).copied()
            .unwrap_or(f64::NAN);
        let failure = |error: String| ConductionError::MaterialQuery {
            property: SURFACE_EMISSIVITY_PROPERTY.to_owned(),
            temperature, upstream: error,
        };
        if surface_name.trim().is_empty() {
            return Err(radiation_error("<unnamed>", "surface name is blank",
                "name the radiating surface"));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        if !temperature.is_finite() || temperature <= 0.0
            || binding.context.axis_quantities().get(binding.temperature_axis)
                .is_some_and(|q| *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS))
        {
            return Err(failure("an existing finite positive absolute-temperature coordinate is required".into()));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, ScalarAdmissibility::StrictlyPositive,
        ).map_err(|e| failure(e.to_string()))?;
        let source = resolve_material_state_point(binding.material, binding.context, &[requirement], selection)
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
            temperature_k: temperature,
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

#[cfg(test)]
mod condition_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};

    fn card() -> MaterialCard {
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve { abscissa: "Ts".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.2), (400.0, 0.8)], dims: EMISSIVITY_DIMS },
            validity: ValidityDomain::unconstrained().with("Ts", 300.0, 400.0)
                .with("pressure", 100_000.0, 110_000.0),
            interpolation: InterpolationPolicy::LinearInside, uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(), provenance: Provenance { source: "synthetic surface fixture".into(),
                license: "CC0: authored test data".into(), artifact: None },
        }).unwrap();
        MaterialCard::assemble(MaterialStateId { chemistry: "synthetic surface".into(),
            phase: "solid".into(), process: "fixture".into(), revision: 0 }, claims, Vec::new()).unwrap()
    }

    #[test]
    fn g1_condition_surface_temperature_and_pressure_share_original_material_owner() {
        let card = card();
        let mut identities = Vec::new();
        for pressure in [101_325.0, 105_000.0] {
            let point = QueryPoint::new().with("Ts", 350.0).unwrap().with("pressure", pressure).unwrap();
            let (emissivity, state) = SurfaceEmissivity::from_condition("surface",
                ThermalMaterialPoint { material: &card, context: &point,
                    temperature_axis: "Ts", liquid_fraction_axis: None },
                MaterialPropertySelection::SingleClaimOnly).unwrap();
            assert!((emissivity.value() - 0.5).abs() < 1e-14);
            assert_eq!(emissivity.temperature_k(), 350.0);
            assert_eq!(state.card_identity(), card.content_hash());
            assert_eq!(state.query_point(), &point);
            identities.push(state.identity());
        }
        assert_ne!(identities[0], identities[1]);
    }

    #[test]
    fn g0_missing_pressure_or_invalid_surface_temperature_cannot_use_a_default() {
        let card = card();
        for point in [QueryPoint::new().with("Ts", 350.0).unwrap(),
            QueryPoint::new().with("Ts", 450.0).unwrap().with("pressure", 101_325.0).unwrap(),
            QueryPoint::new().with("pressure", 101_325.0).unwrap()] {
            assert!(SurfaceEmissivity::from_condition("surface", ThermalMaterialPoint {
                material: &card, context: &point, temperature_axis: "Ts", liquid_fraction_axis: None,
            }, MaterialPropertySelection::SingleClaimOnly).is_err());
        }
    }
}
