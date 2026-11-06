//! Explicit SI coefficient roles for one original dielectric switching fit.
//!
//! Source evidence is retained unchanged; a permittivity, material name or
//! major-loop coercivity alone is not a complete switching model. Only constant
//! coefficients on fixed external source coordinates are admitted. Temperature
//! changes cannot silently replace the fit while retaining incompatible memory.

use super::{
    DielectricHysteresisBranch, FieldStorageError, ScalarDielectricHysteresis,
    UniformHystereticCapacitor,
};
use crate::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS,
    dielectric::reversible::{ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS},
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{EvaluationDecision, MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};
use std::sync::Arc;

/// Exact source property keys for one branch; repeated keys may share a fit.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricHysteresisBranchSources {
    pub susceptibility: PropertyKey,
    pub pinning_field: PropertyKey,
    pub saturation_polarization: PropertyKey,
}

/// All original fitted roles, in stable branch order. No geometry is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricHysteresisSources {
    pub background_permittivity: PropertyKey,
    pub field_limit: PropertyKey,
    pub branches: Vec<DielectricHysteresisBranchSources>,
}

/// Original material card, complete query, selected claims and immutable law.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDielectricHysteresis {
    card: Arc<MaterialCard>,
    material: ResolvedMaterialStatePoint,
    roles: DielectricHysteresisSources,
    law: ScalarDielectricHysteresis,
}

fn invalid(what: &'static str) -> FieldStorageError {
    FieldStorageError::Invalid(what)
}

impl ResolvedDielectricHysteresis {
    /// Resolve a full constant fit using the existing source selection owner.
    /// Numeric time, field, frequency, stress and strain applicability are not
    /// silently frozen at the preparation point. Every role has scalar SI units.
    pub fn resolve(
        card: Arc<MaterialCard>,
        point: &QueryPoint,
        roles: DielectricHysteresisSources,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FieldStorageError> {
        if roles.branches.len() > 1024 {
            return Err(invalid(
                "dielectric switching source exceeds native branch capacity",
            ));
        }
        let mut keys = Vec::<PropertyKey>::new();
        let mut requirements = Vec::new();
        let mut require = |key: &PropertyKey, dims: Dims, admissibility| {
            if key.quantity() != QuantitySpec::dimensional(dims)
                || key.hardness_test().is_some()
                || key.elastic_component().is_some()
                || key.strain_component().is_some()
                || key.stress_component().is_some()
            {
                return Err(invalid(
                    "dielectric switching roles require explicit scalar SI quantities",
                ));
            }
            if !keys.contains(key) {
                requirements.push(ScalarPropertyRequirement::try_with_key(key, admissibility)?);
                keys.push(key.clone());
            }
            Ok(())
        };
        require(
            &roles.background_permittivity,
            ELECTRIC_PERMITTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        require(
            &roles.field_limit,
            ELECTRIC_FIELD_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        for branch in &roles.branches {
            require(
                &branch.susceptibility,
                ELECTRIC_PERMITTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            require(
                &branch.pinning_field,
                ELECTRIC_FIELD_DIMS,
                ScalarAdmissibility::NonNegative,
            )?;
            require(
                &branch.saturation_polarization,
                ELECTRIC_DISPLACEMENT_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?;
        }
        let material = resolve_material_state_point(&card, point, &requirements, selection)?;
        for property in material.properties() {
            if !property.source_axes().is_empty()
                || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
                || !property.is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "dielectric switching memory requires constant fitted coefficients on fixed external source coordinates",
                ));
            }
        }
        let value = |key: &PropertyKey| {
            material
                .property_by_key(key)
                .expect("all original switching roles were resolved")
                .value_si()
        };
        let law = ScalarDielectricHysteresis::new(
            value(&roles.background_permittivity),
            value(&roles.field_limit),
            roles
                .branches
                .iter()
                .map(|branch| DielectricHysteresisBranch {
                    susceptibility_f_m: value(&branch.susceptibility),
                    pinning_field_v_m: value(&branch.pinning_field),
                    saturation_polarization_c_m2: value(&branch.saturation_polarization),
                })
                .collect(),
        )?;
        Ok(Self {
            card,
            material,
            roles,
            law,
        })
    }

    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn sources(&self) -> &DielectricHysteresisSources {
        &self.roles
    }
    pub const fn law(&self) -> &ScalarDielectricHysteresis {
        &self.law
    }

    /// Admit one owned absolute-temperature coordinate without inventing a
    /// temperature-dependent switching law. Every selected coefficient must
    /// explicitly cover the complete interval and remain the same constant fit.
    pub fn admit_temperature_interval(
        &self,
        axis: &str,
        interval_k: [f64; 2],
    ) -> Result<(), FieldStorageError> {
        let [low, high] = interval_k;
        let quantity = self.material.axis_quantities().get(axis).copied();
        let absolute_temperature =
            quantity
                .and_then(QuantitySpec::semantic_type)
                .is_some_and(|ty| {
                    ty.kind() == QuantityKind::AbsoluteTemperature && ty.form() == ValueForm::Static
                })
                || axis == "T"
                    && quantity.is_none_or(|q| q == QuantitySpec::dimensional(Temperature::DIMS));
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || !absolute_temperature
            || !self
                .material
                .query_point()
                .iter()
                .any(|(name, _)| name == axis)
        {
            return Err(invalid(
                "dielectric switching needs one owned positive absolute-temperature interval",
            ));
        }
        for property in self.material.properties() {
            if !property.has_only_fixed_external_conditions(Some(axis), false)
                || !property
                    .answer()
                    .evidence
                    .model
                    .validity
                    .bounds()
                    .contains_key(axis)
            {
                return Err(invalid(
                    "every dielectric coefficient needs the same explicitly source-bounded temperature",
                ));
            }
        }
        for temperature in interval_k {
            let point = self
                .material
                .query()
                .with(axis, temperature)
                .map_err(|_| invalid("dielectric temperature query is invalid"))?;
            self.at_query_point(&point)?;
        }
        Ok(())
    }

    /// Recheck the same selected constant fit without resetting switching memory.
    pub fn at_query_point(&self, point: &QueryPoint) -> Result<Self, FieldStorageError> {
        let material = self.material.at_query_point(&self.card, point)?;
        if self
            .material
            .properties()
            .iter()
            .zip(material.properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "changing dielectric switching coefficients needs a thermal/history law",
            ));
        }
        Ok(Self {
            material,
            ..self.clone()
        })
    }
}

impl UniformHystereticCapacitor {
    /// Bind the original source fit to declared plate area and gap.
    pub fn from_source(
        source: ResolvedDielectricHysteresis,
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, FieldStorageError> {
        let mut capacitor = Self::new(source.law.clone(), area_m2, gap_m)?;
        capacitor.source = Some(source);
        Ok(capacitor)
    }

    /// Original fit and source receipts; None denotes an explicitly authored fit.
    pub const fn source(&self) -> Option<&ResolvedDielectricHysteresis> {
        self.source.as_ref()
    }

    /// Recheck the original constant fit at an explicitly owned source query.
    /// Geometry and the constitutive law remain unchanged; accepted electrical
    /// history belongs to the caller and must not be reinitialized on refresh.
    pub fn at_source_query(&self, point: &QueryPoint) -> Result<Self, FieldStorageError> {
        let source = self
            .source
            .as_ref()
            .ok_or_else(|| invalid("authored dielectric capacitor has no source query"))?
            .at_query_point(point)?;
        Ok(Self {
            source: Some(source),
            ..self.clone()
        })
    }
}
