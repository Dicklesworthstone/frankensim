//! Source-selected coefficients for the original rate-independent play model.
//!
//! A material card supplies the fitted background, pinning, saturation and
//! reversible branch slopes. Explicit keys assign their physical roles; no
//! permeability, coercivity or material name is silently treated as a complete
//! hysteresis fit. This admission preserves the source's evidence, including
//! research provenance, rather than promoting an authored fit to measurements.
//!
//! Coefficients are constant within their declared complete source domain.
//! Changing temperature may move within that domain, but cannot change a fit
//! while retaining incompatible magnetic memory. A temperature-dependent fit
//! needs its own free-energy/entropy and history-transfer law and refuses here.

use super::{FieldStorageError, HysteresisBranch, ScalarHysteresis, UniformHystereticCore};
use crate::field_storage::{
    MAGNETIC_PERMEABILITY_DIMS,
    permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_matdb::{EvaluationDecision, MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{Dims, QuantitySpec};
use std::sync::Arc;

/// Exact material properties for one fitted play branch. Repeated keys across
/// branches deliberately share a coefficient; branch multiplicity stays data.
#[derive(Clone, Debug, PartialEq)]
pub struct HysteresisBranchSources {
    /// Reversible polarization slope dz/dH [H/m], not a dimensionless chi.
    pub susceptibility: PropertyKey,
    /// Irreversible pinning threshold [A/m], not an entire loop's coercivity.
    pub pinning_field: PropertyKey,
    /// Bound on this branch's polarization contribution [T].
    pub saturation_polarization: PropertyKey,
}

/// Complete fitted-law roles. Geometry and initial magnetization remain
/// specimen data. Zero branches admits a linear background with sourced bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct HysteresisSources {
    /// Reversible background contribution [H/m], not total observed mu.
    pub background_permeability: PropertyKey,
    /// Positive bound on actual |H| accepted by the fitted model [A/m].
    pub field_limit: PropertyKey,
    /// Original stable branch order, with at most the native 1024 branches.
    pub branches: Vec<HysteresisBranchSources>,
}

/// The original immutable card/selection, complete physical query and native
/// constitutive owner. All magnetic memory still belongs to ScalarHysteresis.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedScalarHysteresis {
    card: Arc<MaterialCard>,
    material: ResolvedMaterialStatePoint,
    roles: HysteresisSources,
    law: ScalarHysteresis,
}

fn invalid(reason: &'static str) -> FieldStorageError {
    FieldStorageError::Invalid(reason)
}

impl ResolvedScalarHysteresis {
    /// Resolve the complete fitted law once, preserving every selected claim.
    /// The shared state-point capacity bounds distinct properties; many branches
    /// may explicitly share keys. Only constant source models with fixed typed
    /// external applicability are admitted. Dynamic field, frequency, stress or
    /// rate restrictions cannot be frozen at the preparation point.
    pub fn resolve(
        card: Arc<MaterialCard>,
        point: &QueryPoint,
        roles: HysteresisSources,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FieldStorageError> {
        if roles.branches.len() > 1024 {
            return Err(invalid("hysteresis source branches exceed native capacity"));
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
                    "hysteresis roles require their explicit scalar SI quantities",
                ));
            }
            if !keys.contains(key) {
                requirements.push(ScalarPropertyRequirement::try_with_key(key, admissibility)?);
                keys.push(key.clone());
            }
            Ok(())
        };
        require(
            &roles.background_permeability,
            MAGNETIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        require(
            &roles.field_limit,
            FIELD_STRENGTH_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?;
        for branch in &roles.branches {
            require(
                &branch.susceptibility,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?;
            require(
                &branch.pinning_field,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::NonNegative,
            )?;
            require(
                &branch.saturation_polarization,
                FLUX_DENSITY_DIMS,
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
                    "hysteretic memory requires source-constant fitted coefficients without unowned dynamic applicability",
                ));
            }
        }
        let value = |key: &PropertyKey| {
            material
                .property_by_key(key)
                .expect("complete requirements were resolved")
                .value_si()
        };
        let law = ScalarHysteresis::new(
            value(&roles.background_permeability),
            value(&roles.field_limit),
            roles
                .branches
                .iter()
                .map(|branch| HysteresisBranch {
                    susceptibility_h_m: value(&branch.susceptibility),
                    pinning_field_a_m: value(&branch.pinning_field),
                    saturation_polarization_t: value(&branch.saturation_polarization),
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

    /// Original material identity, query quantities and selected source receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Original fitted coefficient roles and branch order.
    pub const fn sources(&self) -> &HysteresisSources {
        &self.roles
    }
    /// The existing history owner, also usable for actual preparation cycles.
    pub const fn law(&self) -> &ScalarHysteresis {
        &self.law
    }

    /// Recheck the SAME constant fit at a new complete physical context. No
    /// source is reselected, constitutive memory changed, or energy reset.
    /// Missing source support or altered coordinate quantities refuses.
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
                "hysteretic coefficient changes require a thermal/history law",
            ));
        }
        Ok(Self {
            material,
            ..self.clone()
        })
    }

    /// Admit heating only over an explicitly source-bounded absolute-temperature
    /// interval. Each selected constant claim must cover it; an unconditioned
    /// coefficient is not silently promoted to a temperature-independent fit.
    /// All nonthermal query coordinates retain their original values/quantities.
    pub fn admit_temperature_interval(
        &self,
        axis: &str,
        interval_k: [f64; 2],
    ) -> Result<(), FieldStorageError> {
        let [low, high] = interval_k;
        if !low.is_finite() || !high.is_finite() || low <= 0.0 || high <= low {
            return Err(invalid(
                "positive increasing core temperature interval required",
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
                    "every magnetic coefficient needs the same owned source-bounded temperature",
                ));
            }
        }
        for temperature in interval_k {
            let point = self
                .material
                .query()
                .with(axis, temperature)
                .map_err(|_| invalid("invalid original core temperature coordinate"))?;
            self.at_query_point(&point)?;
        }
        Ok(())
    }
}

impl UniformHystereticCore {
    /// Bind original source data to the native geometry reduction. All current,
    /// linkage, remanence and pinning heat use the existing play/flux inverse.
    pub fn from_source(
        source: ResolvedScalarHysteresis,
        turns: f64,
        area_m2: f64,
        path_length_m: f64,
    ) -> Result<Self, FieldStorageError> {
        let mut core = Self::new(source.law.clone(), turns, area_m2, path_length_m)?;
        core.source = Some(source);
        Ok(core)
    }

    /// Original source binding, or None for an explicitly authored native fit.
    pub const fn source(&self) -> Option<&ResolvedScalarHysteresis> {
        self.source.as_ref()
    }

    /// Retain this exact geometry and magnetic law while rechecking the original
    /// material at a new context. This does not advance/reset internal memory.
    pub fn at_source_query(&self, point: &QueryPoint) -> Result<Self, FieldStorageError> {
        let source = self
            .source
            .as_ref()
            .ok_or_else(|| invalid("authored core has no original material source"))?
            .at_query_point(point)?;
        Ok(Self {
            source: Some(source),
            ..self.clone()
        })
    }
}
