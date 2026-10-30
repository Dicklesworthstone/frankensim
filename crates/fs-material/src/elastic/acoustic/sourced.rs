//! Isotropic pressure and wave response from original, explicitly selected
//! density/bulk/shear claims. This is a frozen elastic model, not a mechanism
//! for inferring isotropy, frequency response or a thermodynamic convention
//! from property names. Redundant E/nu claims are neither averaged nor used.

use super::ElasticAcousticState;
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

/// Explicit three-coefficient cubic symmetry using the same tensor consumers.
pub mod cubic;

/// Exact property keys in density, bulk-modulus, shear-modulus order. Choosing
/// this type explicitly declares the isotropic model; it is unsuitable for
/// directional crystal constants or a homogenization without a stated law.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotropicBulkShearRequirements {
    properties: [ScalarPropertyRequirement; 3],
}

impl IsotropicBulkShearRequirements {
    /// Preserve the source vocabulary and quantity descriptors verbatim.
    pub fn try_new(
        density: ScalarPropertyRequirement,
        bulk_modulus: ScalarPropertyRequirement,
        shear_modulus: ScalarPropertyRequirement,
    ) -> Result<Self, MaterialStatePointError> {
        for (requirement, dims) in [
            (&density, Density::DIMS),
            (&bulk_modulus, Pressure::DIMS),
            (&shear_modulus, Pressure::DIMS),
        ] {
            require_scalar(requirement, dims)?;
        }
        if density.name() == bulk_modulus.name()
            || density.name() == shear_modulus.name()
            || bulk_modulus.name() == shear_modulus.name()
        {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: bulk_modulus.name().into(),
                reason: "density, bulk modulus and shear modulus require distinct source keys",
            });
        }
        Ok(Self {
            properties: [density, bulk_modulus, shear_modulus],
        })
    }

    /// Canonical dimension-only `density`, `bulk_modulus`, `shear_modulus` keys.
    /// Other source names remain supported through `try_new`, without aliases.
    pub fn canonical() -> Result<Self, MaterialStatePointError> {
        let positive = |name, dims| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        };
        Self::try_new(
            positive("density", Density::DIMS)?,
            positive("bulk_modulus", Pressure::DIMS)?,
            positive("shear_modulus", Pressure::DIMS)?,
        )
    }

    /// Original roles, including exact source keys and admissibility.
    pub const fn properties(&self) -> &[ScalarPropertyRequirement; 3] {
        &self.properties
    }
}

/// Original source bundle and the single admitted elastic tensor used for
/// hydrostatic strain, directional waves and ordinary native solid operators.
/// Source uncertainty/covariance remain on the original receipts; the derived
/// tensor and speeds are not independent observations or certified bounds.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedIsotropicBulkShearState {
    material: ResolvedMaterialStatePoint,
    requirements: IsotropicBulkShearRequirements,
    acoustic: ElasticAcousticState,
}

impl ResolvedIsotropicBulkShearState {
    /// Adopt a complete source bundle without choosing new claims or dropping
    /// additional properties. The explicit strain ceiling limits predictions;
    /// it is not a claim of experimentally measured yield strength.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        requirements: &IsotropicBulkShearRequirements,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, MaterialStatePointError> {
        let value = |required: &ScalarPropertyRequirement| {
            material
                .property(required.name())
                .filter(|p| p.requirement() == required)
                .map(|p| p.value_si())
                .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "bulk/shear response needs the exact original scalar requirement",
                })
        };
        let [density, bulk, shear] = &requirements.properties;
        let values = [value(density)?, value(bulk)?, value(shear)?];
        let mut identity = DomainHasher::new("org.frankensim.fs-material.isotropic-bulk-shear.v1");
        identity.update(material.identity().as_bytes());
        for r in &requirements.properties {
            identity.update(&(r.name().len() as u64).to_le_bytes());
            identity.update(r.name().as_bytes());
        }
        identity.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
        let acoustic = ElasticAcousticState::from_bulk_and_shear(
            values[0],
            values[1],
            values[2],
            identity.finalize(),
            maximum_absolute_mandel_strain,
        )
        .map_err(|_| MaterialStatePointError::InvalidDerived {
            quantity: "finite positive isotropic density/bulk/shear response and strain ceiling",
        })?;
        Ok(Self {
            material,
            requirements: requirements.clone(),
            acoustic,
        })
    }

    /// Complete source/card/query bundle, retaining every selected receipt.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Exact independent source parameters, without a rounded E/nu conversion.
    pub const fn requirements(&self) -> &IsotropicBulkShearRequirements {
        &self.requirements
    }
    /// Frozen law shared by pressure, wave and mechanical consumers.
    pub const fn acoustic(&self) -> &ElasticAcousticState {
        &self.acoustic
    }
    /// Binds the complete original source, parameter roles and strain ceiling.
    pub const fn identity(&self) -> ContentHash {
        self.acoustic.material_state_identity()
    }

    /// Refresh every claim from the same immutable card and original selections.
    /// Unsupported states fail without mutating this state or changing sources.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(
            self.material.at_query_point(card, point)?,
            &self.requirements,
            self.acoustic.maximum_absolute_mandel_strain(),
        )
    }
}

/// Select one density/bulk/shear source bundle at a common complete query.
pub fn resolve_isotropic_bulk_shear_state(
    card: &MaterialCard,
    point: &QueryPoint,
    requirements: &IsotropicBulkShearRequirements,
    maximum_absolute_mandel_strain: f64,
    selection: MaterialPropertySelection,
) -> Result<ResolvedIsotropicBulkShearState, MaterialStatePointError> {
    ResolvedIsotropicBulkShearState::try_from_material(
        resolve_material_state_point(card, point, &requirements.properties, selection)?,
        requirements,
        maximum_absolute_mandel_strain,
    )
}

fn require_scalar(
    requirement: &ScalarPropertyRequirement,
    dims: Dims,
) -> Result<(), MaterialStatePointError> {
    if requirement.quantity() != QuantitySpec::dimensional(dims)
        || requirement.hardness_test().is_some()
        || requirement.elastic_component().is_some()
        || requirement.strain_component().is_some()
        || requirement.stress_component().is_some()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: requirement.name().into(),
            reason: "elastic scalar parameters require the declared scalar quantity, not an addressed tensor or stress component",
        });
    }
    Ok(())
}
