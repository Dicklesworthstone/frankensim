//! Original cubic crystal coefficients shared by waves and solid mechanics.
//!
//! Choosing this resolver declares cubic symmetry and a crystallographic frame;
//! no material name, isotropic averaging or inferred specimen texture selects
//! the law. Four original scalar claims remain four claims. Symmetry-generated
//! matrix entries are algebraic consequences, not additional observations.

use super::{ElasticAcousticState, require_scalar};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Pressure};

/// Exact independent density, C11, C12 and C44 roles in a cubic crystal frame.
/// Source quantity, selection, context and uncertainty are never rewritten.
#[derive(Clone, Debug, PartialEq)]
pub struct CubicElasticRequirements {
    properties: [ScalarPropertyRequirement; 4],
}

impl CubicElasticRequirements {
    pub fn try_new(
        density: ScalarPropertyRequirement,
        c11: ScalarPropertyRequirement,
        c12: ScalarPropertyRequirement,
        c44: ScalarPropertyRequirement,
    ) -> Result<Self, MaterialStatePointError> {
        require_scalar(&density, Density::DIMS)?;
        for r in [&c11, &c12, &c44] {
            require_scalar(r, Pressure::DIMS)?;
        }
        let properties = [density, c11, c12, c44];
        for (i, property) in properties.iter().enumerate() {
            if properties[..i].iter().any(|p| p.name() == property.name()) {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: property.name().into(),
                    reason: "cubic density/C11/C12/C44 must select four distinct source keys",
                });
            }
        }
        Ok(Self { properties })
    }

    /// Canonical scalar keys. C12 is finite, not necessarily positive; the
    /// entire tensor must satisfy the original elastic stability admission.
    pub fn canonical() -> Result<Self, MaterialStatePointError> {
        let positive = |name, dims| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        };
        Self::try_new(
            positive("density", Density::DIMS)?,
            positive("cubic_c11", Pressure::DIMS)?,
            ScalarPropertyRequirement::try_new(
                "cubic_c12",
                Pressure::DIMS,
                ScalarAdmissibility::Finite,
            )?,
            positive("cubic_c44", Pressure::DIMS)?,
        )
    }

    pub const fn properties(&self) -> &[ScalarPropertyRequirement; 4] {
        &self.properties
    }
}

/// The complete original source bundle and its admitted crystal-frame law.
/// A solver rotation consumes this same tensor, density and source identity.
/// Source uncertainty is retained; derived wave speeds are nominal predictions.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedCubicElasticState {
    material: ResolvedMaterialStatePoint,
    requirements: CubicElasticRequirements,
    crystal_frame: ContentHash,
    acoustic: ElasticAcousticState,
}

impl ResolvedCubicElasticState {
    /// Adopt an already selected bundle without discarding extra source roles.
    /// `crystal_frame` identifies the orthonormal [100]/[010]/[001] axes; it
    /// cannot be zero. No world orientation or polycrystal average is inferred.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        requirements: &CubicElasticRequirements,
        crystal_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, MaterialStatePointError> {
        if crystal_frame == ContentHash([0; 32]) {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "nonzero cubic crystallographic frame identity",
            });
        }
        let value = |required: &ScalarPropertyRequirement| {
            material
                .property(required.name())
                .filter(|p| p.requirement() == required)
                .map(|p| p.value_si())
                .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                    property: required.name().into(),
                    reason: "cubic response requires the exact original scalar requirement",
                })
        };
        let [density, c11, c12, c44] = &requirements.properties;
        let values = [value(density)?, value(c11)?, value(c12)?, value(c44)?];
        let mut identity = DomainHasher::new("org.frankensim.fs-material.cubic-elastic.v1");
        identity.update(material.identity().as_bytes());
        identity.update(crystal_frame.as_bytes());
        for r in &requirements.properties {
            identity.update(&(r.name().len() as u64).to_le_bytes());
            identity.update(r.name().as_bytes());
        }
        identity.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
        let acoustic = ElasticAcousticState::from_cubic(
            values[0],
            values[1],
            values[2],
            values[3],
            identity.finalize(),
            maximum_absolute_mandel_strain,
        )
        .map_err(|_| MaterialStatePointError::InvalidDerived {
            quantity: "finite positive-density stable cubic stiffness and strain ceiling",
        })?;
        Ok(Self {
            material,
            requirements: requirements.clone(),
            crystal_frame,
            acoustic,
        })
    }

    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn requirements(&self) -> &CubicElasticRequirements {
        &self.requirements
    }
    pub const fn crystal_frame(&self) -> ContentHash {
        self.crystal_frame
    }
    pub const fn acoustic(&self) -> &ElasticAcousticState {
        &self.acoustic
    }
    pub const fn identity(&self) -> ContentHash {
        self.acoustic.material_state_identity()
    }

    /// Refresh the same complete source selection at an admitted state. Failure
    /// leaves the accepted bundle, physical law and frame identity untouched.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material(
            self.material.at_query_point(card, point)?,
            &self.requirements,
            self.crystal_frame,
            self.acoustic.maximum_absolute_mandel_strain(),
        )
    }
}

/// Resolve the four independent source roles at one original physical query.
pub fn resolve_cubic_elastic_state(
    card: &MaterialCard,
    point: &QueryPoint,
    requirements: &CubicElasticRequirements,
    crystal_frame: ContentHash,
    maximum_absolute_mandel_strain: f64,
    selection: MaterialPropertySelection,
) -> Result<ResolvedCubicElasticState, MaterialStatePointError> {
    ResolvedCubicElasticState::try_from_material(
        resolve_material_state_point(card, point, &requirements.properties, selection)?,
        requirements,
        crystal_frame,
        maximum_absolute_mandel_strain,
    )
}
