//! One material model selection for shape, stress, inertia and vibration.
use super::SolidMaterialError;
use crate::linear3::TetElasticMaterial;
use fs_blake3::ContentHash;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::elastic::OrthotropicElastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};

/// Source-property schema and material axes for the existing solid owner.
/// Thermal expansion independently selects its own source law and fixed axes.
/// Orthotropic stiffness does not imply an anisotropic expansion coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SolidStiffnessModel {
    /// Scalar E and nu, with the original isotropic tensor construction.
    Isotropic,
    /// Nine principal-axis engineering constants, rotated once into mesh axes.
    OrientedOrthotropic {
        /// Proper orthonormal rotation; columns are material axes in mesh coordinates.
        principal_to_mesh: [[f64; 3]; 3],
        /// Nonzero identity of the authored material-axis assignment.
        orientation_identity: ContentHash,
    },
}

/// Thermal consumer name for the shared solid material schema.
pub type ThermalStiffnessModel = SolidStiffnessModel;

impl SolidStiffnessModel {
    fn requirements(self) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
        let mut requirements = Vec::with_capacity(10);
        match self {
            Self::Isotropic => {
                requirements.push(ScalarPropertyRequirement::try_new(
                    YOUNG_MODULUS_PROPERTY,
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )?);
                requirements.push(ScalarPropertyRequirement::try_new(
                    POISSON_RATIO_PROPERTY,
                    Dims::NONE,
                    ScalarAdmissibility::OpenInterval {
                        lower: -1.0,
                        upper: 0.5,
                    },
                )?);
            }
            Self::OrientedOrthotropic { .. } => {
                for property in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES
                    .into_iter()
                    .chain(ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES)
                {
                    requirements.push(ScalarPropertyRequirement::try_new(
                        property,
                        Pressure::DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )?);
                }
                for property in ORTHOTROPIC_POISSON_RATIO_PROPERTIES {
                    requirements.push(ScalarPropertyRequirement::try_new(
                        property,
                        Dims::NONE,
                        ScalarAdmissibility::Finite,
                    )?);
                }
            }
        }
        Ok(requirements)
    }

    /// Resolve the complete elastic law at one current state, without requiring
    /// current density. The original reference inventory owns inertia.
    pub fn resolve(
        self,
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
        let requirements = self.requirements()?;
        // Closed matter uses reference density. Asking the ordinary elastic
        // state resolver for rho(T) would add an unrelated current-density
        // requirement and could silently replace the Lagrangian mass measure.
        resolve_material_state_point(card, point, &requirements, selection)
    }

    /// Resolve density and the chosen elastic schema together at the original
    /// reference state. Callers pin these claims before thermal evolution.
    pub fn resolve_reference(
        self,
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
        let mut requirements = self.requirements()?;
        requirements.push(ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        resolve_material_state_point(card, point, &requirements, selection)
    }

    /// Lower current stiffness into the original solid tensor owner.
    ///
    /// Density is the conserved reference mass divided by reference volume,
    /// never a fresh rho(T) substitution. The caller's nonzero identity binds
    /// that inventory and the current source state; directional axes are then
    /// bound by the solid owner's existing orthotropic admission.
    pub fn material(
        self,
        current: &ResolvedMaterialStatePoint,
        reference_density_kg_m3: f64,
        identity: ContentHash,
        strain_limit: f64,
    ) -> Result<TetElasticMaterial, SolidMaterialError> {
        // This API is shared and public: incomplete or mismatched source
        // receipts must refuse rather than panic on a missing engineering constant.
        for requirement in self.requirements()? {
            let property = current.property(requirement.name()).ok_or(
                MaterialStatePointError::InvalidDerived {
                    quantity: "missing_elastic_source_property",
                },
            )?;
            if property.requirement() != &requirement {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "elastic_source_requirement",
                }
                .into());
            }
        }
        let density = reference_density_kg_m3;
        let value = |key| {
            current
                .property(key)
                .expect("complete matching elastic schema checked above")
                .value_si()
        };
        Ok(match self {
            Self::Isotropic => TetElasticMaterial::try_new(
                density,
                value(YOUNG_MODULUS_PROPERTY),
                value(POISSON_RATIO_PROPERTY),
                identity.0,
            )?,
            Self::OrientedOrthotropic {
                principal_to_mesh,
                orientation_identity,
            } => {
                let law = OrthotropicElastic::new(
                    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.map(value),
                    ORTHOTROPIC_POISSON_RATIO_PROPERTIES.map(value),
                    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.map(value),
                    strain_limit,
                )
                .map_err(|_| MaterialStatePointError::InvalidDerived {
                    quantity: "orthotropic_compliance",
                })?;
                TetElasticMaterial::try_new_oriented_orthotropic(
                    density,
                    &law,
                    principal_to_mesh,
                    identity.0,
                    orientation_identity.0,
                )?
            }
        })
    }
}

