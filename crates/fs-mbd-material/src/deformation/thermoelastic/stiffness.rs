//! One material model selection for thermal shape, stress, inertia and vibration.
use super::*;
use fs_material::elastic::OrthotropicElastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES, POISSON_RATIO_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::{Dims, Pressure};
use fs_solid::linear3::TetElasticMaterial;

/// Source-property schema and material axes for the existing solid owner.
/// Both rungs use the separately sourced isotropic expansion path. Orthotropic
/// stiffness does not imply, or approximate, an anisotropic expansion law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalStiffnessModel {
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

impl ThermalStiffnessModel {
    pub(super) fn resolve(
        self,
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
        let mut requirements = Vec::with_capacity(9);
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
        // Closed matter uses reference density. Asking the ordinary elastic
        // state resolver for rho(T) would add an unrelated current-density
        // requirement and could silently replace the Lagrangian mass measure.
        resolve_material_state_point(card, point, &requirements, selection)
    }

    pub(super) fn material(
        self,
        current: &ResolvedMaterialStatePoint,
        reference: &ResolvedMaterialStatePoint,
        identity: ContentHash,
        strain_limit: f64,
    ) -> Result<TetElasticMaterial, ThermalBodyError> {
        let density = reference
            .property(DENSITY_PROPERTY)
            .expect("reference density resolved")
            .value_si();
        let value = |key| {
            current
                .property(key)
                .expect("elastic requirement resolved")
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
