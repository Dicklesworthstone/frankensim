//! Source-integrated expansion and its one accepted mesh-frame strain.
//!
//! Elastic and thermal symmetry are independent explicit material choices.
//! Principal expansion axes stay fixed over the complete temperature path.

use super::*;
use fs_material::state_point::{
    IntegratedOrthotropicThermalExpansion, integrate_orthotropic_thermal_expansion,
};
use fs_solid::linear3::{
    ElasticTensorOrder, StrainTensorBasis, StrainTensorNotation, TetElasticMaterial,
    ThermalStrainTransformReceipt,
};

/// Source-property schema and fixed material axes for thermal expansion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalExpansionModel {
    /// The original scalar coefficient and isotropic free strain.
    Isotropic,
    /// Three independent coefficients in fixed principal material directions.
    OrientedOrthotropic {
        /// Proper rotation; columns are principal axes in mesh coordinates.
        principal_to_mesh: [[f64; 3]; 3],
        /// Nonzero identity of the principal coordinate frame.
        principal_frame: ContentHash,
    },
}

/// Exact integrated source path and the strain accepted by the solid owner.
///
/// The orthotropic variant retains the complete rotation receipt. All later
/// mechanics consume this accepted strain; they never reconstruct an isotropic
/// approximation from a scalar coefficient.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalExpansionSource {
    /// Original scalar source and its unchanged isotropic lowering.
    Isotropic {
        /// Exact source-integrated scalar expansion.
        path: IntegratedIsotropicThermalExpansion,
        /// Accepted mesh-frame Mandel strain.
        state: TetThermalStrainState,
    },
    /// Independent principal paths transformed once into the mesh frame.
    OrientedOrthotropic {
        /// Exact three-axis source integration.
        path: IntegratedOrthotropicThermalExpansion,
        /// Original solid owner's complete frame-transform receipt.
        transform: ThermalStrainTransformReceipt,
    },
}

impl ThermalExpansionSource {
    /// Accepted state shared by displacement, stress and the vibration tangent.
    #[must_use]
    pub const fn thermal_strain(&self) -> &TetThermalStrainState {
        match self {
            Self::Isotropic { state, .. } => state,
            Self::OrientedOrthotropic { transform, .. } => transform.state(),
        }
    }

    /// Resolved source data at the original reference state.
    #[must_use]
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Isotropic { path, .. } => path.reference().resolved(),
            Self::OrientedOrthotropic { path, .. } => path.reference(),
        }
    }

    /// The same pinned source data at the accepted temperature.
    #[must_use]
    pub const fn current(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Isotropic { path, .. } => path.current().resolved(),
            Self::OrientedOrthotropic { path, .. } => path.current(),
        }
    }

    /// Signed integrated strains in the declared principal material directions.
    #[must_use]
    pub const fn free_principal_strains(&self) -> [f64; 3] {
        match self {
            Self::Isotropic { path, .. } => [path.free_linear_strain(); 3],
            Self::OrientedOrthotropic { path, .. } => *path.free_principal_strains(),
        }
    }

    /// Original integrated-law identity, before consumer orientation.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        match self {
            Self::Isotropic { path, .. } => path.identity(),
            Self::OrientedOrthotropic { path, .. } => path.identity(),
        }
    }

    /// Complete directional frame transformation when one was requested.
    #[must_use]
    pub const fn transform(&self) -> Option<&ThermalStrainTransformReceipt> {
        match self {
            Self::Isotropic { .. } => None,
            Self::OrientedOrthotropic { transform, .. } => Some(transform),
        }
    }
}

impl ThermalExpansionModel {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn resolve(
        self,
        card: &MaterialCard,
        reference: &QueryPoint,
        current: &QueryPoint,
        selection: MaterialPropertySelection,
        material: &TetElasticMaterial,
        mesh_frame: ContentHash,
        strain_limit: f64,
    ) -> Result<ThermalExpansionSource, ThermalBodyError> {
        match self {
            Self::Isotropic => {
                let path = integrate_isotropic_thermal_expansion(
                    card, reference, current, selection,
                )?;
                let strain = path.free_linear_strain();
                let state = TetThermalStrainState::try_new(
                    current.axes()["T"],
                    reference.axes()["T"],
                    [strain, strain, strain, 0.0, 0.0, 0.0],
                    material.material_state_identity,
                    path.identity(),
                    strain_limit,
                )?;
                Ok(ThermalExpansionSource::Isotropic { path, state })
            }
            Self::OrientedOrthotropic { principal_to_mesh, principal_frame } => {
                let path = integrate_orthotropic_thermal_expansion(
                    card, reference, current, selection,
                )?;
                let [first, second, third] = *path.free_principal_strains();
                let transform = TetThermalStrainState::try_from_material_and_oriented_strain(
                    material,
                    mesh_frame,
                    current.axes()["T"],
                    reference.axes()["T"],
                    [first, second, third, 0.0, 0.0, 0.0],
                    StrainTensorBasis {
                        notation: StrainTensorNotation::Tensor,
                        order: ElasticTensorOrder::XxYyZzXyYzZx,
                        frame: principal_frame,
                    },
                    principal_to_mesh,
                    path.identity(),
                    strain_limit,
                )?;
                Ok(ThermalExpansionSource::OrientedOrthotropic { path, transform })
            }
        }
    }
}
