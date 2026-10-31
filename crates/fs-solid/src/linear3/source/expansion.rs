//! Source-integrated expansion and its one accepted mesh-frame strain.
//!
//! Elastic and thermal symmetry are independent explicit material choices.
//! Principal expansion axes stay fixed over the complete temperature path.

use super::SolidMaterialError;
use crate::linear3::TetThermalStrainState;
use crate::linear3::{
    ElasticTensorOrder, StrainTensorBasis, StrainTensorNotation, TetElasticMaterial,
    ThermalStrainTransformReceipt,
};
use fs_blake3::ContentHash;
use fs_matdb::{ClaimSelection, MaterialCard, QueryPoint};
use fs_material::state_point::{
    IntegratedIsotropicThermalExpansion, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    affine_temperature_slope, integrate_isotropic_thermal_expansion,
};
use fs_material::state_point::{
    IntegratedOrthotropicThermalExpansion, ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES,
    integrate_orthotropic_thermal_expansion,
};
use fs_material::visco_tensor::TensorThermalExpansion;

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

    /// Exact first and second free-strain derivatives for a globally affine
    /// expansion coefficient on the complete active temperature interval.
    ///
    /// Every knot remains pinned to the original integrated source. The current
    /// state anchors the derivative, so the stress-free reference may be outside
    /// this interval and need not share its affine extension. Principal rates
    /// use the same accepted mesh rotation as the integrated free strain.
    ///
    /// This does not replace source integration or admit a changed source/frame.
    /// Coefficient kinks, support gaps, unowned conditions, another card and
    /// unrepresentable derivative arithmetic refuse the complete result.
    pub fn affine_temperature_derivatives(
        &self,
        card: &MaterialCard,
        interval_k: [f64; 2],
    ) -> Result<TensorThermalExpansion, SolidMaterialError> {
        let invalid = |quantity| {
            SolidMaterialError::Material(MaterialStatePointError::InvalidDerived { quantity })
        };
        let current = self.current();
        let temperature_k = current
            .state_coordinate("T")
            .ok_or_else(|| invalid("expansion derivatives require the current T coordinate"))?;
        let [low, high] = interval_k;
        if current.card_identity() != card.content_hash()
            || self.reference().card_identity() != card.content_hash()
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || !(low..=high).contains(&temperature_k)
        {
            return Err(invalid(
                "expansion derivatives require the original card and a bounded current interval",
            ));
        }
        let lower = current.query().with("T", low)
            .map_err(|_| invalid("invalid expansion derivative lower query"))?;
        let upper = current.query().with("T", high)
            .map_err(|_| invalid("invalid expansion derivative upper query"))?;
        let coefficient = |name: &str| -> Result<(f64, f64), SolidMaterialError> {
            let property = current.property(name)
                .ok_or_else(|| invalid("original expansion coefficient is absent"))?;
            if !property.has_only_fixed_external_conditions(Some("T"), true) {
                return Err(invalid("expansion derivatives cannot freeze unowned source conditions"));
            }
            let selected = property.answer().receipt.selected;
            let claim = card.claims().claim(selected)
                .ok_or_else(|| invalid("selected expansion source claim is absent"))?;
            let path = card.claims().query_line_typed(
                &claim.key, "T", &lower, &upper, ClaimSelection::Pinned(selected),
            ).map_err(|source| MaterialStatePointError::Support {
                property: name.to_owned(), source: Box::new(source),
            })?;
            let knots: Vec<_> = path.knots.iter()
                .map(|(t, answer)| (*t, answer.evidence.value.value)).collect();
            let alpha = property.value_si();
            let slope = affine_temperature_slope(&knots, temperature_k, alpha)
                .map_err(invalid)?;
            Ok((alpha, slope))
        };
        let (alpha, slope) = match self {
            Self::Isotropic { .. } => {
                let (alpha, slope) = coefficient(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)?;
                ([alpha, alpha, alpha, 0.0, 0.0, 0.0],
                 [slope, slope, slope, 0.0, 0.0, 0.0])
            }
            Self::OrientedOrthotropic { transform, .. } => {
                let mut principal = [(0.0, 0.0); 3];
                for (value, name) in principal.iter_mut()
                    .zip(ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES) {
                    *value = coefficient(name)?;
                }
                let rotation = fs_material::tensor::mandel_rotation(*transform.source_to_target());
                let alpha = core::array::from_fn(|i| {
                    rotation[i][0] * principal[0].0
                        + rotation[i][1] * principal[1].0
                        + rotation[i][2] * principal[2].0
                });
                let slope = core::array::from_fn(|i| {
                    rotation[i][0] * principal[0].1
                        + rotation[i][1] * principal[1].1
                        + rotation[i][2] * principal[2].1
                });
                (alpha, slope)
            }
        };
        TensorThermalExpansion::try_new(temperature_k, alpha, slope, interval_k)
            .map_err(|_| invalid("source-integrated expansion derivatives are unrepresentable"))
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
    /// Integrate the selected source law from the original stress-free state
    /// and admit its one mesh-frame eigenstrain through the existing solid owner.
    /// The principal axes are fixed for the complete temperature path.
    #[allow(clippy::too_many_arguments)]
    pub fn resolve(
        self,
        card: &MaterialCard,
        reference: &QueryPoint,
        current: &QueryPoint,
        selection: MaterialPropertySelection,
        material: &TetElasticMaterial,
        mesh_frame: ContentHash,
        strain_limit: f64,
    ) -> Result<ThermalExpansionSource, SolidMaterialError> {
        match self {
            Self::Isotropic => {
                let path =
                    integrate_isotropic_thermal_expansion(card, reference, current, selection)?;
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
            Self::OrientedOrthotropic {
                principal_to_mesh,
                principal_frame,
            } => {
                let path =
                    integrate_orthotropic_thermal_expansion(card, reference, current, selection)?;
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
