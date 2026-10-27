//! Source-backed in-plane dielectric storage on explicit fixed material axes.
//!
//! Both principal values are absolute SI permittivities. The existing electrical
//! frame and second-order tensor rotation give epsilon=Q diag(e1,e2) Q^T.
//! No conductivity orientation, missing normal response or isotropy is inferred.
//! The physical plane must be decoupled from the omitted normal direction.
//! Static snapshots may use admitted external conditions. Retained charge
//! owners additionally require every selected coefficient to be source-declared
//! constant throughout its validity domain; changing axes/storage needs its
//! own mechanical/thermal parameter-work law.

use super::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, FieldStorageError, coefficient,
    constant_coefficient,
};
use crate::{conductor::planar::PlanarElectricalFrame, state_point::ResolvedMaterialStatePoint};

/// Material-axis order 1, 2, each in absolute farads per metre.
pub const PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES: [&str; 2] =
    ["electric_permittivity_1", "electric_permittivity_2"];

/// An explicit source-law choice, independent of the conductivity's basis/frame.
/// Frame correspondence to the physical mesh belongs to the field owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlanarDielectricLaw {
    Isotropic,
    Orthotropic(PlanarElectricalFrame),
}

/// Actual scalar or symmetric tensor response in mesh coordinates [F/m].
/// A bare coefficient is observational data, not an admissible material source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlanarPermittivity {
    Scalar(f64),
    Tensor([[f64; 2]; 2]),
}
impl PlanarPermittivity {
    pub fn matrix_f_m(self) -> [[f64; 2]; 2] {
        match self {
            Self::Scalar(value) => [[value, 0.0], [0.0, value]],
            Self::Tensor(matrix) => matrix,
        }
    }

    /// D=epsilon E in C/m2. Native owners check representability and integrate
    /// E.D/2 against the original physical cell volumes, never a scalar average.
    pub fn displacement_c_m2(self, field_v_m: [f64; 2]) -> [f64; 2] {
        let [x, y] = field_v_m;
        match self {
            Self::Scalar(value) => [value * x, value * y],
            Self::Tensor(matrix) => [
                matrix[0][0].mul_add(x, matrix[0][1] * y),
                matrix[1][0].mul_add(x, matrix[1][1] * y),
            ],
        }
    }
}

impl PlanarDielectricLaw {
    pub const fn frame(self) -> Option<PlanarElectricalFrame> {
        match self {
            Self::Isotropic => None,
            Self::Orthotropic(frame) => Some(frame),
        }
    }

    /// Every independent material coefficient, with no fallback between laws.
    pub const fn source_property_names(self) -> &'static [&'static str] {
        match self {
            Self::Isotropic => &[ELECTRIC_PERMITTIVITY_PROPERTY],
            Self::Orthotropic(_) => &PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES,
        }
    }

    /// Resolve a fixed-condition snapshot from the complete original bundle.
    /// The shared linear-storage gate retains quantity, source support and
    /// applicability checks, including rejection of unowned frequency/field axes.
    pub fn resolve(
        self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<PlanarPermittivity, FieldStorageError> {
        self.resolve_with(material, false)
    }

    /// Resolve storage suitable for a trajectory with no dielectric parameter
    /// work. Each principal source must explicitly declare validity-wide
    /// constancy; equal values at two temperatures do not establish that fact.
    pub fn resolve_constant(
        self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<PlanarPermittivity, FieldStorageError> {
        self.resolve_with(material, true)
    }

    fn resolve_with(
        self,
        material: &ResolvedMaterialStatePoint,
        constant: bool,
    ) -> Result<PlanarPermittivity, FieldStorageError> {
        let resolve = |name| {
            if constant {
                constant_coefficient(material, name, ELECTRIC_PERMITTIVITY_DIMS)
            } else {
                coefficient(material, name, ELECTRIC_PERMITTIVITY_DIMS)
            }
        };
        match self {
            Self::Isotropic => Ok(PlanarPermittivity::Scalar(resolve(
                ELECTRIC_PERMITTIVITY_PROPERTY,
            )?)),
            Self::Orthotropic(frame) => {
                let [first, second] = PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES;
                let principal = [resolve(first)?, resolve(second)?];
                let rotation = crate::tensor::rotation([0.0, 0.0, 1.0], frame.angle_rad());
                let rotated = crate::tensor::rotate(
                    &[principal[0], principal[1], 0.0, 0.0, 0.0, 0.0],
                    &rotation,
                );
                let matrix = [[rotated[0], rotated[3]], [rotated[3], rotated[1]]];
                let scale = principal[0].max(principal[1]);
                let determinant = (matrix[0][0] / scale).mul_add(
                    matrix[1][1] / scale,
                    -(matrix[0][1] / scale) * (matrix[1][0] / scale),
                );
                if matrix.iter().flatten().any(|value| !value.is_finite())
                    || !determinant.is_finite()
                    || determinant <= 0.0
                {
                    return Err(FieldStorageError::Invalid(
                        "rotated dielectric tensor is not representably positive",
                    ));
                }
                Ok(PlanarPermittivity::Tensor(matrix))
            }
        }
    }
}
