//! Source-backed three-dimensional linear dielectric response.
//!
//! Three absolute SI principal permittivities and fixed orthonormal axes give
//! epsilon = Q diag(epsilon_i) Q^T in mesh coordinates. Every direction has its
//! own source; no conductivity axis, missing normal response or relative
//! permittivity is inferred. This is a fixed-condition lossless snapshot, not
//! a dispersive, nonlinear, moving-material or thermal-work law.

use super::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, FieldStorageError, coefficient,
};
use crate::state_point::ResolvedMaterialStatePoint;

/// Material-axis order 1, 2, 3, each in absolute farads per metre.
pub const PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES: [&str; 3] = [
    "electric_permittivity_1",
    "electric_permittivity_2",
    "electric_permittivity_3",
];

/// Validated principal directions, expressed as rows in physical mesh coordinates.
/// The signs of principal axes are immaterial for a second-order dielectric tensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialDielectricFrame {
    axes: [[f64; 3]; 3],
}
impl SpatialDielectricFrame {
    /// Original orthonormal rows in physical mesh coordinates.
    pub const fn principal_axes(self) -> [[f64; 3]; 3] {
        self.axes
    }
}

/// Explicit source-law choice, independent of the conductor's orientation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpatialDielectricLaw {
    /// One absolute scalar source, shared by all directions.
    Isotropic,
    /// Three principal sources on an admitted, fixed material frame.
    Orthotropic(SpatialDielectricFrame),
}

/// Actual scalar or symmetric tensor response in mesh coordinates [F/m].
/// A numerical coefficient is observational data, not an admissible source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpatialPermittivity {
    /// Isotropic absolute permittivity [F/m].
    Scalar(f64),
    /// Complete symmetric tensor in physical mesh coordinates [F/m].
    Tensor([[f64; 3]; 3]),
}
impl SpatialPermittivity {
    /// Complete coefficient matrix, including every off-diagonal response.
    pub fn matrix_f_m(self) -> [[f64; 3]; 3] {
        match self {
            Self::Scalar(value) => [[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]],
            Self::Tensor(matrix) => matrix,
        }
    }

    /// D = epsilon E [C/m2]. Field owners check finite results and integrate
    /// E.D/2 on the original physical volumes, without a scalar surrogate.
    pub fn displacement_c_m2(self, field_v_m: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Scalar(value) => field_v_m.map(|e| value * e),
            Self::Tensor(matrix) => matrix.map(|row| {
                row[0].mul_add(
                    field_v_m[0],
                    row[1].mul_add(field_v_m[1], row[2] * field_v_m[2]),
                )
            }),
        }
    }
}

impl SpatialDielectricLaw {
    /// Each row is one principal direction in the destination mesh frame.
    /// Axes must already be finite and orthonormal; no repair or normalization
    /// silently changes the caller's material orientation.
    pub fn orthotropic(axes: [[f64; 3]; 3]) -> Result<Self, FieldStorageError> {
        for i in 0..3 {
            for j in 0..3 {
                let dot = axes[i][0].mul_add(
                    axes[j][0],
                    axes[i][1].mul_add(axes[j][1], axes[i][2] * axes[j][2]),
                );
                let expected = if i == j { 1.0 } else { 0.0 };
                if !dot.is_finite() || (dot - expected).abs() > 1e-12 {
                    return Err(FieldStorageError::Invalid(
                        "spatial dielectric axes must be finite and orthonormal",
                    ));
                }
            }
        }
        Ok(Self::Orthotropic(SpatialDielectricFrame { axes }))
    }

    /// Explicit principal frame, or no preferred directions for an isotropic law.
    pub const fn principal_axes(self) -> Option<[[f64; 3]; 3]> {
        match self {
            Self::Isotropic => None,
            Self::Orthotropic(frame) => Some(frame.principal_axes()),
        }
    }

    /// Every required source property; an isotropic value never fills a
    /// missing principal component of an explicitly orthotropic law.
    pub const fn source_property_names(self) -> &'static [&'static str] {
        match self {
            Self::Isotropic => &[ELECTRIC_PERMITTIVITY_PROPERTY],
            Self::Orthotropic(_) => &PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES,
        }
    }

    /// Resolve absolute SI coefficients from the complete original bundle.
    /// The shared storage gate enforces source applicability and quantity
    /// semantics, including refusal of frequency, field, strain and stress axes.
    pub fn resolve(
        self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<SpatialPermittivity, FieldStorageError> {
        let resolve = |name| coefficient(material, name, ELECTRIC_PERMITTIVITY_DIMS);
        match self {
            Self::Isotropic => Ok(SpatialPermittivity::Scalar(resolve(
                ELECTRIC_PERMITTIVITY_PROPERTY,
            )?)),
            Self::Orthotropic(frame) => {
                let [first, second, third] = PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES;
                let principal = [resolve(first)?, resolve(second)?, resolve(third)?];
                // tensor::rotate expects the principal directions as columns.
                let rotation = std::array::from_fn(|i| std::array::from_fn(|j| frame.axes[j][i]));
                let rotated = crate::tensor::rotate(
                    &[principal[0], principal[1], principal[2], 0.0, 0.0, 0.0],
                    &rotation,
                );
                let matrix = [
                    [rotated[0], rotated[3], rotated[5]],
                    [rotated[3], rotated[1], rotated[4]],
                    [rotated[5], rotated[4], rotated[2]],
                ];
                // Positivity of sources is exact data admission; this separate
                // check refuses a tensor whose floating-point rotation lost it.
                // Weighted FEEC retains its own interval-SPD assembly gate.
                let scale = principal.into_iter().fold(0.0, f64::max);
                let n = matrix.map(|row| row.map(|value| value / scale));
                let minor = n[0][0].mul_add(n[1][1], -n[0][1] * n[0][1]);
                let determinant = n[0][0] * (n[1][1] * n[2][2] - n[1][2] * n[1][2])
                    - n[0][1] * (n[0][1] * n[2][2] - n[1][2] * n[0][2])
                    + n[0][2] * (n[0][1] * n[1][2] - n[1][1] * n[0][2]);
                if matrix.iter().flatten().any(|value| !value.is_finite())
                    || n[0][0] <= 0.0
                    || !minor.is_finite()
                    || minor <= 0.0
                    || !determinant.is_finite()
                    || determinant <= 0.0
                {
                    return Err(FieldStorageError::Invalid(
                        "rotated spatial dielectric tensor is not representably positive",
                    ));
                }
                Ok(SpatialPermittivity::Tensor(matrix))
            }
        }
    }
}
