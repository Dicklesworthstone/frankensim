//! Fixed-condition symmetric magnetic permeability in the physical mesh frame.
//!
//! B = mu H, mu = Q diag(mu_i) Q^T, and the reversible density is H.B/2.
//! Reluctivity uses Q diag(1/mu_i) Q^T, never componentwise reciprocals of the
//! rotated tensor. Principal axes and source values are retained so inversion
//! does not subtract nearly equal entries of a high-contrast tensor.
//! Reference: COMSOL, *Magnetostatics, Theory*, "Linear Magnetic Materials":
//! <https://www.comsol.com/multiphysics/magnetostatics-theory>.
//!
//! Every principal value is an absolute H/m source at fixed external conditions.
//! This lossless reciprocal law supplies no saturation, magnetic history,
//! gyrotropy, frequency dispersion, thermal parameter work or mechanical stress.

use super::super::{
    FieldStorageError, MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, coefficient,
    positive,
};
use crate::state_point::ResolvedMaterialStatePoint;

/// Material-axis order 1, 2, 3, each in absolute henries per metre.
pub const PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES: [&str; 3] = [
    "magnetic_permeability_1",
    "magnetic_permeability_2",
    "magnetic_permeability_3",
];

/// Validated principal directions as rows in physical mesh coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialMagneticFrame {
    axes: [[f64; 3]; 3],
}
impl SpatialMagneticFrame {
    /// Original principal directions in material-axis order.
    pub const fn principal_axes(self) -> [[f64; 3]; 3] {
        self.axes
    }
}

/// Magnetic orientation is independent of electrical and dielectric axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpatialMagneticLaw {
    /// One absolute scalar source for all directions.
    Isotropic,
    /// Three original positive principal sources in a fixed material frame.
    Orthotropic(SpatialMagneticFrame),
}

/// Prepared numerical response; original source authority stays with the caller.
/// The retained spectral representation supplies the inverse without a fresh
/// matrix factorization, while B(H) uses the cached matrix or scalar fast path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialPermeability {
    principal_h_m: [f64; 3],
    frame: Option<SpatialMagneticFrame>,
    matrix_h_m: [[f64; 3]; 3],
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn rotated_positive(
    principal: [f64; 3],
    frame: SpatialMagneticFrame,
) -> Result<[[f64; 3]; 3], FieldStorageError> {
    // tensor::rotate stores principal directions as columns.
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
    // Positive source eigenvalues do not guarantee that a rounded Cartesian
    // matrix still represents them. Refuse lost definiteness; FEEC separately
    // retains its interval-SPD admission on the matrix actually assembled.
    let scale = principal.into_iter().fold(0.0_f64, f64::max);
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
            "rotated magnetic tensor is not representably positive",
        ));
    }
    Ok(matrix)
}

impl SpatialPermeability {
    /// Checked numerical scalar, for owners of an already admitted recoil law.
    /// This constructs no material source, property receipt or physical claim.
    pub fn isotropic(value_h_m: f64) -> Result<Self, FieldStorageError> {
        let value = positive(value_h_m, "magnetic permeability must be finite positive")?;
        Ok(Self {
            principal_h_m: [value; 3],
            frame: None,
            matrix_h_m: [[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]],
        })
    }

    /// Scalar fast path only when the original law declared isotropy.
    pub const fn isotropic_h_m(self) -> Option<f64> {
        if self.frame.is_none() {
            Some(self.principal_h_m[0])
        } else {
            None
        }
    }

    /// Complete symmetric coefficient in physical mesh coordinates [H/m].
    pub const fn matrix_h_m(self) -> [[f64; 3]; 3] {
        self.matrix_h_m
    }

    /// Largest source eigenvalue, e.g. for a conservative magnetic Peclet bound.
    pub fn maximum_permeability_h_m(self) -> f64 {
        self.principal_h_m.into_iter().fold(0.0_f64, f64::max)
    }

    /// Full B = mu H [T]. Field owners check representability of the result.
    pub fn flux_density_t(self, field_a_m: [f64; 3]) -> [f64; 3] {
        if let Some(mu) = self.isotropic_h_m() {
            field_a_m.map(|value| mu * value)
        } else {
            self.matrix_h_m.map(|row| dot(row, field_a_m))
        }
    }

    /// Full inverse response in metres per henry, prepared from principal data.
    /// A finite permeability need not have a representable reciprocal; inverse
    /// owners reject that condition without restricting scalar-potential owners.
    pub fn reluctivity_matrix_m_h(self) -> Result<[[f64; 3]; 3], FieldStorageError> {
        let reciprocal = self.principal_h_m.map(|value| value.recip());
        for value in reciprocal {
            positive(value, "magnetic reluctivity must be finite positive")?;
        }
        if let Some(frame) = self.frame {
            rotated_positive(reciprocal, frame)
        } else {
            let value = reciprocal[0];
            Ok([[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]])
        }
    }

    /// H = mu^-1 B [A/m]. B is the induced part for a separate recoil offset.
    /// The scalar path retains division; tensors project onto retained principal
    /// directions instead of numerically inverting the rounded Cartesian matrix.
    pub fn field_strength_a_m(
        self,
        flux_density_t: [f64; 3],
    ) -> Result<[f64; 3], FieldStorageError> {
        let field = if let Some(frame) = self.frame {
            let projected: [f64; 3] =
                std::array::from_fn(|i| dot(frame.axes[i], flux_density_t) / self.principal_h_m[i]);
            std::array::from_fn(|j| {
                frame.axes[0][j].mul_add(
                    projected[0],
                    frame.axes[1][j].mul_add(projected[1], frame.axes[2][j] * projected[2]),
                )
            })
        } else {
            flux_density_t.map(|value| value / self.principal_h_m[0])
        };
        if field.iter().any(|value| !value.is_finite()) {
            return Err(FieldStorageError::Invalid(
                "magnetic inverse response is unrepresentable",
            ));
        }
        Ok(field)
    }

    /// (nu_next - nu_self) B for the same fixed principal frame [A/m].
    /// This is a numerical parameter-work contraction, not source or thermal
    /// admission. Principal differences avoid subtracting rounded inverse
    /// matrices; an unchanged coefficient gives exactly zero in every mode.
    pub fn reluctivity_change_action(
        self,
        next: Self,
        flux_density_t: [f64; 3],
    ) -> Result<[f64; 3], FieldStorageError> {
        if self.frame != next.frame || flux_density_t.iter().any(|value| !value.is_finite()) {
            return Err(FieldStorageError::Invalid(
                "magnetic parameter work requires the same fixed frame and finite flux",
            ));
        }
        let projected = match self.frame {
            Some(frame) => frame.axes.map(|axis| dot(axis, flux_density_t)),
            None => flux_density_t,
        };
        let mut change = [0.0; 3];
        for i in 0..3 {
            let (before, after) = (self.principal_h_m[i], next.principal_h_m[i]);
            if before == after || projected[i] == 0.0 {
                continue;
            }
            let relative = (before - after) / before.max(after);
            let field = projected[i] / before.min(after);
            // |relative| <= 1. If division overflows, multiplying first may
            // still recover a representable action for nearly equal sources.
            change[i] = if field.is_finite() {
                relative * field
            } else {
                (relative * projected[i]) / before.min(after)
            };
            if !change[i].is_finite() || change[i] == 0.0 {
                return Err(FieldStorageError::Invalid(
                    "magnetic reluctivity change is unrepresentable",
                ));
            }
        }
        let action = match self.frame {
            Some(frame) => std::array::from_fn(|j| {
                frame.axes[0][j].mul_add(
                    change[0],
                    frame.axes[1][j].mul_add(change[1], frame.axes[2][j] * change[2]),
                )
            }),
            None => change,
        };
        if action.iter().any(|value| !value.is_finite()) {
            return Err(FieldStorageError::Invalid(
                "magnetic reluctivity-change action is unrepresentable",
            ));
        }
        Ok(action)
    }
}

impl SpatialMagneticLaw {
    /// Each row is a principal direction in the destination mesh frame.
    /// Finite orthonormal input is required; no silent orientation repair occurs.
    pub fn orthotropic(axes: [[f64; 3]; 3]) -> Result<Self, FieldStorageError> {
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                let actual = dot(axes[i], axes[j]);
                if !actual.is_finite() || (actual - expected).abs() > 1e-12 {
                    return Err(FieldStorageError::Invalid(
                        "spatial magnetic axes must be finite and orthonormal",
                    ));
                }
            }
        }
        Ok(Self::Orthotropic(SpatialMagneticFrame { axes }))
    }

    /// Fixed principal frame, or no preferred directions for an isotropic law.
    pub const fn principal_axes(self) -> Option<[[f64; 3]; 3]> {
        match self {
            Self::Isotropic => None,
            Self::Orthotropic(frame) => Some(frame.axes),
        }
    }

    /// An isotropic property never fills a missing principal source.
    pub const fn source_property_names(self) -> &'static [&'static str] {
        match self {
            Self::Isotropic => &[MAGNETIC_PERMEABILITY_PROPERTY],
            Self::Orthotropic(_) => &PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES,
        }
    }

    /// Resolve absolute H/m from the complete original material bundle.
    /// Every axis uses the shared fixed-condition storage admission, including
    /// refusal of field, frequency, strain and stress dependence or bounds.
    pub fn resolve(
        self,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<SpatialPermeability, FieldStorageError> {
        match self {
            Self::Isotropic => SpatialPermeability::isotropic(coefficient(
                material,
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
            )?),
            Self::Orthotropic(frame) => {
                let principal_h_m = [
                    coefficient(
                        material,
                        PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES[0],
                        MAGNETIC_PERMEABILITY_DIMS,
                    )?,
                    coefficient(
                        material,
                        PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES[1],
                        MAGNETIC_PERMEABILITY_DIMS,
                    )?,
                    coefficient(
                        material,
                        PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES[2],
                        MAGNETIC_PERMEABILITY_DIMS,
                    )?,
                ];
                Ok(SpatialPermeability {
                    matrix_h_m: rotated_positive(principal_h_m, frame)?,
                    principal_h_m,
                    frame: Some(frame),
                })
            }
        }
    }
}
