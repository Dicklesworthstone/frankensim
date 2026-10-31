//! Three-dimensional source-bearing Ohmic response in physical mesh coordinates.
//!
//! Three principal sources and an explicit fixed orthonormal frame define
//! sigma = Q diag(sigma_i) Q^T. Resistivity sources are inverted in their
//! principal frame after evaluation. There is no Hall response, electrical
//! memory, inferred missing direction, evolving frame or phase-mixing law.

use super::{ConductorError, ElectricalPropertyBasis, ResolvedConductorState};
use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::PropertyKey;

/// Native coefficient in S/m. Scalar arithmetic remains distinct, so existing
/// isotropic assemblies retain their original multiplication order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpatialConductivity {
    Scalar(f64),
    Tensor([[f64; 3]; 3]),
}
impl SpatialConductivity {
    pub fn matrix_s_m(self) -> [[f64; 3]; 3] {
        match self {
            Self::Scalar(s) => [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]],
            Self::Tensor(s) => s,
        }
    }
    pub fn current_density_a_m2(self, electric_field_v_m: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Scalar(s) => electric_field_v_m.map(|v| v * s),
            Self::Tensor(s) => std::array::from_fn(|i| {
                s[i][0].mul_add(electric_field_v_m[0],
                    s[i][1].mul_add(electric_field_v_m[1], s[i][2] * electric_field_v_m[2]))
            }),
        }
    }
    /// Symmetric-matrix infinity norm bounds the largest conductivity
    /// eigenvalue, including rounding in the actual rotated coefficient.
    pub fn maximum_conductivity_s_m(self) -> f64 {
        match self {
            Self::Scalar(s) => s,
            Self::Tensor(s) => s.into_iter()
                .map(|row| row.into_iter().map(f64::abs).sum::<f64>())
                .fold(0.0, f64::max),
        }
    }
}
mod sealed {
    pub trait Sealed {}
}

/// Shared immutable source state for native 3D current/induction/heat owners.
/// A numerical tensor alone cannot implement this source-admission boundary.
pub trait SpatialConductorMaterial: sealed::Sealed + Clone + PartialEq {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn source_basis(&self) -> ElectricalPropertyBasis;
    fn source_keys(&self) -> impl Iterator<Item = &PropertyKey>;
    /// Rows are principal directions in the same physical frame as E, J and B.
    fn principal_axes(&self) -> Option<[[f64; 3]; 3]>;
    fn conductivity(&self) -> SpatialConductivity;
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError>;
    fn has_admitted_static_sources(&self) -> bool {
        self.source_keys().all(|key| {
            self.material().property_by_key(key)
                .is_some_and(|p| p.is_static_at_fixed_external_conditions())
        })
    }
    fn has_same_electrical_source(&self, other: &Self) -> bool {
        self.source_basis() == other.source_basis()
            && self.source_keys().eq(other.source_keys())
            && self.principal_axes() == other.principal_axes()
    }
}
impl sealed::Sealed for ResolvedConductorState {}
impl SpatialConductorMaterial for ResolvedConductorState {
    fn material(&self) -> &ResolvedMaterialStatePoint { self.material() }
    fn source_basis(&self) -> ElectricalPropertyBasis { self.source_basis() }
    fn source_keys(&self) -> impl Iterator<Item = &PropertyKey> {
        core::iter::once(self.source_key())
    }
    fn principal_axes(&self) -> Option<[[f64; 3]; 3]> { None }
    fn conductivity(&self) -> SpatialConductivity {
        SpatialConductivity::Scalar(self.resistivity_ohm_m().recip())
    }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        ResolvedConductorState::with_material(self, material)
    }
}

/// Complete original multiphysics bundle and three explicitly selected
/// principal electrical sources. Fixed axes are data, not material-name dispatch.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedOrthotropicConductor3 {
    material: ResolvedMaterialStatePoint,
    basis: ElectricalPropertyBasis,
    keys: [PropertyKey; 3],
    axes: [[f64; 3]; 3],
    principal_s_m: [f64; 3],
    coefficient: SpatialConductivity,
}
impl ResolvedOrthotropicConductor3 {
    /// Retain the exact keys and complete source state. Repeated keys are an
    /// explicit equal-principal-value declaration, never an implicit fallback.
    /// Every key must be a positive static scalar source in the chosen basis.
    /// Axes are finite orthonormal rows in the destination physical frame.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: ElectricalPropertyBasis,
        keys: [PropertyKey; 3],
        axes: [[f64; 3]; 3],
    ) -> Result<Self, ConductorError> {
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| axes[i][k] * axes[j][k]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                if !dot.is_finite() || (dot - expected).abs() > 1e-12 {
                    return Err(ConductorError::InvalidInput {
                        quantity: "finite orthonormal spatial electrical axes",
                    });
                }
            }
        }
        let mut principal_s_m = [0.0; 3];
        for (i, key) in keys.iter().enumerate() {
            let scalar = ResolvedConductorState::try_from_material_with_property(
                material.clone(), basis, key,
            )?;
            if !scalar.source_property().is_static_at_fixed_external_conditions() {
                return Err(ConductorError::InvalidInput {
                    quantity: "static spatial principal electrical source",
                });
            }
            principal_s_m[i] = scalar.resistivity_ohm_m().recip();
            if !principal_s_m[i].is_finite() || principal_s_m[i] <= 0.0 {
                return Err(ConductorError::InvalidDerived {
                    quantity: "positive representable spatial principal conductivity",
                });
            }
        }
        let q = std::array::from_fn(|i| std::array::from_fn(|j| axes[j][i]));
        let s = crate::tensor::rotate(
            &[principal_s_m[0], principal_s_m[1], principal_s_m[2], 0.0, 0.0, 0.0], &q,
        );
        let matrix = [[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]];
        let scale = principal_s_m.into_iter().fold(0.0, f64::max);
        let n = matrix.map(|row| row.map(|v| v / scale));
        let minor = n[0][0].mul_add(n[1][1], -n[0][1] * n[0][1]);
        let determinant = n[0][0] * (n[1][1] * n[2][2] - n[1][2] * n[1][2])
            - n[0][1] * (n[0][1] * n[2][2] - n[1][2] * n[0][2])
            + n[0][2] * (n[0][1] * n[1][2] - n[1][1] * n[0][2]);
        let coefficient = SpatialConductivity::Tensor(matrix);
        if matrix.iter().flatten().any(|v| !v.is_finite())
            || n[0][0] <= 0.0 || !minor.is_finite() || minor <= 0.0
            || !determinant.is_finite() || determinant <= 0.0
            || !coefficient.maximum_conductivity_s_m().is_finite()
        {
            return Err(ConductorError::InvalidDerived {
                quantity: "representably positive spatial conductivity tensor",
            });
        }
        Ok(Self { material, basis, keys, axes, principal_s_m, coefficient })
    }
    pub const fn principal_conductivity_s_m(&self) -> &[f64; 3] {
        &self.principal_s_m
    }
}
impl sealed::Sealed for ResolvedOrthotropicConductor3 {}
impl SpatialConductorMaterial for ResolvedOrthotropicConductor3 {
    fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }
    fn source_basis(&self) -> ElectricalPropertyBasis { self.basis }
    fn source_keys(&self) -> impl Iterator<Item = &PropertyKey> { self.keys.iter() }
    fn principal_axes(&self) -> Option<[[f64; 3]; 3]> { Some(self.axes) }
    fn conductivity(&self) -> SpatialConductivity { self.coefficient }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        Self::try_from_material(material, self.basis, self.keys.clone(), self.axes)
    }
}

/// Heterogeneous scalar/orthotropic cells share one native field and heat solve.
#[derive(Clone, Debug, PartialEq)]
pub enum ResolvedSpatialConductor {
    Isotropic(ResolvedConductorState),
    Orthotropic(ResolvedOrthotropicConductor3),
}
impl From<ResolvedConductorState> for ResolvedSpatialConductor {
    fn from(value: ResolvedConductorState) -> Self { Self::Isotropic(value) }
}
impl From<ResolvedOrthotropicConductor3> for ResolvedSpatialConductor {
    fn from(value: ResolvedOrthotropicConductor3) -> Self { Self::Orthotropic(value) }
}
impl sealed::Sealed for ResolvedSpatialConductor {}
impl SpatialConductorMaterial for ResolvedSpatialConductor {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        match self { Self::Isotropic(s) => s.material(), Self::Orthotropic(s) => s.material() }
    }
    fn source_basis(&self) -> ElectricalPropertyBasis {
        match self { Self::Isotropic(s) => s.source_basis(), Self::Orthotropic(s) => s.source_basis() }
    }
    fn source_keys(&self) -> impl Iterator<Item = &PropertyKey> {
        let keys = match self {
            Self::Isotropic(s) => [Some(s.source_key()), None, None],
            Self::Orthotropic(s) => s.keys.each_ref().map(Some),
        };
        keys.into_iter().flatten()
    }
    fn principal_axes(&self) -> Option<[[f64; 3]; 3]> {
        match self { Self::Isotropic(_) => None, Self::Orthotropic(s) => s.principal_axes() }
    }
    fn conductivity(&self) -> SpatialConductivity {
        match self {
            Self::Isotropic(s) => SpatialConductorMaterial::conductivity(s),
            Self::Orthotropic(s) => s.conductivity(),
        }
    }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        match self {
            Self::Isotropic(s) => s.with_material(material).map(Self::Isotropic),
            Self::Orthotropic(s) => s.with_material(material).map(Self::Orthotropic),
        }
    }
}
