//! Oriented finite thermal deformation from the original engineering-length paths.
//!
//! With H=Fth^-1 and Jth=det(Fth), the original-volume response is
//! W0=Jth We(F H), P0=Jth Pe H^T and
//! A0[iJ,kL]=Jth H[Ja] Ae[ia,kb] H[Lb]. The two tangent contractions
//! are separated, and the reference tensors are prepared once. No source
//! queries, interpolation, matrix inversions or allocation occur in mechanics.
//! This is the multiplicative thermal decomposition described in COMSOL's
//! Solid Mechanics theory, "Hyperelastic Materials / Thermal Expansion".

use super::{FiniteSolidError, is_fixed_external_coordinate};
use fs_blake3::ContentHash;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::thermal::transport::expansion::{
    FreeExpansionBasis, FreeExpansionGeometrySample, PreparedFreeExpansion,
};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Orientation of the source's stress-free expansion relative to the mesh.
/// This describes the thermal reference, independently of elastic symmetry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ExpansionThermalFrame {
    /// Allowed only for a path explicitly prepared as isotropic.
    Isotropic,
    /// Columns are the fixed principal directions in the named mesh frame.
    FixedPrincipal {
        /// Must equal the original directional source's frame identity.
        principal_frame: ContentHash,
        /// Nonzero identity of the mesh reference-coordinate frame.
        mesh_frame: ContentHash,
        /// Proper rotation with principal directions as columns.
        principal_to_mesh: [[f64; 3]; 3],
    },
}

/// Held finite thermal reference retaining its complete source path and frame.
///
/// Engineering length changes are normalized by their ORIGINAL datum through
/// `PreparedFreeExpansion`; they are not small-strain CTE integrals. The free
/// lengths define the relaxed constitutive reference even when supports prevent
/// actual expansion. Current geometry and density remain owned by the total
/// displacement map, not by the stress-free volume ratio.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpansionThermalReference {
    source: PreparedFreeExpansion,
    geometry: FreeExpansionGeometrySample,
    frame: ExpansionThermalFrame,
    reference_context: QueryPoint,
    current_context: QueryPoint,
    thermal_gradient: [f64; 9],
    inverse: [[f64; 3]; 3],
    cofactor: [[f64; 3]; 3],
}

impl ExpansionThermalReference {
    /// Prepare one endpoint on the original, source-admitted temperature path.
    /// The source card, directional frame and every fixed physical condition
    /// must match. Stress, strain, history and other evolving axes are refused.
    /// Negative expansion is allowed; nonpositive lengths already refuse in
    /// the original path owner. The rotation must be finite, proper and
    /// orthonormal within 1e-12; no orthogonalization or frame inference occurs.
    pub fn try_new(
        card: &MaterialCard,
        source: PreparedFreeExpansion,
        temperature_k: f64,
        frame: ExpansionThermalFrame,
    ) -> Result<Self, FiniteSolidError> {
        if source.material_card_identity() != card.content_hash() {
            return Err(FiniteSolidError::Invalid(
                "finite expansion requires its original material card",
            ));
        }
        let q = match (source.basis(), frame) {
            (FreeExpansionBasis::Isotropic, ExpansionThermalFrame::Isotropic) => IDENTITY,
            (
                FreeExpansionBasis::FixedPrincipal { frame: expected },
                ExpansionThermalFrame::FixedPrincipal {
                    principal_frame,
                    mesh_frame,
                    principal_to_mesh,
                },
            ) if expected == principal_frame && mesh_frame != ContentHash([0; 32]) => {
                admit_rotation(principal_to_mesh)?;
                if principal_frame == mesh_frame
                    && principal_to_mesh
                        .iter()
                        .flatten()
                        .zip(IDENTITY.iter().flatten())
                        .any(|(a, b)| (a - b).abs() > 1.0e-12)
                {
                    return Err(FiniteSolidError::Invalid(
                        "identical thermal source and mesh frames require the identity map",
                    ));
                }
                principal_to_mesh
            }
            _ => {
                return Err(FiniteSolidError::Invalid(
                    "finite expansion requires the original source basis and named principal/mesh frames",
                ));
            }
        };
        let axis = source.temperature_axis();
        let reference = source.reference_query_point();
        if axis != "T" && reference.axes().contains_key("T") {
            return Err(FiniteSolidError::Invalid(
                "finite expansion has an ambiguous extra canonical temperature coordinate",
            ));
        }
        let admit_axis = |name: &str| -> Result<(), FiniteSolidError> {
            if name != axis
                && !is_fixed_external_coordinate(
                    name,
                    reference.axis_quantities().get(name).copied(),
                    reference.axes().get(name).copied(),
                )
            {
                return Err(FiniteSolidError::Source(format!(
                    "finite expansion axis {name} requires an evolving physical owner",
                )));
            }
            Ok(())
        };
        for name in reference.axes().keys() {
            admit_axis(name)?;
        }
        for line in source.source_lines() {
            let selected = line.knots[0].1.receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(FiniteSolidError::Invalid(
                    "finite expansion lost its original selected source",
                ))?;
            for name in claim.validity.bounds().keys() {
                admit_axis(name)?;
            }
            for (name, _) in claim.value.axes() {
                admit_axis(name)?;
            }
        }
        let geometry = source
            .at_temperature(temperature_k)
            .map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        let context = |temperature: f64| -> Result<QueryPoint, FiniteSolidError> {
            let mut point = QueryPoint::new();
            for (name, value) in reference.axes() {
                if name == axis {
                    continue;
                }
                point = match reference.axis_quantities().get(name) {
                    Some(quantity) => point.with_quantity(name, *quantity, *value),
                    None => point.with(name, *value),
                }
                .map_err(|error| FiniteSolidError::Source(error.to_string()))?;
            }
            point
                .with_quantity("T", reference.axis_quantities()[axis], temperature)
                .map_err(|error| FiniteSolidError::Source(error.to_string()))
        };
        let reference_context = context(reference.axes()[axis])?;
        let current_context = context(temperature_k)?;
        let s = geometry.principal_stretches();
        let reciprocal = s.map(|value| 1.0 / value);
        let areas = [s[1] * s[2], s[0] * s[2], s[0] * s[1]];
        if reciprocal
            .into_iter()
            .chain(areas)
            .any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(FiniteSolidError::Invalid(
                "finite expansion inverse or area ratio is unrepresentable",
            ));
        }
        let rotate = |principal: [f64; 3]| -> [[f64; 3]; 3] {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| (0..3).map(|a| (q[i][a] * principal[a]) * q[j][a]).sum())
            })
        };
        let gradient = rotate(s);
        let inverse = rotate(reciprocal);
        let cofactor = rotate(areas);
        if gradient
            .iter()
            .chain(&inverse)
            .chain(&cofactor)
            .flatten()
            .any(|x| !x.is_finite())
        {
            return Err(FiniteSolidError::Invalid(
                "oriented thermal reference is unrepresentable",
            ));
        }
        let thermal_gradient = std::array::from_fn(|i| gradient[i / 3][i % 3]);
        Ok(Self {
            source,
            geometry,
            frame,
            reference_context,
            current_context,
            thermal_gradient,
            inverse,
            cofactor,
        })
    }

    /// Original engineering length curves, selections, reference and fixed context.
    pub const fn source(&self) -> &PreparedFreeExpansion {
        &self.source
    }
    /// Original normalized principal lengths and stress-free volume ratio.
    pub const fn geometry(&self) -> &FreeExpansionGeometrySample {
        &self.geometry
    }
    /// Explicit source-to-mesh orientation used by all mechanical operations.
    pub const fn frame(&self) -> ExpansionThermalFrame {
        self.frame
    }
    /// Original reference-density context with the temperature normalized to T.
    pub const fn reference_context(&self) -> &QueryPoint {
        &self.reference_context
    }
    /// Current held-temperature context of the underlying elastic law.
    pub const fn current_context(&self) -> &QueryPoint {
        &self.current_context
    }
    /// Row-major stress-free Fth in the mesh frame; not the actual constrained F.
    pub const fn thermal_gradient(&self) -> &[f64; 9] {
        &self.thermal_gradient
    }

    pub(super) fn elastic_deformation(&self, f: &[f64; 9]) -> [f64; 9] {
        if self.frame == ExpansionThermalFrame::Isotropic {
            return f.map(|value| value / self.geometry.principal_stretches()[0]);
        }
        std::array::from_fn(|index| {
            let (i, j) = (index / 3, index % 3);
            (0..3).map(|a| f[3 * i + a] * self.inverse[a][j]).sum()
        })
    }

    pub(super) fn pullback_piola(&self, elastic: [f64; 9]) -> [f64; 9] {
        if self.frame == ExpansionThermalFrame::Isotropic {
            let s = self.geometry.principal_stretches()[0];
            return elastic.map(|value| (s * s) * value);
        }
        std::array::from_fn(|index| {
            let (i, j) = (index / 3, index % 3);
            (0..3)
                .map(|a| elastic[3 * i + a] * self.cofactor[j][a])
                .sum()
        })
    }

    pub(super) fn pullback_tangent(&self, elastic: [[f64; 9]; 9]) -> [[f64; 9]; 9] {
        if self.frame == ExpansionThermalFrame::Isotropic {
            let s = self.geometry.principal_stretches()[0];
            return elastic.map(|row| row.map(|value| s * value));
        }
        let right: [[f64; 9]; 9] = std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                let (k, l) = (column / 3, column % 3);
                (0..3)
                    .map(|b| elastic[row][3 * k + b] * self.inverse[l][b])
                    .sum()
            })
        });
        std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                let (i, j) = (row / 3, row % 3);
                (0..3)
                    .map(|a| self.cofactor[j][a] * right[3 * i + a][column])
                    .sum()
            })
        })
    }
}

fn admit_rotation(q: [[f64; 3]; 3]) -> Result<(), FiniteSolidError> {
    if q.iter().flatten().any(|x| !x.is_finite()) {
        return Err(FiniteSolidError::Invalid(
            "finite thermal material axes required",
        ));
    }
    for i in 0..3 {
        for j in 0..3 {
            let dot: f64 = (0..3).map(|k| q[k][i] * q[k][j]).sum();
            if (dot - if i == j { 1.0 } else { 0.0 }).abs() > 1.0e-12 {
                return Err(FiniteSolidError::Invalid(
                    "orthonormal thermal material axes required",
                ));
            }
        }
    }
    let determinant = q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
        - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
        + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
    if (determinant - 1.0).abs() > 1.0e-12 {
        return Err(FiniteSolidError::Invalid(
            "proper thermal material rotation required",
        ));
    }
    Ok(())
}
