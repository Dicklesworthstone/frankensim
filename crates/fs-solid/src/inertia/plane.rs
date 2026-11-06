//! Consistent two-component mass in the original thin-solid configuration.
//!
//! Density is explicitly constant on each original element and thickness is
//! the original physical thickness. Neither current deformation, natural-volume
//! changes nor a plane-stress thickness update changes the closed-matter mass.
//! P1 uses the original three-point degree-two rule, not the constitutive
//! centroid. Q1 uses four Gauss points: shape products times the bilinear map's
//! affine Jacobian determinant are integrated exactly for constant cell density.

use super::{InertiaError, ReferenceMass, checkpoint, finite, invalid};
use crate::mesh2::{Mesh2, checked_reference_element, quad_points, shapes_at};
use fs_exec::Cx;
use fs_sparse::Coo;

impl ReferenceMass {
    /// Build immutable P1/Q1 reference mass with two components per mesh node.
    ///
    /// Exactly one positive original density is required per element. All
    /// original corner Jacobians must be positive, including Q1 corners outside
    /// the integration points. The explicit byte cap bounds conservative COO,
    /// sorting and CSR assembly storage before allocation. Original shapes and
    /// quadrature are accumulated in element/point order, with identical bits
    /// in each symmetric pair and each uncoupled spatial component.
    pub(crate) fn from_plane(
        mesh: &Mesh2,
        reference_thickness_m: f64,
        reference_density_kg_m3: &[f64],
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, InertiaError> {
        checkpoint(cx)?;
        if mesh.nodes.is_empty()
            || mesh.elems.is_empty()
            || reference_density_kg_m3.len() != mesh.elems.len()
            || !reference_thickness_m.is_finite()
            || reference_thickness_m <= 0.0
        {
            return Err(invalid(
                "nonempty plane mesh, original positive thickness and exact element density coverage required",
            ));
        }
        let ndofs = mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("plane reference mass displacement count overflows"))?;
        let mut entries = 0usize;
        for (element, conn) in mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            if !(conn.len() == 3 || conn.len() == 4) {
                return Err(invalid(
                    "plane reference mass requires P1 triangles or Q1 quads",
                ));
            }
            let rho = reference_density_kg_m3[element];
            if !rho.is_finite() || rho <= 0.0 {
                return Err(invalid(
                    "one finite positive original density per plane element required",
                ));
            }
            entries = entries
                .checked_add(2 * conn.len() * conn.len())
                .ok_or_else(|| invalid("plane reference mass entry count overflows"))?;
        }
        // Two independent component blocks. The same conservative per-entry
        // allowance as Tet10 includes COO growth, stable-sort scratch and CSR.
        // One fixed local tile, centered coordinates and native shape/rule
        // vectors are the only quadrature workspace.
        let estimated = entries
            .checked_mul(128)
            .and_then(|bytes| {
                ndofs
                    .checked_add(1)
                    .and_then(|rows| rows.checked_mul(core::mem::size_of::<usize>()))
                    .and_then(|rows| bytes.checked_add(rows))
            })
            .and_then(|bytes| bytes.checked_add(1024))
            .ok_or_else(|| invalid("plane reference mass assembly size overflows"))?;
        if estimated > maximum_assembly_bytes {
            return Err(invalid(
                "plane reference mass exceeds its assembly byte cap",
            ));
        }
        for (node, point) in mesh.nodes.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            if point.iter().any(|value| !value.is_finite()) {
                return Err(invalid("plane reference coordinates must be finite"));
            }
        }
        let mut matrix = Coo::new(ndofs, ndofs);
        let mut reference_mass_kg = 0.0;
        for (element, conn) in mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            let reference = checked_reference_element(&mesh.nodes, conn).ok_or_else(|| {
                invalid("plane reference element has invalid connectivity or orientation")
            })?;
            let arity = conn.len();
            let mut local = [[0.0; 4]; 4];
            for (xi, eta, weight) in quad_points(arity) {
                checkpoint(cx)?;
                let (shapes, gradients, jacobian) =
                    shapes_at(&reference, &[0, 1, 2, 3][..arity], xi, eta);
                let volume = finite(weight * jacobian * reference_thickness_m)?;
                if volume <= 0.0
                    || gradients.iter().flatten().any(|v| !v.is_finite())
                    || shapes
                        .iter()
                        .any(|n| !n.is_finite() || *n < 0.0 || *n > 1.0)
                    || (shapes.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
                {
                    return Err(invalid(
                        "positive plane reference quadrature and finite shapes required",
                    ));
                }
                let point_mass = finite(reference_density_kg_m3[element] * volume)?;
                if point_mass <= 0.0 {
                    return Err(invalid("plane reference point mass underflowed"));
                }
                reference_mass_kg = finite(reference_mass_kg + point_mass)?;
                for a in 0..arity {
                    for b in a..arity {
                        local[a][b] = (point_mass * shapes[a]).mul_add(shapes[b], local[a][b]);
                    }
                }
            }
            for a in 0..arity {
                for b in 0..arity {
                    let value = finite(local[a.min(b)][a.max(b)])?;
                    if value <= 0.0 {
                        return Err(invalid("consistent plane mass entry underflowed"));
                    }
                    for component in 0..2 {
                        matrix.push(2 * conn[a] + component, 2 * conn[b] + component, value);
                    }
                }
            }
        }
        Self::from_matrix(matrix.assemble(), reference_mass_kg, ndofs, cx)
    }
}

#[cfg(test)]
mod tests;
