//! Consistent two-component mass in the original thin-solid configuration.
//!
//! Density is either cellwise constant or an original P1/Q1 nodal field;
//! thickness is the original physical thickness. Neither current deformation,
//! natural-volume changes nor a plane-stress thickness update changes the mass.
//! P1 uses the original three-point degree-two rule, not the constitutive
//! centroid. Q1 uses four Gauss points: shape products times the bilinear map's
//! affine Jacobian determinant are integrated exactly for constant cell density.
//! Graded P1 uses exact cubic simplex moments; graded Q1 uses a 3-by-3 rule.

use super::{InertiaError, ReferenceMass, checkpoint, finite, invalid};
use crate::mesh2::{Mesh2, checked_reference_element, quad_points, shapes_at};
use fs_exec::Cx;
use fs_sparse::Coo;

#[derive(Clone, Copy)]
enum PlaneDensity<'a> {
    Elements(&'a [f64]),
    Nodes(&'a [f64]),
}

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
        Self::from_plane_density(
            mesh,
            reference_thickness_m,
            PlaneDensity::Elements(reference_density_kg_m3),
            maximum_assembly_bytes,
            cx,
        )
        .map(|(mass, _)| mass)
    }

    /// Exact consistent inertia for a positive density interpolated on the
    /// original mesh. Volume-weighted cell means are reporting data only;
    /// the consistent operator retains every nodal density contribution.
    pub(crate) fn from_plane_nodal_density(
        mesh: &Mesh2,
        reference_thickness_m: f64,
        reference_density_kg_m3: &[f64],
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<f64>), InertiaError> {
        Self::from_plane_density(
            mesh,
            reference_thickness_m,
            PlaneDensity::Nodes(reference_density_kg_m3),
            maximum_assembly_bytes,
            cx,
        )
    }

    fn from_plane_density(
        mesh: &Mesh2,
        reference_thickness_m: f64,
        density: PlaneDensity<'_>,
        maximum_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<(Self, Vec<f64>), InertiaError> {
        checkpoint(cx)?;
        let (values, count) = match density {
            PlaneDensity::Elements(values) => (values, mesh.elems.len()),
            PlaneDensity::Nodes(values) => (values, mesh.nodes.len()),
        };
        let nodal = matches!(density, PlaneDensity::Nodes(_));
        if mesh.nodes.is_empty()
            || mesh.elems.is_empty()
            || values.len() != count
            || !reference_thickness_m.is_finite()
            || reference_thickness_m <= 0.0
        {
            return Err(invalid(
                "nonempty plane mesh, original positive thickness and exact density coverage required",
            ));
        }
        let ndofs = mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("plane reference mass displacement count overflows"))?;
        let mut entries = 0usize;
        for conn in &mesh.elems {
            checkpoint(cx)?;
            if !(conn.len() == 3 || conn.len() == 4) {
                return Err(invalid(
                    "plane reference mass requires P1 triangles or Q1 quads",
                ));
            }
            entries = entries
                .checked_add(2 * conn.len() * conn.len())
                .ok_or_else(|| invalid("plane reference mass entry count overflows"))?;
        }
        for (i, &rho) in values.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            if !rho.is_finite() || rho <= 0.0 {
                return Err(invalid("finite positive original plane densities required"));
            }
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
            .and_then(|bytes| {
                let means = if nodal { mesh.elems.len() } else { 0 };
                means
                    .checked_mul(8)
                    .and_then(|extra| bytes.checked_add(extra))
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
        let mut means = Vec::with_capacity(if nodal { mesh.elems.len() } else { 0 });
        let triangle_rule = quad_points(3);
        let quad_rule = quad_points(4);
        let g = (3.0_f64 / 5.0).sqrt();
        let rule = [(-g, 5.0 / 9.0), (0.0, 8.0 / 9.0), (g, 5.0 / 9.0)];
        let graded_quad_rule: [(f64, f64, f64); 9] = std::array::from_fn(|i| {
            let (xi, wx) = rule[i % 3];
            let (eta, wy) = rule[i / 3];
            (xi, eta, wx * wy)
        });
        for (element, conn) in mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            let reference = checked_reference_element(&mesh.nodes, conn).ok_or_else(|| {
                invalid("plane reference element has invalid connectivity or orientation")
            })?;
            let arity = conn.len();
            let rho: [f64; 4] = match density {
                PlaneDensity::Elements(values) => [values[element]; 4],
                PlaneDensity::Nodes(values) => {
                    std::array::from_fn(|a| values[conn[a.min(arity - 1)]])
                }
            };
            let uniform = rho[..arity].iter().all(|r| r.to_bits() == rho[0].to_bits());
            let scale = rho.into_iter().fold(0.0_f64, f64::max);
            let normalized = if uniform {
                [0.0; 4]
            } else {
                rho.map(|r| r / scale)
            };
            let mut local = [[0.0; 4]; 4];
            // On P1, integral(rho*N_a*N_b) is cubic. The constitutive
            // centroid and the old edge-midpoint rule both underintegrate it.
            if arity == 3 && !uniform {
                let (_, _, jacobian) = shapes_at(&reference, &[0, 1, 2], 0.0, 0.0);
                let volume = finite(0.5 * jacobian * reference_thickness_m)?;
                let total = normalized[..3].iter().sum::<f64>();
                let mean = finite(scale * (total / 3.0))?;
                let mass = finite(volume * mean)?;
                if volume <= 0.0 || mass <= 0.0 {
                    return Err(invalid("graded plane reference mass underflowed"));
                }
                for a in 0..3 {
                    for b in a..3 {
                        // Diagonal V*(sum(rho)+2*rho_a)/30;
                        // off-diagonal V*(sum(rho)+rho_a+rho_b)/60.
                        let fraction = (normalized[a] + normalized[b]) / total;
                        local[a][b] = (mass / if a == b { 10.0 } else { 20.0 }) * (1.0 + fraction);
                    }
                }
                reference_mass_kg = finite(reference_mass_kg + mass)?;
                means.push(mean);
            } else {
                let points = if arity == 3 {
                    triangle_rule.as_slice()
                } else if uniform {
                    quad_rule.as_slice()
                } else {
                    &graded_quad_rule
                };
                let mut element_volume = 0.0;
                let mut element_mass = 0.0;
                for &(xi, eta, weight) in points {
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
                    let point_density = if uniform {
                        rho[0]
                    } else {
                        scale
                            * shapes
                                .iter()
                                .zip(normalized)
                                .map(|(n, r)| n * r)
                                .sum::<f64>()
                    };
                    let point_mass = finite(point_density * volume)?;
                    if point_mass <= 0.0 {
                        return Err(invalid("plane reference point mass underflowed"));
                    }
                    reference_mass_kg = finite(reference_mass_kg + point_mass)?;
                    if !uniform {
                        element_volume = finite(element_volume + volume)?;
                        element_mass = finite(element_mass + point_mass)?;
                    }
                    for a in 0..arity {
                        for b in a..arity {
                            local[a][b] = (point_mass * shapes[a]).mul_add(shapes[b], local[a][b]);
                        }
                    }
                }
                if nodal {
                    means.push(if uniform {
                        rho[0]
                    } else {
                        finite(element_mass / element_volume)?
                    });
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
        Self::from_matrix(matrix.assemble(), reference_mass_kg, ndofs, cx).map(|mass| (mass, means))
    }
}

#[cfg(test)]
mod tests;
