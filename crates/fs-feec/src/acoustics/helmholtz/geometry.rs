//! Quadratic geometry shared by acoustic volume, boundary and receiver operators.
//!
//! One physical midside point belongs to each original canonical mesh edge.
//! Adjacent elements therefore have the same quadratic face map. Geometry does
//! not change the hierarchical pressure basis or its orientation convention.
//! Admission uses the existing rigorous TET10 query over the entire closed
//! simplex, never only the eventual integration points. A negative stored
//! corner orientation is reordered for that proof alone; original pressure
//! indices and barycentric coordinates are retained everywhere else.
//!
//! The proof establishes strictly nonzero local orientation. Global injectivity,
//! nonoverlap of different elements, and a geometric approximation error relative
//! to an acquired surface remain separate mesh/source obligations. Quadrature
//! and inverse-point tolerances remain independent convergence axes.

use super::AcousticFeError;
use crate::whitney::ElementGeometry;
use fs_exec::Cx;
use fs_math::det;
use fs_query::tetrahedron::quadratic::{
    QuadraticTetPathError, QuadraticTetPathReport, QuadraticTetPathVerdict,
    quadratic_tetrahedron_orientation,
};
use fs_rep_mesh::TetComplex;

/// Explicit work, storage and determinant-quality limits of the original
/// rigorous whole-element orientation query. The same nodes are supplied at
/// both path endpoints because this acoustic geometry is stationary.
pub use fs_query::tetrahedron::quadratic::QuadraticTetPathControls as QuadraticGeometryControls;

// Match the original fs-query TET10 owner: 01,12,02,03,13,23.
const EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]];
const DL: [[f64; 3]; 4] = [
    [-1.0, -1.0, -1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// The original reason a quadratic element was refused. An unresolved proof
/// carries its bounded material domain rather than being relabelled a fold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QuadraticGeometryRefusal {
    /// The query refused its work, storage, input or numerical range.
    Query(QuadraticTetPathError),
    /// Whole-element positive orientation could not be established.
    Unresolved(QuadraticTetPathReport),
    /// Geometry tables or an element request are inconsistent.
    InvalidInput(&'static str),
    /// A sampled physical point, Jacobian or inverse is unrepresentable.
    NumericalRange,
}

/// Physical geometry at one point of the original reference tetrahedron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryPoint {
    /// Interpolated position [m].
    pub position_m: [f64; 3],
    /// Positive volume measure `abs(det(dX/d(lambda1,lambda2,lambda3)))` [m³].
    /// Multiply a reference tetrahedron weight, whose sum is 1/6, by this.
    pub determinant_m3: f64,
    /// Original barycentric gradients `grad_X(lambda_a)` [1/m].
    pub barycentric_gradients: [[f64; 3]; 4],
}

/// Pointwise physical surface geometry on an original tetrahedron face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryFacePoint {
    /// The same volume map evaluated at this face point.
    pub point: GeometryPoint,
    /// Unit normal pointing out of the selected element, regardless of winding.
    pub outward_normal: [f64; 3],
    /// Surface measure [m²] for a triangle rule whose weights sum to 1/2.
    pub surface_jacobian_m2: f64,
}

impl GeometryPoint {
    /// Cofactor surface metric for the face `lambda_omitted = 0`.
    ///
    /// `-abs(det J) grad(lambda_omitted)` is the outward area vector for
    /// reference triangle coordinates. Using the same inverse as the volume
    /// operator preserves the physical normal even for curved shared faces.
    pub fn face(self, omitted: usize) -> Result<GeometryFacePoint, AcousticFeError> {
        let gradient =
            self.barycentric_gradients
                .get(omitted)
                .ok_or(AcousticFeError::InvalidInput {
                    what: "surface geometry requires one omitted local vertex in 0..4",
                })?;
        let scale = gradient.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        if !scale.is_finite() || scale <= 0.0 {
            return Err(AcousticFeError::InvalidInput {
                what: "surface geometry has an unrepresentable normal",
            });
        }
        let scaled = gradient.map(|v| v / scale);
        let length = det::sqrt(dot(scaled, scaled));
        let surface_jacobian_m2 = (self.determinant_m3 * scale) * length;
        let outward_normal = scaled.map(|v| -v / length);
        if !surface_jacobian_m2.is_finite()
            || surface_jacobian_m2 <= 0.0
            || outward_normal.iter().any(|v| !v.is_finite())
        {
            return Err(AcousticFeError::InvalidInput {
                what: "surface geometry has an unrepresentable area measure",
            });
        }
        Ok(GeometryFacePoint {
            point: self,
            outward_normal,
            surface_jacobian_m2,
        })
    }
}

#[derive(Debug, Clone)]
struct ElementMap {
    edges: [usize; 6],
    orientation: f64,
    // The original query's exact power-of-two normalization. Reuse it in
    // point evaluation instead of cubing dimensioned Jacobian entries.
    normalization_length_m: f64,
}

/// Only HelmholtzFe can own this proof together with its original coordinates;
/// no public sampler accepts a replacement, unproved coordinate table.
#[derive(Debug, Clone)]
pub(super) struct QuadraticGeometry {
    edge_midpoints_m: Vec<[f64; 3]>,
    elements: Vec<ElementMap>,
}

impl QuadraticGeometry {
    pub(super) fn new(
        complex: &TetComplex,
        positions: &[[f64; 3]],
        edge_midpoints_m: &[[f64; 3]],
        controls: QuadraticGeometryControls,
        cx: &Cx<'_>,
    ) -> Result<Self, AcousticFeError> {
        cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
        controls
            .validate()
            .map_err(|_| AcousticFeError::InvalidInput {
                what: "quadratic geometry requires valid explicit orientation-query limits",
            })?;
        if positions.len() != complex.vertex_count
            || positions.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(AcousticFeError::InvalidPositions {
                vertices: complex.vertex_count,
                positions: positions.len(),
            });
        }
        if edge_midpoints_m.len() != complex.edges.len()
            || edge_midpoints_m.iter().flatten().any(|v| !v.is_finite())
            || complex.edges.windows(2).any(|w| w[0] >= w[1])
            || complex
                .edges
                .iter()
                .any(|e| e[0] >= e[1] || e[1] as usize >= positions.len())
        {
            return Err(AcousticFeError::InvalidInput {
                what: "quadratic geometry requires one finite midside point per canonical mesh edge",
            });
        }
        let mut elements = Vec::with_capacity(complex.tets.len());
        for (element, tet) in complex.tets.iter().enumerate() {
            cx.checkpoint().map_err(|_| AcousticFeError::Cancelled)?;
            let mut sorted = *tet;
            sorted.sort_unstable();
            if sorted.windows(2).any(|w| w[0] == w[1])
                || sorted.iter().any(|&v| v as usize >= positions.len())
            {
                return Err(refusal(
                    element,
                    QuadraticGeometryRefusal::InvalidInput(
                        "quadratic geometry requires distinct original tetrahedron vertices",
                    ),
                ));
            }
            let corners = tet.map(|v| positions[v as usize]);
            let orientation = corner_orientation(corners)
                .ok_or(AcousticFeError::DegenerateTet { tet: element })?;
            let mut edges = [0; 6];
            let mut nodes = [[0.0; 3]; 10];
            nodes[..4].copy_from_slice(&corners);
            for (local, &[a, b]) in EDGES.iter().enumerate() {
                let mut key = [tet[a], tet[b]];
                key.sort_unstable();
                edges[local] = complex.edges.binary_search(&key).map_err(|_| {
                    refusal(
                        element,
                        QuadraticGeometryRefusal::InvalidInput(
                            "original tetrahedron edge missing from canonical edge table",
                        ),
                    )
                })?;
                nodes[4 + local] = edge_midpoints_m[edges[local]];
            }
            if orientation < 0.0 {
                // Swap local corners 1 and 2 for proof only, also permuting
                // their six midsides. Original H1 coordinates never change.
                nodes = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8].map(|i| nodes[i]);
            }
            let report =
                quadratic_tetrahedron_orientation(nodes, nodes, controls, cx).map_err(|error| {
                    match error {
                        QuadraticTetPathError::Cancelled => AcousticFeError::Cancelled,
                        _ => refusal(element, QuadraticGeometryRefusal::Query(error)),
                    }
                })?;
            if !matches!(report.verdict, QuadraticTetPathVerdict::Positive { .. }) {
                return Err(refusal(
                    element,
                    QuadraticGeometryRefusal::Unresolved(report),
                ));
            }
            elements.push(ElementMap {
                edges,
                orientation,
                normalization_length_m: report.normalization_length,
            });
        }
        Ok(Self {
            edge_midpoints_m: edge_midpoints_m.to_vec(),
            elements,
        })
    }

    /// Evaluate with the same original geometry admitted above. Coordinates
    /// may extrapolate during a bounded inverse solve, but only the closed
    /// simplex has whole-element orientation evidence. Nonfinite values and
    /// changed determinant sign always refuse.
    pub(super) fn point(
        &self,
        complex: &TetComplex,
        positions: &[[f64; 3]],
        element: usize,
        lambda: [f64; 4],
    ) -> Result<GeometryPoint, AcousticFeError> {
        validate_lambda(lambda)?;
        let (tet, map) = self.element(complex, element)?;
        let corners = tet.map(|v| positions[v as usize]);
        let origin = corners[0];
        let scale = map.normalization_length_m;
        let offsets =
            corners.map(|p| core::array::from_fn::<_, 3, _>(|i| (p[i] - origin[i]) / scale));
        let mut x = [0.0; 3];
        let mut jacobian = [[0.0; 3]; 3];
        for i in 0..3 {
            for a in 1..4 {
                x[i] += lambda[a] * offsets[a][i];
                jacobian[i][a - 1] = offsets[a][i];
            }
        }
        // X = X_affine + 4 sum_e lambda_a lambda_b delta_e. This form
        // removes the partition-of-unity cancellation in a ten-node sum,
        // and is exactly affine when every midpoint offset is zero.
        for (edge, &[a, b]) in EDGES.iter().enumerate() {
            let midpoint = self.edge_midpoints_m[map.edges[edge]];
            let shape = 4.0 * lambda[a] * lambda[b];
            let reference_gradient = core::array::from_fn::<_, 3, _>(|j| {
                4.0 * (lambda[a] * DL[b][j] + lambda[b] * DL[a][j])
            });
            for i in 0..3 {
                let midpoint_offset = (midpoint[i] - origin[i]) / scale;
                let delta = midpoint_offset - 0.5 * (offsets[a][i] + offsets[b][i]);
                x[i] += shape * delta;
                for j in 0..3 {
                    jacobian[i][j] += reference_gradient[j] * delta;
                }
            }
        }
        let cofactor = cofactors(jacobian);
        let determinant = dot(jacobian[0], cofactor[0]);
        let determinant_m3 = ((determinant.abs() * scale) * scale) * scale;
        let inverse_denominator = determinant * scale;
        if !(determinant * map.orientation > 0.0)
            || !determinant_m3.is_finite()
            || determinant_m3 <= 0.0
            || !inverse_denominator.is_finite()
            || inverse_denominator == 0.0
        {
            return Err(refusal(element, QuadraticGeometryRefusal::NumericalRange));
        }
        let mut barycentric_gradients = [[0.0; 3]; 4];
        for a in 1..4 {
            for i in 0..3 {
                barycentric_gradients[a][i] = cofactor[i][a - 1] / inverse_denominator;
            }
        }
        barycentric_gradients[0] = core::array::from_fn(|i| {
            -(barycentric_gradients[1][i]
                + barycentric_gradients[2][i]
                + barycentric_gradients[3][i])
        });
        let position_m = core::array::from_fn(|i| origin[i] + scale * x[i]);
        if position_m
            .iter()
            .chain(barycentric_gradients.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(refusal(element, QuadraticGeometryRefusal::NumericalRange));
        }
        Ok(GeometryPoint {
            position_m,
            determinant_m3,
            barycentric_gradients,
        })
    }

    pub(super) fn face_point(
        &self,
        complex: &TetComplex,
        positions: &[[f64; 3]],
        element: usize,
        omitted: usize,
        lambda: [f64; 4],
    ) -> Result<GeometryFacePoint, AcousticFeError> {
        self.point(complex, positions, element, lambda)?
            .face(omitted)
    }

    /// Conservative admission for legacy APIs which promise a single normal
    /// and an affine three-corner triangle. A planar yet nonlinear face map is
    /// deliberately not classified as affine.
    pub(super) fn face_is_affine(
        &self,
        complex: &TetComplex,
        positions: &[[f64; 3]],
        element: usize,
        omitted: usize,
    ) -> Result<bool, AcousticFeError> {
        let (tet, map) = self.element(complex, element)?;
        if omitted >= 4 {
            return Err(AcousticFeError::InvalidInput {
                what: "surface geometry requires one omitted local vertex in 0..4",
            });
        }
        for (edge, &[a, b]) in EDGES.iter().enumerate() {
            if a == omitted || b == omitted {
                continue;
            }
            let midpoint = self.edge_midpoints_m[map.edges[edge]];
            let pa = positions[tet[a] as usize];
            let pb = positions[tet[b] as usize];
            if (0..3).any(|i| midpoint[i] != f64::midpoint(pa[i], pb[i])) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn element<'a>(
        &'a self,
        complex: &TetComplex,
        element: usize,
    ) -> Result<([u32; 4], &'a ElementMap), AcousticFeError> {
        match (complex.tets.get(element), self.elements.get(element)) {
            (Some(tet), Some(map)) => Ok((*tet, map)),
            _ => Err(AcousticFeError::InvalidInput {
                what: "geometry requires an original tetrahedron index",
            }),
        }
    }
}

/// Existing affine metrics are passed through unchanged. This is the geometry
/// hook for ordinary HelmholtzFe constructors, not a fresh inverse per point.
pub(super) fn affine_point(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    element: usize,
    lambda: [f64; 4],
    metrics: &ElementGeometry,
) -> Result<GeometryPoint, AcousticFeError> {
    validate_lambda(lambda)?;
    let tet = complex
        .tets
        .get(element)
        .ok_or(AcousticFeError::InvalidInput {
            what: "geometry requires an original tetrahedron index",
        })?;
    let mut anchor = 0;
    for local in 1..4 {
        if lambda[local] > lambda[anchor] {
            anchor = local;
        }
    }
    let origin = positions[tet[anchor] as usize];
    let mut position_m = origin;
    // Match the original pressure probe exactly: inverse-map residuals and
    // the eventual material query must refer to the very same floating point.
    // A contributing anchor and skipped zero weights retain exact vertices
    // and constant face coordinates without extrapolating a boundary source.
    for (local, (&vertex, &bary)) in tet.iter().zip(&lambda).enumerate() {
        if local == anchor || bary == 0.0 {
            continue;
        }
        for (axis, position) in position_m.iter_mut().enumerate() {
            let offset = positions[vertex as usize][axis] - origin[axis];
            *position = bary.mul_add(offset, *position);
        }
    }
    let determinant_m3 = 6.0 * metrics.vol_signed[element].abs();
    if position_m.iter().any(|v| !v.is_finite())
        || !determinant_m3.is_finite()
        || determinant_m3 <= 0.0
    {
        return Err(AcousticFeError::DegenerateTet { tet: element });
    }
    Ok(GeometryPoint {
        position_m,
        determinant_m3,
        barycentric_gradients: metrics.grads[element],
    })
}

fn validate_lambda(lambda: [f64; 4]) -> Result<(), AcousticFeError> {
    if lambda.iter().any(|v| !v.is_finite())
        || (lambda.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(AcousticFeError::InvalidInput {
            what: "geometry requires finite barycentric coordinates summing to one",
        });
    }
    Ok(())
}

fn refusal(tet: usize, reason: QuadraticGeometryRefusal) -> AcousticFeError {
    AcousticFeError::CurvedGeometry { tet, reason }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cofactors(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            m[1][1] * m[2][2] - m[1][2] * m[2][1],
            m[1][2] * m[2][0] - m[1][0] * m[2][2],
            m[1][0] * m[2][1] - m[1][1] * m[2][0],
        ],
        [
            m[0][2] * m[2][1] - m[0][1] * m[2][2],
            m[0][0] * m[2][2] - m[0][2] * m[2][0],
            m[0][1] * m[2][0] - m[0][0] * m[2][1],
        ],
        [
            m[0][1] * m[1][2] - m[0][2] * m[1][1],
            m[0][2] * m[1][0] - m[0][0] * m[1][2],
            m[0][0] * m[1][1] - m[0][1] * m[1][0],
        ],
    ]
}

fn corner_orientation(corners: [[f64; 3]; 4]) -> Option<f64> {
    let edges: [[f64; 3]; 3] =
        core::array::from_fn(|a| core::array::from_fn(|i| corners[a + 1][i] - corners[0][i]));
    let scale = edges.iter().flatten().fold(0.0_f64, |m, x| m.max(x.abs()));
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    let normalized = edges.map(|e| e.map(|v| v / scale));
    let determinant = dot(normalized[0], cofactors(normalized)[0]);
    if !determinant.is_finite() || determinant == 0.0 {
        None
    } else {
        Some(if determinant > 0.0 { 1.0 } else { -1.0 })
    }
}

#[cfg(test)]
mod tests;
