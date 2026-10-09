//! Material-weighted Whitney mass on a genuine two-dimensional complex.
//!
//! The topology remains 2-D. `Metric2` supplies either a physical planar
//! thickness or the axisymmetric volume measure `angle * radius * dA`.
//! Barycentric moments integrate that affine radius exactly; multiplying a
//! planar mass matrix by the centroid radius would not do so. Coefficients
//! are constant on each cell, in the stored Cartesian/meridian frame.
//!
//! These are bulk Galerkin pairings, not interface laws or a complete
//! axisymmetric vector PDE. In particular, no cylindrical connection terms,
//! boundary conditions, gauge, or stability theorem are inferred here.

use crate::weighted::{
    WeightedAssemblyLimits, WeightedError, WeightedMass, assemble_triplets, poll, validate_tensor,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};

/// A cellwise material coefficient in the complex's two-coordinate frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CellWeight2 {
    /// Finite, strictly positive scalar; valid for degrees 0, 1, and 2.
    Scalar(f64),
    /// Symmetric positive-definite tensor; valid for degree-1 vector proxies.
    /// Its positivity must be established by outward interval principal minors.
    Tensor([[f64; 2]; 2]),
}

impl CellWeight2 {
    fn validate(self, cell: usize, degree: u8) -> Result<(), WeightedError> {
        let refuse = |reason| WeightedError::InvalidCoefficient { cell, reason };
        match self {
            Self::Scalar(value) => {
                if !value.is_finite() || value <= 0.0 {
                    return Err(refuse("scalar weight must be finite and positive"));
                }
            }
            Self::Tensor(tensor) => {
                if degree != 1 {
                    return Err(refuse("2-D tensor weights require degree 1"));
                }
                validate_tensor(&tensor, cell)?;
            }
        }
        Ok(())
    }

    fn scalar(self) -> f64 {
        match self {
            Self::Scalar(value) => value,
            Self::Tensor(_) => unreachable!("tensor admitted only for vector forms"),
        }
    }

    fn pairing(self, left: [f64; 2], right: [f64; 2]) -> f64 {
        match self {
            Self::Scalar(value) => value * left[0].mul_add(right[0], left[1] * right[1]),
            Self::Tensor(tensor) => left[0].mul_add(
                tensor[0][0].mul_add(right[0], tensor[0][1] * right[1]),
                left[1] * tensor[1][0].mul_add(right[0], tensor[1][1] * right[1]),
            ),
        }
    }

    fn normalized(self) -> (Self, f64) {
        match self {
            Self::Scalar(value) => (Self::Scalar(1.0), value),
            Self::Tensor(tensor) => {
                let scale = tensor.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
                (
                    Self::Tensor(tensor.map(|row| row.map(|value| value / scale))),
                    scale,
                )
            }
        }
    }
}

/// Assemble a weighted Galerkin mass/Hodge pairing on `TriComplex2`.
///
/// Degrees 0/1/2 use vertex hats, canonically oriented edge Whitney forms,
/// and face-integral-normalized constant 2-forms. Coordinates and planar
/// thickness are in metres. The entry dimension is
/// `coefficient_dims * metre^(3 - 2*degree)` for both admitted metrics.
/// This volume pairing does not turn the complex into a 3-D complex.
///
/// All size/coefficient admission precedes allocation. Each cell is one
/// bounded numerical tile; sorting, reduction and CSR publication also poll
/// `cancelled`. A caller-owned `Cx` can be captured by that callback. Refusal
/// returns no partial matrix. The result is reproducible for the same ordered
/// complex/coefficient input; arbitrary cell permutations may change rounding.
pub fn weighted_mass_matrix_2d(
    complex: &TriComplex2,
    degree: u8,
    weights: &[CellWeight2],
    coefficient_dims: Dims,
    limits: WeightedAssemblyLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WeightedMass, WeightedError> {
    poll(cancelled)?;
    let local_count: usize = match degree {
        0 | 1 => 3,
        2 => 1,
        _ => return Err(WeightedError::InvalidDegree(degree)),
    };
    let cells = complex.faces().len();
    let dofs = match degree {
        0 => complex.vertices().len(),
        1 => complex.edges().len(),
        _ => cells,
    };
    let triplets = cells
        .checked_mul(local_count * local_count)
        .ok_or(WeightedError::SizeOverflow)?;
    limits.check(cells, dofs, triplets)?;
    if weights.len() != cells {
        return Err(WeightedError::CoefficientCount {
            expected: cells,
            actual: weights.len(),
        });
    }
    let entry_dims = coefficient_dims
        .checked_plus(Dims([3 - 2 * degree as i8, 0, 0, 0, 0, 0]))
        .ok_or(WeightedError::DimensionOverflow)?;
    for (cell, &weight) in weights.iter().enumerate() {
        poll(cancelled)?;
        weight.validate(cell, degree)?;
    }
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(triplets)
        .map_err(|_| WeightedError::Allocation)?;
    for (cell, (&face, &weight)) in complex.faces().iter().zip(weights).enumerate() {
        poll(cancelled)?;
        let points = face.map(|v| complex.vertices()[v as usize]);
        let [a, b, c] = points;
        let determinant = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        let area = 0.5 * determinant.abs();
        if !area.is_finite() || area <= 0.0 {
            return Err(WeightedError::InvalidGeometry {
                cell,
                reason: "triangle area is not representable",
            });
        }
        let moments = barycentric_moments(complex.metric(), points, area);
        let mut emit = |row, column, value: f64| -> Result<(), WeightedError> {
            if !value.is_finite() {
                return Err(WeightedError::NonFiniteEntry);
            }
            if row == column && value <= 0.0 {
                return Err(WeightedError::InvalidCoefficient {
                    cell,
                    reason: "weighted diagonal is not representably positive",
                });
            }
            entries.push((row, column, value));
            Ok(())
        };
        match degree {
            0 => {
                for p in 0..3 {
                    for q in p..3 {
                        let value = weight.scalar() * moments[p][q];
                        emit(face[p] as usize, face[q] as usize, value)?;
                        if p != q {
                            emit(face[q] as usize, face[p] as usize, value)?;
                        }
                    }
                }
            }
            1 => {
                let (weight, scale) = weight.normalized();
                let gradients = [
                    [(b[1] - c[1]) / determinant, (c[0] - b[0]) / determinant],
                    [(c[1] - a[1]) / determinant, (a[0] - c[0]) / determinant],
                    [(a[1] - b[1]) / determinant, (b[0] - a[0]) / determinant],
                ];
                if gradients.iter().flatten().any(|value| !value.is_finite()) {
                    return Err(WeightedError::InvalidGeometry {
                        cell,
                        reason: "triangle gradients are not representable",
                    });
                }
                let mut edges = [(0_usize, 0_usize, 0_usize); 3];
                for (index, (p, q)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                    let (p, q) = if face[p] < face[q] { (p, q) } else { (q, p) };
                    let edge = complex.edge_index(face[p], face[q]).ok_or(
                        WeightedError::InvalidGeometry {
                            cell,
                            reason: "triangle edge is absent from the complex",
                        },
                    )?;
                    edges[index] = (edge, p, q);
                }
                for i in 0..3 {
                    let (ei, p, q) = edges[i];
                    for &(ej, r, s) in &edges[i..] {
                        // (lambda_p g_q - lambda_q g_p)^T W
                        // (lambda_r g_s - lambda_s g_r), exactly integrated.
                        let value = scale
                            * (moments[p][r] * weight.pairing(gradients[q], gradients[s])
                                - moments[p][s] * weight.pairing(gradients[q], gradients[r])
                                - moments[q][r] * weight.pairing(gradients[p], gradients[s])
                                + moments[q][s] * weight.pairing(gradients[p], gradients[r]));
                        emit(ei, ej, value)?;
                        if ei != ej {
                            emit(ej, ei, value)?;
                        }
                    }
                }
            }
            _ => {
                // The top-form dof is integral over the oriented embedded
                // triangle, not over its thickness/radial volume measure.
                let measure = complex
                    .face_measure(cell)
                    .ok_or(WeightedError::InvalidGeometry {
                        cell,
                        reason: "triangle measure is absent",
                    })?;
                emit(cell, cell, ((measure / area) / area) * weight.scalar())?;
            }
        }
    }
    let matrix = assemble_triplets(dofs, entries, cancelled)?;
    Ok(WeightedMass::from_matrix(
        matrix,
        degree,
        coefficient_dims,
        entry_dims,
    ))
}

/// Integral of lambda_p lambda_q against the declared physical measure.
/// For a triangle, integral(prod lambda_i^a_i) = 2 A prod(a_i!)/(2+sum a_i)!.
/// A linear radial factor requires cubic moments, not centroid quadrature.
fn barycentric_moments(metric: Metric2, points: [[f64; 2]; 3], area: f64) -> [[f64; 3]; 3] {
    let mut moments = [[0.0; 3]; 3];
    for p in 0..3 {
        for q in p..3 {
            let value = match metric {
                Metric2::Planar { thickness } => area * thickness / if p == q { 6.0 } else { 12.0 },
                Metric2::Axisymmetric { angular_span } => {
                    let mut mean = 0.0;
                    for (r, point) in points.iter().enumerate() {
                        let denominator = if p == q && q == r {
                            10.0
                        } else if p == q || p == r || q == r {
                            30.0
                        } else {
                            60.0
                        };
                        mean += point[0] / denominator;
                    }
                    (area * angular_span) * mean
                }
            };
            moments[p][q] = value;
            moments[q][p] = value;
        }
    }
    moments
}
