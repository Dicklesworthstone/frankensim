//! Explicit sampled support for scalar material laws of two state coordinates.
//!
//! Source triangles, not an inferred grid or convex hull, define interpolation
//! support. Exact predicates own membership and topology; floating arithmetic
//! evaluates the admitted affine law without claiming a numerical enclosure.

use std::collections::BTreeSet;

use fs_ivl::{Sign, orient2d};

use crate::{EvaluationDecision, InterpolationPolicy, MatDbError};

/// Maximum samples in one two-coordinate material response.
pub const MAX_SURFACE_SAMPLES: usize = 4_096;
/// Maximum source-declared triangles. Pairwise admission is bounded by this cap.
pub const MAX_SURFACE_TRIANGLES: usize = 1_024;

// The exact-predicate owner's no-underflow domain is explicit. In this range,
// coordinate differences and degree-two products, including expansion tails,
// remain normal and finite. Zero is admitted separately. These limits apply to
// source and query coordinates alike, without silently rescaling source bits.
const MIN_COORDINATE: f64 = f64::from_bits((1023 - 200) << 52);
const MAX_COORDINATE: f64 = f64::from_bits((1023 + 200) << 52);

fn coordinate_admitted(value: f64) -> bool {
    value.is_finite() && (value == 0.0 || (MIN_COORDINATE..=MAX_COORDINATE).contains(&value.abs()))
}

fn invalid(reason: &'static str) -> MatDbError {
    MatDbError::MalformedSurface { reason }
}

pub(crate) fn validate(
    axes: &[String; 2],
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
) -> Result<(), MatDbError> {
    if axes
        .iter()
        .any(|axis| axis.trim().is_empty() || axis.len() > crate::MAX_PROPERTY_USAGE_AXIS_BYTES)
        || axes[0] == axes[1]
    {
        return Err(invalid(
            "surface axes must be distinct, nonblank bounded names",
        ));
    }
    if samples.is_empty() || samples.len() > MAX_SURFACE_SAMPLES {
        return Err(invalid("surface sample count must be in 1..=4096"));
    }
    if triangles.len() > MAX_SURFACE_TRIANGLES {
        return Err(invalid("surface triangle count exceeds 1024"));
    }
    let mut coordinates = BTreeSet::new();
    for (point, value) in samples {
        if !point.iter().copied().all(coordinate_admitted) {
            return Err(invalid(
                "surface coordinates must be zero or have magnitude in [2^-200, 2^200]",
            ));
        }
        if !value.is_finite() {
            return Err(invalid("surface ordinates must be finite"));
        }
        let key = point.map(|v| if v == 0.0 { 0 } else { v.to_bits() });
        if !coordinates.insert(key) {
            return Err(invalid("surface sample coordinates must be unique"));
        }
    }
    for triangle in triangles {
        if triangle
            .iter()
            .any(|&index| index as usize >= samples.len())
            || triangle[0] == triangle[1]
            || triangle[1] == triangle[2]
            || triangle[2] == triangle[0]
        {
            return Err(invalid(
                "surface triangle requires three distinct existing samples",
            ));
        }
        let points = vertices(samples, *triangle);
        if orient2d(points[0], points[1], points[2]) == Sign::Zero
            || determinant(points[0], points[1], points[2]) == 0.0
        {
            return Err(invalid(
                "surface triangle is degenerate or numerically unresolved",
            ));
        }
        if samples.iter().enumerate().any(|(index, (point, _))| {
            !triangle.contains(&(index as u32)) && contains(points, *point)
        }) {
            return Err(invalid(
                "a source sample inside or on a cell must be one of its declared vertices",
            ));
        }
    }
    validate_cell_intersections(samples, triangles)
}

fn validate_cell_intersections(
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
) -> Result<(), MatDbError> {
    for (index, triangle) in triangles.iter().enumerate() {
        let points = vertices(samples, *triangle);
        for other in &triangles[..index] {
            let other_points = vertices(samples, *other);
            if boxes_disjoint(points, other_points) {
                continue;
            }
            let shared = triangle
                .iter()
                .filter(|vertex| other.contains(vertex))
                .count();
            if shared == 3 {
                return Err(invalid("surface triangles must not repeat the same cell"));
            }
            if shared == 2 {
                let edge: Vec<_> = triangle
                    .iter()
                    .filter(|vertex| other.contains(vertex))
                    .collect();
                let a = samples[*edge[0] as usize].0;
                let b = samples[*edge[1] as usize].0;
                let c = samples[*triangle
                    .iter()
                    .find(|vertex| !other.contains(vertex))
                    .unwrap() as usize]
                    .0;
                let d = samples[*other
                    .iter()
                    .find(|vertex| !triangle.contains(vertex))
                    .unwrap() as usize]
                    .0;
                if orient2d(a, b, c) == orient2d(a, b, d) {
                    return Err(invalid("surface triangles overlap across a shared edge"));
                }
                continue;
            }
            // Source admission already excludes any foreign vertex inside/on a
            // cell. Only crossing interiors can remain for non-shared edges.
            for edge in 0..3 {
                for other_edge in 0..3 {
                    if cross(
                        points[edge],
                        points[(edge + 1) % 3],
                        other_points[other_edge],
                        other_points[(other_edge + 1) % 3],
                    ) {
                        return Err(invalid("surface triangle edges cross"));
                    }
                }
            }
        }
    }
    Ok(())
}

fn vertices(samples: &[([f64; 2], f64)], triangle: [u32; 3]) -> [[f64; 2]; 3] {
    triangle.map(|index| samples[index as usize].0)
}

fn determinant(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])
}

fn opposite(a: Sign, b: Sign) -> bool {
    matches!(
        (a, b),
        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
    )
}

fn cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    opposite(orient2d(a, b, c), orient2d(a, b, d)) && opposite(orient2d(c, d, a), orient2d(c, d, b))
}

fn boxes_disjoint(a: [[f64; 2]; 3], b: [[f64; 2]; 3]) -> bool {
    (0..2).any(|axis| {
        let amin = a
            .iter()
            .map(|point| point[axis])
            .fold(f64::INFINITY, f64::min);
        let amax = a
            .iter()
            .map(|point| point[axis])
            .fold(f64::NEG_INFINITY, f64::max);
        let bmin = b
            .iter()
            .map(|point| point[axis])
            .fold(f64::INFINITY, f64::min);
        let bmax = b
            .iter()
            .map(|point| point[axis])
            .fold(f64::NEG_INFINITY, f64::max);
        amax < bmin || bmax < amin
    })
}

fn contains(triangle: [[f64; 2]; 3], point: [f64; 2]) -> bool {
    if boxes_disjoint(triangle, [point; 3]) {
        return false;
    }
    let winding = orient2d(triangle[0], triangle[1], triangle[2]);
    (0..3).all(|edge| {
        let sign = orient2d(triangle[edge], triangle[(edge + 1) % 3], point);
        sign == Sign::Zero || sign == winding
    })
}

pub(crate) fn evaluate(
    axes: &[String; 2],
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
    point: [f64; 2],
    policy: InterpolationPolicy,
) -> Result<(f64, EvaluationDecision), MatDbError> {
    if policy == InterpolationPolicy::ConstantWithinValidity {
        return Err(MatDbError::UnsupportedEvaluation {
            reason: "a sampled surface cannot be evaluated as a validity-wide constant",
        });
    }
    if !point.into_iter().all(coordinate_admitted) {
        return Err(MatDbError::UnsupportedEvaluation {
            reason: "surface query coordinates must be zero or have magnitude in [2^-200, 2^200]",
        });
    }
    if let Some((sample, (_, value))) = samples.iter().enumerate().find(|(_, (at, _))| {
        at.iter()
            .zip(point)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    }) {
        return Ok((
            *value,
            EvaluationDecision::ExactSurfaceSample {
                sample: sample as u32,
            },
        ));
    }
    if policy == InterpolationPolicy::LinearInside {
        for (index, &triangle) in triangles.iter().enumerate() {
            let positions = vertices(samples, triangle);
            if !contains(positions, point) {
                continue;
            }
            let mut weights = [0.0; 3];
            for vertex in 0..3 {
                let a = positions[(vertex + 1) % 3];
                let b = positions[(vertex + 2) % 3];
                if orient2d(a, b, point) != Sign::Zero {
                    weights[vertex] = determinant(a, b, point).abs();
                    if weights[vertex] == 0.0 {
                        return Err(MatDbError::UnsupportedEvaluation {
                            reason: "surface interpolation weights are numerically unresolved",
                        });
                    }
                }
            }
            let total: f64 = weights.iter().sum();
            let value: f64 = weights
                .iter()
                .zip(triangle)
                .map(|(weight, vertex)| (weight / total) * samples[vertex as usize].1)
                .sum();
            if !value.is_finite() || !total.is_finite() || total <= 0.0 {
                return Err(MatDbError::UnsupportedEvaluation {
                    reason: "surface interpolation is not representable as a finite scalar",
                });
            }
            return Ok((
                value,
                EvaluationDecision::LinearInsideTriangle {
                    triangle: index as u32,
                },
            ));
        }
    }
    Err(MatDbError::OutsideSurfaceSupport {
        axes: axes.clone(),
        requested: point,
    })
}
