//! Explicit sampled support for scalar material laws of two state coordinates.
//!
//! Source triangles, not an inferred grid or convex hull, define interpolation
//! support. Exact predicates own membership and topology; floating arithmetic
//! evaluates the admitted affine law without claiming a numerical enclosure.

use std::collections::{BTreeMap, BTreeSet};

use fs_ivl::expansion::{diff_expansion, estimate, expansion_diff, expansion_product};
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

// A genuine near-edge area can disappear in the rounded product difference
// even though the exact membership predicate resolves its nonzero sign. The
// original expansion owner retains BOTH coordinate-subtraction and product
// tails; rounding that exact determinant recovers the weight without moving
// the query, clipping it to an edge, or changing source-cell membership.
// Source/query coordinate admission already bounds this degree-two arithmetic
// inside the finite, no-underflow expansion domain. This is a cold fallback.
fn expanded_determinant(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let acx = diff_expansion(a[0], c[0]);
    let acy = diff_expansion(a[1], c[1]);
    let bcx = diff_expansion(b[0], c[0]);
    let bcy = diff_expansion(b[1], c[1]);
    estimate(&expansion_diff(
        &expansion_product(&acx, &bcy),
        &expansion_product(&acy, &bcx),
    ))
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

/// Source-cell continuity on a fixed-coordinate line. Only an already admitted
/// claim enters here. Connectivity is combinatorial: two cut cells join only
/// at their SAME declared vertex/edge intersecting the line. Thus rounding two
/// unrelated cuts to nearby numbers cannot fill a real source-support hole.
pub(crate) fn line_knots(
    samples: &[([f64; 2], f64)],
    triangles: &[[u32; 3]],
    varying: usize,
    lower: [f64; 2],
    upper_coordinate: f64,
) -> Result<Vec<f64>, MatDbError> {
    let fixed = 1 - varying;
    let level = lower[fixed];
    let mut upper = lower;
    upper[varying] = upper_coordinate;
    let refusal = |reason| MatDbError::UnsupportedEvaluation { reason };
    let active: Vec<_> = triangles
        .iter()
        .map(|triangle| {
            triangle
                .iter()
                .any(|&i| samples[i as usize].0[fixed] <= level)
                && triangle
                    .iter()
                    .any(|&i| samples[i as usize].0[fixed] >= level)
        })
        .collect();
    let mut visited = vec![false; triangles.len()];
    let mut pending = Vec::new();
    for (index, &triangle) in triangles.iter().enumerate() {
        if active[index] && contains(vertices(samples, triangle), lower) {
            visited[index] = true;
            pending.push(index);
        }
    }
    while let Some(index) = pending.pop() {
        for (other, triangle) in triangles.iter().enumerate() {
            if visited[other] || !active[other] {
                continue;
            }
            let mut shared = [0_u32; 2];
            let mut count = 0;
            for &vertex in &triangles[index] {
                if triangle.contains(&vertex) {
                    // Duplicate cells were refused at source admission.
                    shared[count] = vertex;
                    count += 1;
                }
            }
            let touches = match count {
                1 => samples[shared[0] as usize].0[fixed] == level,
                2 => {
                    let a = samples[shared[0] as usize].0[fixed];
                    let b = samples[shared[1] as usize].0[fixed];
                    a.min(b) <= level && level <= a.max(b)
                }
                _ => false,
            };
            if touches {
                visited[other] = true;
                pending.push(other);
            }
        }
    }
    if !triangles
        .iter()
        .enumerate()
        .any(|(index, &triangle)| visited[index] && contains(vertices(samples, triangle), upper))
    {
        return Err(refusal(
            "declared surface cells do not continuously cover the requested fixed-coordinate line",
        ));
    }

    // A vertex cut uses (i,i); an interior edge cut uses sorted vertex ids.
    // Shared cuts are computed once in that canonical orientation. Distinct
    // cuts collapsing to the same float refuse; no proximity tolerance exists.
    let mut cuts = BTreeMap::new();
    for (index, triangle) in triangles.iter().enumerate() {
        if !visited[index] {
            continue;
        }
        for edge in 0..3 {
            let i = triangle[edge].min(triangle[(edge + 1) % 3]);
            let j = triangle[edge].max(triangle[(edge + 1) % 3]);
            let a = samples[i as usize].0;
            let b = samples[j as usize].0;
            for (vertex, at) in [(i, a), (j, b)] {
                if at[fixed] == level
                    && lower[varying] <= at[varying]
                    && at[varying] <= upper_coordinate
                {
                    cuts.insert((vertex, vertex), at[varying]);
                }
            }
            if !(a[fixed].min(b[fixed]) < level && level < a[fixed].max(b[fixed]))
                || cuts.contains_key(&(i, j))
            {
                continue;
            }
            let fraction = (level - a[fixed]) / (b[fixed] - a[fixed]);
            let x = (b[varying] - a[varying]).mul_add(fraction, a[varying]);
            if !coordinate_admitted(x)
                || x < a[varying].min(b[varying])
                || x > a[varying].max(b[varying])
            {
                return Err(refusal(
                    "surface line intersection is numerically unresolved",
                ));
            }
            if x < lower[varying] || x > upper_coordinate {
                continue;
            }
            if (x == lower[varying] && orient2d(a, b, lower) != Sign::Zero)
                || (x == upper_coordinate && orient2d(a, b, upper) != Sign::Zero)
            {
                return Err(refusal(
                    "surface line intersection collapses onto a distinct endpoint",
                ));
            }
            cuts.insert((i, j), x);
        }
    }
    let mut ordered: Vec<_> = cuts.into_values().collect();
    ordered.sort_by(f64::total_cmp);
    if ordered.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(refusal(
            "distinct source line intersections collapse to one coordinate",
        ));
    }
    ordered.retain(|x| lower[varying] < *x && *x < upper_coordinate);
    Ok(ordered)
}

// Equal ordinates along a coordinate-aligned edge prove that this affine
// cell depends only on the other coordinate. Canonical endpoint ordering
// makes the same one-coordinate law identical in adjacent source cells.
fn axis_independent_value(
    positions: [[f64; 2]; 3],
    values: [f64; 3],
    point: [f64; 2],
) -> Option<f64> {
    for axis in 0..2 {
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            if positions[a][axis].to_bits() != positions[b][axis].to_bits()
                || values[a].to_bits() != values[b].to_bits()
            {
                continue;
            }
            let (low, high) = if positions[a][axis] < positions[c][axis] {
                (a, c)
            } else {
                (c, a)
            };
            let fraction = (point[axis] - positions[low][axis])
                / (positions[high][axis] - positions[low][axis]);
            let value = values[low] + fraction * (values[high] - values[low]);
            if value.is_finite() {
                return Some(value);
            }
        }
    }
    // Opposite extreme ordinates can overflow their difference even when
    // the original weighted sum is representable; retain that evaluator.
    None
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
                        weights[vertex] = expanded_determinant(a, b, point).abs();
                    }
                    if weights[vertex] == 0.0 {
                        return Err(MatDbError::UnsupportedEvaluation {
                            reason: "surface interpolation weights are numerically unresolved",
                        });
                    }
                }
            }
            let total: f64 = weights.iter().sum();
            if !total.is_finite() || total <= 0.0 {
                return Err(MatDbError::UnsupportedEvaluation {
                    reason: "surface interpolation is not representable as a finite scalar",
                });
            }
            // A constant cell or edge must retain its source value exactly:
            // rounded normalized weights need not sum to one. Only vertices
            // with exactly nonzero weights participate, so an off-edge vertex
            // cannot perturb a constant edge or hide variation inside a cell.
            let mut active_values = weights
                .iter()
                .zip(triangle)
                .filter(|(weight, _)| **weight != 0.0)
                .map(|(_, vertex)| samples[vertex as usize].1);
            let constant = active_values
                .next()
                .filter(|first| active_values.all(|value| value.to_bits() == first.to_bits()));
            let value = constant
                .or_else(|| {
                    axis_independent_value(
                        positions,
                        triangle.map(|v| samples[v as usize].1),
                        point,
                    )
                })
                .unwrap_or_else(|| {
                    weights
                        .iter()
                        .zip(triangle)
                        .map(|(weight, vertex)| (weight / total) * samples[vertex as usize].1)
                        .sum::<f64>()
                });
            if !value.is_finite() {
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
