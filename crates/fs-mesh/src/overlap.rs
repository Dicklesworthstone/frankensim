//! Measured tetrahedron/Cartesian-cell overlap and polynomial volume moments.
//!
//! Six axis-aligned half-space clips construct the actual convex intersection;
//! callers need not supply guessed overlap fractions. Moments integrate affine
//! mass/momentum fields and their quadratic kinetic/angular moments. This is
//! f64 geometry, not a certified intersection or a constitutive-history remap.

use fs_exec::Cx;
use fs_ivl::{Sign, orient3d};

#[cfg(test)]
mod tests;

type Point = [f64; 3];
type Face = Vec<Point>;
// A tetrahedron clipped by a box has at most ten supporting faces. The slack
// admits retained collinear edge points without an unbounded geometry path.
const MAX_FACES: usize = 16;
const MAX_FACE_VERTICES: usize = 16;

/// Integrals over the intersection, using the caller's moment origin.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VolumeMoments {
    /// Integral of one [length³]. Zero denotes an empty/zero-volume overlap.
    pub volume: f64,
    /// Integral of `(x - moment_origin)` [length⁴].
    pub first: [f64; 3],
    /// Integral of `(x - moment_origin)(x - moment_origin)^T` [length⁵].
    pub second: [[f64; 3]; 3],
}

/// Invalid geometry, unrepresentable metric arithmetic, or cancellation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OverlapError {
    InvalidInput,
    DegenerateTetrahedron,
    Unrepresentable,
    GeometryLimit,
    Cancelled,
}
impl core::fmt::Display for OverlapError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedron/box overlap: {self:?}")
    }
}
impl std::error::Error for OverlapError {}

/// Intersect a tetrahedron with a closed Cartesian box and integrate moments.
///
/// Both tetrahedron orientations are admitted. Coordinates share one length
/// unit; the moment origin is explicit and should be near the specimen. Local
/// normalization bounds predicate arguments and avoids world-origin volume
/// cancellation. No tolerance expands/shrinks the supplied geometry. Exact
/// signs apply to the normalized f64 tetrahedron; clipping intersections and
/// moments remain measured f64 quantities. Callers must check aggregate
/// coverage/conservation against their own error allowance.
///
/// Work is bounded by six clips, sixteen faces and sixteen vertices per face.
/// Cancellation is polled at each clip and each integrated face. The function
/// publishes a complete result or an error and never mutates the source mesh.
///
/// # Errors
/// Refuses non-finite inputs, unordered box bounds, degenerate tetrahedra,
/// unrepresentable local geometry/positive volume/moments and cancellation.
pub fn tetrahedron_box_moments(
    vertices: [[f64; 3]; 4],
    lower: [f64; 3],
    upper: [f64; 3],
    moment_origin: [f64; 3],
    cx: &Cx<'_>,
) -> Result<VolumeMoments, OverlapError> {
    poll(cx)?;
    if vertices
        .iter()
        .flatten()
        .chain(&lower)
        .chain(&upper)
        .chain(&moment_origin)
        .any(|x| !x.is_finite())
        || (0..3).any(|i| lower[i] >= upper[i])
    {
        return Err(OverlapError::InvalidInput);
    }
    let base = vertices[0];
    let mut local = vertices.map(|p| core::array::from_fn(|i| p[i] - base[i]));
    let scale = local
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() {
        return Err(OverlapError::Unrepresentable);
    }
    if scale == 0.0 {
        return Err(OverlapError::DegenerateTetrahedron);
    }
    for point in &mut local {
        for x in point {
            *x /= scale;
        }
    }
    let orientation = orient3d(local[0], local[1], local[2], local[3]);
    if orientation == Sign::Zero {
        return Err(OverlapError::DegenerateTetrahedron);
    }
    let source_det = determinant(
        sub(local[1], local[0]),
        sub(local[2], local[0]),
        sub(local[3], local[0]),
    );
    if source_det == 0.0
        || !source_det.is_finite()
        || (source_det > 0.0) != (orientation == Sign::Negative)
    {
        return Err(OverlapError::Unrepresentable);
    }
    let source_volume = ((source_det.abs() / 6.0 * scale) * scale) * scale;
    if !source_volume.is_finite() || source_volume <= 0.0 {
        return Err(OverlapError::Unrepresentable);
    }
    // Tighten to the tetrahedron's bounding box before subtraction/division:
    // a far-away but finite box must not overflow local-coordinate arithmetic.
    let mut lo = [0.0; 3];
    let mut hi = [0.0; 3];
    for i in 0..3 {
        let min = vertices.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min);
        let max = vertices
            .iter()
            .map(|p| p[i])
            .fold(f64::NEG_INFINITY, f64::max);
        let a = lower[i].max(min);
        let b = upper[i].min(max);
        if a >= b {
            return Ok(VolumeMoments::default());
        }
        lo[i] = (a - base[i]) / scale;
        hi[i] = (b - base[i]) / scale;
        if !lo[i].is_finite() || !hi[i].is_finite() || lo[i] >= hi[i] {
            return Err(OverlapError::Unrepresentable);
        }
    }
    let mut faces: Vec<Face> = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
        .map(|indices| indices.map(|i| local[i]).to_vec())
        .into();
    // orient3d has the opposite sign to det(v1-v0,v2-v0,v3-v0).
    if orientation == Sign::Positive {
        for face in &mut faces {
            face.reverse();
        }
    }
    for axis in 0..3 {
        for (bound, keep_below) in [(lo[axis], false), (hi[axis], true)] {
            poll(cx)?;
            faces = clip(faces, axis, bound, keep_below)?;
            if faces.is_empty() {
                return Ok(VolumeMoments::default());
            }
        }
    }
    integrate(&faces, base, scale, moment_origin, cx)
}

fn poll(cx: &Cx<'_>) -> Result<(), OverlapError> {
    cx.checkpoint().map_err(|_| OverlapError::Cancelled)
}
fn sub(a: Point, b: Point) -> Point {
    core::array::from_fn(|i| a[i] - b[i])
}
fn determinant(a: Point, b: Point, c: Point) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1])
        + a[1] * (b[2] * c[0] - b[0] * c[2])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}
fn order(a: &Point, b: &Point) -> core::cmp::Ordering {
    a[0].total_cmp(&b[0])
        .then(a[1].total_cmp(&b[1]))
        .then(a[2].total_cmp(&b[2]))
}

fn intersection(
    mut a: Point,
    mut b: Point,
    axis: usize,
    bound: f64,
) -> Result<Point, OverlapError> {
    // Neighboring faces visit an edge in opposite directions. A canonical
    // endpoint order gives the identical shared intersection in both faces.
    if order(&a, &b).is_gt() {
        core::mem::swap(&mut a, &mut b);
    }
    if a[axis] == bound {
        return Ok(a);
    }
    if b[axis] == bound {
        return Ok(b);
    }
    let t = (bound - a[axis]) / (b[axis] - a[axis]);
    if !t.is_finite() || !(0.0..=1.0).contains(&t) {
        return Err(OverlapError::Unrepresentable);
    }
    let mut p = core::array::from_fn(|i| a[i] + t * (b[i] - a[i]));
    p[axis] = bound;
    Ok(p)
}

fn clip(faces: Vec<Face>, axis: usize, bound: f64, below: bool) -> Result<Vec<Face>, OverlapError> {
    let inside = |p: Point| {
        if below {
            p[axis] <= bound
        } else {
            p[axis] >= bound
        }
    };
    if faces.iter().flatten().copied().all(inside) {
        return Ok(faces);
    }
    let mut output = Vec::new();
    let mut cap = Vec::new();
    for face in faces {
        let mut polygon = Vec::new();
        for i in 0..face.len() {
            let a = face[i];
            let b = face[(i + 1) % face.len()];
            if inside(a) {
                polygon.push(a);
            }
            if inside(a) != inside(b) {
                polygon.push(intersection(a, b, axis, bound)?);
            }
        }
        polygon.dedup();
        if polygon.len() > 1 && polygon.first() == polygon.last() {
            polygon.pop();
        }
        if polygon.len() > MAX_FACE_VERTICES {
            return Err(OverlapError::GeometryLimit);
        }
        for &point in &polygon {
            if point[axis] == bound && !cap.contains(&point) {
                cap.push(point);
            }
        }
        if polygon.len() >= 3 {
            output.push(polygon);
        }
    }
    if cap.len() > MAX_FACE_VERTICES {
        return Err(OverlapError::GeometryLimit);
    }
    if cap.len() >= 3 && !output.is_empty() {
        let center: Point =
            core::array::from_fn(|i| cap.iter().map(|p| p[i]).sum::<f64>() / cap.len() as f64);
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        cap.sort_by(|a, b| {
            fs_math::det::atan2(a[v] - center[v], a[u] - center[u])
                .total_cmp(&fs_math::det::atan2(b[v] - center[v], b[u] - center[u]))
                .then(order(a, b))
        });
        if !below {
            cap.reverse();
        }
        output.push(cap);
    }
    if output.len() > MAX_FACES {
        return Err(OverlapError::GeometryLimit);
    }
    Ok(output)
}

fn integrate(
    faces: &[Face],
    base: Point,
    scale: f64,
    origin: Point,
    cx: &Cx<'_>,
) -> Result<VolumeMoments, OverlapError> {
    let count = faces.iter().map(Vec::len).sum::<usize>();
    let center: Point =
        core::array::from_fn(|i| faces.iter().flatten().map(|p| p[i]).sum::<f64>() / count as f64);
    let offset = sub(base, origin);
    if offset.iter().any(|x| !x.is_finite()) {
        return Err(OverlapError::Unrepresentable);
    }
    let mut sums = [0.0; 13];
    let mut corrections = [0.0; 13];
    for face in faces {
        poll(cx)?;
        for i in 1..face.len() - 1 {
            let points = [center, face[0], face[i], face[i + 1]];
            let det = determinant(
                sub(points[1], center),
                sub(points[2], center),
                sub(points[3], center),
            );
            if det < 0.0 {
                return Err(OverlapError::Unrepresentable);
            }
            if det == 0.0 {
                continue;
            }
            let volume = ((det / 6.0 * scale) * scale) * scale;
            if !volume.is_finite() || volume <= 0.0 {
                return Err(OverlapError::Unrepresentable);
            }
            let r = points.map(|p| core::array::from_fn::<_, 3, _>(|i| offset[i] + scale * p[i]));
            let sum: Point = core::array::from_fn(|i| r.iter().map(|p| p[i]).sum());
            let mut values = [0.0; 13];
            values[0] = volume;
            for j in 0..3 {
                values[1 + j] = volume * (sum[j] / 4.0);
                for k in 0..3 {
                    values[4 + 3 * j + k] = volume
                        * ((r.iter().map(|p| p[j] * p[k]).sum::<f64>() + sum[j] * sum[k]) / 20.0);
                }
            }
            for j in 0..13 {
                if !values[j].is_finite() {
                    return Err(OverlapError::Unrepresentable);
                }
                let adjusted = values[j] - corrections[j];
                let next = sums[j] + adjusted;
                corrections[j] = (next - sums[j]) - adjusted;
                sums[j] = next;
            }
        }
    }
    poll(cx)?;
    if sums.iter().any(|x| !x.is_finite()) {
        return Err(OverlapError::Unrepresentable);
    }
    if sums[0] > 0.0 && (0..3).any(|i| sums[4 + 4 * i] <= 0.0) {
        return Err(OverlapError::Unrepresentable);
    }
    Ok(VolumeMoments {
        volume: sums[0],
        first: [sums[1], sums[2], sums[3]],
        second: core::array::from_fn(|i| core::array::from_fn(|j| sums[4 + 3 * i + j])),
    })
}
