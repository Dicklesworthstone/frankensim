//! Deterministic triangulation of one simple, planar polygon.
//!
//! Vertex IDs and winding are retained; no vertex is moved or inserted. This
//! is a bounded floating-point materializer, not an exact-predicate topology
//! certificate. Near-degenerate/ambiguous input may refuse. Callers admitting
//! untrusted meshes must still run their usual topology and repair boundary.

use fs_geom::Point3;
use std::fmt;

/// Hard bound on one polygon's input vertices.
pub const MAX_POLYGON_VERTICES: usize = 4096;
/// Hard bound on edge-pair and ear-containment work for one polygon.
pub const MAX_POLYGON_WORK: usize = 16_000_000;

/// Why a polygon could not be materialized without guessing its surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolygonError {
    /// Too few vertices, a bad index, or a non-finite coordinate.
    Input(&'static str),
    /// A vertex-count, work, or allocation bound was reached.
    Resource(&'static str),
    /// The polygon is degenerate or numerically ambiguous at this scale.
    Degenerate,
    /// The vertices do not define a plane to the admitted roundoff tolerance.
    NonPlanar,
    /// Edges cross, overlap, or touch other than at consecutive endpoints.
    NonSimple,
}

impl fmt::Display for PolygonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input(what) | Self::Resource(what) => f.write_str(what),
            Self::Degenerate => f.write_str("polygon is degenerate or numerically ambiguous"),
            Self::NonPlanar => f.write_str("polygon is not planar; supply an explicit triangulation"),
            Self::NonSimple => f.write_str("polygon has crossing, touching, or overlapping edges"),
        }
    }
}
impl std::error::Error for PolygonError {}

fn reserved<T>(count: usize) -> Result<Vec<T>, PolygonError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count)
        .map_err(|_| PolygonError::Resource("polygon allocation failed"))?;
    Ok(values)
}

fn spend(work: &mut usize) -> Result<(), PolygonError> {
    if *work == MAX_POLYGON_WORK {
        return Err(PolygonError::Resource("polygon triangulation work bound exceeded"));
    }
    *work += 1;
    Ok(())
}

// Normalized coordinates keep these products bounded. Ambiguous turns are
// treated conservatively: they cannot establish a convex ear, and they do
// not establish that two potentially intersecting edges are disjoint.
fn turn(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> i8 {
    let left = (b[0] - a[0]) * (c[1] - a[1]);
    let right = (b[1] - a[1]) * (c[0] - a[0]);
    let det = left - right;
    let uncertainty = 16.0 * f64::EPSILON * (left.abs() + right.abs())
        + f64::MIN_POSITIVE;
    if det > uncertainty { 1 } else if det < -uncertainty { -1 } else { 0 }
}

fn may_intersect(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    for axis in 0..2 {
        if a[axis].max(b[axis]) < c[axis].min(d[axis])
            || c[axis].max(d[axis]) < a[axis].min(b[axis])
        {
            return false;
        }
    }
    turn(a, b, c) * turn(a, b, d) <= 0 && turn(c, d, a) * turn(c, d, b) <= 0
}

/// Triangulate a simple planar polygon, retaining all boundary vertices.
///
/// Output uses the original IDs, follows the input winding, and contains
/// exactly `indices.len() - 2` nondegenerate triangles. Convex faces retain
/// the conventional first-vertex fan. Collinear boundary vertices are not
/// silently dropped, which would introduce cracks against neighboring faces.
///
/// Limits are [`MAX_POLYGON_VERTICES`] and [`MAX_POLYGON_WORK`]. There is no
/// hole, self-intersection, nonplanar-face, or exact-topology admission claim.
/// Planarity is tested after translation and scaling, with a roundoff-scale
/// tolerance of `64 * EPSILON * vertex_count`; coordinates are never repaired.
/// This synchronous bounded API does not claim cancellation support.
///
/// # Errors
/// Returns [`PolygonError`] instead of guessing a triangulation for malformed,
/// ambiguous, nonsimple, nonplanar, or over-budget input. No partial output is
/// published on failure.
pub fn triangulate_polygon(
    positions: &[Point3],
    indices: &[u32],
) -> Result<Vec<[u32; 3]>, PolygonError> {
    let n = indices.len();
    if n < 3 { return Err(PolygonError::Input("polygon needs at least three vertices")); }
    if n > MAX_POLYGON_VERTICES {
        return Err(PolygonError::Resource("polygon vertex count exceeds the face cap"));
    }
    let mut points = reserved(n)?;
    for &index in indices {
        let p = positions.get(index as usize)
            .ok_or(PolygonError::Input("polygon vertex index out of range"))?;
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return Err(PolygonError::Input("polygon contains a non-finite coordinate"));
        }
        points.push([p.x, p.y, p.z]);
    }
    // Subtract first to retain small features far from the origin. Only halve
    // globally when a subtraction would overflow; halving is an exact binary
    // rescale for ordinary coordinates and avoids an infinite extent.
    let origin = points[0];
    let halve = points.iter().any(|p| (0..3).any(|a| !(p[a] - origin[a]).is_finite()));
    let mut extent = 0.0_f64;
    for p in &mut points {
        for axis in 0..3 {
            p[axis] = if halve { 0.5 * p[axis] - 0.5 * origin[axis] }
                else { p[axis] - origin[axis] };
            extent = extent.max(p[axis].abs());
        }
    }
    if !(extent.is_finite() && extent > 0.0) { return Err(PolygonError::Degenerate); }
    for p in &mut points { for value in p { *value /= extent; } }

    // Newell's area normal chooses a stable cyclic projection. Its magnitude
    // also refuses zero-area bow ties and faces below the roundoff scale.
    let mut normal = [0.0_f64; 3];
    let mut compensation = [0.0_f64; 3];
    for i in 0..n {
        let a = points[i];
        let b = points[(i + 1) % n];
        for axis in 0..3 {
            let j = (axis + 1) % 3;
            let k = (axis + 2) % 3;
            let term = a[j] * b[k] - a[k] * b[j] - compensation[axis];
            let sum = normal[axis] + term;
            compensation[axis] = (sum - normal[axis]) - term;
            normal[axis] = sum;
        }
    }
    let mut drop_axis = 0;
    for axis in 1..3 {
        if normal[axis].abs() > normal[drop_axis].abs() { drop_axis = axis; }
    }
    let tolerance = 64.0 * f64::EPSILON * n as f64;
    let magnitude = normal[drop_axis].abs();
    if magnitude <= tolerance { return Err(PolygonError::Degenerate); }
    let winding = if normal[drop_axis] > 0.0 { 1 } else { -1 };
    for component in &mut normal { *component /= magnitude; }
    for p in &points {
        let distance = p[0] * normal[0] + p[1] * normal[1] + p[2] * normal[2];
        if distance.abs() > tolerance { return Err(PolygonError::NonPlanar); }
    }
    let mut projected = reserved(n)?;
    for p in &points { projected.push([p[(drop_axis + 1) % 3], p[(drop_axis + 2) % 3]]); }

    let mut work = 0;
    for i in 0..n {
        let a = projected[i];
        let b = projected[(i + 1) % n];
        if a == b { return Err(PolygonError::NonSimple); }
        // Adjacent collinear edges may continue straight, but must not double
        // back. The latter overlap is excluded from the nonadjacent-pair scan.
        let c = projected[(i + 2) % n];
        if turn(a, b, c) == 0
            && ((a[0] - b[0]) * (c[0] - b[0]) + (a[1] - b[1]) * (c[1] - b[1])) > 0.0
        {
            return Err(PolygonError::NonSimple);
        }
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) { continue; }
            spend(&mut work)?;
            if may_intersect(a, b, projected[j], projected[(j + 1) % n]) {
                return Err(PolygonError::NonSimple);
            }
        }
    }

    let mut remaining = reserved(n)?;
    remaining.extend(0..n);
    let mut triangles = reserved(n - 2)?;
    while remaining.len() > 3 {
        let m = remaining.len();
        let mut ear = None;
        // Starting at vertex 1 retains existing fan output for convex input.
        for offset in 1..=m {
            spend(&mut work)?;
            let at = offset % m;
            let a = remaining[(at + m - 1) % m];
            let b = remaining[at];
            let c = remaining[(at + 1) % m];
            if turn(projected[a], projected[b], projected[c]) != winding { continue; }
            let mut blocked = false;
            for &p in &remaining {
                if p == a || p == b || p == c { continue; }
                spend(&mut work)?;
                if turn(projected[a], projected[b], projected[p]) * winding >= 0
                    && turn(projected[b], projected[c], projected[p]) * winding >= 0
                    && turn(projected[c], projected[a], projected[p]) * winding >= 0
                {
                    blocked = true;
                    break;
                }
            }
            if !blocked { ear = Some((at, [indices[a], indices[b], indices[c]])); break; }
        }
        let Some((at, triangle)) = ear else { return Err(PolygonError::Degenerate); };
        triangles.push(triangle);
        remaining.remove(at);
    }
    let [a, b, c] = [remaining[0], remaining[1], remaining[2]];
    if turn(projected[a], projected[b], projected[c]) != winding {
        return Err(PolygonError::Degenerate);
    }
    triangles.push([indices[a], indices[b], indices[c]]);
    Ok(triangles)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn points(xy: &[[f64; 2]]) -> Vec<Point3> {
        xy.iter().map(|p| Point3::new(p[0], p[1], 0.0)).collect()
    }
    fn twice_area(p: &[Point3], t: [u32; 3]) -> f64 {
        let [a, b, c] = t.map(|i| p[i as usize]);
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }

    #[test]
    fn concave_notch_has_no_outside_fan_triangles_in_either_winding() {
        let p = points(&[[0., 0.], [3., 0.], [3., 3.], [2., 3.],
            [2., 1.], [1., 1.], [1., 3.], [0., 3.]]);
        for reverse in [false, true] {
            let mut ids: Vec<u32> = (0..8).collect();
            if reverse { ids.reverse(); }
            let triangles = triangulate_polygon(&p, &ids).unwrap();
            assert_eq!(triangles.len(), 6);
            let sign = if reverse { -1.0 } else { 1.0 };
            assert!(triangles.iter().all(|&t| twice_area(&p, t) * sign > 0.0));
            assert_eq!(triangles.iter().map(|&t| twice_area(&p, t)).sum::<f64>(), sign * 14.0);
            for t in &triangles {
                let x = t.iter().map(|&i| p[i as usize].x).sum::<f64>() / 3.0;
                let y = t.iter().map(|&i| p[i as usize].y).sum::<f64>() / 3.0;
                assert!(!(x > 1.0 && x < 2.0 && y > 1.0));
            }
            assert_eq!(triangles, triangulate_polygon(&p, &ids).unwrap());
        }
    }

    #[test]
    fn convex_fan_and_collinear_boundary_edges_are_preserved() {
        let p = points(&[[0., 0.], [1., 0.], [1., 1.], [0., 1.]]);
        assert_eq!(triangulate_polygon(&p, &[0, 1, 2, 3]).unwrap(), vec![[0, 1, 2], [0, 2, 3]]);
        let p = points(&[[0., 0.], [1., 0.], [2., 0.], [2., 2.], [0., 2.]]);
        let t = triangulate_polygon(&p, &[0, 1, 2, 3, 4]).unwrap();
        assert_eq!(t.len(), 3);
        for edge in [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]] {
            assert_eq!(t.iter().filter(|t| (0..3).any(|i| [t[i], t[(i + 1) % 3]] == edge)).count(), 1);
        }
        assert_eq!(t.iter().map(|&t| twice_area(&p, t)).sum::<f64>(), 8.0);
    }

    #[test]
    fn crossed_touching_repeated_and_nonplanar_polygons_refuse() {
        for p in [
            points(&[[0., 0.], [2., 2.], [0., 2.], [2., 0.]]),
            points(&[[0., 0.], [2., 0.], [1., 0.], [2., 2.], [0., 2.]]),
            points(&[[0., 0.], [2., 0.], [2., 2.], [0., 0.], [0., 2.]]),
        ] {
            let ids: Vec<u32> = (0..p.len() as u32).collect();
            assert!(triangulate_polygon(&p, &ids).is_err());
        }
        let mut p = points(&[[0., 0.], [1., 0.], [1., 1.], [0., 1.]]);
        p[2].z = 0.1;
        assert_eq!(triangulate_polygon(&p, &[0, 1, 2, 3]), Err(PolygonError::NonPlanar));
    }

    #[test]
    fn arbitrary_plane_and_extreme_finite_scales() {
        let xy = [[0., 0.], [3., 0.], [3., 3.], [1., 1.], [0., 3.]];
        for scale in [1e-200, 1.0, 1e200] {
            let p: Vec<_> = xy.iter().map(|&[x, y]| Point3::new(scale * x, scale * y, scale * (x + y))).collect();
            assert_eq!(triangulate_polygon(&p, &[0, 1, 2, 3, 4]).unwrap().len(), 3);
        }
        let p = points(&[[-1e308, -1e308], [1e308, -1e308], [1e308, 1e308], [-1e308, 1e308]]);
        assert_eq!(triangulate_polygon(&p, &[0, 1, 2, 3]).unwrap().len(), 2);
    }

    #[test]
    fn malformed_and_resource_inputs_refuse() {
        let p = points(&[[0., 0.], [1., 0.], [0., 1.]]);
        assert!(triangulate_polygon(&p, &[0, 1]).is_err());
        assert!(triangulate_polygon(&p, &[0, 1, 3]).is_err());
        assert!(matches!(triangulate_polygon(&p, &vec![0; MAX_POLYGON_VERTICES + 1]), Err(PolygonError::Resource(_))));
        let mut p = p;
        p[0].x = f64::NAN;
        assert!(triangulate_polygon(&p, &[0, 1, 2]).is_err());
        let mut work = MAX_POLYGON_WORK;
        assert!(matches!(spend(&mut work), Err(PolygonError::Resource(_))));
    }
}
