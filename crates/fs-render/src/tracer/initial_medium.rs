//! Launch-point medium ownership shared by camera, probe and light paths.
//!
//! A first back-facing hit is not an inside test: an open sheet or an overlap
//! can produce the same hit. Meshes first need opposite paired edges, then the
//! complete oriented solid angle decides membership. Nested launch regions
//! additionally need a conservative containment witness; no arbitrary order is
//! inferred from the direction in which a particular ray happens to travel.

use super::*;
use fs_evidence::NumericalKind;
use fs_geom::{Aabb, TraceStepClaim};

#[derive(Default)]
pub(super) struct InitialMedia<'scene> {
    boundaries: Vec<Boundary<'scene>>,
    fixed_eye: OnceLock<InitialMediumIds>,
    animated: bool,
}

#[derive(Clone, Copy)]
struct InitialMediumIds {
    boundaries: [usize; MAX_MEDIUM_STACK_DEPTH],
    len: usize,
}

impl InitialMediumIds {
    fn materialize(self, scene: &Scene) -> Result<MediumStack<'_>, TracerError> {
        let mut stack = MediumStack::new();
        for &boundary_primitive in &self.boundaries[..self.len] {
            let Material::Dielectric { glass, .. } = &scene.primitives[boundary_primitive].material
            else {
                return Err(TracerError::InvalidInput);
            };
            stack.push(MediumEntry {
                boundary_primitive,
                glass,
            })?;
        }
        Ok(stack)
    }
}

struct Boundary<'scene> {
    primitive: usize,
    geometry: Geometry<'scene>,
    bounds: Aabb,
    placement: Placement<'scene>,
}

enum Geometry<'scene> {
    Mesh { mesh: &'scene TriMesh, closed: bool },
    Chart(&'scene dyn Chart),
}

enum Placement<'scene> {
    Fixed(RigidTransform),
    Animated(&'scene AnimatedGeometryInstance),
}

fn refusal(primitive: usize, reason: &'static str) -> TracerError {
    TracerError::InitialMedium {
        boundary_primitive: primitive,
        reason,
    }
}

impl<'scene> InitialMedia<'scene> {
    pub(super) fn try_new(scene: &'scene Scene, cx: &Cx<'_>) -> Result<Self, TracerError> {
        let mut result = Self::default();
        for (primitive, item) in scene.primitives.iter().enumerate() {
            cx.checkpoint()?;
            let Material::Dielectric { .. } = &item.material else {
                continue;
            };
            let (geometry, placement) = match &item.shape {
                Shape::Mesh(mesh) => (
                    Geometry::mesh(mesh, cx)?,
                    Placement::Fixed(RigidTransform::identity()),
                ),
                Shape::Chart(chart) => (
                    Geometry::Chart(chart.as_ref()),
                    Placement::Fixed(RigidTransform::identity()),
                ),
                Shape::Instance(instance) => (
                    Geometry::shared(instance.geometry(), cx)?,
                    Placement::Fixed(instance.transform()),
                ),
                Shape::AnimatedInstance(instance) => {
                    result.animated = true;
                    (
                        Geometry::shared(instance.geometry(), cx)?,
                        Placement::Animated(instance),
                    )
                }
            };
            let bounds = geometry.bounds();
            result
                .boundaries
                .try_reserve(1)
                .map_err(|_| TracerError::InvalidInput)?;
            result.boundaries.push(Boundary {
                primitive,
                geometry,
                bounds,
                placement,
            });
        }
        Ok(result)
    }

    /// Query one original boundary's closure and interior ownership without
    /// importing unrelated media at this diagnostic point into a path stack.
    #[cfg(feature = "photometry")]
    pub(super) fn boundary_contains(
        &self,
        cx: &Cx<'_>,
        primitive: usize,
        origin: Point3,
        time: Option<&PathTime>,
    ) -> Result<bool, TracerError> {
        cx.checkpoint()?;
        let boundary = self
            .boundaries
            .iter()
            .find(|boundary| boundary.primitive == primitive)
            .ok_or_else(|| {
                refusal(
                    primitive,
                    "selected primitive has no dielectric medium boundary",
                )
            })?;
        let local = boundary
            .placement(cx, time)?
            .inverse()?
            .transform_point(origin);
        boundary.geometry.contains(cx, local, primitive)
    }

    pub(super) fn stack<'path>(
        &self,
        scene: &'path Scene,
        cx: &Cx<'_>,
        origin: Point3,
        time: Option<&PathTime>,
    ) -> Result<MediumStack<'path>, TracerError> {
        cx.checkpoint()?;
        if self.boundaries.is_empty() {
            return Ok(MediumStack::new());
        }
        let fixed_eye = !self.animated && origin == scene.camera.eye;
        if fixed_eye && let Some(stack) = self.fixed_eye.get() {
            return stack.materialize(scene);
        }
        if [origin.x, origin.y, origin.z]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(TracerError::InvalidInput);
        }
        let mut containing: [Option<(&Boundary<'scene>, RigidTransform)>; MAX_MEDIUM_STACK_DEPTH] =
            [None; MAX_MEDIUM_STACK_DEPTH];
        let mut count = 0;
        for boundary in &self.boundaries {
            cx.checkpoint()?;
            let placement = boundary.placement(cx, time)?;
            let local = placement.inverse()?.transform_point(origin);
            if !boundary.bounds.contains(local) {
                continue;
            }
            if boundary.geometry.contains(cx, local, boundary.primitive)? {
                let slot = containing
                    .get_mut(count)
                    .ok_or(TracerError::MediumStackOverflow)?;
                *slot = Some((boundary, placement));
                count += 1;
            }
        }
        // Rank by proven enclosing supports, not nearest surface distance.
        // A conservative witness can refuse a valid complicated nesting; it
        // must never silently choose a medium for an unsupported overlap.
        let mut ranks = [0usize; MAX_MEDIUM_STACK_DEPTH];
        for inner in 0..count {
            for outer in 0..inner {
                let a = containing[inner].expect("occupied launch prefix");
                let b = containing[outer].expect("occupied launch prefix");
                if b.0.encloses(cx, b.1, a.0, a.1)? {
                    ranks[inner] += 1;
                } else if a.0.encloses(cx, a.1, b.0, b.1)? {
                    ranks[outer] += 1;
                } else {
                    return Err(refusal(
                        a.0.primitive,
                        "containing boundaries overlap or lack a strict nesting witness",
                    ));
                }
            }
        }
        let mut stack = InitialMediumIds {
            boundaries: [0; MAX_MEDIUM_STACK_DEPTH],
            len: count,
        };
        for rank in 0..count {
            let index = ranks[..count]
                .iter()
                .position(|&value| value == rank)
                .ok_or_else(|| {
                    refusal(
                        containing[0].unwrap().0.primitive,
                        "inconsistent containing-boundary order",
                    )
                })?;
            stack.boundaries[rank] = containing[index].unwrap().0.primitive;
        }
        if fixed_eye {
            let _ = self.fixed_eye.set(stack);
        }
        stack.materialize(scene)
    }
}

impl Boundary<'_> {
    fn placement(
        &self,
        cx: &Cx<'_>,
        time: Option<&PathTime>,
    ) -> Result<RigidTransform, TracerError> {
        match self.placement {
            Placement::Fixed(transform) => Ok(transform),
            Placement::Animated(instance) => {
                let time = time.ok_or(TracerError::MissingRayTime)?;
                if let Some(cached) = time
                    .cached_animated
                    .iter()
                    .flatten()
                    .find(|entry| entry.primitive_index == self.primitive)
                {
                    return Ok(cached.instance.transform());
                }
                Ok(instance
                    .instance_at(cx, time.interval.time_at(time.normalized))?
                    .transform())
            }
        }
    }

    fn encloses(
        &self,
        cx: &Cx<'_>,
        outer_pose: RigidTransform,
        inner: &Boundary<'_>,
        inner_pose: RigidTransform,
    ) -> Result<bool, TracerError> {
        let to_outer = outer_pose.inverse()?;
        let corners = box_corners(inner.bounds)
            .map(|point| to_outer.transform_point(inner_pose.transform_point(point)));
        if corners
            .iter()
            .any(|point| [point.x, point.y, point.z].iter().any(|v| !v.is_finite()))
        {
            return Ok(false);
        }
        match self.geometry {
            Geometry::Mesh { mesh, closed: true } => {
                // The intersection of all inward half-spaces is contained in
                // a closed, consistently oriented mesh. It is conservative
                // for concave hosts, and exact for a convex polyhedral host.
                for (index, triangle) in mesh.triangles.iter().enumerate() {
                    if index % 128 == 0 {
                        cx.checkpoint()?;
                    }
                    let [a, b, c] = triangle.map(|i| point(mesh.vertices[i as usize]));
                    let normal = cross(b.delta_from(a), c.delta_from(a));
                    for corner in corners {
                        let delta = corner.delta_from(a);
                        let products = [normal.x * delta.x, normal.y * delta.y, normal.z * delta.z];
                        let distance = products.iter().sum::<f64>();
                        let roundoff =
                            64.0 * f64::EPSILON * products.iter().map(|v| v.abs()).sum::<f64>();
                        if !distance.is_finite() || distance >= -roundoff {
                            return Ok(false);
                        }
                    }
                }
                Ok(true)
            }
            Geometry::Chart(chart) => {
                let center = Point3::new(
                    corners.iter().map(|p| p.x / 8.0).sum(),
                    corners.iter().map(|p| p.y / 8.0).sum(),
                    corners.iter().map(|p| p.z / 8.0).sum(),
                );
                let radius = corners
                    .iter()
                    .map(|p| p.delta_from(center).norm())
                    .fold(0.0, f64::max);
                let sample = chart.eval(center, cx);
                let enclosure = chart.trace_value_enclosure(center, &sample, cx);
                let lipschitz = match chart.trace_step_claim() {
                    TraceStepClaim::ExactDistance => 1.0,
                    TraceStepClaim::LipschitzImplicit => sample.lipschitz.unwrap_or(f64::INFINITY),
                    TraceStepClaim::NoClaim => return Ok(false),
                };
                Ok(valid_enclosure(enclosure)
                    && radius.is_finite()
                    && lipschitz.is_finite()
                    && lipschitz > 0.0
                    && enclosure.hi + (radius * lipschitz).next_up() < 0.0)
            }
            Geometry::Mesh { closed: false, .. } => Ok(false),
        }
    }
}

impl<'scene> Geometry<'scene> {
    fn shared(geometry: &'scene SharedGeometry, cx: &Cx<'_>) -> Result<Self, TracerError> {
        match geometry {
            SharedGeometry::Mesh(mesh) => Self::mesh(mesh.as_ref(), cx),
            SharedGeometry::Chart(chart) => Ok(Self::Chart(chart.as_ref())),
        }
    }

    fn mesh(mesh: &'scene TriMesh, cx: &Cx<'_>) -> Result<Self, TracerError> {
        let closed = mesh_closed(&mesh.vertices, &mesh.triangles, || {
            cx.checkpoint().map_err(Into::into)
        })?;
        Ok(Self::Mesh { mesh, closed })
    }

    fn bounds(&self) -> Aabb {
        match self {
            Self::Chart(chart) => chart.support(),
            Self::Mesh { mesh, .. } => {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for vertex in &mesh.vertices {
                    for axis in 0..3 {
                        lo[axis] = lo[axis].min(vertex[axis]);
                        hi[axis] = hi[axis].max(vertex[axis]);
                    }
                }
                Aabb {
                    min: point(lo),
                    max: point(hi),
                }
            }
        }
    }

    fn contains(&self, cx: &Cx<'_>, origin: Point3, primitive: usize) -> Result<bool, TracerError> {
        match self {
            Self::Chart(chart) => {
                let sample = chart.eval(origin, cx);
                let enclosure = chart.trace_value_enclosure(origin, &sample, cx);
                cx.checkpoint()?;
                if !valid_enclosure(enclosure) {
                    return Err(refusal(
                        primitive,
                        "chart has no finite certified sign at the launch point",
                    ));
                }
                if enclosure.hi < 0.0 {
                    Ok(true)
                } else if enclosure.lo > 0.0 {
                    Ok(false)
                } else {
                    Err(refusal(
                        primitive,
                        "launch point is on or unresolved against a dielectric boundary",
                    ))
                }
            }
            Self::Mesh { closed: false, .. } => Err(refusal(
                primitive,
                "launch point intersects the bounds of an open or inconsistently oriented dielectric mesh",
            )),
            Self::Mesh { mesh, closed: true } => {
                mesh_contains(&mesh.vertices, &mesh.triangles, origin, primitive, || {
                    cx.checkpoint().map_err(Into::into)
                })
            }
        }
    }
}

fn valid_enclosure(value: fs_evidence::NumericalCertificate) -> bool {
    matches!(value.kind, NumericalKind::Exact | NumericalKind::Enclosure)
        && value.lo.is_finite()
        && value.hi.is_finite()
        && value.lo <= value.hi
}

fn point(value: [f64; 3]) -> Point3 {
    Point3::new(value[0], value[1], value[2])
}

fn box_corners(bounds: Aabb) -> [Point3; 8] {
    std::array::from_fn(|i| {
        Point3::new(
            if i & 1 == 0 {
                bounds.min.x
            } else {
                bounds.max.x
            },
            if i & 2 == 0 {
                bounds.min.y
            } else {
                bounds.max.y
            },
            if i & 4 == 0 {
                bounds.min.z
            } else {
                bounds.max.z
            },
        )
    })
}

const UNIT_ROUNDOFF: f64 = f64::EPSILON * 0.5;
/// `det::atan2`'s declared budget in ulps (fs-math CONTRACT.md: `atan` 4 plus 1 for the composition; externally
/// audited by the MPFR lane). An audited budget, not a proof: it is the one assumed stage of the bound below.
const ATAN2_ULPS: f64 = 5.0;
/// The compensated sum's error `(2u + O(N u^2)) sum |Omega_i|` (Higham 2nd ed. (4.8); Goldberg 1991, Theorem 8) is taken
/// as `3u sum |Omega_i|` for at most this many triangles. The sources state the second-order term only as `O(N u^2)`, so
/// the range is stated, not its constant: `3u` covers `2u + c N u^2` whenever `c N u <= 1`, which holds for any `c <= 2^23`
/// at `N <= 2^30`. Above it the plain recursive bound `gamma_{N-1} sum |Omega_i|` (Higham section 4.2) is used.
const COMPENSATED_SUM_MAX_TERMS: f64 = 1_073_741_824.0;
/// A winding whose bound reaches this is not resolved at all: nearest-integer rounding needs it below 1/2, and the
/// bound is kept to half of that.
const WINDING_BOUND_LIMIT: f64 = 0.25;

/// `k u / (1 - k u)` for u = 2^-53 (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., Lemma 3.1:
/// `prod (1 + d_i)^(+-1) = 1 + theta_k`, `|theta_k| <= gamma_k`). The caller keeps `k u < 1`.
fn gamma(k: f64) -> f64 {
    k * UNIT_ROUNDOFF / (1.0 - k * UNIT_ROUNDOFF)
}

/// Absolute error `e0` of the computed pair `(D, N)` of one triangle, the same for every triangle.
///
/// A triangle's solid angle is `Omega = 2 atan2(N, D)` with `N = a.(b x c)` and `D = 1 + a.b + b.c + c.a` for the unit
/// vectors a, b, c towards its corners (Van Oosterom & Strackee, IEEE TBME 30(2), 1983). Near a face edge `N` and `D`
/// both go to 0, the origin being almost in the triangle's plane, and the conditioning of the angle is `1 / |(D, N)|`, so
/// the bound is computed per triangle from its own `(D, N)`. A fixed multiple of epsilon cannot serve.
///
/// With u = 2^-53 and gamma_k as above:
/// 1. A computed unit vector. `delta = vertex - origin` has relative error u per component; `|delta|` is
///    `sqrt(x*x + y*y + z*z)` summed left to right, relative gamma_3 / 2 + u; the reciprocal and the product add u each.
///    The computed vector is within `eps = gamma_6` of the exact unit vector (5.5 u) and has norm at most `m = 1 + eps`.
/// 2. `N` and `D` of the exact unit vectors against those of the computed ones. `N` is trilinear and `D` is bilinear
///    in the vectors, so the difference telescopes: `|dN| <= eps (|b x c| + |a x c| + |a x b|) <= 3 eps m^2` and
///    `|dD| <= eps (|b| + |a| + |c| + |b| + |c| + |a|) <= 6 eps m`.
/// 3. Rounding in `N`. Each component of `fl(b x c)` is within `gamma_2 (|b_i c_j| + |b_j c_i|) <= gamma_2 |b| |c|`
///    (Cauchy-Schwarz), so the vector is within `sqrt(3) gamma_2 m^2`; the dot with `a` adds `m` times that and its own
///    `gamma_3 |a| |w|` (Higham 3.5): `m^3 (sqrt(3) gamma_2 + gamma_3 (1 + sqrt(3) gamma_2))`.
/// 4. Rounding in `D`. Three dots, each within `gamma_3 m^2`, and the sum `((1 + ab) + bc) + ca`, whose three
///    roundings are at most `u (|s_1| + |s_2| + |s_3|) <= u (3 + 6 m^2) <= 9 u m^2`.
///
/// The components of the pair are `e0 = hypot(3 eps m^2 + [3.], 6 eps m + [4.])`.
fn triangle_pair_error() -> f64 {
    let eps = gamma(6.0);
    let m = 1.0 + eps;
    let root3 = 3.0_f64.sqrt();
    let n_error = 3.0 * eps * m * m
        + m.powi(3) * (root3 * gamma(2.0) + gamma(3.0) * (1.0 + root3 * gamma(2.0)));
    let d_error = 6.0 * eps * m + 3.0 * gamma(3.0) * m * m + 9.0 * UNIT_ROUNDOFF * m * m;
    n_error.hypot(d_error)
}

/// Oriented solid-angle winding of a closed mesh at `origin` and a bound `E` on its rounding error.
///
/// `E` assembles, over the triangles, `Omega_i`'s error and divides by 4 pi. The computed pair `(D, N)` lies within
/// `e0` of the exact pair (`triangle_pair_error`); the angle between two vectors of the plane, one of them known to lie in
/// the disc of radius `e0` about the other of length `|(D, N)|`, is at most `asin(rho)` with `rho = e0 / |(D, N)|`
/// (rigorous: no first-order expansion), and `asin(rho) <= rho / sqrt(1 - rho^2)`. So `Omega_i` is within
/// `2 rho / sqrt(1 - rho^2)` for `rho < 1` (otherwise `E` is infinite), plus the `det::atan2` budget `5 ulp(Omega_i)
/// <= 5 * 2^-52 |Omega_i|`. The sum is compensated (Kahan), error `3u sum |Omega_i|` up to `COMPENSATED_SUM_MAX_TERMS`
/// terms and `gamma_{N-1} sum |Omega_i|` beyond. The division by `4 pi` (the f64 constant is within 0.35 u of it) adds `2u |w|`. Evaluating
/// `E` itself costs at most `64 + N` roundings of relative size u on positive terms, `gamma_{64 + N}`.
///
/// `E < 1/2` is what nearest-integer rounding of the exact integer winding needs; `WINDING_BOUND_LIMIT` keeps half of
/// that for the one assumed stage.
fn mesh_winding(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    origin: Point3,
    primitive: usize,
    mut checkpoint: impl FnMut() -> Result<(), TracerError>,
) -> Result<(f64, f64), TracerError> {
    let pair_error = triangle_pair_error();
    let mut sum = 0.0;
    let mut correction = 0.0;
    let mut angle_error = 0.0;
    let mut magnitude = 0.0;
    for (index, triangle) in triangles.iter().enumerate() {
        if index % 128 == 0 {
            checkpoint()?;
        }
        let mut directions = [Vec3::new(0.0, 0.0, 0.0); 3];
        for (slot, &vertex) in directions.iter_mut().zip(triangle) {
            let delta = point(vertices[vertex as usize]).delta_from(origin);
            let length = delta.norm();
            if !length.is_finite() || length <= 0.0 {
                return Err(refusal(
                    primitive,
                    "launch point touches a dielectric vertex or has unresolvable scale",
                ));
            }
            *slot = delta.scale(1.0 / length);
        }
        let [a, b, c] = directions;
        // Van Oosterom & Strackee, IEEE TBME 30(2), 1983,
        // doi:10.1109/TBME.1983.325207. Unit vectors avoid cubic overflow.
        let numerator = a.dot(cross(b, c));
        let denominator = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
        if numerator.abs() <= 64.0 * f64::EPSILON && denominator <= 64.0 * f64::EPSILON {
            return Err(refusal(
                primitive,
                "launch point is on or unresolved against a dielectric triangle",
            ));
        }
        let angle = 2.0 * det::atan2(numerator, denominator);
        let rho = pair_error / (numerator * numerator + denominator * denominator).sqrt();
        angle_error += if rho < 1.0 {
            2.0 * rho / (1.0 - rho * rho).sqrt() + ATAN2_ULPS * 2.0_f64.powi(-52) * angle.abs()
        } else {
            f64::INFINITY
        };
        magnitude += angle.abs();
        let adjusted = angle - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    let winding = sum / (4.0 * PI);
    let count = triangles.len() as f64;
    let sum_coefficient = if count <= COMPENSATED_SUM_MAX_TERMS {
        3.0 * UNIT_ROUNDOFF
    } else {
        gamma(count - 1.0)
    };
    let rounding = (angle_error + sum_coefficient * magnitude) / (4.0 * PI)
        + 2.0 * UNIT_ROUNDOFF * winding.abs();
    Ok((winding, rounding * (1.0 + gamma(64.0 + count))))
}

/// Whether `origin` is inside the closed mesh. The exact winding of a closed, consistently oriented, non-intersecting
/// mesh is 1 inside and 0 outside; the computed one is accepted only within its derived bound `E` of that integer.
///
/// - `E` not finite or at least 1/4: refuse, numerical resolution;
/// - the nearest integer `k` is not 0 or 1: refuse, reversed or intersecting;
/// - `k` is 0 or 1 and `|w - k| > E`: refuse, unresolved winding;
/// - otherwise the origin is inside exactly when `k == 1`.
fn mesh_contains(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    origin: Point3,
    primitive: usize,
    checkpoint: impl FnMut() -> Result<(), TracerError>,
) -> Result<bool, TracerError> {
    let (winding, bound) = mesh_winding(vertices, triangles, origin, primitive, checkpoint)?;
    if !winding.is_finite() || !(bound < WINDING_BOUND_LIMIT) {
        return Err(refusal(
            primitive,
            "dielectric winding evaluation exceeds its numerical resolution",
        ));
    }
    let nearest = winding.round();
    if nearest != 0.0 && nearest != 1.0 {
        return Err(refusal(
            primitive,
            "dielectric mesh has reversed or intersecting winding at the launch point",
        ));
    }
    if (winding - nearest).abs() > bound {
        return Err(refusal(
            primitive,
            "dielectric mesh has an unresolved winding at the launch point",
        ));
    }
    Ok(nearest == 1.0)
}

fn mesh_closed(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    mut checkpoint: impl FnMut() -> Result<(), TracerError>,
) -> Result<bool, TracerError> {
    let capacity = triangles
        .len()
        .checked_mul(3)
        .ok_or(TracerError::InvalidInput)?;
    let mut edges = Vec::new();
    edges
        .try_reserve_exact(capacity)
        .map_err(|_| TracerError::InvalidInput)?;
    let mut closed = !triangles.is_empty();
    for (index, &[a, b, c]) in triangles.iter().enumerate() {
        if index % 128 == 0 {
            checkpoint()?;
        }
        if [a, b, c].iter().any(|&i| i as usize >= vertices.len()) || a == b || b == c || c == a {
            closed = false;
            break;
        }
        let [pa, pb, pc] = [a, b, c].map(|i| point(vertices[i as usize]));
        let area = cross(pb.delta_from(pa), pc.delta_from(pa)).norm();
        if !area.is_finite() || area == 0.0 {
            closed = false;
            break;
        }
        for (from, to) in [(a, b), (b, c), (c, a)] {
            edges.push((from.min(to), from.max(to), from < to));
        }
    }
    edges.sort_unstable();
    if edges.len() % 2 != 0 {
        closed = false;
    }
    for edge in edges.chunks_exact(2) {
        if edge[0].0 != edge[1].0 || edge[0].1 != edge[1].1 || edge[0].2 == edge[1].2 {
            closed = false;
            break;
        }
    }
    // Exactly two faces per undirected edge; four paired faces cannot
    // masquerade as two independent manifold edges with the same indices.
    if edges
        .windows(3)
        .any(|e| e[0].0 == e[2].0 && e[0].1 == e[2].1)
    {
        closed = false;
    }
    checkpoint()?;
    Ok(closed)
}

#[cfg(test)]
mod tests {
    use super::*;

    const VERTICES: [[f64; 3]; 8] = [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ];
    const TRIANGLES: [[u32; 3]; 12] = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];

    #[test]
    fn g0_complete_winding_classifies_interior_and_exterior() {
        assert!(mesh_closed(&VERTICES, &TRIANGLES, || Ok(())).unwrap());
        for p in [Point3::new(0.0, 0.0, 0.0), Point3::new(0.6, -0.3, 0.9)] {
            assert!(mesh_contains(&VERTICES, &TRIANGLES, p, 7, || Ok(())).unwrap());
        }
        for p in [Point3::new(0.0, 0.0, 2.0), Point3::new(2.0, 1.0, 0.0)] {
            assert!(!mesh_contains(&VERTICES, &TRIANGLES, p, 7, || Ok(())).unwrap());
        }
    }

    #[test]
    fn g0_open_reversed_and_nonmanifold_boundaries_do_not_mint_media() {
        assert!(!mesh_closed(&VERTICES, &TRIANGLES[..11], || Ok(())).unwrap());
        let mut reversed = TRIANGLES;
        for triangle in &mut reversed {
            triangle.swap(1, 2);
        }
        assert!(
            mesh_contains(&VERTICES, &reversed, Point3::new(0.0, 0.0, 0.0), 7, || Ok(
                ()
            ))
            .is_err()
        );
        let duplicated = [TRIANGLES.as_slice(), TRIANGLES.as_slice()].concat();
        assert!(!mesh_closed(&VERTICES, &duplicated, || Ok(())).unwrap());
        reversed = TRIANGLES;
        reversed[0].swap(1, 2);
        assert!(!mesh_closed(&VERTICES, &reversed, || Ok(())).unwrap());
    }

    #[test]
    fn g0_launch_on_face_edge_or_vertex_refuses() {
        for p in [
            Point3::new(0.2, 0.3, 1.0),
            Point3::new(1.0, 0.3, 1.0),
            Point3::new(1.0, 1.0, 1.0),
        ] {
            assert!(mesh_contains(&VERTICES, &TRIANGLES, p, 7, || Ok(())).is_err());
        }
    }

    // The bdpt fixture's slab: x and y in [-2, 2], z from `bottom` to 1, the cube's triangles.
    fn slab(bottom: f64) -> [[f64; 3]; 8] {
        VERTICES.map(|v| {
            [
                2.0 * v[0],
                2.0 * v[1],
                if v[2] < 0.0 { bottom } else { 1.0 },
            ]
        })
    }

    #[test]
    fn g1_a_launch_point_micrometres_from_a_face_edge_is_resolved_inside_and_outside() {
        for bottom in [-5.0e-5, -0.5e-6] {
            let vertices = slab(bottom);
            // Over the bottom face and its diagonal edge, 4e-7 m above it in the thin case.
            for (x, y) in [(0.0, 0.0), (0.3, 0.3), (-0.7, 0.2)] {
                let origin = Point3::new(x, y, -0.1e-6);
                let (winding, bound) =
                    mesh_winding(&vertices, &TRIANGLES, origin, 7, || Ok(())).unwrap();
                assert!(bound < 1.0e-7, "{bound}");
                assert!((winding - 1.0).abs() <= bound, "{winding} {bound}");
                assert!(mesh_contains(&vertices, &TRIANGLES, origin, 7, || Ok(())).unwrap());
            }
            // The same distance below the bottom face: outside, with the same conditioning.
            let outside = Point3::new(0.0, 0.0, bottom - 4.0e-7);
            let (winding, bound) =
                mesh_winding(&vertices, &TRIANGLES, outside, 7, || Ok(())).unwrap();
            assert!(winding.abs() <= bound, "{winding} {bound}");
            assert!(!mesh_contains(&vertices, &TRIANGLES, outside, 7, || Ok(())).unwrap());
        }
        // The bound is the conditioning, not a constant: it grows as the point nears the face, and the
        // computed winding stays inside it all the way to the on-surface guard.
        let vertices = slab(-0.5e-6);
        let mut previous = 0.0;
        for height in [1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12, 1.0e-13] {
            let origin = Point3::new(0.0, 0.0, -0.5e-6 + height);
            let (winding, bound) =
                mesh_winding(&vertices, &TRIANGLES, origin, 7, || Ok(())).unwrap();
            assert!(
                bound > previous && bound < WINDING_BOUND_LIMIT,
                "{height} {bound}"
            );
            assert!(
                (winding - 1.0).abs() <= bound,
                "{height}: {winding} {bound}"
            );
            assert!(mesh_contains(&vertices, &TRIANGLES, origin, 7, || Ok(())).unwrap());
            previous = bound;
        }
    }

    #[test]
    fn g0_the_derived_bound_still_refuses_reversed_doubled_open_and_boundary_launches() {
        let vertices = slab(-0.5e-6);
        let inside = Point3::new(0.0, 0.0, -0.1e-6);
        let mut reversed = TRIANGLES;
        for triangle in &mut reversed {
            triangle.swap(1, 2);
        }
        let refused = |result: Result<bool, TracerError>, needle: &str| match result {
            Err(TracerError::InitialMedium { reason, .. }) => {
                assert!(reason.contains(needle), "{reason}");
            }
            other => panic!("expected a refusal containing {needle:?}, got {other:?}"),
        };
        // A reversed slab: w = -1.
        refused(
            mesh_contains(&vertices, &reversed, inside, 7, || Ok(())),
            "reversed or intersecting",
        );
        // A doubly covered mesh: w = 2.
        let doubled = [TRIANGLES.as_slice(), TRIANGLES.as_slice()].concat();
        refused(
            mesh_contains(&vertices, &doubled, inside, 7, || Ok(())),
            "reversed or intersecting",
        );
        // An open mesh (the wall triangle [3, 4, 7] removed), 0.01 in front of the gap: w is about 1/2.
        let gap = Point3::new(-1.99, 1.0, 0.9);
        let (winding, bound) =
            mesh_winding(&vertices, &TRIANGLES[..11], gap, 7, || Ok(())).unwrap();
        assert!(
            winding > 0.2 && winding < 0.8 && bound < 1.0e-6,
            "{winding} {bound}"
        );
        refused(
            mesh_contains(&vertices, &TRIANGLES[..11], gap, 7, || Ok(())),
            "unresolved winding",
        );
        // On a face: the on-surface guard, with its own reason.
        refused(
            mesh_contains(&vertices, &TRIANGLES, Point3::new(0.2, 0.3, 1.0), 7, || {
                Ok(())
            }),
            "on or unresolved against a dielectric triangle",
        );
        // So close to the bottom face that the guard refuses before any bound is asked.
        refused(
            mesh_contains(
                &vertices,
                &TRIANGLES,
                Point3::new(0.0, 0.0, -0.5e-6 + 1.0e-15),
                7,
                || Ok(()),
            ),
            "on or unresolved against a dielectric triangle",
        );
        // 1 m outside: not an error.
        assert!(
            !mesh_contains(&vertices, &TRIANGLES, Point3::new(0.0, 0.0, 2.0), 7, || Ok(
                ()
            ))
            .unwrap()
        );
    }

    #[test]
    fn g3_winding_is_invariant_to_translation_and_uniform_scaling() {
        for scale in [1.0e-8, 1.0, 1.0e8] {
            let vertices = VERTICES.map(|v| {
                [
                    scale * (v[0] + 2.0),
                    scale * (v[1] - 4.0),
                    scale * (v[2] + 3.0),
                ]
            });
            let origin = Point3::new(scale * 2.0, scale * -4.0, scale * 3.0);
            assert!(mesh_contains(&vertices, &TRIANGLES, origin, 7, || Ok(())).unwrap());
        }
    }

    #[test]
    fn g4_launch_topology_and_membership_poll_cancellation() {
        assert_eq!(
            mesh_closed(&VERTICES, &TRIANGLES, || Err(TracerError::Cancelled)),
            Err(TracerError::Cancelled)
        );
        assert_eq!(
            mesh_contains(
                &VERTICES,
                &TRIANGLES,
                Point3::new(0.0, 0.0, 0.0),
                7,
                || Err(TracerError::Cancelled)
            ),
            Err(TracerError::Cancelled)
        );
    }
}
