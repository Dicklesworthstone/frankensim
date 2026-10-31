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
            // Pure emitters do not scatter and their material is ignored by
            // the tracer. An emitter's decorative material is not a volume.
            if item.emission.is_some() {
                continue;
            }
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

fn mesh_contains(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    origin: Point3,
    primitive: usize,
    mut checkpoint: impl FnMut() -> Result<(), TracerError>,
) -> Result<bool, TracerError> {
    let mut sum = 0.0;
    let mut correction = 0.0;
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
        let adjusted = angle - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    let winding = sum / (4.0 * PI);
    let tolerance = 512.0 * f64::EPSILON * (triangles.len().max(1) as f64);
    if !winding.is_finite() || tolerance >= 0.01 {
        return Err(refusal(
            primitive,
            "dielectric winding evaluation exceeds its numerical resolution",
        ));
    }
    if winding.abs() <= tolerance {
        Ok(false)
    } else if (winding - 1.0).abs() <= tolerance {
        Ok(true)
    } else {
        Err(refusal(
            primitive,
            "dielectric mesh has reversed, intersecting or unresolved winding at the launch point",
        ))
    }
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
