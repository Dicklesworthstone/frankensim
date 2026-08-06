//! Shared immutable geometry placed by validated proper-rigid transforms.
//!
//! A render instance is visualization state. Its transform never changes the
//! underlying chart identity, mass properties, or contact geometry.

use core::fmt;
use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_exec::{Cancelled, Cx};
use fs_geom::{Chart, Point3, Vec3};

use crate::charts::{Hit, Ray, TraceAudit, TraceTermination, TriMesh, sphere_trace};

const TRANSFORM_DOMAIN: &str = "org.frankensim.render.rigid-transform.v1";
const INSTANCE_DOMAIN: &str = "org.frankensim.render.geometry-instance.v1";
const SCENE_DOMAIN: &str = "org.frankensim.render.instance-scene.v1";
const UNIT_QUATERNION_TOLERANCE: f64 = 1.0e-12;

/// Proper-rigid body-to-world transform: a unit quaternion plus translation.
/// Scale, shear, and reflection are deliberately not representable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidTransform {
    rotation_xyzw: [f64; 4],
    translation_m: [f64; 3],
}

impl RigidTransform {
    /// Identity placement.
    #[must_use]
    pub const fn identity() -> Self {
        Self {
            rotation_xyzw: [0.0, 0.0, 0.0, 1.0],
            translation_m: [0.0; 3],
        }
    }

    /// Admit a finite near-unit quaternion and finite SI translation. Accepted
    /// quaternions are normalized and sign-canonicalized so `q` and `-q` have
    /// one representation.
    pub fn try_new(
        rotation_xyzw: [f64; 4],
        translation_m: [f64; 3],
    ) -> Result<Self, InstanceError> {
        if rotation_xyzw
            .iter()
            .chain(translation_m.iter())
            .any(|value| !value.is_finite())
        {
            return Err(InstanceError::InvalidTransform);
        }
        let norm_squared = rotation_xyzw.iter().map(|value| value * value).sum::<f64>();
        if !norm_squared.is_finite() || (norm_squared - 1.0).abs() > UNIT_QUATERNION_TOLERANCE {
            return Err(InstanceError::InvalidTransform);
        }
        Ok(Self::normalized(rotation_xyzw, translation_m))
    }

    fn normalized(mut rotation_xyzw: [f64; 4], mut translation_m: [f64; 3]) -> Self {
        let inverse_norm = 1.0
            / rotation_xyzw
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
        for value in &mut rotation_xyzw {
            *value *= inverse_norm;
        }
        if rotation_xyzw
            .iter()
            .find(|value| **value != 0.0)
            .is_some_and(|value| *value < 0.0)
        {
            for value in &mut rotation_xyzw {
                *value = -*value;
            }
        }
        for value in rotation_xyzw.iter_mut().chain(&mut translation_m) {
            if *value == 0.0 {
                *value = 0.0;
            }
        }
        Self {
            rotation_xyzw,
            translation_m,
        }
    }

    /// Canonical body-to-world quaternion in `(x, y, z, w)` order.
    #[must_use]
    pub const fn rotation_xyzw(self) -> [f64; 4] {
        self.rotation_xyzw
    }

    /// Body-origin translation in world metres.
    #[must_use]
    pub const fn translation_m(self) -> [f64; 3] {
        self.translation_m
    }

    /// Place a local point in world coordinates.
    #[must_use]
    pub fn transform_point(self, point: Point3) -> Point3 {
        let rotated = self.transform_vector(Vec3::new(point.x, point.y, point.z));
        Point3::new(
            rotated.x + self.translation_m[0],
            rotated.y + self.translation_m[1],
            rotated.z + self.translation_m[2],
        )
    }

    /// Rotate a local vector into world coordinates.
    #[must_use]
    pub fn transform_vector(self, vector: Vec3) -> Vec3 {
        rotate(self.rotation_xyzw, vector)
    }

    /// Inverse world-to-body transform.
    pub fn inverse(self) -> Result<Self, InstanceError> {
        let [x, y, z, w] = self.rotation_xyzw;
        let rotation = [-x, -y, -z, w];
        let translated = rotate(
            rotation,
            Vec3::new(
                -self.translation_m[0],
                -self.translation_m[1],
                -self.translation_m[2],
            ),
        );
        Self::try_new(rotation, [translated.x, translated.y, translated.z])
    }

    /// Compose `self` after `local_to_parent`.
    pub fn compose(self, local_to_parent: Self) -> Result<Self, InstanceError> {
        let rotation = multiply_quaternions(self.rotation_xyzw, local_to_parent.rotation_xyzw);
        let child_translation = self.transform_vector(Vec3::new(
            local_to_parent.translation_m[0],
            local_to_parent.translation_m[1],
            local_to_parent.translation_m[2],
        ));
        Self::try_new(
            rotation,
            [
                child_translation.x + self.translation_m[0],
                child_translation.y + self.translation_m[1],
                child_translation.z + self.translation_m[2],
            ],
        )
    }

    /// Content identity of this canonical placement alone.
    #[must_use]
    pub fn content_identity(self) -> ContentHash {
        let mut bytes = Vec::with_capacity(56);
        for value in self.rotation_xyzw.into_iter().chain(self.translation_m) {
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        hash_domain(TRANSFORM_DOMAIN, &bytes)
    }

    fn world_to_local_ray(self, ray: &Ray) -> Result<Ray, InstanceError> {
        let inverse = self.inverse()?;
        let local = Ray {
            origin: inverse.transform_point(ray.origin),
            dir: inverse.transform_vector(ray.dir),
        };
        if [
            local.origin.x,
            local.origin.y,
            local.origin.z,
            local.dir.x,
            local.dir.y,
            local.dir.z,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            Ok(local)
        } else {
            Err(InstanceError::InvalidTransform)
        }
    }
}

/// Shared immutable local-space geometry.
#[derive(Clone)]
pub enum SharedGeometry {
    /// Certified chart, evaluated without conversion or tessellation.
    Chart(Arc<dyn Chart>),
    /// Native deterministic-BVH mesh.
    Mesh(Arc<TriMesh>),
}

impl SharedGeometry {
    /// Share a concrete chart.
    #[must_use]
    pub fn chart(chart: impl Chart + 'static) -> Self {
        Self::Chart(Arc::new(chart))
    }

    /// Share a mesh.
    #[must_use]
    pub fn mesh(mesh: TriMesh) -> Self {
        Self::Mesh(Arc::new(mesh))
    }

    /// Whether two handles point to the same allocation and backend kind.
    #[must_use]
    pub fn ptr_eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Chart(left), Self::Chart(right)) => Arc::ptr_eq(left, right),
            (Self::Mesh(left), Self::Mesh(right)) => Arc::ptr_eq(left, right),
            _ => false,
        }
    }
}

impl fmt::Debug for SharedGeometry {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Chart(_) => formatter.write_str("SharedGeometry::Chart(..)"),
            Self::Mesh(mesh) => formatter
                .debug_tuple("SharedGeometry::Mesh")
                .field(&mesh.bvh_fingerprint())
                .finish(),
        }
    }
}

/// Immutable geometry plus one object identity and body-to-world placement.
#[derive(Debug, Clone)]
pub struct GeometryInstance {
    object_id: u64,
    geometry_identity: ContentHash,
    geometry: SharedGeometry,
    transform: RigidTransform,
}

impl GeometryInstance {
    /// Bind a nonzero object ID and caller-supplied immutable-geometry identity.
    pub fn try_new(
        object_id: u64,
        geometry_identity: ContentHash,
        geometry: SharedGeometry,
        transform: RigidTransform,
    ) -> Result<Self, InstanceError> {
        if object_id == 0 {
            return Err(InstanceError::InvalidObjectId);
        }
        if geometry_identity.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(InstanceError::InvalidGeometryIdentity);
        }
        Ok(Self {
            object_id,
            geometry_identity,
            geometry,
            transform,
        })
    }

    /// Stable nonzero object ID used for tie-breaking and AOVs.
    #[must_use]
    pub const fn object_id(&self) -> u64 {
        self.object_id
    }

    /// Caller-supplied identity of the immutable local geometry.
    #[must_use]
    pub const fn geometry_identity(&self) -> ContentHash {
        self.geometry_identity
    }

    /// Shared local geometry handle.
    #[must_use]
    pub const fn geometry(&self) -> &SharedGeometry {
        &self.geometry
    }

    /// Body-to-world placement.
    #[must_use]
    pub const fn transform(&self) -> RigidTransform {
        self.transform
    }

    /// Frame identity binds object, immutable geometry, and placement.
    #[must_use]
    pub fn frame_identity(&self) -> ContentHash {
        let mut bytes = Vec::with_capacity(72);
        bytes.extend_from_slice(&self.object_id.to_le_bytes());
        bytes.extend_from_slice(self.geometry_identity.as_bytes());
        bytes.extend_from_slice(self.transform.content_identity().as_bytes());
        hash_domain(INSTANCE_DOMAIN, &bytes)
    }

    /// Intersect one placed object, preserving chart authority/audit.
    pub fn intersect(
        &self,
        cx: &Cx<'_>,
        world_ray: &Ray,
        t_max: f64,
        eps: f64,
    ) -> Result<Option<InstanceHit>, InstanceError> {
        cx.checkpoint()?;
        let ray_values = [
            world_ray.origin.x,
            world_ray.origin.y,
            world_ray.origin.z,
            world_ray.dir.x,
            world_ray.dir.y,
            world_ray.dir.z,
        ];
        let direction_scale = world_ray
            .dir
            .x
            .abs()
            .max(world_ray.dir.y.abs())
            .max(world_ray.dir.z.abs());
        if ray_values.iter().any(|value| !value.is_finite())
            || direction_scale == 0.0
            || !t_max.is_finite()
            || t_max <= 0.0
            || !eps.is_finite()
            || eps <= 0.0
        {
            return Err(InstanceError::InvalidIntersectionInput);
        }
        let local_ray = self.transform.world_to_local_ray(world_ray)?;
        let (hit, audit) = match &self.geometry {
            SharedGeometry::Chart(chart) => {
                let (hit, trace_audit) =
                    sphere_trace(chart.as_ref(), cx, &local_ray, t_max, eps, 1.0);
                if matches!(
                    trace_audit.termination,
                    TraceTermination::Hit
                        | TraceTermination::ResidualLimit
                        | TraceTermination::Miss
                ) && !trace_audit.certified
                {
                    return Err(InstanceError::UncertifiedTrace);
                }
                let hit = match trace_audit.termination {
                    TraceTermination::Cancelled => return Err(InstanceError::Cancelled),
                    TraceTermination::Miss => None,
                    TraceTermination::Hit => {
                        Some(hit.ok_or(InstanceError::BackendFailure(TraceTermination::Hit))?)
                    }
                    termination => return Err(InstanceError::BackendFailure(termination)),
                };
                (
                    hit.map(|hit| (hit, InstanceSurfaceFeature::ChartUnavailable)),
                    InstanceBackendAudit::Chart(trace_audit),
                )
            }
            SharedGeometry::Mesh(mesh) => {
                let hit = mesh
                    .intersect_surface_with_cx(cx, &local_ray)?
                    .map(|mesh_hit| {
                        (
                            mesh_hit.hit,
                            InstanceSurfaceFeature::MeshTriangle {
                                triangle_index: mesh_hit.triangle_index,
                                barycentric: mesh_hit.barycentric,
                            },
                        )
                    });
                (
                    hit,
                    InstanceBackendAudit::Mesh {
                        bvh_fingerprint: mesh.bvh_fingerprint(),
                    },
                )
            }
        };
        let Some((local_hit, surface_feature)) =
            hit.filter(|(hit, _)| hit.t > 0.0 && hit.t <= t_max)
        else {
            return Ok(None);
        };
        let world_hit = transform_hit(self.transform, local_hit)?;
        cx.checkpoint()?;
        Ok(Some(InstanceHit {
            object_id: self.object_id,
            geometry_identity: self.geometry_identity,
            frame_identity: self.frame_identity(),
            local_hit,
            surface_feature,
            hit: world_hit,
            backend_audit: audit,
        }))
    }
}

/// Backend authority retained by an instance hit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InstanceBackendAudit {
    /// Complete certified chart trace audit.
    Chart(TraceAudit),
    /// Mesh layout receipt used for traversal.
    Mesh {
        /// Deterministic scalar-BVH fingerprint.
        bvh_fingerprint: u64,
    },
}

/// Stable local feature witness retained by an instance hit.
///
/// A mesh triangle index is stable for one immutable ordered mesh artifact and
/// its barycentric coordinates reconstruct the same local point. Generic
/// charts do not currently expose a stable surface parameter, so their object
/// and material identities remain usable while temporal correspondence is
/// explicitly unavailable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InstanceSurfaceFeature {
    /// The chart backend has no admitted stable surface parameterization.
    ChartUnavailable,
    /// Original mesh triangle and its local barycentric point.
    MeshTriangle {
        /// Index into the immutable mesh's triangle array.
        triangle_index: u32,
        /// Weights ordered like that triangle's three vertex indices.
        barycentric: [f64; 3],
    },
}

/// Closest world-space hit with stable object and geometry identities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InstanceHit {
    /// Stable object ID.
    pub object_id: u64,
    /// Immutable local geometry identity.
    pub geometry_identity: ContentHash,
    /// Identity of object, geometry, and current transform.
    pub frame_identity: ContentHash,
    /// Complete local-space hit used to map the same rigid material point at
    /// another pose. Local geometric and shading normals remain distinct.
    pub local_hit: Hit,
    /// Backend-specific stable feature witness, or an explicit chart refusal.
    pub surface_feature: InstanceSurfaceFeature,
    /// World-space hit with unchanged ray parameter.
    pub hit: Hit,
    /// Authority/audit from the local backend.
    pub backend_audit: InstanceBackendAudit,
}

/// Deterministically ordered collection of collision-free object IDs.
#[derive(Debug, Clone)]
pub struct InstanceScene {
    instances: Vec<GeometryInstance>,
    identity: ContentHash,
}

impl InstanceScene {
    /// Validate IDs and sort by object ID so exact-distance ties are stable.
    pub fn try_new(mut instances: Vec<GeometryInstance>) -> Result<Self, InstanceError> {
        instances.sort_by_key(GeometryInstance::object_id);
        if instances
            .windows(2)
            .any(|pair| pair[0].object_id == pair[1].object_id)
        {
            return Err(InstanceError::DuplicateObjectId);
        }
        let mut hasher = DomainHasher::new(SCENE_DOMAIN);
        let instance_count =
            u64::try_from(instances.len()).map_err(|_| InstanceError::TooManyInstances)?;
        hasher.update(&instance_count.to_le_bytes());
        for instance in &instances {
            hasher.update(instance.frame_identity().as_bytes());
        }
        let identity = hasher.finalize();
        Ok(Self {
            instances,
            identity,
        })
    }

    /// Ordered instances.
    #[must_use]
    pub fn instances(&self) -> &[GeometryInstance] {
        &self.instances
    }

    /// Identity of the ordered frame placements.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Closest hit. Exact `t` ties select the lowest object ID.
    pub fn intersect(
        &self,
        cx: &Cx<'_>,
        ray: &Ray,
        t_max: f64,
        eps: f64,
    ) -> Result<Option<InstanceHit>, InstanceError> {
        let mut best: Option<InstanceHit> = None;
        for instance in &self.instances {
            cx.checkpoint()?;
            if let Some(candidate) = instance.intersect(cx, ray, t_max, eps)?
                && best
                    .as_ref()
                    .is_none_or(|current| candidate.hit.t < current.hit.t)
            {
                best = Some(candidate);
            }
        }
        Ok(best)
    }
}

/// Fail-closed transform, scene, and backend errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InstanceError {
    /// Non-finite or non-unit rigid transform.
    InvalidTransform,
    /// Object ID zero is reserved as invalid.
    InvalidObjectId,
    /// Immutable geometry identity is absent.
    InvalidGeometryIdentity,
    /// Two instances reuse one object ID.
    DuplicateObjectId,
    /// The instance count cannot be represented by the canonical identity.
    TooManyInstances,
    /// Ray or intersection limits were non-finite, zero, or negative.
    InvalidIntersectionInput,
    /// Execution was cancelled.
    Cancelled,
    /// Chart stopped without a certified hit or clean miss.
    BackendFailure(TraceTermination),
    /// Chart terminal state lacked its typed trace claim.
    UncertifiedTrace,
    /// Shading requires a finite geometric normal.
    MissingNormal,
    /// A backend produced non-finite world-space hit data.
    InvalidHit,
}

impl From<Cancelled> for InstanceError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl fmt::Display for InstanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "invalid render instance: {self:?}")
    }
}

impl std::error::Error for InstanceError {}

fn transform_hit(transform: RigidTransform, hit: Hit) -> Result<Hit, InstanceError> {
    let normal =
        transform_optional_unit(transform, hit.normal)?.ok_or(InstanceError::MissingNormal)?;
    let shading_normal = transform_optional_unit(transform, hit.shading_normal)?;
    let tangent_u = transform_optional_unit(transform, hit.tangent_u)?;
    let tangent_v = transform_optional_unit(transform, hit.tangent_v)?;
    let dp_du = transform_optional_finite(transform, hit.dp_du)?;
    let dp_dv = transform_optional_finite(transform, hit.dp_dv)?;
    let point = transform.transform_point(hit.point);
    if !hit.t.is_finite()
        || hit.t <= 0.0
        || ![point.x, point.y, point.z]
            .iter()
            .all(|value| value.is_finite())
    {
        return Err(InstanceError::InvalidHit);
    }
    Ok(Hit {
        t: hit.t,
        point,
        normal: Some(normal),
        shading_normal,
        tangent_u,
        tangent_v,
        dp_du,
        dp_dv,
        steps: hit.steps,
    })
}

fn transform_optional_unit(
    transform: RigidTransform,
    vector: Option<Vec3>,
) -> Result<Option<Vec3>, InstanceError> {
    vector
        .map(|value| normalize(transform.transform_vector(value)).ok_or(InstanceError::InvalidHit))
        .transpose()
}

fn transform_optional_finite(
    transform: RigidTransform,
    vector: Option<Vec3>,
) -> Result<Option<Vec3>, InstanceError> {
    vector
        .map(|value| {
            let transformed = transform.transform_vector(value);
            [transformed.x, transformed.y, transformed.z]
                .iter()
                .all(|component| component.is_finite())
                .then_some(transformed)
                .ok_or(InstanceError::InvalidHit)
        })
        .transpose()
}

fn normalize(vector: Vec3) -> Option<Vec3> {
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    let scaled = vector.scale(1.0 / scale);
    let norm = scaled.norm();
    (norm.is_finite() && norm > 0.0).then(|| scaled.scale(1.0 / norm))
}

fn rotate([x, y, z, w]: [f64; 4], vector: Vec3) -> Vec3 {
    let q = Vec3::new(x, y, z);
    let twice_cross = cross(q, vector).scale(2.0);
    let weighted_cross = twice_cross.scale(w);
    let double_cross = cross(q, twice_cross);
    Vec3::new(
        vector.x + weighted_cross.x + double_cross.x,
        vector.y + weighted_cross.y + double_cross.y,
        vector.z + weighted_cross.z + double_cross.z,
    )
}

fn cross(left: Vec3, right: Vec3) -> Vec3 {
    Vec3::new(
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    )
}

fn multiply_quaternions([ax, ay, az, aw]: [f64; 4], [bx, by, bz, bw]: [f64; 4]) -> [f64; 4] {
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}
