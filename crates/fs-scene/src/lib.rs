//! fs-scene — the static environment a robot experiment moves through.
//!
//! # Why this crate exists
//!
//! Every robot experiment in this workspace had its own private notion of
//! "the world": the walking owner carried a `Vec<ObstacleBox>` and a bespoke
//! sphere/box test, the manipulation owner carried exactly one axis-aligned
//! keep-out box, and the surfaces the robots actually stand and work on — the
//! floor, the table — were not represented in any physics at all. They existed
//! only as meshes in the renderer.
//!
//! That gap is not academic. A renderer is free to draw a foundation slab, or
//! a work surface, wherever it likes; nothing checks it against the robot. The
//! failure mode is a robot standing *inside* the floor, which no collision
//! guard can catch because the floor was never a body.
//!
//! This crate is the shared answer. It owns:
//!
//! - [`SceneBody`], the geometry an environment is made of (yawed boxes and
//!   half-spaces);
//! - [`BodyRole`], the distinction between a body nothing may enter
//!   ([`BodyRole::KeepOut`]) and a surface things may rest on but not sink
//!   through ([`BodyRole::Support`]);
//! - [`StaticScene`], a validated collection of them;
//! - exact interpenetration queries for the primitive shapes, plus a
//!   [`ConvexSupportMap`](fs_query::ConvexSupportMap) bridge so an
//!   arbitrary convex collider can be handed to `fs-query`'s certified
//!   separation machinery.
//!
//! It knows nothing about any particular robot. A walking humanoid, a fixed
//! arm, and a mobile base all describe their colliders the same way and get
//! the same guarantee.
//!
//! # What "support" means
//!
//! A keep-out body is simple: any overlap is a violation. A support body is
//! not, because standing on the floor *requires* touching it, and compliant
//! contact models legitimately indent a few millimetres. So a support body
//! reports a violation only past a declared skin depth. That is what makes
//! "the robot is 27 cm inside the floor" a detectable, refusable state while
//! "the foot is 3 mm into the floor" is normal contact.
//!
//! # Determinism
//!
//! Every query here is a closed-form arithmetic evaluation over the body list
//! in declaration order. No iteration counts, no tolerances beyond the caller's
//! declared skin, no allocation in the hot path.

#![forbid(unsafe_code)]

use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_query::ConvexOrientedBox;

/// Failure to admit a scene body or a query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SceneError {
    /// A supplied coordinate, extent, or radius was not finite.
    NonFinite {
        /// Which field failed admission.
        field: &'static str,
    },
    /// A box half-extent was not strictly positive, or a skin was negative.
    NonPositive {
        /// Which field failed admission.
        field: &'static str,
    },
    /// A half-space normal was not usable as a direction.
    DegenerateNormal,
}

impl core::fmt::Display for SceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite { field } => write!(f, "scene field `{field}` is not finite"),
            Self::NonPositive { field } => {
                write!(f, "scene field `{field}` must be strictly positive")
            }
            Self::DegenerateNormal => write!(f, "half-space normal has no direction"),
        }
    }
}

impl core::error::Error for SceneError {}

/// What a body means for the things that touch it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BodyRole {
    /// Solid: no collider may overlap it at all. Walls, furniture, a cabinet.
    KeepOut,
    /// A surface things rest on. Contact is expected; sinking past the
    /// caller's skin depth is a violation. Floors, table tops, the ground.
    Support,
}

/// One piece of static environment geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SceneBody {
    /// A box, yawed about the world +Z axis.
    Box {
        /// Box centre in world coordinates \[m\].
        center_m: [f64; 3],
        /// Half extents along the box's own axes before yaw \[m\].
        half_extents_m: [f64; 3],
        /// Rotation about world +Z \[rad\].
        yaw_rad: f64,
    },
    /// The solid half-space on the negative side of a plane.
    ///
    /// `normal_m` points OUT of the solid, into free space. A ground plane at
    /// `z = 0` is `normal_m = [0, 0, 1]`, `offset_m = 0`: everything below
    /// `z = 0` is solid.
    HalfSpace {
        /// Outward unit normal (normalised on admission).
        normal_m: [f64; 3],
        /// Signed plane offset along the normal \[m\].
        offset_m: f64,
    },
}

impl SceneBody {
    /// A ground plane at height `height_m`, solid below it.
    #[must_use]
    pub const fn ground_plane(height_m: f64) -> Self {
        Self::HalfSpace {
            normal_m: [0.0, 0.0, 1.0],
            offset_m: height_m,
        }
    }

    fn validate(&self) -> Result<Self, SceneError> {
        match *self {
            Self::Box {
                center_m,
                half_extents_m,
                yaw_rad,
            } => {
                for value in center_m {
                    if !value.is_finite() {
                        return Err(SceneError::NonFinite { field: "center_m" });
                    }
                }
                for value in half_extents_m {
                    if !value.is_finite() {
                        return Err(SceneError::NonFinite {
                            field: "half_extents_m",
                        });
                    }
                    if value <= 0.0 {
                        return Err(SceneError::NonPositive {
                            field: "half_extents_m",
                        });
                    }
                }
                if !yaw_rad.is_finite() {
                    return Err(SceneError::NonFinite { field: "yaw_rad" });
                }
                Ok(*self)
            }
            Self::HalfSpace { normal_m, offset_m } => {
                for value in normal_m {
                    if !value.is_finite() {
                        return Err(SceneError::NonFinite { field: "normal_m" });
                    }
                }
                if !offset_m.is_finite() {
                    return Err(SceneError::NonFinite { field: "offset_m" });
                }
                let length = (normal_m[0] * normal_m[0]
                    + normal_m[1] * normal_m[1]
                    + normal_m[2] * normal_m[2])
                    .sqrt();
                if !(length > 1.0e-12) {
                    return Err(SceneError::DegenerateNormal);
                }
                Ok(Self::HalfSpace {
                    normal_m: [
                        normal_m[0] / length,
                        normal_m[1] / length,
                        normal_m[2] / length,
                    ],
                    offset_m,
                })
            }
        }
    }

    /// Depth by which a sphere of `radius_m` at `center_m` overlaps this body.
    ///
    /// Zero when the sphere is clear. Positive is real overlap, measured along
    /// the shortest exit. Exact closed form for both shapes: no iteration.
    #[must_use]
    pub fn sphere_overlap_depth(&self, center_m: &[f64; 3], radius_m: f64) -> f64 {
        match *self {
            Self::Box {
                center_m: box_center,
                half_extents_m: half,
                yaw_rad,
            } => {
                let (sin_yaw, cos_yaw) = yaw_rad.sin_cos();
                sphere_box_overlap_depth(center_m, radius_m, &box_center, &half, cos_yaw, sin_yaw)
            }
            Self::HalfSpace { normal_m, offset_m } => {
                let signed = normal_m[0] * center_m[0]
                    + normal_m[1] * center_m[1]
                    + normal_m[2] * center_m[2]
                    - offset_m;
                (radius_m - signed).max(0.0)
            }
        }
    }

    /// Radius of a sphere centred on this body's own centre that contains it,
    /// or `None` for an unbounded body.
    #[must_use]
    pub fn bounding_radius_m(&self) -> Option<f64> {
        match *self {
            Self::Box { half_extents_m, .. } => Some(
                (half_extents_m[0] * half_extents_m[0]
                    + half_extents_m[1] * half_extents_m[1]
                    + half_extents_m[2] * half_extents_m[2])
                    .sqrt(),
            ),
            Self::HalfSpace { .. } => None,
        }
    }

    /// Centre of a bounded body, or `None` for an unbounded one.
    #[must_use]
    pub const fn center_m(&self) -> Option<[f64; 3]> {
        match *self {
            Self::Box { center_m, .. } => Some(center_m),
            Self::HalfSpace { .. } => None,
        }
    }

    /// Convex support map for a bounded body, so an arbitrary convex collider
    /// can be tested against it with `fs-query`'s certified separation.
    ///
    /// Half-spaces are unbounded and have no support map; use
    /// [`Self::sphere_overlap_depth`] or a caller-side plane test for those.
    pub fn convex_support_map(&self) -> Option<ConvexOrientedBox> {
        match *self {
            Self::Box {
                center_m,
                half_extents_m,
                yaw_rad,
            } => {
                let (sin_yaw, cos_yaw) = yaw_rad.sin_cos();
                ConvexOrientedBox::new(
                    Point3::new(center_m[0], center_m[1], center_m[2]),
                    [
                        GeomVec3::new(cos_yaw, sin_yaw, 0.0),
                        GeomVec3::new(-sin_yaw, cos_yaw, 0.0),
                        GeomVec3::new(0.0, 0.0, 1.0),
                    ],
                    GeomVec3::new(half_extents_m[0], half_extents_m[1], half_extents_m[2]),
                )
                .ok()
            }
            Self::HalfSpace { .. } => None,
        }
    }
}

/// Sphere-versus-yawed-box overlap depth with the box's yaw trig supplied by
/// the caller, so a loop over many spheres against one box computes it once.
#[inline]
#[must_use]
pub fn sphere_box_overlap_depth(
    center_m: &[f64; 3],
    radius_m: f64,
    box_center_m: &[f64; 3],
    half_extents_m: &[f64; 3],
    cos_yaw: f64,
    sin_yaw: f64,
) -> f64 {
    let dx = center_m[0] - box_center_m[0];
    let dy = center_m[1] - box_center_m[1];
    let dz = center_m[2] - box_center_m[2];
    // World -> box frame: rotate by -yaw about Z.
    let lx = cos_yaw * dx + sin_yaw * dy;
    let ly = -sin_yaw * dx + cos_yaw * dy;
    let lz = dz;
    let qx = lx.clamp(-half_extents_m[0], half_extents_m[0]);
    let qy = ly.clamp(-half_extents_m[1], half_extents_m[1]);
    let qz = lz.clamp(-half_extents_m[2], half_extents_m[2]);
    let ddx = lx - qx;
    let ddy = ly - qy;
    let ddz = lz - qz;
    let outside_sq = ddx * ddx + ddy * ddy + ddz * ddz;
    if outside_sq > 0.0 {
        // Centre outside the solid: overlap iff the radius reaches the surface.
        (radius_m - outside_sq.sqrt()).max(0.0)
    } else {
        // Centre inside the solid: depth is the radius plus the distance to
        // the nearest face.
        let face = (half_extents_m[0] - lx.abs())
            .min((half_extents_m[1] - ly.abs()).min(half_extents_m[2] - lz.abs()));
        radius_m + face
    }
}

/// One admitted body plus what it means and how far things may sink into it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SceneEntry {
    /// The geometry.
    pub body: SceneBody,
    /// Keep-out or support.
    pub role: BodyRole,
    /// Overlap tolerated before it counts as a violation \[m\]. Keep-out
    /// bodies normally use a small positive skin so grazing contact does not
    /// fire; support bodies use the contact model's expected indentation.
    pub skin_m: f64,
}

/// A violation of the environment: what was hit, and by how much.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScenePenetration {
    /// Index of the offending body in the scene's declaration order.
    pub body_index: usize,
    /// Overlap depth \[m\], always greater than that body's skin.
    pub depth_m: f64,
    /// Depth in excess of the body's skin \[m\].
    pub excess_m: f64,
    /// What kind of body it was.
    pub role: BodyRole,
    /// Index of the offending collider in the caller's slice.
    pub collider_index: usize,
}

/// A scene body with its per-body constants hoisted out of the query loop.
///
/// A scan over many colliders against many bodies otherwise recomputes each
/// box's yaw trigonometry and bounding radius once per PAIR. Preparing the
/// scene once per rollout turns that into once per body, and adds a
/// bounding-sphere reject that skips the exact test for distant pairs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedBody {
    /// The admitted geometry.
    pub body: SceneBody,
    /// Keep-out or support.
    pub role: BodyRole,
    /// Overlap tolerated before it counts as a violation \[m\].
    pub skin_m: f64,
    cos_yaw: f64,
    sin_yaw: f64,
    center_m: [f64; 3],
    bounding_radius_m: f64,
    bounded: bool,
}

impl PreparedBody {
    fn new(entry: &SceneEntry) -> Self {
        let (cos_yaw, sin_yaw) = match entry.body {
            SceneBody::Box { yaw_rad, .. } => {
                let (s, c) = yaw_rad.sin_cos();
                (c, s)
            }
            SceneBody::HalfSpace { .. } => (1.0, 0.0),
        };
        Self {
            body: entry.body,
            role: entry.role,
            skin_m: entry.skin_m,
            cos_yaw,
            sin_yaw,
            center_m: entry.body.center_m().unwrap_or([0.0; 3]),
            bounding_radius_m: entry.body.bounding_radius_m().unwrap_or(0.0),
            bounded: entry.body.center_m().is_some(),
        }
    }

    /// Cheap conservative reject: `false` guarantees the sphere cannot breach
    /// this body's skin, so the exact test can be skipped. Never rejects a
    /// pair that could violate.
    #[inline]
    #[must_use]
    pub fn may_breach(&self, center_m: &[f64; 3], radius_m: f64) -> bool {
        if !self.bounded {
            return true;
        }
        let dx = center_m[0] - self.center_m[0];
        let dy = center_m[1] - self.center_m[1];
        let dz = center_m[2] - self.center_m[2];
        let reach = self.bounding_radius_m + radius_m;
        dx * dx + dy * dy + dz * dz <= reach * reach
    }

    /// Overlap depth, with the body's trigonometry already resolved.
    #[inline]
    #[must_use]
    pub fn sphere_overlap_depth(&self, center_m: &[f64; 3], radius_m: f64) -> f64 {
        match self.body {
            SceneBody::Box {
                center_m: box_center,
                half_extents_m: half,
                ..
            } => sphere_box_overlap_depth(
                center_m,
                radius_m,
                &box_center,
                &half,
                self.cos_yaw,
                self.sin_yaw,
            ),
            SceneBody::HalfSpace { .. } => self.body.sphere_overlap_depth(center_m, radius_m),
        }
    }
}

/// A validated static environment.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct StaticScene {
    entries: Vec<SceneEntry>,
}

impl StaticScene {
    /// An empty environment.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Admit one body. Returns its index.
    pub fn push(
        &mut self,
        body: SceneBody,
        role: BodyRole,
        skin_m: f64,
    ) -> Result<usize, SceneError> {
        if !skin_m.is_finite() {
            return Err(SceneError::NonFinite { field: "skin_m" });
        }
        if skin_m < 0.0 {
            return Err(SceneError::NonPositive { field: "skin_m" });
        }
        let body = body.validate()?;
        self.entries.push(SceneEntry { body, role, skin_m });
        Ok(self.entries.len() - 1)
    }

    /// Admitted bodies, in declaration order.
    #[must_use]
    pub fn entries(&self) -> &[SceneEntry] {
        &self.entries
    }

    /// How many bodies the environment holds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the environment is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Deepest violation over every (collider, body) pair, or `None` when
    /// every collider is within every body's declared skin.
    ///
    /// `colliders` are (centre, radius) spheres. This is the query a robot
    /// rollout runs per step: pass the body's collider set and refuse the step
    /// if it returns `Some`.
    #[must_use]
    pub fn deepest_sphere_penetration(
        &self,
        colliders: &[([f64; 3], f64)],
    ) -> Option<ScenePenetration> {
        let mut worst: Option<ScenePenetration> = None;
        for (body_index, entry) in self.entries.iter().enumerate() {
            // Bounded bodies get a cheap sphere reject before the exact test.
            let bound = entry.body.center_m().zip(entry.body.bounding_radius_m());
            for (collider_index, (center, radius)) in colliders.iter().enumerate() {
                if let Some((body_center, body_radius)) = bound {
                    let dx = center[0] - body_center[0];
                    let dy = center[1] - body_center[1];
                    let dz = center[2] - body_center[2];
                    let reach = body_radius + radius + entry.skin_m;
                    if dx * dx + dy * dy + dz * dz > reach * reach {
                        continue;
                    }
                }
                let depth = entry.body.sphere_overlap_depth(center, *radius);
                if depth <= entry.skin_m {
                    continue;
                }
                let excess = depth - entry.skin_m;
                if worst.is_none_or(|current| excess > current.excess_m) {
                    worst = Some(ScenePenetration {
                        body_index,
                        depth_m: depth,
                        excess_m: excess,
                        role: entry.role,
                        collider_index,
                    });
                }
            }
        }
        worst
    }

    /// Largest overlap depth over every pair, ignoring skins. Reported by
    /// rollouts that want to publish how close they came even when nothing
    /// was violated.
    #[must_use]
    pub fn maximum_sphere_overlap(&self, colliders: &[([f64; 3], f64)]) -> f64 {
        let mut deepest = 0.0_f64;
        for entry in &self.entries {
            for (center, radius) in colliders {
                let depth = entry.body.sphere_overlap_depth(center, *radius);
                if depth > deepest {
                    deepest = depth;
                }
            }
        }
        deepest
    }

    /// The scene with per-body constants hoisted, for a hot scan loop.
    #[must_use]
    pub fn prepare(&self) -> Vec<PreparedBody> {
        self.entries.iter().map(PreparedBody::new).collect()
    }

    /// Convex support maps for the bounded bodies, paired with their index,
    /// for callers that test non-spherical colliders through `fs-query`.
    #[must_use]
    pub fn convex_bodies(&self) -> Vec<(usize, ConvexOrientedBox)> {
        self.entries
            .iter()
            .enumerate()
            .filter_map(|(index, entry)| {
                entry.body.convex_support_map().map(|shape| (index, shape))
            })
            .collect()
    }
}

/// A convex collider's support map, for symmetry with [`StaticScene`] users.
pub use fs_query::ConvexSupportMap as SceneSupportMap;

#[cfg(test)]
mod tests {
    use super::*;

    fn boxed(center: [f64; 3], half: [f64; 3], yaw: f64) -> SceneBody {
        SceneBody::Box {
            center_m: center,
            half_extents_m: half,
            yaw_rad: yaw,
        }
    }

    #[test]
    fn a_sphere_clear_of_a_box_reports_no_overlap() {
        let body = boxed([0.0, 0.0, 0.0], [0.5, 0.5, 0.5], 0.0);
        assert_eq!(body.sphere_overlap_depth(&[2.0, 0.0, 0.0], 0.25), 0.0);
        // Touching exactly: still zero, not a violation.
        assert!(body.sphere_overlap_depth(&[0.75, 0.0, 0.0], 0.25).abs() < 1e-12);
    }

    #[test]
    fn overlap_depth_is_the_shortest_exit_outside_and_inside() {
        let body = boxed([0.0, 0.0, 0.0], [0.5, 0.5, 0.5], 0.0);
        // 0.1 m of the sphere is past the +x face.
        assert!((body.sphere_overlap_depth(&[0.65, 0.0, 0.0], 0.25) - 0.1).abs() < 1e-12);
        // Centre dead inside: radius plus the distance to the nearest face.
        assert!((body.sphere_overlap_depth(&[0.0, 0.0, 0.0], 0.25) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn yaw_rotates_the_box_and_the_query_follows_it() {
        let long = boxed([0.0, 0.0, 0.0], [1.0, 0.1, 0.5], 0.0);
        let turned = boxed(
            [0.0, 0.0, 0.0],
            [1.0, 0.1, 0.5],
            core::f64::consts::FRAC_PI_2,
        );
        // A point 0.5 m along +y clears the unrotated slab and hits the turned one.
        assert_eq!(long.sphere_overlap_depth(&[0.0, 0.5, 0.0], 0.05), 0.0);
        assert!(turned.sphere_overlap_depth(&[0.0, 0.5, 0.0], 0.05) > 0.0);
    }

    #[test]
    fn a_ground_plane_is_solid_below_its_height() {
        let ground = SceneBody::ground_plane(0.0).validate().expect("admits");
        // A sphere resting exactly on the surface just touches.
        assert!(ground.sphere_overlap_depth(&[0.0, 0.0, 0.25], 0.25).abs() < 1e-12);
        // Sunk 0.1 m.
        assert!((ground.sphere_overlap_depth(&[0.0, 0.0, 0.15], 0.25) - 0.1).abs() < 1e-12);
        // Well above it.
        assert_eq!(ground.sphere_overlap_depth(&[0.0, 0.0, 1.0], 0.25), 0.0);
    }

    #[test]
    fn a_support_skin_tolerates_contact_but_not_sinking() {
        let mut scene = StaticScene::new();
        scene
            .push(SceneBody::ground_plane(0.0), BodyRole::Support, 0.01)
            .expect("ground admits");
        // 3 mm of compliant indentation: normal contact.
        assert!(
            scene
                .deepest_sphere_penetration(&[([0.0, 0.0, 0.247], 0.25)])
                .is_none()
        );
        // 27 cm below the floor: the failure this crate exists to catch.
        let sunk = scene
            .deepest_sphere_penetration(&[([0.0, 0.0, -0.02], 0.25)])
            .expect("a body a quarter metre inside the floor is a violation");
        assert_eq!(sunk.role, BodyRole::Support);
        assert!((sunk.depth_m - 0.27).abs() < 1e-12);
        assert!((sunk.excess_m - 0.26).abs() < 1e-12);
        assert_eq!(sunk.collider_index, 0);
    }

    #[test]
    fn the_deepest_violation_wins_and_names_its_collider() {
        let mut scene = StaticScene::new();
        scene
            .push(
                boxed([0.0, 0.0, 0.0], [0.5, 0.5, 0.5], 0.0),
                BodyRole::KeepOut,
                0.0,
            )
            .expect("box admits");
        scene
            .push(
                boxed([3.0, 0.0, 0.0], [0.5, 0.5, 0.5], 0.0),
                BodyRole::KeepOut,
                0.0,
            )
            .expect("box admits");
        let hit = scene
            .deepest_sphere_penetration(&[
                ([0.62, 0.0, 0.0], 0.1),  // 0.02 into body 0
                ([3.40, 0.0, 0.0], 0.25), // 0.35 into body 1
            ])
            .expect("both overlap");
        assert_eq!(hit.body_index, 1);
        assert_eq!(hit.collider_index, 1);
        assert!((hit.depth_m - 0.35).abs() < 1e-12);
    }

    #[test]
    fn the_broad_phase_reject_never_changes_an_answer() {
        let mut scene = StaticScene::new();
        for index in 0..40 {
            scene
                .push(
                    boxed(
                        [index as f64 * 0.5, 0.0, 0.0],
                        [0.2, 0.2, 0.2],
                        0.1 * index as f64,
                    ),
                    BodyRole::KeepOut,
                    0.0,
                )
                .expect("admits");
        }
        for step in 0..60 {
            let center = [step as f64 * 0.31, 0.05, 0.02];
            let collider = [(center, 0.15)];
            let fast = scene.deepest_sphere_penetration(&collider);
            // Exhaustive reference: same test with no reject.
            let mut reference: Option<ScenePenetration> = None;
            for (body_index, entry) in scene.entries().iter().enumerate() {
                let depth = entry.body.sphere_overlap_depth(&center, 0.15);
                if depth <= entry.skin_m {
                    continue;
                }
                let excess = depth - entry.skin_m;
                if reference.is_none_or(|c| excess > c.excess_m) {
                    reference = Some(ScenePenetration {
                        body_index,
                        depth_m: depth,
                        excess_m: excess,
                        role: entry.role,
                        collider_index: 0,
                    });
                }
            }
            assert_eq!(fast, reference, "step {step}");
        }
    }

    #[test]
    fn malformed_bodies_are_refused_by_name() {
        let mut scene = StaticScene::new();
        assert_eq!(
            scene.push(
                boxed([0.0; 3], [0.0, 0.1, 0.1], 0.0),
                BodyRole::KeepOut,
                0.0
            ),
            Err(SceneError::NonPositive {
                field: "half_extents_m"
            })
        );
        assert_eq!(
            scene.push(
                boxed([f64::NAN, 0.0, 0.0], [0.1; 3], 0.0),
                BodyRole::KeepOut,
                0.0
            ),
            Err(SceneError::NonFinite { field: "center_m" })
        );
        assert_eq!(
            scene.push(
                SceneBody::HalfSpace {
                    normal_m: [0.0; 3],
                    offset_m: 0.0
                },
                BodyRole::Support,
                0.0
            ),
            Err(SceneError::DegenerateNormal)
        );
        assert_eq!(
            scene.push(SceneBody::ground_plane(0.0), BodyRole::Support, -0.1),
            Err(SceneError::NonPositive { field: "skin_m" })
        );
        assert_eq!(
            scene.push(SceneBody::ground_plane(0.0), BodyRole::Support, f64::NAN),
            Err(SceneError::NonFinite { field: "skin_m" })
        );
        assert!(scene.is_empty(), "nothing malformed was admitted");
    }

    #[test]
    fn a_declared_normal_is_normalised_on_admission() {
        let mut scene = StaticScene::new();
        scene
            .push(
                SceneBody::HalfSpace {
                    normal_m: [0.0, 0.0, 5.0],
                    offset_m: 2.0,
                },
                BodyRole::Support,
                0.0,
            )
            .expect("admits");
        // Sunk 0.1 m below a surface at z = 2.
        let hit = scene
            .deepest_sphere_penetration(&[([0.0, 0.0, 2.15], 0.25)])
            .expect("violates");
        assert!((hit.depth_m - 0.1).abs() < 1e-12);
    }

    #[test]
    fn bounded_bodies_expose_a_convex_support_map_and_planes_do_not() {
        let scene_box = boxed([1.0, 2.0, 3.0], [0.4, 0.5, 0.6], 0.3);
        assert!(scene_box.convex_support_map().is_some());
        assert!(SceneBody::ground_plane(0.0).convex_support_map().is_none());
        assert!(SceneBody::ground_plane(0.0).bounding_radius_m().is_none());
        assert!(
            (scene_box.bounding_radius_m().expect("bounded")
                - (0.4_f64 * 0.4 + 0.5 * 0.5 + 0.6 * 0.6).sqrt())
            .abs()
                < 1e-12
        );
    }
}

#[cfg(test)]
mod prepared_tests {
    use super::*;

    fn scene() -> StaticScene {
        let mut scene = StaticScene::new();
        for index in 0..30 {
            scene
                .push(
                    SceneBody::Box {
                        center_m: [index as f64 * 0.37 - 3.0, (index % 5) as f64 * 0.4, 0.5],
                        half_extents_m: [0.22, 0.18, 0.3],
                        yaw_rad: 0.11 * index as f64,
                    },
                    BodyRole::KeepOut,
                    0.01,
                )
                .expect("admits");
        }
        scene
            .push(SceneBody::ground_plane(0.0), BodyRole::Support, 0.05)
            .expect("admits");
        scene
    }

    /// The prepared form is a pure optimisation: identical depths, and the
    /// reject never hides a breach.
    #[test]
    fn prepared_bodies_agree_with_the_direct_query_everywhere() {
        let scene = scene();
        let prepared = scene.prepare();
        assert_eq!(prepared.len(), scene.len());
        for step in 0..120 {
            let center = [
                (step as f64) * 0.09 - 4.0,
                ((step % 11) as f64) * 0.21 - 1.0,
                ((step % 7) as f64) * 0.15 - 0.2,
            ];
            for radius in [0.05, 0.17, 0.4] {
                for (index, body) in prepared.iter().enumerate() {
                    let direct = scene.entries()[index]
                        .body
                        .sphere_overlap_depth(&center, radius);
                    let fast = body.sphere_overlap_depth(&center, radius);
                    assert!(
                        (direct - fast).abs() < 1e-12,
                        "body {index} disagreed: {direct} vs {fast}"
                    );
                    if !body.may_breach(&center, radius) {
                        // The reject is bounding-sphere based, so a rejected
                        // pair cannot touch the body at all, not merely stay
                        // inside its skin.
                        assert!(
                            direct == 0.0,
                            "the reject discarded a pair with overlap {direct}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn an_unbounded_body_is_never_rejected() {
        let scene = scene();
        let prepared = scene.prepare();
        let ground = prepared.last().expect("ground is last");
        assert!(matches!(ground.body, SceneBody::HalfSpace { .. }));
        assert!(ground.may_breach(&[500.0, -500.0, 900.0], 0.01));
    }
}
