//! Bounded local shooting from a point source to a native receiving endpoint.
//!
//! Each trial uses the native scene's one-pass differentiated path. The two
//! analytic arrival columns give a Newton correction in launch solid angle;
//! bounded steps and backtracking retain the original ordered optical branch.
//! No finite-difference retraces or prescribed interface positions are used.
//!
//! The initial direction must already reach the selected finite receiver on a
//! regular branch. This is a local connection solver, not a global search for
//! all images. Blockers, missing detector support, changed features or media,
//! caustics and numerical contact resolution can prevent convergence. Success
//! checks the full three-dimensional endpoint error without snapping the ray.
//! Backtracking handles completed trials that change branch or miss the target;
//! native field, interface and geometry errors terminate the solve immediately.

use super::super::CurvedPolarizedScene;
use super::{CurvedBundlePathResult, CurvedBundleSceneError, CurvedBundleSettings};
use crate::charts::Ray;
use crate::dielectric::DielectricEvent;
use crate::gradient_index::bundle::{BundleError, PlanarBundleProjection, ProjectionOptions};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{MAX_POLARIZED_PATH_EVENTS, PathTermination};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Deterministic finite-seed discovery of distinct regular optical images.
pub mod images;

/// One fixed source and polarization law while its launch direction is solved.
/// The transverse axis is projected into each trial ray's transverse plane;
/// the normalized Stokes coordinates are defined in that projected frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleConnectionLaunch {
    pub ray: Ray,
    pub transverse_axis: Vec3,
    pub wavelength_nm: f64,
    pub normalized_source_stokes: Stokes,
    pub initial_medium: usize,
}

/// Deterministic optical branch selection, replayed from the beginning for
/// every Newton/backtracking trial. The original source's physical Fresnel
/// power is retained; these choices are not Monte Carlo proposal probabilities.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BundleBranchPolicy<'a> {
    /// Transmit whenever permitted and reflect at total internal reflection.
    #[default]
    TransmitOrReflectTir,
    /// Exactly one decision per encountered dielectric boundary, including
    /// coherent films and authored multilayers. Conductors always reflect and
    /// do not consume decisions. The whole sequence must be consumed when the
    /// receiver is reached. Empty sequences admit paths with no dielectric hit.
    Sequence(&'a [DielectricEvent]),
}

/// Target point on a named native planar receiver or in an original medium,
/// with a section normal toward light.
/// Position tolerance is a numerical endpoint allowance, not an aperture or an
/// error bound on irradiance. The scene supplies a surface receiver's finite
/// aperture; a volume target instead validates the point's original medium.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleConnectionTarget {
    object_id: u64,
    position_m: Point3,
    normal: Vec3,
    basis: [Vec3; 2],
    position_tolerance_m: f64,
    kind: ReceiverKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ReceiverKind {
    Detector,
    Diffuse,
    Volume { medium: usize },
}

impl BundleConnectionTarget {
    /// An original `Interaction::Detector`; ordinary opaque or diffuse
    /// surfaces cannot replace this receiving endpoint.
    pub fn try_new(
        object_id: u64,
        position_m: Point3,
        normal: Vec3,
        position_tolerance_m: f64,
    ) -> Result<Self, BundleConnectionError> {
        if ![position_m.x, position_m.y, position_m.z]
            .iter()
            .all(|x| x.is_finite())
            || !position_tolerance_m.is_finite()
            || position_tolerance_m <= 0.0
        {
            return Err(BundleConnectionError::InvalidTarget);
        }
        let normal = unit(normal).ok_or(BundleConnectionError::InvalidTarget)?;
        let axis = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
            Vec3::new(1.0, 0.0, 0.0)
        } else if normal.y.abs() <= normal.z.abs() {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let first = unit(cross(normal, axis)).ok_or(BundleConnectionError::InvalidTarget)?;
        Ok(Self {
            object_id,
            position_m,
            normal,
            basis: [first, cross(normal, first)],
            position_tolerance_m,
            kind: ReceiverKind::Detector,
        })
    }

    /// Receive incident power at one original planar depolarizing diffuse
    /// surface, before applying its reflectance or sampling a reflected ray.
    /// The native walk validates the arriving medium and finite geometry. Only
    /// this selected object is a receiving endpoint; other diffuse encounters
    /// remain unsupported. Detector targets keep their original semantics.
    pub fn try_diffuse_receiver(
        object_id: u64,
        position_m: Point3,
        normal: Vec3,
        position_tolerance_m: f64,
    ) -> Result<Self, BundleConnectionError> {
        let mut result = Self::try_new(object_id, position_m, normal, position_tolerance_m)?;
        result.kind = ReceiverKind::Diffuse;
        Ok(result)
    }

    /// Receive a source connection at a point in an original scene medium.
    /// `section_normal` defines only the numerical shooting section, oriented
    /// toward the source. The native geometry remains authoritative before it;
    /// no detector, surface interaction, or substituted geometry is created.
    /// The section is active only while the walk occupies `medium`.
    ///
    /// For molecular radiometry use the final ray's transverse area, i.e. a
    /// projection normal parallel to its direction, not this numerical section.
    /// Position tolerance bounds the resolved endpoint residual, not irradiance
    /// error. A point inside an unresolved physical contact collar is refused.
    pub fn try_volume_receiver(
        medium: usize,
        position_m: Point3,
        section_normal: Vec3,
        position_tolerance_m: f64,
    ) -> Result<Self, BundleConnectionError> {
        let mut result = Self::try_new(0, position_m, section_normal, position_tolerance_m)?;
        result.kind = ReceiverKind::Volume { medium };
        Ok(result)
    }

    /// Native surface identity for detector/diffuse targets. Volume targets
    /// have no object and return zero; use `volume_medium_index` for them.
    pub const fn object_id(self) -> u64 {
        self.object_id
    }
    pub const fn position_m(self) -> Point3 {
        self.position_m
    }
    pub const fn normal(self) -> Vec3 {
        self.normal
    }
    pub const fn position_tolerance_m(self) -> f64 {
        self.position_tolerance_m
    }

    /// This target receives incident power at a selected diffuse surface.
    pub const fn is_diffuse_receiver(self) -> bool {
        matches!(self.kind, ReceiverKind::Diffuse)
    }

    pub const fn is_volume_receiver(self) -> bool {
        matches!(self.kind, ReceiverKind::Volume { .. })
    }

    pub const fn volume_medium_index(self) -> Option<usize> {
        match self.kind {
            ReceiverKind::Volume { medium } => Some(medium),
            _ => None,
        }
    }

    pub(crate) fn native_receiver(self) -> Option<super::BundleReceiver> {
        match self.kind {
            ReceiverKind::Detector => None,
            ReceiverKind::Diffuse => Some(super::BundleReceiver::Diffuse(self.object_id)),
            ReceiverKind::Volume { medium } => Some(super::BundleReceiver::Volume {
                medium,
                position_m: self.position_m,
                normal: self.normal,
            }),
        }
    }
}

/// Bounds on local shooting in addition to each native path's work limits.
/// At most `max_trace_evaluations` paths are traced, each with `scene`'s limits.
/// The angular step bounds the actual direction change in radians because
/// normalizing a tangent increment of length `h` rotates by `atan(h) <= h`.
#[derive(Clone, Copy, Debug)]
pub struct BundleConnectionOptions {
    pub scene: CurvedBundleSettings,
    pub projection: ProjectionOptions,
    pub max_iterations: u32,
    pub max_trace_evaluations: u32,
    pub max_backtracks: u32,
    pub max_angular_step_rad: f64,
}

impl Default for BundleConnectionOptions {
    fn default() -> Self {
        Self {
            scene: CurvedBundleSettings::default(),
            projection: ProjectionOptions::default(),
            max_iterations: 16,
            max_trace_evaluations: 96,
            max_backtracks: 8,
            max_angular_step_rad: 0.25,
        }
    }
}

impl BundleConnectionOptions {
    pub(crate) fn validate(self) -> Result<(), BundleConnectionError> {
        self.scene.validate()?;
        if self.max_iterations == 0
            || self.max_iterations > 256
            || self.max_trace_evaluations == 0
            || self.max_trace_evaluations > 4096
            || self.max_backtracks > 32
            || !self.max_angular_step_rad.is_finite()
            || self.max_angular_step_rad <= 0.0
            || self.max_angular_step_rad > 1.0
        {
            return Err(BundleConnectionError::InvalidOptions);
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BundleConnectionError {
    InvalidTarget,
    InvalidOptions,
    InvalidBranchPlan,
    BranchPlanNotConsumed {
        planned: usize,
        consumed: usize,
    },
    ReceiverNotReached(PathTermination),
    WrongDetector,
    NonPlanarReceiver,
    ReceiverPlaneMismatch,
    BackFacingReceiver,
    BranchChanged,
    /// The tangent residual is resolved, but the unsnapped native endpoint is
    /// still farther from the target than its requested numerical tolerance.
    ContactResolution {
        residual_m: f64,
    },
    IterationLimit {
        residual_m: f64,
    },
    TraceBudgetExceeded,
    NoDescent {
        residual_m: f64,
    },
    Cancelled,
    NumericalRange,
    Scene(CurvedBundleSceneError),
    Projection(BundleError),
}

impl From<CurvedBundleSceneError> for BundleConnectionError {
    fn from(error: CurvedBundleSceneError) -> Self {
        Self::Scene(error)
    }
}
impl From<BundleError> for BundleConnectionError {
    fn from(error: BundleError) -> Self {
        Self::Projection(error)
    }
}
impl core::fmt::Display for BundleConnectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native optical connection: {self:?}")
    }
}
impl core::error::Error for BundleConnectionError {}

/// A converged visible branch and actual work from all its shooting trials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BundleConnection {
    path: CurvedBundlePathResult,
    residual_m: f64,
    iterations: u32,
    work: Work,
}

impl BundleConnection {
    pub const fn path(&self) -> &CurvedBundlePathResult {
        &self.path
    }
    pub const fn residual_m(&self) -> f64 {
        self.residual_m
    }
    pub const fn iterations(&self) -> u32 {
        self.iterations
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.work.traces
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.work.fields
    }
    pub const fn query_work(&self) -> u64 {
        self.work.geometry
    }

    // Source-survival replay uses this SAME receiver admission and ordered
    // optical history. Its stochastic subdivisions do not replace the solved
    // connection's geometric Jacobian or physical power transfer.
    pub(crate) fn validate_replay(
        &self,
        replay: CurvedBundlePathResult,
        target: BundleConnectionTarget,
        projection: ProjectionOptions,
    ) -> Result<(), BundleConnectionError> {
        if replay.branch_signature() != self.path.branch_signature() {
            return Err(BundleConnectionError::BranchChanged);
        }
        let replay = arrival(replay, target, projection)?;
        if replay.position_error > target.position_tolerance_m {
            return Err(BundleConnectionError::ContactResolution {
                residual_m: replay.position_error,
            });
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Work {
    traces: u32,
    fields: u64,
    geometry: u64,
}

struct Arrival {
    path: CurvedBundlePathResult,
    projection: PlanarBundleProjection,
    residual: [f64; 2],
    tangent_error: f64,
    normal_error: f64,
    position_error: f64,
}

impl CurvedPolarizedScene<'_> {
    /// Solve the local branch that transmits at dielectric interfaces and
    /// reflects under total internal reflection (and at smooth conductors).
    /// Source origin, wavelength, media and authored polarization law are fixed.
    /// The first direction must reach the selected receiver; every accepted
    /// iterate retains its ordered native geometry/feature/branch signature.
    /// Failure returns no partially connected optical result.
    pub fn connect_bundle_transmitted(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
    ) -> Result<BundleConnection, BundleConnectionError> {
        self.connect_bundle_with_policy(
            cx,
            launch,
            target,
            options,
            BundleBranchPolicy::TransmitOrReflectTir,
        )
    }

    /// Connect one prescribed reflected/transmitted optical branch. A plan is
    /// replayed for every candidate launch, so Newton iteration cannot consume
    /// decisions belonging to a later path. Geometry, local media and Fresnel
    /// transfer remain native; no supplied leg length or contact is accepted.
    /// This shares the transmitted connector's local-basin, visibility,
    /// numerical-resolution and bounded-work requirements.
    /// Exactly dark branches retain zero power while geometric continuation
    /// validates the remaining interfaces and actual receiver; they are not
    /// confused with missing, blocked or numerically unresolved connections.
    pub fn connect_bundle_with_policy(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
    ) -> Result<BundleConnection, BundleConnectionError> {
        self.connect_bundle_with_policy_work(
            cx,
            launch,
            target,
            options,
            policy,
            &mut Work::default(),
        )
    }

    // The caller retains work even when a completed trial misses the receiver
    // or a bounded local solve fails. Image discovery uses this same solver;
    // it neither retraces the first ray nor estimates failed-attempt work.
    fn connect_bundle_with_policy_work(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        work: &mut Work,
    ) -> Result<BundleConnection, BundleConnectionError> {
        options.validate()?;
        if matches!(policy, BundleBranchPolicy::Sequence(events)
            if events.len() > MAX_POLARIZED_PATH_EVENTS as usize)
        {
            return Err(BundleConnectionError::InvalidBranchPlan);
        }
        let first = self.connection_trial(cx, launch, target, options, policy, work)?;
        let signature = first.branch_signature();
        let mut current = arrival(first, target, options.projection)?;
        for iteration in 0..=options.max_iterations {
            cx.checkpoint()
                .map_err(|_| BundleConnectionError::Cancelled)?;
            if current.position_error <= target.position_tolerance_m {
                return Ok(BundleConnection {
                    path: current.path,
                    residual_m: current.position_error,
                    iterations: iteration,
                    work: *work,
                });
            }
            if current.tangent_error <= target.position_tolerance_m
                && current.normal_error > target.position_tolerance_m
            {
                return Err(BundleConnectionError::ContactResolution {
                    residual_m: current.position_error,
                });
            }
            if iteration == options.max_iterations {
                return Err(BundleConnectionError::IterationLimit {
                    residual_m: current.position_error,
                });
            }
            let step = newton_step(&current, target, options.max_angular_step_rad)?;
            let initial = current.path.bundle().propagation();
            let direction = initial.initial_ray().dir;
            let axes = initial.launch_basis();
            let mut damping = 1.0;
            let mut accepted = None;
            for _ in 0..=options.max_backtracks {
                let direction = unit(add(
                    direction,
                    add(
                        axes[0].scale(step[0] * damping),
                        axes[1].scale(step[1] * damping),
                    ),
                ))
                .ok_or(BundleConnectionError::NumericalRange)?;
                let trial_launch = BundleConnectionLaunch {
                    ray: Ray {
                        origin: launch.ray.origin,
                        dir: direction,
                    },
                    ..launch
                };
                let path =
                    self.connection_trial(cx, trial_launch, target, options, policy, work)?;
                // A changed visibility/topology cannot use the current smooth
                // Jacobian. Smaller trials may remain in the admitted basin.
                if path.branch_signature() == signature {
                    match arrival(path, target, options.projection) {
                        Ok(candidate)
                            if candidate.position_error <= target.position_tolerance_m
                                || candidate.tangent_error < current.tangent_error =>
                        {
                            accepted = Some(candidate);
                            break;
                        }
                        Ok(_) => {}
                        Err(
                            BundleConnectionError::ReceiverNotReached(_)
                            | BundleConnectionError::WrongDetector
                            | BundleConnectionError::NonPlanarReceiver
                            | BundleConnectionError::ReceiverPlaneMismatch
                            | BundleConnectionError::BackFacingReceiver
                            | BundleConnectionError::Projection(
                                BundleError::SingularProjection | BundleError::GrazingProjection,
                            ),
                        ) => {}
                        Err(error) => return Err(error),
                    }
                }
                damping *= 0.5;
            }
            current = accepted.ok_or(BundleConnectionError::NoDescent {
                residual_m: current.position_error,
            })?;
        }
        Err(BundleConnectionError::IterationLimit {
            residual_m: current.position_error,
        })
    }

    fn connection_trial(
        &self,
        cx: &Cx<'_>,
        launch: BundleConnectionLaunch,
        target: BundleConnectionTarget,
        options: BundleConnectionOptions,
        policy: BundleBranchPolicy<'_>,
        work: &mut Work,
    ) -> Result<CurvedBundlePathResult, BundleConnectionError> {
        cx.checkpoint()
            .map_err(|_| BundleConnectionError::Cancelled)?;
        if work.traces >= options.max_trace_evaluations {
            return Err(BundleConnectionError::TraceBudgetExceeded);
        }
        work.traces += 1;
        let bundle = self.launch_bundle(
            launch.ray,
            launch.transverse_axis,
            launch.wavelength_nm,
            launch.normalized_source_stokes,
            launch.initial_medium,
        )?;
        let (path, branch_plan) = match policy {
            BundleBranchPolicy::TransmitOrReflectTir => (
                self.trace_bundle_with_checked_branches_to_receiver(
                    cx,
                    bundle,
                    launch.initial_medium,
                    options.scene,
                    true,
                    target.native_receiver(),
                    |_, context| {
                        Ok(if context.total_internal_reflection {
                            DielectricEvent::Reflection
                        } else {
                            DielectricEvent::Transmission
                        })
                    },
                )?,
                None,
            ),
            BundleBranchPolicy::Sequence(events) => {
                let mut consumed = 0;
                let path = self.trace_bundle_with_checked_branches_to_receiver(
                    cx,
                    bundle,
                    launch.initial_medium,
                    options.scene,
                    true,
                    target.native_receiver(),
                    |_, _| {
                        let event = events.get(consumed).copied().ok_or(
                            CurvedBundleSceneError::BranchPlanExhausted {
                                dielectric_event: consumed,
                            },
                        )?;
                        consumed += 1;
                        Ok(event)
                    },
                )?;
                (path, Some((events.len(), consumed)))
            }
        };
        work.fields = work
            .fields
            .checked_add(path.field_evaluations())
            .ok_or(BundleConnectionError::NumericalRange)?;
        work.geometry = work
            .geometry
            .checked_add(path.query_work())
            .ok_or(BundleConnectionError::NumericalRange)?;
        if let Some((planned, consumed)) = branch_plan {
            if (path.termination() == PathTermination::Detected
                || path.is_diffuse_receiver()
                || path.is_volume_receiver())
                && consumed != planned
            {
                return Err(BundleConnectionError::BranchPlanNotConsumed { planned, consumed });
            }
        }
        Ok(path)
    }
}

fn arrival(
    path: CurvedBundlePathResult,
    target: BundleConnectionTarget,
    options: ProjectionOptions,
) -> Result<Arrival, BundleConnectionError> {
    let reached = match target.kind {
        ReceiverKind::Detector => path.termination() == PathTermination::Detected,
        ReceiverKind::Diffuse => {
            path.termination() == PathTermination::Opaque && path.is_diffuse_receiver()
        }
        ReceiverKind::Volume { medium } => {
            path.is_volume_receiver() && path.medium_index() == medium
        }
    };
    if !reached {
        return Err(BundleConnectionError::ReceiverNotReached(
            path.termination(),
        ));
    }
    if !target.is_volume_receiver() {
        let contact = path
            .last_contact()
            .ok_or(BundleConnectionError::WrongDetector)?;
        if contact.object_id != target.object_id {
            return Err(BundleConnectionError::WrongDetector);
        }
        if contact.witness.departure_radius_m.is_none() {
            return Err(BundleConnectionError::NonPlanarReceiver);
        }
        // Native planar normals and the supplied receiver orientation must agree
        // to roundoff-scale angular resolution; either geometric orientation is valid.
        if norm(cross(contact.witness.normal, target.normal)) > 1e-8
            || contact
                .witness
                .surface_point
                .delta_from(target.position_m)
                .dot(target.normal)
                .abs()
                > target.position_tolerance_m
        {
            return Err(BundleConnectionError::ReceiverPlaneMismatch);
        }
    }
    let bundle = path.bundle().propagation();
    let endpoint = bundle.propagation().ray;
    if endpoint.dir.dot(target.normal) >= 0.0 {
        return Err(BundleConnectionError::BackFacingReceiver);
    }
    let delta = endpoint.origin.delta_from(target.position_m);
    let residual = target.basis.map(|axis| delta.dot(axis));
    let tangent_error = det::hypot(residual[0], residual[1]);
    let position_error = norm(delta);
    if !position_error.is_finite() || !tangent_error.is_finite() {
        return Err(BundleConnectionError::NumericalRange);
    }
    Ok(Arrival {
        path,
        projection: bundle.project_at_endpoint(target.normal, options)?,
        residual,
        tangent_error,
        normal_error: delta.dot(target.normal).abs(),
        position_error,
    })
}

fn newton_step(
    arrival: &Arrival,
    target: BundleConnectionTarget,
    max_step: f64,
) -> Result<[f64; 2], BundleConnectionError> {
    let columns = arrival.projection.angular_derivatives();
    let matrix = target
        .basis
        .map(|axis| columns.map(|c| axis.dot(c.position_m_per_rad)));
    let scale = matrix.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return Err(BundleError::SingularProjection.into());
    }
    let [[a, b], [c, d]] = matrix.map(|row| row.map(|x| x / scale));
    let determinant = a * d - b * c;
    if determinant == 0.0 || !determinant.is_finite() {
        return Err(BundleError::SingularProjection.into());
    }
    let [r, s] = arrival.residual.map(|x| x / scale);
    let mut step = [(b * s - d * r) / determinant, (c * r - a * s) / determinant];
    let length = det::hypot(step[0], step[1]);
    if !length.is_finite() {
        return Err(BundleConnectionError::NumericalRange);
    }
    if length > max_step {
        step = step.map(|x| (x / length) * max_step);
    }
    Ok(step)
}

fn norm(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
fn unit(v: Vec3) -> Option<Vec3> {
    let length = norm(v);
    if !length.is_finite() || length <= 0.0 {
        return None;
    }
    Some(Vec3::new(v.x / length, v.y / length, v.z / length))
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod diffuse_tests;

#[cfg(test)]
mod volume_tests;
