//! Automatic sphere/plane geometry for the native force/contact/drift step.
//!
//! Sphere geometry is centred on each original COM; its mass distribution may
//! have arbitrary physical principal inertia inside that sphere. Constant world
//! forces, body-frame torques and gravity use the existing native kick, followed
//! by one simultaneous normal/Coulomb solve and the original free drift.
//!
//! Complete straight COM drift paths are checked analytically, including pairs
//! omitted from the current contact solve. A rejected crossing shortens the
//! timestep; no position is projected. This is geometrically refined FIRST-order
//! stepping, not time-error control or interval-certified arbitrary-shape CCD.
//! Positive activation gaps are an explicitly reported spatial resolution.
//! Impact/contact loss stays unallocated, and force splitting work is not heat.

use super::spheres::ContactSphere;
use super::step::{
    RigidContactPath, RigidContactStepControls, RigidContactStepError, RigidContactStepReport,
    step_rigid_frictional_contacts,
};
use super::{
    RigidFrictionalContact, RigidFrictionalImpactError, RigidImpactCounterpart, RigidImpactPoint,
};
use crate::{DynamicsError, Gravity, MassProperties, RigidBodyState, Vec3, Wrench};
use fs_contact::broad_phase::{
    SweepConfig, SweepError, SweepItem, SweepPair, sweep_candidates, swept_sphere_aabb,
};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_contact::swept_sphere::{SweptSphereInput, first_sphere_contact};
use fs_exec::CancelGate;

/// Infinite plane translating with constant velocity and a fixed normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereScenePlane {
    /// Stable plane identity, unique among planes.
    pub id: u64,
    /// A point on the plane at absolute time zero.
    pub point_world: Vec3,
    /// Normal towards the permitted half-space; normalized at admission.
    pub normal_world: Vec3,
    /// Constant world velocity of the prescribed support.
    pub velocity_world: Vec3,
}

/// Stable source identities for resolving an actual detected contact law.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SphereSceneContactKey {
    /// Two sphere IDs in increasing order.
    Spheres { first: u64, second: u64 },
    /// A dynamic sphere and an independently identified prescribed plane.
    SpherePlane { sphere: u64, plane: u64 },
}

/// Source-resolved ideal point-contact coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereSceneLaw {
    pub friction_coefficient: f64,
    pub restitution: f64,
}

/// Physical gates and bounded geometry/time work.
#[derive(Clone, Copy, Debug)]
pub struct SphereSceneControls {
    /// Existing native force/contact/drift controls, also used for global gates.
    pub step: RigidContactStepControls,
    /// Complete broad-phase candidate/comparison caps, per query.
    pub broad_phase: SweepConfig,
    /// Maximum FIRST-order time interval; the caller controls time accuracy.
    pub maximum_step_s: f64,
    /// Maximum accepted intervals, at most 65,536.
    pub maximum_steps: usize,
    /// Maximum rejected intervals across the whole call, at most 65,536.
    pub maximum_rejections: usize,
    /// Whole-call cap on exact sphere/plane gap and complete-path queries.
    pub maximum_plane_queries: usize,
    /// Positive contact activation distance, at most native geometry tolerance.
    ///
    /// Contact may act before exact touch by this declared spatial resolution.
    pub activation_gap_m: f64,
}

/// Analytic admission of the entire proposed numerical drift.
#[derive(Clone, Copy, Debug)]
pub struct SphereScenePathReport {
    /// Minimum over queried candidate pairs and all plane paths.
    /// Broad-phase exclusions have strictly positive separation; if no exact
    /// query was needed this value is positive infinity.
    pub minimum_path_gap_m: f64,
    pub sphere_pairs_checked: usize,
    pub plane_pairs_checked: usize,
}

/// One accepted native step, including its original physics receipt.
pub struct SphereSceneStep {
    /// End of the accepted interval, relative to this call's start.
    pub elapsed_s: f64,
    /// Stable keys in the native contact-ID order.
    pub contact_keys: Vec<SphereSceneContactKey>,
    /// Largest actual geometric gap among activated contacts.
    pub maximum_contact_gap_m: f64,
    pub native: RigidContactStepReport<SphereScenePathReport>,
}
impl core::fmt::Debug for SphereSceneStep {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SphereSceneStep")
            .field("elapsed_s", &self.elapsed_s)
            .field("contact_keys", &self.contact_keys)
            .field("maximum_contact_gap_m", &self.maximum_contact_gap_m)
            .field("path", &self.native.path_admission)
            .finish_non_exhaustive()
    }
}

/// Whole-call trajectory and work/momentum accounting.
#[derive(Debug)]
pub struct SphereSceneReport {
    pub duration_s: f64,
    pub steps: Vec<SphereSceneStep>,
    pub rejected_steps: usize,
    pub sphere_pair_queries: usize,
    pub broad_phase_comparisons: usize,
    pub plane_queries: usize,
    /// Minimum among all exact complete-path queries; see the path receipt.
    pub minimum_path_gap_m: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub force_impulse_work_j: f64,
    pub force_displacement_work_j: f64,
    /// Signed kick-work minus displacement-work defect, for time refinement.
    pub force_splitting_work_j: f64,
    pub torque_impulse_work_j: f64,
    pub support_work_j: f64,
    pub support_impulse_ns: Vec3,
    pub support_angular_impulse_nms: Vec3,
    /// Total contact loss, including normal impulse/integration losses.
    /// No friction-heat partition is inferred.
    pub unallocated_contact_loss_j: f64,
    pub drift_energy_defect_j: f64,
    pub total_energy_residual_j: f64,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
}

/// Geometry refusal from an initial or complete-path query.
#[derive(Debug)]
pub enum SphereScenePathError {
    Penetration {
        contact: SphereSceneContactKey,
        gap_m: f64,
    },
    /// A complete rejected drift entered a surface. The suggested shorter
    /// duration is only a retry hint: force kick, solve and geometry all rerun.
    Crossing {
        contact: SphereSceneContactKey,
        gap_m: f64,
        suggested_duration_s: Option<f64>,
    },
    Dynamics(DynamicsError),
    BroadPhase(SweepError),
    PlaneQueryBudget,
    Invalid(&'static str),
    Cancelled,
}
impl From<DynamicsError> for SphereScenePathError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<SweepError> for SphereScenePathError {
    fn from(e: SweepError) -> Self {
        Self::BroadPhase(e)
    }
}
impl core::fmt::Display for SphereScenePathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sphere scene geometry: {self:?}")
    }
}
impl std::error::Error for SphereScenePathError {}

/// Every refusal preserves all caller-owned body states.
#[derive(Debug)]
pub enum SphereSceneError {
    Invalid(&'static str),
    MissingLaw(SphereSceneContactKey),
    Dynamics(DynamicsError),
    Impact(RigidFrictionalImpactError),
    Geometry(SphereScenePathError),
    Step(RigidContactStepError<SphereScenePathError>),
    StepBudget,
    RejectionBudget,
    NoProgress,
    Cancelled,
    Audit(&'static str),
}
impl From<DynamicsError> for SphereSceneError {
    fn from(e: DynamicsError) -> Self {
        Self::Dynamics(e)
    }
}
impl From<RigidFrictionalImpactError> for SphereSceneError {
    fn from(e: RigidFrictionalImpactError) -> Self {
        Self::Impact(e)
    }
}
impl From<SphereScenePathError> for SphereSceneError {
    fn from(e: SphereScenePathError) -> Self {
        if matches!(e, SphereScenePathError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Geometry(e)
        }
    }
}
impl core::fmt::Display for SphereSceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sphere scene: {self:?}")
    }
}
impl std::error::Error for SphereSceneError {}

#[derive(Default)]
struct Queries {
    sphere: usize,
    comparisons: usize,
    plane: usize,
}
impl Queries {
    fn plane(&mut self, c: SphereSceneControls) -> Result<(), SphereScenePathError> {
        if self.plane >= c.maximum_plane_queries {
            return Err(SphereScenePathError::PlaneQueryBudget);
        }
        self.plane += 1;
        Ok(())
    }
    fn sphere(&mut self) -> Result<(), SphereScenePathError> {
        self.sphere = self
            .sphere
            .checked_add(1)
            .ok_or(SphereScenePathError::Invalid("query count overflow"))?;
        Ok(())
    }
}
#[derive(Default)]
struct Exchange {
    force_work: f64,
    displacement_work: f64,
    splitting_work: f64,
    torque_work: f64,
    support_work: f64,
    loss: f64,
    drift_defect: f64,
    applied_linear: Vec3,
    applied_angular: Vec3,
    support_linear: Vec3,
    support_angular: Vec3,
}

/// Advance a complete sphere/plane scene through native simultaneous contact.
///
/// The law callback MUST be pure: return the same source-resolved coefficients
/// for the same key and do not publish external history. It is called only for
/// actual contacts within the activation distance, once per staged starting
/// configuration; None refuses rather than disabling collision geometry.
///
/// Restitution targets use the incoming point velocity BEFORE the force kick:
/// target = restitution * max(-incoming_normal_speed, 0). Thus a resting loaded
/// contact has target zero even when its material restitution is nonzero. The
/// shared law still admits against the actual kicked state and its passivity
/// gate; no positive work is silently authorized.
///
/// Initial/endpoint and whole straight drift paths use the true spherical
/// boundaries, independent of spin. An omitted crossing or negative native
/// endpoint retries a shorter interval from the same staged state. A crossing's
/// first-encounter time is only a proposal; every retry reruns the force kick,
/// contact solve and complete path admission. Other
/// refusals propagate. All accepted intervals remain private until the full
/// duration and global native momentum/work gates succeed.
///
/// The native free attitude integrator retains its asymmetric-body accuracy
/// limits. This does not provide arbitrary-shape CCD, pose repair, exact
/// impact-time resolution, finite-patch/rolling drag or an automatic time-error
/// estimate. Refine maximum_step_s and activation_gap_m for convergence.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_sphere_scene(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    planes: &[SphereScenePlane],
    wrenches: &[Wrench],
    gravity: Gravity,
    start_time_s: f64,
    duration_s: f64,
    controls: SphereSceneControls,
    law: impl Fn(SphereSceneContactKey) -> Option<SphereSceneLaw>,
    gate: &CancelGate,
) -> Result<SphereSceneReport, SphereSceneError> {
    poll(gate)?;
    let planes = prepare(
        bodies,
        properties,
        spheres,
        planes,
        wrenches,
        gravity,
        start_time_s,
        duration_s,
        controls,
        gate,
    )?;
    let (before, linear_before, angular_before) = super::totals(bodies, properties, gate)?;
    let mut staged = bodies.to_vec();
    let mut steps = Vec::new();
    let mut queries = Queries::default();
    let mut exchange = Exchange::default();
    let mut elapsed = 0.0;
    let mut rejected_steps = 0usize;
    let mut minimum_gap = f64::INFINITY;
    while elapsed < duration_s {
        poll(gate)?;
        if steps.len() >= controls.maximum_steps {
            return Err(SphereSceneError::StepBudget);
        }
        let time = start_time_s + elapsed;
        let current_planes = planes_at(&planes, time)?;
        let (contacts, keys, maximum_contact_gap) = detected_contacts(
            &staged,
            properties,
            spheres,
            &current_planes,
            controls,
            &law,
            &mut queries,
            gate,
        )?;
        let remaining = duration_s - elapsed;
        let mut h = remaining.min(controls.maximum_step_s);
        let native = loop {
            poll(gate)?;
            if h <= 0.0 || elapsed + h == elapsed || time + h == time {
                return Err(SphereSceneError::NoProgress);
            }
            let original = &staged;
            let (applied_linear, applied_angular) =
                applied_kick(original, properties, wrenches, gravity, h, gate)?;
            let result = step_rigid_frictional_contacts(
                &mut staged,
                properties,
                wrenches,
                gravity,
                &contacts,
                h,
                controls.step,
                |path| admit_path(path, spheres, &current_planes, controls, &mut queries, gate),
                gate,
            );
            match result {
                Ok(report) => {
                    exchange.applied_linear = exchange.applied_linear.add(applied_linear);
                    exchange.applied_angular = exchange.applied_angular.add(applied_angular);
                    break report;
                }
                Err(error) => {
                    if gate.is_requested() {
                        return Err(SphereSceneError::Cancelled);
                    }
                    let retry = matches!(
                        &error,
                        RigidContactStepError::Path(SphereScenePathError::Crossing { .. })
                            | RigidContactStepError::Penetration { endpoint: true, .. }
                    );
                    if !retry {
                        return Err(SphereSceneError::Step(error));
                    }
                    if rejected_steps >= controls.maximum_rejections {
                        return Err(SphereSceneError::RejectionBudget);
                    }
                    rejected_steps += 1;
                    let suggested = match &error {
                        RigidContactStepError::Path(SphereScenePathError::Crossing {
                            suggested_duration_s,
                            ..
                        }) => *suggested_duration_s,
                        _ => None,
                    };
                    h = match suggested.filter(|t| t.is_finite() && *t > 0.0) {
                        Some(t) => t.min(0.9 * h),
                        None => 0.5 * h,
                    };
                }
            }
        };
        elapsed = if h == remaining {
            duration_s
        } else {
            elapsed + h
        };
        if !elapsed.is_finite() || elapsed > duration_s {
            return Err(SphereSceneError::NoProgress);
        }
        minimum_gap = minimum_gap.min(native.path_admission.minimum_path_gap_m);
        exchange.force_work += native.force_impulse_work_j;
        exchange.displacement_work += native.force_displacement_work_j;
        exchange.splitting_work += native.force_splitting_work_j;
        exchange.torque_work += native.torque_impulse_work_j;
        exchange.support_work += native.impact.support_work_j;
        exchange.support_linear = exchange
            .support_linear
            .add(native.impact.support_impulse_ns);
        exchange.support_angular = exchange
            .support_angular
            .add(native.impact.support_angular_impulse_nms);
        exchange.loss += native.impact.contact.dissipated_energy;
        exchange.drift_defect += native.drift_energy_defect_j;
        steps.push(SphereSceneStep {
            elapsed_s: elapsed,
            contact_keys: keys,
            maximum_contact_gap_m: maximum_contact_gap,
            native,
        });
    }
    poll(gate)?;
    let (after, linear_after, angular_after) = super::totals(&staged, properties, gate)?;
    let energy_residual =
        after - before - exchange.force_work - exchange.torque_work - exchange.support_work
            + exchange.loss;
    let linear = linear_after
        .sub(linear_before)
        .sub(exchange.applied_linear)
        .sub(exchange.support_linear);
    let angular = angular_after
        .sub(angular_before)
        .sub(exchange.applied_angular)
        .sub(exchange.support_angular);
    if [
        energy_residual,
        exchange.force_work,
        exchange.displacement_work,
        exchange.splitting_work,
        exchange.torque_work,
        exchange.support_work,
        exchange.loss,
        exchange.drift_defect,
        super::magnitude(linear),
        super::magnitude(angular),
    ]
    .iter()
    .any(|x| !x.is_finite())
        || energy_residual.abs() > controls.step.work_energy_tolerance_j
        || super::magnitude(linear) > controls.step.linear_momentum_tolerance_ns
        || super::magnitude(angular) > controls.step.angular_momentum_tolerance_nms
    {
        return Err(SphereSceneError::Audit(
            "whole-call native momentum or energy",
        ));
    }
    poll(gate)?;
    bodies.copy_from_slice(&staged);
    Ok(SphereSceneReport {
        duration_s,
        steps,
        rejected_steps,
        sphere_pair_queries: queries.sphere,
        broad_phase_comparisons: queries.comparisons,
        plane_queries: queries.plane,
        minimum_path_gap_m: minimum_gap,
        kinetic_before_j: before,
        kinetic_after_j: after,
        force_impulse_work_j: exchange.force_work,
        force_displacement_work_j: exchange.displacement_work,
        force_splitting_work_j: exchange.splitting_work,
        torque_impulse_work_j: exchange.torque_work,
        support_work_j: exchange.support_work,
        support_impulse_ns: exchange.support_linear,
        support_angular_impulse_nms: exchange.support_angular,
        unallocated_contact_loss_j: exchange.loss,
        drift_energy_defect_j: exchange.drift_defect,
        total_energy_residual_j: energy_residual,
        linear_momentum_residual_ns: linear,
        angular_momentum_residual_nms: angular,
    })
}

#[allow(clippy::too_many_arguments)]
fn prepare(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    planes: &[SphereScenePlane],
    wrenches: &[Wrench],
    gravity: Gravity,
    start: f64,
    duration: f64,
    c: SphereSceneControls,
    gate: &CancelGate,
) -> Result<Vec<SphereScenePlane>, SphereSceneError> {
    gravity.validate()?;
    if bodies.is_empty()
        || bodies.len() != properties.len()
        || bodies.len() != spheres.len()
        || bodies.len() != wrenches.len()
        || bodies.len() > c.step.impact.maximum_bodies
        || c.maximum_steps == 0
        || c.maximum_steps > 65_536
        || c.maximum_rejections > 65_536
        || !c.maximum_step_s.is_finite()
        || c.maximum_step_s <= 0.0
        || !c.activation_gap_m.is_finite()
        || c.activation_gap_m <= 0.0
        || c.activation_gap_m > c.step.impact.geometry_tolerance_m
        || !start.is_finite()
        || !duration.is_finite()
        || duration <= 0.0
        || !(start + duration).is_finite()
    {
        return Err(SphereSceneError::Invalid(
            "matching bounded spheres, finite time and positive geometry/step controls required",
        ));
    }
    if start + duration == start {
        return Err(SphereSceneError::NoProgress);
    }
    let mut ids = Vec::with_capacity(spheres.len());
    for i in 0..bodies.len() {
        poll(gate)?;
        bodies[i].validate()?;
        properties[i].validate()?;
        wrenches[i].validate()?;
        let radius = spheres[i].radius_m;
        let inertia = properties[i].principal_inertia_body();
        let bound = (2.0 * properties[i].mass() * radius) * radius;
        if !radius.is_finite()
            || radius <= 0.0
            || !bound.is_finite()
            || inertia.x + inertia.y + inertia.z > bound
        {
            return Err(SphereSceneError::Invalid(
                "physical native inertia inside each sphere required",
            ));
        }
        ids.push(spheres[i].id);
    }
    ids.sort_unstable();
    if ids.windows(2).any(|x| x[0] == x[1]) {
        return Err(SphereSceneError::Invalid("sphere IDs must be distinct"));
    }
    if planes
        .len()
        .checked_mul(bodies.len())
        .is_none_or(|count| count > c.maximum_plane_queries)
    {
        return Err(SphereScenePathError::PlaneQueryBudget.into());
    }
    let mut prepared = planes.to_vec();
    prepared.sort_by_key(|p| p.id);
    for i in 0..prepared.len() {
        poll(gate)?;
        if i > 0 && prepared[i - 1].id == prepared[i].id {
            return Err(SphereSceneError::Invalid("plane IDs must be distinct"));
        }
        if !prepared[i].point_world.is_finite() || !prepared[i].velocity_world.is_finite() {
            return Err(SphereSceneError::Invalid(
                "finite plane point and velocity required",
            ));
        }
        prepared[i].normal_world = prepared[i].normal_world.normalized("scene plane normal")?;
    }
    Ok(prepared)
}

fn planes_at(
    planes: &[SphereScenePlane],
    time: f64,
) -> Result<Vec<SphereScenePlane>, SphereScenePathError> {
    planes
        .iter()
        .map(|p| {
            let point = p.point_world.add(p.velocity_world.scale(time));
            if !point.is_finite() {
                return Err(SphereScenePathError::Invalid("plane trajectory overflow"));
            }
            Ok(SphereScenePlane {
                point_world: point,
                ..*p
            })
        })
        .collect()
}

fn candidates(
    start: &[RigidBodyState],
    end: &[RigidBodyState],
    spheres: &[ContactSphere],
    padding: f64,
    c: SphereSceneControls,
    q: &mut Queries,
    gate: &CancelGate,
) -> Result<Vec<SweepPair>, SphereScenePathError> {
    let mut items = Vec::with_capacity(spheres.len());
    for i in 0..spheres.len() {
        geometry_poll(gate)?;
        items.push(SweepItem {
            id: spheres[i].id,
            bounds: swept_sphere_aabb(
                super::components(start[i].pose().position_world()),
                super::components(end[i].pose().position_world()),
                spheres[i].radius_m,
                padding,
            )?,
        });
    }
    let sweep = sweep_candidates(&items, c.broad_phase, || gate.is_requested())?;
    q.comparisons = q
        .comparisons
        .checked_add(sweep.checked_pairs)
        .ok_or(SphereScenePathError::Invalid("comparison count overflow"))?;
    Ok(sweep.pairs)
}

struct Detected {
    key: SphereSceneContactKey,
    point: RigidImpactPoint,
    other: RigidImpactCounterpart,
    normal: Vec3,
    gap: f64,
    incoming: f64,
}
#[allow(clippy::too_many_arguments)]
fn detected_contacts(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    spheres: &[ContactSphere],
    planes: &[SphereScenePlane],
    c: SphereSceneControls,
    law: &impl Fn(SphereSceneContactKey) -> Option<SphereSceneLaw>,
    q: &mut Queries,
    gate: &CancelGate,
) -> Result<(Vec<RigidFrictionalContact>, Vec<SphereSceneContactKey>, f64), SphereSceneError> {
    let mut detected = Vec::new();
    for pair in candidates(
        bodies,
        bodies,
        spheres,
        0.5 * c.activation_gap_m,
        c,
        q,
        gate,
    )? {
        geometry_poll(gate)?;
        q.sphere()?;
        let (a, b) = (pair.index_a, pair.index_b);
        let key = SphereSceneContactKey::Spheres {
            first: pair.id_a,
            second: pair.id_b,
        };
        let relative = bodies[a]
            .pose()
            .position_world()
            .sub(bodies[b].pose().position_world());
        let gap = checked_gap(
            relative.stable_norm("scene sphere distance")?
                - (spheres[a].radius_m + spheres[b].radius_m),
            key,
        )?;
        if gap > c.activation_gap_m {
            continue;
        }
        let normal = relative.normalized("scene sphere normal")?;
        let point = RigidImpactPoint {
            body: a,
            arm_body: bodies[a]
                .pose()
                .orientation()
                .rotate_world_to_body(normal.scale(-spheres[a].radius_m)),
        };
        let opposite = RigidImpactPoint {
            body: b,
            arm_body: bodies[b]
                .pose()
                .orientation()
                .rotate_world_to_body(normal.scale(spheres[b].radius_m)),
        };
        let va = bodies[a]
            .point_kinematics(properties[a], point.arm_body)?
            .point_velocity_world;
        let vb = bodies[b]
            .point_kinematics(properties[b], opposite.arm_body)?
            .point_velocity_world;
        detected.push(Detected {
            key,
            point,
            other: RigidImpactCounterpart::Body(opposite),
            normal,
            gap,
            incoming: va.sub(vb).dot(normal),
        });
        if detected.len() > c.step.impact.maximum_contacts {
            return Err(SphereSceneError::Invalid(
                "detected contact count exceeds native cap",
            ));
        }
    }
    for (i, sphere) in spheres.iter().enumerate() {
        for plane in planes {
            geometry_poll(gate)?;
            q.plane(c)?;
            let key = SphereSceneContactKey::SpherePlane {
                sphere: sphere.id,
                plane: plane.id,
            };
            let gap = checked_gap(
                bodies[i]
                    .pose()
                    .position_world()
                    .sub(plane.point_world)
                    .dot(plane.normal_world)
                    - sphere.radius_m,
                key,
            )?;
            if gap > c.activation_gap_m {
                continue;
            }
            let point = RigidImpactPoint {
                body: i,
                arm_body: bodies[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(plane.normal_world.scale(-sphere.radius_m)),
            };
            let velocity = bodies[i]
                .point_kinematics(properties[i], point.arm_body)?
                .point_velocity_world;
            detected.push(Detected {
                key,
                point,
                other: RigidImpactCounterpart::PrescribedPlane {
                    point_world: plane.point_world,
                    velocity_world: plane.velocity_world,
                },
                normal: plane.normal_world,
                gap,
                incoming: velocity.sub(plane.velocity_world).dot(plane.normal_world),
            });
            if detected.len() > c.step.impact.maximum_contacts {
                return Err(SphereSceneError::Invalid(
                    "detected contact count exceeds native cap",
                ));
            }
        }
    }
    detected.sort_by_key(|d| d.key);
    let mut contacts = Vec::with_capacity(detected.len());
    let mut keys = Vec::with_capacity(detected.len());
    let mut maximum_gap = 0.0_f64;
    for (id, d) in detected.into_iter().enumerate() {
        poll(gate)?;
        let resolved = law(d.key).ok_or(SphereSceneError::MissingLaw(d.key))?;
        if !resolved.friction_coefficient.is_finite()
            || resolved.friction_coefficient < 0.0
            || !resolved.restitution.is_finite()
            || !(0.0..=1.0).contains(&resolved.restitution)
            || !d.incoming.is_finite()
        {
            return Err(SphereSceneError::Invalid(
                "finite source-resolved contact coefficients and velocity required",
            ));
        }
        let axis = if d.normal.x.abs() < 0.8 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        contacts.push(RigidFrictionalContact {
            id: id as u64,
            point: d.point,
            other: d.other,
            normal_world: d.normal,
            tangent_world: d.normal.cross(axis),
            activation_tolerance_m: c.activation_gap_m,
            friction_coefficient: resolved.friction_coefficient,
            normal_law: NormalVelocityLaw::Target {
                speed: -resolved.restitution * d.incoming.min(0.0),
            },
        });
        keys.push(d.key);
        maximum_gap = maximum_gap.max(d.gap);
    }
    Ok((contacts, keys, maximum_gap))
}

fn admit_path(
    path: &RigidContactPath<'_>,
    spheres: &[ContactSphere],
    planes: &[SphereScenePlane],
    c: SphereSceneControls,
    q: &mut Queries,
    gate: &CancelGate,
) -> Result<SphereScenePathReport, SphereScenePathError> {
    let mut minimum = f64::INFINITY;
    let mut sphere_pairs_checked = 0usize;
    let mut plane_pairs_checked = 0usize;
    for pair in candidates(path.post_contact, path.endpoint, spheres, 0.0, c, q, gate)? {
        geometry_poll(gate)?;
        q.sphere()?;
        sphere_pairs_checked += 1;
        let (a, b) = (pair.index_a, pair.index_b);
        let key = SphereSceneContactKey::Spheres {
            first: pair.id_a,
            second: pair.id_b,
        };
        let r = path.post_contact[a]
            .pose()
            .position_world()
            .sub(path.post_contact[b].pose().position_world());
        let end = path.endpoint[a]
            .pose()
            .position_world()
            .sub(path.endpoint[b].pose().position_world());
        let va = path.post_contact[a].center_of_mass_velocity_world(path.properties[a])?;
        let vb = path.post_contact[b].center_of_mass_velocity_world(path.properties[b])?;
        let d = va.sub(vb).scale(path.duration_s);
        let scale = r
            .stable_norm("scene path separation")?
            .max(d.stable_norm("scene path displacement")?);
        let fraction = if scale == 0.0 {
            0.0
        } else {
            let rs = Vec3::new(r.x / scale, r.y / scale, r.z / scale);
            let ds = Vec3::new(d.x / scale, d.y / scale, d.z / scale);
            let square = ds.norm_squared();
            if square == 0.0 {
                0.0
            } else {
                (-rs.dot(ds) / square).clamp(0.0, 1.0)
            }
        };
        let radius = spheres[a].radius_m + spheres[b].radius_m;
        let interior = r
            .add(d.scale(fraction))
            .stable_norm("scene closest approach")?
            - radius;
        let initial = checked_gap(r.stable_norm("scene initial separation")? - radius, key)?;
        let endpoint = end.stable_norm("scene final separation")? - radius;
        if !interior.is_finite() || !endpoint.is_finite() {
            return Err(SphereScenePathError::Invalid("nonfinite sphere path gap"));
        }
        let gap = initial.min(interior).min(endpoint);
        if gap < 0.0 {
            let suggested = first_sphere_contact(SweptSphereInput {
                relative_position_world: super::components(r),
                relative_velocity_world: super::components(va.sub(vb)),
                radius_sum: radius,
                duration_s: path.duration_s,
            })
            .ok()
            .flatten()
            .map(|hit| hit.time_s.next_down());
            return Err(SphereScenePathError::Crossing {
                contact: key,
                gap_m: gap,
                suggested_duration_s: suggested,
            });
        }
        minimum = minimum.min(gap);
    }
    for (i, sphere) in spheres.iter().enumerate() {
        for plane in planes {
            geometry_poll(gate)?;
            q.plane(c)?;
            plane_pairs_checked += 1;
            let key = SphereSceneContactKey::SpherePlane {
                sphere: sphere.id,
                plane: plane.id,
            };
            let plane_end = plane
                .point_world
                .add(plane.velocity_world.scale(path.duration_s));
            if !plane_end.is_finite() {
                return Err(SphereScenePathError::Invalid("plane path overflow"));
            }
            // True bottom gap is affine along the complete straight COM drift,
            // regardless of the old material witness or sphere spin.
            let initial = checked_gap(
                path.post_contact[i]
                    .pose()
                    .position_world()
                    .sub(plane.point_world)
                    .dot(plane.normal_world)
                    - sphere.radius_m,
                key,
            )?;
            let endpoint = path.endpoint[i]
                .pose()
                .position_world()
                .sub(plane_end)
                .dot(plane.normal_world)
                - sphere.radius_m;
            if !endpoint.is_finite() {
                return Err(SphereScenePathError::Invalid("nonfinite plane path gap"));
            }
            if endpoint < 0.0 {
                let speed = path.post_contact[i]
                    .center_of_mass_velocity_world(path.properties[i])?
                    .sub(plane.velocity_world)
                    .dot(plane.normal_world);
                let target = 0.5 * c.activation_gap_m;
                let suggested = if speed < 0.0 {
                    Some(if initial > target {
                        (initial - target) / (-speed)
                    } else {
                        (initial / (-speed)).next_down()
                    })
                } else {
                    None
                };
                return Err(SphereScenePathError::Crossing {
                    contact: key,
                    gap_m: endpoint,
                    suggested_duration_s: suggested,
                });
            }
            minimum = minimum.min(initial).min(endpoint);
        }
    }
    Ok(SphereScenePathReport {
        minimum_path_gap_m: minimum,
        sphere_pairs_checked,
        plane_pairs_checked,
    })
}

fn checked_gap(gap: f64, contact: SphereSceneContactKey) -> Result<f64, SphereScenePathError> {
    if !gap.is_finite() {
        return Err(SphereScenePathError::Invalid("nonfinite geometric gap"));
    }
    if gap < 0.0 {
        return Err(SphereScenePathError::Penetration {
            contact,
            gap_m: gap,
        });
    }
    Ok(gap)
}
fn applied_kick(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    wrenches: &[Wrench],
    gravity: Gravity,
    h: f64,
    gate: &CancelGate,
) -> Result<(Vec3, Vec3), SphereSceneError> {
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    for ((body, property), wrench) in bodies.iter().zip(properties).zip(wrenches) {
        poll(gate)?;
        let impulse = wrench
            .force_world
            .add(gravity.acceleration_world().scale(property.mass()))
            .scale(h);
        linear = linear.add(impulse);
        angular = angular
            .add(body.pose().position_world().cross(impulse))
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(wrench.torque_body.scale(h)),
            );
    }
    if !linear.is_finite() || !angular.is_finite() {
        return Err(SphereSceneError::Audit("applied kick range"));
    }
    Ok((linear, angular))
}
fn geometry_poll(gate: &CancelGate) -> Result<(), SphereScenePathError> {
    if gate.is_requested() {
        Err(SphereScenePathError::Cancelled)
    } else {
        Ok(())
    }
}
fn poll(gate: &CancelGate) -> Result<(), SphereSceneError> {
    if gate.is_requested() {
        Err(SphereSceneError::Cancelled)
    } else {
        Ok(())
    }
}
