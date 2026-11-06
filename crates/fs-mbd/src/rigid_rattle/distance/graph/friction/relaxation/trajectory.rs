//! Moving articulated mechanisms through the original native pose and friction owners.
//!
//! Each split advances the original constrained rigid graph for half a step,
//! resolves the coupled continuous Coulomb flow at that moved geometry, then
//! advances the graph for the other half. The flow's travel is a work integral;
//! it never becomes a second pose update. Repeated step doubling controls the
//! differences in actual poses, momenta, per-port travel, actions and work.
//! These are numerical estimates, not certified error bounds or event discovery.

use super::{
    ArticulatedJoint, Dd, DistanceGraphError, DistanceStepError, JointFrictionFlowControls,
    JointFrictionFlowError, JointFrictionFlowMotion, JointFrictionFlowPort,
    MAX_JOINT_FRICTION_PORTS, MassProperties, RigidBodyState, Vec3, advance_joint_friction,
};
use crate::rigid_rattle::distance::graph::advance_articulated_graph;
use crate::{DynamicsError, RigidBodyIntegrator, Wrench};
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowError;
use fs_math::det;

/// Absolute differences accumulated from one split versus two half splits.
/// Controls use these same units as a budget over the complete duration. Hinge
/// travel is multiplied by its original angular length; axial impulse is divided
/// by that length. Every body and every selected port participates in the maxima.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct JointFrictionTrajectoryAccuracy {
    pub position_m: f64,
    pub orientation_rad: f64,
    pub linear_momentum_ns: f64,
    /// World spin momentum, compared in a common frame at the actual endpoints.
    pub angular_momentum_nms: f64,
    pub scaled_travel_m: f64,
    pub scaled_impulse_ns: f64,
    pub work_j: f64,
}
impl JointFrictionTrajectoryAccuracy {
    fn values(self) -> [f64; 7] {
        [
            self.position_m,
            self.orientation_rad,
            self.linear_momentum_ns,
            self.angular_momentum_nms,
            self.scaled_travel_m,
            self.scaled_impulse_ns,
            self.work_j,
        ]
    }
    fn from_values(x: [f64; 7]) -> Self {
        Self {
            position_m: x[0],
            orientation_rad: x[1],
            linear_momentum_ns: x[2],
            angular_momentum_nms: x[3],
            scaled_travel_m: x[4],
            scaled_impulse_ns: x[5],
            work_j: x[6],
        }
    }
}

/// Original native gates, bounded refinement and explicit numerical work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionTrajectoryControls {
    pub mechanical: JointFrictionFlowControls,
    /// Maximum compared interval, clipped at each authored source breakpoint.
    pub maximum_step_s: f64,
    /// Smallest individual split duration, including the accepted half splits.
    pub minimum_step_s: f64,
    /// Number of accepted individual splits; each accepted comparison uses two.
    pub maximum_accepted_steps: usize,
    /// All split attempts, including discarded coarse paths and failed trials.
    pub maximum_split_trials: usize,
    /// Native workspace plus this owner's endpoint buffers and final port table.
    /// Borrowed inputs, fixed stack banks and callback-owned memory are excluded.
    pub maximum_live_bytes: usize,
    /// Raw coarse/fine differences, with no assumed-order Richardson divisor.
    /// An interval h receives h/total_duration of each strictly positive budget.
    pub accuracy: JointFrictionTrajectoryAccuracy,
    /// Absolute complete-trajectory energy residual, never converted to heat.
    pub energy_residual_tolerance_j: f64,
    /// Sum of absolute native graph/flow constraint work over accepted splits.
    pub constraint_work_tolerance_j: f64,
}

/// Integrals along the successive, recomputed axial phase paths. Their directions
/// follow the current physical hinge or moving rail, so these scalar actions are
/// not world-vector impulse totals. Pose/travel agreement is obtained by refinement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JointFrictionTrajectoryMotion {
    Revolute {
        angular_length_m: f64,
        rotation_rad: f64,
        angular_distance_rad: f64,
        friction_impulse_nms: f64,
        drive_impulse_nms: f64,
    },
    Prismatic {
        displacement_m: f64,
        distance_m: f64,
        friction_impulse_ns: f64,
        drive_impulse_ns: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JointFrictionTrajectoryPort {
    /// Original joint index, in increasing order.
    pub joint_index: usize,
    pub motion: JointFrictionTrajectoryMotion,
    pub drive_work_j: f64,
    /// Nonpositive resolved mechanical work; no heat fraction is inferred.
    pub friction_work_j: f64,
}

/// One atomic candidate; borrowed accepted bodies and sources are never mutated.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct JointFrictionTrajectoryReport {
    pub states_after: Vec<RigidBodyState>,
    pub start_time_s: f64,
    pub target_time_s: f64,
    pub duration_s: f64,
    pub ports: Vec<JointFrictionTrajectoryPort>,
    /// Actual kinetic plus original uniform-gravity potential change.
    pub mechanical_energy_change_j: f64,
    /// Original free-receipt force/displacement and body torque/midpoint work.
    pub wrench_work_j: f64,
    pub drive_work_j: f64,
    pub dissipated_j: f64,
    pub constraint_work_j: f64,
    pub absolute_constraint_work_j: f64,
    /// Mechanical change - wrench work - drive work + loss - constraint work.
    pub work_energy_residual_j: f64,
    pub estimated_error: JointFrictionTrajectoryAccuracy,
    pub accepted_steps: usize,
    pub split_trials: usize,
    pub rejected_intervals: usize,
    pub planned_live_bytes: usize,
}

#[derive(Debug)]
pub enum JointFrictionTrajectoryError<E> {
    Invalid(&'static str),
    Native(DistanceGraphError),
    Friction(JointFrictionFlowError),
    Source(E),
    Limit(&'static str),
    Allocation,
    Unresolved(&'static str),
    Cancelled,
}
impl<E> From<DynamicsError> for JointFrictionTrajectoryError<E> {
    fn from(error: DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl<E> From<DistanceGraphError> for JointFrictionTrajectoryError<E> {
    fn from(error: DistanceGraphError) -> Self {
        Self::Native(error)
    }
}
impl<E> From<JointFrictionFlowError> for JointFrictionTrajectoryError<E> {
    fn from(error: JointFrictionFlowError) -> Self {
        Self::Friction(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for JointFrictionTrajectoryError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native articulated friction trajectory: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for JointFrictionTrajectoryError<E> {}

#[derive(Clone, Copy)]
struct Port {
    joint: usize,
    // Zero denotes a linear port; positive is the original hinge length.
    length: f64,
    // Signed travel, absolute travel, friction impulse, drive impulse, drive work, friction work.
    sums: [Dd; 6],
}
impl Port {
    const EMPTY: Self = Self {
        joint: 0,
        length: 0.0,
        sums: [Dd::ZERO; 6],
    };
}
#[derive(Clone, Copy)]
struct Account {
    ports: [Port; MAX_JOINT_FRICTION_PORTS],
    count: usize,
    wrench: Dd,
    constraint: Dd,
    absolute_constraint: Dd,
}
impl Account {
    const EMPTY: Self = Self {
        ports: [Port::EMPTY; MAX_JOINT_FRICTION_PORTS],
        count: 0,
        wrench: Dd::ZERO,
        constraint: Dd::ZERO,
        absolute_constraint: Dd::ZERO,
    };
    fn merge<E>(&mut self, other: Self) -> Result<(), JointFrictionTrajectoryError<E>> {
        if self.count == 0 {
            self.ports = other.ports;
            self.count = other.count;
        } else {
            if self.count != other.count {
                return Err(JointFrictionTrajectoryError::Invalid(
                    "changing friction ports",
                ));
            }
            for (a, b) in self.ports[..self.count]
                .iter_mut()
                .zip(&other.ports[..other.count])
            {
                if a.joint != b.joint || a.length != b.length {
                    return Err(JointFrictionTrajectoryError::Invalid(
                        "changing friction port identity",
                    ));
                }
                for (x, y) in a.sums.iter_mut().zip(b.sums) {
                    *x = *x + y;
                }
            }
        }
        self.wrench = self.wrench + other.wrench;
        self.constraint = self.constraint + other.constraint;
        self.absolute_constraint = self.absolute_constraint + other.absolute_constraint;
        Ok(())
    }
    fn work(self) -> [Dd; 4] {
        let mut values = [self.wrench, Dd::ZERO, Dd::ZERO, self.constraint];
        for p in &self.ports[..self.count] {
            values[1] = values[1] + p.sums[4];
            values[2] = values[2] - p.sums[5];
        }
        values
    }
}
struct Candidate {
    states: Vec<RigidBodyState>,
    account: Account,
}

fn checkpoint<E>(
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), JointFrictionTrajectoryError<E>> {
    if cancelled() {
        Err(JointFrictionTrajectoryError::Cancelled)
    } else {
        Ok(())
    }
}
fn norm(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |s, (a, b)| s + Dd::from_f64(a) * Dd::from_f64(b))
}

#[allow(clippy::too_many_arguments)]
fn smooth<E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: JointFrictionFlowControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, JointFrictionTrajectoryError<E>> {
    let out = advance_articulated_graph(
        states,
        properties,
        joints,
        integrator,
        wrenches,
        duration_s,
        controls.graph,
        &mut *cancelled,
    )?;
    let mut account = Account::EMPTY;
    account.constraint = Dd::from_f64(out.constraint_work_j);
    account.absolute_constraint = Dd::from_f64(out.constraint_work_j.abs());
    for ((step, properties), wrench) in out.step.free_steps.iter().zip(properties).zip(wrenches) {
        checkpoint(cancelled)?;
        let before = step.state_before;
        let after = step.state_after;
        let linear = components(after.pose().position_world())
            .into_iter()
            .zip(components(before.pose().position_world()))
            .zip(components(wrench.force_world))
            .fold(Dd::ZERO, |sum, ((a, b), force)| {
                sum + (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(force)
            });
        let before_w = properties.angular_velocity_body_checked(before.angular_momentum_body())?;
        let after_w = properties.angular_velocity_body_checked(after.angular_momentum_body())?;
        let rotational = (dot(wrench.torque_body, before_w) + dot(wrench.torque_body, after_w))
            * Dd::from_f64(0.5 * step.duration_seconds);
        account.wrench = account.wrench + linear + rotational;
    }
    if !account.wrench.to_f64().is_finite() {
        return Err(JointFrictionTrajectoryError::Unresolved(
            "finite native wrench work",
        ));
    }
    Ok(Candidate {
        states: out.step.states_after,
        account,
    })
}

#[allow(clippy::too_many_arguments)]
fn split<E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[JointFrictionFlowPort],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    start_s: f64,
    duration_s: f64,
    controls: JointFrictionTrajectoryControls,
    resolve: &impl Fn(f64, &[RigidBodyState], &mut [JointFrictionFlowPort]) -> Result<(), E>,
    trials: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, JointFrictionTrajectoryError<E>> {
    checkpoint(cancelled)?;
    if *trials >= controls.maximum_split_trials {
        return Err(JointFrictionTrajectoryError::Limit("split trials"));
    }
    *trials += 1;
    let half = 0.5 * duration_s;
    let midpoint = start_s + half;
    let end = start_s + duration_s;
    if duration_s < controls.minimum_step_s
        || !(half > 0.0 && midpoint > start_s && midpoint < end && end.is_finite())
    {
        return Err(JointFrictionTrajectoryError::Limit(
            "representable minimum split duration",
        ));
    }
    let first = smooth(
        states,
        properties,
        joints,
        integrator,
        wrenches,
        half,
        controls.mechanical,
        cancelled,
    )?;
    let mut authored = [sources[0]; MAX_JOINT_FRICTION_PORTS];
    authored[..sources.len()].copy_from_slice(sources);
    checkpoint(cancelled)?;
    resolve(midpoint, &first.states, &mut authored[..sources.len()])
        .map_err(JointFrictionTrajectoryError::Source)?;
    checkpoint(cancelled)?;
    for (&original, &trial) in sources.iter().zip(&authored) {
        if original.joint_index() != trial.joint_index()
            || core::mem::discriminant(&original) != core::mem::discriminant(&trial)
        {
            return Err(JointFrictionTrajectoryError::Invalid(
                "source changed joint identity",
            ));
        }
    }
    let (flow_states, mut account) = {
        let flow = advance_joint_friction(
            &first.states,
            properties,
            joints,
            &authored[..sources.len()],
            duration_s,
            controls.mechanical,
            &mut *cancelled,
        )?;
        let mut account = Account::EMPTY;
        account.count = flow.ports.len();
        account.wrench = first.account.wrench;
        account.constraint = first.account.constraint + Dd::from_f64(flow.constraint_work_j);
        account.absolute_constraint =
            first.account.absolute_constraint + Dd::from_f64(flow.constraint_work_j.abs());
        for (slot, p) in account.ports.iter_mut().zip(&flow.ports) {
            let (length, travel, distance, friction, drive) = match p.motion {
                JointFrictionFlowMotion::Revolute {
                    angular_length_m,
                    rotation_rad,
                    angular_distance_rad,
                    friction_impulse_nms,
                    drive_impulse_nms,
                    ..
                } => (
                    angular_length_m,
                    rotation_rad,
                    angular_distance_rad,
                    friction_impulse_nms,
                    drive_impulse_nms,
                ),
                JointFrictionFlowMotion::Prismatic {
                    displacement_m,
                    distance_m,
                    friction_impulse_ns,
                    drive_impulse_ns,
                    ..
                } => (
                    0.0,
                    displacement_m,
                    distance_m,
                    friction_impulse_ns,
                    drive_impulse_ns,
                ),
            };
            *slot = Port {
                joint: p.joint_index,
                length,
                sums: [
                    travel,
                    distance,
                    friction,
                    drive,
                    p.drive_work_j,
                    p.friction_work_j,
                ]
                .map(Dd::from_f64),
            };
        }
        (flow.states_after, account)
    };
    drop(first);
    let second = smooth(
        &flow_states,
        properties,
        joints,
        integrator,
        wrenches,
        half,
        controls.mechanical,
        cancelled,
    )?;
    account.wrench = account.wrench + second.account.wrench;
    account.constraint = account.constraint + second.account.constraint;
    account.absolute_constraint = account.absolute_constraint + second.account.absolute_constraint;
    checkpoint(cancelled)?;
    Ok(Candidate {
        states: second.states,
        account,
    })
}

fn difference<E>(
    coarse: &Candidate,
    fine: &Candidate,
) -> Result<[f64; 7], JointFrictionTrajectoryError<E>> {
    let mut error = [0.0_f64; 7];
    for (&a, &b) in coarse.states.iter().zip(&fine.states) {
        error[0] = error[0].max(norm(
            a.pose().position_world().sub(b.pose().position_world()),
        ));
        let qa = a.pose().orientation().components();
        let qb = b.pose().orientation().components();
        let mut minus = 0.0;
        let mut plus = 0.0;
        for i in 0..4 {
            minus = det::hypot(minus, qa[i] - qb[i]);
            plus = det::hypot(plus, qa[i] + qb[i]);
        }
        error[1] = error[1].max(4.0 * det::asin(0.5 * minus.min(plus)));
        error[2] = error[2].max(norm(
            a.linear_momentum_world().sub(b.linear_momentum_world()),
        ));
        let world = |s: RigidBodyState| {
            s.pose()
                .orientation()
                .rotate_body_to_world(s.angular_momentum_body())
        };
        error[3] = error[3].max(norm(world(a).sub(world(b))));
    }
    if coarse.account.count != fine.account.count {
        return Err(JointFrictionTrajectoryError::Invalid(
            "changing friction ports",
        ));
    }
    for (a, b) in coarse.account.ports[..coarse.account.count]
        .iter()
        .zip(&fine.account.ports)
    {
        if a.joint != b.joint || a.length != b.length {
            return Err(JointFrictionTrajectoryError::Invalid(
                "changing friction port identity",
            ));
        }
        let length = if a.length > 0.0 { a.length } else { 1.0 };
        for i in 0..6 {
            let gap = (a.sums[i] - b.sums[i]).to_f64().abs();
            let (slot, scaled) = if i < 2 {
                (4, gap * length)
            } else if i < 4 {
                (5, gap / length)
            } else {
                (6, gap)
            };
            error[slot] = error[slot].max(scaled);
        }
    }
    for (a, b) in coarse.account.work().into_iter().zip(fine.account.work()) {
        error[6] = error[6].max((a - b).to_f64().abs());
    }
    if error.iter().any(|x| !x.is_finite()) {
        return Err(JointFrictionTrajectoryError::Unresolved(
            "finite split comparison",
        ));
    }
    Ok(error)
}
fn retryable<E>(error: &JointFrictionTrajectoryError<E>) -> bool {
    matches!(
        error,
        JointFrictionTrajectoryError::Native(DistanceGraphError::Step(
            DistanceStepError::LocalBranch | DistanceStepError::NotConverged
        )) | JointFrictionTrajectoryError::Friction(JointFrictionFlowError::Flow(
            CoupledFrictionFlowError::ImpulseBudget { .. }
                | CoupledFrictionFlowError::SegmentLimit { .. }
        ))
    )
}

/// Advance native articulated poses and momenta with changing axial Coulomb loads.
///
/// `resolve` receives the physical midpoint time, actual moved native states and
/// a fixed stack copy of `sources`. It may change capacities and drives, retaining
/// each port's joint index and kind. Trials repeat and may be rejected: sources
/// must be deterministic, stage their work, and publish no external state. Loads
/// are authored or obtained from an existing physical source; aggregate bilateral
/// reactions do not determine a bearing's normal pressure/preload distribution.
///
/// External world forces, body torques and original uniform gravity are held over
/// this transaction. Axial loads that must participate in exact static balance
/// belong in the drive ports: splitting a separate body wrench against dry
/// stiction can give first-order endpoint velocity errors and exhaust refinement. The source drives must not duplicate those wrenches. Known
/// source discontinuities MUST appear in `source_breakpoints_s`, as strictly
/// increasing absolute times inside the interval; no comparison crosses them.
/// The caller must bound the maximum step to resolve other source features.
/// Step doubling cannot detect an unknown pulse skipped by all source samples.
///
/// The represented endpoint `start_time_s + duration_s` fixes the actual elapsed
/// duration reported and integrated; midpoint times must remain distinct.
///
/// Native half steps update poses once in total; friction path travel is used
/// only for work and exported travel. Stopping can lower the splitting order,
/// so the raw one-versus-two difference is used without an assumed-order factor.
/// Only the two-half-step candidate is accumulated. All refusals, including a
/// late source failure or cancellation, preserve the original accepted bodies.
/// Mechanical loss and numerical constraint/balance work remain separate. The
/// supplied axial model may declare distinct static and kinetic capacities; it
/// is an ideal discontinuous force law, not compliant presliding history. Each
/// phase holds only an exact resting input or an isolated zero event. Finite
/// native roundoff or source-created nonzero velocity is not silently classified
/// as sticking. Persistent geometric locking needs its constraint/history owner.
/// No thermal partition, wear or unprovided friction history is inferred.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_joint_friction_trajectory<E>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    sources: &[JointFrictionFlowPort],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    start_time_s: f64,
    duration_s: f64,
    source_breakpoints_s: &[f64],
    controls: JointFrictionTrajectoryControls,
    resolve: impl Fn(f64, &[RigidBodyState], &mut [JointFrictionFlowPort]) -> Result<(), E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<JointFrictionTrajectoryReport, JointFrictionTrajectoryError<E>> {
    checkpoint(&mut cancelled)?;
    let target = start_time_s + duration_s;
    if !start_time_s.is_finite()
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !target.is_finite()
        || target <= start_time_s
        || states.is_empty()
        || states.len() != properties.len()
        || states.len() != wrenches.len()
        || states.len() > controls.mechanical.graph.maximum_bodies
        || sources.is_empty()
        || sources.len() > MAX_JOINT_FRICTION_PORTS
        || !controls.minimum_step_s.is_finite()
        || controls.minimum_step_s <= 0.0
        || !controls.maximum_step_s.is_finite()
        || controls.maximum_step_s < 2.0 * controls.minimum_step_s
        || controls.maximum_accepted_steps < 2
        || controls.maximum_split_trials < 3
        || controls
            .accuracy
            .values()
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        || !controls.energy_residual_tolerance_j.is_finite()
        || controls.energy_residual_tolerance_j < 0.0
        || !controls.constraint_work_tolerance_j.is_finite()
        || controls.constraint_work_tolerance_j < 0.0
    {
        return Err(JointFrictionTrajectoryError::Invalid(
            "trajectory controls, topology or times",
        ));
    }
    // The native clock is f64: the represented endpoint fixes the physical
    // interval, including every force integral and local refinement fraction.
    let duration_s = target - start_time_s;
    if source_breakpoints_s.len() >= controls.maximum_accepted_steps / 2 {
        return Err(JointFrictionTrajectoryError::Limit(
            "source breakpoint steps",
        ));
    }
    let mut previous = start_time_s;
    for &time in source_breakpoints_s {
        checkpoint(&mut cancelled)?;
        if !time.is_finite() || time <= previous || time >= target {
            return Err(JointFrictionTrajectoryError::Invalid(
                "ordered interior source breakpoints",
            ));
        }
        previous = time;
    }
    // At most one native graph/flow workspace is live at once. The factor two
    // covers their existing independently bounded payloads; eight body vectors
    // cover retained coarse/fine/current endpoints and native call inputs.
    let bytes = 2 * controls.mechanical.graph.maximum_live_bytes as u128
        + 8 * states.len() as u128 * core::mem::size_of::<RigidBodyState>() as u128
        + sources.len() as u128 * core::mem::size_of::<JointFrictionTrajectoryPort>() as u128;
    if bytes > controls.maximum_live_bytes as u128 {
        return Err(JointFrictionTrajectoryError::Limit("trajectory live bytes"));
    }
    let mut current = Vec::new();
    current
        .try_reserve_exact(states.len())
        .map_err(|_| JointFrictionTrajectoryError::Allocation)?;
    current.extend_from_slice(states);
    let mut account = Account::EMPTY;
    let mut estimated = [Dd::ZERO; 7];
    let budget = controls.accuracy.values();
    let mut elapsed = 0.0;
    let mut breakpoint = 0;
    let mut proposed = controls.maximum_step_s;
    let mut trials = 0;
    let mut accepted = 0;
    let mut rejected = 0;
    while elapsed < duration_s {
        checkpoint(&mut cancelled)?;
        if accepted > controls.maximum_accepted_steps - 2 {
            return Err(JointFrictionTrajectoryError::Limit(
                "accepted trajectory steps",
            ));
        }
        let boundary = source_breakpoints_s
            .get(breakpoint)
            .map_or(duration_s, |t| *t - start_time_s);
        let remaining = boundary - elapsed;
        let mut h = proposed.min(remaining);
        let minimum_comparison = 2.0 * controls.minimum_step_s;
        // Do not strand a sub-minimum tail, including a decimal endpoint's
        // accumulation residue. Integrate the entire remainder when admitted,
        // otherwise divide it into two admissible intervals. Nothing is snapped
        // or omitted, and a rejected sub-minimum proposal is never enlarged.
        if h >= minimum_comparison && remaining > h && remaining - h < minimum_comparison {
            h = if remaining <= controls.maximum_step_s {
                remaining
            } else {
                0.5 * remaining
            };
        }
        if !h.is_finite() || h < minimum_comparison || elapsed + h <= elapsed {
            return Err(JointFrictionTrajectoryError::Limit(
                "minimum comparison duration",
            ));
        }
        let time = start_time_s + elapsed;
        let result = (|| {
            let coarse = split(
                &current,
                properties,
                joints,
                sources,
                integrator,
                wrenches,
                time,
                h,
                controls,
                &resolve,
                &mut trials,
                &mut cancelled,
            )?;
            let first = split(
                &current,
                properties,
                joints,
                sources,
                integrator,
                wrenches,
                time,
                0.5 * h,
                controls,
                &resolve,
                &mut trials,
                &mut cancelled,
            )?;
            let mut fine = split(
                &first.states,
                properties,
                joints,
                sources,
                integrator,
                wrenches,
                time + 0.5 * h,
                0.5 * h,
                controls,
                &resolve,
                &mut trials,
                &mut cancelled,
            )?;
            let mut combined = first.account;
            combined.merge(fine.account)?;
            fine.account = combined;
            let error = difference(&coarse, &fine)?;
            Ok::<_, JointFrictionTrajectoryError<E>>((fine, error))
        })();
        let (fine, error) = match result {
            Ok(candidate) => candidate,
            Err(error) if retryable(&error) => {
                rejected += 1;
                proposed = 0.5 * h;
                continue;
            }
            Err(error) => return Err(error),
        };
        let fraction = h / duration_s;
        if fraction <= 0.0 || !fraction.is_finite() {
            return Err(JointFrictionTrajectoryError::Invalid(
                "representable refinement fraction",
            ));
        }
        let mut ratio = 0.0_f64;
        for i in 0..7 {
            let allowance = budget[i] * fraction;
            if allowance <= 0.0 || !allowance.is_finite() {
                return Err(JointFrictionTrajectoryError::Invalid(
                    "representable local error allowance",
                ));
            }
            ratio = ratio.max(error[i] / allowance);
        }
        let constraint_limit = controls.constraint_work_tolerance_j * fraction;
        if ratio > 1.0 || fine.account.absolute_constraint.to_f64() > constraint_limit {
            rejected += 1;
            proposed = 0.5 * h;
            continue;
        }
        for i in 0..7 {
            estimated[i] = estimated[i] + Dd::from_f64(error[i]);
        }
        account.merge(fine.account)?;
        current = fine.states;
        accepted += 2;
        if h == remaining {
            elapsed = boundary;
            if breakpoint < source_breakpoints_s.len() {
                breakpoint += 1;
            }
        } else {
            elapsed += h;
        }
        proposed = if ratio < 0.0625 {
            (2.0 * h).min(controls.maximum_step_s)
        } else {
            h
        };
    }
    let mut energy = Dd::ZERO;
    for ((&before, &after), &p) in states.iter().zip(&current).zip(properties) {
        checkpoint(&mut cancelled)?;
        energy = energy + Dd::from_f64(integrator.diagnostics(after, p)?.mechanical_energy)
            - Dd::from_f64(integrator.diagnostics(before, p)?.mechanical_energy);
    }
    let work = account.work();
    let residual = (energy - work[0] - work[1] + work[2] - work[3]).to_f64();
    let estimates = estimated.map(Dd::to_f64);
    let absolute_constraint = account.absolute_constraint.to_f64();
    if [energy.to_f64(), residual, absolute_constraint]
        .iter()
        .any(|x| !x.is_finite())
        || work.iter().any(|x| !x.to_f64().is_finite())
        || work[2].to_f64() < 0.0
        || residual.abs() > controls.energy_residual_tolerance_j
        || absolute_constraint > controls.constraint_work_tolerance_j
        || estimates
            .iter()
            .zip(budget)
            .any(|(&e, b)| !e.is_finite() || e > b)
    {
        return Err(JointFrictionTrajectoryError::Unresolved(
            "original trajectory work or error budgets",
        ));
    }
    let mut ports = Vec::new();
    ports
        .try_reserve_exact(account.count)
        .map_err(|_| JointFrictionTrajectoryError::Allocation)?;
    for p in &account.ports[..account.count] {
        checkpoint(&mut cancelled)?;
        let [travel, distance, friction, drive, drive_work, friction_work] = p.sums.map(Dd::to_f64);
        if p.sums.iter().any(|x| !x.to_f64().is_finite()) || distance < 0.0 || friction_work > 0.0 {
            return Err(JointFrictionTrajectoryError::Unresolved(
                "finite physical port accounts",
            ));
        }
        let motion = if p.length > 0.0 {
            JointFrictionTrajectoryMotion::Revolute {
                angular_length_m: p.length,
                rotation_rad: travel,
                angular_distance_rad: distance,
                friction_impulse_nms: friction,
                drive_impulse_nms: drive,
            }
        } else {
            JointFrictionTrajectoryMotion::Prismatic {
                displacement_m: travel,
                distance_m: distance,
                friction_impulse_ns: friction,
                drive_impulse_ns: drive,
            }
        };
        ports.push(JointFrictionTrajectoryPort {
            joint_index: p.joint,
            motion,
            drive_work_j: drive_work,
            friction_work_j: friction_work,
        });
    }
    checkpoint(&mut cancelled)?;
    Ok(JointFrictionTrajectoryReport {
        states_after: current,
        start_time_s,
        target_time_s: target,
        duration_s,
        ports,
        mechanical_energy_change_j: energy.to_f64(),
        wrench_work_j: work[0].to_f64(),
        drive_work_j: work[1].to_f64(),
        dissipated_j: work[2].to_f64(),
        constraint_work_j: work[3].to_f64(),
        absolute_constraint_work_j: absolute_constraint,
        work_energy_residual_j: residual,
        estimated_error: JointFrictionTrajectoryAccuracy::from_values(estimates),
        accepted_steps: accepted,
        split_trials: trials,
        rejected_intervals: rejected,
        planned_live_bytes: bytes as usize,
    })
}
