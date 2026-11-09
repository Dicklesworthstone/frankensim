//! Private held-axis impulse phase shared by the two bonded moment laws.
//!
//! Public rolling and torsion owners select their physical axis and sealed
//! original source. This owner shares the actual native inertia, source-checked
//! midpoint solve, turn admission and independent momentum/work checks.

use super::super::hertz::{RigidHertzError, points, relative_spin};
use super::super::{
    ResultantImpulse, RigidContactBinding, RigidFrictionError, add_components, magnitude, poll,
    rounded_components,
};
use crate::{
    DynamicsError, MassProperties, MomentumImpulseReceipt, PointKinematics, Pose, RigidBodyState,
    Vec3, Wrench,
};
use fs_contact::normal_patch::adhesive::{rolling::*, torsion::*};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Explicit physical-step limits and independent dimensional acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveMomentControls {
    pub maximum_bodies: usize,
    pub maximum_step_s: f64,
    pub maximum_angle_increment_rad: f64,
    /// At most 128 bisections, two bracket evaluations and one elastic predictor.
    pub maximum_iterations: usize,
    pub angle_tolerance_rad: f64,
    pub geometry_tolerance_m: f64,
    pub angular_velocity_tolerance_rad_s: f64,
    pub linear_momentum_tolerance_ns: f64,
    pub angular_momentum_tolerance_nms: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct HeldGeometry {
    indices: [usize; 2],
    arms: [Vec3; 2],
    poses: [Pose; 2],
    properties: [MassProperties; 2],
    frame: TangentFrame,
}
pub(super) struct NativePair {
    pub geometry: HeldGeometry,
    before: [RigidBodyState; 2],
    initial: [PointKinematics; 2],
    axis: Vec3,
    mobility: f64,
}

#[derive(Clone, Copy)]
pub(super) struct MomentWork {
    mean: f64,
    storage_before: f64,
    storage_change: f64,
    loss: f64,
    applied_work: f64,
    plastic: bool,
}

// Only the two original source owners implement this private interface. It
// exports diagnostics, never a surrogate source or writable accepted history.
pub(super) trait MomentSource: Clone + core::fmt::Debug + PartialEq {
    type State: Copy;
    type Step: Copy + core::fmt::Debug + PartialEq;
    type Error: core::fmt::Debug;
    fn check(&self, state: &Self::State) -> Result<(), Self::Error>;
    fn advance(&self, state: &Self::State, angle: f64) -> Result<Self::Step, Self::Error>;
    fn critical_angle(&self) -> f64;
    fn critical_moment(&self) -> f64;
    fn elastic_fraction(&self, state: &Self::State) -> f64;
    fn successor(step: &Self::Step) -> Self::State;
    fn work(step: &Self::Step) -> MomentWork;
}
impl MomentSource for SourceJkrRollingContact {
    type State = SourceJkrRollingState;
    type Step = SourceJkrRollingStep;
    type Error = SourceJkrRollingError;
    fn check(&self, state: &Self::State) -> Result<(), Self::Error> {
        self.stored_energy_j(state).map(|_| ())
    }
    fn advance(&self, state: &Self::State, angle: f64) -> Result<Self::Step, Self::Error> {
        self.advance(state, angle)
    }
    fn critical_angle(&self) -> f64 {
        self.law().critical_rolling_angle_rad()
    }
    fn critical_moment(&self) -> f64 {
        self.law().critical_moment_n_m()
    }
    fn elastic_fraction(&self, state: &Self::State) -> f64 {
        state.elastic_shift_m() / self.law().parameters().critical_rolling_displacement_m
    }
    fn successor(step: &Self::Step) -> Self::State {
        step.state()
    }
    fn work(step: &Self::Step) -> MomentWork {
        let r = step.rolling();
        MomentWork {
            mean: r.mean_moment_n_m,
            storage_before: r.storage_before_j,
            storage_change: r.storage_change_j,
            loss: r.irreversible_loss_j,
            applied_work: r.applied_work_j,
            plastic: r.plastic_rolling_increment_m != 0.0,
        }
    }
}
impl MomentSource for SourceJkrTorsionContact {
    type State = SourceJkrTorsionState;
    type Step = SourceJkrTorsionStep;
    type Error = SourceJkrTorsionError;
    fn check(&self, state: &Self::State) -> Result<(), Self::Error> {
        self.stored_energy_j(state).map(|_| ())
    }
    fn advance(&self, state: &Self::State, angle: f64) -> Result<Self::Step, Self::Error> {
        self.advance(state, angle)
    }
    fn critical_angle(&self) -> f64 {
        self.law().parameters().critical_twist_angle_rad
    }
    fn critical_moment(&self) -> f64 {
        self.law().critical_moment_n_m()
    }
    fn elastic_fraction(&self, state: &Self::State) -> f64 {
        state.elastic_twist_rad() / self.critical_angle()
    }
    fn successor(step: &Self::Step) -> Self::State {
        step.state()
    }
    fn work(step: &Self::Step) -> MomentWork {
        let r = step.torsion();
        MomentWork {
            mean: r.mean_moment_n_m,
            storage_before: r.storage_before_j,
            storage_change: r.storage_change_j,
            loss: r.irreversible_loss_j,
            applied_work: r.applied_work_j,
            plastic: r.plastic_twist_increment_rad != 0.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct MomentReport<S: MomentSource> {
    pub contact: S::Step,
    pub angle_increment_rad: f64,
    pub rate_before_rad_s: f64,
    pub rate_after_rad_s: f64,
    pub rotational_mobility_per_kg_m2: f64,
    pub angular_impulse_on_a_world_nms: Vec3,
    pub external_wrenches: [Wrench; 2],
    pub body_displacement_world_m: [Vec3; 2],
    pub angular_displacement_body_rad: [Vec3; 2],
    pub kinetic_energy_change_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub external_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub angle_residual_rad: f64,
    pub angular_velocity_residual_rad_s: f64,
    pub body_angular_impulse_residual_nms: [Vec3; 2],
    pub body_linear_impulse_residual_ns: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub linear_momentum_residual_ns: Vec3,
    pub angular_momentum_residual_nms: Vec3,
    pub evaluations: usize,
}
pub(super) struct MomentCandidate<S: MomentSource> {
    pub bodies_after: [RigidBodyState; 2],
    pub source_after: S::State,
    pub report: MomentReport<S>,
}

#[derive(Debug)]
pub(super) enum MomentError<S: MomentSource> {
    Invalid(&'static str),
    Source(S::Error),
    Mechanics(RigidHertzError),
    GeometryMismatch,
    StepLimit,
    Unresolved {
        evaluations: usize,
        angle_residual_rad: f64,
    },
    Reversal,
    Audit(Box<MomentReport<S>>),
}
impl<S: MomentSource> From<RigidHertzError> for MomentError<S> {
    fn from(error: RigidHertzError) -> Self {
        Self::Mechanics(error)
    }
}
impl<S: MomentSource> From<DynamicsError> for MomentError<S> {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error.into())
    }
}
impl<S: MomentSource> From<RigidFrictionError> for MomentError<S> {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanics(error.into())
    }
}

/// Original native free-body phase when a bonded moment junction is absent.
/// Both lifecycle owners retain the same actual momentum/work and turn gates.
pub(super) struct DetachedPhase {
    pub impulses: [MomentumImpulseReceipt; 2],
    pub angle: f64,
    pub kinetic: f64,
    pub work: f64,
    pub defect: f64,
    pub linear_residual: Vec3,
    pub angular_residual: Vec3,
}

/// Called only after the original pair/geometry, controls and duration admit.
/// The selected source type carries errors; no constitutive source acts while
/// detached, and the measured angular travel creates no spring memory.
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_detached<S: MomentSource>(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    axis: Vec3,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveMomentControls,
    gate: &CancelGate,
) -> Result<DetachedPhase, MomentError<S>> {
    let apply = |i: usize| -> Result<_, MomentError<S>> {
        poll(gate)?;
        wrenches[i].validate()?;
        Ok(before[i].apply_momentum_impulse(
            properties[i],
            wrenches[i].force_world.scale(dt_s),
            wrenches[i].torque_body.scale(dt_s),
        )?)
    };
    // Each original native receipt independently retains actual momentum and
    // kinetic change, commanded midpoint work and their compensated defect.
    let impulses = [apply(0)?, apply(1)?];
    let mut kinetic = Dd::ZERO;
    let mut work = Dd::ZERO;
    let mut defect = Dd::ZERO;
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    for receipt in impulses {
        if magnitude(receipt.linear_residual_world) > controls.linear_momentum_tolerance_ns
            || magnitude(receipt.angular_residual_body) > controls.angular_momentum_tolerance_nms
        {
            return Err(MomentError::Invalid(
                "detached applied impulse exceeds native momentum tolerance",
            ));
        }
        kinetic = kinetic + Dd::from_f64(receipt.work.kinetic_energy_change_j);
        work = work + Dd::from_f64(receipt.work.impulse_work_j);
        defect = defect + Dd::from_f64(receipt.work.work_energy_residual_j.abs());
        add_components(&mut linear, receipt.linear_residual_world);
        let pose = receipt.state_before.pose();
        add_components(
            &mut angular,
            pose.position_world().cross(receipt.linear_residual_world),
        );
        add_components(
            &mut angular,
            pose.orientation()
                .rotate_body_to_world(receipt.angular_residual_body),
        );
    }
    let linear_residual = rounded_components(linear);
    let angular_residual = rounded_components(angular);
    // Keep independently recovered per-body work defects as well as the summed
    // reported work balance; neither opposite body defects nor large totals can
    // hide a rejected impulse. The original owner computes all kinetic work.
    defect = defect + (kinetic - work).abs();
    if [kinetic.to_f64(), work.to_f64(), defect.to_f64()]
        .iter()
        .any(|v| !v.is_finite())
        || defect.to_f64() > controls.energy_tolerance_j
        || !linear_residual.is_finite()
        || !angular_residual.is_finite()
        || magnitude(linear_residual) > controls.linear_momentum_tolerance_ns
        || magnitude(angular_residual) > controls.angular_momentum_tolerance_nms
    {
        return Err(MomentError::Invalid(
            "detached applied-load work or momentum balance exceeds tolerance",
        ));
    }
    let initial = points(before, properties, arms)?;
    let final_points = points(impulses.map(|r| r.state_after), properties, arms)?;
    let rate_before = relative_spin(initial).dot(axis);
    let rate_after = relative_spin(final_points).dot(axis);
    let angle = ((Dd::from_f64(rate_before) + Dd::from_f64(rate_after)) * Dd::from_f64(0.5 * dt_s))
        .to_f64();
    if !angle.is_finite() || angle.abs() > controls.maximum_angle_increment_rad {
        return Err(MomentError::StepLimit);
    }
    if (rate_before > 0.0 && rate_after < 0.0) || (rate_before < 0.0 && rate_after > 0.0) {
        // Constant held acceleration can turn between equal-angle endpoints.
        // Admit the actual quadratic turning excursion, not just net travel.
        let scale = rate_before.abs().max(rate_after.abs());
        let first = rate_before.abs() / scale;
        let fraction = first / (first + rate_after.abs() / scale);
        let excursion =
            (Dd::from_f64(rate_before) * Dd::from_f64(0.5 * dt_s) * Dd::from_f64(fraction))
                .to_f64();
        if !excursion.is_finite() || excursion.abs() > controls.maximum_angle_increment_rad {
            return Err(MomentError::StepLimit);
        }
    }
    poll(gate)?;
    Ok(DetachedPhase {
        impulses,
        angle,
        kinetic: kinetic.to_f64(),
        work: work.to_f64(),
        defect: defect.to_f64(),
        linear_residual,
        angular_residual,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn capture<S: MomentSource>(
    binding: RigidContactBinding,
    frame: TangentFrame,
    axis: Vec3,
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    controls: RigidAdhesiveMomentControls,
    gate: &CancelGate,
) -> Result<NativePair, MomentError<S>> {
    poll(gate)?;
    validate::<S>(controls)?;
    if bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
    {
        return Err(MomentError::Invalid(
            "distinct bounded finite body pair required",
        ));
    }
    let before = [bodies[binding.body_a], bodies[binding.body_b]];
    let properties = [properties[binding.body_a], properties[binding.body_b]];
    let arms = [binding.arm_a_body, binding.arm_b_body];
    let initial = points(before, properties, arms)?;
    if magnitude(initial[0].point_world.sub(initial[1].point_world)) > controls.geometry_tolerance_m
    {
        return Err(MomentError::GeometryMismatch);
    }
    let mut mobility = 0.0;
    for i in 0..2 {
        let body_axis = before[i].pose().orientation().rotate_world_to_body(axis);
        let factor = properties[i].angular_mobility_factor(body_axis);
        mobility += factor.dot(factor);
    }
    if !mobility.is_finite() || mobility <= 0.0 {
        return Err(MomentError::Invalid("unrepresentable rotational mobility"));
    }
    Ok(NativePair {
        geometry: HeldGeometry {
            indices: [binding.body_a, binding.body_b],
            arms,
            poses: before.map(|b| b.pose()),
            properties,
            frame,
        },
        before,
        initial,
        axis,
        mobility,
    })
}

fn validate<S: MomentSource>(c: RigidAdhesiveMomentControls) -> Result<(), MomentError<S>> {
    if c.maximum_bodies < 2
        || c.maximum_iterations == 0
        || c.maximum_iterations > 128
        || [
            c.maximum_step_s,
            c.maximum_angle_increment_rad,
            c.angle_tolerance_rad,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || [
            c.geometry_tolerance_m,
            c.angular_velocity_tolerance_rad_s,
            c.linear_momentum_tolerance_ns,
            c.angular_momentum_tolerance_nms,
            c.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(MomentError::Invalid(
            "finite positive bonded moment budgets required",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn prepare<S: MomentSource>(
    native: NativePair,
    geometry: HeldGeometry,
    source: &S,
    state: &S::State,
    wrenches: [Wrench; 2],
    dt_s: f64,
    controls: RigidAdhesiveMomentControls,
    gate: &CancelGate,
) -> Result<MomentCandidate<S>, MomentError<S>> {
    if native.geometry != geometry {
        return Err(MomentError::GeometryMismatch);
    }
    if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > controls.maximum_step_s {
        return Err(MomentError::StepLimit);
    }
    source.check(state).map_err(MomentError::Source)?;
    let rate_before = finite(relative_spin(native.initial).dot(native.axis))?;
    let mut accelerations = [Vec3::ZERO; 2];
    for i in 0..2 {
        wrenches[i].validate()?;
        accelerations[i] = native.before[i].pose().orientation().rotate_body_to_world(
            native.geometry.properties[i].angular_velocity_body_checked(wrenches[i].torque_body)?,
        );
    }
    let free_acceleration = finite(accelerations[0].sub(accelerations[1]).dot(native.axis))?;
    let free = finite(
        (Dd::from_f64(dt_s)
            * (Dd::from_f64(rate_before)
                + Dd::from_f64(0.5 * dt_s) * Dd::from_f64(free_acceleration)))
        .to_f64(),
    )?;
    let compliance = finite((0.5 * dt_s) * (dt_s * native.mobility))?;
    if compliance <= 0.0 {
        return Err(MomentError::Invalid("unrepresentable angular compliance"));
    }
    let (angle, step, evaluations) = solve(source, state, free, compliance, controls, gate)?;
    let work = S::work(&step);
    let impulse = finite(dt_s * work.mean)?;
    let angular_impulse = native.axis.scale(impulse);
    let mut after = native.before;
    let mut kinetic_change = Dd::ZERO;
    let mut linear_change = [Dd::ZERO; 3];
    let mut angular_change = [Dd::ZERO; 3];
    let mut external_linear = [Dd::ZERO; 3];
    let mut external_angular = [Dd::ZERO; 3];
    let mut body_residual = [Vec3::ZERO; 2];
    let mut body_linear_residual = [Vec3::ZERO; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        poll(gate)?;
        let body = native.before[i];
        let orientation = body.pose().orientation();
        let external_impulse = wrenches[i].force_world.scale(dt_s);
        let external_torque = wrenches[i].torque_body.scale(dt_s);
        let mut resultant = ResultantImpulse::default();
        resultant.add(
            body,
            native.geometry.arms[i],
            Vec3::ZERO,
            angular_impulse.scale(sign),
        );
        add_components(&mut resultant.linear, external_impulse);
        add_components(&mut resultant.angular, external_torque);
        // Apply the simultaneous contact and applied loads once per body, so
        // rounding an intermediate free state cannot change the contact law.
        let receipt = body.apply_momentum_impulse(
            native.geometry.properties[i],
            rounded_components(resultant.linear),
            rounded_components(resultant.angular),
        )?;
        after[i] = receipt.state_after;
        body_residual[i] = receipt.angular_residual_body;
        body_linear_residual[i] = receipt.linear_residual_world;
        kinetic_change = kinetic_change + Dd::from_f64(receipt.work.kinetic_energy_change_j);
        add_components(&mut linear_change, receipt.linear_increment_world);
        add_components(&mut external_linear, external_impulse);
        add_components(
            &mut angular_change,
            body.pose()
                .position_world()
                .cross(receipt.linear_increment_world),
        );
        add_components(
            &mut angular_change,
            orientation.rotate_body_to_world(receipt.angular_increment_body),
        );
        add_components(
            &mut external_angular,
            body.pose().position_world().cross(external_impulse),
        );
        add_components(
            &mut external_angular,
            orientation.rotate_body_to_world(external_torque),
        );
    }
    let final_points = points(after, native.geometry.properties, native.geometry.arms)?;
    let rate_after = finite(relative_spin(final_points).dot(native.axis))?;
    if ((rate_before > 0.0 && rate_after < 0.0) || (rate_before < 0.0 && rate_after > 0.0))
        && (wrenches.iter().any(|w| w.torque_body != Vec3::ZERO)
            || !admits_elastic_turn(&native, source, work, rate_before))
    {
        return Err(MomentError::Reversal);
    }
    let midpoint = (Dd::from_f64(rate_before) + Dd::from_f64(rate_after)) * Dd::from_f64(0.5);
    let native_angle = midpoint * Dd::from_f64(dt_s);
    let impulse_work = finite((midpoint * Dd::from_f64(impulse)).to_f64())?;
    let mut body_travel = [Vec3::ZERO; 2];
    let mut angular_travel = [Vec3::ZERO; 2];
    let mut external_work = Dd::ZERO;
    let travel = |a: Vec3, b: Vec3| {
        let component =
            |x, y| ((Dd::from_f64(x) + Dd::from_f64(y)) * Dd::from_f64(0.5 * dt_s)).to_f64();
        Vec3::new(
            component(a.x, b.x),
            component(a.y, b.y),
            component(a.z, b.z),
        )
    };
    let pairing = |a: Vec3, b: Vec3| {
        Dd::from_f64(a.x) * Dd::from_f64(b.x)
            + Dd::from_f64(a.y) * Dd::from_f64(b.y)
            + Dd::from_f64(a.z) * Dd::from_f64(b.z)
    };
    for i in 0..2 {
        body_travel[i] = travel(
            native.initial[i].center_of_mass_velocity_world,
            final_points[i].center_of_mass_velocity_world,
        );
        angular_travel[i] = travel(
            native.initial[i].angular_velocity_body,
            final_points[i].angular_velocity_body,
        );
        if !body_travel[i].is_finite() || !angular_travel[i].is_finite() {
            return Err(MomentError::Invalid("unrepresentable native travel"));
        }
        external_work = external_work
            + pairing(wrenches[i].force_world, body_travel[i])
            + pairing(wrenches[i].torque_body, angular_travel[i]);
    }
    let external_work_j = finite(external_work.to_f64())?;
    let linear_residual = rounded_components(core::array::from_fn(|d| {
        linear_change[d] - external_linear[d]
    }));
    let angular_residual = rounded_components(core::array::from_fn(|d| {
        angular_change[d] - external_angular[d]
    }));
    let report = MomentReport {
        contact: step,
        angle_increment_rad: angle,
        rate_before_rad_s: rate_before,
        rate_after_rad_s: rate_after,
        rotational_mobility_per_kg_m2: native.mobility,
        angular_impulse_on_a_world_nms: angular_impulse,
        external_wrenches: wrenches,
        body_displacement_world_m: body_travel,
        angular_displacement_body_rad: angular_travel,
        kinetic_energy_change_j: finite(kinetic_change.to_f64())?,
        midpoint_impulse_work_j: impulse_work,
        external_work_j,
        work_energy_residual_j: finite(
            (kinetic_change - Dd::from_f64(impulse_work) - external_work).to_f64(),
        )?,
        total_energy_residual_j: finite(
            (kinetic_change + Dd::from_f64(work.storage_change) + Dd::from_f64(work.loss)
                - external_work)
                .to_f64(),
        )?,
        angle_residual_rad: finite((native_angle - Dd::from_f64(angle)).to_f64())?,
        angular_velocity_residual_rad_s: finite(
            (Dd::from_f64(rate_after)
                - Dd::from_f64(rate_before)
                - Dd::from_f64(dt_s) * Dd::from_f64(free_acceleration)
                - Dd::from_f64(native.mobility) * Dd::from_f64(impulse))
            .to_f64(),
        )?,
        body_angular_impulse_residual_nms: body_residual,
        body_linear_impulse_residual_ns: body_linear_residual,
        linear_momentum_change_ns: rounded_components(linear_change),
        angular_momentum_change_nms: rounded_components(angular_change),
        linear_momentum_residual_ns: linear_residual,
        angular_momentum_residual_nms: angular_residual,
        evaluations,
    };
    if report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
        || (impulse_work + work.applied_work).abs() > controls.energy_tolerance_j
        || report.angle_residual_rad.abs() > controls.angle_tolerance_rad
        || report.angular_velocity_residual_rad_s.abs() > controls.angular_velocity_tolerance_rad_s
        || body_residual
            .iter()
            .any(|x| !x.is_finite() || magnitude(*x) > controls.angular_momentum_tolerance_nms)
        || body_linear_residual
            .iter()
            .any(|x| !x.is_finite() || magnitude(*x) > controls.linear_momentum_tolerance_ns)
        || !linear_residual.is_finite()
        || !angular_residual.is_finite()
        || magnitude(linear_residual) > controls.linear_momentum_tolerance_ns
        || magnitude(angular_residual) > controls.angular_momentum_tolerance_nms
    {
        return Err(MomentError::Audit(Box::new(report)));
    }
    poll(gate)?;
    Ok(MomentCandidate {
        bodies_after: after,
        source_after: S::successor(&step),
        report,
    })
}

fn solve<S: MomentSource>(
    source: &S,
    state: &S::State,
    free: f64,
    compliance: f64,
    controls: RigidAdhesiveMomentControls,
    gate: &CancelGate,
) -> Result<(f64, S::Step, usize), MomentError<S>> {
    let span = finite(compliance * source.critical_moment())?;
    let mut low = finite(free - span)?
        .next_down()
        .max(-controls.maximum_angle_increment_rad);
    let mut high = finite(free + span)?
        .next_up()
        .min(controls.maximum_angle_increment_rad);
    if low > high {
        return Err(MomentError::StepLimit);
    }
    let probe = |angle| -> Result<_, MomentError<S>> {
        poll(gate)?;
        let step = source.advance(state, angle).map_err(MomentError::Source)?;
        let residual = finite(
            (Dd::from_f64(angle)
                - Dd::from_f64(free)
                - Dd::from_f64(compliance) * Dd::from_f64(S::work(&step).mean))
            .to_f64(),
        )?;
        Ok((step, residual))
    };
    // Verify the analytic elastic predictor with the original source, including
    // retained low history bits and yield. A failed optional guess only falls back.
    let mut evaluations = 0;
    let denominator = 1.0 + (0.5 * span) / source.critical_angle();
    let predictor = (-span).mul_add(source.elastic_fraction(state), free) / denominator;
    if denominator.is_finite() && predictor.is_finite() && predictor >= low && predictor <= high {
        evaluations += 1;
        match probe(predictor) {
            Ok((step, residual))
                if !S::work(&step).plastic
                    && S::work(&step).loss == 0.0
                    && residual.abs() <= controls.angle_tolerance_rad =>
            {
                return Ok((predictor, step, evaluations));
            }
            _ => poll(gate)?,
        }
    }
    evaluations += 1;
    let (low_step, f_low) = probe(low)?;
    if f_low.abs() <= controls.angle_tolerance_rad {
        return Ok((low, low_step, evaluations));
    }
    evaluations += 1;
    let (high_step, f_high) = probe(high)?;
    if f_high.abs() <= controls.angle_tolerance_rad {
        return Ok((high, high_step, evaluations));
    }
    if f_low > 0.0 || f_high < 0.0 {
        return Err(MomentError::StepLimit);
    }
    let mut residual = f_low;
    for _ in 0..controls.maximum_iterations {
        let mid = 0.5 * low + 0.5 * high;
        evaluations += 1;
        let (step, value) = probe(mid)?;
        residual = value;
        if residual.abs() <= controls.angle_tolerance_rad {
            return Ok((mid, step, evaluations));
        }
        if mid == low || mid == high {
            break;
        }
        if residual < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    Err(MomentError::Unresolved {
        evaluations,
        angle_residual_rad: residual,
    })
}

// At held geometry, scalar impulse j gives omega=omega0+W*j and
// K=K0+omega0*j+W*j^2/2. Thus initial U+omega0^2/(2W) bounds the entire
// conservative scalar orbit. This arithmetic screen is not an interval proof.
fn admits_elastic_turn<S: MomentSource>(
    native: &NativePair,
    source: &S,
    work: MomentWork,
    rate: f64,
) -> bool {
    if work.plastic
        || work.loss != 0.0
        || (work.storage_before != 0.0 && !work.storage_before.is_normal())
    {
        return false;
    }
    // Absolute normwise allowances include fixed quaternion rotation errors,
    // admitted norm drift and inverse-inertia conditioning, even at cancellation.
    // Caller-authored mechanics tolerances cannot relax this physical screen.
    const GUARD: f64 = 4096.0 * f64::EPSILON;
    let (mut inverse_scale, mut rate_scale) = (0.0, 0.0);
    for i in 0..2 {
        let inertia = native.geometry.properties[i].principal_inertia_body();
        let minimum = inertia.x.min(inertia.y).min(inertia.z);
        let inverse = minimum.recip().next_up();
        let momentum = native.before[i].angular_momentum_body();
        let norm = if momentum == Vec3::ZERO {
            0.0
        } else {
            ((momentum.x.abs() + momentum.y.abs()).next_up() + momentum.z.abs()).next_up()
        };
        let amplitude = if norm == 0.0 {
            0.0
        } else {
            (norm / minimum).next_up()
        };
        if !inverse.is_normal() || (norm != 0.0 && (!norm.is_normal() || !amplitude.is_normal())) {
            return false;
        }
        inverse_scale = (inverse_scale + inverse).next_up();
        if amplitude > 0.0 {
            rate_scale = (rate_scale + amplitude).next_up();
        }
    }
    let mobility_error = (GUARD * inverse_scale).next_up();
    let rate_error = (GUARD * rate_scale).next_up();
    let mobility_lower = (native.mobility - mobility_error).next_down();
    let rate_upper = (rate.abs() + rate_error).next_up();
    if !mobility_error.is_normal()
        || !rate_error.is_normal()
        || !mobility_lower.is_normal()
        || mobility_lower <= 0.0
        || !rate_upper.is_normal()
    {
        return false;
    }
    let moment = source.critical_moment();
    let angle = source.critical_angle();
    if !moment.is_normal() || !angle.is_normal() {
        return false;
    }
    let yield_lower = (Dd::from_f64(0.5) * Dd::from_f64(moment) * Dd::from_f64(angle))
        .to_f64()
        .next_down();
    if !yield_lower.is_normal() || yield_lower <= 0.0 {
        return false;
    }
    let reduced_rate = Dd::from_f64(rate_upper) / Dd::from_f64(mobility_lower.sqrt());
    let energy_scale = moment.sqrt() * angle.sqrt();
    if !reduced_rate.to_f64().is_normal() || !energy_scale.is_normal() {
        return false;
    }
    let speed = reduced_rate / Dd::from_f64(energy_scale);
    let fraction = Dd::from_f64(work.storage_before) / Dd::from_f64(yield_lower) + speed * speed;
    let upper = fraction.to_f64().next_up();
    // Covers ordinary coefficient rounding, including the original law's order
    // of multiplication when preparing yield storage from physical parameters.
    upper.is_finite() && upper >= 0.0 && upper < 1.0 - GUARD
}

fn finite<S: MomentSource>(value: f64) -> Result<f64, MomentError<S>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MomentError::Invalid(
            "unrepresentable native moment arithmetic",
        ))
    }
}
