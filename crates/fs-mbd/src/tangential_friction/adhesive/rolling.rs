//! Source-owned adhesive rolling on two finite rigid bodies at held geometry.
//!
//! The longitudinal tangent is one fixed bending axis. Its actual inverse-
//! inertia response couples the source's exact work-mean resisting moment to
//! a discrete-gradient midpoint impulse phase. Equal/opposite free couples
//! change the original momenta; compensated native kinetic change independently
//! checks recoverable rolling storage and irreversible rolling loss.
//!
//! This supplies no force reaction enforcing no-slip rolling. The normal source,
//! attached area, frame, material arms, poses and mass properties remain fixed.
//! Full inertia can change angular rates on other axes; those are not clamped.
//! Constitutive path integration is exact, while the inertial midpoint phase is
//! a time discretization. Every nonzero endpoint rate reversal refuses for
//! refinement: its turning point could yield even when the net step is elastic.
//! Pose evolution, normal-area evolution, history transport, arbitrary vector
//! rolling and heat allocation remain with their physical owners.

use super::super::hertz::{RigidHertzError, points, prepare_pair_wrench, relative_spin};
use super::super::{RigidContactBinding, RigidFrictionError, magnitude, poll, vector};
use crate::{DynamicsError, MassProperties, PointKinematics, Pose, RigidBodyState, Vec3};
use fs_contact::normal_patch::adhesive::rolling::{
    SourceJkrRollingContact, SourceJkrRollingError, SourceJkrRollingState, SourceJkrRollingStep,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_tribo::partial_slip::TangentFrame;

/// Physical pair and coincident reaction markers; longitudinal is the rolling axis.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveRollingContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceJkrRollingContact,
}

/// Explicit physical-step limits and independent dimensional acceptance gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidAdhesiveRollingControls {
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
struct HeldGeometry {
    indices: [usize; 2],
    arms: [Vec3; 2],
    poses: [Pose; 2],
    properties: [MassProperties; 2],
    frame: TangentFrame,
}

/// Sealed source memory bound to the original pair, frame and held geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidAdhesiveRollingState {
    source: SourceJkrRollingState,
    geometry: HeldGeometry,
}
impl RigidAdhesiveRollingState {
    pub const fn source_state(&self) -> SourceJkrRollingState {
        self.source
    }
    pub fn accumulated_angle_rad(&self) -> f64 {
        self.source.accumulated_angle_rad()
    }
    pub fn elastic_shift_m(&self) -> f64 {
        self.source.elastic_shift_m()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingReport {
    pub contact: SourceJkrRollingStep,
    pub angle_increment_rad: f64,
    pub rolling_rate_before_rad_s: f64,
    pub rolling_rate_after_rad_s: f64,
    /// Axis dot (I_A^-1 + I_B^-1) axis, with each actual world inertia [1/(kg m2)].
    pub rotational_mobility_per_kg_m2: f64,
    pub angular_impulse_on_a_world_nms: Vec3,
    /// Native actual-momentum change, without subtracting large kinetic totals.
    pub kinetic_energy_change_j: f64,
    pub midpoint_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub angle_residual_rad: f64,
    pub angular_velocity_residual_rad_s: f64,
    pub body_angular_impulse_residual_nms: [Vec3; 2],
    pub linear_momentum_change_ns: Vec3,
    pub angular_momentum_change_nms: Vec3,
    pub evaluations: usize,
}

/// Complete staged proposal; the A/B ordering is the original binding ordering.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveRollingCandidate {
    pub bodies_after: [RigidBodyState; 2],
    pub state_after: RigidAdhesiveRollingState,
    pub report: RigidAdhesiveRollingReport,
}

#[derive(Debug)]
pub enum RigidAdhesiveRollingError {
    Invalid(&'static str),
    Source(SourceJkrRollingError),
    Mechanics(RigidHertzError),
    GeometryMismatch,
    StepLimit,
    Unresolved {
        evaluations: usize,
        angle_residual_rad: f64,
    },
    /// A net-angle trial cannot certify the maximum shift of a reversing path.
    UnresolvedRollingReversal,
    Audit(Box<RigidAdhesiveRollingReport>),
}
impl From<SourceJkrRollingError> for RigidAdhesiveRollingError {
    fn from(x: SourceJkrRollingError) -> Self {
        Self::Source(x)
    }
}
impl From<RigidHertzError> for RigidAdhesiveRollingError {
    fn from(x: RigidHertzError) -> Self {
        Self::Mechanics(x)
    }
}
impl From<DynamicsError> for RigidAdhesiveRollingError {
    fn from(x: DynamicsError) -> Self {
        Self::Mechanics(x.into())
    }
}
impl From<RigidFrictionError> for RigidAdhesiveRollingError {
    fn from(x: RigidFrictionError) -> Self {
        Self::Mechanics(x.into())
    }
}
impl core::fmt::Display for RigidAdhesiveRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native adhesive rolling: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveRollingError {}

struct NativePair {
    geometry: HeldGeometry,
    before: [RigidBodyState; 2],
    initial: [PointKinematics; 2],
    axis: Vec3,
    mobility: f64,
}
impl RigidAdhesiveRollingContact<'_> {
    /// Start unstrained rolling on this exact pair and accepted normal source.
    pub fn initial_state(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        controls: RigidAdhesiveRollingControls,
        gate: &CancelGate,
    ) -> Result<RigidAdhesiveRollingState, RigidAdhesiveRollingError> {
        let native = self.capture(bodies, properties, controls, gate)?;
        Ok(RigidAdhesiveRollingState {
            source: self.source.initial_state(),
            geometry: native.geometry,
        })
    }
    fn capture(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        controls: RigidAdhesiveRollingControls,
        gate: &CancelGate,
    ) -> Result<NativePair, RigidAdhesiveRollingError> {
        poll(gate)?;
        controls.validate()?;
        let binding = self.binding;
        if bodies.len() != properties.len()
            || bodies.len() > controls.maximum_bodies
            || binding.body_a >= bodies.len()
            || binding.body_b >= bodies.len()
            || binding.body_a == binding.body_b
        {
            return Err(RigidAdhesiveRollingError::Invalid(
                "distinct bounded finite body pair required",
            ));
        }
        let before = [bodies[binding.body_a], bodies[binding.body_b]];
        let properties = [properties[binding.body_a], properties[binding.body_b]];
        let arms = [binding.arm_a_body, binding.arm_b_body];
        let initial = points(before, properties, arms)?;
        if magnitude(initial[0].point_world.sub(initial[1].point_world))
            > controls.geometry_tolerance_m
        {
            return Err(RigidAdhesiveRollingError::GeometryMismatch);
        }
        let axis = vector(self.frame.longitudinal());
        let mut mobility = 0.0;
        for i in 0..2 {
            let body_axis = before[i].pose().orientation().rotate_world_to_body(axis);
            let factor = properties[i].angular_mobility_factor(body_axis);
            mobility += factor.dot(factor);
        }
        if !mobility.is_finite() || mobility <= 0.0 {
            return Err(RigidAdhesiveRollingError::Invalid(
                "unrepresentable rotational mobility",
            ));
        }
        Ok(NativePair {
            geometry: HeldGeometry {
                indices: [binding.body_a, binding.body_b],
                arms,
                poses: before.map(|b| b.pose()),
                properties,
                frame: self.frame,
            },
            before,
            initial,
            axis,
            mobility,
        })
    }
}
impl RigidAdhesiveRollingControls {
    fn validate(self) -> Result<(), RigidAdhesiveRollingError> {
        if self.maximum_bodies < 2
            || self.maximum_iterations == 0
            || self.maximum_iterations > 128
            || [
                self.maximum_step_s,
                self.maximum_angle_increment_rad,
                self.angle_tolerance_rad,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || [
                self.geometry_tolerance_m,
                self.angular_velocity_tolerance_rad_s,
                self.linear_momentum_tolerance_ns,
                self.angular_momentum_tolerance_nms,
                self.energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(RigidAdhesiveRollingError::Invalid(
                "finite positive rolling budgets required",
            ));
        }
        Ok(())
    }
}

/// Prepare original-source rolling and actual native recoil without publication.
#[allow(clippy::too_many_lines)]
pub fn prepare_rigid_adhesive_rolling(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    state: &RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingCandidate, RigidAdhesiveRollingError> {
    let native = contact.capture(bodies, properties, controls, gate)?;
    if native.geometry != state.geometry {
        return Err(RigidAdhesiveRollingError::GeometryMismatch);
    }
    if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > controls.maximum_step_s {
        return Err(RigidAdhesiveRollingError::StepLimit);
    }
    contact.source.stored_energy_j(&state.source)?;
    let rate_before = finite(relative_spin(native.initial).dot(native.axis))?;
    let free = finite(dt_s * rate_before)?;
    let compliance = finite((0.5 * dt_s) * (dt_s * native.mobility))?;
    if compliance <= 0.0 {
        return Err(RigidAdhesiveRollingError::Invalid(
            "unrepresentable angular compliance",
        ));
    }
    let (angle, step, evaluations) = solve(
        contact.source,
        &state.source,
        free,
        compliance,
        controls,
        gate,
    )?;
    let impulse = finite(dt_s * step.rolling().mean_moment_n_m)?;
    let angular_impulse = native.axis.scale(impulse);
    let phase = prepare_pair_wrench(
        native.before,
        native.geometry.properties,
        native.geometry.arms,
        Vec3::ZERO,
        angular_impulse,
        gate,
    )?;
    let after = phase.bodies_after;
    let final_points = points(after, native.geometry.properties, native.geometry.arms)?;
    let rate_after = finite(relative_spin(final_points).dot(native.axis))?;
    if (rate_before > 0.0 && rate_after < 0.0) || (rate_before < 0.0 && rate_after > 0.0) {
        return Err(RigidAdhesiveRollingError::UnresolvedRollingReversal);
    }
    let midpoint = (Dd::from_f64(rate_before) + Dd::from_f64(rate_after)) * Dd::from_f64(0.5);
    let native_angle = midpoint * Dd::from_f64(dt_s);
    let impulse_work = finite((midpoint * Dd::from_f64(impulse)).to_f64())?;
    let mut body_residual = [Vec3::ZERO; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let expected = native.before[i]
            .pose()
            .orientation()
            .rotate_world_to_body(angular_impulse.scale(sign));
        let a = after[i].angular_momentum_body();
        let b = native.before[i].angular_momentum_body();
        let difference = |x, y, z| (Dd::from_f64(x) - Dd::from_f64(y) - Dd::from_f64(z)).to_f64();
        body_residual[i] = Vec3::new(
            difference(a.x, b.x, expected.x),
            difference(a.y, b.y, expected.y),
            difference(a.z, b.z, expected.z),
        );
    }
    let rolling = step.rolling();
    let report = RigidAdhesiveRollingReport {
        contact: step,
        angle_increment_rad: angle,
        rolling_rate_before_rad_s: rate_before,
        rolling_rate_after_rad_s: rate_after,
        rotational_mobility_per_kg_m2: native.mobility,
        angular_impulse_on_a_world_nms: angular_impulse,
        kinetic_energy_change_j: phase.kinetic_change_j,
        midpoint_impulse_work_j: impulse_work,
        work_energy_residual_j: finite(
            (Dd::from_f64(phase.kinetic_change_j) - Dd::from_f64(impulse_work)).to_f64(),
        )?,
        total_energy_residual_j: finite(
            (Dd::from_f64(phase.kinetic_change_j)
                + Dd::from_f64(rolling.storage_change_j)
                + Dd::from_f64(rolling.irreversible_loss_j))
            .to_f64(),
        )?,
        angle_residual_rad: finite((native_angle - Dd::from_f64(angle)).to_f64())?,
        angular_velocity_residual_rad_s: finite(
            (Dd::from_f64(rate_after)
                - Dd::from_f64(rate_before)
                - Dd::from_f64(native.mobility) * Dd::from_f64(impulse))
            .to_f64(),
        )?,
        body_angular_impulse_residual_nms: body_residual,
        linear_momentum_change_ns: phase.linear_change_ns,
        angular_momentum_change_nms: phase.angular_change_nms,
        evaluations,
    };
    if report.work_energy_residual_j.abs() > controls.energy_tolerance_j
        || report.total_energy_residual_j.abs() > controls.energy_tolerance_j
        || (impulse_work + rolling.applied_work_j).abs() > controls.energy_tolerance_j
        || report.angle_residual_rad.abs() > controls.angle_tolerance_rad
        || report.angular_velocity_residual_rad_s.abs() > controls.angular_velocity_tolerance_rad_s
        || body_residual
            .iter()
            .any(|x| !x.is_finite() || magnitude(*x) > controls.angular_momentum_tolerance_nms)
        || magnitude(phase.linear_change_ns) > controls.linear_momentum_tolerance_ns
        || magnitude(phase.angular_change_nms) > controls.angular_momentum_tolerance_nms
    {
        return Err(RigidAdhesiveRollingError::Audit(Box::new(report)));
    }
    poll(gate)?;
    Ok(RigidAdhesiveRollingCandidate {
        bodies_after: after,
        state_after: RigidAdhesiveRollingState {
            source: step.state(),
            geometry: state.geometry,
        },
        report,
    })
}

/// Publish both native momenta and original rolling memory only after all gates.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_rolling(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    state: &mut RigidAdhesiveRollingState,
    contact: RigidAdhesiveRollingContact<'_>,
    dt_s: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveRollingReport, RigidAdhesiveRollingError> {
    let candidate =
        prepare_rigid_adhesive_rolling(bodies, properties, state, contact, dt_s, controls, gate)?;
    poll(gate)?;
    bodies[contact.binding.body_a] = candidate.bodies_after[0];
    bodies[contact.binding.body_b] = candidate.bodies_after[1];
    *state = candidate.state_after;
    Ok(candidate.report)
}

fn solve(
    source: &SourceJkrRollingContact,
    state: &SourceJkrRollingState,
    free: f64,
    compliance: f64,
    controls: RigidAdhesiveRollingControls,
    gate: &CancelGate,
) -> Result<(f64, SourceJkrRollingStep, usize), RigidAdhesiveRollingError> {
    let span = finite(compliance * source.law().critical_moment_n_m())?;
    let mut low = finite(free - span)?
        .next_down()
        .max(-controls.maximum_angle_increment_rad);
    let mut high = finite(free + span)?
        .next_up()
        .min(controls.maximum_angle_increment_rad);
    if low > high {
        return Err(RigidAdhesiveRollingError::StepLimit);
    }
    let probe = |angle| -> Result<_, RigidAdhesiveRollingError> {
        poll(gate)?;
        let step = source.advance(state, angle)?;
        let residual = finite(
            (Dd::from_f64(angle)
                - Dd::from_f64(free)
                - Dd::from_f64(compliance) * Dd::from_f64(step.rolling().mean_moment_n_m))
            .to_f64(),
        )?;
        Ok((step, residual))
    };
    // Ordinary elastic midpoint motion has an analytic solution. Use the
    // original source to verify it, including low history bits and yield;
    // this predictor never supplies a replacement constitutive response.
    let mut evaluations = 0;
    let critical_angle = source.law().critical_rolling_angle_rad();
    let fraction =
        state.elastic_shift_m() / source.law().parameters().critical_rolling_displacement_m;
    let denominator = 1.0 + (0.5 * span) / critical_angle;
    let predictor = (-span).mul_add(fraction, free) / denominator;
    if denominator.is_finite() && predictor.is_finite() && predictor >= low && predictor <= high {
        evaluations += 1;
        match probe(predictor) {
            Ok((step, residual))
                if step.rolling().plastic_rolling_increment_m == 0.0
                    && step.rolling().irreversible_loss_j == 0.0
                    && residual.abs() <= controls.angle_tolerance_rad =>
            {
                return Ok((predictor, step, evaluations));
            }
            // An optional failed predictor cannot add a refusal to the
            // original bracketed solve. Cancellation still takes precedence.
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
        return Err(RigidAdhesiveRollingError::StepLimit);
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
    Err(RigidAdhesiveRollingError::Unresolved {
        evaluations,
        angle_residual_rad: residual,
    })
}
fn finite(value: f64) -> Result<f64, RigidAdhesiveRollingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RigidAdhesiveRollingError::Invalid(
            "unrepresentable native rolling arithmetic",
        ))
    }
}
