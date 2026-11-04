//! Prescribed radial and angular actuators in the native articulated graph.

use super::{
    AngularDrive, AngularDriveStep, ArticulatedJoint, Dd, DistanceDrive, DistanceGraphControls,
    DistanceGraphError, DistanceGraphStep, DistanceLink, DistanceStepError, FreeFlow,
    JointReactionImpulse, JointRows, MAX_DISTANCE_ROWS, MassProperties, PointPlaneJoint,
    PrismaticJoint, RevoluteJoint, RevoluteMotorJoint, RigidBodyIntegrator, RigidBodyState, Row,
    Rows, SphericalJoint, Vec3, Wrench, advance_between_with, angular, driven, limit,
    planned_bytes_for_motion, poll, positive, reserved, revolute,
};

/// Finite passive joints and prescribed radial/angular actuators in one graph.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DrivenJoint {
    /// An ideal internal actuator with a held quadratic range command. A
    /// default drive supplies a fixed positive-length rod in the same graph.
    Distance {
        /// Finite body endpoints and the initial positive range [m].
        link: DistanceLink,
        /// Initial range rate and held range acceleration.
        drive: DistanceDrive,
    },
    /// Three Cartesian coincidence conditions, with free relative attitude.
    Spherical(SphericalJoint),
    /// Coincidence and directed axis alignment, with free axial rotation.
    Revolute(RevoluteJoint),
    /// A moving material point on a finite freely recoiling plane carrier.
    PointPlane(PointPlaneJoint),
    /// A freely sliding rail coordinate with locked relative joint frames.
    Prismatic(PrismaticJoint),
    /// A finite hinge with quadratic relative angle and signed source work.
    /// Zero commanded spin holds the relative angle; it is still a six-row
    /// actuator, not the five-row passive [`Self::Revolute`] joint.
    RevoluteMotor {
        /// Coincident material points and the current commanded joint frames.
        joint: RevoluteMotorJoint,
        /// Initial relative spin and held relative angular acceleration.
        drive: AngularDrive,
    },
}
impl DrivenJoint {
    pub(super) fn row_count(self) -> usize {
        match self {
            Self::Distance { .. } | Self::PointPlane(_) => 1,
            Self::Spherical(_) => 3,
            Self::Revolute(_) | Self::Prismatic(_) => 5,
            Self::RevoluteMotor { .. } => 6,
        }
    }

    pub(super) fn write_rows(self, rows: &mut [Row]) {
        match self {
            Self::Distance { link, .. } => ArticulatedJoint::Distance(link).write_rows(rows),
            Self::Spherical(joint) => ArticulatedJoint::Spherical(joint).write_rows(rows),
            Self::Revolute(joint) => ArticulatedJoint::Revolute(joint).write_rows(rows),
            Self::PointPlane(joint) => ArticulatedJoint::PointPlane(joint).write_rows(rows),
            Self::Prismatic(joint) => ArticulatedJoint::Prismatic(joint).write_rows(rows),
            Self::RevoluteMotor { joint, .. } => joint.write_rows(rows),
        }
    }
}

/// Continuation, scalar-row location and work of one authored joint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DrivenJointStep {
    /// Actuator endpoints with their represented range/rate or current angular
    /// frame/rate, or the unchanged passive joint. Continue with this value.
    pub joint_after: DrivenJoint,
    /// Row span in the native step, including every angular lock row.
    pub rows: JointRows,
    /// Represented endpoint range minus the commanded polynomial [m].
    /// Zero for angular actuators and passive joints.
    pub length_representation_residual_m: f64,
    /// Represented final range rate minus the commanded rate [m/s].
    /// Zero for angular actuators and passive joints.
    pub rate_representation_residual_mps: f64,
    /// Actual angular command representation; absent for radial/passive joints.
    pub angular: Option<AngularDriveStep>,
    /// Signed actuator source work [J]. Passive joints contribute zero.
    pub actuator_work_j: f64,
    /// Ordered native force/couple work minus actuator source work [J].
    /// Passive joint work is retained here; this quantity is not heat.
    pub constraint_work_defect_j: f64,
}

/// Atomic driven articulated motion with the original native reaction events.
#[derive(Debug, Clone, PartialEq)]
pub struct DrivenJointGraphStep<FreeStep = crate::StepReceipt> {
    /// One simultaneous native step, in scalar-row order. Point and pure-couple
    /// events retain their actual finite-body momentum and work residuals.
    pub step: DistanceGraphStep<JointReactionImpulse, FreeStep>,
    /// One bounded vector in authored joint order; no separate row-map copy.
    pub joints: Vec<DrivenJointStep>,
    /// Sum of signed prescribed actuator work, accumulated before rounding [J].
    pub actuator_work_j: f64,
    /// Total ordered native reaction work minus prescribed source work [J].
    /// Includes passive joint work and is never projected out of the bodies.
    pub constraint_work_defect_j: f64,
    /// Independently paired total reaction work [J], retaining small relative
    /// work under a common translational boost.
    pub constraint_work_j: f64,
    /// Requested world angular impulse of all reactions [N m s]. Pure couples
    /// cancel internally; finite point-closure errors remain in this account.
    pub constraint_angular_impulse_world_nms: Vec3,
    /// Logical owned payload, including all continuation receipts [bytes].
    pub planned_live_bytes: usize,
}

/// Advance quadratic distance/angular actuators and finite passive joints.
///
/// Every row enters the same full-inertia force/couple mobility and bounded
/// nonlinear solve. Actuator forces can therefore drive attached mechanisms
/// and transmit hinge couples while all bodies retain finite recoil. There is
/// no sequential solve of separately constrained subsystems. Native free steps,
/// local-branch admission, trial reuse, endpoint impulses and cancellation are
/// supplied by the original graph owner.
///
/// Each distance follows `length(t) = length(0) + rate*t + acceleration*t²/2`.
/// The existing distance command owner checks the entire interval's positive
/// range, endpoint representability and separate COM/spin range rates. Only
/// radial rows receive those commands; point/plane and passive angular rows
/// retain their original geometry. Passive joints supply no source work. Their
/// ordered point/couple work remains in the independently reported defect, so
/// it cannot be mislabeled as actuator energy or silently erased.
///
/// A revolute motor advances its current first joint frame about joint-local Z
/// by `rate*t + acceleration*t²/2`. The original full-frame angular Jacobian
/// subtracts this commanded motion before hidden-rate admission; each actual
/// pure-couple impulse is paired with that same world axis and command rate.
/// The complete interval excursion, including interior turning points, and
/// measured frame/rate representation errors must fit the native local bounds.
/// Continuation carries the advanced frame, avoiding absolute-angle wrapping.
///
/// The 32-scalar-row bound counts one row per rod or point/plane, three per ball
/// and five per passive hinge or slider, or six per revolute motor. Two fixed
/// stack banks and one continuation vector
/// extend the native buffers, with no per-trial allocation. Every refusal or
/// cancellation leaves all borrowed inputs unchanged. Command errors identify
/// the radial row or first motor angular row in native residual/event order.
///
/// This remains an ideal bilateral actuator with an external energy source and
/// smooth finite-body joints. Finite motor power/storage, efficiency, friction,
/// unilateral contact, global branches and exact energy remain separate models.
#[allow(clippy::too_many_arguments)]
pub fn advance_driven_joint_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[DrivenJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    cancelled: impl FnMut() -> bool,
) -> Result<DrivenJointGraphStep, DistanceGraphError> {
    advance_with(
        states, properties, joints, integrator, wrenches, duration_s, controls, 0, cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance_with<F: FreeFlow>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[DrivenJoint],
    integrator: RigidBodyIntegrator,
    flow: F,
    duration_s: f64,
    controls: DistanceGraphControls,
    extra_live_bytes: u128,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DrivenJointGraphStep<F::Receipt>, DistanceGraphError> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || !positive(controls.position_tolerance_m)
        || !positive(controls.velocity_tolerance_mps)
    {
        return Err(DistanceStepError::Invalid("prescribed joint graph controls").into());
    }
    limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    let maximum = controls.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut rows = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut targets = [driven::Target::default(); MAX_DISTANCE_ROWS];
    let mut row_count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let width = joint.row_count();
        limit("constraints", (row_count + width) as u128, maximum as u128)?;
        joint.write_rows(&mut rows[row_count..row_count + width]);
        match joint {
            DrivenJoint::Distance { link, drive } => {
                targets[row_count] =
                    driven::Target::prepare(row_count, link, drive, duration_s, controls)?;
            }
            DrivenJoint::RevoluteMotor { joint, drive } => {
                let target = angular::prepare(row_count + 3, joint, drive, duration_s, controls)?;
                targets[row_count + 3..row_count + 6].fill(target);
            }
            _ => {}
        }
        row_count += width;
    }
    let extra =
        joints.len() as u128 * core::mem::size_of::<DrivenJointStep>() as u128 + extra_live_bytes;
    let bytes =
        planned_bytes_for_motion::<JointReactionImpulse, F::Receipt>(states.len(), row_count)
            + extra;
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut native_controls = controls;
    native_controls.maximum_live_bytes -= extra as usize;
    let (step, work) = advance_between_with(
        states,
        properties,
        Rows::Mixed(&rows[..row_count]),
        integrator,
        flow,
        duration_s,
        native_controls,
        Some(&targets[..row_count]),
        &mut cancelled,
    )?;
    finish(
        step,
        joints,
        properties,
        controls,
        &targets[..row_count],
        work,
        bytes,
        &mut cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn finish<FreeStep>(
    step: DistanceGraphStep<JointReactionImpulse, FreeStep>,
    joints: &[DrivenJoint],
    properties: &[MassProperties],
    controls: DistanceGraphControls,
    targets: &[driven::Target],
    work: driven::Work,
    planned_live_bytes: u128,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<DrivenJointGraphStep<FreeStep>, DistanceGraphError> {
    let mut continuation = reserved(joints.len())?;
    let mut first_row = 0;
    let mut source = Dd::ZERO;
    let mut defect = Dd::ZERO;
    for &joint in joints {
        poll(cancelled)?;
        let row_count = joint.row_count();
        let mut joint_source = Dd::ZERO;
        let mut joint_defect = Dd::ZERO;
        for account in &work.links[first_row..first_row + row_count] {
            joint_source = joint_source + account.source;
            joint_defect = joint_defect + account.defect;
        }
        source = source + joint_source;
        defect = defect + joint_defect;
        if !driven::finite(joint_source)
            || !driven::finite(joint_defect)
            || !driven::finite(source)
            || !driven::finite(defect)
        {
            return Err(DistanceStepError::Invalid("prescribed joint graph work sum").into());
        }
        let mut result = DrivenJointStep {
            joint_after: joint,
            rows: JointRows {
                first_row,
                row_count,
            },
            length_representation_residual_m: 0.0,
            rate_representation_residual_mps: 0.0,
            angular: None,
            actuator_work_j: joint_source.to_f64(),
            constraint_work_defect_j: joint_defect.to_f64(),
        };
        if let DrivenJoint::Distance { link, drive } = joint {
            let target = targets[first_row];
            result.joint_after = DrivenJoint::Distance {
                link: DistanceLink {
                    length_m: target.length,
                    ..link
                },
                drive: DistanceDrive {
                    range_rate_mps: target.final_rate.to_f64(),
                    ..drive
                },
            };
            result.length_representation_residual_m = target.length_residual.to_f64();
            result.rate_representation_residual_mps =
                (Dd::from_f64(target.final_rate.to_f64()) - target.final_rate).to_f64();
        }
        if let DrivenJoint::RevoluteMotor { joint, drive } = joint {
            let (joint, drive, report) = angular::continuation(
                joint,
                drive,
                targets[first_row + 3],
                joint.bodies().map(|b| step.states_after[b]),
                joint.bodies().map(|b| properties[b]),
                controls,
            )?;
            result.joint_after = DrivenJoint::RevoluteMotor { joint, drive };
            result.angular = Some(report);
        }
        continuation.push(result);
        first_row += row_count;
    }
    let (constraint_work_j, constraint_angular_impulse_world_nms) =
        revolute::reaction_account(&step, cancelled)?;
    poll(cancelled)?;
    Ok(DrivenJointGraphStep {
        step,
        joints: continuation,
        actuator_work_j: source.to_f64(),
        constraint_work_defect_j: defect.to_f64(),
        constraint_work_j,
        constraint_angular_impulse_world_nms,
        planned_live_bytes: planned_live_bytes as usize,
    })
}
