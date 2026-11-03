//! Spherical and distance closures on the original native graph solver.

use super::{
    ActionReactionImpulseReceipt, Dd, DistanceGraphControls, DistanceGraphError, DistanceGraphStep,
    DistanceLink, DistanceStepError, Geometry, MAX_DISTANCE_ROWS, MassProperties,
    RigidBodyIntegrator, RigidBodyState, Row, RowKind, Rows, Vec3, Wrench, advance_between, driven,
    limit, planned_bytes, poll, reserved,
};

/// A ball-and-socket joint: two material points coincide, while relative
/// attitude remains free. Both endpoints belong to finite native rigid bodies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalJoint {
    /// Distinct indices into the native body/property/wrench slices.
    pub bodies: [usize; 2],
    /// Fixed material points measured from the corresponding body COMs [m].
    pub points_body_m: [Vec3; 2],
}

/// Smooth bilateral closures that can share arbitrary finite rigid bodies.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GraphJoint {
    /// One positive radial distance, retaining the original rod geometry owner.
    Distance(DistanceLink),
    /// Three independent Cartesian coincidence conditions; not a zero-length rod.
    Spherical(SphericalJoint),
}

/// Location of one authored joint's events and residuals in the native step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JointRows {
    /// First index into `step.position_impulses`, `velocity_impulses` and `residuals`.
    pub first_row: usize,
    /// One for distance, three for spherical, five for revolute. Each point
    /// closure uses world x/y/z; a hinge then uses its two body-relative axes.
    pub row_count: usize,
}

/// Atomic motion of a mixed independent-body joint graph.
#[derive(Debug, Clone, PartialEq)]
pub struct JointGraphStep<E = ActionReactionImpulseReceipt> {
    /// Original solver receipt, in scalar-row order. Every original free step
    /// and every actual native reaction remain available for inspection.
    pub step: DistanceGraphStep<E>,
    /// Row spans in authored joint order; no persistent state is inferred.
    pub joint_rows: Vec<JointRows>,
    /// Signed reaction work [J], paired from native relative COM/spin components
    /// before rounding their point velocities. This retains physical small work
    /// under common boosts; original event diagnostics remain in `step`.
    pub constraint_work_j: f64,
    /// Requested angular impulse of all reaction events about the world origin.
    /// Coincident spherical endpoints have no internal couple. Accepted finite
    /// gap residuals can leave a small couple, which is reported rather than
    /// erased. Actual momentum increments remain in the native pair receipts.
    pub constraint_angular_impulse_world_nms: Vec3,
    /// Logical owned payload, including the authored-joint row map [bytes].
    pub planned_live_bytes: usize,
}

/// Advance mixed spherical and fixed-distance joints using one coupled solve.
///
/// A spherical joint contributes the three world Cartesian components of
/// `x_second(point)-x_first(point)` and their native relative velocities. The
/// existing full-inertia `J M^-1 J^T` includes every cross response between
/// axes and between shared bodies. All position corrections, native free
/// steps, backtracking and endpoint velocity reactions use the original graph
/// owner. A zero-length radial norm and its undefined direction are never used.
///
/// `maximum_constraints` counts scalar rows, so at most ten spherical joints
/// plus two rods fit the 32-row solve. The position/speed tolerances bound each
/// spherical joint's Euclidean closure norm; the impulse limit bounds its
/// complete reaction vector at either event. These bounds do not depend on
/// the authored world-axis orientation. Rows apply at unchanged native poses
/// in joint order, retaining shared-body work cross terms.
///
/// Both bodies remain finite and free. Relative joint rotation is unconstrained:
/// this is not a weld, contact/impact law or kinematic support. The local smooth
/// branch, independent-row admission and time-discretization limits are those
/// of the original graph solver. No exact energy, global-root or symplectic
/// claim follows from a small residual. Finite closure gaps can produce a
/// residual angular impulse; it remains in the returned account.
///
/// One fixed-size stack row bank and one bounded row-map vector extend the
/// original reused trial buffers. No per-trial allocation is introduced.
/// Cancellation and every late refusal leave all input states unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_joint_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[GraphJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<JointGraphStep, DistanceGraphError> {
    poll(&mut cancelled)?;
    let maximum = controls.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut rows = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut row_count = 0;
    for joint in joints {
        poll(&mut cancelled)?;
        let width = match joint {
            GraphJoint::Distance(_) => 1,
            GraphJoint::Spherical(_) => 3,
        };
        limit("constraints", (row_count + width) as u128, maximum as u128)?;
        match joint {
            GraphJoint::Distance(link) => rows[row_count] = Row::distance(*link),
            GraphJoint::Spherical(joint) => {
                for axis in 0..3 {
                    rows[row_count + axis] = Row {
                        bodies: joint.bodies,
                        points_body_m: joint.points_body_m,
                        kind: RowKind::Cartesian(axis),
                    };
                }
            }
        }
        row_count += width;
    }
    let extra = joints.len() as u128 * core::mem::size_of::<JointRows>() as u128;
    let bytes = planned_bytes(states.len(), row_count) + extra;
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut native_controls = controls;
    native_controls.maximum_live_bytes -= extra as usize;
    let (step, _) = advance_between::<ActionReactionImpulseReceipt>(
        states,
        properties,
        Rows::Mixed(&rows[..row_count]),
        integrator,
        wrenches,
        duration_s,
        native_controls,
        None,
        &mut cancelled,
    )?;
    let mut joint_rows = reserved(joints.len())?;
    let mut first_row = 0;
    for joint in joints {
        poll(&mut cancelled)?;
        let row_count = match joint {
            GraphJoint::Distance(_) => 1,
            GraphJoint::Spherical(_) => 3,
        };
        joint_rows.push(JointRows {
            first_row,
            row_count,
        });
        first_row += row_count;
    }
    let mut work = Dd::ZERO;
    let mut angular = [Dd::ZERO; 3];
    for event in step.position_impulses.iter().chain(&step.velocity_impulses) {
        poll(&mut cancelled)?;
        let before = driven::relative_velocity(
            event.first.point_kinematics_before,
            event.second.point_kinematics_before,
        );
        let after = driven::relative_velocity(
            event.first.point_kinematics_after,
            event.second.point_kinematics_after,
        );
        let impulse = event.second.impulse_world;
        for ((before, after), j) in before
            .into_iter()
            .zip(after)
            .zip([impulse.x, impulse.y, impulse.z])
        {
            work = work + (before + after) * Dd::from_f64(0.5) * Dd::from_f64(j);
        }
        let native = event.angular_impulse_world;
        for (sum, value) in angular.iter_mut().zip([native.x, native.y, native.z]) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    let constraint_work_j = work.to_f64();
    let constraint_angular_impulse_world_nms = Vec3::new(
        angular[0].to_f64(),
        angular[1].to_f64(),
        angular[2].to_f64(),
    );
    if !constraint_work_j.is_finite() || !constraint_angular_impulse_world_nms.is_finite() {
        return Err(DistanceStepError::Invalid("joint graph reaction account").into());
    }
    poll(&mut cancelled)?;
    Ok(JointGraphStep {
        step,
        joint_rows,
        constraint_work_j,
        constraint_angular_impulse_world_nms,
        planned_live_bytes: bytes as usize,
    })
}

pub(super) fn geometry(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    points: [Vec3; 2],
    axis: usize,
) -> Result<Geometry, DistanceStepError> {
    let a = states[0].point_kinematics(properties[0], points[0])?;
    let b = states[1].point_kinematics(properties[1], points[1])?;
    let component = |v: Vec3| [v.x, v.y, v.z][axis];
    let gap = (Dd::from_f64(component(states[1].pose().position_world()))
        - Dd::from_f64(component(states[0].pose().position_world()))
        + Dd::from_f64(component(b.arm_world))
        - Dd::from_f64(component(a.arm_world)))
    .to_f64();
    let rate = driven::relative_velocity(a, b)[axis].to_f64();
    if !gap.is_finite() || !rate.is_finite() {
        return Err(DistanceStepError::Invalid("Cartesian point closure"));
    }
    let normal = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ][axis];
    Ok(Geometry { normal, gap, rate })
}
