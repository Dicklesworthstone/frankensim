//! Finite-body hinges with native force and pure-couple reactions.

use super::{
    ActionReactionImpulseReceipt, Dd, DistanceGraphControls, DistanceGraphError, DistanceGraphStep,
    DistanceLink, DistanceStepError, Geometry, JointGraphStep, JointRows, MAX_DISTANCE_ROWS,
    MassProperties, NativeReaction, RigidBodyIntegrator, RigidBodyState, Row, RowKind, Rows,
    SphericalJoint, Vec3, Wrench, advance_between, driven, limit, planned_bytes_for, poll,
    positive, reserved,
};
use crate::MomentumImpulseReceipt;

/// A frictionless bilateral hinge between two finite native rigid bodies.
///
/// The material points coincide and the directed axes align, leaving relative
/// rotation about that common axis free. The positive-alignment branch is
/// explicit: oppositely directed axes are not silently reinterpreted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RevoluteJoint {
    bodies: [usize; 2],
    points_body_m: [Vec3; 2],
    axes_body: [Vec3; 2],
    transverse_second_body: [Vec3; 2],
    angular_length_m: f64,
}
impl RevoluteJoint {
    /// Prepare unit axes and one orthonormal transverse basis in body two.
    ///
    /// `angular_length_m = L` is an explicit dimensional scale for the two
    /// alignment equations. With graph tolerances `(p, v)`, admission requires
    /// `sin(axis_angle) <= p/L`, `axis_first.dot(axis_second) > 0`, and the
    /// Euclidean norm of the two alignment derivatives at most `v/L` [1/s].
    /// Thus a scale change also changes these authored tolerances. At an
    /// exactly aligned state, that derivative norm is the relative tilt speed.
    /// The complete reaction couple is bounded by `L * maximum_impulse_ns`.
    ///
    /// # Errors
    /// Refuses equal body indices, nonfinite points, zero/nonfinite axes, and
    /// a nonpositive/nonfinite length scale. The graph additionally requires
    /// `L > position_tolerance_m` and existing distinct body indices.
    pub fn new(
        bodies: [usize; 2],
        points_body_m: [Vec3; 2],
        axes_body: [Vec3; 2],
        angular_length_m: f64,
    ) -> Result<Self, DistanceStepError> {
        if bodies[0] == bodies[1]
            || points_body_m.iter().any(|point| !point.is_finite())
            || !positive(angular_length_m)
        {
            return Err(DistanceStepError::Invalid(
                "finite hinge geometry and length scale",
            ));
        }
        let axes_body = [
            axes_body[0].normalized("hinge first axis")?,
            axes_body[1].normalized("hinge second axis")?,
        ];
        let axis = axes_body[1];
        let seed = if axis.x.abs() <= axis.y.abs().min(axis.z.abs()) {
            Vec3::new(1.0, 0.0, 0.0)
        } else if axis.y.abs() <= axis.z.abs() {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let first = axis.cross(seed).normalized("hinge transverse basis")?;
        let second = axis.cross(first).normalized("hinge transverse basis")?;
        Ok(Self {
            bodies,
            points_body_m,
            axes_body,
            transverse_second_body: [first, second],
            angular_length_m,
        })
    }
    /// Body indices in the supplied graph slices.
    #[must_use]
    pub const fn bodies(self) -> [usize; 2] {
        self.bodies
    }
    /// Coincident material points, measured from their native COMs [m].
    #[must_use]
    pub const fn points_body_m(self) -> [Vec3; 2] {
        self.points_body_m
    }
    /// Prepared unit axes in the respective native body frames.
    #[must_use]
    pub const fn axes_body(self) -> [Vec3; 2] {
        self.axes_body
    }
    /// Authored scale converting alignment equations to metres.
    #[must_use]
    pub const fn angular_length_m(self) -> f64 {
        self.angular_length_m
    }
}

/// Smooth bilateral joints in one shared force-and-couple graph.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArticulatedJoint {
    /// One radial closure through the original rod owner.
    Distance(DistanceLink),
    /// Three Cartesian point closures, with free relative attitude.
    Spherical(SphericalJoint),
    /// Three point closures and two body-relative axis constraints.
    Revolute(RevoluteJoint),
}
impl ArticulatedJoint {
    pub(super) fn row_count(self) -> usize {
        match self {
            Self::Distance(_) => 1,
            Self::Spherical(_) => 3,
            Self::Revolute(_) => 5,
        }
    }

    // Callers bound the row count before supplying exactly this joint's span.
    // Keep one geometry owner for static, actuated and spatial-load graphs.
    pub(super) fn write_rows(self, rows: &mut [Row]) {
        let (bodies, points) = match self {
            Self::Distance(link) => {
                rows[0] = Row::distance(link);
                return;
            }
            Self::Spherical(joint) => (joint.bodies, joint.points_body_m),
            Self::Revolute(joint) => (joint.bodies, joint.points_body_m),
        };
        for axis in 0..3 {
            rows[axis] = Row {
                bodies,
                points_body_m: points,
                kind: RowKind::Cartesian(axis),
            };
        }
        if let Self::Revolute(joint) = self {
            for component in 0..2 {
                rows[3 + component] = Row {
                    bodies,
                    points_body_m: [Vec3::ZERO; 2],
                    kind: RowKind::Revolute(AngularConstraint {
                        axes_body: joint.axes_body,
                        transverse_second_body: joint.transverse_second_body[component],
                        length_m: joint.angular_length_m,
                        component,
                    }),
                };
            }
        }
    }
}

/// Native updates for one equal-and-opposite central angular impulse pair.
/// No artificial point, arm, or force is introduced for a pure couple.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoupleImpulseReceipt {
    /// Actual first-body momentum update, including representability defects.
    pub first: MomentumImpulseReceipt,
    /// Actual second-body momentum update, including representability defects.
    pub second: MomentumImpulseReceipt,
    /// World angular impulse on body two [N m s]; body one receives its negative.
    pub couple_on_second_world_nms: Vec3,
    /// Actual summed world angular momentum increment [N m s]. The commanded
    /// total is zero; this field retains finite transform/update error.
    pub angular_momentum_change_world_nms: Vec3,
    /// Paired actual kinetic-energy change from native momentum increments [J].
    pub kinetic_energy_change_j: f64,
    /// Commanded world couple dotted with midpoint relative angular velocity [J].
    pub impulse_work_j: f64,
    /// Actual paired kinetic change minus commanded couple work [J].
    pub work_energy_residual_j: f64,
}

/// A scalar row's actual native reaction, in authored row order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointReactionImpulse {
    /// Original equal-and-opposite material-point impulses, unchanged.
    Point(ActionReactionImpulseReceipt),
    /// Equal-and-opposite COM couples through the native momentum owner.
    Couple(CoupleImpulseReceipt),
}

/// Native graph receipt plus authored-joint row spans and paired work account.
pub type ArticulatedGraphStep<FreeStep = crate::StepReceipt> =
    JointGraphStep<JointReactionImpulse, FreeStep>;

/// Advance distance, ball and revolute joints in one coupled full-inertia solve.
///
/// Each hinge contributes five independent rows. Its angular equations are
/// `g_i = L * a_first.dot(t_second_i)`, with two body-fixed transverse directions
/// on the second body. Their exact world torque Jacobian on body two is
/// `L * (t_second_i cross a_first)`, and its negative acts on body one. Both axes
/// move with the native bodies. Cross responses between all force/couple rows
/// and shared bodies enter the original `J M^-1 J^T` factorization.
///
/// The original bounded impulse/free-step/impulse owner supplies all trials,
/// backtracking, native free evolution and endpoint reactions. Reactions apply
/// in row order at unchanged poses, retaining shared-body work cross terms.
/// Five scalar rows per hinge count toward the original 32-row limit; no
/// dependent sixth alignment equation is added. The complete point reaction
/// and complete hinge couple have separate vector bounds described by
/// [`RevoluteJoint::new`]. Existing distance/spherical-only entry points retain
/// their original native receipt type and heap payload.
///
/// This is a local smooth bilateral hinge, with finite bodies on both sides.
/// It does not add joint friction, motors, stops, backlash, impacts, a fixed
/// support, or a global/energy/symplectic certificate. Native work and momentum
/// defects remain visible. Cancellation or refusal publishes no partial state.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_articulated_graph(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    integrator: RigidBodyIntegrator,
    wrenches: &[Wrench],
    duration_s: f64,
    controls: DistanceGraphControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedGraphStep, DistanceGraphError> {
    poll(&mut cancelled)?;
    let maximum = controls.maximum_constraints.min(MAX_DISTANCE_ROWS);
    limit("joints", joints.len() as u128, maximum as u128)?;
    let mut rows = [Row::EMPTY; MAX_DISTANCE_ROWS];
    let mut row_count = 0;
    for &joint in joints {
        poll(&mut cancelled)?;
        let width = joint.row_count();
        limit("constraints", (row_count + width) as u128, maximum as u128)?;
        joint.write_rows(&mut rows[row_count..row_count + width]);
        row_count += width;
    }
    let extra = joints.len() as u128 * core::mem::size_of::<JointRows>() as u128;
    let bytes = planned_bytes_for::<JointReactionImpulse>(states.len(), row_count) + extra;
    limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut native_controls = controls;
    native_controls.maximum_live_bytes -= extra as usize;
    let (step, _) = advance_between(
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
    finish(step, joints, bytes, &mut cancelled)
}

pub(super) fn finish<FreeStep>(
    step: DistanceGraphStep<JointReactionImpulse, FreeStep>,
    joints: &[ArticulatedJoint],
    planned_live_bytes: u128,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedGraphStep<FreeStep>, DistanceGraphError> {
    let mut joint_rows = reserved(joints.len())?;
    let mut first_row = 0;
    for &joint in joints {
        poll(cancelled)?;
        let row_count = joint.row_count();
        joint_rows.push(JointRows {
            first_row,
            row_count,
        });
        first_row += row_count;
    }
    let (constraint_work_j, constraint_angular_impulse_world_nms) =
        reaction_account(&step, cancelled)?;
    poll(cancelled)?;
    Ok(JointGraphStep {
        step,
        joint_rows,
        constraint_work_j,
        constraint_angular_impulse_world_nms,
        planned_live_bytes: planned_live_bytes as usize,
    })
}

pub(super) fn reaction_account<FreeStep>(
    step: &DistanceGraphStep<JointReactionImpulse, FreeStep>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, Vec3), DistanceGraphError> {
    let mut work = Dd::ZERO;
    let mut angular = [Dd::ZERO; 3];
    for event in step.position_impulses.iter().chain(&step.velocity_impulses) {
        poll(cancelled)?;
        match event {
            JointReactionImpulse::Point(event) => {
                let before = driven::relative_velocity(
                    event.first.point_kinematics_before,
                    event.second.point_kinematics_before,
                );
                let after = driven::relative_velocity(
                    event.first.point_kinematics_after,
                    event.second.point_kinematics_after,
                );
                work = work + relative_work(before, after, event.second.impulse_world);
                let native = event.angular_impulse_world;
                for (sum, value) in angular.iter_mut().zip(components(native)) {
                    *sum = *sum + Dd::from_f64(value);
                }
            }
            JointReactionImpulse::Couple(event) => {
                work = work + Dd::from_f64(event.impulse_work_j);
            }
        }
    }
    let constraint_work_j = work.to_f64();
    let constraint_angular_impulse_world_nms = vector(angular.map(Dd::to_f64));
    if !constraint_work_j.is_finite() || !constraint_angular_impulse_world_nms.is_finite() {
        return Err(DistanceStepError::Invalid("articulated graph reaction account").into());
    }
    Ok((constraint_work_j, constraint_angular_impulse_world_nms))
}

#[derive(Clone, Copy)]
pub(super) struct AngularConstraint {
    axes_body: [Vec3; 2],
    transverse_second_body: Vec3,
    pub length_m: f64,
    pub component: usize,
}

pub(super) fn geometry(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    row: AngularConstraint,
) -> Result<Geometry, DistanceStepError> {
    let orientation = states.map(|state| state.pose().orientation());
    let a = orientation[0].rotate_body_to_world(row.axes_body[0]);
    let b = orientation[1].rotate_body_to_world(row.axes_body[1]);
    if a.dot(b) <= 0.0 {
        return Err(DistanceStepError::LocalBranch);
    }
    let transverse = orientation[1].rotate_body_to_world(row.transverse_second_body);
    let normal = transverse.cross(a).scale(row.length_m);
    let gap = (components(a)
        .into_iter()
        .zip(components(transverse))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
        * Dd::from_f64(row.length_m))
    .to_f64();
    let angular_velocity = [
        states[0]
            .point_kinematics(properties[0], Vec3::ZERO)?
            .angular_velocity_world,
        states[1]
            .point_kinematics(properties[1], Vec3::ZERO)?
            .angular_velocity_world,
    ];
    let rate = difference(angular_velocity[1], angular_velocity[0])
        .into_iter()
        .zip(components(normal))
        .fold(Dd::ZERO, |sum, (w, n)| sum + w * Dd::from_f64(n))
        .to_f64();
    if !normal.is_finite() || !gap.is_finite() || !rate.is_finite() {
        return Err(DistanceStepError::Invalid("hinge angular geometry"));
    }
    Ok(Geometry { normal, gap, rate })
}

impl NativeReaction for JointReactionImpulse {
    fn point(event: ActionReactionImpulseReceipt) -> Self {
        Self::Point(event)
    }
    fn couple(
        states: [RigidBodyState; 2],
        properties: [MassProperties; 2],
        impulse_world_nms: Vec3,
    ) -> Result<Self, DistanceGraphError> {
        let orientation = states.map(|state| state.pose().orientation());
        let first = states[0].apply_momentum_impulse(
            properties[0],
            Vec3::ZERO,
            orientation[0].rotate_world_to_body(impulse_world_nms.scale(-1.0)),
        )?;
        let second = states[1].apply_momentum_impulse(
            properties[1],
            Vec3::ZERO,
            orientation[1].rotate_world_to_body(impulse_world_nms),
        )?;
        let events = [first, second];
        let mut before = [Vec3::ZERO; 2];
        let mut after = [Vec3::ZERO; 2];
        let mut change = Dd::ZERO;
        let mut angular = [Dd::ZERO; 3];
        for i in 0..2 {
            let old = properties[i]
                .angular_velocity_body_checked(events[i].state_before.angular_momentum_body())?;
            let new = properties[i]
                .angular_velocity_body_checked(events[i].state_after.angular_momentum_body())?;
            change = change
                + crate::midpoint_impulse_pairing(old, new, events[i].angular_increment_body);
            before[i] = orientation[i].rotate_body_to_world(old);
            after[i] = orientation[i].rotate_body_to_world(new);
            let delta = orientation[i].rotate_body_to_world(events[i].angular_increment_body);
            for (sum, value) in angular.iter_mut().zip(components(delta)) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        let work = relative_work(
            difference(before[1], before[0]),
            difference(after[1], after[0]),
            impulse_world_nms,
        );
        let angular_momentum_change_world_nms = vector(angular.map(Dd::to_f64));
        let kinetic_energy_change_j = change.to_f64();
        let impulse_work_j = work.to_f64();
        let work_energy_residual_j = (change - work).to_f64();
        if !angular_momentum_change_world_nms.is_finite()
            || !kinetic_energy_change_j.is_finite()
            || !impulse_work_j.is_finite()
            || !work_energy_residual_j.is_finite()
        {
            return Err(DistanceStepError::Invalid("native hinge couple account").into());
        }
        Ok(Self::Couple(CoupleImpulseReceipt {
            first,
            second,
            couple_on_second_world_nms: impulse_world_nms,
            angular_momentum_change_world_nms,
            kinetic_energy_change_j,
            impulse_work_j,
            work_energy_residual_j,
        }))
    }
    fn states_after(&self) -> [RigidBodyState; 2] {
        match self {
            Self::Point(event) => [event.first.state_after, event.second.state_after],
            Self::Couple(event) => [event.first.state_after, event.second.state_after],
        }
    }
    fn impulse_work_j(&self) -> f64 {
        match self {
            Self::Point(event) => event.impulse_work_j,
            Self::Couple(event) => event.impulse_work_j,
        }
    }
    fn add_drive_work(
        work: &mut driven::Work,
        events: &[Self],
        impulses: &[f64],
        targets: &[driven::Target],
        endpoint: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DistanceGraphError> {
        for (i, ((event, impulse), target)) in events.iter().zip(impulses).zip(targets).enumerate()
        {
            poll(cancelled)?;
            match event {
                Self::Point(event) => work.add_point(i, event, *impulse, *target, endpoint)?,
                Self::Couple(event) => {
                    // Hinges are passive. Only distance rows own prescribed
                    // rates; retain the physical couple work in their defect.
                    work.add_work(i, Dd::ZERO, Dd::from_f64(event.impulse_work_j))?;
                }
            }
        }
        Ok(())
    }
}

fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn difference(a: Vec3, b: Vec3) -> [Dd; 3] {
    let a = components(a);
    let b = components(b);
    core::array::from_fn(|i| Dd::from_f64(a[i]) - Dd::from_f64(b[i]))
}
fn relative_work(before: [Dd; 3], after: [Dd; 3], impulse: Vec3) -> Dd {
    before
        .into_iter()
        .zip(after)
        .zip(components(impulse))
        .fold(Dd::ZERO, |sum, ((a, b), j)| {
            sum + (a + b) * Dd::from_f64(0.5) * Dd::from_f64(j)
        })
}
