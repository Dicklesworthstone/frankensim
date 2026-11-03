//! Supported payload catches: scalar oracles and independent rigid-body sums.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::inventory::supported::{
    SupportedPayloadCapture, SupportedPayloadControls, capture_payload,
};
use fs_mbd::articulated::motion::inventory::{PayloadError, PayloadState};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionState, PrescribedBasePath, step_prescribed_base,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia,
    forward_kinematics,
};
use std::cell::Cell;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-9 * (1.0 + a.abs().max(b.abs())),
        "{a:.16e} != {b:.16e}"
    );
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn rotation(axis: Vec3) -> So3 {
    So3::exp(So3Tangent::new(axis)).unwrap()
}
fn pose(axis: Vec3, origin: Vec3) -> Se3 {
    Se3::from_parts(rotation(axis), origin).unwrap()
}
fn inertia(mass: f64, center: Vec3, principal: Vec3, axes: Vec3) -> SpatialInertia {
    let r = rotation(axes).matrix();
    let diagonal = Mat3 {
        m: [
            principal.x,
            0.,
            0.,
            0.,
            principal.y,
            0.,
            0.,
            0.,
            principal.z,
        ],
    };
    SpatialInertia::new(mass, center, r.compose(diagonal).compose(r.transpose())).unwrap()
}
fn sphere(mass: f64) -> SpatialInertia {
    inertia(
        mass,
        Vec3::default(),
        Vec3::new(1., 1., 1.),
        Vec3::default(),
    )
}
fn model(joint: JointModel, mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "catcher",
        None,
        Se3::identity(),
        joint,
        sphere(mass),
    )])
    .unwrap()
}
fn state(base: BaseState, position: Vec<f64>, velocity: Vec<f64>) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 4.,
        base,
        joint_position: position,
        joint_velocity: velocity,
    }
}
fn payload(mass: f64, velocity: Vec3) -> PayloadState {
    PayloadState {
        inertia: sphere(mass),
        world_from_body: Se3::identity(),
        twist_body: Twist::new(Vec3::default(), velocity),
    }
}
fn controls(tolerances: &[f64]) -> SupportedPayloadControls<'_> {
    SupportedPayloadControls {
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_removed_relative_energy_j: 1000.,
    }
}

struct Physical {
    linear: Vec3,
    angular: Vec3,
    energy: f64,
    joint: Vec<f64>,
}

// Only native FK is shared with the implementation. Momentum and kinetic
// energy use independent COM/spin formulas, and generalized momenta project
// those physical quantities on each ancestor joint's world screw.
fn physical(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    payload: Option<(usize, PayloadState)>,
) -> Physical {
    let kin = forward_kinematics(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let mut bodies: Vec<_> = model
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| (i, link.inertia(), kin.world_from_link[i], kin.body_twist[i]))
        .collect();
    if let Some((target, incoming)) = payload {
        bodies.push((
            target,
            incoming.inertia,
            incoming.world_from_body,
            incoming.twist_body,
        ));
    }
    let mut result = Physical {
        linear: Vec3::default(),
        angular: Vec3::default(),
        energy: 0.,
        joint: vec![0.; model.dof_count()],
    };
    for (target, mass, pose, twist) in bodies {
        let center = pose.transform_point(mass.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(mass.center_of_mass()))
            .unwrap();
        let linear = velocity.scale(mass.mass());
        let spin_body = mass.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        result.linear = result.linear + linear;
        result.angular = result.angular + spin + center.cross(linear);
        result.energy +=
            0.5 * mass.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
        let mut ancestor = Some(target);
        while let Some(index) = ancestor {
            if let Some(dof) = model.dof_index(index) {
                let screw = model.links()[index].joint().motion_subspace();
                let joint_pose = kin.world_from_link[index];
                let angular = joint_pose.rotation().rotate(screw.angular).unwrap();
                let linear_axis = joint_pose.rotation().rotate(screw.linear).unwrap();
                result.joint[dof] += angular
                    .dot(spin + (center - joint_pose.translation()).cross(linear))
                    + linear_axis.dot(linear);
            }
            ancestor = model.links()[index].parent();
        }
    }
    result
}

fn audit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target: usize,
    payload: PayloadState,
    result: &SupportedPayloadCapture,
) {
    let before = physical(model, state, Some((target, payload)));
    let after = physical(&result.model_after, &result.state_after, None);
    let balance = &result.balance;
    vector(balance.momentum_before.linear_momentum_world, before.linear);
    vector(
        balance.momentum_before.angular_momentum_world,
        before.angular,
    );
    close(balance.momentum_before.kinetic_energy, before.energy);
    vector(balance.momentum_after.linear_momentum_world, after.linear);
    vector(balance.momentum_after.angular_momentum_world, after.angular);
    close(balance.momentum_after.kinetic_energy, after.energy);
    vector(
        balance.linear_momentum_change_world,
        after.linear - before.linear,
    );
    vector(
        balance.angular_momentum_change_world,
        after.angular - before.angular,
    );
    vector(
        balance.support_impulse_world.force,
        after.linear - before.linear,
    );
    vector(
        balance.support_impulse_world.torque,
        after.angular - before.angular,
    );
    vector(balance.linear_momentum_residual_world, Vec3::default());
    vector(balance.angular_momentum_residual_world, Vec3::default());
    for ((old, new), residual) in before
        .joint
        .iter()
        .zip(&after.joint)
        .zip(&balance.joint_momentum_residual)
    {
        close(*new, *old);
        close(*residual, new - old);
    }
    close(
        balance.kinetic_energy_change_j,
        after.energy - before.energy,
    );
    close(
        balance.kinetic_energy_change_j + result.removed_relative_energy_j,
        balance.support_work_j,
    );
    close(
        balance.support_work_j,
        balance
            .support_impulse_body
            .torque
            .dot(state.base.twist_body.angular)
            + balance
                .support_impulse_body
                .force
                .dot(state.base.twist_body.linear),
    );
    close(
        balance.accounting_kinetic_energy_change_j,
        balance.kinetic_after_j - balance.kinetic_before_j,
    );
    close(
        balance.accounting_support_work_j,
        balance.support_work_j
            - balance
                .accounting_boost_world
                .dot(balance.support_impulse_world.force),
    );
    close(balance.energy_residual_j, 0.);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    assert_eq!(state.time_s, result.state_after.time_s);
    assert_eq!(state.base, result.state_after.base);
    assert_eq!(state.joint_position, result.state_after.joint_position);
}

#[test]
fn fixed_support_absorbs_incoming_momentum_and_continues_with_loaded_mass() {
    let model = model(JointModel::FIXED, 3.);
    let state = state(BaseState::stationary(Se3::identity()), vec![], vec![]);
    let incoming = payload(2., Vec3::new(3., 0., 0.));
    let result = capture_payload(&model, &state, 0, incoming, controls(&[]), || false).unwrap();
    close(result.model_after.links()[0].inertia().mass(), 5.);
    close(result.removed_relative_energy_j, 9.);
    vector(
        result.balance.support_impulse_world.force,
        Vec3::new(-6., 0., 0.),
    );
    close(result.balance.support_work_j, 0.);
    close(result.balance.kinetic_energy_change_j, -9.);
    assert_eq!(
        result.attachment.retained_inertia(),
        model.links()[0].inertia()
    );
    audit(&model, &state, 0, incoming, &result);
    let next = step_prescribed_base(
        &result.model_after,
        &result.state_after,
        PrescribedBasePath::stationary(Se3::identity()),
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::new(0., -10., 0.),
        },
        0.01,
        || false,
    )
    .unwrap();
    close(next.state_after.time_s, 4.01);
    vector(next.support_wrench_body.force, Vec3::new(0., 50., 0.));
}

#[test]
fn moving_support_does_work_while_capturing_a_stationary_payload() {
    let model = model(JointModel::FIXED, 3.);
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(2., 0., 0.)),
        Twist::zero(),
    );
    let state = state(base, vec![], vec![]);
    let incoming = payload(2., Vec3::default());
    let result = capture_payload(&model, &state, 0, incoming, controls(&[]), || false).unwrap();
    vector(
        result.payload_after.twist_body.linear,
        Vec3::new(2., 0., 0.),
    );
    vector(
        result.balance.support_impulse_world.force,
        Vec3::new(4., 0., 0.),
    );
    close(result.balance.support_work_j, 8.);
    close(result.balance.kinetic_energy_change_j, 4.);
    close(result.removed_relative_energy_j, 4.);
    close(result.balance.accounting_support_work_j, 0.);
    audit(&model, &state, 0, incoming, &result);
}

#[test]
fn sliding_root_uses_relative_joint_rate_with_a_driven_base() {
    let model = model(
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
        3.,
    );
    let state = state(
        BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(2., 0., 0.)),
            Twist::zero(),
        ),
        vec![0.4],
        vec![1.],
    );
    let incoming = PayloadState {
        world_from_body: pose(Vec3::default(), Vec3::new(0.4, 0., 0.)),
        ..payload(2., Vec3::new(-1., 0., 0.))
    };
    let result = capture_payload(&model, &state, 0, incoming, controls(&[2e-9]), || false).unwrap();
    close(result.state_after.joint_velocity[0], -0.6);
    close(result.payload_after.twist_body.linear.x, 1.4);
    vector(result.balance.support_impulse_world.force, Vec3::default());
    close(result.removed_relative_energy_j, 9.6);
    audit(&model, &state, 0, incoming, &result);
}

#[test]
fn rotating_catch_retains_full_offset_payload_tensor_and_support_couple() {
    let retained = inertia(
        3.,
        Vec3::new(-0.2, 0.3, -0.1),
        Vec3::new(3., 4., 5.),
        Vec3::new(0.3, -0.2, 0.4),
    );
    let model = ArticulatedModel::new(vec![Link::new(
        "hinged catcher",
        None,
        Se3::identity(),
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        retained,
    )])
    .unwrap();
    let state = state(BaseState::stationary(Se3::identity()), vec![0.], vec![0.7]);
    let incoming = PayloadState {
        inertia: inertia(
            1.8,
            Vec3::new(0.2, -0.1, 0.3),
            Vec3::new(0.4, 0.5, 0.7),
            Vec3::new(-0.2, 0.3, 0.1),
        ),
        world_from_body: pose(Vec3::new(0.5, 0.4, -0.3), Vec3::new(0.9, -0.6, 0.5)),
        twist_body: Twist::new(Vec3::new(0.4, -0.7, 0.6), Vec3::new(1.3, -0.8, 0.2)),
    };
    let before = physical(&model, &state, Some((0, incoming)));
    let r = incoming.world_from_body.rotation().matrix();
    let rotated = r
        .compose(incoming.inertia.inertia_com())
        .compose(r.transpose());
    let c0 = retained.center_of_mass();
    let cp = incoming
        .world_from_body
        .transform_point(incoming.inertia.center_of_mass())
        .unwrap();
    let j0 = retained.inertia_com().m[8] + retained.mass() * (c0.x * c0.x + c0.y * c0.y);
    let jp = rotated.m[8] + incoming.inertia.mass() * (cp.x * cp.x + cp.y * cp.y);
    let expected_rate = before.angular.z / (j0 + jp);
    let result = capture_payload(&model, &state, 0, incoming, controls(&[2e-9]), || false).unwrap();
    close(result.state_after.joint_velocity[0], expected_rate);
    close(result.balance.support_impulse_body.torque.z, 0.);
    assert!(result.balance.support_impulse_body.torque.x.abs() > 1e-4);
    let combined = result.model_after.links()[0].inertia();
    let center =
        (c0.scale(retained.mass()) + cp.scale(incoming.inertia.mass())).scale(1. / combined.mass());
    vector(combined.center_of_mass(), center);
    let old_arm = c0 - center;
    let new_arm = cp - center;
    let arrays = |v: Vec3| [v.x, v.y, v.z];
    let a = arrays(old_arm);
    let b = arrays(new_arm);
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1. } else { 0. };
            let expected = retained.inertia_com().m[3 * i + j]
                + rotated.m[3 * i + j]
                + retained.mass() * (delta * old_arm.dot(old_arm) - a[i] * a[j])
                + incoming.inertia.mass() * (delta * new_arm.dot(new_arm) - b[i] * b[j]);
            close(combined.inertia_com().m[3 * i + j], expected);
        }
    }
    audit(&model, &state, 0, incoming, &result);
}

#[test]
fn branched_tree_capture_retains_joint_momentum_with_rotating_support() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "rail",
            None,
            pose(Vec3::new(0.1, 0.2, -0.1), Vec3::new(0.2, -0.3, 0.4)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            sphere(4.),
        ),
        Link::new(
            "arm",
            Some(0),
            pose(Vec3::new(0.2, -0.3, 0.4), Vec3::new(0.6, 0.3, -0.2)),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            sphere(2.),
        ),
        Link::new(
            "sibling",
            Some(0),
            pose(Vec3::new(-0.2, 0.1, 0.3), Vec3::new(-0.5, 0.6, 0.2)),
            JointModel::helical(Vec3::new(0., 1., 0.), 0.1, None).unwrap(),
            sphere(1.),
        ),
    ])
    .unwrap();
    let state = state(
        BaseState::prescribed(
            pose(Vec3::new(0.3, -0.4, 0.2), Vec3::new(2., -3., 1.)),
            Twist::new(Vec3::new(0.2, -0.3, 0.4), Vec3::new(1.1, -0.6, 0.2)),
            Twist::zero(),
        ),
        vec![0.3, -0.4, 0.2],
        vec![0.7, -0.5, 0.9],
    );
    let kin = forward_kinematics(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let incoming = PayloadState {
        inertia: inertia(
            1.5,
            Vec3::new(0.1, -0.2, 0.3),
            Vec3::new(0.4, 0.5, 0.6),
            Vec3::new(0.2, -0.3, 0.4),
        ),
        world_from_body: kin.world_from_link[1]
            .compose(pose(Vec3::new(0.4, -0.2, 0.1), Vec3::new(0.7, -0.2, 0.3)))
            .unwrap(),
        twist_body: Twist::new(Vec3::new(-0.3, 0.6, 0.2), Vec3::new(0.3, -0.8, 1.2)),
    };
    let result =
        capture_payload(&model, &state, 1, incoming, controls(&[2e-9; 3]), || false).unwrap();
    audit(&model, &state, 1, incoming, &result);
    assert_ne!(
        result.state_after.joint_velocity[0],
        state.joint_velocity[0]
    );
    assert_ne!(
        result.state_after.joint_velocity[1],
        state.joint_velocity[1]
    );
    assert_ne!(
        result.state_after.joint_velocity[2],
        state.joint_velocity[2]
    );
}

#[test]
fn cancellation_and_loss_budget_leave_inventory_and_state_unchanged() {
    let model = model(JointModel::FIXED, 3.);
    let state = state(BaseState::stationary(Se3::identity()), vec![], vec![]);
    let incoming = payload(2., Vec3::new(3., 0., 0.));
    let saved_model = model.clone();
    let saved_state = state.clone();
    let polls = Cell::new(0);
    capture_payload(&model, &state, 0, incoming, controls(&[]), || {
        polls.set(polls.get() + 1);
        false
    })
    .unwrap();
    for stop in [1, polls.get()] {
        let seen = Cell::new(0);
        let error = capture_payload(&model, &state, 0, incoming, controls(&[]), || {
            seen.set(seen.get() + 1);
            seen.get() >= stop
        })
        .unwrap_err();
        assert_eq!(error, PayloadError::Cancelled);
    }
    let limited = SupportedPayloadControls {
        maximum_removed_relative_energy_j: 8.9,
        ..controls(&[])
    };
    assert!(matches!(
        capture_payload(&model, &state, 0, incoming, limited, || false),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert_eq!(model, saved_model);
    assert_eq!(state, saved_state);
}

#[test]
fn endpoint_joint_speed_limits_and_control_admission_remain_native() {
    let joint = JointModel::prismatic(
        Vec3::new(1., 0., 0.),
        Some(JointLimits::new(-2., 2., 0.2, 0.1).unwrap()),
    )
    .unwrap();
    let model = model(joint, 3.);
    let state = state(BaseState::stationary(Se3::identity()), vec![0.], vec![0.]);
    let incoming = payload(2., Vec3::new(2., 0., 0.));
    assert!(matches!(
        capture_payload(&model, &state, 0, incoming, controls(&[2e-9]), || false),
        Err(PayloadError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    let invalid = SupportedPayloadControls {
        linear_momentum_tolerance_ns: 0.,
        ..controls(&[2e-9])
    };
    assert!(capture_payload(&model, &state, 0, incoming, invalid, || false).is_err());
    assert!(capture_payload(&model, &state, 1, incoming, controls(&[2e-9]), || false).is_err());
    let mut invalid_state = state.clone();
    invalid_state.time_s = f64::NAN;
    assert!(
        capture_payload(
            &model,
            &invalid_state,
            0,
            incoming,
            controls(&[2e-9]),
            || false
        )
        .is_err()
    );
    assert_eq!(state.joint_velocity, vec![0.]);
    close(model.links()[0].inertia().mass(), 3.);
}

#[test]
fn common_world_rotation_translation_and_boost_preserve_the_capture() {
    let model = ArticulatedModel::new(vec![Link::new(
        "driven hinge",
        None,
        pose(Vec3::new(0.2, 0.1, -0.3), Vec3::new(-0.2, 0.4, 0.1)),
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        inertia(
            3.,
            Vec3::new(-0.2, 0.3, -0.1),
            Vec3::new(3., 4., 5.),
            Vec3::new(0.3, -0.2, 0.4),
        ),
    )])
    .unwrap();
    let state = state(
        BaseState::prescribed(
            pose(Vec3::new(0.3, -0.4, 0.2), Vec3::new(2., -3., 1.)),
            Twist::new(Vec3::new(0.2, -0.3, 0.4), Vec3::new(1.1, -0.6, 0.2)),
            Twist::zero(),
        ),
        vec![0.35],
        vec![-0.6],
    );
    let kin = forward_kinematics(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let incoming = PayloadState {
        inertia: inertia(
            1.5,
            Vec3::new(0.1, -0.2, 0.3),
            Vec3::new(0.4, 0.5, 0.6),
            Vec3::new(0.2, -0.3, 0.4),
        ),
        world_from_body: kin.world_from_link[0]
            .compose(pose(Vec3::new(0.4, -0.2, 0.1), Vec3::new(0.7, -0.2, 0.3)))
            .unwrap(),
        twist_body: Twist::new(Vec3::new(0.3, 0.2, -0.4), Vec3::new(0.8, -1.1, 0.5)),
    };
    let original =
        capture_payload(&model, &state, 0, incoming, controls(&[2e-9]), || false).unwrap();
    let world_change = pose(Vec3::new(-0.4, 0.6, 0.3), Vec3::new(-5., 3., 7.));
    let boost = Vec3::new(17., -9., 6.);
    let shifted_twist = |old: Twist, new_pose: Se3| {
        Twist::new(
            old.angular,
            old.linear + new_pose.rotation().inverse().rotate(boost).unwrap(),
        )
    };
    let shifted_base_pose = world_change.compose(state.base.world_from_base).unwrap();
    let shifted_state = PrescribedBaseMotionState {
        base: BaseState::prescribed(
            shifted_base_pose,
            shifted_twist(state.base.twist_body, shifted_base_pose),
            Twist::new(
                state.base.acceleration_body.angular,
                state.base.acceleration_body.linear
                    - state.base.twist_body.angular.cross(
                        shifted_base_pose
                            .rotation()
                            .inverse()
                            .rotate(boost)
                            .unwrap(),
                    ),
            ),
        ),
        ..state.clone()
    };
    let shifted_payload_pose = world_change.compose(incoming.world_from_body).unwrap();
    let shifted_payload = PayloadState {
        world_from_body: shifted_payload_pose,
        twist_body: shifted_twist(incoming.twist_body, shifted_payload_pose),
        ..incoming
    };
    let shifted = capture_payload(
        &model,
        &shifted_state,
        0,
        shifted_payload,
        controls(&[2e-9]),
        || false,
    )
    .unwrap();
    close(
        shifted.state_after.joint_velocity[0],
        original.state_after.joint_velocity[0],
    );
    close(
        shifted.removed_relative_energy_j,
        original.removed_relative_energy_j,
    );
    vector(
        shifted.balance.support_impulse_body.force,
        original.balance.support_impulse_body.force,
    );
    vector(
        shifted.balance.support_impulse_body.torque,
        original.balance.support_impulse_body.torque,
    );
    let force = world_change
        .rotation()
        .rotate(original.balance.support_impulse_world.force)
        .unwrap();
    let torque = world_change
        .rotation()
        .rotate(original.balance.support_impulse_world.torque)
        .unwrap()
        + world_change.translation().cross(force);
    vector(shifted.balance.support_impulse_world.force, force);
    vector(shifted.balance.support_impulse_world.torque, torque);
    close(
        shifted.balance.support_work_j,
        original.balance.support_work_j + boost.dot(force),
    );
    close(
        shifted.balance.kinetic_energy_change_j,
        original.balance.kinetic_energy_change_j + boost.dot(force),
    );
    vector(
        shifted.balance.accounting_boost_world,
        world_change
            .rotation()
            .rotate(original.balance.accounting_boost_world)
            .unwrap()
            + boost,
    );
    close(
        shifted.balance.kinetic_before_j,
        original.balance.kinetic_before_j,
    );
    close(
        shifted.balance.kinetic_after_j,
        original.balance.kinetic_after_j,
    );
    close(
        shifted.balance.accounting_support_work_j,
        original.balance.accounting_support_work_j,
    );
    audit(&model, &state, 0, incoming, &original);
    audit(&model, &shifted_state, 0, shifted_payload, &shifted);
}
