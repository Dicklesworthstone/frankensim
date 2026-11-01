//! Payload capture checked against scalar mechanics, independent COM momentum
//! sums, and parallel-axis inertia assembly, without a second dynamics solver.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::{
    PayloadCapture, PayloadControls, PayloadError, PayloadState, capture_payload,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Kinematics, Link,
    SpatialInertia, forward_kinematics,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * (1. + actual.abs().max(expected.abs())),
        "{actual:.16e} != {expected:.16e}",
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}
fn pose(axes: Vec3, position: Vec3) -> Se3 {
    Se3::from_parts(rotation(axes), position).unwrap()
}
fn inertia(mass: f64, center: Vec3, diagonal: Vec3, axes: Vec3) -> SpatialInertia {
    let r = rotation(axes).matrix();
    let principal = Mat3 {
        m: [diagonal.x, 0., 0., 0., diagonal.y, 0., 0., 0., diagonal.z],
    };
    SpatialInertia::new(mass, center, r.compose(principal).compose(r.transpose())).unwrap()
}
fn sphere(mass: f64, moment: f64) -> SpatialInertia {
    inertia(
        mass,
        Vec3::default(),
        Vec3::new(moment, moment, moment),
        Vec3::default(),
    )
}
fn controls(tolerances: &[f64], loss_budget: f64) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-10,
        linear_momentum_tolerance_ns: 2e-10,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-9,
        maximum_energy_exchange_j: loss_budget,
        maximum_constraint_work_j: 1e-8,
        minimum_reciprocal_condition: 1e-10,
    }
}

struct Observables {
    linear: Vec3,
    angular: Vec3,
    kinetic: f64,
    generalized: Vec<f64>,
}

fn kinematics(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Kinematics {
    forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap()
}

// The only native operation here is kinematics. Each body's COM velocity,
// spin and orbital momentum are evaluated directly, then projected onto
// ancestor joint screws. The incoming body participates at its target link.
fn physical(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    payload: Option<(usize, &PayloadState)>,
    boost: Vec3,
) -> Observables {
    let kinematics = kinematics(model, state);
    let mut bodies: Vec<_> = model
        .links()
        .iter()
        .enumerate()
        .map(|(i, link)| {
            (
                i,
                link.inertia(),
                kinematics.world_from_link[i],
                kinematics.body_twist[i],
            )
        })
        .collect();
    if let Some((target, incoming)) = payload {
        bodies.push((
            target,
            incoming.inertia,
            incoming.world_from_body,
            incoming.twist_body,
        ));
    }
    let mut linear = Vec3::default();
    let mut angular = Vec3::default();
    let mut about_base = Vec3::default();
    let mut kinetic = 0.;
    let mut momenta = Vec::new();
    for (target, value, pose, twist) in bodies {
        let x = pose.transform_point(value.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(value.center_of_mass()))
            .unwrap()
            - boost;
        let p = velocity.scale(value.mass());
        let spin_body = value.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        linear = linear + p;
        angular = angular + spin + x.cross(p);
        about_base = about_base + spin + (x - state.base.world_from_base.translation()).cross(p);
        kinetic += 0.5 * value.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
        momenta.push((target, x, p, spin));
    }
    let inverse = state.base.world_from_base.rotation().inverse();
    let a = inverse.rotate(about_base).unwrap();
    let p = inverse.rotate(linear).unwrap();
    let mut generalized = vec![a.x, a.y, a.z, p.x, p.y, p.z];
    generalized.resize(6 + model.dof_count(), 0.);
    for (joint, link) in model.links().iter().enumerate() {
        let Some(dof) = model.dof_index(joint) else {
            continue;
        };
        let pose = kinematics.world_from_link[joint];
        let screw = link.joint().motion_subspace();
        let a = pose.rotation().rotate(screw.angular).unwrap();
        let v = pose.rotation().rotate(screw.linear).unwrap();
        for &(target, x, p, spin) in &momenta {
            let mut ancestor = Some(target);
            while let Some(index) = ancestor {
                if index == joint {
                    generalized[6 + dof] +=
                        a.dot(spin + (x - pose.translation()).cross(p)) + v.dot(p);
                    break;
                }
                ancestor = model.links()[index].parent();
            }
        }
    }
    Observables {
        linear,
        angular,
        kinetic,
        generalized,
    }
}

fn check_capture(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target: usize,
    payload: &PayloadState,
    rows: &[JointConstraint],
    result: &PayloadCapture,
) {
    let before = physical(model, state, Some((target, payload)), Vec3::default());
    let next = &result.state_after;
    let after = physical(&result.model_after, next, None, Vec3::default());
    let balance = &result.balance;
    vector(after.linear, before.linear);
    vector(after.angular, before.angular);
    vector(balance.momentum_before.linear_momentum_world, before.linear);
    vector(
        balance.momentum_before.angular_momentum_world,
        before.angular,
    );
    vector(balance.momentum_after.linear_momentum_world, after.linear);
    vector(balance.momentum_after.angular_momentum_world, after.angular);
    assert_eq!(next.base.world_from_base, state.base.world_from_base);
    assert_eq!(next.joint_position, state.joint_position);
    assert_eq!(balance.multipliers.len(), rows.len());
    assert_eq!(
        balance.generalized_constraint_impulse.len(),
        before.generalized.len()
    );
    for i in 0..before.generalized.len() {
        let impulse = if i < 6 {
            0.
        } else {
            rows.iter()
                .zip(&balance.multipliers)
                .map(|(row, lambda)| row.coefficients[i - 6] * lambda)
                .sum()
        };
        close(balance.generalized_constraint_impulse[i], impulse);
        close(after.generalized[i] - before.generalized[i], impulse);
        close(
            balance.generalized_momentum_before[i],
            before.generalized[i],
        );
        close(balance.generalized_momentum_after[i], after.generalized[i]);
    }
    let old_common = physical(
        model,
        state,
        Some((target, payload)),
        balance.accounting_boost_world,
    );
    let new_common = physical(
        &result.model_after,
        next,
        None,
        balance.accounting_boost_world,
    );
    close(balance.kinetic_before_j, old_common.kinetic);
    close(balance.kinetic_after_j, new_common.kinetic);
    close(
        balance.energy_residual_j,
        new_common.kinetic - old_common.kinetic + result.removed_relative_energy_j,
    );
    close(balance.energy_residual_j, 0.);
    close(balance.reaction_work_j, 0.);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    assert!(result.removed_relative_energy_j >= 0.);
    let kinematics = kinematics(&result.model_after, next);
    let target_pose = kinematics.world_from_link[target];
    let target_twist = kinematics.body_twist[target];
    let angular = target_pose.rotation().rotate(target_twist.angular).unwrap();
    let linear = target_pose.rotation().rotate(target_twist.linear).unwrap();
    let attached = &result.payload_after;
    vector(
        attached.world_from_body.translation(),
        payload.world_from_body.translation(),
    );
    vector(
        attached
            .world_from_body
            .rotation()
            .rotate(attached.twist_body.angular)
            .unwrap(),
        angular,
    );
    vector(
        attached
            .world_from_body
            .rotation()
            .rotate(attached.twist_body.linear)
            .unwrap(),
        linear + angular.cross(attached.world_from_body.translation() - target_pose.translation()),
    );
}

fn coax(
    limits: Option<JointLimits>,
    rate: f64,
) -> (ArticulatedModel, FreeFloatingMotionState, PayloadState) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2., 2.),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), limits).unwrap(),
            sphere(1., 1.),
        ),
    ])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0., 0., 0.4), Vec3::new(0.5, 0., 0.)),
        ),
        joint_position: vec![0.],
        joint_velocity: vec![rate],
    };
    let payload = PayloadState {
        inertia: sphere(3., 3.),
        world_from_body: Se3::identity(),
        twist_body: Twist::new(Vec3::new(0., 0., -0.8), Vec3::new(-1., 0., 0.)),
    };
    (model, state, payload)
}

#[test]
fn coax_capture_preserves_momenta_and_resolves_small_loss_under_enormous_common_boost() {
    let (model, state, payload) = coax(None, 1.2);
    let settings = controls(&[2e-10], 4.);
    let result = capture_payload(&model, &state, 1, payload, &[], settings, || false).unwrap();
    check_capture(&model, &state, 1, &payload, &[], &result);
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 0.4),
    );
    vector(
        result.state_after.base.twist_body.linear,
        Vec3::new(-0.25, 0., 0.),
    );
    close(result.state_after.joint_velocity[0], -0.6);
    close(result.payload_after.twist_body.angular.z, -0.2);
    close(result.removed_relative_energy_j, 2.16 + 1.6875);
    assert_eq!(result.model_after.links()[1].inertia(), sphere(4., 4.));
    assert_eq!(result.attachment.target_link(), 1);
    assert_eq!(result.attachment.payload_inertia(), payload.inertia);
    assert_eq!(
        result.attachment.retained_inertia(),
        model.links()[1].inertia()
    );
    assert_eq!(
        result.attachment.loaded_inertia(),
        result.model_after.links()[1].inertia()
    );
    let dt = 1e-3;
    let continued = step_free_floating(
        &result.model_after,
        &result.state_after,
        FreeFloatingMotionLoads {
            generalized_force: &[0.],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        dt,
        || false,
    )
    .unwrap();
    close(continued.state_after.base.twist_body.angular.z, 0.4);
    close(continued.state_after.joint_velocity[0], -0.6);
    close(continued.state_after.joint_position[0], -0.6 * dt);
    assert!(continued.kinetic_work_residual.abs() < 1e-8);

    // All boosted inputs and the expected quarter-unit output are exactly
    // representable; lab kinetic energy exceeds 1e18 and cannot resolve 4 J.
    let boost = Vec3::new(2_f64.powi(30), -2_f64.powi(29), 2_f64.powi(28));
    let mut moving_state = state.clone();
    moving_state.base.twist_body.linear = moving_state.base.twist_body.linear + boost;
    let mut moving_payload = payload;
    moving_payload.twist_body.linear = moving_payload.twist_body.linear + boost;
    let moving = capture_payload(
        &model,
        &moving_state,
        1,
        moving_payload,
        &[],
        settings,
        || false,
    )
    .unwrap();
    vector(
        moving.state_after.base.twist_body.linear - boost,
        result.state_after.base.twist_body.linear,
    );
    vector(
        moving.state_after.base.twist_body.angular,
        result.state_after.base.twist_body.angular,
    );
    close(moving.state_after.joint_velocity[0], -0.6);
    vector(
        moving.payload_after.twist_body.linear - boost,
        result.payload_after.twist_body.linear,
    );
    close(moving.removed_relative_energy_j, 3.8475);
    close(moving.balance.kinetic_before_j, 5.775);
    close(moving.balance.kinetic_after_j, 1.9275);
    assert!(moving.balance.energy_residual_j.abs() < 2e-9);
    assert!(moving.balance.maximum_normalized_momentum_residual <= 1.);
}

fn serial_tree() -> (ArticulatedModel, FreeFloatingMotionState) {
    let centers = [
        Vec3::new(0.12, -0.21, 0.08),
        Vec3::new(-0.08, 0.17, 0.04),
        Vec3::new(0.11, 0.03, -0.09),
        Vec3::new(0.04, -0.12, 0.06),
    ];
    let masses = [2.3, 1.1, 0.8, 0.4];
    let diagonals = [
        Vec3::new(1.2, 1.6, 2.),
        Vec3::new(0.5, 0.7, 0.9),
        Vec3::new(0.3, 0.45, 0.6),
        Vec3::new(0.2, 0.25, 0.3),
    ];
    let inertias: Vec<_> = (0..4)
        .map(|i| {
            inertia(
                masses[i],
                centers[i],
                diagonals[i],
                Vec3::new(0.2 + 0.1 * i as f64, -0.15, 0.3),
            )
        })
        .collect();
    let mounts = [
        pose(Vec3::new(0.2, -0.1, 0.15), Vec3::new(0.17, -0.04, 0.11)),
        pose(Vec3::new(0.1, 0.25, -0.15), Vec3::new(0.6, 0.08, -0.11)),
        pose(Vec3::new(-0.18, 0.11, 0.2), Vec3::new(0.1, -0.35, 0.27)),
        pose(Vec3::new(0.11, -0.12, 0.18), Vec3::new(-0.2, 0.17, 0.33)),
    ];
    let joints = [
        JointModel::FIXED,
        JointModel::revolute(Vec3::new(1., 2., -1.), None).unwrap(),
        JointModel::prismatic(Vec3::new(-1., 1., 2.), None).unwrap(),
        JointModel::helical(Vec3::new(2., -1., 1.), 0.13, None).unwrap(),
    ];
    let model = ArticulatedModel::new(
        (0_usize..4)
            .map(|i| {
                Link::new(
                    ["carrier", "rotor", "slider", "screw"][i],
                    i.checked_sub(1),
                    mounts[i],
                    joints[i],
                    inertias[i],
                )
            })
            .collect(),
    )
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(Vec3::new(0.31, -0.27, 0.19), Vec3::new(2., -1., 0.6)),
            Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(-0.2, 0.6, 0.35)),
        ),
        joint_position: vec![0.2, 0.1, 0.15],
        joint_velocity: vec![0.4, 0.2, 0.3],
    };
    (model, state)
}

fn check_composite(
    old: SpatialInertia,
    payload: &PayloadState,
    relative: Se3,
    actual: SpatialInertia,
) {
    let payload_center = relative
        .transform_point(payload.inertia.center_of_mass())
        .unwrap();
    let mass = old.mass() + payload.inertia.mass();
    let center = (old.center_of_mass().scale(old.mass())
        + payload_center.scale(payload.inertia.mass()))
    .scale(1. / mass);
    close(actual.mass(), mass);
    vector(actual.center_of_mass(), center);
    let r = relative.rotation().matrix();
    let rotated = r
        .compose(payload.inertia.inertia_com())
        .compose(r.transpose());
    let mut expected = [0.; 9];
    for (part_mass, part_center, tensor) in [
        (old.mass(), old.center_of_mass(), old.inertia_com()),
        (payload.inertia.mass(), payload_center, rotated),
    ] {
        let d = part_center - center;
        let components = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in 0..3 {
                expected[3 * i + j] += tensor.m[3 * i + j]
                    + part_mass
                        * (if i == j { d.dot(d) } else { 0. } - components[i] * components[j]);
            }
        }
    }
    for (actual, expected) in actual.inertia_com().m.iter().zip(expected) {
        close(*actual, expected);
    }
}

#[test]
fn offset_full_tensor_capture_obeys_mixed_joint_momenta_and_scene_covariance() {
    let (model, state) = serial_tree();
    let target = 3;
    let settings = controls(&[2e-10; 3], 20.);
    let kinematics = kinematics(&model, &state);
    let relative = pose(Vec3::new(0.25, -0.17, 0.31), Vec3::new(0.28, -0.19, 0.23));
    let payload = PayloadState {
        inertia: inertia(
            0.7,
            Vec3::new(0.08, -0.06, 0.11),
            Vec3::new(0.18, 0.24, 0.3),
            Vec3::new(-0.2, 0.3, 0.15),
        ),
        world_from_body: kinematics.world_from_link[target]
            .compose(relative)
            .unwrap(),
        twist_body: Twist::new(Vec3::new(-0.3, 0.55, 0.2), Vec3::new(0.7, -0.4, 0.15)),
    };
    let result = capture_payload(&model, &state, target, payload, &[], settings, || false).unwrap();
    check_capture(&model, &state, target, &payload, &[], &result);
    check_composite(
        model.links()[target].inertia(),
        &payload,
        relative,
        result.model_after.links()[target].inertia(),
    );
    for i in 0..target {
        assert_eq!(model.links()[i], result.model_after.links()[i]);
    }
    assert!(
        (result.state_after.base.twist_body.angular - state.base.twist_body.angular).norm() > 1e-3
    );
    let scene = pose(Vec3::new(-0.4, 0.2, 0.3), Vec3::new(4., -2.5, 1.7));
    for boost in [Vec3::default(), Vec3::new(3., -2., 4.)] {
        let mut transformed = state.clone();
        transformed.base.world_from_base = scene.compose(state.base.world_from_base).unwrap();
        let base_boost = transformed
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
        transformed.base.twist_body.linear = state.base.twist_body.linear + base_boost;
        let mut incoming = payload;
        incoming.world_from_body = scene.compose(payload.world_from_body).unwrap();
        let payload_boost = incoming
            .world_from_body
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
        incoming.twist_body.linear = payload.twist_body.linear + payload_boost;
        let captured = capture_payload(
            &model,
            &transformed,
            target,
            incoming,
            &[],
            settings,
            || false,
        )
        .unwrap();
        check_capture(&model, &transformed, target, &incoming, &[], &captured);
        check_composite(
            model.links()[target].inertia(),
            &incoming,
            relative,
            captured.model_after.links()[target].inertia(),
        );
        vector(
            captured.state_after.base.twist_body.angular,
            result.state_after.base.twist_body.angular,
        );
        vector(
            captured.state_after.base.twist_body.linear,
            result.state_after.base.twist_body.linear + base_boost,
        );
        for (actual, expected) in captured
            .state_after
            .joint_velocity
            .iter()
            .zip(&result.state_after.joint_velocity)
        {
            close(*actual, *expected);
        }
        vector(
            captured.payload_after.twist_body.angular,
            result.payload_after.twist_body.angular,
        );
        vector(
            captured.payload_after.twist_body.linear,
            result.payload_after.twist_body.linear + payload_boost,
        );
        close(
            captured.removed_relative_energy_j,
            result.removed_relative_energy_j,
        );
    }
}

#[test]
fn locked_capture_admits_the_constrained_endpoint_and_refuses_atomically() {
    let limits = Some(JointLimits::new(-2., 2., 0.5, 100.).unwrap());
    let (model, state, payload) = coax(limits, 0.);
    let retained = (model.clone(), state.clone(), payload);
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![1.],
        position: 0.,
        position_tolerance: 2e-10,
        velocity_tolerance: 2e-10,
    }];
    let settings = controls(&[2e-10], 3.);
    // The free recoil would be omega=.4 and qdot=-.9, beyond the .5
    // speed limit. Only the locked endpoint is a physically valid candidate.
    assert!(capture_payload(&model, &state, 1, payload, &[], settings, || false).is_err());
    let result = capture_payload(&model, &state, 1, payload, &rows, settings, || false).unwrap();
    check_capture(&model, &state, 1, &payload, &rows, &result);
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., -0.2),
    );
    vector(
        result.state_after.base.twist_body.linear,
        Vec3::new(-0.25, 0., 0.),
    );
    close(result.state_after.joint_velocity[0], 0.);
    close(result.balance.multipliers[0], 1.2);
    close(result.balance.generalized_constraint_impulse[6], 1.2);
    close(result.removed_relative_energy_j, 1.08 + 1.6875);
    assert!(matches!(
        capture_payload(
            &model,
            &state,
            1,
            payload,
            &rows,
            controls(&[2e-10], 2.7),
            || false
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert!(
        capture_payload(
            &model,
            &state,
            model.links().len(),
            payload,
            &rows,
            settings,
            || false
        )
        .is_err()
    );
    let mut wrong_position = state.clone();
    wrong_position.joint_position[0] = 0.01;
    let mut wrong_rate = state.clone();
    wrong_rate.joint_velocity[0] = 0.01;
    for invalid in [&wrong_position, &wrong_rate] {
        let snapshot = invalid.clone();
        assert!(matches!(
            capture_payload(&model, invalid, 1, payload, &rows, settings, || { false }),
            Err(PayloadError::Constraint(_))
        ));
        assert_eq!(*invalid, snapshot);
    }
    assert!(matches!(
        capture_payload(&model, &state, 1, payload, &rows, settings, || true),
        Err(PayloadError::Cancelled)
    ));
    assert_eq!((model, state, payload), retained);
}
