//! Physical inertia redistribution at unchanged articulated configuration.
//! The oracles sum individual COM momenta and project their actual world
//! moments onto downstream joint screws; they never use a generalized mass
//! matrix, a production momentum reconstruction, or an independent solver.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::{
    InertiaReconfiguration, ReconfigurationControls, ReconfigurationError, transfer_inertia,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, forward_kinematics,
};

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * (1. + a.abs().max(b.abs())),
        "{a:.16e} != {b:.16e}; difference {}",
        (a - b).abs(),
    );
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
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
fn controls(tolerances: &[f64], work: f64) -> ReconfigurationControls<'_> {
    ReconfigurationControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-11,
        linear_momentum_tolerance_ns: 2e-11,
        joint_momentum_tolerances: tolerances,
        maximum_parameter_work_j: work,
    }
}
fn with_inertias(model: &ArticulatedModel, values: &[SpatialInertia]) -> ArticulatedModel {
    assert_eq!(model.link_count(), values.len());
    ArticulatedModel::new(
        model
            .links()
            .iter()
            .zip(values)
            .map(|(link, value)| {
                Link::new(
                    link.name(),
                    link.parent(),
                    link.parent_from_child_zero(),
                    link.joint(),
                    *value,
                )
            })
            .collect(),
    )
    .unwrap()
}

struct Observables {
    linear: Vec3,
    angular: Vec3,
    center: Vec3,
    kinetic: f64,
    generalized: Vec<f64>,
}

fn physical(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Observables {
    let kinematics = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let mut linear = Vec3::default();
    let mut angular = Vec3::default();
    let mut angular_base_origin = Vec3::default();
    let mut center = Vec3::default();
    let mut total_mass = 0.;
    let mut kinetic = 0.;
    let mut momenta = Vec::new();
    for (i, link) in model.links().iter().enumerate() {
        let value = link.inertia();
        let pose = kinematics.world_from_link[i];
        let twist = kinematics.body_twist[i];
        let x = pose.transform_point(value.center_of_mass()).unwrap();
        let v = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(value.center_of_mass()))
            .unwrap();
        let p = v.scale(value.mass());
        let spin_body = value.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        linear = linear + p;
        angular = angular + spin + x.cross(p);
        angular_base_origin =
            angular_base_origin + spin + (x - state.base.world_from_base.translation()).cross(p);
        kinetic += 0.5 * value.mass() * v.dot(v) + 0.5 * twist.angular.dot(spin_body);
        center = center + x.scale(value.mass());
        total_mass += value.mass();
        momenta.push((x, p, spin));
    }
    let inverse = state.base.world_from_base.rotation().inverse();
    let local_angular = inverse.rotate(angular_base_origin).unwrap();
    let local_linear = inverse.rotate(linear).unwrap();
    let mut generalized = vec![
        local_angular.x,
        local_angular.y,
        local_angular.z,
        local_linear.x,
        local_linear.y,
        local_linear.z,
    ];
    generalized.resize(6 + model.dof_count(), 0.);
    for (joint, link) in model.links().iter().enumerate() {
        let Some(dof) = model.dof_index(joint) else {
            continue;
        };
        let pose = kinematics.world_from_link[joint];
        let screw = link.joint().motion_subspace();
        let angular_axis = pose.rotation().rotate(screw.angular).unwrap();
        let linear_axis = pose.rotation().rotate(screw.linear).unwrap();
        for (body, &(x, p, spin)) in momenta.iter().enumerate() {
            let mut ancestor = Some(body);
            while let Some(index) = ancestor {
                if index == joint {
                    generalized[6 + dof] += angular_axis
                        .dot(spin + (x - pose.translation()).cross(p))
                        + linear_axis.dot(p);
                    break;
                }
                ancestor = model.links()[index].parent();
            }
        }
    }
    Observables {
        linear,
        angular,
        center: center.scale(1. / total_mass),
        kinetic,
        generalized,
    }
}

fn check_receipt(
    previous: &ArticulatedModel,
    next: &ArticulatedModel,
    before: &FreeFloatingMotionState,
    result: &InertiaReconfiguration,
) {
    let old = physical(previous, before);
    let new = physical(next, &result.state_after);
    assert_eq!(
        result.state_after.base.world_from_base,
        before.base.world_from_base
    );
    assert_eq!(result.state_after.joint_position, before.joint_position);
    vector(new.linear, old.linear);
    vector(new.angular, old.angular);
    vector(result.momentum_before.linear_momentum_world, old.linear);
    vector(result.momentum_before.angular_momentum_world, old.angular);
    vector(result.momentum_after.linear_momentum_world, new.linear);
    vector(result.momentum_after.angular_momentum_world, new.angular);
    close(result.momentum_before.kinetic_energy, old.kinetic);
    close(result.momentum_after.kinetic_energy, new.kinetic);
    close(result.parameter_work_j, new.kinetic - old.kinetic);
    assert!(result.maximum_normalized_momentum_residual <= 1.);
    assert_eq!(
        result.generalized_momentum_before.len(),
        old.generalized.len()
    );
    assert_eq!(
        result.generalized_momentum_after.len(),
        new.generalized.len()
    );
    for i in 0..old.generalized.len() {
        close(new.generalized[i], old.generalized[i]);
        close(result.generalized_momentum_before[i], old.generalized[i]);
        close(result.generalized_momentum_after[i], new.generalized[i]);
    }
}

fn wheel(limits: Option<JointLimits>) -> (ArticulatedModel, FreeFloatingMotionState) {
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
            Se3::from_parts(rotation(Vec3::new(0.2, -0.3, 0.1)), Vec3::new(1., -2., 0.5)).unwrap(),
            Twist::new(Vec3::new(0., 0., 0.4), Vec3::new(0.7, -0.2, 0.3)),
        ),
        joint_position: vec![0.23],
        joint_velocity: vec![1.2],
    };
    (model, state)
}

#[test]
fn coaxial_carrier_and_rotor_keep_both_momenta_and_continue_at_new_rates() {
    let (previous, state) = wheel(None);
    let next = with_inertias(&previous, &[sphere(2., 4.), sphere(1., 3.)]);
    let result =
        transfer_inertia(&previous, &next, &state, controls(&[2e-11], 10.), || false).unwrap();
    // p_q = J_r (omega + qdot) = 1.6 and p_b = J_b omega + p_q = 2.4.
    // Thus omega' = (2.4 - 1.6)/4 and qdot' = 1.6/3 - omega'.
    close(result.generalized_momentum_before[2], 2.4);
    close(result.generalized_momentum_before[6], 1.6);
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 0.2),
    );
    vector(
        result.state_after.base.twist_body.linear,
        state.base.twist_body.linear,
    );
    close(result.state_after.joint_velocity[0], 1. / 3.);
    close(result.parameter_work_j, -14. / 15.);
    check_receipt(&previous, &next, &state, &result);

    // Align the frames and use exact powers of two so the common translation
    // contributes exactly zero momentum change. Its kinetic energy dwarfs the
    // rotor work; subtracting the rounded lab-frame totals would return zero.
    let mut translated = state.clone();
    translated.base.world_from_base =
        Se3::from_parts(So3::identity(), state.base.world_from_base.translation()).unwrap();
    translated.joint_position[0] = 0.;
    translated.base.twist_body.linear = Vec3::new(2_f64.powi(30), -2_f64.powi(29), 2_f64.powi(28));
    let large_translation = transfer_inertia(
        &previous,
        &next,
        &translated,
        controls(&[2e-11], 1.),
        || false,
    )
    .unwrap();
    vector(
        large_translation.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 0.2),
    );
    close(large_translation.state_after.joint_velocity[0], 1. / 3.);
    close(large_translation.parameter_work_j, -14. / 15.);
    assert_eq!(
        large_translation.state_after.base.twist_body.linear,
        translated.base.twist_body.linear
    );
    assert_eq!(
        large_translation.state_after.base.world_from_base,
        translated.base.world_from_base
    );
    assert_eq!(
        large_translation.state_after.joint_position,
        translated.joint_position
    );
    assert_eq!(
        large_translation.momentum_after.linear_momentum_world,
        large_translation.momentum_before.linear_momentum_world
    );
    assert!(large_translation.momentum_before.kinetic_energy > 1e18);
    assert_eq!(
        large_translation.momentum_after.kinetic_energy,
        large_translation.momentum_before.kinetic_energy
    );
    assert!(large_translation.maximum_normalized_momentum_residual <= 1.);

    let dt = 0.18;
    let continued = step_free_floating(
        &next,
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
    close(
        continued.state_after.joint_position[0],
        state.joint_position[0] + dt / 3.,
    );
    close(continued.state_after.joint_velocity[0], 1. / 3.);
    vector(
        continued.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 0.2),
    );
    let world_velocity = state
        .base
        .world_from_base
        .rotation()
        .rotate(state.base.twist_body.linear)
        .unwrap();
    vector(
        continued.state_after.base.world_from_base.translation(),
        state.base.world_from_base.translation() + world_velocity.scale(dt),
    );
    let expected_attitude = state
        .base
        .world_from_base
        .rotation()
        .compose(rotation(Vec3::new(0., 0., 0.2 * dt)))
        .unwrap();
    for (actual, expected) in continued
        .state_after
        .base
        .world_from_base
        .rotation()
        .matrix()
        .m
        .iter()
        .zip(expected_attitude.matrix().m)
    {
        close(*actual, expected);
    }
    close(continued.kinetic_work_residual, 0.);
}

fn serial_tree() -> (ArticulatedModel, ArticulatedModel, FreeFloatingMotionState) {
    let centers = [
        Vec3::new(0.12, -0.21, 0.08),
        Vec3::new(-0.08, 0.17, 0.04),
        Vec3::new(0.11, 0.03, -0.09),
        Vec3::new(0.04, -0.12, 0.06),
    ];
    let masses = [2.3, 1.1, 0.8, 0.4];
    let old_diagonals = [
        Vec3::new(1.2, 1.6, 2.),
        Vec3::new(0.5, 0.7, 0.9),
        Vec3::new(0.3, 0.45, 0.6),
        Vec3::new(0.2, 0.25, 0.3),
    ];
    let new_diagonals = [
        Vec3::new(1.7, 2.2, 2.8),
        Vec3::new(0.8, 1.05, 1.25),
        Vec3::new(0.35, 0.55, 0.7),
        Vec3::new(0.45, 0.55, 0.7),
    ];
    let old: Vec<_> = (0..4)
        .map(|i| {
            inertia(
                masses[i],
                centers[i],
                old_diagonals[i],
                Vec3::new(0.2 + 0.1 * i as f64, -0.15, 0.3),
            )
        })
        .collect();
    let new: Vec<_> = (0..4)
        .map(|i| {
            inertia(
                masses[i],
                centers[i],
                new_diagonals[i],
                Vec3::new(-0.2, 0.35, -0.1 - 0.13 * i as f64),
            )
        })
        .collect();
    let pose = |r, x| Se3::from_parts(rotation(r), x).unwrap();
    let previous = ArticulatedModel::new(vec![
        Link::new(
            "offset carrier",
            None,
            pose(Vec3::new(0.2, -0.1, 0.15), Vec3::new(0.17, -0.04, 0.11)),
            JointModel::FIXED,
            old[0],
        ),
        Link::new(
            "oblique rotor",
            Some(0),
            pose(Vec3::new(0.1, 0.25, -0.15), Vec3::new(0.6, 0.08, -0.11)),
            JointModel::revolute(Vec3::new(1., 2., -1.), None).unwrap(),
            old[1],
        ),
        Link::new(
            "oblique slider",
            Some(1),
            pose(Vec3::new(-0.18, 0.11, 0.2), Vec3::new(0.1, -0.35, 0.27)),
            JointModel::prismatic(Vec3::new(-1., 1., 2.), None).unwrap(),
            old[2],
        ),
        Link::new(
            "screw",
            Some(2),
            pose(Vec3::new(0.11, -0.12, 0.18), Vec3::new(-0.2, 0.17, 0.33)),
            JointModel::helical(Vec3::new(2., -1., 1.), 0.13, None).unwrap(),
            old[3],
        ),
    ])
    .unwrap();
    let next = with_inertias(&previous, &new);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(Vec3::new(0.31, -0.27, 0.19), Vec3::new(2., -1., 0.6)),
            Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(-0.2, 0.6, 0.35)),
        ),
        joint_position: vec![0.28, -0.09, 0.17],
        joint_velocity: vec![0.7, -0.2, 0.5],
    };
    (previous, next, state)
}

#[test]
fn full_tensor_serial_transfer_is_reversible_and_frame_and_boost_covariant() {
    let (previous, next, state) = serial_tree();
    let settings = controls(&[2e-11; 3], 100.);
    let result = transfer_inertia(&previous, &next, &state, settings, || false).unwrap();
    check_receipt(&previous, &next, &state, &result);
    assert!(
        (result.state_after.base.twist_body.angular - state.base.twist_body.angular).norm() > 0.01
    );
    assert!((result.state_after.joint_velocity[0] - state.joint_velocity[0]).abs() > 0.01);

    let reversed =
        transfer_inertia(&next, &previous, &result.state_after, settings, || false).unwrap();
    check_receipt(&next, &previous, &result.state_after, &reversed);
    vector(
        reversed.state_after.base.twist_body.angular,
        state.base.twist_body.angular,
    );
    vector(
        reversed.state_after.base.twist_body.linear,
        state.base.twist_body.linear,
    );
    for (actual, expected) in reversed
        .state_after
        .joint_velocity
        .iter()
        .zip(&state.joint_velocity)
    {
        close(*actual, *expected);
    }
    close(result.parameter_work_j + reversed.parameter_work_j, 0.);

    let scene = Se3::from_parts(
        rotation(Vec3::new(-0.4, 0.2, 0.3)),
        Vec3::new(4., -2.5, 1.7),
    )
    .unwrap();
    let before = physical(&previous, &state);
    let total_mass: f64 = previous
        .links()
        .iter()
        .map(|link| link.inertia().mass())
        .sum();
    for boost in [Vec3::default(), Vec3::new(3., -2., 4.)] {
        let mut transformed = state.clone();
        transformed.base.world_from_base = scene.compose(state.base.world_from_base).unwrap();
        let local_boost = transformed
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
        transformed.base.twist_body.linear = state.base.twist_body.linear + local_boost;
        let changed = transfer_inertia(&previous, &next, &transformed, settings, || false).unwrap();
        check_receipt(&previous, &next, &transformed, &changed);
        vector(
            changed.state_after.base.twist_body.angular,
            result.state_after.base.twist_body.angular,
        );
        vector(
            changed.state_after.base.twist_body.linear,
            result.state_after.base.twist_body.linear + local_boost,
        );
        for (actual, expected) in changed
            .state_after
            .joint_velocity
            .iter()
            .zip(&result.state_after.joint_velocity)
        {
            close(*actual, *expected);
        }
        close(changed.parameter_work_j, result.parameter_work_j);
        let rotated_p = scene.rotation().rotate(before.linear).unwrap();
        let boosted_p = boost.scale(total_mass);
        vector(
            changed.momentum_before.linear_momentum_world,
            rotated_p + boosted_p,
        );
        vector(
            changed.momentum_before.angular_momentum_world,
            scene.rotation().rotate(before.angular).unwrap()
                + scene.translation().cross(rotated_p)
                + scene
                    .transform_point(before.center)
                    .unwrap()
                    .cross(boosted_p),
        );
    }

    // Ordinary native flight must accept the transferred state and preserve its
    // physical total momenta while the noncoaxial joints continue moving.
    let continued = step_free_floating(
        &next,
        &result.state_after,
        FreeFloatingMotionLoads {
            generalized_force: &[0.; 3],
            external_body_wrench: &[Wrench::default(); 4],
            gravity_world: Vec3::default(),
        },
        0.015,
        || false,
    )
    .unwrap();
    let after = physical(&next, &continued.state_after);
    vector(after.linear, before.linear);
    vector(after.angular, before.angular);
    assert_ne!(continued.state_after.joint_position, state.joint_position);
}

#[test]
fn identity_inventory_work_speed_and_cancellation_are_atomic() {
    let (previous, state) = wheel(None);
    let next = with_inertias(&previous, &[sphere(2., 4.), sphere(1., 3.)]);
    let retained = (previous.clone(), next.clone(), state.clone());
    let identity = transfer_inertia(&previous, &previous, &state, controls(&[2e-11], 0.), || {
        false
    })
    .unwrap();
    assert_eq!(identity.state_after, state);
    assert_eq!(identity.parameter_work_j, 0.);
    assert_eq!(
        identity.generalized_momentum_before,
        identity.generalized_momentum_after
    );
    let mut moved_links = previous.links().to_vec();
    let root = &previous.links()[0];
    moved_links[0] = Link::new(
        root.name(),
        root.parent(),
        Se3::from_parts(So3::identity(), Vec3::new(0.01, 0., 0.)).unwrap(),
        root.joint(),
        root.inertia(),
    );
    for changed in [
        with_inertias(&previous, &[sphere(2.1, 2.), sphere(1., 1.)]),
        with_inertias(
            &previous,
            &[
                inertia(
                    2.,
                    Vec3::new(0.01, 0., 0.),
                    Vec3::new(2., 2., 2.),
                    Vec3::default(),
                ),
                sphere(1., 1.),
            ],
        ),
        ArticulatedModel::new(moved_links).unwrap(),
    ] {
        assert!(matches!(
            transfer_inertia(&previous, &changed, &state, controls(&[2e-11], 10.), || {
                false
            }),
            Err(ReconfigurationError::Invalid(_))
        ));
    }
    assert!(matches!(
        transfer_inertia(&previous, &next, &state, controls(&[2e-11], 0.1), || false),
        Err(ReconfigurationError::ParameterWork { .. })
    ));
    let accepted =
        transfer_inertia(&previous, &next, &state, controls(&[2e-11], 10.), || false).unwrap();
    assert!(matches!(
        transfer_inertia(
            &next,
            &previous,
            &accepted.state_after,
            controls(&[2e-11], 0.1),
            || false
        ),
        Err(ReconfigurationError::ParameterWork { .. })
    ));
    assert!(transfer_inertia(&previous, &next, &state, controls(&[], 10.), || false).is_err());

    let (limited, limited_state) = wheel(Some(JointLimits::new(-2., 2., 1.5, 100.).unwrap()));
    let contracted = with_inertias(&limited, &[sphere(2., 2.), sphere(1., 0.25)]);
    // Conserved rotor momentum requires qdot' = 1.6/0.25 - 0.4 = 6,
    // beyond the actual 1.5 rad/s bound. Clamping would destroy momentum.
    assert!(matches!(
        transfer_inertia(
            &limited,
            &contracted,
            &limited_state,
            controls(&[2e-11], 100.),
            || false
        ),
        Err(ReconfigurationError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { link: 1, .. }
        ))
    ));

    let mut polls = 0;
    transfer_inertia(&previous, &next, &state, controls(&[2e-11], 10.), || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 2);
    for cancel_at in [1, polls] {
        let mut calls = 0;
        assert!(matches!(
            transfer_inertia(&previous, &next, &state, controls(&[2e-11], 10.), || {
                calls += 1;
                calls == cancel_at
            }),
            Err(ReconfigurationError::Cancelled)
        ));
    }
    assert_eq!((previous, next, state), retained);
}
