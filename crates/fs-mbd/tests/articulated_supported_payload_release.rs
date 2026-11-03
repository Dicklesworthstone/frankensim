//! Supported pickup/drop cycles preserve the original inventory and velocity.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::inventory::supported::release::{
    SupportedPayloadRelease, SupportedPayloadReleaseControls, release_payload,
};
use fs_mbd::articulated::motion::inventory::supported::{
    SupportedPayloadControls, capture_payload,
};
use fs_mbd::articulated::motion::inventory::{PayloadError, PayloadState};
use fs_mbd::articulated::motion::prescribed::{PrescribedBaseMotionState, PrescribedBasePath};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitSpatialConfig, step_prescribed_base_spatial_implicit,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_kinematics,
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

fn pose(axis: Vec3, origin: Vec3) -> Se3 {
    Se3::from_parts(So3::exp(So3Tangent::new(axis)).unwrap(), origin).unwrap()
}

fn inertia(mass: f64, center: Vec3, diagonal: [f64; 3], axes: Vec3) -> SpatialInertia {
    let r = So3::exp(So3Tangent::new(axes)).unwrap().matrix();
    let tensor = Mat3 {
        m: [
            diagonal[0],
            0.,
            0.,
            0.,
            diagonal[1],
            0.,
            0.,
            0.,
            diagonal[2],
        ],
    };
    SpatialInertia::new(mass, center, r.compose(tensor).compose(r.transpose())).unwrap()
}

fn capture_controls(tolerances: &[f64]) -> SupportedPayloadControls<'_> {
    SupportedPayloadControls {
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_removed_relative_energy_j: 1000.,
    }
}

fn release_controls(tolerances: &[f64]) -> SupportedPayloadReleaseControls<'_> {
    SupportedPayloadReleaseControls {
        angular_momentum_tolerance_nms: 2e-9,
        linear_momentum_tolerance_ns: 2e-9,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
    }
}

fn fixture(
    joint: JointModel,
    base: BaseState,
) -> (ArticulatedModel, PrescribedBaseMotionState, PayloadState) {
    let model = ArticulatedModel::new(vec![Link::new(
        "catcher",
        None,
        Se3::identity(),
        joint,
        inertia(
            3.,
            Vec3::new(0.2, -0.3, 0.1),
            [3., 4., 5.],
            Vec3::new(0.3, -0.2, 0.4),
        ),
    )])
    .unwrap();
    let state = PrescribedBaseMotionState {
        time_s: 4.,
        base,
        joint_position: vec![0.2; model.dof_count()],
        joint_velocity: vec![0.7; model.dof_count()],
    };
    let kin =
        forward_kinematics(&model, base, &state.joint_position, &state.joint_velocity).unwrap();
    let payload = PayloadState {
        inertia: inertia(
            1.5,
            Vec3::new(0.1, -0.2, 0.3),
            [0.4, 0.5, 0.6],
            Vec3::new(-0.3, 0.2, 0.1),
        ),
        world_from_body: kin.world_from_link[0]
            .compose(pose(Vec3::new(0.4, -0.2, 0.3), Vec3::new(0.7, -0.4, 0.2)))
            .unwrap(),
        twist_body: Twist::new(Vec3::new(0.2, -0.4, 0.5), Vec3::new(1.2, -0.5, 0.3)),
    };
    (model, state, payload)
}

// Independent COM/spin accounting, with only native geometry reused.
fn physical(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    payload: Option<PayloadState>,
) -> (Vec3, Vec3, f64) {
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
        .map(|(i, link)| (link.inertia(), kin.world_from_link[i], kin.body_twist[i]))
        .collect();
    if let Some(p) = payload {
        bodies.push((p.inertia, p.world_from_body, p.twist_body));
    }
    let mut linear = Vec3::default();
    let mut angular = Vec3::default();
    let mut energy = 0.;
    for (inertia, pose, twist) in bodies {
        let center = pose.transform_point(inertia.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(inertia.center_of_mass()))
            .unwrap();
        let p = velocity.scale(inertia.mass());
        let spin_body = inertia.inertia_com().apply(twist.angular);
        linear = linear + p;
        angular = angular + pose.rotation().rotate(spin_body).unwrap() + center.cross(p);
        energy +=
            0.5 * inertia.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
    }
    (linear, angular, energy)
}

fn audit(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    release: &SupportedPayloadRelease,
) {
    let before = physical(model, state, None);
    let after = physical(
        &release.model_after,
        &release.state_after,
        Some(release.payload_after),
    );
    let balance = &release.balance;
    vector(after.0, before.0);
    vector(after.1, before.1);
    close(after.2, before.2);
    vector(balance.momentum_before.linear_momentum_world, before.0);
    vector(balance.momentum_before.angular_momentum_world, before.1);
    close(balance.momentum_before.kinetic_energy, before.2);
    vector(balance.momentum_after.linear_momentum_world, after.0);
    vector(balance.momentum_after.angular_momentum_world, after.1);
    close(balance.momentum_after.kinetic_energy, after.2);
    assert_eq!(balance.support_impulse_body, Wrench::default());
    assert_eq!(balance.support_impulse_world, Wrench::default());
    assert_eq!(balance.support_work_j, 0.);
    assert_eq!(balance.accounting_support_work_j, 0.);
    vector(balance.linear_momentum_residual_world, after.0 - before.0);
    vector(balance.angular_momentum_residual_world, after.1 - before.1);
    for value in &balance.joint_momentum_residual {
        close(*value, 0.);
    }
    close(balance.energy_residual_j, 0.);
    close(
        balance.energy_residual_j,
        balance.accounting_kinetic_energy_change_j,
    );
    close(balance.kinetic_energy_change_j, after.2 - before.2);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    assert_eq!(release.state_after, *state);
}

#[test]
fn capture_evolve_and_release_restore_tensor_and_current_rigid_velocity() {
    let joint = JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap();
    let base = BaseState::stationary(pose(Vec3::new(0.3, -0.2, 0.4), Vec3::new(1., -2., 3.)));
    let (model, state, payload) = fixture(joint, base);
    let captured = capture_payload(
        &model,
        &state,
        0,
        payload,
        capture_controls(&[2e-9]),
        || false,
    )
    .unwrap();
    // The loaded COM produces a configuration-dependent gravity torque.
    // Carry therefore consumes the new simultaneous nonlinear midpoint path.
    let carried = step_prescribed_base_spatial_implicit(
        &captured.model_after,
        &captured.state_after,
        PrescribedBasePath::stationary(base.world_from_base),
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.3],
                external_body_wrench: &[Wrench::default()],
                gravity_world: Vec3::new(0., -9.81, 0.),
            },
            world_points: &[],
            connections: &[],
        },
        0.001,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    assert!(carried.scaled_residual_norm <= ImplicitSpatialConfig::default().residual_tolerance);
    let evolved = carried.spatial.motion.state_after;
    assert_ne!(evolved.joint_position, captured.state_after.joint_position);
    let release = release_payload(
        &captured.model_after,
        &evolved,
        &captured.attachment,
        release_controls(&[2e-9]),
        || false,
    )
    .unwrap();
    assert_eq!(release.model_after, model);
    assert_eq!(release.payload_after.inertia, payload.inertia);
    let kin = forward_kinematics(
        &captured.model_after,
        evolved.base,
        &evolved.joint_position,
        &evolved.joint_velocity,
    )
    .unwrap();
    let expected_pose = kin.world_from_link[0]
        .compose(captured.attachment.link_from_payload())
        .unwrap();
    let expected_twist =
        kin.body_twist[0].transform_by(&captured.attachment.link_from_payload().inverse().unwrap());
    assert_eq!(release.payload_after.world_from_body, expected_pose);
    assert_eq!(release.payload_after.twist_body, expected_twist);
    assert_ne!(
        release.payload_after.world_from_body,
        captured.payload_after.world_from_body,
    );
    audit(&captured.model_after, &evolved, &release);
}

#[test]
fn moving_support_drop_inherits_offset_velocity_without_support_work() {
    let pose0 = pose(Vec3::new(0.2, -0.4, 0.3), Vec3::new(2., -1., 3.));
    let velocity = Vec3::new(1.2, -0.6, 0.4);
    let spin = Vec3::new(0.3, -0.2, 0.4);
    let (model, state, payload) = fixture(
        JointModel::FIXED,
        BaseState::prescribed(pose0, Twist::new(spin, velocity), Twist::zero()),
    );
    let capture =
        capture_payload(&model, &state, 0, payload, capture_controls(&[]), || false).unwrap();
    let release = release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        release_controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(release.model_after, model);
    let r = capture.attachment.link_from_payload();
    let expected_origin_velocity_in_link = velocity + spin.cross(r.translation());
    vector(
        r.rotation()
            .rotate(release.payload_after.twist_body.linear)
            .unwrap(),
        expected_origin_velocity_in_link,
    );
    vector(
        r.rotation()
            .rotate(release.payload_after.twist_body.angular)
            .unwrap(),
        spin,
    );
    audit(&capture.model_after, &capture.state_after, &release);
}

#[test]
fn nested_inventory_releases_in_reverse_order_and_repeated_release_refuses() {
    let (model, state, payload) =
        fixture(JointModel::FIXED, BaseState::stationary(Se3::identity()));
    let first =
        capture_payload(&model, &state, 0, payload, capture_controls(&[]), || false).unwrap();
    let second_payload = PayloadState {
        world_from_body: pose(Vec3::new(-0.3, 0.2, 0.1), Vec3::new(-0.4, 0.6, -0.2)),
        twist_body: Twist::zero(),
        ..payload
    };
    let second = capture_payload(
        &first.model_after,
        &first.state_after,
        0,
        second_payload,
        capture_controls(&[]),
        || false,
    )
    .unwrap();
    assert!(matches!(
        release_payload(
            &second.model_after,
            &second.state_after,
            &first.attachment,
            release_controls(&[]),
            || false,
        ),
        Err(PayloadError::Invalid(
            "release attachment inventory changed"
        )),
    ));
    let drop_second = release_payload(
        &second.model_after,
        &second.state_after,
        &second.attachment,
        release_controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(drop_second.model_after, first.model_after);
    let drop_first = release_payload(
        &drop_second.model_after,
        &drop_second.state_after,
        &first.attachment,
        release_controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(drop_first.model_after, model);
    assert_eq!(drop_first.state_after, first.state_after);
    assert!(matches!(
        release_payload(
            &drop_first.model_after,
            &drop_first.state_after,
            &first.attachment,
            release_controls(&[]),
            || false,
        ),
        Err(PayloadError::Invalid(
            "release attachment inventory changed"
        )),
    ));
    audit(&second.model_after, &second.state_after, &drop_second);
    audit(
        &drop_second.model_after,
        &drop_second.state_after,
        &drop_first,
    );
}

#[test]
fn cancellation_and_invalid_controls_never_spend_or_modify_attachment() {
    let (model, state, payload) =
        fixture(JointModel::FIXED, BaseState::stationary(Se3::identity()));
    let capture =
        capture_payload(&model, &state, 0, payload, capture_controls(&[]), || false).unwrap();
    let original_model = capture.model_after.clone();
    let original_state = capture.state_after.clone();
    let attachment = capture.attachment.clone();
    let calls = Cell::new(0);
    release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        release_controls(&[]),
        || {
            calls.set(calls.get() + 1);
            false
        },
    )
    .unwrap();
    for stop in [1, calls.get()] {
        let seen = Cell::new(0);
        assert_eq!(
            release_payload(
                &capture.model_after,
                &capture.state_after,
                &capture.attachment,
                release_controls(&[]),
                || {
                    seen.set(seen.get() + 1);
                    seen.get() >= stop
                },
            )
            .unwrap_err(),
            PayloadError::Cancelled,
        );
    }
    let invalid = SupportedPayloadReleaseControls {
        energy_tolerance_j: -1.,
        ..release_controls(&[])
    };
    assert!(
        release_payload(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            invalid,
            || false,
        )
        .is_err()
    );
    assert_eq!(capture.model_after, original_model);
    assert_eq!(capture.state_after, original_state);
    assert_eq!(capture.attachment, attachment);
    let released = release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        release_controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(released.model_after, model);
}
