//! G0/G1/G3/G4: finite payload ejection through original supported dynamics.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::motion::inventory::release::SeparationImpulse;
use fs_mbd::articulated::motion::inventory::supported::release::{
    SupportedPayloadRelease, SupportedPayloadReleaseControls, SupportedPayloadSeparationControls,
    release_payload, release_payload_with_impulse,
};
use fs_mbd::articulated::motion::inventory::supported::{
    SupportedPayloadCapture, SupportedPayloadControls, capture_payload,
};
use fs_mbd::articulated::motion::inventory::{PayloadError, PayloadState};
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
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
fn pose(axis: Vec3, origin: Vec3) -> Se3 {
    Se3::from_parts(So3::exp(So3Tangent::new(axis)).unwrap(), origin).unwrap()
}
fn inertia(mass: f64, center: Vec3, diagonal: [f64; 3], axes: Vec3) -> SpatialInertia {
    let r = So3::exp(So3Tangent::new(axes)).unwrap().matrix();
    SpatialInertia::new(
        mass,
        center,
        r.compose(Mat3 {
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
        })
        .compose(r.transpose()),
    )
    .unwrap()
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
fn controls(tolerances: &[f64]) -> SupportedPayloadSeparationControls<'_> {
    SupportedPayloadSeparationControls {
        balance: SupportedPayloadReleaseControls {
            angular_momentum_tolerance_nms: 2e-9,
            linear_momentum_tolerance_ns: 2e-9,
            joint_momentum_tolerances: tolerances,
            energy_tolerance_j: 2e-8,
        },
        maximum_separation_work_j: 1000.,
    }
}
fn simple(joint: JointModel, speed: f64, rate: f64) -> (ArticulatedModel, SupportedPayloadCapture) {
    let model = ArticulatedModel::new(vec![Link::new(
        "launcher",
        None,
        Se3::identity(),
        joint,
        inertia(3., Vec3::default(), [2., 2., 2.], Vec3::default()),
    )])
    .unwrap();
    let n = model.dof_count();
    let state = PrescribedBaseMotionState {
        time_s: 4.,
        base: BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(speed, 0., 0.)),
            Twist::zero(),
        ),
        joint_position: vec![0.4; n],
        joint_velocity: vec![rate; n],
    };
    let kin = forward_kinematics(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let payload = PayloadState {
        inertia: inertia(2., Vec3::default(), [1., 1., 1.], Vec3::default()),
        world_from_body: kin.world_from_link[0],
        twist_body: kin.body_twist[0],
    };
    let capture = capture_payload(
        &model,
        &state,
        0,
        payload,
        capture_controls(&vec![2e-9; n]),
        || false,
    )
    .unwrap();
    (model, capture)
}
fn axial_impulse() -> SeparationImpulse {
    SeparationImpulse {
        point_link_m: Vec3::default(),
        impulse_on_payload_link_ns: Vec3::new(4., 0., 0.),
        couple_on_payload_link_nms: Vec3::default(),
    }
}

struct Physical {
    linear: Vec3,
    angular: Vec3,
    kinetic: f64,
    joint: Vec<f64>,
}

// Native FK supplies geometry only. COM/spin formulas and each ancestor's
// world screw independently reconstruct physical and joint momentum.
fn physical(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    extra: Option<(usize, PayloadState)>,
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
    if let Some((target, payload)) = extra {
        bodies.push((
            target,
            payload.inertia,
            payload.world_from_body,
            payload.twist_body,
        ));
    }
    let mut result = Physical {
        linear: Vec3::default(),
        angular: Vec3::default(),
        kinetic: 0.,
        joint: vec![0.; model.dof_count()],
    };
    for (target, inertia, pose, twist) in bodies {
        let center = pose.transform_point(inertia.center_of_mass()).unwrap();
        let velocity = pose
            .rotation()
            .rotate(twist.linear + twist.angular.cross(inertia.center_of_mass()))
            .unwrap();
        let linear = velocity.scale(inertia.mass());
        let spin_body = inertia.inertia_com().apply(twist.angular);
        let spin = pose.rotation().rotate(spin_body).unwrap();
        result.linear = result.linear + linear;
        result.angular = result.angular + spin + center.cross(linear);
        result.kinetic +=
            0.5 * inertia.mass() * velocity.dot(velocity) + 0.5 * twist.angular.dot(spin_body);
        let mut ancestor = Some(target);
        while let Some(index) = ancestor {
            if let Some(dof) = model.dof_index(index) {
                let screw = model.links()[index].joint().motion_subspace();
                let joint_pose = kin.world_from_link[index];
                let angular_axis = joint_pose.rotation().rotate(screw.angular).unwrap();
                let linear_axis = joint_pose.rotation().rotate(screw.linear).unwrap();
                result.joint[dof] += angular_axis
                    .dot(spin + (center - joint_pose.translation()).cross(linear))
                    + linear_axis.dot(linear);
            }
            ancestor = model.links()[index].parent();
        }
    }
    result
}
fn audit(
    loaded: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    target: usize,
    release: &SupportedPayloadRelease,
) {
    let before = physical(loaded, state, None);
    let after = physical(
        &release.model_after,
        &release.state_after,
        Some((target, release.payload_after)),
    );
    let b = &release.balance;
    vector(b.momentum_before.linear_momentum_world, before.linear);
    vector(b.momentum_before.angular_momentum_world, before.angular);
    close(b.momentum_before.kinetic_energy, before.kinetic);
    vector(b.momentum_after.linear_momentum_world, after.linear);
    vector(b.momentum_after.angular_momentum_world, after.angular);
    close(b.momentum_after.kinetic_energy, after.kinetic);
    vector(b.linear_momentum_change_world, after.linear - before.linear);
    vector(
        b.angular_momentum_change_world,
        after.angular - before.angular,
    );
    vector(b.support_impulse_world.force, after.linear - before.linear);
    vector(
        b.support_impulse_world.torque,
        after.angular - before.angular,
    );
    vector(b.linear_momentum_residual_world, Vec3::default());
    vector(b.angular_momentum_residual_world, Vec3::default());
    for ((old, new), residual) in before
        .joint
        .iter()
        .zip(&after.joint)
        .zip(&b.joint_momentum_residual)
    {
        close(*old, *new);
        close(*residual, new - old);
    }
    close(b.kinetic_energy_change_j, after.kinetic - before.kinetic);
    close(
        b.kinetic_energy_change_j,
        b.support_work_j + release.supplied_work_j,
    );
    close(
        b.support_work_j,
        b.support_impulse_body.pairing(state.base.twist_body),
    );
    close(
        b.accounting_support_work_j,
        b.support_work_j - b.accounting_boost_world.dot(b.support_impulse_world.force),
    );
    close(
        b.accounting_kinetic_energy_change_j,
        b.kinetic_energy_change_j - b.accounting_boost_world.dot(b.linear_momentum_change_world),
    );
    close(b.energy_residual_j, 0.);
    assert!(b.maximum_normalized_momentum_residual <= 1.);
    assert_eq!(release.state_after.time_s, state.time_s);
    assert_eq!(release.state_after.base, state.base);
    assert_eq!(release.state_after.joint_position, state.joint_position);
}

#[test]
fn grounded_launcher_separates_internal_source_work_from_signed_support_work() {
    for speed in [0., 2., -2.] {
        let (retained, capture) = simple(JointModel::FIXED, speed, 0.);
        let result = release_payload_with_impulse(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            axial_impulse(),
            controls(&[]),
            || false,
        )
        .unwrap();
        assert_eq!(result.model_after, retained);
        assert_eq!(result.state_after, capture.state_after);
        vector(
            result.payload_after.twist_body.linear,
            Vec3::new(speed + 2., 0., 0.),
        );
        vector(
            result.balance.support_impulse_body.force,
            Vec3::new(4., 0., 0.),
        );
        vector(result.balance.support_impulse_body.torque, Vec3::default());
        close(result.supplied_work_j, 4.);
        close(result.balance.support_work_j, 4. * speed);
        close(result.balance.kinetic_energy_change_j, 4. + 4. * speed);
        audit(&capture.model_after, &capture.state_after, 0, &result);
    }
}

#[test]
fn finite_sliding_launcher_recoils_without_support_or_effort_limit_substitution() {
    let joint = JointModel::prismatic(
        Vec3::new(1., 0., 0.),
        Some(JointLimits::new(-2., 2., 10., 1e-8).unwrap()),
    )
    .unwrap();
    let (retained, capture) = simple(joint, 2., 1.);
    let result = release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        axial_impulse(),
        controls(&[2e-9]),
        || false,
    )
    .unwrap();
    assert_eq!(result.model_after, retained);
    close(result.state_after.joint_velocity[0], -1. / 3.);
    close(result.payload_after.twist_body.linear.x, 5.);
    vector(result.balance.support_impulse_world.force, Vec3::default());
    vector(result.balance.support_impulse_world.torque, Vec3::default());
    close(result.supplied_work_j, 20. / 3.);
    close(result.balance.support_work_j, 0.);
    audit(&capture.model_after, &capture.state_after, 0, &result);
}

fn branched() -> SupportedPayloadCapture {
    let sphere = |mass| inertia(mass, Vec3::default(), [2., 2., 2.], Vec3::default());
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
            inertia(
                2.,
                Vec3::new(-0.2, 0.1, 0.3),
                [2., 3., 4.],
                Vec3::new(0.3, -0.2, 0.4),
            ),
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
    let state = PrescribedBaseMotionState {
        time_s: 4.,
        base: BaseState::prescribed(
            pose(Vec3::new(0.3, -0.4, 0.2), Vec3::new(2., -3., 1.)),
            Twist::new(Vec3::new(0.2, -0.3, 0.4), Vec3::new(1.1, -0.6, 0.2)),
            Twist::zero(),
        ),
        joint_position: vec![0.3, -0.4, 0.2],
        joint_velocity: vec![0.7, -0.5, 0.9],
    };
    let kin = forward_kinematics(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let payload = PayloadState {
        inertia: inertia(
            1.5,
            Vec3::new(0.1, -0.2, 0.3),
            [0.4, 0.5, 0.6],
            Vec3::new(0.2, -0.3, 0.4),
        ),
        world_from_body: kin.world_from_link[1]
            .compose(pose(Vec3::new(0.4, -0.2, 0.1), Vec3::new(0.7, -0.2, 0.3)))
            .unwrap(),
        twist_body: Twist::new(Vec3::new(-0.3, 0.6, 0.2), Vec3::new(0.3, -0.8, 1.2)),
    };
    capture_payload(
        &model,
        &state,
        1,
        payload,
        capture_controls(&[2e-9; 3]),
        || false,
    )
    .unwrap()
}
fn spatial_impulse() -> SeparationImpulse {
    SeparationImpulse {
        point_link_m: Vec3::new(0.2, -0.1, 0.3),
        impulse_on_payload_link_ns: Vec3::new(0.6, -0.4, 1.2),
        couple_on_payload_link_nms: Vec3::new(0.3, 0.2, -0.1),
    }
}

// Explicit 3x3 adjugate gives an independent payload COM/angular response,
// unlike the production six-component spatial-inertia Cholesky solve.
fn inverse_apply(a: Mat3, b: Vec3) -> Vec3 {
    let m = a.m;
    let cof = [
        m[4] * m[8] - m[5] * m[7],
        m[2] * m[7] - m[1] * m[8],
        m[1] * m[5] - m[2] * m[4],
        m[5] * m[6] - m[3] * m[8],
        m[0] * m[8] - m[2] * m[6],
        m[2] * m[3] - m[0] * m[5],
        m[3] * m[7] - m[4] * m[6],
        m[1] * m[6] - m[0] * m[7],
        m[0] * m[4] - m[1] * m[3],
    ];
    let determinant = m[0] * cof[0] + m[1] * cof[3] + m[2] * cof[6];
    Mat3 { m: cof }.apply(b).scale(1. / determinant)
}
fn world_point_velocity(pose: Se3, twist: Twist, world_point: Vec3) -> Vec3 {
    let point = pose
        .inverse()
        .unwrap()
        .transform_point(world_point)
        .unwrap();
    pose.rotation()
        .rotate(twist.linear + twist.angular.cross(point))
        .unwrap()
}

#[test]
fn offset_full_tensor_impulses_share_branched_recoil_and_physical_pair_work() {
    let capture = branched();
    let impulse = spatial_impulse();
    let result = release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        impulse,
        controls(&[2e-9; 3]),
        || false,
    )
    .unwrap();
    let r = capture.attachment.link_from_payload();
    let force_body = r
        .rotation()
        .inverse()
        .rotate(impulse.impulse_on_payload_link_ns)
        .unwrap();
    let point_body = r
        .inverse()
        .unwrap()
        .transform_point(impulse.point_link_m)
        .unwrap();
    let free_body = r
        .rotation()
        .inverse()
        .rotate(impulse.couple_on_payload_link_nms)
        .unwrap();
    let center = capture.payload_after.inertia.center_of_mass();
    let torque_com = (point_body - center).cross(force_body) + free_body;
    let angular_delta = inverse_apply(capture.payload_after.inertia.inertia_com(), torque_com);
    let linear_delta =
        force_body.scale(1. / capture.payload_after.inertia.mass()) - angular_delta.cross(center);
    let actual_delta = result
        .payload_after
        .twist_body
        .minus(capture.payload_after.twist_body);
    vector(actual_delta.angular, angular_delta);
    vector(actual_delta.linear, linear_delta);
    vector(result.impulse_on_payload_body.force, force_body);
    vector(
        result.impulse_on_payload_body.torque,
        point_body.cross(force_body) + free_body,
    );
    vector(
        result.impulse_on_retained_link.force,
        impulse.impulse_on_payload_link_ns.scale(-1.),
    );
    vector(
        result.impulse_on_retained_link.torque,
        (impulse
            .point_link_m
            .cross(impulse.impulse_on_payload_link_ns)
            + impulse.couple_on_payload_link_nms)
            .scale(-1.),
    );
    let before = forward_kinematics(
        &capture.model_after,
        capture.state_after.base,
        &capture.state_after.joint_position,
        &capture.state_after.joint_velocity,
    )
    .unwrap();
    let after = forward_kinematics(
        &result.model_after,
        result.state_after.base,
        &result.state_after.joint_position,
        &result.state_after.joint_velocity,
    )
    .unwrap();
    let link_pose = before.world_from_link[1];
    let point_world = link_pose.transform_point(impulse.point_link_m).unwrap();
    let payload_pose = result.payload_after.world_from_body;
    let relative_midpoint =
        (world_point_velocity(payload_pose, capture.payload_after.twist_body, point_world)
            + world_point_velocity(payload_pose, result.payload_after.twist_body, point_world)
            - world_point_velocity(link_pose, before.body_twist[1], point_world)
            - world_point_velocity(link_pose, after.body_twist[1], point_world))
        .scale(0.5);
    let relative_spin = (payload_pose
        .rotation()
        .rotate(capture.payload_after.twist_body.angular + result.payload_after.twist_body.angular)
        .unwrap()
        - link_pose
            .rotation()
            .rotate(before.body_twist[1].angular + after.body_twist[1].angular)
            .unwrap())
    .scale(0.5);
    let force_world = link_pose
        .rotation()
        .rotate(impulse.impulse_on_payload_link_ns)
        .unwrap();
    let free_world = link_pose
        .rotation()
        .rotate(impulse.couple_on_payload_link_nms)
        .unwrap();
    close(
        result.supplied_work_j,
        force_world.dot(relative_midpoint) + free_world.dot(relative_spin),
    );
    assert!(result.supplied_work_j > 0.);
    assert!(result.balance.support_impulse_world.torque.norm() > 1e-3);
    for i in 0..3 {
        assert_ne!(
            result.state_after.joint_velocity[i],
            capture.state_after.joint_velocity[i]
        );
    }
    audit(&capture.model_after, &capture.state_after, 1, &result);
}

#[test]
fn rotating_world_frame_and_inertial_boost_preserve_separation_and_transform_support_work() {
    let capture = branched();
    let impulse = spatial_impulse();
    let result = release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        impulse,
        controls(&[2e-9; 3]),
        || false,
    )
    .unwrap();
    let world = pose(Vec3::new(-0.4, 0.2, 0.5), Vec3::new(3., -2., 1.));
    let boost = Vec3::new(17., -9., 6.);
    let mut moved = capture.state_after.clone();
    moved.base.world_from_base = world.compose(moved.base.world_from_base).unwrap();
    let local_boost = moved
        .base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(boost)
        .unwrap();
    moved.base.twist_body.linear = moved.base.twist_body.linear + local_boost;
    moved.base.acceleration_body.linear =
        moved.base.acceleration_body.linear - moved.base.twist_body.angular.cross(local_boost);
    let transformed = release_payload_with_impulse(
        &capture.model_after,
        &moved,
        &capture.attachment,
        impulse,
        controls(&[2e-9; 3]),
        || false,
    )
    .unwrap();
    assert_eq!(
        transformed.state_after.joint_velocity,
        result.state_after.joint_velocity
    );
    assert_eq!(
        transformed.balance.support_impulse_body,
        result.balance.support_impulse_body
    );
    let expected = result.balance.support_impulse_world.transform_by(&world);
    vector(
        transformed.balance.support_impulse_world.force,
        expected.force,
    );
    vector(
        transformed.balance.support_impulse_world.torque,
        expected.torque,
    );
    close(transformed.supplied_work_j, result.supplied_work_j);
    let source_shift = boost.dot(expected.force);
    close(
        transformed.balance.support_work_j,
        result.balance.support_work_j + source_shift,
    );
    close(
        transformed.balance.kinetic_energy_change_j,
        result.balance.kinetic_energy_change_j + source_shift,
    );
    close(
        transformed.balance.kinetic_before_j,
        result.balance.kinetic_before_j,
    );
    close(
        transformed.balance.kinetic_after_j,
        result.balance.kinetic_after_j,
    );
    let expected_pose = world.compose(result.payload_after.world_from_body).unwrap();
    vector(
        transformed.payload_after.world_from_body.translation(),
        expected_pose.translation(),
    );
    vector(
        transformed.payload_after.twist_body.angular,
        result.payload_after.twist_body.angular,
    );
    vector(
        transformed.payload_after.twist_body.linear,
        result.payload_after.twist_body.linear
            + expected_pose.rotation().inverse().rotate(boost).unwrap(),
    );
    audit(&capture.model_after, &moved, 1, &transformed);
}

#[test]
fn zero_impulse_matches_passive_accounts_and_work_budget_or_cancellation_is_atomic() {
    let (retained, capture) = simple(JointModel::FIXED, 2., 0.);
    let saved_model = capture.model_after.clone();
    let saved_state = capture.state_after.clone();
    let saved_attachment = capture.attachment.clone();
    let c = controls(&[]);
    let passive = release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        c.balance,
        || false,
    )
    .unwrap();
    let zero = release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        SeparationImpulse {
            point_link_m: Vec3::new(0.2, -0.3, 0.1),
            ..Default::default()
        },
        SupportedPayloadSeparationControls {
            maximum_separation_work_j: 0.,
            ..c
        },
        || false,
    )
    .unwrap();
    assert_eq!(zero.model_after, retained);
    assert_eq!(zero.state_after, passive.state_after);
    assert_eq!(zero.payload_after, passive.payload_after);
    assert_eq!(zero.balance, passive.balance);
    assert_eq!(zero.supplied_work_j, 0.);
    let polls = Cell::new(0);
    release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        axial_impulse(),
        c,
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    for stop in [1, polls.get()] {
        let seen = Cell::new(0);
        assert_eq!(
            release_payload_with_impulse(
                &capture.model_after,
                &capture.state_after,
                &capture.attachment,
                axial_impulse(),
                c,
                || {
                    seen.set(seen.get() + 1);
                    seen.get() >= stop
                },
            )
            .unwrap_err(),
            PayloadError::Cancelled
        );
    }
    assert!(matches!(
        release_payload_with_impulse(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            axial_impulse(),
            SupportedPayloadSeparationControls {
                maximum_separation_work_j: 3.9,
                ..c
            },
            || false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert!(
        release_payload_with_impulse(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            axial_impulse(),
            SupportedPayloadSeparationControls {
                maximum_separation_work_j: f64::NAN,
                ..c
            },
            || false,
        )
        .is_err()
    );
    assert!(
        release_payload_with_impulse(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            SeparationImpulse {
                couple_on_payload_link_nms: Vec3::new(f64::INFINITY, 0., 0.),
                ..axial_impulse()
            },
            c,
            || false,
        )
        .is_err()
    );
    assert_eq!(capture.model_after, saved_model);
    assert_eq!(capture.state_after, saved_state);
    assert_eq!(capture.attachment, saved_attachment);
}

#[test]
fn actual_native_joint_speed_limits_still_refuse_ejection() {
    let joint = JointModel::prismatic(
        Vec3::new(1., 0., 0.),
        Some(JointLimits::new(-2., 2., 0.5, 1e-8).unwrap()),
    )
    .unwrap();
    let (_, capture) = simple(joint, 0., 0.);
    let before = capture.state_after.clone();
    assert!(matches!(
        release_payload_with_impulse(
            &capture.model_after,
            &capture.state_after,
            &capture.attachment,
            axial_impulse(),
            controls(&[2e-9]),
            || false,
        ),
        Err(PayloadError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    assert_eq!(capture.state_after, before);
}
