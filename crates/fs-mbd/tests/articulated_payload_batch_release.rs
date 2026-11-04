//! G0/G1/G3/G4: simultaneous retained inertia, shared recoil and signed pair work.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::inventory::bilateral;
use fs_mbd::articulated::motion::inventory::bilateral::release::release_payloads as release_closed;
use fs_mbd::articulated::motion::inventory::release::{
    PayloadBatchRelease, PayloadReleaseInput, SeparationImpulse, release_payload, release_payloads,
};
use fs_mbd::articulated::motion::inventory::{
    PayloadBatchCapture, PayloadCaptureInput, PayloadControls, PayloadError, PayloadState,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::MaterialPoint;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{BilateralLoop, DistanceLoop};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, forward_kinematics,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-9 * (1.0 + expected.abs()),
        "{actual:.16e} != {expected:.16e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn controls(tolerances: &[f64]) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-8,
        linear_momentum_tolerance_ns: 2e-8,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_energy_exchange_j: 100.,
        maximum_constraint_work_j: 1e-8,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn sphere(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [mass, 0., 0., 0., mass, 0., 0., 0., mass],
        },
    )
    .unwrap()
}
fn single_body(inertia: SpatialInertia) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}
fn translating() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    }
}
fn kick(value: f64) -> SeparationImpulse {
    SeparationImpulse {
        impulse_on_payload_link_ns: Vec3::new(value, 0., 0.),
        ..SeparationImpulse::default()
    }
}
fn attach(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    bodies: &[(usize, SpatialInertia, Se3)],
    loops: &[BilateralLoop],
    rows: &[JointConstraint],
) -> PayloadBatchCapture {
    let kin = forward_kinematics(
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
    let inputs: Vec<_> = bodies
        .iter()
        .map(|&(target_link, inertia, local)| PayloadCaptureInput {
            target_link,
            payload: PayloadState {
                inertia,
                world_from_body: kin.world_from_link[target_link].compose(local).unwrap(),
                twist_body: kin.body_twist[target_link].transform_by(&local.inverse().unwrap()),
            },
        })
        .collect();
    let tolerances = vec![2e-8; model.dof_count()];
    let result = bilateral::capture_payloads(
        model,
        state,
        &inputs,
        loops,
        rows,
        controls(&tolerances),
        Default::default(),
        || false,
    )
    .unwrap()
    .capture;
    close(result.removed_relative_energy_j, 0.);
    result
}
fn balance(result: &PayloadBatchRelease) {
    assert!(result.balance.maximum_normalized_momentum_residual <= 1.);
    vector(
        result.balance.momentum_after.linear_momentum_world,
        result.balance.momentum_before.linear_momentum_world,
    );
    vector(
        result.balance.momentum_after.angular_momentum_world,
        result.balance.momentum_before.angular_momentum_world,
    );
    close(result.balance.energy_residual_j, 0.);
    close(
        result.balance.kinetic_after_j - result.balance.kinetic_before_j,
        result.supplied_work_j + result.constraint_work_j,
    );
}

#[test]
fn simultaneous_sources_share_one_retained_recoil_and_keep_cross_work() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(4.),
        ),
        Link::new(
            "left",
            Some(0),
            Se3::identity(),
            JointModel::FIXED,
            sphere(1.),
        ),
        Link::new(
            "right",
            Some(0),
            Se3::identity(),
            JointModel::FIXED,
            sphere(1.),
        ),
    ])
    .unwrap();
    let loaded = attach(
        &model,
        &translating(),
        &[
            (1, sphere(2.), Se3::identity()),
            (2, sphere(3.), Se3::identity()),
        ],
        &[],
        &[],
    );
    let inputs = [
        PayloadReleaseInput {
            attachment: &loaded.attachments[0],
            impulse: Some(kick(1.2)),
        },
        PayloadReleaseInput {
            attachment: &loaded.attachments[1],
            impulse: Some(kick(0.8)),
        },
    ];
    let result = release_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    let v = loaded.state_after.base.twist_body.linear.x;
    assert_eq!(result.model_after, model);
    close(result.state_after.base.twist_body.linear.x, v - 2. / 6.);
    close(
        result.payloads_after[0].payload_after.twist_body.linear.x,
        v + 1.2 / 2.,
    );
    close(
        result.payloads_after[1].payload_after.twist_body.linear.x,
        v + 0.8 / 3.,
    );
    // 1/2 [(J1+J2)^2/M + J1^2/m1 + J2^2/m2] includes J1*J2/M = 0.16 J.
    close(result.supplied_work_j, 0.8);
    close(result.payloads_after[0].supplied_work_j, 0.56);
    close(result.payloads_after[1].supplied_work_j, 0.24);
    close(result.constraint_work_j, 0.);
    balance(&result);

    let permuted = release_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &[inputs[1], inputs[0]],
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(permuted.state_after, result.state_after);
    assert_eq!(permuted.payloads_after[1], result.payloads_after[0]);
    assert_eq!(permuted.payloads_after[0], result.payloads_after[1]);
    close(permuted.supplied_work_j, result.supplied_work_j);
    let delegated = release_closed(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &[],
        &[],
        controls(&[]),
        Default::default(),
        || false,
    )
    .unwrap();
    assert!(delegated.reaction.is_none());
    assert_eq!(delegated.release.state_after, result.state_after);
    assert_eq!(delegated.release.payloads_after, result.payloads_after);
    assert_eq!(delegated.release.balance, result.balance);

    // Two successive releases are a different experiment: the first recoil
    // still carries the second payload. It must not implement this batch.
    let first = release_payload(
        &loaded.model_after,
        &loaded.state_after,
        &loaded.attachments[0],
        Some(kick(1.2)),
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    let second = release_payload(
        &first.model_after,
        &first.state_after,
        &loaded.attachments[1],
        Some(kick(0.8)),
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    close(
        second.state_after.base.twist_body.linear.x,
        v - 1.2 / 9. - 0.8 / 6.,
    );
    assert!(
        (second.state_after.base.twist_body.linear.x - result.state_after.base.twist_body.linear.x)
            .abs()
            > 0.05
    );
}

#[test]
fn repeated_link_inventory_restores_exactly_and_late_refusals_are_atomic() {
    let model = single_body(sphere(5.));
    let loaded = attach(
        &model,
        &translating(),
        &[
            (0, sphere(2.), Se3::identity()),
            (0, sphere(3.), Se3::identity()),
        ],
        &[],
        &[],
    );
    let mut state = loaded.state_after.clone();
    state.base.twist_body.angular.z = -0.0;
    let inputs = [
        PayloadReleaseInput {
            attachment: &loaded.attachments[1],
            impulse: None,
        },
        PayloadReleaseInput {
            attachment: &loaded.attachments[0],
            impulse: Some(kick(0.)),
        },
    ];
    let saved = (
        loaded.model_after.clone(),
        state.clone(),
        loaded.attachments.clone(),
    );
    let polls = Cell::new(0usize);
    let result = release_payloads(
        &loaded.model_after,
        &state,
        &inputs,
        &[],
        controls(&[]),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    assert_eq!(result.model_after, model);
    assert_eq!(result.state_after, state);
    for (actual, original) in result
        .state_after
        .base
        .twist_body
        .to_array()
        .into_iter()
        .zip(state.base.twist_body.to_array())
    {
        assert_eq!(actual.to_bits(), original.to_bits());
    }
    assert_eq!(
        result.payloads_after[0].payload_after,
        loaded.payloads_after[1]
    );
    assert_eq!(
        result.payloads_after[1].payload_after,
        loaded.payloads_after[0]
    );
    assert_eq!(result.supplied_work_j, 0.);
    assert_eq!(result.constraint_work_j, 0.);
    balance(&result);
    for invalid in [
        &[inputs[1], inputs[0]][..],
        &[inputs[0], inputs[0]][..],
        &[][..],
    ] {
        assert!(matches!(
            release_payloads(
                &loaded.model_after,
                &state,
                invalid,
                &[],
                controls(&[]),
                || false,
            ),
            Err(PayloadError::Invalid(_))
        ));
    }
    let total = polls.get();
    assert!(total > 4);
    for threshold in [1, total / 2, total] {
        let calls = Cell::new(0usize);
        assert!(matches!(
            release_payloads(
                &loaded.model_after,
                &state,
                &inputs,
                &[],
                controls(&[]),
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= threshold
                },
            ),
            Err(PayloadError::Cancelled)
        ));
    }
    assert_eq!(
        (
            loaded.model_after.clone(),
            state.clone(),
            loaded.attachments.clone()
        ),
        saved
    );
    let retry = release_payloads(
        &loaded.model_after,
        &state,
        &inputs,
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(retry.state_after, result.state_after);
    assert_eq!(retry.payloads_after, result.payloads_after);
    assert_eq!(retry.balance, result.balance);
}

fn sliders() -> ArticulatedModel {
    let mut links = vec![Link::new(
        "carrier",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(5.),
    )];
    for (i, mass) in [2., 3., 4.].into_iter().enumerate() {
        links.push(Link::new(
            format!("slider {i}"),
            Some(0),
            Se3::identity(),
            JointModel::prismatic(
                Vec3::new(1., 0., 0.),
                Some(JointLimits::new(-10., 10., 1.2, 0.1).unwrap()),
            )
            .unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn slider_state() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
        ..translating()
    }
}
fn closure() -> [BilateralLoop; 1] {
    [BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        length_m: 2.,
    })]
}
fn gear() -> [JointConstraint; 1] {
    [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 2e-3,
    }]
}

#[test]
fn mixed_closed_ejection_uses_the_combined_endpoint_and_signed_constraint_work() {
    let model = sliders();
    let rows = gear();
    let loops = closure();
    let loaded = attach(
        &model,
        &slider_state(),
        &[
            (1, sphere(2.), Se3::identity()),
            (2, sphere(3.), Se3::identity()),
        ],
        &loops,
        &rows,
    );
    let inputs = [
        PayloadReleaseInput {
            attachment: &loaded.attachments[0],
            impulse: Some(kick(-4.)),
        },
        PayloadReleaseInput {
            attachment: &loaded.attachments[1],
            impulse: Some(kick(1.)),
        },
    ];
    let tolerances = [2e-8; 3];
    let controls = PayloadControls {
        // The signed total is below 5 J even though the first pair alone is
        // above 5 J and the sum of absolute pair works is larger still.
        maximum_energy_exchange_j: 5.,
        maximum_constraint_work_j: 1e-2,
        ..controls(&tolerances)
    };
    for defect in [0., 1e-3] {
        let mut state = loaded.state_after.clone();
        state.joint_velocity[2] -= 0.5 * defect;
        let result = release_closed(
            &loaded.model_after,
            &state,
            &inputs,
            &loops,
            &rows,
            controls,
            Default::default(),
            || false,
        )
        .unwrap();
        let released = &result.release;
        let event = result.reaction.as_ref().unwrap();
        // Retained coordinates (V,u), qdot=(2,2,1)u, have M=[[14,14],[14,24]].
        // J1=-4, J2=1: S=-3, du=-S/10, dV=S/35 before finite-slip correction.
        close(
            released.state_after.base.twist_body.linear.x,
            state.base.twist_body.linear.x - 3. / 35. - defect / 7.,
        );
        close(released.state_after.joint_velocity[0], 1.);
        close(released.state_after.joint_velocity[1], 1.);
        close(released.state_after.joint_velocity[2], 0.5);
        close(event.loop_multipliers_ns[0], 104. / 35. + 2. * defect / 7.);
        let lambda = (-3. - 5. * defect) / 7.;
        close(event.joint_multipliers[0], lambda);
        close(released.supplied_work_j, 1037. / 210. - 3. * defect / 14.);
        close(released.constraint_work_j, 0.5 * lambda * defect);
        close(
            released.payloads_after[0].supplied_work_j,
            176. / 35. - 2. * defect / 7.,
        );
        close(
            released.payloads_after[1].supplied_work_j,
            -19. / 210. + defect / 14.,
        );
        assert!(released.payloads_after[1].supplied_work_j < 0.);
        assert_eq!(released.model_after, model);
        for rate in &event.velocity_residual_after_mps {
            close(*rate, 0.);
        }
        for rate in &event.joint_velocity_residuals_after {
            close(*rate, 0.);
        }
        balance(released);
        // The free first-slider candidate would be 2.4, above its 1.2 limit;
        // the simultaneous closed endpoint is 1.0 and is physically admitted.
        if defect > 0. {
            assert!(matches!(
                release_closed(
                    &loaded.model_after,
                    &state,
                    &inputs,
                    &loops,
                    &rows,
                    PayloadControls {
                        maximum_constraint_work_j: 1e-5,
                        ..controls
                    },
                    Default::default(),
                    || false,
                ),
                Err(PayloadError::ConstraintWork { .. })
            ));
        }
    }
}

#[test]
fn cancelling_tree_sources_still_eject_both_payloads_and_obey_the_total_work_cap() {
    let model = single_body(sphere(5.));
    let loaded = attach(
        &model,
        &translating(),
        &[
            (0, sphere(2.), Se3::identity()),
            (0, sphere(3.), Se3::identity()),
        ],
        &[],
        &[],
    );
    let inputs = [
        PayloadReleaseInput {
            attachment: &loaded.attachments[1],
            impulse: Some(kick(2.)),
        },
        PayloadReleaseInput {
            attachment: &loaded.attachments[0],
            impulse: Some(kick(-2.)),
        },
    ];
    let result = release_payloads(
        &loaded.model_after,
        &loaded.state_after,
        &inputs,
        &[],
        controls(&[]),
        || false,
    )
    .unwrap();
    assert_eq!(result.model_after, model);
    assert_eq!(result.state_after, loaded.state_after);
    for (actual, original) in result
        .state_after
        .base
        .twist_body
        .to_array()
        .into_iter()
        .zip(loaded.state_after.base.twist_body.to_array())
    {
        assert_eq!(actual.to_bits(), original.to_bits());
    }
    let v = loaded.state_after.base.twist_body.linear.x;
    close(
        result.payloads_after[0].payload_after.twist_body.linear.x,
        v + 2. / 3.,
    );
    close(
        result.payloads_after[1].payload_after.twist_body.linear.x,
        v - 1.,
    );
    close(result.payloads_after[0].supplied_work_j, 2. / 3.);
    close(result.payloads_after[1].supplied_work_j, 1.);
    close(result.supplied_work_j, 5. / 3.);
    balance(&result);
    assert!(matches!(
        release_payloads(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &[],
            PayloadControls {
                maximum_energy_exchange_j: 1.5,
                ..controls(&[])
            },
            || false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
}

fn solve3(a: Mat3, b: Vec3) -> Vec3 {
    // Independent cofactor inverse: no native articulated or fs-la mass solver.
    let m = a.m;
    let determinant = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    Vec3::new(
        ((m[4] * m[8] - m[5] * m[7]) * b.x
            + (m[2] * m[7] - m[1] * m[8]) * b.y
            + (m[1] * m[5] - m[2] * m[4]) * b.z)
            / determinant,
        ((m[5] * m[6] - m[3] * m[8]) * b.x
            + (m[0] * m[8] - m[2] * m[6]) * b.y
            + (m[2] * m[3] - m[0] * m[5]) * b.z)
            / determinant,
        ((m[3] * m[7] - m[4] * m[6]) * b.x
            + (m[1] * m[6] - m[0] * m[7]) * b.y
            + (m[0] * m[4] - m[1] * m[3]) * b.z)
            / determinant,
    )
}
fn rigid_delta(inertia: SpatialInertia, impulse: Wrench) -> Twist {
    let angular = solve3(
        inertia.inertia_com(),
        impulse.torque - inertia.center_of_mass().cross(impulse.force),
    );
    Twist::new(
        angular,
        impulse.force.scale(1. / inertia.mass()) - angular.cross(inertia.center_of_mass()),
    )
}

#[test]
fn offset_full_tensor_ejections_match_independent_inverses_and_frame_covariance() {
    let retained = SpatialInertia::new(
        5.,
        Vec3::new(0.1, -0.05, 0.07),
        Mat3 {
            m: [4., 0.3, -0.2, 0.3, 5., 0.15, -0.2, 0.15, 6.],
        },
    )
    .unwrap();
    let inertia = [
        SpatialInertia::new(
            2.,
            Vec3::new(-0.1, 0.07, 0.03),
            Mat3 {
                m: [0.9, 0.08, -0.04, 0.08, 1.1, 0.06, -0.04, 0.06, 1.4],
            },
        )
        .unwrap(),
        SpatialInertia::new(
            3.,
            Vec3::new(0.06, -0.02, 0.09),
            Mat3 {
                m: [1.2, -0.07, 0.03, -0.07, 1.6, -0.09, 0.03, -0.09, 1.8],
            },
        )
        .unwrap(),
    ];
    let offsets = [Vec3::new(0.3, -0.2, 0.4), Vec3::new(-0.4, 0.2, 0.1)];
    let local = [
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
            offsets[0],
        )
        .unwrap(),
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(-0.2, 0.1, -0.3))).unwrap(),
            offsets[1],
        )
        .unwrap(),
    ];
    let kicks = [
        SeparationImpulse {
            point_link_m: Vec3::new(0.2, -0.1, 0.25),
            impulse_on_payload_link_ns: Vec3::new(0.7, -0.4, 0.3),
            couple_on_payload_link_nms: Vec3::new(0.08, 0.03, -0.06),
        },
        SeparationImpulse {
            point_link_m: Vec3::new(-0.1, 0.3, -0.2),
            impulse_on_payload_link_ns: Vec3::new(-0.2, 0.5, 0.4),
            couple_on_payload_link_nms: Vec3::new(-0.02, 0.07, 0.05),
        },
    ];
    let mut retained_impulses = Vec::new();
    let mut payload_impulses = Vec::new();
    for i in 0..2 {
        let force = kicks[i].impulse_on_payload_link_ns;
        let couple = kicks[i].couple_on_payload_link_nms;
        retained_impulses.push(Wrench::new(
            (kicks[i].point_link_m.cross(force) + couple).scale(-1.),
            force.scale(-1.),
        ));
        payload_impulses.push(Wrench::new(
            local[i]
                .rotation()
                .inverse()
                .rotate((kicks[i].point_link_m - offsets[i]).cross(force) + couple)
                .unwrap(),
            local[i].rotation().inverse().rotate(force).unwrap(),
        ));
    }
    let total = Wrench::new(
        retained_impulses[0].torque + retained_impulses[1].torque,
        retained_impulses[0].force + retained_impulses[1].force,
    );
    let carrier_delta = rigid_delta(retained, total);
    let payload_delta = [
        rigid_delta(inertia[0], payload_impulses[0]),
        rigid_delta(inertia[1], payload_impulses[1]),
    ];
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(
                So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.3))).unwrap(),
                Vec3::new(1.2, -0.7, 0.4),
            )
            .unwrap(),
            Twist::new(Vec3::new(0.2, -0.1, 0.15), Vec3::new(0.3, -0.4, 0.1)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let frame = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(-0.3, 0.2, 0.1))).unwrap(),
        Vec3::new(-0.4, 0.8, 0.3),
    )
    .unwrap();
    let boost = Vec3::new(3., -1.4, 0.8);
    let mut moved = state.clone();
    moved.base.world_from_base = frame.compose(state.base.world_from_base).unwrap();
    moved.base.twist_body.linear = moved.base.twist_body.linear
        + moved
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
    let model = single_body(retained);
    let mut observations = Vec::new();
    for start in [state, moved] {
        let loaded = attach(
            &model,
            &start,
            &[(0, inertia[0], local[0]), (0, inertia[1], local[1])],
            &[],
            &[],
        );
        let inputs = [
            PayloadReleaseInput {
                attachment: &loaded.attachments[1],
                impulse: Some(kicks[1]),
            },
            PayloadReleaseInput {
                attachment: &loaded.attachments[0],
                impulse: Some(kicks[0]),
            },
        ];
        let result = release_payloads(
            &loaded.model_after,
            &loaded.state_after,
            &inputs,
            &[],
            controls(&[]),
            || false,
        )
        .unwrap();
        let expected = loaded.state_after.base.twist_body.plus(carrier_delta);
        vector(result.state_after.base.twist_body.angular, expected.angular);
        vector(result.state_after.base.twist_body.linear, expected.linear);
        let mut work = 0.;
        for (output, original) in [1usize, 0].into_iter().enumerate() {
            let expected = loaded.payloads_after[original]
                .twist_body
                .plus(payload_delta[original]);
            let actual = &result.payloads_after[output];
            vector(actual.payload_after.twist_body.angular, expected.angular);
            vector(actual.payload_after.twist_body.linear, expected.linear);
            vector(
                actual.impulse_on_payload_body.force,
                payload_impulses[original].force,
            );
            vector(
                actual.impulse_on_payload_body.torque,
                payload_impulses[original].torque,
            );
            let pair_work = 0.5
                * (retained_impulses[original].pairing(carrier_delta)
                    + payload_impulses[original].pairing(payload_delta[original]));
            close(actual.supplied_work_j, pair_work);
            work += pair_work;
        }
        close(result.supplied_work_j, work);
        assert_eq!(result.model_after, model);
        balance(&result);
        observations.push(result);
    }
    let reference = &observations[0];
    let moved = &observations[1];
    close(moved.supplied_work_j, reference.supplied_work_j);
    vector(
        moved.state_after.base.twist_body.angular,
        reference.state_after.base.twist_body.angular,
    );
    vector(
        moved.state_after.base.twist_body.linear,
        reference.state_after.base.twist_body.linear
            + moved
                .state_after
                .base
                .world_from_base
                .rotation()
                .inverse()
                .rotate(boost)
                .unwrap(),
    );
    for (actual, original) in moved.payloads_after.iter().zip(&reference.payloads_after) {
        vector(
            actual.payload_after.twist_body.angular,
            original.payload_after.twist_body.angular,
        );
        vector(
            actual.payload_after.twist_body.linear,
            original.payload_after.twist_body.linear
                + actual
                    .payload_after
                    .world_from_body
                    .rotation()
                    .inverse()
                    .rotate(boost)
                    .unwrap(),
        );
    }
}
