//! G0/G1/G3/G4: supported inventory retains material closure, recoil and source work.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::release::SeparationImpulse;
use fs_mbd::articulated::motion::inventory::supported::bilateral::{
    SupportedBilateralPayloadCapture, capture_payload, release_payload,
    release_payload_with_impulse,
};
use fs_mbd::articulated::motion::inventory::supported::constraints::{
    self, SupportedPayloadConstraintControls,
};
use fs_mbd::articulated::motion::inventory::supported::release::{
    SupportedPayloadReleaseControls, SupportedPayloadSeparationControls,
};
use fs_mbd::articulated::motion::inventory::supported::{
    self, SupportedPayloadBalance, SupportedPayloadControls,
};
use fs_mbd::articulated::motion::inventory::{PayloadError, PayloadState};
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::MaterialPoint;
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopError, DistanceLoop, PointCoincidenceLoop,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia, forward_kinematics,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tolerance {tolerance:.3e})"
    );
}
fn close(actual: f64, expected: f64) {
    near(actual, expected, 3e-9 * (1. + expected.abs()));
}
fn vector(actual: Vec3, expected: Vec3) {
    for (a, e) in [actual.x, actual.y, actual.z]
        .into_iter()
        .zip([expected.x, expected.y, expected.z])
    {
        close(a, e);
    }
}
fn wrench(actual: Wrench, expected: Wrench) {
    vector(actual.torque, expected.torque);
    vector(actual.force, expected.force);
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
fn sliders(maximum_speed: Option<f64>) -> ArticulatedModel {
    let mut links = vec![Link::new(
        "support carrier",
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
                maximum_speed.map(|v| JointLimits::new(-10., 10., v, 0.1).unwrap()),
            )
            .unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 1.25,
        base: BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
            Twist::new(Vec3::default(), Vec3::new(-0.2, 0.1, 0.3)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
    }
}
fn incoming() -> PayloadState {
    PayloadState {
        inertia: sphere(2.),
        world_from_body: Se3::from_parts(So3::identity(), Vec3::new(1., 0., 0.)).unwrap(),
        twist_body: Twist::new(Vec3::default(), Vec3::new(2.3, 0., 0.)),
    }
}
fn loops() -> [BilateralLoop; 1] {
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
fn rows() -> [JointConstraint; 1] {
    [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![0., 1., -2.],
        position: 0.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }]
}
fn controls(tolerances: &[f64]) -> SupportedPayloadControls<'_> {
    SupportedPayloadControls {
        angular_momentum_tolerance_nms: 2e-8,
        linear_momentum_tolerance_ns: 2e-8,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_removed_relative_energy_j: 20.,
    }
}
fn release_controls(tolerances: &[f64]) -> SupportedPayloadReleaseControls<'_> {
    SupportedPayloadReleaseControls {
        angular_momentum_tolerance_nms: 2e-8,
        linear_momentum_tolerance_ns: 2e-8,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
    }
}
fn reaction_controls(maximum_work: f64) -> SupportedPayloadConstraintControls {
    SupportedPayloadConstraintControls {
        minimum_reciprocal_condition: 1e-10,
        maximum_reaction_work_j: maximum_work,
    }
}
fn kick(impulse: f64) -> SeparationImpulse {
    SeparationImpulse {
        impulse_on_payload_link_ns: Vec3::new(impulse, 0., 0.),
        ..Default::default()
    }
}
fn attach(model: &ArticulatedModel) -> SupportedBilateralPayloadCapture {
    capture_payload(
        model,
        &initial(),
        1,
        incoming(),
        &loops(),
        &rows(),
        controls(&[2e-8; 3]),
        reaction_controls(0.1),
        Default::default(),
        || false,
    )
    .unwrap()
}
fn audit(balance: &SupportedPayloadBalance) {
    vector(balance.linear_momentum_residual_world, Vec3::default());
    vector(balance.angular_momentum_residual_world, Vec3::default());
    for &residual in &balance.joint_momentum_residual {
        close(residual, 0.);
    }
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    close(balance.energy_residual_j, 0.);
}

#[test]
fn capture_and_ejection_use_the_combined_endpoint_and_actual_support_work() {
    let model = sliders(Some(0.8));
    let start = initial();
    // The gear-only candidate reaches qdot_1=1, above the native limit.
    assert!(
        constraints::capture_payload(
            &model,
            &start,
            1,
            incoming(),
            &rows(),
            controls(&[2e-8; 3]),
            reaction_controls(0.1),
            || false,
        )
        .is_err()
    );
    let result = attach(&model);
    let capture = &result.capture;
    // With the support fixed in velocity, qdot=(2,2,1)u. The loaded
    // reduced inertia is 32 and the mismatch contributes 4.8 to that mode.
    for (&v, factor) in capture.state_after.joint_velocity.iter().zip([2., 2., 1.]) {
        close(v, factor * 0.35);
    }
    assert_eq!(capture.state_after.base, start.base);
    assert_eq!(capture.state_after.time_s, start.time_s);
    assert_eq!(capture.state_after.joint_position, start.joint_position);
    close(result.reaction.impulse.loop_multipliers_ns[0], 1.2);
    close(result.reaction.impulse.joint_multipliers[0], -0.3);
    assert_eq!(result.reaction.generalized_constraint_impulse.len(), 9);
    for (&actual, expected) in result
        .reaction
        .generalized_constraint_impulse
        .iter()
        .zip([0., 0., 0., 0., 0., 0., -1.2, 0.9, 0.6])
    {
        close(actual, expected);
    }
    close(capture.balance.support_impulse_body.force.x, 0.3);
    close(capture.balance.support_work_j, 0.21);
    close(
        capture.balance.momentum_before.linear_momentum_world.x,
        17.2,
    );
    close(capture.balance.momentum_after.linear_momentum_world.x, 17.5);
    close(capture.balance.momentum_before.kinetic_energy, 11.16);
    close(capture.balance.momentum_after.kinetic_energy, 10.29);
    close(capture.balance.kinetic_energy_change_j, -0.87);
    close(capture.removed_relative_energy_j, 1.08);
    close(result.reaction.midpoint_work_j, 0.);
    audit(&capture.balance);

    let passive = release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        &loops(),
        &rows(),
        release_controls(&[2e-8; 3]),
        reaction_controls(0.1),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(passive.release.model_after, model);
    assert_eq!(passive.release.state_after, capture.state_after);
    close(passive.release.payload_after.twist_body.linear.x, 1.4);
    wrench(
        passive.release.balance.support_impulse_body,
        Wrench::default(),
    );
    close(passive.release.supplied_work_j, 0.);
    audit(&passive.release.balance);

    // Restored inertia is 24. The unreacted first-slider candidate is -1.3;
    // the physical simultaneous endpoint is 1/30 and passes the 0.8 limit.
    let separated = release_payload_with_impulse(
        &capture.model_after,
        &capture.state_after,
        &capture.attachment,
        kick(4.),
        &loops(),
        &rows(),
        SupportedPayloadSeparationControls {
            balance: release_controls(&[2e-8; 3]),
            maximum_separation_work_j: 6.,
        },
        reaction_controls(0.1),
        Default::default(),
        || false,
    )
    .unwrap();
    let release = &separated.release;
    assert_eq!(release.model_after, model);
    assert_eq!(release.state_after.base, start.base);
    for (&v, factor) in release.state_after.joint_velocity.iter().zip([2., 2., 1.]) {
        close(v, factor * (0.35 - 4. / 12.));
    }
    close(separated.reaction.impulse.loop_multipliers_ns[0], -8. / 3.);
    close(separated.reaction.impulse.joint_multipliers[0], 2. / 3.);
    close(release.payload_after.twist_body.linear.x, 3.4);
    close(release.balance.support_impulse_body.force.x, -2. / 3.);
    close(release.balance.support_work_j, -7. / 15.);
    close(release.supplied_work_j, 16. / 3.);
    close(release.balance.kinetic_energy_change_j, 73. / 15.);
    close(separated.reaction.midpoint_work_j, 0.);
    close(separated.reaction.impulse.work_partition_residual_j, 0.);
    audit(&release.balance);
    assert_eq!(start, initial());
}

#[test]
fn finite_tangency_keeps_constraint_work_separate_and_zero_sources_passive() {
    let model = sliders(None);
    let defect = 1. / 1024.;
    let mut state = initial();
    state.joint_velocity[2] -= 0.5 * defect;
    let mut gear = rows();
    gear[0].velocity_tolerance = 2. * defect;
    let capture = capture_payload(
        &model,
        &state,
        1,
        incoming(),
        &loops(),
        &gear,
        controls(&[2e-8; 3]),
        reaction_controls(0.1),
        Default::default(),
        || false,
    )
    .unwrap();
    let lambda = -0.3 - 7. * defect / 8.;
    close(
        capture.reaction.impulse.loop_multipliers_ns[0],
        1.2 + defect / 2.,
    );
    close(capture.reaction.impulse.joint_multipliers[0], lambda);
    close(capture.reaction.midpoint_work_j, 0.5 * lambda * defect);
    close(
        capture.capture.state_after.joint_velocity[2],
        0.35 - defect / 16.,
    );
    let after = &capture.capture.state_after.joint_velocity;
    // Direct physical relative-velocity quadratic, including every retained
    // slider; midpoint constraint work must not be subtracted from this loss.
    let loss = [2., 3., 4.]
        .into_iter()
        .zip(&state.joint_velocity)
        .zip(after)
        .map(|((m, before), after)| 0.5 * m * (after - before).powi(2))
        .sum::<f64>()
        + (0.7 + after[0] - 2.3).powi(2);
    close(capture.capture.removed_relative_energy_j, loss);
    close(
        capture.capture.balance.energy_residual_j,
        capture.reaction.final_rate_pairing_j,
    );
    assert!(matches!(
        capture_payload(
            &model,
            &state,
            1,
            incoming(),
            &loops(),
            &gear,
            controls(&[2e-8; 3]),
            reaction_controls(1e-5),
            Default::default(),
            || false,
        ),
        Err(PayloadError::ConstraintWork { .. })
    ));

    let mut comoving = incoming();
    comoving.twist_body.linear.x = state.base.twist_body.linear.x + state.joint_velocity[0];
    let quiet = capture_payload(
        &model,
        &state,
        1,
        comoving,
        &loops(),
        &gear,
        controls(&[2e-8; 3]),
        reaction_controls(0.),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(quiet.capture.state_after, state);
    assert_eq!(quiet.capture.removed_relative_energy_j, 0.);
    assert_eq!(quiet.reaction.impulse.loop_multipliers_ns, vec![0.]);
    assert_eq!(quiet.reaction.impulse.joint_multipliers, vec![0.]);
    close(
        quiet.reaction.impulse.joint_velocity_residuals_after[0],
        defect,
    );

    let loaded = attach(&model).capture;
    let mut slipping = loaded.state_after.clone();
    slipping.joint_velocity[2] -= 0.5 * defect;
    for source in [None, Some(SeparationImpulse::default())] {
        let result = if let Some(source) = source {
            release_payload_with_impulse(
                &loaded.model_after,
                &slipping,
                &loaded.attachment,
                source,
                &loops(),
                &gear,
                SupportedPayloadSeparationControls {
                    balance: release_controls(&[2e-8; 3]),
                    maximum_separation_work_j: 0.,
                },
                reaction_controls(0.),
                Default::default(),
                || false,
            )
        } else {
            release_payload(
                &loaded.model_after,
                &slipping,
                &loaded.attachment,
                &loops(),
                &gear,
                release_controls(&[2e-8; 3]),
                reaction_controls(0.),
                Default::default(),
                || false,
            )
        }
        .unwrap();
        assert_eq!(result.release.state_after, slipping);
        assert_eq!(result.release.model_after, model);
        assert_eq!(result.release.applied_impulse, source);
        close(
            result.reaction.impulse.joint_velocity_residuals_after[0],
            defect,
        );
        assert_eq!(result.reaction.impulse.loop_multipliers_ns, vec![0.]);
        assert_eq!(result.reaction.impulse.joint_multipliers, vec![0.]);
        assert_eq!(result.release.supplied_work_j, 0.);
        assert_eq!(result.reaction.midpoint_work_j, 0.);
        audit(&result.release.balance);
    }
    let separated = release_payload_with_impulse(
        &loaded.model_after,
        &slipping,
        &loaded.attachment,
        kick(4.),
        &loops(),
        &gear,
        SupportedPayloadSeparationControls {
            balance: release_controls(&[2e-8; 3]),
            maximum_separation_work_j: 6.,
        },
        reaction_controls(0.1),
        Default::default(),
        || false,
    )
    .unwrap();
    let lambda = 4. / 6. - 5. * defect / 6.;
    close(
        separated.reaction.impulse.loop_multipliers_ns[0],
        -8. / 3. + defect / 3.,
    );
    close(separated.reaction.impulse.joint_multipliers[0], lambda);
    close(separated.reaction.midpoint_work_j, 0.5 * lambda * defect);
    close(separated.release.supplied_work_j, 16. / 3. + defect / 3.);
    close(
        separated.release.balance.support_impulse_body.force.x,
        -lambda,
    );
    close(
        separated.release.balance.kinetic_energy_change_j,
        separated.release.balance.support_work_j
            + separated.release.supplied_work_j
            + separated.reaction.midpoint_work_j,
    );
    audit(&separated.release.balance);
    // This source extracts delta^2/192 J from the admitted finite slip.
    // Supported inventory intentionally has no regenerative-source policy.
    assert!(matches!(
        release_payload_with_impulse(
            &loaded.model_after,
            &slipping,
            &loaded.attachment,
            kick(-defect / 8.),
            &loops(),
            &gear,
            SupportedPayloadSeparationControls {
                balance: release_controls(&[2e-8; 3]),
                maximum_separation_work_j: 1.,
            },
            reaction_controls(0.1),
            Default::default(),
            || false,
        ),
        Err(PayloadError::Invalid("negative supported separation work"))
    ));
}

fn solve3(a: Mat3, b: Vec3) -> Vec3 {
    // Independent centroidal cofactor inverse; no articulated mass solver.
    let m = a.m;
    let det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    Vec3::new(
        ((m[4] * m[8] - m[5] * m[7]) * b.x
            + (m[2] * m[7] - m[1] * m[8]) * b.y
            + (m[1] * m[5] - m[2] * m[4]) * b.z)
            / det,
        ((m[5] * m[6] - m[3] * m[8]) * b.x
            + (m[0] * m[8] - m[2] * m[6]) * b.y
            + (m[2] * m[3] - m[0] * m[5]) * b.z)
            / det,
        ((m[3] * m[7] - m[4] * m[6]) * b.x
            + (m[1] * m[6] - m[0] * m[7]) * b.y
            + (m[0] * m[4] - m[1] * m[3]) * b.z)
            / det,
    )
}
fn virtual_work_projection(
    model: &ArticulatedModel,
    state: &PrescribedBaseMotionState,
    body: &[Wrench],
    joint: &[f64],
) -> Vec<f64> {
    // Kinematic virtual work checks every base/joint covector independently
    // of the production zero-rate inverse-dynamics projection.
    let axes = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ];
    (0..6 + model.dof_count())
        .map(|column| {
            let base_twist = if column < 3 {
                Twist::new(axes[column], Vec3::default())
            } else if column < 6 {
                Twist::new(Vec3::default(), axes[column - 3])
            } else {
                Twist::zero()
            };
            let mut velocity = vec![0.; model.dof_count()];
            if column >= 6 {
                velocity[column - 6] = 1.;
            }
            let kin = forward_kinematics(
                model,
                BaseState::prescribed(state.base.world_from_base, base_twist, Twist::zero()),
                &state.joint_position,
                &velocity,
            )
            .unwrap();
            body.iter()
                .zip(kin.body_twist)
                .map(|(j, v)| j.pairing(v))
                .sum::<f64>()
                + joint.iter().zip(velocity).map(|(j, v)| j * v).sum::<f64>()
        })
        .collect()
}

#[test]
fn mobile_root_spherical_closure_and_locks_match_full_tensor_rigid_references() {
    let zero = Vec3::default();
    let base_from_root_zero = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.3))).unwrap(),
        Vec3::new(0.6, -0.3, 0.2),
    )
    .unwrap();
    let mut links = vec![Link::new(
        "mobile root",
        None,
        base_from_root_zero,
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        sphere(5.),
    )];
    for (index, axis) in [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ]
    .into_iter()
    .enumerate()
    {
        links.push(Link::new(
            format!("translation {index}"),
            Some(index),
            Se3::identity(),
            JointModel::prismatic(axis, None).unwrap(),
            sphere(2. + index as f64),
        ));
    }
    links.push(Link::new(
        "payload mount",
        Some(3),
        Se3::identity(),
        JointModel::revolute(Vec3::new(0., 1., 0.), None).unwrap(),
        sphere(1.5),
    ));
    let model = ArticulatedModel::new(links).unwrap();
    let material = [BilateralLoop::PointCoincidence(PointCoincidenceLoop {
        a: MaterialPoint {
            link: 0,
            point_body: zero,
        },
        b: MaterialPoint {
            link: 3,
            point_body: zero,
        },
    })];
    let locks: Vec<_> = [(0, 0.5), (4, 0.)]
        .into_iter()
        .map(|(index, position)| {
            let mut coefficients = vec![0.; 5];
            coefficients[index] = 1.;
            JointConstraint {
                coordinate: ConstraintCoordinate::Angle,
                coefficients,
                position,
                position_tolerance: 1e-9,
                velocity_tolerance: 1e-9,
            }
        })
        .collect();
    let inertia = SpatialInertia::new(
        2.,
        Vec3::new(-0.1, 0.07, 0.03),
        Mat3 {
            m: [0.9, 0.08, -0.04, 0.08, 1.1, 0.06, -0.04, 0.06, 1.4],
        },
    )
    .unwrap();
    let local = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
        Vec3::new(0.3, -0.2, 0.4),
    )
    .unwrap();
    let slip = Twist::new(Vec3::new(0.3, -0.1, 0.2), Vec3::new(-0.4, 0.6, 0.2));
    let source = SeparationImpulse {
        point_link_m: Vec3::new(0.2, -0.1, 0.25),
        impulse_on_payload_link_ns: Vec3::new(0.7, -0.4, 0.3),
        couple_on_payload_link_nms: Vec3::new(0.08, 0.03, -0.06),
    };
    let source_link = Wrench::new(
        source.point_link_m.cross(source.impulse_on_payload_link_ns)
            + source.couple_on_payload_link_nms,
        source.impulse_on_payload_link_ns,
    );
    let source_body = source_link.transform_by(&local.inverse().unwrap());
    let angular = solve3(
        inertia.inertia_com(),
        source_body.torque - inertia.center_of_mass().cross(source_body.force),
    );
    let payload_delta = Twist::new(
        angular,
        source_body.force.scale(1. / inertia.mass()) - angular.cross(inertia.center_of_mass()),
    );
    let reference_source_work = 0.5 * source_body.pairing(payload_delta);
    let frame = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(-0.3, 0.2, 0.1))).unwrap(),
        Vec3::new(-0.4, 0.8, 0.3),
    )
    .unwrap();
    for pose in [Se3::identity(), frame] {
        // Exact immobilization on a stationary support lets every residual
        // joint response remain representable. Moving supports and boosts
        // are covered by the one-mode fixture above and below.
        let base_twist = Twist::zero();
        let state = PrescribedBaseMotionState {
            time_s: 3.,
            base: BaseState::prescribed(pose, base_twist, Twist::zero()),
            joint_position: vec![0.5, 0., 0., 0., 0.],
            joint_velocity: vec![0.; 5],
        };
        let kin = forward_kinematics(
            &model,
            state.base,
            &state.joint_position,
            &state.joint_velocity,
        )
        .unwrap();
        let carried = kin.body_twist[4].transform_by(&local.inverse().unwrap());
        let payload = PayloadState {
            inertia,
            world_from_body: kin.world_from_link[4].compose(local).unwrap(),
            twist_body: carried.plus(slip),
        };
        let base_from_target = pose
            .inverse()
            .unwrap()
            .compose(kin.world_from_link[4])
            .unwrap();
        let mismatch = inertia
            .momentum(slip)
            .unwrap()
            .transform_by(&local)
            .transform_by(&base_from_target);
        let expected_support = Wrench::new(mismatch.torque.scale(-1.), mismatch.force.scale(-1.));
        let captured = capture_payload(
            &model,
            &state,
            4,
            payload,
            &material,
            &locks,
            controls(&[2e-8; 5]),
            reaction_controls(1e-6),
            Default::default(),
            || false,
        )
        .unwrap();
        for &v in &captured.capture.state_after.joint_velocity {
            close(v, 0.);
        }
        assert_eq!(captured.capture.state_after.base, state.base);
        vector(
            captured.capture.payload_after.twist_body.angular,
            carried.angular,
        );
        vector(
            captured.capture.payload_after.twist_body.linear,
            carried.linear,
        );
        wrench(
            captured.capture.balance.support_impulse_body,
            expected_support,
        );
        wrench(
            captured.capture.balance.support_impulse_world,
            expected_support.transform_by(&pose),
        );
        close(
            captured.capture.balance.support_work_j,
            expected_support.pairing(base_twist),
        );
        close(
            captured.capture.removed_relative_energy_j,
            inertia.kinetic_energy(slip).unwrap(),
        );
        let projected = virtual_work_projection(
            &captured.capture.model_after,
            &state,
            &captured.reaction.impulse.constraint_body_impulses,
            &captured.reaction.impulse.constraint_joint_impulses,
        );
        assert_eq!(projected.len(), 11);
        assert_eq!(captured.reaction.generalized_constraint_impulse.len(), 11);
        for (&actual, expected) in captured
            .reaction
            .generalized_constraint_impulse
            .iter()
            .zip(projected)
        {
            close(actual, expected);
        }
        close(captured.reaction.midpoint_work_j, 0.);
        close(captured.reaction.final_rate_pairing_j, 0.);
        audit(&captured.capture.balance);
        let released = release_payload_with_impulse(
            &captured.capture.model_after,
            &captured.capture.state_after,
            &captured.capture.attachment,
            source,
            &material,
            &locks,
            SupportedPayloadSeparationControls {
                balance: release_controls(&[2e-8; 5]),
                maximum_separation_work_j: 2.,
            },
            reaction_controls(1e-6),
            Default::default(),
            || false,
        )
        .unwrap();
        for &v in &released.release.state_after.joint_velocity {
            close(v, 0.);
        }
        assert_eq!(released.release.model_after, model);
        wrench(released.release.impulse_on_payload_body, source_body);
        vector(
            released.release.payload_after.twist_body.angular,
            carried.angular + payload_delta.angular,
        );
        vector(
            released.release.payload_after.twist_body.linear,
            carried.linear + payload_delta.linear,
        );
        let expected_support = source_link.transform_by(&base_from_target);
        wrench(
            released.release.balance.support_impulse_body,
            expected_support,
        );
        close(
            released.release.balance.support_work_j,
            expected_support.pairing(base_twist),
        );
        close(released.release.supplied_work_j, reference_source_work);
        close(
            released.release.balance.kinetic_energy_change_j,
            expected_support.pairing(base_twist) + reference_source_work,
        );
        audit(&released.release.balance);
    }
}

#[test]
fn large_common_boost_preserves_conditioned_capture_loss_and_separation_work() {
    let model = sliders(None);
    for boost in [0., 1_073_741_824.] {
        let mut state = initial();
        state.base.twist_body.linear.x += boost;
        let mut payload = incoming();
        payload.twist_body.linear.x += boost;
        let tolerance = [2e-6; 3];
        let captured = capture_payload(
            &model,
            &state,
            1,
            payload,
            &loops(),
            &rows(),
            SupportedPayloadControls {
                angular_momentum_tolerance_nms: 2e-6,
                linear_momentum_tolerance_ns: 2e-6,
                energy_tolerance_j: 4e-6,
                ..controls(&tolerance)
            },
            reaction_controls(4e-6),
            Default::default(),
            || false,
        )
        .unwrap();
        // Construct the independent mode from the actually representable
        // incoming speed in the support frame, rather than subtracting O(V^2).
        let incoming_relative = payload.twist_body.linear.x - state.base.twist_body.linear.x;
        let mode = (4. * 0.4 + 6. * 0.4 + 4. * 0.2 + 4. * incoming_relative) / 32.;
        let loss = 2. * (2. * mode - 0.4).powi(2) / 2.
            + 3. * (2. * mode - 0.4).powi(2) / 2.
            + 4. * (mode - 0.2).powi(2) / 2.
            + (2. * mode - incoming_relative).powi(2);
        near(captured.capture.state_after.joint_velocity[2], mode, 3e-7);
        near(captured.capture.removed_relative_energy_j, loss, 4e-6);
        assert!(captured.capture.balance.energy_residual_j.abs() <= 4e-6);
        assert!(
            captured
                .capture
                .balance
                .maximum_normalized_momentum_residual
                <= 1.
        );
        near(
            captured.capture.balance.support_work_j,
            captured.capture.balance.support_impulse_body.force.x * state.base.twist_body.linear.x,
            2e-6,
        );
        let released = release_payload_with_impulse(
            &captured.capture.model_after,
            &captured.capture.state_after,
            &captured.capture.attachment,
            kick(4.),
            &loops(),
            &rows(),
            SupportedPayloadSeparationControls {
                balance: SupportedPayloadReleaseControls {
                    angular_momentum_tolerance_nms: 2e-6,
                    linear_momentum_tolerance_ns: 2e-6,
                    energy_tolerance_j: 4e-6,
                    ..release_controls(&tolerance)
                },
                maximum_separation_work_j: 6.,
            },
            reaction_controls(4e-6),
            Default::default(),
            || false,
        )
        .unwrap();
        near(released.release.supplied_work_j, 16. / 3., 4e-6);
        near(
            released.release.balance.support_impulse_body.force.x,
            -2. / 3.,
            2e-6,
        );
        assert!(released.release.balance.energy_residual_j.abs() <= 4e-6);
        assert!(
            released
                .release
                .balance
                .maximum_normalized_momentum_residual
                <= 1.
        );
    }
}

#[test]
fn empty_rows_and_atomic_refusals_retain_native_admission_and_inventory_identity() {
    let model = sliders(None);
    let state = initial();
    let input = incoming();
    let allowances = [2e-8; 3];
    let empty = capture_payload(
        &model,
        &state,
        1,
        input,
        &[],
        &[],
        controls(&allowances),
        reaction_controls(0.),
        Default::default(),
        || false,
    )
    .unwrap();
    let original =
        supported::capture_payload(&model, &state, 1, input, controls(&allowances), || false)
            .unwrap();
    assert_eq!(empty.capture.model_after, original.model_after);
    assert_eq!(empty.capture.state_after, original.state_after);
    close(
        empty.capture.removed_relative_energy_j,
        original.removed_relative_energy_j,
    );
    assert!(empty.reaction.impulse.loop_multipliers_ns.is_empty());
    assert!(empty.reaction.impulse.joint_multipliers.is_empty());
    assert_eq!(empty.reaction.generalized_constraint_impulse, vec![0.; 9]);
    let duplicate = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![-1., 1., 0.],
        position: 2.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }];
    assert!(matches!(
        capture_payload(
            &model,
            &state,
            1,
            input,
            &loops(),
            &duplicate,
            controls(&allowances),
            reaction_controls(0.1),
            Default::default(),
            || false,
        ),
        Err(PayloadError::Bilateral(BilateralLoopError::Mobility { .. }))
    ));
    assert!(matches!(
        capture_payload(
            &model,
            &state,
            1,
            input,
            &loops(),
            &rows(),
            controls(&allowances),
            SupportedPayloadConstraintControls {
                minimum_reciprocal_condition: 0.99,
                ..reaction_controls(0.1)
            },
            Default::default(),
            || false,
        ),
        Err(PayloadError::Constraint(
            ConstraintError::IllConditioned { .. }
        ))
    ));
    assert!(matches!(
        capture_payload(
            &model,
            &state,
            1,
            input,
            &loops(),
            &rows(),
            SupportedPayloadControls {
                maximum_removed_relative_energy_j: 1.,
                ..controls(&allowances)
            },
            reaction_controls(0.1),
            Default::default(),
            || false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    let captured = attach(&model);
    assert!(matches!(
        release_payload(
            &model,
            &captured.capture.state_after,
            &captured.capture.attachment,
            &loops(),
            &rows(),
            release_controls(&allowances),
            reaction_controls(0.1),
            Default::default(),
            || false,
        ),
        Err(PayloadError::Invalid(_))
    ));
    assert!(matches!(
        release_payload_with_impulse(
            &captured.capture.model_after,
            &captured.capture.state_after,
            &captured.capture.attachment,
            kick(4.),
            &loops(),
            &rows(),
            SupportedPayloadSeparationControls {
                balance: release_controls(&allowances),
                maximum_separation_work_j: 5.,
            },
            reaction_controls(0.1),
            Default::default(),
            || false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    let polls = Cell::new(0);
    let expected = capture_payload(
        &model,
        &state,
        1,
        input,
        &loops(),
        &rows(),
        controls(&allowances),
        reaction_controls(0.1),
        Default::default(),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    for limit in [1, polls.get()] {
        let calls = Cell::new(0);
        assert_eq!(
            capture_payload(
                &model,
                &state,
                1,
                input,
                &loops(),
                &rows(),
                controls(&allowances),
                reaction_controls(0.1),
                Default::default(),
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= limit
                },
            ),
            Err(PayloadError::Cancelled)
        );
    }
    assert_eq!(expected, captured);
    assert_eq!(state, initial());
    assert_eq!(input, incoming());
}
