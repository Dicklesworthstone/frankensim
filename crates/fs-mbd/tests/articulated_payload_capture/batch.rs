//! G1/G3/G4 simultaneous capture checked by independent COM/screw sums and
//! scalar mechanics. Native execution is required; parsing is not test success.

use super::*;
use fs_mbd::articulated::motion::inventory::release::release_payload;
use fs_mbd::articulated::motion::inventory::{
    PayloadBatchCapture, PayloadCaptureInput, capture_payloads,
};

fn requests() -> [PayloadCaptureInput; 2] {
    [
        PayloadCaptureInput {
            target_link: 0,
            payload: PayloadState {
                inertia: sphere(1., 1.),
                world_from_body: Se3::identity(),
                twist_body: Twist::new(Vec3::new(0., 0., -0.8), Vec3::new(-1., 0., 0.)),
            },
        },
        PayloadCaptureInput {
            target_link: 1,
            payload: PayloadState {
                inertia: sphere(4., 4.),
                world_from_body: Se3::identity(),
                twist_body: Twist::new(Vec3::new(0., 0., 0.2), Vec3::new(1., 0., 0.)),
            },
        },
    ]
}

fn check_batch(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadCaptureInput],
    rows: &[JointConstraint],
    result: &PayloadBatchCapture,
) {
    let before = physical_many(
        model,
        state,
        inputs
            .iter()
            .map(|input| (input.target_link, &input.payload)),
        Vec3::default(),
    );
    let after = physical(
        &result.model_after,
        &result.state_after,
        None,
        Vec3::default(),
    );
    vector(after.linear, before.linear);
    vector(after.angular, before.angular);
    let balance = &result.balance;
    vector(balance.momentum_before.linear_momentum_world, before.linear);
    vector(
        balance.momentum_before.angular_momentum_world,
        before.angular,
    );
    vector(balance.momentum_after.linear_momentum_world, after.linear);
    vector(balance.momentum_after.angular_momentum_world, after.angular);
    assert_eq!(
        result.state_after.base.world_from_base,
        state.base.world_from_base
    );
    assert_eq!(result.state_after.joint_position, state.joint_position);
    assert_eq!(balance.multipliers.len(), rows.len());
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
    let old = physical_many(
        model,
        state,
        inputs
            .iter()
            .map(|input| (input.target_link, &input.payload)),
        balance.accounting_boost_world,
    );
    let new = physical(
        &result.model_after,
        &result.state_after,
        None,
        balance.accounting_boost_world,
    );
    close(balance.kinetic_before_j, old.kinetic);
    close(balance.kinetic_after_j, new.kinetic);
    close(
        balance.energy_residual_j,
        new.kinetic - old.kinetic + result.removed_relative_energy_j,
    );
    close(balance.energy_residual_j, 0.);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    let kin = kinematics(&result.model_after, &result.state_after);
    assert_eq!(result.attachments.len(), inputs.len());
    assert_eq!(result.payloads_after.len(), inputs.len());
    for ((input, attached), handle) in inputs
        .iter()
        .zip(&result.payloads_after)
        .zip(&result.attachments)
    {
        let target = input.target_link;
        let frame = kin.world_from_link[target];
        let twist = kin.body_twist[target];
        let omega = frame.rotation().rotate(twist.angular).unwrap();
        let origin_velocity = frame.rotation().rotate(twist.linear).unwrap();
        assert_eq!(handle.target_link(), target);
        assert_eq!(handle.payload_inertia(), input.payload.inertia);
        assert_eq!(attached.inertia, input.payload.inertia);
        vector(
            attached.world_from_body.translation(),
            input.payload.world_from_body.translation(),
        );
        for (actual, expected) in attached
            .world_from_body
            .rotation()
            .matrix()
            .m
            .into_iter()
            .zip(input.payload.world_from_body.rotation().matrix().m)
        {
            close(actual, expected);
        }
        vector(
            attached
                .world_from_body
                .rotation()
                .rotate(attached.twist_body.angular)
                .unwrap(),
            omega,
        );
        vector(
            attached
                .world_from_body
                .rotation()
                .rotate(attached.twist_body.linear)
                .unwrap(),
            origin_velocity
                + omega.cross(attached.world_from_body.translation() - frame.translation()),
        );
    }
}

// Direct sum of every central tensor shifted to the independently reconstructed
// common COM. This does not use the production reduced-mass inventory fold.
fn check_loaded_inertias(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    inputs: &[PayloadCaptureInput],
    result: &PayloadBatchCapture,
) {
    let kin = kinematics(model, state);
    for (target, link) in model.links().iter().enumerate() {
        let old = link.inertia();
        let mut pieces = vec![(old.mass(), old.center_of_mass(), old.inertia_com())];
        for input in inputs.iter().filter(|input| input.target_link == target) {
            let p = input.payload;
            let relative = kin.world_from_link[target]
                .inverse()
                .unwrap()
                .compose(p.world_from_body)
                .unwrap();
            let axes = relative.rotation().matrix();
            pieces.push((
                p.inertia.mass(),
                relative
                    .transform_point(p.inertia.center_of_mass())
                    .unwrap(),
                axes.compose(p.inertia.inertia_com())
                    .compose(axes.transpose()),
            ));
        }
        let mass: f64 = pieces.iter().map(|p| p.0).sum();
        let center = pieces
            .iter()
            .fold(Vec3::default(), |sum, p| sum + p.1.scale(p.0))
            .scale(1. / mass);
        let mut tensor = [0.; 9];
        for (mass, x, central) in pieces {
            let r = x - center;
            let a = [r.x, r.y, r.z];
            for i in 0..3 {
                for j in 0..3 {
                    tensor[3 * i + j] += central.m[3 * i + j]
                        + mass * (if i == j { r.dot(r) } else { 0. } - a[i] * a[j]);
                }
            }
        }
        let actual = result.model_after.links()[target].inertia();
        close(actual.mass(), mass);
        vector(actual.center_of_mass(), center);
        for (actual, expected) in actual.inertia_com().m.into_iter().zip(tensor) {
            close(actual, expected);
        }
    }
}

#[test]
fn g1_batch_different_links_share_recoil_and_resolve_loss_under_a_large_boost() {
    let (model, state, _) = coax(None, 1.2);
    let inputs = requests();
    let settings = controls(&[2e-10], 3.);
    let result = capture_payloads(&model, &state, &inputs, &[], settings, || false).unwrap();
    check_batch(&model, &state, &inputs, &[], &result);
    check_loaded_inertias(&model, &state, &inputs, &result);
    // Carrier spin: (2*.4 + 1*(-.8))/3 = 0. Rotor spin:
    // (1*1.6 + 4*.2)/5 = .48. Translation: (3*.5-1+4)/8.
    vector(result.state_after.base.twist_body.angular, Vec3::default());
    vector(
        result.state_after.base.twist_body.linear,
        Vec3::new(9. / 16., 0., 0.),
    );
    close(result.state_after.joint_velocity[0], 12. / 25.);
    close(result.removed_relative_energy_j, 158. / 125. + 103. / 64.);
    let boost = Vec3::new(2_f64.powi(30), -2_f64.powi(29), 2_f64.powi(28));
    let mut moving_state = state.clone();
    moving_state.base.twist_body.linear = moving_state.base.twist_body.linear + boost;
    let moving_inputs = inputs.map(|mut input| {
        input.payload.twist_body.linear = input.payload.twist_body.linear + boost;
        input
    });
    let moving = capture_payloads(&model, &moving_state, &moving_inputs, &[], settings, || {
        false
    })
    .unwrap();
    vector(
        moving.state_after.base.twist_body.linear - boost,
        result.state_after.base.twist_body.linear,
    );
    vector(
        moving.state_after.base.twist_body.angular,
        result.state_after.base.twist_body.angular,
    );
    close(
        moving.state_after.joint_velocity[0],
        result.state_after.joint_velocity[0],
    );
    close(
        moving.removed_relative_energy_j,
        result.removed_relative_energy_j,
    );
    close(
        moving.balance.kinetic_before_j,
        result.balance.kinetic_before_j,
    );
    close(
        moving.balance.kinetic_after_j,
        result.balance.kinetic_after_j,
    );
    assert!(moving.balance.energy_residual_j.abs() <= settings.energy_tolerance_j);
    for (a, b) in moving.payloads_after.iter().zip(&result.payloads_after) {
        vector(a.twist_body.linear - boost, b.twist_body.linear);
        vector(a.twist_body.angular, b.twist_body.angular);
    }
}

#[test]
fn g1_g4_repeated_target_capture_admits_only_the_final_rates_and_releases_exact_inventory() {
    let limits = Some(JointLimits::new(-2., 2., 0.5, 100.).unwrap());
    let (model, mut state, _) = coax(limits, 0.);
    state.base.twist_body = Twist::zero();
    let inputs = [4., -4.].map(|omega| PayloadCaptureInput {
        target_link: 1,
        payload: PayloadState {
            inertia: sphere(1., 1.),
            world_from_body: Se3::identity(),
            twist_body: Twist::new(Vec3::new(0., 0., omega), Vec3::default()),
        },
    });
    let settings = controls(&[2e-10], 17.);
    // Either capture in isolation asks for |qdot|=2, outside the .5 limit.
    for input in inputs {
        assert!(
            capture_payload(
                &model,
                &state,
                input.target_link,
                input.payload,
                &[],
                settings,
                || false
            )
            .is_err()
        );
    }
    let result = capture_payloads(&model, &state, &inputs, &[], settings, || false).unwrap();
    check_batch(&model, &state, &inputs, &[], &result);
    close(result.state_after.joint_velocity[0], 0.);
    vector(result.state_after.base.twist_body.angular, Vec3::default());
    close(result.removed_relative_energy_j, 16.);
    assert_eq!(
        result.attachments[0].retained_inertia(),
        model.links()[1].inertia()
    );
    assert_eq!(
        result.attachments[1].retained_inertia(),
        result.attachments[0].loaded_inertia()
    );
    assert_eq!(
        result.attachments[1].loaded_inertia(),
        result.model_after.links()[1].inertia()
    );
    assert!(
        release_payload(
            &result.model_after,
            &result.state_after,
            &result.attachments[0],
            None,
            &[],
            settings,
            || false,
        )
        .is_err()
    );
    let second = release_payload(
        &result.model_after,
        &result.state_after,
        &result.attachments[1],
        None,
        &[],
        settings,
        || false,
    )
    .unwrap();
    let first = release_payload(
        &second.model_after,
        &second.state_after,
        &result.attachments[0],
        None,
        &[],
        settings,
        || false,
    )
    .unwrap();
    assert_eq!(first.model_after, model);
    assert_eq!(first.state_after, result.state_after);
    close(first.supplied_work_j, 0.);
    close(second.supplied_work_j, 0.);
    assert_eq!(first.payload_after.inertia, inputs[0].payload.inertia);
    assert_eq!(second.payload_after.inertia, inputs[1].payload.inertia);
    assert!(
        release_payload(
            &first.model_after,
            &first.state_after,
            &result.attachments[0],
            None,
            &[],
            settings,
            || false,
        )
        .is_err()
    );
}

fn spatial_inputs(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
) -> Vec<PayloadCaptureInput> {
    let kin = kinematics(model, state);
    [3, 1, 3, 0]
        .into_iter()
        .enumerate()
        .map(|(i, target)| {
            let f = i as f64;
            let mount = pose(
                Vec3::new(0.25 + 0.08 * f, -0.17, 0.31),
                Vec3::new(0.28, -0.19 + 0.1 * f, 0.23),
            );
            PayloadCaptureInput {
                target_link: target,
                payload: PayloadState {
                    inertia: inertia(
                        0.4 + 0.1 * f,
                        Vec3::new(0.08, -0.06, 0.11),
                        Vec3::new(0.18 + 0.02 * f, 0.24, 0.3),
                        Vec3::new(-0.2, 0.3 + 0.04 * f, 0.15),
                    ),
                    world_from_body: kin.world_from_link[target].compose(mount).unwrap(),
                    twist_body: Twist::new(
                        Vec3::new(-0.3 + 0.12 * f, 0.55, 0.2),
                        Vec3::new(0.7, -0.4 + 0.15 * f, 0.15),
                    ),
                },
            }
        })
        .collect()
}

#[test]
fn g3_batch_full_tensors_preserve_mixed_joint_momentum_under_permutation_and_scene_changes() {
    let (model, state) = serial_tree();
    let inputs = spatial_inputs(&model, &state);
    let settings = controls(&[2e-10; 3], 40.);
    let result = capture_payloads(&model, &state, &inputs, &[], settings, || false).unwrap();
    check_batch(&model, &state, &inputs, &[], &result);
    check_loaded_inertias(&model, &state, &inputs, &result);
    let reversed: Vec<_> = inputs.iter().rev().copied().collect();
    let swapped = capture_payloads(&model, &state, &reversed, &[], settings, || false).unwrap();
    check_batch(&model, &state, &reversed, &[], &swapped);
    check_loaded_inertias(&model, &state, &reversed, &swapped);
    vector(
        swapped.state_after.base.twist_body.angular,
        result.state_after.base.twist_body.angular,
    );
    vector(
        swapped.state_after.base.twist_body.linear,
        result.state_after.base.twist_body.linear,
    );
    for (&a, &b) in swapped
        .state_after
        .joint_velocity
        .iter()
        .zip(&result.state_after.joint_velocity)
    {
        close(a, b);
    }
    close(
        swapped.removed_relative_energy_j,
        result.removed_relative_energy_j,
    );
    for (a, b) in swapped
        .payloads_after
        .iter()
        .rev()
        .zip(&result.payloads_after)
    {
        vector(a.twist_body.angular, b.twist_body.angular);
        vector(a.twist_body.linear, b.twist_body.linear);
    }
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
        transformed.base.twist_body.linear = transformed.base.twist_body.linear + base_boost;
        let incoming: Vec<_> = inputs
            .iter()
            .copied()
            .map(|mut input| {
                input.payload.world_from_body =
                    scene.compose(input.payload.world_from_body).unwrap();
                input.payload.twist_body.linear = input.payload.twist_body.linear
                    + input
                        .payload
                        .world_from_body
                        .rotation()
                        .inverse()
                        .rotate(boost)
                        .unwrap();
                input
            })
            .collect();
        let shifted =
            capture_payloads(&model, &transformed, &incoming, &[], settings, || false).unwrap();
        check_batch(&model, &transformed, &incoming, &[], &shifted);
        check_loaded_inertias(&model, &transformed, &incoming, &shifted);
        vector(
            shifted.state_after.base.twist_body.angular,
            result.state_after.base.twist_body.angular,
        );
        vector(
            shifted.state_after.base.twist_body.linear - base_boost,
            result.state_after.base.twist_body.linear,
        );
        for (&a, &b) in shifted
            .state_after
            .joint_velocity
            .iter()
            .zip(&result.state_after.joint_velocity)
        {
            close(a, b);
        }
        close(
            shifted.removed_relative_energy_j,
            result.removed_relative_energy_j,
        );
    }
}

#[test]
fn g0_g1_batch_singleton_and_existing_gears_use_the_same_native_endpoint() {
    let (model, state, payload) = coax(None, 1.2);
    let settings = controls(&[2e-10], 20.);
    let single = capture_payload(&model, &state, 1, payload, &[], settings, || false).unwrap();
    let batch = capture_payloads(
        &model,
        &state,
        &[PayloadCaptureInput {
            target_link: 1,
            payload,
        }],
        &[],
        settings,
        || false,
    )
    .unwrap();
    assert_eq!(single.model_after, batch.model_after);
    assert_eq!(single.state_after, batch.state_after);
    assert_eq!(single.attachment, batch.attachments[0]);
    assert_eq!(single.payload_after, batch.payloads_after[0]);
    assert_eq!(single.balance, batch.balance);
    assert_eq!(
        single.removed_relative_energy_j,
        batch.removed_relative_energy_j
    );

    let limits = Some(JointLimits::new(-2., 2., 0.5, 100.).unwrap());
    let (model, state, _) = coax(limits, 0.);
    let mut inputs = requests();
    inputs[1].payload.twist_body.angular.z = 4.;
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![1.],
        position: 0.,
        position_tolerance: 2e-10,
        velocity_tolerance: 2e-10,
    }];
    assert!(capture_payloads(&model, &state, &inputs, &[], settings, || false).is_err());
    let result = capture_payloads(&model, &state, &inputs, &rows, settings, || false).unwrap();
    check_batch(&model, &state, &inputs, &rows, &result);
    vector(
        result.state_after.base.twist_body.angular,
        Vec3::new(0., 0., 41. / 20.),
    );
    close(result.state_after.joint_velocity[0], 0.);
    close(result.balance.multipliers[0], -123. / 20.);
    close(result.balance.reaction_work_j, 0.);
    close(result.removed_relative_energy_j, 17.359375);
    let duplicates = [rows[0].clone(), rows[0].clone()];
    assert!(capture_payloads(&model, &state, &inputs, &duplicates, settings, || false).is_err());
}

#[test]
fn g4_batch_invalid_late_inventory_work_and_terminal_cancellation_publish_nothing() {
    let (model, state) = serial_tree();
    let inputs = spatial_inputs(&model, &state);
    let snapshot = (model.clone(), state.clone(), inputs.clone());
    let settings = controls(&[2e-10; 3], 40.);
    let mut polls = 0;
    let expected = capture_payloads(&model, &state, &inputs, &[], settings, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 8);
    let mut observed = 0;
    assert!(matches!(
        capture_payloads(&model, &state, &inputs, &[], settings, || {
            observed += 1;
            observed == polls
        }),
        Err(PayloadError::Cancelled)
    ));
    assert_eq!(observed, polls);
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &inputs,
            &[],
            controls(&[2e-10; 3], 0.),
            || false
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    for invalid in [
        PayloadCaptureInput {
            target_link: model.link_count(),
            ..inputs[3]
        },
        PayloadCaptureInput {
            payload: PayloadState {
                twist_body: Twist::new(Vec3::new(f64::NAN, 0., 0.), Vec3::default()),
                ..inputs[3].payload
            },
            ..inputs[3]
        },
        PayloadCaptureInput {
            payload: PayloadState {
                inertia: sphere(1e-30, 1e-30),
                ..inputs[3].payload
            },
            ..inputs[3]
        },
    ] {
        let mut bad = inputs.clone();
        bad[3] = invalid;
        assert!(capture_payloads(&model, &state, &bad, &[], settings, || false).is_err());
    }
    assert!(capture_payloads(&model, &state, &[], &[], settings, || false).is_err());
    assert_eq!(
        (&model, &state, &inputs),
        (&snapshot.0, &snapshot.1, &snapshot.2)
    );
    let retry = capture_payloads(&model, &state, &inputs, &[], settings, || false).unwrap();
    assert_eq!(retry, expected);
}
