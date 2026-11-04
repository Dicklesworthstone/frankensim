//! G0/G1/G3/G4: simultaneous incoming bodies retain the original mixed mechanism.

use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::inventory::bilateral::release::release_payload;
use fs_mbd::articulated::motion::inventory::bilateral::{capture_payload, capture_payloads};
use fs_mbd::articulated::motion::inventory::{
    self, PayloadCaptureInput, PayloadControls, PayloadError, PayloadState,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopError, DistanceLoop, joint,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link, SpatialInertia,
};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn close(actual: f64, expected: f64) {
    near(actual, expected, 3e-9 * (1.0 + expected.abs()));
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
                maximum_speed.map(|speed| JointLimits::new(-10., 10., speed, 0.1).unwrap()),
            )
            .unwrap(),
            sphere(mass),
        ));
    }
    ArticulatedModel::new(links).unwrap()
}
fn initial() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.7, 0., 0.)),
        ),
        joint_position: vec![1., 3., 1.5],
        joint_velocity: vec![0.4, 0.4, 0.2],
    }
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
fn input(target_link: usize, mass: f64, position: f64, velocity: f64) -> PayloadCaptureInput {
    PayloadCaptureInput {
        target_link,
        payload: PayloadState {
            inertia: sphere(mass),
            world_from_body: Se3::from_parts(So3::identity(), Vec3::new(position, 0., 0.)).unwrap(),
            twist_body: Twist::new(Vec3::default(), Vec3::new(velocity, 0., 0.)),
        },
    }
}
fn incoming() -> [PayloadCaptureInput; 2] {
    [input(1, 2., 1., 2.3), input(3, 1., 1.5, -0.4)]
}
fn controls(tolerances: &[f64]) -> PayloadControls<'_> {
    PayloadControls {
        relative_momentum_tolerance: 2e-12,
        angular_momentum_tolerance_nms: 2e-8,
        linear_momentum_tolerance_ns: 2e-8,
        joint_momentum_tolerances: tolerances,
        energy_tolerance_j: 2e-8,
        maximum_energy_exchange_j: 20.,
        maximum_constraint_work_j: 0.1,
        minimum_reciprocal_condition: 1e-10,
    }
}

#[test]
fn simultaneous_capture_transmits_shared_recoil_then_continues_loaded_motion() {
    let model = sliders(Some(0.8));
    let start = initial();
    let inputs = incoming();
    // Without the material closure, even the gear-constrained first slider
    // reaches 1.078. The combined physical endpoint is 0.786 and is admissible.
    assert!(
        inventory::capture_payloads(
            &model,
            &start,
            &inputs,
            &rows(),
            controls(&[2e-8; 3]),
            || false,
        )
        .is_err()
    );
    let result = capture_payloads(
        &model,
        &start,
        &inputs,
        &loops(),
        &rows(),
        controls(&[2e-8; 3]),
        Default::default(),
        || false,
    )
    .unwrap();
    let capture = &result.capture;
    // Independent reduced coordinates V and u have qdot=(2,2,1)*u.
    // All incoming bodies give M=[[17,19],[19,33]], p=(16.8,23.4).
    let velocity = (33. * 16.8 - 19. * 23.4) / 200.;
    let mode = (17. * 23.4 - 19. * 16.8) / 200.;
    close(velocity, 0.549);
    close(mode, 0.393);
    close(capture.state_after.base.twist_body.linear.x, velocity);
    for (rate, factor) in capture.state_after.joint_velocity.iter().zip([2., 2., 1.]) {
        close(*rate, factor * mode);
    }
    assert_eq!(capture.state_after.joint_position, start.joint_position);
    assert_eq!(
        capture.state_after.base.world_from_base,
        start.base.world_from_base
    );
    close(capture.model_after.links()[1].inertia().mass(), 4.);
    close(capture.model_after.links()[3].inertia().mass(), 5.);
    close(
        capture.payloads_after[0].twist_body.linear.x,
        velocity + 2. * mode,
    );
    close(
        capture.payloads_after[1].twist_body.linear.x,
        velocity + mode,
    );
    close(capture.balance.momentum_before.kinetic_energy, 11.24);
    close(
        capture.balance.momentum_after.kinetic_energy,
        0.5 * (16.8 * velocity + 23.4 * mode),
    );
    close(capture.removed_relative_energy_j, 2.0303);
    close(capture.balance.momentum_after.linear_momentum_world.x, 16.8);
    close(capture.balance.energy_residual_j, 0.);
    close(capture.balance.reaction_work_j, 0.);
    let reaction = result.reaction.as_ref().unwrap();
    close(reaction.loop_multipliers_ns[0], 1.46);
    close(reaction.joint_multipliers[0], -0.755);

    let reversed = capture_payloads(
        &model,
        &start,
        &[inputs[1], inputs[0]],
        &loops(),
        &rows(),
        controls(&[2e-8; 3]),
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(reversed.capture.model_after, capture.model_after);
    close(
        reversed.capture.state_after.base.twist_body.linear.x,
        velocity,
    );
    close(reversed.capture.state_after.joint_velocity[2], mode);
    close(
        reversed.capture.removed_relative_energy_j,
        capture.removed_relative_energy_j,
    );
    assert_eq!(reversed.capture.attachments[0], capture.attachments[1]);
    assert_eq!(reversed.capture.attachments[1], capture.attachments[0]);

    let zero = [Wrench::default(); 4];
    let moved = joint::step_with_bushings(
        &capture.model_after,
        &capture.state_after,
        SpatialMotionLoads {
            held: FreeFloatingMotionLoads {
                generalized_force: &[0.; 3],
                external_body_wrench: &zero,
                gravity_world: Vec3::default(),
            },
            world_points: &[],
            connections: &[],
        },
        &loops(),
        &rows(),
        &[],
        0.04,
        Default::default(),
        || false,
    )
    .unwrap();
    for (i, factor) in [2., 2., 1.].into_iter().enumerate() {
        close(moved.bilateral.state_after.joint_velocity[i], factor * mode);
        close(
            moved.bilateral.state_after.joint_position[i],
            start.joint_position[i] + 0.04 * factor * mode,
        );
    }
    close(moved.bilateral.mechanical_work_residual, 0.);
    assert_eq!(start, initial());
    assert_eq!(inputs, incoming());
}

#[test]
fn opposing_same_link_arrivals_preserve_finite_slip_and_release_stack_identity() {
    let model = sliders(None);
    let mut state = initial();
    state.base.twist_body.linear.x = 0.5;
    let defect = 1. / 1024.;
    state.joint_velocity = vec![0.5, 0.5, 0.25 - 0.5 * defect];
    let mut rows = rows();
    rows[0].velocity_tolerance = 2. * defect;
    // Each body has a nonzero mismatch, but their impulses cancel exactly on
    // the same target. The event removes 8 J while the mechanism remains still
    // in velocity space; the existing finite gear defect must not be projected.
    let inputs = [input(1, 2., 1., 3.), input(1, 2., 1., -1.)];
    let result = capture_payloads(
        &model,
        &state,
        &inputs,
        &loops(),
        &rows,
        controls(&[2e-8; 3]),
        Default::default(),
        || false,
    )
    .unwrap();
    let capture = &result.capture;
    assert_eq!(capture.state_after, state);
    close(capture.removed_relative_energy_j, 8.);
    close(capture.balance.energy_residual_j, 0.);
    close(capture.balance.reaction_work_j, 0.);
    close(capture.model_after.links()[1].inertia().mass(), 6.);
    close(capture.attachments[0].retained_inertia().mass(), 2.);
    close(capture.attachments[1].retained_inertia().mass(), 4.);
    let reaction = result.reaction.as_ref().unwrap();
    assert_eq!(reaction.loop_multipliers_ns, vec![0.]);
    assert_eq!(reaction.joint_multipliers, vec![0.]);
    close(reaction.joint_velocity_residuals_after[0], defect);

    assert!(matches!(
        release_payload(
            &capture.model_after,
            &capture.state_after,
            &capture.attachments[0],
            None,
            &loops(),
            &rows,
            controls(&[2e-8; 3]),
            Default::default(),
            || false,
        ),
        Err(PayloadError::Invalid(
            "release attachment inventory changed"
        ))
    ));
    let last = release_payload(
        &capture.model_after,
        &capture.state_after,
        &capture.attachments[1],
        None,
        &loops(),
        &rows,
        controls(&[2e-8; 3]),
        Default::default(),
        || false,
    )
    .unwrap()
    .release;
    close(last.model_after.links()[1].inertia().mass(), 4.);
    assert_eq!(last.state_after, state);
    let first = release_payload(
        &last.model_after,
        &last.state_after,
        &capture.attachments[0],
        None,
        &loops(),
        &rows,
        controls(&[2e-8; 3]),
        Default::default(),
        || false,
    )
    .unwrap()
    .release;
    assert_eq!(first.model_after, model);
    assert_eq!(first.state_after, state);
    close(first.payload_after.twist_body.linear.x, 1.);
    close(last.payload_after.twist_body.linear.x, 1.);
    close(first.supplied_work_j, 0.);
    close(last.supplied_work_j, 0.);
}

#[test]
fn finite_constraint_work_and_large_boost_use_the_complete_incoming_inventory() {
    let model = sliders(None);
    let defect = 1e-3;
    let mut rows = rows();
    rows[0].velocity_tolerance = 2e-3;
    for boost in [0., 1e9] {
        let mut state = initial();
        state.base.twist_body.linear.x += boost;
        state.joint_velocity[2] -= 0.5 * defect;
        let mut inputs = incoming();
        for input in &mut inputs {
            input.payload.twist_body.linear.x += boost;
        }
        let tolerances = [2e-7; 3];
        let controls = PayloadControls {
            energy_tolerance_j: 4e-6,
            ..controls(&tolerances)
        };
        let result = capture_payloads(
            &model,
            &state,
            &inputs,
            &loops(),
            &rows,
            controls,
            Default::default(),
            || false,
        )
        .unwrap();
        // Independent solve in the original base's translating frame. Build
        // its data from the actually representable incoming velocities.
        let v = &state.joint_velocity;
        let incoming_1 = inputs[0].payload.twist_body.linear.x - state.base.twist_body.linear.x;
        let incoming_3 = inputs[1].payload.twist_body.linear.x - state.base.twist_body.linear.x;
        let p = 2. * v[0] + 3. * v[1] + 4. * v[2] + 2. * incoming_1 + incoming_3;
        let h = 4. * v[0] + 6. * v[1] + 4. * v[2] + 4. * incoming_1 + incoming_3;
        let velocity = (33. * p - 19. * h) / 200.;
        let mode = (17. * h - 19. * p) / 200.;
        let kinetic_before = v[0] * v[0]
            + 1.5 * v[1] * v[1]
            + 2. * v[2] * v[2]
            + incoming_1 * incoming_1
            + 0.5 * incoming_3 * incoming_3;
        let kinetic_after = 0.5 * (p * velocity + h * mode);
        let lambda = 0.5 * (4. * v[2] + incoming_3 - 5. * (velocity + mode));
        near(
            result.capture.state_after.base.twist_body.linear.x - state.base.twist_body.linear.x,
            velocity,
            3e-7,
        );
        near(result.capture.state_after.joint_velocity[2], mode, 3e-7);
        near(
            result.capture.removed_relative_energy_j,
            kinetic_before - kinetic_after,
            4e-6,
        );
        near(
            result.reaction.as_ref().unwrap().joint_multipliers[0],
            lambda,
            3e-7,
        );
        near(
            result.capture.balance.reaction_work_j,
            0.5 * lambda * defect,
            2e-9,
        );
        assert!(result.capture.balance.reaction_work_j < -3e-4);
        assert!(result.capture.balance.energy_residual_j.abs() <= controls.energy_tolerance_j);
        assert!(matches!(
            capture_payloads(
                &model,
                &state,
                &inputs,
                &loops(),
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

#[test]
fn singleton_empty_and_atomic_refusals_share_the_original_batch_owner() {
    let model = sliders(None);
    let state = initial();
    let inputs = incoming();
    let tolerances = [2e-8; 3];
    let controls = controls(&tolerances);
    let legacy =
        inventory::capture_payloads(&model, &state, &inputs, &rows(), controls, || false).unwrap();
    let empty = capture_payloads(
        &model,
        &state,
        &inputs,
        &[],
        &rows(),
        controls,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(empty.capture, legacy);
    assert!(empty.reaction.is_none());
    let singleton = capture_payload(
        &model,
        &state,
        inputs[0].target_link,
        inputs[0].payload,
        &loops(),
        &rows(),
        controls,
        Default::default(),
        || false,
    )
    .unwrap();
    let one = capture_payloads(
        &model,
        &state,
        &inputs[..1],
        &loops(),
        &rows(),
        controls,
        Default::default(),
        || false,
    )
    .unwrap();
    assert_eq!(one.capture.model_after, singleton.capture.model_after);
    assert_eq!(one.capture.state_after, singleton.capture.state_after);
    assert_eq!(one.capture.attachments[0], singleton.capture.attachment);
    assert_eq!(
        one.capture.payloads_after[0],
        singleton.capture.payload_after
    );
    assert_eq!(one.capture.balance, singleton.capture.balance);
    assert_eq!(one.reaction, singleton.reaction);

    let redundant = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![-1., 1., 0.],
        position: 2.,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }];
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &inputs,
            &loops(),
            &redundant,
            controls,
            Default::default(),
            || false,
        ),
        Err(PayloadError::Bilateral(BilateralLoopError::Mobility { .. }))
    ));
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &inputs,
            &loops(),
            &rows(),
            PayloadControls {
                minimum_reciprocal_condition: 0.99,
                ..controls
            },
            Default::default(),
            || false,
        ),
        Err(PayloadError::Constraint(
            ConstraintError::IllConditioned { .. }
        ))
    ));
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &inputs,
            &loops(),
            &rows(),
            PayloadControls {
                maximum_energy_exchange_j: 1.,
                ..controls
            },
            Default::default(),
            || false,
        ),
        Err(PayloadError::EnergyExchange { .. })
    ));
    assert!(
        capture_payloads(
            &sliders(Some(0.78)),
            &state,
            &inputs,
            &loops(),
            &rows(),
            controls,
            Default::default(),
            || false,
        )
        .is_err()
    );
    let invalid = [
        inputs[0],
        PayloadCaptureInput {
            target_link: usize::MAX,
            ..inputs[1]
        },
    ];
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &invalid,
            &loops(),
            &rows(),
            controls,
            Default::default(),
            || false,
        ),
        Err(PayloadError::Invalid("payload target link index"))
    ));
    assert!(matches!(
        capture_payloads(
            &model,
            &state,
            &[],
            &loops(),
            &rows(),
            controls,
            Default::default(),
            || false,
        ),
        Err(PayloadError::Invalid("empty payload capture"))
    ));
    let polls = Cell::new(0);
    let expected = capture_payloads(
        &model,
        &state,
        &inputs,
        &loops(),
        &rows(),
        controls,
        Default::default(),
        || {
            polls.set(polls.get() + 1);
            false
        },
    )
    .unwrap();
    let total = polls.get();
    assert!(total > 4);
    for limit in [1, total] {
        let calls = Cell::new(0);
        assert_eq!(
            capture_payloads(
                &model,
                &state,
                &inputs,
                &loops(),
                &rows(),
                controls,
                Default::default(),
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= limit
                },
            ),
            Err(PayloadError::Cancelled)
        );
    }
    assert_eq!(
        capture_payloads(
            &model,
            &state,
            &inputs,
            &loops(),
            &rows(),
            controls,
            Default::default(),
            || false,
        )
        .unwrap(),
        expected
    );
    assert_eq!(state, initial());
    assert_eq!(inputs, incoming());
}
