//! Stiff finite-body drag against independent exponential/work solutions.

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};
use fs_time::lie::viscous::pair::{ViscousPairStep, ViscousPairStepper};

fn diagonal(x: f64, y: f64, z: f64) -> Mat3 {
    Mat3 {
        m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
    }
}

fn sphere(moment: f64) -> InertiaTensor3 {
    InertiaTensor3::from_principal(Vec3::new(moment, moment, moment), So3::identity()).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

fn balances(result: ViscousPairStep, tolerance: f64) {
    assert!(result.dissipated_energy_j >= 0.0);
    assert_eq!(
        result.bodies[0].coupling_angular_impulse_space
            + result.bodies[1].coupling_angular_impulse_space,
        Vec3::default()
    );
    assert!(result.angular_momentum_residual_space.norm() < tolerance);
    close(result.work_energy_residual_j, 0.0, tolerance);
    close(
        result.bodies[0].coupling_work_j + result.bodies[1].coupling_work_j,
        -result.dissipated_energy_j,
        tolerance,
    );
    for body in result.bodies {
        assert!(body.angular_momentum_residual_space.norm() < tolerance);
        close(body.work_energy_residual_j, 0.0, tolerance);
    }
}

#[test]
fn stiff_pair_recoils_to_a_finite_common_spin_without_midpoint_sign_reversal() {
    let inertia = [sphere(2.0), sphere(5.0)];
    let rotations = [So3::identity(); 2];
    let momenta = [Vec3::new(0.0, 0.0, 4.0), Vec3::new(0.0, 0.0, -5.0)];
    let drives = [HeldTorque::Body(Vec3::default()); 2];
    let damping = diagonal(0.0, 0.0, 100.0);
    let dt = 0.2;
    assert!(
        ViscousPairStepper::new(inertia, damping, dt)
            .unwrap()
            .step(rotations, momenta, drives)
            .is_err()
    );
    let result = ViscousPairStepper::new_exponential(inertia, damping, dt)
        .unwrap()
        .step(rotations, momenta, drives)
        .unwrap();
    let lambda = 100.0 * (1.0 / 2.0 + 1.0 / 5.0);
    let slip = 3.0 * (-lambda * dt).exp();
    let mean = -1.0 / 7.0;
    close(
        result.bodies[0].angular_velocity_after_body.z,
        mean + 5.0 / 7.0 * slip,
        2e-14,
    );
    close(
        result.bodies[1].angular_velocity_after_body.z,
        mean - 2.0 / 7.0 * slip,
        2e-14,
    );
    assert!(
        result.bodies[0].angular_velocity_after_body.z
            > result.bodies[1].angular_velocity_after_body.z
    );
    close(
        result.dissipated_energy_j,
        0.5 / 0.7 * 9.0 * (1.0 - (-2.0 * lambda * dt).exp()),
        5e-14,
    );
    balances(result, 8e-14);

    let saturated = ViscousPairStepper::new_exponential(inertia, diagonal(0.0, 0.0, 1e6), dt)
        .unwrap()
        .step(rotations, momenta, drives)
        .unwrap();
    for body in saturated.bodies {
        close(body.angular_velocity_after_body.z, mean, 2e-14);
    }
    close(saturated.dissipated_energy_j, 0.5 / 0.7 * 9.0, 5e-14);
    balances(saturated, 8e-14);
}

#[test]
fn independent_held_drives_have_exact_source_work_and_physical_slip_heat() {
    let inertia = [sphere(2.0), sphere(5.0)];
    let initial = [2.1, -0.4];
    let torque = [0.7, -0.2];
    let damping: f64 = 100.0;
    let duration: f64 = 0.4;
    let lambda = damping * 0.7;
    let terminal = (torque[0] / 2.0 - torque[1] / 5.0) / lambda;
    let amplitude = initial[0] - initial[1] - terminal;
    let decay = (-lambda * duration).exp();
    let relative = terminal + amplitude * decay;
    let relative_angle = terminal * duration + amplitude * (1.0 - decay) / lambda;
    let mean_initial = (2.0 * initial[0] + 5.0 * initial[1]) / 7.0;
    let mean = mean_initial + (torque[0] + torque[1]) * duration / 7.0;
    let mean_angle =
        mean_initial * duration + 0.5 * (torque[0] + torque[1]) / 7.0 * duration * duration;
    let rates = [mean + 5.0 / 7.0 * relative, mean - 2.0 / 7.0 * relative];
    let angles = [
        mean_angle + 5.0 / 7.0 * relative_angle,
        mean_angle - 2.0 / 7.0 * relative_angle,
    ];
    let heat = damping
        * (terminal * terminal * duration
            + 2.0 * terminal * amplitude * (1.0 - decay) / lambda
            + amplitude * amplitude * (1.0 - decay * decay) / (2.0 * lambda));
    let result =
        ViscousPairStepper::new_exponential(inertia, diagonal(0.0, 0.0, damping), duration)
            .unwrap()
            .step(
                [So3::identity(); 2],
                [
                    Vec3::new(0.0, 0.0, 2.0 * initial[0]),
                    Vec3::new(0.0, 0.0, 5.0 * initial[1]),
                ],
                torque.map(|value| HeldTorque::Body(Vec3::new(0.0, 0.0, value))),
            )
            .unwrap();
    for i in 0..2 {
        close(
            result.bodies[i].angular_velocity_after_body.z,
            rates[i],
            2e-14,
        );
        close(result.bodies[i].drive_work_j, torque[i] * angles[i], 2e-14);
    }
    close(result.dissipated_energy_j, heat, 6e-14);
    balances(result, 8e-14);
}

#[test]
fn stiff_full_inertias_and_off_axis_drag_preserve_the_analytic_principal_flow() {
    let axes = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.3, 0.2))).unwrap();
    let axis = axes.rotate(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let inertia = [Vec3::new(2.0, 3.0, 4.0), Vec3::new(3.0, 4.0, 5.0)]
        .map(|moments| InertiaTensor3::from_principal(moments, axes).unwrap());
    let damping = axes
        .matrix()
        .compose(diagonal(1.0, 5.0, 200.0))
        .compose(axes.matrix().transpose());
    let result = ViscousPairStepper::new_exponential(inertia, damping, 0.2)
        .unwrap()
        .step(
            [So3::identity(); 2],
            [axis.scale(8.0), axis.scale(-5.0)],
            [HeldTorque::Body(Vec3::default()); 2],
        )
        .unwrap();
    let slip = 3.0 * (-18.0_f64).exp();
    let mean = 1.0 / 3.0;
    let expected = [
        axis.scale(4.0 * (mean + 5.0 / 9.0 * slip)),
        axis.scale(5.0 * (mean - 4.0 / 9.0 * slip)),
    ];
    for i in 0..2 {
        assert!((result.bodies[i].momentum_after_body - expected[i]).norm() < 3e-12);
    }
    close(
        result.dissipated_energy_j,
        0.5 / 0.45 * 9.0 * (1.0 - (-36.0_f64).exp()),
        3e-12,
    );
    balances(result, 4e-12);
}

#[test]
fn dominant_common_spin_does_not_erase_small_reciprocal_heat() {
    let inertia = [sphere(2.0); 2];
    let initial = [Vec3::new(2.0, 0.0, 1e20), Vec3::new(0.0, 0.0, 1e20)];
    let result = ViscousPairStepper::new_exponential(inertia, diagonal(2e40, 0.0, 0.0), 1e-40)
        .unwrap()
        .step(
            [So3::identity(); 2],
            initial,
            [HeldTorque::Body(Vec3::default()); 2],
        )
        .unwrap();
    close(
        result.bodies[0].momentum_after_body.x,
        1.0 + (-2.0_f64).exp(),
        2e-14,
    );
    close(
        result.bodies[1].momentum_after_body.x,
        1.0 - (-2.0_f64).exp(),
        2e-14,
    );
    for i in 0..2 {
        assert_eq!(result.bodies[i].momentum_after_body.z, initial[i].z);
    }
    close(
        result.dissipated_energy_j,
        0.5 * (1.0 - (-4.0_f64).exp()),
        2e-14,
    );
    balances(result, 3e-14);
}

#[test]
fn weak_and_zero_drag_keep_resolved_impulse_and_native_limits() {
    let inertia = [sphere(2.0); 2];
    let rotations = [So3::identity(); 2];
    let momenta = [Vec3::new(0.0, 0.0, 2.0), Vec3::default()];
    let zero = [HeldTorque::Body(Vec3::default()); 2];
    let weak = ViscousPairStepper::new_exponential(inertia, diagonal(0.0, 0.0, 2e-12), 0.1)
        .unwrap()
        .step(rotations, momenta, zero)
        .unwrap();
    close(
        weak.bodies[1].momentum_after_body.z,
        -(-2e-13_f64).exp_m1(),
        2e-27,
    );
    close(
        weak.dissipated_energy_j,
        -0.5 * (-4e-13_f64).exp_m1(),
        2e-27,
    );
    balances(weak, 2e-15);
    let drives = [
        HeldTorque::Space(Vec3::new(0.1, 0.2, -0.1)),
        HeldTorque::Body(Vec3::new(-0.2, 0.1, 0.3)),
    ];
    let plain = ViscousPairStepper::new(inertia, Mat3::zero(), 0.1)
        .unwrap()
        .step(rotations, momenta, drives)
        .unwrap();
    let exponential = ViscousPairStepper::new_exponential(inertia, Mat3::zero(), 0.1)
        .unwrap()
        .step(rotations, momenta, drives)
        .unwrap();
    assert_eq!(plain, exponential);
    assert!(ViscousPairStepper::new_exponential(inertia, diagonal(-1.0, 0.0, 0.0), 0.1).is_err());
    assert!(ViscousPairStepper::new_exponential(inertia, diagonal(1.0, 0.0, 0.0), 0.0).is_err());
}

#[test]
fn inertia_weighted_heat_does_not_require_representable_raw_spin_squares() {
    let expected_heat = 0.5 / 0.7 * 9.0 * (1.0 - (-28.0_f64).exp());
    let expected_slip = 3.0 * (-14.0_f64).exp();
    for (mass_scale, speed_scale, time_scale) in [
        (1e240_f64, 1e-160_f64, 1.0_f64),
        (1e-240_f64, 1e160_f64, 1e-160_f64),
    ] {
        let inertia = [sphere(2.0 * mass_scale), sphere(5.0 * mass_scale)];
        let step = ViscousPairStepper::new_exponential(
            inertia,
            diagonal(0.0, 0.0, 100.0 * (mass_scale / time_scale)),
            0.2 * time_scale,
        )
        .unwrap()
        .step(
            [So3::identity(); 2],
            [
                Vec3::new(0.0, 0.0, 4.0 * (mass_scale * speed_scale)),
                Vec3::new(0.0, 0.0, -5.0 * (mass_scale * speed_scale)),
            ],
            [HeldTorque::Body(Vec3::default()); 2],
        )
        .unwrap();
        let energy_scale = (mass_scale.sqrt() * speed_scale).powi(2);
        close(
            step.dissipated_energy_j / energy_scale,
            expected_heat,
            2e-12,
        );
        close(step.work_energy_residual_j / energy_scale, 0.0, 2e-12);
        close(
            step.bodies[0].angular_velocity_after_body.z / speed_scale,
            -1.0 / 7.0 + 5.0 / 7.0 * expected_slip,
            2e-13,
        );
        close(
            step.bodies[1].angular_velocity_after_body.z / speed_scale,
            -1.0 / 7.0 - 2.0 / 7.0 * expected_slip,
            2e-13,
        );
    }
}

#[test]
fn stiff_counterrotation_retains_exponential_tails_and_time_partition() {
    // Total momentum is exactly zero. Both finite bodies must keep their
    // opposite exponentially small momenta rather than rounding a nearly
    // complete reciprocal impulse to artificial rest or a residual spin.
    let duration = 0.1;
    for moments in [[2.0_f64, 2.0_f64], [2.0, 5.0]] {
        let inertia = moments.map(sphere);
        let mobility = moments[0].recip() + moments[1].recip();
        for relaxation in [20.0_f64, 80.0, 200.0, 700.0] {
            let damping = diagonal(0.0, 0.0, relaxation / (duration * mobility));
            let whole = ViscousPairStepper::new_exponential(inertia, damping, duration).unwrap();
            let part =
                ViscousPairStepper::new_exponential(inertia, damping, duration / 8.0).unwrap();
            for sign in [-1.0, 1.0] {
                let initial = [
                    Vec3::new(0.0, 0.0, sign * 6.0),
                    Vec3::new(0.0, 0.0, -sign * 6.0),
                ];
                let drives = [HeldTorque::Body(Vec3::default()); 2];
                let result = whole.step([So3::identity(); 2], initial, drives).unwrap();
                let expected = sign * 6.0 * (-relaxation).exp();
                for (body, target) in [expected, -expected].into_iter().enumerate() {
                    close(
                        result.bodies[body].momentum_after_body.z,
                        target,
                        target.abs() * 3e-12,
                    );
                }
                close(
                    result.dissipated_energy_j,
                    18.0 * mobility * -(-2.0 * relaxation).exp_m1(),
                    2e-12,
                );
                balances(result, 2e-12);
                let mut rotations = [So3::identity(); 2];
                let mut momenta = initial;
                for _ in 0..8 {
                    let split = part.step(rotations, momenta, drives).unwrap();
                    rotations = split.bodies.map(|body| body.rotation_after);
                    momenta = split.bodies.map(|body| body.momentum_after_body);
                }
                close(momenta[0].z, expected, expected.abs() * 1e-11);
                close(momenta[1].z, -expected, expected.abs() * 1e-11);
            }
        }
    }
}

#[test]
fn stiff_opposed_drives_keep_the_small_terminal_momentum_from_rest() {
    // The constant source impulse can nearly cancel the integrated drag
    // even when the initial momentum is zero. The endpoint must then use
    // phi_1 directly, not the difference of those much larger impulses.
    let duration = 0.1;
    let moments = [2.0_f64, 5.0_f64];
    let mobility = moments[0].recip() + moments[1].recip();
    for damping in [1e6_f64, 1e20] {
        let rate = damping * mobility;
        let stepper = ViscousPairStepper::new_exponential(
            moments.map(sphere),
            diagonal(0.0, 0.0, damping),
            duration,
        )
        .unwrap();
        for drive in [-1e6_f64, -1e-6, 1e-6, 1e6] {
            let result = stepper
                .step(
                    [So3::identity(); 2],
                    [Vec3::default(); 2],
                    [
                        HeldTorque::Body(Vec3::new(0.0, 0.0, drive)),
                        HeldTorque::Body(Vec3::new(0.0, 0.0, -drive)),
                    ],
                )
                .unwrap();
            let expected = (drive / rate) * -(-rate * duration).exp_m1();
            close(
                result.bodies[0].momentum_after_body.z,
                expected,
                expected.abs() * 3e-12,
            );
            close(
                result.bodies[1].momentum_after_body.z,
                -expected,
                expected.abs() * 3e-12,
            );
            assert_eq!(
                result.bodies[0].momentum_after_body + result.bodies[1].momentum_after_body,
                Vec3::default(),
            );
            let fraction = -(-rate * duration).exp_m1();
            let integrated_momentum = (drive / rate) * (duration - fraction / rate);
            let heat = drive * drive / damping
                * (duration - (2.0 * fraction + 0.5 * (-2.0 * rate * duration).exp_m1()) / rate);
            close(result.dissipated_energy_j, heat, heat.abs() * 3e-12);
            close(result.work_energy_residual_j, 0.0, heat.abs() * 3e-12);
            for (body, moment) in result.bodies.into_iter().zip(moments) {
                let work = (drive / moment) * integrated_momentum;
                close(body.drive_work_j, work, work.abs() * 3e-12);
                close(body.work_energy_residual_j, 0.0, work.abs() * 3e-12);
                assert!(
                    body.angular_momentum_residual_space.norm() < (drive * duration).abs() * 3e-12,
                );
            }
            assert_eq!(
                result.bodies[0].coupling_angular_impulse_space
                    + result.bodies[1].coupling_angular_impulse_space,
                Vec3::default(),
            );
            assert!(
                result.angular_momentum_residual_space.norm() < (drive * duration).abs() * 3e-12,
            );
        }
    }
}

#[test]
fn stiff_coupled_eigenaxis_preserves_common_spin_under_world_and_material_rotation() {
    // A rank-one tensor has a non-coordinate active eigenvector and two exact
    // undamped directions. Independent material/world frames exercise both
    // transforms used to recover each body's momentum from the conserved sum.
    let moments = [2.0_f64, 5.0_f64];
    let duration = 0.1;
    let relaxation = 20.0_f64;
    let mobility = moments[0].recip() + moments[1].recip();
    let coefficient = relaxation / (2.0 * mobility * duration);
    let damping = Mat3 {
        m: [
            coefficient,
            coefficient,
            0.0,
            coefficient,
            coefficient,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    };
    let axis = Vec3::new(1.0, 1.0, 0.0);
    let second_frame = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let second_axis = second_frame.inverse().rotate(axis).unwrap();
    let world = So3::exp(So3Tangent::new(Vec3::new(-0.4, 0.3, 0.2))).unwrap();
    let amplitudes = [6.0_f64, -6.0 + 1e-5];
    let mean = (amplitudes[0] + amplitudes[1]) / (moments[0] + moments[1]);
    let slip = (amplitudes[0] / moments[0] - amplitudes[1] / moments[1]) * (-relaxation).exp();
    let expected = [
        axis.scale(moments[0] * (mean + moments[1] / (moments[0] + moments[1]) * slip)),
        second_axis.scale(moments[1] * (mean - moments[0] / (moments[0] + moments[1]) * slip)),
    ];
    let result = ViscousPairStepper::new_exponential(moments.map(sphere), damping, duration)
        .unwrap()
        .step(
            [world, world.compose(second_frame).unwrap()],
            [axis.scale(amplitudes[0]), second_axis.scale(amplitudes[1])],
            [HeldTorque::Body(Vec3::default()); 2],
        )
        .unwrap();
    for (body, target) in result.bodies.into_iter().zip(expected) {
        assert!((body.momentum_after_body - target).norm() < 5e-13);
    }
    balances(result, 4e-12);
}

#[test]
fn stiff_drag_preserves_exact_driven_cancellation_in_an_undamped_component() {
    // Only z is damped, while a body-held x drive removes body 0's entire
    // initial x momentum. Reconstructing that lossless zero through the
    // finite-precision modal dual basis can invent a residual transverse spin.
    let result = ViscousPairStepper::new_exponential(
        [sphere(2.0), sphere(5.0)],
        diagonal(0.0, 0.0, 1000.0),
        0.125,
    )
    .unwrap()
    .step(
        [So3::identity(); 2],
        [Vec3::new(1.0, 0.0, 6.0), Vec3::new(2.0, 0.0, -6.0)],
        [
            HeldTorque::Body(Vec3::new(-8.0, 0.0, 0.0)),
            HeldTorque::Body(Vec3::default()),
        ],
    )
    .unwrap();
    assert_eq!(result.bodies[0].momentum_after_body.x, 0.0);
    assert_eq!(result.bodies[0].momentum_after_body.y, 0.0);
    close(result.bodies[0].drive_work_j, -0.25, 2e-13);
    assert!(result.dissipated_energy_j > 0.0);
    balances(result, 4e-12);
}
