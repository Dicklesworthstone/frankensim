use fs_ga::{Mat3, So3, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, DryClutchLaw, DryClutchMode, DryClutchPort, RotorBank, RotorState,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn brake(static_capacity_nm: f64, kinetic_capacity_nm: f64) -> DryClutchLaw {
    DryClutchLaw {
        port: DryClutchPort::Brake { rotor: 0 },
        static_capacity_nm,
        kinetic_capacity_nm,
    }
}
fn bank() -> RotorBank<1> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 0.5).unwrap()],
    )
    .unwrap()
}
fn state(bank: &RotorBank<1>, carrier: f64, relative: f64, phase: f64) -> RotorState<1> {
    bank.state_from_velocities(
        So3::identity(),
        Vec3::new(0.0, 0.0, carrier),
        [relative],
        [phase],
    )
    .unwrap()
}

#[test]
fn g1_brake_stops_at_finite_time_and_retains_angle_with_the_exact_heat() {
    let bank = bank();
    let initial = state(&bank, 0.2, 3.0, 0.4);
    let stepper = bank.dry_clutch_stepper(brake(1.2, 0.8), 4.0).unwrap();
    close(stepper.mobility(), 2.5, 0.0);
    let step = stepper
        .step(&stepper.start(initial).unwrap(), [0.0], zero())
        .unwrap();
    assert_eq!(step.mode_after, DryClutchMode::Sticking);
    close(
        step.half_stages[0].synchronization_time_s.unwrap(),
        1.5,
        2e-15,
    );
    assert_eq!(step.half_stages[1].synchronization_time_s, None);
    close(step.half_stages[0].sticking_duration_s, 0.5, 2e-15);
    close(step.half_stages[1].sticking_duration_s, 2.0, 0.0);
    close(step.motion.state_after.relative_angle_rad[0], 2.65, 3e-15);
    close(
        bank.relative_rates(&step.motion.state_after).unwrap()[0],
        0.0,
        1e-15,
    );
    close(
        bank.angular_velocity(&step.motion.state_after).unwrap().z,
        0.8,
        1e-15,
    );
    close(step.friction_impulse_nms, 1.2, 2e-15);
    close(step.sliding_angle_rad, 2.25, 2e-15);
    close(step.dissipated_energy_j, 1.8, 3e-15);
    close(step.friction_work_j[0], -1.8, 3e-15);
    close(step.motion.kinetic_energy_change_j, -1.8, 3e-15);
    close(step.work_energy_residual_j, 0.0, 2e-15);
    assert_eq!(step.drive_work_j, [0.0]);
    assert_eq!(
        step.motion.state_after.total_momentum_body,
        initial.total_momentum_body
    );
    assert!(step.motion.angular_momentum_residual_world.norm() < 2e-15);

    // Carry the proved branch, without repeated inelastic engagement or drift
    // in the stopped phase. A fresh nonzero speed never gets this admission.
    let held = stepper.step(&step.state_after(), [0.0], zero()).unwrap();
    assert_eq!(held.mode_after, DryClutchMode::Sticking);
    assert_eq!(
        held.motion.state_after.relative_angle_rad,
        step.motion.state_after.relative_angle_rad
    );
    assert_eq!(held.dissipated_energy_j, 0.0);
    assert_eq!(held.friction_impulse_nms, 0.0);
    let small = stepper.start(state(&bank, 0.5, 1e-14, 0.0)).unwrap();
    assert_eq!(small.mode(), DryClutchMode::Slipping);
}

#[test]
fn g1_static_capacity_holds_breaks_away_and_carries_an_external_acceleration() {
    let bank = bank();
    let stepper = bank.dry_clutch_stepper(brake(1.2, 0.8), 0.4).unwrap();
    let initial = state(&bank, 0.0, 0.0, 1.25);
    let rest = stepper.start(initial).unwrap();
    for drive in [0.6, 1.2, -0.6, -1.2] {
        let step = stepper.step(&rest, [drive], zero()).unwrap();
        assert_eq!(step.mode_after, DryClutchMode::Sticking);
        assert_eq!(step.motion.state_after, initial);
        assert_eq!(step.dissipated_energy_j, 0.0);
        assert_eq!(step.drive_work_j, [0.0]);
        close(step.friction_impulse_nms, 0.4 * drive, 1e-16);
    }
    for sign in [1.0, -1.0] {
        let step = stepper.step(&rest, [sign * 1.5], zero()).unwrap();
        assert_eq!(step.mode_after, DryClutchMode::Slipping);
        close(
            bank.relative_rates(&step.motion.state_after).unwrap()[0],
            sign * 0.7,
            1e-15,
        );
        close(
            step.motion.state_after.relative_angle_rad[0],
            1.25 + sign * 0.14,
            1e-15,
        );
        close(
            bank.angular_velocity(&step.motion.state_after).unwrap().z,
            -sign * 0.14,
            5e-16,
        );
        close(step.drive_work_j[0], 0.21, 1e-15);
        close(step.dissipated_energy_j, 0.112, 1e-15);
        close(step.motion.kinetic_energy_change_j, 0.098, 1e-15);
        close(step.work_energy_residual_j, 0.0, 1e-15);
    }

    // Static does not mean inertially stationary: the shaft and carrier share
    // alpha = external/(I+J), with reaction drive-J*alpha and zero slip work.
    let moving = state(&bank, 0.3, 0.0, 0.2);
    let loaded = bank.dry_clutch_stepper(brake(1.0, 0.5), 0.8).unwrap();
    let step = loaded
        .step(
            &loaded.start(moving).unwrap(),
            [0.5],
            HeldTorque::Body(Vec3::new(0.0, 0.0, 1.0)),
        )
        .unwrap();
    assert_eq!(step.mode_after, DryClutchMode::Sticking);
    assert_eq!(
        step.motion.state_after.relative_angle_rad,
        moving.relative_angle_rad
    );
    close(
        bank.angular_velocity(&step.motion.state_after).unwrap().z,
        0.62,
        5e-16,
    );
    close(step.friction_impulse_nms, 0.24, 1e-15);
    close(step.motion.external_work_j, 0.368, 1e-15);
    close(step.motion.kinetic_energy_change_j, 0.368, 1e-15);
    assert_eq!(step.dissipated_energy_j, 0.0);
    assert_eq!(step.drive_work_j, [0.0]);
    close(step.work_energy_residual_j, 0.0, 1e-15);
}

#[test]
fn g1_overloaded_brake_reverses_after_synchronization_without_negative_heat() {
    let bank = bank();
    let initial = state(&bank, 0.1, 1.0, 0.2);
    let stepper = bank.dry_clutch_stepper(brake(1.0, 0.5), 1.0).unwrap();
    let step = stepper
        .step(&stepper.start(initial).unwrap(), [-1.5], zero())
        .unwrap();
    // v'=-5 until t=.2, then v'=-2.5 after the static limit is exceeded.
    close(
        step.half_stages[0].synchronization_time_s.unwrap(),
        0.2,
        2e-15,
    );
    assert_eq!(step.half_stages[0].sticking_duration_s, 0.0);
    assert_eq!(step.mode_after, DryClutchMode::Slipping);
    close(
        bank.relative_rates(&step.motion.state_after).unwrap()[0],
        -2.0,
        3e-15,
    );
    close(step.motion.state_after.relative_angle_rad[0], -0.5, 3e-15);
    close(step.sliding_angle_rad, 0.9, 2e-15);
    close(step.dissipated_energy_j, 0.45, 2e-15);
    close(step.friction_work_j[0], -0.45, 2e-15);
    close(step.friction_impulse_nms, -0.3, 1e-15);
    close(step.drive_work_j[0], 1.05, 2e-15);
    close(step.work_energy_residual_j, 0.0, 2e-15);
}

#[test]
fn g1_signed_ratio_clutch_synchronizes_both_shafts_with_full_carrier_recoil() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(4.0, 4.0, 3.0), So3::identity()).unwrap(),
        [0.4, 0.8].map(|j| AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), j).unwrap()),
    )
    .unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.0, 0.0, 0.2),
            [2.0, 0.5],
            [0.1, -0.2],
        )
        .unwrap();
    let law = DryClutchLaw {
        port: DryClutchPort::Shafts {
            first_rotor: 0,
            second_rotor: 1,
            ratio: -2.0,
        },
        static_capacity_nm: 0.75,
        kinetic_capacity_nm: 0.5,
    };
    let stepper = bank.dry_clutch_stepper(law, 2.0).unwrap();
    close(stepper.mobility(), 10.5, 1e-14);
    let step = stepper
        .step(&stepper.start(initial).unwrap(), [0.0; 2], zero())
        .unwrap();
    assert_eq!(step.mode_after, DryClutchMode::Sticking);
    close(
        step.half_stages[0].synchronization_time_s.unwrap(),
        4.0 / 7.0,
        1e-15,
    );
    close(step.friction_impulse_nms, 2.0 / 7.0, 1e-15);
    close(step.motion.axial_momentum_change_nms[0], -2.0 / 7.0, 1e-15);
    close(step.motion.axial_momentum_change_nms[1], -4.0 / 7.0, 1e-15);
    close(
        bank.angular_velocity(&step.motion.state_after).unwrap().z,
        0.2 + 2.0 / 7.0,
        1e-15,
    );
    let rates = bank.relative_rates(&step.motion.state_after).unwrap();
    close(rates[0], 1.0, 2e-15);
    close(rates[1], -0.5, 2e-15);
    close(
        step.motion.state_after.relative_angle_rad[0],
        0.1 + 16.0 / 7.0,
        3e-15,
    );
    close(
        step.motion.state_after.relative_angle_rad[1],
        -0.2 - 5.0 / 7.0,
        3e-15,
    );
    close(step.dissipated_energy_j, 3.0 / 7.0, 2e-15);
    close(step.sliding_angle_rad, 6.0 / 7.0, 2e-15);
    close(step.friction_work_j.iter().sum(), -3.0 / 7.0, 2e-15);
    close(step.work_energy_residual_j, 0.0, 3e-15);
    assert_eq!(
        step.motion.state_after.total_momentum_body,
        initial.total_momentum_body
    );
}

#[test]
fn g3_oblique_sliding_refines_against_original_euler_equations_and_drift_releases_stick() {
    let matrix = Mat3 {
        m: [2.4, 0.3, -0.2, 0.3, 3.1, 0.25, -0.2, 0.25, 4.0],
    };
    let inertia = InertiaTensor3::new(matrix).unwrap();
    let axes = [Vec3::new(1.0, 0.2, 0.0), Vec3::new(0.3, 1.0, 0.4)];
    let polar = [0.4, 0.7];
    let rotors = core::array::from_fn(|i| AxialRotor::new(axes[i], polar[i]).unwrap());
    let bank = RotorBank::from_reduced(inertia, rotors).unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.4, -0.3, 0.7),
            [2.5, -0.8],
            [0.2, -0.1],
        )
        .unwrap();
    let law = DryClutchLaw {
        port: DryClutchPort::Shafts {
            first_rotor: 0,
            second_rotor: 1,
            ratio: 0.6,
        },
        static_capacity_nm: 0.45,
        kinetic_capacity_nm: 0.3,
    };
    let drive = [0.2, -0.1];
    let torque = Vec3::new(0.07, -0.04, 0.03);
    let factor = lu(&matrix.m, 3).unwrap();
    let rhs = |y: [f64; 7]| {
        let l = Vec3::new(y[0], y[1], y[2]);
        let carrier = l - rotors[0].axis_body().scale(y[3]) - rotors[1].axis_body().scale(y[4]);
        let mut w = [carrier.x, carrier.y, carrier.z];
        factor.solve(&mut w);
        let omega = Vec3::new(w[0], w[1], w[2]);
        let s0 = y[3] / polar[0] - rotors[0].axis_body().dot(omega);
        let s1 = y[4] / polar[1] - rotors[1].axis_body().dot(omega);
        assert!(s0 - 0.6 * s1 > 1.0); // fixed positive-slip branch, independently checked
        let dl = l.cross(omega) + torque;
        [
            dl.x,
            dl.y,
            dl.z,
            drive[0] - 0.3,
            drive[1] + 0.6 * 0.3,
            s0,
            s1,
        ]
    };
    let packed = |s: &RotorState<2>| {
        [
            s.total_momentum_body.x,
            s.total_momentum_body.y,
            s.total_momentum_body.z,
            s.axial_momentum_nms[0],
            s.axial_momentum_nms[1],
            s.relative_angle_rad[0],
            s.relative_angle_rad[1],
        ]
    };
    let time = 0.6;
    let mut reference = packed(&initial);
    let dt = time / 4096.0;
    for _ in 0..4096 {
        let k1 = rhs(reference);
        let k2 = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * k1[i]));
        let k3 = rhs(core::array::from_fn(|i| reference[i] + 0.5 * dt * k2[i]));
        let k4 = rhs(core::array::from_fn(|i| reference[i] + dt * k3[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    let mut previous = None;
    for count in [8, 16, 32] {
        let stepper = bank.dry_clutch_stepper(law, time / count as f64).unwrap();
        let mut state = stepper.start(initial).unwrap();
        let mut accounted = 0.0;
        for _ in 0..count {
            let step = stepper
                .step(&state, drive, HeldTorque::Body(torque))
                .unwrap();
            assert_eq!(step.mode_after, DryClutchMode::Slipping);
            assert!(step.dissipated_energy_j > 0.0);
            assert!(step.motion.angular_momentum_residual_world.norm() < 2e-13);
            close(step.work_energy_residual_j, 0.0, 2e-13);
            accounted += step.drive_work_j.iter().sum::<f64>()
                + step.motion.external_work_j
                + step.motion.drift_energy_defect_j
                - step.dissipated_energy_j;
            state = step.state_after();
        }
        close(
            bank.kinetic_energy(state.rotors()).unwrap() - bank.kinetic_energy(&initial).unwrap(),
            accounted,
            2e-12,
        );
        let error = packed(state.rotors())
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>();
        if let Some(old) = previous {
            assert!(
                old / error > 3.8 && old / error < 4.2,
                "refinement ratio {}",
                old / error
            );
        }
        previous = Some(error);
    }

    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 0.5).unwrap()],
    )
    .unwrap();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.5, 0.4, 0.3), [0.0], [0.7])
        .unwrap();
    let stepper = bank.dry_clutch_stepper(brake(0.002, 0.001), 0.1).unwrap();
    let state = stepper.start(initial).unwrap();
    assert_eq!(state.mode(), DryClutchMode::Sticking);
    let step = stepper.step(&state, [0.0], zero()).unwrap();
    assert_eq!(step.half_stages[0].mode_after, DryClutchMode::Sticking);
    assert_eq!(step.half_stages[1].mode_after, DryClutchMode::Slipping);
    assert!(step.dissipated_energy_j > 0.0);
    assert_ne!(
        bank.relative_rates(&step.motion.state_after).unwrap()[0],
        0.0
    );
}

#[test]
fn g0_release_delegates_exactly_and_bad_inputs_refuse_atomically() {
    let bank = bank();
    let initial = state(&bank, 0.3, 0.7, 0.2);
    let open = bank.dry_clutch_stepper(brake(0.0, 0.0), 0.2).unwrap();
    let torque = HeldTorque::Space(Vec3::new(0.1, -0.2, 0.3));
    let free = bank.step(&initial, [0.4], torque, 0.2).unwrap();
    let result = open
        .step(&open.start(initial).unwrap(), [0.4], torque)
        .unwrap();
    assert_eq!(result.motion, free);
    assert_eq!(result.mode_after, DryClutchMode::Open);
    assert_eq!(result.dissipated_energy_j, 0.0);
    assert_eq!(result.friction_work_j, [0.0]);
    for duration in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(bank.dry_clutch_stepper(brake(1.0, 0.5), duration).is_err());
    }
    for (s, k) in [
        (1.0, 1.1),
        (-1.0, 0.0),
        (1.0, -0.1),
        (f64::NAN, 0.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(bank.dry_clutch_stepper(brake(s, k), 0.1).is_err());
    }
    for port in [
        DryClutchPort::Brake { rotor: 1 },
        DryClutchPort::Shafts {
            first_rotor: 0,
            second_rotor: 0,
            ratio: 1.0,
        },
        DryClutchPort::Shafts {
            first_rotor: 0,
            second_rotor: 1,
            ratio: f64::NAN,
        },
    ] {
        assert!(
            bank.dry_clutch_stepper(
                DryClutchLaw {
                    port,
                    ..brake(1.0, 0.5)
                },
                0.1
            )
            .is_err()
        );
    }
    let stepper = bank.dry_clutch_stepper(brake(1.0, 0.5), 0.1).unwrap();
    let retained = stepper.start(initial).unwrap();
    let before = retained;
    assert!(stepper.step(&retained, [f64::NAN], zero()).is_err());
    assert!(
        stepper
            .step(
                &retained,
                [0.0],
                HeldTorque::Body(Vec3::new(f64::INFINITY, 0.0, 0.0))
            )
            .is_err()
    );
    assert_eq!(retained, before);
    // A turnaround's cancelled net phase is admitted at its endpoint-motion
    // scale; this does not excuse erasing a genuine monotone phase increment.
    let coarse_phase = stepper.start(state(&bank, 0.3, 0.7, 1e16)).unwrap();
    assert!(matches!(
        stepper.step(&coarse_phase, [0.0], zero()),
        Err(fs_time::lie::LieStepError::InvalidParameter {
            context: "dry clutch relative phase"
        })
    ));

    let empty = RotorBank::from_reduced(bank.reduced_inertia(), []).unwrap();
    assert!(empty.dry_clutch_stepper(brake(1.0, 0.5), 0.1).is_err());
    let full = RotorBank::from_reduced(bank.reduced_inertia(), [bank.rotors()[0]; 64]).unwrap();
    let port = DryClutchPort::Brake { rotor: 63 };
    let prepared = full
        .dry_clutch_stepper(
            DryClutchLaw {
                port,
                ..brake(1.0, 0.5)
            },
            0.1,
        )
        .unwrap();
    let rest = full
        .state_from_velocities(So3::identity(), Vec3::default(), [0.0; 64], [0.0; 64])
        .unwrap();
    let step = prepared
        .step(&prepared.start(rest).unwrap(), [0.0; 64], zero())
        .unwrap();
    assert_eq!(step.motion.state_after, rest);
}
