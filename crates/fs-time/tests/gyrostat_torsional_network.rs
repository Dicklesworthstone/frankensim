use fs_ga::{Mat3, So3, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, MAX_TORSIONAL_COUPLINGS, RotorBank, RotorState, TorsionalCoupling,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn coupling(ratio: f64, stiffness: f64, damping: f64) -> TorsionalCoupling {
    TorsionalCoupling {
        first_rotor: 0,
        second_rotor: 1,
        ratio,
        stiffness_nm_per_rad: stiffness,
        damping_nms_per_rad: damping,
        rest_angle_rad: 0.1,
    }
}
fn coaxial() -> RotorBank<2> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [0.4, 0.6].map(|j| AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), j).unwrap()),
    )
    .unwrap()
}

#[test]
fn g1_two_shaft_normal_mode_has_the_physical_frequency_energy_and_no_carrier_recoil() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.0, 0.0, 0.3),
            [0.6, -0.4],
            [0.3, -0.2],
        )
        .unwrap();
    let k = 3.0_f64;
    let frequency = (k * (1.0 / 0.4 + 1.0 / 0.6)).sqrt();
    let period = 2.0 * core::f64::consts::PI / frequency;
    let mut previous = None;
    for count in [24, 48, 96] {
        let stepper = bank
            .torsional_network_stepper([coupling(1.0, k, 0.0)], period / count as f64)
            .unwrap();
        close(
            stepper.angular_frequency_bound_rad_s() / frequency,
            1.0,
            1e-12,
        );
        let initial_energy =
            bank.kinetic_energy(&initial).unwrap() + stepper.potential_energy(&initial).unwrap();
        let mut state = initial;
        for _ in 0..count {
            let step = stepper.step(&state, [0.0; 2], zero()).unwrap();
            assert_eq!(step.dissipated_energy_j, [0.0]);
            assert_eq!(step.drive_work_j, [0.0; 2]);
            close(step.work_energy_residual_j, 0.0, 2e-13);
            state = step.motion.state_after;
            close(bank.angular_velocity(&state).unwrap().z, 0.3, 2e-14);
            close(
                state.axial_momentum_nms.iter().sum(),
                initial.axial_momentum_nms.iter().sum(),
                2e-14,
            );
            close(
                bank.kinetic_energy(&state).unwrap() + stepper.potential_energy(&state).unwrap(),
                initial_energy,
                5e-13,
            );
        }
        let rates = bank.relative_rates(&state).unwrap();
        let error = (state.relative_angle_rad[0] - state.relative_angle_rad[1] - 0.5).abs()
            + (rates[0] - rates[1] - 1.0).abs();
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2);
        }
        previous = Some(error);
    }
}

#[test]
fn g1_signed_ratio_loaded_damped_motion_closes_drive_heat_and_external_work() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.0, 0.0, 0.3),
            [0.5, -0.2],
            [0.2, -0.1],
        )
        .unwrap();
    let mut c = coupling(-2.0, 1.2, 0.1);
    c.rest_angle_rad = 0.05;
    let drive = [0.4, -0.2];
    let external = 0.3;
    let time = 0.6;
    // The exact scalar strain equation retains carrier recoil and both drives:
    // x'' + H*d*x' + H*k*x = B*M*drive - (1-ratio)*external/I.
    let h = 1.0 / 0.4 + c.ratio * c.ratio / 0.6 + (1.0 - c.ratio).powi(2) / 2.0;
    let forcing = drive[0] / 0.4 - c.ratio * drive[1] / 0.6
        + (1.0 - c.ratio) * (drive[0] + drive[1] - external) / 2.0;
    let equilibrium = forcing / (h * c.stiffness_nm_per_rad);
    let decay = 0.5 * h * c.damping_nms_per_rad;
    let wd = (h * c.stiffness_nm_per_rad - decay * decay).sqrt();
    let a = initial.relative_angle_rad[0]
        - c.ratio * initial.relative_angle_rad[1]
        - c.rest_angle_rad
        - equilibrium;
    let b = (0.5 - c.ratio * -0.2 + decay * a) / wd;
    let expected_x =
        equilibrium + (-decay * time).exp() * (a * (wd * time).cos() + b * (wd * time).sin());
    let expected_rate = (-decay * time).exp()
        * ((-decay * a + wd * b) * (wd * time).cos() + (-decay * b - wd * a) * (wd * time).sin());
    let mut previous = None;
    for count in [16, 32, 64] {
        let stepper = bank
            .torsional_network_stepper([c], time / count as f64)
            .unwrap();
        let initial_energy =
            bank.kinetic_energy(&initial).unwrap() + stepper.potential_energy(&initial).unwrap();
        let mut state = initial;
        let (mut work, mut heat, mut external_work, mut drift) = ([0.0; 2], 0.0, 0.0, 0.0);
        for _ in 0..count {
            let step = stepper
                .step(
                    &state,
                    drive,
                    HeldTorque::Body(Vec3::new(0.0, 0.0, external)),
                )
                .unwrap();
            assert!(step.dissipated_energy_j[0] >= 0.0);
            close(step.work_energy_residual_j, 0.0, 2e-13);
            for i in 0..2 {
                work[i] += step.drive_work_j[i];
            }
            heat += step.dissipated_energy_j[0];
            external_work += step.motion.external_work_j;
            drift += step.motion.drift_energy_defect_j;
            state = step.motion.state_after;
        }
        assert!(heat > 0.0);
        for i in 0..2 {
            close(
                work[i],
                drive[i] * (state.relative_angle_rad[i] - initial.relative_angle_rad[i]),
                5e-13,
            );
        }
        close(
            state.total_momentum_body.z - initial.total_momentum_body.z,
            external * time,
            5e-13,
        );
        close(
            bank.kinetic_energy(&state).unwrap() + stepper.potential_energy(&state).unwrap()
                - initial_energy,
            work.iter().sum::<f64>() + external_work - heat + drift,
            5e-13,
        );
        let rates = bank.relative_rates(&state).unwrap();
        let error = (state.relative_angle_rad[0]
            - c.ratio * state.relative_angle_rad[1]
            - c.rest_angle_rad
            - expected_x)
            .abs()
            + (rates[0] - c.ratio * rates[1] - expected_rate).abs();
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2);
        }
        previous = Some(error);
    }
}

#[test]
fn g1_coupled_oblique_network_refines_against_original_euler_equations() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::new(Mat3 {
            m: [2.4, 0.3, -0.2, 0.3, 3.1, 0.25, -0.2, 0.25, 4.0],
        })
        .unwrap(),
        [
            AxialRotor::new(Vec3::new(1.0, 0.2, 0.0), 0.4).unwrap(),
            AxialRotor::new(Vec3::new(0.3, 1.0, 0.4), 0.7).unwrap(),
            AxialRotor::new(Vec3::new(-0.2, 0.1, 1.0), 0.9).unwrap(),
        ],
    )
    .unwrap();
    let couplings = [
        coupling(2.0, 0.7, 0.1),
        TorsionalCoupling {
            first_rotor: 1,
            second_rotor: 2,
            ratio: 0.5,
            stiffness_nm_per_rad: 1.1,
            damping_nms_per_rad: 0.2,
            rest_angle_rad: 0.05,
        },
        TorsionalCoupling {
            first_rotor: 0,
            second_rotor: 2,
            ratio: 1.0,
            stiffness_nm_per_rad: 0.4,
            damping_nms_per_rad: 0.05,
            rest_angle_rad: -0.08,
        },
    ];
    // B's third row is the first plus twice the second: dependent loops are
    // still a valid elastic network, without replacing it by rigid constraints.
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.4, -0.3, 0.2),
            [0.6, -0.4, 0.2],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let drive = [0.15, -0.12, 0.07];
    let external = Vec3::new(0.1, 0.05, -0.02);
    let factor = lu(&bank.reduced_inertia().matrix().m, 3).unwrap();
    let pack = |state: &RotorState<3>| {
        [
            state.total_momentum_body.x,
            state.total_momentum_body.y,
            state.total_momentum_body.z,
            state.axial_momentum_nms[0],
            state.axial_momentum_nms[1],
            state.axial_momentum_nms[2],
            state.relative_angle_rad[0],
            state.relative_angle_rad[1],
            state.relative_angle_rad[2],
        ]
    };
    let derivative = |y: [f64; 9]| {
        let total = Vec3::new(y[0], y[1], y[2]);
        let mut carrier = total;
        for i in 0..3 {
            carrier = carrier - bank.rotors()[i].axis_body().scale(y[3 + i]);
        }
        let mut solve = [carrier.x, carrier.y, carrier.z];
        factor.solve(&mut solve);
        let omega = Vec3::new(solve[0], solve[1], solve[2]);
        let torque = total.cross(omega) + external;
        let mut out = [
            torque.x, torque.y, torque.z, drive[0], drive[1], drive[2], 0.0, 0.0, 0.0,
        ];
        for i in 0..3 {
            out[6 + i] = y[3 + i] / bank.rotors()[i].axial_inertia_kg_m2()
                - bank.rotors()[i].axis_body().dot(omega);
        }
        for c in couplings {
            let strain = y[6 + c.first_rotor] - c.ratio * y[6 + c.second_rotor] - c.rest_angle_rad;
            let rate = out[6 + c.first_rotor] - c.ratio * out[6 + c.second_rotor];
            let effort = c.stiffness_nm_per_rad * strain + c.damping_nms_per_rad * rate;
            out[3 + c.first_rotor] -= effort;
            out[3 + c.second_rotor] += c.ratio * effort;
        }
        out
    };
    let time = 0.4;
    let dt = time / 4096.0;
    let mut reference = pack(&initial);
    for _ in 0..4096 {
        let k1 = derivative(reference);
        let k2 = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| reference[i] + dt * k3[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        });
    }
    let mut previous = None;
    for count in [20, 40, 80] {
        let stepper = bank
            .torsional_network_stepper(couplings, time / count as f64)
            .unwrap();
        let mut state = initial;
        for _ in 0..count {
            let step = stepper
                .step(&state, drive, HeldTorque::Body(external))
                .unwrap();
            assert!(step.dissipated_energy_j.iter().all(|heat| *heat >= 0.0));
            close(step.work_energy_residual_j, 0.0, 2e-13);
            assert!(step.motion.angular_momentum_residual_world.norm() < 2e-13);
            state = step.motion.state_after;
        }
        let error = pack(&state)
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        if let Some(old) = previous {
            assert!(old / error > 3.7 && old / error < 4.3, "{old:e}/{error:e}");
        }
        previous = Some(error);
    }
}

#[test]
fn g0_zero_anchor_parallel_edge_limits_and_bounded_atomic_refusals() {
    let bank = coaxial();
    let state = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.0, 0.0, 0.2),
            [0.5, -0.1],
            [0.3, -0.2],
        )
        .unwrap();
    let drive = [0.1, -0.2];
    let dt = 0.01;
    let held = bank.step(&state, drive, zero(), dt).unwrap();
    assert_eq!(
        bank.torsional_network_stepper([], dt)
            .unwrap()
            .step(&state, drive, zero())
            .unwrap()
            .motion,
        held
    );
    assert_eq!(
        bank.torsional_network_stepper([coupling(1.0, 0.0, 0.0)], dt)
            .unwrap()
            .step(&state, drive, zero())
            .unwrap()
            .motion,
        held
    );
    let anchor = bank
        .torsional_network_stepper([coupling(0.0, 2.0, 0.3)], dt)
        .unwrap()
        .step(&state, drive, zero())
        .unwrap();
    let existing = bank
        .torsional_stepper([2.0, 0.0], [0.1, 0.0], [0.3, 0.0], dt)
        .unwrap()
        .step(&state, drive, zero())
        .unwrap();
    for i in 0..2 {
        close(
            anchor.motion.state_after.axial_momentum_nms[i],
            existing.motion.state_after.axial_momentum_nms[i],
            2e-14,
        );
        close(
            anchor.motion.state_after.relative_angle_rad[i],
            existing.motion.state_after.relative_angle_rad[i],
            2e-14,
        );
    }
    let pair = bank
        .torsional_network_stepper([coupling(1.0, 2.0, 0.2), coupling(1.0, 1.0, 0.1)], dt)
        .unwrap()
        .step(&state, drive, zero())
        .unwrap();
    let combined = bank
        .torsional_network_stepper([coupling(1.0, 3.0, 0.3)], dt)
        .unwrap()
        .step(&state, drive, zero())
        .unwrap();
    for i in 0..2 {
        close(
            pair.motion.state_after.axial_momentum_nms[i],
            combined.motion.state_after.axial_momentum_nms[i],
            2e-14,
        );
    }
    assert!(
        bank.torsional_network_stepper([coupling(1.0, 0.001, 0.001); MAX_TORSIONAL_COUPLINGS], dt)
            .is_ok()
    );
    assert!(
        bank.torsional_network_stepper(
            [coupling(1.0, 0.001, 0.001); MAX_TORSIONAL_COUPLINGS + 1],
            dt
        )
        .is_err()
    );
    for bad in [
        TorsionalCoupling {
            ratio: f64::NAN,
            ..coupling(1.0, 1.0, 0.0)
        },
        TorsionalCoupling {
            second_rotor: 0,
            ..coupling(1.0, 1.0, 0.0)
        },
        coupling(1.0, -1.0, 0.0),
        coupling(1.0, 1.0, -1.0),
    ] {
        assert!(bank.torsional_network_stepper([bad], dt).is_err());
    }
    assert!(
        bank.torsional_network_stepper([coupling(1.0, 1e6, 1e6)], 1.0)
            .is_err()
    );
    let stepper = bank
        .torsional_network_stepper([coupling(1.0, 1.0, 0.2)], dt)
        .unwrap();
    let original = state;
    assert!(stepper.step(&state, [f64::NAN, 0.0], zero()).is_err());
    assert_eq!(state, original);
}
