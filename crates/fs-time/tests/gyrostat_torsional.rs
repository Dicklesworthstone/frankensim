use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{AxialRotor, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn zero_torque() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
fn coaxial() -> RotorBank<1> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap()],
    )
    .unwrap()
}
fn oblique() -> RotorBank<3> {
    RotorBank::from_reduced(
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
    .unwrap()
}

#[test]
fn free_torsional_period_energy_phase_and_carrier_recoil_match_the_harmonic_solution() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.3, 0.1)),
            Vec3::new(0.0, 0.0, 0.2),
            [0.7],
            [0.3],
        )
        .unwrap();
    let stiffness = 4.0_f64;
    let frequency = (1.5 * stiffness).sqrt(); // Full carrier recoil: 1/J + 1/I_z.
    let period = 2.0 * core::f64::consts::PI / frequency;
    let mut previous = None;
    for count in [24, 48, 96] {
        let stepper = bank
            .torsional_stepper([stiffness], [0.1], [0.0], period / count as f64)
            .unwrap();
        close(
            stepper.angular_frequency_bound_rad_s() / frequency,
            1.0,
            1e-12,
        );
        let energy0 =
            bank.kinetic_energy(&initial).unwrap() + stepper.potential_energy(&initial).unwrap();
        let mut state = initial;
        let mut maximum_energy_error = 0.0_f64;
        for i in 1..=count {
            let step = stepper.step(&state, [0.0], zero_torque()).unwrap();
            assert_eq!(step.dissipated_energy_j, [0.0]);
            assert_eq!(step.drive_work_j, [0.0]);
            close(step.work_energy_residual_j, 0.0, 3e-14);
            close(step.motion.drift_energy_defect_j, 0.0, 1e-15);
            state = step.motion.state_after;
            let energy =
                bank.kinetic_energy(&state).unwrap() + stepper.potential_energy(&state).unwrap();
            maximum_energy_error = maximum_energy_error.max((energy - energy0).abs());
            // A changing spring force must exchange real rotor/carrier motion.
            if i == count / 4 {
                assert!((bank.angular_velocity(&state).unwrap().z - 0.2).abs() > 0.1);
            }
        }
        let expected_rotation = initial
            .rotation
            .compose(rotation(Vec3::new(
                0.0,
                0.0,
                initial.total_momentum_body.z * period / 3.0,
            )))
            .unwrap();
        let error = (state.relative_angle_rad[0] - initial.relative_angle_rad[0]).abs()
            + (bank.relative_rates(&state).unwrap()[0] - 0.7).abs()
            + state
                .rotation
                .body_minus(expected_rotation)
                .unwrap()
                .angular
                .norm();
        eprintln!(
            "torsional harmonic period: steps={count}, error={error:e}, energy_drift={maximum_energy_error:e}"
        );
        assert!(maximum_energy_error < 5e-13);
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2, "{old:e}/{error:e}");
        }
        previous = Some(error);
        assert!(
            (state.rotation.rotate(state.total_momentum_body).unwrap()
                - initial
                    .rotation
                    .rotate(initial.total_momentum_body)
                    .unwrap())
            .norm()
                < 3e-14
        );
    }
}

#[test]
fn damped_driven_spring_refines_to_analytic_motion_and_distinguishes_heat_from_storage() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.8], [0.1])
        .unwrap();
    let time = 0.6_f64;
    let stiffness = 3.0_f64;
    let damping = 0.4_f64;
    let drive = 0.3_f64;
    let rest = 0.2_f64;
    let equilibrium = rest + drive / stiffness;
    let decay = 0.5 * 1.5 * damping;
    let wd = (1.5 * stiffness - decay * decay).sqrt();
    let a = initial.relative_angle_rad[0] - equilibrium;
    let b = (0.8 + decay * a) / wd;
    let expected_q =
        equilibrium + (-decay * time).exp() * (a * (wd * time).cos() + b * (wd * time).sin());
    let expected_s = (-decay * time).exp()
        * ((-decay * a + wd * b) * (wd * time).cos() + (-decay * b - wd * a) * (wd * time).sin());
    let mut previous = None;
    for count in [16, 32, 64] {
        let stepper = bank
            .torsional_stepper([stiffness], [rest], [damping], time / count as f64)
            .unwrap();
        let energy0 =
            bank.kinetic_energy(&initial).unwrap() + stepper.potential_energy(&initial).unwrap();
        let mut state = initial;
        let mut heat = 0.0;
        let mut work = 0.0;
        let mut drift = 0.0;
        for _ in 0..count {
            let step = stepper.step(&state, [drive], zero_torque()).unwrap();
            assert!(step.dissipated_energy_j[0] >= 0.0);
            close(step.work_energy_residual_j, 0.0, 3e-14);
            close(
                step.elastic_energy_change_j[0],
                stepper.potential_energy(&step.motion.state_after).unwrap()
                    - stepper.potential_energy(&state).unwrap(),
                1e-14,
            );
            heat += step.dissipated_energy_j[0];
            work += step.drive_work_j[0];
            drift += step.motion.drift_energy_defect_j;
            state = step.motion.state_after;
        }
        assert!(heat > 0.0 && work > 0.0);
        close(
            work,
            drive * (state.relative_angle_rad[0] - initial.relative_angle_rad[0]),
            2e-14,
        );
        close(
            bank.kinetic_energy(&state).unwrap() + stepper.potential_energy(&state).unwrap()
                - energy0,
            work - heat + drift,
            4e-14,
        );
        let error = (state.relative_angle_rad[0] - expected_q).abs()
            + (bank.relative_rates(&state).unwrap()[0] - expected_s).abs();
        eprintln!("damped driven torsion: steps={count}, error={error:e}");
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2, "{old:e}/{error:e}");
        }
        previous = Some(error);
    }
}

#[test]
fn oblique_full_tensor_springs_refine_against_independent_original_frame_euler_equations() {
    let bank = oblique();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.3, 0.1)),
            Vec3::new(0.6, -0.3, 0.8),
            [1.0, -0.5, 0.7],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let stiffness = [0.7, 1.1, 0.5];
    let damping = [0.2, 0.1, 0.3];
    let rest = [0.1, 0.2, -0.1];
    let drive = [0.1, -0.2, 0.4];
    let external = Vec3::new(0.1, 0.2, -0.1);
    let factor = lu(&bank.reduced_inertia().matrix().m, 3).unwrap();
    // Independent dense Euler/canonical ODE. No native principal basis,
    // gyrostat split, midpoint stage or factor elimination enters this oracle.
    let derivative = |y: [f64; 9]| {
        let total = Vec3::new(y[0], y[1], y[2]);
        let mut carrier = total;
        for i in 0..3 {
            carrier = carrier - bank.rotors()[i].axis_body().scale(y[3 + i]);
        }
        let mut omega = [carrier.x, carrier.y, carrier.z];
        factor.solve(&mut omega);
        let omega = Vec3::new(omega[0], omega[1], omega[2]);
        let rate = total.cross(omega) + external;
        let mut out = [rate.x, rate.y, rate.z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for i in 0..3 {
            out[6 + i] = y[3 + i] / bank.rotors()[i].axial_inertia_kg_m2()
                - bank.rotors()[i].axis_body().dot(omega);
            out[3 + i] = drive[i] - stiffness[i] * (y[6 + i] - rest[i]) - damping[i] * out[6 + i];
        }
        out
    };
    let pack = |s: &RotorState<3>| {
        [
            s.total_momentum_body.x,
            s.total_momentum_body.y,
            s.total_momentum_body.z,
            s.axial_momentum_nms[0],
            s.axial_momentum_nms[1],
            s.axial_momentum_nms[2],
            s.relative_angle_rad[0],
            s.relative_angle_rad[1],
            s.relative_angle_rad[2],
        ]
    };
    let time = 0.7;
    let dt = time / 8192.0;
    let mut reference = pack(&initial);
    for _ in 0..8192 {
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
            .torsional_stepper(stiffness, rest, damping, time / count as f64)
            .unwrap();
        let mut state = initial;
        for _ in 0..count {
            let step = stepper
                .step(&state, drive, HeldTorque::Body(external))
                .unwrap();
            close(step.work_energy_residual_j, 0.0, 5e-14);
            assert!(step.dissipated_energy_j.iter().all(|value| *value >= 0.0));
            assert!(step.motion.angular_momentum_residual_world.norm() < 3e-14);
            state = step.motion.state_after;
        }
        let error = pack(&state)
            .into_iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        eprintln!("oblique torsion Euler reference: steps={count}, error={error:e}");
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2, "{old:e}/{error:e}");
        }
        previous = Some(error);
    }
}

#[test]
fn zero_stiffness_delegates_exactly_and_unresolved_phase_or_mechanics_refuses_atomically() {
    let bank = coaxial();
    let state = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.8], [0.1])
        .unwrap();
    let external = HeldTorque::Space(Vec3::new(0.1, -0.2, 0.3));
    for damping in [0.0, 0.4] {
        let expected = bank
            .viscous_stepper([damping], 0.05)
            .unwrap()
            .step(&state, [0.3], external)
            .unwrap();
        let actual = bank
            .torsional_stepper([0.0], [123.0], [damping], 0.05)
            .unwrap()
            .step(&state, [0.3], external)
            .unwrap();
        assert_eq!(actual.motion, expected.motion);
        assert_eq!(actual.drive_work_j, expected.drive_work_j);
        assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
        assert_eq!(
            actual.work_energy_residual_j,
            expected.work_energy_residual_j
        );
        assert_eq!(actual.elastic_energy_change_j, [0.0]);
    }
    let before = state;
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(bank.torsional_stepper([bad], [0.0], [0.0], 0.1).is_err());
        assert!(bank.torsional_stepper([1.0], [0.0], [bad], 0.1).is_err());
    }
    assert!(
        bank.torsional_stepper([1.0], [f64::NAN], [0.0], 0.1)
            .is_err()
    );
    assert!(bank.torsional_stepper([1e12], [0.0], [0.0], 1.0).is_err());
    let resolved = bank.torsional_stepper([1e12], [0.0], [0.0], 1e-6).unwrap();
    assert!(resolved.half_stage_resolution_number() < 2.0);
    for rate in [1.0, 5.0] {
        let coarse = bank
            .state_from_velocities(So3::identity(), Vec3::default(), [rate], [1e16])
            .unwrap();
        let saved = coarse;
        // At this stored phase, the first update either vanishes or is only
        // partially representable; a small local strain must not hide that loss.
        assert!(
            bank.torsional_stepper([1.0], [1e16], [0.0], 1.0)
                .unwrap()
                .step(&coarse, [0.0], zero_torque())
                .is_err()
        );
        assert_eq!(coarse, saved);
    }
    let ordinary = bank.torsional_stepper([1.0], [0.0], [0.1], 0.1).unwrap();
    assert!(ordinary.step(&state, [f64::NAN], zero_torque()).is_err());
    assert_eq!(state, before);
    assert!(ordinary.step(&state, [0.0], zero_torque()).is_ok());
}
