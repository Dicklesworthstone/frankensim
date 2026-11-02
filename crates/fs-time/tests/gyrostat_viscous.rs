use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{AxialRotor, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}
fn zero_torque() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tol {tolerance:e})"
    );
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
fn coaxial_drive_drag_phase_and_heat_refine_to_the_closed_exponential() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.3, 0.1)),
            Vec3::new(0.0, 0.0, 0.3),
            [2.0],
            [0.2],
        )
        .unwrap();
    let time: f64 = 0.8;
    let damping: f64 = 0.7;
    let drive: f64 = 0.4;
    let mobility: f64 = 1.5; // 1/J + 1/I_z, including carrier recoil.
    let lambda = damping * mobility;
    let steady = drive / damping;
    let decayed = (-lambda * time).exp();
    let relative = steady + (2.0 - steady) * decayed;
    let phase = steady * time + (2.0 - steady) * (1.0 - decayed) / lambda;
    let p = (relative + initial.total_momentum_body.z / 2.0) / mobility;
    let carrier_angle = initial.total_momentum_body.z * time / 3.0 - phase / 3.0;
    let expected_rotation = initial
        .rotation
        .compose(rotation(Vec3::new(0.0, 0.0, carrier_angle)))
        .unwrap();
    let exact_heat = damping
        * (steady * steady * time
            + 2.0 * steady * (2.0 - steady) * (1.0 - decayed) / lambda
            + (2.0 - steady).powi(2) * (1.0 - decayed * decayed) / (2.0 * lambda));
    let mut previous = None;
    for count in [8, 16, 32] {
        let stepper = bank
            .viscous_stepper([damping], time / count as f64)
            .unwrap();
        let mut state = initial;
        let mut heat = 0.0;
        let mut drive_work = 0.0;
        let mut drift = 0.0;
        for _ in 0..count {
            let step = stepper.step(&state, [drive], zero_torque()).unwrap();
            assert!(step.dissipated_energy_j[0] >= 0.0);
            close(step.work_energy_residual_j, 0.0, 3e-14);
            heat += step.dissipated_energy_j[0];
            drive_work += step.drive_work_j[0];
            drift += step.motion.drift_energy_defect_j;
            state = step.motion.state_after;
        }
        let error = (state.axial_momentum_nms[0] - p).abs()
            + (state.relative_angle_rad[0] - initial.relative_angle_rad[0] - phase).abs()
            + state
                .rotation
                .body_minus(expected_rotation)
                .unwrap()
                .angular
                .norm()
            + (heat - exact_heat).abs();
        eprintln!("coaxial exponential/phase/heat: steps={count}, error={error:e}");
        if let Some(old) = previous {
            assert!(
                old / error > 3.8 && old / error < 4.2,
                "{old:e} / {error:e}"
            );
        }
        previous = Some(error);
        close(
            drive_work,
            drive * (state.relative_angle_rad[0] - initial.relative_angle_rad[0]),
            2e-14,
        );
        close(
            bank.kinetic_energy(&state).unwrap() - bank.kinetic_energy(&initial).unwrap(),
            drive_work - heat + drift,
            3e-14,
        );
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
fn oblique_full_tensor_motion_refines_against_independent_dense_euler_dynamics() {
    let bank = oblique();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.3, 0.1)),
            Vec3::new(0.6, -0.3, 0.8),
            [1.0, -0.5, 0.7],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let damping = [0.3, 0.7, 0.2];
    let drive = [0.1, -0.2, 0.4];
    let external = Vec3::new(0.1, 0.2, -0.1);
    // Independent full body-coordinate Euler equations. This reference uses
    // original dense I and a native LU solve, no principal basis, gyrostat
    // drift, rank-three elimination or midpoint shaft owner.
    let factor = lu(&bank.reduced_inertia().matrix().m, 3).unwrap();
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
            out[3 + i] = drive[i] - damping[i] * out[6 + i];
        }
        out
    };
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
    let mut reference = pack(&initial);
    let dt = 1.0 / 8192.0;
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
        let stepper = bank.viscous_stepper(damping, 1.0 / count as f64).unwrap();
        let mut state = initial;
        for _ in 0..count {
            let step = stepper
                .step(&state, drive, HeldTorque::Body(external))
                .unwrap();
            assert!(step.dissipated_energy_j.iter().all(|heat| *heat >= 0.0));
            close(step.work_energy_residual_j, 0.0, 3e-14);
            assert!(step.motion.angular_momentum_residual_world.norm() < 2e-14);
            state = step.motion.state_after;
        }
        let error = pack(&state)
            .into_iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        eprintln!("oblique dense Euler reference: steps={count}, error={error:e}");
        if let Some(old) = previous {
            assert!(
                old / error > 3.8 && old / error < 4.2,
                "{old:e} / {error:e}"
            );
        }
        previous = Some(error);
    }
}

#[test]
fn stiff_drag_exposes_its_timescale_and_preserves_passivity_when_resolved() {
    let bank = coaxial();
    let damping = 1e10;
    let lambda = 1.5 * damping;
    close(
        bank.viscous_relaxation_rate_bound([damping]).unwrap() / lambda,
        1.0,
        1e-12,
    );
    assert!(bank.viscous_stepper([damping], 10.0 / lambda).is_err());
    let h = 0.4 / lambda;
    let stepper = bank.viscous_stepper([damping], h).unwrap();
    close(stepper.half_stage_relaxation_number(), 0.2, 1e-12);
    close(stepper.duration_s(), h, 0.0);
    assert_eq!(*stepper.damping_nms(), [damping]);
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [3.0], [0.0])
        .unwrap();
    let half_decay = (1.0 - 0.25 * h * lambda) / (1.0 + 0.25 * h * lambda);
    let mut state = initial;
    let mut energy = bank.kinetic_energy(&state).unwrap();
    for index in 1..=20 {
        let step = stepper.step(&state, [0.0], zero_torque()).unwrap();
        assert!(step.dissipated_energy_j[0] > 0.0);
        close(step.drive_work_j[0], 0.0, 0.0);
        close(step.work_energy_residual_j, 0.0, 3e-14);
        state = step.motion.state_after;
        close(
            bank.relative_rates(&state).unwrap()[0],
            3.0 * half_decay.powi(2 * index),
            2e-14,
        );
        let next_energy = bank.kinetic_energy(&state).unwrap();
        assert!(next_energy <= energy);
        close(energy - next_energy, step.dissipated_energy_j[0], 3e-14);
        energy = next_energy;
    }
    // A motor holding relative speed supplies exactly the viscous heat; zero
    // net shaft torque does not mean either power source vanishes.
    let equilibrium = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.5), [1.0], [0.0])
        .unwrap();
    let steady = bank
        .viscous_stepper([1.0], 0.5)
        .unwrap()
        .step(&equilibrium, [1.0], zero_torque())
        .unwrap();
    close(steady.drive_work_j[0], 0.5, 1e-14);
    close(steady.dissipated_energy_j[0], 0.5, 1e-14);
    close(steady.motion.kinetic_energy_change_j, 0.0, 1e-14);
    let braking = bank
        .viscous_stepper([0.5], 0.1)
        .unwrap()
        .step(&equilibrium, [-0.2], zero_torque())
        .unwrap();
    assert!(braking.drive_work_j[0] < 0.0 && braking.dissipated_energy_j[0] > 0.0);
}

#[test]
fn unresolved_stage_work_refuses_even_when_the_rounded_shaft_law_closes() {
    let bank = coaxial();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 1e16), [-1e16], [0.0])
        .unwrap();
    let before = initial;
    // The shaft receives exactly 0.5, but the simultaneous external transfer
    // leaves required carrier Delta L=3.5, rounded to 4 on its huge baseline.
    // The old native allowance admits that partial recoil. The drag equation
    // still closes exactly, so checking only the shaft law misses 5e15 J.
    let raw = bank
        .apply_impulse(&initial, [0.5], Vec3::new(0.0, 0.0, 4.0))
        .unwrap();
    close(raw.axial_momentum_change_nms[0], 0.5, 0.0);
    close(-0.5e-16 * raw.mean_relative_rate_rad_s[0], 0.5, 1e-15);
    assert!(raw.work_energy_residual_j.abs() > 1e15);
    let stepper = bank.viscous_stepper([1e-16], 1.0).unwrap();
    assert!(
        stepper
            .step(&initial, [0.0], HeldTorque::Body(Vec3::new(0.0, 0.0, 8.0)))
            .is_err()
    );
    assert_eq!(initial, before);
}

#[test]
fn zero_drag_preserves_native_motion_and_invalid_or_unresolved_stages_are_atomic() {
    let bank = oblique();
    let state = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.3, 0.1)),
            Vec3::new(-0.4, 0.3, 0.2),
            [0.7, -0.2, 1.1],
            [0.1, 0.0, -0.3],
        )
        .unwrap();
    let before = state;
    let drive = [0.3, -0.7, 0.2];
    for torque in [
        HeldTorque::Body(Vec3::new(0.1, 0.4, -0.2)),
        HeldTorque::Space(Vec3::new(0.1, 0.4, -0.2)),
    ] {
        let plain = bank.step(&state, drive, torque, 0.04).unwrap();
        let zero = bank
            .viscous_stepper([0.0; 3], 0.04)
            .unwrap()
            .step(&state, drive, torque)
            .unwrap();
        assert_eq!(zero.motion, plain);
        assert_eq!(zero.dissipated_energy_j, [0.0; 3]);
        close(zero.work_energy_residual_j, 0.0, 3e-14);
    }
    for interval in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(bank.viscous_stepper([0.1; 3], interval).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(bank.viscous_stepper([0.1, bad, 0.2], 0.1).is_err());
    }
    let stepper = bank.viscous_stepper([0.1; 3], 0.1).unwrap();
    assert!(
        stepper
            .step(&state, [0.0, f64::NAN, 0.0], zero_torque())
            .is_err()
    );
    let mut phase = state;
    phase.relative_angle_rad[0] = 1e30;
    assert!(stepper.step(&phase, [0.0; 3], zero_torque()).is_err());
    assert_eq!(state, before);
    assert_eq!(phase.relative_angle_rad[0], 1e30);
    // Exact zero-drag compatibility also keeps the older held-impulse lane's
    // explicitly reported coarse rounding, rather than changing its admission.
    let axis_bank = coaxial();
    let coarse = axis_bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16], [0.0])
        .unwrap();
    let native = axis_bank.step(&coarse, [2.5], zero_torque(), 1.0).unwrap();
    let delegated = axis_bank
        .viscous_stepper([0.0], 1.0)
        .unwrap()
        .step(&coarse, [2.5], zero_torque())
        .unwrap();
    assert_eq!(delegated.motion, native);
    assert_eq!(delegated.drive_work_j, native.shaft_work_j);
    // Retain the finite-transfer aggregate-h refusal under the implicit law.
    let rotor = AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let large = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap(),
        [rotor; 2],
    )
    .unwrap();
    let spin = large
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16, 0.0], [0.0; 2])
        .unwrap();
    let saved = spin;
    assert!(
        large
            .viscous_stepper([0.0, 1.0], 1.0)
            .unwrap()
            .step(&spin, [0.0, 2.0], zero_torque())
            .is_err()
    );
    assert_eq!(spin, saved);
    let empty = RotorBank::from_reduced(bank.reduced_inertia(), []).unwrap();
    let empty_state = empty
        .state_from_velocities(state.rotation, Vec3::new(0.3, 0.2, -0.1), [], [])
        .unwrap();
    assert_eq!(
        empty
            .viscous_stepper([], 0.1)
            .unwrap()
            .step(&empty_state, [], zero_torque())
            .unwrap()
            .motion,
        empty.step(&empty_state, [], zero_torque(), 0.1).unwrap()
    );
}
