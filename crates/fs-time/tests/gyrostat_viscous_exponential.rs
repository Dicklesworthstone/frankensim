//! G1/G3: original shaft laws, independent continuous work and Euler dynamics.

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{AxialRotor, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tol={tolerance:e}"
    );
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn coaxial<const N: usize>() -> RotorBank<N> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap(); N],
    )
    .unwrap()
}

// Integrate r=steady+(initial-steady)*exp(-rate*t), independently of native
// response coefficients. exp_m1 retains small losses in the weak-drag case.
fn closed(initial: f64, steady: f64, rate: f64, time: f64) -> (f64, f64, f64) {
    let q = initial - steady;
    let first = -(-rate * time).exp_m1() / rate;
    let second = -(-2.0 * rate * time).exp_m1() / (2.0 * rate);
    (
        steady + q * (-rate * time).exp(),
        steady * time + q * first,
        steady * steady * time + 2.0 * steady * q * first + q * q * second,
    )
}

#[test]
fn stiff_and_weak_scalar_decay_integrate_original_phase_and_heat() {
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [3.0], [0.0])
        .unwrap();
    for damping in [2.0, 1e3, 1e100] {
        let time = 0.8;
        let rate = 1.5 * damping;
        let prepared = bank.viscous_exponential_stepper([damping], time).unwrap();
        assert_eq!(prepared.damped_rotor_count(), 1);
        let step = prepared
            .step(&initial, [0.0], zero())
            .unwrap_or_else(|error| panic!("damping={damping:e}: {error:?}"));
        let (relative, phase, square_integral) = closed(3.0, 0.0, rate, time);
        close(
            step.state_after.axial_momentum_nms[0],
            1.0 + relative / 1.5,
            3e-14,
        );
        close(
            step.state_after.relative_angle_rad[0],
            phase,
            phase.abs() * 3e-14,
        );
        close(
            step.dissipated_energy_j,
            damping * square_integral,
            3e-14 * damping * square_integral,
        );
        close(step.work_energy_residual_j, 0.0, 4e-14);
        assert!(step.dissipated_energy_j > 0.0);
        if damping >= 1e3 {
            assert!(bank.viscous_stepper([damping], time).is_err());
        }
    }
    // A weak drag still produces a resolved, separately integrated heat even
    // when it barely changes a driven trajectory. Use independent positive
    // Simpson integration of the scalar exponential for this tiny loss.
    let damping = 1e-8;
    let drive = 0.5;
    let time = 0.8;
    let rate = 1.5 * damping;
    let rest = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [0.0], [0.0])
        .unwrap();
    let step = bank
        .viscous_exponential_stepper([damping], time)
        .unwrap()
        .step(&rest, [drive], zero())
        .unwrap();
    let mut square_integral = 0.0;
    let mut phase = 0.0;
    for i in 0..=2048 {
        let t = time * i as f64 / 2048.0;
        let relative = -(drive / damping) * (-rate * t).exp_m1();
        let weight = if i == 0 || i == 2048 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        square_integral += weight * relative * relative;
        phase += weight * relative;
    }
    let dt_over_three = time / (3.0 * 2048.0);
    close(
        step.dissipated_energy_j,
        damping * square_integral * dt_over_three,
        2e-22,
    );
    close(
        step.state_after.relative_angle_rad[0],
        phase * dt_over_three,
        3e-14,
    );
    close(step.work_energy_residual_j, 0.0, 3e-14);
}

#[test]
fn held_drive_and_external_carrier_torque_have_independent_continuous_work() {
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.1, 0.3)),
            Vec3::default(),
            [3.0],
            [0.0],
        )
        .unwrap();
    let damping = 1e4;
    let drive = 0.4;
    let external = 0.3;
    let time = 0.8;
    let (relative, phase, square_integral) =
        closed(3.0, (drive - external / 3.0) / damping, 1.5 * damping, time);
    let step = bank
        .viscous_exponential_stepper([damping], time)
        .unwrap()
        .step(
            &initial,
            [drive],
            HeldTorque::Body(Vec3::new(0.0, 0.0, external)),
        )
        .unwrap();
    let total = initial.total_momentum_body.z + external * time;
    close(
        step.state_after.axial_momentum_nms[0],
        (relative + total / 2.0) / 1.5,
        2e-14,
    );
    close(step.state_after.relative_angle_rad[0], phase, 2e-15);
    close(step.drive_work_j[0], drive * phase, 2e-15);
    close(step.dissipated_energy_j, damping * square_integral, 4e-14);
    let integrated_carrier =
        (initial.total_momentum_body.z * time + 0.5 * external * time * time - phase) / 3.0;
    close(step.external_work_j, external * integrated_carrier, 3e-14);
    close(
        step.net_shaft_work_j,
        drive * phase - damping * square_integral,
        4e-14,
    );
    close(step.work_energy_residual_j, 0.0, 4e-14);

    // At exact driven equilibrium the drag and motor remain nonzero, while
    // their net impulse vanishes. No exp/phi cancellation may invent a kick.
    let equilibrium = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.5), [1.5], [0.0])
        .unwrap();
    let steady = bank
        .viscous_exponential_stepper([4.0], 100.0)
        .unwrap()
        .step(&equilibrium, [6.0], zero())
        .unwrap();
    assert_eq!(
        steady.state_after.axial_momentum_nms,
        equilibrium.axial_momentum_nms
    );
    close(steady.state_after.relative_angle_rad[0], 150.0, 3e-12);
    close(steady.drive_work_j[0], 900.0, 2e-11);
    close(steady.dissipated_energy_j, 900.0, 2e-11);
    close(steady.work_energy_residual_j, 0.0, 2e-11);
}

#[test]
fn coupled_parallel_shafts_retain_their_distinct_collective_decay_rates() {
    let bank = coaxial::<2>();
    let state = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [3.0, 1.0], [0.0; 2])
        .unwrap();
    let damping = 20.0;
    let time = 0.8;
    let common = (-2.0_f64 * damping * time).exp();
    let difference = (-damping * time).exp();
    let step = bank
        .viscous_exponential_stepper([damping; 2], time)
        .unwrap()
        .step(&state, [0.0; 2], zero())
        .unwrap();
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        close(
            step.state_after.axial_momentum_nms[i],
            1.0 + common + sign * difference,
            3e-14,
        );
        close(
            step.state_after.relative_angle_rad[i],
            ((1.0 - common) + sign * (1.0 - difference)) / damping,
            3e-14,
        );
    }
    let heat = 3.0 - 2.0 * common * common - difference * difference;
    close(step.dissipated_energy_j, heat, 5e-14);
    close(step.work_energy_residual_j, 0.0, 5e-14);
}

#[test]
fn undamped_driven_shafts_still_force_the_damped_carrier_coupling() {
    let bank = coaxial::<2>();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.0, 0.0, 0.1),
            [2.0, -0.3],
            [0.0; 2],
        )
        .unwrap();
    let time = 0.3;
    let damping = 50.0;
    let drive = [0.2, 0.5];
    let external = 0.1;
    let (relative, phase_first, square_integral) =
        closed(2.0, 1.0 / (3.0 * damping), 1.5 * damping, time);
    let prepared = bank
        .viscous_exponential_stepper([damping, 0.0], time)
        .unwrap();
    assert_eq!(prepared.damped_rotor_count(), 1);
    let step = prepared
        .step(
            &initial,
            drive,
            HeldTorque::Body(Vec3::new(0.0, 0.0, external)),
        )
        .unwrap();
    let p2 = initial.axial_momentum_nms[1] + drive[1] * time;
    let total = initial.total_momentum_body.z + external * time;
    let p1 = (relative + 0.5 * total - 0.5 * p2) / 1.5;
    let integral_p2 = initial.axial_momentum_nms[1] * time + 0.5 * drive[1] * time * time;
    let integral_total = initial.total_momentum_body.z * time + 0.5 * external * time * time;
    let phase_second = (phase_first + 4.0 * integral_p2 - integral_total) / 3.0;
    close(step.state_after.axial_momentum_nms[0], p1, 3e-14);
    close(step.state_after.axial_momentum_nms[1], p2, 3e-14);
    close(step.state_after.relative_angle_rad[0], phase_first, 3e-14);
    close(step.state_after.relative_angle_rad[1], phase_second, 3e-14);
    close(step.drive_work_j[0], drive[0] * phase_first, 3e-14);
    close(step.drive_work_j[1], drive[1] * phase_second, 3e-14);
    close(step.dissipated_energy_j, damping * square_integral, 3e-14);
    close(
        step.external_work_j,
        external * (integral_total - integral_p2 - phase_first) / 3.0,
        3e-14,
    );
    close(step.work_energy_residual_j, 0.0, 3e-14);
}

fn matrix() -> Mat3 {
    Mat3 {
        m: [2.4, 0.3, -0.2, 0.3, 3.1, 0.25, -0.2, 0.25, 4.0],
    }
}
fn oblique(scale: f64, frame: So3) -> RotorBank<3> {
    let r = frame.matrix();
    let transformed = r.compose(matrix()).compose(r.transpose());
    RotorBank::from_reduced(
        InertiaTensor3::new(Mat3 {
            m: transformed.m.map(|v| v * scale),
        })
        .unwrap(),
        [
            AxialRotor::new(frame.rotate(Vec3::new(1.0, 0.2, 0.0)).unwrap(), 0.4 * scale).unwrap(),
            AxialRotor::new(frame.rotate(Vec3::new(0.3, 1.0, 0.4)).unwrap(), 0.7 * scale).unwrap(),
            AxialRotor::new(
                frame.rotate(Vec3::new(-0.2, 0.1, 1.0)).unwrap(),
                0.9 * scale,
            )
            .unwrap(),
        ],
    )
    .unwrap()
}
fn inverse(m: Mat3, v: Vec3) -> Vec3 {
    let [a, b, c, _, d, e, _, _, f] = m.m;
    let adj = Mat3 {
        m: [
            d * f - e * e,
            c * e - b * f,
            b * e - c * d,
            c * e - b * f,
            a * f - c * c,
            b * c - a * e,
            b * e - c * d,
            b * c - a * e,
            a * d - b * b,
        ],
    };
    let determinant = a * adj.m[0] + b * adj.m[1] + c * adj.m[2];
    adj.apply(v).scale(1.0 / determinant)
}
type Reference = [f64; 18]; // Quaternion, total momentum, shaft momenta/phases, heat, drive works, external work.

fn pack(state: &RotorState<3>) -> Reference {
    let q = state.rotation.as_quat();
    let mut y = [0.0; 18];
    y[..4].copy_from_slice(&[q.w, q.x, q.y, q.z]);
    y[4..7].copy_from_slice(&[
        state.total_momentum_body.x,
        state.total_momentum_body.y,
        state.total_momentum_body.z,
    ]);
    y[7..10].copy_from_slice(&state.axial_momentum_nms);
    y[10..13].copy_from_slice(&state.relative_angle_rad);
    y
}
fn derivative(bank: &RotorBank<3>, state: Reference) -> Reference {
    let total = Vec3::new(state[4], state[5], state[6]);
    let mut carrier = total;
    for i in 0..3 {
        carrier = carrier - bank.rotors()[i].axis_body().scale(state[7 + i]);
    }
    let omega = inverse(matrix(), carrier);
    // Independent quaternion transport for an external torque held in space.
    let imaginary = Vec3::new(state[1], state[2], state[3]);
    let world = Vec3::new(0.1, 0.2, -0.1);
    let external = world.scale(state[0] * state[0] - imaginary.dot(imaginary))
        + imaginary.scale(2.0 * imaginary.dot(world))
        - imaginary.cross(world).scale(2.0 * state[0]);
    let qdot = omega.scale(state[0]) + imaginary.cross(omega);
    let pdot = total.cross(omega) + external;
    let mut out = [0.0; 18];
    out[0] = -0.5 * imaginary.dot(omega);
    out[1..4].copy_from_slice(&[0.5 * qdot.x, 0.5 * qdot.y, 0.5 * qdot.z]);
    out[4..7].copy_from_slice(&[pdot.x, pdot.y, pdot.z]);
    for i in 0..3 {
        let rotor = bank.rotors()[i];
        let rate = state[7 + i] / rotor.axial_inertia_kg_m2() - rotor.axis_body().dot(omega);
        let drive = [0.1, -0.2, 0.4][i];
        let damping = [0.8, 1.7, 0.4][i];
        out[7 + i] = drive - damping * rate;
        out[10 + i] = rate;
        out[13] += damping * rate * rate;
        out[14 + i] = drive * rate;
    }
    out[17] = external.dot(omega);
    out
}

#[test]
fn full_tensor_rotation_phase_and_physical_work_refine_against_independent_euler_rk4() {
    let bank = oblique(1.0, So3::identity());
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.3, 0.1)),
            Vec3::new(0.6, -0.3, 0.8),
            [1.0, -0.5, 0.7],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let time = 0.7;
    let dt = time / 8192.0;
    let mut reference = pack(&initial);
    for _ in 0..8192 {
        let a = derivative(&bank, reference);
        let b = derivative(
            &bank,
            core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]),
        );
        let c = derivative(
            &bank,
            core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]),
        );
        let d = derivative(&bank, core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let stepper = bank
            .viscous_exponential_stepper([0.8, 1.7, 0.4], time / count as f64)
            .unwrap();
        let mut state = initial;
        let mut accounts = [0.0; 5];
        for _ in 0..count {
            let step = stepper
                .step(
                    &state,
                    [0.1, -0.2, 0.4],
                    HeldTorque::Space(Vec3::new(0.1, 0.2, -0.1)),
                )
                .unwrap();
            assert!(step.dissipated_energy_j >= 0.0);
            close(step.work_energy_residual_j, 0.0, 3e-13);
            assert!(step.angular_momentum_residual_world.norm() < 3e-13);
            accounts[0] += step.dissipated_energy_j;
            for i in 0..3 {
                accounts[1 + i] += step.drive_work_j[i];
            }
            accounts[4] += step.external_work_j;
            state = step.state_after;
        }
        let mut actual = pack(&state);
        actual[13..].copy_from_slice(&accounts);
        let error = actual
            .iter()
            .zip(reference)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        errors.push(error);
    }
    eprintln!("independent Euler/quaternion/heat refinement: {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2);
    }
}

#[test]
fn original_material_frame_and_uniform_inertia_units_preserve_the_coupled_response() {
    let frame = rotation(Vec3::new(0.2, -0.4, 0.3));
    let base = oblique(1.0, So3::identity());
    let attitude = rotation(Vec3::new(-0.3, 0.1, 0.2));
    let omega = Vec3::new(0.4, -0.6, 0.2);
    let relative = [0.7, -0.2, 0.4];
    let phase = [0.0; 3];
    let damping = [1.0, 0.5, 0.2];
    let drive = [0.1, -0.2, 0.3];
    let external = Vec3::new(-0.1, 0.2, 0.3);
    let reference = base
        .viscous_exponential_stepper(damping, 0.05)
        .unwrap()
        .step(
            &base
                .state_from_velocities(attitude, omega, relative, phase)
                .unwrap(),
            drive,
            HeldTorque::Body(external),
        )
        .unwrap();
    for scale in [1e-150, 1.0, 1e150] {
        let bank = oblique(scale, frame);
        let state = bank
            .state_from_velocities(
                attitude.compose(frame.inverse()).unwrap(),
                frame.rotate(omega).unwrap(),
                relative,
                phase,
            )
            .unwrap();
        let result = bank
            .viscous_exponential_stepper(damping.map(|x| scale * x), 0.05)
            .unwrap()
            .step(
                &state,
                drive.map(|x| scale * x),
                HeldTorque::Body(frame.rotate(external).unwrap().scale(scale)),
            )
            .unwrap();
        for i in 0..3 {
            close(
                result.state_after.axial_momentum_nms[i] / scale,
                reference.state_after.axial_momentum_nms[i],
                4e-13,
            );
            close(
                result.state_after.relative_angle_rad[i],
                reference.state_after.relative_angle_rad[i],
                4e-13,
            );
            close(
                result.drive_work_j[i] / scale,
                reference.drive_work_j[i],
                4e-13,
            );
        }
        close(
            result.dissipated_energy_j / scale,
            reference.dissipated_energy_j,
            4e-13,
        );
        close(
            result.external_work_j / scale,
            reference.external_work_j,
            4e-13,
        );
        close(result.work_energy_residual_j / scale, 0.0, 4e-13);
        let world = result
            .state_after
            .rotation
            .rotate(result.state_after.total_momentum_body)
            .unwrap()
            .scale(1.0 / scale);
        let expected_world = reference
            .state_after
            .rotation
            .rotate(reference.state_after.total_momentum_body)
            .unwrap();
        assert!((world - expected_world).norm() < 4e-13);
    }
}

#[test]
fn zero_drag_capacity_and_unresolved_inputs_keep_native_admission() {
    let bank = oblique(1.0, So3::identity());
    let state = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(0.3, -0.2, 0.1),
            [1.0, -0.4, 0.2],
            [0.0; 3],
        )
        .unwrap();
    let drive = [0.1, 0.2, -0.3];
    let external = HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3));
    let original = bank.step(&state, drive, external, 0.2).unwrap();
    let zero_stepper = bank.viscous_exponential_stepper([0.0; 3], 0.2).unwrap();
    assert_eq!(zero_stepper.damped_rotor_count(), 0);
    let result = zero_stepper.step(&state, drive, external).unwrap();
    assert_eq!(result.state_after, original.state_after);
    assert_eq!(
        result.axial_momentum_change_nms,
        original.axial_momentum_change_nms
    );
    assert_eq!(result.drive_work_j, original.shaft_work_j);
    assert_eq!(
        result.external_work_j.to_bits(),
        original.external_work_j.to_bits()
    );
    assert_eq!(
        result.work_energy_residual_j.to_bits(),
        original.work_energy_residual_j.to_bits()
    );
    assert_eq!(result.dissipated_energy_j.to_bits(), 0.0_f64.to_bits());
    for time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(bank.viscous_exponential_stepper([1.0; 3], time).is_err());
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            bank.viscous_exponential_stepper([bad, 0.0, 0.1], 0.1)
                .is_err()
        );
    }
    let prepared = bank.viscous_exponential_stepper([1.0; 3], 0.1).unwrap();
    assert!(prepared.step(&state, [f64::NAN; 3], zero()).is_err());
    let mut unresolved = state;
    unresolved.relative_angle_rad[0] = 1e30;
    let saved = unresolved;
    assert!(prepared.step(&unresolved, [0.0; 3], zero()).is_err());
    assert_eq!(unresolved, saved);
    assert!(prepared.step(&state, [0.0; 3], zero()).is_ok());

    // The new physical work account cannot bypass the original lost-recoil
    // check when a large stored rotor erases another shaft's finite transfer.
    let large = coaxial::<2>();
    let large_state = large
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16, 0.0], [0.0; 2])
        .unwrap();
    assert!(
        large
            .viscous_exponential_stepper([0.0, 1.0], 1.0)
            .unwrap()
            .step(&large_state, [0.0, 2.0], zero())
            .is_err()
    );
    let single = coaxial::<1>();
    let fine = single
        .state_from_velocities(So3::identity(), Vec3::default(), [3.0], [0.0])
        .unwrap();
    assert!(
        single
            .viscous_exponential_stepper([1e-20], 0.8)
            .unwrap()
            .step(&fine, [0.0], zero())
            .is_err()
    );

    let empty = coaxial::<0>();
    let empty_state = empty
        .state_from_velocities(So3::identity(), Vec3::new(0.3, 0.2, -0.1), [], [])
        .unwrap();
    assert_eq!(
        empty
            .viscous_exponential_stepper([], 0.1)
            .unwrap()
            .step(&empty_state, [], zero())
            .unwrap()
            .state_after,
        empty
            .step(&empty_state, [], zero(), 0.1)
            .unwrap()
            .state_after
    );
    let largest = coaxial::<64>();
    let largest_state = largest
        .state_from_velocities(So3::identity(), Vec3::default(), [1.0; 64], [0.0; 64])
        .unwrap();
    let largest_stepper = largest
        .viscous_exponential_stepper([0.03; 64], 0.01)
        .unwrap();
    assert_eq!(largest_stepper.damped_rotor_count(), 64);
    assert!(
        largest_stepper
            .step(&largest_state, [0.0; 64], zero())
            .unwrap()
            .dissipated_energy_j
            > 0.0
    );
    assert!(RotorBank::from_reduced(largest.reduced_inertia(), [largest.rotors()[0]; 65]).is_err());
}
