//! Independent mechanics for finite ideal gearing and carrier momentum.

use fs_ga::{Quat, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{AxialRotor, GearConstraint, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn turn(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn rotation_error(a: So3, b: So3) -> f64 {
    [v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)]
        .into_iter()
        .map(|axis| (a.rotate(axis).unwrap() - b.rotate(axis).unwrap()).norm())
        .fold(0.0, f64::max)
}
fn relation(first: usize, second: usize, ratio: f64) -> GearConstraint {
    GearConstraint {
        first_rotor: first,
        second_rotor: second,
        ratio,
    }
}
fn coaxial() -> RotorBank<2> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(v(2.0, 3.0, 4.0), So3::identity()).unwrap(),
        [
            AxialRotor::new(v(0.0, 0.0, 1.0), 0.3).unwrap(),
            AxialRotor::new(v(0.0, 0.0, 1.0), 0.7).unwrap(),
        ],
    )
    .unwrap()
}

#[test]
fn opposed_gears_match_independent_two_coordinate_driven_solution() {
    let bank = coaxial();
    let train = bank.gear_train([relation(0, 1, -2.0)]).unwrap();
    // Independent original kinetic matrix on carrier angle and q_1:
    // [[I+J0+J1, -2J0+J1],[-2J0+J1,4J0+J1]] = [[5,.1],[.1,1.9]].
    let determinant = 5.0 * 1.9 - 0.1 * 0.1;
    let torque = [1.1, -0.2];
    let generalized_drive = -2.0 * torque[0] + torque[1];
    for external in [0.0, 0.3] {
        let omega = 0.4;
        let rate = 3.0;
        let initial = bank
            .state_from_velocities(
                So3::identity(),
                v(0.0, 0.0, omega),
                [-2.0 * rate, rate],
                [0.7, -0.4],
            )
            .unwrap();
        let dt = 0.4;
        let omega_dot = (1.9 * external - 0.1 * generalized_drive) / determinant;
        let rate_dot = (5.0 * generalized_drive - 0.1 * external) / determinant;
        let result = train
            .step(
                &initial,
                torque,
                HeldTorque::Space(v(0.0, 0.0, external)),
                dt,
            )
            .unwrap();
        let after = result.motion.state_after;
        let expected_omega = omega + omega_dot * dt;
        let expected_rate = rate + rate_dot * dt;
        vector_near(
            bank.angular_velocity(&after).unwrap(),
            v(0.0, 0.0, expected_omega),
            3e-14,
        );
        let rates = bank.relative_rates(&after).unwrap();
        near(rates[0], -2.0 * expected_rate, 3e-14);
        near(rates[1], expected_rate, 3e-14);
        near(
            after.axial_momentum_nms[0],
            0.3 * (expected_omega - 2.0 * expected_rate),
            3e-14,
        );
        near(
            after.axial_momentum_nms[1],
            0.7 * (expected_omega + expected_rate),
            3e-14,
        );
        let angle = rate * dt + 0.5 * rate_dot * dt * dt;
        near(
            after.relative_angle_rad[0] - initial.relative_angle_rad[0],
            -2.0 * angle,
            3e-14,
        );
        near(
            after.relative_angle_rad[1] - initial.relative_angle_rad[1],
            angle,
            3e-14,
        );
        let carrier_angle = omega * dt + 0.5 * omega_dot * dt * dt;
        assert!(rotation_error(after.rotation, turn(v(0.0, 0.0, carrier_angle))) < 3e-14);
        for i in 0..2 {
            near(
                result.drive_work_j[i],
                torque[i] * (after.relative_angle_rad[i] - initial.relative_angle_rad[i]),
                4e-14,
            );
        }
        near(
            result.motion.external_work_j,
            external * carrier_angle,
            4e-14,
        );
        near(result.reaction_work_j, 0.0, 5e-14);
        near(result.work_energy_residual_j, 0.0, 2e-13);
        let mesh = result.constraint_impulse_nms[0];
        near(
            result.motion.axial_momentum_change_nms[0],
            torque[0] * dt + mesh,
            4e-14,
        );
        near(
            result.motion.axial_momentum_change_nms[1],
            torque[1] * dt + 2.0 * mesh,
            4e-14,
        );
        vector_near(
            after.rotation.rotate(after.total_momentum_body).unwrap()
                - initial
                    .rotation
                    .rotate(initial.total_momentum_body)
                    .unwrap(),
            v(0.0, 0.0, external * dt),
            5e-14,
        );
    }
}

#[test]
fn engagement_keeps_generalized_spin_and_accounts_for_only_the_actual_capture_loss() {
    let bank = coaxial();
    let train = bank.gear_train([relation(0, 1, -2.0)]).unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.0, 0.0, 0.4),
            [-3.0, 5.0],
            [0.7, -0.4],
        )
        .unwrap();
    let total = initial.total_momentum_body.z;
    let generalized = -2.0 * initial.axial_momentum_nms[0] + initial.axial_momentum_nms[1];
    let det = 5.0 * 1.9 - 0.01;
    let omega = (1.9 * total - 0.1 * generalized) / det;
    let rate = (5.0 * generalized - 0.1 * total) / det;
    let captured = train.engage(&initial).unwrap();
    let after = captured.impulse.state_after;
    assert_eq!(after.rotation, initial.rotation);
    assert_eq!(after.total_momentum_body, initial.total_momentum_body);
    assert_eq!(after.relative_angle_rad, initial.relative_angle_rad);
    near(bank.angular_velocity(&after).unwrap().z, omega, 2e-14);
    let rates = bank.relative_rates(&after).unwrap();
    near(rates[0], -2.0 * rate, 2e-14);
    near(rates[1], rate, 2e-14);
    near(
        -2.0 * after.axial_momentum_nms[0] + after.axial_momentum_nms[1],
        generalized,
        3e-14,
    );
    // Evaluate the original quadratic energy, independently of receipt helpers.
    let energy_before = 0.5 * 4.0 * 0.4_f64.powi(2)
        + initial.axial_momentum_nms[0].powi(2) / (2.0 * 0.3)
        + initial.axial_momentum_nms[1].powi(2) / (2.0 * 0.7);
    let energy_after = 0.5 * 5.0 * omega * omega + 0.1 * omega * rate + 0.5 * 1.9 * rate * rate;
    near(
        captured.dissipated_energy_j,
        energy_before - energy_after,
        1e-13,
    );
    assert!(captured.dissipated_energy_j > 0.0);
    near(captured.work_energy_residual_j, 0.0, 1e-13);
    near(
        captured.impulse.axial_momentum_change_nms[1],
        2.0 * captured.impulse.axial_momentum_change_nms[0],
        5e-14,
    );
    let repeated = train.engage(&after).unwrap();
    assert_eq!(repeated.impulse.state_after, after);
    assert_eq!(repeated.dissipated_energy_j, 0.0);
}

fn spatial_bank() -> RotorBank<3> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(v(2.4, 3.2, 4.1), turn(v(0.2, -0.3, 0.4))).unwrap(),
        [
            AxialRotor::new(v(1.0, 0.2, 0.1), 0.25).unwrap(),
            AxialRotor::new(v(-0.3, 1.0, 0.2), 0.4).unwrap(),
            AxialRotor::new(v(0.2, -0.1, 1.0), 0.35).unwrap(),
        ],
    )
    .unwrap()
}
fn spatial_state(bank: &RotorBank<3>) -> RotorState<3> {
    bank.state_from_velocities(
        turn(v(-0.15, 0.25, 0.3)),
        v(0.2, -0.35, 0.4),
        [4.0, -2.0, -4.0],
        [0.0; 3],
    )
    .unwrap()
}

#[test]
fn nonparallel_gears_preserve_spatial_momentum_and_reverse_driven_precession() {
    let bank = spatial_bank();
    let train = bank
        .gear_train([relation(0, 1, -2.0), relation(1, 2, 0.5)])
        .unwrap();
    let initial = spatial_state(&bank);
    let spin = initial
        .rotation
        .rotate(initial.total_momentum_body)
        .unwrap();
    let generalized = initial.axial_momentum_nms[0]
        - 0.5 * initial.axial_momentum_nms[1]
        - initial.axial_momentum_nms[2];
    let mut state = initial;
    for _ in 0..100 {
        let next = train.drift(&state, 0.01).unwrap();
        state = next.state_after;
        vector_near(
            state.rotation.rotate(state.total_momentum_body).unwrap(),
            spin,
            3e-12,
        );
        near(
            state.axial_momentum_nms[0]
                - 0.5 * state.axial_momentum_nms[1]
                - state.axial_momentum_nms[2],
            generalized,
            2e-12,
        );
        let rates = bank.relative_rates(&state).unwrap();
        near(rates[0], -2.0 * rates[1], 3e-12);
        near(rates[1], 0.5 * rates[2], 3e-12);
        near(
            state.relative_angle_rad[0],
            -2.0 * state.relative_angle_rad[1],
            3e-12,
        );
        near(
            state.relative_angle_rad[1],
            0.5 * state.relative_angle_rad[2],
            3e-12,
        );
        near(next.work_energy_residual_j, 0.0, 2e-13);
    }
    for _ in 0..100 {
        state = train.drift(&state, -0.01).unwrap().state_after;
    }
    assert!(rotation_error(state.rotation, initial.rotation) < 4e-12);
    vector_near(
        state.total_momentum_body,
        initial.total_momentum_body,
        4e-12,
    );
    for i in 0..3 {
        near(
            state.axial_momentum_nms[i],
            initial.axial_momentum_nms[i],
            4e-12,
        );
        near(state.relative_angle_rad[i], 0.0, 4e-12);
    }
    let drive = [0.4, -0.2, 0.1];
    for torque in [
        HeldTorque::Body(v(0.1, -0.2, 0.3)),
        HeldTorque::Space(v(0.1, -0.2, 0.3)),
    ] {
        let result = train.step(&initial, drive, torque, 0.04).unwrap();
        vector_near(
            result.motion.angular_momentum_residual_world,
            Vec3::default(),
            2e-13,
        );
        near(result.reaction_work_j, 0.0, 2e-13);
        near(result.work_energy_residual_j, 0.0, 2e-13);
        for i in 0..3 {
            near(
                result.drive_work_j[i],
                drive[i] * result.motion.state_after.relative_angle_rad[i],
                2e-13,
            );
        }
        let back = train
            .step(&result.motion.state_after, drive, torque, -0.04)
            .unwrap()
            .motion
            .state_after;
        assert!(rotation_error(back.rotation, initial.rotation) < 3e-12);
        vector_near(back.total_momentum_body, initial.total_momentum_body, 3e-12);
        for i in 0..3 {
            near(
                back.axial_momentum_nms[i],
                initial.axial_momentum_nms[i],
                3e-12,
            );
            near(back.relative_angle_rad[i], 0.0, 3e-12);
        }
    }
}

// Independent four-coordinate Euler equations. The constant mass matrix is
// assembled from the ORIGINAL carrier tensor and each absolute rotor kinetic
// energy; the production retained-inertia formula and gear solver are not used.
fn inverse4(mut a: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut inverse = [[0.0; 4]; 4];
    for i in 0..4 {
        inverse[i][i] = 1.0;
    }
    for col in 0..4 {
        let pivot = (col..4)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .unwrap();
        a.swap(col, pivot);
        inverse.swap(col, pivot);
        let d = a[col][col];
        for j in 0..4 {
            a[col][j] /= d;
            inverse[col][j] /= d;
        }
        for row in 0..4 {
            if row == col {
                continue;
            }
            let scale = a[row][col];
            for j in 0..4 {
                a[row][j] -= scale * a[col][j];
                inverse[row][j] -= scale * inverse[col][j];
            }
        }
    }
    inverse
}
fn qproduct(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}
fn original_mass(bank: &RotorBank<3>) -> [[f64; 4]; 4] {
    let mut m = [[0.0; 4]; 4];
    let original = bank.reduced_inertia().matrix().m;
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = original[3 * i + j];
        }
    }
    for (rotor, c) in bank.rotors().iter().zip([1.0, -0.5, -1.0]) {
        let a = rotor.axis_body();
        let row = [a.x, a.y, a.z, c];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += rotor.axial_inertia_kg_m2() * row[i] * row[j];
            }
        }
    }
    m
}
fn ode(y: [f64; 9], mass: [[f64; 4]; 4], inverse: [[f64; 4]; 4]) -> [f64; 9] {
    let state = [y[0], y[1], y[2], y[3]];
    let momentum = mass.map(|row| row.into_iter().zip(state).map(|(a, b)| a * b).sum::<f64>());
    let omega = v(y[0], y[1], y[2]);
    let euler = v(momentum[0], momentum[1], momentum[2]).cross(omega) + v(0.1, -0.2, 0.3);
    let rhs = [euler.x, euler.y, euler.z, 0.4 - 0.5 * (-0.2) - 0.1];
    let acceleration = inverse.map(|row| row.into_iter().zip(rhs).map(|(a, b)| a * b).sum::<f64>());
    let qdot = qproduct([y[4], y[5], y[6], y[7]], [0.0, y[0], y[1], y[2]]);
    [
        acceleration[0],
        acceleration[1],
        acceleration[2],
        acceleration[3],
        0.5 * qdot[0],
        0.5 * qdot[1],
        0.5 * qdot[2],
        0.5 * qdot[3],
        y[3],
    ]
}
fn reference(bank: &RotorBank<3>, initial: RotorState<3>, duration: f64) -> [f64; 9] {
    let mass = original_mass(bank);
    let inverse = inverse4(mass);
    let omega = bank.angular_velocity(&initial).unwrap();
    let rate = bank.relative_rates(&initial).unwrap()[0];
    let q = initial.rotation.as_quat();
    let mut y = [omega.x, omega.y, omega.z, rate, q.w, q.x, q.y, q.z, 0.0];
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let a = ode(y, mass, inverse);
        let b = ode(
            core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]),
            mass,
            inverse,
        );
        let c = ode(
            core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]),
            mass,
            inverse,
        );
        let d = ode(core::array::from_fn(|i| y[i] + dt * c[i]), mass, inverse);
        for i in 0..9 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn nonparallel_loaded_gears_converge_to_independent_original_kinetic_equations() {
    let bank = spatial_bank();
    let train = bank
        .gear_train([relation(0, 1, -2.0), relation(1, 2, 0.5)])
        .unwrap();
    let initial = spatial_state(&bank);
    let reference = reference(&bank, initial, 1.0);
    let expected = So3::try_from_quat(Quat {
        w: reference[4],
        x: reference[5],
        y: reference[6],
        z: reference[7],
    })
    .unwrap();
    let mut errors = [0.0; 3];
    for (sample, count) in [20, 40, 80].into_iter().enumerate() {
        let mut state = initial;
        for _ in 0..count {
            state = train
                .step(
                    &state,
                    [0.4, -0.2, 0.1],
                    HeldTorque::Body(v(0.1, -0.2, 0.3)),
                    1.0 / count as f64,
                )
                .unwrap()
                .motion
                .state_after;
        }
        let omega = bank.angular_velocity(&state).unwrap();
        let rate = bank.relative_rates(&state).unwrap()[0];
        errors[sample] = (omega - v(reference[0], reference[1], reference[2])).norm()
            + (rate - reference[3]).abs()
            + rotation_error(state.rotation, expected)
            + (state.relative_angle_rad[0] - reference[8]).abs();
    }
    for i in 0..2 {
        let ratio = errors[i] / errors[i + 1];
        assert!((3.7..4.3).contains(&ratio), "errors {errors:?}");
    }
}

#[test]
fn signed_ratio_orientation_and_constraint_order_preserve_the_same_physics() {
    let bank = spatial_bank();
    let forward = bank
        .gear_train([relation(0, 1, -2.0), relation(1, 2, 0.5)])
        .unwrap();
    let reordered = bank
        .gear_train([relation(2, 1, 2.0), relation(1, 0, -0.5)])
        .unwrap();
    let initial = spatial_state(&bank);
    let a = forward
        .step(
            &initial,
            [0.2, -0.1, 0.4],
            HeldTorque::Space(v(0.1, 0.2, -0.3)),
            0.03,
        )
        .unwrap();
    let b = reordered
        .step(
            &initial,
            [0.2, -0.1, 0.4],
            HeldTorque::Space(v(0.1, 0.2, -0.3)),
            0.03,
        )
        .unwrap();
    assert!(rotation_error(a.motion.state_after.rotation, b.motion.state_after.rotation) < 1e-13);
    vector_near(
        a.motion.state_after.total_momentum_body,
        b.motion.state_after.total_momentum_body,
        1e-13,
    );
    for i in 0..3 {
        near(
            a.motion.state_after.axial_momentum_nms[i],
            b.motion.state_after.axial_momentum_nms[i],
            1e-13,
        );
        near(a.drive_work_j[i], b.drive_work_j[i], 1e-13);
    }
}

#[test]
fn invalid_or_redundant_transmissions_refuse_and_empty_constraints_keep_native_motion() {
    let bank = spatial_bank();
    for ratio in [0.0, f64::INFINITY, f64::NAN] {
        assert!(bank.gear_train([relation(0, 1, ratio)]).is_err());
    }
    assert!(bank.gear_train([relation(0, 0, 1.0)]).is_err());
    assert!(bank.gear_train([relation(0, 3, 1.0)]).is_err());
    assert!(
        bank.gear_train([relation(0, 1, 2.0), relation(1, 0, 0.5)])
            .is_err()
    );
    for ratio in [0.25, 0.3] {
        assert!(
            bank.gear_train([
                relation(0, 1, 2.0),
                relation(1, 2, 2.0),
                relation(2, 0, ratio)
            ])
            .is_err()
        );
    }
    let state = spatial_state(&bank);
    let empty = bank.gear_train([]).unwrap();
    let native = bank
        .step(
            &state,
            [0.1, -0.2, 0.3],
            HeldTorque::Body(v(0.2, 0.1, -0.1)),
            0.01,
        )
        .unwrap();
    let result = empty
        .step(
            &state,
            [0.1, -0.2, 0.3],
            HeldTorque::Body(v(0.2, 0.1, -0.1)),
            0.01,
        )
        .unwrap();
    assert_eq!(result.motion, native);
    let train = bank.gear_train([relation(0, 1, 1.0)]).unwrap();
    let before = state;
    assert!(
        train
            .step(&state, [0.0; 3], HeldTorque::Body(Vec3::default()), 0.01)
            .is_err()
    );
    assert_eq!(state, before);
    let train = bank
        .gear_train([relation(0, 1, -2.0), relation(1, 2, 0.5)])
        .unwrap();
    for dt in [0.0, f64::INFINITY, f64::NAN] {
        assert!(
            train
                .step(&state, [0.0; 3], HeldTorque::Body(Vec3::default()), dt)
                .is_err()
        );
    }
}

#[test]
fn equal_counterrotating_gears_do_not_invent_carrier_spin_or_external_momentum() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(v(2.0, 3.0, 4.0), So3::identity()).unwrap(),
        [AxialRotor::new(v(0.0, 0.0, 1.0), 0.5).unwrap(); 2],
    )
    .unwrap();
    let train = bank.gear_train([relation(0, 1, -1.0)]).unwrap();
    near(train.retained_inertia().matrix().m[8], 5.0, 1e-14);
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            Vec3::default(),
            [8.0, -8.0],
            [0.0; 2],
        )
        .unwrap();
    assert_eq!(initial.total_momentum_body, Vec3::default());
    let result = train
        .step(
            &initial,
            [2.0, -2.0],
            HeldTorque::Body(Vec3::default()),
            0.25,
        )
        .unwrap();
    let after = result.motion.state_after;
    assert_eq!(after.total_momentum_body, Vec3::default());
    assert_eq!(after.rotation, initial.rotation);
    assert_eq!(bank.angular_velocity(&after).unwrap(), Vec3::default());
    near(after.axial_momentum_nms[0], 4.5, 2e-14);
    near(after.axial_momentum_nms[1], -4.5, 2e-14);
    near(after.relative_angle_rad[0], 2.125, 2e-14);
    near(after.relative_angle_rad[1], -2.125, 2e-14);
    near(result.drive_work_j[0] + result.drive_work_j[1], 8.5, 3e-14);
    near(result.constraint_impulse_nms[0], 0.0, 2e-14);
    near(result.work_energy_residual_j, 0.0, 3e-14);
}
