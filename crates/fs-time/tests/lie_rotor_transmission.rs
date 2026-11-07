//! Independent original-coordinate mechanics for multi-shaft transmissions.

use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, GearConstraint, RotorBank, RotorState, ShaftConstraint,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}
fn xyz(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
fn turn(a: Vec3) -> So3 {
    So3::exp(So3Tangent::new(a)).unwrap()
}
fn row<const N: usize>(coefficients: [f64; N]) -> ShaftConstraint<N> {
    ShaftConstraint { coefficients }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}, tolerance {tol:e}"
    );
}
fn vector_near(a: Vec3, b: Vec3, tol: f64) {
    for (a, b) in xyz(a).into_iter().zip(xyz(b)) {
        near(a, b, tol);
    }
}
fn rotation_error(a: So3, b: So3) -> f64 {
    [v(1., 0., 0.), v(0., 1., 0.), v(0., 0., 1.)]
        .into_iter()
        .map(|e| (a.rotate(e).unwrap() - b.rotate(e).unwrap()).norm())
        .fold(0., f64::max)
}
fn bank3(coaxial: bool) -> RotorBank<3> {
    let axes = if coaxial {
        [v(0., 0., 1.); 3]
    } else {
        [v(1., 0., 0.), v(0., 1., 0.), v(0., 0., 1.)]
    };
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(v(4., 5., 7.), So3::identity()).unwrap(),
        core::array::from_fn(|i| AxialRotor::new(axes[i], [0.7, 1.1, 1.5][i]).unwrap()),
    )
    .unwrap()
}
fn bank4() -> RotorBank<4> {
    RotorBank::from_reduced(
        InertiaTensor3::new(Mat3 {
            m: [4., 0.2, -0.1, 0.2, 5., 0.3, -0.1, 0.3, 7.],
        })
        .unwrap(),
        [
            AxialRotor::new(v(1., 0., 0.), 0.8).unwrap(),
            AxialRotor::new(v(0., 1., 0.), 1.3).unwrap(),
            AxialRotor::new(v(0., 0., 1.), 1.9).unwrap(),
            AxialRotor::new(v(1., 1., 1.), 0.6).unwrap(),
        ],
    )
    .unwrap()
}

// Independent fixture solve: partial-pivot elimination of the ORIGINAL kinetic
// matrix in [carrier omega (3), every relative shaft rate (N), row multipliers].
// It never constructs or uses the production reduced tensor or reaction basis.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        assert!(a[pivot][k].abs() > 1e-14);
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..n {
            let ratio = a[i][k] / a[k][k];
            for j in k + 1..n {
                a[i][j] -= ratio * a[k][j];
            }
            b[i] -= ratio * b[k];
        }
    }
    let mut x = vec![0.; n];
    for i in (0..n).rev() {
        x[i] = (b[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
fn original_kkt<const N: usize, const M: usize>(
    bank: &RotorBank<N>,
    rows: &[ShaftConstraint<N>; M],
) -> Vec<Vec<f64>> {
    let mut a = vec![vec![0.; 3 + N + M]; 3 + N + M];
    let inertia = bank.reduced_inertia().matrix().m;
    for r in 0..3 {
        for c in 0..3 {
            a[r][c] = inertia[3 * r + c];
        }
    }
    for (i, rotor) in bank.rotors().iter().enumerate() {
        let axis = xyz(rotor.axis_body());
        let j = rotor.axial_inertia_kg_m2();
        for r in 0..3 {
            for c in 0..3 {
                a[r][c] += j * axis[r] * axis[c];
            }
            a[r][3 + i] = j * axis[r];
            a[3 + i][r] = j * axis[r];
        }
        a[3 + i][3 + i] = j;
        for e in 0..M {
            a[3 + i][3 + N + e] = -rows[e].coefficients[i];
            a[3 + N + e][3 + i] = rows[e].coefficients[i];
        }
    }
    a
}
fn original_delta<const N: usize, const M: usize>(
    bank: &RotorBank<N>,
    rows: &[ShaftConstraint<N>; M],
    drives: [f64; N],
    external: Vec3,
) -> Vec<f64> {
    let mut rhs = vec![0.; 3 + N + M];
    rhs[..3].copy_from_slice(&xyz(external));
    rhs[3..3 + N].copy_from_slice(&drives);
    solve(original_kkt(bank, rows), rhs)
}
fn state_near<const N: usize>(a: &RotorState<N>, b: &RotorState<N>, tol: f64) {
    vector_near(a.total_momentum_body, b.total_momentum_body, tol);
    assert!(rotation_error(a.rotation, b.rotation) < tol);
    for i in 0..N {
        near(a.axial_momentum_nms[i], b.axial_momentum_nms[i], tol);
        near(a.relative_angle_rad[i], b.relative_angle_rad[i], tol);
    }
}

#[test]
fn differential_retains_two_free_output_modes_and_matches_original_kinetic_system() {
    let bank = bank3(true);
    let rows = [row([2., -1., -1.])];
    let transmission = bank.transmission(rows).unwrap();
    let omega = v(0., 0., 0.4);
    let rates = [3., 5., 1.];
    let initial = bank
        .state_from_velocities(So3::identity(), omega, rates, [0.7, -0.4, 0.2])
        .unwrap();
    let drives = [1.2, -0.4, 0.7];
    let external = v(0., 0., 0.3);
    let expected = original_delta(&bank, &rows, drives, external);
    let result = transmission
        .apply_impulse(&initial, drives, external)
        .unwrap();
    let after = result.impulse.state_after;
    vector_near(
        bank.angular_velocity(&after).unwrap(),
        omega + v(expected[0], expected[1], expected[2]),
        5e-14,
    );
    let final_rates = bank.relative_rates(&after).unwrap();
    for i in 0..3 {
        near(final_rates[i], rates[i] + expected[3 + i], 5e-14);
        near(
            result.impulse.axial_momentum_change_nms[i],
            drives[i] + rows[0].coefficients[i] * expected[6],
            5e-14,
        );
    }
    near(result.constraint_impulse_nms[0], expected[6], 3e-14);
    near(
        2. * final_rates[0] - final_rates[1] - final_rates[2],
        0.,
        2e-14,
    );
    assert!((final_rates[1] - final_rates[2]).abs() > 3.);
    assert!((bank.angular_velocity(&after).unwrap().z - omega.z).abs() > 0.1);
    near(result.reaction_work_j, 0., 1e-13);
    near(result.work_energy_residual_j, 0., 2e-13);

    // Coaxial motion has constant acceleration: carrier/shaft phases and work
    // have closed forms, including differential output speed asymmetry.
    let dt = 0.4;
    let step = transmission
        .step(&initial, drives, HeldTorque::Space(external), dt)
        .unwrap();
    let after = step.motion.state_after;
    for i in 0..3 {
        let angle = rates[i] * dt + 0.5 * expected[3 + i] * dt * dt;
        near(
            after.relative_angle_rad[i] - initial.relative_angle_rad[i],
            angle,
            5e-14,
        );
        near(step.drive_work_j[i], drives[i] * angle, 5e-14);
        near(
            bank.relative_rates(&after).unwrap()[i],
            rates[i] + dt * expected[3 + i],
            5e-14,
        );
    }
    let body_angle = omega.z * dt + 0.5 * expected[2] * dt * dt;
    assert!(rotation_error(after.rotation, turn(v(0., 0., body_angle))) < 5e-14);
    near(step.motion.external_work_j, external.z * body_angle, 5e-14);
    near(step.reaction_work_j, 0., 1e-13);
    vector_near(
        step.motion.angular_momentum_residual_world,
        v(0., 0., 0.),
        1e-13,
    );
}

#[test]
fn inelastic_engagement_has_analytic_loss_and_preserves_both_free_momenta() {
    let bank = bank3(false);
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.1, 0.3)),
            v(0.3, -0.2, 0.1),
            [4., -1., 2.],
            [0.7, -0.4, 0.2],
        )
        .unwrap();
    assert!(
        transmission
            .apply_impulse(&initial, [0.; 3], Vec3::default())
            .is_err()
    );
    let mobility = 4. / 0.7 + 1. / 1.1 + 1. / 1.5 + 4. / 4. + 1. / 5. + 1. / 7.;
    let lambda = -(2. * 4. - (-1.) - 2.) / mobility;
    let result = transmission.engage(&initial).unwrap();
    let after = result.impulse.state_after;
    assert_eq!(initial.rotation, after.rotation);
    assert_eq!(initial.total_momentum_body, after.total_momentum_body);
    assert_eq!(initial.relative_angle_rad, after.relative_angle_rad);
    near(result.constraint_impulse_nms[0], lambda, 2e-14);
    near(
        result.dissipated_energy_j,
        0.5 * mobility * lambda * lambda,
        3e-14,
    );
    near(
        bank.kinetic_energy(&initial).unwrap() - bank.kinetic_energy(&after).unwrap(),
        result.dissipated_energy_j,
        5e-14,
    );
    for free in [[1., 2., 0.], [1., 0., 2.]] {
        let change = (0..3)
            .map(|i| free[i] * (after.axial_momentum_nms[i] - initial.axial_momentum_nms[i]))
            .sum::<f64>();
        near(change, 0., 2e-14);
    }
    let rate = bank.relative_rates(&after).unwrap();
    near(2. * rate[0] - rate[1] - rate[2], 0., 2e-14);
    assert_eq!(
        transmission.engage(&after).unwrap().impulse.state_after,
        after
    );
}

#[test]
fn overlapping_rows_and_author_units_match_independent_kkt() {
    let bank = bank4();
    let rows = [row([2., -1., -1., 0.]), row([0., 3., -1., -2.])];
    let transmission = bank.transmission(rows).unwrap();
    let omega = v(0.3, -0.2, 0.6);
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, 0.3, -0.1)),
            omega,
            [0., 1., -1., 2.],
            [0.4, 0.1, -0.2, 0.7],
        )
        .unwrap();
    let drives = [0.2, -0.4, 0.7, 0.5];
    let external = v(0.3, -0.2, 0.1);
    let expected = original_delta(&bank, &rows, drives, external);
    let result = transmission
        .apply_impulse(&initial, drives, external)
        .unwrap();
    vector_near(
        bank.angular_velocity(&result.impulse.state_after).unwrap(),
        omega + v(expected[0], expected[1], expected[2]),
        2e-13,
    );
    let rates = bank.relative_rates(&result.impulse.state_after).unwrap();
    for i in 0..4 {
        near(rates[i], [0., 1., -1., 2.][i] + expected[3 + i], 2e-13);
    }
    for e in 0..2 {
        near(result.constraint_impulse_nms[e], expected[7 + e], 2e-13);
    }
    near(result.reaction_work_j, 0., 3e-13);
    near(result.work_energy_residual_j, 0., 3e-13);

    // Reverse row order, reverse one row orientation, and rescale by 360 orders
    // of magnitude. Motion is invariant; authored multiplier units change.
    let scaled = bank
        .transmission([
            row(rows[1].coefficients.map(|x| x * 1e-180)),
            row(rows[0].coefficients.map(|x| x * (-1e180))),
        ])
        .unwrap();
    let other = scaled.apply_impulse(&initial, drives, external).unwrap();
    state_near(
        &other.impulse.state_after,
        &result.impulse.state_after,
        3e-13,
    );
    near(
        other.constraint_impulse_nms[0] * 1e-180,
        result.constraint_impulse_nms[1],
        2e-13,
    );
    near(
        other.constraint_impulse_nms[1] * (-1e180),
        result.constraint_impulse_nms[0],
        2e-13,
    );
    for (a, b) in scaled
        .retained_inertia()
        .matrix()
        .m
        .into_iter()
        .zip(transmission.retained_inertia().matrix().m)
    {
        near(a, b, 3e-13);
    }
}

// Continuous original-coordinate Euler equations plus quaternion kinematics.
// The only transmission information is the authored row in the original KKT
// system; no production projection, retained inertia or phase map is reused.
fn oracle_rhs<const N: usize, const M: usize>(
    y: &[f64],
    bank: &RotorBank<N>,
    matrix: &[Vec<f64>],
    drives: [f64; N],
    external: Vec3,
) -> Vec<f64> {
    let omega = v(y[0], y[1], y[2]);
    let mut total = bank.reduced_inertia().momentum(omega).unwrap();
    for i in 0..N {
        let r = bank.rotors()[i];
        let a = r.axis_body();
        total = total + a.scale(r.axial_inertia_kg_m2() * (a.dot(omega) + y[3 + i]));
    }
    let mut rhs = vec![0.; 3 + N + M];
    rhs[..3].copy_from_slice(&xyz(total.cross(omega) + external));
    rhs[3..3 + N].copy_from_slice(&drives);
    let acceleration = solve(matrix.to_vec(), rhs);
    let mut dy = vec![0.; y.len()];
    dy[..3 + N].copy_from_slice(&acceleration[..3 + N]);
    let q = 3 + N;
    let (w, x, y_axis, z) = (y[q], y[q + 1], y[q + 2], y[q + 3]);
    dy[q] = -0.5 * (x * omega.x + y_axis * omega.y + z * omega.z);
    dy[q + 1] = 0.5 * (w * omega.x + y_axis * omega.z - z * omega.y);
    dy[q + 2] = 0.5 * (w * omega.y + z * omega.x - x * omega.z);
    dy[q + 3] = 0.5 * (w * omega.z + x * omega.y - y_axis * omega.x);
    for i in 0..N {
        dy[q + 4 + i] = y[3 + i];
    }
    for e in 0..M {
        dy[q + 4 + N + e] = acceleration[3 + N + e];
    }
    dy
}
fn oracle<const N: usize, const M: usize>(
    initial: &RotorState<N>,
    bank: &RotorBank<N>,
    rows: &[ShaftConstraint<N>; M],
    drives: [f64; N],
    external: Vec3,
    duration: f64,
    steps: usize,
) -> Vec<f64> {
    let q = initial.rotation.as_quat();
    let mut y = Vec::new();
    y.extend(xyz(bank.angular_velocity(initial).unwrap()));
    y.extend(bank.relative_rates(initial).unwrap());
    y.extend([q.w, q.x, q.y, q.z]);
    y.extend(initial.relative_angle_rad);
    y.extend([0.; M]);
    let matrix = original_kkt(bank, rows);
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let a = oracle_rhs::<N, M>(&y, bank, &matrix, drives, external);
        let stage = |d: &[f64], fraction: f64| {
            y.iter()
                .zip(d)
                .map(|(y, d)| y + fraction * dt * d)
                .collect::<Vec<_>>()
        };
        let b = oracle_rhs::<N, M>(&stage(&a, 0.5), bank, &matrix, drives, external);
        let c = oracle_rhs::<N, M>(&stage(&b, 0.5), bank, &matrix, drives, external);
        let d = oracle_rhs::<N, M>(&stage(&c, 1.), bank, &matrix, drives, external);
        for i in 0..y.len() {
            y[i] += dt * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.;
        }
    }
    y
}

#[test]
fn precessing_differential_converges_to_original_euler_equations() {
    let bank = bank3(false);
    let rows = [row([2., -1., -1.])];
    let transmission = bank.transmission(rows).unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [3., 5., 1.],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let drives = [0.4, -0.3, 0.2];
    let external = v(0.1, -0.2, 0.3);
    let duration = 0.3;
    let reference = oracle(&initial, &bank, &rows, drives, external, duration, 2000);
    let expected_rotation = So3::try_from_quat(Quat {
        w: reference[6],
        x: reference[7],
        y: reference[8],
        z: reference[9],
    })
    .unwrap();
    let mut errors = Vec::new();
    for steps in [30, 60] {
        let mut state = initial;
        let mut reaction = 0.;
        for _ in 0..steps {
            let result = transmission
                .step(
                    &state,
                    drives,
                    HeldTorque::Body(external),
                    duration / steps as f64,
                )
                .unwrap();
            near(result.reaction_work_j, 0., 5e-13);
            near(result.work_energy_residual_j, 0., 5e-12);
            vector_near(
                result.motion.angular_momentum_residual_world,
                v(0., 0., 0.),
                1e-12,
            );
            reaction += result.constraint_impulse_nms[0];
            state = result.motion.state_after;
        }
        let omega = bank.angular_velocity(&state).unwrap();
        let mut error = (omega - v(reference[0], reference[1], reference[2])).norm();
        error = error.max(rotation_error(state.rotation, expected_rotation));
        for i in 0..3 {
            error = error.max((bank.relative_rates(&state).unwrap()[i] - reference[3 + i]).abs());
            error = error.max((state.relative_angle_rad[i] - reference[10 + i]).abs());
        }
        error = error.max((reaction - reference[13]).abs());
        errors.push(error);
    }
    assert!(errors[0] < 1e-4, "errors {errors:?}");
    assert!(
        errors[0] > 3.5 * errors[1],
        "second-order errors {errors:?}"
    );
}

#[test]
fn free_precession_retains_world_momentum_and_each_differential_free_momentum() {
    let bank = bank3(false);
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [3., 5., 1.],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let momentum = initial
        .rotation
        .rotate(initial.total_momentum_body)
        .unwrap();
    let mut state = initial;
    let mut defect = 0.;
    for _ in 0..100 {
        let result = transmission.drift(&state, 0.003).unwrap();
        vector_near(result.angular_momentum_residual_world, v(0., 0., 0.), 1e-12);
        near(result.work_energy_residual_j, 0., 2e-12);
        defect += result.drift_energy_defect_j;
        state = result.state_after;
    }
    vector_near(
        state.rotation.rotate(state.total_momentum_body).unwrap(),
        momentum,
        3e-12,
    );
    for free in [[1., 2., 0.], [1., 0., 2.]] {
        near(
            (0..3)
                .map(|i| free[i] * (state.axial_momentum_nms[i] - initial.axial_momentum_nms[i]))
                .sum(),
            0.,
            1e-12,
        );
    }
    let delta_phase: [f64; 3] =
        core::array::from_fn(|i| state.relative_angle_rad[i] - initial.relative_angle_rad[i]);
    near(
        2. * delta_phase[0] - delta_phase[1] - delta_phase[2],
        0.,
        2e-12,
    );
    assert!((state.axial_momentum_nms[0] - initial.axial_momentum_nms[0]).abs() > 1e-3);
    near(
        bank.kinetic_energy(&state).unwrap() - bank.kinetic_energy(&initial).unwrap(),
        defect,
        3e-12,
    );

    let driven = transmission
        .step(
            &initial,
            [0.3, -0.2, 0.1],
            HeldTorque::Space(v(0.1, 0.2, -0.1)),
            0.02,
        )
        .unwrap();
    let reversed = transmission
        .step(
            &driven.motion.state_after,
            [0.3, -0.2, 0.1],
            HeldTorque::Space(v(0.1, 0.2, -0.1)),
            -0.02,
        )
        .unwrap();
    state_near(&reversed.motion.state_after, &initial, 2e-12);
}

#[test]
fn unconstrained_pairwise_and_fully_locked_limits_keep_native_mechanics() {
    let bank = bank3(false);
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [2., -1., 0.4],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let free = bank.transmission([]).unwrap();
    let drive = [0.2, 0.3, -0.1];
    let external = HeldTorque::Body(v(0.1, -0.1, 0.2));
    assert_eq!(
        free.step(&initial, drive, external, 0.03).unwrap().motion,
        bank.step(&initial, drive, external, 0.03).unwrap()
    );
    let gears = bank
        .gear_train([GearConstraint {
            first_rotor: 0,
            second_rotor: 1,
            ratio: -2.,
        }])
        .unwrap();
    let relation = bank.transmission([row([1., 2., 0.])]).unwrap();
    let geared = gears.step(&initial, drive, external, 0.03).unwrap();
    let general = relation.step(&initial, drive, external, 0.03).unwrap();
    state_near(
        &geared.motion.state_after,
        &general.motion.state_after,
        2e-12,
    );
    near(
        geared.constraint_impulse_nms[0],
        general.constraint_impulse_nms[0],
        1e-12,
    );

    let locked = bank.lock_mode([true; 3]).unwrap();
    let all = bank
        .transmission([row([1., 0., 0.]), row([0., 1., 0.]), row([0., 0., 1.])])
        .unwrap();
    let expected = locked.engage(&initial).unwrap();
    let actual = all.engage(&initial).unwrap();
    state_near(
        &expected.impulse.state_after,
        &actual.impulse.state_after,
        2e-12,
    );
    near(
        expected.dissipated_energy_j,
        actual.dissipated_energy_j,
        3e-12,
    );
    let expected = locked.drift(&expected.impulse.state_after, 0.02).unwrap();
    let actual = all.drift(&actual.impulse.state_after, 0.02).unwrap();
    state_near(&expected.state_after, &actual.state_after, 3e-12);
    assert_eq!(
        actual.state_after.relative_angle_rad,
        initial.relative_angle_rad
    );
}

#[test]
fn rank_state_and_resolution_refusals_are_atomic() {
    let bank = bank3(false);
    for rows in [
        [row([2., -1., -1.]), row([4., -2., -2.])],
        [row([2., -1., -1.]), row([2., -1., -1. + 1e-15])],
        [row([0.; 3]), row([1., 0., 0.])],
        [row([f64::NAN, 0., 0.]), row([1., 0., 0.])],
    ] {
        assert!(bank.transmission(rows).is_err());
    }
    assert!(bank.transmission([row([1., 0., 0.]); 4]).is_err());
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            v(0.5, -0.4, 0.7),
            [3., 5., 1.],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let original = initial;
    assert!(
        transmission
            .step(&initial, [0.; 3], HeldTorque::Body(Vec3::default()), 0.)
            .is_err()
    );
    assert!(
        transmission
            .apply_impulse(&initial, [f64::NAN, 0., 0.], Vec3::default())
            .is_err()
    );
    assert!(
        transmission
            .apply_impulse(&initial, [1e-100, 0., 0.], Vec3::default())
            .is_err()
    );
    assert_eq!(initial, original);
}

#[test]
fn full_rank_drive_reactions_have_no_free_work_and_stationary_carrier_stays_at_rest() {
    let bank = bank3(false);
    let all = bank
        .transmission([row([1., 0., 0.]), row([0., 1., 0.]), row([0., 0., 1.])])
        .unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [0.; 3],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let drives = [0.2, -0.3, 0.1];
    let held = all
        .apply_impulse(&initial, drives, Vec3::default())
        .unwrap();
    assert_eq!(held.impulse.state_after, initial);
    assert_eq!(held.drive_work_j, [0.; 3]);
    assert_eq!(held.reaction_work_j, 0.);
    for i in 0..3 {
        near(held.constraint_impulse_nms[i], -drives[i], 3e-14);
    }
    let external = v(0.2, -0.1, 0.3);
    let result = all.apply_impulse(&initial, drives, external).unwrap();
    let expected = all
        .retained_inertia()
        .velocity(initial.total_momentum_body + external)
        .unwrap();
    vector_near(
        bank.angular_velocity(&result.impulse.state_after).unwrap(),
        expected,
        2e-13,
    );
    assert_eq!(result.drive_work_j, [0.; 3]);
    near(
        result.impulse.kinetic_energy_change_j,
        result.impulse.external_work_j,
        2e-13,
    );

    let differential = bank.transmission([row([2., -1., -1.])]).unwrap();
    let rest = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            Vec3::default(),
            [3., 5., 1.],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let drift = differential.drift(&rest, 0.1).unwrap();
    assert_eq!(drift.state_after.rotation, rest.rotation);
    assert_eq!(
        drift.state_after.total_momentum_body,
        rest.total_momentum_body
    );
    assert_eq!(
        drift.state_after.axial_momentum_nms,
        rest.axial_momentum_nms
    );
    near(
        drift.state_after.relative_angle_rad[1] - rest.relative_angle_rad[1],
        0.5,
        2e-14,
    );
}

#[test]
fn nearly_parallel_rows_cannot_hide_a_forbidden_physical_rate() {
    let bank = bank3(false);
    let transmission = bank
        .transmission([row([1., 1., 0.]), row([1., 1. + 1e-10, 0.])])
        .unwrap();
    // These rows imply s0=s1=0. Their tiny algebraic residual alone is within
    // the absolute-spin rounding scale; orthonormal physical modes still see
    // the finite forbidden differential speed and require explicit capture.
    let initial = bank
        .state_from_velocities(So3::identity(), v(1e3, 1e3, 1e3), [1., -1., 2.], [0.; 3])
        .unwrap();
    assert!(
        transmission
            .apply_impulse(&initial, [0.; 3], Vec3::default())
            .is_err()
    );
    assert!(transmission.drift(&initial, 0.01).is_err());
}

#[test]
fn fixed_ring_planetary_relation_keeps_the_locked_original_phase() {
    let bank = bank3(true);
    let rows = [row([2., 5., -7.]), row([0., 1., 0.])];
    let transmission = bank.transmission(rows).unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            v(0., 0., 0.4),
            [7., 0., 2.],
            [0.7, -0.4, 0.2],
        )
        .unwrap();
    let drives = [1.2, -0.4, 0.7];
    let external = v(0., 0., 0.3);
    let expected = original_delta(&bank, &rows, drives, external);
    let dt = 0.4;
    let result = transmission
        .step(&initial, drives, HeldTorque::Body(external), dt)
        .unwrap();
    let after = result.motion.state_after;
    assert_eq!(after.relative_angle_rad[1], initial.relative_angle_rad[1]);
    assert_eq!(result.drive_work_j[1], 0.0);
    for i in 0..3 {
        near(
            bank.relative_rates(&after).unwrap()[i],
            [7., 0., 2.][i] + dt * expected[3 + i],
            3e-13,
        );
        near(
            after.relative_angle_rad[i] - initial.relative_angle_rad[i],
            [7., 0., 2.][i] * dt + 0.5 * expected[3 + i] * dt * dt,
            3e-13,
        );
    }
    for e in 0..2 {
        near(
            result.constraint_impulse_nms[e],
            dt * expected[6 + e],
            3e-13,
        );
    }
    near(result.work_energy_residual_j, 0., 3e-13);
}
