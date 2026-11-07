//! Original-coordinate mechanics for passive loads inside ideal transmissions.

use fs_ga::{Mat3, Quat, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, RotorBank, RotorState, ShaftConstraint, TorsionalCoupling,
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
fn load(first: usize, second: usize, ratio: f64, k: f64, d: f64, rest: f64) -> TorsionalCoupling {
    TorsionalCoupling {
        first_rotor: first,
        second_rotor: second,
        ratio,
        stiffness_nm_per_rad: k,
        damping_nms_per_rad: d,
        rest_angle_rad: rest,
    }
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
        .map(|axis| (a.rotate(axis).unwrap() - b.rotate(axis).unwrap()).norm())
        .fold(0., f64::max)
}
fn state_near<const N: usize>(a: &RotorState<N>, b: &RotorState<N>, tol: f64) {
    vector_near(a.total_momentum_body, b.total_momentum_body, tol);
    assert!(rotation_error(a.rotation, b.rotation) < tol);
    for i in 0..N {
        near(a.axial_momentum_nms[i], b.axial_momentum_nms[i], tol);
        near(a.relative_angle_rad[i], b.relative_angle_rad[i], tol);
    }
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
fn port<const N: usize>(c: TorsionalCoupling) -> [f64; N] {
    let mut b = [0.; N];
    b[c.first_rotor] = 1.;
    b[c.second_rotor] = -c.ratio;
    b
}
fn mapped(c: TorsionalCoupling, x: &[f64]) -> f64 {
    x[c.first_rotor] - c.ratio * x[c.second_rotor]
}

// Independent original carrier/shaft kinetic matrix, augmented by the authored
// ideal constraints. This does not use the production projection or retained
// inertia. Partial-pivot elimination is confined to the small test oracle.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        assert!(a[p][k].abs() > 1e-14);
        a.swap(p, k);
        b.swap(p, k);
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
    let mut matrix = vec![vec![0.; 3 + N + M]; 3 + N + M];
    let inertia = bank.reduced_inertia().matrix().m;
    for r in 0..3 {
        for c in 0..3 {
            matrix[r][c] = inertia[3 * r + c];
        }
    }
    for (i, rotor) in bank.rotors().iter().enumerate() {
        let axis = xyz(rotor.axis_body());
        let j = rotor.axial_inertia_kg_m2();
        for r in 0..3 {
            for c in 0..3 {
                matrix[r][c] += j * axis[r] * axis[c];
            }
            matrix[r][3 + i] = j * axis[r];
            matrix[3 + i][r] = j * axis[r];
        }
        matrix[3 + i][3 + i] = j;
        for e in 0..M {
            matrix[3 + i][3 + N + e] = -rows[e].coefficients[i];
            matrix[3 + N + e][3 + i] = rows[e].coefficients[i];
        }
    }
    matrix
}
struct ReferenceStage<const N: usize, const M: usize, const E: usize> {
    omega: Vec3,
    rates: [f64; N],
    phase: [f64; N],
    constraint: [f64; M],
    coupling: [f64; E],
    heat: [f64; E],
}
#[allow(clippy::too_many_arguments)]
fn original_midpoint<const N: usize, const M: usize, const E: usize>(
    bank: &RotorBank<N>,
    rows: &[ShaftConstraint<N>; M],
    couplings: &[TorsionalCoupling; E],
    omega: Vec3,
    rates: [f64; N],
    phase: [f64; N],
    drive: [f64; N],
    external: Vec3,
    delta: f64,
) -> ReferenceStage<N, M, E> {
    let mut matrix = original_kkt(bank, rows);
    let mut rhs = vec![0.; 3 + N + M];
    rhs[..3].copy_from_slice(&xyz(external.scale(delta)));
    for i in 0..N {
        rhs[3 + i] = delta * drive[i];
    }
    for &c in couplings {
        let b = port::<N>(c);
        let effective = c.damping_nms_per_rad + 0.5 * delta * c.stiffness_nm_per_rad;
        let old_effort = c.stiffness_nm_per_rad * (mapped(c, &phase) - c.rest_angle_rad)
            + effective * mapped(c, &rates);
        for i in 0..N {
            rhs[3 + i] -= delta * b[i] * old_effort;
            for j in 0..N {
                matrix[3 + i][3 + j] += 0.5 * delta * b[i] * effective * b[j];
            }
        }
    }
    let change = solve(matrix, rhs);
    let mean: [f64; N] = core::array::from_fn(|i| rates[i] + 0.5 * change[3 + i]);
    let after = core::array::from_fn(|i| phase[i] + delta * mean[i]);
    ReferenceStage {
        omega: omega + v(change[0], change[1], change[2]),
        rates: core::array::from_fn(|i| rates[i] + change[3 + i]),
        phase: after,
        constraint: core::array::from_fn(|e| change[3 + N + e]),
        coupling: core::array::from_fn(|e| {
            let c = couplings[e];
            delta
                * (0.5
                    * c.stiffness_nm_per_rad
                    * (mapped(c, &phase) + mapped(c, &after) - 2. * c.rest_angle_rad)
                    + c.damping_nms_per_rad * mapped(c, &mean))
        }),
        heat: core::array::from_fn(|e| {
            let c = couplings[e];
            delta * c.damping_nms_per_rad * mapped(c, &mean).powi(2)
        }),
    }
}

#[test]
fn unequal_differential_loads_match_original_midpoint_kinetics_and_work() {
    let bank = bank3(true);
    let rows = [row([2., -1., -1.])];
    let couplings = [
        load(1, 2, 1., 2.5, 0.15, 0.3),
        load(1, 0, 0., 0., 0.7, 0.),
        load(2, 0, 0., 0., 1.1, 0.),
    ];
    let transmission = bank.transmission(rows).unwrap();
    let dt = 0.08;
    let loaded = transmission
        .torsional_network_stepper(couplings, dt)
        .unwrap();
    let omega = v(0., 0., 0.4);
    let rates = [3., 5., 1.];
    let phase = [0.7, -0.4, 0.2];
    let initial = bank
        .state_from_velocities(So3::identity(), omega, rates, phase)
        .unwrap();
    let drive = [1.2, -0.1, 0.2];
    let external = v(0., 0., 0.3);
    let first = original_midpoint(
        &bank,
        &rows,
        &couplings,
        omega,
        rates,
        phase,
        drive,
        external,
        0.5 * dt,
    );
    let second = original_midpoint(
        &bank,
        &rows,
        &couplings,
        first.omega,
        first.rates,
        first.phase,
        drive,
        external,
        0.5 * dt,
    );
    let result = loaded
        .step(&initial, drive, HeldTorque::Space(external))
        .unwrap();
    let after = result.motion.state_after;
    vector_near(bank.angular_velocity(&after).unwrap(), second.omega, 2e-13);
    assert!(rotation_error(after.rotation, turn(first.omega.scale(dt))) < 2e-13);
    for i in 0..3 {
        near(
            bank.relative_rates(&after).unwrap()[i],
            second.rates[i],
            2e-13,
        );
        near(after.relative_angle_rad[i], second.phase[i], 2e-13);
        near(
            result.drive_work_j[i],
            drive[i] * (after.relative_angle_rad[i] - phase[i]),
            2e-13,
        );
        let load_impulse = (0..3)
            .map(|e| port::<3>(couplings[e])[i] * result.coupling_impulse_nms[e])
            .sum::<f64>();
        near(
            result.motion.axial_momentum_change_nms[i],
            dt * drive[i] - load_impulse
                + rows[0].coefficients[i] * result.constraint_impulse_nms[0],
            3e-13,
        );
    }
    near(
        result.constraint_impulse_nms[0],
        first.constraint[0] + second.constraint[0],
        2e-13,
    );
    for e in 0..3 {
        near(
            result.coupling_impulse_nms[e],
            first.coupling[e] + second.coupling[e],
            2e-13,
        );
        near(
            result.dissipated_energy_j[e],
            first.heat[e] + second.heat[e],
            2e-13,
        );
        assert!(result.dissipated_energy_j[e] >= 0.);
    }
    near(
        loaded.potential_energy(&after).unwrap() - loaded.potential_energy(&initial).unwrap(),
        result.elastic_energy_change_j.iter().sum(),
        3e-13,
    );
    near(result.reaction_work_j, 0., 3e-13);
    near(result.work_energy_residual_j, 0., 3e-13);
    vector_near(
        result.motion.angular_momentum_residual_world,
        Vec3::default(),
        3e-13,
    );
    assert!((second.rates[1] - second.rates[2]).abs() > 3.);
}

#[test]
fn viscous_outputs_transfer_momentum_to_the_finite_carrier_and_heat_once() {
    let bank = bank3(true);
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let loaded = transmission
        .torsional_network_stepper(
            [load(1, 0, 0., 0., 0.7, 0.), load(2, 0, 0., 0., 1.1, 0.)],
            0.05,
        )
        .unwrap();
    let initial = bank
        .state_from_velocities(So3::identity(), v(0., 0., 0.4), [3., 5., 1.], [0.; 3])
        .unwrap();
    let initial_energy = bank.kinetic_energy(&initial).unwrap();
    let momentum = initial.total_momentum_body.z;
    let limit = momentum / (7. + 0.7 + 1.1 + 1.5);
    let mut state = initial;
    let mut heat = 0.;
    let mut defect = 0.;
    for iteration in 0..100 {
        let result = loaded
            .step(&state, [0.; 3], HeldTorque::Body(Vec3::default()))
            .unwrap_or_else(|error| {
                panic!(
                    "damping iteration {iteration}, rates {:?}: {error:?}",
                    bank.relative_rates(&state).unwrap()
                )
            });
        for q in result.dissipated_energy_j {
            assert!(q >= 0.);
            heat += q;
        }
        defect += result.motion.drift_energy_defect_j;
        near(result.work_energy_residual_j, 0., 3e-12);
        state = result.motion.state_after;
    }
    let rates = bank.relative_rates(&state).unwrap();
    assert!(
        rates.into_iter().map(f64::abs).fold(0., f64::max) < 0.2,
        "{rates:?}"
    );
    // The right output crosses zero during the admitted interval; phase
    // addition roundoff must not masquerade as a midpoint-law error.
    assert!(rates[2] < 0.0);
    let omega = bank.angular_velocity(&state).unwrap().z;
    assert!(omega > 1.2 && (omega - limit).abs() < 0.05);
    near(state.total_momentum_body.z, momentum, 3e-12);
    near(
        bank.kinetic_energy(&state).unwrap() + heat - initial_energy - defect,
        0.,
        2e-11,
    );
    assert!(heat > 10.);
}

#[test]
fn conservative_differential_spring_has_the_original_frequency_and_no_heat() {
    let bank = bank3(true);
    let rows = [row([2., -1., -1.])];
    let transmission = bank.transmission(rows).unwrap();
    let coupling = load(1, 2, 1., 5., 0., 0.1);
    let dt = 0.03;
    let loaded = transmission
        .torsional_network_stepper([coupling], dt)
        .unwrap();
    let mut rhs = vec![0.; 7];
    rhs[4] = 1.;
    rhs[5] = -1.;
    let response = solve(original_kkt(&bank, &rows), rhs);
    let frequency = (5. * (response[4] - response[5])).sqrt();
    near(loaded.angular_frequency_bound_rad_s(), frequency, 2e-13);
    let initial = bank
        .state_from_velocities(So3::identity(), v(0., 0., 0.4), [0.; 3], [0., 0.4, -0.2])
        .unwrap();
    let energy =
        bank.kinetic_energy(&initial).unwrap() + loaded.potential_energy(&initial).unwrap();
    let strain0 = mapped(coupling, &initial.relative_angle_rad) - coupling.rest_angle_rad;
    let phase_step = 4. * (0.25 * dt * frequency).atan();
    let mut state = initial;
    for n in 1..=80 {
        let result = loaded
            .step(&state, [0.; 3], HeldTorque::Body(Vec3::default()))
            .unwrap();
        assert_eq!(result.dissipated_energy_j, [0.]);
        near(result.work_energy_residual_j, 0., 3e-12);
        state = result.motion.state_after;
        near(
            mapped(coupling, &state.relative_angle_rad) - coupling.rest_angle_rad,
            strain0 * (n as f64 * phase_step).cos(),
            3e-11,
        );
        near(
            mapped(coupling, &bank.relative_rates(&state).unwrap()),
            -frequency * strain0 * (n as f64 * phase_step).sin(),
            3e-11,
        );
    }
    near(
        bank.kinetic_energy(&state).unwrap() + loaded.potential_energy(&state).unwrap(),
        energy,
        3e-11,
    );
    near(bank.angular_velocity(&state).unwrap().z, 0.4, 2e-12);
}

#[test]
fn dependent_parallel_and_reversed_load_rows_keep_their_physical_shares() {
    let bank = bank3(true);
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let dt = 0.04;
    let single = transmission
        .torsional_network_stepper([load(1, 2, 1., 6., 0.9, 0.2)], dt)
        .unwrap();
    let split = transmission
        .torsional_network_stepper(
            [load(1, 2, 1., 2., 0.3, 0.2), load(2, 1, 1., 4., 0.6, -0.2)],
            dt,
        )
        .unwrap();
    // Forty dependent ports also exercise the native blocked Cholesky/GEMM
    // path past its 32-column panel boundary.
    let many = transmission
        .torsional_network_stepper([load(1, 2, 1., 6. / 40., 0.9 / 40., 0.2); 40], dt)
        .unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            v(0., 0., 0.4),
            [3., 5., 1.],
            [0.7, -0.4, 0.2],
        )
        .unwrap();
    let drive = [0.8, 0.2, -0.1];
    let external = HeldTorque::Body(v(0., 0., 0.1));
    let a = single.step(&initial, drive, external).unwrap();
    let b = split.step(&initial, drive, external).unwrap();
    let c = many.step(&initial, drive, external).unwrap();
    state_near(&a.motion.state_after, &b.motion.state_after, 3e-12);
    state_near(&a.motion.state_after, &c.motion.state_after, 3e-12);
    near(
        b.coupling_impulse_nms[0],
        a.coupling_impulse_nms[0] / 3.,
        2e-12,
    );
    near(
        b.coupling_impulse_nms[1],
        -2. * a.coupling_impulse_nms[0] / 3.,
        2e-12,
    );
    near(
        b.dissipated_energy_j.iter().sum(),
        a.dissipated_energy_j[0],
        2e-12,
    );
    near(
        c.dissipated_energy_j.iter().sum(),
        a.dissipated_energy_j[0],
        2e-12,
    );
}

#[test]
fn blocked_and_fixed_ring_loads_preserve_phases_and_report_static_reactions() {
    let bank = bank3(false);
    let transmission = bank
        .transmission([row([1., 0., 0.]), row([0., 1., 0.]), row([0., 0., 1.])])
        .unwrap();
    let dt = 0.04;
    let couplings = [load(0, 1, -1., 3., 0.2, -0.1), load(2, 0, 0., 1., 0.1, 0.4)];
    let loaded = transmission
        .torsional_network_stepper(couplings, dt)
        .unwrap();
    near(loaded.half_stage_resolution_number(), 0., 0.);
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [0.; 3],
            [0.3, -0.2, 0.1],
        )
        .unwrap();
    let drive = [0.2, -0.3, 0.1];
    let result = loaded
        .step(&initial, drive, HeldTorque::Body(Vec3::default()))
        .unwrap();
    let expected = transmission.drift(&initial, dt).unwrap();
    state_near(&result.motion.state_after, &expected.state_after, 3e-12);
    assert_eq!(
        result.motion.state_after.relative_angle_rad,
        initial.relative_angle_rad
    );
    assert_eq!(result.drive_work_j, [0.; 3]);
    assert_eq!(result.dissipated_energy_j, [0.; 2]);
    assert_eq!(result.elastic_energy_change_j, [0.; 2]);
    for e in 0..2 {
        near(
            result.coupling_impulse_nms[e],
            dt * couplings[e].stiffness_nm_per_rad
                * (mapped(couplings[e], &initial.relative_angle_rad) - couplings[e].rest_angle_rad),
            2e-13,
        );
    }
    for i in 0..3 {
        let passive = (0..2)
            .map(|e| port::<3>(couplings[e])[i] * result.coupling_impulse_nms[e])
            .sum::<f64>();
        near(
            result.constraint_impulse_nms[i],
            result.motion.axial_momentum_change_nms[i] - dt * drive[i] + passive,
            2e-12,
        );
    }

    let bank = bank3(true);
    let rows = [row([2., 5., -7.]), row([0., 1., 0.])];
    let transmission = bank.transmission(rows).unwrap();
    let couplings = [load(1, 0, 0., 2.5, 0.3, 0.1), load(2, 0, 0., 0., 0.6, 0.)];
    let loaded = transmission
        .torsional_network_stepper(couplings, dt)
        .unwrap();
    let omega = v(0., 0., 0.4);
    let rates = [7., 0., 2.];
    let phase = [0.7, -0.4, 0.2];
    let initial = bank
        .state_from_velocities(So3::identity(), omega, rates, phase)
        .unwrap();
    let drive = [1.2, 0.4, 0.];
    let external = v(0., 0., 0.3);
    let first = original_midpoint(
        &bank,
        &rows,
        &couplings,
        omega,
        rates,
        phase,
        drive,
        external,
        0.5 * dt,
    );
    let second = original_midpoint(
        &bank,
        &rows,
        &couplings,
        first.omega,
        first.rates,
        first.phase,
        drive,
        external,
        0.5 * dt,
    );
    let result = loaded
        .step(&initial, drive, HeldTorque::Body(external))
        .unwrap();
    assert_eq!(result.motion.state_after.relative_angle_rad[1], phase[1]);
    assert_eq!(result.drive_work_j[1], 0.);
    assert_eq!(result.dissipated_energy_j[0], 0.);
    vector_near(
        bank.angular_velocity(&result.motion.state_after).unwrap(),
        second.omega,
        3e-12,
    );
    for i in 0..3 {
        near(
            bank.relative_rates(&result.motion.state_after).unwrap()[i],
            second.rates[i],
            3e-12,
        );
    }
    for e in 0..2 {
        near(
            result.constraint_impulse_nms[e],
            first.constraint[e] + second.constraint[e],
            3e-12,
        );
    }
    near(result.work_energy_residual_j, 0., 3e-12);
}

// Continuous original Euler/KKT equations with actual spring/damper forces.
// Quaternion, original shaft phases, physical heats/work and ideal row impulses
// are integrated independently by a refined RK4 reference.
fn rhs<const N: usize, const M: usize, const E: usize>(
    y: &[f64],
    bank: &RotorBank<N>,
    matrix: &[Vec<f64>],
    couplings: &[TorsionalCoupling; E],
    drive: [f64; N],
    external: Vec3,
) -> Vec<f64> {
    let omega = v(y[0], y[1], y[2]);
    let q = 3 + N;
    let phase = q + 4;
    let heat = phase + N;
    let work = heat + E;
    let mut momentum = bank.reduced_inertia().momentum(omega).unwrap();
    for i in 0..N {
        let rotor = bank.rotors()[i];
        let axis = rotor.axis_body();
        momentum =
            momentum + axis.scale(rotor.axial_inertia_kg_m2() * (axis.dot(omega) + y[3 + i]));
    }
    let mut forcing = vec![0.; 3 + N + M];
    forcing[..3].copy_from_slice(&xyz(momentum.cross(omega) + external));
    for i in 0..N {
        forcing[3 + i] = drive[i];
    }
    let mut dy = vec![0.; y.len()];
    for (e, &c) in couplings.iter().enumerate() {
        let speed = mapped(c, &y[3..3 + N]);
        let effort = c.stiffness_nm_per_rad * (mapped(c, &y[phase..phase + N]) - c.rest_angle_rad)
            + c.damping_nms_per_rad * speed;
        for (i, b) in port::<N>(c).into_iter().enumerate() {
            forcing[3 + i] -= b * effort;
        }
        dy[heat + e] = c.damping_nms_per_rad * speed * speed;
    }
    let acceleration = solve(matrix.to_vec(), forcing);
    dy[..3 + N].copy_from_slice(&acceleration[..3 + N]);
    let (w, x, z_y, z) = (y[q], y[q + 1], y[q + 2], y[q + 3]);
    dy[q] = -0.5 * (x * omega.x + z_y * omega.y + z * omega.z);
    dy[q + 1] = 0.5 * (w * omega.x + z_y * omega.z - z * omega.y);
    dy[q + 2] = 0.5 * (w * omega.y + z * omega.x - x * omega.z);
    dy[q + 3] = 0.5 * (w * omega.z + x * omega.y - z_y * omega.x);
    for i in 0..N {
        dy[phase + i] = y[3 + i];
        dy[work + i] = drive[i] * y[3 + i];
    }
    dy[work + N] = external.dot(omega);
    for e in 0..M {
        dy[work + N + 1 + e] = acceleration[3 + N + e];
    }
    dy
}
#[allow(clippy::too_many_arguments)]
fn oracle<const N: usize, const M: usize, const E: usize>(
    bank: &RotorBank<N>,
    rows: &[ShaftConstraint<N>; M],
    couplings: &[TorsionalCoupling; E],
    initial: &RotorState<N>,
    drive: [f64; N],
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
    y.extend([0.; E]);
    y.extend([0.; N]);
    y.push(0.);
    y.extend([0.; M]);
    let matrix = original_kkt(bank, rows);
    let dt = duration / steps as f64;
    for _ in 0..steps {
        let a = rhs::<N, M, E>(&y, bank, &matrix, couplings, drive, external);
        let stage = |d: &[f64], fraction: f64| {
            y.iter()
                .zip(d)
                .map(|(y, d)| y + fraction * dt * d)
                .collect::<Vec<_>>()
        };
        let b = rhs::<N, M, E>(&stage(&a, 0.5), bank, &matrix, couplings, drive, external);
        let c = rhs::<N, M, E>(&stage(&b, 0.5), bank, &matrix, couplings, drive, external);
        let d = rhs::<N, M, E>(&stage(&c, 1.), bank, &matrix, couplings, drive, external);
        for i in 0..y.len() {
            y[i] += dt * (a[i] + 2. * b[i] + 2. * c[i] + d[i]) / 6.;
        }
    }
    y
}

#[test]
fn loaded_nonparallel_modes_converge_to_original_euler_kinetics_and_heat() {
    let bank = bank4();
    let rows = [row([2., -1., -1., 0.]), row([0., 3., -1., -2.])];
    let couplings = [
        load(0, 1, 0.5, 1.2, 0.15, 0.1),
        load(2, 3, -0.4, 0.8, 0.2, -0.2),
        load(1, 3, 0., 0., 0.3, 0.),
    ];
    let transmission = bank.transmission(rows).unwrap();
    let initial = bank
        .state_from_velocities(
            turn(v(0.2, -0.3, 0.1)),
            v(0.5, -0.4, 0.7),
            [2.2, 1.3, 3.1, 0.4],
            [0.2, -0.1, 0.3, 0.4],
        )
        .unwrap();
    let drive = [0.2, -0.1, 0.15, 0.05];
    let external = v(0.1, -0.2, 0.3);
    let duration = 0.3;
    let reference = oracle(
        &bank, &rows, &couplings, &initial, drive, external, duration, 2400,
    );
    let expected_rotation = So3::try_from_quat(Quat {
        w: reference[7],
        x: reference[8],
        y: reference[9],
        z: reference[10],
    })
    .unwrap();
    let mut errors = Vec::new();
    for steps in [30, 60] {
        let loaded = transmission
            .torsional_network_stepper(couplings, duration / steps as f64)
            .unwrap();
        let mut state = initial;
        let mut heat = [0.; 3];
        let mut work = [0.; 4];
        let mut ext = 0.;
        let mut reactions = [0.; 2];
        for _ in 0..steps {
            let result = loaded
                .step(&state, drive, HeldTorque::Body(external))
                .unwrap();
            for e in 0..3 {
                heat[e] += result.dissipated_energy_j[e];
            }
            for i in 0..4 {
                work[i] += result.drive_work_j[i];
            }
            for e in 0..2 {
                reactions[e] += result.constraint_impulse_nms[e];
            }
            ext += result.motion.external_work_j;
            near(result.reaction_work_j, 0., 4e-12);
            near(result.work_energy_residual_j, 0., 4e-12);
            vector_near(
                result.motion.angular_momentum_residual_world,
                Vec3::default(),
                2e-12,
            );
            state = result.motion.state_after;
        }
        let omega = bank.angular_velocity(&state).unwrap();
        let mut error = (omega - v(reference[0], reference[1], reference[2])).norm();
        error = error.max(rotation_error(state.rotation, expected_rotation));
        for i in 0..4 {
            error = error.max((bank.relative_rates(&state).unwrap()[i] - reference[3 + i]).abs());
            error = error.max((state.relative_angle_rad[i] - reference[11 + i]).abs());
            error = error.max((work[i] - reference[18 + i]).abs());
        }
        for e in 0..3 {
            error = error.max((heat[e] - reference[15 + e]).abs());
        }
        for e in 0..2 {
            error = error.max((reactions[e] - reference[23 + e]).abs());
        }
        error = error.max((ext - reference[22]).abs());
        errors.push(error);
    }
    assert!(errors[0] < 2e-4, "errors {errors:?}");
    assert!(
        errors[0] > 3.4 * errors[1],
        "second-order errors {errors:?}"
    );
}

#[test]
fn native_limits_and_unresolved_loaded_steps_refuse_atomically() {
    let bank = bank3(true);
    let transmission = bank.transmission([row([2., -1., -1.])]).unwrap();
    let initial = bank
        .state_from_velocities(
            So3::identity(),
            v(0., 0., 0.4),
            [3., 5., 1.],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let drive = [0.2, 0.1, -0.3];
    let torque = HeldTorque::Space(v(0.1, -0.2, 0.3));
    let dt = 0.02;
    let native = transmission.step(&initial, drive, torque, dt).unwrap();
    let empty = transmission
        .torsional_network_stepper([], dt)
        .unwrap()
        .step(&initial, drive, torque)
        .unwrap();
    let zero = transmission
        .torsional_network_stepper([load(1, 2, -0.3, 0., 0., 0.4)], dt)
        .unwrap()
        .step(&initial, drive, torque)
        .unwrap();
    assert_eq!(empty.motion, native.motion);
    assert_eq!(zero.motion, native.motion);
    assert_eq!(empty.constraint_impulse_nms, native.constraint_impulse_nms);
    assert_eq!(zero.drive_work_j, native.drive_work_j);

    let laws = [load(1, 2, -0.3, 1., 0.2, 0.1)];
    let free = bank
        .transmission([])
        .unwrap()
        .torsional_network_stepper(laws, dt)
        .unwrap();
    let original = bank.torsional_network_stepper(laws, dt).unwrap();
    let actual = free.step(&initial, drive, torque).unwrap();
    let expected = original.step(&initial, drive, torque).unwrap();
    state_near(
        &actual.motion.state_after,
        &expected.motion.state_after,
        3e-12,
    );
    near(
        actual.dissipated_energy_j[0],
        expected.dissipated_energy_j[0],
        2e-12,
    );

    for law in [
        load(0, 3, 1., 1., 0., 0.),
        load(1, 1, 1., 1., 0., 0.),
        load(1, 2, 1., -1., 0., 0.),
        load(1, 2, 1., 0., -1., 0.),
        load(1, 2, f64::NAN, 1., 0., 0.),
        load(1, 2, 1., 1e12, 0., 0.),
        load(1, 2, 1., 0., 1e12, 0.),
    ] {
        assert!(transmission.torsional_network_stepper([law], dt).is_err());
    }
    assert!(transmission.torsional_network_stepper(laws, 0.).is_err());
    assert!(
        transmission
            .torsional_network_stepper([laws[0]; 65], dt)
            .is_err()
    );
    let loaded = transmission.torsional_network_stepper(laws, dt).unwrap();
    let original = initial;
    assert!(loaded.step(&initial, [f64::NAN, 0., 0.], torque).is_err());
    let mut off_manifold = initial;
    off_manifold.axial_momentum_nms[0] += 0.1;
    let unchanged = off_manifold;
    assert!(loaded.step(&off_manifold, drive, torque).is_err());
    assert_eq!(off_manifold, unchanged);
    assert_eq!(initial, original);
}

#[test]
fn exact_torque_neutral_rows_do_not_erase_neighboring_nonzero_recoil() {
    let bank = bank3(true);
    let rows = [row([2., -1., -1. + 1e-8])];
    let transmission = bank.transmission(rows).unwrap();
    let initial = bank
        .state_from_velocities(So3::identity(), v(0., 0., 0.4), [0.; 3], [0.; 3])
        .unwrap();
    let drive = [0., -1., 1.];
    let mut rhs = vec![0.; 7];
    rhs[4] = -1.;
    rhs[5] = 1.;
    let expected = solve(original_kkt(&bank, &rows), rhs);
    let result = transmission
        .apply_impulse(&initial, drive, Vec3::default())
        .unwrap();
    let recoil = bank
        .angular_velocity(&result.impulse.state_after)
        .unwrap()
        .z
        - 0.4;
    assert!(recoil != 0.);
    near(recoil, expected[2], 3e-14);

    // At this carrier scale a physically nonzero recoil is erased by ordinary
    // f64 canonical storage. Neither a nearly torque-neutral row nor a nearly
    // balanced authored source may use the exact-zero exception.
    let large = bank
        .state_from_velocities(So3::identity(), v(0., 0., 1e8), [0.; 3], [0.; 3])
        .unwrap();
    let near_row = bank.transmission([row([2., -1., -1. + 1e-12])]).unwrap();
    let error = near_row
        .apply_impulse(&large, drive, Vec3::default())
        .unwrap_err();
    assert!(
        matches!(
            error,
            fs_time::lie::LieStepError::InvalidParameter {
                context: "unresolved carrier recoil against rotor momentum"
            }
        ),
        "{error:?}"
    );
    let exact = bank.transmission([row([2., -1., -1.])]).unwrap();
    let error = exact
        .apply_impulse(&large, [0., -1., 1. + 1e-12], Vec3::default())
        .unwrap_err();
    assert!(
        matches!(
            error,
            fs_time::lie::LieStepError::InvalidParameter {
                context: "unresolved carrier recoil against rotor momentum"
            }
        ),
        "{error:?}"
    );
}
