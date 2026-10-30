//! G0/G1/G3: source gain, signed work, wind-axis covariance and nonlinear motion.
use fs_airscrew::dynamic_inflow::first_harmonic::skew::SkewInflowControls;
use fs_airscrew::dynamic_inflow::first_harmonic::{FirstHarmonicInflow, HarmonicInduction};

const R: f64 = 1.3;
const RHO: f64 = 1.225;
fn state(u: [f64; 3]) -> HarmonicInduction {
    HarmonicInduction {
        mean_mps: u[0],
        first_mps: [u[1], u[2]],
    }
}
fn model() -> FirstHarmonicInflow {
    FirstHarmonicInflow::new(R, RHO).unwrap()
}
fn controls() -> SkewInflowControls {
    SkewInflowControls {
        maximum_iterations: 64,
        absolute_tolerance_mps: 1e-14,
        relative_tolerance: 1e-14,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a:e} != {b:e}, defect={:e}", a - b);
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn near_vec(a: [f64; 3], b: [f64; 3], tol: f64) {
    for i in 0..3 {
        near(a[i], b[i], tol);
    }
}
fn wind_transform(incident: [f64; 3], x: [f64; 3]) -> [f64; 3] {
    let magnitude = incident[1].hypot(incident[2]);
    let d = [-incident[1] / magnitude, -incident[2] / magnitude];
    [x[0], d[0] * x[1] + d[1] * x[2], -d[1] * x[1] + d[0] * x[2]]
}
// Independent source GAIN, before inversion. The production owner builds the
// inverse transport by a different closed-form algebra and then rotates it.
fn source_gain(mean: f64, incident: [f64; 3]) -> [[f64; 3]; 3] {
    let transverse = incident[1].hypot(incident[2]);
    let lambda = incident[0] + mean;
    let speed = (transverse * transverse + lambda * lambda).sqrt();
    let vm = (transverse * transverse + lambda * (lambda + mean)) / speed;
    let chi = transverse.atan2(lambda);
    let a = 15.0 * std::f64::consts::PI / 64.0;
    [
        [0.5 / speed, -a * (0.5 * chi).tan() / vm, 0.0],
        [
            a * (0.5 * chi).tan() / speed,
            4.0 * chi.cos() / (vm * (1.0 + chi.cos())),
            0.0,
        ],
        [0.0, 0.0, 4.0 / (vm * (1.0 + chi.cos()))],
    ]
}

#[test]
fn actual_transport_inverts_the_primary_gain_with_both_off_diagonal_signs() {
    for incident in [[6.0, -4.0, 2.0], [1.0, 7.0, -3.0], [0.0, -20.0, 1.0]] {
        for mean in [0.0, 0.8, 2.0] {
            let u = state([mean, 0.2, -0.1]);
            let transport = model().skew_transport(mean, incident).unwrap();
            let forces = wind_transform(incident, transport.force_n(u));
            let expected = wind_transform(incident, u.modes());
            let gain = source_gain(mean, incident);
            let rho_a = RHO * std::f64::consts::PI * R * R;
            near_vec(gain.map(|row| dot(row, forces) / rho_a), expected, 2e-12);
        }
    }
    let m = model()
        .skew_transport(2.0, [6.0, -4.0, 0.0])
        .unwrap()
        .operator_ns_per_m;
    assert!(m[0][1] > 0.0 && m[1][0] < 0.0);
    assert!((m[0][1] + m[1][0]).abs() > 0.1);
}

#[test]
fn uniform_thrust_equilibrium_has_the_source_downstream_gradient() {
    let incident = [6.0, -4.0, 0.0];
    let mean = 2.0;
    let speed = 4.0_f64.hypot(6.0 + mean);
    let half = 4.0 / (speed + 6.0 + mean);
    let h = 2.0 * (15.0 * std::f64::consts::PI / 64.0) * half * mean;
    let thrust = 2.0 * RHO * std::f64::consts::PI * R * R * speed * mean;
    let equilibrium = state([mean, h, 0.0]);
    let rate = model()
        .rate_skewed(equilibrium, incident, [thrust, 0.0, 0.0], 0.0)
        .unwrap();
    near_vec(rate.modal.acceleration_mps2, [0.0; 3], 8e-14);
    let step = model()
        .step_skewed_held(
            equilibrium,
            incident,
            [thrust, 0.0, 0.0],
            0.0,
            0.03,
            controls(),
        )
        .unwrap();
    near_vec(step.modal.induction.modes(), equilibrium.modes(), 3e-15);
    assert!(h > 0.0); // high induction is downstream, not on the advancing side.
    let reverse = model()
        .rate_skewed(
            state([mean, -h, 0.0]),
            [6.0, 4.0, 0.0],
            [thrust, 0.0, 0.0],
            0.0,
        )
        .unwrap();
    near_vec(reverse.modal.acceleration_mps2, [0.0; 3], 8e-14);
}

#[test]
fn zero_and_vanishing_crossflow_keep_the_existing_normal_solution() {
    let before = state([1.0, 0.2, -0.1]);
    let force = [140.0, 1.0, -2.0];
    let old = model().step_held(before, 5.0, force, 0.6, 0.003).unwrap();
    let exact = model()
        .step_skewed_held(before, [5.0, 0.0, 0.0], force, 0.6, 0.003, controls())
        .unwrap();
    assert_eq!(old, exact.modal);
    assert_eq!(exact.iterations, 0);
    for tiny in [1e-7, -1e-7, 1e-10, -1e-10] {
        let step = model()
            .step_skewed_held(
                before,
                [5.0, tiny, tiny * 0.7],
                force,
                0.6,
                0.003,
                controls(),
            )
            .unwrap();
        near_vec(
            step.modal.induction.modes(),
            old.induction.modes(),
            tiny.abs() * 0.2 + 2e-12,
        );
    }
}

// Small dense elimination is confined to the independent test oracle. It
// solves source H f=rho*A*u, rather than copying the production inverse.
fn solve(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> [f64; 3] {
    for k in 0..3 {
        let pivot = (k..3)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..3 {
            let ratio = a[i][k] / a[k][k];
            for j in k..3 {
                a[i][j] -= ratio * a[k][j];
            }
            b[i] -= ratio * b[k];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        x[i] = (b[i] - (i + 1..3).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
fn independent_rate(u: [f64; 3], incident: [f64; 3], q: [f64; 3], omega: f64) -> [f64; 3] {
    let rho_a = RHO * std::f64::consts::PI * R * R;
    let uw = wind_transform(incident, u);
    let f = solve(source_gain(u[0], incident), uw.map(|x| x * rho_a));
    let d = wind_transform(incident, [0.0, 1.0, 0.0]);
    // The columns of the orthogonal disk-to-wind transform give its inverse.
    let e = wind_transform(incident, [0.0, 0.0, 1.0]);
    let force = [f[0], d[1] * f[1] + d[2] * f[2], e[1] * f[1] + e[2] * f[2]];
    let masses = [
        8.0 / 3.0 * RHO * R.powi(3),
        16.0 / 45.0 * RHO * R.powi(3),
        16.0 / 45.0 * RHO * R.powi(3),
    ];
    [
        (q[0] - force[0]) / masses[0],
        (q[1] - force[1]) / masses[1] + omega * u[2],
        (q[2] - force[2]) / masses[2] - omega * u[1],
    ]
}

#[test]
fn nonlinear_midpoint_refines_to_independent_continuous_skew_dynamics() {
    let before = state([1.2, 0.1, -0.2]);
    let incident = [5.0, -4.0, 2.0];
    let force = [230.0, 4.0, -3.0];
    let omega = 0.7;
    let horizon = 0.08;
    let mut reference = before.modes();
    let dt = horizon / 8192.0;
    for _ in 0..8192 {
        let k1 = independent_rate(reference, incident, force, omega);
        let k2 = independent_rate(
            core::array::from_fn(|i| reference[i] + 0.5 * dt * k1[i]),
            incident,
            force,
            omega,
        );
        let k3 = independent_rate(
            core::array::from_fn(|i| reference[i] + 0.5 * dt * k2[i]),
            incident,
            force,
            omega,
        );
        let k4 = independent_rate(
            core::array::from_fn(|i| reference[i] + dt * k3[i]),
            incident,
            force,
            omega,
        );
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    let run = |n: usize| {
        let mut u = before;
        let mut max_iterations = 0;
        for _ in 0..n {
            let step = model()
                .step_skewed_held(u, incident, force, omega, horizon / n as f64, controls())
                .unwrap();
            max_iterations = max_iterations.max(step.iterations);
            assert!(step.modal.energy_residual_j.abs() < 3e-12);
            assert!(
                step.modal
                    .momentum_residual_ns
                    .into_iter()
                    .all(|x| x.abs() < 2e-12)
            );
            near(
                step.diagonal_transport_work_j + step.cross_transport_work_j,
                step.modal.transport_work_j,
                2e-14,
            );
            near(step.modal.frame_work_j, 0.0, 1e-17);
            assert!(step.cross_transport_work_j < 0.0);
            u = step.modal.induction;
        }
        let error = (0..3)
            .map(|i| (u.modes()[i] - reference[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        (error, max_iterations)
    };
    let results = [run(16), run(32), run(64)];
    println!("independent skew midpoint (error, max iterations): {results:?}");
    for pair in results.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0].0 / pair[1].0)));
    }
}

fn rebase(x: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    [x[0], c * x[1] + s * x[2], -s * x[1] + c * x[2]]
}
#[test]
fn turning_wind_changes_operator_without_reorienting_memory_and_is_frame_covariant() {
    let angle = 0.62;
    let mut a = state([1.0, 0.3, -0.2]);
    let mut b = state(rebase(a.modes(), angle));
    for tick in 0..60 {
        let t = (tick as f64 + 0.5) * 0.001;
        let incident = [5.0, 4.0 * (30.0 * t).cos(), 4.0 * (30.0 * t).sin()];
        let force = [150.0, 2.0, -1.0];
        let aa = model()
            .step_skewed_held(a, incident, force, 0.3, 0.001, controls())
            .unwrap();
        let bb = model()
            .step_skewed_held(
                b,
                rebase(incident, angle),
                rebase(force, angle),
                0.3,
                0.001,
                controls(),
            )
            .unwrap();
        near_vec(
            rebase(aa.modal.induction.modes(), angle),
            bb.modal.induction.modes(),
            2e-13,
        );
        near(aa.cross_transport_work_j, bb.cross_transport_work_j, 1e-16);
        a = aa.modal.induction;
        b = bb.modal.induction;
    }
    for vx in [-4.0, 4.0] {
        let s = model()
            .step_skewed_held(a, [5.0, vx, 0.0], [150.0, 2.0, -1.0], 0.3, 1e-9, controls())
            .unwrap();
        near_vec(a.modes(), s.modal.induction.modes(), 3e-8);
    }
}

#[test]
fn exhausted_nonlinear_budget_and_invalid_flow_never_return_partial_state() {
    let before = state([1.0, 0.2, -0.1]);
    let tiny = SkewInflowControls {
        maximum_iterations: 1,
        ..controls()
    };
    let error = model()
        .step_skewed_held(before, [5.0, -4.0, 1.0], [400.0, 3.0, 0.0], 0.4, 0.02, tiny)
        .unwrap_err();
    assert_eq!(error.code, "skew-inflow-unresolved");
    for n in [0, 65] {
        assert!(
            SkewInflowControls {
                maximum_iterations: n,
                ..controls()
            }
            .admit()
            .is_err()
        );
    }
    assert!(
        model()
            .rate_skewed(
                state([1.0, 10.0, 0.0]),
                [5.0, 2.0, 1.0],
                [100.0, 0.0, 0.0],
                0.0
            )
            .is_err()
    );
    assert!(
        model()
            .step_skewed_held(
                before,
                [5.0, -4.0, 1.0],
                [100.0, 1e7, 0.0],
                0.0,
                0.01,
                controls()
            )
            .is_err()
    );
    assert!(
        model()
            .skew_transport(1.0, [5.0, f64::INFINITY, 0.0])
            .is_err()
    );
    let run = || {
        model()
            .step_skewed_held(
                before,
                [5.0, -4.0, 1.0],
                [140.0, 3.0, 0.0],
                0.4,
                0.002,
                controls(),
            )
            .unwrap()
    };
    assert_eq!(run(), run());
}
