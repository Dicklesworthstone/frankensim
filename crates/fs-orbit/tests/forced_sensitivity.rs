use fs_math::{c64::C64, det};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales, IslandStateJacobian, OrbitProblem,
    sensitivity::{
        ForcedSensitivityControls, ForcedSensitivityDirection, ForcedSensitivityError,
        PreparedForcedSensitivity,
    },
    solve_hb_scaled_seeded,
};
use std::cell::Cell;

#[path = "forced_sensitivity/split.rs"]
mod split;

fn controls(samples: usize) -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: samples,
        maximum_unknowns: 256,
        maximum_payload_bytes: 16 * 1024 * 1024,
        maximum_work: 500_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_condition_estimate: 1e9,
        maximum_relative_residual: 1e-12,
        maximum_derivative_residual: 1e-10,
    }
}
fn direction(scale: f64) -> ForcedSensitivityDirection {
    ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: scale,
        angular_frequency_derivative: 0.0,
    }
}
fn scales(dimension: usize) -> HbScales {
    HbScales {
        states: vec![1.0; dimension],
        equations: vec![1.0; dimension],
        angular_frequency: 1.0,
    }
}

// Dissipative nonlinear coordinate coupled in BOTH directions to a relaxing
// material coordinate. Its odd harmonics and their parameter directions cannot
// be obtained by differentiating independent linear modal amplitudes.
struct Coupled {
    force: f64,
    beta: f64,
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}
impl Coupled {
    fn new(force: f64, beta: f64) -> Self {
        Self {
            force,
            beta,
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }
    fn rhs(&self, t: f64, state: &[f64], out: &mut [f64]) {
        out[0] = -1.2 * state[0] - self.beta * state[0].powi(3) - 0.4 * state[1]
            + self.force * det::cos(0.9 * t);
        out[1] = 0.3 * state[0] - 0.8 * state[1];
    }
    fn tangent(&self, _: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        out.copy_from_slice(&[-1.2 - 3.0 * self.beta * state[0].powi(2), -0.4, 0.3, -0.8]);
        true
    }
}
impl OrbitProblem for Coupled {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, state: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        self.rhs(t, state, out);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::tangent)
    }
}
fn solve(problem: &Coupled, harmonics: usize) -> HbOrbit {
    solve_hb_scaled_seeded(
        problem,
        HbAnchor::Forced { omega: 0.9 },
        vec![0.0; 2 * (2 * harmonics + 1)],
        &HbBudget {
            harmonics,
            max_newton: 40,
            tolerance: 1e-12,
            deflate_equilibrium: false,
        },
        &scales(2),
    )
    .unwrap()
    .orbit
}
fn amplitude_derivative(problem: &Coupled, harmonics: usize) -> (HbOrbit, Vec<Vec<C64>>) {
    let orbit = solve(problem, harmonics);
    let mut prepared =
        PreparedForcedSensitivity::new(problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    let source_calls = problem.source_calls.get();
    let tangent_calls = problem.tangent_calls.get();
    let response = prepared
        .solve_island_direction(
            direction(1.0),
            |t, _, out| {
                out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    // A second different physical direction reuses the factor and cached state.
    let material = prepared
        .solve_island_direction(
            direction(0.25),
            |_, x, out| {
                out.copy_from_slice(&[-x[0].powi(3), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    assert!(material.coeffs[0][3].abs() > 1e-4);
    assert_eq!(problem.source_calls.get(), source_calls);
    assert_eq!(problem.tangent_calls.get(), tangent_calls);
    (orbit, response.coeffs)
}

#[test]
fn coupled_nonlinear_harmonics_differentiate_force_and_material_without_resolving() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 11);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    let force = prepared
        .solve_island_direction(
            direction(1.0),
            |t, _, out| {
                out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    let material = prepared
        .solve_island_direction(
            direction(0.25),
            |_, x, out| {
                out.copy_from_slice(&[-x[0].powi(3), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    assert!(force.coeffs[0][3].abs() > 1e-3);
    assert!(force.coeffs[1][3].abs() > 1e-4);
    let step = 2e-5;
    for (analytical, low, high) in [
        (
            &force.coeffs,
            Coupled::new(1.1 - step, 1.0),
            Coupled::new(1.1 + step, 1.0),
        ),
        (
            &material.coeffs,
            Coupled::new(1.1, 1.0 - step),
            Coupled::new(1.1, 1.0 + step),
        ),
    ] {
        let low = solve(&low, 11);
        let high = solve(&high, 11);
        for state in 0..2 {
            for harmonic in 0..=11 {
                let finite_difference =
                    (high.coeffs[state][harmonic] - low.coeffs[state][harmonic]).scale(0.5 / step);
                assert!(
                    (analytical[state][harmonic] - finite_difference).abs() < 3e-10,
                    "state {state}, harmonic {harmonic}: {:?} versus {:?}",
                    analytical[state][harmonic],
                    finite_difference
                );
            }
        }
    }
    let (_, repeated) = amplitude_derivative(&problem, 11);
    assert_eq!(force.coeffs, repeated);
}

fn integrated_force_direction() -> Vec<Vec<C64>> {
    // Independent continuous-time RK4 on the original state AND its exact
    // variational ODE. No Fourier residual/Jacobian is used by this check.
    fn rhs(t: f64, z: &[f64; 4]) -> [f64; 4] {
        [
            -1.2 * z[0] - z[0].powi(3) - 0.4 * z[1] + 1.1 * det::cos(0.9 * t),
            0.3 * z[0] - 0.8 * z[1],
            (-1.2 - 3.0 * z[0] * z[0]) * z[2] - 0.4 * z[3] + det::cos(0.9 * t),
            0.3 * z[2] - 0.8 * z[3],
        ]
    }
    let samples = 2048;
    let dt = core::f64::consts::TAU / (0.9 * samples as f64);
    let mut z = [0.0; 4];
    let mut result = vec![vec![C64::ZERO; 16]; 2];
    for step in 0..(41 * samples) {
        let time = step as f64 * dt;
        if step >= 40 * samples {
            let phase = core::f64::consts::TAU * (step % samples) as f64 / samples as f64;
            for state in 0..2 {
                for harmonic in 0..16 {
                    let (sin, cos) = det::sin_cos(harmonic as f64 * phase);
                    result[state][harmonic] = result[state][harmonic]
                        + C64::new(cos, -sin).scale(z[state + 2] / samples as f64);
                }
            }
        }
        let k1 = rhs(time, &z);
        let k2 = rhs(
            time + 0.5 * dt,
            &core::array::from_fn(|i| z[i] + 0.5 * dt * k1[i]),
        );
        let k3 = rhs(
            time + 0.5 * dt,
            &core::array::from_fn(|i| z[i] + 0.5 * dt * k2[i]),
        );
        let k4 = rhs(time + dt, &core::array::from_fn(|i| z[i] + dt * k3[i]));
        for i in 0..4 {
            z[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    result
}

#[test]
fn harmonic_refinement_converges_to_independent_time_domain_variation() {
    let problem = Coupled::new(1.1, 1.0);
    let reference = integrated_force_direction();
    let mut errors = Vec::new();
    for harmonics in [3, 7, 11, 15] {
        let (_, derivative) = amplitude_derivative(&problem, harmonics);
        let mut maximum: f64 = 0.0;
        for state in 0..2 {
            for harmonic in 0..=harmonics {
                maximum =
                    maximum.max((derivative[state][harmonic] - reference[state][harmonic]).abs());
            }
        }
        errors.push(maximum);
    }
    assert!(errors[1] < errors[0] * 0.02, "{errors:?}");
    assert!(errors[2] < errors[1] * 0.02, "{errors:?}");
    assert!(errors[3] < errors[2] * 0.02, "{errors:?}");
    assert!(errors[3] < 3e-11, "{errors:?}");
}

struct FrequencyPort {
    omega: f64,
    forcing: Vec<C64>,
}
impl OrbitProblem for FrequencyPort {
    fn dim(&self) -> usize {
        1
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, _: &[f64], out: &mut [f64]) {
        out[0] = self.forcing[0].re;
        for (n, f) in self.forcing.iter().enumerate().skip(1) {
            let (s, c) = det::sin_cos(n as f64 * self.omega * t);
            out[0] += 2.0 * (f.re * c - f.im * s);
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![C64::from_re(2.0) + s.scale(0.4) + (s * s).scale(1.3)]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|_, _, _, out| {
            out[0] = 0.0;
            true
        })
    }
}
fn frequency_orbit(problem: &FrequencyPort) -> HbOrbit {
    HbOrbit {
        omega: problem.omega,
        coeffs: vec![
            problem
                .forcing
                .iter()
                .enumerate()
                .map(|(n, f)| *f / problem.port(C64::new(0.0, n as f64 * problem.omega))[0])
                .collect(),
        ],
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}

#[test]
fn complete_column_includes_frequency_port_clock_and_all_complex_harmonics() {
    let problem = FrequencyPort {
        omega: 0.8,
        forcing: vec![
            C64::from_re(0.3),
            C64::new(0.5, -0.2),
            C64::ZERO,
            C64::new(-0.1, 0.07),
        ],
    };
    let orbit = frequency_orbit(&problem);
    let prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(64), &scales(1), || false)
            .unwrap();
    let derivative = prepared
        .solve(
            ForcedSensitivityDirection {
                parameter: 0.8,
                parameter_scale: 0.1,
                angular_frequency_derivative: 1.0,
            },
            |point, out| {
                assert_eq!(point.parameter, 0.8);
                out.fill(0.0);
                for n in 1..=point.harmonics {
                    let x = C64::new(point.coefficients[2 * n - 1], point.coefficients[2 * n]);
                    let value = C64::new(
                        -2.6 * n as f64 * n as f64 * point.angular_frequency,
                        0.4 * n as f64,
                    ) * x;
                    out[2 * n - 1] = value.re;
                    out[2 * n] = value.im;
                }
                true
            },
            || false,
        )
        .unwrap();
    for n in 1..=3 {
        let z = problem.port(C64::new(0.0, n as f64 * problem.omega))[0];
        let expected = (C64::new(2.6 * (n * n) as f64 * problem.omega, -0.4 * n as f64)
            * problem.forcing[n])
            / (z * z);
        assert!((derivative.coeffs[0][n] - expected).abs() < 2e-13);
    }
    assert_eq!(derivative.coeffs[0][0], C64::ZERO);
    assert_eq!(derivative.angular_frequency_derivative, 1.0);
}

#[test]
fn parameter_scale_changes_conditioning_not_physical_derivative() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 7);
    let mut responses = Vec::new();
    for scale in [1e-12_f64, 1.0, 1e12] {
        // Absolute residual allowance follows the declared parameter increment.
        // The core returns derivatives per unit, never per scaled increment.
        let mut prepared = PreparedForcedSensitivity::new(
            &problem,
            &orbit,
            ForcedSensitivityControls {
                maximum_derivative_residual: 1e-10 * scale.max(1.0),
                ..controls(128)
            },
            &scales(2),
            || false,
        )
        .unwrap();
        responses.push(
            prepared
                .solve_island_direction(
                    direction(scale),
                    |t, _, out| {
                        out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                        true
                    },
                    || false,
                )
                .unwrap()
                .coeffs,
        );
    }
    for response in &responses[1..] {
        for (state, base) in response.iter().zip(&responses[0]) {
            for (a, b) in state.iter().zip(base) {
                assert!((*a - *b).abs() < 2e-14);
            }
        }
    }
}

#[test]
fn invalid_budgets_stale_orbits_incomplete_directions_and_cancellation_refuse() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 7);
    let calls = problem.source_calls.get();
    let refusal = PreparedForcedSensitivity::new(
        &problem,
        &orbit,
        ForcedSensitivityControls {
            maximum_payload_bytes: 1,
            ..controls(128)
        },
        &scales(2),
        || false,
    );
    assert!(matches!(
        refusal,
        Err(ForcedSensitivityError::Budget { .. })
    ));
    assert_eq!(calls, problem.source_calls.get());
    let mut stale = orbit.clone();
    stale.coeffs[0][1].re += 0.01;
    assert!(matches!(
        PreparedForcedSensitivity::new(&problem, &stale, controls(128), &scales(2), || false),
        Err(ForcedSensitivityError::BaseUnresolved { .. })
    ));
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    assert!(matches!(
        prepared.solve(
            direction(1.0),
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        prepared.solve_island_direction(
            direction(1.0),
            |_, _, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::IslandDirectionRefused { .. })
    ));
    let source_called = Cell::new(false);
    assert!(matches!(
        prepared.solve(
            direction(1.0),
            |_, _| {
                source_called.set(true);
                true
            },
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert!(!source_called.get());
    assert!(matches!(
        prepared.solve_island_direction(
            ForcedSensitivityDirection {
                angular_frequency_derivative: 1.0,
                ..direction(1.0)
            },
            |_, _, _| panic!(),
            || false
        ),
        Err(ForcedSensitivityError::Invalid(_))
    ));
    // A callback that changes its physical direction between grids must fail
    // independent differentiated-equation acceptance even after an exact solve.
    assert!(matches!(
        prepared.solve(
            direction(1.0),
            |point, out| {
                out.fill(if point.aft_samples == 128 { 1.0 } else { 2.0 });
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
}

#[test]
fn disparate_physical_state_and_equation_units_preserve_coupled_derivatives() {
    struct Units<'a>(&'a Coupled, [f64; 2]);
    impl OrbitProblem for Units<'_> {
        fn dim(&self) -> usize {
            2
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            self.0.rhs(t, &[x[0] / self.1[0], x[1] / self.1[1]], out);
            for (value, unit) in out.iter_mut().zip(self.1) {
                *value *= unit;
            }
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|p, t, x, out| {
                p.0.tangent(t, &[x[0] / p.1[0], x[1] / p.1[1]], out);
                for row in 0..2 {
                    for col in 0..2 {
                        out[2 * row + col] *= p.1[row] / p.1[col];
                    }
                }
                true
            })
        }
    }
    let original = Coupled::new(1.1, 1.0);
    let (mut orbit, expected) = amplitude_derivative(&original, 7);
    let units = Units(&original, [1e-6, 1e5]);
    for (row, unit) in orbit.coeffs.iter_mut().zip(units.1) {
        for value in row {
            *value = value.scale(unit);
        }
    }
    let mut prepared = PreparedForcedSensitivity::new(
        &units,
        &orbit,
        controls(128),
        &HbScales {
            states: units.1.to_vec(),
            equations: units.1.to_vec(),
            angular_frequency: 1.0,
        },
        || false,
    )
    .unwrap();
    let actual = prepared
        .solve_island_direction(
            direction(1.0),
            |t, _, out| {
                out.copy_from_slice(&[units.1[0] * det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    for state in 0..2 {
        for harmonic in 0..=7 {
            assert!(
                (actual.coeffs[state][harmonic].scale(1.0 / units.1[state])
                    - expected[state][harmonic])
                    .abs()
                    < 3e-14
            );
        }
    }
}

#[test]
fn singular_and_authored_condition_thresholds_refuse_without_regularization() {
    let problem = FrequencyPort {
        omega: 0.8,
        forcing: vec![C64::ZERO, C64::ZERO],
    };
    let orbit = frequency_orbit(&problem);
    assert!(matches!(
        PreparedForcedSensitivity::new(
            &problem,
            &orbit,
            ForcedSensitivityControls {
                maximum_condition_estimate: 1.0,
                ..controls(32)
            },
            &scales(1),
            || false
        ),
        Err(ForcedSensitivityError::IllConditioned { .. })
    ));
    struct Singular;
    impl OrbitProblem for Singular {
        fn dim(&self) -> usize {
            1
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, _: f64, _: &[f64], out: &mut [f64]) {
            out[0] = 0.0;
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, _, _, out| {
                out[0] = 0.0;
                true
            })
        }
    }
    assert!(matches!(
        PreparedForcedSensitivity::new(&Singular, &orbit, controls(32), &scales(1), || false),
        Err(ForcedSensitivityError::Singular)
    ));
}

// A physical observation mixes the two coupled coordinates BEFORE measuring
// mean square, retaining their coherent cross terms and all generated harmonics.
fn observation_mean_square(orbit: &HbOrbit) -> f64 {
    orbit.coeffs[0]
        .iter()
        .zip(&orbit.coeffs[1])
        .enumerate()
        .map(|(n, (&x, &z))| {
            let y = x + z.scale(0.7);
            (if n == 0 { 1.0 } else { 2.0 }) * (y.re * y.re + y.im * y.im)
        })
        .sum()
}

fn observation_gradient(orbit: &HbOrbit) -> Vec<f64> {
    let n = orbit.coeffs[0].len() - 1;
    let mut gradient = vec![0.0; 2 * (2 * n + 1)];
    for h in 0..=n {
        let y = orbit.coeffs[0][h] + orbit.coeffs[1][h].scale(0.7);
        for (state, weight) in [1.0, 0.7].into_iter().enumerate() {
            if h == 0 {
                gradient[state] = 2.0 * weight * y.re;
            } else {
                gradient[(2 * h - 1) * 2 + state] = 4.0 * weight * y.re;
                gradient[2 * h * 2 + state] = 4.0 * weight * y.im;
            }
        }
    }
    gradient
}

fn objective_direction(gradient: &[f64], derivative: &[Vec<C64>]) -> f64 {
    let mut total = 0.0;
    for h in 0..derivative[0].len() {
        for state in 0..2 {
            if h == 0 {
                total += gradient[state] * derivative[state][0].re;
            } else {
                total += gradient[(2 * h - 1) * 2 + state] * derivative[state][h].re
                    + gradient[2 * h * 2 + state] * derivative[state][h].im;
            }
        }
    }
    total
}

#[test]
fn forced_objective_adjoint_reuses_one_factor_for_nonlinear_force_and_material_gradients() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 11);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    let force = prepared
        .solve_island_direction(
            direction(1.0),
            |t, _, out| {
                out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    let material = prepared
        .solve_island_direction(
            direction(0.25),
            |_, x, out| {
                out.copy_from_slice(&[-x[0].powi(3), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    let gradient = observation_gradient(&orbit);
    let counts = (problem.source_calls.get(), problem.tangent_calls.get());
    let mut adjoint = prepared.prepare_adjoint(&gradient, 1.0, || false).unwrap();
    assert_eq!(adjoint.residual_weights().len(), gradient.len());
    assert!(adjoint.relative_residual < 1e-12 && adjoint.dual_residual < 1e-10);
    // q includes explicit terms 0.25*force + 0.5*beta^2 as well as the
    // complete observation's mean square. They are not state-chain terms.
    let d_force = adjoint
        .gradient_island_direction(
            direction(1.0),
            0.25,
            |t, _, out| {
                out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    let d_material = adjoint
        .gradient_island_direction(
            direction(0.25),
            1.0,
            |_, x, out| {
                out.copy_from_slice(&[-x[0].powi(3), 0.0]);
                true
            },
            || false,
        )
        .unwrap();
    assert!(
        (d_force.derivative - 0.25 - objective_direction(&gradient, &force.coeffs)).abs() < 3e-13
    );
    assert!(
        (d_material.derivative - 1.0 - objective_direction(&gradient, &material.coeffs)).abs()
            < 3e-13
    );
    // An independently authored coefficient direction gives the same answer
    // with two O(N) contractions and without any waveform sampling.
    let direct = adjoint
        .gradient(
            direction(1.0),
            0.25,
            |_, out| {
                out.fill(0.0);
                out[2] = -0.5;
                true
            },
            || false,
        )
        .unwrap();
    assert!((direct.derivative - d_force.derivative).abs() < 3e-13);
    assert_eq!(
        counts,
        (problem.source_calls.get(), problem.tangent_calls.get())
    );
    let step = 2e-5;
    for (got, low, high) in [
        (
            d_force.derivative,
            Coupled::new(1.1 - step, 1.0),
            Coupled::new(1.1 + step, 1.0),
        ),
        (
            d_material.derivative,
            Coupled::new(1.1, 1.0 - step),
            Coupled::new(1.1, 1.0 + step),
        ),
    ] {
        let objective = |p: &Coupled| {
            observation_mean_square(&solve(p, 11)) + 0.25 * p.force + 0.5 * p.beta * p.beta
        };
        let independently_solved = (objective(&high) - objective(&low)) / (2.0 * step);
        assert!(
            (got - independently_solved).abs() < 5e-10,
            "{got} vs {independently_solved}"
        );
    }
}

#[test]
fn objective_adjoint_retains_frequency_port_and_explicit_clock_dependence() {
    let problem = FrequencyPort {
        omega: 0.8,
        forcing: vec![
            C64::from_re(0.3),
            C64::new(0.5, -0.2),
            C64::ZERO,
            C64::new(-0.1, 0.07),
        ],
    };
    let orbit = frequency_orbit(&problem);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(64), &scales(1), || false)
            .unwrap();
    // Mean squared physical velocity: 2*sum(n^2*omega^2*|X_n|^2).
    let mut gradient = vec![0.0; 7];
    let mut explicit = 0.0;
    for n in 1..=3 {
        let x = orbit.coeffs[0][n];
        let scale = 4.0 * (n * n) as f64 * orbit.omega.powi(2);
        gradient[2 * n - 1] = scale * x.re;
        gradient[2 * n] = scale * x.im;
        explicit += 4.0 * (n * n) as f64 * orbit.omega * (x.re * x.re + x.im * x.im);
    }
    let adjoint = prepared.prepare_adjoint(&gradient, 1.0, || false).unwrap();
    let result = adjoint
        .gradient(
            ForcedSensitivityDirection {
                parameter: 0.8,
                parameter_scale: 0.1,
                angular_frequency_derivative: 1.0,
            },
            explicit,
            |point, out| {
                out.fill(0.0);
                for n in 1..=point.harmonics {
                    let x = C64::new(point.coefficients[2 * n - 1], point.coefficients[2 * n]);
                    let derivative = C64::new(
                        -2.6 * (n * n) as f64 * point.angular_frequency,
                        0.4 * n as f64,
                    ) * x;
                    out[2 * n - 1] = derivative.re;
                    out[2 * n] = derivative.im;
                }
                true
            },
            || false,
        )
        .unwrap();
    let objective = |w: f64| {
        let orbit = frequency_orbit(&FrequencyPort {
            omega: w,
            forcing: problem.forcing.clone(),
        });
        (1..=3)
            .map(|n| {
                let x = orbit.coeffs[0][n];
                2.0 * (n * n) as f64 * w * w * (x.re * x.re + x.im * x.im)
            })
            .sum::<f64>()
    };
    let fd = (objective(0.8 + 1e-5) - objective(0.8 - 1e-5)) / 2e-5;
    assert!(
        (result.derivative - fd).abs() < 3e-9,
        "{} vs {fd}",
        result.derivative
    );
}

#[test]
fn objective_adjoint_preserves_independent_state_equation_objective_and_parameter_units() {
    struct Units<'a>(&'a Coupled, [f64; 2], [f64; 2]);
    impl OrbitProblem for Units<'_> {
        fn dim(&self) -> usize {
            2
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            self.0.rhs(t, &[x[0] / self.1[0], x[1] / self.1[1]], out);
            for (v, scale) in out.iter_mut().zip(self.2) {
                *v *= scale;
            }
        }
        fn port(&self, s: C64) -> Vec<C64> {
            vec![
                s.scale(self.2[0] / self.1[0]),
                C64::ZERO,
                C64::ZERO,
                s.scale(self.2[1] / self.1[1]),
            ]
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|p, t, x, out| {
                p.0.tangent(t, &[x[0] / p.1[0], x[1] / p.1[1]], out);
                for row in 0..2 {
                    for col in 0..2 {
                        out[2 * row + col] *= p.2[row] / p.1[col];
                    }
                }
                true
            })
        }
    }
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 11);
    let gradient = observation_gradient(&orbit);
    let mut original =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    let mut adjoint = original.prepare_adjoint(&gradient, 1.0, || false).unwrap();
    let expected = adjoint
        .gradient_island_direction(
            direction(1.0),
            0.0,
            |t, _, out| {
                out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
                true
            },
            || false,
        )
        .unwrap()
        .derivative;
    let units = Units(&problem, [1e-6, 1e5], [1e3, 1e-2]);
    let mut converted = orbit.clone();
    for (row, unit) in converted.coeffs.iter_mut().zip(units.1) {
        for value in row {
            *value = value.scale(unit);
        }
    }
    let gradient: Vec<_> = gradient
        .iter()
        .enumerate()
        .map(|(i, v)| v / units.1[i % 2])
        .collect();
    for (objective_scale, parameter_scale) in [(0.01, 100.0), (1.0, 1.0), (100.0, 0.01)] {
        let mut prepared = PreparedForcedSensitivity::new(
            &units,
            &converted,
            controls(128),
            &HbScales {
                states: units.1.to_vec(),
                equations: units.2.to_vec(),
                angular_frequency: 1.0,
            },
            || false,
        )
        .unwrap();
        let mut adjoint = prepared
            .prepare_adjoint(&gradient, objective_scale, || false)
            .unwrap();
        let actual = adjoint
            .gradient_island_direction(
                direction(parameter_scale),
                0.0,
                |t, _, out| {
                    out.copy_from_slice(&[units.2[0] * det::cos(0.9 * t), 0.0]);
                    true
                },
                || false,
            )
            .unwrap();
        assert!((actual.derivative - expected).abs() < 3e-13);
    }
}

#[test]
fn objective_adjoint_refuses_incomplete_inconsistent_and_cancelled_directions() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 7);
    let gradient = observation_gradient(&orbit);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    assert!(matches!(
        prepared.prepare_adjoint(&gradient[..3], 1.0, || false),
        Err(ForcedSensitivityError::Invalid(_))
    ));
    assert!(matches!(
        prepared.prepare_adjoint(&gradient, 0.0, || false),
        Err(ForcedSensitivityError::Invalid(_))
    ));
    assert!(matches!(
        prepared.prepare_adjoint(&gradient, 1.0, || true),
        Err(ForcedSensitivityError::Cancelled)
    ));
    let mut adjoint = prepared.prepare_adjoint(&gradient, 1.0, || false).unwrap();
    assert!(matches!(
        adjoint.gradient(
            direction(1.0),
            0.0,
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        adjoint.gradient(
            direction(1.0),
            0.0,
            |point, out| {
                out.fill(if point.aft_samples == 128 { 1.0 } else { 2.0 });
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    let calls = Cell::new(0);
    assert!(matches!(
        adjoint.gradient(
            direction(1.0),
            0.0,
            |_, out| {
                calls.set(calls.get() + 1);
                out.fill(0.0);
                true
            },
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        adjoint.gradient_island_direction(
            ForcedSensitivityDirection {
                angular_frequency_derivative: 1.0,
                ..direction(1.0)
            },
            0.0,
            |_, _, _| panic!("clock refusal precedes callback"),
            || false
        ),
        Err(ForcedSensitivityError::Invalid(_))
    ));
    let zero = adjoint
        .gradient(
            direction(1.0),
            0.25,
            |_, out| {
                out.fill(0.0);
                true
            },
            || false,
        )
        .unwrap();
    assert_eq!(zero.derivative, 0.25);
}

#[test]
fn objective_adjoint_rechecks_original_dual_equations_on_both_grids() {
    struct PeriodicCoefficient;
    impl OrbitProblem for PeriodicCoefficient {
        fn dim(&self) -> usize {
            1
        }
        fn autonomous(&self) -> bool {
            false
        }
        fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = -(1.0 + 0.001 * det::cos(16.0 * t)) * x[0];
        }
        fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
            Some(|_, t, _, out| {
                out[0] = -(1.0 + 0.001 * det::cos(16.0 * t));
                true
            })
        }
    }
    let orbit = HbOrbit {
        omega: 1.0,
        coeffs: vec![vec![C64::ZERO; 4]],
        residual_trace: vec![],
        residual: 0.0,
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    // The source's 16th harmonic aliases to DC on the 16-point grid. Deliberately
    // admit that loose tangent comparison: the separate physical dual gate must
    // still refuse a finely solved lambda that fails the coarse equations.
    let mut prepared = PreparedForcedSensitivity::new(
        &PeriodicCoefficient,
        &orbit,
        ForcedSensitivityControls {
            maximum_tangent_sampling_defect: 0.1,
            ..controls(16)
        },
        &scales(1),
        || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.prepare_adjoint(&[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1.0, || false),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
}
