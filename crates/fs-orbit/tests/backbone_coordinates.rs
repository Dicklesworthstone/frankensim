//! Original-coordinate backbone physics, independent of state/row ordering.

use core::cell::Cell;
use fs_math::c64::C64;
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbRefinementBudget, HbScales, IslandStateJacobian, OrbitError,
    OrbitProblem, solve_hb_backbone_coordinate_adaptive_scaled_seeded,
    solve_hb_backbone_coordinate_scaled_seeded, solve_hb_scaled_seeded,
};

// A static first coordinate and a conservative nonlinear oscillator in the
// remaining coordinates. A dense, nonsingular equation multiplier makes both
// rows AND columns of the original frequency port/tangent matter.
const MASS: [f64; 9] = [2.0, 0.3, 0.1, 0.0, 3.0, 0.2, 0.4, 0.0, 5.0];
struct HiddenMode {
    beta: f64,
    exact: bool,
    refuse_tangent: bool,
    autonomous: bool,
    calls: Cell<usize>,
    tangents: Cell<usize>,
}
impl HiddenMode {
    fn new(beta: f64) -> Self {
        Self {
            beta,
            exact: true,
            refuse_tangent: false,
            autonomous: true,
            calls: Cell::new(0),
            tangents: Cell::new(0),
        }
    }
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        self.tangents.set(self.tangents.get() + 1);
        if self.refuse_tangent {
            return false;
        }
        let jac = [
            -7.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            -4.0 - 3.0 * self.beta * x[1] * x[1],
            0.0,
        ];
        for i in 0..3 {
            for j in 0..3 {
                out[3 * i + j] = (0..3).map(|k| MASS[3 * i + k] * jac[3 * k + j]).sum();
            }
        }
        true
    }
}
impl OrbitProblem for HiddenMode {
    fn dim(&self) -> usize {
        3
    }
    fn autonomous(&self) -> bool {
        self.autonomous
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        let flow = [
            -7.0 * (x[0] - 4.0),
            x[2],
            -4.0 * x[1] - self.beta * x[1].powi(3),
        ];
        for i in 0..3 {
            out[i] = (0..3).map(|j| MASS[3 * i + j] * flow[j]).sum();
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        MASS.iter().map(|&entry| s.scale(entry)).collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.exact.then_some(Self::jacobian)
    }
}

fn budget(order: usize) -> HbBudget {
    HbBudget {
        harmonics: order,
        max_newton: 50,
        tolerance: 2e-11,
        deflate_equilibrium: true,
    }
}
fn scales() -> HbScales {
    HbScales {
        states: vec![4.0, 0.4, 1.0],
        equations: vec![2.0, 3.0, 10.0],
        angular_frequency: 2.0,
    }
}
fn seed(order: usize, amplitude: f64, omega: f64, phase: f64) -> Vec<f64> {
    let mut u = vec![0.0; 3 * (2 * order + 1) + 1];
    u[0] = 4.0;
    u[4] = 0.5 * amplitude * phase.cos();
    u[7] = 0.5 * amplitude * phase.sin();
    // q=A*cos(theta+phase), v=-omega*A*sin(theta+phase).
    u[5] = -0.5 * omega * amplitude * phase.sin();
    u[8] = 0.5 * omega * amplitude * phase.cos();
    *u.last_mut().unwrap() = omega;
    u
}

fn independent_frequency(orbit: &HbOrbit, beta: f64) -> f64 {
    let q = orbit.sample(1, 0.0);
    let v = orbit.sample(2, 0.0);
    let energy = 2.0 * q * q + 0.25 * beta * q.powi(4) + 0.5 * v * v;
    let a_squared = 4.0 * energy / (4.0 + (16.0 + 4.0 * beta * energy).sqrt());
    // Energy quadrature dt=dphi/sqrt(omega0²+beta*A²*(1+sin²(phi))/2).
    let panels = 4096;
    let step = core::f64::consts::FRAC_PI_2 / panels as f64;
    let mut sum = 0.0;
    for j in 0..=panels {
        let sine = (j as f64 * step).sin();
        let value = 1.0 / (4.0 + 0.5 * beta * a_squared * (1.0 + sine * sine)).sqrt();
        let weight = if j == 0 || j == panels {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * value;
    }
    core::f64::consts::TAU / (4.0 * step * sum / 3.0)
}

#[test]
fn nonlinear_mode_survives_static_first_state_dense_port_and_phase_rotation() {
    for exact in [false, true] {
        for phase in [0.0, 0.73] {
            let mut problem = HiddenMode::new(8.0);
            problem.exact = exact;
            let b = budget(9);
            let result = solve_hb_backbone_coordinate_scaled_seeded(
                &problem,
                1,
                0.4,
                seed(b.harmonics, 0.4, 2.1, phase),
                &b,
                &scales(),
            )
            .unwrap();
            let orbit = result.orbit;
            assert!((orbit.coeffs[0][0].re - 4.0).abs() < 1e-10);
            assert!(orbit.first_harmonic_amplitude(0) < 1e-11);
            assert!((orbit.first_harmonic_amplitude(1) - 0.4).abs() < 2e-11);
            assert!(orbit.coeffs[1][3].abs() > 1e-4);
            let reference = independent_frequency(&orbit, problem.beta);
            assert!(
                (orbit.omega / reference - 1.0).abs() < 2e-8,
                "frequency {} vs original energy period {reference}",
                orbit.omega
            );
            for n in 1..=b.harmonics {
                let expected_velocity = C64::new(0.0, n as f64 * orbit.omega) * orbit.coeffs[1][n];
                assert!((orbit.coeffs[2][n] - expected_velocity).abs() < 2e-10);
            }
            let mut minimum = f64::INFINITY;
            let mut maximum = 0.0_f64;
            for j in 0..257 {
                let phase = core::f64::consts::TAU * j as f64 / 257.0;
                let q = orbit.sample(1, phase);
                let v = orbit.sample(2, phase);
                let energy = 2.0 * q * q + 0.25 * problem.beta * q.powi(4) + 0.5 * v * v;
                minimum = minimum.min(energy);
                maximum = maximum.max(energy);
            }
            assert!((maximum - minimum) / maximum < 2e-8);
            assert_eq!(problem.tangents.get() > 0, exact);
        }
    }
}

#[test]
fn adaptive_orders_keep_the_chosen_physical_coordinate_and_full_equations() {
    let problem = HiddenMode::new(8.0);
    let b = budget(1);
    let limits = HbRefinementBudget {
        max_harmonics: 15,
        residual_tolerances: vec![2e-7, 3e-7, 1e-6],
        state_tolerances: vec![1e-8, 1e-7, 3e-7],
        frequency_tolerance: 1e-7,
    };
    let result = solve_hb_backbone_coordinate_adaptive_scaled_seeded(
        &problem,
        1,
        0.4,
        seed(1, 0.4, 2.1, 0.31),
        &b,
        &limits,
        &scales(),
    )
    .unwrap();
    assert!(result.refinement_trace.len() >= 3);
    assert!(result.orbit.coeffs[0].len() > 4);
    assert!(result.orbit.first_harmonic_amplitude(0) < 1e-11);
    assert!((result.orbit.first_harmonic_amplitude(1) - 0.4).abs() < 2e-11);
    let last = result.refinement_trace.last().unwrap();
    assert!(last.total_residual_ratio <= 1.0);
    assert!(last.state_change_ratio.unwrap() <= 1.0);
    assert!(last.relative_frequency_change.unwrap() <= limits.frequency_tolerance);
    assert!((result.orbit.omega / independent_frequency(&result.orbit, 8.0) - 1.0).abs() < 2e-8);
}

#[test]
fn last_coordinate_can_pin_physical_velocity_without_moving_output_rows() {
    let problem = HiddenMode::new(0.0);
    let b = budget(3);
    let result = solve_hb_backbone_coordinate_scaled_seeded(
        &problem,
        2,
        0.8,
        seed(3, 0.4, 1.9, 0.63),
        &b,
        &scales(),
    )
    .unwrap();
    assert!((result.orbit.omega - 2.0).abs() < 1e-9);
    assert!((result.orbit.first_harmonic_amplitude(2) - 0.8).abs() < 1e-10);
    assert!((result.orbit.first_harmonic_amplitude(1) - 0.4).abs() < 1e-10);
    assert!((result.orbit.coeffs[0][0].re - 4.0).abs() < 1e-10);
}

struct FirstMode;
impl OrbitProblem for FirstMode {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = -4.0 * x[0];
    }
}

#[test]
fn original_first_coordinate_path_is_bitwise_unchanged() {
    let problem = FirstMode;
    let b = budget(3);
    let s = HbScales {
        states: vec![0.4, 0.8],
        equations: vec![0.8, 1.6],
        angular_frequency: 2.0,
    };
    let mut u = vec![0.0; 2 * (2 * b.harmonics + 1) + 1];
    u[2] = 0.2;
    u[5] = 0.38;
    *u.last_mut().unwrap() = 1.9;
    let expected = solve_hb_scaled_seeded(
        &problem,
        HbAnchor::Backbone {
            amplitude: 0.4,
            omega_guess: 1.9,
        },
        u.clone(),
        &b,
        &s,
    )
    .unwrap();
    let actual = solve_hb_backbone_coordinate_scaled_seeded(&problem, 0, 0.4, u, &b, &s).unwrap();
    assert_eq!(actual.orbit.coeffs, expected.orbit.coeffs);
    assert_eq!(actual.orbit.omega.to_bits(), expected.orbit.omega.to_bits());
    assert_eq!(actual.orbit.residual_trace, expected.orbit.residual_trace);
    assert_eq!(
        actual.scaled_residual.to_bits(),
        expected.scaled_residual.to_bits()
    );
}

#[test]
fn invalid_charts_refuse_before_source_work_and_exact_tangent_refusal_survives() {
    let b = budget(3);
    let problem = HiddenMode::new(8.0);
    for coordinate in [3, usize::MAX] {
        assert!(matches!(
            solve_hb_backbone_coordinate_scaled_seeded(
                &problem,
                coordinate,
                0.4,
                seed(3, 0.4, 2.1, 0.0),
                &b,
                &scales(),
            ),
            Err(OrbitError::BadParameter { .. })
        ));
    }
    for amplitude in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            solve_hb_backbone_coordinate_scaled_seeded(
                &problem,
                1,
                amplitude,
                seed(3, 0.4, 2.1, 0.0),
                &b,
                &scales(),
            )
            .is_err()
        );
    }
    assert!(
        solve_hb_backbone_coordinate_scaled_seeded(&problem, 1, 0.4, vec![], &b, &scales(),)
            .is_err()
    );
    let mut malformed = scales();
    malformed.equations.pop();
    assert!(
        solve_hb_backbone_coordinate_scaled_seeded(
            &problem,
            1,
            0.4,
            seed(3, 0.4, 2.1, 0.0),
            &b,
            &malformed,
        )
        .is_err()
    );
    assert_eq!(problem.calls.get(), 0);
    assert_eq!(problem.tangents.get(), 0);
    let mut problem = HiddenMode::new(8.0);
    problem.autonomous = false;
    assert!(
        solve_hb_backbone_coordinate_scaled_seeded(
            &problem,
            1,
            0.4,
            seed(3, 0.4, 2.1, 0.0),
            &b,
            &scales(),
        )
        .is_err()
    );
    assert_eq!(problem.calls.get(), 0);
    problem.autonomous = true;
    problem.refuse_tangent = true;
    assert!(matches!(
        solve_hb_backbone_coordinate_scaled_seeded(
            &problem,
            1,
            0.4,
            seed(3, 0.4, 2.1, 0.0),
            &b,
            &scales(),
        ),
        Err(OrbitError::TangentRefused { sample: 0 })
    ));
    assert_eq!(problem.tangents.get(), 1);
}
