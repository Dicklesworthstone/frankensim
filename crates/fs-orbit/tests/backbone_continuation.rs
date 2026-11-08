//! Conservative continuation checked against original energy and kinematics.

use core::cell::Cell;
use fs_math::c64::C64;
use fs_orbit::{
    BackboneContinuationBudget, BackboneContinuationError, HbOrbit, HbScales, IslandStateJacobian,
    OrbitError, OrbitProblem, continue_backbone_branch,
};

const MIX: [f64; 9] = [2.0, 0.3, 0.1, 0.0, 3.0, 0.2, 0.4, 0.0, 5.0];

// Static/nodal state zero, a physical oscillator in states 1/2, dense equation
// mixing, and independently changeable original state and equation units.
struct Polynomial {
    k: f64,
    cubic: f64,
    quintic: f64,
    damping: f64,
    exact: bool,
    refuse: bool,
    state_units: [f64; 3],
    row_units: [f64; 3],
    calls: Cell<usize>,
}
impl Polynomial {
    fn new(k: f64, cubic: f64, quintic: f64) -> Self {
        Self {
            k,
            cubic,
            quintic,
            damping: 0.0,
            exact: true,
            refuse: false,
            state_units: [1.0; 3],
            row_units: [1.0; 3],
            calls: Cell::new(0),
        }
    }
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        if self.refuse {
            return false;
        }
        let q = x[1] / self.state_units[1];
        let derivative = -self.k - 3.0 * self.cubic * q * q - 5.0 * self.quintic * q.powi(4);
        let jac = [
            -7.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            derivative,
            -self.damping,
        ];
        for i in 0..3 {
            for j in 0..3 {
                out[3 * i + j] = self.row_units[i] / self.state_units[j]
                    * (0..3).map(|k| MIX[3 * i + k] * jac[3 * k + j]).sum::<f64>();
            }
        }
        true
    }
}
impl OrbitProblem for Polynomial {
    fn dim(&self) -> usize {
        3
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        let q = x[1] / self.state_units[1];
        let v = x[2] / self.state_units[2];
        let flow = [
            -7.0 * (x[0] / self.state_units[0] - 4.0),
            v,
            -self.k * q - self.cubic * q.powi(3) - self.quintic * q.powi(5) - self.damping * v,
        ];
        for i in 0..3 {
            out[i] = self.row_units[i] * (0..3).map(|j| MIX[3 * i + j] * flow[j]).sum::<f64>();
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        (0..9)
            .map(|i| s.scale(self.row_units[i / 3] * MIX[i] / self.state_units[i % 3]))
            .collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.exact.then_some(Self::jacobian)
    }
}

fn budget(problem: &Polynomial, harmonics: usize, steps: usize) -> BackboneContinuationBudget {
    BackboneContinuationBudget {
        harmonics,
        steps,
        max_newton: 35,
        max_attempts_per_step: 8,
        step_sizes: [0.07, 0.002, 0.09],
        tolerance: 2e-10,
        minimum_oscillation: 1e-8,
        maximum_unknowns: 512,
        maximum_payload_bytes: 128 * 1024 * 1024,
        maximum_residual_evaluations: 40_000,
        maximum_condition_estimate: 1e8,
        scales: HbScales {
            states: problem.state_units.to_vec(),
            equations: (0..3).map(|i| 5.0 * problem.row_units[i]).collect(),
            angular_frequency: 1.0,
        },
    }
}
fn orbit(coeffs: Vec<Vec<C64>>, omega: f64) -> HbOrbit {
    // Deliberately false metadata: the continuation must independently correct.
    HbOrbit {
        coeffs,
        omega,
        residual: 0.0,
        residual_trace: vec![f64::NAN],
        aft_samples: 0,
        sampling_defect: 0.0,
    }
}
fn seed(problem: &Polynomial, order: usize, amplitude: f64, omega: f64, phase: f64) -> HbOrbit {
    let mut coeffs = vec![vec![C64::ZERO; order + 1]; 3];
    coeffs[0][0] = C64::from_re(4.0 * problem.state_units[0]);
    let q = C64::new(phase.cos(), phase.sin()).scale(0.5 * amplitude);
    coeffs[1][1] = q.scale(problem.state_units[1]);
    coeffs[2][1] = (C64::new(0.0, omega) * q).scale(problem.state_units[2]);
    orbit(coeffs, omega)
}

fn potential(k: f64, cubic: f64, quintic: f64, q: f64) -> f64 {
    0.5 * k * q * q + 0.25 * cubic * q.powi(4) + quintic * q.powi(6) / 6.0
}
fn energy_period(k: f64, cubic: f64, quintic: f64, q: f64, v: f64) -> f64 {
    let energy = 0.5 * v * v + potential(k, cubic, quintic, q);
    let (mut lo, mut hi) = (0.0, 1.0);
    while potential(k, cubic, quintic, hi) < energy {
        hi *= 2.0;
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if potential(k, cubic, quintic, mid) < energy {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let a = 0.5 * (lo + hi);
    let panels = 4096;
    let step = core::f64::consts::FRAC_PI_2 / panels as f64;
    let mut sum = 0.0;
    for j in 0..=panels {
        let s2 = (step * j as f64).sin().powi(2);
        let value = 1.0
            / (k + 0.5 * cubic * a * a * (1.0 + s2)
                + quintic * a.powi(4) * (1.0 + s2 + s2 * s2) / 3.0)
                .sqrt();
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

fn check_physics(
    problem: &Polynomial,
    point: &fs_orbit::BackboneContinuationPoint,
    tolerance: f64,
) {
    let orbit = &point.orbit;
    assert!((orbit.coeffs[0][0].re / problem.state_units[0] - 4.0).abs() < 1e-9);
    assert!(orbit.first_harmonic_amplitude(0) / problem.state_units[0] < 1e-10);
    let q = orbit.sample(1, 0.0) / problem.state_units[1];
    let v = orbit.sample(2, 0.0) / problem.state_units[2];
    let reference = energy_period(problem.k, problem.cubic, problem.quintic, q, v);
    assert!(
        (orbit.omega / reference - 1.0).abs() < tolerance,
        "frequency {} vs independent energy period {reference}",
        orbit.omega
    );
    for n in 1..orbit.coeffs[1].len() {
        let expected = C64::new(0.0, n as f64 * orbit.omega)
            * orbit.coeffs[1][n].scale(1.0 / problem.state_units[1]);
        let actual = orbit.coeffs[2][n].scale(1.0 / problem.state_units[2]);
        assert!((actual - expected).abs() < 2e-9);
    }
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    for j in 0..129 {
        let theta = core::f64::consts::TAU * j as f64 / 129.0;
        let q = orbit.sample(1, theta) / problem.state_units[1];
        let v = orbit.sample(2, theta) / problem.state_units[2];
        let energy = 0.5 * v * v + potential(problem.k, problem.cubic, problem.quintic, q);
        minimum = minimum.min(energy);
        maximum = maximum.max(energy);
    }
    assert!((maximum - minimum) / maximum < 10.0 * tolerance);
    assert!(point.scaled_residual <= 2e-10);
    assert!(point.scaled_sampling_defect <= 2e-10);
    assert!(orbit.aft_samples > 0);
    assert!(point.jacobian_condition_estimate.is_finite());
}

#[test]
fn full_state_dense_port_backbones_keep_phase_units_energy_and_every_velocity_harmonic() {
    for exact in [false, true] {
        for rescaled in [false, true] {
            let mut problem = Polynomial::new(4.0, 8.0, 0.0);
            problem.exact = exact;
            if rescaled {
                problem.state_units = [1e6, 1e-8, 1e8];
                problem.row_units = [1e-7, 1e5, 3.0];
            }
            let b = budget(&problem, 9, 5);
            let a = seed(&problem, 9, 0.35, 2.1, 0.83);
            let c = seed(&problem, 9, 0.4, 2.2, -1.21);
            let branch = continue_backbone_branch(
                &problem,
                1,
                [
                    (0.35 * problem.state_units[1], &a),
                    (0.4 * problem.state_units[1], &c),
                ],
                &b,
                || false,
            )
            .unwrap();
            assert_eq!(branch.points.len(), 7);
            assert_eq!(branch.coordinate, 1);
            assert_eq!(branch.tangent_preparations > 0, exact);
            assert!(branch.residual_evaluations < b.maximum_residual_evaluations);
            for point in &branch.points {
                check_physics(&problem, point, 3e-8);
            }
            for pair in branch.points.windows(2) {
                assert!(pair[1].first_harmonic_amplitude > pair[0].first_harmonic_amplitude);
                assert!(pair[1].orbit.omega > pair[0].orbit.omega);
            }
            assert_eq!(branch.points[0].step, 0.0);
            assert_eq!(branch.points[1].step, 0.0);
            assert!(branch.points[2].step > 0.0);
        }
    }
}

#[test]
fn true_arclength_crosses_softening_hardening_frequency_turn() {
    let problem = Polynomial::new(1.0, -0.4, 0.1);
    // The long hardening side needs H21: H15 balances its retained rows but
    // misses the energy-period oracle by 2.9e-7 at the last point.
    let b = budget(&problem, 21, 12);
    let a = seed(&problem, 21, 1.2, 0.84, 0.47);
    let c = seed(&problem, 21, 1.25, 0.82, -0.93);
    let branch =
        continue_backbone_branch(&problem, 1, [(1.2, &a), (1.25, &c)], &b, || false).unwrap();
    let mut softening = false;
    let mut hardening = false;
    for point in &branch.points {
        check_physics(&problem, point, 2e-7);
    }
    for pair in branch.points.windows(2) {
        let difference = pair[1].orbit.omega - pair[0].orbit.omega;
        softening |= difference < -1e-3;
        hardening |= difference > 1e-3;
        assert!(pair[1].first_harmonic_amplitude > pair[0].first_harmonic_amplitude);
    }
    assert!(
        softening && hardening,
        "frequencies: {:?}",
        branch
            .points
            .iter()
            .map(|p| p.orbit.omega)
            .collect::<Vec<_>>()
    );
    assert!(
        branch
            .points
            .iter()
            .any(|p| p.orbit.coeffs[1][3].abs() > 1e-3)
    );
}

struct Coupled {
    identical: bool,
}
impl Coupled {
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        let coupling = if self.identical {
            0.0
        } else {
            0.3 + 1.8 * (x[0] - x[2]).powi(2)
        };
        out[1] = 1.0;
        out[4] = -1.0 - coupling;
        out[6] = coupling;
        out[11] = 1.0;
        out[12] = coupling;
        out[14] = -1.0 - coupling;
        true
    }
}
impl OrbitProblem for Coupled {
    fn dim(&self) -> usize {
        4
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        let difference = x[0] - x[2];
        let force = if self.identical {
            0.0
        } else {
            0.3 * difference + 0.6 * difference.powi(3)
        };
        out.copy_from_slice(&[x[1], -x[0] - force, x[3], -x[2] + force]);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}
fn coupled_seed(order: usize, amplitude: f64, omega: f64) -> HbOrbit {
    let mut coeffs = vec![vec![C64::ZERO; order + 1]; 4];
    coeffs[0][1] = C64::from_re(0.5 * amplitude);
    coeffs[1][1] = C64::new(0.0, 0.5 * amplitude * omega);
    coeffs[2][1] = coeffs[0][1].scale(-1.0);
    coeffs[3][1] = coeffs[1][1].scale(-1.0);
    orbit(coeffs, omega)
}

#[test]
fn coupled_two_mass_family_retains_reciprocal_nonlinear_force() {
    let problem = Coupled { identical: false };
    let mut b = budget(&Polynomial::new(1.0, 0.0, 0.0), 9, 4);
    b.scales.states = vec![1.0; 4];
    b.scales.equations = vec![1.0; 4];
    let a = coupled_seed(9, 0.2, 1.3);
    let c = coupled_seed(9, 0.25, 1.35);
    let result =
        continue_backbone_branch(&problem, 0, [(0.2, &a), (0.25, &c)], &b, || false).unwrap();
    for point in result.points {
        let orbit = point.orbit;
        for n in 0..=9 {
            assert!((orbit.coeffs[0][n] + orbit.coeffs[2][n]).abs() < 1e-9);
            assert!((orbit.coeffs[1][n] + orbit.coeffs[3][n]).abs() < 1e-9);
        }
        let reference = energy_period(1.6, 4.8, 0.0, orbit.sample(0, 0.0), orbit.sample(1, 0.0));
        assert!((orbit.omega / reference - 1.0).abs() < 3e-8);
    }
}

#[test]
fn exact_residual_seeds_still_refuse_additional_neutral_directions() {
    let problem = Coupled { identical: true };
    let mut b = budget(&Polynomial::new(1.0, 0.0, 0.0), 3, 0);
    b.scales.states = vec![1.0; 4];
    b.scales.equations = vec![1.0; 4];
    let a = coupled_seed(3, 0.2, 1.0);
    let c = coupled_seed(3, 0.25, 1.0);
    assert!(matches!(
        continue_backbone_branch(&problem, 0, [(0.2, &a), (0.25, &c)], &b, || false),
        Err(BackboneContinuationError::Numerical(
            OrbitError::SingularJacobian
        )) | Err(BackboneContinuationError::IllConditioned { .. })
    ));
}

#[test]
fn nonzero_damped_seed_is_not_turned_into_a_conservative_family() {
    let mut problem = Polynomial::new(4.0, 8.0, 0.0);
    problem.damping = 0.2;
    let b = budget(&problem, 7, 2);
    let a = seed(&problem, 7, 0.3, 2.1, 0.0);
    let c = seed(&problem, 7, 0.35, 2.2, 0.0);
    assert!(matches!(
        continue_backbone_branch(&problem, 1, [(0.3, &a), (0.35, &c)], &b, || false),
        Err(BackboneContinuationError::Numerical(
            OrbitError::NewtonStalled { .. }
        ))
    ));
}

#[test]
fn admission_cancellation_evaluation_limits_and_authored_tangent_refusal_are_enforced() {
    let mut problem = Polynomial::new(4.0, 8.0, 0.0);
    let original = budget(&problem, 3, 1);
    let a = seed(&problem, 3, 0.3, 2.1, 0.0);
    let c = seed(&problem, 3, 0.35, 2.2, 0.0);
    let seeds = [(0.3, &a), (0.35, &c)];
    let mut b = original.clone();
    b.maximum_payload_bytes = 1;
    assert!(matches!(
        continue_backbone_branch(&problem, 1, seeds, &b, || false),
        Err(BackboneContinuationError::Invalid(_))
    ));
    b = original.clone();
    b.maximum_unknowns = 1;
    assert!(matches!(
        continue_backbone_branch(&problem, 1, seeds, &b, || false),
        Err(BackboneContinuationError::Invalid(_))
    ));
    assert!(matches!(
        continue_backbone_branch(&problem, 1, seeds, &original, || true),
        Err(BackboneContinuationError::Cancelled)
    ));
    assert_eq!(problem.calls.get(), 0);
    b = original.clone();
    b.maximum_residual_evaluations = 1;
    assert_eq!(
        continue_backbone_branch(&problem, 1, seeds, &b, || false).unwrap_err(),
        BackboneContinuationError::EvaluationLimit { evaluations: 1 }
    );
    assert_eq!(problem.calls.get(), 16);
    problem.refuse = true;
    assert!(matches!(
        continue_backbone_branch(&problem, 1, seeds, &original, || false),
        Err(BackboneContinuationError::Numerical(
            OrbitError::TangentRefused { sample: 0 }
        ))
    ));
    problem.refuse = false;
    let mut polls = 0;
    assert!(matches!(
        continue_backbone_branch(&problem, 1, seeds, &original, || {
            polls += 1;
            polls > 10
        }),
        Err(BackboneContinuationError::Cancelled)
    ));
}
