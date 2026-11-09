//! Conservative QR refinement against continuous energy and original equations.

use fs_math::c64::C64;
use fs_orbit::{
    AdaptiveBackboneOrbit, BackboneContinuationError, BackboneRefinementControls, HbOrbit,
    HbRefinementBudget, HbScales, IslandStateJacobian, OrbitError, OrbitProblem,
    solve_backbone_adaptive,
};
use std::cell::Cell;

const MIX: [f64; 9] = [2.0, 0.3, 0.1, 0.0, 3.0, 0.2, 0.4, 0.0, 5.0];
struct Polynomial {
    cubic: f64,
    quintic: f64,
    offset: f64,
    state: [f64; 3],
    equation: [f64; 3],
    exact: bool,
    refused: bool,
    calls: Cell<usize>,
}
impl Polynomial {
    fn new() -> Self {
        Self {
            cubic: 80.0,
            quintic: 200.0,
            offset: 4.0,
            state: [1.0; 3],
            equation: [1.0; 3],
            exact: true,
            refused: false,
            calls: Cell::new(0),
        }
    }
    fn scales(&self) -> HbScales {
        HbScales {
            states: self.state.to_vec(),
            equations: self.equation.iter().map(|v| 5.0 * v).collect(),
            angular_frequency: 4.0,
        }
    }
    fn refinement(&self, maximum: usize) -> HbRefinementBudget {
        HbRefinementBudget {
            max_harmonics: maximum,
            residual_tolerances: self.equation.iter().map(|v| 2e-8 * v).collect(),
            state_tolerances: self.state.iter().map(|v| 2e-8 * v).collect(),
            frequency_tolerance: 2e-9,
        }
    }
    fn seed(&self) -> HbOrbit {
        let amplitude: f64 = 0.42;
        let omega = (4.0
            + 0.75 * self.cubic * amplitude.powi(2)
            + 0.625 * self.quintic * amplitude.powi(4))
        .sqrt();
        let phase: f64 = 0.73;
        let q = C64::new(phase.cos(), phase.sin()).scale(0.5 * amplitude);
        HbOrbit {
            coeffs: vec![
                vec![C64::from_re(self.offset * self.state[0]), C64::ZERO],
                vec![C64::ZERO, q.scale(self.state[1])],
                vec![C64::ZERO, (C64::new(0.0, omega) * q).scale(self.state[2])],
            ],
            omega,
            // None of this caller-authored convergence metadata is trusted.
            residual: 0.0,
            residual_trace: vec![f64::NAN],
            aft_samples: usize::MAX,
            sampling_defect: 0.0,
        }
    }
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        if self.refused {
            return false;
        }
        let q = x[1] / self.state[1];
        let slope = -4.0 - 3.0 * self.cubic * q * q - 5.0 * self.quintic * q.powi(4);
        let jac = [-7.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, slope, 0.0];
        for row in 0..3 {
            for col in 0..3 {
                out[3 * row + col] = self.equation[row] / self.state[col]
                    * (0..3)
                        .map(|k| MIX[3 * row + k] * jac[3 * k + col])
                        .sum::<f64>();
            }
        }
        true
    }
    fn solve(&self) -> AdaptiveBackboneOrbit {
        solve_backbone_adaptive(
            self,
            1,
            0.42 * self.state[1],
            &self.seed(),
            &controls(),
            &self.refinement(31),
            &self.scales(),
            || false,
        )
        .unwrap()
    }
    fn potential(&self, q: f64) -> f64 {
        2.0 * q * q + 0.25 * self.cubic * q.powi(4) + self.quintic * q.powi(6) / 6.0
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
        let q = x[1] / self.state[1];
        let v = x[2] / self.state[2];
        let g = [
            -7.0 * (x[0] / self.state[0] - self.offset),
            v,
            -4.0 * q - self.cubic * q.powi(3) - self.quintic * q.powi(5),
        ];
        for row in 0..3 {
            out[row] = self.equation[row] * (0..3).map(|k| MIX[3 * row + k] * g[k]).sum::<f64>();
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        (0..9)
            .map(|i| s.scale(self.equation[i / 3] * MIX[i] / self.state[i % 3]))
            .collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        self.exact.then_some(Self::jacobian)
    }
}
fn controls() -> BackboneRefinementControls {
    BackboneRefinementControls {
        max_newton: 70,
        tolerance: 2e-12,
        minimum_oscillation: 1e-8,
        maximum_unknowns: 256,
        maximum_payload_bytes: 64 * 1024 * 1024,
        maximum_residual_evaluations: 30_000,
        maximum_condition_estimate: 1e10,
    }
}

// Continuous energy quadrature is independent of harmonic balance and uses
// the ORIGINAL polynomial potential, including the generated fifth harmonic.
fn energy_frequency(model: &Polynomial, orbit: &HbOrbit) -> f64 {
    let q = orbit.sample(1, 0.0) / model.state[1];
    let v = orbit.sample(2, 0.0) / model.state[2];
    let energy = model.potential(q) + 0.5 * v * v;
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if model.potential(mid) < energy {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let a = 0.5 * (lo + hi);
    let h = std::f64::consts::FRAC_PI_2 / 4096.0;
    let mut sum = 0.0;
    for j in 0..=4096 {
        let s2 = (j as f64 * h).sin().powi(2);
        let inverse_speed = 1.0
            / (4.0
                + 0.5 * model.cubic * a * a * (1.0 + s2)
                + model.quintic * a.powi(4) * (1.0 + s2 + s2 * s2) / 3.0)
                .sqrt();
        let weight = if j == 0 || j == 4096 {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * inverse_speed;
    }
    std::f64::consts::TAU / (4.0 * h * sum / 3.0)
}

#[test]
fn original_rectangular_balance_refines_generated_modes_and_matches_continuous_period() {
    let model = Polynomial::new();
    let result = model.solve();
    let orbit = &result.point.orbit;
    assert!(result.refinement_trace.len() >= 4);
    let first = &result.refinement_trace[0];
    assert_eq!(first.harmonics, 1);
    assert!(first.retained_residual_ratio < 1.0);
    assert!(
        first.omitted_residual_ratio > 1e4,
        "a solved short truncation must not imply physical resolution"
    );
    let last = result.refinement_trace.last().unwrap();
    assert!(last.total_residual_ratio <= 1.0 && last.omitted_residual_ratio <= 1.0);
    assert!(last.state_change_ratio.unwrap() <= 1.0);
    assert!(last.relative_frequency_change.unwrap() <= model.refinement(31).frequency_tolerance);
    assert!(last.sampling_defect_ratio <= 1.0);
    assert!((orbit.first_harmonic_amplitude(1) - 0.42).abs() < 2e-12);
    assert!(orbit.coeffs[1][3].abs() > 1e-3 && orbit.coeffs[1][5].abs() > 1e-5);
    assert!((orbit.omega / energy_frequency(&model, orbit) - 1.0).abs() < 3e-9);
    assert!(result.point.phase_residual.abs() <= controls().tolerance);
    assert!(result.point.jacobian_condition_estimate.is_finite());
    assert!(
        result.tangent_preparations > 0
            && result.residual_evaluations > result.tangent_preparations
    );
    let (mut low, mut high) = (f64::INFINITY, 0.0_f64);
    for j in 0..257 {
        let theta = std::f64::consts::TAU * j as f64 / 257.0;
        let q = orbit.sample(1, theta);
        let v = orbit.sample(2, theta);
        let energy = model.potential(q) + 0.5 * v * v;
        low = low.min(energy);
        high = high.max(energy);
    }
    assert!((high - low) / high < 2e-8);
    for n in 1..orbit.coeffs[1].len() {
        assert!(
            (orbit.coeffs[2][n] - C64::new(0.0, n as f64 * orbit.omega) * orbit.coeffs[1][n]).abs()
                < 2e-10
        );
    }
}

#[test]
fn tiny_physical_motion_large_static_preload_and_independent_units_keep_absolute_resolution() {
    let original = Polynomial::new().solve();
    let mut converted = Polynomial::new();
    converted.offset = 4e12;
    converted.state = [1e6, 1e-11, 1e8];
    converted.equation = [1e-7, 1e5, 3.0];
    let result = converted.solve();
    assert_eq!(
        result.refinement_trace.len(),
        original.refinement_trace.len()
    );
    assert!((result.point.orbit.omega - original.point.orbit.omega).abs() < 2e-10);
    for state in 1..3 {
        for harmonic in 0..result.point.orbit.coeffs[state].len() {
            assert!(
                (result.point.orbit.coeffs[state][harmonic].scale(1.0 / converted.state[state])
                    - original.point.orbit.coeffs[state][harmonic])
                    .abs()
                    < 2e-10
            );
        }
    }
    let failure = solve_backbone_adaptive(
        &converted,
        1,
        0.42 * converted.state[1],
        &converted.seed(),
        &controls(),
        &converted.refinement(3),
        &converted.scales(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(failure, BackboneContinuationError::Numerical(
        OrbitError::HarmonicsUnresolved { check }
    ) if check.harmonics == 3 && check.omitted_residual_ratio > 1.0));
}

#[test]
fn global_admission_cancellation_and_retained_accuracy_cannot_be_bypassed_by_order_growth() {
    let model = Polynomial::new();
    let mut small = controls();
    small.maximum_unknowns = 20;
    assert!(
        solve_backbone_adaptive(
            &model,
            1,
            0.42,
            &model.seed(),
            &small,
            &model.refinement(31),
            &model.scales(),
            || false
        )
        .is_err()
    );
    assert_eq!(model.calls.get(), 0);
    small = controls();
    small.maximum_payload_bytes = 16;
    assert!(
        solve_backbone_adaptive(
            &model,
            1,
            0.42,
            &model.seed(),
            &small,
            &model.refinement(31),
            &model.scales(),
            || false
        )
        .is_err()
    );
    assert_eq!(model.calls.get(), 0);
    small = controls();
    small.maximum_residual_evaluations = 1;
    assert!(matches!(
        solve_backbone_adaptive(
            &model,
            1,
            0.42,
            &model.seed(),
            &small,
            &model.refinement(31),
            &model.scales(),
            || false
        ),
        Err(BackboneContinuationError::EvaluationLimit { evaluations: 1 })
    ));
    let before = model.calls.get();
    assert!(matches!(
        solve_backbone_adaptive(
            &model,
            1,
            0.42,
            &model.seed(),
            &controls(),
            &model.refinement(31),
            &model.scales(),
            || model.calls.get() > before + 200
        ),
        Err(BackboneContinuationError::Cancelled)
    ));
    let mut inaccurate = model.seed();
    inaccurate.omega *= 1.001;
    small = controls();
    small.tolerance = 0.1;
    assert!(
        matches!(solve_backbone_adaptive(&model, 1, 0.42, &inaccurate, &small,
        &model.refinement(31), &model.scales(), || false),
        Err(BackboneContinuationError::Numerical(OrbitError::HarmonicsUnresolved { check }))
        if check.harmonics == 1 && check.retained_residual_ratio > 1.0)
    );
    let mut refused = Polynomial::new();
    refused.refused = true;
    assert!(matches!(
        solve_backbone_adaptive(
            &refused,
            1,
            0.42,
            &refused.seed(),
            &controls(),
            &refused.refinement(31),
            &refused.scales(),
            || false
        ),
        Err(BackboneContinuationError::Numerical(
            OrbitError::TangentRefused { .. }
        ))
    ));
}

struct AdditionalFamily;
impl OrbitProblem for AdditionalFamily {
    fn dim(&self) -> usize {
        4
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        out.copy_from_slice(&[x[1], -4.0 * x[0], x[3], -4.0 * x[2]]);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|_, _, _, out| {
            out.fill(0.0);
            out[1] = 1.0;
            out[4] = -4.0;
            out[11] = 1.0;
            out[14] = -4.0;
            true
        })
    }
}
#[test]
fn exact_residual_seed_with_an_additional_neutral_family_is_not_regularized_into_acceptance() {
    let seed = HbOrbit {
        coeffs: vec![
            vec![C64::ZERO, C64::from_re(0.2)],
            vec![C64::ZERO, C64::new(0.0, 0.4)],
            vec![C64::ZERO; 2],
            vec![C64::ZERO; 2],
        ],
        omega: 2.0,
        residual: 0.0,
        residual_trace: vec![],
        aft_samples: 0,
        sampling_defect: 0.0,
    };
    let failure = solve_backbone_adaptive(
        &AdditionalFamily,
        0,
        0.4,
        &seed,
        &controls(),
        &HbRefinementBudget {
            max_harmonics: 3,
            residual_tolerances: vec![1e-8; 4],
            state_tolerances: vec![1e-8; 4],
            frequency_tolerance: 1e-8,
        },
        &HbScales {
            states: vec![1.0; 4],
            equations: vec![1.0; 4],
            angular_frequency: 2.0,
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        failure,
        BackboneContinuationError::Numerical(OrbitError::SingularJacobian)
            | BackboneContinuationError::IllConditioned { .. }
    ));
}
