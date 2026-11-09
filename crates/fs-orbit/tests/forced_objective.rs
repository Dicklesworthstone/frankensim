//! G1/G3 cycle-mean observations, real Fourier covectors and physical adjoints.

use fs_math::{c64::C64, det};
use fs_orbit::{
    HbAnchor, HbBudget, HbOrbit, HbScales, IslandStateJacobian, OrbitProblem,
    sensitivity::{
        ForcedMeanObjectiveControls, ForcedSensitivityControls, ForcedSensitivityDirection,
        ForcedSensitivityError, PreparedForcedSensitivity,
    },
    solve_hb_scaled_seeded,
};
use std::cell::Cell;

const OMEGA: f64 = 0.9;

#[path = "forced_objective/mean_partial.rs"]
mod mean_partial;

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
fn objective_controls(objective_scale: f64) -> ForcedMeanObjectiveControls {
    ForcedMeanObjectiveControls {
        objective_scale,
        maximum_value_sampling_defect: 1e-11,
        maximum_gradient_sampling_defect: 1e-11,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() < tolerance * (1.0 + expected.abs()),
        "{actual:.17e} vs {expected:.17e}"
    );
}

// An exact, two-coordinate linear compliance driven in every retained complex
// quadrature. Independent state/equation units exercise the original HB chart.
struct Compliance {
    orbit: HbOrbit,
    states: [f64; 2],
    equations: [f64; 2],
    calls: [Cell<usize>; 3],
}
impl Compliance {
    fn new() -> Self {
        Self {
            orbit: HbOrbit {
                omega: OMEGA,
                coeffs: vec![
                    vec![
                        C64::from_re(0.3),
                        C64::new(0.2, -0.13),
                        C64::new(0.08, 0.11),
                        C64::new(0.01, -0.02),
                    ],
                    vec![
                        C64::from_re(-0.15),
                        C64::new(-0.15, 0.07),
                        C64::new(0.04, -0.05),
                        C64::new(-0.03, 0.015),
                    ],
                ],
                residual_trace: vec![],
                residual: 0.0,
                aft_samples: 0,
                sampling_defect: 0.0,
            },
            states: [1.0; 2],
            equations: [1.0; 2],
            calls: std::array::from_fn(|_| Cell::new(0)),
        }
    }
    fn counts(&self) -> [usize; 3] {
        std::array::from_fn(|i| self.calls[i].get())
    }
    fn prepare(&self, samples: usize) -> PreparedForcedSensitivity {
        let mut orbit = self.orbit.clone();
        for (row, unit) in orbit.coeffs.iter_mut().zip(self.states) {
            for coefficient in row {
                *coefficient = coefficient.scale(unit);
            }
        }
        PreparedForcedSensitivity::new(
            self,
            &orbit,
            controls(samples),
            &HbScales {
                states: self.states.to_vec(),
                equations: self.equations.to_vec(),
                angular_frequency: 1.0,
            },
            || false,
        )
        .unwrap()
    }
    fn square(&self, x: &[f64], scale: f64, value: &mut f64, gradient: &mut [f64]) -> bool {
        let observed = x[0] / self.states[0] + 0.7 * x[1] / self.states[1];
        *value = scale * observed * observed;
        gradient[0] = 2.0 * scale * observed / self.states[0];
        gradient[1] = 1.4 * scale * observed / self.states[1];
        true
    }
}
impl OrbitProblem for Compliance {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, _: &[f64], out: &mut [f64]) {
        self.calls[0].set(self.calls[0].get() + 1);
        for (i, value) in out.iter_mut().enumerate() {
            *value = self.equations[i] * self.orbit.sample(i, OMEGA * t);
        }
    }
    fn port(&self, _: C64) -> Vec<C64> {
        self.calls[1].set(self.calls[1].get() + 1);
        vec![
            C64::from_re(self.equations[0] / self.states[0]),
            C64::ZERO,
            C64::ZERO,
            C64::from_re(self.equations[1] / self.states[1]),
        ]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|p, _, _, out| {
            p.calls[2].set(p.calls[2].get() + 1);
            out.fill(0.0);
            true
        })
    }
}

#[test]
fn mean_square_matches_parseval_and_every_real_coordinate_without_new_source_calls() {
    let problem = Compliance::new();
    let mut prepared = problem.prepare(32);
    let baseline = problem.counts();
    let calls = Cell::new(0);
    let result = prepared
        .mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| {
                calls.set(calls.get() + 1);
                problem.square(x, 1.0, value, gradient)
            },
            || false,
        )
        .unwrap();
    assert_eq!(calls.get(), 32 + 65);
    assert_eq!(problem.counts(), baseline);
    let mut expected = 0.0;
    for n in 0..=3 {
        let y = problem.orbit.coeffs[0][n] + problem.orbit.coeffs[1][n].scale(0.7);
        expected += if n == 0 {
            y.re * y.re
        } else {
            2.0 * (y.re * y.re + y.im * y.im)
        };
        for (state, weight) in [1.0, 0.7].into_iter().enumerate() {
            if n == 0 {
                close(result.gradient[state], 2.0 * weight * y.re, 2e-14);
            } else {
                close(
                    result.gradient[(2 * n - 1) * 2 + state],
                    4.0 * weight * y.re,
                    2e-14,
                );
                close(
                    result.gradient[2 * n * 2 + state],
                    4.0 * weight * y.im,
                    2e-14,
                );
            }
        }
    }
    close(result.value, expected, 2e-14);
    let replay = prepared
        .mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| problem.square(x, 1.0, value, gradient),
            || false,
        )
        .unwrap();
    assert_eq!(result, replay);
}

#[test]
fn objective_and_state_unit_changes_preserve_the_physical_covector() {
    let original = Compliance::new();
    let reference = original
        .prepare(32)
        .mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| original.square(x, 1.0, value, gradient),
            || false,
        )
        .unwrap();
    let mut converted = Compliance::new();
    converted.states = [1e-5, 1e4];
    converted.equations = [1e2, 1e-3];
    let units = 1e6;
    let result = converted
        .prepare(32)
        .mean_objective(
            objective_controls(units),
            |_, x, value, gradient| converted.square(x, units, value, gradient),
            || false,
        )
        .unwrap();
    close(result.value / units, reference.value, 3e-14);
    for (i, value) in result.gradient.iter().enumerate() {
        close(
            value * converted.states[i % 2] / units,
            reference.gradient[i],
            3e-14,
        );
    }
}

struct NonlinearDamper {
    damping: f64,
    force: f64,
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}
impl NonlinearDamper {
    fn new(damping: f64, force: f64) -> Self {
        Self {
            damping,
            force,
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }
    fn solve(&self, harmonics: usize) -> HbOrbit {
        solve_hb_scaled_seeded(
            self,
            HbAnchor::Forced { omega: OMEGA },
            vec![0.0; 2 * (2 * harmonics + 1)],
            &HbBudget {
                harmonics,
                max_newton: 40,
                tolerance: 1e-12,
                deflate_equilibrium: false,
            },
            &unit_scales(),
        )
        .unwrap()
        .orbit
    }
}
impl OrbitProblem for NonlinearDamper {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        out[0] = x[1];
        out[1] = -1.7 * x[0] - 0.3 * x[0].powi(3) - 0.4 * x[1] - self.damping * x[1].powi(3)
            + 0.15
            + self.force * det::cos(OMEGA * t);
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|p, _, x, out| {
            p.tangent_calls.set(p.tangent_calls.get() + 1);
            out.copy_from_slice(&[
                0.0,
                1.0,
                -1.7 - 0.9 * x[0] * x[0],
                -0.4 - 3.0 * p.damping * x[1] * x[1],
            ]);
            true
        })
    }
}
fn unit_scales() -> HbScales {
    HbScales {
        states: vec![1.0; 2],
        equations: vec![1.0; 2],
        angular_frequency: 1.0,
    }
}
fn resolved_mean(orbit: &HbOrbit, function: impl Fn(f64) -> f64) -> f64 {
    (0..1024)
        .map(|i| function(orbit.sample(1, core::f64::consts::TAU * i as f64 / 1024.0)))
        .sum::<f64>()
        / 1024.0
}

#[test]
fn nonlinear_dissipated_power_adjoints_match_independently_resolved_parameter_changes() {
    let problem = NonlinearDamper::new(0.18, 0.8);
    let orbit = problem.solve(11);
    assert!(orbit.coeffs[1][3].abs() > 1e-3);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &unit_scales(), || false)
            .unwrap();
    let baseline = (problem.source_calls.get(), problem.tangent_calls.get());
    let objective = prepared
        .mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| {
                *value = 0.4 * x[1].powi(2) + problem.damping * x[1].powi(4);
                gradient.copy_from_slice(&[0.0, 0.8 * x[1] + 4.0 * problem.damping * x[1].powi(3)]);
                true
            },
            || false,
        )
        .unwrap();
    close(
        objective.value,
        resolved_mean(&orbit, |v| 0.4 * v * v + problem.damping * v.powi(4)),
        2e-13,
    );
    let mut adjoint = prepared
        .prepare_adjoint(&objective.gradient, 1.0, || false)
        .unwrap();
    let direction = ForcedSensitivityDirection {
        parameter: problem.damping,
        parameter_scale: 0.1,
        angular_frequency_derivative: 0.0,
    };
    // The dissipation law itself changes with nonlinear damping, independently
    // of the periodic state; its explicit <v^4> partial must be added once.
    let material = adjoint
        .gradient_island_direction(
            direction,
            resolved_mean(&orbit, |v| v.powi(4)),
            |_, x, out| {
                out.copy_from_slice(&[0.0, -x[1].powi(3)]);
                true
            },
            || false,
        )
        .unwrap();
    let force = adjoint
        .gradient_island_direction(
            ForcedSensitivityDirection {
                parameter: problem.force,
                parameter_scale: 1.0,
                ..direction
            },
            0.0,
            |t, _, out| {
                out.copy_from_slice(&[0.0, det::cos(OMEGA * t)]);
                true
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        baseline,
        (problem.source_calls.get(), problem.tangent_calls.get())
    );
    let step = 2e-5;
    for (actual, low, high) in [
        (
            material.derivative,
            NonlinearDamper::new(0.18 - step, 0.8),
            NonlinearDamper::new(0.18 + step, 0.8),
        ),
        (
            force.derivative,
            NonlinearDamper::new(0.18, 0.8 - step),
            NonlinearDamper::new(0.18, 0.8 + step),
        ),
    ] {
        let independent = |p: &NonlinearDamper| {
            resolved_mean(&p.solve(15), |v| 0.4 * v * v + p.damping * v.powi(4))
        };
        close(
            actual,
            (independent(&high) - independent(&low)) / (2.0 * step),
            2e-8,
        );
    }
}

#[test]
fn equal_objective_values_cannot_hide_an_aliased_gradient_or_the_reverse() {
    let mut problem = Compliance::new();
    for row in &mut problem.orbit.coeffs {
        row.fill(C64::ZERO);
    }
    let mut prepared = problem.prepare(32);
    for gradient_alias in [true, false] {
        let refusal = prepared.mean_objective(
            objective_controls(1.0),
            |time, x, value, gradient| {
                let alias = det::cos(32.0 * OMEGA * time);
                *value = if gradient_alias { alias * x[0] } else { alias };
                gradient.copy_from_slice(&[if gradient_alias { alias } else { 0.0 }, 0.0]);
                true
            },
            || false,
        );
        let Err(ForcedSensitivityError::ObjectiveSamplingUnresolved {
            value_defect,
            gradient_defect,
            ..
        }) = refusal
        else {
            panic!("unresolved objective accepted: {refusal:?}");
        };
        if gradient_alias {
            assert_eq!(value_defect, 0.0);
            assert!(gradient_defect > 0.9);
        } else {
            assert!(value_defect > 0.9);
            assert_eq!(gradient_defect, 0.0);
        }
    }
}

#[test]
fn incomplete_callbacks_and_cancellation_refuse_then_reuse_the_unchanged_waveform() {
    let problem = Compliance::new();
    let mut prepared = problem.prepare(32);
    let baseline = problem.counts();
    for bad in [
        ForcedMeanObjectiveControls {
            objective_scale: 0.0,
            ..objective_controls(1.0)
        },
        ForcedMeanObjectiveControls {
            maximum_value_sampling_defect: f64::NAN,
            ..objective_controls(1.0)
        },
        ForcedMeanObjectiveControls {
            maximum_gradient_sampling_defect: -1.0,
            ..objective_controls(1.0)
        },
    ] {
        assert!(matches!(
            prepared.mean_objective(bad, |_, _, _, _| panic!("invalid controls"), || false),
            Err(ForcedSensitivityError::Invalid(_))
        ));
    }
    assert!(matches!(
        prepared.mean_objective(
            objective_controls(1.0),
            |_, _, _, gradient| {
                gradient.fill(0.0);
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::ObjectiveRefused { sample: 0, .. })
    ));
    let calls = Cell::new(0);
    assert!(matches!(
        prepared.mean_objective(
            objective_controls(1.0),
            |_, _, value, gradient| {
                *value = 0.0;
                if calls.get() == 0 {
                    gradient.fill(0.0);
                } else {
                    gradient[0] = 0.0;
                }
                calls.set(calls.get() + 1);
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::ObjectiveRefused { sample: 1, .. })
    ));
    assert!(matches!(
        prepared.mean_objective(
            objective_controls(1.0),
            |_, _, _, _| panic!("cancelled before callback"),
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    calls.set(0);
    assert!(matches!(
        prepared.mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| {
                calls.set(calls.get() + 1);
                problem.square(x, 1.0, value, gradient)
            },
            || calls.get() == 3
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert_eq!(calls.get(), 3);
    let recovered = prepared
        .mean_objective(
            objective_controls(1.0),
            |_, x, value, gradient| problem.square(x, 1.0, value, gradient),
            || false,
        )
        .unwrap();
    let mut adjoint = prepared
        .prepare_adjoint(&recovered.gradient, 1.0, || false)
        .unwrap();
    let derivative = adjoint
        .gradient_island_direction(
            ForcedSensitivityDirection {
                parameter: 0.0,
                parameter_scale: 1.0,
                angular_frequency_derivative: 0.0,
            },
            0.0,
            |_, _, out| {
                out.fill(0.0);
                true
            },
            || false,
        )
        .unwrap();
    assert_eq!(derivative.derivative, 0.0);
    assert_eq!(problem.counts(), baseline);
}
