//! G1/G3: physical cycle means include the changing autonomous period.

use fs_math::{c64::C64, det};
use fs_orbit::sensitivity::{
    AutonomousSensitivityPoint, ForcedMeanObjectiveControls, ForcedSensitivityControls,
    ForcedSensitivityError, PreparedAutonomousSensitivity,
};
use fs_orbit::{HbAnchor, HbOrbit, HbScales, IslandStateJacobian, OrbitProblem};
use std::cell::Cell;
use std::f64::consts::TAU;

// An independently known nonisochronous attracting cycle: r=sqrt(alpha),
// omega=1.3+0.4*alpha. Unit changes act on states, equations, time, parameters,
// and the observation independently. No harmonic solve supplies this oracle.
struct Cycle {
    alpha: f64,
    units: [f64; 2],
    equations: [f64; 2],
    clock: f64,
    parameter_unit: f64,
    objective_unit: f64,
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}

impl Cycle {
    fn new(converted: bool) -> Self {
        Self {
            alpha: 0.8,
            units: if converted { [1e5, 1e-7] } else { [1.0; 2] },
            equations: if converted { [1e-4, 1e6] } else { [1.0; 2] },
            clock: if converted { 1e3 } else { 1.0 },
            parameter_unit: if converted { 1e-6 } else { 1.0 },
            objective_unit: if converted { 1e9 } else { 1.0 },
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }

    fn orbit(&self, phase: f64) -> HbOrbit {
        let z = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * det::sqrt(self.alpha));
        let mut coeffs = vec![vec![C64::ZERO; 4]; 2];
        coeffs[0][1] = z.scale(self.units[0]);
        coeffs[1][1] = (C64::new(0.0, -1.0) * z).scale(self.units[1]);
        HbOrbit {
            omega: self.clock * (1.3 + 0.4 * self.alpha),
            coeffs,
            residual_trace: vec![],
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }

    fn prepare(&self, orbit: &HbOrbit) -> PreparedAutonomousSensitivity {
        let reference = orbit.packed_seed(&HbAnchor::Autonomous {
            omega_guess: orbit.omega,
        });
        PreparedAutonomousSensitivity::new(
            self,
            orbit,
            &reference,
            ForcedSensitivityControls {
                samples_per_period: 32,
                maximum_unknowns: 256,
                maximum_payload_bytes: 16 * 1024 * 1024,
                maximum_work: 500_000_000,
                maximum_base_residual: 1e-10,
                maximum_tangent_sampling_defect: 1e-10,
                maximum_condition_estimate: 1e9,
                maximum_relative_residual: 1e-12,
                maximum_derivative_residual: 1e-9,
            },
            &HbScales {
                states: self.units.to_vec(),
                equations: self.equations.iter().map(|e| e * self.clock).collect(),
                angular_frequency: 2.0 * self.clock,
            },
            |point, out| self.frequency_column(point, out),
            || false,
        )
        .unwrap()
    }

    fn controls(&self) -> ForcedMeanObjectiveControls {
        ForcedMeanObjectiveControls {
            objective_scale: self.objective_unit,
            maximum_value_sampling_defect: 1e-10,
            maximum_gradient_sampling_defect: 1e-10,
        }
    }

    fn frequency_column(&self, point: AutonomousSensitivityPoint<'_>, out: &mut [f64]) -> bool {
        out.fill(0.0);
        for n in 1..=point.harmonics {
            for j in 0..2 {
                let re = 2 + (n - 1) * 4 + j;
                let im = re + 2;
                let scale = n as f64 * self.equations[j] / self.units[j];
                out[re] = -scale * point.coefficients[im];
                out[im] = scale * point.coefficients[re];
            }
        }
        true
    }

    fn tangent(&self, _: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        let (x, y) = (state[0] / self.units[0], state[1] / self.units[1]);
        let r2 = x * x + y * y;
        let radial = self.alpha - r2;
        let w = 1.3 + 0.4 * r2;
        out.copy_from_slice(&[
            radial - 2.0 * x * x - 0.8 * x * y,
            -2.0 * x * y - w - 0.8 * y * y,
            w + 0.8 * x * x - 2.0 * x * y,
            radial + 0.8 * x * y - 2.0 * y * y,
        ]);
        for row in 0..2 {
            for col in 0..2 {
                out[2 * row + col] *= self.clock * self.equations[row] / self.units[col];
            }
        }
        true
    }

    fn energy(
        &self,
        omega: f64,
        state: &[f64],
        value: &mut f64,
        gradient: &mut [f64],
        frequency: &mut f64,
    ) -> bool {
        let (x, y) = (state[0] / self.units[0], state[1] / self.units[1]);
        let period = TAU * self.clock / omega;
        let scale = period * self.objective_unit;
        *value = scale * (0.4 * (x * x + y * y) + 0.2 * x.powi(4));
        gradient[0] = scale * (0.8 * x + 0.8 * x.powi(3)) / self.units[0];
        gradient[1] = scale * 0.8 * y / self.units[1];
        *frequency = -*value / omega;
        true
    }
}

impl OrbitProblem for Cycle {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, state: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        let (x, y) = (state[0] / self.units[0], state[1] / self.units[1]);
        let r2 = x * x + y * y;
        let w = 1.3 + 0.4 * r2;
        out[0] = self.clock * self.equations[0] * ((self.alpha - r2) * x - w * y);
        out[1] = self.clock * self.equations[1] * (w * x + (self.alpha - r2) * y);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![
            s.scale(self.equations[0] / self.units[0]),
            C64::ZERO,
            C64::ZERO,
            s.scale(self.equations[1] / self.units[1]),
        ]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::tangent)
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn energy_per_cycle_includes_waveform_period_and_independent_units() {
    for converted in [false, true] {
        for phase in [0.0, 0.47, -1.2] {
            let problem = Cycle::new(converted);
            let orbit = problem.orbit(phase);
            let mut prepared = problem.prepare(&orbit);
            let calls = (problem.source_calls.get(), problem.tangent_calls.get());
            let objective = prepared
                .mean_objective(
                    problem.controls(),
                    |_, x, value, gradient, frequency| {
                        problem.energy(orbit.omega, x, value, gradient, frequency)
                    },
                    || false,
                )
                .unwrap();
            let a = problem.alpha;
            let w = 1.3 + 0.4 * a;
            let power = 0.4 * a + 0.075 * a * a;
            let energy = TAU * power / w;
            close(objective.value / problem.objective_unit, energy);
            close(
                objective.angular_frequency_partial * problem.clock / problem.objective_unit,
                -energy / w,
            );
            // The cubic local derivative contributes a third-harmonic covector,
            // even though this exact base cycle has only a fundamental.
            assert!(objective.gradient[10..14].iter().any(|v| v.abs() > 1e-5));
            let mut adjoint = prepared
                .prepare_adjoint(
                    &objective.gradient,
                    objective.angular_frequency_partial,
                    problem.objective_unit,
                    || false,
                )
                .unwrap();
            let result = adjoint
                .gradient_island_direction(
                    problem.parameter_unit,
                    0.0,
                    |_, x, out| {
                        for j in 0..2 {
                            out[j] = problem.clock * problem.equations[j] * x[j]
                                / problem.units[j]
                                / problem.parameter_unit;
                        }
                        true
                    },
                    || false,
                )
                .unwrap();
            let expected = TAU * ((0.4 + 0.15 * a) / w - 0.4 * power / (w * w));
            close(
                result.derivative * problem.parameter_unit / problem.objective_unit,
                expected,
            );
            assert_eq!(
                calls,
                (problem.source_calls.get(), problem.tangent_calls.get())
            );
        }
    }
}

#[test]
fn state_only_mean_has_exact_parseval_covector_and_no_clock_partial() {
    let problem = Cycle::new(false);
    let orbit = problem.orbit(0.63);
    let mut prepared = problem.prepare(&orbit);
    let mut calls = 0;
    let observation = prepared
        .mean_objective(
            problem.controls(),
            |_, x, value, gradient, frequency| {
                calls += 1;
                *value = x[0] * x[0] + x[1] * x[1];
                gradient.copy_from_slice(&[2.0 * x[0], 2.0 * x[1]]);
                *frequency = 0.0;
                true
            },
            || false,
        )
        .unwrap();
    assert_eq!(calls, 32 + 65);
    close(observation.value, problem.alpha);
    assert_eq!(observation.angular_frequency_partial, 0.0);
    for j in 0..2 {
        close(observation.gradient[j], 2.0 * orbit.coeffs[j][0].re);
        for n in 1..=3 {
            let re = 2 + (n - 1) * 4 + j;
            close(observation.gradient[re], 4.0 * orbit.coeffs[j][n].re);
            close(observation.gradient[re + 2], 4.0 * orbit.coeffs[j][n].im);
        }
    }
}

#[test]
fn unresolved_frequency_covector_refuses_even_when_value_and_state_gradients_agree() {
    let problem = Cycle::new(false);
    let orbit = problem.orbit(0.2);
    let mut prepared = problem.prepare(&orbit);
    // Local jet of (omega-omega_base)*cos(32*phase): value and state
    // derivatives vanish at the base, while its frequency derivative aliases.
    let error = prepared
        .mean_objective(
            problem.controls(),
            |time, _, value, gradient, frequency| {
                *value = 0.0;
                gradient.fill(0.0);
                *frequency = det::cos(32.0 * orbit.omega * time);
                true
            },
            || false,
        )
        .unwrap_err();
    match error {
        ForcedSensitivityError::ObjectiveSamplingUnresolved {
            value_defect,
            gradient_defect,
            ..
        } => {
            assert_eq!(value_defect, 0.0);
            assert!(gradient_defect > 1.0);
        }
        other => panic!("unexpected refusal: {other:?}"),
    }
}

#[test]
fn complete_clock_declaration_and_cancellation_preserve_cached_reuse() {
    let problem = Cycle::new(false);
    let orbit = problem.orbit(0.1);
    let mut prepared = problem.prepare(&orbit);
    let calls = (problem.source_calls.get(), problem.tangent_calls.get());
    let missing = prepared
        .mean_objective(
            problem.controls(),
            |_, _, value, gradient, _| {
                *value = 0.0;
                gradient.fill(0.0);
                true
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        missing,
        ForcedSensitivityError::ObjectiveRefused { sample: 0, .. }
    ));
    let mut polls = 0;
    let cancelled = prepared
        .mean_objective(
            problem.controls(),
            |_, x, value, gradient, frequency| {
                problem.energy(orbit.omega, x, value, gradient, frequency)
            },
            || {
                polls += 1;
                polls > 4
            },
        )
        .unwrap_err();
    assert_eq!(cancelled, ForcedSensitivityError::Cancelled);
    let first = prepared
        .mean_objective(
            problem.controls(),
            |_, x, value, gradient, frequency| {
                problem.energy(orbit.omega, x, value, gradient, frequency)
            },
            || false,
        )
        .unwrap();
    let second = prepared
        .mean_objective(
            problem.controls(),
            |_, x, value, gradient, frequency| {
                problem.energy(orbit.omega, x, value, gradient, frequency)
            },
            || false,
        )
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(
        calls,
        (problem.source_calls.get(), problem.tangent_calls.get())
    );
}
