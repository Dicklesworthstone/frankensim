//! G1/G3: self-excited frequency, full memory waveform and a fixed phase gauge.

use fs_math::{c64::C64, det};
use fs_orbit::sensitivity::{
    AutonomousOrbitSensitivity, AutonomousSensitivityPoint, ForcedSensitivityControls,
    ForcedSensitivityError, PreparedAutonomousSensitivity,
};
use fs_orbit::{HbAnchor, HbOrbit, HbScales, IslandStateJacobian, OrbitError, OrbitProblem};
use std::cell::Cell;

// The radius satisfies a=p, and the two lag states feed their squared magnitude
// into angular speed: W=1.2+.3*a+.4*a^2*g^2/(g^2+4*W^2). Their quadratic
// forcing produces a SECOND harmonic with a distinct physical lag. The mass/time
// factor tau=1+.15*p also changes with the parameter. This scalar cycle oracle
// is independent of harmonic balance and of its exact state Jacobian.
struct MemoryCycle {
    p: f64,
    parameter_unit: f64,
    states: [f64; 5],
    equations: [f64; 5],
    gain: f64,
    island_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}
impl MemoryCycle {
    fn new(converted: bool) -> Self {
        Self {
            p: 0.8,
            parameter_unit: if converted { 1e-4 } else { 1.0 },
            states: if converted {
                [3e4, 1e-9, 2e3, 4e-7, 1e2]
            } else {
                [1.0; 5]
            },
            equations: if converted {
                [7.0, 1e2, 5e-6, 3e-8, 2e4]
            } else {
                [1.0; 5]
            },
            gain: 1.0,
            island_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }
    fn scales(&self) -> HbScales {
        HbScales {
            states: self.states.to_vec(),
            equations: self.equations.to_vec(),
            angular_frequency: 2.0,
        }
    }
    fn internal_frequency(&self) -> f64 {
        let (mut lo, mut hi) = (
            1.2 + 0.3 * self.p,
            1.2 + 0.3 * self.p + 0.4 * self.p * self.p,
        );
        for _ in 0..70 {
            let w = 0.5 * (lo + hi);
            if w - 1.2 - 0.3 * self.p - 0.4 * self.p * self.p * 0.49 / (0.49 + 4.0 * w * w) > 0.0 {
                hi = w;
            } else {
                lo = w;
            }
        }
        0.5 * (lo + hi)
    }
    fn orbit(&self, phase: f64) -> HbOrbit {
        let w = self.internal_frequency();
        let x = C64::new(det::cos(phase), det::sin(phase)).scale(0.5 * det::sqrt(self.p));
        let y = C64::new(0.0, -1.0) * x;
        let lag = C64::from_re(0.7) / C64::new(0.7, 2.0 * w);
        let drive = (x * x).scale(2.0);
        let mut coeffs = vec![vec![C64::ZERO; 4]; 5];
        coeffs[0][0] = C64::from_re((0.3 + 0.07 * self.p) * self.states[0]);
        for (j, value) in [x, y, lag * drive, C64::new(0.0, -1.0) * lag * drive]
            .into_iter()
            .enumerate()
        {
            coeffs[j + 1][if j < 2 { 1 } else { 2 }] = value.scale(self.states[j + 1]);
        }
        HbOrbit {
            omega: w / (1.0 + 0.15 * self.p),
            coeffs,
            residual_trace: vec![],
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }
    fn frequency_column(&self, point: AutonomousSensitivityPoint<'_>, out: &mut [f64]) -> bool {
        out.fill(0.0);
        for k in 1..=point.harmonics {
            for j in 0..5 {
                let re = 5 + (k - 1) * 10 + j;
                let im = re + 5;
                let slope = k as f64 * (1.0 + 0.15 * self.p) * self.equations[j] / self.states[j];
                out[re] = -slope * point.coefficients[im];
                out[im] = slope * point.coefficients[re];
            }
        }
        true
    }
    fn parameter_column(&self, point: AutonomousSensitivityPoint<'_>, out: &mut [f64]) -> bool {
        out.fill(0.0);
        out[0] = -0.07 * self.equations[0] / self.parameter_unit;
        for k in 0..=point.harmonics {
            for j in 0..5 {
                let re = if k == 0 { j } else { 5 + (k - 1) * 10 + j };
                let unit = self.equations[j] / self.states[j] / self.parameter_unit;
                if k > 0 {
                    let im = re + 5;
                    let clock = 0.15 * k as f64 * point.angular_frequency * unit;
                    out[re] -= clock * point.coefficients[im];
                    out[im] += clock * point.coefficients[re];
                    if j == 1 || j == 2 {
                        out[im] -= self.gain * unit * point.coefficients[im];
                    }
                }
                if j == 1 || j == 2 {
                    out[re] -= self.gain * unit * point.coefficients[re];
                }
            }
        }
        true
    }
    fn exact_tangent(&self, _: f64, physical: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        let x: [f64; 5] = std::array::from_fn(|j| physical[j] / self.states[j]);
        let a = x[1] * x[1] + x[2] * x[2];
        let b = x[3] * x[3] + x[4] * x[4];
        let radial = self.gain * (self.p - a);
        let w = 1.2 + 0.3 * a + 0.4 * b;
        out.fill(0.0);
        out[0] = -1.0;
        out[6] = radial - 2.0 * self.gain * x[1] * x[1] - 0.6 * x[1] * x[2];
        out[7] = -2.0 * self.gain * x[1] * x[2] - w - 0.6 * x[2] * x[2];
        out[11] = w + 0.6 * x[1] * x[1] - 2.0 * self.gain * x[1] * x[2];
        out[12] = radial - 2.0 * self.gain * x[2] * x[2] + 0.6 * x[1] * x[2];
        for j in [3, 4] {
            out[5 + j] = -0.8 * x[j] * x[2];
            out[10 + j] = 0.8 * x[j] * x[1];
        }
        out[16] = 1.4 * x[1];
        out[17] = -1.4 * x[2];
        out[18] = -0.7;
        out[21] = 1.4 * x[2];
        out[22] = 1.4 * x[1];
        out[24] = -0.7;
        for row in 0..5 {
            for col in 0..5 {
                out[row * 5 + col] *= self.equations[row] / self.states[col];
            }
        }
        true
    }
    fn verify(&self, orbit: &HbOrbit, result: &AutonomousOrbitSensitivity, material: bool) {
        let w = self.internal_frequency();
        let denominator = 0.49 + 4.0 * w * w;
        let correction = 1.0 + 3.2 * self.p * self.p * 0.49 * w / (denominator * denominator);
        let dw = if material {
            (0.3 + 0.8 * self.p * 0.49 / denominator) / correction
        } else {
            1.0 / correction
        };
        let tau = 1.0 + 0.15 * self.p;
        let omega_derivative = if material {
            (dw * tau - 0.15 * w) / (tau * tau)
        } else {
            dw / tau
        };
        assert!(
            (result.angular_frequency_derivative * self.parameter_unit - omega_derivative).abs()
                < 2e-11
        );
        let x = orbit.coeffs[1][1].scale(1.0 / self.states[1]);
        let y = orbit.coeffs[2][1].scale(1.0 / self.states[2]);
        let lag = C64::from_re(0.7) / C64::new(0.7, 2.0 * w);
        let lag_derivative =
            C64::new(0.0, -1.4 * dw) / (C64::new(0.7, 2.0 * w) * C64::new(0.7, 2.0 * w));
        let radial = if material { 0.5 / self.p } else { 0.0 };
        let drive = (x * x).scale(2.0);
        let drive_y = C64::new(0.0, -1.0) * drive;
        let base = [x, y, lag * drive, lag * drive_y];
        let mut derivatives = [
            x.scale(radial),
            y.scale(radial),
            lag_derivative * drive + (lag * drive).scale(2.0 * radial),
            lag_derivative * drive_y + (lag * drive_y).scale(2.0 * radial),
        ];
        let mut dot = 0.0;
        let mut norm = 0.0;
        for (j, (z, dz)) in base.iter().zip(&derivatives).enumerate() {
            let harmonic = if j < 2 { 1.0 } else { 2.0 };
            dot += harmonic * (-z.im * dz.re + z.re * dz.im);
            norm += harmonic * harmonic * (z.re * z.re + z.im * z.im);
        }
        let phase_derivative = -dot / norm;
        for j in 0..4 {
            derivatives[j] = derivatives[j]
                + C64::new(0.0, phase_derivative * if j < 2 { 1.0 } else { 2.0 }) * base[j];
        }
        assert!(
            (result.coeffs[0][0].re * self.parameter_unit / self.states[0]
                - if material { 0.07 } else { 0.0 })
            .abs()
                < 2e-11
        );
        for j in 1..5 {
            let harmonic = if j < 3 { 1 } else { 2 };
            let actual = result.coeffs[j][harmonic].scale(self.parameter_unit / self.states[j]);
            assert!(
                (actual - derivatives[j - 1]).abs() < 3e-11,
                "state {j}: {actual:?} expected {:?}",
                derivatives[j - 1]
            );
            assert!(result.coeffs[j][0].abs() * self.parameter_unit / self.states[j] < 2e-11);
            assert!(
                result.coeffs[j]
                    .iter()
                    .enumerate()
                    .all(|(k, z)| k == harmonic
                        || z.abs() * self.parameter_unit / self.states[j] < 2e-11)
            );
        }
        assert!(result.derivative_residual < 1e-10);
    }
}
impl OrbitProblem for MemoryCycle {
    fn dim(&self) -> usize {
        5
    }
    fn autonomous(&self) -> bool {
        true
    }
    fn island(&self, _: f64, physical: &[f64], out: &mut [f64]) {
        self.island_calls.set(self.island_calls.get() + 1);
        let x: [f64; 5] = std::array::from_fn(|j| physical[j] / self.states[j]);
        let a = x[1] * x[1] + x[2] * x[2];
        let b = x[3] * x[3] + x[4] * x[4];
        let radial = self.gain * (self.p - a);
        let w = 1.2 + 0.3 * a + 0.4 * b;
        out.copy_from_slice(&[
            -(x[0] - 0.3 - 0.07 * self.p),
            radial * x[1] - w * x[2],
            w * x[1] + radial * x[2],
            0.7 * (x[1] * x[1] - x[2] * x[2] - x[3]),
            0.7 * (2.0 * x[1] * x[2] - x[4]),
        ]);
        for (value, unit) in out.iter_mut().zip(self.equations) {
            *value *= unit;
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let mut port = vec![C64::ZERO; 25];
        for j in 0..5 {
            port[j * 5 + j] = s.scale((1.0 + 0.15 * self.p) * self.equations[j] / self.states[j]);
        }
        port
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::exact_tangent)
    }
}
fn controls() -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: 64,
        maximum_unknowns: 128,
        maximum_payload_bytes: 16 << 20,
        maximum_work: 100_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_condition_estimate: 1e9,
        maximum_relative_residual: 1e-12,
        maximum_derivative_residual: 1e-10,
    }
}
fn reference(orbit: &HbOrbit) -> Vec<f64> {
    orbit.packed_seed(&HbAnchor::Autonomous {
        omega_guess: orbit.omega,
    })
}

#[test]
fn memory_frequency_mass_and_complete_waveform_match_independent_cycle_in_physical_units() {
    for (converted, phase) in [(false, 0.37), (true, -1.1)] {
        let problem = MemoryCycle::new(converted);
        let orbit = problem.orbit(phase);
        let mut prepared = PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &reference(&orbit),
            controls(),
            &problem.scales(),
            |point, out| problem.frequency_column(point, out),
            || false,
        )
        .unwrap();
        let calls = (problem.island_calls.get(), problem.tangent_calls.get());
        let material = prepared
            .solve(
                problem.parameter_unit,
                |point, out| problem.parameter_column(point, out),
                || false,
            )
            .unwrap();
        problem.verify(&orbit, &material, true);
        let frequency = prepared
            .solve_island_direction(
                problem.parameter_unit,
                |_, x, out| {
                    out.fill(0.0);
                    out[1] =
                        -x[2] / problem.states[2] * problem.equations[1] / problem.parameter_unit;
                    out[2] =
                        x[1] / problem.states[1] * problem.equations[2] / problem.parameter_unit;
                    true
                },
                || false,
            )
            .unwrap();
        problem.verify(&orbit, &frequency, false);
        assert_eq!(
            calls,
            (problem.island_calls.get(), problem.tangent_calls.get())
        );
        assert_eq!(prepared.dimension(), 5);
        assert_eq!(prepared.harmonics(), 3);
        assert!(prepared.condition_estimate().is_finite());
    }
}

#[test]
fn authored_phase_reference_and_full_unknown_budget_are_enforced_before_source_calls() {
    let problem = MemoryCycle::new(false);
    let orbit = problem.orbit(0.0);
    let mut budget = controls();
    budget.maximum_unknowns = 35;
    assert!(matches!(
        PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &reference(&orbit),
            budget,
            &problem.scales(),
            |point, out| problem.frequency_column(point, out),
            || false
        ),
        Err(ForcedSensitivityError::Budget { .. })
    ));
    assert!(matches!(
        PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &reference(&problem.orbit(0.2)),
            controls(),
            &problem.scales(),
            |point, out| problem.frequency_column(point, out),
            || false
        ),
        Err(ForcedSensitivityError::BaseUnresolved { .. })
    ));
    let mut zero = reference(&orbit);
    zero[5..35].fill(0.0);
    assert!(matches!(
        PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &zero,
            controls(),
            &problem.scales(),
            |point, out| problem.frequency_column(point, out),
            || false
        ),
        Err(ForcedSensitivityError::Orbit(OrbitError::TrivialCollapse))
    ));
    assert_eq!(problem.island_calls.get(), 0);
}

#[test]
fn complete_frequency_declaration_and_isolated_cycle_are_required() {
    let mut problem = MemoryCycle::new(false);
    let orbit = problem.orbit(0.4);
    assert!(matches!(
        PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &reference(&orbit),
            controls(),
            &problem.scales(),
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::FrequencyDirectionRefused { .. })
    ));
    problem.gain = 0.0;
    assert!(matches!(
        PreparedAutonomousSensitivity::new(
            &problem,
            &orbit,
            &reference(&orbit),
            controls(),
            &problem.scales(),
            |point, out| problem.frequency_column(point, out),
            || false
        ),
        Err(ForcedSensitivityError::Singular | ForcedSensitivityError::IllConditioned { .. })
    ));
}

#[test]
fn refused_unresolved_and_cancelled_directions_leave_original_factor_reusable() {
    let problem = MemoryCycle::new(false);
    let orbit = problem.orbit(-0.3);
    let prepared = PreparedAutonomousSensitivity::new(
        &problem,
        &orbit,
        &reference(&orbit),
        controls(),
        &problem.scales(),
        |point, out| problem.frequency_column(point, out),
        || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.solve(
            1.0,
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        prepared.solve(
            1.0,
            |point, out| {
                problem.parameter_column(point, out);
                if point.aft_samples == 64 {
                    out[0] += 0.1;
                }
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    assert!(matches!(
        prepared.solve(
            1.0,
            |point, out| problem.parameter_column(point, out),
            || true
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    let result = prepared
        .solve(
            1.0,
            |point, out| problem.parameter_column(point, out),
            || false,
        )
        .unwrap();
    problem.verify(&orbit, &result, true);
}
