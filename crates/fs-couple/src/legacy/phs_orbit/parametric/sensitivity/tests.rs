//! Nonlinear original pumped physics, with independent motion and power oracles.

use super::*;
use crate::phs_orbit::PhsForcedParameterDirection;
use fs_orbit::sensitivity::{
    ForcedMeanObjectiveControls, ForcedSensitivityControls, ForcedSensitivityDirection,
    ForcedSensitivityError,
};
use fs_phs::PhsStructureParameterDirection;

#[path = "loss_tests.rs"]
mod loss;

#[derive(Clone, Copy, Default)]
struct Family {
    j: f64,
    loss: f64,
    g: f64,
    k: f64,
    beta: f64,
    mu: f64,
    force: f64,
    drive_phase: f64,
    pump: f64,
    pump_phase: f64,
    cross: f64,
    cross_phase: f64,
    omega: f64,
}
impl Family {
    fn source() -> Self {
        Self {
            j: 1.0,
            loss: 1.0,
            g: 1.0,
            k: 1.8,
            beta: 0.6,
            mu: 0.8,
            force: 0.2,
            drive_phase: -0.2,
            pump: 0.25,
            pump_phase: 0.4,
            cross: 0.12,
            cross_phase: 0.3,
            omega: 0.8,
        }
    }
    fn slope() -> Self {
        Self {
            j: 0.11,
            loss: 0.13,
            g: 0.07,
            k: 0.09,
            beta: -0.04,
            mu: 0.05,
            force: 0.08,
            drive_phase: 0.17,
            pump: 0.06,
            pump_phase: -0.13,
            cross: 0.03,
            cross_phase: 0.08,
            omega: 0.12,
        }
    }
    fn shifted(self, p: Self, h: f64) -> Self {
        Self {
            j: self.j + h * p.j,
            loss: self.loss + h * p.loss,
            g: self.g + h * p.g,
            k: self.k + h * p.k,
            beta: self.beta + h * p.beta,
            mu: self.mu + h * p.mu,
            force: self.force + h * p.force,
            drive_phase: self.drive_phase + h * p.drive_phase,
            pump: self.pump + h * p.pump,
            pump_phase: self.pump_phase + h * p.pump_phase,
            cross: self.cross + h * p.cross,
            cross_phase: self.cross_phase + h * p.cross_phase,
            omega: self.omega + h * p.omega,
        }
    }
    fn matrices(self) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            vec![0.0, self.j, 0.0, -self.j, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![
                0.03 * self.loss,
                0.0,
                0.0,
                0.0,
                0.32 * self.loss,
                0.0,
                0.0,
                0.0,
                0.7 * self.loss,
            ],
            vec![0.1 * self.g, self.g, 0.2 * self.g],
        )
    }
    fn system(self) -> PortHamiltonian {
        let (j, r, g) = self.matrices();
        PortHamiltonian::new(3, 1, j, r, g, Box::new(Energy(self))).unwrap()
    }
    fn input(self) -> PeriodicPortInput {
        PeriodicPortInput::new(
            self.omega,
            vec![vec![
                C64::from_re(0.2 * self.force),
                complex(0.5 * self.force, self.drive_phase),
                C64::ZERO,
                complex(0.1 * self.force, 2.0 * self.drive_phase + 0.7),
            ]],
        )
        .unwrap()
    }
    fn pump(self) -> PeriodicQuadraticModulation {
        let mut q = vec![vec![C64::ZERO; 9]; 3];
        q[0][0] = C64::from_re(0.28 * self.pump);
        q[1][2] = complex(0.5 * self.cross, self.cross_phase);
        q[1][6] = q[1][2];
        q[2][0] = complex(0.5 * self.pump, self.pump_phase);
        PeriodicQuadraticModulation::new(self.omega, 3, q).unwrap()
    }
    fn input_partial(self, p: Self) -> PeriodicPortInput {
        PeriodicPortInput::new(
            self.omega,
            vec![vec![
                C64::from_re(0.2 * p.force),
                coefficient_partial(
                    0.5 * self.force,
                    self.drive_phase,
                    0.5 * p.force,
                    p.drive_phase,
                ),
                C64::ZERO,
                coefficient_partial(
                    0.1 * self.force,
                    2.0 * self.drive_phase + 0.7,
                    0.1 * p.force,
                    2.0 * p.drive_phase,
                ),
            ]],
        )
        .unwrap()
    }
    fn pump_partial(self, p: Self) -> PeriodicQuadraticModulation {
        let mut q = vec![vec![C64::ZERO; 9]; 3];
        q[0][0] = C64::from_re(0.28 * p.pump);
        q[1][2] = coefficient_partial(
            0.5 * self.cross,
            self.cross_phase,
            0.5 * p.cross,
            p.cross_phase,
        );
        q[1][6] = q[1][2];
        q[2][0] = coefficient_partial(0.5 * self.pump, self.pump_phase, 0.5 * p.pump, p.pump_phase);
        PeriodicQuadraticModulation::new(self.omega, 3, q).unwrap()
    }
    // Explicit original physical equations and three independent accounts.
    // This oracle calls no production field, power, modulation or tangent owner.
    fn physical(self, time: f64, x: &[f64; 3]) -> ([f64; 3], [f64; 3]) {
        let theta = self.omega * time;
        let q = self.pump * (0.28 + (2.0 * theta + self.pump_phase).cos());
        let cross = self.cross * (theta + self.cross_phase).cos();
        let e = [
            self.k * x[0]
                + self.beta * x[0].powi(3)
                + self.mu * (x[0] - x[2])
                + q * x[0]
                + cross * x[2],
            x[1] / 1.2,
            self.mu * (x[2] - x[0]) + cross * x[0],
        ];
        let u = self.force
            * (0.2
                + (theta + self.drive_phase).cos()
                + 0.2 * (3.0 * theta + 2.0 * self.drive_phase + 0.7).cos());
        let pump = -self.omega * self.pump * (2.0 * theta + self.pump_phase).sin() * x[0] * x[0]
            - self.omega * self.cross * (theta + self.cross_phase).sin() * x[0] * x[2];
        (
            [
                self.j * e[1] - 0.03 * self.loss * e[0] + 0.1 * self.g * u,
                -self.j * e[0] - 0.32 * self.loss * e[1] + self.g * u,
                -0.7 * self.loss * e[2] + 0.2 * self.g * u,
            ],
            [
                pump,
                self.g * u * (0.1 * e[0] + e[1] + 0.2 * e[2]),
                self.loss * (0.03 * e[0] * e[0] + 0.32 * e[1] * e[1] + 0.7 * e[2] * e[2]),
            ],
        )
    }
    fn solve(self) -> HbOrbit {
        let (system, input, pump) = (self.system(), self.input(), self.pump());
        hb(
            &ParametricHamiltonian::new(&system, &input, &pump).unwrap(),
            3,
            9,
        )
    }
    fn means(self, orbit: &HbOrbit) -> [f64; 3] {
        let mut means = [0.0; 3];
        for sample in 0..1024 {
            let theta = TAU * sample as f64 / 1024.0;
            let x = core::array::from_fn(|state| {
                orbit.coeffs[state][0].re
                    + orbit.coeffs[state]
                        .iter()
                        .enumerate()
                        .skip(1)
                        .map(|(n, z)| {
                            2.0 * (z.re * (n as f64 * theta).cos()
                                - z.im * (n as f64 * theta).sin())
                        })
                        .sum::<f64>()
            });
            let powers = self.physical(theta / self.omega, &x).1;
            for account in 0..3 {
                means[account] += powers[account] / 1024.0;
            }
        }
        means
    }
}

struct Energy(Family);
impl Storage for Energy {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let p = self.0;
        0.5 * (p.k * x[0] * x[0] + x[1] * x[1] / 1.2 + p.mu * (x[2] - x[0]).powi(2))
            + 0.25 * p.beta * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let p = self.0;
        out[0] = p.k * x[0] + p.beta * x[0].powi(3) + p.mu * (x[0] - x[2]);
        out[1] = x[1] / 1.2;
        out[2] = p.mu * (x[2] - x[0]);
    }
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }
}
impl ContinuousStorageTangent for Energy {
    fn apply(&self, x: &[f64], v: &[f64], out: &mut [f64]) -> bool {
        let p = self.0;
        out[0] = (p.k + p.mu + 3.0 * p.beta * x[0] * x[0]) * v[0] - p.mu * v[2];
        out[1] = v[1] / 1.2;
        out[2] = p.mu * (v[2] - v[0]);
        true
    }
}
fn coefficient_partial(a: f64, phase: f64, ap: f64, phasep: f64) -> C64 {
    complex(ap, phase) + C64::new(0.0, phasep) * complex(a, phase)
}
fn controls() -> ForcedSensitivityControls {
    ForcedSensitivityControls {
        samples_per_period: 128,
        maximum_unknowns: 256,
        maximum_payload_bytes: 16 << 20,
        maximum_work: 500_000_000,
        maximum_base_residual: 1e-9,
        maximum_tangent_sampling_defect: 1e-9,
        maximum_condition_estimate: 1e10,
        maximum_relative_residual: 1e-11,
        maximum_derivative_residual: 1e-9,
    }
}
fn objectives() -> ForcedMeanObjectiveControls {
    ForcedMeanObjectiveControls {
        objective_scale: 1.0,
        maximum_value_sampling_defect: 1e-10,
        maximum_gradient_sampling_defect: 1e-9,
    }
}
fn chart(omega_partial: f64) -> ForcedSensitivityDirection {
    ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 0.4,
        angular_frequency_derivative: omega_partial,
    }
}
fn empty_direction<'a>() -> ParametricParameterDirection<'a> {
    ParametricParameterDirection {
        physical: PhsForcedParameterDirection {
            structure: Default::default(),
            storage: None,
            input: None,
        },
        modulation: None,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn complete_pumped_waveform_and_accounts_follow_original_neighboring_systems() {
    let f = Family::source();
    let p = Family::slope();
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let source = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let orbit = hb(&source, 3, 9);
    let mut prepared = source
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega), || false)
        .unwrap();
    let (jp, rp, gp) = p.matrices();
    let (up, qp) = (f.input_partial(p), f.pump_partial(p));
    let hp = |x: &[f64], out: &mut [f64]| {
        out[0] = p.k * x[0] + p.beta * x[0].powi(3) + p.mu * (x[0] - x[2]);
        out[1] = 0.0;
        out[2] = p.mu * (x[2] - x[0]);
        true
    };
    let direction = || ParametricParameterDirection {
        physical: PhsForcedParameterDirection {
            structure: PhsStructureParameterDirection {
                interconnection: Some(&jp),
                dissipation: Some(&rp),
                ports: Some(&gp),
            },
            storage: Some(&hp),
            input: Some(&up),
        },
        modulation: Some(&qp),
    };
    let response = prepared
        .solve_parameter(chart(p.omega), direction(), || false)
        .unwrap();
    let epsilon = 2e-5;
    let (lower, upper) = (f.shifted(p, -epsilon), f.shifted(p, epsilon));
    let (lo, hi) = (lower.solve(), upper.solve());
    for ((actual, a), b) in response
        .coeffs
        .iter()
        .flatten()
        .zip(lo.coeffs.iter().flatten())
        .zip(hi.coeffs.iter().flatten())
    {
        assert!((*actual - (*b - *a).scale(0.5 / epsilon)).abs() < 3e-7);
    }
    // Independent original-ODE RK4 motion has no HB/parameter implementation
    // in common. The nonlinear relaxing state and converted harmonics remain.
    let timed_lo = time_coefficients(lower.omega, 1024, 9, |time, x| lower.physical(time, x).0);
    let timed_hi = time_coefficients(upper.omega, 1024, 9, |time, x| upper.physical(time, x).0);
    for ((actual, a), b) in response
        .coeffs
        .iter()
        .flatten()
        .zip(timed_lo.iter().flatten())
        .zip(timed_hi.iter().flatten())
    {
        assert!((*actual - (*b - *a).scale(0.5 / epsilon)).abs() < 3e-6);
    }
    assert!(orbit.coeffs[2][1].abs() > 1e-3);
    assert!(orbit.coeffs[0][3].abs() > 1e-4);
    let (base, low, high) = (f.means(&orbit), lower.means(&lo), upper.means(&hi));
    assert!(base[0].abs() > 1e-7);
    assert!((base[1] - base[2]).abs() > 1e-7);
    for (index, objective) in [
        ParametricPowerObjective::MeanPumpPower,
        ParametricPowerObjective::MeanPortPower,
        ParametricPowerObjective::MeanResistancePower,
        ParametricPowerObjective::PumpEnergyPerCycle,
        ParametricPowerObjective::PortEnergyPerCycle,
        ParametricPowerObjective::ResistanceEnergyPerCycle,
    ]
    .into_iter()
    .enumerate()
    {
        let cycle = index >= 3;
        let mut adjoint = prepared
            .prepare_power_adjoint(objective, objectives(), || false)
            .unwrap();
        close(
            adjoint.value(),
            base[index % 3] * if cycle { TAU / f.omega } else { 1.0 },
            2e-12,
        );
        let derivative = adjoint
            .gradient_parameter(chart(p.omega), direction(), || false)
            .unwrap();
        let expected = if cycle {
            TAU * (high[index % 3] / upper.omega - low[index % 3] / lower.omega) / (2.0 * epsilon)
        } else {
            (high[index % 3] - low[index % 3]) / (2.0 * epsilon)
        };
        close(derivative.derivative, expected, 3e-7);
    }
}

#[test]
fn pump_clock_partial_and_cycle_work_have_the_correct_distinct_chains() {
    let f = Family::source();
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        Energy(f).apply(x, v, out)
    };
    let source = ParametricHamiltonian::new_analytic(&system, &input, &pump, &hessian).unwrap();
    let orbit = hb(&source, 3, 9);
    let mut prepared = source
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega), || false)
        .unwrap();
    let clock = 0.12;
    let response = prepared
        .solve_parameter(chart(clock), empty_direction(), || false)
        .unwrap();
    let mut packed = vec![0.0; 57];
    for state in 0..3 {
        packed[state] = response.coeffs[state][0].re;
    }
    for n in 1..=9 {
        for state in 0..3 {
            let re = 3 + (n - 1) * 6 + state;
            packed[re] = response.coeffs[state][n].re;
            packed[re + 3] = response.coeffs[state][n].im;
        }
    }
    let p = Family {
        omega: clock,
        ..Default::default()
    };
    let epsilon = 2e-5;
    let (lo, hi) = (f.shifted(p, -epsilon), f.shifted(p, epsilon));
    let (low, high) = (lo.means(&lo.solve()), hi.means(&hi.solve()));
    for cycle in [false, true] {
        let objective = if cycle {
            ParametricPowerObjective::PumpEnergyPerCycle
        } else {
            ParametricPowerObjective::MeanPumpPower
        };
        let observation = prepared
            .power_objective(objective, objectives(), || false)
            .unwrap();
        let chain: f64 = observation
            .gradient
            .iter()
            .zip(&packed)
            .map(|(g, x)| g * x)
            .sum();
        let explicit = if cycle {
            0.0
        } else {
            clock / f.omega * observation.value
        };
        let mut adjoint = prepared
            .prepare_power_adjoint(objective, objectives(), || false)
            .unwrap();
        let before = calls.get();
        let result = adjoint
            .gradient_parameter(chart(clock), empty_direction(), || false)
            .unwrap();
        assert_eq!(calls.get(), before);
        close(result.derivative - chain, explicit, 3e-12);
        let expected = if cycle {
            TAU * (high[0] / hi.omega - low[0] / lo.omega) / (2.0 * epsilon)
        } else {
            (high[0] - low[0]) / (2.0 * epsilon)
        };
        close(result.derivative, expected, 2e-7);
    }
    assert!(calls.get() > 0);
}

#[test]
fn pumped_direction_admission_budget_refusal_and_cancellation_leave_source_reusable() {
    let f = Family::source();
    let (system, input, pump) = (f.system(), f.input(), f.pump());
    let source = ParametricHamiltonian::new(&system, &input, &pump).unwrap();
    let orbit = hb(&source, 3, 9);
    let mut tiny = controls();
    tiny.maximum_payload_bytes = 1;
    assert!(matches!(
        source.prepare_sensitivity(&orbit, tiny, &scales(3, f.omega), || false),
        Err(ParametricError::Sensitivity(
            ForcedSensitivityError::Budget {
                quantity: "parametric physical direction payload",
                ..
            }
        ))
    ));
    let mut prepared = source
        .prepare_sensitivity(&orbit, controls(), &scales(3, f.omega), || false)
        .unwrap();
    let wrong =
        PeriodicQuadraticModulation::new(f.omega * 1.1, 3, vec![vec![C64::ZERO; 9]]).unwrap();
    assert!(matches!(
        prepared.solve_parameter(
            chart(0.0),
            ParametricParameterDirection {
                modulation: Some(&wrong),
                ..empty_direction()
            },
            || false
        ),
        Err(ParametricError::Invalid(_))
    ));
    let mut too_wide = vec![vec![C64::ZERO; 9]; 11];
    too_wide[10][0] = C64::from_re(0.1);
    let too_wide = PeriodicQuadraticModulation::new(f.omega, 3, too_wide).unwrap();
    assert!(matches!(
        prepared.solve_parameter(
            chart(0.0),
            ParametricParameterDirection {
                modulation: Some(&too_wide),
                ..empty_direction()
            },
            || false
        ),
        Err(ParametricError::Invalid(_))
    ));
    // A new derivative harmonic above the BASE pump bandwidth is valid when
    // retained by the original orbit; its sampler has its own admitted scratch.
    let mut new_band = vec![vec![C64::ZERO; 9]; 4];
    new_band[3][0] = complex(0.03, 0.2);
    let new_band = PeriodicQuadraticModulation::new(f.omega, 3, new_band).unwrap();
    let direction = || ParametricParameterDirection {
        modulation: Some(&new_band),
        ..empty_direction()
    };
    assert!(
        prepared
            .solve_parameter(chart(0.0), direction(), || false)
            .unwrap()
            .coeffs
            .iter()
            .flatten()
            .any(|z| z.abs() > 1e-7)
    );
    let mut adjoint = prepared
        .prepare_power_adjoint(
            ParametricPowerObjective::PumpEnergyPerCycle,
            objectives(),
            || false,
        )
        .unwrap();
    let baseline = adjoint
        .gradient_parameter(chart(0.0), direction(), || false)
        .unwrap();
    let incomplete = |_: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    let mut bad = direction();
    bad.physical.storage = Some(&incomplete);
    assert!(matches!(
        adjoint.gradient_parameter(chart(0.0), bad, || false),
        Err(ParametricError::SourceRefused)
    ));
    let mut polls = 0;
    assert!(matches!(
        adjoint.gradient_parameter(chart(0.0), direction(), || {
            polls += 1;
            polls > 100
        }),
        Err(ParametricError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let retry = adjoint
        .gradient_parameter(chart(0.0), direction(), || false)
        .unwrap();
    assert_eq!(baseline, retry);
}
