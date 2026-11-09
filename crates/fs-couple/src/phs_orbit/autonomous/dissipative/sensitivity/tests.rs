//! G1/G3: independent loss-aware cycles, all physical parameter chains and work.
use super::*;
use crate::phs_orbit::autonomous::{
    DissipativeAutonomousHamiltonian, DissipativeAutonomousParameterDirection,
    DissipativeAutonomousPowerObjective as Objective, DissipativeAutonomousSensitivityError,
};
use fs_phs::{PhsStructureParameterDirection, harmonic::equilibrium::NonlinearDissipation};

// The existing original physical feedback and four stores are reused. Internal
// radial resistance and memory resistance are NEW, independently parameterized
// losses; the feedback does not cancel or adapt to either one.
#[derive(Clone, Copy)]
struct Loss {
    radial: f64,
    state: f64,
    memory: f64,
}
impl Loss {
    fn source() -> Self {
        Self {
            radial: 0.7,
            state: 0.18,
            memory: 0.12,
        }
    }
    fn zero() -> Self {
        Self {
            radial: 0.0,
            state: 0.0,
            memory: 0.0,
        }
    }
    fn shifted(self, p: Self, h: f64) -> Self {
        Self {
            radial: self.radial + h * p.radial,
            state: self.state + h * p.state,
            memory: self.memory + h * p.memory,
        }
    }
    fn flow(self, x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        let a = x[0] * x[0] + x[1] * x[1];
        let b = e[0] * e[0] + e[1] * e[1];
        let d = self.radial * (1.0 + self.state * a) * b;
        out.copy_from_slice(&[d * e[0], d * e[1], self.memory * e[2], self.memory * e[3]]);
        true
    }
    fn tangent(self, x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
        let a = x[0] * x[0] + x[1] * x[1];
        let da = 2.0 * (x[0] * dx[0] + x[1] * dx[1]);
        let b = e[0] * e[0] + e[1] * e[1];
        let db = 2.0 * (e[0] * de[0] + e[1] * de[1]);
        let d = self.radial * (1.0 + self.state * a) * b;
        let dd = self.radial * (self.state * da * b + (1.0 + self.state * a) * db);
        out.copy_from_slice(&[
            dd * e[0] + d * de[0],
            dd * e[1] + d * de[1],
            self.memory * de[2],
            self.memory * de[3],
        ]);
        true
    }
    fn parameter(self, p: Self, x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        let a = x[0] * x[0] + x[1] * x[1];
        let b = e[0] * e[0] + e[1] * e[1];
        let d = (p.radial * (1.0 + self.state * a) + self.radial * p.state * a) * b;
        out.copy_from_slice(&[d * e[0], d * e[1], p.memory * e[2], p.memory * e[3]]);
        true
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() < tolerance * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
fn zero_slope() -> Cycle {
    Cycle {
        j: 0.0,
        r: 0.0,
        g: 0.0,
        k: 0.0,
        eta: 0.0,
        mu: 0.0,
        h: 0.0,
        gain: 0.0,
        saturation: 0.0,
        shear: 0.0,
        beta: 0.0,
        alpha: 0.0,
    }
}

// Independent circular solution of the original four-state ODE. Radius is a
// positive quadratic root; memory's original pole determines its amplitude and
// phase. Memory back-reaction in the unchanged feedback determines frequency.
fn physical_orbit(c: &Cycle, d: Loss) -> HbOrbit {
    let linear = c.saturation * c.g.powi(4) + d.radial;
    let target = (c.g * c.g * c.gain - c.r) / (c.k * c.k);
    let a = 2.0 * target / (linear + (linear * linear + 4.0 * d.radial * d.state * target).sqrt());
    let gamma = (c.eta + d.memory) * c.mu;
    let k = c.k * c.j + c.k.powi(3) * c.g.powi(4) * c.shear * a;
    let b = c.k * c.g * c.g * c.beta * a * (c.h * c.alpha).powi(2);
    let (mut lo, mut hi) = (k, k + b / (gamma * gamma));
    for _ in 0..70 {
        let w = 0.5 * (lo + hi);
        if w - k - b / (gamma * gamma + w * w) > 0.0 {
            hi = w;
        } else {
            lo = w;
        }
    }
    let w = 0.5 * (lo + hi);
    let x = C64::new(0.43_f64.cos(), 0.43_f64.sin()).scale(0.5 * a.sqrt());
    let y = C64::new(0.0, -1.0) * x;
    let lag = C64::from_re(c.h * c.alpha) / C64::new(gamma, w);
    HbOrbit {
        omega: w,
        coeffs: [x, y, lag * x, lag * y]
            .into_iter()
            .map(|z| vec![C64::ZERO, z, C64::ZERO, C64::ZERO])
            .collect(),
        residual_trace: vec![],
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
fn aligned(base: &HbOrbit, mut candidate: HbOrbit) -> HbOrbit {
    let dot = base
        .coeffs
        .iter()
        .zip(&candidate.coeffs)
        .fold(C64::ZERO, |s, (a, b)| s + a[1].conj() * b[1]);
    let rotation = dot.conj().scale(1.0 / dot.abs());
    for row in &mut candidate.coeffs {
        row[1] = row[1] * rotation;
    }
    candidate
}
fn physical_power(c: &Cycle, d: Loss, orbit: &HbOrbit) -> [f64; 3] {
    let a = 4.0 * (orbit.coeffs[0][1].re.powi(2) + orbit.coeffs[0][1].im.powi(2));
    let b = 4.0 * (orbit.coeffs[2][1].re.powi(2) + orbit.coeffs[2][1].im.powi(2));
    let q = c.g * c.k;
    let supplied =
        (c.gain - c.saturation * q * q * a) * q * q * a + (c.eta + d.memory) * c.mu * c.mu * b;
    let resistance = c.r * c.k * c.k * a + c.eta * c.mu * c.mu * b;
    let nonlinear =
        d.radial * (1.0 + d.state * a) * c.k.powi(4) * a * a + d.memory * c.mu * c.mu * b;
    close(supplied, resistance + nonlinear, 2e-14);
    [supplied, resistance, nonlinear]
}

fn exercise(slope: Cycle, loss_slope: Loss, hold_loss: bool) {
    let c = Cycle::source();
    let d = Loss::source();
    let system = c.system();
    let orbit = physical_orbit(&c, d);
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| d.flow(x, e, out);
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        d.tangent(x, e, dx, de, out)
    };
    let jacobian = |x: &[f64], y: &[f64], out: &mut [f64]| c.jacobian(x, y, out);
    let source = DissipativeAutonomousHamiltonian::new(
        &system,
        &c,
        &jacobian,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 4,
        },
    )
    .unwrap();
    let mut prepared = source
        .prepare_sensitivity(&orbit, &reference(&orbit), controls(), &scales(), || false)
        .unwrap();
    let mut jp = vec![0.0; 16];
    jp[1] = -slope.j;
    jp[4] = slope.j;
    let rp = diagonal([slope.r, slope.r, slope.eta, slope.eta]);
    let gp = diagonal([slope.g, slope.g, slope.h, slope.h]);
    let storage = |x: &[f64], out: &mut [f64]| {
        for i in 0..4 {
            out[i] = x[i] * if i < 2 { slope.k } else { slope.mu };
        }
        true
    };
    let feedback =
        |x: &[f64], y: &[f64], yp: &[f64], out: &mut [f64]| c.parameter(&slope, x, y, yp, out);
    let loss_partial = |x: &[f64], e: &[f64], out: &mut [f64]| d.parameter(loss_slope, x, e, out);
    let direction = || DissipativeAutonomousParameterDirection {
        physical: PhsAutonomousParameterDirection {
            structure: PhsStructureParameterDirection {
                interconnection: Some(&jp),
                dissipation: Some(&rp),
                ports: Some(&gp),
            },
            storage: Some(&storage),
            feedback: &feedback,
        },
        dissipation: if hold_loss { None } else { Some(&loss_partial) },
    };
    let response = prepared
        .solve_parameter(0.2, direction(), || false)
        .unwrap();
    let step = 1e-5;
    let lower_c = c.shifted(&slope, -step);
    let upper_c = c.shifted(&slope, step);
    let lower_d = d.shifted(loss_slope, -step);
    let upper_d = d.shifted(loss_slope, step);
    let lower = aligned(&orbit, physical_orbit(&lower_c, lower_d));
    let upper = aligned(&orbit, physical_orbit(&upper_c, upper_d));
    let expected_frequency = (upper.omega - lower.omega) / (2.0 * step);
    assert!(
        expected_frequency.abs() > 1e-4,
        "loss/material frequency axis must be live"
    );
    close(
        response.angular_frequency_derivative,
        expected_frequency,
        2e-9,
    );
    let mut coefficient_norm = 0.0;
    for i in 0..4 {
        for k in 0..=3 {
            let expected = (upper.coeffs[i][k] - lower.coeffs[i][k]).scale(0.5 / step);
            close(response.coeffs[i][k].re, expected.re, 3e-9);
            close(response.coeffs[i][k].im, expected.im, 3e-9);
            coefficient_norm += expected.abs();
        }
    }
    assert!(coefficient_norm > 1e-3);
    let lower_power = physical_power(&lower_c, lower_d, &lower);
    let upper_power = physical_power(&upper_c, upper_d, &upper);
    let base_power = physical_power(&c, d, &orbit);
    // Build the original real coefficient derivative explicitly.
    let mut packed = vec![0.0; 28];
    for i in 0..4 {
        packed[i] = response.coeffs[i][0].re;
    }
    for k in 1..=3 {
        for i in 0..4 {
            packed[4 + 8 * (k - 1) + i] = response.coeffs[i][k].re;
            packed[8 + 8 * (k - 1) + i] = response.coeffs[i][k].im;
        }
    }
    for (index, mean, work) in [
        (
            0,
            Objective::MeanSuppliedPower,
            Objective::SuppliedEnergyPerCycle,
        ),
        (
            1,
            Objective::MeanResistancePower,
            Objective::ResistanceEnergyPerCycle,
        ),
        (
            2,
            Objective::MeanNonlinearPower,
            Objective::NonlinearEnergyPerCycle,
        ),
    ] {
        for (objective, per_cycle) in [(mean, false), (work, true)] {
            let observation = prepared
                .power_objective(objective, objective_controls(), || false)
                .unwrap();
            let period = if per_cycle {
                core::f64::consts::TAU / orbit.omega
            } else {
                1.0
            };
            close(observation.value, period * base_power[index], 3e-12);
            close(
                observation.angular_frequency_partial,
                if per_cycle {
                    -observation.value / orbit.omega
                } else {
                    0.0
                },
                3e-13,
            );
            let expected = if per_cycle {
                (core::f64::consts::TAU * upper_power[index] / upper.omega
                    - core::f64::consts::TAU * lower_power[index] / lower.omega)
                    / (2.0 * step)
            } else {
                (upper_power[index] - lower_power[index]) / (2.0 * step)
            };
            let mut dual = prepared
                .prepare_power_adjoint(objective, objective_controls(), || false)
                .unwrap();
            let result = dual.gradient_parameter(0.2, direction(), || false).unwrap();
            close(result.derivative, expected, 4e-9);
            assert!(result.gradient_sampling_defect < 1e-10);
            assert!(dual.prepared().dual_residual() < 1e-10);
            // Independent direct integrand partial at unchanged original state;
            // neighbor systems evaluate the ORIGINAL physical work, without
            // using the production parameter action or adjoint callback.
            let mut direct = 0.0;
            for j in 0..128 {
                let theta = core::f64::consts::TAU * j as f64 / 128.0;
                let x: Vec<_> = (0..4).map(|i| orbit.sample(i, theta)).collect();
                let value = |c: &Cycle, d: Loss| {
                    let e = [c.k * x[0], c.k * x[1], c.mu * x[2], c.mu * x[3]];
                    let y = [c.g * e[0], c.g * e[1], c.h * e[2], c.h * e[3]];
                    let a = x[0] * x[0] + x[1] * x[1];
                    let b = e[0] * e[0] + e[1] * e[1];
                    match index {
                        0 => {
                            (c.gain - c.saturation * (y[0] * y[0] + y[1] * y[1]))
                                * (y[0] * y[0] + y[1] * y[1])
                                + c.alpha * (x[0] * y[2] + x[1] * y[3])
                        }
                        1 => c.r * b + c.eta * (e[2] * e[2] + e[3] * e[3]),
                        _ => {
                            d.radial * (1.0 + d.state * a) * b * b
                                + d.memory * (e[2] * e[2] + e[3] * e[3])
                        }
                    }
                };
                direct += period * (value(&upper_c, upper_d) - value(&lower_c, lower_d))
                    / (2.0 * step * 128.0);
            }
            let forward = direct
                + observation
                    .gradient
                    .iter()
                    .zip(&packed)
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
                + observation.angular_frequency_partial * response.angular_frequency_derivative;
            close(result.derivative, forward, 4e-9);
        }
    }
}

#[test]
fn internal_loss_parameter_changes_live_frequency_waveform_and_independent_power_work() {
    exercise(
        zero_slope(),
        Loss {
            radial: 1.0,
            state: 0.0,
            memory: 0.0,
        },
        false,
    );
}
#[test]
fn simultaneous_storage_feedback_geometry_and_internal_loss_keep_all_chains_once() {
    exercise(
        Cycle::slope(),
        Loss {
            radial: 0.08,
            state: -0.025,
            memory: 0.04,
        },
        false,
    );
    exercise(Cycle::slope(), Loss::zero(), true);
}

#[test]
fn incomplete_loss_parameter_cancellation_and_resource_refusal_leave_factor_reusable() {
    let c = Cycle::source();
    let d = Loss::source();
    let system = c.system();
    let orbit = physical_orbit(&c, d);
    let flow = |x: &[f64], e: &[f64], out: &mut [f64]| d.flow(x, e, out);
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        d.tangent(x, e, dx, de, out)
    };
    let jacobian = |x: &[f64], y: &[f64], out: &mut [f64]| c.jacobian(x, y, out);
    let source = DissipativeAutonomousHamiltonian::new(
        &system,
        &c,
        &jacobian,
        NonlinearDissipation {
            flow: &flow,
            tangent: &tangent,
            max_tangent_actions: 4,
        },
    )
    .unwrap();
    let mut small = controls();
    small.maximum_payload_bytes = 8 * (6 * 4 + 4 * 4) - 1;
    assert!(matches!(
        source.prepare_sensitivity(&orbit, &reference(&orbit), small, &scales(), || false),
        Err(DissipativeAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Budget { .. }
        ))
    ));
    assert!(matches!(
        source.prepare_sensitivity(&orbit, &reference(&orbit), controls(), &scales(), || true),
        Err(DissipativeAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let mut prepared = source
        .prepare_sensitivity(&orbit, &reference(&orbit), controls(), &scales(), || false)
        .unwrap();
    let feedback = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let incomplete = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    let direction = || DissipativeAutonomousParameterDirection {
        physical: PhsAutonomousParameterDirection {
            structure: Default::default(),
            storage: None,
            feedback: &feedback,
        },
        dissipation: Some(&incomplete),
    };
    assert!(
        prepared
            .solve_parameter(1.0, direction(), || false)
            .is_err()
    );
    let mut dual = prepared
        .prepare_power_adjoint(
            Objective::NonlinearEnergyPerCycle,
            objective_controls(),
            || false,
        )
        .unwrap();
    assert!(dual.gradient_parameter(1.0, direction(), || false).is_err());
    let zero = || DissipativeAutonomousParameterDirection {
        physical: PhsAutonomousParameterDirection {
            structure: Default::default(),
            storage: None,
            feedback: &feedback,
        },
        dissipation: None,
    };
    assert!(matches!(
        dual.gradient_parameter(1.0, zero(), || true),
        Err(DissipativeAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    close(
        dual.gradient_parameter(1.0, zero(), || false)
            .unwrap()
            .derivative,
        0.0,
        1e-14,
    );
}
