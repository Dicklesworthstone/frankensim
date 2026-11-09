use super::*;
use fs_math::{c64::C64, det};
use fs_orbit::HbAnchor;
use fs_phs::QuadraticStorage;

#[path = "power/automatic_tests.rs"]
mod automatic_power;
#[path = "power/tests.rs"]
mod power;
// Reuse this independent physical cycle fixture for the complete original D law.
#[path = "../dissipative/sensitivity/tests.rs"]
mod nonlinear_loss;

// An active saturating physical-port source feeds a passive rotating pair and
// two passive relaxation stores. Neither material resistance nor stored energy
// is replaced by a fitted amplitude/decay law. The circular solution is an
// independent manufactured oracle, not a production model reduction.
struct Cycle {
    j: f64,
    r: f64,
    g: f64,
    k: f64,
    eta: f64,
    mu: f64,
    h: f64,
    gain: f64,
    saturation: f64,
    shear: f64,
    beta: f64,
    alpha: f64,
}
impl Cycle {
    fn source() -> Self {
        Self {
            j: 1.1,
            r: 0.2,
            g: 0.9,
            k: 1.3,
            eta: 0.6,
            mu: 0.9,
            h: 0.8,
            gain: 1.2,
            saturation: 0.8,
            shear: 0.15,
            beta: 0.25,
            alpha: 0.7,
        }
    }
    fn slope() -> Self {
        Self {
            j: 0.11,
            r: 0.07,
            g: 0.05,
            k: 0.09,
            eta: 0.02,
            mu: -0.03,
            h: 0.04,
            gain: 0.13,
            saturation: 0.03,
            shear: -0.02,
            beta: 0.06,
            alpha: 0.08,
        }
    }
    fn system(&self) -> PortHamiltonian {
        let mut j = vec![0.0; 16];
        j[1] = -self.j;
        j[4] = self.j;
        PortHamiltonian::new(
            4,
            4,
            j,
            diagonal([self.r, self.r, self.eta, self.eta]),
            diagonal([self.g, self.g, self.h, self.h]),
            Box::new(
                QuadraticStorage::new(diagonal([self.k, self.k, self.mu, self.mu]), 4).unwrap(),
            ),
        )
        .unwrap()
    }
    fn radius_squared(&self) -> f64 {
        (self.g * self.g * self.gain - self.r)
            / (self.saturation * self.g.powi(4) * self.k * self.k)
    }
    fn frequency(&self) -> f64 {
        let a = self.radius_squared();
        let gamma = self.eta * self.mu;
        let k = self.k * self.j + self.k.powi(3) * self.g.powi(4) * self.shear * a;
        let b = self.k * self.g * self.g * self.beta * a * (self.h * self.alpha).powi(2);
        let (mut lo, mut hi) = (k, k + b / (gamma * gamma));
        for _ in 0..70 {
            let w = 0.5 * (lo + hi);
            if w - k - b / (gamma * gamma + w * w) > 0.0 {
                hi = w;
            } else {
                lo = w;
            }
        }
        0.5 * (lo + hi)
    }
    fn orbit(&self) -> HbOrbit {
        let w = self.frequency();
        let x =
            C64::new(det::cos(0.43), det::sin(0.43)).scale(0.5 * det::sqrt(self.radius_squared()));
        let y = C64::new(0.0, -1.0) * x;
        let lag = C64::from_re(self.h * self.alpha) / C64::new(self.eta * self.mu, w);
        let coeffs = [x, y, lag * x, lag * y]
            .into_iter()
            .map(|z| vec![C64::ZERO, z, C64::ZERO, C64::ZERO])
            .collect();
        HbOrbit {
            omega: w,
            coeffs,
            residual_trace: vec![],
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }
    fn jacobian(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        let a = flows[0] * flows[0] + flows[1] * flows[1];
        let b = x[2] * x[2] + x[3] * x[3];
        let radial = self.gain - self.saturation * a;
        let angular = self.shear * a + self.beta * b;
        let q = self.g * self.k;
        out.fill(0.0);
        for column in 0..4 {
            let da = if column < 2 {
                2.0 * flows[column] * q
            } else {
                0.0
            };
            let dr = -self.saturation * da;
            let dw = self.shear * da
                + if column >= 2 {
                    2.0 * self.beta * x[column]
                } else {
                    0.0
                };
            let dy0 = if column == 0 { q } else { 0.0 };
            let dy1 = if column == 1 { q } else { 0.0 };
            out[column] = dr * flows[0] + radial * dy0 - dw * flows[1] - angular * dy1;
            out[4 + column] = dw * flows[0] + angular * dy0 + dr * flows[1] + radial * dy1;
        }
        out[8] = self.alpha;
        out[13] = self.alpha;
        true
    }
    fn parameter(
        &self,
        slope: &Self,
        x: &[f64],
        flows: &[f64],
        dy: &[f64],
        out: &mut [f64],
    ) -> bool {
        for j in 0..4 {
            let expected = if j < 2 {
                (slope.g * self.k + self.g * slope.k) * x[j]
            } else {
                (slope.h * self.mu + self.h * slope.mu) * x[j]
            };
            assert!((dy[j] - expected).abs() < 1e-13);
        }
        let a = flows[0] * flows[0] + flows[1] * flows[1];
        let da = 2.0 * (flows[0] * dy[0] + flows[1] * dy[1]);
        let b = x[2] * x[2] + x[3] * x[3];
        let radial = self.gain - self.saturation * a;
        let angular = self.shear * a + self.beta * b;
        let dr = slope.gain - slope.saturation * a - self.saturation * da;
        let dw = slope.shear * a + self.shear * da + slope.beta * b;
        out[0] = dr * flows[0] + radial * dy[0] - dw * flows[1] - angular * dy[1];
        out[1] = dw * flows[0] + angular * dy[0] + dr * flows[1] + radial * dy[1];
        out[2] = slope.alpha * x[0];
        out[3] = slope.alpha * x[1];
        true
    }
    fn verify(&self, slope: &Self, orbit: &HbOrbit, result: &AutonomousOrbitSensitivity) {
        let a = self.radius_squared();
        let n = self.g * self.g * self.gain - self.r;
        let dn = 2.0 * self.g * slope.g * self.gain + self.g * self.g * slope.gain - slope.r;
        let da = a
            * (dn / n
                - slope.saturation / self.saturation
                - 4.0 * slope.g / self.g
                - 2.0 * slope.k / self.k);
        let gamma = self.eta * self.mu;
        let dg = slope.eta * self.mu + self.eta * slope.mu;
        let shear = self.k.powi(3) * self.g.powi(4) * self.shear * a;
        let dk = slope.k * self.j
            + self.k * slope.j
            + shear
                * (3.0 * slope.k / self.k
                    + 4.0 * slope.g / self.g
                    + slope.shear / self.shear
                    + da / a);
        let b = self.k * self.g * self.g * self.beta * a * (self.h * self.alpha).powi(2);
        let db = b
            * (slope.k / self.k
                + 2.0 * slope.g / self.g
                + slope.beta / self.beta
                + da / a
                + 2.0 * slope.h / self.h
                + 2.0 * slope.alpha / self.alpha);
        let denominator = gamma * gamma + orbit.omega * orbit.omega;
        let dw = (dk + db / denominator - 2.0 * b * gamma * dg / (denominator * denominator))
            / (1.0 + 2.0 * b * orbit.omega / (denominator * denominator));
        assert!((result.angular_frequency_derivative - dw).abs() < 2e-11);
        let pole = C64::new(gamma, orbit.omega);
        let lag = C64::from_re(self.h * self.alpha) / pole;
        let dlag = C64::from_re(slope.h * self.alpha + self.h * slope.alpha) / pole
            - C64::new(dg, dw).scale(self.h * self.alpha) / (pole * pole);
        let base: Vec<_> = orbit.coeffs.iter().map(|row| row[1]).collect();
        let radial = 0.5 * da / a;
        let mut derivative = [
            base[0].scale(radial),
            base[1].scale(radial),
            dlag * base[0] + (lag * base[0]).scale(radial),
            dlag * base[1] + (lag * base[1]).scale(radial),
        ];
        let dot: f64 = base
            .iter()
            .zip(&derivative)
            .map(|(z, dz)| -z.im * dz.re + z.re * dz.im)
            .sum();
        let norm: f64 = base.iter().map(|z| z.re * z.re + z.im * z.im).sum();
        for j in 0..4 {
            derivative[j] = derivative[j] + C64::new(0.0, -dot / norm) * base[j];
            assert!(
                (result.coeffs[j][1] - derivative[j]).abs() < 3e-11,
                "state {j}"
            );
            assert!(result.coeffs[j][0].abs() < 1e-11);
            assert!(result.coeffs[j][2..].iter().all(|z| z.abs() < 1e-11));
        }
        assert!(result.derivative_residual < 1e-10);
    }
}
impl AutonomousPortFeedback for Cycle {
    fn port_count(&self) -> usize {
        4
    }
    fn input(&self, x: &[f64], flows: &[f64], out: &mut [f64]) -> bool {
        let a = flows[0] * flows[0] + flows[1] * flows[1];
        let b = x[2] * x[2] + x[3] * x[3];
        let radial = self.gain - self.saturation * a;
        let angular = self.shear * a + self.beta * b;
        out[0] = radial * flows[0] - angular * flows[1];
        out[1] = angular * flows[0] + radial * flows[1];
        out[2] = self.alpha * x[0];
        out[3] = self.alpha * x[1];
        true
    }
}
fn diagonal(values: [f64; 4]) -> Vec<f64> {
    let mut out = vec![0.0; 16];
    for i in 0..4 {
        out[5 * i] = values[i];
    }
    out
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
fn scales() -> HbScales {
    HbScales {
        states: vec![1.0; 4],
        equations: vec![1.0; 4],
        angular_frequency: 2.0,
    }
}
fn reference(orbit: &HbOrbit) -> Vec<f64> {
    orbit.packed_seed(&HbAnchor::Autonomous {
        omega_guess: orbit.omega,
    })
}

#[test]
fn physical_structure_storage_feedback_and_memory_slopes_match_independent_cycle() {
    let cycle = Cycle::source();
    let slope = Cycle::slope();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let feedback = |x: &[f64], flows: &[f64], dy: &[f64], out: &mut [f64]| {
        cycle.parameter(&slope, x, flows, dy, out)
    };
    let storage = |x: &[f64], out: &mut [f64]| {
        for j in 0..4 {
            out[j] = x[j] * if j < 2 { slope.k } else { slope.mu };
        }
        true
    };
    let mut dj = vec![0.0; 16];
    dj[1] = -slope.j;
    dj[4] = slope.j;
    let dr = diagonal([slope.r, slope.r, slope.eta, slope.eta]);
    let dg = diagonal([slope.g, slope.g, slope.h, slope.h]);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let result = prepared
        .solve_parameter(
            1.0,
            PhsAutonomousParameterDirection {
                structure: PhsStructureParameterDirection {
                    interconnection: Some(&dj),
                    dissipation: Some(&dr),
                    ports: Some(&dg),
                },
                storage: Some(&storage),
                feedback: &feedback,
            },
            || false,
        )
        .unwrap();
    cycle.verify(&slope, &orbit, &result);
    assert!(prepared.prepared().base_residual() < 1e-11);
}

#[test]
fn declared_storage_hessian_reuses_the_same_autonomous_owner() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let hessian = |_: &[f64], direction: &[f64], out: &mut [f64]| {
        for j in 0..4 {
            out[j] = direction[j] * if j < 2 { cycle.k } else { cycle.mu };
        }
        true
    };
    let mut native = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let mut analytic = PreparedPhsAutonomousSensitivity::new_analytic(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        &hessian,
        || false,
    )
    .unwrap();
    // A direct drive-strength change with fixed material matrices/storage.
    let direction = |_: &[f64], flows: &[f64], dy: &[f64], out: &mut [f64]| {
        assert!(dy.iter().all(|v| *v == 0.0));
        out.fill(0.0);
        out[0] = flows[0];
        out[1] = flows[1];
        true
    };
    let make = || PhsAutonomousParameterDirection {
        structure: Default::default(),
        storage: None,
        feedback: &direction,
    };
    let a = native.solve_parameter(0.2, make(), || false).unwrap();
    let b = analytic.solve_parameter(0.2, make(), || false).unwrap();
    assert!((a.angular_frequency_derivative - b.angular_frequency_derivative).abs() < 1e-12);
    for (x, y) in a.coeffs.iter().flatten().zip(b.coeffs.iter().flatten()) {
        assert!((*x - *y).abs() < 1e-12);
    }
}

#[test]
fn refused_feedback_direction_and_cancellation_leave_the_factor_reusable() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let incomplete = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    let make = || PhsAutonomousParameterDirection {
        structure: Default::default(),
        storage: None,
        feedback: &incomplete,
    };
    assert!(matches!(
        prepared.solve_parameter(1.0, make(), || false),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::IslandDirectionRefused { .. }
        ))
    ));
    assert!(matches!(
        prepared.solve_parameter(1.0, make(), || true),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let zero = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let result = prepared
        .solve_parameter(
            1.0,
            PhsAutonomousParameterDirection {
                structure: Default::default(),
                storage: None,
                feedback: &zero,
            },
            || false,
        )
        .unwrap();
    assert_eq!(result.angular_frequency_derivative, 0.0);
    assert!(result.coeffs.iter().flatten().all(|z| *z == C64::ZERO));
}

fn objective_controls() -> ForcedMeanObjectiveControls {
    ForcedMeanObjectiveControls {
        objective_scale: 1.0,
        maximum_value_sampling_defect: 1e-10,
        maximum_gradient_sampling_defect: 1e-10,
    }
}

impl Cycle {
    fn shifted(&self, slope: &Self, h: f64) -> Self {
        Self {
            j: self.j + h * slope.j,
            r: self.r + h * slope.r,
            g: self.g + h * slope.g,
            k: self.k + h * slope.k,
            eta: self.eta + h * slope.eta,
            mu: self.mu + h * slope.mu,
            h: self.h + h * slope.h,
            gain: self.gain + h * slope.gain,
            saturation: self.saturation + h * slope.saturation,
            shear: self.shear + h * slope.shear,
            beta: self.beta + h * slope.beta,
            alpha: self.alpha + h * slope.alpha,
        }
    }

    // Independent closed-form physical observations, including both relaxing
    // stores and their frequency-dependent phase/amplitude.
    fn energy_and_cycle_loss(&self) -> (f64, f64) {
        let a = self.radius_squared();
        let w = self.frequency();
        let b = a * (self.h * self.alpha).powi(2) / ((self.eta * self.mu).powi(2) + w * w);
        (
            0.5 * (self.k * a + self.mu * b),
            core::f64::consts::TAU / w
                * (self.r * self.k * self.k * a + self.eta * self.mu * self.mu * b),
        )
    }
}

#[test]
fn physical_objectives_include_all_material_feedback_and_changing_period_chains() {
    let cycle = Cycle::source();
    let slope = Cycle::slope();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let feedback = |x: &[f64], flows: &[f64], dy: &[f64], out: &mut [f64]| {
        cycle.parameter(&slope, x, flows, dy, out)
    };
    let storage = |x: &[f64], out: &mut [f64]| {
        for j in 0..4 {
            out[j] = x[j] * if j < 2 { slope.k } else { slope.mu };
        }
        true
    };
    let mut dj = vec![0.0; 16];
    dj[1] = -slope.j;
    dj[4] = slope.j;
    let dr = diagonal([slope.r, slope.r, slope.eta, slope.eta]);
    let dg = diagonal([slope.g, slope.g, slope.h, slope.h]);
    let direction = || PhsAutonomousParameterDirection {
        structure: PhsStructureParameterDirection {
            interconnection: Some(&dj),
            dissipation: Some(&dr),
            ports: Some(&dg),
        },
        storage: Some(&storage),
        feedback: &feedback,
    };
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let w = orbit.omega;
    let period = core::f64::consts::TAU / w;
    let a = cycle.radius_squared();
    let b = a * (cycle.h * cycle.alpha).powi(2) / ((cycle.eta * cycle.mu).powi(2) + w * w);
    let h = 2e-5;
    let lo = cycle.shifted(&slope, -h).energy_and_cycle_loss();
    let hi = cycle.shifted(&slope, h).energy_and_cycle_loss();
    let base = cycle.energy_and_cycle_loss();
    for per_cycle in [false, true] {
        let objective = if per_cycle {
            prepared
                .mean_objective(
                    objective_controls(),
                    |_, x, value, gradient, frequency| {
                        let effort = [
                            cycle.k * x[0],
                            cycle.k * x[1],
                            cycle.mu * x[2],
                            cycle.mu * x[3],
                        ];
                        *value = period * system.resistance_power_from_effort(&effort).unwrap();
                        for j in 0..4 {
                            let coefficient = if j < 2 {
                                cycle.r * cycle.k * cycle.k
                            } else {
                                cycle.eta * cycle.mu * cycle.mu
                            };
                            gradient[j] = 2.0 * period * coefficient * x[j];
                        }
                        *frequency = -*value / w;
                        true
                    },
                    || false,
                )
                .unwrap()
        } else {
            prepared
                .mean_stored_energy(objective_controls(), || false)
                .unwrap()
        };
        let explicit = if per_cycle {
            period
                * ((slope.r * cycle.k * cycle.k + 2.0 * cycle.r * cycle.k * slope.k) * a
                    + (slope.eta * cycle.mu * cycle.mu + 2.0 * cycle.eta * cycle.mu * slope.mu) * b)
        } else {
            0.5 * (slope.k * a + slope.mu * b)
        };
        let mut adjoint = prepared
            .prepare_adjoint(
                &objective.gradient,
                objective.angular_frequency_partial,
                1.0,
                || false,
            )
            .unwrap();
        let derivative = adjoint
            .gradient_parameter(0.2, explicit, direction(), || false)
            .unwrap();
        let (expected, reference) = if per_cycle {
            (base.1, (hi.1 - lo.1) / (2.0 * h))
        } else {
            (base.0, (hi.0 - lo.0) / (2.0 * h))
        };
        assert!((objective.value - expected).abs() < 3e-12);
        assert!(
            (derivative.derivative - reference).abs() < 3e-9,
            "per_cycle={per_cycle}: {} versus {reference}",
            derivative.derivative
        );
        assert!(adjoint.prepared().dual_residual() < 1e-10);
    }
}

#[test]
fn physical_frequency_objective_uses_the_original_bordered_cycle() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let feedback = |_: &[f64], flows: &[f64], dy: &[f64], out: &mut [f64]| {
        assert!(dy.iter().all(|v| *v == 0.0));
        out.copy_from_slice(&[flows[0], flows[1], 0.0, 0.0]);
        true
    };
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    let objective = prepared
        .mean_objective(
            objective_controls(),
            |_, _, value, gradient, frequency| {
                *value = orbit.omega;
                gradient.fill(0.0);
                *frequency = 1.0;
                true
            },
            || false,
        )
        .unwrap();
    let mut adjoint = prepared
        .prepare_adjoint(
            &objective.gradient,
            objective.angular_frequency_partial,
            1.0,
            || false,
        )
        .unwrap();
    let result = adjoint
        .gradient_parameter(
            0.2,
            0.0,
            PhsAutonomousParameterDirection {
                structure: Default::default(),
                storage: None,
                feedback: &feedback,
            },
            || false,
        )
        .unwrap();
    let h = 1e-5;
    let mut lo = Cycle::source();
    lo.gain -= h;
    let mut hi = Cycle::source();
    hi.gain += h;
    assert!((result.derivative - (hi.frequency() - lo.frequency()) / (2.0 * h)).abs() < 1e-9);
    assert!(result.derivative.abs() > 1e-2);
}

#[test]
fn physical_autonomous_objective_refuses_missing_clock_or_feedback_and_recovers() {
    let cycle = Cycle::source();
    let system = cycle.system();
    let orbit = cycle.orbit();
    let jacobian = |x: &[f64], flows: &[f64], out: &mut [f64]| cycle.jacobian(x, flows, out);
    let mut prepared = PreparedPhsAutonomousSensitivity::new(
        &system,
        &cycle,
        &jacobian,
        &orbit,
        &reference(&orbit),
        controls(),
        &scales(),
        || false,
    )
    .unwrap();
    assert!(
        prepared
            .mean_objective(
                objective_controls(),
                |_, _, value, gradient, _| {
                    *value = 0.0;
                    gradient.fill(0.0);
                    true
                },
                || false
            )
            .is_err()
    );
    let objective = prepared
        .mean_stored_energy(objective_controls(), || false)
        .unwrap();
    let mut adjoint = prepared
        .prepare_adjoint(
            &objective.gradient,
            objective.angular_frequency_partial,
            1.0,
            || false,
        )
        .unwrap();
    let incomplete = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 0.0;
        true
    };
    assert!(
        adjoint
            .gradient_parameter(
                1.0,
                0.0,
                PhsAutonomousParameterDirection {
                    structure: Default::default(),
                    storage: None,
                    feedback: &incomplete,
                },
                || false
            )
            .is_err()
    );
    let zero = |_: &[f64], _: &[f64], _: &[f64], out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let direction = || PhsAutonomousParameterDirection {
        structure: Default::default(),
        storage: None,
        feedback: &zero,
    };
    let mut polls = 0;
    assert!(matches!(
        adjoint.gradient_parameter(1.0, 0.0, direction(), || {
            polls += 1;
            polls > 8
        }),
        Err(PhsAutonomousSensitivityError::Sensitivity(
            ForcedSensitivityError::Cancelled
        ))
    ));
    let result = adjoint
        .gradient_parameter(1.0, 0.7, direction(), || false)
        .unwrap();
    assert_eq!(result.derivative, 0.7);
}
