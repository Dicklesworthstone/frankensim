use super::*;
use fs_math::{c64::C64, det};
use fs_orbit::HbAnchor;
use fs_phs::QuadraticStorage;

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
