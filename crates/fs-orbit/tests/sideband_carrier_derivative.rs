//! Exact original-port carrier derivatives, checked without frequency differencing.

use fs_math::c64::C64;
use fs_orbit::sideband::{
    PreparedPeriodicLinearization, SidebandControls, SidebandError, SidebandScales,
    sensitivity::SidebandCarrierDerivativeControls,
};
use fs_orbit::{HbOrbit, IslandStateJacobian, OrbitProblem};
use std::cell::Cell;

fn controls(harmonics: usize) -> SidebandControls {
    SidebandControls {
        harmonics,
        samples_per_period: 128,
        maximum_unknowns: 128,
        maximum_payload_bytes: 32 << 20,
        maximum_work: 1_000_000_000,
        maximum_base_residual: 1e-11,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_relative_residual: 1e-12,
        maximum_omitted_residual: 1e-9,
    }
}
fn derivative_controls() -> SidebandCarrierDerivativeControls {
    SidebandCarrierDerivativeControls {
        maximum_payload_bytes: 8 << 20,
        maximum_work: 1_000_000_000,
        maximum_relative_residual: 1e-12,
        maximum_omitted_residual: 1e-9,
    }
}
fn base(omega: f64) -> HbOrbit {
    HbOrbit {
        omega,
        coeffs: vec![vec![C64::ZERO; 2]; 2],
        residual_trace: vec![],
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:?} expected={expected:?} error={:?}",
        (actual - expected).abs()
    );
}
fn inverse_pair(a: &[C64], ap: &[C64], f: [C64; 2], fp: [C64; 2]) -> ([C64; 2], [C64; 2]) {
    let det = a[0] * a[3] - a[1] * a[2];
    let detp = ap[0] * a[3] + a[0] * ap[3] - ap[1] * a[2] - a[1] * ap[2];
    let numerators = [a[3] * f[0] - a[1] * f[1], a[0] * f[1] - a[2] * f[0]];
    let slopes = [
        ap[3] * f[0] + a[3] * fp[0] - ap[1] * f[1] - a[1] * fp[1],
        ap[0] * f[1] + a[0] * fp[1] - ap[2] * f[0] - a[2] * fp[0],
    ];
    let x = numerators.map(|z| z / det);
    (
        x,
        [
            (slopes[0] - x[0] * detp) / det,
            (slopes[1] - x[1] * detp) / det,
        ],
    )
}

struct DispersivePort;
impl DispersivePort {
    fn jacobian(&self, _: f64, _: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        true
    }
    fn derivative(&self, w: f64) -> Vec<C64> {
        let i = C64::new(0.0, 1.0);
        let relaxation = C64::ONE + i.scale(0.4 * w);
        let memory_derivative = i.scale(1.2) / (relaxation * relaxation);
        let mass = [2.0, 0.3, 0.3, 1.5];
        let loss_and_gyro = [0.3, 0.47, -0.33, 0.5];
        let memory = [1.0, -0.5, -0.5, 0.25];
        (0..4)
            .map(|j| {
                C64::from_re(-2.0 * w * mass[j])
                    + i.scale(loss_and_gyro[j])
                    + memory_derivative.scale(memory[j])
            })
            .collect()
    }
}
impl OrbitProblem for DispersivePort {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, _: f64, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let relaxation = s.scale(1.2) / (C64::ONE + s.scale(0.4));
        let mass = [2.0, 0.3, 0.3, 1.5];
        let loss_and_gyro = [0.3, 0.47, -0.33, 0.5];
        let stiffness = [4.0, 0.2, 0.2, 2.0];
        let memory = [1.0, -0.5, -0.5, 0.25];
        (0..4)
            .map(|j| {
                (s * s).scale(mass[j])
                    + s.scale(loss_and_gyro[j])
                    + C64::from_re(stiffness[j])
                    + relaxation.scale(memory[j])
            })
            .collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}

fn source_at(w: f64, rotation: C64) -> ([C64; 2], [C64; 2]) {
    let i = C64::new(0.0, 1.0);
    let a = rotation * C64::new(0.2, -0.1);
    let b = rotation * C64::new(-0.04, 0.15);
    let pole = C64::ONE + i.scale(0.2 * w);
    (
        [a * (C64::ONE + i.scale(0.3 * w)), b / pole],
        [a * i.scale(0.3), -b * i.scale(0.2) / (pole * pole)],
    )
}

#[test]
fn full_dispersive_gyroscopic_derivative_and_varying_excitation_match_closed_inverse() {
    let problem = DispersivePort;
    let orbit = base(1.1);
    let scales = SidebandScales {
        states: vec![1e-3, 1e2],
        equations: vec![7.0, 1e-2],
    };
    let prepared =
        PreparedPeriodicLinearization::new(&problem, &orbit, controls(1), &scales, || false)
            .unwrap();
    for carrier in [-1.4, 0.73] {
        let factor = prepared.factor(carrier, || false).unwrap();
        let calls = Cell::new(0);
        let derivative = factor
            .prepare_carrier_derivative(
                derivative_controls(),
                |w, out| {
                    calls.set(calls.get() + 1);
                    out.copy_from_slice(&problem.derivative(w));
                    true
                },
                || false,
            )
            .unwrap();
        assert_eq!(calls.get(), 3);
        assert_eq!(derivative.dimension(), 2);
        assert_eq!(derivative.harmonics(), 1);
        assert_eq!(derivative.carrier_angular_frequency(), carrier);
        for rotation in [C64::ONE, C64::new(-0.3, 0.7)] {
            let mut input = Vec::new();
            let mut input_slope = Vec::new();
            for k in -1..=1 {
                let (f, fp) = source_at(carrier + k as f64 * orbit.omega, rotation);
                input.extend(f);
                input_slope.extend(fp);
            }
            let pair = derivative.solve(&input, &input_slope, || false).unwrap();
            for band in 0..3 {
                let w = carrier + (band as f64 - 1.0) * orbit.omega;
                let (f, fp) = source_at(w, rotation);
                let (expected, slope) = inverse_pair(
                    &problem.port(C64::new(0.0, w)),
                    &problem.derivative(w),
                    f,
                    fp,
                );
                for row in 0..2 {
                    near(
                        pair.response.coefficients[band * 2 + row],
                        expected[row],
                        2e-11,
                    );
                    near(
                        pair.angular_derivative.coefficients[band * 2 + row],
                        slope[row],
                        2e-10,
                    );
                }
            }
            assert_eq!(
                calls.get(),
                3,
                "derivative blocks must be cached across excitations"
            );
            assert_eq!(pair.response.omitted_residual, 0.0);
            assert_eq!(pair.angular_derivative.omitted_residual, 0.0);
            assert!(pair.angular_derivative.relative_residual < 1e-12);
        }
    }
}

struct Conversion {
    omega: f64,
    phase: f64,
    coupling: f64,
}
impl Conversion {
    fn jacobian(&self, t: f64, _: &[f64], out: &mut [f64]) -> bool {
        out.copy_from_slice(&[
            -0.2,
            0.0,
            self.coupling * (self.omega * t + self.phase).cos(),
            -0.7,
        ]);
        true
    }
    fn port_slope(_: f64, out: &mut [C64]) -> bool {
        out.fill(C64::ZERO);
        out[0] = C64::new(0.0, 1.0);
        out[3] = out[0];
        true
    }
}
impl OrbitProblem for Conversion {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = -0.2 * x[0];
        out[1] = self.coupling * (self.omega * t + self.phase).cos() * x[0] - 0.7 * x[1];
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![s, C64::ZERO, C64::ZERO, s]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}

#[test]
fn converted_sideband_slopes_keep_full_coupling_phase_and_signed_frequencies() {
    let source = Conversion {
        omega: 1.3,
        phase: 0.37,
        coupling: 0.8,
    };
    let orbit = base(source.omega);
    let units = SidebandScales {
        states: vec![2.0, 0.03],
        equations: vec![0.4, 7.0],
    };
    let prepared =
        PreparedPeriodicLinearization::new(&source, &orbit, controls(1), &units, || false).unwrap();
    let force = C64::new(0.3, -0.2);
    let input = [C64::ZERO, C64::ZERO, force, C64::ZERO, C64::ZERO, C64::ZERO];
    let i = C64::new(0.0, 1.0);
    for carrier in [-0.43, 0.29] {
        let factor = prepared.factor(carrier, || false).unwrap();
        let derivative = factor
            .prepare_carrier_derivative(derivative_controls(), Conversion::port_slope, || false)
            .unwrap();
        let pair = derivative.solve(&input, &[C64::ZERO; 6], || false).unwrap();
        let rate = C64::new(0.2, carrier);
        let x = force / rate;
        let xp = -i * x / rate;
        for k in -1_isize..=1 {
            let row = (k + 1) as usize * 2;
            near(
                pair.response.coefficients[row],
                if k == 0 { x } else { C64::ZERO },
                2e-12,
            );
            near(
                pair.angular_derivative.coefficients[row],
                if k == 0 { xp } else { C64::ZERO },
                2e-11,
            );
            let (y, yp) = if k == 0 {
                (C64::ZERO, C64::ZERO)
            } else {
                let phase = k as f64 * source.phase;
                let coupling = C64::new(phase.cos(), phase.sin()).scale(0.5 * source.coupling);
                let rate = C64::new(0.7, carrier + k as f64 * source.omega);
                (
                    coupling * x / rate,
                    coupling * (xp / rate - i * x / (rate * rate)),
                )
            };
            near(pair.response.coefficients[row + 1], y, 2e-12);
            near(pair.angular_derivative.coefficients[row + 1], yp, 2e-11);
        }
    }
}

#[test]
fn derivative_has_its_own_omitted_equation_allowance() {
    let source = Conversion {
        omega: 1.3,
        phase: 0.37,
        coupling: 0.1,
    };
    let orbit = base(source.omega);
    let units = SidebandScales {
        states: vec![1.0; 2],
        equations: vec![1.0; 2],
    };
    let mut c = controls(0);
    c.maximum_omitted_residual = 1.0;
    let prepared =
        PreparedPeriodicLinearization::new(&source, &orbit, c, &units, || false).unwrap();
    let factor = prepared.factor(0.3, || false).unwrap();
    let input = [C64::ONE, C64::ZERO];
    factor.solve(&input, || false).unwrap();
    let mut dc = derivative_controls();
    dc.maximum_omitted_residual = 1e-10;
    let derivative = factor
        .prepare_carrier_derivative(dc, Conversion::port_slope, || false)
        .unwrap();
    assert!(matches!(derivative.solve(&input,&[C64::ZERO;2],||false),
        Err(SidebandError::OmittedSidebands{measured,limit}) if measured>0.1 && limit==1e-10));
}

#[test]
fn derivative_admission_refusal_and_cancelled_reuse_do_not_publish_partial_pairs() {
    let source = DispersivePort;
    let orbit = base(1.1);
    let units = SidebandScales {
        states: vec![1.0; 2],
        equations: vec![1.0; 2],
    };
    let prepared =
        PreparedPeriodicLinearization::new(&source, &orbit, controls(1), &units, || false).unwrap();
    let factor = prepared.factor(0.7, || false).unwrap();
    let calls = Cell::new(0);
    let mut tiny = derivative_controls();
    tiny.maximum_payload_bytes = 0;
    assert!(matches!(
        factor.prepare_carrier_derivative(
            tiny,
            |_, _| {
                calls.set(calls.get() + 1);
                true
            },
            || false
        ),
        Err(SidebandError::Budget { .. })
    ));
    tiny = derivative_controls();
    tiny.maximum_work = 0;
    assert!(matches!(
        factor.prepare_carrier_derivative(
            tiny,
            |_, _| {
                calls.set(calls.get() + 1);
                true
            },
            || false
        ),
        Err(SidebandError::Budget { .. })
    ));
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        factor.prepare_carrier_derivative(derivative_controls(), |_, _| false, || false),
        Err(SidebandError::CarrierDerivativeRefused { sideband: -1, .. })
    ));
    assert!(matches!(
        factor.prepare_carrier_derivative(
            derivative_controls(),
            |_, out| {
                out[0] = C64::ZERO;
                true
            },
            || false
        ),
        Err(SidebandError::CarrierDerivativeRefused { .. })
    ));
    let derivative = factor
        .prepare_carrier_derivative(
            derivative_controls(),
            |w, out| {
                out.copy_from_slice(&source.derivative(w));
                true
            },
            || false,
        )
        .unwrap();
    let forcing = [C64::new(0.1, 0.2); 6];
    assert!(matches!(
        derivative.solve(&forcing, &[C64::ZERO; 5], || false),
        Err(SidebandError::Invalid(_))
    ));
    let mut bad = [C64::ZERO; 6];
    bad[4] = C64::new(f64::NAN, 0.0);
    assert!(matches!(
        derivative.solve(&forcing, &bad, || false),
        Err(SidebandError::Invalid(_))
    ));
    let checkpoints = Cell::new(0);
    assert!(matches!(
        derivative.solve(&forcing, &[C64::ZERO; 6], || {
            checkpoints.set(checkpoints.get() + 1);
            checkpoints.get() > 20
        }),
        Err(SidebandError::Cancelled)
    ));
    let first = derivative
        .solve(&forcing, &[C64::ZERO; 6], || false)
        .unwrap();
    let second = derivative
        .solve(&forcing, &[C64::ZERO; 6], || false)
        .unwrap();
    assert_eq!(first.response.coefficients, second.response.coefficients);
    assert_eq!(
        first.angular_derivative.coefficients,
        second.angular_derivative.coefficients
    );
}
