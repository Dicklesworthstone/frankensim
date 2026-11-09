//! Native mixed port observations and physical cycle-mean work.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyParameterControls, HarmonicFrequencyParameterResponse, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 271,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-11,
    }
}

fn limits() -> HarmonicFrequencyParameterControls {
    HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}

fn port_limits() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}",
    );
}

fn close_real(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}",
    );
}

/// Closed rational quotient oracle; this never calls a response solver.
#[derive(Clone, Copy)]
struct Jet {
    v: C64,
    w: C64,
    ww: C64,
    t: C64,
    wt: C64,
}

impl Jet {
    const ZERO: Self = Self {
        v: C64::ZERO,
        w: C64::ZERO,
        ww: C64::ZERO,
        t: C64::ZERO,
        wt: C64::ZERO,
    };

    fn add(self, other: Self) -> Self {
        Self {
            v: self.v + other.v,
            w: self.w + other.w,
            ww: self.ww + other.ww,
            t: self.t + other.t,
            wt: self.wt + other.wt,
        }
    }

    fn mul(self, other: Self) -> Self {
        Self {
            v: self.v * other.v,
            w: self.w * other.v + self.v * other.w,
            ww: self.ww * other.v + (self.w * other.w).scale(2.0) + self.v * other.ww,
            t: self.t * other.v + self.v * other.t,
            wt: self.wt * other.v + self.w * other.t + self.t * other.w + self.v * other.wt,
        }
    }

    fn div(self, denominator: Self) -> Self {
        let d2 = denominator.v * denominator.v;
        let d3 = d2 * denominator.v;
        Self {
            v: self.v / denominator.v,
            w: self.w / denominator.v - self.v * denominator.w / d2,
            ww: self.ww / denominator.v
                - (self.w * denominator.w).scale(2.0) / d2
                - self.v * denominator.ww / d2
                + (self.v * denominator.w * denominator.w).scale(2.0) / d3,
            t: self.t / denominator.v - self.v * denominator.t / d2,
            wt: self.wt / denominator.v
                - (self.w * denominator.t + self.t * denominator.w + self.v * denominator.wt) / d2
                + (self.v * denominator.w * denominator.t).scale(2.0) / d3,
        }
    }
}

const M: f64 = 2.0;
const C: f64 = 0.3;
const K: f64 = 8.0;
const MT: f64 = 0.07;
const CT: f64 = 0.13;
const KT: f64 = -0.11;
const B: f64 = 1.2;
const BT: f64 = 0.18;
const KB: f64 = 2.1;
const TAU: f64 = 0.8;
const TAUT: f64 = 0.23;

fn source(w: f64) -> Jet {
    let (f0, f1, f2) = (
        C64::new(1.7, -0.6),
        C64::new(-0.2, 0.08),
        C64::new(0.03, -0.04),
    );
    let (g0, g1) = (C64::new(0.1, 0.7), C64::new(-0.06, 0.02));
    Jet {
        v: f0 + f1.scale(w) + f2.scale(0.5 * w * w),
        w: f1 + f2.scale(w),
        ww: f2,
        t: g0 + g1.scale(w),
        wt: g1,
    }
}

fn load(w: f64) -> Jet {
    let s = C64::new(1.0, -w * TAU);
    let s2 = s * s;
    let s3 = s2 * s;
    Jet {
        v: C64::from_re(KB * TAU) / s,
        w: C64::new(0.0, KB * TAU * TAU) / s2,
        ww: C64::from_re(-2.0 * KB * TAU * TAU * TAU) / s3,
        t: C64::from_re(KB * TAUT) / s2,
        wt: C64::new(0.0, 2.0 * KB * TAU * TAUT) / s3,
    }
}

fn expected(w: f64) -> (Jet, Jet, Jet, Jet) {
    let mechanical = Jet {
        v: C64::new(K - w * w * M, -w * C),
        w: C64::new(-2.0 * w * M, -C),
        ww: C64::from_re(-2.0 * M),
        t: C64::new(KT - w * w * MT, -w * CT),
        wt: C64::new(-2.0 * w * MT, -CT),
    };
    let s = Jet {
        v: C64::new(1.0, -w * TAU),
        w: C64::new(0.0, -TAU),
        t: C64::new(0.0, -w * TAUT),
        wt: C64::new(0.0, -TAUT),
        ..Jet::ZERO
    };
    let minus_i_w = Jet {
        v: C64::new(0.0, -w),
        w: C64::new(0.0, -1.0),
        ..Jet::ZERO
    };
    let load_numerator = Jet {
        v: C64::from_re(B * B * KB * TAU),
        t: C64::from_re(KB * (2.0 * B * BT * TAU + B * B * TAUT)),
        ..Jet::ZERO
    };
    // The original Maxwell memory equation after clearing 1-i*w*tau:
    // q = F*(1-i*w*tau) / ((K-w^2*M-i*w*C)*(1-i*w*tau)-i*w*B^2*Kb*tau).
    let q = source(w)
        .mul(s)
        .div(mechanical.mul(s).add(minus_i_w.mul(load_numerator)));
    let v = minus_i_w.mul(q);
    let attachment = Jet {
        v: C64::from_re(B),
        t: C64::from_re(BT),
        ..Jet::ZERO
    };
    let flow = attachment.mul(v);
    let effort = load(w).mul(flow);
    (q, v, flow, effort)
}

fn solve_fixture(w: f64, units: f64, cx: &Cx<'_>) -> HarmonicFrequencyParameterResponse {
    let f = source(w);
    let z = load(w);
    let inverse_square_units = 1.0 / (units * units);
    let owner = PreparedHarmonicResponse::try_with_port_impedance(
        &[M],
        &[C],
        &[K],
        1,
        &[B * units],
        1,
        &[z.v.scale(inverse_square_units)],
        w,
        controls(),
        port_limits(),
        cx,
    )
    .unwrap();
    owner
        .solve_with_frequency_parameter_and_port_derivative(
            &[f.v],
            QuadraticParameterDirection {
                mass: &[MT],
                damping: &[CT],
                stiffness: &[KT],
            },
            Some(&[f.w]),
            Some(&[f.ww]),
            Some(&[f.t]),
            Some(&[f.wt]),
            &[z.w.scale(inverse_square_units)],
            &[z.ww.scale(inverse_square_units)],
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    force_map: &[BT * units],
                    impedance: &[z.t.scale(inverse_square_units)],
                },
                impedance_angular_parameter_derivative: &[z.wt.scale(inverse_square_units)],
            },
            limits(),
            cx,
        )
        .unwrap()
}

fn real_inner(a: C64, b: C64) -> f64 {
    a.re * b.re + a.im * b.im
}

fn expected_power(w: f64, velocity: Jet, flow: Jet) -> (f64, f64) {
    let viscous = C * (real_inner(velocity.wt, velocity.v) + real_inner(velocity.w, velocity.t))
        + CT * real_inner(velocity.w, velocity.v);
    // Differentiate the passive branch's real constitutive work 0.5*Re(Z)*|u|^2,
    // independently of the implementation's four complex effort/flow products.
    let d = 1.0 + w * w * TAU * TAU;
    let resistance = KB * TAU / d;
    let resistance_w = -2.0 * KB * TAU * TAU * TAU * w / (d * d);
    let resistance_t = KB * TAUT * (1.0 - w * w * TAU * TAU) / (d * d);
    let resistance_wt = -2.0 * KB * w * TAU * TAU * TAUT * (3.0 - w * w * TAU * TAU) / (d * d * d);
    let port = 0.5 * resistance_wt * flow.v.norm_sq()
        + resistance_w * real_inner(flow.v, flow.t)
        + resistance_t * real_inner(flow.v, flow.w)
        + resistance * (real_inner(flow.wt, flow.v) + real_inner(flow.w, flow.t));
    (viscous, port)
}

fn check_fixture(result: &HarmonicFrequencyParameterResponse, w: f64, units: f64) {
    let (q, v, flow, effort) = expected(w);
    let first = &result.frequency_response.frequency_response;
    close(first.response.displacement[0], q.v);
    close(first.displacement_derivative[0], q.w);
    close(
        result.frequency_response.displacement_second_derivative[0],
        q.ww,
    );
    close(result.parameter_response.displacement_derivative[0], q.t);
    close(result.displacement_mixed_derivative[0], q.wt);
    close(result.velocity_mixed_derivative[0], v.wt);

    let original_port = first.response.port_load.as_ref().unwrap();
    let first_port = first.port_load.as_ref().unwrap();
    let curvature_port = result.frequency_response.port_load.as_ref().unwrap();
    let parameter_port = result.parameter_response.port_load.as_ref().unwrap();
    let port = result.port_load.as_ref().unwrap();
    for (actual, wanted) in [
        (original_port.flow[0], flow.v),
        (first_port.flow_derivative[0], flow.w),
        (curvature_port.flow_second_derivative[0], flow.ww),
        (parameter_port.flow_derivative[0], flow.t),
        (port.flow_mixed_derivative[0], flow.wt),
    ] {
        close(actual.scale(1.0 / units), wanted);
    }
    for (actual, wanted) in [
        (original_port.effort[0], effort.v),
        (first_port.effort_derivative[0], effort.w),
        (curvature_port.effort_second_derivative[0], effort.ww),
        (parameter_port.effort_derivative[0], effort.t),
        (port.effort_mixed_derivative[0], effort.wt),
    ] {
        close(actual.scale(units), wanted);
    }
    let (viscous, port_work) = expected_power(w, v, flow);
    close_real(result.average_viscous_power_mixed_derivative, viscous);
    close_real(port.average_absorbed_power_mixed_derivative, port_work);
    close_real(
        result.average_input_power_mixed_derivative,
        viscous + port_work,
    );
    let envelope = viscous.abs().max(port_work.abs()).max(1.0);
    assert!(result.power_balance_residual_mixed_derivative.abs() < 1e-11 * envelope);
    assert!(result.relative_mixed_derivative_residual < 2e-12);
}

#[test]
fn changing_geometry_and_maxwell_memory_preserve_complete_port_jets_and_mixed_work() {
    with_cx(|cx| {
        for w in [0.0, 0.2, 1.5, 1.99, 2.7, 8.0] {
            check_fixture(&solve_fixture(w, 1.0, cx), w, 1.0);
        }
    });
}

#[test]
fn original_port_units_preserve_mixed_motion_effort_flow_and_signed_work() {
    with_cx(|cx| {
        for exponent in [-400, 0, 400] {
            let units = f64::from_bits(((1023 + exponent) as u64) << 52);
            for w in [0.0, 1.99, 2.7] {
                check_fixture(&solve_fixture(w, units, cx), w, units);
            }
        }
    });
}

#[test]
fn reactive_mixed_source_work_cancels_before_overflowing_products_are_formed() {
    with_cx(|cx| {
        let amplitude = f64::from_bits(((1023 + 600) as u64) << 52);
        let owner =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.0], &[2.0], 1, 1.0, controls(), cx)
                .unwrap();
        let zero = [0.0];
        // A pure source-phase parameter has no physical loss. At w=1,
        // q=F, q_w=2F, q_theta=iF, q_wtheta=2iF and v_wtheta=3F.
        // The two nonzero source-work contractions are exactly +/-1.5*F^2;
        // neither is representable on its own, while their sum is zero.
        let result = owner
            .solve_with_frequency_parameter_derivative(
                &[C64::from_re(amplitude)],
                QuadraticParameterDirection {
                    mass: &zero,
                    damping: &zero,
                    stiffness: &zero,
                },
                None,
                None,
                Some(&[C64::new(0.0, amplitude)]),
                None,
                limits(),
                cx,
            )
            .unwrap();
        assert_eq!(
            result.velocity_mixed_derivative,
            vec![C64::from_re(3.0 * amplitude)]
        );
        assert_eq!(result.average_input_power_mixed_derivative, 0.0);
        assert_eq!(result.average_viscous_power_mixed_derivative, 0.0);
        assert_eq!(result.power_balance_residual_mixed_derivative, 0.0);
        assert!(result.port_load.is_none());
    });
}

#[test]
fn changing_gyroscope_changes_motion_without_creating_direct_viscous_work() {
    with_cx(|cx| {
        let mass = [1.0, 0.0, 0.0, 2.0];
        let damping = [0.0, 0.7, -0.7, 0.0];
        let stiffness = [5.0, -0.6, -0.6, 9.0];
        let zero = [0.0; 4];
        let damping_direction = [0.0, 0.11, -0.11, 0.0];
        for w in [0.0, 0.5, 1.4, 2.6] {
            let owner = PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &stiffness,
                2,
                w,
                controls(),
                cx,
            )
            .unwrap();
            let f = source(w);
            let result = owner
                .solve_with_frequency_parameter_derivative(
                    &[f.v, C64::new(-0.3, 0.7)],
                    QuadraticParameterDirection {
                        mass: &zero,
                        damping: &damping_direction,
                        stiffness: &zero,
                    },
                    Some(&[f.w, C64::new(0.2, -0.1)]),
                    Some(&[f.ww, C64::new(-0.04, 0.03)]),
                    Some(&[f.t, C64::new(0.1, 0.4)]),
                    Some(&[f.wt, C64::new(-0.07, 0.05)]),
                    limits(),
                    cx,
                )
                .unwrap();
            assert_eq!(result.average_viscous_power_mixed_derivative, 0.0);
            assert!(result.average_input_power_mixed_derivative.abs() < 1e-10);
            assert!(result.power_balance_residual_mixed_derivative.abs() < 1e-10);
            assert!(result.port_load.is_none());
        }
    });
}
