//! Authored channel-current source: cutoff, saturation and asymmetric reversal.

use super::*;
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelPoint, GateChannelResponse,
};

// This bounded manufactured source has a square-law knee, exact cutoff and
// a reverse scale of 0.4. It owns channel current only; body diodes are separate.
#[derive(Clone, Debug)]
struct SaturatingChannel;
impl GateChannelResponse for SaturatingChannel {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-1.0, 4.0]
    }
    fn evaluate(&self, vd: f64, vg: f64) -> Result<GateChannelPoint, PhsError> {
        if !vd.is_finite() || !vg.is_finite() || vd.abs() > 20.0 || !(-1.0..=4.0).contains(&vg) {
            return Err(PhsError::Dimension {
                what: "authored current source domain",
            });
        }
        let u = (vg - 1.0).max(0.0);
        let a = vd.abs();
        let scale = if vd < 0.0 { 0.4 } else { 1.0 };
        let (magnitude, drain, gate) = if a < u {
            (2.0 * (2.0 * u * a - a * a), 4.0 * (u - a), 4.0 * a)
        } else {
            (2.0 * u * u, 0.0, 4.0 * u)
        };
        Ok(GateChannelPoint {
            current_a: vd.signum() * scale * magnitude,
            drain_derivative_s: scale * drain,
            gate_derivative_s: vd.signum() * scale * gate,
        })
    }
}

#[test]
fn original_current_law_retains_cutoff_saturation_reverse_asymmetry_and_both_partials() {
    let law = GateChannelLaw::new(SaturatingChannel).unwrap();
    for (vd, vg, expected) in [
        (0.2, 1.5, [0.32, 1.2, 0.8]),
        (1.0, 1.5, [0.5, 0.0, 2.0]),
        (-1.0, 1.5, [-0.2, 0.0, -0.8]),
        (-0.2, 1.5, [-0.128, 0.48, -0.32]),
    ] {
        let point = law.evaluate(vd, vg).unwrap();
        for (actual, expected) in [
            point.current_a,
            point.drain_derivative_s,
            point.gate_derivative_s,
        ]
        .into_iter()
        .zip(expected)
        {
            close(actual, expected, 5e-15);
        }
        let h = 1e-6;
        close(
            (law.evaluate(vd + h, vg).unwrap().current_a
                - law.evaluate(vd - h, vg).unwrap().current_a)
                / (2.0 * h),
            point.drain_derivative_s,
            3e-10,
        );
        close(
            (law.evaluate(vd, vg + h).unwrap().current_a
                - law.evaluate(vd, vg - h).unwrap().current_a)
                / (2.0 * h),
            point.gate_derivative_s,
            3e-10,
        );
        assert!(point.current_a * vd > 0.0);
    }
    for vd in [-20.0, -0.3, 0.0, 0.3, 20.0] {
        let point = law.evaluate(vd, 0.4).unwrap();
        assert_eq!(point.current_a, 0.0);
        assert_eq!(point.drain_derivative_s, 0.0);
        assert_eq!(point.gate_derivative_s, 0.0);
    }
    let zero = law.evaluate(0.0, 1.5).unwrap();
    assert_eq!(zero.current_a, 0.0);
    close(zero.drain_derivative_s, 2.0, 0.0);
    assert_eq!(zero.gate_derivative_s, 0.0);
    let old = GateConductanceCurve::try_new(vec![(-1.0, 0.1), (4.0, 2.1)]).unwrap();
    let current = GateChannelResponse::evaluate(&old, -0.3, 1.5).unwrap();
    close(current.current_a, -0.33, 2e-15);
    close(current.drain_derivative_s, 1.1, 2e-15);
    close(current.gate_derivative_s, -0.12, 2e-15);
}

#[derive(Debug)]
struct GivenPoint {
    range: [f64; 2],
    point: GateChannelPoint,
}
impl GateChannelResponse for GivenPoint {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        self.range
    }
    fn evaluate(&self, _: f64, _: f64) -> Result<GateChannelPoint, PhsError> {
        Ok(self.point)
    }
}

#[test]
fn source_support_active_work_zero_drain_and_nonfinite_responses_refuse() {
    let law = GateChannelLaw::new(SaturatingChannel).unwrap();
    for (vd, vg) in [
        (20.01, 2.0),
        (-20.01, 2.0),
        (1.0, -1.01),
        (1.0, 4.01),
        (f64::NAN, 2.0),
        (1.0, f64::INFINITY),
    ] {
        assert!(law.evaluate(vd, vg).is_err());
    }
    for range in [
        [1.0, 1.0],
        [2.0, -1.0],
        [f64::NAN, 1.0],
        [-1.0, f64::INFINITY],
    ] {
        assert!(
            GateChannelLaw::new(GivenPoint {
                range,
                point: GateChannelPoint {
                    current_a: 0.0,
                    drain_derivative_s: 0.0,
                    gate_derivative_s: 0.0,
                }
            })
            .is_err()
        );
    }
    for (vd, point) in [
        (1.0, [f64::NAN, 0.0, 0.0]),
        (1.0, [1.0, f64::INFINITY, 0.0]),
        (1.0, [1.0, 0.0, f64::NAN]),
        (1.0, [-1.0, 1.0, 0.0]),
        (0.0, [1.0, 1.0, 0.0]),
        (f64::MAX, [2.0, 0.0, 0.0]),
        (1e-300, [1e-300, 0.0, 0.0]),
    ] {
        let result = GateChannelLaw::new(GivenPoint {
            range: [-1.0, 4.0],
            point: GateChannelPoint {
                current_a: point[0],
                drain_derivative_s: point[1],
                gate_derivative_s: point[2],
            },
        })
        .and_then(|law| law.evaluate(vd, 2.0));
        assert!(
            result.is_err(),
            "unphysical or unrepresentable source response {point:?}"
        );
    }
}

#[derive(Debug)]
struct PassiveNegativeSlope;
impl GateChannelResponse for PassiveNegativeSlope {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-1.0, 4.0]
    }
    fn evaluate(&self, v: f64, _: f64) -> Result<GateChannelPoint, PhsError> {
        Ok(GateChannelPoint {
            current_a: v / (1.0 + v * v),
            drain_derivative_s: (1.0 - v * v) / (1.0 + v * v).powi(2),
            gate_derivative_s: 0.0,
        })
    }
}

#[test]
fn passive_source_owned_negative_differential_slope_is_not_replaced_by_conductance() {
    let law = GateChannelLaw::new(PassiveNegativeSlope).unwrap();
    for vd in [-2.0, 2.0] {
        let p = law.evaluate(vd, 0.0).unwrap();
        close(p.current_a, vd / 5.0, 0.0);
        close(p.drain_derivative_s, -3.0 / 25.0, 0.0);
        let h = 1e-5;
        close(
            (law.evaluate(vd + h, 0.0).unwrap().current_a
                - law.evaluate(vd - h, 0.0).unwrap().current_a)
                / (2.0 * h),
            p.drain_derivative_s,
            2e-11,
        );
        assert!(vd * p.current_a > 0.0);
    }
}
