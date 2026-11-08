//! Passive series loss `dp = K U |U|`, with no lagged-flow approximation.
//!
//! The coefficient is supplied by the physical geometry/material owner. This
//! numerical primitive does not infer an edge shape, discharge coefficient or
//! acoustic validity domain. It solves `U + K/R U|U| = U_linear` against the
//! SAME implicit load history and reports `K |U|^3 dt` as irreversible loss.

use crate::waveguide::WaveguideError;

/// Immutable dimensional loss coefficient and admitted magnitude of solved flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticLossSpec {
    coefficient_pa_s2_m6: f64,
    maximum_flow_m3_s: f64,
}

impl QuadraticLossSpec {
    /// `coefficient` is nonnegative [Pa s²/m⁶]; `maximum_flow` is positive [m³/s].
    /// Zero coefficient is an exact linear comparison with the same flow limit.
    pub fn new(coefficient: f64, maximum_flow: f64) -> Result<Self, WaveguideError> {
        if !coefficient.is_finite() || coefficient < 0.0
            || !maximum_flow.is_finite() || maximum_flow <= 0.0 {
            return Err(WaveguideError("quadratic loss needs finite K>=0 and a positive flow limit"));
        }
        Ok(Self { coefficient_pa_s2_m6: coefficient, maximum_flow_m3_s: maximum_flow })
    }

    #[must_use]
    pub const fn coefficient_pa_s2_m6(self) -> f64 { self.coefficient_pa_s2_m6 }
    #[must_use]
    pub const fn maximum_flow_m3_s(self) -> f64 { self.maximum_flow_m3_s }

    /// Signed drop and nonnegative dissipated power at the actually solved flow.
    ///
    /// Past the flow domain, only a non-finite drop or power is refused. With a positive `K`, a
    /// finite nonzero flow so small that `K U|U|` or `K |U|^3` falls below the
    /// smallest subnormal returns that product as the zero it rounds to: zero
    /// is its true value to double precision, and the work balance stays
    /// closed (the power is the drop times the flow either way). Every other
    /// flow returns exactly the arithmetic it always did.
    pub fn at_flow(self, flow: f64) -> Result<(f64, f64), WaveguideError> {
        if !flow.is_finite() || flow.abs() > self.maximum_flow_m3_s {
            return Err(WaveguideError("quadratic terminal solved flow exceeds its declared domain"));
        }
        let drop = (self.coefficient_pa_s2_m6 * flow) * flow.abs();
        let power = drop * flow;
        if !drop.is_finite() || !power.is_finite() {
            return Err(WaveguideError("quadratic terminal loss is not representable"));
        }
        Ok((drop, power))
    }

    pub(crate) fn solve(self, linear_flow: f64, inverse_resistance: f64)
        -> Result<f64, WaveguideError> {
        // Rationalized root: no subtraction of nearly equal square roots, no
        // division by the flow, and exact odd symmetry across reversal.
        let k = self.coefficient_pa_s2_m6 * inverse_resistance;
        let discriminant = k * linear_flow.abs();
        if !linear_flow.is_finite() || !inverse_resistance.is_finite()
            || inverse_resistance <= 0.0 || !k.is_finite() || !discriminant.is_finite() {
            return Err(WaveguideError("quadratic terminal discrete equation is not representable"));
        }
        let flow = if self.coefficient_pa_s2_m6 == 0.0 { linear_flow } else {
            linear_flow / (0.5 + fs_math::det::sqrt(0.25 + discriminant))
        };
        self.at_flow(flow)?;
        Ok(flow)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNREPRESENTABLE: WaveguideError = WaveguideError("quadratic terminal loss is not representable");
    const DOMAIN: WaveguideError = WaveguideError("quadratic terminal solved flow exceeds its declared domain");

    /// `at_flow` as it was before underflow became a zero instead of a refusal, kept as the reference.
    fn previous(spec: QuadraticLossSpec, flow: f64) -> Result<(f64, f64), WaveguideError> {
        if !flow.is_finite() || flow.abs() > spec.maximum_flow_m3_s {
            return Err(DOMAIN);
        }
        let drop = (spec.coefficient_pa_s2_m6 * flow) * flow.abs();
        let power = drop * flow;
        if !drop.is_finite() || !power.is_finite()
            || (spec.coefficient_pa_s2_m6 > 0.0 && flow != 0.0 && (drop == 0.0 || power == 0.0)) {
            return Err(UNREPRESENTABLE);
        }
        Ok((drop, power))
    }

    /// Both signs of every decade from the smallest subnormal up to `limit` (five exactly parsed
    /// mantissas each), the smallest normal, the limit itself and both zeros.
    fn flows(limit: f64) -> Vec<f64> {
        let mut magnitudes = vec![f64::from_bits(1), f64::MIN_POSITIVE, limit];
        for exponent in -324..=2 {
            for mantissa in ["1", "1.7", "3.1", "5.5", "9.3"] {
                let u: f64 = format!("{mantissa}e{exponent}").parse().unwrap();
                if u > 0.0 && u <= limit { magnitudes.push(u); }
            }
        }
        let mut all = vec![0.0, -0.0];
        for u in magnitudes { all.extend([u, -u]); }
        all
    }

    #[test]
    fn every_flow_the_old_arithmetic_admitted_returns_the_same_bits() {
        let (mut same, mut underflow) = (0_usize, 0_usize);
        for coefficient in [0.0, 1.0, 1.0e3, 1.83e6, 1.0e9, 1.0e12] {
            for limit in [1.0e-3, 0.1, 10.0] {
                let spec = QuadraticLossSpec::new(coefficient, limit).unwrap();
                for flow in flows(limit) {
                    let now = spec.at_flow(flow);
                    match previous(spec, flow) {
                        Ok((drop, power)) => {
                            let (d, p) = now.unwrap();
                            assert_eq!((d.to_bits(), p.to_bits()), (drop.to_bits(), power.to_bits()),
                                "K {coefficient:e} U {flow:e}");
                            same += 1;
                        },
                        Err(UNREPRESENTABLE) => {
                            // Finite K, U and limit: only an underflowed drop or power reaches here.
                            let (d, p) = now.unwrap();
                            assert!(coefficient > 0.0 && flow != 0.0, "K {coefficient:e} U {flow:e}");
                            assert_eq!(d.to_bits(), ((coefficient * flow) * flow.abs()).to_bits(), "K {coefficient:e} U {flow:e}");
                            assert_eq!(p.to_bits(), 0, "power is +0, K {coefficient:e} U {flow:e}");
                            underflow += 1;
                        },
                        Err(other) => assert_eq!(now, Err(other), "K {coefficient:e} U {flow:e}"),
                    }
                }
            }
        }
        assert!(same > 1000 && underflow > 1000, "same {same} underflow {underflow}");
    }

    #[test]
    fn a_flow_near_1e_minus_200_gives_zero_loss_and_power_instead_of_a_refusal() {
        // An open pipe end of radius 13.5 mm in air with C_c 1: K = rho / (2 A^2), Mach 0.3 limit.
        let area = core::f64::consts::PI * 0.0135 * 0.0135;
        let spec = QuadraticLossSpec::new(1.2 / (2.0 * area * area), 0.3 * 343.0 * area).unwrap();
        for flow in [1.0e-200, -1.0e-200, 3.7e-201, -3.7e-201] {
            assert_eq!(previous(spec, flow), Err(UNREPRESENTABLE), "U {flow:e}");
            let (drop, power) = spec.at_flow(flow).unwrap();
            assert_eq!(drop.abs().to_bits(), 0, "U {flow:e}");
            assert_eq!(drop.is_sign_negative(), flow < 0.0, "the zero drop keeps the flow's sign, U {flow:e}");
            assert_eq!(power.to_bits(), 0, "U {flow:e}");
            // The implicit solve that feeds this flow now admits it, unchanged.
            assert_eq!(spec.solve(flow, 1.0e-3).unwrap().to_bits(), flow.to_bits(), "U {flow:e}");
        }
    }

    #[test]
    fn non_finite_values_are_still_refused() {
        let spec = QuadraticLossSpec::new(1.0e300, 1.0e10).unwrap();
        assert_eq!(spec.at_flow(1.0e10), Err(UNREPRESENTABLE));
        assert_eq!(spec.at_flow(-1.0e10), Err(UNREPRESENTABLE));
        assert_eq!(spec.at_flow(f64::NAN), Err(DOMAIN));
        assert_eq!(spec.at_flow(f64::INFINITY), Err(DOMAIN));
        assert_eq!(spec.at_flow(2.0e10), Err(DOMAIN));
    }
}
