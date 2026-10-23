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
    pub fn at_flow(self, flow: f64) -> Result<(f64, f64), WaveguideError> {
        if !flow.is_finite() || flow.abs() > self.maximum_flow_m3_s {
            return Err(WaveguideError("quadratic terminal solved flow exceeds its declared domain"));
        }
        let drop = (self.coefficient_pa_s2_m6 * flow) * flow.abs();
        let power = drop * flow;
        if !drop.is_finite() || !power.is_finite()
            || (self.coefficient_pa_s2_m6 > 0.0 && flow != 0.0 && (drop == 0.0 || power == 0.0)) {
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
