//! Analytic positive/negative quadrature of one affine volume-flow history.
//!
//! These are interval-average rates. Multiplying by the interval duration
//! gives the actual forward/reverse swept volumes of the affine interpolant.
//! Upwind enthalpy must be integrated against each sign separately: a zero
//! signed mean does not imply zero transport between unequal-temperature cells.

/// Original mean plus the two nonnegative directional interval-average rates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AffineVolumeRates {
    pub mean_m3_s: f64,
    pub forward_m3_s: f64,
    pub reverse_m3_s: f64,
}

impl AffineVolumeRates {
    /// Retain the native owner's supplied mean exactly. It must agree with
    /// the endpoint trapezoid to arithmetic roundoff. With no sign crossing,
    /// the single nonzero rate is exactly that original mean's magnitude.
    ///
    /// At a crossing, the smaller swept lobe is a triangle:
    /// `small^2 / (2 * (large + small))`. Scaled evaluation avoids squaring a
    /// physical rate or overflowing the sum of opposite endpoint magnitudes.
    /// An unrepresentable nonzero lobe refuses; no speed floor is introduced.
    pub fn try_new(initial: f64, mean: f64, final_flow: f64) -> Result<Self, &'static str> {
        if ![initial, mean, final_flow].into_iter().all(f64::is_finite) {
            return Err("finite affine endpoint and mean rates required");
        }
        let trapezoid = 0.5 * initial + 0.5 * final_flow;
        let scale = initial.abs().max(final_flow.abs()).max(mean.abs());
        if (mean - trapezoid).abs() > 64.0 * f64::EPSILON * scale {
            return Err("native average differs from the affine endpoint mean");
        }
        let crossing = initial < 0.0 && final_flow > 0.0 || initial > 0.0 && final_flow < 0.0;
        let counter = if crossing {
            let small = initial.abs().min(final_flow.abs());
            let large = initial.abs().max(final_flow.abs());
            let ratio = small / large;
            let value = (0.5 * small) * ratio / (1.0 + ratio);
            if !value.is_finite() || value <= 0.0 {
                return Err("affine reversed swept volume is unrepresentable");
            }
            value
        } else {
            if (initial >= 0.0 && final_flow >= 0.0 && mean < 0.0)
                || (initial <= 0.0 && final_flow <= 0.0 && mean > 0.0)
                || (scale > 0.0 && mean == 0.0)
            {
                return Err("affine mean lost its original flow direction");
            }
            0.0
        };
        let forward = mean.max(0.0) + counter;
        let reverse = (-mean).max(0.0) + counter;
        if !forward.is_finite() || !reverse.is_finite() {
            return Err("affine directional volume rate is unrepresentable");
        }
        Ok(Self {
            mean_m3_s: mean,
            forward_m3_s: forward,
            reverse_m3_s: reverse,
        })
    }

    /// The reverse lobe in addition to the original signed mean. Zero denotes
    /// no within-interval reversal, including startup and exact shutdown.
    pub fn counterflow_m3_s(self) -> f64 {
        self.forward_m3_s.min(self.reverse_m3_s)
    }
}

#[cfg(test)]
mod tests;
