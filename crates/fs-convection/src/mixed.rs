//! Manning--Qureshi's laminar aiding isothermal vertical-plate average.
//!
//! Journal of Heat Transfer 116(3), 777--780 (1994), Eqs. (4), (8), (10),
//! (11), doi:10.1115/1.2910940. Admission and directional semantics belong
//! to the card and dimensional owner; this module receives admitted magnitudes.

use fs_math::det;

pub(super) struct AidingVerticalPlate {
    reynolds: f64,
    rayleigh: f64,
    prandtl: f64,
    forced: f64,
    natural: f64,
    exponent: f64,
    forced_correction: f64,
    natural_correction: f64,
}

impl AidingVerticalPlate {
    pub(super) fn new(reynolds: f64, rayleigh: f64, prandtl: f64) -> Self {
        let forced_correction = 0.100 * det::pow(prandtl, -0.75);
        let natural_correction = 0.670 * det::pow(prandtl, -9.0 / 16.0);
        let forced = 2.0 * 0.339 * det::sqrt(reynolds) * det::pow(prandtl, 1.0 / 3.0)
            / det::pow(1.0 + forced_correction, 2.0 / 9.0);
        // Cancel Pr^(1/4) (Ra/Pr)^(1/4) algebraically. Neither Gr nor
        // xi = Gr/Re^2 is formed: both would spoil the zero-flow limits and
        // can underflow/overflow even when the desired Nu remains finite.
        let natural = (4.0 * 0.503 / 3.0) * det::pow(rayleigh, 0.25)
            / det::pow(1.0 + natural_correction, 4.0 / 9.0);
        Self {
            reynolds,
            rayleigh,
            prandtl,
            forced,
            natural,
            exponent: 3.5 * det::pow(prandtl, 0.075),
            forced_correction,
            natural_correction,
        }
    }

    /// Scale before exponentiation, so both powers lie in [0,1]. At least
    /// one is exactly one because the admission gate refuses both drivers zero.
    fn normalized_powers(&self) -> (f64, f64, f64) {
        let scale = self.forced.max(self.natural);
        (
            scale,
            det::pow(self.forced / scale, self.exponent),
            det::pow(self.natural / scale, self.exponent),
        )
    }

    pub(super) fn nusselt(&self) -> f64 {
        // Preserve the exact component at either pure-convection limit.
        if self.forced == 0.0 {
            return self.natural;
        }
        if self.natural == 0.0 {
            return self.forced;
        }
        let (scale, forced_power, natural_power) = self.normalized_powers();
        scale * det::pow(forced_power + natural_power, 1.0 / self.exponent)
    }

    pub(super) fn elasticity(&self, axis: &str) -> f64 {
        let (_, forced_power, natural_power) = self.normalized_powers();
        let sum = forced_power + natural_power;
        let forced_weight = forced_power / sum;
        let natural_weight = natural_power / sum;
        match axis {
            "Re" => 0.5 * forced_weight,
            "Ra" => 0.25 * natural_weight,
            "Pr" => {
                let forced =
                    1.0 / 3.0 + self.forced_correction / (6.0 * (1.0 + self.forced_correction));
                let natural = self.natural_correction / (4.0 * (1.0 + self.natural_correction));
                // n varies with Pr. Differentiating the power sum adds
                // 0.075/n * sum(w_i ln w_i), with 0 ln 0 defined by its limit.
                // Omitting this term gives the wrong Pr derivative in the
                // mixed regime, even when both component derivatives are right.
                let weighted_log = |weight: f64| {
                    if weight == 0.0 {
                        0.0
                    } else {
                        weight * det::ln(weight)
                    }
                };
                forced_weight * forced
                    + natural_weight * natural
                    + 0.075 / self.exponent
                        * (weighted_log(forced_weight) + weighted_log(natural_weight))
            }
            _ => 0.0,
        }
    }

    /// Absolute derivative at an interior input. Log arithmetic prevents a
    /// tiny component's power weight from rounding to zero before division
    /// by a tiny driver can restore a representable derivative.
    pub(super) fn partial(&self, axis: &str) -> f64 {
        let nu = self.nusselt();
        let (component, input, power) = match axis {
            "Re" => (self.forced, self.reynolds, 0.5),
            "Ra" => (self.natural, self.rayleigh, 0.25),
            "Pr" => return nu * self.elasticity(axis) / self.prandtl,
            _ => return 0.0,
        };
        let (scale, forced_power, natural_power) = self.normalized_powers();
        let log_weight =
            self.exponent * det::ln(component / scale) - det::ln(forced_power + natural_power);
        det::exp(det::ln(nu) + det::ln(power) + log_weight - det::ln(input))
    }
}
