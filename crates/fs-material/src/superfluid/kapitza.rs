//! Authored Kapitza heat transfer across a solid/helium interface.
//!
//! `Q = area * sigma * (T_from^n - T_to^n)`, positive from the first side
//! into the second. The common phonon-radiation form has n=4; measured
//! interfaces can require other exponents and surface-dependent coefficients.
//! See L. Bottura, *Cable Stability*, CERN-2014-005, equations (28)-(29):
//! <https://cds.cern.ch/record/1974064/files/CERN-2014-005-p401.pdf>.
//!
//! This is an interface law, not bulk helium conductivity. The caller supplies
//! the coefficient, actual contact area, exponent and valid absolute-temperature
//! interval. No material-name calibration, helium phase boundary, boiling,
//! surface coating, or transient interface storage is inferred. A thermal owner
//! must debit and credit the same heat once. The associated entropy production
//! is not an additional heat source.

use super::{invalid, positive};
use crate::MaterialError;
use fs_math::det;

/// Explicit interface parameters. The coefficient uses the power-law exponent
/// shown here; it is not the differential conductance in W/(m² K).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KapitzaParameters {
    /// Actual heat-transfer area [m²], positive.
    pub area_m2: f64,
    /// Power-law coefficient sigma [W/(m² K^n)], positive.
    pub coefficient_w_m2_k_n: f64,
    /// Authored exponent n >= 1. The common quartic law has n=4.
    pub exponent: f64,
    /// Inclusive lower validity temperature for both sides [K], positive.
    pub min_temperature_k: f64,
    /// Inclusive upper validity temperature for both sides [K].
    pub max_temperature_k: f64,
}

/// Immutable, passive interface with an explicitly admitted temperature domain.
/// Evaluation is allocation-free and uses the shared deterministic math owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KapitzaInterface {
    parameters: KapitzaParameters,
    coefficient_w_k_n: f64,
    minimum_tangent_w_k: f64,
}

impl KapitzaInterface {
    /// Admit a finite positive law whose endpoint tangents and maximum heat
    /// transfer are representable. No extrapolation beyond the range is allowed.
    pub fn try_new(parameters: KapitzaParameters) -> Result<Self, MaterialError> {
        let p = parameters;
        if ![
            p.area_m2,
            p.coefficient_w_m2_k_n,
            p.min_temperature_k,
            p.max_temperature_k,
        ]
        .into_iter()
        .all(positive)
            || !p.exponent.is_finite()
            || p.exponent < 1.0
            || p.min_temperature_k >= p.max_temperature_k
        {
            return Err(invalid(
                "Kapitza law needs positive area/coefficient, n >= 1 and ordered positive temperature bounds",
            ));
        }
        let coefficient_w_k_n = p.area_m2 * p.coefficient_w_m2_k_n;
        if !positive(coefficient_w_k_n) {
            return Err(invalid(
                "Kapitza total interface coefficient is unrepresentable",
            ));
        }
        let mut law = Self {
            parameters,
            coefficient_w_k_n,
            minimum_tangent_w_k: 0.0,
        };
        law.minimum_tangent_w_k = law.tangent_from_w_k(p.min_temperature_k)?;
        law.tangent_from_w_k(p.max_temperature_k)?;
        let span = law.heat_flow_w(p.max_temperature_k, p.min_temperature_k)?;
        if !positive(span) {
            return Err(invalid(
                "Kapitza full-range heat transfer is unrepresentable",
            ));
        }
        Ok(law)
    }

    /// Original declared inputs, including area and the source validity range.
    #[must_use]
    pub const fn parameters(&self) -> KapitzaParameters {
        self.parameters
    }

    /// Inclusive [minimum, maximum] absolute-temperature validity bounds [K].
    #[must_use]
    pub const fn temperature_bounds_k(&self) -> [f64; 2] {
        [
            self.parameters.min_temperature_k,
            self.parameters.max_temperature_k,
        ]
    }

    /// Positive lower bound of each side's differential conductance throughout
    /// the declared range [W/K]. For n >= 1 it occurs at the lower endpoint.
    /// This bound is useful for monotone implicit finite-reservoir solves.
    #[must_use]
    pub const fn minimum_tangent_w_k(&self) -> f64 {
        self.minimum_tangent_w_k
    }

    /// `dQ/dT_from = -dQ/dT_to` when evaluated at the same temperature:
    /// `area * sigma * n * T^(n-1)` [W/K]. Away from equilibrium, evaluate
    /// each side separately; a single secant conductance is not either tangent.
    pub fn tangent_from_w_k(&self, temperature_k: f64) -> Result<f64, MaterialError> {
        self.admit_temperature(temperature_k)?;
        let n = self.parameters.exponent;
        let power = match n {
            1.0 => 1.0,
            2.0 => temperature_k,
            3.0 => temperature_k * temperature_k,
            4.0 => temperature_k * temperature_k * temperature_k,
            _ => det::pow(temperature_k, n - 1.0),
        };
        let tangent = (self.coefficient_w_k_n * power) * n;
        if positive(tangent) {
            Ok(tangent)
        } else {
            Err(invalid(
                "Kapitza differential conductance is unrepresentable",
            ))
        }
    }

    /// Signed total heat power [W], positive from `from_k` into `to_k`.
    /// Integer powers use factored differences; other exponents use a corrected
    /// logarithm and expm1. A one-ULP temperature difference is not subtracted
    /// from two nearly equal powers. Reversing the sides negates the same flux.
    pub fn heat_flow_w(&self, from_k: f64, to_k: f64) -> Result<f64, MaterialError> {
        self.admit_temperature(from_k)?;
        self.admit_temperature(to_k)?;
        self.heat_from_difference_w(from_k, to_k, from_k - to_k)
    }

    /// Heat power between temperatures represented as a common reference plus
    /// separate perturbations [K]. This preserves differences smaller than one
    /// ULP of the absolute temperature, as carried by small-signal wave owners.
    /// The reference must be positive finite; it need not itself be in range
    /// if both resulting temperatures are. An error-free sum also checks the
    /// source bounds, so an out-of-range offset cannot round onto a valid edge.
    pub fn heat_flow_about_reference_w(
        &self,
        reference_k: f64,
        from_perturbation_k: f64,
        to_perturbation_k: f64,
    ) -> Result<f64, MaterialError> {
        if !positive(reference_k)
            || !from_perturbation_k.is_finite()
            || !to_perturbation_k.is_finite()
        {
            return Err(invalid(
                "Kapitza reference and perturbations must be finite; reference > 0",
            ));
        }
        let from = self.temperature_about_reference(reference_k, from_perturbation_k)?;
        let to = self.temperature_about_reference(reference_k, to_perturbation_k)?;
        self.heat_from_difference_w(from, to, from_perturbation_k - to_perturbation_k)
    }

    fn heat_from_difference_w(
        &self,
        from_k: f64,
        to_k: f64,
        difference_k: f64,
    ) -> Result<f64, MaterialError> {
        if difference_k == 0.0 {
            return Ok(0.0);
        }
        let hot = from_k.max(to_k);
        let cold = from_k.min(to_k);
        let delta = difference_k.abs();
        let ratio = cold / hot;
        let n = self.parameters.exponent;
        // Divide the secant by the hot-side tangent: this factor is in [1/n,1].
        // It avoids forming the potentially much larger absolute T^n powers.
        let fraction = match n {
            1.0 => 1.0,
            2.0 => 0.5 * (1.0 + ratio),
            3.0 => (1.0 + ratio + ratio * ratio) / 3.0,
            4.0 => 0.25 * (1.0 + ratio) * (1.0 + ratio * ratio),
            _ => {
                let relative_delta = delta / hot;
                if relative_delta == 0.0 {
                    // A subnormal offset can underflow when divided by T0;
                    // the secant still tends to the finite differential law.
                    1.0
                } else if relative_delta >= 1.0 {
                    1.0 / n
                } else {
                    // log1p(-relative_delta), corrected for rounding in 1-d.
                    // If 1-d rounds to 1, the limiting logarithm is -d.
                    let one_minus = 1.0 - relative_delta;
                    let log_ratio = if one_minus == 1.0 {
                        -relative_delta
                    } else {
                        det::ln(one_minus) * (-relative_delta / (one_minus - 1.0))
                    };
                    (-det::expm1(n * log_ratio) / n) / relative_delta
                }
            }
        };
        let magnitude = (self.tangent_from_w_k(hot)? * fraction) * delta;
        if !magnitude.is_finite() || magnitude < 0.0 {
            return Err(invalid("Kapitza heat transfer is unrepresentable"));
        }
        Ok(if difference_k > 0.0 {
            magnitude
        } else {
            -magnitude
        })
    }

    /// Nonnegative irreversible interface entropy production [W/K],
    /// `Q*(1/T_to - 1/T_from)`. This is an account, not additional deposited heat.
    pub fn entropy_production_w_k(&self, from_k: f64, to_k: f64) -> Result<f64, MaterialError> {
        let power = self.heat_flow_w(from_k, to_k)?;
        let entropy = (power / from_k) * ((from_k - to_k) / to_k);
        if entropy.is_finite() && entropy >= 0.0 {
            Ok(entropy)
        } else {
            Err(invalid("Kapitza entropy production is unrepresentable"))
        }
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), MaterialError> {
        if temperature_k.is_finite()
            && temperature_k >= self.parameters.min_temperature_k
            && temperature_k <= self.parameters.max_temperature_k
        {
            Ok(())
        } else {
            Err(invalid(
                "Kapitza interface temperature is outside its declared range",
            ))
        }
    }

    fn temperature_about_reference(
        &self,
        reference_k: f64,
        perturbation_k: f64,
    ) -> Result<f64, MaterialError> {
        let temperature = reference_k + perturbation_k;
        self.admit_temperature(temperature)?;
        // TwoSum retains the rounding residual even for a sub-ULP offset at
        // an exact source boundary. Only the range check needs this low part.
        let recovered_offset = temperature - reference_k;
        let residual =
            (reference_k - (temperature - recovered_offset)) + (perturbation_k - recovered_offset);
        if (temperature == self.parameters.min_temperature_k && residual < 0.0)
            || (temperature == self.parameters.max_temperature_k && residual > 0.0)
        {
            Err(invalid(
                "Kapitza interface offset is outside its declared range",
            ))
        } else {
            Ok(temperature)
        }
    }
}
