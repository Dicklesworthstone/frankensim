//! Finite vacuum-wavelength blackbody bands and their temperature tangents.
//!
//! With `x = c2/(lambda T)`, `C = 15/pi^4`, and the original conduction
//! Stefan--Boltzmann constant, the hemispherical power is
//! `sigma T^4 C integral x^3/(exp(x)-1) dx`. At FIXED wavelength endpoints its
//! temperature derivative is the positive integral
//! `sigma T^3 C integral x^4 exp(x)/(exp(x)-1)^2 dx`. It is not `4 E_band/T`.
//!
//! Planck's vacuum radiance convention is NIST SP 250-1, section 2:
//! <https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication250-1.pdf>.
//! `c2` uses the exact defining SI h, c and k_B values listed by NIST:
//! <https://physics.nist.gov/cuu/Constants/Table/allascii.txt>.
//! Normalization deliberately retains the existing conduction owner's sigma.
//! The low-x expansion is the Bernoulli generating function, DLMF 24.2.1:
//! <https://dlmf.nist.gov/24.2.E1>.
//!
//! Work and storage are bounded independently of wavelength/temperature: a
//! degree-24 Bernoulli expansion below x=1, at most 48 positive exponential
//! tail terms above it, and a positive 16-point rule for narrow intervals that
//! must not subtract adjacent tails. Near-equilibrium differences integrate
//! the analytic temperature derivative instead of subtracting emitted powers.
//! These are floating-point evaluations, not certified interval error bounds.
//! Spectral tails are coverage diagnostics; omitted wavelengths receive no heat.

use super::super::{STEFAN_BOLTZMANN_W_M2_K4, radiation_error, require_temperature};
use crate::ConductionError;

const C2_M_K: f64 = 6.626_070_15e-34 * 299_792_458.0 / 1.380_649e-23;
const NORMALIZATION: f64 = 15.0
    / (core::f64::consts::PI
        * core::f64::consts::PI
        * core::f64::consts::PI
        * core::f64::consts::PI);

/// One finite, strictly increasing vacuum-wavelength interval in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlackbodyBand {
    lower_wavelength_m: f64,
    upper_wavelength_m: f64,
}

/// Hemispherical band emission, its fixed-endpoint tangent, and spectral coverage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlackbodyBandEmission {
    /// Integrated emitted power per physical surface area [W/m²].
    pub power_w_m2: f64,
    /// Analytic temperature derivative with both wavelengths held fixed [W/m²/K].
    pub temperature_derivative_w_m2_k: f64,
    /// Fraction of the blackbody spectrum inside the declared interval.
    pub in_band_fraction: f64,
    /// Fraction at wavelengths shorter than the lower endpoint.
    pub short_wave_fraction: f64,
    /// Fraction at wavelengths longer than the upper endpoint.
    pub long_wave_fraction: f64,
}

#[derive(Clone, Copy)]
struct Moments {
    energy: f64,
    derivative: f64,
}

impl BlackbodyBand {
    /// Admit explicit finite positive endpoints; zero/infinity do not stand for
    /// omitted-tail assumptions. Construct a broad finite window when needed.
    pub fn new(lower_wavelength_m: f64, upper_wavelength_m: f64) -> Result<Self, ConductionError> {
        if !lower_wavelength_m.is_finite()
            || lower_wavelength_m <= 0.0
            || !upper_wavelength_m.is_finite()
            || upper_wavelength_m <= lower_wavelength_m
        {
            return Err(invalid(
                "finite positive increasing vacuum wavelengths are required",
            ));
        }
        Ok(Self {
            lower_wavelength_m,
            upper_wavelength_m,
        })
    }

    /// Short-wavelength endpoint [m].
    pub const fn lower_wavelength_m(&self) -> f64 {
        self.lower_wavelength_m
    }
    /// Long-wavelength endpoint [m].
    pub const fn upper_wavelength_m(&self) -> f64 {
        self.upper_wavelength_m
    }

    /// Hemispherical blackbody power density per unit log vacuum wavelength
    /// [W/m²]. Integrating this quantity over `d(ln lambda)` reproduces the
    /// finite-band power returned by [`Self::evaluate`]. This is the shared
    /// spectral primitive for solid and participating-gas radiation owners.
    pub fn spectral_power_per_log_wavelength(
        wavelength_m: f64,
        temperature_k: f64,
    ) -> Result<f64, ConductionError> {
        require_temperature("blackbody-spectrum", "temperature", temperature_k)?;
        if !wavelength_m.is_finite() || wavelength_m <= 0.0 {
            return Err(invalid("finite positive vacuum wavelength is required"));
        }
        let (x, log_x) = dimensionless(wavelength_m, temperature_k);
        if !x.is_finite() || x <= 0.0 {
            return Err(invalid("blackbody wavelength/temperature is unrepresentable"));
        }
        let log_shape = if x <= 1.0 {
            3.0 * log_x + (x / x.exp_m1()).ln()
        } else {
            4.0 * log_x - x - (-(-x).exp_m1()).ln()
        };
        let log_power = (STEFAN_BOLTZMANN_W_M2_K4 * NORMALIZATION).ln()
            + 4.0 * temperature_k.ln()
            + log_shape;
        finite_nonnegative(log_power.exp())
    }

    /// Signed difference of blackbody power densities per unit log vacuum
    /// wavelength [W/m²]. Equal temperatures return exactly zero and nearby
    /// temperatures use `expm1` rather than subtracting two large spectra.
    pub fn spectral_difference_per_log_wavelength(
        wavelength_m: f64,
        temperature_k: f64,
        reference_temperature_k: f64,
    ) -> Result<f64, ConductionError> {
        require_temperature("blackbody-spectrum", "temperature", temperature_k)?;
        require_temperature(
            "blackbody-spectrum",
            "reference temperature",
            reference_temperature_k,
        )?;
        if !wavelength_m.is_finite() || wavelength_m <= 0.0 {
            return Err(invalid("finite positive vacuum wavelength is required"));
        }
        if temperature_k == reference_temperature_k {
            return Ok(0.0);
        }
        let (hot, cold, sign) = if temperature_k > reference_temperature_k {
            (temperature_k, reference_temperature_k, 1.0)
        } else {
            (reference_temperature_k, temperature_k, -1.0)
        };
        let hot_emission = Self::spectral_power_per_log_wavelength(wavelength_m, hot)?;
        let x_hot = C2_M_K / wavelength_m / hot;
        let x_cold = C2_M_K / wavelength_m / cold;
        if !x_hot.is_finite() || !x_cold.is_finite() || x_hot <= 0.0 || x_cold <= 0.0 {
            return Err(invalid("blackbody wavelength/temperature is unrepresentable"));
        }
        let exponent_difference = x_hot * ((hot - cold) / cold);
        let value = sign * hot_emission * (-exponent_difference).exp_m1() / (-x_cold).exp_m1();
        if value.is_finite() {
            Ok(value)
        } else {
            Err(invalid("blackbody spectral difference is unrepresentable"))
        }
    }

    /// Resolve emission at an absolute temperature. Logarithmic scaling retains
    /// representable Wien-tail powers even when the corresponding fraction
    /// underflows; no overflowing full-spectrum T^4 intermediate is necessary.
    pub fn evaluate(&self, temperature_k: f64) -> Result<BlackbodyBandEmission, ConductionError> {
        require_temperature("blackbody-band", "temperature", temperature_k)?;
        let (a, log_a) = dimensionless(self.upper_wavelength_m, temperature_k);
        let (b, log_b) = dimensionless(self.lower_wavelength_m, temperature_k);
        let moments = self.moments(a, log_a, b, log_b);
        let (power, derivative) = dimensional(moments, temperature_k)?;
        let short = if b <= 1.0 {
            -(NORMALIZATION.ln() + low_moments(b, log_b, 0.0, 1.0).energy).exp_m1()
        } else {
            (NORMALIZATION.ln() + tail_moments(b).energy).exp()
        };
        let long = if a <= 1.0 {
            (NORMALIZATION.ln() + low_moments(a, log_a, 0.0, 1.0).energy).exp()
        } else {
            -(NORMALIZATION.ln() + tail_moments(a).energy).exp_m1()
        };
        Ok(BlackbodyBandEmission {
            power_w_m2: power,
            temperature_derivative_w_m2_k: derivative,
            in_band_fraction: fraction((NORMALIZATION.ln() + moments.energy).exp())?,
            short_wave_fraction: fraction(short)?,
            long_wave_fraction: fraction(long)?,
        })
    }

    /// Signed `E_band(T) - E_band(reference)` without equilibrium cancellation.
    /// Swapping temperatures negates the identical calculation; equality is
    /// exactly zero. This contains no finite-difference temperature derivative.
    pub fn difference(
        &self,
        temperature_k: f64,
        reference_temperature_k: f64,
    ) -> Result<f64, ConductionError> {
        require_temperature("blackbody-band", "temperature", temperature_k)?;
        require_temperature(
            "blackbody-band",
            "reference temperature",
            reference_temperature_k,
        )?;
        if temperature_k == reference_temperature_k {
            return Ok(0.0);
        }
        let (cold, hot, sign) = if temperature_k > reference_temperature_k {
            (reference_temperature_k, temperature_k, 1.0)
        } else {
            (temperature_k, reference_temperature_k, -1.0)
        };
        let delta = hot - cold;
        let (a, _) = dimensionless(self.upper_wavelength_m, cold);
        let variation = (delta / hot) * (a + 5.0);
        let difference = if variation <= 1e-7 {
            // The midpoint remainder is quadratic in this very small spectral
            // variation. This also makes adjacent-float temperatures inexpensive.
            delta * self.values(cold + 0.5 * delta)?.1
        } else if variation <= 0.5 {
            let mut average = 0.0;
            for (node, weight) in GL8 {
                for direction in [-1.0, 1.0] {
                    let t = (0.5 * (1.0 + direction * node)).mul_add(delta, cold);
                    average += 0.5 * weight * self.values(t)?.1;
                }
            }
            delta * average
        } else {
            self.values(hot)?.0 - self.values(cold)?.0
        };
        Ok(sign * finite_nonnegative(difference)?)
    }

    fn values(&self, temperature: f64) -> Result<(f64, f64), ConductionError> {
        let (a, log_a) = dimensionless(self.upper_wavelength_m, temperature);
        let (b, log_b) = dimensionless(self.lower_wavelength_m, temperature);
        dimensional(self.moments(a, log_a, b, log_b), temperature)
    }

    fn moments(&self, a: f64, log_a: f64, b: f64, log_b: f64) -> Moments {
        if b <= 1.0 {
            let width_ratio =
                (self.upper_wavelength_m - self.lower_wavelength_m) / self.upper_wavelength_m;
            return low_moments(
                b,
                log_b,
                self.lower_wavelength_m / self.upper_wavelength_m,
                width_ratio,
            );
        }
        if a >= 1.0 {
            // Form the width from wavelength separation, not rounded x_b-x_a.
            // Adjacent floating-point wavelengths therefore retain positive heat.
            let width =
                a * ((self.upper_wavelength_m - self.lower_wavelength_m) / self.lower_wavelength_m);
            return high_interval(a, b, width);
        }
        let low = low_moments(1.0, 0.0, a, -log_a.exp_m1());
        let high = high_interval(1.0, b, log_b.exp_m1());
        Moments {
            energy: log_add(low.energy, high.energy),
            derivative: log_add(low.derivative, high.derivative),
        }
    }
}

fn invalid(what: &'static str) -> ConductionError {
    radiation_error(
        "blackbody-band",
        what,
        "supply finite absolute temperatures and an explicit vacuum-wavelength window",
    )
}
fn finite_nonnegative(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(invalid(
            "blackbody band power or tangent is unrepresentable",
        ))
    }
}
fn fraction(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() && (0.0..=1.0 + 32.0 * f64::EPSILON).contains(&value) {
        Ok(value.min(1.0))
    } else {
        Err(invalid("blackbody spectral coverage is unrepresentable"))
    }
}
fn dimensionless(wavelength: f64, temperature: f64) -> (f64, f64) {
    let x = (C2_M_K / wavelength) / temperature;
    if x.is_finite() && x > 0.0 {
        (x, x.ln())
    } else {
        let logarithm = C2_M_K.ln() - wavelength.ln() - temperature.ln();
        (logarithm.exp(), logarithm)
    }
}
fn dimensional(moments: Moments, temperature: f64) -> Result<(f64, f64), ConductionError> {
    let scale = (STEFAN_BOLTZMANN_W_M2_K4 * NORMALIZATION).ln();
    let log_t = temperature.ln();
    Ok((
        finite_nonnegative((scale + 4.0 * log_t + moments.energy).exp())?,
        finite_nonnegative((scale + 3.0 * log_t + moments.derivative).exp())?,
    ))
}
fn log_add(a: f64, b: f64) -> f64 {
    let larger = a.max(b);
    if larger == f64::NEG_INFINITY {
        larger
    } else {
        larger + (a.min(b) - larger).exp().ln_1p()
    }
}
fn log_subtract(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY || b == f64::NEG_INFINITY {
        a
    } else {
        a + (-(b - a).exp_m1()).ln()
    }
}

// Integral coefficients of x^3/(exp(x)-1), through the B_22 term. The
// omitted B_24 contribution below x=1 is < 6e-21 in the unnormalized integral.
const LOW: [(usize, f64); 13] = [
    (3, 1.0 / 3.0),
    (4, -1.0 / 8.0),
    (5, 1.0 / 60.0),
    (7, -1.0 / 5040.0),
    (9, 1.0 / 272160.0),
    (11, -1.0 / 13305600.0),
    (13, 1.0 / 622702080.0),
    (15, -691.0 / 19615115520000.0),
    (17, 1.0 / 1270312243200.0),
    (19, -3617.0 / 202741834014720000.0),
    (21, 43867.0 / 107290978560589824000.0),
    (23, -174611.0 / 18465726242060697600000.0),
    (25, 854513.0 / 3877802510832746496000000.0),
];

// x_low = ratio*x_high. Factor the width and x_high^2 before summing, avoiding
// both adjacent-antiderivative cancellation and small-x cubic underflow.
fn low_moments(high: f64, log_high: f64, ratio: f64, width_ratio: f64) -> Moments {
    let mut energy = 0.0;
    let mut derivative = 0.0;
    let mut sum_powers = 0.0;
    let mut ratio_power = 1.0;
    let mut coefficient = 0;
    for power in 1..=25 {
        sum_powers += ratio_power;
        ratio_power *= ratio;
        if coefficient < LOW.len() && LOW[coefficient].0 == power {
            let term = LOW[coefficient].1 * high.powi(power as i32 - 3) * sum_powers;
            energy += term;
            derivative += (4.0 - power as f64) * term;
            coefficient += 1;
        }
    }
    let scale = 3.0 * log_high + width_ratio.ln();
    Moments {
        energy: scale + energy.ln(),
        derivative: scale + derivative.ln(),
    }
}

// Positive exponential tails, with exp(-x)*x^3 and exp(-x)*x^4 factored out.
// At x>=1 the uncomputed terms after n=48 are below 1e-18 of either tail.
fn tail_moments(x: f64) -> Moments {
    if x > 8192.0 {
        // Even sigma*T^4 for the largest finite T cannot make these powers
        // representable. This guard also avoids infinity*zero arithmetic.
        return Moments {
            energy: f64::NEG_INFINITY,
            derivative: f64::NEG_INFINITY,
        };
    }
    let inverse = 1.0 / x;
    let attenuation = (-x).exp();
    let mut weight = 1.0;
    let mut energy = 0.0;
    let mut derivative = 0.0;
    for n in 1..=48 {
        let inv_n = 1.0 / f64::from(n);
        let q = inverse * inv_n;
        let e = weight * inv_n * (1.0 + q * (3.0 + q * (6.0 + 6.0 * q)));
        let d = weight * (1.0 + q * (4.0 + q * (12.0 + q * (24.0 + 24.0 * q))));
        energy += e;
        derivative += d;
        if e <= 1e-18 * energy && d <= 1e-18 * derivative {
            break;
        }
        weight *= attenuation;
    }
    let logarithm = x.ln();
    Moments {
        energy: -x + 3.0 * logarithm + energy.ln(),
        derivative: -x + 4.0 * logarithm + derivative.ln(),
    }
}

fn high_interval(low: f64, high: f64, width: f64) -> Moments {
    if width <= 1.0 {
        let mut energy = 0.0;
        let mut derivative = 0.0;
        for (node, weight) in GL16 {
            for sign in [-1.0, 1.0] {
                let offset = 0.5 * width * (1.0 + sign * node);
                let relative = 1.0 + offset / low;
                let denominator = -(-(low + offset)).exp_m1();
                let shared = 0.5 * weight * (-offset).exp() / denominator;
                energy += shared * relative.powi(3);
                derivative += shared * relative.powi(4) / denominator;
            }
        }
        let scale = -low + width.ln();
        let logarithm = low.ln();
        Moments {
            energy: scale + 3.0 * logarithm + energy.ln(),
            derivative: scale + 4.0 * logarithm + derivative.ln(),
        }
    } else {
        let a = tail_moments(low);
        let b = tail_moments(high);
        Moments {
            energy: log_subtract(a.energy, b.energy),
            derivative: log_subtract(a.derivative, b.derivative),
        }
    }
}

const GL8: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];
const GL16: [(f64, f64); 8] = [
    (0.09501250983763744, 0.18945061045506850),
    (0.28160355077925891, 0.18260341504492359),
    (0.45801677765722739, 0.16915651939500254),
    (0.61787624440264375, 0.14959598881657673),
    (0.75540440835500303, 0.12462897125553387),
    (0.86563120238783174, 0.09515851168249278),
    (0.94457502307323258, 0.06225352393864789),
    (0.98940093499164993, 0.02715245941175409),
];

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    // Independent physical-wavelength composite Simpson integration. The
    // dimensional NIST radiance formula, not production cumulative moments,
    // owns this oracle; integration is in log(lambda) for the broad windows.
    fn oracle(band: BlackbodyBand, temperature: f64) -> (f64, f64) {
        let h: f64 = 6.62607015e-34;
        let c: f64 = 299792458.0;
        let k: f64 = 1.380649e-23;
        let sigma_si =
            2.0 * core::f64::consts::PI.powi(5) * k.powi(4) / (15.0 * h.powi(3) * c.powi(2));
        let normalization = STEFAN_BOLTZMANN_W_M2_K4 / sigma_si;
        let lo = band.lower_wavelength_m.ln();
        let length = band.upper_wavelength_m.ln() - lo;
        let panels = 16384;
        let mut energy = 0.0;
        let mut derivative = 0.0;
        for i in 0..=panels {
            let wavelength = (lo + length * i as f64 / panels as f64).exp();
            let x = h * c / (k * wavelength * temperature);
            let value = if x > 700.0 {
                0.0
            } else {
                normalization * 2.0 * core::f64::consts::PI * h * c * c
                    / (wavelength.powi(4) * x.exp_m1())
            };
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            energy += weight * value;
            derivative += weight * value * x / (temperature * -(-x).exp_m1());
        }
        (
            energy * length / (3.0 * panels as f64),
            derivative * length / (3.0 * panels as f64),
        )
    }

    #[test]
    fn band_power_and_fixed_wavelength_tangent_match_independent_planck_quadrature() {
        // Independently generated 75-digit Decimal integration with 64/128
        // recomputed Legendre nodes, using exact SI h,c,k and original sigma.
        for (lo, hi, temperature, power, tangent) in [
            (
                5e-6,
                50e-6,
                300.0,
                4.39129759958399795e2,
                5.84277160565876642,
            ),
            (
                1e-7,
                1e-3,
                1200.0,
                1.17580873596434766e5,
                3.91936271172470134e2,
            ),
            (
                2e-7,
                2e-6,
                3000.0,
                3.38866878961543647e6,
                5.29770109715799902e3,
            ),
            (
                1e-3,
                1e-2,
                50.0,
                3.88035434980384928e-4,
                8.62441903023964855e-6,
            ),
            (
                1e-7,
                1.01e-7,
                4000.0,
                1.04766192433856520e-5,
                9.37177839667846831e-8,
            ),
        ] {
            let band = BlackbodyBand::new(lo, hi).unwrap();
            let result = band.evaluate(temperature).unwrap();
            let reference = oracle(band, temperature);
            close(result.power_w_m2, reference.0, 3e-12);
            close(result.temperature_derivative_w_m2_k, reference.1, 3e-12);
            close(result.power_w_m2, power, 6e-14);
            close(result.temperature_derivative_w_m2_k, tangent, 6e-14);
            close(
                result.in_band_fraction + result.short_wave_fraction + result.long_wave_fraction,
                1.0,
                3e-14,
            );
        }
    }

    #[test]
    fn complete_window_recovers_original_sigma_and_band_partition_is_additive() {
        for temperature in [1.0_f64, 300.0, 5000.0] {
            let band = BlackbodyBand::new(
                C2_M_K / (1000.0 * temperature),
                C2_M_K / (1e-6 * temperature),
            )
            .unwrap();
            let result = band.evaluate(temperature).unwrap();
            close(
                result.power_w_m2,
                STEFAN_BOLTZMANN_W_M2_K4 * temperature.powi(4),
                3e-14,
            );
            close(
                result.temperature_derivative_w_m2_k,
                4.0 * STEFAN_BOLTZMANN_W_M2_K4 * temperature.powi(3),
                3e-14,
            );
        }
        let endpoints = [1e-7, 1e-6, 2.5e-6, 10e-6, 1e-4];
        for temperature in [120.0, 900.0, 3500.0] {
            let total = BlackbodyBand::new(endpoints[0], endpoints[4])
                .unwrap()
                .evaluate(temperature)
                .unwrap();
            let parts: Vec<_> = endpoints
                .windows(2)
                .map(|v| {
                    BlackbodyBand::new(v[0], v[1])
                        .unwrap()
                        .evaluate(temperature)
                        .unwrap()
                })
                .collect();
            close(
                parts.iter().map(|p| p.power_w_m2).sum(),
                total.power_w_m2,
                4e-14,
            );
            close(
                parts.iter().map(|p| p.temperature_derivative_w_m2_k).sum(),
                total.temperature_derivative_w_m2_k,
                4e-14,
            );
        }
    }

    #[test]
    fn narrow_bands_and_near_equilibrium_heat_keep_their_sign_and_magnitude() {
        let lo: f64 = 4e-6;
        let hi = f64::from_bits(lo.to_bits() + 4);
        let temperature = 700.0;
        let band = BlackbodyBand::new(lo, hi).unwrap();
        let x = C2_M_K / (lo * temperature);
        let expected = STEFAN_BOLTZMANN_W_M2_K4 * temperature.powi(4) * NORMALIZATION * x.powi(4)
            / x.exp_m1()
            * ((hi - lo) / lo);
        close(
            band.evaluate(temperature).unwrap().power_w_m2,
            expected,
            5e-14,
        );
        for band in [
            band,
            BlackbodyBand::new(1e-6, 20e-6).unwrap(),
            BlackbodyBand::new(1e-3, 2e-3).unwrap(),
        ] {
            let hot = f64::from_bits(temperature.to_bits() + 8);
            let reference = band
                .evaluate(temperature)
                .unwrap()
                .temperature_derivative_w_m2_k
                * (hot - temperature);
            let heat = band.difference(hot, temperature).unwrap();
            assert!(heat > 0.0);
            close(heat, reference, 8e-14);
            assert_eq!(heat, -band.difference(temperature, hot).unwrap());
            assert_eq!(band.difference(temperature, temperature).unwrap(), 0.0);
        }
    }

    #[test]
    fn wien_fraction_underflow_does_not_erase_representable_dimensional_heat() {
        let temperature = 1e70;
        let band = BlackbodyBand::new(
            C2_M_K / (801.0 * temperature),
            C2_M_K / (800.0 * temperature),
        )
        .unwrap();
        let result = band.evaluate(temperature).unwrap();
        assert_eq!(result.in_band_fraction, 0.0);
        assert!(result.power_w_m2 > 0.0 && result.power_w_m2.is_finite());
        assert!(result.temperature_derivative_w_m2_k > 0.0);
    }

    #[test]
    fn physical_wavelength_temperature_scaling_and_invalid_inputs() {
        let base = BlackbodyBand::new(1e-6, 30e-6)
            .unwrap()
            .evaluate(900.0)
            .unwrap();
        let scaled = BlackbodyBand::new(1e-5, 300e-6)
            .unwrap()
            .evaluate(90.0)
            .unwrap();
        close(base.in_band_fraction, scaled.in_band_fraction, 2e-14);
        close(base.power_w_m2 * 1e-4, scaled.power_w_m2, 3e-14);
        close(
            base.temperature_derivative_w_m2_k * 1e-3,
            scaled.temperature_derivative_w_m2_k,
            3e-14,
        );
        for (a, b) in [
            (0.0, 1.0),
            (-1.0, 1.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (f64::NAN, 1.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(BlackbodyBand::new(a, b).is_err());
        }
        let band = BlackbodyBand::new(1e-6, 10e-6).unwrap();
        for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(band.evaluate(temperature).is_err());
            assert!(band.difference(300.0, temperature).is_err());
        }
    }
}
