use super::{FusionErrorV1, FusionReactionV1, FusionReactivityV1, KELVIN_PER_KEV};
use fs_math::det::{exp, pow, sqrt};
use fs_qty::Temperature;

// Literal table VII from Bosch & Hale, Nucl. Fusion 32 (1992) 611–631,
// equation (14), printed p.625. BG and m_r c² retain the paper's keV basis.
// C1 produces cm³/s; convert once after evaluation. Zero coefficients are
// published zeros, not missing observations. No fitted parameter is inferred.
struct Coefficients {
    bg: f64,
    mrc2: f64,
    c: [f64; 7],
}

impl FusionReactionV1 {
    fn coefficients(self) -> Coefficients {
        match self {
            Self::DeuteriumTritium => Coefficients {
                bg: 34.3827,
                mrc2: 1_124_656.0,
                c: [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,
                    1.35000e-2, -1.06750e-4, 1.36600e-5],
            },
            Self::DeuteriumDeuteriumNeutron => Coefficients {
                bg: 31.3970,
                mrc2: 937_814.0,
                c: [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0],
            },
            Self::DeuteriumDeuteriumProton => Coefficients {
                bg: 31.3970,
                mrc2: 937_814.0,
                c: [5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0],
            },
            Self::DeuteriumHelium3 => Coefficients {
                bg: 68.7508,
                mrc2: 1_124_572.0,
                c: [5.51036e-10, 6.41918e-3, -2.02896e-3, -1.91080e-5,
                    1.35776e-4, 0.0, 0.0],
            },
        }
    }
}

pub(super) fn evaluate(
    reaction: FusionReactionV1,
    ion_temperature: Temperature,
) -> Result<FusionReactivityV1, FusionErrorV1> {
    let temperature_k = ion_temperature.value();
    if !temperature_k.is_finite() || temperature_k <= 0.0 {
        return Err(FusionErrorV1::InvalidTemperature);
    }
    let [lower, upper] = reaction.temperature_range();
    // Compare in the API's Kelvin basis, so converted source endpoints remain
    // admitted even if division by KELVIN_PER_KEV rounds one keV ULP outward.
    if temperature_k < lower.value() || temperature_k > upper.value() {
        return Err(FusionErrorV1::OutsideTemperatureDomain {
            reaction, temperature_k, range_k: [lower.value(), upper.value()],
        });
    }
    let t = temperature_k / KELVIN_PER_KEV;
    let Coefficients { bg, mrc2, c: [c1, c2, c3, c4, c5, c6, c7] } = reaction.coefficients();
    let numerator = t * c6.mul_add(t, c4).mul_add(t, c2);
    let denominator = c7.mul_add(t, c5).mul_add(t, c3).mul_add(t, 1.0);
    let correction = 1.0 - numerator / denominator;
    let theta = t / correction;
    let xi = pow(bg * bg / (4.0 * theta), 1.0 / 3.0);
    let value_m3_s = 1.0e-6 * c1 * theta * sqrt(xi / (mrc2 * t * t * t)) * exp(-3.0 * xi);

    // A=T(C2+T(C4+T C6)), B=1+T(C3+T(C5+T C7)), theta=T/(1-A/B).
    // d ln(theta)/dT = 1/T + (A'B-AB')/[B²(1-A/B)].
    // Since xi ∝ theta^(-1/3), d ln(rate)/dT = (5/6+xi)theta'/theta-3/(2T).
    let numerator_derivative = (3.0 * c6).mul_add(t, 2.0 * c4).mul_add(t, c2);
    let denominator_derivative = (3.0 * c7).mul_add(t, 2.0 * c5).mul_add(t, c3);
    let d_ln_theta = 1.0 / t
        + numerator_derivative.mul_add(denominator, -numerator * denominator_derivative)
            / (denominator * denominator * correction);
    let d_ln_rate_kev = (5.0 / 6.0 + xi).mul_add(d_ln_theta, -1.5 / t);
    let logarithmic_temperature_derivative = t * d_ln_rate_kev;
    let derivative_m3_s_k = value_m3_s * d_ln_rate_kev / KELVIN_PER_KEV;
    if !value_m3_s.is_finite() || value_m3_s <= 0.0
        || !derivative_m3_s_k.is_finite()
        || !logarithmic_temperature_derivative.is_finite()
        || logarithmic_temperature_derivative <= -1.5
    {
        return Err(FusionErrorV1::Arithmetic);
    }
    Ok(FusionReactivityV1 {
        ion_temperature, value_m3_s, derivative_m3_s_k, logarithmic_temperature_derivative,
    })
}
