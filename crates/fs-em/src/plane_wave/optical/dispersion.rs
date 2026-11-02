//! Homogeneous visible dispersion and intensity-attenuation conversions.
//!
//! Source/quantity admission stays with the material bridge. Wavelengths are
//! vacuum wavelengths; alpha is power absorption per physical metre. Therefore
//! N=n-i*k for exp(+i omega t), k=alpha*lambda_vac/(4*pi), and
//! epsilon=N²/(mu*c0²). Neither interface reflection nor scattering is alpha.
//! The Sellmeier equation is SCHOTT TIE-29, section 2.3, expressed in SI metres:
//! <https://www.schott.com/en-us/products/optical-glass>.

use fs_math::{c64::C64, dd::Dd, det};

pub(super) const LIGHT_SPEED_M_S: f64 = 299_792_458.0;
pub(super) const FREQUENCY_WAVELENGTH_NM: f64 = core::f64::consts::TAU * LIGHT_SPEED_M_S * 1e9;
pub(super) const VISIBLE_BAND_NM: [f64; 2] = [380.0, 780.0];

pub(super) fn visible_wavelength_nm(omega: f64) -> Option<f64> {
    let [low, high] = VISIBLE_BAND_NM;
    if !omega.is_finite()
        || omega < FREQUENCY_WAVELENGTH_NM / high
        || omega > FREQUENCY_WAVELENGTH_NM / low
    {
        return None;
    }
    // Admission precedes rounding correction, so one adjacent frequency
    // outside the supported band cannot be clamped into a usable source.
    Some((FREQUENCY_WAVELENGTH_NM / omega).clamp(low, high))
}

pub(super) fn cauchy_index([a, b_m2, c_m4]: [f64; 3], wavelength_nm: f64) -> Option<f64> {
    if ![a, b_m2, c_m4].iter().all(|x| x.is_finite() && *x >= 0.0)
        || a == 0.0
        || !wavelength_nm.is_finite()
        || wavelength_nm <= 0.0
    {
        return None;
    }
    let wavelength_m = wavelength_nm * 1e-9;
    let square_m2 = wavelength_m * wavelength_m;
    let n = a + b_m2 / square_m2 + (c_m4 / square_m2) / square_m2;
    (n.is_finite() && n > 0.0).then_some(n)
}

pub(super) fn sellmeier_index(
    strengths: [f64; 3],
    resonances_m2: [f64; 3],
    wavelength_nm: f64,
) -> Option<f64> {
    if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
        return None;
    }
    let wavelength_m = wavelength_nm * 1e-9;
    let square_m2 = wavelength_m * wavelength_m;
    let mut n_squared = Dd::from_f64(1.0);
    for (b, c) in strengths.into_iter().zip(resonances_m2) {
        if !b.is_finite() || b < 0.0 || !c.is_finite() || c < 0.0 {
            return None;
        }
        if b != 0.0 {
            if square_m2 == c {
                return None;
            }
            let term = b * (square_m2 / (square_m2 - c));
            if !term.is_finite() {
                return None;
            }
            n_squared = n_squared + Dd::from_f64(term);
        }
    }
    let n_squared = n_squared.to_f64();
    (n_squared.is_finite() && n_squared > 0.0).then(|| det::sqrt(n_squared))
}

/// A pole between source knots must not hide behind finite sampled indices.
/// Positive strengths make n² monotone between poles; both band endpoints
/// and exclusion of every active pole admit the complete visible interval.
pub(super) fn admits_sellmeier_band(strengths: [f64; 3], resonances_m2: [f64; 3]) -> bool {
    let bounds = VISIBLE_BAND_NM.map(|x| {
        let m = x * 1e-9;
        m * m
    });
    !strengths
        .iter()
        .zip(resonances_m2)
        .any(|(&b, c)| b != 0.0 && c >= bounds[0] && c <= bounds[1])
        && VISIBLE_BAND_NM
            .into_iter()
            .all(|wavelength| sellmeier_index(strengths, resonances_m2, wavelength).is_some())
}

pub(super) fn attenuating_index(n: f64, alpha_per_m: f64, wavelength_nm: f64) -> Option<C64> {
    if !n.is_finite()
        || n <= 0.0
        || !alpha_per_m.is_finite()
        || alpha_per_m < 0.0
        || !wavelength_nm.is_finite()
        || wavelength_nm <= 0.0
    {
        return None;
    }
    let k = alpha_per_m * (wavelength_nm * 1e-9 / (4.0 * core::f64::consts::PI));
    if !k.is_finite() || (alpha_per_m > 0.0 && k == 0.0) {
        return None;
    }
    Some(C64::new(n, -k))
}

pub(super) fn permittivity(index: C64, index_to_sqrt_permittivity: f64) -> Option<C64> {
    // Scale before squaring; the factored difference retains near-zero real
    // epsilon without subtracting two separately rounded squares.
    let n = index.re * index_to_sqrt_permittivity;
    let k = -index.im * index_to_sqrt_permittivity;
    let epsilon = C64::new((n - k) * (n + k), -2.0 * n * k);
    if !epsilon.re.is_finite()
        || !epsilon.im.is_finite()
        || epsilon == C64::ZERO
        || epsilon.im > 0.0
        || (index.im < 0.0 && epsilon.im == 0.0)
    {
        return None;
    }
    Some(epsilon)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_absorption_has_exactly_one_factor_of_two_in_field_decay() {
        for wavelength in [380.0, 532.1, 780.0] {
            for alpha in [0.0, 1e-9, 1e3, 2e7] {
                let index = attenuating_index(1.7, alpha, wavelength).unwrap();
                let q = index.scale(core::f64::consts::TAU / (wavelength * 1e-9));
                assert!(((-2.0 * q.im) - alpha).abs() <= 3e-15 * alpha);
                let d = 1e-8;
                let amplitude = det::exp(q.im * d);
                assert!((amplitude * amplitude - det::exp(-alpha * d)).abs() < 4e-16);
                let mu = 1.3e-6;
                let epsilon = permittivity(index, 1.0 / (LIGHT_SPEED_M_S * det::sqrt(mu))).unwrap();
                let expected =
                    (index * index).scale(1.0 / (mu * LIGHT_SPEED_M_S * LIGHT_SPEED_M_S));
                assert!((epsilon - expected).abs() < 5e-15 * expected.abs());
            }
        }
    }

    #[test]
    fn source_band_rejects_adjacent_outside_frequencies_and_hidden_sellmeier_poles() {
        let lo = FREQUENCY_WAVELENGTH_NM / 780.0;
        let hi = FREQUENCY_WAVELENGTH_NM / 380.0;
        assert_eq!(visible_wavelength_nm(lo), Some(780.0));
        assert_eq!(visible_wavelength_nm(hi), Some(380.0));
        for omega in [lo.next_down(), hi.next_up(), 0.0, f64::NAN, f64::INFINITY] {
            assert!(visible_wavelength_nm(omega).is_none());
        }
        let hidden_pole = (543.21e-9_f64).powi(2);
        assert!(!admits_sellmeier_band(
            [0.5, 0.0, 0.0],
            [hidden_pole, 0.0, 0.0]
        ));
        assert!(admits_sellmeier_band([0.0; 3], [hidden_pole, 0.0, 0.0]));
        assert!(attenuating_index(1.5, -1.0, 550.0).is_none());
        assert!(attenuating_index(1.5, f64::from_bits(1), 550.0).is_none());
    }

    #[test]
    fn cauchy_and_sellmeier_use_metres_and_admit_complete_transparent_band() {
        let a = 1.4;
        let b = 5e-15;
        let c = 2e-28;
        for wavelength in [380.0, 550.0, 780.0] {
            let micrometres: f64 = wavelength / 1000.0;
            let reference = a + 0.005 / micrometres.powi(2) + 0.0002 / micrometres.powi(4);
            assert!((cauchy_index([a, b, c], wavelength).unwrap() - reference).abs() < 4e-16);
        }
        let strengths = [0.8, 0.2, 0.5];
        let resonances = [0.01e-12, 0.04e-12, 81e-12];
        assert!(admits_sellmeier_band(strengths, resonances));
        for wavelength in [380.0, 550.0, 780.0] {
            let u2 = (wavelength / 1000.0_f64).powi(2);
            let reference = det::sqrt(
                1.0 + 0.8 * u2 / (u2 - 0.01) + 0.2 * u2 / (u2 - 0.04) + 0.5 * u2 / (u2 - 81.0),
            );
            assert!(
                (sellmeier_index(strengths, resonances, wavelength).unwrap() - reference).abs()
                    < 4e-16
            );
        }
        assert!(!admits_sellmeier_band([10.0, 0.0, 0.0], [1e-12, 0.0, 0.0]));
    }
}
