//! Analytic frequency jets of the original finite-band optical laws.
//!
//! Wavelength interpolation is differentiated in its ORIGINAL coordinate.
//! For k_wave=omega*N/c0 and lambda_vac proportional to 1/omega,
//! c0*k_wave' = N-lambda*N_lambda and
//! c0*omega*k_wave'' = lambda²*N_lambda_lambda.
//! These combinations preserve exactly zero GVD of a linear n/k segment,
//! and zero attenuation dispersion of constant INTENSITY absorption.
//! Endpoints have no two-sided source neighbourhood. Interior table knots
//! require matching slopes; no smoothing, extrapolation or pole fit is added.

use super::{
    OpticalConductorResponse, OpticalMaterialResponse, ReflectionError,
    SpectralOpticalResponse, SpectralOpticalSource, VisibleDielectricAbsorption,
    VisibleOpticalStatePoint, WaterOpticalResponse, VACUUM_SPEED_M_S,
    FREQUENCY_WAVELENGTH_NM, VISIBLE_DIELECTRIC_WAVELENGTHS_NM, dispersion,
};
use fs_math::{c64::C64, dd::Dd};

type Result<T> = core::result::Result<T, ReflectionError>;

pub(in crate::plane_wave) struct OpticalFrequencyJet {
    pub(in crate::plane_wave) epsilon: C64,
    pub(in crate::plane_wave) epsilon_derivatives: [C64; 2],
    pub(in crate::plane_wave) first_s_m: C64,
    pub(in crate::plane_wave) second_scaled_s_m: C64,
}

fn invalid(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}

fn interior(omega: f64, constant: f64, bounds: [f64; 2]) -> Result<f64> {
    let wavelength = dispersion::wavelength_from_frequency(omega, constant, bounds)
        .ok_or_else(|| invalid("optical dispersion is outside its original source band"))?;
    if omega <= constant / bounds[1] || omega >= constant / bounds[0] {
        return Err(invalid("optical frequency derivatives require a two-sided source neighbourhood"));
    }
    Ok(wavelength)
}

struct Line<const N: usize> {
    slope: [f64; N],
    intercept: [f64; N],
}

// Search the frequency images of ORIGINAL knots, before rounding lambda.
// A canonical knot frequency must not accidentally select one cusp side.
fn line<const N: usize>(
    omega: f64,
    constant: f64,
    count: usize,
    sample: impl Fn(usize) -> (f64, [f64; N]),
) -> Result<Line<N>> {
    if count < 2 {
        return Err(invalid("a monochromatic optical sample has no frequency neighbourhood"));
    }
    interior(omega, constant, [sample(0).0, sample(count - 1).0])?;
    let mut low = 0;
    let mut high = count;
    while low < high {
        let middle = low + (high - low) / 2;
        if omega < constant / sample(middle).0 {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    let upper = low;
    let segment = |left: usize| {
        let (a, va) = sample(left);
        let (b, vb) = sample(left + 1);
        let a = Dd::from_f64(a);
        let b = Dd::from_f64(b);
        let span = b - a;
        let slope = core::array::from_fn(|j| {
            ((Dd::from_f64(vb[j]) - Dd::from_f64(va[j])) / span).to_f64()
        });
        let intercept = core::array::from_fn(|j| {
            ((b * Dd::from_f64(va[j]) - a * Dd::from_f64(vb[j])) / span).to_f64()
        });
        Line { slope, intercept }
    };
    let selected = segment(upper - 1);
    if omega == constant / sample(upper).0
        && selected.slope != segment(upper).slope
    {
        return Err(invalid("an optical interpolation cusp has no two-sided frequency derivative"));
    }
    if !selected.slope.into_iter().chain(selected.intercept).all(f64::is_finite) {
        return Err(invalid("optical interpolation derivatives are unrepresentable"));
    }
    Ok(selected)
}

fn finish(
    index: C64,
    scale: f64,
    wavelength_first: C64,
    wavelength_second: C64,
    group_index: C64,
) -> Result<OpticalFrequencyJet> {
    let epsilon = dispersion::permittivity(index, scale)
        .ok_or_else(|| invalid("original optical permittivity is unrepresentable"))?;
    // epsilon=N²/(mu*c0²), with scaling BEFORE complex products.
    let n = index.scale(scale);
    let first = wavelength_first.scale(scale);
    let second = wavelength_second.scale(scale);
    let epsilon_derivatives = [
        -(n * first).scale(2.0),
        (first * first + n * (second + first.scale(2.0))).scale(2.0),
    ];
    let first_s_m = group_index.scale(VACUUM_SPEED_M_S.recip());
    let second_scaled_s_m = wavelength_second.scale(VACUUM_SPEED_M_S.recip());
    if !epsilon_derivatives.into_iter().chain([first_s_m, second_scaled_s_m])
        .all(super::super::finite)
    {
        return Err(invalid("optical frequency derivatives are unrepresentable"));
    }
    Ok(OpticalFrequencyJet {
        epsilon, epsilon_derivatives, first_s_m, second_scaled_s_m,
    })
}

// The tuple is (lambda*n_lambda, lambda²*n_lambda_lambda, group index).
fn cauchy(coefficients: [f64; 3], wavelength_m: f64) -> (f64, f64, f64) {
    let [a, b, c] = coefficients;
    let square = wavelength_m * wavelength_m;
    let b = b / square;
    let c = (c / square) / square;
    (-2.0 * b - 4.0 * c, 6.0 * b + 20.0 * c, a + 3.0 * b + 5.0 * c)
}

fn sellmeier(
    n: f64,
    strengths: [f64; 3],
    resonances_m2: [f64; 3],
    wavelength_m: f64,
) -> (f64, f64, f64) {
    let square = wavelength_m * wavelength_m;
    let mut first = Dd::from_f64(0.0);
    let mut second = Dd::from_f64(0.0);
    for (b, c) in strengths.into_iter().zip(resonances_m2) {
        if b == 0.0 || c == 0.0 {
            continue;
        }
        let denominator = square - c;
        let term = b * (square / denominator);
        let ratio = c / denominator;
        let moment = term * ratio;
        first = first + Dd::from_f64(-2.0 * moment);
        second = second + Dd::from_f64(moment * (6.0 + 8.0 * ratio));
    }
    let first = first.to_f64() / (2.0 * n);
    let second = second.to_f64() / (2.0 * n) - first * (first / n);
    (first, second, n - first)
}

// Alpha is intensity absorption. Its linear wavelength interpolation gives
// N_im=-alpha*lambda/(4*pi), so absorption contributes a quadratic index term.
// Form the group loss directly: subtracting N_im-lambda*N_im' loses weak slopes.
fn dielectric(
    index: C64,
    scale: f64,
    phase: (f64, f64, f64),
    wavelength_m: f64,
    scaled_alpha_slope: f64,
) -> Result<OpticalFrequencyJet> {
    let group_loss = scaled_alpha_slope * (wavelength_m / (4.0 * core::f64::consts::PI));
    finish(
        index,
        scale,
        C64::new(phase.0, index.im - group_loss),
        C64::new(phase.1, -2.0 * group_loss),
        C64::new(phase.2, group_loss),
    )
}

fn visible_alpha(
    absorption: VisibleDielectricAbsorption,
    omega: f64,
    wavelength_nm: f64,
) -> Result<f64> {
    let VisibleDielectricAbsorption::Spectral { absorption_per_m } = absorption else {
        return Err(invalid("RGB absorption has no monochromatic frequency derivative"));
    };
    let values = line(omega, FREQUENCY_WAVELENGTH_NM, 9, |i| {
        (VISIBLE_DIELECTRIC_WAVELENGTHS_NM[i], [absorption_per_m[i]])
    })?;
    Ok(wavelength_nm * values.slope[0])
}

impl OpticalConductorResponse<'_> {
    pub(in crate::plane_wave) fn frequency_jet(self, omega: f64) -> Result<OpticalFrequencyJet> {
        let samples = self.source.samples();
        let bounds = [samples[0].wavelength_nm, samples[samples.len() - 1].wavelength_nm];
        let wavelength = interior(omega, FREQUENCY_WAVELENGTH_NM, bounds)?;
        let values = line(omega, FREQUENCY_WAVELENGTH_NM, samples.len(), |i| {
            let sample = samples[i];
            (sample.wavelength_nm, [sample.eta, sample.k])
        })?;
        finish(
            self.complex_index(omega)?,
            self.index_to_sqrt_permittivity,
            C64::new(wavelength * values.slope[0], -wavelength * values.slope[1]),
            C64::ZERO,
            C64::new(values.intercept[0], -values.intercept[1]),
        )
    }
}

impl OpticalMaterialResponse<'_> {
    pub(in crate::plane_wave) fn frequency_jet(self, omega: f64) -> Result<OpticalFrequencyJet> {
        let wavelength = interior(omega, FREQUENCY_WAVELENGTH_NM, dispersion::VISIBLE_BAND_NM)?;
        if let VisibleOpticalStatePoint::Conductor(source) = self.source {
            return OpticalConductorResponse {
                source,
                permeability_h_m: self.permeability_h_m,
                index_to_sqrt_permittivity: self.index_to_sqrt_permittivity,
            }.frequency_jet(omega);
        }
        let index = self.complex_index(omega)?;
        let wavelength_m = wavelength * 1e-9;
        let (phase, alpha_slope) = match self.source {
            VisibleOpticalStatePoint::Dielectric(source) => (
                cauchy(source.cauchy_coefficients_si(), wavelength_m),
                visible_alpha(source.absorption(), omega, wavelength)?,
            ),
            VisibleOpticalStatePoint::TabulatedDielectric(source) => {
                let samples = source.samples();
                let values = line(omega, FREQUENCY_WAVELENGTH_NM, samples.len(), |i| {
                    (VISIBLE_DIELECTRIC_WAVELENGTHS_NM[i],
                     [samples[i].phase_index, samples[i].absorption_per_m])
                })?;
                ((wavelength * values.slope[0], 0.0, values.intercept[0]),
                 wavelength * values.slope[1])
            }
            VisibleOpticalStatePoint::DielectricSellmeier(source) => {
                let (strengths, resonances) = source.sellmeier_coefficients_si();
                (sellmeier(index.re, strengths, resonances, wavelength_m),
                 visible_alpha(source.absorption(), omega, wavelength)?)
            }
            VisibleOpticalStatePoint::Conductor(_) => unreachable!("conductor returned above"),
        };
        dielectric(index, self.index_to_sqrt_permittivity, phase, wavelength_m, alpha_slope)
    }
}

pub(super) fn spectral(
    response: SpectralOpticalResponse<'_>,
    omega: f64,
    scale: f64,
) -> Result<OpticalFrequencyJet> {
    let source = response.source();
    let complex = match source {
        SpectralOpticalSource::ComplexIndex(value) => Some(value),
        SpectralOpticalSource::Conductor(value) => Some(value.complex_index()),
        _ => None,
    };
    if let Some(source) = complex {
        let wavelength = interior(omega, FREQUENCY_WAVELENGTH_NM, source.wavelength_bounds_nm())?;
        let samples = source.samples();
        let values = line(omega, FREQUENCY_WAVELENGTH_NM, samples.len(), |i| {
            (samples[i].wavelength_nm, [samples[i].eta, samples[i].k])
        })?;
        return finish(
            response.complex_index(omega)?,
            scale,
            C64::new(wavelength * values.slope[0], -wavelength * values.slope[1]),
            C64::ZERO,
            C64::new(values.intercept[0], -values.intercept[1]),
        );
    }
    let constant = core::f64::consts::TAU * VACUUM_SPEED_M_S;
    let index = response.complex_index(omega)?;
    match source {
        SpectralOpticalSource::TabulatedDielectric(source) => {
            let (low, high) = source.wavelength_bounds_m();
            let wavelength = interior(omega, constant, [low, high])?;
            let samples = source.samples();
            let values = line(omega, constant, samples.len(), |i| {
                (samples[i].wavelength_m, [samples[i].phase_index, samples[i].absorption_per_m])
            })?;
            dielectric(
                index, scale,
                (wavelength * values.slope[0], 0.0, values.intercept[0]),
                wavelength, wavelength * values.slope[1],
            )
        }
        SpectralOpticalSource::DielectricSellmeier(source) => {
            let (low, high) = source.wavelength_bounds_m();
            let wavelength = interior(omega, constant, [low, high])?;
            let samples = source.absorption_samples();
            let values = line(omega, constant, samples.len(), |i| {
                (samples[i].wavelength_m, [samples[i].absorption_per_m])
            })?;
            let (strengths, resonances) = source.sellmeier_coefficients_si();
            dielectric(
                index, scale, sellmeier(index.re, strengths, resonances, wavelength),
                wavelength, wavelength * values.slope[0],
            )
        }
        _ => unreachable!("complex-index sources returned above"),
    }
}

pub(super) fn water(
    response: WaterOpticalResponse<'_>,
    omega: f64,
    scale: f64,
) -> Result<OpticalFrequencyJet> {
    let wavelength = interior(
        omega, FREQUENCY_WAVELENGTH_NM, fs_material::water_refraction::WAVELENGTH_RANGE_NM,
    )?;
    let sample = response.refractive_sample(omega)?;
    finish(
        C64::from_re(sample.phase_index),
        scale,
        C64::from_re(wavelength * sample.wavelength_derivative_per_nm),
        C64::from_re((wavelength * wavelength) * sample.wavelength_second_derivative_per_nm2),
        C64::from_re(sample.group_index),
    )
}
