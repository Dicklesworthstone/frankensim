//! Continuous-frequency finite-path amplification of the native Benney film.
//!
//! A steady path preserves inlet frequency. With local k=omega/c and c=3*q/h,
//! every completed path prefix has log gain A*omega²-B*omega⁴, B>=0. Its exact
//! model maximum is at a band endpoint or its one positive stationary point.
//! Within each constant section, gain is linear in distance, so section ends
//! also suffice for the maximum over the complete piecewise path.

use fs_lbm::film::{FilmError, LongWaveLimits};
use fs_math::det;

use super::{LaminarFilmError, LaminarPourConditions, LaminarTransitSection, section_film};

/// Finite-path transport over every inlet frequency in a closed SI band.
///
/// This uses the same weak-growth, locally parallel, slowly varying model as
/// `LaminarFilmTransitStudy`. It maximizes actual accumulated amplification,
/// rather than summing independently maximized local growth rates. Each trial
/// disturbance keeps its frequency as the local wavelength changes.
///
/// Every section must admit the WHOLE band under the native kh/Re*kh limits
/// and the declared |temporal growth / frequency| budget. Zero frequency uses
/// the continuous zero-growth limit. These are model restrictions, not error
/// enclosures, nonlinear saturation or a turbulence certificate.
#[derive(Debug, Clone, Copy)]
pub struct LaminarFilmTransitBandStudy<'a> {
    /// Positive-length sections in path order, approximating gradual variation.
    pub sections: &'a [LaminarTransitSection],
    /// Inclusive, finite, nonnegative [minimum, maximum] frequencies in rad/s.
    /// A point band is allowed; [0,0] is the exactly neutral conserved mode.
    pub angular_frequencies_rad_s: [f64; 2],
    /// Native long-wave budgets applied to every local wavelength in the band.
    pub limits: LongWaveLimits,
    /// Maximum |s/omega| anywhere in the band, strictly between zero and one.
    pub max_growth_frequency_ratio: f64,
}

/// Independent worst outlet and upstream amplification over the complete band.
/// Logarithmic amplitudes are never exponentiated or clipped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaminarTransitBandReport {
    /// Maximum natural log of outlet relative amplitude / inlet amplitude.
    /// Can be negative when the entire band decays at the outlet.
    pub max_outlet_log_amplification: f64,
    /// Inlet frequency attaining the outlet maximum; the lower wins exact ties.
    pub outlet_angular_frequency_rad_s: f64,
    /// Largest log gain anywhere along the path, including the inlet's zero.
    pub peak_log_amplification: f64,
    /// Completed sections at the first peak; zero denotes the inlet.
    pub peak_after_sections: usize,
    /// Inlet frequency attaining the upstream peak. None when the inlet wins:
    /// every inlet frequency has zero gain there, with no distinguished mode.
    pub peak_angular_frequency_rad_s: Option<f64>,
    /// Kinematic-wave travel time along the path, in seconds.
    pub wave_transit_time_s: f64,
    /// Mean liquid residence time, three times the wave travel time, in seconds.
    pub liquid_residence_time_s: f64,
    /// Smallest/largest local wavenumbers over all sections and frequencies.
    pub wavenumber_range_rad_m: [f64; 2],
    /// Minimum/maximum mean-speed film Reynolds numbers along the path.
    pub reynolds_range: [f64; 2],
    /// Largest local kh over the complete band and path.
    pub max_kh: f64,
    /// Largest local Re*kh over the complete band and path.
    pub max_re_kh: f64,
    /// Largest |s/omega|, including interior cubic extrema in the band.
    pub max_growth_frequency_ratio: f64,
}

impl LaminarFilmTransitBandStudy<'_> {
    /// Evaluate a width candidate in O(sections) work and constant storage.
    /// There is no frequency grid, search tolerance or eigensolver per trial.
    pub fn evaluate(
        &self,
        width_scale: f64,
        conditions: LaminarPourConditions,
    ) -> Result<LaminarTransitBandReport, LaminarFilmError> {
        if self.sections.is_empty() {
            return Err(LaminarFilmError::EmptyPath);
        }
        if !width_scale.is_finite() || width_scale <= 0.0 {
            return Err(LaminarFilmError::InvalidWidthScale);
        }
        if !conditions.volume_flow_m3_s.is_finite() || conditions.volume_flow_m3_s <= 0.0 {
            return Err(LaminarFilmError::InvalidVolumeFlow);
        }
        let [lower, upper] = self.angular_frequencies_rad_s;
        let budget = self.max_growth_frequency_ratio;
        if !lower.is_finite()
            || !upper.is_finite()
            || lower < 0.0
            || upper < lower
            || !budget.is_finite()
            || budget <= 0.0
            || budget >= 1.0
        {
            return Err(LaminarFilmError::InvalidTransitInput);
        }
        // Normalize by the upper frequency before forming the polynomial.
        // Powers of a dimensional frequency never enter the accumulated gain.
        let lower_fraction = if upper == 0.0 { 0.0 } else { lower / upper };
        if lower > 0.0 && lower_fraction == 0.0 {
            return Err(LaminarFilmError::InvalidTransitInput);
        }
        let frequency = |fraction: f64| {
            if fraction == lower_fraction {
                lower
            } else if fraction == 1.0 {
                upper
            } else {
                (fraction * upper).clamp(lower, upper)
            }
        };
        let mut quadratic = Sum::default();
        let mut quartic = Sum::default();
        let mut report = LaminarTransitBandReport {
            max_outlet_log_amplification: 0.0,
            outlet_angular_frequency_rad_s: lower,
            peak_log_amplification: 0.0,
            peak_after_sections: 0,
            peak_angular_frequency_rad_s: None,
            wave_transit_time_s: 0.0,
            liquid_residence_time_s: 0.0,
            wavenumber_range_rad_m: [f64::INFINITY, 0.0],
            reynolds_range: [f64::INFINITY, 0.0],
            max_kh: 0.0,
            max_re_kh: 0.0,
            max_growth_frequency_ratio: 0.0,
        };
        for (section, segment) in self.sections.iter().enumerate() {
            if !segment.length_m.is_finite() || segment.length_m <= 0.0 {
                return Err(LaminarFilmError::InvalidLength { section });
            }
            let film = section_film(segment.film, width_scale, conditions, section)?;
            let failure = |source| LaminarFilmError::Film { section, source };
            let speed = represented(3.0 * film.mean_velocity_m_s(), true).map_err(failure)?;
            let k_low = represented(lower / speed, lower > 0.0).map_err(failure)?;
            let k_high = represented(upper / speed, upper > 0.0).map_err(failure)?;
            let liquid_time =
                represented(segment.length_m / film.mean_velocity_m_s(), true).map_err(failure)?;
            let wave_time = represented(liquid_time / 3.0, true).map_err(failure)?;
            let band = film
                .long_wave_band([k_low, k_high], self.limits)
                .map_err(|source| LaminarFilmError::Band { section, source })?;
            let (quadratic_rate, quartic_rate, ratio) = if upper == 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                let [a, b] = film.long_wave_coefficients().map_err(failure)?;
                let a_high = represented((a * k_high) * k_high, a != 0.0).map_err(failure)?;
                let b_high = represented((((b * k_high) * k_high) * k_high) * k_high, b != 0.0)
                    .map_err(failure)?;
                let ratio = maximum_cubic_absolute(a_high / upper, b_high / upper, lower_fraction)
                    .map_err(failure)?;
                (a_high, b_high, ratio)
            };
            if ratio > budget {
                return Err(LaminarFilmError::OutsideSpatialGrowthRegime { section });
            }
            let a = quadratic
                .add(
                    represented(quadratic_rate * wave_time, quadratic_rate != 0.0)
                        .map_err(failure)?,
                )
                .map_err(failure)?;
            let b = quartic
                .add(represented(quartic_rate * wave_time, quartic_rate != 0.0).map_err(failure)?)
                .map_err(failure)?;
            let (gain, fraction) = maximum_quartic(a, b, lower_fraction).map_err(failure)?;
            let omega = frequency(fraction);
            report.max_outlet_log_amplification = gain;
            report.outlet_angular_frequency_rad_s = omega;
            if gain > report.peak_log_amplification {
                report.peak_log_amplification = gain;
                report.peak_after_sections = section + 1;
                report.peak_angular_frequency_rad_s = Some(omega);
            }
            report.wave_transit_time_s =
                represented(report.wave_transit_time_s + wave_time, true).map_err(failure)?;
            report.liquid_residence_time_s =
                represented(report.liquid_residence_time_s + liquid_time, true).map_err(failure)?;
            report.wavenumber_range_rad_m[0] = report.wavenumber_range_rad_m[0].min(k_low);
            report.wavenumber_range_rad_m[1] = report.wavenumber_range_rad_m[1].max(k_high);
            report.reynolds_range[0] = report.reynolds_range[0].min(film.reynolds());
            report.reynolds_range[1] = report.reynolds_range[1].max(film.reynolds());
            report.max_kh = report.max_kh.max(band.max_kh);
            report.max_re_kh = report.max_re_kh.max(band.max_re_kh);
            report.max_growth_frequency_ratio = report.max_growth_frequency_ratio.max(ratio);
        }
        Ok(report)
    }
}

fn represented(value: f64, nonzero: bool) -> Result<f64, FilmError> {
    if !value.is_finite() || (nonzero && value == 0.0) {
        Err(FilmError::NonRepresentable)
    } else {
        Ok(value)
    }
}

// Signed quadratic contributions can cancel between growing and damping
// sections. Keep that cancellation accurate without retaining the path.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) -> Result<f64, FilmError> {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
        represented(next + self.correction, false)
    }
}

fn quartic_value(a: f64, b: f64, x: f64) -> Result<f64, FilmError> {
    let driving = (-(b * x)).mul_add(x, a);
    represented((driving * x) * x, driving != 0.0 && x != 0.0)
}

// x is inlet frequency / upper frequency, so 0 <= lower <= x <= 1.
// Endpoints and the positive stationary point exhaust the quartic's maxima.
fn maximum_quartic(a: f64, b: f64, lower: f64) -> Result<(f64, f64), FilmError> {
    let mut maximum = (quartic_value(a, b, lower)?, lower);
    if a > 0.0 && b > 0.0 && 0.5 * a < b {
        // Separate square roots retain an interior root even if a/b underflows.
        let stationary = (det::sqrt(a) / det::sqrt(b)) * core::f64::consts::FRAC_1_SQRT_2;
        if stationary > lower && stationary < 1.0 {
            let gain = quartic_value(a, b, stationary)?;
            if gain > maximum.0 {
                maximum = (gain, stationary);
            }
        }
    }
    let gain = quartic_value(a, b, 1.0)?;
    if gain > maximum.0 {
        maximum = (gain, 1.0);
    }
    Ok(maximum)
}

// The weak-growth condition is |alpha*x-beta*x³|, not the quartic gain.
// Its interior maximum can violate the budget even when both endpoints pass.
fn maximum_cubic_absolute(linear: f64, cubic: f64, lower: f64) -> Result<f64, FilmError> {
    represented(linear, false)?;
    represented(cubic, false)?;
    let value = |x: f64| represented(((-(cubic * x)).mul_add(x, linear) * x).abs(), false);
    let mut maximum = value(lower)?.max(value(1.0)?);
    if linear > 0.0 && cubic > 0.0 && linear / 3.0 < cubic {
        let stationary = (det::sqrt(linear) / det::sqrt(cubic)) / det::sqrt(3.0);
        if stationary > lower && stationary < 1.0 {
            maximum = maximum.max(value(stationary)?);
        }
    }
    Ok(maximum)
}
