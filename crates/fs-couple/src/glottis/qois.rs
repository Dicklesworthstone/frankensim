//! Descriptive flow-record measurements, not clinical voice metrics.
//!
//! Pitch uses the fixed-window difference and cumulative mean normalization
//! of de Cheveigne and Kawahara, YIN (JASA 2002, DOI 10.1121/1.1458024),
//! equations 6–8 and section II.E. Unlike unconditional pitch tracking, this
//! observer reports no identified pitch when no trough passes admission.
//! A bounded local frequency refinement then minimizes the joint DC/harmonic
//! residual, so coarse pitch error does not masquerade as weak partials or noise.

use fs_fft::harmonic::{
    HarmonicFitBudget, HarmonicFitError, HarmonicFrequencyBudget, HarmonicFrequencyError,
    refine_harmonic_frequency,
};
use fs_fft::{C64, Fft};
use fs_math::det;

/// Explicit limits for one uniformly sampled glottal-flow record.
#[derive(Debug, Clone, Copy)]
pub struct GlottalQoiOptions {
    /// Lowest admitted fundamental frequency [Hz].
    pub min_f0_hz: f64,
    /// Highest admitted fundamental [Hz], strictly below Nyquist.
    pub max_f0_hz: f64,
    /// Largest admitted cumulative-mean normalized difference, in (0, 1).
    /// This is an estimator threshold, not a calibrated voicing probability.
    pub periodicity_threshold: f64,
    /// Absolute local frequency-correction tolerance [Hz]. This controls
    /// numerical convergence, not physical pitch uncertainty.
    pub frequency_tolerance_hz: f64,
    /// Maximum joint fits during local frequency refinement of a voiced record.
    /// Exhaustion is an explicit error, not acceptance of an unfinished fit.
    pub max_frequency_evaluations: usize,
    /// Conservative limit on temporary numeric storage, excluding inputs
    /// and allocator metadata. Both FFT and joint-fit storage are admitted.
    pub max_workspace_bytes: usize,
}

impl Default for GlottalQoiOptions {
    fn default() -> Self {
        Self {
            min_f0_hz: 40.0,
            max_f0_hz: 500.0,
            periodicity_threshold: 0.1,
            frequency_tolerance_hz: 1e-6,
            max_frequency_evaluations: 24,
            max_workspace_bytes: 256 * 1024 * 1024,
        }
    }
}

/// A refused record produces no partially valid measurement.
#[derive(Debug, Clone, PartialEq)]
pub enum GlottalQoiError {
    /// Invalid samples, shape, clock, analysis range, or overflowing size.
    Invalid {
        /// The failed admission condition.
        what: &'static str,
    },
    /// Numeric storage would exceed the caller's limit.
    Workspace {
        /// Conservative required bytes.
        required_bytes: usize,
        /// Supplied limit.
        max_bytes: usize,
    },
    /// The harmonic fit could not separate the declared columns reliably.
    HarmonicFit(HarmonicFitError),
    /// The local frequency estimate was unidentifiable or failed to converge
    /// within its declared range, evaluation limit, or numeric budget.
    FrequencyRefinement(HarmonicFrequencyError),
}

impl core::fmt::Display for GlottalQoiError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "glottal QoI: {what}"),
            Self::Workspace {
                required_bytes,
                max_bytes,
            } => {
                write!(
                    f,
                    "glottal QoI requires {required_bytes} bytes, budget is {max_bytes}"
                )
            }
            Self::HarmonicFit(error) => error.fmt(f),
            Self::FrequencyRefinement(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for GlottalQoiError {}

/// Descriptive QoIs of the supplied record, without a voice-standard claim.
#[derive(Debug, Clone, Copy)]
pub struct GlottalQois {
    /// Fraction of samples with positive gap; a never-closing cycle has 1.0.
    pub open_quotient: f64,
    /// Regression slope [dB/octave] of resolved members of harmonics 1..6.
    /// NaN if fewer than two levels exceed the numerical and residual floor.
    pub spectral_slope_db_oct: f64,
    /// Fractional mean-crossing period spread (std/mean). NaN when fewer than
    /// three periods or ambiguous/multiple crossings prevent cycle labeling.
    pub jitter: f64,
    /// Identified fundamental [Hz]; zero for constant or unvoiced records.
    pub f0_hz: f64,
    /// Jointly fitted peak amplitudes in input flow units, orders 1..6.
    /// NaN above Nyquist or without an identified periodic model.
    pub harmonic_peak_amplitudes: [f64; 6],
    /// RMS unexplained by the fitted DC and admitted harmonics, in input flow
    /// units. For unvoiced records this is the DC-only residual. It is not
    /// an uncertainty bound or a separation of noise from omitted physics.
    pub harmonic_residual_rms: f64,
}

impl GlottalQois {
    fn unvoiced(open: f64, residual: f64) -> Self {
        Self {
            open_quotient: open,
            spectral_slope_db_oct: f64::NAN,
            jitter: f64::NAN,
            f0_hz: 0.0,
            harmonic_peak_amplitudes: [f64::NAN; 6],
            harmonic_residual_rms: residual,
        }
    }
}

/// Measure a flow record using the default 40–500 Hz analysis band.
/// Invalid inputs return NaNs; use [`glottal_qois_with_options`] for errors.
#[must_use]
pub fn glottal_qois(flow: &[f64], gaps: &[f64], rate: f64) -> GlottalQois {
    glottal_qois_with_options(flow, gaps, rate, GlottalQoiOptions::default()).unwrap_or_else(|_| {
        let mut invalid = GlottalQois::unvoiced(f64::NAN, f64::NAN);
        invalid.f0_hz = f64::NAN;
        invalid
    })
}

/// Measure finite, equally long flow/gap records sampled at `rate` Hz.
/// At least two periods at an admitted frequency must fit the record.
/// Complexity is O(N log N) for coarse pitch and O(E N) for the fixed six-tone
/// refinement, where E is bounded by `options.max_frequency_evaluations`.
/// The entire record is assumed approximately stationary; chirps, aspiration,
/// subharmonics and transients can remain in the exposed fit residual.
pub fn glottal_qois_with_options(
    flow: &[f64],
    gaps: &[f64],
    rate: f64,
    options: GlottalQoiOptions,
) -> Result<GlottalQois, GlottalQoiError> {
    let invalid = |what| GlottalQoiError::Invalid { what };
    let n = flow.len();
    if n < 8 || n != gaps.len() || flow.iter().chain(gaps).any(|x| !x.is_finite()) {
        return Err(invalid(
            "need equal finite records with at least eight samples",
        ));
    }
    if !(rate.is_finite()
        && rate > 0.0
        && options.min_f0_hz.is_finite()
        && options.min_f0_hz > 0.0
        && options.max_f0_hz.is_finite()
        && options.max_f0_hz > options.min_f0_hz
        && options.max_f0_hz < 0.5 * rate
        && options.periodicity_threshold.is_finite()
        && options.periodicity_threshold > 0.0
        && options.periodicity_threshold < 1.0
        && options.frequency_tolerance_hz.is_finite()
        && options.frequency_tolerance_hz > 0.0
        && options.max_frequency_evaluations > 0)
    {
        return Err(invalid(
            "invalid clock, sub-Nyquist band, periodicity threshold, or refinement limits",
        ));
    }
    let min_period = rate / options.max_f0_hz;
    let max_period = (rate / options.min_f0_hz).min(n as f64 / 2.0);
    if min_period > max_period {
        return Err(invalid("record is too short for the declared pitch band"));
    }
    // Include integer centers bracketing the continuous band edges; ceiling
    // the lower edge can skip the actual trough and select its double period.
    let lo = (min_period.floor() as usize).max(2);
    let hi = max_period.ceil() as usize;
    let open = gaps.iter().filter(|&&g| g > 0.0).count() as f64 / n as f64;
    if flow.iter().all(|&x| x == flow[0]) {
        return Ok(GlottalQois::unvoiced(open, 0.0));
    }
    let fft_len = n
        .checked_next_power_of_two()
        .ok_or_else(|| invalid("FFT size overflow"))?;
    // 64L covers three complex FFT arrays, twiddles and subplan headroom.
    // 256N covers the six-tone fit and real/prefix/lag arrays and crossing
    // marks. The refiner separately admits its entire numeric peak after
    // FFT scratch and the AC record are dropped.
    let required_bytes = n
        .checked_mul(256)
        .and_then(|bytes| {
            fft_len
                .checked_mul(64)
                .and_then(|extra| bytes.checked_add(extra))
        })
        .and_then(|bytes| bytes.checked_add(4096))
        .ok_or_else(|| invalid("workspace size overflow"))?;
    if required_bytes > options.max_workspace_bytes {
        return Err(GlottalQoiError::Workspace {
            required_bytes,
            max_bytes: options.max_workspace_bytes,
        });
    }
    // Subtract the preload before scaling whenever representable; scaling the
    // full DC first can erase small but representable physical variations.
    let mut ac: Vec<f64> = flow.iter().map(|x| x - flow[0]).collect();
    let mut scale = ac.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    if !scale.is_finite() {
        scale = flow.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        for (x, &value) in ac.iter_mut().zip(flow) {
            *x = value / scale;
        }
    } else {
        for x in &mut ac {
            *x /= scale;
        }
    }
    let mean = ac.iter().sum::<f64>() / n as f64;
    for x in &mut ac {
        *x -= mean;
    }
    let rms = det::sqrt(ac.iter().map(|x| x * x).sum::<f64>() / n as f64) * scale;
    if !rms.is_finite() {
        return Err(invalid("flow residual RMS is not representable"));
    }
    let mut result = GlottalQois::unvoiced(open, rms);
    let Some((period, jitter)) = yin_period(
        &ac,
        fft_len,
        lo,
        hi,
        min_period,
        max_period,
        options.periodicity_threshold,
    ) else {
        return Ok(result);
    };
    drop(ac);
    result.jitter = jitter;
    let initial_f0 = rate / period;
    // This is a local correction of an admitted periodic model, not an octave
    // search. One inverse record duration or 2% pitch bounds the neighborhood.
    let half_width = (0.02 * initial_f0).min(rate / (n - 1) as f64);
    let lower = (initial_f0 - half_width).max(options.min_f0_hz);
    let upper = (initial_f0 + half_width).min(options.max_f0_hz);
    // Keep one fixed set of columns throughout refinement. If the refined
    // pitch admits one more order, complete refinement with that order inside
    // a tightened Nyquist bound and the remaining evaluation budget.
    let all_orders = [1, 2, 3, 4, 5, 6];
    let order_count = all_orders
        .iter()
        .take_while(|&&h| upper * (h as f64) < 0.5 * rate)
        .count();
    let orders = &all_orders[..order_count];
    let fit_budget = HarmonicFitBudget {
        max_workspace_bytes: options.max_workspace_bytes - 6 * core::mem::size_of::<usize>(),
        min_reciprocal_condition: 1e-10,
    };
    let refined = refine_harmonic_frequency(
        flow,
        1.0 / rate,
        initial_f0,
        orders,
        HarmonicFrequencyBudget {
            frequency_bounds_hz: [lower, upper],
            absolute_tolerance_hz: options.frequency_tolerance_hz,
            max_evaluations: options.max_frequency_evaluations,
            fit: fit_budget,
        },
    )
    .map_err(GlottalQoiError::FrequencyRefinement)?;
    result.f0_hz = refined.fundamental_hz;
    let final_count = all_orders
        .iter()
        .take_while(|&&h| result.f0_hz * (h as f64) < 0.5 * rate)
        .count();
    let final_orders = &all_orders[..final_count];
    let fit = if final_orders.len() == orders.len() {
        refined.fit
    } else {
        let used = refined.evaluations;
        let remaining = options.max_frequency_evaluations.saturating_sub(used);
        if remaining == 0 {
            return Err(GlottalQoiError::FrequencyRefinement(
                HarmonicFrequencyError::NotConverged {
                    fundamental_hz: result.f0_hz,
                    residual_rms: refined.fit.residual_rms,
                    evaluations: used,
                    last_correction_hz: refined.last_correction_hz,
                },
            ));
        }
        drop(refined.fit);
        // At most one order can cross Nyquist in a +/-2% interval for orders
        // 1..6. The second solve uses every order admissible at its seed.
        let upper = upper.min((0.5 * rate / final_orders.len() as f64).next_down());
        let complete = refine_harmonic_frequency(
            flow,
            1.0 / rate,
            result.f0_hz,
            final_orders,
            HarmonicFrequencyBudget {
                frequency_bounds_hz: [lower, upper],
                absolute_tolerance_hz: options.frequency_tolerance_hz,
                max_evaluations: remaining,
                fit: fit_budget,
            },
        )
        .map_err(|error| {
            GlottalQoiError::FrequencyRefinement(match error {
                HarmonicFrequencyError::NotConverged {
                    fundamental_hz,
                    residual_rms,
                    evaluations,
                    last_correction_hz,
                } => HarmonicFrequencyError::NotConverged {
                    fundamental_hz,
                    residual_rms,
                    evaluations: used + evaluations,
                    last_correction_hz,
                },
                other => other,
            })
        })?;
        result.f0_hz = complete.fundamental_hz;
        complete.fit
    };
    for (amplitude, phasor) in result.harmonic_peak_amplitudes.iter_mut().zip(&fit.phasors) {
        *amplitude = phasor.re.hypot(phasor.im);
        if !amplitude.is_finite() {
            return Err(invalid("harmonic peak amplitude is not representable"));
        }
    }
    result.harmonic_residual_rms = fit.residual_rms;
    result.spectral_slope_db_oct =
        harmonic_slope(&result.harmonic_peak_amplitudes, fit.residual_rms);
    Ok(result)
}

// With W=N-max_lag, j+lag never wraps for the requested positive lags, so
// inverse(conj(FFT(x[0..W]))*FFT(x)) needs padding only to next_power_of_two(N).
fn fixed_window_difference(x: &[f64], fft_len: usize, max_lag: usize) -> Vec<f64> {
    let width = x.len() - max_lag;
    let fft = Fft::new(fft_len);
    let mut a = vec![C64::new(0.0, 0.0); fft_len];
    let mut b = a.clone();
    let mut scratch = a.clone();
    for (value, &sample) in a.iter_mut().zip(&x[..width]) {
        value.re = sample;
    }
    for (value, &sample) in b.iter_mut().zip(x) {
        value.re = sample;
    }
    fft.forward(&mut a, &mut scratch);
    fft.forward(&mut b, &mut scratch);
    for (a, b) in a.iter_mut().zip(&b) {
        *a = C64::new(
            a.re.mul_add(b.re, a.im * b.im),
            a.re.mul_add(b.im, -a.im * b.re),
        );
    }
    fft.inverse(&mut a, &mut scratch);
    let mut energy = vec![0.0; x.len() + 1];
    for (i, &sample) in x.iter().enumerate() {
        energy[i + 1] = sample.mul_add(sample, energy[i]);
    }
    (0..=max_lag)
        .map(|lag| (energy[width] + energy[width + lag] - energy[lag] - 2.0 * a[lag].re).max(0.0))
        .collect()
}

fn parabola(values: &[f64], i: usize) -> (f64, f64) {
    let curvature = values[i - 1] - 2.0 * values[i] + values[i + 1];
    let offset = if curvature > 0.0 {
        (0.5 * (values[i - 1] - values[i + 1]) / curvature).clamp(-0.5, 0.5)
    } else {
        0.0
    };
    (
        offset,
        values[i] - 0.25 * (values[i - 1] - values[i + 1]) * offset,
    )
}

fn yin_period(
    x: &[f64],
    fft_len: usize,
    lo: usize,
    hi: usize,
    min_period: f64,
    max_period: f64,
    threshold: f64,
) -> Option<(f64, f64)> {
    let difference = fixed_window_difference(x, fft_len, hi + 1);
    let mut normalized = vec![1.0; difference.len()];
    let mut sum = 0.0;
    for lag in 1..difference.len() {
        sum += difference[lag];
        if sum > 0.0 {
            normalized[lag] = difference[lag] * lag as f64 / sum;
        }
    }
    for lag in lo..=hi {
        if normalized[lag] <= normalized[lag - 1]
            && normalized[lag] <= normalized[lag + 1]
            && parabola(&normalized, lag).1 < threshold
        {
            // Select with d', but interpolate raw d to avoid normalization bias.
            let mut period = lag as f64 + parabola(&difference, lag).0;
            let (jitter, crossing_period) = crossing_periods(x, period);
            if jitter < 0.01 && (crossing_period / period - 1.0).abs() < 0.05 {
                period = crossing_period;
            }
            // A trough just outside the requested continuous band must not
            // become a clamped boundary value or its in-band double period.
            return (min_period..=max_period)
                .contains(&period)
                .then_some((period, jitter));
        }
    }
    None
}

fn crossing_periods(ac: &[f64], period: f64) -> (f64, f64) {
    let marks: Vec<f64> = ac
        .windows(2)
        .enumerate()
        .filter_map(|(i, w)| (w[0] <= 0.0 && w[1] > 0.0).then(|| i as f64 - w[0] / (w[1] - w[0])))
        .collect();
    if marks.len() < 4
        || marks
            .windows(2)
            .any(|w| !(0.6 * period..=1.4 * period).contains(&(w[1] - w[0])))
    {
        return (f64::NAN, f64::NAN);
    }
    let mean = (marks[marks.len() - 1] - marks[0]) / (marks.len() - 1) as f64;
    let variance = marks
        .windows(2)
        .map(|w| (w[1] - w[0] - mean).powi(2))
        .sum::<f64>()
        / (marks.len() - 1) as f64;
    (det::sqrt(variance) / mean, mean)
}

fn harmonic_slope(amplitudes: &[f64; 6], residual: f64) -> f64 {
    let peak = amplitudes
        .iter()
        .filter(|a| a.is_finite())
        .fold(0.0_f64, |a, &b| a.max(b));
    let floor = (64.0 * f64::EPSILON).max(core::f64::consts::SQRT_2 * (residual / peak));
    let levels: Vec<(f64, f64)> = amplitudes
        .iter()
        .enumerate()
        .filter_map(|(i, &a)| {
            let relative = a / peak;
            (relative > floor).then(|| {
                (
                    det::ln((i + 1) as f64) / core::f64::consts::LN_2,
                    20.0 * det::ln(relative) / core::f64::consts::LN_10,
                )
            })
        })
        .collect();
    if levels.len() < 2 {
        return f64::NAN;
    }
    let mx = levels.iter().map(|&(x, _)| x).sum::<f64>() / levels.len() as f64;
    let my = levels.iter().map(|&(_, y)| y).sum::<f64>() / levels.len() as f64;
    levels
        .iter()
        .map(|&(x, y)| (x - mx) * (y - my))
        .sum::<f64>()
        / levels.iter().map(|&(x, _)| (x - mx).powi(2)).sum::<f64>()
}

#[cfg(test)]
#[path = "qois/tests.rs"]
mod tests;
