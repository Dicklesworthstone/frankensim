//! Joint Welch cross spectra as positive-semidefinite factors.
//!
//! Every channel shares the same sample clock, segment boundaries, window and
//! detrending. A bin retains one complex channel vector per segment, rather
//! than discarding relative phase or building a dense channel covariance.
//! With the negative-time phasor convention, `S_ab = sum_r B_ar conj(B_br)`
//! is the one-sided cross-spectral density per Hz. Here
//! `B_ar = conj(FFT(w*x_a)_r) * sqrt(d*dt/(R*sum(w²)))`, where d=2 for
//! interior bins and d=1 at DC/Nyquist. These are NOT coherent-gain tone fits.
//!
//! This is Welch's averaged modified periodogram, extended by retaining its
//! outer-product factors: <https://doi.org/10.1109/TAU.1967.1161901>.
//! Overlapping segments need not be statistically independent for this
//! algebraic factorization. It supplies no degrees of freedom, confidence
//! interval, stationarity guarantee or cross-frequency covariance. A finite
//! record's leakage and resolution limits remain. In particular, one segment
//! has rank one and cannot establish physical coherence by itself.

use crate::{C64, Fft};
use fs_math::{dd::Dd, det};

#[path = "welch/stream.rs"]
mod stream;
pub use stream::WelchStream;

/// Bounds the uninterrupted work of one existing FFT. Cancellation is also
/// checked during input admission and between segment/channel transforms.
pub const MAX_WELCH_SEGMENT_SAMPLES: usize = 65_536;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WelchWindow {
    Rectangular,
    /// Periodic Hann: `w[j]=(1-cos(2*pi*j/N))/2`, j=0..N-1.
    PeriodicHann,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WelchDetrend {
    None,
    /// Subtract each channel's unweighted mean on each segment before tapering.
    SegmentMean,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WelchConfig {
    /// Power of two in 4..=MAX_WELCH_SEGMENT_SAMPLES. No zero padding.
    pub segment_samples: usize,
    /// Segment advance in 1..=segment_samples. Incomplete tails are reported.
    pub hop_samples: usize,
    pub window: WelchWindow,
    pub detrend: WelchDetrend,
}

/// Limits cover input scanning, FFT count and all retained/temporary buffers.
/// Caller-owned input storage and allocator metadata are excluded.
#[derive(Clone, Copy, Debug)]
pub struct WelchBudget {
    pub maximum_channels: usize,
    pub maximum_samples_per_channel: usize,
    pub maximum_transforms: usize,
    pub maximum_workspace_bytes: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WelchError {
    Invalid(&'static str),
    BudgetExceeded(&'static str),
    NonFinite,
    Unrepresentable,
    Allocation,
    Cancelled,
}

impl core::fmt::Display for WelchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(why) => write!(f, "invalid Welch observation: {why}"),
            Self::BudgetExceeded(what) => write!(f, "Welch budget exceeded: {what}"),
            Self::NonFinite => f.write_str("Welch input contains a nonfinite sample"),
            Self::Unrepresentable => f.write_str("Welch result is not representable"),
            Self::Allocation => f.write_str("Welch buffer allocation failed"),
            Self::Cancelled => f.write_str("Welch observation cancelled"),
        }
    }
}
impl std::error::Error for WelchError {}

/// Immutable joint spectrum, with factors ordered bin, segment, channel.
/// Storage is O(bins * segments * channels), never O(bins * channels²).
#[derive(Clone, Debug)]
pub struct WelchCrossSpectrum {
    config: WelchConfig,
    timestep_s: f64,
    bin_width_hz: f64,
    channels: usize,
    segments: usize,
    input_samples: usize,
    used_samples: usize,
    factors: Vec<C64>,
}

/// Borrowed interior bin admitted at a physical observer's exact frequency.
///
/// The factors remain in input-units / sqrt(Hz). Apply [`Self::peak_factor`]
/// once to each coefficient before a peak-phasor observer: its quadratic
/// quantities then represent this bin's power, including the observer's 1/2
/// time average. This is not a coherent-gain tone fit or a cross-frequency
/// covariance. The caller still binds channel order, units and source identity.
#[derive(Clone, Copy, Debug)]
pub struct WelchPeakBin<'a> {
    factors: &'a [C64],
    channels: usize,
    segments: usize,
    frequency_hz: f64,
    peak_scale: f64,
}

impl<'a> WelchPeakBin<'a> {
    pub const fn channel_count(self) -> usize {
        self.channels
    }
    pub const fn segment_count(self) -> usize {
        self.segments
    }
    pub const fn frequency_hz(self) -> f64 {
        self.frequency_hz
    }

    /// Borrow the original **unscaled** input-units / sqrt(Hz) coefficients.
    /// Each slice is one segment's factor, in the permanent channel order.
    /// No coefficient is copied, allocated, phase-rotated or normalized here.
    pub fn factors(self) -> impl ExactSizeIterator<Item = &'a [C64]> {
        self.factors.chunks_exact(self.channels)
    }

    /// Convert one raw coefficient to peak-equivalent units for this bin.
    /// Multiplication by sqrt(2*df) supplies bin width and the peak convention
    /// exactly once. A nonzero component that rounds to zero, overflow, and
    /// nonfinite input are refused; genuine zero components remain zero.
    pub fn peak_factor(self, coefficient: C64) -> Result<C64, WelchError> {
        let convert = |value: f64| {
            if !value.is_finite() {
                return Err(WelchError::NonFinite);
            }
            let result = Dd::from_f64(value) * Dd::from_f64(self.peak_scale);
            if !result.hi.is_finite()
                || !result.lo.is_finite()
                || (value != 0.0 && result.to_f64() == 0.0)
            {
                Err(WelchError::Unrepresentable)
            } else {
                Ok(result.to_f64())
            }
        };
        Ok(C64::new(convert(coefficient.re)?, convert(coefficient.im)?))
    }
}

impl WelchCrossSpectrum {
    /// Estimate synchronized real channels in their original physical units.
    /// Channels must be nonempty and have equal lengths. No interpolation,
    /// padding, extrapolation, phase synthesis or missing-data fill is applied.
    pub fn estimate(
        channels: &[&[f64]],
        timestep_s: f64,
        config: WelchConfig,
        budget: WelchBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        let samples = channels.first().map_or(0, |c| c.len());
        if channels.iter().any(|c| c.len() != samples) {
            return Err(WelchError::Invalid("channel sample counts differ"));
        }
        Self::estimate_with(
            channels.len(),
            samples,
            |c, j| channels[c][j],
            timestep_s,
            config,
            budget,
            cancelled,
        )
    }

    /// Observe an existing Cartesian record without transposing or copying it.
    pub fn estimate_interleaved<const D: usize>(
        samples: &[[f64; D]],
        timestep_s: f64,
        config: WelchConfig,
        budget: WelchBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        Self::estimate_with(
            D,
            samples.len(),
            |c, j| samples[j][c],
            timestep_s,
            config,
            budget,
            cancelled,
        )
    }

    /// Estimate a dynamically shaped record without transposing or copying it.
    /// `sample(channel, time_index)` must return the same physical sample on
    /// every call; it is read during admission and again during transforms.
    /// Channel order and the declared clock remain fixed for the whole record.
    /// The same original shape, work, allocation and cancellation bounds apply.
    #[allow(clippy::too_many_arguments)]
    pub fn estimate_with(
        channels: usize,
        samples: usize,
        sample: impl Fn(usize, usize) -> f64,
        timestep_s: f64,
        config: WelchConfig,
        budget: WelchBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        poll(&mut cancelled)?;
        let shape = admit(channels, samples, timestep_s, config, budget, 0)?;
        // Validate the entire supplied record, including its explicitly unused
        // tail, before any estimate or allocation can be published.
        validate_samples(channels, samples, &sample, &mut cancelled)?;
        let mut transform =
            SegmentTransform::new(timestep_s, config, shape.segments, &mut cancelled)?;
        let mut factors = allocate(shape.factor_count, C64::default())?;
        for segment in 0..shape.segments {
            let start = segment * config.hop_samples;
            for channel in 0..channels {
                transform.channel(
                    |j| sample(channel, start + j),
                    |bin, value| {
                        factors[(bin * shape.segments + segment) * channels + channel] = value
                    },
                    &mut cancelled,
                )?;
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            config,
            timestep_s,
            bin_width_hz: shape.bin_width_hz,
            channels,
            segments: shape.segments,
            input_samples: samples,
            used_samples: (shape.segments - 1) * config.hop_samples + config.segment_samples,
            factors,
        })
    }

    pub const fn config(&self) -> WelchConfig {
        self.config
    }
    pub const fn timestep_s(&self) -> f64 {
        self.timestep_s
    }
    pub const fn bin_width_hz(&self) -> f64 {
        self.bin_width_hz
    }
    pub const fn channel_count(&self) -> usize {
        self.channels
    }
    pub const fn segment_count(&self) -> usize {
        self.segments
    }
    pub const fn bin_count(&self) -> usize {
        self.config.segment_samples / 2 + 1
    }
    pub const fn used_sample_count(&self) -> usize {
        self.used_samples
    }
    pub const fn discarded_tail_samples(&self) -> usize {
        self.input_samples - self.used_samples
    }

    pub fn frequency_hz(&self, bin: usize) -> Result<f64, WelchError> {
        self.check_bin(bin)?;
        Ok(bin as f64 * self.bin_width_hz)
    }

    /// One algebraic factor per segment, each in permanent channel order.
    /// Units are input-units / sqrt(Hz), in the negative-time convention.
    pub fn factors(&self, bin: usize) -> Result<impl ExactSizeIterator<Item = &[C64]>, WelchError> {
        self.check_bin(bin)?;
        let count = self.channels * self.segments;
        Ok(self.factors[bin * count..(bin + 1) * count].chunks_exact(self.channels))
    }

    /// Admit an interior bin at an observer's original frequency in Hz.
    /// Exact equality is required; no interpolation or frequency snapping is
    /// inferred. DC and Nyquist do not represent the positive harmonic used by
    /// peak-phasor propagation and are refused before any factor is exposed.
    pub fn peak_bin_hz(
        &self,
        bin: usize,
        expected_frequency_hz: f64,
    ) -> Result<WelchPeakBin<'_>, WelchError> {
        let view = self.interior_peak_bin(bin)?;
        if view.frequency_hz.to_bits() != expected_frequency_hz.to_bits() {
            return Err(WelchError::Invalid(
                "observer frequency must exactly match the Welch bin in Hz",
            ));
        }
        Ok(view)
    }

    /// Admit an interior bin at an observer's original angular frequency.
    /// Match exactly `TAU * bin_frequency_hz`, without dividing the observer's
    /// angular frequency back to Hz: that round trip can move a valid bin by
    /// one ULP. An unrepresentable angular frequency is refused.
    pub fn peak_bin_rad_s(
        &self,
        bin: usize,
        expected_angular_frequency_rad_s: f64,
    ) -> Result<WelchPeakBin<'_>, WelchError> {
        let view = self.interior_peak_bin(bin)?;
        let omega = core::f64::consts::TAU * view.frequency_hz;
        if !omega.is_finite() || omega.to_bits() != expected_angular_frequency_rad_s.to_bits() {
            return Err(WelchError::Invalid(
                "observer angular frequency must exactly match TAU times the Welch bin",
            ));
        }
        Ok(view)
    }

    fn interior_peak_bin(&self, bin: usize) -> Result<WelchPeakBin<'_>, WelchError> {
        if bin == 0 || bin >= self.bin_count() - 1 {
            return Err(WelchError::Invalid(
                "peak-phasor observation requires an interior positive Welch bin",
            ));
        }
        let count = self.channels * self.segments;
        Ok(WelchPeakBin {
            factors: &self.factors[bin * count..(bin + 1) * count],
            channels: self.channels,
            segments: self.segments,
            frequency_hz: bin as f64 * self.bin_width_hz,
            peak_scale: core::f64::consts::SQRT_2 * det::sqrt(self.bin_width_hz),
        })
    }

    /// One-sided `S_ab` per Hz; the diagonal is channel power density.
    /// Does not invent phase/covariance between DIFFERENT frequency bins.
    pub fn cross_spectrum(&self, bin: usize, a: usize, b: usize) -> Result<C64, WelchError> {
        if a >= self.channels || b >= self.channels {
            return Err(WelchError::Invalid(
                "channel index outside recorded spectrum",
            ));
        }
        let mut re = Dd::ZERO;
        let mut im = Dd::ZERO;
        for factor in self.factors(bin)? {
            let ar = Dd::from_f64(factor[a].re);
            let ai = Dd::from_f64(factor[a].im);
            let br = Dd::from_f64(factor[b].re);
            let bi = Dd::from_f64(factor[b].im);
            re = re + ar * br + ai * bi;
            im = im + ai * br - ar * bi;
        }
        let result = C64::new(re.to_f64(), im.to_f64());
        if result.re.is_finite() && result.im.is_finite() {
            Ok(result)
        } else {
            Err(WelchError::Unrepresentable)
        }
    }

    fn check_bin(&self, bin: usize) -> Result<(), WelchError> {
        if bin < self.bin_count() {
            Ok(())
        } else {
            Err(WelchError::Invalid("bin outside recorded spectrum"))
        }
    }
}

struct WelchShape {
    bin_width_hz: f64,
    segments: usize,
    factor_count: usize,
}

fn admit(
    channels: usize,
    samples: usize,
    timestep_s: f64,
    config: WelchConfig,
    budget: WelchBudget,
    input_buffers: usize,
) -> Result<WelchShape, WelchError> {
    let n = config.segment_samples;
    if channels == 0
        || samples < n
        || n < 4
        || !n.is_power_of_two()
        || n > MAX_WELCH_SEGMENT_SAMPLES
        || config.hop_samples == 0
        || config.hop_samples > n
    {
        return Err(WelchError::Invalid(
            "nonempty channels and valid segment/hop required",
        ));
    }
    if !timestep_s.is_finite() || timestep_s <= 0.0 || !(1.0 / timestep_s).is_finite() {
        return Err(WelchError::Invalid("positive finite sample rate required"));
    }
    let bin_width_hz = (1.0 / timestep_s) / n as f64;
    if bin_width_hz == 0.0 {
        return Err(WelchError::Unrepresentable);
    }
    let segments = 1 + (samples - n) / config.hop_samples;
    let transforms = channels
        .checked_mul(segments)
        .ok_or(WelchError::BudgetExceeded("transform count overflow"))?;
    if channels > budget.maximum_channels
        || samples > budget.maximum_samples_per_channel
        || transforms > budget.maximum_transforms
    {
        return Err(WelchError::BudgetExceeded(
            "channels, samples or transforms",
        ));
    }
    let factor_count = (n / 2 + 1)
        .checked_mul(transforms)
        .ok_or(WelchError::BudgetExceeded("factor count overflow"))?;
    let mut workspace = factor_count
        .checked_mul(core::mem::size_of::<C64>())
        .and_then(|v| v.checked_add(n * (8 + 2 * core::mem::size_of::<C64>())))
        .ok_or(WelchError::BudgetExceeded("workspace size overflow"))?;
    if input_buffers != 0 {
        let input_bytes = n
            .checked_mul(channels)
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|v| v.checked_add(core::mem::size_of::<Vec<f64>>()))
            .and_then(|v| v.checked_mul(input_buffers))
            .ok_or(WelchError::BudgetExceeded("stream input storage overflow"))?;
        workspace = workspace
            .checked_add(input_bytes)
            .ok_or(WelchError::BudgetExceeded("workspace size overflow"))?;
    }
    // Include the original FFT plan and its optional six-step subplans.
    let mut plan_n = n;
    loop {
        workspace = workspace
            .checked_add(plan_n / 2 * core::mem::size_of::<C64>())
            .and_then(|v| v.checked_add(core::mem::size_of::<Fft>()))
            .ok_or(WelchError::BudgetExceeded("FFT plan size overflow"))?;
        if !Fft::takes_sixstep(plan_n) {
            break;
        }
        plan_n = 1usize << (plan_n.trailing_zeros() / 2);
    }
    if workspace > budget.maximum_workspace_bytes {
        return Err(WelchError::BudgetExceeded(
            "retained factors, FFT workspace and stream input storage",
        ));
    }
    Ok(WelchShape {
        bin_width_hz,
        segments,
        factor_count,
    })
}

fn validate_samples(
    channels: usize,
    samples: usize,
    sample: &impl Fn(usize, usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), WelchError> {
    for c in 0..channels {
        for j in 0..samples {
            if j % 256 == 0 {
                poll(cancelled)?;
            }
            if !sample(c, j).is_finite() {
                return Err(WelchError::NonFinite);
            }
        }
    }
    Ok(())
}

// The batch and chunked owners share this exact segment implementation and
// final-record gain. Chunk boundaries never become new window boundaries.
struct SegmentTransform {
    config: WelchConfig,
    window: Vec<f64>,
    gain: f64,
    fft: Fft,
    buffer: Vec<C64>,
    scratch: Vec<C64>,
}

impl SegmentTransform {
    fn new(
        timestep_s: f64,
        config: WelchConfig,
        segments: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        let n = config.segment_samples;
        let mut window = allocate(n, 0.0)?;
        let mut window_energy = Dd::ZERO;
        for (j, w) in window.iter_mut().enumerate() {
            if j % 256 == 0 {
                poll(cancelled)?;
            }
            *w = match config.window {
                WelchWindow::Rectangular => 1.0,
                WelchWindow::PeriodicHann => {
                    0.5 - 0.5 * det::cos(core::f64::consts::TAU * j as f64 / n as f64)
                }
            };
            let value = Dd::from_f64(*w);
            window_energy = window_energy + value * value;
        }
        let gain =
            det::sqrt(timestep_s) / det::sqrt(segments as f64) / det::sqrt(window_energy.to_f64());
        Ok(Self {
            config,
            window,
            gain,
            fft: Fft::new(n),
            buffer: allocate(n, C64::default())?,
            scratch: allocate(n, C64::default())?,
        })
    }

    fn channel(
        &mut self,
        sample: impl Fn(usize) -> f64,
        mut write: impl FnMut(usize, C64),
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WelchError> {
        poll(cancelled)?;
        let n = self.config.segment_samples;
        let mut scale = 0.0f64;
        for j in 0..n {
            if j % 256 == 0 {
                poll(cancelled)?;
            }
            scale = scale.max(sample(j).abs());
        }
        if scale == 0.0 {
            // A cancelled chunk may have touched these future factor slots.
            // A subsequently accepted zero channel must overwrite all of them.
            for bin in 0..=n / 2 {
                if bin % 256 == 0 {
                    poll(cancelled)?;
                }
                write(bin, C64::default());
            }
            return Ok(());
        }
        let mut mean = Dd::ZERO;
        if self.config.detrend == WelchDetrend::SegmentMean {
            for j in 0..n {
                if j % 256 == 0 {
                    poll(cancelled)?;
                }
                mean = mean + normalized(sample(j), scale)?;
            }
            mean = mean / Dd::from_f64(n as f64);
        }
        for j in 0..n {
            if j % 256 == 0 {
                poll(cancelled)?;
            }
            let centered = normalized(sample(j), scale)? - mean;
            let tapered = (centered * Dd::from_f64(self.window[j])).to_f64();
            if centered.hi != 0.0 && self.window[j] != 0.0 && tapered == 0.0 {
                return Err(WelchError::Unrepresentable);
            }
            self.buffer[j] = C64::new(tapered, 0.0);
        }
        self.fft.forward(&mut self.buffer, &mut self.scratch);
        poll(cancelled)?;
        for (bin, value) in self.buffer.iter().take(n / 2 + 1).enumerate() {
            if bin % 256 == 0 {
                poll(cancelled)?;
            }
            let endpoint = bin == 0 || bin == n / 2;
            let gain = if endpoint {
                self.gain
            } else {
                self.gain * core::f64::consts::SQRT_2
            };
            write(
                bin,
                C64::new(
                    scaled(value.re, gain, scale)?,
                    if endpoint {
                        0.0
                    } else {
                        scaled(-value.im, gain, scale)?
                    },
                ),
            );
        }
        Ok(())
    }
}

fn allocate<T: Clone>(count: usize, value: T) -> Result<Vec<T>, WelchError> {
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| WelchError::Allocation)?;
    out.resize(count, value);
    Ok(out)
}

fn scaled(value: f64, gain: f64, scale: f64) -> Result<f64, WelchError> {
    // DD adds precision, not exponent range. Pair the smallest and largest
    // magnitudes before the middle one: a tiny FFT residual times a tiny
    // spectral gain must not vanish before the physical source scale restores
    // its representable amplitude.
    let mut factors = [value, gain, scale];
    if factors[0].abs() > factors[1].abs() {
        factors.swap(0, 1);
    }
    if factors[1].abs() > factors[2].abs() {
        factors.swap(1, 2);
    }
    if factors[0].abs() > factors[1].abs() {
        factors.swap(0, 1);
    }
    let result =
        ((Dd::from_f64(factors[0]) * Dd::from_f64(factors[2])) * Dd::from_f64(factors[1])).to_f64();
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(WelchError::Unrepresentable)
    } else {
        Ok(result)
    }
}

fn normalized(value: f64, scale: f64) -> Result<Dd, WelchError> {
    // Keep the division's low part until AFTER detrending and tapering.
    // Rounding x/max(|x|) to f64 first can change a representable one-ULP
    // fluctuation on a DC offset by tens of percent; subsequent DD summation
    // cannot recover it. Original sample precision must survive mean removal.
    let result = Dd::from_f64(value) / Dd::from_f64(scale);
    if !result.hi.is_finite() || !result.lo.is_finite() || (value != 0.0 && result.to_f64() == 0.0)
    {
        // A record wider than one f64 FFT's exponent range needs a different
        // representation. Refuse instead of erasing a quiet source component.
        Err(WelchError::Unrepresentable)
    } else {
        Ok(result)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WelchError> {
    if cancelled() {
        Err(WelchError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "welch/tests.rs"]
mod tests;
