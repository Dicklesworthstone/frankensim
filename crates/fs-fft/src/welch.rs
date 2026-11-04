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

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn estimate_with(
        channels: usize,
        samples: usize,
        sample: impl Fn(usize, usize) -> f64,
        timestep_s: f64,
        config: WelchConfig,
        budget: WelchBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WelchError> {
        poll(&mut cancelled)?;
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
        let bin_count = n / 2 + 1;
        let factor_count = bin_count
            .checked_mul(transforms)
            .ok_or(WelchError::BudgetExceeded("factor count overflow"))?;
        let mut workspace = factor_count
            .checked_mul(core::mem::size_of::<C64>())
            .and_then(|v| v.checked_add(n * (8 + 2 * core::mem::size_of::<C64>())))
            .ok_or(WelchError::BudgetExceeded("workspace size overflow"))?;
        // Account for the actual FFT plan, including optional six-step subplans.
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
                "retained factors and FFT workspace",
            ));
        }
        // Validate the entire supplied record, including its explicitly unused
        // tail, before any estimate or allocation can be published.
        for c in 0..channels {
            for j in 0..samples {
                if j % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if !sample(c, j).is_finite() {
                    return Err(WelchError::NonFinite);
                }
            }
        }
        let mut window = allocate(n, 0.0)?;
        let mut window_energy = Dd::ZERO;
        for (j, w) in window.iter_mut().enumerate() {
            if j % 256 == 0 {
                poll(&mut cancelled)?;
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
        let fft = Fft::new(n);
        let mut buffer = allocate(n, C64::default())?;
        let mut scratch = allocate(n, C64::default())?;
        let mut factors = allocate(factor_count, C64::default())?;
        for segment in 0..segments {
            let start = segment * config.hop_samples;
            for channel in 0..channels {
                poll(&mut cancelled)?;
                let mut scale = 0.0f64;
                for j in 0..n {
                    if j % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    scale = scale.max(sample(channel, start + j).abs());
                }
                if scale == 0.0 {
                    continue;
                }
                let mut mean = Dd::ZERO;
                if config.detrend == WelchDetrend::SegmentMean {
                    for j in 0..n {
                        if j % 256 == 0 {
                            poll(&mut cancelled)?;
                        }
                        mean = mean + normalized(sample(channel, start + j), scale)?;
                    }
                    mean = mean / Dd::from_f64(n as f64);
                }
                for j in 0..n {
                    if j % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let centered = normalized(sample(channel, start + j), scale)? - mean;
                    let tapered = (centered * Dd::from_f64(window[j])).to_f64();
                    if centered.hi != 0.0 && window[j] != 0.0 && tapered == 0.0 {
                        return Err(WelchError::Unrepresentable);
                    }
                    buffer[j] = C64::new(tapered, 0.0);
                }
                fft.forward(&mut buffer, &mut scratch);
                poll(&mut cancelled)?;
                for (bin, value) in buffer.iter().take(bin_count).enumerate() {
                    if bin % 256 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let endpoint = bin == 0 || bin == n / 2;
                    let factor_gain = if endpoint {
                        gain
                    } else {
                        gain * core::f64::consts::SQRT_2
                    };
                    let offset = (bin * segments + segment) * channels + channel;
                    factors[offset] = C64::new(
                        scaled(value.re, factor_gain, scale)?,
                        if endpoint {
                            0.0
                        } else {
                            scaled(-value.im, factor_gain, scale)?
                        },
                    );
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            config,
            timestep_s,
            bin_width_hz,
            channels,
            segments,
            input_samples: samples,
            used_samples: (segments - 1) * config.hop_samples + n,
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
