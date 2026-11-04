//! Sample the original settled harmonic field on an explicit observer clock.
//!
//! Sampling neither changes propagation nor adds an attack envelope. All
//! physical lines must lie strictly below the output Nyquist frequency. A
//! coarser output clock refuses instead of silently aliasing or dropping tones.

use super::*;
use fs_fft::harmonic_window::phase_hz_from_clock_dd;

/// Uniform output samples at `epoch_s + index / sample_rate_hz`.
///
/// The epoch is on the same physical clock as every input phasor. It is not a
/// propagation delay: each source-to-receiver transfer already includes that.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSampleClock {
    epoch_s: Dd,
    sample_rate_hz: u32,
}

impl HarmonicSampleClock {
    /// Preserve a short sample interval even on a large source-clock epoch.
    pub fn new(epoch_s: Dd, sample_rate_hz: u32) -> Result<Self, AeroacError> {
        if sample_rate_hz == 0 || !epoch_s.hi.is_finite() || !epoch_s.lo.is_finite() {
            return Err(invalid(
                "harmonic sampling needs a finite epoch and positive rate",
            ));
        }
        let clock = Self {
            epoch_s: Dd::from_pair(epoch_s.hi, epoch_s.lo),
            sample_rate_hz,
        };
        clock.validate_index(0)?;
        Ok(clock)
    }

    pub const fn epoch_s(self) -> Dd {
        self.epoch_s
    }
    pub const fn sample_rate_hz(self) -> u32 {
        self.sample_rate_hz
    }

    /// Source-clock time without converting a long integer sample count to f64.
    pub fn time_at(self, index: u64) -> Result<Dd, AeroacError> {
        let integer = Dd::from_f64((index >> 32) as f64) * Dd::from_f64(4_294_967_296.0)
            + Dd::from_f64((index & 0xffff_ffff) as f64);
        let time = self.epoch_s + integer / Dd::from_f64(f64::from(self.sample_rate_hz));
        if !time.hi.is_finite() || !time.lo.is_finite() {
            return Err(nonfinite("harmonic sample clock"));
        }
        Ok(time)
    }

    fn validate_index(self, index: u64) -> Result<(), AeroacError> {
        let time = self.time_at(index)?;
        let step = Dd::ONE / Dd::from_f64(f64::from(self.sample_rate_hz));
        let next = time + step;
        if !next.hi.is_finite() || !next.lo.is_finite() || !time.lt(next) {
            return Err(invalid(
                "harmonic sample interval is unresolved on the source clock",
            ));
        }
        Ok(())
    }
}

/// Per-request work limits; output storage is supplied by the caller.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSampleBudget {
    pub max_samples: usize,
    /// Maximum sample count times the number of distinct physical lines.
    pub max_tone_samples: usize,
}

impl<const D: usize> PreparedHarmonicField<D> {
    /// Admit an entire sample range before a streaming writer changes output.
    /// Empty requests are allowed; an out-of-band line refuses even at a null.
    pub fn validate_sample_range(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        sample_count: u64,
    ) -> Result<(), AeroacError> {
        let end = first_sample
            .checked_add(sample_count)
            .ok_or_else(|| invalid("harmonic sample counter overflow"))?;
        clock.validate_index(first_sample)?;
        clock.validate_index(end)?;
        for line in &self.lines {
            if line.frequency_hz >= 0.5 * f64::from(clock.sample_rate_hz) {
                return Err(invalid(
                    "harmonic field contains a line at or above output Nyquist",
                ));
            }
            // Bounding both ends covers the absolute phase of every sample.
            phase_hz_from_clock_dd(line.frequency_hz, clock.time_at(first_sample)?)
                .map_err(window_error)?;
            phase_hz_from_clock_dd(line.frequency_hz, clock.time_at(end)?).map_err(window_error)?;
        }
        Ok(())
    }

    /// Real pressure [Pa] on the original source clock, with coherent sums
    /// retained in DD until the final scalar. No propagation is recomputed.
    pub fn pressure_at_clock_dd(&self, time_s: Dd) -> Result<f64, AeroacError> {
        self.pressure_at_clock_dd_cancellable(time_s, || false)
    }

    pub fn pressure_at_clock_dd_cancellable(
        &self,
        time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AeroacError> {
        poll(&mut cancelled)?;
        if !time_s.hi.is_finite() || !time_s.lo.is_finite() {
            return Err(nonfinite("harmonic pressure source clock"));
        }
        let mut pressure = Dd::ZERO;
        for (index, line) in self.lines.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let [real, imaginary] =
                phase_hz_from_clock_dd(line.frequency_hz, time_s).map_err(window_error)?;
            pressure = pressure + line.pressure.re * real - line.pressure.im * imaginary;
        }
        poll(&mut cancelled)?;
        let pressure_pa = pressure.to_f64();
        if !pressure.hi.is_finite() || !pressure.lo.is_finite() || !pressure_pa.is_finite() {
            return Err(nonfinite("sampled harmonic pressure"));
        }
        Ok(pressure_pa)
    }

    /// Allocation-free O(samples * distinct lines) observation. Every sample is
    /// evaluated from its absolute integer index, so callback partitioning
    /// cannot reset phase or accumulate oscillator drift. A refused/cancelled
    /// request may have written a prefix into `output`; discard that buffer.
    pub fn sample_pressure_into(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [f64],
        budget: HarmonicSampleBudget,
    ) -> Result<(), AeroacError> {
        self.sample_pressure_into_cancellable(clock, first_sample, output, budget, || false)
    }

    pub fn sample_pressure_into_cancellable(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [f64],
        budget: HarmonicSampleBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(&mut cancelled)?;
        if output.len() > budget.max_samples
            || output
                .len()
                .checked_mul(self.lines.len())
                .is_none_or(|work| work > budget.max_tone_samples)
        {
            return Err(invalid(
                "harmonic pressure sampling exceeds the work budget",
            ));
        }
        let count = u64::try_from(output.len())
            .map_err(|_| invalid("harmonic sample count exceeds u64"))?;
        self.validate_sample_range(clock, first_sample, count)?;
        for (offset, value) in output.iter_mut().enumerate() {
            let time = clock.time_at(first_sample + offset as u64)?;
            *value = self.pressure_at_clock_dd_cancellable(time, &mut cancelled)?;
        }
        poll(&mut cancelled)
    }
}

#[cfg(test)]
mod tests;
