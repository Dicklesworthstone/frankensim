//! Endpoint values and changes on the original source clock.
//!
//! A storage increment is an endpoint difference, not a normalized detector
//! average or an interval integral. Constant storage cancels algebraically;
//! short durations survive even when start + duration rounds back to start.

use super::*;

impl PreparedHarmonicProductSpectrum {
    /// Evaluate the real observable without rounding its source clock or result
    /// to one f64. Units are those of the original quadratic observable.
    pub fn evaluate_from_clock_dd(&self, time_s: Dd) -> Result<Dd, ProductSpectrumError> {
        self.evaluate_from_clock_dd_cancellable(time_s, || false)
    }

    pub fn evaluate_from_clock_dd_cancellable(
        &self,
        time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        let time = source_time(time_s)?;
        let mut sum = self.dc;
        for (index, tone) in self.tones.iter().enumerate() {
            if index % 64 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let angle = checked(tone.angular_frequency_rad_s * time)?;
            let (real, imaginary) = proper_phase(angle);
            sum =
                checked(sum + checked(tone.amplitude[0] * real - tone.amplitude[1] * imaginary)?)?;
        }
        spectrum_poll(&mut cancelled)?;
        Ok(sum)
    }

    /// Compute f(start + duration) - f(start), in the original observable units.
    /// Duration must be finite and nonnegative. Zero duration is exactly zero.
    ///
    /// Each tone uses exp(-i omega start) * exp(-i omega duration/2) and
    /// -2i sin(omega duration/2). It never subtracts two rounded endpoint
    /// values or adds duration to an f64 absolute time. DC cancels exactly.
    /// This is allocation-free, O(distinct tones), and uses the original
    /// deterministic DD phase reduction. It is not an error enclosure.
    pub fn increment_dd(&self, start_s: f64, duration_s: f64) -> Result<Dd, ProductSpectrumError> {
        self.increment_from_clock_dd(Dd::from_f64(start_s), duration_s)
    }

    pub fn increment_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
    ) -> Result<Dd, ProductSpectrumError> {
        self.increment_from_clock_dd_cancellable(start_s, duration_s, || false)
    }

    /// Preserve a DD source epoch and cancellation before returning any result.
    /// A nonzero phase increment that underflows refuses rather than returning
    /// a false constant. DD retains cancellation within its finite range.
    pub fn increment_from_clock_dd_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        let start = source_time(start_s)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(WindowError::InvalidWindow.into());
        }
        if duration_s == 0.0 {
            return Ok(Dd::ZERO);
        }
        let duration = Dd::from_f64(duration_s);
        let half = Dd::from_f64(0.5);
        let mut change = Dd::ZERO;
        for (index, tone) in self.tones.iter().enumerate() {
            if index % 64 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let omega = tone.angular_frequency_rad_s;
            let phase = checked(omega * start)?;
            let width = checked(omega * duration)?;
            let half_width = checked(width * half)?;
            if is_zero(width) || is_zero(half_width) {
                return Err(WindowError::Unrepresentable.into());
            }
            let (real, imaginary) = proper_phase(phase);
            let rotated_real = checked(tone.amplitude[0] * real - tone.amplitude[1] * imaginary)?;
            let rotated_imag = checked(tone.amplitude[0] * imaginary + tone.amplitude[1] * real)?;
            let (half_real, half_imag) = proper_phase(half_width);
            let midpoint_imag = checked(rotated_real * half_imag + rotated_imag * half_real)?;
            let increment = checked(Dd::from_f64(-2.0) * half_imag * midpoint_imag)?;
            change = checked(change + increment)?;
        }
        spectrum_poll(&mut cancelled)?;
        Ok(change)
    }
}

fn source_time(time_s: Dd) -> Result<Dd, ProductSpectrumError> {
    if !time_s.hi.is_finite() || !time_s.lo.is_finite() {
        return Err(WindowError::InvalidWindow.into());
    }
    Ok(checked(Dd::from_pair(time_s.hi, time_s.lo))?)
}

#[cfg(test)]
mod tests;
