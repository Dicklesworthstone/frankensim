//! Analytic, phase-continuous excitation on the scenario's physical clock.

use crate::ScenarioError;
use fs_math::{dd::Dd, det};
use fs_qty::QtyAny;

// Mathematical 2*pi, including its binary64 tail. Frequencies are authored in
// hertz, so whole cycles must disappear before this conversion to radians.
const TAU: Dd = Dd {
    hi: core::f64::consts::TAU,
    lo: 2.449_293_598_294_706_4e-16,
};

/// A continuous cosine excitation with linearly changing instantaneous hertz.
///
/// `y(t) = offset + amplitude*cos(phase(t))`. On the authored sweep interval,
/// phase integrates the linear frequency exactly in real arithmetic. Outside
/// it, oscillation continues at the respective endpoint frequency with the
/// same accumulated phase; neither the waveform nor its clock is held/reset.
/// Equal endpoint frequencies give an ordinary continuous tone, including DC.
///
/// This is a continuous source, independent of the solver or recording clock.
/// Callers must admit their sample rate using [`Self::max_frequency_hz`] and
/// resolve the changing phase; evaluation does not apply an antialias filter.
/// Evaluation is deterministic, allocation-free and constant-work. Double-
/// double phase products and cycle reduction avoid secular error from a rounded
/// `2*pi` multiplied by a large clock. This is not an interval error bound.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearChirp {
    t_start_s: f64,
    t_end_s: f64,
    frequency_start_hz: f64,
    frequency_end_hz: f64,
    phase_at_start_rad: f64,
    offset: QtyAny,
    amplitude: QtyAny,
}

impl LinearChirp {
    /// Author one finite nonempty sweep, with finite nonnegative frequencies,
    /// finite phase, and finite offset/amplitude of the same dimensions.
    /// Signed amplitudes are permitted. The initial phase belongs to `t_start_s`.
    ///
    /// # Errors
    /// [`ScenarioError::Dimensions`] if offset and amplitude have different
    /// dimensions; [`ScenarioError::Evaluate`] for invalid or nonfinite inputs,
    /// including an unrepresentable sweep duration.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t_start_s: f64,
        t_end_s: f64,
        frequency_start_hz: f64,
        frequency_end_hz: f64,
        phase_at_start_rad: f64,
        offset: QtyAny,
        amplitude: QtyAny,
    ) -> Result<Self, ScenarioError> {
        if !t_start_s.is_finite()
            || !t_end_s.is_finite()
            || t_start_s >= t_end_s
            || !(t_end_s - t_start_s).is_finite()
        {
            return Err(invalid("chirp needs a finite positive sweep duration"));
        }
        if ![frequency_start_hz, frequency_end_hz]
            .into_iter()
            .all(|frequency| frequency.is_finite() && frequency >= 0.0)
            || !phase_at_start_rad.is_finite()
            || !offset.value.is_finite()
            || !amplitude.value.is_finite()
        {
            return Err(invalid(
                "chirp frequencies must be finite and nonnegative; phase, offset and amplitude must be finite",
            ));
        }
        if offset.dims != amplitude.dims {
            return Err(ScenarioError::Dimensions {
                context: "chirp offset and amplitude".to_string(),
                expected: offset.dims.0,
                got: amplitude.dims.0,
            });
        }
        Ok(Self {
            t_start_s,
            t_end_s,
            frequency_start_hz,
            frequency_end_hz,
            phase_at_start_rad,
            offset,
            amplitude,
        })
    }

    /// Start of the linear frequency sweep, in seconds on the authored clock.
    #[must_use]
    pub const fn t_start_s(&self) -> f64 {
        self.t_start_s
    }
    /// End of the linear frequency sweep, in seconds on the same clock.
    #[must_use]
    pub const fn t_end_s(&self) -> f64 {
        self.t_end_s
    }
    /// Instantaneous frequency at and before the sweep start, in hertz.
    #[must_use]
    pub const fn frequency_start_hz(&self) -> f64 {
        self.frequency_start_hz
    }
    /// Instantaneous frequency at and after the sweep end, in hertz.
    #[must_use]
    pub const fn frequency_end_hz(&self) -> f64 {
        self.frequency_end_hz
    }
    /// Authored cosine phase at [`Self::t_start_s`], in radians, without wrapping.
    #[must_use]
    pub const fn phase_at_start_rad(&self) -> f64 {
        self.phase_at_start_rad
    }
    /// Constant offset in coherent SI, retaining the signal's dimensions.
    #[must_use]
    pub const fn offset(&self) -> QtyAny {
        self.offset
    }
    /// Signed peak cosine amplitude in the same dimensions as the offset.
    #[must_use]
    pub const fn amplitude(&self) -> QtyAny {
        self.amplitude
    }

    /// Instantaneous physical frequency in hertz at a time in seconds,
    /// including both constant-frequency tails.
    ///
    /// # Errors
    /// [`ScenarioError::Evaluate`] for nonfinite time or unrepresentable
    /// frequency interpolation.
    pub fn instantaneous_frequency_hz(&self, time: f64) -> Result<f64, ScenarioError> {
        finite_time(time)?;
        if time <= self.t_start_s {
            return Ok(self.frequency_start_hz);
        }
        if time >= self.t_end_s {
            return Ok(self.frequency_end_hz);
        }
        let elapsed = Dd::from_f64(time) - Dd::from_f64(self.t_start_s);
        let duration = Dd::from_f64(self.t_end_s) - Dd::from_f64(self.t_start_s);
        let delta = Dd::from_f64(self.frequency_end_hz) - Dd::from_f64(self.frequency_start_hz);
        finite((Dd::from_f64(self.frequency_start_hz) + delta * (elapsed / duration)).to_f64())
    }

    /// Maximum instantaneous frequency over a finite closed interval.
    /// A point interval is admitted; reversed endpoints refuse. The linear
    /// sweep and constant tails are monotone, so endpoints determine the maximum.
    /// Both query endpoints are seconds on the authored clock.
    ///
    /// # Errors
    /// [`ScenarioError::Evaluate`] for nonfinite or reversed query endpoints,
    /// or unrepresentable frequency interpolation.
    pub fn max_frequency_hz(&self, start: f64, end: f64) -> Result<f64, ScenarioError> {
        if start > end {
            return Err(invalid("chirp frequency interval must be ordered"));
        }
        Ok(self
            .instantaneous_frequency_hz(start)?
            .max(self.instantaneous_frequency_hz(end)?))
    }

    /// Sample the same continuous excitation at a finite physical time.
    /// Time is in seconds; the returned value retains the offset's dimensions.
    ///
    /// # Errors
    /// [`ScenarioError::Evaluate`] for nonfinite time, unrepresentable phase
    /// arithmetic, or an unrepresentable output value.
    pub fn eval(&self, time: f64) -> Result<QtyAny, ScenarioError> {
        finite_time(time)?;
        if self.amplitude.value == 0.0 {
            return Ok(self.offset);
        }
        let (_, cosine) = self.phase_sin_cos(time)?;
        let value = finite(
            (Dd::from_f64(self.offset.value) + Dd::from_f64(self.amplitude.value) * cosine)
                .to_f64(),
        )?;
        Ok(QtyAny::new(value, self.offset.dims))
    }

    /// Return value, first derivative and second derivative in coherent SI.
    ///
    /// Derivatives divide the value's dimensions by seconds and seconds squared.
    /// At the sweep endpoints the frequency is continuous but its slope jumps.
    /// For a nonconstant sweep this operation refuses those exact knots instead
    /// of assigning one side's second derivative. Value-only evaluation remains
    /// valid there. No numerical differentiation or extra sample clock is used.
    /// Zero amplitude is the identically constant signal and has zero
    /// derivatives even at the sweep endpoints.
    ///
    /// # Errors
    /// [`ScenarioError::Evaluate`] for nonfinite time, a frequency-slope knot
    /// of a nonzero-amplitude sweep, or unrepresentable phase/value/derivatives.
    pub fn value_derivatives(&self, time: f64) -> Result<(f64, f64, f64), ScenarioError> {
        finite_time(time)?;
        if self.amplitude.value == 0.0 {
            return Ok((self.offset.value, 0.0, 0.0));
        }
        let changing = self.frequency_start_hz != self.frequency_end_hz;
        if changing && (time == self.t_start_s || time == self.t_end_s) {
            return Err(invalid(
                "chirp second derivative is undefined at a frequency-slope knot",
            ));
        }
        let (sine, cosine) = self.phase_sin_cos(time)?;
        let amplitude = Dd::from_f64(self.amplitude.value);
        let omega = checked(Dd::from_f64(self.instantaneous_frequency_hz(time)?) * TAU)?;
        let slope = if changing && time > self.t_start_s && time < self.t_end_s {
            checked(
                ((Dd::from_f64(self.frequency_end_hz) - Dd::from_f64(self.frequency_start_hz))
                    / (Dd::from_f64(self.t_end_s) - Dd::from_f64(self.t_start_s)))
                    * TAU,
            )?
        } else {
            Dd::ZERO
        };
        let value = finite((Dd::from_f64(self.offset.value) + amplitude * cosine).to_f64())?;
        let first = finite((-(amplitude * sine) * omega).to_f64())?;
        let second =
            finite((-(amplitude * cosine) * omega * omega - amplitude * sine * slope).to_f64())?;
        Ok((value, first, second))
    }

    fn phase_sin_cos(&self, time: f64) -> Result<(Dd, Dd), ScenarioError> {
        let f0 = Dd::from_f64(self.frequency_start_hz);
        let f1 = Dd::from_f64(self.frequency_end_hz);
        let cycles = if time <= self.t_start_s || self.frequency_start_hz == self.frequency_end_hz {
            // DC has no clock product to overflow, even at remote finite times.
            if self.frequency_start_hz == 0.0 {
                Dd::ZERO
            } else {
                reduced(checked(
                    f0 * (Dd::from_f64(time) - Dd::from_f64(self.t_start_s)),
                )?)
            }
        } else {
            let duration = Dd::from_f64(self.t_end_s) - Dd::from_f64(self.t_start_s);
            if time < self.t_end_s {
                let elapsed = Dd::from_f64(time) - Dd::from_f64(self.t_start_s);
                // Integral of f0 + (f1-f0)*u/T. Form u/T first so no
                // unnecessary squared timestamp or unbounded chirp rate occurs.
                let mean_frequency = f0 + (f1 - f0) * (elapsed / duration) * Dd::from_f64(0.5);
                reduced(checked(mean_frequency * elapsed)?)
            } else {
                // Reduce each interval before adding it to the other. A long
                // endpoint tail must not erase the sweep's retained phase.
                let mean_frequency = f0 + (f1 - f0) * Dd::from_f64(0.5);
                let sweep = reduced(checked(mean_frequency * duration)?);
                let tail = if self.frequency_end_hz == 0.0 {
                    Dd::ZERO
                } else {
                    reduced(checked(
                        f1 * (Dd::from_f64(time) - Dd::from_f64(self.t_end_s)),
                    )?)
                };
                reduced(sweep + tail)
            }
        };
        let phase = checked(Dd::from_f64(self.phase_at_start_rad) + cycles * TAU)?;
        let (sine, cosine) = det::sin_cos_dd(phase);
        Ok((checked(sine)?, checked(cosine)?))
    }
}

// Both DD components may contain whole cycles. Removing only hi.round()
// leaves a large low component for very long clocks and reintroduces drift.
fn reduced(cycles: Dd) -> Dd {
    let fraction =
        Dd::from_f64(cycles.hi - cycles.hi.round()) + Dd::from_f64(cycles.lo - cycles.lo.round());
    fraction - Dd::from_f64(fraction.hi.round())
}

fn checked(value: Dd) -> Result<Dd, ScenarioError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(invalid("chirp phase or derivative is not representable"))
    }
}

fn finite(value: f64) -> Result<f64, ScenarioError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("chirp value or derivative is not representable"))
    }
}

fn finite_time(time: f64) -> Result<(), ScenarioError> {
    if time.is_finite() {
        Ok(())
    } else {
        Err(invalid("chirp evaluation time must be finite"))
    }
}

fn invalid(what: &str) -> ScenarioError {
    ScenarioError::Evaluate {
        what: what.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::Dims;

    fn chirp(start: f64, end: f64, f0: f64, f1: f64, phase: f64) -> LinearChirp {
        LinearChirp::new(
            start,
            end,
            f0,
            f1,
            phase,
            QtyAny::dimensionless(0.0),
            QtyAny::dimensionless(1.0),
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_chirp_integrates_frequency_and_continues_phase_in_both_tails() {
        let source = chirp(2.0, 5.0, 1.25, 3.5, 0.31);
        for t in [-1.0, 1.875, 2.0, 2.125, 3.375, 5.0, 5.125, 9.25] {
            let elapsed = t - 2.0;
            let cycles = if t < 2.0 {
                1.25 * elapsed
            } else if t <= 5.0 {
                1.25 * elapsed + 0.375 * elapsed * elapsed
            } else {
                7.125 + 3.5 * (t - 5.0)
            };
            close(
                source.eval(t).unwrap().value,
                (0.31 + core::f64::consts::TAU * cycles).cos(),
                2e-14,
            );
            let expected_hz = 1.25 + 0.75 * elapsed.clamp(0.0, 3.0);
            close(
                source.instantaneous_frequency_hz(t).unwrap(),
                expected_hz,
                1e-15,
            );
        }
        close(source.eval(2.0).unwrap().value, 0.31_f64.cos(), 2e-16);
        assert_ne!(
            source.eval(5.125).unwrap().value,
            source.eval(5.0).unwrap().value
        );
        assert_eq!(source.max_frequency_hz(-8.0, 2.0).unwrap(), 1.25);
        assert_eq!(source.max_frequency_hz(0.0, 9.0).unwrap(), 3.5);
        let descending = chirp(2.0, 5.0, 3.5, 1.25, 0.0);
        assert_eq!(descending.max_frequency_hz(3.0, 9.0).unwrap(), 2.75);
    }

    #[test]
    fn g3_chirp_is_independent_of_sample_chunks_and_absolute_epoch() {
        let source = chirp(2.0, 5.0, 1.25, 3.5, -0.4);
        let epoch = 2.0_f64.powi(36);
        let shifted = chirp(2.0 + epoch, 5.0 + epoch, 1.25, 3.5, -0.4);
        for t in [-1.25, 2.125, 4.75, 5.0, 9.25] {
            assert_eq!(source.eval(t).unwrap(), shifted.eval(t + epoch).unwrap());
        }
        let before = source.eval(3.25).unwrap();
        let _ = source.eval(-100.0).unwrap();
        let _ = source.eval(1e9).unwrap();
        assert_eq!(source.eval(3.25).unwrap(), before);
    }

    #[test]
    fn g1_tone_limit_preserves_fractional_cycles_on_a_long_clock() {
        let frequency = 1_000_000.25;
        let source = chirp(0.0, 1.0, frequency, frequency, 0.0);
        let t = 2.0_f64.powi(36) + 0.125;
        // The large integer clock contributes only whole cycles; the exact
        // fractional sample contributes 125000 + 1/32 cycles.
        close(
            source.eval(t).unwrap().value,
            (core::f64::consts::TAU / 32.0).cos(),
            3e-16,
        );
        assert_eq!(source.eval(2.0_f64.powi(36)).unwrap().value, 1.0);
        let dc = chirp(-1e308, -0.5e308, 0.0, 0.0, 0.37);
        close(dc.eval(1e308).unwrap().value, 0.37_f64.cos(), 2e-16);
        let constant_tone = chirp(1.0, 3.0, 0.5, 0.5, 0.0);
        for t in [-3.0, 1.0, 1.5, 3.0, 5.75] {
            close(
                constant_tone.eval(t).unwrap().value,
                (core::f64::consts::PI * (t - 1.0)).cos(),
                3e-15,
            );
            assert!(constant_tone.value_derivatives(t).is_ok());
        }
    }

    #[test]
    fn g1_analytic_derivatives_match_the_integrated_phase_law() {
        let source = chirp(-1.0, 2.0, 0.25, 1.75, 0.27);
        for t in [-2.0_f64, -0.75, 0.125, 1.5, 3.0] {
            let (value, first, second) = source.value_derivatives(t).unwrap();
            assert_eq!(source.eval(t).unwrap().value, value);
            let h = 1e-5;
            let before = source.eval(t - h).unwrap().value;
            let after = source.eval(t + h).unwrap().value;
            close(first, (after - before) / (2.0 * h), 3e-8);
            close(second, (after - 2.0 * value + before) / (h * h), 2e-5);
        }
        assert!(source.value_derivatives(-1.0).is_err());
        assert!(source.value_derivatives(2.0).is_err());
        assert!(source.eval(-1.0).is_ok());
        assert!(source.eval(2.0).is_ok());
    }

    #[test]
    fn g0_zero_amplitude_is_constant_at_sweep_knots_and_remote_times() {
        let source = LinearChirp::new(
            0.0,
            1.0,
            f64::MAX,
            0.0,
            f64::MAX,
            QtyAny::dimensionless(7.0),
            QtyAny::dimensionless(0.0),
        )
        .unwrap();
        for time in [-f64::MAX, 0.0, 0.25, 1.0, f64::MAX] {
            assert_eq!(source.eval(time).unwrap(), QtyAny::dimensionless(7.0));
            assert_eq!(source.value_derivatives(time).unwrap(), (7.0, 0.0, 0.0));
        }
        // Frequency queries retain the authored source clock even when its
        // amplitude is zero; they do not infer a different physical waveform.
        assert_eq!(source.instantaneous_frequency_hz(0.0).unwrap(), f64::MAX);
        assert!(source.eval(f64::NAN).is_err());
        assert!(source.value_derivatives(f64::INFINITY).is_err());
    }

    #[test]
    fn g0_chirp_admits_units_and_refuses_invalid_clocks_or_values() {
        let force = Dims([1, 1, -2, 0, 0, 0]);
        let source = LinearChirp::new(
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            QtyAny::new(3.0, force),
            QtyAny::new(-2.0, force),
        )
        .unwrap();
        assert_eq!(source.eval(0.0).unwrap(), QtyAny::new(1.0, force));
        assert!(
            LinearChirp::new(
                0.0,
                1.0,
                0.0,
                1.0,
                0.0,
                QtyAny::new(3.0, force),
                QtyAny::dimensionless(2.0),
            )
            .is_err()
        );
        for [a, b, f0, f1, phase] in [
            [0.0, 0.0, 1.0, 2.0, 0.0],
            [-f64::MAX, f64::MAX, 1.0, 2.0, 0.0],
            [0.0, 1.0, -1.0, 2.0, 0.0],
            [0.0, 1.0, 1.0, f64::INFINITY, 0.0],
            [0.0, 1.0, 1.0, 2.0, f64::NAN],
        ] {
            assert!(
                LinearChirp::new(
                    a,
                    b,
                    f0,
                    f1,
                    phase,
                    QtyAny::dimensionless(0.0),
                    QtyAny::dimensionless(1.0)
                )
                .is_err()
            );
        }
        assert!(source.eval(f64::NAN).is_err());
        assert!(source.max_frequency_hz(2.0, 1.0).is_err());
        assert!(source.max_frequency_hz(0.0, f64::INFINITY).is_err());
        assert_eq!(source.max_frequency_hz(0.5, 0.5).unwrap(), 0.5);
        let unrepresentable = chirp(0.0, 1.0, f64::MAX, f64::MAX, 0.0);
        assert!(unrepresentable.eval(2.0).is_err());
    }
}
