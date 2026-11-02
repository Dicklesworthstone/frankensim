//! A physical input-energy excursion bound for unresolved temporal forcing.
//!
//! The spectrum's primitive is compared with its endpoint chord. This is the
//! input error made by replacing a source by its exact interval-average value;
//! an endpoint energy match alone cannot reveal oscillations between endpoints.

use super::*;

impl PreparedHarmonicProductSpectrum {
    /// Bound the deviation of accumulated input from constant-mean input.
    ///
    /// For F(s)=integral_0^s f(t0+u) du, h=duration and 0<=s<=h, returns B
    /// such that |F(s)-(s/h)F(h)| <= B for the analytic retained spectrum,
    /// independently of the source-clock start t0. Units are f-units * seconds.
    /// At h=0 the bound is zero. DC contributes exactly zero.
    ///
    /// One positive-frequency tone with amplitude envelope A and frequency w
    /// contributes at most min(A*w*h²/8, A*h/2, 2*A/w). The first follows from
    /// the linear interpolation remainder for the primitive (|F''|<=A*w);
    /// the others use bounded input and the harmonic primitive, respectively.
    /// See NIST DLMF 3.3.E5 for the interpolation remainder.
    ///
    /// The envelope |Re amplitude|+|Im amplitude| is conservatively bounded
    /// from all DD components. Every positive arithmetic operation rounds
    /// outward by one f64 successor. The bound therefore encloses this exact
    /// retained analytic spectrum, with no sampled extrema or phase search.
    /// It does not enclose source/propagation/model errors, and does not by
    /// itself certify a nonlinear state integration or its admissibility.
    ///
    /// Work is O(distinct tones), allocation-free, with bounded cancellation.
    pub fn input_ripple_bound(&self, duration_s: f64) -> Result<f64, ProductSpectrumError> {
        self.input_ripple_bound_cancellable(duration_s, || false)
    }

    pub fn input_ripple_bound_cancellable(
        &self,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(WindowError::InvalidWindow.into());
        }
        if duration_s == 0.0 {
            return Ok(0.0);
        }
        let half_duration = mul_up(duration_s, 0.5);
        let eighth_duration = mul_up(duration_s, 0.125);
        let mut bound = 0.0;
        for (index, tone) in self.tones.iter().enumerate() {
            if index % 32 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let frequency = tone.angular_frequency_dd();
            let low = frequency.lo.abs();
            let upper_frequency = add_up(frequency.hi, low);
            let lower_frequency = if low == 0.0 {
                frequency.hi
            } else {
                (frequency.hi - low).next_down().max(0.0)
            };
            let phase_upper = mul_up(upper_frequency, duration_s);
            // If w*h >=4, the derivative estimate cannot improve h/2.
            // This branch avoids overflow in a discarded h² intermediate.
            let derivative = if phase_upper < 4.0 {
                mul_up(phase_upper, eighth_duration)
            } else {
                f64::INFINITY
            };
            let primitive = if lower_frequency > 0.0 {
                div_up(2.0, lower_frequency)
            } else {
                f64::INFINITY
            };
            let gain = derivative.min(half_duration).min(primitive);
            // Scale each component first: their envelope could overflow even
            // when its product with a very short primitive gain is finite.
            for amplitude in tone.amplitude {
                bound = add_up(bound, mul_up(amplitude.hi.abs(), gain));
                bound = add_up(bound, mul_up(amplitude.lo.abs(), gain));
            }
        }
        spectrum_poll(&mut cancelled)?;
        if !bound.is_finite() {
            return Err(WindowError::Unrepresentable.into());
        }
        Ok(bound)
    }
}

// All operands are nonnegative. Exact zero shortcuts avoid inventing a
// positive bound for a genuinely dark term; positive underflow rounds upward
// to the smallest subnormal instead of silently returning zero.
fn add_up(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        (a + b).next_up()
    }
}

fn mul_up(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        (a * b).next_up()
    }
}

fn div_up(a: f64, b: f64) -> f64 {
    (a / b).next_up()
}

#[cfg(test)]
mod tests;
