//! Causal integration of analytic interference products for linear dynamics.
//!
//! This is an unnormalized physical time integral, not another detector
//! profile. It retains the source clock and duration separately, so a short
//! accepted step need not be representable as start + duration in one f64.

use super::*;

/// Two integrals over the same accepted interval, in observable-units * s.
///
/// For f(t), source start t0, duration h and nonnegative decay lambda:
/// integral = integral_0^h f(t0+s) ds;
/// weighted = integral_0^h exp(-lambda*(h-s)) f(t0+s) ds.
///
/// A thermal RC node, for example, uses weighted/C in its temperature
/// response. The caller owns physical units and whether f represents heat,
/// force or another applied quantity. A general product is not dissipation.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicCausalIntegral {
    integral: Dd,
    weighted_integral: Dd,
    decay: Dd,
    one_minus_decay: Dd,
}

impl HarmonicCausalIntegral {
    /// Unweighted signed input integral, retaining low components.
    pub const fn integral_dd(self) -> Dd {
        self.integral
    }

    /// Causal exponential convolution, retaining low components.
    pub const fn weighted_integral_dd(self) -> Dd {
        self.weighted_integral
    }

    /// exp(-lambda*h), with a retained low component near unity.
    pub const fn decay_dd(self) -> Dd {
        self.decay
    }

    /// 1-exp(-lambda*h), evaluated without a near-unity subtraction.
    pub const fn one_minus_decay_dd(self) -> Dd {
        self.one_minus_decay
    }
}

impl PreparedHarmonicProductSpectrum {
    /// Integrate the complete spectrum and its causal exponential response.
    ///
    /// All arguments must be finite; duration and decay must be nonnegative.
    /// The zero-duration result is exactly zero with unit decay. Zero decay
    /// makes both integrals identical. Work is O(distinct tones), with no
    /// allocation, time grid, Nyquist assumption or dropped carrier terms.
    ///
    /// Closed forms and local divided-difference series avoid subtraction of
    /// nearly equal exponentials and phases. DD algebra retains weak residuals
    /// within its finite dynamic range; transcendental accuracy remains that
    /// of fs_math::det. This is not an interval or null-relative certificate.
    pub fn causal_integral(
        &self,
        start_s: f64,
        duration_s: f64,
        decay_rate_per_s: f64,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        self.causal_integral_cancellable(start_s, duration_s, decay_rate_per_s, || false)
    }

    pub fn causal_integral_cancellable(
        &self,
        start_s: f64,
        duration_s: f64,
        decay_rate_per_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        self.causal_integral_from_clock_dd_cancellable(
            Dd::from_f64(start_s),
            duration_s,
            decay_rate_per_s,
            cancelled,
        )
    }

    /// Retain a physical source epoch plus accepted elapsed time in DD.
    /// This avoids erasing later substeps when an f64 absolute clock cannot
    /// resolve their offset. Duration and decay semantics are unchanged.
    pub fn causal_integral_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
        decay_rate_per_s: f64,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        self.causal_integral_from_clock_dd_cancellable(
            start_s,
            duration_s,
            decay_rate_per_s,
            || false,
        )
    }

    pub fn causal_integral_from_clock_dd_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        decay_rate_per_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        self.causal_integral_from_clock_dd_scaled_cancellable(
            start_s,
            duration_s,
            decay_rate_per_s,
            1.0,
            cancelled,
        )
    }

    /// Integrate after an explicit finite positive output-unit conversion.
    ///
    /// The scale multiplies each DC/tone amplitude BEFORE the existing time
    /// kernels. A small observable in large physical units therefore need not
    /// lose a representable physical integral through an unscaled underflow.
    /// Results are in converted-output-units * s. Decay factors are unchanged.
    ///
    /// Unrepresentable converted amplitudes or nonzero final duration products
    /// that underflow to zero refuse; they are never published as darkness.
    /// Signed cancellation between representable contributions remains valid.
    /// This is allocation-free and does not widen the DD dynamic range or
    /// change the existing kernel/phase accuracy guarantees.
    pub fn causal_integral_from_clock_dd_scaled(
        &self,
        start_s: Dd,
        duration_s: f64,
        decay_rate_per_s: f64,
        output_scale: f64,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        self.causal_integral_from_clock_dd_scaled_cancellable(
            start_s,
            duration_s,
            decay_rate_per_s,
            output_scale,
            || false,
        )
    }

    pub fn causal_integral_from_clock_dd_scaled_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        decay_rate_per_s: f64,
        output_scale: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicCausalIntegral, ProductSpectrumError> {
        spectrum_poll(&mut cancelled)?;
        if !start_s.hi.is_finite()
            || !start_s.lo.is_finite()
            || !duration_s.is_finite()
            || duration_s < 0.0
            || !decay_rate_per_s.is_finite()
            || decay_rate_per_s < 0.0
            || !output_scale.is_finite()
            || output_scale <= 0.0
        {
            return Err(WindowError::InvalidWindow.into());
        }
        let start_s = checked(Dd::from_pair(start_s.hi, start_s.lo))?;
        if duration_s == 0.0 {
            return Ok(HarmonicCausalIntegral {
                integral: Dd::ZERO,
                weighted_integral: Dd::ZERO,
                decay: Dd::ONE,
                one_minus_decay: Dd::ZERO,
            });
        }
        let duration = Dd::from_f64(duration_s);
        let alpha = checked(Dd::from_f64(decay_rate_per_s) * duration)?;
        if decay_rate_per_s > 0.0 && is_zero(alpha) {
            return Err(WindowError::Unrepresentable.into());
        }
        let (decay, one_minus_decay) = decay_parts(alpha)?;
        let zero_decay = decay_rate_per_s == 0.0;
        let dc_gain = if zero_decay {
            duration
        } else {
            checked(duration * checked(one_minus_decay / alpha)?)?
        };
        let output_scale = Dd::from_f64(output_scale);
        let dc = checked_source_product(self.dc, output_scale)?;
        let mut integral = checked_source_product(dc, duration)?;
        let mut weighted_integral = checked_source_product(dc, dc_gain)?;
        for (index, tone) in self.tones.iter().enumerate() {
            if index % 32 == 0 {
                spectrum_poll(&mut cancelled)?;
            }
            let omega = tone.angular_frequency_dd();
            let phase_start = checked(omega * start_s)?;
            let beta = checked(omega * duration)?;
            let (real, imaginary) = proper_phase(phase_start);
            let scaled = [
                checked_source_product(tone.amplitude[0], output_scale)?,
                checked_source_product(tone.amplitude[1], output_scale)?,
            ];
            let amplitude = [
                checked(scaled[0] * real - scaled[1] * imaginary)?,
                checked(scaled[0] * imaginary + scaled[1] * real)?,
            ];
            let (plain, weighted) = integration_kernels(alpha, beta, one_minus_decay)?;
            integral = checked(integral + integrate_source_tone(duration, amplitude, plain)?)?;
            if !zero_decay {
                weighted_integral = checked(
                    weighted_integral + integrate_source_tone(duration, amplitude, weighted)?,
                )?;
            }
        }
        spectrum_poll(&mut cancelled)?;
        if zero_decay {
            weighted_integral = integral;
        }
        Ok(HarmonicCausalIntegral {
            integral,
            weighted_integral,
            decay,
            one_minus_decay,
        })
    }
}

// Prefer physical time gains, whose natural scale is 1/omega or 1/lambda
// for long/stiff intervals. A very short interval can instead underflow one
// gain before two normalized components meet; retain that original ordering
// as a whole-tone fallback. Neither ordering publishes a failed product.
fn integrate_source_tone(
    duration: Dd,
    amplitude: [Dd; 2],
    kernel: [Dd; 2],
) -> Result<Dd, WindowError> {
    let gain_first = || -> Result<Dd, WindowError> {
        let component = |index: usize| -> Result<Dd, WindowError> {
            if is_zero(amplitude[index]) {
                Ok(Dd::ZERO)
            } else {
                checked_source_product(
                    amplitude[index],
                    checked_source_product(duration, kernel[index])?,
                )
            }
        };
        checked(component(0)? - component(1)?)
    };
    match gain_first() {
        Ok(value) => Ok(value),
        Err(WindowError::Unrepresentable) => {
            let normalized = checked(
                checked_source_product(amplitude[0], kernel[0])?
                    - checked_source_product(amplitude[1], kernel[1])?,
            )?;
            checked_source_product(duration, normalized)
        }
        Err(error) => Err(error),
    }
}

// Guard unit conversion and the final time factors, not signed sums: exact
// destructive cancellation between nonzero contributions is still permitted.
fn checked_source_product(left: Dd, right: Dd) -> Result<Dd, WindowError> {
    let product = checked(left * right)?;
    if !is_zero(left) && !is_zero(right) && is_zero(product) {
        Err(WindowError::Unrepresentable)
    } else {
        Ok(product)
    }
}

// Both returned components describe the SAME decay, including its decrement
// below one. The small-argument series keeps alpha² terms lost by f64 expm1.
fn decay_parts(alpha: Dd) -> Result<(Dd, Dd), WindowError> {
    if is_zero(alpha) {
        return Ok((Dd::ONE, Dd::ZERO));
    }
    if alpha.hi <= 0.5 {
        let mut term = -alpha;
        let mut minus_one = term;
        for n in 2..=32 {
            term = checked(term * (-alpha) / Dd::from_f64(f64::from(n)))?;
            minus_one = checked(minus_one + term)?;
        }
        return Ok((checked(Dd::ONE + minus_one)?, -minus_one));
    }
    // Beyond this point the exponential is below the f64 dynamic range.
    // Do not form exp(-alpha.lo) if a huge alpha has a huge low component.
    if alpha.hi > 746.0 {
        return Ok((Dd::ZERO, Dd::ONE));
    }
    let high = Dd::from_f64(det::exp(-alpha.hi));
    let decay = checked(high * checked(Dd::ONE + Dd::from_f64(det::expm1(-alpha.lo)))?)?;
    Ok((decay, checked(Dd::ONE - decay)?))
}

// Kernels normalized by duration:
// plain = integral_0^1 exp(-i beta u) du
// weighted = integral_0^1 exp[-alpha(1-u)-i beta u] du.
fn integration_kernels(
    alpha: Dd,
    beta: Dd,
    one_minus_decay: Dd,
) -> Result<([Dd; 2], [Dd; 2]), WindowError> {
    if alpha.hi <= 0.5 && beta.hi.abs() <= 0.5 {
        // The exponential divided difference has no singularity at alpha=0
        // or beta=0. Expand complete homogeneous powers directly:
        // sum_{n>=0} sum_{k=0}^n (-alpha)^(n-k)(-i beta)^k/(n+1)!.
        // 32 terms suffice well below DD roundoff over this bounded square.
        let mut power = [Dd::ONE, Dd::ZERO];
        let mut term = power;
        let mut plain = power;
        let mut weighted = power;
        for n in 1..=32 {
            let divisor = Dd::from_f64(f64::from(n + 1));
            power = [
                checked(power[1] * beta / divisor)?,
                checked(-power[0] * beta / divisor)?,
            ];
            term = [
                checked(-alpha * term[0] / divisor + power[0])?,
                checked(-alpha * term[1] / divisor + power[1])?,
            ];
            plain = [checked(plain[0] + power[0])?, checked(plain[1] + power[1])?];
            weighted = [
                checked(weighted[0] + term[0])?,
                checked(weighted[1] + term[1])?,
            ];
        }
        return Ok((plain, weighted));
    }
    let half_beta = beta * Dd::from_f64(0.5);
    let (half_real, half_imaginary) = proper_phase(half_beta);
    let sinc_half = sinc(half_beta);
    let plain = [
        checked(half_real * sinc_half)?,
        checked(half_imaginary * sinc_half)?,
    ];
    if is_zero(alpha) {
        return Ok((plain, plain));
    }
    let (phase_real, phase_imaginary) = proper_phase(beta);
    // exp(-i beta)-exp(-alpha), without losing either near-unit decrement.
    let numerator = [
        checked(checked(phase_real - Dd::ONE)? + one_minus_decay)?,
        phase_imaginary,
    ];
    // Divide by alpha-i beta using a scaled quotient, never alpha²+beta².
    let weighted = if alpha.hi >= beta.hi.abs() {
        let ratio = checked(beta / alpha)?;
        let denominator = checked(Dd::ONE + ratio * ratio)?;
        let real = checked(numerator[0] / alpha)?;
        let imaginary = checked(numerator[1] / alpha)?;
        [
            checked(checked(real - ratio * imaginary)? / denominator)?,
            checked(checked(imaginary + ratio * real)? / denominator)?,
        ]
    } else {
        let ratio = checked(alpha / beta)?;
        let denominator = checked(Dd::ONE + ratio * ratio)?;
        let real = checked(numerator[0] / beta)?;
        let imaginary = checked(numerator[1] / beta)?;
        [
            checked(checked(ratio * real - imaginary)? / denominator)?,
            checked(checked(ratio * imaginary + real)? / denominator)?,
        ]
    };
    Ok((plain, weighted))
}

#[cfg(test)]
mod tests;
