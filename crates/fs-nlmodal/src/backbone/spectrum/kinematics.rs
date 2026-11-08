//! Physical derivatives of the same conservative Jacobi spectrum.

use super::{DuffingBackboneError, DuffingSpectrum, binary_parts, denominator, prepare_spectrum};
use fs_math::qd::Qd;

/// Displacement, velocity and acceleration harmonics at one physical turning point.
#[derive(Debug, Clone, PartialEq)]
pub struct DuffingKinematicSpectrum {
    /// The existing displacement spectrum, including its frequency and tail.
    pub displacement: DuffingSpectrum,
    /// Signed peak sine amplitudes in original displacement units per second.
    /// Entry `j` multiplies `sin((2*j+1)*omega*t)` and equals
    /// `-(2*j+1)*omega` times the exact displacement cosine coefficient.
    pub velocity_sine_amplitudes: Vec<f64>,
    /// Signed peak cosine amplitudes in original displacement units per second².
    /// Entry `j` equals `-((2*j+1)*omega)²` times the exact displacement
    /// cosine coefficient. These are physical acceleration, not phase derivatives.
    pub acceleration_cosine_amplitudes: Vec<f64>,
    /// Uniform omitted-velocity absolute-sum bound in physical units per second.
    /// This estimates harmonic truncation, excluding floating-point evaluation
    /// error; it is not a validated interval enclosure.
    pub velocity_tail_bound: f64,
    /// Corresponding omitted-acceleration bound in physical units per second².
    /// Resolving displacement alone does not establish either derivative bound.
    pub acceleration_tail_bound: f64,
}

/// Physical kinematic Fourier seed of the original conservative Duffing orbit.
///
/// This shares the displacement query's original parameters, exact Jacobi
/// series, bounded work budget and invariant phase convention. The scalar
/// frequency/nome preparation runs once. Velocity is a sine series and
/// acceleration a cosine series; both have the sign opposite the displacement
/// coefficients. For an `exp(+i*h*omega*t)` convention, complex velocity is
/// `-i/2` times the reported sine amplitude.
///
/// Derivatives are formed from binary-scaled source amplitudes before rounding
/// to displacement coefficients. Thus a representable acceleration harmonic
/// survives even if its displacement or velocity rounds below the subnormal
/// range, and no intermediate `omega²` or `h*omega` overflow erases a valid
/// result. A derivative coefficient or bound that cannot fit `f64` refuses.
/// The existing displacement-only query retains its wider admission range.
///
/// Omitted tails use closed sums of `h*q^j` and `h²*q^j`, with the same
/// denominator bound as the displacement owner. A nonlinear tail that rounds
/// below the subnormal range is reported as the least positive subnormal.
/// Retaining more harmonics does not renormalize any coefficient. These are
/// conservative reference/solver seeds, not damping or coupled stability claims.
pub fn duffing_kinematic_spectrum_exact(
    omega0: f64,
    beta: f64,
    amplitude: f64,
    odd_harmonics: usize,
) -> Result<DuffingKinematicSpectrum, DuffingBackboneError> {
    let prepared = prepare_spectrum(omega0, beta, amplitude, odd_harmonics)?;
    let displacement = prepared.displacement(odd_harmonics)?;
    let omega = prepared.angular_frequency;
    let mut numerator = ScaledPositive::new(amplitude).times(prepared.first_ratio);
    let mut velocity_sine_amplitudes = Vec::with_capacity(odd_harmonics);
    let mut acceleration_cosine_amplitudes = Vec::with_capacity(odd_harmonics);
    for j in 0..odd_harmonics {
        let order = 2 * j + 1;
        let inverse_denominator = 1.0 / denominator(prepared.decay, order, prepared.softening);
        let velocity = numerator.times(omega).times(order as f64);
        let acceleration = velocity.times(omega).times(order as f64);
        let sign = if prepared.softening && j % 2 != 0 {
            1.0
        } else {
            -1.0
        };
        velocity_sine_amplitudes.push(sign * velocity.times(inverse_denominator).value()?);
        acceleration_cosine_amplitudes
            .push(sign * acceleration.times(inverse_denominator).value()?);
        numerator = numerator
            .times(prepared.root_nome)
            .times(prepared.root_nome);
    }
    let (velocity_tail_bound, acceleration_tail_bound) = if prepared.root_nome == 0.0 {
        (0.0, 0.0)
    } else {
        let first = (2 * odd_harmonics + 1) as f64;
        let nome = prepared.root_nome * prepared.root_nome;
        let gap = prepared.one_minus_nome;
        let tail_denominator = if prepared.softening {
            denominator(prepared.decay, 2 * odd_harmonics + 1, true)
        } else {
            1.0
        };
        // Weighted geometric sums: E[h] = first + 2*q/(1-q), and
        // E[h²] = E[h]² + 4*q/(1-q)² for j starting at the first omitted term.
        let mean_order = first + 2.0 * nome / gap;
        let mean_square = mean_order * mean_order + 4.0 * nome / (gap * gap);
        let tail = numerator.times(1.0 / gap).times(1.0 / tail_denominator);
        (
            tail.times(omega).times(mean_order).upper_bound()?,
            tail.times(omega)
                .times(omega)
                .times(mean_square)
                .upper_bound()?,
        )
    };
    Ok(DuffingKinematicSpectrum {
        displacement,
        velocity_sine_amplitudes,
        acceleration_cosine_amplitudes,
        velocity_tail_bound,
        acceleration_tail_bound,
    })
}

/// Carry the exponent through every harmonic before choosing a physical unit.
/// All factors come from admitted positive source/series data, or exact zero.
#[derive(Clone, Copy)]
struct ScaledPositive {
    mantissa: f64,
    exponent: i32,
}

impl ScaledPositive {
    fn new(value: f64) -> Self {
        if value == 0.0 {
            Self {
                mantissa: 0.0,
                exponent: 0,
            }
        } else {
            let (mantissa, exponent) = binary_parts(value);
            Self { mantissa, exponent }
        }
    }

    fn times(mut self, factor: f64) -> Self {
        if self.mantissa == 0.0 || factor == 0.0 {
            return Self::new(0.0);
        }
        let (mantissa, exponent) = binary_parts(factor);
        self.mantissa *= mantissa;
        self.exponent += exponent;
        if self.mantissa >= 2.0 {
            self.mantissa *= 0.5;
            self.exponent += 1;
        }
        self
    }

    fn value(self) -> Result<f64, DuffingBackboneError> {
        if self.mantissa == 0.0 || self.exponent < -1075 {
            return Ok(0.0);
        }
        if self.exponent > 1023 {
            return Err(DuffingBackboneError::NumericalRange);
        }
        let value = Qd::from_f64(self.mantissa)
            .scale_power_of_two(self.exponent)
            .to_f64();
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DuffingBackboneError::NumericalRange)
        }
    }

    fn upper_bound(self) -> Result<f64, DuffingBackboneError> {
        let value = self.value()?.next_up();
        if value.is_finite() {
            Ok(value)
        } else {
            Err(DuffingBackboneError::NumericalRange)
        }
    }
}

#[cfg(test)]
mod tests;
