//! Odd harmonics of the original conservative cubic spring.
//!
//! Hardening motion is `A cn(lambda*t,k)`, with
//! `lambda^2 = omega0^2 + beta*A^2`, `k^2 = beta*A^2/(2*lambda^2)`.
//! Softening motion inside the well is `A cd(lambda*t,k)`, with
//! `lambda^2 = omega0^2 + beta*A^2/2`, `k^2 = -beta*A^2/(2*lambda^2)`.
//! NIST DLMF 22.11.2 and 22.11.4 give their cosine series; 22.2.1 defines
//! the nome. The softening coefficients alternate in sign. Both represent
//! a physical turning point at phase zero, not a prescribed fundamental.

use super::{
    DuffingBackboneError, arithmetic_geometric_mean, binary_parts, duffing_backbone_exact,
    scaled_parameters,
};
use fs_math::{dd::Dd, det, qd::Qd};

mod kinematics;
pub use kinematics::{DuffingKinematicSpectrum, duffing_kinematic_spectrum_exact};

/// Maximum number of retained odd harmonics; this is a work/storage limit.
pub const MAX_DUFFING_HARMONICS: usize = 4096;

/// A finite spectral image of the full conservative Duffing orbit.
#[derive(Debug, Clone, PartialEq)]
pub struct DuffingSpectrum {
    /// Fundamental angular frequency [rad/s], from the full energy integral.
    pub angular_frequency: f64,
    /// Signed peak cosine amplitudes in the original displacement units.
    /// Entry `j` belongs to order `2*j+1`; DC and every even order are zero.
    /// Thus `q(t) = sum_j a_j*cos((2*j+1)*omega*t)` plus the omitted tail.
    /// These are peak amplitudes: the positive-frequency complex Fourier
    /// coefficient used by harmonic balance is `a_j/2`, not `a_j`.
    pub cosine_amplitudes: Vec<f64>,
    /// Analytic upper bound on the omitted displacement's absolute value at
    /// every phase, evaluated numerically in the original displacement units.
    /// It bounds truncation only, excluding floating-point evaluation error;
    /// this is not a validated interval enclosure or a coupled-model bound.
    /// A positive tail smaller than the least subnormal is reported as that
    /// subnormal, never as exact zero. No retained coefficients are rescaled
    /// to hide the missing turning displacement.
    pub displacement_tail_bound: f64,
}

/// Exact-formula Fourier amplitudes at an authored physical turning point.
///
/// The same original equation, admissible softening well, binary scaling and
/// near-separatrix cancellation policy as [`duffing_backbone_exact`] apply.
/// `odd_harmonics` explicitly retains orders `1,3,...,2*odd_harmonics-1`;
/// it must be in `1..=MAX_DUFFING_HARMONICS`. Refining it changes neither the
/// frequency nor any already retained coefficient. The fundamental generally
/// differs from `amplitude`; this difference is physical, not normalization.
///
/// Two bounded AGMs evaluate the elliptic nome. A separately converged nome
/// sum normalizes the *infinite* series at the turning point, preserving the
/// linear limit without dividing two tiny elliptic prefactors. The returned
/// finite series is never normalized to its retained order. Work and storage
/// are linear in the retained count after this scalar preparation; there is
/// no time stepping, frequency-bin sampling, or nonlinear solve.
///
/// Useful as a physical harmonic-balance seed and independent conservative
/// reference. It adds no damping, forcing, mode coupling or stability claim,
/// and is not an instrument-render replacement. For velocity seeds, multiply
/// each complex displacement coefficient by `i*(2*j+1)*angular_frequency`.
/// [`duffing_kinematic_spectrum_exact`] performs derivative scaling before
/// rounding and reports the separate physical velocity/acceleration tails.
///
/// References: <https://dlmf.nist.gov/22.11.E2>,
/// <https://dlmf.nist.gov/22.11.E4>, <https://dlmf.nist.gov/22.2.E1>.
pub fn duffing_spectrum_exact(
    omega0: f64,
    beta: f64,
    amplitude: f64,
    odd_harmonics: usize,
) -> Result<DuffingSpectrum, DuffingBackboneError> {
    prepare_spectrum(omega0, beta, amplitude, odd_harmonics)?.displacement(odd_harmonics)
}

/// Scalar preparation shared by displacement and physical derivative queries.
struct PreparedSpectrum {
    angular_frequency: f64,
    amplitude: f64,
    first_ratio: f64,
    root_nome: f64,
    decay: f64,
    one_minus_nome: f64,
    softening: bool,
}

fn prepare_spectrum(
    omega0: f64,
    beta: f64,
    amplitude: f64,
    odd_harmonics: usize,
) -> Result<PreparedSpectrum, DuffingBackboneError> {
    if !(1..=MAX_DUFFING_HARMONICS).contains(&odd_harmonics) {
        return Err(DuffingBackboneError::InvalidInput);
    }
    let angular_frequency = duffing_backbone_exact(omega0, beta, amplitude)?;
    if beta == 0.0 || amplitude == 0.0 {
        return Ok(PreparedSpectrum {
            angular_frequency,
            amplitude,
            first_ratio: 1.0,
            root_nome: 0.0,
            decay: f64::INFINITY,
            one_minus_nome: 1.0,
            softening: false,
        });
    }
    let parameters = scaled_parameters(omega0, beta, amplitude)?;
    let softening = beta < 0.0;
    let (carrier, other) = if softening {
        (parameters.middle, parameters.endpoint)
    } else {
        (parameters.endpoint, parameters.middle)
    };
    // Compute k before squaring its small numerator. A dimensionless k^2
    // (or nome) may underflow while a physical third harmonic is representable.
    let (b, eb) = binary_parts(beta.abs());
    let (a, ea) = binary_parts(amplitude);
    let exponent = eb + 2 * ea - 2 * parameters.rate_exponent;
    let small_modulus = Qd::from_f64(
        det::sqrt(b) * a / det::sqrt(2.0 * carrier)
            * if exponent.rem_euclid(2) == 0 {
                1.0
            } else {
                core::f64::consts::SQRT_2
            },
    )
    .scale_power_of_two(exponent.div_euclid(2))
    .to_f64();
    let complement = det::sqrt(other / carrier);
    // Near the separatrix, separate rounded sqrt(beta) and amplitude products
    // can put k above one. The independently retained positive well margin
    // supplies k' instead; this form is well-conditioned for k near one.
    let modulus = if complement < core::f64::consts::FRAC_1_SQRT_2 {
        det::sqrt((1.0 - complement) * (1.0 + complement))
    } else {
        small_modulus
    };
    if !(modulus > 0.0 && modulus <= 1.0 && complement > 0.0 && complement <= 1.0) {
        return Err(DuffingBackboneError::NumericalRange);
    }
    let decay = core::f64::consts::PI * arithmetic_geometric_mean(1.0, complement)?
        / arithmetic_geometric_mean(1.0, modulus)?;
    let root_nome = det::exp(-0.5 * decay);
    if !(root_nome > 0.0 && root_nome < 1.0) {
        return Err(DuffingBackboneError::NumericalRange);
    }
    let nome = root_nome * root_nome;
    let one_minus_nome = -det::expm1(-decay);
    let normalization = infinite_normalization(nome, decay, softening, one_minus_nome)?;
    Ok(PreparedSpectrum {
        angular_frequency,
        amplitude,
        first_ratio: denominator(decay, 1, softening) / normalization,
        root_nome,
        decay,
        one_minus_nome,
        softening,
    })
}

impl PreparedSpectrum {
    fn displacement(&self, odd_harmonics: usize) -> Result<DuffingSpectrum, DuffingBackboneError> {
        if self.root_nome == 0.0 {
            let mut cosine_amplitudes = vec![0.0; odd_harmonics];
            cosine_amplitudes[0] = self.amplitude;
            return Ok(DuffingSpectrum {
                angular_frequency: self.angular_frequency,
                cosine_amplitudes,
                displacement_tail_bound: 0.0,
            });
        }
        let mut numerator = self.amplitude * self.first_ratio;
        let mut cosine_amplitudes = Vec::with_capacity(odd_harmonics);
        for j in 0..odd_harmonics {
            let magnitude = numerator / denominator(self.decay, 2 * j + 1, self.softening);
            if !magnitude.is_finite() {
                return Err(DuffingBackboneError::NumericalRange);
            }
            cosine_amplitudes.push(if self.softening && j % 2 != 0 {
                -magnitude
            } else {
                magnitude
            });
            // Split the nome multiplication so a small dimensionless nome does
            // not erase a representable physical harmonic before multiplication.
            numerator = (numerator * self.root_nome) * self.root_nome;
        }
        // If a_j = C q^j/(1 +/- q^(2j+1)), the omitted absolute sum is
        // <= C q^N/(1-q) for cn; for cd also divide by 1-q^(2N+1).
        let tail_denominator = if self.softening {
            denominator(self.decay, 2 * odd_harmonics + 1, true)
        } else {
            1.0
        };
        let tail = (numerator / self.one_minus_nome) / tail_denominator;
        let displacement_tail_bound = if tail == 0.0 {
            f64::from_bits(1)
        } else {
            tail.next_up()
        };
        if !displacement_tail_bound.is_finite() {
            return Err(DuffingBackboneError::NumericalRange);
        }
        Ok(DuffingSpectrum {
            angular_frequency: self.angular_frequency,
            cosine_amplitudes,
            displacement_tail_bound,
        })
    }
}

fn denominator(decay: f64, order: usize, softening: bool) -> f64 {
    let exponent = -(order as f64) * decay;
    if softening {
        -det::expm1(exponent)
    } else {
        1.0 + det::exp(exponent)
    }
}

/// Normalize ratios to the first coefficient, independently of retained N.
/// cd is an alternating decreasing series; cn has a geometric tail bound.
fn infinite_normalization(
    nome: f64,
    decay: f64,
    softening: bool,
    one_minus_nome: f64,
) -> Result<f64, DuffingBackboneError> {
    let first_denominator = denominator(decay, 1, softening);
    let mut power = nome;
    let mut sum = Dd::ONE;
    for j in 1..=MAX_DUFFING_HARMONICS {
        let magnitude = power * first_denominator / denominator(decay, 2 * j + 1, softening);
        let remainder = if softening {
            magnitude
        } else {
            power * first_denominator / one_minus_nome
        };
        if remainder <= 0.125 * f64::EPSILON * sum.to_f64().abs() {
            let normalization = sum.to_f64();
            return if normalization.is_finite() && normalization > 0.0 {
                Ok(normalization)
            } else {
                Err(DuffingBackboneError::NumericalRange)
            };
        }
        sum = sum
            + Dd::from_f64(if softening && j % 2 != 0 {
                -magnitude
            } else {
                magnitude
            });
        power *= nome;
    }
    Err(DuffingBackboneError::NumericalRange)
}

#[cfg(test)]
mod tests;
