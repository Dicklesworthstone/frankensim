//! Conservative Duffing frequency at a physical turning-point amplitude.
//!
//! Energy separation gives
//! `T = 4 integral_0^(pi/2) [omega0^2 + beta*A^2*(1+sin(theta)^2)/2]^(-1/2) dtheta`.
//! The complete elliptic integral identity (NIST DLMF 19.8.5) therefore gives
//! `omega = AGM(sqrt(omega0^2+beta*A^2), sqrt(omega0^2+beta*A^2/2))`.
//! This evaluates the original scalar conservative orbit, with no harmonic
//! truncation, time integration, quadrature, or small-amplitude expansion.

use fs_math::{det, qd::Qd};

/// Refusal of the scalar conservative periodic-orbit query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuffingBackboneError {
    /// Frequencies/amplitudes must be finite and nonnegative; beta finite.
    InvalidInput,
    /// No finite-period orbit about zero: no restoring force, or a softening
    /// amplitude at or beyond the separatrix `|beta|*A^2 = omega0^2`.
    NoPeriodicOrbit,
    /// An active stress channel drives a different displacement coordinate.
    NonInvariantAxis,
    /// The physical frequency/derivative cannot be represented, or bounded
    /// AGM refinement did not resolve its floating-point limit.
    NumericalRange,
}

impl core::fmt::Display for DuffingBackboneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-NLMODAL-DUFFING: {self:?}")
    }
}
impl std::error::Error for DuffingBackboneError {}

/// A direction in the ORIGINAL Duffing parameters, per caller-owned parameter.
/// These are partial input changes, not logarithmic derivatives. In particular,
/// a nonzero `beta` direction is meaningful at a linear baseline `beta = 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DuffingBackboneDirection {
    pub omega0: f64,
    pub beta: f64,
    pub amplitude: f64,
}

/// Conservative turning-amplitude frequency and its analytic local derivative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DuffingBackboneDerivative {
    pub angular_frequency_rad_s: f64,
    /// rad/s per unit of the parameter defining the supplied direction.
    pub angular_frequency_derivative: f64,
}

/// Differentiate the complete conservative period, including large-amplitude
/// hardening and the steep softening slope approaching the separatrix.
///
/// This differentiates the same energy-integral AGM as [`duffing_backbone_exact`],
/// with no differencing step, weak-amplitude expansion or numerical quadrature.
/// It does not infer a damped/forced peak or a coupled orbit's stability. Source
/// and direction units are unchanged. Frequency and derivative products are
/// independently binary-scaled before quad-double contraction; only the final
/// physical derivative must fit `f64`. An unrepresentable derivative refuses.
///
/// Each of the frequency and tangent passes is bounded by 32 AGM iterations
/// and uses constant memory. The frequency uses the original owner unchanged.
/// At zero amplitude the limit is differentiated from the admitted nonnegative
/// side; a free baseline refuses because it has no finite-period orbit.
pub fn duffing_backbone_directional_derivative(
    omega0: f64,
    beta: f64,
    amplitude: f64,
    direction: DuffingBackboneDirection,
) -> Result<DuffingBackboneDerivative, DuffingBackboneError> {
    if [direction.omega0, direction.beta, direction.amplitude]
        .iter()
        .any(|value| !value.is_finite())
    {
        return Err(DuffingBackboneError::InvalidInput);
    }
    let frequency = duffing_backbone_exact(omega0, beta, amplitude)?;
    if amplitude == 0.0 {
        return Ok(DuffingBackboneDerivative {
            angular_frequency_rad_s: frequency,
            angular_frequency_derivative: direction.omega0,
        });
    }
    let product = |factors: &[f64]| {
        let mut value = Qd::from_f64(1.0);
        let mut exponent = 0;
        for &factor in factors {
            if factor == 0.0 {
                return (Qd::ZERO, -8192);
            }
            let (mantissa, power) = binary_parts(factor.abs());
            value = value * mantissa.copysign(factor);
            exponent += power;
        }
        (value, exponent)
    };
    let linear = product(&[omega0, omega0]);
    let nonlinear = product(&[beta, amplitude, amplitude]);
    let rate_exponent = linear.1.max(nonlinear.1).div_euclid(2);
    let scaled = |term: (Qd, i32), exponent| term.0.scale_power_of_two(term.1 - exponent);
    let linear = scaled(linear, 2 * rate_exponent);
    let nonlinear = scaled(nonlinear, 2 * rate_exponent);
    let mut a = det::sqrt((linear + nonlinear).to_f64());
    let mut b = det::sqrt((linear + nonlinear * 0.5).to_f64());
    let linear_direction = product(&[2.0, omega0, direction.omega0]);
    let beta_direction = product(&[direction.beta, amplitude, amplitude]);
    let amplitude_direction = product(&[2.0, beta, amplitude, direction.amplitude]);
    let derivative_exponent = linear_direction
        .1
        .max(beta_direction.1)
        .max(amplitude_direction.1);
    let dl = scaled(linear_direction, derivative_exponent);
    let dn = scaled(beta_direction, derivative_exponent)
        + scaled(amplitude_direction, derivative_exponent);
    let mut da = (dl + dn) / (2.0 * a);
    let mut db = (dl + dn * 0.5) / (2.0 * b);
    for _ in 0..32 {
        let mean = 0.5 * (a + b);
        let mean_direction = (da + db) * 0.5;
        // At coincident primal iterates the two AGM partial derivatives are
        // exactly 1/2, even when the two incoming tangent values differ.
        if (a - b).abs() <= 4.0 * f64::EPSILON * mean {
            let derivative = mean_direction
                .scale_power_of_two(derivative_exponent - rate_exponent)
                .to_f64();
            if !derivative.is_finite() || (derivative == 0.0 && !mean_direction.is_zero()) {
                return Err(DuffingBackboneError::NumericalRange);
            }
            return Ok(DuffingBackboneDerivative {
                angular_frequency_rad_s: frequency,
                angular_frequency_derivative: derivative,
            });
        }
        let geometric = det::sqrt(a * b);
        let geometric_direction = (da * b + db * a) / (2.0 * geometric);
        a = mean;
        b = geometric;
        da = mean_direction;
        db = geometric_direction;
    }
    Err(DuffingBackboneError::NumericalRange)
}

/// Exact conservative angular frequency [rad/s] for
/// `q'' + omega0^2*q + beta*q^3 = 0` with `q(0)=amplitude`, `q'(0)=0`.
///
/// Hardening, pure quartic motion, and softening strictly inside the well are
/// supported. Zero amplitude returns the small-amplitude limit when `omega0`
/// is positive. A free coordinate or the softening separatrix refuses instead
/// of presenting a nonperiodic trajectory as a zero-frequency resonance.
///
/// `amplitude` is the physical turning displacement, **not** the fundamental
/// Fourier amplitude pinned by a harmonic-balance solver. The result is useful
/// as that solver's frequency seed or as an independent frequency reference;
/// it does not establish a forced-response peak, damping, or coupled-mode
/// stability. A coupled structure needs an invariant scalar displacement axis
/// before this frequency describes its actual motion.
///
/// Binary scaling avoids squared-frequency overflow/underflow. The existing
/// quad-double owner retains the cancellation of the original input products
/// near the softening separatrix; taking `sqrt(|beta|)*amplitude` first would
/// lose that separation. The final AGM uses deterministic square roots, at
/// most 32 iterations, constant memory, and no allocation. This is numerical
/// evaluation of an exact formula, not a validated interval enclosure.
///
/// Reference: <https://dlmf.nist.gov/19.8.E5> and the energy integral above.
pub fn duffing_backbone_exact(
    omega0: f64,
    beta: f64,
    amplitude: f64,
) -> Result<f64, DuffingBackboneError> {
    if !omega0.is_finite()
        || omega0 < 0.0
        || !beta.is_finite()
        || !amplitude.is_finite()
        || amplitude < 0.0
    {
        return Err(DuffingBackboneError::InvalidInput);
    }
    if beta == 0.0 || amplitude == 0.0 {
        return if omega0 > 0.0 {
            Ok(omega0)
        } else {
            Err(DuffingBackboneError::NoPeriodicOrbit)
        };
    }
    let (b, eb) = binary_parts(beta.abs());
    let (a, ea) = binary_parts(amplitude);
    let (w, ew) = if omega0 == 0.0 {
        (0.0, -4096)
    } else {
        binary_parts(omega0)
    };
    let nonlinear_exponent = eb + 2 * ea;
    // An even exponent makes the final frequency recovery a power of two.
    let rate_exponent = (2 * ew).max(nonlinear_exponent).div_euclid(2);
    let square_exponent = 2 * rate_exponent;
    let linear = (Qd::from_f64(w) * w).scale_power_of_two(2 * ew - square_exponent);
    let nonlinear = (Qd::from_f64(b.copysign(beta)) * a * a)
        .scale_power_of_two(nonlinear_exponent - square_exponent);
    let endpoint = (linear + nonlinear).to_f64();
    if endpoint <= 0.0 {
        return Err(DuffingBackboneError::NoPeriodicOrbit);
    }
    let middle = (linear + nonlinear * 0.5).to_f64();
    let endpoint_rate = det::sqrt(endpoint);
    let middle_rate = det::sqrt(middle);
    let mut high = endpoint_rate.max(middle_rate);
    let mut low = endpoint_rate.min(middle_rate);
    for _ in 0..32 {
        let mean = 0.5 * (high + low);
        if high - low <= 4.0 * f64::EPSILON * high {
            let frequency = Qd::from_f64(mean)
                .scale_power_of_two(rate_exponent)
                .to_f64();
            return if frequency.is_finite() && frequency > 0.0 {
                Ok(frequency)
            } else {
                Err(DuffingBackboneError::NumericalRange)
            };
        }
        low = det::sqrt(high * low);
        high = mean;
    }
    Err(DuffingBackboneError::NumericalRange)
}

/// Exact `x = mantissa * 2^exponent`, with `1 <= mantissa < 2`.
fn binary_parts(x: f64) -> (f64, i32) {
    let (value, correction) = if x < f64::MIN_POSITIVE {
        (x * 18_014_398_509_481_984.0, -54)
    } else {
        (x, 0)
    };
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
    (mantissa, exponent)
}

/// Source adapter's scalar reduction, after original storage admission. Each
/// active channel must preserve the axis; no cancellation or coupling cutoff
/// manufactures a decoupled model. The borrowed source arrays remain intact.
pub(crate) fn invariant_axis_backbone<'a>(
    omegas: &[f64],
    channels: impl Iterator<Item = (f64, &'a [f64])>,
    mode: usize,
    amplitude: f64,
) -> Result<f64, DuffingBackboneError> {
    let (omega, beta) = invariant_axis_parameters(omegas, channels, mode)?;
    duffing_backbone_exact(omega, beta, amplitude)
}

/// Source-bound amplitude derivative, with the same exact invariant-axis
/// admission as the frequency itself.
pub(crate) fn invariant_axis_backbone_amplitude_derivative<'a>(
    omegas: &[f64],
    channels: impl Iterator<Item = (f64, &'a [f64])>,
    mode: usize,
    amplitude: f64,
) -> Result<DuffingBackboneDerivative, DuffingBackboneError> {
    let (omega, beta) = invariant_axis_parameters(omegas, channels, mode)?;
    duffing_backbone_directional_derivative(
        omega,
        beta,
        amplitude,
        DuffingBackboneDirection {
            omega0: 0.0,
            beta: 0.0,
            amplitude: 1.0,
        },
    )
}

fn invariant_axis_parameters<'a>(
    omegas: &[f64],
    channels: impl Iterator<Item = (f64, &'a [f64])>,
    mode: usize,
) -> Result<(f64, f64), DuffingBackboneError> {
    let n = omegas.len();
    let omega = *omegas.get(mode).ok_or(DuffingBackboneError::InvalidInput)?;
    let square = n.checked_mul(n).ok_or(DuffingBackboneError::InvalidInput)?;
    let mut beta = 0.0;
    for (coefficient, coupling) in channels {
        if coupling.len() != square || !coefficient.is_finite() || coefficient < 0.0 {
            return Err(DuffingBackboneError::InvalidInput);
        }
        let diagonal = coupling[mode * n + mode];
        if coefficient != 0.0 && diagonal != 0.0 {
            if (0..n).any(|j| {
                j != mode && (coupling[mode * n + j] != 0.0 || coupling[j * n + mode] != 0.0)
            }) {
                return Err(DuffingBackboneError::NonInvariantAxis);
            }
            let term = coefficient * diagonal * diagonal;
            if !term.is_finite() || term <= 0.0 {
                return Err(DuffingBackboneError::NumericalRange);
            }
            beta += term;
        }
    }
    Ok((omega, beta))
}

#[cfg(test)]
mod tests {
    use super::{
        DuffingBackboneDirection, DuffingBackboneError, duffing_backbone_directional_derivative,
        duffing_backbone_exact, invariant_axis_backbone,
        invariant_axis_backbone_amplitude_derivative,
    };

    fn close(got: f64, expected: f64, relative: f64) {
        assert!(
            (got - expected).abs() <= relative * expected.abs(),
            "frequency {got:.17e}, expected {expected:.17e}"
        );
    }

    // Independent composite Simpson integration of the energy-separated
    // quarter period. This reference uses no elliptic integral or AGM.
    fn period_frequency(omega0: f64, beta: f64, amplitude: f64) -> f64 {
        let count = 16_384;
        let step = core::f64::consts::FRAC_PI_2 / f64::from(count);
        let mut sum = 0.0;
        for i in 0..=count {
            let sine = (f64::from(i) * step).sin();
            let weight = if i == 0 || i == count {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight
                / (omega0 * omega0 + 0.5 * beta * amplitude * amplitude * (1.0 + sine * sine))
                    .sqrt();
        }
        core::f64::consts::TAU / (4.0 * step * sum / 3.0)
    }

    #[test]
    fn g1_full_amplitude_backbone_matches_independent_energy_periods() {
        for ratio in [-0.999, -0.9, -0.1, -1e-4, 0.0, 1e-4, 0.1, 1.0, 10.0, 1e6] {
            let (omega, amplitude) = (7.0, 0.4);
            let beta = ratio * omega * omega / (amplitude * amplitude);
            close(
                duffing_backbone_exact(omega, beta, amplitude).unwrap(),
                period_frequency(omega, beta, amplitude),
                2e-12,
            );
        }
    }

    #[test]
    fn g1_linear_weakly_nonlinear_and_pure_quartic_limits() {
        for omega in [1e-300, 1.0, 1e300] {
            assert_eq!(duffing_backbone_exact(omega, 0.0, 1e200).unwrap(), omega);
            assert_eq!(duffing_backbone_exact(omega, -1e300, 0.0).unwrap(), omega);
        }
        for ratio in [-1e-3, 1e-3] {
            let expansion = 1.0 + 3.0 * ratio / 8.0 - 21.0 * ratio * ratio / 256.0;
            close(
                duffing_backbone_exact(1.0, ratio, 1.0).unwrap(),
                expansion,
                4e-11,
            );
        }
        // pi/(2 K(1/sqrt(2))), independently tabulated lemniscatic constant.
        close(
            duffing_backbone_exact(0.0, 1.0, 1.0).unwrap(),
            0.847_213_084_793_979_1,
            8e-16,
        );
    }

    #[test]
    fn g1_near_separatrix_keeps_cancellation_before_the_elliptic_limit() {
        let t = 2.0_f64.powi(-53);
        let amplitude = 1.0 - t;
        let beta = -(1.0 + 2.0 * t);
        // 1 - (1+2t)*(1-t)^2 = 3t^2 - 2t^3. Binary64 multiplication
        // alone erases the positive well margin and invents a separatrix.
        assert_eq!(1.0 + beta * (amplitude * amplitude), 0.0);
        let margin = (3.0 - 2.0 * t) * t * t;
        let high = ((1.0 + margin) * 0.5).sqrt();
        let low = margin.sqrt();
        // K(k) = log(4/k') + O(k'^2 log(1/k')); here k'^2 < 1e-30.
        let asymptotic = core::f64::consts::FRAC_PI_2 * high / (4.0 * high / low).ln();
        close(
            duffing_backbone_exact(1.0, beta, amplitude).unwrap(),
            asymptotic,
            2e-15,
        );
        let closer = duffing_backbone_exact(1.0, -1.0, 1.0_f64.next_down()).unwrap();
        assert!(closer < duffing_backbone_exact(1.0, -1.0, 0.999).unwrap());
        assert!(duffing_backbone_exact(1.0, beta, amplitude).unwrap() < closer);
    }

    #[test]
    fn g3_time_and_displacement_units_survive_squared_frequency_overflow_and_underflow() {
        for ratio in [-0.75, 3.0] {
            let reference = duffing_backbone_exact(1.0, ratio, 1.0).unwrap();
            for (time_exp, amplitude_exp) in [(400, 0), (-400, 0), (800, 400), (-800, -400)] {
                let time = 2.0_f64.powi(time_exp);
                let amplitude = 2.0_f64.powi(amplitude_exp);
                let beta = ratio * 2.0_f64.powi(2 * (time_exp - amplitude_exp));
                close(
                    duffing_backbone_exact(time, beta, amplitude).unwrap() / time,
                    reference,
                    2e-15,
                );
            }
        }
        // sqrt(beta)*A itself overflows, but the physical frequency does not.
        close(
            duffing_backbone_exact(0.0, 4.0, 1e308).unwrap(),
            1.694_426_169_587_958_2e308,
            8e-16,
        );
        let smallest = f64::from_bits(1);
        let rate = 2.0_f64.powi(163);
        close(
            duffing_backbone_exact(rate, smallest, 2.0_f64.powi(700)).unwrap() / rate,
            duffing_backbone_exact(1.0, 1.0, 1.0).unwrap(),
            2e-15,
        );
    }

    #[test]
    fn g0_nonperiodic_and_invalid_sources_never_return_a_fabricated_resonance() {
        for (omega, beta, amplitude) in [
            (0.0, 0.0, 1.0),
            (0.0, 1.0, 0.0),
            (0.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, -1.0, 1.01),
        ] {
            assert_eq!(
                duffing_backbone_exact(omega, beta, amplitude),
                Err(DuffingBackboneError::NoPeriodicOrbit)
            );
        }
        for (omega, beta, amplitude) in [
            (-1.0, 1.0, 1.0),
            (1.0, 1.0, -1.0),
            (f64::NAN, 1.0, 1.0),
            (1.0, f64::INFINITY, 1.0),
            (1.0, 1.0, f64::NAN),
        ] {
            assert_eq!(
                duffing_backbone_exact(omega, beta, amplitude),
                Err(DuffingBackboneError::InvalidInput)
            );
        }
        assert_eq!(
            duffing_backbone_exact(f64::MAX, f64::MAX, f64::MAX),
            Err(DuffingBackboneError::NumericalRange)
        );
    }

    #[test]
    fn g1_source_axis_uses_all_original_stress_coefficients() {
        let omegas = [7.0, 11.0];
        let first = [2.0, 0.0, 0.0, -3.0];
        let second = [-4.0, 0.0, 0.0, 5.0];
        for mode in 0..2 {
            let beta =
                2.0 * first[3 * mode] * first[3 * mode] + 3.0 * second[3 * mode] * second[3 * mode];
            assert_eq!(
                invariant_axis_backbone(
                    &omegas,
                    [(2.0, &first[..]), (3.0, &second[..])].into_iter(),
                    mode,
                    0.7
                )
                .unwrap(),
                duffing_backbone_exact(omegas[mode], beta, 0.7).unwrap(),
            );
            let response = invariant_axis_backbone_amplitude_derivative(
                &omegas,
                [(2.0, &first[..]), (3.0, &second[..])].into_iter(),
                mode,
                0.7,
            )
            .unwrap();
            let direction = DuffingBackboneDirection {
                omega0: 0.0,
                beta: 0.0,
                amplitude: 1.0,
            };
            close(
                response.angular_frequency_derivative,
                integral_derivative(omegas[mode], beta, 0.7, direction),
                3e-12,
            );
        }
    }

    #[test]
    fn g0_source_axis_never_discards_a_coupled_force_or_invents_an_axis() {
        let omegas = [7.0, 11.0];
        let coupled = [2.0, 0.5, 0.5, -3.0];
        let zero_along_axis = [0.0, 0.5, 0.5, -3.0];
        for mode in 0..2 {
            assert_eq!(
                invariant_axis_backbone(&omegas, [(1.0, &coupled[..])].into_iter(), mode, 0.7),
                Err(DuffingBackboneError::NonInvariantAxis),
            );
            assert_eq!(
                invariant_axis_backbone_amplitude_derivative(
                    &omegas,
                    [(1.0, &coupled[..])].into_iter(),
                    mode,
                    0.7,
                ),
                Err(DuffingBackboneError::NonInvariantAxis),
            );
            assert_eq!(
                invariant_axis_backbone(&omegas, [(0.0, &coupled[..])].into_iter(), mode, 0.7)
                    .unwrap(),
                omegas[mode],
            );
        }
        // A channel with q^T E q identically zero on the selected axis has
        // zero restoring force there, even with nonzero transverse tangents.
        assert_eq!(
            invariant_axis_backbone(&omegas, [(1.0, &zero_along_axis[..])].into_iter(), 0, 0.7)
                .unwrap(),
            omegas[0]
        );
        assert_eq!(
            invariant_axis_backbone(&omegas, [].into_iter(), 2, 0.7),
            Err(DuffingBackboneError::InvalidInput)
        );
        assert_eq!(
            invariant_axis_backbone(&omegas, [(1.0, &coupled[..3])].into_iter(), 0, 0.7),
            Err(DuffingBackboneError::InvalidInput)
        );
        // A positive cubic coefficient that underflows must not silently turn
        // the source into a linear spring; large amplitude can reveal it.
        assert_eq!(
            invariant_axis_backbone(&[1.0], [(1e-200, &[1e-200][..])].into_iter(), 0, 1e300),
            Err(DuffingBackboneError::NumericalRange),
        );
    }

    // Independently differentiate the energy-separated quarter-period integral.
    // This uses ordinary real quadrature and no differentiated AGM recurrence.
    fn integral_derivative(w: f64, beta: f64, a: f64, d: DuffingBackboneDirection) -> f64 {
        let count = 16_384;
        let step = core::f64::consts::FRAC_PI_2 / f64::from(count);
        let mut period = 0.0;
        let mut derivative = 0.0;
        for i in 0..=count {
            let sine = (f64::from(i) * step).sin();
            let shape = 0.5 * (1.0 + sine * sine);
            let stiffness = w * w + beta * a * a * shape;
            let tangent =
                2.0 * w * d.omega0 + (d.beta * a * a + 2.0 * beta * a * d.amplitude) * shape;
            let weight = if i == 0 || i == count {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            period += weight / stiffness.sqrt();
            derivative -= 0.5 * weight * tangent / (stiffness * stiffness.sqrt());
        }
        -core::f64::consts::TAU * derivative / (4.0 * step * period * period / 3.0)
    }

    #[test]
    fn g1_exact_backbone_tuning_matches_independent_period_integral() {
        for ratio in [-0.99, -0.7, -0.1, 0.0, 0.01, 1.0, 100.0, 1e6] {
            let (w, a) = (7.0, 0.4);
            let beta = ratio * w * w / (a * a);
            for d in [
                DuffingBackboneDirection {
                    omega0: 1.0,
                    beta: 0.0,
                    amplitude: 0.0,
                },
                DuffingBackboneDirection {
                    omega0: 0.0,
                    beta: 1.0,
                    amplitude: 0.0,
                },
                DuffingBackboneDirection {
                    omega0: 0.0,
                    beta: 0.0,
                    amplitude: 1.0,
                },
                DuffingBackboneDirection {
                    omega0: 0.3,
                    beta: -0.2,
                    amplitude: 0.1,
                },
            ] {
                let actual = duffing_backbone_directional_derivative(w, beta, a, d).unwrap();
                assert_eq!(
                    actual.angular_frequency_rad_s,
                    duffing_backbone_exact(w, beta, a).unwrap()
                );
                let expected = integral_derivative(w, beta, a, d);
                if expected == 0.0 {
                    assert_eq!(actual.angular_frequency_derivative, 0.0);
                } else {
                    close(actual.angular_frequency_derivative, expected, 3e-12);
                }
            }
        }
    }

    #[test]
    fn g1_linear_pure_quartic_and_zero_amplitude_tuning_limits() {
        let d = DuffingBackboneDirection {
            omega0: 0.0,
            beta: 1.0,
            amplitude: 0.0,
        };
        close(
            duffing_backbone_directional_derivative(7.0, 0.0, 0.4, d)
                .unwrap()
                .angular_frequency_derivative,
            3.0 * 0.4 * 0.4 / (8.0 * 7.0),
            2e-15,
        );
        let d = DuffingBackboneDirection {
            omega0: 0.2,
            beta: -2.0,
            amplitude: 3.0,
        };
        assert_eq!(
            duffing_backbone_directional_derivative(7.0, -30.0, 0.0, d)
                .unwrap()
                .angular_frequency_derivative,
            0.2,
        );
        let d = DuffingBackboneDirection {
            omega0: 0.0,
            beta: 0.0,
            amplitude: 1.0,
        };
        let response = duffing_backbone_directional_derivative(0.0, 20.0, 0.4, d).unwrap();
        close(
            response.angular_frequency_derivative,
            response.angular_frequency_rad_s / 0.4,
            2e-15,
        );
    }

    #[test]
    fn g3_tuning_retains_homogeneity_and_cancellation_across_units() {
        for (rate_exp, amplitude_exp) in [(0, 0), (800, 400), (-800, -400)] {
            let w = 2.0_f64.powi(rate_exp);
            let a = 2.0_f64.powi(amplitude_exp);
            let beta = 2.0_f64.powi(2 * (rate_exp - amplitude_exp));
            let d = DuffingBackboneDirection {
                omega0: w,
                beta: 2.0 * beta,
                amplitude: 0.0,
            };
            let response = duffing_backbone_directional_derivative(w, beta, a, d).unwrap();
            close(
                response.angular_frequency_derivative,
                response.angular_frequency_rad_s,
                3e-15,
            );
            let d = DuffingBackboneDirection {
                omega0: 0.0,
                beta: -2.0 * beta,
                amplitude: a,
            };
            assert_eq!(
                duffing_backbone_directional_derivative(w, beta, a, d)
                    .unwrap()
                    .angular_frequency_derivative,
                0.0
            );
        }
    }

    #[test]
    fn g1_softening_tuning_diverges_toward_the_physical_separatrix() {
        let t = 2.0_f64.powi(-53);
        let a = 1.0 - t;
        let d = DuffingBackboneDirection {
            omega0: 0.0,
            beta: 0.0,
            amplitude: 1.0,
        };
        for (beta, margin) in [
            (-1.0, (1.0 - a) * (1.0 + a)),
            // Ordinary binary64 stiffness products erase this well margin.
            (-(1.0 + 2.0 * t), (3.0 - 2.0 * t) * t * t),
        ] {
            let response = duffing_backbone_directional_derivative(1.0, beta, a, d).unwrap();
            let high = ((1.0 + margin) * 0.5).sqrt();
            let logarithm = (4.0 * high / margin.sqrt()).ln();
            let log_high_derivative = beta * a / (1.0 + margin);
            let expected = response.angular_frequency_rad_s
                * (log_high_derivative - (log_high_derivative - beta * a / margin) / logarithm);
            close(response.angular_frequency_derivative, expected, 2e-14);
            assert!(response.angular_frequency_derivative < -1e12);
        }
    }

    #[test]
    fn g0_tuning_refuses_invalid_or_unrepresentable_results() {
        let bad = DuffingBackboneDirection {
            omega0: 0.0,
            beta: f64::NAN,
            amplitude: 0.0,
        };
        assert_eq!(
            duffing_backbone_directional_derivative(1.0, 1.0, 1.0, bad),
            Err(DuffingBackboneError::InvalidInput)
        );
        let d = DuffingBackboneDirection {
            omega0: 0.0,
            beta: 1.0,
            amplitude: 0.0,
        };
        assert_eq!(
            duffing_backbone_directional_derivative(1.0, -1.0, 1.0, d),
            Err(DuffingBackboneError::NoPeriodicOrbit)
        );
        assert_eq!(
            duffing_backbone_directional_derivative(1e-300, 0.0, 1e300, d),
            Err(DuffingBackboneError::NumericalRange)
        );
    }
}
