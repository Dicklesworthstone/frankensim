//! Factored Maxwell spectra and frequency slopes, before physical scaling.
//!
//! Keeping the branch stiffness outside these binary-scaled coefficients avoids
//! losing a representable response when omega*tau or its square cannot fit.

#[derive(Clone, Copy)]
struct Scaled {
    fraction: f64,
    exponent: i32,
}

fn split(value: f64) -> (f64, i32) {
    let (value, correction) = if value.abs() < f64::MIN_POSITIVE {
        (value * 18_014_398_509_481_984.0, -54)
    } else {
        (value, 0)
    };
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
    let fraction = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (fraction, exponent)
}

impl Scaled {
    fn ratio(numerator: &[f64], denominator: &[f64]) -> Self {
        if numerator.contains(&0.0) {
            return Self {
                fraction: 0.0,
                exponent: 0,
            };
        }
        let (mut fraction, mut exponent) = (1.0, 0);
        for &value in numerator {
            let (f, e) = split(value);
            fraction *= f;
            exponent += e;
        }
        for &value in denominator {
            let (f, e) = split(value);
            fraction /= f;
            exponent -= e;
        }
        let (fraction, adjustment) = split(fraction);
        Self {
            fraction,
            exponent: exponent + adjustment,
        }
    }

    fn apply(self, stiffness: f64) -> Option<f64> {
        if !stiffness.is_finite() {
            return None;
        }
        if self.fraction == 0.0 || stiffness == 0.0 {
            return Some(0.0);
        }
        let (c, ce) = split(stiffness);
        let (fraction, adjustment) = split(self.fraction * c);
        let exponent = self.exponent + ce + adjustment;
        let result = if exponent > 1023 {
            return None;
        } else if exponent < -1075 {
            0.0_f64.copysign(fraction)
        } else if exponent < -1022 {
            // The first multiplication is exact normal binary scaling; only
            // the final multiplication rounds into the subnormal range.
            let scale = f64::from_bits(((exponent + 1074 + 1023) as u64) << 52);
            (fraction * scale) * f64::from_bits(1)
        } else {
            fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
        };
        result.is_finite().then_some(result)
    }
}

/// One original branch, including the derivative of C''/omega directly.
/// Primal and tangent representability are admitted independently.
pub(crate) struct PronySpectrumWeights {
    storage: Scaled,
    loss: Scaled,
    storage_derivative: Scaled,
    loss_derivative: Scaled,
    viscosity_derivative: Scaled,
}

impl PronySpectrumWeights {
    pub(crate) fn new(omega: f64, tau: f64) -> Option<Self> {
        if !omega.is_finite() || !tau.is_finite() || tau <= 0.0 {
            return None;
        }
        let w = omega.abs();
        let x = w * tau;
        if x <= 1.0 {
            let den = x.mul_add(x, 1.0);
            let error = w.mul_add(tau, -x);
            let difference = (-x).mul_add(x, 1.0) - 2.0 * x * error - error * error;
            Some(Self {
                storage: Scaled::ratio(&[w, tau, w, tau], &[den]),
                loss: Scaled::ratio(&[omega, tau], &[den]),
                storage_derivative: Scaled::ratio(&[2.0, omega, tau, tau], &[den, den]),
                loss_derivative: Scaled::ratio(&[tau, difference], &[den, den]),
                viscosity_derivative: Scaled::ratio(&[-2.0, omega, tau, tau, tau], &[den, den]),
            })
        } else {
            let u = 1.0 / x;
            let den = u.mul_add(u, 1.0);
            let difference = if x <= 2.0 {
                let error = w.mul_add(tau, -x);
                ((-x).mul_add(x, 1.0) - 2.0 * x * error - error * error) / (x * x)
            } else {
                u.mul_add(u, -1.0)
            };
            let sign = omega.signum();
            Some(Self {
                storage: Scaled::ratio(&[1.0], &[den]),
                loss: Scaled::ratio(&[sign], &[w, tau, den]),
                storage_derivative: Scaled::ratio(&[2.0 * sign], &[w, w, w, tau, tau, den, den]),
                loss_derivative: Scaled::ratio(&[difference], &[w, w, tau, den, den]),
                viscosity_derivative: Scaled::ratio(&[-2.0 * sign], &[w, w, w, tau, den, den]),
            })
        }
    }

    pub(crate) fn primal(&self, stiffness: f64) -> Option<[f64; 2]> {
        Some([self.storage.apply(stiffness)?, self.loss.apply(stiffness)?])
    }

    pub(crate) fn modulus_derivative(&self, stiffness: f64) -> Option<[f64; 2]> {
        Some([
            self.storage_derivative.apply(stiffness)?,
            self.loss_derivative.apply(stiffness)?,
        ])
    }

    /// [d storage/d omega, d loss/d omega, d(loss/omega)/d omega].
    pub(crate) fn derivative(&self, stiffness: f64) -> Option<[f64; 3]> {
        Some([
            self.storage_derivative.apply(stiffness)?,
            self.loss_derivative.apply(stiffness)?,
            self.viscosity_derivative.apply(stiffness)?,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::PronySpectrumWeights;

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 3e-14 * b.abs().max(1e-310),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn analytic_branch_slopes_and_frequency_parity() {
        for w in [0.0_f64, 0.125, 0.5, 1.0, 2.0, 8.0, 64.0] {
            let (tau, c) = (2.0, 3.0);
            let x = w * tau;
            let den = 1.0 + x * x;
            let expected = [
                2.0 * c * tau * x / den.powi(2),
                c * tau * (1.0 - x * x) / den.powi(2),
                -2.0 * c * tau * tau * x / den.powi(2),
            ];
            for sign in [-1.0, 1.0] {
                let weights = PronySpectrumWeights::new(sign * w, tau).unwrap();
                let primal = weights.primal(c).unwrap();
                close(primal[0], c * x * x / den);
                close(primal[1], sign * c * x / den);
                for (i, value) in weights.derivative(c).unwrap().into_iter().enumerate() {
                    close(value, expected[i] * if i == 1 { 1.0 } else { sign });
                }
            }
        }
    }

    #[test]
    fn overflowing_rate_time_keeps_representable_loss_and_slope() {
        let weights = PronySpectrumWeights::new(1e100, 1e250).unwrap();
        let primal = weights.primal(1e150).unwrap();
        close(primal[0], 1e150);
        close(primal[1], 1e-200);
        close(weights.derivative(1e150).unwrap()[1], -1e-300);
    }

    #[test]
    fn tiny_rate_time_retains_physical_storage_and_derivatives() {
        let weights = PronySpectrumWeights::new(1e-100, 1e-100).unwrap();
        let primal = weights.primal(1e150).unwrap();
        close(primal[0], 1e-250);
        close(primal[1], 1e-50);
        let derivative = weights.derivative(1e150).unwrap();
        close(derivative[0], 2e-150);
        close(derivative[1], 1e50);
        close(derivative[2], -2e-250);
    }

    #[test]
    fn subnormal_relaxation_time_is_not_discarded() {
        let weights = PronySpectrumWeights::new(2.0_f64.powi(500), f64::from_bits(1)).unwrap();
        let c = 2.0_f64.powi(500);
        close(weights.primal(c).unwrap()[1], 2.0_f64.powi(-74));
        close(weights.derivative(c).unwrap()[1], 2.0_f64.powi(-574));
    }

    #[test]
    fn finite_dc_primal_survives_an_unrepresentable_tangent() {
        let weights = PronySpectrumWeights::new(0.0, 1e200).unwrap();
        assert_eq!(weights.primal(1e150), Some([0.0, 0.0]));
        assert!(weights.derivative(1e150).is_none());
        let crossover = PronySpectrumWeights::new(1e-200, 1e200).unwrap();
        assert!(crossover.modulus_derivative(1.0).is_some());
        assert!(crossover.derivative(1.0).is_none());
    }

    #[test]
    fn invalid_frequency_time_or_stiffness_refuses() {
        for (w, tau) in [
            (f64::NAN, 1.0),
            (f64::INFINITY, 1.0),
            (1.0, 0.0),
            (1.0, -1.0),
        ] {
            assert!(PronySpectrumWeights::new(w, tau).is_none());
        }
        let weights = PronySpectrumWeights::new(1.0, 1.0).unwrap();
        assert!(weights.primal(f64::NAN).is_none());
        assert!(weights.derivative(f64::INFINITY).is_none());
    }
}
