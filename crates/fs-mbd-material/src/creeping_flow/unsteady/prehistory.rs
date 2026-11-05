//! Explicit preparation of the original translational Basset modes.

/// Physical slip history preceding the initial native body state at time zero.
///
/// The endpoint slip is `(p - m*U)/m`, in world axes. Geometry and material
/// transport are held fixed during preparation. The initial body state is
/// supplied, not integrated backwards; its holding forces and past work are
/// outside the subsequent simulation's accounting window. Neither option
/// supplies rotational prehistory or changes the native body momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SlipPrehistory {
    /// Comoving before `-duration_s`, then constant relative acceleration to
    /// the current slip at time zero. Duration must be finite and positive.
    /// Each original branch retains its exact response to that linear ramp;
    /// the preparation has finite continuum disturbance storage.
    LinearRamp { duration_s: f64 },
    /// The current slip has been maintained throughout the infinite past.
    /// Every finite spectral branch starts at that slip, with zero history
    /// force and nonzero stored energy. Absolute branch storage depends on the
    /// slow cutoff and diverges in the unbounded steady Stokes continuum;
    /// this is a finite spectral reference, not finite total exterior energy.
    SteadySlip,
}

impl SlipPrehistory {
    /// u(0)/w(0) for u_dot = rate*(w-u), with the declared preparation.
    /// None denotes an invalid or unrepresentable preparation timescale.
    pub(super) fn memory_fraction(self, rate_per_s: f64) -> Option<f64> {
        if !rate_per_s.is_finite() || rate_per_s <= 0.0 {
            return None;
        }
        let Self::LinearRamp { duration_s } = self else {
            return Some(1.0);
        };
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return None;
        }
        let alpha = rate_per_s * duration_s;
        if !alpha.is_finite() || alpha == 0.0 {
            return None;
        }
        // Exact convolution gives 1 + expm1(-alpha)/alpha. Near zero its
        // subtraction loses the retained acceleration history. The convergent
        // series starts at alpha/2; at alpha < 1/2 the first omitted term is
        // below (1/2)^18/19!, well below a binary64 rounding unit.
        let fraction = if alpha < 0.5 {
            let mut term = 0.5 * alpha;
            let mut sum = term;
            for divisor in 3..=18 {
                term *= -alpha / f64::from(divisor);
                sum += term;
            }
            sum
        } else {
            1.0 + fs_math::det::expm1(-alpha) / alpha
        };
        (fraction.is_finite() && fraction > 0.0 && fraction <= 1.0).then_some(fraction)
    }
}

#[cfg(test)]
mod tests {
    use super::SlipPrehistory;

    #[test]
    fn ramp_memory_matches_an_independent_convolution_integral() {
        for alpha in [1e-12, 1e-6, 0.1, 0.499, 0.5, 2.0, 10.0, 1000.0, 1e8] {
            // Integrate rate*exp(-rate*(T-s))*s/T, not the closed form.
            // The scaled variable resolves fast modes without a narrow layer.
            let upper = f64::min(alpha, 40.0);
            let n = 16384;
            let dx = upper / f64::from(n);
            let integrand = |y: f64| (1.0 - y / alpha) * (-y).exp();
            let mut sum = integrand(0.0) + integrand(upper);
            for i in 1..n {
                sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(f64::from(i) * dx);
            }
            let expected = sum * dx / 3.0;
            let actual = SlipPrehistory::LinearRamp { duration_s: 1.0 }
                .memory_fraction(alpha)
                .unwrap();
            assert!((actual / expected - 1.0).abs() < 3e-12, "alpha={alpha:e}");
        }
    }

    #[test]
    fn slow_modes_retain_tiny_acceleration_history_and_fast_modes_equilibrate() {
        for alpha in [1e-300, 1e-100, 1e-20, 1e-16] {
            let actual = SlipPrehistory::LinearRamp { duration_s: 1.0 }
                .memory_fraction(alpha)
                .unwrap();
            assert_eq!(actual, 0.5 * alpha);
            assert!(actual > 0.0);
        }
        let mut previous = 0.0;
        for alpha in [1e-10, 0.01, 0.49, 0.5, 1.0, 100.0, 1e10, 1e100] {
            let actual = SlipPrehistory::LinearRamp { duration_s: 1.0 }
                .memory_fraction(alpha)
                .unwrap();
            assert!(actual > previous && actual <= 1.0);
            previous = actual;
        }
        assert_eq!(previous, 1.0);
        assert_eq!(
            SlipPrehistory::SteadySlip.memory_fraction(1e-200),
            Some(1.0)
        );
    }

    #[test]
    fn invalid_or_unrepresentable_preparations_refuse() {
        for duration_s in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(
                SlipPrehistory::LinearRamp { duration_s }
                    .memory_fraction(1.0)
                    .is_none()
            );
        }
        for rate in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(SlipPrehistory::SteadySlip.memory_fraction(rate).is_none());
        }
        for (duration_s, rate) in [(1e300, 1e300), (1e-300, 1e-300)] {
            assert!(
                SlipPrehistory::LinearRamp { duration_s }
                    .memory_fraction(rate)
                    .is_none()
            );
        }
        let least = f64::from_bits(1);
        let ramp = SlipPrehistory::LinearRamp { duration_s: 1.0 };
        assert_eq!(ramp.memory_fraction(least), None);
        assert_eq!(ramp.memory_fraction(2.0 * least), Some(least));
    }
}
