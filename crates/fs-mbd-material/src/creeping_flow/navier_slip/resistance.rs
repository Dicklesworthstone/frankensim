//! Scalar closed-form resistance used by the real Navier spherical boundary.
//!
//! Ichiki, Kobryn & Kovalenko, arXiv:1302.0461, Eqs. (44), (49), (54):
//! Gamma_m,n = (1 + m*b/a)/(1 + n*b/a). Translation uses Gamma_2,3;
//! rotation uses Gamma_0,3. See <https://arxiv.org/abs/1302.0461>.
//! Both factors multiply the respective no-slip Stokes resistance.

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ResistanceFactors {
    pub translation: f64,
    pub rotation: f64,
}

impl ResistanceFactors {
    /// Finite physical b >= 0 and radius a > 0. Taking the smaller ratio keeps
    /// b/a, 3*b and a+3*b from overflowing. A finite-slip rotary resistance
    /// must remain representably positive; zero is not invented as shear-free.
    pub fn new(radius_m: f64, slip_length_m: f64) -> Option<Self> {
        if !radius_m.is_finite()
            || radius_m <= 0.0
            || !slip_length_m.is_finite()
            || slip_length_m < 0.0
        {
            return None;
        }
        let result = if slip_length_m <= radius_m {
            let slip_ratio = slip_length_m / radius_m;
            let denominator = 1.0 + 3.0 * slip_ratio;
            Self {
                translation: (1.0 + 2.0 * slip_ratio) / denominator,
                rotation: 1.0 / denominator,
            }
        } else {
            let radius_ratio = radius_m / slip_length_m;
            let denominator = radius_ratio + 3.0;
            Self {
                translation: (radius_ratio + 2.0) / denominator,
                rotation: radius_ratio / denominator,
            }
        };
        (result.rotation > 0.0).then_some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs(),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn zero_slip_is_exactly_the_original_stokes_boundary() {
        for radius in [f64::MIN_POSITIVE, 0.01, 1.0, f64::MAX] {
            assert_eq!(
                ResistanceFactors::new(radius, 0.0),
                Some(ResistanceFactors {
                    translation: 1.0,
                    rotation: 1.0
                })
            );
        }
    }

    #[test]
    fn finite_slip_changes_translation_and_spin_by_different_exact_factors() {
        let same = ResistanceFactors::new(0.01, 0.01).unwrap();
        assert_eq!(same.translation, 0.75);
        assert_eq!(same.rotation, 0.25);
        let third = ResistanceFactors::new(3.0, 1.0).unwrap();
        close(third.translation, 5.0 / 6.0);
        assert_eq!(third.rotation, 0.5);
        let twice = ResistanceFactors::new(1.0, 2.0).unwrap();
        close(twice.translation, 5.0 / 7.0);
        close(twice.rotation, 1.0 / 7.0);
    }

    #[test]
    fn large_slip_approaches_shear_free_translation_and_vanishing_torque() {
        let factors = ResistanceFactors::new(1.0, 1e150).unwrap();
        assert_eq!(factors.translation, 2.0 / 3.0);
        close(factors.rotation, 1e-150 / 3.0);
    }

    #[test]
    fn scaled_evaluation_avoids_overflow_and_preserves_the_physical_ratio() {
        for scale in [f64::MIN_POSITIVE, 1.0, f64::MAX / 4.0] {
            let factors = ResistanceFactors::new(scale, 2.0 * scale).unwrap();
            close(factors.translation, 5.0 / 7.0);
            close(factors.rotation, 1.0 / 7.0);
        }
        let factors = ResistanceFactors::new(0.01, f64::MAX).unwrap();
        assert_eq!(factors.translation, 2.0 / 3.0);
        assert!(factors.rotation > 0.0);
    }

    #[test]
    fn resistance_and_power_stay_passive_and_decrease_with_slip() {
        let mut previous = ResistanceFactors::new(1.0, 0.0).unwrap();
        for slip in [1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6] {
            let factors = ResistanceFactors::new(1.0, slip).unwrap();
            assert!((2.0 / 3.0..=1.0).contains(&factors.translation));
            assert!(factors.rotation > 0.0 && factors.rotation <= 1.0);
            assert!(factors.translation < previous.translation);
            assert!(factors.rotation < previous.rotation);
            previous = factors;
        }
    }

    #[test]
    fn thermal_growth_recomputes_the_ratio_without_rescaling_the_slip_length() {
        let cold = ResistanceFactors::new(1.0, 1.0).unwrap();
        let hot = ResistanceFactors::new(2.0, 1.0).unwrap();
        close(hot.translation, 0.8);
        close(hot.rotation, 0.4);
        assert_ne!(hot, cold);
    }

    #[test]
    fn invalid_and_unrepresentable_finite_boundaries_refuse() {
        for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ResistanceFactors::new(radius, 0.0).is_none());
        }
        for slip in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(ResistanceFactors::new(1.0, slip).is_none());
        }
        assert!(ResistanceFactors::new(f64::MIN_POSITIVE, f64::MAX).is_none());
    }
}
