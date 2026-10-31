//! A restricted bound on transport omitted by a finite diffuse path depth.

use fs_ivl::Interval;

/// A finite-depth remainder for an admitted ideal diffuse vacuum model.
///
/// This metadata is available only for uniform Lambertian calibrated area
/// emitters and passive tabulated [`crate::spectral::SpectralReflectance`]
/// materials in vacuum, with no other source or transport mechanism. If every
/// represented reflectance is at most `q < 1`, the positive diffuse transport
/// operator contracts the maximum spectral radiance by at most `q` per
/// reflection. A nonnegative spectral source envelope with photopic ceiling
/// `S` therefore bounds the omitted luminance by `S * q^D / (1 - q)`, where
/// `D` is the native maximum path depth. The illuminance bound is `pi` times
/// that ceiling. Using order `D` is conservative for native estimators whose
/// terminal direct-light proposal also includes the next reflection order.
///
/// The envelope uses the admitted uniform-emitter radiance coefficients,
/// retained source spectra and observer weights. The maximum reflectance is
/// taken over every knot of every admitted diffuse material, so the same
/// contraction applies to colored as well as gray materials.
///
/// These bounds concern only omitted transport in that ideal represented
/// model. They exclude sampling error, floating-point transport arithmetic,
/// geometry and source-normalization error, calibration and input uncertainty,
/// and disagreement between the represented model and a physical fixture.
/// They are not bounds on an individual sample or the sampled estimate, and
/// do not strengthen its numerical certificate or its unknown total error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiffuseVacuumRemainder {
    max_depth: u32,
    max_reflectance: f64,
    source_luminance_ceiling_cd_m2: f64,
    luminance_bound_cd_m2: f64,
    illuminance_bound_lux: f64,
}

impl DiffuseVacuumRemainder {
    /// The native maximum path depth used in the remainder calculation.
    #[must_use]
    pub const fn max_depth(self) -> u32 {
        self.max_depth
    }

    /// Upper bound on all represented diffuse reflectances, strictly below one.
    #[must_use]
    pub const fn max_reflectance(self) -> f64 {
        self.max_reflectance
    }

    /// Photopic ceiling of the admitted spectral source envelope in cd/m².
    ///
    /// This includes the retained spectra's represented observer integrals;
    /// it does not assume that their floating-point normalization is exact.
    #[must_use]
    pub const fn source_luminance_ceiling_cd_m2(self) -> f64 {
        self.source_luminance_ceiling_cd_m2
    }

    /// Upper bound in cd/m² on the ideal model's omitted diffuse transport.
    #[must_use]
    pub const fn luminance_bound_cd_m2(self) -> f64 {
        self.luminance_bound_cd_m2
    }

    /// Upper bound in lux on the ideal model's omitted hemisphere integral.
    #[must_use]
    pub const fn illuminance_bound_lux(self) -> f64 {
        self.illuminance_bound_lux
    }

    /// Evaluate a previously admitted scene envelope with outward arithmetic.
    ///
    /// Scene-wide source, material and vacuum eligibility must be checked by
    /// the caller. `source_luminance` encloses the nonnegative photopic source
    /// envelope, including the retained spectral normalization. Refuse invalid
    /// input or a bound whose radiance or illuminance cannot be represented by
    /// a finite `f64`; absence of this metadata does not reject a render job.
    pub(crate) fn from_envelope(
        max_depth: u32,
        max_reflectance: f64,
        source_luminance: Interval,
    ) -> Option<Self> {
        if max_depth == 0
            || !max_reflectance.is_finite()
            || !(0.0..1.0).contains(&max_reflectance)
            || !source_luminance.lo().is_finite()
            || !source_luminance.hi().is_finite()
            || source_luminance.lo() < 0.0
        {
            return None;
        }

        // Exact model zeros need no rounded arithmetic. Canonicalize signed
        // zero as well, so metadata has a deterministic positive-zero result.
        let max_reflectance = if max_reflectance == 0.0 {
            0.0
        } else {
            max_reflectance
        };
        let source_ceiling = if source_luminance.hi() == 0.0 {
            0.0
        } else {
            source_luminance.hi()
        };
        let (luminance_bound, illuminance_bound) =
            if max_reflectance == 0.0 || source_ceiling == 0.0 {
                (0.0, 0.0)
            } else {
                let power = reflectance_power(max_reflectance, max_depth);
                let denominator = Interval::point(1.0) - Interval::point(max_reflectance);
                if denominator.lo() <= 0.0 {
                    return None;
                }
                let luminance = (source_luminance * power) / denominator;
                let upper = luminance.hi();
                // In this branch the envelope's worst-case tail is positive.
                // A positive underflow must retain an outward subnormal
                // ceiling, never become an exact-zero model claim.
                if !upper.is_finite() || upper <= 0.0 {
                    return None;
                }
                let pi = Interval::new(
                    fs_math::next_down(std::f64::consts::PI),
                    fs_math::next_up(std::f64::consts::PI),
                );
                let lux = (Interval::point(upper) * pi).hi();
                if !lux.is_finite() || lux <= 0.0 {
                    return None;
                }
                (upper, lux)
            };

        Some(Self {
            max_depth,
            max_reflectance,
            source_luminance_ceiling_cd_m2: source_ceiling,
            luminance_bound_cd_m2: luminance_bound,
            illuminance_bound_lux: illuminance_bound,
        })
    }
}

/// At most 32 exponent bits even for the largest admitted native path depth.
fn reflectance_power(reflectance: f64, mut exponent: u32) -> Interval {
    let mut result = Interval::point(1.0);
    let mut base = Interval::point(reflectance);
    while exponent > 0 {
        if exponent & 1 != 0 {
            result = unit_product(result, base);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = unit_product(base, base);
        }
    }
    result
}

fn unit_product(left: Interval, right: Interval) -> Interval {
    let product = left * right;
    // Both factors represent powers of q in [0, 1]. Intersecting with this
    // known range is safe, and removes negative lower endpoints introduced
    // when a positive product rounds to zero. The positive outward upper
    // endpoint is retained even below the smallest representable subnormal.
    Interval::new(product.lo().max(0.0), product.hi().min(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dyadic_geometric_tail_uses_the_source_envelope_upper_endpoint() {
        let bound =
            DiffuseVacuumRemainder::from_envelope(3, 0.5, Interval::new(8.0, 16.0)).unwrap();
        // The worst admitted S is 16, so S * (1/2)^3 / (1 - 1/2) = 4.
        assert_eq!(bound.max_depth(), 3);
        assert_eq!(bound.max_reflectance(), 0.5);
        assert_eq!(bound.source_luminance_ceiling_cd_m2(), 16.0);
        assert!(bound.luminance_bound_cd_m2() >= 4.0);
        assert!(bound.luminance_bound_cd_m2() <= 4.0 * (1.0 + 1e-13));
        assert!(bound.illuminance_bound_lux() >= 4.0 * std::f64::consts::PI);
        assert!(bound.illuminance_bound_lux() <= 4.0 * std::f64::consts::PI * (1.0 + 1e-13));
    }

    #[test]
    fn reflectance_immediately_below_one_keeps_a_finite_enclosing_tail() {
        let q = fs_math::next_down(1.0);
        let bound = DiffuseVacuumRemainder::from_envelope(1, q, Interval::point(1.0)).unwrap();
        // q = 1 - 2^-53, so q / (1-q) = 2^53 - 1 exactly.
        let exact = 9_007_199_254_740_991.0;
        assert!(bound.luminance_bound_cd_m2() >= exact);
        assert!(bound.luminance_bound_cd_m2() <= exact * (1.0 + 1e-13));
        assert!(bound.illuminance_bound_lux().is_finite());
    }

    #[test]
    fn zero_reflectance_or_zero_source_has_an_exact_model_tail_zero() {
        for (q, source) in [(0.0, f64::MAX), (-0.0, 2.0), (0.75, 0.0), (0.75, -0.0)] {
            let bound = DiffuseVacuumRemainder::from_envelope(u32::MAX, q, Interval::point(source))
                .unwrap();
            assert_eq!(bound.luminance_bound_cd_m2().to_bits(), 0);
            assert_eq!(bound.illuminance_bound_lux().to_bits(), 0);
            if q == 0.0 {
                assert_eq!(bound.max_reflectance().to_bits(), 0);
            }
            if source == 0.0 {
                assert_eq!(bound.source_luminance_ceiling_cd_m2().to_bits(), 0);
            }
        }
    }

    #[test]
    fn positive_underflow_keeps_a_nonzero_upper_bound() {
        let tiny = f64::from_bits(1);
        for (depth, q, source) in [(2000, 0.5, 1.0), (1, 0.5, tiny), (2, tiny, 1.0)] {
            let bound =
                DiffuseVacuumRemainder::from_envelope(depth, q, Interval::point(source)).unwrap();
            assert!(bound.luminance_bound_cd_m2() > 0.0);
            assert!(bound.luminance_bound_cd_m2() < f64::MIN_POSITIVE);
            assert!(bound.illuminance_bound_lux() > 0.0);
            assert!(bound.illuminance_bound_lux() < f64::MIN_POSITIVE);
        }
    }

    #[test]
    fn maximum_native_depth_requires_only_logarithmic_work() {
        let bound =
            DiffuseVacuumRemainder::from_envelope(u32::MAX, 0.5, Interval::point(1.0)).unwrap();
        assert_eq!(bound.max_depth(), u32::MAX);
        assert!(bound.luminance_bound_cd_m2() > 0.0);
        assert!(bound.luminance_bound_cd_m2() < f64::MIN_POSITIVE);
        assert!(bound.illuminance_bound_lux().is_finite());
    }

    #[test]
    fn invalid_or_nonrepresentable_envelopes_are_not_reported_as_bounds() {
        assert!(DiffuseVacuumRemainder::from_envelope(0, 0.5, Interval::point(1.0)).is_none());
        for q in [-1.0, 1.0, 2.0, f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            assert!(DiffuseVacuumRemainder::from_envelope(1, q, Interval::point(1.0)).is_none());
        }
        for source in [
            Interval::new(-1.0, 1.0),
            Interval::new(-2.0, -1.0),
            Interval::new(0.0, f64::INFINITY),
            Interval::new(f64::NEG_INFINITY, 0.0),
            Interval::point(f64::INFINITY),
        ] {
            assert!(DiffuseVacuumRemainder::from_envelope(1, 0.5, source).is_none());
        }
        // The first envelope overflows radiance; the second still has finite
        // radiance but its hemisphere integral cannot fit in a finite f64.
        for (q, source) in [(0.9, f64::MAX / 4.0), (0.5, f64::MAX / 2.0)] {
            assert!(DiffuseVacuumRemainder::from_envelope(1, q, Interval::point(source)).is_none());
        }
    }
}
