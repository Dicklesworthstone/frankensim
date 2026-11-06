//! Newtonian viscous capillary instability of a slender cylindrical jet.
//!
//! This is the linearized one-dimensional model of Eggers & Dupont (1994),
//! equations (27)-(32), not the full viscous Rayleigh dispersion:
//! <https://people.maths.bris.ac.uk/~majge/ED94.pdf>.
//!
//! In the translating rest frame,
//! s^2 + 3 nu k^2 s = sigma k^2 (1 - R^2 k^2) / (2 rho R).
//! The factor 3 is the Newtonian extensional-viscosity contribution. The
//! coefficient assumptions are incompressible Newtonian liquid, constant
//! density/viscosity/tension, a uniform circular cylinder, small axisymmetric
//! disturbances, and a passive exterior with negligible inertia and viscosity.
//! Stretching, nonlinear necking, gas drag and non-Newtonian stresses are absent.
//!
//! The caller must prescribe a maximum kR in (0, 1]. Values near one retain
//! the model's neutral cutoff but need not satisfy a stringent long-wave
//! accuracy target. Setting a budget does not certify the approximation.
//! The zero-viscosity limit remains the slender approximation: its maximum
//! is kR=1/sqrt(2), distinct from the exact inviscid cylindrical maximum.
//! All evaluations are allocation-free and require constant work.

use super::{CapillaryJet, JetBandMaximum, JetError, JetMode};

/// Constant-property viscous slender cylinder with an explicit wavelength budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousSlenderJet {
    coefficients: CapillaryJet,
    kinematic_viscosity_m2_s: f64,
    maximum_kr: f64,
}

impl ViscousSlenderJet {
    /// Prescribe SI coefficients and the largest admitted dimensionless kR.
    ///
    /// Radius and density must be positive; tension and kinematic viscosity
    /// must be nonnegative; all inputs must be finite. Zero tension explicitly
    /// removes the capillary drive. No fluid-name lookup or hidden viscosity
    /// inference is performed.
    pub fn new(
        radius_m: f64,
        density_kg_m3: f64,
        surface_tension_n_m: f64,
        kinematic_viscosity_m2_s: f64,
        maximum_kr: f64,
    ) -> Result<Self, JetError> {
        if !kinematic_viscosity_m2_s.is_finite() || kinematic_viscosity_m2_s < 0.0 {
            return Err(JetError::InvalidInput("kinematic viscosity"));
        }
        if !maximum_kr.is_finite() || maximum_kr <= 0.0 || maximum_kr > 1.0 {
            return Err(JetError::InvalidInput("maximum kR"));
        }
        // Reuse validated physical coefficients and the capillary time scale,
        // not the inviscid owner's dispersion relation.
        let coefficients =
            CapillaryJet::inviscid(radius_m, density_kg_m3, surface_tension_n_m)?;
        Ok(Self {
            coefficients,
            kinematic_viscosity_m2_s,
            maximum_kr,
        })
    }

    /// Unperturbed circular radius [m].
    #[must_use]
    pub const fn radius_m(self) -> f64 {
        self.coefficients.radius_m()
    }

    /// Prescribed liquid density [kg/m^3].
    #[must_use]
    pub const fn density_kg_m3(self) -> f64 {
        self.coefficients.density_kg_m3()
    }

    /// Prescribed liquid/exterior tension [N/m].
    #[must_use]
    pub const fn surface_tension_n_m(self) -> f64 {
        self.coefficients.surface_tension_n_m()
    }

    /// Newtonian liquid kinematic viscosity [m^2/s].
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(self) -> f64 {
        self.kinematic_viscosity_m2_s
    }

    /// Caller-authored upper bound on kR, at most one.
    #[must_use]
    pub const fn maximum_kr(self) -> f64 {
        self.maximum_kr
    }

    /// Inverse capillary time sqrt(sigma/(rho R^3)) [1/s].
    #[must_use]
    pub const fn capillary_frequency_per_s(self) -> f64 {
        self.coefficients.capillary_frequency_per_s()
    }

    /// Evaluate the most amplified shape branch at physical k [rad/m].
    ///
    /// The admitted branch is growing for 0<kR<1 and neutral at zero,
    /// at kR=1, or when tension is zero. Its frequency is zero. The second,
    /// decaying velocity branch is not returned. Wavenumbers outside the
    /// declared long-wave budget are refused even when tension is zero.
    pub fn mode(self, wavenumber_rad_m: f64) -> Result<JetMode, JetError> {
        if !wavenumber_rad_m.is_finite() || wavenumber_rad_m < 0.0 {
            return Err(JetError::InvalidInput("wavenumber"));
        }
        let kr = wavenumber_rad_m * self.radius_m();
        if !kr.is_finite() || (wavenumber_rad_m > 0.0 && kr == 0.0) {
            return Err(JetError::NonRepresentable);
        }
        if kr > self.maximum_kr {
            return Err(JetError::InvalidInput("wavenumber exceeds maximum kR"));
        }
        let mut mode = JetMode {
            kr,
            growth_per_s: 0.0,
            angular_frequency_rad_s: 0.0,
        };
        if kr == 0.0 || kr == 1.0 || self.surface_tension_n_m() == 0.0 {
            return Ok(mode);
        }
        // g^2 is the inviscid slender drive. Factoring 1-x^2 retains the
        // small distance to the neutral cutoff without cancellation.
        let g = self.capillary_frequency_per_s()
            * kr
            * (0.5 * (1.0 - kr) * (1.0 + kr)).sqrt();
        positive(g)?;
        if self.kinematic_viscosity_m2_s == 0.0 {
            mode.growth_per_s = g;
            return Ok(mode);
        }
        let damping = ((self.kinematic_viscosity_m2_s * wavenumber_rad_m)
            * wavenumber_rad_m)
            * 1.5;
        // sqrt(a^2+g^2)-a = g^2/(sqrt(a^2+g^2)+a). Scale before
        // squaring and never subtract nearby numbers at high viscosity.
        let growth = if damping <= g {
            let q = damping / g;
            g / ((1.0 + q * q).sqrt() + q)
        } else {
            let q = g / damping;
            (g * q) / ((1.0 + q * q).sqrt() + 1.0)
        };
        mode.growth_per_s = if growth.is_finite() && growth > 0.0 {
            growth
        } else {
            // Rare exponent-range fallback: an overflowing damping scale
            // may still leave a representable small positive growth rate.
            // The normal physical range uses only arithmetic and sqrt.
            let log_g = fs_math::det::ln(g);
            let log_a = fs_math::det::ln(1.5)
                + fs_math::det::ln(self.kinematic_viscosity_m2_s)
                + 2.0 * fs_math::det::ln(wavenumber_rad_m);
            let log_growth = if log_a <= log_g {
                let q = fs_math::det::exp(log_a - log_g);
                log_g - fs_math::det::ln((1.0 + q * q).sqrt() + q)
            } else {
                let q = fs_math::det::exp(log_g - log_a);
                2.0 * log_g - log_a
                    - fs_math::det::ln((1.0 + q * q).sqrt() + 1.0)
            };
            positive(fs_math::det::exp(log_growth))?
        };
        Ok(mode)
    }

    /// Maximize growth over every k in a closed physical interval [rad/m].
    ///
    /// The dispersion is unimodal. For Oh=nu*sqrt(rho/(sigma*R)), its
    /// stationary point is x*=1/sqrt(2+3*sqrt(2)*Oh). Two endpoint
    /// evaluations and at most one interior evaluation give the maximum,
    /// with lower k winning equal computed growth. Both band endpoints
    /// must lie within the declared long-wave budget.
    pub fn fastest_mode(
        self,
        wavenumbers_rad_m: [f64; 2],
    ) -> Result<JetBandMaximum, JetError> {
        let [lower, upper] = wavenumbers_rad_m;
        if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || upper < lower {
            return Err(JetError::InvalidInput("wavenumber band"));
        }
        let low = self.mode(lower)?;
        let high = self.mode(upper)?;
        let mut maximum = JetBandMaximum {
            wavenumber_rad_m: lower,
            mode: low,
        };
        if self.surface_tension_n_m() > 0.0 && lower < upper {
            let critical = self.critical_kr()?;
            if low.kr < critical && critical < high.kr {
                let k = positive(critical / self.radius_m())?.clamp(lower, upper);
                let mode = self.mode(k)?;
                if mode.growth_per_s > maximum.mode.growth_per_s {
                    maximum = JetBandMaximum {
                        wavenumber_rad_m: k,
                        mode,
                    };
                }
            }
        }
        // Evaluate candidates in increasing k so an equal-growth plateau
        // cannot replace a lower-wavenumber candidate.
        if high.growth_per_s > maximum.mode.growth_per_s {
            maximum = JetBandMaximum {
                wavenumber_rad_m: upper,
                mode: high,
            };
        }
        Ok(maximum)
    }

    fn critical_kr(self) -> Result<f64, JetError> {
        if self.kinematic_viscosity_m2_s == 0.0 {
            return Ok(core::f64::consts::FRAC_1_SQRT_2);
        }
        // Oh=nu/(Omega*R^2), consistent with the cached frequency used by
        // mode(). Avoid squaring R and avoid overflowing 3*sqrt(2)*Oh.
        let oh = (self.kinematic_viscosity_m2_s / self.radius_m())
            / (self.capillary_frequency_per_s() * self.radius_m());
        let c = 3.0 * core::f64::consts::SQRT_2;
        if oh.is_finite() && oh > 0.0 {
            return positive(if oh <= 1.0 {
                1.0 / (2.0 + c * oh).sqrt()
            } else {
                (1.0 / oh.sqrt()) / (c + 2.0 / oh).sqrt()
            });
        }
        // A finite peak can survive an unrepresentable Oh intermediate.
        let log_oh = fs_math::det::ln(self.kinematic_viscosity_m2_s)
            - fs_math::det::ln(self.capillary_frequency_per_s())
            - 2.0 * fs_math::det::ln(self.radius_m());
        let log_term = fs_math::det::ln(c) + log_oh;
        let log_two = core::f64::consts::LN_2;
        let log_denominator = if log_term > log_two {
            log_term + fs_math::det::ln(1.0 + fs_math::det::exp(log_two - log_term))
        } else {
            log_two + fs_math::det::ln(1.0 + fs_math::det::exp(log_term - log_two))
        };
        positive(fs_math::det::exp(-0.5 * log_denominator))
    }
}

fn positive(value: f64) -> Result<f64, JetError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(JetError::NonRepresentable)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 5e-13 * expected.abs().max(f64::MIN_POSITIVE),
            "{actual:e} != {expected:e}"
        );
    }

    fn dimensionless(oh: f64, maximum_kr: f64) -> ViscousSlenderJet {
        ViscousSlenderJet::new(1.0, 1.0, 1.0, oh, maximum_kr).unwrap()
    }

    #[test]
    fn g1_satisfies_extensional_dispersion_and_neutral_cutoff() {
        for nu in [0.0, 1e-6, 1.0, 1e6] {
            let jet = dimensionless(nu, 1.0);
            for x in [0.0001, 0.05, 0.4, 0.9, 0.999_999] {
                let mode = jet.mode(x).unwrap();
                let s = mode.growth_per_s;
                near(
                    s * s + 3.0 * nu * x * x * s,
                    0.5 * x * x * (1.0 - x) * (1.0 + x),
                );
                assert!(s > 0.0);
                assert_eq!(mode.angular_frequency_rad_s, 0.0);
            }
            for x in [0.0, 1.0] {
                let mode = jet.mode(x).unwrap();
                assert_eq!(mode.growth_per_s, 0.0);
                assert_eq!(mode.angular_frequency_rad_s, 0.0);
            }
        }
    }

    #[test]
    fn g1_closed_form_continuous_maximum() {
        for oh in [0.0, 0.001, 0.1, 1.0, 10.0, 1e6, 1e200] {
            let jet = dimensionless(oh, 1.0);
            let peak = jet.fastest_mode([0.0, 1.0]).unwrap();
            let expected_x = 1.0 / (2.0 + 3.0 * core::f64::consts::SQRT_2 * oh).sqrt();
            let expected_growth = 1.0 / (8.0_f64.sqrt() + 6.0 * oh);
            near(peak.mode.kr, expected_x);
            near(peak.mode.growth_per_s, expected_growth);
            for fraction in [0.01, 0.1, 0.5, 0.99, 1.01, 1.5] {
                let x = peak.mode.kr * fraction;
                if x <= 1.0 {
                    assert!(
                        jet.mode(x).unwrap().growth_per_s
                            <= peak.mode.growth_per_s * (1.0 + 2e-14)
                    );
                }
            }
        }
    }

    #[test]
    fn g1_viscosity_slowdown_and_large_oh_avoid_cancellation() {
        let low = dimensionless(0.0, 1.0).mode(0.5).unwrap();
        let middle = dimensionless(1.0, 1.0).mode(0.5).unwrap();
        let high = dimensionless(1e200, 1.0).mode(0.5).unwrap();
        assert!(low.growth_per_s > middle.growth_per_s);
        assert!(middle.growth_per_s > high.growth_per_s);
        near(high.growth_per_s, 1.25e-201);
        // The slender zero-viscosity limit is explicitly distinct from the
        // exact cylindrical inviscid owner, while retaining the same cutoff.
        let exact = CapillaryJet::inviscid(1.0, 1.0, 1.0).unwrap();
        assert!(low.growth_per_s > exact.mode(0.5).unwrap().growth_per_s);
    }

    #[test]
    fn g1_band_edges_and_explicit_long_wave_budget() {
        let jet = dimensionless(1.0, 1.0);
        for (band, expected) in [
            ([0.0, 0.1], 0.1),
            ([0.9, 1.0], 0.9),
            ([0.4, 0.4], 0.4),
        ] {
            assert_eq!(jet.fastest_mode(band).unwrap().wavenumber_rad_m, expected);
        }
        let restricted = dimensionless(0.0, 0.3);
        assert_eq!(
            restricted.fastest_mode([0.0, 0.3]).unwrap().wavenumber_rad_m,
            0.3
        );
        assert!(restricted.mode(0.300_001).is_err());
        assert!(restricted.fastest_mode([0.0, 0.5]).is_err());
        assert!(jet.mode(1.000_001).is_err());
    }

    #[test]
    fn g3_fixed_oh_preserves_radius_density_and_tension_scalings() {
        let base = ViscousSlenderJet::new(0.001, 1000.0, 0.072, 0.001, 1.0).unwrap();
        let rate = base.mode(500.0).unwrap().growth_per_s;
        for (radius, density, tension, viscosity, k, factor) in [
            (0.004, 1000.0, 0.072, 0.002, 125.0, 0.125),
            (0.001, 4000.0, 0.072, 0.0005, 500.0, 0.5),
            (0.001, 1000.0, 0.648, 0.003, 500.0, 3.0),
        ] {
            let jet = ViscousSlenderJet::new(radius, density, tension, viscosity, 1.0).unwrap();
            near(jet.mode(k).unwrap().growth_per_s, factor * rate);
        }
    }

    #[test]
    fn g1_zero_tension_is_neutral_and_selects_lower_band_edge() {
        let jet = ViscousSlenderJet::new(1.0, 1.0, 0.0, f64::MAX, 1.0).unwrap();
        let maximum = jet.fastest_mode([0.2, 0.8]).unwrap();
        assert_eq!(maximum.wavenumber_rad_m, 0.2);
        assert_eq!(maximum.mode.growth_per_s, 0.0);
        assert_eq!(maximum.mode.angular_frequency_rad_s, 0.0);
        assert!(jet.mode(1.1).is_err());
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_refuse() {
        for value in [f64::NAN, f64::INFINITY, -1.0] {
            assert!(ViscousSlenderJet::new(1.0, 1.0, 1.0, value, 1.0).is_err());
        }
        for budget in [0.0, -1.0, 1.000_001, f64::NAN, f64::INFINITY] {
            assert!(ViscousSlenderJet::new(1.0, 1.0, 1.0, 0.1, budget).is_err());
        }
        for (radius, density, tension) in [
            (0.0, 1.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, -1.0),
            (f64::NAN, 1.0, 1.0),
            (1.0, f64::INFINITY, 1.0),
        ] {
            assert!(ViscousSlenderJet::new(radius, density, tension, 0.1, 1.0).is_err());
        }
        let jet = dimensionless(1.0, 1.0);
        for k in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(jet.mode(k).is_err());
        }
        assert!(jet.fastest_mode([0.8, 0.2]).is_err());
        assert!(jet.fastest_mode([0.0, f64::INFINITY]).is_err());
        let tiny = ViscousSlenderJet::new(1e-100, 1.0, 1.0, 0.0, 1.0).unwrap();
        assert_eq!(tiny.mode(1e-300), Err(JetError::NonRepresentable));
        let viscous = dimensionless(f64::MAX, 1.0);
        let x = 0.9;
        let expected = ((1.0 - x) * (1.0 + x) / f64::MAX) / 6.0;
        near(viscous.mode(x).unwrap().growth_per_s, expected);
        // This strictly unstable mode has a rate smaller than f64 can hold.
        let almost_one = f64::from_bits(1.0_f64.to_bits() - 1);
        assert_eq!(viscous.mode(almost_one), Err(JetError::NonRepresentable));
    }
}
