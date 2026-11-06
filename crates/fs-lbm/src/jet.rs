//! Axisymmetric capillary modes of a uniform inviscid cylindrical liquid jet.
//!
//! In its translating rest frame, an infinitesimal radius perturbation has
//! s^2 = sigma/(rho R^3) * x I1(x)/I0(x) * (1-x^2), where x=kR.
//! The azimuthal curvature is essential: a planar liquid strip is not this
//! problem. Modes 0<x<1 grow; x=0 and x=1 are neutral; x>1 oscillates.
//! Gas inertia/viscosity, liquid viscosity, gravity-induced stretching,
//! surfactants, finite-length boundary modes and nonlinear pinch-off are absent.
//! Callers must justify the inviscid approximation; the constructor does not
//! infer it from a liquid name or a viscosity hidden elsewhere in a scenario.
//!
//! Source: J. W. M. Bush, MIT 18.357, Lecture 11, equation (11.30):
//! <https://ocw.mit.edu/courses/18-357-interfacial-phenomena-fall-2010/d77da8d5f1bb8b69a62682bc0dbc5845_MIT18_357F10_Lecture11.pdf>.
//! These are deterministic floating-point model evaluations, not error
//! enclosures or experimental calibration. All queries are allocation-free.

/// Bind the original liquid source and explicit viscosity applicability.
pub mod material;

/// Newtonian viscous slender-jet model with an explicit long-wave budget.
pub mod viscous;

/// Fixed-property cylindrical model consumed by physical free-flight studies.
/// Implementations retain their own inviscid or long-wave admission boundaries.
/// An admitted model has a finite positive radius and finite nonnegative growth;
/// the returned maximum covers the entire requested band or refuses it.
pub trait JetDispersion: Copy {
    /// Unperturbed cylindrical radius in metres.
    fn radius_m(self) -> f64;

    /// Continuous maximum over the entire admitted physical wavenumber band.
    fn fastest_mode(self, wavenumbers_rad_m: [f64; 2]) -> Result<JetBandMaximum, JetError>;
}

impl JetDispersion for CapillaryJet {
    fn radius_m(self) -> f64 {
        CapillaryJet::radius_m(self)
    }

    fn fastest_mode(self, wavenumbers_rad_m: [f64; 2]) -> Result<JetBandMaximum, JetError> {
        CapillaryJet::fastest_mode(self, wavenumbers_rad_m)
    }
}

impl JetDispersion for viscous::ViscousSlenderJet {
    fn radius_m(self) -> f64 {
        viscous::ViscousSlenderJet::radius_m(self)
    }

    fn fastest_mode(self, wavenumbers_rad_m: [f64; 2]) -> Result<JetBandMaximum, JetError> {
        viscous::ViscousSlenderJet::fastest_mode(self, wavenumbers_rad_m)
    }
}

/// Most amplified capillary shape branch, measured in the jet's rest frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JetMode {
    /// Wavenumber times the unperturbed cylindrical radius.
    pub kr: f64,
    /// Temporal amplitude growth in 1/s; zero for neutral/oscillatory modes.
    pub growth_per_s: f64,
    /// Positive angular frequency of the oscillatory pair, in rad/s.
    /// Zero on the unstable and neutral branches; advection is not included.
    pub angular_frequency_rad_s: f64,
}

/// Continuous-band maximum of temporal amplitude growth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JetBandMaximum {
    /// Maximizing physical wavenumber [rad/m]; lower k wins equal growth.
    pub wavenumber_rad_m: f64,
    /// Native cylindrical mode at the maximizing wavenumber.
    pub mode: JetMode,
}

/// Admission failure; no clipped or substituted physics is returned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JetError {
    /// A named input lies outside its declared SI domain.
    InvalidInput(&'static str),
    /// A required nonzero value underflowed, or a derived value overflowed.
    NonRepresentable,
}

impl core::fmt::Display for JetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(name) => write!(f, "invalid capillary-jet input: {name}"),
            Self::NonRepresentable => f.write_str("capillary-jet state is not representable"),
        }
    }
}

impl std::error::Error for JetError {}

/// Admitted constant-property inviscid cylinder in SI units.
///
/// Construction computes its capillary frequency once. Queries neither
/// allocate nor evolve a grid; they evaluate the exact inviscid linear
/// dispersion with floating-point Bessel ratios.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryJet {
    radius_m: f64,
    density_kg_m3: f64,
    surface_tension_n_m: f64,
    capillary_frequency_per_s: f64,
}

impl CapillaryJet {
    /// Construct an explicitly inviscid jet with positive radius and density.
    ///
    /// Zero surface tension explicitly neglects capillarity and makes every
    /// mode neutral. A nonzero coefficient must give a representable capillary
    /// frequency. This API accepts prescribed coefficients; it supplies no
    /// material-source authority and does not silently ignore a viscosity input.
    pub fn inviscid(
        radius_m: f64,
        density_kg_m3: f64,
        surface_tension_n_m: f64,
    ) -> Result<Self, JetError> {
        for (name, value, positive) in [
            ("radius", radius_m, true),
            ("density", density_kg_m3, true),
            ("surface tension", surface_tension_n_m, false),
        ] {
            if !value.is_finite() || value < 0.0 || (positive && value == 0.0) {
                return Err(JetError::InvalidInput(name));
            }
        }
        // Avoid cubing R or forming sigma/rho before taking the square root.
        let capillary_frequency_per_s = if surface_tension_n_m == 0.0 {
            0.0
        } else {
            fs_math::det::exp(
                0.5 * (fs_math::det::ln(surface_tension_n_m)
                    - fs_math::det::ln(density_kg_m3)
                    - 3.0 * fs_math::det::ln(radius_m)),
            )
        };
        represented(capillary_frequency_per_s, surface_tension_n_m > 0.0)?;
        Ok(Self {
            radius_m,
            density_kg_m3,
            surface_tension_n_m,
            capillary_frequency_per_s,
        })
    }

    /// Unperturbed cylindrical radius [m].
    #[must_use]
    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }

    /// Prescribed liquid density [kg/m^3].
    #[must_use]
    pub const fn density_kg_m3(self) -> f64 {
        self.density_kg_m3
    }

    /// Prescribed liquid/exterior surface tension [N/m].
    #[must_use]
    pub const fn surface_tension_n_m(self) -> f64 {
        self.surface_tension_n_m
    }

    /// Inverse capillary time sqrt(sigma/(rho R^3)) [1/s].
    #[must_use]
    pub const fn capillary_frequency_per_s(self) -> f64 {
        self.capillary_frequency_per_s
    }

    /// Evaluate a nonnegative physical wavenumber in rad/m.
    ///
    /// A stable mode returns its oscillation frequency instead of taking a
    /// square root of negative growth or falsely reporting exponential decay.
    /// The unstable result is the positive member of the real eigenvalue pair.
    pub fn mode(self, wavenumber_rad_m: f64) -> Result<JetMode, JetError> {
        if !wavenumber_rad_m.is_finite() || wavenumber_rad_m < 0.0 {
            return Err(JetError::InvalidInput("wavenumber"));
        }
        let x = wavenumber_rad_m * self.radius_m;
        represented(x, wavenumber_rad_m > 0.0)?;
        let mut mode = JetMode {
            kr: x,
            growth_per_s: 0.0,
            angular_frequency_rad_s: 0.0,
        };
        if x == 0.0 || x == 1.0 || self.surface_tension_n_m == 0.0 {
            return Ok(mode);
        }
        let ratio = bessel_i1_over_i0(x);
        represented(ratio, true)?;
        // Split the roots instead of forming x^3 at large k, and factor
        // 1-x^2 near the neutral cutoff to retain its small signed distance.
        let rate = self.capillary_frequency_per_s
            * x.sqrt()
            * ratio.sqrt()
            * (1.0 - x).abs().sqrt()
            * (1.0 + x).sqrt();
        represented(rate, true)?;
        if x < 1.0 {
            mode.growth_per_s = rate;
        } else {
            mode.angular_frequency_rad_s = rate;
        }
        Ok(mode)
    }

    /// Maximize growth over EVERY wavenumber in a closed SI interval.
    ///
    /// The inviscid dispersion is unimodal on 0<kR<1. Its dimensionless
    /// stationary point is independent of R, rho and sigma, so two endpoint
    /// evaluations and at most one interior evaluation replace a sampled scan.
    /// All stable modes have zero growth; the lower endpoint wins that tie.
    /// Both endpoints must be numerically representable, including a stable
    /// tail that does not contain the maximum.
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
        if high.growth_per_s > maximum.mode.growth_per_s {
            maximum = JetBandMaximum {
                wavenumber_rad_m: upper,
                mode: high,
            };
        }
        // Rounded 70-decimal positive-series root of
        // (1-x^2)*(1-r^2)-2*x*r=0, r=I1(x)/I0(x). This is the derivative
        // condition, not a fitted experimental wavelength or a grid sample.
        const PEAK_KR: f64 = 0.697_018_898_306_328_6;
        if self.surface_tension_n_m > 0.0 && low.kr < PEAK_KR && PEAK_KR < high.kr {
            let k = PEAK_KR / self.radius_m;
            represented(k, true)?;
            let mode = self.mode(k)?;
            if mode.growth_per_s > maximum.mode.growth_per_s
                || (mode.growth_per_s == maximum.mode.growth_per_s
                    && k < maximum.wavenumber_rad_m)
            {
                maximum = JetBandMaximum {
                    wavenumber_rad_m: k,
                    mode,
                };
            }
        }
        Ok(maximum)
    }
}

// Evaluate the ratio directly, never the exponentially large Bessel functions.
// The unstable branch needs only the short positive power series. A backward
// continued fraction covers 1<x<=64; above it eight asymptotic terms have
// ~3e-15 truncation at the switch and decreasing error toward infinity.
fn bessel_i1_over_i0(x: f64) -> f64 {
    if x <= 1.0 {
        let z = 0.25 * x * x;
        let (mut a, mut b, mut sum_a, mut sum_b) = (1.0, 1.0, 1.0, 1.0);
        for n in 1..=12 {
            let n = f64::from(n);
            a *= z / (n * n);
            b *= z / (n * (n + 1.0));
            sum_a += a;
            sum_b += b;
        }
        return (0.5 * x) * (sum_b / sum_a);
    }
    if x <= 64.0 {
        let square = x * x;
        let mut tail = 0.0;
        for n in (2..=128).rev() {
            tail = square / (2.0 * f64::from(n) + tail);
        }
        return x / (2.0 + tail);
    }
    // Coefficients follow from r'=1-r/x-r^2 with r -> 1 as x -> infinity.
    const C: [f64; 8] = [
        -0.5,
        -0.125,
        -0.125,
        -25.0 / 128.0,
        -13.0 / 32.0,
        -1073.0 / 1024.0,
        -103.0 / 32.0,
        -375_733.0 / 32_768.0,
    ];
    let inverse = 1.0 / x;
    let polynomial = C
        .iter()
        .rev()
        .fold(0.0, |tail, &coefficient| inverse.mul_add(tail, coefficient));
    inverse.mul_add(polynomial, 1.0)
}

fn represented(value: f64, nonzero: bool) -> Result<(), JetError> {
    if !value.is_finite() || (nonzero && value == 0.0) {
        Err(JetError::NonRepresentable)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-13 * expected.abs().max(f64::MIN_POSITIVE),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g1_cylindrical_cutoff_and_stable_oscillations() {
        let jet = CapillaryJet::inviscid(1.0, 1.0, 1.0).unwrap();
        for k in [0.0, 1.0] {
            let mode = jet.mode(k).unwrap();
            assert_eq!(mode.growth_per_s, 0.0);
            assert_eq!(mode.angular_frequency_rad_s, 0.0);
        }
        let growing = jet.mode(0.5).unwrap();
        // Independent positive-series I1(0.5)/I0(0.5) reference.
        near(
            growing.growth_per_s,
            (0.5_f64 * 0.242_499_612_580_801_95 * 0.75).sqrt(),
        );
        assert_eq!(growing.angular_frequency_rad_s, 0.0);
        let stable = jet.mode(2.0).unwrap();
        near(
            stable.angular_frequency_rad_s,
            (6.0_f64 * 0.697_774_657_964_008).sqrt(),
        );
        assert_eq!(stable.growth_per_s, 0.0);
        // The planar high-k limit omega^2 -> sigma*k^3/rho remains finite
        // here even though forming k^3 would overflow.
        near(jet.mode(1e100).unwrap().angular_frequency_rad_s, 1e150);
    }

    #[test]
    fn g1_continuous_maximum_and_band_edges() {
        let jet = CapillaryJet::inviscid(1.0, 1.0, 1.0).unwrap();
        let peak = jet.fastest_mode([0.0, 2.0]).unwrap();
        near(peak.mode.kr, 0.697_018_898_306_328_6);
        near(peak.mode.growth_per_s, 0.343_338_701_861_299_3);
        let x = peak.mode.kr;
        let r = bessel_i1_over_i0(x);
        assert!(((1.0 - x * x) * (1.0 - r * r) - 2.0 * x * r).abs() < 2e-15);
        for k in [0.0, 0.1, 0.3, 0.6, 0.69, 0.7, 0.8, 1.0, 2.0] {
            assert!(jet.mode(k).unwrap().growth_per_s <= peak.mode.growth_per_s);
        }
        for (band, expected_k) in [
            ([0.1, 0.3], 0.3),
            ([0.8, 1.2], 0.8),
            ([1.2, 3.0], 1.2),
            ([0.4, 0.4], 0.4),
        ] {
            assert_eq!(jet.fastest_mode(band).unwrap().wavenumber_rad_m, expected_k);
        }
    }

    #[test]
    fn g3_radius_density_and_tension_scalings() {
        let reference = CapillaryJet::inviscid(0.001, 1000.0, 0.072).unwrap();
        let base = reference.mode(500.0).unwrap().growth_per_s;
        for (radius, density, tension, k, factor) in [
            (0.004, 1000.0, 0.072, 125.0, 0.125),
            (0.001, 4000.0, 0.072, 500.0, 0.5),
            (0.001, 1000.0, 0.648, 500.0, 3.0),
        ] {
            let jet = CapillaryJet::inviscid(radius, density, tension).unwrap();
            near(jet.mode(k).unwrap().growth_per_s, factor * base);
        }
        let none = CapillaryJet::inviscid(1.0, 1.0, 0.0).unwrap();
        let neutral = none.fastest_mode([0.2, 3.0]).unwrap();
        assert_eq!(neutral.wavenumber_rad_m, 0.2);
        assert_eq!(neutral.mode.growth_per_s, 0.0);
        assert_eq!(neutral.mode.angular_frequency_rad_s, 0.0);
    }

    #[test]
    fn g1_bessel_ratio_across_numerical_ranges() {
        for (x, expected) in [
            (0.1, 0.049_937_603_987_938_92),
            (1.0, 0.446_389_965_896_534_5),
            (8.0, 0.935_235_493_529_438_6),
            (32.0, 0.984_248_915_536_987),
            (64.0, 0.992_156_493_548_811_2),
            (65.0, 0.992_277_640_036_529_7),
            (128.0, 0.996_086_060_261_160_6),
        ] {
            near(bessel_i1_over_i0(x), expected);
        }
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_refuse() {
        for (radius, density, tension) in [
            (0.0, 1.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, -1.0),
            (f64::NAN, 1.0, 1.0),
            (1.0, f64::INFINITY, 1.0),
        ] {
            assert!(CapillaryJet::inviscid(radius, density, tension).is_err());
        }
        assert_eq!(
            CapillaryJet::inviscid(1e-250, 1.0, 1.0),
            Err(JetError::NonRepresentable)
        );
        let jet = CapillaryJet::inviscid(1.0, 1.0, 1.0).unwrap();
        assert!(jet.mode(-1.0).is_err());
        assert!(jet.mode(f64::NAN).is_err());
        assert!(jet.fastest_mode([2.0, 1.0]).is_err());
        assert!(jet.fastest_mode([0.0, f64::INFINITY]).is_err());
    }
}
