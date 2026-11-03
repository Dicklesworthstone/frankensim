//! Shared conservative convection weights for energy and staggered momentum.

use fs_math::det;

/// Convection discretization used by the finite-volume energy and momentum
/// equations. Each option preserves the same paired face-flux assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConvectionScheme {
    /// Patankar's inexpensive approximation to exponential fitting:
    /// `A(P) = max(0, (1 - 0.1 |P|)^5)`.
    #[default]
    PowerLaw,
    /// First-order upwind: `A(P) = 1`.
    Upwind,
    /// Conservative limited linear reconstruction of source specific enthalpy.
    ///
    /// Available only with an `AdvectiveEnthalpyField`. The positive upwind
    /// tangent matrix retains molecular diffusion; a frozen minmod face
    /// correction is refreshed by the original nonlinear material owner.
    /// Graded cell-centre distances define the slopes. Open boundaries and
    /// cells without two fluid neighbours along an axis retain donor transport.
    /// This is an energy scheme, not a momentum or temperature-only scheme.
    /// It does not establish a multidimensional discrete maximum principle.
    EnthalpyMinmod,
    /// Exponential fitting: `A(P) = |P| / expm1(|P|)`.
    ///
    /// Reproduces the source-free, constant-coefficient one-dimensional
    /// convection-diffusion profile between neighbouring unknowns, including
    /// graded centre distances and boundary half cells. The coefficients stay
    /// nonnegative at every Peclet number. This is not a higher-order claim
    /// for multidimensional crosswind transport or varying coefficients.
    /// More expensive than `PowerLaw`: each nonzero face evaluates an exponential.
    Exponential,
}

impl ConvectionScheme {
    pub(crate) fn weight(self, peclet: f64) -> f64 {
        match self {
            Self::Upwind | Self::EnthalpyMinmod => 1.0,
            Self::PowerLaw => {
                let t = 0.1_f64.mul_add(-peclet.abs(), 1.0).max(0.0);
                let t2 = t * t;
                t2 * t2 * t
            }
            Self::Exponential => {
                let p = peclet.abs();
                if p == 0.0 {
                    1.0
                } else if p == f64::INFINITY {
                    // Finite F/D can overflow when diffusion is negligible.
                    // Its fitted diffusion limit is zero, not inf * 0 = NaN.
                    0.0
                } else if p > 700.0 {
                    // Form the product in log space so exp(-p) does not
                    // underflow before multiplication by p can rescue it.
                    det::exp(det::ln(p) - p)
                } else if p > 50.0 {
                    // exp(-p) is below rounding at one, so 1-exp(-p) rounds
                    // to one. Avoid both positive-exponent overflow and a
                    // second transcendental evaluation.
                    p * det::exp(-p)
                } else {
                    // exp_m1 retains the central/conductive limit without
                    // subtracting two nearly equal floating-point numbers.
                    p / det::expm1(p)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::ConvectionScheme;

    #[test]
    fn exponential_weight_resolves_small_and_extreme_peclet_numbers() {
        let scheme = ConvectionScheme::Exponential;
        assert_eq!(scheme.weight(0.0), 1.0);
        assert_eq!(scheme.weight(-0.0), 1.0);
        for p in [
            f64::MIN_POSITIVE,
            1e-14,
            1e-8,
            1e-4,
            0.1,
            1.0,
            10.0,
            50.0,
            100.0,
            710.0,
            740.0,
            746.0,
            1000.0,
            f64::MAX,
            f64::INFINITY,
        ] {
            let weight = scheme.weight(p);
            assert!(weight.is_finite() && (0.0..=1.0).contains(&weight), "p={p}");
            assert_eq!(weight, scheme.weight(-p));
        }
        assert!((scheme.weight(1e-8) - (1.0 - 0.5e-8)).abs() < 2e-16);
        assert!((scheme.weight(1.0) - 0.581_976_706_869_326_5).abs() < 2e-16);
        assert!((scheme.weight(10.0) - 0.000_454_019_910_096_877_64).abs() < 2e-19);
        assert!(scheme.weight(710.0) > 0.0);
        assert!(scheme.weight(740.0) > 0.0);
        assert!(scheme.weight(746.0) > 0.0);
        assert_eq!(scheme.weight(1000.0), 0.0);
        assert!(scheme.weight(f64::NAN).is_nan());
    }

    #[test]
    fn exponential_faces_reproduce_the_constant_flux_on_a_graded_line() {
        // Independent exact solution of T'' - F T' = 0, k*A = L = 1,
        // T(0)=2, T(1)=7. Include boundary half cells and both flow signs.
        let points = [0.0, 0.017, 0.08, 0.31, 0.47, 0.89, 1.0];
        for flow in [-30.0_f64, -4.0, -1e-7, 0.0, 1e-7, 4.0, 30.0] {
            let temperature = |x: f64| {
                if flow == 0.0 {
                    2.0 + 5.0 * x
                } else {
                    2.0 + 5.0 * (flow * x).exp_m1() / flow.exp_m1()
                }
            };
            let exact_flux = if flow == 0.0 {
                -5.0
            } else {
                flow * (2.0 - 5.0 / flow.exp_m1())
            };
            for edge in points.windows(2) {
                let diffusion = 1.0 / (edge[1] - edge[0]);
                let weight = ConvectionScheme::Exponential.weight(flow / diffusion);
                let neighbour = diffusion.mul_add(weight, (-flow).max(0.0));
                let flux =
                    (neighbour + flow) * temperature(edge[0]) - neighbour * temperature(edge[1]);
                assert!(
                    (flux - exact_flux).abs() < 4e-13 * (1.0 + exact_flux.abs()),
                    "flow={flow}, edge={edge:?}: {flux} != {exact_flux}"
                );
            }
        }
    }
}
