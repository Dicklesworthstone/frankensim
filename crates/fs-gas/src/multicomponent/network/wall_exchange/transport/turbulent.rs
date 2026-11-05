//! Developed turbulent circular-pipe friction on the original gas-wall owner.

use super::{Length, WallExchangeError, admit, positive};
use fs_math::det;

/// Coupled smooth-pipe friction and Gnielinski turbulent heat transfer.
///
/// Bulk Reynolds and Prandtl numbers are recomputed from the original gas
/// composition, calorics, viscosity and conductivity at every source interval.
/// This explicitly selected correlation requires fully developed flow in a
/// smooth circular pipe, modest wall-to-bulk property variation and low Mach
/// number. It does not infer development from Reynolds number, model roughness,
/// or interpolate through transition. The original branch must be at least ten
/// diameters long; that geometric check alone does not establish development.
///
/// The supported range is deliberately restricted to Re in [1e4, 1e6] and Pr
/// in [0.5, 200]. Callers may narrow both ranges. The same Darcy factor drives
/// momentum and heat; no constant Nusselt input is needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularGnielinskiWall {
    /// Positive diameter of the original constant-area circular branch [m].
    pub diameter: Length,
    /// Inclusive declared Reynolds interval, contained in [1e4, 1e6].
    pub reynolds_number: [f64; 2],
    /// Inclusive declared Prandtl interval, contained in [0.5, 200].
    pub prandtl_number: [f64; 2],
    /// Positive upper wall-relative Mach bound, no greater than 0.3.
    pub maximum_mach_number: f64,
}

impl CircularGnielinskiWall {
    pub(super) fn validate(self) -> Result<(), WallExchangeError> {
        let [re_min, re_max] = self.reynolds_number;
        let [pr_min, pr_max] = self.prandtl_number;
        if !positive(self.diameter.value())
            || !re_min.is_finite()
            || !re_max.is_finite()
            || re_min < 1e4
            || re_min > re_max
            || re_max > 1e6
            || !pr_min.is_finite()
            || !pr_max.is_finite()
            || pr_min < 0.5
            || pr_min > pr_max
            || pr_max > 200.0
            || !positive(self.maximum_mach_number)
            || self.maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "Gnielinski circular wall geometry or applicability",
            ));
        }
        Ok(())
    }

    /// Return `(Darcy friction factor, Nusselt number)` from the same state.
    ///
    /// `f = (0.79 ln(Re) - 1.64)^-2`,
    /// `Nu = (f/8)(Re-1000)Pr / [1 + 12.7 sqrt(f/8)(Pr^(2/3)-1)]`.
    /// Reference: Gnielinski (1976), International Chemical Engineering 16(2),
    /// 359--368, turbulent pipe-flow heat-transfer correlation.
    /// This allocation-free evaluation shares the smooth-pipe relation used by
    /// `fs-convection` without introducing its FEM dependencies into gas flow.
    pub fn coefficients(
        self,
        reynolds: f64,
        prandtl: f64,
    ) -> Result<(f64, f64), WallExchangeError> {
        self.validate()?;
        admit(
            "Reynolds number",
            reynolds,
            self.reynolds_number[0],
            self.reynolds_number[1],
        )?;
        admit(
            "Prandtl number",
            prandtl,
            self.prandtl_number[0],
            self.prandtl_number[1],
        )?;
        let inverse = 0.79_f64.mul_add(det::ln(reynolds), -1.64);
        let friction = 1.0 / (inverse * inverse);
        let eighth = friction / 8.0;
        let nusselt = eighth * (reynolds - 1000.0) * prandtl
            / (1.0 + 12.7 * det::sqrt(eighth) * (det::pow(prandtl, 2.0 / 3.0) - 1.0));
        if !positive(friction) || !positive(nusselt) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable Gnielinski coefficients",
            ));
        }
        Ok((friction, nusselt))
    }
}

/// Colebrook-White Darcy friction with explicit equivalent sand-grain roughness.
///
/// The original mixture supplies viscosity, density and wall-relative speed at
/// every source interval. This engineering closure admits developed circular
/// pipe flow only: Re in [4000, 1e8], equivalent roughness/diameter in [0, 0.05],
/// and a caller-declared Mach ceiling no greater than 0.3. The caller may narrow
/// these limits. A transition out of that domain refuses the complete step;
/// no laminar/turbulent blend or roughness inferred from a material name is used.
///
/// The Nusselt number remains an independent, authored heat-transfer closure.
/// Friction does not supply a heat-transfer correlation or resolve entrance,
/// oscillatory, separated, compressible boundary-layer or non-Newtonian flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularTurbulentWall {
    /// Positive diameter of the original constant-area circular branch [m].
    pub diameter: Length,
    /// Nonnegative equivalent sand-grain roughness [m], not an RMS height.
    pub equivalent_sand_roughness: Length,
    /// Nonnegative declared `h D / k`; zero disables radial heat exchange.
    pub nusselt_number: f64,
    /// Inclusive supported Reynolds interval, contained in [4000, 1e8].
    pub reynolds_number: [f64; 2],
    /// Positive upper wall-relative Mach bound, no greater than 0.3.
    pub maximum_mach_number: f64,
}

impl CircularTurbulentWall {
    pub(super) fn validate(self) -> Result<(), WallExchangeError> {
        let [minimum, maximum] = self.reynolds_number;
        if !positive(self.diameter.value())
            || !self.equivalent_sand_roughness.value().is_finite()
            || self.equivalent_sand_roughness.value() < 0.0
            || !self.nusselt_number.is_finite()
            || self.nusselt_number < 0.0
            || !minimum.is_finite()
            || !maximum.is_finite()
            || minimum < 4000.0
            || minimum > maximum
            || maximum > 1e8
            || !positive(self.maximum_mach_number)
            || self.maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "turbulent circular wall geometry or applicability",
            ));
        }
        // Compare lengths before division so a boundary authored as 0.05*D
        // is not spuriously outside its own limit after one rounded quotient.
        if self.equivalent_sand_roughness.value() > 0.05 * self.diameter.value() {
            return Err(WallExchangeError::TransportDomain {
                quantity: "equivalent sand roughness / diameter",
                value: self.equivalent_sand_roughness.value() / self.diameter.value(),
                minimum: 0.0,
                maximum: 0.05,
            });
        }
        Ok(())
    }

    /// Solve the original implicit Colebrook relation for the Darcy factor.
    ///
    /// `1/sqrt(f) = -2 log10(epsilon/(3.7 D) + 2.51/(Re sqrt(f)))`.
    /// Reference: Colebrook (1939), DOI 10.1680/ijoti.1939.13150;
    /// NASA KSC-STD-Z-0017 Rev A, section 7.4.1, equation 4:
    /// <https://standards.nasa.gov/sites/default/files/standards/KSC/A/0/KSC-STD-Z-0017_RevA_Final.pdf>.
    ///
    /// The reciprocal square-root variable gives a strictly increasing scalar
    /// residual. Safeguarded Newton uses deterministic elementary functions,
    /// no allocation, a finite bracket and at most 32 iterations. The numerical
    /// solve is much tighter than the accuracy of the empirical correlation;
    /// it does not promote that correlation to a resolved turbulence model.
    pub fn darcy_factor(self, reynolds: f64) -> Result<f64, WallExchangeError> {
        self.validate()?;
        admit(
            "Reynolds number",
            reynolds,
            self.reynolds_number[0],
            self.reynolds_number[1],
        )?;
        let a = (self.equivalent_sand_roughness.value() / self.diameter.value()) / 3.7;
        let b = 2.51 / reynolds;
        let scale = 2.0 * core::f64::consts::LOG10_E;
        // Over the admitted domain, F(1)<0 and F(32)>0, including smooth pipe.
        let (mut lower, mut upper) = (1.0, 32.0);
        let mut x = 7.0;
        for _ in 0..32 {
            let argument = a + b * x;
            let residual = x + scale * det::ln(argument);
            if residual.abs() <= 64.0 * f64::EPSILON * x {
                return Ok(1.0 / (x * x));
            }
            if residual < 0.0 {
                lower = x;
            } else {
                upper = x;
            }
            let derivative = 1.0 + scale * b / argument;
            let newton = x - residual / derivative;
            x = if newton > lower && newton < upper {
                newton
            } else {
                0.5 * (lower + upper)
            };
        }
        Err(WallExchangeError::Invalid(
            "Colebrook friction solve did not converge",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn heat_model() -> CircularGnielinskiWall {
        CircularGnielinskiWall {
            diameter: Length::new(0.02),
            reynolds_number: [1e4, 1e6],
            prandtl_number: [0.5, 200.0],
            maximum_mach_number: 0.3,
        }
    }

    #[test]
    fn gnielinski_g1_friction_heat_relation_and_prandtl_dependence() {
        for re in [1e4_f64, 1e5, 1e6] {
            let expected_f = (0.79 * re.ln() - 1.64).powi(-2);
            let mut previous_nu = 0.0;
            for pr in [0.5_f64, 1.0, 7.0, 200.0] {
                let (f, nu) = heat_model().coefficients(re, pr).unwrap();
                // Independent platform functions and the original equation;
                // at Pr=1 the denominator is exactly one (Reynolds analogy).
                let expected_nu = (expected_f / 8.0) * (re - 1000.0) * pr
                    / (1.0 + 12.7 * (expected_f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
                assert!((f / expected_f - 1.0).abs() < 2e-13);
                assert!((nu / expected_nu - 1.0).abs() < 2e-13);
                assert!(nu > previous_nu);
                previous_nu = nu;
                if pr == 1.0 {
                    assert_eq!(nu, (f / 8.0) * (re - 1000.0));
                }
            }
        }
        let (_, slow) = heat_model().coefficients(1e4, 0.7).unwrap();
        let (_, fast) = heat_model().coefficients(1e5, 0.7).unwrap();
        assert!(fast > slow);
    }

    #[test]
    fn gnielinski_g0_requires_authored_supported_thermal_and_flow_regimes() {
        for re in [9999.0, 1e6 + 1.0, f64::NAN, f64::INFINITY] {
            assert!(heat_model().coefficients(re, 1.0).is_err());
        }
        for pr in [0.0, 0.49, 200.01, f64::NAN, f64::INFINITY] {
            assert!(heat_model().coefficients(1e5, pr).is_err());
        }
        let narrowed = CircularGnielinskiWall {
            reynolds_number: [2e4, 2e5],
            prandtl_number: [0.7, 7.0],
            ..heat_model()
        };
        assert!(narrowed.coefficients(1e4, 1.0).is_err());
        assert!(narrowed.coefficients(1e5, 0.5).is_err());
        assert!(narrowed.coefficients(2e4, 0.7).is_ok());
        for invalid in [
            CircularGnielinskiWall {
                diameter: Length::new(0.0),
                ..heat_model()
            },
            CircularGnielinskiWall {
                reynolds_number: [4000.0, 1e6],
                ..heat_model()
            },
            CircularGnielinskiWall {
                prandtl_number: [0.5, 201.0],
                ..heat_model()
            },
            CircularGnielinskiWall {
                maximum_mach_number: 0.4,
                ..heat_model()
            },
        ] {
            assert!(invalid.coefficients(1e5, 1.0).is_err());
        }
    }

    fn model(relative_roughness: f64) -> CircularTurbulentWall {
        CircularTurbulentWall {
            diameter: Length::new(1.0),
            equivalent_sand_roughness: Length::new(relative_roughness),
            nusselt_number: 0.0,
            reynolds_number: [4000.0, 1e8],
            maximum_mach_number: 0.3,
        }
    }

    #[test]
    fn colebrook_g1_original_equation_and_rough_pipe_limit() {
        for roughness in [0.0, 1e-8, 1e-4, 0.01, 0.05] {
            let pipe = model(roughness);
            let mut previous = f64::INFINITY;
            for reynolds in [4000.0, 1e4, 1e5, 1e6, 1e8] {
                let factor = pipe.darcy_factor(reynolds).unwrap();
                // Independent platform-log check in the original f variable.
                let inverse = 1.0 / factor.sqrt();
                let residual =
                    inverse + 2.0 * (roughness / 3.7 + 2.51 * inverse / reynolds).log10();
                assert!(residual.abs() < 2e-12, "{roughness} {reynolds}: {residual}");
                assert!(factor > 0.0 && factor < previous);
                previous = factor;
            }
        }
        let rough = model(0.05).darcy_factor(1e8).unwrap();
        let fully_rough = 1.0 / (-2.0 * (0.05_f64 / 3.7).log10()).powi(2);
        assert!((rough / fully_rough - 1.0).abs() < 5e-6);
        let smooth = model(0.0).darcy_factor(1e5).unwrap();
        assert!((smooth - 0.017_989_773_084_273_84).abs() < 2e-14);
        assert!(model(1e-4).darcy_factor(1e5).unwrap() > smooth);
        for diameter in [0.2, 1.0, 10.0] {
            let scaled = CircularTurbulentWall {
                diameter: Length::new(diameter),
                equivalent_sand_roughness: Length::new(0.05 * diameter),
                ..model(0.05)
            };
            assert!((scaled.darcy_factor(1e8).unwrap() - rough).abs() < 2e-15);
        }
    }

    #[test]
    fn colebrook_g0_refuses_transition_extrapolation_and_invalid_surface() {
        for reynolds in [0.0, 2300.0, 3999.0, 1e8 + 1.0, f64::NAN, f64::INFINITY] {
            assert!(model(0.001).darcy_factor(reynolds).is_err());
        }
        for roughness in [-0.001, 0.051, f64::NAN, f64::INFINITY] {
            assert!(model(roughness).darcy_factor(1e5).is_err());
        }
        let mut pipe = model(0.001);
        pipe.reynolds_number = [1e5, 1e6];
        assert!(pipe.darcy_factor(1e4).is_err());
        assert!(pipe.darcy_factor(1e5).is_ok());
        pipe.maximum_mach_number = 0.4;
        assert!(pipe.darcy_factor(1e5).is_err());
    }
}
