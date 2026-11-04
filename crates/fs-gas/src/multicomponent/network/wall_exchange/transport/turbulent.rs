//! Developed turbulent circular-pipe friction on the original gas-wall owner.

use super::{Length, WallExchangeError, admit, positive};
use fs_math::det;

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
