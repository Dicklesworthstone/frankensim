//! Quasistatic supports on the existing layered uniaxial constitutive law.
use super::{ConcentricAxialResponse, ConcentricThermoelasticSection, checked_product, finite};
use crate::state_point::MaterialStatePointError;
use super::super::invalid;

/// Physical axial support data, independent of the material being installed.
/// These constraints describe equilibrium, not support dynamics or buckling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AxialConstraint {
    /// Fixed total engineering strain in the original stress-free configuration.
    FixedStrain(f64),
    /// Prescribed resultant [N], positive in tension. Zero is a freely expanding
    /// member; differential expansion can still leave nonzero layer stresses.
    FixedForce(f64),
    /// A linear elastic support with resultant
    /// `N = preload + stiffness * reference_length * (anchor_strain - strain)`.
    /// The anchor is stationary during this equilibrium evaluation.
    LinearSpring {
        /// Support stiffness [N/m]; explicit zero is a force-controlled limit.
        stiffness_n_m: f64,
        /// Common stress-free member length [m] defining engineering strain.
        reference_length_m: f64,
        /// Anchor displacement divided by that reference length.
        anchor_strain: f64,
        /// Signed force bias [N], distinct from spring stiffness.
        preload_n: f64,
    },
}

impl ConcentricThermoelasticSection {
    /// Solve the selected support equation, then evaluate EVERY layer through
    /// the original `response_at_axial_strain` law and its strain-domain checks.
    /// No effective expansion coefficient or material-specific support rule is
    /// introduced. Force/spring control differs physically from fixed strain.
    ///
    /// This is one-way quasistatic thermoelastic response at a supplied material
    /// state. Its stored elastic energy is not dissipated heat or a closed
    /// thermodynamic balance; changing supports may require external work.
    ///
    /// # Errors
    /// Invalid support data, unrepresentable arithmetic and any layer outside
    /// the declared linear regime refuse without returning a partial response.
    pub fn response_with_axial_constraint(
        &self, constraint: AxialConstraint, linear_strain_limit: f64,
    ) -> Result<ConcentricAxialResponse, MaterialStatePointError> {
        let ea = self.section().axial_stiffness_n();
        let thermal = self.thermal_force_n();
        let strain = match constraint {
            AxialConstraint::FixedStrain(strain) => finite(strain)?,
            AxialConstraint::FixedForce(force) => {
                let numerator = finite(finite(force)? + thermal)?;
                quotient(numerator, ea)?
            }
            AxialConstraint::LinearSpring {
                stiffness_n_m, reference_length_m, anchor_strain, preload_n,
            } => {
                if !stiffness_n_m.is_finite() || stiffness_n_m < 0.0
                    || !reference_length_m.is_finite() || reference_length_m <= 0.0
                {
                    return Err(invalid("axial support needs nonnegative stiffness and positive reference length"));
                }
                finite(anchor_strain)?;
                finite(preload_n)?;
                let spring = checked_product(stiffness_n_m, reference_length_m)?;
                let denominator = finite(ea + spring)?;
                // Divide before combining: a stiff support must not overflow
                // merely because spring*anchor_strain is formed prematurely.
                finite(quotient(thermal, denominator)? + quotient(preload_n, denominator)?
                    + checked_product(quotient(spring, denominator)?, anchor_strain)?)?
            }
        };
        self.response_at_axial_strain(strain, linear_strain_limit)
    }
}

fn quotient(a: f64, b: f64) -> Result<f64, MaterialStatePointError> {
    let value = finite(a / b)?;
    if a != 0.0 && value == 0.0 {
        return Err(invalid("axial support strain contribution underflowed"));
    }
    Ok(value)
}
