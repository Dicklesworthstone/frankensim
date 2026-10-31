//! Source-retained cycle-equivalent bending loss for circular modal members.
//!
//! The same section loss-angle and circular thermoelastic owners used by piano
//! strings supply every retained mode's intrinsic dissipation. This wrapper
//! preserves their original material states, selectors and frequency support;
//! it does not infer a viscosity or a causal broadband relaxation spectrum.
//! The modal rate is H(omega) k^4 / (mu omega), with the actual wave number,
//! angular frequency and dynamic line mass supplied by the mechanical owner.
//!
//! Thermoelastic loss alone omits unknown anelastic excess. An explicitly
//! selected excess must belong to exactly the same elastic section and remain
//! inside its source band; it is never silently dropped outside that band.

use super::loss_angle::ConcentricLossAngleSection;
use super::ted::ConcentricThermoelasticBending;
use super::{ConcentricElasticSection, MaterialStatePointError, invalid, positive};
use crate::visco::circular_ted::circular_relaxation_shape;

/// One original section's loss-angle or computed thermal bending dissipation.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricModalLossSection {
    law: ModalBendingLoss,
    omega_band_rad_s: (f64, f64),
}

#[derive(Clone, Debug, PartialEq)]
enum ModalBendingLoss {
    LossAngle(ConcentricLossAngleSection),
    Thermoelastic {
        thermal: ConcentricThermoelasticBending,
        excess: Option<ConcentricLossAngleSection>,
    },
}

impl ConcentricModalLossSection {
    /// Retain every layer's measured loss tangent and intersected source support.
    #[must_use]
    pub fn from_loss_angle(section: ConcentricLossAngleSection) -> Self {
        Self { omega_band_rad_s: section.omega_band_rad_s(),
            law: ModalBendingLoss::LossAngle(section) }
    }

    /// Retain computed circular thermoelastic loss and optional measured excess.
    ///
    /// With no excess, this is the intrinsic thermoelastic contribution only,
    /// not a claim that all other material loss vanishes.
    ///
    /// # Errors
    /// Excess attached to another geometry/source state or disjoint source bands.
    pub fn from_thermoelastic(
        thermal: ConcentricThermoelasticBending,
        excess: Option<ConcentricLossAngleSection>,
    ) -> Result<Self, MaterialStatePointError> {
        let mut band = thermal.omega_band_rad_s();
        if let Some(loss) = &excess {
            if loss.elastic() != thermal.elastic() {
                return Err(invalid("thermoelastic excess must retain exactly the same elastic section and material states"));
            }
            let interval = loss.omega_band_rad_s();
            band = (band.0.max(interval.0), band.1.min(interval.1));
        }
        if band.1 < band.0 {
            return Err(invalid("thermoelastic and selected excess source bands do not intersect"));
        }
        Ok(Self { law: ModalBendingLoss::Thermoelastic { thermal, excess },
            omega_band_rad_s: band })
    }

    /// Original section, including each layer's geometry and source receipts.
    #[must_use]
    pub fn elastic(&self) -> &ConcentricElasticSection {
        match &self.law {
            ModalBendingLoss::LossAngle(section) => section.elastic(),
            ModalBendingLoss::Thermoelastic { thermal, .. } => thermal.elastic(),
        }
    }

    /// Original selected loss tangent, either complete or explicitly excess.
    #[must_use]
    pub fn loss_angle(&self) -> Option<&ConcentricLossAngleSection> {
        match &self.law {
            ModalBendingLoss::LossAngle(section) => Some(section),
            ModalBendingLoss::Thermoelastic { excess, .. } => excess.as_ref(),
        }
    }

    /// Original circular thermoelastic law and its thermal material coefficients.
    #[must_use]
    pub fn thermoelastic(&self) -> Option<&ConcentricThermoelasticBending> {
        match &self.law {
            ModalBendingLoss::Thermoelastic { thermal, .. } => Some(thermal),
            ModalBendingLoss::LossAngle(_) => None,
        }
    }

    /// Intersection of every selected coefficient's frequency support [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }

    /// Imaginary bending stiffness H(omega) [N m2] at an admitted frequency.
    ///
    /// # Errors
    /// Nonpositive/out-of-band frequency or unrepresentable summed loss.
    pub fn hysteretic_stiffness_n_m2(
        &self, omega_rad_s: f64,
    ) -> Result<f64, MaterialStatePointError> {
        let omega = positive(omega_rad_s, "modal bending angular frequency")?;
        let (lo, hi) = self.omega_band_rad_s;
        if omega < lo || omega > hi {
            return Err(invalid("modal bending frequency is outside the original selected source support"));
        }
        let excess = self.loss_angle().map_or(0.0,
            ConcentricLossAngleSection::hysteretic_stiffness_n_m2);
        let thermal = self.thermoelastic().map_or(0.0, |law|
            law.hysteretic_scale_n_m2()
                * circular_relaxation_shape(omega, law.base_time_s()));
        let value = excess + thermal;
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("modal hysteretic bending stiffness is unrepresentable"));
        }
        Ok(value)
    }

    /// Cycle-equivalent rate 2*zeta*omega [1/s], with tension dilution retained.
    ///
    /// Use each actual mode's frequency and total dynamic line mass, including
    /// explicitly bound fluid inertia. This does not return Kelvin-Voigt
    /// viscosity or alter elastic bending stiffness/bridge shear.
    ///
    /// # Errors
    /// Unsupported frequency, nonpositive wave number/mass, or unrepresentable loss.
    pub fn damping_rate_per_s(
        &self, omega_rad_s: f64, wave_number_per_m: f64, line_density_kg_m: f64,
    ) -> Result<f64, MaterialStatePointError> {
        let wave_number = positive(wave_number_per_m, "modal bending wave number")?;
        let mass = positive(line_density_kg_m, "modal bending dynamic line density")?;
        let stiffness = self.hysteretic_stiffness_n_m2(omega_rad_s)?;
        let rate = stiffness * wave_number.powi(4) / mass / omega_rad_s;
        if !rate.is_finite() || rate < 0.0 || (stiffness > 0.0 && rate == 0.0) {
            return Err(invalid("modal bending dissipation rate is unrepresentable"));
        }
        Ok(rate)
    }
}
