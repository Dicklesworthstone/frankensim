//! Zener thermoelastic bending loss of a homogeneous solid circular member.
//!
//! The loss is computed, not sourced: E, density, linear thermal expansion,
//! specific heat capacity and thermal conductivity of the member's one
//! material state, and the absolute temperature of that state's query point,
//! give `Delta = E alpha^2 T / (rho c_p)` and the circular-section series of
//! [`crate::visco::circular_ted`]. The imaginary part of the bending stiffness
//! at angular frequency omega is `EI Delta S(omega)`, where S is that series'
//! shape. Intrinsic (anelastic) loss beyond thermoelasticity is not included:
//! a consumer adds a measured excess inside its source band, or reports it as
//! unknown, so this law alone is a lower bound on bending loss.
//!
//! A layered section has heat flow across its interfaces and is refused, as
//! are hollow sections, a missing or ambiguous absolute temperature, and
//! state-varying coefficients.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid, positive};
use crate::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use crate::visco::{ThermoelasticZener, circular_ted::circular_relaxation_shape};
use fs_qty::{QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};

/// A homogeneous solid circular section with its computed Zener bending loss.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricThermoelasticBending {
    elastic: ConcentricElasticSection,
    law: ThermoelasticZener,
    base_time_s: f64,
    hysteretic_scale_n_m2: f64,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricElasticSection {
    /// Bind Zener thermoelastic bending loss computed from the section's own
    /// material state.
    ///
    /// # Errors
    /// Hollow or layered sections, missing/wrong-quantity/state-varying thermal
    /// coefficients, a missing or ambiguous positive absolute temperature axis,
    /// and unrepresentable derived values refuse atomically.
    pub fn with_thermoelastic_bending(self) -> Result<ConcentricThermoelasticBending, MaterialStatePointError> {
        if self.inner_radius_m() != 0.0 {
            return Err(invalid("thermoelastic bending needs a solid circular section"));
        }
        let state = &self.layers()[0].material;
        if self.layers().iter().any(|layer| layer.material.identity() != state.identity()) {
            return Err(invalid("thermoelastic bending needs one homogeneous material; layered heat flow is not modelled"));
        }
        let mut band = (0.0_f64, f64::MAX);
        let mut coefficient = |name: &str, dims| -> Result<f64, MaterialStatePointError> {
            let property = state.property(name)
                .ok_or_else(|| invalid("thermoelastic bending needs density, modulus, expansion, heat capacity and conductivity"))?;
            if property.requirement().quantity().dims() != dims {
                return Err(invalid("thermoelastic coefficient has the wrong dimensions"));
            }
            if !property.is_constant_at_fixed_environment() {
                return Err(invalid("thermoelastic coefficients require the original held material environment"));
            }
            if let Some((lo, hi)) = property.angular_frequency_band_rad_s()? {
                band = (band.0.max(lo), band.1.min(hi));
            }
            Ok(property.value_si())
        };
        let absolute = QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static));
        let mut temperatures = state.query_point().iter()
            .filter(|(axis, _)| state.axis_quantities().get(axis) == Some(&absolute));
        let temperature = match (temperatures.next(), temperatures.next()) {
            (Some((_, t)), None) => positive(*t, "positive absolute temperature")?,
            _ => return Err(invalid("thermoelastic bending needs exactly one typed absolute-temperature axis")),
        };
        let law = ThermoelasticZener {
            e: positive(coefficient(YOUNG_MODULUS_PROPERTY, fs_qty::Pressure::DIMS)?, "positive modulus")?,
            rho: positive(coefficient(DENSITY_PROPERTY, fs_qty::Density::DIMS)?, "positive density")?,
            alpha_t: coefficient(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS)?,
            cp: positive(coefficient(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS)?, "positive heat capacity")?,
            conductivity: positive(coefficient(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS)?, "positive conductivity")?,
            t0: temperature,
        };
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("thermoelastic coefficient frequency bands do not intersect"));
        }
        let delta = law.relaxation_strength();
        if !delta.is_finite() || delta < 0.0 {
            return Err(invalid("thermoelastic relaxation strength is unrepresentable"));
        }
        let base_time_s = positive(law.circular_base_time(self.outer_radius_m()), "thermal relaxation time")?;
        let hysteretic_scale_n_m2 = self.bending_stiffness_n_m2() * delta;
        if !hysteretic_scale_n_m2.is_finite() {
            return Err(invalid("thermoelastic hysteretic stiffness is unrepresentable"));
        }
        Ok(ConcentricThermoelasticBending { elastic: self, law, base_time_s,
            hysteretic_scale_n_m2, omega_band_rad_s: band })
    }
}

impl ConcentricThermoelasticBending {
    /// Original elastic section and its single source bundle.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Resolved Zener law, including the state's absolute temperature.
    #[must_use]
    pub const fn law(&self) -> ThermoelasticZener { self.law }
    /// `rho c_p R^2 / k` [s]; mode n relaxes in this time over `xi_n^2`.
    #[must_use]
    pub const fn base_time_s(&self) -> f64 { self.base_time_s }
    /// `EI Delta` [N m2]: multiplies the series shape to give Im(EI*).
    #[must_use]
    pub const fn hysteretic_scale_n_m2(&self) -> f64 { self.hysteretic_scale_n_m2 }
    /// Intersection of all selected coefficients' source support [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
    /// Thermoelastic bending loss tangent at angular frequency `omega`.
    /// Consumers must admit their actual frequency against `omega_band_rad_s`.
    #[must_use]
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        self.law.relaxation_strength() * circular_relaxation_shape(omega, self.base_time_s)
    }
}
