//! Source-bound loss-angle (hysteretic) bending loss for bonded sections.
//!
//! A layer whose complex modulus is `E_i (1 + i tan(delta_i))` over its source
//! band gives the section the complex bending stiffness `EI + i H`, with
//! `H = sum(E_i I_i tan(delta_i))`. Layers keep their own loss; nothing is
//! averaged. The claim is the loss tangent `Im(E)/Re(E)` (the loss angle phi
//! of the low-loss literature, the loss factor eta); a quality factor or a
//! logarithmic decrement is a different quantity and is never converted here.
//!
//! A modal time integrator carries viscous rates only. For a mode at angular
//! frequency omega with wavenumber k, the rate that dissipates the same energy
//! per cycle is `H k^4 / (mu omega)` (a Kelvin-Voigt coefficient `H / omega`).
//! The consumer evaluates it at each actual mode's own frequency, which must
//! lie inside the intersected source band. This is frozen-state small-strain
//! bending loss, not thermoelastic relaxation computed from thermal data.

use super::{ConcentricElasticSection, MaterialStatePointError, invalid};
use crate::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use fs_qty::{Dims, QuantitySpec};

/// Bonded section with per-layer sourced loss tangents in bending.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricLossAngleSection {
    elastic: ConcentricElasticSection,
    loss_property: String,
    layer_hysteretic_stiffness_n_m2: Vec<f64>,
    hysteretic_stiffness_n_m2: f64,
    omega_band_rad_s: (f64, f64),
}

impl ConcentricElasticSection {
    /// Bind a named dimensionless loss-tangent claim from EVERY layer.
    ///
    /// Density, E and the loss tangent are evaluated at one held material
    /// environment under the shared fixed-environment policy. Temperature,
    /// pressure, composition and moisture stay at their resolved values; no
    /// strain, stress, field or frequency response is frozen. Each tangent
    /// declares a frequency band;
    /// all coefficient bands intersect. Explicit zero is lossless, but absent
    /// loss data is an error, never a lossless layer.
    ///
    /// # Errors
    /// Missing/ambiguous/wrong-quantity/negative data, state-varying
    /// coefficients, disjoint support and unrepresentable integrals refuse.
    pub fn with_loss_angle_bending(
        self,
        loss_property: &str,
    ) -> Result<ConcentricLossAngleSection, MaterialStatePointError> {
        let mut band = (0.0_f64, f64::MAX);
        let mut total = 0.0;
        let mut layers = Vec::with_capacity(self.layers().len());
        for (layer, geometry) in self.layers().iter().zip(self.layer_integrals()) {
            let loss = layer.material.property(loss_property)
                .ok_or_else(|| invalid("every layer needs an unambiguous bending loss tangent"))?;
            let tangent = loss.value_si();
            if loss.requirement().quantity() != QuantitySpec::dimensional(Dims::NONE)
                || !tangent.is_finite() || tangent < 0.0 {
                return Err(invalid("layer loss tangent must be a nonnegative dimensionless ratio"));
            }
            let interval = loss.angular_frequency_band_rad_s()?
                .ok_or_else(|| invalid("each loss tangent needs an explicit frequency validity band"))?;
            band = (band.0.max(interval.0), band.1.min(interval.1));
            for key in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, loss_property] {
                let coefficient = layer.material.property(key)
                    .ok_or_else(|| invalid("layer bending law needs density, modulus and loss tangent"))?;
                if !coefficient.is_constant_at_fixed_environment() {
                    return Err(invalid("layer bending coefficients require the original held material environment"));
                }
                if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s()? {
                    band = (band.0.max(lo), band.1.min(hi));
                }
            }
            let value = tangent * geometry.bending_stiffness_n_m2();
            if !value.is_finite() || (tangent > 0.0 && value == 0.0) {
                return Err(invalid("layer hysteretic flexural stiffness is unrepresentable"));
            }
            total += value;
            if !total.is_finite() { return Err(invalid("summed hysteretic flexural stiffness is unrepresentable")); }
            layers.push(value);
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("layer bending frequency validity bands do not intersect"));
        }
        Ok(ConcentricLossAngleSection {
            elastic: self, loss_property: loss_property.into(),
            layer_hysteretic_stiffness_n_m2: layers, hysteretic_stiffness_n_m2: total,
            omega_band_rad_s: band,
        })
    }
}

impl ConcentricLossAngleSection {
    /// Original elastic section and its source bundles.
    #[must_use]
    pub const fn elastic(&self) -> &ConcentricElasticSection { &self.elastic }
    /// Name of the selected loss-tangent property.
    #[must_use]
    pub fn loss_property(&self) -> &str { &self.loss_property }
    /// `E_i I_i tan(delta_i)` per layer [N m2], in spatial layer order.
    #[must_use]
    pub fn layer_hysteretic_stiffness_n_m2(&self) -> &[f64] { &self.layer_hysteretic_stiffness_n_m2 }
    /// Imaginary part H of the complex bending stiffness [N m2].
    #[must_use]
    pub const fn hysteretic_stiffness_n_m2(&self) -> f64 { self.hysteretic_stiffness_n_m2 }
    /// Intersection of all coefficients' source support [rad/s].
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band_rad_s }
}
