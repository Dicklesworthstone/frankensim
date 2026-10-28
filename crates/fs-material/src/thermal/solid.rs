//! Phase-restricted access to the original equilibrium caloric chart.
//!
//! Register this borrowed law on a whole-solid heat/sublimation timeline.
//! The unchanged curve still supports ordinary melting simulations elsewhere.

use crate::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// The finite all-solid enthalpy interval of one immutable material chart.
///
/// Temperature and relative entropy delegate to the original curve. Admission
/// uses its authored zero-liquid-fraction knots, not temperature: an isothermal
/// melting plateau is not an all-solid inventory. The melting-onset endpoint
/// belongs to the solid branch; any larger enthalpy refuses.
///
/// Register this value as the cell's `ThermalLaw` when using the original
/// single-surface, shared-chamber, or external solid/gas exchange APIs. Their
/// initial, heat, nonlinear trial and accepted caloric queries then share this
/// phase restriction. A rejected nonlinear trial may still backtrack to a valid
/// solid endpoint. No phase fractions, latent heat or densities are invented.
///
/// This is an opt-in model-domain restriction, not mixed-phase sublimation,
/// melting competition, a vapor-law compatibility check or a topology update.
#[derive(Clone, Copy, Debug)]
pub struct PureSolidThermalLaw<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    range_j_kg: [f64; 2],
    maximum_slope: f64,
}

impl<'a> PureSolidThermalLaw<'a> {
    /// Retain the complete source curve and admit its positive-length solid span.
    /// Preparation scans only the solid prefix once; caloric queries allocate
    /// nothing and use the original logarithmic-time interpolation owner.
    pub fn new(curve: &'a EquilibriumEnthalpyPhaseCurve) -> Result<Self, PhaseStateError> {
        let knots = curve.knots();
        let count = knots.partition_point(|knot| knot.liquid_mass_fraction == 0.0);
        if count < 2 {
            return Err(PhaseStateError::InvalidCurve {
                what: "whole-solid thermal transport requires a nonempty all-solid enthalpy interval",
            });
        }
        let mut maximum_slope = 0.0_f64;
        for pair in knots[..count].windows(2) {
            let delta_t = pair[1].temperature_k - pair[0].temperature_k;
            let slope = delta_t
                / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg);
            if !slope.is_finite() || (delta_t > 0.0 && slope == 0.0) {
                return Err(PhaseStateError::UnrepresentableTemperatureDerivative);
            }
            maximum_slope = maximum_slope.max(slope);
        }
        Ok(Self {
            curve,
            range_j_kg: [knots[0].specific_enthalpy_j_kg, knots[count - 1].specific_enthalpy_j_kg],
            maximum_slope,
        })
    }

    /// Original unmodified curve, including its full source and curve identities.
    #[must_use]
    pub const fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.curve
    }

    /// Inclusive solid-only domain in the original specific-enthalpy datum.
    #[must_use]
    pub const fn enthalpy_range_j_kg(&self) -> [f64; 2] {
        self.range_j_kg
    }

    fn admit(&self, h: f64) -> Result<(), DepositionError> {
        if h.is_finite() && (self.range_j_kg[0]..=self.range_j_kg[1]).contains(&h) {
            Ok(())
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}

impl EquilibriumEnthalpyPhaseCurve {
    /// Borrow this chart's all-solid branch for the existing thermal timeline.
    pub fn pure_solid_thermal_law(&self) -> Result<PureSolidThermalLaw<'_>, PhaseStateError> {
        PureSolidThermalLaw::new(self)
    }
}

impl ThermalLaw for PureSolidThermalLaw<'_> {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.admit(h)?;
        self.curve.temperature_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self, initial: f64, final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.admit(initial)?;
        self.admit(final_h)?;
        self.curve.specific_entropy_change_j_kg_k(initial, final_h)
            .map(Some).map_err(|_| DepositionError::MaterialDomain)
    }
}
