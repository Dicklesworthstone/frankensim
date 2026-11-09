//! Liquid dilation from original Cp(T), engineering free lengths and density.
//!
//! Liquid applicability is an explicit physical declaration by the caller,
//! not phase-stability evidence inferred from a material name. The same
//! source card, pressure, composition, caloric datum and source receipts are
//! retained. No density or expansion chart is resampled for the fluid solver.

use super::{IsobaricFreeExpansion, invalid};
use crate::phase::SolidLiquidPhase;
use crate::thermal::transport::{
    ThermalTransportError,
    expansion::{FreeExpansionBasis, FreeExpansionDensitySample},
};
use fs_blake3::ContentHash;
use fs_thermal::enthalpy::ThermalLaw;

/// A single liquid isobar, borrowing its original caloric and density owner.
///
/// Register `caloric()` on a native thermal timeline. This object only admits
/// liquid use and evaluates its instantaneous volume response; it is not a
/// second thermal law, EOS, phase-transition model or time integrator.
#[derive(Clone, Debug)]
pub struct LiquidIsobaricExpansion<'a> {
    caloric: &'a IsobaricFreeExpansion<'a>,
    // Exact original expansion knots mapped through the original Cp integral.
    // This locates directional branches without approximating any property.
    knots: Vec<(f64, f64)>,
}

/// Original source state and signed material-volume response at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidIsobaricDilation {
    pub caloric_identity: ContentHash,
    pub state: FreeExpansionDensitySample,
    pub specific_enthalpy_j_kg: f64,
    pub specific_heat_cp_j_kg_k: f64,
    /// Material derivative Dh/Dt [W/kg], with the source isobar maintained.
    pub material_specific_enthalpy_rate_w_kg: f64,
    /// Exact directional d(1/rho)/dh [m³/J].
    pub specific_volume_derivative_m3_j: f64,
    /// rho * d(1/rho)/dh * Dh/Dt [1/s]; contraction remains negative.
    pub dilation_per_s: f64,
}

impl<'a> LiquidIsobaricExpansion<'a> {
    /// Declare that the original source path physically describes one liquid
    /// throughout its supported interval. This declaration is not new measured
    /// phase evidence; callers retain responsibility for phase applicability.
    /// Solid/coexistence mobility and fixed-principal expansion are unsupported.
    /// Material names are never parsed to guess the phase.
    pub fn try_new(
        caloric: &'a IsobaricFreeExpansion<'a>,
        phase: SolidLiquidPhase,
    ) -> Result<Self, ThermalTransportError> {
        if phase != SolidLiquidPhase::Liquid
            || caloric.density().basis() != FreeExpansionBasis::Isotropic
        {
            return Err(invalid(
                "liquid dilation requires an explicitly liquid isotropic source path",
            ));
        }
        let knots = caloric
            .density()
            .expansion()
            .temperatures_k()
            .map(|t| caloric.caloric().specific_enthalpy_j_kg(t).map(|h| (h, t)))
            .collect::<Result<Vec<_>, _>>()?;
        if knots.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
            return Err(invalid(
                "original expansion knots must resolve distinct caloric coordinates",
            ));
        }
        Ok(Self { caloric, knots })
    }

    /// Exact original density/enthalpy owner, also used by expanding conduction.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsobaricFreeExpansion<'a> {
        self.caloric
    }

    /// Evaluate the original source with the tangent appropriate to the actual
    /// material enthalpy rate. Positive Cp makes its direction the temperature
    /// direction. Cooling at an interior source knot uses the exact LEFT
    /// source slope; heating uses the right. Outward endpoint rates refuse.
    /// A zero rate uses the ordinary source tangent and returns zero dilation.
    /// Negative expansion is neither clipped nor converted to positive flow.
    ///
    /// Source knots are located in exact source enthalpy coordinates, so a
    /// rounded Cp inverse cannot select the wrong side of a genuine kink.
    /// Preparation allocates once; evaluation performs bounded binary searches
    /// and no allocation, source query, perturbation or density differencing.
    pub fn dilation_at_specific_enthalpy(
        &self,
        h: f64,
        material_specific_enthalpy_rate_w_kg: f64,
    ) -> Result<LiquidIsobaricDilation, ThermalTransportError> {
        let rate = material_specific_enthalpy_rate_w_kg;
        if !rate.is_finite() || !h.is_finite() {
            return Err(invalid(
                "finite source enthalpy and material heat rate required",
            ));
        }
        let index = self.knots.partition_point(|&(knot_h, _)| knot_h < h);
        if index == self.knots.len() || (index == 0 && h < self.knots[0].0) {
            return Err(invalid("enthalpy is outside the original liquid source"));
        }
        let at_knot = self.knots[index].0 == h;
        let temperature = if at_knot {
            self.knots[index].1
        } else {
            self.caloric
                .temperature_k(h)
                .map_err(|_| invalid("enthalpy is outside the original liquid source"))?
        };
        let expansion = self.caloric.density().expansion();
        let tangent = if at_knot && rate == 0.0 {
            expansion.at_temperature_with_tangent(temperature)?
        } else {
            let increasing = if at_knot {
                rate > 0.0
            } else {
                // A rounded inverse at either end of the containing interval
                // still uses that interval's original derivative.
                let (low, high) = (self.knots[index - 1].1, self.knots[index].1);
                if temperature < low || temperature > high {
                    return Err(invalid(
                        "inverse temperature leaves its original source interval",
                    ));
                }
                temperature == low || (temperature != high && rate >= 0.0)
            };
            expansion.at_temperature_with_directional_tangent(temperature, increasing)?
        };
        let cp = self
            .caloric
            .caloric()
            .specific_heat_cp_j_kg_k(temperature)?;
        let slope = self.caloric.specific_volume_derivative_from_tangent(
            tangent.relative_volume_derivative_per_k(),
            cp,
        )?;
        let state = self.caloric.density().at_temperature(temperature)?;
        let coefficient = state.density_kg_m3() * slope;
        let dilation = coefficient * rate;
        if !coefficient.is_finite()
            || !dilation.is_finite()
            || (slope != 0.0 && coefficient == 0.0)
            || (coefficient != 0.0 && rate != 0.0 && dilation == 0.0)
        {
            return Err(invalid("original liquid dilation is unrepresentable"));
        }
        Ok(LiquidIsobaricDilation {
            caloric_identity: self.caloric.identity(),
            state,
            specific_enthalpy_j_kg: h,
            specific_heat_cp_j_kg_k: cp,
            material_specific_enthalpy_rate_w_kg: rate,
            specific_volume_derivative_m3_j: slope,
            dilation_per_s: dilation,
        })
    }
}
