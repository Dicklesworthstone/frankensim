//! Native accepted-temperature phase and absorption for water and silica.

use super::{
    DielectricResponse, MaterialCard, MaterialStatePointError, QueryPoint, ReflectionError,
    ResolvedMaterialStatePoint, ThermoOpticalSource,
};
use fs_em::plane_wave::{SilicaOpticalResponse, WaterOpticalResponse};
use fs_material::state_point::SpectralDielectricAbsorptionStatePoint;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum AnalyticPhase {
    IapwsWater,
    Corning7980Silica,
}

/// An explicitly selected original analytic phase law and same-state loss.
///
/// The source retains the canonical nine vacuum wavelength/alpha samples,
/// their original selected claims and every additional requirement. Water also
/// requires bulk density; both models require actual permeability. Temperature
/// is read from the original thermal cell owner and refreshes this complete
/// snapshot once. The analytic phase law is evaluated at that same temperature
/// and water density; no reference-temperature phase claims are needed hot.
///
/// Existing monochromatic or mutually incoherent broadband heating integrates
/// the resulting coherent complex-index fields. Alpha is intensity absorption,
/// with no additional Ohmic heat or reservoir credit. The interval retains the
/// existing quasi-steady thermal split, fixed geometry and lossless exterior
/// requirements. This source supplies no transient constitutive storage,
/// scattering, fluid flow, thermal expansion or automatic phase classification.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticThermoOpticalSource {
    source: SpectralDielectricAbsorptionStatePoint,
    phase: AnalyticPhase,
    temperature_axis: Arc<str>,
    maximum_knot_span_m: f64,
}

impl AnalyticThermoOpticalSource {
    /// Explicitly select original IAPWS ordinary-water phase applicability.
    /// The supplied absorption snapshot also owns T, bulk density and scalar mu.
    pub fn try_iapws_water(
        source: SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        Self::try_new(source, AnalyticPhase::IapwsWater, temperature_axis, maximum_knot_span_m)
    }

    /// Explicitly select original CHARMS Corning 7980 phase applicability.
    /// Only the original 30--300 K phase domain is supported; the current
    /// absorption source can impose a narrower temperature or wavelength range.
    pub fn try_corning_7980_silica(
        source: SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        Self::try_new(source, AnalyticPhase::Corning7980Silica, temperature_axis, maximum_knot_span_m)
    }

    fn try_new(
        source: SpectralDielectricAbsorptionStatePoint,
        phase: AnalyticPhase,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        let value = Self {
            source,
            phase,
            temperature_axis: Arc::from(temperature_axis),
            maximum_knot_span_m,
        };
        // The same native adapter admits original source state, analytic phase,
        // permeability and the complete phase/absorption intersection.
        match phase {
            AnalyticPhase::IapwsWater => {
                WaterOpticalResponse::try_with_current_absorption(
                    &value.source, &value.temperature_axis, value.maximum_knot_span_m,
                )?;
            }
            AnalyticPhase::Corning7980Silica => {
                SilicaOpticalResponse::try_with_current_absorption(
                    &value.source, &value.temperature_axis, value.maximum_knot_span_m,
                )?;
            }
        }
        Ok(value)
    }

    /// Complete original current-state loss and its retained source identity.
    pub const fn source(&self) -> &SpectralDielectricAbsorptionStatePoint {
        &self.source
    }

    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }

    /// Retained positive finite vacuum-wavelength interpolation limit [m].
    pub const fn maximum_knot_span_m(&self) -> f64 {
        self.maximum_knot_span_m
    }
}

impl ThermoOpticalSource for AnalyticThermoOpticalSource {
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.source.resolved()
    }

    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Ok(Self {
            source: self.source.at_query_point(card, point)?,
            phase: self.phase,
            temperature_axis: Arc::clone(&self.temperature_axis),
            maximum_knot_span_m: self.maximum_knot_span_m,
        })
    }

    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError> {
        // The native thermal owner calls this after each accepted-temperature
        // source refresh and BEFORE publishing either applied or endpoint states.
        // Moving source knots and changing density/mu retain their original
        // claims; phase or spectral support exits abort the joint transaction.
        match self.phase {
            AnalyticPhase::IapwsWater => {
                let response = WaterOpticalResponse::try_with_current_absorption(
                    &self.source, &self.temperature_axis, self.maximum_knot_span_m,
                )?;
                response.permittivity_f_m(omega_rad_s)?;
                Ok(DielectricResponse::WaterOptical(response))
            }
            AnalyticPhase::Corning7980Silica => {
                let response = SilicaOpticalResponse::try_with_current_absorption(
                    &self.source, &self.temperature_axis, self.maximum_knot_span_m,
                )?;
                response.permittivity_f_m(omega_rad_s)?;
                Ok(DielectricResponse::SilicaOptical(response))
            }
        }
    }
}
