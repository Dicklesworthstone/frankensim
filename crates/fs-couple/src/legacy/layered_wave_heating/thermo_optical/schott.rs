//! Analytic catalog refraction on the native cellwise wave/heat timeline.

use super::{
    DielectricResponse, MaterialCard, MaterialStatePointError, QueryPoint, ReflectionError,
    ResolvedMaterialStatePoint, SpectralOpticalResponse, SpectralOpticalStatePoint,
    ThermoOpticalSource,
};
use fs_em::plane_wave::ThermoOpticOpticalResponse;
use fs_material::{
    state_point::SpectralDielectricAbsorptionStatePoint, thermo_optic::SchottThermoOpticLaw,
};
use std::sync::Arc;

#[derive(Debug, PartialEq)]
struct Reference {
    optical: SpectralOpticalStatePoint,
    law: SchottThermoOpticLaw,
    temperature_axis: String,
    maximum_knot_span_m: f64,
    retain_reference_absorption: bool,
}

#[derive(Clone, Debug, PartialEq)]
enum Actual {
    FixedLoss(ResolvedMaterialStatePoint),
    CurrentAbsorption {
        source: SpectralDielectricAbsorptionStatePoint,
        maximum_knot_span_m: f64,
    },
}

/// One immutable reference spectrum and the current cell's material state.
///
/// Native monochromatic and broadband thermo-optical heating refresh only the
/// selected requirements of `actual`: its permeability and every additional
/// supplied requirement remain source-pinned. The original reference spectrum
/// stays at the thermal law's reference temperature, so a T0-only optical card
/// can supply refraction at hot cells without extrapolating its optical claims.
/// The Maxwell owner applies the SCHOTT correction from that original reference
/// on every call; repeated heat/cool cycles never accumulate index increments.
///
/// Reference data/law are shared across cells and endpoint snapshots. The
/// actual state belongs to the same immutable card and carries the original
/// static absolute-temperature coordinate and semantics. The current-absorption
/// mode refreshes original alpha claims alongside permeability without querying
/// reference phase claims hot. No pressure, stress, expansion, invented
/// absorption-temperature fit or transient EM law is added.
#[derive(Clone, Debug, PartialEq)]
pub struct SchottThermoOpticalSource {
    reference: Arc<Reference>,
    actual: Actual,
}

impl SchottThermoOpticalSource {
    /// Select the original lossless reference and the actual permeability state.
    /// Nonzero reference absorption is refused rather than silently discarded.
    pub fn try_lossless(
        optical: SpectralOpticalStatePoint,
        law: SchottThermoOpticLaw,
        actual: ResolvedMaterialStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        Self::try_new(
            optical,
            law,
            Actual::FixedLoss(actual),
            temperature_axis,
            maximum_knot_span_m,
            false,
        )
    }

    /// Retain the original absorption spectrum while changing real refraction.
    ///
    /// This explicitly holds alpha(lambda) fixed at its original reference
    /// temperature. The existing Maxwell field owner computes local absorbed
    /// watts from the complete complex permittivity, including the corrected
    /// real index; no Beer-law heat or extra attenuation term is added here.
    pub fn try_with_reference_absorption(
        optical: SpectralOpticalStatePoint,
        law: SchottThermoOpticLaw,
        actual: ResolvedMaterialStatePoint,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        Self::try_new(
            optical,
            law,
            Actual::FixedLoss(actual),
            temperature_axis,
            maximum_knot_span_m,
            true,
        )
    }

    /// Refresh independently sourced alpha(T,lambda) on each native cell state.
    ///
    /// The current source contains its own original wavelength/alpha claims and
    /// same-state permeability. Reference phase remains at T0; current loss is
    /// evaluated once by the existing full complex Maxwell response. Its table
    /// has a separately retained gap admission, which is rechecked after refresh.
    pub fn try_with_current_absorption(
        optical: SpectralOpticalStatePoint,
        law: SchottThermoOpticLaw,
        absorption: SpectralDielectricAbsorptionStatePoint,
        temperature_axis: &str,
        maximum_reference_knot_span_m: f64,
        maximum_absorption_knot_span_m: f64,
    ) -> Result<Self, ReflectionError> {
        Self::try_new(
            optical,
            law,
            Actual::CurrentAbsorption {
                source: absorption,
                maximum_knot_span_m: maximum_absorption_knot_span_m,
            },
            temperature_axis,
            maximum_reference_knot_span_m,
            false,
        )
    }

    fn try_new(
        optical: SpectralOpticalStatePoint,
        law: SchottThermoOpticLaw,
        actual: Actual,
        temperature_axis: &str,
        maximum_knot_span_m: f64,
        retain_reference_absorption: bool,
    ) -> Result<Self, ReflectionError> {
        let source = Self {
            reference: Arc::new(Reference {
                optical,
                law,
                temperature_axis: temperature_axis.into(),
                maximum_knot_span_m,
                retain_reference_absorption,
            }),
            actual,
        };
        // Source identity, exact reference temperature, original interpolation
        // gaps, same-state mu, and intersected T/lambda support all remain the
        // native Maxwell owner's admission, shared with its other consumers.
        source.response()?;
        Ok(source)
    }

    pub fn reference_source(&self) -> &SpectralOpticalStatePoint {
        &self.reference.optical
    }
    pub fn law(&self) -> &SchottThermoOpticLaw {
        &self.reference.law
    }
    pub fn temperature_axis(&self) -> &str {
        &self.reference.temperature_axis
    }
    pub fn maximum_knot_span_m(&self) -> f64 {
        self.reference.maximum_knot_span_m
    }
    pub fn retains_reference_absorption(&self) -> bool {
        self.reference.retain_reference_absorption
    }

    /// Original independently sourced loss at the current accepted cell state.
    pub fn current_absorption(&self) -> Option<&SpectralDielectricAbsorptionStatePoint> {
        match &self.actual {
            Actual::CurrentAbsorption { source, .. } => Some(source),
            Actual::FixedLoss(_) => None,
        }
    }

    /// Native complex response at this accepted cell state, before selecting a
    /// frequency. The source band's original endpoint admission is retained.
    pub fn response(&self) -> Result<ThermoOpticOpticalResponse<'_>, ReflectionError> {
        let reference = SpectralOpticalResponse::try_new(
            (&self.reference.optical).into(),
            self.reference.maximum_knot_span_m,
        )?;
        if let Actual::CurrentAbsorption {
            source,
            maximum_knot_span_m,
        } = &self.actual
        {
            return ThermoOpticOpticalResponse::try_with_current_absorption(
                reference,
                &self.reference.law,
                source,
                &self.reference.temperature_axis,
                *maximum_knot_span_m,
            );
        }
        if self.reference.retain_reference_absorption {
            ThermoOpticOpticalResponse::try_with_reference_absorption(
                reference,
                &self.reference.law,
                self.resolved(),
                &self.reference.temperature_axis,
            )
        } else {
            ThermoOpticOpticalResponse::try_lossless(
                reference,
                &self.reference.law,
                self.resolved(),
                &self.reference.temperature_axis,
            )
        }
    }
}

impl ThermoOpticalSource for SchottThermoOpticalSource {
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        match &self.actual {
            Actual::FixedLoss(source) => source,
            Actual::CurrentAbsorption { source, .. } => source.resolved(),
        }
    }

    fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        Ok(Self {
            reference: Arc::clone(&self.reference),
            actual: match &self.actual {
                Actual::FixedLoss(source) => Actual::FixedLoss(source.at_query_point(card, point)?),
                Actual::CurrentAbsorption {
                    source,
                    maximum_knot_span_m,
                } => Actual::CurrentAbsorption {
                    source: source.at_query_point(card, point)?,
                    maximum_knot_span_m: *maximum_knot_span_m,
                },
            },
        })
    }

    fn wave_response(&self, omega_rad_s: f64) -> Result<DielectricResponse<'_>, ReflectionError> {
        let response = self.response()?;
        response.permittivity_f_m(omega_rad_s)?;
        Ok(DielectricResponse::ThermoOpticOptical(response))
    }
}
