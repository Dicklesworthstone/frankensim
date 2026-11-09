//! Finite-band material-card optics routed to the original native transport.

use super::*;
use crate::thermal_surface::{CoatingEmissionError, CoatingThermalEmission};
use fs_material::state_point::{
    SpectralConductorStatePoint, SpectralOpticalStatePoint, SpectralSellmeierStatePoint,
    SpectralTabulatedDielectricStatePoint,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Refusal before a material state becomes a slab's physical thermal source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialSlabEmissionError {
    /// The supplied state is not the exact accepted state of this binding.
    StateMismatch,
    /// Incoherent slab emission requires a smooth, uncoated dielectric.
    Surface,
    /// An existing typed absolute-temperature coordinate is required.
    TemperatureAxis,
    /// Original bulk optics, source support or thermal emission refused.
    Emission(CoatingEmissionError),
}

impl core::fmt::Display for MaterialSlabEmissionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "material slab emission refused: {self:?}")
    }
}

impl core::error::Error for MaterialSlabEmissionError {}

impl MaterialOpticalBinding {
    /// Observe finite glass thickness using the same accepted material optics
    /// and its explicitly typed absolute-temperature coordinate. No emissivity
    /// is fitted: the original directional slab absorptance supplies Kirchhoff
    /// emission into the selected vacuum port, per wavelength and polarization.
    ///
    /// `thickness_m` and the opposite medium's original phase-index law declare
    /// the physical slab geometry. The binding must have smooth uncoated faces.
    /// The immutable result feeds native spectral observation and exchange;
    /// observing it does not debit heat. A changed material state requires the
    /// existing state refresh and a new emission binding. The original support
    /// and incoherent weak-absorption approximation remain binding.
    pub fn incoherent_slab_emission(
        &self,
        accepted: &ResolvedMaterialStatePoint,
        temperature_axis: &str,
        thickness_m: f64,
        opposite_ior: DielectricIor,
    ) -> Result<CoatingThermalEmission, MaterialSlabEmissionError> {
        if accepted.identity() != self.accepted_state_identity() {
            return Err(MaterialSlabEmissionError::StateMismatch);
        }
        let OpticalAppearance::Dielectric { glass, surface } = self.appearance else {
            return Err(MaterialSlabEmissionError::Surface);
        };
        if !surface.is_delta() || surface.coating().is_some() {
            return Err(MaterialSlabEmissionError::Surface);
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let point = accepted.query();
        if point.axis_quantities().get(temperature_axis) != Some(&absolute) {
            return Err(MaterialSlabEmissionError::TemperatureAxis);
        }
        let temperature = point
            .axes()
            .get(temperature_axis)
            .copied()
            .ok_or(MaterialSlabEmissionError::TemperatureAxis)?;
        CoatingThermalEmission::try_incoherent_slab(
            glass,
            thickness_m,
            temperature,
            opposite_ior,
            self.identity(),
        )
        .map_err(MaterialSlabEmissionError::Emission)
    }

    /// Bind an original analytic Sellmeier fit on its own source wavelength
    /// domain. Its independent intensity-absorption spectrum is retained, and
    /// only their common support is available to native transport. Exact phase
    /// and group dispersion are never replaced by a sampled index table.
    /// `None` selects a smooth surface; source state and surface remain explicit.
    pub fn try_spectral_sellmeier(
        accepted: &ResolvedMaterialStatePoint,
        optical: &SpectralSellmeierStatePoint,
        roughness_alpha: Option<f64>,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        Self::validate_state(accepted, optical.resolved(), surface_state_identity)?;
        let glass = bulk_optics::spectral_sellmeier(optical)?;
        let surface = match roughness_alpha {
            Some(alpha) => DielectricSurface::try_rough(alpha)?,
            None => DielectricSurface::SMOOTH,
        };
        let mut binding = Self::new(
            accepted,
            optical.resolved(),
            surface_state_identity,
            OpticalAppearance::Dielectric { glass, surface },
        );
        binding.optical_state_identity = optical.identity();
        Ok(binding)
    }

    /// Bind nine source-selected conductor knots on their original finite band.
    /// The shared complex-index owner already carries converted vacuum nm. The
    /// explicit positive maximum knot gap is an interpolation admission limit,
    /// not an error certificate. Zero roughness selects a smooth mirror.
    /// Accepted material and optics must share one exact card and physical state.
    pub fn try_spectral_conductor(
        accepted: &ResolvedMaterialStatePoint,
        optical: &SpectralConductorStatePoint,
        max_knot_span_nm: f64,
        roughness_alpha: f64,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        Self::validate_state(accepted, optical.resolved(), surface_state_identity)?;
        let optics = bulk_optics::spectral_conductor(optical, max_knot_span_nm)?;
        let surface = if roughness_alpha == 0.0 {
            ConductorSurface::SMOOTH
        } else {
            ConductorSurface::try_rough(roughness_alpha)?
        };
        let mut binding = Self::new(
            accepted,
            optical.resolved(),
            surface_state_identity,
            OpticalAppearance::Conductor { optics, surface },
        );
        // The source schema is part of the identity, even if an adopted bundle
        // also contains the complete legacy visible or alternate-family schema.
        binding.optical_state_identity = optical.identity();
        Ok(binding)
    }

    /// Bind original finite-band phase-index and intensity-absorption knots.
    /// This preserves the native weak-absorption model, including transparent
    /// high-index sources; it does not infer complex Fresnel or broadband emissivity.
    /// Interpolation stays between source knots. `None` selects a smooth surface.
    pub fn try_spectral_tabulated_dielectric(
        accepted: &ResolvedMaterialStatePoint,
        optical: &SpectralTabulatedDielectricStatePoint,
        max_knot_span_nm: f64,
        roughness_alpha: Option<f64>,
        surface_state_identity: ContentHash,
    ) -> Result<Self, MaterialOpticalError> {
        Self::validate_state(accepted, optical.resolved(), surface_state_identity)?;
        let glass = bulk_optics::spectral_tabulated(optical, max_knot_span_nm)?;
        let surface = match roughness_alpha {
            Some(alpha) => DielectricSurface::try_rough(alpha)?,
            None => DielectricSurface::SMOOTH,
        };
        let mut binding = Self::new(
            accepted,
            optical.resolved(),
            surface_state_identity,
            OpticalAppearance::Dielectric { glass, surface },
        );
        binding.optical_state_identity = optical.identity();
        Ok(binding)
    }

    /// Refresh a finite-band binding at the new accepted material endpoint.
    /// Exact source claims, family, interpolation-gap limit, surface finish and
    /// surface identity are retained. Every new source wavelength and value is
    /// readmitted before either result is published. The original inputs remain
    /// unchanged on refusal; no support along the intervening path is inferred.
    pub fn at_spectral_material_state(
        &self,
        optical: &SpectralOpticalStatePoint,
        card: &MaterialCard,
        accepted: &ResolvedMaterialStatePoint,
    ) -> Result<(Self, SpectralOpticalStatePoint), MaterialOpticalRefreshError> {
        if optical.identity() != self.optical_state_identity
            || accepted.card_identity() != self.material_card_identity
        {
            return Err(MaterialOpticalRefreshError::SourceMismatch);
        }
        let refreshed = optical
            .at_query_point(card, &accepted.query())
            .map_err(MaterialOpticalRefreshError::Source)?;
        let binding = match (self.appearance, &refreshed) {
            (
                OpticalAppearance::Conductor { optics, surface },
                SpectralOpticalStatePoint::Conductor(source),
            ) => Self::try_spectral_conductor(
                accepted,
                source,
                optics.max_knot_span_nm(),
                surface.roughness_alpha(),
                self.surface_state_identity,
            ),
            (
                OpticalAppearance::Dielectric { glass, surface },
                SpectralOpticalStatePoint::TabulatedDielectric(source),
            ) => {
                let DielectricIor::Tabulated(table) = glass.ior() else {
                    return Err(MaterialOpticalRefreshError::SourceMismatch);
                };
                Self::try_spectral_tabulated_dielectric(
                    accepted,
                    source,
                    table.max_knot_span_nm(),
                    surface.roughness_alpha(),
                    self.surface_state_identity,
                )
            }
            (
                OpticalAppearance::Dielectric { glass, surface },
                SpectralOpticalStatePoint::DielectricSellmeier(source),
            ) => {
                if !matches!(glass.ior(), DielectricIor::Sellmeier(_)) {
                    return Err(MaterialOpticalRefreshError::SourceMismatch);
                }
                Self::try_spectral_sellmeier(
                    accepted,
                    source,
                    surface.roughness_alpha(),
                    self.surface_state_identity,
                )
            }
            _ => return Err(MaterialOpticalRefreshError::SourceMismatch),
        }
        .map_err(MaterialOpticalRefreshError::Binding)?;
        let binding = self
            .retain_surface_finish(binding)
            .map_err(MaterialOpticalRefreshError::Binding)?;
        Ok((binding, refreshed))
    }
}

#[cfg(test)]
mod sellmeier_tests;
#[cfg(test)]
mod tests;
