//! Finite-band material-card optics routed to the original native transport.

use super::*;
use fs_material::state_point::{
    SpectralConductorStatePoint, SpectralOpticalStatePoint, SpectralTabulatedDielectricStatePoint,
};

impl MaterialOpticalBinding {
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
    /// This preserves the native weak-absorption model and its `[1, 3]` phase
    /// index domain; it does not infer complex Fresnel or broadband emissivity.
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
            _ => return Err(MaterialOpticalRefreshError::SourceMismatch),
        }
        .map_err(MaterialOpticalRefreshError::Binding)?;
        Ok((binding, refreshed))
    }
}

#[cfg(test)]
mod tests;
