//! Static spectral rendering through an accepted absorbing/emitting grid.
//!
//! This entry uses the original surface path integrator and random streams.
//! Cell transport is deterministic Beer--Lambert/LTE integration; there is no
//! separate image-space fog filter, scattering approximation, or RGB source.

use super::*;
use crate::volumes::transport::AbsorbingGrid;

/// Progressive film bound to one exact medium snapshot and render policy.
///
/// The accumulated film is read-only outside this wrapper. Every append
/// validates the snapshot identity and all settings, including the total sample
/// ceiling. As with the ordinary render-range API, the caller must retain the
/// same scene geometry, camera, materials and lights between appends; this is
/// not a serialized scene checkpoint or a scene-content hash.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorbingFilm {
    film: Film,
    medium_identity: ContentHash,
    settings: Settings,
}

impl AbsorbingFilm {
    /// Allocate an empty film and bind its physical medium and sample policy.
    pub fn try_new(medium: &AbsorbingGrid<'_>, settings: &Settings) -> Result<Self, TracerError> {
        if settings.spp == 0 || settings.max_depth == 0 {
            return Err(TracerError::InvalidInput);
        }
        Ok(Self {
            film: Film::try_new(settings.width, settings.height)?,
            medium_identity: medium.content_identity(),
            settings: *settings,
        })
    }

    /// Read the accumulated XYZ sums and completed sample prefix.
    #[must_use]
    pub const fn film(&self) -> &Film {
        &self.film
    }

    /// Finish progressive accumulation and take the ordinary output film.
    #[must_use]
    pub fn into_film(self) -> Film {
        self.film
    }

    pub(super) fn checked_append(
        &mut self,
        medium: &AbsorbingGrid<'_>,
        settings: &Settings,
        from: u32,
        to: u32,
    ) -> Result<&mut Film, TracerError> {
        if to < from {
            return Err(TracerError::InvalidRange { from, to });
        }
        if *settings != self.settings
            || medium.content_identity() != self.medium_identity
            || self.film.spp_done != from
            || to > settings.spp
        {
            return Err(TracerError::InvalidInput);
        }
        Ok(&mut self.film)
    }
}

/// Render a static scene through a physical, non-scattering index-one medium.
///
/// Density is kg/m³, mass absorption is m²/kg, world distance is metres and
/// local-equilibrium emission is W/(m² sr nm). The same volume affects camera,
/// reflected and direct-light shadow segments. A closed dielectric displaces
/// the ambient grid and retains its own homogeneous
/// spectral absorption. Reflected/refracted continuations carry the same grid
/// after exit. Nonzero grids use complete interface random walks rather than
/// the surface-only specialized refractive NEE proposals. A vacuum grid retains
/// the ordinary tracer's exact pixel stream and results.
pub fn render_absorbing(
    scene: &Scene,
    medium: &AbsorbingGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = AbsorbingFilm::try_new(medium, settings)?;
    render_absorbing_range(scene, medium, cx, settings, &mut film, 0, settings.spp)?;
    Ok(film.into_film())
}

/// Append `[from, to)` without exposing partially traced pixels on failure.
///
/// `from` must equal the completed prefix and `to` must not exceed the bound
/// `settings.spp`. Changing any settings or the medium content identity returns
/// `InvalidInput`. Cancellation, geometry refusal, volume traversal-budget
/// exhaustion and nonfinite radiance leave every accumulated bit unchanged.
#[allow(clippy::too_many_arguments)]
pub fn render_absorbing_range(
    scene: &Scene,
    medium: &AbsorbingGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut AbsorbingFilm,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    cx.checkpoint()?;
    let film = film.checked_append(medium, settings, from, to)?;
    scattering::render_medium_range(
        scene,
        cx,
        settings,
        film,
        (!medium.is_vacuum()).then_some(PathMedium::Absorbing(medium)),
        medium.has_thermal_emission(),
        from,
        to,
    )
}
