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
}

/// Render a static scene through a physical, non-scattering index-one medium.
///
/// Density is kg/m³, mass absorption is m²/kg, world distance is metres and
/// local-equilibrium emission is W/(m² sr nm). The same volume affects camera,
/// reflected and direct-light shadow segments. A nonzero medium combined with
/// any dielectric material currently refuses: curved/reflected dielectric
/// connectors need their own complete segment integration. A vacuum grid
/// retains the ordinary tracer's exact pixel stream and results.
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
    if to < from {
        return Err(TracerError::InvalidRange { from, to });
    }
    if *settings != film.settings
        || medium.content_identity() != film.medium_identity
        || film.film.spp_done != from
        || to > settings.spp
    {
        return Err(TracerError::InvalidInput);
    }
    if !medium.is_vacuum() {
        for primitive in &scene.primitives {
            cx.checkpoint()?;
            if matches!(primitive.material, Material::Dielectric { .. }) {
                return Err(TracerError::VolumeDielectricUnsupported);
            }
        }
    }
    let lighting =
        validate_scene_with_medium_emission(scene, cx, None, medium.has_thermal_emission())?;
    if from == to {
        return Ok(());
    }
    let mut staged = Vec::new();
    staged
        .try_reserve_exact(film.film.xyz.len())
        .map_err(|_| TracerError::InvalidInput)?;
    for chunk in film.film.xyz.chunks(4096) {
        cx.checkpoint()?;
        staged.extend_from_slice(chunk);
    }
    let key = [
        (settings.seed & 0xffff_ffff) as u32,
        (settings.seed >> 32) as u32,
    ];
    let sobol = pixel_sobol(settings.sampler, settings.seed);
    let path_sobol = path_sobol(settings.sampler);
    let kn = 1.0 / y_integral();
    let render_cx = cx.with_stream_seed(settings.seed);
    let medium = (!medium.is_vacuum()).then_some(medium);
    for (pixel, slot) in staged.iter_mut().enumerate() {
        let pixel = u32::try_from(pixel).map_err(|_| TracerError::InvalidInput)?;
        for sample in from..to {
            cx.checkpoint()?;
            let (jx, jy, ul) = pixel_dims(settings, sobol.as_ref(), key, pixel, sample)?;
            let xyz = trace_path_with_medium(
                scene,
                &lighting,
                &render_cx,
                settings,
                kn,
                pixel,
                sample,
                jx,
                jy,
                ul,
                None,
                CameraPath::Legacy,
                false,
                path_sobol.as_ref(),
                None,
                medium,
            )?
            .xyz;
            for channel in 0..3 {
                slot[channel] += xyz[channel];
                if !slot[channel].is_finite() {
                    return Err(TracerError::InvalidInput);
                }
            }
        }
    }
    cx.checkpoint()?;
    film.film.xyz = staged;
    film.film.spp_done = to;
    film.film.time_mode = FilmTimeMode::Static;
    Ok(())
}
