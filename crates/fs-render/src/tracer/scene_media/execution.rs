//! Existing serial/tile kernels with immutable physical ownership and atomic publication.

use super::*;

/// Progressive film bound to the admitted scene, medium owners and settings.
#[derive(Clone, Debug, PartialEq)]
pub struct SceneMediaFilm {
    film: Film,
    binding_identity: ContentHash,
    settings: Settings,
}

impl SceneMediaFilm {
    /// Allocate an empty film for one accepted physical scene.
    pub fn try_new(media: &SceneMedia<'_, '_>, settings: &Settings) -> Result<Self, TracerError> {
        validate_settings(settings)?;
        Ok(Self {
            film: Film::try_new(settings.width, settings.height)?,
            binding_identity: media.content_identity(),
            settings: *settings,
        })
    }

    /// Read accumulated XYZ sums and completed absolute sample prefix.
    #[must_use]
    pub const fn film(&self) -> &Film {
        &self.film
    }

    /// Finish accumulation and take the ordinary render film.
    #[must_use]
    pub fn into_film(self) -> Film {
        self.film
    }

    fn check(
        &self,
        media: &SceneMedia<'_, '_>,
        settings: &Settings,
        from: u32,
        to: u32,
    ) -> Result<(), TracerError> {
        if to < from {
            return Err(TracerError::InvalidRange { from, to });
        }
        if self.settings != *settings
            || self.binding_identity != media.content_identity()
            || self.film.spp_done != from
            || to > settings.spp
        {
            return Err(TracerError::InvalidInput);
        }
        Ok(())
    }
}

fn validate_settings(settings: &Settings) -> Result<(), TracerError> {
    if settings.spp == 0 || settings.max_depth == 0 {
        return Err(TracerError::InvalidInput);
    }
    Ok(())
}

/// Render physical absorption, thermal emission and scattering in glass bodies.
///
/// All interfaces and field regions share the existing spectral path walk.
/// `max_depth` bounds surface and scattering vertices together; grids retain
/// their admitted per-segment traversal limits. The camera must start outside
/// every declared closed dielectric.
pub fn render_scene_media(
    media: &SceneMedia<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = SceneMediaFilm::try_new(media, settings)?;
    render_scene_media_range(media, cx, settings, &mut film, 0, settings.spp)?;
    Ok(film.into_film())
}

/// Append `[from,to)` atomically; failures preserve the previous complete film.
pub fn render_scene_media_range(
    media: &SceneMedia<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut SceneMediaFilm,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    cx.checkpoint()?;
    film.check(media, settings, from, to)?;
    render_range_with_scene_exposure_impl(
        media.scene,
        cx,
        settings,
        &mut film.film,
        from,
        to,
        None,
        CameraPath::Legacy,
        media.ambient_path().map(MediumInput::Static),
        media.has_thermal_emission(),
        None,
        Some(media),
    )
}

/// Render through the existing bounded native tile executor.
///
/// Every worker borrows the same scene/field owners and prepared light
/// proposals. Pixel/sample identities and summation order match serial output.
pub fn render_scene_media_with_execution(
    media: &SceneMedia<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    cx.checkpoint()?;
    validate_settings(settings)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_scene_exposure_impl(
        media.scene,
        cx,
        settings,
        None,
        None,
        0,
        settings.spp,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
        media.ambient_path().map(MediumInput::Static),
        media.has_thermal_emission(),
        None,
        Some(media),
    )?;
    let publication_started = Instant::now();
    let film = Film {
        width: settings.width,
        height: settings.height,
        xyz: result.xyz.ok_or(RenderExecutionError::Internal(
            "physical scene media omitted film",
        ))?,
        spp_done: settings.spp,
        time_mode: result.requested_mode,
    };
    let mut report = result.report;
    report.publication_ns = elapsed_ns(publication_started);
    Ok(RenderExecutionOutput { film, report })
}

/// Append a complete range using deterministic tiles and atomic publication.
pub fn render_scene_media_range_with_execution(
    media: &SceneMedia<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut SceneMediaFilm,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    cx.checkpoint()?;
    film.check(media, settings, from, to)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_scene_exposure_impl(
        media.scene,
        cx,
        settings,
        Some(&film.film),
        None,
        from,
        to,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
        media.ambient_path().map(MediumInput::Static),
        media.has_thermal_emission(),
        None,
        Some(media),
    )?;
    let mut report = result.report;
    if let Some(xyz) = result.xyz {
        let publication_started = Instant::now();
        film.film.xyz = xyz;
        film.film.spp_done = to;
        film.film.time_mode = result.requested_mode;
        report.publication_ns = elapsed_ns(publication_started);
    }
    Ok(report)
}
