//! Transactional film accumulation over once-prepared accepted scene states.

use super::*;

/// Progressive film bound to the complete accepted scene, camera and medium
/// history, shutter policy and render settings.
///
/// Each append keeps the same absolute sample addresses as a fresh render.
/// Changed accepted states or sample policy are refused before publication;
/// cancellation and failed paths preserve the entire previous film.
#[derive(Clone, Debug, PartialEq)]
pub struct SceneExposureFilm {
    film: Film,
    binding_identity: ContentHash,
    settings: Settings,
    time_mode: FilmTimeMode,
}

impl SceneExposureFilm {
    /// Allocate an empty film bound to all admitted exposure inputs.
    pub fn try_new(exposure: &SceneExposure<'_>, settings: &Settings) -> Result<Self, TracerError> {
        validate_settings(settings)?;
        Ok(Self {
            film: Film::try_new(settings.width, settings.height)?,
            binding_identity: exposure.binding_identity(),
            settings: *settings,
            time_mode: time_mode(exposure, settings)?,
        })
    }

    /// Read the accumulated XYZ sums and completed sample prefix.
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
        exposure: &SceneExposure<'_>,
        settings: &Settings,
        from: u32,
        to: u32,
    ) -> Result<(), TracerError> {
        if to < from {
            return Err(TracerError::InvalidRange { from, to });
        }
        if self.settings != *settings
            || self.binding_identity != exposure.binding_identity()
            || self.time_mode != time_mode(exposure, settings)?
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

fn time_mode(
    exposure: &SceneExposure<'_>,
    settings: &Settings,
) -> Result<FilmTimeMode, TracerError> {
    requested_time_mode(Some(exposure.shutter), settings.seed, exposure.camera)
}

/// Render held accepted scenes with their matching source proposals at each
/// sampled shutter time. Optional camera and medium histories use that time.
pub fn render_scene_exposure(
    exposure: &SceneExposure<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = SceneExposureFilm::try_new(exposure, settings)?;
    render_scene_exposure_range(exposure, cx, settings, &mut film, 0, settings.spp)?;
    Ok(film.into_film())
}

/// Append `[from,to)` atomically using the admitted scene history. Scene
/// selection happens once per camera path, before geometry or light sampling.
pub fn render_scene_exposure_range(
    exposure: &SceneExposure<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut SceneExposureFilm,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    cx.checkpoint()?;
    film.check(exposure, settings, from, to)?;
    render_range_with_scene_exposure_impl(
        exposure.initial_scene(),
        cx,
        settings,
        &mut film.film,
        from,
        to,
        Some(exposure.shutter),
        exposure.camera,
        #[cfg(feature = "volumes")]
        exposure.medium,
        true,
        Some(exposure),
    )
}

/// Render a complete scene exposure through the existing bounded tile pool.
/// Workers borrow the prepared history and a fresh render allocates one film.
pub fn render_scene_exposure_with_execution(
    exposure: &SceneExposure<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    cx.checkpoint()?;
    validate_settings(settings)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_scene_exposure_impl(
        exposure.initial_scene(),
        cx,
        settings,
        None,
        None,
        0,
        settings.spp,
        Some(exposure.shutter),
        exposure.camera,
        execution,
        &pool,
        #[cfg(feature = "volumes")]
        exposure.medium,
        true,
        Some(exposure),
    )?;
    let publication_started = Instant::now();
    let film = Film {
        width: settings.width,
        height: settings.height,
        xyz: result.xyz.ok_or(RenderExecutionError::Internal(
            "accepted scene exposure omitted film",
        ))?,
        spp_done: settings.spp,
        time_mode: result.requested_mode,
    };
    let mut report = result.report;
    report.publication_ns = elapsed_ns(publication_started);
    Ok(RenderExecutionOutput { film, report })
}

/// Append scene exposure samples with deterministic bounded tile execution.
/// A failed or cancelled tile never publishes a partial sample prefix.
pub fn render_scene_exposure_range_with_execution(
    exposure: &SceneExposure<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut SceneExposureFilm,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    cx.checkpoint()?;
    film.check(exposure, settings, from, to)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_scene_exposure_impl(
        exposure.initial_scene(),
        cx,
        settings,
        Some(&film.film),
        None,
        from,
        to,
        Some(exposure.shutter),
        exposure.camera,
        execution,
        &pool,
        #[cfg(feature = "volumes")]
        exposure.medium,
        true,
        Some(exposure),
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
