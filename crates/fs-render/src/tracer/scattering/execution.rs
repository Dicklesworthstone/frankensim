//! Physical media use the existing bounded tile executor and publication path.

use super::*;

/// Render a borrowed participating medium with the ordinary tile worker pool.
///
/// The field is shared immutably by workers. Pixel/sample random identities and
/// accumulation order match [`render_scattering`] exactly, independently of
/// tile shape and worker count. The existing execution policy bounds staging
/// and tile scratch; no second image-sized buffer is allocated for a fresh film.
pub fn render_scattering_with_execution(
    scene: &Scene,
    medium: &ScatteringGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    render_medium_with_execution(
        scene,
        cx,
        settings,
        scattering_medium(medium),
        medium.has_thermal_emission(),
        execution,
    )
}

/// Append a complete scattering sample range using the bounded tile executor.
///
/// The same field/settings binding and completed-prefix rules apply as for
/// [`render_scattering_range`]. Worker failure, cancellation, exhausted cell
/// traversal, memory refusal or nonfinite radiance leaves the film unchanged.
#[allow(clippy::too_many_arguments)]
pub fn render_scattering_range_with_execution(
    scene: &Scene,
    medium: &ScatteringGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut ScatteringFilm,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    cx.checkpoint()?;
    let film = film.checked_append(medium, settings, from, to)?;
    render_medium_range_with_execution(
        scene,
        cx,
        settings,
        film,
        scattering_medium(medium),
        medium.has_thermal_emission(),
        from,
        to,
        execution,
    )
}

/// Render a borrowed absorbing/emitting field with bounded parallel execution.
///
/// Beer--Lambert attenuation and LTE emission retain the serial estimator's
/// exact pixel stream. Workers borrow the same accepted field without copying
/// its density or temperature arrays. A fresh render allocates one film payload.
pub fn render_absorbing_with_execution(
    scene: &Scene,
    medium: &AbsorbingGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    render_medium_with_execution(
        scene,
        cx,
        settings,
        (!medium.is_vacuum()).then_some(PathMedium::Absorbing(medium)),
        medium.has_thermal_emission(),
        execution,
    )
}

/// Append an absorbing sample range transactionally with the ordinary executor.
///
/// Range, field and settings checks match [`render_absorbing_range`]. All tile
/// results remain private until successful drain and the final checkpoint.
#[allow(clippy::too_many_arguments)]
pub fn render_absorbing_range_with_execution(
    scene: &Scene,
    medium: &AbsorbingGrid<'_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut AbsorbingFilm,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    cx.checkpoint()?;
    let film = film.checked_append(medium, settings, from, to)?;
    render_medium_range_with_execution(
        scene,
        cx,
        settings,
        film,
        (!medium.is_vacuum()).then_some(PathMedium::Absorbing(medium)),
        medium.has_thermal_emission(),
        from,
        to,
        execution,
    )
}

fn render_medium_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    medium: Option<PathMedium<'_, '_>>,
    has_emission: bool,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    cx.checkpoint()?;
    if settings.spp == 0 || settings.max_depth == 0 {
        return Err(TracerError::InvalidInput.into());
    }
    validate_medium_scene(scene, cx, medium)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_medium_impl(
        scene,
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
        medium,
        has_emission,
    )?;
    let publication_started = Instant::now();
    let film = Film {
        width: settings.width,
        height: settings.height,
        xyz: result.xyz.ok_or(RenderExecutionError::Internal(
            "fresh physical render omitted staging film",
        ))?,
        spp_done: settings.spp,
        time_mode: result.requested_mode,
    };
    let mut report = result.report;
    report.publication_ns = elapsed_ns(publication_started);
    Ok(RenderExecutionOutput { film, report })
}

#[allow(clippy::too_many_arguments)]
fn render_medium_range_with_execution(
    scene: &Scene,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut Film,
    medium: Option<PathMedium<'_, '_>>,
    has_emission: bool,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    validate_medium_scene(scene, cx, medium)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_medium_impl(
        scene,
        cx,
        settings,
        Some(film),
        None,
        from,
        to,
        None,
        CameraPath::Legacy,
        execution,
        &pool,
        medium,
        has_emission,
    )?;
    let mut report = result.report;
    if let Some(xyz) = result.xyz {
        let publication_started = Instant::now();
        film.xyz = xyz;
        film.spp_done = to;
        if film.time_mode == FilmTimeMode::Uninitialized {
            film.time_mode = result.requested_mode;
        }
        report.publication_ns = elapsed_ns(publication_started);
    }
    Ok(report)
}
