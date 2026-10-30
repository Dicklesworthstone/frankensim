//! One accepted physical state for every ray in a sampled camera exposure.

use super::*;
use crate::volumes::transport::timeline::{MediumTimeline, PhysicalMedium, TimelineError};

/// Shutter, physical history and optional cinematic camera admitted together.
///
/// Every path selects one held snapshot using the same absolute time as its
/// camera and animated geometry. Density/temperature discontinuities are never
/// blended; this is quasi-static exposure integration, not light travel time.
#[derive(Clone, Copy)]
pub struct MediumExposure<'view, 'field> {
    history: &'view MediumTimeline<'view, 'field>,
    shutter: ShutterInterval,
    camera: CameraPath<'view>,
    camera_identity: Option<ContentHash>,
}

impl<'view, 'field> MediumExposure<'view, 'field> {
    /// Admit complete history coverage for a scene-camera exposure.
    pub fn try_new(
        cx: &Cx<'_>,
        history: &'view MediumTimeline<'view, 'field>,
        shutter: ShutterInterval,
    ) -> Result<Self, TracerError> {
        history.admit_shutter(cx, shutter)?;
        Ok(Self {
            history,
            shutter,
            camera: CameraPath::Legacy,
            camera_identity: None,
        })
    }

    /// Use a thin-lens/keyframed camera at the same absolute physical time.
    /// Existing camera cut and exposure-coverage rules remain binding.
    pub fn with_camera(
        mut self,
        cx: &Cx<'_>,
        camera: &'view AnimatedCamera,
        cut_side: CutSide,
    ) -> Result<Self, TracerError> {
        let exposure = camera.admit_shutter(cx, self.shutter, cut_side)?;
        self.camera_identity = Some(camera_identity(cx, camera)?);
        self.camera = CameraPath::Cinematic { camera, exposure };
        Ok(self)
    }

    /// Accepted exposure interval and its temporal sampling policy.
    #[must_use]
    pub const fn shutter(&self) -> ShutterInterval {
        self.shutter
    }

    /// Numerical identity of the accepted physical history.
    #[must_use]
    pub fn history_identity(&self) -> ContentHash {
        self.history.content_identity()
    }

    fn time_mode(&self, settings: &Settings) -> Result<FilmTimeMode, TracerError> {
        requested_time_mode(Some(self.shutter), settings.seed, self.camera)
    }
}

/// Progressive exposure bound to the physical history, cinematic camera and shutter.
/// Scene geometry/materials/lights and the legacy scene camera must remain fixed
/// between appends; animated owners are evaluated at the sampled physical times.
#[derive(Clone, Debug, PartialEq)]
pub struct MediumExposureFilm {
    film: Film,
    history_identity: ContentHash,
    camera_identity: Option<ContentHash>,
    settings: Settings,
    time_mode: FilmTimeMode,
}

impl MediumExposureFilm {
    /// Allocate an empty film bound to the complete exposure/sample policy.
    pub fn try_new(
        exposure: &MediumExposure<'_, '_>,
        settings: &Settings,
    ) -> Result<Self, TracerError> {
        validate_settings(settings)?;
        Ok(Self {
            film: Film::try_new(settings.width, settings.height)?,
            history_identity: exposure.history_identity(),
            camera_identity: exposure.camera_identity,
            settings: *settings,
            time_mode: exposure.time_mode(settings)?,
        })
    }

    /// Read accumulated XYZ sums and completed sample prefix.
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
        exposure: &MediumExposure<'_, '_>,
        settings: &Settings,
        from: u32,
        to: u32,
    ) -> Result<(), TracerError> {
        if to < from {
            return Err(TracerError::InvalidRange { from, to });
        }
        if self.settings != *settings
            || self.history_identity != exposure.history_identity()
            || self.camera_identity != exposure.camera_identity
            || self.time_mode != exposure.time_mode(settings)?
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

// Cache the numeric optics once, not just a caller-chosen shot identifier.
// Reuse the same physical-camera hashing as cinematic progressive rendering.
fn camera_identity(cx: &Cx<'_>, camera: &AnimatedCamera) -> Result<ContentHash, TracerError> {
    let mut hash = DomainHasher::new("org.frankensim.render.medium-exposure-camera.v1");
    hash.update(&(camera.shots().len() as u64).to_le_bytes());
    for shot in camera.shots() {
        cx.checkpoint()?;
        hash.update(&shot.shot_id().to_le_bytes());
        update_continuity_f64(&mut hash, shot.start_s());
        update_continuity_f64(&mut hash, shot.end_s());
        hash.update(&(shot.keyframes().len() as u64).to_le_bytes());
        for keyframe in shot.keyframes() {
            cx.checkpoint()?;
            update_continuity_f64(&mut hash, keyframe.absolute_time_s());
            update_physical_camera_continuity(&mut hash, keyframe.camera());
            match keyframe.focus() {
                KeyframeFocus::AxialDistance => hash.update(&[0]),
                KeyframeFocus::WorldPoint(point) => {
                    hash.update(&[1]);
                    update_continuity_point(&mut hash, point);
                }
            }
        }
    }
    cx.checkpoint()?;
    Ok(hash.finalize())
}

/// Render accepted physical history with the scene or admitted cinematic camera.
pub fn render_medium_exposure(
    scene: &Scene,
    exposure: &MediumExposure<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
) -> Result<Film, TracerError> {
    cx.checkpoint()?;
    let mut film = MediumExposureFilm::try_new(exposure, settings)?;
    render_medium_exposure_range(scene, exposure, cx, settings, &mut film, 0, settings.spp)?;
    Ok(film.into_film())
}

/// Append `[from,to)` atomically, refusing changed history, shutter or camera.
#[allow(clippy::too_many_arguments)]
pub fn render_medium_exposure_range(
    scene: &Scene,
    exposure: &MediumExposure<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut MediumExposureFilm,
    from: u32,
    to: u32,
) -> Result<(), TracerError> {
    cx.checkpoint()?;
    film.check(exposure, settings, from, to)?;
    exposure.history.admit_shutter(cx, exposure.shutter)?;
    render_range_with_medium_impl(
        scene,
        cx,
        settings,
        &mut film.film,
        from,
        to,
        Some(exposure.shutter),
        exposure.camera,
        Some(MediumInput::Timeline(exposure.history)),
        exposure.history.has_thermal_emission(),
    )
}

/// Render a complete physical exposure through the existing bounded tile pool.
/// Workers share immutable snapshots; a fresh render allocates one film payload.
pub fn render_medium_exposure_with_execution(
    scene: &Scene,
    exposure: &MediumExposure<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionOutput, RenderExecutionError> {
    cx.checkpoint()?;
    validate_settings(settings)?;
    exposure
        .history
        .admit_shutter(cx, exposure.shutter)
        .map_err(TracerError::from)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_medium_impl(
        scene,
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
        Some(MediumInput::Timeline(exposure.history)),
        exposure.history.has_thermal_emission(),
    )?;
    let publication_started = Instant::now();
    let film = Film {
        width: settings.width,
        height: settings.height,
        xyz: result.xyz.ok_or(RenderExecutionError::Internal(
            "physical exposure omitted film",
        ))?,
        spp_done: settings.spp,
        time_mode: result.requested_mode,
    };
    let mut report = result.report;
    report.publication_ns = elapsed_ns(publication_started);
    Ok(RenderExecutionOutput { film, report })
}

/// Append physical exposure samples with deterministic bounded tile execution.
#[allow(clippy::too_many_arguments)]
pub fn render_medium_exposure_range_with_execution(
    scene: &Scene,
    exposure: &MediumExposure<'_, '_>,
    cx: &Cx<'_>,
    settings: &Settings,
    film: &mut MediumExposureFilm,
    from: u32,
    to: u32,
    execution: &RenderExecutionConfig,
) -> Result<RenderExecutionReport, RenderExecutionError> {
    cx.checkpoint()?;
    film.check(exposure, settings, from, to)?;
    exposure
        .history
        .admit_shutter(cx, exposure.shutter)
        .map_err(TracerError::from)?;
    let pool = build_render_pool(execution, cx.mode(), settings.seed);
    let result = render_range_parallel_with_medium_impl(
        scene,
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
        Some(MediumInput::Timeline(exposure.history)),
        exposure.history.has_thermal_emission(),
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

#[derive(Clone, Copy)]
pub(super) enum MediumInput<'view, 'field> {
    Static(PathMedium<'view, 'field>),
    Timeline(&'view MediumTimeline<'view, 'field>),
}

impl<'view, 'field> MediumInput<'view, 'field> {
    pub(super) fn at_time(
        self,
        time: Option<f64>,
    ) -> Result<Option<PathMedium<'view, 'field>>, TracerError> {
        match self {
            Self::Static(medium) => Ok(Some(medium)),
            Self::Timeline(history) => {
                let medium = history.at_time(time.ok_or(TracerError::MissingRayTime)?)?;
                Ok(match medium {
                    PhysicalMedium::Absorbing(grid) => {
                        (!grid.is_vacuum()).then_some(PathMedium::Absorbing(grid))
                    }
                    PhysicalMedium::Scattering(grid) if grid.is_vacuum() => None,
                    PhysicalMedium::Scattering(grid) if grid.has_scattering() => {
                        Some(PathMedium::Scattering(grid))
                    }
                    PhysicalMedium::Scattering(grid) => {
                        Some(PathMedium::Absorbing(grid.absorption_grid()))
                    }
                })
            }
        }
    }
}

impl From<TimelineError> for TracerError {
    fn from(error: TimelineError) -> Self {
        match error {
            TimelineError::Cancelled => Self::Cancelled,
            other => Self::MediumTimeline(other),
        }
    }
}
