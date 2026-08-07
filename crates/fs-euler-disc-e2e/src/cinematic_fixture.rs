//! Direct, deterministic Euler-disc critique-clip producer.
//!
//! This is a product-facing preview path: it runs the source-bound reduced mechanics,
//! renders the resulting trajectory, derives dry sound from mechanics control
//! channels, and optionally muxes the image sequence and WAV with `ffmpeg`.
//! It is deliberately labeled as an analytical, timestep-refined trajectory
//! with uncalibrated acoustics rather than experimental ground truth.

use core::{fmt, num::NonZeroUsize};
use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};

use crate::modal_synthesis::{EulerModalParameterSet, EulerModalParameterSetInput};
use crate::{
    AUDIO_EXCITATION_ALGORITHM_VERSION, AUDIO_RECONSTRUCTION_FILTER_VERSION,
    AUDIO_RESAMPLING_ALGORITHM_VERSION, AudioArtifactBudget, AudioDryMixSpec,
    AudioEventFractionalDelay, AudioExcitationBudget, AudioExcitationMapper,
    AudioExcitationModelInput, AudioExcitationReduction, AudioMasterSource,
    AudioReconstructionFilterSpec, AudioResampler, AudioResamplingBoundaryPolicy,
    AudioResamplingBudget, AudioResamplingModelInput, ContactModeShape, ContactParticipationPolicy,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerControlStream, EulerRenderTrajectoryArtifact,
    ListenerPose as SpatialListenerPose, ListenerPoseTrack, MicrophoneDirectivity,
    ModalPresetAuthority, ModalSynthesisBudget, ModalSynthesisModelInput,
    ModeContactParticipationRule, OfflineSpatializer, RenderTrajectory, RenderTrajectoryCodecBudget,
    RepresentativeDiscMaterial, SoundWavArtifact, SourcePositionTrack, SpatialAudioAuthority,
    SpatialAudioBudget, SpatialAudioConfig, SpatialAudioOutput, SpatialAudioRenderInput,
    SpatialAudioSource, SpatialDelayPolicy, SpatialMonoSignal, SpatialOutputHorizon,
    SpatialStemComponent, StemGainPan, StereoSample, WavMetadata, WavSampleEncoding, measure_audio,
    mix_dry_modal_stems,
    reduced_decay::{
        ReducedDecayRun, RefinementEvidence, Thorne2026SteelGlassBenchmark,
        thorne_2026_refinement_evidence,
    },
    render_scene_bridge::{
        EulerCinematicScene, EulerEnvironmentStyle, EulerFrameRequest, EulerMaterialStyle,
        EulerRectLightSpec, EulerSceneConfig, EulerStudioEnvironmentSpec, EulerTessellationConfig,
        euler_scene_smoke_settings,
    },
    representative_modal_preset,
    timeline_resampling::ExposureEventPolicy,
};
use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_evidence::{
    cinematic::{CinematicClock, CinematicClockDomain, SoundAuthority},
    cinematic_config::{CinematicComponentRef, CinematicComponentRole},
    cinematic_sound::{
        ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_SYNTHESIS_SCHEMA_VERSION,
        SoundAmplitudeReference, SoundChannelLayout, SoundExcitationChannel,
        SoundExcitationControl, SoundModalComponent, SoundModelAssumption, SoundRoomResponse,
        SoundSynthesisConfig, SoundSynthesisInput, SoundTerminalPolicy, SoundTrajectoryDisposition,
    },
};
use fs_exec::{Cx, RunId};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_img::{
    CinematicColorConfig, CinematicColorLimits, PngColor, PreviewDither,
    TEMPORAL_DENOISE_PIPELINE_VERSION, TemporalDenoiseConfig, TemporalDenoiseInput,
    TemporalDenoiseLimits, TemporalDenoisedFrame, TemporalFrameBoundary, temporal_denoise_rgb,
    transform_cinematic_preview, write_png16,
};
use fs_math::det;
use fs_mbd::Vec3;
use fs_render::{
    aov::{CinematicAovConfig, CinematicAovLimits, CinematicAovProfile, CinematicAovProvenance},
    camera::{AnimatedCamera, Aperture, CameraProjection, CutSide, PhysicalCamera},
    conductor::{ConductorOptics, ConductorSurface},
    motion::{ShutterConvention, ShutterDistribution},
    tracer::{
        MAX_RENDER_TILE_EDGE, MAX_RENDER_WORKERS, RenderExecutionConfig, RenderExecutionReport,
        RenderWorkerPool, film_to_exr,
    },
};

/// Fixed master frame rate used by the cinematic sound contract.
pub const CRITIQUE_FPS: u32 = 24;
/// Minimum admitted cinematic duration: 192 frames = 8 seconds.
pub const CRITIQUE_FRAMES: u32 = 192;
/// Five-millisecond deterministic taper applied at the censored soundtrack end.
const TERMINAL_FADE_SAMPLE_FRAMES: u32 = 240;
/// Twenty-millisecond presentation fade suppressing the analytical crop boundary.
const INITIAL_FADE_SAMPLE_FRAMES: u32 = 960;
/// Conservative reconstruction ceiling for the sub-100 Hz trajectory-derived
/// harmonic-one contact modulation and its slowly varying rolling envelope.
const CRITIQUE_DECLARED_AUDIO_SOURCE_BANDWIDTH_HZ: f64 = 256.0;
/// Stable caller-ledgered root for animation frame render jobs.
const CRITIQUE_RENDER_RUN: RunId = RunId(0x4555_4c45_5252_454e);
/// Stable placement-only seed for the reusable parked render crew.
const CRITIQUE_RENDER_SCHEDULER_SEED: u64 = 0x5354_5544_494f_5631;
/// Bit-affecting version of the absolute-frame render-seed schedule.
const CRITIQUE_FRAME_SEED_SCHEDULE_VERSION: u16 = 1;
/// World-space camera/listener origin shared by picture and spatial sound.
const CRITIQUE_CAMERA_EYE_M: [f64; 3] = [0.24, -0.30, 0.18];
/// World-space look target shared by picture and spatial sound.
const CRITIQUE_CAMERA_TARGET_M: [f64; 3] = [0.0, 0.0, 0.025];
/// Artistic near-field reference distance for the desk-scale stereo preview.
///
/// The spatializer's distance gain is `reference / max(distance, reference)`.
/// Keeping this just below the camera-to-subject distance preserves a modest
/// distance cue without needlessly spending roughly 18 dB of the bounded
/// mechanics-to-digital mastering range on a five-centimetre reference.
const CRITIQUE_SPATIAL_REFERENCE_DISTANCE_M: f64 = 0.4;

/// Absolute frame selection for an affordable still or contiguous lookdev shot.
///
/// Mechanics and audio retain the complete cinematic horizon. A partial window
/// only limits expensive image rendering and is never eligible for muxing as a
/// complete movie.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CinematicFrameWindow {
    /// Render every frame in the configured cinematic horizon.
    #[default]
    Full,
    /// Render `frame_count` frames beginning at absolute `first_frame`.
    Range {
        /// Zero-based absolute frame index.
        first_frame: u32,
        /// Number of contiguous frames to render.
        frame_count: u32,
    },
}

/// Bounded settings for one watchable critique artifact.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicFixtureConfig {
    /// Display and raw-master width in pixels.
    pub width: u32,
    /// Display and raw-master height in pixels.
    pub height: u32,
    /// Exact 24 Hz frame count. This source-bound fixture is fixed to 192.
    pub frames: u32,
    /// Complete sequence or an absolute contiguous render subset.
    pub frame_window: CinematicFrameWindow,
    /// Uniform path-tracing samples per pixel.
    pub samples_per_pixel: u32,
    /// Caller-selected independent scramble salt for replicated raw renders.
    pub render_seed_salt: u64,
    /// Maximum path depth, including dielectric traversal.
    pub max_depth: u32,
    /// Physical exposure angle at 24 Hz. `180` means a 1/48-second exposure.
    pub shutter_angle_degrees: u16,
    /// Reused tile-render worker count. This affects throughput, not image bits.
    pub render_workers: usize,
    /// Logical render-tile width in pixels.
    pub tile_width: u32,
    /// Logical render-tile height in pixels.
    pub tile_height: u32,
    /// Per-frame renderer operation-memory ceiling.
    pub render_memory_limit_bytes: u64,
    /// Whether previews use the explicitly biased animation-aware denoiser.
    pub denoise_previews: bool,
    /// Persist the full DailyCore AOV EXR rather than three-channel raw beauty.
    /// Rendering retains DailyCore in memory whenever denoising needs guides.
    pub retain_full_aov_exr: bool,
    /// Whether mechanics-derived dry stems use the bounded spatial-audio path.
    pub spatialize_audio: bool,
    /// Whether to ask `ffmpeg` for a non-authoritative convenience movie.
    pub mux_with_ffmpeg: bool,
    /// `ffmpeg` executable name or path.
    pub ffmpeg_executable: PathBuf,
}

impl Default for CinematicFixtureConfig {
    fn default() -> Self {
        Self {
            width: 320,
            height: 180,
            frames: CRITIQUE_FRAMES,
            frame_window: CinematicFrameWindow::Full,
            samples_per_pixel: 1,
            render_seed_salt: 0,
            max_depth: 6,
            shutter_angle_degrees: 180,
            render_workers: default_render_workers(),
            tile_width: 32,
            tile_height: 32,
            render_memory_limit_bytes: 4 * 1024 * 1024 * 1024,
            denoise_previews: true,
            retain_full_aov_exr: true,
            spatialize_audio: true,
            mux_with_ffmpeg: true,
            ffmpeg_executable: PathBuf::from("ffmpeg"),
        }
    }
}

impl CinematicFixtureConfig {
    /// Check raster, cinematic-clock, tracer, and mux admission bounds.
    pub fn validate(&self) -> Result<(), CinematicFixtureError> {
        if self.width == 0 || self.height == 0 || self.width > 3_840 || self.height > 2_160 {
            return Err(CinematicFixtureError::InvalidConfig(
                "dimensions must be nonzero and no larger than 3840x2160",
            ));
        }
        if self.frames != CRITIQUE_FRAMES {
            return Err(CinematicFixtureError::InvalidConfig(
                "the source-bound eight-second fixture requires exactly 192 frames",
            ));
        }
        let range = self.render_frame_range()?;
        if self.mux_with_ffmpeg && range.len() != self.frames as usize {
            return Err(CinematicFixtureError::InvalidConfig(
                "partial frame windows cannot be muxed as complete movies",
            ));
        }
        if self.samples_per_pixel == 0 || self.samples_per_pixel > 4_096 {
            return Err(CinematicFixtureError::InvalidConfig(
                "samples_per_pixel must be in 1..=4096",
            ));
        }
        if self.max_depth == 0 || self.max_depth > 64 {
            return Err(CinematicFixtureError::InvalidConfig(
                "max_depth must be in 1..=64",
            ));
        }
        if self.shutter_angle_degrees > 360 {
            return Err(CinematicFixtureError::InvalidConfig(
                "shutter_angle_degrees must be in 0..=360",
            ));
        }
        if self.render_workers == 0 || self.render_workers > MAX_RENDER_WORKERS {
            return Err(CinematicFixtureError::InvalidConfig(
                "render_workers must be in 1..=256",
            ));
        }
        if self.tile_width == 0
            || self.tile_height == 0
            || self.tile_width > MAX_RENDER_TILE_EDGE
            || self.tile_height > MAX_RENDER_TILE_EDGE
        {
            return Err(CinematicFixtureError::InvalidConfig(
                "tile dimensions must be in 1..=4096",
            ));
        }
        if self.render_memory_limit_bytes == 0 {
            return Err(CinematicFixtureError::InvalidConfig(
                "render_memory_limit_bytes must be nonzero",
            ));
        }
        if self.mux_with_ffmpeg && (self.width % 2 != 0 || self.height % 2 != 0) {
            return Err(CinematicFixtureError::InvalidConfig(
                "muxed 4:2:0 video requires even width and height",
            ));
        }
        Ok(())
    }

    fn render_frame_range(&self) -> Result<core::ops::Range<u32>, CinematicFixtureError> {
        match self.frame_window {
            CinematicFrameWindow::Full => Ok(0..self.frames),
            CinematicFrameWindow::Range {
                first_frame,
                frame_count,
            } => {
                if frame_count == 0 {
                    return Err(CinematicFixtureError::InvalidConfig(
                        "frame window count must be nonzero",
                    ));
                }
                let end = first_frame.checked_add(frame_count).ok_or(
                    CinematicFixtureError::InvalidConfig("frame window end overflows u32"),
                )?;
                if first_frame >= self.frames || end > self.frames {
                    return Err(CinematicFixtureError::InvalidConfig(
                        "frame window must lie inside the configured sequence",
                    ));
                }
                Ok(first_frame..end)
            }
        }
    }
}

fn default_render_workers() -> usize {
    std::thread::available_parallelism()
        .map(NonZeroUsize::get)
        .unwrap_or(1)
        .min(MAX_RENDER_WORKERS)
}

/// Successful fixture paths and the optional convenience movie.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicFixtureReport {
    /// Atomically published output root.
    pub output_directory: PathBuf,
    /// Deterministic top-level critique manifest.
    pub manifest_path: PathBuf,
    /// Verified float32 stereo master.
    pub wav_path: PathBuf,
    /// First display-referred frame, useful for quick inspection.
    pub first_preview_path: PathBuf,
    /// Convenience movie when muxing was requested and succeeded.
    pub movie_path: Option<PathBuf>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum MuxOutcome {
    Disabled,
    Unavailable(String),
    Failed(i32),
    Written(PathBuf),
}

struct FixtureAudio {
    artifact: SoundWavArtifact,
    modal_parameter_set_identity: ContentHash,
    modal_parameter_set_disclosure: String,
    chirp_start_hz: f64,
    chirp_end_hz: f64,
    pre_master_peak_fs: f64,
    master_gain_db: f64,
    spatialization: Option<FixtureSpatialAudioEvidence>,
}

#[derive(Clone, Debug, PartialEq)]
struct FixtureSpatialAudioEvidence {
    config_identity: ContentHash,
    input_identity: ContentHash,
    raw_output_identity: ContentHash,
    mastered_output_identity: ContentHash,
    authority: SpatialAudioAuthority,
    maximum_distance_m: f64,
    maximum_delay_frames: f64,
    discarded_tail_frames: u64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct FixtureRenderEvidence {
    frames: u32,
    maximum_effective_workers: usize,
    maximum_peak_memory_bytes: u64,
    setup_ns: u128,
    traversal_ns: u128,
    tile_compute_ns: u128,
    tile_merge_ns: u128,
    publication_ns: u128,
    idle_worker_ns: u128,
}

impl FixtureRenderEvidence {
    fn observe(&mut self, report: &RenderExecutionReport) {
        self.frames += 1;
        self.maximum_effective_workers = self.maximum_effective_workers.max(report.workers);
        self.maximum_peak_memory_bytes =
            self.maximum_peak_memory_bytes.max(report.memory.peak_bytes);
        self.setup_ns += u128::from(report.setup_ns);
        self.traversal_ns += u128::from(report.traversal_ns);
        self.tile_compute_ns += u128::from(report.tile_compute_ns);
        self.tile_merge_ns += u128::from(report.tile_merge_ns);
        self.publication_ns += u128::from(report.publication_ns);
        self.idle_worker_ns += u128::from(report.idle_worker_ns);
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct FixtureDenoiseEvidence {
    applied_frames: u32,
    maximum_retained_bytes: u64,
    maximum_history_frames: u16,
}

impl FixtureDenoiseEvidence {
    fn observe(&mut self, frame: &TemporalDenoisedFrame) {
        self.applied_frames += 1;
        self.maximum_retained_bytes = self.maximum_retained_bytes.max(frame.retained_bytes());
        self.maximum_history_frames = self
            .maximum_history_frames
            .max(frame.history_length().iter().copied().max().unwrap_or(0));
    }
}

/// Fail-closed fixture refusal; a failed run never publishes the requested root.
#[derive(Debug)]
pub enum CinematicFixtureError {
    /// A bounded public configuration rule was violated.
    InvalidConfig(&'static str),
    /// A requested final or incomplete staging path already exists.
    OutputAlreadyExists(PathBuf),
    /// Filesystem publication failed.
    Io(std::io::Error),
    /// A typed mechanics, rendering, color, or sound stage refused.
    Pipeline(String),
    /// The caller cancelled the fixture at a bounded checkpoint.
    Cancelled,
}

impl fmt::Display for CinematicFixtureError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig(message) => write!(formatter, "invalid fixture config: {message}"),
            Self::OutputAlreadyExists(path) => {
                write!(
                    formatter,
                    "output directory already exists: {}",
                    path.display()
                )
            }
            Self::Io(error) => write!(formatter, "fixture I/O failed: {error}"),
            Self::Pipeline(error) => write!(formatter, "fixture pipeline failed: {error}"),
            Self::Cancelled => formatter.write_str("cinematic fixture cancelled"),
        }
    }
}

impl std::error::Error for CinematicFixtureError {}

impl From<std::io::Error> for CinematicFixtureError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

/// Produce one complete critique clip. The concrete staged implementation is
/// kept here rather than hidden behind a workflow engine.
pub fn run_cinematic_fixture(
    config: &CinematicFixtureConfig,
    output_directory: &Path,
    cx: &Cx<'_>,
    mut progress: impl FnMut(&str),
) -> Result<CinematicFixtureReport, CinematicFixtureError> {
    config.validate()?;
    let render_frame_range = config.render_frame_range()?;
    if output_directory.exists() {
        return Err(CinematicFixtureError::OutputAlreadyExists(
            output_directory.to_path_buf(),
        ));
    }
    let output_parent = output_directory
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_parent)?;
    let output_name = output_directory
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or(CinematicFixtureError::InvalidConfig(
            "output directory must have a UTF-8 final component",
        ))?;
    let staging_directory =
        output_parent.join(format!(".{output_name}.incomplete-{}", std::process::id()));
    if staging_directory.exists() {
        return Err(CinematicFixtureError::OutputAlreadyExists(
            staging_directory,
        ));
    }
    let duration_s = f64::from(config.frames) / f64::from(CRITIQUE_FPS);

    progress("stage=mechanics begin");
    let benchmark = Thorne2026SteelGlassBenchmark::ambient().map_err(pipeline)?;
    let profile = benchmark.resolve_specimen(cx).map_err(pipeline)?;
    // Reuse the coarse member of the dt/dt/2 comparison as the picture and
    // sound source. This avoids executing the identical coarse run twice while
    // ensuring every published clip carries current encoded-model refinement
    // evidence for both source-bound loss channels.
    let refinement = thorne_2026_refinement_evidence(&benchmark).map_err(pipeline)?;
    let run = &refinement.coarse;
    let trajectory =
        RenderTrajectory::from_reduced_decay_run(run, &profile, cx).map_err(pipeline)?;
    let retained_time_s = trajectory
        .samples()
        .last()
        .expect("admitted trajectory is nonempty")
        .input()
        .time_s;
    if (retained_time_s - duration_s).abs() > 16.0 * f64::EPSILON * duration_s {
        return Err(CinematicFixtureError::Pipeline(format!(
            "source-bound render tail is {retained_time_s:.17e}s, expected {duration_s:.17e}s"
        )));
    }
    progress("stage=mechanics complete");

    let mut campaign_hasher =
        DomainHasher::new("org.frankensim.euler-critique.literature-campaign.v2");
    campaign_hasher.update(trajectory.metadata().configuration_identity.as_bytes());
    campaign_hasher.update(&config.frames.to_le_bytes());
    campaign_hasher.update(&CRITIQUE_FPS.to_le_bytes());
    let campaign_identity = campaign_hasher.finalize();
    let trajectory_artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
        campaign_identity,
        trajectory,
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .map_err(pipeline)?;
    progress(&trajectory_diagnostics(&trajectory_artifact));

    fs::create_dir(&staging_directory)?;
    let trajectory_directory = staging_directory.join("trajectory");
    let raw_directory = staging_directory.join("raw");
    let preview_directory = staging_directory.join("preview");
    let sound_directory = staging_directory.join("sound");
    for directory in [
        &trajectory_directory,
        &raw_directory,
        &preview_directory,
        &sound_directory,
    ] {
        fs::create_dir(directory)?;
    }
    let trajectory_path = trajectory_directory.join("euler-trajectory.fset");
    let mut trajectory_file = create_new_file(&trajectory_path)?;
    let trajectory_receipt = trajectory_artifact
        .write_to(
            &mut trajectory_file,
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .map_err(pipeline)?;
    trajectory_file.flush()?;

    progress("stage=render begin");
    let camera = critique_camera(duration_s).map_err(pipeline)?;
    let mut scene_config = EulerSceneConfig::reference(camera);
    scene_config.tessellation = EulerTessellationConfig {
        azimuthal_segments: 64,
        arc_subdivisions_per_arc: 8,
    };
    scene_config.show_spin_fiducial = true;
    scene_config.disc_material = EulerMaterialStyle::Conductor {
        optics: ConductorOptics::representative_stainless_steel(),
        surface: ConductorSurface::try_rough(0.12).map_err(pipeline)?,
    };
    // Keep the emitter above the camera frustum while retaining a broad,
    // downward-facing studio source. The reference light otherwise appears as
    // a distracting white bar across the top of this particular composition.
    scene_config.light = EulerRectLightSpec {
        corner_world_m: Point3::new(-0.09, 0.06, 0.24),
        edge_u_world_m: GeomVec3::new(0.18, 0.0, 0.0),
        edge_v_world_m: GeomVec3::new(0.0, -0.12, 0.0),
        linear_rgb: [1.0, 0.96, 0.90],
        radiance_scale: 48.0,
    };
    scene_config.environment =
        EulerEnvironmentStyle::StudioGradient(EulerStudioEnvironmentSpec::SOFT_NEUTRAL);
    let scene = EulerCinematicScene::try_build(&trajectory_artifact, &profile, scene_config, cx)
        .map_err(pipeline)?;
    let mut base_render_settings = euler_scene_smoke_settings(config.width, config.height);
    base_render_settings.spp = config.samples_per_pixel;
    base_render_settings.max_depth = config.max_depth;
    let composition_identity = composition_identity(config, scene.scene_identity());
    let mut raw_sequence = DomainHasher::new("org.frankensim.euler-critique.raw-sequence.v1");
    let mut preview_sequence =
        DomainHasher::new("org.frankensim.euler-critique.preview-sequence.v1");
    let mut over_range_channels = 0_u64;
    let mut gamut_mapped_pixels = 0_u64;
    let mut render_evidence = FixtureRenderEvidence::default();
    let mut denoise_evidence = FixtureDenoiseEvidence::default();
    let mut denoise_history: Option<TemporalDenoisedFrame> = None;
    let denoise_config = TemporalDenoiseConfig::default();
    let aov_profile = if config.denoise_previews || config.retain_full_aov_exr {
        CinematicAovProfile::DailyCore
    } else {
        CinematicAovProfile::BeautyOnly
    };
    let raster_width = usize::try_from(config.width)
        .map_err(|_| CinematicFixtureError::Pipeline("raster width exceeds usize".into()))?;
    let raster_height = usize::try_from(config.height)
        .map_err(|_| CinematicFixtureError::Pipeline("raster height exceeds usize".into()))?;
    let exposure_duration_s =
        f64::from(config.shutter_angle_degrees) / 360.0 / f64::from(CRITIQUE_FPS);
    let pool_execution = RenderExecutionConfig::try_new(
        config.tile_width,
        config.tile_height,
        config.render_workers,
        config.render_memory_limit_bytes,
        CRITIQUE_RENDER_RUN,
    )
    .map_err(pipeline)?;
    let render_pool =
        RenderWorkerPool::new(&pool_execution, cx.mode(), CRITIQUE_RENDER_SCHEDULER_SEED);
    let trajectory_samples = trajectory_artifact.trajectory().samples();
    let trajectory_start_s = trajectory_samples[0].input().time_s;
    let trajectory_end_s = trajectory_samples[trajectory_samples.len() - 1]
        .input()
        .time_s;
    render_pool.with_parked_crew_local(|renderer| -> Result<(), CinematicFixtureError> {
        for frame in render_frame_range.clone() {
            let mut render_settings = base_render_settings;
            render_settings.seed =
                frame_render_seed(base_render_settings.seed, config.render_seed_salt, frame);
            let (frame_time_s, previous_time_s, next_time_s) =
                frame_reference_times(frame, config.frames, trajectory_start_s, trajectory_end_s);
            let prepared = scene
                .prepare_frame(EulerFrameRequest {
                    frame_time_s,
                    exposure_duration_s,
                    // Treat each image as the exposure ending at its video
                    // frame boundary. The final shutter therefore closes on
                    // the analytical validity cutoff at exactly 8 s without
                    // retiming or extrapolating the mechanics.
                    convention: ShutterConvention::BackLoaded,
                    distribution: ShutterDistribution::StratifiedCounterV1 {
                        strata: config.samples_per_pixel,
                    },
                    event_policy: ExposureEventPolicy::Refuse,
                    cut_side: CutSide::After,
                })
                .map_err(pipeline)?;
            if prepared.segments().len() != 1 {
                return Err(CinematicFixtureError::Pipeline(format!(
                    "frame {frame} exposure unexpectedly resolved to {} segments",
                    prepared.segments().len()
                )));
            }
            let provenance = CinematicAovProvenance::try_new(
                u64::from(frame),
                frame_time_s,
                previous_time_s,
                next_time_s,
                trajectory_receipt.artifact_identity(),
                scene.scene_identity(),
                composition_identity,
            )
            .map_err(pipeline)?;
            let frame_execution = RenderExecutionConfig::try_new(
                config.tile_width,
                config.tile_height,
                config.render_workers,
                config.render_memory_limit_bytes,
                CRITIQUE_RENDER_RUN.derive(
                    "org.frankensim.euler-critique.render-frame.v1",
                    u64::from(frame),
                ),
            )
            .map_err(pipeline)?;
            let output = renderer
                .render_cinematic_with_aovs(
                    scene.scene(),
                    scene.camera(),
                    prepared.cut_side(),
                    cx,
                    &render_settings,
                    prepared.segments()[0].shutter(),
                    CinematicAovConfig::new(aov_profile, provenance, CinematicAovLimits::default()),
                    &frame_execution,
                )
                .map_err(pipeline)?;
            progress(&format!(
                concat!(
                    "stage=render frame={}/{} workers={} tiles={} ",
                    "traversal_ms={:.3} compute_ms={:.3} merge_ms={:.3} peak_mib={:.3}"
                ),
                frame - render_frame_range.start + 1,
                render_frame_range.len(),
                output.report.workers,
                output.report.layout.tile_count(),
                output.report.traversal_ns as f64 / 1.0e6,
                output.report.tile_compute_ns as f64 / 1.0e6,
                output.report.tile_merge_ns as f64 / 1.0e6,
                output.report.memory.peak_bytes as f64 / (1024.0 * 1024.0),
            ));
            render_evidence.observe(&output.report);
            let film = output.film;
            let exr = if config.retain_full_aov_exr {
                film.to_exr().map_err(pipeline)?
            } else {
                film_to_exr(film.beauty()).map_err(pipeline)?
            };
            raw_sequence.update(hash_domain("frame", &exr).as_bytes());
            write_new(&raw_directory.join(format!("frame-{frame:06}.exr")), &exr)?;
            drop(exr);

            let [red, green, blue] = film.beauty().to_linear_srgb();
            let mut color = CinematicColorConfig::reference_srgb_16();
            color.exposure_ev = 1;
            color.dither = PreviewDither::Disabled;
            let preview = if config.denoise_previews {
                let guides = film.denoise_guides().map_err(pipeline)?;
                let denoise_boundary = if denoise_history.is_none() && frame != 0 {
                    TemporalFrameBoundary::Cut
                } else {
                    TemporalFrameBoundary::Continuous
                };
                let denoised = temporal_denoise_rgb(
                    TemporalDenoiseInput {
                        frame_index: u64::from(frame),
                        width: raster_width,
                        height: raster_height,
                        red: &red,
                        green: &green,
                        blue: &blue,
                        motion_prev_x: guides.motion_prev_x(),
                        motion_prev_y: guides.motion_prev_y(),
                        axial_depth_m: guides.axial_depth_m(),
                        normal_x: guides.normal_x(),
                        normal_y: guides.normal_y(),
                        normal_z: guides.normal_z(),
                        primary_coverage: guides.primary_coverage(),
                        variance_luminance: guides.variance_luminance(),
                        object_ids: None,
                        material_ids: None,
                    },
                    denoise_history.as_ref(),
                    denoise_boundary,
                    denoise_config,
                    TemporalDenoiseLimits::reference_4k(),
                )
                .map_err(pipeline)?;
                let [denoised_red, denoised_green, denoised_blue] = denoised.linear_rgb();
                let preview = transform_cinematic_preview(
                    config.width,
                    config.height,
                    [denoised_red, denoised_green, denoised_blue],
                    color,
                    CinematicColorLimits::reference_4k(),
                )
                .map_err(pipeline)?;
                denoise_evidence.observe(&denoised);
                denoise_history = Some(denoised);
                preview
            } else {
                transform_cinematic_preview(
                    config.width,
                    config.height,
                    [&red, &green, &blue],
                    color,
                    CinematicColorLimits::reference_4k(),
                )
                .map_err(pipeline)?
            };
            over_range_channels += preview.metadata().over_range_linear_channels();
            gamut_mapped_pixels += preview.metadata().gamut_mapped_pixels();
            let samples = preview.samples().as_u16().ok_or_else(|| {
                CinematicFixtureError::Pipeline(
                    "16-bit color pipeline returned 8-bit samples".into(),
                )
            })?;
            let png = write_png16(config.width, config.height, PngColor::Rgb, samples)
                .map_err(pipeline)?;
            preview_sequence.update(hash_domain("frame", &png).as_bytes());
            write_new(
                &preview_directory.join(format!("frame-{frame:06}.png")),
                &png,
            )?;
        }
        Ok(())
    })?;
    let raw_sequence_identity = raw_sequence.finalize();
    let preview_sequence_identity = preview_sequence.finalize();
    progress("stage=render complete");

    progress("stage=audio begin");
    let audio = build_audio(&trajectory_artifact, config, cx)?;
    let wav_path = sound_directory.join("master.float32.wav");
    write_new(&wav_path, audio.artifact.wav_bytes())?;
    let audio_manifest_path = sound_directory.join("master.manifest.json");
    write_new(
        &audio_manifest_path,
        audio.artifact.manifest().to_manifest_json().as_bytes(),
    )?;
    audio
        .artifact
        .verify(AudioArtifactBudget::DEFAULT, cx)
        .map_err(pipeline)?;
    progress("stage=audio complete");

    let movie_path = staging_directory.join("euler-disc-critique.mkv");
    let mux = if config.mux_with_ffmpeg {
        progress("stage=mux begin");
        mux_movie(config, &staging_directory, &movie_path)
    } else {
        MuxOutcome::Disabled
    };
    let completed_movie = match &mux {
        MuxOutcome::Written(path) => Some(path.clone()),
        _ => None,
    };
    let manifest_path = staging_directory.join("critique-manifest.json");
    let manifest = fixture_manifest(
        config,
        duration_s,
        run,
        &refinement,
        &trajectory_artifact,
        raw_sequence_identity,
        preview_sequence_identity,
        audio.artifact.manifest().wav().wav_identity(),
        audio.modal_parameter_set_identity,
        &audio.modal_parameter_set_disclosure,
        audio.chirp_start_hz,
        audio.chirp_end_hz,
        audio.pre_master_peak_fs,
        audio.master_gain_db,
        over_range_channels,
        gamut_mapped_pixels,
        &render_evidence,
        &denoise_evidence,
        audio.spatialization.as_ref(),
        &mux,
    );
    write_new(&manifest_path, manifest.as_bytes())?;
    if output_directory.exists() {
        return Err(CinematicFixtureError::OutputAlreadyExists(
            output_directory.to_path_buf(),
        ));
    }
    fs::rename(&staging_directory, output_directory)?;
    progress("stage=complete");
    Ok(CinematicFixtureReport {
        output_directory: output_directory.to_path_buf(),
        manifest_path: output_directory.join("critique-manifest.json"),
        wav_path: output_directory.join("sound/master.float32.wav"),
        first_preview_path: output_directory
            .join(format!("preview/frame-{:06}.png", render_frame_range.start)),
        movie_path: completed_movie.map(|_| output_directory.join("euler-disc-critique.mkv")),
    })
}

fn pipeline(error: impl fmt::Display) -> CinematicFixtureError {
    CinematicFixtureError::Pipeline(error.to_string())
}

fn create_new_file(path: &Path) -> Result<File, CinematicFixtureError> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(CinematicFixtureError::Io)
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), CinematicFixtureError> {
    let mut file = create_new_file(path)?;
    file.write_all(bytes)?;
    file.flush()?;
    Ok(())
}

fn frame_reference_times(
    frame: u32,
    frames: u32,
    trajectory_start_s: f64,
    trajectory_end_s: f64,
) -> (f64, f64, f64) {
    debug_assert!(frames > 0 && frame < frames);
    let frame_duration_s = 1.0 / f64::from(CRITIQUE_FPS);
    let frame_time_s = (f64::from(frame) + 1.0) * frame_duration_s;
    // AOV motion vectors evaluate both references, so keep them inside the
    // producer's admitted trajectory without loosening the renderer's
    // no-extrapolation contract.
    let previous_time_s = (frame_time_s - frame_duration_s).max(trajectory_start_s);
    let next_time_s = (frame_time_s + frame_duration_s).min(trajectory_end_s);
    (frame_time_s, previous_time_s, next_time_s)
}

fn composition_identity(config: &CinematicFixtureConfig, scene: ContentHash) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.composition.v2");
    hasher.update(scene.as_bytes());
    hasher.update(&config.width.to_le_bytes());
    hasher.update(&config.height.to_le_bytes());
    hasher.update(&config.frames.to_le_bytes());
    hasher.update(&CRITIQUE_FPS.to_le_bytes());
    hasher.update(&config.samples_per_pixel.to_le_bytes());
    hasher.update(&config.max_depth.to_le_bytes());
    hasher.update(&config.shutter_angle_degrees.to_le_bytes());
    hasher.update(
        &euler_scene_smoke_settings(config.width, config.height)
            .seed
            .to_le_bytes(),
    );
    hasher.update(&config.render_seed_salt.to_le_bytes());
    hasher.update(&CRITIQUE_FRAME_SEED_SCHEDULE_VERSION.to_le_bytes());
    if config.denoise_previews {
        let identity = TemporalDenoiseConfig::default()
            .identity()
            .expect("default temporal denoiser configuration is valid");
        hasher.update(identity.as_bytes());
    } else {
        hasher.update(b"temporal-denoise-disabled");
    }
    hasher.update(
        b"back-loaded-frame-boundary-stratified-shutter-v1;daily-core-aov;aces-srgb16;exposure-ev-plus-1",
    );
    hasher.finalize()
}

fn frame_render_seed(base_seed: u64, seed_salt: u64, absolute_frame: u32) -> u64 {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.frame-seed.v1");
    hasher.update(&CRITIQUE_FRAME_SEED_SCHEDULE_VERSION.to_le_bytes());
    hasher.update(&base_seed.to_le_bytes());
    hasher.update(&seed_salt.to_le_bytes());
    hasher.update(&absolute_frame.to_le_bytes());
    let digest = hasher.finalize();
    let mut seed = [0_u8; 8];
    seed.copy_from_slice(&digest.as_bytes()[..8]);
    u64::from_le_bytes(seed)
}

fn trajectory_diagnostics(artifact: &EulerRenderTrajectoryArtifact) -> String {
    let samples = artifact.trajectory().samples();
    let first = samples
        .first()
        .expect("admitted trajectory is nonempty")
        .input();
    let last = samples
        .last()
        .expect("admitted trajectory is nonempty")
        .input();
    let mut minimum = first.center_of_mass_world_m;
    let mut maximum = first.center_of_mass_world_m;
    let mut maximum_radius_m = first
        .center_of_mass_world_m
        .x
        .hypot(first.center_of_mass_world_m.y);
    let mut base_min_m = first.base_mode.map_or(0.0, |base| base.displacement_m);
    let mut base_max_m = base_min_m;
    for sample in samples {
        let input = sample.input();
        let center = input.center_of_mass_world_m;
        minimum.x = minimum.x.min(center.x);
        minimum.y = minimum.y.min(center.y);
        minimum.z = minimum.z.min(center.z);
        maximum.x = maximum.x.max(center.x);
        maximum.y = maximum.y.max(center.y);
        maximum.z = maximum.z.max(center.z);
        maximum_radius_m = maximum_radius_m.max(center.x.hypot(center.y));
        if let Some(base) = input.base_mode {
            base_min_m = base_min_m.min(base.displacement_m);
            base_max_m = base_max_m.max(base.displacement_m);
        }
    }
    format!(
        concat!(
            "stage=mechanics diagnostics ",
            "com_min_m=[{:.6},{:.6},{:.6}] com_max_m=[{:.6},{:.6},{:.6}] ",
            "max_com_radius_m={:.6} base_range_m=[{:.6},{:.6}] ",
            "first_qoi=[inclination:{:.6},precession:{:.6},spin:{:.6}] ",
            "last_qoi=[inclination:{:.6},precession:{:.6},spin:{:.6}]"
        ),
        minimum.x,
        minimum.y,
        minimum.z,
        maximum.x,
        maximum.y,
        maximum.z,
        maximum_radius_m,
        base_min_m,
        base_max_m,
        first.qois.inclination_rad,
        first.qois.precession_rad_per_s,
        first.qois.spin_rad_per_s,
        last.qois.inclination_rad,
        last.qois.precession_rad_per_s,
        last.qois.spin_rad_per_s,
    )
}

fn critique_camera(duration_s: f64) -> Result<AnimatedCamera, fs_render::camera::CameraError> {
    let eye = Point3::new(
        CRITIQUE_CAMERA_EYE_M[0],
        CRITIQUE_CAMERA_EYE_M[1],
        CRITIQUE_CAMERA_EYE_M[2],
    );
    let target = Point3::new(
        CRITIQUE_CAMERA_TARGET_M[0],
        CRITIQUE_CAMERA_TARGET_M[1],
        CRITIQUE_CAMERA_TARGET_M[2],
    );
    let physical = PhysicalCamera::try_look_at(
        eye,
        target,
        GeomVec3::new(0.0, 0.0, 1.0),
        CameraProjection::try_half_tangent(0.25)?,
        target.delta_from(eye).norm(),
        Aperture::try_circular(0.0)?,
    )?;
    AnimatedCamera::try_static(1, 0.0, duration_s, physical)
}

fn component(
    role: CinematicComponentRole,
    identity: ContentHash,
    version: u32,
) -> Result<CinematicComponentRef, CinematicFixtureError> {
    CinematicComponentRef::try_new(role, identity, version).map_err(pipeline)
}

fn listener_identity(listener: ListenerPose) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.microphone.v1");
    hasher.update(&SOUND_SYNTHESIS_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&[listener.frame as u8]);
    for value in listener
        .position_m
        .into_iter()
        .chain(listener.forward)
        .chain(listener.up)
    {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

fn dry_room_identity() -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.room.v1");
    hasher.update(&SOUND_SYNTHESIS_SCHEMA_VERSION.to_le_bytes());
    hasher.update(b"sound-room-response:dry");
    hasher.finalize()
}

fn trajectory_body_contact_chirp_bounds(
    trajectory: &EulerRenderTrajectoryArtifact,
) -> Result<(f64, f64), CinematicFixtureError> {
    let mut frequencies = trajectory.trajectory().samples().iter().map(|sample| {
        let qois = sample.input().qois;
        qois.precession_rad_per_s * det::cos(qois.inclination_rad) / core::f64::consts::TAU
    });
    let first = frequencies.next().ok_or_else(|| {
        CinematicFixtureError::Pipeline("audio trajectory has no retained samples".into())
    })?;
    if !(first.is_finite() && first > 0.0) {
        return Err(CinematicFixtureError::Pipeline(
            "trajectory-derived body-contact chirp has an invalid first frequency".into(),
        ));
    }
    let mut previous = first;
    for frequency in frequencies {
        if !(frequency.is_finite() && frequency > 0.0)
            || frequency + 32.0 * f64::EPSILON * previous.abs() < previous
        {
            return Err(CinematicFixtureError::Pipeline(
                "trajectory-derived body-contact chirp must be finite, positive, and monotone"
                    .into(),
            ));
        }
        previous = frequency;
    }
    if previous <= first {
        return Err(CinematicFixtureError::Pipeline(
            "trajectory-derived body-contact chirp must rise over the retained horizon".into(),
        ));
    }
    Ok((first, previous))
}

fn build_audio(
    trajectory: &EulerRenderTrajectoryArtifact,
    config: &CinematicFixtureConfig,
    cx: &Cx<'_>,
) -> Result<FixtureAudio, CinematicFixtureError> {
    let audio_frame_count = u64::from(config.frames) * 2_000;
    let controls = EulerControlStream::try_derive(trajectory.trajectory(), cx).map_err(pipeline)?;
    let (chirp_start_hz, chirp_end_hz) = trajectory_body_contact_chirp_bounds(trajectory)?;
    if chirp_end_hz >= CRITIQUE_DECLARED_AUDIO_SOURCE_BANDWIDTH_HZ {
        return Err(CinematicFixtureError::Pipeline(format!(
            "trajectory-derived chirp reaches {chirp_end_hz:.6} Hz, outside the declared reconstruction bandwidth"
        )));
    }
    let preset = representative_modal_preset(RepresentativeDiscMaterial::StainlessSteel);
    let spatial_rules = preset
        .modes()
        .iter()
        .map(|mode| ModeContactParticipationRule {
            mode_id: mode.mode_id,
            shape: if mode.component == SoundModalComponent::Disc {
                ContactModeShape::AzimuthalCosine {
                    harmonic: 1,
                    phase_rad: 0.0,
                }
            } else {
                ContactModeShape::Uniform
            },
        })
        .collect();
    let acoustic_rig_identity = hash_domain(
        "org.frankensim.euler-critique.representative-acoustic-rig.v1",
        trajectory
            .trajectory()
            .metadata()
            .base_model_identity
            .as_bytes(),
    );
    let modal_parameters = EulerModalParameterSet::try_admit(
        EulerModalParameterSetInput {
            authority: ModalPresetAuthority::RepresentativeUncalibrated,
            specimen_identity: trajectory
                .trajectory()
                .metadata()
                .specimen_profile_identity,
            rig_identity: acoustic_rig_identity,
            disclosure: "Representative stainless-steel disc, thick-glass, and base modes; not measured for the Thorne specimen or support"
                .to_owned(),
            calibration: None,
            model: ModalSynthesisModelInput {
            sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
            modes: preset.modes().to_vec(),
            budget: ModalSynthesisBudget::reference_film(audio_frame_count),
        },
        },
        cx,
    )
    .map_err(pipeline)?;
    let modal_parameter_set_identity = modal_parameters.identity();
    let modal_parameter_set_disclosure = modal_parameters.disclosure().to_owned();
    let modal = modal_parameters.into_model();
    let mappings = vec![SoundExcitationControl {
        channel: SoundExcitationChannel::RollingSignedWorkRate,
        target_component: SoundModalComponent::Disc,
        source_scale: 2.0,
    }];
    let interval_count = controls.audio().len();
    let mapper = AudioExcitationMapper::try_new(
        trajectory,
        &controls,
        &modal,
        AudioExcitationModelInput {
            mappings: mappings.clone(),
            reduction: AudioExcitationReduction::RawIntervals,
            spatial_policy: ContactParticipationPolicy::ContactCoordinates {
                rules: spatial_rules,
            },
            artistic_texture: None,
            budget: AudioExcitationBudget::reference_film(interval_count),
        },
        cx,
    )
    .map_err(pipeline)?;
    let selected_count = mapper.grid().interval_count;
    let mut mapped_intervals = Vec::new();
    mapped_intervals
        .try_reserve_exact(selected_count)
        .map_err(|_| CinematicFixtureError::Pipeline("audio interval allocation refused".into()))?;
    let mut map_checkpoint = mapper.initial_checkpoint(cx).map_err(pipeline)?;
    while map_checkpoint.next_interval_index() < selected_count {
        let remaining = selected_count - map_checkpoint.next_interval_index();
        let chunk_size = remaining.min(65_536);
        let chunk = mapper
            .map_next_chunk(
                &map_checkpoint,
                NonZeroUsize::new(chunk_size).expect("positive remaining interval count"),
                cx,
            )
            .map_err(pipeline)?;
        mapped_intervals.extend(chunk.intervals);
        map_checkpoint = chunk.successor;
    }
    if mapped_intervals.len() != selected_count {
        return Err(CinematicFixtureError::Pipeline(format!(
            "audio mapping retained {} of {selected_count} source intervals",
            mapped_intervals.len()
        )));
    }

    let video_clock = CinematicClock::try_new(
        CinematicClockDomain::Video,
        CRITIQUE_FPS,
        1,
        0,
        i64::from(config.frames),
    )
    .map_err(pipeline)?;
    let audio_clock = CinematicClock::try_new(
        CinematicClockDomain::Audio,
        SOUND_MASTER_SAMPLE_RATE_HZ,
        1,
        0,
        i64::try_from(audio_frame_count)
            .map_err(|_| CinematicFixtureError::Pipeline("audio clock overflow".into()))?,
    )
    .map_err(pipeline)?;
    let resampler = AudioResampler::try_new(
        &mapper,
        &modal,
        mapped_intervals,
        AudioResamplingModelInput {
            video_clock,
            audio_clock,
            declared_source_bandwidth_hz: CRITIQUE_DECLARED_AUDIO_SOURCE_BANDWIDTH_HZ,
            filter: AudioReconstructionFilterSpec {
                passband_edge_hz: 2_000.0,
                stopband_edge_hz: 4_800.0,
                half_length: 128,
                maximum_passband_ripple_db: 0.1,
                minimum_stopband_attenuation_db: 80.0,
                response_grid_intervals: 8_192,
            },
            boundary_policy: AudioResamplingBoundaryPolicy::HalfSampleEvenReflectionV1,
            event_fractional_delay: AudioEventFractionalDelay::LinearTwoBoundaryV1,
            budget: AudioResamplingBudget::reference_film(),
        },
        cx,
    )
    .map_err(pipeline)?;
    let timeline_identity = {
        let mut hasher = DomainHasher::new("org.frankensim.euler-critique.master-clocks.v1");
        hasher.update(&CRITIQUE_FPS.to_le_bytes());
        hasher.update(&SOUND_MASTER_SAMPLE_RATE_HZ.to_le_bytes());
        hasher.update(&0_i64.to_le_bytes());
        hasher.update(&config.frames.to_le_bytes());
        hasher.update(&0_i64.to_le_bytes());
        hasher.update(&audio_frame_count.to_le_bytes());
        hasher.finalize()
    };
    let listener = ListenerPose {
        frame: ListenerFrame::AnimatedCamera,
        position_m: [0.0, 0.0, 0.0],
        forward: [0.0, 0.0, -1.0],
        up: [0.0, 1.0, 0.0],
    };
    let sound = SoundSynthesisConfig::try_admit(SoundSynthesisInput {
        schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
        authority: SoundAuthority::PhysicallyInformed,
        trajectory: component(
            CinematicComponentRole::Trajectory,
            trajectory.receipt().artifact_identity(),
            u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION),
        )?,
        excitation: component(
            CinematicComponentRole::AudioExcitation,
            mapper.identity(),
            AUDIO_EXCITATION_ALGORITHM_VERSION,
        )?,
        sound_model: component(
            CinematicComponentRole::SoundModel,
            modal.identity(),
            crate::MODAL_SYNTHESIS_ALGORITHM_VERSION,
        )?,
        microphone: component(
            CinematicComponentRole::Microphone,
            listener_identity(listener),
            1,
        )?,
        room: component(CinematicComponentRole::Room, dry_room_identity(), 1)?,
        timeline: component(CinematicComponentRole::Timeline, timeline_identity, 1)?,
        video_clock,
        audio_clock,
        channel_layout: SoundChannelLayout::Stereo,
        listener,
        excitation_controls: mappings,
        modes: modal.modes().to_vec(),
        room_response: SoundRoomResponse::Dry,
        amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db: 6.0 },
        trajectory_disposition: SoundTrajectoryDisposition::HorizonCensored,
        terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
            fade_sample_frames: TERMINAL_FADE_SAMPLE_FRAMES,
        },
        resampler_identity: resampler.identity(),
        resampler_version: AUDIO_RESAMPLING_ALGORITHM_VERSION,
        filter_identity: resampler.filter_identity(),
        filter_version: AUDIO_RECONSTRUCTION_FILTER_VERSION,
        assumptions: vec![
            SoundModelAssumption::LinearModalSuperposition,
            SoundModelAssumption::TimeInvariantDamping,
            SoundModelAssumption::DeclaredExcitationCompleteness,
            SoundModelAssumption::DeclaredRoomResponse,
        ],
        calibration: None,
    })
    .map_err(pipeline)?;
    mapper
        .validate_sound_configuration(&sound)
        .map_err(pipeline)?;
    modal
        .validate_sound_configuration(&sound)
        .map_err(pipeline)?;
    resampler
        .validate_sound_configuration(&sound)
        .map_err(pipeline)?;

    let mut stems = Vec::new();
    stems
        .try_reserve_exact(audio_frame_count as usize)
        .map_err(|_| CinematicFixtureError::Pipeline("audio stem allocation refused".into()))?;
    let mut resampling_checkpoint = resampler.initial_checkpoint(cx).map_err(pipeline)?;
    let mut modal_checkpoint = modal.initial_checkpoint(cx).map_err(pipeline)?;
    while resampling_checkpoint.next_audio_frame_offset() < audio_frame_count {
        let remaining = audio_frame_count - resampling_checkpoint.next_audio_frame_offset();
        let chunk_frames = usize::try_from(remaining.min(65_536))
            .map_err(|_| CinematicFixtureError::Pipeline("audio chunk length overflow".into()))?;
        let resampled = resampler
            .resample_next_chunk(
                &sound,
                &resampling_checkpoint,
                NonZeroUsize::new(chunk_frames).expect("positive remaining audio frames"),
                cx,
            )
            .map_err(pipeline)?;
        let synthesized = resampled
            .synthesize_modal(&modal, &modal_checkpoint, cx)
            .map_err(pipeline)?;
        stems.extend(synthesized.stem_frames);
        resampling_checkpoint = resampled.successor;
        modal_checkpoint = synthesized.successor;
    }
    if stems.len() as u64 != audio_frame_count {
        return Err(CinematicFixtureError::Pipeline(format!(
            "modal synthesis emitted {} of {audio_frame_count} frames",
            stems.len()
        )));
    }
    apply_initial_fade(&mut stems, INITIAL_FADE_SAMPLE_FRAMES)?;
    if !config.spatialize_audio {
        apply_terminal_fade(&mut stems, TERMINAL_FADE_SAMPLE_FRAMES)?;
    }
    let mut mix = AudioDryMixSpec {
        disc: StemGainPan {
            gain_db: -6.0,
            pan: -0.08,
        },
        glass_plate: StemGainPan {
            gain_db: -4.0,
            pan: 0.08,
        },
        base_assembly: StemGainPan {
            gain_db: -8.0,
            pan: 0.0,
        },
        master_gain_db: 0.0,
    };
    let spatialized = config
        .spatialize_audio
        .then(|| {
            spatialize_fixture_audio(
                trajectory,
                &stems,
                sound.receipt().configuration_identity,
                cx,
            )
        })
        .transpose()?;
    let spatial_premaster = spatialized
        .as_ref()
        .map(|output| {
            apply_spatial_terminal_fade(output.samples(), TERMINAL_FADE_SAMPLE_FRAMES, cx)
        })
        .transpose()?;
    let provisional_meters = if let Some(samples) = &spatial_premaster {
        measure_audio(samples, AudioArtifactBudget::DEFAULT, cx).map_err(pipeline)?
    } else {
        let provisional =
            mix_dry_modal_stems(&stems, mix, AudioArtifactBudget::DEFAULT, cx).map_err(pipeline)?;
        measure_audio(&provisional, AudioArtifactBudget::DEFAULT, cx).map_err(pipeline)?
    };
    let provisional_peak = provisional_meters
        .sample_peak_fs
        .max(provisional_meters.true_peak_estimate_fs);
    const TARGET_PEAK_FS: f64 = 0.45;
    if provisional_peak <= f64::MIN_POSITIVE {
        return Err(CinematicFixtureError::Pipeline(
            "mechanics-derived sound is silent and cannot be mastered".into(),
        ));
    }
    let master_gain_db =
        20.0 * det::ln(TARGET_PEAK_FS / provisional_peak) / core::f64::consts::LN_10;
    if !(-120.0..=120.0).contains(&master_gain_db) {
        return Err(CinematicFixtureError::Pipeline(format!(
            "mechanics-derived sound needs {:.3} dB mastering gain, beyond the admitted range",
            master_gain_db
        )));
    }
    mix.master_gain_db = master_gain_db;
    let (artifact, spatialization) = if let (Some(output), Some(premaster)) =
        (spatialized, spatial_premaster)
    {
        let mastered = apply_stereo_master_gain(&premaster, master_gain_db, cx)?;
        let mastered_output_identity = mastered_spatial_output_identity(&output, master_gain_db);
        let artifact = SoundWavArtifact::try_build(
            &sound,
            AudioMasterSource::SpatializedStereo {
                frames: &mastered,
                spatialization_identity: mastered_output_identity,
                source_synthesis: sound.receipt(),
            },
            WavSampleEncoding::Float32,
            WavMetadata::try_new(Some(
                "FrankenSim Euler-disc mechanics-driven spatial critique preview; uncalibrated"
                    .to_owned(),
            ))
            .map_err(pipeline)?,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .map_err(pipeline)?;
        let diagnostics = output.diagnostics();
        let evidence = FixtureSpatialAudioEvidence {
            config_identity: output.config_identity(),
            input_identity: output.input_identity(),
            raw_output_identity: output.output_identity(),
            mastered_output_identity,
            authority: output.authority(),
            maximum_distance_m: diagnostics.maximum_distance_m,
            maximum_delay_frames: diagnostics.maximum_delay_frames,
            discarded_tail_frames: diagnostics.discarded_tail_frames,
        };
        (artifact, Some(evidence))
    } else {
        let artifact = SoundWavArtifact::try_build(
            &sound,
            AudioMasterSource::DryModalStems {
                frames: &stems,
                mix,
                source_synthesis: sound.receipt(),
            },
            WavSampleEncoding::Float32,
            WavMetadata::try_new(Some(
                "FrankenSim Euler-disc mechanics-driven dry critique preview; uncalibrated"
                    .to_owned(),
            ))
            .map_err(pipeline)?,
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .map_err(pipeline)?;
        (artifact, None)
    };
    Ok(FixtureAudio {
        artifact,
        modal_parameter_set_identity,
        modal_parameter_set_disclosure,
        chirp_start_hz,
        chirp_end_hz,
        pre_master_peak_fs: provisional_peak,
        master_gain_db,
        spatialization,
    })
}

fn spatialize_fixture_audio(
    trajectory: &EulerRenderTrajectoryArtifact,
    stems: &[crate::ModalStemFrame],
    sound_configuration_identity: ContentHash,
    cx: &Cx<'_>,
) -> Result<SpatialAudioOutput, CinematicFixtureError> {
    let (disc_positions, glass_positions) = fixture_audio_positions(trajectory, stems.len(), cx)?;
    let listener = spatial_listener_pose();
    let spatializer = OfflineSpatializer::try_new(
        SpatialAudioConfig {
            sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
            speed_of_sound_m_per_s: 343.0,
            minimum_distance_m: CRITIQUE_SPATIAL_REFERENCE_DISTANCE_M,
            delay_policy: SpatialDelayPolicy::LinearFloorCeil,
            output_horizon: SpatialOutputHorizon::ClampToInputFrames,
            microphone_directivity: MicrophoneDirectivity::Cardioid {
                rear_floor_gain: 0.15,
            },
            authority: SpatialAudioAuthority::Artistic,
            budget: SpatialAudioBudget::PREVIEW,
        },
        cx,
    )
    .map_err(pipeline)?;
    let sources = [
        SpatialAudioSource {
            source_identity: spatial_stem_identity(sound_configuration_identity, b"disc"),
            signal: SpatialMonoSignal::ModalStemFrames {
                frames: stems,
                component: SpatialStemComponent::Disc,
            },
            positions: SourcePositionTrack::PerFrame(&disc_positions),
            // Exact linear equivalents of the established -6/-4/-8 dB mix.
            gain_linear: 0.501_187_233_627_272_2,
            authority: SpatialAudioAuthority::Artistic,
        },
        SpatialAudioSource {
            source_identity: spatial_stem_identity(sound_configuration_identity, b"glass-plate"),
            signal: SpatialMonoSignal::ModalStemFrames {
                frames: stems,
                component: SpatialStemComponent::GlassPlate,
            },
            positions: SourcePositionTrack::PerFrame(&glass_positions),
            gain_linear: 0.630_957_344_480_193_2,
            authority: SpatialAudioAuthority::Artistic,
        },
        SpatialAudioSource {
            source_identity: spatial_stem_identity(sound_configuration_identity, b"base-assembly"),
            signal: SpatialMonoSignal::ModalStemFrames {
                frames: stems,
                component: SpatialStemComponent::BaseAssembly,
            },
            positions: SourcePositionTrack::Static([0.0, 0.0, -0.02]),
            gain_linear: 0.398_107_170_553_497_2,
            authority: SpatialAudioAuthority::Artistic,
        },
    ];
    spatializer
        .spatialize(
            SpatialAudioRenderInput {
                sources: &sources,
                listener: ListenerPoseTrack::Static(listener),
                room_ir: None,
            },
            cx,
        )
        .map_err(pipeline)
}

fn spatial_stem_identity(
    sound_configuration_identity: ContentHash,
    component: &[u8],
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.spatial-stem-source.v1");
    hasher.update(sound_configuration_identity.as_bytes());
    hasher.update(component);
    hasher.finalize()
}

fn spatial_listener_pose() -> SpatialListenerPose {
    let eye = CRITIQUE_CAMERA_EYE_M;
    let target = CRITIQUE_CAMERA_TARGET_M;
    let toward: [f64; 3] = [target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]];
    let forward_norm =
        det::sqrt(toward[0] * toward[0] + toward[1] * toward[1] + toward[2] * toward[2]);
    let forward_unit = toward.map(|value| value / forward_norm);
    let unnormalized_right = [forward_unit[1], -forward_unit[0], 0.0];
    let right_norm = det::sqrt(
        unnormalized_right[0] * unnormalized_right[0]
            + unnormalized_right[1] * unnormalized_right[1],
    );
    SpatialListenerPose {
        position_m: eye,
        forward_unit,
        right_unit: unnormalized_right.map(|value| value / right_norm),
    }
}

fn fixture_audio_positions(
    trajectory: &EulerRenderTrajectoryArtifact,
    audio_frames: usize,
    cx: &Cx<'_>,
) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), CinematicFixtureError> {
    let samples = trajectory.trajectory().samples();
    let metadata = trajectory.trajectory().metadata();
    let first = samples.first().ok_or_else(|| {
        CinematicFixtureError::Pipeline("spatial audio received an empty trajectory".into())
    })?;
    let last = samples
        .last()
        .expect("nonempty trajectory has a last sample");
    let mut disc_positions = Vec::new();
    disc_positions
        .try_reserve_exact(audio_frames)
        .map_err(|_| {
            CinematicFixtureError::Pipeline("disc spatial position allocation refused".into())
        })?;
    let mut glass_positions = Vec::new();
    glass_positions
        .try_reserve_exact(audio_frames)
        .map_err(|_| {
            CinematicFixtureError::Pipeline("glass spatial position allocation refused".into())
        })?;
    let mut upper = 1_usize;
    for frame in 0..audio_frames {
        if frame % 4_096 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicFixtureError::Cancelled)?;
        }
        let time_s = frame as f64 / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
        while upper < samples.len() && samples[upper].input().time_s < time_s {
            upper += 1;
        }
        let (disc_position, glass_position) = if time_s <= first.input().time_s {
            fixture_audio_position_pair(first.input(), metadata)?
        } else if upper >= samples.len() {
            fixture_audio_position_pair(last.input(), metadata)?
        } else {
            let left = samples[upper - 1].input();
            let right = samples[upper].input();
            let interval_s = right.time_s - left.time_s;
            if !interval_s.is_finite() || interval_s <= 0.0 {
                return Err(CinematicFixtureError::Pipeline(
                    "spatial trajectory time interval is not positive".into(),
                ));
            }
            let alpha = ((time_s - left.time_s) / interval_s).clamp(0.0, 1.0);
            let (left_disc, left_glass) = fixture_audio_position_pair(left, metadata)?;
            let (right_disc, right_glass) = fixture_audio_position_pair(right, metadata)?;
            (
                left_disc.scale(1.0 - alpha).add(right_disc.scale(alpha)),
                left_glass.scale(1.0 - alpha).add(right_glass.scale(alpha)),
            )
        };
        disc_positions.push([disc_position.x, disc_position.y, disc_position.z]);
        glass_positions.push([glass_position.x, glass_position.y, glass_position.z]);
    }
    cx.checkpoint()
        .map_err(|_| CinematicFixtureError::Cancelled)?;
    Ok((disc_positions, glass_positions))
}

fn fixture_audio_position_pair(
    input: &crate::RenderTrajectorySampleInput,
    metadata: &crate::RenderTrajectoryMetadata,
) -> Result<(Vec3, Vec3), CinematicFixtureError> {
    let base = input.base_mode.ok_or_else(|| {
        CinematicFixtureError::Pipeline("spatial audio trajectory omitted base state".into())
    })?;
    let base_offset_world = metadata
        .base_frame
        .orientation_base_to_world
        .rotate_body_to_world(Vec3::new(0.0, 0.0, base.displacement_m));
    Ok((
        input.center_of_mass_world_m,
        metadata.base_frame.origin_world_m.add(base_offset_world),
    ))
}

fn apply_spatial_terminal_fade(
    samples: &[StereoSample],
    fade_sample_frames: u32,
    cx: &Cx<'_>,
) -> Result<Vec<StereoSample>, CinematicFixtureError> {
    let fade_frames = usize::try_from(fade_sample_frames)
        .map_err(|_| CinematicFixtureError::Pipeline("spatial fade length overflow".into()))?;
    if fade_frames < 2 || fade_frames > samples.len() {
        return Err(CinematicFixtureError::Pipeline(format!(
            "spatial fade length {fade_frames} is incompatible with {} frames",
            samples.len()
        )));
    }
    let fade_start = samples.len() - fade_frames;
    let denominator = (fade_frames - 1) as f64;
    let mut faded = Vec::new();
    faded
        .try_reserve_exact(samples.len())
        .map_err(|_| CinematicFixtureError::Pipeline("spatial fade allocation refused".into()))?;
    for (index, sample) in samples.iter().copied().enumerate() {
        if index % 4_096 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicFixtureError::Cancelled)?;
        }
        let gain = if index < fade_start {
            1.0
        } else {
            (samples.len() - 1 - index) as f64 / denominator
        };
        faded.push(StereoSample {
            left_fs: sample.left_fs * gain,
            right_fs: sample.right_fs * gain,
        });
    }
    cx.checkpoint()
        .map_err(|_| CinematicFixtureError::Cancelled)?;
    Ok(faded)
}

fn apply_stereo_master_gain(
    samples: &[StereoSample],
    gain_db: f64,
    cx: &Cx<'_>,
) -> Result<Vec<StereoSample>, CinematicFixtureError> {
    let gain = det::exp(gain_db * core::f64::consts::LN_10 / 20.0);
    if !gain.is_finite() || gain <= 0.0 {
        return Err(CinematicFixtureError::Pipeline(
            "spatial master gain is not finite and positive".into(),
        ));
    }
    let mut mastered = Vec::new();
    mastered
        .try_reserve_exact(samples.len())
        .map_err(|_| CinematicFixtureError::Pipeline("spatial master allocation refused".into()))?;
    for (index, sample) in samples.iter().copied().enumerate() {
        if index % 4_096 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicFixtureError::Cancelled)?;
        }
        mastered.push(StereoSample {
            left_fs: sample.left_fs * gain,
            right_fs: sample.right_fs * gain,
        });
    }
    cx.checkpoint()
        .map_err(|_| CinematicFixtureError::Cancelled)?;
    Ok(mastered)
}

fn mastered_spatial_output_identity(
    output: &SpatialAudioOutput,
    master_gain_db: f64,
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.spatial-master.v1");
    hasher.update(output.output_identity().as_bytes());
    hasher.update(&master_gain_db.to_bits().to_le_bytes());
    hasher.update(&TERMINAL_FADE_SAMPLE_FRAMES.to_le_bytes());
    hasher.update(b"post-spatial-linear-terminal-fade;fs-math-exp-linear-gain;no-limiter");
    hasher.finalize()
}

fn apply_terminal_fade(
    stems: &mut [crate::ModalStemFrame],
    fade_sample_frames: u32,
) -> Result<(), CinematicFixtureError> {
    let fade_frames = usize::try_from(fade_sample_frames)
        .map_err(|_| CinematicFixtureError::Pipeline("terminal fade length overflow".into()))?;
    if fade_frames < 2 || fade_frames > stems.len() {
        return Err(CinematicFixtureError::Pipeline(format!(
            "terminal fade length {fade_frames} is incompatible with {} stem frames",
            stems.len()
        )));
    }
    let fade_start = stems.len() - fade_frames;
    let denominator = (fade_frames - 1) as f64;
    for (index, frame) in stems[fade_start..].iter_mut().enumerate() {
        let gain = (fade_frames - 1 - index) as f64 / denominator;
        frame.disc_fs *= gain;
        frame.glass_plate_fs *= gain;
        frame.base_assembly_fs *= gain;
    }
    Ok(())
}

fn apply_initial_fade(
    stems: &mut [crate::ModalStemFrame],
    fade_sample_frames: u32,
) -> Result<(), CinematicFixtureError> {
    let fade_frames = usize::try_from(fade_sample_frames)
        .map_err(|_| CinematicFixtureError::Pipeline("initial fade length overflow".into()))?;
    if fade_frames < 2 || fade_frames > stems.len() {
        return Err(CinematicFixtureError::Pipeline(format!(
            "initial fade length {fade_frames} is incompatible with {} stem frames",
            stems.len()
        )));
    }
    let denominator = (fade_frames - 1) as f64;
    for (index, frame) in stems[..fade_frames].iter_mut().enumerate() {
        let gain = index as f64 / denominator;
        frame.disc_fs *= gain;
        frame.glass_plate_fs *= gain;
        frame.base_assembly_fs *= gain;
    }
    Ok(())
}

fn mux_movie(
    config: &CinematicFixtureConfig,
    output_directory: &Path,
    movie_path: &Path,
) -> MuxOutcome {
    let Some(movie_name) = movie_path.file_name() else {
        return MuxOutcome::Unavailable("movie path has no file name".into());
    };
    let status = Command::new(&config.ffmpeg_executable)
        .current_dir(output_directory)
        .args([
            "-hide_banner",
            "-loglevel",
            "warning",
            "-n",
            "-framerate",
            "24",
            "-i",
            "preview/frame-%06d.png",
            "-i",
            "sound/master.float32.wav",
            "-c:v",
            "libsvtav1",
            "-crf",
            "24",
            "-preset",
            "6",
            "-pix_fmt",
            "yuv420p10le",
            "-c:a",
            "libopus",
            "-b:a",
            "192k",
            "-shortest",
        ])
        .arg(movie_name)
        .status();
    match status {
        Ok(status) if status.success() && movie_path.is_file() => {
            MuxOutcome::Written(PathBuf::from("euler-disc-critique.mkv"))
        }
        Ok(status) => MuxOutcome::Failed(status.code().unwrap_or(-1)),
        Err(error) => MuxOutcome::Unavailable(error.to_string()),
    }
}

#[allow(clippy::too_many_arguments)]
fn fixture_manifest(
    config: &CinematicFixtureConfig,
    duration_s: f64,
    run: &ReducedDecayRun,
    refinement: &RefinementEvidence,
    trajectory: &EulerRenderTrajectoryArtifact,
    raw_sequence_identity: ContentHash,
    preview_sequence_identity: ContentHash,
    wav_identity: ContentHash,
    modal_parameter_set_identity: ContentHash,
    modal_parameter_set_disclosure: &str,
    chirp_start_hz: f64,
    chirp_end_hz: f64,
    audio_pre_master_peak_fs: f64,
    audio_master_gain_db: f64,
    over_range_channels: u64,
    gamut_mapped_pixels: u64,
    render: &FixtureRenderEvidence,
    denoise: &FixtureDenoiseEvidence,
    spatialization: Option<&FixtureSpatialAudioEvidence>,
    mux: &MuxOutcome,
) -> String {
    let rendered_frames = config
        .render_frame_range()
        .expect("validated fixture configuration retains a valid frame window");
    let rendered_frame_count = rendered_frames.len();
    let complete_sequence =
        rendered_frames.start == 0 && rendered_frame_count == config.frames as usize;
    let first_source_sample = run
        .samples
        .first()
        .expect("completed reduced-decay run retains a first sample");
    let last_source_sample = run
        .samples
        .last()
        .expect("completed reduced-decay run retains a last sample");
    let trajectory_samples = trajectory.trajectory().samples();
    let first_visual_sample = trajectory_samples
        .first()
        .expect("admitted trajectory retains a first sample")
        .input();
    let last_visual_sample = trajectory_samples
        .last()
        .expect("admitted trajectory retains a last sample")
        .input();
    let relative_energy_defect = run.energy_closure_residual_j.abs()
        / first_source_sample.energy_j.abs().max(f64::MIN_POSITIVE);
    let trajectory_identity = trajectory.receipt().artifact_identity();
    let mux_json = match mux {
        MuxOutcome::Disabled => "{\"status\":\"disabled\"}".to_owned(),
        MuxOutcome::Unavailable(message) => format!(
            "{{\"status\":\"unavailable\",\"detail\":\"{}\"}}",
            json_escape(message)
        ),
        MuxOutcome::Failed(code) => {
            format!("{{\"status\":\"failed\",\"exit_code\":{code}}}")
        }
        MuxOutcome::Written(path) => format!(
            "{{\"status\":\"written\",\"path\":\"{}\"}}",
            json_escape(&path.display().to_string())
        ),
    };
    let spatialization_json = spatialization.map_or_else(
        || "{\"enabled\":false,\"path\":\"canonical-dry-stereo\"}".to_owned(),
        |spatial| {
            format!(
                concat!(
                    "{{\"enabled\":true,\"path\":\"direct-point-source-stereo\",",
                    "\"config_identity\":\"{}\",\"input_identity\":\"{}\",",
                    "\"raw_output_identity\":\"{}\",\"mastered_output_identity\":\"{}\",",
                    "\"authority\":\"{}\",\"room_response\":\"dry\",",
                    "\"output_horizon\":\"clamp-to-input-frames\",",
                    "\"post_spatial_terminal_fade_sample_frames\":{},",
                    "\"maximum_distance_m\":{:.17e},\"maximum_delay_frames\":{:.17e},",
                    "\"discarded_tail_frames\":{}}}"
                ),
                spatial.config_identity.to_hex(),
                spatial.input_identity.to_hex(),
                spatial.raw_output_identity.to_hex(),
                spatial.mastered_output_identity.to_hex(),
                spatial.authority.code(),
                TERMINAL_FADE_SAMPLE_FRAMES,
                spatial.maximum_distance_m,
                spatial.maximum_delay_frames,
                spatial.discarded_tail_frames,
            )
        },
    );
    let specimen = run
        .provenance
        .literature_specimen
        .as_ref()
        .expect("source-bound fixture run retains its specimen declaration");
    let raw_profile = if config.retain_full_aov_exr {
        "daily-core-aov-float"
    } else {
        "linear-srgb-beauty-float"
    };
    let denoise_pipeline = if config.denoise_previews {
        TEMPORAL_DENOISE_PIPELINE_VERSION
    } else {
        "disabled"
    };
    let modal_disclosure = json_escape(modal_parameter_set_disclosure);
    format!(
        concat!(
            "{{\n",
            "  \"schema\": \"frankensim-euler-cinematic-critique-v2\",\n",
            "  \"authority\": \"source-bound analytical simulation visualization; physically informed but uncalibrated synthesis; artistic spatial presentation\",\n",
            "  \"video\": {{\"width\": {width}, \"height\": {height}, \"sequence_frames\": {sequence_frames}, \"rendered_first_frame\": {rendered_first}, \"rendered_frame_count\": {rendered_count}, \"complete_sequence\": {complete}, \"fps\": {fps}, \"duration_s\": {duration:.9}, \"spp\": {spp}, \"render_seed_salt\": {seed_salt}, \"max_depth\": {max_depth}, \"shutter_angle_degrees\": {shutter_angle}, \"shutter_duration_s\": {shutter_duration:.17e}, \"shutter_convention\": \"back-loaded-frame-boundary\", \"terminal_cutoff_in_final_shutter\": true, \"shutter_distribution\": \"stratified-counter-v1\", \"frame_seed_schedule_version\": {seed_version}, \"denoise_requested\": {denoise_requested}, \"raw_exr_profile\": \"{raw_profile}\", \"exposure_ev\": 1, \"raw_sequence_identity\": \"{raw_sequence}\", \"preview_sequence_identity\": \"{preview_sequence}\", \"over_range_linear_channels\": {over_range}, \"gamut_mapped_pixels\": {gamut_mapped}}},\n",
            "  \"render_execution\": {{\"policy\": \"deterministic-parked-crew-tile-v1\", \"requested_workers\": {requested_workers}, \"maximum_effective_workers\": {effective_workers}, \"tile_width\": {tile_width}, \"tile_height\": {tile_height}, \"memory_limit_bytes\": {memory_limit}, \"maximum_peak_memory_bytes\": {peak_memory}, \"measured_frames\": {measured_frames}, \"timing_ns\": {{\"setup\": {setup_ns}, \"traversal\": {traversal_ns}, \"tile_compute_sum\": {compute_ns}, \"tile_merge_sum\": {merge_ns}, \"publication\": {publication_ns}, \"idle_worker_capacity\": {idle_ns}}}}},\n",
            "  \"denoise\": {{\"requested\": {denoise_requested}, \"applied_frames\": {denoised_frames}, \"pipeline\": \"{denoise_pipeline}\", \"authority\": \"biased-display-derivative\", \"maximum_retained_bytes\": {denoise_bytes}, \"maximum_history_frames\": {history_frames}}},\n",
            "  \"mechanics\": {{\"model\": \"Thorne-2026-small-angle-rolling-plus-Bildsten-boundary-layer\", \"source_id\": \"{source_id}\", \"model_authority\": \"{model_authority}\", \"physical_validation\": \"{physical_validation}\", \"specimen\": {{\"diameter_m\": {diameter:.17e}, \"thickness_m\": {thickness:.17e}, \"mass_kg\": {mass:.17e}, \"outer_fillet_radius_m\": {fillet:.17e}}}, \"integration\": {{\"coarse_timestep_s\": {coarse_dt:.17e}, \"fine_timestep_s\": {fine_dt:.17e}, \"source_sample_count\": {source_samples}, \"source_duration_s\": {source_duration:.17e}, \"retained_tail_sample_count\": {tail_samples}, \"retained_tail_duration_s\": {tail_duration:.17e}, \"terminal\": \"{terminal:?}\", \"positive_validity_cutoff_rad\": {cutoff:.17e}}}, \"refinement\": {{\"terminal_time_difference_s\": {refine_time:.17e}, \"total_work_difference_j\": {refine_work:.17e}, \"claim\": \"single dt/dt2 consistency pair for the encoded analytical model; not experimental validation or an asymptotic-order certificate\"}}, \"channels\": {{\"rolling_coefficient_mu\": {rolling_mu:.17e}, \"rolling_work_j\": {rolling_work:.17e}, \"boundary_layer_work_j\": {gas_work:.17e}}}, \"first_retained_qoi\": {{\"inclination_rad\": {first_theta:.17e}, \"precession_rad_per_s\": {first_precession:.17e}, \"spin_rad_per_s\": {first_spin:.17e}}}, \"last_retained_qoi\": {{\"inclination_rad\": {last_theta:.17e}, \"precession_rad_per_s\": {last_precession:.17e}, \"spin_rad_per_s\": {last_spin:.17e}}}, \"energy\": {{\"initial_j\": {initial_energy:.17e}, \"final_j\": {final_energy:.17e}, \"closure_residual_j\": {energy_residual:.17e}, \"relative_abs_residual\": {relative_residual:.17e}}}, \"trajectory_identity\": \"{trajectory_identity}\"}},\n",
            "  \"audio\": {{\"sample_rate_hz\": {audio_rate}, \"wav_identity\": \"{wav_identity}\", \"authority\": \"physically-informed-uncalibrated\", \"calibrated\": false, \"procedural_texture\": false, \"excitation\": \"published rolling work rate times uncalibrated 2 N/W transfer\", \"contact_phase\": \"body-contact azimuth; harmonic one; instantaneous rate Omega*cos(theta)\", \"chirp_start_hz\": {chirp_start:.17e}, \"chirp_end_hz\": {chirp_end:.17e}, \"declared_source_bandwidth_hz\": {audio_bandwidth:.17e}, \"modal_parameter_set_identity\": \"{modal_identity}\", \"modal_parameter_set_disclosure\": \"{modal_disclosure}\", \"pre_master_peak_fs\": {pre_master_peak:.17e}, \"master_gain_db\": {master_gain:.9}, \"initial_fade_sample_frames\": {initial_fade}, \"terminal_fade_sample_frames\": {terminal_fade}, \"terminal_fade_application\": \"exactly once: dry stems for dry output or post-propagation stereo for spatial output\", \"mix_policy\": \"one content-derived digital mastering gain to 0.45 FS; no limiter\", \"spatialization\": {spatialization}}},\n",
            "  \"mux\": {mux},\n",
            "  \"no_claims\": [\"the analytical model reproduces published equations and fitted rolling coefficient but is not a full fluid-structure-contact solve or a raw measured trajectory\", \"the positive inclination cutoff is horizon censoring, not theta zero, loss of contact, or a resolved terminal impact\", \"the harmonic-one contact shape, modal frequencies, damping, masses, radiation gains, and rolling-power-to-force transfer are representative rather than measured for this specimen and rig\", \"declared excitation completeness means complete only for this authored reduced-channel sonification, not complete physical acoustic forcing\", \"the waveform, loudness, spectral envelope, terminal chatter, microphone, room, HRTF, and sound-pressure level are not experimentally validated\", \"spatial output clamps propagation tails, so listener audio does not claim to contain the exact source cutoff sample\", \"digital mastering is presentation normalization, not a pascal or SPL prediction\", \"the radial spin fiducial is visualization-only and excluded from specimen geometry, contact, and mass\", \"image quality is final only after native-4K sample-rung review and complete-sequence verification\"]\n",
            "}}\n"
        ),
        width = config.width,
        height = config.height,
        sequence_frames = config.frames,
        rendered_first = rendered_frames.start,
        rendered_count = rendered_frame_count,
        complete = complete_sequence,
        fps = CRITIQUE_FPS,
        duration = duration_s,
        spp = config.samples_per_pixel,
        seed_salt = config.render_seed_salt,
        max_depth = config.max_depth,
        shutter_angle = config.shutter_angle_degrees,
        shutter_duration =
            f64::from(config.shutter_angle_degrees) / 360.0 / f64::from(CRITIQUE_FPS),
        seed_version = CRITIQUE_FRAME_SEED_SCHEDULE_VERSION,
        denoise_requested = config.denoise_previews,
        raw_profile = raw_profile,
        raw_sequence = raw_sequence_identity.to_hex(),
        preview_sequence = preview_sequence_identity.to_hex(),
        over_range = over_range_channels,
        gamut_mapped = gamut_mapped_pixels,
        requested_workers = config.render_workers,
        effective_workers = render.maximum_effective_workers,
        tile_width = config.tile_width,
        tile_height = config.tile_height,
        memory_limit = config.render_memory_limit_bytes,
        peak_memory = render.maximum_peak_memory_bytes,
        measured_frames = render.frames,
        setup_ns = render.setup_ns,
        traversal_ns = render.traversal_ns,
        compute_ns = render.tile_compute_ns,
        merge_ns = render.tile_merge_ns,
        publication_ns = render.publication_ns,
        idle_ns = render.idle_worker_ns,
        denoised_frames = denoise.applied_frames,
        denoise_pipeline = denoise_pipeline,
        denoise_bytes = denoise.maximum_retained_bytes,
        history_frames = denoise.maximum_history_frames,
        source_id = run.provenance.small_angle_oracle_source_id,
        model_authority = run.provenance.model_authority,
        physical_validation = run.provenance.physical_validation,
        diameter = specimen.diameter_m,
        thickness = specimen.thickness_m,
        mass = specimen.mass_kg,
        fillet = specimen.outer_fillet_radius_m,
        coarse_dt = run.parameters.timestep_s,
        fine_dt = refinement.fine.parameters.timestep_s,
        source_samples = run.samples.len(),
        source_duration = last_source_sample.time_s,
        tail_samples = trajectory_samples.len(),
        tail_duration = last_visual_sample.time_s,
        terminal = run.terminal,
        cutoff = run.parameters.validity_cutoff_theta_rad,
        refine_time = refinement.terminal_time_difference_s,
        refine_work = refinement.total_work_difference_j,
        rolling_mu = run
            .provenance
            .published_rolling_coefficient_mu
            .expect("source-bound run retains rolling coefficient"),
        rolling_work = last_source_sample.work.published_rolling_j,
        gas_work = last_source_sample.work.bildsten_boundary_layer_j,
        first_theta = first_visual_sample.qois.inclination_rad,
        first_precession = first_visual_sample.qois.precession_rad_per_s,
        first_spin = first_visual_sample.qois.spin_rad_per_s,
        last_theta = last_visual_sample.qois.inclination_rad,
        last_precession = last_visual_sample.qois.precession_rad_per_s,
        last_spin = last_visual_sample.qois.spin_rad_per_s,
        initial_energy = first_source_sample.energy_j,
        final_energy = last_source_sample.energy_j,
        energy_residual = run.energy_closure_residual_j,
        relative_residual = relative_energy_defect,
        trajectory_identity = trajectory_identity.to_hex(),
        audio_rate = SOUND_MASTER_SAMPLE_RATE_HZ,
        wav_identity = wav_identity.to_hex(),
        chirp_start = chirp_start_hz,
        chirp_end = chirp_end_hz,
        audio_bandwidth = CRITIQUE_DECLARED_AUDIO_SOURCE_BANDWIDTH_HZ,
        modal_identity = modal_parameter_set_identity.to_hex(),
        modal_disclosure = modal_disclosure,
        pre_master_peak = audio_pre_master_peak_fs,
        master_gain = audio_master_gain_db,
        initial_fade = INITIAL_FADE_SAMPLE_FRAMES,
        terminal_fade = TERMINAL_FADE_SAMPLE_FRAMES,
        spatialization = spatialization_json,
        mux = mux_json,
    )
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            value if value.is_control() => {
                use core::fmt::Write as _;
                let _ = write!(escaped, "\\u{:04x}", value as u32);
            }
            value => escaped.push(value),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    use super::*;

    fn with_test_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x4555_4c45_525f_5445,
                    kernel_id: 0x4649_5854,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    #[test]
    fn default_is_an_eight_second_practical_preview() {
        let config = CinematicFixtureConfig::default();
        config.validate().unwrap();
        assert_eq!(config.frames, 8 * CRITIQUE_FPS);
        assert_eq!(config.frame_window, CinematicFrameWindow::Full);
        assert_eq!(config.render_seed_salt, 0);
        assert_eq!((config.width, config.height), (320, 180));
        assert_eq!(config.shutter_angle_degrees, 180);
        assert!(config.render_workers > 0);
        assert!(config.denoise_previews);
        assert!(config.retain_full_aov_exr);
        assert!(config.spatialize_audio);
    }

    #[test]
    fn sound_contract_frame_range_is_enforced() {
        let mut config = CinematicFixtureConfig::default();
        config.frames = 191;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(_))
        ));
    }

    #[test]
    fn execution_and_shutter_bounds_are_enforced() {
        let mut config = CinematicFixtureConfig::default();
        config.shutter_angle_degrees = 361;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(_))
        ));
        config.shutter_angle_degrees = 180;
        config.render_workers = 0;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(_))
        ));
        config.render_workers = 1;
        config.tile_width = 0;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(_))
        ));
    }

    #[test]
    fn partial_frame_windows_are_bounded_and_cannot_masquerade_as_movies() {
        let mut config = CinematicFixtureConfig::default();
        config.mux_with_ffmpeg = false;
        config.frame_window = CinematicFrameWindow::Range {
            first_frame: 96,
            frame_count: 1,
        };
        assert_eq!(config.render_frame_range().unwrap(), 96..97);
        config.mux_with_ffmpeg = true;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(
                "partial frame windows cannot be muxed as complete movies"
            ))
        ));

        config.mux_with_ffmpeg = false;
        config.frame_window = CinematicFrameWindow::Range {
            first_frame: CRITIQUE_FRAMES,
            frame_count: 1,
        };
        assert!(config.validate().is_err());
        config.frame_window = CinematicFrameWindow::Range {
            first_frame: 0,
            frame_count: 0,
        };
        assert!(config.validate().is_err());
        config.frame_window = CinematicFrameWindow::Range {
            first_frame: u32::MAX,
            frame_count: 2,
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn absolute_frame_seed_schedule_is_replayable_and_frame_distinct() {
        let base = euler_scene_smoke_settings(320, 180).seed;
        assert_eq!(
            frame_render_seed(base, 0, 96),
            frame_render_seed(base, 0, 96)
        );
        assert_ne!(
            frame_render_seed(base, 0, 95),
            frame_render_seed(base, 0, 96)
        );
        assert_ne!(
            frame_render_seed(base, 0, 96),
            frame_render_seed(base, 0, 97)
        );
        assert_ne!(
            frame_render_seed(base ^ 1, 0, 96),
            frame_render_seed(base, 0, 96)
        );
        assert_ne!(
            frame_render_seed(base, 1, 96),
            frame_render_seed(base, 0, 96)
        );
    }

    #[test]
    fn boundary_reference_times_cover_finite_shutters() {
        let trajectory_start_s = 5.0e-5;
        let trajectory_end_s = 7.999_999_999_988_025;
        let half_shutter_s = 0.25 / f64::from(CRITIQUE_FPS);
        for frame in 0..CRITIQUE_FRAMES {
            let (frame_time_s, previous_time_s, next_time_s) =
                frame_reference_times(frame, CRITIQUE_FRAMES, trajectory_start_s, trajectory_end_s);
            assert!(trajectory_start_s <= previous_time_s);
            assert!(previous_time_s <= frame_time_s - half_shutter_s);
            assert!(frame_time_s + half_shutter_s <= next_time_s);
            assert!(next_time_s <= trajectory_end_s);
            if frame == 0 {
                assert_eq!(previous_time_s, trajectory_start_s);
            }
            if frame + 1 == CRITIQUE_FRAMES {
                assert_eq!(next_time_s, trajectory_end_s);
            }
        }
    }

    #[test]
    fn terminal_fade_preserves_prefix_and_reaches_zero() {
        let mut stems = vec![
            crate::ModalStemFrame {
                disc_fs: 1.0,
                glass_plate_fs: 2.0,
                base_assembly_fs: 3.0,
            };
            4
        ];
        apply_terminal_fade(&mut stems, 3).unwrap();
        assert_eq!(stems[0].disc_fs, 1.0);
        assert_eq!(stems[1].disc_fs, 1.0);
        assert_eq!(stems[2].disc_fs, 0.5);
        assert_eq!(stems[3], crate::ModalStemFrame::default());
    }

    #[test]
    fn initial_fade_starts_at_zero_and_preserves_suffix() {
        let mut stems = vec![
            crate::ModalStemFrame {
                disc_fs: 1.0,
                glass_plate_fs: 2.0,
                base_assembly_fs: 3.0,
            };
            4
        ];
        apply_initial_fade(&mut stems, 3).unwrap();
        assert_eq!(stems[0], crate::ModalStemFrame::default());
        assert_eq!(stems[1].disc_fs, 0.5);
        assert_eq!(stems[2].disc_fs, 1.0);
        assert_eq!(stems[3].disc_fs, 1.0);
    }

    #[test]
    fn spatial_listener_basis_is_finite_unit_and_orthogonal() {
        let listener = spatial_listener_pose();
        let forward_norm = listener
            .forward_unit
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let right_norm = listener
            .right_unit
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let dot = listener
            .forward_unit
            .iter()
            .zip(listener.right_unit)
            .map(|(forward, right)| forward * right)
            .sum::<f64>();
        assert!(listener.position_m.iter().all(|value| value.is_finite()));
        assert!((forward_norm - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert!((right_norm - 1.0).abs() <= 4.0 * f64::EPSILON);
        assert!(dot.abs() <= 4.0 * f64::EPSILON);
        let target_distance_m = CRITIQUE_CAMERA_TARGET_M
            .iter()
            .zip(listener.position_m)
            .map(|(target, listener)| (target - listener) * (target - listener))
            .sum::<f64>()
            .sqrt();
        assert!(CRITIQUE_SPATIAL_REFERENCE_DISTANCE_M > 0.0);
        assert!(CRITIQUE_SPATIAL_REFERENCE_DISTANCE_M < target_distance_m);
    }

    #[test]
    fn post_spatial_fade_preserves_length_and_publishes_exact_zero_endpoint() {
        let samples = vec![
            StereoSample {
                left_fs: 1.0,
                right_fs: -1.0,
            };
            5
        ];
        let faded = with_test_cx(|cx| apply_spatial_terminal_fade(&samples, 3, cx)).unwrap();
        assert_eq!(faded.len(), samples.len());
        assert_eq!(faded[0], samples[0]);
        assert_eq!(faded[1], samples[1]);
        assert_eq!(faded[2], samples[2]);
        assert_eq!(faded[3].left_fs, 0.5);
        assert_eq!(faded[3].right_fs, -0.5);
        assert_eq!(faded[4].left_fs.to_bits(), 0.0_f64.to_bits());
        assert_eq!(faded[4].right_fs, 0.0);
    }
}
