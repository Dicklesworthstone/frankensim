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
    AudioResamplingBudget, AudioResamplingChunk, AudioResamplingCrop, AudioResamplingModelInput,
    ContactModeShape, ContactParticipationPolicy, DerivedEulerQois,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerControlStream, EulerRenderTrajectoryArtifact,
    ListenerPose as SpatialListenerPose, ListenerPoseTrack, MAX_AUDIO_MASTER_GAIN_DB,
    MicrophoneDirectivity, ModalPresetAuthority, ModalStemFrame, ModalSynthesisBudget,
    ModalSynthesisModel, ModalSynthesisModelInput, ModeContactParticipationRule,
    OfflineSpatializer, RenderNormalForceSampling, RenderTrajectory, RenderTrajectoryCodecBudget,
    RepresentativeDiscMaterial, SoundWavArtifact, SourcePositionTrack, SpatialAudioAuthority,
    SpatialAudioBudget, SpatialAudioConfig, SpatialAudioOutput, SpatialAudioRenderInput,
    SpatialAudioSource, SpatialDelayPolicy, SpatialMonoSignal, SpatialOutputHorizon,
    SpatialStemComponent, StemGainPan, StereoSample, WavMetadata, WavSampleEncoding, measure_audio,
    mix_dry_modal_stems,
    reduced_decay::{
        ReducedDecayError, ReducedDecayRun, RefinementEvidence, Thorne2026SteelGlassBenchmark,
        thorne_2026_refinement_evidence,
    },
    render_motion_bridge::EulerRenderMotionBridge,
    render_scene_bridge::{
        EulerCinematicScene, EulerDiscMaterialStateBinding, EulerEnvironmentStyle,
        EulerFrameRequest, EulerMaterialStyle, EulerPreviewMeshReceipt, EulerRectLightSpec,
        EulerSceneConfig, EulerStudioEnvironmentSpec, EulerSupportSurfaceSpec,
        EulerTessellationConfig, MAX_EULER_ARC_SUBDIVISIONS, MAX_EULER_AZIMUTHAL_SEGMENTS,
        euler_scene_smoke_settings,
    },
    representative_modal_preset,
    timeline_resampling::{EventEvaluationSide, ExposureEventPolicy},
};
use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_evidence::{
    ValidityDomain,
    cinematic::{CinematicClock, CinematicClockDomain, SoundAuthority},
    cinematic_config::{CinematicComponentRef, CinematicComponentRole},
    cinematic_sound::{
        ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SOUND_SYNTHESIS_SCHEMA_VERSION,
        SoundAmplitudeReference, SoundChannelLayout, SoundExcitationChannel,
        SoundExcitationControl, SoundModalComponent, SoundModelAssumption, SoundRoomResponse,
        SoundSynthesisConfig, SoundSynthesisInput, SoundTerminalPolicy, SoundTrajectoryDisposition,
    },
};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_exec::{Cx, RunId};
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_img::{
    Channel, CinematicColorConfig, CinematicColorLimits, ExrAttribute, PixelType, PngColor,
    PreviewDither, TEMPORAL_DENOISE_PIPELINE_VERSION, TemporalDenoiseConfig, TemporalDenoiseInput,
    TemporalDenoiseLimits, TemporalDenoisedFrame, TemporalFrameBoundary, temporal_denoise_rgb,
    transform_cinematic_preview, write_exr_with_attributes, write_png16,
};
use fs_math::det;
use fs_material::{
    gas::{GasSpec, GasState},
    state_point::{
        MaterialPropertySelection, VISIBLE_COMPLEX_IOR_ETA_PROPERTIES,
        VISIBLE_COMPLEX_IOR_K_PROPERTIES, resolve_isotropic_elastic_state_point,
        resolve_visible_conductor_state_point,
    },
    visco::RayleighDamping,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_mbd::Vec3;
use fs_modal::SliceOptions;
use fs_qty::{Density, Dims, Pressure};
use fs_render::{
    aov::{
        CinematicAovConfig, CinematicAovLimits, CinematicAovProfile, CinematicAovProvenance,
        cinematic_render_semantics_versions,
    },
    camera::{AnimatedCamera, Aperture, CameraProjection, CutSide, PhysicalCamera},
    conductor::{ConductorOptics, ConductorSurface},
    dielectric::{DielectricGlass, DielectricSurface},
    motion::{
        NormalizedShutterTime, ShotTimeBounds, ShutterConvention, ShutterDistribution,
        ShutterInterval,
    },
    tracer::{
        ADAPTIVE_SAMPLING_SEMANTICS_VERSION, AdaptiveFilm, AdaptiveSamplingConfig,
        MAX_RENDER_TILE_EDGE, MAX_RENDER_WORKERS, RenderExecutionConfig, RenderExecutionReport,
        RenderWorkerPool, Settings, film_to_exr,
    },
};
use fs_solid::TetAssemblyBudget;

/// Fixed master frame rate used by the cinematic sound contract.
pub const CRITIQUE_FPS: u32 = 24;
/// Minimum admitted cinematic duration: 192 frames = 8 seconds.
pub const CRITIQUE_FRAMES: u32 = 192;
/// Five-millisecond deterministic taper applied at the censored soundtrack end.
const TERMINAL_FADE_SAMPLE_FRAMES: u32 = 240;
/// Twenty-millisecond presentation fade at the published clip onset.
const INITIAL_FADE_SAMPLE_FRAMES: u32 = 960;
/// One exact video frame of real source history used to warm the modal state.
const AUDIO_PREROLL_VIDEO_FRAMES: u32 = 1;
/// Exact 48 kHz samples in the source-bound modal warm start.
const AUDIO_PREROLL_SAMPLE_FRAMES: u64 = 2_000;
/// The production film compares the original benchmark's quarter-step against
/// its eighth-step solution. The earlier dt/dt2 pair did not satisfy the
/// preregistered mechanics-to-audio drive convergence gate.
const CINEMATIC_MECHANICS_COARSE_REFINEMENT_FACTOR: u32 = 4;
/// Frozen display exposure shared by in-process and offline critique previews
/// and declared verbatim in the fixture manifest.
pub const CRITIQUE_EXPOSURE_EV: i32 = 0;

/// Exact display finishing used by both in-process and offline previews.
#[must_use]
pub fn critique_color_config() -> CinematicColorConfig {
    let mut color = CinematicColorConfig::reference_srgb_16();
    color.exposure_ev = CRITIQUE_EXPOSURE_EV;
    color.dither = PreviewDither::Disabled;
    color
}
const AUDIO_PREROLL_POLICY_ID: &str =
    "source-bound-one-video-frame-continuous-fir-and-modal-preroll-v3";
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
const CRITIQUE_CAMERA_EYE_M: [f64; 3] = [0.18, -0.24, 0.12];
/// World-space look target shared by picture and spatial sound.
const CRITIQUE_CAMERA_TARGET_M: [f64; 3] = [0.0, 0.0, 0.008];
/// Artistic near-field reference distance for the desk-scale stereo preview.
///
/// The spatializer's distance gain is `reference / max(distance, reference)`.
/// Keeping this just below the camera-to-subject distance preserves a modest
/// distance cue without needlessly spending roughly 18 dB of the bounded
/// mechanics-to-digital mastering range on a five-centimetre reference.
const CRITIQUE_SPATIAL_REFERENCE_DISTANCE_M: f64 = 0.28;

// These output-space gates are fixed before evaluating the dt/dt/2 pair. They
// test whether refinement materially changes the quantities consumed by this
// exact picture-and-sound fixture; they are not experimental validation or an
// asymptotic convergence-order claim.
const OUTPUT_CONVERGENCE_RELATIVE_IMPULSE_LIMIT: f64 = 1.0e-3;
const OUTPUT_CONVERGENCE_COM_LIMIT_M: f64 = 1.0e-5;
const OUTPUT_CONVERGENCE_ORIENTATION_LIMIT_RAD: f64 = 1.0e-3;
const OUTPUT_CONVERGENCE_CHIRP_ABSOLUTE_LIMIT_HZ: f64 = 0.1;
const OUTPUT_CONVERGENCE_CHIRP_RELATIVE_LIMIT: f64 = 1.0e-3;
const OUTPUT_CONVERGENCE_TERMINAL_TIME_ABSOLUTE_LIMIT_S: f64 = 1.0e-4;
const OUTPUT_CONVERGENCE_TERMINAL_TIME_RELATIVE_LIMIT: f64 = 2.0e-5;
const OUTPUT_CONVERGENCE_IMPULSE_ROUNDOFF_PER_INTERVAL: f64 = 64.0;

// These audio-output gates were fixed before evaluating the dt/dt/2 pair.
// Drive is compared before modal synthesis; raw component stems are compared
// before fades, spatialization, mixing, or mastering. They establish only one
// encoded-model output-consistency pair, never acoustic or experimental truth.
const AUDIO_CONVERGENCE_DRIVE_NRMSE_LIMIT: f64 = 1.0e-3;
const AUDIO_CONVERGENCE_DRIVE_NORMALIZED_PEAK_LIMIT: f64 = 5.0e-3;
const AUDIO_CONVERGENCE_STEM_NRMSE_LIMIT: f64 = 1.0e-2;
const AUDIO_CONVERGENCE_STEM_NORMALIZED_PEAK_LIMIT: f64 = 2.0e-2;
const AUDIO_CONVERGENCE_DRIVE_NORMALIZATION_FRACTION: f64 = 1.0e-12;
const AUDIO_CONVERGENCE_STEM_FLOOR_FS: f64 = 1.0e-12;

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

/// Explicit deterministic raw-estimator stopping policy for production renders.
///
/// The thresholds are the per-channel raw-XYZ dispersion proxies defined by
/// [`AdaptiveSamplingConfig`]. They are not confidence intervals, perceptual
/// image-error bounds, or denoiser controls. The hard maximum remains an
/// explicit part of the policy so an adaptive request can never become
/// unbounded work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CinematicAdaptiveSamplingConfig {
    /// First per-pixel sample count at which convergence may be declared.
    pub minimum_samples_per_pixel: u32,
    /// Hard per-pixel path ceiling.
    pub maximum_samples_per_pixel: u32,
    /// Number of additional samples between deterministic decisions.
    pub decision_batch_samples: u32,
    /// Per-channel absolute raw-XYZ dispersion allowance.
    pub absolute_error: f64,
    /// Per-channel relative raw-XYZ dispersion allowance.
    pub relative_error: f64,
    /// Lower raw-XYZ scale used by the relative term for dark channels.
    pub dark_floor: f64,
}

impl CinematicAdaptiveSamplingConfig {
    fn policy(self) -> Result<AdaptiveSamplingConfig, CinematicFixtureError> {
        AdaptiveSamplingConfig::try_new(
            self.minimum_samples_per_pixel,
            self.decision_batch_samples,
            self.absolute_error,
            self.relative_error,
            self.dark_floor,
        )
        .map_err(pipeline)
    }

    fn canonical(self) -> Result<Self, CinematicFixtureError> {
        let policy = self.policy()?;
        Ok(Self {
            minimum_samples_per_pixel: policy.minimum_samples(),
            maximum_samples_per_pixel: self.maximum_samples_per_pixel,
            decision_batch_samples: policy.batch_samples(),
            absolute_error: policy.absolute_error(),
            relative_error: policy.relative_error(),
            dark_floor: policy.dark_floor(),
        })
    }
}

/// Material-state inputs for the physical disc used by picture and sound.
///
/// These are numerical constitutive data, not a material-name selector. A
/// caller can replace every value and provenance string; the same resolved
/// state then drives density/stiffness, structural modes, damping, acoustic
/// radiation, spectral Fresnel response, and surface roughness. The default
/// film instance uses the literature specimen's mass-derived density, a
/// declared room-temperature steel elastic estimate, and room-temperature
/// elemental-iron optical measurements as an explicitly imperfect proxy for
/// the reported but composition-unspecified steel disc.
#[derive(Clone, Debug, PartialEq)]
pub struct CinematicDiscMaterialConfig {
    /// Absolute material/environment query temperature [K].
    pub temperature_k: f64,
    /// Isotropic tangent Young's modulus [Pa].
    pub young_modulus_pa: f64,
    /// Isotropic Poisson ratio.
    pub poisson_ratio: f64,
    /// Mass-proportional Rayleigh damping coefficient [1/s].
    pub rayleigh_alpha_per_s: f64,
    /// Stiffness-proportional Rayleigh damping coefficient [s].
    pub rayleigh_beta_s: f64,
    /// Visible complex-index real parts on FrankenSim's canonical 9-knot grid.
    pub optical_eta: [f64; 9],
    /// Visible complex-index extinction coefficients on the same grid.
    pub optical_k: [f64; 9],
    /// Isotropic GGX microfacet roughness for the explicit surface state.
    pub surface_roughness_alpha: f64,
    /// Identity-only chemistry/specimen description; never dispatches physics.
    pub material_label: String,
    /// Process/surface history identifier; never dispatches physics.
    pub process_label: String,
    /// Citation or explicit estimate disclosure for elastic coefficients.
    pub elastic_source: String,
    /// Citation for the complex optical constants.
    pub optical_source: String,
    /// Citation or explicit estimate disclosure for damping coefficients.
    pub damping_source: String,
}

impl Default for CinematicDiscMaterialConfig {
    fn default() -> Self {
        Self {
            temperature_k: 293.15,
            young_modulus_pa: 200.0e9,
            poisson_ratio: 0.29,
            // A deliberately disclosed, low-loss room-temperature estimate.
            // These are input data, not an outcome-tuned frequency or decay.
            rayleigh_alpha_per_s: 0.35,
            rayleigh_beta_s: 6.0e-8,
            // Linear interpolation of Johnson & Christy's room-temperature
            // elemental-iron n,k table onto 380:50:780 nm. This is visibly and
            // provenance-wise preferable to an artistic steel preset, but it
            // remains a proxy because the paper does not report alloy/finish.
            optical_eta: [
                2.112_307_692_308,
                2.472_777_777_778,
                2.695_2,
                2.888_928_571_429,
                2.940_606_060_606,
                2.892_380_952_381,
                2.892,
                2.865,
                2.895_846_153_846,
            ],
            optical_k: [
                2.494_615_384_615,
                2.706_666_666_667,
                2.841_6,
                2.916_428_571_429,
                2.986_363_636_364,
                3.065_476_190_476,
                3.142,
                3.235,
                3.320_615_384_615,
            ],
            surface_roughness_alpha: 0.12,
            material_label: "Thorne-2026 reported steel; composition unspecified".to_owned(),
            process_label: "machined filleted disc; finish unmeasured".to_owned(),
            elastic_source:
                "FrankenSim room-temperature isotropic steel estimate; not measured on specimen"
                    .to_owned(),
            optical_source: "P. B. Johnson and R. W. Christy, Phys. Rev. B 9, 5056-5070 (1974), DOI 10.1103/PhysRevB.9.5056; elemental-iron proxy, linearly interpolated"
                .to_owned(),
            damping_source:
                "FrankenSim low-loss Rayleigh estimate; not measured on specimen".to_owned(),
        }
    }
}

/// Bounded settings for one watchable critique artifact.
#[derive(Clone, Debug, PartialEq)]
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
    /// Optional deterministic variance-targeted raw sampling policy.
    ///
    /// Beauty, retained AOVs, and denoising guides all observe each pixel's
    /// exact accepted sample prefix. The raw EXR and manifest retain the
    /// authoritative per-pixel counts; denoising remains a biased display
    /// derivative and cannot affect the stopping decision.
    pub adaptive_sampling: Option<CinematicAdaptiveSamplingConfig>,
    /// Caller-selected independent scramble salt for replicated raw renders.
    pub render_seed_salt: u64,
    /// Maximum path depth, including dielectric traversal.
    pub max_depth: u32,
    /// Render-only samples around the disc's complete revolution.
    pub azimuthal_segments: u32,
    /// Render-only equal-angle chords used for each circular meridian arc.
    pub arc_subdivisions_per_arc: u32,
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
    /// Persist the full FinalDiagnostic AOV EXR rather than three-channel raw
    /// beauty. Denoising still retains FinalDiagnostic guides in memory when
    /// disk retention is disabled, so this switch cannot silently weaken edge
    /// identity.
    pub retain_full_aov_exr: bool,
    /// Whether mechanics-derived dry stems use the bounded spatial-audio path.
    pub spatialize_audio: bool,
    /// One parameterized physical material state shared by picture and sound.
    pub disc_material: CinematicDiscMaterialConfig,
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
            adaptive_sampling: None,
            render_seed_salt: 0,
            max_depth: 6,
            // The conservative 4K look-development tier: on the canonical
            // specimen 512 x 64 retains sub-micrometre meridian and azimuthal
            // chord errors. The renderer currently uses geometric facet
            // normals, so the more aggressive 256 x 32 tier is deliberately
            // not the production default.
            azimuthal_segments: 512,
            arc_subdivisions_per_arc: 64,
            shutter_angle_degrees: 180,
            render_workers: default_render_workers(),
            tile_width: default_render_tile_edge(),
            tile_height: default_render_tile_edge(),
            render_memory_limit_bytes: 4 * 1024 * 1024 * 1024,
            denoise_previews: true,
            retain_full_aov_exr: true,
            spatialize_audio: true,
            disc_material: CinematicDiscMaterialConfig::default(),
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
        if let Some(adaptive) = self.adaptive_sampling {
            if adaptive.minimum_samples_per_pixel < 2 || adaptive.minimum_samples_per_pixel > 4_096
            {
                return Err(CinematicFixtureError::InvalidConfig(
                    "adaptive minimum samples per pixel must be in 2..=4096",
                ));
            }
            if adaptive.maximum_samples_per_pixel < adaptive.minimum_samples_per_pixel
                || adaptive.maximum_samples_per_pixel > 4_096
            {
                return Err(CinematicFixtureError::InvalidConfig(
                    "adaptive maximum samples per pixel must be in minimum..=4096",
                ));
            }
            if adaptive.decision_batch_samples == 0 {
                return Err(CinematicFixtureError::InvalidConfig(
                    "adaptive decision batch samples must be nonzero",
                ));
            }
            if [
                adaptive.absolute_error,
                adaptive.relative_error,
                adaptive.dark_floor,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value < 0.0)
            {
                return Err(CinematicFixtureError::InvalidConfig(
                    "adaptive error controls must be finite and nonnegative",
                ));
            }
            adaptive.policy()?;
        }
        if self.max_depth == 0 || self.max_depth > 64 {
            return Err(CinematicFixtureError::InvalidConfig(
                "max_depth must be in 1..=64",
            ));
        }
        if !(8..=MAX_EULER_AZIMUTHAL_SEGMENTS).contains(&self.azimuthal_segments) {
            return Err(CinematicFixtureError::InvalidConfig(
                "azimuthal_segments must be in 8..=4096",
            ));
        }
        if !(1..=MAX_EULER_ARC_SUBDIVISIONS).contains(&self.arc_subdivisions_per_arc) {
            return Err(CinematicFixtureError::InvalidConfig(
                "arc_subdivisions_per_arc must be in 1..=1024",
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
        let material = &self.disc_material;
        if !(material.temperature_k.is_finite()
            && material.temperature_k > 0.0
            && material.young_modulus_pa.is_finite()
            && material.young_modulus_pa > 0.0
            && material.poisson_ratio.is_finite()
            && material.poisson_ratio > -1.0
            && material.poisson_ratio < 0.5
            && material.rayleigh_alpha_per_s.is_finite()
            && material.rayleigh_alpha_per_s >= 0.0
            && material.rayleigh_beta_s.is_finite()
            && material.rayleigh_beta_s >= 0.0
            && material.surface_roughness_alpha.is_finite()
            && material.surface_roughness_alpha >= 1.0e-4
            && material.surface_roughness_alpha <= 1.0)
        {
            return Err(CinematicFixtureError::InvalidConfig(
                "disc material mechanical, damping, or surface scalar is outside its physical domain",
            ));
        }
        if material
            .optical_eta
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || material
                .optical_k
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(CinematicFixtureError::InvalidConfig(
                "disc material optical eta/k table is outside its physical domain",
            ));
        }
        if [
            &material.material_label,
            &material.process_label,
            &material.elastic_source,
            &material.optical_source,
            &material.damping_source,
        ]
        .iter()
        .any(|value| value.trim().is_empty())
        {
            return Err(CinematicFixtureError::InvalidConfig(
                "disc material identity and provenance strings must be nonblank",
            ));
        }
        Ok(())
    }

    fn render_sample_ceiling(&self) -> u32 {
        self.adaptive_sampling
            .map_or(self.samples_per_pixel, |policy| {
                policy.maximum_samples_per_pixel
            })
    }

    fn adaptive_policy(&self) -> Result<Option<AdaptiveSamplingConfig>, CinematicFixtureError> {
        self.adaptive_sampling
            .map(CinematicAdaptiveSamplingConfig::policy)
            .transpose()
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

const fn default_render_tile_edge() -> u32 {
    // Native M4 profiling found 8 x 8 to be the smallest tier before tile
    // overhead erased the load-balance gain. Apply that Apple-aarch64 family
    // policy to M-series hosts while keeping other architectures unchanged;
    // only the measured M4 host carries performance evidence, not M1-M3 or M5.
    if cfg!(all(target_arch = "aarch64", target_os = "macos")) {
        8
    } else {
        32
    }
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
    warm_start_source_identity: ContentHash,
    warm_start_checkpoint_identity: ContentHash,
    source_sound_configuration_identity: ContentHash,
    published_trajectory_identity: ContentHash,
    crop_resampler_identity: ContentHash,
    crop_first_source_audio_frame: u64,
    crop_end_source_audio_frame: u64,
    convergence: FixtureAudioConvergenceEvidence,
    spatialization: Option<FixtureSpatialAudioEvidence>,
}

struct FixturePhysicalDiscState {
    specimen: crate::specimen::ResolvedElasticDiscProfile,
    render_binding: EulerDiscMaterialStateBinding,
    damping_model_identity: ContentHash,
}

fn resolve_fixture_physical_disc(
    reference_profile: &crate::specimen::ResolvedDiscProfile,
    config: &CinematicDiscMaterialConfig,
    cx: &Cx<'_>,
) -> Result<FixturePhysicalDiscState, CinematicFixtureError> {
    let validity = ValidityDomain::unconstrained().with(
        "T",
        config.temperature_k,
        config.temperature_k,
    );
    let mut claims = ClaimSet::new();
    let mut insert_scalar = |name: &str,
                             dims: Dims,
                             value: f64,
                             source: &str|
     -> Result<(), CinematicFixtureError> {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: source.to_owned(),
                    // The fixture's numeric declarations and the CC0 optical
                    // database transcription are redistributable inputs. The
                    // original paper remains the scientific citation.
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .map(|_| ())
            .map_err(pipeline)
    };
    insert_scalar(
        "density",
        Density::DIMS,
        reference_profile.density_kg_per_m3,
        "Thorne 2026 specimen mass and exact FrankenSim profile volume; mass-derived density",
    )?;
    insert_scalar(
        "young_modulus",
        Pressure::DIMS,
        config.young_modulus_pa,
        &config.elastic_source,
    )?;
    insert_scalar(
        "poisson_ratio",
        Dims::NONE,
        config.poisson_ratio,
        &config.elastic_source,
    )?;
    for ((name, value), source) in VISIBLE_COMPLEX_IOR_ETA_PROPERTIES
        .iter()
        .zip(config.optical_eta)
        .zip(core::iter::repeat(config.optical_source.as_str()))
    {
        insert_scalar(name, Dims::NONE, value, source)?;
    }
    for ((name, value), source) in VISIBLE_COMPLEX_IOR_K_PROPERTIES
        .iter()
        .zip(config.optical_k)
        .zip(core::iter::repeat(config.optical_source.as_str()))
    {
        insert_scalar(name, Dims::NONE, value, source)?;
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: config.material_label.clone(),
            phase: "solid".to_owned(),
            process: config.process_label.clone(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .map_err(pipeline)?;
    let point = QueryPoint::new()
        .with("T", config.temperature_k)
        .map_err(pipeline)?;
    let elastic = resolve_isotropic_elastic_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(pipeline)?;
    let optical = resolve_visible_conductor_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .map_err(pipeline)?;
    let specimen = reference_profile
        .spec
        .resolve_with_isotropic_elastic_state(&elastic, cx)
        .map_err(pipeline)?;
    if specimen.profile.content_identities().profile
        != reference_profile.content_identities().profile
    {
        return Err(CinematicFixtureError::Pipeline(
            "resolved elastic material changed the mechanics geometry/density profile".into(),
        ));
    }
    let mut surface =
        DomainHasher::new("org.frankensim.euler-critique.disc-surface-state.v1");
    surface.update(specimen.material_state_identity.as_bytes());
    surface.update(config.process_label.as_bytes());
    surface.update(&config.surface_roughness_alpha.to_bits().to_le_bytes());
    let render_binding = EulerDiscMaterialStateBinding::try_conductor_elastic(
        &elastic,
        &optical,
        config.surface_roughness_alpha,
        surface.finalize(),
    )
    .map_err(pipeline)?;
    let mut damping = DomainHasher::new("org.frankensim.euler-critique.rayleigh-damping.v1");
    damping.update(specimen.material_state_identity.as_bytes());
    damping.update(&config.rayleigh_alpha_per_s.to_bits().to_le_bytes());
    damping.update(&config.rayleigh_beta_s.to_bits().to_le_bytes());
    damping.update(config.damping_source.as_bytes());
    Ok(FixturePhysicalDiscState {
        specimen,
        render_binding,
        damping_model_identity: damping.finalize(),
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FixtureAudioConvergenceEvidence {
    full_audio_frame_count: u64,
    published_audio_frame_count: u64,
    mode_count: usize,
    drive_normalization_floor_n: f64,
    localized_drive_nrmse: f64,
    localized_drive_normalized_peak: f64,
    distributed_drive_nrmse: f64,
    distributed_drive_normalized_peak: f64,
    maximum_drive_nrmse: f64,
    maximum_drive_normalized_peak: f64,
    cropped_stem_nrmse: f64,
    cropped_stem_normalized_peak: f64,
    crop_first_source_audio_frame: u64,
    crop_end_source_audio_frame: u64,
    coarse_crop_identity: ContentHash,
    fine_crop_identity: ContentHash,
    coarse_crop_binding_identity: ContentHash,
    fine_crop_binding_identity: ContentHash,
}

impl FixtureAudioConvergenceEvidence {
    fn diagnostics(self) -> String {
        format!(
            concat!(
                "stage=audio timestep_convergence source=fine-dt2 full_frames={} ",
                "published_frames={} modes={} drive_nrmse={:.6e} ",
                "drive_peak_rel={:.6e} stem_nrmse={:.6e} stem_peak_rel={:.6e}"
            ),
            self.full_audio_frame_count,
            self.published_audio_frame_count,
            self.mode_count,
            self.maximum_drive_nrmse,
            self.maximum_drive_normalized_peak,
            self.cropped_stem_nrmse,
            self.cropped_stem_normalized_peak,
        )
    }

    fn manifest_json(self) -> String {
        format!(
            concat!(
                "{{\"comparison\":\"coarse-dt versus fine-dt2 on one shared continuous 48 kHz preroll clock; fine-dt2 alone is published\",",
                "\"full_audio_frame_count\":{},\"published_audio_frame_count\":{},\"mode_count\":{},",
                "\"drive\":{{\"normalization_floor_n\":{:.17e},",
                "\"localized_nrmse\":{:.17e},\"localized_normalized_peak\":{:.17e},",
                "\"distributed_nrmse\":{:.17e},\"distributed_normalized_peak\":{:.17e},",
                "\"maximum_nrmse\":{:.17e},\"nrmse_limit\":{:.17e},",
                "\"maximum_normalized_peak\":{:.17e},\"normalized_peak_limit\":{:.17e}}},",
                "\"raw_cropped_stems\":{{\"normalization_floor_fs\":{:.17e},",
                "\"nrmse\":{:.17e},\"nrmse_limit\":{:.17e},",
                "\"normalized_peak\":{:.17e},\"normalized_peak_limit\":{:.17e}}},",
                "\"crop_source_range\":{{\"first_audio_frame\":{},\"end_audio_frame\":{}}},",
                "\"identities\":{{\"coarse_crop\":\"{}\",\"fine_crop\":\"{}\",",
                "\"coarse_crop_binding\":\"{}\",\"fine_crop_binding\":\"{}\"}},",
                "\"artistic_localized_impulses\":\"exact all-zero in both members\",",
                "\"comparison_stage\":\"pre-fade, pre-spatialization, pre-mix, pre-master\",",
                "\"claim\":\"one fixed-gate output-consistency pair for the encoded reduced mechanics and declared uncalibrated modal model; not asymptotic convergence, contact-law validation, acoustic calibration, psychoacoustic equivalence, or experiment\"}}"
            ),
            self.full_audio_frame_count,
            self.published_audio_frame_count,
            self.mode_count,
            self.drive_normalization_floor_n,
            self.localized_drive_nrmse,
            self.localized_drive_normalized_peak,
            self.distributed_drive_nrmse,
            self.distributed_drive_normalized_peak,
            self.maximum_drive_nrmse,
            AUDIO_CONVERGENCE_DRIVE_NRMSE_LIMIT,
            self.maximum_drive_normalized_peak,
            AUDIO_CONVERGENCE_DRIVE_NORMALIZED_PEAK_LIMIT,
            AUDIO_CONVERGENCE_STEM_FLOOR_FS,
            self.cropped_stem_nrmse,
            AUDIO_CONVERGENCE_STEM_NRMSE_LIMIT,
            self.cropped_stem_normalized_peak,
            AUDIO_CONVERGENCE_STEM_NORMALIZED_PEAK_LIMIT,
            self.crop_first_source_audio_frame,
            self.crop_end_source_audio_frame,
            self.coarse_crop_identity.to_hex(),
            self.fine_crop_identity.to_hex(),
            self.coarse_crop_binding_identity.to_hex(),
            self.fine_crop_binding_identity.to_hex(),
        )
    }
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

/// Compact copy of the scene builder's render-only preview-mesh receipt.
///
/// Keeping this separate from [`CinematicFixtureConfig`] ensures the manifest
/// records what the scene builder actually admitted and emitted, rather than
/// merely echoing the requested tessellation controls.
#[derive(Clone, Copy, Debug, PartialEq)]
struct FixturePreviewMeshEvidence {
    azimuthal_segments: u32,
    arc_subdivisions_per_arc: u32,
    vertex_count: usize,
    triangle_count: usize,
    maximum_meridian_chord_error_m: f64,
    maximum_azimuthal_chord_error_m: f64,
}

impl FixturePreviewMeshEvidence {
    fn from_receipt(receipt: EulerPreviewMeshReceipt) -> Self {
        Self {
            azimuthal_segments: receipt.azimuthal_segments,
            arc_subdivisions_per_arc: receipt.arc_subdivisions_per_arc,
            vertex_count: receipt.vertex_count,
            triangle_count: receipt.triangle_count,
            maximum_meridian_chord_error_m: receipt.maximum_meridian_chord_error_m,
            maximum_azimuthal_chord_error_m: receipt.maximum_azimuthal_chord_error_m,
        }
    }

    fn manifest_json(self) -> String {
        format!(
            concat!(
                "{{\"authority\":\"render-only chordal approximation\",",
                "\"azimuthal_segments\":{},\"arc_subdivisions_per_arc\":{},",
                "\"vertex_count\":{},\"triangle_count\":{},",
                "\"maximum_meridian_chord_error_m\":{:.17e},",
                "\"maximum_azimuthal_chord_error_m\":{:.17e}}}"
            ),
            self.azimuthal_segments,
            self.arc_subdivisions_per_arc,
            self.vertex_count,
            self.triangle_count,
            self.maximum_meridian_chord_error_m,
            self.maximum_azimuthal_chord_error_m,
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FixtureAdaptiveFrameEvidence {
    frame: u32,
    pixels: u64,
    minimum_samples: u32,
    maximum_samples: u32,
    total_samples: u64,
    converged_pixels: u64,
    maximum_sample_pixels: u64,
    sample_count_identity: ContentHash,
}

/// Exact summaries derived from the published adaptive films.
///
/// The raw EXR retains every per-pixel count in its `samples` channel. This
/// structure keeps compact manifest evidence plus content identities for those
/// row-major maps; it never upgrades the raw dispersion proxy into an image
/// error certificate.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct FixtureAdaptiveSamplingEvidence {
    frames: Vec<FixtureAdaptiveFrameEvidence>,
    pixels: u64,
    minimum_samples: u32,
    maximum_samples: u32,
    total_samples: u64,
    converged_pixels: u64,
    maximum_sample_pixels: u64,
}

impl FixtureAdaptiveSamplingEvidence {
    fn observe(
        &mut self,
        frame: u32,
        film: &AdaptiveFilm,
        policy: CinematicAdaptiveSamplingConfig,
    ) -> Result<ContentHash, CinematicFixtureError> {
        if self
            .frames
            .last()
            .is_some_and(|previous| previous.frame >= frame)
        {
            return Err(CinematicFixtureError::Pipeline(
                "adaptive frame evidence was not observed in strictly increasing order".into(),
            ));
        }
        let admitted_policy = policy.policy()?;
        if film.maximum_samples() != policy.maximum_samples_per_pixel
            || film.policy() != admitted_policy
            || film.semantics_version() != ADAPTIVE_SAMPLING_SEMANTICS_VERSION
        {
            return Err(CinematicFixtureError::Pipeline(format!(
                "adaptive frame {frame} estimator provenance does not match its fixture policy"
            )));
        }
        let summary = film.summary();
        let expected_pixels = u64::from(film.width())
            .checked_mul(u64::from(film.height()))
            .ok_or_else(|| {
                CinematicFixtureError::Pipeline(
                    "adaptive frame pixel-count product overflowed".into(),
                )
            })?;
        if summary.pixels != expected_pixels
            || summary.converged_pixels + summary.maximum_sample_pixels != summary.pixels
            || summary.minimum_samples < policy.minimum_samples_per_pixel
            || summary.maximum_samples > policy.maximum_samples_per_pixel
            || film.sample_counts().len() != usize::try_from(expected_pixels).unwrap_or(usize::MAX)
        {
            return Err(CinematicFixtureError::Pipeline(format!(
                concat!(
                    "adaptive frame {} published inconsistent sample-count evidence: ",
                    "pixels={} expected={} min={} max={} converged={} ceiling={}"
                ),
                frame,
                summary.pixels,
                expected_pixels,
                summary.minimum_samples,
                summary.maximum_samples,
                summary.converged_pixels,
                summary.maximum_sample_pixels,
            )));
        }
        let mut hasher =
            DomainHasher::new("org.frankensim.euler-critique.adaptive-sample-count-frame.v1");
        hasher.update(&frame.to_le_bytes());
        hasher.update(&film.width().to_le_bytes());
        hasher.update(&film.height().to_le_bytes());
        for samples in film.sample_counts() {
            hasher.update(&samples.to_le_bytes());
        }
        let sample_count_identity = hasher.finalize();
        let evidence = FixtureAdaptiveFrameEvidence {
            frame,
            pixels: summary.pixels,
            minimum_samples: summary.minimum_samples,
            maximum_samples: summary.maximum_samples,
            total_samples: summary.total_samples,
            converged_pixels: summary.converged_pixels,
            maximum_sample_pixels: summary.maximum_sample_pixels,
            sample_count_identity,
        };
        self.pixels = self.pixels.checked_add(evidence.pixels).ok_or_else(|| {
            CinematicFixtureError::Pipeline("adaptive rendered-pixel count overflowed".into())
        })?;
        self.total_samples = self
            .total_samples
            .checked_add(evidence.total_samples)
            .ok_or_else(|| {
                CinematicFixtureError::Pipeline("adaptive total path count overflowed".into())
            })?;
        self.converged_pixels = self
            .converged_pixels
            .checked_add(evidence.converged_pixels)
            .ok_or_else(|| {
                CinematicFixtureError::Pipeline("adaptive converged-pixel count overflowed".into())
            })?;
        self.maximum_sample_pixels = self
            .maximum_sample_pixels
            .checked_add(evidence.maximum_sample_pixels)
            .ok_or_else(|| {
                CinematicFixtureError::Pipeline(
                    "adaptive maximum-sample pixel count overflowed".into(),
                )
            })?;
        self.minimum_samples = if self.frames.is_empty() {
            evidence.minimum_samples
        } else {
            self.minimum_samples.min(evidence.minimum_samples)
        };
        self.maximum_samples = self.maximum_samples.max(evidence.maximum_samples);
        self.frames.push(evidence);
        Ok(sample_count_identity)
    }

    fn validate_complete(
        &self,
        config: &CinematicFixtureConfig,
        rendered_frames: &core::ops::Range<u32>,
    ) -> Result<(), CinematicFixtureError> {
        let Some(policy) = config.adaptive_sampling else {
            if self.frames.is_empty() {
                return Ok(());
            }
            return Err(CinematicFixtureError::Pipeline(
                "uniform render unexpectedly retained adaptive sample evidence".into(),
            ));
        };
        if self.frames.len() != rendered_frames.len()
            || self.frames.first().map(|frame| frame.frame) != Some(rendered_frames.start)
            || self.frames.last().map(|frame| frame.frame) != Some(rendered_frames.end - 1)
            || self.converged_pixels + self.maximum_sample_pixels != self.pixels
            || self.minimum_samples < policy.minimum_samples_per_pixel
            || self.maximum_samples > policy.maximum_samples_per_pixel
        {
            return Err(CinematicFixtureError::Pipeline(
                "adaptive sample evidence did not cover the exact rendered frame window".into(),
            ));
        }
        Ok(())
    }

    fn sequence_identity(&self) -> ContentHash {
        let mut hasher =
            DomainHasher::new("org.frankensim.euler-critique.adaptive-sample-count-sequence.v1");
        for frame in &self.frames {
            hasher.update(&frame.frame.to_le_bytes());
            hasher.update(frame.sample_count_identity.as_bytes());
        }
        hasher.finalize()
    }

    fn manifest_json(&self, policy: CinematicAdaptiveSamplingConfig) -> String {
        use core::fmt::Write as _;

        let policy = policy
            .canonical()
            .expect("validated adaptive evidence retains a canonical policy");

        let mut frames = String::from("[");
        for (index, frame) in self.frames.iter().enumerate() {
            if index != 0 {
                frames.push(',');
            }
            let _ = write!(
                frames,
                concat!(
                    "{{\"frame\":{},\"pixels\":{},\"minimum_spp\":{},",
                    "\"maximum_spp\":{},\"total_paths\":{},",
                    "\"error_threshold_pixels\":{},\"maximum_spp_pixels\":{},",
                    "\"sample_count_identity\":\"{}\"}}"
                ),
                frame.frame,
                frame.pixels,
                frame.minimum_samples,
                frame.maximum_samples,
                frame.total_samples,
                frame.converged_pixels,
                frame.maximum_sample_pixels,
                frame.sample_count_identity.to_hex(),
            );
        }
        frames.push(']');
        format!(
            concat!(
                "{{\"mode\":\"adaptive-raw-xyz-dispersion-v1\",",
                "\"minimum_spp\":{},\"maximum_spp\":{},\"decision_batch_spp\":{},",
                "\"absolute_error_xyz\":{:.17e},\"relative_error_xyz\":{:.17e},",
                "\"dark_floor_xyz\":{:.17e},\"sample_count_channel\":\"samples\",",
                "\"rendered_frames\":{},\"rendered_pixels\":{},",
                "\"actual_minimum_spp\":{},\"actual_maximum_spp\":{},",
                "\"total_paths\":{},\"error_threshold_pixels\":{},",
                "\"maximum_spp_pixels\":{},\"sample_count_sequence_identity\":\"{}\",",
                "\"frames\":{},",
                "\"claim\":\"raw per-channel XYZ dispersion stopping; not a confidence interval, perceptual error bound, denoiser decision, or physical-validation claim\"}}"
            ),
            policy.minimum_samples_per_pixel,
            policy.maximum_samples_per_pixel,
            policy.decision_batch_samples,
            policy.absolute_error,
            policy.relative_error,
            policy.dark_floor,
            self.frames.len(),
            self.pixels,
            self.minimum_samples,
            self.maximum_samples,
            self.total_samples,
            self.converged_pixels,
            self.maximum_sample_pixels,
            self.sequence_identity().to_hex(),
            frames,
        )
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

#[derive(Clone, Copy, Debug, PartialEq)]
struct FixtureOutputConvergenceEvidence {
    shutter_pose_sample_count: usize,
    shutter_exact_renderer_sample_count: usize,
    shutter_stratum_boundary_sample_count: usize,
    shutter_interpolation_knot_sample_count: usize,
    shutter_interpolation_midpoint_sample_count: usize,
    terminal_time_difference_s: f64,
    terminal_time_relative_difference: f64,
    coarse_preroll_normal_impulse_n_s: f64,
    fine_preroll_normal_impulse_n_s: f64,
    preroll_normal_impulse_relative_difference: f64,
    coarse_source_normal_impulse_n_s: f64,
    fine_source_normal_impulse_n_s: f64,
    source_normal_impulse_relative_difference: f64,
    coarse_published_normal_impulse_n_s: f64,
    fine_published_normal_impulse_n_s: f64,
    published_normal_impulse_relative_difference: f64,
    maximum_cumulative_normal_impulse_difference_n_s: f64,
    maximum_cumulative_normal_impulse_relative_difference: f64,
    maximum_impulse_identity_residual_n_s: f64,
    maximum_impulse_identity_tolerance_n_s: f64,
    maximum_center_of_mass_difference_m: f64,
    maximum_orientation_difference_rad: f64,
    maximum_chirp_difference_hz: f64,
    maximum_relative_chirp_difference: f64,
}

impl FixtureOutputConvergenceEvidence {
    fn diagnostics(self) -> String {
        format!(
            concat!(
                "stage=mechanics output_convergence source=fine-dt2 ",
                "shutter_queries={} exact_renderer_samples={} stratum_boundaries={} ",
                "interpolation_knots={} interpolation_midpoints={} ",
                "terminal_time_difference_s={:.6e} terminal_time_difference_rel={:.6e} ",
                "preroll_impulse_rel={:.6e} source_impulse_rel={:.6e} ",
                "published_impulse_rel={:.6e} cumulative_impulse_rel={:.6e} ",
                "impulse_identity_residual_n_s={:.6e} impulse_identity_tolerance_n_s={:.6e} ",
                "com_difference_m={:.6e} orientation_difference_rad={:.6e} ",
                "chirp_difference_hz={:.6e} chirp_difference_rel={:.6e}"
            ),
            self.shutter_pose_sample_count,
            self.shutter_exact_renderer_sample_count,
            self.shutter_stratum_boundary_sample_count,
            self.shutter_interpolation_knot_sample_count,
            self.shutter_interpolation_midpoint_sample_count,
            self.terminal_time_difference_s,
            self.terminal_time_relative_difference,
            self.preroll_normal_impulse_relative_difference,
            self.source_normal_impulse_relative_difference,
            self.published_normal_impulse_relative_difference,
            self.maximum_cumulative_normal_impulse_relative_difference,
            self.maximum_impulse_identity_residual_n_s,
            self.maximum_impulse_identity_tolerance_n_s,
            self.maximum_center_of_mass_difference_m,
            self.maximum_orientation_difference_rad,
            self.maximum_chirp_difference_hz,
            self.maximum_relative_chirp_difference,
        )
    }

    fn manifest_json(self) -> String {
        format!(
            concat!(
                "{{\"comparison\":\"coarse-dt versus fine-dt2; fine-dt2 is the published picture-and-sound source\",",
                "\"shutter_pose_query_count\":{},",
                "\"shutter_exact_renderer_sample_count\":{},",
                "\"shutter_stratum_boundary_sample_count\":{},",
                "\"shutter_interpolation_knot_sample_count\":{},",
                "\"shutter_interpolation_midpoint_sample_count\":{},",
                "\"shutter_sampling\":\"one exact pixel-zero renderer sample per temporal stratum and all stratum boundaries, plus exhaustive coarse/fine union-knot and adjacent midpoint coverage inside every exposure; not exact all-ray jitter coverage\",",
                "\"terminal_time\":{{\"absolute_difference_s\":{:.17e},",
                "\"absolute_limit_s\":{:.17e},\"relative_difference\":{:.17e},",
                "\"relative_limit\":{:.17e}}},",
                "\"normal_impulse\":{{\"coarse_preroll_n_s\":{:.17e},\"fine_preroll_n_s\":{:.17e},",
                "\"preroll_relative_difference\":{:.17e},\"coarse_source_n_s\":{:.17e},",
                "\"fine_source_n_s\":{:.17e},\"source_relative_difference\":{:.17e},",
                "\"coarse_published_n_s\":{:.17e},",
                "\"fine_published_n_s\":{:.17e},\"published_relative_difference\":{:.17e},",
                "\"maximum_cumulative_difference_n_s\":{:.17e},",
                "\"maximum_cumulative_relative_difference\":{:.17e},",
                "\"maximum_identity_residual_n_s\":{:.17e},",
                "\"maximum_identity_tolerance_n_s\":{:.17e},",
                "\"relative_difference_limit\":{:.17e}}},",
                "\"motion\":{{\"maximum_center_of_mass_difference_m\":{:.17e},",
                "\"center_of_mass_limit_m\":{:.17e},",
                "\"maximum_orientation_geodesic_difference_rad\":{:.17e},",
                "\"orientation_limit_rad\":{:.17e}}},",
                "\"chirp\":{{\"maximum_difference_hz\":{:.17e},",
                "\"absolute_limit_hz\":{:.17e},\"maximum_relative_difference\":{:.17e},",
                "\"relative_limit\":{:.17e}}},",
                "\"claim\":\"one output-consistency pair for this encoded reduced model; not an asymptotic-order certificate, contact-law validation, acoustic calibration, or experimental validation\"}}"
            ),
            self.shutter_pose_sample_count,
            self.shutter_exact_renderer_sample_count,
            self.shutter_stratum_boundary_sample_count,
            self.shutter_interpolation_knot_sample_count,
            self.shutter_interpolation_midpoint_sample_count,
            self.terminal_time_difference_s,
            OUTPUT_CONVERGENCE_TERMINAL_TIME_ABSOLUTE_LIMIT_S,
            self.terminal_time_relative_difference,
            OUTPUT_CONVERGENCE_TERMINAL_TIME_RELATIVE_LIMIT,
            self.coarse_preroll_normal_impulse_n_s,
            self.fine_preroll_normal_impulse_n_s,
            self.preroll_normal_impulse_relative_difference,
            self.coarse_source_normal_impulse_n_s,
            self.fine_source_normal_impulse_n_s,
            self.source_normal_impulse_relative_difference,
            self.coarse_published_normal_impulse_n_s,
            self.fine_published_normal_impulse_n_s,
            self.published_normal_impulse_relative_difference,
            self.maximum_cumulative_normal_impulse_difference_n_s,
            self.maximum_cumulative_normal_impulse_relative_difference,
            self.maximum_impulse_identity_residual_n_s,
            self.maximum_impulse_identity_tolerance_n_s,
            OUTPUT_CONVERGENCE_RELATIVE_IMPULSE_LIMIT,
            self.maximum_center_of_mass_difference_m,
            OUTPUT_CONVERGENCE_COM_LIMIT_M,
            self.maximum_orientation_difference_rad,
            OUTPUT_CONVERGENCE_ORIENTATION_LIMIT_RAD,
            self.maximum_chirp_difference_hz,
            OUTPUT_CONVERGENCE_CHIRP_ABSOLUTE_LIMIT_HZ,
            self.maximum_relative_chirp_difference,
            OUTPUT_CONVERGENCE_CHIRP_RELATIVE_LIMIT,
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
struct FixtureImpulseAudit {
    cumulative_at_queries_n_s: Vec<f64>,
    total_n_s: f64,
    identity_residual_n_s: f64,
    identity_tolerance_n_s: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.correction;
        let next = self.sum + corrected;
        self.correction = (next - self.sum) - corrected;
        self.sum = next;
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct SymmetricErrorAccumulator {
    sample_count: u64,
    coarse_squared: CompensatedSum,
    fine_squared: CompensatedSum,
    difference_squared: CompensatedSum,
    coarse_peak: f64,
    fine_peak: f64,
    difference_peak: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct RelativeErrorMetrics {
    nrmse: f64,
    normalized_peak: f64,
}

impl SymmetricErrorAccumulator {
    fn observe(&mut self, coarse: f64, fine: f64) -> Result<(), CinematicFixtureError> {
        if !(coarse.is_finite() && fine.is_finite()) {
            return Err(CinematicFixtureError::Pipeline(
                "audio convergence encountered a non-finite sample".into(),
            ));
        }
        let difference = coarse - fine;
        let coarse_squared = coarse * coarse;
        let fine_squared = fine * fine;
        let difference_squared = difference * difference;
        if !(difference.is_finite()
            && coarse_squared.is_finite()
            && fine_squared.is_finite()
            && difference_squared.is_finite())
        {
            return Err(CinematicFixtureError::Pipeline(
                "audio convergence metric arithmetic overflowed".into(),
            ));
        }
        self.sample_count = self.sample_count.checked_add(1).ok_or_else(|| {
            CinematicFixtureError::Pipeline("audio convergence sample count overflow".into())
        })?;
        self.coarse_squared.add(coarse_squared);
        self.fine_squared.add(fine_squared);
        self.difference_squared.add(difference_squared);
        self.coarse_peak = self.coarse_peak.max(coarse.abs());
        self.fine_peak = self.fine_peak.max(fine.abs());
        self.difference_peak = self.difference_peak.max(difference.abs());
        Ok(())
    }

    fn metrics(self, floor: f64) -> Result<RelativeErrorMetrics, CinematicFixtureError> {
        if self.sample_count == 0 {
            return Err(CinematicFixtureError::Pipeline(
                "audio convergence metric has no samples".into(),
            ));
        }
        if !(floor.is_finite() && floor > 0.0) {
            return Err(CinematicFixtureError::Pipeline(
                "audio convergence normalization floor must be finite and positive".into(),
            ));
        }
        let count = self.sample_count as f64;
        let coarse_rms = (self.coarse_squared.sum / count).sqrt();
        let fine_rms = (self.fine_squared.sum / count).sqrt();
        let difference_rms = (self.difference_squared.sum / count).sqrt();
        let rms_scale = coarse_rms.max(fine_rms).max(floor);
        let peak_scale = self.coarse_peak.max(self.fine_peak).max(floor);
        let metrics = RelativeErrorMetrics {
            nrmse: difference_rms / rms_scale,
            normalized_peak: self.difference_peak / peak_scale,
        };
        if !(metrics.nrmse.is_finite()
            && metrics.nrmse >= 0.0
            && metrics.normalized_peak.is_finite()
            && metrics.normalized_peak >= 0.0)
        {
            return Err(CinematicFixtureError::Pipeline(
                "audio convergence produced an invalid normalized metric".into(),
            ));
        }
        Ok(metrics)
    }
}

struct FixtureAudioCandidate {
    resampler: AudioResampler,
    source_sound: SoundSynthesisConfig,
    excitation_identity: ContentHash,
    source_trajectory_identity: ContentHash,
}

struct FixtureAudioPairOutput {
    fine_stems: Vec<ModalStemFrame>,
    fine_sound: SoundSynthesisConfig,
    fine_crop: AudioResamplingCrop,
    convergence: FixtureAudioConvergenceEvidence,
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

fn cinematic_thorne_2026_refinement_evidence(
    benchmark: &Thorne2026SteelGlassBenchmark,
) -> Result<RefinementEvidence, ReducedDecayError> {
    let mut refined = benchmark.clone();
    refined.decay.timestep_s /= f64::from(CINEMATIC_MECHANICS_COARSE_REFINEMENT_FACTOR);
    refined.decay.maximum_steps = refined
        .decay
        .maximum_steps
        .checked_mul(CINEMATIC_MECHANICS_COARSE_REFINEMENT_FACTOR)
        .ok_or(ReducedDecayError::RefinementStepBudgetOverflow)?;
    thorne_2026_refinement_evidence(&refined)
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
    let refinement = cinematic_thorne_2026_refinement_evidence(&benchmark).map_err(pipeline)?;
    let audio_preroll_duration_s = f64::from(AUDIO_PREROLL_VIDEO_FRAMES) / f64::from(CRITIQUE_FPS);
    let (coarse_audio_preroll_trajectory, coarse_trajectory) =
        RenderTrajectory::from_reduced_decay_run_with_causal_preroll(
            &refinement.coarse,
            &profile,
            audio_preroll_duration_s,
            duration_s,
            cx,
        )
        .map_err(pipeline)?;
    let (audio_preroll_trajectory, trajectory) =
        RenderTrajectory::from_reduced_decay_run_with_causal_preroll(
            &refinement.fine,
            &profile,
            audio_preroll_duration_s,
            duration_s,
            cx,
        )
        .map_err(pipeline)?;
    let output_convergence = fixture_output_convergence_evidence(
        &coarse_audio_preroll_trajectory,
        &coarse_trajectory,
        &audio_preroll_trajectory,
        &trajectory,
        config,
        &refinement,
        audio_preroll_duration_s,
        duration_s,
        cx,
    )?;
    progress(&output_convergence.diagnostics());
    // Once the fixed mechanics gate passes, retain the coarse full-preroll
    // source only through the pre-master audio gate. The coarse published
    // picture trajectory has served its sole comparison purpose.
    drop(coarse_trajectory);
    let run = &refinement.fine;
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
    let audio_preroll_artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
        hash_domain(
            "org.frankensim.euler-critique.audio-preroll-source.v1",
            trajectory_artifact.receipt().artifact_identity().as_bytes(),
        ),
        audio_preroll_trajectory,
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .map_err(pipeline)?;
    let coarse_audio_preroll_artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
        hash_domain(
            "org.frankensim.euler-critique.coarse-audio-convergence-source.v1",
            coarse_audio_preroll_trajectory
                .metadata()
                .configuration_identity
                .as_bytes(),
        ),
        coarse_audio_preroll_trajectory,
        Vec::new(),
        RenderTrajectoryCodecBudget::DEFAULT,
        cx,
    )
    .map_err(pipeline)?;
    progress(&trajectory_diagnostics(&trajectory_artifact));

    // Sound construction is much cheaper than path tracing and depends only
    // on the admitted trajectory. Fail here before spending render-hours if a
    // synthesis, mastering, or artifact contract is unsatisfied.
    progress("stage=audio begin");
    let audio = build_audio(
        &trajectory_artifact,
        &coarse_audio_preroll_artifact,
        &audio_preroll_artifact,
        run.parameters.gravity_m_per_s2,
        config,
        cx,
    )?;
    progress(&audio.convergence.diagnostics());
    drop(coarse_audio_preroll_artifact);
    audio
        .artifact
        .verify(AudioArtifactBudget::DEFAULT, cx)
        .map_err(pipeline)?;
    progress("stage=audio complete");

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
    let wav_path = sound_directory.join("master.float32.wav");
    write_new(&wav_path, audio.artifact.wav_bytes())?;
    let audio_manifest_path = sound_directory.join("master.manifest.json");
    write_new(
        &audio_manifest_path,
        audio.artifact.manifest().to_manifest_json().as_bytes(),
    )?;

    progress("stage=render begin");
    let camera = critique_camera(duration_s).map_err(pipeline)?;
    let mut scene_config = EulerSceneConfig::reference(camera);
    // The cinematic defaults increase the actual intersected render geometry
    // so glossy highlights follow matching geometric normals instead of
    // relying on an uncorrected shading-normal approximation. Explicit config
    // controls allow profiling lower-cost render-only approximations without
    // altering the specimen, mechanics, or production defaults.
    scene_config.tessellation = EulerTessellationConfig {
        azimuthal_segments: config.azimuthal_segments,
        arc_subdivisions_per_arc: config.arc_subdivisions_per_arc,
    };
    scene_config.show_spin_fiducial = true;
    scene_config.disc_material = EulerMaterialStyle::Conductor {
        optics: ConductorOptics::representative_stainless_steel(),
        surface: ConductorSurface::try_rough(0.12).map_err(pipeline)?,
    };
    // The product plate is polished glass, not a visibly frosted surface. The
    // reference scene's rough-GGX convenience preset creates rare, extremely
    // bright microfacet paths at practical cinematic sample counts. Use the
    // physically defensible ideal-polished limit here; this changes only the
    // visual material model and never the mechanics/contact base.
    scene_config.plate_material = EulerMaterialStyle::Dielectric {
        glass: DielectricGlass::representative_crown(),
        surface: DielectricSurface::SMOOTH,
    };
    scene_config.support_surface = Some(EulerSupportSurfaceSpec {
        // Cover the full two-metre camera range so the finite tabletop edge
        // cannot cut a black chevron through the background of the hero shot.
        // Scaling this box does not add triangles or alter the mechanical base.
        width_m: 4.0,
        depth_m: 4.0,
        thickness_m: 0.02,
        gap_below_housing_m: 0.0,
        material: EulerMaterialStyle::Lambertian {
            linear_rgb: [0.07, 0.055, 0.045],
        },
    });
    // Keep the emitter above the camera frustum while retaining a broad,
    // downward-facing studio source. The reference light otherwise appears as
    // a distracting white bar across the top of this particular composition.
    scene_config.light = EulerRectLightSpec {
        corner_world_m: Point3::new(-0.175, 0.075, 0.22),
        edge_u_world_m: GeomVec3::new(0.24, 0.0, 0.0),
        edge_v_world_m: GeomVec3::new(0.0, -0.18, 0.0),
        linear_rgb: [1.0, 0.96, 0.90],
        radiance_scale: 24.0,
    };
    scene_config.environment = EulerEnvironmentStyle::StudioGradient(EulerStudioEnvironmentSpec {
        overhead_linear_rgb: [0.12, 0.16, 0.24],
        horizon_linear_rgb: [0.22, 0.15, 0.10],
        floor_linear_rgb: [0.012, 0.009, 0.007],
        radiance_scale: 0.35,
    });
    let scene = EulerCinematicScene::try_build(&trajectory_artifact, &profile, scene_config, cx)
        .map_err(pipeline)?;
    let preview_mesh_evidence =
        FixturePreviewMeshEvidence::from_receipt(scene.preview_mesh_receipt());
    let render_sample_ceiling = config.render_sample_ceiling();
    let adaptive_policy = config.adaptive_policy()?;
    let mut base_render_settings = euler_scene_smoke_settings(config.width, config.height);
    base_render_settings.spp = render_sample_ceiling;
    base_render_settings.max_depth = config.max_depth;
    let composition_identity = composition_identity(config, scene.scene_identity());
    let mut raw_sequence = DomainHasher::new("org.frankensim.euler-critique.raw-sequence.v1");
    let mut preview_sequence =
        DomainHasher::new("org.frankensim.euler-critique.preview-sequence.v1");
    let mut over_range_channels = 0_u64;
    let mut gamut_mapped_pixels = 0_u64;
    let mut render_evidence = FixtureRenderEvidence::default();
    let mut adaptive_sampling_evidence = FixtureAdaptiveSamplingEvidence::default();
    let mut denoise_evidence = FixtureDenoiseEvidence::default();
    let mut denoise_history: Option<TemporalDenoisedFrame> = None;
    let denoise_config = TemporalDenoiseConfig::default();
    let aov_profile = if config.denoise_previews || config.retain_full_aov_exr {
        CinematicAovProfile::FinalDiagnostic
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
            let frame_times =
                frame_timeline_times(frame, config.frames, trajectory_start_s, trajectory_end_s);
            let prepared = scene
                .prepare_frame(EulerFrameRequest {
                    frame_time_s: frame_times.shutter_close_time_s,
                    exposure_duration_s,
                    // Treat each image as the exposure ending at its video
                    // interval's mechanics boundary. Encoded presentation time
                    // remains at the interval start, while the final shutter
                    // closes on the analytical validity cutoff at exactly 8 s.
                    convention: ShutterConvention::BackLoaded,
                    distribution: ShutterDistribution::StratifiedCounterV1 {
                        strata: render_sample_ceiling,
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
                frame_times.presentation_time_s,
                frame_times.previous_presentation_time_s,
                frame_times.next_presentation_time_s,
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
            let color = critique_color_config();
            let (exr, preview, report, sampling_progress) = if let Some(policy) = adaptive_policy {
                let adaptive_config = config
                    .adaptive_sampling
                    .expect("validated adaptive policy retains its fixture configuration");
                let output = renderer
                    .render_cinematic_adaptive_with_aovs(
                        scene.scene(),
                        scene.camera(),
                        prepared.cut_side(),
                        cx,
                        &render_settings,
                        policy,
                        prepared.segments()[0].shutter(),
                        CinematicAovConfig::new(
                            aov_profile,
                            provenance,
                            CinematicAovLimits::default(),
                        ),
                        &frame_execution,
                    )
                    .map_err(pipeline)?;
                let film = output.film;
                let summary = film.beauty().summary();
                let sample_count_identity =
                    adaptive_sampling_evidence.observe(frame, film.beauty(), adaptive_config)?;
                let [red, green, blue] = film.beauty().to_linear_srgb();
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
                            samples_per_pixel: render_settings.spp,
                            sample_counts_per_pixel: Some(film.beauty().sample_counts()),
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
                            object_ids: guides.object_palette_indices(),
                            material_ids: guides.material_palette_indices(),
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
                let exr = if config.retain_full_aov_exr {
                    film.to_exr().map_err(pipeline)?
                } else {
                    adaptive_beauty_to_exr(
                        film.beauty(),
                        [red, green, blue],
                        adaptive_config,
                        render_settings,
                        provenance,
                        sample_count_identity,
                    )?
                };
                let sampling_progress = format!(
                    concat!(
                        "sampling=adaptive min_spp={} max_spp={} total_paths={} ",
                        "error_threshold_pixels={} maximum_spp_pixels={}"
                    ),
                    summary.minimum_samples,
                    summary.maximum_samples,
                    summary.total_samples,
                    summary.converged_pixels,
                    summary.maximum_sample_pixels,
                );
                (exr, preview, output.report, sampling_progress)
            } else {
                let output = renderer
                    .render_cinematic_with_aovs(
                        scene.scene(),
                        scene.camera(),
                        prepared.cut_side(),
                        cx,
                        &render_settings,
                        prepared.segments()[0].shutter(),
                        CinematicAovConfig::new(
                            aov_profile,
                            provenance,
                            CinematicAovLimits::default(),
                        ),
                        &frame_execution,
                    )
                    .map_err(pipeline)?;
                let film = output.film;
                let exr = if config.retain_full_aov_exr {
                    film.to_exr().map_err(pipeline)?
                } else {
                    film_to_exr(film.beauty()).map_err(pipeline)?
                };
                let [red, green, blue] = film.beauty().to_linear_srgb();
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
                            samples_per_pixel: config.samples_per_pixel,
                            sample_counts_per_pixel: None,
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
                            object_ids: guides.object_palette_indices(),
                            material_ids: guides.material_palette_indices(),
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
                (
                    exr,
                    preview,
                    output.report,
                    format!("sampling=uniform spp={}", config.samples_per_pixel),
                )
            };
            progress(&format!(
                concat!(
                    "stage=render frame={}/{} workers={} tiles={} ",
                    "traversal_ms={:.3} compute_ms={:.3} merge_ms={:.3} peak_mib={:.3} {}"
                ),
                frame - render_frame_range.start + 1,
                render_frame_range.len(),
                report.workers,
                report.layout.tile_count(),
                report.traversal_ns as f64 / 1.0e6,
                report.tile_compute_ns as f64 / 1.0e6,
                report.tile_merge_ns as f64 / 1.0e6,
                report.memory.peak_bytes as f64 / (1024.0 * 1024.0),
                sampling_progress,
            ));
            render_evidence.observe(&report);
            raw_sequence.update(hash_domain("frame", &exr).as_bytes());
            write_new(&raw_directory.join(format!("frame-{frame:06}.exr")), &exr)?;
            drop(exr);
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
    adaptive_sampling_evidence.validate_complete(config, &render_frame_range)?;
    let raw_sequence_identity = raw_sequence.finalize();
    let preview_sequence_identity = preview_sequence.finalize();
    progress("stage=render complete");

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
        &output_convergence,
        &audio.convergence,
        &trajectory_artifact,
        raw_sequence_identity,
        preview_sequence_identity,
        audio.artifact.manifest().wav().wav_identity(),
        audio.modal_parameter_set_identity,
        &audio.modal_parameter_set_disclosure,
        audio.warm_start_source_identity,
        audio.warm_start_checkpoint_identity,
        audio.source_sound_configuration_identity,
        audio.published_trajectory_identity,
        audio.crop_resampler_identity,
        audio.crop_first_source_audio_frame,
        audio.crop_end_source_audio_frame,
        audio.chirp_start_hz,
        audio.chirp_end_hz,
        audio.pre_master_peak_fs,
        audio.master_gain_db,
        over_range_channels,
        gamut_mapped_pixels,
        &render_evidence,
        preview_mesh_evidence,
        &adaptive_sampling_evidence,
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

fn adaptive_beauty_to_exr(
    film: &AdaptiveFilm,
    rgb: [Vec<f32>; 3],
    policy: CinematicAdaptiveSamplingConfig,
    settings: Settings,
    provenance: CinematicAovProvenance,
    sample_count_identity: ContentHash,
) -> Result<Vec<u8>, CinematicFixtureError> {
    let policy = policy.canonical()?;
    let pixel_count = film.sample_counts().len();
    if rgb.iter().any(|plane| plane.len() != pixel_count) {
        return Err(CinematicFixtureError::Pipeline(
            "adaptive RGB and sample-count planes have different shapes".into(),
        ));
    }
    let mut samples = Vec::new();
    samples.try_reserve_exact(pixel_count).map_err(|_| {
        CinematicFixtureError::Pipeline(
            "adaptive EXR sample-count channel allocation refused".into(),
        )
    })?;
    for count in film.sample_counts() {
        // The fixture ceiling is 4096, so every count is exactly representable
        // in the float channel used by the frozen fs-img EXR subset.
        #[allow(clippy::cast_precision_loss)]
        samples.push(*count as f32);
    }
    let [red, green, blue] = rgb;
    let channels = [
        Channel {
            name: "R".to_owned(),
            ty: PixelType::Float,
            data: red,
        },
        Channel {
            name: "G".to_owned(),
            ty: PixelType::Float,
            data: green,
        },
        Channel {
            name: "B".to_owned(),
            ty: PixelType::Float,
            data: blue,
        },
        Channel {
            name: "samples".to_owned(),
            ty: PixelType::Float,
            data: samples,
        },
    ];
    let string_attribute = |name: &str, value: String| ExrAttribute {
        name: name.to_owned(),
        ty: "string".to_owned(),
        value: value.into_bytes(),
    };
    let attributes = [
        string_attribute("frankensim.aov.authority", "raw-estimate".to_owned()),
        string_attribute(
            "frankensim.aov.profile",
            "fixture-beauty-plus-samples-v1".to_owned(),
        ),
        string_attribute(
            "frankensim.frame.index",
            provenance.frame_index().to_string(),
        ),
        string_attribute(
            "frankensim.frame.timeSeconds",
            format!("0x{:016x}", provenance.frame_time_s().to_bits()),
        ),
        string_attribute(
            "frankensim.frame.previousTimeS",
            format!("0x{:016x}", provenance.previous_frame_time_s().to_bits()),
        ),
        string_attribute(
            "frankensim.frame.nextTimeS",
            format!("0x{:016x}", provenance.next_frame_time_s().to_bits()),
        ),
        string_attribute(
            "frankensim.source.trajectory",
            provenance.source_trajectory_identity().to_hex(),
        ),
        string_attribute(
            "frankensim.source.sceneHash",
            provenance.scene_identity().to_hex(),
        ),
        string_attribute(
            "frankensim.source.composition",
            provenance.composition_identity().to_hex(),
        ),
        string_attribute("frankensim.render.sampleMode", "adaptive".to_owned()),
        string_attribute("frankensim.render.seed", settings.seed.to_string()),
        string_attribute(
            "frankensim.render.sampler",
            match settings.sampler {
                fs_render::tracer::Sampler::Iid => "iid-philox",
                fs_render::tracer::Sampler::OwenSobol => "owen-sobol",
            }
            .to_owned(),
        ),
        string_attribute(
            "frankensim.render.strategy",
            match settings.strategy {
                fs_render::tracer::DirectStrategy::NeeOnly => "nee-only",
                fs_render::tracer::DirectStrategy::BsdfOnly => "bsdf-only",
                fs_render::tracer::DirectStrategy::Mis => "mis",
            }
            .to_owned(),
        ),
        string_attribute("frankensim.render.maxDepth", settings.max_depth.to_string()),
        string_attribute("frankensim.render.spp", "per-pixel-channel".to_owned()),
        string_attribute(
            "frankensim.render.sppCeiling",
            policy.maximum_samples_per_pixel.to_string(),
        ),
        string_attribute(
            "frankensim.render.adaptive",
            format!(
                concat!(
                    "version={};minimum={};batch={};",
                    "absolute=0x{:016x};relative=0x{:016x};darkFloor=0x{:016x}"
                ),
                ADAPTIVE_SAMPLING_SEMANTICS_VERSION,
                policy.minimum_samples_per_pixel,
                policy.decision_batch_samples,
                policy.absolute_error.to_bits(),
                policy.relative_error.to_bits(),
                policy.dark_floor.to_bits(),
            ),
        ),
        string_attribute(
            "frankensim.render.sampleCounts",
            sample_count_identity.to_hex(),
        ),
        string_attribute(
            "frankensim.render.versions",
            cinematic_render_semantics_versions(),
        ),
    ];
    write_exr_with_attributes(film.width(), film.height(), &channels, &attributes).map_err(pipeline)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FixtureFrameTimelineTimes {
    presentation_time_s: f64,
    previous_presentation_time_s: f64,
    next_presentation_time_s: f64,
    shutter_close_time_s: f64,
}

fn frame_timeline_times(
    frame: u32,
    frames: u32,
    trajectory_start_s: f64,
    trajectory_end_s: f64,
) -> FixtureFrameTimelineTimes {
    debug_assert!(frames > 0 && frame < frames);
    let fps = f64::from(CRITIQUE_FPS);
    let presentation_time_s = f64::from(frame) / fps;
    let shutter_close_time_s = (f64::from(frame) + 1.0) / fps;
    // Motion-vector references use encoded presentation timestamps. The next
    // presentation boundary also encloses this frame's back-loaded shutter.
    let previous_presentation_time_s = if frame == 0 {
        trajectory_start_s
    } else {
        (f64::from(frame) - 1.0) / fps
    }
    .max(trajectory_start_s);
    let next_presentation_time_s = shutter_close_time_s.min(trajectory_end_s);
    debug_assert!(trajectory_start_s <= presentation_time_s);
    debug_assert!(shutter_close_time_s <= trajectory_end_s);
    FixtureFrameTimelineTimes {
        presentation_time_s,
        previous_presentation_time_s,
        next_presentation_time_s,
        shutter_close_time_s,
    }
}

fn composition_identity(config: &CinematicFixtureConfig, scene: ContentHash) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.composition.v4");
    hasher.update(scene.as_bytes());
    hasher.update(&config.width.to_le_bytes());
    hasher.update(&config.height.to_le_bytes());
    hasher.update(&config.frames.to_le_bytes());
    hasher.update(&CRITIQUE_FPS.to_le_bytes());
    hasher.update(&config.render_sample_ceiling().to_le_bytes());
    if let Some(policy) = config.adaptive_sampling {
        let admitted = policy
            .policy()
            .expect("validated fixture configuration retains an admitted adaptive policy");
        hasher.update(b"adaptive-raw-xyz-dispersion-v1");
        hasher.update(&admitted.minimum_samples().to_le_bytes());
        hasher.update(&admitted.batch_samples().to_le_bytes());
        hasher.update(&admitted.absolute_error().to_bits().to_le_bytes());
        hasher.update(&admitted.relative_error().to_bits().to_le_bytes());
        hasher.update(&admitted.dark_floor().to_bits().to_le_bytes());
    }
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
        b"encoded-pts-frame-start-v1;back-loaded-shutter-close-at-frame-end-v2;cinematic-aov;display-color-config-canonical-v1",
    );
    hasher.update(
        &critique_color_config()
            .canonical_bytes()
            .expect("the frozen critique display configuration is valid"),
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

fn fixture_output_convergence_evidence(
    coarse_source: &RenderTrajectory,
    coarse_published: &RenderTrajectory,
    fine_source: &RenderTrajectory,
    fine_published: &RenderTrajectory,
    config: &CinematicFixtureConfig,
    refinement: &RefinementEvidence,
    preroll_duration_s: f64,
    published_duration_s: f64,
    cx: &Cx<'_>,
) -> Result<FixtureOutputConvergenceEvidence, CinematicFixtureError> {
    let gravity_m_per_s2 = refinement.fine.parameters.gravity_m_per_s2;
    let coarse_dt_s = coarse_published.metadata().timestep_s;
    let fine_dt_s = fine_published.metadata().timestep_s;
    let timestep_scale = coarse_dt_s.abs().max(f64::MIN_POSITIVE);
    if coarse_source.metadata().timestep_s.to_bits() != coarse_dt_s.to_bits()
        || fine_source.metadata().timestep_s.to_bits() != fine_dt_s.to_bits()
        || refinement.coarse.parameters.timestep_s.to_bits() != coarse_dt_s.to_bits()
        || refinement.fine.parameters.timestep_s.to_bits() != fine_dt_s.to_bits()
        || (coarse_dt_s - 2.0 * fine_dt_s).abs() > 8.0 * f64::EPSILON * timestep_scale
    {
        return Err(CinematicFixtureError::Pipeline(format!(
            "output convergence requires one exact dt/dt2 source pair; coarse={coarse_dt_s:.17e}s fine={fine_dt_s:.17e}s"
        )));
    }
    if coarse_published.metadata().mass_properties.properties
        != fine_published.metadata().mass_properties.properties
        || coarse_source.metadata().mass_properties.properties
            != fine_source.metadata().mass_properties.properties
        || coarse_source.metadata().mass_properties.properties
            != coarse_published.metadata().mass_properties.properties
    {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence trajectories do not share exact mass properties".into(),
        ));
    }
    if !(gravity_m_per_s2.is_finite()
        && gravity_m_per_s2 > 0.0
        && refinement.coarse.parameters.gravity_m_per_s2.to_bits() == gravity_m_per_s2.to_bits()
        && preroll_duration_s.is_finite()
        && preroll_duration_s > 0.0
        && published_duration_s.is_finite()
        && published_duration_s > 0.0)
    {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence horizons and gravity must be finite and positive".into(),
        ));
    }
    let coarse_terminal_time_s = refinement
        .coarse
        .samples
        .last()
        .ok_or_else(|| {
            CinematicFixtureError::Pipeline(
                "output convergence coarse run has no terminal sample".into(),
            )
        })?
        .time_s;
    let fine_terminal_time_s = refinement
        .fine
        .samples
        .last()
        .ok_or_else(|| {
            CinematicFixtureError::Pipeline(
                "output convergence fine run has no terminal sample".into(),
            )
        })?
        .time_s;
    if !(coarse_terminal_time_s.is_finite()
        && coarse_terminal_time_s > 0.0
        && fine_terminal_time_s.is_finite()
        && fine_terminal_time_s > 0.0)
    {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence terminal times must be finite and positive".into(),
        ));
    }
    let terminal_time_difference_s = (coarse_terminal_time_s - fine_terminal_time_s).abs();
    if terminal_time_difference_s.to_bits() != refinement.terminal_time_difference_s.to_bits() {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence terminal-time evidence does not match its retained endpoints"
                .into(),
        ));
    }
    let terminal_time_relative_difference = terminal_time_difference_s
        / coarse_terminal_time_s
            .abs()
            .max(fine_terminal_time_s.abs())
            .max(f64::MIN_POSITIVE);
    let source_duration_s = preroll_duration_s + published_duration_s;
    for (label, trajectory, expected_duration_s) in [
        ("coarse source", coarse_source, source_duration_s),
        ("coarse published", coarse_published, published_duration_s),
        ("fine source", fine_source, source_duration_s),
        ("fine published", fine_published, published_duration_s),
    ] {
        require_trajectory_horizon(label, trajectory, expected_duration_s)?;
    }
    if coarse_source
        .samples()
        .last()
        .expect("validated trajectory is nonempty")
        .input()
        .disposition
        != fine_source
            .samples()
            .last()
            .expect("validated trajectory is nonempty")
            .input()
            .disposition
        || coarse_published
            .samples()
            .last()
            .expect("validated trajectory is nonempty")
            .input()
            .disposition
            != fine_published
                .samples()
                .last()
                .expect("validated trajectory is nonempty")
                .input()
                .disposition
    {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence trajectories end with different disposition classes".into(),
        ));
    }

    let frame_boundary_times_s = (0..=config.frames)
        .map(|frame| f64::from(frame) / f64::from(CRITIQUE_FPS))
        .collect::<Vec<_>>();
    let coarse_source_impulse =
        fixture_impulse_audit(coarse_source, &[preroll_duration_s], gravity_m_per_s2, cx)?;
    let fine_source_impulse =
        fixture_impulse_audit(fine_source, &[preroll_duration_s], gravity_m_per_s2, cx)?;
    let coarse_published_impulse = fixture_impulse_audit(
        coarse_published,
        &frame_boundary_times_s,
        gravity_m_per_s2,
        cx,
    )?;
    let fine_published_impulse = fixture_impulse_audit(
        fine_published,
        &frame_boundary_times_s,
        gravity_m_per_s2,
        cx,
    )?;
    for (label, audit) in [
        ("coarse source", &coarse_source_impulse),
        ("fine source", &fine_source_impulse),
        ("coarse published", &coarse_published_impulse),
        ("fine published", &fine_published_impulse),
    ] {
        if audit.identity_residual_n_s > audit.identity_tolerance_n_s {
            return Err(CinematicFixtureError::Pipeline(format!(
                "{label} normal-impulse identity residual {:.17e} N s exceeds its {:.17e} N s roundoff allowance",
                audit.identity_residual_n_s, audit.identity_tolerance_n_s
            )));
        }
    }

    let mass_kg = fine_published.metadata().mass_properties.properties.mass();
    let preroll_impulse_scale_n_s = (mass_kg * gravity_m_per_s2 * preroll_duration_s).abs();
    let source_impulse_scale_n_s = (mass_kg * gravity_m_per_s2 * source_duration_s).abs();
    let published_impulse_scale_n_s = (mass_kg * gravity_m_per_s2 * published_duration_s).abs();
    let coarse_preroll_normal_impulse_n_s = coarse_source_impulse.cumulative_at_queries_n_s[0];
    let fine_preroll_normal_impulse_n_s = fine_source_impulse.cumulative_at_queries_n_s[0];
    let preroll_normal_impulse_relative_difference = scaled_absolute_difference(
        coarse_preroll_normal_impulse_n_s,
        fine_preroll_normal_impulse_n_s,
        preroll_impulse_scale_n_s,
    );
    let source_normal_impulse_relative_difference = scaled_absolute_difference(
        coarse_source_impulse.total_n_s,
        fine_source_impulse.total_n_s,
        source_impulse_scale_n_s,
    );
    let published_normal_impulse_relative_difference = scaled_absolute_difference(
        coarse_published_impulse.total_n_s,
        fine_published_impulse.total_n_s,
        published_impulse_scale_n_s,
    );
    let maximum_cumulative_normal_impulse_difference_n_s = coarse_published_impulse
        .cumulative_at_queries_n_s
        .iter()
        .zip(&fine_published_impulse.cumulative_at_queries_n_s)
        .map(|(coarse, fine)| (coarse - fine).abs())
        .fold(0.0_f64, f64::max);
    let maximum_cumulative_normal_impulse_relative_difference =
        maximum_cumulative_normal_impulse_difference_n_s
            / fine_published_impulse
                .total_n_s
                .abs()
                .max(published_impulse_scale_n_s)
                .max(f64::MIN_POSITIVE);

    let coarse_motion = EulerRenderMotionBridge::new(coarse_published);
    let fine_motion = EulerRenderMotionBridge::new(fine_published);
    let coarse_last_time_s = coarse_published
        .samples()
        .last()
        .expect("validated trajectory is nonempty")
        .input()
        .time_s;
    let fine_last_time_s = fine_published
        .samples()
        .last()
        .expect("validated trajectory is nonempty")
        .input()
        .time_s;
    let pose_query_horizon_tolerance_s = 32.0 * f64::EPSILON * published_duration_s.max(1.0);
    let mass = fine_published.metadata().mass_properties.properties;
    let exposure_duration_s =
        f64::from(config.shutter_angle_degrees) / 360.0 / f64::from(CRITIQUE_FPS);
    let shot = ShotTimeBounds::try_new(0.0, published_duration_s).map_err(pipeline)?;
    let base_render_seed = euler_scene_smoke_settings(config.width, config.height).seed;
    let interpolation_knot_capacity = coarse_published
        .samples()
        .len()
        .checked_add(fine_published.samples().len())
        .ok_or_else(|| {
            CinematicFixtureError::Pipeline("shutter interpolation knot count overflow".into())
        })?;
    let mut interpolation_knots_s = Vec::new();
    interpolation_knots_s
        .try_reserve_exact(interpolation_knot_capacity)
        .map_err(|_| {
            CinematicFixtureError::Pipeline("shutter interpolation knot allocation refused".into())
        })?;
    interpolation_knots_s.extend(
        coarse_published
            .samples()
            .iter()
            .map(|sample| sample.input().time_s),
    );
    interpolation_knots_s.extend(
        fine_published
            .samples()
            .iter()
            .map(|sample| sample.input().time_s),
    );
    interpolation_knots_s.sort_by(f64::total_cmp);
    interpolation_knots_s.dedup_by(|left, right| left.to_bits() == right.to_bits());
    let mut maximum_center_of_mass_difference_m = 0.0_f64;
    let mut maximum_orientation_difference_rad = 0.0_f64;
    let mut maximum_chirp_difference_hz = 0.0_f64;
    let mut maximum_relative_chirp_difference = 0.0_f64;
    let mut shutter_pose_sample_count = 0_usize;
    let mut shutter_exact_renderer_sample_count = 0_usize;
    let mut shutter_stratum_boundary_sample_count = 0_usize;
    let mut shutter_interpolation_knot_sample_count = 0_usize;
    let mut shutter_interpolation_midpoint_sample_count = 0_usize;
    let mut compare_at_time = |time_s: f64| -> Result<(), CinematicFixtureError> {
        let coarse_time_s = canonicalize_terminal_query_time_s(
            time_s,
            coarse_last_time_s,
            pose_query_horizon_tolerance_s,
        );
        let fine_time_s = canonicalize_terminal_query_time_s(
            time_s,
            fine_last_time_s,
            pose_query_horizon_tolerance_s,
        );
        let coarse = coarse_motion
            .sample_at_time(coarse_time_s, EventEvaluationSide::RightLimit)
            .map_err(pipeline)?;
        let fine = fine_motion
            .sample_at_time(fine_time_s, EventEvaluationSide::RightLimit)
            .map_err(pipeline)?;
        let coarse_pose = coarse.timeline_sample().state.pose();
        let fine_pose = fine.timeline_sample().state.pose();
        let position_delta = coarse_pose.position_world().sub(fine_pose.position_world());
        let center_of_mass_difference_m = position_delta.norm_squared().sqrt();
        let orientation_difference_rad = quaternion_geodesic_difference_rad(
            coarse_pose.orientation().components(),
            fine_pose.orientation().components(),
        );
        if !(center_of_mass_difference_m.is_finite()
            && center_of_mass_difference_m >= 0.0
            && orientation_difference_rad.is_finite()
            && orientation_difference_rad >= 0.0)
        {
            return Err(CinematicFixtureError::Pipeline(
                "output convergence encountered a non-finite or negative pose difference".into(),
            ));
        }
        maximum_center_of_mass_difference_m =
            maximum_center_of_mass_difference_m.max(center_of_mass_difference_m);
        maximum_orientation_difference_rad =
            maximum_orientation_difference_rad.max(orientation_difference_rad);
        let coarse_chirp_hz = body_contact_chirp_hz(coarse.timeline_sample().state, mass)?;
        let fine_chirp_hz = body_contact_chirp_hz(fine.timeline_sample().state, mass)?;
        let chirp_difference_hz = (coarse_chirp_hz - fine_chirp_hz).abs();
        maximum_chirp_difference_hz = maximum_chirp_difference_hz.max(chirp_difference_hz);
        maximum_relative_chirp_difference = maximum_relative_chirp_difference
            .max(chirp_difference_hz / fine_chirp_hz.abs().max(f64::MIN_POSITIVE));
        shutter_pose_sample_count = shutter_pose_sample_count.checked_add(1).ok_or_else(|| {
            CinematicFixtureError::Pipeline("shutter convergence query count overflow".into())
        })?;
        Ok(())
    };
    let render_sample_ceiling = config.render_sample_ceiling();
    for frame in 0..config.frames {
        cx.checkpoint()
            .map_err(|_| CinematicFixtureError::Cancelled)?;
        let close_s = f64::from(frame + 1) / f64::from(CRITIQUE_FPS);
        let shutter = ShutterInterval::resolve(
            close_s,
            exposure_duration_s,
            ShutterConvention::BackLoaded,
            ShutterDistribution::StratifiedCounterV1 {
                strata: render_sample_ceiling,
            },
            shot,
        )
        .map_err(pipeline)?;
        let frame_seed = frame_render_seed(base_render_seed, config.render_seed_salt, frame);
        let strata = usize::try_from(render_sample_ceiling).map_err(|_| {
            CinematicFixtureError::Pipeline("shutter stratum count exceeds usize".into())
        })?;
        let mut visited_strata = vec![false; strata];
        // For any fixed pixel, consecutive sample identities 0..SPP form a
        // permutation of every temporal stratum. Pixel zero is therefore a
        // bounded exact subset of renderer timestamps that covers the complete
        // partition; the explicit boundaries below close each stratum.
        for sample in 0..render_sample_ceiling {
            let normalized = shutter.sample_for_stream(frame_seed, 0, u64::from(sample));
            let stratum = (normalized.value() * f64::from(render_sample_ceiling)).floor();
            if !(stratum.is_finite() && stratum >= 0.0 && stratum < strata as f64) {
                return Err(CinematicFixtureError::Pipeline(format!(
                    "frame {frame} renderer shutter sample escaped its temporal strata"
                )));
            }
            let stratum = stratum as usize;
            if core::mem::replace(&mut visited_strata[stratum], true) {
                return Err(CinematicFixtureError::Pipeline(format!(
                    "frame {frame} renderer shutter samples did not form a stratum permutation"
                )));
            }
            compare_at_time(shutter.time_at(normalized))?;
            shutter_exact_renderer_sample_count = shutter_exact_renderer_sample_count
                .checked_add(1)
                .ok_or_else(|| {
                    CinematicFixtureError::Pipeline(
                        "exact renderer shutter sample count overflow".into(),
                    )
                })?;
        }
        if visited_strata.iter().any(|visited| !visited) {
            return Err(CinematicFixtureError::Pipeline(format!(
                "frame {frame} renderer shutter samples left a temporal stratum uncovered"
            )));
        }
        for boundary in 0..=render_sample_ceiling {
            let normalized = NormalizedShutterTime::try_new(
                f64::from(boundary) / f64::from(render_sample_ceiling),
            )
            .map_err(pipeline)?;
            compare_at_time(shutter.time_at(normalized))?;
            shutter_stratum_boundary_sample_count = shutter_stratum_boundary_sample_count
                .checked_add(1)
                .ok_or_else(|| {
                    CinematicFixtureError::Pipeline(
                        "shutter stratum boundary count overflow".into(),
                    )
                })?;
        }
        let first_internal_knot =
            interpolation_knots_s.partition_point(|time_s| *time_s <= shutter.open_s());
        let end_internal_knot =
            interpolation_knots_s.partition_point(|time_s| *time_s < shutter.close_s());
        let internal_knots = &interpolation_knots_s[first_internal_knot..end_internal_knot];
        for time_s in internal_knots.iter().copied() {
            compare_at_time(time_s)?;
            shutter_interpolation_knot_sample_count = shutter_interpolation_knot_sample_count
                .checked_add(1)
                .ok_or_else(|| {
                    CinematicFixtureError::Pipeline(
                        "shutter interpolation knot sample count overflow".into(),
                    )
                })?;
        }
        let mut left_time_s = shutter.open_s();
        for right_time_s in internal_knots
            .iter()
            .copied()
            .chain(core::iter::once(shutter.close_s()))
        {
            if right_time_s > left_time_s {
                let midpoint_s = 0.5_f64.mul_add(right_time_s - left_time_s, left_time_s);
                compare_at_time(midpoint_s)?;
                shutter_interpolation_midpoint_sample_count =
                    shutter_interpolation_midpoint_sample_count
                        .checked_add(1)
                        .ok_or_else(|| {
                            CinematicFixtureError::Pipeline(
                                "shutter interpolation midpoint sample count overflow".into(),
                            )
                        })?;
            }
            left_time_s = right_time_s;
        }
    }

    let evidence = FixtureOutputConvergenceEvidence {
        shutter_pose_sample_count,
        shutter_exact_renderer_sample_count,
        shutter_stratum_boundary_sample_count,
        shutter_interpolation_knot_sample_count,
        shutter_interpolation_midpoint_sample_count,
        terminal_time_difference_s,
        terminal_time_relative_difference,
        coarse_preroll_normal_impulse_n_s,
        fine_preroll_normal_impulse_n_s,
        preroll_normal_impulse_relative_difference,
        coarse_source_normal_impulse_n_s: coarse_source_impulse.total_n_s,
        fine_source_normal_impulse_n_s: fine_source_impulse.total_n_s,
        source_normal_impulse_relative_difference,
        coarse_published_normal_impulse_n_s: coarse_published_impulse.total_n_s,
        fine_published_normal_impulse_n_s: fine_published_impulse.total_n_s,
        published_normal_impulse_relative_difference,
        maximum_cumulative_normal_impulse_difference_n_s,
        maximum_cumulative_normal_impulse_relative_difference,
        maximum_impulse_identity_residual_n_s: [
            coarse_source_impulse.identity_residual_n_s,
            fine_source_impulse.identity_residual_n_s,
            coarse_published_impulse.identity_residual_n_s,
            fine_published_impulse.identity_residual_n_s,
        ]
        .into_iter()
        .fold(0.0_f64, f64::max),
        maximum_impulse_identity_tolerance_n_s: [
            coarse_source_impulse.identity_tolerance_n_s,
            fine_source_impulse.identity_tolerance_n_s,
            coarse_published_impulse.identity_tolerance_n_s,
            fine_published_impulse.identity_tolerance_n_s,
        ]
        .into_iter()
        .fold(0.0_f64, f64::max),
        maximum_center_of_mass_difference_m,
        maximum_orientation_difference_rad,
        maximum_chirp_difference_hz,
        maximum_relative_chirp_difference,
    };
    admit_fixture_output_convergence(evidence)?;
    Ok(evidence)
}

fn require_trajectory_horizon(
    label: &str,
    trajectory: &RenderTrajectory,
    expected_duration_s: f64,
) -> Result<(), CinematicFixtureError> {
    let samples = trajectory.samples();
    let first_time_s = samples
        .first()
        .expect("validated trajectory is nonempty")
        .input()
        .time_s;
    let last_time_s = samples
        .last()
        .expect("validated trajectory is nonempty")
        .input()
        .time_s;
    let tolerance_s = 32.0 * f64::EPSILON * expected_duration_s.max(1.0);
    if first_time_s.abs() > tolerance_s || (last_time_s - expected_duration_s).abs() > tolerance_s {
        return Err(CinematicFixtureError::Pipeline(format!(
            "{label} horizon [{first_time_s:.17e}, {last_time_s:.17e}] s does not match [0, {expected_duration_s:.17e}] s"
        )));
    }
    Ok(())
}

fn fixture_impulse_audit(
    trajectory: &RenderTrajectory,
    query_times_s: &[f64],
    gravity_m_per_s2: f64,
    cx: &Cx<'_>,
) -> Result<FixtureImpulseAudit, CinematicFixtureError> {
    let controls = EulerControlStream::try_derive(trajectory, cx).map_err(pipeline)?;
    let intervals = controls.audio();
    let mut interval_measures_n_s = Vec::new();
    interval_measures_n_s
        .try_reserve_exact(intervals.len())
        .map_err(|_| {
            CinematicFixtureError::Pipeline("normal-impulse audit allocation refused".into())
        })?;
    for interval in intervals {
        let full_contact_mean_available = interval.channels.contact.available().is_some();
        let declared_scalar_is_interval_mean = matches!(
            interval.normal_force_sampling,
            RenderNormalForceSampling::IntervalMean
                | RenderNormalForceSampling::AppliedSubstepZeroOrderHold
        );
        if !full_contact_mean_available && !declared_scalar_is_interval_mean {
            return Err(CinematicFixtureError::Pipeline(format!(
                "normal-impulse audit requires a full-contact duration mean, IntervalMean scalar, or applied zero-order hold at source sample {}",
                interval.source_sample_index
            )));
        }
        let mean_normal_force_n = interval.mean_base_normal_contact_force_n.ok_or_else(|| {
            CinematicFixtureError::Pipeline(format!(
                "normal-impulse audit is missing the mean normal load at source sample {}",
                interval.source_sample_index
            ))
        })?;
        let measure_n_s = mean_normal_force_n * interval.duration_s;
        if !(mean_normal_force_n.is_finite()
            && mean_normal_force_n >= 0.0
            && measure_n_s.is_finite()
            && measure_n_s >= 0.0)
        {
            return Err(CinematicFixtureError::Pipeline(format!(
                "normal-impulse audit found an invalid interval measure at source sample {}",
                interval.source_sample_index
            )));
        }
        interval_measures_n_s.push(measure_n_s);
    }
    let first_time_s = controls
        .visualization()
        .first()
        .expect("validated control stream is nonempty")
        .time_s;
    let last_time_s = controls
        .visualization()
        .last()
        .expect("validated control stream is nonempty")
        .time_s;
    let horizon_tolerance_s = 32.0 * f64::EPSILON * (last_time_s - first_time_s).abs().max(1.0);
    let mut canonical_query_times_s = Vec::new();
    canonical_query_times_s
        .try_reserve_exact(query_times_s.len())
        .map_err(|_| {
            CinematicFixtureError::Pipeline("normal-impulse query allocation refused".into())
        })?;
    let mut previous_query_s = None;
    for (query_index, query_time_s) in query_times_s.iter().copied().enumerate() {
        let canonical_query_time_s =
            canonicalize_terminal_query_time_s(query_time_s, last_time_s, horizon_tolerance_s);
        if !canonical_query_time_s.is_finite()
            || canonical_query_time_s < first_time_s
            || canonical_query_time_s > last_time_s
            || previous_query_s.is_some_and(|previous| canonical_query_time_s < previous)
        {
            return Err(CinematicFixtureError::Pipeline(format!(
                "normal-impulse query {query_index} at {query_time_s:.17e}s is outside the ordered trajectory horizon"
            )));
        }
        canonical_query_times_s.push(canonical_query_time_s);
        previous_query_s = Some(canonical_query_time_s);
    }

    let mut completed = CompensatedSum::default();
    let mut interval_index = 0_usize;
    let mut cumulative_at_queries_n_s = Vec::new();
    cumulative_at_queries_n_s
        .try_reserve_exact(query_times_s.len())
        .map_err(|_| {
            CinematicFixtureError::Pipeline("cumulative-impulse audit allocation refused".into())
        })?;
    for query_time_s in &canonical_query_times_s {
        while interval_index < intervals.len()
            && intervals[interval_index].end_time_s <= *query_time_s
        {
            completed.add(interval_measures_n_s[interval_index]);
            interval_index += 1;
        }
        let mut cumulative_n_s = completed.sum;
        if let Some(interval) = intervals.get(interval_index) {
            if interval.start_time_s < *query_time_s {
                let overlap_s = (*query_time_s).min(interval.end_time_s) - interval.start_time_s;
                let mean_normal_force_n = interval
                    .mean_base_normal_contact_force_n
                    .expect("validated interval retains a mean normal load");
                cumulative_n_s += mean_normal_force_n * overlap_s;
            }
        }
        cumulative_at_queries_n_s.push(cumulative_n_s);
    }

    let mut total = CompensatedSum::default();
    for measure_n_s in &interval_measures_n_s {
        total.add(*measure_n_s);
    }
    let first_velocity_z_m_per_s = controls
        .visualization()
        .first()
        .expect("validated control stream is nonempty")
        .center_of_mass_velocity_world_m_per_s
        .z;
    let last_velocity_z_m_per_s = controls
        .visualization()
        .last()
        .expect("validated control stream is nonempty")
        .center_of_mass_velocity_world_m_per_s
        .z;
    let mass_kg = trajectory.metadata().mass_properties.properties.mass();
    let horizon_s = last_time_s - first_time_s;
    let expected_n_s = mass_kg
        * (gravity_m_per_s2 * horizon_s + last_velocity_z_m_per_s - first_velocity_z_m_per_s);
    let identity_residual_n_s = (total.sum - expected_n_s).abs();
    let identity_scale_n_s = total
        .sum
        .abs()
        .max(expected_n_s.abs())
        .max((mass_kg * gravity_m_per_s2 * horizon_s).abs())
        .max(f64::MIN_POSITIVE);
    let identity_tolerance_n_s = OUTPUT_CONVERGENCE_IMPULSE_ROUNDOFF_PER_INTERVAL
        * f64::EPSILON
        * (intervals.len().max(1) as f64)
        * identity_scale_n_s;
    Ok(FixtureImpulseAudit {
        cumulative_at_queries_n_s,
        total_n_s: total.sum,
        identity_residual_n_s,
        identity_tolerance_n_s,
    })
}

fn canonicalize_terminal_query_time_s(
    query_time_s: f64,
    last_time_s: f64,
    horizon_tolerance_s: f64,
) -> f64 {
    if query_time_s > last_time_s && query_time_s - last_time_s <= horizon_tolerance_s {
        last_time_s
    } else {
        query_time_s
    }
}

fn scaled_absolute_difference(coarse: f64, fine: f64, physical_scale: f64) -> f64 {
    (coarse - fine).abs() / fine.abs().max(physical_scale.abs()).max(f64::MIN_POSITIVE)
}

fn quaternion_geodesic_difference_rad(coarse_wxyz: [f64; 4], fine_wxyz: [f64; 4]) -> f64 {
    let absolute_dot = coarse_wxyz
        .iter()
        .zip(fine_wxyz)
        .map(|(coarse, fine)| coarse * fine)
        .sum::<f64>()
        .abs()
        .clamp(0.0, 1.0);
    2.0 * det::acos(absolute_dot)
}

fn body_contact_chirp_hz(
    state: fs_mbd::RigidBodyState,
    mass: fs_mbd::MassProperties,
) -> Result<f64, CinematicFixtureError> {
    let qois = DerivedEulerQois::from_state(state, mass, 0.0).map_err(pipeline)?;
    let frequency_hz =
        qois.precession_rad_per_s * det::cos(qois.inclination_rad) / core::f64::consts::TAU;
    if !(frequency_hz.is_finite() && frequency_hz > 0.0) {
        return Err(CinematicFixtureError::Pipeline(
            "output convergence produced an invalid body-contact chirp frequency".into(),
        ));
    }
    Ok(frequency_hz)
}

fn admit_fixture_output_convergence(
    evidence: FixtureOutputConvergenceEvidence,
) -> Result<(), CinematicFixtureError> {
    let finite_nonnegative = [
        evidence.terminal_time_difference_s,
        evidence.terminal_time_relative_difference,
        evidence.coarse_preroll_normal_impulse_n_s,
        evidence.fine_preroll_normal_impulse_n_s,
        evidence.preroll_normal_impulse_relative_difference,
        evidence.coarse_source_normal_impulse_n_s,
        evidence.fine_source_normal_impulse_n_s,
        evidence.source_normal_impulse_relative_difference,
        evidence.coarse_published_normal_impulse_n_s,
        evidence.fine_published_normal_impulse_n_s,
        evidence.published_normal_impulse_relative_difference,
        evidence.maximum_cumulative_normal_impulse_difference_n_s,
        evidence.maximum_cumulative_normal_impulse_relative_difference,
        evidence.maximum_impulse_identity_residual_n_s,
        evidence.maximum_impulse_identity_tolerance_n_s,
        evidence.maximum_center_of_mass_difference_m,
        evidence.maximum_orientation_difference_rad,
        evidence.maximum_chirp_difference_hz,
        evidence.maximum_relative_chirp_difference,
    ];
    if finite_nonnegative
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || evidence.shutter_pose_sample_count == 0
        || evidence
            .shutter_exact_renderer_sample_count
            .checked_add(evidence.shutter_stratum_boundary_sample_count)
            .and_then(|count| count.checked_add(evidence.shutter_interpolation_knot_sample_count))
            .and_then(|count| {
                count.checked_add(evidence.shutter_interpolation_midpoint_sample_count)
            })
            != Some(evidence.shutter_pose_sample_count)
    {
        return Err(CinematicFixtureError::Pipeline(
            "output-space dt/dt2 convergence evidence is non-finite, negative, empty, or internally inconsistent"
                .into(),
        ));
    }
    let impulse_maximum = evidence
        .preroll_normal_impulse_relative_difference
        .max(evidence.source_normal_impulse_relative_difference)
        .max(evidence.published_normal_impulse_relative_difference)
        .max(evidence.maximum_cumulative_normal_impulse_relative_difference);
    if evidence.terminal_time_difference_s > OUTPUT_CONVERGENCE_TERMINAL_TIME_ABSOLUTE_LIMIT_S
        || evidence.terminal_time_relative_difference
            > OUTPUT_CONVERGENCE_TERMINAL_TIME_RELATIVE_LIMIT
        || impulse_maximum > OUTPUT_CONVERGENCE_RELATIVE_IMPULSE_LIMIT
        || evidence.maximum_center_of_mass_difference_m > OUTPUT_CONVERGENCE_COM_LIMIT_M
        || evidence.maximum_orientation_difference_rad > OUTPUT_CONVERGENCE_ORIENTATION_LIMIT_RAD
        || evidence.maximum_chirp_difference_hz > OUTPUT_CONVERGENCE_CHIRP_ABSOLUTE_LIMIT_HZ
        || evidence.maximum_relative_chirp_difference > OUTPUT_CONVERGENCE_CHIRP_RELATIVE_LIMIT
    {
        return Err(CinematicFixtureError::Pipeline(format!(
            concat!(
                "output-space dt/dt2 convergence refused: terminal_s={terminal:.17e} ",
                "(limit {terminal_limit:.17e}), terminal_rel={terminal_relative:.17e} ",
                "(limit {terminal_relative_limit:.17e}), impulse_rel={impulse_maximum:.17e} ",
                "(limit {impulse_limit:.17e}), com_m={com:.17e} (limit {com_limit:.17e}), ",
                "orientation_rad={orientation:.17e} (limit {orientation_limit:.17e}), ",
                "chirp_hz={chirp:.17e} (limit {chirp_limit:.17e}), ",
                "chirp_rel={chirp_relative:.17e} (limit {chirp_relative_limit:.17e})"
            ),
            terminal = evidence.terminal_time_difference_s,
            terminal_limit = OUTPUT_CONVERGENCE_TERMINAL_TIME_ABSOLUTE_LIMIT_S,
            terminal_relative = evidence.terminal_time_relative_difference,
            terminal_relative_limit = OUTPUT_CONVERGENCE_TERMINAL_TIME_RELATIVE_LIMIT,
            impulse_maximum = impulse_maximum,
            impulse_limit = OUTPUT_CONVERGENCE_RELATIVE_IMPULSE_LIMIT,
            com = evidence.maximum_center_of_mass_difference_m,
            com_limit = OUTPUT_CONVERGENCE_COM_LIMIT_M,
            orientation = evidence.maximum_orientation_difference_rad,
            orientation_limit = OUTPUT_CONVERGENCE_ORIENTATION_LIMIT_RAD,
            chirp = evidence.maximum_chirp_difference_hz,
            chirp_limit = OUTPUT_CONVERGENCE_CHIRP_ABSOLUTE_LIMIT_HZ,
            chirp_relative = evidence.maximum_relative_chirp_difference,
            chirp_relative_limit = OUTPUT_CONVERGENCE_CHIRP_RELATIVE_LIMIT,
        )));
    }
    Ok(())
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
        CameraProjection::try_half_tangent(0.28)?,
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

fn map_all_audio_intervals(
    mapper: &AudioExcitationMapper<'_, '_>,
    cx: &Cx<'_>,
) -> Result<Vec<crate::AudioExcitationInterval>, CinematicFixtureError> {
    let selected_count = mapper.grid().interval_count;
    let mut mapped = Vec::new();
    mapped
        .try_reserve_exact(selected_count)
        .map_err(|_| CinematicFixtureError::Pipeline("audio interval allocation refused".into()))?;
    let mut checkpoint = mapper.initial_checkpoint(cx).map_err(pipeline)?;
    while checkpoint.next_interval_index() < selected_count {
        let remaining = selected_count - checkpoint.next_interval_index();
        let chunk_size = remaining.min(65_536);
        let chunk = mapper
            .map_next_chunk(
                &checkpoint,
                NonZeroUsize::new(chunk_size).expect("positive remaining interval count"),
                cx,
            )
            .map_err(pipeline)?;
        mapped.extend(chunk.intervals);
        checkpoint = chunk.successor;
    }
    if mapped.len() != selected_count {
        return Err(CinematicFixtureError::Pipeline(format!(
            "audio mapping retained {} of {selected_count} source intervals",
            mapped.len()
        )));
    }
    Ok(mapped)
}

fn fixture_resampling_input(
    video_clock: CinematicClock,
    audio_clock: CinematicClock,
) -> AudioResamplingModelInput {
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
    }
}

fn fixture_timeline_identity(
    frames: u32,
    audio_frame_count: u64,
    warm_start_checkpoint_identity: Option<ContentHash>,
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.master-clocks.v3");
    hasher.update(&CRITIQUE_FPS.to_le_bytes());
    hasher.update(&SOUND_MASTER_SAMPLE_RATE_HZ.to_le_bytes());
    hasher.update(&0_i64.to_le_bytes());
    hasher.update(&frames.to_le_bytes());
    hasher.update(&0_i64.to_le_bytes());
    hasher.update(&audio_frame_count.to_le_bytes());
    if let Some(identity) = warm_start_checkpoint_identity {
        hasher.update(AUDIO_PREROLL_POLICY_ID.as_bytes());
        hasher.update(&AUDIO_PREROLL_SAMPLE_FRAMES.to_le_bytes());
        hasher.update(identity.as_bytes());
    } else {
        hasher.update(b"zero-state-origin");
    }
    hasher.finalize()
}

#[allow(clippy::too_many_arguments)]
fn admit_fixture_sound(
    trajectory: &EulerRenderTrajectoryArtifact,
    excitation_identity: ContentHash,
    modal: &ModalSynthesisModel,
    resampler_identity: ContentHash,
    filter_identity: ContentHash,
    video_clock: CinematicClock,
    audio_clock: CinematicClock,
    timeline_identity: ContentHash,
    mappings: Vec<SoundExcitationControl>,
) -> Result<SoundSynthesisConfig, CinematicFixtureError> {
    let listener = ListenerPose {
        frame: ListenerFrame::AnimatedCamera,
        position_m: [0.0, 0.0, 0.0],
        forward: [0.0, 0.0, -1.0],
        up: [0.0, 1.0, 0.0],
    };
    SoundSynthesisConfig::try_admit(SoundSynthesisInput {
        schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
        authority: SoundAuthority::PhysicallyInformed,
        trajectory: component(
            CinematicComponentRole::Trajectory,
            trajectory.receipt().artifact_identity(),
            u32::from(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION),
        )?,
        excitation: component(
            CinematicComponentRole::AudioExcitation,
            excitation_identity,
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
        timeline: component(CinematicComponentRole::Timeline, timeline_identity, 3)?,
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
        resampler_identity,
        resampler_version: AUDIO_RESAMPLING_ALGORITHM_VERSION,
        filter_identity,
        filter_version: AUDIO_RECONSTRUCTION_FILTER_VERSION,
        assumptions: vec![
            SoundModelAssumption::LinearModalSuperposition,
            SoundModelAssumption::TimeInvariantDamping,
            SoundModelAssumption::DeclaredExcitationCompleteness,
            SoundModelAssumption::DeclaredRoomResponse,
        ],
        calibration: None,
    })
    .map_err(pipeline)
}

#[allow(clippy::too_many_arguments)]
fn prepare_fixture_audio_candidate(
    source: &EulerRenderTrajectoryArtifact,
    modal: &ModalSynthesisModel,
    mappings: &[SoundExcitationControl],
    spatial_rules: &[ModeContactParticipationRule],
    video_clock: CinematicClock,
    audio_clock: CinematicClock,
    video_frame_count: u32,
    audio_frame_count: u64,
    cx: &Cx<'_>,
) -> Result<FixtureAudioCandidate, CinematicFixtureError> {
    let controls = EulerControlStream::try_derive(source.trajectory(), cx).map_err(pipeline)?;
    let interval_count = controls.audio().len();
    let mapper = AudioExcitationMapper::try_new(
        source,
        &controls,
        modal,
        AudioExcitationModelInput {
            mappings: mappings.to_vec(),
            reduction: AudioExcitationReduction::RawIntervals,
            spatial_policy: ContactParticipationPolicy::ContactCoordinates {
                rules: spatial_rules.to_vec(),
            },
            // The convergence comparison must not contain a stochastic or
            // authored impulse bank that can obscure mechanics refinement.
            artistic_texture: None,
            budget: AudioExcitationBudget::reference_film(interval_count),
        },
        cx,
    )
    .map_err(pipeline)?;
    let intervals = map_all_audio_intervals(&mapper, cx)?;
    let resampler = AudioResampler::try_new(
        &mapper,
        modal,
        intervals,
        fixture_resampling_input(video_clock, audio_clock),
        cx,
    )
    .map_err(pipeline)?;
    if resampler.total_audio_frames() != audio_frame_count {
        return Err(CinematicFixtureError::Pipeline(format!(
            "audio candidate admitted {} frames, expected {audio_frame_count}",
            resampler.total_audio_frames()
        )));
    }
    let excitation_identity = mapper.identity();
    let source_sound = admit_fixture_sound(
        source,
        excitation_identity,
        modal,
        resampler.identity(),
        resampler.filter_identity(),
        video_clock,
        audio_clock,
        fixture_timeline_identity(video_frame_count, audio_frame_count, None),
        mappings.to_vec(),
    )?;
    mapper
        .validate_sound_configuration(&source_sound)
        .map_err(pipeline)?;
    modal
        .validate_sound_configuration(&source_sound)
        .map_err(pipeline)?;
    resampler
        .validate_sound_configuration(&source_sound)
        .map_err(pipeline)?;
    Ok(FixtureAudioCandidate {
        resampler,
        source_sound,
        excitation_identity,
        source_trajectory_identity: source.receipt().artifact_identity(),
    })
}

#[allow(clippy::too_many_arguments)]
fn validate_audio_pair_progress(
    stage: &'static str,
    expected_start: u64,
    expected_end: u64,
    coarse_start: u64,
    coarse_end: u64,
    coarse_successor: u64,
    fine_start: u64,
    fine_end: u64,
    fine_successor: u64,
) -> Result<(), CinematicFixtureError> {
    if coarse_start != expected_start
        || coarse_end != expected_end
        || coarse_successor != expected_end
        || fine_start != expected_start
        || fine_end != expected_end
        || fine_successor != expected_end
    {
        return Err(CinematicFixtureError::Pipeline(format!(
            "{stage} coarse/fine checkpoint range mismatch: expected [{expected_start},{expected_end}), coarse [{coarse_start},{coarse_end}) -> {coarse_successor}, fine [{fine_start},{fine_end}) -> {fine_successor}"
        )));
    }
    Ok(())
}

fn observe_resampled_drive_pair(
    coarse: &AudioResamplingChunk,
    fine: &AudioResamplingChunk,
    mode_count: usize,
    localized: &mut [SymmetricErrorAccumulator],
    distributed: &mut [SymmetricErrorAccumulator; 3],
) -> Result<(), CinematicFixtureError> {
    if mode_count == 0 || localized.len() != mode_count {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence requires one accumulator per canonical mode".into(),
        ));
    }
    if coarse.drive_frames.len() != fine.drive_frames.len() {
        return Err(CinematicFixtureError::Pipeline(
            "coarse/fine resampled drive frame counts differ".into(),
        ));
    }
    let expected_localized = coarse
        .drive_frames
        .len()
        .checked_mul(mode_count)
        .ok_or_else(|| CinematicFixtureError::Pipeline("localized drive length overflow".into()))?;
    if coarse.preparticipated_localized_force_n.len() != expected_localized
        || fine.preparticipated_localized_force_n.len() != expected_localized
        || coarse.preparticipated_localized_impulse_n_s.len() != expected_localized
        || fine.preparticipated_localized_impulse_n_s.len() != expected_localized
    {
        return Err(CinematicFixtureError::Pipeline(format!(
            "coarse/fine localized drive arrays do not match expected row-major length {expected_localized}"
        )));
    }
    if coarse
        .preparticipated_localized_impulse_n_s
        .iter()
        .chain(&fine.preparticipated_localized_impulse_n_s)
        .any(|value| !value.is_finite() || *value != 0.0)
    {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence requires exact all-zero localized impulse banks in both members"
                .into(),
        ));
    }
    for (index, (&coarse_force, &fine_force)) in coarse
        .preparticipated_localized_force_n
        .iter()
        .zip(&fine.preparticipated_localized_force_n)
        .enumerate()
    {
        localized[index % mode_count].observe(coarse_force, fine_force)?;
    }
    for (coarse_frame, fine_frame) in coarse.drive_frames.iter().zip(&fine.drive_frames) {
        let coarse_values = coarse_frame.distributed_generalized_force_n;
        let fine_values = fine_frame.distributed_generalized_force_n;
        for (accumulator, (coarse_value, fine_value)) in distributed.iter_mut().zip([
            (coarse_values.disc, fine_values.disc),
            (coarse_values.glass_plate, fine_values.glass_plate),
            (coarse_values.base_assembly, fine_values.base_assembly),
        ]) {
            accumulator.observe(coarse_value, fine_value)?;
        }
    }
    Ok(())
}

fn observe_stem_pair(
    coarse: &[ModalStemFrame],
    fine: &[ModalStemFrame],
    stems: &mut [SymmetricErrorAccumulator; 3],
) -> Result<(), CinematicFixtureError> {
    if coarse.len() != fine.len() {
        return Err(CinematicFixtureError::Pipeline(
            "coarse/fine raw modal stem frame counts differ".into(),
        ));
    }
    for (coarse_frame, fine_frame) in coarse.iter().zip(fine) {
        for (accumulator, (coarse_value, fine_value)) in stems.iter_mut().zip([
            (coarse_frame.disc_fs, fine_frame.disc_fs),
            (coarse_frame.glass_plate_fs, fine_frame.glass_plate_fs),
            (coarse_frame.base_assembly_fs, fine_frame.base_assembly_fs),
        ]) {
            accumulator.observe(coarse_value, fine_value)?;
        }
    }
    Ok(())
}

fn maximum_error_metrics(
    accumulators: &[SymmetricErrorAccumulator],
    floor: f64,
) -> Result<RelativeErrorMetrics, CinematicFixtureError> {
    if accumulators.is_empty() {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence metric bank is empty".into(),
        ));
    }
    let mut maximum = RelativeErrorMetrics {
        nrmse: 0.0,
        normalized_peak: 0.0,
    };
    for accumulator in accumulators {
        let metric = accumulator.metrics(floor)?;
        maximum.nrmse = maximum.nrmse.max(metric.nrmse);
        maximum.normalized_peak = maximum.normalized_peak.max(metric.normalized_peak);
    }
    Ok(maximum)
}

fn enforce_audio_convergence_thresholds(
    localized: RelativeErrorMetrics,
    distributed: RelativeErrorMetrics,
    stems: RelativeErrorMetrics,
) -> Result<(), CinematicFixtureError> {
    let drive_nrmse = localized.nrmse.max(distributed.nrmse);
    let drive_peak = localized.normalized_peak.max(distributed.normalized_peak);
    if drive_nrmse > AUDIO_CONVERGENCE_DRIVE_NRMSE_LIMIT
        || drive_peak > AUDIO_CONVERGENCE_DRIVE_NORMALIZED_PEAK_LIMIT
        || stems.nrmse > AUDIO_CONVERGENCE_STEM_NRMSE_LIMIT
        || stems.normalized_peak > AUDIO_CONVERGENCE_STEM_NORMALIZED_PEAK_LIMIT
    {
        return Err(CinematicFixtureError::Pipeline(format!(
            concat!(
                "pre-master audio dt/dt2 convergence refused: drive NRMSE {drive_nrmse:.17e} ",
                "(limit {drive_nrmse_limit:.17e}), drive normalized peak {drive_peak:.17e} ",
                "(limit {drive_peak_limit:.17e}), raw cropped stem NRMSE {stem_nrmse:.17e} ",
                "(limit {stem_nrmse_limit:.17e}), raw cropped stem normalized peak ",
                "{stem_peak:.17e} (limit {stem_peak_limit:.17e})"
            ),
            drive_nrmse = drive_nrmse,
            drive_nrmse_limit = AUDIO_CONVERGENCE_DRIVE_NRMSE_LIMIT,
            drive_peak = drive_peak,
            drive_peak_limit = AUDIO_CONVERGENCE_DRIVE_NORMALIZED_PEAK_LIMIT,
            stem_nrmse = stems.nrmse,
            stem_nrmse_limit = AUDIO_CONVERGENCE_STEM_NRMSE_LIMIT,
            stem_peak = stems.normalized_peak,
            stem_peak_limit = AUDIO_CONVERGENCE_STEM_NORMALIZED_PEAK_LIMIT,
        )));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn continuous_audio_crop_binding_identity(
    source_trajectory_identity: ContentHash,
    published_trajectory_identity: ContentHash,
    source_sound_configuration_identity: ContentHash,
    full_resampler_identity: ContentHash,
    crop: AudioResamplingCrop,
    first_chunk_identity: ContentHash,
    crop_start_resampling_checkpoint_identity: ContentHash,
    crop_start_modal_checkpoint_identity: ContentHash,
    end_resampling_checkpoint_identity: ContentHash,
    end_modal_checkpoint_identity: ContentHash,
) -> ContentHash {
    let mut hasher =
        DomainHasher::new("org.frankensim.euler-critique.continuous-audio-crop-binding.v3");
    hasher.update(AUDIO_PREROLL_POLICY_ID.as_bytes());
    hasher.update(source_trajectory_identity.as_bytes());
    hasher.update(published_trajectory_identity.as_bytes());
    hasher.update(source_sound_configuration_identity.as_bytes());
    hasher.update(full_resampler_identity.as_bytes());
    hasher.update(crop.identity().as_bytes());
    hasher.update(&crop.first_source_audio_frame().to_le_bytes());
    hasher.update(&crop.end_source_audio_frame().to_le_bytes());
    hasher.update(first_chunk_identity.as_bytes());
    hasher.update(crop_start_resampling_checkpoint_identity.as_bytes());
    hasher.update(crop_start_modal_checkpoint_identity.as_bytes());
    hasher.update(end_resampling_checkpoint_identity.as_bytes());
    hasher.update(end_modal_checkpoint_identity.as_bytes());
    hasher.finalize()
}

#[allow(clippy::too_many_arguments)]
fn synthesize_audio_convergence_pair(
    coarse_source: &EulerRenderTrajectoryArtifact,
    fine_source: &EulerRenderTrajectoryArtifact,
    fine_published: &EulerRenderTrajectoryArtifact,
    coarse: &FixtureAudioCandidate,
    fine: &FixtureAudioCandidate,
    modal: &ModalSynthesisModel,
    mappings: &[SoundExcitationControl],
    full_audio_frame_count: u64,
    published_video_frame_count: u32,
    published_audio_frame_count: u64,
    output_video_clock: CinematicClock,
    output_audio_clock: CinematicClock,
    drive_normalization_floor_n: f64,
    cx: &Cx<'_>,
) -> Result<FixtureAudioPairOutput, CinematicFixtureError> {
    if !(drive_normalization_floor_n.is_finite() && drive_normalization_floor_n > 0.0) {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence force floor must be finite and positive".into(),
        ));
    }
    if full_audio_frame_count.checked_sub(published_audio_frame_count)
        != Some(AUDIO_PREROLL_SAMPLE_FRAMES)
    {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence full and published horizons do not share the exact preroll".into(),
        ));
    }
    let mode_count = modal.modes().len();
    if mode_count == 0 {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence modal bank is empty".into(),
        ));
    }
    let mut localized = vec![SymmetricErrorAccumulator::default(); mode_count];
    let mut distributed = [SymmetricErrorAccumulator::default(); 3];
    let mut stem_errors = [SymmetricErrorAccumulator::default(); 3];

    let coarse_initial_resampling = coarse.resampler.initial_checkpoint(cx).map_err(pipeline)?;
    let fine_initial_resampling = fine.resampler.initial_checkpoint(cx).map_err(pipeline)?;
    let coarse_initial_modal = modal.initial_checkpoint(cx).map_err(pipeline)?;
    let fine_initial_modal = modal.initial_checkpoint(cx).map_err(pipeline)?;
    let preroll_chunk_frames = NonZeroUsize::new(
        usize::try_from(AUDIO_PREROLL_SAMPLE_FRAMES)
            .expect("audio preroll frame count is representable as usize"),
    )
    .expect("audio preroll is nonzero");
    let coarse_preroll = coarse
        .resampler
        .resample_next_chunk(
            &coarse.source_sound,
            &coarse_initial_resampling,
            preroll_chunk_frames,
            cx,
        )
        .map_err(pipeline)?;
    let fine_preroll = fine
        .resampler
        .resample_next_chunk(
            &fine.source_sound,
            &fine_initial_resampling,
            preroll_chunk_frames,
            cx,
        )
        .map_err(pipeline)?;
    validate_audio_pair_progress(
        "resampling preroll",
        0,
        AUDIO_PREROLL_SAMPLE_FRAMES,
        coarse_preroll.diagnostics.start_audio_frame_offset,
        coarse_preroll.diagnostics.end_audio_frame_offset,
        coarse_preroll.successor.next_audio_frame_offset(),
        fine_preroll.diagnostics.start_audio_frame_offset,
        fine_preroll.diagnostics.end_audio_frame_offset,
        fine_preroll.successor.next_audio_frame_offset(),
    )?;
    observe_resampled_drive_pair(
        &coarse_preroll,
        &fine_preroll,
        mode_count,
        &mut localized,
        &mut distributed,
    )?;
    let coarse_warmed = coarse_preroll
        .synthesize_modal(modal, &coarse_initial_modal, cx)
        .map_err(pipeline)?;
    let fine_warmed = fine_preroll
        .synthesize_modal(modal, &fine_initial_modal, cx)
        .map_err(pipeline)?;
    validate_audio_pair_progress(
        "modal preroll",
        0,
        AUDIO_PREROLL_SAMPLE_FRAMES,
        coarse_warmed.diagnostics.start_sample_frame,
        coarse_warmed.diagnostics.end_sample_frame,
        coarse_warmed.successor.next_sample_frame(),
        fine_warmed.diagnostics.start_sample_frame,
        fine_warmed.diagnostics.end_sample_frame,
        fine_warmed.successor.next_sample_frame(),
    )?;
    let coarse_crop_start_resampling_checkpoint_identity = coarse_preroll.successor.identity();
    let fine_crop_start_resampling_checkpoint_identity = fine_preroll.successor.identity();
    let coarse_crop_start_modal_checkpoint_identity = coarse_warmed.successor.identity();
    let fine_crop_start_modal_checkpoint_identity = fine_warmed.successor.identity();
    let coarse_first_chunk_identity = coarse_preroll.identity;
    let fine_first_chunk_identity = fine_preroll.identity;
    let mut coarse_resampling_checkpoint = coarse_preroll.successor;
    let mut fine_resampling_checkpoint = fine_preroll.successor;
    let mut coarse_modal_checkpoint = coarse_warmed.successor;
    let mut fine_modal_checkpoint = fine_warmed.successor;
    let published_capacity = usize::try_from(published_audio_frame_count)
        .map_err(|_| CinematicFixtureError::Pipeline("audio stem length overflow".into()))?;
    let mut fine_stems = Vec::new();
    fine_stems
        .try_reserve_exact(published_capacity)
        .map_err(|_| CinematicFixtureError::Pipeline("audio stem allocation refused".into()))?;

    while fine_resampling_checkpoint.next_audio_frame_offset() < full_audio_frame_count {
        let expected_start = fine_resampling_checkpoint.next_audio_frame_offset();
        if coarse_resampling_checkpoint.next_audio_frame_offset() != expected_start
            || coarse_modal_checkpoint.next_sample_frame() != expected_start
            || fine_modal_checkpoint.next_sample_frame() != expected_start
        {
            return Err(CinematicFixtureError::Pipeline(
                "coarse/fine audio checkpoints lost lockstep before a chunk".into(),
            ));
        }
        let remaining = full_audio_frame_count - expected_start;
        let chunk_frames = usize::try_from(remaining.min(65_536))
            .map_err(|_| CinematicFixtureError::Pipeline("audio chunk length overflow".into()))?;
        let expected_end = expected_start
            .checked_add(chunk_frames as u64)
            .ok_or_else(|| CinematicFixtureError::Pipeline("audio chunk end overflow".into()))?;
        let maximum_frames =
            NonZeroUsize::new(chunk_frames).expect("positive remaining audio frame count");
        let coarse_resampled = coarse
            .resampler
            .resample_next_chunk(
                &coarse.source_sound,
                &coarse_resampling_checkpoint,
                maximum_frames,
                cx,
            )
            .map_err(pipeline)?;
        let fine_resampled = fine
            .resampler
            .resample_next_chunk(
                &fine.source_sound,
                &fine_resampling_checkpoint,
                maximum_frames,
                cx,
            )
            .map_err(pipeline)?;
        validate_audio_pair_progress(
            "resampling",
            expected_start,
            expected_end,
            coarse_resampled.diagnostics.start_audio_frame_offset,
            coarse_resampled.diagnostics.end_audio_frame_offset,
            coarse_resampled.successor.next_audio_frame_offset(),
            fine_resampled.diagnostics.start_audio_frame_offset,
            fine_resampled.diagnostics.end_audio_frame_offset,
            fine_resampled.successor.next_audio_frame_offset(),
        )?;
        observe_resampled_drive_pair(
            &coarse_resampled,
            &fine_resampled,
            mode_count,
            &mut localized,
            &mut distributed,
        )?;
        let coarse_synthesized = coarse_resampled
            .synthesize_modal(modal, &coarse_modal_checkpoint, cx)
            .map_err(pipeline)?;
        let fine_synthesized = fine_resampled
            .synthesize_modal(modal, &fine_modal_checkpoint, cx)
            .map_err(pipeline)?;
        validate_audio_pair_progress(
            "modal synthesis",
            expected_start,
            expected_end,
            coarse_synthesized.diagnostics.start_sample_frame,
            coarse_synthesized.diagnostics.end_sample_frame,
            coarse_synthesized.successor.next_sample_frame(),
            fine_synthesized.diagnostics.start_sample_frame,
            fine_synthesized.diagnostics.end_sample_frame,
            fine_synthesized.successor.next_sample_frame(),
        )?;
        observe_stem_pair(
            &coarse_synthesized.stem_frames,
            &fine_synthesized.stem_frames,
            &mut stem_errors,
        )?;
        fine_stems.extend(fine_synthesized.stem_frames);
        coarse_resampling_checkpoint = coarse_resampled.successor;
        fine_resampling_checkpoint = fine_resampled.successor;
        coarse_modal_checkpoint = coarse_synthesized.successor;
        fine_modal_checkpoint = fine_synthesized.successor;
    }
    if fine_stems.len() != published_capacity {
        return Err(CinematicFixtureError::Pipeline(format!(
            "fine modal synthesis retained {} of {published_audio_frame_count} published frames",
            fine_stems.len()
        )));
    }
    validate_audio_pair_progress(
        "full-horizon checkpoints",
        full_audio_frame_count,
        full_audio_frame_count,
        coarse_resampling_checkpoint.next_audio_frame_offset(),
        coarse_modal_checkpoint.next_sample_frame(),
        coarse_modal_checkpoint.next_sample_frame(),
        fine_resampling_checkpoint.next_audio_frame_offset(),
        fine_modal_checkpoint.next_sample_frame(),
        fine_modal_checkpoint.next_sample_frame(),
    )?;

    let coarse_crop = coarse
        .resampler
        .try_crop(
            AUDIO_PREROLL_SAMPLE_FRAMES,
            full_audio_frame_count,
            output_video_clock,
            output_audio_clock,
        )
        .map_err(pipeline)?;
    let fine_crop = fine
        .resampler
        .try_crop(
            AUDIO_PREROLL_SAMPLE_FRAMES,
            full_audio_frame_count,
            output_video_clock,
            output_audio_clock,
        )
        .map_err(pipeline)?;
    if coarse_crop.first_source_audio_frame() != fine_crop.first_source_audio_frame()
        || coarse_crop.end_source_audio_frame() != fine_crop.end_source_audio_frame()
        || coarse_crop.output_video_clock() != fine_crop.output_video_clock()
        || coarse_crop.output_audio_clock() != fine_crop.output_audio_clock()
    {
        return Err(CinematicFixtureError::Pipeline(
            "coarse/fine typed audio crops do not share an exact range and output clocks".into(),
        ));
    }
    let coarse_end_resampling_checkpoint_identity = coarse_resampling_checkpoint.identity();
    let fine_end_resampling_checkpoint_identity = fine_resampling_checkpoint.identity();
    let coarse_end_modal_checkpoint_identity = coarse_modal_checkpoint.identity();
    let fine_end_modal_checkpoint_identity = fine_modal_checkpoint.identity();
    let coarse_crop_binding_identity = continuous_audio_crop_binding_identity(
        coarse.source_trajectory_identity,
        coarse.source_trajectory_identity,
        coarse.source_sound.receipt().configuration_identity,
        coarse.resampler.identity(),
        coarse_crop,
        coarse_first_chunk_identity,
        coarse_crop_start_resampling_checkpoint_identity,
        coarse_crop_start_modal_checkpoint_identity,
        coarse_end_resampling_checkpoint_identity,
        coarse_end_modal_checkpoint_identity,
    );
    let fine_crop_binding_identity = continuous_audio_crop_binding_identity(
        fine.source_trajectory_identity,
        fine_published.receipt().artifact_identity(),
        fine.source_sound.receipt().configuration_identity,
        fine.resampler.identity(),
        fine_crop,
        fine_first_chunk_identity,
        fine_crop_start_resampling_checkpoint_identity,
        fine_crop_start_modal_checkpoint_identity,
        fine_end_resampling_checkpoint_identity,
        fine_end_modal_checkpoint_identity,
    );
    let coarse_crop_sound = admit_fixture_sound(
        coarse_source,
        coarse.excitation_identity,
        modal,
        coarse_crop.identity(),
        coarse.resampler.filter_identity(),
        output_video_clock,
        output_audio_clock,
        fixture_timeline_identity(
            published_video_frame_count,
            published_audio_frame_count,
            Some(coarse_crop_binding_identity),
        ),
        mappings.to_vec(),
    )?;
    let fine_crop_sound = admit_fixture_sound(
        fine_source,
        fine.excitation_identity,
        modal,
        fine_crop.identity(),
        fine.resampler.filter_identity(),
        output_video_clock,
        output_audio_clock,
        fixture_timeline_identity(
            published_video_frame_count,
            published_audio_frame_count,
            Some(fine_crop_binding_identity),
        ),
        mappings.to_vec(),
    )?;
    for sound in [&coarse_crop_sound, &fine_crop_sound] {
        modal
            .validate_sound_configuration(sound)
            .map_err(pipeline)?;
    }
    coarse
        .resampler
        .validate_cropped_sound_configuration(&coarse_crop, &coarse_crop_sound)
        .map_err(pipeline)?;
    fine.resampler
        .validate_cropped_sound_configuration(&fine_crop, &fine_crop_sound)
        .map_err(pipeline)?;

    let localized_metrics = maximum_error_metrics(&localized, drive_normalization_floor_n)?;
    let distributed_metrics = maximum_error_metrics(&distributed, drive_normalization_floor_n)?;
    let stem_metrics = maximum_error_metrics(&stem_errors, AUDIO_CONVERGENCE_STEM_FLOOR_FS)?;
    enforce_audio_convergence_thresholds(localized_metrics, distributed_metrics, stem_metrics)?;
    let convergence = FixtureAudioConvergenceEvidence {
        full_audio_frame_count,
        published_audio_frame_count,
        mode_count,
        drive_normalization_floor_n,
        localized_drive_nrmse: localized_metrics.nrmse,
        localized_drive_normalized_peak: localized_metrics.normalized_peak,
        distributed_drive_nrmse: distributed_metrics.nrmse,
        distributed_drive_normalized_peak: distributed_metrics.normalized_peak,
        maximum_drive_nrmse: localized_metrics.nrmse.max(distributed_metrics.nrmse),
        maximum_drive_normalized_peak: localized_metrics
            .normalized_peak
            .max(distributed_metrics.normalized_peak),
        cropped_stem_nrmse: stem_metrics.nrmse,
        cropped_stem_normalized_peak: stem_metrics.normalized_peak,
        crop_first_source_audio_frame: fine_crop.first_source_audio_frame(),
        crop_end_source_audio_frame: fine_crop.end_source_audio_frame(),
        coarse_crop_identity: coarse_crop.identity(),
        fine_crop_identity: fine_crop.identity(),
        coarse_crop_binding_identity,
        fine_crop_binding_identity,
    };
    Ok(FixtureAudioPairOutput {
        fine_stems,
        fine_sound: fine_crop_sound,
        fine_crop,
        convergence,
    })
}

fn build_audio(
    trajectory: &EulerRenderTrajectoryArtifact,
    coarse_preroll_trajectory: &EulerRenderTrajectoryArtifact,
    fine_preroll_trajectory: &EulerRenderTrajectoryArtifact,
    gravity_m_per_s2: f64,
    config: &CinematicFixtureConfig,
    cx: &Cx<'_>,
) -> Result<FixtureAudio, CinematicFixtureError> {
    let (audio_frame_count, video_clock, audio_clock) = fixture_master_clocks(config.frames)?;
    let preroll_video_frames = config
        .frames
        .checked_add(AUDIO_PREROLL_VIDEO_FRAMES)
        .ok_or_else(|| CinematicFixtureError::Pipeline("audio preroll frame overflow".into()))?;
    let (preroll_audio_frame_count, preroll_video_clock, preroll_audio_clock) =
        fixture_master_clocks(preroll_video_frames)?;
    if preroll_audio_frame_count.checked_sub(audio_frame_count) != Some(AUDIO_PREROLL_SAMPLE_FRAMES)
    {
        return Err(CinematicFixtureError::Pipeline(
            "audio preroll is not exactly one 24 Hz frame at 48 kHz".into(),
        ));
    }
    let (chirp_start_hz, chirp_end_hz) = trajectory_body_contact_chirp_bounds(trajectory)?;
    if chirp_end_hz >= CRITIQUE_DECLARED_AUDIO_SOURCE_BANDWIDTH_HZ {
        return Err(CinematicFixtureError::Pipeline(format!(
            "trajectory-derived chirp reaches {chirp_end_hz:.6} Hz, outside the declared reconstruction bandwidth"
        )));
    }
    let preset = representative_modal_preset(RepresentativeDiscMaterial::StainlessSteel);
    let spatial_rules: Vec<ModeContactParticipationRule> = preset
        .modes()
        .iter()
        .map(|mode| ModeContactParticipationRule {
            mode_id: mode.mode_id,
            shape: if matches!(
                mode.component,
                SoundModalComponent::Disc | SoundModalComponent::GlassPlate
            ) {
                // The moving reaction point sweeps both bodies. This compact
                // harmonic-one ansatz preserves that forcing frequency for
                // the disc and plate mode families; the base modes retain
                // their declared spatially distributed participation.
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
                // One model admits both the warm-start source horizon and the
                // rebased published crop. Its identity therefore binds the
                // longer of those two exact clocks.
                budget: ModalSynthesisBudget::reference_film(preroll_audio_frame_count),
            },
        },
        cx,
    )
    .map_err(pipeline)?;
    let modal_parameter_set_identity = modal_parameters.identity();
    let modal_parameter_set_disclosure = modal_parameters.disclosure().to_owned();
    let modal = modal_parameters.into_model();
    // Drive the two contacting bodies from the retained SI normal reaction.
    // Equal magnitudes and opposite signs encode action/reaction without the
    // previous arbitrary watts-to-newtons transfer. Absolute acoustic output
    // remains uncalibrated because modal radiation and the listening rig are
    // representative rather than measured.
    let mappings = vec![
        SoundExcitationControl {
            channel: SoundExcitationChannel::ContactNormalForce,
            target_component: SoundModalComponent::Disc,
            source_scale: 1.0,
        },
        SoundExcitationControl {
            channel: SoundExcitationChannel::ContactNormalForce,
            target_component: SoundModalComponent::GlassPlate,
            source_scale: -1.0,
        },
    ];
    let coarse_candidate = prepare_fixture_audio_candidate(
        coarse_preroll_trajectory,
        &modal,
        &mappings,
        &spatial_rules,
        preroll_video_clock,
        preroll_audio_clock,
        preroll_video_frames,
        preroll_audio_frame_count,
        cx,
    )?;
    let fine_candidate = prepare_fixture_audio_candidate(
        fine_preroll_trajectory,
        &modal,
        &mappings,
        &spatial_rules,
        preroll_video_clock,
        preroll_audio_clock,
        preroll_video_frames,
        preroll_audio_frame_count,
        cx,
    )?;
    let mass_kg = trajectory
        .trajectory()
        .metadata()
        .mass_properties
        .properties
        .mass();
    if !(mass_kg.is_finite()
        && mass_kg > 0.0
        && gravity_m_per_s2.is_finite()
        && gravity_m_per_s2 > 0.0)
    {
        return Err(CinematicFixtureError::Pipeline(
            "audio convergence mass and gravity reference must be finite and positive".into(),
        ));
    }
    let drive_normalization_floor_n =
        AUDIO_CONVERGENCE_DRIVE_NORMALIZATION_FRACTION * mass_kg * gravity_m_per_s2;
    let pair = synthesize_audio_convergence_pair(
        coarse_preroll_trajectory,
        fine_preroll_trajectory,
        trajectory,
        &coarse_candidate,
        &fine_candidate,
        &modal,
        &mappings,
        preroll_audio_frame_count,
        config.frames,
        audio_frame_count,
        video_clock,
        audio_clock,
        drive_normalization_floor_n,
        cx,
    )?;
    let source_sound_configuration_identity =
        fine_candidate.source_sound.receipt().configuration_identity;
    let warm_start_source_identity = fine_candidate.source_trajectory_identity;
    let FixtureAudioPairOutput {
        fine_stems: mut stems,
        fine_sound: sound,
        fine_crop: crop,
        convergence,
    } = pair;
    let warm_start_checkpoint_identity = convergence.fine_crop_binding_identity;
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
    if !(-MAX_AUDIO_MASTER_GAIN_DB..=MAX_AUDIO_MASTER_GAIN_DB).contains(&master_gain_db) {
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
        warm_start_source_identity,
        warm_start_checkpoint_identity,
        source_sound_configuration_identity,
        published_trajectory_identity: trajectory.receipt().artifact_identity(),
        crop_resampler_identity: crop.identity(),
        crop_first_source_audio_frame: crop.first_source_audio_frame(),
        crop_end_source_audio_frame: crop.end_source_audio_frame(),
        convergence,
        spatialization,
    })
}

fn fixture_master_clocks(
    frames: u32,
) -> Result<(u64, CinematicClock, CinematicClock), CinematicFixtureError> {
    if SOUND_MASTER_SAMPLE_RATE_HZ % CRITIQUE_FPS != 0 {
        return Err(CinematicFixtureError::Pipeline(
            "audio/video master rates do not have an integral frame ratio".into(),
        ));
    }
    let audio_frames_per_video_frame = SOUND_MASTER_SAMPLE_RATE_HZ / CRITIQUE_FPS;
    let audio_frame_count = u64::from(frames)
        .checked_mul(u64::from(audio_frames_per_video_frame))
        .ok_or_else(|| CinematicFixtureError::Pipeline("audio frame count overflow".into()))?;
    let video_clock = CinematicClock::try_new(
        CinematicClockDomain::Video,
        CRITIQUE_FPS,
        1,
        0,
        i64::from(frames),
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
    Ok((audio_frame_count, video_clock, audio_clock))
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
    output_convergence: &FixtureOutputConvergenceEvidence,
    audio_convergence: &FixtureAudioConvergenceEvidence,
    trajectory: &EulerRenderTrajectoryArtifact,
    raw_sequence_identity: ContentHash,
    preview_sequence_identity: ContentHash,
    wav_identity: ContentHash,
    modal_parameter_set_identity: ContentHash,
    modal_parameter_set_disclosure: &str,
    audio_warm_start_source_identity: ContentHash,
    audio_warm_start_checkpoint_identity: ContentHash,
    audio_source_sound_configuration_identity: ContentHash,
    audio_published_trajectory_identity: ContentHash,
    audio_crop_resampler_identity: ContentHash,
    audio_crop_first_source_frame: u64,
    audio_crop_end_source_frame: u64,
    chirp_start_hz: f64,
    chirp_end_hz: f64,
    audio_pre_master_peak_fs: f64,
    audio_master_gain_db: f64,
    over_range_channels: u64,
    gamut_mapped_pixels: u64,
    render: &FixtureRenderEvidence,
    preview_mesh: FixturePreviewMeshEvidence,
    adaptive_sampling: &FixtureAdaptiveSamplingEvidence,
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
    let raw_profile = if config.adaptive_sampling.is_some() && config.retain_full_aov_exr {
        "adaptive-final-diagnostic-aov-float"
    } else if config.adaptive_sampling.is_some() {
        "adaptive-linear-srgb-beauty-plus-sample-count-float"
    } else if config.retain_full_aov_exr {
        "final-diagnostic-aov-float"
    } else {
        "linear-srgb-beauty-float"
    };
    let sampling_json = config.adaptive_sampling.map_or_else(
        || {
            format!(
                "{{\"mode\":\"uniform\",\"samples_per_pixel\":{}}}",
                config.samples_per_pixel
            )
        },
        |policy| adaptive_sampling.manifest_json(policy),
    );
    let denoise_pipeline = if config.denoise_previews {
        TEMPORAL_DENOISE_PIPELINE_VERSION
    } else {
        "disabled"
    };
    let modal_disclosure = json_escape(modal_parameter_set_disclosure);
    let output_convergence_json = output_convergence.manifest_json();
    let audio_convergence_json = audio_convergence.manifest_json();
    let preview_mesh_json = preview_mesh.manifest_json();
    format!(
        concat!(
            "{{\n",
            "  \"schema\": \"frankensim-euler-cinematic-critique-v6\",\n",
            "  \"authority\": \"source-bound analytical simulation visualization; physically informed but uncalibrated synthesis; artistic spatial presentation\",\n",
            "  \"video\": {{\"width\": {width}, \"height\": {height}, \"sequence_frames\": {sequence_frames}, \"rendered_first_frame\": {rendered_first}, \"rendered_frame_count\": {rendered_count}, \"complete_sequence\": {complete}, \"fps\": {fps}, \"duration_s\": {duration:.9}, \"spp\": {spp}, \"sampling\": {sampling}, \"render_seed_salt\": {seed_salt}, \"max_depth\": {max_depth}, \"shutter_angle_degrees\": {shutter_angle}, \"shutter_duration_s\": {shutter_duration:.17e}, \"shutter_convention\": \"back-loaded-frame-boundary\", \"final_shutter_closes_at_cutoff\": true, \"first_shutter_opens_at_s\": {first_shutter_open:.17e}, \"shutter_distribution\": \"stratified-counter-v1\", \"frame_seed_schedule_version\": {seed_version}, \"denoise_requested\": {denoise_requested}, \"raw_exr_profile\": \"{raw_profile}\", \"exposure_ev\": {exposure_ev}, \"raw_sequence_identity\": \"{raw_sequence}\", \"preview_sequence_identity\": \"{preview_sequence}\", \"over_range_linear_channels\": {over_range}, \"gamut_mapped_pixels\": {gamut_mapped}}},\n",
            "  \"timeline\": {{\"video_pts_convention\": \"encoded-frame-start\", \"video_first_pts_s\": 0.00000000000000000e0, \"video_final_pts_s\": {video_final_pts:.17e}, \"video_frame_interval_s\": {frame_interval:.17e}, \"audio_first_sample_time_s\": 0.00000000000000000e0, \"audio_video_clock_origins_aligned\": true, \"shutter_close_convention\": \"mechanics-frame-end\", \"first_shutter_close_time_s\": {first_shutter_close:.17e}, \"final_shutter_close_time_s\": {final_shutter_close:.17e}, \"note\": \"encoded PTS is distinct from the mechanics time at which each back-loaded shutter closes\"}},\n",
            "  \"render_execution\": {{\"policy\": \"deterministic-parked-crew-tile-v1\", \"requested_workers\": {requested_workers}, \"maximum_effective_workers\": {effective_workers}, \"tile_width\": {tile_width}, \"tile_height\": {tile_height}, \"memory_limit_bytes\": {memory_limit}, \"maximum_peak_memory_bytes\": {peak_memory}, \"measured_frames\": {measured_frames}, \"timing_ns\": {{\"setup\": {setup_ns}, \"traversal\": {traversal_ns}, \"tile_compute_sum\": {compute_ns}, \"tile_merge_sum\": {merge_ns}, \"publication\": {publication_ns}, \"idle_worker_capacity\": {idle_ns}}}}},\n",
            "  \"preview_mesh\": {preview_mesh},\n",
            "  \"denoise\": {{\"requested\": {denoise_requested}, \"applied_frames\": {denoised_frames}, \"pipeline\": \"{denoise_pipeline}\", \"authority\": \"biased-display-derivative\", \"maximum_retained_bytes\": {denoise_bytes}, \"maximum_history_frames\": {history_frames}}},\n",
            "  \"mechanics\": {{\"model\": \"Thorne-2026-small-angle-rolling-plus-Bildsten-boundary-layer\", \"source_id\": \"{source_id}\", \"model_authority\": \"{model_authority}\", \"physical_validation\": \"{physical_validation}\", \"specimen\": {{\"diameter_m\": {diameter:.17e}, \"thickness_m\": {thickness:.17e}, \"mass_kg\": {mass:.17e}, \"outer_fillet_radius_m\": {fillet:.17e}}}, \"inputs\": {{\"gravity_m_per_s2\": {gravity:.17e}, \"source_initial_inclination_rad\": {source_initial_theta:.17e}, \"air_density_kg_per_m3\": {air_density:.17e}, \"air_dynamic_viscosity_pa_s\": {air_viscosity:.17e}, \"bildsten_dimensionless_prefactor\": {bildsten_prefactor:.17e}, \"maximum_steps\": {maximum_steps}}}, \"integration\": {{\"coarse_timestep_s\": {coarse_dt:.17e}, \"fine_timestep_s\": {fine_dt:.17e}, \"published_source_timestep_s\": {source_dt:.17e}, \"source_sample_count\": {source_samples}, \"source_duration_s\": {source_duration:.17e}, \"retained_tail_sample_count\": {tail_samples}, \"retained_tail_duration_s\": {tail_duration:.17e}, \"terminal\": \"{terminal:?}\", \"positive_validity_cutoff_rad\": {cutoff:.17e}}}, \"refinement\": {{\"terminal_time_difference_s\": {refine_time:.17e}, \"total_work_difference_j\": {refine_work:.17e}, \"output_consistency\": {output_convergence}, \"claim\": \"single admitted dt/dt2 consistency pair for terminal time, encoded interval impulse, renderer-stratum pose, and body-contact chirp; not experimental validation or an asymptotic-order certificate\"}}, \"channels\": {{\"rolling_coefficient_mu\": {rolling_mu:.17e}, \"rolling_work_j\": {rolling_work:.17e}, \"boundary_layer_work_j\": {gas_work:.17e}}}, \"first_retained_qoi\": {{\"inclination_rad\": {first_theta:.17e}, \"precession_rad_per_s\": {first_precession:.17e}, \"spin_rad_per_s\": {first_spin:.17e}}}, \"last_retained_qoi\": {{\"inclination_rad\": {last_theta:.17e}, \"precession_rad_per_s\": {last_precession:.17e}, \"spin_rad_per_s\": {last_spin:.17e}}}, \"energy\": {{\"initial_j\": {initial_energy:.17e}, \"final_j\": {final_energy:.17e}, \"closure_residual_j\": {energy_residual:.17e}, \"relative_abs_residual\": {relative_residual:.17e}}}, \"trajectory_identity\": \"{trajectory_identity}\"}},\n",
            "  \"audio\": {{\"sample_rate_hz\": {audio_rate}, \"wav_identity\": \"{wav_identity}\", \"authority\": \"physically-informed-uncalibrated\", \"calibrated\": false, \"procedural_texture\": false, \"excitation\": \"kinematically implied interval-mean SI normal reaction applied as unit action/reaction to disc and glass modes\", \"contact_phase\": \"disc modes use body-contact harmonic-one azimuth with rate Omega*cos(theta); glass modes use the independently retained base-frame contact azimuth\", \"chirp_start_hz\": {chirp_start:.17e}, \"chirp_end_hz\": {chirp_end:.17e}, \"assumed_reconstruction_ceiling_hz\": {audio_bandwidth:.17e}, \"modal_parameter_set_identity\": \"{modal_identity}\", \"modal_parameter_set_binding_scope\": \"outer fixture manifest; WAV sound config binds the prepared modal model identity\", \"modal_parameter_set_disclosure\": \"{modal_disclosure}\", \"continuous_crop_policy\": \"{warm_start_policy}\", \"continuous_crop_first_source_audio_frame\": {crop_first_source_frame}, \"continuous_crop_end_source_audio_frame\": {crop_end_source_frame}, \"continuous_crop_source_trajectory_identity\": \"{warm_start_source_identity}\", \"continuous_crop_published_trajectory_identity\": \"{published_trajectory_identity}\", \"continuous_crop_source_sound_configuration_identity\": \"{source_sound_configuration_identity}\", \"continuous_crop_resampler_identity\": \"{crop_resampler_identity}\", \"continuous_crop_binding_identity\": \"{warm_start_checkpoint_identity}\", \"timestep_convergence\": {audio_convergence}, \"pre_master_peak_fs\": {pre_master_peak:.17e}, \"master_gain_db\": {master_gain:.9}, \"initial_fade_sample_frames\": {initial_fade}, \"terminal_fade_sample_frames\": {terminal_fade}, \"terminal_fade_application\": \"exactly once: dry stems for dry output or post-propagation stereo for spatial output\", \"mix_policy\": \"one content-derived digital mastering gain to 0.45 FS; no limiter\", \"spatialization\": {spatialization}}},\n",
            "  \"mux\": {mux},\n",
            "  \"no_claims\": [\"the analytical model reproduces published equations and fitted rolling coefficient but is not a full fluid-structure-contact solve or a raw measured trajectory\", \"the output-space dt/dt2 gate establishes one encoded-model consistency pair for terminal time, interval impulse, renderer-stratum pose, body-contact chirp, resampled modal drive, and unmastered modal stems; it is not an asymptotic convergence-order certificate, contact-law validation, acoustic calibration, psychoacoustic equivalence, or experimental validation\", \"the positive inclination cutoff is horizon censoring, not theta zero, loss of contact, or a resolved terminal impact\", \"the rendered 10 mm glass plate, housing, and support are representative studio apparatus; mechanics uses a static rigid plane and does not simulate the displayed base's compliance, wobble, thickness, mounting, or the paper's 3.2 mm glass-on-butyl rig\", \"the harmonic-one contact shapes, modal frequencies, damping, masses, and radiation gains are representative rather than measured for this specimen and rig\", \"declared excitation completeness means complete only for this authored reduced-channel sonification, not complete physical acoustic forcing\", \"the 256 Hz reconstruction ceiling is a conservative authored assumption, not a certified bandlimit of the interval controls\", \"the centered FIR and modal state both continue through one source-bound video frame (2000 samples) of real prehistory and are then cropped without a new reflection boundary; earlier acoustic history is not reconstructed and the initial fade remains a presentation taper\", \"the waveform, loudness, spectral envelope, terminal chatter, microphone, room, HRTF, and sound-pressure level are not experimentally validated\", \"spatial output clamps propagation tails, so listener audio does not claim to contain the exact source cutoff sample\", \"digital mastering is presentation normalization, not a pascal or SPL prediction\", \"the radial spin fiducial is visualization-only and excluded from specimen geometry, contact, and mass\", \"image quality is final only after native-4K sample-rung review and complete-sequence verification\"]\n",
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
        video_final_pts = f64::from(config.frames - 1) / f64::from(CRITIQUE_FPS),
        frame_interval = 1.0 / f64::from(CRITIQUE_FPS),
        first_shutter_close = 1.0 / f64::from(CRITIQUE_FPS),
        final_shutter_close = f64::from(config.frames) / f64::from(CRITIQUE_FPS),
        spp = config.render_sample_ceiling(),
        sampling = sampling_json,
        seed_salt = config.render_seed_salt,
        max_depth = config.max_depth,
        shutter_angle = config.shutter_angle_degrees,
        shutter_duration =
            f64::from(config.shutter_angle_degrees) / 360.0 / f64::from(CRITIQUE_FPS),
        first_shutter_open = 1.0 / f64::from(CRITIQUE_FPS)
            - f64::from(config.shutter_angle_degrees) / 360.0 / f64::from(CRITIQUE_FPS),
        seed_version = CRITIQUE_FRAME_SEED_SCHEDULE_VERSION,
        denoise_requested = config.denoise_previews,
        raw_profile = raw_profile,
        exposure_ev = CRITIQUE_EXPOSURE_EV,
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
        preview_mesh = preview_mesh_json,
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
        gravity = run.parameters.gravity_m_per_s2,
        source_initial_theta = run.parameters.initial_theta_rad,
        air_density = run
            .provenance
            .bildsten_density_kg_per_m3
            .expect("source-bound run retains air density"),
        air_viscosity = run
            .provenance
            .bildsten_dynamic_viscosity_pa_s
            .expect("source-bound run retains air viscosity"),
        bildsten_prefactor = run
            .provenance
            .bildsten_dimensionless_prefactor
            .expect("source-bound run retains Bildsten prefactor"),
        maximum_steps = run.parameters.maximum_steps,
        coarse_dt = refinement.coarse.parameters.timestep_s,
        fine_dt = refinement.fine.parameters.timestep_s,
        source_dt = run.parameters.timestep_s,
        source_samples = run.samples.len(),
        source_duration = last_source_sample.time_s,
        tail_samples = trajectory_samples.len(),
        tail_duration = last_visual_sample.time_s,
        terminal = run.terminal,
        cutoff = run.parameters.validity_cutoff_theta_rad,
        refine_time = refinement.terminal_time_difference_s,
        refine_work = refinement.total_work_difference_j,
        output_convergence = output_convergence_json,
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
        warm_start_policy = AUDIO_PREROLL_POLICY_ID,
        warm_start_source_identity = audio_warm_start_source_identity.to_hex(),
        warm_start_checkpoint_identity = audio_warm_start_checkpoint_identity.to_hex(),
        source_sound_configuration_identity = audio_source_sound_configuration_identity.to_hex(),
        published_trajectory_identity = audio_published_trajectory_identity.to_hex(),
        crop_resampler_identity = audio_crop_resampler_identity.to_hex(),
        crop_first_source_frame = audio_crop_first_source_frame,
        crop_end_source_frame = audio_crop_end_source_frame,
        pre_master_peak = audio_pre_master_peak_fs,
        master_gain = audio_master_gain_db,
        initial_fade = INITIAL_FADE_SAMPLE_FRAMES,
        terminal_fade = TERMINAL_FADE_SAMPLE_FRAMES,
        spatialization = spatialization_json,
        audio_convergence = audio_convergence_json,
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

    fn stereo_rms(samples: &[StereoSample]) -> f64 {
        assert!(!samples.is_empty());
        let sum_of_squares = samples.iter().fold(0.0, |sum, sample| {
            sample.left_fs.mul_add(
                sample.left_fs,
                sample.right_fs.mul_add(sample.right_fs, sum),
            )
        });
        (sum_of_squares / (2 * samples.len()) as f64).sqrt()
    }

    #[test]
    fn default_is_an_eight_second_practical_preview() {
        let config = CinematicFixtureConfig::default();
        config.validate().unwrap();
        assert_eq!(config.frames, 8 * CRITIQUE_FPS);
        assert_eq!(config.frame_window, CinematicFrameWindow::Full);
        assert_eq!(config.render_seed_salt, 0);
        assert_eq!((config.width, config.height), (320, 180));
        assert_eq!(config.azimuthal_segments, 512);
        assert_eq!(config.arc_subdivisions_per_arc, 64);
        assert_eq!(config.shutter_angle_degrees, 180);
        assert!(config.render_workers > 0);
        assert_eq!(config.tile_width, default_render_tile_edge());
        assert_eq!(config.tile_height, default_render_tile_edge());
        assert!(config.denoise_previews);
        assert!(config.retain_full_aov_exr);
        assert!(config.spatialize_audio);
        assert_eq!(config.adaptive_sampling, None);
        assert_eq!(config.render_sample_ceiling(), config.samples_per_pixel);
    }

    #[test]
    fn render_tessellation_controls_are_bounded_without_changing_production_defaults() {
        let mut config = CinematicFixtureConfig::default();
        config.azimuthal_segments = 8;
        config.arc_subdivisions_per_arc = 1;
        config.validate().unwrap();
        config.azimuthal_segments = MAX_EULER_AZIMUTHAL_SEGMENTS;
        config.arc_subdivisions_per_arc = MAX_EULER_ARC_SUBDIVISIONS;
        config.validate().unwrap();

        config.azimuthal_segments = 7;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(
                "azimuthal_segments must be in 8..=4096"
            ))
        ));
        config.azimuthal_segments = 8;
        config.arc_subdivisions_per_arc = 0;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(
                "arc_subdivisions_per_arc must be in 1..=1024"
            ))
        ));
    }

    #[test]
    fn preview_mesh_manifest_reports_admitted_resolution_counts_and_errors() {
        let manifest = FixturePreviewMeshEvidence {
            azimuthal_segments: 512,
            arc_subdivisions_per_arc: 64,
            vertex_count: 33_280,
            triangle_count: 65_536,
            maximum_meridian_chord_error_m: 1.25e-7,
            maximum_azimuthal_chord_error_m: 2.5e-7,
        }
        .manifest_json();
        assert!(manifest.contains("\"authority\":\"render-only chordal approximation\""));
        assert!(manifest.contains("\"azimuthal_segments\":512"));
        assert!(manifest.contains("\"arc_subdivisions_per_arc\":64"));
        assert!(manifest.contains("\"vertex_count\":33280"));
        assert!(manifest.contains("\"triangle_count\":65536"));
        assert!(manifest.contains(&format!(
            "\"maximum_meridian_chord_error_m\":{:.17e}",
            1.25e-7
        )));
        assert!(manifest.contains(&format!(
            "\"maximum_azimuthal_chord_error_m\":{:.17e}",
            2.5e-7
        )));
    }

    fn adaptive_fixture_config() -> CinematicFixtureConfig {
        let mut config = CinematicFixtureConfig::default();
        config.denoise_previews = false;
        config.retain_full_aov_exr = false;
        config.adaptive_sampling = Some(CinematicAdaptiveSamplingConfig {
            minimum_samples_per_pixel: 8,
            maximum_samples_per_pixel: 64,
            decision_batch_samples: 4,
            absolute_error: 1.0e-4,
            relative_error: 0.02,
            dark_floor: 1.0e-5,
        });
        config
    }

    #[test]
    fn adaptive_policy_is_explicit_bounded_and_canonical() {
        let config = adaptive_fixture_config();
        config.validate().unwrap();
        assert_eq!(config.render_sample_ceiling(), 64);
        let policy = config.adaptive_policy().unwrap().unwrap();
        assert_eq!(policy.minimum_samples(), 8);
        assert_eq!(policy.batch_samples(), 4);
        assert_eq!(policy.absolute_error().to_bits(), 1.0e-4_f64.to_bits());
        assert_eq!(policy.relative_error().to_bits(), 0.02_f64.to_bits());
        assert_eq!(policy.dark_floor().to_bits(), 1.0e-5_f64.to_bits());

        let mut negative_zero = config;
        let adaptive = negative_zero.adaptive_sampling.as_mut().unwrap();
        adaptive.absolute_error = -0.0;
        adaptive.relative_error = -0.0;
        adaptive.dark_floor = -0.0;
        negative_zero.validate().unwrap();
        let canonical = negative_zero.adaptive_policy().unwrap().unwrap();
        assert_eq!(canonical.absolute_error().to_bits(), 0.0_f64.to_bits());
        assert_eq!(canonical.relative_error().to_bits(), 0.0_f64.to_bits());
        assert_eq!(canonical.dark_floor().to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn adaptive_policy_accepts_final_diagnostic_aovs_and_count_aware_denoising() {
        let mut final_diagnostic = adaptive_fixture_config();
        final_diagnostic.retain_full_aov_exr = true;
        final_diagnostic.validate().unwrap();

        let mut denoised = adaptive_fixture_config();
        denoised.denoise_previews = true;
        denoised.validate().unwrap();

        let mut final_diagnostic_and_denoised = adaptive_fixture_config();
        final_diagnostic_and_denoised.retain_full_aov_exr = true;
        final_diagnostic_and_denoised.denoise_previews = true;
        final_diagnostic_and_denoised.validate().unwrap();
    }

    #[test]
    fn adaptive_policy_refuses_invalid_controls() {
        let mut below_minimum = adaptive_fixture_config();
        below_minimum
            .adaptive_sampling
            .as_mut()
            .unwrap()
            .maximum_samples_per_pixel = 7;
        assert!(below_minimum.validate().is_err());

        let mut non_finite = adaptive_fixture_config();
        non_finite
            .adaptive_sampling
            .as_mut()
            .unwrap()
            .relative_error = f64::NAN;
        assert!(non_finite.validate().is_err());
    }

    #[test]
    fn adaptive_raw_artifact_exports_exact_counts_and_manifest_evidence() {
        use fs_render::{
            lighting::EnvironmentMap,
            tracer::{
                Camera, DirectStrategy, RenderExecutionConfig, Sampler, Scene, Settings,
                render_adaptive_with_execution,
            },
        };

        let scene = Scene {
            primitives: Vec::new(),
            lights: Vec::new(),
            environment: Some(
                EnvironmentMap::try_from_linear_srgb(4, 2, vec![[0.25, 0.5, 0.75]; 8], 1.0)
                    .unwrap(),
            ),
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 0.0),
                forward: GeomVec3::new(1.0, 0.0, 0.0),
                up: GeomVec3::new(0.0, 1.0, 0.0),
                half_tan: 0.5,
            },
        };
        let adaptive = CinematicAdaptiveSamplingConfig {
            minimum_samples_per_pixel: 2,
            maximum_samples_per_pixel: 4,
            decision_batch_samples: 1,
            absolute_error: 1.0e30,
            relative_error: 0.0,
            dark_floor: 0.0,
        };
        let settings = Settings {
            width: 2,
            height: 1,
            spp: adaptive.maximum_samples_per_pixel,
            max_depth: 1,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed: 0x4144_4150_5449_5645,
        };
        let execution =
            RenderExecutionConfig::try_new(1, 1, 1, 64 * 1024 * 1024, RunId(0x4144_4150_5449_5645))
                .unwrap();
        let output = with_test_cx(|cx| {
            render_adaptive_with_execution(
                &scene,
                cx,
                &settings,
                adaptive.policy().unwrap(),
                &execution,
            )
            .unwrap()
        });
        assert_eq!(output.film.sample_counts(), [2, 2]);

        let mut evidence = FixtureAdaptiveSamplingEvidence::default();
        let count_identity = evidence.observe(17, &output.film, adaptive).unwrap();
        let mut config = adaptive_fixture_config();
        config.adaptive_sampling = Some(adaptive);
        evidence.validate_complete(&config, &(17..18)).unwrap();

        let provenance = CinematicAovProvenance::try_new(
            17,
            17.0 / 24.0,
            16.0 / 24.0,
            18.0 / 24.0,
            hash_domain("adaptive-test-trajectory", b"trajectory"),
            hash_domain("adaptive-test-scene", b"scene"),
            hash_domain("adaptive-test-composition", b"composition"),
        )
        .unwrap();
        let rgb = output.film.to_linear_srgb();
        let exr = adaptive_beauty_to_exr(
            &output.film,
            rgb,
            adaptive,
            settings,
            provenance,
            count_identity,
        )
        .unwrap();
        let decoded = fs_img::read_exr(&exr).unwrap();
        let samples = decoded
            .channels
            .iter()
            .find(|channel| channel.name == "samples")
            .unwrap();
        assert_eq!(samples.data, [2.0, 2.0]);
        assert!(decoded.attributes.iter().any(|attribute| {
            attribute.name == "frankensim.render.sampleCounts"
                && attribute.value == count_identity.to_hex().as_bytes()
        }));

        let manifest = evidence.manifest_json(adaptive);
        assert!(manifest.contains("\"mode\":\"adaptive-raw-xyz-dispersion-v1\""));
        assert!(manifest.contains("\"actual_minimum_spp\":2"));
        assert!(manifest.contains("\"actual_maximum_spp\":2"));
        assert!(manifest.contains("\"total_paths\":4"));
        assert!(manifest.contains(&count_identity.to_hex()));

        let mut mismatched_policy = adaptive;
        mismatched_policy.maximum_samples_per_pixel = 3;
        assert!(
            FixtureAdaptiveSamplingEvidence::default()
                .observe(17, &output.film, mismatched_policy)
                .is_err()
        );
    }

    #[test]
    fn impulse_queries_only_canonicalize_terminal_roundoff() {
        let nominal_terminal_s = 8.0_f64;
        let retained_terminal_s = f64::from_bits(nominal_terminal_s.to_bits() - 1);
        let tolerance_s = 32.0 * f64::EPSILON * nominal_terminal_s.max(1.0);
        assert_eq!(
            canonicalize_terminal_query_time_s(
                nominal_terminal_s,
                retained_terminal_s,
                tolerance_s,
            )
            .to_bits(),
            retained_terminal_s.to_bits()
        );

        let interior_query_s = retained_terminal_s - 0.25;
        assert_eq!(
            canonicalize_terminal_query_time_s(interior_query_s, retained_terminal_s, tolerance_s,)
                .to_bits(),
            interior_query_s.to_bits()
        );
        let out_of_horizon_query_s = retained_terminal_s + 2.0 * tolerance_s;
        assert_eq!(
            canonicalize_terminal_query_time_s(
                out_of_horizon_query_s,
                retained_terminal_s,
                tolerance_s,
            )
            .to_bits(),
            out_of_horizon_query_s.to_bits()
        );
    }

    #[test]
    fn source_bound_fixture_requires_exact_eight_second_frame_count() {
        let mut config = CinematicFixtureConfig::default();
        config.frames = 191;
        assert!(matches!(
            config.validate(),
            Err(CinematicFixtureError::InvalidConfig(_))
        ));
        config.frames = 193;
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
    fn encoded_pts_and_back_loaded_shutter_use_distinct_times_in_one_frame_interval() {
        let trajectory_start_s = 0.0;
        let trajectory_end_s = 8.0;
        let exposure_s = 0.5 / f64::from(CRITIQUE_FPS);
        for frame in 0..CRITIQUE_FRAMES {
            let times =
                frame_timeline_times(frame, CRITIQUE_FRAMES, trajectory_start_s, trajectory_end_s);
            let expected_pts_s = f64::from(frame) / f64::from(CRITIQUE_FPS);
            let expected_shutter_close_s = (f64::from(frame) + 1.0) / f64::from(CRITIQUE_FPS);
            let shutter_open_s = times.shutter_close_time_s - exposure_s;
            assert_eq!(times.presentation_time_s, expected_pts_s);
            assert_eq!(times.shutter_close_time_s, expected_shutter_close_s);
            assert!(trajectory_start_s <= times.previous_presentation_time_s);
            assert!(times.previous_presentation_time_s <= times.presentation_time_s);
            assert!(times.presentation_time_s <= shutter_open_s);
            assert!(shutter_open_s <= times.shutter_close_time_s);
            assert!(times.shutter_close_time_s <= times.next_presentation_time_s);
            assert!(times.next_presentation_time_s <= trajectory_end_s);
            if frame == 0 {
                assert_eq!(times.presentation_time_s, 0.0);
                assert_eq!(times.previous_presentation_time_s, 0.0);
                assert_eq!(shutter_open_s, 1.0 / 48.0);
                assert_eq!(times.shutter_close_time_s, 1.0 / 24.0);
                assert_eq!(times.next_presentation_time_s, 1.0 / 24.0);
            }
            if frame + 1 == CRITIQUE_FRAMES {
                assert_eq!(times.presentation_time_s, 191.0 / 24.0);
                assert_eq!(times.shutter_close_time_s, trajectory_end_s);
                assert_eq!(times.next_presentation_time_s, trajectory_end_s);
            }
        }
    }

    #[test]
    fn audio_and_encoded_video_master_clocks_share_the_zero_origin() {
        let (audio_frames, video_clock, audio_clock) =
            fixture_master_clocks(CRITIQUE_FRAMES).unwrap();
        assert_eq!(video_clock.start_tick(), 0);
        assert_eq!(audio_clock.start_tick(), 0);
        assert_eq!(video_clock.end_tick_exclusive(), i64::from(CRITIQUE_FRAMES));
        assert_eq!(audio_frames, u64::from(CRITIQUE_FRAMES) * 2_000);
        assert_eq!(
            audio_clock.end_tick_exclusive(),
            i64::try_from(audio_frames).unwrap()
        );
    }

    fn accumulated_metrics(pairs: &[(f64, f64)], floor: f64) -> RelativeErrorMetrics {
        let mut accumulator = SymmetricErrorAccumulator::default();
        for &(coarse, fine) in pairs {
            accumulator.observe(coarse, fine).unwrap();
        }
        accumulator.metrics(floor).unwrap()
    }

    #[test]
    fn audio_convergence_identical_zero_signals_pass_fixed_gates() {
        let identical = accumulated_metrics(&[(0.0, 0.0), (0.0, 0.0)], 1.0e-12);
        assert_eq!(identical.nrmse.to_bits(), 0.0_f64.to_bits());
        assert_eq!(identical.normalized_peak.to_bits(), 0.0_f64.to_bits());
        enforce_audio_convergence_thresholds(identical, identical, identical).unwrap();
    }

    #[test]
    fn audio_convergence_rejects_equal_integral_different_drive_envelopes() {
        let different_envelope = accumulated_metrics(&[(2.0, 0.0), (0.0, 2.0)], 1.0e-12);
        let identical = accumulated_metrics(&[(1.0, 1.0)], 1.0e-12);
        assert!(different_envelope.nrmse > AUDIO_CONVERGENCE_DRIVE_NRMSE_LIMIT);
        assert!(
            enforce_audio_convergence_thresholds(different_envelope, identical, identical,)
                .is_err()
        );
    }

    #[test]
    fn audio_convergence_rejects_raw_stem_amplitude_error_before_mastering() {
        let amplitude_error = accumulated_metrics(&[(2.0, 1.0), (2.0, 1.0)], 1.0e-12);
        let identical = accumulated_metrics(&[(1.0, 1.0)], 1.0e-12);
        assert!(amplitude_error.nrmse > AUDIO_CONVERGENCE_STEM_NRMSE_LIMIT);
        assert!(
            enforce_audio_convergence_thresholds(identical, identical, amplitude_error).is_err()
        );
    }

    #[test]
    fn audio_convergence_rejects_checkpoint_range_mismatch() {
        assert!(
            validate_audio_pair_progress(
                "test", 2_000, 4_000, 2_000, 4_000, 4_000, 2_000, 4_001, 4_001
            )
            .is_err()
        );
    }

    #[test]
    fn audio_convergence_metrics_are_bit_exact_on_replay() {
        let pairs = [(0.25, 0.20), (-0.5, -0.45), (0.75, 0.70)];
        let first = accumulated_metrics(&pairs, 1.0e-12);
        let second = accumulated_metrics(&pairs, 1.0e-12);
        assert_eq!(first.nrmse.to_bits(), second.nrmse.to_bits());
        assert_eq!(
            first.normalized_peak.to_bits(),
            second.normalized_peak.to_bits()
        );
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

    #[test]
    fn source_bound_audio_is_audible_and_admitted_before_rendering() {
        with_test_cx(|cx| {
            let benchmark = Thorne2026SteelGlassBenchmark::ambient().unwrap();
            let profile = benchmark.resolve_specimen(cx).unwrap();
            let refinement = cinematic_thorne_2026_refinement_evidence(&benchmark).unwrap();
            let preroll_duration_s =
                f64::from(AUDIO_PREROLL_VIDEO_FRAMES) / f64::from(CRITIQUE_FPS);
            let published_duration_s = f64::from(CRITIQUE_FRAMES) / f64::from(CRITIQUE_FPS);
            let (coarse_preroll_trajectory, coarse_trajectory) =
                RenderTrajectory::from_reduced_decay_run_with_causal_preroll(
                    &refinement.coarse,
                    &profile,
                    preroll_duration_s,
                    published_duration_s,
                    cx,
                )
                .unwrap();
            let (preroll_trajectory, trajectory) =
                RenderTrajectory::from_reduced_decay_run_with_causal_preroll(
                    &refinement.fine,
                    &profile,
                    preroll_duration_s,
                    published_duration_s,
                    cx,
                )
                .unwrap();
            let source_duration_s = preroll_duration_s + published_duration_s;
            for source in [&coarse_preroll_trajectory, &preroll_trajectory] {
                assert_eq!(
                    source
                        .samples()
                        .last()
                        .expect("source trajectory endpoint")
                        .input()
                        .time_s
                        .to_bits(),
                    source_duration_s.to_bits()
                );
            }
            for published in [&coarse_trajectory, &trajectory] {
                assert_eq!(
                    published
                        .samples()
                        .last()
                        .expect("published trajectory endpoint")
                        .input()
                        .time_s
                        .to_bits(),
                    published_duration_s.to_bits()
                );
                let motion = EulerRenderMotionBridge::new(published);
                motion
                    .sample_at_time(published_duration_s, EventEvaluationSide::RightLimit)
                    .expect("the exact exclusive film boundary is an admitted shutter endpoint");
                assert!(
                    motion
                        .sample_at_time(
                            published_duration_s + 1.0e-9,
                            EventEvaluationSide::RightLimit,
                        )
                        .is_err(),
                    "the strict motion timeline must still refuse real extrapolation"
                );
            }
            let mut convergence_config = CinematicFixtureConfig::default();
            convergence_config.samples_per_pixel = 4;
            let convergence = fixture_output_convergence_evidence(
                &coarse_preroll_trajectory,
                &coarse_trajectory,
                &preroll_trajectory,
                &trajectory,
                &convergence_config,
                &refinement,
                preroll_duration_s,
                published_duration_s,
                cx,
            )
            .unwrap();
            eprintln!("{}", convergence.diagnostics());
            assert_eq!(
                convergence.shutter_exact_renderer_sample_count,
                CRITIQUE_FRAMES as usize * 4
            );
            assert_eq!(
                convergence.shutter_stratum_boundary_sample_count,
                CRITIQUE_FRAMES as usize * 5
            );
            assert!(convergence.shutter_interpolation_knot_sample_count > CRITIQUE_FRAMES as usize);
            assert_eq!(
                convergence.shutter_interpolation_midpoint_sample_count,
                convergence.shutter_interpolation_knot_sample_count + CRITIQUE_FRAMES as usize
            );
            assert_eq!(
                convergence.shutter_pose_sample_count,
                convergence.shutter_exact_renderer_sample_count
                    + convergence.shutter_stratum_boundary_sample_count
                    + convergence.shutter_interpolation_knot_sample_count
                    + convergence.shutter_interpolation_midpoint_sample_count
            );
            assert!(
                convergence.terminal_time_difference_s
                    <= OUTPUT_CONVERGENCE_TERMINAL_TIME_ABSOLUTE_LIMIT_S
            );
            assert!(
                convergence.terminal_time_relative_difference
                    <= OUTPUT_CONVERGENCE_TERMINAL_TIME_RELATIVE_LIMIT
            );
            assert!(
                convergence.maximum_impulse_identity_residual_n_s
                    <= convergence.maximum_impulse_identity_tolerance_n_s
            );
            assert!(
                convergence.maximum_cumulative_normal_impulse_relative_difference
                    <= OUTPUT_CONVERGENCE_RELATIVE_IMPULSE_LIMIT
            );
            assert!(
                convergence.maximum_center_of_mass_difference_m <= OUTPUT_CONVERGENCE_COM_LIMIT_M
            );
            assert!(
                convergence.maximum_orientation_difference_rad
                    <= OUTPUT_CONVERGENCE_ORIENTATION_LIMIT_RAD
            );
            assert!(
                convergence.maximum_chirp_difference_hz
                    <= OUTPUT_CONVERGENCE_CHIRP_ABSOLUTE_LIMIT_HZ
            );
            assert!(
                convergence.maximum_relative_chirp_difference
                    <= OUTPUT_CONVERGENCE_CHIRP_RELATIVE_LIMIT
            );
            let mut non_finite = convergence;
            non_finite.maximum_center_of_mass_difference_m = f64::NAN;
            assert!(admit_fixture_output_convergence(non_finite).is_err());
            drop(coarse_trajectory);
            let source_crop_boundary = preroll_trajectory
                .samples()
                .iter()
                .find(|sample| sample.input().time_s.to_bits() == preroll_duration_s.to_bits())
                .expect("causal source retains the exact picture/audio crop boundary");
            let published_initial = &trajectory.samples()[0];
            assert_eq!(source_crop_boundary.state(), published_initial.state());
            assert_eq!(
                source_crop_boundary.input().qois,
                published_initial.input().qois
            );
            assert_eq!(
                published_initial.input().time_s.to_bits(),
                0.0_f64.to_bits()
            );
            assert_eq!(
                published_initial.input().interval_start_time_s.to_bits(),
                0.0_f64.to_bits()
            );
            let initial_channels = published_initial.input().channels;
            let initial_work_j = initial_channels.gravity.work_j
                + initial_channels.contact.work_j
                + initial_channels.rolling.work_j
                + initial_channels.base.work_j
                + initial_channels.gas.work_j;
            assert_eq!(initial_work_j.to_bits(), 0.0_f64.to_bits());
            let artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
                hash_domain(
                    "org.frankensim.euler-critique.audio-regression.v1",
                    b"source-bound-eight-second-tail",
                ),
                trajectory,
                Vec::new(),
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            )
            .unwrap();
            let preroll_artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
                hash_domain(
                    "org.frankensim.euler-critique.audio-regression-preroll.v1",
                    artifact.receipt().artifact_identity().as_bytes(),
                ),
                preroll_trajectory,
                Vec::new(),
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            )
            .unwrap();
            let coarse_preroll_artifact = EulerRenderTrajectoryArtifact::try_from_trajectory(
                hash_domain(
                    "org.frankensim.euler-critique.audio-regression-coarse-preroll.v1",
                    coarse_preroll_trajectory
                        .metadata()
                        .configuration_identity
                        .as_bytes(),
                ),
                coarse_preroll_trajectory,
                Vec::new(),
                RenderTrajectoryCodecBudget::DEFAULT,
                cx,
            )
            .unwrap();

            let mut expected_warm_start = None;
            for spatialize_audio in [true, false] {
                let mut config = CinematicFixtureConfig::default();
                config.spatialize_audio = spatialize_audio;
                let audio = build_audio(
                    &artifact,
                    &coarse_preroll_artifact,
                    &preroll_artifact,
                    refinement.fine.parameters.gravity_m_per_s2,
                    &config,
                    cx,
                )
                .unwrap();
                assert!(audio.pre_master_peak_fs.is_finite());
                assert!(audio.pre_master_peak_fs > 0.0);
                assert_eq!(
                    audio.warm_start_source_identity,
                    preroll_artifact.receipt().artifact_identity()
                );
                assert_eq!(
                    audio.published_trajectory_identity,
                    artifact.receipt().artifact_identity()
                );
                assert_ne!(
                    audio.warm_start_source_identity, audio.published_trajectory_identity,
                    "full source and rebased picture trajectory must remain distinct"
                );
                assert_eq!(
                    audio.crop_first_source_audio_frame,
                    AUDIO_PREROLL_SAMPLE_FRAMES
                );
                assert_eq!(
                    audio.crop_end_source_audio_frame,
                    u64::from(CRITIQUE_FRAMES + AUDIO_PREROLL_VIDEO_FRAMES) * 2_000
                );
                let synthesis = audio.artifact.manifest().synthesis();
                assert_eq!(
                    synthesis.trajectory_identity, audio.warm_start_source_identity,
                    "WAV provenance must name the full source trajectory that generated the FIR/modal samples"
                );
                assert_ne!(
                    synthesis.configuration_identity, audio.source_sound_configuration_identity,
                    "published crop configuration must be distinct from the full-horizon source configuration"
                );
                assert_ne!(
                    audio.crop_resampler_identity,
                    audio.source_sound_configuration_identity
                );
                let warm_start = (
                    audio.warm_start_source_identity,
                    audio.warm_start_checkpoint_identity,
                );
                if let Some(expected) = expected_warm_start {
                    assert_eq!(warm_start, expected);
                } else {
                    expected_warm_start = Some(warm_start);
                }
                // Content-derived presentation normalization remains
                // explicitly non-SPL. The source itself is now an SI contact
                // reaction, so no arbitrary lower bound on mastering gain is
                // physically meaningful.
                assert!(audio.master_gain_db.is_finite());
                assert!(audio.master_gain_db <= MAX_AUDIO_MASTER_GAIN_DB);
                let meters = audio.artifact.manifest().meters();
                let peak = meters.sample_peak_fs.max(meters.true_peak_estimate_fs);
                assert!(peak > 0.40 && peak <= 0.46, "stored peak {peak}");
                assert_eq!(
                    audio.artifact.manifest().wav().sample_frame_count(),
                    u64::from(CRITIQUE_FRAMES) * 2_000
                );
                let decoded = crate::decode_stereo_wav(
                    audio.artifact.wav_bytes(),
                    AudioArtifactBudget::DEFAULT,
                    cx,
                )
                .unwrap();
                assert_eq!(decoded.samples.len(), CRITIQUE_FRAMES as usize * 2_000);
                assert!(
                    decoded.samples.iter().all(|sample| {
                        sample.left_fs.is_finite() && sample.right_fs.is_finite()
                    })
                );
                let startup_rms = stereo_rms(&decoded.samples[..480]);
                let post_fade_rms = stereo_rms(&decoded.samples[960..4_800]);
                let startup_ratio = startup_rms / post_fade_rms;
                eprintln!(
                    "spatialize_audio={spatialize_audio} modal_warm_start_frames={AUDIO_PREROLL_SAMPLE_FRAMES} startup_0_10ms_rms_fs={startup_rms:.9e} post_fade_20_100ms_rms_fs={post_fade_rms:.9e} startup_ratio={startup_ratio:.9}"
                );
                assert!(startup_rms > 0.0);
                assert!(post_fade_rms > 0.0);
                assert!(
                    startup_rms <= post_fade_rms,
                    "warm-started onset overshot the 20-100 ms reference: {startup_rms:.9e} > {post_fade_rms:.9e} FS"
                );
                let terminal = decoded.samples.last().unwrap();
                assert_eq!(terminal.left_fs, 0.0);
                assert_eq!(terminal.right_fs, 0.0);
                audio
                    .artifact
                    .verify(AudioArtifactBudget::DEFAULT, cx)
                    .unwrap();
            }
        });
    }
}
