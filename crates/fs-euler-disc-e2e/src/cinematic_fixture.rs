//! Direct, deterministic Euler-disc critique-clip producer.
//!
//! This is a product-facing preview path: it runs the reduced mechanics,
//! renders the resulting trajectory, derives dry sound from mechanics control
//! channels, and optionally muxes the image sequence and WAV with `ffmpeg`.
//! It is deliberately labeled uncalibrated and timestep-unconverged.

use core::{fmt, num::NonZeroUsize};
use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    process::Command,
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
use fs_exec::Cx;
use fs_geom::{Point3, Vec3 as GeomVec3};
use fs_img::{
    CinematicColorConfig, CinematicColorLimits, PngColor, PreviewDither,
    transform_cinematic_preview, write_png16,
};
use fs_math::det;
use fs_mbd::{UnitQuaternion, Vec3};
use fs_render::{
    aov::{CinematicAovConfig, CinematicAovLimits, CinematicAovProfile, CinematicAovProvenance},
    camera::{AnimatedCamera, Aperture, CameraProjection, CutSide, PhysicalCamera},
    conductor::{ConductorOptics, ConductorSurface},
    motion::{ShutterConvention, ShutterDistribution},
    tracer::render_cinematic_with_aovs,
};
use fs_rep_frep::SquatDiscEdgeTreatment;

use crate::{
    AUDIO_EXCITATION_ALGORITHM_VERSION, AUDIO_RECONSTRUCTION_FILTER_VERSION,
    AUDIO_RESAMPLING_ALGORITHM_VERSION, AudioArtifactBudget, AudioDryMixSpec,
    AudioEventFractionalDelay, AudioExcitationBudget, AudioExcitationMapper,
    AudioExcitationModelInput, AudioExcitationReduction, AudioMasterSource,
    AudioReconstructionFilterSpec, AudioResampler, AudioResamplingBoundaryPolicy,
    AudioResamplingBudget, AudioResamplingModelInput, ContactParticipationPolicy,
    EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerControlStream, EulerRenderTrajectoryArtifact,
    ModalSynthesisBudget, ModalSynthesisModel, ModalSynthesisModelInput, RenderBaseFrame,
    RenderBaseModeState, RenderChannelAvailability, RenderMassProperties, RenderTrajectory,
    RenderTrajectoryAuthority, RenderTrajectoryCodecBudget, RenderTrajectoryMetadata,
    RenderUnitSystem, RenderWorldFrame, RepresentativeDiscMaterial, SoundWavArtifact, StemGainPan,
    WavMetadata, WavSampleEncoding,
    coupled_runner::{
        CoupledChannelFactors, CoupledControls, CoupledInitialState, CoupledTerminal,
        run_closed_profile_reduced,
    },
    measure_audio, mix_dry_modal_stems,
    render_scene_bridge::{
        EulerCinematicScene, EulerFrameRequest, EulerMaterialStyle, EulerSceneConfig,
        EulerTessellationConfig, euler_scene_smoke_settings,
    },
    representative_modal_preset,
    specimen::DiscProfileSpec,
    timeline_resampling::ExposureEventPolicy,
};

/// Fixed master frame rate used by the cinematic sound contract.
pub const CRITIQUE_FPS: u32 = 24;
/// Minimum admitted cinematic duration: 192 frames = 8 seconds.
pub const CRITIQUE_FRAMES: u32 = 192;
/// Five-millisecond deterministic taper applied at the censored soundtrack end.
const TERMINAL_FADE_SAMPLE_FRAMES: u32 = 240;

/// Bounded settings for one watchable critique artifact.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicFixtureConfig {
    /// Display and raw-master width in pixels.
    pub width: u32,
    /// Display and raw-master height in pixels.
    pub height: u32,
    /// Exact 24 Hz frame count; the sound contract admits 192 through 288.
    pub frames: u32,
    /// Uniform path-tracing samples per pixel.
    pub samples_per_pixel: u32,
    /// Maximum path depth, including dielectric traversal.
    pub max_depth: u32,
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
            samples_per_pixel: 1,
            max_depth: 6,
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
        if !(192..=288).contains(&self.frames) {
            return Err(CinematicFixtureError::InvalidConfig(
                "frames must be in the cinematic sound-contract range 192..=288",
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
        if self.mux_with_ffmpeg && (self.width % 2 != 0 || self.height % 2 != 0) {
            return Err(CinematicFixtureError::InvalidConfig(
                "muxed 4:2:0 video requires even width and height",
            ));
        }
        Ok(())
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
    // The 50 us preview rung halves the eight-second cumulative energy defect
    // seen at 100 us for this configuration (about 2.23% versus 4.52%). It is
    // still explicitly unconverged; the manifest publishes the actual defect.
    let timestep_s = 5.0e-5;
    let maximum_steps = (duration_s / timestep_s).round() as u32;

    progress("stage=mechanics begin");
    let profile = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    }
    .resolve(19_250.0, cx)
    .map_err(pipeline)?;
    let channels = CoupledChannelFactors {
        gravity_m_per_s2: 9.806_65,
        sliding_friction_coefficient: 0.42,
        rolling_resistance_m: 4.0e-5,
        contact_stiffness_n_per_m: 8.0e4,
        contact_damping_n_s_per_m: 3.0,
        base_effective_mass_kg: 0.25,
        base_stiffness_n_per_m: 4.0e4,
        base_damping_n_s_per_m: 4.0,
        gas_rotational_damping_n_m_s: 2.0e-7,
        gas_translation_damping_n_s_per_m: 4.0e-4,
    };
    let controls = CoupledControls {
        timestep_s,
        maximum_steps,
        terminal_inclination_rad: 1.0e-6,
        reimpact_limit: 128,
    };
    // Use an opposite-sign, gravity-scale initial twist. The former same-sign
    // `(precession=16, spin=120)` twist forced the no-slip initializer to give
    // the centre of mass roughly 5 m/s of lateral speed, so the simulated disc
    // left the 180 mm plate. The 36 rad/s candidate is the rounded scale
    // `sqrt(4g/(R*theta))`; it was not fitted to a desired trajectory. The
    // present 6 mm squat profile is outside the thin-disc oracle's admitted
    // geometry, and these rates are only initial conditions: they neither
    // constrain the subsequent coupled motion nor upgrade its authority.
    let inclination_rad = 0.8;
    let precession_rad_per_s = 36.0;
    let initial = CoupledInitialState {
        inclination_rad,
        precession_rad_per_s,
        spin_rad_per_s: -precession_rad_per_s * det::cos(inclination_rad),
    };
    let run = run_closed_profile_reduced(&profile, channels, controls, initial, None, cx)
        .map_err(pipeline)?;
    if !matches!(run.terminal, CoupledTerminal::HorizonReached) {
        return Err(CinematicFixtureError::Pipeline(format!(
            "eight-second preview requires a complete source horizon; mechanics ended as {:?} at {:.9}s",
            run.terminal, run.checkpoint.time_s
        )));
    }
    progress("stage=mechanics complete");

    let profile_identities = profile.content_identities();
    let configuration_identity = mechanics_configuration_identity(
        config,
        profile_identities.profile,
        &channels,
        &controls,
        &initial,
    );
    let trajectory = RenderTrajectory::from_coupled_run(
        RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: profile_identities.profile,
            specimen_chart_identity: profile_identities.chart,
            mass_properties: RenderMassProperties {
                identity: profile_identities.mass_properties,
                properties: run.mass_properties,
            },
            initial_state: run.configuration_initial_state,
            initial_base_mode: RenderBaseModeState {
                displacement_m: run.configuration_initial_base_deflection_m,
                velocity_m_per_s: run.configuration_initial_base_velocity_m_per_s,
            },
            base_model_identity: hash_domain(
                "org.frankensim.euler-critique.base-model.v1",
                b"one-mode-kelvin-voigt-glass-base",
            ),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: hash_domain(
                "org.frankensim.euler-critique.mechanics-model.v1",
                b"closed-profile-reduced-coupled-runner",
            ),
            channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
            configuration_identity,
            configuration_fingerprint: run.checkpoint.configuration_fingerprint,
            timestep_s,
            producer_version: "euler-cinematic-fixture-v1".to_owned(),
            applicability: "deterministic eight-second visual and auditory critique preview"
                .to_owned(),
            no_claims: vec![
                "preview timestep has not been convergence-qualified".to_owned(),
                "mechanics and acoustic parameters have not been calibrated to experiment"
                    .to_owned(),
                "radial spin fiducial is visualization-only and excluded from specimen, contact, and mass mechanics"
                    .to_owned(),
            ],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        },
        &run,
    )
    .map_err(pipeline)?;
    let campaign_identity = hash_domain(
        "org.frankensim.euler-critique.campaign.v1",
        configuration_identity.as_bytes(),
    );
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
        optics: ConductorOptics::representative_tungsten(),
        surface: ConductorSurface::try_rough(0.12).map_err(pipeline)?,
    };
    let scene = EulerCinematicScene::try_build(&trajectory_artifact, &profile, scene_config, cx)
        .map_err(pipeline)?;
    let mut render_settings = euler_scene_smoke_settings(config.width, config.height);
    render_settings.spp = config.samples_per_pixel;
    render_settings.max_depth = config.max_depth;
    let composition_identity = composition_identity(config, scene.scene_identity());
    let mut raw_sequence = DomainHasher::new("org.frankensim.euler-critique.raw-sequence.v1");
    let mut preview_sequence =
        DomainHasher::new("org.frankensim.euler-critique.preview-sequence.v1");
    let mut over_range_channels = 0_u64;
    let mut gamut_mapped_pixels = 0_u64;
    for frame in 0..config.frames {
        if frame % CRITIQUE_FPS == 0 {
            progress(&format!("stage=render frame={frame}/{}", config.frames));
        }
        let frame_time_s = (f64::from(frame) + 0.5) / f64::from(CRITIQUE_FPS);
        let previous_time_s = if frame == 0 {
            frame_time_s
        } else {
            frame_time_s - 1.0 / f64::from(CRITIQUE_FPS)
        };
        let next_time_s = if frame + 1 == config.frames {
            frame_time_s
        } else {
            frame_time_s + 1.0 / f64::from(CRITIQUE_FPS)
        };
        let prepared = scene
            .prepare_frame(EulerFrameRequest {
                frame_time_s,
                exposure_duration_s: 0.0,
                convention: ShutterConvention::Centered,
                distribution: ShutterDistribution::UniformCounterV1,
                event_policy: ExposureEventPolicy::Refuse,
                cut_side: CutSide::After,
            })
            .map_err(pipeline)?;
        if prepared.segments().len() != 1 {
            return Err(CinematicFixtureError::Pipeline(format!(
                "zero-width frame {frame} unexpectedly resolved to {} segments",
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
        let film = render_cinematic_with_aovs(
            scene.scene(),
            scene.camera(),
            prepared.cut_side(),
            cx,
            &render_settings,
            prepared.segments()[0].shutter(),
            CinematicAovConfig::new(
                CinematicAovProfile::DailyCore,
                provenance,
                CinematicAovLimits::default(),
            ),
        )
        .map_err(pipeline)?;
        let exr = film.to_exr().map_err(pipeline)?;
        raw_sequence.update(hash_domain("frame", &exr).as_bytes());
        write_new(&raw_directory.join(format!("frame-{frame:06}.exr")), &exr)?;

        let [red, green, blue] = film.beauty().to_linear_srgb();
        let mut color = CinematicColorConfig::reference_srgb_16();
        color.exposure_ev = -2;
        color.dither = PreviewDither::Disabled;
        let preview = transform_cinematic_preview(
            config.width,
            config.height,
            [&red, &green, &blue],
            color,
            CinematicColorLimits::reference_4k(),
        )
        .map_err(pipeline)?;
        over_range_channels += preview.metadata().over_range_linear_channels();
        gamut_mapped_pixels += preview.metadata().gamut_mapped_pixels();
        let samples = preview.samples().as_u16().ok_or_else(|| {
            CinematicFixtureError::Pipeline("16-bit color pipeline returned 8-bit samples".into())
        })?;
        let png =
            write_png16(config.width, config.height, PngColor::Rgb, samples).map_err(pipeline)?;
        preview_sequence.update(hash_domain("frame", &png).as_bytes());
        write_new(
            &preview_directory.join(format!("frame-{frame:06}.png")),
            &png,
        )?;
    }
    let raw_sequence_identity = raw_sequence.finalize();
    let preview_sequence_identity = preview_sequence.finalize();
    progress("stage=render complete");

    progress("stage=audio begin");
    let audio = build_audio(&trajectory_artifact, config, cx)?;
    let wav_path = sound_directory.join("master.float32.wav");
    write_new(&wav_path, audio.wav_bytes())?;
    let audio_manifest_path = sound_directory.join("master.manifest.json");
    write_new(
        &audio_manifest_path,
        audio.manifest().to_manifest_json().as_bytes(),
    )?;
    audio
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
        timestep_s,
        &run,
        trajectory_receipt.artifact_identity(),
        raw_sequence_identity,
        preview_sequence_identity,
        audio.manifest().wav().wav_identity(),
        over_range_channels,
        gamut_mapped_pixels,
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
        first_preview_path: output_directory.join("preview/frame-000000.png"),
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

fn mechanics_configuration_identity(
    config: &CinematicFixtureConfig,
    profile_identity: ContentHash,
    channels: &CoupledChannelFactors,
    controls: &CoupledControls,
    initial: &CoupledInitialState,
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.mechanics-config.v1");
    hasher.update(profile_identity.as_bytes());
    hasher.update(&config.frames.to_le_bytes());
    for value in [
        channels.gravity_m_per_s2,
        channels.sliding_friction_coefficient,
        channels.rolling_resistance_m,
        channels.contact_stiffness_n_per_m,
        channels.contact_damping_n_s_per_m,
        channels.base_effective_mass_kg,
        channels.base_stiffness_n_per_m,
        channels.base_damping_n_s_per_m,
        channels.gas_rotational_damping_n_m_s,
        channels.gas_translation_damping_n_s_per_m,
        controls.timestep_s,
        controls.terminal_inclination_rad,
        initial.inclination_rad,
        initial.precession_rad_per_s,
        initial.spin_rad_per_s,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.update(&controls.maximum_steps.to_le_bytes());
    hasher.update(&controls.reimpact_limit.to_le_bytes());
    hasher.finalize()
}

fn composition_identity(config: &CinematicFixtureConfig, scene: ContentHash) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-critique.composition.v1");
    hasher.update(scene.as_bytes());
    hasher.update(&config.width.to_le_bytes());
    hasher.update(&config.height.to_le_bytes());
    hasher.update(&config.frames.to_le_bytes());
    hasher.update(&CRITIQUE_FPS.to_le_bytes());
    hasher.update(&config.samples_per_pixel.to_le_bytes());
    hasher.update(&config.max_depth.to_le_bytes());
    hasher.update(b"zero-width-shutter;daily-core-aov;aces-srgb16;exposure-ev-minus-2");
    hasher.finalize()
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
    let eye = Point3::new(0.24, -0.30, 0.18);
    let target = Point3::new(0.0, 0.0, 0.025);
    let physical = PhysicalCamera::try_look_at(
        eye,
        target,
        GeomVec3::new(0.0, 0.0, 1.0),
        CameraProjection::try_half_tangent(0.38)?,
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

fn build_audio(
    trajectory: &EulerRenderTrajectoryArtifact,
    config: &CinematicFixtureConfig,
    cx: &Cx<'_>,
) -> Result<SoundWavArtifact, CinematicFixtureError> {
    let audio_frame_count = u64::from(config.frames) * 2_000;
    let controls = EulerControlStream::try_derive(trajectory.trajectory(), cx).map_err(pipeline)?;
    let preset = representative_modal_preset(RepresentativeDiscMaterial::Tungsten);
    let modal = ModalSynthesisModel::try_new(
        ModalSynthesisModelInput {
            sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
            modes: preset.modes().to_vec(),
            budget: ModalSynthesisBudget::reference_film(audio_frame_count),
        },
        cx,
    )
    .map_err(pipeline)?;
    let mappings = vec![
        SoundExcitationControl {
            channel: SoundExcitationChannel::ContactNormalForce,
            target_component: SoundModalComponent::Disc,
            source_scale: 0.02,
        },
        SoundExcitationControl {
            channel: SoundExcitationChannel::RollingSignedWorkRate,
            target_component: SoundModalComponent::Disc,
            source_scale: 2.0,
        },
        SoundExcitationControl {
            channel: SoundExcitationChannel::BaseDampingSignedWorkRate,
            target_component: SoundModalComponent::BaseAssembly,
            source_scale: 1.0,
        },
    ];
    let interval_count = controls.audio().len();
    let mapper = AudioExcitationMapper::try_new(
        trajectory,
        &controls,
        &modal,
        AudioExcitationModelInput {
            mappings: mappings.clone(),
            reduction: AudioExcitationReduction::RawIntervals,
            spatial_policy: ContactParticipationPolicy::DeclaredStatic,
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
            declared_source_bandwidth_hz: 1_500.0,
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
    apply_terminal_fade(&mut stems, TERMINAL_FADE_SAMPLE_FRAMES)?;
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
    let provisional =
        mix_dry_modal_stems(&stems, mix, AudioArtifactBudget::DEFAULT, cx).map_err(pipeline)?;
    let provisional_meters =
        measure_audio(&provisional, AudioArtifactBudget::DEFAULT, cx).map_err(pipeline)?;
    let provisional_peak = provisional_meters
        .sample_peak_fs
        .max(provisional_meters.true_peak_estimate_fs);
    drop(provisional);
    const TARGET_PEAK_FS: f64 = 0.45;
    if provisional_peak > TARGET_PEAK_FS {
        mix.master_gain_db =
            20.0 * det::ln(TARGET_PEAK_FS / provisional_peak) / core::f64::consts::LN_10;
        if mix.master_gain_db < -120.0 {
            return Err(CinematicFixtureError::Pipeline(format!(
                "mechanics-derived sound needs {:.3} dB attenuation, beyond the admitted mix range",
                mix.master_gain_db
            )));
        }
    }
    SoundWavArtifact::try_build(
        &sound,
        AudioMasterSource::DryModalStems {
            frames: &stems,
            mix,
            source_synthesis: sound.receipt(),
        },
        WavSampleEncoding::Float32,
        WavMetadata::try_new(Some(
            "FrankenSim Euler-disc mechanics-driven critique preview; uncalibrated".to_owned(),
        ))
        .map_err(pipeline)?,
        AudioArtifactBudget::DEFAULT,
        cx,
    )
    .map_err(pipeline)
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
    timestep_s: f64,
    run: &crate::coupled_runner::CoupledRun,
    trajectory_identity: ContentHash,
    raw_sequence_identity: ContentHash,
    preview_sequence_identity: ContentHash,
    wav_identity: ContentHash,
    over_range_channels: u64,
    gamut_mapped_pixels: u64,
    mux: &MuxOutcome,
) -> String {
    let first_sample = run
        .samples
        .first()
        .expect("completed fixture run retains a first sample");
    let last_sample = run
        .samples
        .last()
        .expect("completed fixture run retains a last sample");
    let relative_energy_defect = last_sample.energy_defect_j
        / run
            .checkpoint
            .initial_total_energy_j
            .abs()
            .max(f64::MIN_POSITIVE);
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
    format!(
        concat!(
            "{{\n",
            "  \"schema\": \"frankensim-euler-cinematic-critique-v1\",\n",
            "  \"authority\": \"simulation-derived-visualization-and-physically-informed-sound\",\n",
            "  \"video\": {{\"width\": {}, \"height\": {}, \"frames\": {}, \"fps\": {}, \"duration_s\": {:.9}, \"spp\": {}, \"max_depth\": {}, \"raw_sequence_identity\": \"{}\", \"preview_sequence_identity\": \"{}\", \"over_range_linear_channels\": {}, \"gamut_mapped_pixels\": {}}},\n",
            "  \"mechanics\": {{\"model\": \"closed-profile-reduced-coupled-runner\", \"timestep_s\": {:.9e}, \"sample_count\": {}, \"retained_time_s\": {:.9}, \"terminal\": \"{:?}\", \"initial_rate_selection\": \"rounded gravity-scale estimate; not fitted or calibrated\", \"first_qoi\": {{\"inclination_rad\": {:.17e}, \"precession_rad_per_s\": {:.17e}, \"spin_rad_per_s\": {:.17e}}}, \"last_qoi\": {{\"inclination_rad\": {:.17e}, \"precession_rad_per_s\": {:.17e}, \"spin_rad_per_s\": {:.17e}}}, \"energy\": {{\"initial_total_j\": {:.17e}, \"final_total_j\": {:.17e}, \"defect_j\": {:.17e}, \"relative_defect\": {:.17e}}}, \"trajectory_identity\": \"{}\"}},\n",
            "  \"audio\": {{\"sample_rate_hz\": {}, \"wav_identity\": \"{}\", \"calibrated\": false, \"procedural_texture\": false, \"terminal_fade_sample_frames\": {}, \"mix_policy\": \"single static content-derived gain, peak target 0.45 FS, no limiter\"}},\n",
            "  \"mux\": {},\n",
            "  \"no_claims\": [\"preview timestep and endpoint phase have not been convergence-qualified; inspect the published energy defect\", \"initial rates use a thin-disc gravity scale outside that oracle's admitted geometry for this squat specimen\", \"reduced tangential contact has gross sliding but no static-stick solve\", \"mechanics and acoustic parameters have not been calibrated to experiment\", \"radial spin fiducial is visualization-only and excluded from specimen, contact, and mass mechanics\", \"one-sample-per-pixel preview is intended for motion and composition critique, not final image quality\"]\n",
            "}}\n"
        ),
        config.width,
        config.height,
        config.frames,
        CRITIQUE_FPS,
        duration_s,
        config.samples_per_pixel,
        config.max_depth,
        raw_sequence_identity.to_hex(),
        preview_sequence_identity.to_hex(),
        over_range_channels,
        gamut_mapped_pixels,
        timestep_s,
        run.samples.len(),
        run.checkpoint.time_s,
        run.terminal,
        first_sample.inclination_rad,
        first_sample.precession_rad_per_s,
        first_sample.spin_rad_per_s,
        last_sample.inclination_rad,
        last_sample.precession_rad_per_s,
        last_sample.spin_rad_per_s,
        run.checkpoint.initial_total_energy_j,
        last_sample.mechanical_energy_j,
        last_sample.energy_defect_j,
        relative_energy_defect,
        trajectory_identity.to_hex(),
        SOUND_MASTER_SAMPLE_RATE_HZ,
        wav_identity.to_hex(),
        TERMINAL_FADE_SAMPLE_FRAMES,
        mux_json,
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
    use super::*;

    #[test]
    fn default_is_an_eight_second_practical_preview() {
        let config = CinematicFixtureConfig::default();
        config.validate().unwrap();
        assert_eq!(config.frames, 8 * CRITIQUE_FPS);
        assert_eq!((config.width, config.height), (320, 180));
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
}
