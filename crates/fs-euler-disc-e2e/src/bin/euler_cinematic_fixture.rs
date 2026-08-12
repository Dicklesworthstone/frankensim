//! Command-line producer for a directly watchable Euler-disc critique clip.
//!
//! It emits replayable intermediate artifacts plus an optional convenience
//! movie; scientific and acoustic no-claims live in the output manifest.

use std::path::PathBuf;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_euler_disc_e2e::cinematic_fixture::{
    CinematicAdaptiveSamplingConfig, CinematicCanonicalMediaSource, CinematicFixtureConfig,
    CinematicFrameWindow, CinematicIndependentPilotSamplingConfig, CinematicInitialMotionConfig,
    CinematicPeriodicSurfaceConfig, CinematicRenderIntegrator, CinematicSurfaceSpectrumConfig,
    run_cinematic_fixture,
};
use fs_euler_disc_e2e::render_scene_bridge::{
    MAX_EULER_ARC_SUBDIVISIONS, MAX_EULER_AZIMUTHAL_SEGMENTS,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_render::tracer::Sampler;
use fs_tribo::{InputAuthority, surface_excitation::SelfAffinePeriodicProfileSpectrum};

fn main() {
    if let Err(error) = run() {
        eprintln!("status=error message={error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let mut config = CinematicFixtureConfig::default();
    let mut output = PathBuf::from("target/euler-cinematic-fixture");
    let mut frame_start = None;
    let mut frame_count = None;
    let mut uniform_spp_seen = false;
    let mut canonical_trajectory_path = None;
    let mut canonical_trajectory_identity = None;
    let mut canonical_listening_master_path = None;
    let mut canonical_listening_master_identity = None;
    let mut adaptive = AdaptiveCliOptions::default();
    let mut pilot = IndependentPilotCliOptions::default();
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--output" => output = PathBuf::from(next_value(&mut args, "--output")?),
            "--width" => config.width = parse(&next_value(&mut args, "--width")?, "width")?,
            "--height" => config.height = parse(&next_value(&mut args, "--height")?, "height")?,
            "--frames" => config.frames = parse(&next_value(&mut args, "--frames")?, "frames")?,
            "--frame-start" => {
                frame_start = Some(parse(
                    &next_value(&mut args, "--frame-start")?,
                    "frame-start",
                )?)
            }
            "--frame-count" => {
                frame_count = Some(parse(
                    &next_value(&mut args, "--frame-count")?,
                    "frame-count",
                )?)
            }
            "--spp" => {
                uniform_spp_seen = true;
                config.samples_per_pixel = parse(&next_value(&mut args, "--spp")?, "spp")?;
            }
            "--adaptive" => adaptive.enabled = true,
            "--adaptive-min-spp" => {
                adaptive.touched = true;
                adaptive.minimum_samples_per_pixel = Some(parse(
                    &next_value(&mut args, "--adaptive-min-spp")?,
                    "adaptive-min-spp",
                )?);
            }
            "--adaptive-max-spp" => {
                adaptive.touched = true;
                adaptive.maximum_samples_per_pixel = Some(parse(
                    &next_value(&mut args, "--adaptive-max-spp")?,
                    "adaptive-max-spp",
                )?);
            }
            "--adaptive-batch-spp" => {
                adaptive.touched = true;
                adaptive.decision_batch_samples = Some(parse(
                    &next_value(&mut args, "--adaptive-batch-spp")?,
                    "adaptive-batch-spp",
                )?);
            }
            "--adaptive-abs-error" => {
                adaptive.touched = true;
                adaptive.absolute_error = Some(parse(
                    &next_value(&mut args, "--adaptive-abs-error")?,
                    "adaptive-abs-error",
                )?);
            }
            "--adaptive-rel-error" => {
                adaptive.touched = true;
                adaptive.relative_error = Some(parse(
                    &next_value(&mut args, "--adaptive-rel-error")?,
                    "adaptive-rel-error",
                )?);
            }
            "--adaptive-dark-floor" => {
                adaptive.touched = true;
                adaptive.dark_floor = Some(parse(
                    &next_value(&mut args, "--adaptive-dark-floor")?,
                    "adaptive-dark-floor",
                )?);
            }
            "--independent-pilot" => pilot.enabled = true,
            "--pilot-spp" => {
                pilot.touched = true;
                pilot.pilot_samples_per_pixel =
                    Some(parse(&next_value(&mut args, "--pilot-spp")?, "pilot-spp")?);
            }
            "--pilot-min-spp" => {
                pilot.touched = true;
                pilot.minimum_samples_per_pixel = Some(parse(
                    &next_value(&mut args, "--pilot-min-spp")?,
                    "pilot-min-spp",
                )?);
            }
            "--pilot-max-spp" => {
                pilot.touched = true;
                pilot.maximum_samples_per_pixel = Some(parse(
                    &next_value(&mut args, "--pilot-max-spp")?,
                    "pilot-max-spp",
                )?);
            }
            "--pilot-abs-error" => {
                pilot.touched = true;
                pilot.absolute_error = Some(parse(
                    &next_value(&mut args, "--pilot-abs-error")?,
                    "pilot-abs-error",
                )?);
            }
            "--pilot-rel-error" => {
                pilot.touched = true;
                pilot.relative_error = Some(parse(
                    &next_value(&mut args, "--pilot-rel-error")?,
                    "pilot-rel-error",
                )?);
            }
            "--pilot-dark-floor" => {
                pilot.touched = true;
                pilot.dark_floor = Some(parse(
                    &next_value(&mut args, "--pilot-dark-floor")?,
                    "pilot-dark-floor",
                )?);
            }
            "--pilot-safety-factor" => {
                pilot.touched = true;
                pilot.safety_factor = Some(parse(
                    &next_value(&mut args, "--pilot-safety-factor")?,
                    "pilot-safety-factor",
                )?);
            }
            "--render-seed-salt" => {
                config.render_seed_salt = parse(
                    &next_value(&mut args, "--render-seed-salt")?,
                    "render-seed-salt",
                )?
            }
            "--sampler" => {
                config.render_sampler = parse_sampler(&next_value(&mut args, "--sampler")?)?
            }
            "--integrator" => {
                config.render_integrator =
                    parse_integrator(&next_value(&mut args, "--integrator")?)?
            }
            "--max-depth" => {
                config.max_depth = parse(&next_value(&mut args, "--max-depth")?, "max-depth")?
            }
            "--initial-inclination-rad" => {
                config.mechanics.initial_motion =
                    parse_initial_inclination(&next_value(&mut args, "--initial-inclination-rad")?)?
            }
            "--support-plate-thickness-m" => {
                config.support_plate.thickness_m = parse_support_plate_thickness(&next_value(
                    &mut args,
                    "--support-plate-thickness-m",
                )?)?
            }
            "--mechanics-sample-rate" => {
                config.mechanics.sample_rate_hz = parse(
                    &next_value(&mut args, "--mechanics-sample-rate")?,
                    "mechanics-sample-rate",
                )?
            }
            "--mechanics-preroll-frames" => {
                config.mechanics_preroll_video_frames = parse(
                    &next_value(&mut args, "--mechanics-preroll-frames")?,
                    "mechanics-preroll-frames",
                )?
            }
            "--canonical-trajectory" => {
                canonical_trajectory_path = Some(PathBuf::from(next_value(
                    &mut args,
                    "--canonical-trajectory",
                )?));
            }
            "--canonical-trajectory-identity" => {
                canonical_trajectory_identity = Some(parse_content_hash(
                    &next_value(&mut args, "--canonical-trajectory-identity")?,
                    "canonical-trajectory-identity",
                )?);
            }
            "--canonical-listening-master" => {
                canonical_listening_master_path = Some(PathBuf::from(next_value(
                    &mut args,
                    "--canonical-listening-master",
                )?));
            }
            "--canonical-listening-master-identity" => {
                canonical_listening_master_identity = Some(parse_content_hash(
                    &next_value(&mut args, "--canonical-listening-master-identity")?,
                    "canonical-listening-master-identity",
                )?);
            }
            "--listening-gain-fs-per-pa" => {
                config.listening_gain_fs_per_pa = parse(
                    &next_value(&mut args, "--listening-gain-fs-per-pa")?,
                    "listening-gain-fs-per-pa",
                )?
            }
            "--disc-surface-self-affine" => apply_self_affine_surface(
                &mut config.mechanics.disc_surface,
                "disc",
                &next_value(&mut args, "--disc-surface-self-affine")?,
            )?,
            "--support-surface-self-affine" => apply_self_affine_surface(
                &mut config.mechanics.support_surface,
                "support",
                &next_value(&mut args, "--support-surface-self-affine")?,
            )?,
            "--disc-surface-samples" => {
                config.mechanics.disc_surface.sample_count = parse(
                    &next_value(&mut args, "--disc-surface-samples")?,
                    "disc-surface-samples",
                )?
            }
            "--support-surface-samples" => {
                config.mechanics.support_surface.sample_count = parse(
                    &next_value(&mut args, "--support-surface-samples")?,
                    "support-surface-samples",
                )?
            }
            "--azimuthal-segments" => {
                let value = next_value(&mut args, "--azimuthal-segments")?;
                apply_tessellation_control(&mut config, "--azimuthal-segments", &value)?;
            }
            "--arc-subdivisions" => {
                let value = next_value(&mut args, "--arc-subdivisions")?;
                apply_tessellation_control(&mut config, "--arc-subdivisions", &value)?;
            }
            "--shutter-angle" => {
                config.shutter_angle_degrees =
                    parse(&next_value(&mut args, "--shutter-angle")?, "shutter-angle")?
            }
            "--workers" => {
                config.render_workers = parse(&next_value(&mut args, "--workers")?, "workers")?
            }
            "--tile-width" => {
                config.tile_width = parse(&next_value(&mut args, "--tile-width")?, "tile-width")?
            }
            "--tile-height" => {
                config.tile_height = parse(&next_value(&mut args, "--tile-height")?, "tile-height")?
            }
            "--render-memory-mib" => {
                let mib: u64 = parse(
                    &next_value(&mut args, "--render-memory-mib")?,
                    "render-memory-mib",
                )?;
                config.render_memory_limit_bytes = mib
                    .checked_mul(1024 * 1024)
                    .ok_or_else(|| "render-memory-mib overflows bytes".to_owned())?;
            }
            "--no-denoise" => config.denoise_previews = false,
            "--beauty-only-exr" => config.retain_full_aov_exr = false,
            "--dry-audio" => config.spatialize_audio = false,
            "--no-mux" => config.mux_with_ffmpeg = false,
            "--ffmpeg" => {
                config.ffmpeg_executable = PathBuf::from(next_value(&mut args, "--ffmpeg")?)
            }
            "--help" | "-h" => {
                println!(
                    "Usage: euler_cinematic_fixture [--output DIR] [--width PX] [--height PX] \
                     [--frames 1..=192] [--frame-start N --frame-count N --no-mux] \
                     [--spp N] [--sampler iid|owen-pixel|owen-full-path] \
                     [--integrator camera-path|bdpt|bdpt-camera-connected] \
                     [--render-seed-salt N] [--max-depth N] [--shutter-angle 0..360] \
                     [--initial-inclination-rad 0<THETA<=0.35] \
                     [--support-plate-thickness-m POSITIVE_METRES] \
                     [--mechanics-sample-rate HZ] [--mechanics-preroll-frames 1..240] \
                     [--canonical-trajectory FSET --canonical-trajectory-identity BLAKE3 \
                      --canonical-listening-master WAV \
                      --canonical-listening-master-identity BLAKE3] \
                     [--listening-gain-fs-per-pa POSITIVE] \
                     [--disc-surface-self-affine RMS_M,HURST,MIN_CYCLE,MAX_CYCLE,SEED] \
                     [--support-surface-self-affine RMS_M,HURST,MIN_CYCLE,MAX_CYCLE,SEED] \
                     [--disc-surface-samples N] [--support-surface-samples N] \
                     [--azimuthal-segments 8..4096] [--arc-subdivisions 1..1024] \
                     [--adaptive --adaptive-min-spp N --adaptive-max-spp N \
                      --adaptive-batch-spp N --adaptive-abs-error X --adaptive-rel-error X \
                      --adaptive-dark-floor X] \
                     [--independent-pilot --pilot-spp N --pilot-min-spp N --pilot-max-spp N \
                      --pilot-abs-error X --pilot-rel-error X --pilot-dark-floor X \
                      --pilot-safety-factor X] \
                     [--workers N] [--tile-width PX] [--tile-height PX] \
                     [--render-memory-mib MIB] [--no-denoise] [--beauty-only-exr] [--dry-audio] \
                     [--no-mux] [--ffmpeg PATH]"
                );
                return Ok(());
            }
            _ => return Err(format!("unknown argument: {argument}")),
        }
    }
    adaptive.apply(&mut config, uniform_spp_seen)?;
    pilot.apply(&mut config, uniform_spp_seen)?;
    config.canonical_media_source = match (
        canonical_trajectory_path,
        canonical_trajectory_identity,
        canonical_listening_master_path,
        canonical_listening_master_identity,
    ) {
        (None, None, None, None) => None,
        (
            Some(trajectory_path),
            Some(trajectory_identity),
            Some(listening_master_path),
            Some(listening_master_identity),
        ) => Some(CinematicCanonicalMediaSource {
            trajectory_path,
            trajectory_identity,
            listening_master_path,
            listening_master_identity,
        }),
        _ => {
            return Err(concat!(
                "canonical replay requires all four flags: --canonical-trajectory, ",
                "--canonical-trajectory-identity, --canonical-listening-master, and ",
                "--canonical-listening-master-identity"
            )
            .to_owned());
        }
    };
    match (frame_start, frame_count) {
        (None, None) => {}
        (Some(first_frame), Some(frame_count)) => {
            config.frame_window = CinematicFrameWindow::Range {
                first_frame,
                frame_count,
            };
        }
        _ => {
            return Err("--frame-start and --frame-count must be supplied together".to_owned());
        }
    }

    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_434c,
                kernel_id: 0x4649_5854,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let report = run_cinematic_fixture(&config, &output, &cx, |message| {
            eprintln!("status=progress {message}");
        })
        .map_err(|error| error.to_string())?;
        println!(
            "status=complete manifest={}",
            report.manifest_path.display()
        );
        println!("wav={}", report.wav_path.display());
        if let Some(movie) = report.movie_path {
            println!("movie={}", movie.display());
        }
        Ok(())
    })
}

#[derive(Default)]
struct AdaptiveCliOptions {
    enabled: bool,
    touched: bool,
    minimum_samples_per_pixel: Option<u32>,
    maximum_samples_per_pixel: Option<u32>,
    decision_batch_samples: Option<u32>,
    absolute_error: Option<f64>,
    relative_error: Option<f64>,
    dark_floor: Option<f64>,
}

impl AdaptiveCliOptions {
    fn apply(
        self,
        config: &mut CinematicFixtureConfig,
        uniform_spp_seen: bool,
    ) -> Result<(), String> {
        if !self.enabled {
            if self.touched {
                return Err("adaptive controls require the explicit --adaptive opt-in".to_owned());
            }
            return Ok(());
        }
        if uniform_spp_seen {
            return Err(
                "--spp selects uniform sampling and cannot be combined with --adaptive; use --adaptive-max-spp"
                    .to_owned(),
            );
        }
        config.adaptive_sampling = Some(CinematicAdaptiveSamplingConfig {
            minimum_samples_per_pixel: required_adaptive(
                self.minimum_samples_per_pixel,
                "--adaptive-min-spp",
            )?,
            maximum_samples_per_pixel: required_adaptive(
                self.maximum_samples_per_pixel,
                "--adaptive-max-spp",
            )?,
            decision_batch_samples: required_adaptive(
                self.decision_batch_samples,
                "--adaptive-batch-spp",
            )?,
            absolute_error: required_adaptive(self.absolute_error, "--adaptive-abs-error")?,
            relative_error: required_adaptive(self.relative_error, "--adaptive-rel-error")?,
            dark_floor: required_adaptive(self.dark_floor, "--adaptive-dark-floor")?,
        });
        Ok(())
    }
}

#[derive(Default)]
struct IndependentPilotCliOptions {
    enabled: bool,
    touched: bool,
    pilot_samples_per_pixel: Option<u32>,
    minimum_samples_per_pixel: Option<u32>,
    maximum_samples_per_pixel: Option<u32>,
    absolute_error: Option<f64>,
    relative_error: Option<f64>,
    dark_floor: Option<f64>,
    safety_factor: Option<f64>,
}

impl IndependentPilotCliOptions {
    fn apply(
        self,
        config: &mut CinematicFixtureConfig,
        uniform_spp_seen: bool,
    ) -> Result<(), String> {
        if !self.enabled {
            if self.touched {
                return Err(
                    "pilot controls require the explicit --independent-pilot opt-in".to_owned(),
                );
            }
            return Ok(());
        }
        if uniform_spp_seen {
            return Err(
                "--spp selects uniform sampling and cannot be combined with --independent-pilot; use --pilot-max-spp"
                    .to_owned(),
            );
        }
        if config.adaptive_sampling.is_some() {
            return Err("--adaptive and --independent-pilot are mutually exclusive".to_owned());
        }
        config.independent_pilot_sampling = Some(CinematicIndependentPilotSamplingConfig {
            pilot_samples_per_pixel: required_pilot(self.pilot_samples_per_pixel, "--pilot-spp")?,
            minimum_samples_per_pixel: required_pilot(
                self.minimum_samples_per_pixel,
                "--pilot-min-spp",
            )?,
            maximum_samples_per_pixel: required_pilot(
                self.maximum_samples_per_pixel,
                "--pilot-max-spp",
            )?,
            absolute_error: required_pilot(self.absolute_error, "--pilot-abs-error")?,
            relative_error: required_pilot(self.relative_error, "--pilot-rel-error")?,
            dark_floor: required_pilot(self.dark_floor, "--pilot-dark-floor")?,
            safety_factor: required_pilot(self.safety_factor, "--pilot-safety-factor")?,
        });
        Ok(())
    }
}

fn required_adaptive<T>(value: Option<T>, flag: &str) -> Result<T, String> {
    value.ok_or_else(|| format!("--adaptive requires {flag}"))
}

fn required_pilot<T>(value: Option<T>, flag: &str) -> Result<T, String> {
    value.ok_or_else(|| format!("--independent-pilot requires {flag}"))
}

fn apply_tessellation_control(
    config: &mut CinematicFixtureConfig,
    flag: &str,
    value: &str,
) -> Result<(), String> {
    match flag {
        "--azimuthal-segments" => {
            config.azimuthal_segments =
                parse_bounded_u32(value, "azimuthal-segments", 8, MAX_EULER_AZIMUTHAL_SEGMENTS)?;
        }
        "--arc-subdivisions" => {
            config.arc_subdivisions_per_arc =
                parse_bounded_u32(value, "arc-subdivisions", 1, MAX_EULER_ARC_SUBDIVISIONS)?;
        }
        _ => return Err(format!("unknown tessellation control: {flag}")),
    }
    Ok(())
}

fn parse_bounded_u32(value: &str, name: &str, minimum: u32, maximum: u32) -> Result<u32, String> {
    let parsed = parse(value, name)?;
    if !(minimum..=maximum).contains(&parsed) {
        return Err(format!("{name} must be in {minimum}..={maximum}: {value}"));
    }
    Ok(parsed)
}

fn parse_sampler(value: &str) -> Result<Sampler, String> {
    match value {
        "iid" => Ok(Sampler::Iid),
        "owen-pixel" => Ok(Sampler::OwenSobol),
        "owen-full-path" => Ok(Sampler::OwenSobolFullPath),
        _ => Err(format!(
            "invalid sampler: {value}; expected iid, owen-pixel, or owen-full-path"
        )),
    }
}

fn parse_integrator(value: &str) -> Result<CinematicRenderIntegrator, String> {
    match value {
        "camera-path" => Ok(CinematicRenderIntegrator::CameraPathMis),
        "bdpt" => Ok(CinematicRenderIntegrator::Bidirectional),
        "bdpt-camera-connected" => Ok(CinematicRenderIntegrator::BidirectionalCameraConnected),
        _ => Err(format!(
            "invalid integrator: {value}; expected camera-path, bdpt, or bdpt-camera-connected"
        )),
    }
}

fn parse_initial_inclination(value: &str) -> Result<CinematicInitialMotionConfig, String> {
    let inclination_rad: f64 = parse(value, "initial-inclination-rad")?;
    if !(inclination_rad.is_finite() && inclination_rad > 0.0 && inclination_rad <= 0.35) {
        return Err(format!(
            "initial-inclination-rad must be finite and in (0, 0.35]: {value}"
        ));
    }
    Ok(CinematicInitialMotionConfig::SmallAngleNoSlip { inclination_rad })
}

fn parse_support_plate_thickness(value: &str) -> Result<f64, String> {
    let thickness_m: f64 = parse(value, "support-plate-thickness-m")?;
    if !(thickness_m.is_finite() && thickness_m > 0.0) {
        return Err(format!(
            "support-plate-thickness-m must be finite and positive: {value}"
        ));
    }
    Ok(thickness_m)
}

fn apply_self_affine_surface(
    surface: &mut CinematicPeriodicSurfaceConfig,
    surface_name: &str,
    value: &str,
) -> Result<(), String> {
    let fields = value.split(',').collect::<Vec<_>>();
    if fields.len() != 5 {
        return Err(format!(
            "invalid {surface_name} self-affine surface: expected RMS_M,HURST,MIN_CYCLE,MAX_CYCLE,SEED: {value}"
        ));
    }
    let rms_height_m = parse(fields[0], &format!("{surface_name}-surface-rms-height-m"))?;
    let hurst_exponent = parse(fields[1], &format!("{surface_name}-surface-hurst"))?;
    let minimum_cycles_per_track = parse(fields[2], &format!("{surface_name}-surface-min-cycle"))?;
    let maximum_cycles_per_track = parse(fields[3], &format!("{surface_name}-surface-max-cycle"))?;
    let phase_seed = parse(fields[4], &format!("{surface_name}-surface-seed"))?;
    let spectrum = SelfAffinePeriodicProfileSpectrum::new(
        rms_height_m,
        hurst_exponent,
        minimum_cycles_per_track,
        maximum_cycles_per_track,
        phase_seed,
    )
    .map_err(|error| format!("invalid {surface_name} self-affine surface: {error}"))?;
    surface.spectrum = CinematicSurfaceSpectrumConfig::SelfAffine(spectrum);
    surface.authority = InputAuthority::CallerDeclared;
    surface.source_id = format!(
        "cli:caller-declared-self-affine:{surface_name}:rms={rms_height_m:.17e}:hurst={hurst_exponent:.17e}:cycles={minimum_cycles_per_track}..={maximum_cycles_per_track}:seed={phase_seed}"
    );
    Ok(())
}

fn next_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}"))
}

fn parse<T: core::str::FromStr>(value: &str, name: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}

fn parse_content_hash(value: &str, name: &str) -> Result<ContentHash, String> {
    ContentHash::from_hex(value)
        .ok_or_else(|| format!("invalid {name}: expected exactly 64 hexadecimal characters"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_media_identities_are_exact_content_hashes() {
        let hex = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
        assert_eq!(parse_content_hash(hex, "trajectory").unwrap().to_hex(), hex);
        assert_eq!(
            parse_content_hash("0123", "trajectory").unwrap_err(),
            "invalid trajectory: expected exactly 64 hexadecimal characters"
        );
        assert_eq!(
            parse_content_hash(
                "g123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
                "trajectory"
            )
            .unwrap_err(),
            "invalid trajectory: expected exactly 64 hexadecimal characters"
        );
    }

    #[test]
    fn integrator_cli_names_are_explicit_and_fail_closed() {
        assert_eq!(
            parse_integrator("camera-path").unwrap(),
            CinematicRenderIntegrator::CameraPathMis
        );
        assert_eq!(
            parse_integrator("bdpt").unwrap(),
            CinematicRenderIntegrator::Bidirectional
        );
        assert_eq!(
            parse_integrator("bdpt-camera-connected").unwrap(),
            CinematicRenderIntegrator::BidirectionalCameraConnected
        );
        assert!(parse_integrator("path").is_err());
    }

    fn complete_adaptive_options() -> AdaptiveCliOptions {
        AdaptiveCliOptions {
            enabled: true,
            touched: true,
            minimum_samples_per_pixel: Some(8),
            maximum_samples_per_pixel: Some(64),
            decision_batch_samples: Some(4),
            absolute_error: Some(1.0e-4),
            relative_error: Some(0.02),
            dark_floor: Some(1.0e-5),
        }
    }

    fn complete_pilot_options() -> IndependentPilotCliOptions {
        IndependentPilotCliOptions {
            enabled: true,
            touched: true,
            pilot_samples_per_pixel: Some(16),
            minimum_samples_per_pixel: Some(16),
            maximum_samples_per_pixel: Some(256),
            absolute_error: Some(1.0e-4),
            relative_error: Some(0.02),
            dark_floor: Some(1.0e-5),
            safety_factor: Some(2.0),
        }
    }

    #[test]
    fn adaptive_cli_requires_explicit_complete_controls() {
        let mut missing = complete_adaptive_options();
        missing.relative_error = None;
        let error = missing
            .apply(&mut CinematicFixtureConfig::default(), false)
            .unwrap_err();
        assert_eq!(error, "--adaptive requires --adaptive-rel-error");

        let mut unmarked = complete_adaptive_options();
        unmarked.enabled = false;
        let error = unmarked
            .apply(&mut CinematicFixtureConfig::default(), false)
            .unwrap_err();
        assert_eq!(
            error,
            "adaptive controls require the explicit --adaptive opt-in"
        );
    }

    #[test]
    fn adaptive_cli_refuses_ambiguous_uniform_spp_and_builds_exact_policy() {
        let error = complete_adaptive_options()
            .apply(&mut CinematicFixtureConfig::default(), true)
            .unwrap_err();
        assert!(error.contains("--spp selects uniform sampling"));

        let mut config = CinematicFixtureConfig::default();
        complete_adaptive_options()
            .apply(&mut config, false)
            .unwrap();
        config.validate().unwrap();
        assert!(config.denoise_previews);
        assert!(config.retain_full_aov_exr);
        assert_eq!(
            config.adaptive_sampling,
            Some(CinematicAdaptiveSamplingConfig {
                minimum_samples_per_pixel: 8,
                maximum_samples_per_pixel: 64,
                decision_batch_samples: 4,
                absolute_error: 1.0e-4,
                relative_error: 0.02,
                dark_floor: 1.0e-5,
            })
        );
    }

    #[test]
    fn independent_pilot_cli_requires_opt_in_and_complete_controls() {
        let mut missing = complete_pilot_options();
        missing.safety_factor = None;
        let error = missing
            .apply(&mut CinematicFixtureConfig::default(), false)
            .unwrap_err();
        assert_eq!(error, "--independent-pilot requires --pilot-safety-factor");

        let mut unmarked = complete_pilot_options();
        unmarked.enabled = false;
        let error = unmarked
            .apply(&mut CinematicFixtureConfig::default(), false)
            .unwrap_err();
        assert_eq!(
            error,
            "pilot controls require the explicit --independent-pilot opt-in"
        );
    }

    #[test]
    fn independent_pilot_cli_refuses_competing_modes_and_builds_exact_policy() {
        let error = complete_pilot_options()
            .apply(&mut CinematicFixtureConfig::default(), true)
            .unwrap_err();
        assert!(error.contains("--spp selects uniform sampling"));

        let mut adaptive = CinematicFixtureConfig::default();
        adaptive.adaptive_sampling = Some(CinematicAdaptiveSamplingConfig {
            minimum_samples_per_pixel: 8,
            maximum_samples_per_pixel: 64,
            decision_batch_samples: 4,
            absolute_error: 1.0e-4,
            relative_error: 0.02,
            dark_floor: 1.0e-5,
        });
        let error = complete_pilot_options()
            .apply(&mut adaptive, false)
            .unwrap_err();
        assert_eq!(
            error,
            "--adaptive and --independent-pilot are mutually exclusive"
        );

        let mut config = CinematicFixtureConfig::default();
        complete_pilot_options().apply(&mut config, false).unwrap();
        config.validate().unwrap();
        assert_eq!(
            config.independent_pilot_sampling,
            Some(CinematicIndependentPilotSamplingConfig {
                pilot_samples_per_pixel: 16,
                minimum_samples_per_pixel: 16,
                maximum_samples_per_pixel: 256,
                absolute_error: 1.0e-4,
                relative_error: 0.02,
                dark_floor: 1.0e-5,
                safety_factor: 2.0,
            })
        );
    }

    #[test]
    fn tessellation_cli_controls_parse_apply_and_enforce_scene_bounds() {
        let mut config = CinematicFixtureConfig::default();
        apply_tessellation_control(&mut config, "--azimuthal-segments", "512").unwrap();
        apply_tessellation_control(&mut config, "--arc-subdivisions", "64").unwrap();
        assert_eq!(config.azimuthal_segments, 512);
        assert_eq!(config.arc_subdivisions_per_arc, 64);
        config.validate().unwrap();

        assert_eq!(
            apply_tessellation_control(&mut config, "--azimuthal-segments", "7").unwrap_err(),
            "azimuthal-segments must be in 8..=4096: 7"
        );
        assert_eq!(
            apply_tessellation_control(&mut config, "--arc-subdivisions", "1025").unwrap_err(),
            "arc-subdivisions must be in 1..=1024: 1025"
        );
        assert_eq!(
            apply_tessellation_control(&mut config, "--arc-subdivisions", "not-a-number")
                .unwrap_err(),
            "invalid arc-subdivisions: not-a-number"
        );
    }

    #[test]
    fn sampler_cli_names_are_explicit_and_reject_unknown_estimators() {
        assert_eq!(parse_sampler("iid").unwrap(), Sampler::Iid);
        assert_eq!(parse_sampler("owen-pixel").unwrap(), Sampler::OwenSobol);
        assert_eq!(
            parse_sampler("owen-full-path").unwrap(),
            Sampler::OwenSobolFullPath
        );
        assert_eq!(
            parse_sampler("sobol").unwrap_err(),
            "invalid sampler: sobol; expected iid, owen-pixel, or owen-full-path"
        );
    }

    #[test]
    fn initial_inclination_cli_is_an_explicit_bounded_initial_condition() {
        assert_eq!(
            parse_initial_inclination("0.02").unwrap(),
            CinematicInitialMotionConfig::SmallAngleNoSlip {
                inclination_rad: 0.02,
            }
        );
        assert_eq!(
            parse_initial_inclination("0").unwrap_err(),
            "initial-inclination-rad must be finite and in (0, 0.35]: 0"
        );
        assert_eq!(
            parse_initial_inclination("NaN").unwrap_err(),
            "initial-inclination-rad must be finite and in (0, 0.35]: NaN"
        );
    }

    #[test]
    fn support_plate_thickness_cli_changes_the_shared_physical_plate() {
        assert_eq!(parse_support_plate_thickness("0.020").unwrap(), 0.020);
        assert_eq!(
            parse_support_plate_thickness("0").unwrap_err(),
            "support-plate-thickness-m must be finite and positive: 0"
        );
        assert_eq!(
            parse_support_plate_thickness("NaN").unwrap_err(),
            "support-plate-thickness-m must be finite and positive: NaN"
        );
    }

    #[test]
    fn self_affine_surface_cli_changes_upstream_contact_geometry() {
        let mut config = CinematicFixtureConfig::default();
        apply_self_affine_surface(
            &mut config.mechanics.disc_surface,
            "disc",
            "4.5e-10,0.8,3,96,24301",
        )
        .unwrap();
        assert_eq!(
            config.mechanics.disc_surface.authority,
            InputAuthority::CallerDeclared
        );
        let CinematicSurfaceSpectrumConfig::SelfAffine(spectrum) =
            config.mechanics.disc_surface.spectrum
        else {
            panic!("CLI must configure self-affine contact geometry");
        };
        assert_eq!(spectrum.rms_height_m(), 4.5e-10);
        assert_eq!(spectrum.hurst_exponent(), 0.8);
        assert_eq!(spectrum.minimum_cycles_per_track(), 3);
        assert_eq!(spectrum.maximum_cycles_per_track(), 96);
        assert_eq!(spectrum.phase_seed(), 24_301);
        config.validate().unwrap();

        assert!(
            apply_self_affine_surface(
                &mut config.mechanics.support_surface,
                "support",
                "4.5e-10,0.8,3,96",
            )
            .unwrap_err()
            .contains("expected RMS_M,HURST,MIN_CYCLE,MAX_CYCLE,SEED")
        );
        assert!(
            apply_self_affine_surface(
                &mut config.mechanics.support_surface,
                "support",
                "4.5e-10,1.2,3,96,24302",
            )
            .unwrap_err()
            .contains("invalid support self-affine surface")
        );
    }
}
