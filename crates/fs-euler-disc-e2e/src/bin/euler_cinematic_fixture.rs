//! Command-line producer for a directly watchable Euler-disc critique clip.
//!
//! It emits replayable intermediate artifacts plus an optional convenience
//! movie; scientific and acoustic no-claims live in the output manifest.

use std::path::PathBuf;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_euler_disc_e2e::cinematic_fixture::{
    CinematicAdaptiveSamplingConfig, CinematicFixtureConfig, CinematicFrameWindow,
    run_cinematic_fixture,
};
use fs_euler_disc_e2e::render_scene_bridge::{
    MAX_EULER_ARC_SUBDIVISIONS, MAX_EULER_AZIMUTHAL_SEGMENTS,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

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
    let mut adaptive = AdaptiveCliOptions::default();
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
            "--render-seed-salt" => {
                config.render_seed_salt = parse(
                    &next_value(&mut args, "--render-seed-salt")?,
                    "render-seed-salt",
                )?
            }
            "--max-depth" => {
                config.max_depth = parse(&next_value(&mut args, "--max-depth")?, "max-depth")?
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
                     [--frames 192] [--frame-start N --frame-count N --no-mux] \
                     [--spp N] [--render-seed-salt N] [--max-depth N] [--shutter-angle 0..360] \
                     [--azimuthal-segments 8..4096] [--arc-subdivisions 1..1024] \
                     [--adaptive --adaptive-min-spp N --adaptive-max-spp N \
                      --adaptive-batch-spp N --adaptive-abs-error X --adaptive-rel-error X \
                      --adaptive-dark-floor X] \
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

fn required_adaptive<T>(value: Option<T>, flag: &str) -> Result<T, String> {
    value.ok_or_else(|| format!("--adaptive requires {flag}"))
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

fn next_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}"))
}

fn parse<T: core::str::FromStr>(value: &str, name: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
