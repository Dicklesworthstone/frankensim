//! Command-line producer for a directly watchable Euler-disc critique clip.
//!
//! It emits replayable intermediate artifacts plus an optional convenience
//! movie; scientific and acoustic no-claims live in the output manifest.

use std::path::PathBuf;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_euler_disc_e2e::cinematic_fixture::{
    CinematicFixtureConfig, CinematicFrameWindow, run_cinematic_fixture,
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
            "--spp" => config.samples_per_pixel = parse(&next_value(&mut args, "--spp")?, "spp")?,
            "--render-seed-salt" => {
                config.render_seed_salt = parse(
                    &next_value(&mut args, "--render-seed-salt")?,
                    "render-seed-salt",
                )?
            }
            "--max-depth" => {
                config.max_depth = parse(&next_value(&mut args, "--max-depth")?, "max-depth")?
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
                     [--workers N] [--tile-width PX] [--tile-height PX] \
                     [--render-memory-mib MIB] [--no-denoise] [--beauty-only-exr] [--dry-audio] \
                     [--no-mux] [--ffmpeg PATH]"
                );
                return Ok(());
            }
            _ => return Err(format!("unknown argument: {argument}")),
        }
    }
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

fn next_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    args.next()
        .ok_or_else(|| format!("missing value for {flag}"))
}

fn parse<T: core::str::FromStr>(value: &str, name: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid {name}: {value}"))
}
