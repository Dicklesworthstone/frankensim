//! Fixed-gain transform for canonical FrankenSim stereo WAV artifacts.
//!
//! This command changes only presentation level. It admits the crate's strict
//! WAV subset, multiplies every decoded sample by one caller-declared linear
//! gain, refuses clipping, preserves metadata, and re-encodes through the same
//! canonical codec. It performs no normalization, EQ, dynamics processing, or
//! synthesis.

use std::{
    ffi::OsString,
    fs::{self, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_euler_disc_e2e::{
    AudioArtifactBudget, StereoSample, WavSampleEncoding, decode_stereo_wav, encode_stereo_wav,
    measure_audio,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn main() {
    if let Err(error) = run() {
        eprintln!("status=error message={error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let (input, output, gain) = parse_cli(std::env::args().skip(1))?;
    require_absent(&output, "output")?;
    if input == output {
        return Err("input and output paths must differ".to_owned());
    }

    let input_bytes =
        fs::read(&input).map_err(|error| format!("failed to read {}: {error}", input.display()))?;
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_4741,
                kernel_id: 0x494e_5f47_4149_4e31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let decoded = decode_stereo_wav(&input_bytes, AudioArtifactBudget::DEFAULT, &cx)
            .map_err(|error| format!("input WAV refused: {error:?}"))?;
        if decoded.receipt.encoding() != WavSampleEncoding::Pcm24 {
            return Err(format!(
                "input must be canonical PCM24, decoded {:?}",
                decoded.receipt.encoding()
            ));
        }
        let samples = apply_gain(&decoded.samples, gain)?;
        let meters = measure_audio(&samples, AudioArtifactBudget::DEFAULT, &cx)
            .map_err(|error| format!("output metering refused: {error:?}"))?;
        let (output_bytes, receipt) = encode_stereo_wav(
            &samples,
            decoded.receipt.sample_rate_hz(),
            WavSampleEncoding::Pcm24,
            &decoded.metadata,
            AudioArtifactBudget::DEFAULT,
            &cx,
        )
        .map_err(|error| format!("output WAV refused: {error:?}"))?;
        publish_new(&output, &output_bytes)?;
        println!("status=complete output={}", output.display());
        println!("input_wav_identity={}", decoded.receipt.wav_identity());
        println!("output_wav_identity={}", receipt.wav_identity());
        println!("gain_linear={gain:.17e}");
        println!("sample_peak_fs={:.17e}", meters.sample_peak_fs);
        println!(
            "true_peak_estimate_fs={:.17e}",
            meters.true_peak_estimate_fs
        );
        println!("stereo_rms_fs={:.17e}", meters.stereo_rms_fs);
        if let Some(loudness) = meters.integrated_loudness_lufs {
            println!("integrated_loudness_lufs={loudness:.9}");
        }
        Ok(())
    })
}

fn parse_cli(mut args: impl Iterator<Item = String>) -> Result<(PathBuf, PathBuf, f64), String> {
    let mut input = None;
    let mut output = None;
    let mut gain = None;
    while let Some(argument) = args.next() {
        let value = |flag: &str, args: &mut dyn Iterator<Item = String>| {
            args.next()
                .ok_or_else(|| format!("missing value for {flag}"))
        };
        match argument.as_str() {
            "--input" => input = Some(PathBuf::from(value("--input", &mut args)?)),
            "--output" => output = Some(PathBuf::from(value("--output", &mut args)?)),
            "--gain-linear" => {
                gain = Some(
                    value("--gain-linear", &mut args)?
                        .parse::<f64>()
                        .map_err(|_| "--gain-linear must be a number".to_owned())?,
                )
            }
            "--help" | "-h" => {
                return Err(
                    "usage: euler_cinematic_audio_gain --input WAV --output WAV --gain-linear POSITIVE"
                        .to_owned(),
                );
            }
            _ => return Err(format!("unknown argument: {argument}")),
        }
    }
    let gain = gain.ok_or_else(|| "missing --gain-linear".to_owned())?;
    if !gain.is_finite() || gain <= 0.0 {
        return Err("--gain-linear must be finite and positive".to_owned());
    }
    Ok((
        input.ok_or_else(|| "missing --input".to_owned())?,
        output.ok_or_else(|| "missing --output".to_owned())?,
        gain,
    ))
}

fn apply_gain(samples: &[StereoSample], gain: f64) -> Result<Vec<StereoSample>, String> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(samples.len())
        .map_err(|_| "output sample allocation refused".to_owned())?;
    for (index, sample) in samples.iter().copied().enumerate() {
        let scaled = StereoSample {
            left_fs: sample.left_fs * gain,
            right_fs: sample.right_fs * gain,
        };
        if !scaled.left_fs.is_finite() || !scaled.right_fs.is_finite() {
            return Err(format!(
                "gain produced a non-finite sample at frame {index}"
            ));
        }
        if !(-1.0..=1.0).contains(&scaled.left_fs) || !(-1.0..=1.0).contains(&scaled.right_fs) {
            return Err(format!("gain would clip at frame {index}"));
        }
        output.push(scaled);
    }
    Ok(output)
}

fn publish_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let file_name = path.file_name().ok_or_else(|| {
        format!(
            "output must name a file rather than a filesystem root: {}",
            path.display()
        )
    })?;
    let mut staging_name = OsString::from(".");
    staging_name.push(file_name);
    staging_name.push(format!(".incomplete-{}", std::process::id()));
    let staging = path.with_file_name(staging_name);
    require_absent(&staging, "staging output")?;
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&staging)
        .map_err(|error| format!("failed to create {}: {error}", staging.display()))?;
    file.write_all(bytes)
        .and_then(|()| file.sync_all())
        .map_err(|error| {
            format!(
                "failed to complete {}; incomplete staging file preserved at {}: {error}",
                path.display(),
                staging.display()
            )
        })?;
    drop(file);
    if fs::symlink_metadata(path).is_ok() {
        return Err(format!(
            "output appeared before publication and will not be overwritten; complete staging file preserved at {}",
            staging.display()
        ));
    }
    fs::rename(&staging, path).map_err(|error| {
        format!(
            "failed to publish {} atomically; complete staging file preserved at {}: {error}",
            path.display(),
            staging.display()
        )
    })
}

fn require_absent(path: &Path, label: &str) -> Result<(), String> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(format!(
            "{label} already exists and will not be overwritten: {}",
            path.display()
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(format!(
            "could not verify that {label} is absent at {}: {error}",
            path.display()
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_gain_scales_both_channels_and_refuses_clipping() {
        let samples = [StereoSample {
            left_fs: 0.25,
            right_fs: -0.5,
        }];
        assert_eq!(
            apply_gain(&samples, 2.0).unwrap(),
            vec![StereoSample {
                left_fs: 0.5,
                right_fs: -1.0,
            }]
        );
        assert!(apply_gain(&samples, 2.000_000_1).is_err());
    }

    #[test]
    fn cli_refuses_nonpositive_and_nonfinite_gain() {
        for gain in ["0", "-1", "NaN", "inf"] {
            assert!(
                parse_cli(
                    [
                        "--input",
                        "in.wav",
                        "--output",
                        "out.wav",
                        "--gain-linear",
                        gain
                    ]
                    .into_iter()
                    .map(str::to_owned)
                )
                .is_err()
            );
        }
    }
}
