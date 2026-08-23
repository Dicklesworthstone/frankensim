//! Listening receipt for the bowed-string fixture
//! (bead frankensim-music-v8-root-3ez8g.7.5).
//!
//! Renders one open-string bow stroke through the plate-body configuration,
//! encodes it as 16-bit mono WAV at the card sample rate, and emits a JSONL
//! receipt binding the gesture, the WAV bytes by BLAKE3, and the honest
//! no-claim: the automated gates cover the MECHANISM; whether the render
//! "reads as bowed, not plucked" is listener judgment and is deliberately
//! NOT claimed here.
//!
//! Usage:
//! ```text
//! bowed_listen --wav PATH --receipt PATH
//! ```

use std::io::Write as _;
use std::path::PathBuf;

use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination, run_bowed,
};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_scenario::RadiatingPlate;

const SCHEMA: &str = "frankensim.bowed-listen.v1";
const SAMPLE_RATE_HZ: u32 = 48_000;

fn parse_args() -> Result<(PathBuf, PathBuf), String> {
    let mut wav = None;
    let mut receipt = None;
    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--wav" => wav = Some(iter.next().ok_or("--wav needs a path")?),
            "--receipt" => receipt = Some(iter.next().ok_or("--receipt needs a path")?),
            other => return Err(format!("unknown argument {other}")),
        }
    }
    Ok((
        PathBuf::from(wav.ok_or("--wav PATH is required")?),
        PathBuf::from(receipt.ok_or("--receipt PATH is required")?),
    ))
}

fn main() -> std::process::ExitCode {
    let (wav_path, receipt_path) = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("bowed_listen: refusal: {e}");
            return std::process::ExitCode::from(2);
        }
    };
    for dir in [wav_path.parent(), receipt_path.parent()] {
        if let Some(d) = dir {
            if let Err(e) = std::fs::create_dir_all(d) {
                eprintln!("bowed_listen: cannot create {}: {e}", d.display());
                return std::process::ExitCode::from(2);
            }
        }
    }

    let card = BowedStringCard {
        length_m: 0.65,
        tension_n: 60.0,
        linear_density_kg_m: 6.0e-4,
        mode_count: 16,
        zetas: (0..16).map(|k| 1.0e-3 * (1.0 + 0.55 * k as f64)).collect(),
        sample_rate_hz: SAMPLE_RATE_HZ,
    };
    let rosin = StribeckFriction {
        mu_static: 0.8,
        mu_dynamic: 0.4,
        stiction_m_s: 0.04,
    };
    let gesture = match BowGesture::admit(0.45, 3.9, 0.11) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("bowed_listen: gesture refused: {e:?}");
            return std::process::ExitCode::from(2);
        }
    };
    let body_spec = RadiatingPlate {
        area_m2: 3.0e-3,
        mass_kg: 0.15,
        frequency_hz: 280.0,
        damping_ratio: 0.02,
    };
    let body = match CompactBody::from_radiator(body_spec) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("bowed_listen: body refused: {e:?}");
            return std::process::ExitCode::from(2);
        }
    };
    let cfg = BowedRunConfig {
        card,
        island: FrictionIsland::Stribeck(rosin),
        gesture,
        steps: 24_000, // 500 ms open-string stroke
        subsamples: 16,
        termination: Termination::PlateOnePort(Box::new(body)),
        listener_m: 1.0,
    };
    let log = match run_bowed(&cfg) {
        Ok(l) => l,
        Err(e) => {
            eprintln!("bowed_listen: run refused: {e}");
            return std::process::ExitCode::from(2);
        }
    };

    // Peak-normalize the radiated pressure to a modest monitoring level so
    // the WAV is a LISTENING AID, not a calibrated acoustic estimate.
    let peak = log
        .radiated_pressure_pa
        .iter()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    if !(peak > 0.0) {
        eprintln!("bowed_listen: silent radiation; refusing empty render");
        return std::process::ExitCode::from(2);
    }
    let full_scale_pa = peak / 0.3; // -10 dBFS peak
    let (bytes, frames) = match encode_pcm16_wav(&log.radiated_pressure_pa, SAMPLE_RATE_HZ, full_scale_pa)
    {
        Ok(r) => r,
        Err(e) => {
            eprintln!("bowed_listen: encode refused: {e:?}");
            return std::process::ExitCode::from(2);
        }
    };
    if let Err(e) = std::fs::write(&wav_path, &bytes) {
        eprintln!("bowed_listen: write failed: {e}");
        return std::process::ExitCode::from(2);
    }
    let hash = fs_blake3::hash_bytes(&bytes).to_hex();

    let mut sink_file = match std::fs::File::create(&receipt_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("bowed_listen: receipt create failed: {e}");
            return std::process::ExitCode::from(2);
        }
    };
    let row = format!(
        concat!(
            r#"{{"schema":"{SCHEMA}","kind":"listening-receipt","#,
            r#""gesture":{{"v_bow_m_s":0.45,"normal_force_n":3.9,"station_fraction":0.11}},"#,
            r#""termination":"plate-one-port","frames":{frames},"sample_rate_hz":48000,"#,
            r#""peak_normalization_dbfs":-10.0,"wav_sha256":"{hash}","#,
            r#""no_claim":"automated gates cover the emergent mechanism; whether the stroke reads as bowed rather than plucked is listener judgment and is not claimed here"}}"#
        ),
        SCHEMA = SCHEMA,
        frames = frames,
        hash = hash,
    );
    if writeln!(sink_file, "{row}").is_err() {
        eprintln!("bowed_listen: receipt write failed");
        return std::process::ExitCode::from(2);
    }
    eprintln!(
        "bowed_listen: {} frames -> {}; receipt {}",
        frames,
        wav_path.display(),
        receipt_path.display()
    );
    std::process::ExitCode::SUCCESS
}
