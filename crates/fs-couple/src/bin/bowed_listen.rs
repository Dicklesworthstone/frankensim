//! Render physical bow performances through the existing plate-body fixture.
//! The default remains a 500 ms constant open-string stroke. A canonical gesture
//! file adds ramps, reversals, station changes, release and re-entry, with one
//! retained string/body state and callback-independent sample timing.
//!
//! The peak-normalized WAV is a listening aid, not calibrated acoustic output.
//! The existing first-order stiction and one-way bridge/body approximation remain;
//! perceived realism and experimental validity are not asserted by this command.
//!
//! ```text
//! bowed_listen --wav out.wav --receipt out.jsonl
//! bowed_listen --gesture crates/fs-couple/examples/bow-phrases.gesture \
//!   --track bow --samples 24000 --block 256 --wav out.wav --receipt out.jsonl
//! ```

use std::io::{Read as _, Write as _};
use std::path::{Path, PathBuf};

use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
};
use fs_couple::bowed_string::runtime::BowedStringState;
use fs_couple::bowed_string::runtime::schedule::ScheduledBowedRenderer;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::RadiatingPlate;
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const SCHEMA: &str = "frankensim.bowed-listen.v2";
const SAMPLE_RATE_HZ: u32 = 48_000;
const MAX_SAMPLES: usize = 2_880_000; // Explicit 60-second offline render bound.
const MAX_GESTURE_BYTES: usize = 262_144;
const MAX_WORK: u64 = 10_000_000;
const MAX_EVENTS: usize = 1_000_000;

#[derive(Debug)]
struct Options {
    wav: PathBuf,
    receipt: PathBuf,
    gesture: Option<PathBuf>,
    track: String,
    samples: usize,
    block: usize,
}

fn positive_limit(value: &str, name: &str, max: usize) -> Result<usize, String> {
    let n = value.parse::<usize>().map_err(|_| format!("{name} requires an integer"))?;
    if n == 0 || n > max {
        return Err(format!("{name} must be in 1..={max}"));
    }
    Ok(n)
}

fn parse_args(mut iter: impl Iterator<Item = String>) -> Result<Options, String> {
    let mut wav = None;
    let mut receipt = None;
    let mut gesture = None;
    let mut track = "bow".to_string();
    let mut samples = 24_000;
    let mut block = 256;
    while let Some(arg) = iter.next() {
        let value = iter.next().ok_or_else(|| format!("{arg} needs a value"))?;
        match arg.as_str() {
            "--wav" => wav = Some(PathBuf::from(value)),
            "--receipt" => receipt = Some(PathBuf::from(value)),
            "--gesture" => gesture = Some(PathBuf::from(value)),
            "--track" => track = value,
            "--samples" => samples = positive_limit(&value, "--samples", MAX_SAMPLES)?,
            "--block" => block = positive_limit(&value, "--block", 8192)?,
            other => return Err(format!("unknown argument {other}")),
        }
    }
    let wav = wav.ok_or("--wav PATH is required")?;
    let receipt = receipt.ok_or("--receipt PATH is required")?;
    if wav == receipt || gesture.as_ref() == Some(&wav) || gesture.as_ref() == Some(&receipt) {
        return Err("input gesture, WAV and receipt must use distinct paths".to_string());
    }
    if track.trim().is_empty() {
        return Err("--track must not be empty".to_string());
    }
    Ok(Options { wav, receipt, gesture, track, samples, block })
}

fn decode_schedule(bytes: &[u8]) -> Result<GestureSchedule, String> {
    if bytes.len() > MAX_GESTURE_BYTES {
        return Err("gesture input exceeds 256 KiB".to_string());
    }
    let text = std::str::from_utf8(bytes).map_err(|_| "gesture input is not UTF-8".to_string())?;
    let line_count = text.lines().count();
    // This command has one physical voice; never silently drop other tracks.
    // Preflight declared counts before the existing decoder reserves its vectors.
    if text.lines().nth(2) != Some("tracks\t1") {
        return Err("bowed_listen requires exactly one gesture track".to_string());
    }
    for line in text.lines() {
        if let Some(count) = line.strip_prefix("events\t") {
            let count = count.parse::<usize>().map_err(|_| "invalid event count".to_string())?;
            if count > line_count {
                return Err("declared events exceed the available source lines".to_string());
            }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(bytes).map_err(|e| e.to_string())?;
    if !matches!(schedule.tracks()[0].target, GestureTarget::BowStroke { .. }) {
        return Err("the gesture track must target a bow stroke".to_string());
    }
    Ok(schedule)
}

fn load_schedule(path: Option<&Path>, track: &str, initial: BowGesture) -> Result<GestureSchedule, String> {
    let schedule = if let Some(path) = path {
        let mut bytes = Vec::new();
        std::fs::File::open(path).map_err(|e| format!("gesture open: {e}"))?
            .take(MAX_GESTURE_BYTES as u64 + 1).read_to_end(&mut bytes)
            .map_err(|e| format!("gesture read: {e}"))?;
        decode_schedule(&bytes)?
    } else {
        GestureSchedule::try_new(700, vec![GestureTrack {
            id: track.to_string(), target: GestureTarget::BowStroke { string: 0 },
            initial: GestureValue::Bow {
                velocity_m_per_s: initial.v_bow_m_s,
                normal_force_n: initial.normal_force_n,
                station: initial.station_fraction,
            },
            events: Vec::new(),
        }]).map_err(|e| e.to_string())?
    };
    if schedule.tracks()[0].id != track {
        return Err(format!("gesture has no selected track {track:?}"));
    }
    Ok(schedule)
}

fn run(options: Options) -> Result<(), String> {
    let card = BowedStringCard {
        length_m: 0.65,
        tension_n: 60.0,
        linear_density_kg_m: 6.0e-4,
        bending_stiffness_n_m2: 0.0,
        viscous_bending_n_m2_s: 0.0,
        mode_count: 16,
        zetas: (0..16).map(|k| 1.0e-3 * (1.0 + 0.55 * f64::from(k))).collect(),
        sample_rate_hz: SAMPLE_RATE_HZ,
    };
    let rosin = StribeckFriction { mu_static: 0.8, mu_dynamic: 0.4, stiction_m_s: 0.04 };
    let gesture = BowGesture::admit(0.45, 3.9, 0.11).map_err(|e| format!("gesture: {e:?}"))?;
    let schedule = load_schedule(options.gesture.as_deref(), &options.track, gesture)?;
    let body = CompactBody::from_radiator(RadiatingPlate {
        area_m2: 3.0e-3, mass_kg: 0.15, frequency_hz: 280.0, damping_ratio: 0.02,
    }).map_err(|e| format!("body: {e:?}"))?;
    let ambient = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)
        .map_err(|e| format!("ambient: {e}"))?;
    let cfg = BowedRunConfig {
        card, island: FrictionIsland::Stribeck(rosin), gesture,
        steps: options.samples, subsamples: 16,
        termination: Termination::PlateOnePort { body: Box::new(body), ambient },
        listener_m: 1.0,
    };
    let state = BowedStringState::new(&cfg, options.block).map_err(|e| format!("voice: {e}"))?;
    let mut render = ScheduledBowedRenderer::new(
        state, &schedule, &options.track, options.samples as u64, MAX_WORK, MAX_EVENTS,
    ).map_err(|e| e.to_string())?;
    let mut pressure = vec![0.0; options.samples];
    render.pressure_under_gate(&CancelGate::new(), &mut pressure, options.block)
        .map_err(|e| e.to_string())?;

    // Normalize only after the physical run: monitoring gain never feeds physics.
    let peak = pressure.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    if !peak.is_finite() || peak <= 0.0 {
        return Err("silent or non-finite radiation; refusing empty render".to_string());
    }
    let full_scale_pa = peak / 0.3;
    let (bytes, clipped) = encode_pcm16_wav(&pressure, SAMPLE_RATE_HZ, full_scale_pa)
        .map_err(|e| format!("encode: {e:?}"))?;
    let hash = fs_blake3::hash_bytes(&bytes).to_hex();
    let source_hash = schedule.content_hash().to_hex();
    let row = format!(
        concat!(
            r#"{{"schema":"{schema}","kind":"listening-receipt","#,
            r#""gesture_schedule_blake3":"{source_hash}","track_index":0,"#,
            r#""control_rate_hz":{control_rate},"applied_bow_controls":{controls},"#,
            r#""termination":"plate-one-port","frames":{frames},"clipped_samples":{clipped},"#,
            r#""sample_rate_hz":48000,"peak_normalization_factor":0.3,"full_scale_pa":{full_scale_pa},"#,
            r#""wav_blake3":"{hash}","no_claim":"peak-normalized listening aid; first-order stiction and one-way bridge-body fixture, not calibrated or experimentally validated instrument audio"}}"#
        ),
        schema = SCHEMA, source_hash = source_hash, control_rate = schedule.control_rate_hz,
        controls = render.applied_controls().len(), frames = pressure.len(), clipped = clipped,
        full_scale_pa = full_scale_pa, hash = hash,
    );
    // No output files are touched until the complete performance and encoding succeed.
    for dir in [options.wav.parent(), options.receipt.parent()].into_iter().flatten() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
        }
    }
    std::fs::write(&options.wav, bytes).map_err(|e| format!("WAV write: {e}"))?;
    let mut sink = std::fs::File::create(&options.receipt).map_err(|e| format!("receipt create: {e}"))?;
    writeln!(sink, "{row}").map_err(|e| format!("receipt write: {e}"))?;
    eprintln!("bowed_listen: {} frames, {} bow controls -> {}; receipt {}",
        pressure.len(), render.applied_controls().len(), options.wav.display(), options.receipt.display());
    Ok(())
}

fn main() -> std::process::ExitCode {
    match parse_args(std::env::args().skip(1)).and_then(run) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("bowed_listen: refusal: {error}");
            std::process::ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options(args: &[&str]) -> Result<Options, String> {
        parse_args(args.iter().map(|s| s.to_string()))
    }

    #[test]
    fn legacy_invocation_and_explicit_performance_sizes_are_admitted() {
        let old = options(&["--wav", "out.wav", "--receipt", "out.jsonl"]).unwrap();
        assert_eq!((old.samples, old.block, old.gesture), (24000, 256, None));
        let new = options(&["--wav", "out.wav", "--receipt", "out.jsonl",
            "--gesture", "bow.gesture", "--track", "bow", "--samples", "12345", "--block", "37"]).unwrap();
        assert_eq!((new.samples, new.block), (12345, 37));
        for (flag, value) in [("--samples", "0"), ("--samples", "2880001"),
            ("--block", "0"), ("--block", "8193"), ("--samples", "1.5")] {
            assert!(options(&["--wav", "a", "--receipt", "b", flag, value]).is_err());
        }
        assert!(options(&["--wav", "a", "--receipt", "a"]).is_err());
        assert!(options(&["--wav", "a", "--receipt", "b", "--gesture", "a"]).is_err());
    }

    #[test]
    fn shipped_performance_contains_physical_reversal_release_and_reentry() {
        let source = decode_schedule(include_bytes!("../../examples/bow-phrases.gesture")).unwrap();
        assert_eq!(source.control_rate_hz, 700);
        for (tick, negative_speed, lifted) in [(140, true, false), (210, true, true),
            (259, false, false), (329, false, true)] {
            let GestureValue::Bow { velocity_m_per_s, normal_force_n, .. } =
                source.sample_value("bow", tick).unwrap() else { panic!("bow") };
            assert_eq!(velocity_m_per_s < 0.0, negative_speed);
            assert_eq!(normal_force_n == 0.0, lifted);
        }
        let mut invalid = source.to_canonical_bytes();
        invalid.resize(MAX_GESTURE_BYTES + 1, b' ');
        assert!(decode_schedule(&invalid).is_err());
        assert!(decode_schedule(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t999999999\n").is_err());
        assert!(decode_schedule(b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\nevents\t999999999\n").is_err());
    }
}
