//! Music-lane render CLI (bead `frankensim-music-t-out-render-ib15w`,
//! program root `frankensim-music-v8-root-3ez8g`): the first binary that
//! puts music-lane audio ON DISK.
//!
//! ```text
//! music_render <fixture> <out.wav> [--seconds S] [--block N] [--full-scale-pa P] [--schedule FILE]
//! Reed runs also accept `--temperature-k K` and `--ambient-pressure-pa PA`.
//! `--relative-humidity RH` selects the shared moist-air mixture (fraction
//! in [0, 1], default 0). Its saturation-fit and vapor-fraction limits apply;
//! condensation and moisture uptake by solids are not modeled.
//! `--duct-length-m M` and `--duct-radius-m M` set its cylindrical air column.
//! `--duct-outlet-radius-m M` optionally makes that column a linear cone;
//! `--duct-radius-m` then names the inlet radius.
//! These set the shared gas state, not a pitch shift or output gain.
//! They do not model structural heating or phase change of the reed or duct.
//! GasState's declared validity limits apply; a gas-model admission is not proof
//! that the fixed solid specimen remains valid at that ambient temperature.
//! ```
//!
//! `music_render modal INPUT.performance OUT.wav [--block N]` renders supplied
//! mass-normalized structural images and physical actuator histories instead of
//! either fixture. The file owns sample count, physical scale and numerical
//! budgets; the command requires its explicit 48 kHz clock, without resampling.
//! See `examples/MODAL_PERFORMANCES.md` for the format and scope.
//!
//! `music_render plate INPUT.performance OUT.wav [--block N]` derives modes
//! from a supplied flat triangle mesh, per-triangle material/thickness sections,
//! supports and force footprint. See `examples/PLATE_PERFORMANCES.md`.
//!
//! `--schedule performance.gesture` loads the existing canonical
//! `GestureSchedule` format. The reed fixture accepts exactly one explicitly
//! typed blowing-pressure track, bound to voice zero. A new command interrupts
//! a ramp at its current value. Unsupported targets refuse; no controls are
//! silently discarded. The source
//! clock is independent of `--block`, and the schedule digest is retained in
//! the provenance sidecar. This is not an assembly loader or a MIDI renderer.
//!
//! Fixtures are PINNED compositions of gated machinery (`reed`: the
//! massless-reed 2.2 mm characteristic-line voice; `string`: a plucked
//! three-mode exact-ZOH modal string). Rendering goes through the block
//! render API (`fs_couple::render`) — the same path the budget lane
//! measures — and encoding goes through the ONE pascals→PCM owner,
//! `fs_couple::pcm_wav::encode_pcm16_wav`: mono PCM16, physically scaled
//! by the declared full-scale, NEVER peak-normalized (normalization would
//! hide a material or temperature change), clips COUNTED and reported,
//! never rewritten. WAV output uses the incremental encoder with block-sized
//! pressure/PCM staging, and hashes the finalized file through its open handle.
//! The fixture's baked physics storage is unchanged; this is offline file I/O,
//! not a device callback or a whole-program constant-memory claim.
//!
//! Seam decision (recorded here and on beads ib15w + h7xu5.7.8): the
//! music lane's encoder owner is `fs_couple::pcm_wav`; the cinematic
//! stack's receipt-hashed stereo encoder stays cinematic; the 7.8 adapter
//! consumes music-lane pressure through its own admission contract when
//! it lands. No third RIFF writer exists in the music lane.
//!
//! Output contract (campaign-script discipline): the output path and its
//! `.provenance.json` sidecar are REFUSED if they already exist; every
//! emitted JSON-lines field is deterministic — same arguments produce
//! bit-identical WAV bytes and provenance on any single host (the WAV
//! content hash is the replay check). Sample rate is pinned at 48 kHz to
//! keep the ecosystem coherent (fs-psycho refuses other rates).

#[path = "music_render/modal_input.rs"]
mod modal_input;
#[path = "music_render/plate_input.rs"]
mod plate_input;
#[path = "music_render/stream_output.rs"]
mod stream_output;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::{
    PressureGestureBinding, ScheduledControl, ScheduledRenderer, compile_pressure_gestures,
};
use fs_couple::render::{ModalStringVoice, ReedBoreVoice, RenderContext, RenderVoice};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
#[cfg(test)]
use fs_material::gas::GasSpec;
use fs_material::gas::GasState;
use fs_scenario::BeatingReed;
use fs_scenario::gesture::GestureSchedule;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const RATE: u32 = 48_000;
const WAV_HASH_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";

// Input and compilation budgets are separate: a small, densely scheduled file
// must not turn into an unbounded control-clock scan or callback log.
const MAX_SCHEDULE_BYTES: usize = 1 << 20;
const MAX_SCHEDULE_ITEMS: usize = 16_384;
const MAX_SCHEDULE_WORK: u64 = 1_000_000;

fn json_string(value: &str) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => {
                write!(&mut out, "\\u{:04x}", u32::from(ch)).expect("writing to String");
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn fail(what: &str) -> ! {
    // Paths, track ids and nested errors may contain quotes or newlines.
    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"refused\",\"what\":{}}}",
        json_string(what)
    );
    std::process::exit(1)
}

fn decode_schedule(bytes: &[u8]) -> Result<GestureSchedule, String> {
    if bytes.len() > MAX_SCHEDULE_BYTES {
        return Err("schedule exceeds the 1 MiB input budget".to_string());
    }
    let text = core::str::from_utf8(bytes).map_err(|_| "schedule is not UTF-8".to_string())?;
    let lines = text.lines().count();
    // The existing decoder reserves from the declared counts. Bound EVERY
    // count before calling it, including counts in malformed/trailing records.
    for line in text.lines() {
        if let Some(count) = line
            .strip_prefix("tracks\t")
            .or_else(|| line.strip_prefix("events\t"))
        {
            let count = count
                .parse::<usize>()
                .map_err(|_| "invalid schedule item count".to_string())?;
            if count > MAX_SCHEDULE_ITEMS || count > lines {
                return Err("schedule item count exceeds the input/16384-item budget".to_string());
            }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(bytes).map_err(|e| e.to_string())?;
    // Reject ignored suffixes, extra fields and lossy decoder aliases. The
    // accepted file is exactly the artifact whose digest will be published.
    if schedule.to_canonical_bytes() != bytes {
        return Err(
            "schedule must be canonical bytes without trailing or ignored fields".to_string(),
        );
    }
    Ok(schedule)
}

fn load_schedule(path: &Path) -> Result<GestureSchedule, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("schedule open failed: {e}"))?;
    let mut bytes = Vec::new();
    file.take((MAX_SCHEDULE_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("schedule read failed: {e}"))?;
    decode_schedule(&bytes)
}

fn scheduled_controls(
    fixture: &str,
    path: Option<&Path>,
    samples: usize,
) -> Result<(Vec<ScheduledControl>, String), String> {
    let Some(path) = path else {
        return Ok((Vec::new(), String::new()));
    };
    if fixture != "reed" {
        return Err(
            "--schedule currently requires the reed fixture and a blowing-pressure track"
                .to_string(),
        );
    }
    let schedule = load_schedule(path)?;
    let [track] = schedule.tracks() else {
        return Err(
            "the reed fixture requires exactly one gesture track bound to voice zero".to_string(),
        );
    };
    let bindings = [PressureGestureBinding {
        track: track.id.clone(),
        voice: 0,
    }];
    let controls = compile_pressure_gestures(
        &schedule,
        &bindings,
        RATE,
        samples as u64,
        MAX_SCHEDULE_WORK,
    )
    .map_err(|e| e.to_string())?;
    // Source paths are deliberately absent: relocation cannot change replay
    // identity. v1 fixture-only sidecars remain byte-for-byte unchanged.
    #[allow(clippy::format_collect)]
    let schedule_hash: String = schedule
        .content_hash()
        .0
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    let provenance = format!(
        ",\"gesture_schedule\":{{\"blake3\":\"{}\",\"control_rate_hz\":{},\
         \"track\":{},\"voice\":0,\"compiled_controls\":{},\
         \"clock_policy\":\"ceil-control-tick-to-audio-sample\"}}",
        schedule_hash,
        schedule.control_rate_hz,
        json_string(&track.id),
        controls.len()
    );
    Ok((controls, provenance))
}

// create_new closes the exists-check race without replacing evidence. An I/O
// failure may leave incomplete NEW files; it always refuses and never claims
// a rendered artifact. No existing file is truncated and no path is deleted.
fn create_outputs(out: &Path, sidecar: &Path) -> Result<(std::fs::File, std::fs::File), String> {
    // Read access is needed only for bounded hashing of the finalized WAV.
    let audio_file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(out)
        .map_err(|e| format!("wav create refused: {e}"))?;
    let sidecar_file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(sidecar)
        .map_err(|e| format!("sidecar create refused: {e}"))?;
    Ok((audio_file, sidecar_file))
}

#[cfg(test)]
fn cylindrical_duct(length: f64, radius: f64) -> Result<Duct, &'static str> {
    axial_duct(length, radius, radius)
}

fn axial_duct(length: f64, inlet: f64, outlet: f64) -> Result<Duct, &'static str> {
    if [length, inlet, outlet]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err("duct length and endpoint radii must be finite positive metres");
    }
    let segment = if inlet == outlet {
        Segment::Cylinder {
            radius: inlet,
            length,
        }
    } else {
        Segment::Cone {
            inlet_radius: inlet,
            outlet_radius: outlet,
            length,
        }
    };
    Ok(Duct {
        segments: vec![segment],
    })
}

fn reed_context(samples: usize, block: usize, air: &GasState, duct: &Duct) -> RenderContext {
    let reed = BeatingReed {
        rest_opening_m: 4.0e-4,
        width_m: 0.013,
        closing_pressure_pa: 6_000.0,
        blowing_pressure_pa: 2_800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    };
    let voice = ReedBoreVoice::new(
        duct,
        air,
        reed,
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        RATE,
        samples,
        None,
    )
    .unwrap_or_else(|e| fail(&format!("reed voice refused: {e:?}")));
    RenderContext::new(vec![RenderVoice::ReedBore(voice)], block)
}

fn string_context(block: usize) -> RenderContext {
    // A plucked steel-ish string: three modes of the 0.65 m / 60 N /
    // 0.6 g/m card (the bake-off fixture's card, reused deliberately so
    // artifacts cross-reference).
    let wave_speed = (60.0f64 / 6.0e-4).sqrt();
    let modes = (1..=3)
        .map(|k| ModalAcousticMode {
            angular_frequency_rad_s: f64::from(k) * core::f64::consts::PI * wave_speed / 0.65,
            damping_ratio: 1.0e-3 * f64::from(k),
            pressure_per_modal_velocity: fs_math::c64::C64::new(2.0, 0.0),
        })
        .collect::<Vec<_>>();
    let mut model =
        ModalAcousticTimeModel::try_new(RATE, modes, ModalAcousticTimeBudget::audible_reference())
            .unwrap_or_else(|e| fail(&format!("modal model refused: {e:?}")));
    let pluck: Vec<ModalAcousticState> = (1..=3)
        .map(|k| ModalAcousticState {
            displacement_m_sqrt_kg: 1.0e-3 / f64::from(k),
            velocity_m_sqrt_kg_per_s: 0.0,
        })
        .collect();
    model
        .restore_states(&pluck)
        .unwrap_or_else(|e| fail(&format!("pluck refused: {e:?}")));
    let voice = ModalStringVoice::new(model, vec![0.0; 3])
        .unwrap_or_else(|e| fail(&format!("string voice refused: {e:?}")));
    RenderContext::new(vec![RenderVoice::ModalString(voice)], block)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().is_some_and(|arg| arg == "plate") {
        plate_input::run(&args[1..]).unwrap_or_else(|e| fail(&e));
        return;
    }
    if args.first().is_some_and(|arg| arg == "modal") {
        modal_input::run(&args[1..]).unwrap_or_else(|e| fail(&e));
        return;
    }
    let mut positional = Vec::new();
    let mut seconds = 1.0f64;
    let mut block = 512usize;
    let mut full_scale_pa = 200.0f64;
    let mut schedule_path: Option<PathBuf> = None;
    let mut temperature_k = 293.15;
    let mut ambient_pressure_pa = 101_325.0;
    let mut relative_humidity = 0.0;
    let mut ambient_explicit = false;
    let mut duct_length_m = 0.50;
    let mut duct_radius_m = 0.0022;
    let mut duct_outlet_radius_m = None;
    let mut geometry_explicit = false;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--relative-humidity" => {
                relative_humidity = iter
                    .next()
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or_else(|| fail("--relative-humidity requires a fraction in [0, 1]"));
                ambient_explicit = true;
            }
            "--duct-length-m" | "--duct-radius-m" | "--duct-outlet-radius-m" => {
                let value = iter
                    .next()
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or_else(|| fail("duct dimensions require a number in metres"));
                if arg == "--duct-length-m" {
                    duct_length_m = value;
                } else if arg == "--duct-radius-m" {
                    duct_radius_m = value;
                } else {
                    duct_outlet_radius_m = Some(value);
                }
                geometry_explicit = true;
            }
            "--temperature-k" | "--ambient-pressure-pa" => {
                let value = iter
                    .next()
                    .and_then(|v| v.parse::<f64>().ok())
                    .unwrap_or_else(|| {
                        fail("ambient options require a number in Kelvin or Pascal")
                    });
                if arg == "--temperature-k" {
                    temperature_k = value;
                } else {
                    ambient_pressure_pa = value;
                }
                ambient_explicit = true;
            }
            "--seconds" => {
                seconds = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--seconds needs a positive number"));
            }
            "--block" => {
                block = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--block needs a positive integer"));
            }
            "--full-scale-pa" => {
                full_scale_pa = iter
                    .next()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| fail("--full-scale-pa needs a positive number"));
            }
            "--schedule" => {
                if schedule_path.is_some() {
                    fail("--schedule may only be specified once");
                }
                let path = iter
                    .next()
                    .unwrap_or_else(|| fail("--schedule needs a file path"));
                schedule_path = Some(PathBuf::from(path.as_str()));
            }
            other if other.starts_with('-') => fail(&format!("unknown option: {other}")),
            other => positional.push(other.to_string()),
        }
    }
    let [fixture, out_path] = positional.as_slice() else {
        fail(
            "usage: music_render <reed|string> <out.wav> [--seconds S] [--block N] [--full-scale-pa P] [--schedule FILE] [--temperature-k K] [--ambient-pressure-pa PA] [--relative-humidity RH] [--duct-length-m M] [--duct-radius-m M] [--duct-outlet-radius-m M] (ambient and duct options require reed)",
        );
    };
    if !(seconds > 0.0 && seconds <= 600.0) {
        fail("--seconds must be in (0, 600]");
    }
    if block == 0 || block > 1 << 16 {
        fail("--block must be in 1..=65536");
    }
    if !full_scale_pa.is_finite() || full_scale_pa <= 0.0 {
        fail("--full-scale-pa must be positive and finite");
    }
    let out = Path::new(out_path);
    let sidecar = out.with_extension("provenance.json");
    if out == sidecar.as_path() {
        fail("output and sidecar must have distinct file paths");
    }
    if out.exists() || sidecar.exists() {
        fail("output or sidecar already exists; this lane refuses to overwrite evidence");
    }

    let samples = (seconds * f64::from(RATE)).round() as usize;
    if samples == 0 {
        fail("--seconds must round to at least one audio sample");
    }
    let (controls, schedule_provenance) =
        scheduled_controls(fixture, schedule_path.as_deref(), samples).unwrap_or_else(|e| fail(&e));
    if ambient_explicit && fixture != "reed" {
        fail("ambient gas options require reed; the string fixture has no gas coupling");
    }
    if geometry_explicit && fixture != "reed" {
        fail("duct geometry options require reed; the string fixture has no duct");
    }
    let mut ambient_provenance = String::new();
    let context = match fixture.as_str() {
        "reed" => {
            let outlet = duct_outlet_radius_m.unwrap_or(duct_radius_m);
            let duct = axial_duct(duct_length_m, duct_radius_m, outlet).unwrap_or_else(|e| fail(e));
            let air =
                GasState::try_new_moist_air(temperature_k, ambient_pressure_pa, relative_humidity)
                    .unwrap_or_else(|e| fail(&format!("ambient gas state refused: {e:?}")));
            let gas_spec = if relative_humidity == 0.0 {
                "dry_air_ussa1976"
            } else {
                "moist_air_ussa1976_water_vapor_nist"
            };
            let geometry = if outlet == duct_radius_m {
                format!(
                    "{{\"shape\":\"cylinder\",\"length_m\":{duct_length_m:e},\"radius_m\":{duct_radius_m:e}}}"
                )
            } else {
                format!(
                    "{{\"shape\":\"cone\",\"length_m\":{duct_length_m:e},\"inlet_radius_m\":{duct_radius_m:e},\"outlet_radius_m\":{outlet:e}}}"
                )
            };
            ambient_provenance = format!(
                ",\"ambient_gas\":{{\"spec\":\"{gas_spec}\",\"temperature_k\":{temperature_k:e},\"pressure_pa\":{ambient_pressure_pa:e},\"relative_humidity\":{relative_humidity:e}}},\"duct\":{geometry}"
            );
            reed_context(samples, block, &air, &duct)
        }
        "string" => string_context(block),
        _ => fail("fixture must be `reed` or `string`"),
    };

    let mut renderer = ScheduledRenderer::new(context, controls, MAX_SCHEDULE_WORK as usize)
        .unwrap_or_else(|e| fail(&format!("schedule admission refused: {e}")));
    // Reserve both paths before advancing the physical performance. An error
    // leaves only new incomplete files, never a success sidecar or overwritten
    // evidence. Waveforms are streamed; no full pressure/PCM history is staged.
    let (mut audio_file, mut sidecar_file) =
        create_outputs(out, &sidecar).unwrap_or_else(|e| fail(&e));
    let rendered = stream_output::render_waveform(
        &mut renderer,
        &mut audio_file,
        samples,
        block,
        full_scale_pa,
    )
    .unwrap_or_else(|e| fail(&e));
    let clipped = rendered.clipped;
    let peak = rendered.peak_pa;
    let rms = rendered.rms_pa;
    let hash_hex = rendered.hash.to_hex();

    // Deterministic provenance sidecar: everything a replayer needs. No
    // wall-clock, no commit stamp (the git history of committed artifacts
    // carries those); the WAV content hash is the replay check.
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"{fixture}\",\
         \"sample_rate_hz\":{RATE},\"samples\":{samples},\"block\":{block},\
         \"full_scale_pa\":{full_scale_pa:e},\"clipped_samples\":{clipped},\
         \"peak_pa\":{peak:e},\"rms_pa\":{rms:e},\"wav_blake3\":\"{hash_hex}\",\
         \"encoder\":\"fs_couple::pcm_wav (mono PCM16, never peak-normalized)\"{schedule_provenance}{ambient_provenance}}}"
    );
    writeln!(sidecar_file, "{provenance}")
        .and_then(|_| sidecar_file.flush())
        .unwrap_or_else(|e| fail(&format!("sidecar write failed: {e}")));

    println!(
        "{{\"suite\":\"music-render\",\"verdict\":\"rendered\",\"fixture\":\"{fixture}\",\
         \"wav\":{},\"samples\":{samples},\"clipped\":{clipped},\"peak_pa\":{peak:.3},\
         \"rms_pa\":{rms:.3},\"wav_blake3\":\"{hash_hex}\"}}",
        json_string(&out.display().to_string())
    );
}

#[cfg(test)]
mod schedule_input_tests {
    use super::*;
    use fs_scenario::gesture::{GestureTarget, GestureTrack, GestureValue};

    #[test]
    fn tapered_ducts_render_from_endpoint_geometry() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(axial_duct(0.5, 0.002, bad).is_err());
        }
        assert!(matches!(
            axial_duct(0.5, 0.002, 0.002).unwrap().segments[0],
            Segment::Cylinder { .. }
        ));
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let render = |outlet, block| {
            let duct = axial_duct(0.5, 0.002, outlet).unwrap();
            assert_eq!(duct.segments[0].inlet_radius(), 0.002);
            assert_eq!(duct.segments[0].outlet_radius(), outlet);
            let mut context = reed_context(4096, block, &gas, &duct);
            let mut output = vec![0.0; 4096];
            for chunk in output.chunks_mut(block) {
                context.block(chunk).unwrap();
            }
            assert!(output.iter().all(|p| p.is_finite()));
            assert!(output.iter().any(|p| p.abs() > 1.0));
            output
        };
        let cylinder = render(0.002, 37);
        for outlet in [0.0015, 0.0022] {
            let cone = render(outlet, 37);
            assert!(
                cone.iter()
                    .zip(&cylinder)
                    .any(|(a, b)| a.to_bits() != b.to_bits())
            );
            assert!(
                cone.iter()
                    .zip(render(outlet, 512))
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
    }

    #[test]
    fn declared_duct_dimensions_drive_audio_and_refuse_invalid_geometry() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(cylindrical_duct(bad, 0.0022).is_err());
            assert!(cylindrical_duct(0.5, bad).is_err());
        }
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let render = |length, radius, block| {
            let duct = cylindrical_duct(length, radius).unwrap();
            let mut context = reed_context(4096, block, &gas, &duct);
            let mut output = vec![0.0; 4096];
            for chunk in output.chunks_mut(block) {
                context.block(chunk).unwrap();
            }
            assert!(output.iter().all(|x| x.is_finite()));
            assert!(output.iter().any(|x| x.abs() > 1.0));
            output
        };
        let baseline = render(0.5, 0.0022, 37);
        for (length, radius) in [(0.6, 0.0022), (0.5, 0.0020)] {
            let changed = render(length, radius, 37);
            assert!(
                changed
                    .iter()
                    .zip(&baseline)
                    .any(|(a, b)| a.to_bits() != b.to_bits())
            );
            let other_partition = render(length, radius, 512);
            assert!(
                changed
                    .iter()
                    .zip(other_partition)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
            );
        }
    }

    #[test]
    fn ambient_gas_changes_physical_reed_audio_without_callback_dependence() {
        let render = |temperature, pressure, block| {
            let gas =
                GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
            let duct = cylindrical_duct(0.50, 0.0022).unwrap();
            let mut context = reed_context(4096, block, &gas, &duct);
            let mut output = vec![0.0; 4096];
            for chunk in output.chunks_mut(block) {
                context.block(chunk).unwrap();
            }
            assert!(output.iter().all(|x| x.is_finite()));
            assert!(output.iter().any(|x| x.abs() > 1.0));
            output
        };
        let baseline = render(293.15, 101_325.0, 37);
        for (temperature, pressure) in [(313.15, 101_325.0), (293.15, 80_000.0)] {
            let changed = render(temperature, pressure, 37);
            assert!(
                changed
                    .iter()
                    .zip(&baseline)
                    .any(|(a, b)| a.to_bits() != b.to_bits())
            );
            assert_eq!(changed, render(temperature, pressure, 512));
        }
    }

    fn canonical_schedule() -> Vec<u8> {
        GestureSchedule::try_new(
            200,
            vec![GestureTrack {
                id: "blow".to_string(),
                target: GestureTarget::BlowingPressure,
                initial: GestureValue::PressurePa(2800.0),
                events: Vec::new(),
            }],
        )
        .unwrap()
        .to_canonical_bytes()
    }

    #[test]
    fn canonical_input_round_trips_but_ignored_suffixes_refuse() {
        let bytes = canonical_schedule();
        let schedule = decode_schedule(&bytes).unwrap();
        assert_eq!(schedule.to_canonical_bytes(), bytes);
        let mut trailing = bytes.clone();
        trailing.extend_from_slice(b"ignored\n");
        assert!(
            decode_schedule(&trailing)
                .unwrap_err()
                .contains("canonical")
        );
        let extra = String::from_utf8(bytes)
            .unwrap()
            .replace("blowing-pressure\n", "blowing-pressure\textra\n");
        assert!(decode_schedule(extra.as_bytes()).is_err());
    }

    #[test]
    fn untrusted_declared_sizes_are_rejected_before_decoder_reservation() {
        for bytes in [
            b"frankensim-gesture-schedule-v1\ncontrol_rate_hz\t200\ntracks\t18446744073709551615\n"
                .as_slice(),
            b"events\t9999999999999999999999999999999999999999\n".as_slice(),
            b"tracks\t16385\n".as_slice(),
            b"events\t100\n".as_slice(),
        ] {
            assert!(decode_schedule(bytes).unwrap_err().contains("count"));
        }
        assert!(decode_schedule(&vec![b'x'; MAX_SCHEDULE_BYTES + 1]).is_err());
        assert!(decode_schedule(&[0xff]).is_err());
    }

    #[test]
    fn structured_diagnostics_escape_untrusted_strings() {
        assert_eq!(json_string("a\"b"), "\"a\\\"b\"");
        assert_eq!(json_string("a\\b"), "\"a\\\\b\"");
        assert_eq!(json_string("\n\r\t"), "\"\\n\\r\\t\"");
        assert_eq!(json_string("\0"), "\"\\u0000\"");
        assert_eq!(json_string("é"), "\"é\"");
    }

    #[test]
    fn no_schedule_does_not_add_provenance_or_controls() {
        let (controls, provenance) = scheduled_controls("string", None, 480).unwrap();
        assert!(controls.is_empty());
        assert!(provenance.is_empty());
    }
}
