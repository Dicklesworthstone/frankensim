//! Real file-driven contact and plate physics through the ensemble command.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::render::plate::file::PlatePerformance;
use fs_couple::render::schedule::force::file::ModalPerformance;
use fs_material::gas::GasState;
use fs_couple::brass_loop::performance::LipDescription;
use fs_couple::voice_performance::VoiceDescription;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_material::gas_species::input::ThermochemicalSources;

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
#[allow(dead_code)]
mod gas_sources;

const STRIKER: &str = include_str!("../examples/free-striker-192k.performance");
const PLATE: &str = include_str!("../examples/plate-mesh.performance");
const REED_OUTLET: &str = include_str!("../examples/reed-baffled-mixture.performance");
static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-ensemble-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap(); path
}
fn run(output: &Path, inputs: &[(&str, &Path)], block: usize, decimate: bool) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
    command.arg("ensemble").arg(output).args(["--full-scale-pa", "1", "--block", &block.to_string()]);
    for &(kind, path) in inputs { command.arg(kind).arg(path); }
    if decimate { command.arg("--decimate"); }
    command.output().unwrap()
}
fn succeeded(output: Output) {
    assert!(output.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
}
fn reference(striker: &str, plate: &str) -> Vec<f64> {
    let performance = ModalPerformance::from_bytes(striker.as_bytes(), 37).unwrap();
    let mut input = vec![0.0; performance.info().samples as usize];
    let mut renderer = performance.into_renderer();
    for block in input.chunks_mut(37) { renderer.block(block).unwrap(); }
    let mut filter = Decimator::new(4, 1).unwrap();
    let delay = filter.delay_output_frames() as usize;
    let mut pressure: Vec<_> = input.chunks_exact(4).map(|group| {
        let p = filter.preview(group).unwrap()[0]; filter.commit(); p
    }).collect();
    let plate = PlatePerformance::from_bytes(plate.as_bytes(), 37).unwrap();
    assert_eq!(plate.info().samples as usize, pressure.len());
    let mut plate_pressure = vec![0.0; pressure.len()];
    let mut renderer = plate.into_renderer();
    for block in plate_pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    for (p, q) in pressure[delay..].iter_mut().zip(plate_pressure) { *p += q; }
    pressure
}

#[test]
fn command_mixes_actual_contact_and_geometry_with_aligned_pcm_across_partitions() {
    let dir = directory();
    let striker = dir.join("striker.performance"); let plate = dir.join("plate.performance");
    std::fs::write(&striker, STRIKER).unwrap(); std::fs::write(&plate, PLATE).unwrap();
    let pressure = reference(STRIKER, PLATE);
    assert!(pressure.iter().any(|p| p.abs() > 0.01));
    let (expected, clips) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    let hash = fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1", &expected).to_hex();
    for block in [1, 37, 512] {
        let out = dir.join(format!("block-{block}.wav"));
        succeeded(run(&out, &[("--modal", &striker), ("--plate", &plate)], block, true));
        assert_eq!(std::fs::read(&out).unwrap(), expected);
        let sidecar = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        for token in ["\"samples\":4801", "\"mechanics_samples\":19204", "\"ratio\":4",
            "\"common_delay_output_samples\":44", "\"alignment_delay_samples\":44",
            "\"alignment_delay_samples\":0", "\"source_pcm_scale_applied\":false",
            "\"delay_compensated\":false", "\"tail\":\"no-flush-declared-window\""] {
            assert!(sidecar.contains(token), "{token}: {sidecar}");
        }
        assert!(sidecar.contains(&hash));
        assert!(sidecar.contains(&format!("\"clipped_samples\":{clips}")));
    }
    let relocated = dir.join("relocated.performance"); std::fs::write(&relocated, STRIKER).unwrap();
    let replay = dir.join("replay.wav");
    succeeded(run(&replay, &[("--modal", &relocated), ("--plate", &plate)], 37, true));
    assert_eq!(std::fs::read(replay.with_extension("provenance.json")).unwrap(),
        std::fs::read(dir.join("block-37.provenance.json")).unwrap());
    assert!(!run(&replay, &[("--modal", &striker), ("--plate", &plate)], 37, true).status.success());
    assert_eq!(std::fs::read(&replay).unwrap(), expected);
}

#[test]
fn a_material_change_reaches_the_mixed_waveform_not_just_its_source_hash() {
    let dir = directory(); let striker = dir.join("striker.performance");
    std::fs::write(&striker, STRIKER).unwrap();
    let mut waves = Vec::new();
    for (i, text) in [PLATE.to_string(), PLATE.replace("11000000000.0", "22000000000.0")].iter().enumerate() {
        let plate = dir.join(format!("plate-{i}.performance")); std::fs::write(&plate, text).unwrap();
        let out = dir.join(format!("material-{i}.wav"));
        succeeded(run(&out, &[("--modal", &striker), ("--plate", &plate)], 37, true));
        let (expected, _) = encode_pcm16_wav(&reference(STRIKER, text), 48_000, 1.0).unwrap();
        let actual = std::fs::read(&out).unwrap(); assert_eq!(actual, expected); waves.push(actual);
    }
    assert_ne!(&waves[0][44..], &waves[1][44..]);
}

fn modal(rate: u32, samples: u64, scale: f64) -> String {
    format!("frankensim-modal-performance-v1\nsample_rate_hz {rate}\nsamples {samples}\nfull_scale_pa {scale}\n\
        limits 0.9 1 1000 1000 1000\ncompile_limits 10 20\nvoices 1\nvoice retain-state 1 1\n\
        mode 1000 0.01 1 0.1 0.0001 0\nport 1 1\nevents 1\nforce 37 0 0 0\n")
}

#[test]
fn source_pcm_scales_are_not_hidden_gains_on_the_same_clock_path() {
    let dir = directory(); let a = dir.join("a.performance"); let b = dir.join("b.performance");
    std::fs::write(&a, modal(48_000, 257, 1.0)).unwrap();
    let mut waves = Vec::new();
    for (i, scale) in [1.0, 200.0].into_iter().enumerate() {
        std::fs::write(&b, modal(48_000, 257, scale)).unwrap();
        let out = dir.join(format!("scale-{i}.wav"));
        succeeded(run(&out, &[("--modal", &a), ("--modal", &b)], 37, false));
        waves.push(std::fs::read(&out).unwrap());
        let metadata = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(metadata.contains("\"common_delay_output_samples\":0"));
    }
    assert_eq!(waves[0], waves[1]);
    let mut renderer = ModalPerformance::from_bytes(modal(48_000, 257, 1.0).as_bytes(), 257).unwrap().into_renderer();
    let mut p = vec![0.0; 257]; renderer.block(&mut p).unwrap();
    for value in &mut p { *value = *value + *value; }
    assert_eq!(waves[0], encode_pcm16_wav(&p, 48_000, 1.0).unwrap().0);
}

#[test]
fn all_sources_and_clocks_are_admitted_before_either_output_artifact_is_created() {
    let dir = directory(); let a = dir.join("a.performance"); let b = dir.join("b.performance");
    std::fs::write(&a, modal(48_000, 257, 1.0)).unwrap();
    for (i, (source, decimate)) in [
        (modal(48_000, 258, 1.0), false), // unequal physical duration
        (modal(96_000, 515, 1.0), true), // incomplete output interval
        (modal(96_000, 514, 1.0), false), // conversion was not requested
        (modal(100_000, 514, 1.0), true), // noninteger conversion
        (modal(44_100, 257, 1.0), true), // no upsampling
        (modal(48_000, 257, 1.0), true), // redundant conversion request
        ("unrecognized model\n".into(), false),
    ].into_iter().enumerate() {
        std::fs::write(&b, source).unwrap(); let out = dir.join(format!("refused-{i}.wav"));
        assert!(!run(&out, &[("--modal", &a), ("--modal", &b)], 37, decimate).status.success());
        assert!(!out.exists()); assert!(!out.with_extension("provenance.json").exists());
    }
}

#[test]
fn a_part_heard_in_another_gas_refuses_before_either_output_exists() {
    // The exact mix reed_gas_mixture::outlet once summed: one baffled reed heard
    // in its 100 kPa terminal gas, the same reed in the same mixture at 50 kPa.
    let dir = directory(); let (a, b) = (dir.join("100kPa.performance"), dir.join("50kPa.performance"));
    std::fs::write(&a, REED_OUTLET).unwrap();
    assert!(REED_OUTLET.contains("mixture 293.15 100000"));
    std::fs::write(&b, REED_OUTLET.replace("mixture 293.15 100000", "mixture 293.15 50000")).unwrap();
    let out = dir.join("two-gases.wav");
    let refused = run(&out, &[("--reed", &a), ("--reed", &b)], 37, false);
    assert!(!refused.status.success());
    let message = String::from_utf8_lossy(&refused.stdout);
    assert!(message.contains("ensemble part 1 (reed) is heard in a different gas than part 0 (reed)")
        && message.contains("one ensemble has one listener gas"), "{message}");
    assert!(!out.exists()); assert!(!out.with_extension("provenance.json").exists());
    // Either reed alone, or twice, is one listener gas.
    succeeded(run(&dir.join("one-gas.wav"), &[("--reed", &b), ("--reed", &b)], 37, false));
}

#[test]
fn the_receipt_names_the_one_listener_gas_and_each_unconstrained_part() {
    let dir = directory(); let striker = dir.join("striker.performance");
    let dry = dir.join("dry.performance"); let thin = dir.join("thin.performance");
    std::fs::write(&striker, STRIKER).unwrap(); std::fs::write(&dry, PLATE).unwrap();
    let row = "ambient 293.15 101325 0\n";
    assert_eq!(PLATE.matches(row).count(), 1);
    std::fs::write(&thin, PLATE.replace(row, "ambient 293.15 70927.5 0\n")).unwrap();
    let medium = |gas: GasState, first: usize, rest: &str| format!(
        "\"listener_medium\":{{\"density_kg_m3\":{:e},\"sound_speed_m_s\":{:e},\"set_by_part\":{first},\"unconstrained_parts\":[{rest}]}}",
        gas.density, gas.sound_speed);
    let receipt = |out: &Path| std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
    // One gas, bit for bit, set by the first part.
    let same = dir.join("one-gas.wav");
    succeeded(run(&same, &[("--plate", &dry), ("--plate", &dry)], 37, false));
    let air = GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap();
    assert!(receipt(&same).contains(&medium(air, 0, "")), "{}", receipt(&same));
    // An authored modal transfer does not constrain the listener, and says why.
    let mixed = dir.join("unconstrained.wav");
    succeeded(run(&mixed, &[("--modal", &striker), ("--plate", &thin)], 37, true));
    let thin_air = GasState::try_new_moist_air(293.15, 70_927.5, 0.0).unwrap();
    let modal_part = "{\"part\":0,\"kind\":\"modal\",\"reason\":\"authored transfer; no radiation block\"}";
    assert!(receipt(&mixed).contains(&medium(thin_air, 1, modal_part)), "{}", receipt(&mixed));
    // No part declares a listener gas.
    let a = dir.join("a.performance"); std::fs::write(&a, modal(48_000, 257, 1.0)).unwrap();
    let none = dir.join("none.wav");
    succeeded(run(&none, &[("--modal", &a), ("--modal", &a)], 37, false));
    assert!(receipt(&none).contains(&format!("\"listener_medium\":{{\"declared\":\"none\",\"unconstrained_parts\":[{modal_part},{}]}}",
        modal_part.replace("\"part\":0", "\"part\":1"))), "{}", receipt(&none));
}

// Synthetic geometry/excitation fixtures: the two existing source owners
// determine all flow, mechanics, propagation and pressure, with no timbre preset.
fn lip_source(gas: &str) -> String {
    format!("frankensim-lip-performance-v1\n\
        audio 48000 1200\nfull_scale_pa 30000\nambient {gas}\n\
        lip_response 200 7 0.1 0.012 0.0005 0.0001\n\
        lip_source test fixture, not a measurement\n\
        bore_source test fixture, not a measurement\n\
        cup_volume_m3 1.5e-6\ntermination flanged\nmodes 3 1\n\
        cylinder 0.006 0.30\ncone 0.006 0.012 0.25\n\
        microphone 0 0 0.3 4 16 3000\nschedule\n\
        frankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\ntracks\t1\n\
        track\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t1\n\
        event\t0e0\t1e-2\tpressure-pa\t2.8e3\n")
}

fn voice_source(gas: &str) -> String {
    let mut text = format!("frankensim-voice-performance-v1\n\
        audio 48000 1200 2\nfull_scale_pa 2\nambient {gas}\n\
        fold two_mass 1.25e-4 2.5e-5 128 12.8 40 0.1 3.57e-4 1.4e-2 2.5e-3\n\
        termination flanged\ntract_source test-fixture\ntract_licence test fixture, not data\n\
        microphone 0 0 0.3 4 16 4000\nsubglottal 0 0\nsubglottal 0.02 900\n");
    for area in [0.5e-4, 0.6e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4] {
        text.push_str(&format!("section {area:e} 0.02\n"));
    }
    text
}

fn lip_voice_reference(
    lip: &str,
    voice: &str,
    sources: Option<&ThermochemicalSources<'_>>,
) -> (Vec<f64>, f64) {
    let mut lip = LipDescription::from_bytes_with_gas_sources(lip.as_bytes(), sources).unwrap().build().unwrap();
    let mut voice = VoiceDescription::from_bytes_with_gas_sources(voice.as_bytes(), sources).unwrap().build().unwrap();
    let delays = [lip.filter_delay_samples(), voice.filter_delay_samples()];
    let common = delays[0].max(delays[1]);
    let mut a = vec![0.0; 1200]; let mut b = a.clone();
    for block in a.chunks_mut(37) { lip.block(block).unwrap(); }
    for block in b.chunks_mut(37) { voice.block(block).unwrap(); }
    assert!(a.iter().any(|p| p.abs() > 1e-12));
    assert!(b.iter().any(|p| p.abs() > 1e-12));
    let mut pressure = vec![0.0; a.len()];
    for (source, delay) in [a, b].iter().zip(delays) {
        let extra = common - delay;
        assert_eq!(extra.fract(), 0.0);
        for (p, q) in pressure.iter_mut().skip(extra as usize).zip(source) { *p += q; }
    }
    (pressure, common)
}

#[test]
fn lip_and_voice_keep_their_material_medium_substeps_and_original_microphone_latency() {
    let dir = directory(); let lip = dir.join("lip.performance"); let voice = dir.join("voice.performance");
    let dry = dir.join("modal.performance");
    let gas = "mixture 300 100000 mole 2 N2 0.8 O2 0.2";
    let a = lip_source(gas); let b = voice_source(gas);
    std::fs::write(&lip, &a).unwrap(); std::fs::write(&voice, &b).unwrap();
    let (mut pressure, delay) = lip_voice_reference(&a, &b, None);
    assert!(delay > 0.0, "the differentiating microphones already contribute latency");
    // This third physical source has no microphone filter. Its pressure must
    // receive the other owners' existing latency exactly once before summing.
    let dry_text = modal(48_000, 1200, 1.0);
    std::fs::write(&dry, &dry_text).unwrap();
    let mut renderer = ModalPerformance::from_bytes(dry_text.as_bytes(), 37).unwrap().into_renderer();
    let mut dry_pressure = vec![0.0; 1200];
    for block in dry_pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    for (p, q) in pressure.iter_mut().skip(delay as usize).zip(dry_pressure) { *p += q; }
    let expected = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap().0;
    for block in [37, 113] {
        let out = dir.join(format!("lip-voice-{block}.wav"));
        succeeded(run(&out, &[("--lip", &lip), ("--voice", &voice), ("--modal", &dry)], block, false));
        assert_eq!(std::fs::read(&out).unwrap(), expected);
        let receipt = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        for field in ["\"kind\":\"lip\"", "\"kind\":\"voice\"", "\"mechanics_sample_rate_hz\":96000",
            "\"mechanics_samples\":2400", "\"pressure_samples\":1200", "\"part\":2,\"kind\":\"modal\"",
            "\"modes\":4,\"compiled_controls\":2", "\"modes\":2,\"compiled_controls\":2",
            "original multirate baffled exterior receiver"] {
            assert!(receipt.contains(field), "{field}: {receipt}");
        }
        assert!(receipt.contains(&format!("\"common_delay_output_samples\":{delay}")));
    }
    // A common listener cannot inhabit two pressure/composition states. Refuse
    // both conflicts before publishing a WAV or success receipt.
    for (i, other) in ["mixture 300 50000 mole 2 N2 0.8 O2 0.2",
        "mixture 300 100000 mole 2 N2 0.6 CO2 0.4"].into_iter().enumerate() {
        std::fs::write(&voice, voice_source(other)).unwrap();
        let out = dir.join(format!("lip-voice-conflict-{i}.wav"));
        let result = run(&out, &[("--lip", &lip), ("--voice", &voice)], 37, false);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stdout).contains("one ensemble has one listener gas"));
        assert!(!out.exists() && !out.with_extension("provenance.json").exists());
    }
}

#[test]
fn arbitrary_source_pack_changes_reach_lip_and_voice_ensemble_pressure() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let dir = directory(); let lip = dir.join("lip.performance"); let voice = dir.join("voice.performance");
    let gas = "thermochemical 300 100000 mole 1 ensemble-source-gas 1";
    let a = lip_source(gas); let b = voice_source(gas);
    std::fs::write(&lip, &a).unwrap(); std::fs::write(&voice, &b).unwrap();
    let mut pressures = Vec::new();
    for (index, molar_mass) in [0.028, 0.042].into_iter().enumerate() {
        let packs = gas_sources::packs("ensemble-source-gas", molar_mass, 3.5, 200.0);
        let bindings = [packs.binding()];
        let registry = ThermochemicalSources { species: &bindings, water_species: None, phase: ThermochemicalPhase::FrozenGas };
        let (pressure, _) = lip_voice_reference(&a, &b, Some(&registry));
        let paths = [dir.join(format!("species-{index}.pack")), dir.join(format!("nasa9-{index}.pack")),
            dir.join(format!("transport-{index}.pack"))];
        for (path, bytes) in paths.iter().zip([packs.species.to_bytes(), packs.nasa9.to_bytes(), packs.transport.to_bytes()]) {
            std::fs::write(path, bytes).unwrap();
        }
        let out = dir.join(format!("sourced-lip-voice-{index}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("ensemble").arg(&out).args(["--full-scale-pa", "1000", "--block", "37"])
            .arg("--lip").arg(&lip).arg("--voice").arg(&voice)
            .arg("--gas-component").args(&paths).args(["--gas-water", "none", "--gas-phase", "frozen"])
            .output().unwrap();
        succeeded(result);
        assert_eq!(std::fs::read(&out).unwrap(), encode_pcm16_wav(&pressure, 48_000, 1000.0).unwrap().0);
        let receipt = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(receipt.contains("ensemble-source-gas"));
        pressures.push(pressure);
    }
    let energy: f64 = pressures[0].iter().map(|p| p*p).sum();
    let difference: f64 = pressures[0].iter().zip(&pressures[1]).map(|(a,b)| (a-b).powi(2)).sum();
    assert!(energy > 0.0 && difference > 1e-8 * energy,
        "changing only source material data must alter the mixed physical pressure");
}
