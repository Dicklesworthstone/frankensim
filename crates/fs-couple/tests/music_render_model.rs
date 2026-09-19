//! Actual-command tests: custom input -> existing physics -> streamed WAV.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_blake3::hash_domain;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::schedule::force::file::MODAL_PERFORMANCE_HASH_DOMAIN;
use fs_math::c64::C64;

const INPUT: &str = include_str!("../examples/modal-ports.performance");
const WAV_DOMAIN: &str = "org.frankensim.fs-couple.music-render-wav.v1";
static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let time = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-modal-input-{}-{time}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path
}
fn invoke(input: &Path, output: &Path, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_music_render")).arg("modal").arg(input).arg(output)
        .args(extra).output().unwrap()
}
fn model(modes: &[(f64, f64, f64)]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(48000, modes.iter().map(|&(omega, damping, transfer)| ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: damping,
        pressure_per_modal_velocity: C64::new(transfer, 0.0),
    }).collect(), ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn direct_wav() -> Vec<u8> {
    // Independent input wiring: do not call the importer or schedule compiler.
    let mut a = model(&[(1000.0, 0.03, 1.0), (2000.0, 0.04, 0.75)]);
    a.initialize_static_equilibrium(&[0.4375, -0.75]).unwrap();
    let mut b = model(&[(1500.0, 0.025, 0.5)]);
    b.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 0.0,
        velocity_m_sqrt_kg_per_s: 0.01 }]).unwrap();
    let pressure: Vec<_> = (0..4801).map(|sample| {
        let p = if sample < 13 { 0.5 } else { 0.0 };
        let q = if sample < 71 { -0.25 } else if sample < 120 { 0.5 } else { 0.0 };
        a.step(&[p + 0.25*q, -0.5*p + 2.0*q]).unwrap().observer_pressure_pa
            + b.step(&[if sample < 2400 { 0.0 } else { 0.125 }]).unwrap().observer_pressure_pa
    }).collect();
    encode_pcm16_wav(&pressure, 48000, 0.02).unwrap().0
}

#[test]
fn command_uses_imported_states_ports_and_short_tail_at_every_partition() {
    let dir = directory();
    let input = dir.join("model.performance");
    std::fs::write(&input, INPUT).unwrap();
    let expected = direct_wav();
    let input_hash = hash_domain(MODAL_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()).to_hex();
    let wav_hash = hash_domain(WAV_DOMAIN, &expected).to_hex();
    for block in ["1", "37", "128"] {
        let output = dir.join(format!("{block}.wav"));
        let result = invoke(&input, &output, &["--block", block]);
        assert!(result.status.success(), "{} {}", String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
        let actual = std::fs::read(&output).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 44 + 2*4801);
        assert!(actual[44..].iter().any(|&x| x != 0));
        let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(meta.contains(&input_hash));
        assert!(meta.contains(&wav_hash));
        assert!(meta.contains("\"voices\":2,\"modes\":3,\"force_events\":4"));
        assert!(meta.contains("authored reduced model; no physical-validation claim"));
    }
}

#[test]
fn input_relocation_replays_but_physical_mode_changes_change_audio() {
    let dir = directory();
    let a = dir.join("a.performance");
    let b = dir.join("relocated.performance");
    std::fs::write(&a, INPUT).unwrap();
    std::fs::write(&b, INPUT).unwrap();
    let out_a = dir.join("a.wav");
    let out_b = dir.join("b.wav");
    assert!(invoke(&a, &out_a, &[]).status.success());
    assert!(invoke(&b, &out_b, &[]).status.success());
    assert_eq!(std::fs::read(&out_a).unwrap(), std::fs::read(&out_b).unwrap());
    assert_eq!(std::fs::read(out_a.with_extension("provenance.json")).unwrap(),
        std::fs::read(out_b.with_extension("provenance.json")).unwrap());
    let changed = dir.join("changed.performance");
    std::fs::write(&changed, INPUT.replace("mode 1500", "mode 1800")).unwrap();
    let out_changed = dir.join("changed.wav");
    assert!(invoke(&changed, &out_changed, &[]).status.success());
    assert_ne!(std::fs::read(&out_a).unwrap(), std::fs::read(&out_changed).unwrap());
}

#[test]
fn invalid_model_and_conflicting_overrides_create_no_output_and_existing_files_survive() {
    let dir = directory();
    for (i, text) in [
        INPUT.replace("sample_rate_hz 48000", "sample_rate_hz 44100"),
        INPUT.replace("mode 1500", "mode 1000000"),
        INPUT.replace("force 2400", "force 4801"),
        INPUT.replace("port 0 1.25", "port NaN 1.25"),
        format!("{INPUT}ignored 1\n"),
    ].into_iter().enumerate() {
        let input = dir.join(format!("bad{i}.performance"));
        let output = dir.join(format!("bad{i}.wav"));
        std::fs::write(&input, text).unwrap();
        let result = invoke(&input, &output, &[]);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stdout).contains("\"verdict\":\"refused\""));
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
    let input = dir.join("valid.performance");
    std::fs::write(&input, INPUT).unwrap();
    for flag in ["--seconds", "--schedule", "--full-scale-pa"] {
        let output = dir.join(format!("{flag}.wav"));
        assert!(!invoke(&input, &output, &[flag, "1"]).status.success());
        assert!(!output.exists());
    }
    let existing = dir.join("existing.wav");
    std::fs::write(&existing, b"existing audio").unwrap();
    assert!(!invoke(&input, &existing, &[]).status.success());
    assert_eq!(std::fs::read(&existing).unwrap(), b"existing audio");
    let output = dir.join("reserved.wav");
    let sidecar = output.with_extension("provenance.json");
    std::fs::write(&sidecar, b"existing metadata").unwrap();
    assert!(!invoke(&input, &output, &[]).status.success());
    assert!(!output.exists());
    assert_eq!(std::fs::read(&sidecar).unwrap(), b"existing metadata");
}
