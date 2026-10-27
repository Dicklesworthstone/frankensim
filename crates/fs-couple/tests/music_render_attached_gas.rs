//! Attached atmosphere sources through the actual native CLI and WAV owner.
//! Source fixtures are synthetic; these are integration tests, not realism grades.
#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod support;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);
const PATHS: [&str; 3] = ["gas/species.pack", "gas/nasa9.pack", "gas/transport.pack"];
const RECORD: &str = "thermochemical 300 101325 mole 1 custom-alpha 1";

fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
        .unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-attached-gas-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    std::fs::create_dir(path.join("gas")).unwrap();
    path
}
fn write_sources(dir: &Path, mass: f64, legacy: bool) -> String {
    let packs = support::packs("custom-alpha", mass, 3.5, 250.0);
    for (name, bytes) in PATHS.into_iter().zip([
        packs.species.to_bytes(), packs.nasa9.to_bytes(),
        if legacy { packs.transport.claims_pack().to_bytes() } else { packs.transport.to_bytes() },
    ]) { std::fs::write(dir.join(name), bytes).unwrap(); }
    packs.transport.content_hash().to_hex()
}
fn command(dir: &Path, output: &str, decimate: bool, attached: Option<&[&str]>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_music_render"));
    cmd.current_dir(dir).args(["wind", "input.performance", output, "--block", "37",
        "--gas-component", PATHS[0], PATHS[1], PATHS[2],
        "--gas-water", "none", "--gas-phase", "frozen"]);
    if decimate { cmd.arg("--decimate"); }
    if let Some(names) = attached {
        cmd.arg("--self-contained");
        for name in names { cmd.args(["--attached-source", name]); }
    }
    cmd.output().unwrap()
}
fn diagnostic(out: &Output) -> String {
    format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr))
}
fn success(out: Output) { assert!(out.status.success(), "{}", diagnostic(&out)); }
fn reed() -> String {
    include_str!("../examples/reed-duct.performance")
        .replace("ambient 293.15 101325 0", &format!("ambient {RECORD}"))
}

#[test]
fn attached_and_ordinary_gas_packs_produce_identical_native_reed_and_flue_audio() {
    let flue = include_str!("../examples/flue-dg.performance")
        .replace("gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2", &format!("gas {RECORD}"));
    for (source, decimate) in [(reed(), false), (flue, true)] {
        for legacy in [false, true] {
            let dir = directory();
            let identity = write_sources(&dir, 0.028, legacy);
            std::fs::write(dir.join("input.performance"), &source).unwrap();
            success(command(&dir, "ordinary.wav", decimate, None));
            success(command(&dir, "attached.wav", decimate, Some(&PATHS)));
            let audio = std::fs::read(dir.join("attached.wav")).unwrap();
            assert_eq!(audio, std::fs::read(dir.join("ordinary.wav")).unwrap());
            assert!(audio[44..].iter().any(|byte| *byte != 0));
            let metadata = std::fs::read(dir.join("attached.provenance.json")).unwrap();
            assert_eq!(metadata, std::fs::read(dir.join("ordinary.provenance.json")).unwrap());
            if !legacy { assert!(String::from_utf8(metadata).unwrap().contains(&identity)); }
        }
    }
}

#[test]
fn a_source_only_atmosphere_swap_changes_actual_pcm_without_rewriting_the_performance() {
    let dir = directory();
    let source = reed();
    std::fs::write(dir.join("input.performance"), &source).unwrap();
    write_sources(&dir, 0.028, false);
    success(command(&dir, "first.wav", false, Some(&PATHS)));
    write_sources(&dir, 0.032, false);
    success(command(&dir, "second.wav", false, Some(&PATHS)));
    assert_ne!(std::fs::read(dir.join("first.wav")).unwrap()[44..],
        std::fs::read(dir.join("second.wav")).unwrap()[44..]);
    assert_eq!(std::fs::read_to_string(dir.join("input.performance")).unwrap(), source);
}

#[test]
fn existing_but_unlisted_gas_sources_refuse_without_creating_audio_or_sidecars() {
    let dir = directory();
    write_sources(&dir, 0.028, false);
    std::fs::write(dir.join("input.performance"), reed()).unwrap();
    for (index, names) in [PATHS[..0].to_vec(), PATHS[..2].to_vec()].iter().enumerate() {
        let output = format!("missing-{index}.wav");
        let out = command(&dir, &output, false, Some(names));
        assert!(!out.status.success());
        assert!(diagnostic(&out).contains("--self-contained refuses external gas source files"));
        assert!(!dir.join(&output).exists());
        assert!(!dir.join(&output).with_extension("provenance.json").exists());
    }
}
