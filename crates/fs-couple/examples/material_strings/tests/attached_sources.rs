//! The installed performer and web source packet use the same decoded physics.
use super::*;
use std::{ffi::OsString, path::{Path, PathBuf}, sync::atomic::{AtomicU64, Ordering}};
#[path = "../../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

static NEXT: AtomicU64 = AtomicU64::new(0);

fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-string-sources-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path
}

fn fixture(dir: &Path) {
    let source = include_str!("../pluck.material-strings")
        .replace("48000 4096 0.01 128", "48000 256 0.001 37");
    std::fs::write(dir.join("input.performance"), source).unwrap();
    for name in ["wire.fsmatpk", "body.fsmatpk"] {
        std::fs::write(dir.join(name), pack(name, 7800.0, 200_000.0, None).to_bytes()).unwrap();
    }
}

fn args(dir: &Path, output: &str, attached: Option<&[&str]>) -> Vec<OsString> {
    let mut args = vec!["material_strings".into(), dir.join("input.performance").into_os_string(),
        dir.join(output).into_os_string()];
    if let Some(names) = attached {
        args.push("--self-contained".into());
        for name in names { args.extend(["--attached-source".into(), (*name).into()]); }
    }
    args
}

#[test]
fn g3_attached_material_sources_match_ordinary_pcm_and_preserve_data_only_swaps() {
    let dir = directory(); fixture(&dir);
    run(&args(&dir, "ordinary.wav", None)).unwrap();
    let names = ["wire.fsmatpk", "body.fsmatpk"];
    run(&args(&dir, "attached.wav", Some(&names))).unwrap();
    let original = std::fs::read(dir.join("attached.wav")).unwrap();
    assert_eq!(original, std::fs::read(dir.join("ordinary.wav")).unwrap());
    assert!(original[44..].iter().any(|byte| *byte != 0));
    let source = std::fs::read(dir.join("input.performance")).unwrap();
    std::fs::write(dir.join("wire.fsmatpk"), pack("wire.fsmatpk", 15600.0, 200_000.0, None).to_bytes()).unwrap();
    run(&args(&dir, "replacement.wav", Some(&names))).unwrap();
    assert_ne!(original[44..], std::fs::read(dir.join("replacement.wav")).unwrap()[44..]);
    assert_eq!(source, std::fs::read(dir.join("input.performance")).unwrap());
}

#[test]
fn g3_attached_gas_sources_share_the_same_string_and_cavity_state() {
    let dir = directory(); fixture(&dir);
    let source = std::fs::read_to_string(dir.join("input.performance")).unwrap().replace(
        "gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
        "gas thermochemical 293.15 101325 mole 1 bow-atmosphere 1");
    std::fs::write(dir.join("input.performance"), source).unwrap();
    let packs = gas_support::packs("bow-atmosphere", 0.028, 3.5, 200.0);
    for (name, bytes) in ["species.pack", "nasa9.pack", "transport.pack"].into_iter().zip([
        packs.species.to_bytes(), packs.nasa9.to_bytes(), packs.transport.to_bytes(),
    ]) { std::fs::write(dir.join(name), bytes).unwrap(); }
    for restricted in [false, true] {
        let names = ["wire.fsmatpk", "body.fsmatpk", "species.pack", "nasa9.pack", "transport.pack"];
        let mut command = args(&dir, if restricted { "attached.wav" } else { "ordinary.wav" },
            restricted.then_some(names.as_slice()));
        command.extend(["--gas-component".into(),
            if restricted { "species.pack".into() } else { dir.join("species.pack").into_os_string() },
            if restricted { "nasa9.pack".into() } else { dir.join("nasa9.pack").into_os_string() },
            if restricted { "transport.pack".into() } else { dir.join("transport.pack").into_os_string() },
            "--gas-water".into(), "none".into(), "--gas-phase".into(), "frozen".into()]);
        run(&command).unwrap();
    }
    assert_eq!(std::fs::read(dir.join("attached.wav")).unwrap(),
        std::fs::read(dir.join("ordinary.wav")).unwrap());
}

#[test]
fn g3_decimated_attached_strings_preserve_output_clock_and_callback_identity() {
    let dir = directory(); fixture(&dir);
    let original = std::fs::read_to_string(dir.join("input.performance")).unwrap();
    let names = ["wire.fsmatpk", "body.fsmatpk"];
    for (block, output) in [(37, "prime.wav"), (128, "wide.wav")] {
        let source = original.replace("48000 256 0.001 37", &format!("96000 512 0.001 {block}"));
        std::fs::write(dir.join("input.performance"), source).unwrap();
        let mut command = args(&dir, output, Some(&names));
        command.push("--decimate".into());
        run(&command).unwrap();
    }
    let wav = std::fs::read(dir.join("prime.wav")).unwrap();
    assert_eq!(wav, std::fs::read(dir.join("wide.wav")).unwrap());
    assert_eq!(u32::from_le_bytes(wav[24..28].try_into().unwrap()), 48_000);
    assert_eq!(u32::from_le_bytes(wav[40..44].try_into().unwrap()), 256 * 2);
    assert_eq!(wav.len(), 44 + 256 * 2);
    assert!(wav[44..].iter().any(|byte| *byte != 0));

    // An incomplete observation interval must not silently shorten the gesture
    // horizon or publish an output with an unrendered mechanical tail.
    std::fs::write(dir.join("input.performance"),
        original.replace("48000 256 0.001 37", "96000 513 0.001 37")).unwrap();
    let mut command = args(&dir, "partial.wav", Some(&names));
    command.push("--decimate".into());
    assert!(run(&command).unwrap_err().to_string().contains("complete output intervals"));
    assert!(!dir.join("partial.wav").exists());
}

#[test]
fn g0_unlisted_material_profiles_and_cache_sources_refuse_before_output() {
    let dir = directory(); fixture(&dir);
    let mut missing = args(&dir, "missing.wav", Some(&["body.fsmatpk"]));
    assert!(run(&missing).unwrap_err().to_string().contains("unlisted source"));
    assert!(!dir.join("missing.wav").exists());
    missing = args(&dir, "profile.wav", Some(&["body.fsmatpk", "wire.fsmatpk"]));
    std::fs::write(dir.join("body.profile"), b"not admitted").unwrap();
    missing.extend(["--body-profile".into(), "body.profile".into()]);
    assert!(run(&missing).unwrap_err().to_string().contains("unlisted source"));
    assert!(!dir.join("profile.wav").exists());
    let mut cached = args(&dir, "cache.wav", Some(&[]));
    cached.extend(["--observer-store".into(), dir.clone().into_os_string()]);
    assert!(run(&cached).unwrap_err().to_string().contains("refuses --observer-store"));
    assert!(!dir.join("cache.wav").exists());
}
