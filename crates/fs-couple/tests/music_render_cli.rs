//! E2E for the music-lane render CLI (bead
//! `frankensim-music-t-out-render-ib15w`): drive the REAL binary; assert
//! determinism (same args → bit-identical WAV + provenance), the
//! never-overwrite refusal, the fixture/argument refusals, RIFF shape,
//! provenance completeness, and the never-peak-normalize law (a hotter
//! full-scale yields QUIETER samples of the SAME physics, and clipping is
//! counted, not hidden).

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_music_render")
}

fn scratch(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("music-render-{}-{name}", std::process::id()));
    if dir.exists() {
        std::fs::remove_dir_all(&dir).expect("clean scratch");
    }
    std::fs::create_dir_all(&dir).expect("mkdir scratch");
    dir
}

fn run(args: &[&str]) -> (bool, String) {
    let output = Command::new(bin()).args(args).output().expect("spawn");
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

#[test]
fn renders_deterministically_with_provenance() {
    let dir = scratch("determinism");
    let a = dir.join("a.wav");
    let b = dir.join("b.wav");
    for path in [&a, &b] {
        let (ok, stdout) = run(&["string", path.to_str().expect("utf8"), "--seconds", "0.25"]);
        assert!(ok, "render must succeed:\n{stdout}");
        assert!(stdout.contains("\"verdict\":\"rendered\""), "{stdout}");
        assert!(stdout.contains("\"wav_blake3\":\""), "{stdout}");
    }
    let wav_a = std::fs::read(&a).expect("read a");
    let wav_b = std::fs::read(&b).expect("read b");
    assert_eq!(wav_a, wav_b, "same args must produce bit-identical WAVs");
    assert_eq!(&wav_a[0..4], b"RIFF", "WAV container shape");
    let prov_a = std::fs::read_to_string(a.with_extension("provenance.json")).expect("sidecar");
    let prov_b = std::fs::read_to_string(b.with_extension("provenance.json")).expect("sidecar");
    // Sidecars differ only in nothing — fully deterministic.
    assert_eq!(prov_a, prov_b);
    for field in [
        "\"schema\":\"frankensim-music-render-provenance-v1\"",
        "\"sample_rate_hz\":48000",
        "\"clipped_samples\":",
        "\"wav_blake3\":\"",
        "never peak-normalized",
    ] {
        assert!(
            prov_a.contains(field),
            "sidecar missing {field:?}:\n{prov_a}"
        );
    }
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn refuses_to_overwrite_evidence() {
    let dir = scratch("overwrite");
    let out = dir.join("once.wav");
    let (ok, _) = run(&["string", out.to_str().expect("utf8"), "--seconds", "0.05"]);
    assert!(ok);
    let before = std::fs::read(&out).expect("read");
    let (ok, stdout) = run(&["string", out.to_str().expect("utf8"), "--seconds", "0.05"]);
    assert!(!ok, "re-render onto an existing path must refuse");
    assert!(stdout.contains("refuses to overwrite evidence"), "{stdout}");
    assert_eq!(
        std::fs::read(&out).expect("read"),
        before,
        "the refusal must leave the artifact untouched"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn refusal_arms_are_typed() {
    let dir = scratch("refusals");
    let out = dir.join("x.wav");
    let out_str = out.to_str().expect("utf8");
    for (args, needle) in [
        (vec!["kazoo", out_str], "fixture must be"),
        (vec!["string"], "usage:"),
        (
            vec!["string", out_str, "--seconds", "0"],
            "--seconds must be",
        ),
        (vec!["string", out_str, "--block", "0"], "--block must be"),
    ] {
        let (ok, stdout) = run(&args);
        assert!(!ok, "args {args:?} must refuse");
        assert!(
            stdout.contains(needle),
            "args {args:?}: wrong refusal:\n{stdout}"
        );
    }
    assert!(!out.exists(), "refusals must write nothing");
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn full_scale_is_physics_not_normalization() {
    // The never-peak-normalize law, observable: rendering the same
    // physics at 10x the full-scale (200 kPa vs 20 kPa: the reed peaks ~10.7 kPa, so both scales are clip-free) must produce the SAME peak_pa in
    // provenance (physics unchanged) and quieter PCM (samples scale
    // down); a tiny full-scale must CLIP and say so.
    let dir = scratch("fullscale");
    let quiet = dir.join("quiet.wav");
    let loud = dir.join("loud.wav");
    let clipped = dir.join("clipped.wav");
    let (ok, _) = run(&[
        "reed",
        quiet.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "200000",
    ]);
    assert!(ok);
    let (ok, _) = run(&[
        "reed",
        loud.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "20000",
    ]);
    assert!(ok);
    let (ok, stdout) = run(&[
        "reed",
        clipped.to_str().expect("utf8"),
        "--seconds",
        "0.1",
        "--full-scale-pa",
        "10",
    ]);
    assert!(ok, "clipping renders (counted, not refused):\n{stdout}");

    let read_prov =
        |p: &Path| std::fs::read_to_string(p.with_extension("provenance.json")).expect("sidecar");
    let peak_of = |prov: &str| -> f64 {
        let start = prov.find("\"peak_pa\":").expect("peak field") + "\"peak_pa\":".len();
        let rest = &prov[start..];
        let end = rest.find(',').expect("comma");
        rest[..end].parse().expect("peak parse")
    };
    let clips_of = |prov: &str| -> u64 {
        let start =
            prov.find("\"clipped_samples\":").expect("clip field") + "\"clipped_samples\":".len();
        let rest = &prov[start..];
        let end = rest.find(',').expect("comma");
        rest[..end].parse().expect("clip parse")
    };
    let quiet_prov = read_prov(&quiet);
    let loud_prov = read_prov(&loud);
    let clipped_prov = read_prov(&clipped);
    // Same physics: identical peak pascals across full-scale choices.
    assert!(
        (peak_of(&quiet_prov) - peak_of(&loud_prov)).abs() < 1.0e-12,
        "full-scale must not change the physics"
    );
    assert_eq!(clips_of(&quiet_prov), 0, "200 kPa full-scale must not clip");
    assert_eq!(
        clips_of(&loud_prov),
        0,
        "20 kPa full-scale must not clip either"
    );
    assert!(
        clips_of(&clipped_prov) > 0,
        "10 Pa full-scale must clip a screaming reed and SAY so"
    );
    // Quieter mapping: the 200-kPa WAV's peak sample magnitude is ~10x
    // smaller than the 20-kPa one's (both clip-free, so the ratio is exact).
    let peak_sample = |path: &Path| -> i32 {
        let bytes = std::fs::read(path).expect("wav");
        bytes[44..]
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| i32::from(i16::from_le_bytes(*c)).abs())
            .max()
            .unwrap_or(0)
    };
    let ratio = f64::from(peak_sample(&loud)) / f64::from(peak_sample(&quiet)).max(1.0);
    assert!(
        (8.0..=12.0).contains(&ratio),
        "PCM peaks must scale ~10x with a 10x full-scale change (got {ratio})"
    );
    std::fs::remove_dir_all(&dir).expect("cleanup");
}
