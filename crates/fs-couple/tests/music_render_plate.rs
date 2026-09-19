//! The actual command consumes geometry/materials, renders and hashes a WAV.
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::render::plate::file::PlatePerformance;

const SOURCE: &[u8] = include_bytes!("../examples/plate-mesh.performance");
static NEXT: AtomicU64 = AtomicU64::new(0);
fn scratch() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-music-plate-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path // Retained on failure/success for diagnostics; existing paths never removed.
}
fn run(input: &Path, output: &Path, block: usize) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("plate").arg(input).arg(output).arg("--block").arg(block.to_string())
        .output().unwrap()
}
fn rendered(input: &Path, output: &Path, block: usize) -> Vec<u8> {
    let result = run(input, output, block);
    assert!(result.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
    assert!(String::from_utf8_lossy(&result.stdout).contains("\"verdict\":\"rendered\""));
    std::fs::read(output).unwrap()
}

#[test]
fn command_matches_public_physics_and_pcm_across_partitions_and_relocation() {
    let root = scratch();
    let source = root.join("mesh.performance");
    let moved = root.join("same-mesh-elsewhere.performance");
    std::fs::write(&source, SOURCE).unwrap();
    std::fs::write(&moved, SOURCE).unwrap();
    let mut performance = PlatePerformance::from_bytes(SOURCE, 4801).unwrap().into_renderer();
    let mut pressure = vec![0.0; 4801];
    performance.block(&mut pressure).unwrap();
    let (expected, clips) = encode_pcm16_wav(&pressure, 48_000, 1.0).unwrap();
    assert!(expected[44..].iter().any(|&b| b!=0));
    let hash = fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1", &expected).to_hex();
    let a = root.join("a.wav");
    let b = root.join("b.wav");
    assert_eq!(rendered(&source, &a, 37), expected);
    assert_eq!(rendered(&moved, &b, 571), expected);
    assert_eq!(expected.len(), 44+4801*2);
    for output in [&a,&b] {
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for fragment in ["\"nodes\":30", "\"triangles\":40", "\"sections\":2",
            "\"requested_modes\":3", "\"force_events\":3", "\"samples\":4801"] {
            assert!(sidecar.contains(fragment));
        }
        assert!(sidecar.contains(&format!("\"clipped_samples\":{clips}")));
        assert!(sidecar.contains(&hash));
        assert!(!sidecar.contains(&source.display().to_string()));
    }
    let previous = std::fs::read(&a).unwrap();
    assert!(!run(&source, &a, 64).status.success());
    assert_eq!(std::fs::read(&a).unwrap(), previous);
}

#[test]
fn modifying_a_constitutive_region_changes_the_produced_sound() {
    let root = scratch();
    let a = root.join("first.performance");
    let b = root.join("stiffer.performance");
    std::fs::write(&a, SOURCE).unwrap();
    let changed = std::str::from_utf8(SOURCE).unwrap().replace("11000000000.0", "22000000000.0");
    std::fs::write(&b, changed).unwrap();
    let first = rendered(&a, &root.join("first.wav"), 64);
    let second = rendered(&b, &root.join("stiffer.wav"), 64);
    assert_ne!(&first[44..], &second[44..], "material assignment must reach the audio, not only metadata");
}

#[test]
fn invalid_geometry_and_conflicting_options_leave_no_output_artifact() {
    let root = scratch();
    let input = root.join("invalid.performance");
    let out = root.join("invalid.wav");
    let invalid = std::str::from_utf8(SOURCE).unwrap().replace("triangle 0 1 7 0", "triangle 0 1 7 50");
    std::fs::write(&input, invalid).unwrap();
    assert!(!run(&input, &out, 37).status.success());
    assert!(!out.exists());
    assert!(!out.with_extension("provenance.json").exists());
    std::fs::write(&input, SOURCE).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
        .arg("plate").arg(&input).arg(&out).args(["--seconds", "1"]).output().unwrap();
    assert!(!result.status.success());
    assert!(!out.exists());
}
