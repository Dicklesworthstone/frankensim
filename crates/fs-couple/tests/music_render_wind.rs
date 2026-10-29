//! Supplied wind performances through the actual CLI, observer and sole PCM owner.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::{encode_pcm16_wav, decimate::Decimator};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::plate::file::PlatePerformance;
use fs_couple::render::schedule::reed::ReedPerformance;

#[path = "music_render_wind/resolved.rs"]
mod resolved;

const REED: &str = include_str!("../examples/reed-duct.performance");
const PLATE: &str = include_str!("../examples/plate-mesh.performance");
static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-wind-file-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&dir).unwrap(); dir
}
fn high_rate() -> String {
    // Explicit geometry change keeps the existing low-ka load inside its domain
    // all the way to the higher clock's Nyquist frequency. No hidden retuning.
    REED.replace("audio 48000 4801", "audio 96000 9602")
        .replace("cylinder 0.0022", "cylinder 0.0011")
}
fn reference(text: &str) -> (Vec<f64>, usize) {
    let mut source = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    let info = source.info();
    let mut physical = vec![0.0; info.samples as usize];
    for block in physical.chunks_mut(37) { source.block(block).unwrap(); }
    let ratio = info.sample_rate_hz as usize / 48_000;
    let mut filter = Decimator::new(ratio, 1).unwrap();
    let delay = filter.delay_output_frames() as usize;
    let observed = physical.chunks_exact(ratio).map(|group| {
        let value = filter.preview(group).unwrap()[0]; filter.commit(); value
    }).collect();
    (observed, delay)
}
fn run(input: &Path, output: &Path, block: usize, decimate: bool) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
    command.arg("wind").arg(input).arg(output).args(["--block", &block.to_string()]);
    if decimate { command.arg("--decimate"); }
    command.output().unwrap()
}
fn success(output: Output) {
    assert!(output.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
}
// music_render reports a refusal as one JSON line on stdout; this is its decoded `what`.
fn refusal_what(output: &Output) -> String {
    let line = String::from_utf8_lossy(&output.stdout);
    let what = line.split_once("\"verdict\":\"refused\",\"what\":\"").map_or("", |(_, w)| w);
    let (mut text, mut chars) = (String::new(), what.chars());
    while let Some(c) = chars.next() {
        match c {
            '"' => break,
            '\\' => match chars.next() { Some('n') => text.push('\n'), Some(e) => text.push(e), None => break },
            _ => text.push(c),
        }
    }
    text
}

#[test]
fn wind_command_writes_exact_physical_pcm_on_native_and_explicitly_decimated_clocks() {
    let dir = directory();
    for (i, text) in [REED.to_string(), high_rate()].iter().enumerate() {
        let input = dir.join(format!("input-{i}.performance")); std::fs::write(&input, text).unwrap();
        let (pressure, delay) = reference(text);
        let (expected, clips) = encode_pcm16_wav(&pressure, 48_000, 20000.0).unwrap();
        let hash = fs_blake3::hash_domain("org.frankensim.fs-couple.music-render-wav.v1", &expected).to_hex();
        for block in [1, 37, 512] {
            let output = dir.join(format!("wind-{i}-{block}.wav"));
            success(run(&input, &output, block, i == 1));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["frankensim-reed-performance-v1", "\"gesture_events\":5",
                "\"samples\":4801", "not a calibrated exterior microphone", "\"tail\":\"no-flush-declared-window\""] {
                assert!(sidecar.contains(token), "{token}: {sidecar}");
            }
            assert!(sidecar.contains(&hash));
            assert!(sidecar.contains(&format!("\"delay_output_samples\":{delay}")));
            assert!(sidecar.contains(&format!("\"clipped_samples\":{clips}")));
        }
        let relocated = dir.join(format!("relocated-{i}.performance")); std::fs::write(&relocated, text).unwrap();
        let replay = dir.join(format!("replay-{i}.wav")); success(run(&relocated, &replay, 37, i == 1));
        assert_eq!(std::fs::read(replay.with_extension("provenance.json")).unwrap(),
            std::fs::read(dir.join(format!("wind-{i}-37.provenance.json"))).unwrap());
        assert!(!run(&input, &replay, 37, i == 1).status.success());
        assert_eq!(std::fs::read(&replay).unwrap(), expected);
    }
}

#[test]
fn reed_and_mesh_mix_with_independent_clocks_without_hidden_source_pcm_gains() {
    let dir = directory();
    let input = dir.join("reed.performance"); let plate_path = dir.join("plate.performance");
    std::fs::write(&plate_path, PLATE).unwrap();
    let mut plate = PlatePerformance::from_bytes(PLATE.as_bytes(), 37).unwrap().into_renderer();
    let mut plate_pressure = vec![0.0; 4801];
    for block in plate_pressure.chunks_mut(37) { plate.block(block).unwrap(); }
    let mut outputs = Vec::new();
    for (i, text) in [high_rate(), high_rate().replace("9602 20000", "9602 1"),
        high_rate().replace("ambient 293.15", "ambient 303.15")].iter().enumerate() {
        std::fs::write(&input, text).unwrap();
        let (mut expected, delay) = reference(text);
        for (p, q) in expected[delay..].iter_mut().zip(&plate_pressure) { *p += q; }
        let (wav, _) = encode_pcm16_wav(&expected, 48_000, 20000.0).unwrap();
        let output = dir.join(format!("mixed-{i}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--reed").arg(&input).arg("--plate").arg(&plate_path)
            .args(["--block", "37", "--decimate", "--full-scale-pa", "20000"]).output().unwrap());
        let actual = std::fs::read(&output).unwrap(); assert_eq!(actual, wav); outputs.push(actual);
        let sidecar = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"kind\":\"reed\"", "\"kind\":\"plate\"", "\"mechanics_samples\":9602",
            "\"common_delay_output_samples\":40", "\"alignment_delay_samples\":40",
            "\"source_pcm_scale_applied\":false", "bore-pressure-plus-compact-jet-proxy"] {
            assert!(sidecar.contains(token), "{token}: {sidecar}");
        }
    }
    assert_eq!(outputs[0], outputs[1], "source PCM scale is not a gain");
    assert_ne!(outputs[0][44..], outputs[2][44..], "gas changes must reach physical PCM, not only metadata");
}

#[test]
fn invalid_wind_sources_clocks_and_overrides_refuse_before_creating_output_files() {
    let dir = directory(); let input = dir.join("input.performance");
    for (i, (text, decimate)) in [
        (high_rate(), false),
        (REED.to_string(), true),
        (high_rate().replace("9602 20000", "9603 20000"), true),
        (high_rate().replace("audio 96000", "audio 100000"), true),
        (REED.replace("ambient 293.15", "ambient NaN"), false),
        (format!("{REED}ignored\n"), false),
    ].iter().enumerate() {
        std::fs::write(&input, text).unwrap(); let output = dir.join(format!("refused-{i}.wav"));
        assert!(!run(&input, &output, 37, *decimate).status.success());
        assert!(!output.exists()); assert!(!output.with_extension("provenance.json").exists());
    }
    std::fs::write(&input, REED).unwrap();
    for (i, flag) in ["--seconds", "--schedule", "--full-scale-pa", "--gain", "--frequency"].iter().enumerate() {
        let output = dir.join(format!("override-{i}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render")).arg("wind").arg(&input).arg(&output)
            .args([*flag, "1"]).output().unwrap();
        assert!(!result.status.success()); assert!(!output.exists());
    }
}

#[test]
fn self_contained_reed_and_flues_keep_the_same_physical_pcm_and_provenance() {
    let dir = directory();
    for (index, (source, decimate)) in [
        (REED, false),
        (include_str!("../examples/flue-dg.performance"), true),
        (include_str!("../examples/flue-nolle.performance"), false),
    ].iter().enumerate() {
        let input = dir.join(format!("inline-{index}.performance"));
        std::fs::write(&input, source).unwrap();
        let ordinary = dir.join(format!("ordinary-{index}.wav"));
        success(run(&input, &ordinary, 37, *decimate));
        let restricted = dir.join(format!("self-contained-{index}.wav"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
        command.arg("wind").arg(&input).arg(&restricted)
            .args(["--block", "37", "--self-contained"]);
        if *decimate { command.arg("--decimate"); }
        success(command.output().unwrap());
        assert_eq!(std::fs::read(&restricted).unwrap(), std::fs::read(&ordinary).unwrap());
        assert_eq!(std::fs::read(restricted.with_extension("provenance.json")).unwrap(),
            std::fs::read(ordinary.with_extension("provenance.json")).unwrap());
    }
}

#[test]
fn self_contained_dg_wall_refuses_before_opening_missing_or_existing_source_packs() {
    let dir = directory();
    let input = dir.join("external-wall.performance");
    let source = include_str!("../examples/flue-dg.performance").replace("\nflue ",
        "\ndg_wall isotropic wall.pack 293.15 1000 0.001 stress-free 50 2000 elastic\nflue ");
    std::fs::write(&input, source).unwrap();
    for present in [false, true] {
        if present {
            std::fs::write(dir.join("wall.pack"), b"not a material-card pack").unwrap();
        }
        let output = dir.join(format!("restricted-{present}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("wind").arg(&input).arg(&output)
            .args(["--decimate", "--self-contained"]).output().unwrap();
        assert!(!result.status.success());
        let refusal = refusal_what(&result);
        assert!(refusal.contains("--self-contained refuses external material-card pack \"wall.pack\""),
            "{refusal}");
        assert!(!refusal.contains("open failed"), "{refusal}");
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());

        // The ordinary CLI still resolves the same declared artifact. A missing
        // file reports open failure; an existing malformed file reaches decoding.
        let ordinary = dir.join(format!("ordinary-wall-{present}.wav"));
        let result = run(&input, &ordinary, 37, true);
        assert!(!result.status.success());
        let refusal = refusal_what(&result);
        assert!(refusal.contains("material-card pack \"wall.pack\""), "{refusal}");
        assert!(!refusal.contains("--self-contained"), "{refusal}");
        assert_eq!(refusal.contains("open failed"), !present, "{refusal}");
        assert!(!ordinary.exists());
        assert!(!ordinary.with_extension("provenance.json").exists());
    }
}

#[test]
fn self_contained_duplicate_and_external_gas_options_refuse_before_source_file_access() {
    let dir = directory();
    let input = dir.join("inline.performance");
    std::fs::write(&input, REED).unwrap();
    for (index, (flags, expected)) in [
        (vec!["--self-contained", "--self-contained"],
            "--self-contained may only be supplied once"),
        (vec!["--self-contained", "--gas-record", "missing.gas"],
            "--self-contained refuses external gas source files"),
        (vec!["--self-contained", "--gas-component", "missing-species.pack",
            "missing-nasa.pack", "missing-transport.pack", "--gas-water", "none",
            "--gas-phase", "frozen"],
            "--self-contained refuses external gas source files"),
    ].iter().enumerate() {
        let output = dir.join(format!("refused-self-contained-{index}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("wind").arg(&input).arg(&output).args(flags).output().unwrap();
        assert!(!result.status.success());
        let refusal = refusal_what(&result);
        assert!(refusal.contains(*expected), "{refusal}");
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[path = "music_render_wind/valve.rs"]
mod valve;

#[path = "music_render_wind/gas_mixture.rs"]
mod gas_mixture;

#[path = "music_render_wind/flue.rs"]
mod flue;

#[path = "music_render_wind/material_modal.rs"]
mod material_modal;

#[path = "music_render_wind/dg.rs"]
mod dg;

#[path = "music_render_wind/flue_dg.rs"]
mod flue_dg;


#[test]
fn attached_materials_are_explicit_and_cannot_enable_ambient_source_callbacks() {
    let dir = directory();
    let input = dir.join("attached-wall.performance");
    let source = include_str!("../examples/flue-dg.performance").replace("\nflue ",
        "\ndg_wall isotropic wall.pack 293.15 1000 0.001 stress-free 50 2000 elastic\nflue ");
    std::fs::write(&input, source).unwrap();
    std::fs::write(dir.join("wall.pack"), b"original malformed pack, for loader admission only").unwrap();
    std::fs::write(dir.join("other.pack"), b"unrelated source").unwrap();
    for (index, flags, expected) in [
        (0, vec!["--attached-source", "wall.pack"], "--attached-source requires --self-contained"),
        (1, vec!["--self-contained", "--attached-source", "other.pack"],
            "--self-contained refuses external material-card pack"),
        (2, vec!["--self-contained", "--attached-source", "wall.pack"], "material-card pack"),
        (3, vec!["--self-contained", "--attached-source", "../wall.pack"], "relative path"),
        (4, vec!["--self-contained", "--attached-source", "wall.pack",
            "--attached-source", "wall.pack"], "duplicate --attached-source"),
    ] {
        let output = dir.join(format!("attached-refusal-{index}.wav"));
        let result = Command::new(env!("CARGO_BIN_EXE_music_render"))
            .arg("wind").arg(&input).arg(&output).arg("--decimate").args(flags).output().unwrap();
        assert!(!result.status.success());
        let refusal = refusal_what(&result);
        assert!(refusal.contains(expected), "{refusal}");
        if index == 2 {
            assert!(!refusal.contains("--self-contained refuses"), "{refusal}");
            assert!(!refusal.contains("open failed"), "{refusal}");
        }
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}
