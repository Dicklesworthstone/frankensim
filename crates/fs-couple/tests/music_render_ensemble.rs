//! Real file-driven contact and plate physics through the ensemble command.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::render::plate::file::PlatePerformance;
use fs_couple::render::schedule::force::file::ModalPerformance;

const STRIKER: &str = include_str!("../examples/free-striker-192k.performance");
const PLATE: &str = include_str!("../examples/plate-mesh.performance");
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
