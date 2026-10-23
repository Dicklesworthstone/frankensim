//! Named mixtures through the actual TMM reed CLI, not a second audio model.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::render::schedule::reed::ReedPerformance;

const BASE: &str = include_str!("../examples/reed-mixture.performance");
const CO2: &str = include_str!("../examples/reed-mixture-co2.performance");
const LOW: &str = include_str!("../examples/reed-mixture-low-pressure.performance");
const AMBIENT: &str = "ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3";
static NEXT: AtomicU64 = AtomicU64::new(0);
fn directory() -> PathBuf {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-reed-mixture-{}-{stamp}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&path).unwrap();
    path
}
fn success(result: Output) {
    assert!(result.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&result.stdout), String::from_utf8_lossy(&result.stderr));
}
fn wind(input: &Path, output: &Path, block: usize, decimate: bool) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_music_render"));
    command.arg("wind").arg(input).arg(output).args(["--block", &block.to_string()]);
    if decimate { command.arg("--decimate"); }
    command.output().unwrap()
}
fn pressure(text: &str) -> Vec<f64> {
    let mut source = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    let info = source.info();
    let mut raw = vec![0.0; info.samples as usize];
    for block in raw.chunks_mut(37) { source.block(block).unwrap(); }
    assert!(raw.iter().all(|p| p.is_finite()));
    assert!(raw.iter().any(|p| p.abs() > 1e-8));
    let ratio = info.sample_rate_hz as usize / 48000;
    let mut filter = Decimator::new(ratio, 1).unwrap();
    raw.chunks_exact(ratio).map(|block| {
        let value = filter.preview(block).unwrap()[0];
        filter.commit();
        value
    }).collect()
}
fn high_rate(text: &str) -> String {
    // This is an explicit common geometry change for the higher Nyquist band,
    // applied to EVERY composition, never an automatic material-dependent fix.
    text.replace("audio 48000 4801", "audio 96000 9602")
        .replace("cylinder 0.0015", "cylinder 0.00075")
}

#[test]
fn atmosphere_files_change_only_the_medium_and_reach_exact_native_and_decimated_pcm() {
    let dir = directory();
    for high in [false, true] {
        let mut output_bytes = Vec::new();
        for (index, fixture) in [BASE, CO2, LOW].into_iter().enumerate() {
            let lines: Vec<_> = fixture.lines().collect();
            for (line, baseline) in BASE.lines().enumerate() {
                if line != 2 { assert_eq!(lines[line], baseline); }
            }
            let text = if high { high_rate(fixture) } else { fixture.to_string() };
            let input = dir.join(format!("input-{high}-{index}.performance"));
            std::fs::write(&input, &text).unwrap();
            let (expected, _) = encode_pcm16_wav(&pressure(&text), 48000, 20000.0).unwrap();
            assert!(expected[44..].iter().any(|byte| *byte != 0));
            for block in [37, 113] {
                let output = dir.join(format!("render-{high}-{index}-{block}.wav"));
                success(wind(&input, &output, block, high));
                assert_eq!(std::fs::read(&output).unwrap(), expected);
                let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
                assert!(meta.contains("frankensim-reed-performance-v1"));
                assert!(meta.contains("not a calibrated exterior microphone"));
            }
            output_bytes.push(expected);
        }
        assert_ne!(output_bytes[0][44..], output_bytes[1][44..], "composition must change PCM");
        assert_ne!(output_bytes[0][44..], output_bytes[2][44..], "pressure must change PCM");
    }
}

#[test]
fn five_species_mass_basis_and_mixture_ensemble_keep_the_declared_physical_pressure() {
    let dir = directory();
    let all = BASE.replace(AMBIENT,
        "ambient mixture 300 100000 mass 5 nitrogen 0.7 oxygen 0.2 argon 0.05 carbon-dioxide 0.04 water-vapor 0.01");
    let a = dir.join("all-species.performance");
    let b = dir.join("co2.performance");
    std::fs::write(&a, &all).unwrap();
    std::fs::write(&b, CO2).unwrap();
    let mut mixed = pressure(&all);
    for (p, q) in mixed.iter_mut().zip(pressure(CO2)) { *p += q; }
    let (expected, _) = encode_pcm16_wav(&mixed, 48000, 40000.0).unwrap();
    for (index, text) in [all.clone(), all.replace("4801 20000", "4801 1")].iter().enumerate() {
        std::fs::write(&a, text).unwrap();
        let output = dir.join(format!("ensemble-{index}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--reed").arg(&a).arg("--reed").arg(&b)
            .args(["--block", "37", "--full-scale-pa", "40000"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(), expected, "source PCM scales are not source gains");
    }
}

#[test]
fn unsupported_materials_refuse_before_creating_either_cli_output() {
    let dir = directory();
    let input = dir.join("input.performance");
    for (index, gas) in ["mixture 300 100000 mole 1 Xe 1", "mixture 300 100000 mole 2 N2 0.5 nitrogen 0.5",
        "mixture 300 100000 mole 2 N2 70 O2 30", "mixture 300 100000 mole 1 H2O 1",
        "mixture 500 100000 mole 1 N2 1", "mixture 300 0 mole 1 N2 1",
        "mixture 300 100000 mole 1 N2 NaN", "mixture 300 100000 volume 1 N2 1",
        "mixture 300 100000 mole 1 N2 1 extra"].into_iter().enumerate() {
        std::fs::write(&input, BASE.replace(AMBIENT, &format!("ambient {gas}"))).unwrap();
        let output = dir.join(format!("refused-{index}.wav"));
        assert!(!wind(&input, &output, 37, false).status.success(), "accepted {gas}");
        assert!(!output.exists());
        assert!(!output.with_extension("provenance.json").exists());
    }
}

#[path = "reed_gas_mixture/regional.rs"]
mod regional;

#[path = "reed_gas_mixture/outlet.rs"]
mod outlet;
