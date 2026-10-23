//! Composition and static pressure through the actual CLI and shared PCM owner.
use super::{directory, run, success};
use fs_couple::bernoulli_aperture::performance::{CoupledAperture, file::PlateValvePerformance};
use fs_couple::pcm_wav::{encode_pcm16_wav, decimate::Decimator};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;
use std::process::Command;

const BASE: &str = include_str!("../../examples/plate-valve-mixture.performance");
const CO2: &str = include_str!("../../examples/plate-valve-mixture-co2.performance");
const LOW: &str = include_str!("../../examples/plate-valve-mixture-low-pressure.performance");

fn load(text: &str) -> PlateValvePerformance {
    PlateValvePerformance::from_bytes(text.as_bytes(), 37, &CancelGate::new()).unwrap()
}
fn reference(text: &str) -> Vec<f64> {
    let source = load(text);
    let info = source.info();
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; info.samples as usize];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    let ratio = info.sample_rate_hz as usize / 48000;
    let mut filter = Decimator::new(ratio, 1).unwrap();
    pressure.chunks_exact(ratio).map(|group| {
        let value = filter.preview(group).unwrap()[0]; filter.commit(); value
    }).collect()
}

#[test]
fn mixture_and_pressure_files_reach_exact_cli_pcm_and_resolved_material_metadata() {
    let dir = directory();
    let mut results = Vec::new();
    for (index, text) in [BASE, CO2, LOW].iter().enumerate() {
        // Deliberately isolate atmosphere: no geometry, pressure-phrase or gain edits.
        assert_eq!(BASE.lines().skip(3).collect::<Vec<_>>(), text.lines().skip(3).collect::<Vec<_>>());
        let input = dir.join(format!("mixture-{index}.performance"));
        std::fs::write(&input, text).unwrap();
        let source = load(text);
        let gas = source.ambient_gas();
        let (expected, clips) = encode_pcm16_wav(&reference(text), 48000, source.info().full_scale_pa).unwrap();
        assert!(expected[44..].iter().any(|b| *b != 0));
        for block in [37, 512] {
            let output = dir.join(format!("mixture-{index}-{block}.wav"));
            success(run(&input, &output, block, true));
            assert_eq!(std::fs::read(&output).unwrap(), expected);
            let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
            for token in ["\"ambient_gas\":{", "\"viscothermal\":{", "\"samples\":4801"] {
                assert!(metadata.contains(token), "{token}: {metadata}");
            }
            for (name, value) in [("static_pressure_pa", gas.pressure), ("temperature_k", gas.temperature),
                ("density_kg_m3", gas.density), ("sound_speed_m_s", gas.sound_speed),
                ("specific_heat_cp_j_kg_k", gas.specific_heat_cp), ("thermal_conductivity_w_m_k", gas.thermal_conductivity)] {
                assert!(metadata.contains(&format!("\"{name}\":{value:e}")), "missing {name}: {metadata}");
            }
            assert!(metadata.contains(&format!("\"clipped_samples\":{clips}")));
        }
        results.push(expected);
    }
    assert_ne!(results[0][44..], results[1][44..], "composition must change sound, not only provenance");
    assert_ne!(results[0][44..], results[2][44..], "static pressure must change sound without pitch retuning");
}

#[test]
fn five_species_regional_outlet_retains_its_own_material_in_cli_sound_and_sidecar() {
    let dir = directory();
    let text = BASE.replace("duct_section 1 2 0.125 0.007 0.002 viscothermal",
        "duct_section 1 2 0.125 0.007 0.002 gas mixture 303.15 100000 mole 5 N2 0.70 O2 0.20 Ar 0.05 CO2 0.04 H2O 0.01 viscothermal");
    assert_ne!(text, BASE);
    let input = dir.join("regional-mixture.performance"); std::fs::write(&input, &text).unwrap();
    let output = dir.join("regional-mixture.wav"); success(run(&input, &output, 37, true));
    let (expected, _) = encode_pcm16_wav(&reference(&text), 48000, 0.01).unwrap();
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    assert_ne!(expected, encode_pcm16_wav(&reference(BASE), 48000, 0.01).unwrap().0);
    let p = load(&text);
    let CoupledAperture::Network(n) = p.renderer().system() else { panic!("regional graph required") };
    let local = n.section_gases().unwrap().iter().find(|g| g.water_mole_fraction > 0.0).unwrap();
    assert_ne!(local.density, p.ambient_gas().density);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for (name, value) in [("water_mole_fraction", local.water_mole_fraction),
        ("specific_heat_cp_j_kg_k", local.specific_heat_cp), ("density_kg_m3", local.density)] {
        assert!(metadata.contains(&format!("\"{name}\":{value:e}")), "missing local {name}");
    }
}

#[test]
fn invalid_composition_and_pressure_refuse_before_cli_creates_artifacts() {
    let dir = directory(); let input = dir.join("invalid.performance");
    // He became a supported species in ad216794e; Xe is still unsupported (as in exterior.rs).
    for (i, text) in [BASE.replace("N2 0.7", "Xe 0.7"), BASE.replace("O2 0.3", "O2 0.2"),
        BASE.replace("N2 0.7 O2 0.3", "N2 1 O2 0"), BASE.replace("293.15 100000", "400 100000"),
        BASE.replace("mole 2 N2 0.7 O2 0.3", "mole 2 N2 0.9 H2O 0.1"),
        BASE.replace("0.002 viscothermal", "0.002 gas mixture 300 90000 mole 1 N2 1 viscothermal")]
        .iter().enumerate() {
        std::fs::write(&input, text).unwrap();
        let output = dir.join(format!("invalid-{i}.wav"));
        let result = run(&input, &output, 37, true);
        assert!(!result.status.success(), "accepted {text}");
        assert!(!output.exists()); assert!(!output.with_extension("provenance.json").exists());
        // music_render's refusal is a JSON verdict on stdout (fail()).
        assert!(String::from_utf8_lossy(&result.stdout).contains("line"));
    }
}

#[test]
fn mixture_ensemble_uses_the_same_physical_sources_without_per_source_pcm_gain() {
    let dir = directory(); let base = dir.join("base.performance"); let enriched = dir.join("co2.performance");
    std::fs::write(&base, BASE).unwrap();
    let mut combined = reference(BASE);
    for (out, value) in combined.iter_mut().zip(reference(CO2)) { *out += value; }
    let (expected, _) = encode_pcm16_wav(&combined, 48000, 0.02).unwrap();
    for (index, text) in [CO2.to_string(), CO2.replace("9602 0.01", "9602 100")].iter().enumerate() {
        std::fs::write(&enriched, text).unwrap();
        let output = dir.join(format!("mixture-ensemble-{index}.wav"));
        success(Command::new(env!("CARGO_BIN_EXE_music_render")).arg("ensemble").arg(&output)
            .arg("--valve").arg(&base).arg("--valve").arg(&enriched)
            .args(["--decimate", "--block", "37", "--full-scale-pa", "0.02"]).output().unwrap());
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let meta = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert_eq!(meta.matches("\"ambient_gas\":{").count(), 2);
        assert!(meta.contains("\"source_pcm_scale_applied\":false"));
    }
}

#[path = "gas_mixture/exterior.rs"]
mod exterior;
