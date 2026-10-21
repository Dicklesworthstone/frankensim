//! Material-card spectra reach PCM through the same public wind entry point.
use super::*;

fn expected(text: &str, material: &NormalizedMaterialCardPack) -> Vec<u8> {
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(material.clone()),
    )
    .unwrap();
    let mut renderer = loaded.into_renderer();
    let mut pressure = vec![0.0; 256];
    for block in pressure.chunks_mut(37) {
        renderer.block(block).unwrap();
    }
    encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0
}

#[test]
fn original_material_spectrum_changes_actual_cli_pcm_with_full_modal_histories() {
    let dir = directory();
    let text = include_str!("../../../examples/material_modal_valve/paired.performance")
        .replace(
            "audio 768000 48000 2400 0.005",
            "audio 768000 48000 256 0.005",
        )
        .replace(
            "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000",
            "loss proportional-maxwell relaxed 1 maxwell_modulus maxwell_time",
        );
    let input = dir.join("relaxing.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let reference = expected(&text, &material);
    assert!(reference[44..].iter().any(|byte| *byte != 0));
    for block in [1, 37] {
        let output = dir.join(format!("relaxing-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), reference);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in [
            "\"retained_modes\":6",
            "\"full_tensor_viscosity\":false",
            "\"proportional_material_relaxation\":{",
            "\"initial_memory\":\"relaxed\"",
            "\"live_arms\":6",
            "\"regional_spectra\":[",
            "\"modulus_pa\":1e8",
            "\"relaxation_time_s\":1e-3",
            "proportional Maxwell internal strains",
        ] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
        assert!(metadata.contains(&material.content_hash().to_hex().to_string()));
    }
    let changed = pack_with_relaxation(900.0, 1e9, 0.0002);
    std::fs::write(dir.join("blade.matpack"), changed.to_bytes()).unwrap();
    let output = dir.join("changed-material-spectrum.wav");
    success(run(&input, &output, 37, false));
    let changed_expected = expected(&text, &changed);
    assert_ne!(changed_expected[44..], reference[44..]);
    assert_eq!(std::fs::read(output).unwrap(), changed_expected);
}
