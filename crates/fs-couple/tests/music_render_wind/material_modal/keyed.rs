//! The native material key route reaches the actual CLI and sole PCM encoder.
use super::*;

const KEYED: &str = include_str!("../../../examples/material_modal_valve/paired-keyed-sampled.performance");

fn expected(text: &str, key: &NormalizedMaterialCardPack) -> Vec<u8> {
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path|
            Ok(if path == "key.matpack" { key.clone() } else { pack() })).unwrap();
    assert_eq!(loaded.boundary_records().len(), 2);
    assert_eq!(loaded.info().modes, 6);
    let frames = loaded.info().samples as usize;
    let mut source = loaded.into_renderer();
    let mut pressure = vec![0.0; frames];
    for block in pressure.chunks_mut(37) { source.block(block).unwrap(); }
    encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0
}

#[test]
fn native_keys_material_replacement_and_subframe_controls_reach_cli_pcm() {
    let dir = directory();
    let text = KEYED.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
        .replace("boundary_plate blade.matpack", "boundary_plate key.matpack");
    let input = dir.join("keyed.performance");
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), pack().to_bytes()).unwrap();
    let key = pack();
    std::fs::write(dir.join("key.matpack"), key.to_bytes()).unwrap();
    let original = expected(&text, &key);
    assert!(original[44..].iter().any(|byte| *byte != 0));
    for block in [1, 37, 128] {
        let output = dir.join(format!("keyed-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), original);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"retained_modes\":6", "\"full_tensor_viscosity\":true", "\"moving_boundaries\":[",
            "\"side\":\"upstream\",\"source\":{\"node\":3", "\"side\":\"downstream\",\"source\":{\"node\":3",
            "\"boundary_control_clock_hz\":768000",
            "\"mechanical_sample\":193,\"closing_force_n\":2e-3,\"exterior_pressure_pa\":1e-1",
            "\"mechanical_sample\":321,\"closing_force_n\":3e-3,\"exterior_pressure_pa\":3e-1",
            "\"pack_path\":\"key.matpack\"",
            "\"source\":\"local-chimney-gas\"", "\"force_parameter\":\"closing_force\""] {
            assert!(metadata.contains(token), "missing physical source {token}: {metadata}");
        }
        assert!(metadata.contains(&format!("\"pack_hash\":\"{}\"", key.content_hash().to_hex())));
    }
    let heavier = pack_with_density(1800.0);
    std::fs::write(dir.join("key.matpack"), heavier.to_bytes()).unwrap();
    let modified = expected(&text, &heavier);
    let output = dir.join("keyed-new-material.wav");
    success(run(&input, &output, 37, false));
    assert_eq!(std::fs::read(&output).unwrap(), modified);
    assert_ne!(modified, original, "a key-only material substitution must survive quantization");
    let malformed = text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 321 2 0.003 0.3");
    std::fs::write(&input, malformed).unwrap();
    let output = dir.join("invalid-key-address.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
