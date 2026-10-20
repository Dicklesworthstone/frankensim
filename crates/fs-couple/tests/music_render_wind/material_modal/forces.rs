//! Original physical plate loads survive native input, the receiver and PCM.
use super::*;

const FORCED: &str = include_str!("../../../examples/material_modal_valve/paired-forced.performance");

fn expected(text: &str, material: &NormalizedMaterialCardPack) -> Vec<u8> {
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).unwrap();
    assert_eq!(loaded.info().modes, 6);
    assert_eq!(loaded.boundary_records().len(), 2);
    assert_eq!(loaded.renderer().pending_plate_force_controls(0).len(), 4);
    assert_eq!(loaded.renderer().pending_plate_force_controls(1).len(), 2);
    let mut renderer = loaded.into_renderer();
    let mut pressure = vec![0.0; 256];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0
}

#[test]
fn source_plate_force_programs_reach_actual_cli_wav_and_source_records() {
    let dir = directory();
    let text = FORCED.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005");
    let input = dir.join("source-forced.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let reference = expected(&text, &material);
    for block in [1, 37, 128] {
        let output = dir.join(format!("source-forced-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), reference);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"source_force_programs\":[", "\"source_force_clock_hz\":768000",
            "\"node_force_weights\":[", "\"modal_force_projection\":[",
            "\"mechanical_sample\":97,\"force_n\":-1e-4",
            "\"mechanical_sample\":225,\"force_n\":-1.5e-4",
            "\"moving_boundaries\":[", "no tongue or lip contact law"] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
    }
    let heavier = pack_with_density(1150.0);
    std::fs::write(dir.join("blade.matpack"), heavier.to_bytes()).unwrap();
    let output = dir.join("source-forced-new-material.wav");
    success(run(&input, &output, 37, false));
    let changed = expected(&text, &heavier);
    assert_ne!(changed, reference);
    assert_eq!(std::fs::read(output).unwrap(), changed);
    std::fs::write(&input, text.replace("source_force 97 0", "source_force 97 2")).unwrap();
    let output = dir.join("invalid-source-force.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
