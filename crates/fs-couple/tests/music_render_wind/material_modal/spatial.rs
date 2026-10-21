//! Source-mesh closure must run through the ordinary wind command and PCM owner.
use super::*;

fn expected(text: &str, material: &NormalizedMaterialCardPack, block: usize) -> Vec<u8> {
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(material.clone()),
    )
    .unwrap();
    let mut source = loaded.into_renderer();
    let mut pressure = vec![0.0; 256];
    for values in pressure.chunks_mut(block) {
        source.block(values).unwrap();
    }
    encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0
}

#[test]
fn spatial_material_reeds_reach_actual_cli_pcm_without_a_second_physics_path() {
    let dir = directory();
    let text = include_str!("../../../examples/material_modal_valve/paired-spatial.performance")
        .replace(
            "audio 768000 48000 2400 0.005",
            "audio 768000 48000 256 0.005",
        );
    let input = dir.join("spatial.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let reference = expected(&text, &material, 37);
    assert!(reference[44..].iter().any(|byte| *byte != 0));
    for block in [1, 37] {
        let output = dir.join(format!("spatial-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), reference);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        assert!(metadata.contains(
            "source-mesh distributed contact and positive-P1 slit flow on every retained mode"
        ));
        assert!(metadata.contains("\"retained_modes\":6"));
        assert!(metadata.contains(&material.content_hash().to_hex().to_string()));
    }
    // The same geometry, source program and receiver consume a data-only swap.
    let altered = pack_with_density(1150.0);
    std::fs::write(dir.join("blade.matpack"), altered.to_bytes()).unwrap();
    let changed = expected(&text, &altered, 37);
    assert_ne!(reference[44..], changed[44..]);
    let output = dir.join("changed-spatial-material.wav");
    success(run(&input, &output, 37, false));
    assert_eq!(std::fs::read(output).unwrap(), changed);
    let broken = text.replace("lay_triangle 23", "lay_triangle 22");
    std::fs::write(&input, broken).unwrap();
    let output = dir.join("invalid-spatial.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
