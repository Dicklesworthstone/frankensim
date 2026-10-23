//! Original pressure-face geometry through the actual CLI and PCM owner.
use super::*;

const WETTED: &str = include_str!("../../../examples/material_modal_valve/paired-wetted.performance");

fn expected(text: &str, material: &NormalizedMaterialCardPack) -> (Vec<u8>, String) {
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).unwrap();
    let fs_couple::bernoulli_aperture::two_sided::LoadedAperture::Modal(valve) =
        loaded.renderer().source().system().valve() else { panic!("material modal owner"); };
    let port = valve.material_plates()[0].wetted_surface().unwrap();
    assert_eq!(port.surface().upstream_facets().len(), 12);
    let specimen = port.surface().specimen_identity().to_hex();
    let mut renderer = loaded.into_renderer();
    let mut pressure = vec![0.0; 256];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    (encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0, specimen)
}

#[test]
fn wetted_geometry_and_material_replacement_reach_actual_cli_sound_and_source_records() {
    let dir = directory();
    let text = WETTED.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005");
    let input = dir.join("wetted.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let (reference, specimen) = expected(&text, &material);
    for block in [1, 37] {
        let output = dir.join(format!("wetted-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), reference);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"wetted_surfaces\":[", "manufactured-exposed-tip", "\"triangle\":6",
            "\"upstream_face\":\"positive\"", "\"modal_pressure_areas_m2\":[", specimen.as_str()] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
    }
    let flipped = text.replace(" positive\n", " negative\n");
    std::fs::write(&input, &flipped).unwrap();
    let output = dir.join("wetted-opposite-face.wav");
    success(run(&input, &output, 37, false));
    let (opposite, _) = expected(&flipped, &material);
    assert_eq!(std::fs::read(&output).unwrap(), opposite);
    assert_ne!(opposite[44..], reference[44..]);
    std::fs::write(&input, &text).unwrap();
    let heavier = pack_with_density(1150.0);
    std::fs::write(dir.join("blade.matpack"), heavier.to_bytes()).unwrap();
    let output = dir.join("wetted-new-material.wav");
    success(run(&input, &output, 37, false));
    let (changed, _) = expected(&text, &heavier);
    assert_eq!(std::fs::read(&output).unwrap(), changed);
    assert_ne!(changed[44..], reference[44..]);
    std::fs::write(&input, text.replace("wetted_facet 7 positive", "wetted_facet 6 negative")).unwrap();
    let output = dir.join("duplicate-wetted-face.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}
