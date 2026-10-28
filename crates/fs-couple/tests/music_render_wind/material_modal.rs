//! Original material modes and imported exterior through the actual CLI.
use super::{directory, run, success};
use fs_couple::bernoulli_aperture::performance::file::material_modal::MaterialModalValvePerformance;
use fs_couple::pcm_wav::{encode_pcm16_wav, observation::PressureRenderer};
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const INPUT: &str = include_str!("../../examples/material_modal_valve/paired-mutual.performance");
const ARTIFACT: &str = include_str!("../../examples/material_modal_valve/mutual.radiation");

#[test]
fn circular_supply_reaches_cli_pcm_and_geometry_metadata_without_an_instrument_solver() {
    let dir = directory();
    let text = include_str!("../../examples/material_modal_valve/paired.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
        .replace("supply none", "supply circular 0.02 0.0005 32");
    let input = dir.join("circular-supply.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let source = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone())).unwrap();
    let mut source = source.into_renderer();
    let mut pressure = vec![0.0; 512];
    for block in pressure.chunks_mut(37) { source.block(block).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, 0.005).unwrap().0;
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    for block in [1, 37] {
        let output = dir.join(format!("circular-supply-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"kind\":\"circular\"", "\"radius_m\":5e-4", "\"radial_cells\":32"] {
            assert!(metadata.contains(token), "missing {token}: {metadata}");
        }
    }
}

#[test]
fn material_key_file_reaches_actual_wav_through_the_existing_reed_renderer() {
    let dir = directory();
    let text = include_str!("../../examples/material_modal_valve/paired-keyed.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005");
    let input = dir.join("keyed.performance");
    let blade = pack();
    let key = pack_with_density(1800.0);
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), blade.to_bytes()).unwrap();
    std::fs::write(dir.join("key.matpack"), key.to_bytes()).unwrap();
    let loaded = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None,
        |path| Ok(if path == "key.matpack" { key.clone() } else { blade.clone() })).unwrap();
    assert_eq!(loaded.boundary_records().len(), 1);
    let mut renderer = loaded.into_renderer();
    let mut pressure = vec![0.0; 512];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, 0.005).unwrap().0;
    for block in [1, 37] {
        let output = dir.join(format!("keys-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"moving_boundaries\":[", "\"side\":\"downstream\"", "\"node\":3",
            "\"pack_path\":\"key.matpack\"", "\"force_parameter\":\"closing_force\"", "local-chimney-gas"] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
        assert!(metadata.contains(&format!("\"pack_hash\":\"{}\"", key.content_hash().to_hex())));
    }
    // Only the key material artifact changes; the source file, reed cards,
    // both ducts, gesture program and microphone remain the original values.
    std::fs::write(dir.join("key.matpack"), blade.to_bytes()).unwrap();
    let changed = dir.join("changed-key.wav");
    success(run(&input, &changed, 37, false));
    assert_ne!(std::fs::read(&changed).unwrap()[44..], expected[44..]);
}

fn pack() -> NormalizedMaterialCardPack {
    pack_with_density(900.0)
}

#[test]
fn native_material_dg_observer_gas_preserves_source_and_pcm() {
    use fs_couple::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance;
    use fs_couple::pcm_wav::observation::DecimatedRenderer;
    let text = include_str!("../../examples/material_modal_valve/paired-dg-chimney.performance")
        .replace("audio 96000 96000 4800 20", "audio 96000 96000 512 0.02");
    let material = pack();
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(), 37, &CancelGate::new_clock_free(), None, |_| Ok(material.clone()),
    ).unwrap();
    let report = source.renderer().observer_gas().unwrap();
    assert_eq!(report.gas, *source.renderer().system().load().radiation_record()
        .unwrap().gas_record().state());
    assert_eq!(report.surface_filter_indices.len(), 2*8*32);
    let mut renderer = DecimatedRenderer::new(source.into_renderer(), 96000, 48000, 37).unwrap();
    let mut pressure = vec![0.0;256];
    for block in pressure.chunks_mut(37) { renderer.block(block).unwrap(); }
    let expected = encode_pcm16_wav(&pressure, 48000, 0.02).unwrap().0;
    let dir = directory();
    let input = dir.join("reed-observer-gas.performance");
    std::fs::write(&input, text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    let output = dir.join("reed-observer-gas.wav");
    success(run(&input, &output, 37, true));
    assert_eq!(std::fs::read(&output).unwrap(), expected);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    for token in ["\"observer_gas\":{\"law\":\"stokes-kirchhoff\"",
        "\"retained_modes\":6", "\"molecular_relaxation_omitted\":true"] {
        assert!(metadata.contains(token), "{token}");
    }
}
fn pack_with_density(density: f64) -> NormalizedMaterialCardPack {
    pack_with_relaxation(density, 1e8, 0.001)
}
fn pack_with_relaxation(density: f64, modulus: f64, time: f64) -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let artifact = fs_blake3::hash_domain("org.frankensim.test.material-modal-cli", b"manufactured-data");
    let provenance = Provenance { source: "manufactured material-modal CLI regression".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured material-modal valve".into(), method: "synthetic integration fixture".into(),
        artifact, caveats: "not a measured material or instrument".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [
        ("density", density, Density::DIMS), ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE), ("eta11", 22000.0, DynViscosity::DIMS),
        ("eta22", 35000.0, DynViscosity::DIMS), ("eta66", 11000.0, DynViscosity::DIMS),
        ("eta12", 4500.0, DynViscosity::DIMS), ("eta16", 1000.0, DynViscosity::DIMS),
        ("eta26", -1500.0, DynViscosity::DIMS),
        ("maxwell_modulus", modulus, Pressure::DIMS),
        ("maxwell_time", time, Dims([0, 0, 1, 0, 0, 0])),
    ] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: fs_evidence::ValidityDomain::unconstrained().with("T", 280.0, 320.0)
                .with("omega", 0.0, 100000.0), uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: provenance.clone(), observations: vec![observation.clone()],
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "manufactured".into(), revision: 0 }, claims, vec![]).unwrap();
    let normalized = NormalizedPack::new("synthetic-reed", "synthetic-v1", artifact,
        "CC0-1.0", card.claims().clone(), vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(card.id().clone(), normalized).unwrap()
}

#[path = "material_modal/keyed.rs"]
mod keyed;
#[path = "material_modal/forces.rs"]
mod forces;
#[path = "material_modal/relaxation.rs"]
mod relaxation;

#[test]
fn native_material_mimo_cli_matches_original_performer_and_preserves_source_authority() {
    let dir = directory();
    let text = INPUT.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
        .replace("network 2 1 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_section 1 0 0.12 0.007 0.002",
            "network 3 2 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_node junction\nduct_section 1 2 0.06 0.007 0.002\nduct_section 2 0 0.06 0.007 0.002")
        .replace("valve 0.0004\n", "walls 1 20 1000\nwall upstream 2 isotropic blade.matpack 293.15 100 0.0002 stress-free elastic\nvalve 0.0004\n");
    let input = dir.join("source.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    std::fs::write(dir.join("mutual.radiation"), ARTIFACT).unwrap();
    let loaded = MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone()),
        |_| Ok(ARTIFACT.as_bytes().to_vec())).unwrap();
    let hash = loaded.mimo_radiation_records()[0].1.artifact_hash().to_hex();
    assert_eq!(loaded.wall_records().len(), 1);
    assert_eq!(loaded.wall_records()[0].pack(), &material);
    let mut source = loaded.into_renderer();
    let mut pressure = vec![0.0; 512];
    for block in pressure.chunks_mut(37) { source.block(block).unwrap(); }
    let (expected, _) = encode_pcm16_wav(&pressure, 48000, 0.005).unwrap();
    assert!(expected[44..].iter().any(|byte| *byte != 0));
    for block in [1, 37] {
        let output = dir.join(format!("mutual-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        for token in ["\"retained_modes\":6", "\"full_tensor_viscosity\":true",
            "\"radiation_terminals\":2", "\"mimo_radiation_groups\":[",
            "\"side\":\"downstream\"", "\"path\":\"mutual.radiation\"",
            "\"material_duct_walls\":[{\"side\":\"upstream\",\"node\":2", "fixed-state local isotropic membrane",
            "\"ordered_apertures\":[", "signed terminal exchange", "no cross-side exterior coupling"] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
        assert!(metadata.contains(&format!("\"artifact_blake3\":\"{hash}\"")));
        assert!(metadata.contains(&format!("\"declaration_identity\":\"{}\"", "11".repeat(32))));
        assert!(metadata.contains(&format!("\"realization_identity\":\"{}\"", "22".repeat(32))));
    }
    // A stale operator is refused before the CLI creates WAV or sidecar, even
    // though the complete instrument's atmosphere was changed coherently.
    std::fs::write(&input, text.replace("101325", "65000")).unwrap();
    let output = dir.join("stale-atmosphere.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
    std::fs::write(&input, text.replace("mutual.radiation", "missing.radiation")).unwrap();
    let output = dir.join("missing-artifact.wav");
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists());
    assert!(!output.with_extension("provenance.json").exists());
}

#[path = "material_modal/spatial.rs"]
mod spatial;
#[path = "material_modal/wetted.rs"]
mod wetted;
