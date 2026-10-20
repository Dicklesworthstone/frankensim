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

fn pack() -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let artifact = fs_blake3::hash_domain("org.frankensim.test.material-modal-cli", b"manufactured-data");
    let provenance = Provenance { source: "manufactured material-modal CLI regression".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured material-modal valve".into(), method: "synthetic integration fixture".into(),
        artifact, caveats: "not a measured material or instrument".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [
        ("density", 900.0, Density::DIMS), ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE), ("eta11", 22000.0, DynViscosity::DIMS),
        ("eta22", 35000.0, DynViscosity::DIMS), ("eta66", 11000.0, DynViscosity::DIMS),
        ("eta12", 4500.0, DynViscosity::DIMS), ("eta16", 1000.0, DynViscosity::DIMS),
        ("eta26", -1500.0, DynViscosity::DIMS),
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

#[test]
fn native_material_mimo_cli_matches_original_performer_and_preserves_source_authority() {
    let dir = directory();
    let text = INPUT.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005");
    let input = dir.join("source.performance");
    let material = pack();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("blade.matpack"), material.to_bytes()).unwrap();
    std::fs::write(dir.join("mutual.radiation"), ARTIFACT).unwrap();
    let loaded = MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(material.clone()),
        |_| Ok(ARTIFACT.as_bytes().to_vec())).unwrap();
    let hash = loaded.mimo_radiation_records()[0].1.artifact_hash().to_hex();
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
