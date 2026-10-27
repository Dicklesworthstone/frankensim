//! Material-card substitutions through jet feedback, wall storage and the CLI.
use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

#[path = "material/boundaries.rs"]
mod boundaries;
#[path = "material/mimo.rs"]
mod mimo;
#[path = "material/wall_memory.rs"]
mod wall_memory;
#[path = "material/dg.rs"]
mod dg;

fn pack(density: f64, young: f64, omit: Option<&str>, upper_omega: f64)
    -> NormalizedMaterialCardPack
{
    pack_with_time(density, young, omit, upper_omega, 0.001)
}
fn pack_with_time(density: f64, young: f64, omit: Option<&str>, upper_omega: f64, time_s: f64)
    -> NormalizedMaterialCardPack
{
    let mut claims = ClaimSet::new();
    let provenance = Provenance { source: "synthetic compliant cylinder; not measured flute material".into(),
        license: "CC0-1.0".into(), artifact: Some(ContentHash([11; 32])) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic isotropic membrane".into(), method: "manufactured physical regression".into(),
        artifact: ContentHash([11; 32]), caveats: "not a measured material or instrument".into(),
        provenance: provenance.clone(),
    }).unwrap();
    for (name, dims, value) in [("density", Density::DIMS, density),
        ("young_modulus", Pressure::DIMS, young), ("poisson_ratio", Dims::NONE, 0.3),
        ("kelvin_voigt_young_viscosity", DynViscosity::DIMS, 2.0),
        ("boundary_modulus", Pressure::DIMS, 0.5 * young),
        ("boundary_time", Dims([0, 0, 1, 0, 0, 0]), time_s)] {
        if omit == Some(name) { continue; }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 330.0)
                .with("omega", 0.0, upper_omega),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()], provenance: provenance.clone() }).unwrap();
    }
    let claims = NormalizedPack::new("synthetic-wall", "executable-test-v1",
        ContentHash([11; 32]), "synthetic redistribution permitted", claims, vec![], vec![]).unwrap();
    let pack = NormalizedMaterialCardPack::new(MaterialStateId {
        chemistry: "synthetic isotropic wall".into(), phase: "solid".into(),
        process: "manufactured test state".into(), revision: 0,
    }, claims).unwrap();
    NormalizedMaterialCardPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn source() -> String {
    short().replace("tube 0.3 0.006 1e-10 -0.9",
        "network 3 2\nduct_node inlet\nduct_node junction\nduct_node reflection -0.9\nduct_section 0 1 0.15 0.006 5e-11\nduct_section 1 2 0.15 0.006 5e-11")
        .replace("\nflue ", "\nwalls 1 20 1000\nwall 1 isotropic wall.fsmcdpk 293.15 1000 0.0002 stress-free kelvin-voigt\nradiation compact-baffled\nflue ")
        .replace("viscothermal 8 300 500",
            "viscothermal sections 2\nloss_section 0 8 300 500\nloss_section 1 8 300 500")
        .replace("\ncontrols ", "\noutlets 2 0\ncontrols ")
}
fn loaded(text: &str, pack: &NormalizedMaterialCardPack, block: usize) -> FlueFilePerformance {
    FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |name| {
            assert_eq!(name, "wall.fsmcdpk");
            Ok(pack.clone())
        }).unwrap()
}
struct Trace {
    pressure: Vec<f64>,
    incoming: Vec<f64>,
    maximum_wall_energy: f64,
    maximum_terminal_energy: f64,
}
fn trace(p: FlueFilePerformance) -> Trace {
    let mut renderer = p.into_renderer();
    let mut out = Trace { pressure: Vec::new(), incoming: Vec::new(),
        maximum_wall_energy: 0.0, maximum_terminal_energy: 0.0 };
    while renderer.remaining_samples() != 0 {
        let mut one = [0.0]; renderer.block(&mut one).unwrap();
        out.pressure.push(one[0]);
        let network = renderer.voice().network().unwrap();
        let frame = network.last_frame().unwrap();
        out.incoming.push(frame.incoming_pressure_pa);
        let scale = (frame.stored_energy_j + frame.storage_change_j.abs()
            + frame.inlet_work_j.abs() + frame.interior_loss_j + frame.terminal_loss_j).max(1e-30);
        assert!(frame.balance_residual_j().abs() < 1e-10 * scale);
        let wall = network.wall_state(1).unwrap();
        out.maximum_wall_energy = out.maximum_wall_energy.max(wall.stored_energy_j);
        let terminal = network.node_frame(renderer.config().outlet_node).unwrap();
        out.maximum_terminal_energy = out.maximum_terminal_energy.max(terminal.stored_energy_j);
    }
    assert!(out.pressure.iter().all(|p| p.is_finite()));
    assert!(out.pressure.iter().any(|p| p.abs() > 1e-10));
    assert!(out.maximum_wall_energy > 0.0);
    out
}

#[test]
fn a_data_only_material_swap_changes_wall_feedback_and_actual_microphone_pressure() {
    let text = source();
    let original = pack(1000.0, 1e6, None, 100_000.0);
    let p = loaded(&text, &original, 37);
    assert_eq!(p.wall_records()[0].pack(), &original);
    assert_eq!(p.wall_records()[0].node(), 1);
    let n = p.renderer().voice().network().unwrap();
    assert_eq!(n.source().graph.sections.len(), 2);
    assert_eq!(n.material_walls().len(), 1);
    assert_eq!(n.material_walls()[0].length_m(), 0.15);
    assert_eq!(n.material_walls()[0].source().material.card_identity(), original.card().content_hash());
    assert_eq!(n.material_walls()[0].patch().wall.surface_density, 0.2);
    // Wall source addresses and physical coverage survive loss/fractional lowering.
    assert_eq!(n.source().viscothermal.len(), 2);
    assert!(n.graph().sections.len() > n.source().graph.sections.len());
    assert_eq!(n.material_walls()[0].source().node, 1);
    assert_eq!(n.graph().nodes[1], n.material_walls()[0].shunt());
    let reference = trace(p);
    assert!(reference.maximum_terminal_energy > 0.0);
    for changed in [pack(1000.0, 4e6, None, 100_000.0), pack(2000.0, 1e6, None, 100_000.0)] {
        let changed = trace(loaded(&text, &changed, 37));
        assert_ne!(changed.incoming, reference.incoming, "material must change the jet's returning bore wave");
        assert_ne!(changed.pressure, reference.pressure, "material must reach actual exterior pressure");
        assert_ne!(changed.maximum_wall_energy, reference.maximum_wall_energy);
    }
}

#[test]
fn exterior_selection_drives_terminal_reaction_and_receiver_through_the_same_medium() {
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let base = source();
    let warm = base.replace("exterior same", "exterior 323.15 101325 0");
    let p = loaded(&warm, &material, 37);
    let renderer = p.renderer();
    let exterior = renderer.radiation_input().unwrap().exterior;
    let binding = &p.radiation_bindings()[0];
    assert_eq!(binding.opening().node, renderer.config().outlet_node);
    assert_eq!(binding.load().density_kg_m3(), exterior.density);
    assert_eq!(binding.load().sound_speed_m_s(), exterior.sound_speed);
    assert_eq!(binding.interior_gas(), &renderer.voice().network().unwrap().source().gases[1]);
    assert_ne!(binding.interior_gas().density, exterior.density);
    assert_eq!(renderer.voice().network().unwrap().graph().nodes[binding.opening().node],
        binding.load().termination());
    let a = trace(loaded(&base, &material, 37));
    let b = trace(p);
    assert_ne!(a.incoming, b.incoming, "exterior material changes the coupled source when reaction is selected");
    assert_ne!(a.pressure, b.pressure);
    // The legacy prescribed reflection deliberately leaves exterior reaction out.
    let a = trace(loaded(&base.replace("compact-baffled", "prescribed"), &material, 37));
    let b = trace(loaded(&warm.replace("compact-baffled", "prescribed"), &material, 37));
    assert_eq!(a.incoming, b.incoming);
    assert_eq!(a.maximum_wall_energy, b.maximum_wall_energy);
    assert_ne!(a.pressure, b.pressure);
}

#[test]
fn normalized_material_pack_reaches_the_cli_and_preserves_pressure_across_callbacks() {
    let dir = directory();
    let input = dir.join("material-flue.performance");
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let text = source();
    std::fs::write(&input, &text).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
    let reference = trace(loaded(&text, &material, 37));
    let expected = encode_pcm16_wav(&reference.pressure, 48_000, 1.0).unwrap().0;
    for block in [1, 37, 512] {
        let output = dir.join(format!("material-flue-{block}.wav"));
        success(run(&input, &output, block, false));
        assert_eq!(std::fs::read(&output).unwrap(), expected);
        let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
        let pack_hash = material.content_hash().to_hex();
        let card_hash = material.card().content_hash().to_hex();
        for token in ["isotropic-local-ring-membrane", "\"wall_strips\":1",
            "\"covered_length_m\":1.5e-1", "\"radiation_feedback_added\":true",
            "\"kind\":\"relaxation\"", "compact-baffled-piston-positive-real-v1",
            "\"composition_authority\":\"EstimateOnly\"", pack_hash.as_str(), card_hash.as_str()] {
            assert!(metadata.contains(token), "{token}: {metadata}");
        }
    }
    // The performance bytes and all mechanics stay fixed; only the pack changes.
    let substituted = pack(2000.0, 4e6, None, 100_000.0);
    std::fs::write(dir.join("wall.fsmcdpk"), substituted.to_bytes()).unwrap();
    let output = dir.join("substituted.wav");
    success(run(&input, &output, 37, false));
    let expected = trace(loaded(&text, &substituted, 37));
    assert_ne!(expected.pressure, reference.pressure);
    assert_eq!(std::fs::read(&output).unwrap(),
        encode_pcm16_wav(&expected.pressure, 48_000, 1.0).unwrap().0);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains(&substituted.content_hash().to_hex()));
}

#[test]
fn unsupported_walls_missing_sources_and_out_of_band_reaction_refuse_before_output() {
    let text = source();
    let gate = CancelGate::new_clock_free();
    assert!(FlueFilePerformance::from_bytes(text.as_bytes(), 37, &gate).is_err());
    let material = pack(1000.0, 1e6, None, 100_000.0);
    for changed in [text.replace("wall.fsmcdpk", "../wall.fsmcdpk"),
        text.replace("0.0002 stress-free", "0.001 stress-free"),
        text.replace("walls 1 20 1000", "walls 65 20 1000"),
        text.replace("walls 1 20 1000", "walls 1 20 500"),
        text.replace("kelvin-voigt", "guess"),
        text.replace("exterior same", "exterior 293.15 80000 0"),
        text.replace("fractional 1000", "fractional 6000")
            .replace("16 1000", "16 6000").replace("walls 1 20 1000", "walls 1 20 6000")] {
        assert!(FlueFilePerformance::from_bytes_with_sources(changed.as_bytes(), 37, &gate,
            None, |_| Ok(material.clone())).is_err(), "accepted {changed}");
    }
    for missing in [pack(1000.0, 1e6, Some("poisson_ratio"), 100_000.0),
        pack(1000.0, 1e6, Some("kelvin_voigt_young_viscosity"), 100_000.0),
        pack(1000.0, 1e6, None, 2000.0)] {
        assert!(FlueFilePerformance::from_bytes_with_sources(text.as_bytes(), 37, &gate,
            None, |_| Ok(missing.clone())).is_err());
    }
    let elastic = pack(1000.0, 1e6, Some("kelvin_voigt_young_viscosity"), 100_000.0);
    assert!(FlueFilePerformance::from_bytes_with_sources(text.replace("kelvin-voigt", "elastic").as_bytes(),
        37, &gate, None, |_| Ok(elastic.clone())).is_ok());
    let legacy = short();
    assert_eq!(physical(&legacy), physical(&legacy.replace("\nflue ", "\nwalls none\nradiation prescribed\nflue ")));
    let dir = directory(); let input = dir.join("missing-material.performance");
    let output = dir.join("refused.wav"); std::fs::write(&input, text).unwrap();
    assert!(!run(&input, &output, 37, false).status.success());
    assert!(!output.exists()); assert!(!output.with_extension("provenance.json").exists());
}

#[test]
fn walls_and_reaction_preserve_branched_geometry_section_gases_and_all_outlets() {
    let branched = source().replace(
        "network 3 2\nduct_node inlet\nduct_node junction\nduct_node reflection -0.9\nduct_section 0 1 0.15 0.006 5e-11\nduct_section 1 2 0.15 0.006 5e-11",
        "network 5 4\nduct_node inlet\nduct_node junction\nduct_node junction\nduct_node reflection -0.9\nduct_node reflection -0.7\nduct_section 0 1 0.08 0.006 5e-11\nduct_section 1 2 0.08 0.006 5e-11\nduct_section 2 3 0.15 0.006 5e-11\nduct_section 2 4 0.05 0.003 5e-11 gas 320 101325 0")
        .replace("outlets 2 0", "outlets 3 1\noutlet 4 0.2 0")
        // The 0.08 m sections at 96 kHz admit at most 5 loss cells: each of the 4*CELLS
        // quarters must span one transit sample (343 m/s / 96 kHz = 3.6 mm).
        .replace("loss_section 0 8 300 500\nloss_section 1 8 300 500",
            "loss_section 0 4 300 500\nloss_section 1 4 300 500");
    assert_eq!(branched.matches("loss_section 0 4 300 500").count(), 1);
    let material = pack(1000.0, 1e6, None, 100_000.0);
    let p = loaded(&branched, &material, 37);
    assert_eq!(p.radiation_bindings().len(), 2);
    assert_eq!(p.renderer().outlet_nodes(), &[3,4]);
    let n = p.renderer().voice().network().unwrap();
    assert_eq!(n.source().graph.sections.len(), 4);
    assert_eq!(n.source().gases[3].temperature, 320.0);
    assert_eq!(n.material_walls()[0].length_m(), 0.08);
    assert_eq!(p.radiation_bindings()[1].interior_gas(), &n.source().gases[3]);
    let a = trace(p);
    let b = trace(loaded(&branched, &pack(1500.0, 4e6, None, 100_000.0), 37));
    assert_ne!(a.incoming, b.incoming);
    assert_ne!(a.pressure, b.pressure);
    let dir = directory(); let input = dir.join("branched-material.performance");
    let output = dir.join("branched-material.wav");
    std::fs::write(&input, branched).unwrap();
    std::fs::write(dir.join("wall.fsmcdpk"), material.to_bytes()).unwrap();
    success(run(&input, &output, 37, false));
    assert_eq!(std::fs::read(&output).unwrap(), encode_pcm16_wav(&a.pressure, 48_000, 1.0).unwrap().0);
    let metadata = std::fs::read_to_string(output.with_extension("provenance.json")).unwrap();
    assert!(metadata.contains("\"node\":4,\"original_reflection\":"));
}
