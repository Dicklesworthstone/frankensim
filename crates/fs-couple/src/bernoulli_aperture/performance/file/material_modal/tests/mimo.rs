use super::*;

const PERFORMANCE: &str = include_str!("../../../../../../examples/material_modal_valve/paired-mutual.performance");
const OPERATOR: &str = include_str!("../../../../../../examples/material_modal_valve/mutual.radiation");

fn mimo_input() -> String {
    PERFORMANCE.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
}
fn build_mimo(text: &str, artifact: &str, block: usize)
    -> Result<MaterialModalValvePerformance, PlateValveInputError>
{
    MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None,
        |path| Ok(pack(if path == "altered.matpack" { 1150.0 } else { 900.0 }, false)),
        |path| { assert_eq!(path, "mutual.radiation"); Ok(artifact.as_bytes().to_vec()) })
}

#[test]
fn native_mimo_retains_original_material_modes_gases_order_and_unselected_side_loads() {
    let text = mimo_input().replace("supply none", "supply slit 0.00002 0.01 0.001 16");
    let loaded = build_mimo(&text, OPERATOR, 37).unwrap();
    assert!(loaded.radiation_bindings().is_empty());
    assert_eq!(loaded.mimo_radiation_records().len(), 1);
    let (side, record) = &loaded.mimo_radiation_records()[0];
    assert_eq!(*side, ValveSide::Downstream);
    assert_eq!(record.path(), "mutual.radiation");
    assert_eq!(record.declaration_identity(), [0x11; 32]);
    assert_eq!(record.realization_identity(), [0x22; 32]);
    assert_ne!(record.artifact_hash().0, record.realization_identity());
    assert_eq!(valve(&loaded).modal_states().len(), 6);
    assert_eq!(valve(&loaded).viscosity_matrix_n_s_m().unwrap().len(), 36);
    assert_eq!(loaded.material_records().len(), 4);
    let system = loaded.renderer().source().system();
    assert_eq!(system.source(ValveSide::Upstream).graph.nodes[0],
        crate::bernoulli_aperture::network::NetworkNode::Termination { reflection: 0.7 });
    assert!(system.mimo_environment(ValveSide::Upstream).is_none());
    let environment = system.mimo_environment(ValveSide::Downstream).unwrap();
    assert_eq!(environment.source().apertures, record.apertures());
    assert_eq!(environment.source().exterior, *record.gas_record().state());
    assert_eq!(environment.radiation().openings.iter().map(|o| (o.network,o.node)).collect::<Vec<_>>(),
        vec![(0,2),(0,3)]);
    assert_eq!(system.source(ValveSide::Downstream).gases[1].temperature, 303.15);
    assert_eq!(system.viscothermal_sections(ValveSide::Downstream).len(), 1);
    assert!(system.fractional_sections(ValveSide::Downstream).is_some());
    assert_eq!(system.supply_channel().unwrap().spec().length_m, 0.00002);
    assert_eq!(loaded.renderer().source().mechanical_samples_rendered(), 0);
}

#[test]
fn mutual_terms_and_material_substitution_reach_native_pcm_with_callback_parity() {
    let text = mimo_input();
    let (pressure, rendered) = render(build_mimo(&text, OPERATOR, 37).unwrap());
    let (partitioned, other) = render(build_mimo(&text, OPERATOR, 1).unwrap());
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(pressure, partitioned);
    assert_eq!(encode_pcm16_wav(&pressure, 48000, 0.005).unwrap(),
        encode_pcm16_wav(&partitioned, 48000, 0.005).unwrap());
    let states = rendered.source().system().mimo_pole_states(ValveSide::Downstream).unwrap();
    assert!(states.iter().any(|state| state.re != 0.0 || state.im != 0.0));
    assert_eq!(states, other.source().system().mimo_pole_states(ValveSide::Downstream).unwrap());
    assert_eq!(rendered.source().mechanical_samples_rendered(), 512 * 16);
    assert!(rendered.pending_controls().is_empty());
    let independent = OPERATOR.replace("1000000 0 1000000 0", "0 0 0 0")
        .replace("direct 2000 300 300 2000", "direct 2000 0 0 2000");
    let (uncoupled, _) = render(build_mimo(&text, &independent, 37).unwrap());
    assert!(pressure.iter().zip(uncoupled).any(|(a,b)| (a-b).abs() > 1e-15));
    let material = text.replacen("section blade.matpack", "section altered.matpack", 1);
    let (changed_material, _) = render(build_mimo(&material, OPERATOR, 37).unwrap());
    assert!(pressure.iter().zip(changed_material).any(|(a,b)| (a-b).abs() > 1e-15));
}

#[test]
fn changed_atmosphere_requires_its_own_independently_encoded_operator_gas() {
    let text = mimo_input().replace("101325", "65000");
    assert!(matches!(build_mimo(&text, OPERATOR, 37), Err(PlateValveInputError::Physics(_))));
    // This second manufactured artifact explicitly declares its own gas. It is
    // integration data, not an automatically refitted or calibrated operator.
    let artifact = OPERATOR.replace("101325", "65000");
    let loaded = build_mimo(&text, &artifact, 37).unwrap();
    let record = &loaded.mimo_radiation_records()[0].1;
    assert_eq!(record.gas_record().state().pressure, 65000.0);
    assert_eq!(loaded.renderer().source().system().mimo_environment(ValveSide::Downstream)
        .unwrap().source().exterior.pressure, 65000.0);
    let (low, _) = render(loaded);
    let (standard, _) = render(build_mimo(&mimo_input(), OPERATOR, 37).unwrap());
    assert!(low.iter().zip(standard).any(|(a,b)| (a-b).abs() > 1e-15));
}

#[test]
fn native_selection_reuses_runtime_geometry_band_passivity_and_resource_refusals() {
    for artifact in [
        OPERATOR.replace("aperture 0.03 0", "aperture 0.04 0"),
        OPERATOR.replace("0.007", "0.008"),
        OPERATOR.replace("reciprocal true", "reciprocal false"),
        OPERATOR.replace("direct 2000 300 300 2000", "direct 2000 300 301 2000"),
        OPERATOR.replace("direct 2000 300 300 2000", "direct -2000 0 0 -2000"),
        OPERATOR.replace("improper 0 0 0 0", "improper 1 0 0 1"),
        OPERATOR.replace("pole real -7000", "pole real 7000"),
        OPERATOR.replace("band 1 1000000", "band 1 1000"),
    ] {
        assert!(matches!(build_mimo(&mimo_input(), &artifact, 37), Err(PlateValveInputError::Physics(_))));
    }
    for text in [
        mimo_input().replace("256 8388608 65", "1 8388608 65"),
        mimo_input().replace("256 8388608 65", "256 1 65"),
        mimo_input().replace("opening downstream 3 0.03 0", "opening downstream 3 0.04 0"),
    ] { assert!(build_mimo(&text, OPERATOR, 37).is_err()); }
}

#[test]
fn missing_artifact_invalid_group_and_cancellation_never_return_a_partial_performer() {
    let text = mimo_input();
    assert!(matches!(MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false))),
        Err(PlateValveInputError::RadiationSource { .. })));
    for altered in [
        text.replace("radiation mimo 1", "radiation mimo 0"),
        text.replace("radiation_group downstream", "radiation_group upstream"),
        text.replace("mutual.radiation", "../mutual.radiation"),
        text.replace("radiation mimo 1", "radiation mimo 2\nradiation_group downstream mutual.radiation 256 8388608 65"),
    ] {
        assert!(MaterialModalValvePerformance::from_bytes_with_all_sources(altered.as_bytes(), 37,
            &CancelGate::new_clock_free(), None,
            |_| panic!("invalid radiation selector must refuse before material loading"),
            |_| panic!("invalid radiation selector must refuse before artifact loading")).is_err());
    }
    let gate = CancelGate::new_clock_free();
    let refused = MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), 37,
        &gate, None, |_| Ok(pack(900.0, false)), |_| {
            gate.request(); Ok(OPERATOR.as_bytes().to_vec())
        });
    assert!(matches!(refused, Err(PlateValveInputError::RadiationSource { .. })));
}

#[test]
fn explicit_compact_selector_and_legacy_omission_preserve_identical_pcm() {
    let legacy = input();
    let explicit = legacy.replacen("valve 0.0004", "radiation compact-baffled\nvalve 0.0004", 1);
    let compact = MaterialModalValvePerformance::from_bytes_with_all_sources(explicit.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)),
        |_| panic!("compact radiation never loads an external MIMO source")).unwrap();
    assert!(compact.mimo_radiation_records().is_empty());
    assert_eq!(render(compact).0, render(build(&legacy, 1)).0);
}
