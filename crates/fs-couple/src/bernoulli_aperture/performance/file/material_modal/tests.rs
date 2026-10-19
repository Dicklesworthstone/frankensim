use super::*;
use crate::bernoulli_aperture::two_sided::{LoadedAperture, ValveSide};
use crate::pcm_wav::{encode_pcm16_wav, observation::PressureRenderer};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedPack, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const INPUT: &str = include_str!("../../../../../examples/material_modal_valve/paired.performance");

fn pack(rho: f64, missing: bool) -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", rho, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("eta11", 22000.0, DynViscosity::DIMS),
        ("eta22", 35000.0, DynViscosity::DIMS),
        ("eta66", 11000.0, DynViscosity::DIMS),
        ("eta12", 4500.0, DynViscosity::DIMS),
        ("eta16", 1000.0, DynViscosity::DIMS),
        ("eta26", -1500.0, DynViscosity::DIMS),
    ] {
        if missing && name == "eta26" { continue; }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0).with("omega", 0.0, 100000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "manufactured material-modal fixture".into(),
                license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(),
        process: "manufactured".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let normalized = NormalizedPack::new("synthetic-reed", "synthetic-v1",
        hash_domain("org.frankensim.test.material-modal", &rho.to_le_bytes()),
        "CC0-1.0", card.claims().clone(), vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(card.id().clone(), normalized).unwrap()
}

fn input() -> String {
    INPUT.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
}
fn build(text: &str, block: usize) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |path| {
            Ok(pack(if path == "altered.matpack" { 1150.0 } else { 900.0 }, false))
        }).unwrap()
}
fn valve(p: &MaterialModalValvePerformance) -> &ModalAperture {
    match p.renderer().source().system().valve() {
        LoadedAperture::Modal(a) => a,
        _ => panic!("the file must retain its actual multi-mode valve"),
    }
}
fn render(p: MaterialModalValvePerformance) -> (Vec<f64>, ScheduledTwoSidedPressure) {
    let samples = p.info().samples as usize;
    let mut source = p.into_renderer();
    let mut pressure = vec![0.0; samples];
    for block in pressure.chunks_mut(source.max_block_len()) { source.block(block).unwrap(); }
    (pressure, source)
}

#[test]
fn g0_original_material_modes_gases_loss_sections_and_radiation_are_retained() {
    let mut loads = 0;
    // Manufactured very short passage: the original full loaded-clock gate
    // must still admit its inertance; ordinary long supplies may require a
    // different coupled discretization and are never silently weakened here.
    let supplied = input().replace("supply none", "supply slit 0.00002 0.01 0.001 16");
    let p = MaterialModalValvePerformance::from_bytes_with_sources(supplied.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        }).unwrap();
    assert_eq!(loads, 1); // both blades and every region share one original artifact
    assert_eq!((p.info().blades, p.info().modes, p.info().nodes, p.info().triangles), (2, 6, 42, 48));
    assert_eq!(p.material_records().len(), 4);
    let a = valve(&p);
    assert_eq!(a.material_plates()[0].mode_range(), 0..3);
    assert_eq!(a.material_plates()[1].mode_range(), 3..6);
    assert_eq!(a.viscosity_matrix_n_s_m().unwrap().len(), 36);
    for blade in a.material_plates() {
        assert_eq!(blade.source().regions().len(), 2);
        assert!(blade.viscosity().is_some());
        assert!(blade.source().bending_viscosity().is_none());
        assert_eq!(blade.displacement_shapes().len(), 21);
    }
    let system = p.renderer().source().system();
    assert_eq!(system.source(ValveSide::Upstream).graph.inlet_section_index().unwrap(), 0);
    assert_eq!(system.source(ValveSide::Upstream).graph.sections[0].nodes, [1, 0]);
    assert_eq!(system.source(ValveSide::Downstream).viscothermal.len(), 1);
    assert_eq!(system.source(ValveSide::Downstream).viscothermal[0].section, 1);
    assert_eq!(system.source(ValveSide::Downstream).gases[1].temperature, 303.15);
    assert!(system.fractional_sections(ValveSide::Downstream).is_some());
    assert_eq!(p.radiation_bindings().len(), 1);
    assert_eq!(p.radiation_bindings()[0].opening().node, 2);
    assert_eq!(p.renderer().source().mechanics_per_output(), 16);
    assert_eq!(p.renderer().source().mechanical_samples_rendered(), 0);
    assert_eq!(system.supply_channel().unwrap().spec().length_m, 0.00002);
    assert_eq!(system.supply_channel().unwrap().gas(), &system.source(ValveSide::Upstream).gases[0]);
    let mut renderer = p.into_renderer();
    renderer.block(&mut [0.0; 8]).unwrap();
    assert_eq!(renderer.source().mechanical_samples_rendered(), 128);
}

#[test]
fn g3_data_only_material_and_pressure_changes_reach_the_actual_microphone() {
    let text = input();
    let changed = text.replacen("section blade.matpack", "section altered.matpack", 1);
    let original = build(&text, 37);
    let changed = build(&changed, 37);
    assert_ne!(valve(&original).material_plates()[0].source().specimen_identity(),
        valve(&changed).material_plates()[0].source().specimen_identity());
    assert_eq!(valve(&original).material_plates()[1].source().specimen_identity(),
        valve(&changed).material_plates()[1].source().specimen_identity());
    let (pressure, original) = render(original);
    let (material, _) = render(changed);
    let (low_pressure, _) = render(build(&text.replace("101325", "65000"), 37));
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert!(pressure.iter().zip(&material).any(|(a,b)| (a-b).abs() > 1e-15));
    assert!(pressure.iter().zip(&low_pressure).any(|(a,b)| (a-b).abs() > 1e-15));
    assert_eq!(original.source().mechanical_samples_rendered(), 512 * 16);
    assert!(original.pending_controls().is_empty());
    assert!(original.source().system().stored_energy_j() > 0.0);
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    assert_eq!(encode_pcm16_wav(&pressure, 48000, 0.005).unwrap(),
        encode_pcm16_wav(&partitioned, 48000, 0.005).unwrap());
}

#[test]
fn g0_missing_material_coefficients_and_incomplete_physical_inputs_refuse() {
    let gate = CancelGate::new_clock_free();
    let text = input();
    let missing = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &gate, None, |_| Ok(pack(900.0, true))).err().unwrap();
    assert!(matches!(missing, PlateValveInputError::Material { .. }));
    for altered in [
        text.replacen("blades 2", "blades 3", 1),
        text.replacen("blade 3", "blade 33", 1),
        text.replacen("768000 48000", "800000 48000", 1),
        text.replacen("blade.matpack", "../blade.matpack", 1),
        text.replacen("duct_node reflection 0.7", "duct_node baffled-low-ka 1000", 1),
        text.replacen("303.15 101325", "303.15 65000", 1),
    ] {
        assert!(MaterialModalValvePerformance::from_bytes_with_sources(altered.as_bytes(), 37,
            &gate, None, |_| Ok(pack(900.0, false))).is_err());
    }
    gate.request();
    let mut requested = false;
    assert!(MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &gate, None, |_| { requested = true; Ok(pack(900.0, false)) }).is_err());
    assert!(!requested);
}
