use super::*;
use crate::bernoulli_aperture::two_sided::{LoadedAperture, ValveSide};
use crate::pcm_wav::{encode_pcm16_wav, observation::PressureRenderer};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const INPUT: &str = include_str!("../../../../../examples/material_modal_valve/paired.performance");

mod mimo;
mod boundaries;
mod keyed;
mod walls;
mod forces;
mod pressure;
mod wetted;
mod relaxation;
mod thickness;
mod solid;
mod observer_gas;

fn pack(rho: f64, missing: bool) -> NormalizedMaterialCardPack {
    pack_with_relaxation(rho, missing, 1e8, 0.001)
}
fn pack_with_relaxation(rho: f64, missing: bool, modulus: f64, time: f64) -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let identity = [rho.to_le_bytes(), modulus.to_le_bytes(), time.to_le_bytes()].concat();
    let artifact = hash_domain("org.frankensim.test.material-modal", &identity);
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured material-modal reed".into(),
        method: "manufactured physical regression".into(), artifact,
        caveats: "not a measured material".into(),
        provenance: Provenance { source: "manufactured material-modal fixture".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
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
        ("boundary_modulus", modulus, Pressure::DIMS),
        ("boundary_time", time, Dims([0, 0, 1, 0, 0, 0])),
        ("fast_modulus", 3e7, Pressure::DIMS),
        ("fast_time", 0.0001, Dims([0, 0, 1, 0, 0, 0])),
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
            observations: vec![observation.clone()],
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(),
        process: "manufactured".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let normalized = NormalizedPack::new("synthetic-reed", "synthetic-v1", artifact,
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
    // The exact structural image still retains the original passage inertance
    // and its separate no-slip flow state; selecting a clock never edits geometry.
    let supplied = input().replace("supply none", "supply slit 0.00002 0.01 0.001 16");
    let p = MaterialModalValvePerformance::from_bytes_with_sources(supplied.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        }).unwrap();
    assert_eq!(loads, 1); // both blades and every region share one original artifact
    assert_eq!(p.info().propagator, BladePropagator::ExactZoh);
    assert_eq!((p.info().blades, p.info().modes, p.info().nodes, p.info().triangles), (2, 6, 42, 48));
    assert_eq!(p.material_records().len(), 4);
    let a = valve(&p);
    assert_eq!(a.blade_propagator(), BladePropagator::ExactZoh);
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
    assert_eq!(system.supply_channel().unwrap().slit().unwrap().spec().length_m, 0.00002);
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
fn g3_circular_supply_geometry_gas_and_reed_material_reach_the_same_microphone() {
    let text = input().replace("supply none", "supply circular 0.02 0.0005 32");
    let source = build(&text, 37);
    let system = source.renderer().source().system();
    let channel = system.supply_channel().unwrap().circular().unwrap();
    assert_eq!(channel.spec().radial_cells, 32);
    assert_eq!(channel.spec().length_m, 0.02);
    assert_eq!(channel.spec().radius_m, 0.0005);
    assert_eq!(channel.gas(), &system.source(ValveSide::Upstream).gases[0]);
    let (pressure, source) = render(source);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    let channel = source.source().system().supply_channel().unwrap();
    assert!(channel.last_step().kinetic_energy_j > 0.0);
    assert!(channel.last_step().viscous_dissipation_j > 0.0);
    assert!(channel.velocities_m_s().windows(2).any(|v| v[0] != v[1]));
    let states = |source: &ScheduledTwoSidedPressure| match source.source().system().valve() {
        LoadedAperture::Modal(valve) => valve.modal_states().to_vec(),
        _ => panic!("original modal material source must remain active"),
    };
    for changed in [text.replace("0.02 0.0005 32", "0.02 0.0007 32"),
        text.replace("101325", "65000"),
        text.replacen("section blade.matpack", "section altered.matpack", 1)] {
        let (altered_pressure, altered) = render(build(&changed, 37));
        assert!(pressure.iter().zip(&altered_pressure).any(|(a, b)| (a - b).abs() > 1e-15));
        assert_ne!(states(&source), states(&altered),
            "geometry, atmosphere and material must alter actual reed mechanics, not only audio gain");
    }
    let (partitioned, partitioned_source) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    assert_eq!(states(&source), states(&partitioned_source));
    assert_eq!(channel.last_step(), partitioned_source.source().system().supply_channel().unwrap().last_step());
}

#[test]
fn g0_circular_supply_malformed_geometry_refuses_before_material_loading() {
    for supply in ["circular 0.02 -0.0005 32", "circular 0.02 0.0005 1",
        "circular 0.02 0.0005 4097", "circular 0.02 0.0005", "circular NaN 0.0005 32"] {
        let text = input().replace("supply none", &format!("supply {supply}"));
        assert!(MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
            &CancelGate::new_clock_free(), None,
            |_| panic!("invalid supply geometry must refuse before material loading")).is_err());
    }
}

#[test]
fn g3_omitted_and_explicit_midpoint_keep_identical_original_physics() {
    let legacy = input().replace("propagator exact-zoh\n", "");
    let explicit = legacy.replacen("ambient ", "propagator midpoint\nambient ", 1);
    let a = build(&legacy, 37);
    let b = build(&explicit, 1);
    assert_eq!(a.info().propagator, BladePropagator::Midpoint);
    assert_eq!(b.info().propagator, BladePropagator::Midpoint);
    assert_eq!(valve(&a).spec(), valve(&b).spec());
    assert_eq!(valve(&a).viscosity_matrix_n_s_m(), valve(&b).viscosity_matrix_n_s_m());
    let (pa, a) = render(a);
    let (pb, b) = render(b);
    assert_eq!(pa, pb);
    let states = |source: &ScheduledTwoSidedPressure| match source.source().system().valve() {
        LoadedAperture::Modal(a) => a.modal_states().to_vec(),
        _ => panic!("original modal source"),
    };
    assert_eq!(states(&a), states(&b));
}

#[test]
fn g3_file_selected_exact_structure_admits_high_modes_and_materials_reach_the_microphone() {
    // A deliberately coarser numerical rung exercises the selector. This is
    // not a claim that nonlinear contact/channel time error converged at 96 kHz.
    let text = input().replace("audio 768000 48000", "audio 96000 48000")
        .replace("blade 3 0.2 3000 0.04", "blade 3 0.2 3000 0.25");
    let midpoint = text.replace("propagator exact-zoh", "propagator midpoint");
    let refused = MaterialModalValvePerformance::from_bytes_with_sources(midpoint.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false))).err().unwrap();
    assert!(matches!(refused, PlateValveInputError::Physics(
        crate::acoustic_realize::AcousticRealizeError::TimeResolution(_))));
    let source = build(&text, 37);
    assert_eq!(source.info().mechanics_rate_hz, 96000);
    assert_eq!(source.info().propagator, BladePropagator::ExactZoh);
    assert_eq!(valve(&source).modal_states().len(), 6);
    assert_eq!(valve(&source).viscosity_matrix_n_s_m().unwrap().len(), 36);
    let (pressure, rendered) = render(source);
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert_eq!(rendered.source().mechanical_samples_rendered(), 512 * 2);
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
    let changed = text.replacen("section blade.matpack", "section altered.matpack", 1);
    let (material, _) = render(build(&changed, 37));
    assert!(pressure.iter().zip(material).any(|(a,b)| (a-b).abs() > 1e-15));
    let elastic = text.replace("loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000", "loss elastic");
    let elastic = build(&elastic, 37);
    assert_eq!(valve(&elastic).blade_propagator(), BladePropagator::ExactZoh);
    assert!(valve(&elastic).viscosity_matrix_n_s_m().is_none());
}

#[test]
fn g0_invalid_structural_selector_refuses_before_loading_materials() {
    for token in ["unknown", "frequency-matched", "exact-zoh extra", ""] {
        let text = input().replace("propagator exact-zoh", &format!("propagator {token}"));
        let mut loads = 0;
        let error = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |_| { loads += 1; Ok(pack(900.0, false)) })
            .err().unwrap();
        assert!(matches!(error, PlateValveInputError::Input { line: 3, .. }));
        assert_eq!(loads, 0);
    }
    let duplicated = input().replace("propagator exact-zoh\n",
        "propagator exact-zoh\npropagator midpoint\n");
    assert!(MaterialModalValvePerformance::from_bytes_with_sources(duplicated.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |_| panic!("duplicate selector must refuse before sources"))
        .is_err());
}

#[test]
fn g0_missing_material_coefficients_and_incomplete_physical_inputs_refuse() {
    let gate = CancelGate::new_clock_free();
    let text = input();
    let missing = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &gate, None, |_| Ok(pack(900.0, true))).err().unwrap();
    assert!(matches!(missing, PlateValveInputError::Material { .. }));
    // A section mixture record carries its own pressure (the legacy section
    // record inherits it); at the shared pressure it is admitted, so the
    // 65000 Pa refusal below is the static-pressure check.
    let mixture = |pa: &str| text.replacen("gas 303.15 0 ",
        &format!("gas mixture 303.15 {pa} mole 2 N2 0.79 O2 0.21 "), 1);
    assert!(MaterialModalValvePerformance::from_bytes_with_sources(mixture("101325").as_bytes(),
        37, &gate, None, |_| Ok(pack(900.0, false))).is_ok());
    for altered in [
        text.replacen("blades 2", "blades 3", 1),
        text.replacen("blade 3", "blade 33", 1),
        text.replacen("768000 48000", "800000 48000", 1),
        text.replacen("blade.matpack", "../blade.matpack", 1),
        text.replacen("duct_node reflection 0.7", "duct_node baffled-low-ka 1000", 1),
        text.replacen("exterior 293.15 101325", "exterior 293.15 65000", 1),
        mixture("65000"),
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

mod spatial;


#[test]
fn g3_local_opening_gas_and_geometry_reach_modal_valve_and_microphone() {
    // Explicit manufactured closure, not a measured instrument coefficient.
    let reference = input();
    let text = reference.replace("valve 0.0004",
        "opening_losses 1\nopening_loss downstream 2 0.6 0.2 authored-comparison\nvalve 0.0004");
    assert_ne!(reference, text);
    let source = build(&text, 37);
    let system = source.renderer().source().system();
    let binding = &system.opening_losses(ValveSide::Downstream)[0];
    assert_eq!(system.opening_losses(ValveSide::Upstream).len(), 0);
    assert_eq!(binding.input().node, 2);
    assert_eq!(binding.radius_m(), 0.007);
    assert_eq!(binding.gas(), &system.source(ValveSide::Downstream).gases[1]);
    assert_eq!(binding.gas().temperature, 303.15);
    let coefficient = binding.loss().coefficient_pa_s2_m6();
    let density = binding.gas().density;
    let mixed = text.replace("gas 303.15 0 viscothermal",
        "gas mixture 303.15 101325 mole 2 N2 0.6 O2 0.4 viscothermal");
    assert_ne!(mixed, text);
    let mixture = build(&mixed, 37);
    let mixed_binding = &mixture.renderer().source().system()
        .opening_losses(ValveSide::Downstream)[0];
    assert_ne!(mixed_binding.gas().density, density);
    assert!((mixed_binding.loss().coefficient_pa_s2_m6()/coefficient
        - mixed_binding.gas().density/density).abs() < 2e-12);
    let narrow = text.replace("duct_section 1 2 0.15 0.007",
        "duct_section 1 2 0.15 0.006");
    assert_ne!(narrow, text);
    let narrow = build(&narrow, 37);
    let narrow_binding = &narrow.renderer().source().system()
        .opening_losses(ValveSide::Downstream)[0];
    assert!((narrow_binding.loss().coefficient_pa_s2_m6()/coefficient
        - (0.007_f64/0.006).powi(4)).abs() < 2e-12);
    let (pressure, rendered) = render(source);
    let node = rendered.source().system().node_frame(ValveSide::Downstream, 2).unwrap();
    assert!(node.quadratic_loss_j > 0.0);
    assert!(node.absorbed_energy_j >= node.quadratic_loss_j);
    for alternative in [build(&reference, 37), mixture, narrow] {
        let (altered, _) = render(alternative);
        assert!(pressure.iter().zip(&altered).any(|(a,b)| a != b),
            "declared opening physics must reach the actual microphone");
    }
    let (partitioned, _) = render(build(&text, 1));
    assert_eq!(pressure, partitioned);
}

mod dg;

// Manufactured source laws exercise the original material query, not an
// empirical wet-cane fit or an inferred equilibrium with the surrounding gas.
fn environment_pack(base: NormalizedMaterialCardPack, typed: bool) -> NormalizedMaterialCardPack {
    use fs_qty::{QuantitySpec, semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm}};
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static));
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "manufactured held-moisture material".into(),
        method: "declared source interpolation".into(), artifact:base.source_artifact(),
        caveats: "not measured cane, wall tissue or an instrument".into(),
        provenance: Provenance {source:"manufactured environmental law".into(),
            license:"CC0-1.0".into(), artifact:None},
    }).unwrap();
    for (_, original) in base.card().claims().claims_ordered() {
        let PropertyValue::Scalar {value, dims} = &original.value else { panic!("scalar test source") };
        let name = original.key.name();
        let factor = if name == "density" { 1.1 }
            else if name == "young_modulus" || name.contains("modulus") || name.starts_with("branch_") { 0.8 }
            else if name.starts_with("eta") || name.contains("viscosity") { 1.5 }
            else if name.ends_with("time") { 1.2 } else { 1.0 };
        let mut claim = original.clone();
        claim.value = PropertyValue::Curve {abscissa:"moisture_content".into(),
            abscissa_dims:Dims::NONE, knots:vec![(0.0,*value),(0.2,*value*factor)], dims:*dims};
        claim.validity = if typed {
            claim.validity.with_quantity("moisture_content", moisture, 0.0, 0.2)
        } else { claim.validity.with("moisture_content", 0.0, 0.2) };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claim.observations = vec![observation.clone()];
        claims.insert_claim(claim).unwrap();
    }
    let normalized = NormalizedPack::new("held-environment-fixture", "manufactured-v1",
        base.source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(base.card().id().clone(), normalized).unwrap()
}
fn at_environment(text: &str, prefix: &str, moisture: f64) -> String {
    text.lines().map(|line| if line.starts_with(prefix) {
        format!("{line} environment 1 moisture_content {moisture}")
    } else {line.to_owned()}).collect::<Vec<_>>().join("\n") + "\n"
}

#[test]
fn g3_declared_moisture_reaches_original_reed_mass_stiffness_loss_and_microphone() {
    let pack = environment_pack(pack(900.0, false), true);
    let tensor = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";
    for loss in [tensor, "loss regional-maxwell relaxed 1 boundary_modulus boundary_time"] {
        let text = input().replace(tensor, loss);
        let load = |moisture, block| MaterialModalValvePerformance::from_bytes_with_sources(
            at_environment(&text, "section ", moisture).as_bytes(), block,
            &CancelGate::new_clock_free(), None, |_| Ok(pack.clone())).unwrap();
        let dry = load(0.0, 37);
        let wet = load(0.2, 37);
        let region = |p: &MaterialModalValvePerformance| {
            valve(p).material_plates()[0].source().regions()[0].input.material.clone()
        };
        let a = region(&dry);
        let b = region(&wet);
        assert_eq!(a.card_identity(), b.card_identity(), "one unchanged original source card");
        assert_ne!(a.identity(), b.identity());
        assert_eq!(a.property("density").unwrap().value_si(), 900.0);
        assert_eq!(b.property("density").unwrap().value_si(), 900.0 * 1.1);
        assert_eq!(b.property("young_modulus").unwrap().value_si(), 4e9 * 0.8);
        assert_ne!(valve(&dry).spec().modes, valve(&wet).spec().modes);
        if loss == tensor {
            assert_ne!(valve(&dry).viscosity_matrix_n_s_m(), valve(&wet).viscosity_matrix_n_s_m());
        } else {
            assert!(valve(&wet).material_relaxation().is_some());
        }
        let (reference, _) = render(dry);
        let (changed, accepted) = render(wet);
        assert!(changed.iter().all(|value| value.is_finite()));
        assert!(reference.iter().zip(&changed).any(|(a,b)| (a-b).abs() > 1e-15));
        let frame = accepted.source().last_frame().unwrap();
        let scale = frame.stored_energy_j + frame.storage_change_j.abs()
            + frame.source_work_j.abs() + frame.aperture.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() <= 1e-16 + 1e-7 * scale);
        assert_eq!(changed, render(load(0.2, 1)).0);
    }
}

#[test]
fn g0_native_material_environments_refuse_missing_untyped_and_unsupported_source_states() {
    let text = input();
    let typed = environment_pack(pack(900.0, false), true);
    let untyped = environment_pack(pack(900.0, false), false);
    let check = |text: &str, pack: &NormalizedMaterialCardPack| {
        MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |_| Ok(pack.clone())).is_err()
    };
    assert!(check(&text, &typed), "source response axis must not receive a guessed state");
    assert!(check(&at_environment(&text, "section ", 0.1), &untyped));
    assert!(check(&at_environment(&text, "section ", 0.3), &typed));
    let valid = at_environment(&text, "section ", 0.1);
    for suffix in ["environment 0", "environment 9 moisture_content 0.1",
        "environment 1 moisture_content NaN", "environment 1 T 300",
        "environment 1 omega 2000", "environment 2 moisture_content 0.1 moisture_content 0.2",
        "environment 1 unknown_coordinate 0.1", "environment 1 moisture_content 0.1 extra"] {
        assert!(check(&valid.replace("environment 1 moisture_content 0.1", suffix), &typed), "{suffix}");
    }
}
