use super::*;
use crate::bernoulli_aperture::network::NetworkNode;

const MUTUAL: &str = include_str!("../../../../../../examples/material_modal_valve/paired-mutual.performance");
const OPERATOR: &str = include_str!("../../../../../../examples/material_modal_valve/mutual.radiation");

fn input_walls() -> String {
    include_str!("../../../../../../examples/material_modal_valve/paired-walls.performance")
        .replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")
}

fn upstream(text: &str) -> String {
    text.replace("network 2 1 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_section 1 0 0.12 0.007 0.002",
        "network 3 2 1048576\nduct_node reflection 0.7\nduct_node inlet\nduct_node junction\nduct_section 1 2 0.06 0.007 0.002\nduct_section 2 0 0.06 0.007 0.002 gas 303.15 0")
}
fn with_walls(text: &str, rows: &str) -> String {
    upstream(text).replacen("valve 0.0004\n", &format!("{rows}valve 0.0004\n"), 1)
}
fn wall_pack(young: f64, omit: Option<&str>, max_omega: f64) -> NormalizedMaterialCardPack {
    wall_pack_with_memory(young, omit, max_omega, None)
}
fn wall_pack_with_memory(young: f64, omit: Option<&str>, max_omega: f64,
    memory: Option<(f64, Dims, f64)>) -> NormalizedMaterialCardPack
{
    let mut claims = ClaimSet::new();
    let mut identity = young.to_le_bytes().to_vec();
    if let Some((time, _, upper_omega)) = memory {
        identity.extend_from_slice(&time.to_le_bytes());
        identity.extend_from_slice(&upper_omega.to_le_bytes());
    }
    let artifact = hash_domain("org.frankensim.test.modal-duct-wall", &identity);
    let provenance = Provenance { source: "manufactured local membrane, not measured airway tissue".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic duct wall".into(), method: "manufactured regression".into(), artifact,
        caveats: "not a measured material or instrument".into(), provenance: provenance.clone(),
    }).unwrap();
    let mut properties: Vec<_> = [("density", 1000.0, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS), ("poisson_ratio", 0.3, Dims::NONE),
        ("kelvin_voigt_young_viscosity", 2.0, DynViscosity::DIMS)]
        .into_iter().map(|(name, value, dims)| (name, value, dims, max_omega)).collect();
    if let Some((time, dims, upper_omega)) = memory {
        properties.extend([("wall_modulus_1", 8e5, Pressure::DIMS, max_omega),
            ("wall_modulus_2", 4e5, Pressure::DIMS, max_omega),
            ("wall_time", time, dims, upper_omega)]);
    }
    for (name, value, dims, upper_omega) in properties {
        if omit == Some(name) { continue; }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, upper_omega),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: provenance.clone(), observations: vec![observation.clone()] }).unwrap();
    }
    let normalized = NormalizedPack::new("synthetic-wall", "integration-v1", artifact,
        "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(MaterialStateId { chemistry: "synthetic".into(),
        phase: "solid".into(), process: "manufactured membrane".into(), revision: 0 }, normalized).unwrap()
}
fn loaded(text: &str, wall: &NormalizedMaterialCardPack, block: usize)
    -> Result<MaterialModalValvePerformance, PlateValveInputError>
{
    MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None,
        |path| Ok(if path == "wall.matpack" { wall.clone() } else { pack(900.0, false) }),
        |_| Ok(OPERATOR.as_bytes().to_vec()))
}

fn directional_wall_pack() -> NormalizedMaterialCardPack {
    let base = wall_pack(1e6, None, 100000.0);
    let mut claims = base.card().claims().clone();
    let template = claims.claims_for("young_modulus")[0].1.clone();
    for (name, dims, value) in [("young_modulus_1", Pressure::DIMS, 1e6),
        ("young_modulus_2", Pressure::DIMS, 4e6), ("poisson_ratio_12", Dims::NONE, 0.1),
        ("wall_viscosity_time", Dims([0,0,1,0,0,0]), 2e-6)] {
        let mut claim = template.clone();
        claim.key = PropertyKey::new(name, dims);
        claim.value = PropertyValue::Scalar { value, dims };
        claims.insert_claim(claim).unwrap();
    }
    let normalized = NormalizedPack::new("synthetic-directional-wall", "integration-v1",
        base.source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(base.card().id().clone(), normalized).unwrap()
}

fn normal_viscous_wall_pack(cross: f64) -> NormalizedMaterialCardPack {
    let base = directional_wall_pack();
    let mut claims = base.card().claims().clone();
    let template = claims.claims_for("young_modulus_1")[0].1.clone();
    for (name, value) in [("wall_eta_hh", 400.0), ("wall_eta_ha", cross), ("wall_eta_aa", 90.0)] {
        let mut claim = template.clone();
        claim.key = PropertyKey::new(name, DynViscosity::DIMS);
        claim.value = PropertyValue::Scalar { value, dims:DynViscosity::DIMS };
        claims.insert_claim(claim).unwrap();
    }
    let normalized = NormalizedPack::new("synthetic-normal-viscous-wall", "integration-v1",
        base.source_artifact(), "CC0-1.0", claims, vec![], vec![]).unwrap();
    NormalizedMaterialCardPack::new(base.card().id().clone(), normalized).unwrap()
}

#[test]
fn g3_normal_viscosity_keeps_axial_memory_in_both_original_reed_ducts() {
    use crate::bernoulli_aperture::wall::material::{WallMaterialLoss, WallMaterialModel};
    let text = input_walls().replace("isotropic wall.matpack", "orthotropic wall.matpack")
        .replace("strain-fixed", "stress-free")
        .replace("kelvin-voigt", "1 2 normal-kelvin-voigt 1 0.8 wall_eta_hh wall_eta_ha wall_eta_aa");
    let mut pressures = Vec::new();
    let mut patches = Vec::new();
    let mut reeds = Vec::new();
    for cross in [-60.0, 60.0] {
        let material = normal_viscous_wall_pack(cross);
        let p = loaded(&text, &material, 37).unwrap();
        reeds.push(valve(&p).material_plates()[0].source().specimen_identity());
        let system = p.renderer().source().system();
        let mut original_patches = Vec::new();
        for side in [ValveSide::Upstream, ValveSide::Downstream] {
            let wall = &system.material_walls(side)[0];
            assert!(matches!(wall.source().model, WallMaterialModel::OrthotropicNormal { axes:[0,1], .. }));
            let WallMaterialLoss::Maxwell(spectrum) = &wall.source().loss else { panic!("lost axial memory") };
            assert_eq!(spectrum.material().terms.len(), 1);
            assert_eq!(wall.load().compliance_terms().len(), 1);
            assert!(wall.patch().wall.resistance > 0.0);
            original_patches.push(wall.patch());
        }
        patches.push(original_patches);
        let mut renderer = p.into_renderer();
        let mut pressure = Vec::new();
        let mut maximum_memory = [0.0_f64; 2];
        while renderer.remaining_samples() != 0 {
            let mut sample = [0.0]; renderer.block(&mut sample).unwrap();
            pressure.push(sample[0]);
            let system = renderer.source().system();
            let frame = renderer.source().last_frame().unwrap();
            let scale = frame.stored_energy_j + frame.storage_change_j.abs()
                + frame.source_work_j.abs() + frame.aperture.mechanical_work_j.abs();
            assert!(frame.balance_residual_j().abs() <= 1e-16 + 1e-7 * scale);
            for (i, side, node) in [(0, ValveSide::Upstream, 2), (1, ValveSide::Downstream, 1)] {
                let energy = system.wall_relaxation_energy_j(side, node).unwrap();
                assert!(energy.is_finite() && energy >= 0.0);
                maximum_memory[i] = maximum_memory[i].max(energy);
            }
        }
        assert!(maximum_memory.iter().all(|energy| *energy > 0.0));
        assert!(pressure.iter().all(|pressure| pressure.is_finite()));
        assert_eq!(pressure, render(loaded(&text, &material, 37).unwrap()).0);
        pressures.push(pressure);
    }
    assert_eq!(reeds[0], reeds[1]);
    assert_eq!(patches[0], patches[1]);
    assert!(pressures[0].iter().zip(&pressures[1]).any(|(a,b)| (a-b).abs() > 1e-15),
        "a data-only change of axial memory must reach the actual exterior receiver");
}

#[test]
fn g3_directional_duct_walls_share_the_original_reed_and_external_pressure_owner() {
    use crate::bernoulli_aperture::wall::material::WallMaterialModel;
    let text = input_walls().replace("isotropic wall.matpack", "orthotropic wall.matpack")
        .replace("stress-free kelvin-voigt", "stress-free 1 2 proportional-kelvin-voigt wall_viscosity_time")
        .replace("strain-fixed kelvin-voigt", "strain-fixed 1 2 proportional-kelvin-voigt wall_viscosity_time");
    let material = directional_wall_pack();
    let p = loaded(&text, &material, 37).unwrap();
    let reed = valve(&p).material_plates()[0].source().specimen_identity();
    for side in [ValveSide::Upstream, ValveSide::Downstream] {
        let wall = &p.renderer().source().system().material_walls(side)[0];
        assert_eq!(wall.source().model, WallMaterialModel::Orthotropic { axes:[0,1],
            proportional_viscosity_time_property:Some("wall_viscosity_time".into()) });
        assert!(wall.source().material.property("young_modulus").is_none());
        assert_eq!(wall.source().material.card_identity(), material.card().content_hash());
        assert_eq!(wall.patch().wall.surface_density, 0.2);
    }
    let reference = render(p).0;
    assert!(reference.iter().any(|pressure| pressure.abs()>1e-15));
    assert_eq!(reference, render(loaded(&text, &material, 1).unwrap()).0);
    let rotated = text.replace("strain-fixed 1 2", "strain-fixed 2 1");
    let p = loaded(&rotated, &material, 37).unwrap();
    assert_eq!(valve(&p).material_plates()[0].source().specimen_identity(), reed);
    let wall = &p.renderer().source().system().material_walls(ValveSide::Downstream)[0];
    let expected = 4e6/0.96 * 0.0002 / 0.007f64.powi(2);
    assert!((wall.patch().wall.stiffness_per_area/expected-1.0).abs()<1e-14);
    assert!(reference.iter().zip(render(p).0).any(|(a,b)| (a-b).abs()>1e-15));
    // Even a complete isotropic bundle cannot replace the missing directional claims.
    assert!(loaded(&text, &wall_pack(1e6, None, 100000.0), 37).is_err());
    for invalid in [text.replace("1 2 proportional", "1 1 proportional"),
        text.replace("proportional-kelvin-voigt wall_viscosity_time", "kelvin-voigt")] {
        assert!(loaded(&invalid, &material, 37).is_err());
    }
}

#[test]
fn g0_material_walls_share_original_sources_geometry_gases_and_both_radiation_paths() {
    let material = wall_pack(1e6, None, 100000.0);
    let text = input_walls();
    let mut loads = [0; 2];
    let p = MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            if path == "wall.matpack" { loads[0] += 1; Ok(material.clone()) }
            else { loads[1] += 1; Ok(pack(900.0, false)) }
        }).unwrap();
    assert_eq!(loads, [1, 1]);
    assert_eq!(p.wall_records().len(), 2);
    assert_eq!(p.wall_records()[0].pack(), &material);
    let system = p.renderer().source().system();
    for (side, node, length, axial) in [(ValveSide::Upstream, 2, 0.06, 1.0),
        (ValveSide::Downstream, 1, 0.125, 1.0 / (1.0 - 0.3f64.powi(2)))] {
        let wall = &system.material_walls(side)[0];
        assert_eq!(wall.source().node, node);
        assert_eq!(wall.length_m(), length);
        assert_eq!(wall.patch().wall.surface_density, 0.2);
        assert!((wall.patch().wall.stiffness_per_area / (1e6 * 0.0002 / 0.007f64.powi(2) * axial) - 1.0).abs() < 1e-14);
        assert_eq!(system.source(side).graph.nodes[node], NetworkNode::Junction);
        assert_eq!(system.graph(side).nodes[node], wall.shunt());
        assert_eq!(system.source(side).gases[1].temperature, 303.15);
        assert_eq!(wall.source().material.card_identity(), material.card().content_hash());
    }
    assert_eq!(p.radiation_bindings().len(), 1);
    assert_eq!(valve(&p).material_plates().len(), 2);
    let up_only = "walls 1 20 1000\nwall upstream 2 isotropic wall.matpack 293.15 100 0.0002 stress-free kelvin-voigt\n";
    for original in [input(), MUTUAL.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 512 0.005")] {
        for supply in ["supply none", "supply slit 0.00002 0.01 0.001 16"] {
            let text = with_walls(&original.replace("supply none", supply), up_only);
            let p = loaded(&text, &material, 37).unwrap();
            let system = p.renderer().source().system();
            assert_eq!(system.material_walls(ValveSide::Upstream).len(), 1);
            assert_eq!(system.supply_channel().is_some(), supply != "supply none");
            let (pressure, _) = render(p);
            assert!(pressure.iter().any(|p| p.abs() > 1e-15));
        }
    }
}

#[test]
fn g3_tract_or_bore_material_changes_accepted_motion_and_pressure_without_changing_the_reed() {
    let text = input_walls();
    let material = wall_pack(1e6, None, 100000.0);
    let p = loaded(&text, &material, 37).unwrap();
    let reed = valve(&p).material_plates()[0].source().specimen_identity();
    let mut renderer = p.into_renderer();
    let mut pressure = Vec::new();
    let mut wall_energy = [0.0_f64; 2];
    while renderer.remaining_samples() != 0 {
        let mut sample = [0.0]; renderer.block(&mut sample).unwrap(); pressure.push(sample[0]);
        let system = renderer.source().system();
        let frame = renderer.source().last_frame().unwrap();
        let scale = frame.stored_energy_j + frame.storage_change_j.abs()
            + frame.source_work_j.abs() + frame.aperture.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() <= 1e-16 + 1e-7 * scale);
        assert_eq!(frame.stored_energy_j, system.stored_energy_j());
        for (i, side, node) in [(0, ValveSide::Upstream, 2), (1, ValveSide::Downstream, 1)] {
            wall_energy[i] = wall_energy[i].max(system.wall_state(side, node).unwrap().stored_energy_j);
        }
    }
    assert!(wall_energy.iter().all(|e| *e > 0.0));
    assert_eq!(pressure, render(loaded(&text, &material, 37).unwrap()).0);
    for side in ["upstream 2", "downstream 1"] {
        let changed = text.replace(&format!("wall {side} isotropic wall.matpack"),
            &format!("wall {side} isotropic altered-wall.matpack"));
        let p = MaterialModalValvePerformance::from_bytes_with_sources(changed.as_bytes(), 37,
            &CancelGate::new_clock_free(), None, |path| Ok(match path {
                "wall.matpack" => material.clone(),
                "altered-wall.matpack" => wall_pack(4e6, None, 100000.0),
                _ => pack(900.0, false),
            })).unwrap();
        assert_eq!(valve(&p).material_plates()[0].source().specimen_identity(), reed);
        let (changed_pressure, changed) = render(p);
        assert!(pressure.iter().zip(changed_pressure).any(|(a,b)| (a-b).abs() > 1e-15));
        assert_ne!(renderer.source().last_frame().unwrap().aperture,
            changed.source().last_frame().unwrap().aperture);
    }
    let legacy = upstream(&input());
    assert_eq!(render(build(&legacy, 37)).0, render(build(&with_walls(&input(), "walls none\n"), 1)).0);
}

#[test]
fn g0_missing_out_of_band_and_incompatible_wall_sources_refuse() {
    let text = input_walls();
    for source in [wall_pack(1e6, Some("kelvin_voigt_young_viscosity"), 100000.0),
        wall_pack(1e6, None, 1000.0)] {
        assert!(loaded(&text, &source, 37).is_err());
    }
    let material = wall_pack(1e6, None, 100000.0);
    for changed in [
        text.replace("walls 2 20 1000", "walls 129 20 1000"),
        text.replace("walls 2 20 1000", "walls 2 20 999"),
        text.replace("wall downstream 1", "wall upstream 2"),
        text.replace("wall downstream 1", "wall downstream 2"),
        text.replace("wall.matpack", "../wall.matpack"),
        text.replace("wall.matpack 303.15", "wall.matpack 500"),
        text.replace("0.0002 strain-fixed", "0.001 strain-fixed"),
    ] { assert!(loaded(&changed, &material, 37).is_err()); }
}

const WALL_MAXWELL: &str =
    "maxwell 2 1 0.1 wall_modulus_1 wall_time wall_modulus_2 wall_time";
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

#[test]
fn g3_source_wall_relaxation_reaches_both_ducts_reed_motion_and_microphone() {
    let text = input_walls().replace("kelvin-voigt", WALL_MAXWELL)
        .replace("wall upstream 2 isotropic wall.matpack", "wall upstream 2 isotropic tract-wall.matpack")
        .replace("wall downstream 1 isotropic wall.matpack", "wall downstream 1 isotropic bore-wall.matpack");
    let original = wall_pack_with_memory(1e6, None, 100000.0, Some((0.001, TIME, 100000.0)));
    let substituted = wall_pack_with_memory(1e6, None, 100000.0, Some((0.0001, TIME, 100000.0)));
    // The performance bytes, gas states, equilibrium stiffness, geometry and
    // reed sources stay fixed. Only one wall card's relaxation time changes.
    let make = |changed: Option<&str>, block| {
        MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), block,
            &CancelGate::new_clock_free(), None, |path| Ok(match path {
                "tract-wall.matpack" | "bore-wall.matpack" => {
                    if Some(path) == changed { substituted.clone() } else { original.clone() }
                }
                _ => pack(900.0, false),
            })).unwrap()
    };
    let p = make(None, 37);
    let reed = valve(&p).material_plates()[0].source().specimen_identity();
    for record in p.wall_records() { assert_eq!(record.pack(), &original); }
    let system = p.renderer().source().system();
    for side in [ValveSide::Upstream, ValveSide::Downstream] {
        let wall = &system.material_walls(side)[0];
        let crate::bernoulli_aperture::wall::material::WallMaterialLoss::Maxwell(spectrum)
            = &wall.source().loss else { panic!("source spectrum was not retained"); };
        assert_eq!(spectrum.material().e_inf, 1e6);
        assert_eq!(spectrum.material().terms, vec![(8e5, 0.001), (4e5, 0.001)]);
        assert_eq!(wall.load().compliance_terms().len(), 2);
        assert_eq!(wall.patch().wall.resistance, 0.0, "Maxwell must not also add Kelvin-Voigt loss");
    }
    let mut renderer = p.into_renderer();
    let mut pressure = Vec::new();
    let mut maximum_memory_energy = [0.0_f64; 2];
    while renderer.remaining_samples() != 0 {
        let mut sample = [0.0];
        renderer.block(&mut sample).unwrap();
        pressure.push(sample[0]);
        let system = renderer.source().system();
        let frame = renderer.source().last_frame().unwrap();
        let scale = frame.stored_energy_j + frame.storage_change_j.abs()
            + frame.source_work_j.abs() + frame.aperture.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() <= 1e-16 + 1e-7 * scale);
        for (i, side, node) in [(0, ValveSide::Upstream, 2), (1, ValveSide::Downstream, 1)] {
            let memory = system.wall_relaxation_energy_j(side, node).unwrap();
            assert!(memory.is_finite() && memory >= 0.0);
            maximum_memory_energy[i] = maximum_memory_energy[i].max(memory);
        }
    }
    assert!(maximum_memory_energy.iter().all(|e| *e > 0.0));
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    let (partitioned, chunked) = render(make(None, 37));
    assert_eq!(pressure, partitioned);
    assert_eq!(renderer.source().last_frame().unwrap().aperture,
        chunked.source().last_frame().unwrap().aperture);
    for (side, node) in [(ValveSide::Upstream, 2), (ValveSide::Downstream, 1)] {
        assert_eq!(renderer.source().system().wall_relaxation_energy_j(side, node).unwrap(),
            chunked.source().system().wall_relaxation_energy_j(side, node).unwrap());
    }
    for source in ["tract-wall.matpack", "bore-wall.matpack"] {
        let p = make(Some(source), 37);
        assert_eq!(valve(&p).material_plates()[0].source().specimen_identity(), reed);
        let (changed_pressure, changed) = render(p);
        assert!(pressure.iter().zip(changed_pressure).any(|(a, b)| (a - b).abs() > 1e-15));
        assert_ne!(renderer.source().last_frame().unwrap().aperture,
            changed.source().last_frame().unwrap().aperture);
    }
}

#[test]
fn g0_wall_relaxation_requires_source_dimensions_band_and_resolved_times() {
    let text = input_walls().replace("kelvin-voigt", WALL_MAXWELL);
    for material in [
        wall_pack_with_memory(1e6, Some("wall_modulus_1"), 100000.0, Some((0.001, TIME, 100000.0))),
        wall_pack_with_memory(1e6, Some("wall_time"), 100000.0, Some((0.001, TIME, 100000.0))),
        wall_pack_with_memory(1e6, None, 100000.0, Some((0.001, Pressure::DIMS, 100000.0))),
        wall_pack_with_memory(1e6, None, 100000.0, Some((0.001, TIME, 1000.0))),
        wall_pack_with_memory(1e6, None, 100000.0, Some((0.0, TIME, 100000.0))),
        wall_pack_with_memory(1e6, None, 100000.0, Some((1e-7, TIME, 100000.0))),
    ] {
        assert!(loaded(&text, &material, 37).is_err(), "invalid hereditary source was accepted");
    }
    let material = wall_pack_with_memory(1e6, None, 100000.0, Some((0.001, TIME, 100000.0)));
    for changed in [
        text.replace("maxwell 2", "maxwell 0"),
        text.replace("maxwell 2", "maxwell 9"),
        text.replace("wall_modulus_2 wall_time", "wall_modulus_1 wall_time"),
        text.replace("maxwell 2 1 0.1", "maxwell 2 1 0.0001"),
    ] {
        assert!(loaded(&changed, &material, 37).is_err(), "invalid spectrum declaration was accepted");
    }
}

#[test]
fn g3_held_wall_environment_reaches_both_ducts_and_the_same_coupled_reed() {
    let material = environment_pack(wall_pack(1e6, None, 100000.0), true);
    let dry_text = at_environment(&input_walls(), "wall ", 0.0);
    let wet_text = at_environment(&input_walls(), "wall ", 0.2);
    let dry = loaded(&dry_text, &material, 37).unwrap();
    let wet = loaded(&wet_text, &material, 37).unwrap();
    assert_eq!(valve(&dry).material_plates()[0].source().specimen_identity(),
        valve(&wet).material_plates()[0].source().specimen_identity());
    for side in [ValveSide::Upstream, ValveSide::Downstream] {
        let a = &dry.renderer().source().system().material_walls(side)[0];
        let b = &wet.renderer().source().system().material_walls(side)[0];
        assert_eq!(a.source().material.card_identity(), b.source().material.card_identity());
        assert_ne!(a.source().material.identity(), b.source().material.identity());
        assert_ne!(a.patch(), b.patch());
    }
    let (reference, dry) = render(dry);
    let (changed, wet) = render(wet);
    assert!(changed.iter().all(|p| p.is_finite()));
    assert!(reference.iter().zip(&changed).any(|(a,b)| (a-b).abs() > 1e-15));
    let motion = |source: &ScheduledTwoSidedPressure| match source.source().system().valve() {
        LoadedAperture::Modal(valve) => valve.modal_states().to_vec(),
        _ => panic!("original reed motion"),
    };
    assert_ne!(motion(&dry), motion(&wet), "wall state changes actual acoustic back-reaction");
    assert_eq!(changed, render(loaded(&wet_text, &material, 1).unwrap()).0);

    // Uniform and regional DG adapters call this same original wall input.
    // Geometry, constitutive selection and use band remain fixed.
    let load_dg = |moisture| {
        let text = format!("isotropic wall.matpack 293.15 1000 0.001 strain-fixed 50 2000 kelvin-voigt environment 1 moisture_content {moisture}");
        crate::dg_port::input::wall::DgMaterialWallInput::read(1, &mut text.split_ascii_whitespace())
            .unwrap().unwrap().resolve(std::sync::Arc::new(material.clone()),
                &CancelGate::new_clock_free()).unwrap()
    };
    let a = load_dg(0.0);
    let b = load_dg(0.2);
    assert_eq!(a.pack().content_hash(), b.pack().content_hash());
    assert_ne!(a.wall().material().identity(), b.wall().material().identity());
    assert_eq!(a.wall().thickness_m(), b.wall().thickness_m());
    assert_ne!(a.wall().specific_impedance(0.007, 1000.0).unwrap(),
        b.wall().specific_impedance(0.007, 1000.0).unwrap());
}
