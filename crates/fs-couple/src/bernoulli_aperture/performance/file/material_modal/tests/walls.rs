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
    let mut claims = ClaimSet::new();
    let artifact = hash_domain("org.frankensim.test.modal-duct-wall", &young.to_le_bytes());
    let provenance = Provenance { source: "manufactured local membrane, not measured airway tissue".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic duct wall".into(), method: "manufactured regression".into(), artifact,
        caveats: "not a measured material or instrument".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, value, dims) in [("density", 1000.0, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS), ("poisson_ratio", 0.3, Dims::NONE),
        ("kelvin_voigt_young_viscosity", 2.0, DynViscosity::DIMS)] {
        if omit == Some(name) { continue; }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, max_omega),
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
