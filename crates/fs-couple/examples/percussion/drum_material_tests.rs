//! Source properties enter the original head pencil and accepted drum motion.
use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack,
    ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_qty::{Density, Dims, Pressure, Temperature};

fn input() -> String {
    "frankensim-drum-spec-v1\ngeometry,0.1703,0.1651,0.1778\n\
     head_material,batter,0.000254,3000,stated-unknown,film.pack synthetic T=280\n\
     head_material,resonant,0.0000762,1500,stated-unknown,film.pack synthetic T=280\n\
     mesh,2,8\nband_hz,80,500\nambient,293.15 101325 0.65\n".into()
}

fn pack(omit: Option<&str>, frequency_bound: bool) -> NormalizedMaterialCardPack {
    let artifact = ContentHash([27;32]);
    let provenance = Provenance { source: "synthetic head material; not measured film".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact) };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: "synthetic isotropic film".into(), method: "independent scaling fixture".into(),
        artifact, caveats: "not calibrated material data".into(), provenance: provenance.clone(),
    }).unwrap();
    for (name, dims, first) in [("density", Density::DIMS, 1390.0),
        ("young_modulus", Pressure::DIMS, 4e9), ("poisson_ratio", Dims::NONE, 0.38)] {
        if omit == Some(name) { continue; }
        let constant = name == "poisson_ratio";
        let value = if constant { PropertyValue::Scalar { value: first, dims } }
            else { PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, first), (320.0, 2.0*first)], dims } };
        let mut validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        if frequency_bound { validity = validity.with("omega", 0.0, 100_000.0); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims), value, validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if constant { InterpolationPolicy::ConstantWithinValidity }
                else { InterpolationPolicy::LinearInside },
            observations: vec![observation.clone()], provenance: provenance.clone() }).unwrap();
    }
    NormalizedMaterialCardPack::new(MaterialStateId { chemistry: "synthetic-film".into(),
        phase: "solid".into(), process: "synthetic".into(), revision: 0 },
        NormalizedPack::new("synthetic-film", "synthetic-head-v1", artifact, "CC0-1.0",
            claims, vec![], vec![]).unwrap()).unwrap()
}

fn read(text: &str, pack: &NormalizedMaterialCardPack) -> Spec {
    Spec::read_with_all_sources(text, None, |_| Ok(pack.clone())).unwrap()
}

#[test]
fn source_coefficients_preserve_the_numeric_head_pencil_and_keep_actual_sources() {
    let pack = pack(None, false); let mut reads = 0;
    let sourced = Spec::read_with_all_sources(&input(), None, |name| {
        reads += 1; assert_eq!(name, "film.pack"); Ok(pack.clone())
    }).unwrap();
    assert_eq!(reads, 1, "the same original pack supplies both heads");
    let mut direct = sourced.clone(); direct.head_materials = [None, None];
    for i in 0..2 {
        let a = sourced.head(i).unwrap(); let b = direct.head(i).unwrap();
        assert_eq!(a.model.m, b.model.m); assert_eq!(a.model.k, b.model.k);
        assert_eq!(a.mesh, b.mesh); assert_eq!(a.mass_kg, b.mass_kg);
        let source = sourced.head_materials[i].as_ref().unwrap();
        assert_eq!(source.pack().content_hash(), pack.content_hash());
        assert_eq!(source.state().card_identity(), pack.card().content_hash());
        assert_eq!(source.record().axes, [("T".to_string(), 280.0)]);
        assert!(source.receipt().contains(&pack.content_hash().to_hex()));
    }
    // The entire execution remains bitwise equal to explicitly entering the
    // admitted coefficients, including the same external gas and installed load.
    let build = |s| crate::drum_with_spec(128, 2e-6, false, false, None, true,
        crate::Stroke { speed_m_s: 2.0, position_m: Some([0.06, 0.01]) }, false, None, Some(s)).unwrap();
    let mut a = build(sourced); let mut b = build(direct);
    assert_eq!(a.observer_a, b.observer_a); assert_eq!(a.force, b.force);
    assert_eq!(a.head_materials[0].as_ref().unwrap().pack().content_hash(), pack.content_hash());
    let gate = CancelGate::new_clock_free();
    for _ in 0..128 {
        let x = a.system.step(&a.force, &gate).unwrap(); let y = b.system.step(&b.force, &gate).unwrap();
        assert_eq!(x.stored_energy_j.to_bits(), y.stored_energy_j.to_bits());
        assert_eq!(a.system.state(), b.system.state());
    }
}

#[test]
fn same_card_temperature_changes_mass_bending_stretching_and_physical_pressure() {
    let pack = pack(None, false);
    let cold = read(&input(), &pack); let warm = read(&input().replace("T=280", "T=320"), &pack);
    for i in 0..2 {
        let a = cold.head(i).unwrap(); let b = warm.head(i).unwrap();
        assert_eq!(a.mesh, b.mesh);
        assert_eq!(a.spec.thickness_m, b.spec.thickness_m);
        assert_eq!(a.spec.tension_n_m, b.spec.tension_n_m);
        assert_eq!(b.mass_kg, 2.0*a.mass_kg);
        // Independently assemble zero-tension bending. The original pretension
        // stays fixed: doubling E adds K_b, while doubling density adds M.
        let rim: Vec<_> = (0..a.mesh.nodes.len()).filter(|n| a.model.dof_map[3*n].is_none()).collect();
        let bending = fs_plate::assemble(&a.mesh, &a.section, &rim, &[], &fs_plate::AssemblyOptions {
            pretension: 0.0, support: fs_plate::EdgeSupport::SimplySupported }).unwrap();
        for row in 0..a.model.free { for col in 0..a.model.free {
            assert_eq!(b.model.m.get(row, col), 2.0*a.model.m.get(row, col));
            let expected = a.model.k.get(row, col) + bending.k.get(row, col);
            assert!((b.model.k.get(row, col)-expected).abs() <= 1e-10*(1.0+expected.abs()));
        }}
    }
    let build = |s| crate::drum_with_spec(128, 2e-6, false, false, None, true,
        crate::Stroke { speed_m_s: 2.0, position_m: Some([0.06, 0.01]) }, false, None, Some(s)).unwrap();
    let mut first = build(cold.clone()); let mut changed = build(warm); let mut replay = build(cold);
    assert_eq!(first.ambient, changed.ambient);
    assert_eq!(first.pressure.as_ref().unwrap().volume_m3, changed.pressure.as_ref().unwrap().volume_m3);
    let gate = CancelGate::new_clock_free(); let mut changed_pressure = 0.0_f64; let mut stretching = 0.0_f64;
    for _ in 0..128 {
        for e in [&mut first, &mut changed, &mut replay] {
            let frame = e.system.step(&e.force, &gate).unwrap();
            assert!(frame.balance_residual_j.abs() < 2e-7);
            stretching = stretching.max(e.system.membrane_observation(1).unwrap().stretching_energy_j);
        }
        assert_eq!(first.system.state(), replay.system.state());
        let p = crate::cavity_pressure(first.pressure.as_ref().unwrap(), first.system.state());
        let q = crate::cavity_pressure(changed.pressure.as_ref().unwrap(), changed.system.state());
        changed_pressure = changed_pressure.max((p-q).abs());
    }
    assert!(stretching > 0.0); assert!(changed_pressure > 1e-9);
}

#[test]
fn original_measured_pack_loads_relative_to_the_drum_and_retains_its_identity() {
    let claims = NormalizedPack::from_bytes(include_bytes!(
        "../../tests/fixtures/matdb/aluminum-2024-t3-nasa-tn-d6448.fsmatpk")).unwrap();
    let pack = NormalizedMaterialCardPack::new(MaterialStateId { chemistry: "aluminum-2024-t3".into(),
        phase: "solid".into(), process: "as-tabulated".into(), revision: 0 }, claims).unwrap();
    let dir = std::env::temp_dir().join(format!("fs-drum-material-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("film.pack"), pack.to_bytes()).unwrap();
    let text = input().replace("synthetic T=280", "as-tabulated T=300");
    std::fs::write(dir.join("source.fsd"), &text).unwrap();
    let file = Spec::load(dir.join("source.fsd").to_str().unwrap()).unwrap();
    let direct = read(&text, &pack);
    assert_eq!(file, direct);
    for source in file.head_materials.iter().flatten() {
        assert_eq!(source.pack().to_bytes(), pack.to_bytes());
    }
    assert!(Spec::read(&text).is_err(), "data-only reader cannot silently open source files");
}

#[test]
fn malformed_foreign_state_missing_dynamic_and_overwritten_sources_refuse() {
    let good = pack(None, false);
    for bad in [input().replace("synthetic T=280", "foreign T=280"), input().replace("T=280", "T=350"),
        input().replace("film.pack", "../film.pack"), input().replace("T=280", "T=280 T=290"),
        input().replace("stated-unknown", "0.001"),
        format!("{}head,batter,0.000254,4000000000,0.38,1390,3000,stated-unknown\n", input()),
        format!("{}head_material,batter,0.000254,3000,stated-unknown,film.pack synthetic T=280\n", input())] {
        assert!(Spec::read_with_all_sources(&bad, None, |_| Ok(good.clone())).is_err());
    }
    assert!(Spec::read_with_all_sources(&input(), None, |_| Err("missing pack".into())).is_err());
    let missing = pack(Some("poisson_ratio"), false);
    assert!(Spec::read_with_all_sources(&input(), None, |_| Ok(missing.clone())).is_err());
    let dynamic = pack(None, true);
    assert!(Spec::read_with_all_sources(&input().replace("T=280", "T=280 omega=1000"), None,
        |_| Ok(dynamic.clone())).is_err(), "one frequency sample is not a frozen elastic material law");
    let mut altered = read(&input(), &good); altered.heads[0].density_kg_m3 *= 2.0;
    assert!(altered.head(0).is_err(), "source identity cannot authorize overwritten inertia");
    let relaxation = read(&input().replace("stated-unknown", "relaxation"), &good);
    assert!(relaxation.admit_head_loss(false).is_err());
    assert!(relaxation.admit_head_loss(true).is_ok());
    let law = crate::head_relaxation::Spec::read(include_str!("estimated-head-relaxation.fshr")).unwrap();
    law.admit(&relaxation).unwrap();
    let changed = read(&input().replace("stated-unknown", "relaxation").replace("T=280", "T=320"), &good);
    assert!(law.admit(&changed).is_err(), "source stiffness cannot silently retune a distinct Maxwell law");
    let original = include_str!("estimated_drum.fsd"); let mut reads = 0;
    let legacy = Spec::read_with_all_sources(original, None, |_| { reads += 1; Err("unexpected read".into()) }).unwrap();
    assert_eq!(legacy, Spec::read(original).unwrap()); assert_eq!(reads, 0);
}
