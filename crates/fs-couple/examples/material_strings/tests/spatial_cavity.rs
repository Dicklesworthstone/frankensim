//! Authored original volume, shell, atmosphere and material reach the same
//! mechanical performance. This manufactured closed shell is not an instrument.
use super::*;

#[path = "spatial_cavity/openings.rs"]
mod openings;
#[path = "spatial_cavity/thermal.rs"]
mod thermal;
#[path = "spatial_cavity/regional.rs"]
mod regional;
#[path = "spatial_cavity/molecular.rs"]
mod molecular;

const BODY: &str = include_str!("../closed-shell-cavity.body");
const PERFORMANCE: &str = include_str!("../closed-shell-cavity.material-strings");
const RESOLUTION: &str = "interior-fe-v1 2 1100 24 256 80 10000000 1e-8";

fn shell_pack(name: &str, density: f64) -> NormalizedPack {
    let provenance = Provenance {
        source: "synthetic closed-shell cavity integration".into(),
        license: "CC0-1.0".into(), artifact: Some(fs_blake3::ContentHash([7;32])),
    };
    let mut claims = ClaimSet::new();
    let observation = claims.register_observation(ObservationDataset {
        specimen: name.into(), method: "manufactured coefficients".into(),
        artifact: fs_blake3::ContentHash([7;32]), caveats: "not measured instrument data".into(),
        provenance: provenance.clone(),
    }).unwrap();
    let mut entries = vec![
        ("density", Density::DIMS, density),
        ("young_modulus", Pressure::DIMS, 100_000.0),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("shear_modulus_13", Pressure::DIMS, 100_000.0/2.6),
        ("shear_modulus_23", Pressure::DIMS, 100_000.0/2.6),
    ];
    for key in [
        "plane_stress_viscosity_11", "plane_stress_viscosity_22",
        "plane_stress_viscosity_66", "plane_stress_viscosity_12",
        "plane_stress_viscosity_16", "plane_stress_viscosity_26",
        "transverse_shear_viscosity_13", "transverse_shear_viscosity_23",
        "transverse_shear_viscosity_13_23",
    ] { entries.push((key, DynViscosity::DIMS, 0.0)); }
    for (key,dims,value) in entries {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key,dims), value: PropertyValue::Scalar {value,dims},
            validity: ValidityDomain::unconstrained().with("T",280.0,320.0)
                .with("omega",0.0,200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation.clone()], provenance: provenance.clone(),
        }).unwrap();
    }
    let pack = NormalizedPack::new(name,"synthetic-closed-shell-v1",
        fs_blake3::ContentHash([7;32]),"synthetic redistribution permitted",
        claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

fn source(text: &str) -> input::Spec {
    input::Spec::parse(text.replace("48000 1024", "48000 256").as_bytes()).unwrap()
}
fn ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, body: &str, density: f64) -> Ready {
    spec.compile_with_body(cx,gate,Some(body::Profile::parse(body.as_bytes()).unwrap()),
        |name| Ok(if name == "closed-shell.fsmatpk" { shell_pack(name,density) }
            else { pack(name,7800.0,200_000.0,None) })).unwrap()
}
fn uniform_stiffness(ready: &Ready) -> f64 {
    let source = ready.performance.bank().source().bank();
    let count = source.spatial_cavity().unwrap().coupling().cavity_modes();
    let links = source.body_connections();
    let spring = &links[links.len()-count];
    spring.stiffness_n_m * spring.left.shapes[0].powi(2)
}

#[test]
fn authored_volume_reaches_shared_body_and_string_motion_before_preload_and_pcm() {
    scoped(|cx,gate| {
        let mut loaded = ready(cx,gate,source(PERFORMANCE),BODY,700.0);
        let mut replay = ready(cx,gate,source(PERFORMANCE),BODY,700.0);
        let cavity = loaded.performance.bank().source().bank().spatial_cavity().unwrap();
        assert_eq!(*cavity.gas(),loaded.spec.gas);
        assert_eq!(cavity.order(),2);
        assert_eq!(cavity.positions().len(),43);
        assert_eq!(cavity.tetrahedra().len(),80);
        assert_eq!(cavity.pressure_shapes()[0].len(),205);
        assert!(cavity.coupling().cavity_modes() > 8);
        let retained = cavity.coupling().cavity_modes();
        assert_eq!(loaded.performance.bank().system().mode_count(),4+1+retained-1);
        assert_eq!(cavity.coupling().total_modes(),cavity.coupling().structural_modes()+retained-1);
        let original = loaded.performance.bank().source().bank().material_shell_image().unwrap();
        let expected = original.cavity_volume_shapes_m2_per_sqrt_kg()[0]/cavity.volume_m3().sqrt();
        assert!((cavity.overlap_rows()[0]-expected).abs() < 1e-10*expected.abs());
        assert!(loaded.performance.preload_energy_j().unwrap() > 0.0);
        let audio = pressure(&mut loaded,1);
        assert!(audio.iter().any(|p| p.abs() > 1e-14));
        assert_eq!(audio,pressure(&mut replay,37));
        assert!(loaded.performance.bank().system().components()[2].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        let without_air = BODY.split_once("interior-fe-v1").unwrap().0;
        let mut bare = ready(cx,gate,source(PERFORMANCE),without_air,700.0);
        let bare_audio = pressure(&mut bare,128);
        assert!(audio.iter().zip(bare_audio).any(|(p,q)| (p-q).abs() > 1e-14));
        assert_ne!(loaded.performance.bank().system().components()[0].states(),
            bare.performance.bank().system().components()[0].states());
        let frame = loaded.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn source_pressure_composition_and_solid_substitutions_rebuild_physical_storage() {
    scoped(|cx,gate| {
        let ordinary = ready(cx,gate,source(PERFORMANCE),BODY,700.0);
        let doubled = ready(cx,gate,source(&PERFORMANCE.replace("293.15 101325", "293.15 202650")),BODY,700.0);
        let dense = ready(cx,gate,source(PERFORMANCE),BODY,1400.0);
        let mixed = ready(cx,gate,source(&PERFORMANCE.replace("mole 2 N2 0.8 O2 0.2",
            "mole 3 N2 0.7 O2 0.2 CO2 0.1")),BODY,700.0);
        let a = ordinary.performance.bank().source().bank().spatial_cavity().unwrap();
        let b = doubled.performance.bank().source().bank().spatial_cavity().unwrap();
        let d = dense.performance.bank().source().bank().spatial_cavity().unwrap();
        let m = mixed.performance.bank().source().bank().spatial_cavity().unwrap();
        assert_eq!(a.basis_identity(),b.basis_identity());
        assert_eq!(a.basis_identity(),m.basis_identity());
        assert_ne!(a.basis_identity(),d.basis_identity());
        assert_eq!(a.positions(),d.positions());
        assert_eq!(a.tetrahedra(),m.tetrahedra());
        assert!((uniform_stiffness(&doubled)/uniform_stiffness(&ordinary)-2.0).abs() < 1e-9);
        assert!((uniform_stiffness(&dense)/uniform_stiffness(&ordinary)-0.5).abs() < 1e-8);
        assert_eq!(*m.gas(),mixed.spec.gas);
        let ratio = m.modal_report().frequency_hz[1]/a.modal_report().frequency_hz[1];
        assert!((ratio-mixed.spec.gas.sound_speed/ordinary.spec.gas.sound_speed).abs() < 2e-7);
        let bulk_ratio = mixed.spec.gas.density*mixed.spec.gas.sound_speed.powi(2)
            /(ordinary.spec.gas.density*ordinary.spec.gas.sound_speed.powi(2));
        assert!((uniform_stiffness(&mixed)/uniform_stiffness(&ordinary)-bulk_ratio).abs() < 1e-8);
    });
}

#[test]
fn uniform_compression_and_distributed_air_keep_the_complete_resource_and_fluid_counts() {
    let full = body::Profile::parse(BODY.as_bytes()).unwrap();
    let spec = source(PERFORMANCE);
    assert_eq!(full.additional_body_connections(),24);
    assert_eq!(full.additional_body_modes(),23);
    assert_eq!(spec.network_dimensions_with_body_modes(true,24,23).unwrap(),(28,26,28*27*27));
    let uniform = BODY.replace(RESOLUTION,"interior-fe-v1 2 1 1 256 80 10000000 1e-8");
    let profile = body::Profile::parse(uniform.as_bytes()).unwrap();
    assert_eq!(profile.additional_body_connections(),1);
    assert_eq!(profile.additional_body_modes(),0);
    assert_eq!(spec.network_dimensions_with_body_modes(true,1,0).unwrap(),(5,3,5*4*4));
    scoped(|cx,gate| {
        let compressed = ready(cx,gate,spec,&uniform,700.0);
        let cavity = compressed.performance.bank().source().bank().spatial_cavity().unwrap();
        assert_eq!(cavity.coupling().cavity_modes(),1);
        assert_eq!(cavity.coupling().total_modes(),1);
        assert_eq!(compressed.performance.bank().system().components().len(),2);
        assert!(uniform_stiffness(&compressed) > 0.0);
        let flow = PERFORMANCE.replace("cavity none", "cylinder-gas-loading 0.01 4 2\ncavity none");
        let mut fluid = ready(cx,gate,source(&flow),BODY,700.0);
        let source = fluid.performance.bank().source().bank();
        let cavity = source.spatial_cavity().unwrap();
        assert_eq!(source.cylinder_gas_loading().unwrap().strings()[0].flow().gas(),cavity.gas());
        let count = cavity.coupling().cavity_modes();
        assert_eq!(fluid.performance.bank().system().components().iter()
            .map(|c| c.modes().len()).collect::<Vec<_>>(),vec![4,1,count-1,8]);
        assert_eq!(fluid.performance.bank().system().mode_count(),4+count+8);
        pressure(&mut fluid,37);
        let frame = fluid.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn strict_volume_and_whole_network_admission_precede_material_reads() {
    for (from,to) in [
        (RESOLUTION,"interior-fe-v1 0 1100 24 256 80 10000000 1e-8"),
        (RESOLUTION,"interior-fe-v1 2 1100 24 204 80 10000000 1e-8"),
        (RESOLUTION,"interior-fe-v1 2 1100 0 256 80 10000000 1e-8"),
        (RESOLUTION,"interior-fe-v1 2 1100 24 256 80 10000000 0"),
        ("air-tet 42 0 12 14","air-tet 42 0 12 43"),
        ("air-tet 42 0 12 14","air-tet 42 0 12 12"),
        ("air-tet 42 11 13 12","air-tet 42 0 12 14"),
    ] {
        assert!(body::Profile::parse(BODY.replace(from,to).as_bytes()).is_err(),"{to}");
    }
    assert!(body::Profile::parse(format!("{BODY}ignored\n").as_bytes()).is_err());
    scoped(|cx,gate| {
        let mut budget = source(PERFORMANCE);
        budget.prepared_linear = Some(PreparedModalBudget {
            max_modes: 27, max_dense_bytes: 67_108_864, max_multiply_terms: 200_000_000_000,
        });
        let duplicate = source(&PERFORMANCE.replace("cavity none","cavity neck 0.003 0.01 0.02 100 1"));
        for spec in [budget,duplicate] {
            let mut reads = 0;
            assert!(spec.compile_with_body::<NormalizedPack>(cx,gate,Some(body::Profile::parse(BODY.as_bytes()).unwrap()),
                |_| { reads += 1; Err("material read after inadmissible source".into()) }).is_err());
            assert_eq!(reads,0);
        }
        let open = BODY.replace("pressure-surfaces 80 80","pressure-surfaces 80 79")
            .replace("cavity-facet 79 positive\n","");
        let error = source(PERFORMANCE).compile_with_body(cx,gate,
            Some(body::Profile::parse(open.as_bytes()).unwrap()),|name| Ok(
                if name == "closed-shell.fsmatpk" { shell_pack(name,700.0) }
                else { pack(name,7800.0,200_000.0,None) }))
            .err().expect("an open cavity must not acquire an invented closing wall");
        assert!(error.to_string().contains("material shell cavity"));
    });
}


