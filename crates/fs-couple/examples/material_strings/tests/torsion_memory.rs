//! Source data drives torsional history through the existing performer and PCM.
use super::*;

const MEMORY: &str = "torsion-memory 2\ntorsion-layer shear_modulus 0\ntorsion-layer shear_modulus 2 coating_delta_1 coating_tau_1 coating_delta_2 coating_tau_2\n";
fn base() -> String {
    include_str!("../longitudinal-cylinder-bow.material-strings")
        .replace("48000 4096 0.01 128", "48000 129 0.01 128")
        .replace("prepared-linear 64 67108864 200000000000\n", "")
}
fn selected() -> String {
    base().replace("cylinder-gas-loading ", &format!("{MEMORY}cylinder-gas-loading "))
}
fn coating(name: &str, tau: f64) -> NormalizedPack {
    let mut claims = pack(name,7800.0,200_000.0,None).claims().clone();
    let artifact = fs_blake3::ContentHash([61;32]);
    let provenance = Provenance { source:"manufactured torsional spectrum, not measured wire data".into(),
        license:"CC0-1.0".into(), artifact:Some(artifact) };
    let observation = claims.register_observation(ObservationDataset {
        specimen:name.into(), method:"synthetic shear relaxation".into(), artifact,
        caveats:"not measured".into(), provenance:provenance.clone(),
    }).unwrap();
    for (name,value,dims) in [
        ("coating_delta_1",12e9,Pressure::DIMS),
        ("coating_tau_1",0.0001,Dims([0,0,1,0,0,0])),
        ("coating_delta_2",28e9,Pressure::DIMS),
        ("coating_tau_2",tau,Dims([0,0,1,0,0,0])),
    ] {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar { value,dims },
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
            uncertainty:UncertaintyModel::Unstated, interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![observation], provenance:provenance.clone(),
        }).unwrap();
    }
    let pack = NormalizedPack::new(name,"synthetic-torsion-memory-v1",artifact,
        "synthetic redistribution permitted",claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn build(cx:&Cx<'_>, gate:&CancelGate, text:&str, tau:f64) -> Result<Ready> {
    input::Spec::parse(text.as_bytes())?.compile(cx,gate,|name| Ok(match name {
        "sleeve.fsmatpk" => coating(name,tau),
        // Neither the elastic core nor the body needs coating memory properties.
        _ => pack(name,7800.0,200_000.0,None),
    }))
}

#[test]
fn coating_shear_memory_reaches_shared_dynamics_pressure_and_partition_independent_pcm() {
    scoped(|cx,gate| {
        let text = selected();
        let mut a = build(cx,gate,&text,0.0004).unwrap();
        let mut b = build(cx,gate,&text,0.0012).unwrap();
        let source = a.performance.bank().source().bank();
        assert!(source.cylinder_gas_loading().is_some());
        assert!(source.longitudinal_waves().is_some());
        assert!(source.material_cavity().is_some());
        // One angular plane: transverse polarization must not duplicate memory.
        assert_eq!(a.performance.bank().system().material_memory().len(),2);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_eq!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_ne!(a.catalog["sleeve.fsmatpk"].pack.content_hash(),b.catalog["sleeve.fsmatpk"].pack.content_hash());
        let x = pressure(&mut a,37); let y = pressure(&mut b,128);
        assert!(x.iter().any(|p|p.abs()>1e-15));
        assert!(x.iter().zip(&y).any(|(x,y)|(x-y).abs()>1e-18));
        assert!(a.performance.bank().system().material_memory().iter().any(|x|*x != 0.0));
        assert_eq!(y,pressure(&mut build(cx,gate,&text,0.0012).unwrap(),1));
        let frame = b.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        let (wav,summary) = build(cx,gate,&text,0.0004).unwrap()
            .write(Cursor::new(Vec::new()),gate).unwrap();
        assert_eq!(wav.into_inner(),fs_couple::pcm_wav::encode_pcm16_wav(&x,48000,0.01).unwrap().0);
        assert_eq!(summary.samples,129);
    });
}

#[test]
fn explicitly_elastic_torsional_layers_preserve_existing_arithmetic_without_dummy_data() {
    scoped(|cx,gate| {
        let elastic = selected().replace(MEMORY,
            "torsion-memory 2\ntorsion-layer shear_modulus 0\ntorsion-layer shear_modulus 0\n");
        let mut a = build(cx,gate,&base(),0.0004).unwrap();
        let mut b = input::Spec::parse(elastic.as_bytes()).unwrap()
            .compile(cx,gate,|name|Ok(pack(name,7800.0,200_000.0,None))).unwrap();
        assert!(b.performance.bank().system().material_memory().is_empty());
        assert_eq!(pressure(&mut a,37),pressure(&mut b,1));
    });
}

#[test]
fn torsional_memory_structure_and_joint_history_budget_refuse_before_material_io() {
    let text = selected();
    for record in [
        "torsion-memory 0\n",
        "torsion-memory 1\ntorsion-layer shear_modulus 0\n",
        "torsion-memory 2\ntorsion-layer shear_modulus 0\ntorsion-layer shear_modulus 1 coating_delta_1 shear_modulus\n",
        "torsion-memory 2\ntorsion-layer shear_modulus 0\ntorsion-layer shear_modulus 2 coating_delta_1 coating_tau_1 coating_delta_1 coating_tau_1\n",
    ] {
        assert!(input::Spec::parse(text.replace(MEMORY,record).as_bytes()).is_err());
    }
    assert!(input::Spec::parse(text.replace("tension 20 1 1 0.8","tension 20 1 0 0.8").as_bytes()).is_err());
    assert!(input::Spec::parse(text.replace("longitudinal-waves ",
        "prepared-linear 64 67108864 200000000000\nlongitudinal-waves ").as_bytes()).is_ok());
    scoped(|cx,gate| {
        let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
        spec.strings[0].twist_modes=128;
        spec.strings[0].torsion_memory.as_mut().unwrap()[0].branches
            .push(("core_delta".into(),"core_tau".into()));
        let mut reads=0;
        assert!(spec.compile(cx,gate,|name| { reads+=1; Ok(pack(name,7800.0,200_000.0,None)) }).is_err());
        assert_eq!(reads,0);
        let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
        spec.prepared_linear=Some(PreparedModalBudget {max_modes:1,..MEMORY_PREPARATION});
        let mut reads=0;
        assert!(spec.compile(cx,gate,|name| { reads+=1; Ok(pack(name,7800.0,200_000.0,None)) }).is_err());
        assert_eq!(reads,0);
    });
}

#[test]
fn missing_or_invalid_shear_memory_has_no_elastic_fallback() {
    scoped(|cx,gate| {
        assert!(build(cx,gate,&selected().replace("coating_tau_2","absent_time"),0.0004).is_err());
        assert!(build(cx,gate,&selected(),0.0).is_err());
        assert!(build(cx,gate,&selected(),-0.001).is_err());
    });
}

#[test]
fn prepared_torsional_memory_keeps_the_original_coupled_history_and_pressure() {
    scoped(|cx,gate| {
        let text=selected();
        let mut spec=input::Spec::parse(text.as_bytes()).unwrap();
        spec.prepared_linear=Some(MEMORY_PREPARATION);
        let prepared=spec.compile(cx,gate,|name|Ok(match name {
            "sleeve.fsmatpk"=>coating(name,0.0004),
            _=>pack(name,7800.0,200_000.0,None),
        })).unwrap();
        assert_eq!(prepared.performance.bank().system().material_memory().len(),2);
        assert_prepared_memory_parity(build(cx,gate,&text,0.0004).unwrap(),prepared);
    });
}
