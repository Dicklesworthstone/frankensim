//! Actual pack -> material tensor -> shared instrument -> pressure/PCM tests.
//! All coefficients below are synthetic; they are not measured wood or metal.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,ObservationDataset,PropertyClaim,PropertyKey,PropertyValue,
    Provenance,InterpolationPolicy,UncertaintyModel};

fn pack(name:&str, viscosity:f64, omit:Option<&str>) -> NormalizedPack {
    let body=name=="body.fsmatpk";
    let mut values=vec![("density",Density::DIMS,if body {600.0}else{7800.0}),
        ("young_modulus",Pressure::DIMS,if body {8e9}else{200e9}),
        ("poisson_ratio",Dims::NONE,0.2),(ETA,DynViscosity::DIMS,10.0),
        (SHEAR,Pressure::DIMS,80e9),(SHEAR_ETA,DynViscosity::DIMS,1000.0),
        ("young_modulus_1",Pressure::DIMS,8e9),("young_modulus_2",Pressure::DIMS,2e9),
        ("poisson_ratio_12",Dims::NONE,0.2),("shear_modulus_12",Pressure::DIMS,0.6e9)];
    for (key,value) in VISCOSITY.into_iter().zip([3e5,6e4,4e4,1e4,4e3,-3e3]) {
        values.push((key,DynViscosity::DIMS,value*viscosity));
    }
    let provenance=Provenance {source:"synthetic tensor body test".into(),license:"CC0-1.0".into(),
        artifact:Some(fs_blake3::ContentHash([8;32]))};
    let mut claims=ClaimSet::new();
    let obs=claims.register_observation(ObservationDataset {specimen:name.into(),method:"synthetic".into(),
        artifact:fs_blake3::ContentHash([8;32]),caveats:"not measured".into(),provenance:provenance.clone()}).unwrap();
    for (key,dims,value) in values {
        if body && omit==Some(key) {continue;}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,200_000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![obs.clone()],provenance:provenance.clone()}).unwrap();
    }
    let p=NormalizedPack::new(name,"synthetic-tensor-v1",fs_blake3::ContentHash([8;32]),
        "synthetic redistribution permitted",claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&p.to_bytes()).unwrap()
}
fn scoped<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|a|f(&Cx::new(&gate,a,fs_exec::StreamKey {seed:0,kernel_id:0,tile:0,iteration:0},
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn spec()->input::Spec {
    let mut s=input::Spec::parse(include_bytes!("../pluck.material-strings")).unwrap();
    s.samples=256;s.body.modes=3;s.body.thickness=0.002;s.block=128;s
}
fn profile()->Profile {Profile::parse(include_bytes!("../orthotropic.body")).unwrap()}
fn compile(cx:&Cx<'_>,gate:&CancelGate,p:Option<Profile>,scale:f64)->Ready {
    spec().compile_with_body(cx,gate,p,|name|Ok(pack(name,scale,None))).unwrap()
}
fn pressure(r:&mut Ready,n:usize)->Vec<f64> {
    let mut p=vec![0.0;r.spec.samples as usize];
    for block in p.chunks_mut(n) {r.performance.block(block).unwrap();}p
}

#[test]
fn no_profile_keeps_the_original_material_performance_path() {
    scoped(|cx,gate| {
        let mut a=spec().compile(cx,gate,|n|Ok(pack(n,1.0,None))).unwrap();
        let mut b=compile(cx,gate,None,1.0);
        assert!(a.body_profile.is_none());
        assert!(a.performance.bank().source().bank().material_plate().is_some());
        assert_eq!(pressure(&mut a,37),pressure(&mut b,128));
    });
}
#[test]
fn directional_body_and_full_source_loss_reach_actual_pressure_and_pcm() {
    scoped(|cx,gate| {
        let mut a=compile(cx,gate,Some(profile()),1.0);
        let source=a.performance.bank().source().bank();
        let chart=source.material_chart().unwrap();
        assert_eq!(chart.specimen.regions()[0].input.model,PlateMaterialModel::Orthotropic12);
        assert!(chart.specimen.bending_viscosity().is_none());
        let image=source.material_chart_viscosity().unwrap();
        let c=image.damping().matrix_per_s();let n=image.damping().mode_count();
        assert_eq!(n,3);
        assert!(c.iter().enumerate().any(|(i,v)|i/n!=i%n && v.abs()>1e-8));
        assert!(image.radiators().iter().all(|r|r.zeta==0.0));
        assert!(source.body_connections().iter().any(|l|l.stiffness_n_m==0.0 && l.damping_n_s_m>0.0));
        assert_eq!(a.body_profile.as_ref().unwrap().bytes,profile().bytes);
        let p=pressure(&mut a,37);assert!(p.iter().any(|v|v.abs()>1e-12));
        let mut b=compile(cx,gate,Some(profile()),1.0);
        let expected=fs_couple::pcm_wav::encode_pcm16_wav(&p,b.spec.rate,b.spec.full_scale).unwrap().0;
        let (actual,report)=b.write(std::io::Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(actual.into_inner(),expected);assert_eq!(report.samples,256);
    });
}
#[test]
fn grain_and_viscosity_are_independent_data_only_changes_not_retuning() {
    scoped(|cx,gate| {
        let mut a=compile(cx,gate,Some(profile()),1.0);
        let turned=Profile::parse(include_str!("../orthotropic.body").replace("grain 0.17","grain 0.87").as_bytes()).unwrap();
        let mut b=compile(cx,gate,Some(turned),1.0);
        let mut c=compile(cx,gate,Some(profile()),8.0);
        assert_eq!(a.performance.bank().cards()[0].fundamental_hz(),b.performance.bank().cards()[0].fundamental_hz());
        let modes=|r:&Ready|r.performance.bank().system().components()[1].modes().to_vec();
        assert_ne!(modes(&a),modes(&b));assert_eq!(modes(&a),modes(&c));
        let pa=pressure(&mut a,1);let pb=pressure(&mut b,128);let pc=pressure(&mut c,37);
        assert!(pa.iter().zip(pb).any(|(a,b)|(a-b).abs()>1e-12));
        assert!(pa.iter().zip(pc).any(|(a,b)|(a-b).abs()>1e-12));
        for r in [&a,&b,&c] {
            let f=r.performance.bank().system().last_frame().unwrap();
            assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
        }
    });
}
#[test]
fn absent_tensor_coefficients_refuse_even_with_legacy_bending_loss_present() {
    scoped(|cx,gate| {
        for key in VISCOSITY {
            assert!(spec().compile_with_body(cx,gate,Some(profile()),|n|Ok(pack(n,1.0,Some(key)))).is_err());
        }
        let mut p=profile();p.window=(1e12,2e12);
        assert!(spec().compile_with_body(cx,gate,Some(p),|n|Ok(pack(n,1.0,None))).is_err());
        let stopped=CancelGate::new_clock_free();stopped.request();let mut calls=0;
        assert!(spec().compile_with_body(cx,&stopped,Some(profile()),|n| {calls+=1;Ok(pack(n,1.0,None))}).is_err());
        assert_eq!(calls,0);
    });
}
#[test]
fn profile_rejects_unbounded_mesh_silent_fields_and_fixed_boundary_bridge() {
    let text=include_str!("../orthotropic.body");
    for (from,to) in [("mesh 4 4","mesh 99999 4"),("mesh 4 4","mesh 1 4"),
        ("grain 0.17","grain NaN"),("window 1 100000000","window 3 2"),
        ("bridge 2 2","bridge 0 2"),("bridge 2 2","bridge 2 4"),
        ("elastic orthotropic12","elastic wood"),("pretension 0","pretension -1"),
        ("grain 0.17","grain 0.17 ignored")] {
        assert!(Profile::parse(text.replace(from,to).as_bytes()).is_err(),"accepted {to}");
    }
    assert!(Profile::parse(format!("{text}extra\n").as_bytes()).is_err());
}
