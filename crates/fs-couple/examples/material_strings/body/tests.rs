//! Actual pack -> material tensor -> shared instrument -> pressure/PCM tests.
//! All coefficients below are synthetic; they are not measured wood or metal.
use super::*;
#[path = "rayleigh_tests.rs"]
mod rayleigh_tests;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,ObservationDataset,PropertyClaim,PropertyKey,PropertyValue,
    Provenance,InterpolationPolicy,UncertaintyModel};

fn pack(name:&str, viscosity:f64, omit:Option<&str>) -> NormalizedPack {
    pack_with_shear(name, viscosity, 1.0, omit)
}
fn pack_with_shear(name:&str, viscosity:f64, shear_scale:f64, omit:Option<&str>) -> NormalizedPack {
    let body=name=="body.fsmatpk";
    let mut values=vec![("density",Density::DIMS,if body {600.0}else{7800.0}),
        ("young_modulus",Pressure::DIMS,if body {8e9}else{200e9}),
        ("poisson_ratio",Dims::NONE,0.2),(ETA,DynViscosity::DIMS,10.0),
        (SHEAR,Pressure::DIMS,80e9),(SHEAR_ETA,DynViscosity::DIMS,1000.0),
        ("young_modulus_1",Pressure::DIMS,8e9),("young_modulus_2",Pressure::DIMS,2e9),
        ("poisson_ratio_12",Dims::NONE,0.2),("shear_modulus_12",Pressure::DIMS,0.6e9),
        ("shear_modulus_13",Pressure::DIMS,0.5e9*shear_scale),
        ("shear_modulus_23",Pressure::DIMS,0.25e9*shear_scale)];
    for (key,value) in fs_couple::shear_shell::TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES
        .into_iter().zip([9000.0,7000.0,0.0]) {
        values.push((key,DynViscosity::DIMS,value*viscosity));
    }
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
        let mut p=profile();
        let BodyModel::Plate(plate) = &mut p.model else { panic!("expected plate profile"); };
        plate.window=(1e12,2e12);
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

const ROTATIONAL_BRIDGE_ROW: &str = "rotational-bridge 0 0.8 0 1 0.001 1e-7 0";
const ROTATIONAL_BRIDGE_BLOCK: &str = "rotational-bridges 1\nrotational-bridge 0 0.8 0 1 0.001 1e-7 0\n";
fn rotational_input() -> String {
    include_str!("../torsional-bridge.material-strings")
        .replace("audio 48000 4096 0.01 128", "audio 48000 256 0.01 128")
}
fn rotational_profile() -> Profile {
    Profile::parse(include_bytes!("../torsional-bridge.body")).unwrap()
}

#[test]
fn rotational_bridge_file_drives_released_body_with_original_material_loss_and_pcm() {
    scoped(|cx, gate| {
        // The only initial energy is this explicitly declared small-angle
        // spring preload. Every bow is lifted and there are no external loads.
        let unforced = rotational_input().replace(
            "forces preload 1 1\nport 0 0.2 surface 0.002\nforce 0 0 0",
            "forces retain 0 0",
        );
        let active = unforced.replace(ROTATIONAL_BRIDGE_ROW,
            "rotational-bridge 0 0.8 0 1 0.001 1e-7 0.002");
        let omitted = unforced.replace(ROTATIONAL_BRIDGE_BLOCK, "");
        let s = input::Spec::parse(active.as_bytes()).unwrap();
        assert!(s.forces.is_none());
        let p = rotational_profile();
        assert_eq!(p.bridge_separation_m(&s, 0, 1).unwrap(), 0.2);
        let mut connected = s.compile_with_body(cx, gate, Some(p.clone()),
            |name| Ok(pack(name, 1.0, None))).unwrap();
        let mut silent = input::Spec::parse(omitted.as_bytes()).unwrap()
            .compile_with_body(cx, gate, Some(p.clone()), |name| Ok(pack(name, 1.0, None))).unwrap();
        assert_eq!(connected.spec.bytes, active.as_bytes());
        assert_eq!(connected.rotational_bridge_baselines_m, vec![0.2]);
        let source = connected.performance.bank().source().bank();
        let original = silent.performance.bank().source().bank();
        let joint = &source.rotational_body_connections()[0];
        assert_eq!((joint.string, joint.station_fraction, joint.body.component), (0, 0.8, 0));
        assert_eq!((joint.stiffness_n_m_per_rad, joint.damping_n_m_s_per_rad, joint.rest_angle_rad),
            (0.001, 1e-7, 0.002));
        let a = &source.strings()[0].connections[0].body.shapes;
        let b = &source.strings()[1].connections[0].body.shapes;
        assert_eq!(joint.body.angular_shapes, a.iter().zip(b).map(|(a,b)| (a-b)/0.2).collect::<Vec<_>>());
        assert_eq!(source.body_connections(), original.body_connections());
        assert_eq!(source.material_chart_viscosity().unwrap().damping().matrix_per_s(),
            original.material_chart_viscosity().unwrap().damping().matrix_per_s());
        assert_eq!(connected.performance.bank().system().components()[2].modes(),
            silent.performance.bank().system().components()[2].modes());

        let actual_pressure = pressure(&mut connected, 37);
        assert!(actual_pressure.iter().any(|p| p.abs() > 1e-16));
        assert!(pressure(&mut silent, 128).iter().all(|p| *p == 0.0));
        let body = &connected.performance.bank().system().components()[2];
        assert!(body.states().iter().any(|s|
            s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0));
        let frame = connected.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

        // The ordinary file performer's streaming writer consumes exactly the
        // same coupled pressure; the joint is not an output-stage sound effect.
        let mut replay = input::Spec::parse(active.as_bytes()).unwrap()
            .compile_with_body(cx, gate, Some(p), |name| Ok(pack(name, 1.0, None))).unwrap();
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(
            &actual_pressure, replay.spec.rate, replay.spec.full_scale).unwrap().0;
        let (wav, report) = replay.write(std::io::Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(report.samples, 256);
    });
}

#[test]
fn rotational_bridge_records_are_optional_bounded_and_version_independent() {
    let text = rotational_input();
    let s = input::Spec::parse(text.as_bytes()).unwrap();
    assert_eq!(s.rotational_bridges.len(), 1);
    assert!(input::Spec::parse(text.replace(ROTATIONAL_BRIDGE_BLOCK, "").as_bytes())
        .unwrap().rotational_bridges.is_empty());
    assert!(input::Spec::parse(text.replace(ROTATIONAL_BRIDGE_BLOCK, "rotational-bridges 0\n").as_bytes())
        .unwrap().rotational_bridges.is_empty());
    let v2 = text.replace("frankensim-material-strings-v1", "frankensim-material-strings-v2")
        .replace("wire.fsmatpk 0.5 0.0003 20", "wire.fsmatpk 0.5 0.0003 tension 20")
        .replace("wire.fsmatpk 0.5 0.0003 22", "wire.fsmatpk 0.5 0.0003 tension 22")
        .replace("0.8 3000 0.1\n", "0.8 3000 0.1\nsection homogeneous\n");
    assert_eq!(input::Spec::parse(v2.as_bytes()).unwrap().rotational_bridges, s.rotational_bridges);
    for record in [
        "rotational-bridge 2 0.8 0 1 0.001 1e-7 0",
        "rotational-bridge 0 0.8 2 1 0.001 1e-7 0",
        "rotational-bridge 0 0.8 0 2 0.001 1e-7 0",
        "rotational-bridge 0 0 0 1 0.001 1e-7 0",
        "rotational-bridge 0 1 0 1 0.001 1e-7 0",
        "rotational-bridge 0 0.8 0 0 0.001 1e-7 0",
        "rotational-bridge 0 0.8 0 1 -0.001 1e-7 0",
        "rotational-bridge 0 0.8 0 1 0.001 -1e-7 0",
        "rotational-bridge 0 0.8 0 1 0 0 0",
        "rotational-bridge 0 0.8 0 1 0.001 1e-7 NaN",
        "rotational-bridge 0 0.8 0 1 0.001 1e-7 0 ignored",
    ] {
        assert!(input::Spec::parse(text.replace(ROTATIONAL_BRIDGE_ROW, record).as_bytes()).is_err(),
            "accepted {record}");
    }
    for (from, to) in [
        ("rotational-bridges 1", "rotational-bridges 17"),
        ("rotational-bridges 1", "rotational-bridges 2"),
        ("20 4 3 0.8", "20 4 0 0.8"),
    ] {
        assert!(input::Spec::parse(text.replace(from, to).as_bytes()).is_err());
    }
}

#[test]
fn rotational_bridge_geometry_refuses_before_loading_materials() {
    scoped(|cx, gate| {
        let s = input::Spec::parse(rotational_input().as_bytes()).unwrap();
        let text = include_str!("../torsional-bridge.body");
        let coincident = Profile::parse(text.replace("bridge 3 2", "bridge 1 2").as_bytes()).unwrap();
        let missing = Profile::parse(text.replace("bridges 2\nbridge 1 2\nbridge 3 2",
            "bridges 1\nbridge 1 2").as_bytes()).unwrap();
        for profile in [None, Some(profile()), Some(coincident), Some(missing)] {
            let mut loads = 0;
            assert!(s.clone().compile_with_body(cx, gate, profile, |name| {
                loads += 1; Ok(pack(name, 1.0, None))
            }).is_err());
            assert_eq!(loads, 0);
        }
        let p = rotational_profile();
        assert_eq!(p.bridge_separation_m(&s, 0, 1).unwrap(), p.bridge_separation_m(&s, 1, 0).unwrap());
        assert!(p.bridge_separation_m(&s, 0, 2).is_err());
        let mut collapsed = s;
        collapsed.body.length = f64::from_bits(1);
        assert!(p.bridge_separation_m(&collapsed, 0, 1).is_err());
    });
}

#[path="regions.rs"]
mod regions;

#[path="fluid.rs"]
mod fluid;

#[path = "shell/tests.rs"]
mod shell_input;

#[path = "shell/axial_bearing_tests.rs"]
mod axial_bearing_tests;

#[path = "shell/assembly_tests.rs"]
mod assembly_input;

#[path = "relaxation_tests.rs"]
mod material_memory;

#[path = "tensor_memory_tests.rs"]
mod tensor_memory;
