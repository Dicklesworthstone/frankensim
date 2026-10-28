//! Native tests execute the real pack decoder, material reducers and pressure writer.
use super::*;
use std::io::Cursor;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,ObservationDataset,PropertyClaim,PropertyKey,PropertyValue,Provenance,InterpolationPolicy,UncertaintyModel};
use fs_scenario::gesture::{GestureTrack,GestureEvent};

#[path = "tests/normal_sources.rs"]
mod normal_sources;

#[path = "tests/material_cards.rs"]
mod material_cards;

#[path = "tests/fibers.rs"]
mod fibers;

#[path = "tests/bow_motion.rs"]
mod bow_motion;

#[path = "tests/attached_sources.rs"]
mod attached_sources;

#[path = "tests/shell_radiation.rs"]
mod shell_radiation;

#[path = "tests/observer.rs"]
mod shell_observer;

#[path = "tests/one_way_observer.rs"]
mod one_way_observer;

fn pack(name:&str,rho:f64,band:f64,omit:Option<&str>)->NormalizedPack {
    let mut claims=ClaimSet::new();
    let provenance=Provenance {source:"synthetic executable test".into(),license:"CC0-1.0".into(),artifact:Some(fs_blake3::ContentHash([7;32]))};
    let obs=claims.register_observation(ObservationDataset {
        specimen:name.into(),method:"synthetic executable test".into(),
        artifact:fs_blake3::ContentHash([7;32]),caveats:"synthetic test fixture".into(),
        provenance:provenance.clone(),
    }).unwrap();
    for (key,dims,value) in [("density",Density::DIMS,rho),("young_modulus",Pressure::DIMS,200e9),
        ("poisson_ratio",Dims::NONE,0.3),(ETA,DynViscosity::DIMS,10.0),
        (SHEAR,Pressure::DIMS,80e9),(SHEAR_ETA,DynViscosity::DIMS,1000.0)] {
        if omit==Some(key) {continue;}
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0).with("omega",0.0,band),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![obs.clone()],provenance:provenance.clone()}).unwrap();
    }
    let pack=NormalizedPack::new(name,"synthetic-test-v1",fs_blake3::ContentHash([7;32]),
        "synthetic redistribution permitted",claims,vec![],vec![]).unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}
fn scoped<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&Cx::new(&gate,arena,fs_exec::StreamKey {seed:0,kernel_id:0,tile:0,iteration:0},
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn spec()->input::Spec {
    let bytes=include_str!("pluck.material-strings").replace("48000 4096 0.01 128","48000 256 0.01 128");
    input::Spec::parse(bytes.as_bytes()).unwrap()
}
fn compile(cx:&Cx<'_>,gate:&CancelGate,s:input::Spec,rho:f64)->Ready {
    s.compile(cx,gate,|name|Ok(pack(name,if name=="wire.fsmatpk"{rho}else{7800.0},200_000.0,None))).unwrap()
}
fn pressure(ready:&mut Ready,block:usize)->Vec<f64> {
    let mut out=vec![0.0;ready.spec.samples as usize];
    for chunk in out.chunks_mut(block) {ready.performance.block(chunk).unwrap();}out
}
#[test]
fn material_pack_to_physical_pressure_and_wav_uses_the_same_sample_history() {
    scoped(|cx,gate| {
        let mut direct=compile(cx,gate,spec(),7800.0);
        assert_eq!(direct.catalog.len(),2);
        assert!(direct.performance.bank().source().bank().material_cavity().is_some());
        assert!(direct.performance.preload_energy_j().unwrap()>0.0);
        let expected=pressure(&mut direct,37);assert!(expected.iter().any(|p|p.abs()>1e-10));
        let (encoded,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&expected,48000,0.01).unwrap();
        let mut actual=compile(cx,gate,spec(),7800.0);
        let (output,summary)=actual.write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        assert_eq!(output.into_inner(),encoded);assert_eq!(summary.samples,256);
        assert_eq!(summary.clipped_samples,clips as u64);
        assert_eq!(actual.performance.samples_rendered(),256);
    });
}
#[test]
fn source_density_swap_changes_modes_and_actual_pressure_without_changing_geometry_or_body() {
    scoped(|cx,gate| {
        let mut a=compile(cx,gate,spec(),7800.0);let mut b=compile(cx,gate,spec(),15600.0);
        assert_eq!(a.spec.bytes,b.spec.bytes);
        assert_ne!(a.catalog["wire.fsmatpk"].pack.content_hash(),b.catalog["wire.fsmatpk"].pack.content_hash());
        assert_eq!(a.performance.bank().system().components()[1].modes(),b.performance.bank().system().components()[1].modes());
        assert!(b.performance.bank().cards()[0].fundamental_hz()<a.performance.bank().cards()[0].fundamental_hz());
        let x=pressure(&mut a,128);let y=pressure(&mut b,37);
        assert!(x.iter().zip(y).any(|(x,y)|(x-y).abs()>1e-12));
    });
}
#[test]
fn gas_pressure_is_shared_by_strings_cavity_and_actual_output_not_a_post_effect() {
    scoped(|cx,gate| {
        let changed=String::from_utf8(spec().bytes).unwrap().replace(
            "gas mixture 293.15 101325", "gas mixture 293.15 60000");
        let changed=input::Spec::parse(changed.as_bytes()).unwrap();
        let mut a=compile(cx,gate,spec(),7800.0);let mut b=compile(cx,gate,changed,7800.0);
        assert_eq!(b.performance.bank().source().bank().gas(),b.performance.bank().source().bank().material_cavity().unwrap().cavity.gas());
        assert_ne!(a.performance.bank().cards()[0].zetas,b.performance.bank().cards()[0].zetas);
        let x=pressure(&mut a,128);let y=pressure(&mut b,128);
        assert!(x.iter().zip(y).any(|(x,y)|(x-y).abs()>1e-12));
    });
}
#[test]
fn canonical_double_stop_ramps_and_pluck_share_the_original_joint_owner() {
    scoped(|cx,gate| {
        let mut s=spec();s.strings.push(s.strings[0].clone());s.strings[1].prestress=StringPrestress::FixedTension(22.0);
        let track=|id:&str,string,station|GestureTrack {id:id.into(),target:GestureTarget::BowStroke {string},
            initial:GestureValue::Bow {velocity_m_per_s:0.0,normal_force_n:0.0,station},events:vec![
                GestureEvent {time_s:0.0,transition_s:0.002,value:GestureValue::Bow {velocity_m_per_s:0.01,normal_force_n:0.0001,station}},
                GestureEvent {time_s:0.003,transition_s:0.0,value:GestureValue::Bow {velocity_m_per_s:0.0,normal_force_n:0.0,station}}]};
        s.schedule=GestureSchedule::try_new(1000,vec![track("first",0,0.2),track("second",1,0.31)]).unwrap();
        // Exercise the public source grammar too; canonical bytes are generated
        // by its sole encoder, not a handwritten competing schedule serializer.
        let header=String::from_utf8(s.bytes.clone()).unwrap().split_once("\nschedule\n").unwrap().0
            .replace("strings 1\nstring wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1",
                "strings 2\nstring wire.fsmatpk 0.5 0.0003 20 4 3 0.8 3000 0.1\nstring wire.fsmatpk 0.5 0.0003 22 4 3 0.8 3000 0.1");
        let mut bytes=format!("{header}\nschedule\n").into_bytes();bytes.extend(s.schedule.to_canonical_bytes());
        let mut a=compile(cx,gate,input::Spec::parse(&bytes).unwrap(),7800.0);
        let mut b=compile(cx,gate,input::Spec::parse(&bytes).unwrap(),7800.0);
        assert_eq!(a.catalog.len(),2,"same immutable pack should be loaded once");
        assert_eq!(a.performance.bank().source().bows().len(),2);
        let cols=a.performance.projected_string_forces().unwrap().initial_forces()[0].clone();
        assert!(cols[7..].iter().all(|f|*f==0.0),"other string/body must not acquire a duplicate direct load");
        assert_eq!(pressure(&mut a,1),pressure(&mut b,128));
        let last=a.performance.bank().system().last_frame().unwrap();
        assert!(last.energy_residual_j.abs()<=last.energy_tolerance_j);
        assert!(last.network_energy_j>0.0);
    });
}
#[test]
fn missing_loss_band_and_cancel_refuse_instead_of_constructing_a_fallback_material() {
    scoped(|cx,gate| {
        for (band,omit) in [(200_000.0,Some(ETA)),(200_000.0,Some(SHEAR_ETA)),(2000.0,None)] {
            assert!(spec().compile(cx,gate,|name|Ok(pack(name,7800.0,band,omit))).is_err());
        }
        let stop=CancelGate::new_clock_free();stop.request();let mut calls=0;
        assert!(spec().compile(cx,&stop,|name|{calls+=1;Ok(pack(name,7800.0,200_000.0,None))}).is_err());
        assert_eq!(calls,0);
    });
}
#[test]
fn malformed_or_unobserved_controls_and_unbounded_decoder_counts_refuse() {
    let text=include_str!("pluck.material-strings");
    for (from,to) in [("strings 1","strings 999999999"),("4 3 0.8","128 128 0.8"),
        ("wire.fsmatpk","../wire.fsmatpk"),("force 0 0 0","force 4096 0 0"),
        ("force 0 0 0","force 0 1 0"),("events\t0","events\t18446744073709551615"),
        ("tracks\t1","tracks\t0"),("initial\tbow\t0e0\t0e0\t2e-1","initial\tbow\t0e0\t0e0\t2e-1\textra"),
        ("O2 0.2","unknown 0.2"),("audio 48000","audio 0")] {
        assert!(input::Spec::parse(text.replace(from,to).as_bytes()).is_err(),"accepted {to}");
    }
    assert!(input::Spec::parse(format!("{text}unconsumed\n").as_bytes()).is_err());
    assert!(input::Spec::parse(include_bytes!("pluck-and-bow.material-strings")).is_ok());
}

#[path = "tests/layered.rs"]
mod layered;

#[path = "tests/tension_loss.rs"]
mod tension_loss;

#[path = "tests/bow_trace.rs"]
mod bow_trace;

#[path = "tests/residual.rs"]
mod residual;

#[path = "tests/polarized.rs"]
mod polarized;

#[path = "stops/tests.rs"]
mod stopped;

#[path = "tests/fluid.rs"]
mod fluid;

#[path = "tests/prepared.rs"]
mod prepared;

#[path = "stops/motion_tests.rs"]
mod moving_stops;

#[path = "stops/contact_load_tests.rs"]
mod contact_loaded_bow;

#[path = "tests/nonlinear.rs"]
mod nonlinear;

#[path = "tests/longitudinal.rs"]
mod longitudinal;

#[path = "tests/bow_interfaces.rs"]
mod bow_interfaces;

#[path = "tests/bridge_interfaces.rs"]
pub(crate) mod bridge_interfaces;

#[path = "tests/thermochemical.rs"]
mod thermochemical_gas;

#[path = "tests/cavity.rs"]
mod physical_cavity;

#[path = "tests/spatial_cavity.rs"]
mod spatial_cavity;

#[path = "stops/station_tests.rs"]
mod migrating_stops;

#[path = "tests/thermal_support.rs"]
mod thermal_support;

#[path = "tests/relaxation.rs"]
mod relaxation;

#[path = "tests/torsion_memory.rs"]
mod torsion_memory;

#[path = "tests/wound.rs"]
mod wound;

#[path = "tests/shear_layers.rs"]
mod shear_layers;
