//! Source-profile -> solved normal load -> existing pressure/PCM history.
use super::*;
use crate::stops::Profile;

fn input_text()->String {
    include_str!("../contact-bow.material-strings").replace("48000 4096", "48000 192")
}
fn profile_text()->String {
    include_str!("../contact-bow.stops").replace("knot 480", "knot 32")
        .replace("knot 960", "knot 64").replace("knot 1440", "knot 128")
}
fn ready(cx:&Cx<'_>,gate:&CancelGate,profile:&str,rho:f64)->Ready {
    input::Spec::parse(input_text().as_bytes()).unwrap().compile_with_profiles(
        cx,gate,None,Some(Profile::parse(profile.as_bytes()).unwrap()),
        |name|Ok(pack(name,if name=="wire.fsmatpk"{rho}else{7800.0},200000.0,None))).unwrap()
}

#[test]
fn file_driven_normal_contact_sets_bow_traction_and_original_pressure_pcm() {
    scoped(|cx,gate| {
        let text=profile_text();let mut a=ready(cx,gate,&text,7800.0);
        let binding=a.stop_profile.as_ref().unwrap().bow_contacts()[0];
        assert_eq!(binding.bow,0);assert_eq!(binding.stop,0);
        assert_eq!(a.performance.bank().source().bow_contact_loads(),&[binding]);
        assert_eq!(a.stop_profile.as_ref().unwrap().bytes,text.as_bytes());
        assert!(a.performance.bank().source().bank().material_cavity().is_some());
        let (mut active,mut separated)=(0,0);let mut out=Vec::new();
        for _ in 0..192 {
            let mut sample=[0.0];a.performance.block(&mut sample).unwrap();out.push(sample[0]);
            let frame=a.performance.bank().system().last_frame().unwrap();
            let n=a.performance.bank().system().last_obstacle_frames().unwrap()[0].normal_force_n;
            assert_eq!(frame.ports[0].normal_force_n,n);
            if n==0.0 {separated+=1;assert_eq!(frame.ports[0].traction_n,0.0);} else {active+=1;}
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        }
        assert!(active>0 && separated>0 && out.iter().any(|p|p.abs()>1e-12));
        let mut b=ready(cx,gate,&text,7800.0);
        assert_eq!(pressure(&mut b,37),out);
        let mut c=ready(cx,gate,&text,7800.0);
        let (bytes,summary)=c.write(Cursor::new(Vec::new()),gate).unwrap();
        let (expected,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&out,c.spec.rate,c.spec.full_scale).unwrap();
        assert_eq!(bytes.into_inner(),expected);assert_eq!(summary.clipped_samples,clips as u64);
    });
}

#[test]
fn normal_binding_and_material_data_change_actual_sound_not_the_observer_gain() {
    scoped(|cx,gate| {
        let text=profile_text();
        let mut a=ready(cx,gate,&text,7800.0);
        let mut b=ready(cx,gate,&text.replace("bow-contact 0 bow_surface\n",""),7800.0);
        let mut c=ready(cx,gate,&text,15600.0);
        for (a,b) in a.performance.bank().system().components().iter().zip(b.performance.bank().system().components()) {
            assert_eq!(a.modes(),b.modes());
        }
        assert_eq!(a.performance.bank().system().acceleration_pressure_weights(),
            b.performance.bank().system().acceleration_pressure_weights());
        assert_eq!(a.performance.bank().system().components()[1].modes(),
            c.performance.bank().system().components()[1].modes());
        let pa=pressure(&mut a,128);let pb=pressure(&mut b,128);let pc=pressure(&mut c,128);
        assert!(pa.iter().zip(&pb).any(|(a,b)|(a-b).abs()>1e-12));
        assert!(pa.iter().zip(&pc).any(|(a,c)|(a-c).abs()>1e-12));
        assert_ne!(a.catalog["wire.fsmatpk"].pack.content_hash(),c.catalog["wire.fsmatpk"].pack.content_hash());
    });
}

#[test]
fn conflicting_or_unmatched_bow_contacts_refuse_before_pack_io() {
    let text=profile_text();
    for replacement in ["bow-contact 0 missing", "bow-contact 16 bow_surface",
        "bow-contact 0 bow_surface\nbow-contact 0 bow_surface"] {
        assert!(Profile::parse(text.replace("bow-contact 0 bow_surface",replacement).as_bytes()).is_err());
    }
    scoped(|cx,gate| {
        for (profile,input) in [
            (text.replace("bow-contact 0", "bow-contact 3"),input_text()),
            (text.replace("0.2 normal", "0.31 normal"),input_text()),
            (text.clone(),input_text().replace("3e-1\t0e0", "3e-1\t1e-3")),
        ] {
            let mut calls=0;
            let result=input::Spec::parse(input.as_bytes()).unwrap().compile_with_profiles(
                cx,gate,None,Some(Profile::parse(profile.as_bytes()).unwrap()),|name| {
                    calls+=1;Ok(pack(name,7800.0,200000.0,None))
                });
            assert!(result.is_err());assert_eq!(calls,0);
        }
        // Physical orientation is checked again against the realized source,
        // never repaired by flipping the existing torsional tangent map.
        assert!(input::Spec::parse(input_text().as_bytes()).unwrap().compile_with_profiles(
            cx,gate,None,Some(Profile::parse(text.replace("normal 0 -1", "normal 0 1").as_bytes()).unwrap()),
            |name|Ok(pack(name,7800.0,200000.0,None))).is_err());
    });
}
