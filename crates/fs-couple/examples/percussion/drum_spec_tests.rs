use super::*;
use crate::{Stroke, drum_with_air, drum_with_spec, mechanics::Mechanics};
use fs_exec::CancelGate;

fn text(s: Spec) -> String {
    let mut out=format!("{HEADER}\ngeometry,{:.17e},{:.17e},{:.17e}\n",
        s.radius_m,s.depth_m,s.outer_radius_m);
    for (name,h) in ["batter","resonant"].into_iter().zip(s.heads) {
        out.push_str(&format!("head,{name},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}\n",
            h.thickness_m,h.young_pa,h.poisson,h.density_kg_m3,h.tension_n_m,h.damping_ratio));
    }
    out.push_str(&format!("mesh,{},{}\nband_hz,{:.17e},{:.17e}\n",
        s.radial_intervals,s.azimuths,s.band_hz[0],s.band_hz[1]));
    out
}
fn small() -> Spec { Spec {radial_intervals:2,azimuths:16,..Spec::reference()} }
fn stroke() -> Stroke {Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])}}

#[test]
fn complete_materials_and_mesh_round_trip_without_defaults() {
    let mut spec=small();
    spec.heads[0].tension_n_m=2456.7;spec.heads[1].young_pa=3.2e9;
    spec.heads[1].damping_ratio=0.0;spec.depth_m=0.20;
    let file=text(spec);
    assert_eq!(Spec::read(&file).unwrap(),spec);
    let mut records:Vec<_>=file.lines().skip(1).collect();records.reverse();
    assert_eq!(Spec::read(&format!("# inputs, not certified measurements\n{HEADER}\n{}\n",records.join("\n"))).unwrap(),spec);
    let film=spec.head(1).unwrap();
    assert_eq!(film.spec.young_pa,spec.heads[1].young_pa);
    assert_eq!(film.spec.thickness_m,spec.heads[1].thickness_m);
    assert_eq!(film.spec.density_kg_m3,spec.heads[1].density_kg_m3);
    assert_eq!(film.spec.poisson,spec.heads[1].poisson);
    assert_eq!(film.spec.tension_n_m,spec.heads[1].tension_n_m);
    assert_eq!(film.spec.radius_m,spec.radius_m);
    assert_eq!(film.mesh.nodes.len(),1+spec.radial_intervals*spec.azimuths);
}

#[test]
fn incomplete_duplicate_nonfinite_and_nonphysical_inputs_refuse() {
    let valid=text(small());
    for (i,_) in valid.lines().enumerate() {
        let missing=valid.lines().enumerate().filter(|(j,_)|*j!=i)
            .map(|(_,s)|s).collect::<Vec<_>>().join("\n");
        assert!(Spec::read(&missing).is_err(),"accepted missing record {i}");
    }
    for row in valid.lines().skip(1) {
        assert!(Spec::read(&format!("{valid}{row}\n")).is_err());
    }
    for bad in [String::new(),valid.replace(HEADER,"wrong-version"),
        format!("{valid}pitch,440\n"),format!("{valid}{}","#".repeat(MAX_BYTES)),
        valid.replace("mesh,2,16","mesh,2.5,16"),valid.replace("mesh,2,16","mesh,0,16"),
        valid.replace("mesh,2,16","mesh,33,16"),valid.replace("mesh,2,16","mesh,2,129"),
        valid.replace("head,batter,","head,unknown,"),
        valid.replacen("geometry,","geometry,NaN,",1)] {
        assert!(Spec::read(&bad).is_err(),"accepted malformed declaration");
    }
    for bad in [f64::NAN,f64::INFINITY,-1.0,0.0] {
        let mut s=small();s.heads[0].tension_n_m=bad;
        assert!(Spec::read(&text(s)).is_err());
        let mut s=small();s.radius_m=bad;
        assert!(Spec::read(&text(s)).is_err());
    }
    let mut s=small();s.outer_radius_m=s.radius_m;assert!(s.validate().is_err());
    let mut s=small();s.heads[0].damping_ratio=-0.01;assert!(s.validate().is_err());
    let mut s=small();s.heads[1].young_pa=0.0;assert!(s.validate().is_err());
    let mut s=small();s.heads[1].poisson=1.0;assert!(s.validate().is_err());
}

#[test]
fn option_composes_with_playing_and_existing_physical_images() {
    let mut args:Vec<_>=["drum-stretch-mic","128","--drum-spec","instrument.fsd",
        "--cavity-modes","--prepared-nonlinear","--strike-speed-m-s","4"]
        .into_iter().map(String::from).collect();
    assert_eq!(option(&mut args).unwrap().as_deref(),Some("instrument.fsd"));
    assert!(crate::cavity::option(&mut args).unwrap());
    assert!(crate::mechanics::prepared_option(&mut args).unwrap());
    let (args,stroke)=crate::playing::parse(args).unwrap();
    assert_eq!(args,["drum-stretch-mic","128"]);assert_eq!(stroke.speed_m_s,4.0);
    for command in ["drum","drum-mic","drum-stretch-wav","drum-modal","snare","snare-off-mic"] {
        assert!(admit_command(Some("instrument.fsd"),command).is_ok());
    }
    for command in ["splash","splash-mic","unknown"] {
        assert!(admit_command(Some("instrument.fsd"),command).is_err());
    }
    for raw in ["--drum-spec", "--drum-spec --cavity-modes", "--drum-spec a --drum-spec b"] {
        let mut args:Vec<_>=raw.split_whitespace().map(String::from).collect();
        let saved=args.clone();assert!(option(&mut args).is_err());assert_eq!(args,saved);
    }
}

#[test]
fn clocks_acoustic_band_and_fixed_snare_geometry_are_not_silently_changed() {
    let spec=small();
    for dt in [0.0,-1.0,f64::NAN,0.001] {assert!(spec.admit_clock(dt,false).is_err());}
    let wide=Spec {band_hz:[20.0,2000.0],..spec};
    assert!(wide.admit_clock(2e-6,false).is_ok());
    assert!(wide.admit_clock(crate::acoustics::MECHANICAL_DT,true).is_err());
    let wires=crate::snare::SnareSet::reference(false);
    assert!(spec.admit_snare(wires).is_ok());
    assert!(Spec {radius_m:0.08,outer_radius_m:0.09,..spec}.admit_snare(wires).is_err());
}

#[test]
fn supplied_reference_keeps_the_original_contact_air_and_accepted_trajectory() {
    let supplied=Spec::read(&text(Spec::reference())).unwrap();
    let mut original=drum_with_air(64,2e-6,false,false,None,false,stroke(),false,None).unwrap();
    let mut imported=drum_with_spec(64,2e-6,false,false,None,false,stroke(),false,None,Some(supplied)).unwrap();
    assert_eq!(original.observer_a,imported.observer_a);
    assert_eq!(original.observer_b,imported.observer_b);
    assert_eq!(original.pressure.as_ref().unwrap().volume_m3,imported.pressure.as_ref().unwrap().volume_m3);
    assert_eq!(original.pressure.as_ref().unwrap().areas,imported.pressure.as_ref().unwrap().areas);
    assert_eq!(original.system.state(),imported.system.state());
    let gate=CancelGate::new_clock_free();
    for _ in 0..64 {
        let a=original.system.step(&original.force,&gate).unwrap();
        let b=imported.system.step(&imported.force,&gate).unwrap();
        assert_eq!(a.stored_energy_j.to_bits(),b.stored_energy_j.to_bits());
        assert_eq!(original.system.state(),imported.system.state());
    }
}

#[test]
fn tension_and_elasticity_change_the_actual_pencil_and_resonances() {
    let a=small();let mut b=a;
    for h in &mut b.heads {h.young_pa*=4.0;h.tension_n_m*=4.0;}
    b.band_hz=a.band_hz.map(|f|2.0*f);
    let (fa,ma)=a.prepare(2e-6,false).unwrap();
    let (fb,mb)=b.prepare(2e-6,false).unwrap();
    for i in 0..2 {
        assert_eq!(fa[i].mass_kg,fb[i].mass_kg);
        assert_eq!(fa[i].mesh.nodes,fb[i].mesh.nodes);
        assert_eq!(ma[i].len(),mb[i].len());
        let n=fa[i].model.free;
        for r in 0..n {for c in 0..n {
            assert_eq!(fa[i].model.m.get(r,c),fb[i].model.m.get(r,c));
            let expected=4.0*fa[i].model.k.get(r,c);
            assert!((fb[i].model.k.get(r,c)-expected).abs()<1e-10*(1.0+expected.abs()));
        }}
        for (x,y) in ma[i].iter().zip(&mb[i]) {
            assert!((y.lambda/(4.0*x.lambda)-1.0).abs()<1e-5,"physical stiffness scaling must double frequency");
        }
    }
    let mut tuned=a;for h in &mut tuned.heads {h.tension_n_m*=2.0;}
    let (_,mt)=tuned.prepare(2e-6,false).unwrap();
    let lowest=|m:&[ModePair]|m.iter().map(|p|p.lambda).fold(f64::INFINITY,f64::min);
    for i in 0..2 {assert!(lowest(&mt[i])>lowest(&ma[i]));}
}

#[test]
fn depth_changes_real_air_feedback_without_retuning_the_heads() {
    let a=small();let b=Spec {depth_m:1.5*a.depth_m,..a};
    let build=|s|drum_with_spec(128,2e-6,false,false,None,false,stroke(),false,None,Some(s)).unwrap();
    let mut shallow=build(a);let mut deep=build(b);
    assert_eq!(shallow.observer_a,deep.observer_a);
    assert_eq!(shallow.system.state(),deep.system.state());
    let va=shallow.pressure.as_ref().unwrap();let vb=deep.pressure.as_ref().unwrap();
    assert_eq!(va.areas,vb.areas);
    assert!((vb.volume_m3/va.volume_m3-1.5).abs()<1e-14);
    let gate=CancelGate::new_clock_free();let mut changed=0.0_f64;
    for _ in 0..128 {
        for e in [&mut shallow,&mut deep] {
            let f=e.system.step(&e.force,&gate).unwrap();
            assert!(f.stored_energy_j.is_finite());assert!(f.balance_residual_j.abs()<2e-7);
        }
        for (x,y) in shallow.system.state().iter().zip(deep.system.state()) {changed=changed.max((x-y).abs());}
    }
    assert!(changed>1e-14,"changed cavity volume must feed back into the accepted head/contact motion");
}

#[test]
fn larger_supplied_heads_reach_their_matching_exterior_not_the_stock_shell() {
    let s=Spec {radius_m:0.20,outer_radius_m:0.21,depth_m:0.14,..small()};
    let e=drum_with_spec(2,crate::acoustics::MECHANICAL_DT,true,false,None,false,
        stroke(),false,None,Some(s)).unwrap();
    assert!(e.acoustics.is_some());
    assert_eq!(e.pressure.as_ref().unwrap().volume_m3,s.volume_m3());
    // The old 0.1778m outer surface cannot enclose these 0.20m heads.
    let (films,modes)=s.prepare(2e-6,false).unwrap();
    assert!(crate::acoustics::Boundary::drum(&films,&modes,s.depth_m,0.1778).is_err());
}

#[test]
fn supplied_drum_retains_stretching_distributed_air_and_prepared_nonlinearity() {
    let s=Spec {depth_m:0.18,..small()};
    let mut e=drum_with_spec(64,2e-6,false,false,None,true,stroke(),true,None,Some(s)).unwrap();
    assert!(e.air.is_some());
    e.system=e.system.into_prepared_nonlinear().unwrap();
    assert!(matches!(&e.system,Mechanics::Nonlinear(_)));
    let gate=CancelGate::new_clock_free();let mut stretching=0.0_f64;
    for _ in 0..64 {
        let f=e.system.step(&e.force,&gate).unwrap();assert!(f.balance_residual_j.abs()<2e-7);
        stretching=stretching.max(e.system.membrane_observation(1).unwrap().stretching_energy_j);
        let (a,b)=e.air.as_ref().unwrap().points(e.system.state()).unwrap();
        assert!(a.is_finite() && b.is_finite());
    }
    assert!(stretching>0.0);
}
