use super::*;
use fs_couple::render::plate::impact::{ImpactBody,BodyPotential,ImpactSystem,radiation::Pole};

// Manufactured breathing coordinate on a real triangulated closed BEM surface.
// This is an integration oracle, NOT a cymbal eigenmode or a calibrated body.
fn boundary()->Boundary {
    let s=SpherePanels::icosphere(0.02,0).unwrap();
    Boundary{triangles:s.triangles().unwrap().to_vec(),weights:vec![vec![100.;s.areas().len()]],state_modes:vec![0]}
}
fn spec()->radiation_spec::Spec {
    radiation_spec::Spec{band_hz:[40.,400.],training_intervals:16,max_order:8,subdivisions:0,
        max_panels:80,max_dense_work:100_000_000}
}
fn experiment(frames:u64)->Experiment {
    let mut body=ImpactBody::free_mass(1.,0.,0.0001).unwrap().0;
    body.potential=BodyPotential::Linear(vec![1000.]);
    let s=ImpactSystem::new(vec![body],vec![],vec![],vec![],crate::config(frames*SUBSTEPS as u64,MECHANICAL_DT)).unwrap();
    Experiment{flexible_sticks:[None,None],mute:None,system:Mechanics::Reference(s),force:vec![0.],stick_weight:1.,second_stick:None,
        observer_a:vec![100.],observer_b:vec![0.],pressure:None,acoustics:Some(boundary()),air:None}
}
fn energy(s:&Mechanics)->f64 {
    match s{Mechanics::Reference(s)=>s.stored_energy_j(),Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(),Mechanics::Driven{inner,..}=>energy(inner),_=>panic!()}
}
#[test]
fn acceleration_bem_projection_matches_independent_velocity_driven_boundary_work() {
    let mut b=boundary();b.weights.push(b.weights[0].iter().map(|x|-0.4*x).collect());
    let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let w=500.;let medium=Medium::air();
    let a=acceleration_fields(&b.weights,w);let ar:Vec<_>=a.iter().map(Vec::as_slice).collect();
    let acc=solve_radiation_batch(&surface,w/medium.sound_speed,medium,&ar,Formulation::PlainCbie).unwrap();
    let actual=project(&surface,&b.weights,&acc,w).unwrap();
    let v:Vec<Vec<_>>=b.weights.iter().map(|row|row.iter().map(|x|C64::from_re(*x)).collect()).collect();
    let vr:Vec<_>=v.iter().map(Vec::as_slice).collect();
    let vel=solve_radiation_batch(&surface,w/medium.sound_speed,medium,&vr,Formulation::PlainCbie).unwrap();
    for i in 0..2 {for j in 0..2 {
        let expected=vel[j].pressure.iter().zip(surface.areas()).zip(&b.weights[i])
            .fold(C64::ZERO,|a,((p,area),b)|a+p.scale(area*b));
        assert!((actual[2*i+j]-expected).abs()<1e-11*expected.abs());
    }}
    assert!(actual[0].re>0. && actual[0].im<0.);
    assert!(actual[1].re<0.,"signed off-diagonal radiation is not a negative physical loss");
    assert!((actual[1]-actual[2]).abs()<1e-11*actual[0].abs());
}
#[test]
fn full_power_form_rejects_active_cross_coupling_even_with_positive_diagonals() {
    let r=C64::from_re;
    for z in [vec![r(1.),r(2.),r(2.),r(1.)],
        vec![r(1.),C64::new(0.,2.),C64::new(0.,-2.),r(1.)]] {
        assert!(admit_power(&z,2).is_err());
    }
    assert!(admit_power(&[r(1.),r(-0.5),r(-0.5),r(1.)],2).is_ok());
    assert!(admit_power(&[C64::new(f64::NAN,0.)],1).is_err());
    assert!(admit_power(&[r(1.)],2).is_err());
}
#[test]
fn real_bem_load_changes_played_motion_and_both_receivers_share_its_history() {
    let gate=CancelGate::new_clock_free();let receivers=[Receiver::FinitePoint([0.,0.,0.12]),Receiver::FinitePoint([0.1,0.,0.08])];
    let original=experiment(96);let prefix=original.system.state().to_vec();
    let (mut e,prepared)=prepare(original,96,20.,&receivers,spec(),&gate).unwrap();
    assert_eq!(&e.system.state()[..prefix.len()],prefix);assert_eq!(e.force,vec![0.]);
    assert!(e.system.state()[prefix.len()..].iter().all(|v|*v==0.));
    assert!(observation(&e.system).is_some());let initial=energy(&e.system);
    // A manual one-way trajectory uses the IDENTICAL admitted observer bank.
    let mut unreacted=experiment(96);unreacted.system=unreacted.system.into_analytic_nonlinear().unwrap();
    let one_way=render_baked_with_gate(&mut unreacted,96,20.,&prepared.baked,&gate).unwrap();
    e.system=e.system.into_analytic_nonlinear().unwrap();
    let wav=prepared.render(&mut e,96,20.,&gate).unwrap();
    assert_eq!(&wav[..4],b"RIFF");assert_eq!(wav.len(),44+96*2*2);assert_ne!(wav,one_way);
    assert_ne!(&e.system.state()[..prefix.len()],unreacted.system.state());
    assert!(observation(&e.system).unwrap().stored_energy_j>0.);
    assert!(energy(&e.system)<initial);
}
#[test]
fn acoustic_reaction_composes_with_real_stick_head_cavity_and_material_memory() {
    let mut drum=crate::drum_spec::Spec{radial_intervals:2,azimuths:8,..crate::drum_spec::Spec::reference()};
    for h in &mut drum.heads{h.damping_ratio=0.;}
    let material=crate::head_relaxation::Spec::read(include_str!("estimated-head-relaxation.fshr")).unwrap();
    let make=||crate::drum_with_material(256,2e-6,true,false,None,true,
        crate::Stroke{speed_m_s:2.,position_m:Some([0.06,0.01])},true,None,Some(drum.clone()),None,&[],20.,None,false,Some(&material)).unwrap();
    let mut original=make();let mut loaded=make();let prefix=original.system.state().len();
    let sources=loaded.acoustics.as_ref().unwrap().state_modes.clone();
    // Declared test load, not a BEM-derived specimen claim; the previous test
    // separately executes BEM fitting and loaded pressure rendering.
    let model=Model{ports:sources.len(),poles:vec![Pole{omega:1200.,zeta:0.2,
        coupling:sources.iter().enumerate().map(|(i,_)|if i%2==0{150.}else{-90.}).collect()}]};
    let Mechanics::Reference(s)=loaded.system else{panic!()};
    loaded.system=Mechanics::Reference(s.with_radiation_load(&model,&sources,1).unwrap());
    assert_eq!(&loaded.system.state()[..prefix],original.system.state());assert_eq!(loaded.force,original.force);
    assert_eq!(loaded.observer_a,original.observer_a);assert_eq!(loaded.observer_b,original.observer_b);
    assert_eq!(loaded.acoustics.as_ref().unwrap().state_modes,original.acoustics.as_ref().unwrap().state_modes);
    let initial=energy(&loaded.system);let mut loss=0.;let mut changed=0.0_f64;
    loaded.system=loaded.system.into_analytic_nonlinear().unwrap();original.system=original.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..256 {
        let f=loaded.system.step(&loaded.force,&gate).unwrap();original.system.step(&original.force,&gate).unwrap();
        loss+=f.dissipated_energy_j;assert!(f.balance_residual_j.abs()<1e-7);
        for (a,b) in loaded.system.state().iter().zip(original.system.state()){changed=changed.max((a-b).abs());}
    }
    assert!(changed>0. && observation(&loaded.system).unwrap().stored_energy_j>0.);
    assert!(crate::head_relaxation::observation(&loaded.system).stored_energy_j>0.);
    assert!((energy(&loaded.system)+loss-initial).abs()<1e-6);
}
#[test]
fn feedback_admission_never_hides_a_different_model_or_consumes_playing_time() {
    let mut args=vec!["hihat-mic".into(),"pair.fshh".into(),"--radiation-feedback".into()];
    assert!(option(&mut args).unwrap());assert_eq!(args,["hihat-mic","pair.fshh"]);
    let mut dup=vec!["--radiation-feedback".into();2];let before=dup.clone();assert!(option(&mut dup).is_err());assert_eq!(dup,before);
    for name in ["splash-mic","drum-stretch-wav","snare-mic","hihat-wav"]{admit_command(true,name,false).unwrap();}
    for name in ["splash","drum","drum-modal-mic","hihat"]{assert!(admit_command(true,name,false).is_err());}
    assert!(admit_command(true,"drum-mic",true).is_err());
    let e=experiment(4);let before=e.system.state().to_vec();let gate=CancelGate::new_clock_free();gate.request();
    assert!(bake_scene_with_spec(e.acoustics.as_ref().unwrap(),&[Receiver::FinitePoint([0.,0.,0.12])],spec(),&gate,true).is_err());
    assert_eq!(e.system.state(),before);
    let mut b=boundary();b.weights=vec![b.weights[0].clone();33];b.state_modes=(0..33).collect();assert!(admit_boundary(&b,spec()).is_err());
    assert!(admit_boundary(&boundary(),radiation_spec::Spec{training_intervals:65,..spec()}).is_err());
    let mut e=experiment(4);e.system=e.system.into_analytic_nonlinear().unwrap();assert!(admit_instrument(&e).is_err());
    let e=crate::drum(64,MECHANICAL_DT,true,true).unwrap();assert!(admit_instrument(&e).is_err());
}
