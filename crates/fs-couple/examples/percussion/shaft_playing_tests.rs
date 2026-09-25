use super::*;
use crate::{compliant_mute, drum_spec, drum_with_mallets, drum_with_shafts,
    mechanics::{self, Mechanics}, snare, specimen, splash_with_mallets, splash_with_shafts};
use fs_couple::render::plate::impact::{ImpactSubstepConfig, supported::TranslatingSupport};
use fs_couple::render::plate::impact::striker::RadiusStation;
use fs_exec::CancelGate;
use fs_plate::shell::stiffened::beam::RoundBeamSpec;

fn shaft() -> FlexibleStriker {
    FlexibleStriker::new(&[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)|
        RadiusStation{position_m,radius_m}),RoundBeamSpec{young_pa:12e9,density_kg_m3:800.0,
        pivot_m:0.1,contact_m:0.39,hand_m:0.16,subdivisions:8,maximum_hz:3000.0,maximum_modes:17},0.001).unwrap()
}
fn selection() -> Selection {Selection{first:Some(shaft()),second:Some(shaft())}}
fn first() -> Stroke {Stroke{speed_m_s:0.8,position_m:Some([0.06,0.01])}}
fn second() -> Stroke {Stroke{speed_m_s:0.6,position_m:Some([-0.05,0.02])}}
fn drum() -> drum_spec::Spec {drum_spec::Spec{radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}}
fn shell() -> specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn drum_build(s:&Selection,audio:bool,wires:Option<snare::SnareSet>,air:bool,stretch:bool) -> Experiment {
    drum_with_shafts(256,2e-6,audio,false,wires,stretch,first(),air,None,Some(drum()),Some(second()),
        &[],if air{25.0}else{0.0},None,false,None,&mallets::Selection::default(),s).unwrap()
}
fn shell_build(s:&Selection,audio:bool) -> Experiment {
    splash_with_shafts(256,2e-6,audio,first(),Some(shell()),&[],Some(second()),None,
        &mallets::Selection::default(),s).unwrap()
}
fn energy(s:&Mechanics) -> f64 {
    match s {Mechanics::Reference(s)=>s.stored_energy_j(),Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(),Mechanics::Driven{inner,..}=>energy(inner),
        Mechanics::Prepared(_)=>panic!("nonlinear-capable fixture")}
}
fn prepared(mut e:Experiment) -> Experiment {
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();e
}

#[test]
fn selection_refuses_missing_hands_overlapping_tip_models_and_clock_aliasing() {
    let both=selection();let no_mallet=mallets::Selection::default();
    for command in ["splash","drum","drum-stretch-mic","snare-off-wav","drum-modal"] {
        both.admit(command,Some(second()),&no_mallet).unwrap();
    }
    assert!(both.admit("unknown",Some(second()),&no_mallet).is_err());
    assert!(both.admit("drum",None,&no_mallet).is_err());
    let felt=mallets::Spec::parse(include_str!("estimated-felt-mallet.fsmallet")).unwrap();
    let first_mallet=mallets::Selection{first:Some(felt),second:None};
    assert!(both.admit("splash",Some(second()),&first_mallet).is_err());
    Selection{first:None,second:Some(shaft())}.admit("drum",Some(second()),&first_mallet).unwrap();
    assert!(both.build(4,2,first(),Some(second()),0.01).is_err());
    assert!(both.build(usize::MAX,2,first(),Some(second()),2e-6).is_err());
    assert!(both.build(4,4,first(),Some(second()),2e-6).is_err(),"rigid coordinates cannot overlap appended flexure");
}

#[test]
fn unselected_shafts_preserve_original_trajectories_and_output_columns() {
    let mut old=drum_with_mallets(256,2e-6,false,false,None,false,first(),false,None,Some(drum()),
        Some(second()),&[],0.0,None,false,None,&mallets::Selection::default()).unwrap();
    let mut new=drum_build(&Selection::default(),false,None,false,false);
    let gate=CancelGate::new_clock_free();
    for _ in 0..8 {let a=old.system.step(&old.force,&gate).unwrap();let b=new.system.step(&new.force,&gate).unwrap();
        assert_eq!(old.system.state(),new.system.state());assert_eq!(a.stored_energy_j,b.stored_energy_j);}
    let old=splash_with_mallets(256,2e-6,false,first(),Some(shell()),&[],Some(second()),None,&mallets::Selection::default()).unwrap();
    let new=shell_build(&Selection::default(),false);
    assert_eq!(old.system.state(),new.system.state());assert_eq!(old.observer_b,new.observer_b);
    let mut h=Vec::new();let mut r=Vec::new();header(&new,&mut h).unwrap();row(&new,&mut r).unwrap();
    assert!(h.is_empty()&&r.is_empty());
}

#[test]
fn both_resonators_keep_the_same_acoustic_sources_and_rigid_only_launch() {
    for e in [shell_build(&selection(),true),drum_build(&selection(),true,None,false,true)] {
        let sources=e.acoustics.as_ref().unwrap().state_modes();
        let p=e.flexible_sticks[0].as_ref().unwrap();let q=e.flexible_sticks[1].as_ref().unwrap();
        assert_eq!(p.rigid_coordinate(),0);assert_eq!(q.rigid_coordinate(),e.second_stick.unwrap().coordinate);
        assert_eq!(p.elastic_start(),q.rigid_coordinate()+1);
        assert_eq!(q.elastic_start(),p.elastic_start()+p.elastic_modes());
        assert_eq!(e.force.len(),q.elastic_start()+q.elastic_modes());
        for (hand,port) in [p,q].into_iter().enumerate() {
            assert!(!sources.contains(&port.rigid_coordinate()));
            for k in port.elastic_start()..port.elastic_start()+port.elastic_modes() {
                assert!(!sources.contains(&k));assert_eq!(e.observer_a[k],0.0);
                if let Some(v)=&e.pressure {assert_eq!(v.areas[k],0.0);}
            }
            let o=port.observe(e.system.state()).unwrap();assert_eq!(o.flexural_energy_j,0.0);
            assert!((o.tip_displacement_m+0.0002).abs()<1e-16);
            let speed=if hand==0{first().speed_m_s}else{second().speed_m_s};
            assert!((o.tip_velocity_m_s-speed).abs()<1e-14);
            assert_ne!(port.tip_row(e.force.len()).unwrap(),port.hand_row(e.force.len()).unwrap());
        }
        let mut h=Vec::new();let mut r=Vec::new();header(&e,&mut h).unwrap();row(&e,&mut r).unwrap();
        assert_eq!(String::from_utf8(h).unwrap().split(',').count(),11);
        assert_eq!(String::from_utf8(r).unwrap().split(',').count(),11);
    }
}

#[test]
fn actual_drum_and_shell_contacts_excite_both_shafts_without_invented_work() {
    for mut e in [prepared(shell_build(&selection(),false)),prepared(drum_build(&selection(),false,None,false,true))] {
        let initial=energy(&e.system);let gate=CancelGate::new_clock_free();
        let (mut loss,mut bending)=(0.0,[0.0_f64;2]);
        for tick in 0..256 {
            if tick==160 {let old=e.system.state().to_vec();let mut bad=e.force.clone();bad[0]=f64::NAN;
                assert!(e.system.step(&bad,&gate).is_err());assert_eq!(e.system.state(),old);}
            let f=e.system.step(&e.force,&gate).unwrap();loss+=f.dissipated_energy_j;
            assert_eq!(f.supplied_work_j,0.0);assert!(f.balance_residual_j.abs()<1e-7);
            for (i,p) in e.flexible_sticks.iter().enumerate() {
                bending[i]=bending[i].max(p.as_ref().unwrap().observe(e.system.state()).unwrap().flexural_energy_j);
            }
        }
        assert!(bending.iter().all(|e|*e>1e-12),"both real impacts must excite bending: {bending:?}");
        assert!((energy(&e.system)+loss-initial).abs()<1e-6);
    }
}

#[test]
fn physical_hand_force_matches_manual_full_row_projection_and_retry() {
    let s=Selection{first:Some(shaft()),second:None};
    let mut manual=prepared(drum_build(&s,false,None,false,false));
    let mut driven=prepared(drum_build(&s,false,None,false,false));
    let source="0,0\n0.00008,0.1\n0.00016,0\n0.00032,-0.05\n0.0005,0";
    let programs=||[Some(Program::parse(source).unwrap()),None];
    let (scalar,spatial)=player_inputs(&driven,programs()).unwrap();assert!(scalar.is_empty());
    assert_eq!(spatial.len(),1);
    let expected=manual.flexible_sticks[0].as_ref().unwrap().hand_row(manual.force.len()).unwrap();
    assert_eq!(spatial[0].weights,expected);
    driven.system=driven.system.with_player_drives(scalar,spatial,2e-6,256,driven.force.len()).unwrap();
    let (scalar,spatial)=player_inputs(&manual,programs()).unwrap();
    let mut force=mechanics::drive::StickDrive::new_mixed(scalar,spatial,2e-6,256,manual.force.len()).unwrap();
    let gate=CancelGate::new_clock_free();let initial=energy(&driven.system);let mut net=0.0;
    for tick in 0..256 {
        if tick==140 {let before=driven.system.state().to_vec();let mut bad=driven.force.clone();bad[0]=1e7;
            assert!(driven.system.step(&bad,&gate).is_err());assert_eq!(driven.system.state(),before);}
        let a=manual.system.step(force.forces(&manual.force).unwrap(),&gate).unwrap();force.accept();
        let b=driven.system.step(&driven.force,&gate).unwrap();
        assert_eq!(manual.system.state(),driven.system.state());assert_eq!(a.supplied_work_j,b.supplied_work_j);
        net+=b.supplied_work_j-b.dissipated_energy_j;
    }
    assert!((energy(&driven.system)-initial-net).abs()<1e-6);
}

#[test]
fn full_snare_carrier_and_cavity_precede_no_hidden_direct_shaft_pressure() {
    let source=snare::SnareSet::reference(false);
    let carrier=TranslatingSupport{mass_kg:0.02,stiffness_n_m:12.0,damping_n_s_m:0.02,
        initial_position_m:0.0,initial_velocity_m_s:0.0,maximum_travel_m:0.01,maximum_slope:0.2};
    let wires=snare::SnareSet{carrier:Some(carrier),..source};
    let e=drum_build(&selection(),true,Some(wires),true,true);
    let p=e.flexible_sticks[0].as_ref().unwrap();let q=e.flexible_sticks[1].as_ref().unwrap();
    let start=e.second_stick.unwrap().coordinate+1;
    assert_eq!(p.elastic_start(),start+wires.mode_count().unwrap());
    let (index,_)=snare::carrier::force_port(&e.system,4).unwrap();assert_eq!(index,start+160);
    let air=e.air.as_ref().unwrap();assert_eq!(air.coupling.structural_modes(),q.elastic_start()+q.elastic_modes());
    assert!(e.acoustics.as_ref().unwrap().state_modes().iter().all(|k|*k<start-1));
    let mut displaced=e.system.state().to_vec();
    let before=air.points(&displaced).unwrap();
    for port in [p,q] {for k in port.elastic_start()..port.elastic_start()+port.elastic_modes(){displaced[2*k]=0.0001;}}
    assert_eq!(air.points(&displaced).unwrap(),before);assert_eq!(air.uniform_pressure(&displaced).unwrap(),0.0);
}

#[test]
fn mute_jaws_and_opposite_hand_mallet_keep_their_own_inertia_and_history() {
    let s=Selection{first:None,second:Some(shaft())};
    let selected=mallets::Selection{first:Some(mallets::Spec::parse(include_str!("estimated-felt-mallet.fsmallet")).unwrap()),second:None};
    let mute=compliant_mute::Spec::parse("frankensim-compliant-mute-v1\nsurface,shell\n\
        site,0.073,-0.005,0.0001\nsite,0.083,0.005,0.0001\n\
        jaw,above,0.025,0.15,0.00002,0.004,30000,0.2,2.2,3,0.15,0.7,0\n\
        creep,above,4000,12\nforce,above,0,0\nforce,above,0.0005,0").unwrap();
    let e=splash_with_shafts(256,2e-6,true,first(),Some(shell()),&[],Some(second()),Some(&mute),&selected,&s).unwrap();
    let p=e.flexible_sticks[1].as_ref().unwrap();let jaws=e.mute.as_ref().unwrap();
    assert!(e.flexible_sticks[0].is_none());assert_eq!(p.elastic_start(),jaws.ports[0].coordinate+1);
    assert_eq!(jaws.first_pad,10,"six stand plus four independent mallet patches");
    let sources=e.acoustics.as_ref().unwrap().state_modes();
    assert!(sources.iter().all(|k|*k<e.second_stick.unwrap().coordinate));
    // All finite-area felt histories remain private trailing states. Bending
    // cannot be mistaken for a stand, mallet or moving-jaw creep coordinate.
    let Mechanics::Reference(system)=&e.system else {panic!()};
    for i in 0..12 {assert!(system.felt_history(i).is_some());}
    assert_eq!(p.observe(e.system.state()).unwrap().flexural_energy_j,0.0);
}

#[test]
fn linear_modal_image_retains_the_selected_elastic_shaft_bodies() {
    let s=selection();
    let e=drum_with_shafts(4,2e-6,true,true,None,false,first(),false,None,Some(drum()),
        Some(second()),&[],0.0,None,false,None,&mallets::Selection::default(),&s).unwrap();
    let Mechanics::Prepared(system)=&e.system else {panic!("explicit modal image must remain modal")};
    assert_eq!(system.contact_count(),2);
    assert_eq!(system.state().len(),2*e.force.len());
    for p in e.flexible_sticks.iter().flatten() {
        assert!(p.elastic_modes()>0);assert_eq!(p.observe(system.state()).unwrap().flexural_energy_j,0.0);
        assert!(e.acoustics.as_ref().unwrap().state_modes().iter().all(|i|*i<p.elastic_start()));
    }
}
