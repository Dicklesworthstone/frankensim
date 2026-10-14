//! Existing instrument geometry, contacts and material memory with joint kinematic/air elimination.
use crate::*;
use fs_couple::render::plate::impact::{PreparedImpactSystem,radiation::{Model,Pole}};

struct Run {system:PreparedImpactSystem,force:Vec<f64>,sources:Vec<usize>,observers:[Vec<f64>;2]}
fn attach(e:Experiment,condensed:bool)->Run {
    let sources=e.acoustics.as_ref().unwrap().state_modes().to_vec();
    let model=Model{ports:sources.len(),poles:(0..16).map(|i|Pole {
        omega:600.+40.*i as f64,zeta:0.2,
        coupling:sources.iter().enumerate().map(|(j,_)|if (i+j)%2==0{0.8}else{-0.4}).collect(),
    }).collect()};
    // This manufactured passive load isolates execution equivalence. Existing
    // radiation-feedback tests separately execute actual BEM/load/receiver fits.
    let Mechanics::Reference(s)=e.system else{panic!("unprepared physical constructor")};
    let prefix=s.state().to_vec();let before=s.stored_energy_j();
    let mut system=s.with_radiation_load(&model,&sources,16).unwrap().prepare_analytic().unwrap();
    assert_eq!(&system.state()[..prefix.len()],prefix);assert_eq!(system.stored_energy_j(),before);
    assert_eq!(system.condensed_newton_dimension(),Some(prefix.len()-e.force.len()+1));
    system.set_radiation_condensation(condensed).unwrap();
    Run{system,force:e.force,sources,observers:[e.observer_a,e.observer_b]}
}
fn compare(mut fast:Run,mut dense:Run,ticks:usize)->Run {
    assert_eq!(fast.force,dense.force);assert_eq!(fast.sources,dense.sources);assert_eq!(fast.observers,dense.observers);
    assert_eq!(fast.system.state(),dense.system.state());
    let stopped=CancelGate::new_clock_free();stopped.request();let before=fast.system.state().to_vec();
    assert!(fast.system.step(&fast.force,&stopped).is_err());assert_eq!(fast.system.state(),before);
    let gate=CancelGate::new_clock_free();let initial=fast.system.stored_energy_j();
    let (mut work,mut loss,mut condensed_solves,mut full_solves,mut motion)=(0.,0.,0,0,0.0_f64);
    for tick in 0..ticks {
        fast.force[0]=if tick<40 {0.04}else{0.};dense.force.copy_from_slice(&fast.force);
        let f=fast.system.step(&fast.force,&gate).unwrap();let d=dense.system.step(&dense.force,&gate).unwrap();
        assert_eq!(f.time_s,d.time_s);assert_eq!(f.sample,d.sample);
        for (a,b) in fast.system.state().iter().zip(dense.system.state()){assert!((a-b).abs()<2e-8);}
        assert!(f.balance_residual_j.abs()<1e-7);work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
        // The owner explicitly retains full LU when elimination is numerically
        // unsuitable. That safe path must obey the SAME state/work comparison,
        // not fail this physical regression solely for taking one fallback.
        let calls=fast.system.newton_linear_solve_counts();condensed_solves+=calls.0;full_solves+=calls.1;
        assert_eq!(dense.system.newton_linear_solve_counts().0,0);
        let displacement=fast.observers[0].iter().enumerate().map(|(i,b)|b*fast.system.state()[2*i]).sum::<f64>();
        motion=motion.max(displacement.abs());
    }
    assert!(condensed_solves>full_solves && motion>1e-10,"condensation must still do the majority of this fixture's solves");
    assert!((fast.system.stored_energy_j()-initial+loss-work).abs()<1e-6);
    assert!(fast.system.radiation_observation().unwrap().stored_energy_j>0.);
    fast
}

#[test]
fn condensed_radiation_keeps_two_stick_curved_shell_motion_and_stand_history() {
    let make=||{
        let mut spec=specimen::Specimen::reference();spec.azimuths=8;
        splash_with_sticks(256,2e-6,true,Stroke{speed_m_s:0.8,position_m:Some([0.06,0.01])},
            Some(spec),&[],Some(Stroke{speed_m_s:0.6,position_m:Some([-0.05,0.02])})).unwrap()
    };
    let fast=compare(attach(make(),true),attach(make(),false),220);
    for i in 0..6 {assert!(fast.system.felt_history(i).is_some());}
    assert!(fast.sources.iter().all(|&i|i>0 && i<fast.force.len()-1));
}

#[test]
fn condensed_radiation_keeps_stretching_heads_cavity_and_hereditary_material() {
    let mut spec=drum_spec::Spec{radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()};
    for head in &mut spec.heads {head.damping_ratio=0.;}
    let material=head_relaxation::Spec::read(include_str!("estimated-head-relaxation.fshr")).unwrap();
    let make=||drum_with_material(256,2e-6,true,false,None,true,
        Stroke{speed_m_s:2.,position_m:Some([0.06,0.01])},true,None,Some(spec),None,&[],20.,None,false,Some(&material)).unwrap();
    let fast=compare(attach(make(),true),attach(make(),false),180);
    assert!(fast.system.relaxation_observation().unwrap().stored_energy_j>0.);
    assert!(fast.system.membrane_observation(1).unwrap().stretching_energy_j>0.);
    assert!(fast.sources.len()<fast.force.len()-1,"cavity coordinates must not become radiation sources");
}


#[test]
fn kinematic_radiation_solve_keeps_active_wire_contacts_nonlinear_heads_and_cavity() {
    let wires=snare::SnareSet {strands:2,modes_per_strand:2,contact_cells:4,
        clearance_m:-2e-6,..snare::SnareSet::reference(false)};
    let make=||drum_with_cavity_loss(256,2e-6,true,false,Some(wires),true,
        Stroke {speed_m_s:0.0,position_m:Some([0.06,0.01])},true,None,
        Some(drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}),None,
        &[muffling::Muffler {surface:muffling::Surface::Batter,position_m:[0.07,0.01],resistance_n_s_m:0.1}],25.0).unwrap();
    // Interference is deliberate initial stored contact energy, NOT an inferred
    // equilibrium or a calibrated snare setting. Both numerical paths see it.
    let first=make();let start=first.air.as_ref().unwrap().coupling.structural_modes()-4;
    let fast=compare(attach(first,true),attach(make(),false),128);
    assert!(fast.system.state()[2*start..2*(start+4)].iter().any(|v|v.abs()>1e-10));
    assert!(!fast.sources.iter().any(|&i|(start..start+4).contains(&i)));
    assert!(fast.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
    assert!(fast.system.membrane_observation(2).unwrap().stretching_energy_j>0.0);
}

#[test]
fn complete_twenty_strand_bank_retains_all_coordinates_with_a_smaller_loaded_border() {
    let e=drum_with_cavity_loss(2,2e-6,true,false,Some(snare::SnareSet::reference(false)),true,
        Stroke {speed_m_s:0.0,position_m:Some([0.06,0.01])},true,None,
        Some(drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}),None,&[],25.0).unwrap();
    let m=e.force.len();let physical=e.system.state().len();
    let source_modes=e.acoustics.as_ref().unwrap().state_modes().to_vec();
    assert!(m>160);let fast=attach(e,true);
    assert_eq!(fast.system.state().len(),physical+32);
    assert_eq!(fast.system.condensed_newton_dimension(),Some(physical-m+1));
    assert_eq!(fast.sources,source_modes);assert_eq!(fast.force.len(),m);
    // No need to run a large idle trajectory to prove preparation kept the
    // bank. The active small-bank test above checks the full nonlinear solve.
}
