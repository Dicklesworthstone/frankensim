use super::*;
use crate::{BodyPotential, ImpactSystem, Stroke, config, mallets, shaft_playing, shell_prepare, specimen};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_exec::CancelGate;

fn specimen()->specimen::Specimen {
    let mut s=specimen::Specimen::reference();s.azimuths=8;s
}

#[test]
fn rocking_is_explicit_bounded_and_only_admitted_for_cymbals() {
    let mut args=vec!["splash".into(),"--cymbal-rocking".into(),"0.03".into()];
    assert_eq!(option(&mut args).unwrap(),Some(0.03));assert_eq!(args,["splash"]);
    for command in ["splash","splash-wav","splash-mic","hihat","hihat-wav","hihat-mic"] {
        admit_command(Some(0.03),command).unwrap();
    }
    for command in ["drum","snare","drum-modal"] {assert!(admit_command(Some(0.03),command).is_err());}
    for text in ["0","-0.01","0.101","NaN","inf"] {
        let mut args=vec!["--cymbal-rocking".into(),text.into()];let before=args.clone();
        assert!(option(&mut args).is_err());assert_eq!(args,before);
    }
}

#[test]
fn excessive_rigid_rotation_refuses_initial_and_every_time_image_without_consuming_input() {
    let dt=2e-6;let cap=1e-5;
    let (_,shell)=shell_prepare::prepare_rocking(&specimen(),dt,cap).unwrap();
    let n=shell.mode_count();let mut probe=vec![0.;n];probe[1]=1e-10;
    let rotation=shell.rigid_rotation(&probe).unwrap().unwrap();
    let radians_per_q=rotation.iter().map(|v|v*v).sum::<f64>().sqrt()/probe[1];
    assert!(radians_per_q>0.);
    let outside=2.*cap/radians_per_q;
    let make_body=|| {
        let mut body=crate::zero_body(BodyPotential::Shell(shell.clone()),shell.omegas());
        body.damping_per_s.fill(0.);body
    };
    let mut bad=make_body();bad.initial[1].displacement_m_sqrt_kg=outside;
    assert!(ImpactSystem::new(vec![bad],vec![],vec![],vec![],config(2,dt)).is_err());
    for image in 0..4 {
        let mut body=make_body();body.initial[1].velocity_m_sqrt_kg_per_s=outside/dt;
        let reference=ImpactSystem::new(vec![body],vec![],vec![],vec![],config(2,dt)).unwrap();
        assert_eq!(reference.shell_rotation(0).unwrap(),[0.;3]);
        let mut system=Mechanics::Reference(reference);
        if image>0 {system=if image==1 {system.into_prepared_nonlinear().unwrap()}
            else {system.into_analytic_nonlinear().unwrap()};}
        if image==3 {system=system.with_impact_substeps(ImpactSubstepConfig{max_depth:3,max_attempts:15}).unwrap();}
        system=system.with_stick_drives(vec![crate::mechanics::drive::Input {
            program:crate::mechanics::drive::Program::parse("0,0\n0.000001,0.001\n0.000004,0").unwrap(),
            coordinate:0,tip_weight:1.,
        }],dt,2,n).unwrap();
        let before=system.state().to_vec();let force=vec![0.;n];let gate=CancelGate::new_clock_free();
        let staged=match &mut system {Mechanics::Driven{drive,..}=>drive.forces(&force).unwrap().to_vec(),
            _=>unreachable!()};
        for _ in 0..2 {
            let error=system.step(&force,&gate).err().expect("excessive rotation must refuse");
            assert!(error.to_string().contains("linearized rigid-rotation limit"),"{error}");
            assert_eq!(system.state(),before,"an angle-domain refusal cannot retain a partial step");
            assert_eq!(observation(&system,0).unwrap(),[0.;3]);
            if let Mechanics::Driven{drive,inner}=&mut system {
                assert_eq!(drive.forces(&force).unwrap(),staged);
                let samples=match inner.as_ref() {
                    Mechanics::Reference(s)=>s.samples(),Mechanics::Nonlinear(s)=>s.samples(),
                    Mechanics::Substepped(s)=>s.samples(),_=>unreachable!(),
                };
                assert_eq!(samples,0);
            }
        }
    }
}

#[test]
fn an_off_center_hit_rocks_the_cymbal_against_all_six_original_stand_contacts() {
    let dt=2e-6;let steps=1024;
    let stroke=Stroke{speed_m_s:0.4,position_m:Some([0.068,0.01])};
    let make=|rocking|crate::splash_with_rocking(steps,dt,false,stroke,Some(specimen()),&[],None,
        None,&mallets::Selection::default(),&shaft_playing::Selection::default(),None,rocking).unwrap();
    let fixed=make(None);let mut e=make(Some(0.03));
    assert_eq!(e.force.len(),fixed.force.len()+2);
    assert!(observation(&fixed.system,1).is_none());
    if let Mechanics::Reference(system)=&e.system {
        for i in 0..6 {assert!(system.felt_history(i).is_some());}
        assert!(system.felt_history(6).is_none());
    }
    let initial=match &e.system {Mechanics::Reference(s)=>s.stored_energy_j(),_=>unreachable!()};
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
    let gate=CancelGate::new_clock_free();let mut peak=0.0_f64;let mut loss=0.;
    for _ in 0..steps {
        let frame=e.system.step(&e.force,&gate).unwrap();loss+=frame.dissipated_energy_j;
        assert!((frame.stored_energy_j+loss-initial).abs()<1e-6);
        let r=observation(&e.system,1).unwrap();
        peak=peak.max(r.iter().map(|v|v*v).sum::<f64>().sqrt());
    }
    assert!(peak>1e-8 && peak<=0.03,"the off-center impulse must excite real bounded angular inertia");
    assert!(loss>0.,"the same stand felt retains physical loss while rocking");
}
