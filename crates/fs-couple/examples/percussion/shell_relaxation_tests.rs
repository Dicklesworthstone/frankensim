use super::*;
use crate::{Experiment,Mechanics,Stroke,head_relaxation,mallets,shaft_playing,specimen};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_exec::CancelGate;

const CARD:&str="frankensim-shell-relaxation-v1\nintrinsic_loss,replace\ninitial,relaxed\nband_hz,0,5000\nshell,single\nbranch,single,0.05,0.0002\n";
fn shell()->specimen::Specimen {let mut s=specimen::Specimen::reference();s.azimuths=8;s}
fn energy(system:&Mechanics)->f64 {
    match system {Mechanics::Reference(s)=>s.stored_energy_j(),Mechanics::Nonlinear(s)=>s.stored_energy_j(),
        Mechanics::Substepped(s)=>s.stored_energy_j(),Mechanics::Driven{inner,..}=>energy(inner),
        Mechanics::Prepared(_)=>panic!("shell material requires its existing joint owner")}
}
fn build(material:Option<&Spec>)->Experiment {
    crate::splash_with_material(512,2e-6,true,Stroke{speed_m_s:0.8,position_m:Some([0.06,0.01])},
        Some(shell()),&[],Some(Stroke{speed_m_s:0.,position_m:Some([-0.05,0.02])}),None,
        &mallets::Selection::default(),&shaft_playing::Selection::default(),material).unwrap()
}

#[test]
fn explicit_loss_replacement_and_bounded_physical_input_are_required_on_every_shell_command() {
    let material=Spec::read(CARD,false).unwrap();material.admit_windows(&[[50.,1200.]]).unwrap();
    assert!(material.admit_windows(&[[50.,6000.]]).is_err());
    assert!(material.admit_windows(&[[50.,1200.];2]).is_err());
    for command in ["splash","splash-wav","splash-mic","hihat","hihat-wav","hihat-mic"] {
        admit_command(true,command).unwrap();
        let mut args=vec![command.into(),"--shell-relaxation".into(),"supplied.fssr".into(),"512".into()];
        assert_eq!(option(&mut args).unwrap().as_deref(),Some("supplied.fssr"));
        assert_eq!(args,[command,"512"]);
    }
    for command in ["drum","drum-modal-mic","snare","unknown"] {assert!(admit_command(true,command).is_err());}
    for bad in [CARD.replace("intrinsic_loss,replace\n",""),CARD.replace("replace","add"),
        CARD.replace("initial,relaxed\n",""),CARD.replace("shell,single\n",""),
        CARD.replace("0.05","-1"),CARD.replace("0.0002","NaN"),CARD.replace("0.0002","0"),
        format!("{CARD}shell,single\n"),format!("{CARD}intrinsic_loss,replace\n"),
        format!("{CARD}{}","branch,single,0.01,0.001\n".repeat(8))] {
        assert!(Spec::read(&bad,false).is_err());
    }
    assert!(Spec::read(CARD,true).is_err());
    let mut missing=vec!["--shell-relaxation".into()];assert!(option(&mut missing).is_err());
    let mut duplicate=vec!["--shell-relaxation".into(),"a".into(),"--shell-relaxation".into(),"b".into()];
    assert!(option(&mut duplicate).is_err());
}

#[test]
fn actual_free_shell_translation_has_no_material_stiffness_or_memory_coordinate() {
    let (_,reduction)=crate::shell_prepare::prepare(&shell(),2e-6).unwrap();
    let n=reduction.mode_count();let bending=reduction.bending_stiffness(&vec![1.;reduction.facet_count()],2_000_000).unwrap();
    assert_eq!(reduction.omegas()[0],0.);
    assert!(bending[..n].iter().all(|v|*v==0.));
    assert!((0..n).all(|i|bending[i*n]==0.));
    let material=Spec::read(CARD,false).unwrap();let e=build(Some(&material));
    let sources=e.acoustics.as_ref().unwrap().state_modes();
    assert_eq!(sources.len(),n);
    assert_eq!(head_relaxation::observation(&e.system).states,n-1);
    assert_eq!(head_relaxation::observation(&e.system).stored_energy_j,0.);
    assert_eq!(e.system.state().len(),2*e.force.len()+6+n-1);
    assert_eq!(sources,(1..e.second_stick.unwrap().coordinate).collect::<Vec<_>>());
}

#[test]
fn contact_excites_physical_shell_memory_with_complete_energy_and_retry_accounting() {
    let material=Spec::read(CARD,false).unwrap();
    let elastic=Spec::read(&CARD.replace("branch,single,0.05,0.0002\n",""),false).unwrap();
    // Both inputs explicitly replace intrinsic modal loss. Their only material
    // difference is the supplied bending-memory branch, not a damping retune.
    let mut relaxing=build(Some(&material));let mut bare=build(Some(&elastic));
    let prefix=bare.system.state().len();assert_eq!(&relaxing.system.state()[..prefix],bare.system.state());
    assert_eq!(relaxing.force,bare.force);assert_eq!(relaxing.observer_a,bare.observer_a);
    let initial=energy(&relaxing.system);
    let prepare=|s:Mechanics|s.into_analytic_nonlinear().unwrap().with_impact_substeps(
        ImpactSubstepConfig{max_depth:8,max_attempts:511}).unwrap();
    relaxing.system=prepare(relaxing.system);bare.system=prepare(bare.system);
    let gate=CancelGate::new_clock_free();let mut loss=0.;let (mut changed,mut memory_peak)=(0.0_f64,0.0_f64);
    for tick in 0..512 {
        if tick==192 {
            let state=relaxing.system.state().to_vec();let memory=head_relaxation::observation(&relaxing.system);
            let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(relaxing.system.step(&relaxing.force,&cancel).is_err());
            assert_eq!(relaxing.system.state(),state);
            assert_eq!(head_relaxation::observation(&relaxing.system),memory);
        }
        let frame=relaxing.system.step(&relaxing.force,&gate).unwrap();bare.system.step(&bare.force,&gate).unwrap();
        loss+=frame.dissipated_energy_j;assert!(frame.balance_residual_j.abs()<1e-7);
        for (a,b) in relaxing.system.state().iter().zip(bare.system.state()) {changed=changed.max((a-b).abs());}
        let memory=head_relaxation::observation(&relaxing.system);
        memory_peak=memory_peak.max(memory.stored_energy_j);assert!(memory.dissipated_power_w>=0.);
    }
    assert!(changed>0. && memory_peak>0. && loss>0.);
    assert!((energy(&relaxing.system)+loss-initial).abs()<1e-6);
}
