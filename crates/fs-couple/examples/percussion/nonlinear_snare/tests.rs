use super::*;
use super::super::{drum_with_cavity_loss,drum_spec,snare::SnareSet,Stroke,Mechanics,
    Experiment,acoustics,cavity,mechanics,muffling};
use fs_exec::CancelGate;
use fs_couple::render::plate::impact::ImpactSubstepConfig;

fn energy(system: &Mechanics) -> f64 {
    match system {
        Mechanics::Reference(s) => s.stored_energy_j(),
        Mechanics::Nonlinear(s) => s.stored_energy_j(),
        Mechanics::Substepped(s) => s.stored_energy_j(),
        Mechanics::Driven { inner, .. } => energy(inner),
        Mechanics::Prepared(_) => panic!("this fixture requires nonlinear mechanics"),
    }
}

fn specimen()->drum_spec::Spec {
    drum_spec::Spec {radial_intervals:2,azimuths:8,..drum_spec::Spec::reference()}
}
// Deliberate installed interference in this small runtime fixture. It stores
// initial contact energy; no equilibrium, preload calibration or stock setting
// is inferred. Production retains its unchanged positive 20-micrometre gap.
fn wires()->SnareSet {SnareSet {strands:2,modes_per_strand:2,contact_cells:4,
    clearance_m:-2e-6,..SnareSet::reference(false)}}
fn build(spec:SnareSet,air:bool,second:bool,audio:bool)->Experiment {
    let dt=if audio {acoustics::MECHANICAL_DT}else{2e-6};
    drum_with_cavity_loss(256,dt,audio,false,Some(spec),true,
        Stroke {speed_m_s:0.0,position_m:Some([0.06,0.01])},air,None,Some(specimen()),
        second.then_some(Stroke {speed_m_s:0.0,position_m:Some([-0.05,0.02])}),
        &[muffling::Muffler {surface:muffling::Surface::Batter,position_m:[0.07,0.01],resistance_n_s_m:0.1}],
        if air {25.0}else{0.0}).unwrap()
}

#[test]
fn nonlinear_snare_selection_is_explicit_and_composes_with_existing_front_doors() {
    let mut args=vec!["snare-mic".into(),"48".into(),"--head-stretching".into(),"--analytic-newton".into()];
    assert!(option(&mut args).unwrap());assert_eq!(args,["snare-mic","48","--analytic-newton"]);
    assert!(!option(&mut args).unwrap());assert!(option(&mut vec!["--head-stretching".into();2]).is_err());
    for command in ["snare","snare-off","snare-wav","snare-off-wav","snare-mic","snare-off-mic"] {
        admit_command(true,command).unwrap();admit_prepared_command(true,true,command).unwrap();
        assert!(admit_prepared_command(true,false,command).is_err());
        cavity::admit_command(true,command).unwrap();
        if command.ends_with("-mic") {acoustics::stereo::admit_command(Some([0.3,0.0,0.5]),command).unwrap();}
    }
    for command in ["drum-modal","splash","drum-stretch","unknown"] {assert!(admit_command(true,command).is_err());}
    assert!(admit_image(true,true,true).is_err());assert!(admit_image(false,true,false).is_err());
    admit_image(false,true,true).unwrap();admit_image(true,true,false).unwrap();
}

#[test]
fn nonlinear_heads_retain_the_complete_twenty_strand_bank_and_acoustic_geometry() {
    let e=build(SnareSet::reference(false),true,true,true);
    assert!(matches!(&e.system,Mechanics::Reference(_)));assert!(e.acoustics.is_some());
    assert!(e.system.membrane_observation(1).is_some() && e.system.membrane_observation(2).is_some());
    let air=e.air.as_ref().unwrap();let structural=air.coupling.structural_modes();
    let second=e.second_stick.unwrap();assert_eq!(structural,second.coordinate+1+160);
    let start=second.coordinate+1;
    assert!(e.system.state()[2*start..].iter().all(|v|*v==0.0));
    assert!(e.observer_a[start..].iter().chain(&e.observer_b[start..]).all(|v|*v==0.0));
    assert!(e.pressure.as_ref().unwrap().areas[start..].iter().all(|v|*v==0.0));
    assert!(air.coupling.total_modes()>structural);assert_eq!(e.force.len(),air.coupling.total_modes());
    // The 256-coordinate work ceiling retains the bank, not a truncated proxy.
    assert!(structural>160 && e.force.len()<=256);
}

#[test]
fn wire_contacts_exchange_energy_with_both_stretching_heads_in_the_same_state() {
    let mut active=build(wires(),false,false,false);
    let mut inactive=build(SnareSet {contact_stiffness_per_length:0.0,..wires()},false,false,false);
    assert_eq!(active.system.state(),inactive.system.state());
    let initial=energy(&active.system);let start=active.force.len()-4;
    active.system=active.system.into_analytic_nonlinear().unwrap();
    inactive.system=inactive.system.into_analytic_nonlinear().unwrap();
    let gate=CancelGate::new_clock_free();let mut losses=0.0;let mut peak=0.0_f64;let mut stretch=0.0_f64;
    for _ in 0..128 {
        let f=active.system.step(&active.force,&gate).unwrap();inactive.system.step(&inactive.force,&gate).unwrap();
        losses+=f.dissipated_energy_j;assert!(f.balance_residual_j.abs()<1e-7);
        peak=peak.max(active.system.state()[2*start..].iter().map(|v|v.abs()).fold(0.0,f64::max));
        for head in [1,2] {
            let observation=active.system.membrane_observation(head).unwrap();
            assert!(observation.maximum_slope<=0.2);stretch=stretch.max(observation.stretching_energy_j);
        }
    }
    assert!(peak>1e-10 && stretch>0.0 && losses>0.0);
    assert!(inactive.system.state()[2..].iter().all(|v|*v==0.0));
    assert!((energy(&active.system)-initial+losses).abs()<1e-6);
    assert!(active.system.state()[2..2*start].iter().any(|v|v.abs()>1e-12));
}

#[test]
fn nonlinear_snare_keeps_cavity_two_player_ports_refinement_and_exact_retry() {
    let prepare=|| {
        let mut e=build(wires(),true,true,false);let second=e.second_stick.unwrap();
        e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(
            ImpactSubstepConfig {max_depth:4,max_attempts:31}).unwrap();
        let inputs=vec![mechanics::drive::Input {program:mechanics::drive::Program::parse("0,0.02\n0.000512,0.02").unwrap(),coordinate:0,tip_weight:e.stick_weight},
            mechanics::drive::Input {program:mechanics::drive::Program::parse("0,-0.01\n0.000512,-0.01").unwrap(),coordinate:second.coordinate,tip_weight:second.weight}];
        e.system=e.system.with_stick_drives(inputs,2e-6,256,e.force.len()).unwrap();e
    };
    let mut e=prepare();let mut clean=prepare();let initial=energy(&e.system);
    let gate=CancelGate::new_clock_free();let cancel=CancelGate::new_clock_free();cancel.request();
    let before=e.system.state().to_vec();assert!(e.system.step(&e.force,&cancel).is_err());
    assert_eq!(e.system.state(),before);let mut net=0.0;let mut work=0.0;let mut pressure=0.0_f64;
    for tick in 1..=32 {
        let f=e.system.step(&e.force,&gate).unwrap();clean.system.step(&clean.force,&gate).unwrap();
        assert_eq!(e.system.state(),clean.system.state());assert_eq!(f.time_s,tick as f64*2e-6);
        net+=f.supplied_work_j-f.dissipated_energy_j;work+=f.supplied_work_j.abs();
        pressure=pressure.max(e.air.as_ref().unwrap().uniform_pressure(e.system.state()).unwrap().abs());
        assert!((f.stored_energy_j-initial-net).abs()<1e-6);
    }
    assert!(work>0.0 && pressure>0.0);
    // The new head law is not an aperture-radiation implementation.
    let neck=cavity::NeckOptions {radius_m:0.005,effective_length_m:0.012,resistance_pa_s_m3:5000.0,
        azimuth_rad:0.4,axial_position_m:0.08};
    assert!(drum_with_cavity_loss(2,acoustics::MECHANICAL_DT,true,false,Some(wires()),true,
        Stroke::default(),true,Some(neck),Some(specimen()),None,&[],0.0).is_err());
}
