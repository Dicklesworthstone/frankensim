use super::*;
use super::super::{drive, config};
use super::super::super::{drum, drum_with_cavity_loss, Stroke, acoustics, muffling};
use fs_couple::render::plate::impact::{ImpactBody, ImpactError, ImpactSystem};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

fn bounds() -> ImpactSubstepConfig { ImpactSubstepConfig { max_depth:8, max_attempts:511 } }

#[test]
fn substep_controls_preserve_physical_inputs_and_do_not_convert_modal_snare_physics() {
    let mut args = ["splash-mic","128","20","--impact-substeps","8","511",
        "--analytic-newton","--strike-speed-m-s","4"].map(String::from).to_vec();
    assert_eq!(option(&mut args).unwrap(),Some(bounds()));
    assert!(super::super::analytic_option(&mut args).unwrap());
    let (positional, stroke)=super::super::super::playing::parse(args).unwrap();
    assert_eq!(positional,["splash-mic","128","20"]);assert_eq!(stroke.speed_m_s,4.0);
    for bad in ["--impact-substeps", "--impact-substeps 8", "--impact-substeps 11 10",
        "--impact-substeps 8 0", "--impact-substeps 8 2048", "--impact-substeps -1 10",
        "--impact-substeps 8 NaN", "--impact-substeps 8 4 --impact-substeps 8 4"] {
        let mut args=bad.split_whitespace().map(String::from).collect();
        assert!(option(&mut args).is_err(),"{bad}");
    }
    assert!(drum(1,2e-6,false,true).unwrap().system.with_impact_substeps(bounds()).is_err());
    // Method order is a numerical choice only. Select the analytic Jacobian
    // after refinement and after an already accepted sample without a reset.
    let (body,_)=ImpactBody::free_mass(0.02,-0.001,0.8).unwrap();
    let mut s=Mechanics::Reference(ImpactSystem::new(vec![body],vec![],vec![],vec![],config(4,2e-6)).unwrap())
        .with_impact_substeps(bounds()).unwrap();
    s.step(&[0.0],&CancelGate::new_clock_free()).unwrap();let before=s.state().to_vec();
    s=s.into_prepared_nonlinear().unwrap().into_analytic_nonlinear().unwrap();
    assert_eq!(s.state(),before);
    let Mechanics::Substepped(s)=s else {panic!("must retain refinement")};assert_eq!(s.samples(),1);
}

fn driven_collision(iterations:usize) -> Mechanics {
    let dt=128e-6;let (body,w)=ImpactBody::free_mass(0.02,-0.0001,1.0).unwrap();
    let contact=Obstacle::new(vec![w],1,1,vec![0.0],vec![1.0],1e8,1.5,
        "physical drive retry regression".into()).unwrap();
    let mut p=ImpactSystem::new(vec![body],vec![contact],vec![],vec![],config(4,dt)).unwrap().prepare_analytic().unwrap();
    p.set_iteration_limit(iterations).unwrap();
    Mechanics::Nonlinear(p).with_stick_drive(
        drive::Program::parse("0,0\n0.000128,0\n0.000256,0.01\n0.000384,0").unwrap(),dt,4,w,1).unwrap()
        .with_impact_substeps(ImpactSubstepConfig {max_depth:8,max_attempts:2}).unwrap()
}

#[test]
fn failed_partial_tick_consumes_neither_the_player_schedule_nor_instrument_time() {
    let gate=CancelGate::new_clock_free();let mut retried=driven_collision(1);let mut clean=driven_collision(50);
    let before=retried.state().to_vec();
    let error=match retried.step(&[0.0],&gate) {Err(e)=>e,Ok(_)=>panic!("bounded partial tick must refuse")};
    assert!(matches!(error,ImpactError::SubstepBudget {accepted_substeps:1,..}),"{error:?}");
    assert_eq!(retried.state(),before);
    let Mechanics::Driven {inner,..}=&mut retried else {panic!("player retained")};
    let Mechanics::Substepped(s)=inner.as_mut() else {panic!("refinement retained")};
    assert_eq!(s.samples(),0);s.set_iteration_limit(50).unwrap();
    let mut work=0.0;
    for tick in 1..=4 {
        let a=retried.step(&[0.0],&gate).unwrap();let b=clean.step(&[0.0],&gate).unwrap();
        assert_eq!(retried.state(),clean.state());assert_eq!(a.time_s,tick as f64*128e-6);
        assert_eq!(a.supplied_work_j,b.supplied_work_j);work+=a.supplied_work_j;
        assert!((a.stored_energy_j-0.01-work).abs()<1e-7);
    }
    assert!(work!=0.0,"the later physical force must actually execute");
    assert!(matches!(retried.step(&[0.0],&gate),Err(ImpactError::Budget)));
}

#[test]
fn refined_two_stick_stretching_drum_retains_cavity_loss_and_pressure_geometry() {
    let dt=acoustics::MECHANICAL_DT;
    let first=Stroke {speed_m_s:4.0,position_m:Some([0.06,0.01])};
    let second=Stroke {speed_m_s:2.5,position_m:Some([-0.05,0.02])};
    let pad=muffling::Muffler {surface:muffling::Surface::Batter,position_m:[0.08,0.01],resistance_n_s_m:0.4};
    let mut e=drum_with_cavity_loss(128,dt,true,false,None,true,first,true,None,None,Some(second),&[pad],20.0).unwrap();
    let before=e.system.state().to_vec();let a=e.observer_a.clone();let b=e.observer_b.clone();
    e.system=e.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();
    let Mechanics::Substepped(s)=&e.system else {panic!("prepared nonlinear")};let energy=s.stored_energy_j();
    assert_eq!(s.state(),before);assert!(e.acoustics.is_some());
    let port=e.second_stick.unwrap();let air=e.air.as_ref().unwrap();
    let solid=air.coupling.structural_modes();assert!(solid>port.coordinate);
    assert!(e.force[solid..].iter().all(|v|*v==0.0));
    let program=|force| drive::Program::parse(&format!("0,0\n{},0\n{},{force}\n{},0",32.0*dt,64.0*dt,96.0*dt)).unwrap();
    let inputs=vec![drive::Input {program:program(0.05),coordinate:0,tip_weight:e.stick_weight},
        drive::Input {program:program(-0.03),coordinate:port.coordinate,tip_weight:port.weight}];
    e.system=e.system.with_stick_drives(inputs,dt,128,e.force.len()).unwrap();
    let gate=CancelGate::new_clock_free();let mut net=0.0;let mut peak=0.0_f64;
    for tick in 1..=128 {
        let f=e.system.step(&e.force,&gate).unwrap();net+=f.supplied_work_j-f.dissipated_energy_j;
        assert_eq!(f.time_s,tick as f64*dt);
        assert!(f.balance_residual_j.abs()<1e-7 && (f.stored_energy_j-energy-net).abs()<1e-6);
        peak=peak.max(air.uniform_pressure(e.system.state()).unwrap().abs());
    }
    assert!(peak>0.0 && e.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
    assert_eq!(e.observer_a,a);assert_eq!(e.observer_b,b);
}
