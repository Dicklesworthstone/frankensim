use crate::{splash_with_sticks,splash_with_mufflers,Stroke,mechanics::Mechanics,
    muffling::{Muffler,Surface},specimen::Specimen};
use crate::mechanics::drive::{Input,Program,StickDrive};
use fs_couple::render::plate::impact::ImpactSubstepConfig;
use fs_exec::CancelGate;

fn first() -> Stroke {Stroke {speed_m_s:0.8,position_m:Some([0.06,0.01])}}
fn second() -> Stroke {Stroke {speed_m_s:0.8,position_m:Some([-0.05,0.02])}}
fn bounds() -> ImpactSubstepConfig {ImpactSubstepConfig {max_depth:8,max_attempts:511}}

#[test]
fn two_impacts_share_one_curved_shell_and_preserve_stand_memory_and_observer_addresses() {
    let dt=2e-6;let count=192;
    let pad=Muffler {surface:Surface::Shell,position_m:[0.075,0.0],resistance_n_s_m:0.2};
    let mut one=splash_with_mufflers(count,dt,true,first(),None,&[pad]).unwrap();
    let mut two=splash_with_sticks(count,dt,true,first(),None,&[pad],Some(second())).unwrap();
    let port=two.second_stick.unwrap();let old_n=one.force.len();
    assert_eq!(port.coordinate,old_n);assert_eq!(two.force.len(),old_n+1);
    assert_eq!(&two.system.state()[..2*old_n],&one.system.state()[..2*old_n]);
    // The new mechanical pair precedes the Kelvin tail; tail VALUES do not move.
    assert_eq!(&two.system.state()[2*(old_n+1)..],&one.system.state()[2*old_n..]);
    for (a,b) in [(&one.observer_a,&two.observer_a),(&one.observer_b,&two.observer_b)] {
        assert_eq!(&b[..old_n],a.as_slice());assert_eq!(b[old_n],0.0);
    }
    assert!(one.acoustics.is_some() && two.acoustics.is_some());
    let Mechanics::Reference(a)=&one.system else {unreachable!()};
    let Mechanics::Reference(b)=&two.system else {unreachable!()};
    for i in 0..6 {assert_eq!(a.felt_history(i),b.felt_history(i));assert_eq!(a.felt_observation(i),b.felt_observation(i));}
    let initial=b.stored_energy_j();
    assert!((initial-a.stored_energy_j()-0.5*(second().speed_m_s/port.weight).powi(2)).abs()<1e-12);
    one.system=one.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();
    two.system=two.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();
    let gate=CancelGate::new_clock_free();let mut loss=0.0;let mut difference=0.0_f64;
    for tick in 1..=count {
        one.system.step(&one.force,&gate).unwrap();let f=two.system.step(&two.force,&gate).unwrap();
        loss+=f.dissipated_energy_j;
        assert_eq!(f.time_s,tick as f64*dt);
        assert!((f.stored_energy_j+loss-initial).abs()<1e-6 && f.balance_residual_j.abs()<1e-7);
        for (a,b) in one.system.state()[2..2*old_n].iter().zip(&two.system.state()[2..2*old_n]) {
            difference=difference.max((a-b).abs());
        }
    }
    assert!(difference>1e-12,"second contact must alter the same shell, not an independently mixed voice");
    assert!((two.system.state()[2*old_n+1]*port.weight-second().speed_m_s).abs()>1e-8,
        "the surface must react on the second stick");
    let Mechanics::Substepped(s)=&two.system else {panic!("one nonlinear joint owner")};
    for i in 0..6 {assert!(s.felt_history(i).is_some());}
    assert_eq!(s.samples(),count);
}

#[test]
fn supplied_shell_accepts_independent_signed_player_inputs_without_consuming_failed_ticks() {
    let dt=2e-6;let count=192;
    let make=||splash_with_sticks(count,dt,false,first(),Some(Specimen::reference()),&[],Some(second())).unwrap();
    let mut driven=make();let mut manual=make();let n=driven.force.len();let port=driven.second_stick.unwrap();
    // Distinct schedules; source geometry/stand material are never reconstructed during playback.
    let inputs=||vec![Input {program:Program::parse("0,0\n0.000128,0.02\n0.000384,0").unwrap(),
        coordinate:0,tip_weight:driven.stick_weight},Input {
        program:Program::parse("0,0\n0.000096,-0.01\n0.000256,0").unwrap(),
        coordinate:port.coordinate,tip_weight:port.weight}];
    let mut staging=StickDrive::new_inputs(inputs(),dt,count,n).unwrap();
    let programs=inputs();
    driven.system=driven.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap()
        .with_stick_drives(programs,dt,count,n).unwrap();
    manual.system=manual.system.into_analytic_nonlinear().unwrap().with_impact_substeps(bounds()).unwrap();
    let gate=CancelGate::new_clock_free();let mut work=0.0_f64;
    for tick in 0..count {
        if tick==48 {
            let before=driven.system.state().to_vec();let mut bad=vec![0.0;n];bad[port.coordinate]=1e7;
            assert!(driven.system.step(&bad,&gate).is_err());assert_eq!(driven.system.state(),before);
        }
        let f=driven.system.step(&driven.force,&gate).unwrap();
        let expected=manual.system.step(staging.forces(&manual.force).unwrap(),&gate).unwrap();staging.accept();
        assert_eq!(driven.system.state(),manual.system.state());
        assert_eq!(f.supplied_work_j,expected.supplied_work_j);work+=f.supplied_work_j.abs();
        assert_eq!(f.time_s,(tick+1) as f64*dt);assert!(f.balance_residual_j.abs()<1e-7);
    }
    assert!(work>0.0);
}

#[test]
fn absent_second_stick_preserves_single_shell_and_hole_or_outside_stations_refuse() {
    let mut old=splash_with_mufflers(4,2e-6,false,first(),None,&[]).unwrap();
    let mut new=splash_with_sticks(4,2e-6,false,first(),None,&[],None).unwrap();
    assert!(new.second_stick.is_none());assert_eq!(old.force,new.force);
    let gate=CancelGate::new_clock_free();
    for _ in 0..4 {
        let a=old.system.step(&old.force,&gate).unwrap();let b=new.system.step(&new.force,&gate).unwrap();
        assert_eq!(old.system.state(),new.system.state());assert_eq!(a.stored_energy_j,b.stored_energy_j);
    }
    for position in [None,Some([0.0,0.0]),Some([0.2,0.0])] {
        assert!(splash_with_sticks(1,2e-6,false,first(),None,&[],Some(Stroke {
            speed_m_s:0.8,position_m:position })).is_err(),"{position:?}");
    }
}
