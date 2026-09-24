use super::*;
use crate::{Stroke,drum_with_spec,drum_with_sticks,mechanics::drive};
use fs_exec::CancelGate;

fn source()->String {include_str!("estimated_drum.fsd").replace("mesh,5,32","mesh,2,16")}
fn varied()->Spec {
    Spec::read(&format!("{}\ntension_variation,batter,100,-50,20,200,-100,150,100\n\
        tension_variation,resonant,-30,60,-10,-100,50,100,-50\n",source())).unwrap()
}
fn stroke()->Stroke {Stroke {speed_m_s:2.0,position_m:Some([0.06,0.01])}}
fn observed(e:&crate::Experiment)->f64 {
    e.observer_a.iter().enumerate().map(|(i,g)|g*e.system.state()[2*i+1]).sum()
}

#[test]
fn tension_records_are_independent_complete_si_tensors_not_pitches_or_defaults() {
    let base=Spec::read(&source()).unwrap();let s=varied();
    assert_eq!(s.heads,base.heads);assert_eq!(s.volume_m3(),base.volume_m3());
    assert_eq!(s.tension_variation[0].constant_n_m,[100.0,-50.0,20.0]);
    assert_eq!(s.tension_variation[1].gradient_n_m2,[-100.0,50.0,100.0,-50.0]);
    let zero=Spec::read(&format!("{}\ntension_variation,batter,0,0,0,0,0,0,0\n",source())).unwrap();
    assert_eq!(zero,base);
    for row in ["tension_variation,unknown,0,0,0,0,0,0,0",
        "tension_variation,batter,1,2,3", "tension_variation,batter,0,0,0,0,0,0,0,0",
        "tension_variation,batter,NaN,0,0,0,0,0,0",
        "tension_variation,batter,0,0,0,0,0,0,inf",
        "tension_variation,batter,0,0,0,0,0,0,0\ntension_variation,batter,1,0,0,0,0,0,0"] {
        assert!(Spec::read(&format!("{}\n{row}\n",source())).is_err(),"{row}");
    }
    let mut invalid=s;invalid.tension_variation[0].gradient_n_m2[0]=1e6;
    // Center tensor is positive. The actual rim must still refuse before modes.
    assert!(invalid.prepare(2e-6,false).is_err());
}

#[test]
fn selected_tension_prepares_actual_modes_and_the_same_reciprocal_air_boundary() {
    let s=varied();let plain=Spec::read(&source()).unwrap();
    let (films,modes)=s.prepare(2e-6,true).unwrap();
    for i in 0..2 {
        let original=plain.head(i).unwrap();
        assert_eq!(films[i].model.m,original.model.m);
        assert_eq!(films[i].model.dof_map,original.model.dof_map);
        assert_ne!(films[i].model.k,original.model.k);
        assert_eq!(films[i].tension_variation,s.tension_variation[i]);
        for mode in &modes[i] {assert!(films[i].modal_area(&mode.phi).unwrap().is_finite());}
    }
    let boundary=crate::acoustics::Boundary::drum(&films,&modes,s.depth_m,s.outer_radius_m).unwrap();
    assert_eq!(boundary.state_modes(),&(1..1+modes.iter().map(Vec::len).sum::<usize>()).collect::<Vec<_>>());
    assert_eq!(films[0].mesh.nodes,films[1].mesh.nodes);
}

#[test]
fn real_tuned_head_contact_changes_motion_and_preserves_energy_in_the_prepared_host() {
    let make=|s|drum_with_spec(256,2e-6,false,true,None,false,stroke(),false,None,Some(s)).unwrap();
    let mut tuned=make(varied());let mut plain=make(Spec::read(&source()).unwrap());
    assert_eq!(&tuned.system.state()[..2],&plain.system.state()[..2]);
    assert_eq!(tuned.stick_weight,plain.stick_weight);
    let initial=0.5*tuned.system.state()[1].powi(2);let mut loss=0.0;let mut changed=0.0_f64;
    let gate=CancelGate::new_clock_free();
    for _ in 0..256 {
        let f=tuned.system.step(&tuned.force,&gate).unwrap();plain.system.step(&plain.force,&gate).unwrap();
        loss+=f.dissipated_energy_j;
        assert!((f.stored_energy_j+loss-initial).abs()<1e-7);
        changed=changed.max((observed(&tuned)-observed(&plain)).abs());
    }
    assert!(changed>1e-10,"the field must change physical head motion, not only input metadata");
    assert!((tuned.system.state()[1]*tuned.stick_weight-stroke().speed_m_s).abs()>1e-6);
}

#[test]
fn nonuniform_two_head_stretching_keeps_two_hand_drive_cavity_history_and_exact_retry() {
    let second=Stroke {speed_m_s:1.0,position_m:Some([-0.05,0.02])};
    let make=|| {
        let mut e=drum_with_sticks(192,2e-6,false,false,None,true,stroke(),true,None,Some(varied()),Some(second)).unwrap();
        let other=e.second_stick.unwrap();
        e.system=e.system.into_analytic_nonlinear().unwrap();
        e.system=e.system.with_stick_drives(vec![
            drive::Input {program:drive::Program::parse("0,0\n0.0001,0.02\n0.0003,0").unwrap(),coordinate:0,tip_weight:e.stick_weight},
            drive::Input {program:drive::Program::parse("0,0\n0.00015,-0.01\n0.0003,0").unwrap(),coordinate:other.coordinate,tip_weight:other.weight},
        ],2e-6,192,e.force.len()).unwrap();e
    };
    let mut a=make();let mut b=make();let gate=CancelGate::new_clock_free();let mut stretch=0.0_f64;
    for tick in 0..192 {
        if tick==96 {
            let old=a.system.state().to_vec();let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(a.system.step(&a.force,&cancel).is_err());assert_eq!(a.system.state(),old);
        }
        let f=a.system.step(&a.force,&gate).unwrap();b.system.step(&b.force,&gate).unwrap();
        assert_eq!(a.system.state(),b.system.state());assert!(f.balance_residual_j.abs()<2e-7);
        stretch=stretch.max(a.system.membrane_observation(1).unwrap().stretching_energy_j);
        let (top,bottom)=a.air.as_ref().unwrap().points(a.system.state()).unwrap();
        assert!(top.is_finite() && bottom.is_finite());
    }
    assert!(stretch>0.0);assert!(observed(&a).abs()>0.0);
}
