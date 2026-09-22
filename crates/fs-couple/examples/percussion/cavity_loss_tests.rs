//! Actual drum compositions exercise the same front door as the CLI.
use super::*;
use super::super::{Mechanics, Stroke, drum_with_cavity_loss, drum_with_mufflers, muffling, snare};
use fs_couple::render::plate::impact::ImpactError;

fn neck() -> NeckOptions {
    NeckOptions { radius_m:0.005, effective_length_m:0.012, resistance_pa_s_m3:5000.0,
        azimuth_rad:0.4, axial_position_m:0.08 }
}
fn stroke() -> Stroke { Stroke { speed_m_s:2.0, position_m:Some([0.06,0.01]) } }

#[test]
fn explicit_acoustic_loss_controls_preserve_other_inputs_and_refuse_wrong_images() {
    let mut args:Vec<String>=["snare","256","--cavity-modes","--cavity-drag-per-s","50",
        "--cavity-neck","0.005","0.012","5000","0.4","0.08",
        "--strike-speed-m-s","2"].map(String::from).to_vec();
    let drag=drag_option(&mut args).unwrap(); assert_eq!(drag,Some(50.0));
    let vent=neck_option(&mut args).unwrap(); assert_eq!(vent,Some(neck()));
    assert!(option(&mut args).unwrap());
    let (positional,playing)=super::super::playing::parse(args).unwrap();
    assert_eq!(positional,["snare","256"]); assert_eq!(playing.speed_m_s,2.0);
    for command in ["drum","drum-stretch","drum-modal","snare","snare-off"] {
        assert!(admit_drag_command(drag,true,command).is_ok());
        assert!(admit_neck_command(vent,true,command).is_ok());
        for suffix in ["-wav","-mic"] {
            let command=format!("{command}{suffix}");
            assert!(admit_drag_command(drag,true,&command).is_ok());
            assert!(admit_neck_command(vent,true,&command).is_err());
        }
    }
    assert!(admit_drag_command(Some(0.0),false,"snare").is_err());
    assert!(admit_drag_command(drag,true,"splash").is_err());
    for text in ["--cavity-drag-per-s", "--cavity-drag-per-s -1", "--cavity-drag-per-s NaN",
        "--cavity-drag-per-s inf", "--cavity-drag-per-s 100001",
        "--cavity-drag-per-s 1 --cavity-drag-per-s 2"] {
        let mut args:Vec<String>=text.split_whitespace().map(String::from).collect();
        let original=args.clone(); assert!(drag_option(&mut args).is_err()); assert_eq!(args,original);
    }
    // Public composition refuses before eigensolve, not just at CLI parsing.
    assert!(drum_with_cavity_loss(1,2e-6,false,true,None,false,stroke(),false,
        None,None,None,&[],50.0).is_err());
    assert!(drum_with_cavity_loss(1,2e-6,true,true,Some(snare::SnareSet::reference(false)),
        false,stroke(),true,Some(neck()),None,None,&[],50.0).is_err());
}

#[test]
fn vented_twenty_strand_two_stick_snare_retains_mufflers_and_real_pressure_work() {
    let wires=snare::SnareSet::reference(false);
    let second=Stroke { speed_m_s:1.6, position_m:Some([-0.05,0.02]) };
    let pad=muffling::Muffler { surface:muffling::Surface::Batter,
        position_m:[0.07,0.01], resistance_n_s_m:0.1 };
    let make=|resistance| drum_with_cavity_loss(385,2e-6,false,true,Some(wires),false,
        stroke(),true,Some(NeckOptions {resistance_pa_s_m3:resistance,..neck()}),
        None,Some(second),&[pad],500.0).unwrap();
    let mut ideal=make(0.0); let mut lossy=make(5000.0);
    let Mechanics::Prepared(system)=&lossy.system else {panic!("must retain prepared snare solver");};
    assert_eq!(system.contact_count(),2+wires.strands*wires.contact_cells);
    let initial=system.frame().stored_energy_j;
    let solid=lossy.air.as_ref().unwrap().coupling.structural_modes();
    assert!(solid>160 && solid<256);
    assert_eq!(lossy.system.state(),ideal.system.state());
    let second=lossy.second_stick.unwrap();
    assert_eq!(solid-(second.coordinate+1),wires.strands*wires.modes_per_strand);
    assert!(lossy.system.state()[2*(second.coordinate+1)..].iter().all(|v| *v==0.0));
    assert!(lossy.observer_a[second.coordinate..].iter().all(|v| *v==0.0));
    assert!(lossy.force[solid..].iter().all(|v| *v==0.0));
    let gate=CancelGate::new_clock_free(); let mut loss=0.0; let mut changed=0.0_f64;
    let mut peak_flow=0.0_f64; let mut neck_power=0.0_f64;
    for _ in 0..384 {
        ideal.system.step(&ideal.force,&gate).unwrap();
        let f=lossy.system.step(&lossy.force,&gate).unwrap();
        loss+=f.dissipated_energy_j;
        assert_eq!(f.supplied_work_j,0.0);
        assert!(f.balance_residual_j.abs()<1e-7);
        assert!((f.stored_energy_j+loss-initial).abs()<1e-6);
        let x=lossy.system.state(); let air=lossy.air.as_ref().unwrap();
        let vent=air.coupling.neck_observation(x,0).unwrap();
        assert_eq!(vent.coordinate,lossy.force.len()-1);
        assert!((vent.resistive_pressure_drop_pa-5000.0*vent.volume_flow_m3_s).abs()<1e-12);
        peak_flow=peak_flow.max(vent.volume_flow_m3_s.abs());
        neck_power=neck_power.max(vent.dissipated_power_w);
        let volume=lossy.pressure.as_ref().unwrap();
        let expected=super::super::cavity_pressure(volume,x)
            -volume.bulk_modulus_pa/volume.volume_m3*vent.displaced_volume_m3;
        assert!((air.uniform_pressure(x).unwrap()-expected).abs()<1e-7*(1.0+expected.abs()));
        changed=changed.max(x[..2*second.coordinate].iter().zip(ideal.system.state())
            .map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max));
    }
    assert!(peak_flow>1e-10 && neck_power>0.0 && changed>1e-14);
    let before=lossy.system.state().to_vec();
    let cancelled=CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(lossy.system.step(&lossy.force,&cancelled),Err(ImpactError::Cancelled)));
    let mut invalid=lossy.force.clone(); invalid[solid]=f64::NAN;
    assert!(lossy.system.step(&invalid,&gate).is_err()); assert_eq!(lossy.system.state(),before);
    let f=lossy.system.step(&lossy.force,&gate).unwrap();
    assert_eq!(f.time_s,385.0*2e-6);
}

#[test]
fn acoustic_drag_changes_actual_sealed_head_motion_in_both_execution_images() {
    let gate=CancelGate::new_clock_free();
    for prepared in [false,true] {
        let make=|drag| drum_with_cavity_loss(256,2e-6,false,prepared,None,!prepared,
            stroke(),true,None,None,None,&[],drag).unwrap();
        let mut zero=make(0.0); let mut lossy=make(2000.0);
        assert_eq!(zero.system.state(),lossy.system.state());
        assert_eq!(zero.observer_a,lossy.observer_a);
        let solid=lossy.air.as_ref().unwrap().coupling.structural_modes();
        let mut changed=0.0_f64; let mut pressure_change=0.0_f64;
        for _ in 0..256 {
            zero.system.step(&zero.force,&gate).unwrap();
            let f=lossy.system.step(&lossy.force,&gate).unwrap();
            assert!(f.balance_residual_j.abs()<1e-7);
            changed=changed.max(zero.system.state()[..2*solid].iter().zip(lossy.system.state())
                .map(|(a,b)|(a-b).abs()).fold(0.0_f64,f64::max));
            let a=zero.air.as_ref().unwrap().points(zero.system.state()).unwrap().0;
            let b=lossy.air.as_ref().unwrap().points(lossy.system.state()).unwrap().0;
            pressure_change=pressure_change.max((a-b).abs());
        }
        assert!(changed>1e-12 && pressure_change>1e-6,"loss must act on physical gas/head motion");
        if !prepared {assert!(lossy.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);}
    }
}

#[test]
fn explicit_zero_drag_keeps_the_original_prepared_trajectory_and_sealed_audio_is_admitted() {
    let mut old=drum_with_mufflers(192,2e-6,false,true,None,false,stroke(),true,
        None,None,None,&[]).unwrap();
    let mut new=drum_with_cavity_loss(192,2e-6,false,true,None,false,stroke(),true,
        None,None,None,&[],0.0).unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..192 {
        let a=old.system.step(&old.force,&gate).unwrap();
        let b=new.system.step(&new.force,&gate).unwrap();
        assert_eq!(old.system.state(),new.system.state());
        assert_eq!(a.stored_energy_j,b.stored_energy_j);
        assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
    }
    // Construct the real undeformed boundary, but do not claim that this test
    // exercises the expensive BEM fit or produces audible output.
    let audio=drum_with_cavity_loss(16,super::super::acoustics::MECHANICAL_DT,true,true,
        None,false,stroke(),true,None,None,None,&[],50.0).unwrap();
    assert!(audio.acoustics.is_some()); assert!(audio.air.is_some());
    assert!(audio.force.iter().all(|v| *v==0.0));
}
