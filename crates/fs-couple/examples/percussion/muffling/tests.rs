use super::*;
use super::super::{Stroke, drum_with_mufflers, drum_with_sticks, splash_with_mufflers,
    drum_spec, mechanics::Mechanics, snare};
use fs_exec::CancelGate;
fn args(text: &str) -> Vec<String> { text.split_whitespace().map(str::to_owned).collect() }
fn batter(resistance: f64) -> Muffler {
    Muffler { surface: Surface::Batter, position_m: [0.06,0.01], resistance_n_s_m: resistance }
}
fn stroke() -> Stroke { Stroke { speed_m_s: 4.0, position_m: Some([0.06,0.01]) } }

#[test]
fn repeated_si_controls_preserve_playing_inputs_and_refuse_wrong_surfaces() {
    let mut a = args("drum-modal 128 --muffler batter 0.06 0.01 2 --strike-speed-m-s 4 --muffler resonant -0.05 0.02 0");
    let m = options(&mut a).unwrap();
    assert_eq!(m[0],batter(2.0)); assert_eq!(m[1].surface,Surface::Resonant);
    assert_eq!(m[1].resistance_n_s_m,0.0);
    let (rest,first) = super::super::playing::parse(a).unwrap();
    assert_eq!(rest,["drum-modal","128"]); assert_eq!(first.speed_m_s,4.0);
    assert!(admit_command(&m,"snare-off-mic").is_ok());
    assert!(admit_command(&m,"splash").is_err());
    let shell = Muffler { surface:Surface::Shell,..batter(1.0) };
    assert!(admit_command(&[shell],"splash-mic").is_ok());
    assert!(admit_command(&[shell],"drum-mic").is_err());
    for text in ["--muffler batter 0 0", "--muffler rim 0 0 1", "--muffler shell NaN 0 1",
        "--muffler batter 0 0 -1", "--muffler resonant 0 0 inf"] {
        assert!(options(&mut args(text)).is_err(),"{text}");
    }
    assert!(options(&mut args(&"--muffler batter 0 0 1 ".repeat(17))).is_err());
}

#[test]
fn actual_head_rows_preserve_physical_scaling_and_all_unattached_coordinates() {
    let (films,modes) = drum_spec::Spec::reference().prepare(2e-6,false).unwrap();
    let end = 1+modes.iter().map(Vec::len).sum::<usize>();
    let specs = [batter(2.0),Muffler { surface:Surface::Resonant,..batter(0.7) }];
    let dampers = head_ports(&specs,&films,&modes,end+17).unwrap();
    for (head,port) in dampers.iter().enumerate() {
        let expected = film_shapes(&films[head],&modes[head],&[specs[head].position_m]).unwrap().remove(0);
        let start = if head==0 {1} else {1+modes[0].len()};
        assert_eq!(&port.weights[start..start+expected.len()],expected);
        assert!(port.weights[..start].iter().chain(&port.weights[start+expected.len()..]).all(|b|*b==0.0));
        assert_eq!(port.damping_n_s_m,specs[head].resistance_n_s_m);
    }
    assert!(head_ports(&[Muffler {position_m:[1.0,0.0],..batter(1.0)}],&films,&modes,end).is_err());
    assert!(head_ports(&[batter(1.0)],&films,&modes,end-1).is_err());
}

#[test]
fn a_muffler_changes_real_head_motion_while_the_existing_energy_balance_closes() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false,true] {
        let mut free = drum_with_sticks(256,2e-6,false,prepared,None,false,stroke(),false,None,None,None).unwrap();
        let mut held = drum_with_mufflers(256,2e-6,false,prepared,None,false,stroke(),false,None,None,None,&[batter(2.0)]).unwrap();
        assert_eq!(free.system.state(),held.system.state());
        assert_eq!(free.observer_a,held.observer_a); assert_eq!(free.observer_b,held.observer_b);
        let initial = 0.5*(stroke().speed_m_s/held.stick_weight).powi(2);
        let mut loss = 0.0; let mut changed = 0.0_f64;
        for _ in 0..256 {
            free.system.step(&free.force,&gate).unwrap();
            let f = held.system.step(&held.force,&gate).unwrap();
            loss += f.dissipated_energy_j;
            assert!(f.dissipated_energy_j>=0.0 && f.supplied_work_j==0.0);
            assert!((f.stored_energy_j+loss-initial).abs()<1e-6);
            for (&a,&b) in free.system.state()[2..].iter().zip(&held.system.state()[2..]) {
                changed = changed.max((a-b).abs());
            }
        }
        assert!(changed>1e-10 && loss>0.0,"actual head/contact motion must respond to the attached resistance");
    }
}

#[test]
fn mufflers_compose_with_two_sticks_snare_air_and_the_unchanged_audio_boundary() {
    let second = Some(Stroke {speed_m_s:2.5,position_m:Some([-0.05,0.02])});
    let wires = snare::SnareSet::reference(false);
    let e = drum_with_mufflers(4,super::super::acoustics::MECHANICAL_DT,true,true,
        Some(wires),false,stroke(),true,None,None,second,&[batter(0.4)]).unwrap();
    let Mechanics::Prepared(system) = &e.system else {panic!("prepared snare/cavity image")};
    let stick = e.second_stick.unwrap(); let air = e.air.as_ref().unwrap();
    assert_eq!(system.contact_count(),2+wires.strands*wires.contact_cells);
    assert_eq!(air.coupling.structural_modes(),stick.coordinate+1+wires.mode_count().unwrap());
    assert_eq!(e.force.len(),air.coupling.total_modes());
    assert!(e.observer_a[stick.coordinate..].iter().all(|b|*b==0.0));
    assert!(e.observer_b[stick.coordinate..].iter().all(|b|*b==0.0));
    assert!(e.acoustics.is_some());
    // The nonlinear cavity wrapper must keep the same structural port too.
    let mut e = drum_with_mufflers(4,2e-6,false,false,None,false,stroke(),true,None,None,
        second,&[batter(0.4)]).unwrap();
    let gate = CancelGate::new_clock_free();
    e.system = e.system.into_prepared_nonlinear().unwrap();
    for _ in 0..4 {e.system.step(&e.force,&gate).unwrap();}
}

#[test]
fn cymbal_mufflers_keep_felt_storage_and_prepare_the_same_nonlinear_shell() {
    let shell = Muffler {surface:Surface::Shell,position_m:[0.075,0.0],resistance_n_s_m:0.5};
    let mut e = splash_with_mufflers(4,2e-6,false,stroke(),None,&[shell]).unwrap();
    // Kelvin felt memory is still behind the original mechanical prefix.
    assert!(e.system.state().len()>2*e.force.len());
    e.system = e.system.into_prepared_nonlinear().unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..4 {
        let f = e.system.step(&e.force,&gate).unwrap();
        assert!(f.stored_energy_j.is_finite() && f.dissipated_energy_j>=0.0);
    }
    assert!(splash_with_mufflers(4,2e-6,false,stroke(),None,&[batter(1.0)]).is_err());
}
