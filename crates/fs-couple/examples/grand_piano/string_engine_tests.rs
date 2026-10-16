use super::*;
use super::super::{radiation,hammer_footprint,ShankGeometry};
use super::super::super::{board,geometry};

fn courses()->Vec<geometry::Course> {
    let all=geometry::demonstration_scale().unwrap();vec![all[48],all[51]]
}
fn spec(courses:&[geometry::Course],ea:Option<f64>,slope:f64)->string_stretching::Specification {
    let mut text=format!("{}\n",string_stretching::HEADER);
    for c in courses {text.push_str(&match ea {
        Some(ea)=>format!("stretch,{},{ea},{slope}\n",c.midi),None=>format!("linear,{}\n",c.midi),
    });}
    string_stretching::Specification::read(&text,courses).unwrap()
}
fn piano(ea:Option<f64>)->Instrument {
    let courses=courses();let mut faces=format!("{}\n",hammer_footprint::HEADER);
    for c in &courses {faces.push_str(&format!("span,{},0.012,4\n",c.midi));}
    let faces=hammer_footprint::Specification::read(&faces,&courses).unwrap();
    let selection=spec(&courses,ea,0.2);
    let mut p=Instrument::new_with_footprints(courses,&board::demonstration(),48000,4,12,true,&faces).unwrap();
    p.configure_string_stretching(&selection).unwrap();p
}
fn balance(p:&Instrument){assert!((p.energy_j()+p.accounting.dissipated_j()-p.accounting.input_work_j).abs()<1e-7);}
fn same(a:&Instrument,b:&Instrument) {
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
    assert_eq!(a.accounting.input_work_j,b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(),b.accounting.dissipated_j());
    for (a,b) in a.contacts.iter().zip(&b.contacts) {
        assert_eq!(a.overlap,b.overlap);assert_eq!(a.force,b.force);assert_eq!(a.memory,b.memory);
        assert_eq!(a.state.eps_max,b.state.eps_max);assert_eq!(a.state.sig_max,b.state.sig_max);
    }
    for (a,b) in a.hammers.iter().zip(&b.hammers) {
        assert_eq!(a.jack.elapsed_s,b.jack.elapsed_s);assert_eq!(a.jack.peak_n,b.jack.peak_n);
        assert_eq!(a.motion.q,b.motion.q);assert_eq!(a.motion.v,b.motion.v);
    }
}

#[test]
fn played_finite_hammer_chord_changes_actual_tension_and_board_motion_with_one_work_ledger() {
    let mut nonlinear=piano(Some(150000.));let mut linear=piano(None);
    assert_eq!(nonlinear.bank.q,linear.bank.q);
    assert_eq!(nonlinear.hammer_contact_count(),linear.hammer_contact_count());
    assert_eq!(nonlinear.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>(),
        linear.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>());
    let rest=nonlinear.bank.string_stretching_observation(0).unwrap().tension_n;
    for p in [&mut nonlinear,&mut linear] {p.note_on(69,2.).unwrap();p.note_on(72,1.4).unwrap();}
    let mut stretch=0.0_f64;let mut tension=rest;let mut motion=0.0_f64;
    for tick in 0..1800 {
        if tick==600 {for p in [&mut nonlinear,&mut linear]{p.set_sustain(0.3).unwrap();p.note_off(69).unwrap();p.note_off(72).unwrap();}}
        nonlinear.step().unwrap();linear.step().unwrap();
        for si in 0..nonlinear.bank.strings.len() {
            let o=nonlinear.bank.string_stretching_observation(si).unwrap();
            assert!(o.slope_bound<=0.2);stretch=stretch.max(o.stretching_energy_j);
        }
        tension=tension.max(nonlinear.bank.string_stretching_observation(0).unwrap().tension_n);
        motion=motion.max(nonlinear.bank.v[nonlinear.bank.modes.len()..].iter()
            .zip(&linear.bank.v[linear.bank.modes.len()..]).map(|(a,b)|(a-b).abs()).fold(0.,f64::max));
        balance(&nonlinear);
    }
    assert!(stretch>0. && tension>rest && motion>1e-10);
    assert!(nonlinear.accounting.felt_loss_j>0. && nonlinear.accounting.felt_relaxation_loss_j>0.);
    assert!(nonlinear.accounting.damper_loss_j>0.);
}

#[test]
fn nonlinear_trial_failure_restores_string_felt_pedal_and_acoustic_history_for_exact_retry() {
    let mut a=piano(Some(150000.));let mut b=piano(Some(150000.));
    let model=radiation::Model{ports:a.bank.board_count,poles:vec![radiation::Pole{
        omega:std::f64::consts::TAU*300.,zeta:0.02,coupling:vec![40.;a.bank.board_count]}]};
    for p in [&mut a,&mut b]{p.configure_radiation(&model).unwrap();p.set_una_corda(true);
        p.note_on(69,2.).unwrap();p.note_on(72,1.4).unwrap();}
    for _ in 0..300{a.step().unwrap();b.step().unwrap();}
    assert!(a.bank.string_stretching_observation(0).unwrap().stretching_energy_j>0.);
    assert!(a.radiation_energy_j()>0.);
    let n=a.contacts.len();let index=(0..n).find(|&i|a.bank.strings[a.bank.contact_strings[i]].course==1).unwrap();
    let saved=a.contact_h[index*n+index];a.contact_h[index*n+index]=f64::NAN;
    assert!(a.step().is_err());a.contact_h[index*n+index]=saved;same(&a,&b);
    for _ in 0..100{assert_eq!(a.step().unwrap(),b.step().unwrap());}same(&a,&b);balance(&a);
    let selection=spec(&a.courses,Some(150000.),0.2);
    assert!(a.configure_string_stretching(&selection).is_err());same(&a,&b);
}

#[test]
fn slope_refusal_does_not_publish_a_clipped_string_or_consume_a_sample() {
    let c=courses()[0];let mut p=Instrument::new(vec![c],&board::demonstration(),48000,4,12,true).unwrap();
    p.configure_string_stretching(&spec(&[c],Some(150000.),1e-10)).unwrap();p.note_on(69,2.).unwrap();
    let mut refused=false;
    for _ in 0..600 {
        let q=p.bank.q.clone();let v=p.bank.v.clone();let energy=p.energy_j();let loss=p.accounting.dissipated_j();
        if p.step().is_err() {
            assert_eq!(q,p.bank.q);assert_eq!(v,p.bank.v);assert_eq!(p.energy_j(),energy);
            assert_eq!(p.accounting.dissipated_j(),loss);
            assert!(p.step().is_err());assert_eq!(q,p.bank.q);assert_eq!(v,p.bank.v);
            refused=true;break;
        }
    }
    assert!(refused,"this explicit domain must reject rather than clip the struck string");
}

#[test]
fn all_linear_selection_preserves_jack_clock_and_original_material_dynamics_bitwise() {
    let c=courses()[0];let law=super::super::felt::demonstration_law().unwrap();
    let build=||Instrument::new_with_course_shanks(vec![c],&board::demonstration(),48000,4,12,true,
        vec![(law.clone(),super::super::relaxation::demonstration_prony())],ShankGeometry::published()).unwrap();
    let mut original=build();let mut selected=build();
    selected.configure_string_stretching(&spec(&[c],None,0.2)).unwrap();
    assert!(!selected.bank.has_string_stretching());
    original.jack_on(69,70.,0.007).unwrap();selected.jack_on(69,70.,0.007).unwrap();
    for _ in 0..1200{assert_eq!(original.step().unwrap(),selected.step().unwrap());}
    same(&original,&selected);balance(&selected);assert!(selected.accounting.felt_loss_j>0.);
}
