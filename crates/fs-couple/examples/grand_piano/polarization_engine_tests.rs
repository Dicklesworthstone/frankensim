//! The original nonlinear hammer/Prony/pedal/acoustic loop with vector strings.
use super::*;
use super::super::{board,geometry};
fn course()->Course {
    Course {unison:2,duplex_length_m:0.12,detune_cents:0.8,..geometry::demonstration_scale().unwrap()[48]}
}
fn build(side:f64,ratio:f64,spatial:bool,stretch:bool)->Instrument {
    let c=course();let modes=board::demonstration();
    let row=vec![(0..modes.len()).map(|i|side*(0.3+i as f64).cos()).collect::<Vec<_>>()];
    let faces=hammer_footprint::Specification::read("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n",&[c]).unwrap();
    let mut p=Instrument::new_with_transverse_contact_geometry(vec![c],&modes,48_000,4,12,true,
        vec![(felt::demonstration_law().unwrap(),relaxation::demonstration_prony())],
        Some(ShankGeometry::published()),Some(&faces),Some((&row,&[ratio]))).unwrap();
    if spatial {p.configure_dampers(&dampers::Specification::estimated(&[c]).unwrap()).unwrap();}
    if stretch {p.configure_string_stretching(&super::super::linear::string_stretching::Specification::read(
        "frankensim-piano-string-stretching-v1\nstretch,69,100000,0.2\n",&[c]).unwrap()).unwrap();}
    p
}
#[test]
fn actual_finite_hammer_drives_both_directions_without_doubling_area_or_launch_work() {
    let c=course();let modes=board::demonstration();
    let faces=hammer_footprint::Specification::read("frankensim-hammer-footprints-v1\nspan,69,0.008,2\n",&[c]).unwrap();
    let mut bare=Instrument::new_with_contact_geometry(vec![c],&modes,48_000,4,12,true,
        vec![(felt::demonstration_law().unwrap(),relaxation::demonstration_prony())],
        Some(ShankGeometry::published()),Some(&faces)).unwrap();
    let mut vector=build(0.06,0.3,true,true);
    assert_eq!(bare.hammers.len(),vector.hammers.len());
    assert_eq!(bare.contact_areas,vector.contact_areas);
    assert_eq!(bare.hammer_contact_count(),vector.hammer_contact_count());
    assert_eq!(vector.bank.strings.len(),2*bare.bank.strings.len());
    assert!(vector.bank.contact_strings.iter().all(|&s|vector.bank.strings[s].polarization==0));
    bare.note_on(69,1.).unwrap();vector.note_on(69,1.).unwrap();
    assert_eq!(bare.accounting.input_work_j,vector.accounting.input_work_j);
    for n in 0..2400 {
        if n==1200 {vector.note_off(69).unwrap();}
        vector.step().unwrap();
    }
    assert!(vector.accounting.felt_loss_j>0. && vector.accounting.shank_loss_j>0.);
    assert!(vector.accounting.damper_loss_j>0.);
    assert!(vector.bank.strings.iter().filter(|s|s.polarization==1).any(|s|
        vector.bank.q[s.modes.clone()].iter().any(|q|q.abs()>1e-12)));
    for (i,s) in vector.bank.strings.iter().enumerate().filter(|(_,s)|s.polarization==0) {
        let j=vector.bank.strings.iter().position(|p|p.course==s.course && p.member==s.member
            && p.duplex==s.duplex && p.polarization==1).unwrap();
        assert_eq!(vector.bank.string_stretching_observation(i).unwrap().tension_n,
            vector.bank.string_stretching_observation(j).unwrap().tension_n);
    }
    assert!((vector.accounting.input_work_j-vector.energy_j()-vector.accounting.dissipated_j()).abs()<1e-7);
}
#[test]
fn explicit_lateral_pads_follow_pedals_and_never_damp_duplex_coordinates() {
    for spatial in [false,true] {for ratio in [0.,0.25] {
        let mut p=build(0.,ratio,spatial,false);
        for s in p.bank.strings.iter().filter(|s|s.polarization==1) {
            p.bank.v[s.modes.clone()].fill(0.003);
        }
        let initial=p.bank.v.clone();let energy=p.energy_j();
        p.set_sustain(1.).unwrap();assert_eq!(p.damp(0.001).unwrap(),0.);
        assert_eq!(p.bank.v,initial);p.set_sustain(0.).unwrap();
        let loss=p.damp(0.001).unwrap();
        assert!((p.energy_j()+loss-energy).abs()<1e-12);
        if ratio==0. {assert_eq!(loss,0.);assert_eq!(p.bank.v,initial);}
        else {
            assert!(loss>0.);
            for s in &p.bank.strings {
                if s.polarization==1 && !s.duplex {
                    assert_ne!(&p.bank.v[s.modes.clone()],&initial[s.modes.clone()]);
                } else {assert_eq!(&p.bank.v[s.modes.clone()],&initial[s.modes.clone()]);}
            }
        }
    }}
}
#[test]
fn reacted_vector_string_samples_roll_back_both_directions_and_share_one_trace_clock() {
    let load=|p:&mut Instrument| {
        let r=p.bank.board_count;
        p.configure_radiation(&radiation::Model {ports:r,poles:vec![radiation::Pole {
            omega:700.,zeta:0.2,coupling:vec![100.;r],
        }]}).unwrap();p.note_on(69,1.).unwrap();
    };
    let mut a=build(0.06,0.3,false,true);let mut b=build(0.06,0.3,false,true);
    load(&mut a);load(&mut b);let mut trace=vec![0.;a.board_trace_len()];let r=a.bank.board_count;
    for _ in 0..500 {
        a.step_with_board_trace(&mut trace).unwrap();
        for sub in 0..b.substeps {
            b.mechanics_step().unwrap();
            assert_eq!(&trace[sub*r..(sub+1)*r],&b.bank.v[b.bank.modes.len()..]);
        }
    }
    let q=a.bank.q.clone();let v=a.bank.v.clone();let energy=a.energy_j();
    let air=a.radiation_energy_j();let work=a.accounting.input_work_j;
    a.damper_drag_ns_m=f64::NAN;assert!(a.step_with_board_trace(&mut trace).is_err());
    assert_eq!(a.bank.q,q);assert_eq!(a.bank.v,v);assert_eq!(a.energy_j(),energy);
    assert_eq!(a.radiation_energy_j(),air);assert_eq!(a.accounting.input_work_j,work);
    a.damper_drag_ns_m=0.4;a.step().unwrap();b.step().unwrap();
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    assert_eq!(a.radiation_energy_j(),b.radiation_energy_j());
    assert!((a.accounting.input_work_j-a.energy_j()-a.accounting.dissipated_j()).abs()<1e-7);
}
