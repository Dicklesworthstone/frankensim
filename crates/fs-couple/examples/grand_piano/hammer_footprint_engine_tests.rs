//! The real felt/Prony/hammer owner, with several contacts on each string.
use super::*;
use super::super::{board,geometry};

fn instrument(span:bool,flexible:bool)->Instrument {
    let c=geometry::demonstration_scale().unwrap()[48];
    let text=if span {format!("{}\nspan,69,{},4",hammer_footprint::HEADER,0.06*c.length_m)}
        else {format!("{}\npoint,69",hammer_footprint::HEADER)};
    let geometry=hammer_footprint::Specification::read(&text,&[c]).unwrap();
    Instrument::new_with_contact_geometry(vec![c],&board::demonstration(),48_000,4,12,true,
        vec![(felt::demonstration_law().unwrap(),relaxation::demonstration_prony())],
        flexible.then(ShankGeometry::published),Some(&geometry)).unwrap()
}
fn balance(p:&Instrument){
    let defect=p.accounting.input_work_j-p.accounting.dissipated_j()-p.energy_j();
    assert!(defect.abs()<1e-7,"complete contact work defect {defect:e}");
}
fn same_history(a:&Instrument,b:&Instrument){
    assert_eq!(a.bank.q,b.bank.q);assert_eq!(a.bank.v,b.bank.v);
    for (a,b) in a.contacts.iter().zip(&b.contacts){
        assert_eq!(a.state.eps_max,b.state.eps_max);assert_eq!(a.state.sig_max,b.state.sig_max);
        assert_eq!(a.memory,b.memory);assert_eq!(a.overlap,b.overlap);assert_eq!(a.force,b.force);
        assert_eq!(a.enabled,b.enabled);
    }
    assert_eq!(a.accounting.input_work_j,b.accounting.input_work_j);
    assert_eq!(a.accounting.dissipated_j(),b.accounting.dissipated_j());
}

#[test]
fn quadrature_splits_area_and_prony_stiffness_without_multiplying_the_hammer(){
    let p=instrument(true,false);let original=instrument(false,false);
    assert_eq!(p.hammers.len(),1);assert_eq!(p.bank.q,original.bank.q);
    assert_eq!(p.hammer_models[0].compliance(),original.hammer_models[0].compliance());
    let total:f64=p.contact_areas.iter().sum();
    assert!((total-p.courses[0].felt_area_m2).abs()<1e-18);
    let reference_compliance=original.creep[0].compliance();
    for i in 0..p.contacts.len(){
        let fraction=p.bank.contact_area_fraction(i);
        assert!((p.creep[i].compliance()*fraction/reference_compliance-1.0).abs()<1e-14);
        // Equal stress gives equal creep motion and area-proportional energy.
        let (site,loss)=p.creep[i].advance(&relaxation::Memory::default(),fraction*3.0);
        let (whole,reference_loss)=original.creep[0].advance(&relaxation::Memory::default(),3.0);
        assert!((p.creep[i].deformation(&site)-original.creep[0].deformation(&whole)).abs()<1e-15);
        assert!((loss-fraction*reference_loss).abs()<1e-15);
    }
}

#[test]
fn actual_finite_contact_keeps_distinct_felt_memory_and_changes_board_motion(){
    let mut spatial=instrument(true,false);let mut point=instrument(false,false);
    spatial.note_on(69,2.0).unwrap();point.note_on(69,2.0).unwrap();
    let mut changed=0.0_f64;let mut peak=0.0_f64;
    for _ in 0..2400{
        peak=peak.max(spatial.step().unwrap().abs());point.step().unwrap();
        for (a,b) in spatial.bank.q.iter().zip(&point.bank.q){changed=changed.max((a-b).abs());}
    }
    assert!(peak>1e-10 && changed>1e-12);
    assert!(spatial.contacts.iter().all(|c|c.state.eps_max>0.0));
    assert!(spatial.contacts[..4].windows(2).any(|c|c[0].state.eps_max!=c[1].state.eps_max));
    assert!(spatial.accounting.felt_relaxation_loss_j>0.0);balance(&spatial);balance(&point);
}

#[test]
fn una_corda_selects_whole_strings_not_a_subset_of_one_strings_sites(){
    let mut p=instrument(true,false);assert_eq!(p.courses[0].unison,3);
    p.set_una_corda(true);p.note_on(69,1.0).unwrap();
    assert_eq!(p.contacts.iter().filter(|c|c.enabled).count(),8);
    for sites in p.contacts.chunks_exact(4){assert!(sites.iter().all(|c|c.enabled==sites[0].enabled));}
    let active_area:f64=p.contacts.iter().zip(&p.contact_areas).filter(|(c,_)|c.enabled).map(|(_,a)|a).sum();
    assert!((active_area/p.courses[0].felt_area_m2-2.0/3.0).abs()<1e-15);
    p.set_sostenuto(true);p.note_off(69).unwrap();assert!(p.hammers[0].latched);
    for _ in 0..1200{p.step().unwrap();}
    assert!(p.contacts[..8].iter().any(|c|c.state.eps_max>0.0));
    assert!(p.contacts[8..].iter().all(|c|c.state.eps_max==0.0));
    assert_eq!(p.accounting.damper_loss_j,0.0);
    p.set_sostenuto(false);for _ in 0..128{p.step().unwrap();}
    assert!(p.accounting.damper_loss_j>0.0);balance(&p);
}

#[test]
fn finite_contact_refusal_and_restrike_keep_every_site_history(){
    let mut a=instrument(true,false);let mut b=instrument(true,false);
    // Compose with the concurrently added passive multiport radiation owner;
    // this declared pole is a runtime fixture, not a fitted acoustic specimen.
    let model=radiation::Model{ports:a.bank.board_count,poles:vec![radiation::Pole{
        omega:std::f64::consts::TAU*300.0,zeta:0.02,coupling:vec![40.0;a.bank.board_count],
    }]};
    a.configure_radiation(&model).unwrap();b.configure_radiation(&model).unwrap();
    a.note_on(69,2.0).unwrap();b.note_on(69,2.0).unwrap();
    for _ in 0..1200{
        a.step().unwrap();b.step().unwrap();
        if a.contacts.iter().any(|c|c.force>0.0){break;}
    }
    assert!(a.contacts.iter().any(|c|c.force>0.0));
    a.note_off(69).unwrap();b.note_off(69).unwrap();
    // Exercise the existing refusal after the damping/prediction prefix. No
    // production failpoint or altered physical tolerance is introduced.
    let air=a.radiation_energy_j();
    let diagonal=a.contact_h[0];a.contact_h[0]=f64::NAN;
    assert!(a.step().is_err());a.contact_h[0]=diagonal;same_history(&a,&b);
    assert_eq!(a.radiation_energy_j(),air);
    for _ in 0..8000{a.step().unwrap();b.step().unwrap();if !a.hammers[0].active{break;}}
    assert!(!a.hammers[0].active);same_history(&a,&b);
    let history:Vec<_>=a.contacts.iter().map(|c|c.state.eps_max).collect();
    assert!(history.iter().any(|e|*e>0.0));
    a.note_on(69,1.0).unwrap();b.note_on(69,1.0).unwrap();
    assert_eq!(history,a.contacts.iter().map(|c|c.state.eps_max).collect::<Vec<_>>());
    for _ in 0..1200{assert_eq!(a.step().unwrap(),b.step().unwrap());}
    same_history(&a,&b);balance(&a);
}

#[test]
fn point_footprints_are_bit_identical_and_finite_sites_accept_physical_jack_drive(){
    let c=geometry::demonstration_scale().unwrap()[48];
    let mut legacy=Instrument::new(vec![c],&board::demonstration(),48_000,4,12,true).unwrap();
    let mut point=instrument(false,false);legacy.note_on(69,2.0).unwrap();point.note_on(69,2.0).unwrap();
    for _ in 0..1200{assert_eq!(legacy.step().unwrap(),point.step().unwrap());}
    same_history(&legacy,&point);
    let mut flexible=instrument(true,true);flexible.jack_on(69,70.0,0.007).unwrap();
    for _ in 0..2400{flexible.step().unwrap();}
    assert!(flexible.contacts.iter().any(|c|c.state.eps_max>0.0));
    assert!(flexible.accounting.input_work_j>0.0 && flexible.accounting.shank_loss_j>0.0);
    assert_eq!(flexible.hammers[0].jack.peak_n,0.0);balance(&flexible);
}
