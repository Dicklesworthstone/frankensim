//! Original radiative, interface, charge, Maxwell and beam owners in one solve.
#[path="bridge_heat_network/common.rs"]
mod common;
#[path="bridge_heat_network/continuation.rs"]
mod continuation;
use common::*;
use fs_couple::strip_enthalpy::bridge::electrical::*;

#[test]
fn endpoint_heat_sources_and_beam_replay_the_original_owners() {
    with_case(0.7,0.5,true,false,vec![],|m,cx,_,network,t,b,c| {
        let prepared=network.prepare(cx,t).unwrap();
        let next=step_heat_network_candidate(cx,t,b,&map(m,false),c,&prepared,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let exact=prepared.resolve(cx,next.thermal()).unwrap();
        let r=next.report();
        let interface=&r.material.thermal_contacts[0];
        near(interface.conductance_w_k,
            1.0/(1.0/interface.bulk_conductance_w_k+interface.interface.value_m2_k_per_w()/0.001),1e-12);
        let rad=r.material.radiation.as_ref().unwrap();
        assert!(rad.sources_w()[&CELLS[0]]<0.0 && rad.sources_w()[&CELLS[1]]>0.0);
        assert_eq!(rad.sources_w(),exact.bulk.source.radiation.as_ref().unwrap().sources_w());
        for surface in rad.surfaces() {
            near(surface.phase.temperature_k(),next.thermal().temperature_k(surface.cell).unwrap(),1e-10);
        }
        let original=c.step_with_contacts_candidate(exact.bulk.resistors,exact.contacts,&[],DT,
            controls().circuit,||false).unwrap();
        assert_eq!(next.circuit(),original.circuit());
        let temperatures=(0..2).map(|i|((20,LAYERS[i]),next.thermal().temperature_k(CELLS[i]).unwrap())).collect();
        let original=b.advance_controlled(DT,&temperatures,1.0,&[0.02],mechanics(),&mut ||false).unwrap();
        assert_eq!(next.bridge(),&original.0);
        near(r.other_external_heat_j,0.0,1e-12);
        near(r.closed_radiation_net_j,DT*rad.net_heat_w(),1e-15);
        near(r.thermal.external_heat_j,r.deposited_joule_j+r.deposited_mechanical_heat_j
            +r.closed_radiation_net_j,controls().bridge.feedback.total_energy_tolerance_j);
        near(next.thermal().energy_j()-t.energy_j()+next.bridge().stored_energy_j()-b.stored_energy_j()
            +r.electrical.stored_after_j-r.electrical.stored_before_j+r.mechanical.energy_j.numerical_dissipation
            -r.mechanical.load_work_j-r.electrical.native.supplied,0.0,config().total_energy_tolerance_j);
    });
}

#[test]
fn thermal_interface_or_emissivity_only_changes_reach_current_temperature_and_sag() {
    let run=|e,rt|with_case(e,rt,true,false,vec![],|m,cx,_,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();
        let next=step_heat_network_candidate(cx,t,b,&map(m,false),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        [next.thermal().temperature_k(CELLS[0]).unwrap(),next.circuit().charges_c()[0],
            next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],next.bridge().mass_kg()]
    });
    let baseline=run(0.2,0.1);
    for changed in [run(0.9,0.1),run(0.2,2.0)] {
        assert!((baseline[0]-changed[0]).abs()>1e-5);
        assert!((baseline[1]-changed[1]).abs()>1e-10);
        assert!((baseline[2]-changed[2]).abs()>1e-11);
        assert_eq!(baseline[3],changed[3]);
    }
}

#[test]
fn absent_extra_transfer_paths_keep_contact_only_arithmetic() {
    with_case(0.7,0.5,false,false,vec![],|m,cx,_,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();let contact=n.electrical.prepare(t,||false).unwrap();
        let next=step_heat_network_candidate(cx,t,b,&map(m,false),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let original=step_with_contacts_candidate(t,b,&map(m,false),c,&contact,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        assert_eq!(next.circuit(),original.circuit());assert_eq!(next.bridge(),original.bridge());
        assert_eq!(next.thermal().cells(),original.thermal().cells());
        assert_eq!(next.report().other_external_heat_j,original.report().other_external_heat_j);
        assert_eq!(next.report().closed_radiation_net_j,0.0);
        assert!(next.report().material.radiation.is_none());
        assert!(next.report().material.thermal_contacts.is_empty());
    });
}

#[test]
fn budget_and_scope_cancellation_preserve_all_accepted_histories() {
    with_case(0.7,0.5,true,false,vec![],|m,cx,gate,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();let before_t=t.cells().clone();let before_b=b.clone();let before_c=c.clone();
        let good=step_heat_network_candidate(cx,t,b,&map(m,false),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let mut limited=controls();limited.bridge.feedback.maximum_structural_trials=good.report().structural_trials-1;
        assert!(step_heat_network_candidate(cx,t,b,&map(m,false),c,&p,&[],DT,
            1.0,&[0.02],limited,&[0.0;2]).is_err());
        gate.request();
        assert!(matches!(step_heat_network_candidate(cx,t,b,&map(m,false),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]),Err(ElectricalError::Cancelled)));
        assert_eq!(t.cells(),&before_t);assert_eq!(b,&before_b);assert_eq!(c,&before_c);
    });
}
