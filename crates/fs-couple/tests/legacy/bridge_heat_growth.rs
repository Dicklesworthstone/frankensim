//! Radiation/contact heat and scheduled spatial births publish together.
#[path="bridge_heat_network/common.rs"]
mod common;
#[path="bridge_heat_growth/continuation.rs"]
mod continuation;
use common::*;
use fs_couple::strip_enthalpy::bridge::{BridgeThermalBinding, electrical};
use fs_couple::strip_enthalpy::bridge::growth::{GrowingBridge, TipBirth};
use fs_couple::strip_enthalpy::bridge::growth::electrical::ElectricalGrowthError;
use fs_couple::strip_enthalpy::LayerThermalBinding;
use fs_solid::bridge::growth::CellBirth;
use fs_thermal::enthalpy::{CellId, Event, StepConfig};

fn future(m:&Material,thermal:u64,axial:u64,length:f64)->TipBirth<'_> {
    TipBirth {binding:BridgeThermalBinding {axial_cell:axial,section:LayerThermalBinding {
        layer:LAYERS[0],cell:CellId(thermal),material:&m.mechanical,caloric:&m.caloric}},length_m:length}
}

#[test]
fn full_heat_interval_and_birth_replay_original_owners_and_retain_sources() {
    with_case(0.7,0.5,true,true,vec![arrival(77,99,LENGTH,0.02,1)],|m,cx,_,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();
        let (old,_)=GrowingBridge::try_new(t.clone(),b.clone(),1,map(m,true),
            vec![future(&m[1],77,99,0.02)],2,mechanics(),1e-9,||false).unwrap();
        let left=electrical::step_heat_network_candidate(cx,t,b,old.bindings(),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let manual_t=left.thermal().advance_to(DT,StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
        let manual_b=left.bridge().append_cantilever_controlled(&[CellBirth {id:99,length_m:0.02,
            mass_kg:manual_t.cells()[&CellId(77)].mass_kg,temperature_k:manual_t.temperature_k(CellId(77)).unwrap(),
            law:m[1].mechanical.law()}],mechanics(),&mut ||false).unwrap().0;
        let (next,circuit,r)=old.step_heat_network_candidate(cx,c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        assert_eq!(next.thermal().cells(),manual_t.cells());assert_eq!(next.bridge(),&manual_b);
        assert_eq!(&circuit,left.circuit());assert_eq!(r.births.thermal.substeps,0);
        assert_eq!(r.births.mechanical.len(),1);assert!(r.birth_equilibrium_pending);
        assert_eq!(r.interval.material.radiation.as_ref().unwrap().sources_w(),
            left.report().material.radiation.as_ref().unwrap().sources_w());
        assert_eq!(r.interval.material.thermal_contacts.len(),1);
        assert_eq!(next.bridge().cells().len(),2);assert_eq!(next.thermal().cells().len(),3);
        assert_eq!(next.bridge().cells()[1].section().stored_energy_j(),0.0);
        assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
        assert!(p.resolve(cx,next.thermal()).is_err(),"a grown topology requires complete new material bindings");
        assert_eq!(old.thermal().time_s(),0.0);assert_eq!(c.charges_c(),&[CHARGE]);
    });
}

#[test]
fn late_bad_birth_rejects_charge_radiative_receiver_and_earlier_append_together() {
    let mut bad=arrival(78,101,0.12,0.015,1);
    if let Event::Birth(cell)=&mut bad.event {cell.mass_kg*=2.0;}
    with_case(0.7,0.5,true,true,vec![arrival(77,99,LENGTH,0.02,1),bad],|m,cx,_,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();
        let (old,_)=GrowingBridge::try_new(t.clone(),b.clone(),1,map(m,true),
            vec![future(&m[1],77,99,0.02),future(&m[1],78,101,0.015)],3,mechanics(),1e-9,||false).unwrap();
        let left=electrical::step_heat_network_candidate(cx,t,b,old.bindings(),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        assert!(left.circuit().charges_c()[0]<CHARGE);
        assert!(left.thermal().temperature_k(CELLS[1]).unwrap()>290.0);
        let saved=old.clone();let charge=c.clone();
        assert!(matches!(old.step_heat_network_candidate(cx,c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]),Err(ElectricalGrowthError::Birth(_))));
        assert_eq!(old.thermal().cells(),saved.thermal().cells());
        assert_eq!(old.thermal().links(),saved.thermal().links());
        assert_eq!(old.bridge(),saved.bridge());assert_eq!(c,&charge);
        assert_eq!(old.remaining_births(),2);
    });
}

#[test]
fn passive_arrival_retains_incoming_enthalpy_without_an_extra_circuit_step() {
    let mut passive=arrival(77,99,LENGTH,0.02,1);
    if let Event::Birth(cell)=&mut passive.event {cell.body_id=2;cell.specific_h_j_kg=20.0;}
    with_case(0.7,0.5,true,true,vec![passive],|m,cx,_,n,t,b,c| {
        let p=n.prepare(cx,t).unwrap();
        let (old,_)=GrowingBridge::try_new(t.clone(),b.clone(),1,map(m,true),
            vec![],1,mechanics(),1e-9,||false).unwrap();
        let left=electrical::step_heat_network_candidate(cx,t,b,old.bindings(),c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        let (next,circuit,r)=old.step_heat_network_candidate(cx,c,&p,&[],DT,
            1.0,&[0.02],controls(),&[0.0;2]).unwrap();
        assert_eq!(&circuit,left.circuit());assert_eq!(next.bridge(),left.bridge());
        assert_eq!(next.thermal().cells().len(),3);assert!(r.births.mechanical.is_empty());
        assert!(!r.birth_equilibrium_pending);
        near(r.births.thermal.incoming_energy_j,DENSITY*AREA*0.02*20.0,1e-15);
        near(next.thermal().energy_j()-left.thermal().energy_j(),r.births.thermal.incoming_energy_j,1e-12);
        assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
    });
}
