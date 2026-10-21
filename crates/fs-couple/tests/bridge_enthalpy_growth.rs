//! Original thermal events, bridge birth geometry and nonlinear feedback together.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::{LayerThermalBinding, StripEnthalpyError};
use fs_couple::strip_enthalpy::bridge::{self as coupled, BridgeEnthalpyConfig, BridgeEnthalpyError, BridgeThermalBinding};
use fs_couple::strip_enthalpy::bridge::growth::{GrowingBridge, GrowthError, TipBirth};
use fs_couple::circuit_enthalpy::material::{ConductionLink,MaterialCell,MaterialNetwork};
use fs_material::state_point::{MaterialPropertySelection,THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{PropertyKey,QueryPoint};
use fs_solid::bridge::{self as solid,Bridge,Support,growth::CellBirth};
use fs_thermal::enthalpy::*;

fn controls()->BridgeEnthalpyConfig { BridgeEnthalpyConfig {feedback:config(),maximum_deflection_m:0.009} }
fn mechanics()->solid::Config {solid::Config {section:config().mechanics,maximum_deflection_m:0.009}}
fn binding(m:&Material,id:u64,axial:u64)->BridgeThermalBinding<'_> {
    BridgeThermalBinding {axial_cell:axial,section:LayerThermalBinding {
        layer:LAYERS[0],cell:CellId(id),material:&m.mechanical,caloric:&m.caloric}}
}
fn future(m:&Material,id:u64,axial:u64,len:f64)->TipBirth<'_> {
    TipBirth {binding:binding(m,id,axial),length_m:len}
}
fn arrival(id:u64,axial:u64,start:f64,len:f64,material:u64)->TimedEvent {
    let z=geometry(0).centroid_z_m;
    TimedEvent {time_s:0.2,event:Event::Birth(DepositedCell {id:CellId(id),segment_id:axial,
        body_id:1,material_id:MaterialId(material),mass_kg:DENSITY*AREA*len,specific_h_j_kg:0.0,
        geometry:CellGeometry {start_m:[start,0.0,z],end_m:[start+len,0.0,z],
            volume_m3:AREA*len,surface_area_m2:0.002},boundary:boundary()})}
}
fn seed(m:&[Material;2],later:Vec<TimedEvent>)->(DepositionSimulation<'_>,Bridge<'_>) {
    let mut first=arrival(90,20,0.0,LENGTH,0);first.time_s=0.0;
    let mut events=vec![first];events.extend(later);
    let heat=DepositionSimulation::new(events,&[
        (MaterialId(0),m[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),m[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:20,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let beam=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
        solid::Budget {max_cells:4,max_layers_per_cell:1,max_branch_states:8}).unwrap()
        .deposit_bonded_material_layer(LAYERS[0],geometry(0),&m[0].mechanical,&m[0].caloric,
            &BTreeMap::from([(20,300.0)]),&mut ||false).unwrap();
    (heat,beam)
}
fn materials()->[Material;2] { [5.0,10.0].map(|cp|material(cp,3e8,0.05,0.0)) }
fn transport_cell<'a>(m:&'a Material,context:&'a QueryPoint,key:&'a PropertyKey,id:u64,material:u64)->MaterialCell<'a> {
    MaterialCell {cell:CellId(id),material_id:MaterialId(material),curve:m.caloric.curve(),
        point:ThermalMaterialPoint {material:&m.card,context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))}
}

#[test]
fn real_tip_birth_matches_original_owners_then_changes_the_next_coupled_sag() {
    let m=materials();
    let events=vec![arrival(99,3,LENGTH,0.02,1),TimedEvent {time_s:0.2,event:Event::Contact(
        ThermalLink {id:5,a:CellId(90),b:CellId(99),area_m2:AREA,conductance_w_k:9.0})}];
    let (heat,beam)=seed(&m,events);let map=vec![binding(&m[0],90,20)];
    let (old,initial)=GrowingBridge::try_new(heat.clone(),beam.clone(),1,map.clone(),
        vec![future(&m[1],99,3,0.02)],2,mechanics(),1e-9,||false).unwrap();
    assert_eq!(initial.thermal.events,0);
    let context=point(300.0);let key=PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=[transport_cell(&m[0],&context,&key,90,0)];
    let network=MaterialNetwork {cells:&cells,resistors:&[],links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();
    let left=coupled::step_transport_candidate(&heat,&beam,&map,&prepared,0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    let (manual_heat,_) = left.thermal().advance_to(0.2,
        StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap();
    let born=&manual_heat.cells()[&CellId(99)];
    let (manual_beam,_) = left.bridge().append_cantilever_controlled(&[CellBirth {id:3,length_m:0.02,
        mass_kg:born.mass_kg,temperature_k:manual_heat.temperature_k(CellId(99)).unwrap(),law:m[1].mechanical.law()}],
        mechanics(),&mut ||false).unwrap();
    let (next,r)=old.step_transport_candidate(&prepared,0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    assert_eq!(next.bridge(),&manual_beam);assert_eq!(next.thermal().cells(),manual_heat.cells());
    assert_eq!(r.births.thermal.events,2);assert_eq!(r.births.thermal.substeps,0);
    assert_eq!(r.births.mechanical.len(),1);assert!(r.birth_equilibrium_pending);
    assert_eq!(next.bridge().cells()[0],left.bridge().cells()[0]);
    assert!(!r.interval.layer_energy_j.contains_key(&(3,LAYERS[0])));
    assert_eq!(next.bridge().cells()[1].section().stored_energy_j(),0.0);
    assert!(prepared.resolve(next.thermal(),||false).is_err());
    let expanded=[transport_cell(&m[0],&context,&key,90,0),transport_cell(&m[1],&context,&key,99,1)];
    let links=[ConductionLink {link:5,distance_a_m:0.05,distance_b_m:0.01}];
    let successor=MaterialNetwork {cells:&expanded,resistors:&[],links:&links};
    let continued=prepared.extend_thermal(&successor,next.thermal(),||false).unwrap();
    let (after,r2)=next.step_transport_candidate(&continued,0.4,1.0,&[0.02;2],controls(),&[0.0;2],||false).unwrap();
    assert!(r2.interval.layer_energy_j[&(3,LAYERS[0])].dissipated>0.0);
    assert_eq!(after.bridge().cells()[1].section().layers()[&LAYERS[0]].birth_time_s,0.2);
    assert!(!r2.birth_equilibrium_pending);
    assert!(after.bridge().deformed_position_m(0.12,0.0).unwrap()[1]
        <next.bridge().deformed_position_m(0.12,0.0).unwrap()[1]);
    near(after.bridge().mass_kg(),MASS+DENSITY*AREA*0.02,1e-15);
    assert!(r2.total_energy_residual_j.abs()<controls().feedback.total_energy_tolerance_j);
    assert_eq!(old.remaining_births(),1);assert_eq!(old.thermal().time_s(),0.0);
}

#[test]
fn event_order_not_cell_id_order_determines_spatial_tip_order() {
    let m=materials();let (heat,beam)=seed(&m,vec![arrival(99,10,0.1,0.02,1),arrival(2,1,0.12,0.015,0)]);
    let (old,_) = GrowingBridge::try_new(heat,beam,1,vec![binding(&m[0],90,20)],
        vec![future(&m[0],2,1,0.015),future(&m[1],99,10,0.02)],3,mechanics(),1e-9,||false).unwrap();
    let (next,r)=old.step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||false).unwrap();
    assert_eq!(next.bridge().cells().iter().map(|c|c.id()).collect::<Vec<_>>(),vec![20,10,1]);
    near(next.bridge().length_m(),0.135,1e-15);assert_eq!(r.births.mechanical.len(),2);
    for i in 0..2 {near(r.births.mechanical[0].tip_m[i],r.births.mechanical[1].previous_tip_m[i],1e-14);}
    assert_eq!(next.remaining_births(),0);
}

#[test]
fn final_invalid_birth_discards_the_interval_and_earlier_successful_appends() {
    let m=materials();
    for wrong_body in [false,true] {
        let mut bad=arrival(2,1,0.12,0.015,0);
        if let Event::Birth(c)=&mut bad.event {if wrong_body {c.body_id=2;} else {c.mass_kg*=2.0;}}
        let (heat,beam)=seed(&m,vec![arrival(99,10,0.1,0.02,1),bad]);
        let (old,_) = GrowingBridge::try_new(heat,beam,1,vec![binding(&m[0],90,20)],
            vec![future(&m[1],99,10,0.02),future(&m[0],2,1,0.015)],3,mechanics(),1e-9,||false).unwrap();
        let saved=old.clone();let calls=std::cell::Cell::new(0);
        let result=old.step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],|t|{calls.set(calls.get()+1);coefficients(t)},||false);
        assert!(calls.get()>0);assert!(matches!(result,Err(GrowthError::Birth(_))));
        assert_eq!(old.bridge(),saved.bridge());assert_eq!(old.thermal().cells(),saved.thermal().cells());
        assert_eq!(old.remaining_births(),2);assert_eq!(old.thermal().time_s(),0.0);
    }
}

#[test]
fn final_cancellation_and_checkpoint_replay_are_atomic() {
    let m=materials();let (heat,beam)=seed(&m,vec![arrival(99,3,0.1,0.02,1)]);
    let (old,_) = GrowingBridge::try_new(heat,beam,1,vec![binding(&m[0],90,20)],
        vec![future(&m[1],99,3,0.02)],2,mechanics(),1e-9,||false).unwrap();
    assert!(old.step_candidate(0.3,1.0,&[0.02],controls(),&[0.0],coefficients,||false).is_err());
    let mut polls=0;
    let (next,_) = old.step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||{polls+=1;false}).unwrap();
    let mut current=0;
    assert!(matches!(old.step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||{current+=1;current==polls}),
        Err(GrowthError::Birth(BridgeEnthalpyError::Coupled(StripEnthalpyError::Cancelled)))));
    let (replay,_) = old.clone().step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||false).unwrap();
    assert_eq!(next.bridge(),replay.bridge());assert_eq!(next.thermal().cells(),replay.thermal().cells());
    assert_eq!(old.thermal().time_s(),0.0);assert_eq!(old.remaining_births(),1);
}

#[test]
fn pending_endpoint_initialization_and_passive_birth_keep_original_incoming_energy() {
    let m=materials();let mut passive=arrival(55,6,0.1,0.01,1);
    if let Event::Birth(c)=&mut passive.event {c.body_id=2;c.specific_h_j_kg=20.0;}
    let (heat,beam)=seed(&m,vec![arrival(99,3,0.1,0.02,1),passive]);
    let map=vec![binding(&m[0],90,20)];
    let left=coupled::step_candidate(&heat,&beam,&map,0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||false).unwrap();
    let (t,b,_) = left.into_parts();
    let (next,audit) = GrowingBridge::try_new(t,b,1,map,
        vec![future(&m[1],99,3,0.02)],2,mechanics(),1e-9,||false).unwrap();
    assert_eq!(audit.thermal.events,2);assert_eq!(audit.mechanical.len(),1);
    near(audit.thermal.incoming_energy_j,DENSITY*AREA*0.01*20.0,1e-15);
    assert_eq!(next.thermal().cells().len(),3);assert_eq!(next.bridge().cells().len(),2);
    assert!(next.birth_equilibrium_pending());
    near(next.thermal().temperature_k(CellId(55)).unwrap(),302.0,1e-10);
}

#[test]
fn bounded_complete_registry_is_required_without_reinitializing_material() {
    let m=materials();let (heat,beam)=seed(&m,vec![arrival(99,3,0.1,0.02,1)]);
    let map=vec![binding(&m[0],90,20)];
    assert!(GrowingBridge::try_new(heat.clone(),beam.clone(),1,map.clone(),
        vec![future(&m[1],99,3,0.02)],1,mechanics(),1e-9,||false).is_err());
    assert!(GrowingBridge::try_new(heat.clone(),beam.clone(),2,map.clone(),vec![],1,mechanics(),1e-9,||false).is_err());
    let (unregistered,_) = GrowingBridge::try_new(heat,beam,1,map,vec![],1,mechanics(),1e-9,||false).unwrap();
    assert!(matches!(unregistered.step_candidate(0.2,1.0,&[0.02],controls(),&[0.0],coefficients,||false),Err(GrowthError::Birth(_))));
}
