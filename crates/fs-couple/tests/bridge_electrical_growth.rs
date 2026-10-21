//! Finite electrical storage survives refused deposition and advances only once.
#[path = "strip_enthalpy/common.rs"]
mod common;
#[path = "electrical_contacts/source.rs"]
mod contact_source;
#[path = "bridge_electrical_growth/contacts.rs"]
mod contact_growth;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::{LayerThermalBinding, StripEnthalpyError};
use fs_couple::strip_enthalpy::bridge::{BridgeEnthalpyConfig, BridgeEnthalpyError, BridgeThermalBinding};
use fs_couple::strip_enthalpy::bridge::electrical::{self, ElectricalConfig};
use fs_couple::strip_enthalpy::bridge::growth::{GrowingBridge, TipBirth};
use fs_couple::strip_enthalpy::bridge::growth::electrical::ElectricalGrowthError;
use fs_couple::circuit_enthalpy::ResistorCell;
use fs_couple::circuit_enthalpy::material::{CellResistor, ConductionLink, MaterialCell, MaterialNetwork, PreparedMaterialNetwork};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::{PropertyKey, QueryPoint};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solid::bridge::{self as solid, Bridge, Support, growth::CellBirth};
use fs_thermal::enthalpy::*;

const Q0: f64 = 0.1;
const CAP: f64 = 0.01;
fn controls() -> ElectricalConfig {
    ElectricalConfig {bridge:BridgeEnthalpyConfig {feedback:config(),maximum_deflection_m:0.009},
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-10,power_tolerance_w:1e-9},maximum_circuit_trials:1000}
}
fn mechanics() -> solid::Config { solid::Config {section:config().mechanics,maximum_deflection_m:0.009} }
fn binding(m:&Material,id:u64,axial:u64) -> BridgeThermalBinding<'_> {
    BridgeThermalBinding {axial_cell:axial,section:LayerThermalBinding {layer:LAYERS[0],cell:CellId(id),
        material:&m.mechanical,caloric:&m.caloric}}
}
fn birth(m:&Material,id:u64,axial:u64,length:f64) -> TipBirth<'_> {
    TipBirth {binding:binding(m,id,axial),length_m:length}
}
fn arrival(id:u64,axial:u64,start:f64,length:f64,material:u64) -> TimedEvent {
    TimedEvent {time_s:0.2,event:Event::Birth(DepositedCell {id:CellId(id),segment_id:axial,body_id:1,
        material_id:MaterialId(material),mass_kg:DENSITY*AREA*length,specific_h_j_kg:0.0,
        geometry:CellGeometry {start_m:[start,0.0,geometry(0).centroid_z_m],
            end_m:[start+length,0.0,geometry(0).centroid_z_m],volume_m3:AREA*length,surface_area_m2:0.002},
        boundary:boundary()})}
}
fn seed<'a>(m:&'a [Material;2],later:Vec<TimedEvent>,future:Vec<TipBirth<'a>>) -> GrowingBridge<'a,'a> {
    let mut first=arrival(90,20,0.0,LENGTH,0); first.time_s=0.0;
    let mut events=vec![first]; events.extend(later);
    let heat=DepositionSimulation::new(events,&[(MaterialId(0),m[0].caloric.curve() as &dyn ThermalLaw),
        (MaterialId(1),m[1].caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let beam=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
        solid::Budget {max_cells:4,max_layers_per_cell:1,max_branch_states:8}).unwrap()
        .deposit_bonded_material_layer(LAYERS[0],geometry(0),&m[0].mechanical,&m[0].caloric,
            &BTreeMap::from([(20,300.0)]),&mut ||false).unwrap();
    GrowingBridge::try_new(heat,beam,1,vec![binding(&m[0],90,20)],future,4,mechanics(),1e-9,||false).unwrap().0
}
fn cells<'a>(m:&'a [Material;2],context:&'a QueryPoint,key:&'a PropertyKey,ids:&[(u64,usize)]) -> Vec<MaterialCell<'a>> {
    ids.iter().map(|&(id,i)|MaterialCell {cell:CellId(id),material_id:MaterialId(i as u64),curve:m[i].caloric.curve(),
        point:ThermalMaterialPoint {material:&m[i].card,context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((key,MaterialPropertySelection::SingleClaimOnly))}).collect()
}
fn resistors() -> [CellResistor;1] {
    [CellResistor {destination:ResistorCell {branch:1,cell:CellId(90)},selection:MaterialPropertySelection::SingleClaimOnly,
        length_m:LENGTH,area_m2:AREA}]
}
fn circuit(p:&PreparedMaterialNetwork<'_,'_>,t:&DepositionSimulation<'_>) -> MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![(1,0,Branch::Capacitor {farads:CAP}),
        (1,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},p.resolve(t,||false).unwrap().resistors,
        vec![],vec![Q0],&[]).unwrap()
}
fn materials() -> [Material;2] { [5.0,10.0].map(|cp|material(cp,3e8,0.05,0.0)) }
fn key() -> PropertyKey { PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS) }

#[test]
fn electrical_interval_and_real_birth_match_original_owners_then_continue_with_finite_supply() {
    let m=materials();
    let old=seed(&m,vec![arrival(99,3,LENGTH,0.02,1),TimedEvent {time_s:0.2,event:Event::Contact(
        ThermalLink {id:5,a:CellId(90),b:CellId(99),area_m2:AREA,conductance_w_k:9.0})}],vec![birth(&m[1],99,3,0.02)]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]); let resistors=resistors();
    let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let p=network.prepare(old.thermal(),||false).unwrap(); let circuit=circuit(&p,old.thermal());
    let left=electrical::step_candidate(old.thermal(),old.bridge(),old.bindings(),&circuit,&p,&[],
        0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    let manual_heat=left.thermal().advance_to(0.2,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let manual_beam=left.bridge().append_cantilever_controlled(&[CellBirth {id:3,length_m:0.02,
        mass_kg:manual_heat.cells()[&CellId(99)].mass_kg,temperature_k:manual_heat.temperature_k(CellId(99)).unwrap(),
        law:m[1].mechanical.law()}],mechanics(),&mut ||false).unwrap().0;
    let (born,next_circuit,r)=old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    assert_eq!(born.thermal().cells(),manual_heat.cells()); assert_eq!(born.bridge(),&manual_beam);
    assert_eq!(&next_circuit,left.circuit()); assert!(next_circuit.charges_c()[0]<Q0);
    assert_eq!(r.births.thermal.events,2); assert_eq!(r.births.thermal.substeps,0);
    assert_eq!(r.births.mechanical.len(),1); assert!(r.birth_equilibrium_pending);
    assert!(!r.interval.layer_energy_j.contains_key(&(3,LAYERS[0])));
    assert_eq!(born.bridge().cells()[1].section().stored_energy_j(),0.0);
    assert!(p.resolve(born.thermal(),||false).is_err());
    let expanded=cells(&m,&context,&key,&[(90,0),(99,1)]);
    let links=[ConductionLink {link:5,distance_a_m:0.05,distance_b_m:0.01}];
    let successor=MaterialNetwork {cells:&expanded,resistors:&resistors,links:&links};
    // Explicit new selection epoch; retain the returned circuit, never rebuild its charge.
    let next_p=successor.prepare(born.thermal(),||false).unwrap();
    let (after,final_circuit,r2)=born.step_electrical_candidate(&next_circuit,&next_p,&[],0.4,1.0,&[0.02;2],
        controls(),&[0.0;2],||false).unwrap();
    assert!(final_circuit.charges_c()[0]<next_circuit.charges_c()[0]);
    assert_eq!(final_circuit.resistors().len(),1,"new material does not imply a new circuit branch");
    assert!(r2.interval.layer_energy_j[&(3,LAYERS[0])].dissipated>0.0);
    assert_eq!(after.bridge().cells()[1].section().layers()[&LAYERS[0]].birth_time_s,0.2);
    assert!(!r2.birth_equilibrium_pending);
    let total=after.thermal().energy_j()+after.bridge().stored_energy_j()+r2.interval.electrical.stored_after_j
        +r.interval.mechanical.energy_j.numerical_dissipation+r2.interval.mechanical.energy_j.numerical_dissipation
        -r.interval.mechanical.load_work_j-r2.interval.mechanical.load_work_j;
    near(total,Q0*Q0/(2.0*CAP),2.0*config().total_energy_tolerance_j);
    assert_eq!(circuit.charges_c(),&[Q0]); assert_eq!(old.thermal().time_s(),0.0);
}

#[test]
fn final_bad_birth_rolls_back_the_electrical_debit_and_earlier_valid_append() {
    let m=materials(); let mut bad=arrival(2,1,0.12,0.015,0);
    if let Event::Birth(c)=&mut bad.event { c.mass_kg*=2.0; }
    let old=seed(&m,vec![arrival(99,3,LENGTH,0.02,1),bad],vec![birth(&m[1],99,3,0.02),birth(&m[0],2,1,0.015)]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]); let resistors=resistors();
    let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let p=network.prepare(old.thermal(),||false).unwrap(); let circuit=circuit(&p,old.thermal());
    let left=electrical::step_candidate(old.thermal(),old.bridge(),old.bindings(),&circuit,&p,&[],
        0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    assert!(left.circuit().charges_c()[0]<Q0,"the interval alone really consumes stored energy");
    let before=old.clone(); let original=circuit.clone();
    assert!(matches!(old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false),
        Err(ElectricalGrowthError::Birth(_))));
    assert_eq!(circuit,original); assert_eq!(old.bridge(),before.bridge());
    assert_eq!(old.thermal().cells(),before.thermal().cells()); assert_eq!(old.remaining_births(),2);
}

#[test]
fn final_cancellation_preserves_charge_and_replays_the_complete_growth_transaction() {
    let m=materials(); let old=seed(&m,vec![arrival(99,3,LENGTH,0.02,1)],vec![birth(&m[1],99,3,0.02)]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]); let resistors=resistors();
    let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let p=network.prepare(old.thermal(),||false).unwrap(); let circuit=circuit(&p,old.thermal());
    assert!(old.step_electrical_candidate(&circuit,&p,&[],0.3,1.0,&[0.02],controls(),&[0.0],||false).is_err());
    let mut polls=0;
    let (next,c,_) = old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||{actual+=1;actual==polls}),
        Err(ElectricalGrowthError::Birth(BridgeEnthalpyError::Coupled(StripEnthalpyError::Cancelled)))));
    let (replay,replayed_c,_) = old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    assert_eq!(next.bridge(),replay.bridge()); assert_eq!(next.thermal().cells(),replay.thermal().cells());
    assert_eq!(c,replayed_c); assert_eq!(circuit.charges_c(),&[Q0]); assert_eq!(old.remaining_births(),1);
}

#[test]
fn passive_incoming_enthalpy_is_not_charged_to_the_finite_electrical_supply() {
    let m=materials(); let mut passive=arrival(55,6,LENGTH,0.01,1);
    if let Event::Birth(c)=&mut passive.event {c.body_id=2;c.specific_h_j_kg=20.0;}
    let old=seed(&m,vec![passive],vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]); let resistors=resistors();
    let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let p=network.prepare(old.thermal(),||false).unwrap(); let circuit=circuit(&p,old.thermal());
    let left=electrical::step_candidate(old.thermal(),old.bridge(),old.bindings(),&circuit,&p,&[],
        0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    let (next,c,r)=old.step_electrical_candidate(&circuit,&p,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
    assert_eq!(&c,left.circuit()); assert_eq!(next.bridge(),left.bridge());
    assert!(r.births.mechanical.is_empty()); assert!(!r.birth_equilibrium_pending);
    near(r.births.thermal.incoming_energy_j,DENSITY*AREA*0.01*20.0,1e-15);
    assert_eq!(next.thermal().cells().len(),2); assert_eq!(next.bridge().cells().len(),1);
    assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
}
