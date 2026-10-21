//! G0/G1/G3: original charge, Maxwell, beam and enthalpy owners in one residual.
#[path = "strip_enthalpy/common.rs"]
mod common;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::bridge::{BridgeEnthalpyConfig, BridgeThermalBinding};
use fs_couple::strip_enthalpy::bridge::electrical::*;
use fs_couple::strip_enthalpy::LayerThermalBinding;
use fs_couple::circuit_enthalpy::{self, CircuitEnthalpyInput, CircuitEnthalpySolve, ResistorCell};
use fs_couple::circuit_enthalpy::material::{CellResistor, MaterialCell, MaterialNetwork, PreparedMaterialNetwork};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_matdb::QueryPoint;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solid::bridge::{self as solid, Bridge, Support};
use fs_thermal::enthalpy::*;

const CAP: f64 = 0.01;
const CHARGE: f64 = 0.1;
const ID: CellId = CellId(90);
fn controls()->ElectricalConfig {
    ElectricalConfig {bridge:BridgeEnthalpyConfig {feedback:config(),maximum_deflection_m:0.009},
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-10,power_tolerance_w:1e-9},
        maximum_circuit_trials:1000}
}
fn material_r(branch:f64, resistivity:Vec<(f64,f64)>)->Material {
    material_with_transport(5.0,branch,0.05,0.0,vec![(280.0,0.002),(320.0,0.002)],None,resistivity)
}
fn seed(m:&Material,heater:f64)->(DepositionSimulation<'_>,Bridge<'_>) {
    let mut c=cell(ID,0); c.segment_id=20; c.boundary.source_w=heater;
    let t=DepositionSimulation::new(vec![TimedEvent {time_s:0.0,event:Event::Birth(c)}],
        &[(MaterialId(0),m.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {max_cells:8,max_links:8,max_events:16,max_substeps:100}).unwrap()
        .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
    let b=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
        solid::Budget {max_cells:1,max_layers_per_cell:1,max_branch_states:2}).unwrap()
        .deposit_bonded_material_layer(LAYERS[0],geometry(0),&m.mechanical,&m.caloric,
            &BTreeMap::from([(20,300.0)]),&mut ||false).unwrap();
    (t,b)
}
fn bindings(m:&Material)->[BridgeThermalBinding<'_>;1] {
    [BridgeThermalBinding {axial_cell:20,section:LayerThermalBinding {layer:LAYERS[0],
        cell:ID,material:&m.mechanical,caloric:&m.caloric}}]
}
fn cells<'a>(m:&'a Material,context:&'a QueryPoint)->[MaterialCell<'a>;1] {
    // Deliberately electrical-only: source compatibility must not depend on k.
    [MaterialCell {cell:ID,material_id:MaterialId(0),curve:m.caloric.curve(),
        point:ThermalMaterialPoint {material:&m.card,context,temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:None}]
}
fn resistors()->[CellResistor;1] {
    [CellResistor {destination:ResistorCell {branch:1,cell:ID},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}]
}
fn circuit(prepared:&PreparedMaterialNetwork<'_,'_>,thermal:&DepositionSimulation<'_>,charge:f64)->MaterialCircuit {
    MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:CAP}),(1,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        prepared.resolve(thermal,||false).unwrap().resistors,vec![],vec![charge],&[]).unwrap()
}

#[test]
fn finite_capacitor_heat_matches_midpoint_rc_without_thermalizing_numerical_beam_loss() {
    let m=material_r(0.0,vec![(280.0,0.1),(320.0,0.1)]);let (heat,beam)=seed(&m,0.0);
    let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();
    assert_eq!(prepared.resistor_destinations(),vec![ResistorCell {branch:1,cell:ID}]);
    for sign in [-1.0,1.0] {
        let old=circuit(&prepared,&heat,sign*CHARGE);
        let next=step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
        let eta=0.2/(2.0*1000.0*CAP);let expected_q=sign*CHARGE*(1.0-eta)/(1.0+eta);
        let expected_heat=(CHARGE*CHARGE-expected_q*expected_q)/(2.0*CAP);
        near(next.circuit().charges_c()[0],expected_q,1e-12);
        near(next.thermal().temperature_k(ID).unwrap(),300.0+expected_heat/(MASS*5.0),1e-7);
        near(next.report().deposited_joule_j,expected_heat,1e-10);
        assert_eq!(next.report().deposited_mechanical_heat_j,0.0);
        assert!(next.report().mechanical.energy_j.numerical_dissipation>0.0);
        assert!(next.report().total_energy_residual_j.abs()<config().total_energy_tolerance_j);
        assert_eq!(old.charges_c(),&[sign*CHARGE]);assert_eq!(heat.time_s(),0.0);
    }
}

#[test]
fn electrical_resistivity_only_substitution_changes_current_temperature_and_sag() {
    let run=|rho| {
        let m=material_r(3e8,vec![(280.0,rho),(320.0,rho)]);let (heat,beam)=seed(&m,0.0);
        let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
        let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
        let prepared=network.prepare(&heat,||false).unwrap();let old=circuit(&prepared,&heat,CHARGE);
        let next=step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).unwrap();
        assert!(next.report().deposited_mechanical_heat_j>0.0);
        (next.circuit().charges_c()[0],next.thermal().temperature_k(ID).unwrap(),
            next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1])
    };
    let low=run(0.1);let high=run(0.2);
    assert!(low.0<high.0 && low.1>high.1);
    assert!((low.2-high.2).abs()>1e-8,"electrical data reaches spatial mechanics through actual heat");
}

#[test]
fn endpoint_resistance_feedback_replays_original_owners_and_preserves_finite_supply() {
    let m=material_r(3e8,vec![(280.0,0.1),(300.0,0.1),(320.0,0.5)]);
    let (mut heat,mut beam)=seed(&m,0.001);
    let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();let mut old=circuit(&prepared,&heat,CHARGE);
    let mut work=0.0;let mut other=0.0;let mut numerical=0.0;
    for i in 1..=4 {
        let g=if i<=2 {1.0}else{0.0};let loads=[if i<=2 {0.02}else{0.0}];
        let run=||step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],i as f64*0.2,
            g,&loads,controls(),&[0.0],||false).unwrap();
        let next=run();let replay=run();
        assert_eq!(next.circuit(),replay.circuit());assert_eq!(next.bridge(),replay.bridge());
        assert_eq!(next.thermal().cells(),replay.thermal().cells());
        let t=next.thermal().temperature_k(ID).unwrap();let r=next.report();
        near(r.electrical.resistors[0].material.conductor.resistance_ohm(),(0.1+0.02*(t-300.0))*LENGTH/AREA,1e-7);
        let actual=prepared.resolve(next.thermal(),||false).unwrap();
        let electrical=old.step_candidate(actual.resistors,&[],i as f64*0.2-heat.time_s(),controls().circuit,||false).unwrap();
        assert_eq!(next.circuit(),electrical.circuit());
        let manual=beam.advance_controlled(i as f64*0.2,&BTreeMap::from([((20,LAYERS[0]),t)]),g,&loads,
            solid::Config {section:config().mechanics,maximum_deflection_m:0.009},&mut ||false).unwrap();
        assert_eq!(next.bridge(),&manual.0);
        near(r.other_external_heat_j,0.0002,1e-12);
        near(r.thermal.external_heat_j,r.deposited_joule_j+r.deposited_mechanical_heat_j+r.other_external_heat_j,1e-10);
        work+=r.mechanical.load_work_j;other+=r.other_external_heat_j;numerical+=r.mechanical.energy_j.numerical_dissipation;
        near(next.bridge().stored_energy_j()+next.thermal().energy_j()+r.electrical.stored_after_j
            +numerical-work-other,CHARGE*CHARGE/(2.0*CAP),i as f64*config().total_energy_tolerance_j);
        let (t,b,c,_) = next.into_parts();heat=t;beam=b;old=c;
    }
    assert!(beam.cells()[0].section().layers()[&LAYERS[0]].material_points()[0].branch_strains()[0].abs()>1e-10);
}

#[test]
fn unloaded_lossless_bridge_matches_the_existing_circuit_thermal_solver() {
    let m=material_r(0.0,vec![(280.0,0.1),(320.0,0.1)]);let (heat,beam)=seed(&m,0.0);
    let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();let old=circuit(&prepared,&heat,CHARGE);
    let next=step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,0.0,&[0.0],controls(),&[0.0],||false).unwrap();
    let c=config();let routes=prepared.resistor_destinations();
    let manual=circuit_enthalpy::step_candidate(&heat,0.2,CircuitEnthalpyInput {circuit:&old,ports:&[],destinations:&routes},
        CircuitEnthalpySolve {maximum_cells:c.maximum_cells,maximum_circuit_trials:1000,
            maximum_newton_iterations:c.maximum_newton_iterations,difference_step_j:c.difference_step_j,
            cell_energy_tolerance_j:c.cell_energy_tolerance_j,thermal_energy_tolerance_j:c.thermal_energy_tolerance_j,
            total_energy_tolerance_j:c.total_energy_tolerance_j,circuit:controls().circuit,newton:c.newton},
        &[0.0],|t|prepared.resolve(t,||false),||false).unwrap();
    assert_eq!(next.circuit(),manual.circuit());assert_eq!(next.thermal().cells(),manual.thermal().cells());
    assert_eq!(next.bridge().stored_energy_j(),0.0);
}

#[test]
fn electrical_only_context_mismatch_and_expired_source_refuse_all_histories() {
    for narrow in [false,true] {
        let m=material_r(3e8,vec![(280.0,0.1),(if narrow {300.02}else{320.0},0.1)]);
        let (heat,beam)=seed(&m,0.0);
        let context=if narrow {point(300.0)}else{point(300.0).with("lot",12.0).unwrap()};
        let cells=cells(&m,&context);let resistors=resistors();
        let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
        let prepared=network.prepare(&heat,||false).unwrap();let old=circuit(&prepared,&heat,CHARGE);
        let saved=old.clone();let h=heat.cells().clone();let b=beam.clone();
        assert!(step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],controls(),&[0.0],||false).is_err());
        assert_eq!(old,saved);assert_eq!(heat.cells(),&h);assert_eq!(beam,b);
    }
}

#[test]
fn final_cancellation_and_structural_budget_leave_three_owners_unchanged() {
    let m=material_r(3e8,vec![(280.0,0.1),(320.0,0.1)]);let (heat,beam)=seed(&m,0.0);
    let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();let old=circuit(&prepared,&heat,CHARGE);
    let mut polls=0;
    let good=step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],controls(),&[0.0],||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],controls(),&[0.0],
        ||{actual+=1;actual==polls}),Err(ElectricalError::Cancelled)));
    let mut limited=controls();limited.bridge.feedback.maximum_structural_trials=good.report().structural_trials-1;
    assert!(step_candidate(&heat,&beam,&bindings(&m),&old,&prepared,&[],0.2,1.0,&[0.02],limited,&[0.0],||false).is_err());
    assert_eq!(old.charges_c(),&[CHARGE]);assert_eq!(heat.time_s(),0.0);assert_eq!(beam.time_s(),0.0);
    assert_eq!(beam.stored_energy_j(),0.0);
}
