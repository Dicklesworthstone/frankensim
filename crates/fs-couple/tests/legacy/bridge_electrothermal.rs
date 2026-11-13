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
fn finite_electrical_heat_bends_bonded_layers_from_expansion_data_alone() {
    let run=|alpha:f64| {
        let slope=-0.005;
        let expansion=[alpha,0.0];
        let materials=expansion.map(|a|material_with_affine_scale(5.0,0.0,1.0,40.0*a,
            vec![(280.0,0.002),(320.0,0.002)],None,
            vec![(280.0,0.1),(320.0,0.1)],slope));
        let mut events:Vec<_>=CELLS.into_iter().enumerate().map(|(i,id)| {
            let mut c=cell(id,i);c.segment_id=20;
            TimedEvent {time_s:0.0,event:Event::Birth(c)}
        }).collect();
        events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
            id:5,a:CELLS[0],b:CELLS[1],area_m2:0.001,conductance_w_k:0.002})});
        let heat=DepositionSimulation::new(events,&[
            (MaterialId(0),materials[0].caloric.curve() as &dyn ThermalLaw),
            (MaterialId(1),materials[1].caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {max_cells:2,max_links:1,max_events:3,max_substeps:100}).unwrap()
            .advance_to(0.0,StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
        let mut beam=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
            solid::Budget {max_cells:1,max_layers_per_cell:2,max_branch_states:4}).unwrap();
        for i in 0..2 {
            beam=beam.deposit_bonded_material_layer(LAYERS[i],geometry(i),
                &materials[i].mechanical,&materials[i].caloric,
                &BTreeMap::from([(20,300.0)]),&mut ||false).unwrap();
        }
        let map=[0,1].map(|i|BridgeThermalBinding {axial_cell:20,
            section:LayerThermalBinding {layer:LAYERS[i],cell:CELLS[i],
                material:&materials[i].mechanical,caloric:&materials[i].caloric}});
        let context=point(300.0);
        let k=fs_matdb::PropertyKey::new("conductivity",fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS);
        let cells=[0,1].map(|i|MaterialCell {cell:CELLS[i],material_id:MaterialId(i as u64),
            curve:materials[i].caloric.curve(),point:ThermalMaterialPoint {
                material:&materials[i].card,context:&context,temperature_axis:"T",liquid_fraction_axis:None},
            conductivity:Some((&k,MaterialPropertySelection::SingleClaimOnly))});
        let resistors=[CellResistor {destination:ResistorCell {branch:1,cell:CELLS[0]},
            selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}];
        let links=[fs_couple::circuit_enthalpy::material::ConductionLink {
            link:5,distance_a_m:0.0005,distance_b_m:0.0005}];
        let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
        let prepared=network.prepare(&heat,||false).unwrap();
        let old=circuit(&prepared,&heat,CHARGE);
        let next=step_candidate(&heat,&beam,&map,&old,&prepared,&[],0.2,0.0,&[0.0],
            controls(),&[0.0;2],||false).unwrap();
        let t=CELLS.map(|id|next.thermal().temperature_k(id).unwrap());
        let s=t.map(|v|1.0+slope*(v-300.0));
        let e=s.map(|v|1e8*v);
        let f=[0,1].map(|i|expansion[i]*(t[i]-300.0));
        let z=[0,1].map(|i|geometry(i).centroid_z_m);
        let inertia=AREA*0.001_f64.powi(2)/12.0;
        // Independent cross-section equilibrium: integrate E*(eps-z*kappa-f)
        // and its moment through each rectangular layer, with zero applied load.
        let k00=AREA*(e[0]+e[1]);let k01=-AREA*(e[0]*z[0]+e[1]*z[1]);
        let k11=e[0]*(AREA*z[0]*z[0]+inertia)+e[1]*(AREA*z[1]*z[1]+inertia);
        let r0=AREA*(e[0]*f[0]+e[1]*f[1]);
        let r1=-AREA*(e[0]*z[0]*f[0]+e[1]*z[1]*f[1]);
        let det=k00*k11-k01*k01;
        let strain=(r0*k11-r1*k01)/det;let curvature=(k00*r1-k01*r0)/det;
        let section=next.bridge().cells()[0].section();
        near(section.axial_strain(),strain,1e-12);
        near(section.curvature_per_m(),curvature,1e-9);
        let tip=next.bridge().deformed_position_m(LENGTH,0.0).unwrap();
        near(tip[0],LENGTH*(1.0+strain),1e-12);
        near(tip[1],0.5*curvature*LENGTH*LENGTH,1e-11);
        let eta=0.2/(2.0*1000.0*CAP);let charge=CHARGE*(1.0-eta)/(1.0+eta);
        let joule=(CHARGE*CHARGE-charge*charge)/(2.0*CAP);
        near(next.circuit().charges_c()[0],charge,1e-12);
        let mut u=0.0;let mut entropy_energy=0.0;let mut psi_sum=0.0;let mut reversible=0.0;
        for i in 0..2 {
            let elastic=strain-z[i]*curvature-f[i];
            let psi=0.5*e[i]*LENGTH*(AREA*elastic*elastic+inertia*curvature*curvature);
            let integrated_stress=VOLUME*e[i]*elastic;
            let entropy=expansion[i]*integrated_stress-slope*psi/s[i];
            let rev=integrated_stress*f[i]-t[i]*entropy;
            u+=psi+t[i]*entropy;entropy_energy+=t[i]*entropy;psi_sum+=psi;reversible+=rev;
            let conduction=0.2*0.002*(t[1-i]-t[i]);
            near(MASS*5.0*(t[i]-300.0),(if i==0 {joule}else{0.0})+conduction+rev,
                config().cell_energy_tolerance_j);
        }
        let report=next.report();
        near(report.mechanical_internal_energy_change_j,u,1e-11);
        near(report.mechanical_entropy_energy_change_j,entropy_energy,1e-11);
        near(report.reversible_heat_j,reversible,1e-11);
        near(report.mechanical.energy_j.numerical_dissipation,psi_sum,1e-11);
        assert_eq!(report.mechanical.energy_j.dissipated,0.0);
        near(report.deposited_mechanical_heat_j,reversible,1e-11);
        near(report.other_external_heat_j,0.0,1e-14);
        near(u+next.thermal().energy_j()+report.electrical.stored_after_j+psi_sum,
            CHARGE*CHARGE/(2.0*CAP),config().total_energy_tolerance_j);
        assert_eq!(old.charges_c(),&[CHARGE]);assert_eq!(beam.time_s(),0.0);
        (t,tip)
    };
    let fixed=run(0.0);let expanding=run(1e-4);let contracting=run(-1e-4);
    near(fixed.1[1],0.0,1e-14);
    assert!(expanding.1[1]>0.001 && expanding.1[0]>LENGTH);
    near(expanding.1[1],-contracting.1[1],1e-10);
    near(expanding.1[0]-LENGTH,LENGTH-contracting.1[0],1e-12);
    for i in 0..2 {near(expanding.0[i],contracting.0[i],1e-7);}
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
    let slope=-0.005;
    let m=material_with_affine_scale(5.0,3e8,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],None,
        vec![(280.0,0.1),(300.0,0.1),(320.0,0.5)],slope);
    let (mut heat,mut beam)=seed(&m,0.001);
    let context=point(300.0);let cells=cells(&m,&context);let resistors=resistors();
    let network=MaterialNetwork {cells:&cells,resistors:&resistors,links:&[]};
    let prepared=network.prepare(&heat,||false).unwrap();let mut old=circuit(&prepared,&heat,CHARGE);
    let mut work=0.0;let mut other=0.0;let mut numerical=0.0;
    let mut mechanical_u=0.0;let mut entropy_energy=0.0;let mut saw_coefficient_work=false;
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
        near(r.deposited_mechanical_heat_j,r.mechanical.energy_j.dissipated+r.reversible_heat_j,1e-12);
        near(r.reversible_heat_j+r.mechanical.energy_j.modulus_work+r.mechanical_entropy_energy_change_j,0.0,1e-12);
        assert!(r.mechanical.energy_j.dissipated>=0.0);
        if i==1 {assert!(r.reversible_heat_j<0.0);}
        saw_coefficient_work|=r.mechanical.energy_j.modulus_work.abs()>1e-12;
        work+=r.mechanical.load_work_j;other+=r.other_external_heat_j;numerical+=r.mechanical.energy_j.numerical_dissipation;
        mechanical_u+=r.mechanical_internal_energy_change_j;entropy_energy+=r.mechanical_entropy_energy_change_j;
        // Use the authored affine law and retained native point histories as an
        // independent endpoint oracle, including the previously loaded branches.
        let xi_j=next.bridge().cells()[0].section().layers()[&LAYERS[0]].material_points()
            .iter().map(|p|0.5*VOLUME*p.stored_energy_j_m3()/(1.0+slope*(t-300.0))).sum::<f64>();
        near(mechanical_u,(1.0-300.0*slope)*xi_j,1e-12);
        near(entropy_energy,-slope*t*xi_j,1e-12);
        near(mechanical_u+next.thermal().energy_j()+r.electrical.stored_after_j
            +numerical-work-other,CHARGE*CHARGE/(2.0*CAP),i as f64*config().total_energy_tolerance_j);
        let (t,b,c,_) = next.into_parts();heat=t;beam=b;old=c;
    }
    assert!(saw_coefficient_work,"heating must recouple the preloaded material, not a fresh spring");
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
    let m=material_with_affine_scale(5.0,3e8,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],None,vec![(280.0,0.1),(320.0,0.1)],-0.005);
    let (heat,beam)=seed(&m,0.0);
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
