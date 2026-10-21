//! Interface-only changes reach both finite bodies and the actual beam response.
#[path = "strip_enthalpy/common.rs"]
mod common;
#[path = "electrical_contacts/source.rs"]
mod contact_source;
use common::*;
use std::collections::BTreeMap;
use fs_couple::strip_enthalpy::bridge::{BridgeEnthalpyConfig, BridgeThermalBinding};
use fs_couple::strip_enthalpy::bridge::electrical::*;
use fs_couple::strip_enthalpy::LayerThermalBinding;
use fs_couple::circuit_enthalpy::ResistorCell;
use fs_couple::circuit_enthalpy::material::{CellResistor, MaterialCell, MaterialNetwork};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_couple::circuit_enthalpy::contacts::material::{
    ElectricalContactNetwork, MaterialElectricalContact, PreparedElectricalContactNetwork,
};
use fs_couple::circuit_enthalpy::contacts::ContactHeatCells;
use fs_couple::circuit_enthalpy::material::ConductionLink;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::{interface::ThermalInterfacePoint, source::ThermalMaterialPoint};
use fs_matdb::PropertyKey;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solid::bridge::{self as solid, Bridge, Support};
use fs_thermal::enthalpy::*;

const CAP: f64 = 0.01;
const CHARGE: f64 = 0.1;
const LINK_AREA: f64 = 0.001;
const DT: f64 = 0.2;

fn controls() -> ElectricalConfig {
    ElectricalConfig {bridge:BridgeEnthalpyConfig {feedback:config(),maximum_deflection_m:0.009},
        circuit:MaterialCircuitSolve {energy_tolerance_j:1e-10,power_tolerance_w:1e-9},
        maximum_circuit_trials:1000}
}

fn map(m: &[Material;2]) -> [BridgeThermalBinding<'_>;2] {
    [0,1].map(|i|BridgeThermalBinding {axial_cell:20,section:LayerThermalBinding {
        layer:LAYERS[i],cell:CELLS[i],material:&m[i].mechanical,caloric:&m[i].caloric}})
}

// All source objects outlive every step performed by the test's callback.
// Reuse existing material, thermal and ordered-interface fixtures unchanged.
fn with_fixture<R>(branch: f64, contact_ohm: f64, slope: f64, fraction: f64, charge: f64,
    run: impl FnOnce(&[Material;2], &DepositionSimulation<'_>, &Bridge<'_>,
        &MaterialCircuit, &PreparedElectricalContactNetwork<'_,'_>) -> R) -> R {
    let m=[5.0,10.0].map(|cp|material_with_transport(cp,branch,0.05,0.0,
        vec![(280.0,0.0),(320.0,0.0)],None,vec![(280.0,0.1),(320.0,0.1)]));
    let (heat,_)=common::initial(&m,9.0,vec![]);
    let mut beam=Bridge::new(&[(20,LENGTH)],Support::CantileverLeft,
        solid::Budget {max_cells:1,max_layers_per_cell:2,max_branch_states:4}).unwrap();
    for i in 0..2 {
        beam=beam.deposit_bonded_material_layer(LAYERS[i],geometry(i),&m[i].mechanical,&m[i].caloric,
            &BTreeMap::from([(20,300.0)]),&mut ||false).unwrap();
    }
    let context=point(300.0); let key=PropertyKey::new("conductivity",THERMAL_CONDUCTIVITY_DIMS);
    let cells=[0,1].map(|i|MaterialCell {cell:CELLS[i],material_id:MaterialId(i as u64),
        curve:m[i].caloric.curve(),point:ThermalMaterialPoint {material:&m[i].card,context:&context,
            temperature_axis:"T",liquid_fraction_axis:None},
        conductivity:Some((&key,MaterialPropertySelection::SingleClaimOnly))});
    let resistors=[CellResistor {destination:ResistorCell {branch:1,cell:CELLS[0]},
        selection:MaterialPropertySelection::SingleClaimOnly,length_m:LENGTH,area_m2:AREA}];
    let links=[ConductionLink {link:5,distance_a_m:0.0005,distance_b_m:0.0005}];
    let bulk=MaterialNetwork {cells:&cells,resistors:&resistors,links:&links};
    let card=contact_source::card(&m[0].card,&m[1].card,contact_ohm*LINK_AREA,slope);
    let interface_context=contact_source::point();
    let contact=[MaterialElectricalContact {link:5,branch:1,point:ThermalInterfacePoint {
        interface:&card,materials:[&m[0].card,&m[1].card],context:&interface_context,
        temperature_axes:["Ta","Tb"],liquid_fraction_axes:[Some("melt_a"),Some("melt_b")],
        mean_temperature_axis:None},selection:MaterialPropertySelection::SingleClaimOnly,
        heat_fraction_to_first:fraction}];
    let bases=[(1,ElectricalPropertyBasis::Resistivity)];
    let network=ElectricalContactNetwork {bulk:&bulk,contacts:&contact,electrical_bases:&bases};
    let prepared=network.prepare(&heat,||false).unwrap();
    let initial=prepared.resolve(&heat,||false).unwrap();
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:CAP}),(1,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        initial.bulk.resistors,vec![],vec![charge],&[]).unwrap()
        .rebind_contacts_candidate(initial.contacts,&[]).unwrap();
    run(&m,&heat,&beam,&circuit,&prepared)
}

#[test]
fn contact_heat_enters_both_materials_with_the_explicit_partition_and_charge_sign() {
    for fraction in [0.0,0.25,1.0] { for sign in [-1.0,1.0] {
        with_fixture(0.0,1000.0,0.0,fraction,sign*CHARGE,|m,heat,beam,circuit,prepared| {
            let next=step_with_contacts_candidate(heat,beam,&map(m),circuit,prepared,&[],DT,
                1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
            let eta=DT/(2.0*2000.0*CAP);
            let q1=sign*CHARGE*(1.0-eta)/(1.0+eta);
            let current=0.5*(sign*CHARGE+q1)/(CAP*2000.0);
            let loss=DT*current*current*1000.0;
            near(next.circuit().charges_c()[0],q1,1e-12);
            near(next.thermal().temperature_k(CELLS[0]).unwrap(),300.0+loss*(1.0+fraction)/(MASS*5.0),1e-7);
            near(next.thermal().temperature_k(CELLS[1]).unwrap(),300.0+loss*(1.0-fraction)/(MASS*10.0),1e-7);
            let r=next.report();
            near(r.electrical.contacts[0].heat_j,loss,1e-10);
            near(r.deposited_joule_j,2.0*loss,1e-10);
            assert_eq!(r.contact_heat_destinations,vec![ContactHeatCells {id:5,cells:CELLS,
                heat_fraction_to_first:fraction}]);
            assert_eq!(r.deposited_mechanical_heat_j,0.0);
            assert!(r.mechanical.energy_j.numerical_dissipation>0.0);
            assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
            assert_eq!(r.material.links[0].conductance_w_k,0.0);
            assert_eq!(heat.time_s(),0.0); assert_eq!(beam.time_s(),0.0);
        });
    }}
}

#[test]
fn changing_only_the_interface_resistance_changes_current_heating_and_sag() {
    let run=|resistance|with_fixture(3e8,resistance,0.0,0.25,CHARGE,|m,t,b,c,p| {
        let next=step_with_contacts_candidate(t,b,&map(m),c,p,&[],DT,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        (next.circuit().charges_c()[0],next.thermal().temperature_k(CELLS[0]).unwrap(),
            next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1],next.bridge().mass_kg())
    });
    let low=run(200.0); let high=run(2000.0);
    assert!(low.0<high.0 && low.1>high.1);
    assert!((low.2-high.2).abs()>1e-10);
    assert_eq!(low.3,high.3);
}

#[test]
fn partition_alone_changes_deformation_without_changing_the_finite_electrical_supply() {
    let run=|fraction|with_fixture(3e8,1000.0,0.0,fraction,CHARGE,|m,t,b,c,p| {
        let next=step_with_contacts_candidate(t,b,&map(m),c,p,&[],DT,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        (next.circuit().charges_c()[0],next.report().deposited_joule_j,
            next.bridge().deformed_position_m(LENGTH,0.0).unwrap()[1])
    });
    let a=run(0.0); let b=run(0.75);
    near(a.0,b.0,1e-12); near(a.1,b.1,1e-10);
    assert!((a.2-b.2).abs()>1e-10);
}

#[test]
fn actual_endpoint_contact_feedback_replays_the_original_circuit_and_beam() {
    with_fixture(3e8,1000.0,0.005,0.25,CHARGE,|m,t,b,c,p| {
        let mut heat=t.clone(); let mut beam=b.clone(); let mut circuit=c.clone();
        let mut work=0.0; let mut numerical=0.0;
        for n in 1..=4 {
            let target=f64::from(n)*DT;
            let g=if n<=2 {1.0}else{0.0}; let load=if n<=2 {0.02}else{0.0};
            let next=step_with_contacts_candidate(&heat,&beam,&map(m),&circuit,p,&[],target,
                g,&[load],controls(),&[0.0;2],||false).unwrap();
            let r=next.report(); let temperature=next.thermal().temperature_k(CELLS[0]).unwrap();
            near(r.electrical.contacts[0].contact.resistance_ohm(),1000.0*(1.0+0.005*(temperature-300.0)),1e-8);
            let exact=p.resolve(next.thermal(),||false).unwrap();
            let replay=circuit.step_with_contacts_candidate(exact.bulk.resistors,exact.contacts,&[],
                target-heat.time_s(),controls().circuit,||false).unwrap();
            assert_eq!(next.circuit(),replay.circuit());
            let temperatures=[0,1].into_iter().map(|i|((20,LAYERS[i]),next.thermal().temperature_k(CELLS[i]).unwrap())).collect();
            let original=beam.advance_controlled(target,&temperatures,g,&[load],
                solid::Config {section:config().mechanics,maximum_deflection_m:0.009},&mut ||false).unwrap();
            assert_eq!(next.bridge(),&original.0);
            near(r.thermal.external_heat_j,r.deposited_joule_j+r.deposited_mechanical_heat_j,1e-9);
            work+=r.mechanical.load_work_j; numerical+=r.mechanical.energy_j.numerical_dissipation;
            near(next.thermal().energy_j()+next.bridge().stored_energy_j()+r.electrical.stored_after_j
                +numerical-work,CHARGE*CHARGE/(2.0*CAP),f64::from(n)*config().total_energy_tolerance_j);
            let (t,b,c,_)=next.into_parts(); heat=t; beam=b; circuit=c;
        }
    });
}

#[test]
fn source_switch_and_late_cancellation_cannot_publish_any_of_the_three_histories() {
    with_fixture(3e8,1000.0,0.005,0.25,CHARGE,|m,t,b,c,p| {
        let saved_t=t.cells().clone(); let saved_b=b.clone(); let saved_c=c.clone();
        let other=contact_source::card(&m[0].card,&m[1].card,2.0,0.005);
        let changed=c.rebind_contacts_candidate(vec![contact_source::contact(5,1,&other,LINK_AREA)],&[]).unwrap();
        assert!(step_with_contacts_candidate(t,b,&map(m),&changed,p,&[],DT,1.0,&[0.02],controls(),&[0.0;2],||false).is_err());
        let mut polls=0;
        let next=step_with_contacts_candidate(t,b,&map(m),c,p,&[],DT,1.0,&[0.02],controls(),&[0.0;2],||{polls+=1;false}).unwrap();
        let mut actual=0;
        assert!(matches!(step_with_contacts_candidate(t,b,&map(m),c,p,&[],DT,1.0,&[0.02],controls(),&[0.0;2],
            ||{actual+=1;actual==polls}),Err(ElectricalError::Cancelled)));
        let mut limited=controls(); limited.maximum_circuit_trials=next.report().circuit_trials-1;
        assert!(step_with_contacts_candidate(t,b,&map(m),c,p,&[],DT,1.0,&[0.02],limited,&[0.0;2],||false).is_err());
        assert_eq!(t.cells(),&saved_t); assert_eq!(b,&saved_b); assert_eq!(c,&saved_c);
    });
}