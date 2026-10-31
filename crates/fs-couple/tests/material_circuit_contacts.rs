//! Source-bound contact drops in the actual descriptor, not a replacement RC solver.
#[path="circuit_enthalpy/common.rs"] mod common;
#[path="electrical_contacts/source.rs"] mod source;
use common::*;
use fs_couple::material_circuit::{MaterialCircuit,MaterialCircuitError,MaterialCircuitSolve};
use fs_couple::material_circuit::contacts::{SeriesElectricalContact,MAX_SERIES_CONTACTS};
use fs_matdb::MaterialCard;

fn fixture() -> (MaterialCircuit,[MaterialCard;2]) {
    let (a,ca)=material(1.0,1.0,0.0);let (b,cb)=material(3.0,2.0,0.0);
    let cards=[a,b];let curves=[ca,cb];let t=initial(&curves,0.0,None);
    (circuit(&t,&cards,&curves),cards)
}
fn controls() -> MaterialCircuitSolve { MaterialCircuitSolve {energy_tolerance_j:1e-8,power_tolerance_w:1e-8} }
fn close(a:f64,b:f64) {assert!((a-b).abs()<1e-8,"{a} != {b}");}

#[test]
fn source_contacts_change_discharge_and_retain_separate_bulk_and_interface_heat() {
    let (base,cards)=fixture();
    let card=source::card(&cards[0],&cards[1],2.0,0.0);
    let contact=source::contact(7,1,&card,1.0);
    let circuit=base.rebind_contacts_candidate(vec![contact.clone()],&[]).unwrap();
    let saved=circuit.clone();let dt=0.25;
    let proposal=circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),vec![contact.clone()],&[],dt,controls(),||false).unwrap();
    // Independent implicit-midpoint RC solution: C=1, total R=1+2+3=6.
    let alpha=dt/12.0;let charge=10.0*(1.0-alpha)/(1.0+alpha);
    let current=0.5*(10.0+charge)/6.0;
    close(proposal.circuit().charges_c()[0],charge);
    let report=proposal.step();
    close(report.resistors[0].heat_j,dt*current*current);
    close(report.resistors[1].heat_j,dt*current*current*3.0);
    close(report.contacts[0].heat_j,dt*current*current*2.0);
    close(report.contacts[0].voltage_v,2.0*current);
    close(report.contacts[0].current_a,current);
    close(report.resistors[0].voltage_v,current);
    close(report.heat_j,0.5*(100.0-charge*charge));
    assert_eq!(&report.contacts[0].contact,&contact);
    assert_eq!(report.resistors[0].material,base.resistors()[0]);
    assert_eq!(circuit,saved);
}

#[test]
fn area_zero_resistance_and_signed_current_have_distinct_physical_meanings() {
    let (base,cards)=fixture();
    let run=|r,area,charge| {
        let card=source::card(&cards[0],&cards[1],r,0.0);
        let contact=source::contact(7,1,&card,area);
        let circuit=MaterialCircuit::new_with_contacts(base.graph().clone(),base.resistors().to_vec(),vec![contact.clone()],
            vec![],vec![charge],&[]).unwrap();
        circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),vec![contact],&[],0.1,controls(),||false).unwrap()
    };
    let a=run(2.0,1.0,10.0);let b=run(2.0,0.5,10.0);let negative=run(2.0,1.0,-10.0);
    assert!(b.circuit().charges_c()[0]>a.circuit().charges_c()[0]);
    close(a.step().heat_j,negative.step().heat_j);
    assert!(negative.step().contacts[0].current_a<0.0 && negative.step().contacts[0].voltage_v<0.0);
    let zero=run(0.0,1.0,10.0);
    let original=base.step_candidate(base.resistors().to_vec(),&[],0.1,controls(),||false).unwrap();
    close(zero.circuit().charges_c()[0],original.circuit().charges_c()[0]);
    assert_eq!(zero.step().contacts[0].heat_j,0.0);
    let no_contacts=base.step_with_contacts_candidate(base.resistors().to_vec(),vec![],&[],0.1,controls(),||false).unwrap();
    assert_eq!(no_contacts.circuit(),original.circuit());
    assert_eq!(no_contacts.step().resistors,original.step().resistors);
}

#[test]
fn missing_destinations_topology_and_late_cancel_never_silently_drop_a_contact() {
    let (base,cards)=fixture();let card=source::card(&cards[0],&cards[1],2.0,0.0);
    let contact=source::contact(7,1,&card,1.0);
    let circuit=base.rebind_contacts_candidate(vec![contact.clone()],&[]).unwrap();
    let saved=circuit.clone();
    assert!(circuit.step_candidate(circuit.resistors().to_vec(),&[],0.1,controls(),||false).is_err());
    assert!(circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),vec![],&[],0.1,controls(),||false).is_err());
    assert!(base.rebind_contacts_candidate(vec![contact.clone(),contact.clone()],&[]).is_err());
    assert!(base.rebind_contacts_candidate(vec![contact.clone();MAX_SERIES_CONTACTS+1],&[]).is_err());
    assert!(base.rebind_contacts_candidate(vec![source::contact(7,0,&card,1.0)],&[]).is_err());
    for area in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(SeriesElectricalContact::new(7,1,contact.source().clone(),area).is_err());
    }
    let mut polls=0;
    circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),vec![contact.clone()],&[],0.1,controls(),||{polls+=1;false}).unwrap();
    let mut actual=0;
    assert!(matches!(circuit.step_with_contacts_candidate(circuit.resistors().to_vec(),vec![contact],&[],0.1,controls(),||{
        actual+=1;actual==polls
    }),Err(MaterialCircuitError::Cancelled)));
    assert_eq!(circuit,saved);
    assert_eq!(circuit.rebind_contacts_candidate(vec![],&[]).unwrap(),base);
}

#[test]
fn source_bound_storage_rebinding_preserves_contact_sources_and_exact_charge() {
    use fs_material::field_storage::*;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_couple::material_circuit::storage::{MaterialStorage,MaterialStorageElement};
    use fs_matdb::QueryPoint;
    let (base,cards)=fixture();
    let storage=|value| {let card=source::scalar_card(ELECTRIC_PERMITTIVITY_PROPERTY,ELECTRIC_PERMITTIVITY_DIMS,value);
        MaterialStorage {branch:0,element:MaterialStorageElement::Capacitor(ResolvedCapacitor::resolve(&card,&QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,1.0,1.0).unwrap())}};
    let circuit=MaterialCircuit::new_with_storage(base.graph().clone(),base.resistors().to_vec(),vec![storage(1.0)],
        vec![],vec![10.0],&[]).unwrap();
    let card=source::card(&cards[0],&cards[1],2.0,0.0);
    let contact=source::contact(7,1,&card,1.0);
    let circuit=circuit.rebind_contacts_candidate(vec![contact.clone()],&[]).unwrap();
    let changed=circuit.rebind_storage_candidate(vec![storage(2.0)],&[],||false).unwrap();
    assert_eq!(changed.circuit().contacts(),&[contact]);
    assert_eq!(changed.circuit().charges_c(),&[10.0]);
    close(changed.parameter_work_j(),-25.0);
    assert_eq!(circuit.storage().len(),1);
}

#[test]
fn contact_heat_and_mutual_field_energy_survive_contact_reconfiguration() {
    use fs_phs::circuit::{Branch,CircuitGraph,MutualInductance};
    let (base,cards)=fixture();
    let graph=CircuitGraph {node_count:3,branches:vec![
        (1,0,Branch::Inductor {henries:2.0}),
        (1,0,Branch::Resistor {ohms:999.0}),
        (2,0,Branch::Resistor {ohms:999.0}),
        (2,0,Branch::Inductor {henries:3.0}),
    ],transformers:vec![]};
    let mutual=vec![MutualInductance {branches:(0,3),henries:1.0}];
    let original=MaterialCircuit::new_with_mutual_inductance(graph,base.resistors().to_vec(),
        mutual.clone(),vec![2.0,1.0],vec![],&[]).unwrap(); // i=(1,0), H=1 J
    let card=source::card(&cards[0],&cards[1],2.0,0.0);
    let contact=source::contact(7,1,&card,1.0);
    let bound=original.rebind_contacts_candidate(vec![contact.clone()],&[]).unwrap();
    assert_eq!(bound.mutual_inductance(),mutual);
    assert_eq!(bound.rebind_contacts_candidate(vec![],&[]).unwrap(),original);
    let step=bound.step_with_contacts_candidate(bound.resistors().to_vec(),vec![contact],
        &[],0.25,controls(),||false).unwrap();
    // The primary loop has 1 ohm bulk + 2 ohm interface. The induced secondary
    // current heats only its own 3 ohm bulk; no shared field energy is heat.
    let aa=2.375;let dd=3.375;
    let i1=(dd*1.625-1.0)/(aa*dd-1.0);
    let i2=(aa-1.625)/(aa*dd-1.0);
    let mid=0.5*(1.0+i1);
    close(step.step().resistors[0].heat_j,0.25*mid*mid);
    close(step.step().contacts[0].heat_j,0.5*mid*mid);
    close(step.step().resistors[1].heat_j,0.75*(0.5*i2)*(0.5*i2));
    close(step.step().stored_after_j+step.step().heat_j,1.0);
    assert_eq!(step.circuit().mutual_inductance(),mutual);
}
