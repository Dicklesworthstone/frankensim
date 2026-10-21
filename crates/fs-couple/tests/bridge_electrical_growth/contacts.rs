//! The original contact-current owner participates in the same birth transaction.
use super::*;
use fs_couple::circuit_enthalpy::contacts::{ContactHeatCells, material::{
    ElectricalContactNetwork, MaterialElectricalContact, PreparedElectricalContactNetwork,
}};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::thermal::interface::ThermalInterfacePoint;

// The contact joins the original single-layer beam to a passive thermal body.
// It does not turn that body into a mechanical support or a future beam segment.
fn with_contact<R>(bad_last: bool,
    run: impl FnOnce(&[Material;2], &GrowingBridge<'_,'_>, &MaterialCircuit,
        &PreparedElectricalContactNetwork<'_,'_>) -> R) -> R {
    let m=materials();
    let mut passive=arrival(55,6,0.0,LENGTH,1);
    passive.time_s=0.0;
    if let Event::Birth(c)=&mut passive.event { c.body_id=2; }
    let mut events=vec![passive,TimedEvent {time_s:0.0,event:Event::Contact(
        ThermalLink {id:7,a:CellId(90),b:CellId(55),area_m2:0.001,conductance_w_k:9.0})},
        arrival(99,3,LENGTH,0.02,1)];
    let mut future=vec![birth(&m[1],99,3,0.02)];
    if bad_last {
        let mut bad=arrival(2,1,0.12,0.015,0);
        if let Event::Birth(c)=&mut bad.event { c.mass_kg*=2.0; }
        events.push(bad);
        future.push(birth(&m[0],2,1,0.015));
    }
    let old=seed(&m,events,future);
    let context=point(300.0); let key=key();
    let bound=cells(&m,&context,&key,&[(90,0),(55,1)]);
    let resistors=resistors();
    let links=[ConductionLink {link:7,distance_a_m:0.0005,distance_b_m:0.0005}];
    let bulk=MaterialNetwork {cells:&bound,resistors:&resistors,links:&links};
    let card=contact_source::card(&m[0].card,&m[1].card,1.0,0.005);
    let interface_context=contact_source::point();
    let contacts=[MaterialElectricalContact {link:7,branch:1,
        point:ThermalInterfacePoint {interface:&card,materials:[&m[0].card,&m[1].card],
            context:&interface_context,temperature_axes:["Ta","Tb"],
            liquid_fraction_axes:[Some("melt_a"),Some("melt_b")],mean_temperature_axis:None},
        selection:MaterialPropertySelection::SingleClaimOnly,heat_fraction_to_first:0.25}];
    let bases=[(1,ElectricalPropertyBasis::Resistivity)];
    let network=ElectricalContactNetwork {bulk:&bulk,contacts:&contacts,electrical_bases:&bases};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    let initial=prepared.resolve(old.thermal(),||false).unwrap();
    let circuit=MaterialCircuit::new(CircuitGraph {node_count:2,branches:vec![
        (1,0,Branch::Capacitor {farads:CAP}),(1,0,Branch::Resistor {ohms:999.0})],transformers:vec![]},
        initial.bulk.resistors,vec![],vec![Q0],&[]).unwrap()
        .rebind_contacts_candidate(initial.contacts,&[]).unwrap();
    run(&m,&old,&circuit,&prepared)
}

#[test]
fn contact_interval_and_tip_birth_replay_original_owners_without_an_extra_electrical_step() {
    with_contact(false,|m,old,circuit,p| {
        let left=electrical::step_with_contacts_candidate(old.thermal(),old.bridge(),old.bindings(),
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        let manual_heat=left.thermal().advance_to(0.2,
            StepConfig {max_dt_s:0.1,energy_tolerance_j:1e-9},||Ok(())).unwrap().0;
        let manual_beam=left.bridge().append_cantilever_controlled(&[CellBirth {id:3,length_m:0.02,
            mass_kg:manual_heat.cells()[&CellId(99)].mass_kg,
            temperature_k:manual_heat.temperature_k(CellId(99)).unwrap(),law:m[1].mechanical.law()}],
            mechanics(),&mut ||false).unwrap().0;
        let (born,next_circuit,r)=old.step_electrical_with_contacts_candidate(
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        assert_eq!(&next_circuit,left.circuit());
        assert_eq!(born.thermal().cells(),manual_heat.cells());
        assert_eq!(born.bridge(),&manual_beam);
        assert!(r.interval.electrical.contacts[0].heat_j>0.0);
        assert_eq!(r.interval.contact_heat_destinations,vec![ContactHeatCells {
            id:7,cells:[CellId(90),CellId(55)],heat_fraction_to_first:0.25}]);
        near(r.interval.deposited_joule_j,left.report().deposited_joule_j,1e-12);
        assert!(born.thermal().temperature_k(CellId(55)).unwrap()>300.0);
        assert_eq!(born.thermal().cells().len(),3);
        assert_eq!(born.bridge().cells().len(),2,"the passive contact receiver is not a beam cell");
        assert_eq!(born.bridge().cells()[1].section().stored_energy_j(),0.0);
        assert_eq!(r.births.thermal.events,1); assert_eq!(r.births.thermal.substeps,0);
        assert!(r.birth_equilibrium_pending);
        assert!(r.total_energy_residual_j.abs()<config().total_energy_tolerance_j);
        assert_eq!(old.thermal().time_s(),0.0); assert_eq!(circuit.charges_c(),&[Q0]);
    });
}

#[test]
fn rejected_last_birth_preserves_contact_receivers_finite_supply_and_old_beam() {
    with_contact(true,|_,old,circuit,p| {
        let left=electrical::step_with_contacts_candidate(old.thermal(),old.bridge(),old.bindings(),
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        assert!(left.circuit().charges_c()[0]<Q0);
        assert!(left.thermal().temperature_k(CellId(55)).unwrap()>300.0);
        let before=old.clone(); let original=circuit.clone();
        assert!(matches!(old.step_electrical_with_contacts_candidate(
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||false),
            Err(ElectricalGrowthError::Birth(_))));
        assert_eq!(circuit,&original); assert_eq!(old.bridge(),before.bridge());
        assert_eq!(old.thermal().cells(),before.thermal().cells());
        assert_eq!(old.thermal().links(),before.thermal().links());
        assert_eq!(old.remaining_births(),2);
    });
}

#[test]
fn final_contact_growth_cancellation_preserves_all_histories_and_replays() {
    with_contact(false,|_,old,circuit,p| {
        let mut polls=0;
        let (next,c,_) = old.step_electrical_with_contacts_candidate(
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||{polls+=1;false}).unwrap();
        let mut actual=0;
        assert!(matches!(old.step_electrical_with_contacts_candidate(
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||{actual+=1;actual==polls}),
            Err(ElectricalGrowthError::Birth(BridgeEnthalpyError::Coupled(StripEnthalpyError::Cancelled)))));
        let (replay,rc,_) = old.step_electrical_with_contacts_candidate(
            circuit,p,&[],0.2,1.0,&[0.02],controls(),&[0.0;2],||false).unwrap();
        assert_eq!(next.bridge(),replay.bridge()); assert_eq!(next.thermal().cells(),replay.thermal().cells());
        assert_eq!(c,rc); assert_eq!(circuit.charges_c(),&[Q0]);
        assert_eq!(old.thermal().time_s(),0.0); assert_eq!(old.remaining_births(),1);
    });
}
