//! Existing electrical source choices survive actual thermal topology growth.
use super::*;
use fs_couple::circuit_enthalpy::material::MaterialNetworkError;

fn competing_material(expired: bool) -> Material {
    material_with_transport_alternatives(5.0,3e8,0.05,0.0,
        vec![(280.0,0.002),(320.0,0.002)],None,
        (vec![(280.0,0.1),(if expired {300.02}else{320.0},0.1)],
            Some(vec![(300.01,0.3),(320.0,0.3)])))
}

#[test]
fn electrical_only_cell_keeps_its_claim_when_growth_adds_its_first_heat_contact() {
    let m=[competing_material(false),material(10.0,3e8,0.05,0.0)];
    let old=seed(&m,vec![arrival(99,3,LENGTH,0.02,1),TimedEvent {time_s:0.2,event:Event::Contact(
        ThermalLink {id:5,a:CellId(90),b:CellId(99),area_m2:AREA,conductance_w_k:9.0})}],
        vec![birth(&m[1],99,3,0.02)]);
    let context=point(300.0); let key=key(); let mut bound=cells(&m,&context,&key,&[(90,0)]);
    bound[0].conductivity=None;
    let resistors=resistors();
    let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    let initial=prepared.resolve(old.thermal(),||false).unwrap();
    let claim=initial.resistors[0].conductor.material().property("electrical_resistivity").unwrap().answer().receipt.selected;
    assert!(initial.source.cells[0].conductivity.is_none());
    let circuit=circuit(&prepared,old.thermal());
    let (born,charged,r)=old.step_electrical_candidate(&circuit,&prepared,&[],0.2,1.0,&[0.02],
        controls(),&[0.0],||false).unwrap();
    assert!(born.thermal().temperature_k(CellId(90)).unwrap()>300.01);
    let expanded=cells(&m,&context,&key,&[(90,0),(99,1)]);
    let links=[ConductionLink {link:5,distance_a_m:0.05,distance_b_m:0.01}];
    let successor=MaterialNetwork {cells:&expanded,resistors:&resistors,links:&links};
    assert!(successor.prepare(born.thermal(),||false).is_err(),"two electrical fits are now eligible");
    assert!(prepared.extend_thermal(&successor,born.thermal(),||false).is_err());
    let continued=prepared.extend_fixed_circuit(&successor,born.thermal(),||false).unwrap();
    let actual=continued.resolve(born.thermal(),||false).unwrap();
    assert_eq!(actual.resistors[0].conductor.material().property("electrical_resistivity").unwrap().answer().receipt.selected,claim);
    near(actual.resistors[0].conductor.resistance_ohm(),0.1*LENGTH/AREA,1e-9);
    assert!(actual.source.cells.iter().all(|s|s.conductivity.is_some()));
    let (after,final_circuit,r2)=born.step_electrical_candidate(&charged,&continued,&[],0.4,1.0,&[0.02;2],
        controls(),&[0.0;2],||false).unwrap();
    assert!(final_circuit.charges_c()[0]<charged.charges_c()[0]);
    assert!(r2.interval.layer_energy_j[&(3,LAYERS[0])].dissipated>0.0);
    near(after.thermal().energy_j()+after.bridge().stored_energy_j()+r2.interval.electrical.stored_after_j
        +r.interval.mechanical.energy_j.numerical_dissipation+r2.interval.mechanical.energy_j.numerical_dissipation
        -r.interval.mechanical.load_work_j-r2.interval.mechanical.load_work_j,
        Q0*Q0/(2.0*CAP),2.0*config().total_energy_tolerance_j);
    assert_eq!(circuit.charges_c(),&[Q0]); assert_eq!(old.thermal().time_s(),0.0);
}

#[test]
fn expired_electrical_source_refuses_instead_of_selecting_the_new_hot_fit() {
    let m=[competing_material(true),material(10.0,3e8,0.05,0.0)];
    let old=seed(&m,vec![TimedEvent {time_s:0.0,event:Event::Boundary(CellId(90),
        Boundary {source_w:0.01,..boundary()})}],vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]);
    let resistors=resistors(); let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    // Warm outside the electrical lane to challenge its original validity domain.
    let (warm,_) = old.thermal().advance_to(0.1,
        StepConfig {max_dt_s:0.05,energy_tolerance_j:1e-9},||Ok(())).unwrap();
    assert!(warm.temperature_k(CellId(90)).unwrap()>300.02);
    assert!(network.prepare(&warm,||false).is_ok());
    assert!(prepared.extend_fixed_circuit(&network,&warm,||false).is_err());
    assert!(prepared.resolve(old.thermal(),||false).is_ok());
    assert_eq!(old.thermal().time_s(),0.0);
}

#[test]
fn extension_cannot_rewire_heat_or_change_resistor_geometry_and_material_conditions() {
    let m=materials(); let old=seed(&m,vec![],vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0)]);
    let original=resistors(); let network=MaterialNetwork {cells:&bound,resistors:&original,links:&[]};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    for mode in 0..3 {
        let mut changed=resistors();
        match mode {
            0=>changed[0].destination.branch=2,
            1=>changed[0].destination.cell=CellId(99),
            _=>{changed[0].length_m*=2.0;changed[0].area_m2*=0.5;},
        }
        let successor=MaterialNetwork {cells:&bound,resistors:&changed,links:&[]};
        assert!(matches!(prepared.extend_fixed_circuit(&successor,old.thermal(),||false),
            Err(MaterialNetworkError::Invalid(_))));
    }
    let removed=MaterialNetwork {cells:&bound,resistors:&[],links:&[]};
    assert!(prepared.extend_fixed_circuit(&removed,old.thermal(),||false).is_err());
    let [first]=resistors(); let [second]=resistors(); let extra=[first,second];
    let added=MaterialNetwork {cells:&bound,resistors:&extra,links:&[]};
    assert!(prepared.extend_fixed_circuit(&added,old.thermal(),||false).is_err());
    let changed_context=context.clone().with("lot",12.0).unwrap();
    let changed=cells(&m,&changed_context,&key,&[(90,0)]);
    let retargeted=MaterialNetwork {cells:&changed,resistors:&original,links:&[]};
    assert!(prepared.extend_fixed_circuit(&retargeted,old.thermal(),||false).is_err());
    assert_eq!(prepared.resistor_destinations(),vec![ResistorCell {branch:1,cell:CellId(90)}]);
}

#[test]
fn late_cancellation_leaves_old_sources_usable_and_reordering_replays() {
    let m=materials(); let mut passive=arrival(55,6,0.0,LENGTH,1);
    passive.time_s=0.0; if let Event::Birth(c)=&mut passive.event {c.body_id=2;}
    let old=seed(&m,vec![passive],vec![]);
    let context=point(300.0); let key=key(); let bound=cells(&m,&context,&key,&[(90,0),(55,1)]);
    let resistors=resistors(); let network=MaterialNetwork {cells:&bound,resistors:&resistors,links:&[]};
    let prepared=network.prepare(old.thermal(),||false).unwrap();
    let before=prepared.resolve(old.thermal(),||false).unwrap();
    let reversed=cells(&m,&context,&key,&[(55,1),(90,0)]);
    let successor=MaterialNetwork {cells:&reversed,resistors:&resistors,links:&[]};
    let mut polls=0;
    let next=prepared.extend_fixed_circuit(&successor,old.thermal(),||{polls+=1;false}).unwrap();
    let actual=next.resolve(old.thermal(),||false).unwrap();
    assert_eq!(actual.source,before.source); assert_eq!(actual.resistors,before.resistors);
    let mut count=0;
    assert!(matches!(prepared.extend_fixed_circuit(&successor,old.thermal(),||{count+=1;count==polls}),
        Err(MaterialNetworkError::Thermal(DepositionError::Cancelled))));
    let replay=prepared.extend_fixed_circuit(&successor,old.thermal(),||false).unwrap()
        .resolve(old.thermal(),||false).unwrap();
    assert_eq!(replay.source,before.source); assert_eq!(replay.resistors,before.resistors);
    assert_eq!(prepared.resolve(old.thermal(),||false).unwrap().source,before.source);
}
