#![cfg(feature = "enthalpy-coupling")]
use std::cell::Cell;
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::*;
use fs_lbm::conjugate::enthalpy::transport::*;
use fs_material::phase::{EnthalpyPhaseKnot,EquilibriumEnthalpyPhaseCurve};

fn curve()->EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xD1;32]),
        [(0.0,300.0,0.0),(20.0,310.0,0.0),(120.0,310.0,1.0),(1120.0,810.0,1.0)]
            .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn config()->StepConfig {StepConfig {max_dt_s:0.01,energy_tolerance_j:1e-8}}
fn setup(v:f64)->ElectricSetup {
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive:ElectrodeDrive::Voltage(v)},
        Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)}]}
}
fn initial(law:&EquilibriumEnthalpyPhaseCurve,h:f64)->DepositionSimulation<'_> {
    let mut events:Vec<_>=(0..3).map(|c|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:CellId(c),segment_id:c,body_id:1,material_id:MaterialId(1),mass_kg:1.0,
        specific_h_j_kg:if c==1 {h} else {0.0},
        geometry:CellGeometry {start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],
            volume_m3:1.0,surface_area_m2:6.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
    })}).collect();
    for c in 0..2 {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:c,a:CellId(c),b:CellId(c+1),area_m2:1.0,conductance_w_k:1.0})});}
    DepositionSimulation::new(events,&[(MaterialId(1),law as &dyn ThermalLaw)],DepositionBudget {
        max_cells:3,max_links:2,max_events:8,max_substeps:10000,
    }).unwrap().advance_to(0.0,config(),||Ok(())).unwrap().0
}
const IDS:[CellId;3]=[CellId(0),CellId(1),CellId(2)];
fn input<'a>(domain:&'a VoxelDomain,e:&'a ElectricSetup,contacts:&'a [ElectricalContact])->JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup:e,contacts,cell_ids:&IDS,max_cells:3,max_current_solves:10001,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn series<F>(p:&ContactPatch,s:[&EnthalpyTransportSample<F>;2])->Result<(f64,()),JouleEnthalpyError> {
    let [a,b]=s.map(|x|x.response.conductivity_w_m_k());
    Ok((if a==0.0 || b==0.0 {0.0} else {p.link.area_m2/(0.5/a+0.5/b)},()))
}
fn near(a:f64,b:f64){assert!((a-b).abs()<1e-7,"{a} != {b}");}

#[test]
fn constant_shared_transport_preserves_the_original_current_and_enthalpy_path() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(10.0);
    let gate=CancelGate::new();let request=input(&d,&e,&[]);
    let (a,ra)=advance_joule_enthalpy(&old,0.2,config(),request,&gate,|_,_,_|Ok(Some(2.0))).unwrap();
    let (b,rb)=advance_joule_enthalpy_with_material_transport(&old,0.2,config(),request,&[0,1],&gate,
        |_,_,_|EnthalpyTransport::try_new((),Some(2.0),1.0),series).unwrap();
    assert_eq!(a.cells(),b.cells());assert_eq!(ra,rb.joule);assert_eq!(a.links(),b.links());
    assert_eq!(rb.last_applied.as_ref().unwrap().time_s,rb.joule.last_applied.as_ref().unwrap().time_s);
    assert_eq!(rb.endpoint.time_s,0.2);
}

#[test]
fn one_nonclone_phase_sample_supplies_both_conduction_and_current_on_the_plateau() {
    struct Payload {fraction:f64,query:usize}
    let law=curve();let old=initial(&law,40.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(10.0);
    let queries=Cell::new(0usize);let links=Cell::new(0usize);
    let (next,r)=advance_joule_enthalpy_with_material_transport(&old,0.2,config(),input(&d,&e,&[]),&[0,1],&CancelGate::new(),
        |_,c,t| {
            let phase=law.state_at_specific_enthalpy(c.specific_h_j_kg).unwrap();near(phase.temperature_k(),t);
            let fraction=phase.liquid_mass_fraction();let query=queries.get();queries.set(query+1);
            EnthalpyTransport::try_new(Payload {fraction,query},Some(1.0+fraction),0.1+2.0*fraction)
        },|p,s| {
            links.set(links.get()+1);
            for a in s {near(a.response.conductivity_w_m_k(),0.1+2.0*a.response.material().fraction);}
            let g=series(p,s)?.0;
            Ok((g,[s[0].response.material().query,s[1].response.material().query]))
        }).unwrap();
    assert_eq!(queries.get(),3*r.joule.current_solves);assert_eq!(links.get(),2*r.joule.current_solves);
    let used=r.last_applied.as_ref().unwrap();let end=&r.endpoint;
    assert_eq!(used.cells[1].temperature_k,end.cells[1].temperature_k);
    assert!(end.cells[1].response.material().fraction>used.cells[1].response.material().fraction);
    assert!(end.cells[1].response.conductivity_w_m_k()>used.cells[1].response.conductivity_w_m_k());
    assert!(r.joule.endpoint.electric.electrodes[0].current_a<r.joule.last_applied.as_ref().unwrap().electric.electrodes[0].current_a);
    near(next.energy_j()-old.energy_j(),r.joule.deposited_joule_j);
    assert_eq!(next.mass_kg(),old.mass_kg());assert_eq!(next.links(),old.links());
}

#[test]
fn thermal_material_data_changes_phase_and_subsequent_current_with_unchanged_electrical_law() {
    let law=curve();let old=initial(&law,40.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(10.0);let gate=CancelGate::new();
    let run=|k:f64|advance_joule_enthalpy_with_material_transport(&old,0.5,config(),input(&d,&e,&[]),&[0,1],&gate,
        |_,c,_| {
            let f=law.state_at_specific_enthalpy(c.specific_h_j_kg).unwrap().liquid_mass_fraction();
            EnthalpyTransport::try_new(f,Some(1.0+f),k)
        },series).unwrap();
    let (a,ra)=run(0.1);let (b,rb)=run(2.0);
    assert!((a.cells()[&IDS[1]].specific_h_j_kg-b.cells()[&IDS[1]].specific_h_j_kg).abs()>0.01);
    assert!((ra.joule.endpoint.electric.electrodes[0].current_a-rb.joule.endpoint.electric.electrodes[0].current_a).abs()>1e-4);
    for (s,r) in [(&a,&ra),(&b,&rb)] {near(s.energy_j()-old.energy_j(),r.joule.deposited_joule_j);}
    let off=setup(0.0);
    let (retained,r)=advance_joule_enthalpy_with_material_transport(&a,0.6,config(),input(&d,&off,&[]),&[0,1],&gate,
        |_,_,_|EnthalpyTransport::try_new((),Some(1.0),0.1),series).unwrap();
    near(retained.energy_j(),a.energy_j());assert_eq!(r.joule.deposited_joule_j,0.0);
    assert!(retained.cells()[&IDS[0]].specific_h_j_kg>a.cells()[&IDS[0]].specific_h_j_kg);
}

#[test]
fn contact_heat_is_deposited_once_and_thermal_insulation_keeps_coverage() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);
    let contacts=[ElectricalContact {cells:[0,1],resistance_ohm_m2:3.0,heat_fraction_to_first:0.25}];
    let (next,r)=advance_joule_enthalpy_with_material_transport(&old,0.1,config(),input(&d,&e,&contacts),&[0,1],&CancelGate::new(),
        |_,_,_|EnthalpyTransport::try_new((),Some(2.0),0.0),series).unwrap();
    for (c,q) in [0.75,4.25,0.0].into_iter().enumerate() {near(next.cells()[&IDS[c]].specific_h_j_kg,0.1*q);}
    near(r.joule.deposited_joule_j,0.5);near(r.joule.thermal.external_heat_j,0.5);
    assert_eq!(next.contact_patch(0),old.contact_patch(0));assert_eq!(next.exposed_area_m2(IDS[1]).unwrap(),4.0);
}

#[test]
fn current_work_and_material_queries_follow_conductivity_limited_substeps() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);
    let (next,r)=advance_joule_enthalpy_with_material_transport(&old,0.1,StepConfig {max_dt_s:1.0,..config()},
        input(&d,&e,&[]),&[0,1],&CancelGate::new(),|_,_,_|EnthalpyTransport::try_new((),Some(2.0),1000.0),series).unwrap();
    assert!(r.joule.thermal.substeps>100);assert_eq!(r.joule.current_solves,r.joule.thermal.substeps+1);
    near(r.joule.deposited_joule_j,1.25);near(next.energy_j(),1.25);
    assert!(next.cells()[&IDS[0]].specific_h_j_kg>0.0);
}

#[test]
fn late_link_refusal_budget_and_cancellation_do_not_publish_a_partial_state() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);
    let gate=CancelGate::new();let input=input(&d,&e,&[]);let cells=old.cells().clone();let links=old.links().clone();
    let failure=ChtError::InvalidInput {field:"thermal.interface",reason:"endpoint unsupported".into()};
    let result=advance_joule_enthalpy_with_material_transport(&old,0.01,config(),input,&[0,1],&gate,
        |_,_,_|EnthalpyTransport::try_new((),Some(2.0),1.0),|p,s| {
            if s.iter().any(|v|v.cell.specific_h_j_kg>0.0) {Err(failure.clone().into())} else {series(p,s)}
        });
    assert!(matches!(result,Err(JouleEnthalpyError::Electric(e)) if e==failure));
    assert!(matches!(advance_joule_enthalpy_with_material_transport(&old,0.01,config(),JouleEnthalpyInput {max_current_solves:1,..input},
        &[0,1],&gate,|_,_,_|EnthalpyTransport::try_new((),Some(2.0),1.0),series),Err(JouleEnthalpyError::CurrentBudget)));
    let result=advance_joule_enthalpy_with_material_transport(&old,0.01,config(),input,&[0,1],&gate,
        |_,c,_|if c.specific_h_j_kg>0.0 {Err(DepositionError::Cancelled.into())}
            else {EnthalpyTransport::try_new((),Some(2.0),1.0)},series);
    assert!(matches!(result,Err(JouleEnthalpyError::Thermal(DepositionError::Cancelled))));
    assert_eq!(old.cells(),&cells);assert_eq!(old.links(),&links);assert_eq!(old.time_s(),0.0);
}

#[test]
fn selected_link_admission_is_explicit_and_zero_duration_checks_support() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(0.0);let gate=CancelGate::new();
    for ids in [&[0,0][..],&[4][..]] {
        assert!(advance_joule_enthalpy_with_material_transport(&old,0.0,config(),input(&d,&e,&[]),ids,&gate,
            |_,_,_|EnthalpyTransport::try_new((),Some(1.0),1.0),series).is_err());
    }
    let (next,r)=advance_joule_enthalpy_with_material_transport(&old,0.0,config(),input(&d,&e,&[]),&[0],&gate,
        |_,_,_|EnthalpyTransport::try_new((),Some(1.0),2.0),series).unwrap();
    assert!(r.last_applied.is_none());assert!(r.joule.last_applied.is_none());assert_eq!(r.joule.current_solves,1);
    assert_eq!(r.endpoint.links.len(),1);assert_eq!(next.cells(),old.cells());
    assert!(EnthalpyTransport::try_new((),None,f64::NAN).is_err());
}
