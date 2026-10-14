#![cfg(feature = "enthalpy-coupling")]
//! Actual current solver -> conservative specific enthalpy -> phase relation.
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

fn curve(latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xA7;32]),
        [(0.0,300.0,0.0),(20.0,310.0,0.0),(20.0+latent,310.0,1.0),(1020.0+latent,560.0,1.0)]
        .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn config() -> StepConfig { StepConfig {max_dt_s:0.125,energy_tolerance_j:1e-8} }
fn setup(n:usize,drive:ElectrodeDrive) -> ElectricSetup {
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive},
        Electrode {lo:[n-1,0,0],hi:[n,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ]}
}
fn initial<'a>(law:&'a EquilibriumEnthalpyPhaseCurve,n:usize,h:f64,links:bool) -> (DepositionSimulation<'a>,Vec<CellId>) {
    let ids:Vec<_>=(0..n).map(|c|CellId(10+7*c as u64)).collect();
    let mut events:Vec<_>=(0..n).map(|c|TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
        id:ids[c],segment_id:c as u64,body_id:1,material_id:MaterialId(1),mass_kg:1.0,
        specific_h_j_kg:if c==1 {h} else {0.0},
        geometry:CellGeometry {start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],
            volume_m3:1.0,surface_area_m2:6.0},
        boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
    })}).collect();
    if links { for c in 1..n { events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:c as u64,a:ids[c-1],b:ids[c],area_m2:1.0,conductance_w_k:100.0,
    })}); } }
    let heat=DepositionSimulation::new(events,&[(MaterialId(1),law as &dyn ThermalLaw)],
        DepositionBudget {max_cells:n,max_links:n,max_events:2*n,max_substeps:10000}).unwrap()
        .advance_to(0.0,config(),||Ok(())).unwrap().0;
    (heat,ids)
}
fn input<'a>(domain:&'a VoxelDomain,setup:&'a ElectricSetup,ids:&'a [CellId],contacts:&'a [ElectricalContact]) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup,contacts,cell_ids:ids,max_cells:16,max_current_solves:10001,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-7,"{a} != {b}");}

#[test]
fn actual_joule_power_crosses_latent_plateau_with_fixed_reference_mass() {
    let law=curve(100.0); let (old,ids)=initial(&law,3,0.0,false);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap(); let electric=setup(3,ElectrodeDrive::Current(4.0));
    let gate=CancelGate::new();
    let (mid,r)=advance_joule_enthalpy(&old,1.0,config(),input(&domain,&electric,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    near(mid.temperature_k(ids[1]).unwrap(),310.0); near(mid.cells()[&ids[1]].specific_h_j_kg,80.0);
    near(law.state_at_specific_enthalpy(80.0).unwrap().liquid_mass_fraction(),0.6);
    near(r.electrode_energy_j,80.0);near(r.thermal.external_heat_j,80.0);
    assert_eq!(r.current_solves,r.thermal.substeps+1);
    let (end,r)=advance_joule_enthalpy(&mid,2.0,config(),input(&domain,&electric,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    near(end.temperature_k(ids[1]).unwrap(),320.0);near(end.energy_j(),160.0);
    near(r.deposited_joule_j,80.0);assert_eq!(end.mass_kg(),old.mass_kg());
    assert_eq!(old.time_s(),0.0);near(old.energy_j(),0.0);
}

#[test]
fn resistance_evolves_with_liquid_fraction_at_identical_kelvin_temperature() {
    let law=curve(100.0);let (old,ids)=initial(&law,3,40.0,false);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let electric=setup(3,ElectrodeDrive::Voltage(10.0));
    let state=|_:usize,cell:&DepositedCell,_:f64| {
        let f=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap().liquid_mass_fraction();
        Ok(Some(1.0+f))
    };
    let (next,r)=advance_joule_enthalpy(&old,0.1,StepConfig {max_dt_s:0.01,..config()},
        input(&domain,&electric,&ids,&[]),&CancelGate::new(),state).unwrap();
    near(next.temperature_k(ids[1]).unwrap(),310.0);
    assert!(r.endpoint.specific_enthalpy_j_kg[1]>40.0);
    assert!(r.endpoint.resistivity_ohm_m[1].unwrap()>1.2);
    assert!(r.endpoint.electric.electrodes[0].current_a<10.0/1.2);
    let applied=r.last_applied.unwrap();
    assert_eq!(applied.temperature_k[1],r.endpoint.temperature_k[1]);
    assert!(applied.specific_enthalpy_j_kg[1]<r.endpoint.specific_enthalpy_j_kg[1]);
    assert!(applied.electric.electrodes[0].current_a>r.endpoint.electric.electrodes[0].current_a);
    near(next.energy_j()-old.energy_j(),r.deposited_joule_j);
}

#[test]
fn contact_partition_reaches_electrode_and_bulk_thermal_cells_exactly_once() {
    let law=curve(100.0);let (old,ids)=initial(&law,4,0.0,false);
    let domain=VoxelDomain::new(4,1,1,1.0).unwrap();let electric=setup(4,ElectrodeDrive::Voltage(9.0));
    let contacts=[ElectricalContact {cells:[0,1],resistance_ohm_m2:3.0,heat_fraction_to_first:0.25}];
    let (next,r)=advance_joule_enthalpy(&old,0.25,config(),input(&domain,&electric,&ids,&contacts),
        &CancelGate::new(),|c,_,_|Ok(Some(if c==1 {2.0} else {4.0}))).unwrap();
    // R=2+4+3=9, I=1. Electrode receives .25*3; all of R_bulk stays local.
    for (c,q) in [0.75,4.25,4.0,0.0].into_iter().enumerate() {
        near(next.cells()[&ids[c]].specific_h_j_kg,q*0.25);
        assert_eq!(next.cells()[&ids[c]].boundary.source_w,0.0);
    }
    near(r.electrode_energy_j,2.25);near(r.thermal.external_heat_j,2.25);
    near(r.electrical_defect_j,0.0);near(r.transfer_defect_j,0.0);
}

#[test]
fn removing_electrical_drive_preserves_latent_energy_and_material_substitution_changes_phase() {
    let a=curve(100.0);let b=curve(200.0);let (old,ids)=initial(&a,3,0.0,false);let (other,_)=initial(&b,3,0.0,false);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let on=setup(3,ElectrodeDrive::Current(4.0));
    let off=setup(3,ElectrodeDrive::Voltage(0.0));let gate=CancelGate::new();
    let (hot,_)=advance_joule_enthalpy(&old,1.0,config(),input(&domain,&on,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    let (long,_)=advance_joule_enthalpy(&other,1.0,config(),input(&domain,&on,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    near(hot.energy_j(),long.energy_j());near(hot.temperature_k(ids[1]).unwrap(),long.temperature_k(ids[1]).unwrap());
    near(a.state_at_specific_enthalpy(hot.cells()[&ids[1]].specific_h_j_kg).unwrap().liquid_mass_fraction(),0.6);
    near(b.state_at_specific_enthalpy(long.cells()[&ids[1]].specific_h_j_kg).unwrap().liquid_mass_fraction(),0.3);
    let (retained,r)=advance_joule_enthalpy(&hot,2.0,config(),input(&domain,&off,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    near(retained.energy_j(),hot.energy_j());near(r.deposited_joule_j,0.0);near(r.electrode_energy_j,0.0);
}

#[test]
fn cfl_limited_substeps_share_electrical_and_thermal_physical_time() {
    let law=curve(100.0);let (old,ids)=initial(&law,3,0.0,true);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let electric=setup(3,ElectrodeDrive::Current(4.0));
    let (next,r)=advance_joule_enthalpy(&old,0.1,config(),input(&domain,&electric,&ids,&[]),
        &CancelGate::new(),|_,_,_|Ok(Some(5.0))).unwrap();
    assert!(r.thermal.substeps>1);assert_eq!(r.current_solves,r.thermal.substeps+1);
    near(next.energy_j(),8.0);near(r.deposited_joule_j,8.0);near(r.thermal.residual_j,0.0);
    assert!(next.cells()[&ids[0]].specific_h_j_kg>0.0); // actual conduction into electrode body
}

#[test]
fn failed_endpoint_law_work_budget_and_cancellation_preserve_original_state() {
    let law=curve(100.0);let (old,ids)=initial(&law,3,0.0,false);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let electric=setup(3,ElectrodeDrive::Current(4.0));
    let failure=ChtError::InvalidInput {field:"material.phase",reason:"endpoint unsupported".into()};
    let gate=CancelGate::new();let request=input(&domain,&electric,&ids,&[]);
    let result=advance_joule_enthalpy(&old,0.125,config(),request,&gate,|_,cell,_| {
        if cell.specific_h_j_kg>0.0 {Err(failure.clone())} else {Ok(Some(5.0))}
    });
    assert!(matches!(result,Err(JouleEnthalpyError::Electric(e)) if e==failure));
    assert!(matches!(advance_joule_enthalpy(&old,0.125,config(),JouleEnthalpyInput {max_current_solves:1,..request},
        &gate,|_,_,_|Ok(Some(5.0))),Err(JouleEnthalpyError::CurrentBudget)));
    let cancelled=CancelGate::new();
    assert!(advance_joule_enthalpy(&old,0.125,config(),request,&cancelled,|_,_,_| {
        cancelled.request();Ok(Some(5.0))
    }).is_err());
    near(old.energy_j(),0.0);assert_eq!(old.time_s(),0.0);
    let (_,r)=advance_joule_enthalpy(&old,0.125,config(),request,&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    near(r.deposited_joule_j,10.0);
}

#[test]
fn invalid_correspondence_refuses_before_material_queries_and_zero_time_only_observes() {
    let law=curve(100.0);let (old,ids)=initial(&law,3,0.0,false);
    let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let electric=setup(3,ElectrodeDrive::Current(4.0));
    let gate=CancelGate::new();
    for map in [vec![ids[0],ids[0],ids[2]],vec![ids[2],ids[1],ids[0]],vec![CellId(0);3]] {
        assert!(advance_joule_enthalpy(&old,1.0,config(),input(&domain,&electric,&map,&[]),&gate,
            |_,_,_|panic!("map admission precedes material queries")).is_err());
    }
    let (next,r)=advance_joule_enthalpy(&old,0.0,config(),input(&domain,&electric,&ids,&[]),&gate,|_,_,_|Ok(Some(5.0))).unwrap();
    assert_eq!(next.cells(),old.cells());assert!(r.last_applied.is_none());
    assert_eq!(r.current_solves,1);near(r.electrode_energy_j,0.0);
}
