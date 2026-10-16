#![cfg(feature = "enthalpy-coupling")]
//! G0/G1/G3/G4: phase-controlled electrical connectivity on ONE thermal history.
//! Synthetic constitutive thresholds, not measured fuse or contact mechanics.
use std::collections::BTreeMap;
use std::cell::Cell;
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_conduction::observers::ObserverCouplingError;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricalContact, ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::enthalpy::transport::{EnthalpyTransport, EnthalpyTransportSample};
use fs_lbm::conjugate::enthalpy::transport::interfaces::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

const IDS: [CellId;3] = [CellId(11),CellId(23),CellId(41)];
const FACE: ElectricalContact = ElectricalContact {
    cells:[0,1], resistance_ohm_m2:3.0, heat_fraction_to_first:0.25,
};
fn config(dt:f64) -> StepConfig { StepConfig {max_dt_s:dt,energy_tolerance_j:1e-8} }
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-8,"{a} != {b}");}
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xD7;32]),
        [(0.0,300.0,0.0),(20.0,310.0,0.0),(120.0,310.0,1.0),(1120.0,810.0,1.0)]
        .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn initial(law:&EquilibriumEnthalpyPhaseCurve, h:[f64;3], link:bool) -> DepositionSimulation<'_> {
    let mut events:Vec<_>=IDS.iter().enumerate().map(|(c,&id)|TimedEvent {time_s:0.0,
        event:Event::Birth(DepositedCell {id,segment_id:c as u64,body_id:1,material_id:MaterialId(7),
            mass_kg:1.0,specific_h_j_kg:h[c],geometry:CellGeometry {
                start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],volume_m3:1.0,surface_area_m2:6.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})}).collect();
    if link {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:17,a:IDS[0],b:IDS[1],area_m2:1.0,conductance_w_k:1.0})});}
    DepositionSimulation::new(events,&[(MaterialId(7),law as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:1,max_events:4,max_substeps:10000}).unwrap()
        .advance_to(0.0,config(0.125),||Ok(())).unwrap().0
}
fn setup(drive:ElectrodeDrive) -> ElectricSetup {
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive},
        Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ]}
}
fn input<'a>(domain:&'a VoxelDomain,setup:&'a ElectricSetup) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup,contacts:std::slice::from_ref(&FACE),cell_ids:&IDS,
        max_cells:3,max_current_solves:10001,power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn bulk(law:&EquilibriumEnthalpyPhaseCurve,cell:&DepositedCell,t:f64)
    -> Result<EnthalpyTransport<EquilibriumPhaseState>,JouleEnthalpyError> {
    let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
    assert_eq!(phase.temperature_k(),t);
    EnthalpyTransport::try_new(phase,Some(2.0),1.0)
}
// Deliberately no Clone: evidence/state payloads are retained by the one epoch.
#[derive(Debug,PartialEq)]
struct InterfacePoint {phase:EquilibriumPhaseState}
fn response(template:&ElectricalContact,sides:[&EnthalpyTransportSample<EquilibriumPhaseState>;2],threshold:f64)
    -> Result<ElectricalInterfaceState<InterfacePoint>,JouleEnthalpyError> {
    let phase=*sides[1].response.material();
    let point=InterfacePoint {phase};
    if phase.specific_enthalpy_j_kg()>=threshold {Ok(ElectricalInterfaceState::Open(point))}
    else {Ok(ElectricalInterfaceResponse::try_new(point,template.resistance_ohm_m2,
        template.heat_fraction_to_first)?.into())}
}
fn run<'a>(old:&DepositionSimulation<'a>,law:&EquilibriumEnthalpyPhaseCurve,
    domain:&VoxelDomain,electric:&ElectricSetup,target:f64,dt:f64,threshold:f64)
    -> (DepositionSimulation<'a>,JouleSwitchingReport<EquilibriumPhaseState,(),InterfacePoint>) {
    let (next,(),r)=advance_joule_enthalpy_with_switching_interfaces_observed(old,target,config(dt),
        input(domain,electric),&[],&CancelGate::new(),&(),|_,cell,t|bulk(law,cell,t),
        |_,_|Ok((0.0,())),|_,template,sides|response(template,sides,threshold)).unwrap();
    (next,r)
}

#[test]
fn material_threshold_changes_current_shutoff_at_equal_temperature_without_erasing_heat() {
    let law=curve();let old=initial(&law,[40.0;3],false);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Voltage(5.0));
    // Thresholds lie strictly between sampling states, avoiding a test that
    // depends on the sign of the current solver's last rounding bit.
    for (threshold,seconds) in [(44.0,1.0),(48.25,2.0)] {
        let (next,r)=run(&old,&law,&d,&e,2.5,0.125,threshold);
        let expected_h=40.0+4.25*seconds;
        near(next.cells()[&IDS[1]].specific_h_j_kg,expected_h);
        near(next.temperature_k(IDS[1]).unwrap(),310.0);
        near(next.energy_j()-old.energy_j(),5.0*seconds);
        near(r.transport.joule.deposited_joule_j,5.0*seconds);
        near(r.transport.joule.electrode_energy_j,5.0*seconds);
        assert!(r.endpoint.contacts[0].response.is_open());
        assert!(r.endpoint.contacts[0].applied().is_none());
        let retained=r.endpoint.contacts[0].response.material().phase;
        near(retained.liquid_mass_fraction(),(expected_h-20.0)/100.0);
        assert_eq!(retained.phase_curve_identity(),law.identity());
        assert_eq!(r.endpoint.contacts[0].authored,FACE);
        assert_eq!(r.transport.joule.endpoint.electric.total_joule_w,0.0);
        assert_eq!(next.mass_kg(),old.mass_kg());
        for id in IDS {assert_eq!(next.cells()[&id].geometry,old.cells()[&id].geometry);}
    }
    near(old.energy_j(),120.0);assert_eq!(old.time_s(),0.0);
}

#[test]
fn endpoint_switch_does_not_retroactively_erase_work_and_reclosing_preserves_history() {
    let law=curve();let old=initial(&law,[40.0;3],false);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Voltage(5.0));
    let (hot,r)=run(&old,&law,&d,&e,1.0,0.125,44.0);
    assert!(!r.last_applied.as_ref().unwrap().contacts[0].response.is_open());
    assert!(r.endpoint.contacts[0].response.is_open());
    near(r.transport.joule.deposited_joule_j,5.0);
    let (retained,r)=run(&hot,&law,&d,&e,2.0,0.125,44.0);
    assert_eq!(retained.cells(),hot.cells());near(r.transport.joule.deposited_joule_j,0.0);
    let (direct,_)=run(&old,&law,&d,&e,2.0,0.125,44.0);
    assert_eq!(direct.cells(),retained.cells());
    // A separate explicitly closed constitutive state, not a history reset.
    let (closed,r)=run(&retained,&law,&d,&e,2.25,0.125,f64::INFINITY);
    near(closed.energy_j()-retained.energy_j(),1.25);
    near(closed.cells()[&IDS[1]].specific_h_j_kg,44.25+4.25*0.25);
    near(r.transport.joule.endpoint.electric.electrodes[0].current_a,1.0);
}

#[test]
fn electrical_opening_does_not_remove_thermal_contact_or_its_covered_surface() {
    let law=curve();let old=initial(&law,[0.0,40.0,0.0],true);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Voltage(5.0));
    let (next,r)=run(&old,&law,&d,&e,0.125,0.125,0.0);
    near(next.cells()[&IDS[0]].specific_h_j_kg,1.25);
    near(next.cells()[&IDS[1]].specific_h_j_kg,38.75);
    near(next.energy_j(),old.energy_j());near(r.transport.joule.deposited_joule_j,0.0);
    assert_eq!(next.contact_patch(17),old.contact_patch(17));
    near(next.exposed_area_m2(IDS[0]).unwrap(),5.0);
    assert_eq!(r.endpoint.contacts[0].authored,FACE);
}

#[test]
fn closed_variant_is_the_identical_material_interface_transaction() {
    let law=curve();let old=initial(&law,[40.0;3],false);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Voltage(5.0));
    let (new,r)=run(&old,&law,&d,&e,0.5,0.125,f64::INFINITY);
    let (legacy,(),reference)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.5,config(0.125),
        input(&d,&e),&[],&CancelGate::new(),&(),|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
        |_,template,sides|ElectricalInterfaceResponse::try_new(InterfacePoint {phase:*sides[1].response.material()},
            template.resistance_ohm_m2,template.heat_fraction_to_first)).unwrap();
    assert_eq!(new.cells(),legacy.cells());assert_eq!(r.transport.joule,reference.transport.joule);
    assert_eq!(r.transport.endpoint,reference.transport.endpoint);
    assert_eq!(r.endpoint.contacts[0].applied(),Some(reference.endpoint.contacts[0].applied()));
}

#[derive(Clone,Debug,PartialEq)]
struct Heat {watts:f64,joules:f64,steps:usize,reject:bool}
impl DepositionObserver for Heat {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>) -> Result<(),Self::Error> {Ok(())}
    fn sources(&mut self,_:&DepositionSimulation<'_>,check:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error> {check()?;Ok(BTreeMap::from([(IDS[1],self.watts)]))}
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        check()?;self.joules+=dt*self.watts;self.steps+=1;
        if self.reject {Err(DepositionError::InvalidInput("independent owner refuses"))} else {Ok(())}
    }
}

#[test]
fn independent_heat_and_zero_current_floating_terminal_share_one_work_account() {
    let law=curve();let old=initial(&law,[40.0;3],false);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Current(0.0));
    let heat=Heat {watts:2.0,joules:0.0,steps:0,reject:false};
    let (next,accepted,r)=advance_joule_enthalpy_with_switching_interfaces_observed(&old,0.5,config(0.125),
        input(&d,&e),&[],&CancelGate::new(),&heat,|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
        |_,template,sides|response(template,sides,0.0)).unwrap();
    near(next.energy_j()-old.energy_j(),1.0);near(accepted.joules,1.0);assert_eq!(accepted.steps,4);
    assert!(r.transport.joule.endpoint.electric.electrodes[0].voltage_v.is_nan());
    near(r.transport.joule.electrode_energy_j,0.0);near(r.transport.joule.deposited_joule_j,0.0);
    near(r.transport.joule.thermal.external_heat_j,1.0);near(r.transport.joule.thermal.residual_j,0.0);
    assert_eq!(heat.steps,0);near(heat.joules,0.0);
}

#[test]
fn late_opening_current_incompatibility_support_budget_and_cancellation_roll_back() {
    let law=curve();let old=initial(&law,[40.0;3],false);let d=VoxelDomain::new(3,1,1,1.0).unwrap();
    let current=setup(ElectrodeDrive::Current(1.0));let gate=CancelGate::new();
    let failure=advance_joule_enthalpy_with_switching_interfaces_observed(&old,1.0,config(0.125),
        input(&d,&current),&[],&gate,&(),|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
        |_,template,sides|response(template,sides,44.0));
    assert!(matches!(failure,Err(ObserverCouplingError::First(JouleEnthalpyError::Electric(
        ChtError::InvalidInput {field:"electric.electrodes",..})))));
    let e=setup(ElectrodeDrive::Voltage(5.0));let heat=Heat {watts:0.0,joules:0.0,steps:0,reject:false};
    for cancelled in [false,true] {
        let gate=CancelGate::new();let calls=Cell::new(0);
        let result=advance_joule_enthalpy_with_switching_interfaces_observed(&old,0.125,config(0.125),
            input(&d,&e),&[],&gate,&heat,|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
            |_,template,sides| {
                calls.set(calls.get()+1);
                if calls.get()==2 {
                    if cancelled {gate.request();}
                    else {return Err(JouleEnthalpyError::Invalid("unsupported endpoint interface"));}
                }
                response(template,sides,44.25)
            });
        assert!(result.is_err());assert_eq!(calls.get(),2);
    }
    let request=JouleEnthalpyInput {max_current_solves:1,..input(&d,&e)};
    let result=advance_joule_enthalpy_with_switching_interfaces_observed(&old,0.125,config(0.125),
        request,&[],&gate,&heat,|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
        |_,template,sides|response(template,sides,44.25));
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::CurrentBudget))));
    let rejecting=Heat {reject:true,..heat.clone()};
    assert!(advance_joule_enthalpy_with_switching_interfaces_observed(&old,0.125,config(0.125),
        input(&d,&e),&[],&gate,&rejecting,|_,cell,t|bulk(&law,cell,t),|_,_|Ok((0.0,())),
        |_,template,sides|response(template,sides,0.0)).is_err());
    assert_eq!(old.time_s(),0.0);near(old.energy_j(),120.0);assert_eq!(heat.steps,0);assert_eq!(rejecting.steps,0);
}

#[test]
fn state_triggered_switch_overshoot_is_bounded_by_the_explicit_sampling_grid() {
    let law=curve();let old=initial(&law,[40.0;3],false);
    let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(ElectrodeDrive::Voltage(5.0));
    let threshold=45.1;
    for dt in [0.25,0.125,0.0625,0.03125] {
        let (next,r)=run(&old,&law,&d,&e,2.0,dt,threshold);
        let h=next.cells()[&IDS[1]].specific_h_j_kg;
        assert!(h>=threshold && h-threshold<=4.25*dt+1e-10);
        assert!(r.endpoint.contacts[0].response.is_open());
        near(next.energy_j()-old.energy_j(),r.transport.joule.electrode_energy_j);
    }
}
