#![cfg(feature = "enthalpy-coupling")]
//! Original current/heat solvers with state-dependent ordered electrical contacts.
use std::cell::Cell;
use std::collections::BTreeMap;
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_conduction::observers::ObserverCouplingError;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::*;
use fs_lbm::conjugate::enthalpy::transport::*;
use fs_lbm::conjugate::enthalpy::transport::interfaces::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};

const IDS: [CellId; 3] = [CellId(10), CellId(20), CellId(30)];
fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xC3;32]),
        [(0.0,300.0,0.0),(20.0,310.0,0.0),(120.0,310.0,1.0),(1120.0,810.0,1.0)]
        .map(|(h,t,f)| EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.0}).to_vec()).unwrap()
}
fn config() -> StepConfig { StepConfig {max_dt_s:1.0/64.0,energy_tolerance_j:1e-8} }
fn initial(law: &EquilibriumEnthalpyPhaseCurve, h: f64) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = IDS.iter().enumerate().map(|(c,&id)| TimedEvent {time_s:0.0,
        event:Event::Birth(DepositedCell { id,segment_id:c as u64,body_id:1,material_id:MaterialId(7),
            mass_kg:1.0,specific_h_j_kg:h,geometry:CellGeometry {
                start_m:[c as f64,0.5,0.5],end_m:[c as f64+1.0,0.5,0.5],volume_m3:1.0,surface_area_m2:6.0},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0}})
    }).collect();
    for c in 0..2 {events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:c as u64,a:IDS[c],b:IDS[c+1],area_m2:1.0,conductance_w_k:1.0})});}
    DepositionSimulation::new(events,&[(MaterialId(7),law as &dyn ThermalLaw)],DepositionBudget {
        max_cells:3,max_links:2,max_events:5,max_substeps:10000}).unwrap()
        .advance_to(0.0,config(),||Ok(())).unwrap().0
}
fn setup(voltage: f64) -> ElectricSetup { ElectricSetup {conductors:vec![],electrodes:vec![
    Electrode {lo:[0,0,0],hi:[1,1,1],drive:ElectrodeDrive::Voltage(voltage)},
    Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)}]}}
fn contact() -> ElectricalContact {
    ElectricalContact {cells:[0,1],resistance_ohm_m2:3.0,heat_fraction_to_first:0.25}
}
fn input<'a>(d:&'a VoxelDomain,e:&'a ElectricSetup,c:&'a [ElectricalContact]) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain:d,setup:e,contacts:c,cell_ids:&IDS,max_cells:3,max_current_solves:10001,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn bulk(_:usize,_:&DepositedCell,_:f64)->Result<EnthalpyTransport<()>,JouleEnthalpyError> {
    EnthalpyTransport::try_new((),Some(2.0),0.0)
}
fn link<F>(_:&ContactPatch,s:[&EnthalpyTransportSample<F>;2])->Result<(f64,()),JouleEnthalpyError> {
    let [a,b]=s.map(|s|s.response.conductivity_w_m_k());
    Ok((if a==0.0 || b==0.0 {0.0} else {1.0/(0.5/a+0.5/b)},()))
}
fn fixed<F>(_:usize,c:&ElectricalContact,_:[&EnthalpyTransportSample<F>;2])
    ->Result<ElectricalInterfaceResponse<()>,JouleEnthalpyError> {
    ElectricalInterfaceResponse::try_new((),c.resistance_ohm_m2,c.heat_fraction_to_first)
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-7,"{a} != {b}");}

#[test]
fn constant_interface_response_preserves_the_original_physical_path_exactly() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();
    let e=setup(5.0);let c=[contact()];let request=input(&d,&e,&c);let gate=CancelGate::new();
    let (a,ra)=advance_joule_enthalpy_with_material_transport(&old,0.125,config(),request,&[0,1],&gate,bulk,link).unwrap();
    let (b,(),rb)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),request,&[0,1],&gate,&(),bulk,link,fixed).unwrap();
    assert_eq!(a.cells(),b.cells());assert_eq!(a.links(),b.links());assert_eq!(ra.joule,rb.transport.joule);
    assert_eq!(ra.endpoint,rb.transport.endpoint);assert_eq!(ra.last_applied,rb.transport.last_applied);
    for (cell,watts) in IDS.into_iter().zip([0.75,4.25,0.0]) {near(b.cells()[&cell].specific_h_j_kg,0.125*watts);}
    assert_eq!(rb.endpoint.contacts[0].applied(),c[0]);assert_eq!(c[0],contact());
}

#[test]
fn phase_dependent_interfaces_reuse_nonclone_bulk_samples_and_change_current_at_fixed_kelvin() {
    struct Material { phase:f64, serial:usize }
    struct Interface { serials:[usize;2] }
    let law=curve();let old=initial(&law,40.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(10.0);
    let c=[contact()];let q=Cell::new(0);let iq=Cell::new(0);let lq=Cell::new(0);
    let (next,(),r)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),input(&d,&e,&c),&[0,1],&CancelGate::new(),&(),
        |_,cell,t| {
            let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();near(t,phase.temperature_k());
            let serial=q.get();q.set(serial+1);
            EnthalpyTransport::try_new(Material {phase:phase.liquid_mass_fraction(),serial},Some(2.0),0.0)
        },|p,s| {
            lq.set(lq.get()+1);assert_eq!(s[0].response.material().serial/3,s[1].response.material().serial/3);
            link(p,s)
        },|index,template,s| {
            iq.set(iq.get()+1);assert_eq!(index,0);assert_eq!(template.cells,[0,1]);
            assert_eq!(s[0].cell.id,IDS[0]);assert_eq!(s[1].cell.id,IDS[1]);
            ElectricalInterfaceResponse::try_new(Interface {serials:s.map(|v|v.response.material().serial)},
                3.0+10.0*s[1].response.material().phase,0.25)
        }).unwrap();
    let j=&r.transport.joule;assert_eq!(q.get(),3*j.current_solves);assert_eq!(iq.get(),j.current_solves);assert_eq!(lq.get(),2*j.current_solves);
    let used=r.last_applied.as_ref().unwrap();let end=&r.endpoint;
    let source=r.transport.last_applied.as_ref().unwrap();
    assert_eq!(used.contacts[0].response.material().serials,[source.cells[0].response.material().serial,source.cells[1].response.material().serial]);
    assert!(end.contacts[0].response.resistance_ohm_m2()>used.contacts[0].response.resistance_ohm_m2());
    assert_eq!(used.time_s,source.time_s);assert_eq!(end.time_s,r.transport.endpoint.time_s);
    near(next.temperature_k(IDS[1]).unwrap(),310.0);near(old.temperature_k(IDS[1]).unwrap(),310.0);
    assert!(j.endpoint.electric.electrodes[0].current_a<j.last_applied.as_ref().unwrap().electric.electrodes[0].current_a);
    near(next.energy_j()-old.energy_j(),j.deposited_joule_j);near(j.electrical_defect_j,0.0);near(j.transfer_defect_j,0.0);
}

#[test]
fn material_only_contact_substitution_and_surface_or_current_reversal_have_physical_effects() {
    let law=curve();let old=initial(&law,40.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let gate=CancelGate::new();
    let run=|r0:f64,reverse_surface:bool,voltage:f64| {
        let e=setup(voltage);let c=[ElectricalContact {cells:if reverse_surface {[1,0]} else {[0,1]},..contact()}];
        advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),input(&d,&e,&c),&[0,1],&gate,&(),bulk,link,
            |_,_,s| {
                let fraction=law.state_at_specific_enthalpy(s[usize::from(!reverse_surface)].cell.specific_h_j_kg).unwrap().liquid_mass_fraction();
                ElectricalInterfaceResponse::try_new((),r0+fraction,if reverse_surface {0.75} else {0.25})
            }).unwrap()
    };
    let (a,(),ra)=run(3.0,false,10.0);let (b,(),rb)=run(8.0,false,10.0);
    assert!(a.cells()[&IDS[1]].specific_h_j_kg>b.cells()[&IDS[1]].specific_h_j_kg);
    assert!(ra.transport.joule.endpoint.electric.electrodes[0].current_a>rb.transport.joule.endpoint.electric.electrodes[0].current_a);
    for (surface,voltage) in [(true,10.0),(false,-10.0),(true,-10.0)] {
        let (other,(),r)=run(3.0,surface,voltage);
        for id in IDS {near(a.cells()[&id].specific_h_j_kg,other.cells()[&id].specific_h_j_kg);}
        near(ra.transport.joule.deposited_joule_j,r.transport.joule.deposited_joule_j);
    }
}

#[derive(Clone,Debug,PartialEq)]
struct Heater {watts:f64,joules:f64,steps:usize,refuse:bool}
impl DepositionObserver for Heater {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn sources(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error>{
        Ok(BTreeMap::from([(IDS[1],self.watts)]))
    }
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{
        self.steps+=1;self.joules+=dt*self.watts;
        if self.refuse {Err(DepositionError::MaterialDomain)} else {Ok(())}
    }
}

#[test]
fn shared_cfl_and_external_heat_preserve_each_account_and_source_off_memory() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);let c=[contact()];let gate=CancelGate::new();
    let h=Heater {watts:2.0,joules:0.0,steps:0,refuse:false};
    let (hot,history,r)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),input(&d,&e,&c),&[0,1],&gate,&h,
        |_,_,_|EnthalpyTransport::try_new((),Some(2.0),1000.0),link,fixed).unwrap();
    let j=&r.transport.joule;assert!(history.steps>8);assert_eq!(history.steps,j.thermal.substeps);assert_eq!(j.current_solves,history.steps+1);
    near(history.joules,0.25);near(j.deposited_joule_j,0.625);near(j.thermal.external_heat_j,0.875);near(hot.energy_j(),0.875);
    assert_eq!(h.joules,0.0);assert_eq!(old.energy_j(),0.0);
    let off=setup(0.0);let h=Heater {watts:0.0,..history};
    let (cold,history,r)=advance_joule_enthalpy_with_material_interfaces_observed(&hot,0.25,config(),input(&d,&off,&c),&[0,1],&gate,&h,
        |_,_,_|EnthalpyTransport::try_new((),Some(2.0),1000.0),link,fixed).unwrap();
    near(cold.energy_j(),hot.energy_j());near(history.joules,0.25);near(r.transport.joule.deposited_joule_j,0.0);
}

#[test]
fn late_interface_or_observer_refusal_budget_and_cancellation_return_no_candidate() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);let c=[contact()];let gate=CancelGate::new();
    let h=Heater {watts:1.0,joules:0.0,steps:0,refuse:false};let config=StepConfig {max_dt_s:0.125,..config()};
    let result=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config,input(&d,&e,&c),&[0,1],&gate,&h,bulk,link,
        |i,c,s| {if s[1].cell.specific_h_j_kg>0.0 {Err(JouleEnthalpyError::Invalid("unsupported endpoint"))} else {fixed(i,c,s)}});
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::Invalid("unsupported endpoint")))));
    let one=JouleEnthalpyInput {max_current_solves:1,..input(&d,&e,&c)};
    assert!(matches!(advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config,one,&[0,1],&gate,&h,bulk,link,fixed),
        Err(ObserverCouplingError::First(JouleEnthalpyError::CurrentBudget))));
    let bad=Heater {refuse:true,..h.clone()};
    assert!(matches!(advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config,input(&d,&e,&c),&[0,1],&gate,&bad,bulk,link,fixed),
        Err(ObserverCouplingError::Second(DepositionError::MaterialDomain))));
    let cancelled=CancelGate::new();
    assert!(advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config,input(&d,&e,&c),&[0,1],&cancelled,&h,bulk,link,
        |i,c,s| {cancelled.request();fixed(i,c,s)}).is_err());
    assert_eq!(old.energy_j(),0.0);assert_eq!(old.time_s(),0.0);assert_eq!(h.joules,0.0);assert_eq!(bad.steps,0);
}

#[test]
fn malformed_contacts_refuse_before_sampling_and_zero_duration_only_observes() {
    let law=curve();let old=initial(&law,0.0);let d=VoxelDomain::new(3,1,1,1.0).unwrap();let e=setup(5.0);let gate=CancelGate::new();
    for cells in [[0,2],[1,1],[3,0]] {
        let c=[ElectricalContact {cells,..contact()}];
        assert!(advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),input(&d,&e,&c),&[0,1],&gate,&(),
            |_,_,_|->Result<EnthalpyTransport<()>,JouleEnthalpyError>{panic!("topology before samples")},link,fixed).is_err());
    }
    let c=[contact(),ElectricalContact {cells:[1,0],..contact()}];
    assert!(advance_joule_enthalpy_with_material_interfaces_observed(&old,0.125,config(),input(&d,&e,&c),&[0,1],&gate,&(),bulk,link,fixed).is_err());
    let c=[contact()];
    let (next,(),r)=advance_joule_enthalpy_with_material_interfaces_observed(&old,0.0,config(),input(&d,&e,&c),&[0,1],&gate,&(),bulk,link,fixed).unwrap();
    assert_eq!(old.cells(),next.cells());assert_eq!(r.transport.joule.current_solves,1);assert!(r.last_applied.is_none());
    assert!(r.transport.last_applied.is_none());near(r.transport.joule.deposited_joule_j,0.0);
    for (resistance,fraction) in [(-1.0,0.5),(f64::NAN,0.5),(f64::INFINITY,0.5),(1.0,-0.1),(1.0,1.1),(1.0,f64::NAN)] {
        assert!(ElectricalInterfaceResponse::try_new((),resistance,fraction).is_err());
    }
}
