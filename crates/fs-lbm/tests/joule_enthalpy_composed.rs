#![cfg(feature = "enthalpy-coupling")]
//! Actual electrical solve + existing phase chart + independently owned heating.
use std::collections::BTreeMap;
use fs_blake3::ContentHash;
use fs_conduction::transient::deposition::*;
use fs_conduction::observers::ObserverCouplingError;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyInput, JouleEnthalpyError};
use fs_lbm::conjugate::enthalpy::transport::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

fn curve(latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0xD9; 32]),
        [(0.0, 300.0, 0.0), (20.0, 310.0, 0.0), (20.0+latent, 310.0, 1.0), (1020.0+latent, 560.0, 1.0)]
        .map(|(h,t,f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg:h, temperature_k:t,
            liquid_mass_fraction:f, bulk_density_kg_m3:1.0 }).to_vec()).unwrap()
}
fn config() -> StepConfig { StepConfig { max_dt_s:0.125, energy_tolerance_j:1e-8 } }
fn initial(law: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events:Vec<_> = (0..3).map(|c| TimedEvent { time_s:0.0, event:Event::Birth(DepositedCell {
        id:CellId(c), segment_id:c, body_id:1, material_id:MaterialId(1), mass_kg:1.0, specific_h_j_kg:0.0,
        geometry:CellGeometry { start_m:[c as f64,0.5,0.5], end_m:[c as f64+1.0,0.5,0.5], volume_m3:1.0,surface_area_m2:6.0 },
        boundary:Boundary { ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0 },
    }) }).collect();
    events.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
        id:7,a:CellId(0),b:CellId(1),area_m2:1.0,conductance_w_k:0.0,
    })});
    DepositionSimulation::new(events, &[(MaterialId(1),law as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:1,max_events:4,max_substeps:10000}).unwrap()
        .advance_to(0.0,config(),||Ok(())).unwrap().0
}
fn setup() -> ElectricSetup {
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo:[0,0,0],hi:[1,1,1],drive:ElectrodeDrive::Current(4.0)},
        Electrode {lo:[2,0,0],hi:[3,1,1],drive:ElectrodeDrive::Voltage(0.0)},
    ]}
}
fn input<'a>(domain:&'a VoxelDomain, setup:&'a ElectricSetup, ids:&'a [CellId]) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {domain,setup,contacts:&[],cell_ids:ids,max_cells:3,max_current_solves:10001,
        power_tolerance_w:1e-8,energy_tolerance_j:1e-8}
}
fn resolve(law:&EquilibriumEnthalpyPhaseCurve, cell:&DepositedCell) -> Result<EnthalpyTransport<EquilibriumPhaseState>,JouleEnthalpyError> {
    let phase=law.state_at_specific_enthalpy(cell.specific_h_j_kg).unwrap();
    EnthalpyTransport::try_new(phase,Some(5.0),2.0)
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-7,"{a} != {b}");}

#[derive(Clone, Debug, Default, PartialEq)]
struct Heater { joules:f64, steps:usize, link:Option<f64>, fail_step:Option<usize> }
impl DepositionObserver for Heater {
    type Error=DepositionError;
    fn sources(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>) -> Result<BTreeMap<CellId,f64>,Self::Error> {
        Ok(BTreeMap::from([(CellId(1),20.0)]))
    }
    fn conductances(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>) -> Result<BTreeMap<u64,f64>,Self::Error> {
        Ok(self.link.map(|g|(7,g)).into_iter().collect())
    }
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self,_:f64,dt:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{
        self.joules+=20.0*dt; self.steps+=1;
        if self.fail_step==Some(self.steps) {return Err(DepositionError::InvalidInput("extra physical owner rejected endpoint"));}
        Ok(())
    }
}

#[test]
fn independent_work_and_joule_heat_change_one_material_phase_history() {
    for latent in [100.0,200.0] {
        let law=curve(latent);let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let setup=setup();
        let ids=[CellId(0),CellId(1),CellId(2)];let extra=Heater::default();
        let (next,accepted,report)=advance_joule_enthalpy_with_material_transport_observed(
            &old,1.0,config(),input(&domain,&setup,&ids),&[],&CancelGate::new(),&extra,
            |_,cell,_|resolve(&law,cell),|_,_|Ok((0.0,()))).unwrap();
        near(next.energy_j(),100.0);near(report.joule.deposited_joule_j,80.0);near(accepted.joules,20.0);
        near(report.joule.thermal.external_heat_j,100.0);
        near(report.endpoint.cells[1].response.material().liquid_mass_fraction(),80.0/latent);
        near(next.temperature_k(ids[1]).unwrap(),310.0);
        assert_eq!(accepted.steps,report.joule.thermal.substeps);
        assert_eq!(report.joule.current_solves,accepted.steps+1);
        assert_eq!(extra,Heater::default());near(old.energy_j(),0.0);
        assert!(next.cells().values().all(|c|c.boundary.source_w==0.0));
    }
}

#[test]
fn no_op_preserves_the_original_material_transport_solution_and_receipts() {
    let law=curve(100.0);let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let setup=setup();
    let ids=[CellId(0),CellId(1),CellId(2)];let gate=CancelGate::new();let request=input(&domain,&setup,&ids);
    let (reference,report)=advance_joule_enthalpy_with_material_transport(&old,0.5,config(),request,&[7],&gate,
        |_,cell,_|resolve(&law,cell),|_,_|Ok((10.0,()))).unwrap();
    let (next,(),observed)=advance_joule_enthalpy_with_material_transport_observed(&old,0.5,config(),request,&[7],&gate,&(),
        |_,cell,_|resolve(&law,cell),|_,_|Ok((10.0,()))).unwrap();
    assert_eq!(next.cells(),reference.cells());assert_eq!(observed.joule,report.joule);
    assert_eq!(observed.endpoint,report.endpoint);assert_eq!(observed.last_applied,report.last_applied);
}

#[test]
fn endpoint_material_failure_or_extra_history_failure_publishes_neither_candidate() {
    let law=curve(100.0);let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let setup=setup();
    let ids=[CellId(0),CellId(1),CellId(2)];let gate=CancelGate::new();let request=input(&domain,&setup,&ids);
    let extra=Heater::default();
    let result=advance_joule_enthalpy_with_material_transport_observed(&old,1.0,config(),request,&[],&gate,&extra,
        |_,cell,_| if cell.specific_h_j_kg>99.0 {Err(JouleEnthalpyError::Invalid("source endpoint unsupported"))} else {resolve(&law,cell)},
        |_,_|Ok((0.0,())));
    assert!(matches!(result,Err(ObserverCouplingError::First(JouleEnthalpyError::Invalid(_)))));
    let failing=Heater {fail_step:Some(3),..Heater::default()};
    let result=advance_joule_enthalpy_with_material_transport_observed(&old,1.0,config(),request,&[],&gate,&failing,
        |_,cell,_|resolve(&law,cell),|_,_|Ok((0.0,())));
    assert!(matches!(result,Err(ObserverCouplingError::Second(DepositionError::InvalidInput(_)))));
    near(old.energy_j(),0.0);assert_eq!(old.time_s(),0.0);near(extra.joules,0.0);near(failing.joules,0.0);
}

#[test]
fn extra_transport_cannot_replace_the_shared_material_epoch_link() {
    let law=curve(100.0);let old=initial(&law);let domain=VoxelDomain::new(3,1,1,1.0).unwrap();let setup=setup();
    let ids=[CellId(0),CellId(1),CellId(2)];let extra=Heater {link:Some(10.0),..Heater::default()};
    let result=advance_joule_enthalpy_with_material_transport_observed(&old,1.0,config(),input(&domain,&setup,&ids),&[7],&CancelGate::new(),&extra,
        |_,cell,_|resolve(&law,cell),|_,_|Ok((10.0,())));
    assert!(matches!(result,Err(ObserverCouplingError::ConflictingConductance(7))));
    near(old.energy_j(),0.0);assert_eq!(old.links()[&7].conductance_w_k,0.0);
}
