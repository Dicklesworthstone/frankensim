//! Shared partners advance once through the original thermal network.
use super::*;
use std::collections::BTreeMap;
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget,
    DepositionError, DepositionObserver, DepositionSimulation, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent};

fn bindings()->Vec<ReservoirCell> {
    [("right",3),("support",77),("left",9)].into_iter().map(|(key,id)|
        ReservoirCell {reservoir:key.into(),cell:CellId(id)}).collect()
}
fn thermal_config(dt:f64)->StepConfig { StepConfig {max_dt_s:dt,energy_tolerance_j:1e-12} }
fn heated(c:&EquilibriumEnthalpyPhaseCurve,chi:Option<f64>)->HeatedImpactSystem<'_> {
    let stores=[("left",0.4),("right",0.1),("support",0.0)].into_iter().map(|(key,h)|
        CaloricReservoir::new(key.into(),c,c.state_at_specific_enthalpy(h).unwrap(),
            0.001,1e-6,0.0006).unwrap()).collect();
    HeatedImpactSystem::new(collision(chi,0.0,true).prepare_analytic().unwrap(),stores,
        vec![route("touch"),route("other")],1e-10).unwrap()
}
fn network<'a>(c:&'a EquilibriumEnthalpyPhaseCurve, extra:Vec<TimedEvent>)->DepositionSimulation<'a> {
    let mut schedule=Vec::new();
    for (id,h) in [(9,0.4),(3,0.1),(77,0.0)] {
        schedule.push(TimedEvent {time_s:0.0,event:Event::Birth(DepositedCell {
            id:CellId(id),segment_id:id,body_id:id,material_id:MaterialId(1),mass_kg:0.001,
            specific_h_j_kg:h,geometry:CellGeometry {start_m:[id as f64,0.0,0.0],
                end_m:[id as f64+0.01,0.0,0.0],volume_m3:1e-6,surface_area_m2:0.0006},
            boundary:Boundary {ambient_k:300.0,htc_w_m2_k:0.0,reservoir:None,source_w:0.0},
        })});
    }
    for (id,a,b) in [(101,9,3),(102,3,77)] {
        schedule.push(TimedEvent {time_s:0.0,event:Event::Contact(ThermalLink {
            id,a:CellId(a),b:CellId(b),area_m2:1e-6,conductance_w_k:0.01,
        })});
    }
    schedule.extend(extra);
    let initial=DepositionSimulation::new(schedule,&[(MaterialId(1),c as &dyn ThermalLaw)],
        DepositionBudget {max_cells:3,max_links:2,max_events:16,max_substeps:16}).unwrap();
    initial.advance_to(0.0,thermal_config(1e-6),||Ok(())).unwrap().0
}
fn unchanged(_: &fs_couple::render::plate::impact::material::MaterialContact,
    _: &[CaloricReservoir<'_>])->Result<Option<SourcePointContact>,ImpactError> {Ok(None)}

#[test]
fn original_network_parity_and_three_partner_simultaneous_refinement() {
    let c=curve(0.02,1.0);let gate=CancelGate::new_clock_free();let mut errors=Vec::new();
    for dt in [2e-6,1e-6,0.5e-6] {
        let mut heated=heated(&c,None);let mut thermal=network(&c,vec![]);
        let mut reference=thermal.clone();let mut observer=();
        for sample in 1..=500 {
            reference=reference.advance_to(sample as f64*config().dt_s,thermal_config(dt),||Ok(())).unwrap().0;
            let frame=heated.step_with_thermal_network(&mut thermal,&mut observer,&bindings(),
                thermal_config(dt),&[0.0;2],&gate,unchanged).unwrap();
            assert_eq!(thermal.cells(),reference.cells());assert_eq!(thermal.time_s(),reference.time_s());
            assert_eq!(frame.network.events,0);assert_eq!(frame.impact.thermal.impact.deposited_heat_j,0.0);
            near(frame.non_contact_heat_j,0.0,1e-18);
            for binding in bindings() {
                let reservoir=heated.reservoirs().iter().find(|r|r.key()==binding.reservoir).unwrap();
                assert_eq!(reservoir.state().specific_enthalpy_j_kg(),thermal.cells()[&binding.cell].specific_h_j_kg);
            }
        }
        // Three equal capacities in a chain: nonuniform eigenmodes decay with
        // rates g/C and 3g/C. No pair may see another pair's partially updated state.
        let time=thermal.time_s();let mean=925.0/3.0;
        let exact=[mean+10.0*(-500.0*time).exp()+(5.0/3.0)*(-1500.0*time).exp(),
            mean-(10.0/3.0)*(-1500.0*time).exp(),
            mean-10.0*(-500.0*time).exp()+(5.0/3.0)*(-1500.0*time).exp()];
        let error=[9,3,77].into_iter().zip(exact).map(|(id,t)|
            (thermal.temperature_k(CellId(id)).unwrap()-t).abs()).fold(0.0,f64::max);
        errors.push(error);
    }
    assert!(errors[0]>1.9*errors[1] && errors[1]>1.9*errors[2]);
}

#[derive(Clone,Default,Debug,PartialEq)]
struct MaterialObserver {steps:usize,source_w:f64,fail_at:Option<usize>}
impl DepositionObserver for MaterialObserver {
    type Error=DepositionError;
    fn event(&mut self,_:&Event,_:&DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn sources(&mut self,_:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<CellId,f64>,Self::Error>{Ok([(CellId(77),self.source_w)].into_iter().collect())}
    fn conductances(&mut self,t:&DepositionSimulation<'_>,_:&mut dyn FnMut()->Result<(),DepositionError>)
        ->Result<BTreeMap<u64,f64>,Self::Error> {
        // Synthetic temperature dependence; no new integration rule.
        let g=0.01*(1.0+0.001*(t.temperature_k(CellId(3))?-300.0));
        Ok([(101,g),(102,g)].into_iter().collect())
    }
    fn step(&mut self,_:f64,_:f64,_:&BTreeMap<CellId,f64>,_:&DepositionSimulation<'_>,
        _:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        self.steps+=1;
        if self.fail_at==Some(self.steps){Err(DepositionError::MaterialDomain)}else{Ok(())}
    }
}

#[test]
fn actual_contact_work_and_other_sources_share_one_network_without_persistent_generated_heat() {
    let c=curve(0.02,1.0);let gate=CancelGate::new_clock_free();
    let mut s=heated(&c,Some(0.4));let mut thermal=network(&c,vec![]);
    let mut observer=MaterialObserver {source_w:0.02,..Default::default()};
    let mut plain=collision(Some(0.4),0.0,true).prepare_analytic().unwrap();
    let initial=plain.stored_energy_j();let mut heat=0.0;let mut supplied=0.0;let mut substeps=0;
    for _ in 0..1800 {
        plain.step(&[0.0;2],&gate).unwrap();
        let frame=s.step_with_thermal_network(&mut thermal,&mut observer,&bindings(),
            thermal_config(5e-7),&[0.0;2],&gate,unchanged).unwrap();
        assert_eq!(plain.state(),s.mechanics().state());
        heat+=frame.impact.thermal.impact.deposited_heat_j;
        supplied+=frame.non_contact_heat_j;substeps+=frame.network.substeps;
        assert!(frame.impact.thermal.impact.combined_energy_residual_j.abs()<1e-9);
        assert!(thermal.cells().values().all(|cell|cell.boundary.source_w==0.0));
        assert!(thermal.links().values().all(|link|link.conductance_w_k==0.01));
    }
    assert!(heat>1e-9);assert_eq!(observer.steps,substeps);assert!(substeps>=7200);
    near(supplied,0.02*thermal.time_s(),1e-15);
    near(s.mechanics().stored_energy_j()+s.reservoirs().iter().map(|r|r.storage_change_j()).sum::<f64>()-supplied,initial,1e-9);
    // Remove the generated source. The accepted heat and network state remain.
    let frame=s.step_with_thermal_network(&mut thermal,&mut (),&bindings(),
        thermal_config(5e-7),&[0.0;2],&gate,unchanged).unwrap();
    assert_eq!(frame.network.external_heat_j,0.0);
}

#[test]
fn late_network_refusal_rolls_back_material_refresh_motion_all_cells_and_observer_history() {
    let c=curve(0.02,1.0);let gate=CancelGate::new_clock_free();let mut s=heated(&c,Some(0.4));
    let mut thermal=network(&c,vec![]);let cells=thermal.cells().clone();let stores=s.reservoirs().to_vec();
    let motion=s.mechanics().state().to_vec();let sources=s.mechanics().material_contacts().to_vec();
    let mut observer=MaterialObserver {fail_at:Some(2),..Default::default()};let saved=observer.clone();
    assert!(s.step_with_thermal_network(&mut thermal,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&gate,
        |_,_|Ok(Some(crate::common::source(3e7,Some(0.4))))).is_err());
    assert_eq!(observer,saved);assert_eq!(thermal.cells(),&cells);assert_eq!(thermal.time_s(),0.0);
    assert_eq!(s.reservoirs(),stores);assert_eq!(s.mechanics().state(),motion);
    assert_eq!(s.mechanics().material_contacts(),sources);assert_eq!(s.mechanics().samples(),0);
    observer.fail_at=None;
    s.step_with_thermal_network(&mut thermal,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).unwrap();
    assert_eq!(s.mechanics().samples(),1);assert_eq!(thermal.time_s(),config().dt_s);
}

#[test]
fn stale_chart_geometry_mapping_clock_pending_events_and_cancellation_do_not_publish() {
    let c=curve(0.02,1.0);let foreign=c.clone();let gate=CancelGate::new_clock_free();
    let mut s=heated(&c,Some(0.4));let stores=s.reservoirs().to_vec();let mut observer=();
    let mut wrong_chart=network(&foreign,vec![]);
    assert!(s.step_with_thermal_network(&mut wrong_chart,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    let mut thermal=network(&c,vec![]);let mut wrong=bindings();wrong[1].cell=wrong[0].cell;
    assert!(s.step_with_thermal_network(&mut thermal,&mut observer,&wrong,thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    let mut stale=thermal.advance_to(1e-6,thermal_config(1e-6),||Ok(())).unwrap().0;
    assert!(s.step_with_thermal_network(&mut stale,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    let mut pending=network(&c,vec![TimedEvent {time_s:config().dt_s,event:Event::Separation(101)}]);
    assert!(s.step_with_thermal_network(&mut pending,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&gate,unchanged).is_err());
    assert_eq!(pending.time_s(),0.0);assert_eq!(pending.links().len(),2);
    let stopped=CancelGate::new_clock_free();stopped.request();
    assert!(s.step_with_thermal_network(&mut thermal,&mut observer,&bindings(),thermal_config(5e-7),&[0.0;2],&stopped,unchanged).is_err());
    assert_eq!(thermal.time_s(),0.0);assert_eq!(s.reservoirs(),stores);assert_eq!(s.mechanics().samples(),0);
}

#[path="network/contact.rs"]
mod contact;
