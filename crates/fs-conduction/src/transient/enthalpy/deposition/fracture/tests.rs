use super::*;
use fs_material::cohesive::BilinearCohesive;
use crate::transient::deposition::{Boundary, CellGeometry, DepositedCell,
    DepositionBudget, MaterialId, ThermalLaw, TimedEvent};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static HEAT: Sensible = Sensible;
fn law() -> BilinearCohesive { BilinearCohesive::new(1e9,2e5,100.0,2e9,0.01).unwrap() }
fn policy() -> ContactPolicy {
    ContactPolicy { open_gap_m: 0.0001, reclose_gap_m: 0.0,
        reclosed_conductance_w_k: 0.05, maximum_slip_m: 0.005 }
}
fn link(id: u64) -> ThermalLink {
    ThermalLink { id, a: CellId(1), b: CellId(2), area_m2: 0.0002, conductance_w_k: 0.2 }
}
fn birth(id: u64, temperature: f64, htc: f64) -> TimedEvent {
    TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
        mass_kg: 0.001, specific_h_j_kg: 2000.0*(temperature-300.0),
        geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],
            end_m: [0.01,0.0,id as f64*0.001], volume_m3: 1e-6, surface_area_m2: 0.001 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: htc, reservoir: None, source_w: 0.0 },
    }) }
}
fn thermal(links: &[ThermalLink], htc: f64, extra: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    let mut events = vec![birth(1,500.0,htc),birth(2,300.0,htc)];
    events.extend(links.iter().map(|l| TimedEvent { time_s: 0.0,event: Event::Contact(*l) }));
    events.extend(extra);
    DepositionSimulation::new(events,&[(MaterialId(1),&HEAT)],DepositionBudget {
        max_cells: 2,max_links: 4,max_events: 16,max_substeps: 10000,
    }).unwrap().advance_to(0.0,config(),|| Ok(())).unwrap().0
}
fn config() -> StepConfig { StepConfig { max_dt_s: 0.01,energy_tolerance_j: 1e-8 } }
fn binding(law: &BilinearCohesive, reference: ThermalLink) -> Interface<'_> {
    Interface { reference,point: law.bonded_point(),closed: true,policy: policy() }
}
fn coupled(law: &BilinearCohesive, htc: f64) -> FractureCooling<'static,'_> {
    FractureCooling::new(thermal(&[link(7)],htc,vec![]),vec![binding(law,link(7))],1,|| Ok(())).unwrap()
}
fn near(a: f64,b: f64,tol: f64) { assert!((a-b).abs() < tol,"{a:e} != {b:e}"); }

// G1: without ambient loss, contact exchange has an independent two-body limit.
#[test]
fn failure_opens_actual_contact_and_changes_subsequent_temperatures() {
    let law = law();
    let state = coupled(&law,0.0);
    let (held,audit) = state.advance_to(1.0,config(),|| Ok(())).unwrap();
    near(held.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-0.2_f64).exp(),0.02);
    near(audit.residual_j,0.0,1e-8);
    let broken = law.bonded_point().trial([law.failure_separation_m(),0.0,0.0],1e-8).unwrap().point;
    let (open,report) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    assert_eq!(report.opened,1);
    near(report.fracture_dissipation_j,0.02,1e-12);
    assert!(open.thermal().links().is_empty());
    assert_eq!(open.thermal().cells(),state.thermal().cells());
    near(open.thermal().exposed_area_m2(CellId(1)).unwrap(),0.001,1e-15);
    let (cooled,audit) = open.advance_to(1.0,config(),|| Ok(())).unwrap();
    assert_eq!(cooled.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(cooled.thermal().temperature_k(CellId(2)).unwrap(),300.0);
    assert_eq!(audit.external_heat_j,0.0);
    assert!(cooled.interfaces()[&7].point.failed());
    assert_eq!(state.thermal().links()[&7],link(7));
}

#[test]
fn reclosure_restores_separate_conductance_but_never_tensile_strength() {
    let law = law();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (open,_) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let compressed = broken.trial([-0.00001,0.0,0.0],1e-8).unwrap().point;
    let (closed,r) = open.synchronize(0.0,&[(7,compressed)],1e-8,|| Ok(())).unwrap();
    assert_eq!(r.reclosed,1);
    assert_eq!(r.fracture_dissipation_j,0.0);
    assert_eq!(closed.thermal().links()[&7].conductance_w_k,0.05);
    assert_eq!(closed.interfaces()[&7].point.damage(),1.0);
    let inside_band = compressed.trial([0.00005,0.0,0.0],1e-8).unwrap().point;
    assert_eq!(inside_band.traction_pa(),[0.0;3]);
    let (still_closed,_) = closed.synchronize(0.0,&[(7,inside_band)],1e-8,|| Ok(())).unwrap();
    assert!(still_closed.interfaces()[&7].closed);
    let (heated,_) = still_closed.advance_to(1.0,config(),|| Ok(())).unwrap();
    near(heated.thermal().temperature_k(CellId(1)).unwrap(),400.0+100.0*(-0.05_f64).exp(),0.002);
    // The same gap stays open when approached from the open branch.
    let approach = broken.trial([0.00005,0.0,0.0],1e-8).unwrap().point;
    let (still_open,_) = open.synchronize(0.0,&[(7,approach)],1e-8,|| Ok(())).unwrap();
    assert!(!still_open.interfaces()[&7].closed);
}

#[test]
fn shear_failure_and_partial_damage_do_not_invent_missing_contact_area() {
    let law = law();
    let state = coupled(&law,0.0);
    let partial = law.bonded_point().trial([0.0005,0.0,0.0],1e-8).unwrap().point;
    let (damaged,_) = state.synchronize(0.0,&[(7,partial)],1e-8,|| Ok(())).unwrap();
    assert!(partial.damage() > 0.0 && !partial.failed());
    assert_eq!(damaged.thermal().links()[&7],link(7));
    let shear = law.bonded_point().trial([0.0,0.001,0.0],1e-8).unwrap().point;
    let (slipped,r) = state.synchronize(0.0,&[(7,shear)],1e-8,|| Ok(())).unwrap();
    assert_eq!(r.changed_conductances,1);
    assert_eq!(r.opened,0);
    assert_eq!(slipped.thermal().links()[&7].area_m2,link(7).area_m2);
    assert_eq!(slipped.thermal().links()[&7].conductance_w_k,0.05);
}

#[test]
fn inconsistent_history_law_ids_time_and_slip_refuse_without_mutation() {
    let law = law();
    let other = law.clone();
    let state = coupled(&law,0.0);
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (open,_) = state.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let excessive_slip = broken.trial([0.0,0.006,0.0],1e-8).unwrap().point;
    for endpoints in [vec![],vec![(8,broken)],vec![(7,law.bonded_point())],
        vec![(7,other.bonded_point())],vec![(7,excessive_slip)]] {
        assert!(open.synchronize(0.0,&endpoints,1e-8,|| Ok(())).is_err());
    }
    assert!(open.synchronize(1.0,&[(7,broken)],1e-8,|| Ok(())).is_err());
    // A larger excursion unseen by this owner cannot be smuggled in at unloading.
    let unload = broken.trial([0.0003,0.0,0.0],1e-8).unwrap().point;
    assert!(state.synchronize(0.0,&[(7,unload)],1e-8,|| Ok(())).is_err());
    assert!(open.interfaces()[&7].point.failed());
    assert!(open.thermal().links().is_empty());
    assert_eq!(open.thermal().energy_j(),state.thermal().energy_j());
}

#[test]
fn cancelled_contact_replacement_discards_the_entire_mechanical_proposal() {
    let law = law();
    let state = coupled(&law,0.0);
    let point = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    for allowed in [0,2,4,6,8] {
        let mut calls = 0;
        assert!(state.synchronize(0.0,&[(7,point)],1e-8,|| {
            calls += 1;
            if calls > allowed { Err(DepositionError::Cancelled) } else { Ok(()) }
        }).is_err());
        assert_eq!(state.interfaces()[&7].point.damage(),0.0);
        assert_eq!(state.thermal().links()[&7],link(7));
    }
}

#[test]
fn same_time_external_changes_cannot_silently_override_mechanical_ownership() {
    let law = law();
    let events = vec![TimedEvent { time_s: 0.5,event: Event::Separation(7) },
        TimedEvent { time_s: 0.5,event: Event::Contact(link(7)) }];
    let state = FractureCooling::new(thermal(&[link(7)],0.0,events),
        vec![binding(&law,link(7))],1,|| Ok(())).unwrap();
    assert!(state.advance_to(1.0,config(),|| Ok(())).is_err());
    assert_eq!(state.thermal().time_s(),0.0);
    assert_eq!(state.thermal().temperature_k(CellId(1)).unwrap(),500.0);
    assert_eq!(state.interfaces()[&7].point.damage(),0.0);
}

#[test]
fn checked_contact_updates_refuse_stale_duplicate_and_invalid_area_batches() {
    let state = thermal(&[link(7)],0.0,vec![]);
    let open = ContactUpdate { id: 7,expected: Some(link(7)),replacement: None };
    assert!(state.replace_contacts(&[open,open],|| Ok(())).is_err());
    assert!(state.replace_contacts(&[ContactUpdate { expected: None,..open }],|| Ok(())).is_err());
    let huge = ThermalLink { area_m2: 0.1,..link(7) };
    assert!(state.replace_contacts(&[ContactUpdate { replacement: Some(huge),..open }],|| Ok(())).is_err());
    assert_eq!(state.links()[&7],link(7));
    let opened = state.replace_contacts(&[open],|| Ok(())).unwrap();
    assert_eq!(state.cells(),opened.cells());
    assert_eq!(state.time_s(),opened.time_s());
    assert!(opened.replace_contacts(&[open],|| Ok(())).is_err());
}

#[test]
fn area_partition_and_same_grid_restart_preserve_heat_and_fracture_work() {
    let law = law();
    let whole = coupled(&law,100.0);
    let a = ThermalLink { area_m2: 0.0001,conductance_w_k: 0.1,..link(7) };
    let b = ThermalLink { id: 8,..a };
    let mut ia = binding(&law,a);
    ia.policy.reclosed_conductance_w_k = 0.025;
    let split = FractureCooling::new(thermal(&[a,b],100.0,vec![]),
        vec![ia,Interface { reference: b,..ia }],2,|| Ok(())).unwrap();
    let broken = law.bonded_point().trial([0.001,0.0,0.0],1e-8).unwrap().point;
    let (whole,wr) = whole.synchronize(0.0,&[(7,broken)],1e-8,|| Ok(())).unwrap();
    let (split,sr) = split.synchronize(0.0,&[(8,broken),(7,broken)],1e-8,|| Ok(())).unwrap();
    assert_eq!(wr.fracture_dissipation_j,sr.fracture_dissipation_j);
    let cfg = StepConfig { max_dt_s: 0.125,..config() };
    let (whole,_) = whole.advance_to(2.0,cfg,|| Ok(())).unwrap();
    let (prefix,_) = split.advance_to(1.0,cfg,|| Ok(())).unwrap();
    let (resumed,_) = prefix.clone().advance_to(2.0,cfg,|| Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(),resumed.thermal().cells());
    assert_eq!(resumed.interfaces()[&7].point,broken);
    assert_eq!(resumed.interfaces()[&8].point,broken);
}
