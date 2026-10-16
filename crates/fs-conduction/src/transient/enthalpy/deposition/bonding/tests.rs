use super::*;
use crate::transient::deposition::{Boundary, CellGeometry, DepositedCell};
use fs_material::healing::{MobilityKnot, MobilityTable};
use fs_material::healing::strength::{StrengthKnot, StrengthMap};

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64,DepositionError> {
        if !(0.0..=600_000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0+h/2000.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.0005 }
}
static THERMAL: Sensible = Sensible;
fn mobility() -> MobilityTable {
    MobilityTable::new("test-mobility".into(),vec![
        MobilityKnot { temperature_k: 300.0, rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 400.0, rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 500.0, rate_per_s: 1.0 },
        MobilityKnot { temperature_k: 600.0, rate_per_s: 3.0 },
    ]).unwrap()
}
fn map() -> StrengthMap {
    StrengthMap::new("test-strength".into(),"test-mobility".into(),vec![
        StrengthKnot { dose: 0.0, peak_traction_pa: 0.0, fracture_energy_j_m2: 0.0 },
        StrengthKnot { dose: 10.0, peak_traction_pa: 2e5, fracture_energy_j_m2: 100.0 },
    ],1e9,2e9,0.01).unwrap()
}
fn boundary(htc: f64,source: f64) -> Boundary {
    Boundary { ambient_k: 300.0, htc_w_m2_k: htc, source_w: source, reservoir: None }
}
fn birth(id: u64,temperature: f64,htc: f64) -> TimedEvent {
    event(0.0,Event::Birth(DepositedCell { id: CellId(id), segment_id: id, body_id: 1,
        material_id: MaterialId(1), mass_kg: 0.001, specific_h_j_kg: 2000.0*(temperature-300.0),
        geometry: CellGeometry { start_m: [0.0,0.0,id as f64*0.001],
            end_m: [0.01,0.0,id as f64*0.001], volume_m3: 1e-6, surface_area_m2: 0.0012 },
        boundary: boundary(htc,0.0) }))
}
fn event(time_s: f64,event: Event) -> TimedEvent { TimedEvent { time_s,event } }
fn contact(id: u64,area: f64,g: f64) -> ThermalLink {
    ThermalLink { id,a: CellId(1),b: CellId(2),area_m2: area,conductance_w_k: g }
}
fn thermal_budget() -> DepositionBudget {
    DepositionBudget { max_cells: 8,max_links: 8,max_events: 64,max_substeps: 100_000 }
}
fn config(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt,energy_tolerance_j: 1e-8 } }
fn simulation<'a>(schedule: Vec<TimedEvent>,law: &'a dyn HealingLaw, ids: &[u64],max_updates: usize)
    -> BondedDeposition<'a> {
    let bindings = ids.iter().map(|&interface_id| BondBinding { interface_id,
        a: CellId(1),b: CellId(2),law,temperature_weight_a: 0.5 }).collect();
    BondedDeposition::new(schedule,&[(MaterialId(1),&THERMAL)],thermal_budget(),bindings,
        BondBudget { max_interfaces: 8,max_updates }).unwrap()
}
fn schedule(temperature: f64,htc: f64) -> Vec<TimedEvent> {
    vec![birth(1,temperature,htc),birth(2,temperature,htc),
        event(0.0,Event::Contact(contact(7,0.0002,0.2)))]
}
fn close(a: f64,b: f64,tol: f64) { assert!((a-b).abs() < tol,"{a:e} != {b:e}"); }

#[test]
fn actual_contact_window_not_birth_age_drives_exposure() {
    let law = mobility();
    let s = simulation(vec![birth(1,500.0,0.0),birth(2,500.0,0.0),
        event(0.25,Event::Contact(contact(7,0.0002,0.2))),
        event(0.75,Event::Separation(7))],&law,&[7,8],1000);
    let (before,_,_) = s.advance_to(0.2,config(0.2),|| Ok(())).unwrap();
    assert!(before.records().is_empty());
    assert!(before.freeze_bond(7,&map()).is_err());
    let (after,_,_) = before.advance_to(1.0,config(0.2),|| Ok(())).unwrap();
    let record = &after.records()[&7];
    close(record.history().dose(),0.5,1e-14);
    close(record.history().contact_time_s(),0.5,1e-14);
    assert_eq!(record.first_contact_s(),0.25);
    assert_eq!(record.separated_at_s(),Some(0.75));
    assert!(!after.records().contains_key(&8));
    assert!(after.freeze_bond(7,&map()).is_err());
    let (later,_,_) = after.advance_to(5.0,config(0.2),|| Ok(())).unwrap();
    assert_eq!(later.records()[&7].history().dose().to_bits(),record.history().dose().to_bits());
}

#[test]
fn observing_preserves_the_exact_thermal_solution_and_energy_audit() {
    let law = mobility();
    let schedule = vec![birth(1,350.0,100.0),birth(2,550.0,100.0),
        event(0.0,Event::Contact(contact(7,0.0002,0.2)))];
    let original = DepositionSimulation::new(schedule.clone(),&[(MaterialId(1),&THERMAL)],thermal_budget()).unwrap();
    let coupled = simulation(schedule,&law,&[7],10000);
    let (a,aa) = original.advance_to(3.0,config(20.0),|| Ok(())).unwrap();
    let (b,ba,updates) = coupled.advance_to(3.0,config(20.0),|| Ok(())).unwrap();
    assert_eq!(a.cells(),b.thermal().cells());
    assert_eq!(a.links(),b.thermal().links());
    assert_eq!(aa,ba);
    assert_eq!(updates,aa.substeps);
    // Repeat with an explicitly CFL-limited interval (dt_max > 0.8*C/(G+hA)).
    let (_,aa) = original.advance_to(20.0,config(20.0),|| Ok(())).unwrap();
    let (_,ba,updates) = coupled.advance_to(20.0,config(20.0),|| Ok(())).unwrap();
    assert!(aa.substeps > 1);
    assert_eq!(aa,ba);
    assert_eq!(updates,aa.substeps);
}

#[test]
fn thermal_history_refinement_converges_to_independent_cooling_integral() {
    let law = mobility();
    let s = simulation(schedule(450.0,100.0),&law,&[7],10000);
    // Equal cells remain isothermal. C=2 J/K, ambient G=0.1 W/K; tau=20 s.
    // r(T)=max((T-400)/100,0), T=300+150 exp(-t/20).
    let exact = 10.0-20.0*1.5_f64.ln();
    let mut errors = Vec::new();
    for dt in [1.0,0.5,0.25] {
        let (result,audit,_) = s.advance_to(20.0,config(dt),|| Ok(())).unwrap();
        close(audit.residual_j,0.0,1e-8);
        errors.push((result.records()[&7].history().dose()-exact).abs());
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0]);
    for pair in errors.windows(2) { assert!((1.7..2.3).contains(&(pair[0]/pair[1])),"{errors:?}"); }
}

#[test]
fn reheating_restarts_mobility_without_resetting_prior_bond_history() {
    let law = mobility();
    let mut events = schedule(450.0,100.0);
    for id in [CellId(1),CellId(2)] { events.push(event(12.0,Event::Boundary(id,boundary(100.0,100.0)))); }
    let s = simulation(events,&law,&[7],10000);
    let (cold,_,_) = s.advance_to(12.0,config(0.25),|| Ok(())).unwrap();
    assert!(cold.thermal().temperature_k(CellId(1)).unwrap() < 400.0);
    let dose = cold.records()[&7].history().dose();
    let (hot,_,_) = cold.advance_to(14.0,config(0.25),|| Ok(())).unwrap();
    assert!(hot.thermal().temperature_k(CellId(1)).unwrap() > 400.0);
    assert!(hot.records()[&7].history().dose() > dose);
    close(hot.records()[&7].history().contact_time_s(),14.0,1e-12);
}

#[test]
fn equal_time_recontact_cannot_erase_or_reuse_a_bond_cohort() {
    let law = mobility();
    let mut events = schedule(500.0,0.0);
    events.push(event(0.5,Event::Separation(7)));
    events.push(event(0.5,Event::Contact(contact(7,0.0001,0.1))));
    let s = simulation(events,&law,&[7],1000);
    let (accepted,_,_) = s.advance_to(0.25,config(0.25),|| Ok(())).unwrap();
    assert!(matches!(accepted.advance_to(0.5,config(0.25),|| Ok(())),Err(BondingError::Invalid(_))));
    assert_eq!(accepted.thermal().time_s(),0.25);
    assert_eq!(accepted.records()[&7].history().dose(),0.25);
    assert_eq!(accepted.records()[&7].separated_at_s(),None);
    assert_eq!(accepted.records()[&7].contact().area_m2,0.0002);
}

#[test]
fn bond_budget_domain_and_cancellation_refuse_the_whole_candidate() {
    let law = mobility();
    let s = simulation(schedule(500.0,0.0),&law,&[7],1);
    assert!(matches!(s.advance_to(1.0,config(0.1),|| Ok(())),Err(BondingError::Budget)));
    assert!(s.thermal().cells().is_empty() && s.records().is_empty());
    let s = simulation(schedule(500.0,0.0),&law,&[7],10000);
    let mut polls = 0;
    assert!(s.advance_to(1.0,config(0.1),|| {
        polls += 1;
        if polls > 25 { Err(DepositionError::Cancelled) } else { Ok(()) }
    }).is_err());
    assert!(s.thermal().cells().is_empty() && s.records().is_empty());
    let short = MobilityTable::new("short-domain".into(),vec![
        MobilityKnot { temperature_k: 300.0,rate_per_s: 0.0 },
        MobilityKnot { temperature_k: 450.0,rate_per_s: 1.0 },
    ]).unwrap();
    let s = simulation(schedule(500.0,0.0),&short,&[7],10000);
    assert!(matches!(s.advance_to(1.0,config(0.1),|| Ok(())),Err(BondingError::Healing { .. })));
    assert!(s.thermal().cells().is_empty() && s.records().is_empty());
}

#[test]
fn area_partition_preserves_dose_and_total_declared_failure_force() {
    let law = mobility();
    let single = simulation(schedule(500.0,0.0),&law,&[7],10000);
    let split = simulation(vec![birth(1,500.0,0.0),birth(2,500.0,0.0),
        event(0.0,Event::Contact(contact(7,0.0001,0.1))),
        event(0.0,Event::Contact(contact(8,0.0001,0.1)))],&law,&[7,8],10000);
    let (single,_,_) = single.advance_to(1.0,config(0.25),|| Ok(())).unwrap();
    let (split,_,_) = split.advance_to(1.0,config(0.25),|| Ok(())).unwrap();
    let a = single.freeze_bond(7,&map()).unwrap();
    let force = a.area_m2*a.bond.law.unwrap().strength_pa();
    let sum: f64 = [7,8].into_iter().map(|id| {
        let b = split.freeze_bond(id,&map()).unwrap();
        assert_eq!(b.bond.dose,1.0);
        b.area_m2*b.bond.law.unwrap().strength_pa()
    }).sum();
    close(force,sum,1e-12);
}

#[test]
fn same_thermal_grid_restart_preserves_every_bond_and_temperature_bit() {
    let law = mobility();
    let s = simulation(schedule(450.0,100.0),&law,&[7],10000);
    let (whole,_,_) = s.advance_to(2.0,config(0.25),|| Ok(())).unwrap();
    let (prefix,_,_) = s.advance_to(1.0,config(0.25),|| Ok(())).unwrap();
    let (resumed,_,_) = prefix.advance_to(2.0,config(0.25),|| Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(),resumed.thermal().cells());
    assert_eq!(whole.records()[&7].history().dose().to_bits(),resumed.records()[&7].history().dose().to_bits());
    assert_eq!(whole.records()[&7].history().contact_time_s().to_bits(),resumed.records()[&7].history().contact_time_s().to_bits());
}
