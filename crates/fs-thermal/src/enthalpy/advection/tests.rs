//! Real timeline regressions; fixture caloric data are not measured materials.
use std::collections::BTreeMap;
use crate::enthalpy::*;
use super::{Advected, MassFlow, MassBoundary};
use crate::observers::CoupledObservers;

struct Caloric { cp: f64, latent: f64 }
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1e6).contains(&h) { return Err(DepositionError::MaterialDomain); }
        let sensible = if h < 100_000.0 { h }
            else if h <= 100_000.0 + self.latent { 100_000.0 }
            else { h - self.latent };
        Ok(300.0 + sensible/self.cp)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0/self.cp }
}
static LAW: Caloric = Caloric { cp: 2000.0, latent: 100_000.0 };
fn cfg(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-7 } }
fn budget() -> DepositionBudget {
    DepositionBudget { max_cells: 10, max_links: 20, max_events: 30, max_substeps: 100_000 }
}
fn cell(id: u64, mass: f64, h: f64) -> DepositedCell {
    DepositedCell { id: CellId(id), segment_id: id+10, body_id: 99, material_id: MaterialId(1),
        mass_kg: mass, specific_h_j_kg: h,
        geometry: CellGeometry { start_m: [id as f64,0.0,0.0], end_m: [id as f64+1.0,0.0,0.0],
            volume_m3: mass/1000.0, surface_area_m2: 1.0 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 } }
}
fn prepared<'a>(law: &'a dyn ThermalLaw, cells: Vec<DepositedCell>, mut extra: Vec<TimedEvent>)
    -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells.into_iter().map(|c| TimedEvent { time_s: 0.0, event: Event::Birth(c) }).collect();
    events.append(&mut extra);
    let s = DepositionSimulation::new(events, &[(MaterialId(1),law),(MaterialId(2),law)], budget()).unwrap();
    s.advance_to(0.0,cfg(1.0),||Ok(())).unwrap().0
}
fn ring(rate: f64) -> Vec<MassFlow> {
    (0..3).map(|id| MassFlow { id, a: CellId(id), b: CellId((id+1)%3), kg_s: rate }).collect()
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} != {b}"); }

#[test]
fn unequal_mass_upwind_step_conserves_heat_without_external_sources() {
    let masses = [1.0,2.0,4.0];
    let h = [100.0,300.0,500.0];
    let s = prepared(&LAW, (0..3).map(|i|cell(i as u64,masses[i],h[i])).collect(), vec![]);
    let observer = Advected::new((),ring(0.25),3).unwrap();
    let (n,o,a) = s.advance_observed(&observer,0.125,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(a.substeps,1);
    assert_eq!(n.time_s(),0.125);
    for i in 0..3 {
        let c = &n.cells()[&CellId(i as u64)];
        close(c.specific_h_j_kg,h[i]+0.125*0.25*(h[(i+2)%3]-h[i])/masses[i],1e-12);
        assert_eq!(c.mass_kg,masses[i]);
        assert_eq!(c.material_id,MaterialId(1));
        assert_eq!(c.segment_id,i as u64+10);
    }
    close(n.energy_j(),s.energy_j(),1e-10);
    assert_eq!(a.external_heat_j,0.0);
    assert_eq!(a.incoming_mass_kg,0.0);
    close(o.audit().transported_mass_kg,3.0*0.125*0.25,0.0);
    close(o.audit().maximum_courant,0.125*0.25,0.0);
    assert_eq!(observer.audit().substeps,0);
}

#[test]
fn latent_enthalpy_moves_even_when_every_temperature_is_identical() {
    let s = prepared(&LAW, vec![cell(0,1.0,110_000.0),cell(1,1.0,190_000.0),cell(2,1.0,150_000.0)], vec![]);
    let o = Advected::new((),ring(0.5),3).unwrap();
    let (n,_,a) = s.advance_observed(&o,0.25,cfg(1.0),||Ok(())).unwrap();
    for id in 0..3 {
        assert_eq!(s.temperature_k(CellId(id)).unwrap(),350.0);
        assert_eq!(n.temperature_k(CellId(id)).unwrap(),350.0);
    }
    close(n.cells()[&CellId(0)].specific_h_j_kg,115_000.0,1e-10);
    close(n.cells()[&CellId(1)].specific_h_j_kg,180_000.0,1e-10);
    assert_eq!(a.external_heat_j,0.0);
    close(a.residual_j,0.0,1e-9);
}

#[test]
fn changing_only_caloric_data_changes_observed_temperature_not_transport_law() {
    let ordinary = Caloric { cp: 2000.0, latent: 0.0 };
    let different = Caloric { cp: 1000.0, latent: 100_000.0 };
    let cells = vec![cell(0,1.0,110_000.0),cell(1,1.0,190_000.0),cell(2,1.0,150_000.0)];
    let a = prepared(&ordinary,cells.clone(),vec![]);
    let b = prepared(&different,cells,vec![]);
    let o = Advected::new((),ring(0.5),3).unwrap();
    let (a,_,_) = a.advance_observed(&o,0.25,cfg(1.0),||Ok(())).unwrap();
    let (b,_,_) = b.advance_observed(&o,0.25,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(a.cells(),b.cells());
    close(a.temperature_k(CellId(0)).unwrap(),357.5,1e-12);
    close(b.temperature_k(CellId(0)).unwrap(),400.0,1e-12);
}

#[derive(Clone, Default)]
struct Witness { durations: Vec<f64>, watts: f64 }
impl DepositionObserver for Witness {
    type Error = DepositionError;
    fn sources(&mut self, _: &DepositionSimulation<'_>, _: &mut dyn FnMut()->Result<(),DepositionError>)
        -> Result<BTreeMap<CellId,f64>,Self::Error> {
        Ok(BTreeMap::from([(CellId(0),self.watts)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>)->Result<(),Self::Error>{Ok(())}
    fn step(&mut self, _: f64, dt: f64, _: &BTreeMap<CellId,f64>, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error>{self.durations.push(dt);Ok(())}
}

#[test]
fn conduction_and_advection_share_a_combined_monotone_step_and_hit_events() {
    let cells = vec![cell(0,1.0,0.0),cell(1,1.0,80_000.0),cell(2,1.0,40_000.0)];
    let extra = vec![TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 9, a: CellId(0), b: CellId(1), area_m2: 0.1, conductance_w_k: 2000.0 }) },
        TimedEvent { time_s: 0.51, event: Event::Boundary(CellId(0),cells[0].boundary) }];
    let s = prepared(&LAW,cells,extra);
    let o = Advected::new(Witness::default(),ring(1.0),3).unwrap();
    let (n,o,a) = s.advance_observed(&o,0.91,cfg(100.0),||Ok(())).unwrap();
    // Conductive max decay = 1/s, advective = 1/s => combined dt <= 0.4s.
    close(o.inner().durations[0],0.4,1e-15);
    close(o.inner().durations[1],0.11,1e-15);
    assert_eq!(a.substeps,3);
    assert_eq!(a.events,1);
    assert_eq!(n.time_s(),0.91);
    for c in n.cells().values() { assert!((0.0..=80_000.0).contains(&c.specific_h_j_kg)); }
    assert_eq!(a.external_heat_j,0.0);
    close(a.residual_j,0.0,1e-8);
}

#[test]
fn zero_flow_preserves_legacy_steps_and_composed_sources_are_counted_once() {
    let s = prepared(&LAW,vec![cell(0,1.0,40_000.0),cell(1,1.0,50_000.0),cell(2,1.0,60_000.0)],vec![]);
    let source = Witness { watts: 10.0, ..Witness::default() };
    let (legacy,lo,la) = s.advance_observed(&source,0.5,cfg(0.125),||Ok(())).unwrap();
    let empty = Advected::new(source.clone(),ring(0.0),3).unwrap();
    let (wrapped,wo,wa) = s.advance_observed(&empty,0.5,cfg(0.125),||Ok(())).unwrap();
    assert_eq!(legacy.cells(),wrapped.cells());
    assert_eq!(lo.durations,wo.inner().durations);
    assert_eq!(la,wa);
    let pair = CoupledObservers::new(Advected::new((),ring(0.25),3).unwrap(),source);
    let (n,o,a) = s.advance_observed(&pair,0.5,cfg(0.125),||Ok(())).unwrap();
    close(a.external_heat_j,5.0,1e-12);
    close(n.energy_j()-s.energy_j(),5.0,1e-7);
    assert_eq!(o.first.audit().substeps,o.second.durations.len());
    assert_ne!(n.cells()[&CellId(1)].specific_h_j_kg,legacy.cells()[&CellId(1)].specific_h_j_kg);
}

#[test]
fn ordering_and_signed_orientation_do_not_change_solution_and_restart_replays() {
    let s = prepared(&LAW,vec![cell(0,1.0,10_000.0),cell(1,2.0,20_000.0),cell(2,4.0,30_000.0)],vec![]);
    let flows = ring(0.25);
    let equivalent: Vec<_> = flows.iter().rev().map(|f| MassFlow { a:f.b,b:f.a,kg_s:-f.kg_s,..*f }).collect();
    let a = Advected::new((),flows,3).unwrap();
    let b = Advected::new((),equivalent,3).unwrap();
    // dt=0.8/(0.8/0.5 + 0.25)=16/37; align the restart with that actual step.
    let dt = 0.5/(1.0+0.5/(0.8/0.25));
    let (whole,_,_) = s.advance_observed(&a,2.0*dt,cfg(0.5),||Ok(())).unwrap();
    let (same,_,_) = s.advance_observed(&b,2.0*dt,cfg(0.5),||Ok(())).unwrap();
    assert_eq!(whole.cells(),same.cells());
    let (mid,observer,_) = s.advance_observed(&a,dt,cfg(0.5),||Ok(())).unwrap();
    let (resumed,_,_) = mid.advance_observed(&observer,2.0*dt,cfg(0.5),||Ok(())).unwrap();
    assert_eq!(whole.cells(),resumed.cells());
    let reverse = Advected::new((),ring(-0.25),3).unwrap();
    let (reversed,_,_) = s.advance_observed(&reverse,2.0*dt,cfg(0.5),||Ok(())).unwrap();
    assert_ne!(reversed.cells(),whole.cells());
}

#[test]
fn flow_and_material_refusals_leave_original_state_and_observer_unchanged() {
    let cells = vec![cell(0,1.0,10_000.0),cell(1,1.0,20_000.0),cell(2,1.0,30_000.0)];
    let s = prepared(&LAW,cells.clone(),vec![]);
    let mut unbalanced = ring(1.0); unbalanced[1].kg_s = 2.0;
    let mut inactive = ring(1.0); inactive[1].b = CellId(9);
    for flows in [unbalanced,inactive,ring(f64::MAX)] {
        let o = Advected::new(Witness::default(),flows,3).unwrap();
        assert!(s.advance_observed(&o,1.0,cfg(1.0),||Ok(())).is_err());
        assert_eq!(s.time_s(),0.0);
        assert!(o.inner().durations.is_empty());
        assert_eq!(s.cells()[&CellId(0)],cells[0]);
    }
    let mut mixed = cells; mixed[1].material_id = MaterialId(2);
    let mixed = prepared(&LAW,mixed,vec![]);
    let o = Advected::new((),ring(1.0),3).unwrap();
    assert!(mixed.advance_observed(&o,1.0,cfg(1.0),||Ok(())).is_err());
    let mut duplicate = ring(1.0); duplicate[1].id = 0;
    assert!(Advected::new((),duplicate,3).is_err());
    assert!(Advected::new((),ring(1.0),2).is_err());
    assert!(Advected::new((),ring(f64::NAN),3).is_err());
}

#[test]
fn cancellation_and_whole_advance_budget_reject_tentative_transport() {
    let s = prepared(&LAW,vec![cell(0,1.0,10_000.0),cell(1,1.0,20_000.0),cell(2,1.0,30_000.0)],vec![]);
    let o = Advected::new(Witness::default(),ring(1.0),3).unwrap();
    let mut polls = 0;
    let cancelled = s.advance_observed(&o,1.0,cfg(0.1),|| {
        polls += 1;
        if polls > 60 { Err(DepositionError::Cancelled) } else { Ok(()) }
    });
    assert!(matches!(cancelled,Err(DepositionError::Cancelled)));
    assert!(matches!(s.advance_observed_limited(&o,1.0,cfg(0.1),1,||Ok(())),
        Err(DepositionError::Budget("thermal substeps"))));
    assert_eq!(s.time_s(),0.0);
    assert!(o.inner().durations.is_empty());
    assert_eq!(o.audit().substeps,0);
}


fn boundaries(cell: u64, rate: f64, h: f64) -> Vec<MassBoundary> {
    vec![MassBoundary::Inlet { id: 10, cell: CellId(cell), kg_s: rate,
        material_id: MaterialId(1), specific_h_j_kg: h },
        MassBoundary::Outlet { id: 11, cell: CellId(cell), kg_s: rate }]
}

#[test]
fn inlet_and_outlet_carry_enthalpy_once_without_changing_stored_mass() {
    let s = prepared(&LAW,vec![cell(0,2.0,10_000.0)],vec![]);
    let flow = Advected::new((),vec![],0).unwrap().with_boundaries(boundaries(0,1.0,80_000.0),2).unwrap();
    let (n,o,a) = s.advance_observed(&flow,0.25,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(a.substeps,1);
    close(n.cells()[&CellId(0)].specific_h_j_kg,18_750.0,0.0);
    close(n.mass_kg(),2.0,0.0);
    close(a.external_heat_j,17_500.0,0.0);
    close(a.stored_energy_change_j,17_500.0,0.0);
    close(o.audit().external_enthalpy_j,a.external_heat_j,0.0);
    close(o.audit().inlet_mass_kg,0.25,0.0);
    close(o.audit().outlet_mass_kg,0.25,0.0);
    close(o.audit().transported_mass_kg,0.5,0.0);
    close(o.audit().absolute_enthalpy_transfer_j,22_500.0,0.0);
    assert_eq!(a.incoming_mass_kg,0.0);
    assert_eq!(a.incoming_energy_j,0.0);
    assert_eq!(n.cells()[&CellId(0)].boundary,s.cells()[&CellId(0)].boundary);
    assert_eq!(flow.audit(),super::AdvectionAudit::default());
}

#[test]
fn through_flow_separates_interior_exchange_from_reservoir_energy() {
    let s = prepared(&LAW,vec![cell(0,1.0,20_000.0),cell(1,2.0,60_000.0)],vec![]);
    let flow = Advected::new(Witness { watts: 16.0, ..Witness::default() },
        vec![MassFlow {id:0,a:CellId(0),b:CellId(1),kg_s:0.5}],1).unwrap()
        .with_boundaries(vec![
            MassBoundary::Inlet {id:1,cell:CellId(0),kg_s:0.5,material_id:MaterialId(1),specific_h_j_kg:100_000.0},
            MassBoundary::Outlet {id:2,cell:CellId(1),kg_s:0.5},
        ],3).unwrap();
    let (n,o,a) = s.advance_observed(&flow,0.125,cfg(1.0),||Ok(())).unwrap();
    close(n.cells()[&CellId(0)].specific_h_j_kg,25_002.0,0.0);
    close(n.cells()[&CellId(1)].specific_h_j_kg,58_750.0,0.0);
    close(o.audit().external_enthalpy_j,2500.0,0.0);
    close(a.external_heat_j,2502.0,0.0);
    close(a.stored_energy_change_j,2502.0,0.0);
    assert_eq!(a.substeps,o.inner().durations.len());
    assert_eq!(a.substeps,o.audit().substeps);
    assert!(n.cells().values().all(|c|c.boundary.source_w==0.0));
    close(n.mass_kg(),s.mass_kg(),0.0);
}

#[test]
fn reservoir_cooling_refines_to_well_mixed_reference_and_preserves_caloric_datum() {
    let initial_h = 80_000.0;
    let inlet_h = 10_000.0;
    let duration = 2.0;
    let exact = inlet_h+(initial_h-inlet_h)*(-1.0_f64).exp(); // q/m=0.5/s.
    let s = prepared(&LAW,vec![cell(0,2.0,initial_h)],vec![]);
    let flow = Advected::new((),vec![],0).unwrap().with_boundaries(boundaries(0,1.0,inlet_h),2).unwrap();
    let mut errors = Vec::new();
    for dt in [0.125,0.0625,0.03125] {
        let (n,o,a) = s.advance_observed(&flow,duration,cfg(dt),||Ok(())).unwrap();
        assert_eq!(n.time_s(),duration);
        let h = n.cells()[&CellId(0)].specific_h_j_kg;
        assert!((inlet_h..initial_h).contains(&h));
        errors.push((h-exact).abs());
        assert!(a.external_heat_j < 0.0);
        close(a.external_heat_j,o.audit().external_enthalpy_j,1e-8);
        close(a.stored_energy_change_j,a.external_heat_j,1e-7);
        close(o.audit().inlet_mass_kg,2.0,1e-12);
        close(o.audit().outlet_mass_kg,2.0,1e-12);
        assert!(o.audit().maximum_courant<=0.8);
    }
    for pair in errors.windows(2) { assert!((1.8..2.2).contains(&(pair[0]/pair[1])),"{errors:?}"); }
    // Shift the law's energy datum as well as ALL state/reservoir enthalpies.
    struct Shifted;
    impl ThermalLaw for Shifted {
        fn temperature_k(&self,h:f64)->Result<f64,DepositionError>{LAW.temperature_k(h+100_000.0)}
        fn maximum_temperature_slope(&self)->f64{LAW.maximum_temperature_slope()}
    }
    let shifted = prepared(&Shifted,vec![cell(0,2.0,initial_h-100_000.0)],vec![]);
    let sf = Advected::new((),vec![],0).unwrap()
        .with_boundaries(boundaries(0,1.0,inlet_h-100_000.0),2).unwrap();
    let (normal,_,a) = s.advance_observed(&flow,0.5,cfg(0.125),||Ok(())).unwrap();
    let (shifted,_,b) = shifted.advance_observed(&sf,0.5,cfg(0.125),||Ok(())).unwrap();
    close(normal.temperature_k(CellId(0)).unwrap(),shifted.temperature_k(CellId(0)).unwrap(),1e-10);
    close(a.external_heat_j,b.external_heat_j,1e-8);
}

#[test]
fn boundary_identity_support_balance_and_limits_refuse_atomically() {
    let s = prepared(&LAW,vec![cell(0,1.0,50_000.0)],vec![]);
    let mut wrong_material = boundaries(0,1.0,0.0);
    wrong_material[0] = MassBoundary::Inlet {id:10,cell:CellId(0),kg_s:1.0,
        material_id:MaterialId(2),specific_h_j_kg:0.0};
    let mut unbalanced = boundaries(0,1.0,0.0);
    unbalanced[1] = MassBoundary::Outlet {id:11,cell:CellId(0),kg_s:2.0};
    for bad in [wrong_material,unbalanced,boundaries(9,1.0,0.0),boundaries(0,1.0,2e6)] {
        let flow = Advected::new(Witness::default(),vec![],0).unwrap().with_boundaries(bad,2).unwrap();
        assert!(s.advance_observed(&flow,1.0,cfg(1.0),||Ok(())).is_err());
        assert_eq!(flow.audit().substeps,0);
        assert!(flow.inner().durations.is_empty());
    }
    let empty = || Advected::new((),vec![],0).unwrap();
    assert!(empty().with_boundaries(boundaries(0,0.0,0.0),2).is_err());
    assert!(empty().with_boundaries(boundaries(0,1.0,f64::NAN),2).is_err());
    assert!(empty().with_boundaries(boundaries(0,1.0,0.0),1).is_err());
    let mut duplicate = boundaries(0,1.0,0.0); duplicate[1]=duplicate[0];
    assert!(empty().with_boundaries(duplicate,2).is_err());
    let interior = Advected::new((),vec![MassFlow{id:10,a:CellId(0),b:CellId(1),kg_s:0.0}],1).unwrap();
    assert!(interior.with_boundaries(boundaries(0,1.0,0.0),3).is_err());
    assert_eq!(s.time_s(),0.0);
    close(s.cells()[&CellId(0)].specific_h_j_kg,50_000.0,0.0);
}

struct DatumShift { offset: f64 }
impl ThermalLaw for DatumShift {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        LAW.temperature_k(h-self.offset)
    }
    fn maximum_temperature_slope(&self) -> f64 { LAW.maximum_temperature_slope() }
}

#[test]
fn large_formation_datum_does_not_corrupt_a_flowing_calorimeter() {
    // G1/G3: one explicit, well-mixed calorimeter step. The independent
    // balance is dh=dt/m * (mass_rate*(h_in-h)+P), regardless of energy datum.
    let mut reference = None;
    for offset in [0.0, 2.0_f64.powi(52), -2.0_f64.powi(52)] {
        let law = DatumShift { offset };
        let s = prepared(&law,vec![cell(0,1.0,offset+100.0)],vec![]);
        let flow = Advected::new(Witness { watts: 4.0/3.0, ..Witness::default() },vec![],0)
            .unwrap().with_boundaries(boundaries(0,1.0/3.0,offset+200.0),2).unwrap();
        let balance = s.power_balance_observed(&flow,||Ok(())).unwrap();
        close(balance.net_power_w[&CellId(0)],(200.0-100.0)/3.0+4.0/3.0,1e-13);
        let (n,o,a) = s.advance_observed(&flow,0.75,cfg(2.0),||Ok(())).unwrap();
        assert_eq!(a.substeps,1);
        close(n.cells()[&CellId(0)].specific_h_j_kg-offset,126.0,1e-12);
        close(a.external_heat_j,26.0,1e-12);
        close(o.audit().external_enthalpy_j,25.0,1e-12);
        close(a.residual_j,0.0,1e-12);
        assert_eq!(o.audit().absolute_mass_balance_roundoff_kg,0.0);
        assert_eq!(o.audit().absolute_datum_energy_roundoff_j,0.0);
        let temperature = n.temperature_k(CellId(0)).unwrap();
        if let Some(expected) = reference { assert_eq!(temperature,expected); }
        else { reference = Some(temperature); }
        // The gross transported-energy diagnostic deliberately retains datum.
        assert!(o.audit().absolute_enthalpy_transfer_j>0.0);
    }
}

#[test]
fn uniform_ring_retains_small_real_power_beside_large_reference_enthalpy() {
    let law = DatumShift { offset: 2.0_f64.powi(52) };
    let s = prepared(&law,(0..3).map(|id|cell(id,1.0,law.offset+100.0)).collect(),vec![]);
    let flow = Advected::new(Witness { watts: 0.125, ..Witness::default() },ring(0.3),3).unwrap();
    let balance = s.power_balance_observed(&flow,||Ok(())).unwrap();
    assert_eq!(balance.external_power_w,0.125);
    assert_eq!(balance.net_power_w[&CellId(0)],0.125);
    assert_eq!(balance.net_power_w[&CellId(1)],0.0);
    assert_eq!(balance.net_power_w[&CellId(2)],0.0);
    let equilibrium = Advected::new((),ring(0.3),3).unwrap();
    let (next,_,audit) = s.advance_observed(&equilibrium,1.0,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(next.cells(),s.cells());
    assert_eq!(audit.residual_j,0.0);
}

#[test]
fn disconnected_latent_loops_retain_local_references_and_canonical_order() {
    let law = DatumShift { offset: 2.0_f64.powi(52) };
    let h = [110_000.0,190_000.0,150_000.0,150_000.0];
    let s = prepared(&law,(0..4).map(|i|cell(i,1.0,law.offset+h[i as usize])).collect(),vec![]);
    let flows = vec![
        MassFlow {id:0,a:CellId(0),b:CellId(1),kg_s:0.25},
        MassFlow {id:1,a:CellId(1),b:CellId(0),kg_s:0.25},
        MassFlow {id:2,a:CellId(2),b:CellId(3),kg_s:0.25},
        MassFlow {id:3,a:CellId(3),b:CellId(2),kg_s:0.25},
    ];
    let reversed = flows.iter().rev().map(|flow| MassFlow {
        a:flow.b,b:flow.a,kg_s:-flow.kg_s,..*flow
    }).collect();
    let a = Advected::new((),flows,4).unwrap();
    let b = Advected::new((),reversed,4).unwrap();
    let (n,_,audit) = s.advance_observed(&a,0.25,cfg(1.0),||Ok(())).unwrap();
    let (same,_,_) = s.advance_observed(&b,0.25,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(n.cells(),same.cells());
    for (id,expected) in [115_000.0,185_000.0,150_000.0,150_000.0].into_iter().enumerate() {
        assert_eq!(n.cells()[&CellId(id as u64)].specific_h_j_kg-law.offset,expected);
        assert_eq!(n.temperature_k(CellId(id as u64)).unwrap(),350.0);
    }
    assert_eq!(audit.external_heat_j,0.0);
    assert_eq!(audit.residual_j,0.0);
    assert_eq!(a.reference_cells[&CellId(0)],CellId(0));
    assert_eq!(a.reference_cells[&CellId(3)],CellId(2));
}

#[test]
fn admitted_mass_roundoff_has_an_explicit_datum_account() {
    let law = DatumShift { offset: 2.0_f64.powi(52) };
    let h = law.offset+100.0;
    let s = prepared(&law,vec![cell(0,1.0,h)],vec![]);
    let inlet = 0.5_f64;
    let outlet = inlet.next_up();
    let flow = Advected::new((),vec![],0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet {id:1,cell:CellId(0),kg_s:inlet,material_id:MaterialId(1),specific_h_j_kg:h},
        MassBoundary::Outlet {id:2,cell:CellId(0),kg_s:outlet},
    ],2).unwrap();
    let (n,o,a) = s.advance_observed(&flow,0.25,cfg(1.0),||Ok(())).unwrap();
    assert_eq!(n.cells(),s.cells());
    assert_eq!(a.external_heat_j,0.0);
    assert_eq!(o.audit().absolute_mass_balance_roundoff_kg,0.25*(outlet-inlet));
    assert_eq!(o.audit().absolute_datum_energy_roundoff_j,0.25*h*(outlet-inlet));
}
