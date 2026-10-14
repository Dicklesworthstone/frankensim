//! Generated sources enter the original mass/enthalpy update, including plateaus.
use std::collections::BTreeMap;
use fs_conduction::transient::deposition::*;

struct Law;
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(if h < 20.0 { 300.0+h/2.0 } else if h <= 120.0 { 310.0 } else { 310.0+(h-120.0)/2.0 })
    }
    fn maximum_temperature_slope(&self) -> f64 { 0.5 }
}
const LAW: Law = Law;
fn boundary(watts: f64) -> Boundary {
    Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: watts }
}
fn config() -> StepConfig { StepConfig { max_dt_s: 0.125, energy_tolerance_j: 1e-10 } }
fn setup(initial_h: f64, authored: f64, mut events: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    let cell = DepositedCell { id: CellId(1), segment_id: 1, body_id: 1,
        material_id: MaterialId(1), mass_kg: 1.0, specific_h_j_kg: initial_h,
        geometry: CellGeometry { start_m: [0.0;3], end_m: [1.0,0.0,0.0], volume_m3: 1.0, surface_area_m2: 6.0 },
        boundary: boundary(authored) };
    events.insert(0, TimedEvent { time_s: 0.0, event: Event::Birth(cell) });
    DepositionSimulation::new(events, &[(MaterialId(1), &LAW as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 8, max_substeps: 1000 })
        .unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}

#[derive(Clone, Default)]
struct Source {
    watts: f64, calls: usize, energy_j: f64, last_w: f64,
    invalid: bool, fail_after: Option<usize>, feedback: bool,
    seen_authored: Vec<f64>,
}
impl DepositionObserver for Source {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(),Self::Error> { Ok(()) }
    fn sources(&mut self, s: &DepositionSimulation<'_>, check: &mut dyn FnMut()->Result<(),DepositionError>)
        -> Result<BTreeMap<CellId,f64>, Self::Error> {
        check()?;
        if self.fail_after == Some(self.calls) { return Err(DepositionError::Cancelled); }
        self.calls += 1;
        self.seen_authored.push(s.cells()[&CellId(1)].boundary.source_w);
        self.last_w = if self.feedback { 2.0 * (s.temperature_k(CellId(1))?-300.0) } else { self.watts };
        Ok(BTreeMap::from([(CellId(if self.invalid { 99 } else { 1 }), self.last_w)]))
    }
    fn step(&mut self, _:f64, dt:f64, _: &BTreeMap<CellId,f64>, s:&DepositionSimulation<'_>,
        _: &mut dyn FnMut()->Result<(),DepositionError>) -> Result<(),Self::Error> {
        self.energy_j += self.last_w*dt;
        // Generated watts must already have been removed from public boundaries.
        assert_eq!(s.cells()[&CellId(1)].boundary.source_w, *self.seen_authored.last().unwrap());
        Ok(())
    }
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<1e-9, "{a} != {b}"); }

#[test]
fn generated_heat_crosses_latent_plateau_without_inventing_heat_capacity() {
    let original=setup(0.0,0.0,vec![]);
    let source=Source { watts:80.0, ..Source::default() };
    let (mid,owner,audit)=original.advance_observed(&source,1.0,config(),||Ok(())).unwrap();
    near(mid.cells()[&CellId(1)].specific_h_j_kg,80.0);
    near(mid.temperature_k(CellId(1)).unwrap(),310.0);
    near(owner.energy_j,80.0); near(audit.external_heat_j,80.0);
    let (end,owner,audit)=mid.advance_observed(&owner,2.0,config(),||Ok(())).unwrap();
    near(end.cells()[&CellId(1)].specific_h_j_kg,160.0);
    near(end.temperature_k(CellId(1)).unwrap(),330.0);
    near(owner.energy_j,160.0); near(audit.external_heat_j,80.0);
    near(original.energy_j(),0.0); assert_eq!(source.calls,0);
    assert_eq!(end.mass_kg(),original.mass_kg());
}

#[test]
fn sources_compose_with_boundary_events_without_accumulation_or_reset() {
    let original=setup(0.0,4.0,vec![TimedEvent {time_s:0.5,event:Event::Boundary(CellId(1),boundary(0.0))}]);
    let source=Source {watts:6.0,..Source::default()};
    let (end,owner,audit)=original.advance_observed(&source,1.0,config(),||Ok(())).unwrap();
    near(end.energy_j(),8.0); near(audit.external_heat_j,8.0); near(owner.energy_j,6.0);
    assert_eq!(owner.seen_authored,vec![4.0,4.0,4.0,4.0,0.0,0.0,0.0,0.0]);
    let (off,_,audit)=end.advance_observed(&Source::default(),2.0,config(),||Ok(())).unwrap();
    near(off.energy_j(),8.0); near(audit.external_heat_j,0.0);
}

#[test]
fn live_source_is_reevaluated_each_substep_and_zero_source_matches_legacy() {
    let original=setup(1.0,0.0,vec![]);
    let source=Source {feedback:true,..Source::default()};
    let (end,owner,_)=original.advance_observed(&source,1.0,config(),||Ok(())).unwrap();
    near(end.energy_j(),1.125_f64.powi(8)); assert_eq!(owner.calls,8);
    let (plain,a)=original.advance_to(1.0,config(),||Ok(())).unwrap();
    let (zero,_,b)=original.advance_observed(&Source::default(),1.0,config(),||Ok(())).unwrap();
    assert_eq!(plain.cells(),zero.cells()); assert_eq!(a,b);
}

#[test]
fn signed_cooling_removes_latent_energy_at_constant_temperature() {
    let original=setup(100.0,0.0,vec![]);
    let (end,owner,audit)=original.advance_observed(&Source {watts:-40.0,..Source::default()},1.0,config(),||Ok(())).unwrap();
    near(end.temperature_k(CellId(1)).unwrap(),310.0); near(end.energy_j(),60.0);
    near(owner.energy_j,-40.0); near(audit.external_heat_j,-40.0);
}

#[test]
fn source_failure_invalid_mapping_budget_and_domain_failure_are_atomic() {
    let original=setup(1.0,0.0,vec![]);
    for source in [Source {invalid:true,..Source::default()},
        Source {watts:f64::NAN,..Source::default()},
        Source {watts:1.0,fail_after:Some(2),..Source::default()},
        Source {watts:1e9,..Source::default()}] {
        assert!(original.advance_observed(&source,1.0,config(),||Ok(())).is_err());
        assert_eq!(source.calls,0); near(original.time_s(),0.0); near(original.energy_j(),1.0);
    }
    assert!(original.advance_observed_limited(&Source {watts:1.0,..Source::default()},1.0,config(),0,||Ok(())).is_err());
    let (_,s,_)=original.advance_observed_limited(&Source::default(),0.0,config(),0,||Ok(())).unwrap();
    assert_eq!(s.calls,0);
}
