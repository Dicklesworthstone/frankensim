use super::*;
use crate::enthalpy::*;

#[derive(Clone, Debug, PartialEq)]
struct Producer {
    cell: CellId,
    watts: f64,
    link: Option<(u64, f64)>,
    source_states: Vec<(f64, f64)>,
    grid: Vec<(f64, f64)>,
    energy_j: f64,
    events: usize,
    fail_step: Option<usize>,
}
impl Producer {
    fn new(cell: u64, watts: f64) -> Self {
        Self { cell: CellId(cell), watts, link: None, source_states: vec![], grid: vec![],
            energy_j: 0.0, events: 0, fail_step: None }
    }
}
impl DepositionObserver for Producer {
    type Error = DepositionError;
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        // Keep invalid-field tests possible without a panic inside the producer.
        let h = thermal.cells().get(&self.cell).map_or(0.0, |c| c.specific_h_j_kg);
        self.source_states.push((thermal.time_s(), h));
        Ok(BTreeMap::from([(self.cell, self.watts)]))
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<u64, f64>, Self::Error> {
        if let Some((time, h)) = self.source_states.last() {
            assert_eq!(*time, thermal.time_s());
            assert_eq!(*h, thermal.cells()[&self.cell].specific_h_j_kg);
        }
        Ok(self.link.into_iter().collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events += 1;
        Ok(())
    }
    fn step(&mut self, start: f64, dt: f64, _: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        self.grid.push((start, dt));
        self.energy_j += self.watts * dt;
        assert_eq!(thermal.cells()[&self.cell].boundary.source_w, 0.0);
        if self.fail_step == Some(self.grid.len()) {
            return Err(DepositionError::InvalidInput("observer rejected endpoint"));
        }
        Ok(())
    }
}
struct Law(f64);
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=10000.0).contains(&h) { return Err(DepositionError::MaterialDomain); }
        Ok(300.0 + h / self.0)
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 / self.0 }
}
static A: Law = Law(2.0);
static B: Law = Law(4.0);
fn config() -> StepConfig { StepConfig { max_dt_s: 0.25, energy_tolerance_j: 1e-9 } }
fn initial(link: bool) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = (1..=2).map(|id| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id), segment_id: id, body_id: 1, material_id: MaterialId(id),
            mass_kg: 1.0, specific_h_j_kg: 0.0,
            geometry: CellGeometry { start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64 + 1.0, 0.0, 0.0], volume_m3: 1.0, surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    if link {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 7, a: CellId(1), b: CellId(2), area_m2: 1.0, conductance_w_k: 0.0,
        }) });
    }
    DepositionSimulation::new(events, &[(MaterialId(1), &A), (MaterialId(2), &B)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 4, max_substeps: 10000 }).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-9, "{a} != {b}"); }

#[test]
fn independent_sources_add_once_and_see_identical_material_history() {
    let old = initial(false);
    let observers = CoupledObservers::new(Producer::new(1, 4.0), Producer::new(1, 6.0));
    let (next, accepted, audit) = old.advance_observed(&observers, 1.0, config(), || Ok(())).unwrap();
    near(next.energy_j(), 10.0);
    near(next.temperature_k(CellId(1)).unwrap(), 305.0);
    near(audit.external_heat_j, accepted.first.energy_j + accepted.second.energy_j);
    assert_eq!(accepted.first.source_states, accepted.second.source_states);
    assert_eq!(accepted.first.grid, accepted.second.grid);
    assert_eq!(accepted.first.events, 2);
    assert_eq!(accepted.second.events, 2);
    assert_eq!(observers.first.grid.len(), 0);
    assert_eq!(old.time_s(), 0.0);
    assert!(old.cells().is_empty());
}

#[test]
fn material_data_changes_response_without_changing_source_or_solver() {
    let pair = CoupledObservers::new(Producer::new(1, 8.0), Producer::new(2, 8.0));
    let (next, _, audit) = initial(false).advance_observed(&pair, 1.0, config(), || Ok(())).unwrap();
    near(next.cells()[&CellId(1)].specific_h_j_kg, 8.0);
    near(next.cells()[&CellId(2)].specific_h_j_kg, 8.0);
    near(next.temperature_k(CellId(1)).unwrap(), 304.0);
    near(next.temperature_k(CellId(2)).unwrap(), 302.0);
    near(audit.external_heat_j, 16.0);
}

#[test]
fn one_transport_owner_controls_cfl_and_both_source_energy_accounts() {
    let mut first = Producer::new(1, 4.0);
    first.link = Some((7, 100.0));
    let pair = CoupledObservers::new(first, Producer::new(1, 6.0));
    let (next, accepted, audit) = initial(true).advance_observed(&pair, 1.0, config(), || Ok(())).unwrap();
    assert!(audit.substeps > 4);
    assert_eq!(accepted.first.grid, accepted.second.grid);
    assert_eq!(accepted.first.grid.len(), audit.substeps);
    near(next.energy_j(), 10.0);
    near(audit.external_heat_j, accepted.first.energy_j + accepted.second.energy_j);
    assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
    assert_eq!(next.links()[&7].conductance_w_k, 0.0); // authored data restored
}

#[test]
fn duplicate_transport_ownership_refuses_even_when_values_match() {
    for second_g in [100.0, 200.0] {
        let mut first = Producer::new(1, 4.0);
        first.link = Some((7, 100.0));
        let mut second = Producer::new(2, 6.0);
        second.link = Some((7, second_g));
        let old = initial(true);
        let pair = CoupledObservers::new(first, second);
        let result = old.advance_observed(&pair, 1.0, config(), || Ok(()));
        assert!(matches!(result, Err(ObserverCouplingError::ConflictingConductance(7))));
        assert!(old.cells().is_empty());
        assert!(pair.first.source_states.is_empty());
    }
}

#[test]
fn nested_owners_support_signed_cooling_without_a_second_thermal_step() {
    let pair = CoupledObservers::new(
        CoupledObservers::new(Producer::new(1, 9.0), Producer::new(1, -2.0)),
        Producer::new(1, 3.0),
    );
    let (next, accepted, audit) = initial(false).advance_observed(&pair, 1.0, config(), || Ok(())).unwrap();
    near(next.energy_j(), 10.0);
    near(audit.external_heat_j, 10.0);
    assert_eq!(accepted.first.first.grid, accepted.first.second.grid);
    assert_eq!(accepted.first.first.grid, accepted.second.grid);
}

#[test]
fn late_observer_refusal_rolls_back_all_histories_births_and_energy() {
    let old = initial(false);
    let mut second = Producer::new(1, 6.0);
    second.fail_step = Some(3);
    let pair = CoupledObservers::new(Producer::new(1, 4.0), second);
    let result = old.advance_observed(&pair, 1.0, config(), || Ok(()));
    assert!(matches!(result, Err(ObserverCouplingError::Second(DepositionError::InvalidInput(_)))));
    assert!(old.cells().is_empty());
    assert_eq!(old.time_s(), 0.0);
    assert!(pair.first.grid.is_empty());
    assert!(pair.second.grid.is_empty());
    assert_eq!(pair.first.energy_j, 0.0);
}

#[test]
fn invalid_sources_cannot_cancel_or_overwrite_each_other_into_validity() {
    let pairs = [
        (Producer::new(1, f64::INFINITY), Producer::new(1, f64::NEG_INFINITY)),
        (Producer::new(1, f64::MAX), Producer::new(1, f64::MAX)),
        (Producer::new(99, 1.0), Producer::new(99, -1.0)),
        (Producer::new(1, 1.0), Producer::new(2, f64::NAN)),
    ];
    for (a, b) in pairs {
        assert!(initial(false).advance_observed(&CoupledObservers::new(a, b), 1.0, config(), || Ok(())).is_err());
    }
}

#[test]
fn invalid_link_override_refuses_before_either_history_commits() {
    for (id, g) in [(8, 1.0), (7, -1.0), (7, f64::NAN), (7, f64::INFINITY)] {
        let mut source = Producer::new(1, 4.0);
        source.link = Some((id, g));
        assert!(initial(true).advance_observed(&CoupledObservers::new((), source), 1.0, config(), || Ok(())).is_err());
    }
}

#[test]
fn no_op_composition_preserves_existing_state_and_audit_exactly() {
    let old = initial(true);
    let source = Producer::new(1, 7.0);
    let (reference, observer, audit) = old.advance_observed(&source, 1.0, config(), || Ok(())).unwrap();
    let (left, pair, left_audit) = old.advance_observed(&CoupledObservers::new(source.clone(), ()), 1.0, config(), || Ok(())).unwrap();
    assert_eq!(left.cells(), reference.cells());
    assert_eq!(left_audit, audit);
    assert_eq!(pair.into_parts().0, observer);
    let (right, pair, right_audit) = old.advance_observed(&CoupledObservers::new((), source), 1.0, config(), || Ok(())).unwrap();
    assert_eq!(right.cells(), reference.cells());
    assert_eq!(right_audit, audit);
    assert_eq!(pair.into_parts().1, observer);
}

#[test]
fn cancellation_at_every_observed_checkpoint_keeps_originals() {
    let old = initial(false);
    let pair = CoupledObservers::new(Producer::new(1, 4.0), Producer::new(2, 6.0));
    let mut polls = 0;
    old.advance_observed(&pair, 0.25, config(), || { polls += 1; Ok(()) }).unwrap();
    for stop in 1..=polls {
        let mut seen = 0;
        let result = old.advance_observed(&pair, 0.25, config(), || {
            seen += 1;
            if seen == stop { Err(DepositionError::Cancelled) } else { Ok(()) }
        });
        assert!(result.is_err(), "checkpoint {stop}/{polls}");
        assert!(old.cells().is_empty());
        assert!(pair.first.grid.is_empty());
        assert!(pair.second.grid.is_empty());
    }
}


// Duration-dependent power must use the one interval accepted by the original
// conduction/advection/event owner, including through nested observer wrappers.
#[derive(Clone, Debug, PartialEq)]
struct IntervalProducer {
    cell: CellId,
    watts: f64,
    slope_w_s: f64,
    source_h: Option<f64>,
    prepared: Vec<(f64, f64, f64, f64)>,
    accepted: usize,
    fail_prepare: Option<usize>,
    fail_accept: Option<usize>,
}
impl IntervalProducer {
    fn new(watts: f64, slope_w_s: f64) -> Self {
        Self { cell: CellId(1), watts, slope_w_s, source_h: None,
            prepared: vec![], accepted: 0, fail_prepare: None, fail_accept: None }
    }
    fn heat_j(&self) -> f64 {
        self.prepared.iter().map(|(_, dt, _, _)| dt * (self.watts + self.slope_w_s * dt)).sum()
    }
}
impl DepositionObserver for IntervalProducer {
    type Error = DepositionError;
    fn requires_step_interval(&self) -> bool { true }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        assert!(self.source_h.is_none(), "source epoch evaluated twice");
        self.source_h = Some(thermal.cells()[&CellId(1)].specific_h_j_kg);
        Ok(BTreeMap::new())
    }
    fn step_sources(&mut self, thermal: &DepositionSimulation<'_>,
        start: f64, dt: f64, end: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        let h = thermal.cells()[&CellId(1)].specific_h_j_kg;
        assert_eq!(self.source_h, Some(h), "thermal state changed before interval preparation");
        assert_eq!(thermal.time_s().to_bits(), start.to_bits());
        assert_eq!(self.prepared.len(), self.accepted, "interval prepared twice");
        self.prepared.push((start, dt, end, h));
        if self.fail_prepare == Some(self.prepared.len()) {
            return Err(DepositionError::InvalidInput("interval preparation refused"));
        }
        Ok(BTreeMap::from([(self.cell, self.watts + self.slope_w_s * dt)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        assert!(self.source_h.is_none(), "event crossed an unaccepted interval");
        Ok(())
    }
    fn step(&mut self, start: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        let (prepared_start, prepared_dt, end, h) = self.prepared[self.accepted];
        assert_eq!(prepared_start.to_bits(), start.to_bits());
        assert_eq!(prepared_dt.to_bits(), dt.to_bits());
        assert_eq!(end.to_bits(), thermal.time_s().to_bits());
        near(old[&CellId(1)], 300.0 + h / 2.0);
        self.accepted += 1;
        self.source_h = None;
        if self.fail_accept == Some(self.accepted) {
            return Err(DepositionError::InvalidInput("interval acceptance refused"));
        }
        Ok(())
    }
}

fn interval_initial() -> DepositionSimulation<'static> {
    let (active, _) = initial(false).advance_to(0.0, config(), || Ok(())).unwrap();
    let mut schedule: Vec<_> = active.cells().values().cloned().map(|mut cell| {
        cell.material_id = MaterialId(1); // one advected material and caloric datum
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    schedule.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 7, a: CellId(1), b: CellId(2), area_m2: 1.0, conductance_w_k: 100.0,
    }) });
    schedule.push(TimedEvent { time_s: 0.035, event: Event::Boundary(
        CellId(1), active.cells()[&CellId(1)].boundary,
    ) });
    DepositionSimulation::new(schedule, &[(MaterialId(1), &A)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 4, max_substeps: 100 }).unwrap()
}

#[test]
fn interval_heat_uses_combined_cfl_and_exact_event_endpoints_once() {
    use crate::enthalpy::advection::{Advected, MassFlow};
    let physical = CoupledObservers::new(
        CoupledObservers::new(Producer::new(1, 4.0), IntervalProducer::new(3.0, 100.0)),
        IntervalProducer::new(-1.0, 50.0),
    );
    let observer = Advected::new(physical, vec![
        MassFlow { id: 10, a: CellId(1), b: CellId(2), kg_s: 10.0 },
        MassFlow { id: 11, a: CellId(2), b: CellId(1), kg_s: 10.0 },
    ], 2).unwrap();
    let old = interval_initial();
    let (next, accepted, audit) = old.advance_observed(&observer, 0.05, config(), || Ok(())).unwrap();
    let owners = accepted.inner();
    let a = &owners.first.second;
    let b = &owners.second;
    assert_eq!(a.prepared, b.prepared);
    assert_eq!(a.accepted, audit.substeps);
    assert_eq!(b.accepted, audit.substeps);
    assert_eq!(audit.substeps, 5);
    // Conductive 0.8/50 s and advective 0.8/10 s combine before clipping.
    near(a.prepared[0].1, 0.8 / 60.0);
    assert!(a.prepared.iter().any(|&(_, _, end, _)| end == 0.035));
    assert_eq!(a.prepared.last().unwrap().2, 0.05);
    for (&(start, dt, _, h), (&grid, &state)) in a.prepared.iter()
        .zip(owners.first.first.grid.iter().zip(&owners.first.first.source_states)) {
        assert_eq!((start, dt), grid);
        assert_eq!((start, h), state);
    }
    near(audit.external_heat_j, owners.first.first.energy_j + a.heat_j() + b.heat_j());
    near(next.energy_j(), audit.external_heat_j);
    assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 0.0);
    assert_eq!(old.time_s(), 0.0);
    assert!(observer.inner().first.second.prepared.is_empty());
}

#[test]
fn interval_preparation_acceptance_and_cancellation_roll_back_together() {
    let old = interval_initial();
    let observer = CoupledObservers::new(Producer::new(1, 4.0), IntervalProducer::new(3.0, 100.0));
    for prepare in [true, false] {
        let mut failed = observer.clone();
        if prepare { failed.second.fail_prepare = Some(2); }
        else { failed.second.fail_accept = Some(2); }
        assert!(old.advance_observed(&failed, 0.05, config(), || Ok(())).is_err());
        assert!(failed.first.grid.is_empty());
        assert!(failed.second.prepared.is_empty());
        assert!(old.cells().is_empty());
    }
    let mut polls = 0;
    old.advance_observed(&observer, 0.05, config(), || { polls += 1; Ok(()) }).unwrap();
    for stop in 1..=polls {
        let mut seen = 0;
        assert!(old.advance_observed(&observer, 0.05, config(), || {
            seen += 1;
            if seen == stop { Err(DepositionError::Cancelled) } else { Ok(()) }
        }).is_err(), "interval checkpoint {stop}/{polls}");
        assert!(observer.second.prepared.is_empty());
        assert!(old.cells().is_empty());
    }
}

#[test]
fn invalid_interval_heat_cannot_be_hidden_by_composition() {
    let old = initial(false);
    for watts in [f64::NAN, f64::INFINITY] {
        assert!(old.advance_observed(&IntervalProducer::new(watts, 0.0),
            0.25, config(), || Ok(())).is_err());
    }
    let mut missing = IntervalProducer::new(0.0, 0.0);
    missing.cell = CellId(99);
    assert!(old.advance_observed(&missing, 0.25, config(), || Ok(())).is_err());
    let too_large = CoupledObservers::new(
        IntervalProducer::new(f64::MAX, 0.0), IntervalProducer::new(f64::MAX, 0.0),
    );
    assert!(old.advance_observed(&too_large, 0.25, config(), || Ok(())).is_err());
    let (active, _) = old.advance_to(0.0, config(), || Ok(())).unwrap();
    let wrapped = crate::enthalpy::advection::Advected::new(
        CoupledObservers::new((), IntervalProducer::new(1.0, 1.0)), vec![], 0,
    ).unwrap();
    assert!(wrapped.requires_step_interval());
    assert!(active.power_balance_observed(&wrapped, || Ok(())).is_err());
    // Ordinary instantaneous observers still expose their complete power field.
    near(active.power_balance_observed(&Producer::new(1, 7.0), || Ok(())).unwrap().external_power_w, 7.0);
}
