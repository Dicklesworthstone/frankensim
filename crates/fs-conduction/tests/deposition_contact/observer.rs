//! The material contact wrapper preserves native physical hooks and source pins.
use super::*;
use fs_matdb::{ClaimSet, SelectionPolicy};
use fs_qty::{Dims, QuantitySpec};

#[derive(Clone, Debug, PartialEq)]
enum BulkError {
    Thermal(DepositionError),
    Stage(&'static str),
}
impl From<DepositionError> for BulkError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Bulk {
    calls: [usize; 9],
    extra_heat: bool,
    refuse_step: bool,
    transition_s: Option<f64>,
    refuse_after_transition: bool,
    intervals: Vec<(f64, f64, f64)>,
}
impl Bulk {
    fn new() -> Self {
        Self {
            calls: [0; 9],
            extra_heat: false,
            refuse_step: false,
            transition_s: None,
            refuse_after_transition: false,
            intervals: Vec::new(),
        }
    }
}
impl DepositionObserver for Bulk {
    type Error = BulkError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.calls[0] += 1;
        Ok(())
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        checkpoint()?;
        if self.refuse_after_transition
            && self.transition_s.is_some_and(|at| thermal.time_s() >= at)
        {
            return Err(BulkError::Stage("original source transition"));
        }
        Ok(self.transition_s.filter(|&at| thermal.time_s() < at))
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.calls[1] += 1;
        Ok(if let Some(at) = self.transition_s {
            BTreeMap::from([(CellId(30), if thermal.time_s() < at { 2.0 } else { 6.0 })])
        } else if self.extra_heat {
            BTreeMap::from([(CellId(30), 7.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.calls[2] += 1;
        Ok(BTreeMap::from([(7, 2.0)]))
    }
    fn reservoir_conductances(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.calls[3] += 1;
        Ok(if self.extra_heat {
            BTreeMap::from([(CellId(10), 3.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn heat_transfers(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.calls[4] += 1;
        Ok(if self.extra_heat {
            vec![HeatTransfer {
                a: CellId(30),
                b: CellId(10),
                power_w: 11.0,
            }]
        } else {
            Vec::new()
        })
    }
    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.calls[5] += 1;
        Ok(proposed * 0.5)
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        assert!(dt > 0.0 && end > start);
        self.calls[6] += 1;
        Ok(if self.extra_heat {
            BTreeMap::from([(CellId(10), 13.0)])
        } else {
            BTreeMap::new()
        })
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.calls[7] += 1;
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.calls[8] += 1;
        self.intervals.push((start_s, duration_s, thermal.time_s()));
        if self.refuse_step {
            Err(BulkError::Stage("original bulk step"))
        } else {
            Ok(())
        }
    }
}

#[test]
fn material_contact_preserves_source_transitions_and_their_heat() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let old = initial(&curves, 1.0);
    let (source, id) = source(&cards, 1.0);
    let context = point();
    let contact = MaterialContactLink {
        link: 7,
        sides: sides(&curves),
        point: binding(&cards, &source, &context),
        selection: ClaimSelection::Pinned(id),
    };
    let mut bulk = Bulk::new();
    bulk.transition_s = Some(0.125);
    let observer = MaterialContactObserver::new(bulk, &old, vec![contact], || Ok(())).unwrap();
    let config = StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-10,
    };
    let (next, accepted, audit) = old
        .advance_observed(&observer, 0.5, config, || Ok(()))
        .unwrap();

    // The contact and the inner half-step limiter must not erase or approach
    // the pulse edge asymptotically. The original owner clips after all CFLs.
    assert_eq!(
        accepted.inner().intervals,
        [(0.0, 0.125, 0.125), (0.125, 0.375, 0.5)]
    );
    assert_eq!(audit.substeps, 2);
    assert_eq!(audit.events, 0);
    let heat_j = 2.0 * 0.125 + 6.0 * 0.375;
    assert!((audit.external_heat_j - heat_j).abs() < 1e-12);
    assert!((next.energy_j() - old.energy_j() - heat_j).abs() < 1e-12);
    assert_eq!(
        accepted.endpoint_sources()[0].interface.receipt().selected,
        id
    );
    assert!(observer.inner().intervals.is_empty());

    let (middle, first, _) = old
        .advance_observed(&observer, 0.125, config, || Ok(()))
        .unwrap();
    let (split, last, _) = middle
        .advance_observed(&first, 0.5, config, || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), split.cells());
    assert_eq!(accepted.inner().intervals, last.inner().intervals);
}

#[test]
fn contact_wrapped_transition_failure_rolls_back_the_entire_heat_history() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let old = initial(&curves, 1.0);
    let (source, id) = source(&cards, 1.0);
    let context = point();
    let contact = MaterialContactLink {
        link: 7,
        sides: sides(&curves),
        point: binding(&cards, &source, &context),
        selection: ClaimSelection::Pinned(id),
    };
    let mut bulk = Bulk::new();
    bulk.transition_s = Some(0.125);
    bulk.refuse_after_transition = true;
    let observer = MaterialContactObserver::new(bulk, &old, vec![contact], || Ok(())).unwrap();
    let snapshot = old.cells().clone();
    let config = StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-10,
    };
    let error = old
        .advance_observed(&observer, 0.5, config, || Ok(()))
        .err()
        .unwrap();
    assert_eq!(
        error,
        MaterialContactObserverError::Bulk(BulkError::Stage("original source transition"))
    );
    assert_eq!(old.cells(), &snapshot);
    assert_eq!(old.time_s(), 0.0);
    assert!(observer.inner().intervals.is_empty());
    assert_eq!(observer.inner().calls, [0; 9]);
}

#[test]
fn all_native_hooks_are_forwarded_and_fresh_bulk_is_composed_once() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let old = initial(&curves, 1.0);
    let (source, id) = source(&cards, 1.0);
    let context = point();
    let contacts = vec![MaterialContactLink {
        link: 7,
        sides: sides(&curves),
        point: binding(&cards, &source, &context),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let mut bulk = Bulk::new();
    bulk.extra_heat = true;
    let mut observer =
        MaterialContactObserver::new(bulk.clone(), &old, contacts, || Ok(())).unwrap();
    assert_eq!(
        observer.inner(),
        &bulk,
        "attachment probes a clone, not retained history"
    );
    assert_eq!(observer.bindings()[0].selection, ClaimSelection::Pinned(id));
    let mut check = || Ok(());
    observer.begin().unwrap();
    assert_eq!(
        observer.sources(&old, &mut check).unwrap()[&CellId(30)],
        7.0
    );
    let first = observer.conductances(&old, &mut check).unwrap();
    let second = observer.conductances(&old, &mut check).unwrap();
    assert_eq!(first, second);
    assert!((first[&7] - 4.0 / 7.0).abs() < 1e-14);
    assert_eq!(observer.last_samples()[0].bulk_conductance_w_k, 2.0);
    assert_eq!(
        observer.last_samples()[0]
            .series
            .as_ref()
            .unwrap()
            .terms()
            .len(),
        2
    );
    assert_eq!(
        observer.reservoir_conductances(&old, &mut check).unwrap()[&CellId(10)],
        3.0
    );
    assert_eq!(
        observer.heat_transfers(&old, &mut check).unwrap()[0].power_w,
        11.0
    );
    assert_eq!(
        observer.limit_timestep_s(&old, 0.4, &mut check).unwrap(),
        0.2
    );
    assert!(observer.requires_step_interval());
    assert_eq!(
        observer
            .step_sources(&old, 0.0, 0.1, 0.1, &mut check)
            .unwrap()[&CellId(10)],
        13.0
    );
    observer
        .event(
            &Event::Boundary(CellId(10), old.cells()[&CellId(10)].boundary),
            &old,
        )
        .unwrap();
    observer
        .step(0.0, 0.1, &BTreeMap::new(), &old, &mut check)
        .unwrap();
    assert_eq!(observer.inner().calls, [1, 1, 2, 1, 1, 1, 1, 1, 1]);
    assert_eq!(
        observer.endpoint_sources()[0].interface.receipt().selected,
        id
    );
}

#[test]
fn real_native_step_retains_actual_endpoint_sources_and_original_bulk_error() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let old = initial(&curves, 1.0);
    let (source, id) = source(&cards, 1.0);
    let context = point();
    let contact = MaterialContactLink {
        link: 7,
        sides: sides(&curves),
        point: binding(&cards, &source, &context),
        selection: ClaimSelection::Pinned(id),
    };
    let observer =
        MaterialContactObserver::new(Bulk::new(), &old, vec![contact], || Ok(())).unwrap();
    let config = StepConfig {
        max_dt_s: 0.2,
        energy_tolerance_j: 1e-10,
    };
    let (next, accepted, audit) = old
        .advance_observed(&observer, 0.1, config, || Ok(()))
        .unwrap();
    let q = 0.1 * (4.0 / 7.0) * 50.0;
    assert!((next.cells()[&CellId(30)].specific_h_j_kg - (25.0 - q)).abs() < 1e-12);
    assert!((next.cells()[&CellId(10)].specific_h_j_kg - (-50.0 + q)).abs() < 1e-12);
    assert!((next.energy_j() - old.energy_j()).abs() < 1e-12);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(
        accepted.last_samples()[0].interface.value_m2_k_per_w(),
        1.25
    );
    assert!(
        (accepted.endpoint_sources()[0].interface.value_m2_k_per_w() - (1.0 + (25.0 - q) / 100.0))
            .abs()
            < 1e-13
    );
    assert_eq!(accepted.inner().calls[8], 1);
    assert_eq!(observer.inner().calls, [0; 9]);
    let mut failing = Bulk::new();
    failing.refuse_step = true;
    let failing = MaterialContactObserver::new(failing, &old, vec![contact], || Ok(())).unwrap();
    let error = old
        .advance_observed(&failing, 0.1, config, || Ok(()))
        .err()
        .unwrap();
    assert_eq!(
        error,
        MaterialContactObserverError::Bulk(BulkError::Stage("original bulk step"))
    );
    assert_eq!(failing.inner().calls, [0; 9]);
    assert_eq!(old.cells()[&CellId(30)].specific_h_j_kg, 25.0);
}

#[test]
fn endpoint_domain_exit_cannot_switch_contact_claims_and_short_retry_is_clean() {
    let cards = [bulk("A"), bulk("B")];
    let curves = [curve(&cards[0]), curve(&cards[1])];
    let old = initial(&curves, 1.0);
    let context = point();
    let mut original = claim("fa", Dims::NONE, 1.0);
    original.validity =
        original
            .validity
            .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.24, 1.0);
    let mut alternate = claim("fa", Dims::NONE, 4.0);
    alternate.validity =
        alternate
            .validity
            .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.0, 0.24);
    let mut claims = ClaimSet::new();
    let original_id = claims.insert_claim(original).unwrap();
    let alternate_id = claims.insert_claim(alternate).unwrap();
    let source = interface(&cards, claims);
    let contact = MaterialContactLink {
        link: 7,
        sides: sides(&curves),
        point: binding(&cards, &source, &context),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    };
    let observer =
        MaterialContactObserver::new(Bulk::new(), &old, vec![contact], || Ok(())).unwrap();
    assert_eq!(
        observer.bindings()[0].selection,
        ClaimSelection::Pinned(original_id)
    );
    assert_ne!(original_id, alternate_id);
    let config = StepConfig {
        max_dt_s: 0.2,
        energy_tolerance_j: 1e-10,
    };
    assert!(matches!(
        old.advance_observed(&observer, 0.1, config, || Ok(())),
        Err(MaterialContactObserverError::Contact(_))
    ));
    assert_eq!(observer.inner().calls, [0; 9]);
    let (next, accepted, _) = old
        .advance_observed(&observer, 0.01, config, || Ok(()))
        .unwrap();
    assert_eq!(
        accepted.endpoint_sources()[0].interface.receipt().selected,
        original_id
    );
    assert!((next.cells()[&CellId(30)].specific_h_j_kg - (25.0 - 2.0 / 7.0)).abs() < 1e-12);
}
