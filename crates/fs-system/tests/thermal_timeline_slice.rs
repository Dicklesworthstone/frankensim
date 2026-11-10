//! Owner slice: the fs-thermal enthalpy timeline on the composition kernel.
//!
//! `fs_thermal::component::ThermalTimeline` puts the explicit
//! `DepositionSimulation` behind the Component contract. Here it passes the
//! conformance check, is compared with the owner driven directly, and takes
//! the heat of the fs-mbd damped network through heat routes: the "rigid
//! damping heat" pairing with both of its real owners, as a declared graph.
//! Run with `cargo test -p fs-system --features slices`.

// Exact float comparisons here assert bitwise agreement with the owners.
#![allow(clippy::float_cmp)]

use std::cell::Cell;
use std::collections::BTreeMap;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::ViscousNetwork;
use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
        network::{ViscousNetworkLimits, ViscousNetworkLink},
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads,
    viscous::network::{ViscousNetworkMotionControls, step_viscous_network_with_loads},
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, Proposal, Receipt};
use fs_system::{Ledger, RouteTemperature, SystemBuilder, heat_sink};
use fs_thermal::component::{
    ThermalTimeline, TimelineBuildError, TimelineRefusal, TimelineState,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, Reservoir, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Constant heat capacity above a datum: `T = 300 K + (e - datum) / cp`. A
/// synthetic law for the wiring, not a material. When a gate is armed, the
/// next temperature query requests cancellation on it; when told to refuse,
/// every query is refused.
struct Sensible {
    datum_j_kg: f64,
    tripwire: Cell<Option<&'static CancelGate>>,
    refuse: Cell<bool>,
}

const LAW_REFUSAL: DepositionError = DepositionError::InvalidInput("the test law refuses");

impl ThermalLaw for Sensible {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        if self.refuse.get() {
            return Err(LAW_REFUSAL);
        }
        if let Some(gate) = self.tripwire.get() {
            gate.request();
        }
        Ok(300.0 + (specific_h_j_kg - self.datum_j_kg) / CP)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}

const CP: f64 = 500.0;
const MASS: f64 = 0.01;
/// One cell's heat capacity [J/K].
const CAPACITY: f64 = MASS * CP;
const LINK_W_K: f64 = 0.2;
const LINK_AREA_M2: f64 = 1.0e-3;
const SURFACE_M2: f64 = 0.01;
const A: CellId = CellId(10);
const B: CellId = CellId(20);

/// The owner's documented explicit bound (`enthalpy.rs`, module docs):
/// `dt <= 0.8 m / (max(dT/de) * sum(outgoing conductance))`.
fn explicit_bound_s(conductance_w_k: f64) -> f64 {
    0.8 * CAPACITY / conductance_w_k
}

fn config(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1.0e-9,
    }
}

fn insulated() -> Boundary {
    Boundary {
        ambient_k: 293.15,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}

/// Two cells at 310 K and 295 K, born at `birth_s` and joined by one link,
/// with `extra` events after that. Each cell has its own material, whose
/// law has its own enthalpy datum.
struct Scene {
    link_w_k: f64,
    boundaries: [Boundary; 2],
    extra: Vec<TimedEvent>,
    datums_j_kg: [f64; 2],
    birth_s: f64,
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            link_w_k: LINK_W_K,
            boundaries: [insulated(); 2],
            extra: Vec::new(),
            datums_j_kg: [0.0; 2],
            birth_s: 0.0,
        }
    }
}

impl Scene {
    /// A cell of material `material` (0 or 1).
    fn cell(
        &self,
        id: CellId,
        material: usize,
        temperature_k: f64,
        boundary: Boundary,
    ) -> DepositedCell {
        let row = id.0 as f64;
        DepositedCell {
            id,
            segment_id: id.0,
            body_id: id.0,
            material_id: MaterialId(material as u64),
            mass_kg: MASS,
            specific_h_j_kg: self.datums_j_kg[material] + CP * (temperature_k - 300.0),
            geometry: CellGeometry {
                start_m: [0.0, row, 0.0],
                end_m: [0.1, row, 0.0],
                volume_m3: 1.0e-6,
                surface_area_m2: SURFACE_M2,
            },
            boundary,
        }
    }

    /// The timeline and the first cell's law. The laws live as long as the
    /// test process: the kernel keeps component states as `'static` values.
    fn build(&self) -> (DepositionSimulation<'static>, &'static Sensible) {
        let [law, other] = self.datums_j_kg.map(|datum_j_kg| {
            let law: &'static Sensible = Box::leak(Box::new(Sensible {
                datum_j_kg,
                tripwire: Cell::new(None),
                refuse: Cell::new(false),
            }));
            law
        });
        let at = |event| TimedEvent {
            time_s: self.birth_s,
            event,
        };
        let mut schedule = vec![
            at(Event::Birth(self.cell(A, 0, 310.0, self.boundaries[0]))),
            at(Event::Birth(self.cell(B, 1, 295.0, self.boundaries[1]))),
            at(Event::Contact(ThermalLink {
                id: 7,
                a: A,
                b: B,
                area_m2: LINK_AREA_M2,
                conductance_w_k: self.link_w_k,
            })),
        ];
        schedule.extend(self.extra.iter().cloned());
        let simulation = DepositionSimulation::new(
            schedule,
            &[
                (MaterialId(0), law as &dyn ThermalLaw),
                (MaterialId(1), other as &dyn ThermalLaw),
            ],
            DepositionBudget {
                max_cells: 8,
                max_links: 8,
                max_events: 16,
                max_substeps: 100_000,
            },
        )
        .expect("timeline");
        (simulation, law)
    }

    fn timeline(&self) -> DepositionSimulation<'static> {
        self.build().0
    }
}

fn boundary_event(time_s: f64, cell: CellId, boundary: Boundary) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Boundary(cell, boundary),
    }
}

fn temperature(state: &TimelineState<'_>, id: CellId) -> f64 {
    state.simulation().temperature_k(id).expect("temperature")
}

/// `state` is at the owner's `native` time with its temperatures. The
/// wrapper hands the owner each port's energy and the direct run gives it
/// rates, so the two agree to rounding, not bit for bit.
fn same_cells(state: &TimelineState<'_>, native: &DepositionSimulation<'_>) {
    assert_eq!(state.simulation().time_s(), native.time_s());
    for id in [A, B] {
        let (wrapped, owner) = (temperature(state, id), native.temperature_k(id).expect("T"));
        assert!((wrapped - owner).abs() <= 1.0e-11, "{id:?}: {wrapped} vs {owner}");
    }
}

/// One step of `component` from `state`, outside the kernel.
fn propose(
    component: &ThermalTimeline<'static>,
    state: &TimelineState<'static>,
    dt_s: f64,
    inputs: &[f64],
    gate: &CancelGate,
) -> Result<(TimelineState<'static>, Receipt, Vec<f64>), TimelineRefusal> {
    let mut next = state.clone();
    let mut outputs = vec![f64::NAN; component.ports().len()];
    let proposal = with_cx(gate, |cx| {
        component.propose(state, dt_s, inputs, &mut next, &mut outputs, cx)
    })?;
    let Proposal::Step(receipt) = proposal else {
        panic!("the timeline takes the step or refuses");
    };
    Ok((next, receipt, outputs))
}

/// The conformance check on one case, which must step.
fn conforms(
    component: &ThermalTimeline<'static>,
    state: &TimelineState<'static>,
    stale: &TimelineState<'static>,
    dt_s: f64,
    inputs: &[f64],
) -> (Receipt, f64, Vec<f64>) {
    let gate = CancelGate::new_clock_free();
    let case = Case {
        state,
        stale,
        dt_s,
        inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(component, &case, cx)).expect("conforms");
    let Observed::Stepped {
        receipt,
        audit,
        outputs,
    } = observed
    else {
        panic!("the timeline takes the step");
    };
    (receipt, audit.port_energy_in_j, outputs.clone())
}

/// A state at another time and other enthalpies, for the stale slot.
fn later(
    component: &ThermalTimeline<'static>,
    simulation: &DepositionSimulation<'static>,
    by_s: f64,
) -> TimelineState<'static> {
    let (advanced, _) = simulation
        .advance_to(simulation.time_s() + by_s, config(0.05), || Ok(()))
        .expect("advance");
    component.start(&advanced).expect("stale state")
}

/// The owner driven directly with held heat, through its own observer
/// interface; it also records the first step bound the owner proposes.
#[derive(Clone)]
struct HeldSources {
    rates_w: BTreeMap<CellId, f64>,
    first_bound_s: Option<f64>,
}

impl HeldSources {
    fn new(rates_w: impl IntoIterator<Item = (CellId, f64)>) -> Self {
        Self {
            rates_w: rates_w.into_iter().collect(),
            first_bound_s: None,
        }
    }
}

impl DepositionObserver for HeldSources {
    type Error = DepositionError;

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        Ok(self.rates_w.clone())
    }

    fn limit_timestep_s(
        &mut self,
        _: &DepositionSimulation<'_>,
        proposed_s: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, DepositionError> {
        self.first_bound_s.get_or_insert(proposed_s);
        Ok(proposed_s)
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), DepositionError> {
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        Ok(())
    }
}

#[test]
fn the_timeline_conforms_and_its_step_is_the_explicit_rule_of_its_owner() {
    // The second scene stores +1e7 J in one cell and -1e7 J in the other
    // against their laws' datums. The total is small and its terms are not:
    // the receipt's tolerance has to cover the rounding of the terms.
    for datum_j_kg in [0.0, 1.0e9] {
        let scene = Scene {
            datums_j_kg: [datum_j_kg, -datum_j_kg],
            ..Scene::default()
        };
        let simulation = scene.timeline();
        let component = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
        let start = component.start(&simulation).expect("start");
        assert_eq!((temperature(&start, A), temperature(&start, B)), (310.0, 295.0));
        assert_eq!(component.ports().len(), 2);
        let stale = later(&component, &simulation, 1.0);
        assert_ne!(start, stale);
        for (dt_s, rates) in [(0.01, [2.0, -0.5]), (0.05, [0.0, 0.0]), (0.002, [-30.0, 45.0])] {
            let (receipt, port_energy_in_j, outputs) =
                conforms(&component, &start, &stale, dt_s, &rates);
            // Explicit Euler by hand: each cell gains its port's heat and
            // the link's heat at the start temperatures.
            let through = LINK_W_K * (310.0 - 295.0);
            let after_a = 310.0 + dt_s * (rates[0] - through) / CAPACITY;
            let after_b = 295.0 + dt_s * (rates[1] + through) / CAPACITY;
            // At the datum the enthalpy carries 1e9 J/kg, so a temperature
            // is known to one ulp of that over cp.
            let resolution_k = 1.0e-12 + 4.0 * f64::EPSILON * datum_j_kg / CP;
            assert!((outputs[0] - 0.5 * (310.0 + after_a)).abs() <= resolution_k, "{}", outputs[0]);
            assert!((outputs[1] - 0.5 * (295.0 + after_b)).abs() <= resolution_k, "{}", outputs[1]);
            let delivered = dt_s * (rates[0] + rates[1]);
            let size = dt_s * (rates[0].abs() + rates[1].abs());
            assert!((port_energy_in_j - delivered).abs() <= 1.0e-14 * size);
            // An enthalpy of 5 kJ/kg above its datum is known to an ulp.
            let rounding_j =
                1.0e-12 * size + 8.0 * f64::EPSILON * MASS * (datum_j_kg + CP * 10.0);
            assert!((receipt.stored_change_j - delivered).abs() <= rounding_j);
            assert!(receipt.numerical_loss_j.abs() <= rounding_j);
        }
    }
}

#[test]
fn the_step_limit_is_the_owners_explicit_bound_and_a_longer_step_is_its_substeps() {
    let gate = CancelGate::new_clock_free();
    // A generous maximum step: the owner's conduction bound is what limits.
    let component = ThermalTimeline::new(vec![A, B], false, config(1.0e3)).expect("component");
    let simulation = Scene::default().timeline();
    let start = component.start(&simulation).expect("start");
    let rates = [1.5, -0.25];
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);

    // The bound is the owner's rule in closed form, and what the owner
    // itself proposes when driven directly.
    let bound = component.limit_dt_s(&start, 1.0e3).expect("limit");
    assert!((bound - explicit_bound_s(LINK_W_K)).abs() <= 1.0e-12 * bound, "{bound}");
    let origin = start.simulation();
    let (_, observer, _) = origin
        .advance_observed(&held, 1.0e3, config(1.0e3), || Ok(()))
        .expect("owner");
    assert_eq!(observer.first_bound_s, Some(bound));
    assert_eq!(component.limit_dt_s(&start, bound).expect("limit"), bound);
    assert_eq!(component.limit_dt_s(&start, 0.25 * bound).expect("limit"), 0.25 * bound);
    // The owner's maximum step is part of its bound.
    let capped = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    assert_eq!(capped.limit_dt_s(&start, 1.0e3).expect("limit"), 0.05);
    assert_eq!(component.next_event_s(&start), None);

    // One substep at the bound, then two and a half bounds in three
    // substeps: the owner's own advance with the rates held.
    for (dt_s, substeps) in [(bound, 1), (2.5 * bound, 3)] {
        let (next, receipt, outputs) =
            propose(&component, &start, dt_s, &rates, &gate).expect("step");
        let (native, _, audit) = origin
            .advance_observed(&held, dt_s, config(1.0e3), || Ok(()))
            .expect("owner");
        assert_eq!(audit.substeps, substeps);
        assert_eq!(next.simulation().time_s(), dt_s);
        same_cells(&next, &native);
        let delivered = dt_s * (rates[0] + rates[1]);
        let rounding_j = 1.0e-12 * delivered;
        assert!((receipt.stored_change_j - audit.stored_energy_change_j).abs() <= rounding_j);
        assert!((receipt.stored_change_j - delivered).abs() <= rounding_j);
        assert!(receipt.numerical_loss_j.abs() <= rounding_j);
        for (port, id) in [A, B].into_iter().enumerate() {
            let mean = 0.5 * (temperature(&start, id) + native.temperature_k(id).expect("T"));
            assert!((outputs[port] - mean).abs() <= 1.0e-11, "{} vs {mean}", outputs[port]);
        }
    }
}

#[test]
fn a_state_is_what_it_holds() {
    let gate = CancelGate::new_clock_free();
    let component = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    let start = component.start(&Scene::default().timeline()).expect("start");
    let step = |rates: [f64; 2]| propose(&component, &start, 0.01, &rates, &gate).expect("step").0;
    // The same heat, the same state; other heat, another.
    assert_eq!(step([1.5, -0.25]), step([1.5, -0.25]));
    assert_ne!(step([1.5, -0.25]), step([0.0, 0.0]));
    // The same cells on another link, or with an event ahead, are other states.
    let stiffer = Scene {
        link_w_k: 0.3,
        ..Scene::default()
    };
    assert_ne!(start, component.start(&stiffer.timeline()).expect("start"));
    let relabel = Boundary {
        ambient_k: 280.0,
        ..insulated()
    };
    let eventful = Scene {
        extra: vec![boundary_event(5.0, A, relabel)],
        ..Scene::default()
    };
    assert_ne!(start, component.start(&eventful.timeline()).expect("start"));
}

#[test]
fn a_step_the_clock_cannot_represent_keeps_the_energy_and_one_it_cannot_take_is_refused() {
    // At 2^30 s one ulp of the owner's clock is 2^-22 s, about 0.24 us.
    let birth_s = 1_073_741_824.0_f64;
    let ulp_s = birth_s * f64::EPSILON;
    let scene = Scene {
        birth_s,
        ..Scene::default()
    };
    let (waiting, _) = scene
        .timeline()
        .advance_to(birth_s, config(0.05), || Ok(()))
        .expect("advance to the births");
    let gate = CancelGate::new_clock_free();
    let rates = [1000.0, -250.0];
    // (the owner's maximum step, the kernel's step, the clock's interval in
    // ulps, the owner's substeps):
    // - 1 us ends four ulps on, 4.6 % short, in one substep the target clips;
    // - 1.1 us with a 1.1 us bound ends five ulps on, 8 % long, in one
    //   substep the bound limits, which integrates 1.1 us;
    // - 10 us with a 0.3 us bound takes 42 substeps of one ulp each, every
    //   one integrating 0.3 us but the last: 26 % more than 10 us.
    for (max_dt_s, dt_s, ulps, substeps) in
        [(0.05, 1.0e-6, 4.0, 1_u32), (1.1e-6, 1.1e-6, 5.0, 1), (3.0e-7, 1.0e-5, 42.0, 42)]
    {
        let component =
            ThermalTimeline::new(vec![A, B], false, config(max_dt_s)).expect("component");
        let start = component.start(&waiting).expect("start");
        assert_eq!(start.simulation().time_s(), birth_s);
        let stale = later(&component, &waiting, 1.0);
        let (receipt, _, _) = conforms(&component, &start, &stale, dt_s, &rates);
        let (next, _, _) = propose(&component, &start, dt_s, &rates, &gate).expect("step");
        assert_eq!(next.simulation().time_s() - birth_s, ulps * ulp_s);
        let (_, audit) = start
            .simulation()
            .advance_to(next.simulation().time_s(), config(max_dt_s), || Ok(()))
            .expect("owner");
        assert_eq!(audit.substeps, substeps as usize);
        // The cells received the kernel's energy, whatever the owner's
        // durations sum to. An enthalpy of 5 kJ/kg is known to an ulp, at
        // every substep.
        let delivered = dt_s * (rates[0] + rates[1]);
        let quantum_j = 8.0 * f64::EPSILON * MASS * (CP * 10.0) * f64::from(substeps);
        assert!(quantum_j < 1.0e-8 * delivered);
        let stored = receipt.stored_change_j;
        assert!((stored - delivered).abs() <= quantum_j, "{dt_s}: {stored} vs {delivered}");
        assert!(receipt.numerical_loss_j.abs() <= quantum_j, "{}", receipt.numerical_loss_j);
        let warmed_k = temperature(&next, A) + temperature(&next, B) - (310.0 + 295.0);
        assert!((CAPACITY * warmed_k - delivered).abs() <= 1.0e-6 * delivered, "{warmed_k}");

        // A step under half an ulp does not move the clock: it is refused,
        // and the step limit has no bound to report for it.
        assert_eq!(
            propose(&component, &start, 0.25 * ulp_s, &rates, &gate).unwrap_err(),
            TimelineRefusal::Interval { dt_s: 0.25 * ulp_s }
        );
        assert_eq!(component.limit_dt_s(&start, 0.25 * ulp_s).expect("limit"), 0.25 * ulp_s);
    }
}

#[test]
fn a_step_cut_to_an_event_ends_on_the_event() {
    let gate = CancelGate::new_clock_free();
    // Find a start time and an event time for which `t + (e - t)` is not
    // `e` in floating point.
    let (start_s, event_s) = (1_u32..400)
        .flat_map(|k| (0_u32..40).map(move |j| (f64::from(k) * 0.013, 1.0 + f64::from(j) * 0.1)))
        .find(|&(t, e)| t < 0.5 * e && t + (e - t) != e)
        .expect("a pair that does not round-trip");
    let scene = Scene {
        extra: vec![boundary_event(
            event_s,
            A,
            Boundary {
                ambient_k: 280.0,
                ..insulated()
            },
        )],
        ..Scene::default()
    };
    let component = ThermalTimeline::new(vec![A, B], false, config(1.0e3)).expect("component");
    let (at_start, _) = scene
        .timeline()
        .advance_to(start_s, config(1.0e3), || Ok(()))
        .expect("advance");
    let state = component.start(&at_start).expect("state");
    let ahead_s = component.next_event_s(&state).expect("an event ahead");
    assert_eq!(ahead_s, event_s - start_s);

    let rates = [0.75, 0.0];
    let (next, _, _) = propose(&component, &state, ahead_s, &rates, &gate).expect("step");
    assert_eq!(next.simulation().time_s(), event_s);
    // The event was applied and nothing is pending.
    assert_eq!(next.simulation().cells()[&A].boundary.ambient_k, 280.0);
    assert_eq!(component.next_event_s(&next), None);
    // The owner driven directly reaches the event in one substep, with the
    // same cells.
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);
    let (native, _, audit) = at_start
        .advance_observed(&held, event_s, config(1.0e3), || Ok(()))
        .expect("owner");
    assert_eq!(audit.substeps, 1);
    same_cells(&next, &native);
    // A shorter step is not moved onto the event.
    let (short, _, _) = propose(&component, &state, 0.5 * ahead_s, &rates, &gate).expect("step");
    assert_eq!(short.simulation().time_s(), start_s + 0.5 * ahead_s);
    assert_eq!(short.simulation().cells()[&A].boundary.ambient_k, 293.15);
}

#[test]
fn the_timeline_refuses_what_it_cannot_carry() {
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        ThermalTimeline::new(Vec::new(), false, config(0.05)).unwrap_err(),
        TimelineBuildError::NoPorts
    );
    assert_eq!(
        ThermalTimeline::new(vec![A, B, A], false, config(0.05)).unwrap_err(),
        TimelineBuildError::Duplicate { cell: 10 }
    );
    // A port on a cell the timeline does not have.
    let missing = ThermalTimeline::new(vec![A, CellId(99)], false, config(0.05)).expect("ports");
    assert_eq!(
        missing.start(&Scene::default().timeline()).unwrap_err(),
        TimelineRefusal::InactiveCell { cell: 99 }
    );

    let component = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    let start = component.start(&Scene::default().timeline()).expect("start");
    assert!(matches!(
        propose(&component, &start, 0.01, &[f64::NAN, 0.0], &gate),
        Err(TimelineRefusal::Input { value }) if value.is_nan()
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        propose(&component, &start, 0.01, &[0.0, 0.0], &cancelled).unwrap_err(),
        TimelineRefusal::Cancelled
    );

    // The owner's own refusal comes through unchanged, from the step limit
    // and from the step.
    let (simulation, law) = Scene::default().build();
    let start = component.start(&simulation).expect("start");
    law.refuse.set(true);
    assert_eq!(
        component.limit_dt_s(&start, 0.01).unwrap_err(),
        TimelineRefusal::Thermal(LAW_REFUSAL)
    );
    assert_eq!(
        propose(&component, &start, 0.01, &[0.0, 0.0], &gate).unwrap_err(),
        TimelineRefusal::Thermal(LAW_REFUSAL)
    );

    // A cell born after the start brings enthalpy with no port to carry it:
    // refused when the step ends on the birth (where the kernel cuts it)
    // and when it crosses it.
    let scene = Scene::default();
    let growing = Scene {
        extra: vec![TimedEvent {
            time_s: 0.005,
            event: Event::Birth(scene.cell(CellId(30), 1, 320.0, insulated())),
        }],
        ..Scene::default()
    };
    let start = component.start(&growing.timeline()).expect("start");
    assert_eq!(component.next_event_s(&start), Some(0.005));
    for dt_s in [0.005, 0.02] {
        assert_eq!(
            propose(&component, &start, dt_s, &[0.0, 0.0], &gate).unwrap_err(),
            TimelineRefusal::Birth {
                cell: 30,
                energy_j: MASS * (CP * 20.0)
            }
        );
    }
}

#[test]
fn cancellation_inside_the_owners_step_is_seen() {
    let gate: &'static CancelGate = Box::leak(Box::new(CancelGate::new_clock_free()));
    let component = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    let (simulation, law) = Scene::default().build();
    let start = component.start(&simulation).expect("start");
    // Armed now, the request arrives when the owner first evaluates a
    // temperature inside the step, after the wrapper's own entry poll.
    law.tripwire.set(Some(gate));
    assert_eq!(
        propose(&component, &start, 0.01, &[1.0, 0.0], gate).unwrap_err(),
        TimelineRefusal::Cancelled
    );
}

#[test]
fn boundary_heat_needs_the_environment_port_whatever_it_sums_to() {
    let gate = CancelGate::new_clock_free();
    let no_port = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    let film = Boundary {
        htc_w_m2_k: 25.0,
        ..insulated()
    };
    let source = |source_w| Boundary {
        source_w,
        ..insulated()
    };
    // A film on the second cell.
    let filmed = Scene {
        boundaries: [insulated(), film],
        ..Scene::default()
    };
    assert_eq!(
        no_port.start(&filmed.timeline()).unwrap_err(),
        TimelineRefusal::Boundary { cell: 20 }
    );
    // A state admitted by a timeline that has the port is refused by one
    // that has not.
    let with_port = ThermalTimeline::new(vec![A, B], true, config(0.05)).expect("component");
    let admitted = with_port.start(&filmed.timeline()).expect("start");
    assert_eq!(
        propose(&no_port, &admitted, 0.01, &[0.0, 0.0], &gate).unwrap_err(),
        TimelineRefusal::Boundary { cell: 20 }
    );
    // A bed under the first cell, at the cell's own temperature.
    let bedded = Scene {
        boundaries: [
            Boundary {
                reservoir: Some(Reservoir {
                    area_m2: 1.0e-3,
                    conductance_w_k: 0.5,
                    temperature_k: 310.0,
                }),
                ..insulated()
            },
            insulated(),
        ],
        ..Scene::default()
    };
    assert_eq!(
        no_port.start(&bedded.timeline()).unwrap_err(),
        TimelineRefusal::Boundary { cell: 10 }
    );
    // Authored sources that cancel: no net heat, and still off the wires.
    let opposed = Scene {
        boundaries: [source(5.0), source(-5.0)],
        ..Scene::default()
    };
    assert_eq!(
        no_port.start(&opposed.timeline()).unwrap_err(),
        TimelineRefusal::Boundary { cell: 10 }
    );
    // A film at ambient temperature exchanges nothing yet; it is refused
    // for what it is.
    let idle_film = Scene {
        boundaries: [
            Boundary {
                ambient_k: 310.0,
                ..film
            },
            insulated(),
        ],
        ..Scene::default()
    };
    assert_eq!(
        no_port.start(&idle_film.timeline()).unwrap_err(),
        TimelineRefusal::Boundary { cell: 10 }
    );
    // An event that gives a cell a film, at the end of a step or inside it.
    let later_film = Scene {
        extra: vec![boundary_event(0.02, B, film)],
        ..Scene::default()
    };
    let start = no_port.start(&later_film.timeline()).expect("insulated at the start");
    for dt_s in [0.02, 0.05] {
        assert_eq!(
            propose(&no_port, &start, dt_s, &[0.0, 0.0], &gate).unwrap_err(),
            TimelineRefusal::Boundary { cell: 20 }
        );
    }
    // An event that only changes an insulated cell's ambient is no exchange.
    let relabelled = Scene {
        extra: vec![boundary_event(
            0.02,
            B,
            Boundary {
                ambient_k: 250.0,
                ..insulated()
            },
        )],
        ..Scene::default()
    };
    let start = no_port.start(&relabelled.timeline()).expect("start");
    let (next, receipt, _) = propose(&no_port, &start, 0.05, &[0.0, 0.0], &gate).expect("step");
    assert_eq!(next.simulation().time_s(), 0.05);
    assert!(receipt.stored_change_j.abs() <= 1.0e-12);
}

#[test]
fn the_environment_port_returns_the_boundaries_heat() {
    let film = Boundary {
        htc_w_m2_k: 25.0,
        ..insulated()
    };
    let source = |source_w| Boundary {
        source_w,
        ..insulated()
    };
    // With the port, boundaries conform and their heat is the port's
    // output: the film's rate at the start temperature over the surface the
    // contact leaves free, plus the authored source.
    let with_port = ThermalTimeline::new(vec![A, B], true, config(0.05)).expect("component");
    let exposed = Scene {
        boundaries: [film, source(-0.4)],
        ..Scene::default()
    };
    let simulation = exposed.timeline();
    let start = with_port.start(&simulation).expect("start");
    let stale = later(&with_port, &simulation, 1.0);
    let (dt_s, inputs) = (0.01, [2.0, -0.5, 293.15]);
    let (receipt, port_energy_in_j, outputs) = conforms(&with_port, &start, &stale, dt_s, &inputs);
    let film_w = 25.0 * (SURFACE_M2 - LINK_AREA_M2) * (293.15 - 310.0);
    let environment_w = film_w - 0.4;
    assert!((outputs[2] - environment_w).abs() <= 1.0e-12 * environment_w.abs(), "{}", outputs[2]);
    let entered = dt_s * (inputs[0] + inputs[1] + environment_w);
    assert!((port_energy_in_j - entered).abs() <= 1.0e-12 * entered.abs());
    assert!((receipt.stored_change_j - entered).abs() <= 1.0e-12 * entered.abs());
    assert!(receipt.numerical_loss_j.abs() <= 1.0e-12 * entered.abs());

    // A timeline with the environment port alone.
    let only = ThermalTimeline::new(Vec::new(), true, config(0.05)).expect("component");
    let start = only.start(&simulation).expect("start");
    let stale = later(&only, &simulation, 1.0);
    let (_, _, outputs) = conforms(&only, &start, &stale, dt_s, &[293.15]);
    assert!((outputs[0] - environment_w).abs() <= 1.0e-12 * environment_w.abs(), "{}", outputs[0]);

    // While no boundary exchanges, the port returns exactly nothing, whatever
    // the heat on the cell ports.
    let gate = CancelGate::new_clock_free();
    let quiet = with_port.start(&Scene::default().timeline()).expect("start");
    for rates in [[1.0e5, -3.0e4], [0.1, 0.7], [-123.456, 7.0e3], [3.0, 1.0 / 3.0]] {
        let inputs = [rates[0], rates[1], 293.15];
        let (_, _, outputs) = propose(&with_port, &quiet, 0.037, &inputs, &gate).expect("step");
        assert_eq!(outputs[2], 0.0);
    }
}

#[test]
fn boundary_heat_is_on_the_environment_wire_and_reaches_a_named_sink() {
    let gate = CancelGate::new_clock_free();
    let film = |htc_w_m2_k| Boundary {
        htc_w_m2_k,
        ..insulated()
    };
    // Fan changes at 0.4 s and 0.7 s: the kernel must not step across them.
    let scene = Scene {
        boundaries: [film(25.0), insulated()],
        extra: vec![boundary_event(0.4, A, film(5.0)), boundary_event(0.7, A, film(40.0))],
        ..Scene::default()
    };
    let component = ThermalTimeline::new(vec![A, B], true, config(0.05)).expect("component");
    let start = component.start(&scene.timeline()).expect("start");
    let origin = start.simulation().clone();
    let mut builder = SystemBuilder::new();
    let cells = builder.add("cells", component, start).expect("cells");
    let port = |builder: &SystemBuilder, name: &str| builder.port(cells, name).expect("port");
    builder
        .boundary(port(&builder, "cell10"), Box::new(|_, _, out: &mut [f64]| out[0] = 0.0))
        .expect("idle");
    builder
        .boundary(port(&builder, "cell20"), Box::new(|_, _, out: &mut [f64]| out[0] = 1.5))
        .expect("heated");
    builder
        .heat_route(
            port(&builder, "environment"),
            &[],
            Some(heat_sink("room", 293.15).expect("sink")),
            RouteTemperature::Fixed(293.15),
        )
        .expect("route");
    let mut system = builder.seal().expect("sealed");
    let (mut steps, mut cut_steps, mut landed) = (0_u32, 0_u32, [false; 2]);
    with_cx(&gate, |cx| {
        while system.time_s() < 1.0 - 1.0e-12 {
            let report = system.step(0.03, cx).expect("step");
            cut_steps += u32::from(report.dt_s < 0.029);
            landed[0] |= (report.time_s - 0.4).abs() <= 1.0e-12;
            landed[1] |= (report.time_s - 0.7).abs() <= 1.0e-12;
            steps += 1;
        }
    });
    // 0.03 does not divide 0.4: thirteen steps reach 0.39 s, the fourteenth
    // is cut to 0.01 s to end on the first change, ten more end on the
    // second (the tenth cut by rounding only), and ten more on 1 s. Stepping
    // across the changes would end at 1.02 s with no step cut.
    assert_eq!(landed, [true; 2], "a step ends on each fan change");
    assert_eq!((steps, cut_steps), (34, 1));
    assert!((system.time_s() - 1.0).abs() <= 1.0e-12, "{}", system.time_s());

    // The owner alone over the same second, with the same held heat.
    let held = HeldSources::new([(A, 0.0), (B, 1.5)]);
    let (native, _, audit) = origin
        .advance_observed(&held, system.time_s(), config(0.03), || Ok(()))
        .expect("owner");
    let reached = system.state::<TimelineState<'static>>(cells).expect("state");
    // The timeline's clock is its own sum of the steps; it agrees with the
    // kernel's to rounding.
    assert!((reached.simulation().time_s() - system.time_s()).abs() <= 1.0e-12);
    for id in [A, B] {
        let (kernel, owner) = (temperature(reached, id), native.temperature_k(id).expect("T"));
        assert!((kernel - owner).abs() <= 1.0e-10, "{kernel} vs {owner}");
    }
    // What the boundary took is in the sink: the owner's own external heat
    // less what the port delivered.
    let ledger = system.ledger();
    let delivered = 1.5 * system.time_s();
    let through_film = audit.external_heat_j - delivered;
    assert!(through_film < -0.1, "the film cools cell A: {through_film} J");
    let room = ledger.sink_heat_j();
    assert!((room + through_film).abs() <= 1.0e-9 * through_film.abs(), "{room} vs {through_film}");
    let stored = ledger.stored_change_j();
    assert!((stored - audit.stored_energy_change_j).abs() <= 1.0e-9 * delivered, "{stored}");
}

// ---------------------------------------------------------------------------
// The damping-heat pairing with both real owners.
// ---------------------------------------------------------------------------

fn controls() -> ViscousNetworkMotionControls {
    ViscousNetworkMotionControls {
        source: ViscousConnectorControls {
            impulse: InternalWrenchImpulseControls {
                angular_origin_world_m: Vec3::new(0.2, -0.3, 0.1),
                linear_tolerance_ns: 3e-10,
                angular_tolerance_nms: 3e-10,
                work_tolerance_j: 3e-10,
            },
            linear_speed_tolerance_m_s: 3e-10,
            angular_speed_tolerance_rad_s: 3e-10,
        },
        source_substeps: 1,
        limits: ViscousNetworkLimits {
            maximum_bodies: 32,
            maximum_connectors: 64,
            maximum_pair_relaxations: 100_000,
            maximum_live_bytes: 1_000_000,
        },
    }
}

fn axial(bodies: [usize; 2], coefficient: f64) -> ViscousNetworkLink {
    let mut c = [[0.0; 6]; 6];
    c[0][0] = coefficient;
    ViscousNetworkLink {
        bodies,
        arms_body_m: [Vec3::ZERO; 2],
        damping: PassiveBushingMatrix::try_new(
            c,
            BushingCoordinateScales {
                translation_m: 1.0,
                rotation_rad: 1.0,
            },
        )
        .expect("passive"),
    }
}

fn bodies() -> (Vec<MassProperties>, Vec<ViscousNetworkLink>, Vec<RigidBodyState>) {
    let properties =
        vec![MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).expect("mass"); 3];
    let states = [1.2, -0.8, 0.2]
        .into_iter()
        .enumerate()
        .map(|(index, speed)| {
            RigidBodyState::new(
                Pose::new(Vec3::new(index as f64 - 1.0, 0.0, 0.0), UnitQuaternion::IDENTITY)
                    .expect("pose"),
                Vec3::new(speed * 2.0, 0.0, 0.0),
                Vec3::ZERO,
            )
            .expect("state")
        })
        .collect();
    (properties, vec![axial([0, 1], 0.6), axial([1, 2], 0.9)], states)
}

const FRACTION_A: f64 = 0.6;
const FRACTION_B: f64 = 0.3;
/// The link between the cells in the pairing [W/K]: stiff enough that the
/// thermal owner's explicit bound, 10 ms, is the system's clock.
const STIFF_LINK_W_K: f64 = 400.0;

#[test]
fn damper_heat_reaches_the_timelines_cells_on_the_thermal_clock() {
    let gate = CancelGate::new_clock_free();
    let (properties, links, start_bodies) = bodies();
    let network =
        ViscousNetwork::new(properties.clone(), links.clone(), Gravity::ZERO, controls())
            .expect("network");
    let scene = Scene {
        link_w_k: STIFF_LINK_W_K,
        ..Scene::default()
    };
    let thermal = ThermalTimeline::new(vec![A, B], false, config(1.0)).expect("timeline");
    let start_cells = thermal.start(&scene.timeline()).expect("start");
    let clock = thermal.limit_dt_s(&start_cells, 1.0).expect("bound");
    assert!((clock - explicit_bound_s(STIFF_LINK_W_K)).abs() <= 1.0e-12 * clock, "{clock}");

    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, start_bodies.clone()).expect("network");
    let cells = builder.add("cells", thermal, start_cells.clone()).expect("cells");
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).expect("port");
    for body in 0..3 {
        let force = port(&builder, moving, &format!("force{body}"));
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .expect("free");
    }
    let (cell_a, cell_b) = (port(&builder, cells, "cell10"), port(&builder, cells, "cell20"));
    builder
        .heat_route(
            port(&builder, moving, "heat0"),
            &[(cell_a, FRACTION_A), (cell_b, FRACTION_B)],
            Some(heat_sink("unallocated", 293.15).expect("sink")),
            RouteTemperature::SinkEquivalent,
        )
        .expect("route 0");
    builder
        .heat_route(
            port(&builder, moving, "heat1"),
            &[(cell_b, 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .expect("route 1");
    let mut system = builder.seal().expect("sealed");
    for name in ["heat0", "heat1"] {
        let heat = system.port(moving, name).expect("heat port");
        assert!(system.guess(heat, &[300.0]));
    }

    // The same protocol by hand, as the hand-written coupling does it: the
    // thermal owner picks each interval, mechanics advances once over it,
    // and connector heat enters the cells as mean power over that interval.
    let idle = [Wrench::ZERO; 3];
    let mut twin_bodies = start_bodies;
    let mut twin_cells = start_cells.into_simulation();
    let steps = 300_usize;
    let (mut heat_0, mut heat_1) = (0.0_f64, 0.0_f64);
    with_cx(&gate, |cx| {
        for step in 0..steps {
            // Ask for more than the bound: the thermal component cuts it.
            let dt_s = system.step(1.0, cx).expect("kernel step").dt_s;
            assert!((dt_s - clock).abs() <= 1.0e-12 * clock, "step {step}: {dt_s} vs {clock}");
            let motion = step_viscous_network_with_loads(
                &twin_bodies,
                &properties,
                &links,
                HeldExternalLoads {
                    wrenches: &idle,
                    gravity: Gravity::ZERO,
                },
                dt_s,
                controls(),
                || false,
            )
            .expect("owner motion");
            let (d0, d1) = (
                motion.connector_dissipation_j[0],
                motion.connector_dissipation_j[1],
            );
            heat_0 += d0;
            heat_1 += d1;
            let held = HeldSources::new([
                (A, FRACTION_A * d0 / dt_s),
                (B, FRACTION_B * d0 / dt_s + d1 / dt_s),
            ]);
            let target_s = twin_cells.time_s() + dt_s;
            let (next, _, audit) = twin_cells
                .advance_observed(&held, target_s, config(1.0), || Ok(()))
                .expect("owner heat");
            assert_eq!(audit.substeps, 1, "step {step}: the interval is one thermal substep");
            twin_cells = next;
            twin_bodies = motion.into_states_after();
        }
    });

    // The motion does not depend on where its heat goes: bit for bit.
    let reached_bodies = system.state::<Vec<RigidBodyState>>(moving).expect("bodies");
    assert_eq!(reached_bodies, &twin_bodies);
    // The cells took the same heat on the same clock; the kernel delivers
    // fraction * rate where the hand protocol delivers fraction * heat / dt.
    let reached_cells = system.state::<TimelineState<'static>>(cells).expect("cells");
    for id in [A, B] {
        let (kernel, hand) = (temperature(reached_cells, id), twin_cells.temperature_k(id).expect("T"));
        assert!((kernel - hand).abs() <= 1.0e-11, "{id:?}: {kernel} vs {hand}");
    }
    assert!((system.time_s() - twin_cells.time_s()).abs() <= 1.0e-12 * twin_cells.time_s());
    // The damper heat warmed the pair: without it the insulated cells would
    // only have shared what they started with.
    let warmed_k = temperature(reached_cells, A) + temperature(reached_cells, B) - (310.0 + 295.0);
    let kept_j = FRACTION_A * heat_0 + FRACTION_B * heat_0 + heat_1;
    assert!((warmed_k - kept_j / CAPACITY).abs() <= 1.0e-9, "{warmed_k} K for {kept_j} J");

    every_joule_is_in_the_cells_or_the_sink(system.ledger(), heat_0, heat_1);
    assert_eq!(system.ledger().steps, steps as u64);
}

/// Connector heat (`heat_0` from the first damper, `heat_1` from the second,
/// measured on the owner network) is in the cells or in the named sink.
fn every_joule_is_in_the_cells_or_the_sink(ledger: &Ledger, heat_0: f64, heat_1: f64) {
    let total = heat_0 + heat_1;
    assert!(total > 0.5, "enough heat to measure: {total} J");
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
    let warmed = account("cells").stored_change_j.value();
    let lost = ledger.sink_heat_j();
    let unallocated = (1.0 - FRACTION_A - FRACTION_B) * heat_0;
    assert!((lost - unallocated).abs() <= 1.0e-11 * total, "{lost} vs {unallocated}");
    assert!((warmed - (total - unallocated)).abs() <= 1.0e-11 * total, "{warmed}");
    assert!((account("network").dissipated_j.value() - total).abs() <= 1.0e-11 * total);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * total);
}
