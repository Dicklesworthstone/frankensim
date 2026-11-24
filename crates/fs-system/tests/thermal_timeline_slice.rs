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
use fs_system::elements::Conductance;
use fs_system::{Ledger, NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{
    ImplicitScheme, ThermalTimeline, TimelineBuildError, TimelineRefusal, TimelineState,
};
use fs_thermal::enthalpy::conduction::implicit::coupled_source::{
    ImplicitCoupledHeatSourceLaw, ImplicitCoupledHeatSourceState,
};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
    ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionReport,
    ImplicitConductionState, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::implicit::source::{
    ImplicitHeatSourceLaw, ImplicitHeatSourceState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
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
    /// For the implicit graph: conductivity `1 + slope (T - 300 K)` W/(m K).
    conductivity_slope_per_k: f64,
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

    /// `cp ln(T1 / T0)`, which the implicit owner's entropy account needs.
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let before_k = self.temperature_k(initial_h_j_kg)?;
        let after_k = self.temperature_k(final_h_j_kg)?;
        Ok(Some(CP * ((after_k - before_k) / before_k).ln_1p()))
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
    /// Whether the two cells are in contact at all. Without the link they
    /// exchange heat only through what is connected to their ports.
    linked: bool,
    /// With `linked`: the conduction graph names the contact as carried
    /// outside. It still covers its area of both cells' faces, and its heat
    /// is whatever is connected to the cells' ports.
    outside: bool,
    boundaries: [Boundary; 2],
    extra: Vec<TimedEvent>,
    datums_j_kg: [f64; 2],
    birth_s: f64,
    conductivity_slope_per_k: f64,
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            link_w_k: LINK_W_K,
            linked: true,
            outside: false,
            boundaries: [insulated(); 2],
            extra: Vec::new(),
            datums_j_kg: [0.0; 2],
            birth_s: 0.0,
            conductivity_slope_per_k: 0.0,
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
        let (simulation, [law, _]) = self.build_laws();
        (simulation, law)
    }

    fn build_laws(&self) -> (DepositionSimulation<'static>, [&'static Sensible; 2]) {
        let [law, other] = self.datums_j_kg.map(|datum_j_kg| {
            let law: &'static Sensible = Box::leak(Box::new(Sensible {
                datum_j_kg,
                tripwire: Cell::new(None),
                refuse: Cell::new(false),
                conductivity_slope_per_k: self.conductivity_slope_per_k,
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
        ];
        if self.linked {
            schedule.push(at(Event::Contact(ThermalLink {
                id: 7,
                a: A,
                b: B,
                area_m2: LINK_AREA_M2,
                conductance_w_k: self.link_w_k,
            })));
        }
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
        (simulation, [law, other])
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

/// What `component` answers to one step from `state`, outside the kernel.
fn proposal(
    component: &ThermalTimeline<'static>,
    state: &TimelineState<'static>,
    dt_s: f64,
    inputs: &[f64],
    gate: &CancelGate,
) -> Result<Proposal, TimelineRefusal> {
    let mut next = state.clone();
    let mut outputs = vec![f64::NAN; component.ports().len()];
    with_cx(gate, |cx| component.propose(state, dt_s, inputs, &mut next, &mut outputs, cx))
}

/// A system of the timeline alone: `rates_w` held on cell ports `cell10`
/// and `cell20`, and the environment port routed to a sink named `room`.
fn heated_alone(
    component: ThermalTimeline<'static>,
    start: TimelineState<'static>,
    rates_w: [f64; 2],
) -> (System, NodeId) {
    let mut builder = SystemBuilder::new();
    let cells = builder.add("cells", component, start).expect("cells");
    for (name, rate_w) in [("cell10", rates_w[0]), ("cell20", rates_w[1])] {
        builder
            .boundary(
                builder.port(cells, name).expect("port"),
                Box::new(move |_, _, out: &mut [f64]| out[0] = rate_w),
            )
            .expect("heated");
    }
    builder
        .heat_route(
            builder.port(cells, "environment").expect("port"),
            &[],
            Some(heat_sink("room", 293.15).expect("sink")),
            RouteTemperature::Fixed(293.15),
        )
        .expect("route");
    (builder.seal().expect("sealed"), cells)
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
fn published_temperatures_are_those_at_the_start_of_the_step() {
    let gate = CancelGate::new_clock_free();
    let simulation = Scene::default().timeline();
    // One heat port, on cell A; both temperatures published, B first.
    let component = ThermalTimeline::new(vec![A], false, config(0.05))
        .expect("component")
        .publishing(vec![B, A])
        .expect("signals");
    let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(names, ["cell10", "temperature20", "temperature10"]);
    let start = component.start(&simulation).expect("start");
    let stale = later(&component, &simulation, 1.0);
    // Whatever heat the step delivers, the signals are the start
    // temperatures; the heat port's output is the mean over the step.
    for rate_w in [0.0, 40.0, -250.0] {
        let (_, _, outputs) = conforms(&component, &start, &stale, 0.01, &[rate_w]);
        assert_eq!(&outputs[1..], [295.0, 310.0]);
        let through = LINK_W_K * (310.0 - 295.0);
        let after_a = 310.0 + 0.01 * (rate_w - through) / CAPACITY;
        assert!((outputs[0] - 0.5 * (310.0 + after_a)).abs() <= 1.0e-12, "{}", outputs[0]);
    }
    // A step of three owner substeps publishes the temperatures it started
    // from, and an environment port before the signals does not move them.
    let film = Boundary {
        htc_w_m2_k: 25.0,
        ..insulated()
    };
    let exposed = Scene {
        boundaries: [film, insulated()],
        ..Scene::default()
    };
    let fine = ThermalTimeline::new(vec![A], true, config(0.004))
        .expect("component")
        .publishing(vec![B, A])
        .expect("signals");
    let state = fine.start(&exposed.timeline()).expect("start");
    let (next, _, outputs) = propose(&fine, &state, 0.01, &[40.0, 293.15], &gate).expect("step");
    assert_eq!(&outputs[2..], [295.0, 310.0]);
    assert!(outputs[1] < 0.0, "the film takes heat: {} W", outputs[1]);
    assert!(temperature(&next, B) > 295.0);

    // After a step they are the temperatures that step ended on.
    let (next, _, first) = propose(&component, &start, 0.01, &[40.0], &gate).expect("step");
    let (_, _, second) = propose(&component, &next, 0.01, &[40.0], &gate).expect("step");
    assert_eq!(second[1], temperature(&next, B));
    assert_eq!(second[2], temperature(&next, A));
    assert!(second[2] > first[2] && second[1] > first[1], "both cells warmed");

    assert_eq!(
        ThermalTimeline::new(vec![A], false, config(0.05))
            .expect("component")
            .publishing(vec![B, B])
            .unwrap_err(),
        TimelineBuildError::Duplicate { cell: 20 }
    );
    assert_eq!(
        ThermalTimeline::new(vec![A], false, config(0.05))
            .expect("component")
            .publishing(vec![B])
            .expect("signals")
            .publishing(vec![A, B])
            .unwrap_err(),
        TimelineBuildError::Duplicate { cell: 20 }
    );
    let absent = ThermalTimeline::new(vec![A], false, config(0.05))
        .expect("component")
        .publishing(vec![CellId(77)])
        .expect("signals");
    assert_eq!(
        absent.start(&simulation).unwrap_err(),
        TimelineRefusal::InactiveCell { cell: 77 }
    );
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
    let (mut system, cells) = heated_alone(component, start, [0.0, 1.5]);
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
// The implicit driver.
// ---------------------------------------------------------------------------

impl ConductivityLaw for Sensible {
    fn conductivity_w_m_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let kelvin = self.temperature_k(specific_h_j_kg)?;
        Ok(1.0 + self.conductivity_slope_per_k * (kelvin - 300.0))
    }
}

impl ImplicitConductionLaw for Sensible {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    /// 250 K to 400 K.
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [self.datum_j_kg - 50.0 * CP, self.datum_j_kg + 100.0 * CP]
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Ok(ImplicitConductionState {
            temperature_k: self.temperature_k(specific_h_j_kg)?,
            temperature_slope_kg_k_j: 1.0 / CP,
            conductivity_w_m_k: self.conductivity_w_m_k(specific_h_j_kg)?,
            conductivity_slope_kg_w_m_k_j: self.conductivity_slope_per_k / CP,
        })
    }
}

impl Scene {
    /// The scene with the events it starts with applied, and its laws: what
    /// a conduction graph is prepared against.
    fn active(&self) -> (DepositionSimulation<'static>, [&'static Sensible; 2]) {
        let (simulation, laws) = self.build_laws();
        let (active, _) = simulation
            .advance_to(self.birth_s, config(0.05), || Ok(()))
            .expect("the scene's own events");
        (active, laws)
    }

    /// The conduction graph of the scene's one link. Its half distances give
    /// the link the scene's conductance where the conductivity is 1 W/(m K),
    /// at 300 K.
    fn graph(
        &self,
        active: &DepositionSimulation<'static>,
        laws: [&'static Sensible; 2],
    ) -> ImplicitConductionGraph<'static> {
        let half_m = LINK_AREA_M2 / (2.0 * self.link_w_k);
        let link = [ConductiveLink {
            link_id: 7,
            half_distances_m: [half_m; 2],
            contact_resistance_m2_k_w: 0.0,
        }];
        let conduction = ConductionObserver::new(
            active,
            &[
                (MaterialId(0), laws[0] as &dyn ConductivityLaw),
                (MaterialId(1), laws[1] as &dyn ConductivityLaw),
            ],
            if self.linked && !self.outside { &link } else { &[] },
            &[],
            || Ok(()),
        )
        .expect("conduction");
        let conduction = if self.linked && self.outside {
            conduction
                .with_links_carried_outside(active, &[7], || Ok(()))
                .expect("carried outside")
        } else {
            conduction
        };
        ImplicitConductionGraph::new(
            active,
            &conduction,
            &[
                (MaterialId(0), laws[0] as &dyn ImplicitConductionLaw),
                (MaterialId(1), laws[1] as &dyn ImplicitConductionLaw),
            ],
            ImplicitConductionLinearControls::default(),
            || Ok(()),
        )
        .expect("graph")
    }
}

fn implicit_controls(maximum_dt_s: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s,
        cell_tolerance_j: 1.0e-11,
        total_tolerance_j: 1.0e-9,
        ..ImplicitConductionControls::default()
    }
}

/// The owner's own implicit solve to `target_s`, with held rates.
fn owner_solve(
    graph: &ImplicitConductionGraph<'static>,
    from: &DepositionSimulation<'static>,
    held: &HeldSources,
    target_s: f64,
    controls: ImplicitConductionControls,
    scheme: ImplicitScheme,
) -> (DepositionSimulation<'static>, ImplicitConductionReport) {
    let mut workspace = ImplicitConductionWorkspace::new(graph.cell_ids().len()).expect("workspace");
    let advanced = match scheme {
        ImplicitScheme::BackwardEuler => {
            graph.advance_observed_to(from, held, target_s, controls, &mut workspace, || Ok(()))
        }
        ImplicitScheme::Sdirk2 => {
            graph.advance_sdirk2_observed_to(from, held, target_s, controls, &mut workspace, || Ok(()))
        }
    };
    let (next, _, report) = advanced.expect("the owner's solve");
    (next, report)
}

#[test]
fn the_implicit_timeline_conforms_and_its_step_is_the_owners_solve() {
    let gate = CancelGate::new_clock_free();
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.5);
    let (dt_s, rates) = (0.25, [2.0, -0.5]);
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);
    let delivered = dt_s * (rates[0] + rates[1]);
    let mut ends = Vec::new();
    for scheme in [ImplicitScheme::BackwardEuler, ImplicitScheme::Sdirk2] {
        let component =
            ThermalTimeline::implicit(vec![A, B], scene.graph(&active, laws), controls, scheme)
                .expect("component");
        let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
        assert_eq!(names, ["cell10", "cell20", "environment"]);
        let start = component.start(&active).expect("start");
        // The limit is the controls' longest interval, or the proposal.
        assert_eq!(component.limit_dt_s(&start, 10.0).expect("limit"), 0.5);
        assert_eq!(component.limit_dt_s(&start, 0.1).expect("limit"), 0.1);
        let stale = later(&component, &active, 1.0);
        let inputs = [rates[0], rates[1], 293.15];
        let (receipt, _, outputs) = conforms(&component, &start, &stale, dt_s, &inputs);
        let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");

        // The owner's own solve with the rates held: one interval, the same
        // cells and the same stored change.
        let twin = scene.graph(&active, laws);
        let (native, report) = owner_solve(&twin, &active, &held, dt_s, controls, scheme);
        assert_eq!((report.energy.substeps, report.rejected_intervals), (1, 0));
        same_cells(&next, &native);
        let stored = receipt.stored_change_j;
        assert!((stored - report.energy.stored_energy_change_j).abs() <= 1.0e-12 * delivered);
        assert!((stored - delivered).abs() <= 1.0e-12 * delivered, "{stored} vs {delivered}");
        assert!(receipt.numerical_loss_j.abs() <= 1.0e-12 * delivered);
        for (port, id) in [A, B].into_iter().enumerate() {
            let mean = 0.5 * (temperature(&start, id) + native.temperature_k(id).expect("T"));
            assert!((outputs[port] - mean).abs() <= 1.0e-11, "{} vs {mean}", outputs[port]);
        }
        // Nothing exchanges here: the port carries only the rounding
        // between the owner's external heat and the port heat it was handed.
        assert!(outputs[2].abs() <= 1.0e-13, "{} W", outputs[2]);
        // The receipt names the tolerance the solve was held to (the
        // conformance check has already held the residual to it).
        assert_eq!(receipt.tolerance, controls.cell_tolerance_j);
        ends.push([temperature(&next, A), temperature(&next, B)]);
    }
    // Backward Euler by hand for two equal cells on one link: the sum of the
    // temperatures takes the delivered heat, and the difference d obeys
    // C (d' - d) = dt (qa - qb) - 2 dt G d'.
    let sum = 310.0 + 295.0 + delivered / CAPACITY;
    let difference = (CAPACITY * 15.0 + dt_s * (rates[0] - rates[1]))
        / (CAPACITY + 2.0 * dt_s * LINK_W_K);
    let by_hand = [0.5 * (sum + difference), 0.5 * (sum - difference)];
    for cell in 0..2 {
        assert!((ends[0][cell] - by_hand[cell]).abs() <= 1.0e-9, "{:?} vs {by_hand:?}", ends[0]);
    }
    // The second scheme is another quadrature of the same step.
    assert!((ends[1][0] - ends[0][0]).abs() > 1.0e-6, "{ends:?}");
}

#[test]
fn an_implicit_timeline_books_its_owners_conduction_entropy() {
    // Cell 10 also loses heat to a 280 K film over the faces the link does
    // not cover; cell 20 is insulated.
    let gate = CancelGate::new_clock_free();
    let film = Boundary {
        ambient_k: 280.0,
        htc_w_m2_k: 5.0,
        ..insulated()
    };
    let scene = Scene {
        boundaries: [film, insulated()],
        ..Scene::default()
    };
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.5);
    let (dt_s, rates) = (0.25, [2.0, -0.5]);
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);
    let inputs = [rates[0], rates[1], 280.0];
    for scheme in [ImplicitScheme::BackwardEuler, ImplicitScheme::Sdirk2] {
        let component =
            ThermalTimeline::implicit(vec![A, B], scene.graph(&active, laws), controls, scheme)
                .expect("component");
        let start = component.start(&active).expect("start");
        let (next, receipt, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
        let twin = scene.graph(&active, laws);
        let (native, report) = owner_solve(&twin, &active, &held, dt_s, controls, scheme);
        same_cells(&next, &native);
        let (booked, owner) = (
            receipt.entropy_generated_j_per_k,
            report.entropy.conduction_generation_j_k,
        );
        assert!(owner > 0.0);
        assert!(
            (booked - owner).abs() <= 1.0e-12 * owner,
            "{booked} vs {owner}"
        );
        if matches!(scheme, ImplicitScheme::BackwardEuler) {
            // Backward Euler's own reading at the end temperatures: the link's
            // G (Ta - Tb)^2 / (Ta Tb) and the film's q (T0 - Ta) / (Ta T0).
            let (ta, tb) = (temperature(&next, A), temperature(&next, B));
            let film_j = dt_s * 5.0 * (SURFACE_M2 - LINK_AREA_M2) * (280.0 - ta);
            let by_hand = dt_s * LINK_W_K * (ta - tb).powi(2) / (ta * tb)
                + film_j * (280.0 - ta) / (ta * 280.0);
            assert!(
                (booked - by_hand).abs() <= 1.0e-12 * by_hand,
                "{booked} vs {by_hand}"
            );
        }
    }
    // The explicit rule keeps no entropy account and books zero.
    let simulation = scene.timeline();
    let explicit = ThermalTimeline::new(vec![A, B], true, config(0.05)).expect("component");
    let start = explicit.start(&simulation).expect("start");
    let (_, receipt, _) = propose(&explicit, &start, 0.01, &inputs, &gate).expect("step");
    assert_eq!(
        receipt.entropy_generated_j_per_k.to_bits(),
        0.0_f64.to_bits()
    );
}

/// Heat `q_w` moved from cell 10 into cell 20 at every query: a closed
/// passive exchange while cell 10 is the warmer, declared with `tolerance_w`.
struct Moved {
    calorics: [&'static dyn ThermalLaw; 2],
    q_w: f64,
    tolerance_w: f64,
}

impl ImplicitCoupledHeatSourceLaw for Moved {
    fn caloric_laws(&self) -> &[&dyn ThermalLaw] {
        &self.calorics
    }

    fn power_dependencies(&self) -> &[(usize, usize)] {
        &[]
    }

    fn state_at_specific_enthalpies_and_time(
        &self,
        _specific_h_j_kg: &[f64],
        _time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitCoupledHeatSourceState, DepositionError> {
        checkpoint()?;
        Ok(ImplicitCoupledHeatSourceState {
            powers_w: vec![-self.q_w, self.q_w],
            power_enthalpy_derivatives_w_kg_j: Vec::new(),
        })
    }

    fn passive_exchange_tolerance_w(&self) -> Option<f64> {
        Some(self.tolerance_w)
    }
}

#[test]
fn an_implicit_timeline_books_a_declared_exchanges_production_and_none_below_zero() {
    // 0.05 W moved inside the graph from cell 10 (310 K) to cell 20 (295 K),
    // declared closed, adds its production to the link's. 1e-3 W moved the
    // other way, inside a declared 1 W, leaves the owner's signed loss in its
    // audit and books the conduction alone.
    let gate = CancelGate::new_clock_free();
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.5);
    let (dt_s, rates) = (0.25, [2.0, -0.5]);
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);
    let inputs = [rates[0], rates[1], 293.15];
    let calorics = laws.map(|law| law as &'static dyn ThermalLaw);
    for (q_w, tolerance_w) in [(0.05, 0.0), (-1.0e-3, 1.0)] {
        let moved: &'static Moved = Box::leak(Box::new(Moved {
            calorics,
            q_w,
            tolerance_w,
        }));
        let graph = || {
            scene
                .graph(&active, laws)
                .with_coupled_heat_source(&active, &[A, B], moved, 0, || Ok(()))
                .expect("declared exchange")
        };
        for scheme in [ImplicitScheme::BackwardEuler, ImplicitScheme::Sdirk2] {
            let component = ThermalTimeline::implicit(vec![A, B], graph(), controls, scheme)
                .expect("component");
            let start = component.start(&active).expect("start");
            let (_, receipt, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
            let (_, report) = owner_solve(&graph(), &active, &held, dt_s, controls, scheme);
            let (conduction, exchange) = (
                report.entropy.conduction_generation_j_k,
                report.entropy.passive_exchange_generation_j_k,
            );
            let (booked, owner) = if q_w > 0.0 {
                assert!(exchange > 0.0, "{exchange}");
                (receipt.entropy_generated_j_per_k, conduction + exchange)
            } else {
                assert!(exchange < 0.0, "{exchange}");
                (receipt.entropy_generated_j_per_k, conduction)
            };
            assert!(
                (booked - owner).abs() <= 1.0e-12 * owner,
                "{booked} vs {owner}"
            );
        }
    }
}

#[test]
fn a_backward_euler_timeline_can_return_the_temperatures_at_the_end_of_the_step() {
    let gate = CancelGate::new_clock_free();
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.5);
    let euler = || {
        let graph = scene.graph(&active, laws);
        ThermalTimeline::implicit(vec![A, B], graph, controls, ImplicitScheme::BackwardEuler)
            .expect("component")
    };
    let (mean, ending) = (euler(), euler().returning_end_temperatures().expect("end"));
    let start = ending.start(&active).expect("start");
    let stale = later(&ending, &active, 1.0);
    let (dt_s, inputs) = (0.25, [2.0, -0.5, 293.15]);
    let (receipt, _, outputs) = conforms(&ending, &start, &stale, dt_s, &inputs);
    let (next, _, _) = propose(&ending, &start, dt_s, &inputs, &gate).expect("step");
    // It is the same step: the same state and the same account. Only what
    // the two cell ports return differs.
    let (same, same_receipt, means) = propose(&mean, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next, same);
    assert_eq!(receipt, same_receipt);
    assert_eq!(outputs[2].to_bits(), means[2].to_bits());
    for (port, id) in [A, B].into_iter().enumerate() {
        let (before, after) = (temperature(&start, id), temperature(&next, id));
        assert_eq!(outputs[port].to_bits(), after.to_bits());
        assert_eq!(means[port].to_bits(), (0.5 * (before + after)).to_bits());
        assert!((after - before).abs() > 0.01, "{before} to {after}");
    }
    // Only the backward-Euler solve is written at its end state alone.
    let explicit = ThermalTimeline::new(vec![A, B], false, config(0.05)).expect("component");
    assert_eq!(
        explicit.returning_end_temperatures().err(),
        Some(TimelineBuildError::EndTemperatures)
    );
    let staged =
        ThermalTimeline::implicit(vec![A, B], scene.graph(&active, laws), controls, ImplicitScheme::Sdirk2)
            .expect("component");
    assert_eq!(
        staged.returning_end_temperatures().err(),
        Some(TimelineBuildError::EndTemperatures)
    );
}

/// The two cells of `apart`, which are not in contact, in a backward-Euler
/// timeline, with the kernel's ideal conductance of `LINK_W_K` between
/// their two heat ports; the ports return the end temperatures if `end`.
fn joined_through_ports(
    apart: &Scene,
    controls: ImplicitConductionControls,
    end: bool,
) -> (System, NodeId) {
    let (alone, laws) = apart.active();
    let timeline = ThermalTimeline::implicit(
        vec![A, B],
        apart.graph(&alone, laws),
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let timeline = if end {
        timeline.returning_end_temperatures().expect("end")
    } else {
        timeline
    };
    let start = timeline.start(&alone).expect("start");
    let mut builder = SystemBuilder::new();
    let cells = builder.add("cells", timeline, start).expect("cells");
    let wall = builder
        .add_with_tangent("wall", Conductance::new(LINK_W_K).expect("wall"), ())
        .expect("wall");
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).expect("port");
    let (a, b) = (port(&builder, wall, "a"), port(&builder, wall, "b"));
    builder
        .common_effort(port(&builder, cells, "cell10"), &[a])
        .expect("hot side");
    builder
        .common_effort(port(&builder, cells, "cell20"), &[b])
        .expect("cold side");
    builder
        .heat_route(
            port(&builder, cells, "environment"),
            &[],
            Some(heat_sink("room", 293.15).expect("sink")),
            RouteTemperature::Fixed(293.15),
        )
        .expect("environment");
    let mut system = builder.seal().expect("sealed");
    // The wall's two temperature inputs are in the loop and need a start.
    assert!(system.guess(a, &[310.0]));
    assert!(system.guess(b, &[295.0]));
    (system, cells)
}

#[test]
fn cells_joined_through_end_temperature_ports_follow_the_owners_implicit_link() {
    // The reference: two cells joined by the owner's link of 0.2 W/K and
    // stepped by the owner's backward-Euler solve. The graph: the same two
    // cells with no link inside, and the kernel's ideal conductance of
    // 0.2 W/K between their ports. Backward Euler writes the link's heat at
    // the end temperatures, so with ports that return those the two are
    // one discrete rule, and they agree step by step to what the two
    // solves leave: the owner's cell tolerance on either side, and the
    // kernel's 1e-12 of the temperatures it solves for, carried through
    // the wall.
    let gate = CancelGate::new_clock_free();
    let controls = implicit_controls(0.5);
    let (dt_s, steps) = (0.5, 40_u32);
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let twin = scene.graph(&active, laws);
    let apart = Scene {
        linked: false,
        ..Scene::default()
    };
    let (mut ending, ending_cells) = joined_through_ports(&apart, controls, true);
    let (mut mean, mean_cells) = joined_through_ports(&apart, controls, false);
    let unheated = HeldSources::new([(A, 0.0), (B, 0.0)]);
    let mut native = active.clone();
    let (mut worst_ending, mut worst_mean) = (0.0_f64, 0.0_f64);
    with_cx(&gate, |cx| {
        for step in 0..steps {
            let target_s = native.time_s() + dt_s;
            let (next, report) = owner_solve(
                &twin,
                &native,
                &unheated,
                target_s,
                controls,
                ImplicitScheme::BackwardEuler,
            );
            assert_eq!((report.energy.substeps, report.rejected_intervals), (1, 0));
            native = next;
            assert_eq!(ending.step(dt_s, cx).expect("kernel step").dt_s, dt_s, "step {step}");
            assert_eq!(mean.step(dt_s, cx).expect("kernel step").dt_s, dt_s, "step {step}");
            let (ended, meant) = (
                ending.state::<TimelineState<'static>>(ending_cells).expect("cells"),
                mean.state::<TimelineState<'static>>(mean_cells).expect("cells"),
            );
            for id in [A, B] {
                let owner = native.temperature_k(id).expect("T");
                worst_ending = worst_ending.max((temperature(ended, id) - owner).abs());
                worst_mean = worst_mean.max((temperature(meant, id) - owner).abs());
            }
        }
    });
    // The cells came most of the way together: 15 K apart at the start.
    let apart_k = native.temperature_k(A).expect("T") - native.temperature_k(B).expect("T");
    assert!(apart_k > 1.0 && apart_k < 5.0, "{apart_k}");
    let a_step_k = (2.0 * controls.cell_tolerance_j + LINK_W_K * dt_s * 1.0e-12 * 310.0) / CAPACITY;
    eprintln!(
        "through end-temperature ports: {worst_ending:.3e} K from the owner's link (allowed \
         {:.3e}); through mean-temperature ports: {worst_mean:.3e} K",
        f64::from(steps) * a_step_k
    );
    assert!(worst_ending <= f64::from(steps) * a_step_k, "{worst_ending}");
    // The default ports return the mean of the two ends, which makes the
    // wall's heat the trapezoidal rule's: another rule, first order away.
    assert!(worst_mean > 1.0e-3, "{worst_mean}");
}

#[test]
fn a_contact_carried_through_the_ports_still_covers_the_cells_faces() {
    // Two cells in contact that also lose heat to the room through a film
    // on their exposed faces. The reference is the owner's backward-Euler
    // solve with the contact's conduction inside. The graph names the
    // contact as carried outside the timeline and puts the kernel's ideal
    // conductance of the same 0.2 W/K between the two cell ports. The
    // contact covers 1e-3 m2 of each cell's 1e-2 m2, so in both the films
    // act on 9e-3 m2, and the two agree step by step as far as two solves
    // of one discrete rule can. A third system, whose cells are not in
    // contact at all, has films on the whole faces and ends elsewhere:
    // that is what losing the covered area would look like.
    let gate = CancelGate::new_clock_free();
    let controls = implicit_controls(0.5);
    let (dt_s, steps) = (0.5, 40_u32);
    let film = Boundary {
        htc_w_m2_k: 5.5,
        ..insulated()
    };
    let scene = Scene {
        boundaries: [film; 2],
        ..Scene::default()
    };
    let (active, laws) = scene.active();
    let twin = scene.graph(&active, laws);
    let carried = Scene {
        outside: true,
        boundaries: [film; 2],
        ..Scene::default()
    };
    let apart = Scene {
        linked: false,
        boundaries: [film; 2],
        ..Scene::default()
    };
    let (mut covered, covered_cells) = joined_through_ports(&carried, controls, true);
    let (mut bare, bare_cells) = joined_through_ports(&apart, controls, true);
    let unheated = HeldSources::new([(A, 0.0), (B, 0.0)]);
    let mut native = active.clone();
    let (mut worst_covered, mut worst_bare, mut to_room_j) = (0.0_f64, 0.0_f64, 0.0_f64);
    with_cx(&gate, |cx| {
        for step in 0..steps {
            let target_s = native.time_s() + dt_s;
            let (next, report) = owner_solve(
                &twin,
                &native,
                &unheated,
                target_s,
                controls,
                ImplicitScheme::BackwardEuler,
            );
            assert_eq!((report.energy.substeps, report.rejected_intervals), (1, 0));
            to_room_j -= report.energy.external_heat_j;
            native = next;
            assert_eq!(covered.step(dt_s, cx).expect("kernel step").dt_s, dt_s, "step {step}");
            assert_eq!(bare.step(dt_s, cx).expect("kernel step").dt_s, dt_s, "step {step}");
            let (kept, lost) = (
                covered.state::<TimelineState<'static>>(covered_cells).expect("cells"),
                bare.state::<TimelineState<'static>>(bare_cells).expect("cells"),
            );
            for id in [A, B] {
                let owner = native.temperature_k(id).expect("T");
                worst_covered = worst_covered.max((temperature(kept, id) - owner).abs());
                worst_bare = worst_bare.max((temperature(lost, id) - owner).abs());
            }
        }
    });
    // The films took about a sixth of the 93.5 J the cells held above the
    // room (5 J/K each, 16.85 K and 1.85 K above it at the start): two
    // films of 0.0495 W/K over 20 s.
    assert!(to_room_j > 15.0 && to_room_j < 18.5, "{to_room_j} J");
    let a_step_k = (2.0 * controls.cell_tolerance_j + LINK_W_K * dt_s * 1.0e-12 * 310.0) / CAPACITY;
    eprintln!(
        "a contact carried through the ports: {worst_covered:.3e} K from the owner's link \
         (allowed {:.3e}), {to_room_j:.3} J to the room; cells not in contact: {worst_bare:.3e} K",
        f64::from(steps) * a_step_k
    );
    assert!(worst_covered <= f64::from(steps) * a_step_k, "{worst_covered}");
    assert!(worst_bare > 1.0e-2, "{worst_bare}");
    // The kernel's account of the same run: what left the cells is in the
    // room, as the owner's own audit has it.
    let ledger = covered.ledger();
    let room_j = ledger.sink_heat_j();
    assert!((room_j - to_room_j).abs() <= 1.0e-9 * to_room_j, "{room_j} vs {to_room_j}");
    assert!(ledger.energy_residual_j().abs() <= 1.0e-9 * to_room_j);
}

#[test]
fn an_implicit_step_is_not_bound_by_conduction() {
    let gate = CancelGate::new_clock_free();
    // A 400 W/K link: the explicit rule may not step past 10 ms.
    let stiff = Scene {
        link_w_k: 400.0,
        ..Scene::default()
    };
    let explicit = ThermalTimeline::new(vec![A, B], false, config(1.0)).expect("component");
    let bounded = explicit.start(&stiff.timeline()).expect("start");
    assert!(explicit.limit_dt_s(&bounded, 1.0).expect("limit") <= 0.0101);
    let (active, laws) = stiff.active();
    let controls = implicit_controls(0.5);
    let component = ThermalTimeline::implicit(
        vec![A, B],
        stiff.graph(&active, laws),
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let start = component.start(&active).expect("start");
    // The stiff link does not move the implicit limit.
    assert_eq!(component.limit_dt_s(&start, 10.0).expect("limit"), 0.5);
    let (dt_s, rates) = (0.5, [3.0, 1.0]);
    let inputs = [rates[0], rates[1], 293.15];
    let (next, receipt, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.simulation().time_s(), 0.5);
    // Fifty explicit bounds in one backward-Euler interval: the difference
    // d obeys C (d' - d) = dt (qa - qb) - 2 dt G d', so d' = 76 / 405 K.
    // Fifty explicit substeps would leave about 0.0025 K.
    let apart_k = temperature(&next, A) - temperature(&next, B);
    assert!((apart_k - 76.0 / 405.0).abs() <= 1.0e-9, "{apart_k}");
    let delivered = dt_s * (rates[0] + rates[1]);
    assert!((receipt.stored_change_j - delivered).abs() <= 1.0e-10);

}

#[test]
fn a_rejected_implicit_interval_is_the_kernels_to_shorten() {
    let gate = CancelGate::new_clock_free();
    // A conductivity that rises with temperature and a solve allowed two
    // Newton iterations: the owner, left to itself, rejects an interval and
    // halves it until two iterations are enough (measured here: 15 accepted
    // intervals and 118 rejected ones in the half second).
    let curved = Scene {
        link_w_k: 400.0,
        conductivity_slope_per_k: 0.01,
        ..Scene::default()
    };
    let (active, laws) = curved.active();
    let strict = ImplicitConductionControls {
        maximum_newton_iterations: 2,
        ..implicit_controls(0.5)
    };
    let (target_s, rates) = (0.5, [3.0, 1.0]);
    let held = HeldSources::new([(A, rates[0]), (B, rates[1])]);
    let twin = curved.graph(&active, laws);
    let (native, report) =
        owner_solve(&twin, &active, &held, target_s, strict, ImplicitScheme::BackwardEuler);
    assert!(report.rejected_intervals > 100, "{report:?}");
    assert!(report.energy.substeps > 10, "{report:?}");

    // The Component does not halve inside a step: it asks for half the
    // step, down to the controls' shortest interval, where the owner's
    // refusal comes through.
    let implicit = |controls| {
        ThermalTimeline::implicit(
            vec![A, B],
            curved.graph(&active, laws),
            controls,
            ImplicitScheme::BackwardEuler,
        )
        .expect("component")
    };
    let component = implicit(strict);
    let start = component.start(&active).expect("start");
    let inputs = [rates[0], rates[1], 293.15];
    assert_eq!(
        proposal(&component, &start, 0.5, &inputs, &gate).expect("answer"),
        Proposal::Shorten { dt_s: 0.25 }
    );
    let not_converged = |answer: Result<Proposal, TimelineRefusal>| {
        matches!(
            answer,
            Err(TimelineRefusal::Implicit(ImplicitConductionError::NotConverged { .. }))
        )
    };
    // A half of exactly the shortest interval is still asked for; under
    // it, the owner's refusal.
    let at_floor = implicit(ImplicitConductionControls {
        minimum_dt_s: 0.25,
        ..strict
    });
    assert_eq!(
        proposal(&at_floor, &start, 0.5, &inputs, &gate).expect("answer"),
        Proposal::Shorten { dt_s: 0.25 }
    );
    let under_floor = implicit(ImplicitConductionControls {
        minimum_dt_s: 0.3,
        ..strict
    });
    assert!(not_converged(proposal(&under_floor, &start, 0.5, &inputs, &gate)));
    // The owner's depth is kept. Allowed one halving of its half-second
    // interval, the owner alone refuses this half second, and so does the
    // Component once a quarter second has been rejected.
    let shallow = ImplicitConductionControls {
        maximum_halvings: 1,
        ..strict
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).expect("workspace");
    assert!(
        twin.advance_observed_to(&active, &held, target_s, shallow, &mut workspace, || Ok(()))
            .is_err()
    );
    let capped = implicit(shallow);
    assert_eq!(
        proposal(&capped, &start, 0.5, &inputs, &gate).expect("answer"),
        Proposal::Shorten { dt_s: 0.25 }
    );
    assert!(not_converged(proposal(&capped, &start, 0.25, &inputs, &gate)));

    // In the kernel the shortened steps are the owner's own intervals: as
    // many steps as it accepts, as many shortenings as it rejects, and the
    // same cells.
    let (mut system, cells) = heated_alone(component, start, rates);
    let (mut steps, mut shortenings) = (0_usize, 0_usize);
    with_cx(&gate, |cx| {
        while system.time_s() < target_s {
            let step = system.step(target_s - system.time_s(), cx).expect("step");
            steps += 1;
            shortenings += step.shortenings as usize;
            assert!(steps < 1000);
        }
    });
    assert_eq!(steps, report.energy.substeps);
    assert_eq!(shortenings, report.rejected_intervals);
    let reached = system.state::<TimelineState<'static>>(cells).expect("state");
    same_cells(reached, &native);
    let delivered = target_s * (rates[0] + rates[1]);
    let stored = system.ledger().stored_change_j();
    assert!((stored - delivered).abs() <= 1.0e-9, "{stored} vs {delivered}");
}

#[test]
fn a_loose_implicit_solve_reports_its_residual_and_books_what_it_does_not_store() {
    // Two Newton iterations are all this solve gets, and a cell tolerance
    // of 1e-4 J accepts what they leave (measured: 3e-5 J in a cell).
    let curved = Scene {
        link_w_k: 400.0,
        conductivity_slope_per_k: 0.01,
        ..Scene::default()
    };
    let (active, laws) = curved.active();
    let loose = ImplicitConductionControls {
        maximum_newton_iterations: 2,
        cell_tolerance_j: 1.0e-4,
        total_tolerance_j: 1.0e-3,
        ..ImplicitConductionControls::default()
    };
    let component = ThermalTimeline::implicit(
        vec![A, B],
        curved.graph(&active, laws),
        ImplicitConductionControls {
            maximum_dt_s: 0.5,
            ..loose
        },
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let start = component.start(&active).expect("start");
    let stale = later(&component, &active, 1.0);
    let (dt_s, inputs) = (0.5, [3.0, 1.0, 293.15]);
    let (receipt, _, outputs) = conforms(&component, &start, &stale, dt_s, &inputs);
    // The receipt says how well the solve closed each cell, and what it
    // was asked.
    assert!(receipt.residual > 1.0e-6, "{}", receipt.residual);
    assert!(receipt.residual <= receipt.tolerance);
    assert_eq!(receipt.tolerance, 1.0e-4);
    // After its Newton steps, what one cell's balance lacks here is heat
    // the other has too much of (the owner's form is conservative and
    // nothing in this scene exchanges), so the total is kept to rounding.
    let delivered = dt_s * (inputs[0] + inputs[1]);
    assert!((receipt.stored_change_j - delivered).abs() <= 1.0e-12 * delivered);
    assert!(receipt.numerical_loss_j.abs() <= 1.0e-12 * delivered, "{}", receipt.numerical_loss_j);
    // The environment wire is the owner's external heat less the port
    // heat: nothing, whatever the solve left.
    assert!(outputs[2].abs() <= 1.0e-13, "{} W", outputs[2]);

    // The owner accepts an interval without iterating when every cell's
    // balance at the start state is already inside its tolerance. On the
    // 0.2 W/K link over 20 us, cell A's 6e-5 J of port heat is what the
    // link would carry away and cell B lacks 8e-5 J: both under 1e-4 J.
    // The cells keep their enthalpy, and the receipt books every joule
    // delivered as numerical loss, where the ledger shows it.
    let slack = Scene::default();
    let (active, laws) = slack.active();
    let component = ThermalTimeline::implicit(
        vec![A, B],
        slack.graph(&active, laws),
        ImplicitConductionControls {
            maximum_dt_s: 0.5,
            ..loose
        },
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let start = component.start(&active).expect("start");
    let stale = later(&component, &active, 1.0);
    let brief_s = 2.0e-5;
    let (receipt, port_energy_in_j, _) = conforms(&component, &start, &stale, brief_s, &inputs);
    let unstored = brief_s * (inputs[0] + inputs[1]);
    assert!((port_energy_in_j - unstored).abs() <= 1.0e-12 * unstored);
    assert_eq!(receipt.stored_change_j, 0.0);
    assert!((receipt.numerical_loss_j - unstored).abs() <= 1.0e-12 * unstored);
    assert!(receipt.residual > 1.0e-5 && receipt.residual <= receipt.tolerance);
}

#[test]
fn an_implicit_run_ends_a_step_on_a_boundary_change_and_refuses_a_birth() {
    let gate = CancelGate::new_clock_free();
    let film = |htc_w_m2_k| Boundary {
        htc_w_m2_k,
        ..insulated()
    };
    // A film on cell A, and a fan change at 0.4 s.
    let scene = Scene {
        boundaries: [film(25.0), insulated()],
        extra: vec![boundary_event(0.4, A, film(5.0))],
        ..Scene::default()
    };
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.15);
    let component = ThermalTimeline::implicit(
        vec![A, B],
        scene.graph(&active, laws),
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let start = component.start(&active).expect("start");
    assert_eq!(component.next_event_s(&start), Some(0.4));
    let (mut system, cells) = heated_alone(component, start, [0.0, 1.5]);
    let (mut steps, mut landed) = (0_usize, false);
    with_cx(&gate, |cx| {
        while system.time_s() < 1.0 - 1.0e-12 {
            let report = system.step(1.0 - system.time_s(), cx).expect("step");
            landed |= (report.time_s - 0.4).abs() <= 1.0e-12;
            steps += 1;
            assert!(steps < 100);
        }
    });
    // 0.15, 0.30, then a step cut to the change at 0.40, then four more.
    assert!(landed, "a step ends on the fan change");
    assert_eq!(steps, 7);

    // The owner alone over the same second with the same held heat takes
    // the same intervals and reaches the same cells.
    let held = HeldSources::new([(A, 0.0), (B, 1.5)]);
    let twin = scene.graph(&active, laws);
    let (native, report) =
        owner_solve(&twin, &active, &held, 1.0, controls, ImplicitScheme::BackwardEuler);
    assert_eq!((report.energy.substeps, report.energy.events), (7, 1));
    let reached = system.state::<TimelineState<'static>>(cells).expect("state");
    assert!((reached.simulation().time_s() - 1.0).abs() <= 1.0e-12);
    assert_eq!(reached.simulation().cells()[&A].boundary, film(5.0));
    for id in [A, B] {
        let (kernel, owner) = (temperature(reached, id), native.temperature_k(id).expect("T"));
        assert!((kernel - owner).abs() <= 1.0e-10, "{kernel} vs {owner}");
    }
    // What the film took is in the room: the owner's external heat less
    // the port heat.
    let through_film = report.energy.external_heat_j - 1.5;
    assert!(through_film < -0.1, "the film cools cell A: {through_film} J");
    let room = system.ledger().sink_heat_j();
    assert!((room + through_film).abs() <= 1.0e-9 * through_film.abs(), "{room} vs {through_film}");

    // A birth changes the graph's bindings: the owner refuses the step
    // that reaches it.
    let growing = Scene {
        extra: vec![TimedEvent {
            time_s: 0.3,
            event: Event::Birth(Scene::default().cell(CellId(30), 1, 320.0, insulated())),
        }],
        ..Scene::default()
    };
    let (active, laws) = growing.active();
    let component = ThermalTimeline::implicit(
        vec![A, B],
        growing.graph(&active, laws),
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .expect("component");
    let start = component.start(&active).expect("start");
    assert_eq!(component.next_event_s(&start), Some(0.3));
    let inputs = [0.0, 0.0, 293.15];
    assert!(propose(&component, &start, 0.1, &inputs, &gate).is_ok());
    let refused = propose(&component, &start, 0.3, &inputs, &gate).unwrap_err();
    assert!(
        matches!(refused, TimelineRefusal::Implicit(ImplicitConductionError::Invalid(_))),
        "{refused}"
    );
}

/// A constant heat source bound to a cell of the conduction graph: heat the
/// graph owns, which no cell's boundary shows.
struct Steady {
    law: &'static Sensible,
    power_w: f64,
}

impl ImplicitHeatSourceLaw for Steady {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.law
    }

    fn state_at_specific_enthalpy(
        &self,
        _: f64,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        Ok(ImplicitHeatSourceState {
            power_w: self.power_w,
            power_enthalpy_derivative_w_kg_j: 0.0,
        })
    }
}

#[test]
fn heat_the_conduction_graph_owns_is_on_the_environment_wire() {
    let gate = CancelGate::new_clock_free();
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let source: &'static Steady = Box::leak(Box::new(Steady {
        law: laws[0],
        power_w: 3.0,
    }));
    let controls = implicit_controls(0.5);
    for scheme in [ImplicitScheme::BackwardEuler, ImplicitScheme::Sdirk2] {
        let graph = scene
            .graph(&active, laws)
            .with_heat_sources(&active, &[(A, source as &dyn ImplicitHeatSourceLaw)], || Ok(()))
            .expect("source");
        let component =
            ThermalTimeline::implicit(vec![B], graph, controls, scheme).expect("component");
        let start = component.start(&active).expect("start");
        let stale = later(&component, &active, 1.0);
        // Both cells are insulated: nothing about them says heat enters.
        assert!(start.simulation().cells().values().all(|cell| cell.boundary == insulated()));
        let (dt_s, rate_w) = (0.25, 0.75);
        let inputs = [rate_w, 293.15];
        let (receipt, port_energy_in_j, outputs) =
            conforms(&component, &start, &stale, dt_s, &inputs);
        // The graph's three watts are the environment port's output, and
        // with the port heat they are what the cells stored.
        assert!((outputs[1] - 3.0).abs() <= 1.0e-9, "{scheme:?}: {} W", outputs[1]);
        let entered = dt_s * (rate_w + 3.0);
        assert!((port_energy_in_j - entered).abs() <= 1.0e-9);
        assert!((receipt.stored_change_j - entered).abs() <= 1.0e-9);
        assert!(receipt.numerical_loss_j.abs() <= 1.0e-9, "{}", receipt.numerical_loss_j);
        let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
        let warmed_k = temperature(&next, A) + temperature(&next, B) - (310.0 + 295.0);
        assert!((CAPACITY * warmed_k - entered).abs() <= 1.0e-9, "{scheme:?}: {warmed_k} K");
    }
}

#[test]
fn the_implicit_timeline_refuses_what_it_cannot_carry() {
    let gate = CancelGate::new_clock_free();
    let scene = Scene::default();
    let (active, laws) = scene.active();
    let controls = implicit_controls(0.5);
    let implicit = |controls| {
        ThermalTimeline::implicit(
            vec![A, B],
            scene.graph(&active, laws),
            controls,
            ImplicitScheme::BackwardEuler,
        )
        .expect("component")
    };
    let component = implicit(controls);
    // The scene before its births are applied: the graph was not prepared
    // against that.
    assert_eq!(
        component.start(&scene.timeline()).unwrap_err(),
        TimelineRefusal::PendingEvent { time_s: 0.0 }
    );
    assert_eq!(
        ThermalTimeline::implicit(
            vec![A, A],
            scene.graph(&active, laws),
            controls,
            ImplicitScheme::Sdirk2
        )
        .unwrap_err(),
        TimelineBuildError::Duplicate { cell: 10 }
    );
    for maximum_dt_s in [0.0, -1.0, f64::INFINITY] {
        assert_eq!(
            ThermalTimeline::implicit(
                vec![A, B],
                scene.graph(&active, laws),
                ImplicitConductionControls {
                    maximum_dt_s,
                    ..controls
                },
                ImplicitScheme::BackwardEuler
            )
            .unwrap_err(),
            TimelineBuildError::Interval { maximum_dt_s }
        );
    }
    assert!(matches!(
        ThermalTimeline::implicit(
            vec![A, B],
            scene.graph(&active, laws),
            ImplicitConductionControls {
                maximum_dt_s: f64::NAN,
                ..controls
            },
            ImplicitScheme::BackwardEuler
        ),
        Err(TimelineBuildError::Interval { maximum_dt_s }) if maximum_dt_s.is_nan()
    ));
    let start = component.start(&active).expect("start");
    let inputs = [1.0, 0.0, 293.15];
    assert!(matches!(
        propose(&component, &start, 0.1, &[f64::NAN, 0.0, 293.15], &gate),
        Err(TimelineRefusal::Input { value }) if value.is_nan()
    ));
    // The owner's refusal of its controls comes through unchanged.
    let unworkable = implicit(ImplicitConductionControls {
        maximum_newton_iterations: 0,
        ..controls
    });
    assert!(matches!(
        propose(&unworkable, &start, 0.1, &inputs, &gate),
        Err(TimelineRefusal::Implicit(ImplicitConductionError::Invalid(_)))
    ));
    // So does the law's refusal, and a cancellation raised inside the
    // owner's call.
    laws[0].refuse.set(true);
    assert_eq!(
        propose(&component, &start, 0.1, &inputs, &gate).unwrap_err(),
        TimelineRefusal::Thermal(LAW_REFUSAL)
    );
    laws[0].refuse.set(false);
    let tripped: &'static CancelGate = Box::leak(Box::new(CancelGate::new_clock_free()));
    laws[0].tripwire.set(Some(tripped));
    assert_eq!(
        propose(&component, &start, 0.1, &inputs, tripped).unwrap_err(),
        TimelineRefusal::Cancelled
    );
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

/// The pairing's system: the network's first connector heats cell A by
/// `FRACTION_A`, cell B by `FRACTION_B` and a sink named `unallocated` with
/// the rest; the second heats cell B (a receiver two routes share). No body
/// is loaded. With `environment`, the timeline's environment port goes to a
/// sink named `room`. Returns the system, the network and the cells.
fn pairing(
    network: ViscousNetwork,
    bodies: Vec<RigidBodyState>,
    thermal: ThermalTimeline<'static>,
    cells: TimelineState<'static>,
    environment: bool,
) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, bodies).expect("network");
    let cells = builder.add("cells", thermal, cells).expect("cells");
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
    if environment {
        builder
            .heat_route(
                port(&builder, cells, "environment"),
                &[],
                Some(heat_sink("room", 293.15).expect("sink")),
                RouteTemperature::Fixed(293.15),
            )
            .expect("environment");
    }
    let mut system = builder.seal().expect("sealed");
    // The receivers' temperature comes back to the network's heat ports;
    // the kernel solves the two together and starts from this.
    for name in ["heat0", "heat1"] {
        let heat = system.port(moving, name).expect("heat port");
        assert!(system.guess(heat, &[300.0]));
    }
    (system, moving, cells)
}

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

    let (mut system, moving, cells) =
        pairing(network, start_bodies.clone(), thermal, start_cells.clone(), false);

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

#[test]
fn damper_heat_reaches_an_implicit_timeline_on_steps_the_explicit_rule_cannot_take() {
    let gate = CancelGate::new_clock_free();
    let (properties, links, start_bodies) = bodies();
    let scene = Scene {
        link_w_k: STIFF_LINK_W_K,
        ..Scene::default()
    };
    let (active, laws) = scene.active();
    // An eighth of a second: twelve and a half of the explicit rule's bounds.
    let solve = implicit_controls(0.125);
    assert!(solve.maximum_dt_s > 12.0 * explicit_bound_s(STIFF_LINK_W_K));
    let idle = [Wrench::ZERO; 3];
    for scheme in [ImplicitScheme::BackwardEuler, ImplicitScheme::Sdirk2] {
        let network =
            ViscousNetwork::new(properties.clone(), links.clone(), Gravity::ZERO, controls())
                .expect("network");
        let thermal =
            ThermalTimeline::implicit(vec![A, B], scene.graph(&active, laws), solve, scheme)
                .expect("timeline");
        let start_cells = thermal.start(&active).expect("start");
        let (mut system, moving, cells) =
            pairing(network, start_bodies.clone(), thermal, start_cells, true);

        // The same protocol by hand with the owners: mechanics over each
        // interval, then the owner's implicit solve fed the mean powers.
        let twin = scene.graph(&active, laws);
        let mut twin_bodies = start_bodies.clone();
        let mut twin_cells = active.clone();
        let steps = 12_usize;
        let (mut heat_0, mut heat_1) = (0.0_f64, 0.0_f64);
        with_cx(&gate, |cx| {
            for step in 0..steps {
                let dt_s = system.step(1.0, cx).expect("kernel step").dt_s;
                assert_eq!(dt_s, 0.125, "step {step}");
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
                let (next, report) =
                    owner_solve(&twin, &twin_cells, &held, target_s, solve, scheme);
                assert_eq!((report.energy.substeps, report.rejected_intervals), (1, 0));
                twin_cells = next;
                twin_bodies = motion.into_states_after();
            }
        });

        let reached_bodies = system.state::<Vec<RigidBodyState>>(moving).expect("bodies");
        assert_eq!(reached_bodies, &twin_bodies, "{scheme:?}");
        let reached_cells = system.state::<TimelineState<'static>>(cells).expect("cells");
        assert_eq!(reached_cells.simulation().time_s(), twin_cells.time_s());
        for id in [A, B] {
            let (kernel, hand) =
                (temperature(reached_cells, id), twin_cells.temperature_k(id).expect("T"));
            assert!((kernel - hand).abs() <= 1.0e-10, "{scheme:?} {id:?}: {kernel} vs {hand}");
        }
        // The stiff link has equalised the pair, and the damper heat warmed it.
        let apart_k = temperature(reached_cells, A) - temperature(reached_cells, B);
        assert!(apart_k.abs() < 0.05, "{scheme:?}: {apart_k} K apart");
        let total = heat_0 + heat_1;
        assert!(total > 0.5, "enough heat to measure: {total} J");
        let ledger = system.ledger();
        let account =
            |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
        let sink = |name: &str| {
            let found = ledger.sinks.iter().find(|sink| sink.name == name).expect("sink");
            found.heat_j.value()
        };
        let unallocated = (1.0 - FRACTION_A - FRACTION_B) * heat_0;
        assert!((sink("unallocated") - unallocated).abs() <= 1.0e-11 * total);
        // Nothing exchanges with the room.
        assert!(sink("room").abs() <= 1.0e-13, "{scheme:?}: {} J", sink("room"));
        // The cells stored what they were given: the solves closed to
        // rounding in the total (their per-cell tolerance is 1e-11 J).
        let closure_j = 1.0e-11 * total;
        let warmed = account("cells").stored_change_j.value();
        assert!((warmed - (total - unallocated)).abs() <= closure_j, "{scheme:?}: {warmed}");
        assert!((account("network").dissipated_j.value() - total).abs() <= 1.0e-11 * total);
        assert!(ledger.energy_residual_j().abs() <= closure_j);
        assert_eq!(ledger.steps, steps as u64);
    }
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

// Durable snapshots of the timeline (bead rc-1109 .5.4). The three kernel
// helpers are those of tests/snapshots.rs.
use fs_blake3::identity::CanonicalLimits;
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotAlgorithmIdV2,
    SnapshotBudgetStateIdV2, SnapshotDeterminismV2, SnapshotExecutionFingerprintIdV2,
    SnapshotLimitsV2, SnapshotPauseRequestIdV2, SnapshotProvenanceIdV2, SnapshotRngCounterIdV2,
};
use fs_exec::{DrainTracker, RunId};
use fs_system::SystemImage;
use fs_thermal::enthalpy::{ContactPatch, ContactPatchUpdate, ThermalStorageBasis};

fn limits() -> SnapshotLimitsV2 {
    SnapshotLimitsV2::new(
        1 << 22,
        4_096,
        CanonicalLimits::new(16_384, 4_096, 32, 32, 64),
        1 << 16,
        1 << 22,
        4_096,
    )
}

fn paused() -> PausedSnapshotBoundaryV2 {
    let gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(RunId(17), &gate);
    let worker = tracker.register_worker().expect("worker registers");
    gate.request();
    drop(worker);
    let report = tracker.finalize().expect("run drained");
    PausedSnapshotBoundaryV2::from_drain_report(
        report,
        SnapshotPauseRequestIdV2::from_bytes([0x5A; 32]),
        1,
    )
}

fn context(system: &System) -> ExpectedResumeContextV2 {
    ExpectedResumeContextV2::for_paused_state::<SystemImage>(
        SnapshotAlgorithmIdV2::from_bytes([0xA1; 32]),
        1,
        system.problem_id(),
        SnapshotRngCounterIdV2::from_bytes([0; 32]),
        SnapshotDeterminismV2::Deterministic,
        SnapshotExecutionFingerprintIdV2::from_bytes([0xE1; 32]),
        SnapshotBudgetStateIdV2::from_bytes([0xB1; 32]),
        SnapshotProvenanceIdV2::from_bytes([0xC1; 32]),
        paused(),
    )
}

/// `Sensible`'s temperature with no datum, under a name, a storage basis and
/// a slope bound of its own. `Sensible` has no caloric identity.
struct Named {
    name: Option<u8>,
    basis: ThermalStorageBasis,
    slope: f64,
}

impl ThermalLaw for Named {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }

    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.name.map(|byte| [byte; 32])
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        Ok(300.0 + specific_h_j_kg / CP)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.slope
    }
}

/// The first law's name, basis and slope bound, and the budget's links.
type Kept = (Option<u8>, ThermalStorageBasis, f64, usize);
const KEPT: Kept = (Some(1), ThermalStorageBasis::Enthalpy, 1.0 / CP, 8);

/// The timeline of `scene` before its start, on two named laws: its two
/// cells born at 0 s and joined by link 7, then the scene's later events.
fn named_timeline(scene: &Scene, first: Kept) -> DepositionSimulation<'static> {
    let (name, basis, slope, max_links) = first;
    let second = Named {
        name: Some(2),
        basis: ThermalStorageBasis::Enthalpy,
        slope: 1.0 / CP,
    };
    let [law, other] = [Named { name, basis, slope }, second].map(|law| &*Box::leak(Box::new(law)));
    let at = |event| TimedEvent { time_s: 0.0, event };
    let mut schedule = vec![
        at(Event::Birth(scene.cell(A, 0, 310.0, scene.boundaries[0]))),
        at(Event::Birth(scene.cell(B, 1, 295.0, scene.boundaries[1]))),
        at(Event::Contact(ThermalLink {
            id: 7,
            a: A,
            b: B,
            area_m2: LINK_AREA_M2,
            conductance_w_k: scene.link_w_k,
        })),
    ];
    schedule.extend(scene.extra.iter().cloned());
    let budget = DepositionBudget {
        max_cells: 8,
        max_links,
        max_events: 16,
        max_substeps: 100_000,
    };
    let laws = [
        (MaterialId(0), law as &dyn ThermalLaw),
        (MaterialId(1), other as &dyn ThermalLaw),
    ];
    DepositionSimulation::new(schedule, &laws, budget).expect("timeline")
}

fn separation(time_s: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Separation(7),
    }
}

fn state_words(simulation: &DepositionSimulation<'_>) -> Vec<u64> {
    let mut words = Vec::new();
    simulation.state_words(&mut words);
    words
}

/// Claims, declared before the first run: the owner's identity words are
/// equal for a timeline built again, and different for another law name,
/// storage basis, slope bound or budget limit, another field of a born cell
/// or of a link, and another time or kind of a later event. A timeline with
/// a law that has no caloric identity has no words and writes none.
#[test]
fn the_owners_identity_words_name_its_laws_its_budget_and_every_event() {
    let timeline = |change: usize| {
        let mut scene = Scene {
            extra: vec![separation(2.0)],
            ..Scene::default()
        };
        let mut first = KEPT;
        match change {
            0 | 1 => {}
            2 => first.0 = Some(3),
            3 => first.1 = ThermalStorageBasis::InternalEnergy,
            4 => first.2 = 2.0 / CP,
            5 => first.3 = 9,
            6 => scene.boundaries[0].source_w = 0.5,
            7 => scene.link_w_k = 0.25,
            8 => scene.extra[0].time_s = 2.5,
            9 => scene.extra[0].event = Event::Boundary(A, insulated()),
            _ => first.0 = None,
        }
        named_timeline(&scene, first)
    };
    let mut seen: Vec<Vec<u64>> = Vec::new();
    for change in 1..10 {
        let mut words = Vec::new();
        assert!(
            timeline(change).identity_words(&mut words),
            "change {change}"
        );
        assert!(
            !seen.contains(&words),
            "change {change} names an earlier timeline"
        );
        seen.push(words);
    }
    let mut again = Vec::new();
    assert!(timeline(0).identity_words(&mut again));
    assert_eq!(again, seen[0]);
    let mut none = vec![7];
    assert!(!timeline(10).identity_words(&mut none));
    assert_eq!(none, [7]);
}

/// Claims, declared before the first run. The owner's state words, read back
/// on a simulation that is still before its start (no cell, no link), give
/// the state that wrote them: its time, cells, links and the areas a contact
/// covers on each side, which the state's equality does not see; two states
/// that differ only in those areas have different words. A state past a
/// separation and a boundary change is read back without the link and with
/// the new boundary, and goes on to the bits of the run that was never
/// saved. Words the owner never writes, and words another timeline's budget
/// does not admit, are refused.
#[test]
fn the_owners_state_words_are_read_back_on_its_timeline_and_no_others_are() {
    let heated = Boundary {
        source_w: 0.25,
        ..insulated()
    };
    let scene = Scene {
        extra: vec![separation(2.0), boundary_event(3.0, B, heated)],
        ..Scene::default()
    };
    let kept = named_timeline(&scene, KEPT);
    let to = |from: &DepositionSimulation<'static>, time_s: f64| {
        from.advance_to(time_s, config(0.5), || Ok(()))
            .expect("advance")
            .0
    };
    let started = to(&kept, 0.0);
    let even = started.contact_patch(7).expect("link");
    let uneven = ContactPatch {
        covered_area_m2: [2.0e-3, 1.5e-3],
        ..even
    };
    let update = ContactPatchUpdate {
        id: 7,
        expected: Some(even),
        replacement: Some(uneven),
    };
    let covered = started
        .replace_contact_patches(&[update], || Ok(()))
        .expect("patch");
    let aged = to(&covered, 1.0);
    let words = state_words(&aged);
    assert_ne!(words, state_words(&to(&started, 1.0)), "the covered areas");

    let read = kept.restored(&words).expect("read back");
    assert!(read.shares_timeline(&kept));
    assert!(!read.shares_timeline(&named_timeline(&scene, KEPT)));
    assert_eq!(state_words(&read), words);
    assert_eq!(
        (read.time_s(), read.cells(), read.links()),
        (1.0, aged.cells(), aged.links())
    );
    assert_eq!(
        (read.contact_patch(7), read.next_event_time_s()),
        (Some(uneven), Some(2.0))
    );
    // A simulation before its start is a state the owner holds.
    assert!(kept.restored(&state_words(&kept)).is_some());

    let later = to(&aged, 3.5);
    let resumed = kept.restored(&state_words(&later)).expect("resumed");
    assert_eq!((resumed.cells().len(), resumed.links().len()), (2, 0));
    assert_eq!(resumed.cells()[&B].boundary, heated);
    assert_eq!(resumed.next_event_time_s(), None);
    // Read back on a simulation that still holds the link: nothing of it is kept.
    let on_started = started
        .restored(&state_words(&later))
        .expect("on a started one");
    assert_eq!(state_words(&on_started), state_words(&later));
    assert_eq!(
        state_words(&to(&resumed, 5.0)),
        state_words(&to(&later, 5.0))
    );

    // Two cells of 21 words after three, then one link of seven after its count.
    assert_eq!(words.len(), 53);
    let refused = |what: &str, change: &dyn Fn(&mut Vec<u64>)| {
        let mut changed = words.clone();
        change(&mut changed);
        assert!(kept.restored(&changed).is_none(), "{what}");
    };
    let (nan, zero, one) = (f64::NAN.to_bits(), 0.0_f64.to_bits(), 1.0_f64.to_bits());
    refused("a word short", &|words| words.truncate(52));
    refused("a word over", &|words| words.push(0));
    refused("a time that is not one", &|words| words[0] = nan);
    refused("a cursor past the schedule", &|words| words[1] = 6);
    refused("a cursor past a later event", &|words| words[1] = 4);
    refused("a count of other cells", &|words| words[2] = 1);
    refused("an unknown material", &|words| words[6] = 9);
    refused("a cell with no mass", &|words| words[7] = zero);
    refused("a reservoir tag never written", &|words| words[19] = 2);
    refused("cells out of order", &|words| {
        let (first, second) = words[3..45].split_at_mut(21);
        first.swap_with_slice(second);
    });
    refused("a link to a missing cell", &|words| words[47] = 99);
    refused("a covered area of nothing", &|words| words[51] = zero);
    refused("contacts beyond a surface", &|words| words[51] = one);
    let no_links = named_timeline(&scene, (KEPT.0, KEPT.1, KEPT.2, 0));
    assert!(no_links.restored(&words).is_none(), "a link over budget");
}

/// A timeline with a heat port on each of `cells`, stepped at most 0.05 s.
fn ported(cells: Vec<CellId>) -> ThermalTimeline<'static> {
    ThermalTimeline::new(cells, false, config(0.05)).expect("component")
}

/// The timeline on a port of cell A that keeps `kept`, and the state one
/// heated step after its start.
fn stepped(
    kept: &DepositionSimulation<'static>,
) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let timeline = ported(vec![A]).restoring(kept);
    let start = timeline.start(kept).expect("start");
    let gate = CancelGate::new_clock_free();
    let (next, ..) = propose(&timeline, &start, 0.05, &[3.0], &gate).expect("step");
    (timeline, next)
}

/// What `component` writes for `state`; a refusal writes nothing.
fn encoded(
    component: &ThermalTimeline<'static>,
    state: &TimelineState<'static>,
) -> Option<Vec<u64>> {
    let mut words = vec![9];
    let wrote = component.encode_state(state, &mut words);
    assert!(wrote || words == [9], "a refusal writes nothing");
    wrote.then(|| words.split_off(1))
}

fn identity(component: &ThermalTimeline<'static>) -> Vec<u64> {
    let mut words = Vec::new();
    component.identity_words(&mut words);
    words
}

/// Claims, declared before the first run. A timeline that keeps a named
/// timeline writes a state on it and reads it back. Its identity is equal
/// wherever the kept simulation stands, and different for other ports,
/// published cells, step limit, environment port and kept timeline. A state
/// on another timeline is not written. A state with an event due at its
/// time, without one of the port cells, or with a boundary that exchanges
/// heat and no environment port, is not read.
#[test]
fn the_timeline_writes_and_reads_a_state_on_the_named_timeline_it_keeps() {
    let scene = Scene {
        extra: vec![separation(0.2)],
        ..Scene::default()
    };
    let kept = named_timeline(&scene, KEPT);
    let (timeline, next) = stepped(&kept);
    let words = encoded(&timeline, &next).expect("a codec");
    assert_eq!(words, state_words(next.simulation()));
    let read = timeline.decode_state(&words).expect("read back");
    assert!(read == next && encoded(&timeline, &read).as_ref() == Some(&words));

    let base = identity(&timeline);
    assert!(!base.is_empty());
    let moved = ported(vec![A]).restoring(next.simulation());
    assert_eq!(identity(&moved), base);
    let other = Scene {
        link_w_k: 0.25,
        extra: vec![separation(0.2)],
        ..Scene::default()
    };
    let other = named_timeline(&other, KEPT);
    let limited = ThermalTimeline::new(vec![A], false, config(0.04)).expect("component");
    let open = ThermalTimeline::new(vec![A], true, config(0.05)).expect("component");
    let published = ported(vec![A]).publishing(vec![B]).expect("published");
    for (what, changed) in [
        ("another port", ported(vec![B]).restoring(&kept)),
        ("a second port", ported(vec![A, B]).restoring(&kept)),
        ("a step limit", limited.restoring(&kept)),
        ("an environment port", open.clone().restoring(&kept)),
        ("another timeline", ported(vec![A]).restoring(&other)),
        ("a published cell", published.restoring(&kept)),
    ] {
        let named = identity(&changed);
        assert!(!named.is_empty() && named != base, "{what}");
    }

    // A state on another timeline, built alike, is not this timeline's.
    let alike = named_timeline(&scene, KEPT);
    let foreign = ported(vec![A]).start(&alike).expect("start");
    assert_eq!(encoded(&timeline, &foreign), None);
    // The same state at the time of the separation still ahead of it: due.
    let mut due = words.clone();
    due[0] = 0.2_f64.to_bits();
    assert!(timeline.decode_state(&due).is_none());
    let elsewhere = ported(vec![CellId(30)]).restoring(&kept);
    assert!(elsewhere.decode_state(&words).is_none());
    let mut exchanging = Scene::default();
    exchanging.boundaries[0].source_w = 0.25;
    let warm = named_timeline(&exchanging, KEPT);
    let (born, _) = warm.advance_to(0.0, config(0.05), || Ok(())).expect("born");
    let warm_words = state_words(&born);
    let closed = ported(vec![A]).restoring(&warm);
    assert!(closed.decode_state(&warm_words).is_none());
    assert!(open.restoring(&warm).decode_state(&warm_words).is_some());
}

/// Claims, declared before the first run: a timeline that keeps no
/// timeline, one that keeps a timeline with a law that has no name, and an
/// implicit one that keeps a named timeline have no identity words, write
/// no state and read none.
#[test]
fn a_timeline_that_keeps_none_or_an_unnamed_law_or_is_implicit_has_no_codec() {
    let scene = Scene::default();
    let kept = named_timeline(&scene, KEPT);
    let (timeline, next) = stepped(&kept);
    let words = encoded(&timeline, &next).expect("a codec");
    let unnamed = named_timeline(&scene, (None, KEPT.1, KEPT.2, KEPT.3));
    let (active, laws) = scene.active();
    let graph = scene.graph(&active, laws);
    let scheme = ImplicitScheme::BackwardEuler;
    let implicit = ThermalTimeline::implicit(vec![A, B], graph, implicit_controls(0.5), scheme);
    for (what, none) in [
        ("no kept timeline", ported(vec![A])),
        ("a law with no name", ported(vec![A]).restoring(&unnamed)),
        (
            "an implicit driver",
            implicit.expect("implicit").restoring(&kept),
        ),
    ] {
        assert!(identity(&none).is_empty(), "{what}");
        assert_eq!(encoded(&none, &next), None, "{what}");
        assert!(none.decode_state(&words).is_none(), "{what}");
    }
}

/// Claims, declared before the first run (condition (a) of 35536). The kept
/// simulation is the one before its start, with no cell and no link. A
/// system stepped past the separation is sealed; a second system, built
/// alike and still holding the start's two cells and one link, resumes it
/// and then holds the two cells and no link, and both go on to the same
/// bits. The same system on a timeline that keeps none is refused a
/// snapshot by its node's name.
#[test]
fn a_snapshot_past_a_separation_resumes_without_the_link_and_goes_on_alike() {
    let gate = CancelGate::new_clock_free();
    let scene = Scene {
        extra: vec![separation(0.2)],
        ..Scene::default()
    };
    let kept = named_timeline(&scene, KEPT);
    let build = |restoring: bool| {
        let component = ThermalTimeline::new(vec![A, B], true, config(0.05)).expect("component");
        let component = if restoring {
            component.restoring(&kept)
        } else {
            component
        };
        let start = component.start(&kept).expect("start");
        heated_alone(component, start, [2.0, -1.0])
    };
    let (mut first, cells) = build(true);
    let held = |system: &System| {
        let state = system
            .state::<TimelineState<'static>>(cells)
            .expect("state");
        let simulation = state.simulation();
        (
            simulation.cells().len(),
            simulation.links().len(),
            state_words(simulation),
        )
    };
    let steps = |system: &mut System, count: usize| {
        with_cx(&gate, |cx| {
            for _ in 0..count {
                system.step(0.05, cx).expect("step");
            }
        });
    };
    steps(&mut first, 6);
    let sealed = first
        .snapshot(&context(&first), limits(), || false)
        .expect("sealed");
    let opened =
        SystemImage::unseal_v2_expected(sealed.bytes(), &sealed.expectation(), limits(), || false)
            .expect("opened");
    let (mut second, _) = build(true);
    assert_eq!(held(&second).1, 1, "the start's link");
    second.resume(&opened).expect("resumed");
    let resumed = held(&second);
    assert_eq!((resumed.0, resumed.1), (2, 0));
    assert_eq!(resumed, held(&first));
    steps(&mut first, 4);
    steps(&mut second, 4);
    assert_eq!(held(&second), held(&first));
    assert_eq!(second.time_s().to_bits(), first.time_s().to_bits());

    let (unnamed, _) = build(false);
    let error = unnamed
        .snapshot(&context(&unnamed), limits(), || false)
        .expect_err("a timeline with no codec");
    assert_eq!(format!("{error:?}"), r#"NoCodec { node: "cells" }"#);
}
