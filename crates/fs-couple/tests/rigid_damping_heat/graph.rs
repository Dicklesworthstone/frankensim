//! The declared graph against the hand-written coupling.
//!
//! `advance_candidate` is this module's explicit damping-heat step: the
//! thermal owner picks each interval, the damped network advances once over
//! it, and connector heat enters the authored cells as mean power. Here the
//! same two owners are kernel Components
//! (`fs_mbd::component::network::ViscousNetwork` and
//! `fs_thermal::component::ThermalTimeline`), the authored heat fractions are
//! heat routes, and the kernel in `fs-system` steps them. No code is written
//! for the combination. The same is done for the implicit thermal entry
//! points and for temperature-dependent damping.
//!
//! Where the module succeeds the two must agree: the bodies bit for bit (to
//! the kernel's solve tolerance when a temperature signal feeds the
//! damping), the cells to the rounding of handing heat over as energy
//! instead of as a rate, and the report's accounts with the kernel's ledger.
//! When the implicit solve rejects an interval, the module re-runs its
//! mechanics over half of it; the graph's timeline answers `Shorten` and the
//! kernel re-steps every component over half of it, so the two take the same
//! intervals, provided the timeline's solve is tried with the heat it will
//! be given. That holds when the network does not depend on what the cells
//! return (heat released at a fixed temperature, as in the module). When the
//! receivers' temperature is fed back, the kernel solves the two together
//! and the timeline can answer on an iterate: the graph's run can then be
//! finer than the module's, and in the fixture below it is.
//!
//! Limits of the comparison. The module advances its mechanics over the
//! thermal clock's `end - start`, the kernel over the step it negotiated.
//! On every step of these fixtures the two are the same number; that is how
//! these times round, not a property of every grid. Where they differ by an
//! ulp of the time, the bodies agree to rounding and not by bits. Not
//! compared: the
//! module refuses a whole advance when its law has no form at the final
//! temperature, where the graph takes that step and refuses the next. The
//! `relaxation` entry points are compared in `relaxation/graph.rs` and the
//! `driven` ones in `driven_graph.rs`.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::{DampingLaw, ViscousNetwork};
use fs_system::{Ledger, NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};

pub(crate) fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
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

/// The kernel keeps component states as `'static` values, so the laws the
/// timeline borrows live as long as the test process.
fn static_laws() -> &'static [Caloric; 2] {
    Box::leak(Box::new(laws()))
}

/// The graph of the pairing: every connector's heat routed to the two cells
/// by its authored fractions, the rest to a sink named `unallocated`; the
/// timeline's own boundary heat to a sink named `environment`. Returns the
/// system and its two nodes, the network and the cells.
fn graph(
    bodies: &Bodies,
    routes: &[BTreeMap<CellId, f64>; 2],
    heat: &DepositionSimulation<'static>,
) -> (System, NodeId, NodeId) {
    let network = ViscousNetwork::new(
        bodies.properties.to_vec(),
        bodies.links.to_vec(),
        bodies.gravity,
        config().mechanical,
    )
    .unwrap();
    // The fixtures' cells sit on a bed (of zero conductance), so the
    // timeline declares the port its boundary heat leaves by.
    let timeline = ThermalTimeline::new(vec![CellId(0), CellId(1)], true, config().thermal).unwrap();
    let start = timeline.start(heat).unwrap();
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, bodies.states.to_vec()).unwrap();
    let cells = builder.add("cells", timeline, start).unwrap();
    // These fixtures hold no external load on any body.
    for body in 0..bodies.loads.len() {
        let force = builder.port(moving, &format!("force{body}")).unwrap();
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    for (connector, fractions) in routes.iter().enumerate() {
        let receivers: Vec<_> = fractions
            .iter()
            .map(|(cell, fraction)| {
                (builder.port(cells, &format!("cell{}", cell.0)).unwrap(), *fraction)
            })
            .collect();
        builder
            .heat_route(
                builder.port(moving, &format!("heat{connector}")).unwrap(),
                &receivers,
                Some(heat_sink("unallocated", 300.0).unwrap()),
                RouteTemperature::SinkEquivalent,
            )
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    for connector in 0..routes.len() {
        let heat_port = system.port(moving, &format!("heat{connector}")).unwrap();
        assert!(system.guess(heat_port, &[300.0]));
    }
    (system, moving, cells)
}

pub(crate) fn sink(ledger: &Ledger, name: &str) -> f64 {
    ledger
        .sinks
        .iter()
        .find(|sink| sink.name == name)
        .map_or(0.0, |sink| sink.heat_j.value())
}

pub(crate) fn node<'a>(ledger: &'a Ledger, name: &str) -> &'a fs_system::ledger::NodeAccount {
    ledger.nodes.iter().find(|node| node.name == name).unwrap()
}

/// Step `system` to `target_s`, asking each time for all that is left: the
/// thermal Component's bound and events cut the step, as the thermal owner's
/// do inside `advance_candidate`.
fn run(system: &mut System, target_s: f64) -> Vec<f64> {
    run_counting(system, target_s).0
}

/// [`run`], and the number of times the kernel shortened a step on the way.
pub(crate) fn run_counting(system: &mut System, target_s: f64) -> (Vec<f64>, usize) {
    let gate = CancelGate::new_clock_free();
    let mut intervals = Vec::new();
    let mut shortenings = 0_usize;
    with_cx(&gate, |cx| {
        while system.time_s() < target_s {
            let report = system.step(target_s - system.time_s(), cx).unwrap();
            intervals.push(report.dt_s);
            shortenings += report.shortenings as usize;
            assert!(intervals.len() < 1000);
        }
    });
    (intervals, shortenings)
}

#[test]
fn the_declared_graph_reproduces_the_explicit_candidate() {
    let laws = static_laws();
    let bodies = Bodies::new([0.7, 1.1]);
    let routes = routes();
    // Without and with the fixtures' thermal events: an authored 2 W source
    // on cell 0 from the start, a contact, and the source raised to 4 W at
    // 0.35 s. Targets on and off the 0.125 s thermal step.
    for (events, target_s) in [(false, 0.5), (false, 0.45), (true, 0.5), (true, 0.45)] {
        let heat = thermal(laws, events, true);
        let candidate = advance_candidate(
            bodies.input(&routes),
            &heat,
            &(),
            target_s,
            config(),
            || false,
        )
        .unwrap();
        let report = candidate.report();

        let (mut system, moving, cells) = graph(&bodies, &routes, &heat);
        let intervals = run(&mut system, target_s);
        let case = format!("events {events}, target {target_s}");

        // The same intervals, so the same mechanics: bit for bit.
        assert_eq!(report.thermal.substeps, intervals.len(), "{case}: {intervals:?}");
        let reached = system.state::<Vec<RigidBodyState>>(moving).unwrap();
        assert_eq!(reached.as_slice(), candidate.states(), "{case}");
        // The same cells, to the rounding of the heat's hand-over.
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        assert_eq!(
            timeline.simulation().time_s().to_bits(),
            candidate.thermal().time_s().to_bits(),
            "{case}: the two thermal clocks"
        );
        for id in [CellId(0), CellId(1)] {
            let kernel = timeline.simulation().temperature_k(id).unwrap();
            let module = candidate.thermal().temperature_k(id).unwrap();
            assert!((kernel - module).abs() <= 1e-11, "{case}: {id:?} {kernel} vs {module}");
            assert!(kernel > 300.0 + 1e-3, "{case}: the damper heat warmed {id:?}");
        }

        // The module's report and the kernel's ledger are the same account.
        let ledger = system.ledger();
        let dissipated = report.dissipated_energy_j;
        assert!(dissipated > 0.5, "{case}: enough heat to compare, {dissipated} J");
        let close = |ours: f64, theirs: f64, what: &str| {
            assert!(
                (ours - theirs).abs() <= 1e-11 * dissipated,
                "{case}: {what}: ledger {ours:e}, report {theirs:e}"
            );
        };
        close(node(ledger, "network").dissipated_j.value(), dissipated, "dissipated");
        close(
            report.connector_heat_j.iter().sum::<f64>(),
            dissipated,
            "the report's connectors sum to its dissipation",
        );
        close(sink(ledger, "unallocated"), report.unallocated_energy_j, "unallocated");
        // The environment sink receives what leaves the cells; the report
        // counts the authored source as heat supplied to them.
        close(-sink(ledger, "environment"), report.external_thermal_heat_j, "boundary heat");
        close(
            node(ledger, "cells").stored_change_j.value(),
            report.deposited_heat_j + report.external_thermal_heat_j,
            "stored in the cells",
        );
        close(
            node(ledger, "cells").stored_change_j.value(),
            report.thermal.stored_energy_change_j,
            "the thermal owner's stored change",
        );
        close(
            node(ledger, "network").stored_change_j.value(),
            report.mechanical_energy_change_j,
            "mechanical energy",
        );
        close(ledger.energy_residual_j(), 0.0, "the ledger closes");
        if events {
            // 2 W to 0.35 s, 4 W after.
            let authored = 2.0 * 0.35 + 4.0 * (target_s - 0.35);
            close(report.external_thermal_heat_j, authored, "the authored source");
        } else {
            close(report.external_thermal_heat_j, 0.0, "no boundary heat");
        }
    }
}

/// The module's own synthetic temperature law, offered to the network
/// Component through that Component's law interface.
#[derive(Debug)]
struct ModuleLaw {
    beta_per_k: f64,
    maximum_k: f64,
}

impl ModuleLaw {
    fn source(&self) -> LinearTemperatureDamping {
        LinearTemperatureDamping {
            beta_per_k: self.beta_per_k,
            maximum_k: self.maximum_k,
        }
    }
}

impl DampingLaw for ModuleLaw {
    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, String> {
        self.source()
            .at_temperature(connector, reference, temperature_k)
            .map_err(str::to_owned)
    }
}

/// The feedback graph: one connector, all of its heat into cell 0, whose
/// start-of-step temperature `timeline` publishes and the network's law
/// reads. `timeline` has a heat port and a published temperature for cell 0
/// and the environment port. Returns the system and its two nodes, the
/// network and the cells.
fn feedback_graph(
    bodies: &Bodies,
    heat: &DepositionSimulation<'static>,
    timeline: ThermalTimeline<'static>,
    mechanical: ViscousNetworkMotionControls,
    law: &ModuleLaw,
) -> (System, NodeId, NodeId) {
    let network = ViscousNetwork::temperature_dependent(
        bodies.properties.to_vec(),
        bodies.links[..1].to_vec(),
        bodies.gravity,
        mechanical,
        std::rc::Rc::new(ModuleLaw {
            beta_per_k: law.beta_per_k,
            maximum_k: law.maximum_k,
        }),
    )
    .unwrap();
    let start = timeline.start(heat).unwrap();
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, bodies.states.to_vec()).unwrap();
    let cells = builder.add("cells", timeline, start).unwrap();
    for body in 0..bodies.loads.len() {
        let force = builder.port(moving, &format!("force{body}")).unwrap();
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(moving, "heat0").unwrap(),
            &[(builder.port(cells, "cell0").unwrap(), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    builder
        .heat_route(
            builder.port(cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    builder
        .signal(
            builder.port(cells, "temperature0").unwrap(),
            builder.port(moving, "temperature0").unwrap(),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    // The kernel solves the two inputs with the rest of the block and
    // starts from these.
    for name in ["heat0", "temperature0"] {
        let port = system.port(moving, name).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, moving, cells)
}

#[test]
fn the_declared_graph_reproduces_the_temperature_feedback_candidate() {
    // The fixture of the module's own feedback tests: one connector between
    // two 2 kg bodies, all of its heat into cell 0 (1 J/K), whose
    // temperature sets the connector's damping: c = 0.7 (1 + 2 (T - 300)).
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
    let bodies = Bodies::new([0.7, 0.0]);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let heat = thermal(laws, false, true);
    let law = ModuleLaw {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    let mut controls = config();
    controls.thermal.max_dt_s = 0.03125;
    for target_s in [0.5, 0.2] {
        let mut input = bodies.input(&routes);
        input.links = &bodies.links[..1];
        let candidate = advance_temperature(
            input,
            &heat,
            &(),
            &law.source(),
            &[CellId(0)],
            target_s,
            controls,
            || false,
        )
        .unwrap();
        let report = candidate.report();

        let timeline = ThermalTimeline::new(vec![CellId(0)], true, controls.thermal)
            .unwrap()
            .publishing(vec![CellId(0)])
            .unwrap();
        let (mut system, moving, cells) =
            feedback_graph(&bodies, &heat, timeline, controls.mechanical, &law);
        let intervals = run(&mut system, target_s);
        let case = format!("target {target_s}");

        assert_eq!(report.thermal.substeps, intervals.len(), "{case}: {intervals:?}");
        // The feedback is real: the cell warmed by more than half a kelvin,
        // so the damping more than doubled on the way.
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        let kernel_k = timeline.simulation().temperature_k(CellId(0)).unwrap();
        let module_k = candidate.thermal().temperature_k(CellId(0)).unwrap();
        assert!(module_k > 300.5, "{case}: {module_k}");
        assert!((kernel_k - module_k).abs() <= 1e-11, "{case}: {kernel_k} vs {module_k}");
        // The kernel solves the temperature signal with the heat in one
        // block, to a relative tolerance of 1e-12 (3e-10 K here), so the law
        // sees the module's temperature to that and the bodies agree to a
        // like tolerance. (In this fixture they come out equal bit for bit;
        // the kernel does not promise that.)
        let reached = system.state::<Vec<RigidBodyState>>(moving).unwrap();
        for (body, (ours, theirs)) in reached.iter().zip(candidate.states()).enumerate() {
            let (ours, theirs) = (
                ours.center_of_mass_velocity_world(bodies.properties[body]).unwrap(),
                theirs.center_of_mass_velocity_world(bodies.properties[body]).unwrap(),
            );
            assert!((ours.x - theirs.x).abs() <= 1e-9, "{case}: body {body}: {ours:?} vs {theirs:?}");
        }
        // The first body started at 1 m/s.
        let slowed = reached[0]
            .center_of_mass_velocity_world(bodies.properties[0])
            .unwrap();
        assert!(slowed.x < 0.9, "{case}: the first body lost speed, {slowed:?}");
        let ledger = system.ledger();
        let dissipated = report.dissipated_energy_j;
        assert!(dissipated > 0.5, "{case}: {dissipated} J");
        let close = |ours: f64, theirs: f64, what: &str| {
            assert!(
                (ours - theirs).abs() <= 1e-11 * dissipated,
                "{case}: {what}: ledger {ours:e}, report {theirs:e}"
            );
        };
        close(node(ledger, "network").dissipated_j.value(), dissipated, "dissipated");
        close(node(ledger, "cells").stored_change_j.value(), report.deposited_heat_j, "stored");
        close(
            node(ledger, "network").stored_change_j.value(),
            report.mechanical_energy_change_j,
            "mechanical energy",
        );
        close(ledger.energy_residual_j(), 0.0, "the ledger closes");
    }
}

/// The stiff pair's graph: the network, the timeline stepped by the
/// owner's implicit solve on a conduction graph prepared as the module's
/// is, and both connectors' heat routed whole to cell 0 (one receiver
/// shared by two routes), released at `released`. Returns the system and
/// its two nodes.
fn implicit_graph(
    bodies: &Bodies,
    laws: &'static [Caloric; 2],
    heat: &DepositionSimulation<'static>,
    config: RigidDampingImplicitConfig,
    scheme: ImplicitScheme,
    released: RouteTemperature,
) -> (System, NodeId, NodeId) {
    let network = ViscousNetwork::new(
        bodies.properties.to_vec(),
        bodies.links.to_vec(),
        bodies.gravity,
        config.mechanical,
    )
    .unwrap();
    let timeline =
        ThermalTimeline::implicit(vec![CellId(0)], stiff_graph(laws, heat), config.thermal, scheme)
            .unwrap();
    let start = timeline.start(heat).unwrap();
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, bodies.states.to_vec()).unwrap();
    let cells = builder.add("cells", timeline, start).unwrap();
    for body in 0..bodies.loads.len() {
        let force = builder.port(moving, &format!("force{body}")).unwrap();
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    for connector in 0..2 {
        builder
            .heat_route(
                builder.port(moving, &format!("heat{connector}")).unwrap(),
                &[(builder.port(cells, "cell0").unwrap(), 1.0)],
                None,
                released,
            )
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(cells, "environment").unwrap(),
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    if matches!(released, RouteTemperature::SinkEquivalent) {
        // The receivers' temperature comes back to the network, so the
        // kernel solves the two together and starts from this.
        for connector in 0..2 {
            let port = system.port(moving, &format!("heat{connector}")).unwrap();
            assert!(system.guess(port, &[300.0]));
        }
    }
    (system, moving, cells)
}

#[test]
fn the_declared_graph_follows_the_module_through_rejected_intervals() {
    // The stiff pair with calorics that curve (e = a dT + b dT^2 / 2) and a
    // solve allowed two Newton iterations. The module's solve rejects
    // intervals; it halves each one and re-runs its mechanics over the half.
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
    let bodies = Bodies::new([0.7, 1.1]);
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(0), 1.0)]),
    ];
    let heat = stiff_thermal(laws, 256);
    for (module_scheme, scheme) in [
        (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
        (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
    ] {
        let mut config = implicit_config(module_scheme);
        config.thermal.maximum_newton_iterations = 2;
        let conduction = stiff_graph(laws, &heat);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (candidate, native) = advance_implicit_candidate(
            bodies.input(&routes),
            &heat,
            &Observe::default(),
            0.5,
            &conduction,
            config,
            &mut workspace,
            || false,
        )
        .unwrap();
        let case = format!("{scheme:?}");
        assert!(native.rejected_intervals > 0, "{case}: {native:?}");
        let module_intervals: Vec<f64> =
            candidate.observer().intervals.iter().map(|&(_, dt_s, _)| dt_s).collect();
        assert_eq!(module_intervals.len(), native.energy.substeps, "{case}");

        // The graph's timeline answers `Shorten` where the module's solve
        // rejects, and the kernel re-steps the network with it: the same
        // accepted intervals, one shortening for every rejection. The heat
        // is released at a fixed temperature here, as in the module, whose
        // mechanics never sees a cell: the network then depends on nothing
        // the cells return, the kernel steps it first, and the timeline's
        // solve is tried with the heat it will be given.
        let (mut system, moving, cells) = implicit_graph(
            &bodies,
            laws,
            &heat,
            config,
            scheme,
            RouteTemperature::Fixed(300.0),
        );
        let (intervals, shortenings) = run_counting(&mut system, 0.5);
        assert_eq!(intervals, module_intervals, "{case}");
        assert_eq!(shortenings, native.rejected_intervals, "{case}");
        let reached = system.state::<Vec<RigidBodyState>>(moving).unwrap();
        assert_eq!(reached.as_slice(), candidate.states(), "{case}");
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        for id in [CellId(0), CellId(1)] {
            let kernel = timeline.simulation().temperature_k(id).unwrap();
            let module = candidate.thermal().temperature_k(id).unwrap();
            assert!((kernel - module).abs() <= 1e-9, "{case}: {id:?} {kernel} vs {module}");
        }
        let dissipated = candidate.report().dissipated_energy_j;
        let stored = node(system.ledger(), "cells").stored_change_j.value();
        assert!((stored - dissipated).abs() <= 1e-11 * dissipated, "{case}: {stored}");

        // Released at the receivers' temperature, the heat's temperature
        // comes back to the network, the kernel solves the two together,
        // and the timeline can answer `Shorten` on an iterate that is not
        // yet the heat it will be given. From one state that can only add
        // shortenings; here the run comes out finer than the module's
        // (measured: 23 steps and 76 shortenings against 20 and 64 for
        // backward Euler, 16 and 44 against 15 and 40 for SDIRK2) and
        // agrees with it as two discretisations do (measured: 6e-6 K).
        let (mut looped, _, cells) = implicit_graph(
            &bodies,
            laws,
            &heat,
            config,
            scheme,
            RouteTemperature::SinkEquivalent,
        );
        let (finer, more) = run_counting(&mut looped, 0.5);
        assert!(finer.len() >= module_intervals.len(), "{case}: {}", finer.len());
        assert!(more >= native.rejected_intervals, "{case}: {more}");
        let timeline = looped.state::<TimelineState<'static>>(cells).unwrap();
        for id in [CellId(0), CellId(1)] {
            let kernel = timeline.simulation().temperature_k(id).unwrap();
            let module = candidate.thermal().temperature_k(id).unwrap();
            assert!((kernel - module).abs() <= 1e-4, "{case}: {id:?} {kernel} vs {module}");
        }
    }
}

#[test]
fn the_declared_graph_reproduces_the_implicit_candidate() {
    // The fixture of the module's own implicit test: two cells on a
    // 1000 W/K link, which the explicit rule cannot step in its substep
    // allowance, and both connectors' heat into cell 0. A quarter-second
    // thermal interval, backward Euler and SDIRK2.
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
    let bodies = Bodies::new([0.7, 1.1]);
    let routes = [
        BTreeMap::from([(CellId(0), 1.0)]),
        BTreeMap::from([(CellId(0), 1.0)]),
    ];
    let heat = stiff_thermal(laws, 8);
    for (module_scheme, scheme) in [
        (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
        (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
    ] {
        let config = implicit_config(module_scheme);
        let conduction = stiff_graph(laws, &heat);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (candidate, native) = advance_implicit_candidate(
            bodies.input(&routes),
            &heat,
            &(),
            0.5,
            &conduction,
            config,
            &mut workspace,
            || false,
        )
        .unwrap();
        let report = candidate.report();
        assert_eq!((native.energy.substeps, native.rejected_intervals), (2, 0));

        let (mut system, moving, cells) = implicit_graph(
            &bodies,
            laws,
            &heat,
            config,
            scheme,
            RouteTemperature::SinkEquivalent,
        );
        let intervals = run(&mut system, 0.5);
        let case = format!("{scheme:?}");

        // Two quarter-second steps, where the explicit rule's bound is
        // under a millisecond: the same intervals as the module's.
        assert_eq!(intervals, [0.25, 0.25], "{case}");
        let reached = system.state::<Vec<RigidBodyState>>(moving).unwrap();
        assert_eq!(reached.as_slice(), candidate.states(), "{case}");
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        assert_eq!(
            timeline.simulation().time_s().to_bits(),
            candidate.thermal().time_s().to_bits(),
            "{case}: the two thermal clocks"
        );
        // Both sides close each cell to 1e-10 J on a capacity of 1 J/K.
        for id in [CellId(0), CellId(1)] {
            let kernel = timeline.simulation().temperature_k(id).unwrap();
            let module = candidate.thermal().temperature_k(id).unwrap();
            assert!((kernel - module).abs() <= 1e-9, "{case}: {id:?} {kernel} vs {module}");
            assert!(kernel > 300.5, "{case}: the damper heat warmed {id:?} to {kernel}");
        }
        let ledger = system.ledger();
        let dissipated = report.dissipated_energy_j;
        assert!(dissipated > 1.0, "{case}: enough heat to compare, {dissipated} J");
        let close = |ours: f64, theirs: f64, allowance: f64, what: &str| {
            assert!(
                (ours - theirs).abs() <= allowance,
                "{case}: {what}: ledger {ours:e}, report {theirs:e}"
            );
        };
        let rounding = 1e-11 * dissipated;
        // The tolerance the owner's solve gates its energy balance with.
        let closure = config.thermal.total_tolerance_j;
        close(node(ledger, "network").dissipated_j.value(), dissipated, rounding, "dissipated");
        // A property of the fixture: its routes leave nothing unallocated.
        close(report.deposited_heat_j, dissipated, rounding, "all of it deposited");
        close(
            node(ledger, "cells").stored_change_j.value(),
            report.thermal.stored_energy_change_j,
            closure,
            "stored in the cells",
        );
        close(node(ledger, "cells").stored_change_j.value(), dissipated, closure, "stored");
        close(sink(ledger, "environment"), 0.0, closure, "no boundary heat");
        close(
            node(ledger, "network").stored_change_j.value(),
            report.mechanical_energy_change_j,
            rounding,
            "mechanical energy",
        );
        close(ledger.energy_residual_j(), 0.0, closure, "the ledger closes");
    }
}

#[test]
fn the_declared_graph_reproduces_the_implicit_temperature_feedback_candidate() {
    // Temperature feedback on the stiff pair: one connector between two
    // 2 kg bodies heats cell 0, which a 1000 W/K link shares with cell 1,
    // and cell 0's temperature sets the connector's damping. Sixteenth-second
    // thermal intervals, backward Euler and SDIRK2.
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
    let bodies = Bodies::new([0.7, 0.0]);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let heat = stiff_thermal(laws, 8);
    let law = ModuleLaw {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    for (module_scheme, scheme) in [
        (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
        (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
    ] {
        let mut config = implicit_config(module_scheme);
        config.thermal.maximum_dt_s = 0.0625;
        let mut input = bodies.input(&routes);
        input.links = &bodies.links[..1];
        let conduction = stiff_graph(laws, &heat);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (candidate, native) =
            fs_couple::rigid_damping_heat::temperature::advance_implicit_candidate(
                input,
                &heat,
                &(),
                &law.source(),
                &[CellId(0)],
                0.5,
                &conduction,
                config,
                &mut workspace,
                || false,
            )
            .unwrap();
        let report = candidate.report();
        assert_eq!((native.energy.substeps, native.rejected_intervals), (8, 0));

        let timeline = ThermalTimeline::implicit(
            vec![CellId(0)],
            stiff_graph(laws, &heat),
            config.thermal,
            scheme,
        )
        .unwrap()
        .publishing(vec![CellId(0)])
        .unwrap();
        let (mut system, moving, cells) =
            feedback_graph(&bodies, &heat, timeline, config.mechanical, &law);
        let intervals = run(&mut system, 0.5);
        let case = format!("{scheme:?}");

        assert_eq!(intervals, [0.0625; 8], "{case}");
        let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
        // The feedback is real: cell 0 warmed, so the damping rose on the way.
        let module_k = candidate.thermal().temperature_k(CellId(0)).unwrap();
        assert!(module_k > 300.2, "{case}: {module_k}");
        for id in [CellId(0), CellId(1)] {
            let kernel = timeline.simulation().temperature_k(id).unwrap();
            let module = candidate.thermal().temperature_k(id).unwrap();
            assert!((kernel - module).abs() <= 1e-9, "{case}: {id:?} {kernel} vs {module}");
        }
        // The law sees the module's temperature to the kernel's solve
        // tolerance and to the closure of the two thermal solves.
        let reached = system.state::<Vec<RigidBodyState>>(moving).unwrap();
        for (body, (ours, theirs)) in reached.iter().zip(candidate.states()).enumerate() {
            let (ours, theirs) = (
                ours.center_of_mass_velocity_world(bodies.properties[body]).unwrap(),
                theirs.center_of_mass_velocity_world(bodies.properties[body]).unwrap(),
            );
            assert!((ours.x - theirs.x).abs() <= 1e-8, "{case}: body {body}: {ours:?} vs {theirs:?}");
        }
        let ledger = system.ledger();
        let dissipated = report.dissipated_energy_j;
        assert!(dissipated > 0.5, "{case}: {dissipated} J");
        let closure = 8.0 * config.thermal.total_tolerance_j;
        let close = |ours: f64, theirs: f64, what: &str| {
            assert!(
                (ours - theirs).abs() <= closure,
                "{case}: {what}: ledger {ours:e}, report {theirs:e}"
            );
        };
        close(node(ledger, "network").dissipated_j.value(), dissipated, "dissipated");
        close(node(ledger, "cells").stored_change_j.value(), report.deposited_heat_j, "stored");
        close(
            node(ledger, "network").stored_change_j.value(),
            report.mechanical_energy_change_j,
            "mechanical energy",
        );
        close(ledger.energy_residual_j(), 0.0, "the ledger closes");
    }
}
