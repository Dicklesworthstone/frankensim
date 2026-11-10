//! The declared graph against the hand-written coupling.
//!
//! `advance_candidate` is this module's explicit damping-heat step: the
//! thermal owner picks each interval, the damped network advances once over
//! it, and connector heat enters the authored cells as mean power. Here the
//! same two owners are kernel Components
//! (`fs_mbd::component::network::ViscousNetwork` and
//! `fs_thermal::component::ThermalTimeline`), the authored heat fractions are
//! heat routes, and the kernel in `fs-system` steps them. No code is written
//! for the combination. Where the module succeeds the two must agree: the
//! bodies bit for bit (to the kernel's solve tolerance when a temperature
//! signal feeds the damping), the cells to the rounding of handing heat over
//! as energy instead of as a rate, and the report's accounts with the
//! kernel's ledger. Not compared: the module refuses a whole advance when
//! its law has no form at the final temperature, where the graph takes that
//! step and refuses the next.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::{DampingLaw, ViscousNetwork};
use fs_system::{Ledger, NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ThermalTimeline, TimelineState};

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

fn sink(ledger: &Ledger, name: &str) -> f64 {
    ledger
        .sinks
        .iter()
        .find(|sink| sink.name == name)
        .map_or(0.0, |sink| sink.heat_j.value())
}

fn node<'a>(ledger: &'a Ledger, name: &str) -> &'a fs_system::ledger::NodeAccount {
    ledger.nodes.iter().find(|node| node.name == name).unwrap()
}

/// Step `system` to `target_s`, asking each time for all that is left: the
/// thermal Component's bound and events cut the step, as the thermal owner's
/// do inside `advance_candidate`.
fn run(system: &mut System, target_s: f64) -> Vec<f64> {
    let gate = CancelGate::new_clock_free();
    let mut intervals = Vec::new();
    with_cx(&gate, |cx| {
        while system.time_s() < target_s {
            let report = system.step(target_s - system.time_s(), cx).unwrap();
            intervals.push(report.dt_s);
            assert!(intervals.len() < 1000);
        }
    });
    intervals
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
/// start-of-step temperature the timeline publishes and the network's law
/// reads. Returns the system and its two nodes, the network and the cells.
fn feedback_graph(
    bodies: &Bodies,
    heat: &DepositionSimulation<'static>,
    controls: RigidDampingThermalConfig,
    law: &ModuleLaw,
) -> (System, NodeId, NodeId) {
    let network = ViscousNetwork::temperature_dependent(
        bodies.properties.to_vec(),
        bodies.links[..1].to_vec(),
        bodies.gravity,
        controls.mechanical,
        std::rc::Rc::new(ModuleLaw {
            beta_per_k: law.beta_per_k,
            maximum_k: law.maximum_k,
        }),
    )
    .unwrap();
    let timeline = ThermalTimeline::new(vec![CellId(0)], true, controls.thermal)
        .unwrap()
        .publishing(vec![CellId(0)])
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

        let (mut system, moving, cells) = feedback_graph(&bodies, &heat, controls, &law);
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
