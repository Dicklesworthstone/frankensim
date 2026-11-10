//! The declared graph against the hand-written coupling.
//!
//! `advance_candidate` is this module's explicit damping-heat step: the
//! thermal owner picks each interval, the damped network advances once over
//! it, and connector heat enters the authored cells as mean power. Here the
//! same two owners are kernel Components
//! (`fs_mbd::component::network::ViscousNetwork` and
//! `fs_thermal::component::ThermalTimeline`), the authored heat fractions are
//! heat routes, and the kernel in `fs-system` steps them. No code is written
//! for the combination. The two must agree: the bodies bit for bit, the
//! cells to the rounding of handing heat over as energy instead of as a
//! rate, and the report's accounts with the kernel's ledger.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::ViscousNetwork;
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
