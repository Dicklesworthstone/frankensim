//! The declared graph against the hand-written Maxwell-memory coupling.
//!
//! `relaxation::advance_candidate` and `advance_implicit_candidate` advance
//! bodies, the complete elastic memory and the thermal cells together. Here
//! the mechanical owner is the kernel Component
//! `fs_mbd::component::relaxation::RelaxationNetwork` (bodies and memory in
//! one state), the thermal owner is `fs_thermal::component::ThermalTimeline`
//! (explicit or implicit), the authored heat fractions are a heat route, and
//! the kernel in `fs-system` steps them. The two must agree: bodies and
//! memory bit for bit (on these steps, where the kernel's step is the
//! thermal clock's `end - start`; see the parent `graph.rs`), the cells to
//! the rounding of handing heat over as energy (in this fixture they come
//! out equal to the bit with either driver), and the module's report with
//! the kernel's ledger.
//!
//! Each comparison is made three times: with the fixture's held internal
//! effort zero and a Component that has no actuator, and with 0.7 N and
//! -0.7 N held on the connector (against the closing motion and with it),
//! where the Component has an actuator (`driven`) and the effort is a
//! boundary on its `drive0` port. The module's signed actuator work is then
//! the energy that entered through that boundary. Both sides call the same
//! owner here, so that agreement shows the port's rate is conjugate to the
//! effort; it is not a second computation of the work. The heat is
//! released at a fixed temperature, as in the module, whose mechanics never
//! sees a cell.

use super::*;
use crate::graph::{node, run_counting, sink};
use fs_mbd::component::relaxation::{DriveFrame, RelaxationNetwork, RelaxationState};
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::RelaxationNetworkLink;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};

/// 0.7 of the connector's heat to cell 0, 0.2 to cell 1, the rest
/// unallocated: the module's own routes for this fixture.
fn routes() -> [BTreeMap<CellId, f64>; 1] {
    [BTreeMap::from([(CellId(0), 0.7), (CellId(1), 0.2)])]
}

/// The fixture's pair as a graph with `timeline` (heat ports on cells 0 and
/// 1, and the environment port). With `held` the connector has an actuator
/// holding that force [N] along the first body's x axis; without, it has
/// none. Returns the system, the network and the cells.
fn memory_graph(
    fixture: &MemoryPair,
    mechanical: RelaxationNetworkMotionControls,
    timeline: ThermalTimeline<'static>,
    heat: &DepositionSimulation<'static>,
    held: Option<f64>,
) -> (System, NodeId, NodeId) {
    let mut network = RelaxationNetwork::new(
        fixture.bodies.properties[..2].to_vec(),
        vec![RelaxationNetworkLink {
            viscous: fixture.links[0].clone(),
            branches: fixture.branches[0].clone(),
        }],
        fixture.bodies.gravity,
        mechanical,
    )
    .unwrap();
    if held.is_some() {
        network = network.driven(0, DriveFrame::MaterialFirst).unwrap();
    }
    let start = RelaxationState {
        bodies: fixture.bodies.states[..2].to_vec(),
        memory: fixture.memory.to_vec(),
    };
    let cells_start = timeline.start(heat).unwrap();
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, start).unwrap();
    let cells = builder.add("cells", timeline, cells_start).unwrap();
    for body in 0..2 {
        let force = builder.port(moving, &format!("force{body}")).unwrap();
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    if let Some(force) = held {
        let drive = builder.port(moving, "drive0").unwrap();
        builder
            .boundary(
                drive,
                Box::new(move |_, _, out: &mut [f64]| out.copy_from_slice(&[force, 0.0, 0.0])),
            )
            .unwrap();
        let couple = builder.port(moving, "couple0").unwrap();
        builder
            .boundary(couple, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    builder
        .heat_route(
            builder.port(moving, "heat0").unwrap(),
            &[
                (builder.port(cells, "cell0").unwrap(), 0.7),
                (builder.port(cells, "cell1").unwrap(), 0.2),
            ],
            Some(heat_sink("unallocated", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
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
    (builder.seal().unwrap(), moving, cells)
}

/// What the module's candidate holds.
struct Module<'a> {
    states: &'a [RigidBodyState],
    memory: &'a [Vec<[f64; 6]>],
    heat: &'a DepositionSimulation<'static>,
    report: &'a fs_couple::rigid_damping_heat::RigidDampingThermalReport,
}

/// The module's candidate against the graph's reached state and ledger.
fn same_run(
    case: &str,
    system: &System,
    (moving, cells): (NodeId, NodeId),
    module: &Module<'_>,
    cell_k: f64,
) {
    let Module {
        states,
        memory,
        heat,
        report,
    } = *module;
    let reached = system.state::<RelaxationState>(moving).unwrap();
    assert_eq!(reached.bodies.as_slice(), states, "{case}: bodies");
    assert_eq!(reached.memory.as_slice(), memory, "{case}: memory");
    let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
    assert_eq!(
        timeline.simulation().time_s().to_bits(),
        heat.time_s().to_bits(),
        "{case}: the two thermal clocks"
    );
    for id in [CellId(0), CellId(1)] {
        let kernel = timeline.simulation().temperature_k(id).unwrap();
        let module = heat.temperature_k(id).unwrap();
        assert!((kernel - module).abs() <= cell_k, "{case}: {id:?} {kernel} vs {module}");
    }
    let ledger = system.ledger();
    let dissipated = report.dissipated_energy_j;
    assert!(dissipated > 0.1, "{case}: enough heat to compare, {dissipated} J");
    let close = |ours: f64, theirs: f64, what: &str| {
        assert!(
            (ours - theirs).abs() <= 1e-11 * dissipated,
            "{case}: {what}: ledger {ours:e}, report {theirs:e}"
        );
    };
    close(node(ledger, "network").dissipated_j.value(), dissipated, "dissipated");
    close(sink(ledger, "unallocated"), report.unallocated_energy_j, "unallocated");
    close(
        node(ledger, "cells").stored_change_j.value(),
        report.deposited_heat_j,
        "stored in the cells",
    );
    // Kinetic and elastic-memory storage together.
    close(
        node(ledger, "network").stored_change_j.value(),
        report.mechanical_energy_change_j,
        "mechanical and memory energy",
    );
    assert!(report.memory_energy_change_j.abs() > 1e-3, "{case}: the memory moved");
    // The actuator's signed work is what entered through its boundary; a
    // graph without an actuator has no such boundary and the module's
    // actuator did nothing.
    let drive = ledger.boundaries.iter().find(|account| account.port == "network.drive0");
    match drive {
        Some(account) => {
            assert!(report.actuator_work_j.abs() > 0.05, "{case}: {} J", report.actuator_work_j);
            close(account.energy_in_j.value(), report.actuator_work_j, "the actuator's work");
        }
        // Exactly: a zero effort does no work.
        None => assert!(report.actuator_work_j.abs() <= 0.0, "{case}"),
    }
    close(ledger.boundary_energy_in_j(), report.actuator_work_j, "all boundary work");
    close(ledger.energy_residual_j(), 0.0, "the ledger closes");
}

/// The fixture's held force and the graph's actuator: none, then 0.7 N
/// and -0.7 N.
const HELD: [(f64, Option<f64>); 3] = [(0.0, None), (0.7, Some(0.7)), (-0.7, Some(-0.7))];

#[test]
fn the_declared_graph_reproduces_the_memory_candidate() {
    // The module's pair: Kelvin viscosity 0.4, one Maxwell branch (factor
    // 1.3, time 0.2 s) preloaded to 0.6 sqrt(J), closing at 2 m/s, without
    // and with a held force. Explicit thermal steps of an eighth of a
    // second.
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
    let routes = routes();
    let heat = thermal(laws, false, true);
    let controls = memory_controls();
    for ((force, held), target_s) in HELD.into_iter().flat_map(|held| [(held, 0.25), (held, 0.3)]) {
        let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, force);
        let candidate = advance_memory(
            fixture.input(&routes),
            &heat,
            &(),
            target_s,
            controls,
            || false,
        )
        .unwrap();
        let timeline =
            ThermalTimeline::new(vec![CellId(0), CellId(1)], true, controls.thermal).unwrap();
        let (mut system, moving, cells) =
            memory_graph(&fixture, controls.mechanical, timeline, &heat, held);
        let (intervals, shortenings) = run_counting(&mut system, target_s);
        let case = format!("explicit, force {force}, target {target_s}");
        assert_eq!(candidate.report().thermal.substeps, intervals.len(), "{case}");
        // Neither Component of an explicit graph shortens a step.
        assert_eq!(shortenings, 0, "{case}");
        let module = Module {
            states: candidate.states(),
            memory: candidate.memory(),
            heat: candidate.thermal(),
            report: candidate.report(),
        };
        same_run(&case, &system, (moving, cells), &module, 1e-11);
        // The fixture is the module's: its candidate is the pair's own
        // equations, as the module's tests have them.
        let expected = reference(0.4, 1.3, 0.2, force, -2.0, 0.6, target_s);
        near(relative(candidate.states(), &fixture.bodies), expected[0], 8e-11);
        near(candidate.memory()[0][0][0], expected[1], 8e-11);
    }
}

/// The implicit comparison's allowance on a cell temperature. Both sides
/// run the same solve on heat that differs by rounding at most; in this
/// fixture the cells come out equal to the bit (measured).
const CELL_K: f64 = 1e-10;

#[test]
fn the_declared_graph_reproduces_the_implicit_memory_candidate() {
    // The same pair on the stiff thermal fixture, by the implicit solve.
    // With two Newton iterations allowed the solve rejects intervals: the
    // module halves each one and restarts bodies and memory from the
    // accepted banks; the graph's timeline answers `Shorten`.
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
    let routes = routes();
    let heat = stiff_thermal(laws, 128);
    let schemes = [
        (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
        (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
    ];
    for ((force, held), (module_scheme, scheme)) in
        HELD.into_iter().flat_map(|held| schemes.map(|scheme| (held, scheme)))
    {
        let fixture = MemoryPair::new(0.4, 1.3, 0.2, 0.6, -2.0, force);
        let mut controls = memory_implicit_controls(module_scheme);
        controls.thermal.maximum_newton_iterations = 2;
        controls.thermal.cell_tolerance_j = 1e-8;
        let conduction = stiff_graph(laws, &heat);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (candidate, native) = advance_memory_implicit(
            fixture.input(&routes),
            &heat,
            &Observe::default(),
            0.25,
            &conduction,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        let case = format!("{scheme:?}, force {force}");
        assert!(native.rejected_intervals > 0, "{case}: {native:?}");
        let module_intervals: Vec<f64> =
            candidate.observer().intervals.iter().map(|&(_, dt_s, _)| dt_s).collect();

        let timeline = ThermalTimeline::implicit(
            vec![CellId(0), CellId(1)],
            stiff_graph(laws, &heat),
            controls.thermal,
            scheme,
        )
        .unwrap();
        let (mut system, moving, cells) =
            memory_graph(&fixture, controls.mechanical, timeline, &heat, held);
        let (intervals, shortenings) = run_counting(&mut system, 0.25);
        assert_eq!(intervals, module_intervals, "{case}");
        assert_eq!(shortenings, native.rejected_intervals, "{case}");
        let module = Module {
            states: candidate.states(),
            memory: candidate.memory(),
            heat: candidate.thermal(),
            report: candidate.report(),
        };
        same_run(&case, &system, (moving, cells), &module, CELL_K);
    }
}
