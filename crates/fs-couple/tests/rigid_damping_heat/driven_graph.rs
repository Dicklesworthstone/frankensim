//! The declared graph against the hand-written driven coupling.
//!
//! `driven::advance_candidate` (one pair) and `driven::network::advance_candidate`
//! (any connector graph) hold an actuator effort on each connector while the
//! connector's viscous loss heats the thermal cells; both have an implicit
//! twin. Here the mechanical owner is the kernel Component
//! `fs_mbd::component::relaxation::RelaxationNetwork` with no Maxwell branch
//! and an actuator on every connector (`driven`), the held efforts are
//! boundaries on its `drive<e>` ports, the thermal owner is
//! `fs_thermal::component::ThermalTimeline`, and the authored heat fractions
//! are heat routes.
//!
//! The two sides enter fs-mbd by different doors. The module's pair entry
//! calls the driven-pair owner and its network entry the driven-network
//! owner; the Component calls the relaxation-network owner, which solves
//! the same connectors (Kelvin viscosity and held effort) as one system
//! with an empty set of branches. The three owners share fs-mbd's connector
//! geometry and rigid drift, and the two network owners the application of
//! their impulses, so a defect there would sit on both sides and this
//! comparison would not see it; what differs is how each assembles and
//! solves its source. The bodies therefore agree to rounding (measured
//! below) and not by bits.
//! The module's report is compared with the kernel's ledger: each
//! actuator's signed work is the energy that entered through its drive
//! boundary, and only the viscous loss is heat.
//!
//! Not compared: the temperature-dependent twins of these entry points
//! (the Component's Kelvin viscosity is fixed).

// Exact float comparisons here assert that nothing entered through an idle
// port.
#![allow(clippy::float_cmp)]

use super::*;
use crate::driven::{controls as pair_controls, effort, pair_input};
use crate::driven_network::{cycle_bodies, cycle_links, cycle_routes, graph_controls, graph_input, held};
use crate::graph::{node, run_counting, sink};
use fs_couple::rigid_damping_heat::RigidDampingThermalReport;
use fs_couple::rigid_damping_heat::driven::network::{
    DrivenNetworkDampingImplicitConfig, advance_candidate as advance_cycle,
    advance_implicit_candidate as advance_cycle_implicit,
};
use fs_couple::rigid_damping_heat::driven::{
    DrivenDampingImplicitConfig, advance_candidate as advance_pair,
    advance_implicit_candidate as advance_pair_implicit,
};
use fs_mbd::component::relaxation::{DriveFrame, RelaxationNetwork, RelaxationState};
use fs_mbd::spatial_loads::internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::{
    RelaxationNetworkLimits, RelaxationNetworkLink,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::RelaxationNetworkMotionControls;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};

/// What one driven graph is made of.
struct Parts<'a> {
    properties: &'a [MassProperties],
    states: &'a [RigidBodyState],
    links: &'a [ViscousNetworkLink],
    /// The force each connector's actuator holds along the first body's x
    /// axis [N].
    forces: &'a [f64],
    routes: &'a [BTreeMap<CellId, f64>],
    mechanical: RelaxationNetworkMotionControls,
}

/// The relaxation owner's controls for a network of this size with no
/// branch, from the connector controls the module's entry points use.
fn mechanical(
    source: ViscousConnectorControls,
    limits: DrivenViscousNetworkLimits,
) -> RelaxationNetworkMotionControls {
    RelaxationNetworkMotionControls {
        source,
        limits: RelaxationNetworkLimits {
            source: limits,
            maximum_branches: 0,
        },
    }
}

/// The graph: bodies and connectors with an actuator each, the held forces
/// as boundaries, every connector's heat routed to the timeline's cells by
/// its authored fractions and the rest to a sink named `unallocated`.
/// Returns the system, the network and the cells.
fn driven_graph(
    parts: &Parts<'_>,
    timeline: ThermalTimeline<'static>,
    heat: &DepositionSimulation<'static>,
) -> (System, NodeId, NodeId) {
    let links = parts
        .links
        .iter()
        .map(|link| RelaxationNetworkLink {
            viscous: link.clone(),
            branches: Vec::new(),
        })
        .collect();
    let mut network =
        RelaxationNetwork::new(parts.properties.to_vec(), links, Gravity::ZERO, parts.mechanical)
            .unwrap();
    for connector in 0..parts.links.len() {
        network = network.driven(connector, DriveFrame::MaterialFirst).unwrap();
    }
    let start = network.relaxed(parts.states.to_vec());
    let cells_start = timeline.start(heat).unwrap();
    let mut builder = SystemBuilder::new();
    let moving = builder.add("network", network, start).unwrap();
    let cells = builder.add("cells", timeline, cells_start).unwrap();
    for body in 0..parts.states.len() {
        let force = builder.port(moving, &format!("force{body}")).unwrap();
        builder
            .boundary(force, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    for (connector, &force) in parts.forces.iter().enumerate() {
        let drive = builder.port(moving, &format!("drive{connector}")).unwrap();
        builder
            .boundary(
                drive,
                Box::new(move |_, _, out: &mut [f64]| out.copy_from_slice(&[force, 0.0, 0.0])),
            )
            .unwrap();
        let couple = builder.port(moving, &format!("couple{connector}")).unwrap();
        builder
            .boundary(couple, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .unwrap();
    }
    for (connector, fractions) in parts.routes.iter().enumerate() {
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
                RouteTemperature::Fixed(300.0),
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
    (builder.seal().unwrap(), moving, cells)
}

/// What the module's candidate holds.
struct Module<'a> {
    states: &'a [RigidBodyState],
    heat: &'a DepositionSimulation<'static>,
    report: &'a RigidDampingThermalReport,
}

/// Allowance on a body's position [m] and momenta [N s, N m s] between the
/// two mechanical owners over a run of half a second. Measured in these
/// fixtures: 1.7e-15 at most (the pair), 2.2e-16 (the cycle).
const BODY: f64 = 1e-13;

/// Allowance on a cell temperature [K]. Measured: the cells come out equal
/// to the bit in every case here.
const CELL_K: f64 = 1e-12;

/// Allowance on an energy account, relative to the dissipation plus the
/// actuators' work in magnitude. Measured: 7.3e-15 at most (the pair),
/// 4.5e-15 (the cycle).
const ACCOUNT: f64 = 1e-12;

fn apart(a: Vec3, b: Vec3) -> f64 {
    (a.x - b.x).abs().max((a.y - b.y).abs()).max((a.z - b.z).abs())
}

/// The module's candidate against the graph's reached state and ledger.
fn same_run(
    case: &str,
    system: &System,
    (moving, cells): (NodeId, NodeId),
    module: &Module<'_>,
    forces: &[f64],
) {
    let Module {
        states,
        heat,
        report,
    } = *module;
    let reached = system.state::<RelaxationState>(moving).unwrap();
    assert_eq!(reached.bodies.len(), states.len(), "{case}");
    for (body, (ours, theirs)) in reached.bodies.iter().zip(states).enumerate() {
        let position = apart(ours.pose().position_world(), theirs.pose().position_world());
        let momentum = apart(ours.linear_momentum_world(), theirs.linear_momentum_world());
        let spin = apart(ours.angular_momentum_body(), theirs.angular_momentum_body());
        assert!(
            position <= BODY && momentum <= BODY && spin <= BODY,
            "{case}: body {body}: position {position:e}, momentum {momentum:e}, spin {spin:e}"
        );
        // Nothing turns in these fixtures, on either side.
        assert_eq!(ours.pose().orientation(), theirs.pose().orientation(), "{case}");
    }
    let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
    assert_eq!(
        timeline.simulation().time_s().to_bits(),
        heat.time_s().to_bits(),
        "{case}: the two thermal clocks"
    );
    for id in [CellId(0), CellId(1)] {
        let kernel = timeline.simulation().temperature_k(id).unwrap();
        let module = heat.temperature_k(id).unwrap();
        assert!((kernel - module).abs() <= CELL_K, "{case}: {id:?} {kernel} vs {module}");
    }

    let ledger = system.ledger();
    let dissipated = report.dissipated_energy_j;
    let driven: f64 = report.connector_actuator_work_j.iter().map(|work| work.abs()).sum();
    assert!(dissipated > 0.05, "{case}: enough heat to compare, {dissipated} J");
    assert!(driven > 0.05, "{case}: enough actuator work to compare, {driven} J");
    let scale = dissipated + driven;
    let close = |ours: f64, theirs: f64, what: &str| {
        assert!(
            (ours - theirs).abs() <= ACCOUNT * scale,
            "{case}: {what}: ledger {ours:e}, report {theirs:e}"
        );
    };
    close(node(ledger, "network").dissipated_j.value(), dissipated, "dissipated");
    close(sink(ledger, "unallocated"), report.unallocated_energy_j, "unallocated");
    close(
        node(ledger, "cells").stored_change_j.value(),
        report.deposited_heat_j + report.external_thermal_heat_j,
        "stored in the cells",
    );
    close(
        node(ledger, "network").stored_change_j.value(),
        report.mechanical_energy_change_j,
        "mechanical energy",
    );
    // Each actuator's signed work entered through its own drive boundary,
    // and nothing came in anywhere else.
    let boundary = |port: String| {
        ledger
            .boundaries
            .iter()
            .find(|account| account.port == port)
            .unwrap_or_else(|| panic!("{case}: boundary {port}"))
            .energy_in_j
            .value()
    };
    assert_eq!(report.connector_actuator_work_j.len(), forces.len(), "{case}");
    for (connector, &work) in report.connector_actuator_work_j.iter().enumerate() {
        close(boundary(format!("network.drive{connector}")), work, "an actuator's work");
        assert_eq!(boundary(format!("network.couple{connector}")), 0.0, "{case}");
    }
    for body in 0..states.len() {
        assert_eq!(boundary(format!("network.force{body}")), 0.0, "{case}");
    }
    close(ledger.boundary_energy_in_j(), report.actuator_work_j, "all actuator work");
    close(ledger.energy_residual_j(), 0.0, "the ledger closes");
}

/// The cycle of the module's network tests: 2 kg bodies moving at 1, 0 and
/// -1 m/s along x, joined 0-1, 1-2 and 0-2 by axial connectors of
/// 0.7 N s/m with held forces `(F, -F / 2, 2F)`, and a fourth body joined
/// to nothing.
fn cycle_parts<'a>(
    bodies: &'a Bodies,
    links: &'a [ViscousNetworkLink; 3],
    forces: &'a [f64; 3],
    routes: &'a [BTreeMap<CellId, f64>; 3],
) -> Parts<'a> {
    let module = graph_controls().mechanical;
    Parts {
        properties: &bodies.properties,
        states: &bodies.states,
        links,
        forces,
        routes,
        mechanical: mechanical(module.source, module.limits),
    }
}

#[test]
fn the_declared_graph_follows_the_driven_network_candidate() {
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = cycle_routes();
    let input = graph_input(&bodies, &links, &routes);
    for force in [-0.35, 0.4] {
        // Three different forces, so that no two connectors do the same
        // work and a swapped account would show.
        let forces = [force, -0.5 * force, 2.0 * force];
        let efforts = [held(forces[0]), held(forces[1]), held(forces[2])];
        let parts = cycle_parts(&bodies, &links, &forces, &routes);

        // Explicit thermal steps of an eighth of a second.
        let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
        let heat = thermal(laws, false, true);
        let controls = graph_controls();
        let candidate =
            advance_cycle(input, &efforts, &heat, &(), 0.5, controls, || false).unwrap();
        let timeline =
            ThermalTimeline::new(vec![CellId(0), CellId(1)], true, controls.thermal).unwrap();
        let (mut system, moving, cells) = driven_graph(&parts, timeline, &heat);
        let (intervals, shortenings) = run_counting(&mut system, 0.5);
        let case = format!("cycle, explicit, force {force}");
        assert_eq!(candidate.report().thermal.substeps, intervals.len(), "{case}");
        assert_eq!(shortenings, 0, "{case}");
        let module = Module {
            states: candidate.states(),
            heat: candidate.thermal(),
            report: candidate.report(),
        };
        same_run(&case, &system, (moving, cells), &module, &forces);

        // The stiff pair of cells by the implicit solve, both schemes.
        let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
        let heat = stiff_thermal(laws, 8);
        for (module_scheme, scheme) in [
            (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
            (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
        ] {
            let common = implicit_config(module_scheme);
            let implicit = DrivenNetworkDampingImplicitConfig {
                mechanical: graph_controls().mechanical,
                thermal: common.thermal,
                scheme: module_scheme,
                maximum_heat_bindings: common.maximum_heat_bindings,
                energy_tolerance_j: common.energy_tolerance_j,
            };
            let conduction = stiff_graph(laws, &heat);
            let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
            let (candidate, native) = advance_cycle_implicit(
                input,
                &efforts,
                &heat,
                &Observe::default(),
                0.5,
                &conduction,
                implicit,
                &mut workspace,
                || false,
            )
            .unwrap();
            let module_intervals: Vec<f64> =
                candidate.observer().intervals.iter().map(|&(_, dt_s, _)| dt_s).collect();
            let timeline = ThermalTimeline::implicit(
                vec![CellId(0), CellId(1)],
                stiff_graph(laws, &heat),
                implicit.thermal,
                scheme,
            )
            .unwrap();
            let (mut system, moving, cells) = driven_graph(&parts, timeline, &heat);
            let (intervals, shortenings) = run_counting(&mut system, 0.5);
            let case = format!("cycle, {scheme:?}, force {force}");
            assert_eq!(intervals, module_intervals, "{case}");
            assert_eq!(shortenings, native.rejected_intervals, "{case}");
            let module = Module {
                states: candidate.states(),
                heat: candidate.thermal(),
                report: candidate.report(),
            };
            same_run(&case, &system, (moving, cells), &module, &forces);
        }
    }
}

#[test]
fn the_declared_graph_follows_the_driven_pair_candidate() {
    // The module's pair: two 2 kg bodies closing at 2 m/s, one axial
    // connector of 0.7 N s/m, a quarter of its heat to cell 0 and half to
    // cell 1.
    let bodies = Bodies::new([0.7, 0.0]);
    let routes = [BTreeMap::from([(CellId(0), 0.25), (CellId(1), 0.5)])];
    let input = pair_input(&bodies, &routes);
    for force in [-1.4, 1.4] {
        let forces = [force];
        let parts = Parts {
            properties: &bodies.properties[..2],
            states: &bodies.states[..2],
            links: &bodies.links[..1],
            forces: &forces,
            routes: &routes,
            mechanical: mechanical(
                pair_controls().mechanical,
                DrivenViscousNetworkLimits {
                    maximum_bodies: 2,
                    maximum_connectors: 1,
                    maximum_live_bytes: 128_000_000,
                    maximum_preparation_terms: 10_000_000_000,
                },
            ),
        };

        let laws: &'static [Caloric; 2] = Box::leak(Box::new(laws()));
        let heat = thermal(laws, false, true);
        let controls = pair_controls();
        let candidate =
            advance_pair(input, effort(force), &heat, &(), 0.5, controls, || false).unwrap();
        let timeline =
            ThermalTimeline::new(vec![CellId(0), CellId(1)], true, controls.thermal).unwrap();
        let (mut system, moving, cells) = driven_graph(&parts, timeline, &heat);
        let (intervals, shortenings) = run_counting(&mut system, 0.5);
        let case = format!("pair, explicit, force {force}");
        assert_eq!(candidate.report().thermal.substeps, intervals.len(), "{case}");
        assert_eq!(shortenings, 0, "{case}");
        let module = Module {
            states: candidate.states(),
            heat: candidate.thermal(),
            report: candidate.report(),
        };
        same_run(&case, &system, (moving, cells), &module, &forces);

        let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
        let heat = stiff_thermal(laws, 8);
        for (module_scheme, scheme) in [
            (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
            (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
        ] {
            let common = implicit_config(module_scheme);
            let implicit = DrivenDampingImplicitConfig {
                mechanical: pair_controls().mechanical,
                thermal: common.thermal,
                scheme: module_scheme,
                maximum_heat_bindings: common.maximum_heat_bindings,
                energy_tolerance_j: common.energy_tolerance_j,
            };
            let conduction = stiff_graph(laws, &heat);
            let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
            let (candidate, native) = advance_pair_implicit(
                input,
                effort(force),
                &heat,
                &Observe::default(),
                0.5,
                &conduction,
                implicit,
                &mut workspace,
                || false,
            )
            .unwrap();
            let module_intervals: Vec<f64> =
                candidate.observer().intervals.iter().map(|&(_, dt_s, _)| dt_s).collect();
            let timeline = ThermalTimeline::implicit(
                vec![CellId(0), CellId(1)],
                stiff_graph(laws, &heat),
                implicit.thermal,
                scheme,
            )
            .unwrap();
            let (mut system, moving, cells) = driven_graph(&parts, timeline, &heat);
            let (intervals, shortenings) = run_counting(&mut system, 0.5);
            let case = format!("pair, {scheme:?}, force {force}");
            assert_eq!(intervals, module_intervals, "{case}");
            assert_eq!(shortenings, native.rejected_intervals, "{case}");
            let module = Module {
                states: candidate.states(),
                heat: candidate.thermal(),
                report: candidate.report(),
            };
            same_run(&case, &system, (moving, cells), &module, &forces);
        }
    }
}
