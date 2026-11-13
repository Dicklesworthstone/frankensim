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
//! The temperature-dependent twins of these entry points
//! (`advance_temperature_candidate` and its implicit form, for the pair and
//! for the network) are compared too: the timeline publishes cell 0's
//! start-of-step temperature and the Component's damping law reads it on a
//! `temperature<e>` signal, one per connector. The law is the module's own
//! synthetic one on both sides, and so is the thermal owner; what the
//! feedback comparison adds is how the temperature reaches the law (a
//! published signal in place of a call) on top of the two mechanical
//! entries.

// Exact float comparisons here assert that nothing entered through an idle
// port.
#![allow(clippy::float_cmp)]

use super::*;
use crate::driven::{controls as pair_controls, effort, pair_input};
use crate::driven_network::{cycle_bodies, cycle_links, cycle_routes, graph_controls, graph_input, held};
use crate::graph::{ModuleLaw, node, run_counting, sink};
use fs_couple::rigid_damping_heat::RigidDampingThermalReport;
use fs_couple::rigid_damping_heat::driven::network::{
    DrivenNetworkDampingImplicitConfig, advance_candidate as advance_cycle,
    advance_implicit_candidate as advance_cycle_implicit,
    advance_temperature_candidate as advance_cycle_temperature,
    advance_temperature_implicit_candidate as advance_cycle_temperature_implicit,
};
use fs_couple::rigid_damping_heat::driven::{
    DrivenDampingImplicitConfig, advance_candidate as advance_pair,
    advance_implicit_candidate as advance_pair_implicit,
    advance_temperature_candidate as advance_pair_temperature,
    advance_temperature_implicit_candidate as advance_pair_temperature_implicit,
};
use fs_mbd::component::relaxation::{DriveFrame, RelaxationNetwork, RelaxationState};
use fs_mbd::spatial_loads::internal_impulse::viscous::network::driven::DrivenViscousNetworkLimits;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::{
    RelaxationNetworkLimits, RelaxationNetworkLink,
};
use fs_mbd::spatial_loads::internal_motion::viscous::network::relaxation::RelaxationNetworkMotionControls;
use fs_couple::rigid_damping_heat::temperature::TemperatureDampingThermalCandidate;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::StepConfig;
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionReport;

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
    /// With a law, every connector's Kelvin form follows it at cell 0's
    /// start-of-step temperature, which the timeline must publish.
    feedback: Option<&'a ModuleLaw>,
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
    if let Some(law) = parts.feedback {
        network = network
            .temperature_dependent(std::rc::Rc::new(ModuleLaw {
                beta_per_k: law.beta_per_k,
                maximum_k: law.maximum_k,
            }))
            .unwrap();
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
    if parts.feedback.is_some() {
        // One published temperature feeds every connector's law.
        for connector in 0..parts.links.len() {
            builder
                .signal(
                    builder.port(cells, "temperature0").unwrap(),
                    builder.port(moving, &format!("temperature{connector}")).unwrap(),
                )
                .unwrap();
        }
    }
    let mut system = builder.seal().unwrap();
    if parts.feedback.is_some() {
        // The kernel solves the signals with the rest of their block and
        // starts from these.
        for connector in 0..parts.links.len() {
            let port = system.port(moving, &format!("temperature{connector}")).unwrap();
            assert!(system.guess(port, &[300.0]));
        }
    }
    (system, moving, cells)
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

/// How closely a graph must follow the module: bodies [m, N s, N m s],
/// cells [K], energy accounts (relative).
#[derive(Clone, Copy)]
struct Allowance {
    body: f64,
    cell_k: f64,
    account: f64,
}

/// With fixed Kelvin forms.
const FIXED: Allowance = Allowance {
    body: BODY,
    cell_k: CELL_K,
    account: ACCOUNT,
};

/// With temperature feedback. Measured: 5.9e-15 for the bodies, 5.7e-14 K
/// for the cells (one ulp of 300 K), 5.2e-15 for the accounts, with cell 0
/// warming by 0.39 to 1.37 K. These bounds hold the kernel to what it does
/// today, delivering a start-of-step signal exactly; its contract is the
/// solve tolerance (1e-12 relative, 3e-10 K), at which the pair's Kelvin
/// form would move by 2e-10 and the bodies by about 5e-11 a step. If the
/// kernel ever uses that tolerance here, these become 1e-9, as in
/// `graph.rs`.
const FEEDBACK: Allowance = Allowance {
    body: 1e-12,
    cell_k: 1e-11,
    account: 1e-12,
};

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
    allowed: Allowance,
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
            position <= allowed.body && momentum <= allowed.body && spin <= allowed.body,
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
        assert!((kernel - module).abs() <= allowed.cell_k, "{case}: {id:?} {kernel} vs {module}");
    }

    let ledger = system.ledger();
    let dissipated = report.dissipated_energy_j;
    let driven: f64 = report.connector_actuator_work_j.iter().map(|work| work.abs()).sum();
    assert!(dissipated > 0.05, "{case}: enough heat to compare, {dissipated} J");
    assert!(driven > 0.05, "{case}: enough actuator work to compare, {driven} J");
    let scale = dissipated + driven;
    let close = |ours: f64, theirs: f64, what: &str| {
        assert!(
            (ours - theirs).abs() <= allowed.account * scale,
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
        feedback: None,
    }
}

/// The relaxation owner's controls for the module's pair.
fn pair_mechanical() -> RelaxationNetworkMotionControls {
    mechanical(
        pair_controls().mechanical,
        DrivenViscousNetworkLimits {
            maximum_bodies: 2,
            maximum_connectors: 1,
            maximum_live_bytes: 128_000_000,
            maximum_preparation_terms: 10_000_000_000,
        },
    )
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
        same_run(&case, &system, (moving, cells), &module, &forces, FIXED);

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
            same_run(&case, &system, (moving, cells), &module, &forces, FIXED);
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
            mechanical: pair_mechanical(),
            feedback: None,
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
        same_run(&case, &system, (moving, cells), &module, &forces, FIXED);

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
            same_run(&case, &system, (moving, cells), &module, &forces, FIXED);
        }
    }
}

/// One temperature-feedback fixture and the module's two entry points for
/// it.
struct FeedbackCase<'a, E, I> {
    name: &'a str,
    parts: Parts<'a>,
    /// Where the explicit run ends: sixteen thermal steps, on a grid whose
    /// steps sum to it exactly.
    target_s: f64,
    /// The module's explicit entry, given the thermal state, the thermal
    /// step limits and the target.
    explicit: E,
    /// The module's implicit entry, run to half a second in two steps.
    implicit: I,
}

/// The graph of `parts` run to `target_s` on `timeline`, compared with the
/// module's candidate. Returns the kernel's intervals and shortenings.
fn follow(
    parts: &Parts<'_>,
    (label, target_s): (&str, f64),
    timeline: ThermalTimeline<'static>,
    heat: &DepositionSimulation<'static>,
    module: &Module<'_>,
) -> (Vec<f64>, usize) {
    let (mut system, moving, cells) = driven_graph(parts, timeline, heat);
    let counted = run_counting(&mut system, target_s);
    // The feedback is real: cell 0 warmed, so the Kelvin forms changed on
    // the way.
    let warmed = module.heat.temperature_k(CellId(0)).unwrap() - 300.0;
    assert!(warmed > 0.05, "{label}: cell 0 warmed by {warmed} K");
    same_run(label, &system, (moving, cells), module, parts.forces, FEEDBACK);
    counted
}

/// `case` against the module: explicit thermal steps, then the stiff pair
/// of cells by the implicit solve in both schemes.
fn follow_with_feedback<E, I>(case: &FeedbackCase<'_, E, I>)
where
    E: Fn(&DepositionSimulation<'static>, StepConfig, f64) -> TemperatureDampingThermalCandidate<'static, ()>,
    I: Fn(
        &DepositionSimulation<'static>,
        RigidDampingImplicitConfig,
        &ImplicitConductionGraph<'_>,
        &mut ImplicitConductionWorkspace,
    ) -> (TemperatureDampingThermalCandidate<'static, Observe>, ImplicitConductionReport),
{
    let laws: &'static [Caloric; 2] = Box::leak(Box::new(linear_calorics()));
    // All the heat goes to cell 0: the timeline has a heat port for it
    // alone, and the kernel wants every declared port connected.
    let cells = || vec![CellId(0)];

    let heat = thermal(laws, false, true);
    let mut step = config().thermal;
    step.max_dt_s = case.target_s / 16.0;
    let candidate = (case.explicit)(&heat, step, case.target_s);
    // The module evaluated its law at more than the starting temperature.
    assert!(candidate.law_evaluations() > 1, "{}: explicit", case.name);
    let timeline = ThermalTimeline::new(cells(), true, step)
        .unwrap()
        .publishing(vec![CellId(0)])
        .unwrap();
    let module = Module {
        states: candidate.states(),
        heat: candidate.thermal(),
        report: candidate.report(),
    };
    let label = format!("{}, explicit", case.name);
    let (intervals, shortenings) =
        follow(&case.parts, (&label, case.target_s), timeline, &heat, &module);
    assert_eq!(candidate.report().thermal.substeps, intervals.len(), "{label}");
    assert_eq!(shortenings, 0, "{label}");

    // Two implicit steps of a quarter of a second, so that the second
    // step's forms are read at a warmed cell.
    let heat = stiff_thermal(laws, 8);
    for (module_scheme, scheme) in [
        (ImplicitThermalScheme::BackwardEuler, ImplicitScheme::BackwardEuler),
        (ImplicitThermalScheme::Sdirk2, ImplicitScheme::Sdirk2),
    ] {
        let common = implicit_config(module_scheme);
        let conduction = stiff_graph(laws, &heat);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (candidate, native) = (case.implicit)(&heat, common, &conduction, &mut workspace);
        assert!(candidate.law_evaluations() > 1, "{}: {scheme:?}", case.name);
        let timeline =
            ThermalTimeline::implicit(cells(), stiff_graph(laws, &heat), common.thermal, scheme)
                .unwrap()
                .publishing(vec![CellId(0)])
                .unwrap();
        let module = Module {
            states: candidate.states(),
            heat: candidate.thermal(),
            report: candidate.report(),
        };
        let label = format!("{}, {scheme:?}", case.name);
        let (intervals, shortenings) = follow(&case.parts, (&label, 0.5), timeline, &heat, &module);
        assert_eq!(intervals, [0.25, 0.25], "{label}");
        assert_eq!(native.energy.substeps, intervals.len(), "{label}");
        assert_eq!(shortenings, native.rejected_intervals, "{label}");
    }
}

#[test]
fn the_declared_graph_follows_the_driven_pair_with_temperature_feedback() {
    // The module's own feedback fixture: 0.7 N s/m, all of its heat into
    // cell 0 (1 J/K), c = 0.7 (1 + 2 (T - 300)), 0.6 N held.
    let bodies = Bodies::new([0.7, 0.0]);
    let routes = [BTreeMap::from([(CellId(0), 1.0)])];
    let law = ModuleLaw {
        beta_per_k: 2.0,
        maximum_k: 303.0,
    };
    let forces = [0.6];
    follow_with_feedback(&FeedbackCase {
        name: "pair",
        parts: Parts {
            properties: &bodies.properties[..2],
            states: &bodies.states[..2],
            links: &bodies.links[..1],
            forces: &forces,
            routes: &routes,
            mechanical: pair_mechanical(),
            feedback: Some(&law),
        },
        target_s: 0.5,
        explicit: |heat: &DepositionSimulation<'static>, step: StepConfig, target_s: f64| {
            let mut controls = pair_controls();
            controls.thermal = step;
            advance_pair_temperature(
                pair_input(&bodies, &routes),
                effort(forces[0]),
                heat,
                &(),
                &law.source(),
                &[CellId(0)],
                target_s,
                controls,
                || false,
            )
            .unwrap()
        },
        implicit: |heat: &DepositionSimulation<'static>,
                   common: RigidDampingImplicitConfig,
                   conduction: &ImplicitConductionGraph<'_>,
                   workspace: &mut ImplicitConductionWorkspace| {
            let implicit = DrivenDampingImplicitConfig {
                mechanical: pair_controls().mechanical,
                thermal: common.thermal,
                scheme: common.scheme,
                maximum_heat_bindings: common.maximum_heat_bindings,
                energy_tolerance_j: common.energy_tolerance_j,
            };
            advance_pair_temperature_implicit(
                pair_input(&bodies, &routes),
                effort(forces[0]),
                heat,
                &Observe::default(),
                &law.source(),
                &[CellId(0)],
                0.5,
                conduction,
                implicit,
                workspace,
                || false,
            )
            .unwrap()
        },
    });
}

#[test]
fn the_declared_graph_follows_the_driven_network_with_temperature_feedback() {
    // The module's cycle with all heat into cell 0 and
    // c = 0.7 (1 + (T - 300) / 2), forces (0.2, -0.1, 0.4) N.
    let bodies = cycle_bodies();
    let links = cycle_links();
    let routes = core::array::from_fn::<_, 3, _>(|_| BTreeMap::from([(CellId(0), 1.0)]));
    let law = ModuleLaw {
        beta_per_k: 0.5,
        maximum_k: 303.0,
    };
    let forces = [0.2, -0.1, 0.4];
    let efforts = [held(forces[0]), held(forces[1]), held(forces[2])];
    follow_with_feedback(&FeedbackCase {
        name: "cycle",
        parts: Parts {
            feedback: Some(&law),
            ..cycle_parts(&bodies, &links, &forces, &routes)
        },
        // The module's own test runs to 0.3 s in sixteen steps. Sixteen
        // steps of 0.3 / 16 s do not sum to 0.3 s in floating point: the
        // module's thermal clock lands on its target and the kernel's sum
        // ends one ulp short (the limit stated in `graph.rs`). A quarter of
        // a second in sixteen steps is exact.
        target_s: 0.25,
        explicit: |heat: &DepositionSimulation<'static>, step: StepConfig, target_s: f64| {
            let mut controls = graph_controls();
            controls.thermal = step;
            advance_cycle_temperature(
                graph_input(&bodies, &links, &routes),
                &efforts,
                heat,
                &(),
                &law.source(),
                &[CellId(0); 3],
                target_s,
                controls,
                || false,
            )
            .unwrap()
        },
        implicit: |heat: &DepositionSimulation<'static>,
                   common: RigidDampingImplicitConfig,
                   conduction: &ImplicitConductionGraph<'_>,
                   workspace: &mut ImplicitConductionWorkspace| {
            let implicit = DrivenNetworkDampingImplicitConfig {
                mechanical: graph_controls().mechanical,
                thermal: common.thermal,
                scheme: common.scheme,
                maximum_heat_bindings: common.maximum_heat_bindings,
                energy_tolerance_j: common.energy_tolerance_j,
            };
            advance_cycle_temperature_implicit(
                graph_input(&bodies, &links, &routes),
                &efforts,
                heat,
                &Observe::default(),
                &law.source(),
                &[CellId(0); 3],
                0.5,
                conduction,
                implicit,
                workspace,
                || false,
            )
            .unwrap()
        },
    });
}
