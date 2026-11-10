//! Owner slice: the fs-mbd damped rigid-body network on the composition kernel.
//!
//! `fs_mbd::component::network::ViscousNetwork` puts the owner's stepper for
//! bodies joined by viscous connectors behind the Component contract, with a
//! heat port per connector. Here it passes the conformance check, is shown to
//! add nothing to the owner's step, and is joined by the kernel to thermal
//! bodies through heat routes with fractions: the "rigid damping heat"
//! pairing as a declared graph, with no code written for the pair.
//! Run with `cargo test -p fs-system --features slices`.

// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::{
    ViscousNetwork, ViscousNetworkBuildError, ViscousNetworkRefusal,
};
use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
        network::{ViscousNetworkLimits, ViscousNetworkLink},
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads,
    viscous::network::{
        ViscousNetworkMotionControls, ViscousNetworkMotionStep, step_viscous_network_with_loads,
    },
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, Proposal};
use fs_system::elements::HeatCapacity;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};

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

fn damping(diagonal: [f64; 6]) -> PassiveBushingMatrix {
    let mut c = [[0.0; 6]; 6];
    for (axis, value) in diagonal.into_iter().enumerate() {
        c[axis][axis] = value;
    }
    PassiveBushingMatrix::try_new(
        c,
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
    )
    .expect("a diagonal non-negative damping form is passive")
}

fn body(p: MassProperties, x: Vec3, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, UnitQuaternion::IDENTITY).expect("pose"),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).expect("spin"),
    )
    .expect("state")
}

/// Three bodies in a row. The first connector damps translation on every
/// axis and rotation about one, through arms that are off the centres of
/// mass; the second damps along the row only.
fn fixture() -> (Vec<MassProperties>, Vec<ViscousNetworkLink>, Vec<RigidBodyState>) {
    fixture_at(Vec3::ZERO)
}

/// The same bodies with every position moved by `shift`.
fn fixture_at(
    shift: Vec3,
) -> (Vec<MassProperties>, Vec<ViscousNetworkLink>, Vec<RigidBodyState>) {
    let properties = vec![
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).expect("properties"),
        MassProperties::new(1.5, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.5)).expect("properties"),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(4.0, 2.0, 3.0)).expect("properties"),
    ];
    let links = vec![
        ViscousNetworkLink {
            bodies: [0, 1],
            arms_body_m: [Vec3::new(0.3, 0.1, 0.0), Vec3::new(-0.2, 0.05, 0.1)],
            damping: damping([0.6, 0.4, 0.3, 0.0, 0.05, 0.0]),
        },
        ViscousNetworkLink {
            bodies: [1, 2],
            arms_body_m: [Vec3::ZERO; 2],
            damping: damping([0.9, 0.0, 0.0, 0.0, 0.0, 0.0]),
        },
    ];
    let states = vec![
        body(
            properties[0],
            Vec3::new(-1.0, 0.0, 0.2).add(shift),
            Vec3::new(1.2, 0.1, 0.0),
            Vec3::new(0.0, 0.3, -0.2),
        ),
        body(
            properties[1],
            Vec3::new(0.0, 0.1, 0.0).add(shift),
            Vec3::new(-0.8, 0.0, 0.2),
            Vec3::new(0.1, 0.0, 0.4),
        ),
        body(
            properties[2],
            Vec3::new(1.0, -0.1, 0.0).add(shift),
            Vec3::new(0.2, -0.3, 0.0),
            Vec3::ZERO,
        ),
    ];
    (properties, links, states)
}

fn network(gravity: Gravity) -> (ViscousNetwork, Vec<RigidBodyState>) {
    let (properties, links, states) = fixture();
    (
        ViscousNetwork::new(properties, links, gravity, controls()).expect("network"),
        states,
    )
}

fn owner_step(
    states: &[RigidBodyState],
    wrenches: &[Wrench],
    gravity: Gravity,
    dt_s: f64,
) -> ViscousNetworkMotionStep {
    let (properties, links, _) = fixture();
    step_viscous_network_with_loads(
        states,
        &properties,
        &links,
        HeldExternalLoads { wrenches, gravity },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step")
}

fn norm3(v: &[f64]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn force(x: f64, y: f64, z: f64) -> Wrench {
    Wrench {
        force_world: Vec3::new(x, y, z),
        torque_body: Vec3::ZERO,
    }
}

#[test]
fn the_damped_network_conforms_and_books_the_owners_work_and_heat() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("gravity");
    let (component, start) = network(gravity);
    assert_eq!((component.bodies(), component.connectors()), (3, 2));
    assert_eq!(component.ports().len(), 5);
    // A stale candidate of another length: the component must write all of it.
    let stale = vec![start[2]];
    for (dt_s, scale) in [(1.0e-3, 1.0), (5.0e-3, -2.0), (2.0e-2, 0.5)] {
        let inputs = [
            1.5 * scale,
            -0.5 * scale,
            2.0 * scale,
            0.0,
            0.7 * scale,
            -0.3 * scale,
            -0.4 * scale,
            0.0,
            0.2 * scale,
            305.0,
            288.0,
        ];
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        let observed =
            with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
        let Observed::Stepped {
            receipt,
            audit,
            outputs,
        } = observed
        else {
            panic!("the network takes the step");
        };
        // The owner's own account of the same step, called directly.
        let wrenches = [
            force(inputs[0], inputs[1], inputs[2]),
            force(inputs[3], inputs[4], inputs[5]),
            force(inputs[6], inputs[7], inputs[8]),
        ];
        let native = owner_step(&start, &wrenches, gravity, dt_s);
        let heat_j: f64 = native.connector_dissipation_j.iter().sum();
        assert!(heat_j > 0.0, "the connectors dissipate");
        let work_j = native.external_loads.applied_force_work_j;
        // Sizes of the terms the kernel sums, which may cancel in the total.
        let paired_j: f64 = (0..3)
            .map(|body| {
                dt_s * norm3(&inputs[3 * body..3 * body + 3]) * norm3(&outputs[3 * body..3 * body + 3])
            })
            .sum();
        let exchanged_j = paired_j + heat_j;
        assert!(work_j.abs() > 0.01 * paired_j, "the forces do net work");
        assert!(
            (audit.port_energy_in_j - (work_j - heat_j)).abs() <= 1.0e-12 * exchanged_j,
            "{} vs {} - {}",
            audit.port_energy_in_j,
            work_j,
            heat_j
        );
        assert_eq!(outputs[9], -native.connector_dissipation_j[0] / dt_s);
        assert_eq!(outputs[10], -native.connector_dissipation_j[1] / dt_s);
        assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
        assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
        let entropy = native.connector_dissipation_j[0] / 305.0
            + native.connector_dissipation_j[1] / 288.0;
        assert!((receipt.entropy_generated_j_per_k - entropy).abs() <= 1.0e-15 * entropy);
        // What the owner's splitting fails to conserve is reported, and small.
        assert!(receipt.numerical_loss_j.abs() <= 1.0e-9 * (exchanged_j + receipt.stored_change_j.abs()));
    }
}

#[test]
fn the_wrapper_adds_nothing_to_the_owners_step() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let (component, start) = network(gravity);
    let inputs = [0.3, 0.0, -0.1, 0.0, 0.0, 0.0, -0.2, 0.4, 0.0, 300.0, 300.0];
    let dt_s = 4.0e-3;
    let mut next = Vec::new();
    let mut outputs = [f64::NAN; 11];
    let proposal = with_cx(&gate, |cx| {
        component.propose(&start, dt_s, &inputs, &mut next, &mut outputs, cx)
    })
    .expect("step");
    assert!(matches!(proposal, Proposal::Step(_)));
    let wrenches = [
        force(0.3, 0.0, -0.1),
        force(0.0, 0.0, 0.0),
        force(-0.2, 0.4, 0.0),
    ];
    let native = owner_step(&start, &wrenches, gravity, dt_s);
    assert_eq!(next.as_slice(), native.states_after());
    for (index, state) in next.iter().enumerate() {
        let (a, b) = (
            start[index].pose().position_world(),
            state.pose().position_world(),
        );
        assert_eq!(outputs[3 * index], (b.x - a.x) / dt_s);
        assert_eq!(outputs[3 * index + 1], (b.y - a.y) / dt_s);
        assert_eq!(outputs[3 * index + 2], (b.z - a.z) / dt_s);
    }
}

#[test]
fn the_network_refuses_what_it_cannot_represent() {
    let gate = CancelGate::new_clock_free();
    let (properties, links, start) = fixture();
    assert_eq!(
        ViscousNetwork::new(Vec::new(), Vec::new(), Gravity::ZERO, controls()).unwrap_err(),
        ViscousNetworkBuildError::NoBodies
    );
    let mut wrong = links.clone();
    wrong[1].bodies = [1, 3];
    assert_eq!(
        ViscousNetwork::new(properties.clone(), wrong, Gravity::ZERO, controls()).unwrap_err(),
        ViscousNetworkBuildError::Link { connector: 1 }
    );
    let mut same = links.clone();
    same[0].bodies = [2, 2];
    assert_eq!(
        ViscousNetwork::new(properties.clone(), same, Gravity::ZERO, controls()).unwrap_err(),
        ViscousNetworkBuildError::Link { connector: 0 }
    );

    let component =
        ViscousNetwork::new(properties, links, Gravity::ZERO, controls()).expect("network");
    let mut inputs = [0.0; 11];
    inputs[9] = 300.0;
    inputs[10] = 0.0;
    let mut next = start.clone();
    let mut outputs = [0.0; 11];
    let refused = with_cx(&gate, |cx| {
        component.propose(&start, 1.0e-3, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(
        refused,
        ViscousNetworkRefusal::Temperature {
            connector: 1,
            kelvin: 0.0
        }
    );
    inputs[10] = 300.0;
    let short = vec![start[0], start[1]];
    let refused = with_cx(&gate, |cx| {
        component.propose(&short, 1.0e-3, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(
        refused,
        ViscousNetworkRefusal::Bodies {
            expected: 3,
            found: 2
        }
    );
    assert!(component.stored_energy_j(&short).is_err());

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let refused = with_cx(&cancelled, |cx| {
        component.propose(&start, 1.0e-3, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(refused, ViscousNetworkRefusal::Cancelled);
}

const FRACTION_A: f64 = 0.6;
const FRACTION_B: f64 = 0.3;
const CAPACITY_A_J_K: f64 = 5.0;
const CAPACITY_B_J_K: f64 = 8.0;

struct DampingHeat {
    system: System,
    network: NodeId,
    cell_a: NodeId,
    cell_b: NodeId,
}

/// The damped network with its first connector's heat split 60 % to cell A,
/// 30 % to cell B and the rest to a named sink, and its second connector's
/// heat all to cell B. Each body's force port is held by a boundary.
fn damping_heat(gravity: Gravity, held: [[f64; 3]; 3]) -> DampingHeat {
    let (component, start) = network(gravity);
    let mut builder = SystemBuilder::new();
    let network = builder.add("network", component, start).expect("network");
    let cell_a = builder
        .add_with_tangent("cell_a", HeatCapacity::new(CAPACITY_A_J_K).expect("cell"), 300.0)
        .expect("cell a");
    let cell_b = builder
        .add_with_tangent("cell_b", HeatCapacity::new(CAPACITY_B_J_K).expect("cell"), 310.0)
        .expect("cell b");
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).expect("port");
    for (index, value) in held.into_iter().enumerate() {
        let target = port(&builder, network, &format!("force{index}"));
        builder
            .boundary(
                target,
                Box::new(move |_t0, _dt, out: &mut [f64]| out.copy_from_slice(&value)),
            )
            .expect("boundary");
    }
    builder
        .heat_route(
            port(&builder, network, "heat0"),
            &[
                (port(&builder, cell_a, "heat"), FRACTION_A),
                (port(&builder, cell_b, "heat"), FRACTION_B),
            ],
            Some(heat_sink("unallocated", 293.15).expect("sink")),
            RouteTemperature::SinkEquivalent,
        )
        .expect("route 0");
    builder
        .heat_route(
            port(&builder, network, "heat1"),
            &[(port(&builder, cell_b, "heat"), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .expect("route 1");
    let mut system = builder.seal().expect("sealed");
    for name in ["heat0", "heat1"] {
        let heat = system.port(network, name).expect("heat port");
        assert!(system.guess(heat, &[305.0]));
    }
    DampingHeat {
        system,
        network,
        cell_a,
        cell_b,
    }
}

fn account<'a>(system: &'a System, name: &str) -> &'a fs_system::ledger::NodeAccount {
    system
        .ledger()
        .nodes
        .iter()
        .find(|node| node.name == name)
        .expect("account")
}

#[test]
fn connector_heat_reaches_the_cells_by_the_declared_fractions() {
    let gate = CancelGate::new_clock_free();
    let mut graph = damping_heat(Gravity::ZERO, [[0.0; 3]; 3]);
    let dt_s = 5.0e-3;
    let steps = 400_usize;

    // The owner marched alone: the motion does not depend on where its heat
    // goes, so the kernel must reproduce it bit for bit.
    let (_, _, mut twin) = fixture();
    let idle = [force(0.0, 0.0, 0.0); 3];
    let (mut heat_0, mut heat_1) = (0.0_f64, 0.0_f64);
    let mut mechanical_j = 0.0_f64;
    for _ in 0..steps {
        let native = owner_step(&twin, &idle, Gravity::ZERO, dt_s);
        heat_0 += native.connector_dissipation_j[0];
        heat_1 += native.connector_dissipation_j[1];
        mechanical_j += native.mechanical_energy_change_j;
        twin = native.into_states_after();
        with_cx(&gate, |cx| graph.system.step(dt_s, cx)).expect("kernel step");
    }
    let reached = graph
        .system
        .state::<Vec<RigidBodyState>>(graph.network)
        .expect("network state");
    assert_eq!(reached, &twin);

    let total = heat_0 + heat_1;
    assert!(total > 0.5, "enough heat to measure: {total} J");
    // No load: the mechanical energy lost is the heat released, to the
    // owner's own splitting residual.
    assert!((mechanical_j + total).abs() <= 1.0e-7 * total, "{mechanical_j} vs {total}");
    let warmed_a = account(&graph.system, "cell_a").stored_change_j.value();
    let warmed_b = account(&graph.system, "cell_b").stored_change_j.value();
    let lost = graph.system.ledger().sink_heat_j();
    let close = |actual: f64, expected: f64| {
        assert!(
            (actual - expected).abs() <= 1.0e-11 * total,
            "{actual} vs {expected}"
        );
    };
    close(warmed_a, FRACTION_A * heat_0);
    close(warmed_b, FRACTION_B * heat_0 + heat_1);
    close(lost, (1.0 - FRACTION_A - FRACTION_B) * heat_0);
    let network = account(&graph.system, "network");
    close(network.dissipated_j.value(), total);
    close(network.stored_change_j.value(), mechanical_j);
    // The cells' states agree with their accounts.
    let cell = |node: NodeId| *graph.system.state::<f64>(node).expect("temperature");
    assert!(((cell(graph.cell_a) - 300.0) * CAPACITY_A_J_K - warmed_a).abs() <= 1.0e-9 * total);
    assert!(((cell(graph.cell_b) - 310.0) * CAPACITY_B_J_K - warmed_b).abs() <= 1.0e-9 * total);

    let ledger = graph.system.ledger();
    assert_eq!(ledger.steps, steps as u64);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * total);
    // Heat released near 305 K: the entropy generated is close to Q / T.
    let entropy = ledger.entropy_generated_j_per_k();
    assert!(entropy > total / 320.0 && entropy < total / 290.0, "{entropy}");
}

#[test]
fn a_held_force_under_gravity_enters_as_boundary_work() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.81)).expect("gravity");
    let held = [[0.5, 0.0, 19.62], [0.0, 0.0, 14.715], [0.0, -0.2, 29.43]];
    let mut graph = damping_heat(gravity, held);
    let dt_s = 2.0e-3;
    let steps = 300_usize;

    let (_, _, mut twin) = fixture();
    let wrenches = held.map(|f| force(f[0], f[1], f[2]));
    let (mut work_j, mut heat_j, mut mechanical_j) = (0.0_f64, 0.0_f64, 0.0_f64);
    for _ in 0..steps {
        let native = owner_step(&twin, &wrenches, gravity, dt_s);
        work_j += native.external_loads.applied_force_work_j;
        heat_j += native.dissipated_energy_j;
        mechanical_j += native.mechanical_energy_change_j;
        twin = native.into_states_after();
        with_cx(&gate, |cx| graph.system.step(dt_s, cx)).expect("kernel step");
    }
    let reached = graph
        .system
        .state::<Vec<RigidBodyState>>(graph.network)
        .expect("network state");
    assert_eq!(reached, &twin);

    let ledger = graph.system.ledger();
    let exchanged_j = work_j.abs() + heat_j;
    assert!(work_j.abs() > 0.1, "the held forces do work: {work_j} J");
    // The kernel booked the owner's applied work from the wire values.
    assert!(
        (ledger.boundary_energy_in_j() - work_j).abs() <= 1.0e-11 * exchanged_j,
        "{} vs {}",
        ledger.boundary_energy_in_j(),
        work_j
    );
    // First law of the whole graph: work in = mechanical change + heat, and
    // every joule of heat is in a cell or in the named sink.
    let stored_heat = account(&graph.system, "cell_a").stored_change_j.value()
        + account(&graph.system, "cell_b").stored_change_j.value()
        + ledger.sink_heat_j();
    assert!((stored_heat - heat_j).abs() <= 1.0e-11 * exchanged_j);
    assert!((work_j - mechanical_j - heat_j).abs() <= 1.0e-8 * exchanged_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * exchanged_j);
}

/// March the network alone, checking every step against the conformance kit
/// (which includes the audit of the endpoint energies).
fn march_checked(
    component: &ViscousNetwork,
    mut state: Vec<RigidBodyState>,
    dt_s: f64,
    steps: usize,
    what: &str,
) -> Vec<RigidBodyState> {
    let gate = CancelGate::new_clock_free();
    let mut inputs = [0.0; 11];
    inputs[9] = 300.0;
    inputs[10] = 300.0;
    let stale = state.clone();
    let mut outputs = [0.0; 11];
    with_cx(&gate, |cx| {
        for step in 0..steps {
            let case = Case {
                state: &state,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            check_component(component, &case, cx)
                .unwrap_or_else(|failure| panic!("{what}, step {step}: {failure}"));
            let mut next = Vec::new();
            component
                .propose(&state, dt_s, &inputs, &mut next, &mut outputs, cx)
                .unwrap_or_else(|refusal| panic!("{what}, step {step}: {refusal}"));
            state = next;
        }
    });
    state
}

#[test]
fn the_network_keeps_its_account_through_free_fall_and_under_oblique_gravity() {
    // Free fall for three seconds: kinetic and potential energy grow against
    // each other, and a tolerance sized by their difference refuses.
    let g = 9.81;
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -g)).expect("gravity");
    let (component, start) = network(gravity);
    let (properties, _, _) = fixture();
    let centre = |states: &[RigidBodyState]| -> (f64, f64) {
        let (mut mass, mut height, mut momentum) = (0.0, 0.0, 0.0);
        for (state, body) in states.iter().zip(&properties) {
            mass += body.mass();
            height += body.mass() * state.pose().position_world().z;
            momentum += state.linear_momentum_world().z;
        }
        (height / mass, momentum / mass)
    };
    let (z0, v0) = centre(&start);
    let fallen = march_checked(&component, start, 1.0e-3, 3_000, "free fall");
    let (z, v) = centre(&fallen);
    let t = 3.0;
    assert!((v - (v0 - g * t)).abs() <= 1.0e-10 * g * t, "{v}");
    assert!((z - (z0 + v0 * t - 0.5 * g * t * t)).abs() <= 1.0e-10 * g * t * t, "{z}");

    // An oblique field, a kilometre from the origin along a level surface:
    // the three terms of m g . x are large and cancel.
    let oblique = Gravity::new(Vec3::new(6.9367, 0.0, -6.9367)).expect("gravity");
    let (properties, links, far) = fixture_at(Vec3::new(1000.0, 0.0, 1000.0));
    let component = ViscousNetwork::new(properties, links, oblique, controls()).expect("network");
    march_checked(&component, far, 1.0e-3, 300, "oblique gravity");
}
