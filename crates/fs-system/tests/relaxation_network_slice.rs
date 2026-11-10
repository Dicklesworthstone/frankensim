//! Owner slice: the fs-mbd rigid-body network with Maxwell memory on the
//! composition kernel.
//!
//! `fs_mbd::component::relaxation::RelaxationNetwork` puts the owner's
//! `step_relaxation_network_with_loads` behind the Component contract: the
//! bodies and the complete elastic memory are one state, and the memory is
//! stored energy. Here it passes the conformance check, equals the owner
//! driven directly, releases its preloaded memory as motion and heat, and
//! heats thermal bodies through routes; with an actuator on a connector its
//! drive ports carry the owner's actuator work, in material and in world
//! axes; and its Kelvin form, or its viscosity and relaxation times
//! together, follow a temperature signal.
//! Run with `cargo test -p fs-system --features slices`.

// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::ViscousNetwork;
use std::rc::Rc;

use fs_mbd::component::relaxation::{
    DampingLaw, DriveFrame, RelaxationNetwork, RelaxationNetworkBuildError,
    RelaxationNetworkRefusal, RelaxationState, ShiftLaw,
};
use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls,
    viscous::{
        BushingCoordinateScales, PassiveBushingMatrix, ViscousConnectorControls,
        network::{
            ViscousNetworkLimits, ViscousNetworkLink,
            driven::DrivenViscousNetworkLimits,
            relaxation::{
                RelaxationNetworkLimits, RelaxationNetworkLink, SpatialRelaxationBranch,
            },
        },
    },
};
use fs_mbd::spatial_loads::internal_motion::{
    HeldExternalLoads,
    viscous::network::{
        ViscousNetworkMotionControls,
        relaxation::{
            HeldInternalEffort, HeldRelaxationKinetics, RelaxationNetworkMotionControls,
            RelaxationNetworkMotionStep, step_relaxation_network_with_loads,
            step_relaxation_network_with_loads_and_kinetics,
        },
    },
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::PortKind;
use fs_ports::component::{Carrier, Component, Proposal, Receipt};
use fs_system::elements::HeatCapacity;
use fs_system::{RouteTemperature, SystemBuilder, heat_sink};

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

fn source() -> ViscousConnectorControls {
    ViscousConnectorControls {
        impulse: InternalWrenchImpulseControls {
            angular_origin_world_m: Vec3::new(0.2, -0.3, 0.1),
            linear_tolerance_ns: 3e-10,
            angular_tolerance_nms: 3e-10,
            work_tolerance_j: 3e-10,
        },
        linear_speed_tolerance_m_s: 3e-10,
        angular_speed_tolerance_rad_s: 3e-10,
    }
}

fn controls() -> RelaxationNetworkMotionControls {
    RelaxationNetworkMotionControls {
        source: source(),
        limits: RelaxationNetworkLimits {
            source: DrivenViscousNetworkLimits {
                maximum_bodies: 8,
                maximum_connectors: 8,
                maximum_live_bytes: 128_000_000,
                maximum_preparation_terms: 10_000_000_000,
            },
            maximum_branches: 8,
        },
    }
}

fn scales() -> BushingCoordinateScales {
    BushingCoordinateScales {
        translation_m: 1.0,
        rotation_rad: 1.0,
    }
}

fn damping(diagonal: [f64; 6]) -> PassiveBushingMatrix {
    let mut c = [[0.0; 6]; 6];
    for (axis, value) in diagonal.into_iter().enumerate() {
        c[axis][axis] = value;
    }
    PassiveBushingMatrix::try_new(c, scales()).expect("a diagonal damping form is passive")
}

/// A Maxwell branch whose factor rows are `roots` on the diagonal, so that
/// its stiffness on axis `i` is `roots[i]^2`, relaxing in `tau_s`.
fn branch(roots: [f64; 6], tau_s: f64) -> SpatialRelaxationBranch {
    let mut rows = [[0.0; 6]; 6];
    for (axis, value) in roots.into_iter().enumerate() {
        rows[axis][axis] = value;
    }
    SpatialRelaxationBranch::from_gram_rows(rows, tau_s).expect("branch")
}

fn body(p: MassProperties, x: Vec3, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, UnitQuaternion::IDENTITY).expect("pose"),
        Vec3::new(p.mass() * v.x, p.mass() * v.y, p.mass() * v.z),
        w,
    )
    .expect("state")
}

fn properties() -> Vec<MassProperties> {
    vec![
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).expect("properties"),
        MassProperties::new(1.5, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.5)).expect("properties"),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(4.0, 2.0, 3.0)).expect("properties"),
    ]
}

fn viscous_links() -> Vec<ViscousNetworkLink> {
    vec![
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
    ]
}

/// Three bodies and two connectors: the first with two Maxwell branches
/// (one of them also rotational), the second with one.
fn links() -> Vec<RelaxationNetworkLink> {
    let [first, second] = <[ViscousNetworkLink; 2]>::try_from(viscous_links()).expect("two links");
    vec![
        RelaxationNetworkLink {
            viscous: first,
            branches: vec![
                branch([1.3, 0.0, 0.0, 0.0, 0.0, 0.0], 0.2),
                branch([0.0, 0.8, 0.0, 0.0, 0.5, 0.0], 0.05),
            ],
        },
        RelaxationNetworkLink {
            viscous: second,
            branches: vec![branch([0.9, 0.0, 0.0, 0.0, 0.0, 0.0], 0.5)],
        },
    ]
}

fn moving() -> Vec<RigidBodyState> {
    let p = properties();
    vec![
        body(p[0], Vec3::new(-1.0, 0.0, 0.2), Vec3::new(1.2, 0.1, 0.0), Vec3::new(0.0, 0.3, -0.2)),
        body(p[1], Vec3::new(0.0, 0.1, 0.0), Vec3::new(-0.8, 0.0, 0.2), Vec3::new(0.1, 0.0, 0.4)),
        body(p[2], Vec3::new(1.0, -0.1, 0.0), Vec3::new(0.2, -0.3, 0.0), Vec3::ZERO),
    ]
}

/// A preloaded memory: 0.045 + 0.02 + 0.0125 + 0.08 = 0.1575 J stored.
fn preloaded() -> Vec<Vec<[f64; 6]>> {
    vec![
        vec![
            [0.3, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.2, 0.0, 0.0, 0.05, 0.15],
        ],
        vec![[-0.4, 0.0, 0.0, 0.0, 0.0, 0.0]],
    ]
}

const PRELOAD_J: f64 = 0.5 * (0.09 + 0.04 + 0.0025 + 0.0225 + 0.16);

/// Two 2 kg bodies on the x axis joined by one axial connector with Kelvin
/// viscosity `PAIR_C` and one Maxwell branch (factor `PAIR_ROOT_K`, time
/// `PAIR_TAU_S`). The owner's source is exact for this geometry, so a run
/// can be checked against the two-variable equations by hand.
const PAIR_C: f64 = 0.4;
const PAIR_ROOT_K: f64 = 1.3;
const PAIR_TAU_S: f64 = 0.2;

fn pair_properties() -> Vec<MassProperties> {
    vec![MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).expect("properties"); 2]
}

fn pair_links() -> Vec<RelaxationNetworkLink> {
    pair_links_with(branch([PAIR_ROOT_K, 0.0, 0.0, 0.0, 0.0, 0.0], PAIR_TAU_S))
}

fn pair_links_with(memory: SpatialRelaxationBranch) -> Vec<RelaxationNetworkLink> {
    vec![RelaxationNetworkLink {
        viscous: ViscousNetworkLink {
            bodies: [0, 1],
            arms_body_m: [Vec3::ZERO; 2],
            damping: damping([PAIR_C, 0.0, 0.0, 0.0, 0.0, 0.0]),
        },
        branches: vec![memory],
    }]
}

fn pair() -> RelaxationNetwork {
    RelaxationNetwork::new(pair_properties(), pair_links(), Gravity::ZERO, controls())
        .expect("network")
}

/// The pair with its branch lossless: an equilibrium channel, which stores
/// and returns energy and never heats.
fn elastic_pair() -> RelaxationNetwork {
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = PAIR_ROOT_K;
    let lossless = SpatialRelaxationBranch::from_equilibrium_gram_rows(rows).expect("branch");
    RelaxationNetwork::new(pair_properties(), pair_links_with(lossless), Gravity::ZERO, controls())
        .expect("network")
}

/// The pair with relative speed `u_m_s` (body 1 minus body 0) and memory
/// amplitude `q` [sqrt(J)].
fn pair_state(u_m_s: f64, q: f64) -> RelaxationState {
    let p = pair_properties();
    RelaxationState {
        bodies: vec![
            body(p[0], Vec3::new(0.0, 0.0, 0.0), Vec3::new(-0.5 * u_m_s, 0.0, 0.0), Vec3::ZERO),
            body(p[1], Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.5 * u_m_s, 0.0, 0.0), Vec3::ZERO),
        ],
        memory: vec![vec![[q, 0.0, 0.0, 0.0, 0.0, 0.0]]],
    }
}

fn relative_speed(state: &RelaxationState) -> f64 {
    let p = pair_properties();
    let speed = |body: usize| {
        state.bodies[body]
            .center_of_mass_velocity_world(p[body])
            .expect("velocity")
            .x
    };
    speed(1) - speed(0)
}

/// The pair's own equations, integrated independently (RK4, 16384 steps).
/// Two 2 kg masses have relative mobility 1 / kg, so with u the relative
/// speed and q the memory amplitude: `u' = -c u - k q`, `q' = k u - q / tau`
/// and heat `c u^2 + q^2 / tau`. Returns `[u, q, heat]` at `t_s`.
fn pair_reference(u: f64, q: f64, t_s: f64) -> [f64; 3] {
    pair_reference_at(1.0 / PAIR_TAU_S, u, q, t_s)
}

/// [`pair_reference`] with the memory relaxing at `rate_per_s` (zero for an
/// equilibrium channel).
fn pair_reference_at(rate_per_s: f64, u: f64, q: f64, t_s: f64) -> [f64; 3] {
    let f = |s: [f64; 3]| {
        [
            -PAIR_C * s[0] - PAIR_ROOT_K * s[1],
            PAIR_ROOT_K * s[0] - rate_per_s * s[1],
            PAIR_C * s[0] * s[0] + rate_per_s * s[1] * s[1],
        ]
    };
    let add = |a: [f64; 3], b: [f64; 3], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [u, q, 0.0];
    let h = t_s / 16384.0;
    for _ in 0..16384 {
        let a = f(state);
        let b = f(add(state, a, 0.5 * h));
        let c = f(add(state, b, 0.5 * h));
        let d = f(add(state, c, h));
        for i in 0..3 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

fn network(gravity: Gravity) -> RelaxationNetwork {
    RelaxationNetwork::new(properties(), links(), gravity, controls()).expect("network")
}

fn owner_step(
    state: &RelaxationState,
    wrenches: &[Wrench],
    gravity: Gravity,
    dt_s: f64,
) -> RelaxationNetworkMotionStep {
    let idle = [HeldInternalEffort::MaterialFirst([0.0; 6]); 2];
    owner_step_driven(state, wrenches, &idle, gravity, dt_s)
}

/// The owner's step of the three-body network with `efforts` held on its
/// two connectors.
fn owner_step_driven(
    state: &RelaxationState,
    wrenches: &[Wrench],
    efforts: &[HeldInternalEffort; 2],
    gravity: Gravity,
    dt_s: f64,
) -> RelaxationNetworkMotionStep {
    step_relaxation_network_with_loads(
        &state.bodies,
        &properties(),
        &links(),
        &state.memory,
        efforts,
        HeldExternalLoads { wrenches, gravity },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step")
}

/// A state of the three-body network other than `start` (a tenth of a
/// second later, under `efforts`), for the conformance check to hand in as
/// the stale content of the state a propose must overwrite.
fn stale_after(
    start: &RelaxationState,
    efforts: &[HeldInternalEffort; 2],
    gravity: Gravity,
) -> RelaxationState {
    let step = owner_step_driven(start, &held_wrenches(), efforts, gravity, 0.1);
    RelaxationState {
        bodies: step.states_after().to_vec(),
        memory: step.memory_after().clone(),
    }
}

fn force(x: f64, y: f64, z: f64) -> Wrench {
    Wrench {
        force_world: Vec3::new(x, y, z),
        torque_body: Vec3::ZERO,
    }
}

/// One step of `component` from `state`, outside the kernel.
fn propose(
    component: &RelaxationNetwork,
    state: &RelaxationState,
    dt_s: f64,
    inputs: &[f64],
    gate: &CancelGate,
) -> Result<(RelaxationState, Receipt, Vec<f64>), RelaxationNetworkRefusal> {
    let mut next = state.clone();
    let width: usize = component
        .ports()
        .iter()
        .map(fs_ports::component::PortSpec::output_width)
        .sum();
    let mut outputs = vec![f64::NAN; width];
    let proposal = with_cx(gate, |cx| {
        component.propose(state, dt_s, inputs, &mut next, &mut outputs, cx)
    })?;
    let Proposal::Step(receipt) = proposal else {
        panic!("the network takes the step or refuses");
    };
    Ok((next, receipt, outputs))
}

const HELD: [f64; 9] = [0.3, 0.0, -0.1, 0.0, 0.0, 0.0, -0.2, 0.4, 0.0];

fn held_wrenches() -> [Wrench; 3] {
    [force(0.3, 0.0, -0.1), force(0.0, 0.0, 0.0), force(-0.2, 0.4, 0.0)]
}

#[test]
fn the_relaxation_network_conforms_and_books_the_owners_work_heat_and_memory() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let component = network(gravity);
    let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(names, ["force0", "force1", "force2", "heat0", "heat1"]);
    assert_eq!((component.bodies(), component.connectors()), (3, 2));
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let stale = {
        let step = owner_step(&start, &held_wrenches(), gravity, 0.1);
        RelaxationState {
            bodies: step.states_after().to_vec(),
            memory: step.memory_after().clone(),
        }
    };
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0]);
    for dt_s in [4.0e-3, 2.0e-2] {
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
        // The owner's own accounts: the kernel books applied work in and
        // connector heat out; the stored change includes the memory.
        let native = owner_step(&start, &held_wrenches(), gravity, dt_s);
        let applied = native.external_loads.applied_work_j;
        let heat = native.dissipated_energy_j;
        let size = applied.abs() + heat + native.mechanical_energy_change_j.abs();
        assert!((audit.port_energy_in_j - (applied - heat)).abs() <= 1.0e-11 * size);
        assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
        assert!(native.memory_energy_change_j.abs() > 1.0e-6, "the memory moved");
        // The idle actuator does exactly no work, and the owner's residual
        // is the receipt's numerical loss, with its sign.
        assert_eq!(native.actuator_work_j, 0.0);
        assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
        assert!(native.work_heat_residual_j != 0.0, "a residual to carry");
        assert!((receipt.dissipated_j - heat).abs() <= 1.0e-14 * heat);
        // Entropy where the heat is released: 300 K and 320 K.
        let entropy =
            native.connector_dissipation_j[0] / 300.0 + native.connector_dissipation_j[1] / 320.0;
        assert!((receipt.entropy_generated_j_per_k - entropy).abs() <= 1.0e-14 * entropy);
        for connector in 0..2 {
            assert_eq!(outputs[9 + connector], -native.connector_dissipation_j[connector] / dt_s);
        }
    }
}

#[test]
fn the_wrapper_adds_nothing_to_the_owners_step() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let component = network(gravity);
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 300.0]);
    let dt_s = 4.0e-3;
    let (next, _, outputs) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    let native = owner_step(&start, &held_wrenches(), gravity, dt_s);
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
    assert_ne!(next.memory, start.memory);
    for (index, state) in next.bodies.iter().enumerate() {
        let (a, b) = (
            start.bodies[index].pose().position_world(),
            state.pose().position_world(),
        );
        assert_eq!(outputs[3 * index], (b.x - a.x) / dt_s);
        assert_eq!(outputs[3 * index + 1], (b.y - a.y) / dt_s);
        assert_eq!(outputs[3 * index + 2], (b.z - a.z) / dt_s);
    }
}

#[test]
fn memory_is_stored_energy_and_its_release_is_motion_and_heat() {
    let gate = CancelGate::new_clock_free();
    // The three-body network at rest stores exactly its preloaded memory,
    // and nothing when relaxed.
    let component = network(Gravity::ZERO);
    let p = properties();
    let rest = vec![
        body(p[0], Vec3::new(-1.0, 0.0, 0.2), Vec3::ZERO, Vec3::ZERO),
        body(p[1], Vec3::new(0.0, 0.1, 0.0), Vec3::ZERO, Vec3::ZERO),
        body(p[2], Vec3::new(1.0, -0.1, 0.0), Vec3::ZERO, Vec3::ZERO),
    ];
    let loaded = RelaxationState {
        bodies: rest.clone(),
        memory: preloaded(),
    };
    let stored = component.stored_energy_j(&loaded).expect("energy");
    assert!((stored - PRELOAD_J).abs() <= 1.0e-15, "{stored} vs {PRELOAD_J}");
    let relaxed = component.relaxed(rest);
    assert_eq!(relaxed.memory, vec![vec![[0.0; 6]; 2], vec![[0.0; 6]; 1]]);
    assert_eq!(component.stored_energy_j(&relaxed).expect("energy"), 0.0);

    // The pair, at rest with a preloaded memory of 0.6 sqrt(J) (0.18 J),
    // nothing applied: the memory sets the bodies moving and heats the
    // connector. Twenty steps to a quarter second, against the pair's own
    // equations integrated independently.
    let pair = pair();
    let (q0, dt_s, steps) = (0.6, 0.0125, 20_u32);
    let mut state = pair_state(0.0, q0);
    assert!((pair.stored_energy_j(&state).expect("energy") - 0.5 * q0 * q0).abs() <= 1.0e-15);
    let inputs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0];
    let mut heat_j = 0.0_f64;
    for _ in 0..steps {
        let (next, receipt, outputs) = propose(&pair, &state, dt_s, &inputs, &gate).expect("step");
        heat_j -= dt_s * outputs[6];
        assert!(receipt.numerical_loss_j.abs() <= 1.0e-10);
        state = next;
    }
    let [u, q, heat] = pair_reference(0.0, q0, dt_s * f64::from(steps));
    assert!(u.abs() > 0.1, "elastic recovery sets the pair in motion: {u} m/s");
    assert!((relative_speed(&state) - u).abs() <= 1.0e-9, "{} vs {u}", relative_speed(&state));
    assert!((state.memory[0][0][0] - q).abs() <= 1.0e-9, "{} vs {q}", state.memory[0][0][0]);
    assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");
    // Storage: relative kinetic energy u^2 / 2 and memory q^2 / 2; what left
    // it is the heat.
    let left = pair.stored_energy_j(&state).expect("energy");
    assert!((left - 0.5 * (u * u + q * q)).abs() <= 1.0e-9);
    assert!((0.5 * q0 * q0 - left - heat_j).abs() <= 1.0e-10);
}

#[test]
fn an_equilibrium_branch_stores_and_returns_energy_and_never_heats() {
    let gate = CancelGate::new_clock_free();
    // The pair with a lossless branch: only the Kelvin viscosity heats.
    let elastic = elastic_pair();
    let (u0, q0, dt_s, steps) = (-2.0, 0.6, 0.0125, 20_u32);
    let mut state = pair_state(u0, q0);
    let stored0 = elastic.stored_energy_j(&state).expect("energy");
    assert!((stored0 - 0.5 * (u0 * u0 + q0 * q0)).abs() <= 1.0e-14);
    let inputs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0];
    let mut heat_j = 0.0_f64;
    for _ in 0..steps {
        let (next, _, outputs) = propose(&elastic, &state, dt_s, &inputs, &gate).expect("step");
        heat_j -= dt_s * outputs[6];
        state = next;
    }
    let [u, q, heat] = pair_reference_at(0.0, u0, q0, dt_s * f64::from(steps));
    assert!((relative_speed(&state) - u).abs() <= 1.0e-9, "{} vs {u}", relative_speed(&state));
    assert!((state.memory[0][0][0] - q).abs() <= 1.0e-9, "{} vs {q}", state.memory[0][0][0]);
    // The heat is the Kelvin viscosity's alone; with a Maxwell branch of the
    // same factor the same run heats more.
    assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");
    let [_, _, maxwell_heat] = pair_reference(u0, q0, dt_s * f64::from(steps));
    assert!(maxwell_heat > heat + 0.05, "{maxwell_heat} vs {heat}");
    // The closing pair unloaded the branch (its amplitude ran from 0.6 down
    // through zero), and what storage lost is the Kelvin heat alone.
    assert!((q - q0).abs() > 0.5, "the memory moved: {q}");
    let stored = elastic.stored_energy_j(&state).expect("energy");
    assert!((stored - 0.5 * (u * u + q * q)).abs() <= 1.0e-9);
    assert!((stored0 - stored - heat_j).abs() <= 1.0e-10);
}

#[test]
fn without_branches_the_network_is_the_viscous_network() {
    let gate = CancelGate::new_clock_free();
    let plain: Vec<RelaxationNetworkLink> = viscous_links()
        .into_iter()
        .map(|viscous| RelaxationNetworkLink {
            viscous,
            branches: Vec::new(),
        })
        .collect();
    let memoryless =
        RelaxationNetwork::new(properties(), plain, Gravity::ZERO, controls()).expect("network");
    let viscous = ViscousNetwork::new(
        properties(),
        viscous_links(),
        Gravity::ZERO,
        ViscousNetworkMotionControls {
            source: source(),
            source_substeps: 1,
            limits: ViscousNetworkLimits {
                maximum_bodies: 32,
                maximum_connectors: 64,
                maximum_pair_relaxations: 100_000,
                maximum_live_bytes: 1_000_000,
            },
        },
    )
    .expect("network");
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 300.0]);
    let dt_s = 4.0e-3;
    let start = memoryless.relaxed(moving());
    assert_eq!(start.memory, vec![Vec::<[f64; 6]>::new(); 2]);
    let (_, _, ours) = propose(&memoryless, &start, dt_s, &inputs, &gate).expect("step");
    let mut next = Vec::new();
    let mut theirs = vec![f64::NAN; 11];
    with_cx(&gate, |cx| viscous.propose(&moving(), dt_s, &inputs, &mut next, &mut theirs, cx))
        .expect("step");
    // Two owners of the same connectors by two methods (one simultaneous
    // solve of the whole graph here, sweeps over the pairs there): the same
    // velocities and the same heat to the difference of the methods at this
    // step (measured: 4e-7 of the first connector's heat rate).
    for (port, (a, b)) in ours.iter().zip(&theirs).enumerate() {
        let size = a.abs().max(b.abs()).max(1.0e-3);
        assert!((a - b).abs() <= 1.0e-5 * size, "port output {port}: {a} vs {b}");
    }
    assert!(ours[9] < -1.0e-3 && ours[10] < -1.0e-3, "both connectors heat");
}

#[test]
fn the_relaxation_network_refuses_what_it_cannot_represent() {
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        RelaxationNetwork::new(Vec::new(), Vec::new(), Gravity::ZERO, controls()).unwrap_err(),
        RelaxationNetworkBuildError::NoBodies
    );
    let mut wrong = links();
    wrong[1].viscous.bodies = [1, 5];
    assert_eq!(
        RelaxationNetwork::new(properties(), wrong, Gravity::ZERO, controls()).unwrap_err(),
        RelaxationNetworkBuildError::Link { connector: 1 }
    );
    let mut same = links();
    same[0].viscous.bodies = [2, 2];
    assert_eq!(
        RelaxationNetwork::new(properties(), same, Gravity::ZERO, controls()).unwrap_err(),
        RelaxationNetworkBuildError::Link { connector: 0 }
    );

    let component = network(Gravity::ZERO);
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let mut inputs = vec![0.0; 9];
    inputs.extend([300.0, 300.0]);
    let step = |state: &RelaxationState, inputs: &[f64], gate: &CancelGate| {
        propose(&component, state, 4.0e-3, inputs, gate).map(|_| ())
    };
    // A memory that is not one amplitude per branch.
    let mut short = start.clone();
    short.memory[0].pop();
    assert_eq!(
        step(&short, &inputs, &gate).unwrap_err(),
        RelaxationNetworkRefusal::Memory { connector: 0 }
    );
    assert_eq!(
        component.stored_energy_j(&short).unwrap_err(),
        RelaxationNetworkRefusal::Memory { connector: 0 }
    );
    let mut missing = start.clone();
    missing.memory.pop();
    assert_eq!(
        step(&missing, &inputs, &gate).unwrap_err(),
        RelaxationNetworkRefusal::Memory { connector: 1 }
    );
    let mut wild = start.clone();
    wild.memory[1][0][3] = f64::NAN;
    assert_eq!(
        step(&wild, &inputs, &gate).unwrap_err(),
        RelaxationNetworkRefusal::Memory { connector: 1 }
    );
    let mut fewer = start.clone();
    fewer.bodies.pop();
    assert_eq!(
        step(&fewer, &inputs, &gate).unwrap_err(),
        RelaxationNetworkRefusal::Bodies {
            expected: 3,
            found: 2
        }
    );
    let mut longer = start.clone();
    longer.memory.push(Vec::new());
    assert_eq!(
        step(&longer, &inputs, &gate).unwrap_err(),
        RelaxationNetworkRefusal::Memory { connector: 2 }
    );
    // A heat port's temperature that is no temperature.
    for kelvin in [0.0, -20.0, f64::INFINITY] {
        let mut cold = inputs.clone();
        cold[10] = kelvin;
        assert_eq!(
            step(&start, &cold, &gate).unwrap_err(),
            RelaxationNetworkRefusal::Temperature {
                connector: 1,
                kelvin
            }
        );
    }
    let mut unknown = inputs.clone();
    unknown[9] = f64::NAN;
    assert!(matches!(
        step(&start, &unknown, &gate),
        Err(RelaxationNetworkRefusal::Temperature { connector: 0, kelvin }) if kelvin.is_nan()
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        step(&start, &inputs, &cancelled).unwrap_err(),
        RelaxationNetworkRefusal::Cancelled
    );
}

#[test]
fn connector_heat_with_memory_reaches_thermal_bodies_and_the_ledger_closes() {
    let gate = CancelGate::new_clock_free();
    // The pair closing at 2 m/s with a preloaded memory; its connector's
    // heat goes 70 % to one thermal body, 20 % to another and the rest to
    // a named sink.
    let start = pair_state(-2.0, 0.6);
    let mut builder = SystemBuilder::new();
    let moving_node = builder.add("network", pair(), start.clone()).expect("network");
    let cell_a = builder
        .add("cell_a", HeatCapacity::new(5.0).expect("capacity"), 300.0_f64)
        .expect("cell a");
    let cell_b = builder
        .add("cell_b", HeatCapacity::new(8.0).expect("capacity"), 300.0_f64)
        .expect("cell b");
    for index in 0..2 {
        let port = builder.port(moving_node, &format!("force{index}")).expect("port");
        builder
            .boundary(port, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .expect("free");
    }
    builder
        .heat_route(
            builder.port(moving_node, "heat0").expect("port"),
            &[
                (builder.port(cell_a, "heat").expect("port"), 0.7),
                (builder.port(cell_b, "heat").expect("port"), 0.2),
            ],
            Some(heat_sink("unallocated", 293.15).expect("sink")),
            RouteTemperature::SinkEquivalent,
        )
        .expect("route");
    let mut system = builder.seal().expect("sealed");
    let heat_port = system.port(moving_node, "heat0").expect("heat port");
    assert!(system.guess(heat_port, &[300.0]));

    // The owner marched alone over the same steps.
    let idle = [Wrench::ZERO; 2];
    let none = [HeldInternalEffort::MaterialFirst([0.0; 6]); 1];
    let mut twin = start.clone();
    let (dt_s, steps) = (0.0125, 20_u32);
    let mut heat_j = 0.0_f64;
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            let report = system.step(dt_s, cx).expect("kernel step");
            assert_eq!(report.dt_s, dt_s);
            let native = step_relaxation_network_with_loads(
                &twin.bodies,
                &pair_properties(),
                &pair_links(),
                &twin.memory,
                &none,
                HeldExternalLoads {
                    wrenches: &idle,
                    gravity: Gravity::ZERO,
                },
                dt_s,
                controls(),
                || false,
            )
            .expect("the owner takes the step");
            heat_j += native.connector_dissipation_j[0];
            twin = RelaxationState {
                bodies: native.states_after().to_vec(),
                memory: native.memory_after().clone(),
            };
        }
    });
    // The motion and the memory do not depend on where the heat goes.
    let reached = system.state::<RelaxationState>(moving_node).expect("state");
    assert_eq!(reached, &twin);
    // And they are the pair's own equations.
    let [u, q, heat] = pair_reference(-2.0, 0.6, dt_s * f64::from(steps));
    assert!((relative_speed(reached) - u).abs() <= 1.0e-9);
    assert!((reached.memory[0][0][0] - q).abs() <= 1.0e-9);
    assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");

    // Every joule: what left kinetic and memory storage is in the two
    // thermal bodies or in the named sink.
    assert!(heat_j > 0.3, "enough heat to measure: {heat_j} J");
    let ledger = system.ledger();
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
    let warmed_a = account("cell_a").stored_change_j.value();
    let warmed_b = account("cell_b").stored_change_j.value();
    assert!((warmed_a - 0.7 * heat_j).abs() <= 1.0e-11 * heat_j, "{warmed_a}");
    assert!((warmed_b - 0.2 * heat_j).abs() <= 1.0e-11 * heat_j, "{warmed_b}");
    assert!((ledger.sink_heat_j() - 0.1 * heat_j).abs() <= 1.0e-11 * heat_j);
    let released = -account("network").stored_change_j.value();
    assert!((released - heat_j).abs() <= 1.0e-10, "{released} vs {heat_j}");
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10);
}

/// The three-body network with an actuator on both connectors: the first
/// held in the material axes of body 0, the second in world axes.
fn driven_network(gravity: Gravity) -> RelaxationNetwork {
    network(gravity)
        .driven(0, DriveFrame::MaterialFirst)
        .expect("actuator on connector 0")
        .driven(1, DriveFrame::World)
        .expect("actuator on connector 1")
}

/// Held efforts, force then couple: connector 0 in the material axes of
/// body 0, connector 1 in world axes.
const DRIVE: [f64; 12] = [
    0.5, -0.2, 0.1, 0.03, 0.0, -0.04, -0.3, 0.25, 0.0, 0.0, 0.02, 0.05,
];

fn drive_efforts() -> [HeldInternalEffort; 2] {
    [
        HeldInternalEffort::MaterialFirst(core::array::from_fn(|c| DRIVE[c])),
        HeldInternalEffort::World(core::array::from_fn(|c| DRIVE[6 + c])),
    ]
}

fn turn(state: RigidBodyState, orientation: UnitQuaternion) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(state.pose().position_world(), orientation).expect("pose"),
        state.linear_momentum_world(),
        state.angular_momentum_body(),
    )
    .expect("state")
}

/// [`moving`] with the first two bodies turned well away from the world
/// axes, so that material and world components differ.
fn turned() -> Vec<RigidBodyState> {
    let orientations = [
        UnitQuaternion::new(0.9, -0.1, 0.3, 0.2).expect("orientation"),
        UnitQuaternion::new(0.8, 0.2, -0.3, 0.4).expect("orientation"),
        UnitQuaternion::IDENTITY,
    ];
    moving()
        .into_iter()
        .zip(orientations)
        .map(|(state, orientation)| turn(state, orientation))
        .collect()
}

/// Energy the kernel books on the two ports of the actuator in `slot` of a
/// network with `bodies` bodies and `connectors` connectors, and the sum of
/// the sizes of its six products.
fn drive_port_energy_j(
    (bodies, connectors, slot): (usize, usize, usize),
    dt_s: f64,
    inputs: &[f64],
    outputs: &[f64],
) -> (f64, f64) {
    let base = 3 * bodies + connectors + 6 * slot;
    drive_energy_at((base, base), dt_s, inputs, outputs)
}

/// The same for an actuator whose six efforts start at `input` in the
/// input vector and whose six rates start at `output` in the output vector.
fn drive_energy_at(
    (input, output): (usize, usize),
    dt_s: f64,
    inputs: &[f64],
    outputs: &[f64],
) -> (f64, f64) {
    let (mut energy, mut size) = (0.0_f64, 0.0_f64);
    for c in 0..6 {
        energy += inputs[input + c] * outputs[output + c];
        size += (inputs[input + c] * outputs[output + c]).abs();
    }
    (dt_s * energy, dt_s * size)
}

#[test]
fn a_driven_network_conforms_and_its_drive_ports_carry_the_owners_actuator_work() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let component = driven_network(gravity);
    let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(
        names,
        ["force0", "force1", "force2", "heat0", "heat1", "drive0", "couple0", "drive1", "couple1"]
    );
    // A held force is conjugate to a velocity, a held couple to an angular
    // velocity; each port is three components wide.
    for (index, port) in component.ports()[5..].iter().enumerate() {
        let kind = if index % 2 == 0 {
            PortKind::MechanicalForceVelocity
        } else {
            PortKind::RotationalTorqueAngularVelocity
        };
        assert_eq!(port.carrier(), Carrier::Power(kind), "{}", port.id().as_str());
        assert_eq!(port.width(), 3);
    }
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let efforts = drive_efforts();
    let stale = stale_after(&start, &efforts, gravity);
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0]);
    inputs.extend(DRIVE);
    for dt_s in [4.0e-3, 2.0e-2] {
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
        let native = owner_step_driven(&start, &held_wrenches(), &efforts, gravity, dt_s);
        // The kernel books applied work and actuator work in and connector
        // heat out.
        let applied = native.external_loads.applied_work_j;
        let heat = native.dissipated_energy_j;
        let supplied = native.actuator_work_j;
        let size = applied.abs() + heat + supplied.abs() + native.mechanical_energy_change_j.abs();
        assert!(
            (audit.port_energy_in_j - (applied + supplied - heat)).abs() <= 1.0e-11 * size,
            "{} vs {}",
            audit.port_energy_in_j,
            applied + supplied - heat
        );
        // Each actuator's two ports carry the owner's work on that
        // connector, in material axes and in world axes alike.
        for connector in 0..2 {
            let (port_j, terms_j) = drive_port_energy_j((3, 2, connector), dt_s, &inputs, &outputs);
            let native_j = native.connector_actuator_work_j[connector];
            assert!(native_j.abs() > 1.0e-5, "connector {connector} does work: {native_j} J");
            assert!(
                (port_j - native_j).abs() <= 1.0e-13 * terms_j,
                "connector {connector}: ports {port_j:e} J, owner {native_j:e} J"
            );
        }
        // With every connector driven nothing is left over: the numerical
        // loss is the owner's residual alone.
        // The allowance is rounding: 64 ulp of the terms (about 13 J of
        // kinetic and potential energy here), measured 1.9e-13 J.
        assert!(
            receipt.energy_tolerance_j <= 1.0e-12,
            "{:e} J",
            receipt.energy_tolerance_j
        );
        assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
        assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
        assert!((receipt.dissipated_j - heat).abs() <= 1.0e-14 * heat);
        for connector in 0..2 {
            assert_eq!(outputs[9 + connector], -native.connector_dissipation_j[connector] / dt_s);
        }
        // The wrapper adds nothing to the driven step either.
        let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
        assert_eq!(next.bodies.as_slice(), native.states_after());
        assert_eq!(&next.memory, native.memory_after());
    }
}

#[test]
fn one_actuator_on_a_network_of_two_connectors_carries_all_the_actuator_work() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let efforts = drive_efforts();
    let stale = stale_after(&start, &efforts, gravity);
    // An actuator on connector 1 only: connector 0 holds no effort.
    let single = network(gravity).driven(1, DriveFrame::World).expect("actuator");
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0]);
    inputs.extend(&DRIVE[6..]);
    let only = [HeldInternalEffort::MaterialFirst([0.0; 6]), efforts[1]];
    let dt_s = 2.0e-2;
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&single, &case, cx)).expect("conforms");
    let Observed::Stepped {
        receipt, outputs, ..
    } = observed
    else {
        panic!("the network takes the step");
    };
    let native = owner_step_driven(&start, &held_wrenches(), &only, gravity, dt_s);
    // A connector without an actuator holds a zero effort, and the owner's
    // work for it is exactly zero: nothing is left outside the ports.
    assert_eq!(native.connector_actuator_work_j[0], 0.0);
    assert!(native.actuator_work_j.abs() > 1.0e-5, "{} J", native.actuator_work_j);
    let (port_j, terms_j) = drive_port_energy_j((3, 2, 0), dt_s, &inputs, &outputs);
    assert!((port_j - native.actuator_work_j).abs() <= 1.0e-13 * terms_j);
    assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
    let (next, _, _) = propose(&single, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
}

/// The driven pair's own equations (RK4, 16384 steps): with a held axial
/// force `f` on the connector, `u' = f - c u - k q`, `q' = k u - q / tau`,
/// heat `c u^2 + q^2 / tau` and actuator work `f u`. Returns
/// `[u, q, heat, work]` at `t_s`.
fn driven_pair_reference(force: f64, u: f64, q: f64, t_s: f64) -> [f64; 4] {
    let f = |s: [f64; 4]| {
        [
            force - PAIR_C * s[0] - PAIR_ROOT_K * s[1],
            PAIR_ROOT_K * s[0] - s[1] / PAIR_TAU_S,
            PAIR_C * s[0] * s[0] + s[1] * s[1] / PAIR_TAU_S,
            force * s[0],
        ]
    };
    let add = |a: [f64; 4], b: [f64; 4], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [u, q, 0.0, 0.0];
    let h = t_s / 16384.0;
    for _ in 0..16384 {
        let a = f(state);
        let b = f(add(state, a, 0.5 * h));
        let c = f(add(state, b, 0.5 * h));
        let d = f(add(state, c, h));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

/// `steps` steps of the pair with a held axial force on its connector, in
/// `frame`. Returns the reached state, the heat and the energy the drive
/// ports carried in.
fn drive_pair(
    frame: DriveFrame,
    start: &RelaxationState,
    force: f64,
    (dt_s, steps): (f64, u32),
    gate: &CancelGate,
) -> (RelaxationState, f64, f64) {
    let pair = pair().driven(0, frame).expect("actuator");
    let mut inputs = vec![0.0; 6];
    inputs.push(300.0);
    inputs.extend([force, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let mut state = start.clone();
    let (mut heat_j, mut work_j) = (0.0_f64, 0.0_f64);
    for _ in 0..steps {
        let (next, receipt, outputs) = propose(&pair, &state, dt_s, &inputs, gate).expect("step");
        heat_j -= dt_s * outputs[6];
        work_j += drive_port_energy_j((2, 1, 0), dt_s, &inputs, &outputs).0;
        assert!(receipt.numerical_loss_j.abs() <= 1.0e-10);
        state = next;
    }
    (state, heat_j, work_j)
}

#[test]
fn a_held_force_on_the_pair_does_the_work_of_its_own_equations() {
    let gate = CancelGate::new_clock_free();
    let pair = pair();
    let run = (0.0125, 20_u32);
    let t_s = run.0 * f64::from(run.1);
    // Pushed apart from rest with a relaxed memory, then driven against its
    // motion (the actuator recovers energy while the connector still heats).
    for (force, u0, q0) in [(3.0, 0.0, 0.0), (-0.7, 2.0, 0.6)] {
        let start = pair_state(u0, q0);
        let (state, heat_j, work_j) = drive_pair(DriveFrame::MaterialFirst, &start, force, run, &gate);
        let [u, q, heat, work] = driven_pair_reference(force, u0, q0, t_s);
        assert!((relative_speed(&state) - u).abs() <= 1.0e-9, "{} vs {u}", relative_speed(&state));
        assert!((state.memory[0][0][0] - q).abs() <= 1.0e-9, "{} vs {q}", state.memory[0][0][0]);
        assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");
        assert!((work_j - work).abs() <= 1.0e-9, "{work_j} vs {work}");
        assert!(work.abs() > 0.05 && heat > 0.01, "work {work} J and heat {heat} J to measure");
        // First law over the run: the work that came in is stored or is heat.
        let stored_j = pair.stored_energy_j(&state).expect("energy")
            - pair.stored_energy_j(&start).expect("energy");
        assert!((work_j - heat_j - stored_j).abs() <= 1.0e-10);
        if force < 0.0 {
            // The actuator takes energy out (-0.30 J) while the connector
            // still heats (0.67 J).
            assert!(
                work_j < -0.25 && heat_j > 0.5,
                "recovery {work_j} J, connector heat {heat_j} J"
            );
        }
        // Nothing turns in this pair, so world axes are material axes.
        let (world, world_heat_j, world_work_j) =
            drive_pair(DriveFrame::World, &start, force, run, &gate);
        assert_eq!(world, state);
        assert_eq!((world_heat_j, world_work_j), (heat_j, work_j));
    }
}

/// The pair's relative rate in world axes, from a state alone: the
/// velocity of the second centre of mass against the point of the first
/// body that is there, `v1 - v0 - w0 x (x1 - x0)`, and the relative angular
/// velocity `w1 - w0`.
fn pair_relative_rate(state: &RelaxationState) -> [f64; 6] {
    let p = pair_properties();
    let velocity = |body: usize| {
        state.bodies[body]
            .center_of_mass_velocity_world(p[body])
            .expect("velocity")
    };
    let spin = |body: usize| {
        let s = state.bodies[body];
        s.pose().orientation().rotate_body_to_world(
            p[body]
                .angular_velocity_body_checked(s.angular_momentum_body())
                .expect("angular velocity"),
        )
    };
    let (v0, v1, w0, w1) = (velocity(0), velocity(1), spin(0), spin(1));
    let (x0, x1) = (
        state.bodies[0].pose().position_world(),
        state.bodies[1].pose().position_world(),
    );
    let r = [x1.x - x0.x, x1.y - x0.y, x1.z - x0.z];
    let carried = [
        w0.y * r[2] - w0.z * r[1],
        w0.z * r[0] - w0.x * r[2],
        w0.x * r[1] - w0.y * r[0],
    ];
    [
        v1.x - v0.x - carried[0],
        v1.y - v0.y - carried[1],
        v1.z - v0.z - carried[2],
        w1.x - w0.x,
        w1.y - w0.y,
        w1.z - w0.z,
    ]
}

#[test]
fn a_world_held_drive_returns_the_relative_velocity_in_world_axes() {
    let gate = CancelGate::new_clock_free();
    // The pair with its first body turned about one radian from the world
    // axes and spinning, its second body spinning otherwise, a metre apart.
    // The connector's attachments are the two centres of mass, and the
    // first body's wrench is carried to the second attachment; so the rate
    // conjugate to a force is the velocity of the second centre against
    // the point of the first body that is there, and the rate conjugate to
    // a couple is the relative angular velocity.
    let p = pair_properties();
    let start = RelaxationState {
        bodies: vec![
            turn(
                body(p[0], Vec3::ZERO, Vec3::new(-0.6, 0.2, 0.1), Vec3::new(0.4, -0.6, 0.9)),
                UnitQuaternion::new(0.85, 0.3, -0.2, 0.4).expect("orientation"),
            ),
            body(p[1], Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.5, -0.1, 0.3), Vec3::new(-0.2, 0.5, 0.3)),
        ],
        memory: vec![vec![[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]]],
    };
    let held = [0.4, -0.3, 0.2, 0.05, 0.1, -0.08];
    let mut inputs = vec![0.0; 6];
    inputs.push(300.0);
    inputs.extend(held);
    // The mean of the start and end relative rates is the conjugate rate
    // to second order in the step: halving the step quarters the
    // difference.
    let gap = |frame: DriveFrame, dt_s: f64| {
        let pair = pair().driven(0, frame).expect("actuator");
        let (next, _, outputs) = propose(&pair, &start, dt_s, &inputs, &gate).expect("step");
        let (before, after) = (pair_relative_rate(&start), pair_relative_rate(&next));
        let mut rate = [0.0; 6];
        rate.copy_from_slice(&outputs[7..13]);
        if frame == DriveFrame::MaterialFirst {
            // The same vector in the first body's axes, which turn during
            // the step: turned to world axes by the start orientation and
            // by the end orientation, the mean of the two is the world
            // rate to second order.
            let material = rate;
            rate = [0.0; 6];
            for state in [&start, &next] {
                let rotation = state.bodies[0].pose().orientation();
                let linear =
                    rotation.rotate_body_to_world(Vec3::new(material[0], material[1], material[2]));
                let angular =
                    rotation.rotate_body_to_world(Vec3::new(material[3], material[4], material[5]));
                for (c, value) in [linear.x, linear.y, linear.z, angular.x, angular.y, angular.z]
                    .into_iter()
                    .enumerate()
                {
                    rate[c] += 0.5 * value;
                }
            }
        }
        let mut worst = 0.0_f64;
        for c in 0..6 {
            let mean = 0.5 * (before[c] + after[c]);
            worst = worst.max((rate[c] - mean).abs());
        }
        worst
    };
    for frame in [DriveFrame::World, DriveFrame::MaterialFirst] {
        // Measured: 2.0e-5 then 5.0e-6 in world axes, 3.2e-5 then 8.2e-6
        // in material axes, of a relative rate of about 1 m/s.
        let (coarse, fine) = (gap(frame, 2.0e-2), gap(frame, 1.0e-2));
        assert!(coarse <= 1.0e-4, "{frame:?}: rate against the relative velocity: {coarse:e}");
        assert!(
            (3.0..=5.0).contains(&(coarse / fine)),
            "{frame:?}: second order: {coarse:e} then {fine:e}"
        );
    }
}

#[test]
fn an_actuator_is_refused_where_it_cannot_be() {
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        network(Gravity::ZERO).driven(2, DriveFrame::World).unwrap_err(),
        RelaxationNetworkBuildError::Drive { connector: 2 }
    );
    assert_eq!(
        network(Gravity::ZERO)
            .driven(1, DriveFrame::World)
            .expect("actuator")
            .driven(1, DriveFrame::MaterialFirst)
            .unwrap_err(),
        RelaxationNetworkBuildError::Drive { connector: 1 }
    );
    let component = driven_network(Gravity::ZERO);
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let mut inputs = vec![0.0; 9];
    inputs.extend([300.0, 300.0]);
    inputs.extend(DRIVE);
    assert!(propose(&component, &start, 4.0e-3, &inputs, &gate).is_ok());
    // A held effort that is not a number, on the force port and on the
    // couple port.
    for (index, connector) in [(11, 0), (16, 0), (17, 1), (22, 1)] {
        for value in [f64::NAN, f64::INFINITY] {
            let mut wild = inputs.clone();
            wild[index] = value;
            assert_eq!(
                propose(&component, &start, 4.0e-3, &wild, &gate).unwrap_err(),
                RelaxationNetworkRefusal::Effort { connector }
            );
        }
    }
}

#[test]
fn a_held_drive_is_a_boundary_and_its_work_is_in_the_ledger() {
    let gate = CancelGate::new_clock_free();
    // The pair at rest and relaxed, pushed apart by 3 N held on its
    // connector; the connector's heat goes to one thermal body.
    let start = pair_state(0.0, 0.0);
    let component = pair().driven(0, DriveFrame::MaterialFirst).expect("actuator");
    let mut builder = SystemBuilder::new();
    let moving_node = builder.add("network", component, start.clone()).expect("network");
    let cell = builder
        .add("cell", HeatCapacity::new(5.0).expect("capacity"), 300.0_f64)
        .expect("cell");
    for index in 0..2 {
        let port = builder.port(moving_node, &format!("force{index}")).expect("port");
        builder
            .boundary(port, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .expect("free");
    }
    builder
        .boundary(
            builder.port(moving_node, "drive0").expect("port"),
            Box::new(|_, _, out: &mut [f64]| out.copy_from_slice(&[3.0, 0.0, 0.0])),
        )
        .expect("held force");
    builder
        .boundary(
            builder.port(moving_node, "couple0").expect("port"),
            Box::new(|_, _, out: &mut [f64]| out.fill(0.0)),
        )
        .expect("no couple");
    builder
        .heat_route(
            builder.port(moving_node, "heat0").expect("port"),
            &[(builder.port(cell, "heat").expect("port"), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .expect("route");
    let mut system = builder.seal().expect("sealed");
    let heat_port = system.port(moving_node, "heat0").expect("heat port");
    assert!(system.guess(heat_port, &[300.0]));

    // The owner marched alone over the same steps with the same effort.
    let idle = [Wrench::ZERO; 2];
    let held = [HeldInternalEffort::MaterialFirst([3.0, 0.0, 0.0, 0.0, 0.0, 0.0])];
    let mut twin = start.clone();
    let (dt_s, steps) = (0.0125, 20_u32);
    let (mut heat_j, mut work_j) = (0.0_f64, 0.0_f64);
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            system.step(dt_s, cx).expect("kernel step");
            let native = step_relaxation_network_with_loads(
                &twin.bodies,
                &pair_properties(),
                &pair_links(),
                &twin.memory,
                &held,
                HeldExternalLoads {
                    wrenches: &idle,
                    gravity: Gravity::ZERO,
                },
                dt_s,
                controls(),
                || false,
            )
            .expect("the owner takes the step");
            heat_j += native.connector_dissipation_j[0];
            work_j += native.connector_actuator_work_j[0];
            twin = RelaxationState {
                bodies: native.states_after().to_vec(),
                memory: native.memory_after().clone(),
            };
        }
    });
    let reached = system.state::<RelaxationState>(moving_node).expect("state");
    assert_eq!(reached, &twin);
    let [u, q, heat, work] = driven_pair_reference(3.0, 0.0, 0.0, dt_s * f64::from(steps));
    assert!((relative_speed(reached) - u).abs() <= 1.0e-9);
    assert!((reached.memory[0][0][0] - q).abs() <= 1.0e-9);
    assert!((heat_j - heat).abs() <= 1.0e-9 && (work_j - work).abs() <= 1.0e-9);

    // The ledger: the actuator's work came in through the drive boundary
    // and nowhere else, and it is stored in the pair or is heat in the cell.
    let ledger = system.ledger();
    let boundary = |port: &str| {
        ledger
            .boundaries
            .iter()
            .find(|account| account.port == port)
            .unwrap_or_else(|| panic!("boundary {port}"))
            .energy_in_j
            .value()
    };
    assert!(work_j > 0.05 && heat_j > 0.01, "work {work_j} J and heat {heat_j} J to measure");
    assert!((boundary("network.drive0") - work_j).abs() <= 1.0e-13 * work_j);
    for idle_port in ["network.couple0", "network.force0", "network.force1"] {
        assert_eq!(boundary(idle_port), 0.0, "{idle_port}");
    }
    assert!((ledger.boundary_energy_in_j() - work_j).abs() <= 1.0e-13 * work_j);
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
    let warmed = account("cell").stored_change_j.value();
    assert!((warmed - heat_j).abs() <= 1.0e-11 * heat_j, "{warmed} vs {heat_j}");
    let stored = account("network").stored_change_j.value();
    assert!((work_j - heat_j - stored).abs() <= 1.0e-10);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10);
}

#[test]
fn actuators_given_out_of_connector_order_keep_their_own_ports() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    // Connector 1 first (world axes), then connector 0 (material axes).
    let component = network(gravity)
        .driven(1, DriveFrame::World)
        .expect("actuator on connector 1")
        .driven(0, DriveFrame::MaterialFirst)
        .expect("actuator on connector 0");
    let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(
        names,
        ["force0", "force1", "force2", "heat0", "heat1", "drive1", "couple1", "drive0", "couple0"]
    );
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let efforts = drive_efforts();
    let dt_s = 2.0e-2;
    let stale = stale_after(&start, &efforts, gravity);
    // The inputs follow the ports: connector 1's effort, then connector 0's.
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0]);
    inputs.extend(&DRIVE[6..]);
    inputs.extend(&DRIVE[..6]);
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
    let Observed::Stepped { outputs, .. } = observed else {
        panic!("the network takes the step");
    };
    let native = owner_step_driven(&start, &held_wrenches(), &efforts, gravity, dt_s);
    for (slot, connector) in [(0, 1), (1, 0)] {
        let (port_j, terms_j) = drive_port_energy_j((3, 2, slot), dt_s, &inputs, &outputs);
        let native_j = native.connector_actuator_work_j[connector];
        assert!(
            (port_j - native_j).abs() <= 1.0e-13 * terms_j,
            "connector {connector}: ports {port_j:e} J, owner {native_j:e} J"
        );
    }
    // The two connectors do different work, so a swap would show.
    let [first, second] = [native.connector_actuator_work_j[0], native.connector_actuator_work_j[1]];
    assert!((first - second).abs() > 1.0e-4, "{first} J and {second} J");
    // The same physics as with the actuators given in connector order.
    let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
    // A refusal names the connector, not the position of its ports.
    for (index, connector) in [(11, 1), (16, 1), (17, 0), (22, 0)] {
        let mut wild = inputs.clone();
        wild[index] = f64::NAN;
        assert_eq!(
            propose(&component, &start, dt_s, &wild, &gate).unwrap_err(),
            RelaxationNetworkRefusal::Effort { connector }
        );
    }
}

#[test]
fn an_actuator_holding_nothing_changes_nothing_and_still_reports_its_rate() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let dt_s = 2.0e-2;
    let mut plain = HELD.to_vec();
    plain.extend([300.0, 320.0]);
    let mut inputs = plain.clone();
    inputs.extend([0.0; 12]);
    let (without, bare, _) =
        propose(&network(gravity), &start, dt_s, &plain, &gate).expect("step");
    let (with, receipt, outputs) =
        propose(&driven_network(gravity), &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(with, without);
    assert_eq!(receipt, bare);
    // Connector 0's rate, in the first body's material axes, is the mean of
    // the two half-step twists the owner reports.
    let native = owner_step(&start, &held_wrenches(), gravity, dt_s);
    let [initial, endpoint] = &native.half_sources;
    for c in 0..6 {
        let mean =
            0.5 * (initial.average_relative_twist[0][c] + endpoint.average_relative_twist[0][c]);
        assert_eq!(outputs[11 + c], mean, "component {c}");
    }
    for slot in 0..2 {
        // No work, and the relative rate of the connector all the same.
        assert_eq!(drive_port_energy_j((3, 2, slot), dt_s, &inputs, &outputs).0, 0.0);
        let rate = &outputs[11 + 6 * slot..17 + 6 * slot];
        assert!(rate.iter().any(|value| value.abs() > 0.05), "{rate:?}");
    }
}

/// The pair's connector with its bodies named the other way round: body 1
/// is the connector's first body.
fn reversed_pair() -> (Vec<RelaxationNetworkLink>, RelaxationState) {
    let mut links = pair_links();
    links[0].viscous.bodies = [1, 0];
    let p = pair_properties();
    let state = RelaxationState {
        bodies: vec![
            turn(
                body(p[0], Vec3::ZERO, Vec3::new(-0.6, 0.2, 0.1), Vec3::new(0.4, -0.6, 0.9)),
                UnitQuaternion::new(0.85, 0.3, -0.2, 0.4).expect("orientation"),
            ),
            turn(
                body(p[1], Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.5, -0.1, 0.3), Vec3::new(-0.2, 0.5, 0.3)),
                UnitQuaternion::new(0.6, -0.5, 0.4, 0.1).expect("orientation"),
            ),
        ],
        memory: vec![vec![[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]]],
    };
    (links, state)
}

#[test]
fn a_world_drive_turns_with_its_connectors_first_body_whichever_that_is() {
    let gate = CancelGate::new_clock_free();
    // Connector 0 with body 1 as its first body: the pose that matters for
    // a world-held effort is body 1's, and both bodies are turned (each
    // differently) and spinning.
    let (links, start) = reversed_pair();
    let held = [0.4, -0.3, 0.2, 0.05, 0.1, -0.08];
    let mut inputs = vec![0.0; 6];
    inputs.push(300.0);
    inputs.extend(held);
    let idle = [Wrench::ZERO; 2];
    for frame in [DriveFrame::World, DriveFrame::MaterialFirst] {
        let component =
            RelaxationNetwork::new(pair_properties(), links.clone(), Gravity::ZERO, controls())
                .expect("network")
                .driven(0, frame)
                .expect("actuator");
        let effort = match frame {
            DriveFrame::World => HeldInternalEffort::World(held),
            DriveFrame::MaterialFirst => HeldInternalEffort::MaterialFirst(held),
        };
        let owner = |dt_s: f64| {
            step_relaxation_network_with_loads(
                &start.bodies,
                &pair_properties(),
                &links,
                &start.memory,
                &[effort],
                HeldExternalLoads {
                    wrenches: &idle,
                    gravity: Gravity::ZERO,
                },
                dt_s,
                controls(),
                || false,
            )
            .expect("the owner takes the step")
        };
        let stale = {
            let step = owner(0.1);
            RelaxationState {
                bodies: step.states_after().to_vec(),
                memory: step.memory_after().clone(),
            }
        };
        let dt_s = 2.0e-2;
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        let observed =
            with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
        let Observed::Stepped { outputs, .. } = observed else {
            panic!("the network takes the step");
        };
        let native = owner(dt_s);
        let (port_j, terms_j) = drive_port_energy_j((2, 1, 0), dt_s, &inputs, &outputs);
        let native_j = native.connector_actuator_work_j[0];
        assert!(native_j.abs() > 1.0e-3, "{frame:?}: the actuator does work: {native_j} J");
        assert!(
            (port_j - native_j).abs() <= 1.0e-13 * terms_j,
            "{frame:?}: ports {port_j:e} J, owner {native_j:e} J"
        );
    }
}

#[test]
fn half_step_works_that_cancel_leave_an_allowance_for_their_rounding() {
    let gate = CancelGate::new_clock_free();
    // The pair with no Kelvin viscosity and a lossless branch is an
    // oscillator of angular frequency 1.3 rad/s. Pushed from rest by a held
    // force for one period, it is back at rest and relaxed: the actuator's
    // work in the first half step is taken back in the second, nothing is
    // stored and nothing is heat. Every total is then zero to rounding,
    // while each half step's work is not.
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = PAIR_ROOT_K;
    let lossless = SpatialRelaxationBranch::from_equilibrium_gram_rows(rows).expect("branch");
    let links = vec![RelaxationNetworkLink {
        viscous: ViscousNetworkLink {
            bodies: [0, 1],
            arms_body_m: [Vec3::ZERO; 2],
            damping: damping([0.0; 6]),
        },
        branches: vec![lossless],
    }];
    let component =
        RelaxationNetwork::new(pair_properties(), links.clone(), Gravity::ZERO, controls())
            .expect("network")
            .driven(0, DriveFrame::MaterialFirst)
            .expect("actuator");
    let force = 0.7;
    let mut inputs = vec![0.0; 6];
    inputs.push(300.0);
    inputs.extend([force, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let period_s = core::f64::consts::TAU / PAIR_ROOT_K;
    let idle = [Wrench::ZERO; 2];
    let stale = pair_state(0.3, -0.2);
    // The last case spins the bodies against each other about the
    // connector's axis at 2 rad/s (2 J of rotation). No couple is held, so
    // the spin does no work on the actuator.
    let p = pair_properties();
    let spinning = RelaxationState {
        bodies: vec![
            body(p[0], Vec3::ZERO, Vec3::ZERO, Vec3::new(-2.0, 0.0, 0.0)),
            body(p[1], Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0)),
        ],
        memory: vec![vec![[0.0; 6]]],
    };
    for (stretch, start) in [
        (1.0, pair_state(0.0, 0.0)),
        (0.999, pair_state(0.0, 0.0)),
        (1.002, pair_state(0.0, 0.0)),
        (1.0, spinning),
    ] {
        let dt_s = stretch * period_s;
        let native = step_relaxation_network_with_loads(
            &start.bodies,
            &pair_properties(),
            &links,
            &start.memory,
            &[HeldInternalEffort::MaterialFirst([force, 0.0, 0.0, 0.0, 0.0, 0.0])],
            HeldExternalLoads {
                wrenches: &idle,
                gravity: Gravity::ZERO,
            },
            dt_s,
            controls(),
            || false,
        )
        .expect("the owner takes the step");
        let [first, second] = [
            native.half_sources[0].connector_actuator_work_j[0],
            native.half_sources[1].connector_actuator_work_j[0],
        ];
        // 2 F^2 / k^2 = 0.58 J in, and the same out.
        assert!(first > 0.5 && second < -0.5, "{first} J then {second} J");
        assert!(native.actuator_work_j.abs() <= 1.0e-4 * first, "{} J", native.actuator_work_j);
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        let observed =
            with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
        let Observed::Stepped { receipt, .. } = observed else {
            panic!("the network takes the step");
        };
        // The allowance is the rounding of what is there: 64 ulp of each
        // half's work, of each half's port product (the same size) and of
        // the energy stored at the two ends, with a quarter to spare. A
        // force is paired with the linear rate only: the spin of the last
        // case adds its stored energy and nothing else. Measured: 3.3e-14 J
        // without the spin and 9.0e-14 J with it.
        let before = component.stored_energy_j(&start).expect("energy");
        let (after, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
        let after = component.stored_energy_j(&after).expect("energy");
        assert!(receipt.energy_tolerance_j >= 64.0 * f64::EPSILON * (first - second));
        assert!(
            receipt.energy_tolerance_j
                <= 80.0 * f64::EPSILON * (2.0 * (first - second) + before + after),
            "{:e} J with {before} J stored",
            receipt.energy_tolerance_j
        );
    }
}

#[test]
fn drives_in_world_axes_with_a_couple_reach_the_ledger_through_the_kernel() {
    let gate = CancelGate::new_clock_free();
    // The three-body network, its first two bodies turned and spinning,
    // with both actuators holding a force and a couple (connector 0 in
    // material axes, connector 1 in world axes), held forces on the bodies,
    // and each connector's heat to its own sink.
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let mut builder = SystemBuilder::new();
    let moving_node = builder
        .add("network", driven_network(Gravity::ZERO), start.clone())
        .expect("network");
    for body in 0..3 {
        let port = builder.port(moving_node, &format!("force{body}")).expect("port");
        let held: [f64; 3] = core::array::from_fn(|axis| HELD[3 * body + axis]);
        builder
            .boundary(port, Box::new(move |_, _, out: &mut [f64]| out.copy_from_slice(&held)))
            .expect("held force");
    }
    for (slot, name) in ["drive0", "couple0", "drive1", "couple1"].into_iter().enumerate() {
        let port = builder.port(moving_node, name).expect("port");
        let held: [f64; 3] = core::array::from_fn(|axis| DRIVE[3 * slot + axis]);
        builder
            .boundary(port, Box::new(move |_, _, out: &mut [f64]| out.copy_from_slice(&held)))
            .expect("held effort");
    }
    for connector in 0..2 {
        builder
            .heat_route(
                builder.port(moving_node, &format!("heat{connector}")).expect("port"),
                &[],
                Some(heat_sink(&format!("sink{connector}"), 300.0).expect("sink")),
                RouteTemperature::Fixed(300.0),
            )
            .expect("route");
    }
    let mut system = builder.seal().expect("sealed");

    // The owner marched alone over the same steps.
    let efforts = drive_efforts();
    let mut twin = start.clone();
    let (dt_s, steps) = (4.0e-3, 5_u32);
    let (mut works_j, mut heats_j, mut applied_j) = ([0.0_f64; 2], [0.0_f64; 2], 0.0_f64);
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            system.step(dt_s, cx).expect("kernel step");
            let native = owner_step_driven(&twin, &held_wrenches(), &efforts, Gravity::ZERO, dt_s);
            for connector in 0..2 {
                works_j[connector] += native.connector_actuator_work_j[connector];
                heats_j[connector] += native.connector_dissipation_j[connector];
            }
            applied_j += native.external_loads.applied_work_j;
            twin = RelaxationState {
                bodies: native.states_after().to_vec(),
                memory: native.memory_after().clone(),
            };
        }
    });
    assert_eq!(system.state::<RelaxationState>(moving_node).expect("state"), &twin);

    let ledger = system.ledger();
    let boundary = |port: &str| {
        ledger
            .boundaries
            .iter()
            .find(|account| account.port == port)
            .unwrap_or_else(|| panic!("boundary {port}"))
            .energy_in_j
            .value()
    };
    let size = works_j[0].abs() + works_j[1].abs() + heats_j[0] + heats_j[1] + applied_j.abs();
    for connector in 0..2 {
        let force_j = boundary(&format!("network.drive{connector}"));
        let couple_j = boundary(&format!("network.couple{connector}"));
        // Both ports of an actuator carry work (-2.1e-2 J and -2.2e-4 J on
        // connector 0, -8.1e-3 J and -1.7e-4 J on connector 1), and together
        // the owner's. Measured: equal to the bit for the actuators and the
        // sinks, 9e-19 J apart for the body forces, a ledger residual of
        // 2e-19 J; `size` is 0.10 J.
        assert!(force_j.abs() > 1.0e-5 && couple_j.abs() > 1.0e-7, "{force_j} J, {couple_j} J");
        assert!(
            (force_j + couple_j - works_j[connector]).abs() <= 1.0e-14 * size,
            "connector {connector}: {} vs {}",
            force_j + couple_j,
            works_j[connector]
        );
        let sink_j = ledger
            .sinks
            .iter()
            .find(|sink| sink.name == format!("sink{connector}"))
            .expect("sink")
            .heat_j
            .value();
        assert!((sink_j - heats_j[connector]).abs() <= 1.0e-14 * size);
    }
    let pushed: f64 = (0..3).map(|body| boundary(&format!("network.force{body}"))).sum();
    assert!((pushed - applied_j).abs() <= 1.0e-14 * size, "{pushed} vs {applied_j}");
    assert!(ledger.energy_residual_j().abs() <= 1.0e-14 * size);
}

/// A synthetic damping law: the built Kelvin form times
/// `1 + beta_e (T - 300 K)`, with its own `beta_e` for each connector,
/// defined from 280 K to 400 K. Not a material.
#[derive(Debug)]
struct Stiffening {
    betas_per_k: Vec<f64>,
}

/// The slopes of the three-body network's two connectors.
const BETAS: [f64; 2] = [0.02, 0.03];

impl Stiffening {
    fn scale(&self, connector: usize, temperature_k: f64) -> f64 {
        1.0 + self.betas_per_k[connector] * (temperature_k - 300.0)
    }
}

impl DampingLaw for Stiffening {
    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, String> {
        if !(280.0..=400.0).contains(&temperature_k) {
            return Err(format!("no form at {temperature_k} K"));
        }
        let scale = self.scale(connector, temperature_k);
        let shifted = reference.native_coefficients().map(|row| row.map(|c| scale * c));
        PassiveBushingMatrix::try_new(shifted, scales()).map_err(|error| format!("{error:?}"))
    }
}

/// A synthetic shift law: `exp(A_e (1 / T - 1 / 300 K))`, with its own
/// `A_e` for each connector, defined from 280 K to 400 K. Not a material.
#[derive(Debug)]
struct Arrhenius {
    activations_k: Vec<f64>,
}

/// The activation temperatures of the three-body network's two connectors.
const ACTIVATIONS: [f64; 2] = [2000.0, 1200.0];

impl ShiftLaw for Arrhenius {
    fn factor(&self, connector: usize, temperature_k: f64) -> Result<f64, String> {
        if !(280.0..=400.0).contains(&temperature_k) {
            return Err(format!("no factor at {temperature_k} K"));
        }
        Ok((self.activations_k[connector] * (1.0 / temperature_k - 1.0 / 300.0)).exp())
    }
}

/// A shift law that gives the same factor whatever it is asked.
#[derive(Debug)]
struct Constant(f64);

impl ShiftLaw for Constant {
    fn factor(&self, _connector: usize, _temperature_k: f64) -> Result<f64, String> {
        Ok(self.0)
    }
}

/// The three-body network's connectors with their Kelvin forms multiplied
/// by `scales`, built by hand.
fn stiffened_links(scales: [f64; 2]) -> Vec<RelaxationNetworkLink> {
    let mut links = links();
    links[0].viscous.damping = damping([0.6, 0.4, 0.3, 0.0, 0.05, 0.0].map(|c| scales[0] * c));
    links[1].viscous.damping = damping([0.9, 0.0, 0.0, 0.0, 0.0, 0.0].map(|c| scales[1] * c));
    links
}

#[test]
fn a_kelvin_form_that_follows_temperature_is_the_owner_on_the_shifted_forms() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let law = || Rc::new(Stiffening { betas_per_k: BETAS.to_vec() }) as Rc<dyn DampingLaw>;
    let component = network(gravity).temperature_dependent(law()).expect("law");
    let names: Vec<&str> = component.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(
        names,
        ["force0", "force1", "force2", "heat0", "heat1", "temperature0", "temperature1"]
    );
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let stale = {
        let step = owner_step(&start, &held_wrenches(), gravity, 0.1);
        RelaxationState {
            bodies: step.states_after().to_vec(),
            memory: step.memory_after().clone(),
        }
    };
    // Heat released at 300 K and 320 K; the damping read at 310 K and 330 K
    // (the two are authored independently), each connector by its own
    // slope: forms 1.2 and 1.9 times the built ones.
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0, 310.0, 330.0]);
    let dt_s = 2.0e-2;
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
    let Observed::Stepped {
        receipt, outputs, ..
    } = observed
    else {
        panic!("the network takes the step");
    };
    let idle = [HeldInternalEffort::MaterialFirst([0.0; 6]); 2];
    let owner = |links: &[RelaxationNetworkLink]| {
        step_relaxation_network_with_loads(
            &start.bodies,
            &properties(),
            links,
            &start.memory,
            &idle,
            HeldExternalLoads {
                wrenches: &held_wrenches(),
                gravity,
            },
            dt_s,
            controls(),
            || false,
        )
        .expect("the owner takes the step")
    };
    let native = owner(&stiffened_links([1.0 + BETAS[0] * 10.0, 1.0 + BETAS[1] * 30.0]));
    let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
    assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
    assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
    // Entropy where the heat is released (300 K and 320 K), not where the
    // damping is read.
    let entropy =
        native.connector_dissipation_j[0] / 300.0 + native.connector_dissipation_j[1] / 320.0;
    assert!((receipt.entropy_generated_j_per_k - entropy).abs() <= 1.0e-14 * entropy);
    let at_signals =
        native.connector_dissipation_j[0] / 310.0 + native.connector_dissipation_j[1] / 330.0;
    assert!((entropy - at_signals).abs() > 1.0e-3 * entropy);
    for connector in 0..2 {
        assert_eq!(outputs[9 + connector], -native.connector_dissipation_j[connector] / dt_s);
    }
    // The law matters: the built forms heat less.
    let built = owner(&links());
    for connector in 0..2 {
        let (warm, cool) = (
            native.connector_dissipation_j[connector],
            built.connector_dissipation_j[connector],
        );
        assert!(warm > 1.02 * cool, "connector {connector}: {warm} J against {cool} J");
    }
    // At the reference temperature the law returns the built form, and the
    // step is the one of the network that follows nothing.
    let mut reference = HELD.to_vec();
    reference.extend([300.0, 320.0, 300.0, 300.0]);
    let (at_reference, _, _) = propose(&component, &start, dt_s, &reference, &gate).expect("step");
    let (plain, _, _) =
        propose(&network(gravity), &start, dt_s, &reference[..11], &gate).expect("step");
    assert_eq!(at_reference, plain);
}

/// One order in which the three-body network was given its actuators (on
/// connector 0 in material axes, on connector 1 in world axes) and its
/// temperature signals: the network, its input vector, where each
/// connector's efforts start in the inputs and its rates in the outputs,
/// and where the signals start in the inputs.
struct Layout {
    component: RelaxationNetwork,
    names: [&'static str; 6],
    inputs: Vec<f64>,
    drives: [(usize, usize); 2],
    signals: usize,
}

/// Three orders: actuators then signals; signals then actuators; and an
/// actuator on each side of the signals, connector 1 first. A signal input
/// has no output, so the rates always start after the two heat flows.
fn layouts(gravity: Gravity) -> [Layout; 3] {
    let law = || Rc::new(Stiffening { betas_per_k: BETAS.to_vec() }) as Rc<dyn DampingLaw>;
    let mut head = HELD.to_vec();
    head.extend([300.0, 320.0]);
    let signals = [310.0, 330.0];
    let join = |parts: [&[f64]; 3]| {
        let mut inputs = head.clone();
        for part in parts {
            inputs.extend(part);
        }
        inputs
    };
    [
        Layout {
            component: driven_network(gravity).temperature_dependent(law()).expect("law"),
            names: ["drive0", "couple0", "drive1", "couple1", "temperature0", "temperature1"],
            inputs: join([&DRIVE[..6], &DRIVE[6..], &signals]),
            drives: [(11, 11), (17, 17)],
            signals: 23,
        },
        Layout {
            component: network(gravity)
                .temperature_dependent(law())
                .expect("law")
                .driven(0, DriveFrame::MaterialFirst)
                .expect("actuator on connector 0")
                .driven(1, DriveFrame::World)
                .expect("actuator on connector 1"),
            names: ["temperature0", "temperature1", "drive0", "couple0", "drive1", "couple1"],
            inputs: join([&signals, &DRIVE[..6], &DRIVE[6..]]),
            drives: [(13, 11), (19, 17)],
            signals: 11,
        },
        Layout {
            component: network(gravity)
                .driven(1, DriveFrame::World)
                .expect("actuator on connector 1")
                .temperature_dependent(law())
                .expect("law")
                .driven(0, DriveFrame::MaterialFirst)
                .expect("actuator on connector 0"),
            names: ["drive1", "couple1", "temperature0", "temperature1", "drive0", "couple0"],
            inputs: join([&DRIVE[6..], &signals, &DRIVE[..6]]),
            drives: [(19, 17), (11, 11)],
            signals: 17,
        },
    ]
}

#[test]
fn temperature_signals_and_actuators_keep_their_places_in_any_order() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let efforts = drive_efforts();
    let dt_s = 2.0e-2;
    let stale = stale_after(&start, &efforts, gravity);
    // The owner on the hand-shifted forms with the same efforts.
    let native = step_relaxation_network_with_loads(
        &start.bodies,
        &properties(),
        &stiffened_links([1.0 + BETAS[0] * 10.0, 1.0 + BETAS[1] * 30.0]),
        &start.memory,
        &efforts,
        HeldExternalLoads {
            wrenches: &held_wrenches(),
            gravity,
        },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step");
    for layout in layouts(gravity) {
        let Layout {
            component,
            names,
            inputs,
            drives,
            signals,
        } = layout;
        let found: Vec<&str> =
            component.ports()[5..].iter().map(|port| port.id().as_str()).collect();
        assert_eq!(found, names);
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        let observed =
            with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
        let Observed::Stepped { outputs, .. } = observed else {
            panic!("the network takes the step");
        };
        assert_eq!(outputs.len(), 23, "{names:?}");
        let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
        assert_eq!(next.bodies.as_slice(), native.states_after(), "{names:?}");
        assert_eq!(&next.memory, native.memory_after(), "{names:?}");
        // Each actuator's ports carry the owner's work on its connector.
        for (connector, offsets) in drives.into_iter().enumerate() {
            let (port_j, terms_j) = drive_energy_at(offsets, dt_s, &inputs, &outputs);
            let native_j = native.connector_actuator_work_j[connector];
            assert!(native_j.abs() > 1.0e-5);
            assert!(
                (port_j - native_j).abs() <= 1.0e-13 * terms_j,
                "{names:?}, connector {connector}: ports {port_j:e} J, owner {native_j:e} J"
            );
        }
        // A wrong signal names its connector.
        let mut wild = inputs.clone();
        wild[signals + 1] = 0.0;
        assert_eq!(
            propose(&component, &start, dt_s, &wild, &gate).unwrap_err(),
            RelaxationNetworkRefusal::DampingTemperature {
                connector: 1,
                kelvin: 0.0
            }
        );
    }
}

/// The pair's equations with Kelvin viscosity `c` and memory relaxing at
/// `rate_per_s` (RK4, 16384 steps): `u' = -c u - k q`, `q' = k u - rate q`,
/// heat `c u^2 + rate q^2`. Returns `[u, q, heat]` at `t_s`.
fn pair_equations(c: f64, rate_per_s: f64, u: f64, q: f64, t_s: f64) -> [f64; 3] {
    let f = |s: [f64; 3]| {
        [
            -c * s[0] - PAIR_ROOT_K * s[1],
            PAIR_ROOT_K * s[0] - rate_per_s * s[1],
            c * s[0] * s[0] + rate_per_s * s[1] * s[1],
        ]
    };
    let add = |a: [f64; 3], b: [f64; 3], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [u, q, 0.0];
    let h = t_s / 16384.0;
    for _ in 0..16384 {
        let a = f(state);
        let b = f(add(state, a, 0.5 * h));
        let c = f(add(state, b, 0.5 * h));
        let d = f(add(state, c, h));
        for i in 0..3 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

/// `steps` steps of `component` (a pair that follows temperature, no
/// actuator) with its signal at `kelvin`. Returns the reached state and the
/// connector's heat.
fn run_pair_at(
    component: &RelaxationNetwork,
    start: &RelaxationState,
    kelvin: f64,
    (dt_s, steps): (f64, u32),
    gate: &CancelGate,
) -> (RelaxationState, f64) {
    let inputs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0, kelvin];
    let mut state = start.clone();
    let mut heat_j = 0.0_f64;
    for _ in 0..steps {
        let (next, _, outputs) = propose(component, &state, dt_s, &inputs, gate).expect("step");
        heat_j -= dt_s * outputs[6];
        state = next;
    }
    (state, heat_j)
}

#[test]
fn kinetics_that_follow_temperature_scale_viscosity_and_relaxation_time_together() {
    let gate = CancelGate::new_clock_free();
    let law =
        || Rc::new(Arrhenius { activations_k: ACTIVATIONS.to_vec() }) as Rc<dyn ShiftLaw>;

    // The three-body network against the owner's held-kinetics entry, with
    // the factors computed here: one per connector, for its Kelvin form and
    // for each of its branches.
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    let component = network(gravity).shifted(law()).expect("law");
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let stale = {
        let step = owner_step(&start, &held_wrenches(), gravity, 0.1);
        RelaxationState {
            bodies: step.states_after().to_vec(),
            memory: step.memory_after().clone(),
        }
    };
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0, 310.0, 330.0]);
    let dt_s = 2.0e-2;
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
    let Observed::Stepped { receipt, .. } = observed else {
        panic!("the network takes the step");
    };
    // Each connector by its own activation temperature, at its own signal.
    let factors = [
        (ACTIVATIONS[0] * (1.0 / 310.0_f64 - 1.0 / 300.0)).exp(),
        (ACTIVATIONS[1] * (1.0 / 330.0_f64 - 1.0 / 300.0)).exp(),
    ];
    assert!(factors[0] < 0.85 && factors[1] < 0.75, "{factors:?}");
    assert!((factors[0] - factors[1]).abs() > 0.05, "{factors:?}");
    let idle = [HeldInternalEffort::MaterialFirst([0.0; 6]); 2];
    let native = step_relaxation_network_with_loads_and_kinetics(
        &start.bodies,
        &properties(),
        &links(),
        &start.memory,
        &idle,
        HeldRelaxationKinetics {
            viscosity_scales: &factors,
            relaxation_time_scales: &[vec![factors[0]; 2], vec![factors[1]; 1]],
        },
        HeldExternalLoads {
            wrenches: &held_wrenches(),
            gravity,
        },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step");
    let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
    assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
    assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
}

#[test]
fn a_shifted_pair_obeys_its_own_equations_with_viscosity_and_time_scaled() {
    let gate = CancelGate::new_clock_free();
    let law = || Rc::new(Arrhenius { activations_k: vec![2000.0] }) as Rc<dyn ShiftLaw>;
    // The pair closing at 2 m/s with a preloaded memory, at 330 K: its
    // viscosity and its relaxation time are both 0.545 of the built ones,
    // against the pair's own equations with those two numbers scaled.
    let factor = (2000.0_f64 * (1.0 / 330.0 - 1.0 / 300.0)).exp();
    let run = (0.0125, 20_u32);
    let t_s = run.0 * f64::from(run.1);
    let start = pair_state(-2.0, 0.6);
    let pair = pair().shifted(law()).expect("law");
    let (state, heat_j) = run_pair_at(&pair, &start, 330.0, run, &gate);
    let [u, q, heat] = pair_equations(factor * PAIR_C, 1.0 / (factor * PAIR_TAU_S), -2.0, 0.6, t_s);
    assert!((relative_speed(&state) - u).abs() <= 1.0e-9, "{} vs {u}", relative_speed(&state));
    assert!((state.memory[0][0][0] - q).abs() <= 1.0e-9, "{} vs {q}", state.memory[0][0][0]);
    assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");
    // The shift is not small: the built pair ends at -1.846 m/s with
    // 0.459 J of heat, the shifted one at -1.913 m/s with 0.332 J.
    let [built_u, built_q, built_heat] = pair_reference(-2.0, 0.6, t_s);
    assert!((u - built_u).abs() > 0.05, "{u} m/s against {built_u} m/s");
    assert!((heat - built_heat).abs() > 0.1, "{heat} J against {built_heat} J");
    // Scaling the viscosity alone (-1.931 m/s), or the relaxation time
    // alone (-1.829 m/s), is not it.
    let [only_c, _, _] = pair_equations(factor * PAIR_C, 1.0 / PAIR_TAU_S, -2.0, 0.6, t_s);
    let [only_tau, _, _] = pair_equations(PAIR_C, 1.0 / (factor * PAIR_TAU_S), -2.0, 0.6, t_s);
    assert!((u - only_c).abs() > 0.01 && (u - only_tau).abs() > 0.05);

    // At the reference temperature the factor is one and the pair is the
    // built pair.
    let (at_reference, reference_heat_j) = run_pair_at(&pair, &start, 300.0, run, &gate);
    assert!((relative_speed(&at_reference) - built_u).abs() <= 1.0e-9);
    assert!((at_reference.memory[0][0][0] - built_q).abs() <= 1.0e-9);
    assert!((reference_heat_j - built_heat).abs() <= 1.0e-9);

    // A lossless branch stays lossless under a shift: only the scaled
    // Kelvin viscosity heats.
    let elastic = elastic_pair().shifted(law()).expect("law");
    let (state, heat_j) = run_pair_at(&elastic, &start, 330.0, run, &gate);
    let [u, q, heat] = pair_equations(factor * PAIR_C, 0.0, -2.0, 0.6, t_s);
    assert!((relative_speed(&state) - u).abs() <= 1.0e-9);
    assert!((state.memory[0][0][0] - q).abs() <= 1.0e-9);
    assert!((heat_j - heat).abs() <= 1.0e-9, "{heat_j} vs {heat}");
}

#[test]
fn without_branches_a_shift_factor_and_a_scaled_kelvin_form_are_the_same_connector() {
    let gate = CancelGate::new_clock_free();
    // The pair with no branch, by the owner's kinetics entry (a factor of
    // 1.3) and by its plain entry on a Kelvin form multiplied by 1.3 (the
    // damping law at 315 K): two ways to the same connector.
    let branchless = || {
        let mut links = pair_links();
        links[0].branches.clear();
        RelaxationNetwork::new(pair_properties(), links, Gravity::ZERO, controls()).expect("network")
    };
    let by_factor = branchless().shifted(Rc::new(Constant(1.3))).expect("law");
    let by_form = branchless()
        .temperature_dependent(Rc::new(Stiffening { betas_per_k: vec![0.02] }))
        .expect("law");
    let start = RelaxationState {
        bodies: pair_state(-2.0, 0.0).bodies,
        memory: vec![Vec::new()],
    };
    let run = (0.0125, 20_u32);
    let (factor_state, factor_heat_j) = run_pair_at(&by_factor, &start, 315.0, run, &gate);
    let (form_state, form_heat_j) = run_pair_at(&by_form, &start, 315.0, run, &gate);
    let (a, b) = (relative_speed(&factor_state), relative_speed(&form_state));
    // Measured: equal to the bit, speed and heat.
    assert!((a - b).abs() <= 1.0e-12, "{a} vs {b}");
    assert!((factor_heat_j - form_heat_j).abs() <= 1.0e-12);
    // And both are the viscous decay `u' = -1.3 c u`.
    let t_s = run.0 * f64::from(run.1);
    let u = -2.0 * (-1.3 * PAIR_C * t_s).exp();
    assert!((a - u).abs() <= 1.0e-9, "{a} vs {u}");
    assert!((factor_heat_j - 0.5 * (4.0 - u * u)).abs() <= 1.0e-9);
}

#[test]
fn temperature_laws_are_refused_where_they_have_nothing_to_give() {
    let gate = CancelGate::new_clock_free();
    let stiffening = || Rc::new(Stiffening { betas_per_k: BETAS.to_vec() }) as Rc<dyn DampingLaw>;
    let arrhenius =
        || Rc::new(Arrhenius { activations_k: ACTIVATIONS.to_vec() }) as Rc<dyn ShiftLaw>;
    // One way of following temperature, not two.
    assert_eq!(
        network(Gravity::ZERO)
            .temperature_dependent(stiffening())
            .expect("law")
            .shifted(arrhenius())
            .unwrap_err(),
        RelaxationNetworkBuildError::Thermal
    );
    assert_eq!(
        network(Gravity::ZERO)
            .shifted(arrhenius())
            .expect("law")
            .shifted(arrhenius())
            .unwrap_err(),
        RelaxationNetworkBuildError::Thermal
    );
    let start = RelaxationState {
        bodies: moving(),
        memory: preloaded(),
    };
    let mut inputs = vec![0.0; 9];
    inputs.extend([300.0, 300.0, 310.0, 330.0]);
    let by_form = network(Gravity::ZERO).temperature_dependent(stiffening()).expect("law");
    let by_factor = network(Gravity::ZERO).shifted(arrhenius()).expect("law");
    for component in [&by_form, &by_factor] {
        assert!(propose(component, &start, 4.0e-3, &inputs, &gate).is_ok());
        // A signal that is no temperature.
        for kelvin in [0.0, -5.0, f64::INFINITY] {
            let mut wild = inputs.clone();
            wild[11] = kelvin;
            assert_eq!(
                propose(component, &start, 4.0e-3, &wild, &gate).unwrap_err(),
                RelaxationNetworkRefusal::DampingTemperature {
                    connector: 0,
                    kelvin
                }
            );
        }
        let mut unknown = inputs.clone();
        unknown[12] = f64::NAN;
        assert!(matches!(
            propose(component, &start, 4.0e-3, &unknown, &gate),
            Err(RelaxationNetworkRefusal::DampingTemperature { connector: 1, kelvin })
                if kelvin.is_nan()
        ));
        // A temperature outside the law's range, on the second connector.
        let mut hot = inputs.clone();
        hot[12] = 450.0;
        assert!(matches!(
            propose(component, &start, 4.0e-3, &hot, &gate),
            Err(RelaxationNetworkRefusal::Law { connector: 1, reason }) if reason.contains("450")
        ));
    }
    // A factor that is not a positive ratio of times is refused as the
    // law's, with or without branches (the owner alone would take a zero
    // factor on a connector that has no branch).
    let branchless = {
        let mut links = pair_links();
        links[0].branches.clear();
        links
    };
    let pair_start = RelaxationState {
        bodies: pair_state(-2.0, 0.0).bodies,
        memory: vec![Vec::new()],
    };
    let pair_inputs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0, 310.0];
    for factor in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let component = network(Gravity::ZERO).shifted(Rc::new(Constant(factor))).expect("law");
        assert!(
            matches!(
                propose(&component, &start, 4.0e-3, &inputs, &gate),
                Err(RelaxationNetworkRefusal::Law { connector: 0, reason })
                    if reason.contains("shift factor")
            ),
            "factor {factor}"
        );
        let component =
            RelaxationNetwork::new(pair_properties(), branchless.clone(), Gravity::ZERO, controls())
                .expect("network")
                .shifted(Rc::new(Constant(factor)))
                .expect("law");
        assert!(
            matches!(
                propose(&component, &pair_start, 4.0e-3, &pair_inputs, &gate),
                Err(RelaxationNetworkRefusal::Law { connector: 0, .. })
            ),
            "factor {factor}, no branch"
        );
    }
}

#[test]
fn a_temperature_held_on_the_signal_reaches_the_kinetics_through_the_kernel() {
    let gate = CancelGate::new_clock_free();
    // The pair closing at 2 m/s with a preloaded memory, its kinetics read
    // at 330 K held on the signal, its connector's heat to one thermal
    // body.
    let start = pair_state(-2.0, 0.6);
    let component = pair()
        .shifted(Rc::new(Arrhenius { activations_k: vec![2000.0] }))
        .expect("law");
    let mut builder = SystemBuilder::new();
    let moving_node = builder.add("network", component, start.clone()).expect("network");
    let cell = builder
        .add("cell", HeatCapacity::new(5.0).expect("capacity"), 300.0_f64)
        .expect("cell");
    for index in 0..2 {
        let port = builder.port(moving_node, &format!("force{index}")).expect("port");
        builder
            .boundary(port, Box::new(|_, _, out: &mut [f64]| out.fill(0.0)))
            .expect("free");
    }
    builder
        .boundary(
            builder.port(moving_node, "temperature0").expect("port"),
            Box::new(|_, _, out: &mut [f64]| out.fill(330.0)),
        )
        .expect("held temperature");
    builder
        .heat_route(
            builder.port(moving_node, "heat0").expect("port"),
            &[(builder.port(cell, "heat").expect("port"), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .expect("route");
    let mut system = builder.seal().expect("sealed");
    let heat_port = system.port(moving_node, "heat0").expect("heat port");
    assert!(system.guess(heat_port, &[300.0]));
    let (dt_s, steps) = (0.0125, 20_u32);
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            system.step(dt_s, cx).expect("kernel step");
        }
    });
    let factor = (2000.0_f64 * (1.0 / 330.0 - 1.0 / 300.0)).exp();
    let [u, q, heat] = pair_equations(
        factor * PAIR_C,
        1.0 / (factor * PAIR_TAU_S),
        -2.0,
        0.6,
        dt_s * f64::from(steps),
    );
    let reached = system.state::<RelaxationState>(moving_node).expect("state");
    assert!((relative_speed(reached) - u).abs() <= 1.0e-9);
    assert!((reached.memory[0][0][0] - q).abs() <= 1.0e-9);
    let ledger = system.ledger();
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).expect("node");
    let warmed = account("cell").stored_change_j.value();
    assert!((warmed - heat).abs() <= 1.0e-9, "{warmed} J against {heat} J");
    // A signal carries no energy: nothing entered through any boundary.
    assert_eq!(ledger.boundary_energy_in_j(), 0.0);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10);
}

#[test]
fn body_torques_sit_with_an_actuator_and_signals_on_the_relaxation_network() {
    let gate = CancelGate::new_clock_free();
    let gravity = Gravity::new(Vec3::new(0.0, -1.0, -9.81)).expect("gravity");
    // A torque port on body 1, an actuator on connector 0, the temperature
    // signals, then a torque port on body 2: the two torque ports sit apart,
    // with ports of every other kind between them.
    let component = network(gravity)
        .torqued(1)
        .expect("torque port on body 1")
        .driven(0, DriveFrame::MaterialFirst)
        .expect("actuator on connector 0")
        .temperature_dependent(Rc::new(Stiffening { betas_per_k: BETAS.to_vec() }))
        .expect("law")
        .torqued(2)
        .expect("torque port on body 2");
    let names: Vec<&str> =
        component.ports()[5..].iter().map(|port| port.id().as_str()).collect();
    assert_eq!(
        names,
        ["torque1", "drive0", "couple0", "temperature0", "temperature1", "torque2"]
    );
    let torques = [[0.03, -0.05, 0.02], [-0.04, 0.01, 0.06]];
    let mut inputs = HELD.to_vec();
    inputs.extend([300.0, 320.0]);
    inputs.extend(torques[0]);
    inputs.extend(&DRIVE[..6]);
    inputs.extend([310.0, 330.0]);
    inputs.extend(torques[1]);
    let start = RelaxationState {
        bodies: turned(),
        memory: preloaded(),
    };
    let efforts = [
        HeldInternalEffort::MaterialFirst(core::array::from_fn(|c| DRIVE[c])),
        HeldInternalEffort::MaterialFirst([0.0; 6]),
    ];
    let mut wrenches = held_wrenches();
    wrenches[1].torque_body = Vec3::new(torques[0][0], torques[0][1], torques[0][2]);
    wrenches[2].torque_body = Vec3::new(torques[1][0], torques[1][1], torques[1][2]);
    let dt_s = 2.0e-2;
    let stale = stale_after(&start, &efforts, gravity);
    let case = Case {
        state: &start,
        stale: &stale,
        dt_s,
        inputs: &inputs,
    };
    let observed = with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
    let Observed::Stepped {
        receipt, outputs, ..
    } = observed
    else {
        panic!("the network takes the step");
    };
    // A signal input has no output: forces 9, heat 2, torque 3, drive 6,
    // torque 3.
    assert_eq!(outputs.len(), 23);
    let native = step_relaxation_network_with_loads(
        &start.bodies,
        &properties(),
        &stiffened_links([1.0 + BETAS[0] * 10.0, 1.0 + BETAS[1] * 30.0]),
        &start.memory,
        &efforts,
        HeldExternalLoads {
            wrenches: &wrenches,
            gravity,
        },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step");
    let (next, _, _) = propose(&component, &start, dt_s, &inputs, &gate).expect("step");
    assert_eq!(next.bodies.as_slice(), native.states_after());
    assert_eq!(&next.memory, native.memory_after());
    assert_eq!(receipt.stored_change_j, native.mechanical_energy_change_j);
    assert_eq!(receipt.numerical_loss_j, -native.work_heat_residual_j);
    // The two torque ports (efforts at 11 and 22, rates at 11 and 20) carry
    // the owner's couple work; the actuator (efforts and rates at 14) its
    // actuator work.
    let energy = |input: usize, output: usize| {
        let (mut energy, mut size) = (0.0_f64, 0.0_f64);
        for axis in 0..3 {
            energy += inputs[input + axis] * outputs[output + axis];
            size += (inputs[input + axis] * outputs[output + axis]).abs();
        }
        (dt_s * energy, dt_s * size)
    };
    let (first_j, first_terms_j) = energy(11, 11);
    let (second_j, second_terms_j) = energy(22, 20);
    let couple_j = native.external_loads.applied_torque_work_j;
    assert!(couple_j.abs() > 1.0e-6, "the torques do work: {couple_j} J");
    assert!(
        (first_j + second_j - couple_j).abs() <= 1.0e-13 * (first_terms_j + second_terms_j),
        "ports {:e} J, owner {couple_j:e} J",
        first_j + second_j
    );
    let (drive_j, drive_terms_j) = drive_energy_at((14, 14), dt_s, &inputs, &outputs);
    let actuator_j = native.connector_actuator_work_j[0];
    assert!(actuator_j.abs() > 1.0e-5);
    assert!((drive_j - actuator_j).abs() <= 1.0e-13 * drive_terms_j);

    assert_eq!(
        network(gravity).torqued(3).unwrap_err(),
        RelaxationNetworkBuildError::Torque { body: 3 }
    );
    assert_eq!(
        network(gravity).torqued(1).expect("torque port").torqued(1).unwrap_err(),
        RelaxationNetworkBuildError::Torque { body: 1 }
    );
}
