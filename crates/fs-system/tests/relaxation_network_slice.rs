//! Owner slice: the fs-mbd rigid-body network with Maxwell memory on the
//! composition kernel.
//!
//! `fs_mbd::component::relaxation::RelaxationNetwork` puts the owner's
//! `step_relaxation_network_with_loads` behind the Component contract: the
//! bodies and the complete elastic memory are one state, and the memory is
//! stored energy. Here it passes the conformance check, equals the owner
//! driven directly, releases its preloaded memory as motion and heat, and
//! heats thermal bodies through routes.
//! Run with `cargo test -p fs-system --features slices`.

// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mbd::component::network::ViscousNetwork;
use fs_mbd::component::relaxation::{
    RelaxationNetwork, RelaxationNetworkBuildError, RelaxationNetworkRefusal, RelaxationState,
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
            HeldInternalEffort, RelaxationNetworkMotionControls, RelaxationNetworkMotionStep,
            step_relaxation_network_with_loads,
        },
    },
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, Proposal, Receipt};
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
    step_relaxation_network_with_loads(
        &state.bodies,
        &properties(),
        &links(),
        &state.memory,
        &idle,
        HeldExternalLoads { wrenches, gravity },
        dt_s,
        controls(),
        || false,
    )
    .expect("the owner takes the step")
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
    let mut outputs = vec![f64::NAN; inputs.len()];
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
