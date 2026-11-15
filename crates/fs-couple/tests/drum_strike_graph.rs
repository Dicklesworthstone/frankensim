//! A stick strikes a drum head as a composition on the kernel (bead rc-1109 .6.4).
//!
//! **Graph.** The stick is a free mass and the head two of its modes, each
//! fs-nlmodal's `ModalStructure` over the owner model the prepared impact
//! image compiles (`ModalAcousticTimeModel`, its exact held-force step, the
//! same physical damping rates), each publishing the end of its flow. The
//! contact between them is fs-dcontact's `ContactInterface` with a Hertz law,
//! reading both ends on the kernel's load line. Each structure's heat goes to
//! a sink at the air's temperature. Nothing here is hand-coupled: the kernel
//! steps the three parts.
//!
//! **Against the module.** `LinearImpactSystem` with the same bodies, the
//! same obstacle and no volume (the cavity off) steps the same owner models
//! and bisects the same reaction bracket. The two differ in rounding only:
//! the module recomputes the closure from the modal displacements each step
//! and holds one network compliance, while the interface carries the closure
//! as its own state and the kernel forms each end's compliance from the
//! structure's tangent. Each side accepts the first trial reaction within
//! its own force tolerance, so where the two trial sequences part, the
//! reactions may differ by up to the two tolerances. The bound is bead
//! .6.4's: each coordinate within -120 dB of its full scale over the run.
//!
//! **The cavity on.** The second test adds the resonant head (two modes,
//! not struck) and the sealed air between the heads, the module's volume
//! spring: `(bulk / V) (a . q)^2 / 2` over every mode's signed surface
//! integral, which couples the two heads. The module compiles it into one
//! exact block over both heads. In the graph the two heads and their air are
//! one exact quadratic owner, fs-system's `ExactPortHamiltonianComponent`,
//! over the same `H`, the same viscous rates and the module's interleaved
//! coordinates, its strike port publishing the end of its flow to the same
//! interface. The reordering, beyond the first test's: the module steps its
//! block by fs-phs's `step_render_cached` on the full generalized-force
//! vector (one input per mode, the contact's `-column R` among them), the
//! component by fs-phs's `PreparedExactZoh` on its one input `R` entering
//! through the strike row's weights; the two prepare the same exponential
//! separately and associate the input differently. The bound is the same.
//! The resonant head moves only through the air, so its coordinates within
//! the bound show the air is carried, not dropped.
//!
//! **Allocations.** A counting allocator, per thread so that the other tests
//! running beside it are not counted, shows that no step of either graph
//! after the first allocates, in contact or out of it. The first step is not
//! counted (the cavity-off graph's allocated 8 times in run i439).
//! `LinearImpactSystem` itself allocates on every step (6 times in run i439)
//! and is not gated here.
//!
//! **Not covered.** Not the live drum: its finite-element heads and its
//! radiation are not in these graphs, the air is the module's small-signal
//! volume spring (no cavity wave modes, leakage or thermal loss), and the
//! tip law is a synthetic one. No timing and no claim about the planner's
//! plan.

use core::num::NonZeroUsize;
use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::render::plate::impact::linear::{
    LinearImpactConfig, LinearImpactSystem, VolumeConnection,
};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, VolumeSpring};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_nlmodal::component::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAttachment, ModalStructure,
};
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_ports::component::{PortRole, PortSpec, PortTangent};
use fs_ports::{PortKind, StableId};
use fs_system::phs_exact::ExactPortHamiltonianComponent;
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};

struct Counting;

// Tests run side by side, so each thread counts its own allocations.
thread_local! {
    static TRACKING: Cell<bool> = const { Cell::new(false) };
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}

fn count() {
    // A thread that is being torn down has no cell to read: it is not
    // counting.
    if TRACKING.try_with(Cell::get).unwrap_or(false) {
        let _ = ALLOCATIONS.try_with(|count| count.set(count.get() + 1));
    }
}

/// How many times `f` allocated on this thread.
fn allocations_of(f: impl FnOnce()) -> usize {
    ALLOCATIONS.with(|count| count.set(0));
    TRACKING.with(|tracking| tracking.set(true));
    f();
    TRACKING.with(|tracking| tracking.set(false));
    ALLOCATIONS.with(Cell::get)
}

// SAFETY: this test-only allocator delegates every operation to the system
// allocator unchanged. The only addition is bookkeeping in two thread-local
// cells with constant initialisers, which neither allocates nor touches the
// memory being delegated.
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        // SAFETY: the caller's valid layout is forwarded.
        unsafe { SystemAllocator.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        count();
        // SAFETY: the caller's valid layout is forwarded.
        unsafe { SystemAllocator.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: pointer and layout came from the delegated allocator.
        unsafe { SystemAllocator.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        count();
        // SAFETY: the caller's pointer, layout and size are forwarded.
        unsafe { SystemAllocator.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

/// The air's temperature, of both heat sinks [K].
const AIR_K: f64 = 293.15;
/// The mechanics clock [Hz] and the run's length in samples (10 ms).
const RATE_HZ: u32 = 48_000;
const SAMPLES: u64 = 480;
/// The stick: mass [kg], start [m] (0.2 mm above the head) and speed
/// towards the head [m/s].
const STICK_KG: f64 = 0.04;
const STICK_START_M: f64 = -2.0e-4;
const STICK_SPEED_M_S: f64 = 0.8;
/// The head's two modes: angular frequencies [rad/s], physical viscous
/// rates [1/s], and unit-mass mode shapes at the strike point [1/sqrt(kg)].
const HEAD_OMEGA: [f64; 2] = [1200.0, 1700.0];
const HEAD_RATE: [f64; 2] = [2.0, 2.0];
const HEAD_SHAPE: [f64; 2] = [10.0, 6.0];
/// The resonant head's two modes (not struck): angular frequencies [rad/s]
/// and physical viscous rates [1/s].
const RESONANT_OMEGA: [f64; 2] = [1400.0, 2100.0];
const RESONANT_RATE: [f64; 2] = [2.0, 2.0];
/// The sealed air between the heads: adiabatic bulk modulus `gamma p0` of
/// air at one atmosphere [Pa], the enclosed volume [m^3], the volume
/// coordinate's reference area [m^2], and each head mode's signed surface
/// integral [m^2/sqrt(kg)], the struck head's two then the resonant head's.
/// Synthetic axisymmetric modes of 0.05 m^2 heads, not a measured drum.
const AIR_BULK_PA: f64 = 1.4 * 101_325.0;
const CAVITY_M3: f64 = 0.01;
const REFERENCE_M2: f64 = 0.05;
const CAVITY_AREAS: [f64; 4] = [0.3, -0.13, 0.25, -0.1];
/// The tip's Hertz law [N/m^1.5], exponent 1.5, no gap.
const TIP_K: f64 = 1.0e7;
const TIP_ALPHA: f64 = 1.5;
const PROVENANCE: &str = "synthetic tip test, not identified wood";

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 11,
            kernel_id: 12,
            tile: 13,
            iteration: 14,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn budget() -> ModalAcousticTimeBudget {
    ModalAcousticTimeBudget {
        maximum_total_energy_j: 10.0,
        ..ModalAcousticTimeBudget::audible_reference()
    }
}

/// The module's contact budgets; the interface's search is given the same.
fn contact_config() -> ModalContactConfig {
    ModalContactConfig {
        max_iterations: 100,
        maximum_force_n: 1.0e4,
        maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1.0e-10,
        force_relative_tolerance: 1.0e-9,
    }
}

fn config() -> LinearImpactConfig {
    LinearImpactConfig {
        sample_rate_hz: RATE_HZ,
        max_steps: SAMPLES,
        maximum_generalized_force: 1.0e6,
        component: budget(),
        coupling: ModalCouplingConfig {
            max_modes: 64,
            max_connections: 8,
            max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 10.0,
            maximum_abs_pressure_pa: 1.0e6,
            maximum_abs_connection_force_n: 1.0e5,
            solve_relative_tolerance: 1.0e-10,
            energy_absolute_tolerance_j: 1.0e-10,
            energy_relative_tolerance: 1.0e-8,
        },
        contact: contact_config(),
        multiple: MultiContactConfig {
            max_contacts: 32,
            max_sweeps: 100,
            max_setup_terms: 100_000,
        },
    }
}

/// The stick as the module's free mass, and its unit-mass coefficient
/// `1/sqrt(m)`.
fn stick_body() -> (ImpactBody, f64) {
    ImpactBody::free_mass(STICK_KG, STICK_START_M, STICK_SPEED_M_S).unwrap()
}

fn head_body() -> ImpactBody {
    ImpactBody {
        potential: BodyPotential::Linear(HEAD_OMEGA.to_vec()),
        initial: vec![ModalAcousticState::default(); 2],
        damping_per_s: HEAD_RATE.to_vec(),
    }
}

fn resonant_body() -> ImpactBody {
    ImpactBody {
        potential: BodyPotential::Linear(RESONANT_OMEGA.to_vec()),
        initial: vec![ModalAcousticState::default(); 2],
        damping_per_s: RESONANT_RATE.to_vec(),
    }
}

/// The contact's row over the module's coordinates (stick, head modes): the
/// closure is the stick's displacement towards the head less the head's at
/// the strike point.
fn column(b: f64) -> [f64; 3] {
    [b, -HEAD_SHAPE[0], -HEAD_SHAPE[1]]
}

/// The strike's obstacle over a module of `count` coordinates: the row, the
/// resonant head's coordinates (if any) not in it.
fn obstacle(count: usize) -> Obstacle {
    let (_, b) = stick_body();
    let mut row = column(b).to_vec();
    row.resize(count, 0.0);
    Obstacle::new(
        row,
        1,
        count,
        vec![0.0],
        vec![1.0],
        TIP_K,
        TIP_ALPHA,
        PROVENANCE.into(),
    )
    .unwrap()
}

fn module() -> LinearImpactSystem {
    let (stick, _) = stick_body();
    LinearImpactSystem::new(
        vec![stick, head_body()],
        vec![obstacle(3)],
        vec![],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

/// The module with the resonant head and the air between the heads: one
/// volume spring over all five coordinates, the stick's area zero.
fn cavity_module() -> LinearImpactSystem {
    let (stick, _) = stick_body();
    let mut areas = vec![0.0];
    areas.extend(CAVITY_AREAS);
    let volume = VolumeConnection {
        spring: VolumeSpring {
            bulk_modulus_pa: AIR_BULK_PA,
            volume_m3: CAVITY_M3,
            areas,
        },
        reference_area_m2: REFERENCE_M2,
    };
    LinearImpactSystem::new(
        vec![stick, head_body(), resonant_body()],
        vec![obstacle(5)],
        vec![volume],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

/// The owner model the module compiles for a body: its frequencies with a
/// zero damping ratio, and its physical rates.
fn owner(body: &ImpactBody) -> ModalAcousticTimeModel {
    let BodyPotential::Linear(omega) = &body.potential else {
        unreachable!("both bodies are linear")
    };
    let modes = omega
        .iter()
        .map(|&w| ModalAcousticMode {
            angular_frequency_rad_s: w,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        })
        .collect();
    ModalAcousticTimeModel::try_new_with_damping_rates(
        RATE_HZ,
        modes,
        &body.damping_per_s,
        budget(),
    )
    .unwrap()
}

/// A structure over `body`'s owner model with one attachment of `weights`,
/// telling the end of its flow.
fn structure(body: &ImpactBody, id: &str, weights: Vec<f64>) -> ModalStructure {
    ModalStructure::new(
        owner(body),
        vec![ModalAttachment {
            id: StableId::new(id).unwrap(),
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::MIN,
            weights,
        }],
    )
    .unwrap()
    .publishing_ends()
    .unwrap()
}

/// Both heads and the air between them as one exact quadratic owner over
/// `x = (q, p)` per mode, interleaved, in the module's mode order (the
/// struck head's two, the resonant head's two): `H = p.p/2 + q.K q/2` with
/// `K = diag(omega^2) + k c c^T`, the module's reference-area spring
/// `k = bulk A (A / V)` on `c = areas / A`; `J` pairs each `q` with its `p`;
/// `R` is each mode's viscous rate on its momentum; the one input, the
/// strike force on port `strike`, enters the struck head's momenta with the
/// strike row's weights. It tells the end of the strike flow.
fn heads_and_air() -> ExactPortHamiltonianComponent {
    let omega = [
        HEAD_OMEGA[0],
        HEAD_OMEGA[1],
        RESONANT_OMEGA[0],
        RESONANT_OMEGA[1],
    ];
    let rate = [
        HEAD_RATE[0],
        HEAD_RATE[1],
        RESONANT_RATE[0],
        RESONANT_RATE[1],
    ];
    let weights = [-HEAD_SHAPE[0], -HEAD_SHAPE[1], 0.0, 0.0];
    let k = (AIR_BULK_PA * REFERENCE_M2) * (REFERENCE_M2 / CAVITY_M3);
    let c = CAVITY_AREAS.map(|area| area / REFERENCE_M2);
    let n = 2 * omega.len();
    let mut storage = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut g = vec![0.0; n];
    for mode in 0..omega.len() {
        let (q, p) = (2 * mode, 2 * mode + 1);
        storage[q * n + q] = omega[mode].powi(2);
        storage[p * n + p] = 1.0;
        j[q * n + p] = 1.0;
        j[p * n + q] = -1.0;
        r[p * n + p] = rate[mode];
        g[p] = weights[mode];
    }
    for a in 0..c.len() {
        for b in 0..c.len() {
            storage[(2 * a) * n + 2 * b] += k * c[a] * c[b];
        }
    }
    let system = PortHamiltonian::new(
        n,
        1,
        j,
        r,
        g,
        Box::new(QuadraticStorage::new(storage, n).unwrap()),
    )
    .unwrap();
    let strike = PortSpec::power(
        StableId::new("strike").unwrap(),
        PortKind::MechanicalForceVelocity,
        PortRole::EffortIn,
        NonZeroUsize::MIN,
    )
    .unwrap();
    let budget = LinearSplitBudget {
        max_modes: 64,
        max_omega_dt: 1.0,
    };
    ExactPortHamiltonianComponent::new(system, vec![strike], budget)
        .unwrap()
        .publishing_ends()
        .unwrap()
}

struct Graph {
    system: System,
    stick: NodeId,
    head: NodeId,
    contact: NodeId,
}

/// The strike as a graph, the head part `head` at `head_start` with a port
/// `strike`. Contact end `a` is on the stick and end `b` on the head. A
/// common flow gives each end minus its structure's flow, so the closure
/// rate `flow_a - flow_b` is `b v_stick - phi . v_head` when the stick's
/// weight is `-b` and the head's are the row's own, `-phi`; the structures
/// receive the ends' efforts `R` and `-R`, generalized forces `-b R` and
/// `phi R`, the module's `-column R`.
fn strike_graph<C>(head: C, head_start: C::State) -> Graph
where
    C: PortTangent + 'static,
    C::State: 'static,
{
    let (stick_body, b) = stick_body();
    let row = column(b);
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let stick = builder
        .add_with_tangent(
            "stick",
            structure(&stick_body, "tip", vec![-row[0]]),
            stick_body.initial.clone(),
        )
        .unwrap();
    let head = builder.add_with_tangent("head", head, head_start).unwrap();
    let law = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        TIP_K,
        TIP_ALPHA,
        PROVENANCE.into(),
    )
    .unwrap();
    let limits = contact_config();
    let solve = ContactSolve {
        max_iterations: limits.max_iterations,
        maximum_force_n: limits.maximum_force_n,
        maximum_penetration_m: limits.maximum_penetration_m,
        force_absolute_tolerance_n: limits.force_absolute_tolerance_n,
        force_relative_tolerance: limits.force_relative_tolerance,
    };
    let interface = ContactInterface::new(law, solve)
        .unwrap()
        .reading_ends([true, true])
        .unwrap();
    // The module's start closure, `column . q`, with the head at rest.
    let start = row[0] * stick_body.initial[0].displacement_m_sqrt_kg;
    let contact = builder
        .add_on_line("contact", interface, ContactState::at(start))
        .unwrap();
    let (a, b_end) = (port(&builder, contact, "a"), port(&builder, contact, "b"));
    builder
        .common_flow(port(&builder, stick, "tip"), &[a])
        .unwrap();
    builder
        .common_flow(port(&builder, head, "strike"), &[b_end])
        .unwrap();
    for (node, sink) in [(stick, "stick_air"), (head, "head_air")] {
        builder
            .heat_route(
                port(&builder, node, "heat"),
                &[],
                Some(heat_sink(sink, AIR_K).unwrap()),
                RouteTemperature::Fixed(AIR_K),
            )
            .unwrap();
    }
    Graph {
        system: builder.seal().unwrap(),
        stick,
        head,
        contact,
    }
}

fn graph() -> Graph {
    let (_, b) = stick_body();
    let row = column(b);
    let head_body = head_body();
    strike_graph(
        structure(&head_body, "strike", vec![row[1], row[2]]),
        head_body.initial.clone(),
    )
}

/// The strike on both heads and their air: the head part is the exact owner
/// over both heads, at rest.
fn cavity_graph() -> Graph {
    strike_graph(heads_and_air(), vec![0.0; 8])
}

/// One sample of both runs: the module's interleaved `q, v` per coordinate,
/// the graph's in the same order, and the interface's reaction.
struct Sample {
    module: Vec<f64>,
    graph: Vec<f64>,
    graph_reaction_n: f64,
}

fn interleaved(states: &[ModalAcousticState]) -> Vec<f64> {
    states
        .iter()
        .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
        .collect()
}

/// Both runs, sample by sample; `head` reads the graph's head part as the
/// module's interleaved `q, v` of the coordinates after the stick's.
fn run(
    mut module: LinearImpactSystem,
    graph: Graph,
    head: impl Fn(&System, NodeId) -> Vec<f64>,
) -> Vec<Sample> {
    let dt = 1.0 / f64::from(RATE_HZ);
    let gate = CancelGate::new_clock_free();
    let external = vec![0.0; module.mode_count()];
    let Graph {
        mut system,
        stick,
        head: head_node,
        contact,
    } = graph;
    let mut samples = Vec::new();
    with_cx(|cx| {
        for _ in 0..SAMPLES {
            module.step(&external, &gate).unwrap();
            let report = system.step(dt, cx).unwrap();
            assert_eq!(report.dt_s.to_bits(), dt.to_bits(), "the kernel cut a step");
            let mut g = interleaved(system.state::<Vec<ModalAcousticState>>(stick).unwrap());
            g.extend(head(&system, head_node));
            assert_eq!(g.len(), module.state().len());
            samples.push(Sample {
                module: module.state().to_vec(),
                graph: g,
                graph_reaction_n: system.state::<ContactState>(contact).unwrap().reaction_n,
            });
        }
    });
    samples
}

/// What both tests read from a run: the first sample in contact, the
/// strikes (maximal runs of samples in contact), per coordinate the
/// module's full scale over the run, the graph's largest difference from it
/// and that difference in dB re the full scale, and the peak reaction.
struct Comparison {
    first: usize,
    strikes: usize,
    scale: Vec<f64>,
    worst: Vec<f64>,
    db: Vec<f64>,
    peak_n: f64,
}

/// Before the first touch the stick takes the same owner steps on both
/// sides: its coordinates are equal to the bit. So are the head's when
/// `same_head_owner` (the same owner model on both sides); otherwise the
/// heads rest on both sides, two exact steps of zero, equal as numbers (a
/// zero's sign is the summation order's).
fn compare(samples: &[Sample], same_head_owner: bool) -> Comparison {
    let touching = |s: &Sample| s.graph_reaction_n > 0.0;
    let first = samples
        .iter()
        .position(touching)
        .expect("the stick must strike the head");
    let mut strikes = 0_usize;
    for n in 0..samples.len() {
        if touching(&samples[n]) && (n == 0 || !touching(&samples[n - 1])) {
            strikes += 1;
        }
    }
    for (n, s) in samples[..first].iter().enumerate() {
        for (k, (g, m)) in s.graph.iter().zip(&s.module).enumerate() {
            if k < 2 || same_head_owner {
                assert_eq!(g.to_bits(), m.to_bits(), "sample {n}, coordinate {k}");
            } else {
                assert!(
                    *g == 0.0 && *m == 0.0,
                    "sample {n}, coordinate {k}: {g} {m}"
                );
            }
        }
    }
    let width = samples[0].module.len();
    let mut scale = vec![0.0_f64; width];
    let mut worst = vec![0.0_f64; width];
    for s in samples {
        for k in 0..width {
            scale[k] = scale[k].max(s.module[k].abs());
            worst[k] = worst[k].max((s.graph[k] - s.module[k]).abs());
        }
    }
    let db = worst
        .iter()
        .zip(&scale)
        .map(|(w, s)| 20.0 * (w / s).log10())
        .collect();
    let peak_n = samples
        .iter()
        .map(|s| s.graph_reaction_n)
        .fold(0.0_f64, f64::max);
    Comparison {
        first,
        strikes,
        scale,
        worst,
        db,
        peak_n,
    }
}

#[test]
fn a_stick_strike_on_the_kernel_follows_the_prepared_impact_image() {
    let samples = run(module(), graph(), |system, node| {
        interleaved(system.state::<Vec<ModalAcousticState>>(node).unwrap())
    });
    // The head is lighter at the strike point than the stick (each mode
    // 1/phi^2 = 10 g and 28 g against 40 g), so it springs back and is
    // struck again.
    let c = compare(&samples, true);
    eprintln!(
        "drum strike: first touch at sample {}, {} strikes, peak {:.6e} N, \
         dB re full scale {:?}",
        c.first, c.strikes, c.peak_n, c.db
    );
    // i434: four strikes (samples 11-46, 92-130, 162-201, 235-272), peak
    // 44.8 N, the two sides' reactions within 5.1e-9 relative (the largest at
    // a grazing entry, where the carried and the recomputed closures differ
    // by an ulp).
    assert!(c.strikes >= 2, "{} strikes", c.strikes);
    // Bead .6.4: equal bitwise, or within a stated bound below -120 dB re
    // full scale with the reordering identified (this file's lead). Each
    // coordinate's difference is measured against its own full scale
    // (i434: -279.8 to -296.6 dB).
    for (k, value) in c.db.iter().enumerate() {
        assert!(
            c.worst[k] <= 1.0e-6 * c.scale[k],
            "coordinate {k}: {value} dB re full scale"
        );
    }
}

#[test]
fn a_stick_strike_on_two_heads_and_their_air_follows_the_prepared_impact_image() {
    let samples = run(cavity_module(), cavity_graph(), |system, node| {
        system.state::<Vec<f64>>(node).unwrap().clone()
    });
    let c = compare(&samples, false);
    eprintln!(
        "drum strike, cavity on: first touch at sample {}, {} strikes, peak {:.6e} N, \
         full scale {:?}, dB re full scale {:?}",
        c.first, c.strikes, c.peak_n, c.scale, c.db
    );
    assert!(c.strikes >= 2, "{} strikes", c.strikes);
    // The air carries the strike to the resonant head (coordinates 6 to 9).
    // A fine Runge-Kutta model of the same equations, written before this
    // run: its modes' displacement full scales are 0.37 and 0.067 of the
    // struck head's first mode's (coordinate 2).
    assert!(c.scale[6] > 0.1 * c.scale[2], "{:?}", c.scale);
    assert!(c.scale[8] > 0.01 * c.scale[2], "{:?}", c.scale);
    // Bead .6.4's bound, per coordinate against its own full scale, with the
    // reordering identified (this file's lead).
    for (k, value) in c.db.iter().enumerate() {
        assert!(
            c.worst[k] <= 1.0e-6 * c.scale[k],
            "coordinate {k}: {value} dB re full scale"
        );
    }
}

/// Steps `graph` over the run and asserts that no step after the first
/// allocated on this thread, in contact or out of it, and that the stick
/// touched the head. The first step is not counted.
fn assert_warm_steps_allocate_nothing(graph: Graph) {
    let dt = 1.0 / f64::from(RATE_HZ);
    let Graph {
        mut system,
        contact,
        ..
    } = graph;
    let mut counts = Vec::new();
    let mut touching = 0_usize;
    with_cx(|cx| {
        system.step(dt, cx).unwrap();
        for _ in 1..SAMPLES {
            counts.push(allocations_of(|| {
                system.step(dt, cx).unwrap();
            }));
            if system.state::<ContactState>(contact).unwrap().reaction_n > 0.0 {
                touching += 1;
            }
        }
    });
    if let Some(n) = counts.iter().position(|&count| count > 0) {
        panic!(
            "step {} allocated {} times; each warm step's count: {counts:?}",
            n + 1,
            counts[n]
        );
    }
    assert!(touching > 0, "the stick never touched the head");
}

/// Bead .6.4: a step of the graph allocates nothing once the first step has
/// been taken (8 allocations in run i439), in contact or out of it.
#[test]
fn a_warm_step_of_the_strike_graph_allocates_nothing() {
    assert_warm_steps_allocate_nothing(graph());
}

/// The same for the graph with the air: both heads and the air as one exact
/// owner, whose duration cache and preparation are filled by the first step.
#[test]
fn a_warm_step_of_the_strike_graph_with_the_air_allocates_nothing() {
    assert_warm_steps_allocate_nothing(cavity_graph());
}
