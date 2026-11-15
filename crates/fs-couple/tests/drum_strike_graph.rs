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
//! **Not covered.** Not the live drum: its finite-element heads, its air
//! cavity (the volume springs between the heads) and its radiation are not
//! in this graph, and the tip law is a synthetic one. No timing, no
//! allocation count, and no claim about the planner's plan.

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::render::plate::impact::linear::{LinearImpactConfig, LinearImpactSystem};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody};
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
use fs_ports::{PortKind, StableId};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};

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

/// The contact's row over the module's coordinates (stick, head modes): the
/// closure is the stick's displacement towards the head less the head's at
/// the strike point.
fn column(b: f64) -> [f64; 3] {
    [b, -HEAD_SHAPE[0], -HEAD_SHAPE[1]]
}

fn module() -> LinearImpactSystem {
    let (stick, b) = stick_body();
    let obstacle = Obstacle::new(
        column(b).to_vec(),
        1,
        3,
        vec![0.0],
        vec![1.0],
        TIP_K,
        TIP_ALPHA,
        PROVENANCE.into(),
    )
    .unwrap();
    LinearImpactSystem::new(
        vec![stick, head_body()],
        vec![obstacle],
        vec![],
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

struct Graph {
    system: System,
    stick: NodeId,
    head: NodeId,
    contact: NodeId,
}

/// The strike as a graph. Contact end `a` is on the stick and end `b` on
/// the head. A common flow gives each end minus its structure's flow, so the
/// closure rate `flow_a - flow_b` is `b v_stick - phi . v_head` when the
/// stick's weight is `-b` and the head's are the row's own, `-phi`; the
/// structures receive the ends' efforts `R` and `-R`, generalized forces
/// `-b R` and `phi R`, the module's `-column R`.
fn graph() -> Graph {
    let (stick_body, b) = stick_body();
    let head_body = head_body();
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
    let head = builder
        .add_with_tangent(
            "head",
            structure(&head_body, "strike", vec![row[1], row[2]]),
            head_body.initial.clone(),
        )
        .unwrap();
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

/// One sample of both runs: the module's interleaved `q, v` per coordinate,
/// the graph's in the same order, and the interface's reaction.
struct Sample {
    module: [f64; 6],
    graph: [f64; 6],
    graph_reaction_n: f64,
}

fn interleaved(states: &[&[ModalAcousticState]]) -> [f64; 6] {
    let mut out = [0.0; 6];
    for (i, s) in states.iter().flat_map(|s| s.iter()).enumerate() {
        out[2 * i] = s.displacement_m_sqrt_kg;
        out[2 * i + 1] = s.velocity_m_sqrt_kg_per_s;
    }
    out
}

fn run() -> Vec<Sample> {
    let dt = 1.0 / f64::from(RATE_HZ);
    let gate = CancelGate::new_clock_free();
    let mut module = module();
    let Graph {
        mut system,
        stick,
        head,
        contact,
    } = graph();
    let mut samples = Vec::new();
    with_cx(|cx| {
        for _ in 0..SAMPLES {
            module.step(&[0.0; 3], &gate).unwrap();
            let report = system.step(dt, cx).unwrap();
            assert_eq!(report.dt_s.to_bits(), dt.to_bits(), "the kernel cut a step");
            let m: [f64; 6] = module.state().try_into().unwrap();
            let g = interleaved(&[
                system
                    .state::<Vec<ModalAcousticState>>(stick)
                    .unwrap()
                    .as_slice(),
                system
                    .state::<Vec<ModalAcousticState>>(head)
                    .unwrap()
                    .as_slice(),
            ]);
            samples.push(Sample {
                module: m,
                graph: g,
                graph_reaction_n: system.state::<ContactState>(contact).unwrap().reaction_n,
            });
        }
    });
    samples
}

#[test]
fn a_stick_strike_on_the_kernel_follows_the_prepared_impact_image() {
    let samples = run();
    let touching = |s: &Sample| s.graph_reaction_n > 0.0;
    let first = samples
        .iter()
        .position(touching)
        .expect("the stick must strike the head");
    // The strikes: maximal runs of samples in contact. The head is lighter at
    // the strike point than the stick (each mode 1/phi^2 = 10 g and 28 g
    // against 40 g), so it springs back and is struck again.
    let mut strikes = 0_usize;
    for n in 0..samples.len() {
        if touching(&samples[n]) && (n == 0 || !touching(&samples[n - 1])) {
            strikes += 1;
        }
    }
    // Before the first touch both sides take the same owner steps on the same
    // states: equal to the bit.
    for (n, s) in samples[..first].iter().enumerate() {
        for k in 0..6 {
            assert_eq!(
                s.graph[k].to_bits(),
                s.module[k].to_bits(),
                "sample {n}, coordinate {k}"
            );
        }
    }
    // Per coordinate: the module's full scale over the run, and the largest
    // difference of the graph from it.
    let mut scale = [0.0_f64; 6];
    let mut worst = [0.0_f64; 6];
    for s in &samples {
        for k in 0..6 {
            scale[k] = scale[k].max(s.module[k].abs());
            worst[k] = worst[k].max((s.graph[k] - s.module[k]).abs());
        }
    }
    let db: Vec<f64> = worst
        .iter()
        .zip(&scale)
        .map(|(w, s)| 20.0 * (w / s).log10())
        .collect();
    let peak = samples
        .iter()
        .map(|s| s.graph_reaction_n)
        .fold(0.0_f64, f64::max);
    eprintln!(
        "drum strike: first touch at sample {first}, {strikes} strikes, peak {peak:.6e} N, \
         dB re full scale {db:?}"
    );
    // i434: four strikes (samples 11-46, 92-130, 162-201, 235-272), peak
    // 44.8 N, the two sides' reactions within 5.1e-9 relative (the largest at
    // a grazing entry, where the carried and the recomputed closures differ
    // by an ulp).
    assert!(strikes >= 2, "{strikes} strikes");
    // Bead .6.4: equal bitwise, or within a stated bound below -120 dB re
    // full scale with the reordering identified (this file's lead). Each
    // coordinate's difference is measured against its own full scale
    // (i434: -279.8 to -296.6 dB).
    for (k, value) in db.iter().enumerate() {
        assert!(
            worst[k] <= 1.0e-6 * scale[k],
            "coordinate {k}: {value} dB re full scale"
        );
    }
}
