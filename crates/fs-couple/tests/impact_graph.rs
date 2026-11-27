//! fs-couple's linear impact description declared on the composition kernel
//! (`render::plate::impact::graph`, bead rc-1109 .6.4), against
//! `LinearImpactSystem`, which compiles the same description.
//!
//! The description: a stick, a free mass, strikes a head of two modes; a
//! second head of two modes is joined to the first by the sealed air between
//! them, the module's volume spring. The values are those of
//! `tests/drum_strike_graph.rs`'s second test, whose graph is written by
//! hand. Here `ImpactGraph::new` declares it from the module's own inputs.
//!
//! Asserted: before the first touch the stick equal to the bit and the heads
//! zero on both sides; at least two strikes; each coordinate within -120 dB
//! of its full scale over the run; the stored energy within 1e-9 of the
//! start's, and the kernel's summed residual within 1e-12 of it; after a
//! catch (both sides rebuilt from their own head states without the stick)
//! the heads within the same bound; no allocation in a warm step; the
//! refusals the image names; and the kernel's static plan: one looped block
//! solved by a Woodbury tear on the contact's two efforts (2 of 6 unknowns
//! kept), and nothing looped without the stick.
//!
//! Not covered: the live drum (engine-drum's own tests), radiation, several
//! contacts, dampers, timing.

use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::cell::Cell;

use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::plate::impact::graph::ImpactGraph;
use fs_couple::render::plate::impact::linear::{
    LinearImpactConfig, LinearImpactSystem, VolumeConnection,
};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactError, VolumeSpring};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_system::TearStrategy;

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

/// The air's temperature, of every heat sink [K].
const AIR_K: f64 = 293.15;
const RATE_HZ: u32 = 48_000;
const SAMPLES: usize = 480;
const AFTER_CATCH: usize = 240;
const STICK_KG: f64 = 0.04;
const STICK_START_M: f64 = -2.0e-4;
const STICK_SPEED_M_S: f64 = 0.8;
const HEAD_OMEGA: [f64; 2] = [1200.0, 1700.0];
const RESONANT_OMEGA: [f64; 2] = [1400.0, 2100.0];
const RATE_PER_S: [f64; 2] = [2.0, 2.0];
/// The struck head's mode shapes at the strike point [1/sqrt(kg)].
const HEAD_SHAPE: [f64; 2] = [10.0, 6.0];
const AIR_BULK_PA: f64 = 1.4 * 101_325.0;
const CAVITY_M3: f64 = 0.01;
const REFERENCE_M2: f64 = 0.05;
/// Signed surface integrals of the four head modes [m^2/sqrt(kg)].
const CAVITY_AREAS: [f64; 4] = [0.3, -0.13, 0.25, -0.1];
const TIP_K: f64 = 1.0e7;
const TIP_ALPHA: f64 = 1.5;
const PROVENANCE: &str = "synthetic tip test, not identified wood";

fn config() -> LinearImpactConfig {
    let budget = ModalAcousticTimeBudget {
        maximum_total_energy_j: 10.0,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    LinearImpactConfig {
        sample_rate_hz: RATE_HZ,
        max_steps: 1 << 20,
        maximum_generalized_force: 1.0e6,
        component: budget,
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
        contact: ModalContactConfig {
            max_iterations: 100,
            maximum_force_n: 1.0e4,
            maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1.0e-10,
            force_relative_tolerance: 1.0e-9,
        },
        multiple: MultiContactConfig {
            max_contacts: 32,
            max_sweeps: 100,
            max_setup_terms: 100_000,
        },
    }
}

fn head(omega: [f64; 2], initial: &[ModalAcousticState]) -> ImpactBody {
    ImpactBody {
        potential: BodyPotential::Linear(omega.to_vec()),
        initial: initial.to_vec(),
        damping_per_s: RATE_PER_S.to_vec(),
    }
}

/// The two heads from `heads` (four modes), behind the stick if `stick`.
fn bodies(stick: bool, heads: &[ModalAcousticState]) -> Vec<ImpactBody> {
    let mut bodies = Vec::new();
    if stick {
        bodies.push(
            ImpactBody::free_mass(STICK_KG, STICK_START_M, STICK_SPEED_M_S)
                .unwrap()
                .0,
        );
    }
    bodies.push(head(HEAD_OMEGA, &heads[..2]));
    bodies.push(head(RESONANT_OMEGA, &heads[2..]));
    bodies
}

fn volume(stick: bool) -> VolumeConnection {
    let mut areas = if stick { vec![0.0] } else { Vec::new() };
    areas.extend(CAVITY_AREAS);
    VolumeConnection {
        spring: VolumeSpring {
            bulk_modulus_pa: AIR_BULK_PA,
            volume_m3: CAVITY_M3,
            areas,
        },
        reference_area_m2: REFERENCE_M2,
    }
}

/// The strike's obstacle over `row` (stick, struck head, resonant head).
fn obstacle(row: [f64; 5]) -> Obstacle {
    Obstacle::new(
        row.to_vec(),
        1,
        5,
        vec![0.0],
        vec![1.0],
        TIP_K,
        TIP_ALPHA,
        PROVENANCE.into(),
    )
    .unwrap()
}

/// The closure is the stick's displacement towards the head less the
/// struck head's at the strike point.
fn strike_row() -> [f64; 5] {
    let b = 1.0 / STICK_KG.sqrt();
    [b, -HEAD_SHAPE[0], -HEAD_SHAPE[1], 0.0, 0.0]
}

fn at_rest() -> Vec<ModalAcousticState> {
    vec![ModalAcousticState::default(); 4]
}

fn module(stick: bool, heads: &[ModalAcousticState]) -> LinearImpactSystem {
    let contacts = if stick {
        vec![obstacle(strike_row())]
    } else {
        Vec::new()
    };
    LinearImpactSystem::new(
        bodies(stick, heads),
        contacts,
        vec![volume(stick)],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn graph(stick: bool, heads: &[ModalAcousticState]) -> ImpactGraph {
    let contacts = if stick {
        vec![obstacle(strike_row())]
    } else {
        Vec::new()
    };
    ImpactGraph::new(
        bodies(stick, heads),
        contacts,
        vec![volume(stick)],
        config(),
        AIR_K,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

/// The head coordinates of an interleaved state, behind `lead` stick
/// coordinates.
fn heads_of(state: &[f64], lead: usize) -> Vec<ModalAcousticState> {
    state[2 * lead..]
        .as_chunks::<2>()
        .0
        .iter()
        .map(|pair| ModalAcousticState {
            displacement_m_sqrt_kg: pair[0],
            velocity_m_sqrt_kg_per_s: pair[1],
        })
        .collect()
}

/// Per coordinate, the module's full scale over `runs` and the largest
/// difference from it.
fn bound(runs: &[(Vec<f64>, Vec<f64>)]) -> (Vec<f64>, Vec<f64>) {
    let width = runs[0].0.len();
    let mut scale = vec![0.0_f64; width];
    let mut worst = vec![0.0_f64; width];
    for (m, g) in runs {
        for k in 0..width {
            scale[k] = scale[k].max(m[k].abs());
            worst[k] = worst[k].max((g[k] - m[k]).abs());
        }
    }
    (scale, worst)
}

fn db(worst: &[f64], scale: &[f64]) -> Vec<f64> {
    worst
        .iter()
        .zip(scale)
        .map(|(w, s)| 20.0 * (w / s).log10())
        .collect()
}

#[test]
fn the_kernel_image_of_a_strike_on_two_heads_and_their_air_follows_the_module() {
    let gate = CancelGate::new_clock_free();
    let mut m = module(true, &at_rest());
    let mut g = graph(true, &at_rest());
    let external = vec![0.0; m.mode_count()];
    let start_j = m.frame().stored_energy_j;
    assert!(
        (g.frame().stored_energy_j - start_j).abs() <= 1.0e-15 * start_j,
        "{} {}",
        g.frame().stored_energy_j,
        start_j
    );
    let mut runs = Vec::with_capacity(SAMPLES);
    let mut touching = Vec::with_capacity(SAMPLES);
    let (mut stored_gap, mut residual, mut heat_g, mut heat_m) = (0.0_f64, 0.0, 0.0, 0.0);
    for _ in 0..SAMPLES {
        let mf = m.step(&external, &gate).unwrap();
        let gf = g.step(&gate).unwrap();
        stored_gap = stored_gap.max((gf.stored_energy_j - mf.stored_energy_j).abs());
        residual += gf.balance_residual_j;
        heat_g += gf.dissipated_energy_j;
        heat_m += mf.dissipated_energy_j;
        touching.push(g.contact().unwrap().reaction_n > 0.0);
        runs.push((m.state().to_vec(), g.state().to_vec()));
    }
    let first = touching
        .iter()
        .position(|&t| t)
        .expect("the stick must strike the head");
    let strikes = (0..SAMPLES)
        .filter(|&n| touching[n] && (n == 0 || !touching[n - 1]))
        .count();
    for (n, (ms, gs)) in runs[..first].iter().enumerate() {
        for (k, (a, b)) in gs.iter().zip(ms).enumerate() {
            if k < 2 {
                assert_eq!(a.to_bits(), b.to_bits(), "sample {n}, coordinate {k}");
            } else {
                assert!(
                    *a == 0.0 && *b == 0.0,
                    "sample {n}, coordinate {k}: {a} {b}"
                );
            }
        }
    }
    let (scale, worst) = bound(&runs);
    eprintln!(
        "impact graph: first touch at sample {first}, {strikes} strikes, stored energy \
         within {stored_gap:.3e} J of the module's (start {start_j:.6e} J), summed residual \
         {residual:.3e} J, heat {heat_g:.9e} J (module {heat_m:.9e} J), dB re full scale {:?}",
        db(&worst, &scale)
    );
    assert!(strikes >= 2, "{strikes} strikes");
    // The air carries the strike to the resonant head (coordinates 6 to 9).
    assert!(scale[6] > 0.1 * scale[2], "{scale:?}");
    for k in 0..scale.len() {
        assert!(
            worst[k] <= 1.0e-6 * scale[k],
            "coordinate {k}: {:?}",
            db(&worst, &scale)
        );
    }
    assert!(stored_gap <= 1.0e-9 * start_j, "{stored_gap:e}");
    assert!(residual.abs() <= 1.0e-12 * start_j, "{residual:e}");
    assert!(
        (heat_g - heat_m).abs() <= 1.0e-9 * start_j,
        "{heat_g:e} {heat_m:e}"
    );

    // The catch: each side rebuilt from its own head states, the stick gone.
    let mut m = module(false, &heads_of(m.state(), 1));
    let mut g = graph(false, &heads_of(g.state(), 1));
    let external = vec![0.0; m.mode_count()];
    let mut after = Vec::with_capacity(AFTER_CATCH);
    for _ in 0..AFTER_CATCH {
        m.step(&external, &gate).unwrap();
        let frame = g.step(&gate).unwrap();
        assert!(frame.dissipated_energy_j >= 0.0);
        after.push((m.state().to_vec(), g.state().to_vec()));
    }
    let (_, worst_after) = bound(&after);
    eprintln!(
        "after the catch: dB re the strike's full scale {:?}",
        db(&worst_after, &scale[2..])
    );
    for k in 0..worst_after.len() {
        assert!(worst_after[k] <= 1.0e-6 * scale[k + 2], "coordinate {k}");
    }
}

/// Bead .6.4: a warm step of the declared graph allocates nothing, in
/// contact or out of it. The first step prepares the exact owner's
/// exponential and is not counted.
#[test]
fn a_warm_step_of_the_kernel_image_allocates_nothing() {
    let gate = CancelGate::new_clock_free();
    let mut g = graph(true, &at_rest());
    g.step(&gate).unwrap();
    let mut counts = Vec::with_capacity(SAMPLES);
    let mut touching = 0_usize;
    for _ in 1..SAMPLES {
        counts.push(allocations_of(|| {
            g.step(&gate).unwrap();
        }));
        if g.contact().unwrap().reaction_n > 0.0 {
            touching += 1;
        }
    }
    if let Some(n) = counts.iter().position(|&count| count > 0) {
        panic!(
            "step {} allocated {} times; each warm step's count: {counts:?}",
            n + 1,
            counts[n]
        );
    }
    assert!(touching > 0, "the stick never touched the head");
}

/// The kernel's static plan for the image (bead .6.4, Do 1). With the stick,
/// the two parts and the contact form one loop on the contact's load line:
/// six unknowns (the stick's and the head's efforts, the contact's flows and
/// ends), solved by a Woodbury tear that keeps the contact's two efforts.
/// Without the stick nothing loops.
#[test]
fn the_kernel_image_plans_the_strike_as_a_woodbury_tear_on_the_contact() {
    let plan = graph(true, &at_rest()).system().execution_plan();
    assert_eq!((plan.block_count, plan.loop_count), (1, 1));
    let block = &plan.blocks[0];
    assert_eq!(block.nodes, ["part0", "part1", "contact"]);
    assert_eq!(block.tear_strategy, TearStrategy::LoadLineWoodbury);
    assert_eq!(block.load_line_node.as_deref(), Some("contact"));
    assert_eq!((block.total_unknowns, block.kept_unknowns), (6, 2));
    assert_eq!(block.tear_variables, ["contact.a", "contact.b"]);

    let free = graph(false, &at_rest()).system().execution_plan();
    assert_eq!((free.block_count, free.loop_count), (1, 0));
    assert_eq!(free.blocks[0].tear_strategy, TearStrategy::Acyclic);
    assert_eq!(free.blocks[0].total_unknowns, 0);
}

fn refusal(result: Result<ImpactGraph, ImpactError>) -> String {
    match result {
        Ok(_) => panic!("admitted"),
        Err(error) => format!("{error:?}"),
    }
}

#[test]
fn the_kernel_image_refuses_what_it_does_not_declare() {
    let gate = CancelGate::new_clock_free();
    let declare = |contacts: Vec<Obstacle>, sink_k: f64| {
        ImpactGraph::new(
            bodies(true, &at_rest()),
            contacts,
            vec![volume(true)],
            config(),
            sink_k,
            &gate,
        )
    };
    let two = refusal(declare(
        vec![obstacle(strike_row()), obstacle(strike_row())],
        AIR_K,
    ));
    assert!(two.contains("one contact point"), "{two}");
    let lossy = obstacle(strike_row()).with_internal_loss(0.5).unwrap();
    let lossy = refusal(declare(vec![lossy], AIR_K));
    assert!(lossy.contains("internal_loss must be zero"), "{lossy}");
    // A row on the heads alone: both heads are one part (the air joins them).
    let inside = refusal(declare(vec![obstacle([0.0, 1.0, 0.0, -1.0, 0.0])], AIR_K));
    assert!(inside.contains("exactly two parts"), "{inside}");
    let cold = refusal(declare(vec![obstacle(strike_row())], 0.0));
    assert!(cold.contains("positive sink temperature"), "{cold}");
}
