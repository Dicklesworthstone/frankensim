//! What one kernel step costs: time and allocations for a damped chain of
//! masses, against a hand-written loop for the same discrete equations.
//!
//! This is its own test binary so that its counting allocator measures
//! nothing else. The measurement is ignored by default; run it optimised:
//! `cargo test --release -p fs-system --test step_cost -- --ignored --nocapture`.

use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::{
    Damper, Force, HeatCapacity, Mass, Motion, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, System, SystemBuilder};

struct Counting;

static TRACKING: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

fn count() {
    if TRACKING.load(Ordering::Relaxed) {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
    }
}

// SAFETY: this test-only allocator delegates every operation to the system
// allocator unchanged. The only addition is relaxed atomic bookkeeping,
// which neither allocates nor touches the memory being delegated.
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

const MASS_KG: f64 = 0.4;
const STIFFNESS_N_M: f64 = 900.0;
const DAMPING_N_S_M: f64 = 1.5;
const FORCE_N: f64 = 2.0;
const BODY_J_K: f64 = 50.0;
const DT_S: f64 = 1.0e-3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 3,
            kernel_id: 9,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// `n` masses in a line, each joined to the next by a spring and a damper
/// and the last to the ground, a constant force on the first, and every
/// damper's heat in one body. `3 n + 2` components in one looped block.
fn chain(n: usize) -> (System, Vec<NodeId>, NodeId) {
    let mut builder = SystemBuilder::new();
    let masses: Vec<NodeId> = (0..n)
        .map(|index| {
            builder
                .add_with_tangent(
                    &format!("mass{index}"),
                    Mass::new(MASS_KG).unwrap(),
                    Motion {
                        position_m: 0.0,
                        velocity_m_s: 0.0,
                    },
                )
                .unwrap()
        })
        .collect();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(FORCE_N).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let mut ends: Vec<Vec<_>> = vec![Vec::new(); n];
    ends[0].push(builder.port(push, "p").unwrap());
    let mut heats = Vec::new();
    for index in 0..n {
        let spring = builder
            .add_with_tangent(
                &format!("spring{index}"),
                Spring::new(STIFFNESS_N_M, 0.0).unwrap(),
                0.0,
            )
            .unwrap();
        let damper = builder
            .add_with_tangent(
                &format!("damper{index}"),
                Damper::new(DAMPING_N_S_M).unwrap(),
                (),
            )
            .unwrap();
        for element in [spring, damper] {
            ends[index].push(builder.port(element, "a").unwrap());
            let far = builder.port(element, "b").unwrap();
            if index + 1 < n {
                ends[index + 1].push(far);
            } else {
                builder
                    .boundary(far, Box::new(|_, _, values| values[0] = 0.0))
                    .unwrap();
            }
        }
        heats.push(builder.port(damper, "heat").unwrap());
    }
    for (mass, members) in masses.iter().zip(&ends) {
        builder
            .common_flow(builder.port(*mass, "p").unwrap(), members)
            .unwrap();
    }
    builder
        .common_effort(builder.port(body, "heat").unwrap(), &heats)
        .unwrap();
    let mut system = builder.seal().unwrap();
    // The dampers' temperature inputs are in the loop and need a start.
    for heat in heats {
        assert!(system.guess(heat, &[300.0]));
    }
    (system, masses, body)
}

/// The same chain by the implicit midpoint rule in its own unknowns: one
/// tridiagonal solve for the masses' mean velocities a step.
struct Hand {
    velocity: Vec<f64>,
    extension: Vec<f64>,
    kelvin: f64,
    // Workspaces of the tridiagonal solve.
    diagonal: Vec<f64>,
    right: Vec<f64>,
    mean: Vec<f64>,
}

impl Hand {
    fn new(n: usize) -> Self {
        Self {
            velocity: vec![0.0; n],
            extension: vec![0.0; n],
            kelvin: 300.0,
            diagonal: vec![0.0; n],
            right: vec![0.0; n],
            mean: vec![0.0; n],
        }
    }

    fn step(&mut self) {
        let n = self.velocity.len();
        // Each spring and damper acts on the mean relative velocity:
        // s (u[i+1] - u[i]) with s = k dt / 2 + c, and u[n] = 0.
        let s = 0.5 * STIFFNESS_N_M * DT_S + DAMPING_N_S_M;
        let inertia = 2.0 * MASS_KG / DT_S;
        for index in 0..n {
            let above = if index > 0 { s } else { 0.0 };
            self.diagonal[index] = inertia + s + above;
            let held = STIFFNESS_N_M * self.extension[index]
                - if index > 0 {
                    STIFFNESS_N_M * self.extension[index - 1]
                } else {
                    0.0
                };
            let pushed = if index == 0 { FORCE_N } else { 0.0 };
            self.right[index] = pushed + held + inertia * self.velocity[index];
        }
        // Thomas algorithm; every off-diagonal entry is -s.
        for index in 1..n {
            let factor = -s / self.diagonal[index - 1];
            self.diagonal[index] += factor * s;
            self.right[index] -= factor * self.right[index - 1];
        }
        for index in (0..n).rev() {
            let below = if index + 1 < n { self.mean[index + 1] } else { 0.0 };
            self.mean[index] = (self.right[index] + s * below) / self.diagonal[index];
        }
        let mut heat_j = 0.0;
        for index in 0..n {
            let below = if index + 1 < n { self.mean[index + 1] } else { 0.0 };
            let rate = below - self.mean[index];
            self.extension[index] += DT_S * rate;
            heat_j += DT_S * DAMPING_N_S_M * rate * rate;
            self.velocity[index] = 2.0 * self.mean[index] - self.velocity[index];
        }
        self.kelvin += heat_j / BODY_J_K;
    }
}

#[test]
fn the_chain_from_elements_is_the_hand_written_midpoint_rule() {
    for n in [1_usize, 4, 9] {
        let (mut system, masses, body) = chain(n);
        let mut hand = Hand::new(n);
        with_cx(|cx| {
            for _ in 0..300 {
                system.step(DT_S, cx).unwrap();
                hand.step();
            }
        });
        for (mass, expected) in masses.iter().zip(&hand.velocity) {
            let velocity = system.state::<Motion>(*mass).unwrap().velocity_m_s;
            assert!(
                (velocity - expected).abs() <= 1.0e-12,
                "{n} masses: {velocity} vs {expected}"
            );
        }
        let kelvin = *system.state::<f64>(body).unwrap();
        assert!(
            (kelvin - hand.kelvin).abs() <= 1.0e-11,
            "{kelvin} vs {}",
            hand.kelvin
        );
        // The chain is moving and has warmed its body.
        assert!(hand.velocity[0].abs() > 1.0e-4 && hand.kelvin > 300.0);
    }
}

#[test]
fn a_loaded_mass_that_has_come_to_rest_keeps_stepping() {
    // Held by a spring against a constant force, the mass moves at the
    // rounding of the two forces that cancel, about 1e-18 m/s, and its
    // damper reports the square of that as heat. The body's heat input is
    // solved for (the damper shares the body's temperature), and that heat
    // has no size of its own to converge relative to: this run used to be
    // refused after about 20 s ("not converged", `body.heat`, a residual of
    // 4.6e-37 W).
    let (mut system, masses, body) = chain(1);
    with_cx(|cx| {
        for step in 0..40_000 {
            system
                .step(DT_S, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
        }
    });
    // At rest the spring holds the force, and what the force did beyond
    // the spring's energy, F x - k x^2 / 2 = F^2 / (2 k), is in the body.
    let motion = system.state::<Motion>(masses[0]).unwrap();
    let rest = FORCE_N / STIFFNESS_N_M;
    assert!((motion.position_m - rest).abs() <= 1.0e-14, "{}", motion.position_m);
    assert!(motion.velocity_m_s.abs() <= 1.0e-15, "{}", motion.velocity_m_s);
    let warmed = system.ledger().nodes[body.index()].stored_change_j.value();
    let expected = FORCE_N * FORCE_N / (2.0 * STIFFNESS_N_M);
    assert!((warmed - expected).abs() <= 1.0e-10 * expected, "{warmed} vs {expected}");
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-12);
}

#[test]
#[ignore = "a measurement: run with --release -- --ignored --nocapture"]
fn cost_of_a_step_against_the_hand_written_loop() {
    println!("masses  components  kernel ns/step  hand ns/step  ratio  allocations/step  iterations/step");
    for (n, steps) in [(1_usize, 20_000_u32), (4, 8_000), (16, 1_500), (48, 150)] {
        let (mut system, _, _) = chain(n);
        let mut hand = Hand::new(n);
        let (kernel_ns, allocations, iterations) = with_cx(|cx| {
            // Warm up, so one-time growth is not counted.
            for _ in 0..20 {
                system.step(DT_S, cx).unwrap();
            }
            let mut iterations = 0_u64;
            ALLOCATIONS.store(0, Ordering::Relaxed);
            TRACKING.store(true, Ordering::Relaxed);
            let start = Instant::now();
            for _ in 0..steps {
                iterations += u64::from(system.step(DT_S, cx).unwrap().iterations);
            }
            let elapsed = start.elapsed();
            TRACKING.store(false, Ordering::Relaxed);
            (
                elapsed.as_nanos() as f64 / f64::from(steps),
                ALLOCATIONS.load(Ordering::Relaxed) as f64 / f64::from(steps),
                iterations as f64 / f64::from(steps),
            )
        });
        let hand_steps = steps * 200;
        let start = Instant::now();
        for _ in 0..hand_steps {
            hand.step();
        }
        let hand_ns = start.elapsed().as_nanos() as f64 / f64::from(hand_steps);
        // Keep the hand loop's result alive.
        assert!(hand.kelvin.is_finite());
        println!(
            "{n:>6}  {:>10}  {kernel_ns:>14.0}  {hand_ns:>12.1}  {:>5.0}  {allocations:>16.1}  {iterations:>15.2}",
            3 * n + 2,
            kernel_ns / hand_ns
        );
    }
}
