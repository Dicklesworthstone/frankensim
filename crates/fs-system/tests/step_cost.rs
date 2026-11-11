//! What one kernel step costs: time and allocations for a damped chain of
//! masses, against a hand-written loop for the same discrete equations.
//!
//! This is its own test binary so that its counting allocator measures
//! nothing else. The measurement is ignored by default; run it optimised:
//! `cargo test --release -p fs-system --test step_cost -- --ignored --nocapture`.

use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::cell::Cell;
use std::time::Instant;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::{
    Damper, Force, HeatCapacity, Mass, Motion, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, System, SystemBuilder};

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
    line(n, true)
}

/// The chain, or with `damped` false the same line without its dampers:
/// masses, springs and the force, every part of it exactly linear. The body
/// is then alone and receives nothing.
fn line(n: usize, damped: bool) -> (System, Vec<NodeId>, NodeId) {
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
        let mut between = vec![spring];
        if damped {
            let damper = builder
                .add_with_tangent(
                    &format!("damper{index}"),
                    Damper::new(DAMPING_N_S_M).unwrap(),
                    (),
                )
                .unwrap();
            between.push(damper);
            heats.push(builder.port(damper, "heat").unwrap());
        }
        for element in between {
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
    }
    for (mass, members) in masses.iter().zip(&ends) {
        builder
            .common_flow(builder.port(*mass, "p").unwrap(), members)
            .unwrap();
    }
    let warmth = builder.port(body, "heat").unwrap();
    if damped {
        builder.common_effort(warmth, &heats).unwrap();
    } else {
        builder
            .boundary(warmth, Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
    }
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
    damping: f64,
    velocity: Vec<f64>,
    extension: Vec<f64>,
    kelvin: f64,
    // Workspaces of the tridiagonal solve.
    diagonal: Vec<f64>,
    right: Vec<f64>,
    mean: Vec<f64>,
}

impl Hand {
    fn new(n: usize, damped: bool) -> Self {
        Self {
            damping: if damped { DAMPING_N_S_M } else { 0.0 },
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
        let s = 0.5 * STIFFNESS_N_M * DT_S + self.damping;
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
            heat_j += DT_S * self.damping * rate * rate;
            self.velocity[index] = 2.0 * self.mean[index] - self.velocity[index];
        }
        self.kelvin += heat_j / BODY_J_K;
    }
}

#[test]
fn the_chain_from_elements_is_the_hand_written_midpoint_rule() {
    for (n, damped) in [(1_usize, true), (4, true), (9, true), (1, false), (6, false)] {
        let (mut system, masses, body) = line(n, damped);
        let mut hand = Hand::new(n, damped);
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
        // The chain is moving, and with dampers it has warmed its body.
        assert!(hand.velocity[0].abs() > 1.0e-4 && (hand.kelvin > 300.0) == damped);
    }
}

#[test]
fn an_accepted_step_allocates_nothing() {
    // The block's workspaces, its lists and the step's rows are kept
    // between steps. A few steps let them reach their size.
    for (n, damped) in [(1_usize, true), (6, true), (6, false)] {
        let (mut system, _, _) = line(n, damped);
        with_cx(|cx| {
            for _ in 0..10 {
                system.step(DT_S, cx).unwrap();
            }
            let allocations = allocations_of(|| {
                for _ in 0..200 {
                    system.step(DT_S, cx).unwrap();
                }
            });
            assert_eq!(allocations, 0, "{n} masses, damped {damped}");
        });
    }
}

#[test]
fn factors_kept_for_an_affine_block_change_no_bit() {
    // Masses, springs and a force are affine throughout, so the block's
    // interface Jacobian depends on the step length only and its factors
    // are kept from one step to the next. A second copy is made to derive
    // them afresh for every step: before each one it takes a step of
    // another length, which replaces the kept factors, and goes back.
    let n = 6;
    let (mut kept, masses, _) = line(n, false);
    let (mut fresh, _, _) = line(n, false);
    with_cx(|cx| {
        for step in 0..400 {
            let report = kept.step(DT_S, cx).unwrap();
            let before = fresh.checkpoint().unwrap();
            fresh.step(0.37 * DT_S, cx).unwrap();
            fresh.restore(&before).unwrap();
            let again = fresh.step(DT_S, cx).unwrap();
            assert_eq!(report.iterations, again.iterations, "step {step}");
            for mass in &masses {
                let a = kept.state::<Motion>(*mass).unwrap();
                let b = fresh.state::<Motion>(*mass).unwrap();
                assert_eq!(a.position_m.to_bits(), b.position_m.to_bits(), "step {step}");
                assert_eq!(a.velocity_m_s.to_bits(), b.velocity_m_s.to_bits(), "step {step}");
            }
        }
    });
    // The line has moved, and the detour left no trace in the copy's clock
    // or accounts.
    assert!(kept.state::<Motion>(masses[0]).unwrap().position_m.abs() > 1.0e-4);
    assert_eq!(kept.time_s().to_bits(), fresh.time_s().to_bits());
    assert_eq!(
        kept.ledger().energy_residual_j().to_bits(),
        fresh.ledger().energy_residual_j().to_bits()
    );
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
    println!("masses  dampers  kernel ns/step  hand ns/step  ratio  allocations/step  iterations/step");
    for (n, damped, steps) in [
        (1_usize, true, 20_000_u32),
        (4, true, 8_000),
        (16, true, 1_500),
        (48, true, 150),
        (1, false, 20_000),
        (4, false, 8_000),
        (16, false, 3_000),
        (48, false, 600),
    ] {
        let (mut system, _, _) = line(n, damped);
        let mut hand = Hand::new(n, damped);
        let (kernel_ns, allocations, iterations) = with_cx(|cx| {
            // Warm up, so one-time growth is not counted.
            for _ in 0..20 {
                system.step(DT_S, cx).unwrap();
            }
            let mut iterations = 0_u64;
            let mut elapsed = std::time::Duration::ZERO;
            let allocations = allocations_of(|| {
                let start = Instant::now();
                for _ in 0..steps {
                    iterations += u64::from(system.step(DT_S, cx).unwrap().iterations);
                }
                elapsed = start.elapsed();
            });
            (
                elapsed.as_nanos() as f64 / f64::from(steps),
                allocations as f64 / f64::from(steps),
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
            "{n:>6}  {:>7}  {kernel_ns:>14.0}  {hand_ns:>12.1}  {:>5.0}  {allocations:>16.1}  {iterations:>15.2}",
            if damped { "yes" } else { "no" },
            kernel_ns / hand_ns
        );
    }
}
