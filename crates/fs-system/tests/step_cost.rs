//! What one kernel step costs: time and allocations for a damped chain of
//! masses, against a hand-written loop for the same discrete equations.
//!
//! This is its own test binary so that its counting allocator measures
//! nothing else. The measurement is ignored by default; run it optimised:
//! `cargo test --release -p fs-system --test step_cost -- --ignored --nocapture`.

use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::{Class, Component, PortSpec, PortTangent, Proposal};
use fs_system::elements::{
    Damper, ElementError, Force, HeatCapacity, Mass, Motion, SourceAccount, Spring, Waveform,
};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, SystemError};

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
    let heat = if damped { Heat::Shared } else { Heat::None };
    let (system, masses, bodies) = build(n, heat, None);
    (system, masses, bodies[0])
}

/// Where the dampers of [`build`] give up their heat.
#[derive(Clone, Copy, PartialEq)]
enum Heat {
    /// No dampers.
    None,
    /// Every damper shares the body's temperature (a common-effort junction).
    Shared,
    /// Every damper's heat is routed half to the body and half to a second
    /// one of a thousandth of its heat capacity that starts 20 K warmer, at
    /// their entropy-equivalent temperature: a driver that is not linear.
    Routed,
}

/// A spring that refuses to give its tangent while its flag is set, and is
/// otherwise the element it wraps: exactly linear, with the same ports,
/// steps and slopes.
struct Balky {
    spring: Spring,
    balk: Arc<AtomicBool>,
}

impl Component for Balky {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        self.spring.ports()
    }

    fn class(&self) -> Class {
        self.spring.class()
    }

    fn stored_energy_j(&self, state: &f64) -> Result<f64, ElementError> {
        self.spring.stored_energy_j(state)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.spring.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl PortTangent for Balky {
    fn jvp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        if self.balk.load(Ordering::Relaxed) {
            return Err(ElementError::Input { value: f64::NAN });
        }
        self.spring.jvp(old, dt_s, inputs, d_inputs, d_outputs, cx)
    }

    fn vjp(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        self.spring.vjp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

/// The line of [`chain`], its dampers' heat going where `heat` says, and
/// the body and (with [`Heat::Routed`]) the second body. With `balky`, the
/// spring of that index is a [`Balky`] one on that flag.
fn build(
    n: usize,
    heat: Heat,
    balky: Option<(usize, &Arc<AtomicBool>)>,
) -> (System, Vec<NodeId>, Vec<NodeId>) {
    let damped = heat != Heat::None;
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
        let name = format!("spring{index}");
        let spring = Spring::new(STIFFNESS_N_M, 0.0).unwrap();
        let spring = match balky {
            Some((at, balk)) if at == index => {
                let balk = Arc::clone(balk);
                builder.add_with_tangent(&name, Balky { spring, balk }, 0.0)
            }
            _ => builder.add_with_tangent(&name, spring, 0.0),
        }
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
    let mut bodies = vec![body];
    match heat {
        Heat::None => builder
            .boundary(warmth, Box::new(|_, _, values| values[0] = 0.0))
            .unwrap(),
        Heat::Shared => builder.common_effort(warmth, &heats).unwrap(),
        Heat::Routed => {
            let other = builder
                .add_with_tangent("other", HeatCapacity::new(1.0e-3 * BODY_J_K).unwrap(), 320.0)
                .unwrap();
            let shares = [(warmth, 0.5), (builder.port(other, "heat").unwrap(), 0.5)];
            for &source in &heats {
                builder
                    .heat_route(source, &shares, None, RouteTemperature::SinkEquivalent)
                    .unwrap();
            }
            bodies.push(other);
        }
    }
    let mut system = builder.seal().unwrap();
    // The dampers' temperature inputs are in the loop and need a start.
    for heat in heats {
        assert!(system.guess(heat, &[300.0]));
    }
    (system, masses, bodies)
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

/// Mass [kg], ratio and starting rate [m/s] of mode `index` of [`star`].
fn mode(index: usize) -> (f64, f64, f64) {
    let at = index as f64;
    (
        MASS_KG * (1.0 + 0.01 * at),
        1.0 / (1.0 + 0.05 * at),
        0.1 * (1.0 + (index % 3) as f64),
    )
}

/// `n` masses that each feel one spring's force through a ratio of their
/// own, while the spring is stretched by the sum of their rates through the
/// same ratios: the modes of a body that meet at one contact. The spring's
/// other end is fixed. `n + 1` components in one looped block.
fn star(n: usize) -> (System, Vec<NodeId>) {
    let mut builder = SystemBuilder::new();
    let contact = builder
        .add_with_tangent("contact", Spring::new(STIFFNESS_N_M, 0.0).unwrap(), 0.0)
        .unwrap();
    let mut members = Vec::new();
    let masses = (0..n)
        .map(|index| {
            let (mass_kg, ratio, velocity_m_s) = mode(index);
            let mass = builder
                .add_with_tangent(
                    &format!("mode{index}"),
                    Mass::new(mass_kg).unwrap(),
                    Motion {
                        position_m: 0.0,
                        velocity_m_s,
                    },
                )
                .unwrap();
            members.push((builder.port(mass, "p").unwrap(), ratio));
            mass
        })
        .collect();
    builder
        .common_effort_ratios(builder.port(contact, "a").unwrap(), &members)
        .unwrap();
    builder
        .boundary(
            builder.port(contact, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    (builder.seal().unwrap(), masses)
}

/// The star by the implicit midpoint rule in one unknown, the spring's mean
/// extension: the modes' compliance at the contact is summed once, and a
/// step is one pass for the rate they bring and one to move them.
struct StarHand {
    modes: Vec<(f64, f64)>,
    velocity: Vec<f64>,
    extension: f64,
    compliance: f64,
}

impl StarHand {
    fn new(n: usize) -> Self {
        let modes: Vec<(f64, f64)> = (0..n).map(|index| (mode(index).0, mode(index).1)).collect();
        Self {
            compliance: modes.iter().map(|(mass, ratio)| ratio * ratio / mass).sum(),
            velocity: (0..n).map(|index| mode(index).2).collect(),
            extension: 0.0,
            modes,
        }
    }

    fn step(&mut self) {
        // The spring is stretched at minus the sum of ratio times rate.
        let rate: f64 = self
            .modes
            .iter()
            .zip(&self.velocity)
            .map(|((_, ratio), velocity)| ratio * velocity)
            .sum();
        let mean = (self.extension - 0.5 * DT_S * rate)
            / (1.0 + 0.25 * STIFFNESS_N_M * DT_S * DT_S * self.compliance);
        let force = STIFFNESS_N_M * mean;
        for ((mass, ratio), velocity) in self.modes.iter().zip(&mut self.velocity) {
            let mean_rate = *velocity + 0.5 * DT_S * ratio * force / mass;
            *velocity = 2.0 * mean_rate - *velocity;
        }
        self.extension = 2.0 * mean - self.extension;
    }
}

#[test]
fn the_plan_solves_for_what_closes_the_loops_and_substitutes_the_copies() {
    // The line is one loop. A force sum for each mass and the body's heat
    // sum are solved together; every end rate of a spring, a damper or the
    // force, and every damper's temperature, is a copy and follows.
    for n in [1_usize, 4, 48] {
        assert_eq!(line(n, true).0.plan_loops(), [(6 * n, n + 1)], "{n} masses, damped");
        assert_eq!(line(n, false).0.plan_loops(), [(3 * n, n)], "{n} masses");
    }
    // Modes that meet at one contact: the contact's rate alone, whatever
    // their number. That is the scalar the hand-written loop solves for.
    for n in [1_usize, 5, 200] {
        let (mut system, masses) = star(n);
        assert_eq!(system.plan_loops(), [(n + 1, 1)], "{n} modes");
        let mut hand = StarHand::new(n);
        with_cx(|cx| {
            for _ in 0..300 {
                assert_eq!(system.step(DT_S, cx).unwrap().iterations, 1);
                hand.step();
            }
        });
        for (mass, expected) in masses.iter().zip(&hand.velocity) {
            let velocity = system.state::<Motion>(*mass).unwrap().velocity_m_s;
            assert!(
                (velocity - expected).abs() <= 1.0e-12,
                "{n} modes: {velocity} vs {expected}"
            );
        }
        // The modes have been turned round by the contact.
        assert!(hand.velocity[0] < 0.0 && hand.extension.abs() > 1.0e-4);
        assert!(system.ledger().energy_residual_j().abs() <= 1.0e-12);
    }
}

#[test]
fn the_chain_from_elements_is_the_hand_written_midpoint_rule() {
    for (n, damped) in [(1_usize, true), (4, true), (9, true), (1, false), (6, false)] {
        let (mut system, masses, body) = line(n, damped);
        let mut hand = Hand::new(n, damped);
        with_cx(|cx| {
            for step in 0..300 {
                // Everything but a damper's heat is linear: one update, and
                // with dampers at most one more for the square of the rate.
                let iterations = system.step(DT_S, cx).unwrap().iterations;
                assert!(
                    if damped { iterations <= 2 } else { iterations == 1 },
                    "{n} masses, step {step}: {iterations} iterations"
                );
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
fn heat_routed_at_an_equivalent_temperature_costs_no_iteration_more() {
    // The dampers' heat goes half to each of two bodies at different
    // temperatures, and each damper is given their entropy-equivalent
    // temperature: a driver whose slopes depend on the bodies and are remade
    // for every assembly. With those slopes right, a step takes exactly the
    // iterations of the line whose dampers share one body's temperature.
    // The second body is small so that a slope left out would show: in one
    // update the heat rate moves by about 1e-4 W, which warms that body by
    // about 5e-7 K within the step, against a tolerance of 3e-10 K on the
    // temperature the dampers are given. (An estimate; no defect was
    // planted to see it fail.)
    for n in [1_usize, 4, 9] {
        let (mut routed, masses, bodies) = build(n, Heat::Routed, None);
        let (mut shared, _, _) = build(n, Heat::Shared, None);
        // Kept: the force sums and the two bodies' heat sums; with one
        // damper, its temperature in place of the two sums it feeds.
        let kept = n + if n == 1 { 1 } else { 2 };
        assert_eq!(routed.plan_loops(), [(6 * n + 1, kept)], "{n} masses");
        let mut hand = Hand::new(n, true);
        with_cx(|cx| {
            for step in 0..300 {
                let iterations = routed.step(DT_S, cx).unwrap().iterations;
                assert_eq!(
                    iterations,
                    shared.step(DT_S, cx).unwrap().iterations,
                    "{n} masses, step {step}"
                );
                hand.step();
            }
        });
        for (mass, expected) in masses.iter().zip(&hand.velocity) {
            let velocity = routed.state::<Motion>(*mass).unwrap().velocity_m_s;
            assert!((velocity - expected).abs() <= 1.0e-12, "{n} masses: {velocity} vs {expected}");
        }
        // Each body has half of the heat the hand-written loop made; the
        // second warms a thousand times as much with it.
        let half = 0.5 * (hand.kelvin - 300.0);
        assert!(half > 0.0);
        for (body, (start, rise)) in bodies.iter().zip([(300.0, half), (320.0, 1.0e3 * half)]) {
            let kelvin = *routed.state::<f64>(*body).unwrap();
            assert!((kelvin - start - rise).abs() <= 1.0e-9 * (1.0 + rise), "{kelvin} vs {start} + {rise}");
        }
        assert!(routed.ledger().energy_residual_j().abs() <= 1.0e-12);
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
fn an_affine_assembly_cut_short_by_a_refusal_leaves_no_factors_behind() {
    // An affine block keeps the factors of its interface Jacobian for one
    // step length. Here the fourth spring of a line of six refuses to give
    // its tangent on demand. Before every third step the line is asked for
    // a step of another length while the spring refuses: that assembly has
    // cleared the kept matrix and written the columns before the spring's
    // when it stops. The step of the usual length that follows must be the
    // one a line that was never interrupted takes, to the bit and in the
    // same single iteration: nothing half-written may be used as factors.
    // (By reading the solve, not by a run with it broken: had the kept
    // step length survived the refusal, that step's first update would
    // come from the old factors over entries of another step length and
    // of none, and a second iteration would be needed.)
    let n = 6;
    let (balk, quiet) = (
        Arc::new(AtomicBool::new(false)),
        Arc::new(AtomicBool::new(false)),
    );
    let (mut cut, masses, _) = build(n, Heat::None, Some((3, &balk)));
    let (mut whole, _, _) = build(n, Heat::None, Some((3, &quiet)));
    assert_eq!(cut.plan_loops(), [(3 * n, n)]);
    let mut refusals = 0_u32;
    with_cx(|cx| {
        for step in 0..90 {
            if step % 3 == 1 {
                let before = (cut.time_s().to_bits(), cut.steps());
                balk.store(true, Ordering::Relaxed);
                let refused = cut.step(0.37 * DT_S, cx);
                balk.store(false, Ordering::Relaxed);
                assert!(
                    matches!(&refused, Err(SystemError::Component { node, .. }) if node == "spring3"),
                    "step {step}: {refused:?}"
                );
                assert_eq!((cut.time_s().to_bits(), cut.steps()), before);
                refusals += 1;
            }
            let interrupted = cut.step(DT_S, cx).unwrap();
            let steady = whole.step(DT_S, cx).unwrap();
            assert_eq!((interrupted.iterations, steady.iterations), (1, 1), "step {step}");
            for mass in &masses {
                let a = cut.state::<Motion>(*mass).unwrap();
                let b = whole.state::<Motion>(*mass).unwrap();
                assert_eq!(a.position_m.to_bits(), b.position_m.to_bits(), "step {step}");
                assert_eq!(a.velocity_m_s.to_bits(), b.velocity_m_s.to_bits(), "step {step}");
            }
        }
    });
    // Thirty assemblies were cut short, and the line moved meanwhile.
    assert_eq!(refusals, 30);
    assert!(cut.state::<Motion>(masses[0]).unwrap().position_m.abs() > 1.0e-4);
    assert_eq!(cut.time_s().to_bits(), whole.time_s().to_bits());
    assert_eq!(
        cut.ledger().energy_residual_j().to_bits(),
        whole.ledger().energy_residual_j().to_bits()
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
    println!(" modes  kernel ns/step  hand ns/step  ratio  allocations/step");
    for (n, steps) in [(1_usize, 20_000_u32), (16, 8_000), (200, 2_000)] {
        let (mut system, _) = star(n);
        let mut hand = StarHand::new(n);
        let (kernel_ns, allocations) = with_cx(|cx| {
            for _ in 0..20 {
                system.step(DT_S, cx).unwrap();
            }
            let mut elapsed = std::time::Duration::ZERO;
            let allocations = allocations_of(|| {
                let start = Instant::now();
                for _ in 0..steps {
                    system.step(DT_S, cx).unwrap();
                }
                elapsed = start.elapsed();
            });
            (
                elapsed.as_nanos() as f64 / f64::from(steps),
                allocations as f64 / f64::from(steps),
            )
        });
        let hand_steps = steps * 200;
        let start = Instant::now();
        for _ in 0..hand_steps {
            hand.step();
        }
        let hand_ns = start.elapsed().as_nanos() as f64 / f64::from(hand_steps);
        assert!(hand.extension.is_finite());
        println!(
            "{n:>6}  {kernel_ns:>14.0}  {hand_ns:>12.1}  {:>5.0}  {allocations:>16.1}",
            kernel_ns / hand_ns
        );
    }
}
