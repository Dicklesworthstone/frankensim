//! The non-smooth solve tier: a component on a load line inside a loop.
//!
//! The component is `elements::Slider`, ideal dry friction with constant
//! coefficients. Each test compares a graph with a closed form of the
//! discrete law or of the motion it approximates.

// Exact float comparisons here assert what the law fixes exactly: a held
// or arrested slider has no slip, and a body that never moved has no
// velocity.
#![allow(clippy::float_cmp)]

use std::collections::HashMap;
use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_load_line};
use fs_ports::component::{Class, Component, LoadLine, PortRole, PortSpec, Proposal, Receipt};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    Damper, ElementError, FlowSource, Force, Gyrator, HeatCapacity, Mass, Motion, Slider,
    SliderEnd, SliderRow, Slip, SourceAccount, Spring, Waveform,
};
use fs_system::{
    AdmissionError, Checkpoint, NodeId, RouteTemperature, SolveConfig, System, SystemBuilder,
    SystemError,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
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

const MASS_KG: f64 = 2.0;
const NORMAL_N: f64 = 19.6;
const HOLDS: f64 = 0.5;
const SLIDES: f64 = 0.3;
const BODY_J_K: f64 = 40.0;

fn slider() -> Slider {
    Slider::new(NORMAL_N, HOLDS, SLIDES).unwrap()
}

fn port(builder: &SystemBuilder, node: NodeId, name: &str) -> fs_system::PortRef {
    builder.port(node, name).unwrap()
}

/// The heat [J] the component named `body` has received, from the kernel's
/// ledger. The body's own temperature agrees with it only to the rounding
/// of adding small rises to 300 K, half a unit in the last place of 300 K
/// a step at worst, and that is checked here.
fn received_j(system: &System, body: NodeId) -> f64 {
    let account = system
        .ledger()
        .nodes
        .iter()
        .find(|account| account.name == "body")
        .expect("a component named body");
    let heat_j = account.port_energy_in_j.value();
    let warmed_j = BODY_J_K * (system.state::<f64>(body).unwrap() - 300.0);
    let rounding_j = BODY_J_K * 300.0 * f64::EPSILON * (system.steps() as f64 + 1.0);
    assert!((warmed_j - heat_j).abs() <= rounding_j, "{warmed_j} vs {heat_j}");
    heat_j
}

/// A mass on the ground: a slider between it and a side held at rest, the
/// slider's loss collected in a body at 300 K, and with `push_n` a constant
/// force on the mass. Returns the system, the mass, the slider and the
/// body.
fn on_the_ground(velocity_m_s: f64, push_n: Option<f64>) -> (System, NodeId, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s,
            },
        )
        .unwrap();
    // A member of a common-flow junction receives minus the flow of the
    // mass that sets it, so that is the slip the slider starts with.
    let rubbing = builder
        .add_on_line("slider", slider(), Slip::sliding(-velocity_m_s))
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let mut members = vec![port(&builder, rubbing, "a")];
    if let Some(push_n) = push_n {
        let push = builder
            .add_with_tangent(
                "push",
                Force::new(Waveform::constant(push_n).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        members.push(port(&builder, push, "p"));
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    (system, mass, rubbing, body)
}

/// What a held mass is joined to besides, in [`pushed_by`]: a second mass
/// of 2 kg with the given velocity [m/s], through a spring [N/m] or a
/// damper [N s/m].
#[derive(Clone, Copy)]
enum Beside {
    Spring(f64, f64),
    Damper(f64, f64),
}

/// A mass at rest on the ground under a push that follows `push`, alone or
/// joined to a second mass: the system, the mass and the slider.
fn pushed_by(push: Waveform, beside: Option<Beside>) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let rubbing = builder.add_on_line("slider", slider(), Slip::stuck()).unwrap();
    let pusher = builder
        .add_with_tangent("push", Force::new(push), SourceAccount::default())
        .unwrap();
    let mut members = vec![port(&builder, rubbing, "a"), port(&builder, pusher, "p")];
    if let Some(beside) = beside {
        let (Beside::Spring(_, velocity_m_s) | Beside::Damper(_, velocity_m_s)) = beside;
        let other = builder
            .add_with_tangent(
                "other",
                Mass::new(2.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s,
                },
            )
            .unwrap();
        let link = match beside {
            Beside::Spring(stiffness, _) => builder
                .add_with_tangent("link", Spring::new(stiffness, 0.0).unwrap(), 0.0)
                .unwrap(),
            Beside::Damper(damping, _) => {
                let link = builder
                    .add_with_tangent("link", Damper::new(damping).unwrap(), ())
                    .unwrap();
                builder
                    .boundary(port(&builder, link, "heat"), Box::new(|_, _, values| values[0] = 300.0))
                    .unwrap();
                link
            }
        };
        members.push(port(&builder, link, "a"));
        builder
            .common_flow(port(&builder, other, "p"), &[port(&builder, link, "b")])
            .unwrap();
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), mass, rubbing)
}

/// What a step that held the mass has left, with `c` the mass's answer to
/// a force, `dt / (2 m)`. Returns the larger of the flow at the slider's
/// port and the mass's mean velocity, as a part of the rounding the audit
/// was told for the wire.
fn held_to_rounding(system: &System, mass: NodeId, rubbing: NodeId, c: f64, at: &str) -> f64 {
    let (at_slider, at_mass) = (system.port(rubbing, "a").unwrap(), system.port(mass, "p").unwrap());
    assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "{at}");
    // What the audit is told for the wire, read back and by name: the
    // rounding of the published flow, 16 eps of the released flow and of
    // what the force takes away (c |force| each), and, where the step was
    // taken with slopes, the wire's own, 16 eps of what its target
    // depends on (the released flow once more). No more than the two.
    let force = system.last_output(at_slider)[0].abs();
    let told = system.last_allowance(at_slider)[0];
    let published = 16.0 * f64::EPSILON * 2.0 * c * force;
    assert!(
        told >= 0.95 * published && told <= 1.05 * 1.5 * published,
        "{at}: told {told}, the published flow's rounding {published}"
    );
    // The solve took the wire at half of that, so the flow that is
    // published and the mass's mean velocity are both well inside it.
    let flow = system.last_input(at_slider)[0].abs();
    let mean = system.last_output(at_mass)[0].abs();
    assert!(flow <= told && mean <= told, "{at}: {flow} and {mean} against {told}");
    if told > 0.0 { flow.max(mean) / told } else { 0.0 }
}

#[test]
fn a_mass_held_under_a_push_that_changes_is_taken_at_every_step() {
    // A slider that holds a mass lets no flow cross, while the released
    // flow, c times the push, is not small. The solve judged the flow's
    // wire against the released flow. Under a push that changes from step
    // to step it made one update and stopped with a flow inside its
    // relative tolerance of the released flow; the audit, which is told
    // the rounding of that flow and nothing of the solve's tolerance,
    // refused the step, and every repetition of it:
    //   4.9 N rising to 4.9196 N in 0.4 s: 1 step taken, then connection 0
    //     left -3.2258e-19 J uncancelled (tolerance 6.3977e-20 J);
    //   nothing rising by 2 N/s: 222 steps, -1.4972e-20 J against
    //     5.2685e-22 J;
    //   4 N sin(2 pi 3.14 t): 35 steps, 8.4391e-19 J against 1.7584e-20 J.
    // (Found by the independent review of ef1631e0d2. The landed test held
    // the mass under a constant push, where the first check takes the step
    // with no update.) A flow on the line is judged by what crosses, and
    // taken at half of what the audit is told for it. Every step is taken;
    // what the audit is told is the named rounding; the flow left at the
    // slider's port and the mass's mean velocity are inside it.
    let dt = 1.0e-3_f64;
    let cases: [(&str, Waveform); 3] = [
        ("4.9 N rising to 4.9196 N", Waveform::table(vec![(0.0, 4.9), (0.4, 4.9196)]).unwrap()),
        ("nothing rising by 2 N/s", Waveform::table(vec![(0.0, 0.0), (1.0, 2.0)]).unwrap()),
        ("4 N sin(2 pi 3.14 t)", Waveform::sine(4.0, 2.0 * PI * 0.5, 0.0, 0.0).unwrap()),
    ];
    for (name, push) in cases {
        let (mut system, mass, rubbing) = pushed_by(push, None);
        let (mut iterations, mut most, mut worst) = (0_u32, 0_u32, 0.0_f64);
        with_cx(|cx| {
            for step in 0..400 {
                let report = system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("{name}, step {step}: {error}"));
                iterations += report.iterations;
                most = most.max(report.iterations);
                let at = format!("{name}, step {step}");
                let c = 0.5 * report.dt_s / MASS_KG;
                worst = worst.max(held_to_rounding(&system, mass, rubbing, c, &at));
            }
        });
        assert_eq!(system.steps(), 400);
        eprintln!(
            "held under {name}: 400 steps taken, {iterations} updates ({most} at most in a step); \
             the flow left at the slider's port and the mass's mean velocity are at most \
             {worst:.3} of the rounding the audit is told"
        );
        // Two updates at most: one for the push that changed, one that
        // closes what the differenced slider left of it.
        assert!(most <= 2, "{name}: {most}");
    }
}

#[test]
fn a_mass_held_under_a_slow_push_is_taken_for_as_long_as_it_lasts() {
    // The solve took a line's wire at exactly the allowance the audit was
    // then told for it. The residual it judges is formed from the released
    // flow and what the force takes away, which swamp the mass's own
    // velocity, so a step taken at 0.995 of that allowance was over it by
    // 0.08 % where the audit adds up the junction, and was refused:
    //   3 N rising to 3.3 N in 20 s: 1728 steps taken, then connection 0
    //     left -2.4416e-20 J uncancelled (tolerance 2.4397e-20 J);
    //   6 N + 2 N sin(2 pi 0.05 t): 14939 steps, -4.2666e-20 J against
    //     4.2640e-20 J.
    // (Found by the independent review of 74529fda7f, whose test ran 400
    // steps.) The solve now takes such a wire at half of what the audit
    // is told.
    let dt = 1.0e-3_f64;
    let c = 0.5 * dt / MASS_KG;
    let cases: [(&str, Waveform, u32); 2] = [
        ("3 N rising to 3.3 N in 20 s", Waveform::table(vec![(0.0, 3.0), (20.0, 3.3)]).unwrap(), 2000),
        ("6 N + 2 N sin(2 pi 0.05 t)", Waveform::sine(2.0, 0.05, 0.0, 6.0).unwrap(), 16000),
    ];
    for (name, push, steps) in cases {
        let (mut system, mass, rubbing) = pushed_by(push, None);
        let mut worst = 0.0_f64;
        with_cx(|cx| {
            for step in 0..steps {
                system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("{name}, step {step}: {error}"));
                let at = format!("{name}, step {step}");
                worst = worst.max(held_to_rounding(&system, mass, rubbing, c, &at));
            }
        });
        eprintln!(
            "held under {name}: {steps} steps taken; the flow at the slider's port and the \
             mass's mean velocity are at most {worst:.3} of the rounding the audit is told"
        );
        // Half the allowance, and the rounding of forming the residual.
        assert!(worst <= 0.6, "{name}: {worst}");
    }
}

#[test]
fn a_mass_held_beside_another_that_moves_is_taken_at_every_step() {
    // The solve judged every flow of a block against the largest flow of
    // its kind in the block. A held mass joined to another that moves was
    // then judged against the other's velocity, and the audit of its own
    // junction, which is at rest, refused what the solve had taken (4.9 N
    // rising to 4.9196 N in 0.4 s on the held mass, a second mass of 2 kg):
    //   at 0.2 m/s on a spring of 100 N/m: 1 step, then connection 1 left
    //     1.7830e-19 J uncancelled (tolerance 6.4675e-20 J);
    //   at 0.2 m/s on a damper of 5 N s/m: refused at the first step,
    //     1.4487e-19 J against 6.1732e-20 J;
    //   at rest on the spring: 68 steps, then the second mass's own
    //     junction, 1.9347e-50 J against 6.8756e-52 J.
    // (Found by the independent review of 74529fda7f.) Every wire of a
    // block with a line is now relative to its own size, and one that is
    // not a flow on the line is allowed the rounding of the largest wire
    // of its kind in the block: the second mass's junction carries
    // 1e-43 J beside one that holds 5 N and is solved together with it.
    let dt = 1.0e-3_f64;
    let cases = [
        ("a mass at 0.2 m/s on a spring of 100 N/m", Beside::Spring(100.0, 0.2)),
        ("a mass at 0.2 m/s on a damper of 5 N s/m", Beside::Damper(5.0, 0.2)),
        ("a mass at rest on a spring of 100 N/m", Beside::Spring(100.0, 0.0)),
    ];
    for (name, beside) in cases {
        let push = Waveform::table(vec![(0.0, 4.9), (0.4, 4.9196)]).unwrap();
        let (mut system, mass, rubbing) = pushed_by(push, Some(beside));
        let (mut updates, mut worst) = (0_u32, 0.0_f64);
        with_cx(|cx| {
            for step in 0..2000 {
                let report = system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("{name}, step {step}: {error}"));
                updates += report.iterations;
                let at = format!("{name}, step {step}");
                let c = 0.5 * report.dt_s / MASS_KG;
                worst = worst.max(held_to_rounding(&system, mass, rubbing, c, &at));
            }
        });
        eprintln!(
            "held beside {name}: 2000 steps taken, {updates} updates; the flow at the slider's \
             port and the held mass's mean velocity are at most {worst:.3} of the rounding the \
             audit is told"
        );
    }
}

/// A small generator of the same numbers on every build (splitmix64).
struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn pick(&mut self, of: &[f64]) -> f64 {
        of[(self.next() % of.len() as u64) as usize]
    }

    fn between(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * (self.next() >> 11) as f64 / 9_007_199_254_740_992.0
    }
}

/// A first mass on the ground through the slider, pushed or not, joined by
/// a spring, a damper or both to a second mass that moves or not.
fn drawn_pair(draw: &mut Draw) -> (System, NodeId, f64) {
    let first_kg = draw.pick(&[0.5, 2.0]);
    let second_kg = draw.pick(&[0.5, 2.0, 30.0]);
    let (spring, damper) = match draw.next() % 3 {
        0 => (Some(draw.pick(&[30.0, 100.0, 800.0])), None),
        1 => (None, Some(draw.pick(&[2.0, 5.0, 60.0]))),
        _ => (Some(draw.pick(&[30.0, 800.0])), Some(draw.pick(&[2.0, 60.0]))),
    };
    let push = match draw.next() % 4 {
        0 => None,
        1 => Some(Waveform::constant(draw.between(-8.0, 8.0)).unwrap()),
        2 => Some(Waveform::sine(draw.between(1.0, 8.0), draw.between(0.2, 9.0), 0.0, 0.0).unwrap()),
        _ => {
            let from = draw.between(-6.0, 6.0);
            let to = from + draw.between(-3.0, 3.0);
            Some(Waveform::table(vec![(0.0, from), (draw.between(0.5, 30.0), to)]).unwrap())
        }
    };
    let dt = draw.pick(&[1.0e-3, 1.0e-4]);
    let second_m_s = draw.pick(&[0.0, 0.0, 0.05, 0.5]);
    let at = |velocity_m_s: f64| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let mut builder = SystemBuilder::new();
    let mass = builder.add_with_tangent("mass", Mass::new(first_kg).unwrap(), at(0.0)).unwrap();
    let rubbing = builder.add_on_line("slider", slider(), Slip::stuck()).unwrap();
    let other = builder
        .add_with_tangent("other", Mass::new(second_kg).unwrap(), at(second_m_s))
        .unwrap();
    let mut members = vec![port(&builder, rubbing, "a")];
    let mut theirs = Vec::new();
    if let Some(push) = push {
        let pusher = builder
            .add_with_tangent("push", Force::new(push), SourceAccount::default())
            .unwrap();
        members.push(port(&builder, pusher, "p"));
    }
    if let Some(stiffness) = spring {
        let link = builder
            .add_with_tangent("spring", Spring::new(stiffness, 0.0).unwrap(), 0.0)
            .unwrap();
        members.push(port(&builder, link, "a"));
        theirs.push(port(&builder, link, "b"));
    }
    if let Some(damping) = damper {
        let link = builder
            .add_with_tangent("damper", Damper::new(damping).unwrap(), ())
            .unwrap();
        builder
            .boundary(port(&builder, link, "heat"), Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
        members.push(port(&builder, link, "a"));
        theirs.push(port(&builder, link, "b"));
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    builder.common_flow(port(&builder, other, "p"), &theirs).unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), rubbing, dt)
}

#[test]
fn no_step_of_160_drawn_pairs_on_the_ground_is_refused() {
    // A standing property: a mass on the ground through the slider, joined
    // to a second mass that moves or rests, under a push that is constant,
    // a sine, a ramp or none, 2500 steps each. Whatever the slider does,
    // hold or slide, no step is refused. (On ef1631e0d2 a junction's audit
    // refused 121 of these 160 graphs, on 74529fda7f 76.)
    let mut draw = Draw(20_261_112);
    let (mut held, mut slid) = (0_u32, 0_u32);
    for case in 0..160 {
        let (mut system, rubbing, dt) = drawn_pair(&mut draw);
        let mut let_go = false;
        with_cx(|cx| {
            for step in 0..2500 {
                system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("pair {case}, step {step}: {error}"));
                let_go |= !system.state::<Slip>(rubbing).unwrap().is_stuck();
            }
        });
        if let_go {
            slid += 1;
        } else {
            held += 1;
        }
    }
    eprintln!("160 drawn pairs, 2500 steps each, none refused: held throughout in {held}, sliding at some step in {slid}");
    // The draw covers both.
    assert!(held >= 8 && slid >= 8, "{held} {slid}");
}

/// A block riding a cart, stuck to it by the slider (`a` on the block, `b`
/// on the cart); the cart is pushed with a constant force and held to a
/// wall by a spring. Returns the system, the block, the cart and the
/// slider.
fn riding_a_cart(push_n: f64, anchor_n_m: f64) -> (System, NodeId, NodeId, NodeId) {
    let at_rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let mut builder = SystemBuilder::new();
    let block = builder.add_with_tangent("block", Mass::new(0.5).unwrap(), at_rest).unwrap();
    let cart = builder.add_with_tangent("cart", Mass::new(2.0).unwrap(), at_rest).unwrap();
    let rubbing = builder.add_on_line("slider", slider(), Slip::stuck()).unwrap();
    let pusher = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(push_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let anchor = builder
        .add_with_tangent("anchor", Spring::new(anchor_n_m, 0.0).unwrap(), 0.0)
        .unwrap();
    builder
        .common_flow(port(&builder, block, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    builder
        .common_flow(
            port(&builder, cart, "p"),
            &[port(&builder, rubbing, "b"), port(&builder, pusher, "p"), port(&builder, anchor, "a")],
        )
        .unwrap();
    builder
        .boundary(port(&builder, anchor, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), block, cart, rubbing)
}

#[test]
fn a_block_riding_a_cart_is_carried_through_the_reversals_of_its_friction() {
    // The slider between two bodies that both answer: a block of 0.5 kg
    // riding a cart of 2 kg, the cart pushed with 5 N against a spring of
    // 100 N/m to a wall. The two swing as one, and the friction that
    // carries the block is its mass times their acceleration: at most
    // 1 N, and through zero twice a swing. There the block's junction has
    // next to nothing of its own beside the cart's 10 N of terms. The
    // solve took the block's force wire by the largest force of the block
    // and the audit refused what that left: after 248 steps, connection 0
    // left 3.7002e-15 J uncancelled (tolerance 1.1119e-15 J). (Found by
    // the review of 0e3df343f2, whose drawn pairs all had the slider to
    // the ground.) Every wire of a block with a line is now relative to
    // its own size.
    let dt = 1.0e-3_f64;
    let (mut system, block, cart, rubbing) = riding_a_cart(5.0, 100.0);
    let at_block = system.port(rubbing, "a").unwrap();
    let (mut least_n, mut most_n, mut before) = (f64::INFINITY, 0.0_f64, 0.0_f64);
    with_cx(|cx| {
        for step in 0..6000 {
            system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
            assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "step {step}");
            let (riding, under) = (
                system.state::<Motion>(block).unwrap().velocity_m_s,
                system.state::<Motion>(cart).unwrap().velocity_m_s,
            );
            assert!((riding - under).abs() <= 1.0e-12, "step {step}: {riding} {under}");
            // The friction is what the block's change of velocity needs.
            let force = system.last_output(at_block)[0];
            assert!((force - 0.5 * (riding - before) / dt).abs() <= 1.0e-9, "step {step}: {force}");
            before = riding;
            least_n = least_n.min(force.abs());
            most_n = most_n.max(force.abs());
        }
    });
    // The two swing about the 5 cm at which the spring carries the push,
    // from rest: the friction starts at a fifth of the push, and passes
    // through zero.
    assert!(most_n > 0.99 && most_n <= 1.0 && least_n < 1.0e-3, "{least_n} {most_n}");
    eprintln!(
        "a block riding a cart: 6000 steps taken, stuck in every one; the friction between \
         {least_n:e} and {most_n:e} N"
    );
}

/// A drawn graph for the standing property below: two or three masses in
/// a row, the slider from the first to the ground or between the first
/// two, springs (a quarter of them differenced) and dampers between
/// neighbours and to a wall, a push on one of the masses. Returns the
/// system, its step and the slider.
fn drawn_graph(draw: &mut Draw) -> (System, f64, NodeId) {
    let count = 2 + (draw.next() % 2) as usize;
    let between = draw.next().is_multiple_of(2);
    let dt = draw.pick(&[1.0e-3, 1.0e-3, 1.0e-4, 3.0e-6, 1.0e-2]);
    let mut builder = SystemBuilder::new();
    let (mut speeds, mut masses) = (Vec::new(), Vec::new());
    let mut members: Vec<Vec<fs_system::PortRef>> = Vec::new();
    for name in ["first", "second", "third"].into_iter().take(count) {
        let velocity_m_s = draw.pick(&[0.0, 0.0, 0.0, 0.05, 0.5, -0.3]);
        let kg = draw.pick(&[0.5, 2.0, 30.0]);
        let state = Motion {
            position_m: 0.0,
            velocity_m_s,
        };
        masses.push(builder.add_with_tangent(name, Mass::new(kg).unwrap(), state).unwrap());
        speeds.push(velocity_m_s);
        members.push(Vec::new());
    }
    // The slider's `a` receives minus the first mass's velocity, its `b`
    // minus the second's or nothing: the slip the two start with.
    let slip = if between { speeds[1] - speeds[0] } else { -speeds[0] };
    let start = if slip == 0.0 { Slip::stuck() } else { Slip::sliding(slip) };
    let rubbing = builder.add_on_line("slider", slider(), start).unwrap();
    members[0].push(port(&builder, rubbing, "a"));
    if between {
        members[1].push(port(&builder, rubbing, "b"));
    } else {
        builder
            .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
    }
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    // A link: a damper, or a spring that is differenced one time in four.
    let mut links = 0_u32;
    let mut link = |builder: &mut SystemBuilder, draw: &mut Draw| {
        links += 1;
        if draw.next().is_multiple_of(3) {
            let damper = builder
                .add_with_tangent(
                    &format!("damper{links}"),
                    Damper::new(draw.pick(&[2.0, 5.0, 60.0])).unwrap(),
                    (),
                )
                .unwrap();
            builder
                .boundary(builder.port(damper, "heat").unwrap(), Box::new(|_, _, values| values[0] = 300.0))
                .unwrap();
            damper
        } else {
            let law = Spring::new(draw.pick(&[30.0, 100.0, 800.0]), 0.0).unwrap();
            let name = format!("spring{links}");
            if draw.next().is_multiple_of(4) {
                builder.add(&name, law, 0.0).unwrap()
            } else {
                builder.add_with_tangent(&name, law, 0.0).unwrap()
            }
        }
    };
    // Between neighbours three times in four; to a wall one time in three,
    // and always for a mass that has nothing yet.
    for index in 0..count - 1 {
        if !draw.next().is_multiple_of(4) {
            let node = link(&mut builder, draw);
            members[index].push(port(&builder, node, "a"));
            members[index + 1].push(port(&builder, node, "b"));
        }
    }
    for own in &mut members {
        if own.is_empty() || draw.next().is_multiple_of(3) {
            let node = link(&mut builder, draw);
            own.push(port(&builder, node, "a"));
            builder
                .boundary(port(&builder, node, "b"), Box::new(|_, _, values| values[0] = 0.0))
                .unwrap();
        }
    }
    // A push that is constant, a sine, a ramp or none, up to 14 N: above
    // the 9.8 N the slider holds.
    let wave = match draw.next() % 4 {
        0 => None,
        1 => Some(Waveform::constant(draw.between(-12.0, 12.0)).unwrap()),
        2 => Some(Waveform::sine(draw.between(1.0, 14.0), draw.between(0.2, 9.0), 0.0, 0.0).unwrap()),
        _ => {
            let from = draw.between(-6.0, 6.0);
            let to = from + draw.between(-8.0, 8.0);
            Some(Waveform::table(vec![(0.0, from), (draw.between(0.05, 3.0), to)]).unwrap())
        }
    };
    if let Some(wave) = wave {
        let pushed = (draw.next() % count as u64) as usize;
        let pusher = builder
            .add_with_tangent("push", Force::new(wave), SourceAccount::default())
            .unwrap();
        members[pushed].push(port(&builder, pusher, "p"));
    }
    for (mass, own) in masses.iter().zip(&members) {
        builder.common_flow(port(&builder, *mass, "p"), own).unwrap();
    }
    (builder.seal().unwrap(), dt, rubbing)
}

#[test]
fn no_step_of_the_drawn_carts_and_graphs_is_refused() {
    // Standing properties, the attack families of the review of
    // 0e3df343f2 and b47fe4f2a8. A later reviewer widens them by the
    // counts and the seeds here.
    //
    // 40 drawn carts: the block riding the cart of the test above under a
    // push of 2 to 8 N and a spring of 50 to 400 N/m, 4000 steps of 1 ms
    // each. The slider is between two bodies and its force passes through
    // zero. (On b47fe4f2a8 a junction's audit refused 33 of the 40.)
    let mut draw = Draw(20_261_114);
    for case in 0..40 {
        let (push_n, anchor_n_m) = (draw.between(2.0, 8.0), draw.between(50.0, 400.0));
        let (mut system, _, _, rubbing) = riding_a_cart(push_n, anchor_n_m);
        with_cx(|cx| {
            for step in 0..4000 {
                system
                    .step(1.0e-3, cx)
                    .unwrap_or_else(|error| panic!("cart {case}, step {step}: {error}"));
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "cart {case}, step {step}");
            }
        });
    }
    // 300 drawn graphs of 1500 steps: two or three masses with drawn
    // velocities, the slider to the ground or between two of them, springs
    // (a quarter of them differenced) and dampers, pushes up to 14 N,
    // steps of 3e-6 to 1e-2 s. Holding, sliding, arrests and break-aways:
    // whatever the slider does, no step is refused.
    let mut draw = Draw(20_261_115);
    let (mut slid, mut updates) = (0_u32, 0_u64);
    for case in 0..300 {
        let (mut system, dt, rubbing) = drawn_graph(&mut draw);
        let mut let_go = false;
        with_cx(|cx| {
            for step in 0..1500 {
                let report = system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("graph {case}, step {step}: {error}"));
                updates += u64::from(report.iterations);
                let_go |= !system.state::<Slip>(rubbing).unwrap().is_stuck();
            }
        });
        slid += u32::from(let_go);
    }
    eprintln!(
        "40 drawn carts of 4000 steps and 300 drawn graphs of 1500, none refused: the slider \
         slides at some step in {slid} of the graphs; {updates} updates in their 450000 steps"
    );
    // The draw covers both.
    assert!((60..=240).contains(&slid), "{slid}");
}

/// How the side of [`carried`] is moved.
enum Moved {
    /// A boundary on the slider's `b`, in the sign that port receives.
    Boundary(Box<dyn Fn(f64, f64) -> f64>),
    /// A `FlowSource` bonded to the slider's `b`: the waveform is the
    /// side's velocity as the mass's is counted.
    Source(Waveform),
    /// The same, publishing, with its span wired to the slider's `span_b`.
    Told(Waveform),
    /// A boundary on the slider's `b` that gives minus `mean(t0, dt)`, and
    /// one on its `span_b` that gives `value` at the step's two ends: the
    /// side's velocity as its mover would return it.
    ToldByABoundary {
        mean: Box<dyn Fn(f64, f64) -> f64>,
        value: Box<dyn Fn(f64) -> f64>,
    },
}

/// A mass at rest stuck to a side that is moved, with an optional push on
/// the mass, the slider's loss collected in a body at 300 K. Returns the
/// system, the mass, the slider and the body.
fn carried(side: Moved, push: Option<Waveform>) -> (System, NodeId, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let law = if matches!(side, Moved::Told(_) | Moved::ToldByABoundary { .. }) {
        slider().told(SliderEnd::B)
    } else {
        slider()
    };
    let rubbing = builder.add_on_line("slider", law, Slip::stuck()).unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let mut members = vec![port(&builder, rubbing, "a")];
    if let Some(push) = push {
        let pusher = builder
            .add_with_tangent("push", Force::new(push), SourceAccount::default())
            .unwrap();
        members.push(port(&builder, pusher, "p"));
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    match side {
        Moved::Told(wave) => {
            let source = builder
                .add_with_tangent("side", FlowSource::new(wave).publishing(), SourceAccount::default())
                .unwrap();
            builder
                .bond(port(&builder, rubbing, "b"), port(&builder, source, "p"))
                .unwrap();
            builder
                .signal(port(&builder, source, "span"), port(&builder, rubbing, "span_b"))
                .unwrap();
        }
        Moved::ToldByABoundary { mean, value } => {
            builder
                .boundary(
                    port(&builder, rubbing, "b"),
                    Box::new(move |t0, dt, values| values[0] = -mean(t0, dt)),
                )
                .unwrap();
            builder
                .boundary(
                    port(&builder, rubbing, "span_b"),
                    Box::new(move |t0, dt, values| {
                        values[0] = value(t0);
                        values[1] = value(t0 + dt);
                    }),
                )
                .unwrap();
        }
        Moved::Boundary(given) => builder
            .boundary(
                port(&builder, rubbing, "b"),
                Box::new(move |t0, dt, values| values[0] = given(t0, dt)),
            )
            .unwrap(),
        Moved::Source(wave) => {
            let source = builder
                .add_with_tangent("side", FlowSource::new(wave), SourceAccount::default())
                .unwrap();
            builder
                .bond(port(&builder, rubbing, "b"), port(&builder, source, "p"))
                .unwrap();
        }
    }
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    (system, mass, rubbing, body)
}

// What the slider does not cover, as numbers: the three tests that follow.
// A prescribed end reaches the slider as a step mean: that is what
// `FlowSource` delivers and what a boundary is handed `(t0, dt)` for. A
// mean does not say where the side ends its step, and the slider takes the
// mean for the end: the side is held over the step and jumps between
// steps. For a side at a constant speed, or one that does jump, that is
// the truth. For a side that accelerates it is off by a dt / 2 at every
// step, and none of what these tests pin gets better with a smaller step.
// (Found by the independent review of ef1631e0d2. The slider is to be told
// the speeds a prescribed end starts and ends its step with.)

/// The acceleration [m/s²] of the side in the first two of them: the mass
/// needs `m a` = 8 N to go along, and the slider holds 9.8 N.
const SIDE_M_S2: f64 = 4.0;

/// The speed of a side that rests until `from_s` and accelerates from then
/// on.
fn leaving_rest_at(from_s: f64) -> Waveform {
    Waveform::table(vec![(from_s, 0.0), (from_s + 100.0, 100.0 * SIDE_M_S2)]).unwrap()
}

#[test]
fn a_side_that_leaves_rest_accelerating_breaks_the_pair_away_within_a_few_steps() {
    // The side leaves rest at t = 0. Its first mean is a dt / 2 and the
    // mass is at rest, so the pair that was declared stuck starts with a
    // slip the slider is not told, and its record is off by a dt / 2 from
    // the first step on. The first step holds on the mean, with the 8 N.
    // What follows hangs on the sign of a rounding residue: either the
    // slider goes on holding on the mean for a while, exactly, or it sees
    // a jump of a dt, twice what there is, gives nothing in that step, is
    // asked for 16 N in the next and breaks away. It breaks away within 27
    // steps at every step length tried, the mass is left behind under the
    // sliding force, at 0.2 s it has about 0.73 of the side's 0.8 m/s, and
    // about 0.12 J of heat are booked for a pair that in truth never slips.
    let (needed_n, slide_n) = (MASS_KG * SIDE_M_S2, SLIDES * NORMAL_N);
    assert!(needed_n < 0.82 * HOLDS * NORMAL_N);
    // At step lengths that are not powers of two apart: 8 N or nothing
    // until it breaks away, then the sliding force. The step it breaks away
    // at is not derived, since it hangs on the sign of a rounding residue,
    // and its bound is read off a run: step index 26 at 1/3 ms and at most
    // 20 at the five others on ts2 chain fsz (2026-11-14), the kernel that
    // releases every input of the line's component before its first ask.
    // The kernel before that broke away at step index 4, 2, 16, 16, 6 and 6
    // (chain fsx). The motion and the heat are asserted at every step
    // length before the bound.
    let latest_break = 26_u32;
    let mut broke = Vec::new();
    for dt in [0.9e-3_f64, 0.7e-3, 0.6e-3, 0.3e-3, 1.5e-3, 1.0e-3 / 3.0] {
        let wave = leaving_rest_at(0.0);
        let (mut system, mass, rubbing, body) = carried(Moved::Source(wave.clone()), None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let steps = (0.2 / dt).round() as u32;
        let mut at = None;
        with_cx(|cx| {
            for step in 0..steps {
                system.step(dt, cx).unwrap();
                let force = system.last_output(at_mass)[0];
                if at.is_none() && !system.state::<Slip>(rubbing).unwrap().is_stuck() {
                    at = Some(step);
                }
                if at.is_some() {
                    assert_eq!(force, slide_n, "dt {dt}, step {step}");
                } else {
                    let holds = (force - needed_n).abs() <= 1.0e-6;
                    assert!(holds || force.abs() <= 1.0e-6, "dt {dt}, step {step}: {force}");
                }
            }
        });
        broke.push(at.expect("it breaks away"));
        let part = system.state::<Motion>(mass).unwrap().velocity_m_s / wave.value_at(system.time_s());
        assert!(part > 0.72 && part < 0.76, "dt {dt}: {part}");
        let heat_j = received_j(&system, body);
        assert!(heat_j > 0.11 && heat_j < 0.14, "dt {dt}: {heat_j}");
    }
    eprintln!(
        "a side that leaves rest at 4 m/s2: broken away at step index {broke:?} at 0.9, 0.7, \
         0.6, 0.3, 1.5 and 1/3 ms"
    );
    assert!(broke.iter().all(|&at| at <= latest_break), "{broke:?}");
    // At 1, 0.5 and 0.25 ms, powers of two apart and so the same
    // arithmetic rescaled: 8 N, nothing, and broken away at the third
    // step, through the source and through a boundary that gives the same
    // means.
    for dt in [1.0e-3_f64, 0.5e-3, 0.25e-3] {
        for through_source in [true, false] {
            let wave = leaving_rest_at(0.0);
            let side = if through_source {
                Moved::Source(wave.clone())
            } else {
                let means = wave.clone();
                // The slider's `b` receives the side's velocity with the
                // sign its `a` receives the mass's: minus.
                Moved::Boundary(Box::new(move |t0, dt| -means.mean(t0, dt)))
            };
            let (mut system, mass, rubbing, body) = carried(side, None);
            let at_mass = system.port(rubbing, "a").unwrap();
            let steps = (0.2 / dt).round() as u32;
            with_cx(|cx| {
                for step in 0..steps {
                    system.step(dt, cx).unwrap();
                    let force = system.last_output(at_mass)[0];
                    let stuck = system.state::<Slip>(rubbing).unwrap().is_stuck();
                    match step {
                        0 => assert!((force - needed_n).abs() <= 1.0e-9 && stuck, "{force}"),
                        1 => assert!(force.abs() <= 1.0e-9 && stuck, "{force}"),
                        _ => assert!(force == slide_n && !stuck, "step {step}: {force}"),
                    }
                }
            });
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let side_m_s = wave.value_at(system.time_s());
            // What the three kinds of step leave: 8 N for one step, nothing
            // for one, the sliding force for the rest.
            let left = dt / MASS_KG * (needed_n + slide_n * f64::from(steps - 2));
            assert!((velocity - left).abs() <= 1.0e-12, "{velocity} against {left}");
            assert!((side_m_s - 0.8).abs() <= 1.0e-12 && velocity < 0.74 * side_m_s);
            let heat_j = received_j(&system, body);
            assert!(heat_j > 0.12 && heat_j < 0.13, "{heat_j}");
            if through_source {
                eprintln!(
                    "a side that leaves rest at 4 m/s2, dt {dt}: 8 N, nothing, then broken away; \
                     at 0.2 s the mass has {velocity} m/s of the side's {side_m_s}, heat {heat_j:e} J"
                );
            }
        }
    }
}

#[test]
fn a_side_that_accelerates_after_a_step_at_rest_is_followed_a_mean_behind_with_heat() {
    // The side rests for one step and then accelerates, so the slider sees
    // every change of the mean as the jump it takes it for. The mass is
    // carried, with the right force after the first step (which gets
    // half), but one mean behind, half a step, and every step books the
    // catch of its jump as heat: m a^2 dt^2 / 2, an eighth of it in the
    // first. A pair that never slips loses nothing; this heat falls with
    // the step only in proportion to it.
    let (accel, needed_n) = (SIDE_M_S2, MASS_KG * SIDE_M_S2);
    let mut heats = Vec::new();
    for dt in [1.0e-3_f64, 0.5e-3, 0.25e-3] {
        let wave = leaving_rest_at(dt);
        let (mut system, mass, rubbing, body) = carried(Moved::Source(wave.clone()), None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let steps = (0.2 / dt).round() as u32;
        with_cx(|cx| {
            for step in 0..steps {
                system.step(dt, cx).unwrap();
                let force = system.last_output(at_mass)[0];
                let wanted = match step {
                    0 => 0.0,
                    1 => 0.5 * needed_n,
                    _ => needed_n,
                };
                assert!((force - wanted).abs() <= 1.0e-9, "step {step}: {force}");
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "step {step}");
            }
        });
        let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
        let side_m_s = wave.value_at(system.time_s());
        assert!((side_m_s - velocity - 0.5 * accel * dt).abs() <= 1.0e-12, "{velocity} {side_m_s}");
        let heat_j = received_j(&system, body);
        let booked_j = MASS_KG * accel * accel * dt * dt * (0.125 + 0.5 * f64::from(steps - 2));
        assert!((heat_j - booked_j).abs() <= 1.0e-12 * booked_j, "{heat_j} against {booked_j}");
        heats.push(heat_j);
    }
    eprintln!(
        "a side that rests one step and then accelerates at 4 m/s2: carried half a step behind, \
         heat {heats:?} J at 1, 0.5 and 0.25 ms for a pair that never slips"
    );
    assert!(heats[0] > 1.98 * heats[1] && heats[1] > 1.98 * heats[2]);
}

/// `mean`, the flow a port is given for the step that starts at `t0`, left
/// as it is or moved by a unit in the last place: a hash of the step and
/// of `seed` picks which. Side 0 is the side as it is.
fn last_digit_moved(seed: u64, t0: f64, dt: f64, mean: f64) -> f64 {
    if seed == 0 || mean == 0.0 {
        return mean;
    }
    let step = (t0 / dt).round() as u64;
    let mut hash = step.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ seed.wrapping_mul(0xD1B5_4A32_D192_ED03);
    hash ^= hash >> 29;
    hash = hash.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    hash ^= hash >> 32;
    match hash % 3 {
        0 => mean,
        1 => f64::from_bits(mean.to_bits() + 1),
        _ => f64::from_bits(mean.to_bits() - 1),
    }
}

#[test]
fn a_side_that_swings_from_rest_breaks_away_at_a_reversal_or_not_by_its_last_digit() {
    // A side that swings from rest, 0.02 (1 - cos(w t)) m at 2 Hz, needs
    // at most m 0.02 w^2 = 6.32 N. The slider's record alternates from the
    // first step on (as for the side that leaves rest above). It survives
    // on holding on the mean, which keeps the mass on the side, until at
    // one of the side's reversals it does not, and the pair breaks away.
    // Whether that happens, and at which reversal, hangs on the last digit
    // of the side's means: the side as it is, and 39 more with every mean
    // left or moved by a unit in the last place, are run for 4 s. On the
    // build this was written on 24 of the 40 break away and 16 never do;
    // the side as it is breaks away after its first reversal, at 0.251 s
    // (moved by a `FlowSource` instead of this boundary: at 0.751 s, and
    // at 0.251 s before the solve judged a line flow by what crosses).
    // Asserted is what does not hang on a digit: some break away, and
    // each of those in the step of a reversal or the one after.
    let omega = 2.0 * PI * 2.0;
    assert!(MASS_KG * 0.02 * omega * omega < 0.65 * HOLDS * NORMAL_N);
    let dt = 1.0e-3_f64;
    let (mut broke, mut as_it_is) = (0_u32, None);
    for seed in 0..40_u64 {
        let wave = Waveform::sine(0.02 * omega, 2.0, 0.0, 0.0).unwrap();
        let side = Moved::Boundary(Box::new(move |t0, dt| {
            last_digit_moved(seed, t0, dt, -wave.mean(t0, dt))
        }));
        let (mut system, _, rubbing, _) = carried(side, None);
        let mut at = None;
        with_cx(|cx| {
            for step in 0..4000_u32 {
                system.step(dt, cx).unwrap();
                if !system.state::<Slip>(rubbing).unwrap().is_stuck() {
                    at = Some(f64::from(step) * dt);
                    break;
                }
            }
        });
        if let Some(time) = at {
            broke += 1;
            // In the step that follows a reversal (every 0.25 s), or in
            // the step of it.
            let since = time - 0.25 * (time / 0.25).round();
            assert!(since > -0.5 * dt && since < 1.5 * dt, "side {seed}: {time}");
        }
        if seed == 0 {
            as_it_is = at;
        }
    }
    eprintln!(
        "a side that swings from rest and needs at most 6.32 N of 9.8: as it is, broken away \
         at {as_it_is:?} s; with its means moved by a unit in the last place, {broke} of 40 \
         sides break away within 4 s, each at a reversal"
    );
    assert!(broke >= 1);
}

// What the slider does with a prescribed side that is told its span: the
// tests that follow. The side's mover says the velocities it has at the
// start and at the end of the step, a `FlowSource::publishing` through its
// signal or a boundary on the slider's `span_b`, and the slider no longer
// takes the mean for the end. The span is in the mover's sign; the flow
// the port receives is the mean, with the other sign.

/// A side moved as `wave` and told its span by boundaries: the flow as
/// the port receives it, minus the mean, and the span as the mover returns
/// it.
fn told_by_a_boundary(wave: &Waveform) -> Moved {
    let (means, values) = (wave.clone(), wave.clone());
    Moved::ToldByABoundary {
        mean: Box::new(move |t0, dt| means.mean(t0, dt)),
        value: Box::new(move |t| values.value_at(t)),
    }
}

/// The slider's record against the slip the mass and its side end the step
/// with: the slider's `a` receives minus the mass's velocity and its `b`
/// minus the side's, which ends the step at `side_m_s`.
fn record_is_the_slip(system: &System, mass: NodeId, rubbing: NodeId, side_m_s: f64, at: &str) {
    let slip = side_m_s - system.state::<Motion>(mass).unwrap().velocity_m_s;
    let record = system.state::<Slip>(rubbing).unwrap().slip_m_s();
    assert!((record - slip).abs() <= 1.0e-12, "{at}: the record {record} against the slip {slip}");
}

#[test]
fn a_side_told_its_span_carries_the_mass_as_it_accelerates() {
    // The sides of the three tests above, 4 m/s2 and 8 N needed of the
    // 9.8 N the slider holds, with the span told. Leaving rest at t = 0;
    // resting one step first; and that under a push rising to 15 N over
    // half a requested step, which leaves the slider -7 N to give. At five
    // step lengths, two of them not powers of two from the others: the force
    // is the needed one at every step, the pair never slips, the
    // slider's record is the slip the two end each step with (none), the
    // mass ends on the side, and no heat is booked.
    let needed_n = MASS_KG * SIDE_M_S2;
    let cases = [
        ("leaves rest at t = 0", 0.0_f64, 0.0_f64),
        ("rests one step first", 1.0, 0.0),
        ("rests one step first, 15 N on the mass from then", 1.0, 15.0),
    ];
    for (name, lead_steps, push_n) in cases {
        let (mut worst_n, mut worst_j) = (0.0_f64, 0.0_f64);
        for dt in [1.0e-3_f64, 0.5e-3, 0.25e-3, 0.9e-3, 0.7e-3] {
            let lead = lead_steps * dt;
            let wave = leaving_rest_at(lead);
            // The scheduled push ramp is shorter than the requested step,
            // so the force assertion uses each accepted interval's mean.
            let push = (push_n != 0.0)
                .then(|| Waveform::table(vec![(lead, 0.0), (lead + 0.5 * dt, push_n)]).unwrap());
            let within_n = if push_n == 0.0 { 1.0e-9 } else { 1.0e-7 };
            let (mut system, mass, rubbing, body) = carried(Moved::Told(wave.clone()), push.clone());
            let at_mass = system.port(rubbing, "a").unwrap();
            with_cx(|cx| {
                let mut step = 0;
                while system.time_s() < 0.2 {
                    assert!(step < 2000, "{name}, dt {dt}: clock did not advance");
                    let start_s = system.time_s();
                    let report = system.step(dt, cx).unwrap();
                    let force = system.last_output(at_mass)[0];
                    let drive = push.as_ref().map_or(0.0, |wave| wave.mean(start_s, report.dt_s));
                    let wanted = if start_s < lead { 0.0 } else { needed_n - drive };
                    assert!((force - wanted).abs() <= within_n, "{name}, dt {dt}, step {step}: {force}");
                    assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "{name}, step {step}");
                    let side_m_s = wave.value_at(system.time_s());
                    record_is_the_slip(&system, mass, rubbing, side_m_s, name);
                    worst_n = worst_n.max((force - wanted).abs());
                    step += 1;
                }
            });
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let side_m_s = wave.value_at(system.time_s());
            assert!((velocity - side_m_s).abs() <= 1.0e-12, "{name}, dt {dt}: {velocity} {side_m_s}");
            let heat_j = received_j(&system, body);
            assert!(heat_j.abs() <= 1.0e-12, "{name}, dt {dt}: {heat_j}");
            worst_j = worst_j.max(heat_j.abs());
        }
        eprintln!(
            "a side told its span that {name}: the needed force to {worst_n:e} N at every step of \
             five step lengths, {worst_j:e} J of heat at most"
        );
    }
}

#[test]
fn a_span_from_the_source_and_one_from_a_boundary_give_the_same_to_the_bit() {
    // The sign of the span, pinned. The same side is told its span twice:
    // by its `FlowSource`, whose signal has the waveform's own values
    // while the port receives minus its mean; and by a boundary on the
    // port that gives minus the mean, with a boundary on the span that
    // gives the waveform's values. A side that swings and reverses, and a
    // side that accelerates under a push on the mass.
    //
    // What the slider receives, the flow at its port and the two
    // velocities of its span, is the same to the bit at every step. The
    // forces are not: the source is a component of the block, so the
    // solve has other unknowns and forms the same numbers in another
    // order, and a held slider's force is the mass's velocity change over
    // c, which magnifies a unit in the last place of a velocity of some
    // m/s to 1e-12 N. They agree to that rounding, as do the mass's
    // velocities (measured: 4.0e-10 N and 1.0e-13 m/s at most after 1000
    // steps under the push), and the slider is stuck or not in the same
    // steps.
    let omega = 2.0 * PI * 2.0;
    let cases = [
        (Waveform::sine(0.02 * omega, 2.0, 0.0, 0.0).unwrap(), None),
        (leaving_rest_at(0.0), Some(Waveform::constant(9.0).unwrap())),
    ];
    for (wave, push) in cases {
        let (mut by_source, mass, rubbing, _) = carried(Moved::Told(wave.clone()), push.clone());
        let (mut by_boundary, _, _, _) = carried(told_by_a_boundary(&wave), push);
        let at_mass = by_source.port(rubbing, "a").unwrap();
        let (at_side, at_span) = (
            by_source.port(rubbing, "b").unwrap(),
            by_source.port(rubbing, "span_b").unwrap(),
        );
        let (mut apart_n, mut apart_m_s, mut differ) = (0.0_f64, 0.0_f64, 0_u32);
        with_cx(|cx| {
            for step in 0..1000 {
                by_source.step(1.0e-3, cx).unwrap();
                by_boundary.step(1.0e-3, cx).unwrap();
                // What the slider receives is the same to the bit: the
                // flow at its port and the two velocities of the span.
                for port in [at_side, at_span] {
                    let (ours, theirs) = (by_source.last_input(port), by_boundary.last_input(port));
                    assert_eq!(ours.len(), theirs.len());
                    for (a, b) in ours.iter().zip(theirs) {
                        assert_eq!(a.to_bits(), b.to_bits(), "step {step}: {a} {b}");
                    }
                }
                let forces = (by_source.last_output(at_mass)[0], by_boundary.last_output(at_mass)[0]);
                let velocities = (
                    by_source.state::<Motion>(mass).unwrap().velocity_m_s,
                    by_boundary.state::<Motion>(mass).unwrap().velocity_m_s,
                );
                differ += u32::from(forces.0.to_bits() != forces.1.to_bits());
                apart_n = apart_n.max((forces.0 - forces.1).abs());
                apart_m_s = apart_m_s.max((velocities.0 - velocities.1).abs());
                let records = (
                    *by_source.state::<Slip>(rubbing).unwrap(),
                    *by_boundary.state::<Slip>(rubbing).unwrap(),
                );
                assert_eq!(records.0.is_stuck(), records.1.is_stuck(), "step {step}");
            }
        });
        eprintln!(
            "a span from the source and one from a boundary: the slider receives the same to the \
             bit; the forces differ in {differ} of 1000 steps, by {apart_n:e} N at most, the \
             mass's velocity by {apart_m_s:e} m/s at most"
        );
        assert!(apart_n <= 1.0e-9 && apart_m_s <= 1.0e-12, "{apart_n} {apart_m_s}");
    }
}

#[test]
fn a_side_that_swings_told_its_span_carries_the_mass_through_its_reversals() {
    // The swinging side of the tests of what is not covered: 0.02 (1 -
    // cos(w t)) m at 2 Hz, at most 6.32 N needed. Told its span, the pair
    // never slips. The force is m times the side's acceleration at the
    // middle of the step to first order in the step (what the midpoint
    // relation makes of a speed that is curved within the step), the
    // mass's end speed is the side's to second order, and the heat booked
    // for a pair that never slips is second order.
    let omega = 2.0 * PI * 2.0;
    let needed_n = MASS_KG * 0.02 * omega * omega;
    let wave = Waveform::sine(0.02 * omega, 2.0, 0.0, 0.0).unwrap();
    let (mut off_n, mut off_m_s, mut heats_j) = (Vec::new(), Vec::new(), Vec::new());
    for dt in [1.0e-3_f64, 0.5e-3, 0.25e-3] {
        let (mut system, mass, rubbing, body) = carried(Moved::Told(wave.clone()), None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let steps = (1.0 / dt).round() as u32;
        let (mut worst_n, mut worst_m_s) = (0.0_f64, 0.0_f64);
        with_cx(|cx| {
            for step in 0..steps {
                system.step(dt, cx).unwrap();
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "dt {dt}, step {step}");
                let wanted = needed_n * (omega * (f64::from(step) + 0.5) * dt).cos();
                worst_n = worst_n.max((system.last_output(at_mass)[0] - wanted).abs());
                let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
                let side_m_s = wave.value_at(system.time_s());
                record_is_the_slip(&system, mass, rubbing, side_m_s, "the swinging side");
                worst_m_s = worst_m_s.max((velocity - side_m_s).abs());
            }
        });
        off_n.push(worst_n);
        off_m_s.push(worst_m_s);
        heats_j.push(received_j(&system, body));
    }
    eprintln!(
        "a side that swings, told its span, at 1, 0.5 and 0.25 ms: never slips; force off m a at \
         mid-step by {off_n:?} N, mass off the side's end speed by {off_m_s:?} m/s, heat \
         {heats_j:?} J"
    );
    assert!(off_n[0] <= 2.0e-2 && off_m_s[0] <= 1.0e-5 && heats_j[0] <= 5.0e-6);
    for finer in 1..3 {
        assert!(off_n[finer] <= 0.55 * off_n[finer - 1], "{off_n:?}");
        assert!(off_m_s[finer] <= 0.3 * off_m_s[finer - 1], "{off_m_s:?}");
        assert!(heats_j[finer] >= 0.0 && heats_j[finer] <= 0.3 * heats_j[finer - 1], "{heats_j:?}");
    }

    // The 40 sides of the test above, the side as it is and 39 with every
    // mean left or moved by a unit in the last place, of which 24 break
    // away when they are taken as held: told their span, none does.
    let dt = 1.0e-3_f64;
    for seed in 0..40_u64 {
        let (means, values) = (wave.clone(), wave.clone());
        let side = Moved::ToldByABoundary {
            mean: Box::new(move |t0, dt| -last_digit_moved(seed, t0, dt, -means.mean(t0, dt))),
            value: Box::new(move |t| values.value_at(t)),
        };
        let (mut system, _, rubbing, _) = carried(side, None);
        with_cx(|cx| {
            for step in 0..4000_u32 {
                system.step(dt, cx).unwrap();
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "side {seed}, step {step}");
            }
        });
    }
}

#[test]
fn a_side_that_speeds_up_told_its_span_breaks_away_in_the_step_that_reaches_the_limit() {
    // The side of the speeding-up test: one step at rest, then the speed
    // (rate / 2) u^2, which needs m rate u. Told its span, the slider
    // holds on the mean, with m rate u at the middle of the step to first
    // order (0.016 N at 1 ms), and breaks away where that passes 9.8 N:
    // in the step that contains u* = mu_s N / (m rate) or in the one
    // after, by that first order. Taken as held, the side is followed one
    // mean behind and breaks away in the first step that starts after u*
    // (the test above).
    let rate = 48.0_f64;
    let (hold_n, slide_n) = (HOLDS * NORMAL_N, SLIDES * NORMAL_N);
    let reaches_s = hold_n / (MASS_KG * rate);
    let mut broke = Vec::new();
    for dt in [1.0e-3_f64, 0.5e-3, 0.25e-3] {
        let side = Moved::ToldByABoundary {
            // The exact mean of the speed over the step that starts at t0.
            mean: Box::new(move |t0, _| {
                let u = t0 - dt;
                if u < -0.5 * dt {
                    0.0
                } else {
                    0.5 * rate * (u * u + u * dt + dt * dt / 3.0)
                }
            }),
            value: Box::new(move |t| {
                let u = (t - dt).max(0.0);
                0.5 * rate * u * u
            }),
        };
        let (mut system, _, rubbing, _) = carried(side, None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let mut at = None;
        with_cx(|cx| {
            for step in 0..(0.12 / dt).round() as u32 {
                system.step(dt, cx).unwrap();
                let force = system.last_output(at_mass)[0];
                let stuck = system.state::<Slip>(rubbing).unwrap().is_stuck();
                if at.is_none() && !stuck {
                    at = Some(step);
                }
                if at.is_some() {
                    assert_eq!(force, slide_n, "dt {dt}, step {step}");
                } else if step >= 1 {
                    let middle = (f64::from(step) - 0.5) * dt;
                    let wanted = MASS_KG * rate * middle;
                    assert!(
                        (force - wanted).abs() <= 2.0e-2 * dt / 1.0e-3,
                        "dt {dt}, step {step}: {force} against {wanted}"
                    );
                }
            }
        });
        let at = at.expect("it breaks away");
        // The step that contains u*, counted after the step at rest, or
        // the next.
        let contains = 1 + (reaches_s / dt).floor() as u32;
        assert!(at == contains || at == contains + 1, "dt {dt}: {at} against {contains}");
        broke.push(at);
    }
    eprintln!("a side that speeds up, told its span: broken away at step {broke:?} at 1, 0.5 and 0.25 ms");
}

#[test]
fn a_change_of_slope_inside_a_step_is_made_up_in_the_next_and_can_break_the_pair_away() {
    // A plain boundary cannot advertise its events. The limit of its
    // told span: a side that rests and then accelerates,
    // the change in the middle of a step. In that step the side's mean is
    // less than the mean of its two ends, and the force that would end
    // the step on the side would move the mass faster than the side on
    // the mean while pulling it along: a negative loss, which the slider
    // does not book. It holds on the mean, with a quarter of the needed
    // force where a half would end on the side, and the mass ends the
    // step behind. The next step catches that up, with a quarter of the
    // needed force more than it needs: m times a quarter of the change of
    // acceleration, for one step. Then the pair is on the side again.
    let dt = 1.0e-3_f64;
    let from_s = 2.5 * dt;
    let run = |accel: f64, steps: u32| {
        let wave = Waveform::table(vec![(from_s, 0.0), (from_s + 100.0, 100.0 * accel)]).unwrap();
        let (mut system, mass, rubbing, _) = carried(told_by_a_boundary(&wave), None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let mut seen = Vec::new();
        with_cx(|cx| {
            for _ in 0..steps {
                system.step(dt, cx).unwrap();
                let side_m_s = wave.value_at(system.time_s());
                record_is_the_slip(&system, mass, rubbing, side_m_s, "a change of slope inside a step");
                let behind = side_m_s - system.state::<Motion>(mass).unwrap().velocity_m_s;
                seen.push((
                    system.last_output(at_mass)[0],
                    system.state::<Slip>(rubbing).unwrap().is_stuck(),
                    behind,
                ));
            }
        });
        seen
    };
    // 3.5 m/s2, 7 N needed: 1.75 N in the step of the change, 8.75 N in
    // the next, 7 N from then on, stuck throughout.
    let seen = run(3.5, 40);
    for (step, (force, stuck, behind)) in seen.iter().enumerate() {
        let wanted = match step {
            0 | 1 => 0.0,
            2 => 1.75,
            3 => 8.75,
            _ => 7.0,
        };
        assert!((force - wanted).abs() <= 1.0e-9 && *stuck, "step {step}: {force}");
        // Behind by a quarter of a step's gain after the step of the
        // change, on the side after the next.
        let lag = if step == 2 { 0.25 * 3.5 * dt } else { 0.0 };
        assert!((behind - lag).abs() <= 1.0e-12, "step {step}: {behind}");
    }
    // 4 m/s2, 8 N needed of the 9.8 N the slider holds: the 10 N the next
    // step wants are beyond the limit, and the pair breaks away although
    // the side never asks for more than 8 N. Stated, not mended: the
    // force that ends the step of the change on the side has a negative
    // loss, and the slider is told no more than the step's two ends.
    let seen = run(4.0, 4);
    assert!((seen[2].0 - 2.0).abs() <= 1.0e-9 && seen[2].1, "{:?}", seen[2]);
    assert!(seen[3].0 == SLIDES * NORMAL_N && !seen[3].1, "{:?}", seen[3]);
}

#[test]
fn source_knots_keep_a_belt_stuck_across_changes_of_acceleration() {
    // Without scheduled knots, the change at 2.5 ms fell inside a 1 ms
    // step. The 2 kg body then needed a spurious 10 N catch-up force and
    // broke away, though the physical force is 8 N and the limit 9.8 N.
    let wave = Waveform::table(vec![
        (0.0025, 0.0),
        (0.1025, 0.4),
        (0.2025, 0.0),
        (0.3025, -0.4),
        (0.4025, 0.0),
    ])
    .unwrap();
    for step_cap in [0.001_f64, 0.003, 0.007] {
        let (mut system, mass, rubbing, body) = carried(Moved::Told(wave.clone()), None);
        let at_mass = system.port(rubbing, "a").unwrap();
        let mut landed = [false; 5];
        with_cx(|cx| {
            let mut steps = 0;
            while system.time_s() < 0.45 {
                assert!(steps < 1000, "clock did not advance at cap {step_cap}");
                let start = system.time_s();
                let requested = step_cap.min(0.45 - start);
                let report = system.step(requested, cx).unwrap();
                assert!(report.dt_s <= requested);
                assert!(
                    report.dt_s > 1.0e-6,
                    "roundoff-only step from {start} at cap {step_cap}: {}",
                    report.dt_s
                );
                for (seen, knot) in landed
                    .iter_mut()
                    .zip([0.0025, 0.1025, 0.2025, 0.3025, 0.4025])
                {
                    assert!(
                        !(start < knot && report.time_s > knot),
                        "crossed {knot} at {step_cap}"
                    );
                    *seen |= report.time_s == knot;
                }
                let speed = wave.value_at(report.time_s);
                let needed = MASS_KG * (speed - wave.value_at(start)) / report.dt_s;
                let force = system.last_output(at_mass)[0];
                assert!(
                    (force - needed).abs() < 1.0e-8,
                    "{start}: {force} vs {needed}"
                );
                assert!(
                    system.state::<Slip>(rubbing).unwrap().is_stuck(),
                    "{start} at {step_cap}"
                );
                record_is_the_slip(&system, mass, rubbing, speed, "source knot");
                assert!(
                    (system.state::<Motion>(mass).unwrap().velocity_m_s - speed).abs() < 1.0e-12
                );
                steps += 1;
            }
        });
        assert!(landed.into_iter().all(|seen| seen));
        assert!(received_j(&system, body).abs() < 1.0e-11);
        assert!(system.ledger().energy_residual_j().abs() < 1.0e-10);
    }
}

/// One step of [`the_friction_belt`]: when it ended, the force the slider
/// gave the block, the slider's record of the slip, whether the slider
/// ended stuck, and the block's own slip, the belt's speed less the
/// block's.
#[derive(Debug, Clone, Copy)]
struct OnTheBelt {
    time_s: f64,
    force_n: f64,
    record_m_s: f64,
    stuck: bool,
    slip_m_s: f64,
}

/// fs-cli's `examples/system/friction_belt.json`, built in the order and
/// with the joins its `system` command builds it: a 2 kg block at rest, a
/// belt whose speed is tabulated (0, 0) (0.25, 1) (0.5, `last_m_s`) (1.5,
/// `last_m_s`) with its span told to the slider's `b`, a slider under 20 N
/// with coefficients 0.5 and 0.2 (it holds 10 N and slides at 4 N), and
/// its loss in a lining of 50 J/K at 300 K. Run for 1.5 s in steps of 5 ms
/// as the command runs it; every step is returned.
///
/// The block follows the first ramp with 8 N, held, and breaks away on the
/// second, which needs 16 N at a last speed of 3 m/s. It then slides at
/// 2 m/s2 from 1 m/s at 0.25 s and catches the belt at 0.25 + (last - 1) / 2
/// s. A step that starts sliding with slip `s0` and no other force on the
/// block is ended at no slip by the force `s0 / (2 c)`, where c = dt / (2 m)
/// = 1.25e-3 m/s per N: 400 s0.
fn the_friction_belt(last_m_s: f64) -> Vec<OnTheBelt> {
    let wave = Waveform::table(vec![
        (0.0, 0.0),
        (0.25, 1.0),
        (0.5, last_m_s),
        (1.5, last_m_s),
    ])
    .unwrap();
    let mut builder = SystemBuilder::new();
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let block = builder
        .add_with_tangent("block", Mass::new(2.0).unwrap(), rest)
        .unwrap();
    let belt = builder
        .add_with_tangent(
            "belt",
            FlowSource::new(wave.clone()).publishing(),
            SourceAccount::default(),
        )
        .unwrap();
    let law = Slider::new(20.0, 0.5, 0.2).unwrap().told(SliderEnd::B);
    let friction = builder.add_on_line("friction", law, Slip::stuck()).unwrap();
    let lining = builder
        .add_with_tangent(
            "lining",
            HeatCapacity::new(50.0).unwrap().publishing(),
            300.0,
        )
        .unwrap();
    builder
        .common_flow(port(&builder, block, "p"), &[port(&builder, friction, "a")])
        .unwrap();
    builder
        .common_flow(port(&builder, belt, "p"), &[port(&builder, friction, "b")])
        .unwrap();
    builder
        .signal(
            port(&builder, belt, "span"),
            port(&builder, friction, "span_b"),
        )
        .unwrap();
    let warmth = port(&builder, friction, "heat");
    builder
        .common_effort(port(&builder, lining, "heat"), &[warmth])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    let at_block = system.port(friction, "a").unwrap();
    let mut steps = Vec::new();
    with_cx(|cx| {
        while 1.5 - system.time_s() > 1.5e-12 {
            assert!(steps.len() < 600, "the clock did not advance");
            let report = system
                .step(0.005_f64.min(1.5 - system.time_s()), cx)
                .unwrap();
            let slip = *system.state::<Slip>(friction).unwrap();
            let block_m_s = system.state::<Motion>(block).unwrap().velocity_m_s;
            steps.push(OnTheBelt {
                time_s: report.time_s,
                force_n: system.last_output(at_block)[0],
                record_m_s: slip.slip_m_s(),
                stuck: slip.is_stuck(),
                slip_m_s: wave.value_at(report.time_s) - block_m_s,
            });
        }
    });
    steps
}

/// Prints steps 248 to 256 of `steps`, about the catches.
fn shown(label: &str, steps: &[OnTheBelt]) {
    use std::fmt::Write as _;
    let mut rows = String::new();
    for (step, number) in steps[247..256].iter().zip(248..) {
        write!(
            rows,
            "\n  step {number}, to {:.6} s: {:e} N, record {:e}, stuck {}, slip {:e} m/s",
            step.time_s, step.force_n, step.record_m_s, step.stuck, step.slip_m_s
        )
        .unwrap();
    }
    eprintln!("the friction belt, {label}:{rows}");
}

/// What the catches of [`the_friction_belt`] share before them: 300 steps,
/// the block held to the belt in steps 1 to 50 (the first ramp needs 8 N of
/// the 10 N the slider holds), and sliding at 4 N from the break-away in
/// step 51 (the second needs 16 N) through step `sliding_to`.
fn carried_then_sliding(steps: &[OnTheBelt], sliding_to: usize) {
    assert_eq!(steps.len(), 300);
    for (index, step) in steps.iter().enumerate() {
        let number = index + 1;
        if number <= 50 {
            assert!(step.stuck, "step {number}: {step:?}");
        } else if number <= sliding_to {
            let sliding = !step.stuck && (step.force_n - 4.0).abs() <= 1.0e-12;
            assert!(sliding, "step {number}: {step:?}");
        }
    }
}

/// Step `number` of `steps` arrests the block: it ends stuck, with the
/// slider's record at exactly zero and a force inside the sliding limit,
/// and every step after it is stuck with the block on the belt.
fn arrested_at(steps: &[OnTheBelt], number: usize) {
    let arrest = steps[number - 1];
    assert!(
        arrest.stuck && arrest.record_m_s == 0.0,
        "step {number}: {arrest:?}"
    );
    assert!(
        arrest.force_n > 0.0 && arrest.force_n <= 4.0,
        "step {number}: {arrest:?}"
    );
    for (index, step) in steps.iter().enumerate().skip(number) {
        assert!(step.stuck, "step {}: {step:?}", index + 1);
        assert!(
            step.slip_m_s.abs() < 1.0e-12,
            "step {}: {step:?}",
            index + 1
        );
    }
}

#[test]
fn a_belt_caught_at_the_end_of_a_step_holds_the_block_from_two_steps_after() {
    // The catch falls at 1.25 s, the end of step 250, where the force that
    // ends the step at no slip is 400 * 0.01 m/s = 4 N, the sliding force,
    // in exact arithmetic: rounding decides whether that step slides or is
    // arrested, and neither is asserted. Sliding, it leaves a slip of
    // rounding size, and the slider's record of it and the block's own
    // differ by their rounding. Step 251 is then arrested; or released,
    // where the block's slip is below half the record, which leaves the
    // record on the other side of it; or turned, where it already was. A
    // step that starts with the record on the other side of the slip is
    // turned, the slip being far inside c mu_k N: it ends with no mean slip
    // under a force inside the limit, a hold, stuck from whichever state
    // it came. So the pair is stuck from step 252 (catch + 2), and held
    // after, to 1.5 s: breaking away takes 12.5 mm/s of released slip.
    // While turned was recorded as sliding, the record and the block's slip
    // kept opposite signs, every step from 252 was turned, and the pair was
    // never stuck again: fs-cli's belt ended with stuck 0.
    let steps = the_friction_belt(3.0);
    shown("caught at the end of step 250", &steps);
    assert!(
        (steps[249].time_s - 1.25).abs() < 1.0e-12,
        "{:?}",
        steps[249]
    );
    carried_then_sliding(&steps, 249);
    for (index, step) in steps.iter().enumerate().skip(251) {
        assert!(step.stuck, "step {}: {step:?}", index + 1);
        assert!(
            step.slip_m_s.abs() < 1.0e-12,
            "step {}: {step:?}",
            index + 1
        );
    }
}

#[test]
fn a_belt_caught_inside_a_step_arrests_the_block_there_at_exactly_no_slip() {
    // A last speed of 3.005 m/s moves the catch to 1.2525 s, half way
    // through step 251. Step 250 starts with 15 mm/s of slip, 6 N to stop
    // it, and slides at 4 N; step 251 starts with 5 mm/s, 2 N, and is
    // arrested there, its record at exactly zero, and every step after it
    // is held from that zero. No step is turned, so it passes whether
    // turned is recorded stuck or sliding.
    let steps = the_friction_belt(3.005);
    shown("caught inside step 251", &steps);
    carried_then_sliding(&steps, 250);
    arrested_at(&steps, 251);
}

#[test]
fn a_belt_caught_just_after_a_step_slides_through_it_and_arrests_in_the_next() {
    // A last speed of 3.000000002 m/s moves the catch 1e-9 s past the end
    // of step 250. That step starts with 0.010000002 m/s of slip, which
    // takes 4.0000008 N to stop: over the sliding limit by 8e-7 N, far
    // above any rounding, so it slides at 4 N. Step 251 starts with 2e-9
    // m/s and is arrested by 8e-7 N, its record at exactly zero. No step is
    // turned, so it passes whether turned is recorded stuck or sliding.
    let steps = the_friction_belt(3.000_000_002);
    shown("caught 1e-9 s after the end of step 250", &steps);
    carried_then_sliding(&steps, 250);
    arrested_at(&steps, 251);
}

#[test]
fn no_drawn_side_told_its_span_breaks_a_pair_away_that_holds() {
    // A standing property. 60 drawn sides, piecewise linear with one to
    // five changes of slope in 0.4 s, anywhere within their steps but at
    // least three steps apart, and slopes within 2 m/s2, under drawn
    // pushes on the mass, such that the friction the pair needs,
    // m V' - P, stays within 0.7 of the static limit. Told its span, the
    // slider never lets go; its force is above the largest needed by at
    // most m times a quarter of the largest change of acceleration, which
    // is what one change alone can ask (changes in consecutive steps add
    // up: the test above); three steps after the side's last change of
    // slope the mass is on the side.
    let (dt, hold_n) = (1.0e-3_f64, HOLDS * NORMAL_N);
    let mut draw = Draw(20_261_113);
    let (mut taken, mut worst_over, mut worst_heat_j) = (0_u32, 0.0_f64, 0.0_f64);
    while taken < 60 {
        let mut times: Vec<f64> = (0..=draw.next() % 5).map(|_| draw.between(0.0, 0.4)).collect();
        times.sort_by(f64::total_cmp);
        if times.windows(2).any(|pair| pair[1] - pair[0] < 3.0 * dt) {
            continue;
        }
        let slopes: Vec<f64> = (0..=times.len()).map(|_| draw.between(-2.0, 2.0)).collect();
        // The side's speed: from rest with the first slope, each change at
        // its time, the last slope held for a second beyond the run.
        let mut knots = vec![(0.0, 0.0)];
        for (index, &slope) in slopes.iter().enumerate() {
            let (from, speed) = knots[index];
            let to = times.get(index).copied().unwrap_or(1.4);
            knots.push((to, speed + slope * (to - from)));
        }
        if knots.windows(2).any(|pair| pair[1].0 <= pair[0].0) {
            continue;
        }
        let side = Waveform::table(knots).unwrap();
        let mut pushed = vec![(0.0, draw.between(-2.0, 2.0))];
        let mut at = 0.0;
        for _ in 0..=draw.next() % 4 {
            at += draw.between(0.01, 0.1);
            pushed.push((at, draw.between(-2.0, 2.0)));
        }
        let push = Waveform::table(pushed.clone()).unwrap();
        // The friction needed is piecewise linear in time between the
        // changes of either, so its largest size is at one of them, on
        // one of its two sides.
        let mut instants = vec![0.0];
        instants.extend(times.iter().copied());
        instants.extend(pushed.iter().map(|knot| knot.0).filter(|&time| time < 0.4));
        instants.push(0.4);
        let slope_at = |time: f64| slopes[times.iter().filter(|&&change| change <= time).count()];
        let needed_n = instants
            .iter()
            .flat_map(|&time| [time, (time - 1.0e-9).max(0.0)])
            .map(|time| (MASS_KG * slope_at(time) - push.value_at(time)).abs())
            .fold(0.0_f64, f64::max);
        if needed_n > 0.7 * hold_n {
            continue;
        }
        let change = slopes.windows(2).map(|pair| (pair[1] - pair[0]).abs()).fold(0.0_f64, f64::max);
        let last_change = times.last().copied().unwrap_or(0.0);
        taken += 1;
        let (mut system, mass, rubbing, body) = carried(Moved::Told(side.clone()), Some(push));
        let at_mass = system.port(rubbing, "a").unwrap();
        with_cx(|cx| {
            let mut step = 0;
            while system.time_s() < 0.4 {
                assert!(step < 600, "side {taken}: clock did not advance");
                system.step(dt, cx).unwrap();
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "side {taken}, step {step}");
                let over = system.last_output(at_mass)[0].abs() - needed_n;
                assert!(over <= 0.25 * MASS_KG * change + 1.0e-9, "side {taken}, step {step}: {over}");
                worst_over = worst_over.max(over / (MASS_KG * change));
                record_is_the_slip(&system, mass, rubbing, side.value_at(system.time_s()), "a drawn side");
                if system.time_s() >= last_change + 3.0 * dt {
                    let behind = side.value_at(system.time_s()) - system.state::<Motion>(mass).unwrap().velocity_m_s;
                    assert!(behind.abs() <= 1.0e-12, "side {taken}, step {step}: {behind}");
                }
                step += 1;
            }
        });
        let heat_j = received_j(&system, body);
        assert!((0.0..=1.0e-4).contains(&heat_j), "side {taken}: {heat_j}");
        worst_heat_j = worst_heat_j.max(heat_j);
    }
    eprintln!(
        "60 drawn sides told their span: none breaks away; the force is above the largest needed \
         by at most {worst_over:.3} of m times the largest change of acceleration; {worst_heat_j:e} J \
         of heat at most"
    );
}

#[test]
fn a_span_on_the_wrong_end_or_of_the_wrong_sign_is_refused_by_name() {
    // A span on an end that answers the slider's force. The mass is at
    // `a`; a slider told the span of `a` is refused at its first step,
    // by name, and nothing is published.
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let rubbing = builder
        .add_on_line("slider", slider().told(SliderEnd::A), Slip::stuck())
        .unwrap();
    builder
        .common_flow(port(&builder, mass, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    for (name, value) in [("b", 0.0), ("heat", 300.0)] {
        builder
            .boundary(port(&builder, rubbing, name), Box::new(move |_, _, values| values[0] = value))
            .unwrap();
    }
    builder
        .boundary(port(&builder, rubbing, "span_a"), Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    let mut system = builder.seal().unwrap();
    let refusal = with_cx(|cx| system.step(1.0e-3, cx)).unwrap_err().to_string();
    assert!(
        refusal.contains("`slider`") && refusal.contains("Span") && refusal.contains("answers the slider's force"),
        "{refusal}"
    );
    assert_eq!(system.steps(), 0);

    // A span of the port's sign. The side moves at 0.2 m/s and faster, so
    // the port receives a flow below zero and the span is above it. One
    // step cannot say that a span on the flow's side of zero is wrong: a
    // side that reverses twice within a step has both its ends on one
    // side and its mean on the other. And the slider may not decide by a
    // step's inputs at all: the kernel also proposes at inputs it has
    // perturbed. So the slider counts the steps taken in a row with the
    // span on the flow's side, and refuses the step after two of them, by
    // name. Given as the mover returns it, every step is taken.
    let wave = Waveform::table(vec![(0.0, 0.2), (100.0, 400.2)]).unwrap();
    for mirrored in [true, false] {
        let (means, values) = (wave.clone(), wave.clone());
        let sign = if mirrored { -1.0 } else { 1.0 };
        let side = Moved::ToldByABoundary {
            mean: Box::new(move |t0, dt| means.mean(t0, dt)),
            value: Box::new(move |t| sign * values.value_at(t)),
        };
        let (mut system, _, _, _) = carried(side, None);
        let outcomes: Vec<_> = with_cx(|cx| (0..4).map(|_| system.step(1.0e-3, cx)).collect());
        if mirrored {
            assert!(outcomes[0].is_ok() && outcomes[1].is_ok(), "{outcomes:?}");
            for refused in &outcomes[2..] {
                let refusal = refused.as_ref().unwrap_err().to_string();
                assert!(
                    refusal.contains("`slider`")
                        && refusal.contains("Span")
                        && refusal.contains("in its last two steps"),
                    "{refusal}"
                );
            }
            assert_eq!(system.steps(), 2);
        } else {
            assert!(outcomes.iter().all(Result::is_ok), "{outcomes:?}");
        }
    }
}

/// A mass riding stuck on a side moved as `wave` and told its span, the
/// mass starting at `mass_m_s`; a publishing `FlowSource` schedules the
/// knots, or raw boundaries leave them unresolved when `scheduled` is false.
/// With `beside`, a second mass of 2 kg at that velocity on a spring of 100 N/m
/// to the first. Returns the system and the slider.
fn riding_a_told_side(
    wave: Waveform,
    mass_m_s: f64,
    beside: Option<f64>,
    scheduled: bool,
) -> (System, NodeId) {
    let at = |velocity_m_s: f64| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let mut builder = SystemBuilder::new();
    let mass = builder.add_with_tangent("mass", Mass::new(MASS_KG).unwrap(), at(mass_m_s)).unwrap();
    let rubbing = builder
        .add_on_line("slider", slider().told(SliderEnd::B), Slip::stuck())
        .unwrap();
    let mut members = vec![port(&builder, rubbing, "a")];
    if let Some(other_m_s) = beside {
        let other = builder.add_with_tangent("other", Mass::new(2.0).unwrap(), at(other_m_s)).unwrap();
        let link = builder
            .add_with_tangent("link", Spring::new(100.0, 0.0).unwrap(), 0.0)
            .unwrap();
        members.push(port(&builder, link, "a"));
        builder
            .common_flow(port(&builder, other, "p"), &[port(&builder, link, "b")])
            .unwrap();
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    if scheduled {
        let source = builder
            .add_with_tangent(
                "side",
                FlowSource::new(wave).publishing(),
                SourceAccount::default(),
            )
            .unwrap();
        builder
            .bond(port(&builder, rubbing, "b"), port(&builder, source, "p"))
            .unwrap();
        builder
            .signal(
                port(&builder, source, "span"),
                port(&builder, rubbing, "span_b"),
            )
            .unwrap();
    } else {
        let means = wave.clone();
        builder
            .boundary(
                port(&builder, rubbing, "b"),
                Box::new(move |t, dt, values| values[0] = -means.mean(t, dt)),
            )
            .unwrap();
        builder
            .boundary(
                port(&builder, rubbing, "span_b"),
                Box::new(move |t, dt, values| {
                    values[0] = wave.value_at(t);
                    values[1] = wave.value_at(t + dt);
                }),
            )
            .unwrap();
    }
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), rubbing)
}

#[test]
fn a_slow_told_side_and_one_that_reverses_twice_in_a_step_are_taken() {
    // The first version of the sign check looked at one step's inputs and
    // refused two sides that were in order (found by the review of
    // b47fe4f2a8).
    //
    // A slow side beside a fast body. A told end that a `FlowSource`
    // moves is an unknown of its block, and the kernel differences the
    // slider in it by sqrt(eps) of the largest flow of the block. With a
    // second mass at 0.2 m/s that is 3e-9 m/s: the flow of a side at
    // 1e-9 m/s changes sign in the perturbed proposal, and the check
    // refused the first step, and every retry of it. Likewise a side that
    // leaves rest at 1e-6 m/s2, at its second step.
    let dt = 1.0e-3;
    let slow: [(Waveform, f64); 2] = [
        (Waveform::constant(1.0e-9).unwrap(), 1.0e-9),
        (Waveform::table(vec![(0.0, 0.0), (1000.0, 1.0e-3)]).unwrap(), 0.0),
    ];
    for (wave, mass_m_s) in slow {
        let (mut system, rubbing) = riding_a_told_side(wave, mass_m_s, Some(0.2), true);
        with_cx(|cx| {
            for step in 0..300 {
                system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "step {step}");
            }
        });
    }

    // A raw boundary leaves the next sides' interior knots unresolved.
    // A side that reverses twice in one step: it starts at 0.1 mm/s,
    // slows at 2 m/s2 for 0.4 ms and speeds up at 2 m/s2 from then, so the
    // first step's ends are 0.1 and 0.5 mm/s and its mean is -0.18 mm/s.
    // The step is taken, and the steps after it.
    let once = Waveform::table(vec![(0.0, 1.0e-4), (0.4e-3, -7.0e-4), (1.0e-3, 5.0e-4), (1.0, 5.0e-4)]).unwrap();
    assert!(once.mean(0.0, dt) < -1.7e-4 && once.value_at(0.0) > 0.0 && once.value_at(dt) > 0.0);
    let (mut system, _) = riding_a_told_side(once, 1.0e-4, None, false);
    with_cx(|cx| {
        for step in 0..20 {
            system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
        }
    });

    // And again five steps later: the count is of steps in a row, and a
    // step in order in between starts it anew.
    let apart = Waveform::table(vec![
        (0.0, 1.0e-4),
        (0.4e-3, -7.0e-4),
        (1.0e-3, 5.0e-4),
        (5.0e-3, 5.0e-4),
        (5.4e-3, -7.0e-4),
        (6.0e-3, 5.0e-4),
        (1.0, 5.0e-4),
    ])
    .unwrap();
    assert!(apart.mean(5.0 * dt, dt) < 0.0 && apart.value_at(6.0 * dt) > 0.0);
    let (mut system, _) = riding_a_told_side(apart, 1.0e-4, None, false);
    with_cx(|cx| {
        for step in 0..20 {
            system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
        }
    });

    // The same twice in a row is a motion the step does not resolve, and
    // it is what a span of the port's sign looks like: the third step is
    // refused by name, as that span is.
    let twice = Waveform::table(vec![
        (0.0, 1.0e-4),
        (0.4e-3, -7.0e-4),
        (1.0e-3, 5.0e-4),
        (1.4e-3, -7.0e-4),
        (2.0e-3, 5.0e-4),
        (1.0, 5.0e-4),
    ])
    .unwrap();
    assert!(twice.mean(dt, dt) < 0.0 && twice.value_at(2.0 * dt) > 0.0);
    let (mut system, _) = riding_a_told_side(twice.clone(), 1.0e-4, None, false);
    let outcomes: Vec<_> = with_cx(|cx| (0..3).map(|_| system.step(dt, cx)).collect());
    assert!(outcomes[0].is_ok() && outcomes[1].is_ok(), "{outcomes:?}");
    let refusal = outcomes[2].as_ref().unwrap_err().to_string();
    assert!(refusal.contains("Span") && refusal.contains("does not resolve the motion"), "{refusal}");
    assert_eq!(system.steps(), 2);

    // The source advertises those same knots, so both reversals are
    // resolved and their admissible acceleration keeps the pair held.
    let (mut system, rubbing) = riding_a_told_side(twice, 1.0e-4, None, true);
    with_cx(|cx| {
        let mut steps = 0;
        while system.time_s() < 0.02 {
            assert!(steps < 40, "scheduled side did not advance");
            system.step(dt, cx).unwrap();
            assert!(system.state::<Slip>(rubbing).unwrap().is_stuck());
            steps += 1;
        }
    });
}

/// The speed of a side that leaves rest with `slopes[0]` and changes to the
/// next slope at each of `changes_s`, the last one held for a second
/// beyond.
fn changing_slope(slopes: &[f64], changes_s: &[f64]) -> Waveform {
    let mut knots = vec![(0.0, 0.0)];
    for (index, &slope) in slopes.iter().enumerate() {
        let (from, speed) = knots[index];
        let to = changes_s.get(index).copied().unwrap_or(from + 1.0);
        knots.push((to, speed + slope * (to - from)));
    }
    Waveform::table(knots).unwrap()
}

#[test]
fn changes_of_slope_in_consecutive_steps_add_up_and_can_break_the_pair_away() {
    // The limit of a raw boundary's told span, further. A step held on the
    // mean leaves a slip in the record; a change of slope in the very next step adds
    // its own before that one is made up, and so on while the changes come
    // in consecutive steps. (b47fe4f2a8 said the force is above the needed
    // one by a quarter of the largest change of acceleration at most, and
    // asserted it over 60 drawn sides, whose changes fell in consecutive
    // steps too rarely to show otherwise. Found by its review.)
    let dt = 1.0e-3;
    let forces = |slopes: &[f64], changes_s: &[f64], push_n: f64, steps: usize| {
        let side = changing_slope(slopes, changes_s);
        let (mut system, _, rubbing, _) = carried(
            told_by_a_boundary(&side),
            Some(Waveform::constant(push_n).unwrap()),
        );
        let at_mass = system.port(rubbing, "a").unwrap();
        with_cx(|cx| {
            (0..steps)
                .map(|_| {
                    system.step(dt, cx).unwrap();
                    (
                        system.last_output(at_mass)[0],
                        system.state::<Slip>(rubbing).unwrap().is_stuck(),
                    )
                })
                .collect::<Vec<_>>()
        })
    };
    // Changes in three consecutive steps, of 1.48, 1.88 and 3.02 m/s2,
    // under a push of -1.92 N: the slider needs 4 N at most, and gives
    // 6.917 N in the fourth step. That is 2.917 N above, 0.483 of m times
    // the largest change, where one change alone gives a quarter.
    let seen = forces(&[-1.58, -0.10, -1.98, 1.04], &[0.52e-3, 1.27e-3, 2.39e-3], -1.92, 8);
    let wanted = [-0.558, 1.194, -2.752, 6.917, 4.0, 4.0, 4.0, 4.0];
    for (step, ((force, stuck), wanted)) in seen.iter().zip(wanted).enumerate() {
        assert!((force - wanted).abs() <= 1.0e-3 && *stuck, "step {step}: {force}");
    }
    // Inside what the property below drew from before it kept its changes
    // apart (slopes within 2 m/s2, a push of 2 N, 6 N needed at most of
    // the 9.8 N the slider holds): changes in four consecutive steps, and
    // the pair breaks away in the sixth.
    let seen = forces(
        &[-1.3186, 0.8725, -1.7886, 2.0, -2.0],
        &[1.6012e-3, 2.4005e-3, 3.0914e-3, 4.463e-3],
        2.0,
        8,
    );
    let wanted = [-4.637, -2.89, -2.168, -1.878, 3.507];
    for (step, ((force, stuck), wanted)) in seen.iter().zip(wanted).enumerate() {
        assert!((force - wanted).abs() <= 1.0e-3 && *stuck, "step {step}: {force}");
    }
    for (step, (force, stuck)) in seen.iter().enumerate().skip(5) {
        assert!(*force == -SLIDES * NORMAL_N && !*stuck, "step {step}: {force}");
    }

    // Measured: 2000 drawn sides with a change of slope in each of two to
    // five consecutive steps (slopes within 2 m/s2, a constant push
    // within 2 N, at most 0.7 of the static limit needed). None of these
    // draws breaks away; the pair above was found by a search.
    let mut draw = Draw(20_261_116);
    let (mut drawn, mut broke, mut worst_of_largest, mut worst_of_sum) = (0_u32, 0_u32, 0.0_f64, 0.0_f64);
    while drawn < 2000 {
        let changes = 2 + (draw.next() % 4) as usize;
        let slopes: Vec<f64> = (0..=changes).map(|_| draw.between(-2.0, 2.0)).collect();
        let times: Vec<f64> = (0..changes)
            .map(|step| (1.0 + step as f64 + draw.between(0.02, 0.98)) * dt)
            .collect();
        let push_n = draw.between(-2.0, 2.0);
        let needed_n = slopes.iter().map(|slope| (MASS_KG * slope - push_n).abs()).fold(0.0_f64, f64::max);
        if needed_n > 0.7 * HOLDS * NORMAL_N {
            continue;
        }
        drawn += 1;
        let steps: Vec<f64> = slopes.windows(2).map(|pair| (pair[1] - pair[0]).abs()).collect();
        let (largest, sum) = (steps.iter().copied().fold(0.0_f64, f64::max), steps.iter().sum::<f64>());
        let seen = forces(&slopes, &times, push_n, changes + 6);
        if seen.iter().any(|(_, stuck)| !stuck) {
            broke += 1;
            continue;
        }
        let over = seen.iter().map(|(force, _)| force.abs()).fold(0.0_f64, f64::max) - needed_n;
        worst_of_largest = worst_of_largest.max(over / (MASS_KG * largest));
        worst_of_sum = worst_of_sum.max(over / (MASS_KG * sum));
    }
    eprintln!(
        "changes of slope in 2 to 5 consecutive steps, 2000 drawn sides: {broke} break away; of \
         the others the force is above the largest needed by at most {worst_of_largest:.3} of m \
         times the largest change of acceleration, {worst_of_sum:.3} of m times their sum"
    );
    // More than the quarter of one change, and within a quarter of their
    // sum: each change adds at most its own quarter. (Asserted for these
    // draws; not proved.)
    assert!(worst_of_largest > 0.3 && worst_of_sum <= 0.25, "{worst_of_largest} {worst_of_sum}");
}

#[test]
fn a_slider_between_two_sides_told_their_spans_ends_with_the_slip_they_end_with() {
    // Both ends prescribed and told their spans: no line, the slider
    // slides at its sliding force, its loss is that force times the mean
    // slip, and its record after each step is the slip the two sides end
    // the step with, not the mean slip a held pair would leave.
    //
    // The slip it is declared to start with is not used. With both ends
    // given, each step's force follows from its mean slip and its record
    // from its mean and its two ends, whatever the step before left:
    // declared sliding with the slip the sides start with, sliding the
    // other way, or stuck, the forces and the records are the same to the
    // bit.
    let (dt, steps) = (1.0e-3_f64, 200_u32);
    // The two sides' speeds as their movers return them: 0.3 + 2 t and
    // 0.1 - t.
    let speeds: [fn(f64) -> f64; 2] = [|t| 0.3 + 2.0 * t, |t| 0.1 - t];
    let run = |declared: Slip| {
        let mut builder = SystemBuilder::new();
        let rubbing = builder
            .add_on_line("slider", slider().told(SliderEnd::A).told(SliderEnd::B), declared)
            .unwrap();
        let body = builder
            .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
            .unwrap();
        for (end, span, speed) in [("a", "span_a", speeds[0]), ("b", "span_b", speeds[1])] {
            builder
                .boundary(
                    port(&builder, rubbing, end),
                    Box::new(move |t0, dt, values| values[0] = -0.5 * (speed(t0) + speed(t0 + dt))),
                )
                .unwrap();
            builder
                .boundary(
                    port(&builder, rubbing, span),
                    Box::new(move |t0, dt, values| {
                        values[0] = speed(t0);
                        values[1] = speed(t0 + dt);
                    }),
                )
                .unwrap();
        }
        let warmth = port(&builder, rubbing, "heat");
        builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(warmth, &[300.0]));
        let at_a = system.port(rubbing, "a").unwrap();
        let (mut loss_j, mut seen) = (0.0, Vec::new());
        with_cx(|cx| {
            for step in 0..steps {
                let t0 = system.time_s();
                system.step(dt, cx).unwrap();
                let ends = speeds[0](t0 + dt) - speeds[1](t0 + dt);
                let record = system.state::<Slip>(rubbing).unwrap().slip_m_s();
                assert!((record + ends).abs() <= 1.0e-15, "step {step}: {record} against {}", -ends);
                let mean = 0.5 * (speeds[0](t0) - speeds[1](t0) + ends);
                loss_j += dt * SLIDES * NORMAL_N * mean;
                seen.push((system.last_output(at_a)[0].to_bits(), record.to_bits()));
            }
        });
        let warmed_j = received_j(&system, body);
        assert!((warmed_j - loss_j).abs() <= 1.0e-12 * loss_j, "{warmed_j} vs {loss_j}");
        seen
    };
    let starts_with = -(speeds[0](0.0) - speeds[1](0.0));
    let seen = run(Slip::sliding(starts_with));
    assert_eq!(run(Slip::sliding(-starts_with)), seen);
    assert_eq!(run(Slip::stuck()), seen);
}

#[test]
fn a_held_mass_is_taken_under_any_push_and_a_planted_flow_is_still_refused() {
    // The flow the kernel publishes for a line input is the released flow
    // less what the efforts take away. For a held slider that is
    // r - c (r / c): zero for some numbers and a unit in the last place of
    // r for others. The junction's audit has to be told that rounding;
    // told nothing, it refused 282 of these 2000 pushes at their first
    // step, on an earlier kernel. Here none is refused, and the mass is at
    // rest to that rounding.
    let dt = 1.0e-3_f64;
    let c = 0.5 * dt / MASS_KG;
    let (mut uneven, mut worst) = (0_u32, 0.0_f64);
    for k in 0..2000_u32 {
        let push = HOLDS * NORMAL_N * (0.05 + 0.000_437 * f64::from(k));
        let (mut held, mass, rubbing, _) = on_the_ground(0.0, Some(push));
        with_cx(|cx| {
            for step in 0..6 {
                held.step(dt, cx)
                    .unwrap_or_else(|error| panic!("push {push} step {step}: {error}"));
            }
        });
        let flow = held.last_input(held.port(rubbing, "a").unwrap())[0].abs();
        let velocity = held.state::<Motion>(mass).unwrap().velocity_m_s.abs();
        // The rounding of the released flow, c * push.
        let rounding = 16.0 * f64::EPSILON * c * push;
        assert!(flow <= rounding && velocity <= rounding, "push {push}: {flow} {velocity}");
        assert!(held.state::<Slip>(rubbing).unwrap().is_stuck());
        uneven += u32::from(flow != 0.0 || velocity != 0.0);
        worst = worst.max(flow).max(velocity);
    }
    eprintln!(
        "held under 2000 pushes: none refused; {uneven} with a flow or a velocity that is not \
         exactly zero, {worst:e} m/s at worst"
    );
    // Since every input of the line's component is released before its
    // first ask, each of these flows and velocities comes out at exactly
    // zero (none of the 2000 is otherwise on ts2 chain fsz, 2026-11-14;
    // 282 were on the kernel before it, chain fsx), so the sweep no longer
    // reaches a push whose flow carries the rounding. That the audit is
    // told the rounding, and no more, rests on the planted half below: the
    // refusal of a hundred times what is allowed, and the two terms read
    // back from the step that takes a quarter of it.
    assert_eq!(uneven, 0, "{worst:e} m/s at worst");

    // The planted negative: what the audit is told is the rounding of that
    // flow and no more. A mass held with half the static limit is brought
    // to rest; then a net force F is planted on the mass's input, where
    // the solve starts. The plant has to get past the solve before the
    // verdict is the audit's, so the solve is told to take a wire that is
    // open by four times its own size: the planted motion, c F on wires
    // whose size it now is, then reaches the audit at the first check.
    // (The solve judges a flow on the line by what crosses, so a tolerance
    // below one no longer takes a planted motion; the end of this test
    // shows what it does with it.) The mass then moves by c F with a force
    // on it that nothing returns: the junction is left with F times that
    // motion, c dt F^2, and the verdict is the audit's.
    let push = 0.5 * HOLDS * NORMAL_N;
    let (mut system, mass, rubbing, _) = on_the_ground(0.0, Some(push));
    with_cx(|cx| {
        for _ in 0..3 {
            system.step(dt, cx).unwrap();
        }
    });
    let (net, at_slider) = (system.port(mass, "p").unwrap(), system.port(rubbing, "a").unwrap());
    let at_rest = system.checkpoint().unwrap();
    let settled_n = system.last_input(net)[0];
    let held_n = system.last_output(at_slider)[0].abs();
    assert!((held_n - push).abs() <= 1.0e-12 * push && settled_n.abs() <= 1.0e-12 * push);
    let loose = SolveConfig {
        relative_tolerance: 4.0,
        ..SolveConfig::default()
    };
    let planted = |system: &mut System, extra_n: f64| {
        system.restore(&at_rest).unwrap();
        system.set_config(loose);
        assert!(system.guess(net, &[settled_n + extra_n]));
        with_cx(|cx| system.step(dt, cx))
    };
    // The audit's tolerance for this junction, read from the refusal of a
    // planted 10 micronewtons, as a flow at the slider's port.
    let tolerance_j = match planted(&mut system, 1.0e-5) {
        Err(SystemError::Interface {
            connection: 0,
            tolerance_j,
            ..
        }) => tolerance_j,
        other => panic!("a planted force is refused by the first connection's audit: {other:?}"),
    };
    let allowed_m_s = tolerance_j / (dt * held_n);
    // Two terms by name, each 16 units of eps. The published flow's own
    // rounding: the released flow and what the effort takes away, c * push
    // each. And the wire's rounding: what its target depends on, the
    // released flow. As a flow at the slider's port, three times
    // 16 eps c push in all. (What the junction's other two wires are told
    // is the rounding of the largest of their kind, times a mass that
    // does not move and times a flow of 2.5e-9 m/s: nothing beside it.)
    let published = 16.0 * f64::EPSILON * (2.0 * c * push);
    let wire = 16.0 * f64::EPSILON * (c * push);
    let within = |got: f64, named: f64| (got - named).abs() <= 0.05 * named;
    // The audit's tolerance is their sum.
    let sum = published + wire;
    assert!(within(allowed_m_s, sum), "{allowed_m_s} against {published} and {wire}");
    eprintln!(
        "the audit allows a held slider's port {allowed_m_s:e} m/s: {:.3} of the published \
         flow's own rounding term",
        allowed_m_s / published
    );
    // A hundred times what is allowed is refused by that junction's audit,
    // by that factor, and no step is counted.
    let leaving = |times: f64| (times * tolerance_j / (c * dt)).sqrt();
    let refusal = planted(&mut system, leaving(100.0)).expect_err("refused");
    let SystemError::Interface {
        connection,
        defect_j,
        tolerance_j,
    } = refusal
    else {
        panic!("expected the audit's refusal, got {refusal}");
    };
    assert_eq!(connection, 0);
    let ratio = defect_j.abs() / tolerance_j;
    assert!(ratio > 80.0 && ratio < 120.0, "{ratio}");
    assert_eq!(system.steps(), 3);
    // A quarter of it is taken, with no update. Of that step, the two
    // terms are what the solve told the audit for the wire, read back;
    // a restore forgets it.
    let report = planted(&mut system, leaving(0.25)).unwrap();
    assert_eq!(report.iterations, 0);
    assert_eq!(system.steps(), 4);
    let told = system.last_allowance(at_slider)[0];
    assert!(within(told, published + wire), "{told} against {published} and {wire}");
    system.restore(&at_rest).unwrap();
    assert_eq!(system.last_allowance(at_slider), [0.0]);

    // At the solve's ordinary tolerance the plant never reaches the audit:
    // the force that was refused a hundred times over is corrected in one
    // update, and the mass is back on its held state, to the rounding of
    // the released flow.
    system.restore(&at_rest).unwrap();
    system.set_config(SolveConfig::default());
    assert!(system.guess(net, &[settled_n + leaving(100.0)]));
    let report = with_cx(|cx| system.step(dt, cx)).unwrap();
    assert_eq!(report.iterations, 1);
    assert_eq!(system.steps(), 4);
    let rounding = 16.0 * f64::EPSILON * c * push;
    let flow = system.last_input(at_slider)[0].abs();
    let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s.abs();
    assert!(flow <= rounding && velocity <= rounding, "{flow} {velocity}");
    assert!((system.last_output(at_slider)[0].abs() - push).abs() <= 1.0e-12 * push);
    assert!(system.state::<Slip>(rubbing).unwrap().is_stuck());
}

/// A mass at rest stuck to a side whose velocity `side(t0, dt)` is
/// prescribed (in the sign the slider's port `b` receives it), the slider's
/// loss collected in a body at 300 K. Returns the system, the mass, the
/// slider and the body.
fn on_a_side(side: impl Fn(f64, f64) -> f64 + 'static) -> (System, NodeId, NodeId, NodeId) {
    on_a_side_at(false, false, side)
}

/// [`on_a_side`], or with the slider's ends exchanged: the mass at its `b`
/// and the side prescribed at its `a`. With `told`, the side's end is also
/// told its span: the side is held over each step, so that is the one
/// value twice, in the mover's sign, minus what the port receives.
fn on_a_side_at(
    exchanged: bool,
    told: bool,
    side: impl Fn(f64, f64) -> f64 + 'static,
) -> (System, NodeId, NodeId, NodeId) {
    let (at_mass, at_side) = if exchanged { ("b", "a") } else { ("a", "b") };
    let side = std::rc::Rc::new(side);
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let law = match (told, exchanged) {
        (false, _) => slider(),
        (true, false) => slider().told(SliderEnd::B),
        (true, true) => slider().told(SliderEnd::A),
    };
    let rubbing = builder.add_on_line("slider", law, Slip::stuck()).unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .common_flow(port(&builder, mass, "p"), &[port(&builder, rubbing, at_mass)])
        .unwrap();
    if told {
        let held = std::rc::Rc::clone(&side);
        builder
            .boundary(
                port(&builder, rubbing, if exchanged { "span_a" } else { "span_b" }),
                Box::new(move |t0, dt, values| values.fill(-held(t0, dt))),
            )
            .unwrap();
    }
    builder
        .boundary(
            port(&builder, rubbing, at_side),
            Box::new(move |t0, dt, values| values[0] = side(t0, dt)),
        )
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    (system, mass, rubbing, body)
}

#[test]
fn a_side_that_jumps_drags_the_mass_under_full_friction_until_it_is_caught() {
    dragged_by_a_jump(false, false);
    // The same with the slider's ends exchanged, the mass at `b` and the
    // side prescribed at `a`: the same motion and the same force on the
    // mass, the record of the slip with the other sign.
    dragged_by_a_jump(true, false);
    // And with the side's end told its span, at either end: a side that
    // is held and jumps says the one value twice, and the step of the
    // jump takes it up where the last one ended.
    dragged_by_a_jump(false, true);
    dragged_by_a_jump(true, true);
}

fn dragged_by_a_jump(exchanged: bool, told: bool) {
    // A mass at rest, stuck to a side at rest. From the third step the side
    // is prescribed at 0.5 m/s. A prescribed value is held over its step, so
    // that step starts with the whole jump as slip, although the slider
    // ended the step before stuck. 0.5 m/s cannot be caught within a step
    // by the 9.8 N the slider holds, so it breaks away and the mass is
    // dragged at the sliding force, 5.88 N, in every step: 2.94 mm/s a
    // step. After 170 such steps 0.2 mm/s are left, less than a step of
    // full friction takes, and the 171st ends at zero slip with the force
    // that needs, 0.4 N.
    //
    // (The first version of the slider kept one slip and advanced it by
    // the midpoint relation as if the side were a mass. After the jump its
    // record alternated between -0.994 and +0.006 m/s, friction acted in
    // every other step only, and every audit passed.)
    let (side_m_s, dt) = (0.5_f64, 1.0e-3_f64);
    let slide_n = SLIDES * NORMAL_N;
    let (mut system, mass, rubbing, body) =
        on_a_side_at(exchanged, told, move |t0, _| if t0 > 1.5e-3 { side_m_s } else { 0.0 });
    let at_mass = system.port(rubbing, if exchanged { "b" } else { "a" }).unwrap();
    // What a step of full friction takes off the slip.
    let per_step = dt * slide_n / MASS_KG;
    let full = (side_m_s / per_step).ceil() as u32 - 1;
    assert_eq!(full, 170);
    with_cx(|cx| {
        for step in 0..2 + full + 40 {
            system.step(dt, cx).unwrap();
            let force = system.last_output(at_mass)[0];
            let state = *system.state::<Slip>(rubbing).unwrap();
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            // The slider's record against the slip the two sides end the
            // step with: the port receives minus the mass's velocity.
            let side_now = if step >= 2 { side_m_s } else { 0.0 };
            let true_slip = if exchanged {
                side_now + velocity
            } else {
                -velocity - side_now
            };
            assert!(
                (state.slip_m_s() - true_slip).abs() <= 1.0e-12,
                "step {step}: {} against {true_slip}",
                state.slip_m_s()
            );
            if step < 2 {
                assert_eq!((force, velocity), (0.0, 0.0), "step {step}");
                assert!(state.is_stuck());
            } else if step < 2 + full {
                assert_eq!(force, -slide_n, "step {step}");
                assert!(!state.is_stuck(), "step {step}");
                let dragged = -per_step * f64::from(step - 1);
                assert!((velocity - dragged).abs() <= 1.0e-12, "step {step}: {velocity}");
            } else if step == 2 + full {
                let left = side_m_s - per_step * f64::from(full);
                let catching_n = -MASS_KG * left / dt;
                assert!((catching_n + 0.4).abs() <= 1.0e-9, "{catching_n}");
                assert!((force - catching_n).abs() <= 1.0e-9, "{force} against {catching_n}");
                assert!(state.is_stuck());
                assert_eq!(state.slip_m_s(), 0.0);
            } else {
                // Caught: it rides at the side's velocity with no force.
                assert!(state.is_stuck(), "step {step}");
                assert!(force.abs() <= 1.0e-8, "step {step}: {force}");
                assert!((velocity + side_m_s).abs() <= 1.0e-12, "step {step}: {velocity}");
            }
        }
    });
    // The side's work is the mass's kinetic energy and the heat, and the
    // two are equal: the side gave the impulse m V at the velocity V, the
    // mass holds m V^2 / 2, the rest was lost in the slip.
    let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
    let kinetic_j = 0.5 * MASS_KG * velocity * velocity;
    let heat_j = received_j(&system, body);
    let work_j = system
        .ledger()
        .boundaries
        .iter()
        .find(|account| account.port == if exchanged { "slider.a" } else { "slider.b" })
        .expect("the side's account")
        .energy_in_j
        .value();
    let given_j = MASS_KG * side_m_s * side_m_s;
    assert!((work_j - given_j).abs() <= 1.0e-12 * given_j, "{work_j}");
    assert!((kinetic_j + heat_j - work_j).abs() <= 1.0e-12 * given_j, "{kinetic_j} {heat_j}");
    assert!((heat_j - 0.5 * given_j).abs() <= 1.0e-12 * given_j, "{heat_j}");
}

#[test]
fn a_side_that_speeds_up_carries_the_mass_up_to_the_static_limit_and_no_further() {
    // The side's speed is (rate / 2) u^2 on its own clock u = t - dt: one
    // step at rest, then a growing acceleration rate * u. What a step is
    // given is the mean of that speed over the step, held, so the side
    // moves in a staircase and every step starts with the jump from the
    // last one as slip. Stuck, the slider gives the force that ends the
    // step at zero slip, m (V_n - V_(n-1)) / dt, which for these means is
    // m * rate * u at the start of the step, up to the static limit. The
    // mass is carried until that force exceeds 9.8 N: in the first step
    // that starts after u* = mu_s N / (m rate).
    let (rate, dt) = (48.0_f64, 1.0e-3_f64);
    let (hold_n, slide_n) = (HOLDS * NORMAL_N, SLIDES * NORMAL_N);
    let held = move |t0: f64| {
        let u = t0 - dt;
        if u < -0.5 * dt {
            0.0
        } else {
            0.5 * rate * (u * u + u * dt + dt * dt / 3.0)
        }
    };
    // The slider's `a` receives minus the mass's velocity, so the side is
    // given with the same sign.
    let (mut system, mass, rubbing, _) = on_a_side(move |t0, _| -held(t0));
    let at_mass = system.port(rubbing, "a").unwrap();
    let mut broke_at = None;
    let mut largest_carried_n = 0.0_f64;
    with_cx(|cx| {
        for step in 0..120_u32 {
            system.step(dt, cx).unwrap();
            let force = system.last_output(at_mass)[0];
            let state = *system.state::<Slip>(rubbing).unwrap();
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            let (now, before) = (
                held(f64::from(step) * dt),
                if step == 0 { 0.0 } else { held(f64::from(step - 1) * dt) },
            );
            let needed_n = MASS_KG * (now - before) / dt;
            if broke_at.is_none() && needed_n <= hold_n {
                // Carried: the force it needs, and the side's velocity.
                assert!(state.is_stuck(), "step {step}");
                assert!((force.abs() - needed_n).abs() <= 1.0e-9, "step {step}: {force}");
                assert!((velocity - now).abs() <= 1.0e-12, "step {step}: {velocity} {now}");
                largest_carried_n = largest_carried_n.max(needed_n);
            } else if broke_at.is_none() {
                assert_eq!(force.abs(), slide_n, "step {step}");
                assert!(!state.is_stuck());
                broke_at = Some(step);
            } else {
                // Left behind: it slides on under the sliding force.
                assert_eq!(force.abs(), slide_n, "step {step}");
                assert!(!state.is_stuck() && velocity < now, "step {step}");
            }
        }
    });
    let broke_s = f64::from(broke_at.expect("it breaks away")) * dt;
    let closed_s = dt + hold_n / (MASS_KG * rate);
    assert!(broke_s > closed_s && broke_s <= closed_s + dt, "{broke_s} against {closed_s}");
    assert_eq!(broke_at, Some(104));
    // It was carried with more than the sliding force could give.
    assert!(largest_carried_n > 1.6 * slide_n && largest_carried_n <= hold_n, "{largest_carried_n}");
}

#[test]
fn at_a_threshold_a_force_stays_inside_its_limit() {
    // A force that is a quotient, released / c, was compared through
    // `released <= c * limit` and could exceed the limit by a unit in the
    // last place: 3.0000000000000004 against 3.0, and there the force fell
    // as the released slip grew.
    let unit = Slider::new(3.0, 1.0, 1.0).unwrap();
    assert_eq!(
        unit.on_line(0.0, 0.300_000_000_000_000_04, 0.1, true),
        (3.0, 0.0, SliderRow::Held)
    );
    assert_eq!(
        unit.on_line(1.0, -0.300_000_000_000_000_04, 0.1, false),
        (-3.0, 0.0, SliderRow::Turned)
    );
    let (force, mean, row) = unit.on_line(1.0, -0.300_000_000_000_000_1, 0.1, false);
    assert_eq!((force, row), (-3.0, SliderRow::Reversed));
    assert!(force * mean >= 0.0);
    // No compliance: the slip is what it is, the force is the sliding one,
    // and none at exactly zero slip. A hair below zero is no compliance.
    assert_eq!(unit.on_line(1.0, -1.0, 0.0, false), (-3.0, -1.0, SliderRow::Slid));
    assert_eq!(unit.on_line(0.0, 2.0, -1.0e-20, true), (3.0, 2.0, SliderRow::Slid));
    assert_eq!(unit.on_line(0.0, 0.0, 0.0, true), (0.0, 0.0, SliderRow::Held));
    // Released exactly zero with a compliance: held with no force.
    assert_eq!(unit.on_line(0.0, 0.0, 0.1, true), (0.0, 0.0, SliderRow::Held));

    // With unequal coefficients the two limits are told apart. Stuck, the
    // slider holds up to the static limit, exactly; a hair beyond it gives
    // the sliding force.
    let pair = Slider::new(3.0, 1.0, 0.5).unwrap();
    let (c, hold, slide) = (0.1_f64, 3.0_f64, 1.5_f64);
    let hair = 1.0 + 4.0 * f64::EPSILON;
    assert_eq!(pair.on_line(0.0, c * hold, c, true), (hold, 0.0, SliderRow::Held));
    let beyond = c * hold * hair;
    assert_eq!(
        pair.on_line(0.0, beyond, c, true),
        (slide, beyond - c * slide, SliderRow::BrokeAway)
    );
    // Sliding, a released slip against the motion stops it within the
    // step up to the sliding limit, exactly, and reverses it beyond.
    assert_eq!(pair.on_line(1.0, -c * slide, c, false), (-slide, 0.0, SliderRow::Turned));
    let (force, mean, row) = pair.on_line(1.0, -c * slide * hair, c, false);
    assert_eq!((force, row), (-slide, SliderRow::Reversed));
    assert!(mean < 0.0 && mean > -1.0e-15);
    // A start slip that 2 N end within the step: a stuck pair carries it,
    // between the two limits; a sliding one cannot, and slides on.
    let start = 2.0e-3;
    let released = 0.5 * start + c * 2.0;
    let (force, mean, row) = pair.on_line(start, released, c, true);
    assert!((force - 2.0).abs() <= 1.0e-12 && force > slide && force < hold);
    assert_eq!((mean, row), (0.5 * start, SliderRow::Arrested));
    assert_eq!(
        pair.on_line(start, released, c, false),
        (slide, released - c * slide, SliderRow::Slid)
    );
}

#[test]
fn stuck_to_a_side_that_moved_the_slider_carries_holds_or_breaks_away() {
    // The three rows of a slider whose last step ended stuck and whose
    // step starts with a slip all the same, because the side it is stuck
    // to was given another value. Over start slips, compliances and a scan
    // of released slips: the loss is not negative and the force is inside
    // the static limit at every point; carried, the step ends at zero slip
    // (mean = start / 2) with a force along the start slip; held, the mean
    // slip is zero; broken away, the sliding force acts along the released
    // slip.
    let law = slider();
    let (hold, slide) = (HOLDS * NORMAL_N, SLIDES * NORMAL_N);
    let mut rows: HashMap<SliderRow, u32> = HashMap::new();
    let mut carried_most_n = 0.0_f64;
    for &start in &[2.0e-3, -2.0e-3, 1.0e-5, -0.4] {
        for &c in &[2.5e-4, 1.0e-3, 4.0e-2] {
            for scan in -400..=400 {
                let released = 5.0e-3 * c * hold * f64::from(scan) + 0.5 * start;
                let (force, mean, row) = law.on_line(start, released, c, true);
                *rows.entry(row).or_default() += 1;
                assert!(force * mean >= 0.0, "{start} {released} {c}: {force} {mean}");
                assert!(force.abs() <= hold, "{force}");
                match row {
                    SliderRow::Arrested => {
                        assert_eq!(mean, 0.5 * start);
                        assert!(force * start >= 0.0);
                        // The force that ends the step at zero slip.
                        assert_eq!(force, (released - 0.5 * start) / c);
                        carried_most_n = carried_most_n.max(force.abs());
                    }
                    SliderRow::Held => {
                        assert_eq!(mean, 0.0);
                        // The force that leaves no slip on the mean.
                        assert_eq!(force, (released / c).clamp(-hold, hold));
                    }
                    SliderRow::BrokeAway => {
                        assert_eq!(force, slide * released.signum());
                        assert!((mean - (released - c * force)).abs() <= 1.0e-15);
                    }
                    other => panic!("a stuck slider took row {other:?}"),
                }
            }
        }
    }
    eprintln!("stuck to a side that moved, rows over the sweep: {rows:?}");
    for row in [SliderRow::Arrested, SliderRow::Held, SliderRow::BrokeAway] {
        assert!(rows.get(&row).copied().unwrap_or(0) >= 1000, "{row:?}: {rows:?}");
    }
    // Carried with more than the sliding force could give, up to the
    // static limit.
    assert!(carried_most_n > 1.6 * slide && carried_most_n <= hold, "{carried_most_n}");
}

#[test]
fn the_slider_s_law_holds_row_by_row_and_conforms_on_its_line() {
    // Every row of the law's table, over start slips, compliances and a
    // scan of released slips. Asserted for each point: the loss is not
    // negative; the force is inside the static limit when stuck and inside
    // the sliding one when sliding; and within the sliding state the force
    // does not decrease as the released slip grows. A third of the points
    // also go through the conformance kit on a matrix that is not
    // symmetric.
    let law = slider();
    let (hold, slide) = (HOLDS * NORMAL_N, SLIDES * NORMAL_N);
    let mut rows: HashMap<SliderRow, u32> = HashMap::new();
    let mut checked = 0_u32;
    for &start in &[0.0, 2.0e-3, -2.0e-3, 0.4, -0.4, 6.0e-2] {
        for &c in &[2.5e-4, 1.0e-3, 4.0e-2] {
            let mut before = f64::NEG_INFINITY;
            for scan in -400..=400 {
                // Released slips across every threshold of this state, out
                // to twice the static limit on either side.
                let released = 5.0e-3 * c * hold * f64::from(scan) + 0.5 * start;
                let (force, mean, row) = law.on_line(start, released, c, false);
                *rows.entry(row).or_default() += 1;
                assert!(force * mean >= 0.0, "{start} {released} {c}: {force} {mean}");
                assert!((mean - (released - c * force)).abs() <= 1.0e-15 + 1.0e-12 * released.abs());
                if start == 0.0 {
                    assert!(force.abs() <= hold, "{force}");
                } else {
                    assert!(force.abs() <= slide, "{force}");
                    assert!(force >= before, "{start} {released} {c}: {force} after {before}");
                    before = force;
                }
                // Each row's own force and where it ends the step.
                let ends = 2.0 * mean - start;
                let way = start.signum();
                match row {
                    SliderRow::Held => {
                        assert!(start == 0.0 && mean == 0.0);
                        assert!((force - released / c).abs() <= 1.0e-12 * hold, "{force}");
                    }
                    SliderRow::BrokeAway => {
                        assert!(start == 0.0 && ends * released > 0.0);
                        assert_eq!(force, slide * released.signum());
                    }
                    SliderRow::Arrested => {
                        assert_eq!(ends, 0.0);
                        assert!(force * way >= 0.0 && force.abs() <= slide);
                    }
                    SliderRow::Slid => {
                        assert!(ends * start > 0.0);
                        assert_eq!(force, slide * way);
                    }
                    SliderRow::Released => assert!(force == 0.0 && ends * start <= 0.0),
                    SliderRow::Turned => {
                        assert!(mean == 0.0 && ends == -start);
                        assert!(force * way <= 0.0 && force.abs() <= slide);
                    }
                    SliderRow::Reversed => {
                        assert!(ends * start < 0.0);
                        assert_eq!(force, -slide * way);
                    }
                }
                if scan % 3 != 0 {
                    continue;
                }
                // On a line over both ends whose matrix gives this `c`.
                let matrix = [0.7 * c, -0.1 * c, 0.05 * c, 0.25 * c];
                assert!((matrix[0] - matrix[1] - matrix[2] + matrix[3] - c).abs() <= 1.0e-12 * c);
                let inputs = [0.3 + released, 0.3, 310.0];
                let state = Slip::sliding(start);
                let stale = Slip::sliding(7.0);
                let case = Case {
                    state: &state,
                    stale: &stale,
                    dt_s: 1.0e-3,
                    inputs: &inputs,
                };
                let seen = with_cx(|cx| check_load_line(&law, &case, &matrix, cx)).unwrap();
                let Observed::Stepped { outputs, receipt, .. } = seen.observed else {
                    panic!("the slider takes the step");
                };
                assert_eq!(outputs[1], -outputs[0]);
                assert!(receipt.dissipated_j >= 0.0);
                assert!(
                    (seen.flows[0] - seen.flows[1] - mean).abs() <= 1.0e-12 * (0.3 + mean.abs()),
                    "{:?} against {mean}",
                    seen.flows
                );
                checked += 1;
            }
        }
    }
    eprintln!("the slider's rows over the sweep: {rows:?}; {checked} points through the kit");
    for row in [
        SliderRow::Held,
        SliderRow::BrokeAway,
        SliderRow::Arrested,
        SliderRow::Slid,
        SliderRow::Released,
        SliderRow::Turned,
        SliderRow::Reversed,
    ] {
        assert!(rows.get(&row).copied().unwrap_or(0) >= 1000, "{row:?}: {rows:?}");
    }
    // What it refuses to be built from.
    assert!(Slider::new(0.0, 0.5, 0.3).is_err());
    assert!(Slider::new(10.0, 0.3, 0.5).is_err());
    assert!(Slider::new(10.0, f64::NAN, 0.3).is_err());
    assert!(Slider::new(10.0, 0.5, 0.0).is_err());
}

#[test]
fn a_thrown_mass_stops_and_its_arrest_time_error_falls_with_the_step() {
    // Thrown at 1.3 m/s against a friction force of 0.3 * 19.6 N, the mass
    // decelerates at 2.94 m/s^2 and stops after v0 / a = 0.442 s, having
    // travelled v0^2 / (2 a). While it slides the force is constant and
    // the midpoint rule is exact; the step in which it would overshoot
    // ends at zero slip, so the arrest falls on the first step boundary at
    // or after the true time and is late by less than a step. That bound
    // is what falls with the step; the lateness itself is the distance
    // from the true time up to the next multiple of the step.
    //
    // The slider's own slip is then exactly zero. The mass's velocity is
    // zero to the solve's tolerance and not to the bit: the slider keeps
    // its own record of the slip, and the two agree only as closely as the
    // loop was closed in each step. What is left flips sign every step
    // with a zero mean, so the position holds to the bit.
    let v0 = 1.3_f64;
    let deceleration = SLIDES * NORMAL_N / MASS_KG;
    let (stops_s, travels_m) = (v0 / deceleration, 0.5 * v0 * v0 / deceleration);
    let mut late = Vec::new();
    let mut left = Vec::new();
    for halvings in 0..3 {
        let dt = 0.01 / f64::from(1_u32 << halvings);
        let (mut system, mass, rubbing, body) = on_the_ground(v0, None);
        let (mut arrested_s, mut rests_at_m) = (f64::NAN, f64::NAN);
        let (mut steps_at_rest, mut iterations_at_rest) = (0_u32, 0_u32);
        with_cx(|cx| {
            for step in 1..=(100_u32 << halvings) {
                let report = system.step(dt, cx).unwrap();
                let motion = *system.state::<Motion>(mass).unwrap();
                let slip = system.state::<Slip>(rubbing).unwrap().slip_m_s();
                if slip == 0.0 {
                    assert!(motion.velocity_m_s.abs() <= 1.0e-12 * v0, "{motion:?}");
                    if arrested_s.is_nan() {
                        (arrested_s, rests_at_m) = (f64::from(step) * dt, motion.position_m);
                    } else {
                        assert_eq!(motion.position_m, rests_at_m, "step {step}");
                        steps_at_rest += 1;
                        iterations_at_rest += report.iterations;
                    }
                } else {
                    assert!(arrested_s.is_nan(), "sliding again at step {step}: {motion:?}");
                    assert!(motion.velocity_m_s > 0.0 && slip < 0.0, "{motion:?} {slip}");
                }
            }
        });
        assert!(steps_at_rest > 10, "{steps_at_rest}");
        let lateness = arrested_s - stops_s;
        assert!(lateness >= 0.0 && lateness < dt, "dt {dt}: {lateness}");
        late.push(lateness);
        left.push(system.state::<Motion>(mass).unwrap().velocity_m_s.abs());
        // It stops within a step's worth of the true distance.
        let off_m = rests_at_m - travels_m;
        assert!(off_m.abs() <= deceleration * dt * dt, "dt {dt}: {off_m}");
        // Every joule it had is in the body.
        let warmed_j = received_j(&system, body);
        let had_j = 0.5 * MASS_KG * v0 * v0;
        assert!((warmed_j - had_j).abs() <= 1.0e-12 * had_j, "{warmed_j} vs {had_j}");
        assert!(system.ledger().energy_residual_j().abs() <= 1.0e-12 * had_j);
        // At rest after a slide the flipping remainder costs an update a
        // step, where a body that never moved costs none.
        assert!(iterations_at_rest <= steps_at_rest + 4, "{iterations_at_rest} {steps_at_rest}");
    }
    eprintln!(
        "a thrown mass: arrest late by {late:?} s at dt, dt/2, dt/4 (true {stops_s} s); velocity \
         left {left:?} m/s"
    );
}

#[test]
fn a_pushed_mass_is_held_up_to_the_static_limit_and_slides_beyond_it() {
    let (holds_n, slides_n) = (HOLDS * NORMAL_N, SLIDES * NORMAL_N);
    let dt = 1.0e-3;
    // Just inside the limit: nothing moves, to the bit, and after the
    // first step nothing is solved.
    let (mut held, mass, rubbing, body) = on_the_ground(0.0, Some(0.999 * holds_n));
    with_cx(|cx| {
        for step in 0..500 {
            let report = held.step(dt, cx).unwrap();
            let motion = *held.state::<Motion>(mass).unwrap();
            assert_eq!((motion.position_m, motion.velocity_m_s), (0.0, 0.0), "step {step}");
            if step > 0 {
                assert_eq!(report.iterations, 0, "step {step}");
            }
        }
    });
    assert_eq!(held.state::<Slip>(rubbing).unwrap().slip_m_s(), 0.0);
    assert_eq!(*held.state::<f64>(body).unwrap(), 300.0);
    // The slider holds what the push asks of it, and what its port is
    // published to carry is the flow on the line, none, not the released
    // flow the solve iterated on (dt F / (2 m), about 2e-3 m/s here).
    let at_mass = held.port(rubbing, "a").unwrap();
    assert!(held.last_input(at_mass)[0].abs() <= 1.0e-15, "{:?}", held.last_input(at_mass));
    assert!((held.last_output(at_mass)[0].abs() - 0.999 * holds_n).abs() <= 1.0e-12 * holds_n);

    // Just outside it: the mass breaks away and accelerates against the
    // sliding force, which is smaller than the one that held it.
    let push_n = 1.001 * holds_n;
    let (mut sliding, mass, _, body) = on_the_ground(0.0, Some(push_n));
    let steps = 500_u32;
    with_cx(|cx| {
        for _ in 0..steps {
            sliding.step(dt, cx).unwrap();
        }
    });
    let elapsed = f64::from(steps) * dt;
    let acceleration = (push_n - slides_n) / MASS_KG;
    let motion = *sliding.state::<Motion>(mass).unwrap();
    let (speed, distance) = (acceleration * elapsed, 0.5 * acceleration * elapsed * elapsed);
    assert!((motion.velocity_m_s - speed).abs() <= 1.0e-12 * speed, "{motion:?} vs {speed}");
    assert!((motion.position_m - distance).abs() <= 1.0e-12 * distance, "{motion:?}");
    // The loss is the sliding force over the distance.
    let warmed_j = received_j(&sliding, body);
    assert!((warmed_j - slides_n * distance).abs() <= 1.0e-12 * slides_n * distance, "{warmed_j}");
}

#[test]
fn a_step_whose_drive_alone_reverses_the_body_has_no_friction_and_no_heat() {
    // The fifth row of the law. A mass sliding at 1 mm/s is pushed back
    // with 1.5 m v / dt, enough to reverse it within the step on its own.
    // A friction force of the old direction would have a negative loss on
    // the step's mean slip, and the mean does not yet call for one of the
    // new direction: the row gives the friction of this one step up to
    // keep the first law exact. No heat reaches the body, and the slider
    // enters the next step sliding the other way.
    let (v, dt) = (1.0e-3_f64, 1.0e-3_f64);
    let push_n = -1.5 * MASS_KG * v / dt;
    let (mut system, mass, rubbing, body) = on_the_ground(v, Some(push_n));
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s(), -v);
    with_cx(|cx| system.step(dt, cx)).unwrap();
    assert_eq!(*system.state::<f64>(body).unwrap(), 300.0);
    assert_eq!(system.last_output(system.port(rubbing, "a").unwrap()), [0.0]);
    let motion = *system.state::<Motion>(mass).unwrap();
    assert!((motion.velocity_m_s + 0.5 * v).abs() <= 1.0e-12 * v, "{motion:?}");
    let slip = system.state::<Slip>(rubbing).unwrap().slip_m_s();
    assert!((slip - 0.5 * v).abs() <= 1.0e-12 * v, "{slip}");
    // From there friction acts again, against the new direction: the
    // push, 3 N, is inside what the slider holds, and the body is brought
    // to rest and kept there.
    with_cx(|cx| {
        for _ in 0..3 {
            system.step(dt, cx).unwrap();
        }
    });
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s(), 0.0);
    assert!(*system.state::<f64>(body).unwrap() > 300.0);
    let rest = *system.state::<Motion>(mass).unwrap();
    assert!(rest.velocity_m_s.abs() <= 1.0e-12 * v, "{rest:?}");
}

#[test]
fn a_mass_thrown_onto_another_shares_its_momentum_and_the_rest_is_heat() {
    // 2 kg at 1.3 m/s on a 3 kg slab at rest, a slider between them, the
    // slab free. The slider's line spans both ends. Their relative speed
    // falls at mu_k N (1/m1 + 1/m2) and they end together at the velocity
    // of their centre of mass; the kinetic energy of the relative motion
    // is the heat.
    let (m1, m2, v0) = (MASS_KG, 3.0_f64, 1.3_f64);
    let mut builder = SystemBuilder::new();
    let add_mass = |builder: &mut SystemBuilder, name: &str, kg: f64, velocity_m_s: f64| {
        builder
            .add_with_tangent(
                name,
                Mass::new(kg).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s,
                },
            )
            .unwrap()
    };
    let top = add_mass(&mut builder, "top", m1, v0);
    let slab = add_mass(&mut builder, "slab", m2, 0.0);
    let rubbing = builder
        .add_on_line("slider", slider(), Slip::sliding(-v0))
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .common_flow(port(&builder, top, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    builder
        .common_flow(port(&builder, slab, "p"), &[port(&builder, rubbing, "b")])
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));

    let dt = 1.0e-3;
    let together = m1 * v0 / (m1 + m2);
    let closing = SLIDES * NORMAL_N * (1.0 / m1 + 1.0 / m2);
    let mut joined_s = f64::NAN;
    with_cx(|cx| {
        for step in 1..=600_u32 {
            system.step(dt, cx).unwrap();
            let (a, b) = (
                system.state::<Motion>(top).unwrap().velocity_m_s,
                system.state::<Motion>(slab).unwrap().velocity_m_s,
            );
            // No outside force: the momentum is the one it was thrown with.
            assert!((m1 * a + m2 * b - m1 * v0).abs() <= 1.0e-12 * m1 * v0, "step {step}");
            if joined_s.is_nan() && system.state::<Slip>(rubbing).unwrap().slip_m_s() == 0.0 {
                joined_s = f64::from(step) * dt;
            }
        }
    });
    let lateness = joined_s - v0 / closing;
    assert!(lateness >= 0.0 && lateness < dt, "{lateness}");
    for mass in [top, slab] {
        let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
        assert!((velocity - together).abs() <= 1.0e-12 * v0, "{velocity} vs {together}");
    }
    let warmed_j = received_j(&system, body);
    let relative_j = 0.5 * m1 * v0 * v0 - 0.5 * (m1 + m2) * together * together;
    assert!((warmed_j - relative_j).abs() <= 1.0e-12 * relative_j, "{warmed_j} vs {relative_j}");
}

/// A block of 1 kg on a belt moving at 0.1 m/s, held to a wall by a spring
/// of 100 N/m, under 10 N with coefficients 0.5 and 0.3. Returns the times
/// [s] at which the block is caught by the belt again, over `steps` steps
/// of `dt`.
fn caught_on_the_belt(dt: f64, steps: u32) -> Vec<f64> {
    let (belt_m_s, stiffness, normal_n) = (0.1_f64, 100.0_f64, 10.0_f64);
    let mut builder = SystemBuilder::new();
    let block = builder
        .add_with_tangent(
            "block",
            Mass::new(1.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: belt_m_s,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(stiffness, 0.0).unwrap(), 0.0)
        .unwrap();
    let rubbing = builder
        .add_on_line("slider", Slider::new(normal_n, 0.5, 0.3).unwrap(), Slip::default())
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .common_flow(
            port(&builder, block, "p"),
            &[port(&builder, spring, "a"), port(&builder, rubbing, "a")],
        )
        .unwrap();
    builder
        .boundary(port(&builder, spring, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    // The belt is the slider's other end, moving with the block at first:
    // as a junction member the slider's `a` receives minus the block's
    // velocity, so the belt is given with the same sign.
    builder
        .boundary(
            port(&builder, rubbing, "b"),
            Box::new(move |_, _, values| values[0] = -belt_m_s),
        )
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    let mut caught = Vec::new();
    let mut stuck = true;
    with_cx(|cx| {
        for step in 1..=steps {
            system.step(dt, cx).unwrap();
            let now = system.state::<Slip>(rubbing).unwrap().slip_m_s() == 0.0;
            if now && !stuck {
                caught.push(f64::from(step) * dt);
            }
            stuck = now;
        }
    });
    // It slid: the body took the loss.
    assert!(*system.state::<f64>(body).unwrap() > 300.0);
    caught
}

#[test]
fn a_block_on_a_belt_sticks_and_slips_with_the_closed_form_s_period() {
    // Stuck, the block rides the belt until the spring pulls with the
    // static limit, at x_s = mu_s N / k. Sliding, it swings about
    // x_k = mu_k N / k as a free oscillator, from x_k + d at the belt's
    // speed (d = x_s - x_k) round to x_k - d at the belt's speed again,
    // where the belt catches it. One period is the ride 2 d / V and the
    // swing 2 (pi - atan(d w / V)) / w.
    //
    // The steps in which the block breaks away and is caught are first
    // order, each off by a part of a step, so the period is off by less
    // than a step: a bound that falls with the step, not a smooth error.
    // (The discrete motion repeats after a whole number of steps, so the
    // measured period is a multiple of the step and its error moves in
    // jumps.)
    let (belt, omega, d) = (0.1_f64, 10.0_f64, (0.5 - 0.3) * 10.0 / 100.0);
    let period = 2.0 * d / belt + 2.0 * (PI - (d * omega / belt).atan()) / omega;
    let mut errors = Vec::new();
    for halvings in 0..5 {
        let dt = 2.0e-3 / f64::from(1_u32 << halvings);
        let caught = caught_on_the_belt(dt, 3_000 << halvings);
        assert!(caught.len() >= 6, "{caught:?}");
        // Over the cycles after the first, which starts from another place.
        let cycles = (caught.len() - 1) as f64;
        let measured = (caught[caught.len() - 1] - caught[0]) / cycles;
        let error = (measured - period).abs();
        assert!(error < dt, "dt {dt}: period {measured} against {period}");
        errors.push(error);
    }
    eprintln!("stick-slip period {period} s: off by {errors:?} s at dt = 2 ms and four halvings");
}

#[test]
fn a_slider_between_given_velocities_warms_a_body_through_a_loop() {
    // Both ends of the slider are prescribed, so there is no line to find;
    // but its heat port shares a body's temperature, which is a loop. The
    // slider is then differenced as a function of that temperature (its
    // loss does not depend on it), and the body takes force times slip.
    let (slip, dt, steps) = (0.2_f64, 1.0e-3_f64, 400_u32);
    let mut builder = SystemBuilder::new();
    let rubbing = builder
        .add_on_line("slider", slider(), Slip::stuck())
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "a"), Box::new(move |_, _, values| values[0] = slip))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });
    let loss_j = SLIDES * NORMAL_N * slip * dt * f64::from(steps);
    let warmed_j = received_j(&system, body);
    assert!((warmed_j - loss_j).abs() <= 1.0e-12 * loss_j, "{warmed_j} vs {loss_j}");
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s(), slip);
}

#[test]
fn a_run_with_a_line_and_its_replay_agree_to_the_bit() {
    // The line of a step is found in that step: the compliance in the
    // drivers and the released flows are set anew before anything is
    // judged. A copy that takes a step of another length and is restored
    // before every step takes the same iterations and reaches the same
    // bits, sliding, at the arrest and at rest.
    let (mut straight, mass, rubbing, body) = on_the_ground(1.3, None);
    let (mut detoured, _, _, _) = on_the_ground(1.3, None);
    let dt = 0.01;
    with_cx(|cx| {
        for step in 0..80 {
            let report = straight.step(dt, cx).unwrap();
            let before = detoured.checkpoint().unwrap();
            detoured.step(0.37 * dt, cx).unwrap();
            detoured.restore(&before).unwrap();
            let again = detoured.step(dt, cx).unwrap();
            assert_eq!(report.iterations, again.iterations, "step {step}");
            let (a, b) = (
                *straight.state::<Motion>(mass).unwrap(),
                *detoured.state::<Motion>(mass).unwrap(),
            );
            assert_eq!(a.velocity_m_s.to_bits(), b.velocity_m_s.to_bits(), "step {step}");
            assert_eq!(a.position_m.to_bits(), b.position_m.to_bits(), "step {step}");
            assert_eq!(
                straight.state::<Slip>(rubbing).unwrap().slip_m_s().to_bits(),
                detoured.state::<Slip>(rubbing).unwrap().slip_m_s().to_bits(),
                "step {step}"
            );
            assert_eq!(
                straight.state::<f64>(body).unwrap().to_bits(),
                detoured.state::<f64>(body).unwrap().to_bits(),
                "step {step}"
            );
        }
    });
    // The run covered the slide, the arrest and steps at rest.
    assert_eq!(straight.state::<Slip>(rubbing).unwrap().slip_m_s(), 0.0);
    assert!(straight.state::<Motion>(mass).unwrap().position_m > 0.2);
}

/// The row of `law` that a step of the kernel is a point of. The slider
/// started the step as `last`, on a line of compliance `answered`, and the
/// step published `force` and `flow` at its port `a`. The released slip is
/// rebuilt from those two, the slip that crossed and what the force took
/// away, and the law must give that force (to 1e-9 N) and that slip (to
/// 1e-12 m/s) for it.
///
/// This checks the rows and their forces: held with no slip and a force
/// inside the limit, sliding with the sliding force along the slip, each
/// reversal row's own force. It does not check `answered`: a released
/// slip rebuilt from a held step satisfies the law for any compliance, and
/// only the thresholds between rows move with it. The compliance is
/// pinned by the probe test, where a mass on a spring is one of the
/// cases.
fn row_of_the_step(law: &Slider, last: Slip, answered: f64, force: f64, flow: f64) -> SliderRow {
    let released = flow + answered * force;
    let (wanted, mean, row) = law.on_line(last.slip_m_s(), released, answered, last.is_stuck());
    assert!(
        (force - wanted).abs() <= 1.0e-9 && (flow - mean).abs() <= 1.0e-12,
        "{row:?}: {force} against {wanted}, {flow} against {mean}"
    );
    row
}

/// Counts of rows by name, in the order of the names.
fn tally(rows: HashMap<SliderRow, u32>) -> Vec<(String, u32)> {
    let mut reached: Vec<_> = rows
        .into_iter()
        .map(|(row, count)| (format!("{row:?}"), count))
        .collect();
    reached.sort();
    reached
}

#[test]
#[allow(clippy::too_many_lines)]
fn a_mass_on_a_spring_loses_two_friction_offsets_a_half_cycle_and_stops_in_the_dead_band() {
    // 1 kg on a spring of 100 N/m, sliding on the ground under 10 N with
    // coefficients 0.5 and 0.3, thrown from the spring's rest point at
    // sqrt(22) m/s. Its first turning point A has k A^2 / 2 + mu_k N A =
    // m v0^2 / 2: A = 0.44 m. In each half cycle it swings about a centre
    // shifted by mu_k N / k = 0.03 m against its motion, so the turning
    // points fall by 0.06 m: -0.38, 0.32, -0.26, 0.20, -0.14, 0.08. At
    // -0.02 the spring pulls with 2 N, inside the 5 N the slider holds,
    // and it stays there. Every turning point is a reversal under a drive
    // that is not friction: the rows of the law that no other test of the
    // kernel reaches.
    let (stiffness, v0) = (100.0_f64, 22.0_f64.sqrt());
    let turns = [0.44, -0.38, 0.32, -0.26, 0.20, -0.14, 0.08, -0.02];
    let mut worst = Vec::new();
    for halvings in 0..2 {
        let dt = 1.0e-3 / f64::from(1_u32 << halvings);
        let mut builder = SystemBuilder::new();
        let block = builder
            .add_with_tangent(
                "block",
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: v0,
                },
            )
            .unwrap();
        let spring = builder
            .add_with_tangent("spring", Spring::new(stiffness, 0.0).unwrap(), 0.0)
            .unwrap();
        let rubbing = builder
            .add_on_line("slider", Slider::new(10.0, 0.5, 0.3).unwrap(), Slip::sliding(-v0))
            .unwrap();
        let body = builder
            .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
            .unwrap();
        builder
            .common_flow(
                port(&builder, block, "p"),
                &[port(&builder, spring, "a"), port(&builder, rubbing, "a")],
            )
            .unwrap();
        for held in [port(&builder, spring, "b"), port(&builder, rubbing, "b")] {
            builder.boundary(held, Box::new(|_, _, values| values[0] = 0.0)).unwrap();
        }
        let warmth = port(&builder, rubbing, "heat");
        builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(warmth, &[300.0]));
        let at_block = system.port(rubbing, "a").unwrap();
        let law = Slider::new(10.0, 0.5, 0.3).unwrap();
        // What the mass and the spring answer a force at the slider's port
        // with: g / (1 + g k dt / 2), g = dt / (2 m).
        let answered = 0.5 * dt / (1.0 + 0.25 * dt * dt * stiffness);
        let mut last = *system.state::<Slip>(rubbing).unwrap();
        let mut rows: HashMap<SliderRow, u32> = HashMap::new();

        // A turning point is where the velocity changes sign; the mass is
        // within a step's travel at a step's speed change of it.
        let mut found = Vec::new();
        let mut before = Motion {
            position_m: 0.0,
            velocity_m_s: v0,
        };
        with_cx(|cx| {
            for step in 0..(4_000_u32 << halvings) {
                system.step(dt, cx).unwrap();
                // Every step is a point of the slider's law on the line
                // its neighbours present, the reversals among them.
                let (force, flow) = (system.last_output(at_block)[0], system.last_input(at_block)[0]);
                let row = row_of_the_step(&law, last, answered, force, flow);
                *rows.entry(row).or_default() += 1;
                last = *system.state::<Slip>(rubbing).unwrap();
                // A released step has no friction in it: the slider ends it
                // sliding, and the next step is judged from sliding.
                assert!(
                    row != SliderRow::Released || !last.is_stuck(),
                    "dt {dt}, step {step}"
                );
                let now = *system.state::<Motion>(block).unwrap();
                if before.velocity_m_s != 0.0 && now.velocity_m_s * before.velocity_m_s <= 0.0 {
                    let turned = if now.position_m.abs() > before.position_m.abs() {
                        now.position_m
                    } else {
                        before.position_m
                    };
                    // Once more only if it has travelled since the last
                    // one: at rest what is left of the velocity flips sign
                    // every step at the level of rounding.
                    if found.last().is_none_or(|last: &f64| (turned - last).abs() > 1.0e-3) {
                        found.push(turned);
                    }
                }
                before = now;
            }
        });
        assert_eq!(found.len(), turns.len(), "{found:?}");
        let off = found
            .iter()
            .zip(turns)
            .map(|(got, want)| (got - want).abs())
            .fold(0.0_f64, f64::max);
        worst.push(off);
        // It stopped, stuck, inside the dead band, and what it was thrown
        // with is in the spring and in the heat.
        let rest = *system.state::<Motion>(block).unwrap();
        assert!(system.state::<Slip>(rubbing).unwrap().is_stuck());
        assert!(rest.velocity_m_s.abs() <= 1.0e-12, "{rest:?}");
        assert!((rest.position_m - turns[7]).abs() <= 2.0 * off + 1.0e-12, "{rest:?}");
        let stretch = *system.state::<f64>(spring).unwrap();
        assert!((stretch.abs() - rest.position_m.abs()).abs() <= 1.0e-12, "{stretch}");
        let stored_j = 0.5 * stiffness * stretch * stretch;
        let thrown_j = 0.5 * v0 * v0;
        let heat_j = received_j(&system, body);
        assert!((heat_j + stored_j - thrown_j).abs() <= 1.0e-12 * thrown_j, "{heat_j} {stored_j}");
        // The rows of a reversal are among them.
        for row in [SliderRow::Released, SliderRow::Turned, SliderRow::Reversed] {
            assert!(rows.get(&row).copied().unwrap_or(0) >= 1, "{row:?}: {rows:?}");
        }
        let reached = tally(rows);
        eprintln!("a mass on a spring with friction, dt {dt}: rows through the kernel {reached:?}");
    }
    eprintln!(
        "a mass on a spring with friction: turning points off the closed form by {worst:?} m at \
         most, at steps of 1 and 0.5 ms"
    );
    // A reversal is a step in which friction of one direction gives way
    // to the other; what the law gets wrong in it is at most the friction
    // of that step, which moves the mass by (mu_k N / m) dt^2. Eight
    // reversals may add up. (Measured: 4.3e-6 and 1.1e-6 m.)
    for (off, dt) in worst.iter().zip([1.0e-3_f64, 0.5e-3]) {
        assert!(*off <= 8.0 * 3.0 * dt * dt, "{worst:?}");
    }
}

#[test]
fn a_small_jump_under_a_push_is_held_on_the_mean_and_caught_in_the_next_step() {
    // A mass pushed with 6 N and held, on a side that jumps by 1 mm/s in
    // the direction of the push. The step of the jump starts with the mass
    // 1 mm/s behind. Ending it at zero slip would leave the slider -4 N to
    // give, a force that holds the mass back while the slip asks for one
    // that pulls it along: a negative loss, which the slider does not
    // book. It holds on the mean instead: the mass moves at the side's
    // speed on the mean of the step, with -2 N, and ends it 1 mm/s ahead.
    // The next step starts with that slip, and now the force that ends it
    // at zero slip, -8 N, is against the slip and inside the 9.8 N: the
    // mass is caught, 4 microjoules are lost, and from then on it is held
    // with -6 N. (The stuck slider's hold on the mean with a start slip:
    // the row of its law that no other test reaches through the kernel.)
    let (push_n, jump_m_s, dt) = (6.0_f64, 1.0e-3_f64, 1.0e-3_f64);
    let side = Moved::Boundary(Box::new(move |t0, _| if t0 > 1.5e-3 { -jump_m_s } else { 0.0 }));
    let (mut system, mass, rubbing, body) =
        carried(side, Some(Waveform::constant(push_n).unwrap()));
    let at_mass = system.port(rubbing, "a").unwrap();
    // What the mass needs to gain the jump in one step: 2 N.
    let catching_n = MASS_KG * jump_m_s / dt;
    let forces = [
        -push_n,
        -push_n,
        2.0 * catching_n - push_n,
        -catching_n - push_n,
        -push_n,
        -push_n,
    ];
    let velocities = [0.0, 0.0, 2.0 * jump_m_s, jump_m_s, jump_m_s, jump_m_s];
    assert_eq!((forces[2], forces[3]), (-2.0, -8.0));
    with_cx(|cx| {
        for (step, (wanted_n, wanted_m_s)) in forces.into_iter().zip(velocities).enumerate() {
            system.step(dt, cx).unwrap();
            let force = system.last_output(at_mass)[0];
            let velocity = system.state::<Motion>(mass).unwrap().velocity_m_s;
            assert!((force - wanted_n).abs() <= 1.0e-9, "step {step}: {force}");
            assert!((velocity - wanted_m_s).abs() <= 1.0e-15, "step {step}: {velocity}");
            let state = *system.state::<Slip>(rubbing).unwrap();
            assert!(state.is_stuck(), "step {step}");
            // The record after the step of the jump: 1 mm/s the other way.
            let recorded = if step == 2 { -jump_m_s } else { 0.0 };
            assert!((state.slip_m_s() - recorded).abs() <= 1.0e-15, "step {step}: {state:?}");
        }
    });
    // The catch: 8 N over half the slip, for one step.
    let heat_j = received_j(&system, body);
    let lost_j = dt * (catching_n + push_n) * 0.5 * jump_m_s;
    assert!((heat_j - lost_j).abs() <= 1.0e-12 * lost_j, "{heat_j} against {lost_j}");
    assert!((lost_j - 4.0e-6).abs() <= 1.0e-18);
}

/// Three masses of 2 kg at rest in a row on two springs of 100 N/m, the
/// first under a push of half the static limit, 4.9 N, and kept where it
/// is by the slider or, as the twin without a load line, by a spring of
/// 4.9 N/m to a wall that starts stretched by the metre that carries that
/// push. Returns the system and the third mass.
fn at_rest_beside(held_by_the_slider: bool) -> (System, NodeId) {
    let push = 0.5 * HOLDS * NORMAL_N;
    let at_rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let mut builder = SystemBuilder::new();
    let masses: Vec<NodeId> = ["mass", "second", "third"]
        .iter()
        .map(|name| builder.add_with_tangent(name, Mass::new(MASS_KG).unwrap(), at_rest).unwrap())
        .collect();
    let pusher = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(push).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let springs: Vec<NodeId> = ["near", "far"]
        .iter()
        .map(|name| builder.add_with_tangent(name, Spring::new(100.0, 0.0).unwrap(), 0.0).unwrap())
        .collect();
    let mut first = vec![port(&builder, pusher, "p"), port(&builder, springs[0], "a")];
    if held_by_the_slider {
        let rubbing = builder.add_on_line("slider", slider(), Slip::stuck()).unwrap();
        first.push(port(&builder, rubbing, "a"));
        builder
            .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        builder
            .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
    } else {
        let wall = builder
            .add_with_tangent("wall", Spring::new(push, 0.0).unwrap(), -1.0)
            .unwrap();
        first.push(port(&builder, wall, "a"));
        builder
            .boundary(port(&builder, wall, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
    }
    // The junctions in the order of the masses: the third mass's is the
    // third connection.
    builder.common_flow(port(&builder, masses[0], "p"), &first).unwrap();
    builder
        .common_flow(
            port(&builder, masses[1], "p"),
            &[port(&builder, springs[0], "b"), port(&builder, springs[1], "a")],
        )
        .unwrap();
    builder
        .common_flow(port(&builder, masses[2], "p"), &[port(&builder, springs[1], "b")])
        .unwrap();
    (builder.seal().unwrap(), masses[2])
}

#[test]
fn a_leak_beside_a_held_slider_is_refused_at_the_rounding_of_the_block_s_largest_force() {
    // In a block with a load line a wire that is not a flow on the line is
    // allowed the rounding of the largest wire of its kind in the block,
    // 16 eps of it, and the audit is told that for it. The planted
    // negative for that term, at the junction farthest from the line,
    // which has nothing of its own: three masses in a row on two springs,
    // the first held by the slider against 4.9 N, the other two at rest.
    // The largest force wire of the block is the held mass's, the push
    // and the friction: 9.8 N of terms.
    //
    // A net force F is planted on the third mass's input, where the solve
    // starts, with the solve told to take a wire that is open by four
    // times its own size. The mass then moves at c F on the mean with a
    // force on it that nothing returns: its junction is left with
    // c dt F^2, and the audit, which multiplies what it is told for a
    // wire by the wire's paired output, allows it the rounding term times
    // that motion, c dt F times 16 eps of 9.8 N. So what is left is F
    // over that term times what is allowed: the verdict is the audit's,
    // and it turns where F is the rounding of the block's largest force.
    let (dt, push) = (1.0e-3_f64, 0.5 * HOLDS * NORMAL_N);
    let c = 0.5 * dt / MASS_KG;
    let loose = SolveConfig {
        relative_tolerance: 4.0,
        ..SolveConfig::default()
    };
    let settled_after = |held_by_the_slider: bool| {
        let (mut system, third) = at_rest_beside(held_by_the_slider);
        with_cx(|cx| {
            for _ in 0..3 {
                system.step(dt, cx).unwrap();
            }
        });
        let far = system.port(third, "p").unwrap();
        let settled = system.checkpoint().unwrap();
        // What is on the third mass is what the rounding of the first
        // one's balance leaves there: nothing to speak of.
        let settled_n = system.last_input(far)[0];
        assert!(settled_n.abs() <= 1.0e-20, "{settled_n}");
        (system, far, settled, settled_n)
    };
    let planted = |system: &mut System, far, settled: &(Checkpoint, f64), force_n: f64| {
        system.restore(&settled.0).unwrap();
        system.set_config(loose);
        assert!(system.guess(far, &[settled.1 + force_n]));
        with_cx(|cx| system.step(dt, cx))
    };

    let (mut system, far, settled, settled_n) = settled_after(true);
    let settled = (settled, settled_n);
    // The term by name, and read back from the last step taken.
    let kind_n = 16.0 * f64::EPSILON * 2.0 * push;
    let told_n = system.last_allowance(far)[0];
    assert!((told_n - kind_n).abs() <= 0.05 * kind_n, "{told_n} against {kind_n}");
    // A hundred times that term is refused by the third junction's audit,
    // by that factor; the tolerance it reports is the term times the
    // motion; and no step is counted.
    let refusal = planted(&mut system, far, &settled, 100.0 * kind_n).expect_err("refused");
    let SystemError::Interface {
        connection,
        defect_j,
        tolerance_j,
    } = refusal
    else {
        panic!("expected the audit's refusal, got {refusal}");
    };
    assert_eq!(connection, 2);
    let named_j = dt * kind_n * c * (100.0 * kind_n);
    assert!((tolerance_j - named_j).abs() <= 0.05 * named_j, "{tolerance_j} against {named_j}");
    let ratio = defect_j.abs() / tolerance_j;
    assert!(ratio > 80.0 && ratio < 120.0, "{ratio}");
    assert_eq!(system.steps(), 3);
    // A quarter of the term is taken, with no update.
    let report = planted(&mut system, far, &settled, 0.25 * kind_n).unwrap();
    assert_eq!(report.iterations, 0);
    assert_eq!(system.steps(), 4);

    // Only in a block with a load line. The twin without one, the first
    // mass kept where it is by a spring to a wall under the same 9.8 N of
    // terms: the third mass's wire is told nothing of the kind, and the
    // quarter that was taken above is refused, a thousand times over.
    let (mut twin, far, settled, settled_n) = settled_after(false);
    let settled = (settled, settled_n);
    assert_eq!(twin.last_allowance(far), [0.0]);
    let refusal = planted(&mut twin, far, &settled, 0.25 * kind_n).expect_err("refused");
    let SystemError::Interface {
        connection: 2,
        defect_j: twin_defect_j,
        tolerance_j: twin_j,
    } = refusal
    else {
        panic!("expected the third junction's refusal, got {refusal}");
    };
    assert!(twin_defect_j.abs() > 1000.0 * twin_j, "{twin_defect_j} against {twin_j}");
    eprintln!(
        "a junction at rest beside a held slider: a leak of 100 times the rounding of the \
         block's largest force is refused by a factor of {ratio:.1}, a quarter of it is taken; \
         beside a spring instead of the slider that quarter is refused by a factor of {:e}",
        twin_defect_j.abs() / twin_j
    );
}

/// A mass at rest on the ground, or on a belt that moves at `belt_m_s`
/// (as the slider's `b` receives it), on a spring to a wall that starts
/// stretched by `stretch_m`; the spring with its tangent declared or
/// left to be differenced. Returns the system, the mass and the slider.
fn on_a_spring_to_a_wall(
    declared: bool,
    stiffness: f64,
    stretch_m: f64,
    belt_m_s: f64,
    slip: Slip,
    push: Option<Waveform>,
) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let rubbing = builder.add_on_line("slider", slider(), slip).unwrap();
    let law = Spring::new(stiffness, 0.0).unwrap();
    let spring = if declared {
        builder.add_with_tangent("spring", law, stretch_m).unwrap()
    } else {
        builder.add("spring", law, stretch_m).unwrap()
    };
    let mut members = vec![port(&builder, rubbing, "a"), port(&builder, spring, "a")];
    if let Some(push) = push {
        let pusher = builder
            .add_with_tangent("push", Force::new(push), SourceAccount::default())
            .unwrap();
        members.push(port(&builder, pusher, "p"));
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(move |_, _, values| values[0] = belt_m_s))
        .unwrap();
    builder
        .boundary(port(&builder, spring, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    (builder.seal().unwrap(), mass, rubbing)
}

#[test]
fn break_away_beside_a_differenced_neighbour_at_rest_is_solved() {
    // The line's matrix is what the neighbours answer at the step's start
    // as the kernel knows it, and a neighbour without a declared tangent
    // is differenced. While the mass was at rest to rounding, the step of
    // that difference, `sqrt(eps)` of a velocity at rounding, moved the
    // spring's force by less than its own rounding: a spring of 400 N/m
    // gave the line no slope, the slider was handed the mass's dt / (2 m)
    // = 0.00125 where the mass and the spring answer with 0.0012484, and
    // the first step that should break away did not converge (the force
    // wire flipped by the static limit less the sliding force, 3.92 N: with
    // a matrix too large there is a band of released slips at break-away
    // with no point to settle on). A stated limit until a differenced
    // column was resolved against the rounding of its outputs. (Found by
    // the independent review of 74529fda7f.)
    //
    // The spring holds 2 N of the push; the push rises from 11.7 N to
    // 11.9 N in 4 s and passes the 11.8 N the slider and the spring hold
    // together in step 400 (5 ms). The spring declared or differenced, the
    // mass breaks away in that step and slides on.
    let dt = 5.0e-3;
    let push = || Waveform::table(vec![(0.0, 11.7), (4.0, 11.9)]).unwrap();
    let mut ends = Vec::new();
    for declared in [true, false] {
        let (mut system, mass, rubbing) =
            on_a_spring_to_a_wall(declared, 400.0, -0.005, 0.0, Slip::stuck(), Some(push()));
        with_cx(|cx| {
            for step in 0..400 {
                system.step(dt, cx).unwrap();
                assert!(system.state::<Slip>(rubbing).unwrap().is_stuck(), "step {step}");
            }
            system.step(dt, cx).unwrap();
            assert!(
                !system.state::<Slip>(rubbing).unwrap().is_stuck(),
                "declared {declared}"
            );
            for _ in 401..1000 {
                system.step(dt, cx).unwrap();
            }
        });
        ends.push(*system.state::<Motion>(mass).unwrap());
    }
    eprintln!("after 1000 steps, the spring declared and differenced: {ends:?}");
}

#[test]
fn a_starting_slip_that_contradicts_the_ports_is_outside_the_claim() {
    // A stated limit. Nothing can check a starting slip against what is
    // attached. A mass at rest on a belt that already moves at 0.2 m/s
    // slips at 0.2 m/s; declared `Slip::stuck()`, the slider's record is
    // off by that slip for ever and alternates in sign, and its force is a
    // difference of flows of 0.2 m/s divided by c, which moves in steps of
    // 1.1e-13 N. Up to ef1631e0d2 such a run went on with its wrong
    // record. Since the solve judges a flow on the line by what crosses,
    // it ends when the mass has come to rest on its spring and those
    // steps are more than its wires may be left with: not converged. (At
    // step 3569 of 6000 on the build this was written on; found by the
    // independent review of 74529fda7f.) With the slip it has, the run is
    // taken to its end.
    let dt = 1.0e-3;
    let run = |slip: Slip| {
        let (mut system, _, _) = on_a_spring_to_a_wall(true, 2000.0, 0.0, -0.2, slip, None);
        let refusal = with_cx(|cx| (0..6000).find_map(|_| system.step(dt, cx).err()));
        (system.steps(), refusal)
    };
    let (steps, refusal) = run(Slip::sliding(0.2));
    assert!(refusal.is_none() && steps == 6000, "{refusal:?}");
    let (steps, refusal) = run(Slip::stuck());
    let refusal = refusal.expect("a record that contradicts the ports does not last");
    assert!(matches!(refusal, SystemError::NotConverged { .. }), "{refusal}");
    eprintln!("a mass at rest declared stuck on a belt at 0.2 m/s: {steps} steps, then {refusal}");
}

/// A component on a load line that returns no effort and keeps the matrix
/// it was handed: what the kernel says its neighbours answer.
struct Probe {
    ports: [PortSpec; 2],
}

impl Probe {
    fn new() -> Self {
        let end = |name: &str| {
            PortSpec::power(
                StableId::new(name).unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                core::num::NonZeroUsize::MIN,
            )
            .unwrap()
        };
        Self {
            ports: [end("a"), end("b")],
        }
    }
}

impl Component for Probe {
    type State = [f64; 4];
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _handed: &[f64; 4]) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &[f64; 4],
        dt_s: f64,
        inputs: &[f64],
        new: &mut [f64; 4],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 4], new, outputs, cx)
    }
}

impl LoadLine for Probe {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn propose_on_line(
        &self,
        _old: &[f64; 4],
        _dt_s: f64,
        _inputs: &[f64],
        compliance: &[f64],
        new: &mut [f64; 4],
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        new.copy_from_slice(compliance);
        outputs.fill(0.0);
        Ok(Proposal::Step(Receipt::lossless(0.0, f64::MIN_POSITIVE)))
    }
}

/// What joins the two masses of [`handed`] besides the probe.
enum Between {
    Nothing,
    Damper(f64),
    Gyrator(f64),
    /// No second mass: the probe's `b` is prescribed.
    GivenSide,
    /// Nothing between the masses; the first is on a spring to a wall.
    SpringToAWall(f64),
}

/// The matrix the kernel hands a probe whose `a` is on a mass of 2 kg and
/// whose `b` is on a mass of 0.5 kg, in a step of `dt`.
fn handed(between: &Between, dt: f64) -> [f64; 4] {
    let mechanical = PortKind::MechanicalForceVelocity;
    let mut builder = SystemBuilder::new();
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let one = builder.add_with_tangent("one", Mass::new(2.0).unwrap(), rest).unwrap();
    let probe = builder.add_on_line("probe", Probe::new(), [f64::NAN; 4]).unwrap();
    let mut at_one = vec![port(&builder, probe, "a")];
    let mut at_two = vec![port(&builder, probe, "b")];
    match between {
        Between::Nothing | Between::GivenSide => {}
        Between::SpringToAWall(stiffness) => {
            let spring = builder
                .add_with_tangent("spring", Spring::new(*stiffness, 0.0).unwrap(), 0.0)
                .unwrap();
            builder
                .boundary(port(&builder, spring, "b"), Box::new(|_, _, values| values[0] = 0.0))
                .unwrap();
            at_one.push(port(&builder, spring, "a"));
        }
        Between::Damper(damping) => {
            let damper = builder
                .add_with_tangent("damper", Damper::new(*damping).unwrap(), ())
                .unwrap();
            builder
                .boundary(port(&builder, damper, "heat"), Box::new(|_, _, values| values[0] = 300.0))
                .unwrap();
            at_one.push(port(&builder, damper, "a"));
            at_two.push(port(&builder, damper, "b"));
        }
        Between::Gyrator(ratio) => {
            let gyrator = builder
                .add_with_tangent("gyrator", Gyrator::new(mechanical, mechanical, *ratio).unwrap(), ())
                .unwrap();
            at_one.push(port(&builder, gyrator, "a"));
            at_two.push(port(&builder, gyrator, "b"));
        }
    }
    builder.common_flow(port(&builder, one, "p"), &at_one).unwrap();
    if matches!(between, Between::GivenSide) {
        builder.boundary(at_two[0], Box::new(|_, _, values| values[0] = 0.3)).unwrap();
    } else {
        let two = builder.add_with_tangent("two", Mass::new(0.5).unwrap(), rest).unwrap();
        builder.common_flow(port(&builder, two, "p"), &at_two).unwrap();
    }
    let mut system = builder.seal().unwrap();
    with_cx(|cx| system.step(dt, cx)).unwrap();
    *system.state::<[f64; 4]>(probe).unwrap()
}

#[test]
fn the_line_a_component_is_handed_is_what_its_neighbours_answer() {
    // Within a row of the component's relation the step a block accepts
    // does not depend on the line's matrix: the component settles where
    // the flow on the line is the flow its neighbours give, whatever
    // positive matrix it was told, and only the path of the iteration
    // differs. (Where the relation jumps, at a slider's break-away, the
    // matrix decides.) So the matrix is read directly, by a component
    // that keeps what it is handed and returns no effort.
    //
    // A mass answers a force with a mean velocity of dt / (2 m) per
    // newton. With g1 and g2 those of the two masses:
    //   nothing between them:      diag(g1, g2);
    //   a damper d:                [[g1 (1 + g2 d), g1 g2 d], [g1 g2 d, g2 (1 + g1 d)]] / (1 + (g1 + g2) d);
    //   a gyrator G (e_a = G f_b, e_b = -G f_a):
    //                              [[g1, -g1 g2 G], [g1 g2 G, g2]] / (1 + g1 g2 G^2), not symmetric;
    //   the second side prescribed: [[g1, 0], [0, 0]];
    //   the first mass on a spring k to a wall, which answers a velocity
    //   with k dt / 2 of force: diag(g1 / (1 + g1 k dt / 2), g2).
    let dt = 1.0e-3_f64;
    let (g1, g2) = (0.5 * dt / 2.0, 0.5 * dt / 0.5);
    let (d, ratio, stiffness) = (700.0_f64, 900.0_f64, 3.0e6_f64);
    let damped = 1.0 + (g1 + g2) * d;
    let turned = 1.0 + g1 * g2 * ratio * ratio;
    let cases = [
        (Between::Nothing, [g1, 0.0, 0.0, g2]),
        (
            Between::Damper(d),
            [
                g1 * (1.0 + g2 * d) / damped,
                g1 * g2 * d / damped,
                g1 * g2 * d / damped,
                g2 * (1.0 + g1 * d) / damped,
            ],
        ),
        (
            Between::Gyrator(ratio),
            [g1 / turned, -g1 * g2 * ratio / turned, g1 * g2 * ratio / turned, g2 / turned],
        ),
        (Between::GivenSide, [g1, 0.0, 0.0, 0.0]),
        (
            Between::SpringToAWall(stiffness),
            [g1 / (1.0 + g1 * stiffness * 0.5 * dt), 0.0, 0.0, g2],
        ),
    ];
    for (case, (between, expected)) in cases.iter().enumerate() {
        let got = handed(between, dt);
        for (entry, (got, want)) in got.iter().zip(expected).enumerate() {
            assert!(
                (got - want).abs() <= 1.0e-12 * g2,
                "case {case} entry {entry}: {got} against {want} ({:?})",
                handed(between, dt)
            );
        }
    }
    // The damper's and the gyrator's off-diagonal entries are not small:
    // the cases tell an entry from its transpose and from zero.
    assert!(g1 * g2 * d / damped > 0.1 * g1 && g1 * g2 * ratio / turned > 0.1 * g1);
    // And the spring takes more than a fifth off the mass's own answer.
    assert!(g1 * stiffness * 0.5 * dt > 0.25);
}

/// The slider, with a record of the step and the two flows it is handed at
/// every call on its line; while `cut` holds a length, it asks to cut any
/// longer step to it.
struct Recorded {
    slider: Slider,
    calls: std::sync::Arc<std::sync::Mutex<Vec<[f64; 3]>>>,
    cut: std::sync::Arc<std::sync::Mutex<Option<f64>>>,
}

impl Component for Recorded {
    type State = Slip;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        self.slider.ports()
    }

    fn class(&self) -> Class {
        self.slider.class()
    }

    fn stored_energy_j(&self, slip: &Slip) -> Result<f64, ElementError> {
        self.slider.stored_energy_j(slip)
    }

    fn propose(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.slider.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl LoadLine for Recorded {
    fn line_ports(&self) -> &[usize] {
        self.slider.line_ports()
    }

    fn propose_on_line(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.calls
            .lock()
            .unwrap()
            .push([dt_s, inputs[0], inputs[1]]);
        if let Some(cut_s) = *self.cut.lock().unwrap()
            && dt_s > cut_s
        {
            return Ok(Proposal::Shorten { dt_s: cut_s });
        }
        self.slider
            .propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn a_component_on_a_line_is_first_asked_at_the_released_flows_of_every_trial() {
    // 2 kg thrown at 1 m/s and pushed with 3 N rubs through the slider on
    // 0.5 kg at rest, and the slider's loss warms a body. The two slide for
    // about 76 steps of 1 ms, then go on together, held. Three steps are
    // cut to half by the slider's asking, two while they slide and one
    // while they are held. Whatever it answers, the slider is first asked
    // in every trial at the released flows of that trial: what each mass
    // does over the trial's dt with no force from the slider, its velocity
    // plus the push's share for the first, and its velocity for the second,
    // each port receiving minus its mass's. The masses are affine, so the
    // update before the first ask gives them to rounding. (Asked at the
    // last step's crossing flows plus the line times the slider's last
    // forces, as before that update, the first ask of the run is off by the
    // 1 m/s of the throw, and later ones by what a step of the slider's
    // force changes.)
    let (dt, push_n) = (1.0e-3_f64, 3.0_f64);
    let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let cut = std::sync::Arc::new(std::sync::Mutex::new(None));
    let mut builder = SystemBuilder::new();
    let moving = |velocity_m_s: f64| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let one = builder
        .add_with_tangent("one", Mass::new(2.0).unwrap(), moving(1.0))
        .unwrap();
    let two = builder
        .add_with_tangent("two", Mass::new(0.5).unwrap(), moving(0.0))
        .unwrap();
    let recorded = Recorded {
        slider: slider(),
        calls: std::sync::Arc::clone(&calls),
        cut: std::sync::Arc::clone(&cut),
    };
    // The slider's `a` receives minus the first mass's velocity.
    let rubbing = builder
        .add_on_line("slider", recorded, Slip::sliding(-1.0))
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(push_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let at_one = [port(&builder, rubbing, "a"), port(&builder, push, "p")];
    let at_two = [port(&builder, rubbing, "b")];
    builder
        .common_flow(port(&builder, one, "p"), &at_one)
        .unwrap();
    builder
        .common_flow(port(&builder, two, "p"), &at_two)
        .unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder
        .common_effort(port(&builder, body, "heat"), &[warmth])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));

    let cut_steps = [10_u32, 50, 150];
    let (mut trials, mut cuts, mut sliding, mut held) = (0_u32, 0_u32, 0_u32, 0_u32);
    // Over every step and over the steps after the first: the largest
    // offset of a first ask from the released flows [m/s], its step, its
    // trial's dt, and what it was handed against the released flows.
    let mut worst = [(0.0_f64, 0_u32, 0.0_f64, [0.0_f64; 2], [0.0_f64; 2]); 2];
    with_cx(|cx| {
        for step in 1..=200_u32 {
            let first = system.state::<Motion>(one).unwrap().velocity_m_s;
            let second = system.state::<Motion>(two).unwrap().velocity_m_s;
            *cut.lock().unwrap() = cut_steps.contains(&step).then_some(0.5 * dt);
            let report = system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
            *cut.lock().unwrap() = None;
            cuts += u32::from(report.dt_s < dt);
            if system.state::<Slip>(rubbing).unwrap().is_stuck() {
                held += 1;
            } else {
                sliding += 1;
            }
            let asked: Vec<[f64; 3]> = core::mem::take(&mut *calls.lock().unwrap());
            for (k, call) in asked.iter().enumerate() {
                // A trial is a run of asks for one dt; its first is checked.
                if k > 0 && asked[k - 1][0] == call[0] {
                    continue;
                }
                trials += 1;
                let released = [-(first + push_n * call[0] / (2.0 * 2.0)), -second];
                let handed = [call[1], call[2]];
                let offset = handed
                    .iter()
                    .zip(&released)
                    .map(|(handed, released)| (handed - released).abs())
                    .fold(0.0_f64, f64::max);
                let spans = if step > 1 { 2 } else { 1 };
                for slot in &mut worst[..spans] {
                    if offset > slot.0 {
                        *slot = (offset, step, call[0], handed, released);
                    }
                }
            }
        }
    });
    let shown = |(offset, step, dt_s, handed, released): (f64, u32, f64, [f64; 2], [f64; 2])| {
        format!(
            "{offset:.3e} m/s (step {step}, a trial of {dt_s:e} s: handed {handed:?}, released \
             {released:?})"
        )
    };
    eprintln!(
        "first asks of {trials} trials ({cuts} cut steps; {sliding} steps end sliding, {held} \
         held): off the released flows by {} over every step, by {} over the steps after the first",
        shown(worst[0]),
        shown(worst[1])
    );
    assert_eq!(cuts, 3);
    assert_eq!(trials, 200 + 3);
    assert!(sliding > 20 && held > 20, "{sliding} sliding, {held} held");
    assert!(
        worst[0].0 <= 1.0e-13,
        "the first ask of a trial is off the released flows by {} over every step, by {} over \
         the steps after the first",
        shown(worst[0]),
        shown(worst[1])
    );
}

/// Two masses that rub on each other through the slider and are joined by
/// `between` besides: 2 kg thrown at 1 m/s at the slider's `a`, 0.5 kg at
/// rest at its `b`, the slider's loss collected in a body at 300 K, and
/// with `push_n` a constant force on the first mass. Returns the system,
/// the two masses, the slider and the body.
fn rubbing_beside(
    between: &Between,
    push_n: Option<f64>,
) -> (System, NodeId, NodeId, NodeId, NodeId) {
    let mechanical = PortKind::MechanicalForceVelocity;
    let mut builder = SystemBuilder::new();
    let moving = |velocity_m_s: f64| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let one = builder.add_with_tangent("one", Mass::new(2.0).unwrap(), moving(1.0)).unwrap();
    let two = builder.add_with_tangent("two", Mass::new(0.5).unwrap(), moving(0.0)).unwrap();
    // The slider's `a` receives minus the first mass's velocity.
    let rubbing = builder.add_on_line("slider", slider(), Slip::sliding(-1.0)).unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    let mut at_one = vec![port(&builder, rubbing, "a")];
    let mut at_two = vec![port(&builder, rubbing, "b")];
    match between {
        Between::Damper(damping) => {
            let damper = builder
                .add_with_tangent("damper", Damper::new(*damping).unwrap(), ())
                .unwrap();
            builder
                .boundary(port(&builder, damper, "heat"), Box::new(|_, _, values| values[0] = 300.0))
                .unwrap();
            at_one.push(port(&builder, damper, "a"));
            at_two.push(port(&builder, damper, "b"));
        }
        Between::Gyrator(ratio) => {
            let gyrator = builder
                .add_with_tangent("gyrator", Gyrator::new(mechanical, mechanical, *ratio).unwrap(), ())
                .unwrap();
            at_one.push(port(&builder, gyrator, "a"));
            at_two.push(port(&builder, gyrator, "b"));
        }
        Between::Nothing | Between::GivenSide | Between::SpringToAWall(_) => {
            panic!("something joins the masses here")
        }
    }
    if let Some(push_n) = push_n {
        let push = builder
            .add_with_tangent(
                "push",
                Force::new(Waveform::constant(push_n).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        at_one.push(port(&builder, push, "p"));
    }
    builder.common_flow(port(&builder, one, "p"), &at_one).unwrap();
    builder.common_flow(port(&builder, two, "p"), &at_two).unwrap();
    let warmth = port(&builder, rubbing, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[warmth]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(warmth, &[300.0]));
    (system, one, two, rubbing, body)
}

#[test]
fn a_slider_gives_its_forces_through_a_line_whose_ends_answer_each_other() {
    // The probe above reads the matrix and returns no effort, so nothing
    // there is taken away through an entry off the diagonal. Here the real
    // slider rubs between the two masses, beside the damper and beside
    // the gyrator of that test: every flow on its line is the released
    // flow less what both of its forces take away, one of them through an
    // entry off the diagonal, and with the gyrator through a matrix that
    // is not symmetric. A matrix used transposed, or an entry dropped,
    // publishes flows that are not the masses' and the junctions' audit
    // refuses the step.
    let kinetic_j = |system: &System, one: NodeId, two: NodeId| {
        let (first, second) = (
            system.state::<Motion>(one).unwrap().velocity_m_s,
            system.state::<Motion>(two).unwrap().velocity_m_s,
        );
        (0.5 * 2.0 * first * first + 0.5 * 0.5 * second * second, 2.0 * first + 0.5 * second)
    };
    let slide_n = SLIDES * NORMAL_N;

    // Beside a damper of 700 N s/m, with 30 N pushing the first mass. The
    // momentum is the 2 kg m/s they start with and the push's impulse.
    // The damper takes their difference of 1 m/s away within a few steps.
    // To go along at 12 m/s2 the second mass needs 6 N, more than the
    // 5.88 N the slider gives while it slides, so the slider never sticks:
    // the two slide on at the difference at which the damper gives the
    // missing 0.12 N, 0.17 mm/s, with the sliding force at every step. The
    // push's work is what the two gain and what the slider and the damper
    // lose.
    let (push_n, dt) = (30.0_f64, 1.0e-3_f64);
    let (mut system, one, two, rubbing, body) = rubbing_beside(&Between::Damper(700.0), Some(push_n));
    let at_one = system.port(rubbing, "a").unwrap();
    let mut rubbed = 0_u32;
    with_cx(|cx| {
        for step in 0..80_u32 {
            system.step(dt, cx).unwrap_or_else(|error| panic!("damper, step {step}: {error}"));
            let (_, momentum) = kinetic_j(&system, one, two);
            let given = 2.0 + push_n * dt * f64::from(step + 1);
            assert!((momentum - given).abs() <= 1.0e-12 * given, "step {step}: {momentum}");
            let force = system.last_output(at_one)[0];
            rubbed += u32::from(force.abs() == slide_n);
            if step >= 20 {
                let (first, second) = (
                    system.state::<Motion>(one).unwrap().velocity_m_s,
                    system.state::<Motion>(two).unwrap().velocity_m_s,
                );
                let missing_n = 0.5 * push_n / 2.5 - slide_n;
                assert!(!system.state::<Slip>(rubbing).unwrap().is_stuck(), "step {step}");
                assert_eq!(force.abs(), slide_n, "step {step}");
                assert!(
                    (first - second - missing_n / 700.0).abs() <= 1.0e-12,
                    "step {step}: {first} {second}"
                );
            }
        }
    });
    let (left_j, _) = kinetic_j(&system, one, two);
    let work_j = push_n * system.state::<Motion>(one).unwrap().position_m;
    let rubbed_j = received_j(&system, body);
    let damped_j = -system
        .ledger()
        .boundaries
        .iter()
        .find(|account| account.port == "damper.heat")
        .expect("the damper's heat")
        .energy_in_j
        .value();
    eprintln!(
        "a slider beside a damper: {rubbed} of 80 steps at the sliding force; the push's \
         {work_j:e} J are {:e} J gained, {rubbed_j:e} J lost in the slider and {damped_j:e} J \
         in the damper",
        left_j - 1.0
    );
    assert!((left_j - 1.0 + rubbed_j + damped_j - work_j).abs() <= 1.0e-12 * work_j);
    assert!(rubbed >= 60 && rubbed_j > 1.0e-3 && damped_j > 1.0e-3, "{rubbed}");

    // Beside a gyrator of 900 N s/m (force on each mass: 900 N s/m times
    // the other's velocity, with opposite signs), which loses nothing:
    // what the two masses lose is the slider's heat, at every step. They
    // swing against each other 143 times a second and rub until they move
    // as one; then the gyrator pushes them apart with 900 N s/m times
    // their velocity, and the slider, stuck, holds exactly that.
    let (ratio, dt) = (900.0_f64, 1.0e-4_f64);
    let (mut system, one, two, rubbing, body) = rubbing_beside(&Between::Gyrator(ratio), None);
    let at_one = system.port(rubbing, "a").unwrap();
    let (mut rubbed, mut stuck_steps, mut heat_before) = (0_u32, 0_u32, 0.0_f64);
    with_cx(|cx| {
        for step in 0..6000 {
            system.step(dt, cx).unwrap_or_else(|error| panic!("gyrator, step {step}: {error}"));
            let (left_j, _) = kinetic_j(&system, one, two);
            let heat_j = received_j(&system, body);
            assert!((left_j + heat_j - 1.0).abs() <= 1.0e-12, "step {step}: {left_j} {heat_j}");
            assert!(heat_j >= heat_before, "step {step}");
            heat_before = heat_j;
            rubbed += u32::from(system.last_output(at_one)[0].abs() == slide_n);
            stuck_steps += u32::from(system.state::<Slip>(rubbing).unwrap().is_stuck());
        }
    });
    let (left_j, _) = kinetic_j(&system, one, two);
    let ends = (
        system.state::<Motion>(one).unwrap().velocity_m_s,
        system.state::<Motion>(two).unwrap().velocity_m_s,
    );
    let force = system.last_output(at_one)[0];
    eprintln!(
        "a slider beside a gyrator: {rubbed} steps at the sliding force, {stuck_steps} stuck; \
         {left_j:e} J left of 1 J, both at {:e} m/s, held with {force:e} N",
        ends.0
    );
    assert!(rubbed >= 1000 && stuck_steps >= 1000, "{rubbed} {stuck_steps}");
    assert!(system.state::<Slip>(rubbing).unwrap().is_stuck());
    assert!((ends.0 - ends.1).abs() <= 1.0e-12, "{ends:?}");
    assert!((force.abs() - ratio * ends.0.abs()).abs() <= 1.0e-9, "{force} {ends:?}");
    assert!(force.abs() > 0.5 * slide_n && force.abs() <= HOLDS * NORMAL_N);
}

#[test]
fn a_refused_step_leaves_nothing_of_its_line_behind() {
    // A mass thrown along the ground on a spring that is added without its
    // tangent, so the kernel differences it. One system is asked for a
    // step with no iteration allowed: the line has been found, the drivers
    // carry its matrix and the inputs hold released flows when the solve
    // gives up. It is also asked for a step with one iteration allowed,
    // which gives up later, after an assembly and an update. The step
    // that follows each of the six refusals must be the step of a system
    // that was never refused, to the bit.
    let build = || {
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent(
                "mass",
                Mass::new(MASS_KG).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 1.3,
                },
            )
            .unwrap();
        let rubbing = builder.add_on_line("slider", slider(), Slip::sliding(-1.3)).unwrap();
        let spring = builder.add("spring", Spring::new(400.0, 0.0).unwrap(), 0.0).unwrap();
        builder
            .common_flow(
                port(&builder, mass, "p"),
                &[port(&builder, rubbing, "a"), port(&builder, spring, "a")],
            )
            .unwrap();
        for held in [port(&builder, rubbing, "b"), port(&builder, spring, "b")] {
            builder.boundary(held, Box::new(|_, _, values| values[0] = 0.0)).unwrap();
        }
        builder
            .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
        (builder.seal().unwrap(), mass, rubbing, spring)
    };
    let (mut tried, mass, rubbing, spring) = build();
    let (mut fresh, _, _, _) = build();
    let dt = 5.0e-3;
    let mut refusals = 0_u32;
    with_cx(|cx| {
        for step in 0..60 {
            if step % 20 == 0 {
                tried.set_config(SolveConfig {
                    max_iterations: 0,
                    ..SolveConfig::default()
                });
                let refusal = tried.step(dt, cx).expect_err("no iteration allowed");
                assert!(matches!(refusal, SystemError::NotConverged { .. }), "{refusal}");
                assert_eq!(tried.steps(), step);
                tried.set_config(SolveConfig::default());
                refusals += 1;
            }
            if step % 20 == 10 {
                // One iteration allowed where two are needed: the solve
                // gives up after an assembly and an update, with the
                // slider differenced and proposed at a perturbed flow.
                tried.set_config(SolveConfig {
                    max_iterations: 1,
                    ..SolveConfig::default()
                });
                let refusal = tried.step(dt, cx).expect_err("one iteration is not enough");
                assert!(
                    matches!(refusal, SystemError::NotConverged { iterations: 1, .. }),
                    "{refusal}"
                );
                assert_eq!(tried.steps(), step);
                tried.set_config(SolveConfig::default());
                refusals += 1;
            }
            let (a, b) = (tried.step(dt, cx).unwrap(), fresh.step(dt, cx).unwrap());
            assert_eq!(a.iterations, b.iterations, "step {step}");
            let (ours, theirs) = (
                *tried.state::<Motion>(mass).unwrap(),
                *fresh.state::<Motion>(mass).unwrap(),
            );
            assert_eq!(ours.velocity_m_s.to_bits(), theirs.velocity_m_s.to_bits(), "step {step}");
            assert_eq!(ours.position_m.to_bits(), theirs.position_m.to_bits(), "step {step}");
            assert_eq!(
                tried.state::<Slip>(rubbing).unwrap(),
                fresh.state::<Slip>(rubbing).unwrap(),
                "step {step}"
            );
            assert_eq!(
                tried.state::<f64>(spring).unwrap().to_bits(),
                fresh.state::<f64>(spring).unwrap().to_bits(),
                "step {step}"
            );
        }
    });
    assert_eq!(refusals, 6);
    // The run was not idle: the mass turned back on the spring.
    assert!(fresh.state::<Motion>(mass).unwrap().velocity_m_s < 0.0);
}

/// A body of negative inertia: pushed, it moves against the push. Nothing
/// physical, and the neighbour that makes a load line unusable.
struct Antimass {
    mass_kg: f64,
    ports: [PortSpec; 1],
}

impl Antimass {
    fn new(mass_kg: f64) -> Self {
        let id = StableId::new("p").unwrap();
        let width = core::num::NonZeroUsize::MIN;
        Self {
            mass_kg,
            ports: [PortSpec::power(
                id,
                PortKind::MechanicalForceVelocity,
                PortRole::EffortIn,
                width,
            )
            .unwrap()],
        }
    }
}

impl Component for Antimass {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, velocity: &f64) -> Result<f64, ElementError> {
        Ok(0.5 * self.mass_kg * velocity * velocity)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        *new = old + dt_s * inputs[0] / self.mass_kg;
        outputs[0] = 0.5 * (old + *new);
        let stored_change_j = self.mass_kg * outputs[0] * (*new - old);
        Ok(Proposal::Step(Receipt::lossless(
            stored_change_j,
            64.0 * f64::EPSILON * stored_change_j.abs() + f64::MIN_POSITIVE,
        )))
    }
}

/// A slider that names these ports as its line.
struct Misplaced(Slider, Vec<usize>);

impl Component for Misplaced {
    type State = Slip;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        self.0.ports()
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, slip: &Slip) -> Result<f64, ElementError> {
        self.0.stored_energy_j(slip)
    }

    fn propose(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.0.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl LoadLine for Misplaced {
    fn line_ports(&self) -> &[usize] {
        &self.1
    }

    fn propose_on_line(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.0.propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

#[test]
fn what_the_tier_cannot_solve_is_refused_by_name() {
    // A line that is not over power ports that receive a flow, in
    // increasing order, is refused when the component is added, and claims
    // no name: on the heat port, over no port, out of order, a port twice,
    // a port the component does not have.
    let mut builder = SystemBuilder::new();
    for line in [vec![0, 2], vec![], vec![1, 0], vec![0, 0], vec![0, 5]] {
        let refusal = builder
            .add_on_line("misplaced", Misplaced(slider(), line.clone()), Slip::default())
            .unwrap_err();
        assert_eq!(
            refusal,
            AdmissionError::LinePorts {
                node: "misplaced".to_owned()
            },
            "{line:?}"
        );
    }
    builder.add_on_line("misplaced", slider(), Slip::default()).unwrap();

    // Two sliders on one mass are two lines in one loop: refused at
    // sealing, both named.
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let mut members = Vec::new();
    for name in ["front", "rear"] {
        let rubbing = builder.add_on_line(name, slider(), Slip::default()).unwrap();
        builder
            .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        builder
            .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
        members.push(port(&builder, rubbing, "a"));
    }
    builder.common_flow(port(&builder, mass, "p"), &members).unwrap();
    let Err(refusal) = builder.seal() else {
        panic!("two lines in one loop are refused");
    };
    assert_eq!(
        refusal,
        AdmissionError::TwoLines {
            first: "front".to_owned(),
            second: "rear".to_owned()
        }
    );

    // A neighbour that answers a force with motion against it gives a
    // line of negative compliance, and one of no inertia answers with a
    // motion that is not finite, so that there is no line at all: the
    // first step is refused, the slider named with the reason, and nothing
    // is published.
    for (inertia_kg, reason) in [(-MASS_KG, "symmetric part is negative"), (0.0, "do not determine")] {
        let mut builder = SystemBuilder::new();
        let wrong = builder.add("wrong", Antimass::new(inertia_kg), 0.0).unwrap();
        let rubbing = builder.add_on_line("slider", slider(), Slip::default()).unwrap();
        builder
            .common_flow(port(&builder, wrong, "p"), &[port(&builder, rubbing, "a")])
            .unwrap();
        builder
            .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        builder
            .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
            .unwrap();
        let mut system = builder.seal().unwrap();
        let refusal = with_cx(|cx| system.step(1.0e-3, cx)).unwrap_err();
        assert!(
            matches!(&refusal, SystemError::LoadLine { node, reason: why }
                if node == "slider" && why.contains(reason)),
            "{inertia_kg} kg: {refusal}"
        );
        assert_eq!(system.steps(), 0);
    }

    // A slider whose two ends are both prescribed is in no loop: there is
    // no line to find and it is stepped by its mandatory tier, sliding at
    // its kinetic force.
    let mut builder = SystemBuilder::new();
    let rubbing = builder
        .add_on_line("slider", slider(), Slip::stuck())
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "a"), Box::new(|_, _, values| values[0] = 0.2))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "b"), Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
    builder
        .boundary(port(&builder, rubbing, "heat"), Box::new(|_, _, values| values[0] = 300.0))
        .unwrap();
    let mut system = builder.seal().unwrap();
    let report = with_cx(|cx| system.step(1.0e-3, cx)).unwrap();
    assert_eq!(report.iterations, 0);
    assert_eq!(
        system.last_output(system.port(rubbing, "a").unwrap()),
        [SLIDES * NORMAL_N]
    );
}


/// Dry friction to the ground with one constant force, which does not
/// arrest inside a step: where its slip would pass zero it asks for the
/// step that ends there, as a product's friction component does (the event
/// path). Port `a` receives a velocity and returns the force; the loss
/// leaves through `heat`. Its state is the slip its last step ended with,
/// zero at rest.
struct CutsAtAStop {
    force_n: f64,
    ports: [PortSpec; 2],
}

/// How near the sliding force the force that ends a step at rest has to
/// be [N] for [`CutsAtAStop`] to end the step there.
const STOP_BAND_N: f64 = 1.0e-9;

impl CutsAtAStop {
    fn new(force_n: f64) -> Self {
        Self {
            force_n,
            ports: [
                PortSpec::power(
                    StableId::new("a").unwrap(),
                    PortKind::MechanicalForceVelocity,
                    PortRole::FlowIn,
                    core::num::NonZeroUsize::MIN,
                )
                .unwrap(),
                PortSpec::heat(StableId::new("heat").unwrap(), PortRole::EffortIn).unwrap(),
            ],
        }
    }
}

impl Component for CutsAtAStop {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _slip: &f64) -> Result<f64, ElementError> {
        Ok(0.0)
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
        self.propose_on_line(old, dt_s, inputs, &[0.0], new, outputs, cx)
    }
}

impl LoadLine for CutsAtAStop {
    fn line_ports(&self) -> &[usize] {
        &[0]
    }

    fn propose_on_line(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let (start, released, c, hold) = (*old, inputs[0], compliance[0], self.force_n);
        let slides = |way: f64| {
            let force = hold * way;
            (force, released - c * force)
        };
        let (force, mean, ends) = if !(c > 0.0) {
            (if released == 0.0 { 0.0 } else { hold * released.signum() }, released, released)
        } else if start == 0.0 {
            if released.abs() <= c * hold {
                ((released / c).clamp(-hold, hold), 0.0, 0.0)
            } else {
                let (force, mean) = slides(released.signum());
                (force, mean, 2.0 * mean)
            }
        } else {
            // `along` is the force that ends the step at zero slip.
            let way = start.signum();
            let along = way * (released - 0.5 * start) / c;
            if along > hold + STOP_BAND_N {
                let (force, mean) = slides(way);
                (force, mean, 2.0 * mean - start)
            } else {
                if along < hold - STOP_BAND_N {
                    // The slip passes zero inside the step: where, by the
                    // line of a body under a held force.
                    let shorter = dt_s * (0.5 * start.abs()) / (c * hold - way * (released - start));
                    if shorter > 0.0 && shorter < dt_s {
                        return Ok(Proposal::Shorten { dt_s: shorter });
                    }
                }
                let force = way * along.max(0.0);
                (force, released - c * force, 0.0)
            }
        };
        *new = ends;
        let loss_w = force * mean;
        outputs[0] = force;
        outputs[1] = -loss_w;
        let dissipated_j = dt_s * loss_w.max(0.0);
        let scale = dt_s * force.abs() * (released.abs() + force.abs() * c.abs()) + dissipated_j;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dissipated_j / inputs[1],
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

/// A mass of 2 kg thrown at 1 m/s along a surface that brakes it with 3 N
/// through [`CutsAtAStop`]: it stops after 2/3 s and 1/3 m. The brake's
/// heat warms a body whose temperature it is handed back, so the block
/// needs a second update in every step. `told`: whether the brake's port is
/// given the flow it starts with.
fn thrown_onto_a_brake(told: bool) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 1.0,
            },
        )
        .unwrap();
    // The port receives minus the mass's velocity.
    let brake = builder.add_on_line("brake", CutsAtAStop::new(3.0), -1.0).unwrap();
    builder
        .common_flow(port(&builder, mass, "p"), &[port(&builder, brake, "a")])
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .heat_route(
            port(&builder, brake, "heat"),
            &[(port(&builder, body, "heat"), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(system.port(brake, "heat").unwrap(), &[300.0]));
    if told {
        let at = system.port(brake, "a").unwrap();
        assert!(system.guess(at, &[-1.0]));
    }
    (system, mass, brake)
}

#[test]
fn a_component_that_cuts_at_a_stop_is_first_asked_at_the_released_flows_of_the_step() {
    // A block's solve starts from the flows of the last step, and before
    // the first step from none. A component that enters a step with a slip
    // of 1 m/s its port shows nothing of sees, at those flows, a stop
    // inside the step, and asks for a shorter one; the shorter trial starts
    // from the same flows and asks again. Before the kernel released the
    // line first, the first step of this mass was refused after 64
    // shortenings (`the step was shortened too often`). Now the other
    // components take one update with the brake's force held, which gives
    // the released flow of the step, and the brake is asked again there:
    // the step is taken whole.
    let dt = 0.05_f64;
    let (mut system, mass, brake) = thrown_onto_a_brake(false);
    let report = with_cx(|cx| system.step(dt, cx)).unwrap();
    assert_eq!((report.shortenings, report.dt_s), (0, dt));
    let after = system.state::<Motion>(mass).unwrap().velocity_m_s;
    assert!((after - (1.0 - 1.5 * dt)).abs() <= 1.0e-15, "{after}");
    assert_eq!(*system.state::<f64>(brake).unwrap(), -after);
    assert_eq!(system.last_output(system.port(brake, "a").unwrap()), [-3.0]);
}

#[test]
fn a_perturbed_proposal_that_asks_to_shorten_does_not_shorten_the_step() {
    // The mass stops at 2/3 s, between two steps of 0.05 s: the step from
    // 0.65 s is cut once, to end there, and ends at rest. At the iterate of
    // that cut step the kernel differences the brake at a flow moved by
    // sqrt(eps), and the moved flow lies beyond the stop, where the brake
    // asks for a shorter step. That proposal is the kernel's own, not the
    // step's. Taken as a request, it cut the step a second time, a hair
    // short of the stop (slip -1.3e-8 m/s in the first product run), and
    // the step after that was refused. The column is now differenced on
    // the other side of the iterate.
    let dt = 0.05_f64;
    for told in [true, false] {
        let (mut system, mass, brake) = thrown_onto_a_brake(told);
        let mut cuts = 0;
        let mut stopped_at = None;
        with_cx(|cx| {
            for step in 0..20 {
                let was_moving = *system.state::<f64>(brake).unwrap() != 0.0;
                let report =
                    system.step(dt, cx).unwrap_or_else(|error| panic!("step {step}: {error}"));
                cuts += report.shortenings;
                let slip = *system.state::<f64>(brake).unwrap();
                if was_moving && slip == 0.0 {
                    // The one cut step, and it ends at the stop.
                    assert_eq!(report.shortenings, 1, "step {step}");
                    assert!(stopped_at.is_none());
                    stopped_at = Some(report.time_s);
                } else {
                    assert_eq!((report.shortenings, report.dt_s), (0, dt), "step {step}");
                    assert_eq!(was_moving, slip != 0.0, "step {step}");
                }
            }
        });
        let stopped_at = stopped_at.expect("the mass stops");
        assert!((stopped_at - 2.0 / 3.0).abs() <= 1.0e-14, "{stopped_at}");
        assert_eq!(cuts, 1);
        // At rest, a third of a metre on, to the rounding of the run.
        let ends = *system.state::<Motion>(mass).unwrap();
        assert!(ends.velocity_m_s.abs() <= 1.0e-14, "{}", ends.velocity_m_s);
        assert!((ends.position_m - 1.0 / 3.0).abs() <= 1.0e-14, "{}", ends.position_m);
        if told {
            println!(
                "a mass thrown onto a brake that cuts at the stop: stopped at {stopped_at:e} s, \
                 {cuts} shortening, velocity {:e} m/s, position {:e} m",
                ends.velocity_m_s, ends.position_m
            );
        }
    }
}

/// A damper of 1 N s/m between two junctions, differenced by the kernel,
/// with nowhere for its heat to go: the loss is booked as lost. With a step
/// above an eighth of a second it asks for half the step where the flow at
/// `b` is within a micrometre per second of zero and not zero: above zero,
/// or (`either_side`) on either side. No iterate of a solve lies there,
/// only the kernel's perturbed proposals about a flow of zero.
struct Touchy {
    either_side: bool,
    ports: [PortSpec; 2],
}

impl Component for Touchy {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let beside_zero = inputs[1] != 0.0 && inputs[1].abs() < 1.0e-6;
        if dt_s > 0.125 && beside_zero && (self.either_side || inputs[1] > 0.0) {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }
        let force = inputs[0] - inputs[1];
        outputs[0] = force;
        outputs[1] = -force;
        let lost_j = dt_s * force * force;
        Ok(Proposal::Step(Receipt {
            numerical_loss_j: lost_j,
            ..Receipt::lossless(0.0, 64.0 * f64::EPSILON * lost_j + f64::MIN_POSITIVE)
        }))
    }
}

#[test]
fn a_differenced_component_shortens_a_step_only_if_it_asks_on_both_sides_of_the_iterate() {
    // No load line here. A mass pushed with 1 N drags a second one through
    // a damper that is differenced. Both start at rest, so the damper is
    // first differenced with no flow at its `b`, and the kernel's perturbed
    // proposals move that flow off zero by 1.5e-8 m/s, either way. That is
    // where the damper asks to shorten. No iterate lies there: the flow is
    // zero before the update and four hundredths of a metre per second
    // after it.
    //
    // Where the damper asks above zero only, the column is differenced
    // below it, and the step of one second is taken whole in one update:
    // the column is the damper's slope, with its sign. (The step was
    // halved three times, down to an eighth of a second.) Where it asks on
    // both sides it cannot be differenced at this step length at all, and
    // the step is halved until it can.
    //
    // The two masses of 2 kg from rest on the midpoint rule: their mean
    // velocities add up to dt / 4 and differ by 1 / (4 / dt + 2).
    for (either_side, shortenings, taken_s) in [(false, 0, 1.0_f64), (true, 3, 0.125)] {
        let at_rest = Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        };
        let end = |name: &str| {
            PortSpec::power(
                StableId::new(name).unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                core::num::NonZeroUsize::MIN,
            )
            .unwrap()
        };
        let mut builder = SystemBuilder::new();
        let pushed = builder.add_with_tangent("pushed", Mass::new(MASS_KG).unwrap(), at_rest).unwrap();
        let dragged = builder.add_with_tangent("dragged", Mass::new(MASS_KG).unwrap(), at_rest).unwrap();
        let pusher = builder
            .add_with_tangent(
                "push",
                Force::new(Waveform::constant(1.0).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        let touchy = builder
            .add(
                "touchy",
                Touchy {
                    either_side,
                    ports: [end("a"), end("b")],
                },
                (),
            )
            .unwrap();
        builder
            .common_flow(
                port(&builder, pushed, "p"),
                &[port(&builder, pusher, "p"), port(&builder, touchy, "a")],
            )
            .unwrap();
        builder
            .common_flow(port(&builder, dragged, "p"), &[port(&builder, touchy, "b")])
            .unwrap();
        let mut system = builder.seal().unwrap();
        let report = with_cx(|cx| system.step(1.0, cx)).unwrap();
        assert_eq!(
            (report.shortenings, report.dt_s, report.iterations),
            (shortenings, taken_s, 1),
            "{either_side}"
        );
        let (sum, difference) = (0.25 * taken_s, 1.0 / (4.0 / taken_s + 2.0));
        let first = system.state::<Motion>(pushed).unwrap().velocity_m_s;
        let second = system.state::<Motion>(dragged).unwrap().velocity_m_s;
        assert!((first - (sum + difference)).abs() <= 1.0e-15, "{either_side}: {first}");
        assert!((second - (sum - difference)).abs() <= 1.0e-15, "{either_side}: {second}");
    }
}

#[test]
fn a_push_that_creeps_up_to_the_static_limit_is_held_either_way() {
    // A stated limit of differencing a relation that jumps (the CONTRACT
    // states it). The slider holds 9.8 N and slides against 5.88 N. The
    // kernel differences it at a released flow moved UP by sqrt(eps) of
    // it, which is worth sqrt(eps) of the force: 1.46e-7 N near the limit.
    // Within that of the limit, with the held force on the side the move
    // goes to, the moved proposal has broken away: the column is the
    // difference across the jump from the held force to the sliding one,
    // not the held row's slope, and an update made with it leaves the
    // force wire with far more than rounding.
    //
    // A mass of 2 kg at rest under a push that rises to 2e-7 N inside the
    // limit in 0.1 s and creeps on to 5e-8 N inside it in 150 steps of
    // 1 ms, one way and the other. Every input of the slider is released
    // before its first ask, so it is first asked at the flow the mass
    // gives it, and this scene no longer reaches the limit: the mass is
    // held to rounding both ways, inside the window and out, to the end.
    // No step is refused. The most updates a held step takes inside the
    // window is printed.
    let dt = 1.0e-3_f64;
    let hold = HOLDS * NORMAL_N;
    let window = f64::EPSILON.sqrt() * hold;
    for way in [1.0_f64, -1.0] {
        let knots = vec![
            (0.0, 0.0),
            (0.1, way * (hold - 2.0e-7)),
            (0.25, way * (hold - 5.0e-8)),
        ];
        let wave = Waveform::table(knots.clone()).unwrap();
        let (mut system, mass, rubbing) = pushed_by(Waveform::table(knots).unwrap(), None);
        let at = system.port(rubbing, "a").unwrap();
        // Outside and inside the window: the most the force is off the
        // push, and the mass's largest speed, while it is held.
        let (mut off, mut speed) = ([0.0_f64; 2], [0.0_f64; 2]);
        let (mut let_go, mut updates) = (None, 0_u32);
        with_cx(|cx| {
            for step in 0..250 {
                let t0 = system.time_s();
                let report = system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("way {way}, step {step}: {error}"));
                let push = wave.mean(t0, dt);
                if !system.state::<Slip>(rubbing).unwrap().is_stuck() {
                    let_go = Some(hold - push.abs());
                    break;
                }
                let inside = usize::from(hold - push.abs() < window);
                off[inside] = off[inside].max((system.last_output(at)[0].abs() - push.abs()).abs());
                speed[inside] =
                    speed[inside].max(system.state::<Motion>(mass).unwrap().velocity_m_s.abs());
                if inside == 1 {
                    updates = updates.max(report.iterations);
                }
            }
        });
        println!(
            "a push creeping to the static limit, way {way}: outside the window of {window:e} N the \
             force is off the push by {:e} N at most and the mass at {:e} m/s at most; inside it \
             {:e} N and {:e} m/s, {updates} updates at most in a step; let go with the push \
             {let_go:?} N inside the limit",
            off[0], speed[0], off[1], speed[1]
        );
        // Held to rounding, either way, outside the window and inside it,
        // to the end: 5e-8 N inside the limit.
        assert!(off[0] <= 1.0e-12 && speed[0] <= 1.0e-15, "way {way}: {off:?} {speed:?}");
        assert!(off[1] <= 1.0e-12 && speed[1] <= 1.0e-15, "way {way}: {off:?} {speed:?}");
        assert_eq!(let_go, None, "way {way}");
    }
}

// The end of a step on the load line.

/// What an [`ExactOscillator`] says of how its flow ends a step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Ends {
    /// It keeps the default promise, the midpoint relation, which it does
    /// not have.
    Promised,
    /// It returns the end of its flow on the signal output `p_end`.
    Told,
    /// It says that it is neither.
    Untold,
}

/// A neighbour that is NOT stepped by the midpoint rule: two undamped
/// oscillators of unit modal mass driven through one port, each stepped
/// exactly under the held force, as a modal structure is. Port `p` receives
/// the force [N] and returns the mean velocity of the step [m/s], which is
/// the displacement the force did its work on; the velocity at the step's
/// end is in the state, and is not twice that mean less the start.
struct ExactOscillator {
    ends: Ends,
    ports: Vec<PortSpec>,
}

/// Angular frequencies [rad/s] and shapes at the port [1/sqrt(kg)].
const OMEGAS: [f64; 2] = [40.0, 170.0];
const SHAPES: [f64; 2] = [0.8, -1.3];

/// Per mode: displacement and velocity.
type Modes = [(f64, f64); 2];

impl ExactOscillator {
    fn new(ends: Ends) -> Self {
        let mut ports = vec![
            PortSpec::power(
                StableId::new("p").unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::EffortIn,
                core::num::NonZeroUsize::MIN,
            )
            .unwrap(),
        ];
        if ends == Ends::Told {
            ports.push(PortSpec::signal_out(
                StableId::new("p_end").unwrap(),
                core::num::NonZeroUsize::MIN,
            ));
        }
        Self { ends, ports }
    }

    /// The velocity at the port [m/s] in `modes`.
    fn velocity(modes: &Modes) -> f64 {
        SHAPES[0] * modes[0].1 + SHAPES[1] * modes[1].1
    }

    /// The port's mean and end velocity over a step of `dt_s` from rest
    /// under a held newton: its two mobilities [m/(N s)].
    fn mobilities(dt_s: f64) -> (f64, f64) {
        let (mut mean, mut end) = (0.0, 0.0);
        for (omega, shape) in OMEGAS.iter().zip(SHAPES) {
            mean += shape * shape * (1.0 - (omega * dt_s).cos()) / (omega * omega * dt_s);
            end += shape * shape * (omega * dt_s).sin() / omega;
        }
        (mean, end)
    }
}

impl Component for ExactOscillator {
    type State = Modes;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, modes: &Modes) -> Result<f64, ElementError> {
        Ok(modes
            .iter()
            .zip(OMEGAS)
            .map(|((q, v), omega)| 0.5 * (v * v + omega * omega * q * q))
            .sum())
    }

    fn end_of_step(&self, _port: usize) -> fs_ports::component::EndOfStep {
        match self.ends {
            Ends::Promised => fs_ports::component::EndOfStep::Midpoint,
            Ends::Told => fs_ports::component::EndOfStep::Told(1),
            Ends::Untold => fs_ports::component::EndOfStep::Untold,
        }
    }

    fn propose(
        &self,
        old: &Modes,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Modes,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let force = inputs[0];
        let (mut moved, mut size) = (0.0_f64, 0.0_f64);
        for mode in 0..2 {
            let (omega, shape, (q, v)) = (OMEGAS[mode], SHAPES[mode], old[mode]);
            let (sin, cos) = (omega * dt_s).sin_cos();
            // The modal force is held: the exact step of an oscillator of
            // unit mass about its shifted rest point. The displacement is
            // advanced by its increment, formed without subtracting the
            // two ends (`1 - cos` as twice the squared sine of half the
            // angle): a difference of ends would round at the size of the
            // displacement, which a port held at rest cannot be solved
            // below.
            let rest = shape * force / (omega * omega);
            let half = (0.5 * omega * dt_s).sin();
            let step = (rest - q) * 2.0 * half * half + v * sin / omega;
            let v1 = -(q - rest) * omega * sin + v * cos;
            new[mode] = (q + step, v1);
            moved += shape * step;
            size += (shape * step).abs();
        }
        outputs[0] = moved / dt_s;
        if self.ends == Ends::Told {
            outputs[1] = Self::velocity(new);
        }
        // The force's work is the change of the modes' energy, exactly.
        let work_j = force * moved;
        let energies = self.stored_energy_j(old)? + self.stored_energy_j(new)?;
        Ok(Proposal::Step(Receipt::lossless(
            work_j,
            64.0 * f64::EPSILON * (force.abs() * size + energies) + f64::MIN_POSITIVE,
        )))
    }
}

/// What a [`Holder`] was handed in its last step, and what it did with it.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Held {
    /// The line's two rows for end `a`: how its mean flow and how its flow
    /// at the step's end answer the effort.
    compliance: [f64; 2],
    /// The effort returned at `a` [N].
    effort: f64,
    /// The slip on the line at that effort: its mean over the step, and at
    /// the step's end.
    mean: f64,
    end: f64,
    /// The work of the effort on the mean slip so far [J], signed.
    work_j: f64,
}

thread_local! {
    /// What every proposal of a [`Holder`] was handed for end `a`: the
    /// released end flow, and the effort it settled on.
    static HANDED: core::cell::RefCell<Vec<(f64, f64)>> =
        const { core::cell::RefCell::new(Vec::new()) };
}

/// The smallest component that needs the end of a step: two ends that it
/// holds together at the END of every step, as a side carries what sticks
/// to it. Ports `a` and `b` receive flows (minus the movers', as every
/// member of a junction does) and return the effort and its opposite;
/// `a_end` is the end input of `a`. End `b` is for a side whose flow is
/// given: it is held over its step, so it ends the step with its mean.
///
/// Held at no END slip, the mean slip of a step is not zero where a mover's
/// mean and end differ, and the effort does a work on it of either sign.
/// That work is a defect of the step, not heat, and is booked as one.
struct Holder {
    ports: [PortSpec; 3],
}

impl Holder {
    fn new() -> Self {
        let end = |name: &str| {
            PortSpec::power(
                StableId::new(name).unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                core::num::NonZeroUsize::MIN,
            )
            .unwrap()
        };
        Self {
            ports: [
                end("a"),
                end("b"),
                PortSpec::signal_in(
                    StableId::new("a_end").unwrap(),
                    core::num::NonZeroUsize::MIN,
                ),
            ],
        }
    }
}

impl Component for Holder {
    type State = Held;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::NonSmooth
    }

    fn stored_energy_j(&self, _held: &Held) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Held,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Held,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        self.propose_on_line(old, dt_s, inputs, &[0.0; 6], new, outputs, cx)
    }
}

impl LoadLine for Holder {
    fn line_ports(&self) -> &[usize] {
        &[0, 1]
    }

    fn line_ends(&self) -> &[Option<usize>] {
        &[Some(2)]
    }

    fn propose_on_line(
        &self,
        old: &Held,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Held,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        // The efforts are (e, -e) on the two ends. The line's square
        // matrix is its first two rows; the third is the end of `a`.
        let (flow_a, flow_b, end_a) = (inputs[0], inputs[1], inputs[2]);
        let c = compliance[0] - compliance[1] - compliance[2] + compliance[3];
        let c_end = compliance[4] - compliance[5];
        // The effort that ends the step with no slip: `b` ends it with its
        // mean.
        let effort = if c_end > 0.0 {
            (end_a - flow_b) / c_end
        } else {
            0.0
        };
        let mean = flow_a - flow_b - c * effort;
        let work_j = dt_s * effort * mean;
        *new = Held {
            compliance: [c, c_end],
            effort,
            mean,
            end: end_a - flow_b - c_end * effort,
            work_j: old.work_j + work_j,
        };
        HANDED.with(|handed| handed.borrow_mut().push((end_a, effort)));
        outputs[0] = effort;
        outputs[1] = -effort;
        let scale = dt_s * effort.abs() * (flow_a.abs() + flow_b.abs() + (c * effort).abs());
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: work_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

/// Both modes moving, so that the port moves at 0.5 m/s.
fn moving() -> Modes {
    [(0.0, 0.3), (0.0, (0.3 * SHAPES[0] - 0.5) / -SHAPES[1])]
}

/// The oscillator on end `a` of a [`Holder`] whose end `b` is a side at
/// `carried_m_s`.
fn held_at_its_end(start: Modes, carried_m_s: f64) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Told), start)
        .unwrap();
    let holder = builder
        .add_on_line("holder", Holder::new(), Held::default())
        .unwrap();
    let held = port(&builder, holder, "a");
    builder
        .common_flow(port(&builder, modes, "p"), &[held])
        .unwrap();
    // A port receives minus the velocity of what moves it.
    builder
        .boundary(
            port(&builder, holder, "b"),
            Box::new(move |_, _, values| values[0] = -carried_m_s),
        )
        .unwrap();
    (builder.seal().unwrap(), modes, holder)
}

#[test]
fn the_end_row_of_the_line_is_the_neighbours_end_mobility_and_carries_it_at_the_steps_end() {
    // A neighbour that is exact for a held force returns the mean velocity
    // of its step, and tells the velocity it ends the step with on a
    // signal. The kernel wires that to the end input of the component on
    // the line, with the sign of the port's own junction, and finds how it
    // answers the effort by the factors that give the line itself: the
    // line's second row. The holder's effort then ends every step with the
    // port at the carried speed, which the mean row alone cannot say.
    //
    // Both rows against the closed form of two exact oscillators. The
    // neighbour is differenced, so they are known to what a forward
    // difference leaves: 1.4e-6 of themselves at worst here (measured on
    // ts2 at 3a42d0c16d), where the block's efforts are of the size of the
    // step it is differenced by. Their ratio is not the 2 of a mass.
    for dt in [1.0e-3_f64, 4.0e-3] {
        let (mut system, modes, holder) = held_at_its_end(moving(), 0.1);
        let (mean_mobility, end_mobility) = ExactOscillator::mobilities(dt);
        let mut most_off = 0.0_f64;
        let mut worst_row = [0.0_f64; 2];
        let (mut iterations, mut most_iterations) = (0_u32, 0_u32);
        let mut before = Held::default();
        with_cx(|cx| {
            for step in 0..400 {
                HANDED.with(|handed| handed.borrow_mut().clear());
                let start = *system.state::<Modes>(modes).unwrap();
                let report = system
                    .step(dt, cx)
                    .unwrap_or_else(|error| panic!("dt {dt}, step {step}: {error}"));
                iterations += report.iterations;
                most_iterations = most_iterations.max(report.iterations);
                // The first proposal of every step, the first step's
                // included, is handed the RELEASED end of this step, as it
                // is handed the released flow: the end the oscillator
                // reaches over the step with no force, in the port's sign,
                // the closed form from the modes the step starts with. The
                // kernel holds the last effort and its end row puts back
                // what that effort takes away, so it is known to what the
                // row's difference leaves of that part (the row is asserted
                // to 1e-4 below) and to rounding.
                let free_end: f64 = (0..2)
                    .map(|mode| {
                        let (omega, (q, v)) = (OMEGAS[mode], start[mode]);
                        let (sin, cos) = (omega * dt).sin_cos();
                        SHAPES[mode] * (-q * omega * sin + v * cos)
                    })
                    .sum();
                let (first, released) = (HANDED.with(|handed| handed.borrow()[0].0), -free_end);
                let taken = before.compliance[1] * before.effort;
                let rounding = 64.0 * f64::EPSILON * (released.abs() + taken.abs());
                assert!(
                    (first - released).abs() <= 1e-4 * taken.abs() + rounding,
                    "dt {dt}, step {step}: handed {first}, the released end is {released}"
                );
                before = *system.state::<Held>(holder).unwrap();
                let held = *system.state::<Held>(holder).unwrap();
                let end = ExactOscillator::velocity(system.state::<Modes>(modes).unwrap());
                most_off = most_off.max((end - 0.1).abs());
                for (row, closed) in [mean_mobility, end_mobility].into_iter().enumerate() {
                    let off = (held.compliance[row] / closed - 1.0).abs();
                    worst_row[row] = worst_row[row].max(off);
                    assert!(
                        off <= 1e-4,
                        "dt {dt}, step {step}: {held:?} against {closed}"
                    );
                }
                // What the port's end input shows after the step is the
                // end flow that crossed, in the port's sign.
                let shown = system.last_input(system.port(holder, "a_end").unwrap())[0];
                assert!(
                    (shown + end).abs() <= 1e-12,
                    "dt {dt}, step {step}: {shown} and {end}"
                );
                assert!(
                    held.end.abs() <= 1e-9,
                    "dt {dt}, step {step}: {held:?}, {end}"
                );
            }
        });
        let held = *system.state::<Held>(holder).unwrap();
        println!(
            "carried at 0.1 m/s at the end of every step, dt {dt}: the port ends its steps within \
             {most_off:.3e} m/s of that (it started at 0.5); the line's rows are {:?} against the \
             closed form's ({mean_mobility}, {end_mobility}), off by {worst_row:?} at most \
             (relative); end over mean {:.6} (a mass: 2); the holder's work on the mean slip, \
             its defect, is {:+.6e} J in all; ledger residual {:.3e} J; {iterations} kernel \
             iterations in 400 steps, {most_iterations} at most in one",
            held.compliance,
            end_mobility / mean_mobility,
            held.work_j,
            system.ledger().energy_residual_j()
        );
        assert!(most_off <= 1e-12, "dt {dt}: {most_off}");
        assert!(
            end_mobility / mean_mobility < 1.999,
            "the neighbour is not a midpoint one"
        );
    }
}

#[test]
fn carried_at_its_port_two_exact_modes_swing_at_the_constrained_frequency() {
    // With the port's velocity given, what is left of two modes is one.
    // The combination `u = s2 q1 - s1 q2` swings at
    // `Omega^2 = (s2^2 w1^2 + s1^2 w2^2) / (s1^2 + s2^2)` about a motion
    // that is linear in time: each mode is drawn out at `a_k = s_k l / w_k^2`
    // with `l = v / sum(s_k^2 / w_k^2)`, so `u` drifts at `s2 a1 - s1 a2`.
    // The discrete motion, carried at the END of every step, has that
    // period to second order in the step: 2.7e-3, 6.8e-4 and 1.7e-4 of it
    // short at steps of 2, 1 and 0.5 ms (measured; asserted as a ratio
    // between 3.5 and 4.5 from one to the next).
    let (s1, s2) = (SHAPES[0], SHAPES[1]);
    let (w1, w2) = (OMEGAS[0] * OMEGAS[0], OMEGAS[1] * OMEGAS[1]);
    let omega = ((s2 * s2 * w1 + s1 * s1 * w2) / (s1 * s1 + s2 * s2)).sqrt();
    let period = 2.0 * PI / omega;
    let carried = 0.1;
    let pull = carried / (s1 * s1 / w1 + s2 * s2 / w2);
    let drift = s2 * (s1 * pull / w1) - s1 * (s2 * pull / w2);
    let mut errors = Vec::new();
    for dt in [2.0e-3_f64, 1.0e-3, 5.0e-4] {
        let (mut system, modes, _) = held_at_its_end(moving(), carried);
        let swing = |system: &System| {
            let state = system.state::<Modes>(modes).unwrap();
            s2 * state[0].1 - s1 * state[1].1 - drift
        };
        // Upward zero crossings of the free combination's velocity about
        // its drift.
        let mut crossings = Vec::new();
        let mut before = (0.0, swing(&system));
        let mut taken = 0_u32;
        with_cx(|cx| {
            while crossings.len() < 12 {
                taken += 1;
                assert!(
                    taken < 20_000,
                    "dt {dt}: {} crossings in {taken} steps",
                    crossings.len()
                );
                let report = system.step(dt, cx).unwrap();
                let now = (report.time_s, swing(&system));
                if before.1 < 0.0 && now.1 >= 0.0 {
                    crossings.push(before.0 + (now.0 - before.0) * before.1 / (before.1 - now.1));
                }
                before = now;
            }
        });
        let measured = (crossings[11] - crossings[1]) / 10.0;
        errors.push((measured - period) / period);
        println!(
            "two exact modes carried at the port, dt {dt}: period {measured:.9} s against \
             {period:.9} s, {:+.3e} of it",
            (measured - period) / period
        );
    }
    assert!(errors[0].abs() < 5e-3, "{errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(ratio > 3.5 && ratio < 4.5, "{errors:?}");
    }
}

#[test]
fn an_end_held_at_exactly_nothing_has_no_size_and_is_solved_to_what_the_solve_is_told() {
    // A stated limit. A wire is solved relative to its size, and to the
    // rounding of the terms its driver sums. The end of a flow that is
    // held at exactly zero has no size, and its rounding is not in those
    // terms: the neighbour forms it as a sum over its modes, which cancel
    // at the port and do not vanish in themselves (0.24 and 0.26 m/s
    // here). The kernel is not told that sum. With the solve's default,
    // no absolute tolerance, the block is refused within a few steps
    // (which step hangs on the rounding of that sum); told an absolute
    // tolerance above the rounding of those terms, it is taken, and the
    // port ends its steps at rest to that tolerance.
    let dt = 1.0e-3;
    let (mut system, _, _) = held_at_its_end(moving(), 0.0);
    let refused = with_cx(|cx| {
        (0..400).find_map(|step| match system.step(dt, cx) {
            Ok(_) => None,
            Err(SystemError::NotConverged {
                worst, residual, ..
            }) => Some((step, worst, residual)),
            Err(other) => panic!("step {step}: {other}"),
        })
    });
    let (step, worst, residual) = refused.expect("refused within 400 steps");
    println!(
        "an end held at exactly nothing, default tolerances: refused at step {step} on `{worst}` \
         with a residual of {residual:.3e} m/s"
    );
    assert_eq!(worst, "holder.a_end");
    assert!(residual.abs() < 1e-15, "{residual}");

    let (mut system, modes, _) = held_at_its_end(moving(), 0.0);
    system.set_config(SolveConfig {
        absolute_tolerance: 1.0e-15,
        ..SolveConfig::default()
    });
    let mut most_end = 0.0_f64;
    with_cx(|cx| {
        for step in 0..400 {
            system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
            let end = ExactOscillator::velocity(system.state::<Modes>(modes).unwrap());
            most_end = most_end.max(end.abs());
        }
    });
    println!(
        "the same with an absolute tolerance of 1e-15: taken for 400 steps, the port ends them \
         at {most_end:.3e} m/s at most"
    );
    assert!(most_end <= 1e-13, "{most_end}");
}

#[test]
fn a_neighbour_that_keeps_the_midpoint_promise_it_does_not_have_is_never_brought_to_rest() {
    // What the end row is for. The exact oscillator that keeps the default
    // promise seals beside the slider, and the slider, which takes the
    // step's end for twice the mean less the start, arrests it and holds
    // it at no MEAN slip. Its record says at rest. The port's own velocity
    // at the end of each step changes sign from step to step and does not
    // decay, and the force that holds it alternates with it.
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Promised), moving())
        .unwrap();
    // Pressed with 10 kN, so that it arrests the port within its first
    // step and holds whatever the modes then ask of it.
    let strong = Slider::new(1.0e4, HOLDS, SLIDES).unwrap();
    let rubbing = builder
        .add_on_line("rubbing", strong, Slip::default())
        .unwrap();
    builder
        .common_flow(port(&builder, modes, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "b"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "heat"),
            Box::new(|_, _, values| values[0] = 300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let dt = 1.0e-3;
    let mut ends = Vec::new();
    with_cx(|cx| {
        for step in 0..200 {
            system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
            ends.push(ExactOscillator::velocity(
                system.state::<Modes>(modes).unwrap(),
            ));
        }
    });
    let flips = ends[20..]
        .windows(2)
        .filter(|pair| pair[0] * pair[1] < 0.0)
        .count();
    let (early, late) = (
        ends[20..40]
            .iter()
            .fold(0.0_f64, |most, end| most.max(end.abs())),
        ends[180..]
            .iter()
            .fold(0.0_f64, |most, end| most.max(end.abs())),
    );
    println!(
        "an exact oscillator that keeps the midpoint promise, beside the slider: the port's end \
         velocity changes sign in {flips} of 179 steps; it is up to {early:.3e} m/s in steps 21 \
         to 40 and {late:.3e} m/s in steps 181 to 200"
    );
    assert!(flips >= 170, "{flips}");
    assert!(late > 0.5 * early && early > 1e-3, "{early} {late}");
}

#[test]
#[allow(clippy::too_many_lines)]
fn a_line_port_and_its_neighbours_must_agree_on_how_a_step_ends() {
    let refused = |builder: SystemBuilder, says: &str| match builder.seal() {
        Err(AdmissionError::LineEnd { port, reason }) => {
            assert!(reason.contains(says), "{port}: {reason}");
        }
        Err(other) => panic!("{other}"),
        Ok(_) => panic!("sealed: expected a refusal that says `{says}`"),
    };
    let on_the_ground = |builder: &mut SystemBuilder, rubbing: NodeId| {
        builder
            .boundary(
                port(builder, rubbing, "b"),
                Box::new(|_, _, values| values[0] = 0.0),
            )
            .unwrap();
        builder
            .boundary(
                port(builder, rubbing, "heat"),
                Box::new(|_, _, values| values[0] = 300.0),
            )
            .unwrap();
    };

    // A neighbour that says it is neither, beside the slider.
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Untold), moving())
        .unwrap();
    let rubbing = builder
        .add_on_line("rubbing", slider(), Slip::default())
        .unwrap();
    builder
        .common_flow(port(&builder, modes, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    on_the_ground(&mut builder, rubbing);
    refused(builder, "does not tell its end");

    // One that tells its end, beside a port that takes none.
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Told), moving())
        .unwrap();
    let rubbing = builder
        .add_on_line("rubbing", slider(), Slip::default())
        .unwrap();
    builder
        .common_flow(port(&builder, modes, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    on_the_ground(&mut builder, rubbing);
    refused(builder, "the port takes no end input");

    // A port that takes an end, beside a mass, which is midpoint.
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.5,
            },
        )
        .unwrap();
    let holder = builder
        .add_on_line("holder", Holder::new(), Held::default())
        .unwrap();
    builder
        .common_flow(port(&builder, mass, "p"), &[port(&builder, holder, "a")])
        .unwrap();
    builder
        .boundary(
            port(&builder, holder, "b"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    refused(builder, "no neighbour tells the end");

    // A port whose flow is the sum of one that tells and one that does not.
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Told), moving())
        .unwrap();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let holder = builder
        .add_on_line("holder", Holder::new(), Held::default())
        .unwrap();
    builder
        .common_effort(
            port(&builder, holder, "a"),
            &[port(&builder, modes, "p"), port(&builder, mass, "p")],
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, holder, "b"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    refused(builder, "others do not");

    // An end input that something is connected to.
    let mut builder = SystemBuilder::new();
    let modes = builder
        .add("modes", ExactOscillator::new(Ends::Told), moving())
        .unwrap();
    let holder = builder
        .add_on_line("holder", Holder::new(), Held::default())
        .unwrap();
    builder
        .boundary(
            port(&builder, holder, "b"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    let body = builder
        .add_with_tangent(
            "body",
            HeatCapacity::new(BODY_J_K).unwrap().publishing(),
            300.0,
        )
        .unwrap();
    builder
        .common_flow(port(&builder, modes, "p"), &[port(&builder, holder, "a")])
        .unwrap();
    builder
        .signal(
            port(&builder, body, "temperature"),
            port(&builder, holder, "a_end"),
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, body, "heat"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    refused(builder, "the kernel makes that wire");
}

#[test]
fn a_component_with_an_end_row_conforms_on_its_line_and_is_handed_the_row() {
    // The conformance check of the load-line tier takes the line's square
    // matrix followed by the rows of the end inputs, as the kernel hands
    // them over, and refuses the square matrix alone by its length.
    let holder = Holder::new();
    // `a` answers at 0.004 for its mean and 0.0078 for its end; `b` is
    // given. Released: `a` at -0.3 (mean) and -0.31 (end), `b` at -0.1.
    let matrix = [0.004, 0.0, 0.0, 0.0, 0.0078, 0.0];
    let inputs = [-0.3, -0.1, -0.31];
    let state = Held::default();
    let stale = Held {
        work_j: 7.0,
        ..Held::default()
    };
    let case = Case {
        state: &state,
        stale: &stale,
        dt_s: 1.0e-3,
        inputs: &inputs,
    };
    let seen = with_cx(|cx| check_load_line(&holder, &case, &matrix, cx)).unwrap();
    let Observed::Stepped { outputs, .. } = seen.observed else {
        panic!("the holder takes the step");
    };
    // No slip at the step's end: `(-0.31 + 0.1) / 0.0078`.
    let effort = (-0.31 + 0.1) / 0.0078;
    assert_eq!(outputs, [effort, -effort]);
    assert!(
        (seen.flows[0] - (-0.3 - 0.004 * effort)).abs() <= 1e-15,
        "{:?}",
        seen.flows
    );
    let short = with_cx(|cx| check_load_line(&holder, &case, &matrix[..4], cx));
    assert!(
        matches!(
            short,
            Err(fs_ports::component::conformance::Failure::Contract(
                fs_ports::component::ContractViolation::VectorLength {
                    expected: 6,
                    actual: 4,
                    ..
                }
            ))
        ),
        "the square matrix alone was not refused by its length"
    );
}

/// The matrix a probe on a line is handed whose `a` is on a mass of 2 kg
/// that is DIFFERENCED (added without its tangent), moving at
/// `velocity_m_s` under a push of `push_n`, in a step of `dt`; `b` is given.
/// With it, the entries the step left unresolved.
fn handed_beside_a_differenced_mass(push_n: f64, velocity_m_s: f64, dt: f64) -> ([f64; 4], u32) {
    let mut builder = SystemBuilder::new();
    let moving = Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let one = builder.add("one", Mass::new(2.0).unwrap(), moving).unwrap();
    let probe = builder
        .add_on_line("probe", Probe::new(), [f64::NAN; 4])
        .unwrap();
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(push_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow(
            port(&builder, one, "p"),
            &[port(&builder, probe, "a"), port(&builder, push, "p")],
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, probe, "b"),
            Box::new(|_, _, values| values[0] = 0.3),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let report = with_cx(|cx| system.step(dt, cx)).unwrap();
    (*system.state::<[f64; 4]>(probe).unwrap(), report.unresolved)
}

#[test]
fn a_differenced_neighbour_is_handed_its_mobility_whatever_the_force_on_it() {
    // A mass that moves needs no force to keep moving, and the kernel
    // differences a component that declares no tangent at a step of
    // `sqrt(eps)` of the force on it, or of the largest force of its block.
    // Under a small force that step moved the mass's mean velocity by less
    // than its rounding: the line was handed no mobility at all at 1 m/s
    // under 1e-6 N and at 100 m/s under 1e-3 N, and 1.3e-2 too much at 1
    // m/s under 1e-3 N; 12 of these 18 cases were off by more than 1e-6,
    // 7 handed nothing. A differenced column is resolved against the
    // rounding of its outputs, and no entry is left unresolved.
    let dt = 1.0e-3_f64;
    let own = 0.5 * dt / 2.0;
    let mut worst = 0.0_f64;
    for velocity in [0.0_f64, 1.0, 100.0] {
        for push in [1.0_f64, 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-12, 0.0] {
            let (matrix, unresolved) = handed_beside_a_differenced_mass(push, velocity, dt);
            let off = (matrix[0] / own - 1.0).abs();
            assert!(
                off <= 1.0e-6,
                "{velocity} m/s under {push:e} N: handed {:e} for {own:e}",
                matrix[0]
            );
            assert_eq!(&matrix[1..], &[0.0; 3][..], "`b` is given");
            assert_eq!(unresolved, 0, "{velocity} m/s under {push:e} N");
            worst = worst.max(off);
        }
    }
    eprintln!(
        "a differenced mass of 2 kg at 0, 1 and 100 m/s under 1 N down to 1e-12 N and none: \
         its mobility handed within {worst:.3e} of dt / (2 m), no entry unresolved"
    );
}

/// What a [`Lopsided`] port `a` moves at under no force [m/s].
const LOPSIDED_SPEED_M_S: f64 = 100.0;
/// How a [`Lopsided`] port `a` answers the force on it, where it does
/// [m/(N s)].
const LOPSIDED_MOBILITY_A: f64 = 2.5e-4;
/// How a [`Lopsided`] port `b` answers the force on it [m/(N s)].
const LOPSIDED_MOBILITY_B: f64 = 1.0;

/// A two-port component that declares no tangent, so that the kernel
/// differences it, whose outputs answer the force at `a` on two scales:
/// `a` moves at [`LOPSIDED_SPEED_M_S`] plus `mobility_a` times that force,
/// plus `kink_mobility` times its excess over `kink_above_n`, plus `cubic`
/// times the cube of its departure from `cubic_about_n` [m/(N^3 s)], and
/// `jump_m_s` faster once the force is above `jump_above_n`, as much slower
/// once it is below `jump_below_n`; `b` at `curve` times its square [m/(N^2
/// s)], plus [`LOPSIDED_MOBILITY_B`] times the force at `b` (without it the
/// line would not be passive). Under a force at `a` above
/// `shortens_above_n` it asks to shorten the step, as a component that
/// locates an event does. What work its ports do it stores.
struct Lopsided {
    ports: [PortSpec; 2],
    mobility_a: f64,
    curve: f64,
    kink_above_n: f64,
    kink_mobility: f64,
    cubic: f64,
    cubic_about_n: f64,
    jump_above_n: f64,
    jump_below_n: f64,
    jump_m_s: f64,
    shortens_above_n: f64,
}

impl Lopsided {
    /// Without a kink, a cube, a jump or an event.
    fn new(mobility_a: f64, curve: f64) -> Self {
        let end = |name: &str| {
            PortSpec::power(
                StableId::new(name).unwrap(),
                PortKind::MechanicalForceVelocity,
                PortRole::EffortIn,
                core::num::NonZeroUsize::MIN,
            )
            .unwrap()
        };
        Self {
            ports: [end("a"), end("b")],
            mobility_a,
            curve,
            kink_above_n: f64::INFINITY,
            kink_mobility: 0.0,
            cubic: 0.0,
            cubic_about_n: 0.0,
            jump_above_n: f64::INFINITY,
            jump_below_n: f64::NEG_INFINITY,
            jump_m_s: 0.0,
            shortens_above_n: f64::INFINITY,
        }
    }
}

impl Component for Lopsided {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, stored_j: &f64) -> Result<f64, ElementError> {
        Ok(*stored_j)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        let (at_a, at_b) = (inputs[0], inputs[1]);
        if at_a > self.shortens_above_n {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }
        let kink_m_s = self.kink_mobility * (at_a - self.kink_above_n).max(0.0);
        let cubic_m_s = self.cubic * (at_a - self.cubic_about_n).powi(3);
        let jump_m_s = if at_a > self.jump_above_n {
            self.jump_m_s
        } else if at_a < self.jump_below_n {
            -self.jump_m_s
        } else {
            0.0
        };
        outputs[0] = LOPSIDED_SPEED_M_S + self.mobility_a * at_a + kink_m_s + cubic_m_s + jump_m_s;
        outputs[1] = self.curve * at_a * at_a + LOPSIDED_MOBILITY_B * at_b;
        let work_j = dt_s * (at_a * outputs[0] + at_b * outputs[1]);
        *new = old + work_j;
        let tolerance_j = 1.0e-12 * (old.abs() + new.abs() + work_j.abs()) + f64::MIN_POSITIVE;
        Ok(Proposal::Step(Receipt::lossless(*new - old, tolerance_j)))
    }
}

/// The matrix a probe on a line is handed beside `lopsided`, its `a` under
/// `push_n`, in the second step of 1 ms, and the entries that step left
/// unresolved. The line is formed where its step starts: in the first step
/// the force at `a` is still the none it started with, in the second it is
/// the push.
fn handed_beside_a_lopsided(lopsided: Lopsided, push_n: f64) -> ([f64; 4], u32) {
    let mut builder = SystemBuilder::new();
    let lopsided = builder.add("lopsided", lopsided, 0.0).unwrap();
    let probe = builder
        .add_on_line("probe", Probe::new(), [f64::NAN; 4])
        .unwrap();
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(push_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow(
            port(&builder, lopsided, "a"),
            &[port(&builder, probe, "a"), port(&builder, push, "p")],
        )
        .unwrap();
    builder
        .common_flow(port(&builder, lopsided, "b"), &[port(&builder, probe, "b")])
        .unwrap();
    let mut system = builder.seal().unwrap();
    let report = with_cx(|cx| {
        system.step(1.0e-3, cx).unwrap();
        system.step(1.0e-3, cx).unwrap()
    });
    (*system.state::<[f64; 4]>(probe).unwrap(), report.unresolved)
}

#[test]
fn each_entry_of_a_differenced_column_keeps_the_step_that_resolved_it() {
    // One input of a component can move one of its outputs past the
    // output's rounding and leave another inside it. Under 1e3 N the
    // `sqrt(eps)` step of the force at `a` of a `Lopsided` moves `b`, 1e-5
    // m/(N^2 s) times its square, at once, and `a`, 100 m/s and a mobility,
    // by 2^18 units in its last place, under the 2^20 of a resolved entry.
    // The column is differenced again at a step 401 times longer, 6.0e-3 N;
    // `b` differenced there would be its secant over that step, 3.0e-6 off
    // its slope, where the first step leaves it 7.5e-9 off. Each entry keeps
    // the difference that first resolved it. The force at `b` moves `b` and
    // not `a`: that entry is nothing at every step, and is the one counted.
    let (push_n, curve) = (1.0e3, 1.0e-5);
    let (matrix, unresolved) =
        handed_beside_a_lopsided(Lopsided::new(LOPSIDED_MOBILITY_A, curve), push_n);
    let slope_b = 2.0 * curve * push_n;
    eprintln!(
        "a lopsided neighbour: handed {matrix:?} for its own [{LOPSIDED_MOBILITY_A:e}, 0, \
         {slope_b:e}, {LOPSIDED_MOBILITY_B:e}], {unresolved} unresolved"
    );
    assert!(
        (matrix[0] / LOPSIDED_MOBILITY_A - 1.0).abs() <= 1.0e-6,
        "{matrix:?}"
    );
    assert_eq!(matrix[1], 0.0, "{matrix:?}");
    assert!((matrix[2] / slope_b - 1.0).abs() <= 1.0e-6, "{matrix:?}");
    assert!(
        (matrix[3] / LOPSIDED_MOBILITY_B - 1.0).abs() <= 1.0e-6,
        "{matrix:?}"
    );
    assert_eq!(unresolved, 1);
}

#[test]
fn a_longer_step_does_not_difference_across_a_jump() {
    // An output that does not answer an input changes under no step, and
    // the kernel cannot tell it from one whose change rounds away: it takes
    // the step longer. Where the output jumps further on, the longer step
    // differences across the jump; where it jumps the other way as far on
    // the other side of the iterate, the step there agrees. (Found while
    // the block solve's columns were resolved too, on a scratch kernel: a
    // sliding slider's force, which does not answer its released flow, was
    // given columns of 120.5 and 25.2 once the longer step carried it to
    // sticking, and the drawn graphs of this file took 13.6 per cent more
    // updates.) A longer step that changes an output by more than 2^-13 of
    // it is not taken for that output. Here `a` does not answer the force on
    // it within 1e-4 N of the push of 1e-6 N, and is 1 m/s faster above and
    // slower below; the longest step, 2^-10 N, crosses. That entry and the
    // nothing of `a` under the force at `b` are left unresolved.
    let (push_n, curve) = (1.0e-6, 10.0);
    let mut lopsided = Lopsided::new(0.0, curve);
    lopsided.jump_above_n = push_n + 1.0e-4;
    lopsided.jump_below_n = push_n - 1.0e-4;
    lopsided.jump_m_s = 1.0;
    let (matrix, unresolved) = handed_beside_a_lopsided(lopsided, push_n);
    eprintln!(
        "a lopsided neighbour that jumps 1e-4 N off its push: handed {matrix:?}, {unresolved} unresolved"
    );
    assert_eq!(matrix[0], 0.0, "{matrix:?}");
    assert!(
        (matrix[2] / (2.0 * curve * push_n) - 1.0).abs() <= 1.0e-6,
        "{matrix:?}"
    );
    assert_eq!(unresolved, 2);
}

#[test]
fn a_longer_step_across_a_kink_is_not_taken() {
    // An output whose slope changes inside the longer step, beyond the
    // first: the difference at the longer step is the secant across the
    // kink, resolved against the output's rounding all the same. It is
    // differenced at that step on the other side of the iterate as well,
    // where it has its slope, and the two do not agree: the entry is left
    // as the first step had it, and counted. Here `a`, 2.5e-4 m/(N s) under
    // a push of 1e-6 N, answers ten times that beyond 1e-4 N more; the
    // first step is 1.5e-14 N, the longest 2^-10 N, and the secant there is
    // 9.08 times the slope. The first step's change rounded away: the
    // entry is nothing. Counted with it, the nothing of `a` under the force
    // at `b`.
    let push_n = 1.0e-6;
    let mut lopsided = Lopsided::new(LOPSIDED_MOBILITY_A, 10.0);
    lopsided.kink_above_n = push_n + 1.0e-4;
    lopsided.kink_mobility = 9.0 * LOPSIDED_MOBILITY_A;
    let (matrix, unresolved) = handed_beside_a_lopsided(lopsided, push_n);
    let widest = 1.0_f64 / 1024.0;
    let secant = LOPSIDED_MOBILITY_A + 9.0 * LOPSIDED_MOBILITY_A * (1.0 - 1.0e-4 / widest);
    eprintln!(
        "a lopsided neighbour with a kink 1e-4 N off its push: handed {matrix:?}, \
         {unresolved} unresolved; its slope {LOPSIDED_MOBILITY_A:e}, the secant {secant:e}"
    );
    assert_eq!(matrix[0], 0.0, "{matrix:?}");
    assert_eq!(unresolved, 2);
}

#[test]
fn a_longer_step_is_no_longer_than_2_10_of_the_input_scale() {
    // An output curved the same way on both sides of the iterate (odd about
    // it) gives two differences that agree, and only the step's length
    // bounds the secant's error. Here `a`, 2.5e-4 m/(N s) under a push of
    // 1e-6 N, has 5e-5 m/(N^3 s) times the cube of the force's departure
    // from the push as well. The step aimed at its rounding is some 6e-3 N,
    // and there the entry would be 7.1e-6 off its slope; no step is longer
    // than 2^-10 of the force's scale or of one newton. At 2^-10 N the cube
    // puts it 1.9e-7 off, and it is handed 1.6e-7 off, the rest the
    // rounding of the 100 m/s that `a` moves at.
    let push_n = 1.0e-6;
    let mut lopsided = Lopsided::new(LOPSIDED_MOBILITY_A, 10.0);
    lopsided.cubic = 5.0e-5;
    lopsided.cubic_about_n = push_n;
    let (matrix, unresolved) = handed_beside_a_lopsided(lopsided, push_n);
    eprintln!(
        "a lopsided neighbour with a cube about its push: handed {matrix:?}, {unresolved} \
         unresolved"
    );
    assert!(
        (matrix[0] / LOPSIDED_MOBILITY_A - 1.0).abs() <= 1.0e-6,
        "{matrix:?}"
    );
    assert_eq!(unresolved, 1);
}

#[test]
fn a_longer_step_is_as_long_as_its_entry_needs() {
    // A longer step is aimed at `2^-26` of the output it resolves; it is not
    // taken as long as it may be. Here `a` of a `Lopsided` under 1e3 N, 100
    // m/s and 2.5e-4 m/(N s), has 1e-8 m/(N^3 s) times the cube of the
    // force's departure from the push as well (odd about it, so the other
    // side agrees). Aimed, the step is some 6.0e-3 N, and the entry 1.4e-9
    // off its slope; at the longest step, 2^-10 of the push, 0.98 N, it
    // would be 3.8e-5 off.
    let push_n = 1.0e3;
    let mut lopsided = Lopsided::new(LOPSIDED_MOBILITY_A, 1.0e-5);
    lopsided.cubic = 1.0e-8;
    lopsided.cubic_about_n = push_n;
    let (matrix, unresolved) = handed_beside_a_lopsided(lopsided, push_n);
    eprintln!(
        "a lopsided neighbour under 1e3 N with a cube about its push: handed {matrix:?}, \
         {unresolved} unresolved"
    );
    assert!(
        (matrix[0] / LOPSIDED_MOBILITY_A - 1.0).abs() <= 1.0e-6,
        "{matrix:?}"
    );
    assert_eq!(unresolved, 1);
}

#[test]
fn a_column_the_component_does_not_take_on_both_sides_is_left_and_counted() {
    // A component that locates an event asks to shorten a step whose input
    // carries it past one, and where the kernel's own move does that the
    // column is differenced on the other side of the iterate (the contract,
    // "Requests to shorten inside a loop"). Its longer steps are taken on
    // that side, but the kernel cannot difference the other to check them:
    // the entry is left as the first step had it, and counted. Here a force
    // at `a` of 1e-6 N asks to shorten anywhere above itself and a
    // femtonewton; the first step's change at `a` rounded away, and the
    // entry is nothing. Counted with it, the nothing of `a` under the force
    // at `b`.
    let push_n = 1.0e-6;
    let mut lopsided = Lopsided::new(LOPSIDED_MOBILITY_A, 10.0);
    lopsided.shortens_above_n = push_n + 1.0e-15;
    let (matrix, unresolved) = handed_beside_a_lopsided(lopsided, push_n);
    eprintln!(
        "a lopsided neighbour that shortens above its push: handed {matrix:?}, {unresolved} \
         unresolved"
    );
    assert_eq!(matrix[0], 0.0, "{matrix:?}");
    assert_eq!(unresolved, 2);
}

/// The matrix a probe on a line is handed, in a step of 1 ms, whose `a` is
/// on a mass of 2 kg at 100 m/s under 1e-6 N, added without its tangent
/// unless `declared`, and whose `b` is on a mass of 3 kg at 1 m/s under
/// 1e-6 N added with its tangent; with the entries left unresolved.
fn handed_beside_two_masses(declared: bool) -> ([f64; 4], u32) {
    let mut builder = SystemBuilder::new();
    let moving = |velocity_m_s| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let one = if declared {
        builder.add_with_tangent("one", Mass::new(2.0).unwrap(), moving(100.0))
    } else {
        builder.add("one", Mass::new(2.0).unwrap(), moving(100.0))
    }
    .unwrap();
    let two = builder
        .add_with_tangent("two", Mass::new(3.0).unwrap(), moving(1.0))
        .unwrap();
    let probe = builder
        .add_on_line("probe", Probe::new(), [f64::NAN; 4])
        .unwrap();
    let mut push = |name: &str| {
        builder
            .add_with_tangent(
                name,
                Force::new(Waveform::constant(1.0e-6).unwrap()),
                SourceAccount::default(),
            )
            .unwrap()
    };
    let (push_one, push_two) = (push("push_one"), push("push_two"));
    builder
        .common_flow(
            port(&builder, one, "p"),
            &[port(&builder, probe, "a"), port(&builder, push_one, "p")],
        )
        .unwrap();
    builder
        .common_flow(
            port(&builder, two, "p"),
            &[port(&builder, probe, "b"), port(&builder, push_two, "p")],
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    let report = with_cx(|cx| system.step(1.0e-3, cx)).unwrap();
    (*system.state::<[f64; 4]>(probe).unwrap(), report.unresolved)
}

#[test]
fn a_neighbour_that_declares_its_tangent_is_handed_to_the_bit_beside_a_differenced_one() {
    // The columns of a neighbour that declares its tangent are its tangent;
    // a differenced neighbour's are resolved, and that touches no other.
    // Beside a differenced mass, a mass added with its tangent hands the
    // line the bits it hands beside the same mass added with its own.
    let (beside, unresolved) = handed_beside_two_masses(false);
    let (declared, none) = handed_beside_two_masses(true);
    eprintln!(
        "beside a differenced mass: handed {beside:?} ({:#018x} at `b`), {unresolved} \
         unresolved; beside a declared one {declared:?} ({:#018x})",
        beside[3].to_bits(),
        declared[3].to_bits()
    );
    assert_eq!(beside[3].to_bits(), declared[3].to_bits());
    assert_eq!(&beside[1..3], &[0.0; 2][..]);
    assert!(
        (beside[0] / (0.5 * 1.0e-3 / 2.0) - 1.0).abs() <= 1.0e-6,
        "{beside:?}"
    );
    assert_eq!((unresolved, none), (0, 0));
}

/// A component, with how many times the kernel asked it to propose (on its
/// line too) and for a tangent.
struct Counted<C> {
    inner: C,
    proposals: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    tangents: std::sync::Arc<std::sync::atomic::AtomicUsize>,
}

impl<C> Counted<C> {
    fn new(inner: C) -> Self {
        Self {
            inner,
            proposals: std::sync::Arc::default(),
            tangents: std::sync::Arc::default(),
        }
    }

    fn counts(&self) -> (CountOf, CountOf) {
        (
            std::sync::Arc::clone(&self.proposals),
            std::sync::Arc::clone(&self.tangents),
        )
    }
}

type CountOf = std::sync::Arc<std::sync::atomic::AtomicUsize>;

fn read(count: &CountOf) -> usize {
    count.load(std::sync::atomic::Ordering::Relaxed)
}

fn tick(count: &CountOf) {
    count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
}

impl Component for Counted<Mass> {
    type State = Motion;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        self.inner.ports()
    }

    fn class(&self) -> Class {
        self.inner.class()
    }

    fn stored_energy_j(&self, state: &Motion) -> Result<f64, ElementError> {
        self.inner.stored_energy_j(state)
    }

    fn propose(
        &self,
        old: &Motion,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Motion,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        tick(&self.proposals);
        self.inner.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl fs_ports::component::PortTangent for Counted<Mass> {
    fn jvp(
        &self,
        old: &Motion,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        tick(&self.tangents);
        fs_ports::component::PortTangent::jvp(
            &self.inner,
            old,
            dt_s,
            inputs,
            d_inputs,
            d_outputs,
            cx,
        )
    }

    fn vjp(
        &self,
        old: &Motion,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        fs_ports::component::PortTangent::vjp(
            &self.inner,
            old,
            dt_s,
            inputs,
            bar_outputs,
            bar_inputs,
            cx,
        )
    }
}

impl Component for Counted<Slider> {
    type State = Slip;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        self.inner.ports()
    }

    fn class(&self) -> Class {
        self.inner.class()
    }

    fn stored_energy_j(&self, slip: &Slip) -> Result<f64, ElementError> {
        self.inner.stored_energy_j(slip)
    }

    fn propose(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        tick(&self.proposals);
        self.inner.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

impl LoadLine for Counted<Slider> {
    fn line_ports(&self) -> &[usize] {
        self.inner.line_ports()
    }

    fn propose_on_line(
        &self,
        old: &Slip,
        dt_s: f64,
        inputs: &[f64],
        compliance: &[f64],
        new: &mut Slip,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        tick(&self.proposals);
        self.inner
            .propose_on_line(old, dt_s, inputs, compliance, new, outputs, cx)
    }
}

#[test]
fn a_held_mass_is_not_proposed_or_differentiated_again_where_finding_the_line_already_did() {
    // The mass of `pushed_by`, held by the slider under 4 N sin(2 pi 3.14 t),
    // a push that changes from step to step. Finding the line proposes the
    // mass and takes its tangent at the step's starting inputs, and only the
    // slider's inputs are written between it and the solve's first pass
    // (the released flows, then `release`). So the first pass takes the
    // mass's proposal from finding the line, and the first assembly its
    // column: a step proposes the mass `1 + n` times, where it updates `n`
    // times (before: `2 + n`), and takes its tangent once for finding the
    // line and once for every assembly after the first (before: once for
    // every assembly too).
    //
    // `release` puts the slider's inputs at what the mass and the push give
    // at its held force, and both are affine, so a held step here is taken
    // at its first check, with no update: at an exact fixed point with no
    // slopes, or with the slopes of one assembly made to judge it. Measured
    // before the change, per step (updates, mass proposals, mass tangents,
    // slider asks): (0, 2, 1, 1) and (0, 2, 2, 3). (A held step of the
    // bowed string takes one update: there the string was proposed 3 times
    // and differentiated twice, bead `.12.5`, C.)
    //
    // The slider is the component whose inputs did change there, and it is
    // asked in every pass, the first included, once more for its one unknown
    // column in every assembly (a forward difference), and once more at the
    // accepted inputs where an assembly is followed by no update
    // (differencing wrote its candidate): 1 or 3 asks with no update, `1 +
    // 2 n` with `n`, before and after. Skipped in the first pass, it would
    // be asked 0 or 2 times.
    let dt = 1.0e-3_f64;
    let mass = Counted::new(Mass::new(MASS_KG).unwrap());
    let rubbing = Counted::new(slider());
    let ((mass_proposals, mass_tangents), (slider_asks, _)) = (mass.counts(), rubbing.counts());
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            mass,
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let rubbing = builder
        .add_on_line("slider", rubbing, Slip::stuck())
        .unwrap();
    let push = Waveform::sine(4.0, 2.0 * PI * 0.5, 0.0, 0.0).unwrap();
    let pusher = builder
        .add_with_tangent("push", Force::new(push), SourceAccount::default())
        .unwrap();
    builder
        .common_flow(
            port(&builder, mass, "p"),
            &[port(&builder, rubbing, "a"), port(&builder, pusher, "p")],
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "b"),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "heat"),
            Box::new(|_, _, values| values[0] = 300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    // Per step: updates, then the mass's proposals and tangents and the
    // slider's asks.
    let mut rows = Vec::new();
    with_cx(|cx| {
        for step in 0..200 {
            let before = (
                read(&mass_proposals),
                read(&mass_tangents),
                read(&slider_asks),
            );
            let report = system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
            assert!(
                system.state::<Slip>(rubbing).unwrap().is_stuck(),
                "step {step}"
            );
            rows.push((
                report.iterations as usize,
                read(&mass_proposals) - before.0,
                read(&mass_tangents) - before.1,
                read(&slider_asks) - before.2,
            ));
        }
    });
    // How many steps of each kind.
    let mut kinds: Vec<((usize, usize, usize, usize), usize)> = Vec::new();
    for &row in &rows {
        match kinds.iter_mut().find(|(kind, _)| *kind == row) {
            Some((_, count)) => *count += 1,
            None => kinds.push((row, 1)),
        }
    }
    kinds.sort_unstable();
    eprintln!(
        "held 200 steps under 4 N sin(2 pi 3.14 t); steps by (updates, mass proposals, mass \
         tangents, slider asks): {kinds:?}"
    );
    let mut assembled = 0;
    for (step, &(n, proposals, tangents, asks)) in rows.iter().enumerate() {
        // The slider first: what it is asked says how many assemblies the
        // step made, and is the same before and after.
        let expected = match (n, asks) {
            // An exact fixed point at the first check: no assembly.
            (0, 1) => (1, 1),
            // One assembly to judge the first check, taken with no update.
            (0, 3) => {
                assembled += 1;
                (1, 1)
            }
            // `n` updates, one assembly before each.
            _ if n >= 1 && asks == 1 + 2 * n => {
                assembled += 1;
                (1 + n, n)
            }
            _ => panic!("step {step}: {n} updates, the slider asked {asks} times"),
        };
        assert_eq!(
            (proposals, tangents),
            expected,
            "step {step}: {n} updates and {asks} asks of the slider; the mass proposed \
             {proposals} times and differentiated {tangents} times"
        );
    }
    assert!(assembled >= 1, "no step of 200 assembled");
    assert_eq!(system.steps(), 200);
}
