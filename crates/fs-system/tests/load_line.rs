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
    Damper, ElementError, Force, Gyrator, HeatCapacity, Mass, Motion, Slider, SliderRow, Slip,
    SourceAccount, Spring, Waveform,
};
use fs_system::{AdmissionError, NodeId, SolveConfig, System, SystemBuilder, SystemError};

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

#[test]
fn a_held_mass_is_taken_under_any_push_and_a_planted_flow_is_still_refused() {
    // The flow the kernel publishes for a line input is the released flow
    // less what the efforts take away. For a held slider that is
    // r - c (r / c): zero for some numbers and a unit in the last place of
    // r for others. The junction's audit has to be told that rounding;
    // told nothing, it refused 282 of these 2000 pushes at their first
    // step. Here none is refused, and the mass is at rest to that rounding.
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
    // The sweep covers the pushes whose flow does not come out at zero.
    assert!(uneven >= 100, "{uneven}");

    // The planted negative: what the audit is told is the rounding of that
    // flow and no more. A mass held with half the static limit is brought
    // to rest; then a net force F is planted on the mass's input, where
    // the solve starts, with the solve loosened to 1e-2 so that it takes
    // the step at its first check. The mass then moves by c F with a force
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
        relative_tolerance: 1.0e-2,
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
    // Two terms, each 16 units of eps: the published flow's own (the
    // released flow and what the effort takes away, c * push each) and the
    // one a step taken at its first check is told for the wire's residual
    // (the released flow). Three times 16 eps c push in all.
    let published = 16.0 * f64::EPSILON * (2.0 * c * push);
    let residual = 16.0 * f64::EPSILON * (c * push);
    assert!(
        allowed_m_s >= published && allowed_m_s <= 1.05 * (published + residual),
        "{allowed_m_s} against {published} and {residual}"
    );
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
    // A quarter of it is taken, with no update.
    let report = planted(&mut system, leaving(0.25)).unwrap();
    assert_eq!(report.iterations, 0);
    assert_eq!(system.steps(), 4);
}

/// A mass at rest stuck to a side whose velocity `side(t0, dt)` is
/// prescribed (in the sign the slider's port `b` receives it), the slider's
/// loss collected in a body at 300 K. Returns the system, the mass, the
/// slider and the body.
fn on_a_side(side: impl Fn(f64, f64) -> f64 + 'static) -> (System, NodeId, NodeId, NodeId) {
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
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(BODY_J_K).unwrap(), 300.0)
        .unwrap();
    builder
        .common_flow(port(&builder, mass, "p"), &[port(&builder, rubbing, "a")])
        .unwrap();
    builder
        .boundary(
            port(&builder, rubbing, "b"),
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
        on_a_side(move |t0, _| if t0 > 1.5e-3 { side_m_s } else { 0.0 });
    let at_mass = system.port(rubbing, "a").unwrap();
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
            let true_slip = -velocity - side_now;
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
        .find(|account| account.port == "slider.b")
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
                    }
                    SliderRow::Held => assert_eq!(mean, 0.0),
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

#[test]
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

        // A turning point is where the velocity changes sign; the mass is
        // within a step's travel at a step's speed change of it.
        let mut found = Vec::new();
        let mut before = Motion {
            position_m: 0.0,
            velocity_m_s: v0,
        };
        with_cx(|cx| {
            for _ in 0..(4_000_u32 << halvings) {
                system.step(dt, cx).unwrap();
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
    // The step a block accepts does not depend on the line's matrix: the
    // component settles where the flow on the line is the flow its
    // neighbours give, whatever positive matrix it was told. Only the
    // path of the iteration does. So the matrix is read directly, by a
    // component that keeps what it is handed and returns no effort.
    //
    // A mass answers a force with a mean velocity of dt / (2 m) per
    // newton. With g1 and g2 those of the two masses:
    //   nothing between them:      diag(g1, g2);
    //   a damper d:                [[g1 (1 + g2 d), g1 g2 d], [g1 g2 d, g2 (1 + g1 d)]] / (1 + (g1 + g2) d);
    //   a gyrator G (e_a = G f_b, e_b = -G f_a):
    //                              [[g1, -g1 g2 G], [g1 g2 G, g2]] / (1 + g1 g2 G^2), not symmetric;
    //   the second side prescribed: [[g1, 0], [0, 0]].
    let dt = 1.0e-3_f64;
    let (g1, g2) = (0.5 * dt / 2.0, 0.5 * dt / 0.5);
    let (d, ratio) = (700.0_f64, 900.0_f64);
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
}

#[test]
fn a_refused_step_leaves_nothing_of_its_line_behind() {
    // A mass thrown along the ground on a spring that is added without its
    // tangent, so the kernel differences it. One system is asked for a
    // step with no iteration allowed: the line has been found, the drivers
    // carry its matrix and the inputs hold released flows when the solve
    // gives up. The step that follows, and a refusal half way, must be the
    // steps of a system that was never refused, to the bit.
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
    assert_eq!(refusals, 3);
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
    for (inertia_kg, reason) in [(-MASS_KG, "symmetric part is negative"), (0.0, "not finite")] {
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
