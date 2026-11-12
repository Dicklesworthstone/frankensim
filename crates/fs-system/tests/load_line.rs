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
    ElementError, Force, HeatCapacity, Mass, Motion, Slider, SliderRow, Slip, SourceAccount,
    Spring, Waveform,
};
use fs_system::{AdmissionError, NodeId, System, SystemBuilder, SystemError};

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
        .add_on_line(
            "slider",
            slider(),
            Slip {
                slip_m_s: -velocity_m_s,
            },
        )
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
                let (force, mean, row) = law.on_line(start, released, c);
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
                let ends = 2.0 * mean - start;
                match row {
                    SliderRow::Held => assert!(start == 0.0 && mean == 0.0),
                    SliderRow::BrokeAway => assert!(start == 0.0 && ends * released > 0.0),
                    SliderRow::Arrested => assert_eq!(ends, 0.0),
                    SliderRow::Slid => assert!(ends * start > 0.0),
                    SliderRow::Released => assert!(force == 0.0 && ends * start <= 0.0),
                    SliderRow::Turned => assert!(mean == 0.0 && ends == -start),
                    SliderRow::Reversed => assert!(ends * start < 0.0),
                }
                if scan % 3 != 0 {
                    continue;
                }
                // On a line over both ends whose matrix gives this `c`.
                let matrix = [0.7 * c, -0.1 * c, 0.05 * c, 0.25 * c];
                assert!((matrix[0] - matrix[1] - matrix[2] + matrix[3] - c).abs() <= 1.0e-12 * c);
                let inputs = [0.3 + released, 0.3, 310.0];
                let state = Slip { slip_m_s: start };
                let stale = Slip { slip_m_s: 7.0 };
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
        assert!(rows.get(&row).copied().unwrap_or(0) >= 20, "{row:?}: {rows:?}");
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
                let slip = system.state::<Slip>(rubbing).unwrap().slip_m_s;
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
    assert_eq!(held.state::<Slip>(rubbing).unwrap().slip_m_s, 0.0);
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
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s, -v);
    with_cx(|cx| system.step(dt, cx)).unwrap();
    assert_eq!(*system.state::<f64>(body).unwrap(), 300.0);
    assert_eq!(system.last_output(system.port(rubbing, "a").unwrap()), [0.0]);
    let motion = *system.state::<Motion>(mass).unwrap();
    assert!((motion.velocity_m_s + 0.5 * v).abs() <= 1.0e-12 * v, "{motion:?}");
    let slip = system.state::<Slip>(rubbing).unwrap().slip_m_s;
    assert!((slip - 0.5 * v).abs() <= 1.0e-12 * v, "{slip}");
    // From there friction acts again, against the new direction: the
    // push, 3 N, is inside what the slider holds, and the body is brought
    // to rest and kept there.
    with_cx(|cx| {
        for _ in 0..3 {
            system.step(dt, cx).unwrap();
        }
    });
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s, 0.0);
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
        .add_on_line("slider", slider(), Slip { slip_m_s: -v0 })
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
            if joined_s.is_nan() && system.state::<Slip>(rubbing).unwrap().slip_m_s == 0.0 {
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
            let now = system.state::<Slip>(rubbing).unwrap().slip_m_s == 0.0;
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
        .add_on_line("slider", slider(), Slip { slip_m_s: slip })
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
    assert_eq!(system.state::<Slip>(rubbing).unwrap().slip_m_s, slip);
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
                straight.state::<Slip>(rubbing).unwrap().slip_m_s.to_bits(),
                detoured.state::<Slip>(rubbing).unwrap().slip_m_s.to_bits(),
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
    assert_eq!(straight.state::<Slip>(rubbing).unwrap().slip_m_s, 0.0);
    assert!(straight.state::<Motion>(mass).unwrap().position_m > 0.2);
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

/// A slider that names its heat port as its line.
struct Misplaced(Slider);

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
        &[0, 2]
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
    // A line on a port that carries no power is refused when the component
    // is added, and claims no name.
    let mut builder = SystemBuilder::new();
    let refusal = builder
        .add_on_line("misplaced", Misplaced(slider()), Slip::default())
        .unwrap_err();
    assert_eq!(
        refusal,
        AdmissionError::LinePorts {
            node: "misplaced".to_owned()
        }
    );
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
    // line of negative compliance: the first step is refused, the slider
    // named, and nothing is published.
    let mut builder = SystemBuilder::new();
    let wrong = builder.add("wrong", Antimass::new(-MASS_KG), 0.0).unwrap();
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
        matches!(&refusal, SystemError::LoadLine { node, .. } if node == "slider"),
        "{refusal}"
    );
    assert_eq!(system.steps(), 0);

    // A slider whose two ends are both prescribed is in no loop: there is
    // no line to find and it is stepped by its mandatory tier, sliding at
    // its kinetic force.
    let mut builder = SystemBuilder::new();
    let rubbing = builder
        .add_on_line("slider", slider(), Slip { slip_m_s: 0.2 })
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
