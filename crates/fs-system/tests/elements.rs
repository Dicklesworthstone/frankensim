//! The library elements: each passes the conformance kit over a sweep of
//! states and inputs, and devices composed from them reproduce closed-form
//! solutions of the discrete equations they are meant to integrate.

// Exact float comparisons here assert bitwise identities of the wiring.
#![allow(clippy::float_cmp)]

use std::f64::consts::PI;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{
    Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};
use fs_system::elements::{
    Conductance, Damper, ElementError, Force, GasContents, GasVolume, Gyrator, HeatCapacity,
    HeatNode, Heater, Mass, Motion, Restriction, SourceAccount, Spring, Waveform,
};
use fs_system::{
    AdmissionError, NodeId, RouteTemperature, SolveConfig, System, SystemBuilder, SystemError,
    heat_sink,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 11,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// Deterministic numbers in `[lo, hi)` (splitmix64), so a sweep is the same
/// on every run.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        lo + (hi - lo) * ((z >> 11) as f64 / (1_u64 << 53) as f64)
    }
}

/// Run the mandatory and tangent checks at one operating point and return
/// the outputs of the step.
fn conforms<C>(element: &C, state: &C::State, stale: &C::State, dt_s: f64, inputs: &[f64]) -> Vec<f64>
where
    C: PortTangent,
    C::State: PartialEq + core::fmt::Debug,
{
    let case = Case {
        state,
        stale,
        dt_s,
        inputs,
    };
    with_cx(|cx| {
        let observed = check_component(element, &case, cx)
            .unwrap_or_else(|failure| panic!("conformance at {inputs:?}: {failure}"));
        let outputs = match observed {
            Observed::Stepped { outputs, .. } => outputs,
            Observed::Shortened { dt_s } => panic!("shortened to {dt_s} at {inputs:?}"),
        };
        // The kit's central difference steps by h = max(1, |inputs|) eps^(1/3).
        // An output that carries a large offset (a temperature) loses
        // about 2 eps |output| / (2 h) of its slope to roundoff there.
        let largest = |values: &[f64]| values.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        let h = largest(inputs) * f64::EPSILON.cbrt();
        let absolute = 1.0e-12_f64.max(8.0 * f64::EPSILON * largest(&outputs) / h);
        check_port_tangent(element, &case, 1.0e-6, absolute, cx)
            .unwrap_or_else(|failure| panic!("tangent at {inputs:?}: {failure}"));
        outputs
    })
}

#[test]
fn every_element_conforms_across_a_sweep_of_states_and_inputs() {
    let mut n = Numbers(20_261_110);
    let mass = Mass::new(0.75).unwrap();
    let linear = Spring::new(350.0, 0.0).unwrap();
    let hardening = Spring::new(350.0, 2.0e6).unwrap();
    let softening = Spring::new(350.0, -1.0e5).unwrap();
    let damper = Damper::new(2.5).unwrap();
    let push = Force::new(Waveform::sine(3.0, 7.0, 0.4, -0.5).unwrap());
    let body = HeatCapacity::new(40.0).unwrap();
    let wall = Conductance::new(0.8).unwrap();
    let heater = Heater::new(Waveform::table(vec![(0.0, 0.0), (1.0, 12.0), (2.0, 4.0)]).unwrap()).unwrap();
    for _ in 0..200 {
        let dt = n.next(1.0e-5, 2.0e-2);
        let (f_a, f_b) = (n.next(-3.0, 3.0), n.next(-3.0, 3.0));
        let (t_a, t_b) = (n.next(250.0, 900.0), n.next(250.0, 900.0));
        let motion = Motion {
            position_m: n.next(-0.1, 0.1),
            velocity_m_s: n.next(-4.0, 4.0),
        };
        let stale = Motion {
            position_m: 9.0,
            velocity_m_s: -9.0,
        };
        let force = n.next(-50.0, 50.0);
        let out = conforms(&mass, &motion, &stale, dt, &[force]);
        // Independent: the mean velocity of a constant force over the step.
        let expected = motion.velocity_m_s + 0.5 * dt * force / 0.75;
        assert!((out[0] - expected).abs() <= 1.0e-14 * (1.0 + expected.abs()));

        let extension = n.next(-0.02, 0.02);
        for spring in [&linear, &hardening, &softening] {
            let out = conforms(spring, &extension, &0.5, dt, &[f_a, f_b]);
            assert_eq!(out[0], -out[1], "a spring pulls both ends equally");
        }
        let out = conforms(&linear, &extension, &0.5, dt, &[f_a, f_b]);
        let mean_extension = extension + 0.5 * dt * (f_a - f_b);
        assert!((out[0] - 350.0 * mean_extension).abs() <= 1.0e-12 * (1.0 + out[0].abs()));

        let out = conforms(&damper, &(), &(), dt, &[f_a, f_b, t_a]);
        assert_eq!(out[0], -out[1]);
        assert_eq!(out[0], 2.5 * (f_a - f_b));
        assert_eq!(out[2], -(out[0] * (f_a - f_b)), "all absorbed power leaves as heat");

        let mut account = SourceAccount::starting_at(n.next(0.0, 3.0));
        account.energy_in_j = n.next(-10.0, 10.0);
        let mut stale_account = SourceAccount::starting_at(-1.0);
        stale_account.energy_in_j = 99.0;
        conforms(&push, &account, &stale_account, dt, &[f_a]);
        let out = conforms(&heater, &account, &stale_account, dt, &[t_a]);
        assert!(out[0] <= 0.0, "a heater only gives heat away");

        let temperature = n.next(250.0, 900.0);
        let heat_rate = n.next(-500.0, 500.0);
        conforms(&body, &temperature, &1.0, dt, &[heat_rate]);

        let out = conforms(&wall, &(), &(), dt, &[t_a, t_b]);
        assert_eq!(out[0], -out[1], "a conductance stores nothing");
        assert_eq!(out[0], 0.8 * (t_a - t_b));
    }
}

#[test]
fn elements_refuse_what_they_cannot_represent() {
    assert!(matches!(
        Mass::new(0.0),
        Err(ElementError::Parameter { name: "mass_kg", .. })
    ));
    assert!(Spring::new(f64::NAN, 0.0).is_err());
    assert!(Spring::new(1.0, f64::INFINITY).is_err());
    assert!(Damper::new(-1.0).is_err());
    assert!(HeatCapacity::new(0.0).is_err());
    assert!(Conductance::new(f64::NAN).is_err());
    assert!(Waveform::constant(f64::NAN).is_err());
    assert!(Waveform::sine(1.0, -1.0, 0.0, 0.0).is_err());
    assert!(matches!(Waveform::table(vec![]), Err(ElementError::Table { .. })));
    assert!(Waveform::table(vec![(0.0, 1.0), (0.0, 2.0)]).is_err());
    assert!(Waveform::table(vec![(0.0, 1.0), (f64::NAN, 2.0)]).is_err());
    // A heater whose power would turn negative is refused when it is built.
    assert!(Heater::new(Waveform::sine(2.0, 1.0, 0.0, 1.0).unwrap()).is_err());
    assert!(Heater::new(Waveform::sine(2.0, 1.0, 0.0, 2.0).unwrap()).is_ok());
    assert!(Heater::new(Waveform::table(vec![(0.0, 1.0), (1.0, -0.1)]).unwrap()).is_err());

    with_cx(|cx| {
        // A temperature that is not positive is refused, not divided by.
        let damper = Damper::new(1.0).unwrap();
        let mut outputs = [0.0; 3];
        let refused = damper.propose(&(), 1.0e-3, &[1.0, 0.0, 0.0], &mut (), &mut outputs, cx);
        assert_eq!(refused, Err(ElementError::Temperature { kelvin: 0.0 }));
        let wall = Conductance::new(1.0).unwrap();
        let mut outputs = [0.0; 2];
        assert!(wall.propose(&(), 1.0e-3, &[300.0, f64::NAN], &mut (), &mut outputs, cx).is_err());
        // A step that would cool a body below absolute zero asks to be shorter.
        let body = HeatCapacity::new(1.0).unwrap();
        let (mut next, mut outputs) = (0.0, [0.0; 1]);
        let answer = body.propose(&300.0, 1.0, &[-400.0], &mut next, &mut outputs, cx);
        assert_eq!(answer, Ok(Proposal::Shorten { dt_s: 0.5 }));
        // A force that is not finite is refused.
        let mass = Mass::new(1.0).unwrap();
        let rest = Motion {
            position_m: 0.0,
            velocity_m_s: 0.0,
        };
        let mut moved = rest;
        let refused = mass.propose(&rest, 1.0e-3, &[f64::INFINITY], &mut moved, &mut [0.0], cx);
        assert!(matches!(refused, Err(ElementError::Input { .. })));
    });
}

/// Composite Simpson rule, for an independent mean.
fn simpson_mean(f: impl Fn(f64) -> f64, t0: f64, dt: f64) -> f64 {
    let intervals = 20_000;
    let h = dt / f64::from(intervals);
    let mut sum = f(t0) + f(t0 + dt);
    for i in 1..intervals {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * f(t0 + h * f64::from(i));
    }
    sum * h / 3.0 / dt
}

#[test]
fn waveform_means_are_exact_averages_over_the_step() {
    let sine = Waveform::sine(2.5, 3.0, 0.7, 1.25).unwrap();
    let table = Waveform::table(vec![(0.5, 2.0), (1.0, -1.0), (1.5, 3.0), (4.0, 3.5)]).unwrap();
    for (t0, dt) in [(0.0, 0.1), (0.3, 1.7), (1.2, 1.0e-6), (0.9, 3.0), (-1.0, 0.25), (5.0, 2.0)] {
        for wave in [&sine, &table] {
            let reference = simpson_mean(|t| wave.value_at(t), t0, dt);
            let mean = wave.mean(t0, dt);
            assert!(
                (mean - reference).abs() <= 1.0e-9 * (1.0 + reference.abs()),
                "{wave:?} over [{t0}, {}]: {mean} vs {reference}",
                t0 + dt
            );
            // Means of adjacent intervals add up to the mean of their union.
            let (first, second) = (wave.mean(t0, 0.3 * dt), wave.mean(t0 + 0.3 * dt, 0.7 * dt));
            let joined = 0.3 * first + 0.7 * second;
            assert!((joined - mean).abs() <= 1.0e-12 * (1.0 + mean.abs()));
        }
    }
    // Hand values: before the first knot and after the last the end values
    // hold, between knots the table is linear.
    assert_eq!(table.value_at(0.0), 2.0);
    assert_eq!(table.value_at(10.0), 3.5);
    assert_eq!(table.value_at(0.75), 0.5);
    assert_eq!(table.mean(0.5, 0.5), 0.5);
    assert_eq!(table.mean(-2.0, 1.0), 2.0);
    // A sine of zero frequency is the constant it starts at.
    let still = Waveform::sine(2.0, 0.0, 0.5, 1.0).unwrap();
    assert_eq!(still.mean(3.0, 0.5), 1.0 + 2.0 * 0.5_f64.sin());
    assert_eq!(still.lowest(), 1.0 + 2.0 * 0.5_f64.sin());
    assert_eq!(table.lowest(), -1.0);
    assert_eq!(sine.lowest(), 1.25 - 2.5);
}

fn step(system: &mut System, dt_s: f64, steps: u32) {
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt_s, cx).expect("step");
        }
    });
}

fn anchored(builder: &mut SystemBuilder, node: NodeId, port: &str) {
    let port = builder.port(node, port).unwrap();
    builder
        .boundary(port, Box::new(|_, _, values| values[0] = 0.0))
        .unwrap();
}

#[test]
fn an_anchored_mass_and_spring_rotate_by_the_trapezoidal_angle_and_keep_their_energy() {
    let (m, k, x0, v0) = (0.4_f64, 900.0_f64, 0.012_f64, -0.35_f64);
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(m).unwrap(),
            Motion {
                position_m: x0,
                velocity_m_s: v0,
            },
        )
        .unwrap();
    // The spring runs from the mass (a) to a fixed point at the origin (b),
    // so its extension is minus the mass's position.
    let spring = builder
        .add_with_tangent("spring", Spring::new(k, 0.0).unwrap(), -x0)
        .unwrap();
    let (p, a) = (builder.port(mass, "p").unwrap(), builder.port(spring, "a").unwrap());
    builder.common_flow(p, &[a]).unwrap();
    anchored(&mut builder, spring, "b");
    let mut system = builder.seal().unwrap();

    let (dt, steps) = (1.0e-3, 5_000_u32);
    step(&mut system, dt, steps);

    // The implicit midpoint rule maps (x, v / w) by a rotation through
    // 2 atan(w dt / 2) per step.
    let omega = (k / m).sqrt();
    let angle = f64::from(steps) * 2.0 * (0.5 * omega * dt).atan();
    let x = x0 * angle.cos() + v0 / omega * angle.sin();
    let v = -x0 * omega * angle.sin() + v0 * angle.cos();
    let motion = *system.state::<Motion>(mass).unwrap();
    assert!((motion.position_m - x).abs() <= 1.0e-10, "{} vs {x}", motion.position_m);
    assert!((motion.velocity_m_s - v).abs() <= 1.0e-8, "{} vs {v}", motion.velocity_m_s);
    // The spring integrates the same flow with the opposite sign; the two
    // agree to the interface tolerance of each step.
    let extension = *system.state::<f64>(spring).unwrap();
    assert!((extension + motion.position_m).abs() <= 1.0e-13, "{extension} vs {}", motion.position_m);

    let energy = |x: f64, v: f64| 0.5 * m * v * v + 0.5 * k * x * x;
    let (start, end) = (energy(x0, v0), energy(motion.position_m, motion.velocity_m_s));
    assert!((end - start).abs() <= 1.0e-12 * start, "{start} -> {end}");
    let ledger = system.ledger();
    assert!(ledger.stored_change_j().abs() <= 1.0e-12 * start);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-12 * start);
    assert_eq!(ledger.entropy_generated_j_per_k(), 0.0);
}

#[test]
fn two_masses_on_a_spring_keep_their_momentum_and_their_centre_of_mass_moves_uniformly() {
    let (m1, m2) = (0.3_f64, 1.1_f64);
    let (v1, v2) = (1.5_f64, -0.2_f64);
    let mut builder = SystemBuilder::new();
    let left = builder
        .add_with_tangent(
            "left",
            Mass::new(m1).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: v1,
            },
        )
        .unwrap();
    let right = builder
        .add_with_tangent(
            "right",
            Mass::new(m2).unwrap(),
            Motion {
                position_m: 0.1,
                velocity_m_s: v2,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(2_000.0, 5.0e5).unwrap(), 0.0)
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (a, b) = (port(&builder, spring, "a"), port(&builder, spring, "b"));
    let (pl, pr) = (port(&builder, left, "p"), port(&builder, right, "p"));
    builder.common_flow(pl, &[a]).unwrap();
    builder.common_flow(pr, &[b]).unwrap();
    let mut system = builder.seal().unwrap();
    let (dt, steps) = (2.0e-4, 4_000_u32);
    step(&mut system, dt, steps);

    let (l, r) = (
        *system.state::<Motion>(left).unwrap(),
        *system.state::<Motion>(right).unwrap(),
    );
    let momentum = m1 * v1 + m2 * v2;
    let now = m1 * l.velocity_m_s + m2 * r.velocity_m_s;
    assert!((now - momentum).abs() <= 1.0e-12 * (m1 * v1.abs() + m2 * v2.abs()));
    let centre = (m1 * l.position_m + m2 * r.position_m) / (m1 + m2);
    let expected = (m2 * 0.1) / (m1 + m2) + momentum / (m1 + m2) * dt * f64::from(steps);
    assert!((centre - expected).abs() <= 1.0e-12, "{centre} vs {expected}");
    // The extension is the separation less the initial separation.
    let extension = *system.state::<f64>(spring).unwrap();
    let separation = r.position_m - l.position_m - 0.1;
    assert!((extension - separation).abs() <= 1.0e-12);
    // The hardening spring was actually exercised.
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-12);
}

#[test]
fn two_bodies_through_a_conductance_relax_by_the_discrete_factor_and_generate_the_entropy_of_mixing() {
    let (c1, c2, g) = (50.0_f64, 120.0_f64, 3.0_f64);
    let (t1, t2) = (400.0_f64, 300.0_f64);
    let mut builder = SystemBuilder::new();
    let hot = builder
        .add_with_tangent("hot", HeatCapacity::new(c1).unwrap(), t1)
        .unwrap();
    let cold = builder
        .add_with_tangent("cold", HeatCapacity::new(c2).unwrap(), t2)
        .unwrap();
    let wall = builder
        .add_with_tangent("wall", Conductance::new(g).unwrap(), ())
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (a, b) = (port(&builder, wall, "a"), port(&builder, wall, "b"));
    let (h, c) = (port(&builder, hot, "heat"), port(&builder, cold, "heat"));
    builder.common_effort(h, &[a]).unwrap();
    builder.common_effort(c, &[b]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(a, &[t1]));
    assert!(system.guess(b, &[t2]));

    let (dt, steps) = (0.05, 2_000_u32);
    step(&mut system, dt, steps);

    let factor = {
        let half = 0.5 * g * dt * (1.0 / c1 + 1.0 / c2);
        (1.0 - half) / (1.0 + half)
    };
    let difference = (t1 - t2) * factor.powi(steps as i32);
    let (h_k, c_k) = (
        *system.state::<f64>(hot).unwrap(),
        *system.state::<f64>(cold).unwrap(),
    );
    assert!(((h_k - c_k) - difference).abs() <= 1.0e-9 * (t1 - t2), "{} vs {difference}", h_k - c_k);
    let energy = c1 * t1 + c2 * t2;
    assert!((c1 * h_k + c2 * c_k - energy).abs() <= 1.0e-12 * energy);

    // Entropy: the ledger's sum of G dt (Ta - Tb)^2 / (Ta Tb) at the step
    // means agrees with the bodies' exact entropy change up to a third-order
    // term per step, which is far below this bound here.
    let exact = c1 * (h_k / t1).ln() + c2 * (c_k / t2).ln();
    let generated = system.ledger().entropy_generated_j_per_k();
    assert!(exact > 0.0);
    assert!((generated - exact).abs() <= 1.0e-6 * exact, "{generated} vs {exact}");
}

#[test]
fn a_heat_node_keeps_its_balance_and_refuses_what_leaves_it_no_temperature() {
    let conductances = [0.8, 2.5, 0.3];
    let total: f64 = conductances.iter().sum();
    let node = HeatNode::new(&conductances).unwrap();
    let ports: Vec<&str> = node.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(ports, ["source", "body0", "body1", "body2"]);
    let mut n = Numbers(0x4EA7_0001);
    for _ in 0..200 {
        let dt = n.next(1.0e-5, 2.0e-2);
        let source = n.next(-40.0, 120.0);
        let bodies = [n.next(250.0, 900.0), n.next(250.0, 900.0), n.next(250.0, 900.0)];
        let inputs = [source, bodies[0], bodies[1], bodies[2]];
        let out = conforms(&node, &(), &(), dt, &inputs);
        // Independent: the node sits where what the source puts in leaves
        // through the conductances, P = sum G_k (T_n - T_k).
        let weighted: f64 = conductances.iter().zip(bodies).map(|(g, t)| g * t).sum();
        let kelvin = (source + weighted) / total;
        assert!((out[0] - kelvin).abs() <= 1.0e-13 * kelvin, "{} vs {kelvin}", out[0]);
        let mut entering = 0.0;
        for ((rate, g), t) in out[1..].iter().zip(conductances).zip(bodies) {
            assert!((rate - g * (t - kelvin)).abs() <= 1.0e-11, "{rate}");
            entering += rate;
        }
        assert!((entering + source).abs() <= 1.0e-10, "{entering} vs {source}");
        // It stores nothing, so the entropy it generates is what leaves
        // through its ports less what enters, and each conductance's share
        // is positive.
        let case = Case {
            state: &(),
            stale: &(),
            dt_s: dt,
            inputs: &inputs,
        };
        let observed = with_cx(|cx| check_component(&node, &case, cx)).unwrap();
        let Observed::Stepped { receipt, audit, .. } = observed else {
            panic!("the node takes the step");
        };
        let generated = receipt.entropy_generated_j_per_k;
        let by_branch: f64 = conductances
            .iter()
            .zip(bodies)
            .map(|(g, t)| dt * g * (t - kelvin) * (t - kelvin) / (t * kelvin))
            .sum();
        assert!(generated > 0.0);
        assert!((generated - by_branch).abs() <= 1.0e-12 * by_branch);
        assert!((generated + audit.port_entropy_in_j_per_k).abs() <= 1.0e-12 * by_branch);
    }

    // Bodies within a nanokelvin of one temperature, and exactly at it: the
    // rates all but vanish while their rounding stays that of the terms
    // G T they are differences of. (The node's first energy tolerance was
    // sized by the rates and refused such a step in a graph.)
    for spread in [1.0e-9, 1.0e-12, 0.0] {
        let out = conforms(&node, &(), &(), 0.5, &[0.0, 329.4, 329.4 + spread, 329.4 - spread]);
        assert!(out[1..].iter().all(|rate| rate.abs() <= 1.0e-8), "{out:?}");
    }

    assert!(HeatNode::new(&[]).is_err());
    assert!(HeatNode::new(&[1.0, 0.0]).is_err());
    assert!(HeatNode::new(&[1.0, -2.0]).is_err());
    assert!(HeatNode::new(&[f64::NAN]).is_err());
    assert!(HeatNode::new(&[f64::MAX, f64::MAX]).is_err());
    with_cx(|cx| {
        let mut outputs = [7.0; 4];
        // A source that takes more than the conductances bring at any
        // temperature above zero; a body at no temperature; a value that is
        // not a number. Nothing is written.
        for inputs in [
            [-total * 300.0, 300.0, 300.0, 300.0],
            [-total * 300.0 - 1.0, 300.0, 300.0, 300.0],
            [0.0, 300.0, 0.0, 300.0],
            [0.0, 300.0, 300.0, f64::NAN],
            [f64::INFINITY, 300.0, 300.0, 300.0],
        ] {
            assert!(node.propose(&(), 1.0e-3, &inputs, &mut (), &mut outputs, cx).is_err(), "{inputs:?}");
            assert_eq!(outputs, [7.0; 4]);
        }
    });
}

/// The pad's heat capacity [J/K], the plane's conductances to the pad and
/// to the held side [W/K], the held temperature [K], the power generated
/// on the plane [W] and the step [s] of [`held_plane`].
const HELD_PAD_J_K: f64 = 50.0;
const HELD_PLANE_W_K: [f64; 2] = [2.0, 6.0];
const HELD_ROOM_K: f64 = 300.0;
const HELD_POWER_W: f64 = 30.0;
const HELD_DT_S: f64 = 0.5;

/// A pad at 400 K on a plane whose other side is held at a fixed
/// temperature, with a heater on the plane. The pad and the plane meet in
/// the first connection made, a junction of two ports. Returns the system,
/// the pad and the pad's heat port.
fn held_plane() -> (System, NodeId, fs_system::PortRef) {
    let room_k = HELD_ROOM_K;
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let pad = builder
        .add_with_tangent("pad", HeatCapacity::new(HELD_PAD_J_K).unwrap(), 400.0)
        .unwrap();
    let plane = builder
        .add_with_tangent("plane", HeatNode::new(&HELD_PLANE_W_K).unwrap(), ())
        .unwrap();
    let rubbing = heater(&mut builder, "rubbing", HELD_POWER_W);
    let (at_pad, at_room) = (port(&builder, plane, "body0"), port(&builder, plane, "body1"));
    let at_source = port(&builder, rubbing, "heat");
    let pad_heat = port(&builder, pad, "heat");
    builder.common_effort(pad_heat, &[at_pad]).unwrap();
    builder
        .boundary(at_room, Box::new(move |_, _, values| values[0] = room_k))
        .unwrap();
    builder
        .common_effort(port(&builder, plane, "source"), &[at_source])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(at_pad, &[400.0]));
    assert!(system.guess(at_source, &[325.0]));
    (system, pad, pad_heat)
}

#[test]
fn a_residual_planted_on_a_junction_at_rest_is_still_refused_at_a_first_check() {
    // The counterpart of the test below. A step taken at its first check
    // tells the junction's audit the rounding of the terms its heat is a
    // difference of, and no more than that. The held plane is brought to
    // rest; then a heat rate is planted on the pad's input, where the
    // solve starts, and the solve is loosened to 1e-6 (of the 30 W the
    // block carries: 3e-5 W) so that it takes the step at its first check
    // all the same. The audit is not told the solve's relative tolerance.
    let dt = HELD_DT_S;
    let (mut system, pad, pad_heat) = held_plane();
    with_cx(|cx| {
        for _ in 0..2_500 {
            system.step(dt, cx).unwrap();
        }
    });
    let rest_k = HELD_ROOM_K + HELD_POWER_W / HELD_PLANE_W_K[1];
    assert!((*system.state::<f64>(pad).unwrap() - rest_k).abs() <= 1.0e-11 * rest_k);
    let at_rest = system.checkpoint().unwrap();
    let settled_w = system.last_input(pad_heat)[0];
    // What the solve holds on the junction at rest is inside its own
    // tolerance for a heat rate of this block, 1e-12 of the 30 W it
    // carries.
    assert!(settled_w.abs() <= 1.0e-10, "{settled_w}");
    let loose = SolveConfig {
        relative_tolerance: 1.0e-6,
        ..SolveConfig::default()
    };
    let planted = |system: &mut System, extra_w: f64| {
        system.restore(&at_rest).unwrap();
        system.set_config(loose);
        assert!(system.guess(pad_heat, &[settled_w + extra_w]));
        with_cx(|cx| system.step(dt, cx))
    };

    // The audit's tolerance for this junction at rest, read from a refusal
    // of a microwatt: the relative part of it is 1e-9 of half a microjoule.
    let Err(SystemError::Interface {
        connection: 0,
        tolerance_j,
        ..
    }) = planted(&mut system, 1.0e-6)
    else {
        panic!("a planted microwatt is refused by the audit of the first connection");
    };
    let allowed_w = tolerance_j / dt;
    // It is the rounding of the terms and nothing looser: 16 units of eps
    // times the sizes the pad's heat depends on through the plane, 1.5 W/K
    // of the pad's temperature and a quarter of the source's power.
    let reach = HELD_PLANE_W_K[0] / (HELD_PLANE_W_K[0] + HELD_PLANE_W_K[1]);
    let terms_w = HELD_PLANE_W_K[1] * reach * rest_k + reach * HELD_POWER_W;
    let estimate_w = 16.0 * f64::EPSILON * terms_w;
    assert!(
        allowed_w >= estimate_w && allowed_w <= 1.5 * estimate_w,
        "{allowed_w} against {estimate_w}"
    );

    // A hundred times what is allowed: the solve takes the step at its
    // first check (the refusal comes from the audit, not from the solve),
    // and the audit refuses it by that factor.
    let refusal = planted(&mut system, 100.0 * allowed_w).expect_err("refused");
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
    assert_eq!(system.steps(), 2_500);

    // A quarter of what is allowed: taken, with no update.
    let report = planted(&mut system, 0.25 * allowed_w).unwrap();
    assert_eq!(report.iterations, 0);
    assert_eq!(system.steps(), 2_501);
}

#[test]
fn a_body_that_has_come_to_its_plane_s_temperature_keeps_stepping() {
    // A pad of 50 J/K at 400 K on a plane of 2 and 6 W/K whose other side
    // is held at 300 K, with 30 W generated on the plane. The pad settles
    // where it exchanges nothing, at 300 K + P / G1 = 305 K, by the
    // midpoint rule's factor for the rate 2 * 6 / (2 + 6) / 50 per second.
    //
    // At rest the heat between pad and plane is a rounding-sized difference
    // of terms G T, while the energy that crosses their junction goes to
    // zero: the junction's audit has nothing left to be relative to and
    // must be told the rounding of those terms, in a step taken at its
    // first check as in one that iterated.
    let (capacity, conductances, room_k, power_w) =
        (HELD_PAD_J_K, HELD_PLANE_W_K, HELD_ROOM_K, HELD_POWER_W);
    let dt = HELD_DT_S;
    let (mut system, pad, _) = held_plane();

    let rest_k = room_k + power_w / conductances[1];
    let rate = conductances[0] * conductances[1] / (conductances[0] + conductances[1]) / capacity;
    let factor = (1.0 - 0.5 * dt * rate) / (1.0 + 0.5 * dt * rate);
    let (mut expected_k, mut taken_at_first_check) = (400.0_f64, 0_u32);
    with_cx(|cx| {
        for step in 0..4_000_u32 {
            let report = system
                .step(dt, cx)
                .unwrap_or_else(|error| panic!("step {step}: {error}"));
            taken_at_first_check += u32::from(report.iterations == 0);
            expected_k = rest_k + factor * (expected_k - rest_k);
            let pad_k = *system.state::<f64>(pad).unwrap();
            assert!(
                (pad_k - expected_k).abs() <= 1.0e-11 * rest_k,
                "step {step}: {pad_k} against {expected_k}"
            );
        }
    });
    // It reached rest, and went on stepping there without an update.
    let pad_k = *system.state::<f64>(pad).unwrap();
    assert!((pad_k - rest_k).abs() <= 1.0e-11 * rest_k, "{pad_k}");
    assert!(taken_at_first_check >= 1_000, "{taken_at_first_check}");
}

#[test]
fn two_bodies_through_a_heat_node_are_joined_by_its_conductances_in_series() {
    // The node with 2 and 6 W/K and nothing on its source is one
    // conductance of 2 * 6 / (2 + 6) = 1.5 W/K between the bodies: the same
    // heat for the same two temperatures, so the two graphs keep step to
    // the rounding of that arithmetic.
    let (c1, c2, t1, t2) = (50.0_f64, 120.0_f64, 400.0_f64, 300.0_f64);
    let bodies = |builder: &mut SystemBuilder| {
        let hot = builder
            .add_with_tangent("hot", HeatCapacity::new(c1).unwrap(), t1)
            .unwrap();
        let cold = builder
            .add_with_tangent("cold", HeatCapacity::new(c2).unwrap(), t2)
            .unwrap();
        (hot, cold)
    };
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    // With `power_w`, a heater on the node's source; without, nothing.
    let through_node = |power_w: Option<f64>| {
        let mut builder = SystemBuilder::new();
        let (hot, cold) = bodies(&mut builder);
        let node = builder
            .add_with_tangent("node", HeatNode::new(&[2.0, 6.0]).unwrap(), ())
            .unwrap();
        let (a, b) = (port(&builder, node, "body0"), port(&builder, node, "body1"));
        builder.common_effort(port(&builder, hot, "heat"), &[a]).unwrap();
        builder.common_effort(port(&builder, cold, "heat"), &[b]).unwrap();
        let source = port(&builder, node, "source");
        let lamp = power_w.map(|power_w| {
            let lamp = heater(&mut builder, "lamp", power_w);
            let heat = port(&builder, lamp, "heat");
            builder.common_effort(source, &[heat]).unwrap();
            heat
        });
        if lamp.is_none() {
            builder
                .boundary(source, Box::new(|_, _, values| values[0] = 0.0))
                .unwrap();
        }
        let mut system = builder.seal().unwrap();
        // Temperatures received inside the loop need a start.
        assert!(system.guess(a, &[t1]));
        assert!(system.guess(b, &[t2]));
        if let Some(heat) = lamp {
            assert!(system.guess(heat, &[t1]));
        }
        (system, hot, cold, source)
    };
    let (mut split, hot, cold, _) = through_node(None);
    let mut builder = SystemBuilder::new();
    let (hot_one, cold_one) = bodies(&mut builder);
    let wall = builder
        .add_with_tangent("wall", Conductance::new(1.5).unwrap(), ())
        .unwrap();
    let (a, b) = (port(&builder, wall, "a"), port(&builder, wall, "b"));
    builder.common_effort(port(&builder, hot_one, "heat"), &[a]).unwrap();
    builder.common_effort(port(&builder, cold_one, "heat"), &[b]).unwrap();
    let mut whole = builder.seal().unwrap();
    assert!(whole.guess(a, &[t1]));
    assert!(whole.guess(b, &[t2]));
    let (dt, steps) = (0.5, 400_u32);
    let mut worst = 0.0_f64;
    with_cx(|cx| {
        for _ in 0..steps {
            split.step(dt, cx).unwrap();
            whole.step(dt, cx).unwrap();
            for (through, direct) in [(hot, hot_one), (cold, cold_one)] {
                let gap = split.state::<f64>(through).unwrap() - whole.state::<f64>(direct).unwrap();
                worst = worst.max(gap.abs());
            }
        }
    });
    // The bodies have come most of the way together: 100 K apart at the
    // start, and the midpoint rule's factor (1 - a) / (1 + a) a step with
    // a = dt G (1/C1 + 1/C2) / 2 leaves 0.0203 K after 400.
    let apart = split.state::<f64>(hot).unwrap() - split.state::<f64>(cold).unwrap();
    assert!((apart - 0.0203).abs() <= 1.0e-4, "{apart}");
    // The two graphs keep step. The kernel accepts a temperature it solves
    // for within 1e-12 of its size; through 1.5 W/K over a step and the
    // smaller heat capacity, and at every step, that allows 2.4e-9 K.
    let allowed = f64::from(steps) * 1.5 * dt * 1.0e-12 * t1 / c1;
    eprintln!("through a heat node against the series conductance: {worst:.3e} K (allowed {allowed:.3e})");
    assert!(worst <= allowed, "{worst}");
    // Each branch generates its own entropy; together, the series
    // conductance's.
    let series = system_entropy(&split);
    assert!(series > 1.0);
    assert!((series - system_entropy(&whole)).abs() <= 1.0e-12 * series, "{series}");

    // With 30 W on the node's source every joule is in the two bodies, and
    // the heat does not divide by a fixed fraction. In the first step the
    // cold body also takes what the hot one gives up through the node,
    // 5.7 times the source; in the end both warm at one rate, the cold body
    // takes C2 / (C1 + C2) = 0.706 of the source, and the node is 4.4 K and
    // 3.5 K above them.
    let power = 30.0;
    let (mut heated, hot, cold, source) = through_node(Some(power));
    let mut shares = Vec::new();
    let mut above = [0.0_f64; 2];
    with_cx(|cx| {
        for _ in 0..steps {
            let before = *heated.state::<f64>(cold).unwrap();
            heated.step(dt, cx).unwrap();
            let after = *heated.state::<f64>(cold).unwrap();
            shares.push(c2 * (after - before) / (dt * power));
            let node_k = heated.last_output(source)[0];
            above = [
                node_k - *heated.state::<f64>(hot).unwrap(),
                node_k - *heated.state::<f64>(cold).unwrap(),
            ];
        }
    });
    let given = power * dt * f64::from(steps);
    let kept = c1 * (heated.state::<f64>(hot).unwrap() - t1) + c2 * (heated.state::<f64>(cold).unwrap() - t2);
    assert!((kept - given).abs() <= 1.0e-12 * given, "{kept} vs {given}");
    assert!(heated.ledger().energy_residual_j().abs() <= 1.0e-12 * given);
    let (first, last) = (shares[0], shares[shares.len() - 1]);
    eprintln!("a heated node: the cold body's share {first:.4} of the source at first, {last:.4} at last; the node {above:?} K above the bodies");
    assert!(first > 5.0 && first < 6.0, "{first}");
    assert!((last - c2 / (c1 + c2)).abs() <= 2.0e-3, "{last}");
    assert!((above[0] - 4.41).abs() <= 0.1 && (above[1] - 3.53).abs() <= 0.1, "{above:?}");
}

fn system_entropy(system: &System) -> f64 {
    system.ledger().entropy_generated_j_per_k()
}

#[test]
fn a_driven_damped_oscillator_heats_its_body_by_exactly_what_the_damper_absorbed() {
    let (m, k, c) = (0.25_f64, 400.0_f64, 1.5_f64);
    let push = Waveform::sine(0.8, 180.0 / (2.0 * PI), 0.0, 0.0).unwrap();
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(m).unwrap(),
            Motion {
                position_m: 2.0e-3,
                velocity_m_s: -0.3,
            },
        )
        .unwrap();
    let spring = builder
        .add_with_tangent("spring", Spring::new(k, 0.0).unwrap(), -2.0e-3)
        .unwrap();
    let damper = builder
        .add_with_tangent("damper", Damper::new(c).unwrap(), ())
        .unwrap();
    let force = builder
        .add_with_tangent("force", Force::new(push.clone()), SourceAccount::default())
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(2.0).unwrap(), 300.0)
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let members = [
        port(&builder, spring, "a"),
        port(&builder, damper, "a"),
        port(&builder, force, "p"),
    ];
    let p = port(&builder, mass, "p");
    builder.common_flow(p, &members).unwrap();
    anchored(&mut builder, spring, "b");
    anchored(&mut builder, damper, "b");
    let (body_heat, damper_heat) = (port(&builder, body, "heat"), port(&builder, damper, "heat"));
    builder.common_effort(body_heat, &[damper_heat]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(damper_heat, &[300.0]));

    let (dt, steps) = (2.0e-4, 10_000_u32);
    step(&mut system, dt, steps);

    // Monolithic midpoint rule for m x'' = -k x - c x' + F, with F the exact
    // mean of the sine over each step.
    let (mut x, mut v, mut loss, mut work) = (2.0e-3_f64, -0.3_f64, 0.0_f64, 0.0_f64);
    let half = 0.5 * dt / m;
    let mut peak = x.abs();
    for i in 0..steps {
        let f = push.mean(f64::from(i) * dt, dt);
        let mean = (v + half * (f - k * x)) / (1.0 + half * c + 0.5 * dt * half * k);
        x += dt * mean;
        v = 2.0 * mean - v;
        loss += dt * c * mean * mean;
        work += dt * f * mean;
        peak = peak.max(x.abs());
    }
    // The force's clock is the kernel's compensated time; the reference
    // multiplies, so the two times differ by an ulp at most.
    assert_eq!(system.state::<SourceAccount>(force).unwrap().time_s(), system.time_s());
    let motion = *system.state::<Motion>(mass).unwrap();
    assert!((motion.position_m - x).abs() <= 1.0e-9 * peak, "{} vs {x}", motion.position_m);
    let speed = peak * (k / m).sqrt();
    assert!((motion.velocity_m_s - v).abs() <= 1.0e-9 * speed, "{} vs {v}", motion.velocity_m_s);

    let ledger = system.ledger();
    let dissipated = ledger.nodes[damper.index()].dissipated_j.value();
    assert!((dissipated - loss).abs() <= 1.0e-9 * loss, "{dissipated} vs {loss}");
    // The body's compensated account matches to roundoff. Its temperature,
    // stored near 300 K, gained 1e4 small increments, each rounded to an ulp
    // of 300 K.
    let account = ledger.nodes[body.index()].stored_change_j.value();
    assert!((account - dissipated).abs() <= 1.0e-12 * dissipated, "{account} vs {dissipated}");
    let warmed = 2.0 * (*system.state::<f64>(body).unwrap() - 300.0);
    let rounding = 2.0 * f64::from(steps) * f64::EPSILON * 300.0;
    assert!((warmed - dissipated).abs() <= rounding, "{warmed} vs {dissipated}");
    // The force's account is minus the work it delivered.
    let delivered = -system.state::<SourceAccount>(force).unwrap().energy_in_j;
    assert!((delivered - work).abs() <= 1.0e-9 * work.abs().max(loss));
    let scale = work.abs() + loss;
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * scale);
    assert!(ledger.stored_change_j().abs() <= 1.0e-9 * scale, "closed apart from the force's account");
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}

#[test]
fn a_heated_body_cooling_to_a_fixed_ambient_approaches_its_steady_temperature_by_the_discrete_factor() {
    let (c, g, power, ambient, start) = (500.0_f64, 2.0_f64, 30.0_f64, 293.15_f64, 293.15_f64);
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(c).unwrap(), start)
        .unwrap();
    let heater = builder
        .add_with_tangent(
            "heater",
            Heater::new(Waveform::constant(power).unwrap()).unwrap(),
            SourceAccount::default(),
        )
        .unwrap();
    let wall = builder
        .add_with_tangent("wall", Conductance::new(g).unwrap(), ())
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (h, a) = (port(&builder, heater, "heat"), port(&builder, wall, "a"));
    let body_heat = port(&builder, body, "heat");
    builder.common_effort(body_heat, &[h, a]).unwrap();
    let outside = port(&builder, wall, "b");
    builder
        .boundary(outside, Box::new(move |_, _, values| values[0] = ambient))
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(h, &[start]));
    assert!(system.guess(a, &[start]));

    let (dt, steps) = (5.0, 400_u32);
    step(&mut system, dt, steps);

    let steady = ambient + power / g;
    let factor = {
        let half = 0.5 * g * dt / c;
        (1.0 - half) / (1.0 + half)
    };
    let expected = steady + (start - steady) * factor.powi(steps as i32);
    let reached = *system.state::<f64>(body).unwrap();
    assert!((reached - expected).abs() <= 1.0e-9 * steady, "{reached} vs {expected}");

    // What the heater gave either stayed in the body or left through the
    // fixed ambient boundary.
    let ledger = system.ledger();
    let given = power * dt * f64::from(steps);
    assert!((-system.state::<SourceAccount>(heater).unwrap().energy_in_j - given).abs() <= 1.0e-9 * given);
    let kept = c * (reached - start);
    let left = -ledger.boundary_energy_in_j();
    assert!((kept + left - given).abs() <= 1.0e-9 * given, "{kept} + {left} vs {given}");
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10 * given);
}

fn heater(builder: &mut SystemBuilder, name: &str, power_w: f64) -> NodeId {
    builder
        .add_with_tangent(
            name,
            Heater::new(Waveform::constant(power_w).unwrap()).unwrap(),
            SourceAccount::default(),
        )
        .unwrap()
}

#[test]
fn two_heat_routes_into_one_body_deliver_the_sum_of_their_shares() {
    let (strong, weak) = (30.0_f64, 12.0_f64);
    let (c_body, c_other) = (40.0_f64, 25.0_f64);
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(c_body).unwrap(), 300.0)
        .unwrap();
    let other = builder
        .add_with_tangent("other", HeatCapacity::new(c_other).unwrap(), 320.0)
        .unwrap();
    let lamp = heater(&mut builder, "lamp", strong);
    let coil = heater(&mut builder, "coil", weak);
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (lamp_heat, coil_heat) = (port(&builder, lamp, "heat"), port(&builder, coil, "heat"));
    let (body_heat, other_heat) = (port(&builder, body, "heat"), port(&builder, other, "heat"));
    // All of the lamp's heat, and half of the coil's, reach the same body.
    builder
        .heat_route(lamp_heat, &[(body_heat, 1.0)], None, RouteTemperature::SinkEquivalent)
        .unwrap();
    builder
        .heat_route(
            coil_heat,
            &[(body_heat, 0.5), (other_heat, 0.25)],
            Some(heat_sink("ambient", 293.15).unwrap()),
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(lamp_heat, &[300.0]));
    assert!(system.guess(coil_heat, &[300.0]));

    let (dt, steps) = (0.5, 200_u32);
    step(&mut system, dt, steps);
    let elapsed = dt * f64::from(steps);
    let given = (strong + weak) * elapsed;

    let ledger = system.ledger();
    let stored = |name: &str| {
        ledger
            .nodes
            .iter()
            .find(|node| node.name == name)
            .unwrap()
            .stored_change_j
            .value()
    };
    let close = |actual: f64, expected: f64| {
        assert!((actual - expected).abs() <= 1.0e-12 * given, "{actual} vs {expected}");
    };
    close(stored("body"), (strong + 0.5 * weak) * elapsed);
    close(stored("other"), 0.25 * weak * elapsed);
    close(ledger.sink_heat_j(), 0.25 * weak * elapsed);
    let reached = *system.state::<f64>(body).unwrap();
    let expected = 300.0 + (strong + 0.5 * weak) * elapsed / c_body;
    assert!((reached - expected).abs() <= 1.0e-10 * expected, "{reached} vs {expected}");

    // The two routes share a receiver, so they are audited in one row: the
    // first holds the pair's defect, the second holds nothing.
    assert_eq!(ledger.connections.len(), 2);
    assert!(ledger.connections[0].energy_defect_j.value().abs() <= 1.0e-12 * given);
    assert_eq!(ledger.connections[1].energy_defect_j.value(), 0.0);
    assert_eq!(ledger.connections[1].entropy_generated_j_per_k.value(), 0.0);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-12 * given);
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}

#[test]
fn a_shared_receiver_takes_route_deliveries_and_nothing_else() {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(40.0).unwrap(), 300.0)
        .unwrap();
    let other = builder
        .add_with_tangent("other", HeatCapacity::new(25.0).unwrap(), 300.0)
        .unwrap();
    let lamp = heater(&mut builder, "lamp", 30.0);
    let coil = heater(&mut builder, "coil", 12.0);
    let third = heater(&mut builder, "third", 5.0);
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (lamp_heat, coil_heat, third_heat) = (
        port(&builder, lamp, "heat"),
        port(&builder, coil, "heat"),
        port(&builder, third, "heat"),
    );
    let (body_heat, other_heat) = (port(&builder, body, "heat"), port(&builder, other, "heat"));
    let whole = RouteTemperature::SinkEquivalent;
    builder.heat_route(lamp_heat, &[(body_heat, 1.0)], None, whole).unwrap();

    let refused = |result: Result<(), AdmissionError>, port: &str| match result {
        Err(AdmissionError::AlreadyConnected { port: named }) => assert_eq!(named, port),
        other => panic!("expected {port} to be refused, got {other:?}"),
    };
    // The same receiver twice in one call.
    refused(
        builder.heat_route(coil_heat, &[(body_heat, 0.5), (body_heat, 0.5)], None, whole),
        "body.heat",
    );
    // A source belongs to one route.
    refused(
        builder.heat_route(lamp_heat, &[(other_heat, 1.0)], None, whole),
        "lamp.heat",
    );
    // A receiver of a route takes no junction.
    refused(builder.common_effort(body_heat, &[third_heat]), "body.heat");
    // A port of a junction is not a route's receiver.
    builder.common_effort(other_heat, &[third_heat]).unwrap();
    refused(
        builder.heat_route(coil_heat, &[(other_heat, 1.0)], None, whole),
        "other.heat",
    );

    // The refused calls changed nothing: the coil is still free to join the
    // lamp at the body, and the description seals and steps.
    builder.heat_route(coil_heat, &[(body_heat, 1.0)], None, whole).unwrap();
    let mut system = builder.seal().unwrap();
    for heat in [lamp_heat, coil_heat, third_heat] {
        assert!(system.guess(heat, &[300.0]));
    }
    step(&mut system, 0.5, 10);
    let warmed = *system.state::<f64>(body).unwrap() - 300.0;
    assert!((warmed - 42.0 * 5.0 / 40.0).abs() <= 1.0e-12 * 300.0, "{warmed}");
}

const ELECTRICAL: PortKind = PortKind::ElectricalVoltageCurrent;
const ROTATION: PortKind = PortKind::RotationalTorqueAngularVelocity;

#[test]
fn the_gyrator_and_the_other_domains_conform_and_exchange_no_net_power() {
    let mut n = Numbers(20_261_111);
    let machine = Gyrator::new(ELECTRICAL, ROTATION, 0.08).unwrap();
    let rotor = Mass::in_domain(ROTATION, 3.0e-4).unwrap();
    let inductor = Mass::in_domain(ELECTRICAL, 2.0e-3).unwrap();
    let capacitor = Spring::in_domain(ELECTRICAL, 1.0 / 4.7e-6, 0.0).unwrap();
    let resistor = Damper::in_domain(ELECTRICAL, 1.2).unwrap();
    let shaft = Spring::in_domain(ROTATION, 40.0, 5.0).unwrap();
    let supply = Force::in_domain(ELECTRICAL, Waveform::sine(12.0, 50.0, 0.3, 1.0).unwrap());
    for _ in 0..100 {
        let dt = n.next(1.0e-6, 1.0e-3);
        let (current, speed) = (n.next(-20.0, 20.0), n.next(-300.0, 300.0));
        let outputs = conforms(&machine, &(), &(), dt, &[current, speed]);
        assert_eq!(outputs, vec![0.08 * speed, -0.08 * current]);
        // The power entering at one port leaves at the other.
        let net = current * outputs[0] + speed * outputs[1];
        assert!(net.abs() <= 4.0 * f64::EPSILON * (0.08 * current * speed).abs(), "{net}");

        let state = Motion {
            position_m: n.next(-1.0, 1.0),
            velocity_m_s: speed,
        };
        let stale = Motion {
            position_m: 9.0,
            velocity_m_s: -9.0,
        };
        let torque = n.next(-2.0, 2.0);
        let mean = conforms(&rotor, &state, &stale, dt, &[torque]);
        let expected = speed + 0.5 * dt * torque / 3.0e-4;
        assert!((mean[0] - expected).abs() <= 1.0e-12 * (1.0 + speed.abs()));
        let state = Motion {
            position_m: 0.0,
            velocity_m_s: current,
        };
        conforms(&inductor, &state, &stale, dt, &[n.next(-24.0, 24.0)]);
        conforms(&capacitor, &n.next(-1.0e-5, 1.0e-5), &0.5, dt, &[current, 0.0]);
        conforms(&shaft, &n.next(-0.2, 0.2), &0.5, dt, &[speed, n.next(-300.0, 300.0)]);
        let heat = conforms(&resistor, &(), &(), dt, &[current, 0.0, n.next(250.0, 400.0)]);
        assert!((heat[2] + 1.2 * current * current).abs() <= 1.0e-12 * (1.0 + current * current));
        conforms(
            &supply,
            &SourceAccount::starting_at(n.next(0.0, 1.0)),
            &SourceAccount::default(),
            dt,
            &[current],
        );
    }
    for ratio in [0.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            Gyrator::new(ELECTRICAL, ROTATION, ratio),
            Err(ElementError::Parameter { name: "ratio", .. })
        ));
    }
}

struct Motor {
    system: System,
    inductor: NodeId,
    rotor: NodeId,
    supply: NodeId,
    case: NodeId,
}

const SUPPLY_V: f64 = 12.0;
const INDUCTANCE_H: f64 = 2.0e-3;
const RESISTANCE_OHM: f64 = 1.2;
const MACHINE_N_M_A: f64 = 0.08;
const ROTOR_KG_M2: f64 = 3.0e-4;
const BEARING_N_M_S: f64 = 2.0e-4;
const CASE_J_K: f64 = 60.0;

/// A DC motor from seven elements: a supply, an inductance and a winding
/// resistance in one current loop; a gyrator to a rotor with a viscous
/// bearing; winding and bearing heat into one thermal body.
fn motor() -> Motor {
    let mut builder = SystemBuilder::new();
    let rest = Motion {
        position_m: 0.0,
        velocity_m_s: 0.0,
    };
    let inductor = builder
        .add_with_tangent("inductor", Mass::in_domain(ELECTRICAL, INDUCTANCE_H).unwrap(), rest)
        .unwrap();
    let rotor = builder
        .add_with_tangent("rotor", Mass::in_domain(ROTATION, ROTOR_KG_M2).unwrap(), rest)
        .unwrap();
    let supply = builder
        .add_with_tangent(
            "supply",
            Force::in_domain(ELECTRICAL, Waveform::constant(SUPPLY_V).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    let winding = builder
        .add_with_tangent("winding", Damper::in_domain(ELECTRICAL, RESISTANCE_OHM).unwrap(), ())
        .unwrap();
    let machine = builder
        .add_with_tangent(
            "machine",
            Gyrator::new(ELECTRICAL, ROTATION, MACHINE_N_M_A).unwrap(),
            (),
        )
        .unwrap();
    let bearing = builder
        .add_with_tangent("bearing", Damper::in_domain(ROTATION, BEARING_N_M_S).unwrap(), ())
        .unwrap();
    let case = builder
        .add_with_tangent("case", HeatCapacity::new(CASE_J_K).unwrap(), 293.15)
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();

    // A rotor port does not belong in the current loop: the kernel's port
    // kinds refuse the wiring mistake, and the refused call changes nothing.
    let wrong = builder.common_flow(
        port(&builder, inductor, "p"),
        &[port(&builder, supply, "p"), port(&builder, bearing, "a")],
    );
    assert!(
        matches!(&wrong, Err(AdmissionError::Mismatch { b, .. }) if b == "bearing.a"),
        "{wrong:?}"
    );

    builder
        .common_flow(
            port(&builder, inductor, "p"),
            &[
                port(&builder, supply, "p"),
                port(&builder, winding, "a"),
                port(&builder, machine, "a"),
            ],
        )
        .unwrap();
    builder
        .common_flow(
            port(&builder, rotor, "p"),
            &[port(&builder, machine, "b"), port(&builder, bearing, "a")],
        )
        .unwrap();
    anchored(&mut builder, winding, "b");
    anchored(&mut builder, bearing, "b");
    let case_heat = port(&builder, case, "heat");
    let (winding_heat, bearing_heat) = (
        port(&builder, winding, "heat"),
        port(&builder, bearing, "heat"),
    );
    for source in [winding_heat, bearing_heat] {
        builder
            .heat_route(source, &[(case_heat, 1.0)], None, RouteTemperature::SinkEquivalent)
            .unwrap();
    }
    let mut system = builder.seal().unwrap();
    assert!(system.guess(winding_heat, &[293.15]));
    assert!(system.guess(bearing_heat, &[293.15]));
    Motor {
        system,
        inductor,
        rotor,
        supply,
        case,
    }
}

#[test]
fn a_dc_motor_composed_from_elements_matches_the_monolithic_midpoint_rule() {
    let mut motor = motor();
    let (dt, steps) = (2.0e-5_f64, 20_000_u32);
    step(&mut motor.system, dt, steps);

    // L di/dt = V - R i - k w,  J dw/dt = k i - b w, by the midpoint rule:
    // two linear equations for the step means of i and w.
    let (l, r, k, j, b) = (
        INDUCTANCE_H,
        RESISTANCE_OHM,
        MACHINE_N_M_A,
        ROTOR_KG_M2,
        BEARING_N_M_S,
    );
    let (a11, a12, a21, a22) = (2.0 * l / dt + r, k, -k, 2.0 * j / dt + b);
    let determinant = a11 * a22 - a12 * a21;
    let (mut current, mut speed) = (0.0_f64, 0.0_f64);
    let (mut supplied, mut joule, mut friction) = (0.0_f64, 0.0_f64, 0.0_f64);
    for _ in 0..steps {
        let (b1, b2) = (SUPPLY_V + 2.0 * l * current / dt, 2.0 * j * speed / dt);
        let mean_i = (b1 * a22 - a12 * b2) / determinant;
        let mean_w = (a11 * b2 - a21 * b1) / determinant;
        current = 2.0 * mean_i - current;
        speed = 2.0 * mean_w - speed;
        supplied += dt * SUPPLY_V * mean_i;
        joule += dt * r * mean_i * mean_i;
        friction += dt * b * mean_w * mean_w;
    }
    let reached_i = motor
        .system
        .state::<Motion>(motor.inductor)
        .unwrap()
        .velocity_m_s;
    let reached_w = motor
        .system
        .state::<Motion>(motor.rotor)
        .unwrap()
        .velocity_m_s;
    // After 0.4 s (many electrical and several mechanical time constants)
    // the machine is near k V / (k^2 + R b) and b V / (k^2 + R b).
    let scale = k * k + r * b;
    let (steady_w, steady_i) = (k * SUPPLY_V / scale, b * SUPPLY_V / scale);
    assert!((reached_w - speed).abs() <= 1.0e-9 * steady_w, "{reached_w} vs {speed}");
    assert!((reached_i - current).abs() <= 1.0e-9 * SUPPLY_V / r, "{reached_i} vs {current}");
    assert!((reached_w - steady_w).abs() <= 0.02 * steady_w, "{reached_w} vs {steady_w}");
    assert!((reached_i - steady_i).abs() <= 0.02 * SUPPLY_V / r, "{reached_i} vs {steady_i}");

    // Electrical energy in = magnetic and kinetic energy stored + all heat,
    // and both heats are in the one thermal body.
    let ledger = motor.system.ledger();
    let delivered = -motor
        .system
        .state::<SourceAccount>(motor.supply)
        .unwrap()
        .energy_in_j;
    assert!((delivered - supplied).abs() <= 1.0e-9 * supplied, "{delivered} vs {supplied}");
    let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).unwrap();
    let warmed = account("case").stored_change_j.value();
    assert!(
        (warmed - joule - friction).abs() <= 1.0e-9 * supplied,
        "{warmed} vs {joule} + {friction}"
    );
    let stored = 0.5 * l * reached_i * reached_i + 0.5 * j * reached_w * reached_w;
    assert!((delivered - stored - warmed).abs() <= 1.0e-9 * supplied);
    // The gyrator converts and keeps nothing.
    let machine = account("machine");
    assert!(machine.port_energy_in_j.value().abs() <= 1.0e-12 * supplied);
    assert_eq!(machine.stored_change_j.value(), 0.0);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * supplied);
    let reached_t = *motor.system.state::<f64>(motor.case).unwrap();
    assert!(((reached_t - 293.15) * CASE_J_K - warmed).abs() <= 1.0e-8 * supplied);
}

// ---------------------------------------------------------------------------
// Pinned after the second independent review (2026-11-10): each test below
// fails on the code as it was before that review's findings were fixed, or
// would have passed with one of the wrong implementations it tried.
// ---------------------------------------------------------------------------

/// The receipt of one conforming step.
fn receipt_of<C>(element: &C, state: &C::State, stale: &C::State, dt_s: f64, inputs: &[f64]) -> Receipt
where
    C: Component,
    C::State: PartialEq + core::fmt::Debug,
{
    let case = Case {
        state,
        stale,
        dt_s,
        inputs,
    };
    with_cx(|cx| match check_component(element, &case, cx) {
        Ok(Observed::Stepped { receipt, .. }) => receipt,
        other => panic!("no step at {inputs:?}: {other:?}"),
    })
}

fn within(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual} vs {expected}"
    );
}

#[test]
fn dissipative_elements_report_the_entropy_of_their_heat_and_sources_their_step_mean() {
    let dt = 0.01;
    let damper = Damper::new(2.5).unwrap();
    let receipt = receipt_of(&damper, &(), &(), dt, &[1.2, -0.3, 310.0]);
    let heat = dt * 2.5 * 1.5 * 1.5;
    within(receipt.dissipated_j, heat, 1.0e-15);
    within(receipt.entropy_generated_j_per_k, heat / 310.0, 1.0e-15);

    let wall = Conductance::new(0.8).unwrap();
    let receipt = receipt_of(&wall, &(), &(), dt, &[350.0, 300.0]);
    within(
        receipt.entropy_generated_j_per_k,
        dt * 0.8 * 50.0 * 50.0 / (350.0 * 300.0),
        1.0e-14,
    );

    // A ramp from 0 to 10 over one second, asked over [0.2, 0.6]: the mean
    // is 4, the value at the start is 2.
    let ramp = || Waveform::table(vec![(0.0, 0.0), (1.0, 10.0)]).unwrap();
    let started = SourceAccount::starting_at(0.2);
    let heater = Heater::new(ramp()).unwrap();
    let outputs = conforms(&heater, &started, &SourceAccount::default(), 0.4, &[320.0]);
    within(-outputs[0], 4.0, 1.0e-15);
    let receipt = receipt_of(&heater, &started, &SourceAccount::default(), 0.4, &[320.0]);
    within(receipt.dissipated_j, 1.6, 1.0e-15);
    within(receipt.entropy_generated_j_per_k, 1.6 / 320.0, 1.0e-15);
    let push = Force::new(ramp());
    let outputs = conforms(&push, &started, &SourceAccount::default(), 0.4, &[0.7]);
    within(outputs[0], 4.0, 1.0e-15);
}

#[test]
fn waveforms_stay_exact_at_large_times_and_refuse_what_would_overflow() {
    // A constant table at a time whose spacing is coarser than the step.
    let flat = Waveform::table(vec![(0.0, 5.0)]).unwrap();
    assert_eq!(flat.mean(1000.0, 1.0e-6), 5.0);
    assert_eq!(flat.mean(1.0e6, 1.0e-3), 5.0);
    assert_eq!(flat.mean(1.0e9, 1.0e-8), 5.0);
    let ramp = Waveform::table(vec![(0.0, 0.0), (2.0e6, 2.0e6)]).unwrap();
    let (t0, dt) = (1.0e6, 1.0e-3);
    within(ramp.mean(t0, dt), 0.5 * (t0 + (t0 + dt)), 4.0 * f64::EPSILON);

    // A time that is not finite gives NaN, never a panic and never zero.
    for wave in [&flat, &ramp, &Waveform::sine(1.0, 2.0, 0.0, 0.0).unwrap()] {
        assert!(wave.value_at(f64::NAN).is_nan());
        assert!(wave.mean(f64::NAN, 0.1).is_nan());
        assert!(wave.mean(f64::INFINITY, 1.0).is_nan());
        assert!(wave.mean(0.0, f64::INFINITY).is_nan());
    }
    // Sources refuse it instead of delivering nothing.
    let lost = SourceAccount::starting_at(f64::NAN);
    let heater = Heater::new(Waveform::constant(30.0).unwrap()).unwrap();
    let mut next = SourceAccount::default();
    let mut outputs = [0.0];
    let refused = with_cx(|cx| {
        let table = Heater::new(Waveform::table(vec![(0.0, 1.0), (1.0, 2.0)]).unwrap()).unwrap();
        (
            table.propose(&lost, 0.1, &[300.0], &mut next, &mut outputs, cx),
            heater.propose(&lost, 0.1, &[300.0], &mut next, &mut outputs, cx).is_err(),
            Force::new(Waveform::sine(1.0, 2.0, 0.0, 0.0).unwrap())
                .propose(&lost, 0.1, &[1.0], &mut next, &mut outputs, cx),
        )
    });
    assert!(matches!(refused.0, Err(ElementError::Input { value }) if value.is_nan()));
    assert!(refused.1, "a clock that is not a time is refused whatever the waveform");
    assert!(matches!(refused.2, Err(ElementError::Input { value }) if value.is_nan()));

    // Parameters whose products or differences would overflow.
    assert!(matches!(
        Waveform::sine(1.0, 1.0e308, 0.0, 100.0),
        Err(ElementError::Parameter { name: "frequency_hz", .. })
    ));
    assert!(Waveform::sine(1.0e200, 1.0, 0.0, 0.0).is_err());
    assert!(Waveform::constant(1.0e200).is_err());
    assert!(matches!(
        Waveform::table(vec![(0.0, -1.0e308), (1.0, 1.0e308)]),
        Err(ElementError::Table { .. })
    ));
    assert!(Waveform::table(vec![(-1.0e308, 0.0), (1.0e308, 1.0)]).is_err());
}

#[test]
fn a_spring_between_two_ends_that_move_together_closes_its_account() {
    // The stretch is tiny beside the common travel of the two ends.
    let hitch = Spring::new(1.0e4, 0.0).unwrap();
    conforms(&hitch, &1.0e-4, &0.5, 1.0e-3, &[30.0, 30.0]);
    conforms(&hitch, &0.0, &0.5, 1.0e-3, &[5.0, 5.000_001]);
    conforms(&hitch, &1.0e-4, &0.5, 1.0e-3, &[-3.0e3, -3.0e3]);
    let shaft = Spring::in_domain(ROTATION, 1.0e4, 0.0).unwrap();
    conforms(&shaft, &1.0e-4, &0.5, 1.0e-3, &[300.0, 300.0]);

    // A car and a trailer at 30 m/s on a stretched hitch: 2,000 steps.
    let mut builder = SystemBuilder::new();
    let moving = |position_m| Motion {
        position_m,
        velocity_m_s: 30.0,
    };
    let car = builder
        .add_with_tangent("car", Mass::new(1.0).unwrap(), moving(0.0))
        .unwrap();
    let trailer = builder
        .add_with_tangent("trailer", Mass::new(1.0).unwrap(), moving(1.0e-4))
        .unwrap();
    let link = builder.add_with_tangent("hitch", hitch, 1.0e-4).unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (a, b) = (port(&builder, link, "a"), port(&builder, link, "b"));
    builder.common_flow(port(&builder, car, "p"), &[a]).unwrap();
    builder.common_flow(port(&builder, trailer, "p"), &[b]).unwrap();
    let mut system = builder.seal().unwrap();
    step(&mut system, 1.0e-3, 2_000);
    let (front, back) = (
        *system.state::<Motion>(car).unwrap(),
        *system.state::<Motion>(trailer).unwrap(),
    );
    // Momentum is kept and the pair has travelled 60 m.
    within(front.velocity_m_s + back.velocity_m_s, 60.0, 1.0e-13);
    within(front.position_m + back.position_m, 120.0 + 1.0e-4, 1.0e-13);
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-12);
}

#[test]
fn a_system_reaching_equilibrium_keeps_stepping() {
    // Two bodies through a wall, marched far past the point where their
    // difference is below the rounding of their temperatures.
    let (c_hot, c_cold, g) = (10.0_f64, 12.0_f64, 0.7_f64);
    let mut builder = SystemBuilder::new();
    let hot = builder
        .add_with_tangent("hot", HeatCapacity::new(c_hot).unwrap(), 400.0)
        .unwrap();
    let cold = builder
        .add_with_tangent("cold", HeatCapacity::new(c_cold).unwrap(), 350.0)
        .unwrap();
    let wall = builder
        .add_with_tangent("wall", Conductance::new(g).unwrap(), ())
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (a, b) = (port(&builder, wall, "a"), port(&builder, wall, "b"));
    builder.common_effort(port(&builder, hot, "heat"), &[a]).unwrap();
    builder.common_effort(port(&builder, cold, "heat"), &[b]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(a, &[400.0]));
    assert!(system.guess(b, &[350.0]));
    step(&mut system, 0.5, 2_000);
    let mixed = (c_hot * 400.0 + c_cold * 350.0) / (c_hot + c_cold);
    within(*system.state::<f64>(hot).unwrap(), mixed, 1.0e-13);
    within(*system.state::<f64>(cold).unwrap(), mixed, 1.0e-13);
    let ledger = system.ledger();
    assert!(ledger.energy_residual_j().abs() <= 1.0e-10);
    // Entropy of mixing, to the accuracy of the midpoint rule's account.
    let expected = c_hot * (mixed / 400.0).ln() + c_cold * (mixed / 350.0).ln();
    within(ledger.entropy_generated_j_per_k(), expected, 1.0e-4);

    // Two masses at 30 and 31 m/s on a damper approach one speed.
    let mut builder = SystemBuilder::new();
    let at = |velocity_m_s| Motion {
        position_m: 0.0,
        velocity_m_s,
    };
    let slow = builder
        .add_with_tangent("slow", Mass::new(2.0).unwrap(), at(30.0))
        .unwrap();
    let fast = builder
        .add_with_tangent("fast", Mass::new(3.0).unwrap(), at(31.0))
        .unwrap();
    let drag = builder
        .add_with_tangent("drag", Damper::new(4.0).unwrap(), ())
        .unwrap();
    let pad = builder
        .add_with_tangent("pad", HeatCapacity::new(5.0).unwrap(), 300.0)
        .unwrap();
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    builder
        .common_flow(port(&builder, slow, "p"), &[port(&builder, drag, "a")])
        .unwrap();
    builder
        .common_flow(port(&builder, fast, "p"), &[port(&builder, drag, "b")])
        .unwrap();
    let heat = port(&builder, drag, "heat");
    builder
        .common_effort(port(&builder, pad, "heat"), &[heat])
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    step(&mut system, 0.01, 6_000);
    let common = (2.0 * 30.0 + 3.0 * 31.0) / 5.0;
    within(system.state::<Motion>(slow).unwrap().velocity_m_s, common, 1.0e-13);
    within(system.state::<Motion>(fast).unwrap().velocity_m_s, common, 1.0e-13);
    // The kinetic energy of the relative motion is in the pad.
    let lost = 0.5 * (2.0 * 3.0 / 5.0) * 1.0;
    let warmed = (*system.state::<f64>(pad).unwrap() - 300.0) * 5.0;
    within(warmed, lost, 1.0e-9);
}

#[test]
fn bodies_in_uniform_motion_on_relaxed_springs_keep_stepping() {
    // Three unequal masses at one speed: every force is zero, and what the
    // solve sees is only the rounding of their equal velocities.
    for speed in [1.0, 0.3, 30.0, -7.5] {
        let mut builder = SystemBuilder::new();
        let at = |position_m| Motion {
            position_m,
            velocity_m_s: speed,
        };
        let masses = [(300.0, "body"), (40.0, "wheel"), (5.0, "hub")].map(|(kg, name)| {
            builder
                .add_with_tangent(name, Mass::new(kg).unwrap(), at(0.0))
                .unwrap()
        });
        let upper = builder
            .add_with_tangent("upper", Spring::new(2.0e4, 0.0).unwrap(), 0.0)
            .unwrap();
        let lower = builder
            .add_with_tangent("lower", Spring::new(1.8e5, 0.0).unwrap(), 0.0)
            .unwrap();
        let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
        builder
            .common_flow(port(&builder, masses[0], "p"), &[port(&builder, upper, "a")])
            .unwrap();
        builder
            .common_flow(
                port(&builder, masses[1], "p"),
                &[port(&builder, upper, "b"), port(&builder, lower, "a")],
            )
            .unwrap();
        builder
            .common_flow(port(&builder, masses[2], "p"), &[port(&builder, lower, "b")])
            .unwrap();
        let mut system = builder.seal().unwrap();
        step(&mut system, 1.0e-3, 1_000);
        for node in masses {
            let motion = *system.state::<Motion>(node).unwrap();
            within(motion.velocity_m_s, speed, 1.0e-12);
            within(motion.position_m, speed, 1.0e-12);
        }
        for node in [upper, lower] {
            assert!(system.state::<f64>(node).unwrap().abs() <= 1.0e-12 * speed.abs());
        }
    }
}

#[test]
fn heat_routes_that_share_receivers_in_a_chain_are_one_group_in_any_order() {
    // lamp -> r1; coil -> r1, r2 and a sink; stove -> r2. The coil links the
    // other two. A junction of its own sits between the routes.
    let powers = [("lamp", 10.0_f64), ("coil", 20.0), ("stove", 40.0)];
    for order in [[0_usize, 2, 1], [1, 0, 2], [2, 1, 0]] {
        let mut builder = SystemBuilder::new();
        let body = |builder: &mut SystemBuilder, name: &str| {
            builder
                .add_with_tangent(name, HeatCapacity::new(50.0).unwrap(), 300.0)
                .unwrap()
        };
        let (r1, r2, r3) = (
            body(&mut builder, "r1"),
            body(&mut builder, "r2"),
            body(&mut builder, "r3"),
        );
        let sources = powers.map(|(name, watts)| heater(&mut builder, name, watts));
        let extra = heater(&mut builder, "extra", 3.0);
        let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
        let heat = |builder: &SystemBuilder, node| port(builder, node, "heat");
        let whole = RouteTemperature::SinkEquivalent;
        for (position, &which) in order.iter().enumerate() {
            let source = heat(&builder, sources[which]);
            match which {
                0 => builder.heat_route(source, &[(heat(&builder, r1), 1.0)], None, whole),
                1 => builder.heat_route(
                    source,
                    &[(heat(&builder, r1), 0.5), (heat(&builder, r2), 0.25)],
                    Some(heat_sink("air", 293.15).unwrap()),
                    whole,
                ),
                _ => builder.heat_route(source, &[(heat(&builder, r2), 1.0)], None, whole),
            }
            .unwrap();
            if position == 0 {
                builder
                    .common_effort(heat(&builder, r3), &[heat(&builder, extra)])
                    .unwrap();
            }
        }
        let mut system = builder.seal().unwrap();
        for node in sources.into_iter().chain([extra]) {
            let source = system.port(node, "heat").unwrap();
            assert!(system.guess(source, &[300.0]));
        }
        step(&mut system, 0.5, 100);
        let elapsed = 50.0;
        let ledger = system.ledger();
        let stored = |name: &str| {
            ledger
                .nodes
                .iter()
                .find(|node| node.name == name)
                .unwrap()
                .stored_change_j
                .value()
        };
        within(stored("r1"), (10.0 + 10.0) * elapsed, 1.0e-13);
        within(stored("r2"), (5.0 + 40.0) * elapsed, 1.0e-13);
        within(stored("r3"), 3.0 * elapsed, 1.0e-13);
        within(ledger.sink_heat_j(), 5.0 * elapsed, 1.0e-13);
        // Four connections: the first route holds the group, the junction
        // its own row, the two later routes nothing.
        assert_eq!(ledger.connections.len(), 4);
        let defects: Vec<f64> = ledger
            .connections
            .iter()
            .map(|row| row.energy_defect_j.value())
            .collect();
        assert!(defects[0].abs() <= 1.0e-12 * 70.0 * elapsed, "{defects:?}");
        assert_eq!((defects[2], defects[3]), (0.0, 0.0), "{order:?}");
        assert!(ledger.energy_residual_j().abs() <= 1.0e-12 * 73.0 * elapsed);
    }
}

#[test]
fn a_sink_named_by_two_routes_is_one_row_at_one_temperature() {
    let mut builder = SystemBuilder::new();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(40.0).unwrap(), 300.0)
        .unwrap();
    let lamp = heater(&mut builder, "lamp", 30.0);
    let coil = heater(&mut builder, "coil", 12.0);
    let third = heater(&mut builder, "third", 5.0);
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let body_heat = port(&builder, body, "heat");
    let whole = RouteTemperature::SinkEquivalent;
    let air = || Some(heat_sink("air", 293.15).unwrap());
    builder
        .heat_route(port(&builder, lamp, "heat"), &[(body_heat, 0.5)], air(), whole)
        .unwrap();
    builder
        .heat_route(port(&builder, coil, "heat"), &[(body_heat, 0.25)], air(), whole)
        .unwrap();
    // The same name at another temperature is another place: refused, and
    // the builder is as it was.
    let refused = builder.heat_route(
        port(&builder, third, "heat"),
        &[(body_heat, 0.5)],
        Some(heat_sink("air", 280.0).unwrap()),
        whole,
    );
    assert!(
        matches!(
            &refused,
            Err(AdmissionError::SinkTemperature { name, declared_k, repeated_k })
                if name == "air" && *declared_k == 293.15 && *repeated_k == 280.0
        ),
        "{refused:?}"
    );
    builder
        .heat_route(port(&builder, third, "heat"), &[(body_heat, 1.0)], None, whole)
        .unwrap();
    let mut system = builder.seal().unwrap();
    for node in [lamp, coil, third] {
        let source = system.port(node, "heat").unwrap();
        assert!(system.guess(source, &[300.0]));
    }
    step(&mut system, 0.5, 20);
    let ledger = system.ledger();
    assert_eq!(ledger.sinks.len(), 1);
    assert_eq!(ledger.sinks[0].name, "air");
    within(ledger.sinks[0].heat_j.value(), (15.0 + 9.0) * 10.0, 1.0e-13);
}

#[test]
fn a_publishing_heat_capacity_reports_the_temperature_its_step_began_with() {
    let plain = HeatCapacity::new(40.0).unwrap();
    let publishing = HeatCapacity::new(40.0).unwrap().publishing();
    let names: Vec<&str> = publishing.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(names, ["heat", "temperature"]);
    // Asking twice adds the signal once.
    assert_eq!(publishing.clone().publishing().ports().len(), 2);
    let mut n = Numbers(20_261_111);
    for _ in 0..50 {
        let (start, heat_w, dt) = (n.next(250.0, 900.0), n.next(-50.0, 50.0), n.next(1.0e-5, 2.0e-2));
        let outputs = conforms(&publishing, &start, &1.0, dt, &[heat_w]);
        // The mean temperature on the heat port, and the start temperature
        // (not the mean, not the end) on the signal.
        let end = start + dt * heat_w / 40.0;
        assert_eq!(outputs.len(), 2);
        assert_eq!(outputs[0].to_bits(), (0.5 * (start + end)).to_bits());
        assert_eq!(outputs[1].to_bits(), start.to_bits());
        // The signal carries nothing: state and receipt are those of the
        // body that does not publish.
        assert_eq!(
            receipt_of(&publishing, &start, &1.0, dt, &[heat_w]),
            receipt_of(&plain, &start, &1.0, dt, &[heat_w])
        );
        assert_eq!(conforms(&plain, &start, &1.0, dt, &[heat_w]), [outputs[0]]);
    }
}

// Streams: gas volumes joined by restrictions.

const AIR_J_KG_K: f64 = 287.0;
const RATIO: f64 = 1.4;

fn air(carried: usize) -> GasVolume {
    GasVolume::new(AIR_J_KG_K, RATIO, carried).unwrap()
}

/// Every entry of the port tangent against a central difference of the
/// port map, input `j` moved on its own by `moves[j]`. An entry must agree
/// to 1e-6 of its own size, or to the rounding of the difference.
fn tangent_entries<C>(element: &C, state: &C::State, dt_s: f64, inputs: &[f64], moves: &[f64])
where
    C: PortTangent,
    C::State: Clone,
{
    let n_out: usize = element.ports().iter().map(PortSpec::output_width).sum();
    with_cx(|cx| {
        let at = |x: &[f64]| {
            let mut slot = state.clone();
            let mut out = vec![0.0; n_out];
            match element.propose(state, dt_s, x, &mut slot, &mut out, cx) {
                Ok(Proposal::Step(_)) => out,
                _ => panic!("no step at {x:?}"),
            }
        };
        let base = at(inputs);
        for (column, by) in moves.iter().enumerate() {
            let mut moved = inputs.to_vec();
            moved[column] = inputs[column] + by;
            let plus = at(&moved);
            moved[column] = inputs[column] - by;
            let minus = at(&moved);
            let mut direction = vec![0.0; inputs.len()];
            direction[column] = 1.0;
            let mut tangent = vec![f64::NAN; n_out];
            assert!(element.jvp(state, dt_s, inputs, &direction, &mut tangent, cx).is_ok());
            for row in 0..n_out {
                let difference = (plus[row] - minus[row]) / (2.0 * by);
                let size = difference.abs().max(tangent[row].abs());
                let rounding = 16.0 * f64::EPSILON * base[row].abs() / by;
                assert!(
                    (difference - tangent[row]).abs() <= 1.0e-6 * size + rounding,
                    "d output {row} / d input {column}: {} against the difference {difference}",
                    tangent[row]
                );
            }
        }
    });
}

#[test]
fn gas_volumes_and_restrictions_conform_and_carry_what_the_upstream_side_holds() {
    let mut n = Numbers(20_261_112);
    let closed = air(2);
    let walled = air(2).with_wall();
    let pipe = Restriction::new(3.0e-8, 2).unwrap();
    let stale = closed.contents(1.0e5, 300.0, 1.0, &[0.5, 0.5]).unwrap();
    for _ in 0..100 {
        let dt = n.next(1.0e-5, 2.0e-4);
        let (pressure, volume) = (n.next(5.0e4, 5.0e5), n.next(0.05, 0.5));
        let fractions = [n.next(0.05, 0.4), n.next(0.05, 0.4)];
        let state = closed
            .contents(pressure, n.next(250.0, 900.0), volume, &fractions)
            .unwrap();
        let rate = n.next(-2.0e-3, 2.0e-3);
        let (enthalpy_w, heat_w) = (rate * n.next(2.0e5, 9.0e5), n.next(-50.0, 50.0));
        let mut inputs = vec![
            rate,
            enthalpy_w,
            rate * n.next(0.0, 1.0),
            rate * n.next(0.0, 1.0),
            heat_w,
        ];
        let out = conforms(&closed, &state, &stale, dt, &inputs);
        // Independent: at a fixed volume p = (ratio - 1) U / V, and U moves
        // by the enthalpy and heat that entered, so the mean pressure is
        let mean = pressure + 0.5 * (RATIO - 1.0) * dt * (enthalpy_w + heat_w) / volume;
        within(out[0], mean, 1.0e-13);
        // The temperature and each fraction returned are the means of their
        // values before and after the step, and the specific enthalpy is
        // c_p times that temperature.
        let capacity = closed.heat_capacity_j_kg_k();
        let mass = state.mass_kg + dt * rate;
        let energy = state.energy_j + dt * (enthalpy_w + heat_w);
        let kelvin =
            0.5 * (state.energy_j / (state.mass_kg * capacity) + energy / (mass * capacity));
        within(out[4], kelvin, 1.0e-13);
        within(out[1], RATIO * capacity * kelvin, 1.0e-13);
        for index in 0..2 {
            let after = (state.carried_kg[index] + dt * inputs[2 + index]) / mass;
            within(out[2 + index], 0.5 * (fractions[index] + after), 1.0e-13);
        }
        // The kit's tangent check takes one norm over all outputs, which
        // leaves the small rows of a port whose coefficients are different
        // quantities loosely pinned. Here every entry is checked on its own,
        // each input moved by 1e-4 of what would double the contents.
        let per_second = 1.0e-4 / dt;
        let mut moves = vec![
            per_second * state.mass_kg,
            per_second * state.energy_j,
            per_second * state.mass_kg,
            per_second * state.mass_kg,
            per_second * state.energy_j,
        ];
        tangent_entries(&closed, &state, dt, &inputs, &moves);

        // Behind a wall the stream and the wall see one pressure.
        inputs.push(n.next(-1.0e-2, 1.0e-2));
        let out = conforms(&walled, &state, &stale, dt, &inputs);
        assert_eq!(out[0], out[5]);
        moves.push(per_second * state.volume_m3);
        tangent_entries(&walled, &state, dt, &inputs, &moves);

        let mut side = || {
            [
                n.next(5.0e4, 5.0e5),
                n.next(2.0e5, 9.0e5),
                n.next(0.0, 1.0),
                n.next(0.0, 1.0),
            ]
        };
        let (a, b) = (side(), side());
        let both: Vec<f64> = a.iter().chain(&b).copied().collect();
        let out = conforms(&pipe, &(), &(), dt, &both);
        // What crosses is the mass rate times the values of the side the
        // mass leaves, and it enters one port as it leaves the other.
        let upstream = if a[0] >= b[0] { a } else { b };
        within(out[0], 3.0e-8 * (a[0] - b[0]), 1.0e-15);
        for index in 1..4 {
            assert_eq!(out[index], out[0] * upstream[index]);
        }
        for index in 0..4 {
            assert_eq!(out[4 + index], -out[index]);
        }
        // The map is linear in each input away from zero flow.
        if (a[0] - b[0]).abs() > 1.0e3 {
            let moves = [100.0, 300.0, 1.0e-3, 1.0e-3, 100.0, 300.0, 1.0e-3, 1.0e-3];
            tangent_entries(&pipe, &(), dt, &both, &moves);
        }
    }
    // At equal pressures nothing crosses, whatever the two sides hold. The
    // port map has a kink there, so only the mandatory checks apply.
    let level = [2.0e5, 3.0e5, 0.2, 0.7, 2.0e5, 8.0e5, 0.9, 0.1];
    let case = Case {
        state: &(),
        stale: &(),
        dt_s: 1.0e-3,
        inputs: &level,
    };
    with_cx(|cx| match check_component(&pipe, &case, cx) {
        Ok(Observed::Stepped { outputs, .. }) => {
            assert!(outputs.iter().all(|value| *value == 0.0), "{outputs:?}");
        }
        other => panic!("no step at equal pressures: {other:?}"),
    });
}

#[test]
fn gas_volumes_and_restrictions_refuse_what_they_cannot_represent() {
    for (gas_constant, ratio) in [(0.0, 1.4), (-287.0, 1.4), (287.0, 1.0), (287.0, 0.9), (f64::NAN, 1.4), (287.0, f64::INFINITY)] {
        assert!(matches!(
            GasVolume::new(gas_constant, ratio, 0),
            Err(ElementError::Parameter { .. })
        ));
    }
    for conductance in [0.0, -1.0e-8, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            Restriction::new(conductance, 0),
            Err(ElementError::Parameter { .. })
        ));
    }
    let vessel = air(1);
    assert!(vessel.contents(1.0e5, 300.0, 0.01, &[]).is_err());
    assert!(vessel.contents(1.0e5, 300.0, 0.01, &[1.5]).is_err());
    assert!(vessel.contents(0.0, 300.0, 0.01, &[0.5]).is_err());
    assert!(vessel.contents(1.0e5, -1.0, 0.01, &[0.5]).is_err());
    assert!(vessel.contents(1.0e5, 300.0, 0.0, &[0.5]).is_err());
    let names: Vec<&str> = vessel.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(names, ["stream", "heat"]);
    // Asking twice adds the wall once.
    assert_eq!(vessel.clone().with_wall().with_wall().ports().len(), 3);

    // A step that would empty the volume, or draw more of a carried
    // quantity than it holds, is too long: the element asks for a shorter
    // one and leaves the candidate alone.
    let state = vessel.contents(1.0e5, 300.0, 0.01, &[0.25]).unwrap();
    let stale = vessel.contents(2.0e5, 400.0, 0.02, &[0.5]).unwrap();
    let enthalpy = RATIO * vessel.heat_capacity_j_kg_k() * 300.0;
    with_cx(|cx| {
        for inputs in [
            [-2.0 * state.mass_kg, -2.0 * state.mass_kg * enthalpy, -0.5 * state.mass_kg, 0.0],
            [-0.5 * state.mass_kg, -0.5 * state.mass_kg * enthalpy, -0.3 * state.mass_kg, 0.0],
            [0.0, 0.0, 0.0, -2.0 * state.energy_j],
        ] {
            let mut next = stale.clone();
            let mut outputs = [7.0; 4];
            let proposal = vessel
                .propose(&state, 1.0, &inputs, &mut next, &mut outputs, cx)
                .unwrap();
            assert!(matches!(proposal, Proposal::Shorten { dt_s } if dt_s == 0.5), "{inputs:?}");
            assert_eq!(next, stale);
            assert_eq!(outputs, [7.0; 4]);
        }
        // Contents that are not a gas are refused, by name.
        let mut next = stale.clone();
        let empty = GasContents { mass_kg: 0.0, ..state.clone() };
        assert!(matches!(
            vessel.propose(&empty, 1.0e-3, &[0.0; 4], &mut next, &mut [0.0; 4], cx),
            Err(ElementError::Parameter { name: "mass_kg", .. })
        ));
        let short = GasContents { carried_kg: vec![], ..state.clone() };
        assert!(vessel.propose(&short, 1.0e-3, &[0.0; 4], &mut next, &mut [0.0; 4], cx).is_err());
        assert!(matches!(
            vessel.propose(&state, 1.0e-3, &[0.0, f64::NAN, 0.0, 0.0], &mut next, &mut [0.0; 4], cx),
            Err(ElementError::Input { .. })
        ));
        assert_eq!(next, stale);
    });
}

/// Two vessels `[m, U, c0, c1]` each, joined by one restriction and heated
/// at `heat(t)` [W] each, integrated by classical RK4 with the heat of each
/// step taken at its midpoint (the heats used here are constant between
/// whole seconds). Independent of the elements: the gas law, the
/// restriction and the donor rule are written out here.
fn vessels_reference(
    mut y: [f64; 8],
    volumes: [f64; 2],
    conductance: f64,
    heat: fn(f64) -> [f64; 2],
    t_end: f64,
    steps: u32,
) -> [f64; 8] {
    let rhs = |q: [f64; 2], y: &[f64; 8]| {
        let pressure = |side: usize| (RATIO - 1.0) * y[4 * side + 1] / volumes[side];
        let rate = conductance * (pressure(0) - pressure(1));
        let from = if rate >= 0.0 { 0 } else { 4 };
        // c_p T = ratio U / m, and each carried mass goes with its fraction.
        let per_mass = [
            1.0,
            RATIO * y[from + 1] / y[from],
            y[from + 2] / y[from],
            y[from + 3] / y[from],
        ];
        let mut d = [0.0; 8];
        for index in 0..4 {
            d[index] = -rate * per_mass[index];
            d[4 + index] = rate * per_mass[index];
        }
        d[1] += q[0];
        d[5] += q[1];
        d
    };
    let h = t_end / f64::from(steps);
    let moved = |y: &[f64; 8], d: &[f64; 8], by: f64| {
        let mut out = *y;
        for (value, slope) in out.iter_mut().zip(d) {
            *value += by * slope;
        }
        out
    };
    for step in 0..steps {
        let q = heat((f64::from(step) + 0.5) * h);
        let k1 = rhs(q, &y);
        let k2 = rhs(q, &moved(&y, &k1, 0.5 * h));
        let k3 = rhs(q, &moved(&y, &k2, 0.5 * h));
        let k4 = rhs(q, &moved(&y, &k3, h));
        for index in 0..8 {
            y[index] += h * (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index]) / 6.0;
        }
    }
    y
}

struct Vessels {
    system: System,
    nodes: [NodeId; 2],
    pipe: NodeId,
}

impl Vessels {
    /// Two vessels of air carrying two quantities, joined by a restriction,
    /// each heated at the step's midpoint value of `heat`.
    fn new(start: [GasContents; 2], conductance: f64, heat: fn(f64) -> [f64; 2]) -> Self {
        let mut builder = SystemBuilder::new();
        let pipe = builder
            .add_with_tangent("pipe", Restriction::new(conductance, 2).unwrap(), ())
            .unwrap();
        let mut nodes = Vec::new();
        for (side, (name, contents)) in ["a", "b"].into_iter().zip(start).enumerate() {
            let node = builder.add_with_tangent(name, air(2), contents).unwrap();
            builder
                .bond(
                    builder.port(node, "stream").unwrap(),
                    builder.port(pipe, name).unwrap(),
                )
                .unwrap();
            let port = builder.port(node, "heat").unwrap();
            builder
                .boundary(
                    port,
                    Box::new(move |t0, dt, values| values[0] = heat(t0 + 0.5 * dt)[side]),
                )
                .unwrap();
            nodes.push(node);
        }
        Self {
            system: builder.seal().unwrap(),
            nodes: [nodes[0], nodes[1]],
            pipe,
        }
    }

    fn contents(&self, side: usize) -> &GasContents {
        self.system.state::<GasContents>(self.nodes[side]).unwrap()
    }

    /// `[m, U, c0, c1]` of both vessels.
    fn numbers(&self) -> [f64; 8] {
        let mut out = [0.0; 8];
        for side in 0..2 {
            let held = self.contents(side);
            out[4 * side..4 * side + 4].copy_from_slice(&[
                held.mass_kg,
                held.energy_j,
                held.carried_kg[0],
                held.carried_kg[1],
            ]);
        }
        out
    }

    /// Mass rate [kg/s] from `a` to `b` during the last step.
    fn rate(&self) -> f64 {
        self.system.last_output(self.system.port(self.pipe, "a").unwrap())[0]
    }
}

/// Largest differences between two sets of vessel numbers: of mass and
/// energy, relative to the reference's value, and of the carried masses,
/// relative to their vessel's mass.
fn vessel_errors(got: &[f64; 8], want: &[f64; 8]) -> [f64; 2] {
    (0..8).fold([0.0_f64; 2], |mut worst, index| {
        let (kind, scale) = if index % 4 >= 2 {
            (1, want[index - index % 4])
        } else {
            (0, want[index])
        };
        worst[kind] = worst[kind].max((got[index] - want[index]).abs() / scale);
        worst
    })
}

fn no_heat(_: f64) -> [f64; 2] {
    [0.0, 0.0]
}

#[test]
fn two_vessels_through_a_restriction_equalize_and_the_donor_keeps_its_composition() {
    let start = || {
        [
            air(2).contents(5.0e5, 350.0, 0.01, &[0.2, 0.0]).unwrap(),
            air(2).contents(1.0e5, 290.0, 0.03, &[0.0, 0.5]).unwrap(),
        ]
    };
    let conductance = 2.0e-8;
    let initial = Vessels::new(start(), conductance, no_heat).numbers();
    let reference = vessels_reference(initial, [0.01, 0.03], conductance, no_heat, 2.0, 2_000);
    let mut errors = Vec::new();
    for steps in [100_u32, 200] {
        let mut vessels = Vessels::new(start(), conductance, no_heat);
        let dt = 2.0 / f64::from(steps);
        with_cx(|cx| {
            for _ in 0..steps {
                assert_eq!(vessels.system.step(dt, cx).unwrap().shortenings, 0);
                // Mass leaves `a` for the whole run, so `a` keeps the
                // composition it started with and receives none of what only
                // `b` holds. Nothing in either element says so.
                assert!(vessels.rate() > 0.0);
                let a = vessels.contents(0);
                within(a.carried_kg[0] / a.mass_kg, 0.2, 1.0e-11);
                assert_eq!(a.carried_kg[1], 0.0);
            }
        });
        let now = vessels.numbers();
        // Mass, each carried quantity and energy: what one vessel lost the
        // other gained.
        for index in 0..4 {
            within(now[index] + now[4 + index], initial[index] + initial[4 + index], 4.0e-15);
        }
        assert!(vessels.system.ledger().energy_residual_j().abs() <= 1.0e-10);
        let [flow, carried] = vessel_errors(&now, &reference);
        errors.push(flow.max(carried));
    }
    // Second order against the independent integration: measured 1.35e-6
    // and 3.37e-7.
    let ratio = errors[0] / errors[1];
    assert!(errors[1] < 5.0e-7 && (3.8..4.2).contains(&ratio), "{errors:?}");

    // Left alone the pressures meet at (p_a V_a + p_b V_b) / (V_a + V_b),
    // which is 2 bar here: the energy is conserved and p V = (ratio - 1) U.
    // The temperatures do not meet: nothing conducts heat between them.
    let mut vessels = Vessels::new(start(), conductance, no_heat);
    step(&mut vessels.system, 0.05, 3_200);
    let volume = air(2);
    for side in 0..2 {
        within(volume.pressure_pa(vessels.contents(side)), 2.0e5, 1.0e-7);
    }
    let (cold, warm) = (
        volume.temperature_k(vessels.contents(0)),
        volume.temperature_k(vessels.contents(1)),
    );
    assert!(cold < 290.0 && warm > 330.0, "{cold} K and {warm} K");
}

/// 800 W into `a` for the first second, then 4 kW into `b` for the next.
fn alternate_heat(t: f64) -> [f64; 2] {
    if t < 1.0 {
        [800.0, 0.0]
    } else if t < 2.0 {
        [0.0, 4_000.0]
    } else {
        [0.0, 0.0]
    }
}

#[test]
fn a_flow_that_turns_round_carries_back_what_the_other_vessel_holds() {
    let start = || {
        [
            air(2).contents(2.0e5, 300.0, 0.01, &[0.3, 0.0]).unwrap(),
            air(2).contents(2.0e5, 300.0, 0.03, &[0.0, 0.6]).unwrap(),
        ]
    };
    let conductance = 2.0e-8;
    let initial = Vessels::new(start(), conductance, alternate_heat).numbers();
    let reference =
        vessels_reference(initial, [0.01, 0.03], conductance, alternate_heat, 3.0, 30_000);
    let mut errors = Vec::new();
    for steps in [300_u32, 600] {
        let mut vessels = Vessels::new(start(), conductance, alternate_heat);
        let dt = 3.0 / f64::from(steps);
        let (mut forward, mut backward, mut turns) = (0_u32, 0_u32, 0_u32);
        let mut last = 0.0_f64;
        with_cx(|cx| {
            for index in 0..steps {
                // The step in which the flow turns is taken whole.
                assert_eq!(vessels.system.step(dt, cx).unwrap().shortenings, 0, "step {index}");
                let rate = vessels.rate();
                if rate > 0.0 {
                    forward += 1;
                    // Until the flow turns, `a` has received nothing.
                    if backward == 0 {
                        assert_eq!(vessels.contents(0).carried_kg[1], 0.0, "step {index}");
                    }
                } else if rate < 0.0 {
                    backward += 1;
                }
                if rate * last < 0.0 {
                    turns += 1;
                }
                last = rate;
            }
        });
        // Out of `a` while it is heated, and back into it once `b` is.
        assert!(forward > steps / 4 && backward > steps / 4 && turns == 1, "{forward} {backward} {turns}");
        let now = vessels.numbers();
        // After the turn `a` holds some of what only `b` had.
        assert!(now[3] > 1.0e-6, "{}", now[3]);
        // Mass and both carried quantities are conserved through the turn;
        // the energy rose by the heat supplied, 800 J and 4 kJ.
        for index in [0, 2, 3] {
            within(now[index] + now[4 + index], initial[index] + initial[4 + index], 4.0e-15);
        }
        within(now[1] + now[5], initial[1] + initial[5] + 4_800.0, 1.0e-13);
        assert!(vessels.system.ledger().energy_residual_j().abs() <= 1.0e-9);
        let [flow, carried] = vessel_errors(&now, &reference);
        // At the turn the supplier of the carried values changes, and they
        // differ between the vessels by their whole size. The step in which
        // that happens costs up to its length squared, depending on where
        // in it the turn falls, so this error is bounded and has no clean
        // ratio. Measured 1.2e-7 and 6.4e-8.
        assert!(carried <= 5.0e-3 * dt * dt, "{carried} at {steps} steps");
        errors.push(flow);
    }
    // Mass and energy stay second order through the turn: measured 6.5e-8
    // and 1.5e-8.
    let ratio = errors[0] / errors[1];
    assert!(errors[1] < 3.0e-8 && (3.8..4.6).contains(&ratio), "{errors:?}");
}

/// A line of equal vessels `[m, U, c]`, each joined to the next by the same
/// restriction, integrated by classical RK4. Independent of the elements.
fn line_reference(
    mut y: Vec<[f64; 3]>,
    volume: f64,
    conductance: f64,
    t_end: f64,
    steps: u32,
) -> Vec<[f64; 3]> {
    let rhs = |y: &[[f64; 3]]| {
        let mut d = vec![[0.0; 3]; y.len()];
        let pressure = |vessel: usize| (RATIO - 1.0) * y[vessel][1] / volume;
        for link in 0..y.len() - 1 {
            let rate = conductance * (pressure(link) - pressure(link + 1));
            let from = if rate >= 0.0 { link } else { link + 1 };
            let per_mass = [
                1.0,
                RATIO * y[from][1] / y[from][0],
                y[from][2] / y[from][0],
            ];
            for index in 0..3 {
                d[link][index] -= rate * per_mass[index];
                d[link + 1][index] += rate * per_mass[index];
            }
        }
        d
    };
    let moved = |y: &[[f64; 3]], d: &[[f64; 3]], by: f64| -> Vec<[f64; 3]> {
        y.iter()
            .zip(d)
            .map(|(value, slope)| core::array::from_fn(|index| value[index] + by * slope[index]))
            .collect()
    };
    let h = t_end / f64::from(steps);
    for _ in 0..steps {
        let k1 = rhs(&y);
        let k2 = rhs(&moved(&y, &k1, 0.5 * h));
        let k3 = rhs(&moved(&y, &k2, 0.5 * h));
        let k4 = rhs(&moved(&y, &k3, h));
        for (vessel, value) in y.iter_mut().enumerate() {
            for index in 0..3 {
                value[index] += h
                    * (k1[vessel][index]
                        + 2.0 * k2[vessel][index]
                        + 2.0 * k3[vessel][index]
                        + k4[vessel][index])
                    / 6.0;
            }
        }
    }
    y
}

/// Three 20 l vessels of air at 300 K in a line, joined by two restrictions
/// of 5e-8 kg/(s Pa); the first holds `fraction` of one carried quantity.
fn vessel_line(pressures: [f64; 3], fraction: f64) -> (System, Vec<NodeId>) {
    let mut builder = SystemBuilder::new();
    let mut vessels = Vec::new();
    for (index, (name, pressure)) in ["first", "second", "third"].into_iter().zip(pressures).enumerate() {
        let held = if index == 0 { fraction } else { 0.0 };
        let contents = air(1).contents(pressure, 300.0, 0.02, &[held]).unwrap();
        let node = builder.add_with_tangent(name, air(1), contents).unwrap();
        anchored(&mut builder, node, "heat");
        vessels.push(node);
    }
    let mut pipes = Vec::new();
    for name in ["upper", "lower"] {
        let pipe = Restriction::new(5.0e-8, 1).unwrap();
        pipes.push(builder.add_with_tangent(name, pipe, ()).unwrap());
    }
    let stream = |builder: &SystemBuilder, index: usize| builder.port(vessels[index], "stream").unwrap();
    let end = |builder: &SystemBuilder, pipe: usize, side: &str| builder.port(pipes[pipe], side).unwrap();
    builder.bond(stream(&builder, 0), end(&builder, 0, "a")).unwrap();
    builder
        .common_effort(stream(&builder, 1), &[end(&builder, 0, "b"), end(&builder, 1, "a")])
        .unwrap();
    builder.bond(end(&builder, 1, "b"), stream(&builder, 2)).unwrap();
    (builder.seal().unwrap(), vessels)
}

#[test]
fn a_quantity_carried_through_a_middle_vessel_arrives_at_second_order() {
    // At 4, 2.5 and 1 bar the flow runs one way down the line, and what the
    // first vessel holds reaches the third only through the second, which
    // receives and gives in the same step. A storage that returned its
    // fractions at the start or the end of the step, and not their mean,
    // would deliver it at first order.
    let pressures = [4.0e5, 2.5e5, 1.0e5];
    let numbers = |system: &System, nodes: &[NodeId]| -> Vec<[f64; 3]> {
        nodes
            .iter()
            .map(|node| {
                let held = system.state::<GasContents>(*node).unwrap();
                [held.mass_kg, held.energy_j, held.carried_kg[0]]
            })
            .collect()
    };
    let (system, nodes) = vessel_line(pressures, 0.4);
    let initial = numbers(&system, &nodes);
    let reference = line_reference(initial.clone(), 0.02, 5.0e-8, 2.0, 2_000);
    let mut errors = Vec::new();
    for steps in [100_u32, 200] {
        let (mut system, nodes) = vessel_line(pressures, 0.4);
        step(&mut system, 2.0 / f64::from(steps), steps);
        let now = numbers(&system, &nodes);
        // Mass, energy and the quantity are conserved down the line.
        for index in 0..3 {
            let total = |all: &[[f64; 3]]| all.iter().map(|vessel| vessel[index]).sum::<f64>();
            within(total(&now), total(&initial), 4.0e-15);
        }
        let (mut flow, mut carried) = (0.0_f64, 0.0_f64);
        for (got, want) in now.iter().zip(&reference) {
            flow = flow
                .max((got[0] - want[0]).abs() / want[0])
                .max((got[1] - want[1]).abs() / want[1]);
            carried = carried.max((got[2] - want[2]).abs() / want[0]);
        }
        errors.push([flow, carried]);
    }
    // The third vessel has received 0.41 g of it. Measured errors: mass and
    // energy 8.8e-7 and 2.2e-7, the quantity 1.6e-7 and 4.1e-8.
    assert!(reference[2][2] > 1.0e-4, "{}", reference[2][2]);
    for kind in 0..2 {
        let ratio = errors[0][kind] / errors[1][kind];
        assert!((3.7..4.3).contains(&ratio), "{errors:?}");
    }
}

#[test]
fn a_vessel_filled_from_a_line_ends_hotter_than_the_line() {
    // m c_v dT = (c_p T_line - c_v T) dm: the gas pushed in does flow work.
    // The energy balance is exact for any step, so the end state is
    // (m0 c_v T0 + m_in c_p T_line) / ((m0 + m_in) c_v) to rounding.
    let vessel = air(1);
    let start = vessel.contents(1.0e5, 300.0, 0.02, &[0.0]).unwrap();
    let (rate, line_k, seconds) = (5.0e-3, 300.0, 4.0);
    let capacity = vessel.heat_capacity_j_kg_k();
    let enthalpy = RATIO * capacity * line_k;
    let mut builder = SystemBuilder::new();
    let node = builder.add_with_tangent("vessel", air(1), start.clone()).unwrap();
    let inlet = builder.port(node, "stream").unwrap();
    // The line is all of the carried quantity.
    builder
        .boundary(
            inlet,
            Box::new(move |_, _, values| values.copy_from_slice(&[rate, rate * enthalpy, rate])),
        )
        .unwrap();
    anchored(&mut builder, node, "heat");
    let mut system = builder.seal().unwrap();
    step(&mut system, 0.01, 400);
    let end = system.state::<GasContents>(node).unwrap();
    // Four hundred additions round by about 1e-14 of the total.
    let added = rate * seconds;
    within(end.mass_kg, start.mass_kg + added, 1.0e-12);
    within(end.carried_kg[0] / end.mass_kg, added / (start.mass_kg + added), 1.0e-12);
    let expected = (start.mass_kg * capacity * 300.0 + added * enthalpy)
        / ((start.mass_kg + added) * capacity);
    within(vessel.temperature_k(end), expected, 1.0e-12);
    // Between the line's temperature and ratio times it, which is the
    // limit of filling an empty vessel.
    assert!(expected > 330.0 && expected < RATIO * line_k, "{expected}");
    within(system.ledger().boundary_energy_in_j(), added * enthalpy, 1.0e-12);
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-9);
}

#[test]
fn a_mass_on_a_gas_volume_rings_at_the_adiabatic_gas_spring_frequency() {
    // A wall of area A and mass M on a volume V at pressure p, the outside
    // pressure carried by a constant force: for small motion
    // w^2 = ratio p A^2 / (M V).
    let (area, mass_kg, pressure, volume_m3) = (0.01, 1.0, 1.0e5, 1.0e-3);
    let gas = air(0).with_wall();
    let start = gas.contents(pressure, 300.0, volume_m3, &[]).unwrap();
    let mut builder = SystemBuilder::new();
    let wall = builder
        .add_with_tangent(
            "wall",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: volume_m3 / area,
                velocity_m_s: 0.01,
            },
        )
        .unwrap();
    let cavity = builder.add_with_tangent("gas", gas.clone(), start.clone()).unwrap();
    let outside = builder
        .add_with_tangent(
            "outside",
            Force::new(Waveform::constant(-pressure * area).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    builder
        .common_flow_ratios(
            builder.port(wall, "p").unwrap(),
            &[
                (builder.port(cavity, "p").unwrap(), area),
                (builder.port(outside, "p").unwrap(), 1.0),
            ],
        )
        .unwrap();
    anchored(&mut builder, cavity, "heat");
    let closed = builder.port(cavity, "stream").unwrap();
    builder
        .boundary(closed, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    let mut system = builder.seal().unwrap();

    let dt = 1.0e-4;
    let omega = (RATIO * pressure * area * area / (mass_kg * volume_m3)).sqrt();
    let mut crossings = Vec::new();
    let mut invariant = 0.0_f64;
    with_cx(|cx| {
        let mut before = 0.01_f64;
        for index in 0..1_400_u32 {
            system.step(dt, cx).unwrap();
            let velocity = system.state::<Motion>(wall).unwrap().velocity_m_s;
            if before < 0.0 && velocity >= 0.0 {
                crossings.push((f64::from(index) + before / (before - velocity)) * dt);
            }
            before = velocity;
            // No heat and no stream: p V^ratio stays what it was.
            let held = system.state::<GasContents>(cavity).unwrap();
            let adiabat = gas.pressure_pa(held) * held.volume_m3.powf(RATIO)
                / (pressure * volume_m3.powf(RATIO));
            invariant = invariant.max((adiabat - 1.0).abs());
            // The wall's position is the gas's length.
            within(area * system.state::<Motion>(wall).unwrap().position_m, held.volume_m3, 1.0e-12);
        }
    });
    assert!(crossings.len() >= 2, "{crossings:?}");
    let period = crossings[1] - crossings[0];
    // The trapezoidal rule turns a linear oscillator by 2 atan(w dt / 2) a
    // step, which lengthens the period by (w dt)^2 / 12, 1.2e-5 here. The
    // measured period is within 7e-8 of that and 1.2e-5 from 2 pi / w.
    within(period, PI * dt / (0.5 * omega * dt).atan(), 1.0e-6);
    within(period, 2.0 * PI / omega, 2.0e-5);
    // Measured 8e-15 and 5e-19 J.
    assert!(invariant <= 1.0e-12, "{invariant}");
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-13);
}

#[test]
fn streams_meet_at_storages_and_other_wirings_are_refused() {
    // Three vessels in a line: the middle one is the node two restrictions
    // meet at.
    let mut builder = SystemBuilder::new();
    let mut vessels = Vec::new();
    for (name, pressure) in [("left", 4.0e5), ("middle", 1.0e5), ("right", 2.5e5)] {
        let contents = air(0).contents(pressure, 300.0, 0.02, &[]).unwrap();
        let node = builder.add_with_tangent(name, air(0), contents).unwrap();
        anchored(&mut builder, node, "heat");
        vessels.push(node);
    }
    let mut pipes = Vec::new();
    for name in ["first", "second"] {
        pipes.push(
            builder
                .add_with_tangent(name, Restriction::new(5.0e-8, 0).unwrap(), ())
                .unwrap(),
        );
    }
    let stream = |builder: &SystemBuilder, index: usize| builder.port(vessels[index], "stream").unwrap();
    let end = |builder: &SystemBuilder, pipe: usize, side: &str| builder.port(pipes[pipe], side).unwrap();
    let wide = builder
        .add_with_tangent("wide", Restriction::new(5.0e-8, 2).unwrap(), ())
        .unwrap();
    let mass = builder
        .add_with_tangent("mass", Mass::new(1.0).unwrap(), Motion { position_m: 0.0, velocity_m_s: 0.0 })
        .unwrap();
    let body = builder.add_with_tangent("body", HeatCapacity::new(10.0).unwrap(), 300.0).unwrap();

    // One flow cannot be handed whole to two storages.
    assert!(matches!(
        builder.common_flow(end(&builder, 0, "a"), &[stream(&builder, 0), stream(&builder, 1)]),
        Err(AdmissionError::StreamCommonFlow { port }) if port == "first.a"
    ));
    // Two storages, or two restrictions, have nothing between them.
    assert!(matches!(
        builder.bond(stream(&builder, 0), stream(&builder, 1)),
        Err(AdmissionError::Role { .. })
    ));
    assert!(matches!(
        builder.bond(end(&builder, 0, "b"), end(&builder, 1, "a")),
        Err(AdmissionError::Role { .. })
    ));
    // Streams that carry different things, and other carriers, do not join.
    assert!(matches!(
        builder.bond(stream(&builder, 0), builder.port(wide, "a").unwrap()),
        Err(AdmissionError::Mismatch { .. })
    ));
    assert!(matches!(
        builder.bond(builder.port(body, "heat").unwrap(), end(&builder, 0, "a")),
        Err(AdmissionError::Mismatch { .. })
    ));
    // A stream is not a power pair: no ratio, transformer or heat route.
    let velocity = builder.port(mass, "p").unwrap();
    assert!(matches!(
        builder.common_flow_ratios(velocity, &[(stream(&builder, 0), 2.0)]),
        Err(AdmissionError::Role { .. })
    ));
    assert!(matches!(
        builder.transformer(velocity, stream(&builder, 0), 2.0),
        Err(AdmissionError::Role { .. })
    ));
    assert!(builder
        .heat_route(end(&builder, 0, "a"), &[(builder.port(body, "heat").unwrap(), 1.0)], None, RouteTemperature::SinkEquivalent)
        .is_err());

    // Every refusal above named these ports, and none was claimed.
    builder.bond(stream(&builder, 0), end(&builder, 0, "a")).unwrap();
    builder
        .common_effort(stream(&builder, 1), &[end(&builder, 0, "b"), end(&builder, 1, "a")])
        .unwrap();
    builder.bond(end(&builder, 1, "b"), stream(&builder, 2)).unwrap();
    for (node, port) in [(wide, "a"), (wide, "b")] {
        let port = builder.port(node, port).unwrap();
        builder
            .boundary(port, Box::new(|_, _, values| values.copy_from_slice(&[1.0e5, 3.0e5, 0.5, 0.5])))
            .unwrap();
    }
    anchored(&mut builder, mass, "p");
    anchored(&mut builder, body, "heat");
    let mut system = builder.seal().unwrap();
    let total = |system: &System| {
        vessels.iter().fold([0.0, 0.0], |sum, node| {
            let held = system.state::<GasContents>(*node).unwrap();
            [sum[0] + held.mass_kg, sum[1] + held.energy_j]
        })
    };
    let before = total(&system);
    step(&mut system, 0.05, 4_000);
    let after = total(&system);
    within(after[0], before[0], 4.0e-15);
    within(after[1], before[1], 1.0e-13);
    // Equal volumes of one gas: the mean of the three pressures.
    for node in &vessels {
        let held = system.state::<GasContents>(*node).unwrap();
        within(air(0).pressure_pa(held), 2.5e5, 1.0e-7);
    }
    assert!(system.ledger().energy_residual_j().abs() <= 1.0e-9);
}

/// A heat source whose power falls with the temperature it is at,
/// `P = P0 exp(-(T - 300 K) / scale)`: friction that weakens as its
/// interface warms. Its one heat port receives the temperature and returns
/// the heat rate into the source, `-P`. Smooth, and added without a
/// tangent, so the kernel differences it. Its state is the energy it has
/// delivered, as a negative stored energy. It counts how often it is asked
/// to propose.
struct FadingSource {
    power_w: f64,
    scale_k: f64,
    ports: [PortSpec; 1],
    proposals: Arc<AtomicU32>,
}

impl FadingSource {
    fn new(power_w: f64, scale_k: f64) -> Self {
        Self {
            power_w,
            scale_k,
            ports: [PortSpec::heat(StableId::new("heat").unwrap(), PortRole::EffortIn).unwrap()],
            proposals: Arc::new(AtomicU32::new(0)),
        }
    }

    fn power_at(&self, kelvin: f64) -> f64 {
        self.power_w * (-(kelvin - 300.0) / self.scale_k).exp()
    }
}

impl Component for FadingSource {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, delivered_j: &f64) -> Result<f64, ElementError> {
        Ok(*delivered_j)
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
        cx.checkpoint()?;
        self.proposals.fetch_add(1, Ordering::Relaxed);
        let kelvin = inputs[0];
        if !(kelvin.is_finite() && kelvin > 0.0) {
            return Err(ElementError::Temperature { kelvin });
        }
        let power_w = self.power_at(kelvin);
        outputs[0] = -power_w;
        *new = old - dt_s * power_w;
        Ok(Proposal::Step(Receipt {
            stored_change_j: -dt_s * power_w,
            dissipated_j: dt_s * power_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * power_w / kelvin,
            energy_tolerance_j: 64.0 * f64::EPSILON * (old.abs() + new.abs() + dt_s * power_w)
                + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

/// Heat capacities of the two bodies of [`warming_plane`] [J/K], their
/// conductances to the plane and to the room [W/K], the source's power at
/// 300 K [W] and the temperature scale it falls with [K], and the room [K].
const WARMING_J_K: [f64; 2] = [1.0; 2];
const WARMING_PLANE_W_K: [f64; 2] = [2.0, 4.0];
const WARMING_FILM_W_K: [f64; 2] = [3.0, 5.0];
const WARMING_POWER_W: f64 = 8.0;
const WARMING_SCALE_K: f64 = 40.0;
const WARMING_ROOM_K: f64 = 300.0;

/// Two bodies on a plane that carries a heat source falling with the
/// plane's temperature, each body also losing heat to a room through a
/// film; everything at the room's temperature. The kernel's default solve
/// settings. Returns the system and the two bodies.
fn warming_plane() -> (System, Vec<NodeId>) {
    let (system, bodies, _, _) = warming_plane_counted();
    (system, bodies)
}

/// [`warming_plane`], and with it the source and the count of the
/// proposals the source has been asked for.
fn warming_plane_counted() -> (System, Vec<NodeId>, NodeId, Arc<AtomicU32>) {
    let room_k = WARMING_ROOM_K;
    let mut builder = SystemBuilder::new();
    let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
    let plane = builder
        .add_with_tangent("plane", HeatNode::new(&WARMING_PLANE_W_K).unwrap(), ())
        .unwrap();
    let mut bodies = Vec::new();
    let mut starts = Vec::new();
    for side in 0..2 {
        let body = builder
            .add_with_tangent(
                &format!("body{side}"),
                HeatCapacity::new(WARMING_J_K[side]).unwrap(),
                room_k,
            )
            .unwrap();
        let film = builder
            .add_with_tangent(
                &format!("film{side}"),
                Conductance::new(WARMING_FILM_W_K[side]).unwrap(),
                (),
            )
            .unwrap();
        let (on_plane, on_film) =
            (port(&builder, plane, &format!("body{side}")), port(&builder, film, "a"));
        builder
            .common_effort(port(&builder, body, "heat"), &[on_plane, on_film])
            .unwrap();
        builder
            .boundary(port(&builder, film, "b"), Box::new(move |_, _, values| values[0] = room_k))
            .unwrap();
        bodies.push(body);
        starts.extend([on_plane, on_film]);
    }
    let law = FadingSource::new(WARMING_POWER_W, WARMING_SCALE_K);
    let proposals = Arc::clone(&law.proposals);
    let source = builder.add("source", law, 0.0).unwrap();
    let at_source = port(&builder, source, "heat");
    builder
        .common_effort(port(&builder, plane, "source"), &[at_source])
        .unwrap();
    starts.push(at_source);
    let mut system = builder.seal().unwrap();
    for start in starts {
        assert!(system.guess(start, &[room_k]));
    }
    (system, bodies, source, proposals)
}

#[test]
fn a_step_taken_at_its_first_check_was_judged_and_its_differenced_part_proposed_again() {
    // The warming plane's source has no declared tangent, so an assembly
    // differences it: one more proposal, at a perturbed temperature, into
    // its candidate. The source counts its proposals, which says what the
    // solve did in a step:
    //   k >= 1 updates: k + 1 checks and k assemblies, 2 k + 1 proposals;
    //   taken at the first check: that check, the assembly that judges
    //   it, and one more at the accepted inputs, 3 proposals;
    //   an exact fixed point, taken without assembling: 1.
    // So a step with no iteration and three proposals was judged by slopes
    // formed in that step, and its source was proposed again afterwards.
    // What the source then publishes is its proposal at the temperature
    // the step accepted, to the bit, not the perturbed one.
    let (mut system, _, source, proposals) = warming_plane_counted();
    let law = FadingSource::new(WARMING_POWER_W, WARMING_SCALE_K);
    let at_source = system.port(source, "heat").unwrap();
    let dt = 0.25_f64;
    let mut taken_at_first_check = 0_u32;
    with_cx(|cx| {
        for step in 0..80_u32 {
            let asked_before = proposals.load(Ordering::Relaxed);
            let delivered_before_j = *system.state::<f64>(source).unwrap();
            let spent = system.step(dt, cx).unwrap().iterations;
            let asked = proposals.load(Ordering::Relaxed) - asked_before;
            if spent == 0 {
                assert_eq!(asked, 3, "step {step}");
                taken_at_first_check += 1;
            } else {
                assert_eq!(asked, 2 * spent + 1, "step {step}");
            }
            let accepted_k = system.last_input(at_source)[0];
            let expected_j = delivered_before_j - dt * law.power_at(accepted_k);
            let delivered_j = *system.state::<f64>(source).unwrap();
            assert_eq!(delivered_j.to_bits(), expected_j.to_bits(), "step {step}");
        }
    });
    assert!(taken_at_first_check >= 10, "{taken_at_first_check}");
}

/// A body whose heat capacity rises with its temperature,
/// `C(T) = c0 + c1 (T - 300 K)`. The temperature it returns on its port,
/// the mean of its two end temperatures, is therefore not linear in the
/// heat rate it receives. Smooth, and added without a tangent.
struct StiffeningBody {
    c0_j_k: f64,
    c1_j_k2: f64,
    ports: [PortSpec; 1],
}

impl StiffeningBody {
    fn new(c0_j_k: f64, c1_j_k2: f64) -> Self {
        Self {
            c0_j_k,
            c1_j_k2,
            ports: [PortSpec::heat(StableId::new("heat").unwrap(), PortRole::FlowIn).unwrap()],
        }
    }

    /// Stored energy above 300 K [J].
    fn energy_j(&self, kelvin: f64) -> f64 {
        let above = kelvin - 300.0;
        self.c0_j_k * above + 0.5 * self.c1_j_k2 * above * above
    }
}

impl Component for StiffeningBody {
    type State = f64;
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, kelvin: &f64) -> Result<f64, ElementError> {
        Ok(self.energy_j(*kelvin))
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
        cx.checkpoint()?;
        let heat_j = dt_s * inputs[0];
        // The end temperature holds the energy the step began with and the
        // heat received: the root of c0 x + c1 x^2 / 2 = E above 300 K.
        let energy_j = self.energy_j(*old) + heat_j;
        let under = self.c0_j_k * self.c0_j_k + 2.0 * self.c1_j_k2 * energy_j;
        if !(under > 0.0) {
            return Ok(Proposal::Shorten { dt_s: 0.5 * dt_s });
        }
        let next = 300.0 + (under.sqrt() - self.c0_j_k) / self.c1_j_k2;
        *new = next;
        outputs[0] = 0.5 * (old + next);
        Ok(Proposal::Step(Receipt {
            stored_change_j: heat_j,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 64.0
                * f64::EPSILON
                * (self.energy_j(*old).abs()
                    + self.energy_j(next).abs()
                    + heat_j.abs()
                    + self.c0_j_k * (old.abs() + next.abs()))
                + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

#[test]
fn a_block_refused_by_the_second_test_after_updates_names_the_open_input() {
    // A body whose heat capacity rises with temperature (1 J/K at 300 K,
    // 0.01 J/K more per kelvin) behind a wall of 1e5 W/K, the far side 3 K
    // warmer, one step of 1 ms from a guess of 300 K at the wall. The
    // body's temperature is not linear in its heat, so Newton needs four
    // updates. Allowed three, the block stands at its fourth check with
    // the wall's temperature input open by about 4e-11 K: inside the first
    // test's 3e-10 K, and through the wall about 4e-6 W of heat the body
    // has not been told of. The refusal comes from the second test, at a
    // check after updates, and names that input and its residual; the
    // body's heat input, first among the block's unknowns, is closed.
    let build = |max_iterations: Option<u32>| {
        let mut builder = SystemBuilder::new();
        let body = builder
            .add("body", StiffeningBody::new(1.0, 0.01), 300.0)
            .unwrap();
        let wall = builder
            .add_with_tangent("wall", Conductance::new(1.0e5).unwrap(), ())
            .unwrap();
        let near = builder.port(wall, "a").unwrap();
        builder
            .common_effort(builder.port(body, "heat").unwrap(), &[near])
            .unwrap();
        builder
            .boundary(
                builder.port(wall, "b").unwrap(),
                Box::new(|_, _, values| values[0] = 303.0),
            )
            .unwrap();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(near, &[300.0]));
        if let Some(max_iterations) = max_iterations {
            system.set_config(SolveConfig {
                max_iterations,
                ..SolveConfig::default()
            });
        }
        (system, body)
    };
    let (mut solved, body) = build(None);
    let report = with_cx(|cx| solved.step(1.0e-3, cx)).unwrap();
    assert_eq!(report.iterations, 4);
    // The step is stiff (dt G / C is 100): the body's mean temperature
    // over it comes within 0.2 K of the far side.
    let mean_k = 0.5 * (300.0 + *solved.state::<f64>(body).unwrap());
    assert!(mean_k > 302.8 && mean_k < 303.0, "{mean_k}");

    let (mut refused, _) = build(Some(3));
    let error = with_cx(|cx| refused.step(1.0e-3, cx)).expect_err("not converged");
    let SystemError::NotConverged {
        iterations,
        worst,
        residual,
        ..
    } = error
    else {
        panic!("expected NotConverged, got {error}");
    };
    assert_eq!((iterations, worst.as_str()), (3, "wall.a"));
    // Open, and inside the first test's 1e-12 of 300 K.
    assert!(residual.abs() > 1.0e-11 && residual.abs() < 2.0e-10, "{residual}");
    assert_eq!(refused.steps(), 0);
}

#[test]
fn what_a_solve_keeps_between_steps_does_not_decide_a_step() {
    // The warming plane settles, and near its steady state steps are taken
    // at their first check, with no Newton update. A copy takes a step of
    // another length and is restored before every step. A checkpoint
    // carries the states, the clock and the accounts, and nothing of the
    // solve's workspace; so if anything the workspace keeps decided whether
    // a step is taken at its first check, the copy would take other
    // iterations and reach other bits.
    let (mut straight, bodies) = warming_plane();
    let (mut detoured, _) = warming_plane();
    let (dt, steps) = (0.25_f64, 80_u32);
    let mut taken_at_first_check = 0_u32;
    with_cx(|cx| {
        for step in 0..steps {
            let report = straight.step(dt, cx).unwrap();
            let before = detoured.checkpoint().unwrap();
            detoured.step(0.37 * dt, cx).unwrap();
            detoured.restore(&before).unwrap();
            let again = detoured.step(dt, cx).unwrap();
            assert_eq!(report.iterations, again.iterations, "step {step}");
            for body in &bodies {
                let (a, b) = (
                    straight.state::<f64>(*body).unwrap(),
                    detoured.state::<f64>(*body).unwrap(),
                );
                assert_eq!(a.to_bits(), b.to_bits(), "step {step}: {a} against {b}");
            }
            taken_at_first_check += u32::from(report.iterations == 0);
        }
    });
    // The comparison covered steps taken at their first check.
    assert!(taken_at_first_check >= 10, "{taken_at_first_check}");
    assert_eq!(straight.time_s().to_bits(), detoured.time_s().to_bits());
    assert_eq!(
        straight.ledger().energy_residual_j().to_bits(),
        detoured.ledger().energy_residual_j().to_bits()
    );
}

/// What [`balance_of_the_warming_plane`] measured over a run: the defect
/// of a body's balance that is largest against its allowance [J], the step
/// it was found at and that allowance [J]; the iterations spent, the
/// steps taken with none, and the bodies' temperatures at the end [K].
struct Balance {
    worst_j: f64,
    worst_step: u32,
    allowed_at_worst_j: f64,
    iterations: u64,
    taken_at_first_check: u32,
    ended_k: [f64; 2],
}

/// Steps [`warming_plane`] `steps` times, step `n` over `length_s(n)`, at
/// the kernel's default solve settings.
///
/// After every step it solves the step's own equations for the states the
/// kernel published and measures, in joules, how far each body is from its
/// balance: with m_k the mean of a body's two end temperatures (what a
/// heat capacity returns on its port) and T the plane temperature that
/// satisfies sum G_k (T - m_k) = P(T),
///   C_k (T_k' - T_k) = dt (G_k (T - m_k) + F_k (300 - m_k)).
/// No wire value of the kernel enters: only the bodies' states. A body is
/// allowed 1e-12 of the energy that crossed the plane in the step, and the
/// rounding of the terms its balance is formed from.
fn balance_of_the_warming_plane(steps: u32, length_s: impl Fn(u32) -> f64) -> Balance {
    let (capacity, plane_w_k, film_w_k) = (WARMING_J_K, WARMING_PLANE_W_K, WARMING_FILM_W_K);
    let (power_w, scale_k, room_k) = (WARMING_POWER_W, WARMING_SCALE_K, WARMING_ROOM_K);
    let law = FadingSource::new(power_w, scale_k);
    let (mut system, bodies) = warming_plane();
    let total_w_k: f64 = plane_w_k.iter().sum();
    let (mut worst_j, mut worst_step, mut allowed_at_worst_j) = (0.0_f64, 0_u32, 0.0_f64);
    let (mut iterations, mut taken_at_first_check) = (0_u64, 0_u32);
    let mut before = [room_k; 2];
    with_cx(|cx| {
        for step in 0..steps {
            let dt = length_s(step);
            let spent = system.step(dt, cx).unwrap().iterations;
            iterations += u64::from(spent);
            taken_at_first_check += u32::from(spent == 0);
            let after = [0, 1].map(|side: usize| *system.state::<f64>(bodies[side]).unwrap());
            let mean = [0, 1].map(|side: usize| 0.5 * (before[side] + after[side]));
            // The plane's balance for these means, by Newton on a function
            // that rises with T: sum G_k (T - m_k) - P(T) = 0.
            let weighted: f64 = (0..2).map(|side| plane_w_k[side] * mean[side]).sum();
            let mut plane_k = weighted / total_w_k;
            for _ in 0..60 {
                let power = law.power_at(plane_k);
                plane_k -= (total_w_k * plane_k - weighted - power) / (total_w_k + power / scale_k);
            }
            let crossed_j = dt * law.power_at(plane_k);
            for side in 0..2 {
                let received_j = dt
                    * (plane_w_k[side] * (plane_k - mean[side])
                        + film_w_k[side] * (room_k - mean[side]));
                let defect_j = capacity[side] * (after[side] - before[side]) - received_j;
                let allowed_j = 1.0e-12 * crossed_j
                    + 64.0
                        * f64::EPSILON
                        * (capacity[side] * after[side]
                            + dt * (plane_w_k[side] + film_w_k[side]) * plane_k);
                if defect_j.abs() / allowed_j > worst_j / allowed_at_worst_j.max(f64::MIN_POSITIVE) {
                    (worst_j, worst_step, allowed_at_worst_j) = (defect_j.abs(), step, allowed_j);
                }
            }
            before = after;
        }
    });
    Balance {
        worst_j,
        worst_step,
        allowed_at_worst_j,
        iterations,
        taken_at_first_check,
        ended_k: before,
    }
}

#[test]
fn a_loop_is_closed_when_the_step_length_changes_at_every_step() {
    // The warming plane from the room's temperature, 240 steps whose length
    // cycles through 0.25, 0.125, 0.0925 and 0.5 s, so no step has the
    // length of the one before it. Each body's balance is measured after
    // every step as in the test below.
    //
    // The slopes that judge a residual are formed in the step that is
    // judged, so a new step length needs nothing from an earlier step. Two
    // things are asked. The balance holds at every step: judged against
    // 300 K alone, a settled step leaves a temperature open and a body
    // about 1e-9 J out here. And a settled step is still taken with no
    // update: a solve that judged by slopes kept from the step before has
    // none at a new length and updates at every one of these steps.
    let steps = 240_u32;
    let balance =
        balance_of_the_warming_plane(steps, |step| 0.25 * [1.0, 0.5, 0.37, 2.0][step as usize % 4]);
    // The run is not trivial: the bodies warmed by most of a kelvin, and
    // the question is asked of steps taken with no update as well as of
    // steps that iterated.
    let Balance {
        worst_j,
        worst_step,
        allowed_at_worst_j,
        iterations,
        taken_at_first_check,
        ended_k,
    } = balance;
    assert!(ended_k[0] > WARMING_ROOM_K + 0.5 && ended_k[1] > WARMING_ROOM_K + 0.5, "{ended_k:?}");
    assert!(taken_at_first_check >= 40 && iterations >= 40, "{taken_at_first_check} {iterations}");
    eprintln!(
        "a loop whose step length changes at every step: a body's balance is off by \
         {worst_j:.3e} J at worst (step {worst_step}, allowed {allowed_at_worst_j:.3e} J); \
         {iterations} iterations in {steps} steps, {taken_at_first_check} steps taken with none"
    );
    assert!(
        worst_j <= allowed_at_worst_j,
        "step {worst_step}: {worst_j:e} J against {allowed_at_worst_j:e} J"
    );
}

#[test]
fn a_block_that_fails_only_the_second_test_names_the_open_input() {
    // Two bodies 1e-10 K apart through a wall of 1e6 W/K. At the start the
    // wall's two temperature inputs hold the guess, 300 K on both sides, so
    // one of them is open by 1e-10 K: inside the first test, which allows
    // 1e-12 of 300 K. Through the wall that is 1e-4 W of heat the wall has
    // not been told of, against bodies that exchange nothing yet. With no
    // iteration allowed the block is refused, and the refusal names the
    // open temperature, not a wire whose own residual is too large (none
    // is).
    let apart_k = 1.0e-10_f64;
    let build = || {
        let mut builder = SystemBuilder::new();
        let port = |builder: &SystemBuilder, node, name: &str| builder.port(node, name).unwrap();
        let wall = builder
            .add_with_tangent("wall", Conductance::new(1.0e6).unwrap(), ())
            .unwrap();
        let mut bodies = Vec::new();
        let mut ends = Vec::new();
        for (name, end, kelvin) in [("cool", "a", 300.0), ("warm", "b", 300.0 + apart_k)] {
            let body = builder
                .add_with_tangent(name, HeatCapacity::new(1.0).unwrap(), kelvin)
                .unwrap();
            let at_wall = port(&builder, wall, end);
            builder
                .common_effort(port(&builder, body, "heat"), &[at_wall])
                .unwrap();
            bodies.push(body);
            ends.push(at_wall);
        }
        let mut system = builder.seal().unwrap();
        for end in ends {
            assert!(system.guess(end, &[300.0]));
        }
        (system, bodies)
    };
    let (mut refused, _) = build();
    refused.set_config(SolveConfig {
        max_iterations: 0,
        ..SolveConfig::default()
    });
    let error = with_cx(|cx| refused.step(1.0e-3, cx)).expect_err("not converged");
    let SystemError::NotConverged {
        iterations,
        worst,
        residual,
        ..
    } = error
    else {
        panic!("expected NotConverged, got {error}");
    };
    assert_eq!((iterations, worst.as_str()), (0, "wall.b"));
    // Its own residual is inside the first test's 1e-12 of 300 K.
    assert!(residual.abs() <= 3.0e-10 && (residual.abs() - apart_k).abs() <= 1.0e-12, "{residual}");
    assert_eq!(refused.steps(), 0);

    // With the default settings the same block is solved: the wall is told
    // the temperatures and the bodies exchange what a millisecond allows.
    let (mut solved, bodies) = build();
    let report = with_cx(|cx| solved.step(1.0e-3, cx)).unwrap();
    assert!(report.iterations >= 1, "{}", report.iterations);
    let gap = solved.state::<f64>(bodies[1]).unwrap() - solved.state::<f64>(bodies[0]).unwrap();
    // The midpoint rule's factor (1 - a) / (1 + a) with a = dt G = 1000:
    // the bodies overshoot to nearly the opposite gap.
    let factor = (1.0 - 1000.0) / (1.0 + 1000.0);
    assert!((gap - factor * apart_k).abs() <= 2.0e-12, "{gap}");
}

#[test]
fn a_loop_is_closed_in_the_heat_its_temperatures_drive_at_the_default_tolerance() {
    // Two bodies of 1 J/K on a plane of 2 and 4 W/K that carries a heat
    // source falling with the plane's temperature, each body also losing
    // heat to a 300 K room through a film (3 and 5 W/K). The kernel's
    // default solve settings, 80 steps of 0.25 s towards the steady state.
    //
    // After every step the balance of each body is measured from the
    // bodies' states alone, by [`balance_of_the_warming_plane`].
    //
    // A loop test that judges a temperature against 300 K lets it stay
    // open by 3e-10 K, which through 4 W/K over 0.25 s is 3e-10 J. Judged
    // by the heat they drive, the residuals leave about 1e-11 J.
    let steps = 80_u32;
    let Balance {
        worst_j,
        worst_step,
        allowed_at_worst_j,
        iterations,
        ended_k,
        ..
    } = balance_of_the_warming_plane(steps, |_| 0.25);
    // The run is not trivial: the bodies warmed by most of a kelvin and
    // the source weakened with them.
    assert!(ended_k[0] > WARMING_ROOM_K + 0.5 && ended_k[1] > WARMING_ROOM_K + 0.5, "{ended_k:?}");
    eprintln!(
        "a loop at the default tolerance: a body's balance is off by {worst_j:.3e} J at worst \
         (step {worst_step}, allowed {allowed_at_worst_j:.3e} J); {iterations} iterations in \
         {steps} steps"
    );
    assert!(
        worst_j <= allowed_at_worst_j,
        "step {worst_step}: {worst_j:e} J against {allowed_at_worst_j:e} J"
    );
}
