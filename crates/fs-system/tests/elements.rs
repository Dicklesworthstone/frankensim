//! The library elements: each passes the conformance kit over a sweep of
//! states and inputs, and devices composed from them reproduce closed-form
//! solutions of the discrete equations they are meant to integrate.

// Exact float comparisons here assert bitwise identities of the wiring.
#![allow(clippy::float_cmp)]

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent, Proposal, Receipt};
use fs_ports::PortKind;
use fs_system::elements::{
    Conductance, Damper, ElementError, Force, Gyrator, HeatCapacity, Heater, Mass, Motion,
    SourceAccount, Spring, Waveform,
};
use fs_system::{
    AdmissionError, NodeId, RouteTemperature, System, SystemBuilder, heat_sink,
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
