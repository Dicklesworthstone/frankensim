//! The library elements: each passes the conformance kit over a sweep of
//! states and inputs, and devices composed from them reproduce closed-form
//! solutions of the discrete equations they are meant to integrate.

// Exact float comparisons here assert bitwise identities of the wiring.
#![allow(clippy::float_cmp)]

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortTangent, Proposal};
use fs_system::elements::{
    Conductance, Damper, ElementError, Force, HeatCapacity, Heater, Mass, Motion, SourceAccount,
    Spring, Waveform,
};
use fs_system::{NodeId, System, SystemBuilder};

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
