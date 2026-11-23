//! The field in a gap of variable width: the conformance kit over a sweep,
//! its efforts against an independent quadrature of the gradient, a held
//! field accelerating its mass uniformly, a gap discharging through a
//! resistance against the hand-solved monolithic discrete-gradient step, and
//! its refusals.
// Exact float comparisons here assert bitwise identities of the wiring.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Component, PortRole, Receipt};
use fs_system::elements::{
    Admittance, Damper, ElementError, Force, GapField, GapFieldState, HeatCapacity, Mass, Motion,
    SourceAccount, Waveform,
};
use fs_system::{NodeId, PortRef, RouteTemperature, System, SystemBuilder};

const ELECTRICAL: PortKind = PortKind::ElectricalVoltageCurrent;
const LIMITS: [f64; 2] = [0.5, 1.5];

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

/// Deterministic numbers in `[lo, hi)` (splitmix64).
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

fn field(inductor: bool, scale: f64) -> GapField {
    if inductor {
        GapField::inductor(scale, LIMITS).unwrap()
    } else {
        GapField::capacitor(scale, LIMITS).unwrap()
    }
}

/// The mandatory and tangent checks at one operating point: the receipt and
/// the outputs of the step.
fn conforms(
    element: &GapField,
    state: &GapFieldState,
    stale: &GapFieldState,
    dt_s: f64,
    inputs: &[f64],
) -> (Receipt, Vec<f64>) {
    let case = Case {
        state,
        stale,
        dt_s,
        inputs,
    };
    with_cx(|cx| {
        let (receipt, outputs) = match check_component(element, &case, cx)
            .unwrap_or_else(|failure| panic!("conformance at {inputs:?}: {failure}"))
        {
            Observed::Stepped {
                receipt, outputs, ..
            } => (receipt, outputs),
            Observed::Shortened { dt_s } => panic!("shortened to {dt_s} at {inputs:?}"),
        };
        check_port_tangent(element, &case, 1.0e-6, 1.0e-12, cx)
            .unwrap_or_else(|failure| panic!("tangent at {inputs:?}: {failure}"));
        (receipt, outputs)
    })
}

#[test]
fn the_field_conforms_and_its_efforts_are_the_mean_gradient_along_the_step() {
    let mut numbers = Numbers(17);
    for case in 0..64 {
        let inductor = case % 2 == 1;
        let s = numbers.next(0.5, 3.0);
        let element = field(inductor, s);
        assert_eq!(
            element.ports()[0].role(),
            if inductor {
                PortRole::EffortIn
            } else {
                PortRole::FlowIn
            }
        );
        let state = GapFieldState {
            coordinate: numbers.next(-2.0, 2.0),
            gap_m: numbers.next(0.75, 1.25),
        };
        let stale = GapFieldState {
            coordinate: numbers.next(-2.0, 2.0),
            gap_m: numbers.next(0.75, 1.25),
        };
        let dt = numbers.next(0.01, 0.2);
        let inputs = [numbers.next(-3.0, 3.0), numbers.next(-1.0, 1.0)];
        let (receipt, outputs) = conforms(&element, &state, &stale, dt, &inputs);

        // The gradient `(x g / s, x^2 / (2 s))` is quadratic along the
        // straight segment, so Simpson's rule gives its mean exactly.
        let (x0, g0) = (state.coordinate, state.gap_m);
        let (x1, g1) = (x0 + dt * inputs[0], g0 + dt * inputs[1]);
        let gradient = |t: f64| {
            let (x, g) = (x0 + t * (x1 - x0), g0 + t * (g1 - g0));
            [x * g / s, 0.5 * x * x / s]
        };
        let (a, m, b) = (gradient(0.0), gradient(0.5), gradient(1.0));
        for k in 0..2 {
            let simpson = (a[k] + 4.0 * m[k] + b[k]) / 6.0;
            let size = a[k].abs() + 4.0 * m[k].abs() + b[k].abs();
            assert!(
                (outputs[k] - simpson).abs() <= 16.0 * f64::EPSILON * size,
                "effort {k}: {} vs {simpson}",
                outputs[k]
            );
        }

        // The step stores exactly the energy difference.
        let energy = |x: f64, g: f64| 0.5 * x * x * g / s;
        let change = energy(x1, g1) - energy(x0, g0);
        let size = energy(x0, g0)
            + energy(x1, g1)
            + (outputs[0] * (x1 - x0)).abs()
            + (outputs[1] * (g1 - g0)).abs();
        assert!(
            (receipt.stored_change_j - change).abs() <= 32.0 * f64::EPSILON * size,
            "{} vs {change}",
            receipt.stored_change_j
        );
        assert_eq!(receipt.dissipated_j, 0.0);
    }
}

/// A field joined to a mass and a held force: the mass sets the velocity at
/// which the gap opens and receives the force less the field's attraction.
struct Device {
    system: System,
    gap: NodeId,
    body: NodeId,
    push: NodeId,
}

fn device(
    element: GapField,
    start: GapFieldState,
    mass_kg: f64,
    velocity_m_s: f64,
    force_n: f64,
    builder: &mut SystemBuilder,
) -> (NodeId, NodeId, NodeId) {
    let gap = builder.add_with_tangent("gap", element, start).unwrap();
    let body = builder
        .add_with_tangent(
            "body",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: start.gap_m,
                velocity_m_s,
            },
        )
        .unwrap();
    let push = builder
        .add_with_tangent(
            "push",
            Force::new(Waveform::constant(force_n).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();
    // The body moves toward a wider gap: the gap receives the body's
    // velocity and the body minus the gap's force.
    builder
        .common_flow_ratios(
            builder.port(body, "p").unwrap(),
            &[
                (builder.port(gap, "gap").unwrap(), -1.0),
                (builder.port(push, "p").unwrap(), 1.0),
            ],
        )
        .unwrap();
    (gap, body, push)
}

fn step(system: &mut System, dt_s: f64) {
    with_cx(|cx| {
        let report = system.step(dt_s, cx).expect("step");
        assert_eq!(report.dt_s.to_bits(), dt_s.to_bits());
    });
}

#[test]
fn a_held_field_accelerates_its_mass_uniformly_and_keeps_its_coordinate() {
    let (s, x0, g0, m, v0, f, dt) = (2.0, 0.8, 1.0, 3.0, 0.02, 0.1, 0.05);
    for inductor in [false, true] {
        let mut builder = SystemBuilder::new();
        let start = GapFieldState {
            coordinate: x0,
            gap_m: g0,
        };
        let (gap, body, push) = device(field(inductor, s), start, m, v0, f, &mut builder);
        // No current into the plates, or no voltage across the winding.
        builder
            .boundary(
                builder.port(gap, "field").unwrap(),
                Box::new(|_, _, values| values[0] = 0.0),
            )
            .unwrap();
        let mut device = Device {
            system: builder.seal().unwrap(),
            gap,
            body,
            push,
        };
        // A constant force: the midpoint rule integrates it exactly.
        let a = (f - 0.5 * x0 * x0 / s) / m;
        for n in 1..=10 {
            step(&mut device.system, dt);
            let t = f64::from(n) * dt;
            let held = *device.system.state::<GapFieldState>(device.gap).unwrap();
            let motion = *device.system.state::<Motion>(device.body).unwrap();
            assert_eq!(held.coordinate, x0, "the coordinate is held");
            assert_eq!(held.gap_m, motion.position_m, "the gap moves with the body");
            let gap = g0 + v0 * t + 0.5 * a * t * t;
            assert!(
                (held.gap_m - gap).abs() <= 1.0e-14,
                "{} vs {gap}",
                held.gap_m
            );
            let velocity = v0 + a * t;
            assert!(
                (motion.velocity_m_s - velocity).abs() <= 1.0e-15,
                "{} vs {velocity}",
                motion.velocity_m_s
            );
        }
        // The force's work is the field's and the body's energy change.
        let held = *device.system.state::<GapFieldState>(device.gap).unwrap();
        let motion = *device.system.state::<Motion>(device.body).unwrap();
        let work = -device
            .system
            .state::<SourceAccount>(device.push)
            .unwrap()
            .energy_in_j;
        let stored = 0.5 * x0 * x0 * (held.gap_m - g0) / s
            + 0.5 * m * (motion.velocity_m_s * motion.velocity_m_s - v0 * v0);
        assert!((work - f * (held.gap_m - g0)).abs() <= 1.0e-16, "{work}");
        assert!((work - stored).abs() <= 1.0e-16, "{work} vs {stored}");
        assert!(device.system.ledger().energy_residual_j().abs() <= 1.0e-16);
    }
}

/// One step of the monolithic discrete-gradient rule for a field whose
/// coordinate decays at `loss` times its effort (an inductor shorted by a
/// resistance `R`, `dx/dt = -R i`, or a capacitor across a conductance `G`,
/// `dx/dt = -G v`), its gap opened by a mass under a held force. Newton on
/// the end coordinate and gap; returns them and the step's heat.
fn monolithic_step(
    [x0, g0, v0]: [f64; 3],
    s: f64,
    loss: f64,
    m: f64,
    f: f64,
    dt: f64,
) -> ([f64; 3], f64) {
    let (mut x1, mut g1) = (x0, g0 + dt * v0);
    for _ in 0..50 {
        let (xm, gm) = (0.5 * (x0 + x1), 0.5 * (g0 + g1));
        let (dx, dg) = (x1 - x0, g1 - g0);
        let effort = xm * gm / s + dx * dg / (12.0 * s);
        let force = 0.5 * (xm * xm + dx * dx / 12.0) / s;
        let v1 = 2.0 * dg / dt - v0;
        let r = [dx + dt * loss * effort, m * (v1 - v0) - dt * (f - force)];
        let by_coordinate = (0.5 * gm + dg / 12.0) / s;
        let mixed = (0.5 * xm + dx / 12.0) / s;
        let j = [
            [1.0 + dt * loss * by_coordinate, dt * loss * mixed],
            [dt * mixed, 2.0 * m / dt],
        ];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let step_x = (r[0] * j[1][1] - j[0][1] * r[1]) / det;
        let step_g = (j[0][0] * r[1] - j[1][0] * r[0]) / det;
        x1 -= step_x;
        g1 -= step_g;
        if step_x.abs() <= 1.0e-17 && step_g.abs() <= 1.0e-17 {
            break;
        }
    }
    let (xm, gm) = (0.5 * (x0 + x1), 0.5 * (g0 + g1));
    let effort = xm * gm / s + (x1 - x0) * (g1 - g0) / (12.0 * s);
    (
        [x1, g1, 2.0 * (g1 - g0) / dt - v0],
        dt * loss * effort * effort,
    )
}

/// A resistance across the field's port, its far end held and its heat
/// routed into a heat capacity `case`: an inductor's winding returns a
/// voltage for the current it receives, a capacitor's plates draw a current
/// for their voltage. Returns the resistance's heat port and the rate at
/// which the coordinate decays per unit effort, `R` or `G`.
fn resistance(
    builder: &mut SystemBuilder,
    gap: NodeId,
    inductor: bool,
    ohms: f64,
) -> (PortRef, f64) {
    let case = builder
        .add_with_tangent("case", HeatCapacity::new(1.0e3).unwrap(), 300.0)
        .unwrap();
    let field = builder.port(gap, "field").unwrap();
    let (node, loss) = if inductor {
        let node = builder
            .add_with_tangent(
                "resistance",
                Damper::in_domain(ELECTRICAL, ohms).unwrap(),
                (),
            )
            .unwrap();
        builder
            .bond(builder.port(node, "a").unwrap(), field)
            .unwrap();
        (node, ohms)
    } else {
        let conductance = Admittance::in_domain(ELECTRICAL, 1.0 / ohms).unwrap();
        let node = builder
            .add_with_tangent("resistance", conductance, ())
            .unwrap();
        builder
            .bond(field, builder.port(node, "a").unwrap())
            .unwrap();
        (node, 1.0 / ohms)
    };
    builder
        .boundary(
            builder.port(node, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    let heat = builder.port(node, "heat").unwrap();
    builder
        .heat_route(
            heat,
            &[(builder.port(case, "heat").unwrap(), 1.0)],
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
    (heat, loss)
}

#[test]
fn a_gap_discharging_through_a_resistance_takes_the_monolithic_discrete_gradient_step() {
    let (s, x0, g0, m, v0, f, dt, ohms) = (1.5, 1.2, 1.0, 2.0, 0.0, 0.05, 0.05, 0.8);
    for inductor in [false, true] {
        let mut builder = SystemBuilder::new();
        let start = GapFieldState {
            coordinate: x0,
            gap_m: g0,
        };
        let (gap, body, push) = device(field(inductor, s), start, m, v0, f, &mut builder);
        let (heat, loss) = resistance(&mut builder, gap, inductor, ohms);
        let mut device = Device {
            system: builder.seal().unwrap(),
            gap,
            body,
            push,
        };
        assert!(device.system.guess(heat, &[300.0]));

        let mut hand = [x0, g0, v0];
        let mut heat_j = 0.0;
        let (mut worst_x, mut worst_g, mut worst_v, mut worst_p) =
            (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..20 {
            step(&mut device.system, dt);
            let (next, joule) = monolithic_step(hand, s, loss, m, f, dt);
            hand = next;
            heat_j += joule;
            let held = *device.system.state::<GapFieldState>(device.gap).unwrap();
            let motion = *device.system.state::<Motion>(device.body).unwrap();
            worst_p = worst_p.max((held.gap_m - motion.position_m).abs());
            worst_x = worst_x.max((held.coordinate - hand[0]).abs());
            worst_g = worst_g.max((held.gap_m - hand[1]).abs());
            worst_v = worst_v.max((motion.velocity_m_s - hand[2]).abs());
        }
        let kind = if inductor { "inductor" } else { "capacitor" };
        // The kernel's interface solve meets the monolithic step to about
        // 1e-15 here, the body's position the gap to an ulp.
        assert!(worst_x <= 1.0e-13, "{kind}: coordinate {worst_x:e}");
        assert!(worst_g <= 1.0e-13, "{kind}: gap {worst_g:e}");
        assert!(worst_v <= 1.0e-13, "{kind}: velocity {worst_v:e}");
        assert!(worst_p <= 1.0e-14, "{kind}: position {worst_p:e}");
        // The field lost over a third of its coordinate (time constants of
        // 1.2 s and 1.9 s), and the gap closed.
        assert!(hand[0].abs() < 0.65 * x0 && hand[1] < g0 - 0.05, "{hand:?}");

        // All the resistance's heat reached the case; the force's work is the
        // field's, the body's and the case's energy change.
        let ledger = device.system.ledger();
        let account = |name: &str| ledger.nodes.iter().find(|node| node.name == name).unwrap();
        let warmed = account("case").stored_change_j.value();
        assert!(
            (warmed - heat_j).abs() <= 2.0e-13,
            "{kind}: {warmed} vs {heat_j}"
        );
        let held = *device.system.state::<GapFieldState>(device.gap).unwrap();
        let motion = *device.system.state::<Motion>(device.body).unwrap();
        let work = -device
            .system
            .state::<SourceAccount>(device.push)
            .unwrap()
            .energy_in_j;
        let stored = 0.5 * held.coordinate * held.coordinate * held.gap_m / s
            - 0.5 * x0 * x0 * g0 / s
            + 0.5 * m * (motion.velocity_m_s * motion.velocity_m_s - v0 * v0);
        assert!(
            (work - stored - warmed).abs() <= 1.0e-12,
            "{kind}: {work} vs {stored} + {warmed}"
        );
        assert!(ledger.energy_residual_j().abs() <= 1.0e-15);
    }
}

#[test]
fn refusals_leave_nothing_behind() {
    // Construction.
    for (scale, limits) in [
        (0.0, LIMITS),
        (-1.0, LIMITS),
        (f64::NAN, LIMITS),
        (f64::INFINITY, LIMITS),
        // Its reciprocal overflows.
        (5.0e-324, LIMITS),
        (1.0, [-0.1, 1.5]),
        (1.0, [f64::NAN, 1.5]),
        (1.0, [0.5, 0.5]),
        (1.0, [0.5, f64::INFINITY]),
    ] {
        for inductor in [false, true] {
            let built = if inductor {
                GapField::inductor(scale, limits)
            } else {
                GapField::capacitor(scale, limits)
            };
            assert!(
                matches!(built, Err(ElementError::Parameter { .. })),
                "{scale} {limits:?}: {built:?}"
            );
        }
    }
    let element = field(false, 2.0);
    assert_eq!(element.coupling_scale(), 2.0);
    assert_eq!(element.gap_limits_m(), LIMITS);
    assert!(!element.is_inductor() && field(true, 2.0).is_inductor());

    // States outside the open range have no energy and no codec.
    for gap_m in [0.5, 1.5, 0.2, f64::NAN] {
        let state = GapFieldState {
            coordinate: 1.0,
            gap_m,
        };
        assert!(matches!(
            element.stored_energy_j(&state),
            Err(ElementError::Law { .. })
        ));
        assert!(!element.encode_state(&state, &mut Vec::new()));
        assert_eq!(
            element.decode_state(&[1.0_f64.to_bits(), gap_m.to_bits()]),
            None
        );
    }
    let state = GapFieldState {
        coordinate: -0.75,
        gap_m: 1.25,
    };
    let mut words = Vec::new();
    assert!(element.encode_state(&state, &mut words));
    assert_eq!(element.decode_state(&words), Some(state));

    // A step whose end leaves the range is refused; so is a non-finite input.
    with_cx(|cx| {
        let mut new = state;
        let mut outputs = [7.0; 2];
        for inputs in [[0.0, 2.6], [0.0, -15.0], [f64::NAN, 0.0]] {
            let refused = element.propose(&state, 0.1, &inputs, &mut new, &mut outputs, cx);
            assert!(
                matches!(
                    refused,
                    Err(ElementError::Law { .. } | ElementError::Input { .. })
                ),
                "{inputs:?}: {refused:?}"
            );
            assert_eq!(new, state);
            assert_eq!(outputs, [7.0; 2]);
        }
    });

    // The kind, the scale and the limits each change the identity.
    let words = |element: &GapField| {
        let mut words = Vec::new();
        element.identity_words(&mut words);
        words
    };
    let base = words(&element);
    assert_ne!(base, words(&field(true, 2.0)));
    assert_ne!(base, words(&field(false, 2.5)));
    assert_ne!(base, words(&GapField::capacitor(2.0, [0.4, 1.5]).unwrap()));
    assert_ne!(base, words(&GapField::capacitor(2.0, [0.5, 1.6]).unwrap()));
}
