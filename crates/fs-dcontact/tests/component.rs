//! `ContactInterface`, the normal power-law contact on the composition
//! contract: a mass bouncing on a Hertz wall against the closed-form impact,
//! the reaction a load line settles on against the linear law's closed form,
//! the conformance kit, and the refusals.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
use fs_dcontact::{DContactError, Obstacle};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Observed, check_component, check_load_line};
use fs_ports::component::{Class, Component, LoadLine, Proposal};
use fs_system::elements::{Mass, Motion};
use fs_system::{NodeId, System, SystemBuilder};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 7,
            kernel_id: 8,
            tile: 9,
            iteration: 10,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// The live drum's Hertz tip, `K = 4/3 E* sqrt(R)` with `E* = 0.8 GPa` and
/// `R = 3 mm` (engine-drum's instrument card), on a one-point unit-opening
/// obstacle with no gap.
const HERTZ_K: f64 = 4.0 / 3.0 * 0.8e9 * 0.054_772_255_750_516_61;

fn law(stiffness: f64, alpha: f64, gap: f64) -> Obstacle {
    Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![gap],
        vec![1.0],
        stiffness,
        alpha,
        "fs-dcontact component test: authored power law".into(),
    )
    .unwrap()
}

fn solve() -> ContactSolve {
    ContactSolve {
        max_iterations: 100,
        maximum_force_n: 1.0e6,
        maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1.0e-9,
        force_relative_tolerance: 1.0e-12,
    }
}

fn contact(stiffness: f64, alpha: f64, gap: f64) -> ContactInterface {
    ContactInterface::new(law(stiffness, alpha, gap), solve()).unwrap()
}

// ---- A mass on a Hertz wall ----------------------------------------------

const MASS_KG: f64 = 0.01;
const SPEED_M_S: f64 = 2.0;
/// The mass starts this far from the wall [m].
const APPROACH_M: f64 = 1.0e-5;

/// A mass moving towards a wall at `SPEED_M_S`, the contact between them:
/// end `a` on the mass (it receives minus the mass's velocity, so a mass
/// moving the negative way closes the contact), end `b` held at rest.
fn mass_on_a_wall() -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let mass = builder
        .add_with_tangent(
            "mass",
            Mass::new(MASS_KG).unwrap(),
            Motion {
                position_m: APPROACH_M,
                velocity_m_s: -SPEED_M_S,
            },
        )
        .unwrap();
    let wall = builder
        .add_on_line(
            "wall",
            contact(HERTZ_K, 1.5, 0.0),
            ContactState::at(-APPROACH_M),
        )
        .unwrap();
    builder
        .bond(
            builder.port(wall, "a").unwrap(),
            builder.port(mass, "p").unwrap(),
        )
        .unwrap();
    builder
        .boundary(
            builder.port(wall, "b").unwrap(),
            Box::new(|_, _, values| values[0] = 0.0),
        )
        .unwrap();
    (builder.seal().unwrap(), mass, wall)
}

/// One bounce: when the contact began and ended [s] (the closure's
/// crossings of zero, interpolated in the step), its deepest sampled
/// closure [m], and the mass's velocity after it [m/s].
struct Bounce {
    entry_s: f64,
    exit_s: f64,
    peak_m: f64,
    exit_m_s: f64,
}

fn bounce(dt_s: f64, cx: &Cx<'_>) -> Bounce {
    let (mut system, mass, wall) = mass_on_a_wall();
    let closure = |system: &System| system.state::<ContactState>(wall).unwrap().closure_m;
    let (mut entry_s, mut exit_s, mut peak_m) = (None, None, 0.0_f64);
    let mut before = closure(&system);
    let mut t = 0.0;
    for _ in 0..1_000_000 {
        system.step(dt_s, cx).unwrap();
        let after = closure(&system);
        let crossing = t + dt_s * (0.0 - before) / (after - before);
        if before <= 0.0 && after > 0.0 && entry_s.is_none() {
            entry_s = Some(crossing);
        }
        if before > 0.0 && after <= 0.0 && entry_s.is_some() {
            exit_s = Some(crossing);
            break;
        }
        peak_m = peak_m.max(after);
        before = after;
        t += dt_s;
    }
    // A few steps apart again: the mass moves freely away.
    for _ in 0..10 {
        system.step(dt_s, cx).unwrap();
    }
    Bounce {
        entry_s: entry_s.expect("the mass reaches the wall"),
        exit_s: exit_s.expect("the mass leaves the wall"),
        peak_m,
        exit_m_s: system.state::<Motion>(mass).unwrap().velocity_m_s,
    }
}

/// `I = integral from 0 to 1 of (1 - s^(alpha + 1))^(-1/2) ds`, by Simpson's
/// rule after `s = 1 - u^2`, which takes the end singularity away. The
/// contact lasts `2 I p_max / v`.
fn duration_integral(alpha: f64) -> f64 {
    let beta = alpha + 1.0;
    let f = |u: f64| {
        if u == 0.0 {
            2.0 / beta.sqrt()
        } else {
            2.0 * u / (1.0 - (1.0 - u * u).powf(beta)).sqrt()
        }
    };
    let n = 200_000;
    let h = 1.0 / f64::from(n);
    let mut sum = f(0.0) + f(1.0);
    for i in 1..n {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * f(f64::from(i) * h);
    }
    sum * h / 3.0
}

#[test]
fn a_mass_on_a_hertz_wall_bounces_at_the_closed_form_depth_time_and_speed() {
    // Energy: m v^2 / 2 = K p^(alpha + 1) / (alpha + 1) at the deepest point.
    let alpha = 1.5;
    let peak = ((alpha + 1.0) * MASS_KG * SPEED_M_S * SPEED_M_S / (2.0 * HERTZ_K))
        .powf(1.0 / (alpha + 1.0));
    let integral = duration_integral(alpha);
    // (2/5) B(2/5, 1/2) for alpha = 3/2.
    assert!((integral - 1.471_639).abs() < 1e-5, "{integral}");
    let duration = 2.0 * integral * peak / SPEED_M_S;
    // The reaction is settled to its tolerance at most, at the deepest
    // force; over the travel in and out that is the most work the settled
    // reactions can do beyond the secants.
    let force = HERTZ_K * peak.powf(alpha);
    let defect_j = (solve().force_absolute_tolerance_n + solve().force_relative_tolerance * force)
        * 2.0
        * peak;
    let kinetic_j = MASS_KG * SPEED_M_S * SPEED_M_S / 2.0;
    with_cx(|cx| {
        let mut errors = Vec::new();
        for steps in [400.0, 800.0] {
            let b = bounce(duration / steps, cx);
            // The contact's secant and the mass's midpoint rule conserve the
            // energy, so the mass leaves at the speed it came with, to the
            // reaction's tolerance (i429: 3.6e-13 and 2.5e-13). The 1e-12
            // beside it, here and on the depth, is for rounding: above 2n
            // ulps of a unit ratio (1.8e-13 at n = 800).
            let speed = b.exit_m_s / SPEED_M_S - 1.0;
            assert!(
                speed.abs() <= defect_j / (2.0 * kinetic_j) + 1e-12,
                "{}",
                b.exit_m_s
            );
            // A step lands within half a step of the turning point, where the
            // closure lies below the peak by at most (alpha + 1) I^2 / (4 n^2)
            // of it for n steps in the contact; it never lies above (i429:
            // -1.96e-6 and -1.98e-6 against bounds of 8.5e-6 and 2.1e-6).
            let depth = b.peak_m / peak - 1.0;
            let sampling = (alpha + 1.0) * integral * integral / (4.0 * steps * steps);
            assert!(
                depth <= 1e-12 && depth >= -sampling,
                "{} vs {peak}: {depth:e}",
                b.peak_m
            );
            // The contact time, from the closure's interpolated crossings of
            // zero (i429: 6.4e-6 and 1.6e-6, a coefficient of 1.02 on 1/n^2).
            // The cap 4/n^2 is not derived: it is about four times that
            // coefficient. The order below is the check of the rule.
            let error = (b.exit_s - b.entry_s) / duration - 1.0;
            assert!(error.abs() <= 4.0 / (steps * steps), "{error:e}");
            errors.push(error.abs());
        }
        // Second order in the step (i429: 1.99995).
        let order = (errors[0] / errors[1]).log2();
        assert!((1.8..=2.2).contains(&order), "{errors:?}: order {order}");
    });
}

// ---- The reaction on a line ----------------------------------------------

fn stepped(observed: Observed) -> (fs_ports::component::Receipt, Vec<f64>) {
    match observed {
        Observed::Stepped {
            receipt, outputs, ..
        } => (receipt, outputs),
        Observed::Shortened { dt_s } => panic!("shortened to {dt_s}"),
    }
}

#[test]
fn a_linear_law_settles_on_its_line_where_the_closed_form_puts_it() {
    // alpha = 1, both ends of the move inside the contact: the secant is
    // K (p0 + p1) / 2 with p1 = free - gap - S R, so
    // R = K (p0 + free - gap) / (2 + K S).
    let (stiffness, gap) = (2.0e5, 1.0e-4);
    let interface = contact(stiffness, 1.0, gap);
    let dt = 1.0e-4;
    let state = ContactState::at(3.0e-4);
    let stale = ContactState::at(-1.0);
    let (flow_a, flow_b) = (0.3, -0.1);
    let (c_aa, c_ab, c_ba, c_bb) = (2.0e-3, 4.0e-4, 3.0e-4, 1.0e-3);
    let slope = dt * ((c_aa - c_ab) - (c_ba - c_bb));
    let free = 3.0e-4 + dt * (flow_a - flow_b);
    let expected = stiffness * ((3.0e-4 - gap) + (free - gap)) / (2.0 + stiffness * slope);
    with_cx(|cx| {
        let inputs = [flow_a, flow_b];
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: dt,
            inputs: &inputs,
        };
        let line = [c_aa, c_ab, c_ba, c_bb];
        check_load_line(&interface, &case, &line, cx).unwrap();
        let mut new = stale;
        let mut outputs = [0.0; 2];
        let Proposal::Step(receipt) = interface
            .propose_on_line(&state, dt, &inputs, &line, &mut new, &mut outputs, cx)
            .unwrap()
        else {
            panic!("a contact step is never shortened");
        };
        let tolerance = 1.0e-9 + 1.0e-12 * expected;
        assert!(
            (outputs[0] - expected).abs() <= 2.0 * tolerance,
            "{} vs {expected}",
            outputs[0]
        );
        // To the bit: the outputs are the reaction and its negative, and the
        // closure is the line's at that reaction, rounded as the search does.
        assert_eq!(outputs[1].to_bits(), (-outputs[0]).to_bits());
        assert_eq!(new.reaction_n.to_bits(), outputs[0].to_bits());
        assert_eq!(
            new.closure_m.to_bits(),
            (free - slope * outputs[0]).to_bits()
        );
        assert!(receipt.residual <= receipt.tolerance);
        // The secant's work is the potential's change.
        let potential = |x: f64| stiffness * (x - gap).max(0.0).powi(2) / 2.0;
        let change = potential(new.closure_m) - potential(state.closure_m);
        assert!((receipt.stored_change_j - change).abs() <= 1e-12 * change.abs());
        assert_eq!(receipt.dissipated_j.to_bits(), 0.0_f64.to_bits());
    });
}

#[test]
fn the_conformance_kit_passes_apart_entering_and_inside_the_contact() {
    let interface = contact(HERTZ_K, 1.5, 0.0);
    assert_eq!(interface.class(), Class::Monotone1d);
    with_cx(|cx| {
        let dt = 1.0e-6;
        let stale = ContactState::at(7.0);
        for (closure, flow_a) in [
            (-1.0e-3, 1.0),
            (-1.0e-6, 2.0),
            (1.0e-4, 0.5),
            (1.0e-4, -0.5),
        ] {
            let state = ContactState::at(closure);
            let inputs = [flow_a, 0.0];
            let case = Case {
                state: &state,
                stale: &stale,
                dt_s: dt,
                inputs: &inputs,
            };
            let (receipt, outputs) = stepped(check_component(&interface, &case, cx).unwrap());
            assert!(
                outputs[0] >= 0.0 && outputs[1].to_bits() == (-outputs[0]).to_bits(),
                "{outputs:?}"
            );
            assert!(receipt.residual <= receipt.tolerance);
            // A mass of 10 g on the line: c = dt / (2 m).
            let line = [dt / (2.0 * MASS_KG), 0.0, 0.0, 0.0];
            check_load_line(&interface, &case, &line, cx).unwrap();
        }
        // Apart and staying apart: no reaction, nothing stored.
        let state = ContactState::at(-1.0e-3);
        let inputs = [1.0, 0.0];
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: dt,
            inputs: &inputs,
        };
        let (receipt, outputs) = stepped(check_component(&interface, &case, cx).unwrap());
        let bits: Vec<u64> = outputs.iter().map(|x| x.to_bits()).collect();
        assert_eq!(bits, [0.0_f64.to_bits(), (-0.0_f64).to_bits()]);
        assert_eq!(receipt.stored_change_j.to_bits(), 0.0_f64.to_bits());
    });
}

#[test]
fn end_inputs_join_told_neighbours_and_do_not_change_the_step() {
    let plain = contact(HERTZ_K, 1.5, 0.0);
    let told = contact(HERTZ_K, 1.5, 0.0)
        .reading_ends([true, true])
        .unwrap();
    let ids: Vec<_> = told
        .ports()
        .iter()
        .map(|port| port.id().as_str().to_owned())
        .collect();
    assert_eq!(ids, ["a", "b", "a.end", "b.end"]);
    assert_eq!(told.line_ends(), [Some(2), Some(3)]);
    let one = contact(HERTZ_K, 1.5, 0.0)
        .reading_ends([false, true])
        .unwrap();
    assert_eq!(one.line_ends(), [None, Some(2)]);
    assert!(matches!(
        contact(HERTZ_K, 1.5, 0.0).reading_ends([false, false]),
        Err(DContactError::Parameter { .. })
    ));
    assert!(matches!(
        told.clone().reading_ends([true, false]),
        Err(DContactError::Parameter { .. })
    ));
    with_cx(|cx| {
        let dt = 1.0e-6;
        let state = ContactState::at(1.0e-4);
        let line = [5.0e-4, 0.0, 0.0, 2.0e-4];
        let (mut a, mut b) = (state, state);
        let (mut out_a, mut out_b) = ([0.0; 2], [0.0; 2]);
        let plain_receipt = plain
            .propose_on_line(&state, dt, &[0.5, -0.2], &line, &mut a, &mut out_a, cx)
            .unwrap();
        let mut told_line = line.to_vec();
        told_line.extend([1.0e-3, 0.0, 0.0, 4.0e-4]);
        let told_receipt = told
            .propose_on_line(
                &state,
                dt,
                &[0.5, -0.2, 0.7, -0.1],
                &told_line,
                &mut b,
                &mut out_b,
                cx,
            )
            .unwrap();
        assert_eq!(plain_receipt, told_receipt);
        assert_eq!((a, out_a), (b, out_b));
    });
}

#[test]
fn contact_interfaces_refuse_shapes_losses_and_budgets_by_name() {
    let shape = |law: Obstacle| {
        matches!(
            ContactInterface::new(law, solve()),
            Err(DContactError::Shape { .. })
        )
    };
    let parameter = |law: Obstacle, solve: ContactSolve| {
        matches!(
            ContactInterface::new(law, solve),
            Err(DContactError::Parameter { .. })
        )
    };
    // Two points, or a point on another map.
    let two = Obstacle::new(
        vec![-1.0, -1.0],
        2,
        1,
        vec![0.0; 2],
        vec![1.0; 2],
        1e8,
        1.5,
        "two".into(),
    )
    .unwrap();
    assert!(shape(two));
    let flipped = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        1e8,
        1.5,
        "flipped".into(),
    )
    .unwrap();
    assert!(shape(flipped));
    // Contact loss needs a heat port.
    assert!(parameter(
        law(1e8, 1.5, 0.0).with_internal_loss(0.2).unwrap(),
        solve()
    ));
    // Budgets.
    for bad in [
        ContactSolve {
            max_iterations: 0,
            ..solve()
        },
        ContactSolve {
            max_iterations: 129,
            ..solve()
        },
        ContactSolve {
            maximum_force_n: 0.0,
            ..solve()
        },
        ContactSolve {
            maximum_penetration_m: f64::INFINITY,
            ..solve()
        },
        ContactSolve {
            force_absolute_tolerance_n: -1.0,
            ..solve()
        },
        ContactSolve {
            force_relative_tolerance: 1.0,
            ..solve()
        },
    ] {
        assert!(parameter(law(1e8, 1.5, 0.0), bad), "{bad:?}");
    }
}

#[test]
fn contact_interface_steps_refuse_what_they_cannot_settle_by_name() {
    with_cx(|cx| {
        let dt = 1.0e-6;
        let mut new = ContactState::at(0.0);
        let mut outputs = [0.0; 2];
        // A reaction above the ceiling.
        let weak = ContactInterface::new(
            law(HERTZ_K, 1.5, 0.0),
            ContactSolve {
                maximum_force_n: 1.0,
                ..solve()
            },
        )
        .unwrap();
        let refused = weak.propose(
            &ContactState::at(1.0e-4),
            dt,
            &[1.0, 0.0],
            &mut new,
            &mut outputs,
            cx,
        );
        assert!(
            matches!(refused, Err(DContactError::Reaction { what, .. }) if what.contains("maximum_force_n")),
            "{refused:?}"
        );
        // A step that ends deeper than the penetration ceiling.
        let shallow = ContactInterface::new(
            law(HERTZ_K, 1.5, 0.0),
            ContactSolve {
                maximum_penetration_m: 1.0e-5,
                ..solve()
            },
        )
        .unwrap();
        let refused = shallow.propose(
            &ContactState::at(9.9e-6),
            dt,
            &[1.0, 0.0],
            &mut new,
            &mut outputs,
            cx,
        );
        assert!(
            matches!(refused, Err(DContactError::Reaction { what, .. }) if what.contains("maximum_penetration_m")),
            "{refused:?}"
        );
        // Too few iterations for the tolerance.
        let hurried = ContactInterface::new(
            law(HERTZ_K, 1.5, 0.0),
            ContactSolve {
                max_iterations: 2,
                ..solve()
            },
        )
        .unwrap();
        let refused = hurried.propose_on_line(
            &ContactState::at(1.0e-4),
            dt,
            &[1.0, 0.0],
            &[5.0e-4, 0.0, 0.0, 0.0],
            &mut new,
            &mut outputs,
            cx,
        );
        assert!(
            matches!(refused, Err(DContactError::Reaction { what, .. }) if what.contains("max_iterations")),
            "{refused:?}"
        );
        // Wrong lengths and nonfinite inputs.
        let interface = contact(HERTZ_K, 1.5, 0.0);
        let state = ContactState::at(0.0);
        assert!(matches!(
            interface.propose(&state, dt, &[1.0], &mut new, &mut outputs, cx),
            Err(DContactError::Shape { .. })
        ));
        assert!(matches!(
            interface.propose(&state, dt, &[f64::NAN, 0.0], &mut new, &mut outputs, cx),
            Err(DContactError::Parameter { .. })
        ));
        assert!(matches!(
            interface.propose(&state, 0.0, &[1.0, 0.0], &mut new, &mut outputs, cx),
            Err(DContactError::Parameter { .. })
        ));
    });
}
