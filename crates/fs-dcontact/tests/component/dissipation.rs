//! G0/G1 checks of normal-impact loss on the existing composition kernel.
//!
//! The scalar test uses the independently solved quadratic for a linear
//! contact. The isolated two-body test checks the exact Hunt–Crossley
//! restitution relation, momentum and physical thermal energy.

use super::*;
use fs_system::elements::HeatCapacity;

fn lossy_contact(stiffness: f64, alpha: f64, chi: f64) -> ContactInterface {
    ContactInterface::new(
        law(stiffness, alpha, 0.0).with_internal_loss(chi).unwrap(),
        solve(),
    )
    .unwrap()
}

#[test]
fn lossy_linear_contact_settles_at_the_closed_form_load_line_root() {
    let stiffness = 2.0e5_f64;
    let chi = 0.5;
    let dt = 1.0e-3;
    let old = ContactState::at(2.0e-3);
    let stale = ContactState::at(-1.0);
    let inputs = [0.3, 0.1, 300.0];
    let line = [2.0e-3, 0.0, 0.0, 4.0e-3];
    let interface = lossy_contact(stiffness, 1.0, chi);

    // While both endpoints are in contact:
    // E(R) = E0 - a R, 1 + chi*u(R) = c - d R.
    // R = (E0 - a R)(c - d R), with the smaller root selected.
    let slope = dt * (line[0] + line[3]);
    let free = old.closure_m + dt * (inputs[0] - inputs[1]);
    let e0 = 0.5 * stiffness * (old.closure_m + free);
    let a = 0.5 * stiffness * slope;
    let c = 1.0 + chi * (inputs[0] - inputs[1]);
    let d = chi * slope / dt;
    let linear = 1.0 + a * c + e0 * d;
    let discriminant = linear * linear - 4.0 * a * d * e0 * c;
    let expected = 2.0 * e0 * c / (linear + discriminant.sqrt());

    with_cx(|cx| {
        let case = Case {
            state: &old,
            stale: &stale,
            dt_s: dt,
            inputs: &inputs,
        };
        check_load_line(&interface, &case, &line, cx).unwrap();
        let mut next = stale;
        let mut outputs = [0.0; 3];
        let Proposal::Step(receipt) = interface
            .propose_on_line(&old, dt, &inputs, &line, &mut next, &mut outputs, cx)
            .unwrap()
        else {
            panic!("the admitted contact step must settle");
        };

        assert!((outputs[0] - expected).abs() <= 1.0e-8);
        assert_eq!(outputs[1].to_bits(), (-outputs[0]).to_bits());
        assert!(next.closure_m > 0.0);

        let moved = next.closure_m - old.closure_m;
        let velocity = moved / dt;
        let elastic = 0.5 * stiffness * (old.closure_m + next.closure_m);
        assert!(1.0 + chi * velocity > 0.0);
        let expected_heat = dt * chi * elastic * velocity * velocity;
        let expected_storage = 0.5 * stiffness * (next.closure_m.powi(2) - old.closure_m.powi(2));

        assert!((-dt * outputs[2] - expected_heat).abs() <= 1.0e-11 * expected_heat);
        assert!((receipt.dissipated_j - expected_heat).abs() <= 1.0e-11 * expected_heat);
        assert!((receipt.stored_change_j - expected_storage).abs() <= 1.0e-12);
        assert!(
            (receipt.entropy_generated_j_per_k - expected_heat / inputs[2]).abs()
                <= 1.0e-11 * expected_heat / inputs[2]
        );
        assert!(receipt.residual <= receipt.tolerance);
    });
}

#[test]
fn weak_loss_survives_force_rounding_and_fast_unloading_is_nonadhesive() {
    with_cx(|cx| {
        let stiffness = 2.0e5;
        let dt = 1.0e-3;
        let old = ContactState::at(2.0e-3);
        let stale = ContactState::at(-1.0);
        let weak = lossy_contact(stiffness, 1.0, 1.0e-20);

        for velocity in [-0.1, 0.1] {
            let inputs = [velocity, 0.0, 300.0];
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s: dt,
                inputs: &inputs,
            };
            let (receipt, outputs) = stepped(check_component(&weak, &case, cx).unwrap());
            let next_closure = old.closure_m + dt * velocity;
            let elastic = 0.5 * stiffness * (old.closure_m + next_closure);
            let expected_power = 1.0e-20 * elastic * velocity * velocity;

            assert!(receipt.dissipated_j > 0.0);
            assert!(-outputs[2] > 0.0);
            assert!((-outputs[2] - expected_power).abs() <= 1.0e-12 * expected_power);
            assert!(
                (receipt.dissipated_j - dt * expected_power).abs() <= 1.0e-12 * dt * expected_power
            );
        }

        // Opening at 3 m/s with chi=1 s/m saturates the normal force at
        // zero. The old potential becomes heat rather than attraction.
        let unloading = lossy_contact(stiffness, 1.0, 1.0);
        let inputs = [-3.0, 0.0, 300.0];
        let case = Case {
            state: &old,
            stale: &stale,
            dt_s: dt,
            inputs: &inputs,
        };
        let (receipt, outputs) = stepped(check_component(&unloading, &case, cx).unwrap());
        let initial_potential = 0.5 * stiffness * old.closure_m.powi(2);

        assert_eq!(outputs[0], 0.0);
        assert_eq!(outputs[1], 0.0);
        assert!((-dt * outputs[2] - initial_potential).abs() <= 1.0e-12);
        assert!((receipt.stored_change_j + initial_potential).abs() <= 1.0e-12);
        assert!((receipt.dissipated_j - initial_potential).abs() <= 1.0e-12);
    });
}

#[test]
fn lossy_contact_end_inputs_follow_the_heat_input_without_changing_the_step() {
    let plain = lossy_contact(HERTZ_K, 1.5, 0.25);
    let told = plain.clone().reading_ends([true, true]).unwrap();
    let names: Vec<_> = told.ports().iter().map(|p| p.id().as_str()).collect();
    assert_eq!(names, ["a", "b", "heat", "a.end", "b.end"]);
    assert_eq!(told.line_ends(), [Some(3), Some(4)]);

    with_cx(|cx| {
        let old = ContactState::at(1.0e-4);
        let line = [5.0e-4, 0.0, 0.0, 2.0e-4];
        let told_line = [5.0e-4, 0.0, 0.0, 2.0e-4, 1.0e-3, 0.0, 0.0, 4.0e-4];
        let (mut next_a, mut next_b) = (old, old);
        let (mut out_a, mut out_b) = ([0.0; 3], [0.0; 3]);
        let receipt_a = plain
            .propose_on_line(
                &old,
                1.0e-6,
                &[0.5, -0.2, 300.0],
                &line,
                &mut next_a,
                &mut out_a,
                cx,
            )
            .unwrap();
        let receipt_b = told
            .propose_on_line(
                &old,
                1.0e-6,
                &[0.5, -0.2, 300.0, 0.7, -0.1],
                &told_line,
                &mut next_b,
                &mut out_b,
                cx,
            )
            .unwrap();
        assert_eq!(receipt_a, receipt_b);
        assert_eq!((next_a, out_a), (next_b, out_b));
    });
}

#[test]
fn zero_loss_keeps_the_elastic_ports_and_bitwise_proposal() {
    let plain = contact(HERTZ_K, 1.5, 0.0);
    let explicit_zero = lossy_contact(HERTZ_K, 1.5, 0.0);
    assert_eq!(explicit_zero.ports().len(), 2);
    with_cx(|cx| {
        let old = ContactState::at(1.0e-4);
        let stale = ContactState::at(-1.0);
        let (mut next_a, mut next_b) = (stale, stale);
        let (mut out_a, mut out_b) = ([0.0; 2], [0.0; 2]);
        let a = plain
            .propose(&old, 1.0e-6, &[0.5, -0.2], &mut next_a, &mut out_a, cx)
            .unwrap();
        let b = explicit_zero
            .propose(&old, 1.0e-6, &[0.5, -0.2], &mut next_b, &mut out_b, cx)
            .unwrap();
        let receipt_bits = |proposal| {
            let Proposal::Step(receipt) = proposal else {
                panic!("the elastic contact must take its step");
            };
            [
                receipt.stored_change_j,
                receipt.dissipated_j,
                receipt.numerical_loss_j,
                receipt.entropy_generated_j_per_k,
                receipt.energy_tolerance_j,
                receipt.residual,
                receipt.tolerance,
            ]
            .map(f64::to_bits)
        };
        assert_eq!(receipt_bits(a), receipt_bits(b));
        assert_eq!(out_a.map(f64::to_bits), out_b.map(f64::to_bits));
        assert_eq!(
            [next_a.closure_m, next_a.reaction_n].map(f64::to_bits),
            [next_b.closure_m, next_b.reaction_n].map(f64::to_bits)
        );
    });
}

#[test]
fn invalid_heat_temperature_refuses_before_publishing_candidate_or_outputs() {
    let interface = lossy_contact(HERTZ_K, 1.5, 0.25);
    let old = ContactState::at(1.0e-4);
    let stale = ContactState::at(-1.0);
    with_cx(|cx| {
        // The positive subnormal also refuses: finite physical heat divided
        // by it cannot be represented as a finite entropy increment.
        for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            let mut next = stale;
            let mut outputs = [11.0, 22.0, 33.0];
            assert!(
                interface
                    .propose(
                        &old,
                        1.0e-6,
                        &[0.5, -0.2, temperature],
                        &mut next,
                        &mut outputs,
                        cx,
                    )
                    .is_err()
            );
            assert_eq!(next, stale);
            assert_eq!(outputs, [11.0, 22.0, 33.0]);
        }
    });
}

const MASS_A: f64 = 0.01;
const MASS_B: f64 = 0.02;
const VELOCITY_A: f64 = -1.5;
const VELOCITY_B: f64 = 0.5;
const CAPACITY: f64 = 10.0;
const TEMPERATURE: f64 = 300.0;
const IMPACT_DT: f64 = 2.0e-7;

struct Impact {
    system: System,
    a: NodeId,
    b: NodeId,
    contact: NodeId,
    thermal: Option<NodeId>,
}

fn impact(chi: f64) -> Impact {
    let mut builder = SystemBuilder::new();
    let a = builder
        .add_with_tangent(
            "body-a",
            Mass::new(MASS_A).unwrap(),
            Motion {
                position_m: APPROACH_M,
                velocity_m_s: VELOCITY_A,
            },
        )
        .unwrap();
    let b = builder
        .add_with_tangent(
            "body-b",
            Mass::new(MASS_B).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: VELOCITY_B,
            },
        )
        .unwrap();
    let contact = builder
        .add_on_line(
            "contact",
            lossy_contact(HERTZ_K, 1.5, chi),
            ContactState::at(-APPROACH_M),
        )
        .unwrap();

    for (end, body) in [("a", a), ("b", b)] {
        builder
            .bond(
                builder.port(contact, end).unwrap(),
                builder.port(body, "p").unwrap(),
            )
            .unwrap();
    }
    let thermal = if chi > 0.0 {
        let body = builder
            .add_with_tangent(
                "thermal-body",
                HeatCapacity::new(CAPACITY).unwrap(),
                TEMPERATURE,
            )
            .unwrap();
        builder
            .bond(
                builder.port(contact, "heat").unwrap(),
                builder.port(body, "heat").unwrap(),
            )
            .unwrap();
        Some(body)
    } else {
        None
    };
    Impact {
        system: builder.seal().unwrap(),
        a,
        b,
        contact,
        thermal,
    }
}

fn heat_j(device: &Impact) -> f64 {
    device.thermal.map_or(0.0, |node| {
        CAPACITY * (*device.system.state::<f64>(node).unwrap() - TEMPERATURE)
    })
}

/// Exact isolated-impact relation, evaluated independently with libm:
/// x - ln(1+x) = -e*x - ln(1-e*x), x = chi * approach speed.
fn restitution(x: f64) -> f64 {
    if x == 0.0 {
        return 1.0;
    }
    let target = x - x.ln_1p();
    let (mut lo, mut hi) = (0.0, (1.0 / x).min(1.0));
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let y = mid * x;
        let value = -y - (-y).ln_1p();
        if value < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn two_body_impact_conserves_momentum_and_converts_the_lost_motion_to_heat() {
    let initial_energy = 0.5 * (MASS_A * VELOCITY_A.powi(2) + MASS_B * VELOCITY_B.powi(2));
    let initial_momentum = MASS_A * VELOCITY_A + MASS_B * VELOCITY_B;
    let approach_speed = VELOCITY_B - VELOCITY_A;

    with_cx(|cx| {
        for chi in [0.0, 0.25] {
            let mut device = impact(chi);
            let (mut entered, mut exited) = (false, false);
            for _ in 0..12_000 {
                device.system.step(IMPACT_DT, cx).unwrap();
                let a = device.system.state::<Motion>(device.a).unwrap();
                let b = device.system.state::<Motion>(device.b).unwrap();
                let state = device.system.state::<ContactState>(device.contact).unwrap();
                let kinetic =
                    0.5 * (MASS_A * a.velocity_m_s.powi(2) + MASS_B * b.velocity_m_s.powi(2));
                let potential = HERTZ_K * state.closure_m.max(0.0).powf(2.5) / 2.5;

                assert!(
                    (MASS_A * a.velocity_m_s + MASS_B * b.velocity_m_s - initial_momentum).abs()
                        < 2.0e-12
                );
                assert!((state.closure_m - (b.position_m - a.position_m)).abs() < 1.0e-12);
                // The absolute-temperature state rounds each small heat
                // increment; this checks that physical inventory directly,
                // independently of the compensated receipt accounts.
                assert!((kinetic + potential + heat_j(&device) - initial_energy).abs() < 4.0e-9);
                assert!(state.reaction_n >= 0.0);

                entered |= state.closure_m > 0.0;
                if entered && state.closure_m <= 0.0 {
                    exited = true;
                    break;
                }
            }
            assert!(
                entered && exited,
                "the bounded run must complete the collision"
            );

            let a = device.system.state::<Motion>(device.a).unwrap();
            let b = device.system.state::<Motion>(device.b).unwrap();
            let observed = (a.velocity_m_s - b.velocity_m_s) / approach_speed;
            let expected = restitution(chi * approach_speed);
            assert!(
                (observed - expected).abs() < 2.0e-4,
                "chi={chi}: restitution {observed}, expected {expected}"
            );

            let ledger = device.system.ledger();
            assert!(ledger.energy_residual_j().abs() < 1.0e-11);
            assert!(ledger.numerical_loss_j().abs() < 1.0e-9);
            assert_eq!(ledger.boundary_energy_in_j(), 0.0);
            assert_eq!(ledger.sink_heat_j(), 0.0);
            if chi > 0.0 {
                assert!(heat_j(&device) > 0.0);
                assert!(ledger.entropy_generated_j_per_k() > 0.0);
                assert!(observed < 1.0);
            }
        }
    });
}

fn observed(device: &Impact) -> (Vec<u64>, fs_system::Ledger) {
    let a = device.system.state::<Motion>(device.a).unwrap();
    let b = device.system.state::<Motion>(device.b).unwrap();
    let contact = device.system.state::<ContactState>(device.contact).unwrap();
    let temperature = device
        .thermal
        .map_or(0.0, |node| *device.system.state::<f64>(node).unwrap());
    (
        vec![
            a.position_m.to_bits(),
            a.velocity_m_s.to_bits(),
            b.position_m.to_bits(),
            b.velocity_m_s.to_bits(),
            contact.closure_m.to_bits(),
            contact.reaction_n.to_bits(),
            temperature.to_bits(),
            device.system.time_s().to_bits(),
            device.system.steps(),
        ],
        device.system.ledger().clone(),
    )
}

#[test]
fn heated_impact_checkpoint_replay_and_cancelled_retry_do_not_duplicate_heat() {
    let mut device = impact(0.25);
    let expected = with_cx(|cx| {
        // Enter contact so the captured state has stored potential, ongoing
        // loss and an already heated destination.
        for _ in 0..150 {
            device.system.step(IMPACT_DT, cx).unwrap();
        }
        assert!(
            device
                .system
                .state::<ContactState>(device.contact)
                .unwrap()
                .closure_m
                > 0.0
        );
        assert!(heat_j(&device) > 0.0);
        let checkpoint = device.system.checkpoint().unwrap();
        let before = observed(&device);

        for _ in 0..17 {
            device.system.step(IMPACT_DT, cx).unwrap();
        }
        let expected = observed(&device);
        device.system.restore(&checkpoint).unwrap();
        assert_eq!(observed(&device), before);

        cx.cancel_gate().request();
        assert!(device.system.step(IMPACT_DT, cx).is_err());
        assert_eq!(observed(&device), before);
        expected
    });

    // A fresh cancellation scope retries against the same accepted banks
    // and disposable candidate storage.
    with_cx(|cx| {
        for _ in 0..17 {
            device.system.step(IMPACT_DT, cx).unwrap();
        }
    });
    assert_eq!(observed(&device), expected);
}
