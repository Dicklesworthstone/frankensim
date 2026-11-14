//! G0/G2: dual power elements admit parallel networks. References are the
//! closed-form midpoint map and an independently reduced two-state circuit.

#![allow(clippy::float_cmp)]

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::PortKind;
use fs_ports::component::conformance::{
    Case, Failure, Observed, Reversal, ReversalAllowance, check_component, check_port_tangent,
    check_time_reversal,
};
use fs_ports::component::{Component, PortTangent};
use fs_system::elements::{
    Admittance, ElementError, FlowSource, HeatCapacity, Mass, Motion, SourceAccount, Storage,
    Waveform,
};
use fs_system::{NodeId, SystemBuilder};

const ELECTRICAL: PortKind = PortKind::ElectricalVoltageCurrent;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey { seed: 19, kernel_id: 31, tile: 0, iteration: 0 };
        f(&Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} vs {expected}; tolerance {tolerance}");
}

fn conforms<C>(element: &C, state: &C::State, stale: &C::State, dt_s: f64, inputs: &[f64]) -> Vec<f64>
where
    C: PortTangent,
    C::State: PartialEq + core::fmt::Debug,
{
    with_cx(|cx| {
        let case = Case { state, stale, dt_s, inputs };
        let observed = check_component(element, &case, cx).expect("component conformance");
        let Observed::Stepped { outputs, .. } = observed else {
            panic!("linear element shortened the step");
        };
        let largest = |values: &[f64]| values.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        let h = largest(inputs) * f64::EPSILON.cbrt();
        let absolute = 1.0e-12_f64.max(8.0 * f64::EPSILON * largest(&outputs) / h);
        check_port_tangent(element, &case, 1.0e-6, absolute, cx).expect("port tangent");
        outputs
    })
}

/// Unit roundoff of f64.
const U: f64 = f64::EPSILON / 2.0;

/// `n u / (1 - n u)`: a product of `n` factors `(1 + d)`, `|d| <= u`, lies
/// within it of 1.
fn gamma(n: u32) -> f64 {
    let nu = f64::from(n) * U;
    nu / (1.0 - nu)
}

/// A storage element there and back under `reverse`. Its state (displacement,
/// twist, charge or volume) is even, so the physical `R` is the identity; the
/// flow is odd and the effort even. The increment `t = fl(dt f)` flips exactly
/// with the flow, so the backward step ends at `q2 = fl(fl(q0 + t) - t)`.
/// Counted from `Storage::propose` and `stored_energy_j` (halving is exact):
/// - State: within `gamma(2) (|q0| + |q1|)` of `q0`.
/// - Effort `fl(fl(q0 / 2 + q1 / 2) / C)`, two roundings each way: within
///   the state's gap over `2 C`, plus `gamma(2) (|q0| + 2 |q1| + |q2|) / 2 C`.
/// - Energy `fl(q / 2 fl(q / C))`, two roundings each: within the state's gap
///   times `(|q0| + |q2|) / 2 C`, plus `gamma(2) (q0^2 + q2^2) / 2 C`.
fn storage_reverses(
    storage: &Storage,
    capacity: f64,
    state: f64,
    dt_s: f64,
    flow: f64,
    reverse: &dyn Fn(&f64) -> f64,
) -> Result<Reversal, Failure<ElementError>> {
    let observe = |state: &f64| vec![*state];
    let allowance = |start: &f64, forward_end: &f64, back_end: &f64| {
        let [q0, q1, q2] = [start, forward_end, back_end].map(|state| state.abs());
        let travel = gamma(2) * (q0 + q1);
        ReversalAllowance {
            outputs: vec![(travel + gamma(2) * (q0 + 2.0 * q1 + q2)) / (2.0 * capacity)],
            observed: vec![travel],
            energy_j: (travel * (q0 + q2) + gamma(2) * (q0 * q0 + q2 * q2)) / (2.0 * capacity),
        }
    };
    let inputs = [flow];
    let case = Case {
        state: &state,
        stale: &9.0,
        dt_s,
        inputs: &inputs,
    };
    with_cx(|cx| check_time_reversal(storage, &case, reverse, &observe, &allowance, cx))
}

#[test]
fn dual_elements_conform_in_each_power_domain_and_match_their_laws() {
    // The largest gap / allowance of the storage reversals, per law, and their count.
    let (mut worst, mut reversals) = ([0.0_f64; 3], 0);
    for kind in [
        PortKind::MechanicalForceVelocity,
        PortKind::RotationalTorqueAngularVelocity,
        ELECTRICAL,
        PortKind::FluidPressureFlux,
    ] {
        let storage = Storage::in_domain(kind, 0.25).unwrap();
        let branch = Admittance::in_domain(kind, 0.4).unwrap();
        let source = FlowSource::in_domain(kind, Waveform::sine(0.7, 3.0, 0.2, -1.0).unwrap());
        for i in 0..100 {
            let x = f64::from(i);
            let dt = 0.0001 + x * 0.0002;
            let state = 0.3 * (0.7 * x).sin();
            let flow = 2.0 * (0.3 * x).cos();
            let out = conforms(&storage, &state, &9.0, dt, &[flow]);
            close(out[0], (state + 0.5 * dt * flow) / 0.25, 1.0e-14);
            let reversal = storage_reverses(&storage, 0.25, state, dt, flow, &|q: &f64| *q)
                .unwrap_or_else(|failure| {
                    panic!("{kind:?} storage reversal at {state}, {dt} s, {flow}: {failure}")
                });
            // A gap of zero is no part of its allowance, also of an allowance of zero.
            let ratio = |gap: f64, bound: f64| if gap > 0.0 { gap / bound } else { 0.0 };
            let allowance = &reversal.allowance;
            worst[0] = worst[0].max(ratio(reversal.output_gaps[0], allowance.outputs[0]));
            worst[1] = worst[1].max(ratio(reversal.observed_gaps[0], allowance.observed[0]));
            worst[2] = worst[2].max(ratio(reversal.energy_gap_j, allowance.energy_j));
            reversals += 1;

            let (a, b, temperature) = (7.0 * (0.4 * x).sin(), 4.0 * (0.6 * x).cos(), 250.0 + x);
            let out = conforms(&branch, &(), &(), dt, &[a, b, temperature]);
            close(out[0], 0.4 * (a - b), 1.0e-14);
            assert_eq!(out[0], -out[1]);
            close(out[2], -0.4 * (a - b).powi(2), 1.0e-13);
            // A common voltage/pressure/force offset changes no branch flow.
            let shifted = conforms(&branch, &(), &(), dt, &[a + 8.0, b + 8.0, temperature]);
            close(shifted[0], out[0], 2.0e-15);
            close(shifted[2], out[2], 1.0e-13);

            let mut account = SourceAccount::starting_at(0.01 * x);
            account.energy_in_j = -0.3 * x;
            let out = conforms(&source, &account, &SourceAccount::default(), dt, &[a]);
            // Integrate the prescribed sine independently of Waveform::mean.
            let omega = 6.0 * PI;
            let phase = omega * account.time_s() + 0.2;
            let expected = -1.0 + 0.7 * (phase.cos() - (phase + omega * dt).cos()) / (omega * dt);
            close(out[0], expected, 2.0e-13);
        }
    }
    // The storage element reversed in time at every point, in every domain.
    println!(
        "storage: {reversals} reversals; largest gap / allowance: outputs {:.3e}, observed {:.3e}, energy {:.3e}",
        worst[0], worst[1], worst[2]
    );
    assert_eq!(reversals, 400);
    assert_eq!(Storage::new(0.25).unwrap().ports()[0].carrier(), Storage::in_domain(PortKind::MechanicalForceVelocity, 0.25).unwrap().ports()[0].carrier());
    assert_eq!(Admittance::new(0.4).unwrap().ports()[0].carrier(), Admittance::in_domain(PortKind::MechanicalForceVelocity, 0.4).unwrap().ports()[0].carrier());
    assert_eq!(FlowSource::new(Waveform::constant(1.0).unwrap()).ports()[0].carrier(), FlowSource::in_domain(PortKind::MechanicalForceVelocity, Waveform::constant(1.0).unwrap()).ports()[0].carrier());
}

#[test]
fn a_storage_reversal_that_negates_the_state_is_refused_at_the_outputs() {
    // Taking the state for odd sends the backward step from -q1 under -f, so
    // it ends near -(q0 + 2 dt f) and its effort has the wrong sign.
    let storage = Storage::new(0.25).unwrap();
    let failure = storage_reverses(&storage, 0.25, 0.2, 0.01, 1.0, &|q: &f64| -*q)
        .expect_err("the state is even");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    println!("refused: {law} ({detail})");
    assert_eq!(
        law,
        "the backward outputs are the forward outputs times their parity"
    );
    assert!(detail.starts_with("output 0:"), "{detail}");
}

#[test]
fn dual_elements_refuse_invalid_coefficients_inputs_and_temperature() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Storage::new(value).is_err());
        assert!(Admittance::new(value).is_err());
    }
    with_cx(|cx| {
        let storage = Storage::new(0.1).unwrap();
        let branch = Admittance::new(0.5).unwrap();
        let source = FlowSource::new(Waveform::constant(-2.0).unwrap());
        assert!(matches!(storage.propose(&1.0, 0.01, &[f64::NAN], &mut 0.0, &mut [0.0], cx), Err(ElementError::Input { .. })));
        assert!(matches!(source.propose(&SourceAccount::default(), 0.01, &[f64::INFINITY], &mut SourceAccount::default(), &mut [0.0], cx), Err(ElementError::Input { .. })));
        for temperature in [0.0, -50.0] {
            assert_eq!(branch.propose(&(), 0.01, &[1.0, 0.0, temperature], &mut (), &mut [0.0; 3], cx), Err(ElementError::Temperature { kelvin: temperature }));
        }
        assert!(matches!(branch.propose(&(), 0.01, &[1.0, f64::NAN, 300.0], &mut (), &mut [0.0; 3], cx), Err(ElementError::Input { .. })));

        cx.cancel_gate().request();
        assert_eq!(storage.propose(&1.0, 0.01, &[1.0], &mut 0.0, &mut [0.0], cx), Err(ElementError::Cancelled));
        assert_eq!(branch.propose(&(), 0.01, &[1.0, 0.0, 300.0], &mut (), &mut [0.0; 3], cx), Err(ElementError::Cancelled));
        assert_eq!(source.propose(&SourceAccount::default(), 0.01, &[1.0], &mut SourceAccount::default(), &mut [0.0], cx), Err(ElementError::Cancelled));
    });
}

fn port(builder: &SystemBuilder, node: NodeId, name: &str) -> fs_system::PortRef {
    builder.port(node, name).unwrap()
}

#[test]
fn parallel_rc_branches_preserve_charge_and_deliver_all_lost_energy_as_heat() {
    let (c1, c2, g1, g2) = (0.04_f64, 0.12_f64, 0.1_f64, 0.2_f64);
    let (v1, v2, heat_capacity, initial_temperature) = (9.0_f64, -2.0_f64, 10.0_f64, 300.0_f64);
    let mut builder = SystemBuilder::new();
    let left = builder.add_with_tangent("left", Storage::in_domain(ELECTRICAL, c1).unwrap(), c1 * v1).unwrap();
    let right = builder.add_with_tangent("right", Storage::in_domain(ELECTRICAL, c2).unwrap(), c2 * v2).unwrap();
    let branch1 = builder.add_with_tangent("branch1", Admittance::in_domain(ELECTRICAL, g1).unwrap(), ()).unwrap();
    let branch2 = builder.add_with_tangent("branch2", Admittance::in_domain(ELECTRICAL, g2).unwrap(), ()).unwrap();
    let body = builder.add_with_tangent("body", HeatCapacity::new(heat_capacity).unwrap(), initial_temperature).unwrap();
    builder.common_effort(port(&builder, left, "p"), &[port(&builder, branch1, "a"), port(&builder, branch2, "a")]).unwrap();
    builder.common_effort(port(&builder, right, "p"), &[port(&builder, branch1, "b"), port(&builder, branch2, "b")]).unwrap();
    let heats = [port(&builder, branch1, "heat"), port(&builder, branch2, "heat")];
    builder.common_effort(port(&builder, body, "heat"), &heats).unwrap();
    let mut system = builder.seal().unwrap();
    for heat in heats {
        assert!(system.guess(heat, &[initial_temperature]));
    }
    let (dt, steps) = (0.001_f64, 600_i32);
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });

    let q1 = *system.state::<f64>(left).unwrap();
    let q2 = *system.state::<f64>(right).unwrap();
    let total_charge = c1 * v1 + c2 * v2;
    close(q1 + q2, total_charge, 1.0e-13);
    let half = 0.5 * dt * (g1 + g2) * (1.0 / c1 + 1.0 / c2);
    let difference = (v1 - v2) * ((1.0 - half) / (1.0 + half)).powi(steps);
    let average = total_charge / (c1 + c2);
    close(q1 / c1, average + c2 / (c1 + c2) * difference, 1.0e-11);
    close(q2 / c2, average - c1 / (c1 + c2) * difference, 1.0e-11);

    let energy_before = 0.5 * (c1 * v1 * v1 + c2 * v2 * v2);
    let energy_after = 0.5 * (q1 * q1 / c1 + q2 * q2 / c2);
    let lost = energy_before - energy_after;
    let ledger = system.ledger();
    close(ledger.nodes[branch1.index()].dissipated_j.value(), lost * g1 / (g1 + g2), 1.0e-11);
    close(ledger.nodes[branch2.index()].dissipated_j.value(), lost * g2 / (g1 + g2), 1.0e-11);
    close(ledger.nodes[body.index()].stored_change_j.value(), lost, 1.0e-11);
    let temperature = *system.state::<f64>(body).unwrap();
    close(heat_capacity * (temperature - initial_temperature), lost, 2.0e-10);
    close(ledger.energy_residual_j(), 0.0, 1.0e-11);
    assert!(ledger.entropy_generated_j_per_k() > 0.0);
}

#[test]
fn parallel_lc_network_has_the_midpoint_frequency_and_conserves_energy() {
    let (c, l, v0, i0) = (0.04_f64, 0.3_f64, 7.0_f64, -1.0_f64);
    let mut builder = SystemBuilder::new();
    let capacitor = builder.add_with_tangent("capacitor", Storage::in_domain(ELECTRICAL, c).unwrap(), c * v0).unwrap();
    let inductor = builder.add_with_tangent("inductor", Mass::in_domain(ELECTRICAL, l).unwrap(), Motion { position_m: 0.0, velocity_m_s: i0 }).unwrap();
    builder.common_effort(port(&builder, capacitor, "p"), &[port(&builder, inductor, "p")]).unwrap();
    let mut system = builder.seal().unwrap();
    let (dt, steps) = (0.002_f64, 4_000_i32);
    with_cx(|cx| {
        for _ in 0..steps {
            system.step(dt, cx).unwrap();
        }
    });
    let omega = 1.0 / (l * c).sqrt();
    let angle = f64::from(steps) * 2.0 * (0.5 * omega * dt).atan();
    let expected_v = v0 * angle.cos() - i0 / (c * omega) * angle.sin();
    let expected_i = i0 * angle.cos() + c * omega * v0 * angle.sin();
    let charge = *system.state::<f64>(capacitor).unwrap();
    let motion = *system.state::<Motion>(inductor).unwrap();
    close(charge / c, expected_v, 1.0e-10);
    close(motion.velocity_m_s, expected_i, 1.0e-10);
    close(charge + motion.position_m, c * v0, 1.0e-13);
    let energy_before = 0.5 * (c * v0 * v0 + l * i0 * i0);
    let energy_after = 0.5 * (charge * charge / c + l * motion.velocity_m_s.powi(2));
    close(energy_after, energy_before, 1.0e-11);
    close(system.ledger().energy_residual_j(), 0.0, 1.0e-11);
    assert_eq!(system.ledger().entropy_generated_j_per_k(), 0.0);
}

#[test]
fn a_current_driven_parallel_rlc_matches_the_monolithic_equations_and_source_work() {
    let (c, l, g, v0, i0) = (0.01_f64, 0.2_f64, 0.05_f64, 2.0_f64, -0.1_f64);
    let mut builder = SystemBuilder::new();
    let capacitor = builder.add_with_tangent("capacitor", Storage::in_domain(ELECTRICAL, c).unwrap(), c * v0).unwrap();
    let inductor = builder.add_with_tangent("inductor", Mass::in_domain(ELECTRICAL, l).unwrap(), Motion { position_m: 0.0, velocity_m_s: i0 }).unwrap();
    let resistor = builder.add_with_tangent("resistor", Admittance::in_domain(ELECTRICAL, g).unwrap(), ()).unwrap();
    let source = builder.add_with_tangent("source", FlowSource::in_domain(ELECTRICAL, Waveform::sine(-0.2, 3.0, 0.3, -0.4).unwrap()), SourceAccount::default()).unwrap();
    let body = builder.add_with_tangent("body", HeatCapacity::new(10.0).unwrap(), 300.0).unwrap();
    builder.common_effort(port(&builder, capacitor, "p"), &[port(&builder, inductor, "p"), port(&builder, resistor, "a"), port(&builder, source, "p")]).unwrap();
    builder.boundary(port(&builder, resistor, "b"), Box::new(|_, _, values| values[0] = 0.0)).unwrap();
    let heat = port(&builder, resistor, "heat");
    builder.common_effort(port(&builder, body, "heat"), &[heat]).unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));

    let (dt, steps) = (0.001_f64, 2_000_u32);
    let (mut voltage, mut current, mut heat_j, mut work_j) = (v0, i0, 0.0_f64, 0.0_f64);
    with_cx(|cx| {
        for step in 0..steps {
            // C V' = J - G V - I; L I' = V. Eliminate midpoint I
            // from these two equations, independently of the port solver.
            let half = 0.5 * dt;
            let phase = 6.0 * PI * f64::from(step) * dt + 0.3;
            let supplied = 0.4 + 0.2 * (phase.cos() - (phase + 6.0 * PI * dt).cos()) / (6.0 * PI * dt);
            let mean_v = (voltage + half / c * (supplied - current))
                / (1.0 + half * g / c + half * half / (c * l));
            voltage = 2.0 * mean_v - voltage;
            current += dt / l * mean_v;
            heat_j += dt * g * mean_v * mean_v;
            work_j += dt * supplied * mean_v;
            system.step(dt, cx).unwrap();
            close(*system.state::<f64>(capacitor).unwrap() / c, voltage, 2.0e-10);
            close(system.state::<Motion>(inductor).unwrap().velocity_m_s, current, 2.0e-11);
        }
    });
    let account = system.state::<SourceAccount>(source).unwrap();
    assert_eq!(account.time_s(), system.time_s());
    close(-account.energy_in_j, work_j, 2.0e-11);
    let ledger = system.ledger();
    close(ledger.nodes[resistor.index()].dissipated_j.value(), heat_j, 2.0e-11);
    close(ledger.nodes[body.index()].stored_change_j.value(), heat_j, 2.0e-11);
    let stored_delta = 0.5 * (c * (voltage * voltage - v0 * v0) + l * (current * current - i0 * i0));
    close(work_j, stored_delta + heat_j, 2.0e-11);
    close(ledger.energy_residual_j(), 0.0, 2.0e-11);
}
