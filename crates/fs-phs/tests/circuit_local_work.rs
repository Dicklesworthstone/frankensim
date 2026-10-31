//! Independent branch work for nonlinear capacitors and reciprocal magnetic storage.
//! Global energy closure alone cannot identify each original component's voltage.
use fs_phs::circuit::capacitance::{
    CapacitorChargeCurve, CapacitorCurve, assemble_circuit_with_capacitor_curves,
};
use fs_phs::circuit::{Branch, CircuitDae, CircuitGraph, MutualInductance};
use fs_phs::{PortHamiltonian, StepWorkspace, discrete_gradient, step};
use std::cell::Cell;

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}
fn law_a(scale: f64) -> CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(
        [
            (-4.0, -7.0),
            (-1.0, -1.0),
            (0.0, 0.0),
            (1.0, 1.0),
            (4.0, 7.0),
        ]
        .into_iter()
        .map(|(q, v)| (q * scale, v))
        .collect(),
    )
    .unwrap()
}
fn law_b(scale: f64) -> CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(
        [
            (-4.0, -5.0),
            (-0.5, -1.5),
            (0.0, 0.0),
            (0.5, 1.5),
            (4.0, 5.0),
        ]
        .into_iter()
        .map(|(q, v)| (q * scale, v))
        .collect(),
    )
    .unwrap()
}
// Independent antiderivatives of the two explicit V(q) fixtures.
fn energy_a(q: f64) -> f64 {
    let q = q.abs();
    if q <= 1.0 {
        0.5 * q * q
    } else {
        q * q - q + 0.5
    }
}
fn energy_b(q: f64) -> f64 {
    let q = q.abs();
    if q <= 0.5 {
        1.5 * q * q
    } else {
        0.5 * q * q + q - 0.25
    }
}
fn two_caps(a: f64, b: f64) -> CircuitDae {
    assemble_circuit_with_capacitor_curves(
        &CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: a }),
                (1, 0, Branch::CurrentSource { port: 0 }),
                (2, 0, Branch::Capacitor { farads: b }),
                (2, 0, Branch::CurrentSource { port: 1 }),
                (3, 0, Branch::Capacitor { farads: 1.0 }),
                (3, 0, Branch::VoltageSource { port: 2 }),
            ],
            transformers: vec![],
        },
        &[
            CapacitorCurve {
                branch: 0,
                curve: law_a(a),
            },
            CapacitorCurve {
                branch: 2,
                curve: law_b(b),
            },
        ],
        &[],
    )
    .unwrap()
}

#[test]
fn two_changing_capacitors_keep_independent_source_work_beside_a_huge_held_bias() {
    for (a, b) in [(1.0, 1.0), (1e-9, 1e-12)] {
        let dae = two_caps(a, b);
        let input = [-1.5 * a, 1.25 * b, 1e8];
        let mut before = vec![0.0; dae.system.state_dim()];
        before[dae.charge_index[0]] = 2.0 * a;
        before[dae.charge_index[1]] = 0.25 * b;
        before[dae.charge_index[2]] = 1e8;
        let before = dae.consistent_initial_state(&before, &input).unwrap();
        let (step, _) = dae.step_audited(&before, &input, 1.0).unwrap();
        close(step.x[dae.charge_index[0]], 0.5 * a, a);
        close(step.x[dae.charge_index[1]], 1.5 * b, b);
        assert_eq!(
            step.x[dae.charge_index[2]].to_bits(),
            before[dae.charge_index[2]].to_bits()
        );
        // Each held current delivers only its own capacitor's integral work.
        // Simultaneous knot crossings have different mean-voltage corrections.
        let delta_a = a * (energy_a(0.5) - energy_a(2.0));
        let delta_b = b * (energy_b(1.5) - energy_b(0.25));
        close(step.y[0], 19.0 / 12.0, 1.0);
        close(step.y[1], 73.0 / 40.0, 1.0);
        close(input[0] * step.y[0], delta_a, a);
        close(input[1] * step.y[1], delta_b, b);
        close(step.delta_h, delta_a + delta_b, a.max(b));
        close(step.supplied, delta_a + delta_b, a.max(b));
        close(step.dissipated, 0.0, 1.0);
        close(step.y[2], 0.0, 1.0);
        if a < 1.0 {
            assert_eq!(
                dae.system.hamiltonian(&before),
                dae.system.hamiltonian(&step.x)
            );
        }
    }
}

#[test]
fn nonlinear_lc_step_preserves_the_capacitor_work_and_original_midpoint_current() {
    let dae = assemble_circuit_with_capacitor_curves(
        &CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (1, 0, Branch::Inductor { henries: 2.0 }),
            ],
            transformers: vec![],
        },
        &[CapacitorCurve {
            branch: 0,
            curve: law_a(1.0),
        }],
        &[],
    )
    .unwrap();
    let mut before = vec![0.0; dae.system.state_dim()];
    before[dae.charge_index[0]] = 2.0;
    let before = dae.consistent_initial_state(&before, &[]).unwrap();
    let dt = 2.0;
    let (step, _) = dae.step_audited(&before, &[], dt).unwrap();
    let q = step.x[dae.charge_index[0]];
    let flux = step.x[dae.flux_index[0]];
    let voltage = f64::midpoint(
        before[dae.node_potential_index[0]],
        step.x[dae.node_potential_index[0]],
    );
    // dq=-dt*(phi0+phi1)/(2L), dphi=dt*Vbar and exact capacitor work
    // imply 3q1²-8q1+3=0 across the knot, with the physical root below one.
    let expected_q = (4.0 - 7.0_f64.sqrt()) / 3.0;
    close(q, expected_q, 1.0);
    close(flux, 2.0 * (2.0 - expected_q), 1.0);
    close((q - 2.0) * voltage, energy_a(q) - energy_a(2.0), 1.0);
    close(q - 2.0, -dt * flux / 4.0, 1.0);
    close(flux, dt * voltage, 1.0);
    close(energy_a(q) + flux * flux / 4.0, 2.5, 1.0);
    close(step.balance_residual(), 0.0, 1.0);
}

#[test]
fn reciprocal_inductor_efforts_do_not_receive_another_capacitors_knot_correction() {
    let dae = assemble_circuit_with_capacitor_curves(
        &CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Inductor { henries: 2.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
                (2, 0, Branch::Inductor { henries: 3.0 }),
                (2, 0, Branch::Resistor { ohms: 1.0 }),
            ],
            transformers: vec![],
        },
        &[CapacitorCurve {
            branch: 1,
            curve: law_a(1.0),
        }],
        &[MutualInductance {
            branches: (0, 2),
            henries: 1.0,
        }],
    )
    .unwrap();
    let mut before = vec![0.0; dae.system.state_dim()];
    let mut after = before.clone();
    for (slot, a, b) in [
        (dae.flux_index[0], 1.0, 2.5),
        (dae.flux_index[1], -0.5, 1.5),
        (dae.charge_index[0], 2.0, 0.5),
    ] {
        before[slot] = a;
        after[slot] = b;
    }
    let effort = discrete_gradient(dae.system.storage(), &before, &after);
    // Inverse of [[2,1],[1,3]] is [[3,-1],[-1,2]]/5. It acts on the
    // true flux midpoint, independently of the simultaneously changing curve.
    close(effort[dae.flux_index[0]], 0.95, 1.0);
    close(effort[dae.flux_index[1]], -0.15, 1.0);
    close(effort[dae.charge_index[0]], 19.0 / 12.0, 1.0);
    let magnetic_work = 1.5 * effort[dae.flux_index[0]] + 2.0 * effort[dae.flux_index[1]];
    close(magnetic_work, 1.125, 1.0);
    close(-1.5 * effort[dae.charge_index[0]], -2.375, 1.0);
}

#[test]
fn zero_and_one_ulp_charge_increments_keep_their_own_constitutive_efforts() {
    let dae = two_caps(1.0, 1.0);
    let next = f64::from_bits(1.0_f64.to_bits() + 1);
    for (a, b, expected) in [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.5),
        (1.0, next, next),
        (next, 1.0, next),
        (-1.0, -next, -next),
    ] {
        let mut before = vec![0.0; dae.system.state_dim()];
        let mut after = before.clone();
        before[dae.charge_index[0]] = a;
        after[dae.charge_index[0]] = b;
        before[dae.charge_index[1]] = 0.25;
        after[dae.charge_index[1]] = 1.5;
        before[dae.charge_index[2]] = 1e8;
        after[dae.charge_index[2]] = 1e8;
        let effort = discrete_gradient(dae.system.storage(), &before, &after);
        assert_eq!(effort[dae.charge_index[0]].to_bits(), expected.to_bits());
        close(effort[dae.charge_index[1]], 73.0 / 40.0, 1.0);
        assert_eq!(effort[dae.charge_index[2]], 1e8);
    }
}

#[test]
fn resolved_branch_work_does_not_require_a_representable_squared_state_increment() {
    let scale = 1e-180;
    let dae = two_caps(scale, scale);
    let mut before = vec![0.0; dae.system.state_dim()];
    let mut after = before.clone();
    before[dae.charge_index[0]] = 2.0 * scale;
    after[dae.charge_index[0]] = 0.5 * scale;
    before[dae.charge_index[1]] = 0.25 * scale;
    after[dae.charge_index[1]] = 1.5 * scale;
    assert_eq!(
        before
            .iter()
            .zip(&after)
            .map(|(a, b)| (b - a).powi(2))
            .sum::<f64>(),
        0.0
    );
    let effort = discrete_gradient(dae.system.storage(), &before, &after);
    close(effort[dae.charge_index[0]], 19.0 / 12.0, 1.0);
    close(effort[dae.charge_index[1]], 73.0 / 40.0, 1.0);
    close(
        (after[dae.charge_index[0]] - before[dae.charge_index[0]]) * effort[dae.charge_index[0]],
        -2.375 * scale,
        scale,
    );
    close(
        (after[dae.charge_index[1]] - before[dae.charge_index[1]]) * effort[dae.charge_index[1]],
        2.28125 * scale,
        scale,
    );
}

#[test]
fn prepared_analytic_and_reference_steps_use_the_original_native_branch_efforts() {
    let dae = two_caps(1.0, 1.0);
    let n = dae.system.state_dim();
    let [a, b, held] = <[usize; 3]>::try_from(dae.charge_index.as_slice()).unwrap();
    let mut resistance = vec![0.0; n * n];
    resistance[a * n + a] = 1.0;
    resistance[b * n + b] = 1.0;
    let mut ports = vec![0.0; 2 * n];
    ports[2 * a] = 1.0;
    ports[2 * b + 1] = 1.0;
    // Reuse the original source-owned circuit storage, with two ordinary
    // resistive current ports. Unused descriptor coordinates remain stationary.
    let system = PortHamiltonian::new(
        n,
        2,
        vec![0.0; n * n],
        resistance,
        ports,
        dae.system.storage().shared_clone().unwrap(),
    )
    .unwrap();
    let mut before = vec![0.0; n];
    before[a] = 2.0;
    before[b] = 0.25;
    let input = [0.0, 3.0];
    let expected = [
        (4.0 - 7.0_f64.sqrt()) / 3.0,
        (10.0 + 55.0_f64.sqrt()) / 12.0,
    ];
    let reference = step(&system, &before, &input, 1.0).unwrap();
    close(reference.x[a], expected[0], 1.0);
    close(reference.x[b], expected[1], 1.0);
    for (analytic, condensed) in [(false, false), (true, false), (true, true)] {
        let mut workspace = StepWorkspace::new(&system).unwrap();
        if condensed {
            workspace.set_condensed_pairs(&[[a, b]]).unwrap();
        }
        let fallback_calls = Cell::new(0);
        let hessian = |x: &[f64], direction: &[f64], out: &mut [f64]| {
            fallback_calls.set(fallback_calls.get() + 1);
            out.fill(0.0);
            out[a] = if x[a].abs() < 1.0 { 1.0 } else { 2.0 } * direction[a];
            out[b] = if x[b].abs() < 0.5 { 3.0 } else { 1.0 } * direction[b];
            out[held] = direction[held];
            true
        };
        let mut after = vec![0.0; n];
        let mut output = [0.0; 2];
        let record = if analytic {
            workspace.step_into_analytic(
                &system,
                &before,
                &input,
                1.0,
                &mut after,
                &mut output,
                &hessian,
            )
        } else {
            workspace.step_into(&system, &before, &input, 1.0, &mut after, &mut output)
        }
        .unwrap();
        for (actual, expected) in after.iter().zip(&reference.x) {
            close(*actual, *expected, 1.0);
        }
        for (actual, expected) in output.iter().zip(&reference.y) {
            close(*actual, *expected, 1.0);
        }
        let effort = workspace.last_step_discrete_effort().unwrap();
        close(
            effort[a] * (after[a] - before[a]),
            energy_a(after[a]) - energy_a(before[a]),
            1.0,
        );
        close(
            effort[b] * (after[b] - before[b]),
            energy_b(after[b]) - energy_b(before[b]),
            1.0,
        );
        close(
            record.dissipated,
            effort[a] * effort[a] + effort[b] * effort[b],
            1.0,
        );
        close(record.supplied, 3.0 * effort[b], 1.0);
        close(record.balance_residual(), 0.0, 1.0);
        if analytic {
            assert_eq!(
                fallback_calls.get(),
                0,
                "native constitutive tangent must replace the Gonzalez/Hessian fallback"
            );
        }
    }
}

#[test]
fn native_endpoint_tangent_matches_each_secant_derivative_and_independent_perturbations() {
    let dae = two_caps(1.0, 1.0);
    let n = dae.system.state_dim();
    let [a, b, held] = <[usize; 3]>::try_from(dae.charge_index.as_slice()).unwrap();
    let rule = dae.system.storage().discrete_gradient_rule().unwrap();
    let mut before = vec![0.0; n];
    let mut after = before.clone();
    before[a] = 2.0;
    after[a] = 0.5;
    before[b] = 0.25;
    after[b] = 1.5;
    let mut direction = vec![0.1; n];
    direction[a] = 0.7;
    direction[b] = -0.3;
    direction[held] = 0.125;
    let mut action = vec![f64::NAN; n];
    assert!(rule.end_jacobian_action(&before, &after, &direction, &mut action));
    // d/db [(H(b)-H(a))/(b-a)] = [V(b)-mean_voltage]/(b-a).
    // These exact rational derivatives use the independent piecewise fixture.
    close(action[a], (13.0 / 18.0) * direction[a], 1.0);
    close(action[b], (27.0 / 50.0) * direction[b], 1.0);
    close(action[held], 0.5 * direction[held], 1.0);
    let h = 1e-6;
    let plus: Vec<_> = after
        .iter()
        .zip(&direction)
        .map(|(x, d)| x + h * d)
        .collect();
    let minus: Vec<_> = after
        .iter()
        .zip(&direction)
        .map(|(x, d)| x - h * d)
        .collect();
    let plus = discrete_gradient(dae.system.storage(), &before, &plus);
    let minus = discrete_gradient(dae.system.storage(), &before, &minus);
    for i in 0..n {
        assert!((action[i] - (plus[i] - minus[i]) / (2.0 * h)).abs() < 2e-9);
        if ![a, b, held].contains(&i) {
            assert_eq!(
                action[i], 0.0,
                "unused algebraic coordinates have no stored effort"
            );
        }
    }
    before[a] = 0.25;
    before[b] = 0.25;
    after.copy_from_slice(&before);
    assert!(rule.end_jacobian_action(&before, &after, &direction, &mut action));
    close(action[a], 0.5 * direction[a], 1.0);
    close(action[b], 1.5 * direction[b], 1.0);
    before[a] = 1.0;
    after[a] = f64::from_bits(1.0_f64.to_bits() + 1);
    assert!(rule.end_jacobian_action(&before, &after, &direction, &mut action));
    close(action[a], direction[a], 1.0);
}
