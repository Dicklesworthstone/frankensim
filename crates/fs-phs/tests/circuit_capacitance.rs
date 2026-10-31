//! G0/G1 nonlinear capacitance: independent charge/voltage, exact piecewise RC
//! transient, and source-free nonlinear LC energy on the original native DAE.

use fs_phs::circuit::capacitance::{
    CapacitorChargeCurve, CapacitorCurve, assemble_circuit_with_capacitor_curves,
};
use fs_phs::circuit::{Branch, CircuitDae, CircuitGraph, MutualInductance, assemble_circuit};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}

fn curve(scale: f64) -> CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(vec![
        (-4.0 * scale, -7.0),
        (-scale, -1.0),
        (0.0, 0.0),
        (scale, 1.0),
        (4.0 * scale, 7.0),
    ])
    .unwrap()
}

fn binding(branch: usize, scale: f64) -> CapacitorCurve {
    CapacitorCurve {
        branch,
        curve: curve(scale),
    }
}

fn initial(dae: &CircuitDae, charges: &[f64], fluxes: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0; dae.system.state_dim()];
    for (&i, &q) in dae.charge_index.iter().zip(charges) {
        x[i] = q;
    }
    for (&i, &f) in dae.flux_index.iter().zip(fluxes) {
        x[i] = f;
    }
    dae.consistent_initial_state(&x, &[]).unwrap()
}

#[test]
fn differential_capacitance_is_not_secant_capacitance_or_half_qv_energy() {
    let law = curve(1.0);
    near(law.voltage_v(2.0).unwrap(), 3.0, 0.0);
    near(law.differential_capacitance_f(2.0).unwrap(), 0.5, 0.0);
    near(law.stored_energy_j(2.0).unwrap(), 2.5, 0.0);
    assert_ne!(law.stored_energy_j(2.0).unwrap(), 0.5 * 2.0 * 3.0);
    for q in [-3.4, -1.0, -1e-20, 0.0, 1e-20, 0.8, 2.0, 4.0] {
        near(
            law.charge_c(law.voltage_v(q).unwrap()).unwrap(),
            q,
            2e-15 * q.abs(),
        );
        near(
            law.stored_energy_j(-q).unwrap(),
            law.stored_energy_j(q).unwrap(),
            2e-15,
        );
    }
    near(law.stored_energy_j(-1e-20).unwrap(), 0.5e-40, 1e-55);
    near(law.energy_change_j(-2.0, 3.0).unwrap(), 4.0, 1e-15);
    near(law.energy_change_j(3.0, -2.0).unwrap(), -4.0, 1e-15);
}

fn discharge(steps: usize, capacitance_scale: f64) -> f64 {
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (
                1,
                0,
                Branch::Resistor {
                    ohms: 1.0 / capacitance_scale,
                },
            ),
            (
                1,
                0,
                Branch::Capacitor {
                    farads: capacitance_scale,
                },
            ),
        ],
        transformers: vec![],
    };
    let dae = assemble_circuit_with_capacitor_curves(&graph, &[binding(1, capacitance_scale)], &[])
        .unwrap();
    let mut x = initial(&dae, &[2.0 * capacitance_scale], &[]);
    let start = dae.system.hamiltonian(&x);
    let mut heat = 0.0;
    for _ in 0..steps {
        let (record, defect) = dae.step_audited(&x, &[], 1.0 / steps as f64).unwrap();
        assert!(record.dissipated > 0.0);
        assert!(record.x[0] < x[0]);
        assert!(defect < 2e-10 * capacitance_scale, "defect={defect:e}");
        heat += record.dissipated;
        x = record.x;
    }
    near(
        dae.system.hamiltonian(&x) + heat,
        start,
        5e-10 * capacitance_scale,
    );
    x[dae.charge_index[0]] / capacitance_scale
}

#[test]
fn nonlinear_rc_crosses_a_capacitance_knot_with_second_order_accuracy() {
    // V(q)=2q-1 above q=1, then V(q)=q. Starting q=2 with R=1,
    // q=.5+1.5*exp(-2t) until t=ln(3)/2, then exp(-(t-ln(3)/2)).
    let exact = 3.0_f64.sqrt() * (-1.0_f64).exp();
    let errors = [32, 64, 128].map(|n| (discharge(n, 1.0) - exact).abs());
    assert!(errors[2] < 1e-5, "{errors:?}");
    assert!(
        errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0,
        "{errors:?}"
    );
}

#[test]
fn realistic_small_capacitance_retains_the_same_dimensionless_rc_response() {
    let expected = discharge(64, 1.0);
    for scale in [1e-9, 1e-12, 1e-15] {
        near(discharge(64, scale), expected, 2e-7);
    }
}

#[test]
fn uncharged_small_capacitors_use_the_actual_finite_charge_support() {
    for scale in [1e-9, 1e-12, 1e-15] {
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: scale }),
                (1, 0, Branch::CurrentSource { port: 0 }),
            ],
            transformers: vec![],
        };
        let dae =
            assemble_circuit_with_capacitor_curves(&graph, &[binding(0, scale)], &[]).unwrap();
        let before = vec![0.0, 0.0, 2.0 * scale];
        let (step, defect) = dae.step_audited(&before, &[2.0 * scale], 1.0).unwrap();
        near(step.x[0] / scale, 2.0, 2e-12);
        near(dae.capacitor_voltages(&step.x).unwrap()[0], 3.0, 3e-12);
        near(step.delta_h / scale, 2.5, 3e-12);
        near(step.supplied / scale, 2.5, 3e-12);
        assert!(defect / scale < 3e-12);
    }
}

#[test]
fn picofarad_knot_crossing_keeps_local_accuracy_beside_a_large_held_bias() {
    let scale = 1e-12;
    let law = CapacitorChargeCurve::try_new(vec![
        (-4.0 * scale, -13.0),
        (-scale, -1.0),
        (0.0, 0.0),
        (scale, 1.0),
        (4.0 * scale, 13.0),
    ])
    .unwrap();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::Capacitor { farads: scale }),
            (2, 0, Branch::Resistor { ohms: 1e6 }),
        ],
        transformers: vec![],
    };
    let dae = assemble_circuit_with_capacitor_curves(
        &graph,
        &[CapacitorCurve {
            branch: 2,
            curve: law,
        }],
        &[],
    )
    .unwrap();
    let mut x = vec![0.0; dae.system.state_dim()];
    x[dae.charge_index[0]] = 1e8;
    x[dae.charge_index[1]] = 2.0 * scale;
    let (record, _) = dae.step_audited(&x, &[1e8], 0.5e-6).unwrap();

    // V(z)=z below the knot and V(z)=4z-3 above it, z=q/scale.
    // H(2)/scale=3.5. Combining the exact voltage integral with RC
    // continuity gives 5*z1^2-16*z1+9=0 across this single crossing.
    let expected_charge = (8.0 - 19.0_f64.sqrt()) / 5.0;
    assert!(expected_charge > 0.0 && expected_charge < 1.0);
    near(
        record.x[dae.charge_index[1]],
        scale * expected_charge,
        1e-24,
    );
    assert_eq!(
        record.x[dae.charge_index[0]].to_bits(),
        x[dae.charge_index[0]].to_bits()
    );
    assert_eq!(
        dae.system.hamiltonian(&x),
        dae.system.hamiltonian(&record.x)
    );
    let expected_change = scale * (0.5 * expected_charge * expected_charge - 3.5);
    near(record.delta_h, expected_change, 1e-23);
    near(record.dissipated, -expected_change, 1e-22);
    near(record.supplied, 0.0, 0.0);
    near(record.balance_residual(), 0.0, 1e-22);
}

#[test]
fn nonlinear_lc_oscillates_across_both_polarities_without_artificial_loss() {
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::Inductor { henries: 1.0 }),
        ],
        transformers: vec![],
    };
    let dae = assemble_circuit_with_capacitor_curves(&graph, &[binding(0, 1.0)], &[]).unwrap();
    let mut x = initial(&dae, &[2.0], &[0.0]);
    let energy = dae.system.hamiltonian(&x);
    let mut minimum = 2.0_f64;
    for _ in 0..200 {
        let (record, defect) = dae.step_audited(&x, &[], 0.05).unwrap();
        near(record.dissipated, 0.0, 0.0);
        near(record.supplied, 0.0, 0.0);
        assert!(defect < 2e-10);
        near(dae.system.hamiltonian(&record.x), energy, 3e-9);
        minimum = minimum.min(record.x[dae.charge_index[0]]);
        x = record.x;
    }
    assert!(minimum < -1.99);
}

#[test]
fn charge_laws_preserve_interleaved_branch_order_and_mutual_inductance() {
    let graph = CircuitGraph {
        node_count: 4,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (3, 0, Branch::Capacitor { farads: 99.0 }),
            (2, 0, Branch::Inductor { henries: 0.5 }),
            (3, 0, Branch::Capacitor { farads: 4.0 }),
            (1, 0, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Resistor { ohms: 1.0 }),
        ],
        transformers: vec![],
    };
    let mutual = [MutualInductance {
        branches: (0, 2),
        henries: 0.5,
    }];
    let dae = assemble_circuit_with_capacitor_curves(&graph, &[binding(1, 1.0)], &mutual).unwrap();
    let mut x = vec![0.0; dae.system.state_dim()];
    assert_eq!(dae.flux_index, [0, 1]);
    assert_eq!(dae.charge_index, [2, 3]);
    x[0] = 1.8;
    x[1] = 0.3;
    x[2] = 2.0;
    x[3] = 12.0;
    let currents = dae.inductor_currents(&x).unwrap();
    near(currents[0], 1.0, 2e-15);
    near(currents[1], -0.4, 2e-15);
    assert_eq!(dae.capacitor_voltages(&x).unwrap(), [3.0, 3.0]);
    near(dae.system.hamiltonian(&x), 0.84 + 2.5 + 18.0, 2e-14);
    assert_eq!(
        dae.system.dirac_j(),
        assemble_circuit(&graph).unwrap().system.dirac_j()
    );
}

#[test]
fn small_energy_change_on_a_large_dc_bias_keeps_the_crossed_knot() {
    let knee = 1e6;
    let law = CapacitorChargeCurve::try_new(vec![
        (-2e6, -2e6),
        (0.0, 0.0),
        (knee, knee),
        (knee + 1.0, knee + 2.0),
    ])
    .unwrap();
    let before = knee - 1e-4;
    let after = knee + 3e-4;
    let left = before - knee;
    let right = after - knee;
    let expected = knee * (right - left) + right * right - 0.5 * left * left;
    near(law.energy_change_j(before, after).unwrap(), expected, 2e-13);
    near(
        law.energy_change_j(after, before).unwrap(),
        -expected,
        2e-13,
    );
    let subtraction = law.stored_energy_j(after).unwrap() - law.stored_energy_j(before).unwrap();
    assert!((subtraction - expected).abs() > 1e-7);
}

#[test]
fn malformed_curves_unknown_branches_and_out_of_support_steps_refuse() {
    for knots in [
        vec![(-1.0, -1.0), (0.0, 0.0), (1.0, -1.0)],
        vec![(-1.0, -1.0), (0.0, 0.0), (1.0, 0.0)],
        vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)],
        vec![(-1.0, -1.0), (0.0, 0.0), (f64::NAN, 1.0)],
    ] {
        assert!(CapacitorChargeCurve::try_new(knots).is_err());
    }
    let law = curve(1.0);
    assert!(law.voltage_v(4.000001).is_err());
    assert!(law.charge_c(7.000001).is_err());
    assert!(law.stored_energy_j(f64::INFINITY).is_err());
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    };
    assert!(assemble_circuit_with_capacitor_curves(&graph, &[binding(1, 1.0)], &[]).is_err());
    assert!(
        assemble_circuit_with_capacitor_curves(&graph, &[binding(0, 1.0), binding(0, 1.0)], &[])
            .is_err()
    );
    let dae = assemble_circuit_with_capacitor_curves(&graph, &[binding(0, 1.0)], &[]).unwrap();
    let mut x = vec![0.0; dae.system.state_dim()];
    x[dae.charge_index[0]] = 3.9;
    let before = x.clone();
    assert!(dae.step_audited(&x, &[100.0], 1.0).is_err());
    assert_eq!(x, before);
    x[dae.charge_index[0]] = 4.01;
    assert!(dae.capacitor_voltages(&x).is_err());
}
