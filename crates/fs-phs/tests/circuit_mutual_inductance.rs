//! G0/G1: real coupled windings retain Faraday polarity, leakage dynamics and
//! the native Kirchhoff storage/supply/heat balance.

use fs_phs::circuit::{
    Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance, assemble_circuit,
    assemble_circuit_with_mutual_inductance,
};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn loops() -> CircuitGraph {
    CircuitGraph {
        node_count: 4,
        // Charge storage and resistors deliberately interleave the windings.
        branches: vec![
            (3, 0, Branch::Capacitor { farads: 4.0 }),
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (2, 0, Branch::Resistor { ohms: 2.0 }),
            (2, 0, Branch::Inductor { henries: 0.5 }),
            (1, 0, Branch::Resistor { ohms: 1.0 }),
            (3, 0, Branch::Resistor { ohms: 3.0 }),
        ],
        transformers: vec![],
    }
}

fn link(m: f64) -> MutualInductance {
    MutualInductance {
        branches: (1, 3),
        henries: m,
    }
}

fn initial(circuit: &CircuitDae, flux: [f64; 2], charge: f64) -> Vec<f64> {
    let mut state = vec![0.0; circuit.system.state_dim()];
    state[circuit.flux_index[0]] = flux[0];
    state[circuit.flux_index[1]] = flux[1];
    if let Some(&i) = circuit.charge_index.first() {
        state[i] = charge;
    }
    circuit.consistent_initial_state(&state, &[]).unwrap()
}

#[test]
fn full_mutual_storage_drives_original_branches_and_capacitors() {
    let circuit = assemble_circuit_with_mutual_inductance(&loops(), &[link(0.5)]).unwrap();
    assert_eq!(circuit.flux_index, [0, 1]);
    assert_eq!(circuit.charge_index, [2]);
    let currents = [1.0, -0.4];
    let charge = 0.8;
    let before = initial(&circuit, [1.8, 0.3], charge); // L * currents
    for (&i, expected) in circuit
        .inductor_currents(&before)
        .unwrap()
        .iter()
        .zip(currents)
    {
        near(i, expected);
    }
    let energy = 0.5 * (2.0 + 2.0 * 0.5 * -0.4 + 0.5 * 0.16) + charge * charge / 8.0;
    near(circuit.system.hamiltonian(&before), energy);

    let dt = 0.1;
    // Independent midpoint KVL equation (L+dt R/2)*i1=(L-dt R/2)*i0.
    let a = 2.0 + dt / 2.0;
    let d = 0.5 + dt;
    let rhs = [1.8 - dt * currents[0] / 2.0, 0.3 - dt * currents[1]];
    let determinant = a * d - 0.25;
    let expected = [
        (d * rhs[0] - 0.5 * rhs[1]) / determinant,
        (a * rhs[1] - 0.5 * rhs[0]) / determinant,
    ];
    let expected_q = charge * (1.0 - dt / 24.0) / (1.0 + dt / 24.0);
    let (record, defect) = circuit.step_audited(&before, &[], dt).unwrap();
    for (&i, expected) in circuit
        .inductor_currents(&record.x)
        .unwrap()
        .iter()
        .zip(expected)
    {
        near(i, expected);
    }
    near(record.x[circuit.charge_index[0]], expected_q);
    let mid_i = [0, 1].map(|i| f64::midpoint(currents[i], expected[i]));
    let cap_voltage = f64::midpoint(charge, expected_q) / 4.0;
    near(
        record.dissipated,
        dt * (mid_i[0].powi(2) + 2.0 * mid_i[1].powi(2) + cap_voltage.powi(2) / 3.0),
    );
    near(record.supplied, 0.0);
    near(record.delta_h, -record.dissipated);
    assert!(defect < 2e-12);
}

#[test]
fn reversing_a_winding_and_its_dot_keeps_the_same_physical_response() {
    let original = assemble_circuit_with_mutual_inductance(&loops(), &[link(0.5)]).unwrap();
    let mut reversed = loops();
    reversed.branches[3].0 = 0;
    reversed.branches[3].1 = 2;
    let reversed = assemble_circuit_with_mutual_inductance(&reversed, &[link(-0.5)]).unwrap();
    let mut a = initial(&original, [1.8, 0.3], 0.8);
    let mut b = initial(&reversed, [1.8, -0.3], 0.8);
    for _ in 0..20 {
        let first = original.step_audited(&a, &[], 0.05).unwrap().0;
        let second = reversed.step_audited(&b, &[], 0.05).unwrap().0;
        let ia = original.inductor_currents(&first.x).unwrap();
        let ib = reversed.inductor_currents(&second.x).unwrap();
        near(ia[0], ib[0]);
        near(ia[1], -ib[1]);
        near(first.delta_h, second.delta_h);
        near(first.dissipated, second.dissipated);
        for &node in &original.node_potential_index {
            near(first.x[node], second.x[node]);
        }
        a = first.x;
        b = second.x;
    }
}

fn transformer(dt: f64, steps: usize) -> [f64; 2] {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Inductor { henries: 4.0 }),
            (2, 0, Branch::Resistor { ohms: 2.0 }),
            (2, 0, Branch::Inductor { henries: 1.0 }),
        ],
        transformers: vec![],
    };
    let circuit = assemble_circuit_with_mutual_inductance(&graph, &[link(1.5)]).unwrap();
    let mut state = circuit
        .consistent_initial_state(&vec![0.0; circuit.system.state_dim()], &[1.0])
        .unwrap();
    let mut source_work = 0.0;
    let mut heat = 0.0;
    for _ in 0..steps {
        let (record, defect) = circuit.step_audited(&state, &[1.0], dt).unwrap();
        assert!(defect < 2e-12);
        source_work += record.supplied;
        heat += record.dissipated;
        state = record.x;
    }
    near(state[circuit.flux_index[0]], dt * steps as f64); // Faraday: dλ1/dt = V.
    near(source_work, circuit.system.hamiltonian(&state) + heat);
    let currents = circuit.inductor_currents(&state).unwrap();
    [currents[0], currents[1]]
}

#[test]
fn loaded_transformer_has_finite_magnetizing_and_leakage_transients() {
    // V1=1, R2=2: tau=(L2-M²/L1)/R2; secondary current is opposed to
    // primary flux growth. An ideal transformer cannot produce this transient.
    let time = 0.2;
    let tau = (1.0 - 1.5 * 1.5 / 4.0) / 2.0;
    let i2 = -1.5 / 8.0 * (1.0 - f64::exp(-time / tau));
    let exact = [(time - 1.5 * i2) / 4.0, i2];
    let mut previous_error = None;
    for steps in [20, 40, 80] {
        let actual = transformer(time / steps as f64, steps);
        let error = (actual[1] - exact[1]).abs();
        if let Some(previous) = previous_error {
            assert!(previous / error > 3.9 && previous / error < 4.1);
        }
        previous_error = Some(error);
        assert!((actual[0] - exact[0]).abs() < 5e-6);
    }
    assert!(previous_error.unwrap() < 1e-6);
}

#[test]
fn every_pair_can_be_admissible_while_the_full_core_has_negative_energy() {
    let graph = CircuitGraph {
        node_count: 4,
        branches: (1..=3)
            .map(|node| (node, 0, Branch::Inductor { henries: 1.0 }))
            .collect(),
        transformers: vec![],
    };
    let links = [
        MutualInductance {
            branches: (0, 1),
            henries: 0.9,
        },
        MutualInductance {
            branches: (0, 2),
            henries: 0.9,
        },
        MutualInductance {
            branches: (1, 2),
            henries: -0.9,
        },
    ];
    // i=(-1,1,1) gives i^T L i = -2.4 despite all three |M| < sqrt(L_i L_j).
    assert!(matches!(
        assemble_circuit_with_mutual_inductance(&graph, &links),
        Err(CircuitError::Invalid {
            what: "full inductance matrix must be positive definite"
        })
    ));
}

#[test]
fn malformed_and_singular_couplings_refuse_instead_of_inventing_leakage() {
    let graph = loops();
    for links in [
        vec![link(1.0)], // |M|=sqrt(2*0.5): dependent flux coordinates.
        vec![link(1.1)],
        vec![link(f64::NAN)],
        vec![link(f64::INFINITY)],
        vec![MutualInductance {
            branches: (1, 1),
            henries: 0.2,
        }],
        vec![MutualInductance {
            branches: (1, 0),
            henries: 0.2,
        }], // capacitor
        vec![MutualInductance {
            branches: (1, 99),
            henries: 0.2,
        }],
        vec![
            link(0.2),
            MutualInductance {
                branches: (3, 1),
                henries: 0.2,
            },
        ],
    ] {
        assert!(assemble_circuit_with_mutual_inductance(&graph, &links).is_err());
    }
    assert!(assemble_circuit_with_mutual_inductance(&graph, &[link(0.99)]).is_ok());
}

#[test]
fn dimensionless_energy_admission_retains_unequal_winding_scales() {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1e-24 }),
            (2, 0, Branch::Inductor { henries: 1e24 }),
        ],
        transformers: vec![],
    };
    let circuit = assemble_circuit_with_mutual_inductance(
        &graph,
        &[MutualInductance {
            branches: (0, 1),
            henries: 0.5,
        }],
    )
    .unwrap();
    // Currents chosen to store order-one energy in each coil; the large
    // difference in turns/inductance must not erase the small winding.
    let mut state = vec![0.0; circuit.system.state_dim()];
    state[0] = 1.5e-12;
    state[1] = 1.5e12;
    let currents = circuit.inductor_currents(&state).unwrap();
    assert!((currents[0] / 1e12 - 1.0).abs() < 1e-12);
    assert!((currents[1] / 1e-12 - 1.0).abs() < 1e-12);
    near(circuit.system.hamiltonian(&state), 1.5);
}

#[test]
fn no_links_preserve_native_circuit_stamps_exactly_and_currents_validate_state() {
    let a = assemble_circuit(&loops()).unwrap();
    let b = assemble_circuit_with_mutual_inductance(&loops(), &[]).unwrap();
    assert_eq!(a.system.structure(), b.system.structure());
    assert_eq!(
        a.system.storage().exact_quadratic_matrix(),
        b.system.storage().exact_quadratic_matrix()
    );
    assert!(a.inductor_currents(&[]).is_err());
    assert!(
        a.inductor_currents(&vec![f64::NAN; a.system.state_dim()])
            .is_err()
    );
    let state = initial(&a, [2.0, -0.2], 0.0);
    assert_eq!(a.inductor_currents(&state).unwrap(), [1.0, -0.4]);
}
