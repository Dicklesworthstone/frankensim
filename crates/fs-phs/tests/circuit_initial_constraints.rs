//! Initial constraints retain stored coordinates even with redundant ideal laws.

use fs_phs::circuit::{Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit};
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, step_descriptor};

fn clamped_rlc() -> CircuitDae {
    assemble_circuit(&CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Capacitor { farads: 1e-6 }),
            (1, 0, Branch::Inductor { henries: 0.05 }),
            (1, 0, Branch::Resistor { ohms: 10.0 }),
        ],
        transformers: vec![],
    })
    .unwrap()
}

#[test]
fn compatible_ideal_voltage_clamp_keeps_charge_and_flux_from_any_multiplier_guess() {
    let circuit = clamped_rlc();
    let voltage = 12.0;
    let flux = 0.02;
    let charge = 1e-6 * voltage;
    let nd = circuit.system.differential_dim();
    let node = circuit.node_potential_index[0];
    for guess in [[0.0, 0.0, 0.0], [-7.0, 2.0, -4.0]] {
        let mut before = vec![0.0; circuit.system.state_dim()];
        before[circuit.flux_index[0]] = flux;
        before[circuit.charge_index[0]] = charge;
        before[nd..].copy_from_slice(&guess);
        let after = circuit
            .consistent_initial_state(&before, &[voltage])
            .unwrap();
        for i in 0..nd {
            assert_eq!(after[i].to_bits(), before[i].to_bits());
        }
        assert_eq!(
            circuit.system.hamiltonian(&after),
            circuit.system.hamiltonian(&before)
        );
        // The existing IC result represents endpoint multipliers, so the
        // actual constraint effort is their midpoint with the caller's guess.
        let v = f64::midpoint(before[node], after[node]);
        let i_cap = f64::midpoint(before[node + 1], after[node + 1]);
        let i_source = f64::midpoint(before[node + 2], after[node + 2]);
        assert!((v - voltage).abs() < 2e-11);
        assert!((v - charge / 1e-6).abs() < 2e-11);
        assert!((i_cap + i_source + flux / 0.05 + v / 10.0).abs() < 2e-11);

        // A real positive-time step determines the formerly free current
        // split: the held voltage fixes capacitor charge, while dphi/dt=V.
        let dt = 1e-3;
        let step = circuit.step_audited(&after, &[voltage], dt).unwrap().0;
        let expected_flux = flux + dt * voltage;
        let expected_source_current = voltage / 10.0 + f64::midpoint(flux, expected_flux) / 0.05;
        assert!((step.x[circuit.flux_index[0]] - expected_flux).abs() < 2e-11);
        assert!((step.x[circuit.charge_index[0]] - charge).abs() < 2e-13);
        assert!((step.y[0] - expected_source_current).abs() < 2e-9);
        assert!((step.dissipated - dt * voltage * voltage / 10.0).abs() < 2e-10);
        assert!(step.balance_residual().abs() < 2e-10);
    }
}

#[test]
fn redundant_voltage_constraint_does_not_project_incompatible_charge() {
    let circuit = clamped_rlc();
    let mut before = vec![0.0; circuit.system.state_dim()];
    before[circuit.flux_index[0]] = 0.02;
    before[circuit.charge_index[0]] = 5e-6;
    let result = circuit.consistent_initial_state(&before, &[12.0]);
    match result {
        Err(CircuitError::InconsistentInitialConditions { residual }) => {
            assert!(residual.is_finite() && residual > 1.0);
        }
        other => panic!("5-volt stored charge cannot satisfy a 12-volt clamp: {other:?}"),
    }
}

#[test]
fn source_free_lc_stamps_have_exact_empty_port_shape_and_rebuild_unchanged() {
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (1, 0, Branch::Capacitor { farads: 3.0 }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let native = &circuit.system;
    let n = native.state_dim();
    let nd = native.differential_dim();
    let (j, r, g) = native.structure();
    assert_eq!(native.port_dim(), 0);
    assert!(g.is_empty());
    let full_q = native.storage().exact_quadratic_matrix().unwrap();
    let mut q = vec![0.0; nd * nd];
    for row in 0..nd {
        q[row * nd..(row + 1) * nd].copy_from_slice(&full_q[row * n..row * n + nd]);
    }
    let rebuilt = DescriptorPortHamiltonian::new(
        n,
        nd,
        0,
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(QuadraticStorage::new(q, nd).unwrap()),
    )
    .unwrap();
    let mut before = vec![0.0; n];
    before[circuit.flux_index[0]] = 0.25;
    before[circuit.charge_index[0]] = 0.4;
    let mut original = circuit.consistent_initial_state(&before, &[]).unwrap();
    let mut copy = original.clone();
    for _ in 0..4 {
        let a = circuit.step_audited(&original, &[], 0.05).unwrap().0;
        let b = step_descriptor(&rebuilt, &copy, &[], 0.05).unwrap();
        for (x, y) in a.x.iter().zip(&b.x) {
            assert!((x - y).abs() < 2e-12);
        }
        assert!(a.y.is_empty() && b.y.is_empty());
        assert_eq!(a.supplied, 0.0);
        assert_eq!(b.supplied, 0.0);
        assert!(a.balance_residual().abs() < 2e-12);
        assert!(b.balance_residual().abs() < 2e-12);
        original = a.x;
        copy = b.x;
    }
}

#[test]
fn finite_voltage_onset_admits_algebraic_currents_without_creating_stored_charge_or_flux() {
    // Original branch order deliberately interleaves C/L, with two unused
    // source ports. A nonzero source can establish immediate capacitor/source
    // currents while every differential coordinate remains exactly zero.
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 4,
        branches: vec![
            (3, 0, Branch::Capacitor { farads: 0.4 }),
            (2, 3, Branch::Inductor { henries: 0.7 }),
            (1, 2, Branch::Resistor { ohms: 0.8 }),
            (1, 0, Branch::VoltageSource { port: 2 }),
            (2, 0, Branch::Capacitor { farads: 0.3 }),
            (3, 0, Branch::Inductor { henries: 4.5 }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let before = vec![0.0; circuit.system.state_dim()];
    let ports = [0.0, 0.0, 0.03];
    let after = circuit.consistent_initial_state(&before, &ports).unwrap();
    for index in 0..circuit.system.differential_dim() {
        assert_eq!(after[index].to_bits(), before[index].to_bits());
    }
    assert_eq!(circuit.system.hamiltonian(&after), 0.0);
    for (&node, expected) in circuit.node_potential_index.iter().zip([0.03, 0.0, 0.0]) {
        assert!((f64::midpoint(before[node], after[node]) - expected).abs() < 2e-12);
    }
    let step = circuit.step_audited(&after, &ports, 0.02).unwrap().0;
    assert!(step.supplied > 0.0 && step.delta_h > 0.0 && step.dissipated > 0.0);
    assert!(step.balance_residual().abs() < 2e-12);
}
