//! Coupled material storage uses the existing descriptor structure and solver.

use fs_phs::circuit::{Branch, CircuitGraph, assemble_circuit};
use fs_phs::{
    DescriptorPortHamiltonian, QuadraticStorage, Storage, discrete_gradient, step_descriptor,
};

struct QuarticCapacitor;

impl Storage for QuarticCapacitor {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        assert_eq!(x.len(), 1, "the law must never receive a multiplier");
        x[0].powi(4) / 4.0
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        assert_eq!((x.len(), out.len()), (1, 1));
        out[0] = x[0].powi(3);
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        assert_eq!((a.len(), b.len(), mid.len()), (1, 1, 1));
        Some(mid[0] * (b[0] - a[0]).powi(3) / 4.0)
    }
}

fn nonlinear_rc(current_scale: f64) -> DescriptorPortHamiltonian {
    // q_dot = i, 0 = voltage_source - q^3 - i. The multiplier is
    // current_scale*i; changing its coordinates must not change physics.
    DescriptorPortHamiltonian::new(
        2,
        1,
        1,
        vec![0.0, 1.0 / current_scale, -1.0 / current_scale, 0.0],
        vec![0.0, 0.0, 0.0, 1.0 / current_scale.powi(2)],
        vec![0.0, 1.0 / current_scale],
        Box::new(QuarticCapacitor),
    )
    .unwrap()
}

#[test]
fn nonlinear_descriptor_keeps_the_complete_energy_remainder_on_charge() {
    let q0 = 0.5_f64;
    let voltage = 2.0;
    let dt = 0.3;
    // Independent scalar reference: the quartic energy divided difference
    // is (q1^3 + q1^2*q0 + q1*q0^2 + q0^3)/4, not midpoint q^3.
    let mut lower = q0;
    let mut upper = 2.0;
    for _ in 0..60 {
        let q1 = f64::midpoint(lower, upper);
        let effort = (q1.powi(3) + q1.powi(2) * q0 + q1 * q0.powi(2) + q0.powi(3)) / 4.0;
        if q1 - q0 - dt * (voltage - effort) > 0.0 {
            upper = q1;
        } else {
            lower = q1;
        }
    }
    let expected_charge = f64::midpoint(lower, upper);
    let expected_current = (expected_charge - q0) / dt;
    for current_scale in [1.0, 8.0] {
        let system = nonlinear_rc(current_scale);
        let before = [q0, current_scale * (voltage - q0.powi(3))];
        let record = step_descriptor(&system, &before, &[voltage], dt).unwrap();
        assert!((record.x[0] - expected_charge).abs() < 2e-12);
        assert!((record.y[0] - expected_current).abs() < 2e-12);
        assert!((record.dissipated - dt * expected_current.powi(2)).abs() < 5e-12);
        assert!((record.supplied - dt * voltage * expected_current).abs() < 5e-12);
        assert!((record.delta_h + record.dissipated - record.supplied).abs() < 5e-12);

        // Reproduce the old defect independently: forming Gonzalez over charge
        // AND the changing current multiplier assigns energy to the latter.
        // Replacing that algebraic effort then discards this nonzero remainder.
        let wrong = discrete_gradient(system.storage(), &before, &record.x);
        let lost = record.delta_h - wrong[0] * (record.x[0] - q0);
        assert!(lost.abs() > 1e-4);
    }
}

#[test]
fn storage_is_differential_only_and_quadratic_metadata_has_zero_algebraic_block() {
    let nonlinear = nonlinear_rc(1.0);
    assert_eq!(nonlinear.hamiltonian(&[2.0, 1e150]), 4.0);
    let mut gradient = [7.0, 7.0];
    nonlinear.storage().gradient(&[2.0, -1e150], &mut gradient);
    assert_eq!(gradient, [8.0, 0.0]);
    let quadratic = DescriptorPortHamiltonian::new(
        2,
        1,
        0,
        vec![0.0; 4],
        vec![0.0; 4],
        vec![],
        Box::new(QuadraticStorage::new(vec![3.0], 1).unwrap()),
    )
    .unwrap();
    let expected = [3.0, 0.0, 0.0, 0.0];
    assert_eq!(
        quadratic.storage().exact_quadratic_matrix().unwrap(),
        &expected
    );
    let part = quadratic.storage().quadratic_part().unwrap();
    assert!(part.exact);
    assert_eq!(part.matrix, expected);
    assert_eq!(quadratic.hamiltonian(&[2.0, 1e150]), 6.0);
}

#[test]
fn reusing_native_rlc_stamps_and_storage_preserves_ports_and_trajectory() {
    let circuit = assemble_circuit(&CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (2, 0, Branch::Capacitor { farads: 3.0 }),
            (1, 2, Branch::Resistor { ohms: 4.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
        ],
        transformers: vec![],
    })
    .unwrap();
    let native = &circuit.system;
    let n = native.state_dim();
    let nd = native.differential_dim();
    let q = native.storage().exact_quadratic_matrix().unwrap();
    let mut differential_q = vec![0.0; nd * nd];
    for row in 0..nd {
        differential_q[row * nd..(row + 1) * nd].copy_from_slice(&q[row * n..row * n + nd]);
    }
    let (j, r, g) = native.structure();
    let rebuilt = DescriptorPortHamiltonian::new(
        n,
        nd,
        native.port_dim(),
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(QuadraticStorage::new(differential_q, nd).unwrap()),
    )
    .unwrap();
    assert_eq!(native.structure(), rebuilt.structure());
    let ports = [0.15];
    let mut state = vec![0.0; n];
    state[circuit.flux_index[0]] = 0.25;
    state[circuit.charge_index[0]] = 0.4;
    let initial = circuit.consistent_initial_state(&state, &ports).unwrap();
    let mut a = initial.clone();
    let mut b = initial;
    for _ in 0..8 {
        let original = circuit.step_audited(&a, &ports, 0.1).unwrap().0;
        let coupled = step_descriptor(&rebuilt, &b, &ports, 0.1).unwrap();
        for (left, right) in original.x.iter().zip(&coupled.x) {
            assert!((left - right).abs() < 2e-12);
        }
        assert!((original.y[0] - coupled.y[0]).abs() < 2e-12);
        assert!((original.dissipated - coupled.dissipated).abs() < 2e-12);
        assert!((coupled.delta_h + coupled.dissipated - coupled.supplied).abs() < 2e-12);
        a = original.x;
        b = coupled.x;
    }
}

#[test]
fn descriptor_constructor_rejects_invalid_partition_and_structure() {
    let make =
        |nd, j, r, g| DescriptorPortHamiltonian::new(2, nd, 1, j, r, g, Box::new(QuarticCapacitor));
    assert!(make(3, vec![0.0; 4], vec![0.0; 4], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0, 1.0, 1.0, 0.0], vec![0.0; 4], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0; 4], vec![0.0, 0.0, 0.0, -1.0], vec![0.0; 2]).is_err());
    assert!(make(1, vec![0.0; 4], vec![0.0; 4], vec![f64::NAN, 0.0]).is_err());
    assert!(
        DescriptorPortHamiltonian::new(
            usize::MAX,
            1,
            1,
            vec![],
            vec![],
            vec![],
            Box::new(QuarticCapacitor),
        )
        .is_err()
    );
}
