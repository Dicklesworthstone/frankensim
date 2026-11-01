//! G0/G1/G3: full reciprocal capacitor matrices on the native Kirchhoff DAE.
//! Closed two-node equations and explicit quadratic energies are independent
//! of the production factorization, state assembly and time-step accounts.

use fs_phs::circuit::capacitance::{
    CapacitanceBlock, CapacitorChargeCurve, CapacitorCurve,
    assemble_circuit_with_capacitance_blocks, assemble_circuit_with_capacitor_curves,
};
use fs_phs::circuit::{Branch, CircuitDae, CircuitGraph, MutualInductance};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn initial(dae: &CircuitDae, charges: &[f64], fluxes: &[f64]) -> Vec<f64> {
    let mut state = vec![0.0; dae.system.state_dim()];
    for (&index, &value) in dae.charge_index.iter().zip(charges) {
        state[index] = value;
    }
    for (&index, &value) in dae.flux_index.iter().zip(fluxes) {
        state[index] = value;
    }
    dae.consistent_initial_state(&state, &[]).unwrap()
}

fn charges(c: &[f64; 4], voltage: [f64; 2]) -> [f64; 2] {
    [
        c[0] * voltage[0] + c[1] * voltage[1],
        c[2] * voltage[0] + c[3] * voltage[1],
    ]
}

fn energy(c: &[f64; 4], voltage: [f64; 2]) -> f64 {
    0.5 * (c[0] * voltage[0].powi(2)
        + 2.0 * c[1] * voltage[0] * voltage[1]
        + c[3] * voltage[1].powi(2))
}

#[test]
fn reciprocal_cross_signs_order_and_prepared_remap_preserve_the_complete_energy() {
    for cross in [-0.75, 0.75] {
        let c = [2.0, cross, cross, 3.0];
        let block = CapacitanceBlock::try_new(vec![5, 2], c.to_vec()).unwrap();
        assert_eq!(block.branches(), &[5, 2]);
        assert_eq!(block.matrix_f(), &c);
        let voltage = [1.2, -0.7];
        let expected_charge = charges(&c, voltage);
        let actual_charge = block.charges_c(&voltage).unwrap();
        for (actual, expected) in actual_charge.iter().zip(expected_charge) {
            near(*actual, expected, 2e-15);
        }
        for (actual, expected) in block
            .voltages_v(&actual_charge)
            .unwrap()
            .iter()
            .zip(voltage)
        {
            near(*actual, expected, 2e-15);
        }
        near(
            block.stored_energy_j(&actual_charge).unwrap(),
            energy(&c, voltage),
            3e-15,
        );
        let remapped = block.with_branches(vec![8, 3]).unwrap();
        assert_eq!(remapped.branches(), &[8, 3]);
        assert!(std::ptr::eq(block.matrix_f(), remapped.matrix_f()));
        assert!(std::ptr::eq(
            block.inverse_matrix_per_f(),
            remapped.inverse_matrix_per_f()
        ));
        assert!(block.with_branches(vec![1]).is_err());
        assert!(block.with_branches(vec![1, 1]).is_err());
        assert!(block.charges_c(&[1.0]).is_err());
        assert!(block.voltages_v(&[0.0, f64::NAN]).is_err());
    }
}

#[test]
fn coupled_rc_uses_the_full_matrix_and_closes_actual_resistor_heat() {
    for cross in [-0.6e-6, 0.6e-6] {
        let c = [2e-6, cross, cross, 3e-6];
        let block = CapacitanceBlock::try_new(vec![2, 0], c.to_vec()).unwrap();
        // Deliberately interleave branch types and reverse block order. Scalar
        // placeholder values must not leak into constitutive voltage or work.
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (2, 0, Branch::Capacitor { farads: 9.0 }),
                (1, 0, Branch::Resistor { ohms: 1000.0 }),
                (1, 0, Branch::Capacitor { farads: 7.0 }),
                (2, 0, Branch::Resistor { ohms: 2200.0 }),
            ],
            transformers: vec![],
        };
        let dae = assemble_circuit_with_capacitance_blocks(&graph, &[block], &[], &[]).unwrap();
        assert!(dae.system.storage().exact_quadratic_matrix().is_some());
        let mut voltage = [1.0, -0.5];
        let q = charges(&c, voltage);
        let mut state = initial(&dae, &[q[1], q[0]], &[]);
        let start = energy(&c, voltage);
        near(dae.system.hamiltonian(&state), start, 2e-20);
        let (mut native_heat, mut physical_heat) = (0.0, 0.0);
        let dt = 1e-4;
        for _ in 0..16 {
            // (C + dt*G/2) v1 = (C - dt*G/2) v0, solved by the
            // explicit two-by-two determinant rather than production LA.
            let a = [c[0] + dt / 2000.0, c[1], c[2], c[3] + dt / 4400.0];
            let rhs = [
                (c[0] - dt / 2000.0) * voltage[0] + c[1] * voltage[1],
                c[2] * voltage[0] + (c[3] - dt / 4400.0) * voltage[1],
            ];
            let determinant = a[0] * a[3] - a[1] * a[2];
            let next_voltage = [
                (rhs[0] * a[3] - a[1] * rhs[1]) / determinant,
                (a[0] * rhs[1] - rhs[0] * a[2]) / determinant,
            ];
            let (step, defect) = dae.step_audited(&state, &[], dt).unwrap();
            let next_charge = charges(&c, next_voltage);
            for (&index, expected) in dae
                .charge_index
                .iter()
                .zip([next_charge[1], next_charge[0]])
            {
                near(step.x[index], expected, 2e-17);
            }
            for (actual, expected) in dae
                .capacitor_voltages(&step.x)
                .unwrap()
                .iter()
                .zip([next_voltage[1], next_voltage[0]])
            {
                near(*actual, expected, 2e-11);
            }
            let mean = [
                0.5 * (voltage[0] + next_voltage[0]),
                0.5 * (voltage[1] + next_voltage[1]),
            ];
            let heat = dt * (mean[0].powi(2) / 1000.0 + mean[1].powi(2) / 2200.0);
            near(step.dissipated, heat, 2e-18);
            near(
                step.delta_h,
                energy(&c, next_voltage) - energy(&c, voltage),
                2e-18,
            );
            assert!(defect < 2e-18);
            native_heat += step.dissipated;
            physical_heat += heat;
            state = step.x;
            voltage = next_voltage;
        }
        near(dae.system.hamiltonian(&state) + native_heat, start, 5e-17);
        near(native_heat, physical_heat, 5e-17);
    }
}

#[test]
fn matrix_storage_coexists_with_nonlinear_charge_and_full_mutual_inductance() {
    let graph = CircuitGraph {
        node_count: 6,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 10.0 }),
            (2, 0, Branch::Capacitor { farads: 20.0 }),
            (1, 0, Branch::Resistor { ohms: 2.0 }),
            (2, 0, Branch::Resistor { ohms: 3.0 }),
            (3, 0, Branch::Capacitor { farads: 30.0 }),
            (3, 0, Branch::Resistor { ohms: 1.0 }),
            (4, 0, Branch::Inductor { henries: 2.0 }),
            (5, 0, Branch::Inductor { henries: 3.0 }),
            (4, 0, Branch::Resistor { ohms: 2.0 }),
            (5, 0, Branch::Resistor { ohms: 3.0 }),
        ],
        transformers: vec![],
    };
    let c = [2.0, 0.5, 0.5, 3.0];
    let block = CapacitanceBlock::try_new(vec![0, 1], c.to_vec()).unwrap();
    let curve = CapacitorChargeCurve::try_new(vec![
        (-3.0, -5.0),
        (-1.0, -1.0),
        (0.0, 0.0),
        (1.0, 1.0),
        (3.0, 5.0),
    ])
    .unwrap();
    let bindings = [CapacitorCurve { branch: 4, curve }];
    let mutual = [MutualInductance {
        branches: (6, 7),
        henries: 0.5,
    }];
    let dae =
        assemble_circuit_with_capacitance_blocks(&graph, &[block], &bindings, &mutual).unwrap();
    let q = charges(&c, [1.0, -0.5]);
    let mut state = initial(&dae, &[q[0], q[1], 1.1], &[0.6, -0.3]);
    let independent = |state: &[f64]| {
        let q = [state[dae.charge_index[0]], state[dae.charge_index[1]]];
        let c_energy = 0.5 * (3.0 * q[0].powi(2) - q[0] * q[1] + 2.0 * q[1].powi(2)) / 5.75;
        let charge = state[dae.charge_index[2]].abs();
        let nonlinear = if charge <= 1.0 {
            0.5 * charge.powi(2)
        } else {
            charge.powi(2) - charge + 0.5
        };
        let flux = [state[dae.flux_index[0]], state[dae.flux_index[1]]];
        let magnetic =
            0.5 * (3.0 * flux[0].powi(2) - flux[0] * flux[1] + 2.0 * flux[1].powi(2)) / 5.75;
        c_energy + nonlinear + magnetic
    };
    let start = independent(&state);
    near(dae.system.hamiltonian(&state), start, 2e-15);
    let mut heat = 0.0;
    for _ in 0..5 {
        let (step, defect) = dae.step_audited(&state, &[], 0.2).unwrap();
        near(
            step.delta_h,
            independent(&step.x) - independent(&state),
            2e-11,
        );
        near(dae.system.hamiltonian(&step.x), independent(&step.x), 2e-14);
        let current = dae.inductor_currents(&step.x).unwrap();
        let flux = [step.x[dae.flux_index[0]], step.x[dae.flux_index[1]]];
        near(current[0], (3.0 * flux[0] - 0.5 * flux[1]) / 5.75, 2e-14);
        near(current[1], (2.0 * flux[1] - 0.5 * flux[0]) / 5.75, 2e-14);
        assert!(defect < 2e-11);
        heat += step.dissipated;
        state = step.x;
    }
    assert!(state[dae.charge_index[2]] < 1.0);
    near(independent(&state) + heat, start, 8e-11);
}

#[test]
fn stable_block_work_retains_a_small_change_beside_unchanged_bias() {
    let block = CapacitanceBlock::try_new(vec![0, 1], vec![1.0, 0.0, 0.0, 1.0]).unwrap();
    let before = [1e8, 1e-6];
    let after = [1e8, 1.1e-6];
    assert_eq!(
        block.stored_energy_j(&before).unwrap(),
        block.stored_energy_j(&after).unwrap()
    );
    let expected = 0.5 * (after[1] - before[1]) * (after[1] + before[1]);
    near(
        block.energy_change_j(&before, &after).unwrap(),
        expected,
        3e-29,
    );
    near(
        block.energy_change_j(&after, &before).unwrap(),
        -expected,
        3e-29,
    );
    assert_eq!(block.energy_change_j(&before, &before).unwrap(), 0.0);
}

#[test]
fn nonreciprocal_indefinite_gauge_and_unresolved_full_blocks_refuse() {
    for matrix in [
        vec![1.0, 0.2, 0.3, 1.0],
        vec![1.0, 1.0, 1.0, 1.0],
        vec![1.0, -1.0, -1.0, 1.0],
        vec![1.0, 2.0, 2.0, 1.0],
        vec![1.0, 1.0 - f64::EPSILON, 1.0 - f64::EPSILON, 1.0],
        vec![1.0, f64::NAN, f64::NAN, 1.0],
        vec![1.0, 0.0, 0.0, -1.0],
    ] {
        assert!(CapacitanceBlock::try_new(vec![0, 1], matrix).is_err());
    }
    // Every two-by-two principal minor is positive, but the full matrix has
    // a negative common mode. Pairwise coupling bounds are insufficient.
    assert!(
        CapacitanceBlock::try_new(
            vec![0, 1, 2],
            vec![1.0, -0.75, -0.75, -0.75, 1.0, -0.75, -0.75, -0.75, 1.0,]
        )
        .is_err()
    );
    assert!(CapacitanceBlock::try_new(vec![0, 0], vec![1.0; 4]).is_err());
    assert!(CapacitanceBlock::try_new(vec![], vec![]).is_err());
    assert!(CapacitanceBlock::try_new((0..129).collect(), vec![]).is_err());
}

#[test]
fn duplicate_storage_owners_and_wrong_branch_types_refuse_without_changing_graph() {
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Capacitor { farads: 2.0 }),
            (1, 0, Branch::Resistor { ohms: 3.0 }),
        ],
        transformers: vec![],
    };
    let before = graph.clone();
    let block = CapacitanceBlock::try_new(vec![0, 1], vec![2.0, 0.5, 0.5, 3.0]).unwrap();
    let curve = CapacitorChargeCurve::try_new(vec![(-1.0, -1.0), (0.0, 0.0), (1.0, 1.0)]).unwrap();
    assert!(
        assemble_circuit_with_capacitance_blocks(&graph, &[block.clone(), block.clone()], &[], &[])
            .is_err()
    );
    assert!(
        assemble_circuit_with_capacitance_blocks(
            &graph,
            &[block.clone()],
            &[CapacitorCurve { branch: 1, curve }],
            &[]
        )
        .is_err()
    );
    for branches in [vec![0, 2], vec![0, 4]] {
        assert!(
            assemble_circuit_with_capacitance_blocks(
                &graph,
                &[block.with_branches(branches).unwrap()],
                &[],
                &[]
            )
            .is_err()
        );
    }
    assert_eq!(graph, before);
    let ordinary = assemble_circuit_with_capacitor_curves(&graph, &[], &[]).unwrap();
    let empty = assemble_circuit_with_capacitance_blocks(&graph, &[], &[], &[]).unwrap();
    assert_eq!(
        ordinary.system.storage().exact_quadratic_matrix(),
        empty.system.storage().exact_quadratic_matrix()
    );
    assert_eq!(ordinary.charge_index, empty.charge_index);
}
