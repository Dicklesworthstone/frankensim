//! Full reciprocal field loads share the original finite C/L circuit step.
use fs_phs::circuit::conductance::{
    ConductanceBlock, MAX_CONDUCTANCE_BLOCK_COORDINATES, MAX_CONDUCTANCE_BLOCKS,
};
use fs_phs::circuit::{
    Branch, CircuitError, CircuitGraph, MutualInductance, assemble_circuit_with_conductance,
    assemble_circuit_with_conductance_polled, assemble_circuit_with_couplings,
    assemble_circuit_with_couplings_polled,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn two_caps() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    }
}

#[test]
fn reciprocal_mutual_currents_keep_sign_reference_and_independent_power() {
    let block = ConductanceBlock::try_new(vec![7, 4, 9], vec![2.0, 1.0, 1.0, 3.0]).unwrap();
    assert_eq!(block.nodes(), &[7, 4, 9]);
    assert_eq!(block.matrix_s(), &[2.0, 1.0, 1.0, 3.0]);
    for voltage in [[10.0, 12.0, 9.0], [-90.0, -88.0, -91.0]] {
        let current = block.terminal_currents_a(&voltage).unwrap();
        assert_eq!(current, vec![-2.0, 3.0, -1.0]);
        close(block.power_w(&voltage).unwrap(), 7.0);
        close(voltage.iter().zip(&current).map(|(v, i)| v * i).sum(), 7.0);
    }
    // Rank-one passive conductance has a true zero-loss mode. Either mutual
    // sign is valid; no positive pair-resistor fit represents this positive one.
    let rank_one = ConductanceBlock::try_new(vec![0, 1, 2], vec![1.0; 4]).unwrap();
    assert_eq!(
        rank_one.terminal_currents_a(&[0.0, 2.0, -2.0]).unwrap(),
        vec![0.0; 3]
    );
    close(rank_one.power_w(&[0.0, 2.0, -2.0]).unwrap(), 0.0);
    let scaled = ConductanceBlock::try_new(vec![0, 1, 2], vec![1e-20, 0.0, 0.0, 1e20]).unwrap();
    close(scaled.power_w(&[0.0, 1e10, 1e-10]).unwrap(), 2.0);
}

#[test]
fn full_mutual_load_changes_both_native_capacitor_modes() {
    let block = ConductanceBlock::try_new(vec![0, 1, 2], vec![2.0, 1.0, 1.0, 2.0]).unwrap();
    let dae = assemble_circuit_with_conductance(&two_caps(), std::slice::from_ref(&block)).unwrap();
    let mut initial = vec![0.0; dae.system.state_dim()];
    initial[dae.charge_index[0]] = 1.0;
    let initial = dae.consistent_initial_state(&initial, &[]).unwrap();
    let dt = 0.25;
    let (step, defect) = dae.step_audited(&initial, &[], dt).unwrap();
    // Exact implicit-midpoint amplitudes for the original rates 3 and 1.
    let fast = 5.0 / 11.0;
    let slow = 7.0 / 9.0;
    close(step.x[dae.charge_index[0]], (fast + slow) / 2.0);
    close(step.x[dae.charge_index[1]], (fast - slow) / 2.0);
    let midpoint = [
        0.0,
        f64::midpoint(
            initial[dae.node_potential_index[0]],
            step.x[dae.node_potential_index[0]],
        ),
        f64::midpoint(
            initial[dae.node_potential_index[1]],
            step.x[dae.node_potential_index[1]],
        ),
    ];
    let analytic_power =
        2.0 * midpoint[1].powi(2) + 2.0 * midpoint[1] * midpoint[2] + 2.0 * midpoint[2].powi(2);
    close(block.power_w(&midpoint).unwrap(), analytic_power);
    close(step.dissipated, dt * analytic_power);
    close(step.delta_h + dt * analytic_power, 0.0);
    close(defect, 0.0);
}

#[test]
fn finite_inductor_and_both_capacitors_use_one_original_kirchhoff_step() {
    let block = ConductanceBlock::try_new(vec![0, 1, 2], vec![2.0, 1.0, 1.0, 2.0]).unwrap();
    let mut graph = two_caps();
    graph
        .branches
        .push((1, 0, Branch::Inductor { henries: 1.0 }));
    let dae = assemble_circuit_with_conductance(&graph, std::slice::from_ref(&block)).unwrap();
    let mut initial = vec![0.0; dae.system.state_dim()];
    initial[dae.flux_index[0]] = 0.4;
    initial[dae.charge_index[0]] = 1.0;
    initial[dae.charge_index[1]] = -0.3;
    let initial = dae.consistent_initial_state(&initial, &[]).unwrap();
    let dt: f64 = 0.2;
    // Independently eliminate the midpoint inductor current from KCL:
    // (1+dt+dt²/4)v1 + (dt/2)v2 = q1 - (dt/2)i0;
    // (dt/2)v1 + (1+dt)v2 = q2.
    let a = 1.0 + dt + dt * dt / 4.0;
    let b = dt / 2.0;
    let d = 1.0 + dt;
    let rhs = 1.0 - dt * 0.4 / 2.0;
    let determinant = a * d - b * b;
    let v1 = (rhs * d + 0.3 * b) / determinant;
    let v2 = (-0.3 * a - b * rhs) / determinant;
    let (step, defect) = dae.step_audited(&initial, &[], dt).unwrap();
    close(step.x[dae.flux_index[0]], 0.4 + dt * v1);
    close(step.x[dae.charge_index[0]], 2.0 * v1 - 1.0);
    close(step.x[dae.charge_index[1]], 2.0 * v2 + 0.3);
    let heat = dt * (2.0 * v1 * v1 + 2.0 * v1 * v2 + 2.0 * v2 * v2);
    close(step.dissipated, heat);
    close(step.delta_h, -heat);
    close(block.power_w(&[0.0, v1, v2]).unwrap() * dt, heat);
    close(defect, 0.0);
}

#[test]
fn coupled_windings_drive_the_full_field_load_with_one_power_balance() {
    let block = ConductanceBlock::try_new(vec![0, 1, 2], vec![2.0, 1.0, 1.0, 3.0]).unwrap();
    for orientation in [1.0, -1.0] {
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 2.0 }),
                if orientation > 0.0 {
                    (2, 0, Branch::Inductor { henries: 1.0 })
                } else {
                    (0, 2, Branch::Inductor { henries: 1.0 })
                },
            ],
            transformers: vec![],
        };
        let links = [MutualInductance {
            branches: (1, 2),
            henries: orientation * 0.5,
        }];
        let dae =
            assemble_circuit_with_couplings(&graph, &links, std::slice::from_ref(&block)).unwrap();
        let mut state = dae
            .consistent_initial_state(&vec![0.0; dae.system.state_dim()], &[1.0])
            .unwrap();
        let (mut flux1, mut i2, mut supplied, mut heat) = (0.0, 0.0, 0.0, 0.0);
        let dt = 0.2;
        for _ in 0..5 {
            // Independent midpoint elimination: V1=1, L=[[2,.5],[.5,1]],
            // i2+V1+3V2=0. The secondary leakage L2-M²/L1 is 7/8 H.
            let leakage = 0.875;
            let h = dt / 2.0;
            let mid_i2 = (i2 - h * (1.0 / 3.0 + 0.25) / leakage) / (1.0 + h / (3.0 * leakage));
            let v2 = -(mid_i2 + 1.0) / 3.0;
            let mid_i1 = (flux1 + h - 0.5 * mid_i2) / 2.0;
            let expected_heat = dt * (2.0 + 2.0 * v2 + 3.0 * v2 * v2);
            let expected_supply = dt * (mid_i1 + 2.0 + v2);
            let (step, defect) = dae.step_audited(&state, &[1.0], dt).unwrap();
            flux1 += dt;
            i2 = 2.0 * mid_i2 - i2;
            let i1 = (flux1 - 0.5 * i2) / 2.0;
            let currents = dae.inductor_currents(&step.x).unwrap();
            close(currents[0], i1);
            close(currents[1], orientation * i2);
            close(step.dissipated, expected_heat);
            close(step.supplied, expected_supply);
            close(block.power_w(&[0.0, 1.0, v2]).unwrap() * dt, expected_heat);
            close(defect, 0.0);
            supplied += step.supplied;
            heat += step.dissipated;
            // Full independently reconstructed magnetic energy, including M.
            let stored = 0.5 * (2.0 * i1 * i1 + i1 * i2 + i2 * i2);
            close(dae.system.hamiltonian(&step.x), stored);
            close(supplied, heat + stored);
            state = step.x;
        }
        let mut refused = links;
        refused[0].henries = 2.0;
        assert!(
            assemble_circuit_with_couplings(&graph, &refused, std::slice::from_ref(&block))
                .is_err()
        );
        assert!(
            assemble_circuit_with_couplings(&graph, &[links[0]; 2], std::slice::from_ref(&block))
                .is_err()
        );
        let mut calls = 0;
        assert!(
            assemble_circuit_with_couplings_polled(
                &graph,
                &links,
                std::slice::from_ref(&block),
                || {
                    calls += 1;
                    if calls == 3 {
                        Err(CircuitError::Invalid {
                            what: "cancelled magnetic preparation",
                        })
                    } else {
                        Ok(())
                    }
                },
            )
            .is_err()
        );
        assert_eq!(calls, 3);
    }
}

#[test]
fn electrode_with_no_lumped_branch_floats_at_zero_net_current() {
    let block = ConductanceBlock::try_new(vec![0, 1, 2], vec![2.0, -1.0, -1.0, 2.0]).unwrap();
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 })],
        transformers: vec![],
    };
    let dae = assemble_circuit_with_conductance(&graph, std::slice::from_ref(&block)).unwrap();
    let mut initial = vec![0.0; dae.system.state_dim()];
    initial[dae.charge_index[0]] = 1.0;
    let initial = dae.consistent_initial_state(&initial, &[]).unwrap();
    close(initial[dae.node_potential_index[1]], 0.5);
    let (step, _) = dae.step_audited(&initial, &[], 0.4).unwrap();
    close(step.x[dae.charge_index[0]], 7.0 / 13.0);
    let voltage = [0.0, 10.0 / 13.0, 5.0 / 13.0];
    let currents = block.terminal_currents_a(&voltage).unwrap();
    close(currents[2], 0.0);
    close(currents[1], 15.0 / 13.0);
    close(step.dissipated, 0.4 * 1.5 * voltage[1].powi(2));

    // A zero block is not a fictitious electrical anchor.
    let zero = ConductanceBlock::try_new(vec![0, 1, 2], vec![0.0; 4]).unwrap();
    assert!(matches!(
        assemble_circuit_with_conductance(&graph, &[zero]),
        Err(CircuitError::FloatingNode { node: 2 })
    ));
}

#[test]
fn arbitrary_reference_and_mapping_stamp_complete_gauge_lift() {
    let block = ConductanceBlock::try_new(vec![2, 1, 0], vec![2.0, 1.0, 1.0, 3.0]).unwrap();
    let dae = assemble_circuit_with_conductance(&two_caps(), &[block]).unwrap();
    let n = dae.system.state_dim();
    let (_, resistance, _) = dae.system.structure();
    let v = &dae.node_potential_index;
    assert_eq!(resistance[v[0] * n + v[0]], 2.0);
    assert_eq!(resistance[v[0] * n + v[1]], -3.0);
    assert_eq!(resistance[v[1] * n + v[0]], -3.0);
    assert_eq!(resistance[v[1] * n + v[1]], 7.0);

    // If the reference current is identically zero, it is not incident.
    let graph = CircuitGraph {
        node_count: 4,
        branches: two_caps().branches,
        transformers: vec![],
    };
    let absent_reference =
        ConductanceBlock::try_new(vec![3, 1, 2], vec![1.0, -1.0, -1.0, 1.0]).unwrap();
    assert!(matches!(
        assemble_circuit_with_conductance(&graph, &[absent_reference]),
        Err(CircuitError::FloatingNode { node: 3 })
    ));
}

#[test]
fn malformed_active_overbudget_and_cancelled_laws_refuse() {
    for (nodes, matrix) in [
        (vec![0, 0], vec![1.0]),
        (vec![0, 1], vec![-1.0]),
        (vec![0, 1], vec![f64::NAN]),
        (vec![0, 1, 2], vec![1.0, 2.0, 2.0, 1.0]),
        (vec![0, 1, 2], vec![2.0, 1.0, 0.0, 2.0]),
    ] {
        assert!(ConductanceBlock::try_new(nodes, matrix).is_err());
    }
    let n = MAX_CONDUCTANCE_BLOCK_COORDINATES + 1;
    assert!(ConductanceBlock::try_new((0..=n).collect(), vec![0.0; n * n]).is_err());
    let refusal = || {
        Err(CircuitError::Invalid {
            what: "cancelled test",
        })
    };
    assert!(ConductanceBlock::try_new_polled(vec![0, 1], vec![1.0], refusal).is_err());
    let block = ConductanceBlock::try_new(vec![0, 1, 2], vec![2.0, 1.0, 1.0, 2.0]).unwrap();
    assert!(block.terminal_currents_a(&[0.0, 1.0]).is_err());
    assert!(block.power_w(&[0.0, f64::INFINITY, 1.0]).is_err());
    assert!(
        assemble_circuit_with_conductance_polled(
            &two_caps(),
            std::slice::from_ref(&block),
            refusal
        )
        .is_err()
    );
    assert!(
        assemble_circuit_with_conductance(&two_caps(), &vec![block; MAX_CONDUCTANCE_BLOCKS + 1])
            .is_err()
    );

    let mut polls = 0;
    assert!(
        ConductanceBlock::try_new_polled(vec![0, 1], vec![1.0], || {
            polls += 1;
            if polls == 2 { refusal() } else { Ok(()) }
        })
        .is_err()
    );
    assert_eq!(polls, 2);
    let source_overflow = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::CurrentSource { port: usize::MAX })],
        transformers: vec![],
    };
    assert!(assemble_circuit_with_conductance(&source_overflow, &[]).is_err());
    let state_overflow = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 }); 128],
        transformers: vec![],
    };
    assert!(assemble_circuit_with_conductance(&state_overflow, &[]).is_err());
    let out_of_range = ConductanceBlock::try_new(vec![0, 3], vec![1.0]).unwrap();
    assert!(assemble_circuit_with_conductance(&two_caps(), &[out_of_range]).is_err());
}
