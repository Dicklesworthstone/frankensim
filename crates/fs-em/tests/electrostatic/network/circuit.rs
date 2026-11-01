//! Full original field matrices reach native charge history without pair fits.
use super::*;
use fs_em::electrostatic::network::floating::FloatingCapacitanceControls;
use fs_phs::circuit::capacitance::assemble_circuit_with_capacitance_blocks;

fn non_pair_problem(scale: f64) -> PlanarElectrostatics {
    hexagon(
        [0.6, 0.0],
        [0.5, 8.0, 0.1, 0.2, 80.0, 1.0].map(|x| x * scale),
    )
}

#[test]
fn positive_mutual_field_matrix_drives_native_rc_at_microfarad_and_picofarad_scales() {
    for scale in [1e-6, 1e-11] {
        let problem = non_pair_problem(scale);
        let matrix = problem
            .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
            .unwrap();
        assert!(matrix.maxwell_matrix_f()[5] > 0.0);
        let volts = [0.0, 1.0, -1.0];
        assert!(
            matrix
                .circuit_at_voltages(&[0, 1, 2], &volts, 1e-15, 1e-15, &mut || false)
                .is_err()
        );
        let bridge = matrix
            .coupled_circuit_at_voltages(0, &[0, 1, 2], &volts, 1e-15, 1e-15, &mut || false)
            .unwrap();
        assert!(std::ptr::eq(bridge.source(), &matrix));
        assert_eq!(bridge.coordinate_terminals(), &[1, 2]);
        assert_eq!(bridge.branches().len(), 2);
        assert_eq!(
            bridge.capacitance_block().matrix_f(),
            matrix.grounded_matrix_f(0).unwrap()
        );
        // Place a resistor BEFORE the source-derived branches to exercise actual
        // graph-index remapping, without refactoring or copying C and its inverse.
        let r1 = 1e-3 / scale;
        let r2 = 2.2e-3 / scale;
        let mut branches = vec![(1, 0, Branch::Resistor { ohms: r1 })];
        branches.extend_from_slice(bridge.branches());
        branches.push((2, 0, Branch::Resistor { ohms: r2 }));
        let block = bridge.block_at_branch_offset(1).unwrap();
        assert_eq!(block.branches(), &[1, 2]);
        assert!(std::ptr::eq(
            block.matrix_f().as_ptr(),
            bridge.capacitance_block().matrix_f().as_ptr()
        ));
        let dae = assemble_circuit_with_capacitance_blocks(
            &CircuitGraph {
                node_count: 3,
                branches,
                transformers: vec![],
            },
            &[block],
            &[],
            &[],
        )
        .unwrap();
        let mut initial = vec![0.0; dae.system.state_dim()];
        for (&slot, &charge) in dae.charge_index.iter().zip(bridge.initial_charges_c()) {
            initial[slot] = charge;
        }
        let initial = dae.consistent_initial_state(&initial, &[]).unwrap();
        near(
            dae.system.hamiltonian(&initial),
            bridge.field_response().energy_j,
        );
        let dt = 1e-4;
        let (next, defect) = dae.step_audited(&initial, &[], dt).unwrap();
        let c = matrix.grounded_matrix_f(0).unwrap();
        // Independent midpoint terminal equations (C/dt + G/2)V1=(C/dt-G/2)V0.
        let a = [
            c[0] / dt + 0.5 / r1,
            c[1] / dt,
            c[2] / dt,
            c[3] / dt + 0.5 / r2,
        ];
        let rhs = [
            (c[0] / dt - 0.5 / r1) * volts[1] + c[1] / dt * volts[2],
            c[2] / dt * volts[1] + (c[3] / dt - 0.5 / r2) * volts[2],
        ];
        let determinant = a[0] * a[3] - a[1] * a[2];
        let after = [
            0.0,
            (rhs[0] * a[3] - a[1] * rhs[1]) / determinant,
            (a[0] * rhs[1] - rhs[0] * a[2]) / determinant,
        ];
        let expected_charge = [
            c[0] * after[1] + c[1] * after[2],
            c[2] * after[1] + c[3] * after[2],
        ];
        for (slot, expected) in dae.charge_index.iter().zip(expected_charge) {
            near(next.x[*slot], expected);
        }
        let recovered = matrix
            .field_at_voltages(&after, 1e-15, 1e-15, &mut || false)
            .unwrap();
        for (&terminal, &slot) in bridge.coordinate_terminals().iter().zip(&dae.charge_index) {
            let weak: f64 = TERMINALS[terminal]
                .iter()
                .map(|&vertex| recovered.nodal_boundary_charge_c()[vertex].unwrap())
                .sum();
            near(next.x[slot], weak);
        }
        near(dae.system.hamiltonian(&next.x), recovered.energy_j());
        let heat = dt * (0.5 * (volts[1] + after[1])).powi(2) / r1
            + dt * (0.5 * (volts[2] + after[2])).powi(2) / r2;
        near(next.dissipated, heat);
        near(
            recovered.energy_j() + heat,
            bridge.field_response().energy_j,
        );
        assert!(defect < 1e-15);
        assert!(bridge.maximum_charge_residual_c() < 1e-15);
        assert!(bridge.energy_residual_j().abs() < 1e-15);
        assert_eq!(bridge.source().problem().materials(), problem.materials());
    }
}

#[test]
fn full_block_respects_reference_permutation_common_voltage_and_physical_charge_order() {
    let problem = non_pair_problem(1e-11);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    let volts = [7.0, 10.0, -1.0];
    let original = matrix.response(&volts, &mut || false).unwrap();
    for reference in 0..3 {
        let nodes: Vec<_> = (0..3)
            .map(|t| {
                if t == reference {
                    0
                } else if t < reference {
                    t + 1
                } else {
                    t
                }
            })
            .collect();
        let bridge = matrix
            .coupled_circuit_at_voltages(reference, &nodes, &volts, 1e-15, 1e-15, &mut || false)
            .unwrap();
        assert_eq!(bridge.reference_terminal(), reference);
        let order: Vec<_> = (0..3).filter(|&t| t != reference).collect();
        assert_eq!(bridge.coordinate_terminals(), order);
        for (coordinate, &terminal) in order.iter().enumerate() {
            assert_eq!(bridge.branches()[coordinate].0, nodes[terminal]);
            assert_eq!(bridge.branches()[coordinate].1, 0);
            near(
                bridge.initial_charges_c()[coordinate],
                original.terminal_charge_c[terminal],
            );
        }
        near(
            bridge
                .capacitance_block()
                .stored_energy_j(bridge.initial_charges_c())
                .unwrap(),
            original.energy_j,
        );
        let shifted = matrix
            .coupled_circuit_at_voltages(
                reference,
                &nodes,
                &volts.map(|v| v + 1e12),
                1e-15,
                1e-15,
                &mut || false,
            )
            .unwrap();
        assert_eq!(shifted.initial_charges_c(), bridge.initial_charges_c());
        assert_eq!(
            shifted.field_response().energy_j,
            bridge.field_response().energy_j
        );
    }
    let permutation = [2, 0, 1];
    let permuted = problem
        .capacitance_matrix(
            &permutation.map(|t| TERMINALS[t]),
            extraction(),
            &mut || false,
        )
        .unwrap();
    let bridge = permuted
        .coupled_circuit_at_voltages(
            1,
            &[2, 0, 1],
            &permutation.map(|t| volts[t]),
            1e-15,
            1e-15,
            &mut || false,
        )
        .unwrap();
    for (&terminal, &charge) in bridge
        .coordinate_terminals()
        .iter()
        .zip(bridge.initial_charges_c())
    {
        near(charge, original.terminal_charge_c[permutation[terminal]]);
    }
    near(
        bridge
            .capacitance_block()
            .stored_energy_j(bridge.initial_charges_c())
            .unwrap(),
        original.energy_j,
    );
}

#[test]
fn charged_floating_electrode_remains_a_native_conserved_coordinate_and_energy_reservoir() {
    let problem = non_pair_problem(1e-6);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    let c = matrix.grounded_matrix_f(0).unwrap();
    let trapped = 2.0 * c[3];
    let floating = matrix
        .prepare_floating(
            &[2],
            FloatingCapacitanceControls {
                maximum_floating_terminals: 3,
                maximum_dense_work: 4096,
                coefficient_tolerance_f: 1e-15,
            },
            &mut || false,
        )
        .unwrap();
    let actual = floating
        .solve(&[0.0, 1.0], &[trapped], 1e-15, 1e-15, &mut || false)
        .unwrap();
    let bridge = matrix
        .coupled_circuit_at_voltages(
            0,
            &[0, 1, 2],
            &actual.terminal_voltage_v,
            1e-15,
            1e-15,
            &mut || false,
        )
        .unwrap();
    assert_eq!(bridge.coordinate_terminals(), &[1, 2]);
    near(bridge.initial_charges_c()[1], trapped);
    let resistance = 1000.0;
    let mut branches = bridge.branches().to_vec();
    branches.push((1, 0, Branch::Resistor { ohms: resistance }));
    let dae = assemble_circuit_with_capacitance_blocks(
        &CircuitGraph {
            node_count: 3,
            branches,
            transformers: vec![],
        },
        &[bridge.capacitance_block().clone()],
        &[],
        &[],
    )
    .unwrap();
    let mut state = vec![0.0; dae.system.state_dim()];
    for (&slot, &charge) in dae.charge_index.iter().zip(bridge.initial_charges_c()) {
        state[slot] = charge;
    }
    state = dae.consistent_initial_state(&state, &[]).unwrap();
    let effective = c[0] - c[1] * c[2] / c[3];
    let dt = 1e-4;
    let voltage = (2.0 * resistance * effective - dt) / (2.0 * resistance * effective + dt);
    let independent = floating
        .solve(&[0.0, voltage], &[trapped], 1e-15, 1e-15, &mut || false)
        .unwrap();
    let (next, defect) = dae.step_audited(&state, &[], dt).unwrap();
    near(next.x[dae.charge_index[1]], trapped);
    near(
        next.x[dae.charge_index[0]],
        independent.terminal_charge_c[1],
    );
    near(
        dae.system.hamiltonian(&next.x),
        independent.field().energy_j(),
    );
    // Trapped charge contributes a nonzero energy floor; it is not Schur-eliminated.
    let expected_energy = 0.5 * effective * voltage * voltage + 0.5 * trapped * trapped / c[3];
    near(dae.system.hamiltonian(&next.x), expected_energy);
    let heat = dt * (0.5 * (1.0 + voltage)).powi(2) / resistance;
    near(next.dissipated, heat);
    near(expected_energy + heat, bridge.field_response().energy_j);
    assert!(defect < 1e-15);
}

#[test]
fn unreferenced_components_node_aliases_work_mismatch_and_cancellation_refuse_atomically() {
    let problem = non_pair_problem(1e-6);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    for (reference, nodes) in [(3, [0, 1, 2]), (0, [0, 0, 2]), (0, [0, usize::MAX, 2])] {
        assert!(
            matrix
                .coupled_circuit_at_voltages(
                    reference,
                    &nodes,
                    &[0.0, 1.0, -1.0],
                    1e-15,
                    1e-15,
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(
        matrix
            .coupled_circuit_at_voltages(
                0,
                &[0, 1, 2],
                &[0.0, 1e20, -1e20],
                1e-30,
                1e-30,
                &mut || false
            )
            .is_err()
    );
    let mut visits = 0;
    let accepted = matrix
        .coupled_circuit_at_voltages(0, &[0, 1, 2], &[0.0, 1.0, -1.0], 1e-15, 1e-15, &mut || {
            visits += 1;
            false
        })
        .unwrap();
    assert!(accepted.block_at_branch_offset(usize::MAX).is_err());
    let mut calls = 0;
    assert!(matches!(
        matrix.coupled_circuit_at_voltages(
            0,
            &[0, 1, 2],
            &[0.0, 1.0, -1.0],
            1e-15,
            1e-15,
            &mut || {
                calls += 1;
                calls == visits
            }
        ),
        Err(ElectrostaticError::Cancelled)
    ));
    let replay = matrix
        .coupled_circuit_at_voltages(0, &[0, 1, 2], &[0.0, 1.0, -1.0], 1e-15, 1e-15, &mut || {
            false
        })
        .unwrap();
    assert_eq!(accepted.initial_charges_c(), replay.initial_charges_c());

    let mesh = TriComplex2::from_indexed_triangles(
        "two independent dielectric components",
        vec![
            [0., 0.],
            [1., 0.],
            [1., 1.],
            [0., 1.],
            [2., 0.],
            [3., 0.],
            [3., 1.],
            [2., 1.],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    let problem =
        PlanarElectrostatics::new(mesh, vec![dielectric(4e-6); 4], limits(), &mut || false)
            .unwrap();
    let disjoint = problem
        .capacitance_matrix(
            &[&[0, 3], &[1, 2], &[4, 7], &[5, 6]],
            extraction(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        disjoint.coupled_circuit_at_voltages(
            0,
            &[0, 1, 2, 3],
            &[0., 1., 0., 2.],
            1e-15,
            1e-15,
            &mut || false
        ),
        Err(ElectrostaticError::Invalid(
            "every coupled dielectric component needs the explicit circuit reference"
        ))
    ));
    // An electrode physically joining pieces supplies a real reference route.
    let joined = problem
        .capacitance_matrix(
            &[&[0, 3], &[1, 2, 4, 7], &[5, 6]],
            extraction(),
            &mut || false,
        )
        .unwrap();
    assert!(
        joined
            .coupled_circuit_at_voltages(0, &[0, 1, 2], &[0., 1., 2.], 1e-15, 1e-15, &mut || false)
            .is_ok()
    );
}
