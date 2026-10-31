//! Independent three-conductor equations and native floating-node RC dynamics.
use super::*;
use fs_em::electrostatic::network::floating::FloatingCapacitanceControls;

fn reduction_controls() -> FloatingCapacitanceControls {
    FloatingCapacitanceControls {
        maximum_floating_terminals: 4,
        maximum_dense_work: 4096,
        coefficient_tolerance_f: 1e-15,
    }
}

#[test]
fn neutral_and_charged_floating_electrodes_keep_native_energy_and_total_charge() {
    for epsilon in [4e-11, 4e-6] {
        let problem = hexagon([0.0; 2], [epsilon; 6]);
        let matrix = problem
            .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
            .unwrap();
        let floating = matrix
            .prepare_floating(&[2], reduction_controls(), &mut || false)
            .unwrap();
        let pair = 7.0 * epsilon * DEPTH / (6.0 * 3.0_f64.sqrt());
        assert!(std::ptr::eq(floating.source(), &matrix));
        assert_eq!(floating.retained_terminals(), &[0, 1]);
        assert_eq!(floating.floating_terminals(), &[2]);
        // The floating route adds two pair capacitances in series to the
        // existing pair. Merely deleting its row would give the WRONG 2*C.
        for (&value, sign) in floating.reduced_matrix_f().iter().zip([1., -1., -1., 1.]) {
            near(value, sign * 1.5 * pair);
        }
        for charge in [0.0, 2.0 * pair, -3.0 * pair] {
            let state = floating
                .solve(&[0., 6.], &[charge], 1e-15, 1e-15, &mut || false)
                .unwrap();
            let expected_float = 3.0 + charge / (2.0 * pair);
            assert_eq!(&state.terminal_voltage_v[..2], &[0., 6.]);
            near(state.terminal_voltage_v[2], expected_float);
            near(state.terminal_charge_c[0], -9.0 * pair - 0.5 * charge);
            near(state.terminal_charge_c[1], 9.0 * pair - 0.5 * charge);
            assert!((state.terminal_charge_c[2] - charge).abs() < 1e-12 * pair);
            near(
                state.field().energy_j(),
                27.0 * pair + charge * charge / (4.0 * pair),
            );
            assert!(state.maximum_floating_charge_residual_c < 1e-12 * pair);
            assert!(state.terminal_energy_residual_j.abs() < 1e-12 * pair);
            assert!(std::ptr::eq(state.field().problem(), &problem));

            let boundary: Vec<_> = TERMINALS
                .iter()
                .zip([0., 6., expected_float])
                .flat_map(|(vertices, value)| vertices.iter().map(move |&v| (v, value)))
                .collect();
            let direct = problem
                .solve(&[0.; 6], &boundary, controls(), &mut || false)
                .unwrap();
            near(state.field().energy_j(), direct.energy_j());
            for (a, b) in state
                .field()
                .electric_field_v_m()
                .iter()
                .zip(direct.electric_field_v_m())
            {
                for axis in 0..2 {
                    assert!((a[axis] - b[axis]).abs() < 1e-12);
                }
            }
            let shifted = floating
                .solve(&[1e12, 1e12 + 6.], &[charge], 1e-15, 1e-15, &mut || false)
                .unwrap();
            assert_eq!(shifted.reference_voltage_v, 1e12);
            assert_eq!(shifted.terminal_voltage_v, state.terminal_voltage_v);
            assert_eq!(shifted.terminal_charge_c, state.terminal_charge_c);
            assert_eq!(
                shifted.field().energy_j().to_bits(),
                state.field().energy_j().to_bits()
            );
        }
        // With all driven voltages zero, trapped charge STILL stores energy.
        let charged = floating
            .solve(&[0.; 2], &[2. * pair], 1e-15, 1e-15, &mut || false)
            .unwrap();
        near(charged.field().energy_j(), pair);
    }
}

#[test]
fn multiple_charge_controls_preserve_order_and_an_explicit_voltage_reference() {
    let epsilon = 4e-6;
    let problem = hexagon([0.0; 2], [epsilon; 6]);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    let pair = 7.0 * epsilon * DEPTH / (6.0 * 3.0_f64.sqrt());
    let floating = matrix
        .prepare_floating(&[2, 0], reduction_controls(), &mut || false)
        .unwrap();
    assert_eq!(floating.retained_terminals(), &[1]);
    assert_eq!(floating.floating_terminals(), &[2, 0]);
    assert!(floating.reduced_matrix_f()[0].abs() < 1e-12 * pair);
    let state = floating
        .solve(&[7.], &[3. * pair, -pair], 1e-15, 1e-15, &mut || false)
        .unwrap();
    assert_eq!(state.reference_terminal, 1);
    assert_eq!(state.reference_voltage_v, 7.);
    near(state.terminal_voltage_v[0], 1.0 / 3.0);
    assert_eq!(state.terminal_voltage_v[1], 0.0);
    near(state.terminal_voltage_v[2], 5.0 / 3.0);
    near(state.terminal_charge_c[0], -pair);
    near(state.terminal_charge_c[1], -2. * pair);
    near(state.terminal_charge_c[2], 3. * pair);
    near(state.field().energy_j(), 7.0 * pair / 3.0);

    let reordered = matrix
        .prepare_floating(&[0, 2], reduction_controls(), &mut || false)
        .unwrap()
        .solve(&[7.], &[-pair, 3. * pair], 1e-15, 1e-15, &mut || false)
        .unwrap();
    for (a, b) in state
        .terminal_voltage_v
        .iter()
        .zip(&reordered.terminal_voltage_v)
    {
        assert!((a - b).abs() < 1e-14);
    }
    let permutation = [2, 0, 1];
    let permuted = problem
        .capacitance_matrix(
            &permutation.map(|i| TERMINALS[i]),
            extraction(),
            &mut || false,
        )
        .unwrap();
    let swapped = permuted
        .prepare_floating(&[0, 1], reduction_controls(), &mut || false)
        .unwrap()
        .solve(&[7.], &[3. * pair, -pair], 1e-15, 1e-15, &mut || false)
        .unwrap();
    for i in 0..3 {
        assert!(
            (swapped.terminal_voltage_v[i] - state.terminal_voltage_v[permutation[i]]).abs()
                < 1e-14
        );
    }
    near(swapped.field().energy_j(), state.field().energy_j());
}

#[test]
fn charged_floating_electrode_survives_native_rc_evolution_and_thermal_loss() {
    let epsilon = 4e-6;
    let problem = hexagon([0.0; 2], [epsilon; 6]);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    let floating = matrix
        .prepare_floating(&[2], reduction_controls(), &mut || false)
        .unwrap();
    let pair = 7.0 * epsilon * DEPTH / (6.0 * 3.0_f64.sqrt());
    let trapped_charge = 1.2 * pair;
    let initial = floating
        .solve(&[0., 3.], &[trapped_charge], 1e-15, 1e-15, &mut || false)
        .unwrap();
    let reduced = matrix
        .circuit_at_voltages(
            &[0, 1, 2],
            &initial.terminal_voltage_v,
            1e-15,
            1e-15,
            &mut || false,
        )
        .unwrap();
    let mut branches = reduced.branches().to_vec();
    let resistance = 1000.0;
    branches.push((1, 0, Branch::Resistor { ohms: resistance }));
    let dae = assemble_circuit(&CircuitGraph {
        node_count: 3,
        branches,
        transformers: vec![],
    })
    .unwrap();
    let mut state = vec![0.; dae.system.state_dim()];
    for (&slot, &charge) in dae.charge_index.iter().zip(reduced.initial_charges_c()) {
        state[slot] = charge;
    }
    state = dae.consistent_initial_state(&state, &[]).unwrap();
    let dt = 1e-4;
    let a = dt / (2.0 * resistance * 1.5 * pair);
    let mut expected_voltage = 3.0;
    let mut heat = 0.0;
    for _ in 0..12 {
        let (next, defect) = dae.step_audited(&state, &[], dt).unwrap();
        expected_voltage *= (1.0 - a) / (1.0 + a);
        heat += next.dissipated;
        let physical = floating
            .solve(
                &[0., expected_voltage],
                &[trapped_charge],
                1e-15,
                1e-15,
                &mut || false,
            )
            .unwrap();
        for (i, &(p, n, branch)) in reduced.branches().iter().enumerate() {
            let Branch::Capacitor { farads } = branch else {
                panic!("capacitor required")
            };
            near(
                next.x[dae.charge_index[i]],
                farads * (physical.terminal_voltage_v[p] - physical.terminal_voltage_v[n]),
            );
        }
        near(dae.system.hamiltonian(&next.x), physical.field().energy_j());
        near(
            physical.field().energy_j() + heat,
            initial.field().energy_j(),
        );
        assert!(defect < 1e-15);
        assert!(next.dissipated > 0.0);
        state = next.x;
    }
}

#[test]
fn positive_mutual_field_matrix_still_supports_charge_constrained_floating_conductors() {
    let problem = hexagon(
        [0.6, 0.0],
        [0.5, 8.0, 0.1, 0.2, 80.0, 1.0].map(|x| x * 1e-6),
    );
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    assert!(matrix.maxwell_matrix_f()[5] > 0.0);
    let state = matrix
        .prepare_floating(&[2], reduction_controls(), &mut || false)
        .unwrap()
        .solve(&[0., 2.], &[1e-7], 1e-14, 1e-14, &mut || false)
        .unwrap();
    let c = matrix.maxwell_matrix_f();
    near(state.terminal_voltage_v[2], (1e-7 - 2.0 * c[7]) / c[8]);
    assert!((state.terminal_charge_c[2] - 1e-7).abs() < 1e-18);
    assert!(state.field().energy_j() > 0.0);
}

#[test]
fn floating_gauges_work_budgets_actual_charge_gates_and_cancellation_refuse() {
    let problem = hexagon([0.0; 2], [4e-6; 6]);
    let matrix = problem
        .capacitance_matrix(&TERMINALS, extraction(), &mut || false)
        .unwrap();
    for indices in [&[][..], &[0, 1, 2], &[1, 1], &[3]] {
        assert!(
            matrix
                .prepare_floating(indices, reduction_controls(), &mut || false)
                .is_err()
        );
    }
    for policy in [
        FloatingCapacitanceControls {
            maximum_floating_terminals: 0,
            ..reduction_controls()
        },
        FloatingCapacitanceControls {
            maximum_dense_work: 1,
            ..reduction_controls()
        },
        FloatingCapacitanceControls {
            coefficient_tolerance_f: f64::NAN,
            ..reduction_controls()
        },
    ] {
        assert!(
            matrix
                .prepare_floating(&[2], policy, &mut || false)
                .is_err()
        );
    }
    let floating = matrix
        .prepare_floating(&[2], reduction_controls(), &mut || false)
        .unwrap();
    assert!(
        floating
            .solve(&[0.], &[0.], 1e-15, 1e-15, &mut || false)
            .is_err()
    );
    assert!(
        floating
            .solve(&[0., 1.], &[f64::NAN], 1e-15, 1e-15, &mut || false)
            .is_err()
    );
    // The actual native field, not just a successful terminal factorization,
    // must pass its charge/energy tolerances at the requested amplitude.
    assert!(
        floating
            .solve(&[0., 1e20], &[1.], 1e-30, 1e-30, &mut || false)
            .is_err()
    );
    let mut visits = 0;
    let expected = floating
        .solve(&[0., 2.], &[1e-7], 1e-15, 1e-15, &mut || {
            visits += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        floating.solve(&[0., 2.], &[1e-7], 1e-15, 1e-15, &mut || {
            count += 1;
            count == visits
        }),
        Err(ElectrostaticError::Cancelled)
    ));
    let replay = floating
        .solve(&[0., 2.], &[1e-7], 1e-15, 1e-15, &mut || false)
        .unwrap();
    assert_eq!(expected.terminal_charge_c, replay.terminal_charge_c);
    assert!(matches!(
        matrix.prepare_floating(&[2], reduction_controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));

    let mesh = TriComplex2::from_indexed_triangles(
        "two disconnected dielectrics",
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
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap();
    let disjoint =
        PlanarElectrostatics::new(mesh, vec![dielectric(4e-6); 4], limits(), &mut || false)
            .unwrap();
    let original = disjoint
        .capacitance_matrix(
            &[&[0, 3], &[1, 2], &[4, 7], &[5, 6]],
            extraction(),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        original.prepare_floating(&[2, 3], reduction_controls(), &mut || false),
        Err(ElectrostaticError::Invalid(
            "every floating dielectric component needs a retained voltage reference"
        ))
    ));
    // A conductor that joins both pieces provides a real capacitive reference
    // route through the first dielectric; topology must include that conductor.
    let joined = disjoint
        .capacitance_matrix(
            &[&[0, 3], &[1, 2, 4, 7], &[5, 6]],
            extraction(),
            &mut || false,
        )
        .unwrap();
    let prepared = joined
        .prepare_floating(&[1, 2], reduction_controls(), &mut || false)
        .unwrap();
    let state = prepared
        .solve(&[0.], &[0., 1e-7], 1e-15, 1e-15, &mut || false)
        .unwrap();
    assert!((state.terminal_charge_c[1]).abs() < 1e-18);
    near(state.terminal_charge_c[2], 1e-7);
}
