//! Floating free-charge constraints against independent capacitor equations.
use super::*;
use fs_em::spatial::electrostatic::network::floating::SpatialFloatingCapacitanceControls;

fn reduction_controls() -> SpatialFloatingCapacitanceControls {
    SpatialFloatingCapacitanceControls {
        maximum_floating_terminals: 4,
        maximum_dense_work: 4096,
        coefficient_tolerance_f: 1e-18,
    }
}

#[test]
fn neutral_and_charged_middle_electrodes_retain_affine_charge_and_trapped_storage() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    let floating = matrix
        .prepare_floating(&[1], reduction_controls(), &mut || false)
        .unwrap();
    assert!(std::ptr::eq(floating.source(), &matrix));
    assert_eq!(floating.retained_terminals(), &[0, 2]);
    assert_eq!(floating.floating_terminals(), &[1]);
    assert_eq!(floating.reference_terminals(), &[0, 0, 0]);
    let total = C1 + C2;
    let series = C1 * C2 / total;
    for (&actual, sign) in floating
        .reduced_matrix_f()
        .iter()
        .zip([1.0, -1.0, -1.0, 1.0])
    {
        near(actual, sign * series);
    }
    for charge in [0.0, 2.0 * total, -0.75 * total] {
        let middle = (charge + 2.0 * C1 - C2) / total;
        let expected_voltage = [0.0, middle - 2.0, -3.0];
        let trapped = 0.5 * charge * charge / total;
        let original = floating
            .solve(&[2.0, -1.0], &[charge], controls(), &mut || false)
            .unwrap();
        for common in [0.0, 2_f64.powi(45)] {
            let state = floating
                .solve(
                    &[common + 2.0, common - 1.0],
                    &[charge],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for (&actual, expected) in state
                .component_relative_voltage_v
                .iter()
                .zip(expected_voltage)
            {
                assert!((actual - expected).abs() < 1e-12);
            }
            assert_eq!(state.component_reference_voltage_v, vec![common + 2.0; 3]);
            for (&actual, expected) in state.terminal_charge_c.iter().zip([
                3.0 * series - C1 * charge / total,
                charge,
                -3.0 * series - C2 * charge / total,
            ]) {
                charge_close(actual, expected);
            }
            for (&actual, expected) in state
                .retained_charge_offset_c
                .iter()
                .zip([-C1 * charge / total, -C2 * charge / total])
            {
                charge_close(actual, expected);
            }
            near(state.trapped_charge_energy_j, trapped);
            near(state.field().energy_j(), 4.5 * series + trapped);
            assert!(state.maximum_floating_charge_residual_c < 1e-18);
            assert!(state.maximum_reduced_charge_residual_c < 1e-18);
            assert!(state.terminal_energy_residual_j.abs() < 1e-18);
            assert!(state.reduced_energy_residual_j.abs() < 1e-18);
            assert!(std::ptr::eq(state.field().problem(), &problem));
            assert!(state.field().report().is_none());
            assert_eq!(
                state.component_relative_voltage_v,
                original.component_relative_voltage_v
            );
            assert_eq!(state.terminal_charge_c, original.terminal_charge_c);
            assert_eq!(
                state.field().energy_j().to_bits(),
                original.field().energy_j().to_bits()
            );
        }
        let direct = direct_field(&problem, &terminals, &[2.0, middle, -1.0]);
        near(direct.energy_j(), original.field().energy_j());
        for (&actual, &expected) in original
            .field()
            .electric_field_v_m()
            .iter()
            .zip(direct.electric_field_v_m())
        {
            near_vector(actual, expected);
        }
        // A neutral floating middle electrode still has nonzero fields and
        // oppositely signed charges on its two physical faces.
        if charge == 0.0 {
            let side_charge: f64 = (8..12)
                .map(|v| original.field().nodal_boundary_charge_c()[v].unwrap())
                .sum();
            near(side_charge, C1 * (middle - 2.0));
            assert!(side_charge.abs() > 1e-12);
            assert!(original.field().energy_j() > 0.0);
        }
    }
    // A reduced matrix with zero retained voltages is not the full charged
    // physical state: the original native field still stores Q_f^2/(2*C_ff).
    let charged = floating
        .solve(&[0.0, 0.0], &[2.0 * total], controls(), &mut || false)
        .unwrap();
    near(charged.component_relative_voltage_v[1], 2.0);
    near(charged.trapped_charge_energy_j, 2.0 * total);
    near(charged.field().energy_j(), 2.0 * total);
}

#[test]
fn coupled_floating_charges_and_input_permutation_retain_the_full_original_energy() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    let q2 = -0.75 * C2;
    let q1 = 1.25 * C1 - q2;
    let expected_energy = 0.5 * (1.25_f64.powi(2) * C1 + 0.75_f64.powi(2) * C2);
    let mut states = Vec::new();
    for (terminals, charge) in [([1, 2], [q1, q2]), ([2, 1], [q2, q1])] {
        let floating = matrix
            .prepare_floating(&terminals, reduction_controls(), &mut || false)
            .unwrap();
        assert_eq!(floating.retained_terminals(), &[0]);
        charge_close(floating.reduced_matrix_f()[0], 0.0);
        let state = floating
            .solve(&[7.0], &charge, controls(), &mut || false)
            .unwrap();
        for (&actual, expected) in state
            .component_relative_voltage_v
            .iter()
            .zip([0.0, 1.25, 0.5])
        {
            assert!((actual - expected).abs() < 1e-12);
        }
        for (&actual, expected) in state.terminal_charge_c.iter().zip([-q1 - q2, q1, q2]) {
            charge_close(actual, expected);
        }
        near(state.field().energy_j(), expected_energy);
        near(state.trapped_charge_energy_j, expected_energy);
        near(state.retained_charge_offset_c[0], -q1 - q2);
        states.push(state);
    }
    for (&a, &b) in states[0]
        .field()
        .electric_field_v_m()
        .iter()
        .zip(states[1].field().electric_field_v_m())
    {
        near_vector(a, b);
    }
    // The original matrix reference may itself become floating; the first
    // retained electrode now supplies the physical voltage origin.
    let other_reference = matrix
        .prepare_floating(&[0, 1], reduction_controls(), &mut || false)
        .unwrap();
    assert_eq!(other_reference.reference_terminals(), &[2, 2, 2]);
    let state = other_reference
        .solve(&[0.5], &[-q1 - q2, q1], controls(), &mut || false)
        .unwrap();
    for (&actual, expected) in state
        .component_relative_voltage_v
        .iter()
        .zip([-0.5, 0.75, 0.0])
    {
        assert!((actual - expected).abs() < 1e-12);
    }
    near(state.field().energy_j(), expected_energy);
}

#[test]
fn independent_large_biases_preserve_small_trapped_charges_in_separate_components() {
    let (problem, terminals) = two_volumes(false);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    let floating = matrix
        .prepare_floating(&[0, 3], reduction_controls(), &mut || false)
        .unwrap();
    assert_eq!(floating.retained_terminals(), &[1, 2]);
    assert_eq!(floating.reference_terminals(), &[1, 1, 2, 2]);
    let charges = [C1 * 1e-6, -2.0 * C2 * 1e-6];
    let original = floating
        .solve(&[0.0, 0.0], &charges, controls(), &mut || false)
        .unwrap();
    let a = 2_f64.powi(55);
    let b = -2_f64.powi(58);
    let biased = floating
        .solve(&[a, b], &charges, controls(), &mut || false)
        .unwrap();
    assert_eq!(biased.component_reference_voltage_v, vec![a, a, b, b]);
    assert_eq!(
        biased.component_relative_voltage_v,
        original.component_relative_voltage_v
    );
    assert_eq!(biased.terminal_charge_c, original.terminal_charge_c);
    assert_eq!(
        biased.field().electric_field_v_m(),
        original.field().electric_field_v_m()
    );
    assert_eq!(
        biased.field().energy_j().to_bits(),
        original.field().energy_j().to_bits()
    );
    near(biased.component_relative_voltage_v[0], 1e-6);
    near(biased.component_relative_voltage_v[3], -2e-6);
    assert_eq!(biased.component_relative_voltage_v[1], 0.0);
    assert_eq!(biased.component_relative_voltage_v[2], 0.0);
    for (&actual, expected) in
        biased
            .terminal_charge_c
            .iter()
            .zip([charges[0], -charges[0], -charges[1], charges[1]])
    {
        near(actual, expected);
    }
    let expected_energy = 0.5 * (C1 + 4.0 * C2) * 1e-12;
    near(biased.field().energy_j(), expected_energy);
    near(biased.trapped_charge_energy_j, expected_energy);
    for (cell, &field) in biased.field().electric_field_v_m().iter().enumerate() {
        near_vector(field, [if cell < 12 { 1e-6 } else { 2e-6 }, 0.0, 0.0]);
    }
    // A reference on one component must never resolve another component's gauge.
    assert!(matches!(
        matrix.prepare_floating(&[0, 1], reduction_controls(), &mut || false),
        Err(SpatialElectrostaticError::Invalid(_))
    ));
}

#[test]
fn floating_work_limits_invalid_charge_and_cancellation_return_no_partial_state() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    for invalid in [&[][..], &[0, 1, 2], &[1, 1], &[3]] {
        assert!(
            matrix
                .prepare_floating(invalid, reduction_controls(), &mut || false)
                .is_err()
        );
    }
    for policy in [
        SpatialFloatingCapacitanceControls {
            maximum_floating_terminals: 0,
            ..reduction_controls()
        },
        SpatialFloatingCapacitanceControls {
            maximum_dense_work: 1,
            ..reduction_controls()
        },
        SpatialFloatingCapacitanceControls {
            coefficient_tolerance_f: f64::NAN,
            ..reduction_controls()
        },
    ] {
        assert!(
            matrix
                .prepare_floating(&[1], policy, &mut || false)
                .is_err()
        );
    }
    let mut preparation_visits = 0;
    let floating = matrix
        .prepare_floating(&[1], reduction_controls(), &mut || {
            preparation_visits += 1;
            false
        })
        .unwrap();
    let mut visits = 0;
    assert!(matches!(
        matrix.prepare_floating(&[1], reduction_controls(), &mut || {
            visits += 1;
            visits == preparation_visits
        }),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    assert!(
        floating
            .solve(&[0.0], &[0.0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        floating
            .solve(&[0.0, 1.0], &[], controls(), &mut || false)
            .is_err()
    );
    assert!(
        floating
            .solve(&[0.0, f64::INFINITY], &[0.0], controls(), &mut || false)
            .is_err()
    );
    assert!(
        floating
            .solve(&[0.0, 1.0], &[f64::NAN], controls(), &mut || false)
            .is_err()
    );
    assert!(
        floating
            .solve(
                &[0.0, 1.0],
                &[0.0],
                ElectrostaticControls {
                    charge_tolerance_c: 0.0,
                    ..controls()
                },
                &mut || false
            )
            .is_err()
    );
    // A successful small terminal factor does not authorize an unrepresentable
    // actual native field energy at a subsequent physical voltage amplitude.
    assert!(
        floating
            .solve(&[0.0, 1e250], &[0.0], controls(), &mut || false)
            .is_err()
    );
    let mut calls = 0;
    let expected = floating
        .solve(&[0.0, 2.0], &[0.3 * (C1 + C2)], controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            floating.solve(&[0.0, 2.0], &[0.3 * (C1 + C2)], controls(), &mut || {
                count += 1;
                count == stop
            }),
            Err(SpatialElectrostaticError::Cancelled)
        ));
    }
    let replay = floating
        .solve(&[0.0, 2.0], &[0.3 * (C1 + C2)], controls(), &mut || false)
        .unwrap();
    assert_eq!(
        expected.component_relative_voltage_v,
        replay.component_relative_voltage_v
    );
    assert_eq!(expected.terminal_charge_c, replay.terminal_charge_c);
    assert_eq!(
        expected.field().energy_j().to_bits(),
        replay.field().energy_j().to_bits()
    );
}
