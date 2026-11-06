//! Three-terminal bulk/contact resistance against independent series algebra.
use super::*;
use fs_em::conduction::matrix::{PlanarConductanceControls, PlanarMixedConductanceControls};

fn response_controls() -> PlanarConductanceControls {
    PlanarConductanceControls {
        maximum_terminals: 4,
        maximum_field_solves: 3,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 1.0,
        coefficient_tolerance_s: 1e-10,
    }
}

fn mixed_controls() -> PlanarMixedConductanceControls {
    PlanarMixedConductanceControls {
        maximum_current_terminals: 2,
        maximum_response_values: 100,
        maximum_dense_work: 100,
    }
}

fn specimen() -> PlanarConduction {
    // First path: two 1 m bodies, sigma=2 S/m, A=.5 m²: R_bulk=2 ohm.
    // Their source R''=.5 ohm m² contributes another 1 ohm. The third
    // body is 1 m long with sigma=2/3 S/m and the same area: R=3 ohm.
    // One middle ideal electrode joins their disconnected exterior traces.
    let material = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let right = conductor(2.0 / 3.0, ElectricalPropertyBasis::Resistivity);
    let first = contacts::mesh(0.5);
    let second = layered(0.5, 1.0, 0.5);
    let mut positions = first.vertices().to_vec();
    positions.extend(second.vertices().iter().map(|p| [p[0] + 2.0, p[1]]));
    let mut faces = first.faces().to_vec();
    faces.extend(second.faces().iter().map(|face| face.map(|v| v + 8)));
    let mesh = TriComplex2::from_indexed_triangles(
        "three-terminal source contact",
        positions,
        faces,
        Metric2::planar(0.5).unwrap(),
    )
    .unwrap();
    let contact = contacts::binding(contacts::interface(&material, &material, 1e5, 1.0));
    PlanarConduction::new_with_contacts(
        mesh,
        (0..8)
            .map(|i| {
                if i < 4 {
                    material.clone()
                } else {
                    right.clone()
                }
            })
            .collect(),
        vec![contact],
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn electrodes() -> [Vec<usize>; 3] {
    [vec![0, 3], vec![5, 6, 8, 11], vec![10, 13]]
}

fn reaction(field: &fs_em::conduction::ConductionField<'_>, terminal: &[usize]) -> f64 {
    terminal
        .iter()
        .map(|&v| field.nodal_boundary_current_a()[v].unwrap())
        .sum()
}

#[test]
fn multi_terminal_matrix_keeps_local_cross_terms_and_original_contact_heat() {
    let problem = specimen();
    let electrodes = electrodes();
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let expected = [
        1.0 / 3.0,
        -1.0 / 3.0,
        0.0,
        -1.0 / 3.0,
        2.0 / 3.0,
        -1.0 / 3.0,
        0.0,
        -1.0 / 3.0,
        1.0 / 3.0,
    ];
    for reference_voltage_v in [0.3, 1.0, 4.0] {
        let matrix = problem
            .conductance_matrix(
                &terminals,
                PlanarConductanceControls {
                    reference_voltage_v,
                    ..response_controls()
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(matrix.field_solves(), 2);
        assert_eq!(matrix.terminal_vertices(), &electrodes);
        assert!(matrix.maximum_coefficient_residual_s() < 1e-10);
        for (&value, exact) in matrix.conductance_s().iter().zip(expected) {
            near(value, exact);
        }
        for sign in [-1.0_f64, 1.0] {
            for bias in [0.0, 2.0_f64.powi(40)] {
                let voltages = [bias, bias + 2.0 * sign, bias + 5.0 * sign];
                let field = matrix
                    .field_at_voltages(&voltages, 1e-10, 1e-9, &mut || false)
                    .unwrap();
                assert!(field.report().is_none());
                for (vertices, exact) in
                    terminals.iter().zip([-2.0 * sign / 3.0, -sign / 3.0, sign])
                {
                    near(reaction(&field, vertices), exact);
                }
                near(field.voltage_v()[0], 0.0);
                near(field.voltage_v()[5], 2.0 * sign);
                near(field.voltage_v()[10], 5.0 * sign);
                near(field.power_w(), 13.0 / 3.0);
                near(field.contact_response()[0].power_w, 4.0 / 9.0);
                for &power in &field.cell_power_w()[..4] {
                    near(power, 2.0 / 9.0);
                }
                for &power in &field.cell_power_w()[4..] {
                    near(power, 0.75);
                }
                // Independent native solve at the same simultaneous controls.
                let boundary: Vec<_> = terminals
                    .iter()
                    .zip(voltages)
                    .flat_map(|(vertices, voltage)| vertices.iter().map(move |&v| (v, voltage)))
                    .collect();
                let direct = problem.solve(&boundary, controls(), &mut || false).unwrap();
                near(direct.power_w(), field.power_w());
                for (&a, &b) in direct.cell_power_w().iter().zip(field.cell_power_w()) {
                    near(a, b);
                }
            }
        }
        // Joule heat is quadratic in the superposed field, with real cross terms.
        let first = matrix
            .field_at_voltages(&[0.0, 2.0, 0.0], 1e-10, 1e-9, &mut || false)
            .unwrap();
        let second = matrix
            .field_at_voltages(&[0.0, 0.0, 5.0], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(first.power_w() + second.power_w() - 13.0 / 3.0, 20.0 / 3.0);
        let zero = matrix
            .field_at_voltages(&[1e200; 3], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(zero.power_w(), 0.0);
    }
}

#[test]
fn reordered_terminals_retain_cross_terms_in_the_contact_itself() {
    let problem = specimen();
    let electrodes = electrodes();
    let terminals = [&electrodes[2][..], &electrodes[0][..], &electrodes[1][..]];
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    let mut powers = Vec::new();
    for (voltages, expected) in [
        ([0.0, -5.0, -3.0], 4.0 / 9.0),
        ([0.0, -5.0, 0.0], 25.0 / 9.0),
        ([0.0, 0.0, -3.0], 1.0),
    ] {
        let field = matrix
            .field_at_voltages(&voltages, 1e-10, 1e-9, &mut || false)
            .unwrap();
        let heat = field.contact_response()[0].power_w;
        near(heat, expected);
        powers.push(heat);
    }
    near(powers[1] + powers[2] - powers[0], 10.0 / 3.0);
}

#[test]
fn total_current_and_mixed_controls_keep_a_floating_electrode_locally_conducting() {
    let problem = specimen();
    let electrodes = electrodes();
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    for sign in [-1.0_f64, 1.0] {
        let field = matrix
            .field_at_currents(&[-sign, 0.0, sign], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(field.voltage_v()[5], 3.0 * sign);
        near(field.voltage_v()[10], 6.0 * sign);
        near(field.power_w(), 6.0);
        near(field.contact_response()[0].power_w, 1.0);
        near(reaction(&field, &[5, 6]), sign);
        near(reaction(&field, &[8, 11]), -sign);
        near(reaction(&field, &terminals[1]), 0.0);
    }
    let mixed = matrix
        .prepare_current_control(&[1], mixed_controls(), &mut || false)
        .unwrap();
    assert_eq!(mixed.voltage_terminals(), &[0, 2]);
    assert_eq!(mixed.current_terminals(), &[1]);
    for current in [-0.5_f64, 0.0, 0.5] {
        let middle = 3.0 + 1.5 * current;
        let field = mixed
            .field_at_controls(&[1e10, 1e10 + 6.0], &[current], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(field.voltage_v()[5], middle);
        near(reaction(&field, terminals[1]), current);
        near(
            field.power_w(),
            middle.powi(2) / 3.0 + (6.0 - middle).powi(2) / 3.0,
        );
    }
    // Terminal zero may itself be current controlled; input order remains authored.
    let reordered = matrix
        .prepare_current_control(&[1, 0], mixed_controls(), &mut || false)
        .unwrap();
    assert_eq!(reordered.voltage_terminals(), &[2]);
    let field = reordered
        .field_at_controls(&[9e20], &[0.0, -1.0], 1e-10, 1e-9, &mut || false)
        .unwrap();
    near(field.voltage_v()[0], 0.0);
    near(field.voltage_v()[5], 3.0);
    near(field.voltage_v()[10], 6.0);
}

#[test]
fn terminal_preparation_refuses_disconnected_controls_and_honors_work_bounds() {
    let problem = specimen();
    let electrodes = electrodes();
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let disconnected: [&[usize]; 4] = [&[0, 3], &[5, 6], &[8, 11], &[10, 13]];
    assert!(
        problem
            .conductance_matrix(&disconnected, response_controls(), controls(), &mut || {
                false
            })
            .is_err()
    );
    let point: [&[usize]; 3] = [&[0], &[5, 6, 8, 11], &[10, 13]];
    assert!(
        problem
            .conductance_matrix(&point, response_controls(), controls(), &mut || false)
            .is_err()
    );
    for bounded in [
        PlanarConductanceControls {
            maximum_terminals: 2,
            ..response_controls()
        },
        PlanarConductanceControls {
            maximum_field_solves: 1,
            ..response_controls()
        },
        PlanarConductanceControls {
            maximum_response_values: 1,
            ..response_controls()
        },
        PlanarConductanceControls {
            maximum_dense_work: 1,
            ..response_controls()
        },
        PlanarConductanceControls {
            reference_voltage_v: 0.0,
            ..response_controls()
        },
        PlanarConductanceControls {
            coefficient_tolerance_s: f64::NAN,
            ..response_controls()
        },
    ] {
        assert!(
            problem
                .conductance_matrix(&terminals, bounded, controls(), &mut || false)
                .is_err()
        );
    }
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    for bad in [&[][..], &[0, 0][..], &[3][..], &[0, 1, 2][..]] {
        assert!(
            matrix
                .prepare_current_control(bad, mixed_controls(), &mut || false)
                .is_err()
        );
    }
    assert!(
        matrix
            .prepare_current_control(
                &[1],
                PlanarMixedConductanceControls {
                    maximum_dense_work: 0,
                    ..mixed_controls()
                },
                &mut || false
            )
            .is_err()
    );
    assert!(
        matrix
            .field_at_currents(&[0.0, 0.0, 1.0], 1e-10, 1e-9, &mut || false)
            .is_err()
    );
    assert!(
        matrix
            .field_at_voltages(&[0.0, f64::NAN, 1.0], 1e-10, 1e-9, &mut || false)
            .is_err()
    );
    assert!(
        matrix
            .field_at_currents(&[-1.0, 0.0, 1.0], 0.0, 1e-9, &mut || false)
            .is_err()
    );
}

#[test]
fn prepared_response_survives_late_cancellation_and_retains_original_sources() {
    let problem = specimen();
    let electrodes = electrodes();
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    let mixed = matrix
        .prepare_current_control(&[1], mixed_controls(), &mut || false)
        .unwrap();
    let mut calls = 0;
    let baseline = mixed
        .field_at_controls(&[0.0, 6.0], &[0.5], 1e-10, 1e-9, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut attempted = 0;
    assert!(matches!(
        mixed.field_at_controls(&[0.0, 6.0], &[0.5], 1e-10, 1e-9, &mut || {
            attempted += 1;
            attempted == calls - 1
        }),
        Err(ConductionError::Cancelled)
    ));
    let retry = mixed
        .field_at_controls(&[0.0, 6.0], &[0.5], 1e-10, 1e-9, &mut || false)
        .unwrap();
    assert_eq!(retry.voltage_v(), baseline.voltage_v());
    assert_eq!(retry.cell_power_w(), baseline.cell_power_w());
    near(
        retry.contact_response()[0].power_w,
        baseline.contact_response()[0].power_w,
    );
    assert!(matches!(
        problem.conductance_matrix(&terminals, response_controls(), controls(), &mut || true),
        Err(ConductionError::Cancelled)
    ));
}
