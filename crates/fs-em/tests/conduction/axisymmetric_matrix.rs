//! Multi-electrode meridional conduction against cylinder and coaxial laws.
//! Synthetic SI fixtures exercise actual swept area, floating conductors and
//! quadratic local heat; these are not measured material qualifications.

use super::*;
use fs_em::conduction::axisymmetric::{
    AxisymmetricConductanceControls, AxisymmetricConduction, AxisymmetricConductionField,
    AxisymmetricMixedConductanceControls,
};

fn response_controls() -> AxisymmetricConductanceControls {
    AxisymmetricConductanceControls {
        maximum_terminals: 4,
        maximum_field_solves: 3,
        maximum_response_values: 100_000,
        maximum_dense_work: 100_000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-9,
    }
}

fn mixed_controls() -> AxisymmetricMixedConductanceControls {
    AxisymmetricMixedConductanceControls {
        maximum_current_terminals: 2,
        maximum_response_values: 100,
        maximum_dense_work: 100,
    }
}

fn join(first: TriComplex2, second: TriComplex2, axial_shift: f64, angle: f64) -> TriComplex2 {
    let offset = first.vertices().len() as u32;
    let mut vertices = first.vertices().to_vec();
    vertices.extend(second.vertices().iter().map(|p| [p[0], p[1] + axial_shift]));
    let mut faces = first.faces().to_vec();
    faces.extend(second.faces().iter().map(|face| face.map(|v| v + offset)));
    TriComplex2::from_indexed_triangles(
        "axisymmetric ideal-electrode connection",
        vertices,
        faces,
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn axial_chain(scale: f64, angle: f64) -> (AxisymmetricConduction, [Vec<usize>; 3], [f64; 2]) {
    let radius = 0.6 * scale;
    let lengths = [0.4 * scale, 0.9 * scale];
    let sigma = [2.0, 0.5];
    let mesh = join(
        axisymmetric::rectangle(0.0, radius, lengths[0], 2, 2, angle),
        axisymmetric::rectangle(0.0, radius, lengths[1], 2, 2, angle),
        lengths[0],
        angle,
    );
    let materials = [
        conductor(sigma[0], ElectricalPropertyBasis::Conductivity),
        conductor(sigma[1], ElectricalPropertyBasis::Resistivity),
    ];
    let problem = AxisymmetricConduction::new(
        mesh,
        (0..16)
            .map(|cell| materials[usize::from(cell >= 8)].clone())
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    // The middle ideal disk electrode joins two originally separate bodies.
    let terminals = [vec![0, 1, 2], vec![6, 7, 8, 9, 10, 11], vec![15, 16, 17]];
    let area = angle * radius * radius / 2.0;
    let resistance = [
        lengths[0] / (sigma[0] * area),
        lengths[1] / (sigma[1] * area),
    ];
    (problem, terminals, resistance)
}

fn reaction(field: &AxisymmetricConductionField<'_>, terminal: &[usize]) -> f64 {
    terminal
        .iter()
        .map(|&v| field.nodal_boundary_current_a()[v].unwrap())
        .sum()
}

fn chain_matrix(resistance: [f64; 2]) -> [f64; 9] {
    let [a, b] = resistance.map(f64::recip);
    [a, -a, 0.0, -a, a + b, -b, 0.0, -b, b]
}

#[test]
fn cylinder_terminal_matrix_retains_swept_geometry_and_simultaneous_heat_cross_terms() {
    for (scale, angle) in [(1.0, core::f64::consts::TAU), (0.2, 0.7)] {
        let (problem, electrodes, resistance) = axial_chain(scale, angle);
        let terminals = electrodes.each_ref().map(Vec::as_slice);
        let matrix = problem
            .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
            .unwrap();
        assert_eq!(matrix.field_solves(), 2);
        assert_eq!(matrix.terminal_vertices(), &electrodes);
        assert_eq!(matrix.problem().conductors(), problem.conductors());
        assert!(matrix.maximum_coefficient_residual_s() < 1e-9);
        for (&actual, expected) in matrix.conductance_s().iter().zip(chain_matrix(resistance)) {
            near(actual, expected);
        }
        let expected_heat = 1.25_f64.powi(2) / resistance[0] + 2.0_f64.powi(2) / resistance[1];
        for sign in [-1.0, 1.0] {
            for bias in [0.0, 2.0_f64.powi(40)] {
                let field = matrix
                    .field_at_voltages(
                        &[bias, bias + 1.25 * sign, bias - 0.75 * sign],
                        1e-9,
                        1e-9,
                        &mut || false,
                    )
                    .unwrap();
                assert!(field.report().is_none());
                let currents = [
                    -1.25 * sign / resistance[0],
                    1.25 * sign / resistance[0] + 2.0 * sign / resistance[1],
                    -2.0 * sign / resistance[1],
                ];
                for (terminal, expected) in terminals.iter().zip(currents) {
                    near(reaction(&field, terminal), expected);
                }
                near(field.voltage_v()[0], 0.0);
                near(field.voltage_v()[6], 1.25 * sign);
                near(field.voltage_v()[15], -0.75 * sign);
                let gradients = [1.25 * sign / (0.4 * scale), -2.0 * sign / (0.9 * scale)];
                for cell in 0..problem.mesh().faces().len() {
                    let region = usize::from(cell >= 8);
                    let sigma = [2.0, 0.5][region];
                    near(field.electric_field_v_m()[cell][0], 0.0);
                    near(field.electric_field_v_m()[cell][1], -gradients[region]);
                    near(field.current_density_a_m2()[cell][0], 0.0);
                    near(
                        field.current_density_a_m2()[cell][1],
                        -sigma * gradients[region],
                    );
                    near(
                        field.cell_power_w()[cell],
                        sigma
                            * gradients[region].powi(2)
                            * problem.mesh().face_measure(cell).unwrap(),
                    );
                }
                near(field.power_w(), expected_heat);
                near(field.boundary_current_a(), 0.0);
                near(field.power_balance_residual_w(), 0.0);
            }
        }
        let first = matrix
            .field_at_voltages(&[0.0, 1.25, 0.0], 1e-9, 1e-9, &mut || false)
            .unwrap();
        let second = matrix
            .field_at_voltages(&[0.0, 0.0, -0.75], 1e-9, 1e-9, &mut || false)
            .unwrap();
        near(
            expected_heat - first.power_w() - second.power_w(),
            2.0 * 1.25 * 0.75 / resistance[1],
        );
        let zero = matrix
            .field_at_voltages(&[1e200; 3], 1e-9, 1e-9, &mut || false)
            .unwrap();
        near(zero.power_w(), 0.0);
    }
}

#[test]
fn current_and_mixed_controls_transfer_current_through_a_floating_disk() {
    let (problem, electrodes, resistance) = axial_chain(1.0, 1.7);
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    for current in [-0.3, 0.0, 0.2] {
        let field = matrix
            .field_at_currents(&[-current, 0.0, current], 1e-9, 1e-9, &mut || false)
            .unwrap();
        assert!(field.report().is_none());
        near(field.voltage_v()[6], current * resistance[0]);
        near(
            field.voltage_v()[15],
            current * resistance.iter().sum::<f64>(),
        );
        near(reaction(&field, &electrodes[1][..3]), current);
        near(reaction(&field, &electrodes[1][3..]), -current);
        near(reaction(&field, terminals[1]), 0.0);
        for (region, &r) in resistance.iter().enumerate() {
            near(
                field.cell_power_w()[region * 8..(region + 1) * 8]
                    .iter()
                    .sum(),
                current * current * r,
            );
        }
    }
    let mixed = matrix
        .prepare_current_control(&[1], mixed_controls(), &mut || false)
        .unwrap();
    assert_eq!(mixed.voltage_terminals(), &[0, 2]);
    assert_eq!(mixed.current_terminals(), &[1]);
    for current in [-0.1, 0.0, 0.1] {
        let middle = (current + 2.0 / resistance[1]) / (1.0 / resistance[0] + 1.0 / resistance[1]);
        let field = mixed
            .field_at_controls(
                &[2.0_f64.powi(40), 2.0_f64.powi(40) + 2.0],
                &[current],
                1e-9,
                1e-9,
                &mut || false,
            )
            .unwrap();
        near(field.voltage_v()[6], middle);
        near(reaction(&field, terminals[1]), current);
        near(
            field.power_w(),
            middle.powi(2) / resistance[0] + (2.0 - middle).powi(2) / resistance[1],
        );
    }
    let reordered = matrix
        .prepare_current_control(&[2, 0], mixed_controls(), &mut || false)
        .unwrap();
    assert_eq!(reordered.voltage_terminals(), &[1]);
    assert_eq!(reordered.current_terminals(), &[2, 0]);
    let mut polls = 0;
    let expected = reordered
        .field_at_controls(&[1e200], &[0.2, -0.2], 1e-9, 1e-9, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    near(expected.voltage_v()[0], 0.0);
    near(expected.voltage_v()[6], 0.2 * resistance[0]);
    near(
        expected.voltage_v()[15],
        0.2 * resistance.iter().sum::<f64>(),
    );
    let mut attempts = 0;
    assert!(matches!(
        reordered.field_at_controls(&[1e200], &[0.2, -0.2], 1e-9, 1e-9, &mut || {
            attempts += 1;
            attempts == polls - 1
        }),
        Err(ConductionError::Cancelled)
    ));
    let retry = reordered
        .field_at_controls(&[1e200], &[0.2, -0.2], 1e-9, 1e-9, &mut || false)
        .unwrap();
    assert_eq!(retry.voltage_v(), expected.voltage_v());
    assert_eq!(retry.cell_power_w(), expected.cell_power_w());
}

#[test]
fn coaxial_ring_matrix_and_floating_voltage_converge_to_independent_logarithmic_law() {
    let (angle, length) = (1.3, 0.8);
    let sigma = [0.25, 0.75];
    let exact_r = [2.0_f64.ln() / sigma[0], 1.5_f64.ln() / sigma[1]].map(|r| r / (angle * length));
    let expected_g = chain_matrix(exact_r);
    let field_controls = ConductionControls {
        linear: fs_solver::CheckedCgConfig {
            tolerance: 1e-13,
            ..config()
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    };
    let mut previous = f64::INFINITY;
    for nr in [8, 16, 32] {
        let nz = 4;
        let first = axisymmetric::rectangle(0.5, 1.0, length, nr, nz, angle);
        let second = axisymmetric::rectangle(1.0, 1.5, length, nr, nz, angle);
        let offset = first.vertices().len();
        let cells = first.faces().len();
        let mesh = join(first, second, 0.0, angle);
        let materials = sigma.map(|s| conductor(s, ElectricalPropertyBasis::Conductivity));
        let p = AxisymmetricConduction::new(
            mesh,
            (0..2 * cells)
                .map(|cell| materials[usize::from(cell >= cells)].clone())
                .collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let inner: Vec<_> = (0..=nz).map(|j| j * (nr + 1)).collect();
        let middle: Vec<_> = (0..=nz)
            .map(|j| j * (nr + 1) + nr)
            .chain((0..=nz).map(|j| offset + j * (nr + 1)))
            .collect();
        let outer: Vec<_> = (0..=nz).map(|j| offset + j * (nr + 1) + nr).collect();
        let terminals = [&inner[..], &middle[..], &outer[..]];
        let matrix = p
            .conductance_matrix(&terminals, response_controls(), field_controls, &mut || {
                false
            })
            .unwrap();
        let mut error = 0.0_f64;
        for (&actual, expected) in matrix.conductance_s().iter().zip(expected_g) {
            if expected == 0.0 {
                near(actual, 0.0);
            } else {
                error = error.max((actual / expected - 1.0).abs());
            }
        }
        assert!(
            error < 0.6 * previous,
            "nr={nr}: error={error}, previous={previous}"
        );
        previous = error;
        let mixed = matrix
            .prepare_current_control(&[1], mixed_controls(), &mut || false)
            .unwrap();
        let f = mixed
            .field_at_controls(&[1.0, 0.0], &[0.0], 1e-9, 1e-9, &mut || false)
            .unwrap();
        let current = 1.0 / exact_r.iter().sum::<f64>();
        // Terminal zero is the inner ring, so all reported voltages are
        // shifted by -1 V; fields, current and heat retain the physical drive.
        let expected_middle = exact_r[1] * current - 1.0;
        assert!((f.voltage_v()[middle[0]] - expected_middle).abs() < 4.0 / (nr * nr) as f64);
        for (region, &r) in exact_r.iter().enumerate() {
            let heat: f64 = f.cell_power_w()[region * cells..(region + 1) * cells]
                .iter()
                .sum();
            assert!((heat / (current * current * r) - 1.0).abs() < 4.0 / (nr * nr) as f64);
        }
        near(reaction(&f, &middle), 0.0);
        assert!(reaction(&f, &middle[..nz + 1]).abs() > 0.0);
        near(f.power_balance_residual_w(), 0.0);
    }
    assert!(previous < 5e-4);
}

#[test]
fn terminal_admission_refuses_zero_area_disconnection_and_unbounded_or_unbalanced_controls() {
    let (problem, electrodes, _) = axial_chain(1.0, 1.7);
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let axis: [&[usize]; 3] = [&[0, 3, 6], &[2, 5, 8], &[11, 14, 17]];
    assert!(
        problem
            .conductance_matrix(&axis, response_controls(), controls(), &mut || false)
            .is_err()
    );
    let disconnected: [&[usize]; 4] = [
        &electrodes[0],
        &electrodes[1][..3],
        &electrodes[1][3..],
        &electrodes[2],
    ];
    assert!(
        problem
            .conductance_matrix(&disconnected, response_controls(), controls(), &mut || {
                false
            })
            .is_err()
    );
    for bounded in [
        AxisymmetricConductanceControls {
            maximum_terminals: 2,
            ..response_controls()
        },
        AxisymmetricConductanceControls {
            maximum_field_solves: 1,
            ..response_controls()
        },
        AxisymmetricConductanceControls {
            maximum_response_values: 1,
            ..response_controls()
        },
        AxisymmetricConductanceControls {
            maximum_dense_work: 1,
            ..response_controls()
        },
    ] {
        assert!(
            problem
                .conductance_matrix(&terminals, bounded, controls(), &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        problem.conductance_matrix(&terminals, response_controls(), controls(), &mut || true),
        Err(ConductionError::Cancelled)
    ));
    let matrix = problem
        .conductance_matrix(&terminals, response_controls(), controls(), &mut || false)
        .unwrap();
    assert!(
        matrix
            .field_at_currents(&[0.0, 0.0, 0.1], 1e-9, 1e-9, &mut || false)
            .is_err()
    );
    assert!(
        matrix
            .field_at_voltages(&[0.0, f64::NAN, 1.0], 1e-9, 1e-9, &mut || false)
            .is_err()
    );
    for current_terminals in [&[][..], &[0, 0][..], &[3][..], &[0, 1, 2][..]] {
        assert!(
            matrix
                .prepare_current_control(current_terminals, mixed_controls(), &mut || false)
                .is_err()
        );
    }
    let floating = matrix
        .prepare_current_control(&[1], mixed_controls(), &mut || false)
        .unwrap();
    assert!(
        floating
            .field_at_controls(&[0.0], &[0.0], 1e-9, 1e-9, &mut || false)
            .is_err()
    );
}
