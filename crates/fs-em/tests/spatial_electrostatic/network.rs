//! Multi-electrode capacitor fields against independent series/block algebra.
#[path = "network/floating.rs"]
mod floating;

use super::*;
use fs_em::spatial::electrostatic::{SpatialCapacitanceControls, SpatialElectricField};

const C1: f64 = 2e-11;
const C2: f64 = 5e-11;

fn matrix_controls() -> SpatialCapacitanceControls {
    SpatialCapacitanceControls {
        maximum_terminals: 8,
        maximum_field_solves: 7,
        maximum_response_values: 10_000,
        maximum_dense_work: 100_000,
        reference_voltage_v: 1.0,
        coefficient_tolerance_f: 1e-18,
    }
}

fn two_volumes(joined: bool) -> (SpatialElectrostatics, Vec<Vec<usize>>) {
    let (mut positions, mut cells) = cube(&[0.0, 0.5, 1.0], 1, 1);
    let offset = positions.len() as u32;
    let second = if joined { 1.0 } else { 3.0 };
    let (other, tets) = cube(&[second, second + 0.5, second + 1.0], 1, 1);
    let split = cells.len();
    positions.extend(other);
    cells.extend(tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    let material = [dielectric(C1), dielectric(C2)];
    let materials = (0..cells.len())
        .map(|i| material[usize::from(i >= split)].clone())
        .collect();
    let problem =
        SpatialElectrostatics::new(positions, cells, materials, limits(), &mut || false).unwrap();
    let xs = if joined {
        vec![0.0, 1.0, 2.0]
    } else {
        vec![0.0, 1.0, 3.0, 4.0]
    };
    let terminals = xs
        .into_iter()
        .map(|x| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0])
                        == [x; 3]
                })
                .collect()
        })
        .collect();
    (problem, terminals)
}

fn charge_close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-18, "{actual} != {expected}");
}

fn direct_field<'a>(
    problem: &'a SpatialElectrostatics,
    terminals: &[Vec<usize>],
    voltage: &[f64],
) -> SpatialElectricField<'a> {
    let boundary = SpatialElectrostaticBoundary {
        voltage: terminals
            .iter()
            .zip(voltage)
            .flat_map(|(faces, &v)| {
                faces.iter().map(move |&face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: [v; 3],
                })
            })
            .collect(),
        displacement: Vec::new(),
    };
    problem
        .solve(
            &vec![0.0; problem.mesh().tets.len()],
            &boundary,
            controls(),
            &mut || false,
        )
        .unwrap()
}

#[test]
fn three_electrode_storage_and_field_match_series_capacitors_without_new_solves() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let expected = [C1, -C1, 0.0, -C1, C1 + C2, -C2, 0.0, -C2, C2];
    for probe in [0.3, 2.0] {
        let matrix = problem
            .capacitance_matrix(
                &faces,
                SpatialCapacitanceControls {
                    reference_voltage_v: probe,
                    ..matrix_controls()
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(matrix.field_solves(), 2);
        assert_eq!(matrix.reference_terminals(), &[0, 0, 0]);
        assert_eq!(matrix.coordinate_terminals(), &[1, 2]);
        for (&actual, expected) in matrix.maxwell_matrix_f().iter().zip(expected) {
            charge_close(actual, expected);
        }
        let baseline = matrix
            .response(&[0.0, 2.0, 5.0], controls(), &mut || false)
            .unwrap();
        for common in [0.0, 2_f64.powi(45)] {
            let voltage = [common, common + 2.0, common + 5.0];
            let response = matrix
                .response(&voltage, controls(), &mut || false)
                .unwrap();
            assert!(response.field().report().is_none());
            for (&actual, expected) in
                response
                    .terminal_charge_c
                    .iter()
                    .zip([-2.0 * C1, 2.0 * C1 - 3.0 * C2, 3.0 * C2])
            {
                charge_close(actual, expected);
            }
            near(response.field().energy_j(), 0.5 * (4.0 * C1 + 9.0 * C2));
            assert_eq!(
                baseline.field().energy_j().to_bits(),
                response.field().energy_j().to_bits()
            );
            let direct = direct_field(&problem, &terminals, &voltage);
            for (cell, (&actual, &native)) in response
                .field()
                .electric_field_v_m()
                .iter()
                .zip(direct.electric_field_v_m())
                .enumerate()
            {
                let exact = if cell < 12 { -2.0 } else { -3.0 };
                near_vector(actual, [exact, 0.0, 0.0]);
                near_vector(native, [exact, 0.0, 0.0]);
                assert_eq!(actual, baseline.field().electric_field_v_m()[cell]);
            }
            assert!(response.maximum_charge_residual_c <= controls().charge_tolerance_c);
            assert!(response.terminal_energy_residual_j.abs() <= controls().energy_tolerance_j);
        }
    }
}

#[test]
fn terminal_permutation_and_face_order_preserve_the_physical_response() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    let permutation = [2, 0, 1];
    let reordered: Vec<_> = permutation
        .iter()
        .map(|&t| terminals[t].iter().copied().rev().collect::<Vec<_>>())
        .collect();
    let reordered_faces: Vec<_> = reordered.iter().map(Vec::as_slice).collect();
    let permuted = problem
        .capacitance_matrix(&reordered_faces, matrix_controls(), controls(), &mut || {
            false
        })
        .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            charge_close(
                permuted.maxwell_matrix_f()[i * 3 + j],
                matrix.maxwell_matrix_f()[permutation[i] * 3 + permutation[j]],
            );
        }
    }
    let voltage = [2.0, -3.0, 5.0];
    let reference = matrix
        .response(&voltage, controls(), &mut || false)
        .unwrap();
    let result = permuted
        .response(&permutation.map(|t| voltage[t]), controls(), &mut || false)
        .unwrap();
    near(result.field().energy_j(), reference.field().energy_j());
    for i in 0..3 {
        charge_close(
            result.terminal_charge_c[i],
            reference.terminal_charge_c[permutation[i]],
        );
    }
    for (&a, &b) in result
        .field()
        .electric_field_v_m()
        .iter()
        .zip(reference.field().electric_field_v_m())
    {
        near_vector(a, b);
    }
}

#[test]
fn disconnected_capacitors_keep_independent_large_voltage_references() {
    let (problem, terminals) = two_volumes(false);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || false)
        .unwrap();
    assert_eq!(matrix.reference_terminals(), &[0, 0, 2, 2]);
    assert_eq!(matrix.coordinate_terminals(), &[1, 3]);
    assert_eq!(matrix.field_solves(), 2);
    let a = 2_f64.powi(45);
    let b = -2_f64.powi(48);
    let plain = matrix
        .response(&[0.0, 2.0, 0.0, -3.0], controls(), &mut || false)
        .unwrap();
    let biased = matrix
        .response(&[a, a + 2.0, b, b - 3.0], controls(), &mut || false)
        .unwrap();
    assert_eq!(
        plain.field().electric_field_v_m(),
        biased.field().electric_field_v_m()
    );
    assert_eq!(
        plain.field().energy_j().to_bits(),
        biased.field().energy_j().to_bits()
    );
    for (&actual, expected) in
        biased
            .terminal_charge_c
            .iter()
            .zip([-2.0 * C1, 2.0 * C1, 3.0 * C2, -3.0 * C2])
    {
        charge_close(actual, expected);
    }
    for i in 0..2 {
        for j in 2..4 {
            assert_eq!(matrix.maxwell_matrix_f()[4 * i + j], 0.0);
        }
    }
}

#[test]
fn multi_electrode_admission_and_cancellation_leave_extraction_reusable() {
    let (problem, terminals) = two_volumes(true);
    let faces: Vec<_> = terminals.iter().map(Vec::as_slice).collect();
    for bad in [
        SpatialCapacitanceControls {
            maximum_field_solves: 1,
            ..matrix_controls()
        },
        SpatialCapacitanceControls {
            maximum_response_values: 1,
            ..matrix_controls()
        },
        SpatialCapacitanceControls {
            maximum_dense_work: 1,
            ..matrix_controls()
        },
        SpatialCapacitanceControls {
            reference_voltage_v: f64::NAN,
            ..matrix_controls()
        },
    ] {
        assert!(
            problem
                .capacitance_matrix(&faces, bad, controls(), &mut || false)
                .is_err()
        );
    }
    assert!(
        problem
            .capacitance_matrix(
                &[faces[0], faces[0]],
                matrix_controls(),
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        problem
            .capacitance_matrix(&[faces[0], &[]], matrix_controls(), controls(), &mut || {
                false
            })
            .is_err()
    );
    let mut calls = 0;
    let matrix = problem
        .capacitance_matrix(&faces, matrix_controls(), controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        problem.capacitance_matrix(&faces, matrix_controls(), controls(), &mut || {
            count += 1;
            count >= calls / 2
        }),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    for voltage in [vec![1.0], vec![0.0, f64::NAN, 1.0]] {
        assert!(
            matrix
                .response(&voltage, controls(), &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        matrix.response(&[0.0, 1.0, 2.0], controls(), &mut || true),
        Err(SpatialElectrostaticError::Cancelled)
    ));
    let replay = matrix
        .response(&[0.0, 1.0, 2.0], controls(), &mut || false)
        .unwrap();
    near(replay.field().energy_j(), 0.5 * (C1 + C2));
}
