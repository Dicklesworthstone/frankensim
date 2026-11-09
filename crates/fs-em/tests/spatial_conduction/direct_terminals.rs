//! G1/G3/G4: direct tetrahedral electrode controls against series/tensor laws.
//! Manufactured sources verify discretization and accounting, not empirical data.
use super::contacts::{conductor, contacts as material_contacts, interface};
use super::*;
use fs_em::conduction::terminals::{ConductionTerminalControl, ConductionTerminalLimits};
use std::collections::BTreeSet;

fn terminal_limits() -> ConductionTerminalLimits {
    ConductionTerminalLimits {
        maximum_terminals: 3,
        maximum_workspace_values: 16_384,
    }
}

fn specimen() -> (SpatialConduction, [Vec<usize>; 3]) {
    // A: unit bulk resistance plus one unit contact resistance. B: 3 ohm.
    // Only the middle ideal electrode joins A and B's disconnected meshes.
    let materials = [conductor(1.0), conductor(3.0)];
    let source = interface(&materials[0], &materials[0], 1e5, 1.0 / 3.0);
    let (mut positions, mut tets) = slab(&[0.0, 0.5]);
    for xs in [&[0.5, 1.0][..], &[1.0, 1.5, 2.0][..]] {
        let (part, cells) = slab(xs);
        let offset = positions.len() as u32;
        positions.extend(part);
        tets.extend(cells.into_iter().map(|tet| tet.map(|v| v + offset)));
    }
    let problem = SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..24)
            .map(|i| materials[usize::from(i >= 12)].clone())
            .collect(),
        material_contacts(source),
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = [0.0, 1.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect()
    });
    (problem, faces)
}

fn reaction(problem: &SpatialConduction, field: &SpatialConductionField, faces: &[usize]) -> f64 {
    faces
        .iter()
        .flat_map(|&face| problem.mesh().faces[face])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|v| field.nodal_boundary_current_a[v as usize].unwrap())
        .sum()
}

#[test]
fn direct_controls_recover_native_contact_heat_and_floating_face_transfer() {
    let (problem, faces) = specimen();
    let terminals = faces.each_ref().map(Vec::as_slice);
    let originals = problem.conductors().to_vec();
    let matrix = problem
        .conductance_matrix(
            &terminals,
            SpatialConductanceControls {
                maximum_terminals: 3,
                maximum_field_solves: 2,
                maximum_response_values: 1000,
                maximum_dense_work: 1000,
                reference_voltage_v: 1.0,
                coefficient_tolerance_s: 1e-10,
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    for sign in [-1.0_f64, 1.0] {
        for injection in [0.0, 0.5] {
            let middle = sign * (2.0 + 1.2 * injection);
            let right = sign * 5.0;
            let left_current = -middle / 2.0;
            let right_current = (right - middle) / 3.0;
            let currents = [left_current, sign * injection, right_current];
            let voltages = [0.0, middle, right];
            let outer = [0.0, right];
            let reordered_currents = [right_current, sign * injection];
            let expected_power = 2.0 * left_current.powi(2) + 3.0 * right_current.powi(2);
            let basis = matrix
                .field_at_voltages(&voltages, 1e-10, 1e-9, &mut || false)
                .unwrap();
            for control in [
                ConductionTerminalControl::Voltages(&voltages),
                ConductionTerminalControl::Currents(&currents),
                ConductionTerminalControl::Mixed {
                    current_terminals: &[1],
                    currents_a: &[sign * injection],
                    voltages_v: &outer,
                },
                ConductionTerminalControl::Mixed {
                    current_terminals: &[2, 1],
                    currents_a: &reordered_currents,
                    voltages_v: &[1e250],
                },
            ] {
                let (field, report) = problem
                    .solve_terminals(
                        &terminals,
                        control,
                        terminal_limits(),
                        controls(),
                        &mut || false,
                    )
                    .unwrap();
                assert!(field.linear.is_some());
                for i in 0..3 {
                    near(report.voltage_v[i], voltages[i]);
                    near(report.current_a[i], currents[i]);
                    near(reaction(&problem, &field, &faces[i]), currents[i]);
                }
                // A floating electrode has zero TOTAL current, while each
                // physical side still transfers a full nonzero branch current.
                near(
                    field.nodal_boundary_current_a[12..16]
                        .iter()
                        .map(|v| v.unwrap())
                        .sum(),
                    -left_current,
                );
                near(
                    field.nodal_boundary_current_a[16..20]
                        .iter()
                        .map(|v| v.unwrap())
                        .sum(),
                    -right_current,
                );
                for (cell, (actual, reference)) in field
                    .cell_power_w
                    .iter()
                    .zip(&basis.cell_power_w)
                    .enumerate()
                {
                    near(*actual, *reference);
                    near(
                        field.current_density_a_m2[cell][0],
                        if cell < 12 {
                            left_current
                        } else {
                            -right_current
                        },
                    );
                    near(field.current_density_a_m2[cell][1], 0.0);
                    near(field.current_density_a_m2[cell][2], 0.0);
                }
                near(field.cell_power_w[..12].iter().sum(), left_current.powi(2));
                near(
                    field.cell_power_w[12..].iter().sum(),
                    3.0 * right_current.powi(2),
                );
                near(
                    field.contacts.iter().map(|contact| contact.power_w).sum(),
                    left_current.powi(2),
                );
                near(field.power_w, expected_power);
                near(field.supplied_power_w, expected_power);
                near(report.supplied_power_w, expected_power);
                near(report.power_residual_w, 0.0);
                assert!(report.maximum_current_residual_a <= controls().current_tolerance_a);
                near(
                    report.supplied_power_w,
                    middle * sign * injection + right * right_current,
                );
            }
        }
    }
    let permutation = [&faces[1][..], &faces[2][..], &faces[0][..]];
    let (_, report) = problem
        .solve_terminals(
            &permutation,
            ConductionTerminalControl::Currents(&[0.0, 1.0, -1.0]),
            terminal_limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&actual, expected) in report.voltage_v.iter().zip([0.0, 3.0, -2.0]) {
        near(actual, expected);
    }
    near(report.supplied_power_w, 5.0);
    assert_eq!(problem.conductors(), originals);
}

#[test]
fn rotated_tensor_terminals_preserve_physical_field_under_huge_common_reference() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    for basis in [
        ElectricalPropertyBasis::Conductivity,
        ElectricalPropertyBasis::Resistivity,
    ] {
        let values = if basis == ElectricalPropertyBasis::Conductivity {
            [2.0, 8.0, 3.0]
        } else {
            [0.5, 0.125, 1.0 / 3.0]
        };
        let source = ResolvedOrthotropicConductor3::try_from_material(
            material(basis, &values),
            basis,
            std::array::from_fn(|i| PropertyKey::new(format!("principal_{i}"), basis.dimensions())),
            axes,
        )
        .unwrap();
        let (mut positions, tets) = slab(&[0.0, 0.5, 1.0]);
        for p in &mut positions {
            *p = [s * (p[0] - p[1]), s * (p[0] + p[1]), p[2]];
        }
        let problem =
            SpatialConduction::new(positions, tets, vec![source; 12], limits(), &mut || false)
                .unwrap();
        let faces = [0_u32, 8].map(|first| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face]
                        .iter()
                        .all(|&v| v >= first && v < first + 4)
                })
                .collect::<Vec<_>>()
        });
        let terminals = faces.each_ref().map(Vec::as_slice);
        for common in [0.0, 1e250] {
            let (field, report) = problem
                .solve_terminals(
                    &terminals,
                    ConductionTerminalControl::Mixed {
                        current_terminals: &[1],
                        currents_a: &[4.0],
                        voltages_v: &[common],
                    },
                    terminal_limits(),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(report.voltage_v[0], 0.0);
            near(report.voltage_v[1], 2.0);
            near(report.current_a[0], -4.0);
            near(report.current_a[1], 4.0);
            near(report.supplied_power_w, 8.0);
            for cell in 0..12 {
                for axis in 0..3 {
                    near(field.electric_field_v_m[cell][axis], -2.0 * axes[0][axis]);
                    near(field.current_density_a_m2[cell][axis], -4.0 * axes[0][axis]);
                }
                near(field.cell_power_w[cell], 2.0 / 3.0);
            }
            near(field.power_w, 8.0);
            near(field.power_residual_w, 0.0);
        }
    }
}

#[test]
fn direct_terminal_admission_and_cancelled_trial_leave_original_reusable() {
    let (problem, faces) = specimen();
    let terminals = faces.each_ref().map(Vec::as_slice);
    let drive = ConductionTerminalControl::Currents(&[-1.0, 0.0, 1.0]);
    let mut polls = 0;
    let (before, before_report) = problem
        .solve_terminals(
            &terminals,
            drive,
            terminal_limits(),
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop_at in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            problem.solve_terminals(
                &terminals,
                drive,
                terminal_limits(),
                controls(),
                &mut || {
                    calls += 1;
                    calls == stop_at
                },
            ),
            Err(SpatialConductionError::Cancelled)
        ));
    }
    for invalid in [
        ConductionTerminalControl::Currents(&[-1.0, 0.0, 2.0]),
        ConductionTerminalControl::Currents(&[-1.0, f64::NAN, 1.0]),
        ConductionTerminalControl::Voltages(&[0.0, 5.0]),
        ConductionTerminalControl::Mixed {
            current_terminals: &[1, 1],
            currents_a: &[0.0, 0.0],
            voltages_v: &[0.0],
        },
    ] {
        assert!(
            problem
                .solve_terminals(
                    &terminals,
                    invalid,
                    terminal_limits(),
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    for bounded in [
        ConductionTerminalLimits {
            maximum_terminals: 2,
            ..terminal_limits()
        },
        ConductionTerminalLimits {
            maximum_workspace_values: 1,
            ..terminal_limits()
        },
    ] {
        assert!(
            problem
                .solve_terminals(&terminals, drive, bounded, controls(), &mut || false)
                .is_err()
        );
    }
    for invalid in [
        [&faces[0][..], &faces[2][..]], // no shared middle electrode joins the bodies
        [&faces[0][..], &faces[0][..]],
        [&[][..], &faces[2][..]],
        [&[usize::MAX][..], &faces[2][..]],
    ] {
        assert!(
            problem
                .solve_terminals(
                    &invalid,
                    ConductionTerminalControl::Voltages(&[0.0, 5.0]),
                    terminal_limits(),
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let (after, after_report) = problem
        .solve_terminals(
            &terminals,
            drive,
            terminal_limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(before_report, after_report);
    assert_eq!(before.voltage_v, after.voltage_v);
    assert_eq!(before.cell_power_w, after.cell_power_w);
    assert_eq!(before.current_density_a_m2, after.current_density_a_m2);
}
