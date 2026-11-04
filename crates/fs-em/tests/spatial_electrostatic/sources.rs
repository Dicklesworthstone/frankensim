//! G1/G3/G4: free interface charge and fixed polarization on the actual 3D
//! dielectric owner. Analytic slab solutions independently fix all load signs.

use super::*;
use fs_em::spatial::{
    electrostatic::{ElectrostaticChargeSheet, SpatialElectrostaticSources},
    force::SpatialMaxwellForceLimits,
};
use fs_ga::Se3;

fn plane_sheets(
    problem: &SpatialElectrostatics,
    x: f64,
    density: f64,
) -> Vec<ElectrostaticChargeSheet> {
    problem
        .mesh()
        .faces
        .iter()
        .enumerate()
        .filter_map(|(face, vertices)| {
            vertices
                .iter()
                .all(|&v| problem.positions_m()[v as usize][0] == x)
                .then_some(ElectrostaticChargeSheet {
                    face,
                    free_charge_density_c_m2: density,
                })
        })
        .collect()
}

#[test]
fn charged_polarized_interface_matches_gauss_energy_and_voltage_gauge() {
    let width = 0.25;
    let epsilon = [2e-11, 5e-11];
    let p = [4e-11, -3e-11];
    let sigma = 6e-11;
    let voltage = 2.5;
    let (positions, tets) = cube(&[0.0, width, 1.0], 1, 1);
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        (0..12).map(|cell| dielectric(epsilon[cell / 6])).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let polarization: Vec<_> = (0..12).map(|cell| [p[cell / 6], 0.0, 0.0]).collect();
    let mut sheets = plane_sheets(&problem, width, sigma);
    assert_eq!(sheets.len(), 2);
    let rhs = -voltage + width * p[0] / epsilon[0] + (1.0 - width) * (p[1] - sigma) / epsilon[1];
    let left_d = rhs / (width / epsilon[0] + (1.0 - width) / epsilon[1]);
    let d = [left_d, left_d + sigma];
    let e = [(d[0] - p[0]) / epsilon[0], (d[1] - p[1]) / epsilon[1]];
    let expected_energy =
        0.5 * (width * epsilon[0] * e[0] * e[0] + (1.0 - width) * epsilon[1] * e[1] * e[1]);
    let expected_polarization_pairing = width * e[0] * p[0] + (1.0 - width) * e[1] * p[1];
    let mut previous = None;
    for bias in [0.0, 2_f64.powi(45)] {
        let field = problem
            .solve_with_sources(
                SpatialElectrostaticSources {
                    free_charge_density_c_m3: &[0.0; 12],
                    impressed_polarization_c_m2: Some(&polarization),
                    sheets: &sheets,
                },
                &end_voltages(&problem, bias, bias + voltage),
                controls(),
                &mut || false,
            )
            .unwrap();
        near(field.sheet_charge_c(), sigma);
        near(field.boundary_charge_c(), -sigma);
        assert_eq!(field.volume_charge_c(), 0.0);
        near(field.energy_j(), expected_energy);
        near(
            field.polarization_pairing_j(),
            expected_polarization_pairing,
        );
        for cell in 0..12 {
            near_vector(field.electric_field_v_m()[cell], [e[cell / 6], 0.0, 0.0]);
            near_vector(field.displacement_c_m2()[cell], [d[cell / 6], 0.0, 0.0]);
        }
        for (v, position) in problem.positions_m().iter().enumerate() {
            let expected = -e[0] * position[0].min(width) - e[1] * (position[0] - width).max(0.0);
            if expected == 0.0 {
                assert!(field.relative_voltage_v()[v].abs() < 1e-14);
            } else {
                near(field.relative_voltage_v()[v], expected);
            }
        }
        let terminals = terminal_faces(&problem);
        for (terminal, expected) in terminals.iter().zip([-d[1], d[0]]) {
            let vertices: std::collections::BTreeSet<_> = terminal
                .iter()
                .flat_map(|&face| problem.mesh().faces[face])
                .collect();
            let reaction: f64 = vertices
                .iter()
                .map(|&v| field.nodal_boundary_charge_c()[v as usize].unwrap())
                .sum();
            near(reaction, expected);
        }
        for sheet in &sheets {
            let jump = field
                .interface_displacement_jumps_c()
                .iter()
                .find(|(face, _)| *face == sheet.face)
                .unwrap()
                .1;
            near(jump, -sigma * problem.face_area_m2(sheet.face).unwrap());
        }
        assert!(
            field
                .interface_displacement_residuals_c()
                .iter()
                .all(|(_, residual)| residual.abs() < 1e-23)
        );
        assert!(field.charge_balance_residual_c().abs() < 1e-23);
        assert!(field.energy_balance_residual_j().abs() < 1e-23);
        assert!(field.hodge_energy_residual_j().abs() < 1e-23);
        assert_eq!(
            field.impressed_polarization_c_m2(),
            Some(polarization.as_slice())
        );
        if let Some((relative, fields, energy, reactions)) = previous {
            assert_eq!(field.relative_voltage_v(), relative);
            assert_eq!(field.electric_field_v_m(), fields);
            assert_eq!(field.energy_j(), energy);
            assert_eq!(field.nodal_boundary_charge_c(), reactions);
        }
        previous = Some((
            field.relative_voltage_v().to_vec(),
            field.electric_field_v_m().to_vec(),
            field.energy_j(),
            field.nodal_boundary_charge_c().to_vec(),
        ));
        // Canonical source assembly must not depend on authored sheet order.
        sheets.reverse();
    }
}

#[test]
fn open_circuit_polarization_retains_depolarizing_storage_and_tensor_rotation() {
    let s = 0.5_f64.sqrt();
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    let principal = [2e-11, 5e-11, 7e-11];
    let properties = std::array::from_fn::<_, 3, _>(|i| {
        (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[i], principal[i])
    });
    let material = source(&properties);
    let electric = [2.0, -1.0, 0.5];
    let epsilon_e = [8.5e-11, -6.5e-11, 3.5e-11];
    let polarization = epsilon_e.map(|v| -v);
    let expected_energy = 0.5 * (0..3).map(|j| electric[j] * epsilon_e[j]).sum::<f64>();
    let rotate = |p: [f64; 3]| [0.6 * p[0] - 0.8 * p[2], p[1], 0.8 * p[0] + 0.6 * p[2]];
    for rotated in [false, true] {
        let (positions, mut tets) = cube(&[0.0, 0.5, 1.0], 2, 2);
        if rotated {
            // Negative stored orientation must not change physical cell volume.
            for tet in &mut tets {
                tet.swap(1, 2);
            }
        }
        let cells = tets.len();
        let law = SpatialDielectricLaw::orthotropic(if rotated { axes.map(rotate) } else { axes })
            .unwrap();
        let problem = SpatialElectrostatics::new_with_laws(
            if rotated {
                positions.into_iter().map(rotate).collect()
            } else {
                positions
            },
            tets,
            vec![material.clone(); cells],
            vec![law; cells],
            limits(),
            &mut || false,
        )
        .unwrap();
        let polarization = vec![
            if rotated {
                rotate(polarization)
            } else {
                polarization
            };
            cells
        ];
        let field = problem
            .solve_with_sources(
                SpatialElectrostaticSources {
                    free_charge_density_c_m3: &vec![0.0; cells],
                    impressed_polarization_c_m2: Some(&polarization),
                    sheets: &[],
                },
                &SpatialElectrostaticBoundary::default(),
                controls(),
                &mut || false,
            )
            .unwrap();
        near(field.energy_j(), expected_energy);
        near(field.polarization_pairing_j(), -2.0 * expected_energy);
        assert_eq!(field.gauge_vertices().len(), 1);
        assert_eq!(field.boundary_charge_c(), 0.0);
        assert!(field.nodal_boundary_charge_c().iter().all(Option::is_none));
        for cell in 0..cells {
            near_vector(
                field.electric_field_v_m()[cell],
                if rotated { rotate(electric) } else { electric },
            );
            assert!(
                field.displacement_c_m2()[cell]
                    .iter()
                    .all(|d| d.abs() < 1e-21)
            );
        }
        // E.D/2 would falsely report a zero store although E is nonzero.
        let wrong_energy: f64 = (0..cells)
            .map(|cell| {
                0.5 * problem.cell_volume_m3(cell).unwrap()
                    * (0..3)
                        .map(|j| {
                            field.electric_field_v_m()[cell][j] * field.displacement_c_m2()[cell][j]
                        })
                        .sum::<f64>()
            })
            .sum();
        assert!(wrong_energy.abs() < 1e-21);
    }
}

#[test]
fn sheet_charge_uses_the_same_pure_neumann_and_mixed_boundary_compatibility() {
    let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
    let epsilon = 2e-11;
    let sigma = 4e-11;
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(epsilon); 12],
        limits(),
        &mut || false,
    )
    .unwrap();
    let sheets = plane_sheets(&problem, 0.5, sigma);
    let polarization = [[epsilon, 0.0, 0.0]; 12];
    let sources = SpatialElectrostaticSources {
        free_charge_density_c_m3: &[0.0; 12],
        impressed_polarization_c_m2: Some(&polarization),
        sheets: &sheets,
    };
    let terminals = terminal_faces(&problem);
    let mut boundary = SpatialElectrostaticBoundary {
        voltage: vec![],
        displacement: terminals[0]
            .iter()
            .map(|&face| ElectrostaticDisplacementFace {
                face,
                outward_displacement_c_m2: sigma,
            })
            .collect(),
    };
    for grounded in [false, true] {
        let field = problem
            .solve_with_sources(sources, &boundary, controls(), &mut || false)
            .unwrap();
        assert_eq!(field.gauge_vertices().len(), usize::from(!grounded));
        near(field.sheet_charge_c(), sigma);
        near(field.boundary_charge_c(), -sigma);
        near(field.energy_j(), 0.5 * epsilon);
        for cell in 0..12 {
            near_vector(
                field.electric_field_v_m()[cell],
                [if cell < 6 { -1.0 } else { 1.0 }, 0.0, 0.0],
            );
        }
        assert!(field.charge_balance_residual_c().abs() < 1e-23);
        boundary.voltage = terminals[1]
            .iter()
            .map(|&face| ElectrostaticVoltageFace {
                face,
                voltage_v: [0.0; 3],
            })
            .collect();
    }
    assert!(
        problem
            .solve_with_sources(
                sources,
                &SpatialElectrostaticBoundary::default(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn invalid_source_support_and_late_cancellation_leave_the_owner_reusable() {
    let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(2e-11); 12],
        limits(),
        &mut || false,
    )
    .unwrap();
    let sheets = plane_sheets(&problem, 0.5, 4e-11);
    let polarization = [[2e-11, 0.0, 0.0]; 12];
    let sources = SpatialElectrostaticSources {
        free_charge_density_c_m3: &[0.0; 12],
        impressed_polarization_c_m2: Some(&polarization),
        sheets: &sheets,
    };
    let boundary = end_voltages(&problem, 0.0, 1.0);
    for invalid in [
        vec![sheets[0], sheets[0]],
        vec![ElectrostaticChargeSheet {
            face: problem.exterior_faces().next().unwrap(),
            ..sheets[0]
        }],
        vec![ElectrostaticChargeSheet {
            face: usize::MAX,
            ..sheets[0]
        }],
        vec![ElectrostaticChargeSheet {
            free_charge_density_c_m2: f64::NAN,
            ..sheets[0]
        }],
    ] {
        assert!(
            problem
                .solve_with_sources(
                    SpatialElectrostaticSources {
                        sheets: &invalid,
                        ..sources
                    },
                    &boundary,
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    for invalid in [&[[0.0; 3]; 11][..], &[[f64::INFINITY, 0.0, 0.0]; 12][..]] {
        assert!(
            problem
                .solve_with_sources(
                    SpatialElectrostaticSources {
                        impressed_polarization_c_m2: Some(invalid),
                        ..sources
                    },
                    &boundary,
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let mut polls = 0;
    let original = problem
        .solve_with_sources(sources, &boundary, controls(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, 8, polls - 1, polls] {
        let mut current = 0;
        assert!(matches!(
            problem.solve_with_sources(sources, &boundary, controls(), &mut || {
                current += 1;
                current == cancel_at
            }),
            Err(SpatialElectrostaticError::Cancelled)
        ));
    }
    let replay = problem
        .solve_with_sources(sources, &boundary, controls(), &mut || false)
        .unwrap();
    assert_eq!(original.relative_voltage_v(), replay.relative_voltage_v());
    assert_eq!(original.energy_j(), replay.energy_j());
}

#[test]
fn maxwell_stress_refuses_polarization_and_charged_sheets_in_its_collar() {
    let (positions, tets) = cube(&[0.0, 0.25, 0.5, 0.75, 1.0], 4, 4);
    let cells = tets.len();
    let problem = SpatialElectrostatics::new(
        positions,
        tets,
        vec![dielectric(2e-11); cells],
        limits(),
        &mut || false,
    )
    .unwrap();
    let selected: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(cell, tet)| {
            tet.iter()
                .all(|&v| {
                    problem.positions_m()[v as usize]
                        .iter()
                        .all(|x| (0.25..=0.75).contains(x))
                })
                .then_some(cell)
        })
        .collect();
    assert_eq!(selected.len(), 48);
    let rho = vec![0.0; cells];
    let boundary = SpatialElectrostaticBoundary::default();
    let force_limits = SpatialMaxwellForceLimits {
        max_cells: cells,
        max_faces: problem.mesh().faces.len(),
    };
    let empty = problem
        .solve(&rho, &boundary, controls(), &mut || false)
        .unwrap();
    empty
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits,
            &mut || false,
        )
        .unwrap();
    let polarization = vec![[2e-11, 0.0, 0.0]; cells];
    let polarized = problem
        .solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: &rho,
                impressed_polarization_c_m2: Some(&polarization),
                sheets: &[],
            },
            &boundary,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        polarized
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits,
                &mut || false
            )
            .is_err()
    );
    let sigma = 4e-11;
    let sheets = plane_sheets(&problem, 0.25, sigma);
    let charged = problem
        .solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: &rho,
                impressed_polarization_c_m2: None,
                sheets: &sheets,
            },
            &end_voltages(&problem, 0.0, 1.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        charged
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits,
                &mut || false
            )
            .is_err()
    );
}
