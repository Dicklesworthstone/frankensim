//! G1/G3: exact layered Gauss-law sources, not empirical material validation.
use super::{controls, dielectric, near, slab};
use fs_em::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, ElectrostaticError,
    PlanarElectrostaticSources,
};

#[test]
fn charged_interface_obeys_displacement_jump_and_electrode_gauss_balance() {
    let problem = slab(dielectric(2e-11), dielectric(5e-11));
    let sigma = 8e-11;
    let sheets = [ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: sigma,
    }];
    let field = problem
        .solve_with_sources(
            PlanarElectrostaticSources {
                charge_sheets: &sheets,
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[(0, 0.0), (3, 0.0), (2, 0.0), (5, 0.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let interface_v = sigma / (2e-11 / 0.5 + 5e-11 / 0.5);
    near(field.relative_voltage_v()[1], interface_v);
    near(field.relative_voltage_v()[4], interface_v);
    let left_e = -interface_v / 0.5;
    let right_e = interface_v / 0.5;
    for cell in 0..4 {
        near(
            field.electric_field_v_m()[cell][0],
            if cell < 2 { left_e } else { right_e },
        );
        assert!(field.electric_field_v_m()[cell][1].abs() < 1e-12);
    }
    near(
        field.displacement_c_m2()[2][0] - field.displacement_c_m2()[0][0],
        sigma,
    );
    near(field.sheet_charge_c(), sigma * 0.2);
    near(field.boundary_charge_c(), -sigma * 0.2);
    near(field.energy_j(), 0.5 * sigma * 0.2 * interface_v);
    assert_eq!(field.volume_charge_c(), 0.0);
    assert!(field.charge_balance_residual_c().abs() < 1e-25);
    assert!(field.energy_balance_residual_j().abs() < 1e-25);
}

#[test]
fn impressed_polarization_with_insulated_exit_has_depolarizing_field() {
    let epsilon = 2e-11;
    let p = 6e-11;
    let problem = slab(dielectric(epsilon), dielectric(epsilon));
    let polarization = [[p, 0.0]; 4];
    let field = problem
        .solve_with_sources(
            PlanarElectrostaticSources {
                impressed_polarization_c_m2: Some(&polarization),
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[(0, 0.0), (3, 0.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    for cell in 0..4 {
        near(field.electric_field_v_m()[cell][0], -p / epsilon);
        assert!(field.displacement_c_m2()[cell][0].abs() < 1e-25);
        assert!(field.displacement_c_m2()[cell][1].abs() < 1e-25);
    }
    near(field.voltage_v()[2], p / epsilon);
    near(field.energy_j(), 0.5 * p * p / epsilon * 0.2);
    near(field.polarization_pairing_j(), -2.0 * field.energy_j());
    assert!(field.boundary_charge_c().abs() < 1e-25);
    assert_eq!(field.impressed_polarization_c_m2().unwrap(), &polarization);
    // D=0 does not mean zero field energy: 1/2 E.D would be incorrect here.
    assert!(field.energy_j() > 0.0);
}

#[test]
fn shorted_polarization_supplies_opposite_free_electrode_charges() {
    let problem = slab(dielectric(2e-11), dielectric(2e-11));
    let p = 6e-11;
    let polarization = [[p, 0.0]; 4];
    let field = problem
        .solve_with_sources(
            PlanarElectrostaticSources {
                impressed_polarization_c_m2: Some(&polarization),
                ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
            },
            &[(0, 9.0), (3, 9.0), (2, 9.0), (5, 9.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(
        field
            .electric_field_v_m()
            .iter()
            .flatten()
            .all(|e| e.abs() < 1e-14)
    );
    let q = field.nodal_boundary_charge_c();
    near(q[0].unwrap() + q[3].unwrap(), p * 0.2);
    near(q[2].unwrap() + q[5].unwrap(), -p * 0.2);
    assert!(field.energy_j() < 1e-35);
    for d in field.displacement_c_m2() {
        near(d[0], p);
    }
}

#[test]
fn mixed_sheet_polarization_and_neumann_sources_keep_gauge_and_physical_units() {
    let problem = slab(dielectric(2e-11), dielectric(5e-11));
    let polarization = [[7e-11, 0.0], [7e-11, 0.0], [-3e-11, 0.0], [-3e-11, 0.0]];
    let sheet = ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: 8e-11,
    };
    let trace = ElectrostaticDisplacementTrace {
        vertices: [2, 5],
        outward_displacement_c_m2: 1.4e-10,
    };
    let run = |reference,
               sheets: &[ElectrostaticChargeSheet],
               traces: &[ElectrostaticDisplacementTrace]| {
        problem
            .solve_with_sources(
                PlanarElectrostaticSources {
                    impressed_polarization_c_m2: Some(&polarization),
                    charge_sheets: sheets,
                    displacement_traces: traces,
                    ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
                },
                &[(0, reference), (3, reference)],
                controls(),
                &mut || false,
            )
            .unwrap()
    };
    let field = run(0.0, &[sheet], &[trace]);
    let shifted = run(
        1e15,
        &[ElectrostaticChargeSheet {
            vertices: [4, 1],
            ..sheet
        }],
        &[ElectrostaticDisplacementTrace {
            vertices: [5, 2],
            ..trace
        }],
    );
    assert_eq!(field.electric_field_v_m(), shifted.electric_field_v_m());
    assert_eq!(field.displacement_c_m2(), shifted.displacement_c_m2());
    assert_eq!(field.energy_j(), shifted.energy_j());
    assert_eq!(field.free_charge_sheets(), shifted.free_charge_sheets());
    assert_eq!(field.displacement_traces(), shifted.displacement_traces());
    let e1 =
        (trace.outward_displacement_c_m2 - sheet.free_charge_density_c_m2 - polarization[0][0])
            / 2e-11;
    let e2 = (trace.outward_displacement_c_m2 - polarization[2][0]) / 5e-11;
    near(field.electric_field_v_m()[0][0], e1);
    near(field.electric_field_v_m()[2][0], e2);
    near(field.voltage_v()[2], -0.5 * (e1 + e2));
    near(
        field.boundary_charge_c(),
        (trace.outward_displacement_c_m2 - sheet.free_charge_density_c_m2) * 0.2,
    );
    near(
        field.outward_displacement_c(),
        trace.outward_displacement_c_m2 * 0.2,
    );
    near(
        field.energy_j(),
        0.5 * 0.1 * (2e-11 * e1 * e1 + 5e-11 * e2 * e2),
    );
    assert!(field.charge_balance_residual_c().abs() < 1e-25);
    assert!(field.energy_balance_residual_j().abs() < 1e-25);
}

#[test]
fn invalid_interface_and_conflicting_boundary_sources_refuse_without_publication() {
    let problem = slab(dielectric(2e-11), dielectric(2e-11));
    let valid_sheet = ElectrostaticChargeSheet {
        vertices: [1, 4],
        free_charge_density_c_m2: 1e-11,
    };
    for sheets in [
        vec![ElectrostaticChargeSheet {
            vertices: [0, 3],
            ..valid_sheet
        }],
        vec![ElectrostaticChargeSheet {
            vertices: [1, usize::MAX],
            ..valid_sheet
        }],
        vec![
            valid_sheet,
            ElectrostaticChargeSheet {
                vertices: [4, 1],
                ..valid_sheet
            },
        ],
        vec![ElectrostaticChargeSheet {
            free_charge_density_c_m2: f64::NAN,
            ..valid_sheet
        }],
    ] {
        assert!(
            problem
                .solve_with_sources(
                    PlanarElectrostaticSources {
                        charge_sheets: &sheets,
                        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
                    },
                    &[(0, 0.0), (3, 0.0)],
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    for vertices in [[0, 3], [1, 4], [2, usize::MAX]] {
        assert!(
            problem
                .solve_with_sources(
                    PlanarElectrostaticSources {
                        displacement_traces: &[ElectrostaticDisplacementTrace {
                            vertices,
                            outward_displacement_c_m2: 1e-11
                        }],
                        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
                    },
                    &[(0, 0.0), (3, 0.0)],
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    for polarization in [vec![[0.0; 2]; 3], vec![[f64::INFINITY, 0.0]; 4]] {
        assert!(
            problem
                .solve_with_sources(
                    PlanarElectrostaticSources {
                        impressed_polarization_c_m2: Some(&polarization),
                        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
                    },
                    &[(0, 0.0), (3, 0.0)],
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
}

#[test]
fn cancellation_during_source_work_preserves_original_problem_for_retry() {
    let problem = slab(dielectric(2e-11), dielectric(2e-11));
    let polarization = [[6e-11, 0.0]; 4];
    let sources = PlanarElectrostaticSources {
        impressed_polarization_c_m2: Some(&polarization),
        ..PlanarElectrostaticSources::volume_charge(&[0.0; 4])
    };
    let boundary = [(0, 0.0), (3, 0.0)];
    let mut total_polls = 0;
    let accepted = problem
        .solve_with_sources(sources, &boundary, controls(), &mut || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, total_polls / 2, total_polls] {
        let mut count = 0;
        assert!(matches!(
            problem.solve_with_sources(sources, &boundary, controls(), &mut || {
                count += 1;
                count >= stop
            }),
            Err(ElectrostaticError::Cancelled)
        ));
    }
    let retry = problem
        .solve_with_sources(sources, &boundary, controls(), &mut || false)
        .unwrap();
    assert_eq!(accepted.electric_field_v_m(), retry.electric_field_v_m());
    assert_eq!(
        accepted.nodal_boundary_charge_c(),
        retry.nodal_boundary_charge_c()
    );
}
