//! Independent nonlinear free-sheet, fixed-polarization and outward-flux oracles.
use super::*;
use fs_em::{
    electrostatic::nonlinear::spatial::SpatialNonlinearElectricField,
    spatial::electrostatic::{
        ElectrostaticChargeSheet, ElectrostaticDisplacementFace, SpatialElectrostaticSources,
    },
};

fn nonlinear_box() -> (Vec<[f64; 3]>, SpatialNonlinearElectrostatics) {
    let (positions, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 2, 2);
    let materials = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells,
        materials,
        limits(),
        &mut || false,
    )
    .unwrap();
    (positions, problem)
}

fn sheets_at(
    problem: &SpatialNonlinearElectrostatics,
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

fn displacement_faces(faces: &[usize], value: f64) -> Vec<ElectrostaticDisplacementFace> {
    faces
        .iter()
        .map(|&face| ElectrostaticDisplacementFace {
            face,
            outward_displacement_c_m2: value,
        })
        .collect()
}

fn refuses_undeclared_source_motion(field: &SpatialNonlinearElectricField<'_>) {
    assert!(field.has_additional_sources());
    assert!(matches!(
        field.shape_forces(&mut || false),
        Err(NonlinearElectrostaticError::Electric(
            ElectrostaticError::Invalid(_)
        ))
    ));
    assert!(matches!(
        field.virtual_work(
            &vec![[0.0; 3]; field.problem().mesh().vertex_count],
            &mut || false
        ),
        Err(NonlinearElectrostaticError::Electric(
            ElectrostaticError::Invalid(_)
        ))
    ));
}

#[test]
fn nonlinear_internal_sheet_retains_the_displacement_jump_storage_and_fixed_source_tangent() {
    let (positions, problem) = nonlinear_box();
    let terminals = [
        faces_at(&problem, &positions, 0.0),
        faces_at(&problem, &positions, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let sheets = sheets_at(&problem, 0.5 * GAP, 0.5 * D0);
    assert!(!sheets.is_empty());
    let sources = SpatialElectrostaticSources {
        free_charge_density_c_m3: &rho,
        impressed_polarization_c_m2: None,
        sheets: &sheets,
    };
    let initial: Vec<_> = positions.iter().map(|p| 2.0 * E0 * p[0]).collect();
    let field = electrodes
        .solve_with_sources(
            sources,
            &[],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(2.0 * V0),
            ],
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.terminal_charge_c()[0], -2.75 * Q0, Q0);
    near(field.terminal_charge_c()[1], 2.25 * Q0, Q0);
    near(field.sheet_charge_c(), 0.5 * Q0, Q0);
    near(field.energy_j(), 1.8125 * U0, U0);
    near(field.coenergy_j(), 3.3125 * U0, U0);
    assert_eq!(field.charge_sheets().len(), sheets.len());
    assert!(field.impressed_polarization_c_m2().is_none());
    for (v, p) in positions.iter().enumerate() {
        if p[0] == 0.5 * GAP {
            near(field.relative_voltage_v()[v], 1.25 * V0, V0);
        }
    }
    for (cell, vertices) in problem.mesh().tets.iter().enumerate() {
        let left = vertices
            .iter()
            .any(|&v| positions[v as usize][0] < 0.5 * GAP);
        let (e, d) = if left {
            (-2.5 * E0, -2.75 * D0)
        } else {
            (-1.5 * E0, -2.25 * D0)
        };
        near_vector(field.electric_field_v_m()[cell], [e, 0.0, 0.0], E0);
        near_vector(field.displacement_c_m2()[cell], [d, 0.0, 0.0], D0);
        near_vector(
            field.constitutive_displacement_c_m2()[cell],
            [d, 0.0, 0.0],
            D0,
        );
    }
    // The internal sheet's potential pairing is essential to the original
    // nonlinear U+W* identity, although it is not electrode charge.
    near(
        field.energy_j() + field.coenergy_j(),
        1.25 * V0 * 0.5 * Q0 + 2.0 * V0 * 2.25 * Q0,
        U0,
    );
    assert!(field.energy_identity_residual_j().abs() <= controls().energy_tolerance_j);
    let differential = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    for (&actual, expected) in differential.matrix_f().iter().zip([0.5, -0.5, -0.5, 0.5]) {
        near(actual, expected * C0, C0);
    }
    let charged = electrodes
        .solve_with_sources(
            sources,
            &[],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(2.25 * Q0),
            ],
            field.relative_voltage_v(),
            controls(),
            &mut || false,
        )
        .unwrap();
    near(charged.terminal_voltage_v()[1], 2.0 * V0, V0);
    near(charged.energy_j(), field.energy_j(), U0);
    refuses_undeclared_source_motion(&field);
}

#[test]
fn fixed_impressed_polarization_has_zero_open_circuit_total_displacement_and_signed_source_work() {
    let (positions, problem) = nonlinear_box();
    let terminals = [
        faces_at(&problem, &positions, 0.0),
        faces_at(&problem, &positions, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let polarization = vec![[2.5 * D0, 0.0, 0.0]; rho.len()];
    let initial: Vec<_> = positions.iter().map(|p| 2.0 * E0 * p[0]).collect();
    let field = electrodes
        .solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: &rho,
                impressed_polarization_c_m2: Some(&polarization),
                sheets: &[],
            },
            &[],
            &[
                ElectrodeControl::Voltage(1e250),
                ElectrodeControl::Charge(0.0),
            ],
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(field.reference_voltage_for_terminal_v(1), Some(1e250));
    near(field.terminal_voltage_v()[1], 2.0 * V0, V0);
    near(field.terminal_charge_c()[0], 0.0, Q0);
    near(field.terminal_charge_c()[1], 0.0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    near(field.polarization_pairing_j(), -5.0 * U0, U0);
    near(
        field.energy_j() + field.coenergy_j() + field.polarization_pairing_j(),
        0.0,
        U0,
    );
    assert_eq!(
        field.impressed_polarization_c_m2(),
        Some(polarization.as_slice())
    );
    for cell in 0..rho.len() {
        near_vector(field.electric_field_v_m()[cell], [-2.0 * E0, 0.0, 0.0], E0);
        near_vector(
            field.constitutive_displacement_c_m2()[cell],
            [-2.5 * D0, 0.0, 0.0],
            D0,
        );
        near_vector(field.displacement_c_m2()[cell], [0.0; 3], D0);
    }
    assert!(field.energy_identity_residual_j().abs() <= controls().energy_tolerance_j);
    assert!(field.polarization_work_residual_j().abs() <= controls().energy_tolerance_j);
    refuses_undeclared_source_motion(&field);
}

#[test]
fn prescribed_outward_displacement_retains_its_physical_sign_and_electrode_charge() {
    let (positions, problem) = nonlinear_box();
    let left = faces_at(&problem, &positions, 0.0);
    let right = faces_at(&problem, &positions, GAP);
    let electrodes = problem.prepare_electrodes(&[&left], &mut || false).unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let flux = displacement_faces(&right, 2.5 * D0);
    let initial: Vec<_> = positions.iter().map(|p| -2.0 * E0 * p[0]).collect();
    let field = electrodes
        .solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: &rho,
                impressed_polarization_c_m2: None,
                sheets: &[],
            },
            &flux,
            &[ElectrodeControl::Voltage(0.0)],
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.terminal_charge_c()[0], 2.5 * Q0, Q0);
    near(field.outward_displacement_c(), 2.5 * Q0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    assert_eq!(field.displacement_faces().len(), flux.len());
    for (v, p) in positions.iter().enumerate() {
        near(field.relative_voltage_v()[v], -2.0 * E0 * p[0], V0);
    }
    for cell in 0..rho.len() {
        near_vector(field.electric_field_v_m()[cell], [2.0 * E0, 0.0, 0.0], E0);
        near_vector(field.displacement_c_m2()[cell], [2.5 * D0, 0.0, 0.0], D0);
    }
    near(
        field.terminal_charge_c()[0] - field.outward_displacement_c(),
        0.0,
        Q0,
    );
    assert!(field.charge_balance_residual_c().abs() <= controls().charge_tolerance_c);
    assert!(field.energy_identity_residual_j().abs() <= controls().energy_tolerance_j);
    refuses_undeclared_source_motion(&field);
}

#[test]
fn combined_tensor_volume_polarization_sheet_and_flux_sources_match_the_native_hodge() {
    let (positions, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], GAP, GAP, 2, 2);
    let names = PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES;
    let material = scalar(&[
        (names[0], 2.0 * EPSILON),
        (names[1], 5.0 * EPSILON),
        (names[2], 7.0 * EPSILON),
    ]);
    let s = 0.5_f64.sqrt();
    let law =
        SpatialDielectricLaw::orthotropic([[s, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, s]]).unwrap();
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells.clone(),
        vec![
            SpatialDielectricMaterial::Linear {
                material: material.clone(),
                law
            };
            cells.len()
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let native = SpatialElectrostatics::new_with_laws(
        positions.clone(),
        cells.clone(),
        vec![material; cells.len()],
        vec![law; cells.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let left = faces_at(&problem, &positions, 0.0);
    let right = faces_at(&problem, &positions, GAP);
    let electrodes = problem.prepare_electrodes(&[&left], &mut || false).unwrap();
    let rho = vec![0.3 * D0 / GAP; cells.len()];
    let polarization: Vec<_> = cells
        .iter()
        .map(|vertices| {
            if vertices
                .iter()
                .any(|&v| positions[v as usize][0] < 0.5 * GAP)
            {
                [0.4 * D0, 0.2 * D0, -0.1 * D0]
            } else {
                [-0.3 * D0, 0.1 * D0, 0.2 * D0]
            }
        })
        .collect();
    let sheets = sheets_at(&problem, 0.5 * GAP, 0.2 * D0);
    let flux = displacement_faces(&right, 0.7 * D0);
    let sources = SpatialElectrostaticSources {
        free_charge_density_c_m3: &rho,
        impressed_polarization_c_m2: Some(&polarization),
        sheets: &sheets,
    };
    let field = electrodes
        .solve_with_sources(
            sources,
            &flux,
            &[ElectrodeControl::Voltage(0.0)],
            &vec![0.0; positions.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    let exact = native
        .solve_with_sources(
            sources,
            &SpatialElectrostaticBoundary {
                voltage: left
                    .iter()
                    .map(|&face| ElectrostaticVoltageFace {
                        face,
                        voltage_v: [0.0; 3],
                    })
                    .collect(),
                displacement: flux,
            },
            ElectrostaticControls {
                linear: linear_controls(),
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    for cell in 0..cells.len() {
        near_vector(
            field.electric_field_v_m()[cell],
            exact.electric_field_v_m()[cell],
            E0,
        );
        near_vector(
            field.displacement_c_m2()[cell],
            exact.displacement_c_m2()[cell],
            D0,
        );
        let total: [f64; 3] = std::array::from_fn(|i| {
            field.constitutive_displacement_c_m2()[cell][i] + polarization[cell][i]
        });
        near_vector(total, field.displacement_c_m2()[cell], D0);
    }
    let energy_scale = E0 * D0 * GAP.powi(3);
    let charge_scale = D0 * GAP * GAP;
    near(field.energy_j(), exact.energy_j(), energy_scale);
    near(field.coenergy_j(), exact.energy_j(), energy_scale);
    near(
        field.polarization_pairing_j(),
        exact.polarization_pairing_j(),
        energy_scale,
    );
    near(field.volume_charge_c(), 0.3 * charge_scale, charge_scale);
    near(field.sheet_charge_c(), 0.2 * charge_scale, charge_scale);
    near(
        field.outward_displacement_c(),
        0.7 * charge_scale,
        charge_scale,
    );
    near(
        field.terminal_charge_c()[0],
        0.2 * charge_scale,
        charge_scale,
    );
    let ground_charge: f64 = electrodes.terminal_vertices()[0]
        .iter()
        .map(|&v| exact.nodal_boundary_charge_c()[v].unwrap())
        .sum();
    near(field.terminal_charge_c()[0], ground_charge, charge_scale);
    assert!(field.maximum_component_charge_residual_c() <= controls().charge_tolerance_c);
    assert!(field.polarization_work_residual_j().abs() <= controls().energy_tolerance_j);
}

#[test]
fn source_trace_conflicts_and_invalid_faces_refuse_while_cancelled_and_zero_sources_remain_reusable()
 {
    let (positions, problem) = nonlinear_box();
    let left = faces_at(&problem, &positions, 0.0);
    let right = faces_at(&problem, &positions, GAP);
    let middle = sheets_at(&problem, 0.5 * GAP, 0.5 * D0);
    let rho = vec![0.0; problem.mesh().tets.len()];
    let zero = vec![0.0; positions.len()];
    let no_sources = SpatialElectrostaticSources {
        free_charge_density_c_m3: &rho,
        impressed_polarization_c_m2: None,
        sheets: &[],
    };
    let metals = problem
        .prepare_electrodes(&[&left, &right], &mut || false)
        .unwrap();
    let conflicting = displacement_faces(&right, 2.5 * D0);
    for control in [
        ElectrodeControl::Voltage(2.0 * V0),
        ElectrodeControl::Charge(2.5 * Q0),
    ] {
        assert!(
            metals
                .solve_with_sources(
                    no_sources,
                    &conflicting,
                    &[ElectrodeControl::Voltage(0.0), control],
                    &zero,
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let grounded = problem.prepare_electrodes(&[&left], &mut || false).unwrap();
    for bad_sheets in [
        vec![ElectrostaticChargeSheet {
            face: right[0],
            free_charge_density_c_m2: 0.5 * D0,
        }],
        vec![middle[0], middle[0]],
    ] {
        assert!(
            grounded
                .solve_with_sources(
                    SpatialElectrostaticSources {
                        sheets: &bad_sheets,
                        ..no_sources
                    },
                    &[],
                    &[ElectrodeControl::Voltage(0.0)],
                    &zero,
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    for bad_flux in [
        displacement_faces(&[middle[0].face], D0),
        displacement_faces(&[right[0], right[0]], D0),
    ] {
        assert!(
            grounded
                .solve_with_sources(
                    no_sources,
                    &bad_flux,
                    &[ElectrodeControl::Voltage(0.0)],
                    &zero,
                    controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let polarization = vec![[2.5 * D0, 0.0, 0.0]; rho.len()];
    let polarized = SpatialElectrostaticSources {
        impressed_polarization_c_m2: Some(&polarization),
        ..no_sources
    };
    assert!(matches!(
        metals.solve_with_sources(
            polarized,
            &[],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(0.0)
            ],
            &zero,
            controls(),
            &mut || true
        ),
        Err(NonlinearElectrostaticError::Electric(
            ElectrostaticError::Cancelled
        ))
    ));
    let retried = metals
        .solve_with_sources(
            polarized,
            &[],
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(0.0),
            ],
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    near(retried.terminal_voltage_v()[1], 2.0 * V0, V0);
    // Explicitly represented zero sources retain geometry validation but do
    // not invent a nonzero transported-source restriction on shape work.
    let zero_p = vec![[0.0; 3]; rho.len()];
    let zero_sheets = sheets_at(&problem, 0.5 * GAP, 0.0);
    let zero_flux = displacement_faces(&right, 0.0);
    let unforced = grounded
        .solve_with_sources(
            SpatialElectrostaticSources {
                free_charge_density_c_m3: &rho,
                impressed_polarization_c_m2: Some(&zero_p),
                sheets: &zero_sheets,
            },
            &zero_flux,
            &[ElectrodeControl::Voltage(0.0)],
            &zero,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(!unforced.has_additional_sources());
    assert!(
        unforced
            .shape_forces(&mut || false)
            .unwrap()
            .nodal_force_n()
            .iter()
            .flatten()
            .all(|&f| f == 0.0)
    );
    assert_eq!(
        unforced
            .virtual_work(&vec![[0.0; 3]; positions.len()], &mut || false)
            .unwrap(),
        0.0
    );
}
