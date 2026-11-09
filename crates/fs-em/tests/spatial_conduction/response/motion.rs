//! Nonlinear original-source closure under prescribed three-dimensional motion.
use super::direct_terminals::terminal_limits;
use super::*;
use fs_em::spatial::conduction::motion::SpatialConductionMotion;
use fs_em::spatial::conduction::response::motion::SpatialMovingConductionResponseDrive;

fn specimen(coordinate: ConductorMagnitude) -> SpatialConductionResponse {
    let source = law_with_basis(
        Some(coordinate),
        ElectricalPropertyBasis::Resistivity,
        1.0,
        0.5,
        10.0,
    );
    let (mut positions, tets) = slab(&[0.0, 1.0]);
    rotate(&mut positions);
    SpatialConductionResponse::new(positions, tets, vec![source; 6], limits(), &mut || false)
        .unwrap()
        .with_relaxation(0.5)
        .unwrap()
}

fn voltage_boundary(faces: &[Vec<usize>], voltage: f64) -> SpatialConductionBoundary {
    SpatialConductionBoundary {
        voltage: faces
            .iter()
            .zip([0.0, voltage])
            .flat_map(|(faces, value)| {
                faces.iter().map(move |&face| ConductionVoltageFace {
                    face,
                    voltage_v: [value; 3],
                })
            })
            .collect(),
        current: Vec::new(),
    }
}

fn current_from_field(coordinate: ConductorMagnitude, field: f64) -> f64 {
    // Independent uniform-bar algebra: rho=1+0.5*|E'| or
    // (1+0.5*|J|)*J=E'. The latter uses its stable quadratic root.
    match coordinate {
        ConductorMagnitude::ElectricField => field / (1.0 + 0.5 * field.abs()),
        ConductorMagnitude::CurrentDensity => {
            2.0 * field / (1.0 + (1.0 + 2.0 * field.abs()).sqrt())
        }
    }
}

#[test]
fn original_field_and_current_laws_close_short_generation_motor_and_open_states() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let model = specimen(coordinate);
        let faces = face_sets(&model, &[0.0, 1.0]);
        let terminals = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
        for sign in [-1.0, 1.0] {
            for voltage in [0.0, 1.0, 2.0, 3.0] {
                let voltage = sign * voltage;
                let motional = 2.0 * sign;
                let conductor_e = motional - voltage;
                let j = current_from_field(coordinate, conductor_e);
                let boundary = voltage_boundary(&faces, voltage);
                let voltages = [0.0, voltage];
                let terminal_drive = ConductionTerminalDrive {
                    terminals: &terminals,
                    control: ConductionTerminalControl::Voltages(&voltages),
                    limits: terminal_limits(),
                };
                for direct in [false, true] {
                    let drive = if direct {
                        SpatialMovingConductionResponseDrive::Terminals(terminal_drive)
                    } else {
                        SpatialMovingConductionResponseDrive::Boundary(&boundary)
                    };
                    let candidate = model
                        .solve_with_motion_candidate(
                            drive,
                            SpatialConductionMotion {
                                velocity_m_s: &[[1.0, 0.0, 0.0]; 6],
                                flux_density_t: &[[0.0, motional, 0.0]; 6],
                            },
                            limits(),
                            controls(),
                            response_controls(),
                            &mut || false,
                        )
                        .unwrap();
                    let report = candidate.report();
                    let field = &report.field;
                    for cell in 0..6 {
                        near(field.electric_field_v_m()[cell][2], -voltage);
                        near(field.conductor_field_v_m()[cell][2], conductor_e);
                        near(field.current_density_a_m2()[cell][2], j);
                        near(field.current_density_a_m2()[cell][0], 0.0);
                        near(field.current_density_a_m2()[cell][1], 0.0);
                        near(field.lorentz_force_density_n_m3()[cell][0], -motional * j);
                        near(field.cell_joule_power_w()[cell], conductor_e * j / 6.0);
                        near(field.cell_mechanical_power_w()[cell], -motional * j / 6.0);
                        let sample = &report.evaluated[cell];
                        let coordinate_value = match coordinate {
                            ConductorMagnitude::ElectricField => conductor_e.abs(),
                            ConductorMagnitude::CurrentDensity => j.abs(),
                        };
                        near(
                            sample.conductor.material().query().axes()["operating"],
                            coordinate_value,
                        );
                        assert!(sample.current_density_defect_a_m2.iter().all(|v| {
                            v.abs() <= response_controls().current_density_tolerance_a_m2
                        }));
                        assert!(sample.power_defect_w.abs() <= controls().power_tolerance_w);
                    }
                    near(field.joule_power_w(), conductor_e * j);
                    near(field.mechanical_power_w(), -motional * j);
                    near(field.terminal_power_w(), -voltage * j);
                    near(
                        field.joule_power_w() + field.mechanical_power_w(),
                        field.terminal_power_w(),
                    );
                    near(field.power_balance_residual_w(), 0.0);
                    near(field.work_pairing_residual_w(), 0.0);
                    assert_eq!(report.field_solves, report.iterations);
                    assert!(report.source_evaluations >= 6 * report.iterations);
                    assert!(report.field_solves <= response_controls().maximum_field_solves);
                    assert_eq!(report.terminals.is_some(), direct);
                    if let Some(terminal) = &report.terminals {
                        near(terminal.voltage_v[1], voltage);
                        near(terminal.current_a[0], j);
                        near(terminal.current_a[1], -j);
                        near(terminal.supplied_power_w, -voltage * j);
                        near(terminal.power_residual_w, 0.0);
                    }
                    if voltage.abs() == 1.0 {
                        assert!(field.terminal_power_w() < 0.0);
                        assert!(field.joule_power_w() > 0.0);
                    } else if voltage.abs() == 3.0 {
                        assert!(field.mechanical_power_w() > 0.0);
                    } else if voltage.abs() == 2.0 {
                        near(field.joule_power_w(), 0.0);
                        near(field.mechanical_power_w(), 0.0);
                    }
                }
            }
        }
        // Every candidate is speculative; the input source still has its
        // originally authored field coordinate.
        assert_eq!(
            model.conductors()[0].material().query().axes()["operating"],
            0.0
        );
    }
}

#[test]
fn direct_current_and_floating_electrodes_keep_the_motional_voltage_offset() {
    let model = specimen(ConductorMagnitude::CurrentDensity);
    let faces = face_sets(&model, &[0.0, 1.0]);
    let terminals = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
    for j in [-0.5_f64, 0.0, 0.5] {
        let currents = [j, -j];
        let current_terminal = [1];
        let right_current = [-j];
        for mixed in [false, true] {
            let drive = ConductionTerminalDrive {
                terminals: &terminals,
                control: if mixed {
                    ConductionTerminalControl::Mixed {
                        current_terminals: &current_terminal,
                        currents_a: &right_current,
                        voltages_v: &[0.0],
                    }
                } else {
                    ConductionTerminalControl::Currents(&currents)
                },
                limits: terminal_limits(),
            };
            let candidate = model
                .solve_with_motion_candidate(
                    drive,
                    SpatialConductionMotion {
                        velocity_m_s: &[[1.0, 0.0, 0.0]; 6],
                        flux_density_t: &[[0.0, 2.0, 0.0]; 6],
                    },
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || false,
                )
                .unwrap();
            let report = candidate.report();
            let rho = 1.0 + 0.5 * j.abs();
            let voltage = 2.0 - rho * j;
            let terminal = report.terminals.as_ref().unwrap();
            near(terminal.voltage_v[1], voltage);
            near(terminal.current_a[0], j);
            near(terminal.current_a[1], -j);
            near(report.field.joule_power_w(), rho * j * j);
            near(report.field.mechanical_power_w(), -2.0 * j);
            near(report.field.terminal_power_w(), -voltage * j);
            for sample in &report.evaluated {
                near(
                    sample.conductor.material().query().axes()["operating"],
                    j.abs(),
                );
            }
        }
    }
}

#[test]
fn zero_motion_parity_cancellation_and_global_trial_budget_preserve_the_owner() {
    let model = specimen(ConductorMagnitude::CurrentDensity);
    let original = model.conductors().to_vec();
    let boundary = voltage_boundary(&face_sets(&model, &[0.0, 1.0]), 1.0);
    let stationary = model
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    let motion = SpatialConductionMotion {
        velocity_m_s: &[[0.0; 3]; 6],
        flux_density_t: &[[0.0, 2.0, 0.0]; 6],
    };
    let mut polls = 0;
    let moved = model
        .solve_with_motion_candidate(
            &boundary,
            motion,
            limits(),
            controls(),
            response_controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(
        moved.report().field.current_density_a_m2(),
        stationary.report().field.current_density_a_m2,
    );
    assert_eq!(
        moved.report().field.cell_joule_power_w(),
        stationary.report().field.cell_power_w,
    );
    assert_eq!(
        moved.report().field_solves,
        stationary.report().field_solves
    );
    near(moved.report().field.mechanical_power_w(), 0.0);
    for cutoff in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            model.solve_with_motion_candidate(
                &boundary,
                motion,
                limits(),
                controls(),
                response_controls(),
                &mut || {
                    count += 1;
                    count == cutoff
                },
            ),
            Err(SpatialConductionResponseError::Cancelled),
        ));
    }
    assert!(matches!(
        model.solve_with_motion_candidate(
            &boundary,
            motion,
            limits(),
            controls(),
            SpatialConductionResponseControls {
                maximum_field_solves: 1,
                ..response_controls()
            },
            &mut || false,
        ),
        Err(SpatialConductionResponseError::FieldSolveBudget {
            completed: 1,
            required: 1,
            limit: 1,
        }),
    ));
    assert_eq!(model.conductors(), original);
    let retry = model
        .solve_with_motion_candidate(
            &boundary,
            motion,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.report().field.voltage_v(),
        moved.report().field.voltage_v()
    );
}

#[test]
fn original_material_support_is_checked_at_conductor_field_not_laboratory_field() {
    let model = specimen(ConductorMagnitude::ElectricField);
    let boundary = voltage_boundary(&face_sets(&model, &[0.0, 1.0]), 0.0);
    // Laboratory E=0 is supported, but the actual constitutive E'=11 is not.
    assert!(matches!(
        model.solve_with_motion_candidate(
            &boundary,
            SpatialConductionMotion {
                velocity_m_s: &[[1.0, 0.0, 0.0]; 6],
                flux_density_t: &[[0.0, 11.0, 0.0]; 6],
            },
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        ),
        Err(SpatialConductionResponseError::Material(_)),
    ));
    assert_eq!(
        model.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}

#[test]
fn stationary_contacts_keep_heat_and_refuse_motion_even_without_motional_emf() {
    let bulk = law_with_basis(
        Some(ConductorMagnitude::CurrentDensity),
        ElectricalPropertyBasis::Resistivity,
        1.0,
        0.5,
        10.0,
    );
    let source = interface(bulk.original(), bulk.original(), 1e5, 0.5);
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (other, extra) = slab(&[1.0, 2.0]);
    positions.extend(other);
    tets.extend(extra.into_iter().map(|tet| tet.map(|vertex| vertex + 8)));
    rotate(&mut positions);
    let model = SpatialConductionResponse::new_with_contacts(
        positions,
        tets,
        vec![bulk; 12],
        contact_faces(source),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap();
    let boundary = voltage_boundary(&face_sets(&model, &[0.0, 2.0]), 2.0);
    let stationary = model
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    let motion = SpatialConductionMotion {
        velocity_m_s: &[[0.0; 3]; 12],
        flux_density_t: &[[1.0, 0.0, 0.0]; 12],
    };
    let moved = model
        .solve_with_motion_candidate(
            &boundary,
            motion,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    near(
        moved.report().field.joule_power_w(),
        stationary.report().field.power_w,
    );
    near(
        moved.report().field.contact_power_w(),
        stationary
            .report()
            .field
            .contacts
            .iter()
            .map(|c| c.power_w)
            .sum(),
    );
    assert!(moved.report().field.contact_power_w() > 0.0);
    // Parallel u and B give zero motional EMF, but the interface is still
    // physically moving and cannot borrow the stationary surface law.
    assert!(matches!(
        model.solve_with_motion_candidate(
            &boundary,
            SpatialConductionMotion {
                velocity_m_s: &[[1.0, 0.0, 0.0]; 12],
                ..motion
            },
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        ),
        Err(SpatialConductionResponseError::Field(
            SpatialConductionError::Invalid(
                "moving electrical contacts require an owned surface work and heat law"
            )
        )),
    ));
}
