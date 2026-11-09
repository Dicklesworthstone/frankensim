//! G1/G3/G4: direct grouped terminals retain the affine motional source.
//! Independent bar balances distinguish mechanical conversion from Joule heat.
use super::*;
use fs_em::conduction::terminals::{ConductionTerminalControl, ConductionTerminalLimits};

fn direct_limits() -> ConductionTerminalLimits {
    ConductionTerminalLimits {
        maximum_terminals: 3,
        maximum_workspace_values: 16_384,
    }
}

#[test]
fn direct_moving_bar_keeps_open_short_generating_and_motoring_states() {
    let problem = scalar_problem(&[0.0, 0.4, 1.0], 5.0);
    let faces = [0.0, 1.0].map(|x| electrode_faces(&problem, x));
    let terminals = faces.each_ref().map(Vec::as_slice);
    let motion = SpatialConductionMotion {
        velocity_m_s: &[[0.0, 0.8, 0.0]; 12],
        flux_density_t: &[[0.0, 0.0, 1.5]; 12],
    };
    let matrix = problem
        .conductance_matrix(&terminals, terminal_controls(), controls(), &mut || false)
        .unwrap();
    for current in [0.0_f64, -6.0, -2.0, 2.0] {
        // I enters the right electrode. EMF=1.2 V, R=0.2 ohm.
        let voltage = 1.2 + current / 5.0;
        let currents = [-current, current];
        let voltages = [0.0, voltage];
        let reversed_current = [-current];
        let reference = matrix
            .field_at_currents_with_motion(&currents, motion, controls(), &mut || false)
            .unwrap();
        for control in [
            ConductionTerminalControl::Currents(&currents),
            ConductionTerminalControl::Voltages(&voltages),
            ConductionTerminalControl::Mixed {
                current_terminals: &[0],
                currents_a: &reversed_current,
                voltages_v: &[1e250],
            },
        ] {
            let (field, report) = problem
                .solve_terminals_with_motion(
                    &terminals,
                    control,
                    direct_limits(),
                    motion,
                    controls(),
                    &mut || false,
                )
                .unwrap();
            assert!(field.solve_report().is_some());
            near(report.voltage_v[0], 0.0);
            near(report.voltage_v[1], voltage);
            near(report.current_a[0], -current);
            near(report.current_a[1], current);
            near(report.supplied_power_w, voltage * current);
            near(report.power_residual_w, 0.0);
            assert!(report.maximum_current_residual_a <= controls().current_tolerance_a);
            near(field.joule_power_w(), current * current / 5.0);
            near(field.mechanical_power_w(), 1.2 * current);
            near(field.electrical_power_w(), voltage * current);
            near(field.terminal_power_w(), report.supplied_power_w);
            near(field.power_balance_residual_w(), 0.0);
            for cell in 0..12 {
                near(field.electric_field_v_m()[cell][0], -voltage);
                near(field.conductor_field_v_m()[cell][0], -current / 5.0);
                near(field.current_density_a_m2()[cell][0], -current);
                near(field.lorentz_force_density_n_m3()[cell][1], 1.5 * current);
                near(
                    field.cell_joule_power_w()[cell],
                    reference.cell_joule_power_w()[cell],
                );
                near(
                    field.cell_mechanical_power_w()[cell],
                    reference.cell_mechanical_power_w()[cell],
                );
            }
            if current == 0.0 {
                assert!(report.voltage_v[1] > 1.0);
                near(field.joule_power_w(), 0.0);
            } else if current == -2.0 {
                assert!(report.supplied_power_w < 0.0);
                assert!(field.joule_power_w() > 0.0);
            } else if current == -6.0 {
                near(report.supplied_power_w, 0.0);
                assert!(field.joule_power_w() > 0.0);
                assert!(field.mechanical_power_w() < 0.0);
            } else {
                assert!(report.supplied_power_w > 0.0);
                assert!(field.mechanical_power_w() > 0.0);
            }
        }
    }
}

#[test]
fn direct_floating_electrode_carries_motional_current_between_disconnected_bodies() {
    let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
    let (other, cells) = slab(&[1.0, 1.5, 2.0]);
    let offset = positions.len() as u32;
    positions.extend(other);
    tets.extend(cells.into_iter().map(|t| t.map(|v| v + offset)));
    let sources = [contacts::conductor(2.0), contacts::conductor(3.0)];
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..24)
            .map(|i| sources[usize::from(i >= 12)].clone())
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = [0.0, 1.0, 2.0].map(|x| electrode_faces(&problem, x));
    let terminals = faces.each_ref().map(Vec::as_slice);
    let velocity: Vec<_> = (0..24)
        .map(|i| [0.0, if i < 12 { 1.2 } else { -0.3 }, 0.0])
        .collect();
    let motion = SpatialConductionMotion {
        velocity_m_s: &velocity,
        flux_density_t: &[[0.0, 0.0, 1.0]; 24],
    };
    let matrix = problem
        .conductance_matrix(&terminals, terminal_controls(), controls(), &mut || false)
        .unwrap();
    let mixed = matrix
        .prepare_current_control(&[1], mixed_controls(), &mut || false)
        .unwrap();
    for injection in [0.0, 0.5] {
        let middle = (injection + 2.0 / 3.0 + 1.2 / 2.0 + 0.3 / 3.0) / (0.5 + 1.0 / 3.0);
        let ja = (1.2 - middle) / 2.0;
        let jb = (middle - 2.0 - 0.3) / 3.0;
        let (field, report) = problem
            .solve_terminals_with_motion(
                &terminals,
                ConductionTerminalControl::Mixed {
                    current_terminals: &[1],
                    currents_a: &[injection],
                    voltages_v: &[0.0, 2.0],
                },
                direct_limits(),
                motion,
                controls(),
                &mut || false,
            )
            .unwrap();
        let reference = mixed
            .field_at_controls_with_motion(
                &[0.0, 2.0],
                &[injection],
                motion,
                controls(),
                &mut || false,
            )
            .unwrap();
        for (i, voltage) in [0.0, middle, 2.0].into_iter().enumerate() {
            near(report.voltage_v[i], voltage);
        }
        for (i, current) in [ja, injection, -jb].into_iter().enumerate() {
            near(report.current_a[i], current);
            near(electrode_current(&problem, &field, &faces[i]), current);
        }
        near(
            field.nodal_boundary_current_a()[8..12]
                .iter()
                .map(|v| v.unwrap())
                .sum(),
            -ja,
        );
        near(
            field.nodal_boundary_current_a()[12..16]
                .iter()
                .map(|v| v.unwrap())
                .sum(),
            jb,
        );
        near(field.joule_power_w(), 2.0 * ja * ja + 3.0 * jb * jb);
        near(field.mechanical_power_w(), -1.2 * ja + 0.3 * jb);
        near(report.supplied_power_w, middle * injection - 2.0 * jb);
        near(report.power_residual_w, 0.0);
        for cell in 0..24 {
            near(
                field.current_density_a_m2()[cell][0],
                if cell < 12 { ja } else { jb },
            );
            near(
                field.cell_joule_power_w()[cell],
                reference.cell_joule_power_w()[cell],
            );
            near(
                field.cell_mechanical_power_w()[cell],
                reference.cell_mechanical_power_w()[cell],
            );
        }
        if injection == 0.0 {
            assert!(ja.abs() > 0.1 && jb.abs() > 0.1);
        }
        let (_, reordered) = problem
            .solve_terminals_with_motion(
                &[terminals[1], terminals[2], terminals[0]],
                ConductionTerminalControl::Currents(&[injection, -jb, ja]),
                direct_limits(),
                motion,
                controls(),
                &mut || false,
            )
            .unwrap();
        for (i, voltage) in [0.0, 2.0 - middle, -middle].into_iter().enumerate() {
            near(reordered.voltage_v[i], voltage);
        }
        near(reordered.supplied_power_w, report.supplied_power_w);
    }
}

#[test]
fn direct_motion_keeps_stationary_contact_heat_and_transactional_refusals() {
    let source = contacts::conductor(1.0);
    let interface = contacts::interface(&source, &source, 1e5, 1.0 / 3.0);
    let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
    let (other, cells) = slab(&[1.0, 2.0]);
    positions.extend(other);
    tets.extend(cells.into_iter().map(|t| t.map(|v| v + 12)));
    let problem = SpatialConduction::new_with_contacts(
        positions,
        tets,
        vec![source; 18],
        vec![
            SpatialConductionContact {
                first_vertices: [8, 9, 11],
                second_vertices: [12, 13, 15],
                source: interface.clone(),
            },
            SpatialConductionContact {
                first_vertices: [8, 10, 11],
                second_vertices: [12, 14, 15],
                source: interface,
            },
        ],
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = [0.0, 2.0].map(|x| electrode_faces(&problem, x));
    let terminals = faces.each_ref().map(Vec::as_slice);
    let mut velocity = [[0.0; 3]; 18];
    velocity[..6].fill([0.0, 0.8, 0.0]);
    let motion = SpatialConductionMotion {
        velocity_m_s: &velocity,
        flux_density_t: &[[0.0, 0.0, 1.0]; 18],
    };
    let drive = ConductionTerminalControl::Currents(&[-1.0, 1.0]);
    let mut polls = 0;
    let (before, before_report) = problem
        .solve_terminals_with_motion(
            &terminals,
            drive,
            direct_limits(),
            motion,
            controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    // The moving half metre gives 0.4 V EMF. Bulk R=2, contact R=1.
    near(before_report.voltage_v[1], 3.4);
    near(before_report.supplied_power_w, 3.4);
    near(before.contact_power_w(), 1.0);
    near(before.cell_joule_power_w().iter().sum(), 2.0);
    near(before.joule_power_w(), 3.0);
    near(before.mechanical_power_w(), 0.4);
    near(before.power_balance_residual_w(), 0.0);
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            problem.solve_terminals_with_motion(
                &terminals,
                drive,
                direct_limits(),
                motion,
                controls(),
                &mut || {
                    seen += 1;
                    seen == stop
                },
            ),
            Err(SpatialConductionError::Cancelled)
        ));
    }
    assert!(
        problem
            .solve_terminals_with_motion(
                &terminals,
                drive,
                ConductionTerminalLimits {
                    maximum_workspace_values: 1,
                    ..direct_limits()
                },
                motion,
                controls(),
                &mut || false,
            )
            .is_err()
    );
    let mut moving_contact = velocity;
    moving_contact[problem.contact_cells()[0][0]] = [0.0, 0.8, 0.0];
    for invalid in [
        SpatialConductionMotion {
            velocity_m_s: &velocity[..17],
            ..motion
        },
        SpatialConductionMotion {
            velocity_m_s: &moving_contact,
            ..motion
        },
        SpatialConductionMotion {
            flux_density_t: &[[0.0, f64::NAN, 1.0]; 18],
            ..motion
        },
    ] {
        assert!(
            problem
                .solve_terminals_with_motion(
                    &terminals,
                    drive,
                    direct_limits(),
                    invalid,
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    let (after, after_report) = problem
        .solve_terminals_with_motion(
            &terminals,
            drive,
            direct_limits(),
            motion,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(before_report, after_report);
    assert_eq!(before.voltage_v(), after.voltage_v());
    assert_eq!(before.current_density_a_m2(), after.current_density_a_m2());
    assert_eq!(before.cell_joule_power_w(), after.cell_joule_power_w());
    assert_eq!(
        before.cell_mechanical_power_w(),
        after.cell_mechanical_power_w()
    );
}
