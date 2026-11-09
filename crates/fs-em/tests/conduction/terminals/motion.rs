//! G1/G3/G4: total moving current, floating transfer and independent work.
use super::*;
use fs_em::conduction::motion::ConductionMotion;

#[test]
fn moving_grouped_electrodes_recover_current_floating_voltage_and_signed_work() {
    let area = 0.4;
    let original = PlanarConduction::new(
        disconnected_layers(Metric2::planar(area).unwrap(), 2),
        conductors(),
        limits(),
        &mut || false,
    )
    .unwrap();
    for sign in [-1.0, 1.0] {
        for (j, speed) in [(0.0, 2.0), (1.0, -1.0), (1.0, 1.0), (1.0, 5.0)] {
            let (j, speed) = (sign * j, sign * speed);
            let velocity = [[speed, 0.0]; 4];
            let motion = ConductionMotion {
                velocity_m_s: &velocity,
                normal_flux_density_t: &[1.0; 4],
            };
            // In each unit length E'_y=-rho*J, E_y=E'_y+u_x*B.
            // Therefore V=[0, J-uB, 4J-2uB], including open-circuit voltage.
            let expected_v = [0.0, j - speed, 4.0 * j - 2.0 * speed];
            let currents = [-area * j, 0.0, area * j];
            for bias in [0.0, 2.0_f64.powi(54), 1e250] {
                for control in [
                    ConductionTerminalControl::Currents(&currents),
                    ConductionTerminalControl::Mixed {
                        current_terminals: &[2, 1],
                        currents_a: &[area * j, 0.0],
                        voltages_v: &[bias],
                    },
                ] {
                    let (field, report) = original
                        .solve_terminals_with_motion(
                            &TERMINALS,
                            control,
                            terminal_limits(),
                            motion,
                            controls(),
                            &mut || false,
                        )
                        .unwrap();
                    for i in 0..3 {
                        near(report.voltage_v[i], expected_v[i]);
                        near(report.current_a[i], currents[i]);
                    }
                    for cell in 0..4 {
                        let rho = if cell < 2 { 1.0 } else { 3.0 };
                        near(field.conductor_field_v_m()[cell][1], -rho * j);
                        near(field.electric_field_v_m()[cell][1], speed - rho * j);
                        near(field.current_density_a_m2()[cell][1], -j);
                        near(field.lorentz_force_density_n_m3()[cell][0], -j);
                        near(field.cell_joule_power_w()[cell], 0.5 * area * rho * j * j);
                    }
                    near(field.joule_power_w(), 4.0 * area * j * j);
                    near(field.mechanical_power_w(), -2.0 * area * speed * j);
                    near(report.supplied_power_w, area * j * expected_v[2]);
                    near(field.terminal_power_w(), report.supplied_power_w);
                    near(report.power_residual_w, 0.0);
                    near(field.power_balance_residual_w(), 0.0);
                    let reaction = |vertices: &[usize]| {
                        vertices
                            .iter()
                            .map(|&v| field.nodal_boundary_current_a()[v].unwrap())
                            .sum::<f64>()
                    };
                    near(reaction(&[2, 3]), area * j);
                    near(reaction(&[4, 5]), -area * j);
                    near(reaction(TERMINALS[1]), 0.0);
                }
            }
            let (field, report) = original
                .solve_terminals_with_motion(
                    &TERMINALS,
                    ConductionTerminalControl::Voltages(&expected_v),
                    terminal_limits(),
                    motion,
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let boundary: Vec<_> = TERMINALS
                .iter()
                .zip(expected_v)
                .flat_map(|(vertices, v)| vertices.iter().map(move |&vertex| (vertex, v)))
                .collect();
            let traced = original
                .solve_with_motion(&boundary, &[], motion, controls(), &mut || false)
                .unwrap();
            for (&a, &b) in field
                .cell_joule_power_w()
                .iter()
                .zip(traced.cell_joule_power_w())
            {
                near(a, b);
            }
            near(report.supplied_power_w, traced.terminal_power_w());
            if j != 0.0 && speed * j > 4.0 * j * j {
                assert!(report.supplied_power_w < 0.0);
                assert!(field.joule_power_w() > 0.0);
            }
        }
    }
    assert_eq!(original.conductors(), conductors());
}

#[test]
fn moving_terminal_budgets_and_cancellation_leave_original_sources_reusable() {
    let original = PlanarConduction::new(
        disconnected_layers(Metric2::planar(1.0).unwrap(), 2),
        conductors(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let current = [-1.0, 0.0, 1.0];
    let motion = ConductionMotion {
        velocity_m_s: &[[2.0, 0.0]; 4],
        normal_flux_density_t: &[1.0; 4],
    };
    let mut calls = 0;
    let (baseline, _) = original
        .solve_terminals_with_motion(
            &TERMINALS,
            ConductionTerminalControl::Currents(&current),
            terminal_limits(),
            motion,
            controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, calls / 2, calls] {
        let mut polls = 0;
        assert!(matches!(
            original.solve_terminals_with_motion(
                &TERMINALS,
                ConductionTerminalControl::Currents(&current),
                terminal_limits(),
                motion,
                controls(),
                &mut || {
                    polls += 1;
                    polls == cutoff
                },
            ),
            Err(ConductionError::Cancelled)
        ));
    }
    for bound in [
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
            original
                .solve_terminals_with_motion(
                    &TERMINALS,
                    ConductionTerminalControl::Currents(&current),
                    bound,
                    motion,
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    assert!(
        original
            .solve_terminals_with_motion(
                &TERMINALS,
                ConductionTerminalControl::Currents(&[-1.0, 0.0, 2.0]),
                terminal_limits(),
                motion,
                controls(),
                &mut || false,
            )
            .is_err()
    );
    let (retry, _) = original
        .solve_terminals_with_motion(
            &TERMINALS,
            ConductionTerminalControl::Currents(&current),
            terminal_limits(),
            motion,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        retry.current_density_a_m2(),
        baseline.current_density_a_m2()
    );
    assert_eq!(retry.cell_joule_power_w(), baseline.cell_joule_power_w());
    assert_eq!(original.conductors(), conductors());
}
