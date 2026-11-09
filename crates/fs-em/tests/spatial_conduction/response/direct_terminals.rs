//! Direct native terminal controls close source laws at the physical operating state.
use super::*;

pub(super) fn terminal_limits() -> ConductionTerminalLimits {
    ConductionTerminalLimits {
        maximum_terminals: 8,
        maximum_workspace_values: 20_000,
    }
}

#[test]
fn direct_voltage_and_conserved_current_controls_stay_inside_actual_source_support() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        // Two physical bodies share one complete middle electrode. That
        // electrode may transfer current between faces while its total is zero.
        let source = law(Some(coordinate), 0.2, 0.05);
        let (mut positions, mut tets) = slab(&[0.0, 1.0]);
        let (other, other_tets) = slab(&[1.0, 2.0]);
        let offset = positions.len() as u32;
        positions.extend(other);
        tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
        rotate(&mut positions);
        let problem = SpatialConductionResponse::new(
            positions,
            tets,
            vec![source; 12],
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = face_sets(&problem, &[0.0, 1.0, 2.0]);
        let terminals = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let operating = 0.02;
        let sigma = 1.0 + 0.2 * operating;
        let voltage_controlled = matches!(coordinate, ConductorMagnitude::ElectricField);
        let (e, j) = if voltage_controlled {
            (operating, sigma * operating)
        } else {
            (operating / sigma, operating)
        };
        for sign in [-1.0, 1.0] {
            let voltages = [0.0, sign * e, sign * 2.0 * e];
            let currents = [-sign * j, 0.0, sign * j];
            let control = if voltage_controlled {
                ConductionTerminalControl::Voltages(&voltages)
            } else {
                ConductionTerminalControl::Currents(&currents)
            };
            let drive = ConductionTerminalDrive {
                terminals: &terminals,
                control,
                limits: terminal_limits(),
            };
            let candidate = problem
                .solve_candidate(
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || false,
                )
                .unwrap();
            let report = candidate.report();
            let terminal = report.terminals.as_ref().unwrap();
            assert!(report.iterations > 1);
            assert_eq!(report.field_solves, report.iterations);
            assert_eq!(report.source_evaluations, 12 * report.iterations);
            for ((&v, &i), (&expected_v, &expected_i)) in terminal
                .voltage_v
                .iter()
                .zip(&terminal.current_a)
                .zip(voltages.iter().zip(&currents))
            {
                near(v, expected_v);
                near(i, expected_i);
            }
            for sample in &report.evaluated {
                near(
                    sample.conductor.material().query().axes()["operating"],
                    operating,
                );
            }
            for (field_e, field_j) in report
                .field
                .electric_field_v_m
                .iter()
                .zip(&report.field.current_density_a_m2)
            {
                near(field_e[2], -sign * e);
                near(field_j[2], -sign * j);
            }
            near(report.field.power_w, 2.0 * e * j);
            near(terminal.supplied_power_w, 2.0 * e * j);
            assert!(terminal.maximum_current_residual_a <= controls().current_tolerance_a);
            assert!(terminal.power_residual_w.abs() <= controls().power_tolerance_w);

            let bounded = SpatialConductionResponseControls {
                maximum_field_solves: 1,
                ..response_controls()
            };
            assert!(matches!(
                problem.solve_candidate(drive, limits(), controls(), bounded, &mut || false),
                Err(SpatialConductionResponseError::FieldSolveBudget {
                    completed: 1,
                    required: 1,
                    limit: 1,
                })
            ));
            assert!(
                problem
                    .solve_candidate(
                        ConductionTerminalDrive {
                            limits: ConductionTerminalLimits {
                                maximum_workspace_values: 0,
                                ..terminal_limits()
                            },
                            ..drive
                        },
                        limits(),
                        controls(),
                        response_controls(),
                        &mut || false,
                    )
                    .is_err()
            );
        }

        // The bound is on actual E/J, not an artificial unit probe. An actual
        // load outside that support must still refuse without changing history.
        let invalid_voltages = [0.0, 0.1, 0.2];
        let invalid_currents = [-0.1, 0.0, 0.1];
        let invalid = ConductionTerminalDrive {
            terminals: &terminals,
            control: if voltage_controlled {
                ConductionTerminalControl::Voltages(&invalid_voltages)
            } else {
                ConductionTerminalControl::Currents(&invalid_currents)
            },
            limits: terminal_limits(),
        };
        assert!(matches!(
            problem.solve_candidate(
                invalid,
                limits(),
                controls(),
                response_controls(),
                &mut || false,
            ),
            Err(SpatialConductionResponseError::Material(_))
        ));
        assert_eq!(
            problem.conductors()[0].material().query().axes()["operating"],
            0.0
        );
    }
}
