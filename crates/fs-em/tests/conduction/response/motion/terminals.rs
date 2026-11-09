//! G1/G3/G4: original nonlinear moving laws under actual electrode controls.
use super::*;
use fs_em::conduction::response::terminals::{
    ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits,
};

fn terminal_limits() -> ConductionTerminalLimits {
    ConductionTerminalLimits {
        maximum_terminals: 4,
        maximum_workspace_values: 8192,
    }
}

#[test]
fn nonlinear_moving_current_and_mixed_electrodes_use_actual_material_frame() {
    for coordinate in [
        ConductorMagnitude::ElectricField,
        ConductorMagnitude::CurrentDensity,
    ] {
        let basis = if coordinate == ConductorMagnitude::ElectricField {
            ElectricalPropertyBasis::Conductivity
        } else {
            ElectricalPropertyBasis::Resistivity
        };
        let law = response(basis, coordinate, 1.0, 0.25, 10.0);
        let original = specimen(law.clone(), law, 0.2)
            .with_relaxation(0.5)
            .unwrap();
        let motion = ConductionMotion {
            velocity_m_s: &[[0.0, 2.0]; 4],
            normal_flux_density_t: &[1.0; 4],
        };
        for j in [-0.8_f64, 0.0, 0.8] {
            let ep = if coordinate == ConductorMagnitude::ElectricField {
                2.0 * j / (1.0 + (1.0 + j.abs()).sqrt())
            } else {
                (1.0 + 0.25 * j.abs()) * j
            };
            for control in [
                ConductionTerminalControl::Currents(&[0.2 * j, -0.2 * j]),
                ConductionTerminalControl::Mixed {
                    current_terminals: &[1],
                    currents_a: &[-0.2 * j],
                    voltages_v: &[2.0_f64.powi(54)],
                },
            ] {
                let candidate = original
                    .solve_terminals_with_motion_candidate(
                        ConductionTerminalDrive {
                            terminals: &[&[0, 3], &[2, 5]],
                            control,
                            limits: terminal_limits(),
                        },
                        motion,
                        limits(),
                        controls(),
                        policy(),
                        256,
                        &mut || false,
                    )
                    .unwrap();
                let report = candidate.report();
                let terminals = report.terminals.as_ref().unwrap();
                near(terminals.voltage_v[1], 2.0 - ep);
                near(terminals.current_a[0], 0.2 * j);
                near(terminals.current_a[1], -0.2 * j);
                near(terminals.supplied_power_w, -0.2 * j * (2.0 - ep));
                near(report.field.power_w, 0.2 * ep * j);
                near(report.motion.mechanical_power_w, -0.4 * j);
                near(report.motion.terminal_power_w, terminals.supplied_power_w);
                near(terminals.power_residual_w, 0.0);
                for cell in 0..4 {
                    near(report.field.current_density_a_m2[cell][0], j);
                    near(report.motion.conductor_field_v_m[cell][0], ep);
                    near(
                        report.field.evaluated[cell]
                            .conductor
                            .material()
                            .query()
                            .axes()["operating"],
                        if coordinate == ConductorMagnitude::ElectricField {
                            ep.abs()
                        } else {
                            j.abs()
                        },
                    );
                }
                if j > 0.0 {
                    assert!(terminals.supplied_power_w < 0.0);
                }
                if j != 0.0 {
                    assert!(report.field.field_solves > 1);
                }
            }
        }
        let drive = ConductionTerminalDrive {
            terminals: &[&[0, 3], &[2, 5]],
            control: ConductionTerminalControl::Currents(&[0.16, -0.16]),
            limits: terminal_limits(),
        };
        assert!(matches!(
            original.solve_terminals_with_motion_candidate(
                drive,
                motion,
                limits(),
                controls(),
                policy(),
                1,
                &mut || false,
            ),
            Err(ConductionResponseError::FieldBudget {
                completed: 1,
                required: 1,
                maximum: 1
            })
        ));
        assert!(matches!(
            original.solve_terminals_with_motion_candidate(
                drive,
                motion,
                limits(),
                controls(),
                policy(),
                0,
                &mut || false,
            ),
            Err(ConductionResponseError::FieldBudget {
                completed: 0,
                required: 1,
                maximum: 0
            })
        ));
        assert!(matches!(
            original.solve_terminals_with_motion_candidate(
                drive,
                motion,
                limits(),
                controls(),
                policy(),
                256,
                &mut || true,
            ),
            Err(ConductionResponseError::Cancelled)
        ));
        assert!(
            original
                .conductors()
                .iter()
                .all(|c| c.material().query().axes()["operating"] == 0.0)
        );
    }
}

#[test]
fn moving_direct_electrodes_keep_stationary_joint_heat_and_velocity_guard() {
    let model = joint_specimen(0.2);
    let velocity = [
        [0.0, 2.0],
        [0.0, 2.0],
        [0.0; 2],
        [0.0; 2],
        [0.0; 2],
        [0.0; 2],
    ];
    let drive = ConductionTerminalDrive {
        terminals: &[&[0, 3], &[7, 8]],
        control: ConductionTerminalControl::Currents(&[0.04, -0.04]),
        limits: terminal_limits(),
    };
    let candidate = model
        .solve_terminals_with_motion_candidate(
            drive,
            ConductionMotion {
                velocity_m_s: &velocity,
                normal_flux_density_t: &[1.0; 6],
            },
            limits(),
            controls(),
            policy(),
            256,
            &mut || false,
        )
        .unwrap();
    // Three unit lengths: J=(1+0.2|E'|)E', plus original R''=0.5.
    let j = 0.2_f64;
    let ep = 2.0 * j / (1.0 + (1.0 + 0.8 * j).sqrt());
    let voltage = 2.0 - 3.0 * ep - 0.5 * j;
    let report = candidate.report();
    near(report.terminals.as_ref().unwrap().voltage_v[1], voltage);
    near(report.field.contact_response[0].power_w, 0.1 * j * j);
    near(report.field.power_w, 0.6 * ep * j + 0.1 * j * j);
    near(report.motion.mechanical_power_w, -0.4 * j);
    near(
        report.terminals.as_ref().unwrap().supplied_power_w,
        -0.2 * j * voltage,
    );
    let mut invalid_velocity = velocity;
    invalid_velocity[2] = [0.0, 1.0];
    assert!(matches!(
        model.solve_terminals_with_motion_candidate(
            drive,
            ConductionMotion {
                velocity_m_s: &invalid_velocity,
                normal_flux_density_t: &[0.0; 6]
            },
            limits(),
            controls(),
            policy(),
            256,
            &mut || false,
        ),
        Err(ConductionResponseError::Field(ConductionError::Invalid(
            "moving electrical contacts require an owned surface work and heat law"
        )))
    ));
}
