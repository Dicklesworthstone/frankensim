//! A stationary resistive joint remains physical beside moving bulk material.
use super::*;

fn moving_bulk_specimen(depth: f64) -> PlanarConduction {
    let left = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let right = conductor(4.0, ElectricalPropertyBasis::Conductivity);
    let source = interface(&left, &right, 1e5, 1.0);
    let mesh = TriComplex2::from_indexed_triangles(
        "moving bulk with stationary electrical joint",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [2.0, 1.0],
        ],
        vec![
            [0, 1, 2],
            [0, 2, 3],
            [1, 4, 5],
            [1, 5, 2],
            [6, 7, 8],
            [6, 8, 9],
        ],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap();
    PlanarConduction::new_with_contacts(
        mesh,
        vec![
            left.clone(),
            left.clone(),
            left.clone(),
            left,
            right.clone(),
            right,
        ],
        vec![ConductionContact {
            first_vertices: [4, 5],
            second_vertices: [6, 9],
            source,
        }],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn moving_bulk_and_stationary_joint_close_short_generator_and_motor_work() {
    let area = 0.2;
    let problem = moving_bulk_specimen(area);
    // Series area-resistance is 2/2 + 1/4 + R'' = 1.75 ohm m².
    // Only x in [0,1] moves; both cells incident on the x=2 joint are fixed.
    for speed in [0.0, 2.0] {
        let velocity = [
            [0.0, speed],
            [0.0, speed],
            [0.0; 2],
            [0.0; 2],
            [0.0; 2],
            [0.0; 2],
        ];
        for voltage in [0.0, 1.0, 3.0] {
            let boundary = [(0, 0.0), (3, 0.0), (7, voltage), (8, voltage)];
            let field = problem
                .solve_with_motion(
                    &boundary,
                    &[],
                    ConductionMotion {
                        velocity_m_s: &velocity,
                        normal_flux_density_t: &[1.0; 6],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let j = (speed - voltage) / 1.75;
            let bulk_heat = area * j * j * 1.25;
            let contact_heat = area * j * j * 0.5;
            let mechanical = -area * j * speed;
            for cell in 0..6 {
                let resistivity = if cell < 4 { 0.5 } else { 0.25 };
                let motional = if cell < 2 { speed } else { 0.0 };
                near(field.current_density_a_m2()[cell][0], j);
                near(field.current_density_a_m2()[cell][1], 0.0);
                near(
                    field.electric_field_v_m()[cell][0],
                    resistivity * j - motional,
                );
                near(field.conductor_field_v_m()[cell][0], resistivity * j);
                near(
                    field.cell_joule_power_w()[cell],
                    0.5 * area * resistivity * j * j,
                );
            }
            let contact = &field.contact_response()[0];
            for jump in contact.voltage_jump_v {
                near(jump, 0.5 * j);
            }
            for current in contact.current_density_a_m2 {
                near(current, j);
            }
            near(contact.current_a, area * j);
            near(contact.power_w, contact_heat);
            near(field.contact_power_w(), contact_heat);
            near(field.joule_power_w(), bulk_heat + contact_heat);
            near(field.mechanical_power_w(), mechanical);
            near(field.electrical_power_w(), bulk_heat + mechanical);
            near(field.terminal_power_w(), -voltage * area * j);
            near(
                field.electrical_power_w() + field.contact_power_w(),
                field.terminal_power_w(),
            );
            near(field.power_balance_residual_w(), 0.0);
            near(field.work_pairing_residual_w(), 0.0);
            near(field.boundary_current_a(), 0.0);
            if speed == 0.0 {
                let stationary = problem.solve(&boundary, controls(), &mut || false).unwrap();
                near(field.joule_power_w(), stationary.power_w());
                assert_eq!(field.contact_response(), stationary.contact_response());
            } else if voltage == 0.0 {
                near(field.terminal_power_w(), 0.0);
                assert!(field.mechanical_power_w() < 0.0);
            } else if voltage < speed {
                assert!(field.terminal_power_w() < 0.0);
            } else {
                assert!(field.mechanical_power_w() > 0.0);
            }
        }
    }
}

#[test]
fn one_reference_drives_current_across_a_stationary_joint_below_voltage_spacing() {
    let problem = moving_bulk_specimen(0.2);
    for j in [-0.8, 0.8] {
        for bias in [0.0, 2.0_f64.powi(54)] {
            let field = problem
                .solve_with_motion(
                    &[(0, bias), (3, bias)],
                    &[ConductionCurrentTrace {
                        vertices: [7, 8],
                        inward_current_density_a_m2: -j,
                    }],
                    ConductionMotion {
                        velocity_m_s: &[
                            [0.0, 2.0],
                            [0.0, 2.0],
                            [0.0; 2],
                            [0.0; 2],
                            [0.0; 2],
                            [0.0; 2],
                        ],
                        normal_flux_density_t: &[1.0; 6],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            near(field.relative_voltage_v()[4], 2.0 - j);
            near(field.relative_voltage_v()[6], 2.0 - 1.5 * j);
            near(field.relative_voltage_v()[7], 2.0 - 1.75 * j);
            near(field.contact_power_w(), 0.1 * j * j);
            near(field.joule_power_w(), 0.35 * j * j);
            near(field.mechanical_power_w(), -0.4 * j);
            near(field.terminal_power_w(), -0.2 * j * (2.0 - 1.75 * j));
            assert!(field.reference_voltage_v().iter().all(|&v| v == bias));
            near(field.power_balance_residual_w(), 0.0);
        }
    }
}

#[test]
fn either_moving_contact_side_refuses_even_without_a_magnetic_field() {
    let problem = moving_bulk_specimen(0.2);
    for cell in [2, 5] {
        let mut velocity = [[0.0; 2]; 6];
        velocity[cell] = [0.0, 1.0];
        assert!(matches!(
            problem.solve_with_motion(
                &[(0, 0.0), (3, 0.0), (7, 1.0), (8, 1.0)],
                &[],
                ConductionMotion {
                    velocity_m_s: &velocity,
                    normal_flux_density_t: &[0.0; 6],
                },
                controls(),
                &mut || false,
            ),
            Err(ConductionError::Invalid(
                "moving electrical contacts require an owned surface work and heat law"
            ))
        ));
    }
}
