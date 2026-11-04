//! The fixed electrical joint retains its own current and heat under bulk motion.
use super::*;

fn specimen() -> SpatialConduction {
    let sources = [0.5, 0.25].map(contacts::conductor);
    let source = contacts::interface(&sources[0], &sources[1], 1e5, 1.0 / 6.0);
    let (mut positions, mut tets) = slab(&[0.0, 1.0, 2.0]);
    let (other, other_tets) = slab(&[2.0, 3.0]);
    positions.extend(other);
    tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + 12)));
    SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..18)
            .map(|cell| sources[usize::from(cell >= 12)].clone())
            .collect(),
        vec![
            SpatialConductionContact {
                first_vertices: [8, 9, 11],
                second_vertices: [12, 13, 15],
                source: source.clone(),
            },
            SpatialConductionContact {
                first_vertices: [8, 10, 11],
                second_vertices: [12, 14, 15],
                source,
            },
        ],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn stationary_contact_closes_native_motional_surface_and_volume_work() {
    let problem = specimen();
    for speed in [0.0, 2.0] {
        let mut velocity = [[0.0; 3]; 18];
        velocity[..6].fill([0.0, speed, 0.0]);
        for voltage in [0.0, 1.0, 3.0] {
            let boundary = ends(&problem, 3.0, Some(voltage), 0.0);
            let field = problem
                .solve_with_motion(
                    &boundary,
                    SpatialConductionMotion {
                        velocity_m_s: &velocity,
                        flux_density_t: &[[0.0, 0.0, 1.0]; 18],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Independent one-dimensional series law, unit cross section.
            // Bulk rho*L=1.25, stationary R''=0.5, motional emf=speed.
            let j = (speed - voltage) / 1.75;
            let bulk_heat = 1.25 * j * j;
            let contact_heat = 0.5 * j * j;
            let mechanical = -j * speed;
            for cell in 0..18 {
                let resistivity = if cell < 12 { 0.5 } else { 0.25 };
                let motional = if cell < 6 { speed } else { 0.0 };
                near(field.current_density_a_m2()[cell][0], j);
                near(field.current_density_a_m2()[cell][1], 0.0);
                near(field.current_density_a_m2()[cell][2], 0.0);
                near(
                    field.electric_field_v_m()[cell][0],
                    resistivity * j - motional,
                );
                near(field.conductor_field_v_m()[cell][0], resistivity * j);
                near(
                    field.cell_joule_power_w()[cell],
                    problem.cell_volume_m3(cell).unwrap() * resistivity * j * j,
                );
            }
            assert_eq!(field.contact_response().len(), 2);
            for contact in field.contact_response() {
                for jump in contact.voltage_jump_v {
                    near(jump, 0.5 * j);
                }
                for current in contact.current_density_a_m2 {
                    near(current, j);
                }
                near(contact.current_a, 0.5 * j);
                near(contact.power_w, 0.5 * contact_heat);
            }
            near(field.contact_power_w(), contact_heat);
            near(field.joule_power_w(), bulk_heat + contact_heat);
            near(field.mechanical_power_w(), mechanical);
            near(field.electrical_power_w(), bulk_heat + mechanical);
            near(field.terminal_power_w(), -voltage * j);
            near(
                field.electrical_power_w() + field.contact_power_w(),
                field.terminal_power_w(),
            );
            near(field.power_balance_residual_w(), 0.0);
            near(field.hodge_power_residual_w(), 0.0);
            near(field.work_pairing_residual_w(), 0.0);
            near(field.boundary_current_a(), 0.0);
            if speed == 0.0 {
                let stationary = problem.solve(&boundary, controls(), &mut || false).unwrap();
                near(field.joule_power_w(), stationary.power_w);
                assert_eq!(field.contact_response(), stationary.contacts.as_slice());
            } else if voltage == 0.0 {
                near(field.terminal_power_w(), 0.0);
                assert!(field.mechanical_power_w() < 0.0);
            } else if voltage < speed {
                assert!(field.terminal_power_w() < 0.0);
            } else {
                assert!(field.mechanical_power_w() > 0.0);
            }
            // The coupled thermal owner consumes this same accepted surface data.
            let contacts_ptr = field.contact_response().as_ptr();
            let snapshot = field.into_snapshot();
            assert_eq!(contacts_ptr, snapshot.contact_response().as_ptr());
            near(snapshot.contact_power_w(), contact_heat);
            near(snapshot.joule_power_w(), bulk_heat + contact_heat);
            near(
                snapshot.electrical_power_w() + snapshot.contact_power_w(),
                snapshot.terminal_power_w(),
            );
        }
    }
}

#[test]
fn moving_either_original_contact_cell_still_requires_a_surface_law() {
    let problem = specimen();
    for &cell in problem.contact_cells().iter().flatten() {
        for field in [[0.0; 3], [0.0, 0.0, 1.0]] {
            let mut velocity = [[0.0; 3]; 18];
            velocity[cell] = [0.0, 1.0, 0.0];
            assert!(matches!(
                problem.solve_with_motion(
                    &ends(&problem, 3.0, Some(1.0), 0.0),
                    SpatialConductionMotion {
                        velocity_m_s: &velocity,
                        flux_density_t: &[field; 18],
                    },
                    controls(),
                    &mut || false,
                ),
                Err(SpatialConductionError::Invalid(
                    "moving electrical contacts require an owned surface work and heat law"
                ))
            ));
        }
    }
}
