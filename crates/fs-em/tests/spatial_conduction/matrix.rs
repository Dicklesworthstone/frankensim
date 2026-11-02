//! Three-terminal reactions and local bulk/contact heat against series algebra.
use super::contacts::{conductor, contacts as material_contacts, interface};
use super::*;
use std::collections::BTreeSet;

fn response_controls() -> SpatialConductanceControls {
    SpatialConductanceControls {
        maximum_terminals: 3,
        maximum_field_solves: 2,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 1.0,
        coefficient_tolerance_s: 1e-10,
    }
}

fn reaction(problem: &SpatialConduction, field: &SpatialConductionField, faces: &[usize]) -> f64 {
    faces
        .iter()
        .flat_map(|&f| problem.mesh().faces[f])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|v| field.nodal_boundary_current_a[v as usize].unwrap())
        .sum()
}

#[test]
fn native_multiterminal_response_retains_cross_terms_contact_work_and_floating_electrode() {
    // Bulk A has rho=1, length=1, area=1 and one R''=1 contact: R_A=2.
    // Bulk B has rho=3, length=1, area=1: R_B=3.
    // The middle ideal electrode joins their otherwise disconnected volumes.
    let materials = [conductor(1.0), conductor(3.0)];
    let source = interface(&materials[0], &materials[0], 1e5, 1.0 / 3.0);
    let (mut positions, mut tets) = slab(&[0.0, 0.5]);
    for xs in [&[0.5, 1.0][..], &[1.0, 1.5, 2.0][..]] {
        let (part, cells) = slab(xs);
        let offset = positions.len() as u32;
        positions.extend(part);
        tets.extend(cells.into_iter().map(|t| t.map(|v| v + offset)));
    }
    let problem = SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..24)
            .map(|i| materials[usize::from(i >= 12)].clone())
            .collect(),
        material_contacts(source),
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = [0.0, 1.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&face| {
                problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect::<Vec<_>>()
    });
    let terminals = faces.each_ref().map(Vec::as_slice);
    let expected_matrix = [
        0.5,
        -0.5,
        0.0,
        -0.5,
        5.0 / 6.0,
        -1.0 / 3.0,
        0.0,
        -1.0 / 3.0,
        1.0 / 3.0,
    ];
    for reference_voltage_v in [0.3, 2.0] {
        let matrix = problem
            .conductance_matrix(
                &terminals,
                SpatialConductanceControls {
                    reference_voltage_v,
                    ..response_controls()
                },
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(matrix.field_solves(), 2);
        assert_eq!(matrix.terminal_faces(), &faces);
        for (&actual, expected) in matrix.conductance_s().iter().zip(expected_matrix) {
            near(actual, expected);
        }
        for sign in [-1.0, 1.0] {
            // Both nonreference terminals are driven simultaneously. B's
            // actual 2 V drop includes the negative cross term between bases.
            let field = matrix
                .field_at_voltages(
                    &[7.0, 7.0 + sign * 2.0, 7.0 + sign * 4.0],
                    1e-10,
                    1e-9,
                    &mut || false,
                )
                .unwrap();
            let current = matrix
                .field_at_currents(
                    &[-sign, sign / 3.0, 2.0 * sign / 3.0],
                    1e-10,
                    1e-9,
                    &mut || false,
                )
                .unwrap();
            for (i, (a, b)) in field.voltage_v.iter().zip(&current.voltage_v).enumerate() {
                near(*a, *b);
                let x = problem.positions_m()[i][0];
                let exact = if i < 8 {
                    x
                } else if i < 16 {
                    x + 1.0
                } else {
                    2.0 * x
                };
                near(*a, sign * exact);
            }
            for (i, density) in field.current_density_a_m2.iter().enumerate() {
                near(density[0], -sign * if i < 12 { 1.0 } else { 2.0 / 3.0 });
                near(density[1], 0.0);
                near(density[2], 0.0);
                near(
                    field.cell_power_w[i],
                    if i < 12 { 1.0 / 12.0 } else { 1.0 / 9.0 },
                );
            }
            for contact in &field.contacts {
                near(contact.current_a, -sign * 0.5);
                near(contact.power_w, 0.5);
            }
            for (terminal, expected) in terminals.iter().zip([-sign, sign / 3.0, 2.0 * sign / 3.0])
            {
                near(reaction(&problem, &field, terminal), expected);
            }
            near(field.power_w, 10.0 / 3.0);
            near(field.supplied_power_w, field.power_w);
            near(current.power_w, field.power_w);

            // I_middle=0 means a solved floating conductor at 2 V, with +1 A
            // and -1 A on its two faces. Setting its local current to zero or
            // grounding it would give a different field and total heat.
            let floating = matrix
                .field_at_currents(&[-sign, 0.0, sign], 1e-10, 1e-9, &mut || false)
                .unwrap();
            near(floating.voltage_v[12], sign * 2.0);
            near(floating.voltage_v[24], sign * 5.0);
            near(reaction(&problem, &floating, &faces[1]), 0.0);
            near(floating.current_density_a_m2[0][0], -sign);
            near(floating.current_density_a_m2[12][0], -sign);
            near(floating.power_w, 5.0);
        }
        let cancellation = matrix
            .field_at_voltages(&[0.0, 3.0, 3.0], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(cancellation.cell_power_w[12..].iter().sum(), 0.0);
        near(cancellation.power_w, 4.5);
        let zero = matrix
            .field_at_currents(&[0.0; 3], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(zero.power_w, 0.0);
        let biased = matrix
            .field_at_voltages(&[1e10, 1e10 + 2.0, 1e10 + 4.0], 1e-10, 1e-9, &mut || false)
            .unwrap();
        near(biased.power_w, 10.0 / 3.0);
        assert!(
            matrix
                .field_at_currents(&[0.0, 0.0, 1.0], 1e-10, 1e-9, &mut || false)
                .is_err()
        );
        assert!(
            matrix
                .field_at_currents(&[0.0, f64::NAN, 0.0], 1e-10, 1e-9, &mut || false)
                .is_err()
        );

        let mut polls = 0;
        matrix
            .field_at_currents(&[-1.0, 0.0, 1.0], 1e-10, 1e-9, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        let mut retry_polls = 0;
        assert!(matches!(
            matrix.field_at_currents(&[-1.0, 0.0, 1.0], 1e-10, 1e-9, &mut || {
                retry_polls += 1;
                retry_polls == polls
            }),
            Err(SpatialConductionError::Cancelled)
        ));
        near(
            matrix
                .field_at_currents(&[-1.0, 0.0, 1.0], 1e-10, 1e-9, &mut || false)
                .unwrap()
                .power_w,
            5.0,
        );
    }
    // Without the shared middle electrode these two terminal components have
    // no electrical connection. A numerical factor must not invent one.
    assert!(
        problem
            .conductance_matrix(
                &[&faces[0], &faces[2]],
                response_controls(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert!(
        problem
            .conductance_matrix(
                &[&faces[0], &faces[0], &faces[2]],
                response_controls(),
                controls(),
                &mut || false,
            )
            .is_err()
    );
    for bounded in [
        SpatialConductanceControls {
            maximum_field_solves: 1,
            ..response_controls()
        },
        SpatialConductanceControls {
            maximum_response_values: 1,
            ..response_controls()
        },
        SpatialConductanceControls {
            maximum_dense_work: 1,
            ..response_controls()
        },
    ] {
        assert!(
            problem
                .conductance_matrix(&terminals, bounded, controls(), &mut || false)
                .is_err()
        );
    }
}
