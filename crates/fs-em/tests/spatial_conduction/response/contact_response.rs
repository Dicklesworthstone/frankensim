//! Actual local interface response on original tetrahedral contact traces.
#[path = "contact_response/quadrature.rs"]
mod quadrature;
use super::*;
#[path = "../../common/electrical_interface_response.rs"]
mod interface_response;
use fs_em::conduction::response::contact::ContactResponseError;

fn specimen() -> SpatialConductionResponse {
    let bulk = law(None, 0., 20.);
    let interface = interface_response::law(bulk.original(), bulk.original(), 0.75);
    let (mut positions, mut tets) = slab(&[0., 1.]);
    let (other, extra) = slab(&[1., 2.]);
    positions.extend(other);
    tets.extend(extra.into_iter().map(|t| t.map(|v| v + 8)));
    rotate(&mut positions);
    let contacts = contact_faces(interface.original().clone());
    assert!(
        SpatialConductionResponse::new_with_contacts(
            positions.clone(),
            tets.clone(),
            vec![bulk.clone(); 12],
            contacts.clone(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    SpatialConductionResponse::new_with_contact_response(
        positions,
        tets,
        vec![bulk; 12],
        contacts,
        vec![interface; 2],
        interface_response::controls(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap()
}

#[test]
fn tetrahedral_contact_current_voltage_and_heat_follow_actual_original_source_states() {
    let adaptation = AdaptiveConductionRelaxation {
        initial_weight: 0.5,
        minimum_weight: 0.02,
        maximum_weight: 1.0,
    };
    let model = specimen().with_adaptive_relaxation(adaptation).unwrap();
    let faces = face_sets(&model, &[0., 2.]);
    for temperature in [300., 350.] {
        let mut contacts = model.contacts().to_vec();
        for (contact, law) in contacts
            .iter_mut()
            .zip(model.contact_response_sources().unwrap().laws())
        {
            contact.source = law
                .at_query_point(
                    &contact
                        .source
                        .resolved()
                        .query()
                        .with("T", temperature)
                        .unwrap(),
                )
                .unwrap();
        }
        for sign in [-1., 1.] {
            let j = interface_response::density(2., temperature, 0.75);
            let rc = 0.5 + 0.002 * (temperature - 300.) + 0.75 * j;
            for current_control in [false, true] {
                let mut boundary = SpatialConductionBoundary {
                    voltage: faces[0]
                        .iter()
                        .map(|&face| ConductionVoltageFace {
                            face,
                            voltage_v: [0.; 3],
                        })
                        .collect(),
                    current: Vec::new(),
                };
                if current_control {
                    boundary.current = faces[1]
                        .iter()
                        .map(|&face| ConductionCurrentFace {
                            face,
                            inward_current_density_a_m2: sign * j,
                        })
                        .collect();
                } else {
                    boundary
                        .voltage
                        .extend(faces[1].iter().map(|&face| ConductionVoltageFace {
                            face,
                            voltage_v: [sign * 2.; 3],
                        }));
                }
                let solved = model
                    .solve_with_materials_and_contacts_candidate(
                        model.conductors().to_vec(),
                        contacts.clone(),
                        &boundary,
                        limits(),
                        controls(),
                        response_controls(),
                        &mut || false,
                    )
                    .unwrap();
                let r = solved.report();
                near(r.field.power_w, 2. * j);
                near(r.field.supplied_power_w, 2. * j);
                near(r.field.cell_power_w.iter().sum::<f64>(), 2. * j * j);
                near(
                    r.field.contacts.iter().map(|c| c.power_w).sum::<f64>(),
                    rc * j * j,
                );
                assert_eq!(r.evaluated_contacts.len(), 2);
                for (heat, sample) in r.field.contacts.iter().zip(&r.evaluated_contacts) {
                    near(heat.current_a, -sign * 0.5 * j);
                    for (&v, &density) in heat.voltage_jump_v.iter().zip(&heat.current_density_a_m2)
                    {
                        near(v, -sign * rc * j);
                        near(density, -sign * j);
                    }
                    for source in &sample.evaluated {
                        near(
                            source
                                .resolved()
                                .state_coordinate("normal_current")
                                .unwrap(),
                            j,
                        );
                        assert_eq!(source.resolved().state_coordinate("T"), Some(temperature));
                    }
                    assert!(
                        sample.maximum_voltage_defect_v
                            <= interface_response::controls().voltage_tolerance_v
                    );
                }
                assert_eq!(solved.problem().relaxation(), 0.5);
                assert_eq!(solved.problem().adaptive_relaxation(), Some(adaptation));
                assert_eq!(
                    solved
                        .problem()
                        .contact_response_sources()
                        .unwrap()
                        .laws()
                        .len(),
                    2
                );
            }
        }
    }
    assert_eq!(
        model.contacts()[0]
            .source
            .resolved()
            .state_coordinate("normal_current"),
        Some(0.)
    );
}

#[test]
fn nonuniform_tetrahedral_contact_and_interrupted_source_closure_publish_no_candidate() {
    let model = specimen();
    let faces = face_sets(&model, &[0., 2.]);
    let mut boundary = SpatialConductionBoundary {
        voltage: faces[0]
            .iter()
            .map(|&face| ConductionVoltageFace {
                face,
                voltage_v: [0.; 3],
            })
            .collect(),
        current: Vec::new(),
    };
    boundary
        .voltage
        .extend(faces[1].iter().map(|&face| ConductionVoltageFace {
            face,
            voltage_v:
                model.mesh().faces[face].map(|v| 1. + 2. * model.positions_m()[v as usize][0]),
        }));
    assert!(matches!(
        model.solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || false
        ),
        Err(SpatialConductionResponseError::Contact(
            ContactResponseError::NonuniformTrace { .. }
        ))
    ));
    for voltage in &mut boundary.voltage {
        if faces[1].contains(&voltage.face) {
            voltage.voltage_v = [2.; 3];
        }
    }
    let mut calls = 0;
    let accepted = model
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut polled = 0;
    assert!(matches!(
        model.solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || {
                polled += 1;
                polled == calls / 2
            }
        ),
        Err(SpatialConductionResponseError::Cancelled)
    ));
    let retry = model
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        accepted.report().field.voltage_v,
        retry.report().field.voltage_v
    );
    assert_eq!(
        accepted.report().field.contacts,
        retry.report().field.contacts
    );
}
