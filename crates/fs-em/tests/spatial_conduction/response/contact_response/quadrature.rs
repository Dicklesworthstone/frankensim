//! Nonuniform original contact laws on the same tetrahedral trace quadrature.
use super::*;
use fs_em::conduction::response::contact::ContactResponseQuadrature;

fn varying_boundary(model: &SpatialConductionResponse) -> SpatialConductionBoundary {
    let faces = face_sets(model, &[0.0, 2.0]);
    SpatialConductionBoundary {
        voltage: faces[0]
            .iter()
            .map(|&face| ConductionVoltageFace {
                face,
                voltage_v: [0.0; 3],
            })
            .chain(
                faces[1].iter().map(|&face| ConductionVoltageFace {
                    face,
                    voltage_v: model.mesh().faces[face]
                        .map(|v| 1.0 + 2.0 * model.positions_m()[v as usize][0]),
                }),
            )
            .collect(),
        current: Vec::new(),
    }
}

#[test]
fn sampled_triangle_contacts_close_every_original_local_law_and_native_work() {
    for subdivisions in [1usize, 2, 4] {
        let points_per_contact = 3 * subdivisions * subdivisions + 3 * subdivisions + 1;
        let model = specimen()
            .with_contact_quadrature(ContactResponseQuadrature {
                subdivisions,
                maximum_points: 2 * points_per_contact,
            })
            .unwrap();
        let boundary = varying_boundary(&model);
        for temperature in [300.0, 350.0] {
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
            let candidate = model
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
            let report = candidate.report();
            let mut total_contact_power = 0.0;
            for (index, (native, sampled)) in report
                .field
                .contacts
                .iter()
                .zip(&report.evaluated_contacts)
                .enumerate()
            {
                assert_eq!(native.quadrature.len(), points_per_contact);
                assert_eq!(sampled.evaluated.len(), points_per_contact);
                let mut area = 0.0;
                let mut current = 0.0;
                let mut heat = 0.0;
                let mut range = [f64::INFINITY, f64::NEG_INFINITY];
                for (point, source) in native.quadrature.iter().zip(&sampled.evaluated) {
                    near(point.basis.iter().sum(), 1.0);
                    assert!(point.area_weight_m2 > 0.0);
                    near(
                        point.voltage_jump_v,
                        point
                            .basis
                            .iter()
                            .zip(native.voltage_jump_v)
                            .map(|(basis, voltage)| basis * voltage)
                            .sum(),
                    );
                    // Independent local quadratic solution: v=(r0+beta*|j|)j.
                    // This is a pointwise law oracle, not an assumed face-average R.
                    let r0 = 0.5 + 0.002 * (temperature - 300.0);
                    let voltage = point.voltage_jump_v;
                    let density = voltage.signum() * 2.0 * voltage.abs()
                        / (r0 + (r0 * r0 + 3.0 * voltage.abs()).sqrt());
                    near(point.current_density_a_m2, density);
                    near(source.resistance_ohm_m2(), r0 + 0.75 * density.abs());
                    near(point.resistance_ohm_m2, source.resistance_ohm_m2());
                    near(
                        source.resolved().state_coordinate("T").unwrap(),
                        temperature,
                    );
                    near(
                        source
                            .resolved()
                            .state_coordinate("normal_current")
                            .unwrap(),
                        density.abs(),
                    );
                    assert_eq!(
                        source.resolved().properties()[0].answer().receipt.selected,
                        contacts[index].source.resolved().properties()[0]
                            .answer()
                            .receipt
                            .selected,
                    );
                    let power = point.area_weight_m2 * voltage * density;
                    near(point.power_w, power);
                    area += point.area_weight_m2;
                    current += point.area_weight_m2 * density;
                    heat += power;
                    range[0] = range[0].min(point.resistance_ohm_m2);
                    range[1] = range[1].max(point.resistance_ohm_m2);
                }
                near(area, 0.5);
                near(native.current_a, current);
                near(native.power_w, heat);
                assert!(range[1] - range[0] > 1e-3);
                assert!(
                    sampled.maximum_voltage_defect_v
                        <= interface_response::controls().voltage_tolerance_v
                );
                total_contact_power += heat;
            }
            near(
                report.field.power_w,
                report.field.cell_power_w.iter().sum::<f64>() + total_contact_power,
            );
            near(report.field.power_w, report.field.supplied_power_w);
            assert_eq!(report.field_solves, report.iterations);
            // Twelve fixed bulk sources are queried once per iteration;
            // every contact point is evaluated and, except on acceptance, relaxed.
            assert_eq!(
                report.source_evaluations,
                12 * report.iterations + 2 * points_per_contact * (2 * report.iterations - 1),
            );
        }
        assert_eq!(
            model.contacts()[0]
                .source
                .resolved()
                .state_coordinate("normal_current"),
            Some(0.0),
        );
    }
}

#[test]
fn sampled_contact_policy_and_cancellation_preserve_continuing_source_ownership() {
    for policy in [
        ContactResponseQuadrature {
            subdivisions: 0,
            maximum_points: 100,
        },
        ContactResponseQuadrature {
            subdivisions: 2,
            maximum_points: 37,
        },
        ContactResponseQuadrature {
            subdivisions: usize::MAX,
            maximum_points: usize::MAX,
        },
    ] {
        assert!(specimen().with_contact_quadrature(policy).is_err());
    }
    let model = specimen()
        .with_contact_quadrature(ContactResponseQuadrature {
            subdivisions: 2,
            maximum_points: 38,
        })
        .unwrap();
    let boundary = varying_boundary(&model);
    let mut polls = 0;
    let first = model
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            model.solve_candidate(
                &boundary,
                limits(),
                controls(),
                response_controls(),
                &mut || {
                    calls += 1;
                    calls == stop
                },
            ),
            Err(SpatialConductionResponseError::Cancelled)
        ));
    }
    let expected = first.report().field.clone();
    let (accepted, _) = first.into_parts();
    let (native, sources) = accepted.into_frozen_parts();
    let resumed =
        SpatialConductionResponse::from_frozen_parts(native, sources, &mut || false).unwrap();
    let retry = resumed
        .solve_candidate(
            &boundary,
            limits(),
            controls(),
            response_controls(),
            &mut || false,
        )
        .unwrap();
    for contact in &retry.report().field.contacts {
        assert_eq!(contact.quadrature.len(), 19);
    }
    near(retry.report().field.power_w, expected.power_w);
    for (actual, expected) in retry
        .report()
        .field
        .voltage_v
        .iter()
        .zip(&expected.voltage_v)
    {
        near(*actual, *expected);
    }
}
