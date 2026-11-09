//! Nonuniform original contact laws on the same tetrahedral trace quadrature.
use super::*;
use fs_em::conduction::response::contact::{
    ContactQuadratureRefinementControls, ContactResponseQuadrature,
};

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

fn refinement_controls() -> ContactQuadratureRefinementControls {
    ContactQuadratureRefinementControls {
        maximum_refinements: 4,
        maximum_field_solves: 1024,
        maximum_source_evaluations: 1_000_000,
        absolute_voltage_tolerance_v: 1e-7,
        absolute_current_tolerance_a: 1e-7,
        absolute_power_tolerance_w: 1e-7,
        absolute_current_density_tolerance_a_m2: 1e-7,
        relative_tolerance: 1e-6,
    }
}

#[test]
fn automatic_triangle_refinement_closes_the_field_and_counts_every_discarded_solve() {
    let policy = refinement_controls();
    let quadrature = ContactResponseQuadrature {
        subdivisions: 1,
        maximum_points: 1634, // two original triangles through n=16
    };
    let model = specimen()
        .with_contact_quadrature_refinement(quadrature, policy)
        .unwrap();
    let boundary = varying_boundary(&model);
    let mut nonlinear = response_controls();
    nonlinear.maximum_field_solves = 1024;
    let mut polls = 0;
    let accepted = model
        .solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let refinement = accepted
        .report()
        .contact_quadrature_refinement
        .as_ref()
        .unwrap();
    assert!(refinement.differences.converged);
    assert!(refinement.refinements >= 1);
    assert_eq!(refinement.initial_subdivisions, 1);
    assert_eq!(
        refinement.final_subdivisions,
        1usize << refinement.refinements
    );
    assert_eq!(refinement.field_solves, accepted.report().field_solves);
    assert_eq!(
        refinement.source_evaluations,
        accepted.report().source_evaluations
    );
    assert!(refinement.field_solves > accepted.report().iterations);
    assert!(refinement.field_solves <= nonlinear.maximum_field_solves);

    // Independently replay the fixed-resolution public solves. Every rung
    // starts from the original source context, so these counts must add up to
    // the adaptive report, including the discarded coarse solutions.
    let mut total_fields = 0;
    let mut total_sources = 0;
    let mut previous = None;
    for level in 0..=refinement.refinements {
        let fixed = specimen()
            .with_contact_quadrature(ContactResponseQuadrature {
                subdivisions: 1usize << level,
                ..quadrature
            })
            .unwrap()
            .solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false)
            .unwrap();
        assert!(fixed.report().contact_quadrature_refinement.is_none());
        total_fields += fixed.report().field_solves;
        total_sources += fixed.report().source_evaluations;
        if level == refinement.refinements {
            assert_eq!(
                fixed.report().field.voltage_v,
                accepted.report().field.voltage_v
            );
            assert_eq!(
                fixed.report().field.contacts,
                accepted.report().field.contacts
            );
            let coarse: &fs_em::spatial::conduction::SpatialConductionField =
                previous.as_ref().unwrap();
            for (old, fine) in coarse.contacts.iter().zip(&fixed.report().field.contacts) {
                assert!(
                    (fine.current_a - old.current_a).abs()
                        <= policy.absolute_current_tolerance_a
                            + policy.relative_tolerance
                                * fine.current_a.abs().max(old.current_a.abs())
                );
                assert!(
                    (fine.power_w - old.power_w).abs()
                        <= policy.absolute_power_tolerance_w
                            + policy.relative_tolerance * fine.power_w.abs().max(old.power_w.abs())
                );
            }
        }
        previous = Some(fixed.report().field.clone());
    }
    assert_eq!(
        (total_fields, total_sources),
        (refinement.field_solves, refinement.source_evaluations)
    );
    near(
        accepted.report().field.power_w,
        accepted.report().field.supplied_power_w,
    );
    let expected_power = accepted.report().field.power_w;
    let mut calls = 0;
    assert!(matches!(
        model.solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || {
            calls += 1;
            calls == polls
        }),
        Err(SpatialConductionResponseError::Cancelled)
    ));

    let (owner, _) = accepted.into_parts();
    let (native, sources) = owner.into_frozen_parts();
    assert_eq!(sources.contacts().unwrap().quadrature(), Some(quadrature));
    assert_eq!(
        sources.contacts().unwrap().quadrature_refinement(),
        Some(policy)
    );
    let resumed =
        SpatialConductionResponse::from_frozen_parts(native, sources, &mut || false).unwrap();
    let next = resumed
        .solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false)
        .unwrap();
    assert!(next.report().contact_quadrature_refinement.is_some());
    near(next.report().field.power_w, expected_power);
    assert_eq!(
        model.contacts()[0]
            .source
            .resolved()
            .state_coordinate("normal_current"),
        Some(0.0)
    );
}

#[test]
fn unresolved_triangle_integration_and_global_work_limits_return_no_candidate() {
    let quadrature = ContactResponseQuadrature {
        subdivisions: 1,
        maximum_points: 1634,
    };
    let mut nonlinear = response_controls();
    nonlinear.maximum_field_solves = 1024;
    let mut exact = refinement_controls();
    exact.maximum_refinements = 1;
    exact.relative_tolerance = f64::EPSILON;
    exact.absolute_voltage_tolerance_v = 0.0;
    exact.absolute_current_tolerance_a = 0.0;
    exact.absolute_power_tolerance_w = 0.0;
    exact.absolute_current_density_tolerance_a_m2 = 0.0;
    let model = specimen()
        .with_contact_quadrature_refinement(quadrature, exact)
        .unwrap();
    let boundary = varying_boundary(&model);
    assert!(matches!(
        model.solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false),
        Err(SpatialConductionResponseError::Contact(
            ContactResponseError::QuadratureNotConverged { refinements: 1, .. }
        ))
    ));
    let mut limited = refinement_controls();
    limited.maximum_source_evaluations = 1;
    let model = specimen()
        .with_contact_quadrature_refinement(quadrature, limited)
        .unwrap();
    assert!(matches!(
        model.solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false),
        Err(SpatialConductionResponseError::Contact(
            ContactResponseError::QuadratureBudget {
                work: "source evaluations",
                completed: 0,
                ..
            }
        ))
    ));
    let model = specimen()
        .with_contact_quadrature_refinement(
            ContactResponseQuadrature {
                maximum_points: 14,
                ..quadrature
            },
            refinement_controls(),
        )
        .unwrap();
    assert!(
        model
            .solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false,)
            .is_err()
    );
    let model = specimen()
        .with_contact_quadrature_refinement(quadrature, refinement_controls())
        .unwrap();
    // The old native cap remains global even when the refinement policy has
    // a larger cap. One nonlinear trial cannot accept this source law.
    nonlinear.maximum_field_solves = 1;
    assert!(
        model
            .solve_candidate(&boundary, limits(), controls(), nonlinear, &mut || false,)
            .is_err()
    );
    assert_eq!(
        model.contact_response_sources().unwrap().quadrature(),
        Some(quadrature)
    );
}
