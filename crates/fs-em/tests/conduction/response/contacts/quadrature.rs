//! G1/G3/G4: a nonuniform original contact law on one complete native field.
use super::*;
#[path = "../../../common/electrical_interface_response.rs"]
mod interface_response;
use fs_em::conduction::{
    ConductionCurrentTrace,
    response::contact::{ContactResponseError, ContactResponseQuadrature},
};

fn specimen() -> PlanarConductionResponse {
    let bulk = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.0,
        0.0,
        20.0,
    );
    let law = interface_response::voltage_jump_law(bulk.original(), bulk.original());
    PlanarConductionResponse::new_with_contact_response(
        mesh(0.4),
        vec![bulk; 4],
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source: law.original().clone(),
        }],
        vec![law],
        interface_response::controls(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap()
}

fn traces(sign: f64) -> (Vec<(usize, f64)>, Vec<ConductionCurrentTrace>) {
    (
        vec![(0, sign * 2.0), (3, sign * 3.0), (5, -sign), (6, -sign)],
        vec![
            ConductionCurrentTrace {
                vertices: [0, 1],
                inward_current_density_a_m2: -sign,
            },
            ConductionCurrentTrace {
                vertices: [2, 3],
                inward_current_density_a_m2: sign,
            },
        ],
    )
}

#[test]
fn sampled_original_contact_closes_nonuniform_fields_current_and_heat() {
    // Exact continuum and P1 fields:
    // A: V=-x+y+2, J=(1,-1); B: V=-x+1, J=(1,0).
    // The interface jump is 1+y, and original R''=|jump| gives J_n=1.
    // Bulk powers are depth*(2+1), contact power is depth*integral_0^1(1+y).
    for adaptive in [false, true] {
        for subdivisions in [1, 2, 4] {
            let quadrature = ContactResponseQuadrature {
                subdivisions,
                maximum_points: 2 * subdivisions + 1,
            };
            let mut model = specimen().with_contact_quadrature(quadrature).unwrap();
            if adaptive {
                model = model
                    .with_adaptive_relaxation(AdaptiveConductionRelaxation {
                        initial_weight: 0.5,
                        minimum_weight: 0.05,
                        maximum_weight: 1.0,
                    })
                    .unwrap();
            }
            for sign in [-1.0, 1.0] {
                let (voltage, current) = traces(sign);
                let solved = model
                    .solve_candidate(
                        &voltage,
                        &current,
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .unwrap();
                let report = solved.report();
                for (&actual, expected) in report
                    .voltage_v
                    .iter()
                    .zip([2.0, 1.0, 2.0, 3.0, 0.0, -1.0, -1.0, 0.0])
                {
                    near(actual, sign * expected);
                }
                for (cell, j) in report.current_density_a_m2.iter().enumerate() {
                    near(j[0], sign);
                    near(j[1], if cell < 2 { -sign } else { 0.0 });
                }
                let native = &report.contact_response[0];
                let closure = &report.evaluated_contacts[0];
                assert_eq!(native.quadrature.len(), quadrature.maximum_points);
                assert_eq!(closure.applied.len(), native.quadrature.len());
                assert_eq!(closure.evaluated.len(), native.quadrature.len());
                near(native.current_a, sign * 0.4);
                near(native.power_w, 0.6);
                near(report.cell_power_w[..2].iter().sum(), 0.8);
                near(report.cell_power_w[2..].iter().sum(), 0.4);
                near(report.power_w, 1.8);
                near(report.power_balance_residual_w, 0.0);
                for ((point, applied), evaluated) in native
                    .quadrature
                    .iter()
                    .zip(&closure.applied)
                    .zip(&closure.evaluated)
                {
                    let expected_jump = 1.0 + point.basis[1];
                    near(point.voltage_jump_v, sign * expected_jump);
                    near(point.current_density_a_m2, sign);
                    assert_eq!(point.resistance_ohm_m2, applied.resistance_ohm_m2());
                    near(evaluated.resistance_ohm_m2(), expected_jump);
                    near(
                        evaluated.resolved().state_coordinate("jump").unwrap(),
                        expected_jump,
                    );
                    near(point.power_w, point.area_weight_m2 * expected_jump);
                }
                near(
                    native.quadrature.iter().map(|p| p.power_w).sum(),
                    native.power_w,
                );
                assert!(closure.absolute_power_defect_w <= controls().power_tolerance_w);
                let minimum_queries = (4 + native.quadrature.len()) * report.field_solves;
                assert!(report.source_evaluations > minimum_queries);
                let original_power = report.power_w;
                let (native, sources) = solved.into_parts().0.into_frozen_parts();
                assert_eq!(sources.contacts().unwrap().quadrature(), Some(quadrature));
                let resumed =
                    PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false)
                        .unwrap();
                let again = resumed
                    .solve_candidate(
                        &voltage,
                        &current,
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .unwrap();
                near(again.report().power_w, original_power);
            }
        }
    }
}

#[test]
fn sampled_contact_bounds_and_cancellation_leave_original_sources_reusable() {
    let base = specimen();
    let (voltage, current) = traces(1.0);
    assert!(matches!(
        base.solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false
        ),
        Err(ConductionResponseError::Contact(
            ContactResponseError::NonuniformTrace { .. }
        ))
    ));
    for quadrature in [
        ContactResponseQuadrature {
            subdivisions: 0,
            maximum_points: 100,
        },
        ContactResponseQuadrature {
            subdivisions: 4,
            maximum_points: 8,
        },
        ContactResponseQuadrature {
            subdivisions: usize::MAX,
            maximum_points: usize::MAX,
        },
    ] {
        assert!(base.clone().with_contact_quadrature(quadrature).is_err());
    }
    let model = base
        .with_contact_quadrature(ContactResponseQuadrature {
            subdivisions: 4,
            maximum_points: 9,
        })
        .unwrap();
    let mut polls = 0;
    let accepted = model
        .solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            model.solve_candidate(
                &voltage,
                &current,
                limits(),
                controls(),
                policy(),
                &mut || {
                    seen += 1;
                    seen == stop
                },
            ),
            Err(ConductionResponseError::Cancelled)
        ));
    }
    let retry = model
        .solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(accepted.report().voltage_v, retry.report().voltage_v);
    assert_eq!(
        accepted.report().evaluated_contacts,
        retry.report().evaluated_contacts
    );
    assert_eq!(
        model.contacts()[0]
            .source
            .resolved()
            .state_coordinate("jump"),
        Some(1.5)
    );
    assert!(matches!(
        model.solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            ConductionResponseControls {
                maximum_iterations: 1,
                ..policy()
            },
            &mut || false,
        ),
        Err(ConductionResponseError::NotConverged { .. })
    ));
    let mut unsupported = voltage;
    unsupported[0].1 = 20.0;
    unsupported[1].1 = 20.0;
    assert!(matches!(
        model.solve_candidate(
            &unsupported,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        ),
        Err(ConductionResponseError::Contact(
            ContactResponseError::Material(_)
        ))
    ));
}

fn refinement_controls() -> fs_em::conduction::response::contact::ContactQuadratureRefinementControls
{
    fs_em::conduction::response::contact::ContactQuadratureRefinementControls {
        maximum_refinements: 5,
        maximum_field_solves: 128,
        maximum_source_evaluations: 100_000,
        absolute_voltage_tolerance_v: 1e-6,
        absolute_current_tolerance_a: 1e-7,
        absolute_power_tolerance_w: 1e-7,
        absolute_current_density_tolerance_a_m2: 1e-6,
        relative_tolerance: 1e-6,
    }
}

#[test]
fn refinement_compares_complete_fields_retains_policy_and_accounts_for_every_trial() {
    let base = specimen();
    let quadrature = ContactResponseQuadrature {
        subdivisions: 1,
        maximum_points: 65,
    };
    let (voltage, current) = traces(1.0);
    let mut expected_work = [0, 0];
    let mut first_fields = 0;
    for subdivisions in [1, 2] {
        let fixed = base
            .clone()
            .with_contact_quadrature(ContactResponseQuadrature {
                subdivisions,
                ..quadrature
            })
            .unwrap()
            .solve_candidate(
                &voltage,
                &current,
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap();
        expected_work[0] += fixed.report().field_solves;
        expected_work[1] += fixed.report().source_evaluations;
        if subdivisions == 1 {
            first_fields = fixed.report().field_solves;
        }
    }
    let model = base
        .clone()
        .with_contact_quadrature_refinement(quadrature, refinement_controls())
        .unwrap();
    let mut polls = 0;
    let accepted = model
        .solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let report = accepted.report();
    let refinement = report.contact_quadrature_refinement.as_ref().unwrap();
    assert_eq!(
        (
            refinement.initial_subdivisions,
            refinement.final_subdivisions
        ),
        (1, 2)
    );
    assert_eq!(refinement.refinements, 1);
    assert!(refinement.differences.converged);
    assert_eq!(
        [refinement.field_solves, refinement.source_evaluations],
        expected_work
    );
    assert_eq!(
        [report.field_solves, report.source_evaluations],
        expected_work
    );
    assert!(report.field_solves <= policy().maximum_iterations);
    assert_eq!(report.contact_response[0].quadrature.len(), 5);
    near(report.power_w, 1.8);
    near(report.contact_response[0].current_a, 0.4);
    near(report.contact_response[0].power_w, 0.6);
    let (native, sources) = accepted.into_parts().0.into_frozen_parts();
    assert_eq!(sources.contacts().unwrap().quadrature(), Some(quadrature));
    assert_eq!(
        sources.contacts().unwrap().quadrature_refinement(),
        Some(refinement_controls())
    );
    let restored =
        PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false).unwrap();
    let again = restored
        .solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    near(again.report().power_w, 1.8);
    assert_eq!(
        again
            .report()
            .contact_quadrature_refinement
            .as_ref()
            .unwrap()
            .initial_subdivisions,
        1
    );
    let fixed = restored.with_contact_quadrature(quadrature).unwrap();
    assert_eq!(
        fixed
            .contact_response_sources()
            .unwrap()
            .quadrature_refinement(),
        None
    );

    // The first rung can converge, but there is no authority to return it when
    // the work cap cannot admit the required independent second resolution.
    let short = base
        .clone()
        .with_contact_quadrature_refinement(
            quadrature,
            fs_em::conduction::response::contact::ContactQuadratureRefinementControls {
                maximum_field_solves: first_fields,
                ..refinement_controls()
            },
        )
        .unwrap();
    assert!(matches!(short.solve_candidate(
        &voltage, &current, limits(), controls(), policy(), &mut || false,
    ), Err(ConductionResponseError::Contact(ContactResponseError::QuadratureBudget {
        work: "field solves", completed, ..
    })) if completed == first_fields));
    let sources = base
        .with_contact_quadrature_refinement(
            quadrature,
            fs_em::conduction::response::contact::ContactQuadratureRefinementControls {
                maximum_source_evaluations: 13, // one complete iteration requires 2*(4+3)
                ..refinement_controls()
            },
        )
        .unwrap();
    assert!(matches!(
        sources.solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        ),
        Err(ConductionResponseError::Contact(
            ContactResponseError::QuadratureBudget {
                work: "source evaluations",
                completed: 0,
                ..
            }
        ))
    ));
    let mut seen = 0;
    assert!(matches!(
        model.solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || {
                seen += 1;
                seen == polls
            },
        ),
        Err(ConductionResponseError::Cancelled)
    ));
    near(
        model
            .solve_candidate(
                &voltage,
                &current,
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap()
            .report()
            .power_w,
        1.8,
    );
}

#[test]
fn unresolved_nonlinear_contact_refines_and_matches_independent_surface_integrals() {
    let bulk = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.0,
        0.0,
        20.0,
    );
    let beta = 0.75;
    let law = interface_response::law(bulk.original(), bulk.original(), beta);
    let base = PlanarConductionResponse::new_with_contact_response(
        mesh(0.4),
        vec![bulk; 4],
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source: law.original().clone(),
        }],
        vec![law],
        interface_response::controls(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap();
    let quadrature = ContactResponseQuadrature {
        subdivisions: 1,
        maximum_points: 65,
    };
    let (voltage, current) = traces(1.0);
    let model = base
        .clone()
        .with_contact_quadrature_refinement(quadrature, refinement_controls())
        .unwrap();
    let accepted = model
        .solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .unwrap();
    let report = accepted.report();
    let refinement = report.contact_quadrature_refinement.as_ref().unwrap();
    assert!(refinement.refinements > 1);
    assert!(refinement.final_subdivisions >= 4);
    assert!(refinement.differences.converged);
    assert!(report.field_solves <= policy().maximum_iterations);
    let contact = &report.contact_response[0];
    assert_eq!(
        contact.quadrature.len(),
        2 * refinement.final_subdivisions + 1
    );

    // Independently integrate J(v) = (sqrt(a²+4*beta*v)-a)/(2*beta)
    // and v*J(v) over the actual affine positive voltage jump. These analytic
    // primitives do not reuse the production quadrature or sampled resistances.
    let [v0, v1] = contact.voltage_jump_v;
    assert!(v0 > 0.0 && v1 - v0 > 0.1);
    let a = 0.5_f64;
    let q0 = (a * a + 4.0 * beta * v0).sqrt();
    let q1 = (a * a + 4.0 * beta * v1).sqrt();
    let delta = v1 - v0;
    let mean_j = ((q1.powi(3) - q0.powi(3)) / (6.0 * beta * delta) - a) / (2.0 * beta);
    let primitive = |q: f64| {
        q.powi(5) / 5.0 - a * q.powi(4) / 4.0 - a * a * q.powi(3) / 3.0 + a * a * a * q * q / 2.0
    };
    let mean_power = (primitive(q1) - primitive(q0)) / (16.0 * beta.powi(3) * delta);
    assert!((contact.current_a - 0.4 * mean_j).abs() < 2e-7);
    assert!((contact.power_w - 0.4 * mean_power).abs() < 2e-7);
    near(report.power_balance_residual_w, 0.0);

    let unresolved = base
        .with_contact_quadrature_refinement(
            quadrature,
            fs_em::conduction::response::contact::ContactQuadratureRefinementControls {
                maximum_refinements: 1,
                ..refinement_controls()
            },
        )
        .unwrap();
    assert!(matches!(
        unresolved.solve_candidate(
            &voltage,
            &current,
            limits(),
            controls(),
            policy(),
            &mut || false,
        ),
        Err(ConductionResponseError::Contact(
            ContactResponseError::QuadratureNotConverged { refinements: 1, .. }
        ))
    ));
}
