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
