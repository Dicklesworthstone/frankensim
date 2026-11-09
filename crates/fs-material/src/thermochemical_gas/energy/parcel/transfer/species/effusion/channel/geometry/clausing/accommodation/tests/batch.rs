use super::*;

#[test]
fn mixed_species_share_one_jointly_converged_geometry_in_authored_order() {
    let coefficients = [1.0, 0.25, 0.0, 0.8, 0.25];
    let (resolved, responses) = geometry(2.0)
        .resolve_thermal_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    let mass_report = resolved.resolution().unwrap();
    assert_eq!(responses.len(), coefficients.len());
    assert_eq!(responses[1], responses[4]);
    assert_eq!(
        responses[0].transmitted_heat_probability(),
        resolved.wall_transmitted_probability()
    );
    assert_eq!(
        responses[0].reflected_heat_probability(),
        resolved.reflected_probability()
    );
    assert_eq!(responses[2].transmitted_heat_probability(), 0.0);
    assert_eq!(responses[2].reflected_heat_probability(), 0.0);
    assert!(responses[0].resolution().is_none());
    assert!(responses[2].resolution().is_none());
    for (&coefficient, response) in coefficients.iter().zip(&responses) {
        assert_eq!(response.coefficient(), coefficient);
        if let Some(report) = response.resolution() {
            assert_eq!(report.rings, mass_report.rings);
            assert_eq!(report.refinements, mass_report.refinements);
            assert_eq!(report.batch_size, coefficients.len());
            assert_eq!(report.thermal_coefficients, 2);
            assert!(report.iterations > mass_report.iterations);
            assert!(
                report.transmitted_error_estimate
                    <= controls().absolute_tolerance
                        + controls().relative_tolerance * response.transmitted_heat_probability()
            );
            assert!(
                report.reflected_error_estimate
                    <= controls().absolute_tolerance
                        + controls().relative_tolerance * response.reflected_heat_probability()
            );
            let (individual, individual_response) = geometry(2.0)
                .resolve_thermal_accommodation(coefficient, controls(), || false)
                .unwrap();
            assert!(mass_report.rings >= individual.resolution().unwrap().rings);
            let individual_report = individual_response.resolution().unwrap();
            assert!(
                (response.transmitted_heat_probability()
                    - individual_response.transmitted_heat_probability())
                .abs()
                    <= report.transmitted_error_estimate
                        + individual_report.transmitted_error_estimate
            );
            assert!(
                (response.reflected_heat_probability()
                    - individual_response.reflected_heat_probability())
                .abs()
                    <= report.reflected_error_estimate + individual_report.reflected_error_estimate
            );
        }
    }
    assert!(
        responses[1].transmitted_heat_probability() < responses[3].transmitted_heat_probability()
    );
    assert!(
        responses[3].transmitted_heat_probability() < responses[0].transmitted_heat_probability()
    );
    assert!(responses[1].reflected_heat_probability() < responses[3].reflected_heat_probability());
    assert!(responses[3].reflected_heat_probability() < responses[0].reflected_heat_probability());
}

#[test]
fn duplicate_species_reuse_actual_solver_work_and_preserve_scalar_results() {
    let mut scalar_polls = 0;
    let (scalar_geometry, scalar) = geometry(1.0)
        .resolve_thermal_accommodation(0.3, controls(), || {
            scalar_polls += 1;
            false
        })
        .unwrap();
    let mut batch_polls = 0;
    let (batch_geometry, batch) = geometry(1.0)
        .resolve_thermal_accommodation_batch(&[0.3; MAX_EFFUSIVE_SPECIES], controls(), || {
            batch_polls += 1;
            false
        })
        .unwrap();
    assert_eq!(scalar_geometry, batch_geometry);
    assert_eq!(scalar_polls, batch_polls);
    let scalar_report = scalar.resolution().unwrap();
    for response in &batch {
        assert_eq!(response, &batch[0]);
        assert_eq!(
            response.transmitted_heat_probability(),
            scalar.transmitted_heat_probability()
        );
        assert_eq!(
            response.reflected_heat_probability(),
            scalar.reflected_heat_probability()
        );
        let report = response.resolution().unwrap();
        assert_eq!(report.batch_size, MAX_EFFUSIVE_SPECIES);
        assert_eq!(report.thermal_coefficients, 1);
        assert_eq!(report.iterations, scalar_report.iterations);
        assert_eq!(report.rings, scalar_report.rings);
        assert_eq!(
            report.transmitted_error_estimate,
            scalar_report.transmitted_error_estimate
        );
        assert_eq!(
            report.reflected_error_estimate,
            scalar_report.reflected_error_estimate
        );
    }
}

#[test]
fn species_permutation_preserves_physics_and_joint_work_deterministically() {
    let coefficients = [0.9, 1e-18, 0.3, 1.0, 0.0, 0.3];
    let original = geometry(1.0)
        .resolve_thermal_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    let repeated = geometry(1.0)
        .resolve_thermal_accommodation_batch(&coefficients, controls(), || false)
        .unwrap();
    assert_eq!(original, repeated);
    let mut reversed = coefficients;
    reversed.reverse();
    let (geometry_reversed, mut responses_reversed) = geometry(1.0)
        .resolve_thermal_accommodation_batch(&reversed, controls(), || false)
        .unwrap();
    responses_reversed.reverse();
    assert_eq!(original, (geometry_reversed, responses_reversed));
    assert!(original.1[1].transmitted_heat_probability() > 0.0);
    assert!(original.1[1].reflected_heat_probability() > 0.0);
}

#[test]
fn endpoint_only_and_zero_length_batches_preserve_supplied_geometry() {
    for original in [geometry(0.0), geometry(2.0)] {
        let (resolved, responses) = original
            .clone()
            .resolve_thermal_accommodation_batch(&[1.0, -0.0, 0.0, 1.0], controls(), || false)
            .unwrap();
        assert_eq!(original, resolved);
        for response in responses {
            assert!(response.resolution().is_none());
            assert_eq!(
                response.transmitted_heat_probability(),
                response.coefficient() * original.wall_transmitted_probability()
            );
            assert_eq!(
                response.reflected_heat_probability(),
                response.coefficient() * original.reflected_probability()
            );
        }
    }
    let original = geometry(0.0);
    let (resolved, responses) = original
        .clone()
        .resolve_thermal_accommodation_batch(&[0.3, 0.7, 1.0], controls(), || false)
        .unwrap();
    assert_eq!(resolved, original);
    assert!(
        responses
            .iter()
            .all(|response| response.transmitted_heat_probability() == 0.0
                && response.reflected_heat_probability() == 0.0
                && response.resolution().is_none())
    );
}

#[test]
fn batches_validate_all_inputs_before_work_and_refuse_exhausted_joint_budgets() {
    for coefficients in [Vec::new(), vec![0.3; MAX_EFFUSIVE_SPECIES + 1]] {
        assert_eq!(
            geometry(1.0).resolve_thermal_accommodation_batch(&coefficients, controls(), || false),
            Err(invalid(
                "accommodation batch exceeds effusive species bounds"
            ))
        );
    }
    for invalid_coefficient in [f64::NAN, f64::INFINITY, -0.1, 1.1] {
        let mut polls = 0;
        let result = geometry(1.0).resolve_thermal_accommodation_batch(
            &[0.3, invalid_coefficient],
            controls(),
            || {
                polls += 1;
                false
            },
        );
        assert_eq!(
            result,
            Err(invalid(
                "thermal accommodation coefficient must be finite and in [0,1]"
            ))
        );
        assert_eq!(polls, 1);
    }
    let exhausted = ClausingControls {
        initial_rings: 4,
        maximum_rings: 16,
        absolute_tolerance: 1e-12,
        relative_tolerance: 0.0,
        ..controls()
    };
    assert_eq!(
        geometry(2.0)
            .resolve_thermal_accommodation_batch(&[0.0, 0.1, 0.9, 1.0], exhausted, || false),
        Err(invalid("accommodation spatial refinement budget exhausted"))
    );
    let original = geometry(1.0);
    let mut polls = 0;
    assert_eq!(
        original
            .clone()
            .resolve_thermal_accommodation_batch(&[0.1, 0.3, 0.9], controls(), || {
                polls += 1;
                polls == 12
            }),
        Err(GasEnergyError::Cancelled)
    );
    assert_eq!(polls, 12);
    assert!(original.resolution().is_none());
}
