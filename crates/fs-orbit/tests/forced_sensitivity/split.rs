use super::*;

#[test]
fn split_physical_columns_share_cached_projection_and_refuse_incomplete_or_aliased_ports() {
    let problem = Coupled::new(1.1, 1.0);
    let orbit = solve(&problem, 7);
    let mut prepared =
        PreparedForcedSensitivity::new(&problem, &orbit, controls(128), &scales(2), || false)
            .unwrap();
    let source_calls = problem.source_calls.get();
    let tangent_calls = problem.tangent_calls.get();
    let source = |t: f64, _: &[f64], out: &mut [f64]| {
        out.copy_from_slice(&[det::cos(0.9 * t), 0.0]);
        true
    };
    let zero_port = |_: fs_orbit::ForcedParameterPoint<'_>, out: &mut [f64]| {
        out.fill(0.0);
        true
    };
    let expected = prepared
        .solve_island_direction(direction(1.0), source, || false)
        .unwrap();
    let actual = prepared
        .solve_split_direction(direction(1.0), source, zero_port, || false)
        .unwrap();
    assert_eq!(actual.coeffs, expected.coeffs);
    assert!(matches!(
        prepared.solve_split_direction(
            direction(1.0),
            source,
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        prepared.solve_split_direction(
            direction(1.0),
            source,
            |point, out| {
                out.fill(0.0);
                out[0] = point.aft_samples as f64;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    let gradient = observation_gradient(&orbit);
    let reference = objective_direction(&gradient, &expected.coeffs);
    let mut adjoint = prepared.prepare_adjoint(&gradient, 1.0, || false).unwrap();
    let result = adjoint
        .gradient_split_direction(direction(1.0), 0.0, source, zero_port, || false)
        .unwrap();
    assert!((result.derivative - reference).abs() < 2e-12);
    assert!(matches!(
        adjoint.gradient_split_direction(
            direction(1.0),
            0.0,
            source,
            |_, out| {
                out[0] = 0.0;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DirectionRefused { .. })
    ));
    assert!(matches!(
        adjoint.gradient_split_direction(
            direction(1.0),
            0.0,
            source,
            |point, out| {
                out.fill(0.0);
                out[0] = point.aft_samples as f64;
                true
            },
            || false
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        adjoint.gradient_split_direction(direction(1.0), 0.0, source, zero_port, || {
            polls += 1;
            polls > 13
        }),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert_eq!(
        result,
        adjoint
            .gradient_split_direction(direction(1.0), 0.0, source, zero_port, || false)
            .unwrap()
    );
    assert_eq!(problem.source_calls.get(), source_calls);
    assert_eq!(problem.tangent_calls.get(), tangent_calls);
}
