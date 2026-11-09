use super::*;

fn controls(dt: f64) -> OutflowMidpointSolve {
    OutflowMidpointSolve {
        duration_s: dt,
        maximum_outgoing_fraction: 0.5,
        mass_tolerance_kg: 1e-14,
        maximum_iterations: 64,
    }
}

#[test]
fn linear_outflow_matches_the_independent_midpoint_solution_with_or_without_a_tangent() {
    for tangent in [Some(0.7), None] {
        let result = solve(2.0, 1.4, controls(0.2), |m, _| {
            Ok::<_, ()>(OutflowRate { mass_flow_kg_s: 0.7 * m, derivative_per_s: tangent })
        }, || false).unwrap();
        let expected = 2.0 * (1.0 - 0.5 * 0.7 * 0.2) / (1.0 + 0.5 * 0.7 * 0.2);
        assert!((result.remaining_mass_kg - expected).abs() < 1e-14);
        assert_eq!(result.expelled_mass_kg, 2.0 - result.remaining_mass_kg);
        assert_eq!(result.midpoint_mass_kg, 2.0 - 0.5 * result.expelled_mass_kg);
        assert!(result.mass_residual_kg.abs() <= controls(0.2).mass_tolerance_kg);
    }
}

#[test]
fn smooth_nonlinear_outflow_converges_at_second_order_over_three_refinements() {
    // dm/dt = -m^(3/2) has m(t) = (1+t/2)^(-2) for m(0)=1.
    let duration = 0.3;
    let exact = 1.0_f64 / (1.0_f64 + 0.5 * duration).powi(2);
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut mass: f64 = 1.0;
        for _ in 0..steps {
            let result = solve(mass, mass * mass.sqrt(), controls(duration / steps as f64),
                |m, _| Ok::<_, ()>(OutflowRate {
                    mass_flow_kg_s: m * m.sqrt(), derivative_per_s: Some(1.5 * m.sqrt()),
                }), || false).unwrap();
            mass = result.remaining_mass_kg;
        }
        errors.push((mass - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.9 && pair[0] / pair[1] < 4.1, "{errors:?}");
    }
}

#[test]
fn an_inventory_limit_returns_the_actual_midpoint_retry_duration_without_clipping() {
    let mut c = controls(1.0);
    c.maximum_outgoing_fraction = 0.1;
    let mut samples = Vec::new();
    let error = solve(1.0, 1.0, c, |m, _| {
        samples.push(m);
        Ok::<_, ()>(OutflowRate { mass_flow_kg_s: m, derivative_per_s: Some(1.0) })
    }, || false).unwrap_err();
    let OutflowMidpointError::InventoryLimit { maximum_duration_s } = error else {
        panic!("{error:?}");
    };
    assert_eq!(samples, [0.95]);
    assert!((maximum_duration_s - (1.0 - 0.9) / 0.95).abs() < 1e-16);
}

#[test]
fn exact_rest_skips_trials_and_cancellation_prevents_final_publication() {
    let rest = solve(1.0, 0.0, controls(1.0), |_, _| -> Result<OutflowRate, ()> {
        panic!("a resting tank has no mass trial");
    }, || false).unwrap();
    assert_eq!(rest.remaining_mass_kg, 1.0);
    assert_eq!(rest.expelled_mass_kg, 0.0);
    assert_eq!(rest.iterations, 0);
    let mut polls = 0;
    let stopped = solve(1.0, 0.1, controls(1.0), |_, _| Ok::<_, ()>(OutflowRate {
        mass_flow_kg_s: 0.1, derivative_per_s: Some(0.0),
    }), || { polls += 1; polls == 3 });
    assert_eq!(stopped, Err(OutflowMidpointError::Cancelled));
    assert_eq!(polls, 3);
}

#[test]
fn actual_rate_work_source_and_unrepresentable_inventory_refusals_are_preserved() {
    let mut c = controls(0.2);
    c.maximum_iterations = 1;
    let stopped = solve(1.0, 1.0, c, |m, _| Ok::<_, ()>(OutflowRate {
        mass_flow_kg_s: m * m, derivative_per_s: Some(2.0 * m),
    }), || false);
    assert!(matches!(stopped, Err(OutflowMidpointError::NoConvergence { iterations: 1, .. })));
    assert_eq!(solve(1.0, 0.1, controls(0.2), |_, _| Err::<OutflowRate, _>(7), || false),
        Err(OutflowMidpointError::Source(7)));
    assert!(matches!(solve(1.0, 1e-30, controls(1.0), |_, _| -> Result<OutflowRate, ()> {
        panic!("an unrepresentable mass change must refuse before chemistry");
    }, || false), Err(OutflowMidpointError::Invalid(_))));
    assert!(matches!(solve(1.0, 0.1, controls(0.2), |_, _| Ok::<_, ()>(OutflowRate {
        mass_flow_kg_s: f64::NAN, derivative_per_s: None,
    }), || false), Err(OutflowMidpointError::Invalid(_))));
}
