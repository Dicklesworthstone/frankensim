use super::*;

fn zero_island(_: f64, _: &[f64], out: &mut [f64]) -> bool {
    out.fill(0.0);
    true
}

#[test]
fn sampled_direct_partial_composes_with_state_response_and_independent_units() {
    for converted in [false, true] {
        let mut problem = Compliance::new();
        if converted {
            problem.states = [1e6, 1e-5];
            problem.equations = [1e-3, 1e4];
        }
        let parameter_unit = if converted { 1e-7 } else { 1.0 };
        let objective_unit = if converted { 1e9 } else { 1.0 };
        let parameter = 0.7;
        let mut prepared = problem.prepare(32);
        let baseline = problem.counts();
        let objective = prepared
            .mean_objective(
                objective_controls(objective_unit),
                |_, x, value, gradient| {
                    problem.square(x, parameter * objective_unit, value, gradient)
                },
                || false,
            )
            .unwrap();
        let mut adjoint = prepared
            .prepare_adjoint(&objective.gradient, objective_unit, || false)
            .unwrap();
        let direct_calls = Cell::new(0);
        let source_calls = Cell::new(0);
        // In a local family x(p) = p*x_base/p_base, q(p) = p*mean(y(p)^2),
        // q' = 3*mean(y_base^2). The explicit and state chains are both nonzero.
        let direct = |_: f64, x: &[f64], out: &mut f64| {
            direct_calls.set(direct_calls.get() + 1);
            let y = x[0] / problem.states[0] + 0.7 * x[1] / problem.states[1];
            *out = objective_unit * y * y / parameter_unit;
            true
        };
        let island = |_: f64, x: &[f64], out: &mut [f64]| {
            source_calls.set(source_calls.get() + 1);
            for j in 0..2 {
                out[j] =
                    problem.equations[j] * x[j] / problem.states[j] / parameter / parameter_unit;
            }
            true
        };
        let mut chart = ForcedSensitivityDirection {
            parameter: parameter * parameter_unit,
            parameter_scale: parameter_unit,
            // This static compliance has no port-frequency dependence; its
            // physical source coefficients are held in the fixed-phase chart.
            angular_frequency_derivative: 0.3 / parameter_unit,
        };
        let split = adjoint
            .gradient_mean_split_direction(
                chart,
                direct,
                island,
                |_, out| {
                    out.fill(0.0);
                    true
                },
                || false,
            )
            .unwrap();
        let expected = 3.0 * objective.value / parameter / objective_unit;
        close(
            split.derivative * parameter_unit / objective_unit,
            expected,
            1e-12,
        );
        assert_eq!(direct_calls.get(), 97);
        assert_eq!(source_calls.get(), 97);
        assert_eq!(problem.counts(), baseline);
        chart.angular_frequency_derivative = 0.0;
        let island_only = adjoint
            .gradient_mean_island_direction(chart, direct, island, || false)
            .unwrap();
        assert_eq!(island_only, split);
    }
}

#[test]
fn direct_and_total_gradient_aliasing_refuse_independently_and_recover() {
    let problem = Compliance::new();
    let mut prepared = problem.prepare(32);
    let mut covector = vec![0.0; 14];
    covector[0] = 1.0;
    let mut adjoint = prepared.prepare_adjoint(&covector, 1.0, || false).unwrap();
    let baseline = problem.counts();
    let chart = ForcedSensitivityDirection {
        parameter: 0.0,
        parameter_scale: 1.0,
        angular_frequency_derivative: 0.0,
    };
    let alias = |time: f64| det::cos(32.0 * OMEGA * time);
    assert!(matches!(
        adjoint.gradient_mean_island_direction(
            chart,
            |t, _, out| {
                *out = alias(t);
                true
            },
            zero_island,
            || false,
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    // Each separately has a 0.75*limit grid change, but together their physical
    // derivative changes by 1.5*limit. Reusing the fine q_p on both grids would
    // incorrectly admit this direction.
    assert!(matches!(
        adjoint.gradient_mean_island_direction(
            chart,
            |t, _, out| {
                *out = 7.5e-11 * alias(t);
                true
            },
            |t, _, out| {
                out.copy_from_slice(&[7.5e-11 * alias(t), 0.0]);
                true
            },
            || false,
        ),
        Err(ForcedSensitivityError::DerivativeUnresolved { .. })
    ));
    assert!(matches!(
        adjoint.gradient_mean_island_direction(chart, |_, _, _| true, zero_island, || false),
        Err(ForcedSensitivityError::ObjectiveRefused {
            sample: 0,
            samples: 32
        })
    ));
    let calls = Cell::new(0);
    assert!(matches!(
        adjoint.gradient_mean_island_direction(
            chart,
            |_, _, out| {
                calls.set(calls.get() + 1);
                *out = 0.3;
                true
            },
            zero_island,
            || calls.get() == 3,
        ),
        Err(ForcedSensitivityError::Cancelled)
    ));
    assert_eq!(calls.get(), 3);
    let run = |adjoint: &mut fs_orbit::sensitivity::PreparedForcedAdjoint<'_>| {
        adjoint
            .gradient_mean_island_direction(
                chart,
                |_, x, out| {
                    *out = 0.3 + x[0];
                    true
                },
                zero_island,
                || false,
            )
            .unwrap()
    };
    let recovered = run(&mut adjoint);
    close(recovered.derivative, 0.6, 1e-14);
    assert_eq!(recovered, run(&mut adjoint));
    assert_eq!(problem.counts(), baseline);
}

#[test]
fn explicit_mean_retains_representable_extreme_dc_without_a_raw_sample_sum() {
    let problem = Compliance::new();
    let mut prepared = problem.prepare(32);
    let covector = vec![0.0; 14];
    let mut adjoint = prepared
        .prepare_adjoint(&covector, 1e308, || false)
        .unwrap();
    let result = adjoint
        .gradient_mean_island_direction(
            ForcedSensitivityDirection {
                parameter: 0.0,
                parameter_scale: 1.0,
                angular_frequency_derivative: 0.0,
            },
            |t, _, out| {
                *out = 1e308 * (1.0 + 0.1 * det::cos(OMEGA * t));
                true
            },
            zero_island,
            || false,
        )
        .unwrap();
    close(result.derivative / 1e308, 1.0, 1e-14);
    assert!(result.gradient_sampling_defect < 1e-14);
}
