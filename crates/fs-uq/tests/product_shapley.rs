//! Physical conditional Gaussian attribution, analytic limits and paid-call resume.
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, ShapleyExecution, UqExecution,
    UqPlan, UqStatus,
};

fn bivariate(rows: usize, rho: f64) -> UqPlan {
    let mut plan = UqPlan::new("response", PropagationMethod::MonteCarlo, rows * 6)
        .with_parameter(ParameterUncertainty::gaussian("load", 100.0, 3.0, "N"))
        .with_parameter(ParameterUncertainty::gaussian("modulus", 20.0, 2.0, "Pa"))
        .with_correlation(CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, rho], vec![rho, 1.0]],
        });
    plan.seed = 83;
    plan
}

fn linear(x: &[f64]) -> Result<f64, &'static str> {
    Ok(2.0 * ((x[0] - 100.0) / 3.0) + (x[1] - 20.0) / 2.0)
}

#[test]
fn dependence_assigns_variance_to_a_physical_input_not_used_by_the_function() {
    for rho in [0.0, 0.6, 1.0, -1.0] {
        let mut run = ShapleyExecution::new(&bivariate(8192, rho)).unwrap();
        let report = run.advance(
            usize::MAX,
            || false,
            |x| Ok::<_, &str>((x[0] - 100.0) / 3.0),
        );
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(report.evaluations_attempted, 8192 * 6);
        let estimate = report.estimate.unwrap();
        // Y=X1 standardized. Var(Y|X1)=0; Var(Y|X2)=1-rho^2.
        // Thus Sh_1=1-rho^2/2 and Sh_2=rho^2/2, despite no direct X2 term.
        assert!((estimate.effects[0].index - (1.0 - rho * rho / 2.0)).abs() < 0.025);
        assert!((estimate.effects[1].index - rho * rho / 2.0).abs() < 0.025);
        assert!(
            (estimate
                .effects
                .iter()
                .map(|effect| effect.index)
                .sum::<f64>()
                - 1.0)
                .abs()
                < 1e-14
        );
        if rho.abs() == 1.0 {
            assert_eq!(estimate.effects[0].index, 0.5);
            assert_eq!(estimate.effects[1].index, 0.5);
        }
        assert_eq!(estimate.effects[0].parameter, "load");
        assert_eq!(estimate.effects[0].unit, "N");
        assert_eq!(estimate.effects[1].parameter, "modulus");
        assert_eq!(estimate.effects[1].unit, "Pa");
        assert!(matches!(report.evidence_color,
            fs_evidence::Color::Estimated { dispersion, .. } if dispersion == f64::INFINITY));
    }
}

#[test]
fn independent_nonlinear_interactions_are_shared_between_their_physical_inputs() {
    let rows = 8192;
    let mut plan = UqPlan::new("interaction", PropagationMethod::MonteCarlo, rows * 14)
        .with_correlation(CorrelationModel::Independent);
    for name in ["a", "b", "c"] {
        plan.parameters
            .push(ParameterUncertainty::gaussian(name, 0.0, 1.0, "1"));
    }
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let estimate = run
        .advance(
            usize::MAX,
            || false,
            |x| Ok::<_, &str>(2.0 * x[0] + x[1] * x[2]),
        )
        .estimate
        .unwrap();
    // Var(2*A)=4 and Var(B*C)=1; the interaction contributes 1/2 to each.
    for (actual, expected) in estimate.effects.iter().zip([0.8, 0.1, 0.1]) {
        assert!((actual.index - expected).abs() < 0.035);
    }
    assert!((estimate.base_output_std_dev.unwrap() - fs_math::det::sqrt(5.0)).abs() < 0.08);
    assert_eq!(estimate.conditional_variance_fractions.len(), 8);
    assert_eq!(estimate.conditional_variance_fractions[0], 1.0);
    assert_eq!(estimate.conditional_variance_fractions[7], 0.0);
}

#[test]
fn conditional_pairs_freeze_exact_physical_coordinates_and_match_base_mc_law() {
    let plan = bivariate(8, 0.6);
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let mut inputs = Vec::new();
    let report = run.advance(
        48,
        || false,
        |x| {
            inputs.push(x.to_vec());
            linear(x)
        },
    );
    assert_eq!(run.plan(), &plan);
    let mut mc_plan = plan.clone();
    mc_plan.budget_max_samples = 16;
    let mut mc = UqExecution::new(&mc_plan).unwrap();
    let mut reference = Vec::new();
    mc.advance(
        16,
        || false,
        |x| {
            reference.push(x.to_vec());
            linear(x)
        },
    );
    for (row, samples) in inputs.chunks_exact(6).enumerate() {
        assert_eq!(samples[0], reference[2 * row]);
        assert_eq!(samples[1], reference[2 * row + 1]);
        assert_eq!(samples[2][0].to_bits(), samples[3][0].to_bits()); // S={load}
        assert_eq!(samples[4][1].to_bits(), samples[5][1].to_bits()); // S={modulus}
    }
    // Independent reconstruction in physical output units, without the
    // production normalization or compensated-reduction helpers.
    let mut costs = [0.0; 4];
    for row in run.observations().chunks_exact(6) {
        for subset in 0..3 {
            costs[subset] += (row[2 * subset] - row[2 * subset + 1]).powi(2) / 16.0;
        }
    }
    let estimate = report.estimate.unwrap();
    for (fraction, cost) in estimate.conditional_variance_fractions.iter().zip(costs) {
        assert!((fraction - cost / costs[0]).abs() < 1e-13);
    }
    let first = (costs[0] - costs[1] + costs[2]) / (2.0 * costs[0]);
    assert!((estimate.effects[0].index - first).abs() < 1e-13);
    assert!((estimate.effects[1].index - (1.0 - first)).abs() < 1e-13);
}

#[test]
fn correlated_linear_effects_preserve_physical_permutations_and_output_affine_changes() {
    let plan = bivariate(8192, 0.6);
    let mut original = ShapleyExecution::new(&plan).unwrap();
    let estimate = original
        .advance(usize::MAX, || false, linear)
        .estimate
        .unwrap();
    // For Y=2*Z1+Z2, V=5+4*rho; conditional variances are 1-rho^2 and 4*(1-rho^2).
    let expected_first = 0.5 * (1.0 + 3.0 * (1.0 - 0.36) / (5.0 + 2.4));
    assert!((estimate.effects[0].index - expected_first).abs() < 0.025);
    let mut permuted = plan.clone();
    permuted.parameters.swap(0, 1);
    // This bivariate unit-diagonal symmetric matrix is unchanged by the swap.
    let mut swapped = ShapleyExecution::new(&permuted).unwrap();
    let other = swapped
        .advance(usize::MAX, || false, |x| linear(&[x[1], x[0]]))
        .estimate
        .unwrap();
    assert_eq!(other.effects[1].parameter, estimate.effects[0].parameter);
    assert!((other.effects[1].index - estimate.effects[0].index).abs() < 0.03);
    let mut affine = ShapleyExecution::new(&plan).unwrap();
    let transformed = affine
        .advance(
            usize::MAX,
            || false,
            |x| Ok::<_, &str>(1e6 - 8.0 * linear(x)?),
        )
        .estimate
        .unwrap();
    for (a, b) in transformed.effects.iter().zip(&estimate.effects) {
        assert!((a.index - b.index).abs() < 1e-9);
    }
}

#[test]
fn paid_half_pairs_survive_cancellation_interruption_chunks_and_clone() {
    let plan = bivariate(4, 0.6);
    let mut full = ShapleyExecution::new(&plan).unwrap();
    let mut expected_inputs = Vec::new();
    let expected = full.advance(
        24,
        || false,
        |x| {
            expected_inputs.push(x.to_vec());
            linear(x)
        },
    );
    let mut split = ShapleyExecution::new(&plan).unwrap();
    let mut paid_inputs = Vec::new();
    let half = split.advance(
        3,
        || false,
        |x| {
            paid_inputs.push(x.to_vec());
            linear(x)
        },
    );
    assert_eq!(half.evaluations_accepted, 3);
    assert_eq!(half.completed_rows, 0);
    assert!(half.estimate.is_none());
    assert_eq!(
        split.advance(
            0,
            || panic!("zero allowance"),
            |_| -> Result<f64, &str> { panic!("zero allowance") }
        ),
        half
    );
    let cancelled = split.advance(1, || true, |_| -> Result<f64, &str> { panic!("cancelled") });
    assert_eq!(cancelled.status, UqStatus::Cancelled);
    let paused = split.advance_interruptible(
        1,
        || false,
        |x| {
            assert_eq!(x, expected_inputs[3].as_slice());
            Ok::<_, &str>(None)
        },
    );
    assert_eq!(paused.evaluations_attempted, 3);
    let mut resumed = split.clone();
    for count in [1, 4, 99] {
        resumed.advance(
            count,
            || false,
            |x| {
                paid_inputs.push(x.to_vec());
                linear(x)
            },
        );
    }
    assert_eq!(paid_inputs, expected_inputs);
    assert_eq!(resumed.report(), expected);
    assert_eq!(resumed.observations(), full.observations());
    assert_eq!(split.observations().len(), 3);
    assert_eq!(
        resumed.advance(
            1,
            || panic!("complete"),
            |_| -> Result<f64, &str> { panic!("complete") }
        ),
        expected
    );
}

#[test]
fn malformed_laws_and_incomplete_design_budgets_refuse_before_physical_work() {
    let mut plan = bivariate(2, 0.6);
    for budget in [0, 6, 11, 13, 1_000_002] {
        plan.budget_max_samples = budget;
        assert!(ShapleyExecution::new(&plan).is_err());
    }
    plan.budget_max_samples = 12;
    for correlation in [
        CorrelationModel::Unknown,
        CorrelationModel::Correlated {
            matrix: vec![vec![1.0, 0.6], vec![0.6, 1.0]],
        },
        CorrelationModel::JointGaussian {
            matrix: vec![vec![1.0, 0.6]],
        },
    ] {
        plan.correlation = correlation;
        assert!(ShapleyExecution::new(&plan).is_err());
    }
    plan.correlation = CorrelationModel::Independent;
    for parameter in [
        ParameterUncertainty::uniform("load", 1.0, 2.0, "N"),
        ParameterUncertainty::gaussian("load", 1.0, 0.0, "N"),
        ParameterUncertainty::interval("load", 1.0, 2.0, "N"),
    ] {
        plan.parameters[0] = parameter;
        assert!(ShapleyExecution::new(&plan).is_err());
    }
    plan = bivariate(2, 0.6);
    plan.method = PropagationMethod::PolynomialChaos;
    assert!(ShapleyExecution::new(&plan).is_err());
    let mut maximum = UqPlan::new("maximum", PropagationMethod::MonteCarlo, 1020)
        .with_correlation(CorrelationModel::Independent);
    for i in 0..8 {
        maximum.parameters.push(ParameterUncertainty::gaussian(
            format!("x{i}"),
            0.0,
            1.0,
            "1",
        ));
    }
    assert!(ShapleyExecution::new(&maximum).is_ok());
    maximum
        .parameters
        .push(ParameterUncertainty::gaussian("ninth", 0.0, 1.0, "1"));
    assert!(ShapleyExecution::new(&maximum).is_err());
}

#[test]
fn finite_sample_effects_are_unclipped_and_undefined_variance_is_not_zero_effect() {
    // Small fixed independent designs can assign a negative effect to a dummy
    // input. Check the exact algebra for each design and retain such a value.
    let mut found_negative = false;
    for seed in 0..32 {
        let mut plan = bivariate(2, 0.0);
        plan.seed = seed;
        let mut run = ShapleyExecution::new(&plan).unwrap();
        let estimate = run
            .advance(12, || false, |x| Ok::<_, &str>(x[0]))
            .estimate
            .unwrap();
        let cost = &estimate.conditional_variance_fractions;
        assert_eq!(estimate.effects[1].index, 0.5 * (1.0 - cost[2]));
        found_negative |= estimate.effects[1].index < 0.0;
    }
    assert!(
        found_negative,
        "sampling noise is not clipped into a physical guarantee"
    );
    let mut constant = ShapleyExecution::new(&bivariate(2, 1.0)).unwrap();
    let report = constant.advance(12, || false, |_| Ok::<_, &str>(f64::MAX));
    assert_eq!(report.status, UqStatus::Complete);
    assert!(report.estimate.is_none() && report.unavailable_reason.is_some());
    let mut extremes = ShapleyExecution::new(&bivariate(128, 0.6)).unwrap();
    let report = extremes.advance(
        usize::MAX,
        || false,
        |x| Ok::<_, &str>(if x[0] >= 100.0 { f64::MAX } else { -f64::MAX }),
    );
    assert_eq!(report.status, UqStatus::Complete);
    assert!(
        report
            .estimate
            .unwrap()
            .effects
            .iter()
            .all(|e| e.index.is_finite())
    );
    for output in [
        Err("physical solver refused"),
        Ok(f64::NAN),
        Ok(f64::INFINITY),
    ] {
        let mut run = ShapleyExecution::new(&bivariate(2, 0.6)).unwrap();
        run.advance(3, || false, linear);
        let failed = run.advance(1, || false, |_| output);
        assert_eq!(failed.status, UqStatus::Refused);
        assert_eq!(failed.evaluations_attempted, 4);
        assert_eq!(failed.evaluations_accepted, 3);
        assert!(failed.estimate.is_none());
        assert_eq!(
            run.advance(
                1,
                || panic!("terminal"),
                |_| -> Result<f64, &str> { panic!("terminal") }
            ),
            failed
        );
    }
}
