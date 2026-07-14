//! G0/G3/G4/G5 battery for checked, cancellation-correct data assimilation
//! (addendum Proposal 11).

use fs_assimilate::{
    AssimError, AssimilatedPosterior, Belief, Color, MAX_DENSE_OBSERVATIONS, MAX_DENSE_STATE_DIM,
    MAX_DENSE_UPDATE_CUBIC_WORK, Observation, assimilate as assimilate_with_cx,
    assimilate_all as assimilate_all_with_cx, assimilate_colored as assimilate_colored_with_cx,
    assimilate_colored_with_shared_poll_quota, misfit as misfit_with_cx, point_sensor,
    scan_observation,
};
use fs_evidence::{MAX_COLOR_IDENTITY_BYTES, color_leaf_identity_reason};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0xA5_51_11_A7_E,
    kernel_id: 0xA551,
    tile: 0,
    iteration: 0,
};

fn with_stream_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    mode: ExecMode,
    stream: StreamKey,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(gate, arena, stream, budget, mode);
        f(&cx)
    })
}

fn with_configured_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    mode: ExecMode,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    with_stream_cx(gate, budget, mode, TEST_STREAM, f)
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_configured_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        ExecMode::Deterministic,
        f,
    )
}

fn belief(mean: Vec<f64>, covariance: Vec<Vec<f64>>) -> Result<Belief, AssimError> {
    with_cx(|cx| Belief::new(mean, covariance, cx))
}

fn diagonal(means: Vec<f64>, variances: &[f64]) -> Result<Belief, AssimError> {
    with_cx(|cx| Belief::diagonal(means, variances, cx))
}

fn validate(belief: &Belief) -> Result<(), AssimError> {
    with_cx(|cx| belief.validate(cx))
}

fn assimilate(prior: &Belief, observation: &Observation) -> Result<Belief, AssimError> {
    with_cx(|cx| assimilate_with_cx(prior, observation, cx))
}

fn assimilate_all(prior: &Belief, observations: &[Observation]) -> Result<Belief, AssimError> {
    with_cx(|cx| assimilate_all_with_cx(prior, observations, cx))
}

fn misfit(prior: &Belief, observations: &[Observation]) -> Result<f64, AssimError> {
    with_cx(|cx| misfit_with_cx(prior, observations, cx))
}

fn assimilate_colored(
    prior: &Belief,
    observations: &[Observation],
    regime_param: &str,
    regime_lo: f64,
    regime_hi: f64,
) -> Result<AssimilatedPosterior, AssimError> {
    with_cx(|cx| {
        assimilate_colored_with_cx(prior, observations, regime_param, regime_lo, regime_hi, cx)
    })
}

fn scalar(mean: f64, variance: f64) -> Belief {
    Belief::scalar(mean, variance).expect("valid scalar fixture")
}

fn sensor(component: usize, dim: usize, value: f64, noise: f64, instrument: &str) -> Observation {
    point_sensor(component, dim, value, noise, instrument).expect("valid sensor fixture")
}

fn estimator_identity(posterior: &AssimilatedPosterior) -> &str {
    match posterior.color() {
        Color::Estimated { estimator, .. } => estimator,
        other => panic!("expected estimated candidate, got {other:?}"),
    }
}

fn assert_same_non_identity_posterior(left: &AssimilatedPosterior, right: &AssimilatedPosterior) {
    assert_eq!(left.belief(), right.belief());
    assert_eq!(left.regime(), right.regime());
    assert_eq!(
        left.misfit_before().to_bits(),
        right.misfit_before().to_bits()
    );
    assert_eq!(
        left.misfit_after().to_bits(),
        right.misfit_after().to_bits()
    );
    match (left.color(), right.color()) {
        (
            Color::Estimated {
                dispersion: left_dispersion,
                ..
            },
            Color::Estimated {
                dispersion: right_dispersion,
                ..
            },
        ) => assert_eq!(left_dispersion.to_bits(), right_dispersion.to_bits()),
        other => panic!("expected matching Estimated colors, got {other:?}"),
    }
}

fn assert_initial_cancelled<T: core::fmt::Debug>(result: Result<T, AssimError>) {
    assert!(matches!(
        result,
        Err(AssimError::Cancelled {
            phase: "initial",
            completed: 0,
            planned,
        }) if planned > 0
    ));
}

#[test]
fn a_measurement_shifts_the_mean_and_shrinks_the_variance() {
    let prior = scalar(0.0, 10.0);
    let obs = sensor(0, 1, 5.0, 1.0, "gauge-A");
    let post = assimilate(&prior, &obs).unwrap();
    // K = 10/11 -> mean ~ 4.545.
    assert!((post.component_mean(0).unwrap() - 5.0 * 10.0 / 11.0).abs() < 1e-9);
    // posterior variance is below BOTH the prior (10) and the measurement (1).
    assert!(post.variance(0).unwrap() < 1.0 && post.variance(0).unwrap() < 10.0);
}

#[test]
fn assimilation_reduces_the_model_data_misfit() {
    let prior = diagonal(vec![0.0, 0.0], &[10.0, 10.0]).unwrap();
    let obs = vec![
        sensor(0, 2, 5.0, 1.0, "gauge-x"),
        sensor(1, 2, -3.0, 1.0, "gauge-y"),
    ];
    let before = misfit(&prior, &obs).unwrap();
    let post = assimilate_all(&prior, &obs).unwrap();
    let after = misfit(&post, &obs).unwrap();
    assert!(after < before, "misfit {after} !< {before}");
    assert!(before > 30.0); // 25 + 9
}

#[test]
fn a_scan_observation_carries_more_noise_than_a_point_sensor() {
    // point sensors have NO registration problem (the R8 fallback).
    let point = sensor(0, 2, 1.0, 0.04, "thermocouple");
    let scan = scan_observation(vec![1.0, 0.0], 1.0, 0.04, 0.25, "CT-scan").unwrap();
    assert!((point.noise_var() - 0.04).abs() < 1e-12);
    // the scan's noise adds the registration variance on top.
    assert!((scan.noise_var() - (0.04 + 0.25)).abs() < 1e-12);
    assert!(scan.noise_var() > point.noise_var());
}

#[test]
fn the_linear_gaussian_posterior_is_order_independent() {
    let prior = scalar(0.0, 10.0);
    let o1 = sensor(0, 1, 5.0, 1.0, "a");
    let o2 = sensor(0, 1, 3.0, 2.0, "b");
    let fwd = assimilate_all(&prior, &[o1.clone(), o2.clone()]).unwrap();
    let rev = assimilate_all(&prior, &[o2.clone(), o1.clone()]).unwrap();
    assert_eq!(fwd, rev);
    // sequential fusion tightens: posterior variance below the single-obs one.
    assert!(
        fwd.variance(0).unwrap()
            < assimilate(&prior, &sensor(0, 1, 5.0, 1.0, "a"))
                .unwrap()
                .variance(0)
                .unwrap()
    );

    // The content identity canonicalizes observation order too.
    let colored_fwd =
        assimilate_colored(&prior, &[o1.clone(), o2.clone()], "Re", 1.0, 2.0).unwrap();
    let colored_rev = assimilate_colored(&prior, &[o2, o1], "Re", 1.0, 2.0).unwrap();
    assert_eq!(
        estimator_identity(&colored_fwd),
        estimator_identity(&colored_rev)
    );
    assert_eq!(colored_fwd, colored_rev);
}

#[test]
fn the_assimilated_posterior_is_an_honest_bounded_estimate() {
    let prior = scalar(20.0, 4.0);
    let obs = vec![
        sensor(0, 1, 22.0, 0.25, "thermocouple-7"),
        sensor(0, 1, 21.5, 0.25, "thermocouple-7"),
    ];
    let out = assimilate_colored(&prior, &obs, "Re", 1e5, 3e5).unwrap();
    assert!(out.misfit_after() <= out.misfit_before());
    assert_eq!(out.regime().bound("Re"), Some((1e5, 3e5)));
    match out.color() {
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            assert!(estimator.starts_with("assimilation-candidate:v3:"));
            assert!(estimator.len() <= MAX_COLOR_IDENTITY_BYTES);
            assert_eq!(color_leaf_identity_reason(estimator), None);
            assert!(dispersion.is_infinite() && dispersion.is_sign_positive());
        }
        other => panic!("expected estimated candidate, got {other:?}"),
    }
}

#[test]
fn belief_construction_rejects_empty_ragged_and_nonfinite_state() {
    assert_eq!(belief(vec![], vec![]), Err(AssimError::EmptyBelief));
    assert_eq!(diagonal(vec![], &[]), Err(AssimError::EmptyBelief));
    assert_eq!(
        diagonal(vec![0.0, 1.0], &[1.0]),
        Err(AssimError::DiagonalDimensionMismatch {
            means: 2,
            variances: 1,
        })
    );
    assert_eq!(
        belief(vec![0.0], vec![]),
        Err(AssimError::CovarianceDimensionMismatch { state: 1, rows: 0 })
    );
    assert_eq!(
        belief(vec![0.0, 0.0], vec![vec![1.0], vec![0.0, 1.0]]),
        Err(AssimError::CovarianceRowDimensionMismatch {
            row: 0,
            expected: 2,
            actual: 1,
        })
    );
    assert_eq!(
        Belief::scalar(f64::NAN, 1.0),
        Err(AssimError::NonFiniteMean { index: 0 })
    );
    assert_eq!(
        Belief::scalar(0.0, f64::INFINITY),
        Err(AssimError::NonFiniteCovariance { row: 0, column: 0 })
    );
}

#[test]
fn g4_empty_diagonal_preflight_precedes_cancellation_and_zero_poll_budget() {
    let cancelled = CancelGate::new();
    cancelled.request();
    let pre_cancelled = with_configured_cx(
        &cancelled,
        Budget::INFINITE,
        ExecMode::Deterministic,
        |cx| Belief::diagonal(Vec::new(), &[], cx),
    );
    assert_eq!(pre_cancelled, Err(AssimError::EmptyBelief));

    let healthy = CancelGate::new();
    let zero_quota = with_configured_cx(
        &healthy,
        Budget::INFINITE.with_poll_quota(0),
        ExecMode::Deterministic,
        |cx| Belief::diagonal(Vec::new(), &[], cx),
    );
    assert_eq!(zero_quota, Err(AssimError::EmptyBelief));
}

#[test]
fn belief_construction_enforces_covariance_semantics() {
    assert_eq!(
        Belief::scalar(0.0, -1.0),
        Err(AssimError::NegativeVariance { index: 0 })
    );
    assert_eq!(
        belief(vec![0.0, 0.0], vec![vec![1.0, 0.5], vec![0.4, 1.0]]),
        Err(AssimError::NonSymmetricCovariance { row: 0, column: 1 })
    );
    assert_eq!(
        belief(vec![0.0, 0.0], vec![vec![1.0, 2.0], vec![2.0, 1.0]]),
        Err(AssimError::CovarianceNotPositiveSemidefinite)
    );
    // A global magnitude must not hide invalid coupling to a zero-variance
    // component in the semidefinite tolerance.
    assert_eq!(
        belief(
            vec![0.0, 0.0, 0.0],
            vec![
                vec![1e300, 0.0, 0.0],
                vec![0.0, 0.0, 1e-20],
                vec![0.0, 1e-20, 1.0],
            ],
        ),
        Err(AssimError::CovarianceNotPositiveSemidefinite)
    );
    // Unit scaling does not reject a valid high-dynamic-range covariance.
    belief(vec![0.0, 0.0], vec![vec![1e300, 0.0], vec![0.0, 1e-300]]).unwrap();
    // Singular positive-semidefinite covariance is valid.
    let semidefinite = belief(vec![0.0, 0.0], vec![vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
    assert_eq!(semidefinite.variance(1), Ok(1.0));
}

#[test]
fn dense_state_dimension_is_bounded_before_quadratic_allocation() {
    let oversized = MAX_DENSE_STATE_DIM + 1;
    let expected = AssimError::StateDimensionLimit {
        dim: oversized,
        max: MAX_DENSE_STATE_DIM,
    };
    assert_eq!(
        diagonal(vec![0.0; oversized], &vec![1.0; oversized]),
        Err(expected.clone())
    );
    assert_eq!(
        point_sensor(0, oversized, 0.0, 1.0, "bounded-sensor"),
        Err(expected.clone())
    );
    assert_eq!(
        Observation::new(vec![1.0; oversized], 0.0, 1.0, "bounded-sensor"),
        Err(expected)
    );
}

#[test]
fn psd_admission_never_tolerance_clamps_negative_curvature() {
    // Eigenvalues are 2 + 1e-14 and -1e-14. The old LDL guard replaced the
    // negative pivot with zero because it fell inside a scale-relative
    // tolerance, laundering a mathematically indefinite matrix into Belief.
    let correlation = 1.0 + 1e-14;
    assert_eq!(
        belief(
            vec![0.0, 0.0],
            vec![vec![1.0, correlation], vec![correlation, 1.0]],
        ),
        Err(AssimError::CovarianceNotPositiveSemidefinite)
    );

    // The determinant is exactly -EPSILON^2 even though both diagonal
    // products round to 1.0 in ordinary f64 arithmetic. Admission compares
    // the underlying binary-rational products exactly before correlation
    // scaling, so this sub-ulp negative 2x2 minor is still refused.
    assert_eq!(
        belief(
            vec![0.0, 0.0],
            vec![vec![1.0 + f64::EPSILON, 1.0], vec![1.0, 1.0 - f64::EPSILON],],
        ),
        Err(AssimError::CovarianceNotPositiveSemidefinite)
    );
}

#[test]
fn indexed_access_and_point_sensor_bounds_are_checked() {
    let prior = scalar(0.0, 1.0);
    assert_eq!(
        prior.component_mean(1),
        Err(AssimError::ComponentOutOfRange {
            component: 1,
            dim: 1,
        })
    );
    assert_eq!(
        prior.variance(usize::MAX),
        Err(AssimError::ComponentOutOfRange {
            component: usize::MAX,
            dim: 1,
        })
    );
    assert_eq!(
        point_sensor(0, 0, 1.0, 1.0, "x"),
        Err(AssimError::EmptyStateDimension)
    );
    assert_eq!(
        point_sensor(2, 2, 1.0, 1.0, "x"),
        Err(AssimError::ComponentOutOfRange {
            component: 2,
            dim: 2,
        })
    );
}

#[test]
fn observation_construction_rejects_malformed_numeric_inputs() {
    assert_eq!(
        Observation::new(vec![], 1.0, 1.0, "x"),
        Err(AssimError::EmptyObservationOperator)
    );
    assert_eq!(
        Observation::new(vec![0.0, -0.0], 1.0, 1.0, "x"),
        Err(AssimError::ZeroObservationOperator)
    );
    assert_eq!(
        Observation::new(vec![f64::NAN], 1.0, 1.0, "x"),
        Err(AssimError::NonFiniteObservationOperator { index: 0 })
    );
    assert_eq!(
        Observation::new(vec![1.0], f64::INFINITY, 1.0, "x"),
        Err(AssimError::NonFiniteObservationValue)
    );
    assert_eq!(
        Observation::new(vec![1.0], 1.0, 0.0, "x"),
        Err(AssimError::NonPositiveNoise)
    );
    assert_eq!(
        Observation::new(vec![1.0], 1.0, f64::NAN, "x"),
        Err(AssimError::NonFiniteNoise)
    );
}

#[test]
fn observation_construction_rejects_unusable_instrument_identities() {
    assert_eq!(
        Observation::new(vec![1.0], 1.0, 1.0, "  "),
        Err(AssimError::EmptyInstrument)
    );
    assert_eq!(
        Observation::new(vec![1.0], 1.0, 1.0, "unknown"),
        Err(AssimError::InvalidIdentity {
            field: "instrument",
            reason: "placeholder",
        })
    );
    assert_eq!(
        Observation::new(vec![1.0], 1.0, 1.0, "derived:v2:forged"),
        Err(AssimError::InvalidIdentity {
            field: "instrument",
            reason: "derived-identity-requires-lineage",
        })
    );
    assert_eq!(
        Observation::new(
            vec![1.0],
            1.0,
            1.0,
            "x".repeat(MAX_COLOR_IDENTITY_BYTES + 1)
        ),
        Err(AssimError::InvalidIdentity {
            field: "instrument",
            reason: "too-long",
        })
    );
}

#[test]
fn scan_noise_validation_is_fail_closed() {
    assert_eq!(
        scan_observation(vec![1.0], 1.0, 0.0, 0.1, "scan"),
        Err(AssimError::NonPositiveNoise)
    );
    assert_eq!(
        scan_observation(vec![1.0], 1.0, 1.0, -0.1, "scan"),
        Err(AssimError::NegativeRegistrationVariance)
    );
    assert_eq!(
        scan_observation(vec![1.0], 1.0, 1.0, f64::INFINITY, "scan"),
        Err(AssimError::NonFiniteRegistrationVariance)
    );
    assert_eq!(
        scan_observation(vec![1.0], 1.0, f64::MAX, f64::MAX, "scan"),
        Err(AssimError::NonFiniteNoise)
    );
}

#[test]
fn assimilation_rejects_dimension_mismatch_and_empty_aggregates() {
    let prior = scalar(0.0, 1.0);
    // operator length != state dim.
    let bad_dim = sensor(0, 2, 1.0, 1.0, "x");
    assert_eq!(
        assimilate(&prior, &bad_dim),
        Err(AssimError::DimMismatch {
            state: 1,
            operator: 2,
        })
    );
    assert_eq!(misfit(&prior, &[]), Err(AssimError::EmptyObservations));
    assert_eq!(
        assimilate_all(&prior, &[]),
        Err(AssimError::EmptyObservations)
    );
    assert_eq!(
        assimilate_colored(&prior, &[], "Re", 1.0, 2.0),
        Err(AssimError::EmptyObservations)
    );
}

#[test]
fn dense_aggregate_work_is_bounded_before_canonicalization_or_updates() {
    let scalar_prior = scalar(0.0, 1.0);
    let scalar_observation = sensor(0, 1, 0.0, 1.0, "bounded-count");
    let too_many = vec![scalar_observation; MAX_DENSE_OBSERVATIONS + 1];
    let expected_count = AssimError::ObservationCountLimit {
        count: too_many.len(),
        max: MAX_DENSE_OBSERVATIONS,
    };
    assert_eq!(
        misfit(&scalar_prior, &too_many),
        Err(expected_count.clone())
    );
    assert_eq!(
        assimilate_all(&scalar_prior, &too_many),
        Err(expected_count.clone())
    );
    assert_eq!(
        assimilate_colored(&scalar_prior, &too_many, "Re", 1.0, 2.0),
        Err(expected_count)
    );

    let dimension = MAX_DENSE_STATE_DIM;
    let dense_prior = diagonal(vec![0.0; dimension], &vec![1.0; dimension])
        .expect("maximum admitted dense state");
    let dense_observation = sensor(0, dimension, 0.0, 1.0, "bounded-work");
    let observations = vec![dense_observation; 5];
    let requested = (dimension as u128).pow(3) * observations.len() as u128;
    assert!(requested > MAX_DENSE_UPDATE_CUBIC_WORK);
    let expected_work = AssimError::AssimilationWorkLimit {
        requested,
        max: MAX_DENSE_UPDATE_CUBIC_WORK,
    };
    assert_eq!(
        assimilate_all(&dense_prior, &observations),
        Err(expected_work.clone())
    );
    assert_eq!(
        assimilate_colored(&dense_prior, &observations, "Re", 1.0, 2.0),
        Err(expected_work)
    );
}

#[test]
fn regime_validation_rejects_empty_malformed_and_invalid_bounds() {
    let prior = scalar(0.0, 1.0);
    let obs = [sensor(0, 1, 1.0, 1.0, "x")];
    assert_eq!(
        assimilate_colored(&prior, &obs, " ", 1.0, 2.0),
        Err(AssimError::EmptyRegime)
    );
    assert_eq!(
        assimilate_colored(&prior, &obs, "derived:v2:axis", 1.0, 2.0),
        Err(AssimError::InvalidIdentity {
            field: "regime_param",
            reason: "derived-identity-requires-lineage",
        })
    );
    assert_eq!(
        assimilate_colored(&prior, &obs, "Re", f64::NAN, 2.0),
        Err(AssimError::NonFiniteRegimeBounds)
    );
    assert_eq!(
        assimilate_colored(&prior, &obs, "Re", 3.0, 2.0),
        Err(AssimError::InvertedRegimeBounds)
    );
}

#[test]
fn candidate_identity_is_unambiguous_and_binds_multiplicity() {
    let prior = scalar(0.0, 2.0);
    let joined_left = vec![sensor(0, 1, 1.0, 1.0, "a+b"), sensor(0, 1, 2.0, 1.0, "c")];
    let joined_right = vec![sensor(0, 1, 1.0, 1.0, "a"), sensor(0, 1, 2.0, 1.0, "b+c")];
    let left = assimilate_colored(&prior, &joined_left, "Re", 1.0, 2.0).unwrap();
    let right = assimilate_colored(&prior, &joined_right, "Re", 1.0, 2.0).unwrap();
    assert_ne!(estimator_identity(&left), estimator_identity(&right));

    let once = assimilate_colored(&prior, &joined_left[..1], "Re", 1.0, 2.0).unwrap();
    let twice = assimilate_colored(
        &prior,
        &[joined_left[0].clone(), joined_left[0].clone()],
        "Re",
        1.0,
        2.0,
    )
    .unwrap();
    assert_ne!(estimator_identity(&once), estimator_identity(&twice));
}

#[test]
fn signed_zero_is_canonicalized_before_hashing_and_computation() {
    let negative = belief(vec![-0.0], vec![vec![-0.0]]).unwrap();
    let positive = belief(vec![0.0], vec![vec![0.0]]).unwrap();
    assert_eq!(negative.mean()[0].to_bits(), 0.0_f64.to_bits());
    assert_eq!(negative.covariance()[0][0].to_bits(), 0.0_f64.to_bits());
    assert_eq!(negative, positive);

    let negative_obs = point_sensor(0, 1, -0.0, 1.0, "x").unwrap();
    let positive_obs = point_sensor(0, 1, 0.0, 1.0, "x").unwrap();
    let negative_out = assimilate_colored(&negative, &[negative_obs], "Re", -0.0, 0.0).unwrap();
    let positive_out = assimilate_colored(&positive, &[positive_obs], "Re", 0.0, 0.0).unwrap();
    assert_eq!(negative_out, positive_out);
}

#[test]
fn finite_inputs_that_overflow_are_rejected_without_panicking() {
    let prior = scalar(0.0, f64::MAX);
    let obs = Observation::new(vec![f64::MAX], 0.0, 1.0, "x").unwrap();
    assert_eq!(
        assimilate(&prior, &obs),
        Err(AssimError::NonFiniteComputation {
            stage: "covariance-times-operator",
        })
    );

    let prior = scalar(f64::MAX, 1.0);
    let obs = sensor(0, 1, -f64::MAX, 1.0, "x");
    assert_eq!(
        misfit(&prior, &[obs]),
        Err(AssimError::NonFiniteComputation {
            stage: "misfit residual",
        })
    );
}

#[test]
fn covariance_update_avoids_an_unnecessary_intermediate_overflow() {
    let prior = scalar(0.0, 1e100);
    let obs = Observation::new(vec![1e100], 0.0, 1.0, "x").unwrap();
    let posterior = assimilate(&prior, &obs).unwrap();
    assert!(posterior.variance(0).unwrap().is_finite());
}

#[test]
fn correlated_update_preserves_exact_covariance_symmetry() {
    let prior = belief(vec![0.0, 0.0], vec![vec![2.0, 0.5], vec![0.5, 1.0]]).unwrap();
    let obs = Observation::new(vec![1.0, 2.0], 1.0, 0.25, "x").unwrap();
    let posterior = assimilate(&prior, &obs).unwrap();
    assert_eq!(
        posterior.covariance()[0][1].to_bits(),
        posterior.covariance()[1][0].to_bits()
    );
    assert!(posterior.variance(0).unwrap() <= prior.variance(0).unwrap());
    assert!(posterior.variance(1).unwrap() <= prior.variance(1).unwrap());
}

#[test]
fn joseph_update_preserves_psd_for_the_cancellation_counterexample() {
    let prior = belief(
        vec![0.0, 0.0],
        vec![
            vec![333_391_946_697.748_96, -472_122_745_149_250.6],
            vec![-472_122_745_149_250.6, 6.685_821_016_386_806e17],
        ],
    )
    .expect("counterexample prior is positive semidefinite");
    let observation = Observation::new(
        vec![-5.327_771_728_981_161, -1.892_343_097_733_501_6],
        0.0,
        1.656_506_156_123_326_4e-8,
        "adversarial-sensor",
    )
    .unwrap();

    let posterior = assimilate(&prior, &observation).expect("Joseph update must stay admissible");
    validate(&posterior).expect("every returned posterior must pass public belief validation");
    assert_eq!(
        posterior.covariance()[0][1].to_bits(),
        posterior.covariance()[1][0].to_bits(),
        "the computed upper triangle must be mirrored exactly"
    );
    let correlation = posterior.covariance()[0][1]
        / posterior.variance(0).unwrap().sqrt()
        / posterior.variance(1).unwrap().sqrt();
    assert!(
        correlation.abs() <= 1.0,
        "posterior correlation escaped [-1, 1]: {correlation}"
    );
}

#[test]
fn randomized_well_conditioned_updates_always_return_valid_beliefs() {
    fn sample(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*seed >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    }

    let mut seed = 0xA551_41A7_E5E0_0001_u64;
    for dimension in 1..=6 {
        for case in 0..128 {
            // Build a strictly positive-definite covariance as L L^T. The
            // bounded off-diagonal entries and diagonals above one keep this
            // randomized battery away from an ambiguous numerical boundary;
            // the dedicated regression above supplies the adversarial case.
            let mut factor = vec![vec![0.0; dimension]; dimension];
            for (row, factor_row) in factor.iter_mut().enumerate() {
                for (column, entry) in factor_row.iter_mut().enumerate().take(row + 1) {
                    *entry = if row == column {
                        1.0 + sample(&mut seed).abs()
                    } else {
                        sample(&mut seed) * 0.25
                    };
                }
            }
            let mut covariance = vec![vec![0.0; dimension]; dimension];
            for row in 0..dimension {
                for column in row..dimension {
                    let mut entry = 0.0;
                    for (left, right) in factor[row].iter().zip(&factor[column]) {
                        entry = left.mul_add(*right, entry);
                    }
                    covariance[row][column] = entry;
                    covariance[column][row] = entry;
                }
            }
            let mean = (0..dimension)
                .map(|_| sample(&mut seed) * 10.0)
                .collect::<Vec<_>>();
            let prior = belief(mean, covariance).expect("constructed SPD prior");
            let mut operator = (0..dimension)
                .map(|_| sample(&mut seed))
                .collect::<Vec<_>>();
            operator[0] += 1.0;
            let observation = Observation::new(
                operator,
                sample(&mut seed) * 10.0,
                0.05 + sample(&mut seed).abs(),
                "seeded-sensor",
            )
            .unwrap();

            let posterior = assimilate(&prior, &observation).unwrap_or_else(|error| {
                panic!("dimension {dimension}, case {case}: update refused: {error}")
            });
            validate(&posterior).unwrap_or_else(|error| {
                panic!("dimension {dimension}, case {case}: invalid posterior: {error}")
            });
            for row in 0..dimension {
                for column in (row + 1)..dimension {
                    assert_eq!(
                        posterior.covariance()[row][column].to_bits(),
                        posterior.covariance()[column][row].to_bits(),
                        "dimension {dimension}, case {case}: asymmetric covariance"
                    );
                }
            }
        }
    }
}

#[test]
fn assimilation_is_deterministic() {
    let prior = diagonal(vec![0.0, 0.0], &[5.0, 5.0]).unwrap();
    let obs = vec![sensor(0, 2, 2.0, 0.5, "a"), sensor(1, 2, 1.0, 0.5, "b")];
    assert_eq!(assimilate_all(&prior, &obs), assimilate_all(&prior, &obs));
    assert_eq!(
        assimilate_colored(&prior, &obs, "Re", 1.0, 2.0),
        assimilate_colored(&prior, &obs, "Re", 1.0, 2.0)
    );
}

#[test]
fn g4_pre_cancelled_entry_points_publish_nothing() {
    let gate = CancelGate::new();
    gate.request();
    let prior = scalar(0.0, 1.0);
    let observations = [sensor(0, 1, 1.0, 0.25, "pre-cancelled")];

    let colored = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        assimilate_colored_with_cx(&prior, &observations, "Re", 1.0, 2.0, cx)
    });
    assert_initial_cancelled(colored);

    let constructor = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        Belief::new(vec![0.0], vec![vec![1.0]], cx)
    });
    assert_initial_cancelled(constructor);

    let diagonal = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        Belief::diagonal(vec![0.0], &[1.0], cx)
    });
    assert_initial_cancelled(diagonal);

    let validation = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        prior.validate(cx)
    });
    assert_initial_cancelled(validation);

    let measured = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        misfit_with_cx(&prior, &observations, cx)
    });
    assert_initial_cancelled(measured);

    let single = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        assimilate_with_cx(&prior, &observations[0], cx)
    });
    assert_initial_cancelled(single);

    let aggregate = with_configured_cx(&gate, Budget::INFINITE, ExecMode::Deterministic, |cx| {
        assimilate_all_with_cx(&prior, &observations, cx)
    });
    assert_initial_cancelled(aggregate);
}

#[test]
fn g4_hostile_maximum_cancels_at_a_bounded_mid_operation_checkpoint() {
    let dimension = MAX_DENSE_STATE_DIM;
    let prior =
        diagonal(vec![0.0; dimension], &vec![1.0; dimension]).expect("maximum admitted prior");
    let observations = (0..4)
        .map(|component| {
            sensor(
                component,
                dimension,
                component as f64,
                1.0,
                "hostile-boundary",
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        (dimension as u128).pow(3) * observations.len() as u128,
        MAX_DENSE_UPDATE_CUBIC_WORK
    );

    let gate = CancelGate::new();
    let result = with_configured_cx(
        &gate,
        Budget::new().with_poll_quota(1),
        ExecMode::Deterministic,
        |cx| assimilate_colored_with_cx(&prior, &observations, "Re", 1.0, 2.0, cx),
    );
    assert!(matches!(
        result,
        Err(AssimError::Cancelled {
            phase: "observation-validation",
            completed,
            planned,
        }) if completed > 0 && completed < planned
    ));
}

#[test]
fn g4_canonical_merge_polls_inside_a_maximum_record_comparison() {
    let dimension = MAX_DENSE_STATE_DIM;
    let prior =
        diagonal(vec![0.0; dimension], &vec![1.0; dimension]).expect("maximum admitted prior");
    let observations = (0..4)
        .map(|_| {
            Observation::new(vec![1.0; dimension], 0.0, 1.0, "merge-equal-")
                .expect("maximum-width comparison fixture")
        })
        .collect::<Vec<_>>();

    let cancelled = with_configured_cx(
        &CancelGate::new(),
        Budget::INFINITE.with_poll_quota(3_899),
        ExecMode::Deterministic,
        |cx| misfit_with_cx(&prior, &observations, cx),
    );
    assert!(matches!(
        cancelled,
        Err(AssimError::Cancelled {
            phase: "canonical-compare",
            completed: 64_376,
            planned: 190_032,
        })
    ));

    let baseline = misfit(&prior, &observations).expect("healthy comparison baseline");
    let replay = misfit(&prior, &observations).expect("cancelled comparison retains no state");
    assert_eq!(baseline.to_bits(), replay.to_bits());
}

#[test]
fn g4_quota_sweep_reaches_update_psd_hash_and_commit_then_replays_cleanly() {
    let prior = diagonal(vec![0.0; 8], &[1.0; 8]).expect("bounded prior");
    let observations = [sensor(3, 8, 1.0, 0.25, "phase-sweep")];
    let baseline =
        assimilate_colored(&prior, &observations, "Re", 1.0, 2.0).expect("baseline assimilation");
    let targets = [
        "joseph-update",
        "posterior-psd",
        "candidate-hash",
        "finalize",
    ];
    let mut seen = Vec::new();
    let mut completed = false;

    for quota in 1..=4_096 {
        let result = with_configured_cx(
            &CancelGate::new(),
            Budget::new().with_poll_quota(quota),
            ExecMode::Deterministic,
            |cx| assimilate_colored_with_cx(&prior, &observations, "Re", 1.0, 2.0, cx),
        );
        match result {
            Err(AssimError::Cancelled {
                phase,
                completed: completed_work,
                planned,
            }) => {
                assert!(completed_work <= planned);
                if targets.contains(&phase) && !seen.contains(&phase) {
                    seen.push(phase);
                    let recovered = assimilate_colored(&prior, &observations, "Re", 1.0, 2.0)
                        .expect("healthy replay after refusal");
                    assert_eq!(recovered, baseline);
                }
            }
            Ok(_) => {
                completed = true;
                break;
            }
            Err(other) => panic!("unexpected phase-sweep failure: {other}"),
        }
    }

    assert!(completed, "a sufficient finite quota must complete");
    for phase in targets {
        assert!(seen.contains(&phase), "quota sweep missed {phase}");
    }
}

#[test]
fn g4_final_checkpoint_can_refuse_publication_after_local_completion() {
    let prior = scalar(0.0, 1.0);
    let observations = [sensor(0, 1, 1.0, 0.25, "final-boundary")];
    let mut final_refusal = None;

    for quota in 1..=512 {
        let gate = CancelGate::new();
        let result = with_configured_cx(
            &gate,
            Budget::new().with_poll_quota(quota),
            ExecMode::Deterministic,
            |cx| assimilate_colored_with_cx(&prior, &observations, "Re", 1.0, 2.0, cx),
        );
        match result {
            Err(
                error @ AssimError::Cancelled {
                    phase: "finalize", ..
                },
            ) => {
                final_refusal = Some(error);
                break;
            }
            Err(AssimError::Cancelled { .. }) | Ok(_) => {}
            Err(other) => panic!("unexpected final-boundary result: {other}"),
        }
    }

    assert!(matches!(
        final_refusal,
        Some(AssimError::Cancelled {
            phase: "finalize",
            completed,
            planned,
        }) if completed <= planned && planned > 0
    ));

    let recovered = assimilate_colored(&prior, &observations, "Re", 1.0, 2.0)
        .expect("healthy replay after final refusal");
    let baseline =
        assimilate_colored(&prior, &observations, "Re", 1.0, 2.0).expect("clean baseline");
    assert_eq!(recovered, baseline);
}

#[test]
fn g4_shared_poll_quota_is_consumed_once_across_nested_assimilations() {
    let prior = scalar(0.0, 1.0);
    let observations = [sensor(0, 1, 1.0, 0.25, "shared-poll-quota")];
    for (ambient, requested) in [(0, 1), (1, 2)] {
        let mut inflated_remaining = requested;
        let inflated = with_configured_cx(
            &CancelGate::new(),
            Budget::INFINITE.with_poll_quota(ambient),
            ExecMode::Deterministic,
            |cx| {
                assimilate_colored_with_shared_poll_quota(
                    &prior,
                    &observations,
                    "Re",
                    1.0,
                    2.0,
                    cx,
                    &mut inflated_remaining,
                )
            },
        );
        assert_eq!(
            inflated,
            Err(AssimError::PollQuotaExceedsAmbient { requested, ambient })
        );
        assert_eq!(inflated_remaining, requested);
    }

    let required = (1..=512)
        .find(|quota| {
            let mut remaining = *quota;
            with_cx(|cx| {
                assimilate_colored_with_shared_poll_quota(
                    &prior,
                    &observations,
                    "Re",
                    1.0,
                    2.0,
                    cx,
                    &mut remaining,
                )
                .is_ok()
            })
        })
        .expect("bounded scalar assimilation has a finite poll requirement");

    let mut remaining = required;
    let (first, second) = with_cx(|cx| {
        let first = assimilate_colored_with_shared_poll_quota(
            &prior,
            &observations,
            "Re",
            1.0,
            2.0,
            cx,
            &mut remaining,
        );
        let second = assimilate_colored_with_shared_poll_quota(
            &prior,
            &observations,
            "Re",
            1.0,
            2.0,
            cx,
            &mut remaining,
        );
        (first, second)
    });
    let shared = first.expect("first nested assimilation consumes the shared slice");
    assert_eq!(remaining, 0);
    assert!(matches!(
        second,
        Err(AssimError::Cancelled {
            phase: "initial",
            completed: 0,
            ..
        })
    ));
    let ambient =
        assimilate_colored(&prior, &observations, "Re", 1.0, 2.0).expect("ambient-budget baseline");
    assert_same_non_identity_posterior(&ambient, &shared);
    assert_ne!(estimator_identity(&ambient), estimator_identity(&shared));
}

#[test]
fn g5_candidate_identity_binds_mode_budget_and_every_stream_field() {
    let prior = scalar(0.0, 1.0);
    let observations = [sensor(0, 1, 1.0, 0.25, "identity-provenance")];
    let base = Budget::new()
        .with_poll_quota(10_000)
        .with_cost_quota(20_000)
        .with_priority(73);
    let deadline = Budget::with_deadline_at_ns(123_456)
        .with_poll_quota(10_000)
        .with_cost_quota(20_000)
        .with_priority(73);
    let variants = [
        (ExecMode::Deterministic, base, TEST_STREAM),
        (ExecMode::Fast, base, TEST_STREAM),
        (ExecMode::Deterministic, deadline, TEST_STREAM),
        (
            ExecMode::Deterministic,
            base.with_poll_quota(10_001),
            TEST_STREAM,
        ),
        (
            ExecMode::Deterministic,
            base.with_cost_quota(20_001),
            TEST_STREAM,
        ),
        (ExecMode::Deterministic, base.with_priority(74), TEST_STREAM),
        (
            ExecMode::Deterministic,
            base,
            StreamKey {
                seed: TEST_STREAM.seed + 1,
                ..TEST_STREAM
            },
        ),
        (
            ExecMode::Deterministic,
            base,
            StreamKey {
                kernel_id: TEST_STREAM.kernel_id + 1,
                ..TEST_STREAM
            },
        ),
        (
            ExecMode::Deterministic,
            base,
            StreamKey {
                tile: TEST_STREAM.tile + 1,
                ..TEST_STREAM
            },
        ),
        (
            ExecMode::Deterministic,
            base,
            StreamKey {
                iteration: TEST_STREAM.iteration + 1,
                ..TEST_STREAM
            },
        ),
    ];
    let variant_count = variants.len();

    let mut identities = Vec::new();
    let mut reference = None;
    for (mode, budget, stream) in variants {
        let gate = CancelGate::new();
        let output = with_stream_cx(&gate, budget, mode, stream, |cx| {
            assimilate_colored_with_cx(&prior, &observations, "Re", 1.0, 2.0, cx)
        })
        .expect("provenance variant must complete");
        if let Some(reference) = &reference {
            assert_same_non_identity_posterior(reference, &output);
        } else {
            reference = Some(output.clone());
        }
        identities.push(estimator_identity(&output).to_owned());
    }
    assert!(
        identities
            .iter()
            .all(|identity| identity.starts_with("assimilation-candidate:v3:"))
    );
    identities.sort();
    identities.dedup();
    assert_eq!(identities.len(), variant_count);
}
