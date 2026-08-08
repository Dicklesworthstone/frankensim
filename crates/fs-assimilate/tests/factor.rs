//! G0/G1/G3/G4/G5 battery for the stable factor-form scalar substrate
//! (bead sj31i.37): Thornton-Bierman UD updates, executable contraction
//! receipts, and the independent dense-Joseph checker lane.

use fs_assimilate::factor::{
    CheckVerdict, ContractionState, FactorBelief, MisfitVerdict, assimilate_belief_scalar_checked,
    assimilate_scalar, scalar_factor_work_estimate, verify_factor_assimilation,
};
use fs_assimilate::{AssimError, Belief, Observation, point_sensor};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x00FA_C708,
    kernel_id: 0xFA07,
    tile: 1,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = fs_exec::VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, TEST_STREAM, budget, ExecMode::Deterministic)
            .with_time_source(&clock);
        f(&cx)
    });
    let stats = pool.stats();
    assert!(
        stats.quiescent(),
        "Cx arena must be quiescent after scope: {}",
        stats.to_json()
    );
    result
}

fn splitmix(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
}

/// Build a dense PSD covariance from a deterministic factor pair and admit
/// it as a validated `Belief`.
fn dense_belief(mean: &[f64], diag: &[f64], upper: &[(usize, usize, f64)], cx: &Cx<'_>) -> Belief {
    let n = mean.len();
    let mut cov = vec![vec![0.0; n]; n];
    for (i, row) in cov.iter_mut().enumerate() {
        row[i] = diag[i];
    }
    for &(i, j, value) in upper {
        cov[i][j] = value;
        cov[j][i] = value;
    }
    // Rebuild exactly from the implied factor so the matrix is PSD by
    // construction: cov = U D U^T with unit upper U from `upper`.
    let mut u = vec![vec![0.0; n]; n];
    for (i, row) in u.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for &(i, j, value) in upper {
        u[i][j] = value;
    }
    let mut exact = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let mut entry = 0.0;
            for k in j..n {
                entry += u[i][k] * diag[k] * u[j][k];
            }
            exact[i][j] = entry;
            exact[j][i] = entry;
        }
    }
    Belief::new(mean.to_vec(), exact, cx).expect("dense belief fixture")
}

fn max_cov_diff(left: &[Vec<f64>], right: &[Vec<f64>]) -> f64 {
    let mut max = 0.0_f64;
    for (lrow, rrow) in left.iter().zip(right) {
        for (l, r) in lrow.iter().zip(rrow) {
            let diff = (l - r).abs();
            if diff > max {
                max = diff;
            }
        }
    }
    max
}

fn max_vec_diff(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(l, r)| (l - r).abs())
        .fold(0.0_f64, f64::max)
}

#[test]
fn g0_single_component_golden_update() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![1.0], &[2.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0], 3.0, 0.5, "golden-scalar").expect("obs");
        let result = assimilate_scalar(&prior, &obs, cx).expect("update");
        // a = 2.5, K = 0.8, mean' = 2.6, var' = 0.4 exactly in binary64.
        assert!((result.belief().mean()[0] - 2.6).abs() <= 1e-15);
        let post_var = result.belief().variance(0).expect("variance");
        assert!((post_var - 0.4).abs() <= 1e-15);
        assert_eq!(result.receipt().state(), ContractionState::Certified);
        assert_eq!(
            result.receipt().misfit_verdict(),
            MisfitVerdict::NonIncreasing
        );
        assert!((result.receipt().innovation_variance() - 2.5).abs() <= 1e-15);
        assert!(result.receipt().max_pivot_ratio() <= 1.0);
        assert!(
            result
                .receipt()
                .identity()
                .starts_with("scalar-contraction:v1:")
        );
    });
}

#[test]
fn g0_two_component_bierman_matches_direct_dense() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        // P = [[1.5, 1.0], [1.0, 2.0]] from U=[[1,0.5],[0,1]], D=diag(1,2).
        let prior = dense_belief(&[0.0, 0.0], &[1.0, 2.0], &[(0, 1, 0.5)], cx);
        let obs = Observation::new(vec![1.0, 1.0], 1.0, 1.0, "bierman-n2").expect("obs");
        let factor_prior = FactorBelief::from_belief(&prior, cx).expect("factor");
        let result = assimilate_scalar(&factor_prior, &obs, cx).expect("update");
        // Direct dense: P' = P - Ph Ph^T / a with Ph = [2.5, 3], a = 6.5.
        let expected = [
            [1.5 - 6.25 / 6.5, 1.0 - 7.5 / 6.5],
            [1.0 - 7.5 / 6.5, 2.0 - 9.0 / 6.5],
        ];
        let dense = result.belief().to_dense_covariance();
        assert!(
            max_cov_diff(&dense, &expected.map(|row| row.to_vec()).to_vec()) <= 1e-12,
            "factor posterior must match the direct dense rank-one downdate: {dense:?}"
        );
        assert_eq!(result.receipt().state(), ContractionState::Certified);
        let dense_check = verify_factor_assimilation(&prior, &obs, &result, cx).expect("check");
        assert_eq!(dense_check.verdict(), CheckVerdict::Verified);
    });
}

#[test]
fn g0_oracle_differential_across_deterministic_priors() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let dims = [1_usize, 2, 3, 8, 17];
        let mut rng = 0xB1E5_0A11_0001_u64;
        for &n in &dims {
            let mean: Vec<f64> = (0..n).map(|_| splitmix(&mut rng) * 2.0 - 1.0).collect();
            let diag: Vec<f64> = (0..n).map(|_| 0.25 + splitmix(&mut rng) * 3.0).collect();
            let mut upper = Vec::new();
            for i in 0..n {
                for j in (i + 1)..n {
                    upper.push((i, j, (splitmix(&mut rng) - 0.5) * 0.4));
                }
            }
            let prior = dense_belief(&mean, &diag, &upper, cx);
            let operator: Vec<f64> = (0..n).map(|_| 0.5 + splitmix(&mut rng)).collect();
            let obs = Observation::new(
                operator,
                splitmix(&mut rng) * 4.0 - 2.0,
                0.1 + splitmix(&mut rng),
                format!("differential-{n}"),
            )
            .expect("obs");
            let factor_prior = FactorBelief::from_belief(&prior, cx).expect("factor");
            let result = assimilate_scalar(&factor_prior, &obs, cx).expect("update");

            let oracle = fs_assimilate::assimilate(&prior, &obs, cx).expect("oracle");
            let factor_dense = result.belief().to_dense_covariance();
            let diff = max_cov_diff(&factor_dense, oracle.covariance());
            let scale = oracle
                .covariance()
                .iter()
                .flat_map(|row| row.iter())
                .fold(0.0_f64, |acc, v| acc.max(v.abs()))
                .max(1.0);
            assert!(
                diff <= 1e-9 * scale,
                "n={n}: factor vs Joseph oracle max diff {diff} (scale {scale})"
            );
            assert!(
                max_vec_diff(result.belief().mean(), oracle.mean()) <= 1e-12 * scale,
                "n={n}: posterior means diverge"
            );
            assert_eq!(
                result.receipt().state(),
                ContractionState::Certified,
                "n={n}: well-conditioned update must certify"
            );
            assert_eq!(
                result.receipt().misfit_verdict(),
                MisfitVerdict::NonIncreasing,
                "n={n}: misfit must not increase"
            );
            let check = verify_factor_assimilation(&prior, &obs, &result, cx).expect("check");
            assert_eq!(
                check.verdict(),
                CheckVerdict::Verified,
                "n={n}: independent checker must verify"
            );
        }
    });
}

#[test]
fn g0_diagonal_variances_never_expand_pointwise() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior =
            FactorBelief::diagonal(vec![0.0, 0.0, 0.0], &[3.0, 1.0, 2.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 1.0, 1.0], 0.5, 0.5, "monotone-diag").expect("obs");
        let prior_variances: Vec<f64> = (0..3).map(|i| prior.variance(i).expect("v")).collect();
        let result = assimilate_scalar(&prior, &obs, cx).expect("update");
        for (i, prior_var) in prior_variances.iter().enumerate() {
            let post_var = result.belief().variance(i).expect("v");
            assert!(
                post_var <= *prior_var,
                "component {i} variance expanded: {prior_var} -> {post_var}"
            );
        }
        assert_eq!(result.receipt().state(), ContractionState::Certified);
    });
}

#[test]
fn g0_rank_singular_prior_stays_exact_zero_in_singular_direction() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        // Zero variance in component 1: that direction is an exact constant
        // and must survive the update untouched.
        let prior = FactorBelief::diagonal(vec![0.0, 5.0], &[2.0, 0.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 1.0], 6.0, 1.0, "singular-prior").expect("obs");
        let result = assimilate_scalar(&prior, &obs, cx).expect("update");
        assert_eq!(result.belief().diag(1), Some(0.0));
        assert_eq!(result.belief().mean()[1], 5.0);
        assert!((result.belief().variance(1).expect("v")).abs() == 0.0);
        // The verdict may certify or be unresolved on a singular prior, but
        // it must never affirmatively refute a contraction that holds.
        assert_ne!(result.receipt().state(), ContractionState::Refuted);
    });
}

#[test]
fn g0_illconditioned_prior_never_falsely_refutes() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior =
            FactorBelief::diagonal(vec![0.0, 0.0], &[1.0e10, 1.0e-6], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 1.0], 0.25, 1.0, "ill-conditioned").expect("obs");
        let result = assimilate_scalar(&prior, &obs, cx).expect("update");
        assert_ne!(
            result.receipt().state(),
            ContractionState::Refuted,
            "a well-posed ill-conditioned update must not be refuted"
        );
        let dense = result.belief().to_dense_covariance();
        assert!(dense[0][0].is_finite() && dense[1][1].is_finite());
    });
}

#[test]
fn g0_factor_round_trip_reconstructs_validated_covariance() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = dense_belief(
            &[1.0, -1.0, 0.5],
            &[2.0, 1.0, 3.0],
            &[(0, 1, 0.25), (0, 2, -0.5), (1, 2, 0.75)],
            cx,
        );
        let factor = FactorBelief::from_belief(&prior, cx).expect("factor");
        let round_trip = factor.to_dense_covariance();
        let diff = max_cov_diff(&round_trip, prior.covariance());
        assert!(
            diff <= 1e-12,
            "factor round trip must reconstruct the admitted covariance: {diff}"
        );
    });
}

#[test]
fn g0_receipt_identity_is_deterministic_and_input_sensitive() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![0.0, 1.0], &[2.0, 1.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 0.5], 1.5, 0.5, "identity-check").expect("obs");
        let first = assimilate_scalar(&prior, &obs, cx).expect("first");
        let second = assimilate_scalar(&prior, &obs, cx).expect("second");
        assert_eq!(first.receipt().identity(), second.receipt().identity());
        assert_eq!(first.belief().mean(), second.belief().mean());
        let changed = Observation::new(vec![1.0, 0.5], 1.5, 0.5, "identity-other").expect("obs2");
        let third = assimilate_scalar(&prior, &changed, cx).expect("third");
        assert_ne!(first.receipt().identity(), third.receipt().identity());
    });
}

#[test]
fn g1_repeated_low_noise_updates_track_analytic_variance_recursion() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let mut belief = FactorBelief::diagonal(vec![0.0], &[4.0], cx).expect("prior");
        let mut analytic = 4.0_f64;
        for step in 0..8 {
            let obs = Observation::new(vec![1.0], 0.25, 1.0, "low-noise").expect("obs");
            let result = assimilate_scalar(&belief, &obs, cx).expect("update");
            analytic = analytic * 1.0 / (analytic + 1.0);
            let computed = result.belief().variance(0).expect("variance");
            assert!(
                (computed - analytic).abs() <= 1e-12,
                "step {step}: variance {computed} vs analytic {analytic}"
            );
            assert_eq!(result.receipt().state(), ContractionState::Certified);
            assert_eq!(
                result.receipt().misfit_verdict(),
                MisfitVerdict::NonIncreasing
            );
            belief = result.belief().clone();
        }
        assert!(analytic < 4.0);
    });
}

#[test]
fn g3_state_permutation_is_metamorphic() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let mean = [0.5, -1.0, 2.0];
        let diag = [2.0, 1.0, 3.0];
        let upper = [(0_usize, 1_usize, 0.5_f64), (0, 2, -0.25), (1, 2, 0.5)];
        let prior = dense_belief(&mean, &diag, &upper, cx);
        let obs = Observation::new(vec![1.0, 0.5, -0.5], 1.25, 0.75, "permutation").expect("obs");
        let base_factor = FactorBelief::from_belief(&prior, cx).expect("factor");
        let base = assimilate_scalar(&base_factor, &obs, cx).expect("base update");

        // Permutation (0->1, 1->2, 2->0): permuted prior/operator must give
        // the permuted posterior. The permuted covariance is built DENSELY
        // (permuting factor entries does not yield an upper factor), and
        // from_belief re-factors it.
        let perm = [1_usize, 2, 0];
        let permute = |values: &[f64]| -> Vec<f64> { perm.iter().map(|&i| values[i]).collect() };
        let base_cov = prior.covariance();
        let n = 3;
        let mut perm_cov = vec![vec![0.0; n]; n];
        for a in 0..n {
            for b in 0..n {
                perm_cov[a][b] = base_cov[perm[a]][perm[b]];
            }
        }
        let perm_mean = permute(&mean);
        let perm_prior = Belief::new(perm_mean, perm_cov, cx).expect("permuted prior");
        let perm_operator = permute(obs.operator());
        let perm_obs = Observation::new(perm_operator, obs.value(), obs.noise_var(), "permutation")
            .expect("perm obs");
        let perm_factor = FactorBelief::from_belief(&perm_prior, cx).expect("perm factor");
        let permuted = assimilate_scalar(&perm_factor, &perm_obs, cx).expect("perm update");

        let base_dense = base.belief().to_dense_covariance();
        let perm_dense = permuted.belief().to_dense_covariance();
        let mut permuted_back = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                permuted_back[perm[i]][perm[j]] = perm_dense[i][j];
            }
        }
        let diff = max_cov_diff(&permuted_back, &base_dense);
        assert!(diff <= 1e-10, "permuted covariance diff {diff}");
        let base_mean = base.belief().mean();
        let base_mean_forward = permute(base_mean);
        assert!(max_vec_diff(&base_mean_forward, permuted.belief().mean()) <= 1e-12);
    });
}

#[test]
fn g4_exhausted_poll_quota_cancels_without_partial_output() {
    let gate = CancelGate::new();
    let (prior, obs) = with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![0.0, 0.0], &[2.0, 1.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 1.0], 1.0, 0.5, "cancel-probe").expect("obs");
        (prior, obs)
    });
    let refusal = with_cx(&gate, Budget::INFINITE.with_poll_quota(0), |cx| {
        assimilate_scalar(&prior, &obs, cx)
            .expect_err("zero poll quota must refuse before publication")
    });
    assert!(matches!(
        refusal,
        AssimError::Cancelled { .. } | AssimError::BudgetRefused(_)
    ));
}

#[test]
fn g4_work_estimate_is_checked_and_monotone() {
    let small = scalar_factor_work_estimate(1).expect("n=1");
    let large = scalar_factor_work_estimate(256).expect("n=256");
    assert!(large.construction > small.construction);
    assert!(large.update_with_receipt > small.update_with_receipt);
    assert!(large.independent_check > small.independent_check);
    assert!(scalar_factor_work_estimate(usize::MAX).is_err());
}

#[test]
fn e2e_point_sensor_to_factor_update_to_independent_check() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = dense_belief(
            &[0.0, 1.0, -1.0],
            &[2.0, 1.5, 1.0],
            &[(0, 1, 0.5), (1, 2, -0.25)],
            cx,
        );
        let obs = point_sensor(1, 3, 2.5, 0.5, "e2e-sensor").expect("point sensor");
        let checked = assimilate_belief_scalar_checked(&prior, &obs, cx).expect("checked path");
        assert_eq!(checked.receipt().state(), ContractionState::Certified);
        assert_eq!(checked.check().verdict(), CheckVerdict::Verified);
        assert!(checked.check().identity_consistent());
        let prior_var = 1.5 + 0.25; // P11 = d1 + u12^2 d2 ... via fixture factor
        let post_var = checked.posterior().variance(1).expect("variance");
        assert!(
            post_var < prior_var,
            "observed component variance must contract: {prior_var} -> {post_var}"
        );
        assert!(
            checked
                .receipt()
                .identity()
                .starts_with("scalar-contraction:v1:"),
            "receipt identity prefix"
        );
    });
}

#[test]
fn g0_deep_subnormal_noise_never_panics_and_stays_decisive_or_unresolved() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![0.0, 0.0], &[1.0, 1.0], cx).expect("prior");
        let obs = Observation::new(vec![1.0, 1.0], 1.0, 1.0e-300, "deep-subnormal").expect("obs");
        let result = assimilate_scalar(&prior, &obs, cx).expect("update");
        // A subnormal-noise update is numerically legal; the receipt must be
        // decisive evidence, never a panic and never a false refutation.
        assert_ne!(result.receipt().state(), ContractionState::Refuted);
        assert!(result.receipt().innovation_variance() > 0.0);
    });
}

#[test]
fn g2_high_dynamic_range_metrology_battery() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        // A metrology-scale state: ten decades of variance across channels,
        // observed at each channel's own scale. Every update must remain
        // stable and every receipt decisive or honestly unresolved.
        let scales: [f64; 4] = [1.0e-6, 1.0e-2, 1.0e2, 1.0e6];
        let variances: Vec<f64> = scales.iter().map(|s| s * s).collect();
        let mut belief = FactorBelief::diagonal(vec![0.0; 4], &variances, cx).expect("prior");
        for (channel, scale) in scales.iter().enumerate() {
            let mut operator = vec![0.0; 4];
            operator[channel] = 1.0;
            let reading = 0.25 * scale;
            let noise = (0.05 * scale) * (0.05 * scale);
            let obs = Observation::new(operator, reading, noise, "metrology").expect("obs");
            let result = assimilate_scalar(&belief, &obs, cx).expect("update");
            assert!(
                result.receipt().state() != ContractionState::Refuted,
                "channel {channel} (scale {scale:e}) falsely refuted"
            );
            let post_var = result.belief().variance(channel).expect("variance");
            let prior_var = scale * scale;
            assert!(
                post_var < prior_var,
                "channel {channel}: variance must contract ({prior_var:e} -> {post_var:e})"
            );
            assert_eq!(
                result.receipt().misfit_verdict(),
                MisfitVerdict::NonIncreasing
            );
            belief = result.belief().clone();
        }
        // The unobserved-scale channels still contract through the shared
        // update only when correlated; here channels are independent, so
        // each channel contracts exactly at its own observation.
        let dense = belief.to_dense_covariance();
        for (i, scale) in scales.iter().enumerate() {
            assert!(dense[i][i].is_finite() && dense[i][i] > 0.0);
            assert!(dense[i][i] < scale * scale);
        }
    });
}

#[test]
fn e2e_log_emission_is_bounded_schema_valid_and_complete() {
    let gate = CancelGate::new();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let prior = dense_belief(
            &[0.0, 1.0, -1.0],
            &[2.0, 1.5, 1.0],
            &[(0, 1, 0.5), (1, 2, -0.25)],
            cx,
        );
        let obs = point_sensor(1, 3, 2.5, 0.5, "log-sensor").expect("point sensor");
        let checked = assimilate_belief_scalar_checked(&prior, &obs, cx).expect("checked path");
        let planned = scalar_factor_work_estimate(3).expect("estimate");
        let mut emitter = fs_obs::Emitter::new(
            fs_assimilate::factor::SCALAR_FACTOR_LOG_SUITE,
            "log-emission-test",
        );
        let event = fs_assimilate::factor::emit_checked_assimilation_log(
            &checked,
            3,
            planned.update_with_receipt,
            &mut emitter,
        )
        .expect("log emission validates");
        let line = event.to_jsonl();
        fs_obs::validate_line(&line).expect("wire schema");
        assert!(line.contains("bierman-ud/v1"));
        assert!(line.contains("\\\"contraction\\\":\\\"certified\\\""));
        assert!(line.contains("\\\"checker\\\":\\\"verified\\\""));
        assert!(line.contains("scalar-contraction:v1:"));
        assert!(line.contains("\\\"dim\\\":3"));
    });
}
