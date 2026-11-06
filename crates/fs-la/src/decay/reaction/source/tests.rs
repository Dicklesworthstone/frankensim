use super::*;

fn edge(parent: usize, daughter: usize, yield_per_removal: f64) -> ReactionEdge {
    ReactionEdge {
        parent,
        daughter,
        yield_per_removal,
    }
}

fn plan(rates: &[f64], edges: &[ReactionEdge], time: f64) -> ReactionSourcePlan {
    ReactionSourcePlan::prepare(rates, edges, time, DecayLimits::default(), || false).unwrap()
}

fn near(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!(
            (actual - expected).abs() <= 5e-12 * expected.abs(),
            "actual {actual:.17e}, expected {expected:.17e}"
        );
    }
}

#[test]
fn uniform_feed_into_a_finite_sink_matches_the_inhomogeneous_scalar_solution() {
    let rate: f64 = 0.7;
    let time = 1.3;
    let phi = -(-rate * time).exp_m1() / (rate * time);
    let prepared = plan(&[rate], &[], time);
    let answer = prepared.apply(&[5.0], &[17.0]).unwrap();
    near(
        answer.populations[0],
        5.0 * (-rate * time).exp() + 17.0 * phi,
    );
    near(
        answer.removals[0],
        5.0 * -(-rate * time).exp_m1() + 17.0 * (1.0 - phi),
    );
    near(prepared.source_population_matrix()[0], phi);
    near(prepared.source_removal_matrix()[0], 1.0 - phi);
    assert!(answer.balance_residuals[0].abs() < 1e-13);
    assert!(answer.population_balance_residual.abs() < 1e-13);
}

#[test]
fn equal_rate_multiplying_daughters_react_during_their_production_interval() {
    let x: f64 = 1.4;
    let parent = -(-x).exp_m1() / x;
    let daughter = 2.0 * (1.0 - (1.0 + x) * (-x).exp()) / x;
    let parent_removals = 1.0 - parent;
    let daughter_removals = 2.0 * parent_removals - daughter;
    let prepared = plan(&[2.0, 2.0, 0.0], &[edge(0, 1, 2.0), edge(1, 2, 0.5)], 0.7);
    let answer = prepared.apply(&[0.0; 3], &[100.0, 0.0, 0.0]).unwrap();
    near(answer.populations[0], 100.0 * parent);
    near(answer.populations[1], 100.0 * daughter);
    near(answer.populations[2], 50.0 * daughter_removals);
    near(answer.removals[0], 100.0 * parent_removals);
    near(answer.removals[1], 100.0 * daughter_removals);
    assert_eq!(answer.removals[2], 0.0);
    // Multiplicity and the finite sink are retained, rather than normalized
    // into a one-daughter conservation model.
    near(
        answer.populations.iter().sum(),
        100.0 + answer.removals[0] - 0.5 * answer.removals[1],
    );
}

#[test]
fn source_self_return_and_multiplication_use_the_net_generator_but_actual_removals() {
    let rate = 0.3;
    let time = 4.0;
    let neutral = plan(&[rate], &[edge(0, 0, 1.0)], time);
    let answer = neutral.apply(&[2.0], &[7.0]).unwrap();
    assert_eq!(answer.populations, [9.0]);
    near(answer.removals[0], rate * time * (2.0 + 7.0 / 2.0));
    assert_eq!(neutral.source_population_matrix(), [1.0]);
    near(neutral.source_removal_matrix()[0], rate * time / 2.0);

    let multiplying = plan(&[rate], &[edge(0, 0, 2.0)], time);
    let growth = (rate * time).exp();
    let phi = (rate * time).exp_m1() / (rate * time);
    let answer = multiplying.apply(&[2.0], &[7.0]).unwrap();
    near(answer.populations[0], 2.0 * growth + 7.0 * phi);
    near(answer.removals[0], 2.0 * (growth - 1.0) + 7.0 * (phi - 1.0));
    assert!(multiplying.source_population_matrix()[0] > 1.0);
    near(multiplying.source_removal_matrix()[0], phi - 1.0);
}

#[test]
fn fed_reversible_cycle_matches_an_independent_affine_forcing_solution() {
    let (a, b, time): (f64, f64, f64) = (0.7, 1.3, 4.3);
    let (u, v) = (3.0, 1.5);
    let total_initial = 14.0;
    let k = a + b;
    let slope = b * (u + v) / k;
    let constant = (b * total_initial + u - slope) / k;
    let departure = 10.0 - constant;
    let first = slope * time + constant + departure * (-k * time).exp();
    let first_integral =
        slope * time * time / 2.0 + constant * time + departure * -(-k * time).exp_m1() / k;
    let total = total_initial + (u + v) * time;
    let total_integral = total_initial * time + (u + v) * time * time / 2.0;
    let prepared = plan(&[a, b], &[edge(0, 1, 1.0), edge(1, 0, 1.0)], time);
    let answer = prepared.apply(&[10.0, 4.0], &[u * time, v * time]).unwrap();
    near(answer.populations[0], first);
    near(answer.populations[1], total - first);
    near(answer.removals[0], a * first_integral);
    near(answer.removals[1], b * (total_integral - first_integral));
    assert!(
        answer.removals[0] > 10.0 + u * time,
        "cycle revisits count as events"
    );
    assert!(answer.population_balance_residual.abs() < 1e-11);
}

#[test]
fn multiplying_cyclic_feed_composes_and_preserves_original_order_and_reuse() {
    let rates = [0.7, 0.2, 0.04, 0.0];
    let edges = [
        edge(0, 1, 1.7),
        edge(1, 0, 0.35),
        edge(1, 2, 1.1),
        edge(2, 3, 1.0),
    ];
    let initial = [11.0, 3.0, 0.0, 0.5];
    let prepared = plan(&rates, &edges, 2.3);
    let total = [4.6, 6.9, 0.0, 2.3];
    let whole = prepared.apply(&initial, &total).unwrap();
    let first = plan(&rates, &edges, 0.9)
        .apply(&initial, &[1.8, 2.7, 0.0, 0.9])
        .unwrap();
    let second = plan(&rates, &edges, 1.4)
        .apply(&first.populations, &[2.8, 4.2, 0.0, 1.4])
        .unwrap();
    let permuted = plan(
        &[0.2, 0.0, 0.7, 0.04],
        &[
            edge(2, 0, 1.7),
            edge(0, 2, 0.35),
            edge(0, 3, 1.1),
            edge(3, 1, 1.0),
        ],
        2.3,
    )
    .apply(&[3.0, 0.5, 11.0, 0.0], &[6.9, 2.3, 4.6, 0.0])
    .unwrap();
    let triple = prepared
        .apply(&initial.map(|n| 3.0 * n), &total.map(|n| 3.0 * n))
        .unwrap();
    for (i, p) in [2, 0, 3, 1].into_iter().enumerate() {
        near(second.populations[i], whole.populations[i]);
        near(first.removals[i] + second.removals[i], whole.removals[i]);
        near(permuted.populations[p], whole.populations[i]);
        near(permuted.removals[p], whole.removals[i]);
        near(triple.populations[i], 3.0 * whole.populations[i]);
    }
    assert_eq!(
        prepared.apply(&initial, &[0.0; 4]).unwrap(),
        prepared.homogeneous_plan().apply(&initial).unwrap(),
    );
}

#[test]
fn stable_and_zero_intervals_keep_exact_limits_without_consuming_extra_nodes() {
    let stable = plan(&[0.0, 0.0], &[edge(0, 1, 2.0)], f64::MAX);
    let answer = stable.apply(&[2.0, 3.0], &[7.0, 11.0]).unwrap();
    assert_eq!(answer.populations, [9.0, 14.0]);
    assert_eq!(answer.removals, [0.0, 0.0]);
    assert_eq!(answer.balance_residuals, [0.0, 0.0]);
    assert_eq!(answer.population_balance_residual, 0.0);
    let zero = plan(&[1.0], &[], 0.0);
    assert_eq!(zero.apply(&[2.0], &[0.0]).unwrap().populations, [2.0]);
    assert!(matches!(
        zero.apply(&[2.0], &[1.0]),
        Err(DecayError::Invalid(_))
    ));
    let maximum = plan(&[0.0; 64], &[], 1.0);
    assert_eq!(maximum.homogeneous_plan().node_count(), 64);
    assert_eq!(
        maximum.apply(&[0.0; 64], &[1.0; 64]).unwrap().populations,
        [1.0; 64]
    );
}

#[test]
fn tiny_fresh_removals_survive_stiffness_and_small_feeds_keep_their_balance_debit() {
    let prepared = plan(
        &[1e-20, 0.0, 1e20, 0.0],
        &[edge(0, 1, 2.0), edge(2, 3, 1.0)],
        1.0,
    );
    let answer = prepared.apply(&[0.0; 4], &[1e30, 0.0, 0.0, 0.0]).unwrap();
    near(answer.removals[0], 5e9);
    near(answer.populations[1], 1e10);
    assert!(prepared.homogeneous_plan().squarings() > 60);
    let unresolved = plan(&[0.0], &[], 1.0).apply(&[1e30], &[1.0]).unwrap();
    assert_eq!(unresolved.populations, [1e30]);
    assert_eq!(unresolved.balance_residuals, [-1.0]);
    assert_eq!(unresolved.population_balance_residual, -1.0);
}

#[test]
fn continuing_feed_cannot_inherit_an_exhausted_initial_inventory_exemption() {
    let prepared = ReactionSourcePlan::prepare(
        &[1e200],
        &[],
        1e200,
        DecayLimits {
            maximum_series_terms: 256,
            ..DecayLimits::default()
        },
        || false,
    )
    .unwrap();
    assert!(prepared.homogeneous_plan().apply(&[1e200]).is_ok());
    assert!(matches!(
        prepared.apply(&[0.0], &[1e200]),
        Err(DecayError::Unrepresentable(_))
    ));
    let small = prepared.apply(&[0.0], &[1.0]).unwrap();
    assert_eq!(small.populations, [0.0]);
    near(small.removals[0], 1.0);
    let finite = plan(&[1e20], &[], 1.0);
    near(
        finite.apply(&[0.0], &[1e200]).unwrap().populations[0],
        1e180,
    );
}

#[test]
fn source_work_input_and_checkpoint_refusals_leave_inputs_and_cached_maps_intact() {
    let rates = [0.7, 1.3];
    let edges = [edge(0, 1, 1.0), edge(1, 0, 1.0)];
    let mut polls = 0;
    let prepared = ReactionSourcePlan::prepare(&rates, &edges, 2.7, DecayLimits::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert_eq!(
            ReactionSourcePlan::prepare(&rates, &edges, 2.7, DecayLimits::default(), || {
                seen += 1;
                seen == stop
            })
            .unwrap_err(),
            DecayError::Cancelled,
        );
    }
    let work = prepared.homogeneous_plan().preparation_work();
    let limits = DecayLimits {
        maximum_work: work,
        ..DecayLimits::default()
    };
    assert!(ReactionSourcePlan::prepare(&rates, &edges, 2.7, limits, || false).is_ok());
    assert!(matches!(
        ReactionSourcePlan::prepare(
            &rates,
            &edges,
            2.7,
            DecayLimits {
                maximum_work: work - 1,
                ..limits
            },
            || false
        ),
        Err(DecayError::Budget { .. }),
    ));
    let initial = [7.0, 2.0];
    let injected = [3.0, 4.0];
    let mut calls = 0;
    let expected = prepared
        .apply_with_cancel(&initial, &injected, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert_eq!(
            prepared.apply_with_cancel(&initial, &injected, || {
                seen += 1;
                seen == stop
            }),
            Err(DecayError::Cancelled)
        );
    }
    for invalid in [
        vec![1.0],
        vec![-1.0, 0.0],
        vec![f64::NAN, 0.0],
        vec![f64::INFINITY, 0.0],
    ] {
        assert!(prepared.apply(&initial, &invalid).is_err());
        assert!(prepared.apply(&invalid, &injected).is_err());
    }
    assert_eq!(initial, [7.0, 2.0]);
    assert_eq!(injected, [3.0, 4.0]);
    assert_eq!(prepared.apply(&initial, &injected).unwrap(), expected);
}
