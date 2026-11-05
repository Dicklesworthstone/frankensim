use super::*;

fn edge(parent: usize, daughter: usize, yield_per_removal: f64) -> ReactionEdge {
    ReactionEdge {
        parent,
        daughter,
        yield_per_removal,
    }
}

fn plan(rates: &[f64], edges: &[ReactionEdge], time: f64) -> ReactionPlan {
    ReactionPlan::prepare(rates, edges, time, DecayLimits::default(), || false).unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    let allowance = relative * expected.abs().max(1.0e-300);
    assert!(
        (actual - expected).abs() <= allowance,
        "actual {actual:.17e}, expected {expected:.17e}, allowance {allowance:.3e}"
    );
}

#[test]
fn competing_products_and_finite_sink_preserve_declared_multiplicity() {
    let initial = 200.0;
    let exposure: f64 = 3.0 * 0.4;
    let removed = initial * -(-exposure).exp_m1();
    let result = plan(
        &[3.0, 0.0, 0.0],
        &[edge(0, 1, 0.05), edge(0, 1, 0.15), edge(0, 2, 1.4)],
        0.4,
    )
    .apply(&[initial, 0.0, 0.0])
    .unwrap();
    near(result.removals[0], removed, 2.0e-14);
    near(result.populations[0], initial * (-exposure).exp(), 2.0e-14);
    near(result.populations[1], 0.2 * removed, 2.0e-14);
    near(result.populations[2], 1.4 * removed, 2.0e-14);
    assert!(result.populations.iter().sum::<f64>() > initial);
    assert!(result.population_balance_residual.abs() < 1.0e-11);
    let sink = plan(&[3.0], &[], 0.4).apply(&[initial]).unwrap();
    near(sink.removals[0], removed, 2.0e-14);
    near(sink.populations[0], initial - removed, 2.0e-14);
}

#[test]
fn equal_rate_products_can_decay_within_their_birth_interval() {
    let n = 50.0;
    let x: f64 = 0.7 * 2.3;
    let result = plan(&[0.7, 0.7, 0.0], &[edge(0, 1, 2.0), edge(1, 2, 1.0)], 2.3)
        .apply(&[n, 0.0, 0.0])
        .unwrap();
    near(result.populations[0], n * (-x).exp(), 2.0e-13);
    near(result.populations[1], 2.0 * n * x * (-x).exp(), 2.0e-13);
    let granddaughter = 2.0 * n * (1.0 - (1.0 + x) * (-x).exp());
    near(result.populations[2], granddaughter, 2.0e-13);
    near(result.removals[0], n * -(-x).exp_m1(), 2.0e-13);
    near(result.removals[1], granddaughter, 2.0e-13);
    assert_eq!(result.removals[2], 0.0);
}

#[test]
fn reversible_cycle_matches_equilibrium_and_integrated_event_oracle() {
    let a = 0.7;
    let b = 1.3;
    let t = 4.3;
    let total = 14.0;
    let equilibrium = total * b / (a + b);
    let departure = 10.0 - equilibrium;
    let attenuation = f64::exp(-(a + b) * t);
    let excess_integral = departure * (1.0 - attenuation) / (a + b);
    let result = plan(&[a, b], &[edge(0, 1, 1.0), edge(1, 0, 1.0)], t)
        .apply(&[10.0, 4.0])
        .unwrap();
    near(
        result.populations[0],
        equilibrium + departure * attenuation,
        2.0e-13,
    );
    near(
        result.populations[1],
        total - result.populations[0],
        2.0e-13,
    );
    near(
        result.removals[0],
        a * (equilibrium * t + excess_integral),
        2.0e-13,
    );
    near(
        result.removals[1],
        b * ((total - equilibrium) * t - excess_integral),
        2.0e-13,
    );
    assert!(result.removals[0] > 10.0, "repeat visits are counted");
}

#[test]
fn stiff_balanced_cycle_never_hides_population_drift_behind_large_counts() {
    for time in [1.0e3, 1.0e8, 1.0e16] {
        let result = ReactionPlan::prepare(
            &[1.0, 1.0],
            &[edge(0, 1, 1.0), edge(1, 0, 1.0)],
            time,
            DecayLimits::default(),
            || false,
        );
        match result {
            Ok(prepared) => {
                let answer = prepared.apply(&[10.0, 0.0]).unwrap();
                near(answer.populations[0], 5.0, 4.0e-12);
                near(answer.populations[1], 5.0, 4.0e-12);
                near(answer.removals[0], 5.0 * time + 2.5, 4.0e-12);
                near(answer.removals[1], 5.0 * time - 2.5, 4.0e-12);
            }
            // The bounded method may refuse unresolved extreme cyclic
            // stiffness; it cannot accept drift against a scale of O(time).
            Err(DecayError::Balance { row: 2, .. }) => {}
            Err(error) => panic!("unexpected cyclic refusal: {error}"),
        }
    }
}

#[test]
fn self_return_and_self_multiplication_keep_their_actual_event_counts() {
    let result = plan(&[0.8], &[edge(0, 0, 1.0)], 2.4).apply(&[7.0]).unwrap();
    assert_eq!(result.populations, vec![7.0]);
    near(result.removals[0], 7.0 * 0.8 * 2.4, 2.0e-14);
    let result = plan(&[0.3], &[edge(0, 0, 2.0)], 4.0).apply(&[2.0]).unwrap();
    near(result.populations[0], 2.0 * 1.2_f64.exp(), 2.0e-14);
    near(result.removals[0], 2.0 * 1.2_f64.exp_m1(), 2.0e-14);
}

#[test]
fn unrelated_stiff_component_does_not_erase_slow_removal() {
    let prepared = plan(
        &[1.0e-20, 0.0, 1.0e20, 0.0],
        &[edge(0, 1, 1.0), edge(2, 3, 1.0)],
        3.0e5,
    );
    let result = prepared.apply(&[1.0e30, 0.0, 4.0, 0.0]).unwrap();
    let slow_removals = 1.0e30 * -(-3.0e-15_f64).exp_m1();
    near(result.removals[0], slow_removals, 2.0e-13);
    near(result.populations[1], slow_removals, 2.0e-13);
    assert_eq!(result.populations[2], 0.0);
    near(result.populations[3], 4.0, 2.0e-13);
    near(result.removals[2], 4.0, 2.0e-13);
    assert!(prepared.squarings() > 80);
}

#[test]
fn tiny_removal_is_integrated_without_subtracting_large_populations() {
    let result = plan(&[1.0e-18, 0.0], &[edge(0, 1, 1.0)], 1.0)
        .apply(&[1.0e20, 0.0])
        .unwrap();
    assert_eq!(result.populations[0], 1.0e20);
    near(result.removals[0], 100.0, 2.0e-14);
    near(result.populations[1], 100.0, 2.0e-14);
}

#[test]
fn cyclic_multiplication_semigroup_and_count_composition_agree() {
    let rates = [0.7, 0.2, 0.04, 0.0];
    let edges = [
        edge(0, 1, 1.7),
        edge(1, 0, 0.35),
        edge(1, 2, 1.1),
        edge(2, 3, 1.0),
    ];
    let initial = [11.0, 3.0, 0.0, 0.5];
    let full = plan(&rates, &edges, 2.3).apply(&initial).unwrap();
    let first = plan(&rates, &edges, 0.9).apply(&initial).unwrap();
    let second = plan(&rates, &edges, 1.4).apply(&first.populations).unwrap();
    for i in 0..4 {
        near(full.populations[i], second.populations[i], 3.0e-13);
        if full.removals[i] > 0.0 {
            near(
                full.removals[i],
                first.removals[i] + second.removals[i],
                3.0e-13,
            );
        }
    }
}

#[test]
fn original_node_order_and_reused_plans_preserve_linearity() {
    let rates = [0.0, 0.8, 1.1];
    let edges = [edge(1, 2, 1.5), edge(2, 0, 1.0)];
    let prepared = plan(&rates, &edges, 0.7);
    let one = prepared.apply(&[0.0, 4.0, 2.0]).unwrap();
    let two = prepared.apply(&[0.0, 12.0, 6.0]).unwrap();
    let permuted = plan(&[0.8, 1.1, 0.0], &[edge(0, 1, 1.5), edge(1, 2, 1.0)], 0.7)
        .apply(&[4.0, 2.0, 0.0])
        .unwrap();
    for (i, p) in [2, 0, 1].into_iter().enumerate() {
        near(two.populations[i], 3.0 * one.populations[i], 3.0e-14);
        near(one.populations[i], permuted.populations[p], 3.0e-14);
    }
    assert_eq!(prepared.duration_s(), 0.7);
    assert_eq!(prepared.node_count(), 3);
    assert_eq!(prepared.transition_matrix().len(), 9);
    assert_eq!(prepared.integrated_removal_matrix().len(), 9);
}

#[test]
fn zero_and_stable_intervals_preserve_inventory_without_events() {
    for prepared in [
        plan(&[1.0, 0.0], &[edge(0, 1, 2.0)], 0.0),
        plan(&[0.0, 0.0], &[edge(0, 1, 2.0)], 3.0),
    ] {
        let result = prepared.apply(&[7.0, 3.0]).unwrap();
        assert_eq!(result.populations, vec![7.0, 3.0]);
        assert_eq!(result.removals, vec![0.0, 0.0]);
        assert_eq!(result.balance_residuals, vec![0.0, 0.0]);
    }
}

#[test]
fn exhausted_tails_are_distinguished_from_large_inventory_amplification() {
    let prepared = plan(
        &[1.0, 1.4, 0.0],
        &[edge(0, 1, 2.0), edge(1, 2, 1.0)],
        1000.0,
    );
    let exhausted = prepared.apply(&[1.0e100, 0.0, 0.0]).unwrap();
    assert_eq!(exhausted.populations[0], 0.0);
    assert_eq!(exhausted.populations[1], 0.0);
    near(exhausted.populations[2], 2.0e100, 3.0e-13);
    assert!(matches!(
        prepared.apply(&[1.0e300, 0.0, 0.0]),
        Err(DecayError::Unrepresentable(_))
    ));
    let limits = DecayLimits {
        maximum_series_terms: 256,
        ..DecayLimits::default()
    };
    let tiny =
        ReactionPlan::prepare(&[1.0, 0.0], &[edge(0, 1, 1.0e-310)], 1.0, limits, || false).unwrap();
    assert!(matches!(
        tiny.apply(&[1.0e100, 0.0]),
        Err(DecayError::Unrepresentable(_))
    ));
}

#[test]
fn invalid_data_and_work_scale_series_budgets_refuse() {
    let prepare = |rates: &[f64], edges: &[ReactionEdge], time, limits| {
        ReactionPlan::prepare(rates, edges, time, limits, || false)
    };
    assert!(matches!(
        prepare(&[], &[], 1.0, DecayLimits::default()),
        Err(DecayError::Invalid(_))
    ));
    assert!(prepare(&[-1.0], &[], 1.0, DecayLimits::default()).is_err());
    assert!(prepare(&[1.0], &[edge(0, 1, 1.0)], 1.0, DecayLimits::default()).is_err());
    assert!(prepare(&[1.0], &[edge(0, 0, f64::NAN)], 1.0, DecayLimits::default()).is_err());
    assert!(prepare(&[1.0], &[], -1.0, DecayLimits::default()).is_err());
    assert!(matches!(
        prepare(
            &[1.0, 0.0],
            &[edge(0, 1, 1.0)],
            1.0,
            DecayLimits {
                maximum_work: 1,
                ..DecayLimits::default()
            }
        ),
        Err(DecayError::Budget { .. })
    ));
    assert!(matches!(
        prepare(
            &[1.0],
            &[],
            4.0,
            DecayLimits {
                maximum_squarings: 0,
                ..DecayLimits::default()
            }
        ),
        Err(DecayError::Budget { .. })
    ));
    assert!(matches!(
        prepare(
            &[1.0],
            &[],
            1.0,
            DecayLimits {
                maximum_series_terms: 1,
                ..DecayLimits::default()
            }
        ),
        Err(DecayError::SeriesNotConverged { .. })
    ));
    let prepared = plan(&[1.0], &[], 0.1);
    assert!(prepared.apply(&[-1.0]).is_err());
    assert!(prepared.apply(&[f64::INFINITY]).is_err());
    assert!(prepared.apply(&[]).is_err());
}

#[test]
fn cancellation_is_bounded_transactional_and_replay_is_exact() {
    let rates = [0.7, 1.3];
    let edges = [edge(0, 1, 1.0), edge(1, 0, 1.0)];
    let mut polls = 0;
    let prepared = ReactionPlan::prepare(&rates, &edges, 2.7, DecayLimits::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        let error = ReactionPlan::prepare(&rates, &edges, 2.7, DecayLimits::default(), || {
            seen += 1;
            seen == stop
        })
        .unwrap_err();
        assert_eq!(error, DecayError::Cancelled);
    }
    let initial = [7.0, 2.0];
    let expected = prepared.apply(&initial).unwrap();
    for stop in 1..=4 {
        let mut seen = 0;
        assert_eq!(
            prepared
                .apply_with_cancel(&initial, || {
                    seen += 1;
                    seen == stop
                })
                .unwrap_err(),
            DecayError::Cancelled
        );
        assert_eq!(initial, [7.0, 2.0]);
    }
    assert_eq!(expected, prepared.apply(&initial).unwrap());
}
