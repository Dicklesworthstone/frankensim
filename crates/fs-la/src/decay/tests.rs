use super::*;

fn edge(parent: usize, daughter: usize, branch_fraction: f64) -> DecayEdge {
    DecayEdge { parent, daughter, branch_fraction }
}
fn prepare(rates: &[f64], edges: &[DecayEdge], time: f64) -> DecayPlan {
    DecayPlan::prepare(rates, edges, time, DecayLimits::default(), || false).unwrap()
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs(),
        "actual {actual:.17e}, expected {expected:.17e}, relative error {:.5e}",
        (actual - expected) / expected);
}

#[test]
fn single_parent_half_life_and_stable_initial_inventory() {
    let plan = prepare(&[core::f64::consts::LN_2, 0.0], &[edge(0, 1, 1.0)], 1.0);
    let answer = plan.apply(&[100.0, 7.0]).unwrap();
    relative(answer.populations[0], 50.0, 2e-14);
    relative(answer.populations[1], 57.0, 2e-14);
    relative(answer.decays[0], 50.0, 2e-14);
    assert_eq!(answer.decays[1], 0.0);
    assert!(answer.population_balance_residual.abs() < 2e-13);
}

#[test]
fn equal_rates_include_daughter_transformations_in_the_same_interval() {
    let time = 0.7;
    let x: f64 = 2.0 * time;
    let parent = (-x).exp();
    let daughter = x * parent;
    let stable = 1.0 - (1.0 + x) * parent;
    let plan = prepare(&[2.0, 2.0, 0.0], &[edge(0, 1, 1.0), edge(1, 2, 1.0)], time);
    let answer = plan.apply(&[1.0, 0.0, 0.0]).unwrap();
    relative(answer.populations[0], parent, 2e-14);
    relative(answer.populations[1], daughter, 2e-14);
    relative(answer.populations[2], stable, 2e-14);
    relative(answer.decays[0], 1.0 - parent, 2e-14);
    relative(answer.decays[1], stable, 2e-14);
    assert_eq!(answer.decays[2], 0.0);
}

#[test]
fn nearly_equal_rates_do_not_subtract_bateman_exponentials() {
    let a = 1.0;
    let b = 1.0 + 1e-11;
    let time: f64 = 2.0;
    let parent = (-a * time).exp();
    let difference = b - a;
    let daughter = a * parent * (-(-difference * time).exp_m1()) / difference;
    let stable = 1.0 - parent - daughter;
    let answer = prepare(&[a, b, 0.0], &[edge(0, 1, 1.0), edge(1, 2, 1.0)], time)
        .apply(&[1.0, 0.0, 0.0]).unwrap();
    relative(answer.populations[1], daughter, 5e-14);
    relative(answer.populations[2], stable, 5e-14);
    relative(answer.decays[1], stable, 5e-14);
}

#[test]
fn arbitrary_node_order_and_branching_retain_original_addresses() {
    let time: f64 = 0.8;
    let p = (-time).exp();
    let d = (-2.0 * time).exp();
    let edges = [edge(3, 1, 0.4), edge(3, 2, 0.6), edge(1, 0, 1.0)];
    let plan = prepare(&[0.0, 2.0, 0.0, 1.0], &edges, time);
    let answer = plan.apply(&[0.0, 7.0, 5.0, 100.0]).unwrap();
    let daughter_decays = 7.0 * (1.0 - d) + 40.0 * (1.0 - 2.0 * p + d);
    relative(answer.populations[0], daughter_decays, 5e-14);
    relative(answer.populations[1], 7.0 * d + 40.0 * (p - d), 5e-14);
    relative(answer.populations[2], 5.0 + 60.0 * (1.0 - p), 5e-14);
    relative(answer.populations[3], 100.0 * p, 5e-14);
    relative(answer.decays[1], daughter_decays, 5e-14);
    relative(answer.decays[3], 100.0 * (1.0 - p), 5e-14);
    assert_eq!((answer.decays[0], answer.decays[2]), (0.0, 0.0));
    assert_eq!(plan.transition_matrix()[3 * 4 + 1], 0.0);
    assert!(plan.transition_matrix()[4 + 3] > 0.0);
}

#[test]
fn distinct_physical_branches_may_share_the_same_daughter() {
    let split = prepare(&[2.0, 0.0], &[edge(0, 1, 0.25), edge(0, 1, 0.75)], 0.4);
    let combined = prepare(&[2.0, 0.0], &[edge(0, 1, 1.0)], 0.4);
    let split = split.apply(&[17.0, 3.0]).unwrap();
    let combined = combined.apply(&[17.0, 3.0]).unwrap();
    relative(split.populations[0], combined.populations[0], 2e-14);
    relative(split.populations[1], combined.populations[1], 2e-14);
    relative(split.decays[0], combined.decays[0], 2e-14);
}

#[test]
fn tiny_transformations_survive_a_rounded_unchanged_parent() {
    for fast_rate in [1e-20, 1e20] {
        // The second disconnected chain changes the preparation scale only.
        let plan = prepare(&[1e-20, 0.0, fast_rate, 0.0],
            &[edge(0, 1, 1.0), edge(2, 3, 1.0)], 1.0);
        let answer = plan.apply(&[1.0, 0.0, 0.0, 0.0]).unwrap();
        assert_eq!(answer.populations[0], 1.0);
        relative(answer.populations[1], 1e-20, 5e-14);
        relative(answer.decays[0], 1e-20, 5e-14);
        assert_eq!(answer.decays[2], 0.0);
    }
}

#[test]
fn stiff_chain_preserves_the_slow_physical_clock() {
    for fast_rate in [1e20, 1e100] {
        let time: f64 = 0.7;
        let plan = prepare(&[fast_rate, 1.0, 0.0],
            &[edge(0, 1, 1.0), edge(1, 2, 1.0)], time);
        assert!(plan.squarings() > 60);
        assert!(plan.series_terms() <= 128);
        let answer = plan.apply(&[1.0, 0.0, 0.0]).unwrap();
        assert_eq!(answer.populations[0], 0.0);
        relative(answer.populations[1], (-time).exp(), 2e-13);
        relative(answer.populations[2], -(-time).exp_m1(), 2e-13);
        relative(answer.decays[0], 1.0, 2e-14);
        relative(answer.decays[1], -(-time).exp_m1(), 2e-13);
    }
}

#[test]
fn long_chain_matches_the_independent_erlang_population_law() {
    let n = MAX_DECAY_NODES;
    let mut rates = vec![1.0; n];
    rates[n - 1] = 0.0;
    let edges: Vec<_> = (0..n - 1).map(|i| edge(i, i + 1, 1.0)).collect();
    let time: f64 = 30.0;
    let plan = prepare(&rates, &edges, time);
    let mut initial = vec![0.0; n]; initial[0] = 1.0;
    let answer = plan.apply(&initial).unwrap();
    let mut poisson = (-time).exp();
    let mut cumulative = 0.0;
    for i in 0..n - 1 {
        if i != 0 { poisson *= time / i as f64; }
        relative(answer.populations[i], poisson, 3e-13);
        cumulative += poisson;
        if i < 25 { relative(answer.decays[i], 1.0 - cumulative, 3e-13); }
    }
    // The final tail is small; an absolute comparison avoids trusting digits
    // lost in this independent oracle's final subtraction from one.
    assert!((answer.populations[n - 1] - (1.0 - cumulative)).abs() < 2e-14);
    assert_eq!(answer.decays[n - 1], 0.0);
}

#[test]
fn independent_intervals_compose_both_populations_and_counts() {
    let rates = [0.3, 0.7, 4.0, 0.0];
    let edges = [edge(0, 1, 0.6), edge(0, 2, 0.4), edge(1, 2, 0.2),
        edge(1, 3, 0.8), edge(2, 3, 1.0)];
    let initial = [8.0, 2.0, 1.0, 4.0];
    let first = prepare(&rates, &edges, 0.7).apply(&initial).unwrap();
    let second = prepare(&rates, &edges, 1.1).apply(&first.populations).unwrap();
    let whole = prepare(&rates, &edges, 1.8).apply(&initial).unwrap();
    for i in 0..4 {
        relative(second.populations[i], whole.populations[i], 8e-14);
        if i < 3 { relative(first.decays[i] + second.decays[i], whole.decays[i], 8e-14); }
        else { assert_eq!(whole.decays[i], 0.0); }
    }
}

#[test]
fn zero_duration_and_all_stable_are_identity_maps() {
    for plan in [prepare(&[1.0, 0.0], &[edge(0, 1, 1.0)], 0.0),
        prepare(&[0.0, 0.0], &[], 1e308)]
    {
        assert_eq!(plan.apply(&[3.0, 8.0]).unwrap(), DecayAdvance {
            populations: vec![3.0, 8.0], decays: vec![0.0, 0.0], population_balance_residual: 0.0,
        });
        assert_eq!(plan.series_terms(), 0);
    }
}

#[test]
fn incomplete_cyclic_invalid_or_unrepresentable_networks_refuse() {
    let check = |rates: &[f64], edges: &[DecayEdge], time| {
        assert!(DecayPlan::prepare(rates, edges, time, DecayLimits::default(), || false).is_err());
    };
    check(&[1.0, 0.0], &[edge(0, 1, 0.5)], 1.0);
    check(&[1.0, 1.0], &[edge(0, 1, 1.0), edge(1, 0, 1.0)], 1.0);
    check(&[1.0], &[edge(0, 0, 1.0)], 1.0);
    check(&[0.0, 0.0], &[edge(0, 1, 1.0)], 1.0);
    check(&[1.0, 0.0], &[edge(0, 2, 1.0)], 1.0);
    check(&[f64::NAN, 0.0], &[edge(0, 1, 1.0)], 1.0);
    check(&[-1.0, 0.0], &[edge(0, 1, 1.0)], 1.0);
    check(&[1.0, 0.0], &[edge(0, 1, 1.0)], -1.0);
    check(&[1e308, 1e-308, 0.0], &[edge(0, 1, 1.0), edge(1, 2, 1.0)], 1.0);
    check(&[1.0, 0.0], &[edge(0, 1, 1e-320), edge(0, 1, 1.0)], 1e-10);
}

#[test]
fn explicit_budgets_and_late_cancellation_refuse_transactionally() {
    let rates = [1.0, 0.0]; let edges = [edge(0, 1, 1.0)];
    for limits in [
        DecayLimits { maximum_nodes: 1, ..DecayLimits::default() },
        DecayLimits { maximum_squarings: 0, ..DecayLimits::default() },
        DecayLimits { maximum_work: 1, ..DecayLimits::default() },
    ] {
        assert!(matches!(DecayPlan::prepare(&rates, &edges, 1.0, limits, || false),
            Err(DecayError::Budget { .. })));
    }
    assert!(matches!(DecayPlan::prepare(&rates, &edges, 1.0,
        DecayLimits { maximum_series_terms: 1, ..DecayLimits::default() }, || false),
        Err(DecayError::SeriesNotConverged { .. })));
    let mut calls = 0;
    assert!(matches!(DecayPlan::prepare(&rates, &edges, 1.0, DecayLimits::default(), || {
        calls += 1; calls == 12
    }), Err(DecayError::Cancelled)));
    let plan = prepare(&rates, &edges, 1.0);
    let initial = [2.0, 3.0]; let mut calls = 0;
    assert_eq!(plan.apply_with_cancel(&initial, || { calls += 1; calls == 3 }),
        Err(DecayError::Cancelled));
    assert_eq!(initial, [2.0, 3.0]);
    assert!(plan.apply(&[-1.0, 1.0]).is_err());
    assert!(plan.apply(&[f64::MAX, f64::MAX]).is_err());
    let again = prepare(&rates, &edges, 1.0);
    assert_eq!(again.transition_matrix(), plan.transition_matrix());
    assert_eq!(again.integrated_decay_matrix(), plan.integrated_decay_matrix());
}

#[test]
fn a_large_inventory_cannot_silently_resurrect_an_underflowed_tail() {
    for time in [730.0, 800.0] {
        let plan = prepare(&[1.0, 0.0], &[edge(0, 1, 1.0)], time);
        assert!(matches!(plan.apply(&[1e200, 0.0]), Err(DecayError::Unrepresentable(_))));
    }
    let chain = prepare(&[1.0, 1.0, 0.0], &[edge(0, 1, 1.0), edge(1, 2, 1.0)], 800.0);
    assert!(matches!(chain.apply(&[1e200, 0.0, 0.0]), Err(DecayError::Unrepresentable(_))));
    // Truly exhausted fast parents remain usable at large, ordinary inventories.
    let exhausted = prepare(&[1e20, 0.0], &[edge(0, 1, 1.0)], 1.0);
    let answer = exhausted.apply(&[1e200, 0.0]).unwrap();
    assert_eq!(answer.populations[0], 0.0);
    relative(answer.populations[1], 1e200, 5e-14);
    relative(answer.decays[0], 1e200, 5e-14);
}

#[test]
fn large_inventory_refuses_a_vanished_multistep_transformation_coefficient() {
    let chain = prepare(&[1.0, 1.0, 1.0, 1.0, 0.0],
        &[edge(0, 1, 1.0), edge(1, 2, 1.0), edge(2, 3, 1.0), edge(3, 4, 1.0)], 1e-100);
    // Four transformations have probability about 1e-400/24, which does not
    // fit in the map, but this population would make their count representable.
    assert!(matches!(chain.apply(&[1e300, 0.0, 0.0, 0.0, 0.0]),
        Err(DecayError::Unrepresentable(_))));
}
