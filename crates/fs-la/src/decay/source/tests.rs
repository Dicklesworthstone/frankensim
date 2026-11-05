use super::*;

fn edge(parent: usize, daughter: usize, branch_fraction: f64) -> DecayEdge {
    DecayEdge {
        parent,
        daughter,
        branch_fraction,
    }
}

fn plan(rates: &[f64], edges: &[DecayEdge], time: f64) -> DecaySourcePlan {
    DecaySourcePlan::prepare(rates, edges, time, DecayLimits::default(), || false).unwrap()
}

fn close(actual: f64, expected: f64) {
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
fn constant_source_single_parent_matches_the_independent_inhomogeneous_solution() {
    let x = core::f64::consts::LN_2;
    let plan = plan(&[x, 0.0], &[edge(0, 1, 1.0)], 1.0);
    let answer = plan.apply(&[20.0, 7.0], &[100.0, 0.0]).unwrap();
    let retained = 0.5 / x;
    close(answer.populations[0], 10.0 + 100.0 * retained);
    close(answer.populations[1], 17.0 + 100.0 * (1.0 - retained));
    close(answer.decays[0], 10.0 + 100.0 * (1.0 - retained));
    assert_eq!(answer.decays[1], 0.0);
    close(plan.source_population_matrix()[0], retained);
    close(plan.source_decay_matrix()[0], 1.0 - retained);
    assert!(answer.population_balance_residual.abs() < 1e-12);
    // The source is neither a dose at the start nor a dose at the end.
    assert!(answer.populations[0] > 60.0 && answer.populations[0] < 110.0);
}

#[test]
fn stable_sources_and_zero_intervals_have_explicit_exact_limits() {
    let stable = plan(&[0.0, 0.0], &[], f64::MAX);
    assert_eq!(
        stable.apply(&[2.0, 3.0], &[7.0, 11.0]).unwrap(),
        DecayAdvance {
            populations: vec![9.0, 14.0],
            decays: vec![0.0, 0.0],
            population_balance_residual: 0.0,
        }
    );
    let zero = plan(&[1.0, 0.0], &[edge(0, 1, 1.0)], 0.0);
    assert_eq!(
        zero.apply(&[2.0, 3.0], &[0.0, 0.0]).unwrap().populations,
        [2.0, 3.0]
    );
    assert!(matches!(
        zero.apply(&[2.0, 3.0], &[1.0, 0.0]),
        Err(DecayError::Invalid(_))
    ));
    assert_eq!(stable.homogeneous_plan().series_terms(), 0);
}

#[test]
fn equal_rates_integrate_ingrowth_and_transformations_of_injected_daughters() {
    let x: f64 = 1.4;
    let survival = (-x).exp();
    let parent = (1.0 - survival) / x;
    let daughter = (1.0 - (1.0 + x) * survival) / x;
    let stable = 1.0 - parent - daughter;
    let plan = plan(&[2.0, 2.0, 0.0], &[edge(0, 1, 1.0), edge(1, 2, 1.0)], 0.7);
    let answer = plan.apply(&[0.0, 0.0, 0.0], &[100.0, 0.0, 0.0]).unwrap();
    close(answer.populations[0], 100.0 * parent);
    close(answer.populations[1], 100.0 * daughter);
    close(answer.populations[2], 100.0 * stable);
    close(answer.decays[0], 100.0 * (1.0 - parent));
    close(answer.decays[1], 100.0 * stable);
    assert_eq!(answer.decays[2], 0.0);
}

#[test]
fn multiple_sources_branching_and_arbitrary_order_match_distinct_rate_integrals() {
    let time: f64 = 2.5;
    let a = 0.3;
    let b = 0.7;
    let phi_a = -(-a * time).exp_m1() / (a * time);
    let phi_b = -(-b * time).exp_m1() / (b * time);
    let parent = 30.0 * phi_a;
    let daughter = 30.0 * 0.4 * a / (b - a) * (phi_a - phi_b) + 5.0 * phi_b;
    let plan = plan(
        &[0.0, b, a],
        &[edge(2, 1, 0.4), edge(2, 0, 0.6), edge(1, 0, 1.0)],
        time,
    );
    let answer = plan.apply(&[0.0; 3], &[7.0, 5.0, 30.0]).unwrap();
    close(answer.populations[2], parent);
    close(answer.populations[1], daughter);
    close(answer.populations[0], 42.0 - parent - daughter);
    close(answer.decays[2], 30.0 - parent);
    close(answer.decays[1], 0.4 * (30.0 - parent) + 5.0 - daughter);
    assert_eq!(answer.decays[0], 0.0);
}

#[test]
fn fixed_rate_injection_composes_across_intervals_with_integrated_counts() {
    let rates = [0.3, 0.7, 0.0];
    let edges = [edge(0, 1, 0.6), edge(0, 2, 0.4), edge(1, 2, 1.0)];
    let initial = [12.0, 7.0, 3.0];
    let first = plan(&rates, &edges, 0.7)
        .apply(&initial, &[2.8, 1.4, 0.7])
        .unwrap();
    let second = plan(&rates, &edges, 1.1)
        .apply(&first.populations, &[4.4, 2.2, 1.1])
        .unwrap();
    let whole = plan(&rates, &edges, 1.8)
        .apply(&initial, &[7.2, 3.6, 1.8])
        .unwrap();
    for i in 0..3 {
        close(second.populations[i], whole.populations[i]);
        close(first.decays[i] + second.decays[i], whole.decays[i]);
    }
    let prepared = plan(&rates, &edges, 0.7);
    let unforced = prepared.homogeneous_plan().apply(&initial).unwrap();
    let zero_feed = prepared.apply(&initial, &[0.0; 3]).unwrap();
    assert_eq!(zero_feed.populations, unforced.populations);
    assert_eq!(zero_feed.decays, unforced.decays);
}

#[test]
fn tiny_source_transformations_keep_the_mean_age_and_ignore_unrelated_stiffness() {
    for slow_rate in [1e-20, 1.0] {
        let prepared = plan(
            &[slow_rate, 0.0, 1e100, 0.0],
            &[edge(0, 1, 1.0), edge(2, 3, 1.0)],
            1.0,
        );
        let answer = prepared.apply(&[0.0; 4], &[1.0, 0.0, 0.0, 0.0]).unwrap();
        let lost = if slow_rate < 1e-10 {
            0.5 * slow_rate
        } else {
            1.0 - (-(-slow_rate).exp_m1() / slow_rate)
        };
        close(answer.populations[0], 1.0 - lost);
        close(answer.populations[1], lost);
        close(answer.decays[0], lost);
        assert_eq!(answer.populations[2..], [0.0, 0.0]);
        assert_eq!(answer.decays[1..], [0.0, 0.0, 0.0]);
    }
}

#[test]
fn continuing_sources_cannot_use_the_initial_inventory_exhaustion_exemption() {
    let prepared = plan(&[1e100, 0.0], &[edge(0, 1, 1.0)], 1e300);
    // An old inventory is exhausted, but continuing injection still has an
    // endpoint parent coefficient about 1e-400. A large feed could resolve it.
    assert!(prepared.homogeneous_plan().apply(&[1e200, 0.0]).is_ok());
    assert!(matches!(
        prepared.apply(&[0.0, 0.0], &[1e200, 0.0]),
        Err(DecayError::Unrepresentable(_))
    ));
    let small = prepared.apply(&[0.0, 0.0], &[1.0, 0.0]).unwrap();
    assert_eq!(small.populations[0], 0.0);
    close(small.populations[1], 1.0);
    close(small.decays[0], 1.0);
    let finite = plan(&[1e20, 0.0], &[edge(0, 1, 1.0)], 1.0);
    close(
        finite.apply(&[0.0; 2], &[1e200, 0.0]).unwrap().populations[0],
        1e180,
    );
}

#[test]
fn source_work_budgets_cancellation_and_input_refusals_publish_no_candidate() {
    let rates = [1.0, 0.0];
    let edges = [edge(0, 1, 1.0)];
    let prepared = plan(&rates, &edges, 1.0);
    let limits = DecayLimits {
        maximum_work: prepared.homogeneous_plan().preparation_work() - 1,
        ..DecayLimits::default()
    };
    assert!(matches!(
        DecaySourcePlan::prepare(&rates, &edges, 1.0, limits, || false),
        Err(DecayError::Budget { .. })
    ));
    assert!(matches!(
        DecaySourcePlan::prepare(&rates, &edges, 1.0, DecayLimits::default(), || true),
        Err(DecayError::Cancelled)
    ));
    let initial = [2.0, 3.0];
    let injected = [7.0, 0.0];
    let mut calls = 0;
    prepared
        .apply_with_cancel(&initial, &injected, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut current = 0;
    assert_eq!(
        prepared.apply_with_cancel(&initial, &injected, || {
            current += 1;
            current == calls
        }),
        Err(DecayError::Cancelled)
    );
    assert_eq!(initial, [2.0, 3.0]);
    assert_eq!(injected, [7.0, 0.0]);
    for bad in [
        vec![1.0],
        vec![-1.0, 0.0],
        vec![f64::NAN, 0.0],
        vec![f64::MAX, f64::MAX],
    ] {
        assert!(prepared.apply(&initial, &bad).is_err());
    }
}
