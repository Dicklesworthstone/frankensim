use super::*;
use crate::decay::reaction::{ReactionEdge, source::ReactionSourcePlan};

fn edges() -> [ReactionEdge; 2] {
    [
        ReactionEdge {
            parent: 0,
            daughter: 1,
            yield_per_removal: 1.0,
        },
        ReactionEdge {
            parent: 1,
            daughter: 0,
            yield_per_removal: 1.0,
        },
    ]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn equilibrium_and_repeated_events_remain_resolved_over_stiff_histories() {
    for (a, b) in [(1.0, 1.0), (0.7, 1.3), (1e-8, 1e8)] {
        for time in [1e3, 1e8, 1e16] {
            let limits = DecayLimits {
                maximum_squarings: 0,
                maximum_series_terms: 1,
                maximum_work: 2_000,
                ..DecayLimits::default()
            };
            let plan = ReactionPlan::prepare(&[a, b], &edges(), time, limits, || false).unwrap();
            let answer = plan.apply(&[10.0, 4.0]).unwrap();
            let equilibrium = 14.0 * b / (a + b);
            let deviation = 10.0 - equilibrium;
            let excess = deviation * -(-(a + b) * time).exp_m1() / (a + b);
            near(answer.populations[0], equilibrium, 3e-15);
            near(answer.populations[1], 14.0 * a / (a + b), 3e-15);
            near(answer.removals[0], a * (equilibrium * time + excess), 3e-15);
            near(
                answer.removals[1],
                b * ((14.0 * a / (a + b)) * time - excess),
                3e-15,
            );
            assert!(answer.population_balance_residual.abs() < 1e-14);
            assert_eq!(plan.squarings(), 0);
            assert_eq!(plan.series_terms(), 0);
        }
    }
}

#[test]
fn tiny_exchange_and_uniform_source_ages_keep_their_first_nonzero_terms() {
    let time = 1e-20;
    let population = 1e25;
    let plan = ReactionSourcePlan::prepare_with_population_integrals(
        &[3.0, 7.0],
        &edges(),
        time,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let initial = plan.apply(&[population, 0.0], &[0.0, 0.0]).unwrap();
    near(initial.populations[1], population * 3.0 * time, 2e-15);
    near(initial.removals[0], population * 3.0 * time, 2e-15);
    near(
        initial.removals[1],
        population * (3.0 * time) * (7.0 * time) / 2.0,
        2e-15,
    );
    let injected = plan.apply(&[0.0, 0.0], &[population, 0.0]).unwrap();
    near(
        injected.populations[1],
        population * 3.0 * time / 2.0,
        2e-15,
    );
    near(injected.removals[0], population * 3.0 * time / 2.0, 2e-15);
    near(
        injected.removals[1],
        population * (3.0 * time) * (7.0 * time) / 6.0,
        2e-15,
    );
    let residence = plan
        .integrate_populations(&[0.0, 0.0], &[population, 0.0])
        .unwrap();
    near(residence[0], population * time / 2.0, 2e-15);
    near(residence[1], population * (3.0 * time) * time / 6.0, 2e-15);
}

#[test]
fn source_and_residence_maps_obey_the_same_semigroup() {
    for time in [0.1, 3.0, 1e12] {
        let prepare = |duration| {
            ReactionSourcePlan::prepare_with_population_integrals(
                &[0.7, 1.3],
                &edges(),
                duration,
                DecayLimits::default(),
                || false,
            )
            .unwrap()
        };
        let whole = prepare(time);
        let half = prepare(time / 2.0);
        let initial = [10.0, 4.0];
        let injection = [3.0, 2.0];
        let half_injection = [1.5, 1.0];
        let answer = whole.apply(&initial, &injection).unwrap();
        let first = half.apply(&initial, &half_injection).unwrap();
        let second = half.apply(&first.populations, &half_injection).unwrap();
        let residence = whole.integrate_populations(&initial, &injection).unwrap();
        let residence1 = half
            .integrate_populations(&initial, &half_injection)
            .unwrap();
        let residence2 = half
            .integrate_populations(&first.populations, &half_injection)
            .unwrap();
        for i in 0..2 {
            near(answer.populations[i], second.populations[i], 3e-14);
            near(
                answer.removals[i],
                first.removals[i] + second.removals[i],
                3e-14,
            );
            near(residence[i], residence1[i] + residence2[i], 3e-14);
            near(residence[i] * [0.7, 1.3][i], answer.removals[i], 3e-14);
        }
        near(residence.iter().sum(), time * (14.0 + 5.0 / 2.0), 3e-14);
    }
}

#[test]
fn overflowing_rate_sum_keeps_a_finite_physical_exposure() {
    let plan = ReactionPlan::prepare(
        &[1e308, 1e308],
        &edges(),
        1e-308,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let answer = plan.apply(&[1.0, 0.0]).unwrap();
    let loss = -(-2.0_f64).exp_m1();
    near(answer.populations[0], 1.0 - 0.5 * loss, 2e-15);
    near(answer.populations[1], 0.5 * loss, 2e-15);
    near(answer.removals[0], 0.5 + 0.25 * loss, 2e-15);
    near(answer.removals[1], 0.5 - 0.25 * loss, 2e-15);
}

#[test]
fn subnormal_rates_are_scaled_before_integrating_repeated_events() {
    let rate = f64::from_bits(1);
    let time = 1e308;
    let population = 1e20;
    let plan = ReactionPlan::prepare(
        &[rate, rate],
        &edges(),
        time,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let answer = plan.apply(&[population, 0.0]).unwrap();
    let exposure = rate * time;
    near(answer.populations[1], population * exposure, 2e-15);
    near(
        answer.removals[1],
        population * exposure * exposure / 2.0,
        2e-15,
    );
}

#[test]
fn exact_yields_work_limits_and_late_cancellation_remain_mandatory() {
    let limits = DecayLimits {
        maximum_work: 100,
        ..DecayLimits::default()
    };
    assert!(matches!(
        ReactionPlan::prepare(&[1.0, 2.0], &edges(), 1e10, limits, || false),
        Err(DecayError::Budget {
            quantity: "preparation work",
            ..
        }),
    ));
    let mut polls = 0;
    assert_eq!(
        ReactionPlan::prepare(&[1.0, 2.0], &edges(), 1e10, DecayLimits::default(), || {
            polls += 1;
            polls == 6
        })
        .unwrap_err(),
        DecayError::Cancelled,
    );
    let mut altered = edges();
    altered[0].yield_per_removal = 1.0 + f64::EPSILON;
    let limits = DecayLimits {
        maximum_squarings: 0,
        ..DecayLimits::default()
    };
    assert!(matches!(
        ReactionPlan::prepare(&[1.0, 2.0], &altered, 1.0, limits, || false),
        Err(DecayError::Budget {
            quantity: "reaction squarings",
            ..
        }),
    ));
}
