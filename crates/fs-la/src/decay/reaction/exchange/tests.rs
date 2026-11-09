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

fn passive_edges(p: f64, q: f64) -> [ReactionEdge; 2] {
    [
        ReactionEdge {
            parent: 0,
            daughter: 1,
            yield_per_removal: p,
        },
        ReactionEdge {
            parent: 1,
            daughter: 0,
            yield_per_removal: q,
        },
    ]
}

#[test]
fn damped_pair_matches_independent_augmented_matrix_exponential() {
    // References were rounded from a 90-decimal-digit exponential of the
    // six-state block matrix [A,0,I; diag(rates),0,0; 0,0,0]. Its four
    // blocks independently produce P, C, Fbar and Gbar, without using the
    // eigenprojectors or scalar moments implemented in the production path.
    let cases = [
        (
            [1e9, 0.4],
            [0.75, 0.4],
            0.5,
            [
                [
                    1.0432298829551144e-10,
                    1.3909731768840133e-10,
                    0.6520186766643813,
                    0.8693582353090881,
                ],
                [
                    1.0559893274713208,
                    0.07465243676752495,
                    0.13997331893910928,
                    0.18663109226655566,
                ],
                [
                    2.1119786549426415e-9,
                    1.493048735350499e-10,
                    0.6998665946955464,
                    0.9331554613327783,
                ],
                [
                    1.0286476571568612,
                    0.03819687902511972,
                    0.07161914817209948,
                    0.09549219793606149,
                ],
            ],
        ),
        (
            [12.0, 2.0],
            [0.4, 0.7],
            0.25,
            [
                [
                    0.07998436328875706,
                    0.08300581706491703,
                    0.2845913727940012,
                    0.6728830566095929,
                ],
                [
                    1.001113438549225,
                    0.2027445045949554,
                    0.11585400262568882,
                    0.4082147452283892,
                ],
                [
                    0.33370447951640836,
                    0.06758150153165181,
                    0.23170800525137764,
                    0.8164294904567784,
                ],
                [
                    0.7001387733439268,
                    0.08460813215083791,
                    0.048347504086193094,
                    0.21741376240355675,
                ],
            ],
        ),
    ];
    for (rates, yields, duration, expected) in cases {
        for swapped in [false, true] {
            let rates = if swapped { [rates[1], rates[0]] } else { rates };
            let yields = if swapped {
                [yields[1], yields[0]]
            } else {
                yields
            };
            let limits = DecayLimits {
                maximum_squarings: 0,
                maximum_series_terms: 1,
                maximum_work: 2_000,
                ..DecayLimits::default()
            };
            let plan = ReactionSourcePlan::prepare_with_population_integrals(
                &rates,
                &passive_edges(yields[0], yields[1]),
                duration,
                limits,
                || false,
            )
            .unwrap();
            let homogeneous = plan.homogeneous_plan();
            let actual = [
                homogeneous.transition_matrix(),
                homogeneous.integrated_removal_matrix(),
                plan.source_population_matrix(),
                plan.source_removal_matrix(),
            ];
            for (actual, expected) in actual.into_iter().zip(expected) {
                for i in 0..4 {
                    near(actual[i], expected[if swapped { 3 - i } else { i }], 3e-14);
                }
            }
            for i in 0..4 {
                near(
                    plan.population_integral_matrix().unwrap()[i] * rates[i / 2],
                    homogeneous.integrated_removal_matrix()[i],
                    3e-14,
                );
                near(
                    plan.source_population_integral_matrix().unwrap()[i] * rates[i / 2],
                    plan.source_removal_matrix()[i],
                    3e-14,
                );
            }
            assert_eq!(homogeneous.squarings(), 0);
            assert_eq!(homogeneous.series_terms(), 0);
        }
    }
}

#[test]
fn damped_pair_retains_tiny_authored_sinks_over_long_histories() {
    let p = 1.0 - 2.0_f64.powi(-27);
    let q = 1.0 - 2.0_f64.powi(-28);
    let plan = ReactionPlan::prepare(
        &[1e9, 0.4],
        &passive_edges(p, q),
        1e8,
        DecayLimits::default(),
        || false,
    )
    .unwrap();
    let answer = plan.apply(&[0.0, 1.0]).unwrap();
    // Independent 90-digit block-exponential reference at the exact binary64
    // input values. Forming 1-p*q would lose the small product of the two
    // authored sinks and visibly change this long-time survival probability.
    near(answer.populations[0], 2.5580865112618294e-10, 3e-14);
    near(answer.populations[1], 0.639521630197861, 3e-14);
    near(answer.removals[0], 3.225505846238208e7, 3e-14);
    near(answer.removals[1], 3.2255058582541537e7, 3e-14);
    let lost = (1.0 - p) * answer.removals[0] + (1.0 - q) * answer.removals[1];
    near(lost, 0.36047836954633034, 3e-14);
    near(answer.populations.iter().sum::<f64>() + lost, 1.0, 3e-15);
    assert!(answer.population_balance_residual.abs() < 3e-15);
}

#[test]
fn damped_source_and_residence_obey_interval_composition() {
    for rates in [[1e9, 0.4], [0.4, 1e9], [12.0, 2.0]] {
        for time in [1e-9, 0.5, 3.0] {
            let prepare = |duration| {
                ReactionSourcePlan::prepare_with_population_integrals(
                    &rates,
                    &passive_edges(0.75, 0.4),
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
                near(answer.populations[i], second.populations[i], 3e-13);
                near(
                    answer.removals[i],
                    first.removals[i] + second.removals[i],
                    3e-13,
                );
                near(residence[i], residence1[i] + residence2[i], 3e-13);
                near(residence[i] * rates[i], answer.removals[i], 3e-13);
            }
            assert!(answer.population_balance_residual.abs() < 5e-14);
        }
    }
}

#[test]
fn damped_pair_preserves_budgets_and_late_cancellation() {
    let rates = [1e9, 0.4];
    let edges = passive_edges(0.75, 0.4);
    assert!(matches!(
        ReactionPlan::prepare(
            &rates,
            &edges,
            0.5,
            DecayLimits {
                maximum_work: 100,
                ..DecayLimits::default()
            },
            || false,
        ),
        Err(DecayError::Budget {
            quantity: "preparation work",
            ..
        }),
    ));
    let mut polls = 0;
    ReactionPlan::prepare(&rates, &edges, 0.5, DecayLimits::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    let cancel_at = polls - 1;
    polls = 0;
    assert_eq!(
        ReactionPlan::prepare(&rates, &edges, 0.5, DecayLimits::default(), || {
            polls += 1;
            polls == cancel_at
        })
        .unwrap_err(),
        DecayError::Cancelled,
    );
}
