use super::*;
use crate::decay::reaction::{ReactionEdge, source::ReactionSourcePlan};

fn edges(crossing_fraction: f64) -> Vec<ReactionEdge> {
    let mut result = vec![
        ReactionEdge {
            parent: 0,
            daughter: 1,
            yield_per_removal: 1.0,
        },
        ReactionEdge {
            parent: 1,
            daughter: 0,
            yield_per_removal: 1.0 - crossing_fraction,
        },
        ReactionEdge {
            parent: 2,
            daughter: 0,
            yield_per_removal: 1.0,
        },
    ];
    if crossing_fraction > 0.0 {
        result.push(ReactionEdge {
            parent: 1,
            daughter: 2,
            yield_per_removal: crossing_fraction,
        });
    }
    result
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn g1_weak_pumping_preserves_nanosecond_exchange_on_the_thermal_clock() {
    let decay = 1e9;
    for pump in [1.0, 10.0, 1e4] {
        for time in [1e-6, 1e-3, 1.0] {
            let plan = ReactionPlan::prepare(
                &[pump, decay, 5e4],
                &edges(0.0),
                time,
                DecayLimits::default(),
                || false,
            )
            .unwrap();
            let answer = plan.apply(&[1.0, 0.0, 0.0]).unwrap();
            let rate = pump + decay;
            let loss = -(-rate * time).exp_m1();
            let excited = pump / rate * loss;
            let repeated = (pump * decay / rate) * time;
            near(answer.populations[0], 1.0 - excited, 3e-14);
            near(answer.populations[1], excited, 3e-14);
            assert_eq!(answer.populations[2], 0.0);
            near(
                answer.removals[0],
                repeated + (pump / rate).powi(2) * loss,
                3e-14,
            );
            near(
                answer.removals[1],
                repeated - (pump / rate) * (decay / rate) * loss,
                3e-14,
            );
            assert_eq!(answer.removals[2], 0.0);
            assert!(answer.population_balance_residual.abs() < 1e-15);
        }
    }
}

#[test]
fn g1_stiff_three_state_source_maps_match_independent_block_exponential() {
    // Rounded from a 90-decimal-digit exponential of the 9x9 block matrix
    // [A,0,I; diag(rates),0,0; 0,0,0] at t=1s. The population, event and
    // source blocks are evaluated independently of the production series.
    let rates = [1.0, 1e9, 5e4];
    let expected_p = [
        0.9999849992250266,
        0.9999849992250266,
        0.9999849992250266,
        9.999849992250267e-10,
        9.999849992250267e-10,
        9.999849992250267e-10,
        1.49997749883754e-5,
        1.49997749883754e-5,
        1.49997749883754e-5,
    ];
    let expected_c = [
        0.9999849995250326,
        0.9999699987500592,
        0.9999649998250482,
        0.9999849985250476,
        1.9999699977500742,
        0.9999649988250632,
        0.7499737491187973,
        1.4999624985375673,
        1.749958749343809,
    ];
    let expected_f = [
        0.9999849995250326,
        0.9999699987500592,
        0.9999649998250482,
        9.999849985250477e-10,
        1.999969997750074e-9,
        9.999649988250632e-10,
        1.4999474982375947e-5,
        2.9999249970751345e-5,
        3.499917498687618e-5,
    ];
    let expected_g = [
        0.4999924999125133,
        0.49997749943754594,
        0.4999725006125168,
        0.49999249891252834,
        1.499977497437576,
        0.4999724996125518,
        0.37497937470941384,
        1.124953123828211,
        1.374944375534427,
    ];
    for permutation in [[0, 1, 2], [2, 0, 1]] {
        let mut permuted_rates = [0.0; N];
        for i in 0..N {
            permuted_rates[permutation[i]] = rates[i];
        }
        let permuted_edges: Vec<_> = edges(0.75)
            .into_iter()
            .map(|edge| ReactionEdge {
                parent: permutation[edge.parent],
                daughter: permutation[edge.daughter],
                yield_per_removal: edge.yield_per_removal,
            })
            .collect();
        let plan = ReactionSourcePlan::prepare_with_population_integrals(
            &permuted_rates,
            &permuted_edges,
            1.0,
            DecayLimits::default(),
            || false,
        )
        .unwrap();
        let homogeneous = plan.homogeneous_plan();
        for (actual, expected) in [
            (homogeneous.transition_matrix(), expected_p),
            (homogeneous.integrated_removal_matrix(), expected_c),
            (plan.source_population_matrix(), expected_f),
            (plan.source_removal_matrix(), expected_g),
        ] {
            for i in 0..N {
                for j in 0..N {
                    near(
                        actual[permutation[i] * N + permutation[j]],
                        expected[i * N + j],
                        3e-14,
                    );
                }
            }
        }
        for i in 0..CELLS {
            near(
                plan.population_integral_matrix().unwrap()[i] * permuted_rates[i / N],
                homogeneous.integrated_removal_matrix()[i],
                3e-14,
            );
            near(
                plan.source_population_integral_matrix().unwrap()[i] * permuted_rates[i / N],
                plan.source_removal_matrix()[i],
                3e-14,
            );
        }
    }
}

#[test]
fn g3_stiff_sources_counts_and_residence_compose_without_population_correction() {
    let rates = [1.0, 1e9, 5e4];
    let edges = edges(0.75);
    let prepare = |time| {
        ReactionSourcePlan::prepare_with_population_integrals(
            &rates,
            &edges,
            time,
            DecayLimits::default(),
            || false,
        )
        .unwrap()
    };
    let whole = prepare(1.0);
    let first_plan = prepare(0.3);
    let second_plan = prepare(0.7);
    let initial = [10.0, 4.0, 3.0];
    let injected = [1.0, 2.0, 5.0];
    let first_injection = injected.map(|value| value * 0.3);
    let second_injection = injected.map(|value| value * 0.7);
    let answer = whole.apply(&initial, &injected).unwrap();
    let first = first_plan.apply(&initial, &first_injection).unwrap();
    let second = second_plan
        .apply(&first.populations, &second_injection)
        .unwrap();
    let residence = whole.integrate_populations(&initial, &injected).unwrap();
    let residence1 = first_plan
        .integrate_populations(&initial, &first_injection)
        .unwrap();
    let residence2 = second_plan
        .integrate_populations(&first.populations, &second_injection)
        .unwrap();
    for i in 0..N {
        near(answer.populations[i], second.populations[i], 5e-14);
        near(
            answer.removals[i],
            first.removals[i] + second.removals[i],
            5e-14,
        );
        near(residence[i], residence1[i] + residence2[i], 5e-14);
    }
    near(answer.populations.iter().sum(), 25.0, 3e-15);
    assert!(answer.population_balance_residual.abs() < 1e-14);
}

#[test]
fn g0_stiff_small_blocks_retain_work_series_and_late_cancellation_refusals() {
    let rates = [1.0, 1e9, 5e4];
    let edges = edges(0.75);
    assert!(matches!(
        ReactionPlan::prepare(
            &rates,
            &edges,
            1.0,
            DecayLimits {
                maximum_work: 1_000,
                ..DecayLimits::default()
            },
            || false
        ),
        Err(DecayError::Budget {
            quantity: "preparation work",
            ..
        }),
    ));
    assert!(matches!(
        ReactionPlan::prepare(
            &rates,
            &edges,
            1.0,
            DecayLimits {
                maximum_series_terms: 1,
                ..DecayLimits::default()
            },
            || false
        ),
        Err(DecayError::SeriesNotConverged { terms: 1 }),
    ));
    let mut polls = 0;
    ReactionPlan::prepare(&rates, &edges, 1.0, DecayLimits::default(), || {
        polls += 1;
        false
    })
    .unwrap();
    let cancel_at = polls - 1;
    polls = 0;
    assert_eq!(
        ReactionPlan::prepare(&rates, &edges, 1.0, DecayLimits::default(), || {
            polls += 1;
            polls == cancel_at
        })
        .unwrap_err(),
        DecayError::Cancelled,
    );
}
