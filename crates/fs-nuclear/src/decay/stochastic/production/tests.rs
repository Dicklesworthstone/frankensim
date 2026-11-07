//! Synthetic immigration/decay moments, source ages, integer balances and replay.

use super::*;
use crate::decay::production::DecayProductionChange;
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_rand::{Stream, StreamCheckpoint, StreamKey};
use fs_thermochem::decay::DecayLawV1;

fn key(seed: u64) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_0002,
        tile: 19,
    }
}

fn budget() -> StochasticDecayProductionBudget {
    StochasticDecayProductionBudget {
        decay: StochasticDecayBudget {
            maximum_events: 10_000,
            energy_tolerance_j: 1e-8,
        },
        maximum_births: 10_000,
        maximum_schedule_changes: 64,
    }
}

fn node(id: &str, rate: f64) -> DecayNuclide {
    DecayNuclide {
        id: id.into(),
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    }
}

fn branch(parent: usize, daughter: usize, probability: f64, deposited_j: f64) -> DecayBranch {
    DecayBranch {
        parent,
        daughter,
        probability,
        energy_per_decay: EnergyChannels {
            deposited_j,
            radiation_j: 0.5 * deposited_j,
            escaped_j: 0.25 * deposited_j,
        },
    }
}

fn model(nodes: Vec<DecayNuclide>, branches: Vec<DecayBranch>) -> Arc<DecayModel> {
    Arc::new(
        DecayModel::new(
            "synthetic-poisson-supply".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn single() -> Arc<DecayModel> {
    model(
        vec![node("parent", 0.7), node("stable", 0.0)],
        vec![branch(0, 1, 1.0, 2.0)],
    )
}

fn schedule(
    model: &Arc<DecayModel>,
    initial: Vec<f64>,
    changes: Vec<(f64, Vec<f64>)>,
) -> Arc<DecayProductionSchedule> {
    Arc::new(
        DecayProductionSchedule::new(
            Arc::clone(model),
            initial,
            changes
                .into_iter()
                .map(|(time_s, rates_nuclei_s)| DecayProductionChange {
                    time_s,
                    rates_nuclei_s,
                })
                .collect(),
            64,
        )
        .unwrap(),
    )
}

fn supplied(
    model: &Arc<DecayModel>,
    initial: Vec<u64>,
    schedule: &Arc<DecayProductionSchedule>,
    seed: u64,
) -> StochasticDecayProduction {
    StochasticDecayProduction::poisson(
        StochasticDecayInventory::new(Arc::clone(model), initial, key(seed)).unwrap(),
        Arc::clone(schedule),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.4e}"
    );
}

fn poisson_moments(sum: f64, squares: f64, samples: u64, expected: f64) {
    let n = samples as f64;
    let mean = sum / n;
    close(mean, expected, 6.0 * (expected / n).sqrt());
    close(
        squares / n - mean * mean,
        expected,
        6.0 * ((expected + 2.0 * expected * expected) / n).sqrt(),
    );
}

#[test]
fn g1_stable_supply_has_poisson_counts_variance_and_independent_nuclides() {
    let model = model(vec![node("stable-a", 0.0), node("stable-b", 0.0)], vec![]);
    let schedule = schedule(&model, vec![3.0, 1.5], vec![]);
    let samples = 8192;
    let mut sums = [0.0; 2];
    let mut squares = [0.0; 2];
    let mut product = 0.0;
    for seed in 0..samples {
        let mut state = supplied(&model, vec![0, 0], &schedule, seed);
        let report = state.advance_to(2.0, budget(), || Ok(())).unwrap();
        assert_eq!(report.produced_nuclei, state.inventory.populations());
        assert_eq!(report.produced_nuclei, state.cumulative_produced_nuclei());
        assert_eq!(report.births, report.produced_nuclei.iter().sum::<u64>());
        assert_eq!(report.decay.events, 0);
        assert_eq!(report.decay.released, EnergyChannels::default());
        for ((sum, square), &count) in sums
            .iter_mut()
            .zip(&mut squares)
            .zip(state.inventory.populations())
        {
            *sum += count as f64;
            *square += (count * count) as f64;
        }
        product += (report.produced_nuclei[0] * report.produced_nuclei[1]) as f64;
    }
    poisson_moments(sums[0], squares[0], samples, 6.0);
    poisson_moments(sums[1], squares[1], samples, 3.0);
    let n = samples as f64;
    close(
        product / n - sums[0] * sums[1] / (n * n),
        0.0,
        6.0 * (18.0 / n).sqrt(),
    );
}

#[test]
fn g1_continuous_arrivals_and_decay_match_poisson_thinning_and_shared_expectations() {
    let model = single();
    let rate = 4.0;
    let end = 2.5;
    let schedule = schedule(&model, vec![rate, 0.0], vec![]);
    let parent = rate / 0.7 * (1.0 - (-0.7_f64 * end).exp());
    let expected = [parent, rate * end - parent];
    let mut deterministic = DecayInventory::new(Arc::clone(&model), vec![0.0, 0.0]).unwrap();
    let deterministic_step = deterministic
        .advance_with_production(
            end,
            &[rate, 0.0],
            DecayBudget {
                propagation: DecayLimits::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    for (&actual, expected) in deterministic.populations().iter().zip(expected) {
        close(actual, expected, 2e-12);
    }
    let samples = 4096;
    let mut sums = [0.0; 2];
    let mut squares = [0.0; 2];
    for seed in 0..samples {
        let mut state = supplied(&model, vec![0, 0], &schedule, seed + 10_000);
        let report = state.advance_to(end, budget(), || Ok(())).unwrap();
        assert_eq!(
            state.inventory.populations().iter().sum::<u64>(),
            report.births
        );
        assert_eq!(report.produced_nuclei, [report.births, 0]);
        assert_eq!(report.decay.events, state.inventory.populations()[1]);
        assert_eq!(
            report.decay.released.deposited_j,
            2.0 * report.decay.events as f64
        );
        assert_eq!(
            report.decay.released.radiation_j,
            report.decay.events as f64
        );
        assert_eq!(
            report.decay.released.escaped_j,
            0.5 * report.decay.events as f64
        );
        for ((sum, square), &count) in sums
            .iter_mut()
            .zip(&mut squares)
            .zip(state.inventory.populations())
        {
            *sum += count as f64;
            *square += (count * count) as f64;
        }
    }
    for ((sum, squares), mean) in sums.into_iter().zip(squares).zip(expected) {
        poisson_moments(sum, squares, samples, mean);
    }
    close(deterministic_step.decay.decays[0], expected[1], 2e-12);
}

#[test]
fn g1_newborn_equal_rate_daughters_decay_with_original_branch_and_energy_correlations() {
    let rate = 1.25;
    let supply = 6.0;
    let end = 1.2;
    let model = model(
        vec![
            node("parent", rate),
            node("daughter", rate),
            node("stable-a", 0.0),
            node("stable-b", 0.0),
        ],
        vec![
            branch(0, 1, 0.25, 3.0),
            branch(0, 2, 0.75, 7.0),
            branch(1, 3, 1.0, 5.0),
        ],
    );
    let schedule = schedule(&model, vec![supply, 0.0, 0.0, 0.0], vec![]);
    let x: f64 = rate * end;
    let parent = supply / rate * (1.0 - (-x).exp());
    let daughter = 0.25 * supply / rate * (1.0 - (1.0 + x) * (-x).exp());
    let parent_decays = supply * end - parent;
    let expected = [
        parent,
        daughter,
        0.75 * parent_decays,
        0.25 * parent_decays - daughter,
    ];
    let samples = 4096;
    let mut sums = [0.0; 4];
    let mut squares = [0.0; 4];
    for seed in 0..samples {
        let mut state = supplied(&model, vec![0; 4], &schedule, seed + 20_000);
        let report = state.advance_to(end, budget(), || Ok(())).unwrap();
        let counts = state.inventory.populations();
        let branches = &report.decay.branch_decays;
        assert_eq!(counts.iter().sum::<u64>(), report.births);
        assert_eq!(counts[1], branches[0] - branches[2]);
        assert_eq!(counts[2], branches[1]);
        assert_eq!(counts[3], branches[2]);
        assert_eq!(
            report.decay.decays,
            [branches[0] + branches[1], branches[2], 0, 0]
        );
        let heat = (3 * branches[0] + 7 * branches[1] + 5 * branches[2]) as f64;
        assert_eq!(report.decay.released.deposited_j, heat);
        assert_eq!(report.decay.released.radiation_j, 0.5 * heat);
        assert_eq!(report.decay.released.escaped_j, 0.25 * heat);
        for ((sum, square), &count) in sums.iter_mut().zip(&mut squares).zip(counts) {
            *sum += count as f64;
            *square += (count * count) as f64;
        }
    }
    for ((sum, squares), mean) in sums.into_iter().zip(squares).zip(expected) {
        poisson_moments(sum, squares, samples, mean);
    }
}

#[test]
fn g1_shutdown_preserves_decay_and_binding_at_later_time_does_not_replay_supply() {
    let model = single();
    let schedule = schedule(&model, vec![8.0, 0.0], vec![(1.0, vec![0.0, 0.0])]);
    let samples = 4096;
    let mut surviving = 0.0;
    for seed in 0..samples {
        let mut state = supplied(&model, vec![0, 0], &schedule, seed + 30_000);
        state.advance_to(1.0, budget(), || Ok(())).unwrap();
        let produced = state.cumulative_produced_nuclei().to_vec();
        let report = state.advance_to(3.0, budget(), || Ok(())).unwrap();
        assert_eq!(report.births, 0);
        assert_eq!(state.cumulative_produced_nuclei(), produced);
        assert_eq!(
            state.inventory.populations().iter().sum::<u64>(),
            produced[0]
        );
        assert_eq!(state.source.clocks, [EventClock::Inactive; 2]);
        surviving += state.inventory.populations()[0] as f64;
    }
    let expected = 8.0 / 0.7 * (1.0 - (-0.7_f64).exp()) * (-1.4_f64).exp();
    close(
        surviving / samples as f64,
        expected,
        6.0 * (expected / samples as f64).sqrt(),
    );

    let resumed_schedule = self::schedule(
        &model,
        vec![8.0, 0.0],
        vec![
            (1.0, vec![0.0, 0.0]),
            (2.0, vec![4.0, 0.0]),
            (3.0, vec![0.0, 0.0]),
        ],
    );
    let inventory =
        StochasticDecayInventory::new_at(Arc::clone(&model), vec![0, 0], key(77), 1.5).unwrap();
    let mut resumed = StochasticDecayProduction::poisson(inventory, resumed_schedule).unwrap();
    assert_eq!(resumed.inventory.rng_checkpoint().index, 0);
    assert_eq!(resumed.next_transition_time_s(), Some(2.0));
    let report = resumed.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(report.births, 0);
    assert_eq!(report.schedule_changes, 1);
    assert_eq!(resumed.inventory.rng_checkpoint().index, 1);
    let report = resumed.advance_to(4.0, budget(), || Ok(())).unwrap();
    assert!(report.births > 0);
    assert_eq!(
        resumed.inventory.populations().iter().sum::<u64>(),
        report.births
    );
}

#[test]
fn g3_partition_clone_and_noop_changes_preserve_complete_future_history() {
    let model = model(
        vec![
            node("parent", 0.8),
            node("daughter", 0.8),
            node("stable", 0.0),
        ],
        vec![
            branch(0, 1, 0.25, 0.1),
            branch(0, 2, 0.75, 0.3),
            branch(1, 2, 1.0, 1e-20),
        ],
    );
    let schedule = schedule(
        &model,
        vec![3.0, 0.5, 0.25],
        vec![
            (0.1, vec![3.0, 0.5, 0.25]),
            (0.3, vec![1.0, 2.0, 0.25]),
            (0.9, vec![0.0; 3]),
            (1.2, vec![3.0, 0.0, 0.0]),
            (2.7, vec![0.0; 3]),
        ],
    );
    let initial = supplied(&model, vec![20, 4, 1], &schedule, 91);
    let mut whole = initial.clone();
    whole.advance_to(8.0, budget(), || Ok(())).unwrap();
    let mut split = initial.clone();
    let idle = split.advance_to(0.0, budget(), || Ok(())).unwrap();
    assert_eq!(split, initial);
    assert_eq!(idle.decay.rng_draws, 0);
    for end in [0.01, 0.1, 0.100_000_01, 0.3, 0.7, 0.9, 1.1] {
        split.advance_to(end, budget(), || Ok(())).unwrap();
    }
    let mut restart = split.clone();
    for end in [1.2, 1.5, 2.7, 4.0, 8.0] {
        assert_eq!(
            split.advance_to(end, budget(), || Ok(())).unwrap(),
            restart.advance_to(end, budget(), || Ok(())).unwrap(),
        );
    }
    assert_eq!(whole, split);
    assert_eq!(whole, restart);

    let constant = self::schedule(&model, vec![3.0, 0.5, 0.25], vec![]);
    let noop = self::schedule(
        &model,
        vec![3.0, 0.5, 0.25],
        vec![(0.5, vec![3.0, 0.5, 0.25])],
    );
    let mut left = supplied(&model, vec![20, 4, 1], &constant, 91);
    let mut right = supplied(&model, vec![20, 4, 1], &noop, 91);
    left.advance_to(2.0, budget(), || Ok(())).unwrap();
    let report = right.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(report.schedule_changes, 1);
    assert_eq!(left.inventory, right.inventory);
    assert_eq!(left.source.clocks, right.source.clocks);
    assert_eq!(left.source.produced_nuclei, right.source.produced_nuclei);
}

#[test]
fn g0_zero_source_preserves_the_original_closed_trajectory_and_report_bits() {
    let model = single();
    let schedule = schedule(
        &model,
        vec![0.0, 0.0],
        vec![(0.25, vec![0.0, 0.0]), (1.0, vec![0.0, 0.0])],
    );
    let mut closed =
        StochasticDecayInventory::new(Arc::clone(&model), vec![20, 2], key(41)).unwrap();
    let mut produced = StochasticDecayProduction::poisson(closed.clone(), schedule).unwrap();
    assert_eq!(produced.inventory, closed);
    for end in [0.0, 0.1, 0.25, 0.7, 1.0, 3.0, 8.0] {
        let expected = closed.advance_to(end, budget().decay, || Ok(())).unwrap();
        let report = produced.advance_to(end, budget(), || Ok(())).unwrap();
        assert_eq!(report.decay, expected);
        assert_eq!(produced.inventory, closed);
        assert_eq!(report.births, 0);
        assert_eq!(produced.cumulative_produced_nuclei(), [0, 0]);
    }
}

#[test]
fn g0_compensated_births_precede_shutdown_and_exact_ties_install_new_rates_first() {
    let start = 1e16_f64;
    let end = start.next_up();
    let model = model(
        vec![
            node("parent", 1e6),
            node("daughter", 1e6),
            node("stable", 0.0),
        ],
        vec![branch(0, 1, 1.0, 2.0), branch(1, 2, 1.0, 3.0)],
    );
    let schedule = schedule(&model, vec![1e-300, 0.0, 0.0], vec![(end, vec![0.0; 3])]);
    let inventory =
        StochasticDecayInventory::new_at(Arc::clone(&model), vec![0; 3], key(13), start).unwrap();
    let mut state = StochasticDecayProduction::poisson(inventory, Arc::clone(&schedule)).unwrap();
    // A known retained sub-ULP birth isolates boundary ordering from random draws.
    state.source.clocks[0] = EventClock::Finite(Accumulator { hi: end, lo: -0.25 });
    assert_eq!(state.next_transition_time_s(), Some(end));
    let original = state.clone();
    state.advance_to(start, budget(), || Ok(())).unwrap();
    assert_eq!(state, original);
    let report = state.advance_to(end, budget(), || Ok(())).unwrap();
    assert_eq!(
        (report.births, report.decay.events, report.schedule_changes),
        (1, 2, 1)
    );
    assert_eq!(state.inventory.populations(), [0, 0, 1]);
    assert_eq!(report.decay.released.deposited_j, 5.0);
    let accepted = state.clone();
    assert_eq!(
        state.advance_to(end, budget(), || Ok(())).unwrap().births,
        0
    );
    assert_eq!(state, accepted);

    let inventory = StochasticDecayInventory::new_at(model, vec![1, 0, 0], key(13), start).unwrap();
    let mut tied = StochasticDecayProduction::poisson(inventory, schedule).unwrap();
    tied.source.clocks[0] = EventClock::Finite(Accumulator::new(end));
    tied.inventory.clocks[0] = EventClock::Finite(Accumulator::new(end));
    let report = tied.advance_to(end, budget(), || Ok(())).unwrap();
    assert_eq!(
        (report.births, report.decay.events, report.schedule_changes),
        (0, 1, 1)
    );
    assert_eq!(tied.inventory.populations(), [0, 1, 0]);
    assert_eq!(tied.next_transition_time_s(), Some(end.next_up()));
}

#[test]
fn g0_beyond_finite_arrivals_remain_pending_until_an_authored_rate_change() {
    let model = single();
    let seed = (0..100)
        .find(|&seed| {
            Stream::resume(StreamCheckpoint::current(key(seed), 0))
                .unwrap()
                .next_exponential()
                > 2.0
        })
        .unwrap();
    let constant = schedule(&model, vec![1e-308, 0.0], vec![]);
    let mut state = supplied(&model, vec![0, 0], &constant, seed);
    assert_eq!(state.source.clocks[0], EventClock::BeyondFiniteTime);
    assert!(state.next_transition_time_s().is_none());
    let report = state.advance_to(f64::MAX, budget(), || Ok(())).unwrap();
    assert_eq!(report.births, 0);
    assert_eq!(report.decay.rng_draws, 0);
    let changed = schedule(&model, vec![1e-308, 0.0], vec![(1.0, vec![1.0, 0.0])]);
    let mut state = supplied(&model, vec![0, 0], &changed, seed);
    let report = state.advance_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!(
        (
            report.births,
            report.schedule_changes,
            report.decay.rng_draws
        ),
        (0, 1, 1)
    );
    assert!(
        state
            .next_transition_time_s()
            .is_some_and(|time| time > 1.0)
    );
}

#[test]
fn g4_work_population_counter_rng_and_callback_refusals_are_atomic() {
    let model = single();
    let schedule = schedule(&model, vec![4.0, 0.0], vec![(0.5, vec![2.0, 0.0])]);
    let state = supplied(&model, vec![10, 0], &schedule, 21);
    for (limit, expected) in [
        (
            StochasticDecayProductionBudget {
                maximum_births: 0,
                ..budget()
            },
            "stochastic production births",
        ),
        (
            StochasticDecayProductionBudget {
                maximum_schedule_changes: 0,
                ..budget()
            },
            "stochastic production schedule changes",
        ),
        (
            StochasticDecayProductionBudget {
                decay: StochasticDecayBudget {
                    maximum_events: 0,
                    ..budget().decay
                },
                ..budget()
            },
            "stochastic decay events",
        ),
    ] {
        let mut trial = state.clone();
        assert_eq!(
            trial.advance_to(3.0, limit, || Ok(())),
            Err(DecayError::Budget(expected))
        );
        assert_eq!(trial, state);
    }
    let mut checkpoints = 0;
    state
        .trial_to(3.0, budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, 3, checkpoints] {
        let mut trial = state.clone();
        let mut calls = 0;
        assert_eq!(
            trial.advance_to(3.0, budget(), || {
                calls += 1;
                if calls == stop {
                    Err(DecayError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(DecayError::Cancelled)
        );
        assert_eq!(trial, state);
    }
    let mut capped = supplied(&model, vec![0, MAX_STOCHASTIC_NUCLEI], &schedule, 21);
    let original = capped.clone();
    assert_eq!(
        capped.advance_to(3.0, budget(), || Ok(())),
        Err(DecayError::Budget("stochastic population ceiling"))
    );
    assert_eq!(capped, original);
    for counter in [false, true] {
        let mut trial = supplied(&model, vec![0, 0], &schedule, 21);
        if counter {
            trial.source.produced_nuclei[0] = u64::MAX;
        } else {
            trial.inventory.rng.index = u64::MAX;
        }
        let original = trial.clone();
        let expected = if counter {
            DecayError::Arithmetic
        } else {
            DecayError::Budget("stochastic random stream exhausted")
        };
        assert_eq!(trial.advance_to(3.0, budget(), || Ok(())), Err(expected));
        assert_eq!(trial, original);
    }
}

#[test]
fn g4_model_endpoint_and_budget_admission_precede_publication() {
    let model = single();
    let schedule = schedule(&model, vec![2.0, 0.0], vec![]);
    let wrong_model = StochasticDecayInventory::new(single(), vec![0, 0], key(1)).unwrap();
    assert!(matches!(
        StochasticDecayProduction::poisson(wrong_model, Arc::clone(&schedule)),
        Err(DecayError::Invalid(_))
    ));
    let mut state = supplied(&model, vec![0, 0], &schedule, 1);
    state.advance_to(1.0, budget(), || Ok(())).unwrap();
    let original = state.clone();
    for end in [-1.0, 0.5, f64::NAN, f64::INFINITY] {
        assert!(state.advance_to(end, budget(), || Ok(())).is_err());
        assert_eq!(state, original);
    }
    for invalid in [
        StochasticDecayProductionBudget {
            maximum_births: MAX_STOCHASTIC_NUCLEI + 1,
            ..budget()
        },
        StochasticDecayProductionBudget {
            decay: StochasticDecayBudget {
                energy_tolerance_j: -1.0,
                ..budget().decay
            },
            ..budget()
        },
    ] {
        assert!(state.advance_to(1.0, invalid, || Ok(())).is_err());
        assert_eq!(state, original);
    }
}
