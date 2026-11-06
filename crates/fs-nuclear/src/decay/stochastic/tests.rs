//! G0/G1/G3/G4/G5: integer conservation, analytic moments, event ordering,
//! interval/restart identity and atomic refusal. All coefficients are synthetic.

use super::*;
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::{Frequency, Time};
use fs_thermochem::decay::DecayLawV1;

fn key(seed: u64) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_0001,
        tile: 17,
    }
}

fn budget() -> StochasticDecayBudget {
    StochasticDecayBudget {
        maximum_events: 10_000,
        energy_tolerance_j: 1e-8,
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
            "synthetic-stochastic-half-life".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn single() -> Arc<DecayModel> {
    model(
        vec![
            DecayNuclide {
                id: "parent".into(),
                law: DecayLawV1::from_half_life(Time::new(2.0)).unwrap(),
            },
            node("stable", 0.0),
        ],
        vec![branch(0, 1, 1.0, 2.0)],
    )
}

fn branching_chain() -> Arc<DecayModel> {
    // Original ordering deliberately differs from topological ordering.
    model(
        vec![
            node("second-stable", 0.0),
            node("parent", 1.0),
            node("first-stable", 0.0),
            node("daughter", 1.0),
        ],
        vec![
            branch(1, 3, 0.25, 3.0),
            branch(3, 0, 1.0, 5.0),
            branch(1, 2, 0.75, 7.0),
        ],
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.4e}"
    );
}

#[test]
fn g1_one_atom_has_bernoulli_survival_and_binomial_variance_at_one_half_life() {
    let model = single();
    let samples = 8192_u64;
    let mut survived = 0_u64;
    for seed in 0..samples {
        let mut state =
            StochasticDecayInventory::new(Arc::clone(&model), vec![1, 0], key(seed)).unwrap();
        let report = state.advance_to(2.0, budget(), || Ok(())).unwrap();
        let parent = state.populations()[0];
        assert!(parent <= 1);
        assert_eq!(state.populations()[1], 1 - parent);
        assert_eq!(report.events, 1 - parent);
        assert_eq!(report.decays, [1 - parent, 0]);
        assert_eq!(report.branch_decays, [1 - parent]);
        assert_eq!(report.released.deposited_j, 2.0 * (1 - parent) as f64);
        assert_eq!(report.released.radiation_j, (1 - parent) as f64);
        assert_eq!(report.released.escaped_j, 0.5 * (1 - parent) as f64);
        survived += parent;
    }
    let mean = survived as f64 / samples as f64;
    close(mean, 0.5, 0.025);
    close(mean * (1.0 - mean), 0.25, 0.001);
}

#[test]
fn g1_branching_chain_moments_match_bateman_and_the_shared_expectation_owner() {
    let model = branching_chain();
    let initial = vec![0.0, 8.0, 0.0, 0.0];
    let duration = 1.3;
    let mut expected = DecayInventory::new(Arc::clone(&model), initial).unwrap();
    let expected_step = expected
        .advance(
            duration,
            DecayBudget {
                propagation: DecayLimits::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    let survival = (-duration).exp();
    let bateman = [
        2.0 * (1.0 - (1.0 + duration) * survival),
        8.0 * survival,
        6.0 * (1.0 - survival),
        2.0 * duration * survival,
    ];
    for (&actual, &analytic) in expected.populations().iter().zip(&bateman) {
        close(actual, analytic, 2e-12);
    }
    let mut populations = [0_u64; 4];
    let mut decays = [0_u64; 4];
    let mut branches = [0_u64; 3];
    let mut energy = [0.0; 3];
    let samples = 4096_u64;
    for seed in 0..samples {
        let mut state =
            StochasticDecayInventory::new(Arc::clone(&model), vec![0, 8, 0, 0], key(seed + 10_000))
                .unwrap();
        let report = state.advance_to(duration, budget(), || Ok(())).unwrap();
        assert_eq!(state.populations().iter().sum::<u64>(), 8);
        assert_eq!(
            report.decays[1],
            report.branch_decays[0] + report.branch_decays[2]
        );
        assert_eq!(report.decays[3], report.branch_decays[1]);
        assert_eq!(state.populations()[0], report.branch_decays[1]);
        assert_eq!(
            state.populations()[3],
            report.branch_decays[0] - report.branch_decays[1]
        );
        let deposited = (3 * report.branch_decays[0]
            + 5 * report.branch_decays[1]
            + 7 * report.branch_decays[2]) as f64;
        assert_eq!(report.released.deposited_j, deposited);
        assert_eq!(report.released.radiation_j, 0.5 * deposited);
        assert_eq!(report.released.escaped_j, 0.25 * deposited);
        for ((population, decay), (&p, &d)) in populations
            .iter_mut()
            .zip(&mut decays)
            .zip(state.populations().iter().zip(&report.decays))
        {
            *population += p;
            *decay += d;
        }
        for (total, count) in branches.iter_mut().zip(&report.branch_decays) {
            *total += count;
        }
        for (total, value) in energy.iter_mut().zip(channel_array(report.released)) {
            *total += value;
        }
    }
    for ((&population, &decay), (&p, &d)) in populations
        .iter()
        .zip(&decays)
        .zip(expected.populations().iter().zip(&expected_step.decays))
    {
        close(population as f64 / samples as f64, p, 0.1);
        close(decay as f64 / samples as f64, d, 0.1);
    }
    for (&count, &expected) in branches.iter().zip(&expected_step.branch_decays) {
        close(count as f64 / samples as f64, expected, 0.1);
    }
    for (&value, expected) in energy.iter().zip(channel_array(expected_step.released)) {
        close(value / samples as f64, expected, 0.75);
    }
}

#[test]
fn g3_partition_restart_and_read_only_observation_retain_the_identical_trajectory() {
    let original = branching_chain();
    // Nonexact and disparate energies expose partition-dependent aggregation:
    // lifetime accounts must receive events, not rounded interval exports.
    let shared = model(
        original.nuclides().to_vec(),
        vec![
            branch(1, 3, 0.25, 0.1),
            branch(3, 0, 1.0, 1e-20),
            branch(1, 2, 0.75, 0.3),
        ],
    );
    let mut whole = StochasticDecayInventory::new(shared, vec![0, 80, 0, 15], key(123)).unwrap();
    let mut split = whole.clone();
    let first = whole.next_event_time_s().unwrap();
    let before = whole.clone();
    let (observed, report) = whole.trial_to(first * 0.5, budget(), || Ok(())).unwrap();
    assert_eq!(whole, before);
    assert_eq!(report.events, 0);
    assert_eq!(observed.rng_checkpoint(), whole.rng_checkpoint());
    assert_eq!(observed.next_event_time_s(), whole.next_event_time_s());
    whole.advance_to(8.0, budget(), || Ok(())).unwrap();
    for end in [0.0, 0.01, 0.1, 0.100_000_01, 0.5, 1.0, 1.3] {
        split.advance_to(end, budget(), || Ok(())).unwrap();
    }
    let mut resumed = split.clone();
    for end in [2.7, 4.0, 8.0] {
        split.advance_to(end, budget(), || Ok(())).unwrap();
        resumed.advance_to(end, budget(), || Ok(())).unwrap();
    }
    assert_eq!(split, whole);
    assert_eq!(resumed, whole);
    assert_eq!(split.rng_checkpoint(), whole.rng_checkpoint());
    assert_eq!(split.cumulative_released(), whole.cumulative_released());
}

#[test]
fn g0_endpoint_event_occurs_once_and_small_waits_survive_a_large_absolute_clock() {
    let model = model(
        vec![
            node("parent", 1e6),
            node("daughter", 1e6),
            node("stable", 0.0),
        ],
        vec![branch(0, 1, 1.0, 1.0), branch(1, 2, 1.0, 2.0)],
    );
    let start = 1e16_f64;
    let mut state = StochasticDecayInventory::new_at(model, vec![1, 0, 0], key(31), start).unwrap();
    let first = state.next_event_time_s().unwrap();
    assert_eq!(first, start.next_up());
    let original = state.clone();
    assert_eq!(
        state.advance_to(start, budget(), || Ok(())).unwrap().events,
        0
    );
    assert_eq!(state, original);
    let report = state.advance_to(first, budget(), || Ok(())).unwrap();
    assert_eq!(report.events, 2);
    assert_eq!(state.populations(), [0, 0, 1]);
    assert_eq!(state.cumulative_decays(), [1, 1, 0]);
    assert_eq!(state.cumulative_branch_decays(), [1, 1]);
    assert_eq!(report.released.deposited_j, 3.0);
    assert!(state.next_event_time_s().is_none());
    let original = state.clone();
    assert_eq!(
        state.advance_to(first, budget(), || Ok(())).unwrap().events,
        0
    );
    assert_eq!(state, original);
}

#[test]
fn g0_stable_empty_and_beyond_finite_horizon_populations_do_not_invent_events() {
    for populations in [vec![0, 0], vec![0, 20]] {
        let mut state = StochasticDecayInventory::new(single(), populations, key(4)).unwrap();
        assert_eq!(state.rng_checkpoint().index, 0);
        assert!(state.next_event_time_s().is_none());
        let report = state.advance_to(f64::MAX, budget(), || Ok(())).unwrap();
        assert_eq!(report.events, 0);
        assert_eq!(report.rng_draws, 0);
        assert_eq!(report.released, EnergyChannels::default());
    }
    let mut late =
        StochasticDecayInventory::new_at(single(), vec![1, 0], key(4), f64::MAX).unwrap();
    assert!(late.next_event_time_s().is_none());
    assert_eq!(late.clocks[0], EventClock::BeyondFiniteTime);
    let report = late.advance_to(f64::MAX, budget(), || Ok(())).unwrap();
    assert_eq!(report.events, 0);
    assert_eq!(late.populations(), [1, 0]);
    assert!(late.model().nuclides()[0].law.half_life().is_some());

    let slow = model(
        vec![node("slow", 1e-308), node("stable", 0.0)],
        vec![branch(0, 1, 1.0, 1.0)],
    );
    let seed = (0..100)
        .find(|&seed| {
            Stream::resume(StreamCheckpoint::current(key(seed), 0))
                .unwrap()
                .next_exponential()
                > 2.0
        })
        .unwrap();
    let mut state = StochasticDecayInventory::new(slow, vec![1, 0], key(seed)).unwrap();
    assert!(state.next_event_time_s().is_none());
    assert_eq!(state.clocks[0], EventClock::BeyondFiniteTime);
    assert_eq!(
        state
            .advance_to(f64::MAX, budget(), || Ok(()))
            .unwrap()
            .events,
        0
    );
    assert_eq!(state.populations(), [1, 0]);
    assert_eq!(state.rng_checkpoint().index, 1);
}

#[test]
fn g4_event_budget_callbacks_rng_and_counter_refusals_retain_all_history() {
    let state = StochasticDecayInventory::new(single(), vec![20, 0], key(7)).unwrap();
    let first = state.next_event_time_s().unwrap();
    let mut trial = state.clone();
    assert_eq!(
        trial.advance_to(
            first,
            StochasticDecayBudget {
                maximum_events: 0,
                ..budget()
            },
            || Ok(())
        ),
        Err(DecayError::Budget("stochastic decay events")),
    );
    assert_eq!(trial, state);
    let mut checkpoints = 0;
    state
        .trial_to(100.0, budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(checkpoints > 2);
    for stop in [1, 3, checkpoints] {
        let mut calls = 0;
        let mut trial = state.clone();
        assert_eq!(
            trial.advance_to(100.0, budget(), || {
                calls += 1;
                if calls == stop {
                    Err(DecayError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(DecayError::Cancelled),
        );
        assert_eq!(trial, state);
    }
    let mut exhausted = state.clone();
    exhausted.rng.index = u64::MAX;
    let before = exhausted.clone();
    assert_eq!(
        exhausted.advance_to(first, budget(), || Ok(())),
        Err(DecayError::Budget("stochastic random stream exhausted")),
    );
    assert_eq!(exhausted, before);
    let mut overflow = state.clone();
    overflow.decays[0] = u64::MAX;
    let before = overflow.clone();
    assert_eq!(
        overflow.advance_to(first, budget(), || Ok(())),
        Err(DecayError::Arithmetic)
    );
    assert_eq!(overflow, before);
}

#[test]
fn g4_population_probability_time_and_energy_admission_are_explicit() {
    let shared = single();
    for populations in [
        vec![1],
        vec![MAX_STOCHASTIC_NUCLEI + 1, 0],
        vec![u64::MAX, 1],
    ] {
        assert!(StochasticDecayInventory::new(Arc::clone(&shared), populations, key(1)).is_err());
    }
    for time in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(
            StochasticDecayInventory::new_at(Arc::clone(&shared), vec![1, 0], key(1), time)
                .is_err()
        );
    }
    let imperfect = model_with_incomplete_categorical_closure();
    assert!(matches!(
        StochasticDecayInventory::new(imperfect, vec![1, 0], key(1)),
        Err(DecayError::Invalid(_)),
    ));
    let mut state = StochasticDecayInventory::new(shared, vec![10, 0], key(1)).unwrap();
    let original = state.clone();
    for end in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(state.advance_to(end, budget(), || Ok(())).is_err());
        assert_eq!(state, original);
    }
    for invalid in [
        StochasticDecayBudget {
            maximum_events: MAX_STOCHASTIC_NUCLEI + 1,
            ..budget()
        },
        StochasticDecayBudget {
            energy_tolerance_j: -1.0,
            ..budget()
        },
        StochasticDecayBudget {
            energy_tolerance_j: f64::NAN,
            ..budget()
        },
    ] {
        assert!(state.advance_to(0.0, invalid, || Ok(())).is_err());
        assert_eq!(state, original);
    }

    let huge_energy = model(
        vec![node("parent", 1.0), node("stable", 0.0)],
        vec![DecayBranch {
            energy_per_decay: EnergyChannels {
                deposited_j: 1e308,
                ..EnergyChannels::default()
            },
            ..branch(0, 1, 1.0, 1.0)
        }],
    );
    let mut overflow = StochasticDecayInventory::new(huge_energy, vec![2, 0], key(1)).unwrap();
    let original = overflow.clone();
    assert_eq!(
        overflow.advance_to(1e300, budget(), || Ok(())),
        Err(DecayError::Arithmetic)
    );
    assert_eq!(overflow, original);

    let tiny = f64::EPSILON * 0.25;
    let mixed = model(
        vec![node("large", 1.0), node("small", 1.0), node("stable", 0.0)],
        vec![branch(0, 2, 1.0, 1.0), branch(1, 2, 1.0, tiny)],
    );
    let mut state = StochasticDecayInventory::new(mixed, vec![1, 1, 0], key(12)).unwrap();
    let original = state.clone();
    assert_eq!(
        state.advance_to(
            100.0,
            StochasticDecayBudget {
                energy_tolerance_j: 0.0,
                ..budget()
            },
            || Ok(()),
        ),
        Err(DecayError::Budget("stochastic emission energy roundoff")),
    );
    assert_eq!(state, original);
    let report = state.advance_to(100.0, budget(), || Ok(())).unwrap();
    assert_eq!(report.events, 2);
    assert!(report.energy_roundoff_j >= tiny);
    assert_eq!(state.released[0].hi, 1.0);
    assert_eq!(state.released[0].lo, tiny);
}

fn model_with_incomplete_categorical_closure() -> Arc<DecayModel> {
    model(
        vec![node("parent", 1.0), node("stable", 0.0)],
        vec![branch(0, 1, 1.0 - f64::EPSILON, 1.0)],
    )
}
