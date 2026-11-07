//! Manufactured integer activation, competing half-lives and accepted histories.
//! Independent exponential/Erlang probabilities supply the physical references;
//! none of these coefficients is evaluated data for the synthetic identities.

use std::sync::Arc;

use super::{
    StochasticActivationBudget, StochasticActivationChange, StochasticActivationSchedule,
    StochasticActivationState,
};
use crate::burnup::{
    BurnupData, BurnupError, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate,
    BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::decay::stochastic::{MAX_STOCHASTIC_NUCLEI, StreamCheckpoint, StreamKey};
use crate::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use crate::{EnergyChannels, FissionError};
use fs_qty::Frequency;
use fs_rand::Stream;
use fs_thermochem::decay::DecayLawV1;

fn limits() -> BurnupLimits {
    BurnupLimits {
        max_neutron_groups: 8,
        max_reactions: 32,
        max_products: 64,
        ..BurnupLimits::default()
    }
}

fn budget() -> StochasticActivationBudget {
    StochasticActivationBudget {
        maximum_events: 100_000,
        maximum_clock_checks: 2_000_000,
        energy_tolerance_j: 1e-8,
    }
}

fn key(seed: u64) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_ac71,
        tile: 0,
    }
}

fn identity(index: u16) -> NuclideId {
    NuclideId::new(10, 25 + index, 0).unwrap()
}

fn energy(deposited_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j: deposited_j * 0.5,
        escaped_j: deposited_j * 0.25,
    }
}

fn capture(parent: usize, daughter: usize, rates: Vec<f64>, heat: f64) -> BurnupReaction {
    BurnupReaction {
        data_id: format!("manufactured-capture-{parent}-{daughter}"),
        parent,
        rate: BurnupRate::Capture {
            cross_sections_m2: rates,
        },
        products: vec![BurnupProduct {
            nuclide: daughter,
            yield_per_event: 1.0,
        }],
        energy_per_event: energy(heat),
    }
}

fn decay(parent: usize, daughter: usize, heat: f64) -> BurnupReaction {
    BurnupReaction {
        data_id: format!("manufactured-decay-{parent}-{daughter}"),
        parent,
        rate: BurnupRate::Decay {
            branch_fraction: 1.0,
        },
        products: vec![BurnupProduct {
            nuclide: daughter,
            yield_per_event: 1.0,
        }],
        energy_per_event: energy(heat),
    }
}

fn model(rates: &[f64], reactions: Vec<BurnupReaction>, groups: usize) -> Arc<BurnupModel> {
    let grid = Arc::new(
        NeutronGroups::new(
            "manufactured-activation-grid".into(),
            (0..=groups).map(|index| index as f64).collect(),
        )
        .unwrap(),
    );
    assemble(rates, reactions, grid)
}

fn assemble(
    rates: &[f64],
    reactions: Vec<BurnupReaction>,
    groups: Arc<NeutronGroups>,
) -> Arc<BurnupModel> {
    Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-discrete-activation-v1".into(),
                groups,
                temperature_min_k: 300.0,
                temperature_max_k: 900.0,
                nuclides: rates
                    .iter()
                    .enumerate()
                    .map(|(index, &rate)| BurnupNuclide {
                        id: identity(index as u16),
                        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
                    })
                    .collect(),
                reactions,
            },
            limits(),
        )
        .unwrap(),
    )
}

fn flux(model: &BurnupModel, values: Vec<f64>) -> GroupFlux {
    GroupFlux::new(Arc::clone(model.groups()), values).unwrap()
}

fn schedule(
    model: &Arc<BurnupModel>,
    initial: Vec<f64>,
    changes: Vec<(f64, Vec<f64>, f64)>,
) -> Arc<StochasticActivationSchedule> {
    Arc::new(
        StochasticActivationSchedule::new(
            Arc::clone(model),
            flux(model, initial),
            300.0,
            changes
                .into_iter()
                .map(
                    |(time_s, values, temperature_k)| StochasticActivationChange {
                        time_s,
                        flux: flux(model, values),
                        temperature_k,
                    },
                )
                .collect(),
            limits(),
            32,
            || Ok(()),
        )
        .unwrap(),
    )
}

fn binomial_reference(actual: u64, total: u64, probability: f64) {
    let expected = total as f64 * probability;
    let allowance = 6.0 * (expected * (1.0 - probability)).sqrt() + 1.0;
    assert!(
        (actual as f64 - expected).abs() < allowance,
        "{actual} versus {expected}, allowance {allowance}"
    );
}

fn same_physics(left: &StochasticActivationState, right: &StochasticActivationState) {
    assert!(Arc::ptr_eq(left.model(), right.model()));
    assert_eq!(left.time_s(), right.time_s());
    assert_eq!(left.populations(), right.populations());
    assert_eq!(
        left.cumulative_reaction_counts(),
        right.cumulative_reaction_counts()
    );
    assert_eq!(left.cumulative_released(), right.cumulative_released());
    assert_eq!(
        left.cumulative_capture_released(),
        right.cumulative_capture_released()
    );
    assert_eq!(
        left.cumulative_decay_released(),
        right.cumulative_decay_released()
    );
    assert_eq!(left.rng_checkpoint(), right.rng_checkpoint());
    assert_eq!(
        left.next_transition_time_s(),
        right.next_transition_time_s()
    );
}

#[test]
fn g1_capture_and_radioactive_removal_compete_without_relabeling_events_or_energy() {
    let lambda = core::f64::consts::LN_2;
    let model = model(
        &[lambda, 0.0, 0.0],
        vec![capture(0, 1, vec![lambda], 2.0), decay(0, 2, 5.0)],
        1,
    );
    let total = 8192;
    let mut state = StochasticActivationState::new(
        schedule(&model, vec![1.0], vec![]),
        vec![total, 0, 0],
        key(101),
    )
    .unwrap();
    let step = state.advance_to(1.0, budget(), || Ok(())).unwrap();
    // At one radioactive half-life, equal capture and decay hazards leave 1/4
    // of the parents and assign each distinct removal cause probability 3/8.
    for (&actual, probability) in state.populations().iter().zip([0.25, 0.375, 0.375]) {
        binomial_reference(actual, total, probability);
    }
    assert_eq!(state.populations().iter().sum::<u64>(), total);
    assert_eq!(step.reaction_counts, [step.captures, step.decays]);
    assert_eq!(step.events, step.captures + step.decays);
    assert_eq!(step.captures, state.populations()[1]);
    assert_eq!(step.decays, state.populations()[2]);
    assert_eq!(step.capture_released, energy(2.0 * step.captures as f64));
    assert_eq!(step.decay_released, energy(5.0 * step.decays as f64));
    assert_eq!(
        step.released,
        energy(2.0 * step.captures as f64 + 5.0 * step.decays as f64)
    );
    assert_eq!(state.cumulative_released(), step.released);
    assert_eq!(state.cumulative_capture_released(), step.capture_released);
    assert_eq!(state.cumulative_decay_released(), step.decay_released);
}

#[test]
fn g1_equal_capture_and_daughter_decay_rates_keep_the_actual_product_birth_time() {
    let lambda = core::f64::consts::LN_2;
    let model = model(
        &[0.0, lambda, 0.0],
        vec![capture(0, 1, vec![lambda], 0.5), decay(1, 2, 1.0)],
        1,
    );
    let total = 8192;
    let mut state = StochasticActivationState::new(
        schedule(&model, vec![1.0], vec![]),
        vec![total, 0, 0],
        key(211),
    )
    .unwrap();
    let step = state.advance_to(1.0, budget(), || Ok(())).unwrap();
    // Independent repeated-root solution: A=e^-kt, B=kt*e^-kt,
    // C=1-(1+kt)*e^-kt. Treating every newborn as present at t=0 fails this.
    for (&actual, probability) in
        state
            .populations()
            .iter()
            .zip([0.5, 0.5 * lambda, 0.5 * (1.0 - lambda)])
    {
        binomial_reference(actual, total, probability);
    }
    assert_eq!(state.populations().iter().sum::<u64>(), total);
    assert_eq!(step.captures, total - state.populations()[0]);
    assert_eq!(step.decays, state.populations()[2]);
    assert_eq!(step.captures - step.decays, state.populations()[1]);
    assert!(
        step.decays > 0,
        "new products must decay in their birth interval"
    );
    assert_eq!(step.capture_released, energy(0.5 * step.captures as f64));
    assert_eq!(step.decay_released, energy(step.decays as f64));
}

#[test]
fn g1_tabulated_multigroup_forcing_changes_rates_and_shutdown_keeps_the_decay_tail() {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-temperature-grid".into(), vec![0.0, 1.0, 3.0]).unwrap(),
    );
    let table = Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: "manufactured-temperature-capture".into(),
                nuclide: identity(0),
                reaction: NeutronReaction::Capture,
                groups: Arc::clone(&groups),
                temperatures_k: vec![300.0, 900.0],
                cross_sections_m2: vec![vec![1.0, 4.0], vec![3.0, 0.0]],
            },
            CrossSectionLimits::default(),
        )
        .unwrap(),
    );
    let lambda = core::f64::consts::LN_2;
    let mut neutron_capture = capture(0, 1, vec![], 2.0);
    neutron_capture.rate = BurnupRate::TabulatedCapture {
        cross_sections: table,
    };
    let model = assemble(
        &[0.0, lambda, 0.0],
        vec![neutron_capture, decay(1, 2, 3.0)],
        groups,
    );
    let calendar = schedule(
        &model,
        vec![0.25, 0.125],
        vec![(1.0, vec![0.125, 0.5], 600.0), (2.0, vec![0.0, 0.0], 900.0)],
    );
    assert_eq!(calendar.rate_constants_at(0.0).unwrap(), &[0.75, lambda]);
    assert_eq!(calendar.rate_constants_at(1.0).unwrap(), &[1.25, lambda]);
    assert_eq!(calendar.rate_constants_at(2.0).unwrap(), &[0.0, lambda]);
    assert_eq!(calendar.temperature_at(1.0).unwrap(), 600.0);
    assert_eq!(calendar.next_after(1.0), Some(2.0));
    let mut state = StochasticActivationState::new(calendar, vec![1024, 0, 0], key(307)).unwrap();
    let irradiated = state.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(irradiated.schedule_changes, 2);
    binomial_reference(state.populations()[0], 1024, (-2.0_f64).exp());
    let shutdown = state.clone();
    let tail = state.advance_to(4.0, budget(), || Ok(())).unwrap();
    assert_eq!(tail.captures, 0);
    assert_eq!(tail.capture_released, EnergyChannels::default());
    assert_eq!(state.populations()[0], shutdown.populations()[0]);
    assert_eq!(state.cumulative_reaction_counts()[0], irradiated.captures);
    assert!(tail.decays > 0);
    assert_eq!(
        tail.decays,
        shutdown.populations()[1] - state.populations()[1]
    );
    binomial_reference(state.populations()[1], shutdown.populations()[1], 0.25);
    assert_eq!(tail.released, tail.decay_released);
    assert_eq!(state.populations().iter().sum::<u64>(), 1024);
}

#[test]
fn g3_noop_conditions_late_rebinding_and_observation_partitions_preserve_future_events() {
    let model = model(
        &[0.0, 0.7, 0.0],
        vec![capture(0, 1, vec![0.4], 0.5), decay(1, 2, 1.0)],
        1,
    );
    let original = StochasticActivationState::new(
        schedule(&model, vec![1.0], vec![]),
        vec![100, 5, 0],
        key(401),
    )
    .unwrap();
    let noop_calendar = schedule(
        &model,
        vec![1.0],
        vec![
            (0.25, vec![1.0], 450.0),
            (0.5, vec![1.0], 600.0),
            (1.0, vec![1.0], 300.0),
        ],
    );
    let noop = original.clone().with_schedule(noop_calendar).unwrap();
    assert_eq!(noop.rng_checkpoint(), original.rng_checkpoint());
    let (plain, plain_step) = original.trial_to(2.5, budget(), || Ok(())).unwrap();
    let (unchanged, noop_step) = noop.trial_to(2.5, budget(), || Ok(())).unwrap();
    same_physics(&plain, &unchanged);
    assert_eq!(plain_step.rng_draws, noop_step.rng_draws);
    assert_eq!(noop_step.schedule_changes, 3);

    let calendar = schedule(
        &model,
        vec![1.0],
        vec![
            (1.0, vec![0.0], 300.0),
            (2.0, vec![2.0], 600.0),
            (3.0, vec![0.0], 900.0),
        ],
    );
    let pulsed = original.with_schedule(Arc::clone(&calendar)).unwrap();
    let (whole, report) = pulsed.trial_to(4.0, budget(), || Ok(())).unwrap();
    let mut split = pulsed.clone();
    let (mut changes, mut captures, mut decays) = (0, 0, 0);
    for end in [0.1, 0.7, 1.0, 1.25, 2.0, 2.25, 3.0, 3.75, 4.0] {
        let step = split.advance_to(end, budget(), || Ok(())).unwrap();
        changes += step.schedule_changes;
        captures += step.captures;
        decays += step.decays;
        if end == 1.25 {
            let (restarted, _) = split.trial_to(4.0, budget(), || Ok(())).unwrap();
            assert_eq!(restarted, whole);
            let replacement = schedule(
                &model,
                vec![5.0],
                vec![
                    (0.5, vec![0.0], 900.0),
                    (2.0, vec![2.0], 600.0),
                    (3.0, vec![0.0], 900.0),
                ],
            );
            let rebound = split.clone().with_schedule(replacement).unwrap();
            same_physics(&split, &rebound);
            let (continued, _) = rebound.trial_to(4.0, budget(), || Ok(())).unwrap();
            same_physics(&whole, &continued);
        }
    }
    assert_eq!(split, whole);
    assert_eq!(
        (changes, captures, decays),
        (report.schedule_changes, report.captures, report.decays)
    );
    assert_eq!(changes, 3);
    assert!(Arc::ptr_eq(split.schedule(), &calendar));
}

#[test]
fn g3_exact_shutdown_ties_and_sub_ulp_events_keep_their_true_order_and_rng_draws() {
    let stream_key = key(503);
    let mut stream = Stream::resume(StreamCheckpoint::current(stream_key, 0)).unwrap();
    let rate = stream.next_exponential();
    assert!(rate > 0.0);
    let tied_model = model(&[0.0, 0.0], vec![capture(0, 1, vec![rate], 0.0)], 1);
    let calendar = schedule(
        &tied_model,
        vec![1.0],
        vec![
            (0.5, vec![1.0], 600.0),
            (1.0, vec![0.0], 900.0),
            (2.0, vec![1.0], 300.0),
        ],
    );
    let mut state = StochasticActivationState::new(calendar, vec![1, 0], stream_key).unwrap();
    // k equals the first exponential draw, so the pending capture is exactly at 1 s.
    assert_eq!(state.rng_checkpoint(), stream.checkpoint());
    let noop = state.advance_to(0.5, budget(), || Ok(())).unwrap();
    assert_eq!((noop.schedule_changes, noop.events), (1, 0));
    assert_eq!(state.rng_checkpoint(), stream.checkpoint());
    let closed = state.advance_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!((closed.schedule_changes, closed.events), (1, 0));
    assert_eq!(state.populations(), &[1, 0]);
    assert_eq!(state.rng_checkpoint(), stream.checkpoint());
    let next_wait = stream.next_exponential() / rate;
    let opened = state.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!((opened.schedule_changes, opened.events), (1, 0));
    assert_eq!(state.rng_checkpoint(), stream.checkpoint());
    let next = state.next_transition_time_s().unwrap();
    let rounded = 2.0 + next_wait;
    assert!(next > 2.0 && next >= rounded && next <= rounded.next_up());
    assert_eq!(
        state
            .advance_to(next, budget(), || Ok(()))
            .unwrap()
            .captures,
        1
    );
    assert_eq!(state.populations(), &[0, 1]);

    let start = 2.0_f64.powi(60);
    let closure = start.next_up();
    let short_model = model(&[0.0, 0.0], vec![capture(0, 1, vec![1.0], 0.0)], 1);
    let calendar = schedule(&short_model, vec![1.0], vec![(closure, vec![0.0], 300.0)]);
    let mut short =
        StochasticActivationState::new_at(calendar, vec![1, 0], key(509), start).unwrap();
    assert_eq!(short.next_transition_time_s(), Some(closure));
    let report = short.advance_to(closure, budget(), || Ok(())).unwrap();
    assert_eq!((report.captures, report.schedule_changes), (1, 1));
    assert_eq!(short.populations(), &[0, 1]);
    assert!(short.next_transition_time_s().is_none());
}

#[test]
fn g4_noninteger_physics_admission_and_trial_refusals_publish_no_partial_history() {
    for kind in 0..5 {
        let mut reaction = capture(0, 1, vec![1.0], 1.0);
        match kind {
            0 => {
                reaction.rate = BurnupRate::Fission {
                    cross_sections_m2: vec![0.0],
                }
            }
            1 => reaction.products[0].yield_per_event = 0.5,
            2 => reaction.products.push(BurnupProduct {
                nuclide: 2,
                yield_per_event: 1.0,
            }),
            3 => reaction.products.clear(),
            _ => reaction.products[0].yield_per_event = 2.0,
        }
        let unsupported = model(&[0.0; 3], vec![reaction], 1);
        assert!(
            StochasticActivationSchedule::new(
                Arc::clone(&unsupported),
                flux(&unsupported, vec![1.0]),
                300.0,
                vec![],
                limits(),
                0,
                || Ok(())
            )
            .is_err()
        );
    }
    // The expectation owner admits small branch closure residuals. Integer
    // events must refuse this incomplete authored distribution without repair.
    let mut incomplete = decay(0, 1, 1.0);
    incomplete.rate = BurnupRate::Decay {
        branch_fraction: 1.0 - f64::EPSILON,
    };
    let incomplete = model(&[1.0, 0.0], vec![incomplete], 1);
    assert!(
        StochasticActivationSchedule::new(
            Arc::clone(&incomplete),
            flux(&incomplete, vec![0.0]),
            300.0,
            vec![],
            limits(),
            0,
            || Ok(())
        )
        .is_err()
    );

    let model = model(
        &[0.0, 0.6, 0.0],
        vec![capture(0, 1, vec![1.0], 0.5), decay(1, 2, 1.0)],
        1,
    );
    let change = |time_s| StochasticActivationChange {
        time_s,
        flux: flux(&model, vec![0.0]),
        temperature_k: 600.0,
    };
    for time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            StochasticActivationSchedule::new(
                Arc::clone(&model),
                flux(&model, vec![1.0]),
                300.0,
                vec![change(time)],
                limits(),
                1,
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(
        StochasticActivationSchedule::new(
            Arc::clone(&model),
            flux(&model, vec![1.0]),
            300.0,
            vec![change(1.0), change(1.0)],
            limits(),
            2,
            || Ok(())
        )
        .is_err()
    );
    assert!(
        StochasticActivationSchedule::new(
            Arc::clone(&model),
            flux(&model, vec![1.0]),
            300.0,
            vec![change(1.0)],
            limits(),
            0,
            || Ok(())
        )
        .is_err()
    );
    assert!(
        StochasticActivationSchedule::new(
            Arc::clone(&model),
            flux(&model, vec![1.0]),
            299.0,
            vec![],
            limits(),
            0,
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(
        StochasticActivationSchedule::new(
            Arc::clone(&model),
            flux(&model, vec![1.0]),
            300.0,
            vec![],
            limits(),
            0,
            || Err(FissionError::Cancelled)
        )
        .unwrap_err(),
        BurnupError::Fission(FissionError::Cancelled)
    );

    let calendar = schedule(&model, vec![1.0], vec![(1.0, vec![0.0], 600.0)]);
    for populations in [
        vec![1, 0],
        vec![MAX_STOCHASTIC_NUCLEI + 1, 0, 0],
        vec![u64::MAX, 1, 0],
    ] {
        assert!(
            StochasticActivationState::new(Arc::clone(&calendar), populations, key(601)).is_err()
        );
    }
    assert!(
        StochasticActivationState::new_at(Arc::clone(&calendar), vec![1, 0, 0], key(601), -1.0)
            .is_err()
    );
    let original =
        StochasticActivationState::new(Arc::clone(&calendar), vec![40, 4, 0], key(607)).unwrap();
    let alien = Arc::new(model.as_ref().clone());
    assert!(
        original
            .clone()
            .with_schedule(schedule(&alien, vec![1.0], vec![]))
            .is_err()
    );
    let mut checks = 0;
    let (accepted, report) = original
        .trial_to(2.0, budget(), || {
            checks += 1;
            Ok(())
        })
        .unwrap();
    assert!(report.events > 0 && report.clock_checks > 0 && checks > 3);
    for stop in [1, checks / 2, checks] {
        let mut attempt = original.clone();
        let mut seen = 0;
        assert_eq!(
            attempt
                .advance_to(2.0, budget(), || {
                    seen += 1;
                    if seen == stop {
                        Err(FissionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .unwrap_err(),
            BurnupError::Fission(FissionError::Cancelled)
        );
        assert_eq!(seen, stop);
        assert_eq!(attempt, original);
    }
    for tight in [
        StochasticActivationBudget {
            maximum_events: report.events - 1,
            ..budget()
        },
        StochasticActivationBudget {
            maximum_clock_checks: report.clock_checks - 1,
            ..budget()
        },
        StochasticActivationBudget {
            energy_tolerance_j: -1.0,
            ..budget()
        },
    ] {
        let mut attempt = original.clone();
        assert!(attempt.advance_to(2.0, tight, || Ok(())).is_err());
        assert_eq!(attempt, original);
    }
    for endpoint in [-1.0, f64::NAN, f64::INFINITY] {
        let mut attempt = original.clone();
        assert!(attempt.advance_to(endpoint, budget(), || Ok(())).is_err());
        assert_eq!(attempt, original);
    }
    let (repeated, _) = original.trial_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(repeated, accepted);
}
