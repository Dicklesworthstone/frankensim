//! G0/G1/G3/G4: authored flow changes, competing physical clocks and replay.
//! Synthetic coefficients; analytic probabilities and existing owner histories
//! provide the references, with no replacement event solver.

use super::StochasticDecayTransportSchedule;
use crate::EnergyChannels;
use crate::burnup::transport::schedule::NuclideTransportChange;
use crate::decay::production::{DecayProductionChange, DecayProductionSchedule};
use crate::decay::stochastic::transport::production::{
    StochasticDecayProductionTransportBudget, StochasticDecayProductionTransportState,
};
use crate::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use crate::decay::{DecayBranch, DecayError, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_rand::Stream;
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn model(rates: &[f64], energies: &[f64]) -> Arc<DecayModel> {
    assert_eq!(rates.len(), energies.len() + 1);
    Arc::new(
        DecayModel::new(
            "synthetic-scheduled-movement".into(),
            rates
                .iter()
                .enumerate()
                .map(|(i, &rate)| DecayNuclide {
                    id: format!("nuclide-{i}"),
                    law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
                })
                .collect(),
            energies
                .iter()
                .enumerate()
                .map(|(parent, &heat)| DecayBranch {
                    parent,
                    daughter: parent + 1,
                    probability: 1.0,
                    energy_per_decay: EnergyChannels {
                        deposited_j: heat,
                        radiation_j: 0.5 * heat,
                        escaped_j: 0.25 * heat,
                    },
                })
                .collect(),
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn paths(forward: &[f64], backward: &[f64], outflow: &[f64]) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![
            NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: forward.to_vec(),
            },
            NuclideTransfer {
                from_compartment: 1,
                to_compartment: 0,
                rates_per_s: backward.to_vec(),
            },
        ],
        outflows: vec![NuclideOutflow {
            compartment: 1,
            rates_per_s: outflow.to_vec(),
        }],
    }
}

fn extraction(rate: f64) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![],
        outflows: vec![NuclideOutflow {
            compartment: 0,
            rates_per_s: vec![rate, 0.0],
        }],
    }
}

fn cohort(
    model: &Arc<DecayModel>,
    populations: Vec<Vec<u64>>,
    transport: NuclideTransport,
    seed: u64,
    time: f64,
) -> StochasticDecayTransportState {
    StochasticDecayTransportState::new(
        populations
            .into_iter()
            .enumerate()
            .map(|(cell, counts)| {
                StochasticDecayInventory::new_at(
                    Arc::clone(model),
                    counts,
                    StreamKey {
                        seed,
                        kernel: 0xdec0_7701,
                        tile: cell as u32,
                    },
                    time,
                )
                .unwrap()
            })
            .collect(),
        transport,
        StreamKey {
            seed,
            kernel: 0xdec0_7702,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap()
}

fn schedule(
    state: &StochasticDecayTransportState,
    initial: NuclideTransport,
    changes: Vec<(f64, NuclideTransport)>,
) -> Arc<StochasticDecayTransportSchedule> {
    Arc::new(
        StochasticDecayTransportSchedule::new(
            state,
            initial,
            changes
                .into_iter()
                .map(|(time_s, transport)| NuclideTransportChange { time_s, transport })
                .collect(),
            64,
        )
        .unwrap(),
    )
}

fn source(
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

fn budget() -> StochasticDecayTransportBudget {
    StochasticDecayTransportBudget {
        decay: StochasticDecayBudget {
            maximum_events: 10_000,
            energy_tolerance_j: 1e-8,
        },
        maximum_movements: 10_000,
        maximum_clock_checks: 1_000_000,
    }
}

fn produced_budget() -> StochasticDecayProductionTransportBudget {
    StochasticDecayProductionTransportBudget {
        transport: budget(),
        maximum_births: 10_000,
        maximum_schedule_changes: 64,
    }
}

fn same_physics(a: &StochasticDecayTransportState, b: &StochasticDecayTransportState) {
    // Inventory equality also compares all retained nuclear clocks/history.
    assert_eq!(a.inventories(), b.inventories());
    assert_eq!(a.transport(), b.transport());
    assert_eq!(a.movement_rng_checkpoint(), b.movement_rng_checkpoint());
    assert_eq!(
        a.cumulative_transfer_counts(),
        b.cumulative_transfer_counts()
    );
    assert_eq!(a.cumulative_outflow_counts(), b.cumulative_outflow_counts());
    assert_eq!(a.next_event_time_s(), b.next_event_time_s());
}

#[test]
fn g1_scheduled_extraction_competes_with_decay_at_the_authored_time_not_an_average_rate() {
    let lambda = core::f64::consts::LN_2;
    let shared = model(&[lambda, 0.0], &[2.0]);
    let (mut survived, mut decayed, mut removed) = (0_u64, 0_u64, 0_u64);
    let samples = 4096_u64;
    for seed in 0..samples {
        let bare = cohort(&shared, vec![vec![1, 0]], extraction(0.0), seed, 0.0);
        let calendar = schedule(&bare, extraction(0.0), vec![(1.0, extraction(lambda))]);
        let mut candidate = bare.with_transport_schedule(calendar).unwrap();
        let step = candidate.advance_to(2.0, budget(), || Ok(())).unwrap();
        let parent = candidate.inventories()[0].populations()[0];
        let daughter = candidate.inventories()[0].populations()[1];
        let outflow = step.outflow_counts[0][0];
        assert_eq!(parent + daughter + outflow, 1);
        assert_eq!(step.decay_events, daughter);
        assert_eq!(step.movement_events, outflow);
        assert_eq!(step.transport_changes, 1);
        assert_eq!(step.decay[0].released.deposited_j, 2.0 * daughter as f64);
        survived += parent;
        decayed += daughter;
        removed += outflow;
    }
    // First half-life: survival1/2 and decay1/2. After opening k=lambda,
    // the remaining half has survival1/4 and equal decay/extraction shares3/8.
    for (actual, probability) in [(survived, 0.125), (decayed, 0.6875), (removed, 0.1875)] {
        let expected = samples as f64 * probability;
        let allowance = 6.0 * (expected * (1.0 - probability)).sqrt();
        assert!(
            (actual as f64 - expected).abs() < allowance,
            "{actual} versus {expected}"
        );
    }
}

#[test]
fn g1_valves_close_reopen_and_reverse_circulation_without_moving_daughter_heat() {
    let shared = model(&[100.0, 100.0, 0.0], &[3.0, 5.0]);
    let initial = paths(&[0.0, 1e8, 0.0], &[0.0; 3], &[0.0; 3]);
    let closed = paths(&[0.0; 3], &[0.0; 3], &[0.0; 3]);
    let reverse = paths(&[0.0; 3], &[0.0, 0.0, 1e6], &[0.0; 3]);
    let forward = paths(&[0.0, 0.0, 1e6], &[0.0; 3], &[0.0; 3]);
    let bare = cohort(
        &shared,
        vec![vec![1, 0, 0], vec![0; 3]],
        initial.clone(),
        17,
        0.0,
    );
    let calendar = schedule(
        &bare,
        initial,
        vec![
            (1.0, closed.clone()),
            (2.0, reverse.clone()),
            (3.0, closed.clone()),
            (4.0, forward.clone()),
            (5.0, closed.clone()),
        ],
    );
    let mut candidate = bare.with_transport_schedule(Arc::clone(&calendar)).unwrap();
    let first = candidate.advance_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!(first.decay[0].branch_decays, [1, 0]);
    assert_eq!(first.decay[1].branch_decays, [0, 1]);
    assert_eq!(candidate.inventories()[1].populations(), [0, 0, 1]);
    assert_eq!(candidate.transport(), &closed);
    let heat: Vec<_> = candidate
        .inventories()
        .iter()
        .map(|cell| cell.cumulative_released())
        .collect();
    assert_eq!((heat[0].deposited_j, heat[1].deposited_j), (3.0, 5.0));
    for (end, active, location) in [
        (2.0, &reverse, 1),
        (3.0, &closed, 0),
        (4.0, &forward, 0),
        (5.0, &closed, 1),
    ] {
        let step = candidate.advance_to(end, budget(), || Ok(())).unwrap();
        assert_eq!(step.transport_changes, 1);
        assert_eq!(step.decay_events, 0);
        assert_eq!(candidate.transport(), active);
        assert_eq!(candidate.inventories()[location].populations()[2], 1);
        assert_eq!(candidate.inventories()[1 - location].populations()[2], 0);
        assert_eq!(
            candidate
                .inventories()
                .iter()
                .map(|cell| cell.cumulative_released())
                .collect::<Vec<_>>(),
            heat
        );
    }
    assert_eq!(
        candidate.cumulative_transfer_counts(),
        [vec![0, 1, 1], vec![0, 0, 1]]
    );
    assert!(Arc::ptr_eq(
        candidate.transport_schedule().unwrap(),
        &calendar
    ));
    assert!(candidate.next_event_time_s().is_none());
}

#[test]
fn g3_constant_noop_and_late_bound_frames_preserve_existing_clocks_and_future_history() {
    let shared = model(&[0.7, 0.0], &[0.1]);
    let frame = paths(&[0.4, 1.0], &[0.3, 0.5], &[0.1, 0.2]);
    let mut original = cohort(
        &shared,
        vec![vec![20, 5], vec![15, 3]],
        frame.clone(),
        43,
        0.0,
    );
    let constant = schedule(&original, frame.clone(), vec![]);
    let noop = schedule(
        &original,
        frame.clone(),
        vec![
            (0.25, frame.clone()),
            (0.5, frame.clone()),
            (1.0, frame.clone()),
        ],
    );
    let mut constant_state = original.clone().with_transport_schedule(constant).unwrap();
    same_physics(&constant_state, &original);
    let mut noop_state = original.clone().with_transport_schedule(noop).unwrap();
    assert_eq!(noop_state.inventories(), original.inventories());
    assert_eq!(
        noop_state.movement_rng_checkpoint(),
        original.movement_rng_checkpoint()
    );
    let plain = original.advance_to(2.0, budget(), || Ok(())).unwrap();
    let fixed = constant_state.advance_to(2.0, budget(), || Ok(())).unwrap();
    let unchanged = noop_state.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(fixed.transport_changes, 0);
    assert_eq!(unchanged.transport_changes, 3);
    assert_eq!(plain.movement_events, unchanged.movement_events);
    same_physics(&constant_state, &original);
    same_physics(&noop_state, &original);
    let past = schedule(
        &original,
        paths(&[0.0; 2], &[0.0; 2], &[0.0; 2]),
        vec![
            (0.5, paths(&[7.0; 2], &[0.0; 2], &[0.0; 2])),
            (1.0, frame.clone()),
        ],
    );
    let mut late = original.clone().with_transport_schedule(past).unwrap();
    same_physics(&late, &original);
    late.advance_to(4.0, budget(), || Ok(())).unwrap();
    original.advance_to(4.0, budget(), || Ok(())).unwrap();
    same_physics(&late, &original);

    // A changed active rate first starts its clock at the actual late binding
    // time, exactly like a fresh current inventory with that rate, not at0.5.
    let stable = model(&[0.0], &[]);
    let closed = paths(&[0.0], &[0.0], &[0.0]);
    let open = paths(&[1.0], &[0.0], &[0.0]);
    let bare = cohort(&stable, vec![vec![1], vec![0]], closed.clone(), 59, 2.0);
    let calendar = schedule(&bare, closed, vec![(0.5, open.clone())]);
    let rebound = bare.with_transport_schedule(calendar).unwrap();
    let expected = cohort(&stable, vec![vec![1], vec![0]], open, 59, 2.0);
    same_physics(&rebound, &expected);
    assert!(rebound.next_event_time_s().unwrap() > 2.0);
    assert_eq!(rebound.cumulative_transfer_counts(), [vec![0], vec![0]]);
}

#[test]
fn g3_flow_and_source_changes_share_partition_and_restart_history_at_common_boundaries() {
    let shared = model(&[0.7, 1.2, 0.0], &[0.1, 0.3]);
    let initial = paths(&[1.0, 0.5, 0.2], &[0.0; 3], &[0.0; 3]);
    let closed = paths(&[0.0; 3], &[0.0; 3], &[0.0; 3]);
    let bare = cohort(
        &shared,
        vec![vec![3, 0, 0], vec![0, 2, 0]],
        initial.clone(),
        89,
        0.0,
    );
    let calendar = schedule(
        &bare,
        initial.clone(),
        vec![
            (0.5, closed.clone()),
            (1.0, paths(&[0.0; 3], &[0.2, 1.0, 0.0], &[0.0; 3])),
            (
                1.5,
                paths(&[0.3, 0.2, 0.0], &[0.7, 0.5, 0.3], &[0.1, 0.05, 0.0]),
            ),
            (2.0, initial),
            (3.0, closed.clone()),
        ],
    );
    let feeds = vec![
        Some(source(
            &shared,
            vec![5.0, 0.0, 0.0],
            vec![
                (0.5, vec![0.0; 3]),
                (1.0, vec![0.0, 3.0, 0.0]),
                (2.0, vec![0.0; 3]),
            ],
        )),
        Some(source(
            &shared,
            vec![0.0; 3],
            vec![
                (1.0, vec![2.0, 0.0, 0.0]),
                (1.5, vec![0.0; 3]),
                (2.0, vec![1.0, 1.0, 0.0]),
            ],
        )),
    ];
    let original = StochasticDecayProductionTransportState::poisson(bare.clone(), feeds.clone())
        .unwrap()
        .with_transport_schedule(Arc::clone(&calendar))
        .unwrap();
    let alternate_order = StochasticDecayProductionTransportState::poisson(
        bare.with_transport_schedule(Arc::clone(&calendar)).unwrap(),
        feeds,
    )
    .unwrap();
    assert_eq!(original, alternate_order);
    let (whole, report) = original
        .trial_to(4.0, produced_budget(), || Ok(()))
        .unwrap();
    let mut split = original.clone();
    let (mut flow_changes, mut source_changes, mut births) = (0, 0, 0);
    for end in [0.1, 0.5, 0.7, 1.0, 1.5, 2.0, 2.25, 3.0, 4.0] {
        let step = split.advance_to(end, produced_budget(), || Ok(())).unwrap();
        flow_changes += step.transport.transport_changes;
        source_changes += step.schedule_changes;
        births += step.births;
        if end == 1.0 {
            let mut restart = split.clone();
            restart
                .advance_to(4.0, produced_budget(), || Ok(()))
                .unwrap();
            assert_eq!(restart, whole);
        }
    }
    assert_eq!(split, whole);
    assert_eq!(
        (flow_changes, source_changes, births),
        (
            report.transport.transport_changes,
            report.schedule_changes,
            report.births
        )
    );
    assert_eq!((flow_changes, source_changes), (5, 6));
    assert_eq!(whole.transport_state().transport(), &closed);
    let retained: u64 = whole
        .inventories()
        .iter()
        .flat_map(|cell| cell.populations())
        .sum();
    let removed: u64 = whole
        .transport_state()
        .cumulative_outflow_counts()
        .iter()
        .flatten()
        .sum();
    assert_eq!(retained + removed, 5 + births);
}

#[test]
fn g3_sub_ulp_movement_precedes_closure_while_exact_ties_and_noop_entries_retain_draw_order() {
    let shared = model(&[0.0], &[]);
    let open = paths(&[1.0], &[0.0], &[0.0]);
    let closed = paths(&[0.0], &[0.0], &[0.0]);
    let start = 2.0_f64.powi(60);
    let close = start.next_up();
    let bare = cohort(&shared, vec![vec![1], vec![0]], open.clone(), 103, start);
    let calendar = schedule(&bare, open, vec![(close, closed.clone())]);
    let mut candidate = bare.with_transport_schedule(calendar).unwrap();
    assert_eq!(candidate.next_event_time_s(), Some(close));
    let report = candidate.advance_to(close, budget(), || Ok(())).unwrap();
    assert_eq!(report.transport_changes, 1);
    assert_eq!(
        report.movement_events, 1,
        "actual arrival is before closure despite a shared rounded endpoint"
    );
    assert_eq!(candidate.inventories()[1].populations(), [1]);

    // With one donor, setting k to its first unit exponential makes the first
    // movement occur exactly at1. A no-op at0.5 must not replace that draw.
    let bare = cohort(&shared, vec![vec![1], vec![0]], closed.clone(), 127, 0.0);
    let mut stream = Stream::resume(bare.movement_rng_checkpoint()).unwrap();
    let rate = stream.next_exponential();
    assert!(rate > 0.0);
    let tied_open = paths(&[rate], &[0.0], &[0.0]);
    let calendar = schedule(
        &bare,
        tied_open.clone(),
        vec![(0.5, tied_open.clone()), (1.0, closed), (2.0, tied_open)],
    );
    let mut tied = bare.with_transport_schedule(calendar).unwrap();
    assert_eq!(tied.movement_rng_checkpoint(), stream.checkpoint());
    let noop = tied.advance_to(0.5, budget(), || Ok(())).unwrap();
    assert_eq!(noop.transport_changes, 1);
    assert_eq!(noop.movement_events, 0);
    assert_eq!(tied.movement_rng_checkpoint(), stream.checkpoint());
    let closure = tied.advance_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!(closure.transport_changes, 1);
    assert_eq!(
        closure.movement_events, 0,
        "the new closed interval owns an exact boundary tie"
    );
    assert_eq!(tied.inventories()[0].populations(), [1]);
    assert_eq!(tied.movement_rng_checkpoint(), stream.checkpoint());
    let next_wait = stream.next_exponential() / rate;
    let reopen = tied.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(reopen.transport_changes, 1);
    assert_eq!(reopen.movement_events, 0);
    assert_eq!(tied.movement_rng_checkpoint(), stream.checkpoint());
    let next = tied.next_event_time_s().unwrap();
    let rounded = 2.0 + next_wait;
    assert!(next > 2.0 && next >= rounded && next <= rounded.next_up());
    assert_eq!(
        tied.advance_to(next, budget(), || Ok(()))
            .unwrap()
            .movement_events,
        1
    );
    assert_eq!(tied.cumulative_transfer_counts(), [vec![1], vec![0]]);
}

#[test]
fn g4_calendar_admission_work_and_early_midway_final_cancellation_preserve_all_history() {
    let shared = model(&[0.0], &[]);
    let initial = paths(&[2.0], &[1.0], &[0.0]);
    let closed = paths(&[0.0], &[0.0], &[0.0]);
    let bare = cohort(&shared, vec![vec![5], vec![5]], initial.clone(), 157, 0.0);
    let change = |time_s| NuclideTransportChange {
        time_s,
        transport: closed.clone(),
    };
    for time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            StochasticDecayTransportSchedule::new(&bare, initial.clone(), vec![change(time)], 4)
                .is_err()
        );
    }
    for times in [[1.0, 1.0], [1.0, 0.5]] {
        assert!(
            StochasticDecayTransportSchedule::new(
                &bare,
                initial.clone(),
                times.map(change).to_vec(),
                4
            )
            .is_err()
        );
    }
    assert!(
        StochasticDecayTransportSchedule::new(&bare, initial.clone(), vec![change(1.0)], 0)
            .is_err()
    );
    let mut changed_topology = closed.clone();
    changed_topology.transfers.swap(0, 1);
    assert!(
        StochasticDecayTransportSchedule::new(
            &bare,
            initial.clone(),
            vec![NuclideTransportChange {
                time_s: 1.0,
                transport: changed_topology,
            }],
            4
        )
        .is_err()
    );
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        let mut invalid = initial.clone();
        invalid.transfers[0].rates_per_s[0] = rate;
        assert!(StochasticDecayTransportSchedule::new(&bare, invalid, vec![], 4).is_err());
    }
    let mut malformed = initial.clone();
    malformed.outflows[0].rates_per_s.clear();
    assert!(StochasticDecayTransportSchedule::new(&bare, malformed, vec![], 4).is_err());
    assert!(matches!(
        StochasticDecayTransportSchedule::new_with_checkpoint(
            &bare,
            initial.clone(),
            vec![change(1.0)],
            4,
            || Err(DecayError::Cancelled)
        ),
        Err(DecayError::Cancelled)
    ));
    let mut admission_calls = 0;
    StochasticDecayTransportSchedule::new_with_checkpoint(
        &bare,
        initial.clone(),
        vec![change(1.0)],
        4,
        || {
            admission_calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        StochasticDecayTransportSchedule::new_with_checkpoint(
            &bare,
            initial.clone(),
            vec![change(1.0)],
            4,
            || {
                seen += 1;
                if seen == admission_calls {
                    Err(DecayError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(DecayError::Cancelled)
    ));
    let calendar = schedule(
        &bare,
        initial.clone(),
        vec![
            (0.5, closed),
            (1.0, paths(&[0.0], &[3.0], &[0.0])),
            (1.5, initial.clone()),
        ],
    );
    let original = bare.with_transport_schedule(Arc::clone(&calendar)).unwrap();
    let mut calls = 0;
    let (_, accepted) = original
        .trial_to(2.0, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(accepted.transport_changes, 3);
    assert!(accepted.movement_events > 0);
    let mut limited = original.clone();
    assert!(
        limited
            .advance_to(
                2.0,
                StochasticDecayTransportBudget {
                    maximum_clock_checks: accepted.clock_checks - 1,
                    ..budget()
                },
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(limited, original);
    for stop in [1, 3, calls / 2, calls] {
        let mut seen = 0;
        let mut candidate = original.clone();
        assert_eq!(
            candidate.advance_to(2.0, budget(), || {
                seen += 1;
                if seen == stop {
                    Err(DecayError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(DecayError::Cancelled)
        );
        assert_eq!(candidate, original);
    }
    let foreign = cohort(
        &model(&[0.0], &[]),
        vec![vec![5], vec![5]],
        initial,
        157,
        0.0,
    );
    assert!(foreign.with_transport_schedule(calendar).is_err());
}

#[test]
fn g4_rebinding_preserves_existing_sources_and_enforces_the_combined_clock_ceiling() {
    let shared = model(&[1.0, 0.0], &[0.1]);
    let feed = source(
        &shared,
        vec![2.0, 0.0],
        vec![(1.0, vec![0.0; 2]), (3.0, vec![1.0, 0.0])],
    );
    let make = |maximum_clocks| {
        let base = cohort(&shared, vec![vec![4, 0]], extraction(0.0), 191, 0.0);
        StochasticDecayTransportState::new(
            base.inventories().to_vec(),
            extraction(0.0),
            base.movement_rng_checkpoint().key,
            StochasticDecayTransportLimits {
                maximum_clocks,
                ..StochasticDecayTransportLimits::default()
            },
        )
        .unwrap()
    };
    // Four original slots + three source slots + one flow-calendar slot.
    let bare = make(8);
    let old = schedule(
        &bare,
        extraction(0.0),
        vec![(1.0, extraction(1.0)), (2.0, extraction(0.0))],
    );
    let replacement = schedule(
        &bare,
        extraction(0.0),
        vec![(1.5, extraction(1.0)), (2.5, extraction(0.0))],
    );
    let unbound =
        StochasticDecayProductionTransportState::poisson(bare, vec![Some(Arc::clone(&feed))])
            .unwrap();
    let mut previous = unbound
        .clone()
        .with_transport_schedule(Arc::clone(&old))
        .unwrap();
    let mut reference = unbound
        .with_transport_schedule(Arc::clone(&replacement))
        .unwrap();
    previous
        .advance_to(0.75, produced_budget(), || Ok(()))
        .unwrap();
    reference
        .advance_to(0.75, produced_budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        previous.clone().with_transport_schedule(old).unwrap(),
        previous
    );
    let mut rebound = previous
        .with_transport_schedule(Arc::clone(&replacement))
        .unwrap();
    assert_eq!(
        rebound, reference,
        "rebinding identical active rates must preserve pending source and nuclear clocks"
    );
    assert!(Arc::ptr_eq(rebound.schedule(0).unwrap(), &feed));
    assert!(Arc::ptr_eq(
        rebound.transport_state().transport_schedule().unwrap(),
        &replacement
    ));
    rebound
        .advance_to(4.0, produced_budget(), || Ok(()))
        .unwrap();
    reference
        .advance_to(4.0, produced_budget(), || Ok(()))
        .unwrap();
    assert_eq!(rebound, reference);

    let tight = make(7);
    let calendar = schedule(&tight, extraction(0.0), vec![(1.0, extraction(1.0))]);
    let produced = StochasticDecayProductionTransportState::poisson(
        tight.clone(),
        vec![Some(Arc::clone(&feed))],
    )
    .unwrap();
    assert!(matches!(
        produced.with_transport_schedule(Arc::clone(&calendar)),
        Err(DecayError::Budget(_))
    ));
    let moving = tight.with_transport_schedule(calendar).unwrap();
    assert!(matches!(
        StochasticDecayProductionTransportState::poisson(moving, vec![Some(feed)]),
        Err(DecayError::Budget(_))
    ));
}
