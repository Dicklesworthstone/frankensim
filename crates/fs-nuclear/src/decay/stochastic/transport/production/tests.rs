//! G0/G1/G3/G4: joint source ages, competing removal, local emission and replay.
//! Every source coefficient is synthetic; expectations below are independent
//! closed-form immigration/transport/decay solutions, not another event loop.

use super::{StochasticDecayProductionTransportBudget, StochasticDecayProductionTransportState};
use crate::EnergyChannels;
use crate::decay::production::{DecayProductionChange, DecayProductionSchedule};
use crate::decay::stochastic::production::{
    StochasticDecayProduction, StochasticDecayProductionBudget,
};
use crate::decay::stochastic::transport::{
    NuclideOutflow, NuclideTransfer, NuclideTransport, StochasticDecayTransportBudget,
    StochasticDecayTransportLimits, StochasticDecayTransportState,
};
use crate::decay::stochastic::{
    MAX_STOCHASTIC_NUCLEI, StochasticDecayBudget, StochasticDecayInventory, StreamKey,
};
use crate::decay::{DecayBranch, DecayError, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_rand::Stream;
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn key(seed: u64, cell: u32) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_8801,
        tile: cell,
    }
}

fn model(rates: &[f64], energies: &[f64]) -> Arc<DecayModel> {
    assert_eq!(energies.len() + 1, rates.len());
    Arc::new(
        DecayModel::new(
            "synthetic-produced-transport".into(),
            rates
                .iter()
                .enumerate()
                .map(|(node, &rate)| DecayNuclide {
                    id: format!("nuclide-{node}"),
                    law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
                })
                .collect(),
            energies
                .iter()
                .enumerate()
                .map(|(parent, &energy)| DecayBranch {
                    parent,
                    daughter: parent + 1,
                    probability: 1.0,
                    energy_per_decay: EnergyChannels {
                        deposited_j: energy,
                        radiation_j: 0.5 * energy,
                        escaped_j: 0.25 * energy,
                    },
                })
                .collect(),
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn schedule(
    model: &Arc<DecayModel>,
    rates: Vec<f64>,
    changes: Vec<(f64, Vec<f64>)>,
) -> Arc<DecayProductionSchedule> {
    Arc::new(
        DecayProductionSchedule::new(
            Arc::clone(model),
            rates,
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

fn cohort(
    model: &Arc<DecayModel>,
    populations: Vec<Vec<u64>>,
    paths: NuclideTransport,
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
                    key(seed, cell as u32),
                    time,
                )
                .unwrap()
            })
            .collect(),
        paths,
        StreamKey {
            seed,
            kernel: 0xdec0_8802,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap()
}

fn budget() -> StochasticDecayProductionTransportBudget {
    StochasticDecayProductionTransportBudget {
        transport: StochasticDecayTransportBudget {
            decay: StochasticDecayBudget {
                maximum_events: 10_000,
                energy_tolerance_j: 1e-8,
            },
            maximum_movements: 10_000,
            maximum_clock_checks: 1_000_000,
        },
        maximum_births: 10_000,
        maximum_schedule_changes: 64,
    }
}

fn mean_close(sum: f64, samples: u64, expected: f64) {
    assert!(expected >= 0.0);
    let actual = sum / samples as f64;
    let allowance = 6.0 * (expected / samples as f64).sqrt() + 2e-12;
    assert!(
        (actual - expected).abs() <= allowance,
        "mean {actual} versus {expected}; allowance {allowance}"
    );
}

fn poisson_moments(sum: f64, squares: f64, samples: u64, expected: f64) {
    mean_close(sum, samples, expected);
    let mean = sum / samples as f64;
    let variance = squares / samples as f64 - mean * mean;
    let allowance = 6.0 * ((expected + 2.0 * expected * expected) / samples as f64).sqrt();
    assert!(
        (variance - expected).abs() < allowance,
        "variance {variance} versus {expected}"
    );
}

#[test]
fn g1_continuous_feed_decay_and_outflow_match_poisson_thinning_and_exact_birth_balance() {
    let (supply, lambda, removal, end) = (5.0, 0.7, 1.3, 2.0);
    let shared = model(&[lambda, 0.0], &[2.0]);
    let feed = schedule(&shared, vec![supply, 0.0], vec![]);
    let hazard: f64 = lambda + removal;
    let parent = supply / hazard * (1.0 - (-hazard * end).exp());
    let transformed_or_removed = supply * end - parent;
    let expected = [
        parent,
        lambda / hazard * transformed_or_removed,
        removal / hazard * transformed_or_removed,
        supply * end,
    ];
    let (mut sums, mut squares) = ([0.0; 4], [0.0; 4]);
    let samples = 4096;
    for seed in 0..samples {
        let transport = cohort(
            &shared,
            vec![vec![0, 0]],
            NuclideTransport {
                transfers: vec![],
                outflows: vec![NuclideOutflow {
                    compartment: 0,
                    rates_per_s: vec![removal, 0.0],
                }],
            },
            seed,
            0.0,
        );
        let mut candidate = StochasticDecayProductionTransportState::poisson(
            transport,
            vec![Some(Arc::clone(&feed))],
        )
        .unwrap();
        let step = candidate.advance_to(end, budget(), || Ok(())).unwrap();
        let counts = candidate.inventories()[0].populations();
        let outflow = step.transport.outflow_counts[0][0];
        assert_eq!(counts.iter().sum::<u64>() + outflow, step.births);
        assert_eq!(step.produced_nuclei, [vec![step.births, 0]]);
        assert_eq!(
            candidate.cumulative_produced_nuclei(0).unwrap(),
            [step.births, 0]
        );
        assert_eq!(step.transport.decay_events, counts[1]);
        assert_eq!(step.transport.movement_events, outflow);
        assert_eq!(
            step.transport.decay[0].released,
            EnergyChannels {
                deposited_j: 2.0 * counts[1] as f64,
                radiation_j: counts[1] as f64,
                escaped_j: 0.5 * counts[1] as f64,
            }
        );
        for (i, value) in [counts[0], counts[1], outflow, step.births]
            .into_iter()
            .enumerate()
        {
            sums[i] += value as f64;
            squares[i] += (value as f64).powi(2);
        }
    }
    // Marking a Poisson arrival by its age and mutually exclusive final fate
    // gives Poisson marginals, including the explicitly untracked exterior.
    for i in 0..4 {
        poisson_moments(sums[i], squares[i], samples, expected[i]);
    }
}

// Unit rates A->B->C in both cells, with only B moving 0->1 at rate1.
// Source feeds A into cell0. For example B0=s/2*(1-exp(-t))^2, while
// B1=s/2*(1-2*t*exp(-t)-exp(-2*t)). C0/C1 integrate these two functions.
fn fed_chain_response(t: f64, supply: f64) -> ([f64; 6], [f64; 4]) {
    let e = (-t).exp();
    let a = supply * (1.0 - e);
    let b0 = 0.5 * supply * (1.0 - e).powi(2);
    let b1 = 0.5 * supply * (1.0 - 2.0 * t * e - e * e);
    let c0 = 0.5 * supply * (t - 2.0 * (1.0 - e) + 0.5 * (1.0 - e * e));
    let c1 = 0.5 * supply * (t - 2.0 * (1.0 - (t + 1.0) * e) - 0.5 * (1.0 - e * e));
    ([a, b0, c0, 0.0, b1, c1], [supply * t - a, c0, 0.0, c1])
}

#[test]
fn g1_fed_daughters_move_and_release_local_energy_before_and_after_shutdown() {
    let (supply, shutdown, end) = (8.0, 1.0, 3.0);
    let shared = model(&[1.0, 1.0, 0.0], &[3.0, 5.0]);
    let feed = schedule(
        &shared,
        vec![supply, 0.0, 0.0],
        vec![(shutdown, vec![0.0; 3])],
    );
    let on = fed_chain_response(shutdown, supply);
    let late = fed_chain_response(end, supply);
    let removed_step = fed_chain_response(end - shutdown, supply);
    let off = (
        std::array::from_fn::<_, 6, _>(|i| late.0[i] - removed_step.0[i]),
        std::array::from_fn::<_, 4, _>(|i| late.1[i] - removed_step.1[i]),
    );
    let (mut populations, mut decays) = ([[0.0; 6]; 2], [[0.0; 4]; 2]);
    let samples = 2048;
    for seed in 0..samples {
        let transport = cohort(
            &shared,
            vec![vec![0; 3], vec![0; 3]],
            NuclideTransport {
                transfers: vec![NuclideTransfer {
                    from_compartment: 0,
                    to_compartment: 1,
                    rates_per_s: vec![0.0, 1.0, 0.0],
                }],
                outflows: vec![],
            },
            seed + 10_000,
            0.0,
        );
        let mut candidate = StochasticDecayProductionTransportState::poisson(
            transport,
            vec![Some(Arc::clone(&feed)), None],
        )
        .unwrap();
        let mut births = 0;
        for (phase, stop) in [shutdown, end].into_iter().enumerate() {
            let step = candidate.advance_to(stop, budget(), || Ok(())).unwrap();
            if phase == 0 {
                births = step.births;
                assert_eq!(step.schedule_changes, 1);
            } else {
                assert_eq!(step.births, 0);
                assert_eq!(step.schedule_changes, 0);
            }
            assert_eq!(
                candidate.cumulative_produced_nuclei(0).unwrap(),
                [births, 0, 0]
            );
            assert_eq!(
                candidate
                    .inventories()
                    .iter()
                    .flat_map(|cell| cell.populations())
                    .sum::<u64>(),
                births
            );
            for (cell, inventory) in candidate.inventories().iter().enumerate() {
                for (node, &count) in inventory.populations().iter().enumerate() {
                    populations[phase][cell * 3 + node] += count as f64;
                }
                for (branch, &count) in inventory.cumulative_branch_decays().iter().enumerate() {
                    decays[phase][cell * 2 + branch] += count as f64;
                }
                let counts = inventory.cumulative_branch_decays();
                let heat = 3.0 * counts[0] as f64 + 5.0 * counts[1] as f64;
                assert_eq!(
                    inventory.cumulative_released(),
                    EnergyChannels {
                        deposited_j: heat,
                        radiation_j: 0.5 * heat,
                        escaped_j: 0.25 * heat,
                    }
                );
            }
        }
    }
    for (phase, expected) in [on, off].into_iter().enumerate() {
        for i in 0..6 {
            mean_close(populations[phase][i], samples, expected.0[i]);
        }
        for i in 0..4 {
            mean_close(decays[phase][i], samples, expected.1[i]);
        }
    }
    assert!(
        decays[0][3] > 0.0,
        "some moved daughters must decay before shutdown"
    );
    assert!(
        decays[1][3] > decays[0][3],
        "daughters keep releasing heat after shutdown"
    );
}

#[test]
fn g3_absent_sources_and_absent_movement_preserve_the_original_owner_histories() {
    let shared = model(&[0.7, 1.2, 0.0], &[0.1, 0.3]);
    let transport = cohort(
        &shared,
        vec![vec![10, 0, 0], vec![0, 5, 0]],
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![0.5, 2.0, 0.1],
            }],
            outflows: vec![NuclideOutflow {
                compartment: 1,
                rates_per_s: vec![0.2, 0.1, 0.0],
            }],
        },
        29,
        0.0,
    );
    let original = transport
        .trial_to(3.0, budget().transport, || Ok(()))
        .unwrap();
    let combined =
        StochasticDecayProductionTransportState::poisson(transport.clone(), vec![None, None])
            .unwrap();
    assert_eq!(combined.transport_state(), &transport);
    let (actual, report) = combined.trial_to(3.0, budget(), || Ok(())).unwrap();
    assert_eq!(actual.transport_state(), &original.0);
    assert_eq!(report.transport, original.1);
    assert_eq!(report.births, 0);
    assert_eq!(report.schedule_changes, 0);
    assert_eq!(report.produced_nuclei, [vec![0; 3], vec![0; 3]]);

    let inventory =
        StochasticDecayInventory::new(Arc::clone(&shared), vec![3, 2, 0], key(61, 0)).unwrap();
    let feed = schedule(
        &shared,
        vec![4.0, 1.0, 0.0],
        vec![
            (0.5, vec![4.0, 1.0, 0.0]),
            (1.0, vec![0.0; 3]),
            (2.0, vec![2.0, 0.0, 0.0]),
        ],
    );
    let mut standalone =
        StochasticDecayProduction::poisson(inventory.clone(), Arc::clone(&feed)).unwrap();
    let bare = StochasticDecayTransportState::new(
        vec![inventory],
        NuclideTransport::default(),
        StreamKey {
            seed: 61,
            kernel: 0xdec0_8802,
            tile: 0,
        },
        StochasticDecayTransportLimits::default(),
    )
    .unwrap();
    let mut combined =
        StochasticDecayProductionTransportState::poisson(bare, vec![Some(Arc::clone(&feed))])
            .unwrap();
    assert!(Arc::ptr_eq(combined.schedule(0).unwrap(), &feed));
    for end in [0.5, 1.0, 1.75, 2.0, 3.0] {
        let simple = standalone
            .advance_to(
                end,
                StochasticDecayProductionBudget {
                    decay: budget().transport.decay,
                    maximum_births: budget().maximum_births,
                    maximum_schedule_changes: budget().maximum_schedule_changes,
                },
                || Ok(()),
            )
            .unwrap();
        let joined = combined.advance_to(end, budget(), || Ok(())).unwrap();
        assert_eq!(&combined.inventories()[0], standalone.inventory());
        assert_eq!(joined.transport.decay[0], simple.decay);
        assert_eq!(joined.produced_nuclei[0], simple.produced_nuclei);
        assert_eq!(joined.births, simple.births);
        assert_eq!(joined.schedule_changes, simple.schedule_changes);
        assert_eq!(
            combined.cumulative_produced_nuclei(0).unwrap(),
            standalone.cumulative_produced_nuclei()
        );
        assert_eq!(
            combined.next_transition_time_s(),
            standalone.next_transition_time_s()
        );
    }
}

#[test]
fn g3_multiple_source_schedules_keep_whole_split_and_restart_trajectories_identical() {
    let shared = model(&[0.7, 1.2, 0.0], &[0.1, 0.3]);
    let schedules = vec![
        Some(schedule(
            &shared,
            vec![3.0, 0.0, 0.0],
            vec![
                (0.75, vec![0.0, 4.0, 0.0]),
                (1.5, vec![0.0; 3]),
                (2.5, vec![2.0, 1.0, 0.0]),
            ],
        )),
        Some(schedule(
            &shared,
            vec![0.0; 3],
            vec![
                (0.5, vec![2.0, 0.0, 0.0]),
                (1.5, vec![2.0, 0.0, 0.0]),
                (3.0, vec![0.0; 3]),
            ],
        )),
    ];
    let original = StochasticDecayProductionTransportState::poisson(
        cohort(
            &shared,
            vec![vec![2, 3, 0], vec![0; 3]],
            NuclideTransport {
                transfers: vec![
                    NuclideTransfer {
                        from_compartment: 0,
                        to_compartment: 1,
                        rates_per_s: vec![0.5, 2.0, 0.1],
                    },
                    NuclideTransfer {
                        from_compartment: 1,
                        to_compartment: 0,
                        rates_per_s: vec![1.0, 0.25, 0.6],
                    },
                ],
                outflows: vec![NuclideOutflow {
                    compartment: 1,
                    rates_per_s: vec![0.15, 0.1, 0.0],
                }],
            },
            113,
            0.0,
        ),
        schedules,
    )
    .unwrap();
    let (identity, zero) = original.trial_to(0.0, budget(), || Ok(())).unwrap();
    assert_eq!(identity, original);
    assert_eq!(
        zero.births + zero.transport.decay_events + zero.transport.movement_events,
        0
    );
    assert_eq!(zero.transport.clock_checks, 0);
    let (whole, report) = original.trial_to(5.0, budget(), || Ok(())).unwrap();
    let mut split = original.clone();
    let (mut births, mut changes, mut decays, mut movements) = (0, 0, 0, 0);
    for end in [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0] {
        let step = split.advance_to(end, budget(), || Ok(())).unwrap();
        births += step.births;
        changes += step.schedule_changes;
        decays += step.transport.decay_events;
        movements += step.transport.movement_events;
        if end == 1.5 {
            let mut restart = split.clone();
            restart.advance_to(5.0, budget(), || Ok(())).unwrap();
            assert_eq!(restart, whole);
        }
    }
    assert_eq!(split, whole);
    assert_eq!(
        (births, changes, decays, movements),
        (
            report.births,
            report.schedule_changes,
            report.transport.decay_events,
            report.transport.movement_events
        )
    );
    assert_eq!(changes, 6);
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
fn g3_sub_ulp_birth_transfer_and_decay_precede_shutdown_but_exact_tied_birth_does_not() {
    let shared = model(&[4.0, 0.0], &[2.0]);
    let paths = NuclideTransport {
        transfers: vec![NuclideTransfer {
            from_compartment: 0,
            to_compartment: 1,
            rates_per_s: vec![1000.0, 0.0],
        }],
        outflows: vec![],
    };
    let start = 2.0_f64.powi(60);
    let shutdown = start.next_up();
    let feed = schedule(&shared, vec![1.0, 0.0], vec![(shutdown, vec![0.0, 0.0])]);
    let original = StochasticDecayProductionTransportState::poisson(
        cohort(
            &shared,
            vec![vec![0, 0], vec![0, 0]],
            paths.clone(),
            73,
            start,
        ),
        vec![Some(feed), None],
    )
    .unwrap();
    assert_eq!(original.next_transition_time_s(), Some(shutdown));
    let (identity, zero) = original.trial_to(start, budget(), || Ok(())).unwrap();
    assert_eq!(identity, original);
    assert_eq!(zero.births, 0);
    let mut split = original.clone();
    let before_shutdown = split.advance_to(shutdown, budget(), || Ok(())).unwrap();
    assert!(before_shutdown.births > 0);
    assert!(before_shutdown.transport.movement_events > 0);
    assert!(before_shutdown.transport.decay_events > 0);
    assert_eq!(before_shutdown.schedule_changes, 1);
    let after_shutdown = split
        .advance_to(shutdown.next_up(), budget(), || Ok(()))
        .unwrap();
    assert_eq!(after_shutdown.births, 0);
    assert_eq!(
        split.cumulative_produced_nuclei(0).unwrap()[0],
        before_shutdown.births
    );
    assert!(
        split
            .inventories()
            .iter()
            .all(|cell| cell.populations()[0] == 0)
    );
    assert_eq!(
        split
            .inventories()
            .iter()
            .map(|cell| cell.populations()[1])
            .sum::<u64>(),
        before_shutdown.births
    );
    assert_eq!(
        split
            .inventories()
            .iter()
            .map(|cell| cell.cumulative_released().deposited_j)
            .sum::<f64>(),
        2.0 * before_shutdown.births as f64
    );
    let whole = original
        .trial_to(shutdown.next_up(), budget(), || Ok(()))
        .unwrap()
        .0;
    assert_eq!(split, whole);

    // Manufacture an exact, fully admitted clock tie using the public stream:
    // unit_wait/source_rate =1, so the first possible birth is exactly t=1.
    let bare = cohort(&shared, vec![vec![0, 0], vec![0, 0]], paths, 97, 0.0);
    let source_rate = Stream::resume(bare.inventories()[0].rng_checkpoint())
        .unwrap()
        .next_exponential();
    assert!(source_rate > 0.0);
    let stop_at_birth = schedule(&shared, vec![source_rate, 0.0], vec![(1.0, vec![0.0, 0.0])]);
    let mut tied =
        StochasticDecayProductionTransportState::poisson(bare, vec![Some(stop_at_birth), None])
            .unwrap();
    assert_eq!(tied.next_transition_time_s(), Some(1.0));
    let report = tied.advance_to(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(report.schedule_changes, 1);
    assert_eq!(
        report.births, 0,
        "old supply is active on the half-open interval [0,1)"
    );
    assert_eq!(
        report.transport.decay_events + report.transport.movement_events,
        0
    );
    assert!(tied.next_transition_time_s().is_none());
}

#[test]
fn g4_global_source_movement_schedule_and_clock_limits_and_cancellation_are_atomic() {
    let shared = model(&[0.0], &[]);
    let original = StochasticDecayProductionTransportState::poisson(
        cohort(
            &shared,
            vec![vec![0], vec![0]],
            NuclideTransport {
                transfers: vec![
                    NuclideTransfer {
                        from_compartment: 0,
                        to_compartment: 1,
                        rates_per_s: vec![5.0],
                    },
                    NuclideTransfer {
                        from_compartment: 1,
                        to_compartment: 0,
                        rates_per_s: vec![4.0],
                    },
                ],
                outflows: vec![],
            },
            151,
            0.0,
        ),
        vec![
            Some(schedule(&shared, vec![20.0], vec![(0.5, vec![20.0])])),
            Some(schedule(&shared, vec![15.0], vec![(0.5, vec![15.0])])),
        ],
    )
    .unwrap();
    let mut calls = 0;
    let (_, reference) = original
        .trial_to(5.0, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let maximum_cell_births = reference
        .produced_nuclei
        .iter()
        .map(|row| row.iter().sum::<u64>())
        .max()
        .unwrap();
    let maximum_path_movements = reference
        .transport
        .transfer_counts
        .iter()
        .flatten()
        .copied()
        .max()
        .unwrap();
    assert!(reference.births > maximum_cell_births);
    assert!(reference.transport.movement_events > maximum_path_movements);
    assert_eq!(reference.schedule_changes, 2);
    for restricted in [
        StochasticDecayProductionTransportBudget {
            maximum_births: maximum_cell_births,
            ..budget()
        },
        StochasticDecayProductionTransportBudget {
            maximum_schedule_changes: 1,
            ..budget()
        },
        StochasticDecayProductionTransportBudget {
            maximum_births: MAX_STOCHASTIC_NUCLEI + 1,
            ..budget()
        },
        StochasticDecayProductionTransportBudget {
            transport: StochasticDecayTransportBudget {
                maximum_movements: maximum_path_movements,
                ..budget().transport
            },
            ..budget()
        },
        StochasticDecayProductionTransportBudget {
            transport: StochasticDecayTransportBudget {
                maximum_clock_checks: 4,
                ..budget().transport
            },
            ..budget()
        },
    ] {
        let mut candidate = original.clone();
        assert!(candidate.advance_to(5.0, restricted, || Ok(())).is_err());
        assert_eq!(candidate, original);
    }
    for end in [-1.0, f64::NAN, f64::INFINITY] {
        let mut candidate = original.clone();
        assert!(candidate.advance_to(end, budget(), || Ok(())).is_err());
        assert_eq!(candidate, original);
    }
    for stop in [1, 3, calls / 2, calls] {
        let mut seen = 0;
        let mut candidate = original.clone();
        assert_eq!(
            candidate.advance_to(5.0, budget(), || {
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
}

#[test]
fn g4_source_admission_and_whole_cohort_population_decay_and_energy_allowances_are_explicit() {
    let shared = model(&[100.0, 0.0], &[0.1]);
    let zero_feed = schedule(&shared, vec![0.0, 0.0], vec![]);
    let original = StochasticDecayProductionTransportState::poisson(
        cohort(
            &shared,
            vec![vec![3, 0], vec![3, 0]],
            NuclideTransport::default(),
            211,
            0.0,
        ),
        vec![Some(Arc::clone(&zero_feed)), Some(Arc::clone(&zero_feed))],
    )
    .unwrap();
    let (_, accepted) = original.trial_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!(accepted.transport.decay_events, 6);
    let individual = accepted
        .transport
        .decay
        .iter()
        .map(|cell| cell.energy_roundoff_j)
        .fold(0.0_f64, f64::max);
    assert!(individual > 0.0 && accepted.transport.energy_roundoff_j > individual);
    for local in [
        StochasticDecayBudget {
            maximum_events: 3,
            ..budget().transport.decay
        },
        StochasticDecayBudget {
            energy_tolerance_j: individual,
            ..budget().transport.decay
        },
        StochasticDecayBudget {
            energy_tolerance_j: f64::NAN,
            ..budget().transport.decay
        },
    ] {
        let mut candidate = original.clone();
        assert!(
            candidate
                .advance_to(
                    1.0,
                    StochasticDecayProductionTransportBudget {
                        transport: StochasticDecayTransportBudget {
                            decay: local,
                            ..budget().transport
                        },
                        ..budget()
                    },
                    || Ok(())
                )
                .is_err()
        );
        assert_eq!(candidate, original);
    }

    let stable = model(&[0.0], &[]);
    let feed = schedule(&stable, vec![1.0], vec![]);
    let mut full = StochasticDecayProductionTransportState::poisson(
        cohort(
            &stable,
            vec![vec![MAX_STOCHASTIC_NUCLEI - 1], vec![1]],
            NuclideTransport::default(),
            241,
            0.0,
        ),
        vec![Some(Arc::clone(&feed)), None],
    )
    .unwrap();
    let unchanged = full.clone();
    let birth = full.next_transition_time_s().unwrap();
    assert!(full.advance_to(birth, budget(), || Ok(())).is_err());
    assert_eq!(full, unchanged);

    let bare = cohort(
        &stable,
        vec![vec![0], vec![0]],
        NuclideTransport::default(),
        251,
        0.0,
    );
    assert!(
        StochasticDecayProductionTransportState::poisson(
            bare.clone(),
            vec![Some(Arc::clone(&feed))]
        )
        .is_err()
    );
    let distinct = model(&[0.0], &[]);
    let wrong = schedule(&distinct, vec![1.0], vec![]);
    assert!(
        StochasticDecayProductionTransportState::poisson(bare.clone(), vec![Some(wrong), None])
            .is_err()
    );
    let tight = StochasticDecayTransportState::new(
        bare.inventories().to_vec(),
        NuclideTransport::default(),
        bare.movement_rng_checkpoint().key,
        StochasticDecayTransportLimits {
            maximum_clocks: 2,
            ..StochasticDecayTransportLimits::default()
        },
    )
    .unwrap();
    assert!(
        StochasticDecayProductionTransportState::poisson(tight.clone(), vec![None, None]).is_ok()
    );
    assert!(
        StochasticDecayProductionTransportState::poisson(tight, vec![Some(feed), None]).is_err()
    );
}
