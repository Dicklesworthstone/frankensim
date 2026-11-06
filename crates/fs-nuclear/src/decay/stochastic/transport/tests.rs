//! G0/G1/G3/G4: physical event competition, exact integer movement,
//! compartment-local emission, replay and atomic refusal. Synthetic data only.

use super::*;
use crate::decay::{DecayBranch, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn key(seed: u64, cell: u32) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_7701,
        tile: cell,
    }
}

fn movement_key(seed: u64) -> StreamKey {
    StreamKey {
        seed,
        kernel: 0xdec0_7702,
        tile: 0,
    }
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

fn model(rates: &[f64], energies: &[f64]) -> Arc<DecayModel> {
    assert_eq!(energies.len() + 1, rates.len());
    Arc::new(
        DecayModel::new(
            "synthetic-stochastic-transport".into(),
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

fn state(
    model: Arc<DecayModel>,
    populations: &[Vec<u64>],
    transport: NuclideTransport,
    seed: u64,
    time: f64,
) -> StochasticDecayTransportState {
    let inventories = populations
        .iter()
        .enumerate()
        .map(|(cell, counts)| {
            StochasticDecayInventory::new_at(
                Arc::clone(&model),
                counts.clone(),
                key(seed, cell as u32),
                time,
            )
            .unwrap()
        })
        .collect();
    StochasticDecayTransportState::new(
        inventories,
        transport,
        movement_key(seed),
        StochasticDecayTransportLimits::default(),
    )
    .unwrap()
}

fn loops() -> NuclideTransport {
    NuclideTransport {
        transfers: vec![
            NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![3.0],
            },
            NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![2.0],
            },
            NuclideTransfer {
                from_compartment: 1,
                to_compartment: 0,
                rates_per_s: vec![5.0],
            },
        ],
        outflows: vec![],
    }
}

#[test]
fn g1_decay_and_outflow_compete_for_the_same_finite_nucleus() {
    let shared = model(&[core::f64::consts::LN_2, 0.0], &[2.0]);
    let samples = 4096_u64;
    let (mut decayed, mut escaped, mut survived) = (0_u64, 0_u64, 0_u64);
    for seed in 0..samples {
        let mut candidate = state(
            Arc::clone(&shared),
            &[vec![1, 0]],
            NuclideTransport {
                transfers: vec![],
                outflows: vec![NuclideOutflow {
                    compartment: 0,
                    rates_per_s: vec![core::f64::consts::LN_2, 0.0],
                }],
            },
            seed,
            0.0,
        );
        let step = candidate.advance_to(1.0, budget(), || Ok(())).unwrap();
        let decay = step.decay[0].events;
        let outflow = step.outflow_counts[0][0];
        let parent = candidate.inventories()[0].populations()[0];
        assert_eq!(decay + outflow + parent, 1);
        assert_eq!(candidate.inventories()[0].populations()[1], decay);
        assert_eq!(step.decay_events, decay);
        assert_eq!(step.movement_events, outflow);
        assert_eq!(step.decay[0].released.deposited_j, 2.0 * decay as f64);
        decayed += decay;
        escaped += outflow;
        survived += parent;
    }
    // Equal competing rates give survival 1/4 at one radioactive half-life,
    // with the other 3/4 divided equally between transformation and outflow.
    for (actual, probability) in [(decayed, 0.375), (escaped, 0.375), (survived, 0.25)] {
        let expected = samples as f64 * probability;
        let allowance = 6.0 * (expected * (1.0 - probability)).sqrt();
        assert!(
            (actual as f64 - expected).abs() < allowance,
            "{actual} versus {expected}"
        );
    }
}

#[test]
fn g1_daughter_moves_at_its_birth_age_and_emits_in_the_receiving_compartment() {
    let mut candidate = state(
        model(&[100.0, 100.0, 0.0], &[3.0, 5.0]),
        &[vec![1, 0, 0], vec![0, 0, 0]],
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![0.0, 1e8, 0.0],
            }],
            outflows: vec![],
        },
        17,
        0.0,
    );
    let step = candidate.advance_to(1.0, budget(), || Ok(())).unwrap();
    assert_eq!(candidate.inventories()[0].populations(), [0, 0, 0]);
    assert_eq!(candidate.inventories()[1].populations(), [0, 0, 1]);
    assert_eq!(step.decay_events, 2);
    assert_eq!(step.movement_events, 1);
    assert_eq!(step.transfer_counts, [vec![0, 1, 0]]);
    assert_eq!(step.decay[0].branch_decays, [1, 0]);
    assert_eq!(step.decay[1].branch_decays, [0, 1]);
    for (cell, heat) in [(0, 3.0), (1, 5.0)] {
        let expected = EnergyChannels {
            deposited_j: heat,
            radiation_j: 0.5 * heat,
            escaped_j: 0.25 * heat,
        };
        assert_eq!(step.decay[cell].released, expected);
        assert_eq!(
            candidate.inventories()[cell].cumulative_released(),
            expected
        );
    }
}

#[test]
fn g1_stable_circulation_retains_each_authored_path_count_without_false_decay() {
    let mut candidate = state(model(&[0.0], &[]), &[vec![1], vec![0]], loops(), 31, 0.0);
    let cell_keys: Vec<_> = candidate
        .inventories()
        .iter()
        .map(|cell| cell.rng_checkpoint())
        .collect();
    let step = candidate.advance_to(20.0, budget(), || Ok(())).unwrap();
    let forward = step.transfer_counts[0][0] + step.transfer_counts[1][0];
    let backward = step.transfer_counts[2][0];
    assert!(step.transfer_counts[0][0] > 0 && step.transfer_counts[1][0] > 0);
    assert!(step.movement_events > 20);
    assert_eq!(step.movement_events, forward + backward);
    assert_eq!(
        forward - backward,
        candidate.inventories()[1].populations()[0]
    );
    assert_eq!(
        candidate.inventories()[0].populations()[0] + candidate.inventories()[1].populations()[0],
        1
    );
    assert_eq!(candidate.cumulative_transfer_counts(), step.transfer_counts);
    assert_eq!(step.decay_events, 0);
    assert_eq!(step.energy_roundoff_j, 0.0);
    for (cell, original) in candidate.inventories().iter().zip(cell_keys) {
        assert_eq!(cell.rng_checkpoint(), original);
        assert_eq!(cell.cumulative_released(), EnergyChannels::default());
    }
}

#[test]
fn g3_observation_partitions_and_whole_cohort_restart_preserve_all_pending_events() {
    let original = state(
        model(&[1.0, 1.0, 0.0], &[0.1, 0.3]),
        &[vec![25, 0, 0], vec![0, 15, 0]],
        NuclideTransport {
            transfers: vec![
                NuclideTransfer {
                    from_compartment: 0,
                    to_compartment: 1,
                    rates_per_s: vec![0.4, 2.0, 0.1],
                },
                NuclideTransfer {
                    from_compartment: 1,
                    to_compartment: 0,
                    rates_per_s: vec![1.2, 0.2, 0.4],
                },
            ],
            outflows: vec![NuclideOutflow {
                compartment: 1,
                rates_per_s: vec![0.1, 0.05, 0.0],
            }],
        },
        71,
        0.0,
    );
    let (unchanged, identity) = original.trial_to(0.0, budget(), || Ok(())).unwrap();
    assert_eq!(unchanged, original);
    assert_eq!(identity.clock_checks, 0);
    assert_eq!(identity.movement_rng_draws, 0);
    assert!(identity.decay.iter().all(|cell| cell.rng_draws == 0));
    let (whole, whole_report) = original.trial_to(8.0, budget(), || Ok(())).unwrap();
    let mut split = original.clone();
    let (mut decays, mut movements) = (0, 0);
    for end in [0.125, 0.5, 1.0, 2.0, 4.0, 8.0] {
        let step = split.advance_to(end, budget(), || Ok(())).unwrap();
        decays += step.decay_events;
        movements += step.movement_events;
        if end == 1.0 {
            let mut restart = split.clone();
            restart.advance_to(8.0, budget(), || Ok(())).unwrap();
            assert_eq!(restart, whole);
        }
    }
    assert_eq!(split, whole);
    assert_eq!(decays, whole_report.decay_events);
    assert_eq!(movements, whole_report.movement_events);
}

#[test]
fn g3_sub_ulp_arrival_precedes_decay_even_when_both_share_one_observation_endpoint() {
    let start = 2.0_f64.powi(60);
    let mut candidate = state(
        model(&[1.0, 0.0], &[2.0]),
        &[vec![0, 0], vec![1, 0]],
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 1,
                to_compartment: 0,
                rates_per_s: vec![1000.0, 0.0],
            }],
            outflows: vec![],
        },
        19,
        start,
    );
    let original = candidate.clone();
    let Some((Event::Movement { .. }, arrival)) = candidate.next_event() else {
        panic!("the faster authored arrival must win this retained clock race");
    };
    let EventClock::Finite(decay) = candidate.inventories[1].clocks[0] else {
        panic!("finite decay");
    };
    assert!(before(arrival, decay));
    assert_eq!(ceiling(arrival), start.next_up());
    assert_eq!(ceiling(decay), start.next_up());
    assert_eq!(candidate.next_event_time_s(), Some(start.next_up()));
    candidate.advance_to(start, budget(), || Ok(())).unwrap();
    assert_eq!(candidate, original);
    let step = candidate
        .advance_to(start.next_up(), budget(), || Ok(()))
        .unwrap();
    assert_eq!(step.movement_events, 1);
    assert_eq!(step.decay_events, 1);
    assert_eq!(step.decay[0].released.deposited_j, 2.0);
    assert_eq!(step.decay[1].released, EnergyChannels::default());
    assert_eq!(candidate.inventories()[0].populations(), [0, 1]);
    assert_eq!(candidate.inventories()[1].populations(), [0, 0]);
}

#[test]
fn g4_cohort_admission_work_energy_rng_and_callback_refusals_are_atomic() {
    let stable = model(&[0.0], &[]);
    let original = state(Arc::clone(&stable), &[vec![1], vec![0]], loops(), 41, 0.0);
    let first = original.next_event_time_s().unwrap();
    for restricted in [
        StochasticDecayTransportBudget {
            maximum_movements: 0,
            ..budget()
        },
        StochasticDecayTransportBudget {
            maximum_clock_checks: 0,
            ..budget()
        },
        StochasticDecayTransportBudget {
            maximum_movements: MAX_STOCHASTIC_NUCLEI + 1,
            ..budget()
        },
    ] {
        let mut candidate = original.clone();
        assert!(candidate.advance_to(first, restricted, || Ok(())).is_err());
        assert_eq!(candidate, original);
    }
    for end in [-1.0, f64::NAN, f64::INFINITY] {
        let mut candidate = original.clone();
        assert!(candidate.advance_to(end, budget(), || Ok(())).is_err());
        assert_eq!(candidate, original);
    }
    let mut calls = 0;
    original
        .trial_to(1.0, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, 3, calls] {
        let mut seen = 0;
        let mut candidate = original.clone();
        assert_eq!(
            candidate.advance_to(1.0, budget(), || {
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
    let mut exhausted = original.clone();
    exhausted.movement_rng.index = u64::MAX;
    let unchanged = exhausted.clone();
    assert_eq!(
        exhausted.advance_to(first, budget(), || Ok(())),
        Err(DecayError::Budget("stochastic random stream exhausted"))
    );
    assert_eq!(exhausted, unchanged);
    let mut overflow = original.clone();
    let Some((Event::Movement { path, node }, _)) = overflow.next_event() else {
        panic!("movement");
    };
    overflow.transfers[path][node] = u64::MAX;
    let unchanged = overflow.clone();
    assert_eq!(
        overflow.advance_to(first, budget(), || Ok(())),
        Err(DecayError::Arithmetic)
    );
    assert_eq!(overflow, unchanged);

    let heated = state(
        model(&[100.0, 0.0], &[0.1]),
        &[vec![3, 0], vec![3, 0]],
        NuclideTransport::default(),
        83,
        0.0,
    );
    let (_, step) = heated.trial_to(1.0, budget(), || Ok(())).unwrap();
    let individual = step
        .decay
        .iter()
        .map(|cell| cell.energy_roundoff_j)
        .fold(0.0_f64, f64::max);
    assert!(individual > 0.0 && step.energy_roundoff_j > individual);
    for allowance in [
        StochasticDecayBudget {
            maximum_events: 3,
            ..budget().decay
        },
        StochasticDecayBudget {
            energy_tolerance_j: individual,
            ..budget().decay
        },
    ] {
        let mut candidate = heated.clone();
        assert!(
            candidate
                .advance_to(
                    1.0,
                    StochasticDecayTransportBudget {
                        decay: allowance,
                        ..budget()
                    },
                    || Ok(())
                )
                .is_err()
        );
        assert_eq!(candidate, heated);
    }

    let cells = original.inventories().to_vec();
    let limits = StochasticDecayTransportLimits::default();
    for restricted in [
        StochasticDecayTransportLimits {
            maximum_compartments: 1,
            ..limits
        },
        StochasticDecayTransportLimits {
            maximum_rate_coefficients: 2,
            ..limits
        },
        StochasticDecayTransportLimits {
            maximum_clocks: 4,
            ..limits
        },
    ] {
        assert!(
            StochasticDecayTransportState::new(
                cells.clone(),
                loops(),
                movement_key(41),
                restricted
            )
            .is_err()
        );
    }
    let duplicate = vec![cells[0].clone(), cells[0].clone()];
    assert!(
        StochasticDecayTransportState::new(duplicate, loops(), movement_key(41), limits).is_err()
    );
    assert!(
        StochasticDecayTransportState::new(
            cells.clone(),
            loops(),
            cells[0].rng_checkpoint().key,
            limits
        )
        .is_err()
    );
    let mut mismatched = cells.clone();
    mismatched[1] = StochasticDecayInventory::new(model(&[0.0], &[]), vec![0], key(41, 1)).unwrap();
    assert!(
        StochasticDecayTransportState::new(mismatched, loops(), movement_key(41), limits).is_err()
    );
    let mut mistimed = cells.clone();
    mistimed[1]
        .advance_to(1.0, budget().decay, || Ok(()))
        .unwrap();
    assert!(
        StochasticDecayTransportState::new(mistimed, loops(), movement_key(41), limits).is_err()
    );
    let too_many = vec![
        StochasticDecayInventory::new(Arc::clone(&stable), vec![MAX_STOCHASTIC_NUCLEI], key(5, 0))
            .unwrap(),
        StochasticDecayInventory::new(stable, vec![1], key(5, 1)).unwrap(),
    ];
    assert!(
        StochasticDecayTransportState::new(too_many, loops(), movement_key(5), limits).is_err()
    );
    for bad_rate in [-1.0, f64::NAN, f64::INFINITY] {
        let mut bad = loops();
        bad.transfers[0].rates_per_s[0] = bad_rate;
        assert!(
            StochasticDecayTransportState::new(cells.clone(), bad, movement_key(41), limits)
                .is_err()
        );
    }
    let mut self_transfer = loops();
    self_transfer.transfers[0].to_compartment = 0;
    assert!(
        StochasticDecayTransportState::new(cells.clone(), self_transfer, movement_key(41), limits)
            .is_err()
    );
    let mut bad_shape = loops();
    bad_shape.transfers[0].rates_per_s.clear();
    assert!(
        StochasticDecayTransportState::new(cells, bad_shape, movement_key(41), limits).is_err()
    );
}
