//! Analytic joint reaction/transport verification; all coefficients manufactured.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{BurnupData, BurnupNuclide, BurnupProduct, BurnupRate, BurnupReaction};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, allowance {tolerance:.3e}"
    );
}

fn limits() -> BurnupTransportLimits {
    BurnupTransportLimits::default()
}

fn budget() -> BurnupTransportBudget {
    BurnupTransportBudget {
        limits: limits(),
        nuclei_tolerance: 1e-6,
        energy_tolerance_j: 1e-6,
    }
}

fn local_budget() -> BurnupBudget {
    budget().local().unwrap()
}

fn integral(rate: f64, duration: f64) -> f64 {
    -(-rate * duration).exp_m1() / rate
}

fn grid() -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("manufactured-transport-grid".into(), vec![0.0, 1.0]).unwrap())
}

fn model(rate: f64) -> Arc<BurnupModel> {
    let mut nuclides = vec![BurnupNuclide {
        id: NuclideId::new(10, 25, 0).unwrap(),
        law: if rate == 0.0 {
            DecayLawV1::stable()
        } else {
            DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap()
        },
    }];
    let mut reactions = Vec::new();
    if rate > 0.0 {
        nuclides.push(BurnupNuclide {
            id: NuclideId::new(11, 25, 0).unwrap(),
            law: DecayLawV1::stable(),
        });
        reactions.push(BurnupReaction {
            data_id: "manufactured-parent-decay".into(),
            parent: 0,
            rate: BurnupRate::Decay {
                branch_fraction: 1.0,
            },
            products: vec![BurnupProduct {
                nuclide: 1,
                yield_per_event: 1.0,
            }],
            energy_per_event: EnergyChannels {
                deposited_j: 3.0,
                radiation_j: 2.0,
                escaped_j: 1.0,
            },
        });
    }
    Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-transport-model".into(),
                groups: grid(),
                temperature_min_k: 200.0,
                temperature_max_k: 1000.0,
                nuclides,
                reactions,
            },
            limits().burnup,
        )
        .unwrap(),
    )
}

fn zero_flux(model: &BurnupModel) -> GroupFlux {
    GroupFlux::new(Arc::clone(model.groups()), vec![0.0]).unwrap()
}

fn state(
    model: Arc<BurnupModel>,
    populations: Vec<Vec<f64>>,
    volumes: Vec<f64>,
    transport: NuclideTransport,
) -> BurnupTransportState {
    BurnupTransportState::new(
        populations
            .into_iter()
            .map(|row| BurnupInventory::new(Arc::clone(&model), row).unwrap())
            .collect(),
        volumes,
        transport,
        limits(),
        || Ok(()),
    )
    .unwrap()
}

fn flow(from: usize, to: usize, rates: Vec<f64>) -> NuclideTransfer {
    NuclideTransfer {
        from_compartment: from,
        to_compartment: to,
        rates_per_s: rates,
    }
}

fn advance(
    initial: &BurnupTransportState,
    duration: f64,
) -> (BurnupTransportState, BurnupTransportStep) {
    let cells = initial.inventories().len();
    let fluxes = vec![zero_flux(initial.model()); cells];
    initial
        .trial(duration, &fluxes, &vec![500.0; cells], budget(), || Ok(()))
        .unwrap()
}

#[test]
fn g1_finite_donor_moves_and_decays_in_both_compartments_without_moving_past_heat() {
    let model = model(1.0);
    let flux = zero_flux(&model);
    let mut donor = BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0]).unwrap();
    donor
        .advance(0.25, &flux, 500.0, local_budget(), || Ok(()))
        .unwrap();
    let before = donor.clone();
    let initial_parent = donor.populations()[0];
    let initial_daughter = donor.populations()[1];
    let receiver = BurnupInventory::new(Arc::clone(&model), vec![0.0, 0.0]).unwrap();
    let state = BurnupTransportState::new(
        vec![donor, receiver],
        vec![1.0, 2.0],
        NuclideTransport {
            transfers: vec![flow(0, 1, vec![1.0, 0.0])],
            outflows: vec![],
        },
        limits(),
        || Ok(()),
    )
    .unwrap();
    let (next, step) = advance(&state, 1.0);
    let e1 = (-1.0_f64).exp();
    let e2 = (-2.0_f64).exp();
    let donor_decays = initial_parent * (1.0 - e2) / 2.0;
    let receiver_decays = initial_parent * (1.0 - e1) - donor_decays;
    near(
        next.inventories()[0].populations()[0],
        initial_parent * e2,
        1e-9,
    );
    near(
        next.inventories()[1].populations()[0],
        initial_parent * (e1 - e2),
        1e-9,
    );
    near(
        next.inventories()[0].populations()[1],
        initial_daughter + donor_decays,
        1e-9,
    );
    near(
        next.inventories()[1].populations()[1],
        receiver_decays,
        1e-9,
    );
    near(step.transfer_counts[0][0], donor_decays, 1e-9);
    assert_eq!(step.transfer_counts[0][1], 0.0);
    near(step.burnup[0].parent_removals[0], donor_decays, 1e-9);
    near(step.burnup[1].parent_removals[0], receiver_decays, 1e-9);
    near(
        step.burnup[0].released.deposited_j,
        3.0 * donor_decays,
        1e-8,
    );
    near(
        step.burnup[1].released.deposited_j,
        3.0 * receiver_decays,
        1e-8,
    );
    near(
        next.inventories()[0].cumulative_released().deposited_j,
        before.cumulative_released().deposited_j + 3.0 * donor_decays,
        1e-8,
    );
    near(
        next.inventories()[1].cumulative_released().deposited_j,
        3.0 * receiver_decays,
        1e-8,
    );
    assert_eq!(step.produced_nuclei, vec![vec![0.0, 0.0]; 2]);
    assert_eq!(
        next.inventories()[0].cumulative_produced_nuclei(),
        vec![0.0, 0.0]
    );
    assert_eq!(state.inventories()[0], before);
    assert!(
        step.cohort_balance_residuals
            .iter()
            .all(|value| value.abs() < 1e-8)
    );
}

#[test]
fn g1_unequal_volume_recirculation_and_exterior_extraction_close_the_finite_inventory() {
    let model = model(0.0);
    let flow_rate = 1.0;
    let transport = NuclideTransport {
        // Q=1 m3/s, donor volumes1 and2m3: equal flow, different rate constants.
        transfers: vec![
            flow(0, 1, vec![flow_rate]),
            flow(1, 0, vec![flow_rate / 2.0]),
        ],
        outflows: vec![],
    };
    let initial = state(
        Arc::clone(&model),
        vec![vec![12.0], vec![0.0]],
        vec![1.0, 2.0],
        transport.clone(),
    );
    let duration = 2.0;
    let mixing = (-1.5_f64 * duration).exp();
    let (next, report) = advance(&initial, duration);
    near(
        next.inventories()[0].populations()[0],
        4.0 + 8.0 * mixing,
        1e-9,
    );
    near(
        next.inventories()[1].populations()[0],
        8.0 - 8.0 * mixing,
        1e-9,
    );
    near(
        report.transfer_counts[0][0],
        4.0 * duration + 8.0 * integral(1.5, duration),
        1e-9,
    );
    near(
        report.transfer_counts[1][0],
        4.0 * duration - 4.0 * integral(1.5, duration),
        1e-9,
    );
    assert!(
        report.transfer_counts[0][0] > 12.0,
        "recirculation counts repeated visits"
    );
    assert_eq!(report.burnup[0].decays, 0.0);
    assert_eq!(report.burnup[1].decays, 0.0);
    assert_eq!(report.burnup[0].released, EnergyChannels::default());

    let removal = 0.2;
    let extracting = NuclideTransport {
        outflows: vec![
            NuclideOutflow {
                compartment: 0,
                rates_per_s: vec![removal],
            },
            NuclideOutflow {
                compartment: 1,
                rates_per_s: vec![removal],
            },
        ],
        ..transport
    };
    let initial = state(
        model,
        vec![vec![12.0], vec![0.0]],
        vec![1.0, 2.0],
        extracting,
    );
    let (next, report) = advance(&initial, duration);
    let survival = (-removal * duration).exp();
    near(
        next.inventories()[0].populations()[0],
        survival * (4.0 + 8.0 * mixing),
        1e-9,
    );
    near(
        next.inventories()[1].populations()[0],
        survival * (8.0 - 8.0 * mixing),
        1e-9,
    );
    let lost_a =
        removal * (4.0 * integral(removal, duration) + 8.0 * integral(removal + 1.5, duration));
    let lost_b =
        removal * (8.0 * integral(removal, duration) - 8.0 * integral(removal + 1.5, duration));
    near(report.outflow_counts[0][0], lost_a, 1e-9);
    near(report.outflow_counts[1][0], lost_b, 1e-9);
    near(lost_a + lost_b, 12.0 * (1.0 - survival), 1e-9);
    near(report.cohort_balance_residuals[0], 0.0, 1e-8);
    assert_eq!(report.produced_nuclei_by_compartment, vec![0.0; 2]);
}

#[test]
fn g1_continuous_feed_moves_and_decays_at_its_actual_age() {
    let model = model(1.0);
    let transport = NuclideTransport {
        transfers: vec![flow(0, 1, vec![1.0, 0.0])],
        outflows: vec![],
    };
    let initial = state(
        model,
        vec![vec![0.0; 2]; 2],
        vec![1.0; 2],
        transport.clone(),
    );
    let fluxes = vec![zero_flux(initial.model()); 2];
    let production = vec![vec![10.0, 0.0], vec![0.0, 0.0]];
    let (next, report) = initial
        .trial_with_inputs(
            1.0,
            &fluxes,
            &[500.0; 2],
            &transport,
            &production,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let parent_a = 10.0 * integral(2.0, 1.0);
    let parent_total = 10.0 * integral(1.0, 1.0);
    let decays_a = 5.0 * (1.0 - integral(2.0, 1.0));
    let decays_b = 10.0 - parent_total - decays_a;
    near(next.inventories()[0].populations()[0], parent_a, 1e-9);
    near(
        next.inventories()[1].populations()[0],
        parent_total - parent_a,
        1e-9,
    );
    near(report.burnup[0].decays, decays_a, 1e-9);
    near(report.burnup[1].decays, decays_b, 1e-9);
    near(report.transfer_counts[0][0], decays_a, 1e-9);
    near(report.burnup[1].released.deposited_j, 3.0 * decays_b, 1e-8);
    assert_eq!(report.produced_nuclei, production);
    assert_eq!(report.produced_nuclei_by_compartment, vec![10.0, 0.0]);
    assert_eq!(
        next.inventories()[0].cumulative_produced_nuclei(),
        vec![10.0, 0.0]
    );
    assert!(next.prepared().unwrap().supports_sources());
}

#[test]
fn g1_new_fission_products_move_and_decay_while_capture_competes() {
    let groups = grid();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-moving-products".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 200.0,
                temperature_max_k: 1000.0,
                nuclides: vec![
                    BurnupNuclide {
                        id: NuclideId::new(92, 235, 0).unwrap(),
                        law: DecayLawV1::stable(),
                    },
                    BurnupNuclide {
                        id: NuclideId::new(40, 100, 0).unwrap(),
                        law: DecayLawV1::from_decay_constant(Frequency::new(1.0)).unwrap(),
                    },
                    BurnupNuclide {
                        id: NuclideId::new(41, 100, 0).unwrap(),
                        law: DecayLawV1::stable(),
                    },
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-fission".into(),
                        parent: 0,
                        rate: BurnupRate::Fission {
                            cross_sections_m2: vec![1.0],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 2.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 3.0,
                            ..EnergyChannels::default()
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-capture".into(),
                        parent: 0,
                        rate: BurnupRate::Capture {
                            cross_sections_m2: vec![1.0],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 1.0,
                            ..EnergyChannels::default()
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-product-decay".into(),
                        parent: 1,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            ..EnergyChannels::default()
                        },
                    },
                ],
            },
            limits().burnup,
        )
        .unwrap(),
    );
    let initial = state(
        model,
        vec![vec![100.0, 0.0, 0.0], vec![0.0; 3]],
        vec![1.0; 2],
        NuclideTransport {
            transfers: vec![flow(0, 1, vec![0.0, 1.0, 0.0])],
            outflows: vec![],
        },
    );
    let flux = GroupFlux::new(groups, vec![1.0]).unwrap();
    let (next, report) = initial
        .trial(1.0, &[flux.clone(), flux], &[500.0; 2], budget(), || Ok(()))
        .unwrap();
    let e1 = (-1.0_f64).exp();
    let e2 = (-2.0_f64).exp();
    let fissions = 50.0 * (1.0 - e2);
    let products_total = 200.0 * (e1 - e2);
    let product_a = 200.0 * e2;
    let decays_a = 50.0 * (1.0 - 3.0 * e2);
    let decays_b = 2.0 * fissions - products_total - decays_a;
    near(next.inventories()[0].populations()[0], 100.0 * e2, 1e-8);
    near(next.inventories()[0].populations()[1], product_a, 1e-8);
    near(
        next.inventories()[1].populations()[1],
        products_total - product_a,
        1e-8,
    );
    near(report.burnup[0].fissions, fissions, 1e-8);
    near(report.burnup[0].captures, fissions, 1e-8);
    near(report.burnup[0].decays, decays_a, 1e-8);
    near(report.burnup[1].decays, decays_b, 1e-8);
    near(report.burnup[0].parent_removals[1], decays_a, 1e-8);
    near(report.transfer_counts[0][1], decays_a, 1e-8);
    near(
        report.burnup[0].released.deposited_j,
        4.0 * fissions + 2.0 * decays_a,
        1e-7,
    );
    near(report.burnup[1].released.deposited_j, 2.0 * decays_b, 1e-7);
    let remaining = next
        .inventories()
        .iter()
        .flat_map(|inventory| inventory.populations())
        .sum::<f64>();
    near(remaining, 100.0 + fissions, 1e-7);
    assert!(
        report
            .cohort_balance_residuals
            .iter()
            .all(|value| value.abs() < 1e-7)
    );
}

#[test]
fn g3_source_changes_reuse_map_and_rate_changes_invalidate_without_losing_histories() {
    let model = model(1.0);
    let transport = NuclideTransport {
        transfers: vec![flow(0, 1, vec![0.5, 0.5]), flow(1, 0, vec![0.25, 0.25])],
        outflows: vec![],
    };
    let initial = state(
        model,
        vec![vec![10.0, 0.0], vec![0.0; 2]],
        vec![1.0, 2.0],
        transport.clone(),
    );
    let fluxes = vec![zero_flux(initial.model()); 2];
    let feed = vec![vec![3.0, 0.0], vec![0.0, 0.0]];
    let (first, _) = initial
        .trial_with_inputs(
            0.25,
            &fluxes,
            &[500.0; 2],
            &transport,
            &feed,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let map = Arc::clone(first.prepared().unwrap());
    let more_feed = vec![vec![9.0, 0.0], vec![0.0, 1.0]];
    let (second, _) = first
        .trial_with_inputs(
            0.25,
            &fluxes,
            &[700.0; 2],
            &transport,
            &more_feed,
            budget(),
            || Ok(()),
        )
        .unwrap();
    assert!(
        Arc::ptr_eq(&map, second.prepared().unwrap()),
        "source magnitudes and equivalent admitted rates reuse one map"
    );
    let (cooling, _) = second
        .trial(0.25, &fluxes, &[600.0; 2], budget(), || Ok(()))
        .unwrap();
    assert!(
        Arc::ptr_eq(&map, cooling.prepared().unwrap()),
        "source-capable map also serves zero feed"
    );
    assert_eq!(
        cooling.inventories()[0].cumulative_produced_nuclei(),
        vec![3.0, 0.0]
    );
    assert_eq!(
        cooling.inventories()[1].cumulative_produced_nuclei(),
        vec![0.0, 0.25]
    );
    let mut changed = transport.clone();
    changed.transfers[0].rates_per_s[0] = 0.8;
    let (different, _) = second
        .trial_with_inputs(
            0.25,
            &fluxes,
            &[500.0; 2],
            &changed,
            &feed,
            budget(),
            || Ok(()),
        )
        .unwrap();
    assert!(!Arc::ptr_eq(&map, different.prepared().unwrap()));
    assert_eq!(different.transport(), &changed);
    assert!(
        different.cumulative_transfer_counts()[0][0] > second.cumulative_transfer_counts()[0][0]
    );
    assert_eq!(second.transport(), &transport);
}

#[test]
fn g3_joint_semigroup_restart_preserves_each_cell_and_transfer_account() {
    let model = model(0.4);
    let initial = state(
        model,
        vec![vec![50.0, 2.0], vec![4.0, 7.0]],
        vec![1.0; 2],
        NuclideTransport {
            transfers: vec![flow(0, 1, vec![0.7, 0.2]), flow(1, 0, vec![0.3, 0.1])],
            outflows: vec![NuclideOutflow {
                compartment: 1,
                rates_per_s: vec![0.1, 0.0],
            }],
        },
    );
    let fluxes = vec![zero_flux(initial.model()); 2];
    let production = vec![vec![1.0, 0.0], vec![2.0, 0.5]];
    let (whole, whole_report) = initial
        .trial_with_inputs(
            1.0,
            &fluxes,
            &[500.0; 2],
            initial.transport(),
            &production,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let (half, _) = initial
        .trial_with_inputs(
            0.5,
            &fluxes,
            &[500.0; 2],
            initial.transport(),
            &production,
            budget(),
            || Ok(()),
        )
        .unwrap();
    let checkpoint = half.clone();
    let (split, split_report) = checkpoint
        .trial_with_inputs(
            0.5,
            &fluxes,
            &[500.0; 2],
            checkpoint.transport(),
            &production,
            budget(),
            || Ok(()),
        )
        .unwrap();
    for (a, b) in split.inventories().iter().zip(whole.inventories()) {
        for (&x, &y) in a.populations().iter().zip(b.populations()) {
            near(x, y, 1e-8);
        }
        for (x, y) in a
            .cumulative_reaction_counts()
            .into_iter()
            .zip(b.cumulative_reaction_counts())
        {
            near(x, y, 1e-8);
        }
        near(
            a.cumulative_released().total_j(),
            b.cumulative_released().total_j(),
            1e-7,
        );
    }
    for (a, b) in split
        .cumulative_transfer_counts()
        .iter()
        .flatten()
        .zip(whole.cumulative_transfer_counts().iter().flatten())
    {
        near(*a, *b, 1e-8);
    }
    for (a, b) in split
        .cumulative_outflow_counts()
        .iter()
        .flatten()
        .zip(whole.cumulative_outflow_counts().iter().flatten())
    {
        near(*a, *b, 1e-8);
    }
    assert!(split_report.population_roundoff_nuclei <= budget().nuclei_tolerance);
    near(
        whole_report.population_roundoff_nuclei,
        whole_report
            .burnup
            .iter()
            .map(|report| report.population_roundoff_nuclei)
            .sum(),
        1e-20,
    );
    let (idle, idle_report) = initial
        .trial(0.0, &fluxes, &[500.0; 2], budget(), || Ok(()))
        .unwrap();
    assert_eq!(idle, initial);
    assert_eq!(idle_report.transfer_counts, vec![vec![0.0; 2]; 2]);
}

#[test]
fn g1_slow_decay_survives_a_million_times_faster_reciprocal_circulation() {
    let lambda = 1e-6;
    let initial = state(
        model(lambda),
        vec![vec![1000.0, 0.0], vec![0.0; 2]],
        vec![1.0; 2],
        NuclideTransport {
            transfers: vec![flow(0, 1, vec![1.0; 2]), flow(1, 0, vec![1.0; 2])],
            outflows: vec![],
        },
    );
    let (next, report) = advance(&initial, 1.0);
    let fraction_a = (1.0 + (-2.0_f64).exp()) / 2.0;
    let parent_total = 1000.0 * (-lambda).exp();
    let total_decays = 1000.0 * -(-lambda).exp_m1();
    let donor_decays = 500.0 * lambda * (integral(lambda, 1.0) + integral(lambda + 2.0, 1.0));
    near(
        next.inventories()[0].populations()[0],
        parent_total * fraction_a,
        1e-8,
    );
    near(
        next.inventories()[1].populations()[0],
        parent_total * (1.0 - fraction_a),
        1e-8,
    );
    near(
        next.inventories()[0].populations()[1],
        total_decays * fraction_a,
        1e-10,
    );
    near(report.burnup[0].decays, donor_decays, 2e-11);
    near(report.burnup[1].decays, total_decays - donor_decays, 2e-11);
    assert!(report.burnup[0].decays > 0.0 && report.burnup[1].decays > 0.0);
    near(
        report.burnup[0].released.deposited_j + report.burnup[1].released.deposited_j,
        3.0 * total_decays,
        6e-11,
    );
}

#[test]
fn g0_work_model_shape_domain_and_callback_refusals_leave_the_whole_cohort_unchanged() {
    let model = model(1.0);
    let transport = NuclideTransport {
        transfers: vec![flow(0, 1, vec![1.0, 0.0])],
        outflows: vec![],
    };
    let initial = state(
        Arc::clone(&model),
        vec![vec![100.0, 0.0], vec![0.0; 2]],
        vec![1.0; 2],
        transport.clone(),
    );
    let (mut accepted, _) = advance(&initial, 0.25);
    let saved = accepted.clone();
    let map = Arc::clone(accepted.prepared().unwrap());
    let fluxes = vec![zero_flux(&model); 2];
    let feed = vec![vec![0.0; 2]; 2];
    let mut wrong = transport.clone();
    wrong.transfers[0].to_compartment = 0;
    assert!(
        accepted
            .advance_with_inputs(0.25, &fluxes, &[500.0; 2], &wrong, &feed, budget(), || Ok(
                ()
            ))
            .is_err()
    );
    assert!(
        accepted
            .advance_with_inputs(
                0.25,
                &fluxes,
                &[100.0; 2],
                &transport,
                &feed,
                budget(),
                || Ok(())
            )
            .is_err()
    );
    assert!(
        accepted
            .advance_with_inputs(
                0.25,
                &fluxes,
                &[500.0; 2],
                &transport,
                &feed[..1],
                budget(),
                || Ok(())
            )
            .is_err()
    );
    let low = BurnupTransportBudget {
        nuclei_tolerance: 0.0,
        ..budget()
    };
    assert!(
        accepted
            .advance_with_inputs(
                0.25,
                &fluxes,
                &[500.0; 2],
                &transport,
                &feed,
                low,
                || Ok(())
            )
            .is_err()
    );
    let work = BurnupTransportBudget {
        limits: BurnupTransportLimits {
            max_rate_coefficients: 0,
            ..limits()
        },
        ..budget()
    };
    assert!(
        accepted
            .advance(0.25, &fluxes, &[500.0; 2], work, || Ok(()))
            .is_err()
    );
    let mut calls = 0;
    assert_eq!(
        accepted
            .advance(0.25, &fluxes, &[500.0; 2], budget(), || {
                calls += 1;
                if calls == 20 {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .unwrap_err(),
        BurnupError::Fission(FissionError::Cancelled)
    );
    assert_eq!(accepted, saved);
    assert!(Arc::ptr_eq(&map, accepted.prepared().unwrap()));
    let too_many = BurnupTransportLimits {
        burnup: BurnupLimits {
            propagation: fs_la::decay::DecayLimits {
                maximum_nodes: 3,
                ..limits().burnup.propagation
            },
            ..limits().burnup
        },
        ..limits()
    };
    assert!(
        BurnupTransportState::new(
            vec![BurnupInventory::new(Arc::clone(&model), vec![1.0, 0.0]).unwrap(); 2],
            vec![1.0; 2],
            transport.clone(),
            too_many,
            || Ok(())
        )
        .is_err()
    );
    let same_data_different_owner = Arc::new((*model).clone());
    assert!(
        BurnupTransportState::new(
            vec![
                BurnupInventory::new(Arc::clone(&model), vec![1.0, 0.0]).unwrap(),
                BurnupInventory::new(same_data_different_owner, vec![0.0, 0.0]).unwrap(),
            ],
            vec![1.0; 2],
            transport,
            limits(),
            || Ok(())
        )
        .is_err()
    );
}
