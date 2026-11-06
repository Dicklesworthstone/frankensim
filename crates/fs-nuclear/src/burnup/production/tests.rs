//! Manufactured source/reaction oracles, not calibrated isotope data.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{BurnupData, BurnupNuclide, BurnupProduct, BurnupRate, BurnupReaction};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}
fn budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}
fn grid() -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("manufactured-feed-grid/v1".into(), vec![0.0, 1.0]).unwrap())
}
fn reaction(
    parent: usize,
    rate: BurnupRate,
    products: &[(usize, f64)],
    q: EnergyChannels,
) -> BurnupReaction {
    BurnupReaction {
        data_id: format!("manufactured-feed-reaction-{parent}/v1"),
        parent,
        rate,
        products: products
            .iter()
            .map(|&(nuclide, yield_per_event)| BurnupProduct {
                nuclide,
                yield_per_event,
            })
            .collect(),
        energy_per_event: q,
    }
}
fn fixture(
    lambdas: &[f64],
    reactions: Vec<BurnupReaction>,
    populations: Vec<f64>,
) -> (BurnupInventory, GroupFlux) {
    let groups = grid();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-production-network/v1".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 250.0,
                temperature_max_k: 1000.0,
                nuclides: lambdas
                    .iter()
                    .enumerate()
                    .map(|(i, &lambda)| BurnupNuclide {
                        id: NuclideId::new(10, 25 + i as u16, 0).unwrap(),
                        law: DecayLawV1::from_decay_constant(Frequency::new(lambda)).unwrap(),
                    })
                    .collect(),
                reactions,
            },
            budget().limits,
        )
        .unwrap(),
    );
    (
        BurnupInventory::new(model, populations).unwrap(),
        GroupFlux::new(groups, vec![1.0]).unwrap(),
    )
}
fn sink(population: f64) -> (BurnupInventory, GroupFlux) {
    fixture(
        &[0.5],
        vec![reaction(
            0,
            BurnupRate::Decay {
                branch_fraction: 1.0,
            },
            &[],
            energy(2.0, 1.0, 0.5),
        )],
        vec![population],
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn close_values(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        close(a, b);
    }
}
fn unchanged(actual: &BurnupInventory, saved: &BurnupInventory) {
    assert_eq!(actual, saved);
    assert_eq!(
        actual.prepared.as_ref().map(Arc::as_ptr),
        saved.prepared.as_ref().map(Arc::as_ptr)
    );
    assert_eq!(
        actual.production_prepared.as_ref().map(Arc::as_ptr),
        saved.production_prepared.as_ref().map(Arc::as_ptr)
    );
}

#[test]
fn g1_supplied_radioactivity_reaches_saturation_and_keeps_decaying_after_shutdown() {
    let (mut inventory, flux) = sink(3.0);
    let report = inventory
        .advance_with_production(4.0, &flux, 300.0, &[12.0], budget(), || Ok(()))
        .unwrap();
    let remaining = 24.0 + (3.0 - 24.0) * (-2.0_f64).exp();
    let decays = 3.0 + 48.0 - remaining;
    close(inventory.populations()[0], remaining);
    close(report.burnup.decays, decays);
    close_values(&report.burnup.parent_removals, &[decays]);
    assert_eq!(report.produced_nuclei, [48.0]);
    assert_eq!(inventory.cumulative_produced_nuclei(), [48.0]);
    close(report.burnup.released.deposited_j, 2.0 * decays);
    close(report.burnup.released.radiation_j, decays);
    close(report.burnup.released.escaped_j, 0.5 * decays);
    assert!(inventory.production_prepared.is_some());

    let later = inventory
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let late_decays = remaining * -(-1.0_f64).exp_m1();
    close(inventory.populations()[0], remaining * (-1.0_f64).exp());
    close(later.released.deposited_j, 2.0 * late_decays);
    assert!(later.released.deposited_j > 0.0);
    assert_eq!(inventory.cumulative_produced_nuclei(), [48.0]);
    assert!(inventory.production_prepared.is_none());
    assert!(inventory.prepared.is_some());
}

#[test]
fn g1_new_fuel_fissions_and_its_newborn_products_decay_in_the_same_interval() {
    let (mut inventory, flux) = fixture(
        &[0.0, 1.0, 0.0],
        vec![
            reaction(
                0,
                BurnupRate::Fission {
                    cross_sections_m2: vec![1.0],
                },
                &[(1, 2.0)],
                energy(3.0, 2.0, 1.0),
            ),
            reaction(
                1,
                BurnupRate::Decay {
                    branch_fraction: 1.0,
                },
                &[(2, 1.0)],
                energy(0.5, 1.0, 2.0),
            ),
        ],
        vec![0.0; 3],
    );
    let report = inventory
        .advance_with_production(1.0, &flux, 300.0, &[10.0, 0.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let survival = (-1.0_f64).exp();
    let fuel = 10.0 * (1.0 - survival);
    let fissions = 10.0 - fuel;
    // Repeated-root convolution: dX/dt=2P-X, P=10(1-exp(-t)).
    let fragments = 20.0 * (1.0 - 2.0 * survival);
    let decays = 2.0 * fissions - fragments;
    close_values(inventory.populations(), &[fuel, fragments, decays]);
    close(report.burnup.fissions, fissions);
    close(report.burnup.decays, decays);
    assert_eq!(report.burnup.captures, 0.0);
    close(inventory.populations().iter().sum(), 10.0 + fissions);
    close(
        report.burnup.released.deposited_j,
        3.0 * fissions + 0.5 * decays,
    );
    close(report.burnup.released.radiation_j, 2.0 * fissions + decays);
    close(report.burnup.released.escaped_j, fissions + 2.0 * decays);
    assert_eq!(report.produced_nuclei, [10.0, 0.0, 0.0]);
}

#[test]
fn g1_feed_retains_cyclic_visits_and_self_return_events_without_invented_losses() {
    let capture = || BurnupRate::Capture {
        cross_sections_m2: vec![1.0],
    };
    let (mut cycle, flux) = fixture(
        &[0.0, 0.0],
        vec![
            reaction(0, capture(), &[(1, 1.0)], energy(1.0, 0.0, 0.0)),
            reaction(1, capture(), &[(0, 1.0)], energy(2.0, 0.0, 0.0)),
        ],
        vec![0.0, 0.0],
    );
    let time = 2.0;
    let report = cycle
        .advance_with_production(time, &flux, 300.0, &[6.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let transient = 1.0 - (-2.0_f64 * time).exp();
    let a = 3.0 * time + 1.5 * transient;
    let b = 3.0 * time - 1.5 * transient;
    let a_events = 1.5 * time * time + 1.5 * (time - 0.5 * transient);
    let b_events = 1.5 * time * time - 1.5 * (time - 0.5 * transient);
    close_values(cycle.populations(), &[a, b]);
    close_values(&report.burnup.reaction_counts, &[a_events, b_events]);
    close(
        report.burnup.released.deposited_j,
        a_events + 2.0 * b_events,
    );
    close(a + b, 6.0 * time);

    let (mut returning, flux) = fixture(
        &[0.0],
        vec![reaction(
            0,
            BurnupRate::Capture {
                cross_sections_m2: vec![2.0],
            },
            &[(0, 1.0)],
            energy(1.0, 0.0, 0.0),
        )],
        vec![1.0],
    );
    let report = returning
        .advance_with_production(2.0, &flux, 300.0, &[3.0], budget(), || Ok(()))
        .unwrap();
    close(returning.populations()[0], 7.0);
    close(report.burnup.captures, 2.0 * (1.0 * 2.0 + 0.5 * 3.0 * 4.0));
    close(report.burnup.released.deposited_j, 16.0);
}

#[test]
fn g3_partition_restart_and_changed_feed_reuse_the_same_reaction_map() {
    let (initial, flux) = sink(2.0);
    let (whole, _) = initial
        .trial_with_production(1.0, &flux, 300.0, &[5.0], budget(), || Ok(()))
        .unwrap();
    let (mut split, _) = initial
        .trial_with_production(0.5, &flux, 300.0, &[5.0], budget(), || Ok(()))
        .unwrap();
    let cache = Arc::clone(split.production_prepared.as_ref().unwrap());
    let mut restart = split.clone();
    let report = split
        .advance_with_production(0.5, &flux, 400.0, &[5.0], budget(), || Ok(()))
        .unwrap();
    let replay = restart
        .advance_with_production(0.5, &flux, 400.0, &[5.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(report, replay);
    unchanged(&split, &restart);
    close_values(split.populations(), whole.populations());
    close_values(
        &split.cumulative_reaction_counts(),
        &whole.cumulative_reaction_counts(),
    );
    assert_eq!(
        split.cumulative_produced_nuclei(),
        whole.cumulative_produced_nuclei()
    );
    close(
        split.cumulative_released().total_j(),
        whole.cumulative_released().total_j(),
    );
    assert!(Arc::ptr_eq(
        &cache,
        split.production_prepared.as_ref().unwrap()
    ));
    split
        .advance_with_production(0.5, &flux, 500.0, &[7.0], budget(), || Ok(()))
        .unwrap();
    assert!(Arc::ptr_eq(
        &cache,
        split.production_prepared.as_ref().unwrap()
    ));
    assert_eq!(split.cumulative_produced_nuclei(), [8.5]);
    assert_eq!(cache.temperature_k(), 300.0);
    assert!(Arc::ptr_eq(cache.model(), initial.model()));

    let prepared = initial
        .model()
        .prepare_production(0.5, &flux, 300.0, budget().limits, || Ok(()))
        .unwrap();
    let (mut explicit, expected) = initial
        .trial_production_prepared(&prepared, &[5.0], budget(), || Ok(()))
        .unwrap();
    let mut other = initial.clone();
    assert_eq!(
        other
            .advance_production_prepared(&prepared, &[5.0], budget(), || Ok(()))
            .unwrap(),
        expected
    );
    assert_eq!(explicit, other);
    explicit
        .advance(0.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert!(explicit.prepared.is_some());
    explicit
        .advance_with_production(0.5, &flux, 300.0, &[0.0], budget(), || Ok(()))
        .unwrap();
    assert!(explicit.prepared.is_none());
    assert!(explicit.production_prepared.is_some());
}

#[test]
fn g0_production_cache_uses_all_temperature_dependent_channel_rates() {
    let table = |reaction, values: Vec<f64>| {
        Arc::new(
            CrossSectionTable::new(
                CrossSectionData {
                    data_id: format!("manufactured-feed-{reaction:?}/v1"),
                    nuclide: NuclideId::new(10, 25, 0).unwrap(),
                    reaction,
                    groups: grid(),
                    temperatures_k: vec![250.0, 400.0, 1000.0],
                    cross_sections_m2: values.into_iter().map(|v| vec![v]).collect(),
                },
                CrossSectionLimits::default(),
            )
            .unwrap(),
        )
    };
    let (mut inventory, flux) = fixture(
        &[0.0, 0.0, 0.0],
        vec![
            reaction(
                0,
                BurnupRate::TabulatedFission {
                    cross_sections: table(NeutronReaction::Fission, vec![0.25, 0.75, 0.5]),
                },
                &[(1, 2.0)],
                energy(3.0, 0.0, 0.0),
            ),
            reaction(
                0,
                BurnupRate::TabulatedCapture {
                    cross_sections: table(NeutronReaction::Capture, vec![0.75, 0.25, 1.5]),
                },
                &[(2, 1.0)],
                energy(1.0, 0.0, 0.0),
            ),
        ],
        vec![1.0, 0.0, 0.0],
    );
    let feed = [2.0, 0.0, 0.0];
    inventory
        .advance_with_production(0.5, &flux, 250.0, &feed, budget(), || Ok(()))
        .unwrap();
    let first = Arc::clone(inventory.production_prepared.as_ref().unwrap());
    assert_eq!(first.reaction_rate_constants(), [0.25, 0.75]);
    let changed = inventory
        .advance_with_production(0.5, &flux, 400.0, &feed, budget(), || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(
        &first,
        inventory.production_prepared.as_ref().unwrap()
    ));
    close(changed.burnup.fissions / changed.burnup.captures, 3.0);
    inventory
        .advance_with_production(0.5, &flux, 250.0, &feed, budget(), || Ok(()))
        .unwrap();
    let equivalent = Arc::clone(inventory.production_prepared.as_ref().unwrap());
    let half_flux = GroupFlux::new(Arc::clone(flux.groups()), vec![0.5]).unwrap();
    inventory
        .advance_with_production(0.5, &half_flux, 1000.0, &[3.0, 0.0, 0.0], budget(), || {
            Ok(())
        })
        .unwrap();
    assert!(Arc::ptr_eq(
        &equivalent,
        inventory.production_prepared.as_ref().unwrap()
    ));
    assert_eq!(equivalent.flux(), &flux);
    let before = inventory.clone();
    assert!(
        inventory
            .advance_with_production(0.5, &half_flux, 1001.0, &feed, budget(), || Ok(()))
            .is_err()
    );
    unchanged(&inventory, &before);
}

#[test]
fn g0_production_schedules_and_preparations_keep_exact_model_and_work_binding() {
    let (inventory, flux) = sink(1.0);
    let change = |time_s, rate| BurnupProductionChange {
        time_s,
        rates_nuclei_s: vec![rate],
    };
    let schedule = BurnupProductionSchedule::new(
        Arc::clone(inventory.model()),
        vec![3.0],
        vec![change(0.5, 0.0), change(1.0, 7.0)],
        2,
    )
    .unwrap();
    assert!(Arc::ptr_eq(schedule.model(), inventory.model()));
    assert_eq!(schedule.changes().len(), 2);
    assert_eq!(schedule.rates_at(0.5_f64.next_down()).unwrap(), [3.0]);
    assert_eq!(schedule.rates_at(0.5).unwrap(), [0.0]);
    assert_eq!(schedule.rates_at(1.0).unwrap(), [7.0]);
    assert_eq!(schedule.next_after(0.0), Some(0.5));
    assert_eq!(schedule.next_after(0.5), Some(1.0));
    assert_eq!(schedule.next_after(1.0), None);
    assert!(schedule.rates_at(f64::INFINITY).is_err());
    assert_eq!(schedule.next_after(f64::NAN), None);
    assert!(
        BurnupProductionSchedule::new(
            Arc::clone(inventory.model()),
            vec![0.0],
            vec![change(1.0, 1.0)],
            0
        )
        .is_err()
    );
    assert!(
        BurnupProductionSchedule::new(
            Arc::clone(inventory.model()),
            vec![0.0],
            vec![change(1.0, 1.0), change(1.0, 0.0)],
            2
        )
        .is_err()
    );
    assert!(
        BurnupProductionSchedule::new(Arc::clone(inventory.model()), vec![-1.0], vec![], 0)
            .is_err()
    );

    let prepared = inventory
        .model()
        .prepare_production(0.5, &flux, 300.0, budget().limits, || Ok(()))
        .unwrap();
    let other_model = Arc::new((**inventory.model()).clone());
    let other = BurnupInventory::new(other_model, vec![1.0]).unwrap();
    assert!(
        other
            .trial_production_prepared(&prepared, &[1.0], budget(), || Ok(()))
            .is_err()
    );
    let mut changed_limits = budget();
    changed_limits.limits.propagation.maximum_work -= 1;
    assert!(
        inventory
            .trial_production_prepared(&prepared, &[1.0], changed_limits, || Ok(()))
            .is_err()
    );
}

#[test]
fn g4_production_budget_range_and_checkpoint_refusals_preserve_histories_and_cache() {
    let (mut inventory, flux) = sink(2.1);
    inventory
        .advance_with_production(0.5, &flux, 300.0, &[3.2], budget(), || Ok(()))
        .unwrap();
    let saved = inventory.clone();
    let mut polls = 0;
    let (_, report) = inventory
        .trial_with_production(0.6, &flux, 300.0, &[3.2], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            inventory.advance_with_production(0.6, &flux, 300.0, &[3.2], budget(), || {
                count += 1;
                if count == stop {
                    Err(FissionError::Budget("original source checkpoint"))
                } else {
                    Ok(())
                }
            }),
            Err(BurnupError::Fission(FissionError::Budget(
                "original source checkpoint"
            )))
        );
        unchanged(&inventory, &saved);
    }
    let mut exhausted = budget();
    exhausted.limits.propagation.maximum_work = 1;
    assert!(
        inventory
            .advance_with_production(0.6, &flux, 300.0, &[3.2], exhausted, || Ok(()))
            .is_err()
    );
    assert!(report.burnup.population_roundoff_nuclei > 0.0);
    let mut strict = budget();
    strict.nuclei_tolerance = 0.0;
    assert!(
        inventory
            .advance_with_production(0.6, &flux, 300.0, &[3.2], strict, || Ok(()))
            .is_err()
    );
    for rates in [vec![], vec![-1.0], vec![f64::NAN], vec![f64::INFINITY]] {
        assert!(
            inventory
                .advance_with_production(0.6, &flux, 300.0, &rates, budget(), || Ok(()))
                .is_err()
        );
    }
    let wrong_groups =
        Arc::new(NeutronGroups::new("different-feed-grid".into(), vec![0.0, 1.0]).unwrap());
    let wrong_flux = GroupFlux::new(wrong_groups, vec![0.0]).unwrap();
    assert!(
        inventory
            .advance_with_production(0.6, &wrong_flux, 300.0, &[0.0], budget(), || Ok(()))
            .is_err()
    );
    unchanged(&inventory, &saved);
}

#[test]
fn g0_finite_produced_counts_are_admitted_even_when_rate_sum_overflows() {
    let (mut stable, flux) = fixture(&[0.0, 0.0], vec![], vec![0.0, 0.0]);
    let report = stable
        .advance_with_production(0.125, &flux, 300.0, &[f64::MAX, f64::MAX], budget(), || {
            Ok(())
        })
        .unwrap();
    let count = f64::MAX * 0.125;
    assert_eq!(stable.populations(), [count, count]);
    assert_eq!(report.produced_nuclei, [count, count]);
    assert_eq!(report.burnup.released, EnergyChannels::default());
    let (mut tiny, flux) = fixture(&[0.0], vec![], vec![0.0]);
    let saved = tiny.clone();
    assert!(matches!(
        tiny.advance_with_production(
            f64::MIN_POSITIVE,
            &flux,
            300.0,
            &[f64::MIN_POSITIVE],
            budget(),
            || Ok(())
        ),
        Err(BurnupError::Arithmetic(_))
    ));
    unchanged(&tiny, &saved);
    let zero = tiny
        .advance_with_production(0.0, &flux, 300.0, &[f64::MAX], budget(), || Ok(()))
        .unwrap();
    assert_eq!(zero.produced_nuclei, [0.0]);
    assert_eq!(tiny, saved);
}
