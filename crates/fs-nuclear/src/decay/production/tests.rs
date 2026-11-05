//! Physical production/decay oracles, original-order accounts and atomic refusal.

use super::{DecayProductionChange, DecayProductionSchedule};
use crate::EnergyChannels;
use crate::decay::{
    DecayBranch, DecayBudget, DecayError, DecayInventory, DecayModel, DecayNuclide,
};
use fs_la::decay::DecayLimits;
use fs_qty::{Frequency, Time};
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn channels(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn node(id: &str, rate: f64) -> DecayNuclide {
    DecayNuclide {
        id: id.into(),
        law: if rate == 0.0 {
            DecayLawV1::stable()
        } else {
            DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap()
        },
    }
}

fn branch(parent: usize, daughter: usize, probability: f64, energy: EnergyChannels) -> DecayBranch {
    DecayBranch {
        parent,
        daughter,
        probability,
        energy_per_decay: energy,
    }
}

fn model(nodes: Vec<DecayNuclide>, branches: Vec<DecayBranch>) -> Arc<DecayModel> {
    Arc::new(
        DecayModel::new(
            "synthetic-production-oracle/revision-1".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn single(rate: f64, energy: EnergyChannels) -> Arc<DecayModel> {
    model(
        vec![node("parent", rate), node("stable", 0.0)],
        vec![branch(0, 1, 1.0, energy)],
    )
}

fn budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}

fn close(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!(
            (actual - expected).abs() <= 3e-11 * expected.abs(),
            "actual {actual:.17e}, expected {expected:.17e}"
        );
    }
}

fn close_values(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        close(actual, expected);
    }
}

fn close_channels(actual: EnergyChannels, expected: EnergyChannels) {
    close(actual.deposited_j, expected.deposited_j);
    close(actual.radiation_j, expected.radiation_j);
    close(actual.escaped_j, expected.escaped_j);
}

#[test]
fn g1_constant_production_builds_toward_equilibrium_and_shutdown_keeps_decay_active() {
    let rate = 0.5;
    let production = 12.0;
    let model = single(rate, channels(2.0, 1.0, 0.5));
    let initial = DecayInventory::new(model, vec![3.0, 2.0]).unwrap();
    let duration = 4.0;
    let (mut next, report) = initial
        .trial_with_production(duration, &[production, 0.0], budget(), || Ok(()))
        .unwrap();

    // Integrating N' = s-lambda*N gives this independent constant-source oracle.
    let equilibrium = production / rate;
    let parent = equilibrium + (3.0 - equilibrium) * (-rate * duration).exp();
    let decays = 3.0 + production * duration - parent;
    close_values(next.populations(), &[parent, 2.0 + decays]);
    close_values(&report.produced_nuclei, &[48.0, 0.0]);
    close_values(&report.decay.decays, &[decays, 0.0]);
    close_values(&report.decay.branch_decays, &[decays]);
    close_channels(
        report.decay.released,
        channels(2.0 * decays, decays, 0.5 * decays),
    );
    assert!(next.populations()[0] > 3.0 && next.populations()[0] < equilibrium);
    assert_eq!(next.cumulative_produced_nuclei(), report.produced_nuclei);
    assert_eq!(initial.populations(), &[3.0, 2.0]);
    assert_eq!(initial.cumulative_produced_nuclei(), [0.0, 0.0]);

    // The ordinary shutdown path retains the production history and releases
    // half the remaining parent inventory over its original physical half-life.
    let births = next.cumulative_produced_nuclei();
    let shutdown = next
        .advance(
            Time::new(core::f64::consts::LN_2 / rate).value(),
            budget(),
            || Ok(()),
        )
        .unwrap();
    close_values(
        next.populations(),
        &[0.5 * parent, 2.0 + decays + 0.5 * parent],
    );
    close_values(&shutdown.decays, &[0.5 * parent, 0.0]);
    close_channels(
        shutdown.released,
        channels(parent, 0.5 * parent, 0.25 * parent),
    );
    assert_eq!(next.cumulative_produced_nuclei(), births);
    close(next.cumulative_decays()[0], decays + 0.5 * parent);
}

#[test]
fn g1_equal_rate_production_credits_daughter_decays_and_energy_inside_the_interval() {
    let rate = 1.25;
    let production = 9.0;
    let duration = 0.8;
    let model = model(
        vec![
            node("parent", rate),
            node("daughter", rate),
            node("stable", 0.0),
        ],
        vec![
            branch(0, 1, 1.0, channels(2.0, 3.0, 5.0)),
            branch(1, 2, 1.0, channels(7.0, 11.0, 13.0)),
        ],
    );
    let initial = DecayInventory::new(model, vec![0.0; 3]).unwrap();
    let (next, report) = initial
        .trial_with_production(duration, &[production, 0.0, 0.0], budget(), || Ok(()))
        .unwrap();

    // Closed forms obtained by integrating the age distribution of uniform
    // production through P -> D -> S, with lambda_P == lambda_D.
    let x = rate * duration;
    let surviving = (-x).exp();
    let scale = production / rate;
    let parent = scale * (1.0 - surviving);
    let daughter = scale * (1.0 - (1.0 + x) * surviving);
    let parent_decays = scale * (x - 1.0 + surviving);
    let daughter_decays = scale * (x - 2.0 + (2.0 + x) * surviving);
    assert!(daughter_decays > 0.0);
    close_values(next.populations(), &[parent, daughter, daughter_decays]);
    close_values(&report.decay.decays, &[parent_decays, daughter_decays, 0.0]);
    close_values(
        &report.decay.branch_decays,
        &[parent_decays, daughter_decays],
    );
    close_values(&report.produced_nuclei, &[production * duration, 0.0, 0.0]);
    close_channels(
        report.decay.released,
        channels(
            2.0 * parent_decays + 7.0 * daughter_decays,
            3.0 * parent_decays + 11.0 * daughter_decays,
            5.0 * parent_decays + 13.0 * daughter_decays,
        ),
    );
    close_channels(next.cumulative_released(), report.decay.released);
    close(next.populations().iter().sum(), production * duration);
}

#[test]
fn g0_stable_feed_and_zero_duration_preserve_the_complete_history() {
    let model = model(vec![node("stable-a", 0.0), node("stable-b", 0.0)], vec![]);
    let mut inventory = DecayInventory::new(Arc::clone(&model), vec![3.0, 7.0]).unwrap();
    let report = inventory
        .advance_with_production(4.0, &[5.0, 2.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(inventory.populations(), &[23.0, 15.0]);
    assert_eq!(report.produced_nuclei, [20.0, 8.0]);
    assert_eq!(inventory.cumulative_produced_nuclei(), [20.0, 8.0]);
    assert_eq!(report.decay.decays, [0.0, 0.0]);
    assert_eq!(report.decay.released, EnergyChannels::default());
    assert_eq!(inventory.activity_bq().unwrap(), [0.0, 0.0]);

    let saved = inventory.clone();
    let idle = inventory
        .advance_with_production(0.0, &[f64::MAX, f64::MAX], budget(), || Ok(()))
        .unwrap();
    assert_eq!(inventory, saved);
    assert_eq!(idle.produced_nuclei, [0.0, 0.0]);
    assert_eq!(idle.decay.decays, [0.0, 0.0]);
    assert_eq!(idle.decay.released, EnergyChannels::default());

    // Individual rates can have an overflowing sum while the physical counts
    // produced over a short admitted interval are all finite and representable.
    let fresh = DecayInventory::new(model, vec![0.0, 0.0]).unwrap();
    let (fed, report) = fresh
        .trial_with_production(
            f64::MIN_POSITIVE,
            &[f64::MAX, f64::MAX],
            budget(),
            || Ok(()),
        )
        .unwrap();
    let produced = f64::MAX * f64::MIN_POSITIVE;
    assert_eq!(fed.populations(), &[produced, produced]);
    assert_eq!(report.produced_nuclei, [produced, produced]);
    assert_eq!(report.decay.released, EnergyChannels::default());
}

#[test]
fn g3_production_partition_and_restart_preserve_birth_decay_and_emission_accounts() {
    let model = model(
        vec![
            node("parent", 0.3),
            node("daughter", 0.7),
            node("stable", 0.0),
        ],
        vec![
            branch(0, 1, 0.6, channels(1.0, 2.0, 0.0)),
            branch(0, 2, 0.4, channels(3.0, 0.0, 4.0)),
            branch(1, 2, 1.0, channels(5.0, 7.0, 11.0)),
        ],
    );
    let rates = [4.0, 0.75, 0.25];
    let initial = DecayInventory::new(Arc::clone(&model), vec![80.0, 2.0, 5.0]).unwrap();
    let (whole, _) = initial
        .trial_with_production(1.8, &rates, budget(), || Ok(()))
        .unwrap();
    let (mut split, first) = initial
        .trial_with_production(0.7, &rates, budget(), || Ok(()))
        .unwrap();
    let checkpoint = split.clone();
    let mut restarted = checkpoint.clone();
    let prepared = model
        .prepare_production(1.1, budget().propagation, || Ok(()))
        .unwrap();
    let second = split
        .advance_production_prepared(&prepared, &rates, budget(), || Ok(()))
        .unwrap();
    let replay = restarted
        .advance_production_prepared(&prepared, &rates, budget(), || Ok(()))
        .unwrap();
    assert_eq!(split, restarted);
    assert_eq!(second, replay);
    assert_eq!(
        checkpoint.cumulative_produced_nuclei(),
        first.produced_nuclei
    );
    close_values(split.populations(), whole.populations());
    close_values(&split.cumulative_decays(), &whole.cumulative_decays());
    close_values(
        &split.cumulative_produced_nuclei(),
        &whole.cumulative_produced_nuclei(),
    );
    close_channels(split.cumulative_released(), whole.cumulative_released());
    for (index, &rate) in rates.iter().enumerate() {
        close(split.cumulative_produced_nuclei()[index], 1.8 * rate);
        close(
            first.produced_nuclei[index] + second.produced_nuclei[index],
            1.8 * rate,
        );
    }
}

#[test]
fn g0_production_schedule_and_rate_admission_preserve_original_order_and_exact_transitions() {
    let model = single(1.0, channels(2.0, 1.0, 0.0));
    let changes = vec![
        DecayProductionChange {
            time_s: 2.0,
            rates_nuclei_s: vec![1.0, 3.0],
        },
        DecayProductionChange {
            time_s: 5.0,
            rates_nuclei_s: vec![0.0, 0.0],
        },
    ];
    let schedule =
        DecayProductionSchedule::new(Arc::clone(&model), vec![4.0, 0.0], changes.clone(), 2)
            .unwrap();
    assert!(Arc::ptr_eq(schedule.model(), &model));
    assert_eq!(schedule.rates_at(0.0).unwrap(), &[4.0, 0.0]);
    assert_eq!(
        schedule
            .rates_at(f64::from_bits(2.0_f64.to_bits() - 1))
            .unwrap(),
        &[4.0, 0.0]
    );
    assert_eq!(schedule.rates_at(2.0).unwrap(), &[1.0, 3.0]);
    assert_eq!(schedule.rates_at(5.0).unwrap(), &[0.0, 0.0]);
    assert_eq!(schedule.next_after(0.0), Some(2.0));
    assert_eq!(schedule.next_after(2.0), Some(5.0));
    assert_eq!(schedule.next_after(5.0), None);
    assert!(matches!(
        DecayProductionSchedule::new(Arc::clone(&model), vec![4.0, 0.0], changes, 1),
        Err(DecayError::Budget(_)),
    ));

    let mut inventory = DecayInventory::new(Arc::clone(&model), vec![7.0, 2.0]).unwrap();
    let saved = inventory.clone();
    for rates in [
        vec![1.0],
        vec![-1.0, 0.0],
        vec![0.0, f64::NAN],
        vec![f64::INFINITY, 0.0],
    ] {
        assert!(
            DecayProductionSchedule::new(Arc::clone(&model), rates.clone(), vec![], 0).is_err()
        );
        assert!(
            DecayProductionSchedule::new(
                Arc::clone(&model),
                vec![0.0, 0.0],
                vec![DecayProductionChange {
                    time_s: 1.0,
                    rates_nuclei_s: rates.clone()
                }],
                1,
            )
            .is_err()
        );
        assert!(
            inventory
                .advance_with_production(1.0, &rates, budget(), || Ok(()))
                .is_err()
        );
        assert_eq!(inventory, saved);
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(schedule.rates_at(bad).is_err());
        assert_eq!(schedule.next_after(bad), None);
        assert!(
            inventory
                .advance_with_production(bad, &[1.0, 0.0], budget(), || Ok(()))
                .is_err()
        );
        assert_eq!(inventory, saved);
    }
    for times in [
        vec![0.0],
        vec![-1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![2.0, 2.0],
        vec![2.0, 1.0],
    ] {
        let changes = times
            .into_iter()
            .map(|time_s| DecayProductionChange {
                time_s,
                rates_nuclei_s: vec![0.0, 0.0],
            })
            .collect();
        assert!(
            DecayProductionSchedule::new(Arc::clone(&model), vec![1.0, 0.0], changes, 2).is_err()
        );
    }
}

#[test]
fn g4_birth_credits_charge_per_node_and_product_roundoff_before_publication() {
    let stable = model(vec![node("stable-a", 0.0), node("stable-b", 0.0)], vec![]);
    let initial = DecayInventory::new(stable, vec![1e16, 1e16 + 2.0]).unwrap();
    let exact = DecayBudget {
        nuclei_tolerance: 0.0,
        ..budget()
    };
    let mut refused = initial.clone();
    assert!(matches!(
        refused.advance_with_production(1.0, &[1.0, 1.0], exact, || Ok(())),
        Err(DecayError::Budget(_)),
    ));
    assert_eq!(refused, initial);
    let (next, report) = initial
        .trial_with_production(
            1.0,
            &[1.0, 1.0],
            DecayBudget {
                nuclei_tolerance: 16.0,
                ..budget()
            },
            || Ok(()),
        )
        .unwrap();
    // One birth rounds down and the other rounds up. Their opposite local
    // defects cannot cancel to satisfy a global population allowance.
    assert_eq!(next.populations(), &[1e16, 1e16 + 4.0]);
    assert!(report.decay.population_roundoff_nuclei >= 2.0);
    assert_eq!(next.cumulative_produced_nuclei(), [1.0, 1.0]);

    let one_stable = model(vec![node("stable", 0.0)], vec![]);
    let mut product_roundoff = DecayInventory::new(one_stable, vec![0.0]).unwrap();
    let saved = product_roundoff.clone();
    assert!(matches!(
        product_roundoff.advance_with_production(0.1, &[0.1], exact, || Ok(())),
        Err(DecayError::Budget(_)),
    ));
    assert_eq!(product_roundoff, saved);
    let admitted = product_roundoff
        .advance_with_production(0.1, &[0.1], budget(), || Ok(()))
        .unwrap();
    assert!(admitted.decay.population_roundoff_nuclei > 0.0);
    assert_eq!(admitted.decay.released, EnergyChannels::default());
    let saved = product_roundoff.clone();
    assert_eq!(
        product_roundoff.advance_with_production(2.0, &[f64::MAX], budget(), || Ok(())),
        Err(DecayError::Arithmetic),
    );
    assert_eq!(product_roundoff, saved);
}

#[test]
fn g4_early_midway_and_publication_cancellation_preserve_production_history_and_plan() {
    let model = single(1.0, channels(2.0, 3.0, 5.0));
    assert!(matches!(
        model.prepare_production(1.0, budget().propagation, || Err(DecayError::Budget(
            "caller checkpoint"
        ))),
        Err(DecayError::Budget("caller checkpoint")),
    ));
    let prepared = model
        .prepare_production(1.0, budget().propagation, || Ok(()))
        .unwrap();
    let initial = DecayInventory::new(model, vec![7.0, 2.0]).unwrap();
    let rates = [5.0, 1.0];
    let mut polls = 0;
    let expected = initial
        .trial_production_prepared(&prepared, &rates, budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(polls > 3);
    for stop in [1, polls / 2, polls] {
        let mut candidate = initial.clone();
        let mut calls = 0;
        let result = candidate.advance_production_prepared(&prepared, &rates, budget(), || {
            calls += 1;
            if calls == stop {
                Err(DecayError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err(DecayError::Cancelled));
        assert_eq!(candidate, initial);
        let retried = initial
            .trial_production_prepared(&prepared, &rates, budget(), || Ok(()))
            .unwrap();
        assert_eq!(retried, expected);
    }
}

#[test]
fn g0_prepared_production_reuses_distinct_rate_vectors_and_rejects_foreign_bindings() {
    let model = single(0.5, channels(2.0, 3.0, 4.0));
    let prepared = model
        .prepare_production(2.0, budget().propagation, || Ok(()))
        .unwrap();
    assert!(Arc::ptr_eq(prepared.model(), &model));
    assert_eq!(prepared.duration_s(), 2.0);
    let initial = DecayInventory::new(Arc::clone(&model), vec![0.0, 0.0]).unwrap();
    let (parent_fed, parent_report) = initial
        .trial_production_prepared(&prepared, &[4.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let (stable_fed, stable_report) = initial
        .trial_production_prepared(&prepared, &[0.0, 4.0], budget(), || Ok(()))
        .unwrap();
    let parent = 8.0 * (1.0 - (-1.0_f64).exp());
    close_values(parent_fed.populations(), &[parent, 8.0 - parent]);
    close_channels(
        parent_report.decay.released,
        channels(
            2.0 * (8.0 - parent),
            3.0 * (8.0 - parent),
            4.0 * (8.0 - parent),
        ),
    );
    assert_eq!(stable_fed.populations(), &[0.0, 8.0]);
    assert_eq!(stable_report.decay.released, EnergyChannels::default());
    assert_eq!(parent_report.produced_nuclei, [8.0, 0.0]);
    assert_eq!(stable_report.produced_nuclei, [0.0, 8.0]);

    let mut foreign = DecayInventory::new(Arc::new((*model).clone()), vec![0.0, 0.0]).unwrap();
    let saved = foreign.clone();
    assert!(matches!(
        foreign.advance_production_prepared(&prepared, &[4.0, 0.0], budget(), || Ok(())),
        Err(DecayError::Invalid(_)),
    ));
    assert_eq!(foreign, saved);
    let mut limits = budget();
    limits.propagation.maximum_work += 1;
    assert!(matches!(
        initial.trial_production_prepared(&prepared, &[4.0, 0.0], limits, || Ok(())),
        Err(DecayError::Invalid(_)),
    ));
    let mut low_work = initial.clone();
    let mut limit = budget();
    limit.propagation.maximum_work = 1;
    assert!(
        low_work
            .advance_with_production(2.0, &[4.0, 0.0], limit, || Ok(()))
            .is_err()
    );
    assert_eq!(low_work, initial);
    assert_eq!(initial.cumulative_produced_nuclei(), [0.0, 0.0]);
}
