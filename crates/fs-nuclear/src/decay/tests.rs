use super::*;
use fs_qty::{Frequency, Time};

fn channels(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn node(id: &str, rate_per_s: f64) -> DecayNuclide {
    DecayNuclide {
        id: id.into(),
        law: if rate_per_s == 0.0 {
            DecayLawV1::stable()
        } else {
            DecayLawV1::from_decay_constant(Frequency::new(rate_per_s)).unwrap()
        },
    }
}

fn branch(
    parent: usize,
    daughter: usize,
    probability: f64,
    energy_per_decay: EnergyChannels,
) -> DecayBranch {
    DecayBranch {
        parent,
        daughter,
        probability,
        energy_per_decay,
    }
}

fn model(nuclides: Vec<DecayNuclide>, branches: Vec<DecayBranch>) -> Arc<DecayModel> {
    Arc::new(
        DecayModel::new(
            "synthetic-decay-oracle/revision-1".into(),
            nuclides,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn single(half_life: Time, energy: EnergyChannels) -> Arc<DecayModel> {
    let parent = DecayNuclide {
        id: "parent".into(),
        law: DecayLawV1::from_half_life(half_life).unwrap(),
    };
    model(
        vec![parent, node("stable-daughter", 0.0)],
        vec![branch(0, 1, 1.0, energy)],
    )
}

fn budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-9,
        energy_tolerance_j: 1e-9,
    }
}

fn close(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!(
            (actual - expected).abs() <= 3e-12 * expected.abs(),
            "actual {actual:.17e}, expected {expected:.17e}"
        );
    }
}

fn close_channels(actual: EnergyChannels, expected: EnergyChannels) {
    close(actual.deposited_j, expected.deposited_j);
    close(actual.radiation_j, expected.radiation_j);
    close(actual.escaped_j, expected.escaped_j);
}

#[test]
fn g1_equal_rates_credit_daughters_born_and_decayed_inside_the_interval() {
    let model = model(
        vec![
            node("parent", 2.0),
            node("daughter", 2.0),
            node("stable", 0.0),
        ],
        vec![
            branch(0, 1, 1.0, channels(1.0, 2.0, 3.0)),
            branch(1, 2, 1.0, channels(5.0, 7.0, 11.0)),
        ],
    );
    let initial = DecayInventory::new(Arc::clone(&model), vec![100.0, 0.0, 3.0]).unwrap();
    let time = Time::new(0.7);
    let (next, report) = initial.trial(time.value(), budget(), || Ok(())).unwrap();
    let x = 2.0 * time.value();
    let survival = (-x).exp();
    let parent_decays = 100.0 * (1.0 - survival);
    let daughter_decays = 100.0 * (1.0 - (1.0 + x) * survival);
    close(next.populations()[0], 100.0 * survival);
    close(next.populations()[1], 100.0 * x * survival);
    close(next.populations()[2], 3.0 + daughter_decays);
    close(report.decays[0], parent_decays);
    close(report.decays[1], daughter_decays);
    assert_eq!(report.decays[2], 0.0);
    close_channels(
        report.released,
        channels(
            parent_decays + 5.0 * daughter_decays,
            2.0 * parent_decays + 7.0 * daughter_decays,
            3.0 * parent_decays + 11.0 * daughter_decays,
        ),
    );
    close_channels(next.cumulative_released(), report.released);
    assert_eq!(next.cumulative_decays(), report.decays);
    assert_eq!(initial.populations(), &[100.0, 0.0, 3.0]);
    assert_eq!(initial.cumulative_released(), EnergyChannels::default());
}

#[test]
fn g1_duplicate_daughter_channels_keep_their_own_energy_and_original_order() {
    let parent = DecayNuclide {
        id: "parent-last".into(),
        law: DecayLawV1::from_half_life(Time::new(4.0)).unwrap(),
    };
    let branches = vec![
        branch(2, 0, 0.25, channels(4.0, 2.0, 0.0)),
        branch(2, 0, 0.5, channels(0.0, 6.0, 2.0)),
        branch(2, 1, 0.25, channels(8.0, 0.0, 4.0)),
    ];
    let model = model(
        vec![
            node("stable-first", 0.0),
            node("stable-second", 0.0),
            parent,
        ],
        branches.clone(),
    );
    assert_eq!(model.branches(), branches);
    assert_eq!(model.data_id(), "synthetic-decay-oracle/revision-1");
    let state = DecayInventory::new(Arc::clone(&model), vec![3.0, 9.0, 80.0]).unwrap();
    let (next, report) = state
        .trial(Time::new(4.0).value(), budget(), || Ok(()))
        .unwrap();
    for (&actual, expected) in next.populations().iter().zip([33.0, 19.0, 40.0]) {
        close(actual, expected);
    }
    for (&actual, expected) in report.branch_decays.iter().zip([10.0, 20.0, 10.0]) {
        close(actual, expected);
    }
    assert_eq!(report.decays[..2], [0.0, 0.0]);
    close(report.decays[2], 40.0);
    close_channels(report.released, channels(120.0, 140.0, 80.0));
    let activity = next.activity_bq().unwrap();
    assert_eq!(activity[..2], [0.0, 0.0]);
    close(activity[2], 40.0 * core::f64::consts::LN_2 / 4.0);
    assert!(Arc::ptr_eq(next.model(), &model));
}

#[test]
fn g1_activity_and_endpoint_power_are_distinct_from_interval_emission() {
    let model = single(Time::new(2.0), channels(2.0, 1.0, 0.5));
    let state = DecayInventory::new(model, vec![100.0, 0.0]).unwrap();
    let rate = core::f64::consts::LN_2 / 2.0;
    close(state.activity_bq().unwrap()[0], 100.0 * rate);
    let initial_power = state.decay_power_w().unwrap();
    close_channels(
        initial_power,
        channels(200.0 * rate, 100.0 * rate, 50.0 * rate),
    );
    let (next, report) = state
        .trial(Time::new(2.0).value(), budget(), || Ok(()))
        .unwrap();
    let endpoint_power = next.decay_power_w().unwrap();
    close_channels(
        endpoint_power,
        channels(100.0 * rate, 50.0 * rate, 25.0 * rate),
    );
    close_channels(report.released, channels(100.0, 50.0, 25.0));
    assert!(endpoint_power.deposited_j * report.duration_s < report.released.deposited_j);
    assert!(initial_power.deposited_j * report.duration_s > report.released.deposited_j);
    assert_eq!(next.activity_bq().unwrap()[1], 0.0);
}

#[test]
fn g0_prepared_reuse_requires_the_same_immutable_model_and_numerical_limits() {
    let model = single(Time::new(2.0), channels(3.0, 1.0, 2.0));
    let prepared = model
        .prepare(Time::new(1.0).value(), budget().propagation, || Ok(()))
        .unwrap();
    let a = DecayInventory::new(Arc::clone(&model), vec![10.0, 0.0]).unwrap();
    let b = DecayInventory::new(Arc::clone(&model), vec![40.0, 0.0]).unwrap();
    let (next_a, report_a) = a.trial_prepared(&prepared, budget(), || Ok(())).unwrap();
    let (next_b, report_b) = b.trial_prepared(&prepared, budget(), || Ok(())).unwrap();
    for (&small, &large) in next_a.populations().iter().zip(next_b.populations()) {
        close(large, 4.0 * small);
    }
    close(
        report_b.released.deposited_j,
        4.0 * report_a.released.deposited_j,
    );
    assert_eq!(a.populations(), &[10.0, 0.0]);
    assert_eq!(b.populations(), &[40.0, 0.0]);
    let mut foreign = DecayInventory::new(Arc::new((*model).clone()), vec![10.0, 0.0]).unwrap();
    let saved = foreign.clone();
    assert!(matches!(
        foreign.advance_prepared(&prepared, budget(), || Ok(())),
        Err(DecayError::Invalid(_))
    ));
    assert_eq!(foreign, saved);
    let mut different_limits = budget();
    different_limits.propagation.maximum_work += 1;
    assert!(matches!(
        a.trial_prepared(&prepared, different_limits, || Ok(())),
        Err(DecayError::Invalid(_))
    ));
    assert_eq!(prepared.duration_s(), 1.0);
    assert!(Arc::ptr_eq(prepared.model(), &model));
}

#[test]
fn g3_interval_partition_and_restart_preserve_population_and_lifetime_accounts() {
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
    let initial = DecayInventory::new(Arc::clone(&model), vec![80.0, 2.0, 5.0]).unwrap();
    let (whole, _) = initial
        .trial(Time::new(1.8).value(), budget(), || Ok(()))
        .unwrap();
    let (mut split, first) = initial
        .trial(Time::new(0.7).value(), budget(), || Ok(()))
        .unwrap();
    let checkpoint = split.clone();
    let mut restarted = checkpoint.clone();
    let prepared = model
        .prepare(Time::new(1.1).value(), budget().propagation, || Ok(()))
        .unwrap();
    let second = split
        .advance_prepared(&prepared, budget(), || Ok(()))
        .unwrap();
    let replay = restarted
        .advance_prepared(&prepared, budget(), || Ok(()))
        .unwrap();
    assert_eq!(replay, second);
    assert_eq!(restarted, split);
    assert_eq!(checkpoint.cumulative_decays(), first.decays);
    for (&actual, &expected) in split.populations().iter().zip(whole.populations()) {
        close(actual, expected);
    }
    for (actual, expected) in split
        .cumulative_decays()
        .into_iter()
        .zip(whole.cumulative_decays())
    {
        close(actual, expected);
    }
    close_channels(split.cumulative_released(), whole.cumulative_released());
    close_channels(
        split.cumulative_released(),
        channels(
            first.released.deposited_j + second.released.deposited_j,
            first.released.radiation_j + second.released.radiation_j,
            first.released.escaped_j + second.released.escaped_j,
        ),
    );
}

#[test]
fn g0_tiny_releases_survive_and_unresolved_parent_debits_consume_the_budget() {
    let model = model(
        vec![node("slow", 1e-20), node("stable", 0.0)],
        vec![branch(0, 1, 1.0, channels(3e5, 2e5, 1e5))],
    );
    let initial = DecayInventory::new(model, vec![1.0, 0.0]).unwrap();
    let allowed = DecayBudget {
        nuclei_tolerance: 1e-19,
        energy_tolerance_j: 1e-28,
        ..budget()
    };
    let (next, report) = initial
        .trial(Time::new(1.0).value(), allowed, || Ok(()))
        .unwrap();
    assert_eq!(next.populations()[0], 1.0);
    close(next.populations()[1], 1e-20);
    close(report.decays[0], 1e-20);
    close_channels(report.released, channels(3e-15, 2e-15, 1e-15));
    assert!(report.population_roundoff_nuclei >= 0.5e-20);
    assert!(report.population_roundoff_nuclei <= allowed.nuclei_tolerance);
    let mut refused = initial.clone();
    assert_eq!(
        refused.advance(
            1.0,
            DecayBudget {
                nuclei_tolerance: 0.0,
                ..allowed
            },
            || Ok(())
        ),
        Err(DecayError::Budget("population accounting roundoff"))
    );
    assert_eq!(refused, initial);
}

#[test]
fn g0_stable_inventories_and_zero_duration_leave_the_complete_history_unchanged() {
    let stable = model(vec![node("stable-a", 0.0), node("stable-b", 0.0)], vec![]);
    let inventory = DecayInventory::new(stable, vec![3.0, 7.0]).unwrap();
    let (same, report) = inventory
        .trial(Time::new(f64::MAX).value(), budget(), || Ok(()))
        .unwrap();
    assert_eq!(same, inventory);
    assert_eq!(report.decays, [0.0, 0.0]);
    assert!(report.branch_decays.is_empty());
    assert_eq!(report.released, EnergyChannels::default());
    assert_eq!(same.activity_bq().unwrap(), [0.0, 0.0]);
    assert_eq!(same.decay_power_w().unwrap(), EnergyChannels::default());
    let mut decayed = DecayInventory::new(
        single(Time::new(2.0), channels(1.0, 2.0, 3.0)),
        vec![4.0, 0.0],
    )
    .unwrap();
    decayed.advance(1.0, budget(), || Ok(())).unwrap();
    let saved = decayed.clone();
    let report = decayed
        .advance(Time::new(0.0).value(), budget(), || Ok(()))
        .unwrap();
    assert_eq!(decayed, saved);
    assert_eq!(report.decays, [0.0, 0.0]);
    assert_eq!(report.released, EnergyChannels::default());
}

#[test]
fn g0_physical_identity_energy_and_inventory_admission_refuse_invalid_data() {
    assert!(
        DecayModel::new(
            " ".into(),
            vec![node("stable", 0.0)],
            vec![],
            DecayLimits::default()
        )
        .is_err()
    );
    assert!(
        DecayModel::new(
            "dataset/v1".into(),
            vec![node("same", 0.0), node("same", 0.0)],
            vec![],
            DecayLimits::default()
        )
        .is_err()
    );
    for energy in [
        channels(-1.0, 0.0, 0.0),
        channels(0.0, f64::NAN, 0.0),
        channels(0.0, 0.0, f64::INFINITY),
    ] {
        assert!(
            DecayModel::new(
                "dataset/v1".into(),
                vec![node("parent", 1.0), node("stable", 0.0)],
                vec![branch(0, 1, 1.0, energy)],
                DecayLimits::default()
            )
            .is_err()
        );
    }
    let model = single(Time::new(1.0), channels(1.0, 0.0, 0.0));
    for populations in [
        vec![1.0],
        vec![1.0, -1.0],
        vec![f64::NAN, 0.0],
        vec![f64::MAX, f64::MAX],
    ] {
        assert!(DecayInventory::new(Arc::clone(&model), populations).is_err());
    }
    let mut inventory = DecayInventory::new(model, vec![7.0, 0.0]).unwrap();
    let saved = inventory.clone();
    assert!(matches!(
        inventory.advance(
            1.0,
            DecayBudget {
                energy_tolerance_j: f64::NAN,
                ..budget()
            },
            || Ok(())
        ),
        Err(DecayError::Invalid(_))
    ));
    assert_eq!(inventory, saved);
}

#[test]
fn g4_early_midway_and_publication_cancellation_preserve_inventory_and_plan() {
    let model = single(Time::new(2.0), channels(3.0, 4.0, 5.0));
    assert!(matches!(
        model.prepare(1.0, budget().propagation, || Err(DecayError::Budget(
            "caller checkpoint"
        ))),
        Err(DecayError::Budget("caller checkpoint"))
    ));
    let prepared = model.prepare(1.0, budget().propagation, || Ok(())).unwrap();
    let transition = prepared.propagator().transition_matrix().to_vec();
    let inventory = DecayInventory::new(model, vec![10.0, 3.0]).unwrap();
    let mut checkpoints = 0;
    inventory
        .trial_prepared(&prepared, budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert!(checkpoints > 3);
    for stop in [1, checkpoints / 2, checkpoints] {
        let mut candidate = inventory.clone();
        let mut calls = 0;
        let result = candidate.advance_prepared(&prepared, budget(), || {
            calls += 1;
            if calls == stop {
                Err(DecayError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err(DecayError::Cancelled));
        assert_eq!(candidate, inventory);
        assert_eq!(prepared.propagator().transition_matrix(), transition);
    }
}

#[test]
fn g4_energy_product_roundoff_and_propagation_range_refuse_without_publication() {
    let first_model = single(Time::new(1.0), channels(0.1, 0.0, 0.0));
    let mut inventory = DecayInventory::new(first_model, vec![14.0, 0.0]).unwrap();
    let saved = inventory.clone();
    assert_eq!(
        inventory.advance(
            1.0,
            DecayBudget {
                energy_tolerance_j: 0.0,
                ..budget()
            },
            || Ok(())
        ),
        Err(DecayError::Budget("emission energy accounting roundoff"))
    );
    assert_eq!(inventory, saved);
    let report = inventory.advance(1.0, budget(), || Ok(())).unwrap();
    assert!(report.energy_roundoff_j > 0.0);
    close(report.released.deposited_j, 0.7);
    let rate_one = model(
        vec![node("parent", 1.0), node("stable", 0.0)],
        vec![branch(0, 1, 1.0, channels(1.0, 0.0, 0.0))],
    );
    let mut huge = DecayInventory::new(rate_one, vec![1e200, 0.0]).unwrap();
    let saved = huge.clone();
    assert!(matches!(
        huge.advance(800.0, budget(), || Ok(())),
        Err(DecayError::Propagation(
            fs_la::decay::DecayError::Unrepresentable(_)
        ))
    ));
    assert_eq!(huge, saved);
}
