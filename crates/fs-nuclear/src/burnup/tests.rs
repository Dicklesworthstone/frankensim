//! Synthetic equation tests: none of these coefficients are evaluated isotope data.

use super::*;
use fs_qty::Frequency;

fn channels(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn limits() -> BurnupLimits {
    BurnupLimits {
        max_neutron_groups: 8,
        max_reactions: 32,
        max_products: 64,
        propagation: DecayLimits::default(),
    }
}

fn budget() -> BurnupBudget {
    BurnupBudget {
        limits: limits(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (allowance {tolerance:.4e})"
    );
}

fn node(index: u16, rate: f64) -> BurnupNuclide {
    BurnupNuclide {
        id: NuclideId::new(10, 25 + index, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    }
}

fn product_to(nuclide: usize, yield_per_event: f64) -> BurnupProduct {
    BurnupProduct {
        nuclide,
        yield_per_event,
    }
}

fn fixture(k_f: f64, k_c: f64, lambda: f64, populations: Vec<f64>) -> (BurnupInventory, GroupFlux) {
    let groups =
        Arc::new(NeutronGroups::new("synthetic-burnup-grid".into(), vec![0.0, 1.0, 3.0]).unwrap());
    let mut reactions = vec![
        BurnupReaction {
            data_id: "synthetic-fission-channel".into(),
            parent: 0,
            rate: BurnupRate::Fission {
                cross_sections_m2: vec![k_f / 2.0, k_f / 4.0],
            },
            products: vec![product_to(1, 2.0)],
            energy_per_event: channels(3.0, 2.0, 1.0),
        },
        BurnupReaction {
            data_id: "synthetic-capture-channel".into(),
            parent: 0,
            rate: BurnupRate::Capture {
                cross_sections_m2: vec![k_c, 0.0],
            },
            products: vec![product_to(3, 1.0)],
            energy_per_event: channels(1.0, 2.0, 0.0),
        },
    ];
    if lambda > 0.0 {
        reactions.push(BurnupReaction {
            data_id: "synthetic-product-decay".into(),
            parent: 1,
            rate: BurnupRate::Decay {
                branch_fraction: 1.0,
            },
            products: vec![product_to(2, 1.0)],
            energy_per_event: channels(0.5, 0.5, 1.0),
        });
    }
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "synthetic-competing-reactions-with-fission-products".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 250.0,
                temperature_max_k: 1000.0,
                nuclides: vec![node(0, 0.0), node(1, lambda), node(2, 0.0), node(3, 0.0)],
                reactions,
            },
            limits(),
        )
        .unwrap(),
    );
    (
        BurnupInventory::new(model, populations).unwrap(),
        GroupFlux::new(groups, vec![1.0, 2.0]).unwrap(),
    )
}

#[test]
fn g1_multigroup_fission_capture_and_product_decay_share_one_interval() {
    let (mut inventory, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 0.0, 0.0, 0.0]);
    let rates = inventory
        .reaction_rates(&flux, 300.0, limits(), || Ok(()))
        .unwrap();
    close(rates[0], 12.0, 1e-14);
    close(rates[1], 8.0, 1e-14);
    assert_eq!(rates[2], 0.0);
    let report = inventory
        .advance(2.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let surviving_parent = 100.0 * (-0.5_f64).exp();
    let fissions = 60.0 * -(-0.5_f64).exp_m1();
    let captures = 40.0 * -(-0.5_f64).exp_m1();
    let fragments = 24.0 / 0.15 * ((-0.5_f64).exp() - (-0.875_f64).exp());
    let decays = 2.0 * fissions - fragments;
    for (&actual, expected) in
        inventory
            .populations()
            .iter()
            .zip([surviving_parent, fragments, decays, captures])
    {
        close(actual, expected, 1e-9);
    }
    close(report.fissions, fissions, 1e-9);
    close(report.captures, captures, 1e-9);
    close(report.decays, decays, 1e-9);
    close(
        report.released.deposited_j,
        3.0 * fissions + captures + 0.5 * decays,
        3e-9,
    );
    close(
        report.released.radiation_j,
        2.0 * fissions + 2.0 * captures + 0.5 * decays,
        3e-9,
    );
    close(report.released.escaped_j, fissions + decays, 3e-9);
    close(
        report.generated_energy_j,
        6.0 * fissions + 3.0 * captures + 2.0 * decays,
        8e-9,
    );
    close(report.generated_energy_j, report.released.total_j(), 1e-12);
    close(report.fission_energy_j, 6.0 * fissions, 3e-9);
    close(report.capture_energy_j, 3.0 * captures, 3e-9);
    close(report.decay_energy_j, 2.0 * decays, 3e-9);
    assert!(inventory.populations().iter().all(|value| *value >= 0.0));
    assert!(
        report.decays > 0.0,
        "products must decay in their birth interval"
    );
    close(sum(inventory.populations()), 100.0 + fissions, 2e-9);
    assert!(
        report
            .population_balance_residuals
            .iter()
            .all(|r| r.abs() <= 1e-8)
    );
}

#[test]
fn g1_equal_parent_and_product_rates_use_the_repeated_root_solution() {
    let (mut inventory, flux) = fixture(0.125, 0.125, 0.25, vec![120.0, 7.0, 0.0, 0.0]);
    let report = inventory
        .advance(3.0, &flux, 500.0, budget(), || Ok(()))
        .unwrap();
    let survival = (-0.75_f64).exp();
    let product = (7.0 + 2.0 * 0.125 * 120.0 * 3.0) * survival;
    close(inventory.populations()[0], 120.0 * survival, 1e-9);
    close(inventory.populations()[1], product, 1e-9);
    close(report.fissions, 60.0 * (1.0 - survival), 1e-9);
    close(report.captures, report.fissions, 1e-12);
    close(report.decays, 7.0 + 2.0 * report.fissions - product, 2e-9);
}

#[test]
fn g1_capture_products_can_capture_again_before_the_interval_ends() {
    let (template, flux) = fixture(0.0, 0.2, 0.0, vec![100.0, 0.0, 0.0, 0.0]);
    let mut data = template.model.data.clone();
    data.reactions[1].products = vec![product_to(1, 1.0)];
    data.reactions.push(BurnupReaction {
        data_id: "second-synthetic-capture".into(),
        parent: 1,
        rate: BurnupRate::Capture {
            cross_sections_m2: vec![0.2, 0.0],
        },
        products: vec![product_to(2, 1.0)],
        energy_per_event: channels(4.0, 0.0, 0.0),
    });
    let mut inventory = BurnupInventory::new(
        Arc::new(BurnupModel::new(data, limits()).unwrap()),
        vec![100.0, 0.0, 0.0, 0.0],
    )
    .unwrap();
    let report = inventory
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let survival = (-0.4_f64).exp();
    let second_captures = 100.0 * (1.0 - 1.4 * survival);
    close(inventory.populations()[0], 100.0 * survival, 1e-9);
    close(inventory.populations()[1], 40.0 * survival, 1e-9);
    close(inventory.populations()[2], second_captures, 1e-9);
    close(report.reaction_counts[2], second_captures, 1e-9);
    close(
        report.released.deposited_j,
        report.reaction_counts[1] + 4.0 * second_captures,
        2e-9,
    );
    assert_eq!(report.fissions, 0.0);
}

#[test]
fn g1_shutdown_keeps_product_history_and_emits_only_remaining_decay_energy() {
    let (mut inventory, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 0.0, 0.0, 0.0]);
    inventory
        .advance(2.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let before = inventory.clone();
    let zero = GroupFlux::new(Arc::clone(flux.groups()), vec![0.0, 0.0]).unwrap();
    let report = inventory
        .advance(4.0, &zero, 300.0, budget(), || Ok(()))
        .unwrap();
    let decays = before.populations()[1] * -(-1.4_f64).exp_m1();
    assert_eq!(report.fissions, 0.0);
    assert_eq!(report.captures, 0.0);
    assert_eq!(inventory.populations()[0], before.populations()[0]);
    assert_eq!(inventory.populations()[3], before.populations()[3]);
    close(report.decays, decays, 1e-9);
    close(
        inventory.populations()[1],
        before.populations()[1] * (-1.4_f64).exp(),
        1e-9,
    );
    close(
        inventory.populations()[2],
        before.populations()[2] + decays,
        1e-9,
    );
    close(report.released.deposited_j, 0.5 * decays, 1e-9);
    close(report.released.radiation_j, 0.5 * decays, 1e-9);
    close(report.released.escaped_j, decays, 1e-9);
    close(
        inventory.cumulative_released().total_j(),
        before.cumulative_released().total_j() + 2.0 * decays,
        3e-9,
    );
}

#[test]
fn g1_partial_group_channels_preserve_energy_dependent_fission_yields() {
    let (template, flux) = fixture(0.0, 0.0, 0.0, vec![100.0, 0.0, 0.0, 0.0]);
    let mut data = template.model.data.clone();
    data.reactions = vec![
        BurnupReaction {
            data_id: "synthetic-low-group-fission".into(),
            parent: 0,
            rate: BurnupRate::Fission {
                cross_sections_m2: vec![0.1, 0.0],
            },
            products: vec![product_to(1, 1.5), product_to(2, 0.5)],
            energy_per_event: channels(2.0, 0.0, 0.0),
        },
        BurnupReaction {
            data_id: "synthetic-high-group-fission".into(),
            parent: 0,
            rate: BurnupRate::Fission {
                cross_sections_m2: vec![0.0, 0.15],
            },
            products: vec![product_to(1, 0.25), product_to(2, 1.75)],
            energy_per_event: channels(3.0, 0.0, 0.0),
        },
    ];
    let mut inventory = BurnupInventory::new(
        Arc::new(BurnupModel::new(data, limits()).unwrap()),
        vec![100.0, 0.0, 0.0, 0.0],
    )
    .unwrap();
    let report = inventory
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let removed = 100.0 * -(-0.8_f64).exp_m1();
    close(report.reaction_counts[0], 0.25 * removed, 1e-9);
    close(report.reaction_counts[1], 0.75 * removed, 1e-9);
    close(
        inventory.populations()[1],
        (0.25 * 1.5 + 0.75 * 0.25) * removed,
        2e-9,
    );
    close(
        inventory.populations()[2],
        (0.25 * 0.5 + 0.75 * 1.75) * removed,
        2e-9,
    );
    close(
        inventory.populations()[1] + inventory.populations()[2],
        2.0 * removed,
        2e-9,
    );
    close(report.released.deposited_j, 2.75 * removed, 3e-9);
}

#[test]
fn g1_spontaneous_branches_retain_their_own_products_and_energy() {
    let (template, flux) = fixture(0.0, 0.0, 0.4, vec![0.0, 40.0, 0.0, 0.0]);
    let mut data = template.model.data.clone();
    data.reactions[2].rate = BurnupRate::Decay {
        branch_fraction: 0.25,
    };
    data.reactions.push(BurnupReaction {
        data_id: "alternate-synthetic-decay".into(),
        parent: 1,
        rate: BurnupRate::Decay {
            branch_fraction: 0.75,
        },
        products: vec![product_to(3, 1.0)],
        energy_per_event: channels(2.0, 3.0, 0.0),
    });
    let mut inventory = BurnupInventory::new(
        Arc::new(BurnupModel::new(data, limits()).unwrap()),
        vec![0.0, 40.0, 0.0, 0.0],
    )
    .unwrap();
    let report = inventory
        .advance(1.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let removed = 40.0 * -(-0.4_f64).exp_m1();
    close(report.decays, removed, 1e-9);
    close(inventory.populations()[2], 0.25 * removed, 1e-9);
    close(inventory.populations()[3], 0.75 * removed, 1e-9);
    close(
        report.released.deposited_j,
        (0.25 * 0.5 + 0.75 * 2.0) * removed,
        2e-9,
    );
    close(
        report.released.radiation_j,
        (0.25 * 0.5 + 0.75 * 3.0) * removed,
        2e-9,
    );
    close(report.released.escaped_j, 0.25 * removed, 1e-9);
}

#[test]
fn g3_constant_operator_partitioning_and_checkpoint_replay_agree() {
    let (initial, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 7.0, 3.0, 1.0]);
    let (whole, whole_step) = initial
        .trial(3.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let (same, same_step) = initial
        .trial(3.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(whole, same);
    assert_eq!(whole_step, same_step);
    let mut split = initial;
    for _ in 0..6 {
        split
            .advance(0.5, &flux, 300.0, budget(), || Ok(()))
            .unwrap();
    }
    for (&actual, &expected) in split.populations().iter().zip(whole.populations()) {
        close(actual, expected, 3e-9);
    }
    for (actual, expected) in split
        .cumulative_reaction_counts()
        .iter()
        .zip(whole.cumulative_reaction_counts())
    {
        close(*actual, expected, 3e-9);
    }
    close(
        split.cumulative_released().deposited_j,
        whole.cumulative_released().deposited_j,
        5e-9,
    );
    close(
        split.cumulative_released().radiation_j,
        whole.cumulative_released().radiation_j,
        5e-9,
    );
    close(
        split.cumulative_released().escaped_j,
        whole.cumulative_released().escaped_j,
        5e-9,
    );
}

#[test]
fn g0_repeated_steps_share_preparation_across_admitted_temperatures() {
    let (initial, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 0.0, 0.0, 0.0]);
    let (first, _) = initial
        .trial(0.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let (second, _) = first.trial(0.5, &flux, 700.0, budget(), || Ok(())).unwrap();
    assert!(Arc::ptr_eq(
        first.prepared.as_ref().unwrap(),
        second.prepared.as_ref().unwrap()
    ));
    let mut other_budget = budget();
    other_budget.limits.max_reactions += 1;
    let (third, _) = second
        .trial(0.5, &flux, 700.0, other_budget, || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(
        second.prepared.as_ref().unwrap(),
        third.prepared.as_ref().unwrap()
    ));
    assert!(
        second
            .trial(0.5, &flux, 1001.0, budget(), || Ok(()))
            .is_err()
    );
    let prepared = initial
        .model
        .prepare(0.5, &flux, 300.0, limits(), || Ok(()))
        .unwrap();
    let (explicit, explicit_step) = initial
        .trial_prepared(&prepared, budget(), || Ok(()))
        .unwrap();
    let (ordinary, ordinary_step) = initial
        .trial(0.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(explicit, ordinary);
    assert_eq!(explicit_step, ordinary_step);
    let different_model = Arc::new(BurnupModel::new(initial.model.data.clone(), limits()).unwrap());
    let rebound = BurnupInventory::new(different_model, initial.populations.clone()).unwrap();
    assert!(
        rebound
            .trial_prepared(&prepared, budget(), || Ok(()))
            .is_err()
    );
}

#[test]
fn g4_cancellation_resource_and_accounting_refusals_are_atomic() {
    let (mut inventory, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 0.0, 0.0, 0.0]);
    let before = inventory.clone();
    let expected = FissionError::Budget("synthetic caller checkpoint");
    let prepared = inventory
        .model
        .prepare(2.0, &flux, 300.0, limits(), || Ok(()))
        .unwrap();
    let mut polls = 0;
    let error = inventory
        .advance_prepared(&prepared, budget(), || {
            polls += 1;
            if polls == 3 {
                Err(expected.clone())
            } else {
                Ok(())
            }
        })
        .unwrap_err();
    assert_eq!(error, BurnupError::Fission(expected));
    assert_eq!(inventory, before);
    let mut insufficient = budget();
    insufficient.limits.max_reactions = 1;
    assert!(
        inventory
            .advance(2.0, &flux, 300.0, insufficient, || Ok(()))
            .is_err()
    );
    assert_eq!(inventory, before);
    let (mut tiny, tiny_flux) = fixture(1e-20, 0.0, 0.0, vec![1e20, 0.0, 0.0, 0.0]);
    let tiny_before = tiny.clone();
    let mut zero_account = budget();
    zero_account.nuclei_tolerance = 0.0;
    assert!(
        tiny.advance(1.0, &tiny_flux, 300.0, zero_account, || Ok(()))
            .is_err()
    );
    assert_eq!(tiny, tiny_before);
    let mut admitted = budget();
    admitted.nuclei_tolerance = 2.0;
    let report = tiny
        .advance(1.0, &tiny_flux, 300.0, admitted, || Ok(()))
        .unwrap();
    close(report.fissions, 1.0, 1e-12);
    close(tiny.populations()[1], 2.0, 1e-12);
    assert!(report.population_roundoff_nuclei >= 1.0);
}

#[test]
fn g0_zero_duration_and_explicitly_stable_zero_flux_are_identities() {
    let (inventory, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 7.0, 3.0, 1.0]);
    let (same, report) = inventory
        .trial(0.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(same, inventory);
    assert_eq!(report.released, EnergyChannels::default());
    assert_eq!(report.reaction_counts, vec![0.0; 3]);
    let (stable, flux) = fixture(0.0, 0.0, 0.0, vec![100.0, 7.0, 3.0, 1.0]);
    let (same, report) = stable
        .trial(1e100, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(same, stable);
    assert_eq!(report.generated_energy_j, 0.0);
}

#[test]
fn g0_identity_branch_data_and_domains_are_admitted_before_evolution() {
    let (inventory, flux) = fixture(0.12, 0.08, 0.35, vec![100.0, 0.0, 0.0, 0.0]);
    let mut data = inventory.model.data.clone();
    data.nuclides[1].id = data.nuclides[0].id;
    assert!(BurnupModel::new(data, limits()).is_err());
    let mut data = inventory.model.data.clone();
    data.reactions[2].rate = BurnupRate::Decay {
        branch_fraction: 0.9,
    };
    assert!(BurnupModel::new(data, limits()).is_err());
    let mut data = inventory.model.data.clone();
    data.reactions[0].products[0].yield_per_event = f64::NAN;
    assert!(BurnupModel::new(data, limits()).is_err());
    let mut data = inventory.model.data.clone();
    data.reactions[0].rate = BurnupRate::Fission {
        cross_sections_m2: vec![0.1],
    };
    assert!(BurnupModel::new(data, limits()).is_err());
    let other_grid = Arc::new(
        NeutronGroups::new(
            "different-data-identity".into(),
            flux.groups().edges_j().to_vec(),
        )
        .unwrap(),
    );
    let wrong_flux = GroupFlux::new(other_grid, vec![0.0, 0.0]).unwrap();
    assert!(
        inventory
            .trial(0.0, &wrong_flux, 300.0, budget(), || Ok(()))
            .is_err()
    );
    assert!(
        inventory
            .trial(0.0, &flux, 249.0, budget(), || Ok(()))
            .is_err()
    );
    assert!(
        inventory
            .trial(-1.0, &flux, 300.0, budget(), || Ok(()))
            .is_err()
    );
    let mut short = budget();
    short.limits.max_products = 1;
    assert!(
        inventory
            .trial(1.0, &flux, 300.0, short, || Ok(()))
            .is_err()
    );
    assert!(BurnupInventory::new(Arc::clone(inventory.model()), vec![1.0]).is_err());
}
