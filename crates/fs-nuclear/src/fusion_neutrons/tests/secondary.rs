//! G1/G3/G4 nonfission neutron production on the original propagators.
//! Coefficients and Q values are manufactured; isotope labels exercise exact
//! baryon/charge bookkeeping, not an evaluated blanket or breeding prediction.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState, ReactionBinding};
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationError, FusionNeutronActivationState,
};
use crate::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryEntry, FusionNeutronHistoryLimits,
    NeutronKineticHistoryBudget,
};
use crate::irradiation::NuclideId;
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;
use fs_thermochem::fusion::neutrons::FusionNeutronIntervalV1;

fn production_model(q_mev: f64, capture_rate: f64) -> Arc<NeutronKineticModel> {
    let mut supplied = data(&[2.0, 14.0], 1);
    supplied.reactions = vec![
        reaction(
            0,
            1,
            0.6,
            DiffusionReactionKind::NeutronProduction {
                neutron_yields: vec![2.0, 0.0],
            },
        ),
        reaction(0, 0, capture_rate, DiffusionReactionKind::Capture),
    ];
    Arc::new(
        NeutronKineticModel::new_with_reaction_energies(
            "manufactured-secondary-energy".into(),
            Arc::new(DiffusionModel::new(supplied, DiffusionLimits::default()).unwrap()),
            vec![2.0 * MEV, 14.0 * MEV],
            vec![0.25, 1.0],
            vec![
                NeutronReactionEnergy {
                    nuclear_q_j: q_mev * MEV,
                    threshold_j: (-q_mev).max(0.0) * 1.1 * MEV,
                },
                NeutronReactionEnergy::default(),
            ],
        )
        .unwrap(),
    )
}

#[test]
fn g1_secondary_neutrons_match_triangular_kinetics_with_signed_prompt_q() {
    for q_mev in [-1.5, 2.0] {
        let owner = production_model(q_mev, 0.3);
        let initial = NeutronKineticState::new(Arc::clone(&owner), vec![5.0, 80.0]).unwrap();
        let time = 0.8_f64;
        let (next, step) = initial
            .trial_with_rates(time, &[0.0; 2], budget(), || Ok(()))
            .unwrap();
        let fast = 80.0 * (-0.6 * time).exp();
        let slow = 5.0 * (-0.3 * time).exp()
            + 2.0 * 0.6 * 80.0 * ((-0.6 * time).exp() - (-0.3 * time).exp()) / (0.3 - 0.6);
        let events = 80.0 - fast;
        let captures = 5.0 + 2.0 * events - slow;
        close(next.neutrons().neutron_counts()[0], slow);
        close(next.neutrons().neutron_counts()[1], fast);
        close(step.diffusion.neutron_producing_reactions, events);
        close(step.diffusion.secondary_neutrons, 2.0 * events);
        close(step.diffusion.captures, captures);
        assert_eq!(step.diffusion.fissions, 0.0);
        assert_eq!(step.diffusion.prompt_neutrons, 0.0);
        assert!(step.diffusion.neutron_balance_residual.abs() < 1e-10);
        let available = (14.0 + q_mev - 2.0 * 2.0) * MEV * events;
        close(
            step.released.deposited_j,
            0.25 * available + 2.0 * MEV * captures,
        );
        close(step.released.radiation_j, 0.75 * available);
        close(step.nuclear_reaction_energy_j, q_mev * MEV * events);
        assert_eq!(
            step.cell_nuclear_reaction_energy_j,
            vec![step.nuclear_reaction_energy_j]
        );
        close(
            next.cell_nuclear_reaction_energy_j(0).unwrap(),
            step.nuclear_reaction_energy_j,
        );
        close(
            step.final_kinetic_energy_j + step.released.total_j(),
            step.initial_kinetic_energy_j + step.nuclear_reaction_energy_j,
        );
        assert!(step.energy_roundoff_j < budget().energy_tolerance_j);

        let (middle, first) = initial
            .trial_with_rates(0.3, &[0.0; 2], budget(), || Ok(()))
            .unwrap();
        let (restarted, second) = middle
            .trial_with_rates(0.5, &[0.0; 2], budget(), || Ok(()))
            .unwrap();
        for (&a, &b) in restarted
            .neutrons()
            .neutron_counts()
            .iter()
            .zip(next.neutrons().neutron_counts())
        {
            close(a, b);
        }
        close(
            first.nuclear_reaction_energy_j + second.nuclear_reaction_energy_j,
            step.nuclear_reaction_energy_j,
        );
        close(
            restarted.cell_released(0).unwrap().total_j(),
            step.released.total_j(),
        );
        assert_eq!(initial.neutrons().neutron_counts(), &[5.0, 80.0]);
    }
}

#[test]
fn g3_secondary_neutron_history_retains_endothermic_accounts_after_source_off() {
    let owner = production_model(-1.5, 0.0);
    let initial = NeutronKineticState::new(Arc::clone(&owner), vec![0.0; 2]).unwrap();
    let history = FusionNeutronHistory::project(
        owner,
        vec![FusionNeutronHistoryEntry {
            cell: 0,
            interval: FusionNeutronIntervalV1::new(
                Time::new(0.0),
                Time::new(0.4),
                emission(0.4, [40.0, 0.0], [14.0, 0.0]),
            )
            .unwrap(),
        }],
        FusionNeutronHistoryLimits::default(),
    )
    .unwrap();
    let (next, step) = initial
        .trial_with_history(
            Time::new(0.0),
            Time::new(1.0),
            &history,
            NeutronKineticHistoryBudget {
                transport: budget(),
                maximum_step_s: 0.25,
                maximum_substeps: 16,
            },
            || Ok(()),
        )
        .unwrap();
    let fast_at_shutdown = 100.0 / 0.6 * -(-0.6_f64 * 0.4).exp_m1();
    let fast = fast_at_shutdown * (-0.6_f64 * 0.6).exp();
    let events = 40.0 - fast;
    close(next.neutrons().neutron_counts()[1], fast);
    close(next.neutrons().neutron_counts()[0], 2.0 * events);
    close(step.nuclear_reaction_energy_j, -1.5 * MEV * events);
    close(
        step.cell_nuclear_reaction_energy_j[0],
        step.nuclear_reaction_energy_j,
    );
    close(step.released.total_j(), 8.5 * MEV * events);
    close(step.injected_kinetic_energy_j, 14.0 * MEV * 40.0);
    close(
        step.final_kinetic_energy_j + step.released.total_j(),
        step.injected_kinetic_energy_j + step.nuclear_reaction_energy_j,
    );
    assert!(step.energy_roundoff_j <= budget().energy_tolerance_j);
    assert!(events > 40.0 - fast_at_shutdown);
}

fn activation_budget() -> FusionNeutronActivationBudget {
    FusionNeutronActivationBudget {
        spatial: DiffusionBurnupBudget {
            diffusion: budget().diffusion,
            burnup: BurnupBudget {
                limits: BurnupLimits::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-19,
            },
            max_iterations: 64,
            reaction_tolerance: 1e-11,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-18,
    }
}

fn finite_parts(
    low_sigma_at_hot: f64,
    inventory_energy: EnergyChannels,
) -> (Arc<NeutronKineticModel>, DiffusionBurnupState) {
    let mut supplied = data(&[2.0, 14.0], 1);
    let groups = Arc::clone(&supplied.groups);
    let nuclide = |z, a| BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
    };
    let target = nuclide(4, 9);
    let table = Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: "manufactured-Be9-two-neutron-table".into(),
                nuclide: target.id,
                reaction: NeutronReaction::TwoNeutronProduction,
                groups: Arc::clone(&groups),
                temperatures_k: vec![300.0, 600.0],
                cross_sections_m2: vec![vec![0.0, 0.01], vec![low_sigma_at_hot, 0.02]],
            },
            CrossSectionLimits::default(),
        )
        .unwrap(),
    );
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-Be9-to-two-He4".into(),
                groups,
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![target, nuclide(2, 4)],
                reactions: vec![BurnupReaction {
                    data_id: "manufactured-explicit-target-removal-and-products".into(),
                    parent: 0,
                    rate: BurnupRate::TabulatedNeutronProduction {
                        cross_sections: table,
                    },
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: 2.0,
                    }],
                    energy_per_event: inventory_energy,
                }],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    supplied.reactions = (0..2)
        .map(|group| {
            reaction(
                0,
                group,
                if group == 0 { 0.0 } else { 1.0 },
                DiffusionReactionKind::NeutronProduction {
                    neutron_yields: vec![2.0, 0.0],
                },
            )
        })
        .collect();
    let diffusion = Arc::new(DiffusionModel::new(supplied, DiffusionLimits::default()).unwrap());
    let model = Arc::new(
        NeutronKineticModel::new_with_reaction_energies(
            "manufactured-finite-secondary-energy".into(),
            Arc::clone(&diffusion),
            vec![2.0 * MEV, 14.0 * MEV],
            vec![1.0; 2],
            vec![
                NeutronReactionEnergy {
                    nuclear_q_j: -1.5 * MEV,
                    threshold_j: 1.7 * MEV
                };
                2
            ],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0, 20.0], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, vec![100.0, 0.0]).unwrap()],
        (0..2)
            .map(|group| ReactionBinding {
                diffusion_reaction: group,
                burnup_reaction: 0,
            })
            .collect(),
    )
    .unwrap();
    (model, spatial)
}

#[test]
fn g1_secondary_finite_targets_products_and_signed_energy_converge_to_independent_solution() {
    let (model, spatial) = finite_parts(0.0, EnergyChannels::default());
    let initial = FusionNeutronActivationState::new(model, spatial).unwrap();
    let duration = 0.7_f64;
    let e = (-80.0 * 0.01 * duration).exp();
    let expected_fast = 80.0 * 20.0 * e / (80.0 + 20.0 * (1.0 - e));
    let expected_events = 20.0 - expected_fast;
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut current = initial.clone();
        let mut produced = 0.0;
        let mut signed_q = 0.0;
        for _ in 0..steps {
            let (next, step) = current
                .trial_with_rates(
                    duration / f64::from(steps),
                    &[300.0],
                    &[0.0; 2],
                    activation_budget(),
                    || Ok(()),
                )
                .unwrap();
            let events = step.spatial.burnup[0].neutron_producing_reactions;
            close(step.spatial.neutrons.neutron_producing_reactions, events);
            close(step.spatial.neutrons.secondary_neutrons, 2.0 * events);
            close(step.kinetic.nuclear_reaction_energy_j, -1.5 * MEV * events);
            close(
                step.nuclear_generated_energy_j,
                step.kinetic.nuclear_reaction_energy_j,
            );
            assert_eq!(step.spatial.burnup[0].captures, 0.0);
            assert_eq!(step.spatial.burnup[0].fissions, 0.0);
            assert_eq!(step.nuclear_released, EnergyChannels::default());
            assert!(step.energy_roundoff_j <= activation_budget().energy_tolerance_j);
            produced += events;
            signed_q += step.nuclear_generated_energy_j;
            current = next;
        }
        let neutron = current.neutrons().neutron_counts();
        let inventory = current.inventories()[0].populations();
        close(100.0 - inventory[0], produced);
        close(inventory[1], 2.0 * produced);
        close(neutron[0], 2.0 * produced);
        close(neutron[0] + neutron[1], 20.0 + produced);
        close(
            9.0 * inventory[0] + 4.0 * inventory[1] + neutron.iter().sum::<f64>(),
            920.0,
        );
        close(4.0 * inventory[0] + 2.0 * inventory[1], 400.0);
        close(current.cell_prompt_nuclear_energy_j(0).unwrap(), signed_q);
        close(
            current.kinetic_energy_j().unwrap() + current.cell_released(0).unwrap().total_j(),
            initial.kinetic_energy_j().unwrap() + signed_q,
        );
        errors.push((produced - expected_events).abs());
    }
    assert!(
        errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 5e-4, "{errors:?}");
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0]);
}

#[test]
fn g0_g4_secondary_sources_refuse_impossible_energy_wrong_bindings_and_late_cancellation() {
    let owner = production_model(-1.5, 0.0);
    let original = owner.diffusion().data().clone();
    let make = |data: DiffusionData, energy: NeutronReactionEnergy| {
        NeutronKineticModel::new_with_reaction_energies(
            "refusal".into(),
            Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
            vec![2.0 * MEV, 14.0 * MEV],
            vec![1.0; 2],
            vec![energy, NeutronReactionEnergy::default()],
        )
    };
    for energy in [
        NeutronReactionEnergy {
            nuclear_q_j: -11.0 * MEV,
            threshold_j: 12.0 * MEV,
        },
        NeutronReactionEnergy {
            nuclear_q_j: -1.5 * MEV,
            threshold_j: 15.0 * MEV,
        },
        NeutronReactionEnergy {
            nuclear_q_j: -1.5 * MEV,
            threshold_j: 0.0,
        },
        NeutronReactionEnergy {
            nuclear_q_j: f64::NAN,
            threshold_j: 0.0,
        },
    ] {
        assert!(make(original.clone(), energy).is_err());
    }
    assert!(
        NeutronKineticModel::new(
            "implicit-production-energy".into(),
            Arc::clone(owner.diffusion()),
            vec![2.0 * MEV, 14.0 * MEV],
            vec![1.0; 2]
        )
        .is_err()
    );

    let (model, spatial) = finite_parts(
        0.0,
        EnergyChannels {
            deposited_j: MEV,
            ..EnergyChannels::default()
        },
    );
    assert!(FusionNeutronActivationState::new(model, spatial).is_err());
    let (model, spatial) = finite_parts(0.01, EnergyChannels::default());
    let initial = FusionNeutronActivationState::new(model, spatial).unwrap();
    assert!(matches!(
        initial.trial_with_rates(0.1, &[600.0], &[0.0; 2], activation_budget(), || Ok(())),
        Err(FusionNeutronActivationError::Kinetic(
            NeutronKineticError::Invalid(_)
        ))
    ));
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0, 20.0]);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0]);

    let (model, spatial) = finite_parts(0.0, EnergyChannels::default());
    for kind in [
        DiffusionReactionKind::Capture,
        DiffusionReactionKind::NeutronProduction {
            neutron_yields: vec![3.0, 0.0],
        },
    ] {
        let mut changed = model.diffusion().data().clone();
        changed.reactions[1].kind = kind;
        let changed = Arc::new(DiffusionModel::new(changed, DiffusionLimits::default()).unwrap());
        assert!(
            DiffusionBurnupState::new(
                DiffusionState::new(changed, vec![0.0, 20.0], vec![]).unwrap(),
                spatial.inventories().to_vec(),
                spatial.bindings().to_vec(),
            )
            .is_err()
        );
    }
    let initial = FusionNeutronActivationState::new(model, spatial).unwrap();
    let mut calls = 0;
    initial
        .trial_with_rates(0.1, &[300.0], &[0.0; 2], activation_budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        initial.trial_with_rates(0.1, &[300.0], &[0.0; 2], activation_budget(), || {
            seen += 1;
            if seen == calls {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(FusionNeutronActivationError::Source(
            FissionError::Cancelled
        ))
    ));
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0, 20.0]);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0]);
    assert_eq!(initial.cell_prompt_nuclear_energy_j(0), Some(0.0));
}
