//! Manufactured prompt-fission and radioactive-product histories. These check
//! equations and energy ownership; isotope labels and Q values are not data.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::diffusion::DelayedGroup;
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState, ReactionBinding};
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationError, FusionNeutronActivationState,
};
use crate::irradiation::NuclideId;
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn multiplying(q_mev: f64) -> Result<NeutronKineticModel, NeutronKineticError> {
    let mut supplied = data(&[2.0], 1);
    supplied.reactions = vec![
        reaction(
            0,
            0,
            0.4,
            DiffusionReactionKind::Fission {
                prompt_yields: vec![2.5],
                delayed_yields: vec![],
            },
        ),
        reaction(0, 0, 0.2, DiffusionReactionKind::Capture),
    ];
    supplied.leakage.push(LeakageBoundary {
        cell: 0,
        conductance_m2: vec![0.1],
    });
    NeutronKineticModel::new_with_reaction_energies(
        "manufactured-prompt-fission-energy".into(),
        Arc::new(DiffusionModel::new(supplied, DiffusionLimits::default()).unwrap()),
        vec![2.0 * MEV],
        vec![0.75, 0.5],
        vec![
            NeutronReactionEnergy {
                nuclear_q_j: q_mev * MEV,
                threshold_j: (-q_mev).max(0.0) * MEV,
            },
            NeutronReactionEnergy::default(),
        ],
    )
}

#[test]
fn g1_prompt_multiplication_capture_and_uniform_feed_share_one_energy_balance() {
    let owner = Arc::new(multiplying(5.0).unwrap());
    let initial = NeutronKineticState::new(Arc::clone(&owner), vec![10.0]).unwrap();
    let duration = 1.3_f64;
    let source = 4.0;
    let growth = 0.3;
    let gained = (growth * duration).exp_m1();
    let expected = 10.0 * (1.0 + gained) + source * gained / growth;
    let residence = 10.0 * gained / growth + source * (gained / growth - duration) / growth;
    let fissions = 0.4 * residence;
    let captures = 0.2 * residence;
    let leaked = 0.1 * residence;
    let (next, step) = initial
        .trial_with_rates(duration, &[source], budget(), || Ok(()))
        .unwrap();
    close(next.neutrons().neutron_counts()[0], expected);
    close(step.diffusion.integrated_fluence_m2[0], residence);
    close(step.diffusion.fissions, fissions);
    close(step.diffusion.prompt_neutrons, 2.5 * fissions);
    close(step.diffusion.captures, captures);
    close(step.diffusion.leaked_neutrons, leaked);
    assert_eq!(step.diffusion.neutron_producing_reactions, 0.0);
    assert_eq!(step.diffusion.secondary_neutrons, 0.0);
    close(step.nuclear_reaction_energy_j, 5.0 * MEV * fissions);
    close(step.released.deposited_j, MEV * (1.5 * fissions + captures));
    close(step.released.radiation_j, MEV * (0.5 * fissions + captures));
    close(step.released.escaped_j, 2.0 * MEV * leaked);
    close(
        next.kinetic_energy_j().unwrap() + step.released.total_j(),
        2.0 * MEV * (10.0 + source * duration) + 5.0 * MEV * fissions,
    );
    assert!(step.energy_roundoff_j <= budget().energy_tolerance_j);
    let (middle, first) = initial
        .trial_with_rates(0.4, &[source], budget(), || Ok(()))
        .unwrap();
    let (restarted, second) = middle
        .trial_with_rates(0.9, &[source], budget(), || Ok(()))
        .unwrap();
    close(restarted.neutrons().neutron_counts()[0], expected);
    close(
        first.diffusion.fissions + second.diffusion.fissions,
        fissions,
    );
    close(
        restarted.cell_released(0).unwrap().total_j(),
        step.released.total_j(),
    );
}

#[test]
fn g1_prompt_fission_retains_signed_conversion_and_distinct_outgoing_energies() {
    // Negative Q is an authored manufactured energy-account case, not a
    // proposed isotope evaluation. The outgoing-energy condition still holds.
    for q_mev in [-1.5, 2.0] {
        let mut supplied = data(&[2.0, 14.0], 1);
        supplied.reactions = vec![
            reaction(
                0,
                1,
                0.6,
                DiffusionReactionKind::Fission {
                    prompt_yields: vec![2.0, 0.0],
                    delayed_yields: vec![],
                },
            ),
            reaction(0, 0, 0.3, DiffusionReactionKind::Capture),
        ];
        let owner = Arc::new(
            NeutronKineticModel::new_with_reaction_energies(
                "manufactured-signed-prompt-fission".into(),
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
        );
        let initial = NeutronKineticState::new(owner, vec![5.0, 80.0]).unwrap();
        let duration = 0.8_f64;
        let (next, step) = initial
            .trial_with_rates(duration, &[0.0; 2], budget(), || Ok(()))
            .unwrap();
        let fast = 80.0 * (-0.6 * duration).exp();
        let slow = 5.0 * (-0.3 * duration).exp()
            + 2.0 * 0.6 * 80.0 * ((-0.6 * duration).exp() - (-0.3 * duration).exp()) / (0.3 - 0.6);
        let fissions = 80.0 - fast;
        let captures = 5.0 + 2.0 * fissions - slow;
        let available = (14.0 + q_mev - 4.0) * MEV * fissions;
        close(next.neutrons().neutron_counts()[0], slow);
        close(next.neutrons().neutron_counts()[1], fast);
        close(step.diffusion.fissions, fissions);
        close(step.diffusion.prompt_neutrons, 2.0 * fissions);
        close(step.diffusion.captures, captures);
        close(
            step.released.deposited_j,
            0.25 * available + 2.0 * MEV * captures,
        );
        close(step.released.radiation_j, 0.75 * available);
        close(step.nuclear_reaction_energy_j, q_mev * MEV * fissions);
        close(
            step.final_kinetic_energy_j + step.released.total_j(),
            step.initial_kinetic_energy_j + step.nuclear_reaction_energy_j,
        );
    }
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
            max_iterations: 96,
            reaction_tolerance: 1e-12,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-18,
    }
}

fn finite_parts(
    source_atoms: f64,
    resident_neutrons: f64,
    prompt_q_mev: f64,
    inventory_prompt_mev: f64,
) -> (Arc<NeutronKineticModel>, DiffusionBurnupState) {
    let mut supplied = data(&[2.0], 1);
    let groups = Arc::clone(&supplied.groups);
    let nuclide = |z, a, rate| BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    };
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-radioactive-source-and-fission-stock".into(),
                groups,
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    nuclide(7, 16, 4.0),
                    nuclide(8, 15, 0.0),
                    nuclide(10, 30, 0.0),
                    nuclide(5, 15, 0.1),
                    nuclide(6, 15, 0.0),
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-source-decay".into(),
                        parent: 0,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 0.25 * MEV,
                            radiation_j: 3.0 * MEV,
                            escaped_j: 0.1 * MEV,
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-finite-fission".into(),
                        parent: 2,
                        rate: BurnupRate::Fission {
                            cross_sections_m2: vec![0.03],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 3,
                            yield_per_event: 2.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: inventory_prompt_mev * MEV,
                            radiation_j: 0.0,
                            escaped_j: 0.0,
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-product-decay".into(),
                        parent: 3,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 4,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 0.5 * MEV,
                            radiation_j: 0.25 * MEV,
                            escaped_j: 0.1 * MEV,
                        },
                    },
                ],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    supplied.reactions.push(reaction(
        0,
        0,
        1.5,
        DiffusionReactionKind::Fission {
            prompt_yields: vec![1.0],
            delayed_yields: vec![],
        },
    ));
    supplied.leakage.push(LeakageBoundary {
        cell: 0,
        conductance_m2: vec![3.0],
    });
    let diffusion = Arc::new(DiffusionModel::new(supplied, DiffusionLimits::default()).unwrap());
    let model = Arc::new(
        NeutronKineticModel::new_with_reaction_energies(
            "manufactured-finite-prompt-fission-energy".into(),
            Arc::clone(&diffusion),
            vec![2.0 * MEV],
            vec![0.8],
            vec![NeutronReactionEnergy {
                nuclear_q_j: prompt_q_mev * MEV,
                threshold_j: 0.0,
            }],
        )
        .unwrap(),
    );
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "manufactured-source-one-neutron-per-decay".into(),
            Arc::clone(&nuclear),
            vec![vec![1.0], vec![0.0], vec![0.0], vec![0.0], vec![0.0]],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![resident_neutrons], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, vec![source_atoms, 0.0, 50.0, 0.0, 0.0]).unwrap()],
        vec![ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: 1,
        }],
    )
    .unwrap()
    .with_decay_neutrons(vec![emission])
    .unwrap();
    (model, spatial)
}

// Independent nonlinear reference: source, neutrons, fuel, radioactive product,
// stable product and leaked neutrons. No production propagator is called.
fn reference(source: f64, resident: f64, steps: usize) -> [f64; 6] {
    fn rhs(x: [f64; 6]) -> [f64; 6] {
        let emission = 4.0 * x[0];
        let fission = 0.03 * x[1] * x[2];
        let product_decay = 0.1 * x[3];
        let leak = 3.0 * x[1];
        [
            -emission,
            emission - leak,
            -fission,
            2.0 * fission - product_decay,
            product_decay,
            leak,
        ]
    }
    let h = 1.0 / steps as f64;
    let mut x = [source, resident, 50.0, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|j| x[j] + 0.5 * h * a[j]));
        let c = rhs(std::array::from_fn(|j| x[j] + 0.5 * h * b[j]));
        let d = rhs(std::array::from_fn(|j| x[j] + h * c[j]));
        x = std::array::from_fn(|j| x[j] + h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0);
    }
    x
}

#[test]
fn g3_finite_fission_from_resident_or_decay_neutrons_retains_product_decay_and_refines() {
    for (source, resident) in [(0.0, 10.0), (10.0, 0.0)] {
        let expected = reference(source, resident, 2048);
        for (coarse, fine) in reference(source, resident, 1024).into_iter().zip(expected) {
            assert!((coarse - fine).abs() < 1e-9);
        }
        let mut errors = Vec::new();
        for steps in [32, 64] {
            let (model, spatial) = finite_parts(source, resident, 5.0, 0.0);
            let mut state = FusionNeutronActivationState::new(model, spatial).unwrap();
            let mut leaked = 0.0;
            for _ in 0..steps {
                let (next, report) = state
                    .trial_with_rates(
                        1.0 / steps as f64,
                        &[300.0],
                        &[0.0],
                        activation_budget(),
                        || Ok(()),
                    )
                    .unwrap();
                let counts = &report.spatial.burnup[0].reaction_counts;
                close(report.spatial.neutrons.fissions, counts[1]);
                close(report.spatial.neutrons.prompt_neutrons, counts[1]);
                close(report.decay_neutrons, counts[0]);
                close(
                    report.nuclear_generated_energy_j,
                    MEV * (3.35 * counts[0] + 5.0 * counts[1] + 0.85 * counts[2]),
                );
                close(
                    report.released.deposited_j,
                    MEV * (0.25 * counts[0] + 4.0 * counts[1] + 0.5 * counts[2]),
                );
                assert_eq!(report.kinetic.diffusion, report.spatial.neutrons);
                assert_eq!(report.external_injected_neutrons, 0.0);
                leaked += report.spatial.neutrons.leaked_neutrons;
                state = next;
            }
            let p = state.inventories()[0].populations();
            let got = [
                p[0],
                state.neutrons().neutron_counts()[0],
                p[2],
                p[3],
                p[4],
                leaked,
            ];
            close(p[1], source - p[0]);
            close(p[3] + p[4], 2.0 * (50.0 - p[2]));
            let total = MEV * (2.0 * resident + 3.35 * p[1] + 5.0 * (50.0 - p[2]) + 0.85 * p[4]);
            close(
                state.kinetic_energy_j().unwrap() + state.cell_released(0).unwrap().total_j(),
                total,
            );
            errors.push(
                got.into_iter()
                    .zip(expected)
                    .map(|(a, b)| (a - b).abs())
                    .sum::<f64>(),
            );
            if steps == 64 {
                // External forcing is absent. The finite radioactive driver
                // and resident neutrons die away; the original fission-product
                // inventory remains and continues releasing its own energy.
                for _ in 0..28 {
                    state = state
                        .trial_with_rates(0.25, &[300.0], &[0.0], activation_budget(), || Ok(()))
                        .unwrap()
                        .0;
                }
                let (_, late) = state
                    .trial_with_rates(0.5, &[300.0], &[0.0], activation_budget(), || Ok(()))
                    .unwrap();
                let counts = &late.spatial.burnup[0].reaction_counts;
                assert!(counts[2] > 1e6 * (counts[0] + counts[1]));
                close(
                    late.released.deposited_j,
                    MEV * (0.25 * counts[0] + 4.0 * counts[1] + 0.5 * counts[2]),
                );
                assert_eq!(late.external_injected_neutrons, 0.0);
            }
        }
        assert!(errors[1] < 0.4 * errors[0], "refinement errors {errors:?}");
        assert!(errors[1] < 5e-3, "refinement errors {errors:?}");
    }
}

#[test]
fn g0_g4_prompt_energy_is_explicit_unique_and_refusal_preserves_all_histories() {
    assert!(multiplying(0.0).is_err());
    let (zero_q, spatial) = finite_parts(10.0, 0.0, 0.0, 1.0);
    assert!(
        NeutronKineticModel::new(
            "no-inferred-fission-Q".into(),
            Arc::clone(zero_q.diffusion()),
            zero_q.group_energies_j().to_vec(),
            zero_q.local_heat_fractions().to_vec(),
        )
        .is_err()
    );
    assert!(matches!(
        FusionNeutronActivationState::new(zero_q, spatial),
        Err(FusionNeutronActivationError::Invalid(_))
    ));
    let (model, spatial) = finite_parts(10.0, 0.0, 5.0, 1.0);
    assert!(FusionNeutronActivationState::new(model, spatial).is_err());

    let mut delayed = data(&[2.0], 1);
    delayed.delayed_groups.push(DelayedGroup {
        data_id: "unfunded-abstract-precursor".into(),
        decay_constant_s_inv: 1.0,
        spectrum: vec![1.0],
    });
    delayed.reactions.push(reaction(
        0,
        0,
        0.1,
        DiffusionReactionKind::Fission {
            prompt_yields: vec![1.0],
            delayed_yields: vec![0.1],
        },
    ));
    assert!(
        NeutronKineticModel::new_with_reaction_energies(
            "no-inferred-precursor-storage".into(),
            Arc::new(DiffusionModel::new(delayed, DiffusionLimits::default()).unwrap()),
            vec![2.0 * MEV],
            vec![1.0],
            vec![NeutronReactionEnergy {
                nuclear_q_j: 5.0 * MEV,
                threshold_j: 0.0
            }],
        )
        .is_err()
    );

    let (model, spatial) = finite_parts(10.0, 0.0, 5.0, 0.0);
    let initial = FusionNeutronActivationState::new(model, spatial).unwrap();
    let (accepted, _) = initial
        .trial_with_rates(0.2, &[300.0], &[0.0], activation_budget(), || Ok(()))
        .unwrap();
    let before = accepted.clone();
    let mut polls = 0;
    accepted
        .trial_with_rates(0.3, &[300.0], &[0.0], activation_budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(
        accepted
            .trial_with_rates(0.3, &[300.0], &[0.0], activation_budget(), || {
                seen += 1;
                if seen == polls - 1 {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(
        accepted.neutrons().neutron_counts(),
        before.neutrons().neutron_counts()
    );
    assert_eq!(
        accepted.neutrons().precursor_counts(),
        before.neutrons().precursor_counts()
    );
    assert_eq!(accepted.inventories(), before.inventories());
    assert_eq!(accepted.cell_released(0), before.cell_released(0));
    assert_eq!(
        accepted.cell_nuclear_released(0),
        before.cell_nuclear_released(0)
    );
    assert_eq!(
        accepted.cell_prompt_nuclear_energy_j(0),
        before.cell_prompt_nuclear_energy_j(0)
    );
}
