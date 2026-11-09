//! Finite induced-fission precursor ownership and future-neutron energy.
//! Manufactured rates and energies isolate coupling, not isotope calibration.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct,
    BurnupReaction,
};
use crate::diffusion::burnup::ReactionBinding;
use crate::diffusion::{
    DelayedGroup, DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionReaction,
    DiffusionReactionKind,
};
use crate::fusion_neutrons::{InducedFissionPrecursorFunding, NeutronReactionEnergy};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance}",
    );
}

fn budget() -> FusionNeutronActivationBudget {
    FusionNeutronActivationBudget {
        spatial: DiffusionBurnupBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-8,
            },
            burnup: BurnupBudget {
                limits: BurnupLimits::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-8,
            },
            max_iterations: 96,
            reaction_tolerance: 1e-11,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-7,
    }
}

fn parts(
    parent_count: f64,
    retained_fragment: bool,
    unrelated_population: f64,
    microscopic: f64,
) -> (Arc<NeutronKineticModel>, DiffusionBurnupState) {
    let groups = Arc::new(
        NeutronGroups::new(
            "manufactured-induced-family-grid".into(),
            vec![0.0, 1.5, 3.0],
        )
        .unwrap(),
    );
    let nuclide = |z, a, lambda| BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(lambda)).unwrap(),
    };
    let inventory_model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-finite-fuel-unrepresented-fission-products".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    nuclide(10, 30, 0.0),
                    nuclide(4, 10, 0.0),
                    nuclide(1, 3, 0.5),
                    nuclide(2, 3, 0.0),
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-finite-induced-fission".into(),
                        parent: 0,
                        rate: BurnupRate::Fission {
                            cross_sections_m2: vec![0.0, microscopic],
                        },
                        products: if retained_fragment {
                            vec![BurnupProduct {
                                nuclide: 1,
                                yield_per_event: 2.0,
                            }]
                        } else {
                            vec![]
                        },
                        energy_per_event: EnergyChannels::default(),
                    },
                    BurnupReaction {
                        data_id: "manufactured-unrelated-radioactive-decay".into(),
                        parent: 2,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 3,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 0.2,
                            radiation_j: 0.3,
                            escaped_j: 0.1,
                        },
                    },
                ],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let mut reactions = (0..2)
        .map(|incident_group| DiffusionReaction {
            data_id: format!("manufactured-bound-fission-group-{incident_group}"),
            cell: 0,
            incident_group,
            macroscopic_cross_section_m_inv: if incident_group == 0 {
                0.0
            } else {
                parent_count * microscopic / 0.01
            },
            kind: DiffusionReactionKind::Fission {
                prompt_yields: vec![0.0, 1.0],
                delayed_yields: vec![0.5],
            },
        })
        .collect::<Vec<_>>();
    reactions.push(DiffusionReaction {
        data_id: "manufactured-slow-neutron-background-capture".into(),
        cell: 0,
        incident_group: 0,
        macroscopic_cross_section_m_inv: 0.3,
        kind: DiffusionReactionKind::Capture,
    });
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-induced-family-transport".into(),
                groups,
                speeds_m_s: vec![1.0; 2],
                cell_volumes_m3: vec![0.01],
                delayed_groups: vec![DelayedGroup {
                    data_id: "manufactured-induced-family-lifetime".into(),
                    decay_constant_s_inv: 1.0,
                    spectrum: vec![1.0, 0.0],
                }],
                reactions,
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let owner = Arc::new(
        NeutronKineticModel::new_with_fission_precursor_funding(
            "manufactured-prompt-and-future-neutron-conversion".into(),
            Arc::clone(&diffusion),
            vec![1.0, 2.0],
            vec![0.25, 0.25, 0.5],
            vec![
                NeutronReactionEnergy {
                    nuclear_q_j: 4.0,
                    threshold_j: 0.0,
                },
                NeutronReactionEnergy {
                    nuclear_q_j: 4.0,
                    threshold_j: 0.0,
                },
                NeutronReactionEnergy::default(),
            ],
            vec![
                Some(
                    InducedFissionPrecursorFunding::AdditionalDelayedConversion {
                        maximum_energy_j_per_fission: 0.6,
                    },
                ),
                Some(
                    InducedFissionPrecursorFunding::AdditionalDelayedConversion {
                        maximum_energy_j_per_fission: 0.6,
                    },
                ),
                None,
            ],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0, 10.0], vec![0.0]).unwrap(),
        vec![
            BurnupInventory::new(
                inventory_model,
                vec![parent_count, 0.0, unrelated_population, 0.0],
            )
            .unwrap(),
        ],
        vec![
            ReactionBinding {
                diffusion_reaction: 0,
                burnup_reaction: 0,
            },
            ReactionBinding {
                diffusion_reaction: 1,
                burnup_reaction: 0,
            },
        ],
    )
    .unwrap();
    (owner, spatial)
}

/// Shared finite source for the original material transaction. The fast group
/// replaces each fission's incident neutron, so its ten-neutron stock remains
/// constant. Fuel has lambda=.2/s; family yield=.5, lambda=1/s, slow capture=.3/s.
/// Group energies=[1,2]J, promptQ=4J, extra future-neutron conversion=.5J/event.
/// Volume=.01m3, fission heat fraction=.25, capture heat fraction=.5.
pub(crate) fn material_induced_precursor_fixture(
    parent_count: f64,
) -> FusionNeutronActivationState {
    let (model, spatial) = parts(parent_count, false, 0.0, 0.0002);
    FusionNeutronActivationState::new_with_unrepresented_fission_precursors(model, spatial).unwrap()
}

fn continuous_history(parent_count: f64, time: f64) -> (f64, f64, f64) {
    let lambda = 0.2;
    let family = 1.0;
    let capture = 0.3;
    let response = |rate: f64| ((-rate * time).exp() - (-capture * time).exp()) / (capture - rate);
    let parent = parent_count * (-lambda * time).exp();
    let precursor = 0.5 * parent_count * lambda * ((-lambda * time).exp() - (-family * time).exp())
        / (family - lambda);
    let slow = 0.5 * parent_count * lambda * family / (family - lambda)
        * (response(lambda) - response(family));
    (parent, precursor, slow)
}

#[test]
fn g1_g3_finite_induced_families_reserve_energy_and_refine_the_original_history() {
    let initial = material_induced_precursor_fixture(100.0);
    let expected = continuous_history(100.0, 1.0);
    let mut errors = Vec::new();
    for intervals in [32, 64] {
        let mut state = initial.clone();
        let mut fissions = 0.0;
        let mut generated = 0.0;
        for _ in 0..intervals {
            let (next, report) = state
                .trial_with_rates(
                    1.0 / intervals as f64,
                    &[300.0],
                    &[0.0; 2],
                    budget(),
                    || Ok(()),
                )
                .unwrap();
            let count = report.spatial.neutrons.fissions;
            close(count, report.spatial.burnup[0].fissions, 1e-9);
            close(report.kinetic.nuclear_reaction_energy_j, 4.0 * count, 1e-9);
            close(report.kinetic.induced_precursor_energy_j, 0.5 * count, 1e-9);
            close(
                report.spatial.neutrons.delayed_precursors_produced,
                0.5 * count,
                1e-9,
            );
            assert_eq!(report.kinetic.injected_precursor_energy_j, 0.0);
            close(report.nuclear_generated_energy_j, 4.5 * count, 1e-9);
            assert_eq!(report.external_injected_neutrons, 0.0);
            assert_eq!(report.external_injected_precursors, 0.0);
            assert_eq!(report.decay_precursors, 0.0);
            assert_eq!(report.reserved_decay_precursor_energy_j, 0.0);
            assert!(report.energy_roundoff_j <= budget().energy_tolerance_j);
            fissions += count;
            generated += report.nuclear_generated_energy_j;
            state = next;
        }
        let parent = state.inventories()[0].populations()[0];
        let precursor = state.neutrons().precursor_counts()[0];
        let slow = state.neutrons().neutron_counts()[0];
        close(parent, expected.0, 1e-8);
        close(state.neutrons().neutron_counts()[1], 10.0, 1e-9);
        close(fissions, 100.0 - parent, 1e-8);
        close(
            state.cumulative_induced_precursor_energy_j().unwrap(),
            0.5 * fissions,
            1e-8,
        );
        close(
            state.cell_induced_precursor_energy_j(0).unwrap(),
            0.5 * fissions,
            1e-8,
        );
        let captures = 0.5 * fissions - precursor - slow;
        let released = state.cell_released(0).unwrap();
        close(released.deposited_j, fissions + 0.5 * captures, 1e-8);
        close(released.radiation_j, 3.0 * fissions + 0.5 * captures, 1e-8);
        assert_eq!(released.escaped_j, 0.0);
        close(
            state.stored_energy_j().unwrap() + released.total_j(),
            20.0 + generated,
            1e-8,
        );
        // Fuel/flux interval means are the existing finite coupling closure.
        // Halving physical intervals must recover the actual continuous
        // parent -> precursor -> captured-neutron history at second order.
        errors.push((precursor - expected.1).abs() + (slow - expected.2).abs());
    }
    assert!(
        errors[1] < 0.3 * errors[0],
        "physical refinement errors {errors:?}"
    );
    assert!(errors[1] < 2e-4, "physical refinement errors {errors:?}");
}

#[test]
fn g0_aggregate_product_boundary_is_explicit_and_applies_to_dormant_channels() {
    let (model, spatial) = parts(100.0, false, 0.0, 0.0002);
    assert!(FusionNeutronActivationState::new(Arc::clone(&model), spatial.clone()).is_err());
    assert!(
        FusionNeutronActivationState::new_with_unrepresented_fission_precursors(model, spatial)
            .is_ok()
    );
    for cross_section in [0.0, 0.0002] {
        let (model, spatial) = parts(100.0, true, 0.0, cross_section);
        assert!(matches!(
            FusionNeutronActivationState::new_with_unrepresented_fission_precursors(model, spatial),
            Err(FusionNeutronActivationError::Invalid(
                "aggregate induced-fission families require every positive retained fission product to be unrepresented"
            )),
        ));
    }
    // Unbound fixed-background fission has no finite reaction whose product
    // boundary can be inspected. It therefore needs the same explicit boundary.
    let (model, spatial) = parts(0.0, false, 0.0, 0.0);
    let empty = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-background-only-inventory".into(),
                groups: Arc::clone(model.diffusion().groups()),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![spatial.inventories()[0].model().nuclides()[0].clone()],
                reactions: vec![],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let background = DiffusionBurnupState::new(
        spatial.neutrons().clone(),
        vec![BurnupInventory::new(empty, vec![0.0]).unwrap()],
        vec![],
    )
    .unwrap();
    assert!(FusionNeutronActivationState::new(Arc::clone(&model), background.clone()).is_err());
    assert!(
        FusionNeutronActivationState::new_with_unrepresented_fission_precursors(model, background)
            .is_ok()
    );
}

#[test]
fn g1_unrelated_radioactive_stock_keeps_its_own_energy_while_induced_families_age() {
    let (model, spatial) = parts(100.0, false, 5.0, 0.0002);
    let initial =
        FusionNeutronActivationState::new_with_unrepresented_fission_precursors(model, spatial)
            .unwrap();
    let time = 0.8_f64;
    let (next, report) = initial
        .trial_with_rates(time, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let decays = 5.0 * -(-0.5 * time).exp_m1();
    close(next.inventories()[0].populations()[2], 5.0 - decays, 1e-9);
    close(next.inventories()[0].populations()[3], decays, 1e-9);
    close(report.nuclear_released.deposited_j, 0.2 * decays, 1e-9);
    close(report.nuclear_released.radiation_j, 0.3 * decays, 1e-9);
    close(report.nuclear_released.escaped_j, 0.1 * decays, 1e-9);
    close(
        report.nuclear_generated_energy_j,
        4.5 * report.spatial.neutrons.fissions + 0.6 * decays,
        1e-8,
    );
    close(
        next.stored_energy_j().unwrap() + report.released.total_j(),
        20.0 + report.nuclear_generated_energy_j,
        1e-8,
    );
    assert_eq!(report.external_injected_precursors, 0.0);
    assert!(report.kinetic.induced_precursor_energy_j > 0.0);
}

#[test]
fn g4_induced_family_continuation_external_feed_and_late_cancel_remain_atomic() {
    let initial = material_induced_precursor_fixture(100.0);
    let (middle, first) = initial
        .trial_with_sources(0.25, &[300.0], &[0.0; 2], &[2.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(first.external_injected_precursors, 0.5);
    close(first.external_injected_precursor_energy_j, 0.5, 1e-12);
    assert!(first.kinetic.induced_precursor_energy_j > 0.0);
    let before_populations = middle.inventories()[0].populations().to_vec();
    let before_precursors = middle.neutrons().precursor_counts().to_vec();
    let before_energy = middle.cumulative_induced_precursor_energy_j().unwrap();
    let mut polls = 0;
    let (accepted, second) = middle
        .trial_with_rates(0.25, &[300.0], &[0.0; 2], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(
        middle
            .trial_with_rates(0.25, &[300.0], &[0.0; 2], budget(), || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(middle.inventories()[0].populations(), before_populations);
    assert_eq!(middle.neutrons().precursor_counts(), before_precursors);
    assert_eq!(
        middle.cumulative_induced_precursor_energy_j().unwrap(),
        before_energy
    );
    let (retry, _) = middle
        .trial_with_rates(0.25, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    close(
        retry.stored_energy_j().unwrap(),
        accepted.stored_energy_j().unwrap(),
        1e-10,
    );
    close(
        accepted.cumulative_induced_precursor_energy_j().unwrap(),
        first.kinetic.induced_precursor_energy_j + second.kinetic.induced_precursor_energy_j,
        1e-9,
    );
    close(
        accepted.stored_energy_j().unwrap() + accepted.cell_released(0).unwrap().total_j(),
        20.5 + first.nuclear_generated_energy_j + second.nuclear_generated_energy_j,
        1e-8,
    );
}
