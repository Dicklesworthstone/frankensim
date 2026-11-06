//! Manufactured coefficients exercise the coupled equations and extensive
//! accounts; none of these cross sections or event energies are isotope data.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct,
    BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::diffusion::burnup::ReactionBinding;
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionReaction, DiffusionReactionKind,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::{Amount, Energy, Frequency, Time};
use fs_thermochem::decay::DecayLawV1;
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL;
use fs_thermochem::fusion::neutrons::FusionNeutronEmissionV1;

const NEUTRON_ENERGY: f64 = 4.0;
const CAPTURE_Q: EnergyChannels = EnergyChannels {
    deposited_j: 4.0,
    radiation_j: 0.5,
    escaped_j: 0.5,
};
const DECAY_Q: EnergyChannels = EnergyChannels {
    deposited_j: 2.0,
    radiation_j: 1.0,
    escaped_j: 0.5,
};

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
            max_iterations: 64,
            reaction_tolerance: 1e-11,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-7,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn nuclide(z: u16, a: u16, lambda: f64) -> BurnupNuclide {
    BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(lambda)).unwrap(),
    }
}

fn fixture(
    neutrons: f64,
    populations: [f64; 3],
    temperature_table: bool,
) -> (Arc<NeutronKineticModel>, DiffusionBurnupState) {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-activation-group".into(), vec![0.0, 8.0]).unwrap(),
    );
    let target = nuclide(10, 25, 0.0);
    let capture = if temperature_table {
        BurnupRate::TabulatedCapture {
            cross_sections: Arc::new(
                CrossSectionTable::new(
                    CrossSectionData {
                        data_id: "manufactured-activation-temperature-table".into(),
                        nuclide: target.id,
                        reaction: NeutronReaction::Capture,
                        groups: Arc::clone(&groups),
                        temperatures_k: vec![300.0, 600.0],
                        cross_sections_m2: vec![vec![0.01], vec![0.02]],
                    },
                    CrossSectionLimits::default(),
                )
                .unwrap(),
            ),
        }
    } else {
        BurnupRate::Capture {
            cross_sections_m2: vec![0.01],
        }
    };
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-target-activation-decay".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![target, nuclide(10, 26, 0.4), nuclide(11, 26, 0.0)],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-capture-q-excludes-incident-kinetic".into(),
                        parent: 0,
                        rate: capture,
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: CAPTURE_Q,
                    },
                    BurnupReaction {
                        data_id: "manufactured-activated-product-decay".into(),
                        parent: 1,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: DECAY_Q,
                    },
                ],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-activation-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                delayed_groups: vec![],
                reactions: vec![
                    DiffusionReaction {
                        data_id: "manufactured-finite-target-capture".into(),
                        cell: 0,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 0.01 * populations[0],
                        kind: DiffusionReactionKind::Capture,
                    },
                    DiffusionReaction {
                        data_id: "manufactured-fixed-within-group-scatter".into(),
                        cell: 0,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 0.2,
                        kind: DiffusionReactionKind::Scatter { outgoing_group: 0 },
                    },
                ],
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let kinetic = Arc::new(
        NeutronKineticModel::new(
            "manufactured-activation-kinetic-partition".into(),
            Arc::clone(&diffusion),
            vec![NEUTRON_ENERGY],
            vec![0.25, 1.0],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![neutrons], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, populations.to_vec()).unwrap()],
        vec![ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: 0,
        }],
    )
    .unwrap();
    (kinetic, spatial)
}

fn state(neutrons: f64, populations: [f64; 3]) -> FusionNeutronActivationState {
    let (model, spatial) = fixture(neutrons, populations, false);
    FusionNeutronActivationState::new(model, spatial).unwrap()
}

fn assert_energy_accounts(report: &FusionNeutronActivationStep) {
    assert_eq!(report.kinetic.diffusion, report.spatial.neutrons);
    let nuclear = &report.spatial.burnup[0];
    let captured_neutrons = report.spatial.neutrons.captures;
    close(report.kinetic.released.deposited_j, captured_neutrons, 1e-9);
    close(
        report.kinetic.released.radiation_j,
        3.0 * captured_neutrons,
        1e-9,
    );
    close(
        report.nuclear_released.deposited_j,
        4.0 * nuclear.captures + 2.0 * nuclear.decays,
        1e-9,
    );
    close(
        report.nuclear_generated_energy_j,
        5.0 * nuclear.captures + 3.5 * nuclear.decays,
        1e-9,
    );
    close(
        report.released.total_j(),
        report.kinetic.released.total_j() + report.nuclear_released.total_j(),
        1e-9,
    );
    close(
        report.kinetic.final_kinetic_energy_j + report.released.total_j(),
        report.kinetic.initial_kinetic_energy_j
            + report.kinetic.injected_kinetic_energy_j
            + report.nuclear_generated_energy_j,
        1e-8,
    );
    assert!(report.energy_roundoff_j <= budget().energy_tolerance_j);
    assert!(report.spatial.reaction_defect <= budget().spatial.reaction_tolerance);
}

#[test]
fn finite_capture_uses_the_same_native_history_and_does_not_double_count_energy() {
    let initial = state(20.0, [100.0, 0.0, 0.0]);
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![0.0; 3]],
    };
    let (native, native_report) = initial
        .spatial()
        .trial_with_sources(0.6, &[300.0], &sources, budget().spatial, || Ok(()))
        .unwrap();
    let (next, report) = initial
        .trial_with_rates(0.6, &[300.0], &[0.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.spatial.neutrons, native_report.neutrons);
    assert_eq!(report.spatial.burnup, native_report.burnup);
    assert_eq!(report.spatial.iterations, native_report.iterations);
    assert_eq!(
        next.neutrons().neutron_counts(),
        native.neutrons().neutron_counts()
    );
    assert_eq!(
        next.inventories()[0].populations(),
        native.inventories()[0].populations()
    );
    assert!(report.spatial.iterations > 1);
    let captures = report.spatial.burnup[0].captures;
    close(
        100.0 - next.inventories()[0].populations()[0],
        captures,
        1e-9,
    );
    close(20.0 - next.neutrons().neutron_counts()[0], captures, 1e-9);
    close(
        next.inventories()[0].populations()[1] + next.inventories()[0].populations()[2],
        captures,
        1e-9,
    );
    assert!(report.spatial.burnup[0].decays > 0.0);
    assert!(next.neutrons().model().reactions()[0].macroscopic_cross_section_m_inv < 1.0);
    assert_eq!(next.model(), initial.model());
    assert_eq!(next.cell_released(0), Some(report.released));
    assert_eq!(next.cell_nuclear_released(0), Some(report.nuclear_released));
    assert_eq!(next.cell_kinetic_released(0), Some(report.kinetic.released));
    assert_energy_accounts(&report);
    assert_eq!(initial.neutrons().neutron_counts(), &[20.0]);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
}

#[test]
fn finite_depletion_converges_to_the_independent_coupled_capture_solution() {
    // dn/dt = dP/dt = -k*n*P; P-n = delta. The exact solution is
    // n(t) = delta*n0*exp(-k*delta*t)/(P0-n0*exp(-k*delta*t)).
    let time = 0.6;
    let decay = (-0.01_f64 * 80.0 * time).exp();
    let expected = 80.0 * 20.0 * decay / (100.0 - 20.0 * decay);
    let evolve = |steps: usize| {
        let mut current = state(20.0, [100.0, 0.0, 0.0]);
        for _ in 0..steps {
            let (next, report) = current
                .trial_with_rates(time / steps as f64, &[300.0], &[0.0], budget(), || Ok(()))
                .unwrap();
            assert_energy_accounts(&report);
            current = next;
        }
        close(
            current.inventories()[0].populations()[0] - current.neutrons().neutron_counts()[0],
            80.0,
            1e-8,
        );
        current.neutrons().neutron_counts()[0]
    };
    let coarse_error = (evolve(6) - expected).abs();
    let fine_error = (evolve(24) - expected).abs();
    assert!(
        fine_error < coarse_error / 8.0,
        "coarse {coarse_error}, fine {fine_error}"
    );
    assert!(fine_error < 5e-4);
}

#[test]
fn activated_products_release_only_actual_decay_heat_after_neutrons_are_gone() {
    let initial = state(0.0, [40.0, 12.0, 5.0]);
    let duration = 2.5;
    let (cooled, report) = initial
        .trial_with_rates(duration, &[300.0], &[0.0], budget(), || Ok(()))
        .unwrap();
    let remaining = 12.0 * (-0.4_f64 * duration).exp();
    let decays = 12.0 - remaining;
    close(cooled.inventories()[0].populations()[1], remaining, 1e-10);
    close(
        cooled.inventories()[0].populations()[2],
        5.0 + decays,
        1e-10,
    );
    close(report.nuclear_generated_energy_j, 3.5 * decays, 1e-9);
    assert_eq!(report.spatial.burnup[0].captures, 0.0);
    assert_eq!(report.kinetic.released, EnergyChannels::default());
    assert_eq!(report.kinetic.injected_kinetic_energy_j, 0.0);
    assert_energy_accounts(&report);
    // Starting an energy adapter at an accepted checkpoint does not replay the
    // checkpoint's earlier nuclear heat or reset its actual product populations.
    let rebound =
        FusionNeutronActivationState::new(Arc::clone(cooled.model()), cooled.spatial().clone())
            .unwrap();
    assert_eq!(rebound.cell_released(0), Some(EnergyChannels::default()));
    assert_eq!(
        rebound.inventories()[0].populations(),
        cooled.inventories()[0].populations()
    );
    assert_eq!(
        rebound.inventories()[0].cumulative_released(),
        cooled.inventories()[0].cumulative_released()
    );
    assert!(rebound.inventories()[0].cumulative_released().deposited_j > 0.0);
}

#[test]
fn original_fusion_source_survives_depletion_and_temperature_dependent_capture() {
    let (model, spatial) = fixture(0.0, [100.0, 0.0, 0.0], true);
    let initial = FusionNeutronActivationState::new(Arc::clone(&model), spatial).unwrap();
    let amounts = [18.0, 6.0].map(|count| Amount::new(count / AVOGADRO_PER_MOL));
    let energies =
        amounts.map(|amount| Energy::new(amount.value() * AVOGADRO_PER_MOL * NEUTRON_ENERGY));
    let emission = FusionNeutronEmissionV1::new(Time::new(0.4), amounts, energies).unwrap();
    let source = FusionNeutronSource::project_interval_mean(model, emission, 0).unwrap();
    let (cold, cold_report) = initial
        .trial_with_fusion(0.4, &[300.0], &[source.clone()], budget(), || Ok(()))
        .unwrap();
    let (hot, hot_report) = initial
        .trial_with_fusion(0.4, &[600.0], &[source.clone()], budget(), || Ok(()))
        .unwrap();
    assert!(hot_report.spatial.burnup[0].captures > cold_report.spatial.burnup[0].captures);
    assert_eq!(hot_report.kinetic.fusion_sources, vec![source.clone()]);
    close(
        hot_report.kinetic.injected_kinetic_energy_j,
        emission.total_kinetic_energy().value(),
        1e-10,
    );
    close(
        hot_report.spatial.neutrons.injected_neutrons,
        emission.total_particle_count(),
        1e-10,
    );
    assert_energy_accounts(&cold_report);
    assert_energy_accounts(&hot_report);
    // Current macroscopic coefficients differ, but the original emission/model
    // still binds for the next interval. It is neither relabelled nor renormalized.
    assert_ne!(hot.neutrons().model(), hot.model().diffusion());
    let (later, report) = hot
        .trial_with_fusion(0.4, &[300.0], &[source.clone()], budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.kinetic.fusion_sources, vec![source]);
    assert_energy_accounts(&report);
    let (_, shutdown) = later
        .trial_with_fusion(0.4, &[300.0], &[], budget(), || Ok(()))
        .unwrap();
    assert_eq!(shutdown.kinetic.injected_kinetic_energy_j, 0.0);
    assert!(shutdown.spatial.burnup[0].decays > 0.0);
    assert_energy_accounts(&shutdown);
    assert!(cold.inventories()[0].populations()[0] < 100.0);
}

#[test]
fn original_topology_and_fixed_backgrounds_are_required_but_bound_capture_can_evolve() {
    let (model, spatial) = fixture(20.0, [100.0, 0.0, 0.0], false);
    for changed in 0..6 {
        let mut data = model.diffusion().data().clone();
        match changed {
            0 => data.data_id.push_str("-different"),
            1 => {
                data.groups = Arc::new(
                    NeutronGroups::new("different-group-identity".into(), vec![0.0, 8.0]).unwrap(),
                )
            }
            2 => data.speeds_m_s[0] = 2.0,
            3 => data.cell_volumes_m3[0] = 2.0,
            4 => data.reactions.swap(0, 1),
            5 => data.reactions[1].macroscopic_cross_section_m_inv = 0.3,
            _ => unreachable!(),
        }
        let changed_model = Arc::new(
            NeutronKineticModel::new(
                model.data_id().into(),
                Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
                model.group_energies_j().to_vec(),
                model.local_heat_fractions().to_vec(),
            )
            .unwrap(),
        );
        assert!(matches!(
            FusionNeutronActivationState::new(changed_model, spatial.clone()),
            Err(FusionNeutronActivationError::Invalid(_))
        ));
    }
    let mut data = model.diffusion().data().clone();
    data.reactions[0].macroscopic_cross_section_m_inv = 0.5;
    let changed_model = Arc::new(
        NeutronKineticModel::new(
            model.data_id().into(),
            Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
            model.group_energies_j().to_vec(),
            model.local_heat_fractions().to_vec(),
        )
        .unwrap(),
    );
    assert!(FusionNeutronActivationState::new(changed_model, spatial).is_ok());
}

#[test]
fn refusal_and_late_cancellation_leave_both_particle_and_energy_histories_unchanged() {
    let initial = state(20.0, [100.0, 0.0, 0.0]);
    let mut checkpoints = 0;
    let (_, successful) = initial
        .trial_with_rates(0.4, &[300.0], &[1.0], budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    assert_energy_accounts(&successful);
    let mut visited = 0;
    assert!(
        initial
            .trial_with_rates(0.4, &[300.0], &[1.0], budget(), || {
                visited += 1;
                if visited == checkpoints {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(visited, checkpoints);
    for rates in [vec![], vec![f64::NAN], vec![-1.0], vec![0.0, 0.0]] {
        assert!(
            initial
                .trial_with_rates(0.0, &[300.0], &rates, budget(), || Ok(()))
                .is_err()
        );
    }
    assert!(
        initial
            .trial_with_rates(0.0, &[299.0], &[0.0], budget(), || Ok(()))
            .is_err()
    );
    let mut invalid_budget = budget();
    invalid_budget.energy_tolerance_j = -1.0;
    assert!(
        initial
            .trial_with_rates(0.4, &[300.0], &[0.0], invalid_budget, || Ok(()))
            .is_err()
    );
    assert_eq!(initial.neutrons().neutron_counts(), &[20.0]);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
    assert_eq!(
        initial.cell_kinetic_released(0),
        Some(EnergyChannels::default())
    );
    assert_eq!(
        initial.cell_nuclear_released(0),
        Some(EnergyChannels::default())
    );
}

#[test]
fn spatial_transfer_retains_each_materials_kinetic_nuclear_and_leakage_accounts() {
    use crate::diffusion::{DiffusionLink, LeakageBoundary};

    let (single_model, single_spatial) = fixture(10.0, [100.0, 0.0, 0.0], false);
    let mut data = single_model.diffusion().data().clone();
    data.cell_volumes_m3 = vec![1.0, 2.0];
    let mut second_reactions = data.reactions.clone();
    for reaction in &mut second_reactions {
        reaction.cell = 1;
        reaction.data_id.push_str("-second-cell");
    }
    second_reactions[0].macroscopic_cross_section_m_inv = 50.0 * 0.01 / 2.0;
    data.reactions.extend(second_reactions);
    data.links = vec![DiffusionLink {
        cell_a: 0,
        cell_b: 1,
        conductance_m2: vec![0.3],
    }];
    data.leakage = vec![LeakageBoundary {
        cell: 1,
        conductance_m2: vec![0.1],
    }];
    let diffusion = Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap());
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-spatial-activation-energy".into(),
            Arc::clone(&diffusion),
            vec![NEUTRON_ENERGY],
            vec![0.25, 1.0, 0.75, 1.0],
        )
        .unwrap(),
    );
    let second = BurnupInventory::new(
        Arc::clone(single_spatial.inventories()[0].model()),
        vec![50.0, 0.0, 0.0],
    )
    .unwrap();
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![10.0, 0.0], vec![]).unwrap(),
        vec![single_spatial.inventories()[0].clone(), second],
        vec![
            ReactionBinding {
                diffusion_reaction: 0,
                burnup_reaction: 0,
            },
            ReactionBinding {
                diffusion_reaction: 2,
                burnup_reaction: 0,
            },
        ],
    )
    .unwrap();
    let initial = FusionNeutronActivationState::new(model, spatial).unwrap();
    let (next, report) = initial
        .trial_with_rates(0.5, &[300.0, 300.0], &[4.0, 0.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.spatial.neutrons, report.kinetic.diffusion);
    assert!(report.spatial.neutrons.link_a_to_b_counts[0] > 0.0);
    assert!(report.spatial.neutrons.leaked_neutrons > 0.0);
    assert!(next.inventories()[1].populations()[1] > 0.0);
    for cell in 0..2 {
        let neutrons_captured = report.spatial.neutrons.reaction_counts[2 * cell];
        let neutron_heat_fraction = if cell == 0 { 0.25 } else { 0.75 };
        let nuclear = &report.spatial.burnup[cell];
        close(
            report.cell_released[cell].deposited_j,
            neutron_heat_fraction * NEUTRON_ENERGY * neutrons_captured
                + 4.0 * nuclear.captures
                + 2.0 * nuclear.decays,
            1e-9,
        );
        assert_eq!(next.cell_released(cell), Some(report.cell_released[cell]));
        assert_eq!(next.cell_nuclear_released(cell), Some(nuclear.released));
    }
    close(
        report.kinetic.released.escaped_j,
        report.spatial.neutrons.leaked_neutrons * NEUTRON_ENERGY,
        1e-9,
    );
    close(
        report.kinetic.final_kinetic_energy_j + report.released.total_j(),
        report.kinetic.initial_kinetic_energy_j
            + report.kinetic.injected_kinetic_energy_j
            + report.nuclear_generated_energy_j,
        1e-8,
    );
    close(
        report.released.total_j(),
        report
            .cell_released
            .iter()
            .map(|energy| energy.total_j())
            .sum(),
        1e-9,
    );
    assert_eq!(next.cell_released(2), None);
    assert!(report.energy_roundoff_j <= budget().energy_tolerance_j);
}
