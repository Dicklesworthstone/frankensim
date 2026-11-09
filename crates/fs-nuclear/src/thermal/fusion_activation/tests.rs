//! Manufactured event energies and cross sections isolate physical coupling;
//! these fixtures are not evaluated isotope or blanket data.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState, ReactionBinding};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionLink, DiffusionModel,
    DiffusionReaction, DiffusionReactionKind, DiffusionState,
};
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationState,
};
use crate::fusion_neutrons::history::{
    FusionNeutronHistory, FusionNeutronHistoryEntry, FusionNeutronHistoryLimits,
};
use crate::fusion_neutrons::{FusionNeutronSource, NeutronKineticModel};
use fs_qty::{Amount, Energy, Frequency, Time};
use fs_thermochem::decay::DecayLawV1;
use fs_thermochem::fusion::neutrons::{FusionNeutronEmissionV1, FusionNeutronIntervalV1};

#[path = "adaptive_tests.rs"]
mod adaptive_tests;

#[path = "decay_tests.rs"]
mod decay_tests;

#[path = "fission_tests.rs"]
mod fission_tests;

#[path = "evaluated_sf_tests.rs"]
mod evaluated_sf_tests;

fn activation_budget() -> FusionNeutronActivationBudget {
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

fn activation_state(populations: &[[f64; 3]]) -> Arc<FusionNeutronActivationState> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-activation-heat-grid".into(), vec![0.0, 8.0]).unwrap(),
    );
    let nuclide = |z, a, lambda| BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(lambda)).unwrap(),
    };
    let target = nuclide(10, 25, 0.0);
    let capture = Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: "manufactured-activation-heat-temperature-capture".into(),
                nuclide: target.id,
                reaction: NeutronReaction::Capture,
                groups: Arc::clone(&groups),
                temperatures_k: vec![300.0, 600.0],
                cross_sections_m2: vec![vec![1e-4], vec![2e-4]],
            },
            CrossSectionLimits::default(),
        )
        .unwrap(),
    );
    let inventory_model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-capture-products-and-delayed-heat".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![target, nuclide(10, 26, 0.4), nuclide(11, 26, 0.0)],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-capture-q-excluding-kinetic".into(),
                        parent: 0,
                        rate: BurnupRate::TabulatedCapture {
                            cross_sections: capture,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 4.0,
                            radiation_j: 0.5,
                            escaped_j: 0.5,
                        },
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
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            radiation_j: 1.0,
                            escaped_j: 0.5,
                        },
                    },
                ],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let count = populations.len();
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-finite-activation-heat-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![0.01; count],
                delayed_groups: vec![],
                reactions: populations
                    .iter()
                    .enumerate()
                    .map(|(cell, population)| DiffusionReaction {
                        data_id: format!("manufactured-bound-capture-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 0.01 * population[0],
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: if count == 2 {
                    vec![DiffusionLink {
                        cell_a: 0,
                        cell_b: 1,
                        conductance_m2: vec![0.002],
                    }]
                } else {
                    vec![]
                },
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-activation-heat-partition".into(),
            Arc::clone(&diffusion),
            vec![4.0],
            vec![0.25; count],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0; count], vec![]).unwrap(),
        populations
            .iter()
            .map(|population| {
                BurnupInventory::new(Arc::clone(&inventory_model), population.to_vec()).unwrap()
            })
            .collect(),
        (0..count)
            .map(|cell| ReactionBinding {
                diffusion_reaction: cell,
                burnup_reaction: 0,
            })
            .collect(),
    )
    .unwrap();
    Arc::new(FusionNeutronActivationState::new(model, spatial).unwrap())
}

fn emission(duration: f64, count: f64) -> FusionNeutronEmissionV1 {
    let amount = Amount::new(count / 6.022_140_76e23);
    FusionNeutronEmissionV1::new(
        Time::new(duration),
        [amount, Amount::new(0.0)],
        [
            Energy::new(4.0 * (amount.value() * 6.022_140_76e23)),
            Energy::new(0.0),
        ],
    )
    .unwrap()
}

fn history(state: &FusionNeutronActivationState) -> Arc<FusionNeutronHistory> {
    let entries = [
        (0, 0.0, 0.25, 4.0),
        (0, 0.25, 0.5, 8.0),
        (1, 0.125, 0.375, 2.0),
    ]
    .into_iter()
    .map(|(cell, start, end, count)| FusionNeutronHistoryEntry {
        cell,
        interval: FusionNeutronIntervalV1::new(
            Time::new(start),
            Time::new(end),
            emission(end - start, count),
        )
        .unwrap(),
    })
    .collect();
    Arc::new(
        FusionNeutronHistory::project(
            Arc::clone(state.model()),
            entries,
            FusionNeutronHistoryLimits::default(),
        )
        .unwrap(),
    )
}

fn material_cells() -> DepositionSimulation<'static> {
    thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0), cell(10, 2, 3.0, 900.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    )
}
fn bindings() -> Vec<(CellId, MaterialId)> {
    vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))]
}

#[test]
fn g1_history_activation_heats_each_original_material_from_one_native_trial() {
    let initial = material_cells();
    let state = activation_state(&[[100.0, 0.0, 0.0], [80.0, 0.0, 0.0]]);
    let history = history(&state);
    let sources = FusionNeutronActivationCellSource::bind_with_history(
        Arc::clone(&state),
        Time::new(0.0),
        bindings(),
        Arc::clone(&history),
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (hot, accepted, _) = heating
        .advance_to(&initial, 0.625, controls(0.125), || Ok(()))
        .unwrap();

    // Direct original-owner evolution uses physical cell order and each
    // material's distinct initial temperature and heat capacity. Thermal IDs
    // deliberately sort in the opposite order.
    let mut direct = (*state).clone();
    let mut temperatures = [300.0, 450.0];
    let mut captures = 0.0;
    let mut decays = 0.0;
    let mut defect_j = 0.0;
    for interval in 0..5 {
        let time = interval as f64 * 0.125;
        let (next, report) = direct
            .trial_with_rates(
                0.125,
                &temperatures,
                history.rates_at(Time::new(time)).unwrap(),
                activation_budget(),
                || Ok(()),
            )
            .unwrap();
        for nuclear in &report.spatial.burnup {
            captures += nuclear.captures;
            decays += nuclear.decays;
        }
        defect_j += report.energy_roundoff_j;
        temperatures[0] = 300.0 + next.cell_released(0).unwrap().deposited_j / 8.0;
        temperatures[1] = 450.0 + next.cell_released(1).unwrap().deposited_j / 6.0;
        direct = next;
    }
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    assert!(Arc::ptr_eq(left.forcing().history().unwrap(), &history));
    let last = accepted.last_step(CellId(20)).unwrap();
    assert!(core::ptr::eq(
        last.activation(),
        accepted.last_step(CellId(10)).unwrap().activation()
    ));
    assert_eq!(last.start(), Time::new(0.5));
    assert_eq!(last.activation().spatial.neutrons.injected_neutrons, 0.0);
    assert!(last.burnup().decays > 0.0);
    assert_eq!(
        last.activation().kinetic.diffusion,
        last.activation().spatial.neutrons
    );
    for index in 0..2 {
        close(
            left.state().neutrons().neutron_counts()[index],
            direct.neutrons().neutron_counts()[index],
            1e-9,
        );
        for (actual, expected) in left.state().inventories()[index]
            .populations()
            .iter()
            .zip(direct.inventories()[index].populations())
        {
            close(*actual, *expected, 1e-9);
        }
        assert!(left.state().inventories()[index].populations()[0] < [100.0, 80.0][index]);
    }
    close(
        hot.temperature_k(CellId(20)).unwrap(),
        temperatures[0],
        1e-9,
    );
    close(
        hot.temperature_k(CellId(10)).unwrap(),
        temperatures[1],
        1e-9,
    );
    let released = accepted.cumulative_released().unwrap();
    close(
        hot.energy_j() - initial.energy_j(),
        released.deposited_j,
        1e-8,
    );
    close(
        accepted.audit().thermal_credit_j,
        released.deposited_j,
        1e-8,
    );
    close(
        accepted.audit().radiation_escaped_j,
        released.radiation_j,
        1e-8,
    );
    close(accepted.audit().captures, captures, 1e-9);
    close(accepted.audit().decays, decays, 1e-9);
    // Shared defects are charged once, not once per material handle.
    close(accepted.audit().nuclear_energy_defect_j, defect_j, 1e-10);
    close(
        left.state().kinetic_energy_j().unwrap() + released.total_j(),
        56.0 + 5.0 * captures + 3.5 * decays,
        1e-7,
    );
    assert_eq!(accepted.audit().fissions, 0.0);
    assert_eq!(accepted.audit().produced_nuclei, 0.0);
    assert_eq!(accepted.audit().substeps, 5);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(state.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
}

#[test]
fn g1_g3_activated_products_heat_after_neutrons_are_gone_and_restart_without_replay() {
    let initial = thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let state = activation_state(&[[40.0, 12.0, 5.0]]);
    let pulse = Arc::new(
        FusionNeutronPulse::new(
            Time::new(0.0),
            FusionNeutronSource::project_interval_mean(
                Arc::clone(state.model()),
                emission(0.25, 0.0),
                0,
            )
            .unwrap(),
        )
        .unwrap(),
    );
    let sources = FusionNeutronActivationCellSource::bind(
        Arc::clone(&state),
        Time::new(0.0),
        vec![(CellId(20), MaterialId(1))],
        Arc::clone(&pulse),
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (whole, end, _) = heating
        .advance_to(&initial, 2.5, controls(0.5), || Ok(()))
        .unwrap();
    let remaining = 12.0 * (-0.4_f64 * 2.5).exp();
    let decays = 12.0 - remaining;
    close(
        end.source(CellId(20)).unwrap().state().inventories()[0].populations()[1],
        remaining,
        1e-10,
    );
    close(whole.energy_j() - initial.energy_j(), 2.0 * decays, 1e-8);
    close(end.audit().decays, decays, 1e-10);
    assert_eq!(end.audit().captures, 0.0);
    assert_eq!(
        end.last_step(CellId(20))
            .unwrap()
            .activation()
            .kinetic
            .released,
        EnergyChannels::default()
    );
    let (middle, running, _) = heating
        .advance_to(&initial, 1.25, controls(0.5), || Ok(()))
        .unwrap();
    let restart_sources = FusionNeutronActivationCellSource::bind(
        Arc::clone(running.source(CellId(20)).unwrap().state()),
        Time::new(middle.time_s()),
        vec![(CellId(20), MaterialId(1))],
        pulse,
    )
    .unwrap();
    let rebound =
        FusionNeutronActivationHeating::new(&middle, restart_sources, activation_budget()).unwrap();
    let (split, continued, _) = rebound
        .advance_to(&middle, 2.5, controls(0.5), || Ok(()))
        .unwrap();
    close(split.energy_j(), whole.energy_j(), 1e-8);
    close(
        continued.cumulative_released().unwrap().total_j(),
        end.cumulative_released().unwrap().total_j(),
        1e-8,
    );
    close(
        continued.audit().thermal_credit_j,
        split.energy_j() - middle.energy_j(),
        1e-8,
    );
    assert!(continued.audit().decays < end.audit().decays);
    assert_eq!(rebound.time_s(), 1.25);
    assert_eq!(state.neutrons().neutron_counts(), &[0.0]);
}

#[test]
fn g0_g4_activation_requires_complete_cohort_clock_and_atomic_final_publication() {
    let initial = material_cells();
    let state = activation_state(&[[100.0, 0.0, 0.0], [80.0, 0.0, 0.0]]);
    let history = history(&state);
    let make = |time, schedule| {
        FusionNeutronActivationCellSource::bind_with_history(
            Arc::clone(&state),
            Time::new(time),
            bindings(),
            schedule,
        )
        .unwrap()
    };
    let sources = make(0.0, Arc::clone(&history));
    let other = make(0.0, Arc::new((*history).clone()));
    assert!(
        FusionNeutronActivationHeating::new(
            &initial,
            vec![sources[0].clone(), other[1].clone()],
            activation_budget()
        )
        .is_err()
    );
    assert!(
        FusionNeutronActivationHeating::new(&initial, sources[..1].to_vec(), activation_budget())
            .is_err()
    );
    assert!(
        FusionNeutronActivationHeating::new(
            &initial,
            make(0.125, Arc::clone(&history)),
            activation_budget()
        )
        .is_err()
    );
    assert!(
        FusionNeutronActivationCellSource::bind_with_history(
            Arc::clone(&state),
            Time::new(0.0),
            vec![bindings()[0]; 2],
            Arc::clone(&history),
        )
        .is_err()
    );
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let mut successful_polls = 0;
    heating
        .advance_to(&initial, 0.625, controls(0.125), || {
            successful_polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(successful_polls > 1);
    let mut polls = 0;
    let refused = heating.advance_to(&initial, 0.625, controls(0.125), || {
        polls += 1;
        if polls == successful_polls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refused,
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(20)).unwrap(), 300.0);
    assert_eq!(initial.temperature_k(CellId(10)).unwrap(), 450.0);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
    assert!(Arc::ptr_eq(
        heating.source(CellId(20)).unwrap().state(),
        &state
    ));
    assert_eq!(state.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
}
