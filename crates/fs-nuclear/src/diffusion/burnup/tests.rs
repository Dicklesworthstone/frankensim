use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::diffusion::{DiffusionData, DiffusionLimits, DiffusionLink, DiffusionReaction};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

pub(super) fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance}"
    );
}
pub(super) fn budget() -> DiffusionBurnupBudget {
    DiffusionBurnupBudget {
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
        reaction_tolerance: 1e-10,
        relaxation: 1.0,
    }
}
fn stable(index: u16) -> BurnupNuclide {
    BurnupNuclide {
        id: NuclideId::new(10, 25 + index, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
    }
}

pub(super) fn fixture(
    cells: usize,
    fission: bool,
    temperature_table: bool,
) -> DiffusionBurnupState {
    let groups =
        Arc::new(NeutronGroups::new("manufactured-spatial-group".into(), vec![0.0, 1.0]).unwrap());
    let sigma = if temperature_table {
        let table = Arc::new(
            CrossSectionTable::new(
                CrossSectionData {
                    data_id: "manufactured-temperature-row".into(),
                    nuclide: stable(0).id,
                    reaction: if fission {
                        NeutronReaction::Fission
                    } else {
                        NeutronReaction::Capture
                    },
                    groups: Arc::clone(&groups),
                    temperatures_k: vec![300.0, 600.0],
                    cross_sections_m2: vec![vec![0.01], vec![0.02]],
                },
                CrossSectionLimits::default(),
            )
            .unwrap(),
        );
        if fission {
            BurnupRate::TabulatedFission {
                cross_sections: table,
            }
        } else {
            BurnupRate::TabulatedCapture {
                cross_sections: table,
            }
        }
    } else if fission {
        BurnupRate::Fission {
            cross_sections_m2: vec![0.01],
        }
    } else {
        BurnupRate::Capture {
            cross_sections_m2: vec![0.01],
        }
    };
    let inventory_model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-spatial-burnup".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![stable(0), stable(1)],
                reactions: vec![BurnupReaction {
                    data_id: "manufactured-neutron-removal".into(),
                    parent: 0,
                    rate: sigma,
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: if fission { 2.0 } else { 1.0 },
                    }],
                    energy_per_event: EnergyChannels {
                        deposited_j: 8.0,
                        radiation_j: 2.0,
                        escaped_j: 1.0,
                    },
                }],
            },
            budget().burnup.limits,
        )
        .unwrap(),
    );
    let inventories = (0..cells)
        .map(|_| BurnupInventory::new(Arc::clone(&inventory_model), vec![100.0, 0.0]).unwrap())
        .collect();
    let model = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-spatial-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0; cells],
                delayed_groups: vec![],
                reactions: (0..cells)
                    .map(|cell| DiffusionReaction {
                        data_id: format!("explicit-channel-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 1.0,
                        kind: if fission {
                            DiffusionReactionKind::Fission {
                                prompt_yields: vec![2.0],
                                delayed_yields: vec![],
                            }
                        } else {
                            DiffusionReactionKind::Capture
                        },
                    })
                    .collect(),
                links: if cells == 2 {
                    vec![DiffusionLink {
                        cell_a: 0,
                        cell_b: 1,
                        conductance_m2: vec![0.5],
                    }]
                } else {
                    vec![]
                },
                leakage: vec![],
            },
            budget().diffusion.limits,
        )
        .unwrap(),
    );
    let neutrons = DiffusionState::new(model, vec![20.0; cells], vec![]).unwrap();
    DiffusionBurnupState::new(
        neutrons,
        inventories,
        (0..cells)
            .map(|cell| ReactionBinding {
                diffusion_reaction: cell,
                burnup_reaction: 0,
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn g1_finite_fuel_fission_multiplication_and_prompt_heat_close_together() {
    let initial = fixture(1, true, false);
    let (next, report) = initial.trial(0.6, &[300.0], budget(), || Ok(())).unwrap();
    let fissions = report.burnup[0].fissions;
    assert!(fissions > 10.0 && fissions < 100.0);
    assert!(report.iterations > 1);
    assert!(report.reaction_defect <= budget().reaction_tolerance);
    close(report.neutrons.fissions, fissions, 1e-10);
    close(next.neutrons().neutron_counts()[0] - 20.0, fissions, 1e-9);
    close(
        100.0 - next.inventories()[0].populations()[0],
        fissions,
        1e-9,
    );
    close(next.inventories()[0].populations()[1], 2.0 * fissions, 1e-9);
    close(report.burnup[0].released.deposited_j, 8.0 * fissions, 1e-8);
    close(report.burnup[0].released.radiation_j, 2.0 * fissions, 1e-8);
    close(report.burnup[0].released.escaped_j, fissions, 1e-8);
    // Original accepted histories remain an independent physical checkpoint.
    assert_eq!(initial.neutrons().neutron_counts(), &[20.0]);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0]);
}

#[test]
fn g3_refinement_approaches_the_independent_nonlinear_fission_solution() {
    // n'=sigma*n*F, F'=-sigma*n*F; n+F=120, for prompt yield two.
    let exact_fuel = 100.0 * 120.0 / (100.0 + 20.0 * (1.2_f64 * 0.6).exp());
    let mut errors = Vec::new();
    for steps in [3, 6, 12] {
        let mut state = fixture(1, true, false);
        for _ in 0..steps {
            state
                .advance(0.6 / steps as f64, &[300.0], budget(), || Ok(()))
                .unwrap();
        }
        close(
            state.neutrons().neutron_counts()[0] + state.inventories()[0].populations()[0],
            120.0,
            2e-8,
        );
        errors.push((state.inventories()[0].populations()[0] - exact_fuel).abs());
    }
    assert!(errors[1] < errors[0] / 1.8, "{errors:?}");
    assert!(errors[2] < errors[1] / 1.8, "{errors:?}");
}

#[test]
fn g1_capture_temperature_and_spatial_current_feed_the_same_product_counts() {
    let initial = fixture(2, false, true);
    let (next, report) = initial
        .trial(0.4, &[300.0, 600.0], budget(), || Ok(()))
        .unwrap();
    let cold = report.burnup[0].captures;
    let hot = report.burnup[1].captures;
    assert!(hot > cold);
    assert!(report.neutrons.link_a_to_b_counts[0] > report.neutrons.link_b_to_a_counts[0]);
    close(report.neutrons.captures, cold + hot, 1e-9);
    close(
        40.0 - next.neutrons().neutron_counts().iter().sum::<f64>(),
        cold + hot,
        1e-9,
    );
    close(next.inventories()[0].populations()[1], cold, 1e-9);
    close(next.inventories()[1].populations()[1], hot, 1e-9);
    for (index, temperature) in [300.0, 600.0].iter().enumerate() {
        let sigma = if *temperature == 300.0 { 0.01 } else { 0.02 };
        let macro_sigma =
            next.neutrons().model().reactions()[index].macroscopic_cross_section_m_inv;
        assert!(
            macro_sigma < 100.0 * sigma
                && macro_sigma > next.inventories()[index].populations()[0] * sigma
        );
    }
}

#[test]
fn g4_unconverged_cancelled_or_out_of_domain_steps_publish_nothing() {
    let mut state = fixture(1, true, true);
    let initial_neutrons = state.neutrons().neutron_counts().to_vec();
    let initial_inventory = state.inventories()[0].clone();
    let mut strict = budget();
    strict.max_iterations = 1;
    assert!(matches!(
        state.advance(0.6, &[300.0], strict, || Ok(())),
        Err(DiffusionBurnupError::Nonconvergence { .. })
    ));
    assert!(state.advance(0.1, &[600.1], budget(), || Ok(())).is_err());
    assert!(matches!(
        state.advance(0.1, &[300.0], budget(), || Err(FissionError::Cancelled)),
        Err(DiffusionBurnupError::Fission(FissionError::Cancelled))
    ));
    assert_eq!(state.neutrons().neutron_counts(), initial_neutrons);
    assert_eq!(state.inventories()[0], initial_inventory);
    let (_, idle) = state.trial(0.0, &[300.0], budget(), || Ok(())).unwrap();
    assert_eq!(idle.burnup[0].fissions, 0.0);
    assert_eq!(idle.reaction_defect, 0.0);
}

#[test]
fn g0_missing_duplicate_and_wrong_kind_bindings_are_refused() {
    let state = fixture(1, true, false);
    assert!(
        DiffusionBurnupState::new(state.neutrons.clone(), state.inventories().to_vec(), vec![])
            .is_err()
    );
    let repeated = vec![
        ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: 0
        };
        2
    ];
    assert!(
        DiffusionBurnupState::new(
            state.neutrons.clone(),
            state.inventories().to_vec(),
            repeated
        )
        .is_err()
    );
    let capture = fixture(1, false, false);
    assert!(
        DiffusionBurnupState::new(
            capture.neutrons,
            state.inventories().to_vec(),
            vec![ReactionBinding {
                diffusion_reaction: 0,
                burnup_reaction: 0
            }]
        )
        .is_err()
    );
}

#[test]
fn g1_small_group_mismatches_survive_large_channel_totals() {
    let large = 2.0_f64.powi(60);
    let mut neutron_count = Accumulator::new(large);
    neutron_count.add(1.0).unwrap();
    // The endpoints round to the same scalar, but one event is unaccounted.
    assert_eq!(neutron_count.value(), large);
    assert_eq!(reaction_count_defect(&neutron_count, large).unwrap(), 1.0);
    let mut defects = Accumulator::new(large);
    defects.add(1.0).unwrap();
    assert_eq!(upper_value(&defects).unwrap(), large.next_up());
}
