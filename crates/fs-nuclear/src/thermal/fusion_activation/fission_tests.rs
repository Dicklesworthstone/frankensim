//! Prompt fission, actual radioactive product ages and material heating.
//! All isotope labels, event energies and rates are manufactured test data.

use super::*;
use crate::diffusion::LeakageBoundary;
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;
use crate::fusion_neutrons::NeutronReactionEnergy;

fn fission_stock() -> Arc<FusionNeutronActivationState> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-fission-material-grid".into(), vec![0.0, 4.0])
            .unwrap(),
    );
    let nuclide = |z, a, rate| BurnupNuclide {
        id: NuclideId::new(z, a, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    };
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-fission-and-delayed-neutron-products".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    nuclide(10, 21, 0.4),
                    nuclide(20, 41, 0.0),
                    nuclide(10, 20, 0.7),
                    nuclide(11, 20, 0.0),
                    nuclide(11, 19, 0.0),
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-driver-neutron-decay".into(),
                        parent: 0,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 1.0,
                            radiation_j: 2.0,
                            escaped_j: 0.5,
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-finite-fuel-fission".into(),
                        parent: 1,
                        rate: BurnupRate::Fission {
                            cross_sections_m2: vec![2e-4],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 2.0,
                        }],
                        // Prompt Q is owned by the neutron energy model.
                        energy_per_event: EnergyChannels::default(),
                    },
                    BurnupReaction {
                        data_id: "manufactured-product-beta-decay".into(),
                        parent: 2,
                        rate: BurnupRate::Decay {
                            branch_fraction: 0.75,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 3,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            radiation_j: 1.0,
                            escaped_j: 0.0,
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-product-delayed-neutron-decay".into(),
                        parent: 2,
                        rate: BurnupRate::Decay {
                            branch_fraction: 0.25,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 4,
                            yield_per_event: 1.0,
                        }],
                        // Parent-average energy is shared by both branches;
                        // the parent-wide mean neutron yield is 0.25.
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            radiation_j: 1.0,
                            escaped_j: 0.0,
                        },
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
                data_id: "manufactured-prompt-fission-material-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![0.01],
                delayed_groups: vec![],
                reactions: vec![DiffusionReaction {
                    data_id: "manufactured-bound-fuel-fission".into(),
                    cell: 0,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 0.8,
                    kind: DiffusionReactionKind::Fission {
                        prompt_yields: vec![2.0],
                        delayed_yields: vec![],
                    },
                }],
                links: vec![],
                leakage: vec![LeakageBoundary {
                    cell: 0,
                    conductance_m2: vec![0.003],
                }],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new_with_reaction_energies(
            "manufactured-fission-prompt-Q-and-heat".into(),
            Arc::clone(&diffusion),
            vec![2.0],
            vec![0.75],
            vec![NeutronReactionEnergy {
                nuclear_q_j: 10.0,
                threshold_j: 0.0,
            }],
        )
        .unwrap(),
    );
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "manufactured-driver-and-product-neutron-emission".into(),
            Arc::clone(&nuclear),
            vec![vec![1.0], vec![0.0], vec![0.25], vec![0.0], vec![0.0]],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0], vec![]).unwrap(),
        vec![
            BurnupInventory::new(nuclear, vec![20.0, 40.0, 0.0, 0.0, 0.0]).unwrap(),
        ],
        vec![ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: 1,
        }],
    )
    .unwrap()
    .with_decay_neutrons(vec![emission])
    .unwrap();
    Arc::new(FusionNeutronActivationState::new(model, spatial).unwrap())
}

#[test]
fn g1_g3_g4_finite_fission_and_delayed_product_neutrons_heat_original_material() {
    let initial = thermal(
        vec![cell(20, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let state = fission_stock();
    let sources = FusionNeutronActivationCellSource::bind_without_external_source(
        Arc::clone(&state),
        Time::new(0.0),
        vec![(CellId(20), MaterialId(1))],
    )
    .unwrap();
    let heating =
        FusionNeutronActivationHeating::new(&initial, sources, activation_budget()).unwrap();
    let (middle, checkpoint, _) = heating
        .advance_to(&initial, 0.2, controls(0.025), || Ok(()))
        .unwrap();
    let mut polls = 0;
    let (warm, accepted, _) = checkpoint
        .advance_to(&middle, 0.4, controls(0.025), || {
            polls += 1;
            Ok(())
        })
        .unwrap();

    for (time, material, heating) in [(0.2_f64, &middle, &checkpoint), (0.4, &warm, &accepted)] {
        let shared = heating.source(CellId(20)).unwrap().state();
        let inventory = &shared.inventories()[0];
        let p = inventory.populations();
        let events = inventory.cumulative_reaction_counts();
        let driver_decays = events[0];
        let fissions = events[1];
        let product_decays = events[2] + events[3];
        close(driver_decays, 20.0 * -(-0.4 * time).exp_m1(), 1e-8);
        close(p[0], 20.0 - driver_decays, 1e-8);
        close(p[1], 40.0 - fissions, 1e-8);
        close(p[2], driver_decays + 2.0 * fissions - product_decays, 1e-8);
        close(p[3], 0.75 * product_decays, 1e-8);
        close(p[4], 0.25 * product_decays, 1e-8);
        assert!(fissions > 0.0 && product_decays > 0.0);
        assert!(Arc::ptr_eq(inventory.model(), state.inventories()[0].model()));

        // E_in + Q - outgoing = 2 + 10 - 2*2 = 8 J per fission.
        // Its prompt heat/radiation is 6/2 J. Product neutron kinetic
        // energy transfers 0.5 J per product decay from its radiation.
        let heat = driver_decays + 6.0 * fissions + 2.0 * product_decays;
        let release = shared.cell_released(0).unwrap();
        close(release.deposited_j, heat, 1e-7);
        close(release.radiation_j, 2.0 * fissions + 0.5 * product_decays, 1e-7);
        close(material.temperature_k(CellId(20)).unwrap(), 300.0 + heat / 8.0, 1e-7);
        close(material.energy_j() - initial.energy_j(), heat, 1e-7);
        assert_eq!(
            material.cells()[&CellId(20)].mass_kg,
            initial.cells()[&CellId(20)].mass_kg
        );
        let leaked = (release.escaped_j - 0.5 * driver_decays) / 2.0;
        assert!(leaked > 0.0);
        close(
            shared.neutrons().neutron_counts()[0] + leaked,
            driver_decays + 0.25 * product_decays + fissions,
            1e-7,
        );
        close(
            shared.kinetic_energy_j().unwrap() + release.total_j(),
            3.5 * driver_decays + 10.0 * fissions + 3.0 * product_decays,
            1e-7,
        );
        close(heating.cumulative_released().unwrap().deposited_j, heat, 1e-7);
        let last = heating.last_step(CellId(20)).unwrap().activation();
        assert_eq!(last.external_injected_neutrons, 0.0);
        let counts = &last.spatial.burnup[0].reaction_counts;
        close(last.decay_neutrons, counts[0] + counts[3], 1e-8);
        close(last.spatial.neutrons.fissions, counts[1], 1e-8);
        assert!(counts[3] > 0.0);
        assert_eq!(last.kinetic.diffusion, last.spatial.neutrons);
    }

    // A late refusal cannot consume finite fuel or publish partial material heat.
    let before = Arc::clone(checkpoint.source(CellId(20)).unwrap().state());
    assert!(polls > 2);
    let mut remaining = polls - 2;
    assert!(
        checkpoint
            .advance_to(&middle, 0.4, controls(0.025), || {
                if remaining == 0 {
                    Err(DepositionError::Cancelled)
                } else {
                    remaining -= 1;
                    Ok(())
                }
            })
            .is_err()
    );
    assert!(Arc::ptr_eq(
        checkpoint.source(CellId(20)).unwrap().state(),
        &before
    ));
    assert_eq!(middle.time_s(), 0.2);
    assert_eq!(state.inventories()[0].populations(), &[20.0, 40.0, 0.0, 0.0, 0.0]);
    let (retry, _, _) = checkpoint
        .advance_to(&middle, 0.4, controls(0.025), || Ok(()))
        .unwrap();
    close(retry.energy_j(), warm.energy_j(), 1e-8);
}
