//! Closed-form transported radioactive births and complete history atomicity.
use std::sync::Arc;

use super::super::tests::close;
use super::BurnupDecayNeutronEmission;
use crate::burnup::transport::{
    BurnupTransportLimits, BurnupTransportState, NuclideOutflow, NuclideTransfer, NuclideTransport,
};
use crate::burnup::{
    BurnupData, BurnupInventory, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate,
    BurnupReaction,
};
use crate::diffusion::burnup::{
    DiffusionBurnupBudget, DiffusionBurnupError, DiffusionBurnupSources, DiffusionBurnupState,
    ReactionBinding,
};
use crate::diffusion::{
    DiffusionData, DiffusionModel, DiffusionReaction, DiffusionReactionKind, DiffusionState,
    LeakageBoundary,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use crate::{EnergyChannels, FissionError};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn budget() -> DiffusionBurnupBudget {
    let mut value = super::super::tests::budget();
    value.max_iterations = 128;
    value.reaction_tolerance = 1e-8;
    value
}

fn fixture(
    rates: &[f64],
    populations: Vec<Vec<f64>>,
    yields: &[f64],
    movement: NuclideTransport,
    neutron_loss: f64,
) -> DiffusionBurnupState {
    fixture_with_capture(rates, populations, yields, movement, neutron_loss, None)
}

fn fixture_with_capture(
    rates: &[f64],
    populations: Vec<Vec<f64>>,
    yields: &[f64],
    movement: NuclideTransport,
    neutron_loss: f64,
    capture: Option<f64>,
) -> DiffusionBurnupState {
    let cells = populations.len();
    let groups = Arc::new(
        NeutronGroups::new("manufactured-moving-emitter-grid".into(), vec![0.0, 4.0]).unwrap(),
    );
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-moving-emitter-chain".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: rates
                    .iter()
                    .enumerate()
                    .map(|(index, &rate)| BurnupNuclide {
                        id: NuclideId::new(10, 25 + index as u16, 0).unwrap(),
                        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
                    })
                    .collect(),
                reactions: {
                    let mut reactions: Vec<_> = rates
                        .iter()
                        .enumerate()
                        .filter(|(_, rate)| **rate > 0.0)
                        .map(|(parent, _)| BurnupReaction {
                            data_id: format!("moving-decay-{parent}"),
                            parent,
                            rate: BurnupRate::Decay {
                                branch_fraction: 1.0,
                            },
                            products: vec![BurnupProduct {
                                nuclide: parent + 1,
                                yield_per_event: 1.0,
                            }],
                            energy_per_event: EnergyChannels {
                                deposited_j: 2.0,
                                radiation_j: 3.0,
                                escaped_j: 1.0,
                            },
                        })
                        .collect();
                    if let Some(sigma) = capture {
                        reactions.push(BurnupReaction {
                            data_id: "moving-parent-capture".into(),
                            parent: 0,
                            rate: BurnupRate::Capture {
                                cross_sections_m2: vec![sigma],
                            },
                            products: vec![BurnupProduct {
                                nuclide: rates.len() - 1,
                                yield_per_event: 1.0,
                            }],
                            energy_per_event: EnergyChannels {
                                deposited_j: 1.0,
                                radiation_j: 0.0,
                                escaped_j: 0.0,
                            },
                        });
                    }
                    reactions
                },
            },
            budget().burnup.limits,
        )
        .unwrap(),
    );
    // Unequal volumes must not rescale extensive inventories or authored
    // donor-proportional movement. Leakage remains neutron_loss per neutron.
    let volumes: Vec<_> = (0..cells).map(|cell| 1.0 + cell as f64).collect();
    let initial_macroscopic: Vec<_> = capture.map_or_else(Vec::new, |sigma| {
        populations
            .iter()
            .zip(&volumes)
            .map(|(row, volume)| sigma * row[0] / volume)
            .collect()
    });
    let limits = BurnupTransportLimits {
        burnup: budget().burnup.limits,
        ..BurnupTransportLimits::default()
    };
    let transported = BurnupTransportState::new(
        populations
            .into_iter()
            .map(|row| BurnupInventory::new(Arc::clone(&model), row).unwrap())
            .collect(),
        volumes.clone(),
        movement,
        limits,
        || Ok(()),
    )
    .unwrap();
    let neutron_model = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-moving-emitter-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: volumes.clone(),
                delayed_groups: vec![],
                reactions: initial_macroscopic
                    .iter()
                    .enumerate()
                    .map(|(cell, &sigma)| DiffusionReaction {
                        data_id: format!("moving-parent-capture-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: sigma,
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![],
                leakage: volumes
                    .iter()
                    .enumerate()
                    .map(|(cell, volume)| LeakageBoundary {
                        cell,
                        conductance_m2: vec![neutron_loss * volume],
                    })
                    .collect(),
            },
            budget().diffusion.limits,
        )
        .unwrap(),
    );
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "manufactured-moving-emitter-yields".into(),
            model,
            yields.iter().map(|&yield_| vec![yield_]).collect(),
        )
        .unwrap(),
    );
    DiffusionBurnupState::new_with_transport(
        DiffusionState::new(neutron_model, vec![0.0; cells], vec![]).unwrap(),
        transported,
        if capture.is_some() {
            (0..cells)
                .map(|cell| ReactionBinding {
                    diffusion_reaction: cell,
                    burnup_reaction: rates.iter().filter(|&&rate| rate > 0.0).count(),
                })
                .collect()
        } else {
            vec![]
        },
        limits,
    )
    .unwrap()
    .with_decay_neutrons(vec![emission; cells])
    .unwrap()
}

fn transfer(rates: Vec<f64>) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![NuclideTransfer {
            from_compartment: 0,
            to_compartment: 1,
            rates_per_s: rates,
        }],
        outflows: vec![],
    }
}

#[test]
fn g1_moving_parent_emits_at_its_actual_location_and_birth_time() {
    let initial = fixture(
        &[1.0, 0.0],
        vec![vec![100.0, 0.0], vec![0.0, 0.0]],
        &[1.0, 0.0],
        transfer(vec![1.0, 0.0]),
        1.0,
    );
    let time = 1.25_f64;
    let e = (-time).exp();
    let e2 = (-2.0 * time).exp();
    let parents = [100.0 * e2, 100.0 * (e - e2)];
    let decays = [50.0 * (1.0 - e2), 100.0 * (1.0 - e) - 50.0 * (1.0 - e2)];
    let neutrons = [100.0 * (e - e2), 100.0 * ((time - 1.0) * e + e2)];
    let (next, report) = initial
        .trial(time, &[300.0; 2], budget(), || Ok(()))
        .unwrap();
    let movement = report.transport.as_ref().unwrap();
    for cell in 0..2 {
        close(
            next.inventories()[cell].populations()[0],
            parents[cell],
            2e-9,
        );
        close(
            next.inventories()[cell].populations()[1],
            decays[cell],
            2e-9,
        );
        close(next.neutrons().neutron_counts()[cell], neutrons[cell], 2e-9);
        close(
            report.decay_neutron_counts.as_ref().unwrap()[cell][0],
            decays[cell],
            2e-9,
        );
        close(report.burnup[cell].reaction_counts[0], decays[cell], 2e-9);
        close(
            report.neutrons.leakage_counts[cell],
            decays[cell] - neutrons[cell],
            2e-9,
        );
        close(
            movement.integrated_populations_s[cell][0],
            decays[cell],
            2e-9,
        );
        close(
            report.burnup[cell].released.deposited_j,
            2.0 * decays[cell],
            4e-9,
        );
        close(
            next.inventories()[cell].cumulative_released().deposited_j,
            2.0 * decays[cell],
            4e-9,
        );
    }
    close(movement.transfer_counts[0][0], decays[0], 2e-9);
    assert_eq!(movement.transfer_counts[0][1], 0.0);
    close(
        next.transport().unwrap().cumulative_transfer_counts()[0][0],
        decays[0],
        2e-9,
    );
    close(
        parents.iter().sum::<f64>() + decays.iter().sum::<f64>(),
        100.0,
        1e-10,
    );
    for residual in &movement.cohort_balance_residuals {
        close(*residual, 0.0, 2e-9);
    }
    assert_eq!(
        initial.transport().unwrap().cumulative_transfer_counts(),
        vec![vec![0.0, 0.0]]
    );
    assert_eq!(
        initial.inventories()[1].cumulative_released(),
        EnergyChannels::default()
    );
}

#[test]
fn g1_daughter_is_born_moved_and_emits_in_the_same_interval() {
    // P_A -> D_A; D moves A->B and decays to S with unit rates.
    // Only D emits a neutron, whose local leakage rate is also one.
    let initial = fixture(
        &[1.0, 1.0, 0.0],
        vec![vec![100.0, 0.0, 0.0], vec![0.0; 3]],
        &[0.0, 1.0, 0.0],
        transfer(vec![0.0, 1.0, 0.0]),
        1.0,
    );
    let time = 1.2_f64;
    let e = (-time).exp();
    let e2 = e * e;
    let parent_decays = 100.0 * (1.0 - e);
    let daughters = [100.0 * (e - e2), 100.0 * ((time - 1.0) * e + e2)];
    let local_decays = 50.0 * (1.0 - e) * (1.0 - e);
    let decays = [
        local_decays,
        100.0 * (1.0 - (1.0 + time) * e) - local_decays,
    ];
    let neutrons = [
        100.0 * ((time - 1.0) * e + e2),
        100.0 * ((0.5 * time * time - time + 1.0) * e - e2),
    ];
    let (next, report) = initial
        .trial(time, &[300.0; 2], budget(), || Ok(()))
        .unwrap();
    close(next.inventories()[0].populations()[0], 100.0 * e, 2e-9);
    assert_eq!(next.inventories()[1].populations()[0], 0.0);
    close(report.burnup[0].reaction_counts[0], parent_decays, 2e-9);
    assert_eq!(report.burnup[1].reaction_counts[0], 0.0);
    for cell in 0..2 {
        close(
            next.inventories()[cell].populations()[1],
            daughters[cell],
            2e-9,
        );
        close(
            next.inventories()[cell].populations()[2],
            decays[cell],
            2e-9,
        );
        close(report.burnup[cell].reaction_counts[1], decays[cell], 2e-9);
        close(
            report.decay_neutron_counts.as_ref().unwrap()[cell][0],
            decays[cell],
            2e-9,
        );
        close(next.neutrons().neutron_counts()[cell], neutrons[cell], 2e-9);
        close(
            report.neutrons.leakage_counts[cell],
            decays[cell] - neutrons[cell],
            2e-9,
        );
    }
    let movement = report.transport.as_ref().unwrap();
    close(movement.transfer_counts[0][1], local_decays, 2e-9);
    assert_eq!(movement.transfer_counts[0][0], 0.0);
    close(
        report.burnup[0].released.deposited_j,
        2.0 * (parent_decays + decays[0]),
        1e-8,
    );
    close(report.burnup[1].released.deposited_j, 2.0 * decays[1], 2e-9);
    assert!(decays[1] > 0.0 && neutrons[1] > 0.0);
}

#[test]
fn g1_outflow_removes_future_retained_emission_without_becoming_neutron_leakage() {
    let initial = fixture(
        &[1.0, 0.0],
        vec![vec![100.0, 0.0]],
        &[1.0, 0.0],
        NuclideTransport {
            transfers: vec![],
            outflows: vec![NuclideOutflow {
                compartment: 0,
                rates_per_s: vec![1.0, 0.0],
            }],
        },
        0.0,
    );
    let (middle, first) = initial.trial(1.0, &[300.0], budget(), || Ok(())).unwrap();
    let (next, second) = middle
        .clone()
        .trial(1.0, &[300.0], budget(), || Ok(()))
        .unwrap();
    let e2 = (-2.0_f64).exp();
    let e4 = (-4.0_f64).exp();
    let first_decays = 50.0 * (1.0 - e2);
    let second_decays = 50.0 * (e2 - e4);
    close(
        first.decay_neutron_counts.as_ref().unwrap()[0][0],
        first_decays,
        2e-9,
    );
    close(
        second.decay_neutron_counts.as_ref().unwrap()[0][0],
        second_decays,
        2e-9,
    );
    close(
        first.transport.as_ref().unwrap().outflow_counts[0][0],
        first_decays,
        2e-9,
    );
    close(
        second.transport.as_ref().unwrap().outflow_counts[0][0],
        second_decays,
        2e-9,
    );
    close(next.inventories()[0].populations()[0], 100.0 * e4, 2e-9);
    let retained_decays = first_decays + second_decays;
    close(
        next.inventories()[0].populations()[1],
        retained_decays,
        2e-9,
    );
    close(next.neutrons().neutron_counts()[0], retained_decays, 2e-9);
    close(
        next.transport().unwrap().cumulative_outflow_counts()[0][0],
        retained_decays,
        2e-9,
    );
    close(
        next.inventories()[0].cumulative_released().deposited_j,
        2.0 * retained_decays,
        4e-9,
    );
    assert_eq!(first.neutrons.leaked_neutrons, 0.0);
    assert_eq!(second.neutrons.leaked_neutrons, 0.0);
    assert!(next.neutrons().neutron_counts()[0] < 50.0);
    close(
        next.inventories()[0].populations().iter().sum::<f64>() + retained_decays,
        100.0,
        4e-9,
    );
}

#[test]
fn g3_continuous_feed_with_zero_flow_reduces_to_the_stationary_source_owner() {
    let initial = fixture(
        &[1.0, 0.0],
        vec![vec![0.0; 2]; 2],
        &[0.5, 0.0],
        NuclideTransport::default(),
        1.0,
    );
    let independent = DiffusionBurnupState::new(
        initial.neutrons().clone(),
        initial.inventories().to_vec(),
        vec![],
    )
    .unwrap()
    .with_decay_neutrons(initial.decay_neutrons().unwrap().to_vec())
    .unwrap();
    let feed = [4.0, 8.0];
    let external = [1.0, 0.0];
    let sources = DiffusionBurnupSources {
        neutron_rates_s: external.to_vec(),
        precursor_rates_s: vec![],
        nuclide_rates_s: feed.iter().map(|&rate| vec![rate, 0.0]).collect(),
    };
    let time = 1.5_f64;
    let (joint, report) = initial
        .trial_with_sources(time, &[300.0; 2], &sources, budget(), || Ok(()))
        .unwrap();
    let (stationary, comparison) = independent
        .trial_with_sources(time, &[300.0; 2], &sources, budget(), || Ok(()))
        .unwrap();
    let e = (-time).exp();
    for cell in 0..2 {
        let parent = feed[cell] * (1.0 - e);
        let decays = feed[cell] * (time - 1.0 + e);
        let neutrons = external[cell] * (1.0 - e) + 0.5 * feed[cell] * (1.0 - (1.0 + time) * e);
        close(joint.inventories()[cell].populations()[0], parent, 2e-9);
        close(joint.inventories()[cell].populations()[1], decays, 2e-9);
        close(joint.neutrons().neutron_counts()[cell], neutrons, 2e-9);
        close(
            report.decay_neutron_counts.as_ref().unwrap()[cell][0],
            0.5 * decays,
            2e-9,
        );
        close(
            report.neutrons.injected_neutron_counts[cell],
            external[cell] * time + 0.5 * decays,
            2e-9,
        );
        close(
            report.neutrons.leakage_counts[cell],
            external[cell] * time + 0.5 * decays - neutrons,
            2e-9,
        );
        close(
            report.produced_nuclei_by_cell[cell],
            feed[cell] * time,
            1e-10,
        );
        close(
            report.transport.as_ref().unwrap().integrated_populations_s[cell][0],
            decays,
            2e-9,
        );
        close(
            joint.neutrons().neutron_counts()[cell],
            stationary.neutrons().neutron_counts()[cell],
            2e-9,
        );
        for (&a, &b) in joint.inventories()[cell]
            .populations()
            .iter()
            .zip(stationary.inventories()[cell].populations())
        {
            close(a, b, 2e-9);
        }
        close(
            report.burnup[cell].released.total_j(),
            comparison.burnup[cell].released.total_j(),
            1e-8,
        );
    }
    assert!(
        report
            .transport
            .as_ref()
            .unwrap()
            .transfer_counts
            .is_empty()
    );
    assert!(report.transport.as_ref().unwrap().outflow_counts.is_empty());
    assert!(comparison.transport.is_none());
}

#[test]
fn g3_moving_emitters_capture_their_neutrons_and_converge_to_an_independent_ode() {
    // [P_A,P_B,n_A,n_B,D_A,D_B,C_A,C_B,escaped neutrons]. The cell
    // volumes are 1 and 2, so equal microscopic sigma gives k_A=.02,k_B=.01.
    fn reference(steps: usize) -> [f64; 9] {
        fn rhs(x: [f64; 9]) -> [f64; 9] {
            let capture_a = 0.02 * x[0] * x[2];
            let capture_b = 0.01 * x[1] * x[3];
            let movement = 0.4 * x[0];
            [
                -x[0] - capture_a - movement,
                movement - x[1] - capture_b,
                x[0] - capture_a - 0.3 * x[2],
                x[1] - capture_b - 0.3 * x[3],
                x[0],
                x[1],
                capture_a,
                capture_b,
                0.3 * (x[2] + x[3]),
            ]
        }
        let h = 2.0 / steps as f64;
        let mut x = [100.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for _ in 0..steps {
            let a = rhs(x);
            let b = rhs(std::array::from_fn(|j| x[j] + 0.5 * h * a[j]));
            let c = rhs(std::array::from_fn(|j| x[j] + 0.5 * h * b[j]));
            let d = rhs(std::array::from_fn(|j| x[j] + h * c[j]));
            x = std::array::from_fn(|j| x[j] + h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0);
        }
        x
    }
    let expected = reference(2048);
    for (a, b) in reference(1024).into_iter().zip(expected) {
        close(a, b, 5e-10);
    }
    let initial = fixture_with_capture(
        &[1.0, 0.0, 0.0],
        vec![vec![100.0, 0.0, 0.0], vec![20.0, 0.0, 0.0]],
        &[1.0, 0.0, 0.0],
        transfer(vec![0.4, 0.0, 0.0]),
        0.3,
        Some(0.02),
    );
    // Moving an unchanged nuclide must not silently change its decay spectrum.
    let mut mismatched = initial.decay_neutrons().unwrap().to_vec();
    mismatched[1] = Arc::new(
        BurnupDecayNeutronEmission::new(
            "different-receiver-decay-spectrum".into(),
            Arc::clone(initial.inventories()[1].model()),
            vec![vec![0.5], vec![0.0], vec![0.0]],
        )
        .unwrap(),
    );
    assert!(initial.clone().with_decay_neutrons(mismatched).is_err());

    let mut errors = Vec::new();
    for steps in [32, 64] {
        let mut state = initial.clone();
        let mut escaped = 0.0;
        let mut births = [0.0; 2];
        let mut iterated = false;
        for _ in 0..steps {
            let report = state
                .advance(2.0 / steps as f64, &[300.0; 2], budget(), || Ok(()))
                .unwrap();
            iterated |= report.iterations > 1;
            assert!(report.reaction_defect <= budget().reaction_tolerance);
            close(
                report.neutrons.captures,
                report.burnup.iter().map(|cell| cell.captures).sum(),
                2e-8,
            );
            escaped += report.neutrons.leaked_neutrons;
            for cell in 0..2 {
                births[cell] += report.decay_neutron_counts.as_ref().unwrap()[cell][0];
            }
        }
        assert!(iterated);
        let a = state.inventories()[0].populations();
        let b = state.inventories()[1].populations();
        let neutrons = state.neutrons().neutron_counts();
        let actual = [
            a[0],
            b[0],
            neutrons[0],
            neutrons[1],
            a[1],
            b[1],
            a[2],
            b[2],
            escaped,
        ];
        let moved = state.transport().unwrap().cumulative_transfer_counts()[0][0];
        close(moved, 0.4 * a[1], 2e-7);
        close(a.iter().sum::<f64>() + moved, 100.0, 2e-7);
        close(b.iter().sum::<f64>() - moved, 20.0, 2e-7);
        close(
            neutrons.iter().sum::<f64>() + a[2] + b[2] + escaped,
            a[1] + b[1],
            8e-6,
        );
        for (cell, populations) in [a, b].into_iter().enumerate() {
            close(births[cell], populations[1], 2e-7);
            close(
                state.inventories()[cell].cumulative_released().deposited_j,
                2.0 * populations[1] + populations[2],
                2e-7,
            );
            assert!(populations[2] > 1.0);
        }
        errors.push(
            actual
                .into_iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    assert!(errors[1] < errors[0] / 2.0, "{errors:?}");
    assert!(errors[1] < 0.01, "{errors:?}");
}

#[test]
fn g4_late_cancellation_and_joint_limits_preserve_nonzero_complete_histories() {
    let mut movement = transfer(vec![1.0, 0.0]);
    movement.outflows.push(NuclideOutflow {
        compartment: 1,
        rates_per_s: vec![0.2, 0.0],
    });
    let mut state = fixture(
        &[1.0, 0.0],
        vec![vec![100.0, 0.0], vec![0.0, 0.0]],
        &[1.0, 0.0],
        movement,
        1.0,
    );
    state
        .advance(0.3, &[300.0; 2], budget(), || Ok(()))
        .unwrap();
    let original = state.clone();
    assert!(original.transport().unwrap().cumulative_transfer_counts()[0][0] > 0.0);
    assert!(original.transport().unwrap().cumulative_outflow_counts()[0][0] > 0.0);
    assert!(original.inventories()[1].cumulative_released().deposited_j > 0.0);
    assert!(original.neutrons().neutron_counts()[1] > 0.0);
    let mut checkpoints = 0;
    state
        .trial(0.7, &[300.0; 2], budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    let cancelled = state.advance(0.7, &[300.0; 2], budget(), || {
        seen += 1;
        if seen == checkpoints {
            Err(FissionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(DiffusionBurnupError::Fission(FissionError::Cancelled))
    ));
    let mut limited = budget();
    // Four nuclide nodes fit their own transport owner; the coupled neutron
    // field adds two more, so the complete six-node solve must be refused.
    limited.burnup.limits.propagation.maximum_nodes = 5;
    limited.diffusion.limits.propagation.maximum_nodes = 5;
    assert!(state.advance(0.7, &[300.0; 2], limited, || Ok(())).is_err());
    assert_eq!(state.transport(), original.transport());
    assert_eq!(state.inventories(), original.inventories());
    assert_eq!(
        state.neutrons().neutron_counts(),
        original.neutrons().neutron_counts()
    );
}
