//! Independent source-age and finite-emitter depletion regressions.
use std::sync::Arc;

use super::super::tests::close;
use super::BurnupDecayNeutronEmission;
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
    let mut budget = super::super::tests::budget();
    budget.max_iterations = 128;
    budget.reaction_tolerance = 1e-8;
    budget
}

// Every unstable node feeds its next daughter. An optional (target, sigma)
// capture feeds the final, distinct stable capture product.
fn unbound(
    rates: &[f64],
    populations: &[f64],
    capture: Option<(usize, f64)>,
    loss: f64,
) -> DiffusionBurnupState {
    let groups =
        Arc::new(NeutronGroups::new("manufactured-emitter-grid".into(), vec![0.0, 4.0]).unwrap());
    let mut reactions: Vec<_> = rates
        .iter()
        .enumerate()
        .filter(|(_, rate)| **rate > 0.0)
        .map(|(parent, _)| BurnupReaction {
            data_id: format!("decay-{parent}"),
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
    let mut bindings = Vec::new();
    let mut neutron_reactions = Vec::new();
    if let Some((parent, sigma)) = capture {
        bindings.push(ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: reactions.len(),
        });
        reactions.push(BurnupReaction {
            data_id: "finite-emitter-capture".into(),
            parent,
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
        neutron_reactions.push(DiffusionReaction {
            data_id: "finite-emitter-capture".into(),
            cell: 0,
            incident_group: 0,
            macroscopic_cross_section_m_inv: sigma * populations[parent],
            kind: DiffusionReactionKind::Capture,
        });
    }
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-finite-radioactive-emitter".into(),
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
                reactions,
            },
            budget().burnup.limits,
        )
        .unwrap(),
    );
    let spatial = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-emitter-transport".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                delayed_groups: vec![],
                reactions: neutron_reactions,
                links: vec![],
                leakage: vec![LeakageBoundary {
                    cell: 0,
                    conductance_m2: vec![loss],
                }],
            },
            budget().diffusion.limits,
        )
        .unwrap(),
    );
    DiffusionBurnupState::new(
        DiffusionState::new(spatial, vec![0.0], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, populations.to_vec()).unwrap()],
        bindings,
    )
    .unwrap()
}

fn bind(state: DiffusionBurnupState, yields: &[f64]) -> DiffusionBurnupState {
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "manufactured-parent-neutron-yields".into(),
            Arc::clone(state.inventories()[0].model()),
            yields.iter().map(|&yield_| vec![yield_]).collect(),
        )
        .unwrap(),
    );
    state.with_decay_neutrons(vec![emission]).unwrap()
}

fn finite_emitter() -> DiffusionBurnupState {
    bind(
        unbound(&[1.0, 0.0, 0.0], &[100.0, 0.0, 0.0], Some((0, 0.02)), 0.3),
        &[1.0, 0.0, 0.0],
    )
}

#[test]
fn g1_radioactive_births_keep_their_ages_without_neutron_capture() {
    let initial = bind(unbound(&[1.0, 0.0], &[100.0, 0.0], None, 2.0), &[0.5, 0.0]);
    let (next, report) = initial.trial(1.0, &[300.0], budget(), || Ok(())).unwrap();
    let decays = 100.0 * -(-1.0_f64).exp_m1();
    let born = 0.5 * decays;
    let neutrons = 50.0 * ((-1.0_f64).exp() - (-2.0_f64).exp());
    close(next.neutrons().neutron_counts()[0], neutrons, 1e-9);
    close(next.inventories()[0].populations()[0], 100.0 - decays, 1e-9);
    close(next.inventories()[0].populations()[1], decays, 1e-9);
    close(report.burnup[0].reaction_counts[0], decays, 1e-9);
    close(
        report.decay_neutron_counts.as_ref().unwrap()[0][0],
        born,
        1e-9,
    );
    close(report.neutrons.injected_neutrons, born, 1e-9);
    close(report.neutrons.leaked_neutrons, born - neutrons, 1e-9);
    close(
        report.neutrons.integrated_fluence_m2[0],
        (born - neutrons) / 2.0,
        1e-9,
    );
    close(report.burnup[0].released.deposited_j, 2.0 * decays, 1e-8);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0]);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0]);
}

// Independent nonlinear oracle: [parent, neutrons, decay daughter,
// capture product, escaped neutrons]. No production propagator is called.
fn reference(steps: usize) -> [f64; 5] {
    fn rhs(x: [f64; 5]) -> [f64; 5] {
        let decay = x[0];
        let capture = 0.02 * x[1] * x[0];
        let loss = 0.3 * x[1];
        [
            -decay - capture,
            decay - capture - loss,
            decay,
            capture,
            loss,
        ]
    }
    let h = 2.0 / steps as f64;
    let mut x = [100.0, 0.0, 0.0, 0.0, 0.0];
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
fn g3_physical_refinement_resolves_capture_of_the_same_radioactive_stock() {
    let expected = reference(2048);
    for (coarse, fine) in reference(1024).into_iter().zip(expected) {
        close(coarse, fine, 2e-10);
    }
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut state = finite_emitter();
        let mut escaped = 0.0;
        let mut born = 0.0;
        for _ in 0..steps {
            let report = state
                .advance(2.0 / steps as f64, &[300.0], budget(), || Ok(()))
                .unwrap();
            assert!(report.reaction_defect <= budget().reaction_tolerance);
            close(report.neutrons.captures, report.burnup[0].captures, 2e-8);
            escaped += report.neutrons.leaked_neutrons;
            born += report.decay_neutron_counts.as_ref().unwrap()[0][0];
        }
        let populations = state.inventories()[0].populations();
        let events = state.inventories()[0].cumulative_reaction_counts();
        let actual = [
            populations[0],
            state.neutrons().neutron_counts()[0],
            populations[1],
            populations[2],
            escaped,
        ];
        close(populations.iter().sum(), 100.0, 1e-7);
        close(events[0], populations[1], 1e-7);
        close(events[1], populations[2], 1e-7);
        close(born, populations[1], 8e-6);
        close(actual[1] + actual[3] + actual[4], actual[2], 8e-6);
        // Integrating P'/P=-1-0.02*n supplies an independent exposure identity.
        close(
            (actual[0] / 100.0).ln() + 2.0 + (0.02 / 0.3) * escaped,
            0.0,
            8e-6,
        );
        assert!(actual.iter().all(|value| *value >= 0.0));
        errors.push(
            actual
                .into_iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    assert!(errors[1] < errors[0] / 2.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 2.0, "{errors:?}");
    assert!(errors[2] < 0.005, "{errors:?}");
    assert!(expected[0] < 0.5 * 100.0 * (-2.0_f64).exp());
}

#[test]
fn g1_fed_parent_and_its_daughter_emit_before_the_interval_ends() {
    let initial = bind(
        unbound(&[1.0, 1.0, 0.0], &[0.0, 0.0, 0.0], None, 1.0),
        &[0.0, 1.0, 0.0],
    );
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![4.0, 0.0, 0.0]],
    };
    let duration = 2.0_f64;
    let (next, report) = initial
        .trial_with_sources(duration, &[300.0], &sources, budget(), || Ok(()))
        .unwrap();
    let tail = (-duration).exp();
    let parent = 4.0 * (1.0 - tail);
    let daughter = 4.0 * (1.0 - (1.0 + duration) * tail);
    let stable = 4.0 * (duration - 2.0 + (2.0 + duration) * tail);
    let neutrons = 4.0 * (1.0 - (1.0 + duration + 0.5 * duration * duration) * tail);
    for (actual, expected) in next.inventories()[0]
        .populations()
        .iter()
        .zip([parent, daughter, stable])
    {
        close(*actual, expected, 1e-9);
    }
    close(next.neutrons().neutron_counts()[0], neutrons, 1e-9);
    close(
        report.decay_neutron_counts.as_ref().unwrap()[0][0],
        stable,
        1e-9,
    );
    close(
        report.burnup[0].reaction_counts[0],
        4.0 * duration - parent,
        1e-9,
    );
    close(report.burnup[0].reaction_counts[1], stable, 1e-9);
    close(report.neutrons.leaked_neutrons, stable - neutrons, 1e-9);
    assert_eq!(report.produced_nuclei, vec![vec![8.0, 0.0, 0.0]]);
    assert_eq!(
        next.inventories()[0].cumulative_produced_nuclei(),
        vec![8.0, 0.0, 0.0]
    );
    assert_eq!(initial.inventories()[0].populations(), &[0.0, 0.0, 0.0]);
}

#[test]
fn g3_newly_created_neutrons_and_capture_targets_do_not_false_converge_at_zero() {
    // P -> D+n; D+n -> C. Both D and n start empty, so zero initial
    // coefficients alone would falsely agree on zero capture events.
    let initial = bind(
        unbound(&[1.0, 0.0, 0.0], &[100.0, 0.0, 0.0], Some((1, 0.02)), 0.0),
        &[1.0, 0.0, 0.0],
    );
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0]);
    assert_eq!(
        initial.neutrons().model().reactions()[0].macroscopic_cross_section_m_inv,
        0.0
    );
    let (next, report) = initial.trial(1.0, &[300.0], budget(), || Ok(())).unwrap();
    let populations = next.inventories()[0].populations();
    let neutrons = next.neutrons().neutron_counts()[0];
    let captures = report.burnup[0].captures;
    let decays = 100.0 * -(-1.0_f64).exp_m1();
    assert!(report.iterations > 1);
    assert!(captures > 1.0, "startup capture was lost: {captures}");
    close(populations[0], 100.0 - decays, 1e-8);
    close(populations[1], neutrons, 2e-8);
    close(populations[2], captures, 1e-8);
    close(report.neutrons.captures, captures, 2e-8);
    close(populations.iter().sum(), 100.0, 1e-8);
    close(neutrons + captures, decays, 2e-8);
    close(
        report.decay_neutron_counts.as_ref().unwrap()[0][0],
        decays,
        1e-8,
    );
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
}

#[test]
fn g4_restart_and_late_refusal_preserve_all_shared_histories() {
    let mut initial = finite_emitter();
    let original = initial.clone();
    let (middle, _) = initial.trial(0.25, &[300.0], budget(), || Ok(())).unwrap();
    let (continued, a) = middle.trial(0.25, &[300.0], budget(), || Ok(())).unwrap();
    let (restarted, b) = middle
        .clone()
        .trial(0.25, &[300.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(continued.inventories(), restarted.inventories());
    assert_eq!(
        continued.neutrons().neutron_counts(),
        restarted.neutrons().neutron_counts()
    );
    assert_eq!(a.decay_neutron_counts, b.decay_neutron_counts);

    let mut checkpoints = 0;
    initial
        .trial(0.25, &[300.0], budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        initial.advance(0.25, &[300.0], budget(), || {
            seen += 1;
            if seen == checkpoints {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DiffusionBurnupError::Fission(FissionError::Cancelled))
    ));
    let mut limited = budget();
    limited.max_iterations = 1;
    assert!(matches!(
        initial.advance(0.25, &[300.0], limited, || Ok(())),
        Err(DiffusionBurnupError::Nonconvergence { .. })
    ));
    let mut limited = budget();
    limited.diffusion.limits.propagation.maximum_nodes = 3;
    limited.burnup.limits.propagation.maximum_nodes = 3;
    assert!(initial.advance(0.25, &[300.0], limited, || Ok(())).is_err());
    assert_eq!(initial.inventories(), original.inventories());
    assert_eq!(
        initial.neutrons().neutron_counts(),
        original.neutrons().neutron_counts()
    );
}

#[test]
fn g0_emission_requires_the_original_cell_models_and_valid_unstable_rows() {
    let initial = unbound(&[1.0, 0.0], &[100.0, 0.0], None, 0.3);
    let model = Arc::clone(initial.inventories()[0].model());
    for yields in [
        vec![vec![1.0]],
        vec![vec![1.0], vec![0.1]],
        vec![vec![-1.0], vec![0.0]],
    ] {
        assert!(
            BurnupDecayNeutronEmission::new("invalid-yield".into(), Arc::clone(&model), yields)
                .is_err()
        );
    }
    assert!(initial.clone().with_decay_neutrons(vec![]).is_err());
    let foreign = unbound(&[1.0, 0.0], &[100.0, 0.0], None, 0.3);
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "foreign-original-model".into(),
            Arc::clone(foreign.inventories()[0].model()),
            vec![vec![1.0], vec![0.0]],
        )
        .unwrap(),
    );
    assert!(initial.clone().with_decay_neutrons(vec![emission]).is_err());
    let (_, unbound_report) = initial.trial(1.0, &[300.0], budget(), || Ok(())).unwrap();
    assert!(unbound_report.decay_neutron_counts.is_none());
    assert_eq!(unbound_report.neutrons.injected_neutrons, 0.0);
}
