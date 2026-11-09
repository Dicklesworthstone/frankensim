//! G1 analytic radioactive production/absorption and G3 history/conservation.
use super::*;
use crate::EnergyChannels;
use crate::decay::{DecayBranch, DecayNuclide};
use crate::diffusion::{
    DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction, DiffusionReactionKind,
    LeakageBoundary,
};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

pub(super) fn budget() -> DecayNeutronBudget {
    DecayNeutronBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        decay: DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    }
}

fn chain(rates: &[f64]) -> Arc<DecayModel> {
    let nodes = rates
        .iter()
        .enumerate()
        .map(|(i, &rate)| DecayNuclide {
            id: format!("manufactured-nuclide-{i}"),
            law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
        })
        .collect();
    let branches = rates
        .iter()
        .enumerate()
        .filter(|(_, rate)| **rate > 0.0)
        .map(|(i, _)| DecayBranch {
            parent: i,
            daughter: i + 1,
            probability: 1.0,
            energy_per_decay: EnergyChannels {
                deposited_j: 2.0,
                radiation_j: 3.0,
                escaped_j: 1.0,
            },
        })
        .collect();
    Arc::new(
        DecayModel::new(
            "manufactured-chain".into(),
            nodes,
            branches,
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn diffusion(capture: f64, leakage: f64) -> DiffusionState {
    let groups = Arc::new(NeutronGroups::new("manufactured-grid".into(), vec![0.0, 4.0]).unwrap());
    let data = DiffusionData {
        data_id: "manufactured-absorber".into(),
        groups,
        speeds_m_s: vec![1.0],
        cell_volumes_m3: vec![1.0],
        delayed_groups: vec![],
        reactions: vec![DiffusionReaction {
            data_id: "capture".into(),
            cell: 0,
            incident_group: 0,
            macroscopic_cross_section_m_inv: capture,
            kind: DiffusionReactionKind::Capture,
        }],
        links: vec![],
        leakage: vec![LeakageBoundary {
            cell: 0,
            conductance_m2: vec![leakage],
        }],
    };
    DiffusionState::new(
        Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
        vec![0.0],
        vec![],
    )
    .unwrap()
}

pub(super) fn state(
    rates: &[f64],
    populations: &[f64],
    yields: &[f64],
    capture: f64,
    leakage: f64,
) -> DecayNeutronState {
    let neutrons = diffusion(capture, leakage);
    let model = chain(rates);
    let emission = Arc::new(
        DecayNeutronEmission::new(
            "manufactured-neutron-yield".into(),
            Arc::clone(&model),
            Arc::clone(neutrons.model().groups()),
            yields.iter().map(|&yield_| vec![yield_]).collect(),
        )
        .unwrap(),
    );
    let inventory = DecayInventory::new(model, populations.to_vec()).unwrap();
    DecayNeutronState::new(
        neutrons,
        vec![DecayNeutronSource {
            cell: 0,
            inventory,
            emission,
        }],
    )
    .unwrap()
}

pub(super) fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_decay_birth_times_determine_same_interval_capture_and_leakage() {
    let initial = state(&[1.0, 0.0], &[100.0, 0.0], &[0.5, 0.0], 1.2, 0.8);
    let (next, step) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let decays = 100.0 * -(-1.0_f64).exp_m1();
    let born = 0.5 * decays;
    let neutrons = 50.0 * ((-1.0_f64).exp() - (-2.0_f64).exp());
    close(next.neutrons().neutron_counts()[0], neutrons);
    close(step.emitted_neutrons[0][0], born);
    close(step.neutrons.injected_neutrons, born);
    close(step.neutrons.captures, 0.6 * (born - neutrons));
    close(step.neutrons.leaked_neutrons, 0.4 * (born - neutrons));
    close(
        step.neutrons.integrated_fluence_m2[0],
        (born - neutrons) / 2.0,
    );
    close(step.decay[0].decays[0], decays);
    close(step.decay[0].released.deposited_j, 2.0 * decays);
    close(
        next.sources()[0].inventory.populations()[0],
        100.0 * (-1.0_f64).exp(),
    );
    let (uniform, _) = initial
        .neutrons()
        .trial_with_sources(1.0, &[born], &[], budget().diffusion, || Ok(()))
        .unwrap();
    assert!(
        (uniform.neutron_counts()[0] - neutrons).abs() > 1.0,
        "an interval-average source must not masquerade as the true birth history"
    );
    close(step.neutrons.neutron_balance_residual, 0.0);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0]);
}

#[test]
fn g1_daughters_emit_and_neutrons_are_absorbed_in_their_birth_interval() {
    // P -> D -> S, lambda_P = lambda_D = capture = 1. D emits one neutron.
    // D=100 t exp(-t), n=50 t^2 exp(-t), including equal-rate limits.
    let initial = state(
        &[1.0, 1.0, 0.0],
        &[100.0, 0.0, 0.0],
        &[0.0, 1.0, 0.0],
        1.0,
        0.0,
    );
    let (next, step) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let tail = (-2.0_f64).exp();
    let born = 100.0 * (1.0 - 3.0 * tail);
    close(next.sources()[0].inventory.populations()[1], 200.0 * tail);
    close(next.neutrons().neutron_counts()[0], 200.0 * tail);
    close(step.emitted_neutrons[0][0], born);
    close(step.neutrons.captures, born - 200.0 * tail);
    close(step.decay[0].decays[1], born);
}

#[test]
fn g3_partition_and_restart_preserve_source_events_residence_and_energy() {
    let initial = state(
        &[0.7, 0.2, 0.0],
        &[40.0, 6.0, 9.0],
        &[0.3, 1.4, 0.0],
        0.6,
        0.1,
    );
    let (whole, total) = initial.trial(3.0, budget(), || Ok(())).unwrap();
    let (middle, first) = initial.trial(1.5, budget(), || Ok(())).unwrap();
    let restart = middle.clone();
    let (split, second) = restart.trial(1.5, budget(), || Ok(())).unwrap();
    assert!(Arc::ptr_eq(
        middle.prepared.as_ref().unwrap(),
        split.prepared.as_ref().unwrap()
    ));
    close(
        split.neutrons().neutron_counts()[0],
        whole.neutrons().neutron_counts()[0],
    );
    for (a, b) in split.sources()[0]
        .inventory
        .populations()
        .iter()
        .zip(whole.sources()[0].inventory.populations())
    {
        close(*a, *b);
    }
    close(
        first.neutrons.captures + second.neutrons.captures,
        total.neutrons.captures,
    );
    close(
        first.neutrons.integrated_fluence_m2[0] + second.neutrons.integrated_fluence_m2[0],
        total.neutrons.integrated_fluence_m2[0],
    );
    close(
        split.sources()[0].inventory.cumulative_released().total_j(),
        whole.sources()[0].inventory.cumulative_released().total_j(),
    );
}

#[test]
fn g3_stable_sources_zero_duration_and_existing_neutrons_keep_their_meanings() {
    let initial = state(&[0.0], &[8.0], &[0.0], 0.3, 0.2);
    let (next, step) = initial.trial(0.0, budget(), || Ok(())).unwrap();
    assert_eq!(next.sources()[0].inventory, initial.sources()[0].inventory);
    assert_eq!(step.emitted_neutrons, vec![vec![0.0]]);
    let mut neutron_state = diffusion(0.3, 0.2);
    neutron_state =
        DiffusionState::new(Arc::clone(neutron_state.model()), vec![20.0], vec![]).unwrap();
    let expected = neutron_state
        .trial(2.0, budget().diffusion, || Ok(()))
        .unwrap();
    let source_free = DecayNeutronState::new(neutron_state, vec![]).unwrap();
    let actual = source_free.trial(2.0, budget(), || Ok(())).unwrap();
    assert_eq!(
        actual.0.neutrons().neutron_counts(),
        expected.0.neutron_counts()
    );
    assert_eq!(actual.1.neutrons, expected.1);
}

#[test]
fn g4_late_cancellation_and_joint_work_limits_preserve_all_histories() {
    let mut initial = state(&[1.0, 0.0], &[100.0, 0.0], &[0.5, 0.0], 1.2, 0.8);
    let original = initial.clone();
    let mut checkpoints = 0;
    initial
        .trial(1.0, budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    let error = initial
        .advance(1.0, budget(), || {
            seen += 1;
            if seen == checkpoints {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .unwrap_err();
    assert_eq!(error, DecayNeutronError::Fission(FissionError::Cancelled));
    assert_eq!(
        initial.sources()[0].inventory,
        original.sources()[0].inventory
    );
    assert_eq!(
        initial.neutrons().neutron_counts(),
        original.neutrons().neutron_counts()
    );
    let mut limited = budget();
    limited.decay.propagation.maximum_nodes = 2;
    limited.diffusion.limits.propagation = limited.decay.propagation;
    assert!(matches!(
        initial.advance(1.0, limited, || Ok(())),
        Err(DecayNeutronError::Budget(_))
    ));
    limited = budget();
    limited.decay.propagation.maximum_work = 1;
    limited.diffusion.limits.propagation = limited.decay.propagation;
    assert!(initial.advance(1.0, limited, || Ok(())).is_err());
    assert_eq!(
        initial.sources()[0].inventory,
        original.sources()[0].inventory
    );
}

#[test]
fn g0_model_grid_and_stable_emission_mismatches_refuse() {
    let initial = state(&[1.0, 0.0], &[100.0, 0.0], &[0.5, 0.0], 1.0, 0.0);
    let source = &initial.sources()[0];
    assert!(
        DecayNeutronEmission::new(
            "bad-stable".into(),
            Arc::clone(source.inventory.model()),
            Arc::clone(source.emission.groups()),
            vec![vec![0.5], vec![1.0]]
        )
        .is_err()
    );
    let mut wrong = source.clone();
    wrong.inventory = DecayInventory::new(chain(&[1.0, 0.0]), vec![100.0, 0.0]).unwrap();
    assert!(DecayNeutronState::new(initial.neutrons().clone(), vec![wrong]).is_err());
    let mut wrong_cell = source.clone();
    wrong_cell.cell = 1;
    assert!(DecayNeutronState::new(initial.neutrons().clone(), vec![wrong_cell]).is_err());
}
