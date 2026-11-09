//! G1/G3 disjoint source, neutron storage, heat and leakage accounts.
use super::super::tests::{budget as spatial_budget, close, state};
use super::*;

fn energy_state() -> DecayNeutronEnergyState {
    let spatial = state(&[1.0, 0.0], &[100.0, 0.0], &[0.5, 0.0], 1.2, 0.8);
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-decay-neutron-energy".into(),
            Arc::clone(spatial.neutrons().model()),
            vec![2.0],
            vec![0.75],
        )
        .unwrap(),
    );
    DecayNeutronEnergyState::new(model, spatial).unwrap()
}
fn budget() -> DecayNeutronEnergyBudget {
    DecayNeutronEnergyBudget {
        spatial: spatial_budget(),
        energy_tolerance_j: 1e-7,
    }
}

#[test]
fn g1_source_radiation_is_transferred_once_then_absorbed_or_leaked() {
    let initial = energy_state();
    let (next, step) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let decays = 100.0 * -(-1.0_f64).exp_m1();
    let neutrons = 50.0 * ((-1.0_f64).exp() - (-2.0_f64).exp());
    let removed = 0.5 * decays - neutrons;
    let capture_energy = 2.0 * 0.6 * removed;
    close(step.transferred_neutron_energy_j, decays);
    close(step.kinetic.injected_kinetic_energy_j, decays);
    close(step.kinetic.final_kinetic_energy_j, 2.0 * neutrons);
    close(
        step.released.deposited_j,
        2.0 * decays + 0.75 * capture_energy,
    );
    close(
        step.released.radiation_j,
        2.0 * decays + 0.25 * capture_energy,
    );
    close(step.released.escaped_j, decays + 2.0 * 0.4 * removed);
    close(step.decay_released.total_j(), 6.0 * decays);
    close(
        step.released.total_j() + step.kinetic.final_kinetic_energy_j,
        step.decay_released.total_j(),
    );
    close(step.energy_balance_residual_j, 0.0);
    close(
        next.cell_released(0).unwrap().total_j(),
        step.released.total_j(),
    );
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
}

#[test]
fn g1_two_sources_transport_and_deposit_energy_in_their_original_cells() {
    use crate::diffusion::{
        DiffusionData, DiffusionLink, DiffusionModel, DiffusionReaction, DiffusionReactionKind,
        DiffusionState,
    };

    let populations = [100.0, 40.0];
    let yields = [0.5, 0.25];
    let heat_fractions = [0.75, 0.25];
    let link_rate = 0.5;
    let templates = [
        state(
            &[1.0, 0.0],
            &[populations[0], 0.0],
            &[yields[0], 0.0],
            1.0,
            0.0,
        ),
        state(
            &[1.0, 0.0],
            &[populations[1], 0.0],
            &[yields[1], 0.0],
            1.0,
            0.0,
        ),
    ];
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-two-cell-decay".into(),
                groups: Arc::clone(templates[0].neutrons().model().groups()),
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0, 1.0],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|cell| DiffusionReaction {
                        data_id: format!("capture-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: 1.0,
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![DiffusionLink {
                    cell_a: 0,
                    cell_b: 1,
                    conductance_m2: vec![link_rate],
                }],
                leakage: vec![],
            },
            spatial_budget().diffusion.limits,
        )
        .unwrap(),
    );
    let sources = templates
        .iter()
        .enumerate()
        .map(|(cell, template)| {
            let mut source = template.sources()[0].clone();
            source.cell = cell;
            source
        })
        .collect();
    let spatial = DecayNeutronState::new(
        DiffusionState::new(Arc::clone(&diffusion), vec![0.0, 0.0], vec![]).unwrap(),
        sources,
    )
    .unwrap();
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-two-cell-decay-energy".into(),
            diffusion,
            vec![2.0],
            heat_fractions.to_vec(),
        )
        .unwrap(),
    );
    let initial = DecayNeutronEnergyState::new(model, spatial).unwrap();
    let duration = 1.3_f64;
    let (next, step) = initial.trial(duration, budget(), || Ok(())).unwrap();

    // With lambda = capture = 1, the symmetric and antisymmetric modes are
    // S = Q_sum t exp(-t), D = Q_diff (exp(-t)-exp(-(1+2k)t))/(2k).
    // Their integrals independently determine captures and both link counts.
    let tail = (-duration).exp();
    let fast_rate = 1.0 + 2.0 * link_rate;
    let fast_tail = (-fast_rate * duration).exp();
    let q_sum = populations[0] * yields[0] + populations[1] * yields[1];
    let q_difference = populations[0] * yields[0] - populations[1] * yields[1];
    let neutron_sum = q_sum * duration * tail;
    let neutron_difference = q_difference * (tail - fast_tail) / (2.0 * link_rate);
    let capture_sum = q_sum * (1.0 - (1.0 + duration) * tail);
    let capture_difference =
        q_difference * ((1.0 - tail) - (1.0 - fast_tail) / fast_rate) / (2.0 * link_rate);
    let captures = [
        0.5 * (capture_sum + capture_difference),
        0.5 * (capture_sum - capture_difference),
    ];
    let neutrons = [
        0.5 * (neutron_sum + neutron_difference),
        0.5 * (neutron_sum - neutron_difference),
    ];
    for cell in 0..2 {
        let decays = populations[cell] * (1.0 - tail);
        let source = &next.spatial().sources()[cell];
        close(source.inventory.populations()[0], populations[cell] * tail);
        close(source.inventory.populations()[1], decays);
        close(source.inventory.cumulative_decays()[0], decays);
        close(
            source.inventory.cumulative_released().total_j(),
            6.0 * decays,
        );
        close(
            step.spatial.emitted_neutrons[cell][0],
            yields[cell] * decays,
        );
        close(
            step.spatial.neutrons.injected_neutron_counts[cell],
            yields[cell] * decays,
        );
        close(
            next.spatial().neutrons().neutron_counts()[cell],
            neutrons[cell],
        );
        close(step.spatial.neutrons.reaction_counts[cell], captures[cell]);
        close(
            step.spatial.neutrons.integrated_fluence_m2[cell],
            captures[cell],
        );

        let heat = 2.0 * decays + heat_fractions[cell] * 2.0 * captures[cell];
        let radiation = (3.0 - 2.0 * yields[cell]) * decays
            + (1.0 - heat_fractions[cell]) * 2.0 * captures[cell];
        close(step.cell_released[cell].deposited_j, heat);
        close(step.cell_released[cell].radiation_j, radiation);
        close(step.cell_released[cell].escaped_j, decays);
        close(next.cell_released(cell).unwrap().deposited_j, heat);
        assert_eq!(
            initial.spatial().sources()[cell].inventory.populations(),
            &[populations[cell], 0.0]
        );
        assert_eq!(initial.cell_released(cell), Some(EnergyChannels::default()));
    }
    close(
        step.spatial.neutrons.link_a_to_b_counts[0],
        link_rate * captures[0],
    );
    close(
        step.spatial.neutrons.link_b_to_a_counts[0],
        link_rate * captures[1],
    );
    close(step.spatial.neutrons.captures, capture_sum);
    close(
        step.transferred_neutron_energy_j,
        2.0 * q_sum * (1.0 - tail),
    );
    close(step.kinetic.final_kinetic_energy_j, 2.0 * neutron_sum);
    close(
        step.decay_released.total_j(),
        6.0 * populations.iter().sum::<f64>() * (1.0 - tail),
    );
    close(
        step.released.total_j() + step.kinetic.final_kinetic_energy_j,
        step.decay_released.total_j(),
    );
    close(step.energy_balance_residual_j, 0.0);
}

#[test]
fn g3_restart_preserves_storage_disjoint_exports_and_original_source_accounts() {
    let initial = energy_state();
    let (whole, total) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let (middle, first) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let (split, second) = middle.clone().trial(1.0, budget(), || Ok(())).unwrap();
    for (a, b) in values(whole.cell_released(0).unwrap())
        .into_iter()
        .zip(values(split.cell_released(0).unwrap()))
    {
        close(a, b);
    }
    close(
        whole.kinetic().kinetic_energy_j().unwrap(),
        split.kinetic().kinetic_energy_j().unwrap(),
    );
    close(
        first.transferred_neutron_energy_j + second.transferred_neutron_energy_j,
        total.transferred_neutron_energy_j,
    );
    close(
        split.spatial().sources()[0]
            .inventory
            .cumulative_released()
            .total_j(),
        total.decay_released.total_j(),
    );
    assert!(
        split.cell_released(0).unwrap().total_j()
            < split.spatial().sources()[0]
                .inventory
                .cumulative_released()
                .total_j()
    );
}

#[test]
fn g0_neutrons_cannot_reuse_local_heat_or_a_different_neutron_model() {
    let excessive = state(&[1.0, 0.0], &[100.0, 0.0], &[2.0, 0.0], 1.2, 0.8);
    let model = Arc::new(
        NeutronKineticModel::new(
            "insufficient-outward-energy".into(),
            Arc::clone(excessive.neutrons().model()),
            vec![2.0],
            vec![1.0],
        )
        .unwrap(),
    );
    // Four joules/neutron source event cannot be funded by three outward
    // joules, even though the total decay Q including heat and escape is six.
    assert!(matches!(
        DecayNeutronEnergyState::new(model, excessive),
        Err(DecayNeutronEnergyError::Invalid(_))
    ));
    let original = energy_state();
    let equivalent_but_separate = energy_state();
    assert!(
        DecayNeutronEnergyState::new(
            Arc::clone(original.model()),
            equivalent_but_separate.spatial().clone()
        )
        .is_err()
    );
}

#[test]
fn g4_energy_refusal_and_late_cancellation_publish_no_physical_history() {
    let mut initial = energy_state();
    let original_source = initial.spatial().sources()[0].inventory.clone();
    let mut limited = budget();
    limited.energy_tolerance_j = 0.0;
    assert!(initial.advance(1.0, limited, || Ok(())).is_err());
    assert_eq!(initial.spatial().sources()[0].inventory, original_source);
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
    let mut checkpoints = 0;
    initial
        .trial(1.0, budget(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut seen = 0;
    let result = initial.advance(1.0, budget(), || {
        seen += 1;
        if seen == checkpoints {
            Err(FissionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(
        result.unwrap_err(),
        DecayNeutronEnergyError::Fission(FissionError::Cancelled)
    );
    assert_eq!(initial.spatial().sources()[0].inventory, original_source);
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
}
