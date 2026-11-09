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

#[test]
fn g1_g3_binding_and_restart_preserve_original_precursor_stock_with_radioactive_sources() {
    use crate::diffusion::{DelayedGroup, DiffusionModel, DiffusionState};
    use crate::fusion_neutrons::NeutronReactionEnergy;

    let template = state(&[1.0, 0.0], &[100.0, 0.0], &[0.5, 0.0], 0.25, 0.0);
    let mut data = template.neutrons().model().data().clone();
    data.delayed_groups = vec![DelayedGroup {
        data_id: "manufactured-resident-family".into(),
        decay_constant_s_inv: 2.0,
        spectrum: vec![1.0],
    }];
    let native = Arc::new(DiffusionModel::new(data, spatial_budget().diffusion.limits).unwrap());
    let model = Arc::new(
        NeutronKineticModel::new_with_precursor_storage(
            "manufactured-preserved-family-storage".into(),
            Arc::clone(&native),
            vec![2.0],
            vec![0.25],
            vec![NeutronReactionEnergy::default()],
        )
        .unwrap(),
    );
    let spatial = DecayNeutronState::new(
        DiffusionState::new(native, vec![0.0], vec![10.0]).unwrap(),
        template.sources().to_vec(),
    )
    .unwrap();
    let initial = DecayNeutronEnergyState::new(Arc::clone(&model), spatial).unwrap();
    assert_eq!(initial.kinetic().neutrons().precursor_counts(), &[10.0]);
    close(initial.kinetic().stored_energy_j().unwrap(), 20.0);
    let time = 1.25_f64;
    let response = |rate: f64| ((-rate * time).exp() - (-0.25 * time).exp()) / (0.25 - rate);
    let decays = 100.0 * -(-time).exp_m1();
    let families = 10.0 * (-2.0 * time).exp();
    let neutrons = 50.0 * response(1.0) + 20.0 * response(2.0);
    let captures = 0.5 * decays + 10.0 - families - neutrons;
    let (whole, step) = initial.trial(time, budget(), || Ok(())).unwrap();
    close(whole.spatial().neutrons().precursor_counts()[0], families);
    close(whole.kinetic().neutrons().precursor_counts()[0], families);
    close(whole.spatial().neutrons().neutron_counts()[0], neutrons);
    close(step.released.deposited_j, 2.0 * decays + 0.5 * captures);
    close(step.released.radiation_j, 2.0 * decays + 1.5 * captures);
    close(step.released.escaped_j, decays);
    close(
        whole.kinetic().stored_energy_j().unwrap() + step.released.total_j(),
        20.0 + 6.0 * decays,
    );
    assert_eq!(step.kinetic.induced_precursor_energy_j, 0.0);
    assert_eq!(step.kinetic.injected_precursor_energy_j, 0.0);

    let (middle, first) = initial.trial(0.5, budget(), || Ok(())).unwrap();
    let rebound = DecayNeutronEnergyState::new(model, middle.spatial().clone()).unwrap();
    assert_eq!(
        rebound.kinetic().neutrons().precursor_counts(),
        middle.spatial().neutrons().precursor_counts(),
    );
    let (split, second) = rebound.trial(0.75, budget(), || Ok(())).unwrap();
    close(
        split.kinetic().stored_energy_j().unwrap(),
        whole.kinetic().stored_energy_j().unwrap(),
    );
    close(
        first.released.total_j() + second.released.total_j(),
        step.released.total_j(),
    );
    close(
        split.spatial().sources()[0].inventory.cumulative_decays()[0],
        decays,
    );
}

/// Shared manufactured coexistence source for the existing material owner.
/// Fast40, initialfamily10, independent radioactive stock60 atlambda.5/s.
/// Group energies[2,4]J; fissionrate1, promptQ6, future reserve1J/fission;
/// familylambda2, slowcapture.25/s, both reaction heatfractions.25, volume.01.
pub(crate) fn material_induced_decay_fixture(outflow_rate_s: f64) -> DecayNeutronEnergyState {
    use crate::diffusion::DiffusionState;
    use crate::diffusion::decay::{DecayNeutronEmission, DecayNeutronSource};

    let native = crate::fusion_neutrons::material_induced_fission_fixture(40.0, outflow_rate_s);
    let original = state(&[0.5, 0.0], &[60.0, 0.0], &[0.5, 0.0], 0.25, 0.0);
    let source = &original.sources()[0];
    let emission = Arc::new(
        DecayNeutronEmission::new(
            "manufactured-separate-radioactive-neutrons".into(),
            Arc::clone(source.inventory.model()),
            Arc::clone(native.model().diffusion().groups()),
            vec![vec![0.5, 0.0], vec![0.0; 2]],
        )
        .unwrap(),
    );
    let spatial = DecayNeutronState::new(
        DiffusionState::new(
            Arc::clone(native.model().diffusion()),
            vec![0.0, 40.0],
            vec![10.0],
        )
        .unwrap(),
        vec![DecayNeutronSource {
            cell: 0,
            inventory: source.inventory.clone(),
            emission,
        }],
    )
    .unwrap();
    DecayNeutronEnergyState::new_with_unrepresented_fission_precursors(
        Arc::clone(native.model()),
        spatial,
    )
    .unwrap()
}

#[test]
fn g1_g3_g4_radioactive_and_funded_background_histories_close_without_replaying_reserves() {
    let initial = material_induced_decay_fixture(0.6);
    assert!(
        DecayNeutronEnergyState::new(Arc::clone(initial.model()), initial.spatial().clone(),)
            .is_err()
    );
    let time = 1.25_f64;
    let removal = 2.6;
    let response = |rate: f64| ((-rate * time).exp() - (-0.25 * time).exp()) / (0.25 - rate);
    let fissions = 40.0 * -(-time).exp_m1();
    let decays = 60.0 * -(-0.5 * time).exp_m1();
    let families = 10.0 * (-removal * time).exp()
        + 20.0 * ((-time).exp() - (-removal * time).exp()) / (removal - 1.0);
    let residence = 10.0 * -(-removal * time).exp_m1() / removal
        + 20.0 / (removal - 1.0) * (-(-time).exp_m1() + (-removal * time).exp_m1() / removal);
    let neutrons = 40.0 * response(1.0)
        + 40.0 / (removal - 1.0) * (response(1.0) - response(removal))
        + 20.0 * response(removal)
        + 15.0 * response(0.5);
    let captures = fissions + 2.0 * residence + 0.5 * decays - neutrons;
    let heat = 2.0 * fissions + 0.5 * captures + 2.0 * decays;
    let radiation = 6.0 * fissions + 1.5 * captures + 2.0 * decays;
    let (whole, step) = initial.trial(time, budget(), || Ok(())).unwrap();
    close(
        whole.spatial().neutrons().neutron_counts()[1],
        40.0 * (-time).exp(),
    );
    close(whole.spatial().neutrons().neutron_counts()[0], neutrons);
    close(whole.kinetic().neutrons().precursor_counts()[0], families);
    close(step.spatial.neutrons.fissions, fissions);
    close(step.spatial.neutrons.captures, captures);
    close(step.transferred_neutron_energy_j, decays);
    close(step.kinetic.induced_precursor_energy_j, fissions);
    close(step.kinetic.nuclear_reaction_energy_j, 6.0 * fissions);
    close(step.kinetic.exported_precursor_energy_j, 1.2 * residence);
    close(step.released.deposited_j, heat);
    close(step.released.radiation_j, radiation);
    close(step.released.escaped_j, decays);
    close(
        whole.kinetic().stored_energy_j().unwrap()
            + step.released.total_j()
            + step.kinetic.exported_precursor_energy_j,
        180.0 + 6.0 * decays + 7.0 * fissions,
    );
    assert_eq!(step.kinetic.injected_precursor_energy_j, 0.0);

    let (middle, first) = initial.trial(0.5, budget(), || Ok(())).unwrap();
    let mut polls = 0;
    let (split, second) = middle
        .trial(0.75, budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    close(
        split.kinetic().stored_energy_j().unwrap(),
        whole.kinetic().stored_energy_j().unwrap(),
    );
    close(
        first.released.total_j() + second.released.total_j(),
        step.released.total_j(),
    );
    close(
        split
            .kinetic()
            .cumulative_induced_precursor_energy_j()
            .unwrap(),
        fissions,
    );
    close(
        split
            .kinetic()
            .cumulative_exported_precursor_energy_j()
            .unwrap(),
        1.2 * residence,
    );
    let before = middle
        .kinetic()
        .cumulative_induced_precursor_energy_j()
        .unwrap();
    let populations = middle.spatial().neutrons().precursor_counts().to_vec();
    let mut seen = 0;
    assert!(
        middle
            .trial(0.75, budget(), || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(
        middle
            .kinetic()
            .cumulative_induced_precursor_energy_j()
            .unwrap(),
        before
    );
    assert_eq!(middle.spatial().neutrons().precursor_counts(), populations);
    let (retry, _) = middle.trial(0.75, budget(), || Ok(())).unwrap();
    close(
        retry.kinetic().stored_energy_j().unwrap(),
        split.kinetic().stored_energy_j().unwrap(),
    );
}
