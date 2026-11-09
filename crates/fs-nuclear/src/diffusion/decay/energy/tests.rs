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
