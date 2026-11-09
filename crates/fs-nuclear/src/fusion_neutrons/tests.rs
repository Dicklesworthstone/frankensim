//! G1 exact constant-source/cascade solutions and G3 source/energy composition.
//! Authored coefficients are manufactured numerical fixtures, not blanket data.

use super::*;
use crate::diffusion::{
    DiffusionData, DiffusionLimits, DiffusionLink, DiffusionReaction, LeakageBoundary,
};
use crate::irradiation::NeutronGroups;
use fs_qty::{Amount, Energy, Time};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL;

#[path = "tests/secondary.rs"]
mod secondary;

#[path = "tests/fission.rs"]
mod fission;

#[path = "tests/precursors.rs"]
mod precursors;

#[path = "tests/induced_fission.rs"]
mod induced_fission;

pub(crate) use induced_fission::material_induced_fission_fixture;

const MEV: f64 = 1.602_176_634e-13;

fn data(energies_mev: &[f64], cells: usize) -> DiffusionData {
    let mut edges = vec![0.0];
    edges.extend(
        energies_mev
            .windows(2)
            .map(|pair| 0.5 * (pair[0] + pair[1]) * MEV),
    );
    edges.push(2.0 * energies_mev[energies_mev.len() - 1] * MEV);
    DiffusionData {
        data_id: "manufactured-fusion-neutron-transport".into(),
        groups: Arc::new(NeutronGroups::new("manufactured-group-edges".into(), edges).unwrap()),
        speeds_m_s: vec![1.0; energies_mev.len()],
        cell_volumes_m3: vec![1.0; cells],
        delayed_groups: vec![],
        reactions: vec![],
        links: vec![],
        leakage: vec![],
    }
}

fn reaction(
    cell: usize,
    group: usize,
    rate: f64,
    kind: DiffusionReactionKind,
) -> DiffusionReaction {
    DiffusionReaction {
        data_id: format!("manufactured-reaction-{cell}-{group}"),
        cell,
        incident_group: group,
        macroscopic_cross_section_m_inv: rate,
        kind,
    }
}

fn model(
    data: DiffusionData,
    energies_mev: &[f64],
    fractions: Vec<f64>,
) -> Arc<NeutronKineticModel> {
    Arc::new(
        NeutronKineticModel::new(
            "manufactured-kinetic-partition".into(),
            Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap()),
            energies_mev.iter().map(|value| value * MEV).collect(),
            fractions,
        )
        .unwrap(),
    )
}

fn emission(duration: f64, counts: [f64; 2], means_mev: [f64; 2]) -> FusionNeutronEmissionV1 {
    FusionNeutronEmissionV1::new(
        Time::new(duration),
        counts.map(|value| Amount::new(value / AVOGADRO_PER_MOL)),
        core::array::from_fn(|i| Energy::new(counts[i] * means_mev[i] * MEV)),
    )
    .unwrap()
}

fn budget() -> NeutronKineticBudget {
    NeutronKineticBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        energy_tolerance_j: 1e-20,
    }
}

fn close(actual: f64, expected: f64) {
    let tolerance = 2e-10 * expected.abs().max(1e-20);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn channel_projection_preserves_original_number_energy_and_duration() {
    let owner = model(data(&[1.0, 5.0, 15.0], 2), &[1.0, 5.0, 15.0], vec![]);
    let original = emission(0.25, [120.0, 80.0], [14.1, 2.45]);
    let source =
        FusionNeutronSource::project_interval_mean(owner.clone(), original.clone(), 1).unwrap();
    assert_eq!(source.emission(), &original);
    assert_eq!(source.cell_index(), 1);
    assert_eq!(&source.neutron_rates_s()[..3], &[0.0; 3]);
    close(source.neutron_rates_s().iter().sum::<f64>() * 0.25, 200.0);
    close(
        source
            .neutron_rates_s()
            .iter()
            .zip(owner.group_energies_j().iter().cycle())
            .map(|(rate, energy)| rate * energy)
            .sum::<f64>()
            * 0.25,
        original.total_kinetic_energy().value(),
    );
    assert!(source.number_projection_residual().abs() < 1e-12);
    assert!(source.energy_projection_residual_j().abs() < 1e-25);
    let initial = NeutronKineticState::new(owner, vec![0.0; 6]).unwrap();
    let (final_state, step) = initial
        .trial_with_fusion(0.25, &[source], budget(), || Ok(()))
        .unwrap();
    close(final_state.neutrons().neutron_counts().iter().sum(), 200.0);
    close(
        final_state.kinetic_energy_j().unwrap(),
        original.total_kinetic_energy().value(),
    );
    assert_eq!(step.released, EnergyChannels::default());
    assert_eq!(step.fusion_sources.len(), 1);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0; 6]);
}

#[test]
fn continuous_fusion_births_have_correct_capture_leakage_and_residence_time() {
    let mut supplied = data(&[14.1], 1);
    supplied
        .reactions
        .push(reaction(0, 0, 0.7, DiffusionReactionKind::Capture));
    supplied.leakage.push(LeakageBoundary {
        cell: 0,
        conductance_m2: vec![0.3],
    });
    let owner = model(supplied, &[14.1], vec![0.25]);
    let initial = NeutronKineticState::new(owner.clone(), vec![0.0]).unwrap();
    let duration = 2.0;
    let source = FusionNeutronSource::project_interval_mean(
        owner,
        emission(duration, [200.0, 0.0], [14.1, 0.0]),
        0,
    )
    .unwrap();
    let (next, step) = initial
        .trial_with_fusion(duration, &[source], budget(), || Ok(()))
        .unwrap();
    let remaining = 100.0 * -(-duration).exp_m1();
    let removed = 200.0 - remaining;
    close(next.neutrons().neutron_counts()[0], remaining);
    close(step.diffusion.captures, 0.7 * removed);
    close(step.diffusion.leaked_neutrons, 0.3 * removed);
    close(step.released.deposited_j, 0.25 * 0.7 * removed * 14.1 * MEV);
    close(step.released.radiation_j, 0.75 * 0.7 * removed * 14.1 * MEV);
    close(step.released.escaped_j, 0.3 * removed * 14.1 * MEV);
    close(
        step.final_kinetic_energy_j + step.released.total_j(),
        200.0 * 14.1 * MEV,
    );
    assert!(step.energy_roundoff_j < budget().energy_tolerance_j);
    // Depositing all 200 births at t=0 would leave only 200 exp(-2), which is wrong.
    assert!((remaining - 200.0 * (-duration).exp()).abs() > 50.0);
    let (after, shutdown) = next
        .trial_with_fusion(3.0, &[], budget(), || Ok(()))
        .unwrap();
    close(
        after.neutrons().neutron_counts()[0],
        remaining * (-3.0_f64).exp(),
    );
    close(shutdown.injected_kinetic_energy_j, 0.0);
    close(
        after.cumulative_cell_released()[0].total_j() + after.kinetic_energy_j().unwrap(),
        200.0 * 14.1 * MEV,
    );
}

#[test]
fn downscattering_and_capture_credit_only_their_actual_kinetic_loss() {
    let mut supplied = data(&[2.0, 14.0], 1);
    supplied.reactions = vec![
        reaction(
            0,
            1,
            2.0,
            DiffusionReactionKind::Scatter { outgoing_group: 0 },
        ),
        reaction(0, 0, 0.5, DiffusionReactionKind::Capture),
    ];
    let owner = model(supplied, &[2.0, 14.0], vec![0.75, 0.25]);
    let initial = NeutronKineticState::new(owner, vec![0.0, 0.0]).unwrap();
    let duration = 1.5;
    let rate = 80.0;
    let (next, step) = initial
        .trial_with_rates(duration, &[0.0, rate], budget(), || Ok(()))
        .unwrap();
    let high = rate * -(-2.0_f64 * duration).exp_m1() / 2.0;
    let low = rate
        * (-(-0.5_f64 * duration).exp_m1() / 0.5
            - ((-2.0_f64 * duration).exp() - (-0.5_f64 * duration).exp()) / (0.5 - 2.0));
    let scattered = rate * duration - high;
    let captured = scattered - low;
    close(next.neutrons().neutron_counts()[1], high);
    close(next.neutrons().neutron_counts()[0], low);
    close(step.diffusion.scatters, scattered);
    close(step.diffusion.captures, captured);
    close(
        step.released.deposited_j,
        (0.75 * 12.0 * scattered + 0.25 * 2.0 * captured) * MEV,
    );
    close(
        step.released.radiation_j,
        (0.25 * 12.0 * scattered + 0.75 * 2.0 * captured) * MEV,
    );
    assert_eq!(step.released.escaped_j, 0.0);
    close(
        next.kinetic_energy_j().unwrap() + step.released.total_j(),
        rate * duration * 14.0 * MEV,
    );
}

#[test]
fn spatial_energy_transfer_is_reciprocal_and_survives_restart() {
    let mut supplied = data(&[3.0], 2);
    supplied.cell_volumes_m3 = vec![1.0, 2.0];
    supplied.links.push(DiffusionLink {
        cell_a: 0,
        cell_b: 1,
        conductance_m2: vec![0.4],
    });
    supplied
        .reactions
        .push(reaction(1, 0, 0.2, DiffusionReactionKind::Capture));
    let owner = model(supplied, &[3.0], vec![1.0]);
    let initial = NeutronKineticState::new(owner, vec![50.0, 0.0]).unwrap();
    let (whole, whole_step) = initial
        .trial_with_rates(2.0, &[10.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let (half, first) = initial
        .trial_with_rates(1.0, &[10.0, 0.0], budget(), || Ok(()))
        .unwrap();
    let (split, second) = half
        .trial_with_rates(1.0, &[10.0, 0.0], budget(), || Ok(()))
        .unwrap();
    for (&actual, &expected) in split
        .neutrons()
        .neutron_counts()
        .iter()
        .zip(whole.neutrons().neutron_counts())
    {
        close(actual, expected);
    }
    close(
        first.released.deposited_j + second.released.deposited_j,
        whole_step.released.deposited_j,
    );
    close(
        split.kinetic_energy_j().unwrap() + split.cumulative_cell_released()[1].deposited_j,
        70.0 * 3.0 * MEV,
    );
    assert_eq!(whole_step.cell_released[0], EnergyChannels::default());
    assert!(whole_step.cell_released[1].deposited_j > 0.0);
}

#[test]
fn capture_radiation_and_neutron_leakage_are_never_automatic_heat() {
    for capture_fraction in [0.0, 1.0] {
        let mut supplied = data(&[2.45], 1);
        supplied
            .reactions
            .push(reaction(0, 0, 2.0, DiffusionReactionKind::Capture));
        let owner = model(supplied, &[2.45], vec![capture_fraction]);
        let initial = NeutronKineticState::new(owner, vec![20.0]).unwrap();
        let (next, step) = initial
            .trial_with_fusion(20.0, &[], budget(), || Ok(()))
            .unwrap();
        close(
            next.neutrons().neutron_counts()[0],
            20.0 * (-40.0_f64).exp(),
        );
        close(
            step.released.deposited_j,
            capture_fraction * 20.0 * 2.45 * MEV,
        );
        close(
            step.released.radiation_j,
            (1.0 - capture_fraction) * 20.0 * 2.45 * MEV,
        );
    }
    let mut supplied = data(&[2.45], 1);
    supplied.leakage.push(LeakageBoundary {
        cell: 0,
        conductance_m2: vec![2.0],
    });
    let owner = model(supplied, &[2.45], vec![]);
    let initial = NeutronKineticState::new(owner, vec![20.0]).unwrap();
    let (_, step) = initial
        .trial_with_fusion(20.0, &[], budget(), || Ok(()))
        .unwrap();
    assert_eq!(step.released.deposited_j, 0.0);
    assert_eq!(step.released.radiation_j, 0.0);
    close(step.released.escaped_j, 20.0 * 2.45 * MEV);
}

#[test]
fn unsupported_models_group_coverage_and_intervals_refuse_atomically() {
    let make =
        |supplied| Arc::new(DiffusionModel::new(supplied, DiffusionLimits::default()).unwrap());
    let mut upscatter = data(&[1.0, 3.0], 1);
    upscatter.reactions.push(reaction(
        0,
        0,
        1.0,
        DiffusionReactionKind::Scatter { outgoing_group: 1 },
    ));
    assert!(
        NeutronKineticModel::new(
            "upscatter".into(),
            make(upscatter),
            vec![MEV, 3.0 * MEV],
            vec![1.0]
        )
        .is_err()
    );
    let mut multiplying = data(&[3.0], 1);
    multiplying.reactions.push(reaction(
        0,
        0,
        1.0,
        DiffusionReactionKind::Fission {
            prompt_yields: vec![2.0],
            delayed_yields: vec![],
        },
    ));
    assert!(
        NeutronKineticModel::new(
            "fission".into(),
            make(multiplying),
            vec![3.0 * MEV],
            vec![1.0]
        )
        .is_err()
    );
    assert!(
        NeutronKineticModel::new(
            "bad-energy".into(),
            make(data(&[3.0], 1)),
            vec![9.0 * MEV],
            vec![]
        )
        .is_err()
    );
    let owner = model(data(&[3.0, 5.0], 1), &[3.0, 5.0], vec![]);
    assert!(
        FusionNeutronSource::project_interval_mean(
            owner.clone(),
            emission(1.0, [1.0, 0.0], [14.1, 0.0]),
            0
        )
        .is_err()
    );
    assert!(
        FusionNeutronSource::project_interval_mean(
            owner.clone(),
            emission(1.0, [1.0, 0.0], [4.0, 0.0]),
            1
        )
        .is_err()
    );
    let source = FusionNeutronSource::project_interval_mean(
        owner.clone(),
        emission(1.0, [1.0, 0.0], [4.0, 0.0]),
        0,
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![1.0, 2.0]).unwrap();
    assert!(
        initial
            .trial_with_fusion(0.5, &[source], budget(), || Ok(()))
            .is_err()
    );
    assert!(
        initial
            .trial_with_rates(1.0, &[0.0], budget(), || Ok(()))
            .is_err()
    );
    assert!(
        initial
            .trial_with_rates(1.0, &[0.0, 0.0], budget(), || Err(FissionError::Cancelled))
            .is_err()
    );
    assert_eq!(initial.neutrons().neutron_counts(), &[1.0, 2.0]);
    assert!(
        initial
            .cumulative_cell_released()
            .iter()
            .all(|channels| *channels == EnergyChannels::default())
    );
}

#[test]
fn zero_emission_and_zero_duration_retain_original_state() {
    let owner = model(data(&[3.0], 1), &[3.0], vec![]);
    let source = FusionNeutronSource::project_interval_mean(
        owner.clone(),
        emission(0.0, [0.0; 2], [0.0; 2]),
        0,
    )
    .unwrap();
    let initial = NeutronKineticState::new(owner, vec![12.0]).unwrap();
    let (next, step) = initial
        .trial_with_fusion(0.0, &[source], budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        next.neutrons().neutron_counts(),
        initial.neutrons().neutron_counts()
    );
    assert_eq!(step.injected_kinetic_energy_j, 0.0);
    assert_eq!(step.released, EnergyChannels::default());
}
