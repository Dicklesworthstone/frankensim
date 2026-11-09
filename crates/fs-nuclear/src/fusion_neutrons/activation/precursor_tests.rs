//! Manufactured finite radioactive families with independent population/energy histories.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct,
    BurnupReaction,
};
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;
use crate::diffusion::{
    DelayedGroup, DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionReaction,
    DiffusionReactionKind,
};
use crate::fusion_neutrons::NeutronReactionEnergy;
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}

fn budget() -> FusionNeutronActivationBudget {
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

pub(crate) fn material_precursor_fixture(
    parent_count: f64,
    decay_radiation_j: f64,
) -> Result<FusionNeutronActivationState, FusionNeutronActivationError> {
    let groups = Arc::new(NeutronGroups::new(
        "manufactured-delayed-energy-grid".into(), vec![0.0, 4.0],
    ).unwrap());
    let family = DelayedGroup {
        data_id: "manufactured-one-second-inverse-family".into(),
        decay_constant_s_inv: 1.0,
        spectrum: vec![1.0],
    };
    let nuclear = Arc::new(BurnupModel::new(
        BurnupData {
            data_id: "manufactured-reserved-decay-stock".into(),
            groups: Arc::clone(&groups),
            temperature_min_k: 250.0,
            temperature_max_k: 600.0,
            nuclides: vec![
                BurnupNuclide {
                    id: NuclideId::new(10, 25, 0).unwrap(),
                    law: DecayLawV1::from_decay_constant(Frequency::new(0.5)).unwrap(),
                },
                BurnupNuclide {
                    id: NuclideId::new(11, 25, 0).unwrap(),
                    law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
                },
            ],
            reactions: vec![BurnupReaction {
                data_id: "manufactured-decay-with-explicit-reservation".into(),
                parent: 0,
                rate: BurnupRate::Decay { branch_fraction: 1.0 },
                products: vec![BurnupProduct { nuclide: 1, yield_per_event: 1.0 }],
                energy_per_event: EnergyChannels {
                    deposited_j: 1.0,
                    radiation_j: decay_radiation_j,
                    escaped_j: 0.5,
                },
            }],
        },
        BurnupLimits::default(),
    ).unwrap());
    let diffusion = Arc::new(DiffusionModel::new(
        DiffusionData {
            data_id: "manufactured-delayed-neutron-capture".into(),
            groups,
            speeds_m_s: vec![1.0],
            cell_volumes_m3: vec![0.01],
            delayed_groups: vec![family.clone()],
            reactions: vec![DiffusionReaction {
                data_id: "manufactured-fixed-background-capture".into(),
                cell: 0,
                incident_group: 0,
                macroscopic_cross_section_m_inv: 0.25,
                kind: DiffusionReactionKind::Capture,
            }],
            links: vec![],
            leakage: vec![],
        },
        DiffusionLimits::default(),
    ).unwrap());
    let kinetic = Arc::new(NeutronKineticModel::new_with_precursor_storage(
        "manufactured-explicit-future-neutron-storage".into(),
        Arc::clone(&diffusion),
        vec![2.0],
        vec![0.25],
        vec![NeutronReactionEnergy::default()],
    ).unwrap());
    let emission = Arc::new(BurnupDecayNeutronEmission::new_with_precursors(
        "manufactured-immediate-and-delayed-births".into(),
        Arc::clone(&nuclear),
        vec![vec![0.5], vec![0.0]],
        vec![family],
        vec![vec![1.0], vec![0.0]],
    ).unwrap());
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0], vec![0.0]).unwrap(),
        vec![BurnupInventory::new(nuclear, vec![parent_count, 0.0]).unwrap()],
        vec![],
    ).unwrap().with_decay_neutrons(vec![emission]).unwrap();
    FusionNeutronActivationState::new(kinetic, spatial)
}

/// Independent scalar Bateman solution for lambda=.5, mu=1, capture=.25.
pub(crate) fn expected_history(parent_count: f64, time: f64) -> (f64, f64, f64, f64) {
    let parent = parent_count * (-0.5 * time).exp();
    let decays = parent_count * -(-0.5 * time).exp_m1();
    let precursors = parent_count * ((-0.5 * time).exp() - (-time).exp());
    let a = ((-0.5 * time).exp() - (-0.25 * time).exp()) / -0.25;
    let b = ((-time).exp() - (-0.25 * time).exp()) / -0.75;
    let neutrons = parent_count * (1.25 * a - b);
    (parent, decays, precursors, neutrons)
}

#[test]
fn g1_g3_reserved_births_decay_and_capture_share_one_finite_energy_history() {
    let initial = material_precursor_fixture(100.0, 4.0).unwrap();
    let (after, step) = initial.trial_with_rates(2.0, &[300.0], &[0.0], budget(), || Ok(())).unwrap();
    let (parent, decays, precursors, neutrons) = expected_history(100.0, 2.0);
    let captures = 1.5 * decays - precursors - neutrons;
    assert!(captures > 0.0 && precursors > 0.0);
    close(after.inventories()[0].populations()[0], parent, 1e-9);
    close(after.inventories()[0].populations()[1], decays, 1e-9);
    close(after.neutrons().precursor_counts()[0], precursors, 1e-9);
    close(after.neutrons().neutron_counts()[0], neutrons, 1e-9);
    close(step.decay_precursors, decays, 1e-9);
    close(step.decay_neutrons, 0.5 * decays, 1e-9);
    close(step.reserved_decay_precursor_energy_j, 2.0 * decays, 1e-8);
    close(step.transferred_decay_neutron_energy_j, decays, 1e-8);
    close(step.kinetic.precursor_conversion_energy_j, 2.0 * (decays - precursors), 1e-8);
    close(step.released.deposited_j, decays + 0.5 * captures, 1e-8);
    close(step.released.radiation_j, decays + 1.5 * captures, 1e-8);
    close(step.released.escaped_j, 0.5 * decays, 1e-8);
    close(after.precursor_energy_j().unwrap(), 2.0 * precursors, 1e-8);
    close(after.stored_energy_j().unwrap() + step.released.total_j(), 5.5 * decays, 1e-8);
    assert_eq!(step.external_injected_precursors, 0.0);
    assert_eq!(step.external_injected_precursor_energy_j, 0.0);
    assert_eq!(step.kinetic.exported_precursor_energy_j, 0.0);
    let (middle, _) = initial.trial_with_rates(0.7, &[300.0], &[0.0], budget(), || Ok(())).unwrap();
    let (resumed, _) = middle.trial_with_rates(1.3, &[300.0], &[0.0], budget(), || Ok(())).unwrap();
    close(resumed.stored_energy_j().unwrap(), after.stored_energy_j().unwrap(), 1e-8);
    close(resumed.cell_released(0).unwrap().total_j(), step.released.total_j(), 1e-8);
    assert!(Arc::ptr_eq(resumed.inventories()[0].model(), initial.inventories()[0].model()));
}

#[test]
fn g0_g4_external_family_energy_is_separate_and_late_refusal_keeps_checkpoint() {
    assert!(matches!(material_precursor_fixture(100.0, 2.9), Err(FusionNeutronActivationError::Invalid(_))));
    let initial = material_precursor_fixture(0.0, 4.0).unwrap();
    let mut polls = 0;
    let (after, report) = initial.trial_with_sources(
        2.0, &[300.0], &[0.0], &[3.0], budget(), || { polls += 1; Ok(()) },
    ).unwrap();
    close(report.external_injected_precursors, 6.0, 1e-12);
    close(report.external_injected_precursor_energy_j, 12.0, 1e-12);
    assert_eq!(report.reserved_decay_precursor_energy_j, 0.0);
    assert_eq!(report.nuclear_generated_energy_j, 0.0);
    close(after.neutrons().precursor_counts()[0], 3.0 * -(-2.0_f64).exp_m1(), 1e-9);
    close(after.stored_energy_j().unwrap() + report.released.total_j(), 12.0, 1e-8);
    let mut remaining = polls - 1;
    assert!(initial.trial_with_sources(2.0, &[300.0], &[0.0], &[3.0], budget(), || {
        if remaining == 0 { Err(FissionError::Cancelled) } else { remaining -= 1; Ok(()) }
    }).is_err());
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0]);
    assert_eq!(initial.neutrons().precursor_counts(), &[0.0]);
    assert_eq!(initial.cell_released(0).unwrap(), EnergyChannels::default());
    let (retry, _) = initial.trial_with_sources(2.0, &[300.0], &[0.0], &[3.0], budget(), || Ok(())).unwrap();
    close(retry.stored_energy_j().unwrap(), after.stored_energy_j().unwrap(), 1e-10);
}
