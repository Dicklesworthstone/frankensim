//! Finite emitters compete for the same stock; energies use actual events.

use super::*;
use crate::diffusion::LeakageBoundary;
use crate::diffusion::burnup::decay::BurnupDecayNeutronEmission;

fn emitter(
    branches: &[(f64, f64)],
    capture_radiation: f64,
) -> Result<FusionNeutronActivationState, FusionNeutronActivationError> {
    emitter_with_yield(branches, capture_radiation, 1.0)
}

fn emitter_with_yield(
    branches: &[(f64, f64)],
    capture_radiation: f64,
    neutron_yield: f64,
) -> Result<FusionNeutronActivationState, FusionNeutronActivationError> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-finite-emitter-grid".into(), vec![0.0, 4.0]).unwrap(),
    );
    let mut reactions = vec![BurnupReaction {
        data_id: "manufactured-emitter-capture".into(),
        parent: 0,
        rate: BurnupRate::Capture {
            cross_sections_m2: vec![0.02],
        },
        products: vec![BurnupProduct {
            nuclide: 2,
            yield_per_event: 1.0,
        }],
        energy_per_event: EnergyChannels {
            deposited_j: 1.0,
            radiation_j: capture_radiation,
            escaped_j: 0.0,
        },
    }];
    reactions.extend(
        branches
            .iter()
            .enumerate()
            .map(|(index, &(fraction, radiation))| BurnupReaction {
                data_id: format!("manufactured-emitting-branch-{index}"),
                parent: 0,
                rate: BurnupRate::Decay {
                    branch_fraction: fraction,
                },
                products: vec![BurnupProduct {
                    nuclide: 1,
                    yield_per_event: 1.0,
                }],
                energy_per_event: EnergyChannels {
                    deposited_j: 0.5,
                    radiation_j: radiation,
                    escaped_j: 0.25,
                },
            }),
    );
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-emitter-stock".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 250.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    nuclide(10, 25, 1.0),
                    nuclide(11, 25, 0.0),
                    nuclide(10, 26, 0.0),
                ],
                reactions,
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-self-irradiating-emitter".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                delayed_groups: vec![],
                reactions: vec![DiffusionReaction {
                    data_id: "manufactured-bound-emitter-capture".into(),
                    cell: 0,
                    incident_group: 0,
                    macroscopic_cross_section_m_inv: 2.0,
                    kind: DiffusionReactionKind::Capture,
                }],
                links: vec![],
                leakage: vec![LeakageBoundary {
                    cell: 0,
                    conductance_m2: vec![0.3],
                }],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new(
            "manufactured-two-joule-neutrons".into(),
            Arc::clone(&diffusion),
            vec![2.0],
            vec![1.0],
        )
        .unwrap(),
    );
    let emission = Arc::new(
        BurnupDecayNeutronEmission::new(
            "manufactured-one-neutron-per-decay".into(),
            Arc::clone(&nuclear),
            vec![vec![neutron_yield], vec![0.0], vec![0.0]],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, vec![100.0, 0.0, 0.0]).unwrap()],
        vec![ReactionBinding {
            diffusion_reaction: 0,
            burnup_reaction: 0,
        }],
    )
    .unwrap()
    .with_decay_neutrons(vec![emission])
    .unwrap();
    FusionNeutronActivationState::new(model, spatial)
}

#[test]
fn g0_rare_emission_charges_neutron_roundoff_in_neutrons_not_parent_decays() {
    let initial = emitter_with_yield(&[(0.1, 2.0), (0.2, 2.0), (0.7, 2.0)], 0.0, 1e-12).unwrap();
    let mut rare_budget = budget();
    rare_budget.spatial.diffusion.population_tolerance = 1e-20;
    for duration in [0.13, 0.27, 0.41] {
        let (_, report) = initial
            .trial_with_rates(duration, &[300.0], &[0.0], rare_budget, || Ok(()))
            .unwrap();
        let d = report.spatial.burnup[0].decays;
        assert!(d > 10.0);
        close(report.decay_neutrons, 1e-12 * d, 1e-23);
        assert!(report.injection_number_roundoff <= 1e-20);
        close(report.transferred_decay_neutron_energy_j, 2e-12 * d, 1e-23);
    }
}

#[test]
fn g1_finite_decay_neutrons_transfer_energy_once_and_keep_external_source_separate() {
    let initial = emitter(&[(1.0, 2.0)], 0.0).unwrap();
    let (next, report) = initial
        .trial_with_rates(0.4, &[300.0], &[3.0], budget(), || Ok(()))
        .unwrap();
    let nuclear = &report.spatial.burnup[0];
    let d = nuclear.decays;
    let c = nuclear.captures;
    let lost = report.spatial.neutrons.leaked_neutrons;
    assert!(d > 0.0 && c > 0.0 && lost > 0.0);
    close(next.inventories()[0].populations()[0] + d + c, 100.0, 1e-8);
    close(report.decay_neutrons, d, 1e-8);
    close(report.transferred_decay_neutron_energy_j, 2.0 * d, 1e-8);
    close(report.external_injected_neutrons, 1.2, 1e-12);
    close(report.external_injected_kinetic_energy_j, 2.4, 1e-12);
    close(
        report.kinetic.injected_kinetic_energy_j,
        2.0 * d + 2.4,
        1e-8,
    );
    close(report.nuclear_released.radiation_j, 2.0 * d, 1e-8);
    close(report.released.deposited_j, 0.5 * d + 3.0 * c, 1e-8);
    close(report.released.radiation_j, 0.0, 1e-8);
    close(report.released.escaped_j, 0.25 * d + 2.0 * lost, 1e-8);
    close(
        next.kinetic_energy_j().unwrap() + report.released.total_j(),
        2.4 + 2.75 * d + c,
        1e-8,
    );
    assert_eq!(report.kinetic.diffusion, report.spatial.neutrons);
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
}

#[test]
fn g0_only_original_rate_weighted_decay_radiation_can_fund_parent_emission() {
    assert!(matches!(
        emitter(&[(1.0, 1.0)], 1e6),
        Err(FusionNeutronActivationError::Invalid(_))
    ));
    let fractions = [(0.25, 0.0), (0.75 + 4.0 * f64::EPSILON, 4.0)];
    let initial = emitter(&fractions, 5.0).unwrap();
    let (_, report) = initial
        .trial_with_rates(0.2, &[300.0], &[0.0], budget(), || Ok(()))
        .unwrap();
    let counts = &report.spatial.burnup[0].reaction_counts;
    let decays = counts[1] + counts[2];
    close(
        report.transferred_decay_neutron_energy_j,
        2.0 * decays,
        1e-8,
    );
    close(
        report.released.radiation_j,
        5.0 * counts[0] + 4.0 * counts[2] - 2.0 * decays,
        1e-8,
    );
    assert!(report.released.radiation_j > 0.0);
}

#[test]
fn g3_original_fusion_source_audit_excludes_radioactive_births_and_restart_keeps_accounts() {
    let initial = emitter(&[(1.0, 2.0)], 0.0).unwrap();
    let amount = Amount::new(1.0 / AVOGADRO_PER_MOL);
    let emission = FusionNeutronEmissionV1::new(
        Time::new(0.2),
        [amount, Amount::new(0.0)],
        [
            Energy::new(2.0 * amount.value() * AVOGADRO_PER_MOL),
            Energy::new(0.0),
        ],
    )
    .unwrap();
    let source =
        FusionNeutronSource::project_interval_mean(Arc::clone(initial.model()), emission, 0)
            .unwrap();
    let (middle, first) = initial
        .trial_with_fusion(0.2, &[300.0], &[source], budget(), || Ok(()))
        .unwrap();
    close(
        first.external_injected_neutrons,
        emission.total_particle_count(),
        1e-12,
    );
    assert!(first.decay_neutrons > first.external_injected_neutrons);
    let (after, last) = middle
        .trial_with_rates(0.2, &[300.0], &[0.0], budget(), || Ok(()))
        .unwrap();
    assert_eq!(last.external_injected_neutrons, 0.0);
    assert!(last.decay_neutrons > 0.0);
    close(
        after.cell_released(0).unwrap().total_j(),
        first.released.total_j() + last.released.total_j(),
        1e-8,
    );
    close(
        after.cell_nuclear_released(0).unwrap().radiation_j,
        first.nuclear_released.radiation_j + last.nuclear_released.radiation_j,
        1e-8,
    );
}

#[test]
fn g4_energy_and_late_checkpoint_refusals_preserve_finite_emitter_and_neutrons() {
    let initial = emitter(&[(1.0, 2.0)], 0.0).unwrap();
    let mut polls = 0;
    initial
        .trial_with_rates(0.2, &[300.0], &[0.0], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut remaining = polls - 2;
    let result = initial.trial_with_rates(0.2, &[300.0], &[0.0], budget(), || {
        if remaining == 0 {
            Err(FissionError::Cancelled)
        } else {
            remaining -= 1;
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(FusionNeutronActivationError::Source(
            FissionError::Cancelled
        ))
    ));
    let mut strict = budget();
    strict.energy_tolerance_j = 0.0;
    assert!(
        initial
            .trial_with_rates(0.2, &[300.0], &[0.0], strict, || Ok(()))
            .is_err()
    );
    assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0]);
    assert_eq!(initial.cell_released(0).unwrap(), EnergyChannels::default());
}
