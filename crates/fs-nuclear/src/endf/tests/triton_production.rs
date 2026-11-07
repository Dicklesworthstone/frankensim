//! G1/G3 exclusive neutron/triton channels from processed MF3 to finite targets.
//! Coefficients, speeds and Q are manufactured. Product decay is disabled in
//! this short-interval fixture; this is not an evaluated Li-7 breeding prediction.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState, ReactionBinding};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind, DiffusionState,
};
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationState,
};
use crate::fusion_neutrons::{NeutronKineticModel, NeutronReactionEnergy};
use fs_thermochem::decay::DecayLawV1;

const MEV: f64 = 1e6 * EV_J;
const VOLUME: f64 = 2e-26;

fn selected(kind: NeutronReaction, temperature: f64) -> EndfSelection {
    EndfSelection {
        material_number: 328,
        nuclide: NuclideId::new(3, 7, 0).unwrap(),
        reaction: kind,
        temperature_k: temperature,
    }
}

fn source(id: EndfSelection, lr: i32) -> String {
    let mat = i32::from(id.material_number);
    let section = (mat, 3, reaction_mt(id.reaction));
    fixture(
        id,
        2,
        &[(0.0, 0.0), (8e6, 0.0), (8e6, 2.0), (20e6, 2.0)],
        &[(4, 1)],
    )
    .replacen(
        &control(3007.0, 233.0, [2, 1, 0, 7], (mat, 1, 451)),
        &control(3007.0, 7.0, [2, 0, 0, 7], (mat, 1, 451)),
        1,
    )
    .replacen(
        &control(3007.0, 233.0, [0; 4], section),
        &control(3007.0, 7.0, [0; 4], section),
        1,
    )
    .replacen(
        &control(2e8, -1e6, [0, 0, 1, 4], section),
        &control(
            -2.4e6,
            if lr == 0 { -2.4e6 } else { -4e6 },
            [0, lr, 1, 4],
            section,
        ),
        1,
    )
}

fn table(curve: &PointwiseCrossSection) -> Arc<CrossSectionTable> {
    Arc::new(
        collapse_temperature_table(
            "manufactured-lithium-neutron-triton-table".into(),
            std::slice::from_ref(curve),
            groups(&[0.0, 8e6, 20e6]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn nuclides(table: Arc<CrossSectionTable>) -> BurnupData {
    BurnupData {
        data_id: "manufactured-explicit-Li7-tritium-helium".into(),
        groups: Arc::clone(table.groups()),
        temperature_min_k: 300.0,
        temperature_max_k: 300.0,
        nuclides: [
            table.nuclide(),
            NuclideId::new(1, 3, 0).unwrap(),
            NuclideId::new(2, 4, 0).unwrap(),
        ]
        .into_iter()
        .map(|id| BurnupNuclide {
            id,
            law: DecayLawV1::stable(),
        })
        .collect(),
        reactions: vec![BurnupReaction {
            data_id: "authored-one-triton-one-alpha".into(),
            parent: 0,
            rate: BurnupRate::TabulatedNeutronProduction {
                cross_sections: table,
            },
            products: vec![
                BurnupProduct {
                    nuclide: 1,
                    yield_per_event: 1.0,
                },
                BurnupProduct {
                    nuclide: 2,
                    yield_per_event: 1.0,
                },
            ],
            // All prompt energy belongs to the explicit neutron energy model.
            energy_per_event: EnergyChannels::default(),
        }],
    }
}

fn finite(table: Arc<CrossSectionTable>) -> FusionNeutronActivationState {
    let nuclear =
        Arc::new(BurnupModel::new(nuclides(Arc::clone(&table)), BurnupLimits::default()).unwrap());
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-Li7-neutron-preserving-breakup".into(),
                groups: Arc::clone(table.groups()),
                speeds_m_s: vec![1.0; 2],
                cell_volumes_m3: vec![VOLUME],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|group| DiffusionReaction {
                        data_id: format!("manufactured-Li7-group-{group}"),
                        cell: 0,
                        incident_group: group,
                        macroscopic_cross_section_m_inv: 100.0
                            * table.cross_section_at_group(300.0, group).unwrap()
                            / VOLUME,
                        kind: DiffusionReactionKind::NeutronProduction {
                            neutron_yields: vec![1.0, 0.0],
                        },
                    })
                    .collect(),
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let model = Arc::new(
        NeutronKineticModel::new_with_reaction_energies(
            "authored-Li7-prompt-energy-and-threshold".into(),
            Arc::clone(&diffusion),
            vec![2.0 * MEV, 14.0 * MEV],
            vec![0.8; 2],
            // Intentionally distinct from both source QM and QI. Source metadata
            // cannot silently become a heat, outgoing spectrum or threshold model.
            vec![
                NeutronReactionEnergy {
                    nuclear_q_j: -2.5 * MEV,
                    threshold_j: 3.0 * MEV
                };
                2
            ],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0, 20.0], vec![]).unwrap(),
        vec![BurnupInventory::new(nuclear, vec![100.0, 0.0, 0.0]).unwrap()],
        (0..2)
            .map(|group| ReactionBinding {
                diffusion_reaction: group,
                burnup_reaction: 0,
            })
            .collect(),
    )
    .unwrap();
    FusionNeutronActivationState::new(model, spatial).unwrap()
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
                energy_tolerance_j: 1e-19,
            },
            max_iterations: 64,
            reaction_tolerance: 1e-11,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-18,
    }
}

fn near(actual: f64, expected: f64) {
    let allowance = 2e-9 * expected.abs().max(1e-20);
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_lithium_neutron_triton_identities_drive_finite_products_without_consuming_neutrons() {
    for (kind, mt, lr) in [
        (NeutronReaction::NeutronTritonProduction, 33, 0),
        (NeutronReaction::ContinuumNeutronTritonBreakup, 91, 33),
    ] {
        let id = selected(kind, 300.0);
        let curve = parse(&source(id, lr), id).unwrap();
        assert_eq!(reaction_mt(curve.reaction()), mt);
        assert_eq!(curve.breakup_flag(), lr);
        assert_eq!(curve.reaction().emitted_neutrons_per_event(), Some(1));
        assert!(!curve.reaction().is_nonfission_absorption());
        close(curve.mass_difference_q_j(), -2.4 * MEV);
        close(
            curve.reaction_q_j(),
            if lr == 0 { -2.4 * MEV } else { -4.0 * MEV },
        );
        let table = table(&curve);
        assert_eq!(table.reaction(), kind);
        assert_eq!(table.cross_section_at_group(300.0, 0).unwrap(), 0.0);
        close(table.cross_section_at_group(300.0, 1).unwrap(), 2e-28);
        let initial = finite(table);
        let duration = 0.7_f64;
        let attenuation = (-80.0 * 0.01 * duration).exp();
        let expected_fast = 80.0 * 20.0 * attenuation / (80.0 + 20.0 * (1.0 - attenuation));
        let expected_events = 20.0 - expected_fast;
        let mut errors = Vec::new();
        for intervals in [8, 16, 32] {
            let mut current = initial.clone();
            let mut events = 0.0;
            for _ in 0..intervals {
                let (next, report) = current
                    .trial_with_rates(
                        duration / f64::from(intervals),
                        &[300.0],
                        &[0.0; 2],
                        budget(),
                        || Ok(()),
                    )
                    .unwrap();
                let count = report.spatial.burnup[0].neutron_producing_reactions;
                near(report.spatial.neutrons.neutron_producing_reactions, count);
                near(report.spatial.neutrons.secondary_neutrons, count);
                near(report.nuclear_generated_energy_j, -2.5 * MEV * count);
                near(report.released.deposited_j, 0.8 * 9.5 * MEV * count);
                near(report.released.radiation_j, 0.2 * 9.5 * MEV * count);
                assert_eq!(report.nuclear_released, EnergyChannels::default());
                assert_eq!(
                    [
                        report.spatial.neutrons.fissions,
                        report.spatial.neutrons.captures
                    ],
                    [0.0; 2]
                );
                assert!(report.energy_roundoff_j <= budget().energy_tolerance_j);
                events += count;
                current = next;
            }
            let neutron = current.neutrons().neutron_counts();
            let nuclei = current.inventories()[0].populations();
            near(neutron[0] + neutron[1], 20.0);
            near(neutron[0], events);
            near(100.0 - nuclei[0], events);
            near(nuclei[1], events);
            near(nuclei[2], events);
            near(
                7.0 * nuclei[0] + 3.0 * nuclei[1] + 4.0 * nuclei[2] + neutron.iter().sum::<f64>(),
                720.0,
            );
            near(3.0 * nuclei[0] + nuclei[1] + 2.0 * nuclei[2], 300.0);
            let q = current.cell_prompt_nuclear_energy_j(0).unwrap();
            near(q, -2.5 * MEV * events);
            near(
                current.kinetic_energy_j().unwrap() + current.cell_released(0).unwrap().total_j(),
                initial.kinetic_energy_j().unwrap() + q,
            );
            errors.push((events - expected_events).abs());
        }
        assert!(
            errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 5e-4, "{errors:?}");
        assert_eq!(initial.neutrons().neutron_counts(), &[0.0, 20.0]);
        assert_eq!(initial.inventories()[0].populations(), &[100.0, 0.0, 0.0]);
    }
}

#[test]
fn g0_lithium_channels_refuse_ambiguous_breakup_absorption_and_wrong_neutron_multiplicity() {
    for (kind, bad_lr) in [
        (NeutronReaction::NeutronTritonProduction, 33),
        (NeutronReaction::NeutronTritonProduction, 1),
        (NeutronReaction::ContinuumNeutronTritonBreakup, 0),
        (NeutronReaction::ContinuumNeutronTritonBreakup, 1),
        (NeutronReaction::ContinuumNeutronTritonBreakup, 22),
    ] {
        let id = selected(kind, 300.0);
        assert!(matches!(
            parse(&source(id, bad_lr), id),
            Err(EndfError::Unsupported(_))
        ));
    }
    let mut curves = Vec::new();
    for (kind, lr, temperature) in [
        (NeutronReaction::NeutronTritonProduction, 0, 300.0),
        (NeutronReaction::ContinuumNeutronTritonBreakup, 33, 600.0),
    ] {
        let id = selected(kind, temperature);
        let curve = parse(&source(id, lr), id).unwrap();
        curves.push(curve);
    }
    assert!(
        collapse_temperature_table(
            "no-automatic-representation-merge".into(),
            &curves,
            groups(&[0.0, 8e6, 20e6]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(()),
        )
        .is_err()
    );
    for curve in curves {
        let table = table(&curve);
        let mut wrong = nuclides(Arc::clone(&table));
        wrong.temperature_min_k = curve.metadata().temperature_k;
        wrong.temperature_max_k = curve.metadata().temperature_k;
        wrong.reactions[0].rate = BurnupRate::TabulatedCapture {
            cross_sections: Arc::clone(&table),
        };
        assert!(BurnupModel::new(wrong, BurnupLimits::default()).is_err());
    }
    let id = selected(NeutronReaction::ContinuumNeutronTritonBreakup, 300.0);
    let initial = finite(table(&parse(&source(id, 33), id).unwrap()));
    for kind in [
        DiffusionReactionKind::Capture,
        DiffusionReactionKind::NeutronProduction {
            neutron_yields: vec![2.0, 0.0],
        },
    ] {
        let mut data = initial.model().diffusion().data().clone();
        data.reactions[1].kind = kind;
        let diffusion = Arc::new(DiffusionModel::new(data, DiffusionLimits::default()).unwrap());
        assert!(
            DiffusionBurnupState::new(
                DiffusionState::new(diffusion, vec![0.0, 20.0], vec![]).unwrap(),
                initial.inventories().to_vec(),
                initial.spatial().bindings().to_vec(),
            )
            .is_err()
        );
    }
}
