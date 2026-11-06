//! G1/G3 bridges from manufactured processed ENDF to the original owners.
//! Cross sections, Q partitions and product stability are numerical fixtures,
//! not evaluated isotope data. Product identities/stoichiometry are explicit.

use super::*;
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
use crate::fusion_neutrons::{FusionNeutronSource, NeutronKineticModel};
use crate::irradiation::{
    GroupFlux, IrradiationBudget, IrradiationModel, IrradiationState, TabulatedIrradiationData,
};
use crate::{EnergyChannels, FissionBudget, FissionModel};
use fs_qty::{Amount, Energy, Frequency, Time};
use fs_thermochem::decay::DecayLawV1;
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL;
use fs_thermochem::fusion::neutrons::FusionNeutronEmissionV1;
use fs_thermochem::fusion::{NuclearSpeciesV1, SPEED_OF_LIGHT_M_S};

const Q: EnergyChannels = EnergyChannels {
    deposited_j: 4e-13,
    radiation_j: 3e-13,
    escaped_j: 2e-13,
};

fn id(z: u16, a: u16) -> NuclideId {
    NuclideId::new(z, a, 0).unwrap()
}

fn stable(id: NuclideId) -> BurnupNuclide {
    BurnupNuclide {
        id,
        law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
    }
}

fn selected(parent: NuclideId, reaction: NeutronReaction, temperature_k: f64) -> EndfSelection {
    EndfSelection {
        material_number: 325,
        nuclide: parent,
        reaction,
        temperature_k,
    }
}

fn source(selection: EndfSelection, sigma_barns: f64) -> String {
    let text = fixture(
        selection,
        2,
        &[(0.0, sigma_barns), (2e7, sigma_barns)],
        &[(2, 2)],
    );
    // Reuse the bounded format fixture while authoring this isotope's synthetic
    // mass and nonfission flag explicitly, instead of retaining the U-235 header.
    let mat = i32::from(selection.material_number);
    let mt = reaction_mt(selection.reaction);
    let za = f64::from(selection.nuclide.atomic_number()) * 1000.0
        + f64::from(selection.nuclide.mass_number());
    let awr = f64::from(selection.nuclide.mass_number());
    text.replacen(
        &control(za, 233.0, [2, 1, 0, 7], (mat, 1, 451)),
        &control(
            za,
            awr,
            [
                2,
                i32::from(selection.reaction == NeutronReaction::Fission),
                0,
                7,
            ],
            (mat, 1, 451),
        ),
        1,
    )
    .replacen(
        &control(za, 233.0, [0; 4], (mat, 3, mt)),
        &control(za, awr, [0; 4], (mat, 3, mt)),
        1,
    )
}

fn table(
    parent: NuclideId,
    reaction: NeutronReaction,
    grid: Arc<NeutronGroups>,
) -> Arc<CrossSectionTable> {
    let curves: Vec<_> = [(300.0, 2.0), (600.0, 6.0)]
        .into_iter()
        .map(|(temperature, sigma)| {
            let selected = selected(parent, reaction, temperature);
            parse(&source(selected, sigma), selected).unwrap()
        })
        .collect();
    Arc::new(
        collapse_temperature_table(
            "manufactured-charged-absorption-temperature-table".into(),
            &curves,
            grid,
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn data(table: Arc<CrossSectionTable>, products: [NuclideId; 2]) -> BurnupData {
    BurnupData {
        data_id: "manufactured-explicit-absorption-products-and-q".into(),
        groups: Arc::clone(table.groups()),
        temperature_min_k: 300.0,
        temperature_max_k: 600.0,
        nuclides: vec![
            stable(table.nuclide()),
            stable(products[0]),
            stable(products[1]),
        ],
        reactions: vec![BurnupReaction {
            data_id: "manufactured-partial-absorption-channel".into(),
            parent: 0,
            rate: BurnupRate::TabulatedCapture {
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
            energy_per_event: Q,
        }],
    }
}

fn burnup_budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-18,
    }
}

#[test]
fn g1_all_charged_channels_preserve_mt_and_reach_original_irradiation_and_burnup() {
    let parent = id(13, 27);
    let grid = groups(&[0.0, 2e7]);
    let flux = GroupFlux::new(Arc::clone(&grid), vec![1e28]).unwrap();
    let expected = 100.0 * -(-1.0_f64).exp_m1(); // 4 b * 1e28/m²/s * 0.25 s.
    for (reaction, mt, products) in [
        (
            NeutronReaction::ProtonProduction,
            103,
            [id(1, 1), id(12, 27)],
        ),
        (
            NeutronReaction::DeuteronProduction,
            104,
            [id(1, 2), id(12, 26)],
        ),
        (
            NeutronReaction::TritonProduction,
            105,
            [id(1, 3), id(12, 25)],
        ),
        (
            NeutronReaction::Helium3Production,
            106,
            [id(2, 3), id(11, 25)],
        ),
        (
            NeutronReaction::AlphaProduction,
            107,
            [id(2, 4), id(11, 24)],
        ),
    ] {
        assert_eq!(reaction_mt(reaction), mt);
        assert!(reaction.is_nonfission_absorption());
        let capture = table(parent, reaction, Arc::clone(&grid));
        assert_eq!(capture.reaction(), reaction);
        assert_eq!(capture.nuclide(), parent);
        close(
            capture.cross_section_at_group(450.0, 0).unwrap(),
            4.0 * BARN_M2,
        );
        assert_eq!(products.iter().map(|p| p.atomic_number()).sum::<u16>(), 13);
        assert_eq!(products.iter().map(|p| p.mass_number()).sum::<u16>(), 28);

        let model = Arc::new(
            BurnupModel::new(
                data(Arc::clone(&capture), products),
                BurnupLimits::default(),
            )
            .unwrap(),
        );
        let original = BurnupInventory::new(model, vec![100.0, 0.0, 0.0]).unwrap();
        let (next, report) = original
            .trial(0.25, &flux, 450.0, burnup_budget(), || Ok(()))
            .unwrap();
        close(report.captures, expected);
        close(next.populations()[0], 100.0 - expected);
        close(next.populations()[1], expected);
        close(next.populations()[2], expected);
        close(report.generated_energy_j, expected * Q.total_j());
        assert_eq!(report.fissions, 0.0);
        assert_eq!(report.decays, 0.0);
        assert_eq!(original.populations(), &[100.0, 0.0, 0.0]);
        // Imported negative QI remains metadata; only explicitly authored energy
        // enters the physical owner. This fixture makes no endothermic claim.
        let selection = selected(parent, reaction, 300.0);
        assert!(
            parse(&source(selection, 2.0), selection)
                .unwrap()
                .reaction_q_j()
                < 0.0
        );

        let fission = Arc::new(
            CrossSectionTable::new(
                CrossSectionData {
                    data_id: "manufactured-zero-fission".into(),
                    nuclide: parent,
                    reaction: NeutronReaction::Fission,
                    groups: Arc::clone(&grid),
                    temperatures_k: vec![300.0, 600.0],
                    cross_sections_m2: vec![vec![0.0], vec![0.0]],
                },
                CrossSectionLimits::default(),
            )
            .unwrap(),
        );
        let irradiation = Arc::new(
            IrradiationModel::new_tabulated(TabulatedIrradiationData {
                data_id: "manufactured-typed-parent-removal".into(),
                fission,
                capture: Arc::clone(&capture),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                fission_model: Arc::new(
                    FissionModel::new("manufactured-inactive-fission".into(), Q, vec![]).unwrap(),
                ),
                capture_energy_per_event: Q,
            })
            .unwrap(),
        );
        let single = IrradiationState::new(irradiation, 100.0, vec![]).unwrap();
        let (single, single_report) = single
            .trial(
                0.25,
                &flux,
                450.0,
                IrradiationBudget {
                    max_neutron_groups: 2,
                    source: FissionBudget {
                        max_decay_groups: 0,
                        energy_tolerance_j: 1e-18,
                        fuel_tolerance_nuclei: 1e-8,
                    },
                },
                || Ok(()),
            )
            .unwrap();
        close(single.remaining_parent_nuclei(), next.populations()[0]);
        close(single_report.captures, expected);
        close(single_report.released.total_j(), report.released.total_j());
    }
    assert!(!NeutronReaction::Fission.is_nonfission_absorption());
}

#[test]
fn g0_absorption_identity_breakup_and_fission_confusion_refuse() {
    let parent = id(3, 6);
    let selection = selected(parent, NeutronReaction::TritonProduction, 300.0);
    let original = source(selection, 2.0);
    let breakup = original.replacen(
        &control(2e8, -1e6, [0, 0, 1, 2], (325, 3, 105)),
        &control(2e8, -1e6, [0, 22, 1, 2], (325, 3, 105)),
        1,
    );
    assert_ne!(original, breakup);
    assert_eq!(
        parse(&breakup, selection),
        Err(EndfError::Unsupported("MF3 breakup flag LR must be zero"))
    );
    assert_eq!(
        parse(
            &original,
            EndfSelection {
                reaction: NeutronReaction::Capture,
                ..selection
            }
        ),
        Err(EndfError::MissingSelection)
    );
    let other = selected(parent, NeutronReaction::AlphaProduction, 600.0);
    let curves = [
        parse(&original, selection).unwrap(),
        parse(&source(other, 6.0), other).unwrap(),
    ];
    let grid = groups(&[0.0, 2e7]);
    assert!(
        collapse_temperature_table(
            "mixed-reaction-refusal".into(),
            &curves,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(())
        )
        .is_err()
    );

    let capture = table(parent, NeutronReaction::TritonProduction, Arc::clone(&grid));
    let mut wrong = data(Arc::clone(&capture), [id(1, 3), id(2, 4)]);
    wrong.reactions[0].rate = BurnupRate::TabulatedFission {
        cross_sections: capture,
    };
    assert!(BurnupModel::new(wrong, BurnupLimits::default()).is_err());
    let fission = table(parent, NeutronReaction::Fission, grid);
    assert!(
        BurnupModel::new(data(fission, [id(1, 3), id(2, 4)]), BurnupLimits::default()).is_err()
    );
}

#[test]
fn g3_li6_channel_breeds_explicit_tritium_and_helium_on_the_same_fusion_history() {
    let grid = groups(&[0.0, 14e6, 20e6]);
    let capture = table(
        id(3, 6),
        NeutronReaction::TritonProduction,
        Arc::clone(&grid),
    );
    let material = Arc::new(
        BurnupModel::new(
            data(Arc::clone(&capture), [id(1, 3), id(2, 4)]),
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let energy = [13e6 * EV_J, 15e6 * EV_J];
    let rest = NuclearSpeciesV1::Neutron.mass_kg() * SPEED_OF_LIGHT_M_S.powi(2);
    let speeds = energy.map(|e| {
        let x = e / rest;
        SPEED_OF_LIGHT_M_S * (x * (2.0 + x)).sqrt() / (1.0 + x)
    });
    let volume = 1e-24;
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-li6-absorption-cells".into(),
                groups: Arc::clone(&grid),
                speeds_m_s: speeds.to_vec(),
                cell_volumes_m3: vec![volume],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|group| DiffusionReaction {
                        data_id: format!("manufactured-li6-MT105-group-{group}"),
                        cell: 0,
                        incident_group: group,
                        macroscopic_cross_section_m_inv: 100.0 / volume * 2.0 * BARN_M2,
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let kinetic = Arc::new(
        NeutronKineticModel::new(
            "manufactured-li6-incident-heat".into(),
            Arc::clone(&diffusion),
            energy.to_vec(),
            vec![1.0; 2],
        )
        .unwrap(),
    );
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(diffusion, vec![0.0; 2], vec![]).unwrap(),
        vec![BurnupInventory::new(material, vec![100.0, 0.0, 0.0]).unwrap()],
        (0..2)
            .map(|diffusion_reaction| ReactionBinding {
                diffusion_reaction,
                burnup_reaction: 0,
            })
            .collect(),
    )
    .unwrap();
    let mut state = FusionNeutronActivationState::new(Arc::clone(&kinetic), spatial).unwrap();
    let allowance = FusionNeutronActivationBudget {
        spatial: DiffusionBurnupBudget {
            diffusion: DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-8,
            },
            burnup: burnup_budget(),
            max_iterations: 64,
            reaction_tolerance: 1e-10,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-18,
    };
    let duration = 2e-8;
    let amount = Amount::new(4.0 / AVOGADRO_PER_MOL);
    let actual_count = amount.value() * AVOGADRO_PER_MOL;
    let emission = FusionNeutronEmissionV1::new(
        Time::new(duration),
        [amount, Amount::new(0.0)],
        [Energy::new(actual_count * 14e6 * EV_J), Energy::new(0.0)],
    )
    .unwrap();
    let source = FusionNeutronSource::project_interval_mean(kinetic, emission, 0).unwrap();
    let mut captures = 0.0;
    let mut at_shutoff = 0.0;
    for step in 0..15 {
        let input = if step < 10 {
            std::slice::from_ref(&source)
        } else {
            &[]
        };
        let (next, report) = state
            .trial_with_fusion(duration, &[450.0], input, allowance, || Ok(()))
            .unwrap();
        assert_eq!(report.kinetic.diffusion, report.spatial.neutrons);
        assert!(report.energy_roundoff_j <= allowance.energy_tolerance_j);
        captures += report.spatial.burnup[0].captures;
        state = next;
        if step == 9 {
            at_shutoff = captures;
        }
    }
    assert!(captures > at_shutoff && captures < 100.0);
    let populations = state.inventories()[0].populations();
    let near = |actual: f64, expected: f64| {
        assert!(
            (actual - expected).abs() < 1e-7,
            "{actual:e} != {expected:e}"
        )
    };
    near(populations[0], 100.0 - captures);
    near(populations[1], captures); // Authored tritium product.
    near(populations[2], captures); // Authored helium-4 product.
    let neutron_count: f64 = state.neutrons().neutron_counts().iter().sum();
    near(neutron_count + captures, 10.0 * actual_count);
    near(
        3.0 * populations[0] + populations[1] + 2.0 * populations[2],
        300.0,
    );
    near(
        6.0 * populations[0] + 3.0 * populations[1] + 4.0 * populations[2] + neutron_count,
        600.0 + 10.0 * actual_count,
    );
    let released = state.cell_released(0).unwrap();
    let nuclear = state.cell_nuclear_released(0).unwrap();
    let incident = state.cell_kinetic_released(0).unwrap();
    assert_eq!(incident.radiation_j, 0.0);
    assert!((nuclear.total_j() - captures * Q.total_j()).abs() < 1e-18);
    assert!((released.deposited_j - nuclear.deposited_j - incident.deposited_j).abs() < 1e-18);
    assert!(
        (state.kinetic_energy_j().unwrap() + released.total_j()
            - 10.0 * emission.total_kinetic_energy().value()
            - captures * Q.total_j())
        .abs()
            < 1e-18
    );
    let BurnupRate::TabulatedCapture { cross_sections } =
        &state.inventories()[0].model().reactions()[0].rate
    else {
        panic!("original typed table lost")
    };
    assert!(Arc::ptr_eq(cross_sections, &capture));
    assert_eq!(cross_sections.reaction(), NeutronReaction::TritonProduction);
    assert!(
        state.neutrons().model().reactions()[0].macroscopic_cross_section_m_inv
            < 100.0 / volume * 4.0 * BARN_M2
    );
}
