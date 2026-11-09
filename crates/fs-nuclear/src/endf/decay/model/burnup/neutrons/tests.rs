//! Manufactured ENDF sources and independent finite-history solutions. The
//! source values exercise dimensional accounting, not isotope calibration.

#[path = "delayed_tests.rs"]
pub(crate) mod delayed;

use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupRate, BurnupReaction,
};
use crate::diffusion::burnup::{DiffusionBurnupBudget, DiffusionBurnupState};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind, DiffusionState,
};
use crate::endf::decay::model::burnup::SpontaneousFissionProductBinding;
use crate::endf::decay::model::{DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry};
use crate::endf::decay::{EndfDecayReadLimits, EndfDecaySelection};
use crate::endf::fission_yields::spontaneous::SpontaneousFissionYieldEvaluation;
use crate::endf::fission_yields::{EndfFissionYieldReadLimits, EndfFissionYieldSelection};
use crate::endf::spontaneous_neutrons::{
    EndfSpontaneousNeutronSelection, SpontaneousFissionSpectrumLaw, SpontaneousNeutronReadLimits,
    tests::{prompt_analytic_fixture, prompt_fixture, prompt_general_evaporation_fixture},
};
use crate::fusion_neutrons::activation::{
    FusionNeutronActivationBudget, FusionNeutronActivationState,
};
use crate::irradiation::NeutronGroups;

const EV_J: f64 = 1.602_176_634e-19;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, allowance {tolerance:.3e}"
    );
}

fn id(z: u16, a: u16) -> NuclideId {
    NuclideId::new(z, a, 0).unwrap()
}
fn parent() -> NuclideId {
    id(80, 202)
}
fn za(id: NuclideId) -> f64 {
    1000.0 * f64::from(id.atomic_number()) + f64::from(id.mass_number())
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let text = format!("{mantissa}{exponent:+}");
        if text.len() <= 11 {
            return format!("{text:>11}");
        }
    }
    panic!("fixture number too wide");
}

fn row(values: [f64; 6], mf: u16, mt: u16) -> String {
    format!("{}{:4}{mf:2}{mt:3}{:5}", values.map(real).concat(), 777, 1)
}

fn control(first: f64, second: f64, integers: [i32; 4], mf: u16, mt: u16) -> String {
    let fields = [
        real(first),
        real(second),
        format!("{:11}", integers[0]),
        format!("{:11}", integers[1]),
        format!("{:11}", integers[2]),
        format!("{:11}", integers[3]),
    ];
    format!("{}{:4}{mf:2}{mt:3}{:5}", fields.concat(), 777, 1)
}

fn evaluation(
    nuclide: NuclideId,
    half_life: Option<f64>,
    energies_ev: [f64; 3],
    modes: &[(f64, f64)],
    neutron_line: Option<(f64, f64)>,
) -> Arc<DecayEvaluation> {
    let mut lines = vec![
        control(
            za(nuclide),
            f64::from(nuclide.mass_number()),
            [
                0,
                0,
                i32::from(half_life.is_none()),
                i32::from(neutron_line.is_some()),
            ],
            8,
            457,
        ),
        control(half_life.unwrap_or(0.0), 0.0, [0, 0, 6, 0], 8, 457),
        row(
            [
                energies_ev[0],
                0.0,
                energies_ev[1],
                0.0,
                energies_ev[2],
                0.0,
            ],
            8,
            457,
        ),
        control(
            0.0,
            1.0,
            [0, 0, 6 * modes.len().max(1) as i32, modes.len() as i32],
            8,
            457,
        ),
    ];
    if modes.is_empty() {
        lines.push(row([0.0; 6], 8, 457));
    } else {
        for &(rtyp, probability) in modes {
            lines.push(row([rtyp, 0.0, 1000.0, 0.0, probability, 0.0], 8, 457));
        }
    }
    if let Some((rtyp, intensity)) = neutron_line {
        // A zero intensity supplies the actual SF average-only STYP=5 shape:
        // no discrete transitions, FD=FC=0, reported average energy retained.
        let transitions = i32::from(intensity > 0.0);
        lines.push(control(0.0, 5.0, [0, 0, 6, transitions], 8, 457));
        lines.push(row(
            [f64::from(transitions), 0.0, 2.0, 0.0, 0.0, 0.0],
            8,
            457,
        ));
        if transitions != 0 {
            lines.push(control(2.0, 0.0, [0, 0, 6, 0], 8, 457));
            lines.push(row([rtyp, 0.0, intensity, 0.0, 0.0, 0.0], 8, 457));
        }
    }
    lines.push(control(0.0, 0.0, [0; 4], 8, 0));
    Arc::new(
        DecayEvaluation::parse(
            format!("manufactured-MF8-{nuclide:?}"),
            EndfDecaySelection {
                material_number: 777,
                nuclide,
            },
            &lines.join("\n"),
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn fission_products() -> Arc<SpontaneousFissionYieldEvaluation> {
    let mut lines = vec![
        control(za(parent()), 202.0, [-1, 0, 1, 3], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 6], 1, 451),
        control(0.0, 0.0, [0, 0, 5, 8], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 0),
        control(za(parent()), 202.0, [1, 0, 0, 0], 8, 454),
        control(0.0, 0.0, [0, 0, 8, 2], 8, 454),
    ];
    for chunk in [
        za(id(39, 99)),
        0.0,
        1.0,
        0.0,
        za(id(41, 101)),
        0.0,
        1.0,
        0.0,
    ]
    .chunks(6)
    {
        lines.push(row(
            std::array::from_fn(|i| chunk.get(i).copied().unwrap_or(0.0)),
            8,
            454,
        ));
    }
    lines.push(control(0.0, 0.0, [0; 4], 8, 0));
    Arc::new(
        SpontaneousFissionYieldEvaluation::parse(
            "manufactured-independent-fragments".into(),
            EndfFissionYieldSelection {
                material_number: 777,
                nuclide: parent(),
            },
            &lines.join("\n"),
            EndfFissionYieldReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn routing(heavy_radiation: f64) -> DecayEnergyRouting {
    DecayEnergyRouting {
        light_particles: DecayEnergyFractions {
            deposited: 1.0,
            radiation: 0.0,
            escaped: 0.0,
        },
        electromagnetic: DecayEnergyFractions {
            deposited: 0.0,
            radiation: 1.0,
            escaped: 0.0,
        },
        heavy_particles: DecayEnergyFractions {
            deposited: 1.0 - heavy_radiation,
            radiation: heavy_radiation,
            escaped: 0.0,
        },
        include_reported_neutrinos: false,
    }
}

fn assembly(
    parent_modes: &[(f64, f64)],
    parent_line: Option<(f64, f64)>,
    fragment_spectrum: bool,
    heavy_radiation: f64,
) -> Arc<EvaluatedBurnupDecayAssembly> {
    let stable = |nuclide| DecayModelEntry {
        evaluation: evaluation(nuclide, None, [0.0; 3], &[], None),
        routing: routing(1.0),
    };
    let mut entries = vec![
        DecayModelEntry {
            evaluation: evaluation(
                parent(),
                Some(2.0),
                [1.0, 2.0, 6.0],
                parent_modes,
                parent_line,
            ),
            routing: routing(heavy_radiation),
        },
        stable(id(78, 198)),
        DecayModelEntry {
            evaluation: evaluation(
                id(39, 99),
                Some(1.0),
                [1.0, 0.0, 3.0],
                &[(1.5, 1.0)],
                fragment_spectrum.then_some((1.5, 1.0)),
            ),
            routing: routing(1.0),
        },
        stable(id(40, 98)),
        stable(id(41, 101)),
    ];
    if parent_modes
        .iter()
        .any(|&(rtyp, probability)| rtyp == 5.0 && probability > 0.0)
    {
        entries.push(stable(id(80, 201)));
    }
    Arc::new(
        EvaluatedBurnupDecayAssembly::new(
            "manufactured-SF-and-delayed-fragment".into(),
            entries,
            vec![SpontaneousFissionProductBinding {
                parent: parent(),
                yields: fission_products(),
                exported_products: vec![],
            }],
            Arc::new(
                NeutronGroups::new(
                    "manufactured-SF-groups".into(),
                    vec![0.0, 2.0 * EV_J, 4.0 * EV_J],
                )
                .unwrap(),
            ),
            [300.0, 600.0],
            BurnupLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn kinetic(groups: Arc<NeutronGroups>) -> Arc<NeutronKineticModel> {
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-fixed-neutron-absorber".into(),
                groups,
                speeds_m_s: vec![1.0; 2],
                cell_volumes_m3: vec![0.01],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|incident_group| DiffusionReaction {
                        data_id: format!("fixed-background-capture-{incident_group}"),
                        cell: 0,
                        incident_group,
                        macroscopic_cross_section_m_inv: 0.3,
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
    Arc::new(
        NeutronKineticModel::new(
            "manufactured-SF-kinetic-model".into(),
            diffusion,
            vec![EV_J, 3.0 * EV_J],
            vec![0.25; 2],
        )
        .unwrap(),
    )
}

fn prompt(
    selection: EndfSpontaneousNeutronSelection,
    probability: f64,
    density: f64,
) -> SpontaneousFissionPromptNeutronBinding {
    let text = prompt_fixture(
        selection,
        2.0,
        &[(probability, vec![(1.0, 0.0), (2.0, density), (3.0, 0.0)])],
    );
    SpontaneousFissionPromptNeutronBinding {
        parent: selection.nuclide,
        multiplicity: Arc::new(
            SpontaneousFissionNeutronYield::parse(
                "manufactured-two-prompt-neutrons".into(),
                selection,
                FissionNeutronKind::Prompt,
                &text,
                SpontaneousNeutronReadLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        ),
        spectrum: Arc::new(
            SpontaneousFissionNeutronSpectrum::parse(
                "manufactured-triangular-prompt-spectrum".into(),
                selection,
                FissionNeutronKind::Prompt,
                &text,
                SpontaneousNeutronReadLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        ),
    }
}

fn selection() -> EndfSpontaneousNeutronSelection {
    EndfSpontaneousNeutronSelection {
        material_number: 777,
        nuclide: parent(),
    }
}

fn general_evaporation(
    temperature_ev: f64,
    restriction_ev: f64,
    shape: &[(f64, f64)],
    law: i32,
) -> SpontaneousFissionPromptNeutronBinding {
    let text = prompt_general_evaporation_fixture(
        selection(),
        2.0,
        1.0,
        temperature_ev,
        restriction_ev,
        shape,
        law,
    );
    let mut source = prompt(selection(), 1.0, 1.0);
    source.spectrum = Arc::new(
        SpontaneousFissionNeutronSpectrum::parse(
            "manufactured-scaled-evaporation-spectrum".into(),
            selection(),
            FissionNeutronKind::Prompt,
            &text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    source
}

fn source_budget() -> EvaluatedBurnupNeutronBudget {
    EvaluatedBurnupNeutronBudget {
        limits: DecayEmissionLimits::default(),
        normalization_tolerance: 1e-12,
    }
}

fn analytic_source(
    law: i32,
    cutoff_ev: f64,
    probability: f64,
) -> SpontaneousFissionPromptNeutronBinding {
    let text = prompt_analytic_fixture(selection(), 2.0, probability, 0.0, 2.0, -cutoff_ev, law);
    let mut source = prompt(selection(), 1.0, 1.0);
    source.spectrum = Arc::new(
        SpontaneousFissionNeutronSpectrum::parse(
            format!("manufactured-finite-analytic-LF{law}"),
            selection(),
            FissionNeutronKind::Prompt,
            &text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    source
}

fn prepare_source(
    assembly: Arc<EvaluatedBurnupDecayAssembly>,
) -> Result<EvaluatedBurnupNeutrons, EvaluatedBurnupNeutronError> {
    let model = kinetic(Arc::clone(assembly.model().groups()));
    EvaluatedBurnupNeutrons::new(
        "manufactured-evaluated-neutron-source".into(),
        assembly,
        model,
        vec![prompt(selection(), 1.0, 1.0)],
        &[],
        source_budget(),
        || Ok(()),
    )
}

fn activation(
    source: &Arc<EvaluatedBurnupNeutrons>,
    parent_count: f64,
) -> FusionNeutronActivationState {
    let bound = source
        .bind_burnup(
            "original-finite-stock-neutrons".into(),
            Arc::clone(source.assembly().model()),
            || Ok(()),
        )
        .unwrap();
    assert!(Arc::ptr_eq(bound.source(), source));
    let mut populations = vec![0.0; source.assembly().nuclides().len()];
    populations[0] = parent_count;
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(
            Arc::clone(source.kinetic_model().diffusion()),
            vec![0.0; 2],
            vec![],
        )
        .unwrap(),
        vec![BurnupInventory::new(Arc::clone(source.assembly().model()), populations).unwrap()],
        vec![],
    )
    .unwrap()
    .with_decay_neutrons(vec![Arc::clone(bound.emission())])
    .unwrap();
    FusionNeutronActivationState::new(Arc::clone(source.kinetic_model()), spatial).unwrap()
}

/// Shared real parser → native finite stock fixture for material calorimetry.
pub(crate) fn material_fixture(
    parent_count: f64,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    material_fixture_with_source(parent_count, prompt(selection(), 1.0, 1.0))
}

pub(crate) fn material_fixture_general_evaporation(
    parent_count: f64,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    material_fixture_with_source(
        parent_count,
        general_evaporation(2.0, -3.0, &[(0.5, 0.0), (1.0, 2.0), (1.5, 0.0)], 2),
    )
}

pub(crate) fn material_fixture_maxwellian(
    parent_count: f64,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    material_fixture_with_source(parent_count, analytic_source(7, 3.0, 1.0))
}

pub(crate) fn material_fixture_evaporation(
    parent_count: f64,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    material_fixture_with_source(parent_count, analytic_source(9, 3.0, 1.0))
}

fn material_fixture_with_source(
    parent_count: f64,
    binding: SpontaneousFissionPromptNeutronBinding,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    let assembly = assembly(&[(6.0, 0.25), (4.0, 0.75)], Some((6.0, 0.0)), true, 0.5);
    let source = Arc::new(
        EvaluatedBurnupNeutrons::new(
            "manufactured-material-neutron-spectrum".into(),
            Arc::clone(&assembly),
            kinetic(Arc::clone(assembly.model().groups())),
            vec![binding],
            &[],
            source_budget(),
            || Ok(()),
        )
        .unwrap(),
    );
    let state = activation(&source, parent_count);
    (source, state)
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
                energy_tolerance_j: 1e-24,
            },
            max_iterations: 64,
            reaction_tolerance: 1e-11,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-24,
    }
}

#[test]
fn g1_evaluated_sf_and_fragment_neutrons_share_actual_births_stock_and_energy() {
    let (source, initial) = material_fixture(100.0);
    let mut polls = 0;
    let (after, report) = initial
        .trial_with_rates(2.0, &[300.0], &[0.0; 2], budget(), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for (&actual, expected) in after.inventories()[0]
        .populations()
        .iter()
        .zip([50.0, 37.5, 6.25, 6.25, 12.5])
    {
        close(actual, expected, 1e-9);
    }
    let lambda = core::f64::consts::LN_2 / 2.0;
    let captures_rate = 0.3;
    let retained = lambda
        * 100.0
        * ((0.5 - (-captures_rate * 2.0_f64).exp()) / (captures_rate - lambda)
            - 0.5 * (0.25 - (-captures_rate * 2.0_f64).exp()) / (captures_rate - 2.0 * lambda));
    let captures = 31.25 - retained;
    for &count in after.neutrons().neutron_counts() {
        close(count, retained * 0.5, 1e-9);
    }
    close(report.decay_neutrons, 31.25, 1e-9);
    close(report.spatial.neutrons.captures, captures, 1e-9);
    close(
        report.transferred_decay_neutron_energy_j,
        62.5 * EV_J,
        1e-27,
    );
    close(
        report.released.deposited_j,
        (206.25 + 0.5 * captures) * EV_J,
        1e-27,
    );
    close(
        report.released.radiation_j,
        (206.25 + 1.5 * captures) * EV_J,
        1e-27,
    );
    close(
        after.kinetic_energy_j().unwrap() + report.released.total_j(),
        475.0 * EV_J,
        1e-27,
    );
    assert_eq!(report.external_injected_neutrons, 0.0);
    let events = &report.spatial.burnup[0].reaction_counts;
    close(report.decay_neutrons, 2.0 * events[0] + events[2], 1e-9);
    assert_eq!(report.kinetic.diffusion, report.spatial.neutrons);
    assert!(Arc::ptr_eq(
        after.inventories()[0].model(),
        source.assembly().model()
    ));
    close(source.projection()[0].particles_per_decay, 0.5, 2e-15);
    close(
        source.projection()[0].kinetic_energy_j_per_decay,
        EV_J,
        1e-33,
    );
    close(source.projection()[2].particles_per_decay, 1.0, 2e-15);
    let normalization = &source.normalization()[0];
    assert_eq!(normalization.spontaneous_branch_sum, 0.25);
    assert_eq!(normalization.total_decay_branch_sum, 1.0);
    assert_eq!(normalization.fission_events_per_decay, 0.25);
    assert_eq!(
        source.coverage(),
        &[
            DecayEmissionCoverage::Reported,
            DecayEmissionCoverage::Stable,
            DecayEmissionCoverage::Reported,
            DecayEmissionCoverage::Stable,
            DecayEmissionCoverage::Stable
        ]
    );

    let (middle, _) = initial
        .trial_with_rates(0.7, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let (resumed, _) = middle
        .trial_with_rates(1.3, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    for (&actual, &expected) in resumed.inventories()[0]
        .populations()
        .iter()
        .zip(after.inventories()[0].populations())
    {
        close(actual, expected, 1e-9);
    }
    for (&actual, &expected) in resumed
        .neutrons()
        .neutron_counts()
        .iter()
        .zip(after.neutrons().neutron_counts())
    {
        close(actual, expected, 1e-9);
    }
    close(
        resumed.cell_released(0).unwrap().total_j(),
        report.released.total_j(),
        1e-27,
    );
    let mut remaining = polls - 1;
    assert!(
        initial
            .trial_with_rates(2.0, &[300.0], &[0.0; 2], budget(), || {
                if remaining == 0 {
                    Err(FissionError::Cancelled)
                } else {
                    remaining -= 1;
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(
        initial.inventories()[0].populations(),
        &[100.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(initial.neutrons().neutron_counts(), &[0.0; 2]);
    assert_eq!(initial.cell_released(0), Some(EnergyChannels::default()));
}

#[test]
fn g0_source_normalization_is_explicit_and_admitted_discrepancies_are_not_repaired() {
    let original = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.5);
    let build = |probability, density, tolerance| {
        EvaluatedBurnupNeutrons::new(
            "declared-normalization-allowance".into(),
            Arc::clone(&original),
            kinetic(Arc::clone(original.model().groups())),
            vec![prompt(selection(), probability, density)],
            &[],
            EvaluatedBurnupNeutronBudget {
                normalization_tolerance: tolerance,
                ..source_budget()
            },
            || Ok(()),
        )
    };
    assert!(matches!(
        build(1.0, 1.001, 1e-12),
        Err(EvaluatedBurnupNeutronError::Normalization {
            component: Some(0),
            ..
        })
    ));
    assert!(matches!(
        build(0.999, 1.0, 1e-12),
        Err(EvaluatedBurnupNeutronError::Normalization {
            component: None,
            ..
        })
    ));
    let admitted = build(1.0, 1.001, 0.002).unwrap();
    let diagnostic = &admitted.normalization()[0];
    assert_eq!(diagnostic.declared_particles_per_fission, 2.0);
    close(diagnostic.integrated_particles_per_fission, 2.002, 2e-15);
    close(diagnostic.multiplicity_residual, 0.002, 2e-15);
    close(admitted.projection()[0].particles_per_decay, 0.5005, 1e-15);
    close(
        admitted.projection()[0].kinetic_energy_j_per_decay,
        1.001 * EV_J,
        1e-33,
    );
    let SpontaneousFissionSpectrumLaw::Tabulated { points, .. } =
        &admitted.spontaneous_sources()[0].spectrum.components()[0].law
    else {
        panic!("original tabulated law changed");
    };
    assert_eq!(points[1].density_per_j, 1.001 / EV_J);
}

#[test]
fn g1_lf5_shape_scaling_matches_original_tabulated_finite_history() {
    let (tabulated, initial) = material_fixture(100.0);
    let (evaporation, scaled) = material_fixture_general_evaporation(100.0);
    assert_eq!(evaporation.projection(), tabulated.projection());
    assert_eq!(evaporation.neutron_yields(), tabulated.neutron_yields());
    let (expected, original) = initial
        .trial_with_rates(2.0, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let (actual, transformed) = scaled
        .trial_with_rates(2.0, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        actual.neutrons().neutron_counts(),
        expected.neutrons().neutron_counts()
    );
    assert_eq!(
        actual.inventories()[0].populations(),
        expected.inventories()[0].populations()
    );
    assert_eq!(transformed.released, original.released);
    assert_eq!(transformed.decay_neutrons, original.decay_neutrons);
    assert_eq!(
        transformed.transferred_decay_neutron_energy_j,
        original.transferred_decay_neutron_energy_j
    );
    let component = &evaporation.spontaneous_sources()[0].spectrum.components()[0];
    assert_eq!(component.restriction_energy_j, Some(-3.0 * EV_J));
    let SpontaneousFissionSpectrumLaw::GeneralEvaporation {
        temperature_j,
        shape,
        ..
    } = &component.law
    else {
        panic!("original LF5 law changed");
    };
    assert_eq!(*temperature_j, 2.0 * EV_J);
    assert_eq!(shape[1].density, 2.0);
    assert_eq!(shape[1].coordinate, 1.0);
}

#[test]
fn g0_lf5_retains_positive_endpoints_normalization_and_finite_cutoff() {
    let original = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.5);
    let build = |temperature_ev: f64, cutoff_ev: f64, shape: &[(f64, f64)]| {
        EvaluatedBurnupNeutrons::new(
            "original-evaporation-normalization".into(),
            Arc::clone(&original),
            kinetic(Arc::clone(original.model().groups())),
            vec![general_evaporation(temperature_ev, -cutoff_ev, shape, 1)],
            &[],
            source_budget(),
            || Ok(()),
        )
    };
    // LF5 histogram endpoints need not vanish. The Jacobian preserves count
    // and first moment under a different theta/shape coordinate pair.
    let a = build(2.0, 3.0, &[(0.5, 1.0), (1.5, 1.0)]).unwrap();
    let b = build(4.0, 3.0, &[(0.25, 2.0), (0.75, 2.0)]).unwrap();
    assert_eq!(a.projection(), b.projection());
    close(a.projection()[0].particles_per_decay, 0.5, 1e-15);
    close(a.projection()[0].kinetic_energy_j_per_decay, EV_J, 1e-33);
    assert!(matches!(
        build(2.0, 3.0, &[(0.5, 0.9), (1.5, 0.9)]),
        Err(EvaluatedBurnupNeutronError::Normalization { .. })
    ));
    assert!(matches!(
        build(2.0, 2.5, &[(0.5, 1.0), (1.5, 1.0)]),
        Err(EvaluatedBurnupNeutronError::Source(_))
    ));
}

#[test]
fn g1_lf7_lf9_moments_follow_actual_fission_births_and_original_finite_stock() {
    // Independent 70-digit incomplete-gamma means for theta=2 eV, R=3 eV.
    for (law, mean_ev) in [(7, 1.4794210037284057), (9, 1.7292098607415109)] {
        let (source, initial) = material_fixture_with_source(100.0, analytic_source(law, 3.0, 1.0));
        let component = &source.spontaneous_sources()[0].spectrum.components()[0];
        assert_eq!(component.incident_energies_j.as_ref(), &[0.0]);
        assert_eq!(component.restriction_energy_j, Some(-3.0 * EV_J));
        close(source.normalization()[0].component_integrals[0], 1.0, 2e-14);
        close(source.projection()[0].particles_per_decay, 0.5, 2e-14);
        close(
            source.projection()[0].kinetic_energy_j_per_decay,
            0.5 * mean_ev * EV_J,
            1e-32,
        );
        let (after, report) = initial
            .trial_with_rates(2.0, &[300.0], &[0.0; 2], budget(), || Ok(()))
            .unwrap();
        for (&actual, expected) in after.inventories()[0]
            .populations()
            .iter()
            .zip([50.0, 37.5, 6.25, 6.25, 12.5])
        {
            close(actual, expected, 1e-9);
        }
        let lambda = core::f64::consts::LN_2 / 2.0;
        let capture_rate = 0.3;
        let a = (0.5 - (-2.0_f64 * capture_rate).exp()) / (capture_rate - lambda);
        let b = (0.25 - (-2.0_f64 * capture_rate).exp()) / (capture_rate - 2.0 * lambda);
        let prompt = 50.0 * lambda * a;
        let delayed = 50.0 * lambda * (a - b);
        let retained = prompt + delayed;
        close(
            after.neutrons().neutron_counts()[0],
            prompt * (3.0 - mean_ev) / 2.0 + delayed / 2.0,
            1e-9,
        );
        close(
            after.neutrons().neutron_counts()[1],
            prompt * (mean_ev - 1.0) / 2.0 + delayed / 2.0,
            1e-9,
        );
        close(report.decay_neutrons, 31.25, 1e-9);
        close(report.spatial.neutrons.captures, 31.25 - retained, 1e-9);
        close(
            report.transferred_decay_neutron_energy_j,
            (25.0 * mean_ev + 12.5) * EV_J,
            1e-27,
        );
        close(
            after.kinetic_energy_j().unwrap(),
            (mean_ev * prompt + 2.0 * delayed) * EV_J,
            1e-27,
        );
        close(
            after.kinetic_energy_j().unwrap() + report.released.total_j(),
            475.0 * EV_J,
            1e-27,
        );
        let events = &report.spatial.burnup[0].reaction_counts;
        close(report.decay_neutrons, 2.0 * events[0] + events[2], 1e-9);
        assert!(Arc::ptr_eq(
            after.inventories()[0].model(),
            source.assembly().model()
        ));
    }
}

#[test]
fn g0_analytic_laws_do_not_repair_probabilities_or_discard_cutoff_support() {
    let original = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.5);
    let prepare = |law, cutoff, probability| {
        EvaluatedBurnupNeutrons::new(
            "original-analytic-source-boundary".into(),
            Arc::clone(&original),
            kinetic(Arc::clone(original.model().groups())),
            vec![analytic_source(law, cutoff, probability)],
            &[],
            source_budget(),
            || Ok(()),
        )
    };
    for law in [7, 9] {
        assert!(matches!(
            prepare(law, 3.0, 0.9),
            Err(EvaluatedBurnupNeutronError::Normalization { .. })
        ));
        assert!(matches!(
            prepare(law, 5.0, 1.0),
            Err(EvaluatedBurnupNeutronError::Source(_))
        ));
    }
}

#[test]
fn g0_sf_and_ordinary_parent_spectra_must_be_complete_and_disjoint() {
    let mixed = assembly(&[(6.0, 0.25), (5.0, 0.75)], Some((5.0, 0.75)), true, 0.5);
    let source = prepare_source(mixed).unwrap();
    // 0.25*2 prompt + parent-wide ordinary yield 0.75; no second branch factor.
    close(source.projection()[0].particles_per_decay, 1.25, 3e-15);
    close(
        source.projection()[0].kinetic_energy_j_per_decay,
        2.5 * EV_J,
        1e-33,
    );
    for line in [None, Some((6.0, 0.0))] {
        assert!(matches!(
            prepare_source(assembly(&[(6.0, 0.25), (5.0, 0.75)], line, true, 0.5)),
            Err(EvaluatedBurnupNeutronError::Projection(
                DecayNeutronBindingError::Unreported { parent: 0, .. }
            ))
        ));
    }
    for line in [(6.0, 0.75), (10.0, 0.75)] {
        assert!(matches!(
            prepare_source(assembly(&[(6.0, 0.25), (5.0, 0.75)], Some(line), true, 0.5)),
            Err(EvaluatedBurnupNeutronError::Invalid(_))
        ));
    }
    let missing_fragment = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, false, 0.5);
    assert!(matches!(
        prepare_source(Arc::clone(&missing_fragment)),
        Err(EvaluatedBurnupNeutronError::Projection(
            DecayNeutronBindingError::Unreported { parent: 2, .. }
        ))
    ));
    assert!(matches!(
        EvaluatedBurnupNeutrons::new(
            "contradictory-zero".into(),
            Arc::clone(&missing_fragment),
            kinetic(Arc::clone(missing_fragment.model().groups())),
            vec![prompt(selection(), 1.0, 1.0)],
            &[2],
            source_budget(),
            || Ok(()),
        ),
        Err(EvaluatedBurnupNeutronError::Invalid(_))
    ));
}

#[test]
fn g0_prompt_identity_and_original_heavy_energy_are_required() {
    let original = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.5);
    let build = |binding| {
        EvaluatedBurnupNeutrons::new(
            "explicit-prompt-source".into(),
            Arc::clone(&original),
            kinetic(Arc::clone(original.model().groups())),
            vec![binding],
            &[],
            source_budget(),
            || Ok(()),
        )
    };
    let mut mismatched = prompt(selection(), 1.0, 1.0);
    let mut different = selection();
    different.material_number = 778;
    mismatched.spectrum = prompt(different, 1.0, 1.0).spectrum;
    assert!(matches!(
        build(mismatched),
        Err(EvaluatedBurnupNeutronError::Invalid(_))
    ));

    let text = prompt_fixture(
        selection(),
        2.0,
        &[(1.0, vec![(1.0, 0.0), (2.0, 1.0), (3.0, 0.0)])],
    );
    let total_text = text
        .lines()
        .map(|line| {
            let mut line = line.to_owned();
            if line[70..72].trim() == "1" && line[72..75].trim() == "456" {
                line.replace_range(72..75, "452");
            }
            line
        })
        .collect::<Vec<_>>()
        .join("\n");
    let mut total = prompt(selection(), 1.0, 1.0);
    total.multiplicity = Arc::new(
        SpontaneousFissionNeutronYield::parse(
            "explicit-total-is-not-prompt".into(),
            selection(),
            FissionNeutronKind::Total,
            &total_text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    assert!(matches!(
        build(total),
        Err(EvaluatedBurnupNeutronError::Invalid(_))
    ));
    // Depositing all heavy energy leaves ample total/gamma energy but none of
    // the designated heavy radiation account can pay for the emitted neutrons.
    assert!(matches!(
        prepare_source(assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.0)),
        Err(EvaluatedBurnupNeutronError::Projection(
            DecayNeutronBindingError::RadiationBudget { parent: 0, .. }
        ))
    ));
}

#[test]
fn g0_near_closed_branching_preserves_actual_sf_event_multiplicity() {
    let source = Arc::new(
        prepare_source(assembly(
            &[(6.0, 0.25), (4.0, 0.75), (4.0, 1e-15)],
            None,
            true,
            0.5,
        ))
        .unwrap(),
    );
    let original = &source.normalization()[0];
    assert!(original.total_decay_branch_sum > 1.0);
    assert_eq!(original.spontaneous_branch_sum, 0.25);
    assert_eq!(
        original.fission_events_per_decay,
        0.25 / original.total_decay_branch_sum
    );
    close(
        source.projection()[0].particles_per_decay,
        2.0 * original.fission_events_per_decay,
        2e-16,
    );
    let initial = activation(&source, 100.0);
    let (_, report) = initial
        .trial_with_rates(1.0, &[300.0], &[0.0; 2], budget(), || Ok(()))
        .unwrap();
    let events = &report.spatial.burnup[0].reaction_counts;
    close(report.decay_neutrons, 2.0 * events[0] + events[3], 1e-12);
}

#[test]
fn g0_binding_retains_provenance_and_only_admits_original_decay_channels() {
    let (source, _) = material_fixture(100.0);
    let original = source.assembly().model();
    let mut reactions = original.reactions().to_vec();
    reactions.reverse();
    reactions.push(BurnupReaction {
        data_id: "extra-finite-neutron-capture".into(),
        parent: 0,
        rate: BurnupRate::Capture {
            cross_sections_m2: vec![0.01; 2],
        },
        products: vec![],
        energy_per_event: EnergyChannels::default(),
    });
    let data = BurnupData {
        data_id: "same-decay-plus-irradiation".into(),
        groups: Arc::clone(original.groups()),
        temperature_min_k: 300.0,
        temperature_max_k: 600.0,
        nuclides: original.nuclides().to_vec(),
        reactions,
    };
    let model = Arc::new(BurnupModel::new(data.clone(), BurnupLimits::default()).unwrap());
    let bound = source
        .bind_burnup("bound-source".into(), Arc::clone(&model), || Ok(()))
        .unwrap();
    assert!(Arc::ptr_eq(bound.source(), &source));
    assert!(Arc::ptr_eq(bound.emission().model(), &model));
    assert_eq!(bound.emission().neutron_yields(), source.neutron_yields());
    let mut changed = data;
    changed.reactions[0].energy_per_event.radiation_j += EV_J;
    let changed = Arc::new(BurnupModel::new(changed, BurnupLimits::default()).unwrap());
    assert!(matches!(
        source.bind_burnup("altered-source".into(), changed, || Ok(())),
        Err(EvaluatedBurnupNeutronError::Assembly(_))
    ));
    assert!(matches!(
        source.bind_burnup("cancelled-source".into(), model, || Err(
            FissionError::Cancelled
        )),
        Err(EvaluatedBurnupNeutronError::Source(_))
    ));
}

#[test]
fn g0_preparation_bounds_and_late_cancellation_publish_no_partial_source() {
    let original = assembly(&[(6.0, 0.25), (4.0, 0.75)], None, true, 0.5);
    let model = kinetic(Arc::clone(original.model().groups()));
    let binding = prompt(selection(), 1.0, 1.0);
    let mut polls = 0;
    let source = EvaluatedBurnupNeutrons::new(
        "bounded-source".into(),
        Arc::clone(&original),
        Arc::clone(&model),
        vec![binding.clone()],
        &[],
        source_budget(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut remaining = polls - 1;
    assert!(matches!(
        EvaluatedBurnupNeutrons::new(
            "late-cancelled-source".into(),
            Arc::clone(&original),
            Arc::clone(&model),
            vec![binding.clone()],
            &[],
            source_budget(),
            || if remaining == 0 {
                Err(FissionError::Cancelled)
            } else {
                remaining -= 1;
                Ok(())
            },
        ),
        Err(EvaluatedBurnupNeutronError::Source(_))
    ));
    for limits in [
        DecayEmissionLimits {
            max_groups: 1,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_group_values: 1,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_work: 1,
            ..DecayEmissionLimits::default()
        },
    ] {
        assert!(
            EvaluatedBurnupNeutrons::new(
                "refused-bounded-source".into(),
                Arc::clone(&original),
                Arc::clone(&model),
                vec![binding.clone()],
                &[],
                EvaluatedBurnupNeutronBudget {
                    limits,
                    ..source_budget()
                },
                || Ok(()),
            )
            .is_err()
        );
    }
    assert!(Arc::ptr_eq(source.assembly(), &original));
    assert!(Arc::ptr_eq(
        &source.spontaneous_sources()[0].multiplicity,
        &binding.multiplicity
    ));
    assert!(Arc::ptr_eq(
        &source.spontaneous_sources()[0].spectrum,
        &binding.spectrum
    ));
}
