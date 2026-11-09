//! Actual parsed MF1/455 and MF5/455 families through finite native histories.
//! Coefficients are manufactured; no isotope calibration is claimed.

use super::*;
use crate::diffusion::DelayedGroup;
use crate::fusion_neutrons::NeutronReactionEnergy;

fn integers(values: [i32; 6], mf: u16, mt: u16) -> String {
    format!(
        "{}{:4}{mf:2}{mt:3}{:5}",
        values.map(|value| format!("{value:11}")).concat(),
        777,
        1
    )
}

fn delayed_records(probabilities: [f64; 2], areas: [f64; 2]) -> String {
    let mut lines = vec![
        control(za(parent()), 202.0, [-1, 0, 1, 3], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 6], 1, 451),
        control(0.0, 0.0, [2, 0, 4, 8], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 0),
        control(za(parent()), 202.0, [0, 1, 0, 0], 1, 455),
        control(0.0, 0.0, [0, 0, 2, 0], 1, 455),
        row([1.0, 2.0, 0.0, 0.0, 0.0, 0.0], 1, 455),
        control(0.0, 0.0, [0, 0, 1, 0], 1, 455),
        row([0.5, 0.0, 0.0, 0.0, 0.0, 0.0], 1, 455),
        control(0.0, 0.0, [0; 4], 1, 0),
        control(za(parent()), 202.0, [0, 0, 2, 0], 5, 455),
    ];
    for (index, (&probability, &area)) in probabilities.iter().zip(&areas).enumerate() {
        lines.push(control(0.0, 0.0, [0, 1, 1, 2], 5, 455));
        lines.push(integers([2, 2, 0, 0, 0, 0], 5, 455));
        lines.push(row([0.0, probability, 1.0, probability, 0.0, 0.0], 5, 455));
        lines.push(control(0.0, 0.0, [0, 0, 1, 1], 5, 455));
        lines.push(integers([1, 2, 0, 0, 0, 0], 5, 455));
        lines.push(control(0.0, 0.0, [0, 0, 1, 3], 5, 455));
        lines.push(integers([3, 2, 0, 0, 0, 0], 5, 455));
        let center = 2.0 + 2.0 * index as f64;
        lines.push(row(
            [center - 1.0, 0.0, center, area, center + 1.0, 0.0],
            5,
            455,
        ));
    }
    lines.push(control(0.0, 0.0, [0; 4], 5, 0));
    lines.join("\n")
}

fn parsed_delayed(
    probabilities: [f64; 2],
    areas: [f64; 2],
) -> (
    Arc<SpontaneousFissionNeutronYield>,
    Arc<SpontaneousFissionNeutronSpectrum>,
) {
    let text = delayed_records(probabilities, areas);
    (
        Arc::new(
            SpontaneousFissionNeutronYield::parse(
                "manufactured-MF1-455-family-multiplicity".into(),
                selection(),
                FissionNeutronKind::Delayed,
                &text,
                SpontaneousNeutronReadLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        ),
        Arc::new(
            SpontaneousFissionNeutronSpectrum::parse(
                "manufactured-MF5-455-family-spectra".into(),
                selection(),
                FissionNeutronKind::Delayed,
                &text,
                SpontaneousNeutronReadLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        ),
    )
}

fn exported_assembly(
    retain_fragment: bool,
    heavy_radiation: f64,
) -> Arc<EvaluatedBurnupDecayAssembly> {
    let mut entries = vec![
        DecayModelEntry {
            evaluation: evaluation(
                parent(),
                Some(2.0),
                [1.0, 2.0, 6.0],
                &[(6.0, 0.25), (4.0, 0.75)],
                Some((6.0, 0.0)),
            ),
            routing: routing(heavy_radiation),
        },
        DecayModelEntry {
            evaluation: evaluation(id(78, 198), None, [0.0; 3], &[], None),
            routing: routing(1.0),
        },
    ];
    let mut exported_products = vec![id(39, 99), id(41, 101)];
    if retain_fragment {
        entries.push(DecayModelEntry {
            evaluation: evaluation(id(39, 99), None, [0.0; 3], &[], None),
            routing: routing(1.0),
        });
        exported_products = vec![id(41, 101)];
    }
    Arc::new(
        EvaluatedBurnupDecayAssembly::new(
            "manufactured-SF-fragments-outside-native-family-cohort".into(),
            entries,
            vec![SpontaneousFissionProductBinding {
                parent: parent(),
                yields: fission_products(),
                exported_products,
            }],
            Arc::new(
                NeutronGroups::new(
                    "manufactured-delayed-SF-groups".into(),
                    vec![0.0, 3.0 * EV_J, 6.0 * EV_J],
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

fn kinetic_with_families(
    groups: Arc<NeutronGroups>,
    families: Vec<DelayedGroup>,
    energies: [f64; 2],
) -> Arc<NeutronKineticModel> {
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-delayed-SF-native-capture".into(),
                groups,
                speeds_m_s: vec![1.0; 2],
                cell_volumes_m3: vec![0.01],
                delayed_groups: families,
                reactions: (0..2)
                    .map(|incident_group| DiffusionReaction {
                        data_id: format!("manufactured-background-capture-{incident_group}"),
                        cell: 0,
                        incident_group,
                        macroscopic_cross_section_m_inv: 0.4,
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
        NeutronKineticModel::new_with_precursor_storage(
            "manufactured-delayed-SF-energy-owner".into(),
            diffusion,
            energies.to_vec(),
            vec![0.25; 2],
            vec![NeutronReactionEnergy::default(); 2],
        )
        .unwrap(),
    )
}

fn prepared(
    groups: Arc<NeutronGroups>,
    probabilities: [f64; 2],
    areas: [f64; 2],
    budget: EvaluatedBurnupNeutronBudget,
) -> Arc<EvaluatedSpontaneousFissionDelayedNeutrons> {
    let (multiplicity, spectrum) = parsed_delayed(probabilities, areas);
    let grid = kinetic_with_families(groups, vec![], [EV_J, 5.0 * EV_J]);
    Arc::new(
        EvaluatedSpontaneousFissionDelayedNeutrons::new(
            "manufactured-prepared-delayed-SF".into(),
            multiplicity,
            spectrum,
            &grid,
            budget,
            || Ok(()),
        )
        .unwrap(),
    )
}

fn immediate_source(
    assembly: Arc<EvaluatedBurnupDecayAssembly>,
    prepared: &EvaluatedSpontaneousFissionDelayedNeutrons,
    families: Vec<DelayedGroup>,
    energies: [f64; 2],
) -> EvaluatedBurnupNeutrons {
    let kinetic = kinetic_with_families(Arc::clone(assembly.model().groups()), families, energies);
    assert_eq!(prepared.families().len(), 2);
    EvaluatedBurnupNeutrons::new(
        "manufactured-original-prompt-plus-delayed-source".into(),
        assembly,
        kinetic,
        vec![prompt(selection(), 1.0, 1.0)],
        &[],
        source_budget(),
        || Ok(()),
    )
    .unwrap()
}

fn binding(
    source: Arc<EvaluatedSpontaneousFissionDelayedNeutrons>,
    cap_ev: f64,
) -> DelayedSpontaneousFissionBinding {
    DelayedSpontaneousFissionBinding {
        parent: parent(),
        source,
        native_family_indices: vec![0, 1],
        energy_funding: DelayedSpontaneousFissionFunding::ParentHeavyRadiation {
            maximum_energy_j_per_fission: cap_ev * EV_J,
        },
    }
}

/// Shared parsed evaluated source for the existing material heating driver.
/// Parent half-life2s, SFbranch.25, prompt2 neutrons/SF at2eV.
/// Delayed nu=.5 with original p=[.4,.6], lambda=[1,2]/s and means[2,4]eV.
/// All SF fragments are exported. Parent categories[1,2,6]eV route to heat,
/// radiation, radiation. Capture=.4/s with .25 local heat; volume=.01m3.
pub(crate) fn material_delayed_sf_fixture(
    parent_count: f64,
) -> (Arc<EvaluatedBurnupNeutrons>, FusionNeutronActivationState) {
    let assembly = exported_assembly(false, 1.0);
    let prepared = prepared(
        Arc::clone(assembly.model().groups()),
        [0.4, 0.6],
        [1.0, 1.0],
        source_budget(),
    );
    let original = immediate_source(
        Arc::clone(&assembly),
        &prepared,
        prepared.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    let source = Arc::new(
        original
            .with_delayed_spontaneous_fission(vec![binding(prepared, 1.7)], source_budget(), || {
                Ok(())
            })
            .unwrap(),
    );
    let bound = source
        .bind_burnup(
            "manufactured-evaluated-finite-delayed-SF".into(),
            Arc::clone(assembly.model()),
            || Ok(()),
        )
        .unwrap();
    let spatial = DiffusionBurnupState::new(
        DiffusionState::new(
            Arc::clone(source.kinetic_model().diffusion()),
            vec![0.0; 2],
            vec![0.0; 2],
        )
        .unwrap(),
        vec![BurnupInventory::new(Arc::clone(assembly.model()), vec![parent_count, 0.0]).unwrap()],
        vec![],
    )
    .unwrap()
    .with_decay_neutrons(vec![Arc::clone(bound.emission())])
    .unwrap();
    let state =
        FusionNeutronActivationState::new(Arc::clone(source.kinetic_model()), spatial).unwrap();
    (source, state)
}

fn native_budget() -> FusionNeutronActivationBudget {
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
            reaction_tolerance: 1e-10,
            relaxation: 1.0,
        },
        energy_tolerance_j: 1e-24,
    }
}

#[test]
fn g1_original_family_order_and_density_discrepancies_survive_preparation() {
    let assembly = exported_assembly(false, 1.0);
    let budget = EvaluatedBurnupNeutronBudget {
        normalization_tolerance: 0.03,
        ..source_budget()
    };
    let source = prepared(
        Arc::clone(assembly.model().groups()),
        [0.4, 0.6],
        [0.98, 1.02],
        budget,
    );
    assert_eq!(source.multiplicity().decay_constants_s_inv(), &[1.0, 2.0]);
    assert_eq!(source.families()[0].decay_constant_s_inv, 1.0);
    assert_eq!(source.families()[1].decay_constant_s_inv, 2.0);
    assert_eq!(source.spectrum().components()[0].probability, 0.4);
    close(source.component_integrals()[0], 0.98, 1e-12);
    close(source.component_integrals()[1], 1.02, 1e-12);
    close(source.component_probability_sum(), 1.0, 1e-12);
    close(source.weighted_spectrum_integral(), 1.004, 1e-12);
    close(source.precursor_yields_per_fission()[0], 0.196, 1e-12);
    close(source.precursor_yields_per_fission()[1], 0.306, 1e-12);
    close(
        source.precursor_yields_per_fission().iter().sum(),
        0.502,
        1e-12,
    );
    close(source.source_energy_j_per_fission(), 1.616 * EV_J, 1e-29);
    for (family, expected_mean_ev) in source.families().iter().zip([2.0, 4.0]) {
        close(family.spectrum.iter().sum(), 1.0, 1e-12);
        close(
            family.spectrum[0] + 5.0 * family.spectrum[1],
            expected_mean_ev,
            1e-12,
        );
    }
    let original = immediate_source(
        Arc::clone(&assembly),
        &source,
        source.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    let attached = original
        .clone()
        .with_delayed_spontaneous_fission(
            vec![binding(Arc::clone(&source), 1.7)],
            budget,
            || Ok(()),
        )
        .unwrap();
    close(attached.precursor_yields()[0][0], 0.049, 1e-12);
    close(attached.precursor_yields()[0][1], 0.0765, 1e-12);
    close(
        attached.delayed_energy_j_per_decay()[0],
        0.404 * EV_J,
        1e-29,
    );
    assert!(Arc::ptr_eq(&attached.delayed_sources()[0].source, &source));
    assert!(
        original
            .with_delayed_spontaneous_fission(
                vec![binding(source, 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
}

fn response(rate: f64, capture: f64, time: f64) -> f64 {
    ((-rate * time).exp() - (-capture * time).exp()) / (capture - rate)
}

#[test]
fn g1_parsed_finite_sf_delays_births_and_reserves_energy_only_once() {
    let (source, initial) = material_delayed_sf_fixture(100.0);
    let time = 1.25_f64;
    let lambda = core::f64::consts::LN_2 / 2.0;
    let decays = 100.0 * -(-lambda * time).exp_m1();
    let prompt = 0.5 * lambda * 100.0 * response(lambda, 0.4, time);
    let mut expected_neutrons = prompt;
    let mut expected_neutron_energy_ev = 2.0 * prompt;
    let mut expected_precursor_energy_ev = 0.0;
    let mut delayed_energy_ev = 0.0;
    let mut expected_precursors = Vec::new();
    for (yield_, beta, mean) in [(0.05, 1.0, 2.0), (0.075, 2.0, 4.0)] {
        let precursor = yield_ * lambda * 100.0 * ((-lambda * time).exp() - (-beta * time).exp())
            / (beta - lambda);
        let neutrons = yield_ * lambda * beta * 100.0 / (beta - lambda)
            * (response(lambda, 0.4, time) - response(beta, 0.4, time));
        expected_neutrons += neutrons;
        expected_neutron_energy_ev += mean * neutrons;
        expected_precursor_energy_ev += mean * precursor;
        delayed_energy_ev += mean * (yield_ * decays - precursor);
        expected_precursors.push(precursor);
    }
    let (after, step) = initial
        .trial_with_rates(time, &[300.0], &[0.0; 2], native_budget(), || Ok(()))
        .unwrap();
    close(
        after.inventories()[0].populations()[0],
        100.0 - decays,
        1e-9,
    );
    close(after.inventories()[0].populations()[1], 0.75 * decays, 1e-9);
    close(
        after.neutrons().neutron_counts().iter().sum(),
        expected_neutrons,
        1e-9,
    );
    for (&got, &expected) in after
        .neutrons()
        .precursor_counts()
        .iter()
        .zip(&expected_precursors)
    {
        close(got, expected, 1e-9);
    }
    close(
        after.kinetic_energy_j().unwrap(),
        expected_neutron_energy_ev * EV_J,
        1e-27,
    );
    close(
        after.precursor_energy_j().unwrap(),
        expected_precursor_energy_ev * EV_J,
        1e-27,
    );
    close(step.decay_neutrons, 0.5 * decays, 1e-9);
    close(step.decay_precursors, 0.125 * decays, 1e-9);
    close(
        step.spatial.decay_precursor_counts.as_ref().unwrap()[0][0],
        0.05 * decays,
        1e-9,
    );
    close(
        step.spatial.decay_precursor_counts.as_ref().unwrap()[0][1],
        0.075 * decays,
        1e-9,
    );
    close(
        step.reserved_decay_precursor_energy_j,
        0.4 * decays * EV_J,
        1e-27,
    );
    close(
        step.transferred_decay_neutron_energy_j,
        decays * EV_J,
        1e-27,
    );
    let captured_energy_ev = decays + delayed_energy_ev - expected_neutron_energy_ev;
    close(
        step.released.deposited_j,
        (decays + 0.25 * captured_energy_ev) * EV_J,
        1e-27,
    );
    close(
        step.released.radiation_j,
        (6.6 * decays + 0.75 * captured_energy_ev) * EV_J,
        1e-27,
    );
    close(
        after.stored_energy_j().unwrap() + step.released.total_j(),
        9.0 * decays * EV_J,
        1e-27,
    );
    close(step.nuclear_generated_energy_j, 9.0 * decays * EV_J, 1e-27);
    assert_eq!(step.external_injected_neutrons, 0.0);
    assert_eq!(step.external_injected_precursors, 0.0);
    assert_eq!(source.assembly().nuclides().len(), 2);
    assert_eq!(initial.neutrons().precursor_counts(), &[0.0; 2]);
    let (middle, first) = initial
        .trial_with_rates(0.5, &[300.0], &[0.0; 2], native_budget(), || Ok(()))
        .unwrap();
    let (split, second) = middle
        .trial_with_rates(0.75, &[300.0], &[0.0; 2], native_budget(), || Ok(()))
        .unwrap();
    close(
        split.stored_energy_j().unwrap(),
        after.stored_energy_j().unwrap(),
        1e-27,
    );
    close(
        first.released.total_j() + second.released.total_j(),
        step.released.total_j(),
        1e-27,
    );
    for (&got, &expected) in split
        .neutrons()
        .precursor_counts()
        .iter()
        .zip(after.neutrons().precursor_counts())
    {
        close(got, expected, 1e-9);
    }
}

#[test]
fn g0_retained_fragments_family_aliasing_grid_and_unfunded_energy_refuse() {
    let assembly = exported_assembly(false, 1.0);
    let source = prepared(
        Arc::clone(assembly.model().groups()),
        [0.4, 0.6],
        [1.0, 1.0],
        source_budget(),
    );
    let original = immediate_source(
        Arc::clone(&assembly),
        &source,
        source.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    for cap in [0.0, 1.0, f64::NAN] {
        assert!(
            original
                .clone()
                .with_delayed_spontaneous_fission(
                    vec![binding(Arc::clone(&source), cap)],
                    source_budget(),
                    || Ok(()),
                )
                .is_err()
        );
    }
    for indices in [vec![0, 0], vec![1, 0], vec![0, 2], vec![0]] {
        let mut wrong = binding(Arc::clone(&source), 1.7);
        wrong.native_family_indices = indices;
        assert!(
            original
                .clone()
                .with_delayed_spontaneous_fission(vec![wrong], source_budget(), || Ok(()),)
                .is_err()
        );
    }
    let retained = exported_assembly(true, 1.0);
    let retained = immediate_source(
        retained,
        &source,
        source.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    assert!(
        retained
            .with_delayed_spontaneous_fission(
                vec![binding(Arc::clone(&source), 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
    // Original prompt neutron energy fits 1.2eV/parent, but the additional
    // 0.4eV reservation cannot borrow local heat or electromagnetic radiation.
    let insufficient = exported_assembly(false, 0.2);
    let insufficient = immediate_source(
        insufficient,
        &source,
        source.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    assert!(
        insufficient
            .with_delayed_spontaneous_fission(
                vec![binding(Arc::clone(&source), 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
    let wrong_energy = immediate_source(
        Arc::clone(&assembly),
        &source,
        source.families().to_vec(),
        [EV_J, 4.5 * EV_J],
    );
    assert!(
        wrong_energy
            .with_delayed_spontaneous_fission(
                vec![binding(Arc::clone(&source), 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
    let other_grid = Arc::new(
        NeutronGroups::new(
            "different-grid-identity".into(),
            vec![0.0, 3.0 * EV_J, 6.0 * EV_J],
        )
        .unwrap(),
    );
    let wrong_grid = prepared(other_grid, [0.4, 0.6], [1.0, 1.0], source_budget());
    assert!(
        original
            .clone()
            .with_delayed_spontaneous_fission(
                vec![binding(wrong_grid, 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
    let mut wrong_family = source.families().to_vec();
    wrong_family[0].decay_constant_s_inv = 1.5;
    let wrong_family = immediate_source(assembly, &source, wrong_family, [EV_J, 5.0 * EV_J]);
    assert!(
        wrong_family
            .with_delayed_spontaneous_fission(
                vec![binding(source, 1.7)],
                source_budget(),
                || Ok(()),
            )
            .is_err()
    );
    assert!(original.delayed_sources().is_empty());
}

#[test]
fn g0_g4_normalization_work_and_late_cancellation_publish_no_partial_binding() {
    let assembly = exported_assembly(false, 1.0);
    let grid = kinetic_with_families(
        Arc::clone(assembly.model().groups()),
        vec![],
        [EV_J, 5.0 * EV_J],
    );
    for (probabilities, areas) in [([0.4, 0.7], [1.0, 1.0]), ([0.4, 0.6], [0.9, 1.1])] {
        let (multiplicity, spectrum) = parsed_delayed(probabilities, areas);
        assert!(
            EvaluatedSpontaneousFissionDelayedNeutrons::new(
                "unrepaired-delayed-normalization".into(),
                multiplicity,
                spectrum,
                &grid,
                source_budget(),
                || Ok(()),
            )
            .is_err()
        );
    }
    let (multiplicity, spectrum) = parsed_delayed([0.4, 0.6], [1.0, 1.0]);
    let mut polls = 0;
    let source = Arc::new(
        EvaluatedSpontaneousFissionDelayedNeutrons::new(
            "manufactured-prepared-delayed-SF".into(),
            Arc::clone(&multiplicity),
            Arc::clone(&spectrum),
            &grid,
            source_budget(),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap(),
    );
    let mut seen = 0;
    assert!(
        EvaluatedSpontaneousFissionDelayedNeutrons::new(
            "manufactured-prepared-delayed-SF".into(),
            Arc::clone(&multiplicity),
            Arc::clone(&spectrum),
            &grid,
            source_budget(),
            || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .is_err()
    );
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
            EvaluatedSpontaneousFissionDelayedNeutrons::new(
                "bounded-delayed-SF".into(),
                Arc::clone(&multiplicity),
                Arc::clone(&spectrum),
                &grid,
                EvaluatedBurnupNeutronBudget {
                    limits,
                    ..source_budget()
                },
                || Ok(()),
            )
            .is_err()
        );
    }
    let original = immediate_source(
        assembly,
        &source,
        source.families().to_vec(),
        [EV_J, 5.0 * EV_J],
    );
    let mut polls = 0;
    let accepted = original
        .clone()
        .with_delayed_spontaneous_fission(
            vec![binding(Arc::clone(&source), 1.7)],
            source_budget(),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(
        original
            .clone()
            .with_delayed_spontaneous_fission(vec![binding(source, 1.7)], source_budget(), || {
                seen += 1;
                if seen == polls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },)
            .is_err()
    );
    assert!(original.delayed_sources().is_empty());
    assert_eq!(accepted.delayed_sources().len(), 1);
    assert!(Arc::ptr_eq(accepted.assembly(), original.assembly()));
    assert!(Arc::ptr_eq(
        accepted.kinetic_model(),
        original.kinetic_model()
    ));
}

#[test]
fn g2_existing_cf252_six_family_lf5_source_prepares_without_normalization_repair() {
    // Existing original ENDF/B-VIII.0 MF1/MF5 fixture; no MF8 inventory or
    // calorimetric funding is inferred for this preparation-only regression.
    let text =
        include_str!("../../../../spontaneous_neutrons/tests/cf252_endfb_viii_0_neutrons.endf");
    let selection = EndfSpontaneousNeutronSelection {
        material_number: 3644,
        nuclide: NuclideId::new(98, 252, 0).unwrap(),
    };
    let multiplicity = Arc::new(
        SpontaneousFissionNeutronYield::parse(
            "ENDF-B-VIII.0-Cf252-delayed-yield".into(),
            selection,
            FissionNeutronKind::Delayed,
            text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let spectrum = Arc::new(
        SpontaneousFissionNeutronSpectrum::parse(
            "ENDF-B-VIII.0-Cf252-delayed-spectrum".into(),
            selection,
            FissionNeutronKind::Delayed,
            text,
            SpontaneousNeutronReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let groups = Arc::new(
        NeutronGroups::new(
            "explicit-Cf252-preparation-grid".into(),
            vec![0.0, 0.5e6 * EV_J, 1.5e6 * EV_J, 3.0e6 * EV_J, 30.0e6 * EV_J],
        )
        .unwrap(),
    );
    let energies = vec![1e-3 * EV_J, 0.75e6 * EV_J, 2.0e6 * EV_J, 29.0e6 * EV_J];
    let native = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "explicit-Cf252-preparation-model".into(),
                groups,
                speeds_m_s: vec![1.0; 4],
                cell_volumes_m3: vec![0.01],
                delayed_groups: vec![],
                reactions: vec![],
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    let grid = NeutronKineticModel::new_with_precursor_storage(
        "explicit-Cf252-preparation-energies".into(),
        native,
        energies.clone(),
        vec![],
        vec![],
    )
    .unwrap();
    let budget = EvaluatedBurnupNeutronBudget {
        normalization_tolerance: 1e-6,
        ..source_budget()
    };
    let source = EvaluatedSpontaneousFissionDelayedNeutrons::new(
        "original-Cf252-six-delayed-families".into(),
        Arc::clone(&multiplicity),
        Arc::clone(&spectrum),
        &grid,
        budget,
        || Ok(()),
    )
    .unwrap();
    // Independent 50-digit integrations of the original decimal histograms,
    // already retained by the original spectra regression, before SI scaling.
    let expected = [
        (1.0000000120373, 525881.5257291035),
        (1.000000035537044, 587883.4445469742),
        (1.000000006764, 380927.25689426),
        (1.00000000646, 449114.1779536),
        (1.00000000055, 421187.90633315),
        (1.00000001174, 493750.5845413),
    ];
    assert_eq!(source.families().len(), 6);
    assert!(Arc::ptr_eq(source.multiplicity(), &multiplicity));
    assert!(Arc::ptr_eq(source.spectrum(), &spectrum));
    let mut expected_total = 0.0;
    let mut expected_energy = 0.0;
    for (index, (area, energy_ev)) in expected.into_iter().enumerate() {
        let probability = spectrum.components()[index].probability;
        let count = multiplicity.yield_per_fission() * probability * area;
        close(source.component_integrals()[index], area, 3e-14);
        assert_eq!(
            source.families()[index].decay_constant_s_inv,
            multiplicity.decay_constants_s_inv()[index]
        );
        close(source.precursor_yields_per_fission()[index], count, 1e-14);
        close(
            source.component_projections()[index].kinetic_energy_j_per_decay / EV_J,
            energy_ev,
            3e-8,
        );
        let conditional_energy_ev: f64 = source.families()[index]
            .spectrum
            .iter()
            .zip(&energies)
            .map(|(weight, energy)| weight * energy / EV_J)
            .sum();
        close(conditional_energy_ev, energy_ev / area, 4e-8);
        expected_total += count;
        expected_energy += multiplicity.yield_per_fission() * probability * energy_ev;
    }
    close(source.component_probability_sum(), 1.00000014, 1e-15);
    close(
        source.precursor_yields_per_fission().iter().sum(),
        expected_total,
        1e-14,
    );
    close(
        source.source_energy_j_per_fission() / EV_J,
        expected_energy,
        1e-8,
    );
    assert!((expected_total - multiplicity.yield_per_fission()).abs() > 1e-12);
    assert!(
        EvaluatedSpontaneousFissionDelayedNeutrons::new(
            "over-tight-Cf252-normalization".into(),
            multiplicity,
            spectrum,
            &grid,
            source_budget(),
            || Ok(()),
        )
        .is_err()
    );
}
