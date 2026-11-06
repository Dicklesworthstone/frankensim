use super::*;
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind, DiffusionState,
};
use crate::endf::spectrum::{FluxSpectrumData, FluxSpectrumLimits};
use crate::irradiation::NuclideId;

fn region(points: usize, interpolation: EndfInterpolation) -> Vec<EndfInterpolationRegion> {
    vec![EndfInterpolationRegion {
        end_point: points,
        interpolation,
    }]
}

fn sigma(energies: &[f64], values: &[f64], law: EndfInterpolation) -> PointwiseCrossSection {
    PointwiseCrossSection {
        data_id: "manufactured-fission-cross-section".into(),
        metadata: EndfMetadata {
            material_number: 1,
            nuclide: NuclideId::new(10, 20, 0).unwrap(),
            temperature_k: 300.0,
            atomic_weight_ratio: 20.0,
            library_id: 0,
            library_version: 1,
            library_release: 0,
            material_revision: 0,
            resonance_flag: 2,
            derivation_flag: 1,
        },
        reaction: NeutronReaction::Fission,
        energies_j: energies.into(),
        cross_sections_m2: values.into(),
        interpolation_regions: region(energies.len(), law).into(),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}

fn spectrum(energies: &[f64], values: &[f64], law: EndfInterpolation) -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-absolute-spectrum".into(),
            energies_j: energies.to_vec(),
            flux_density_m2_s_j: values.to_vec(),
            interpolation_regions: region(energies.len(), law),
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

fn curve_data(energies: &[f64], values: &[f64], law: EndfInterpolation) -> TabulatedYieldData {
    TabulatedYieldData {
        data_id: "manufactured-dimensionless-yield".into(),
        energies_j: energies.to_vec(),
        values: values.to_vec(),
        interpolation_regions: region(energies.len(), law),
    }
}

fn curve(energies: &[f64], values: &[f64], law: EndfInterpolation) -> TabulatedYield {
    TabulatedYield::new(
        curve_data(energies, values, law),
        TabulatedYieldLimits::default(),
    )
    .unwrap()
}

fn groups(edges: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("manufactured-yield-grid".into(), edges.to_vec()).unwrap())
}

fn collapse(
    sigma: &PointwiseCrossSection,
    spectrum: &FluxSpectrum,
    curves: &[TabulatedYield],
    groups: Arc<NeutronGroups>,
) -> ReactionWeightedCollapse {
    collapse_reaction_weighted(
        sigma,
        spectrum,
        curves,
        groups,
        ReactionWeightedLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!(
            (actual - expected).abs() <= 8e-13 * expected.abs(),
            "{actual:.17e} != {expected:.17e}"
        );
    }
}

#[test]
fn cubic_production_integrals_use_reaction_weighting_and_share_original_flux_coefficients() {
    // sigma=phi=y=E on [0,2]: flux=2, reaction=8/3, production=4.
    // The reaction-weighted yield is 3/2; the flux-only mean would be 4/3.
    let cross = sigma(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear);
    let flux = spectrum(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear);
    let curves = [
        curve(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear),
        curve(&[0.0, 2.0], &[2.0, 8.0], EndfInterpolation::Linear),
        curve(&[0.0, 2.0], &[0.0, 0.0], EndfInterpolation::Linear),
    ];
    let grid = groups(&[0.0, 2.0]);
    let answer = collapse(&cross, &flux, &curves, Arc::clone(&grid));
    let original = cross
        .collapse_with_spectrum(
            &flux,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(answer.cross_sections_m2, original.cross_sections_m2);
    assert_eq!(answer.flux, original.flux);
    assert!(Arc::ptr_eq(answer.flux.groups(), &grid));
    near(answer.flux.values_m2_s()[0], 2.0);
    near(answer.cross_sections_m2[0], 4.0 / 3.0);
    for (row, (expected_yield, expected_production)) in
        answer
            .yields
            .iter()
            .zip([(1.5, 4.0), (6.5, 52.0 / 3.0), (0.0, 0.0)])
    {
        near(row.values[0].unwrap(), expected_yield);
        near(
            answer.cross_sections_m2[0] * answer.flux.values_m2_s()[0] * row.values[0].unwrap(),
            expected_production,
        );
        assert_eq!(&*row.data_id, curves[0].data_id());
    }
    assert!((answer.yields[0].values[0].unwrap() - 4.0 / 3.0).abs() > 0.1);
    assert_eq!(&*answer.cross_section_data_id, cross.data_id());
    assert_eq!(&*answer.spectrum_data_id, flux.data_id());
    assert_eq!(answer.metadata, *cross.metadata());
    assert_eq!(answer.reaction, NeutronReaction::Fission);
}

#[test]
fn union_mesh_preserves_histogram_changes_paired_jumps_and_group_boundaries() {
    let cross = sigma(
        &[0.0, 1.0, 1.0, 3.0],
        &[1.0, 1.0, 4.0, 4.0],
        EndfInterpolation::Linear,
    );
    let flux = spectrum(
        &[0.0, 2.0, 3.0],
        &[2.0, 1.0, 1.0],
        EndfInterpolation::Histogram,
    );
    let yield_ = curve(
        &[0.0, 1.5, 1.5, 3.0],
        &[1.0, 4.0, 2.0, 5.0],
        EndfInterpolation::Linear,
    );
    assert_eq!(yield_.value_at_energy(1.5).unwrap(), 2.0);
    near(yield_.value_at_energy(1.5_f64.next_down()).unwrap(), 4.0);
    let fine = collapse(
        &cross,
        &flux,
        std::slice::from_ref(&yield_),
        groups(&[0.0, 1.0, 2.0, 3.0]),
    );
    assert_eq!(fine.flux.values_m2_s(), [2.0, 2.0, 1.0]);
    for group in 0..3 {
        near(fine.cross_sections_m2[group], [1.0, 4.0, 4.0][group]);
        near(
            fine.yields[0].values[group].unwrap(),
            [2.0, 3.0, 4.0][group],
        );
    }
    let coarse = collapse(&cross, &flux, &[yield_], groups(&[0.0, 3.0]));
    near(coarse.cross_sections_m2[0], 14.0 / 5.0);
    near(coarse.yields[0].values[0].unwrap(), 22.0 / 7.0);
    near(
        coarse.cross_sections_m2[0]
            * coarse.flux.values_m2_s()[0]
            * coarse.yields[0].values[0].unwrap(),
        44.0,
    );
}

#[test]
fn undefined_zero_reaction_yields_are_distinct_from_defined_zero_product_yields() {
    let cross = sigma(
        &[0.0, 1.0, 2.0],
        &[0.0, 2.0, 2.0],
        EndfInterpolation::Histogram,
    );
    let flux = spectrum(&[0.0, 2.0], &[1.0, 1.0], EndfInterpolation::Linear);
    let curves = [
        curve(&[0.0, 2.0], &[4.0, 4.0], EndfInterpolation::Linear),
        curve(&[0.0, 2.0], &[0.0, 0.0], EndfInterpolation::Linear),
    ];
    let answer = collapse(&cross, &flux, &curves, groups(&[0.0, 1.0, 2.0]));
    assert_eq!(answer.cross_sections_m2, [0.0, 2.0]);
    assert_eq!(answer.yields[0].values, [None, Some(4.0)]);
    assert_eq!(answer.yields[1].values, [None, Some(0.0)]);
    let zero_flux = spectrum(&[0.0, 2.0], &[0.0, 0.0], EndfInterpolation::Linear);
    assert!(
        collapse_reaction_weighted(
            &cross,
            &zero_flux,
            &curves,
            groups(&[0.0, 2.0]),
            ReactionWeightedLimits::default(),
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn scaled_cubic_terms_keep_narrow_large_products_and_positive_subnormal_totals() {
    let cross = sigma(&[0.0, 1e-200], &[1e200, 1e200], EndfInterpolation::Linear);
    let flux = spectrum(&[0.0, 1e-200], &[1e100, 1e100], EndfInterpolation::Linear);
    let yield_ = curve(&[0.0, 1e-200], &[1e100, 1e100], EndfInterpolation::Linear);
    let answer = collapse(&cross, &flux, &[yield_], groups(&[0.0, 1e-200]));
    near(answer.cross_sections_m2[0], 1e200);
    near(answer.flux.values_m2_s()[0], 1e-100);
    near(answer.yields[0].values[0].unwrap(), 1e100);
    near(
        answer.cross_sections_m2[0]
            * answer.flux.values_m2_s()[0]
            * answer.yields[0].values[0].unwrap(),
        1e200,
    );

    let tiny = f64::from_bits(1);
    let cross = sigma(&[0.0, 1.0], &[tiny, tiny], EndfInterpolation::Linear);
    let flux = spectrum(&[0.0, 1.0], &[1.0, 1.0], EndfInterpolation::Linear);
    let yield_ = curve(&[0.0, 1.0], &[2.0, 2.0], EndfInterpolation::Linear);
    let answer = collapse(&cross, &flux, &[yield_], groups(&[0.0, 1.0]));
    assert_eq!(answer.cross_sections_m2, [tiny]);
    near(answer.yields[0].values[0].unwrap(), 2.0);
    assert_eq!(
        answer.cross_sections_m2[0] * answer.yields[0].values[0].unwrap(),
        f64::from_bits(2)
    );
}

#[test]
fn energy_dependent_neutron_yield_drives_the_original_fission_population_and_event_owner() {
    let cross = sigma(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear);
    let flux = spectrum(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear);
    let yield_ = curve(&[0.0, 2.0], &[0.0, 2.0], EndfInterpolation::Linear);
    let answer = collapse(&cross, &flux, &[yield_], groups(&[0.0, 2.0]));
    let model = DiffusionModel::new(
        DiffusionData {
            data_id: "manufactured-reaction-weighted-fission".into(),
            groups: Arc::clone(answer.flux.groups()),
            speeds_m_s: vec![1.0],
            cell_volumes_m3: vec![1.0],
            delayed_groups: vec![],
            reactions: vec![DiffusionReaction {
                data_id: "manufactured-fission-with-evaluated-yield".into(),
                cell: 0,
                incident_group: 0,
                macroscopic_cross_section_m_inv: answer.cross_sections_m2[0],
                kind: DiffusionReactionKind::Fission {
                    prompt_yields: vec![answer.yields[0].values[0].unwrap()],
                    delayed_yields: vec![],
                },
            }],
            links: vec![],
            leakage: vec![],
        },
        DiffusionLimits::default(),
    )
    .unwrap();
    let initial = DiffusionState::new(Arc::new(model), vec![2.0], vec![]).unwrap();
    let time = 0.7_f64;
    let (next, report) = initial
        .trial(
            time,
            DiffusionBudget {
                limits: DiffusionLimits::default(),
                population_tolerance: 1e-9,
            },
            || Ok(()),
        )
        .unwrap();
    // Frozen authored fuel density 1/m³ gives loss rate 4/3 and neutron
    // production 2/s per neutron: net growth 2/3, not the flux-weighted value.
    near(next.neutron_counts()[0], 2.0 * (2.0 * time / 3.0).exp());
    near(report.fissions, 4.0 * (2.0 * time / 3.0).exp_m1());
    near(report.prompt_neutrons, 6.0 * (2.0 * time / 3.0).exp_m1());
    assert!(report.neutron_balance_residual.abs() < 1e-12);
}

#[test]
fn yield_admission_and_evaluation_preserve_laws_and_refuse_invalid_payloads() {
    let yield_ = TabulatedYield::new(
        TabulatedYieldData {
            data_id: "manufactured-mixed-yield-law".into(),
            energies_j: vec![0.0, 1.0, 2.0],
            values: vec![2.0, 4.0, 8.0],
            interpolation_regions: vec![
                EndfInterpolationRegion {
                    end_point: 2,
                    interpolation: EndfInterpolation::Histogram,
                },
                EndfInterpolationRegion {
                    end_point: 3,
                    interpolation: EndfInterpolation::Linear,
                },
            ],
        },
        TabulatedYieldLimits::default(),
    )
    .unwrap();
    assert_eq!(yield_.value_at_energy(0.5).unwrap(), 2.0);
    assert_eq!(yield_.value_at_energy(1.0).unwrap(), 4.0);
    assert_eq!(yield_.value_at_energy(1.5).unwrap(), 6.0);
    assert_eq!(yield_.value_at_energy(2.0).unwrap(), 8.0);
    for energy in [-1.0, 3.0, f64::NAN, f64::INFINITY] {
        assert!(yield_.value_at_energy(energy).is_err());
    }
    let extreme = curve(&[0.0, 1.0], &[1e300, 0.0], EndfInterpolation::Linear);
    near(
        extreme.value_at_energy(1.0_f64.next_down()).unwrap(),
        1e300 * (1.0 - 1.0_f64.next_down()),
    );
    let valid = curve_data(&[0.0, 1.0], &[0.0, 3.0], EndfInterpolation::Linear);
    let mut bad = Vec::new();
    let mut entry = valid.clone();
    entry.data_id.clear();
    bad.push(entry);
    let mut entry = valid.clone();
    entry.values.pop();
    bad.push(entry);
    let mut entry = valid.clone();
    entry.interpolation_regions.clear();
    bad.push(entry);
    let mut entry = valid.clone();
    entry.interpolation_regions[0].end_point = 1;
    bad.push(entry);
    let mut entry = valid.clone();
    entry.interpolation_regions[0].end_point = 3;
    bad.push(entry);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut entry = valid.clone();
        entry.values[0] = value;
        bad.push(entry);
        let mut entry = valid.clone();
        entry.energies_j[0] = value;
        bad.push(entry);
    }
    bad.push(curve_data(
        &[1.0, 0.0],
        &[1.0, 1.0],
        EndfInterpolation::Linear,
    ));
    bad.push(curve_data(
        &[0.0, 0.0],
        &[1.0, 2.0],
        EndfInterpolation::Linear,
    ));
    bad.push(curve_data(
        &[0.0, 1.0, 1.0, 1.0, 2.0],
        &[1.0; 5],
        EndfInterpolation::Linear,
    ));
    for entry in bad {
        assert!(TabulatedYield::new(entry, TabulatedYieldLimits::default()).is_err());
    }
    for limits in [
        TabulatedYieldLimits {
            max_points: 1,
            max_interpolation_regions: 1,
        },
        TabulatedYieldLimits {
            max_points: 2,
            max_interpolation_regions: 0,
        },
    ] {
        assert!(matches!(
            TabulatedYield::new(valid.clone(), limits),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    assert_eq!(
        TabulatedYield::new_with_checkpoint(valid, TabulatedYieldLimits::default(), || Err(
            FissionError::Cancelled
        ))
        .unwrap_err(),
        EndfError::Physics(FissionError::Cancelled)
    );
}

#[test]
fn shared_work_bounds_domains_range_and_cancellation_refuse_without_partial_results() {
    let cross = sigma(&[0.0, 2.0], &[1.0, 3.0], EndfInterpolation::Linear);
    let flux = spectrum(&[0.0, 2.0], &[1.0, 2.0], EndfInterpolation::Linear);
    let curves = vec![curve(&[0.0, 2.0], &[1.0, 2.0], EndfInterpolation::Linear); 3];
    let grid = groups(&[0.0, 2.0]);
    let limits = ReactionWeightedLimits {
        collapse: EndfCollapseLimits {
            max_neutron_groups: 1,
            max_segment_visits: 5,
        },
        max_curves: 3,
        max_values: 3,
    };
    let expected =
        collapse_reaction_weighted(&cross, &flux, &curves, Arc::clone(&grid), limits, || Ok(()))
            .unwrap();
    // Exactly five segment visits: one shared flux, one shared denominator,
    // and one per yield. A batch does not repeat the first two integrations.
    for strict in [
        ReactionWeightedLimits {
            collapse: EndfCollapseLimits {
                max_segment_visits: 4,
                ..limits.collapse
            },
            ..limits
        },
        ReactionWeightedLimits {
            collapse: EndfCollapseLimits {
                max_neutron_groups: 0,
                ..limits.collapse
            },
            ..limits
        },
        ReactionWeightedLimits {
            max_curves: 2,
            ..limits
        },
        ReactionWeightedLimits {
            max_values: 2,
            ..limits
        },
    ] {
        assert!(matches!(
            collapse_reaction_weighted(
                &cross,
                &flux,
                &curves,
                Arc::clone(&grid),
                strict,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let mut calls = 0;
    collapse_reaction_weighted(&cross, &flux, &curves, Arc::clone(&grid), limits, || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert_eq!(
            collapse_reaction_weighted(&cross, &flux, &curves, Arc::clone(&grid), limits, || {
                seen += 1;
                if seen == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .unwrap_err(),
            EndfError::Physics(FissionError::Cancelled)
        );
    }
    let short_curve = curve(&[0.0, 1.0], &[1.0, 1.0], EndfInterpolation::Linear);
    assert!(
        collapse_reaction_weighted(
            &cross,
            &flux,
            &[short_curve],
            Arc::clone(&grid),
            limits,
            || Ok(())
        )
        .is_err()
    );
    let short_cross = sigma(&[0.0, 1.0], &[1.0, 1.0], EndfInterpolation::Linear);
    assert!(
        collapse_reaction_weighted(
            &short_cross,
            &flux,
            &curves,
            Arc::clone(&grid),
            limits,
            || Ok(())
        )
        .is_err()
    );
    let short_flux = spectrum(&[0.0, 1.0], &[1.0, 1.0], EndfInterpolation::Linear);
    assert!(
        collapse_reaction_weighted(
            &cross,
            &short_flux,
            &curves,
            Arc::clone(&grid),
            limits,
            || Ok(())
        )
        .is_err()
    );
    assert!(
        collapse_reaction_weighted(&cross, &flux, &[], Arc::clone(&grid), limits, || Ok(()))
            .is_err()
    );
    for value in [1e200, 1e-200] {
        let extreme_cross = sigma(&[0.0, 1.0], &[value, value], EndfInterpolation::Linear);
        let unit_flux = spectrum(&[0.0, 1.0], &[1.0, 1.0], EndfInterpolation::Linear);
        let extreme_yield = curve(&[0.0, 1.0], &[value, value], EndfInterpolation::Linear);
        assert!(matches!(
            collapse_reaction_weighted(
                &extreme_cross,
                &unit_flux,
                &[extreme_yield],
                groups(&[0.0, 1.0]),
                limits,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Arithmetic))
        ));
    }
    assert_eq!(
        collapse_reaction_weighted(&cross, &flux, &curves, grid, limits, || Ok(())).unwrap(),
        expected
    );
}
