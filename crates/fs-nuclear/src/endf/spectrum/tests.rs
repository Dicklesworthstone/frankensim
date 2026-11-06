//! Manufactured source records and analytic spectra; no evaluated coefficients.

use super::*;
use crate::cross_sections::NeutronReaction;
use crate::endf::{EV_J, EndfReadLimits, EndfSelection};
use crate::irradiation::NuclideId;

const BARN: f64 = 1e-28;

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("fixture field is too long");
}

fn record(fields: [String; 6], mf: i32, mt: i32) -> String {
    assert!(fields.iter().all(|s| s.len() == 11));
    format!("{}9237{mf:2}{mt:3}{:5}", fields.concat(), 1)
}

fn control(reals: [f64; 2], integers: [i32; 4], mf: i32, mt: i32) -> String {
    record(
        [
            real(reals[0]),
            real(reals[1]),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mf,
        mt,
    )
}

fn curve(
    temperature: f64,
    reaction: NeutronReaction,
    points_ev_barns: &[(f64, f64)],
    regions: &[(usize, EndfInterpolation)],
) -> PointwiseCrossSection {
    let mt = match reaction {
        NeutronReaction::Fission => 18,
        NeutronReaction::Capture => 102,
    };
    let mut lines = vec![
        control([92235.0, 233.0], [2, 1, 0, 7], 1, 451),
        control([0.0, 0.0], [0, 0, 0, 6], 1, 451),
        control([1.0, 2.0e7], [3, 0, 10, 8], 1, 451),
        control([temperature, 0.0], [1, 0, 5, 2], 1, 451),
    ];
    for description in [
        "MANUFACTURED SPECTRUM TEST",
        "NO EVALUATED COEFFICIENTS",
        "ANALYTIC POINTWISE LAW",
        "NO PROCESSING CLAIM",
        "NO NUCLEAR-DATA ACCURACY CLAIM",
    ] {
        lines.push(format!("{description:<66}9237{:2}{:3}{:5}", 1, 451, 1));
    }
    lines.push(control([0.0; 2], [1, 451, 11, 0], 1, 451));
    let count = 2 + regions.len().div_ceil(3) + points_ev_barns.len().div_ceil(3);
    lines.push(control([0.0; 2], [3, mt, count as i32, 0], 1, 451));
    lines.push(control([0.0; 2], [0; 4], 1, 0));
    lines.push(control([0.0; 2], [0; 4], 0, 0));
    lines.push(control([92235.0, 233.0], [0; 4], 3, mt));
    lines.push(control(
        [0.0; 2],
        [0, 0, regions.len() as i32, points_ev_barns.len() as i32],
        3,
        mt,
    ));
    for chunk in regions.chunks(3) {
        let mut fields = std::array::from_fn(|_| format!("{:11}", 0));
        for (index, &(end, interpolation)) in chunk.iter().enumerate() {
            let law = match interpolation {
                EndfInterpolation::Histogram => 1,
                EndfInterpolation::Linear => 2,
            };
            fields[index * 2] = format!("{end:11}");
            fields[index * 2 + 1] = format!("{law:11}");
        }
        lines.push(record(fields, 3, mt));
    }
    for chunk in points_ev_barns.chunks(3) {
        let mut fields = std::array::from_fn(|_| real(0.0));
        for (index, &(energy, sigma)) in chunk.iter().enumerate() {
            fields[index * 2] = real(energy);
            fields[index * 2 + 1] = real(sigma);
        }
        lines.push(record(fields, 3, mt));
    }
    lines.push(control([0.0; 2], [0; 4], 3, 0));
    PointwiseCrossSection::parse(
        format!("manufactured-spectral-xs-{temperature}"),
        EndfSelection {
            material_number: 9237,
            nuclide: NuclideId::new(92, 235, 0).unwrap(),
            reaction,
            temperature_k: temperature,
        },
        &lines.join("\n"),
        EndfReadLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn spectrum_data(
    points_ev: &[(f64, f64)],
    regions: &[(usize, EndfInterpolation)],
    scale: f64,
) -> FluxSpectrumData {
    FluxSpectrumData {
        data_id: format!("manufactured-absolute-spectrum-{scale}"),
        energies_j: points_ev.iter().map(|p| p.0 * EV_J).collect(),
        flux_density_m2_s_j: points_ev.iter().map(|p| p.1 * scale / EV_J).collect(),
        interpolation_regions: regions
            .iter()
            .map(|&(end_point, interpolation)| EndfInterpolationRegion {
                end_point,
                interpolation,
            })
            .collect(),
    }
}

fn spectrum(
    points_ev: &[(f64, f64)],
    regions: &[(usize, EndfInterpolation)],
    scale: f64,
) -> FluxSpectrum {
    FluxSpectrum::new(
        spectrum_data(points_ev, regions, scale),
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

fn groups(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "manufactured-spectral-groups".into(),
            edges_ev.iter().map(|e| e * EV_J).collect(),
        )
        .unwrap(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs() * 5e-13,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn linear_weighting_changes_the_coefficient_and_preserves_the_pointwise_integral() {
    let source = curve(
        300.0,
        NeutronReaction::Fission,
        &[(0.0, 0.0), (4.0, 4.0)],
        &[(2, EndfInterpolation::Linear)],
    );
    let flux = spectrum(
        &[(0.0, 0.0), (4.0, 4.0)],
        &[(2, EndfInterpolation::Linear)],
        1e26,
    );
    let collapsed = source
        .collapse_with_spectrum(
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    // With x=E/eV, σ=x barn and φ=x*1e26/eV: ∫φ=8e26,
    // ∫σφ=(64/3)*1e26 barn/s, and σ_g=(8/3) barn.
    close(collapsed.flux.values_m2_s()[0], 8e26);
    close(collapsed.cross_sections_m2[0], (8.0 / 3.0) * BARN);
    close(
        collapsed.cross_sections_m2[0] * collapsed.flux.values_m2_s()[0],
        (64.0 / 3.0) * 1e26 * BARN,
    );
    let flat = source
        .collapse_flat_in_energy(&groups(&[0.0, 4.0]), EndfCollapseLimits::default(), || {
            Ok(())
        })
        .unwrap();
    close(flat[0], 2.0 * BARN);
    assert!(collapsed.cross_sections_m2[0] > flat[0]);
    let bins = source
        .collapse_with_spectrum(
            &flux,
            groups(&[0.0, 2.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    for (index, expected) in [4.0 / 3.0, 28.0 / 9.0].into_iter().enumerate() {
        close(bins.cross_sections_m2[index], expected * BARN);
    }
    close(bins.flux.values_m2_s()[0], 2e26);
    close(bins.flux.values_m2_s()[1], 6e26);
}

#[test]
fn merged_mesh_integration_resolves_thresholds_mixed_laws_and_both_discontinuities() {
    let source = curve(
        300.0,
        NeutronReaction::Fission,
        &[(0.0, 0.0), (1.0, 0.0), (2.0, 4.0), (2.0, 10.0), (4.0, 0.0)],
        &[
            (2, EndfInterpolation::Histogram),
            (5, EndfInterpolation::Linear),
        ],
    );
    let flux = spectrum(
        &[
            (0.0, 1.0),
            (0.5, 2.0),
            (1.5, 2.0),
            (1.5, 4.0),
            (3.0, 1.0),
            (4.0, 1.0),
        ],
        &[
            (2, EndfInterpolation::Linear),
            (3, EndfInterpolation::Histogram),
            (5, EndfInterpolation::Linear),
            (6, EndfInterpolation::Histogram),
        ],
        1e26,
    );
    close(flux.density_at_energy(1.5 * EV_J).unwrap(), 4e26 / EV_J);
    let collapsed = source
        .collapse_with_spectrum(
            &flux,
            groups(&[0.0, 1.5, 2.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    for (index, (phi, sigma)) in [
        (11.0 / 4.0, 4.0 / 11.0),
        (7.0 / 4.0, 62.0 / 21.0),
        (3.0, 55.0 / 9.0),
    ]
    .into_iter()
    .enumerate()
    {
        close(collapsed.flux.values_m2_s()[index], phi * 1e26);
        close(collapsed.cross_sections_m2[index], sigma * BARN);
    }
    let whole = source
        .collapse_with_spectrum(
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    close(whole.cross_sections_m2[0], (49.0 / 15.0) * BARN);
    close(whole.flux.values_m2_s()[0], 7.5e26);
    let threshold = source
        .collapse_with_spectrum(
            &flux,
            groups(&[0.0, 0.5]),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(threshold.cross_sections_m2[0], 0.0);
    close(threshold.flux.values_m2_s()[0], 0.75e26);
}

#[test]
fn coefficients_are_invariant_to_absolute_flux_scale_and_flat_spectra_match_flat_collapse() {
    let source = curve(
        300.0,
        NeutronReaction::Capture,
        &[(0.0, 0.0), (4.0, 4.0)],
        &[(2, EndfInterpolation::Linear)],
    );
    for scale in [1e-270, 1e-30, 1.0, 1e30, 1e270] {
        let flux = spectrum(
            &[(0.0, 0.0), (4.0, 4.0)],
            &[(2, EndfInterpolation::Linear)],
            scale,
        );
        let collapsed = source
            .collapse_with_spectrum(
                &flux,
                groups(&[0.0, 4.0]),
                EndfCollapseLimits::default(),
                || Ok(()),
            )
            .unwrap();
        close(collapsed.cross_sections_m2[0], (8.0 / 3.0) * BARN);
        close(collapsed.flux.values_m2_s()[0], 8.0 * scale);
    }
    let flux = spectrum(
        &[(0.0, 3.0), (1.0, 3.0), (4.0, 3.0)],
        &[(3, EndfInterpolation::Histogram)],
        5e25,
    );
    let grid = groups(&[0.0, 2.0, 4.0]);
    let weighted = source
        .collapse_with_spectrum(
            &flux,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    let flat = source
        .collapse_flat_in_energy(&grid, EndfCollapseLimits::default(), || Ok(()))
        .unwrap();
    for (&a, b) in weighted.cross_sections_m2.iter().zip(flat) {
        close(a, b);
    }
}

#[test]
fn scaled_positive_integrals_avoid_false_intermediate_overflow_and_keep_subnormal_totals() {
    let source = curve(
        300.0,
        NeutronReaction::Capture,
        &[(0.0, 1e308), (1e-261, 1e308)],
        &[(2, EndfInterpolation::Linear)],
    );
    let width = source.energies_j()[1];
    let grid =
        Arc::new(NeutronGroups::new("narrow-positive-domain".into(), vec![0.0, width]).unwrap());
    let flux = FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-large-density".into(),
            energies_j: vec![0.0, width],
            flux_density_m2_s_j: vec![1e300; 2],
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap();
    assert!((source.cross_sections_m2()[0] * 1e300).is_infinite());
    let collapsed = source
        .collapse_with_spectrum(&flux, grid, EndfCollapseLimits::default(), || Ok(()))
        .unwrap();
    close(
        collapsed.cross_sections_m2[0],
        source.cross_sections_m2()[0],
    );
    close(collapsed.flux.values_m2_s()[0], width * 1e300);

    let unit = curve(
        300.0,
        NeutronReaction::Capture,
        &[(0.0, 1e28), (1e19, 1e28)],
        &[(2, EndfInterpolation::Linear)],
    );
    let grid = Arc::new(NeutronGroups::new("unit-width-domain".into(), vec![0.0, 1.0]).unwrap());
    let data = FluxSpectrumData {
        data_id: "manufactured-subnormal-density".into(),
        energies_j: vec![0.0, 1.0],
        flux_density_m2_s_j: vec![f64::from_bits(1); 2],
        interpolation_regions: vec![EndfInterpolationRegion {
            end_point: 2,
            interpolation: EndfInterpolation::Linear,
        }],
    };
    let tiny = FluxSpectrum::new(data, FluxSpectrumLimits::default()).unwrap();
    let collapsed = unit
        .collapse_with_spectrum(
            &tiny,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(collapsed.flux.values_m2_s()[0], f64::from_bits(1));
    close(collapsed.cross_sections_m2[0], 1.0);
    let unrepresentable = FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-overflow".into(),
            energies_j: vec![0.0, 1.0],
            flux_density_m2_s_j: vec![1e308; 2],
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap();
    let huge_sigma = curve(
        300.0,
        NeutronReaction::Capture,
        &[(0.0, 1e29), (1e19, 1e29)],
        &[(2, EndfInterpolation::Linear)],
    );
    assert_eq!(
        huge_sigma.collapse_with_spectrum(
            &unrepresentable,
            grid,
            EndfCollapseLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Arithmetic))
    );
}

#[test]
fn spectrum_admission_and_undefined_zero_flux_groups_are_refused() {
    let valid = spectrum_data(
        &[(0.0, 1.0), (1.0, 0.0), (4.0, 0.0)],
        &[(3, EndfInterpolation::Histogram)],
        1e26,
    );
    let mut bad_inputs = Vec::new();
    let mut data = valid.clone();
    data.data_id.clear();
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.flux_density_m2_s_j.pop();
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.flux_density_m2_s_j[0] = f64::NAN;
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.flux_density_m2_s_j[0] = -1.0;
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.energies_j = vec![1.0, 0.0, 4.0];
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.energies_j = vec![1.0; 3];
    bad_inputs.push(data);
    let mut data = valid.clone();
    data.interpolation_regions[0].end_point = 2;
    bad_inputs.push(data);
    for data in bad_inputs {
        assert!(FluxSpectrum::new(data, FluxSpectrumLimits::default()).is_err());
    }
    for limits in [
        FluxSpectrumLimits {
            max_points: 2,
            ..FluxSpectrumLimits::default()
        },
        FluxSpectrumLimits {
            max_interpolation_regions: 0,
            ..FluxSpectrumLimits::default()
        },
    ] {
        assert!(matches!(
            FluxSpectrum::new(valid.clone(), limits),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let flux = FluxSpectrum::new(valid, FluxSpectrumLimits::default()).unwrap();
    let source = curve(
        300.0,
        NeutronReaction::Capture,
        &[(0.0, 1.0), (4.0, 1.0)],
        &[(2, EndfInterpolation::Linear)],
    );
    assert!(matches!(
        source.collapse_with_spectrum(
            &flux,
            groups(&[0.0, 1.0, 4.0]),
            EndfCollapseLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::InvalidInput(_)))
    ));
    for points in [[(1.0, 1.0), (4.0, 1.0)], [(5.0, 1.0), (6.0, 1.0)]] {
        let outside = spectrum(&points, &[(2, EndfInterpolation::Linear)], 1e26);
        assert!(
            source
                .collapse_with_spectrum(
                    &outside,
                    groups(&[0.0, 4.0]),
                    EndfCollapseLimits::default(),
                    || Ok(())
                )
                .is_err()
        );
    }
    assert!(flux.density_at_energy(-1.0).is_err());
    assert!(flux.density_at_energy(5.0 * EV_J).is_err());
}

#[test]
fn shared_work_limits_and_early_midway_late_cancellation_return_no_candidate() {
    let curves = [300.0, 600.0].map(|t| {
        curve(
            t,
            NeutronReaction::Fission,
            &[(0.0, 2.0), (4.0, 2.0)],
            &[(2, EndfInterpolation::Linear)],
        )
    });
    let flux = spectrum(
        &[(0.0, 1.0), (4.0, 1.0)],
        &[(2, EndfInterpolation::Linear)],
        1e26,
    );
    let limited = |visits| EndfCollapseLimits {
        max_segment_visits: visits,
        ..EndfCollapseLimits::default()
    };
    assert!(matches!(
        curves[0].collapse_with_spectrum(&flux, groups(&[0.0, 4.0]), limited(1), || Ok(())),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    // Exactly one spectrum interval plus one merged interval per temperature.
    assert!(
        collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &curves,
            &flux,
            groups(&[0.0, 4.0]),
            limited(3),
            CrossSectionLimits::default(),
            || Ok(())
        )
        .is_ok()
    );
    assert!(matches!(
        collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &curves,
            &flux,
            groups(&[0.0, 4.0]),
            limited(2),
            CrossSectionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    assert!(matches!(
        collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &curves,
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            CrossSectionLimits {
                max_values: 1,
                ..CrossSectionLimits::default()
            },
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    let mut calls = 0;
    let baseline = collapse_temperature_table(
        "identified-shared-spectrum".into(),
        &curves,
        &flux,
        groups(&[0.0, 4.0]),
        EndfCollapseLimits::default(),
        CrossSectionLimits::default(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut call = 0;
        let result = collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &curves,
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || {
                call += 1;
                if call == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result, Err(EndfError::Physics(FissionError::Cancelled)));
    }
    assert_eq!(
        collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &curves,
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(())
        )
        .unwrap(),
        baseline
    );
    assert!(
        collapse_temperature_table(
            "identified-shared-spectrum".into(),
            &[curves[1].clone(), curves[0].clone()],
            &flux,
            groups(&[0.0, 4.0]),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(
        FluxSpectrum::new_with_checkpoint(
            spectrum_data(
                &[(0.0, 1.0), (4.0, 1.0)],
                &[(2, EndfInterpolation::Linear)],
                1e26
            ),
            FluxSpectrumLimits::default(),
            || Err(FissionError::Cancelled)
        ),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
}

#[test]
fn weighted_temperature_rows_with_different_meshes_drive_the_original_burnup_integral() {
    use crate::EnergyChannels;
    use crate::burnup::{
        BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
        BurnupProduct, BurnupRate, BurnupReaction,
    };
    use fs_qty::Frequency;
    use fs_thermochem::decay::DecayLawV1;

    let curves = [
        curve(
            300.0,
            NeutronReaction::Fission,
            &[(0.0, 0.0), (4.0, 4.0)],
            &[(2, EndfInterpolation::Linear)],
        ),
        curve(
            600.0,
            NeutronReaction::Fission,
            &[(0.0, 0.0), (1.0, 2.0), (3.0, 6.0), (4.0, 8.0)],
            &[(4, EndfInterpolation::Linear)],
        ),
    ];
    let spectrum = spectrum(
        &[(0.0, 0.0), (4.0, 4.0)],
        &[(2, EndfInterpolation::Linear)],
        1e26,
    );
    let grid = groups(&[0.0, 1.0, 4.0]);
    let (table, flux) = collapse_temperature_table(
        "manufactured-source-spectrum-assembly".into(),
        &curves,
        &spectrum,
        Arc::clone(&grid),
        EndfCollapseLimits::default(),
        CrossSectionLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let separately = curves[1]
        .collapse_with_spectrum(
            &spectrum,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(flux, separately.flux);
    // At 450 K, σ(E)=1.5*(E/eV) barn. With φ=(E/eV)*1e26/eV,
    // the complete pointwise integral is (1.5*64/3)*1e26 barn = 0.32/s.
    close(
        table.rate_constant(&flux, 450.0, 2, || Ok(())).unwrap(),
        0.32,
    );
    let limits = BurnupLimits::default();
    let stable = DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-spectrum-to-burnup".into(),
                groups: grid,
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides: vec![
                    BurnupNuclide {
                        id: curves[0].metadata().nuclide,
                        law: stable,
                    },
                    BurnupNuclide {
                        id: NuclideId::new(40, 95, 0).unwrap(),
                        law: stable,
                    },
                ],
                reactions: vec![BurnupReaction {
                    data_id: "explicit-fragments-and-prompt-energy".into(),
                    parent: 0,
                    rate: BurnupRate::TabulatedFission {
                        cross_sections: Arc::new(table),
                    },
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: 2.0,
                    }],
                    energy_per_event: EnergyChannels {
                        deposited_j: 1e-12,
                        radiation_j: 0.0,
                        escaped_j: 0.0,
                    },
                }],
            },
            limits,
        )
        .unwrap(),
    );
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0]).unwrap();
    let report = inventory
        .advance(
            2.0,
            &flux,
            450.0,
            BurnupBudget {
                limits,
                nuclei_tolerance: 1e-9,
                energy_tolerance_j: 1e-18,
            },
            || Ok(()),
        )
        .unwrap();
    let fissions = 100.0 * -(-0.64_f64).exp_m1();
    close(report.fissions, fissions);
    close(inventory.populations()[0], 100.0 * (-0.64_f64).exp());
    close(inventory.populations()[1], 2.0 * fissions);
    close(report.released.deposited_j, fissions * 1e-12);
}
