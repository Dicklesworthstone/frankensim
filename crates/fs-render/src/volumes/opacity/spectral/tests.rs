use super::*;
use crate::volumes::render_density_transport_box;

fn limits(cells: usize, resolution: [usize; 2]) -> SpectralAbsorptionImageSettings {
    SpectralAbsorptionImageSettings {
        image: AbsorptionImageSettings {
            resolution,
            max_pixels: 128,
            max_grid_cells: cells,
            max_cell_visits: 10_000,
        },
        max_spectral_evaluations: 10_000,
    }
}

fn spectrum<'a>(wavelengths: &'a [f64], coefficients: &'a [f64]) -> DensityAbsorptionSpectrum<'a> {
    DensityAbsorptionSpectrum {
        wavelengths_nm: wavelengths,
        mass_absorption_m2_kg: coefficients,
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn spectral_slabs_preserve_wavelengths_units_and_native_column_attenuation() {
    let grid = VolumeGrid::new([1, 1, 2], &[2.0, 2.0], [0.0; 3], [1.0, 1.0, 0.5]);
    let calibration = DensityAbsorptionSpectrum {
        density_unit_kg_m3: 3.0,
        length_unit_m: 0.2,
        ..spectrum(&[400.0, 550.0, 900.0], &[0.0, 0.25, 1.5])
    };
    let image = render_density_spectral_box(&grid, limits(2, [2, 2]), calibration).unwrap();
    assert_eq!(image.resolution(), [2, 2]);
    assert_eq!(image.wavelengths_nm(), &[400.0, 550.0, 900.0]);
    for pixel in image.samples().chunks_exact(3) {
        for (sample, tau) in pixel.iter().zip([0.0_f64, 0.3, 1.8]) {
            close(sample.transmittance(), (-tau).exp());
            close(sample.opacity(), -(-tau).exp_m1());
            close(sample.optical_depth(), tau);
            close(sample.attenuate(5.0).unwrap(), 5.0 * (-tau).exp());
        }
    }
    let second = VolumeGrid::new([1, 1, 2], &[1.0, 1.0], [0.0; 3], [1.0, 1.0, 0.25]);
    let rescaled = DensityAbsorptionSpectrum {
        density_unit_kg_m3: 6.0,
        length_unit_m: 0.4,
        ..calibration
    };
    let same = render_density_spectral_box(&second, limits(2, [2, 2]), rescaled).unwrap();
    assert_eq!(image, same);
}

#[test]
fn wavelength_dependent_subpixel_coverage_survives_physical_grid_refinement() {
    let density = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let grid = VolumeGrid::new([4, 2, 1], &density, [0.0; 3], [1.0; 3]);
    let coefficients = [1000.0, 0.0, core::f64::consts::LN_2];
    let calibration = spectrum(&[400.0, 550.0, 900.0], &coefficients);
    let image = render_density_spectral_box(&grid, limits(8, [1, 1]), calibration).unwrap();
    for (sample, expected) in image.samples().iter().zip([0.75_f64, 1.0, 0.875]) {
        close(sample.transmittance(), expected);
        close(sample.opacity(), 1.0 - expected);
        close(sample.optical_depth(), -expected.ln());
    }
    let fine: Vec<_> = (0..4)
        .flat_map(|y| (0..8).map(move |x| density[x / 2 + 4 * (y / 2)]))
        .collect();
    let refined = VolumeGrid::new([8, 4, 1], &fine, [1e12, -1e12, 0.0], [0.5, 0.5, 1.0]);
    let second = render_density_spectral_box(&refined, limits(32, [1, 1]), calibration).unwrap();
    assert_eq!(image, second);
}

#[test]
fn mixed_nondivisible_pixels_match_independent_area_integrals_and_monochromatic_owner() {
    let density = [1.0, 0.2, 0.0, 3.0, 2.0, 0.5];
    let grid = VolumeGrid::new([3, 2, 1], &density, [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[450.0, 650.0], &[0.5, 2.0]);
    let request = limits(6, [2, 3]);
    let image = render_density_spectral_box(&grid, request, calibration).unwrap();
    for py in 0..3 {
        for px in 0..2 {
            let pixel = image.pixel(px, py).unwrap();
            for (lane, &coefficient) in calibration.mass_absorption_m2_kg.iter().enumerate() {
                let mut expected = 0.0;
                let mut mean_density = 0.0;
                // Independent physical XY rectangle intersections; pixel area is 1.
                for iy in 0..2 {
                    for ix in 0..3 {
                        let wx = ((px + 1) as f64 * 1.5).min((ix + 1) as f64)
                            - (px as f64 * 1.5).max(ix as f64);
                        let wy = ((py + 1) as f64 * (2.0 / 3.0)).min((iy + 1) as f64)
                            - (py as f64 * (2.0 / 3.0)).max(iy as f64);
                        let area = wx.max(0.0) * wy.max(0.0);
                        expected += area * (-coefficient * density[ix + 3 * iy]).exp();
                        mean_density += area * density[ix + 3 * iy];
                    }
                }
                close(pixel[lane].transmittance(), expected);
                assert!(pixel[lane].transmittance() > (-coefficient * mean_density).exp());
            }
        }
    }
    for (lane, &coefficient) in calibration.mass_absorption_m2_kg.iter().enumerate() {
        let mono = render_density_transport_box(
            &grid,
            request.image,
            DensityAbsorption {
                density_unit_kg_m3: 1.0,
                length_unit_m: 1.0,
                mass_absorption_m2_kg: coefficient,
            },
        )
        .unwrap();
        for (index, value) in mono.iter().enumerate() {
            let sample = image.samples()[2 * index + lane];
            close(sample.transmittance(), value.transmittance());
            close(sample.opacity(), value.opacity());
            close(sample.optical_depth(), value.optical_depth());
        }
    }
}

#[test]
fn scaled_column_mass_preserves_finite_depth_after_the_mass_sum_would_overflow() {
    let grid = VolumeGrid::new([1, 1, 2], &[f64::MAX; 2], [0.0; 3], [1.0; 3]);
    let image = render_density_spectral_box(
        &grid,
        limits(2, [1, 1]),
        spectrum(&[500.0, 600.0], &[1e-308, 0.0]),
    )
    .unwrap();
    let tau = 2.0 * (f64::MAX * 1e-308);
    close(image.samples()[0].optical_depth(), tau);
    close(image.samples()[0].transmittance(), (-tau).exp());
    assert_eq!(image.samples()[1].transmittance(), 1.0);
    assert_eq!(image.samples()[1].opacity(), 0.0);

    let small = VolumeGrid::new([1, 1, 2], &[1e-320; 2], [0.0; 3], [1.0, 1.0, 1e300]);
    let thin =
        render_density_spectral_box(&small, limits(2, [1, 1]), spectrum(&[500.0], &[1.0])).unwrap();
    close(thin.samples()[0].opacity(), (1e-320 * 1e300) * 2.0);
}

#[test]
fn thin_opacity_and_bright_backgrounds_survive_in_separate_spectral_lanes() {
    let grid = VolumeGrid::new([1; 3], &[1.0], [0.0; 3], [1.0; 3]);
    let image = render_density_spectral_box(
        &grid,
        limits(1, [1, 1]),
        spectrum(&[450.0, 650.0], &[1e-20, 800.0]),
    )
    .unwrap();
    let thin = image.samples()[0];
    assert_eq!(thin.transmittance(), 1.0);
    close(thin.opacity(), 1e-20);
    close(thin.optical_depth(), 1e-20);
    let opaque = image.samples()[1];
    assert_eq!(opaque.transmittance(), 0.0);
    close(opaque.optical_depth(), 800.0);
    close(
        opaque.attenuate(1e300).unwrap(),
        (1e300_f64.ln() - 800.0).exp(),
    );
    assert_eq!(opaque.attenuate(0.0).unwrap(), 0.0);
    assert_eq!(
        thin.attenuate(-1.0),
        Err(DvrError::InvalidBackgroundRadiance)
    );
}

#[test]
fn direct_transmission_matches_full_statistics_and_preserves_pixel_band_layout() {
    let grid = VolumeGrid::new([2, 2, 1], &[1.0, 2.0, 3.0, 4.0], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[400.0, 550.0, 700.0], &[0.0, 0.5, 2.0]);
    for resolution in [[2, 2], [3, 1], [1, 3]] {
        let request = limits(4, resolution);
        let full = render_density_spectral_box(&grid, request, calibration).unwrap();
        let direct =
            render_density_spectral_transmittance_box(&grid, request, calibration).unwrap();
        assert_eq!(direct.len(), resolution[0] * resolution[1] * 3);
        for (flat, sample) in direct.iter().zip(full.samples()) {
            assert_eq!(*flat, sample.transmittance());
        }
        assert_eq!(full.pixel(resolution[0], 0), None);
        assert_eq!(full.pixel(0, resolution[1]), None);
        if resolution == [2, 2] {
            for (pixel, density) in direct.chunks_exact(3).zip([1.0_f64, 2.0, 3.0, 4.0]) {
                for (&value, coefficient) in pixel.iter().zip([0.0, 0.5, 2.0]) {
                    close(value, (-density * coefficient).exp());
                }
            }
        }
    }
}

#[test]
fn spectral_work_is_admitted_before_field_scan_and_counts_all_original_overlaps() {
    let grid = VolumeGrid::new([3, 2, 2], &[1.0; 12], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[400.0, 500.0, 600.0], &[0.0, 0.5, 2.0]);
    let exact = SpectralAbsorptionImageSettings {
        image: AbsorptionImageSettings {
            max_cell_visits: 32,
            ..limits(12, [2, 3]).image
        },
        max_spectral_evaluations: 48,
    };
    assert_eq!(
        render_density_spectral_box(&grid, exact, calibration)
            .unwrap()
            .samples()
            .len(),
        18
    );
    let mut invalid = VolumeGrid::new([3, 2, 2], &[f64::NAN; 12], [0.0; 3], [1.0; 3]);
    let short = SpectralAbsorptionImageSettings {
        max_spectral_evaluations: 47,
        ..exact
    };
    assert_eq!(
        render_density_spectral_box(&invalid, short, calibration),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 48,
            limit: 47
        }),
    );
    assert_eq!(
        render_density_spectral_transmittance_box(&invalid, short, calibration),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 48,
            limit: 47
        }),
    );
    assert_eq!(
        render_density_spectral_box(
            &grid,
            SpectralAbsorptionImageSettings {
                image: AbsorptionImageSettings {
                    max_cell_visits: 31,
                    ..exact.image
                },
                ..exact
            },
            calibration
        ),
        Err(DvrError::CellVisitBudgetExceeded {
            requested: 32,
            limit: 31
        }),
    );
    invalid.data = &[-1.0; 12];
    assert_eq!(
        render_density_spectral_box(&invalid, exact, calibration),
        Err(DvrError::NegativeGridAbsorption { index: 0 }),
    );
}

#[test]
fn invalid_spectra_refuse_and_the_maximum_packet_is_bounded() {
    let grid = VolumeGrid::new([1; 3], &[1.0], [0.0; 3], [1.0; 3]);
    for (wavelengths, coefficients) in [
        (&[][..], &[][..]),
        (&[500.0][..], &[][..]),
        (&[500.0, 400.0][..], &[1.0, 1.0][..]),
        (&[500.0, 500.0][..], &[1.0, 1.0][..]),
        (&[0.0][..], &[1.0][..]),
        (&[f64::NAN][..], &[1.0][..]),
        (&[500.0][..], &[-1.0][..]),
        (&[500.0][..], &[f64::INFINITY][..]),
        (&[500.0; 17][..], &[1.0; 17][..]),
    ] {
        assert_eq!(
            render_density_spectral_box(
                &grid,
                limits(1, [1, 1]),
                spectrum(wavelengths, coefficients)
            ),
            Err(DvrError::InvalidAbsorptionSpectrum),
        );
    }
    let wavelengths = std::array::from_fn::<_, MAX_OPACITY_WAVELENGTHS, _>(|i| 400.0 + i as f64);
    let calibration = spectrum(&wavelengths, &[1.0; MAX_OPACITY_WAVELENGTHS]);
    let image = render_density_spectral_box(&grid, limits(1, [1, 1]), calibration).unwrap();
    assert_eq!(image.wavelengths_nm(), wavelengths);
    assert_eq!(image.samples().len(), MAX_OPACITY_WAVELENGTHS);
    assert_eq!(
        render_density_spectral_box(
            &grid,
            limits(1, [1, 1]),
            DensityAbsorptionSpectrum {
                density_unit_kg_m3: 0.0,
                ..calibration
            }
        ),
        Err(DvrError::InvalidDensityAbsorption),
    );
}
