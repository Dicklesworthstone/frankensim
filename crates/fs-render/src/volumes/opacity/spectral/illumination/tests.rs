use super::*;
use crate::{
    spectral::{cie_x, cie_y, cie_z},
    volumes::{AbsorptionImageSettings, render_density_spectral_box},
};

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
        (actual - expected).abs() <= expected.abs().max(1e-300) * 5e-12,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn original_column_lighting_preserves_opacity_correlation_in_each_wavelength() {
    let grid = VolumeGrid::new([2, 1, 1], &[0.0, 1.0], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[450.0, 650.0], &[core::f64::consts::LN_2; 2]);
    let backgrounds = [2.0, 10.0, 20.0, 4.0];
    let image = render_density_spectral_radiance_box(
        &grid,
        limits(2, [1, 1]),
        calibration,
        SpectralIllumination::Columns(&backgrounds),
    )
    .unwrap();
    close(image[0], 6.0);
    close(image[1], 6.0);
    // Averaging source and transmission separately would give 8.25 and 5.25.
    assert!((image[0] - 8.25).abs() > 2.0);
    assert!((image[1] - 5.25).abs() > 0.7);
    let fine = render_density_spectral_radiance_box(
        &grid,
        limits(2, [2, 1]),
        calibration,
        SpectralIllumination::Columns(&backgrounds),
    )
    .unwrap();
    for lane in 0..2 {
        close(image[lane], (fine[lane] + fine[2 + lane]) * 0.5);
    }
}

#[test]
fn indexed_layers_nondivisible_footprints_and_xyz_match_independent_integrals() {
    let density = [0.5, 1.0, 2.0, 3.0, 0.0, 0.25, 2.0, 0.0, 1.0, 0.5, 1.5, 3.0];
    let materials = [0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1];
    let coefficients = [0.25, 1.0, 0.5, 0.75];
    let wavelengths = [450.0, 650.0];
    let backgrounds = [
        1.0, 5.0, 11.0, 2.0, 3.0, 9.0, 7.0, 4.0, 8.0, 12.0, 6.0, 10.0,
    ];
    let weights = [25.0, 100.0];
    let grid = VolumeGrid::new([3, 2, 2], &density, [0.0; 3], [1.0, 1.0, 0.5]);
    let calibration = IndexedDensityAbsorptionSpectrum {
        wavelengths_nm: &wavelengths,
        mass_absorption_m2_kg: &coefficients,
        cell_materials: &materials,
        density_unit_kg_m3: 2.0,
        length_unit_m: 0.5,
    };
    let settings = IndexedSpectralAbsorptionImageSettings {
        spectral: limits(12, [2, 3]),
        max_materials: 2,
    };
    let image = render_density_spectral_indexed_radiance_box(
        &grid,
        settings,
        calibration,
        SpectralIllumination::Columns(&backgrounds),
    )
    .unwrap();
    let xyz = render_density_spectral_indexed_xyz_box(
        &grid,
        settings,
        calibration,
        SpectralIllumination::Columns(&backgrounds),
        &weights,
    )
    .unwrap();
    for py in 0..3 {
        for px in 0..2 {
            let mut expected = [0.0; 2];
            for iy in 0..2 {
                for ix in 0..3 {
                    let wx = (((px + 1) as f64 * 1.5).min((ix + 1) as f64)
                        - (px as f64 * 1.5).max(ix as f64))
                    .max(0.0);
                    let wy = (((py + 1) as f64 * (2.0 / 3.0)).min((iy + 1) as f64)
                        - (py as f64 * (2.0 / 3.0)).max(iy as f64))
                    .max(0.0);
                    let column = ix + 3 * iy;
                    for lane in 0..2 {
                        let tau: f64 = (0..2)
                            .map(|z| {
                                let cell = column + 6 * z;
                                density[cell]
                                    * coefficients[materials[cell] as usize * 2 + lane]
                                    * 0.5
                            })
                            .sum();
                        expected[lane] += wx * wy * backgrounds[column * 2 + lane] * (-tau).exp();
                    }
                }
            }
            let pixel = px + 2 * py;
            for lane in 0..2 {
                close(image[pixel * 2 + lane], expected[lane]);
            }
            for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
                close(
                    xyz[pixel][channel],
                    (0..2)
                        .map(|lane| {
                            expected[lane] * weights[lane] * cmf(wavelengths[lane]) / y_integral()
                        })
                        .sum(),
                );
            }
        }
    }
}

#[test]
fn uniform_illumination_agrees_with_retained_opacity_and_indexed_material_paths() {
    let density = [0.0, 2.0, 1.0, 3.0];
    let grid = VolumeGrid::new([2, 1, 2], &density, [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[450.0, 550.0, 650.0], &[0.1, 0.5, 2.0]);
    let background = [2.0, 5.0, 11.0];
    let weights = [20.0, 100.0, 50.0];
    let settings = limits(4, [3, 2]);
    let opacity = render_density_spectral_box(&grid, settings, calibration).unwrap();
    let radiance = render_density_spectral_radiance_box(
        &grid,
        settings,
        calibration,
        SpectralIllumination::Uniform(&background),
    )
    .unwrap();
    for (index, sample) in opacity.samples().iter().enumerate() {
        close(
            radiance[index],
            sample.attenuate(background[index % 3]).unwrap(),
        );
    }
    let xyz = render_density_spectral_xyz_box(
        &grid,
        settings,
        calibration,
        SpectralIllumination::Uniform(&background),
        &weights,
    )
    .unwrap();
    let reference = opacity
        .approximate_transmitted_xyz(&weights, &background)
        .unwrap();
    for (actual, expected) in xyz.iter().flatten().zip(reference.iter().flatten()) {
        close(*actual, *expected);
    }
    let indexed = IndexedDensityAbsorptionSpectrum {
        wavelengths_nm: calibration.wavelengths_nm,
        mass_absorption_m2_kg: calibration.mass_absorption_m2_kg,
        cell_materials: &[0; 4],
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
    };
    let indexed_image = render_density_spectral_indexed_radiance_box(
        &grid,
        IndexedSpectralAbsorptionImageSettings {
            spectral: settings,
            max_materials: 1,
        },
        indexed,
        SpectralIllumination::Uniform(&background),
    )
    .unwrap();
    assert_eq!(radiance, indexed_image);
}

#[test]
fn projection_applies_all_factors_before_underflow_or_column_overflow() {
    let grid = VolumeGrid::new([2, 1, 1], &[1.0, 0.0], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[550.0], &[1400.0]);
    let backgrounds = [1e200, 0.0];
    let raw = render_density_spectral_radiance_box(
        &grid,
        limits(2, [1, 1]),
        calibration,
        SpectralIllumination::Columns(&backgrounds),
    )
    .unwrap();
    assert_eq!(raw, [0.0]);
    let observed = render_density_spectral_xyz_box(
        &grid,
        limits(2, [1, 1]),
        calibration,
        SpectralIllumination::Columns(&backgrounds),
        &[1e300],
    )
    .unwrap();
    for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
        close(
            observed[0][channel],
            (1e200_f64.ln() + 1e300_f64.ln() - 1400.0 + (cmf(550.0) / y_integral()).ln()
                - 2.0_f64.ln())
            .exp(),
        );
    }
    let clear = VolumeGrid::new([8, 1, 1], &[0.0; 8], [0.0; 3], [1.0; 3]);
    let bright = [1e308, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let pixel = render_density_spectral_xyz_box(
        &clear,
        limits(8, [1, 1]),
        spectrum(&[550.0], &[0.0]),
        SpectralIllumination::Columns(&bright),
        &[1000.0],
    )
    .unwrap();
    assert!((1e308 * 1000.0 * cie_y(550.0) / y_integral()).is_infinite());
    for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
        close(
            pixel[0][channel],
            (1e308 / 8.0) * (1000.0 * cmf(550.0) / y_integral()),
        );
    }
}

#[test]
fn bright_backgrounds_survive_zero_unit_transmission_and_finite_means_do_not_overflow() {
    let grid = VolumeGrid::new([2, 1, 1], &[1.0, 1.0], [0.0; 3], [1.0; 3]);
    let radiance = render_density_spectral_radiance_box(
        &grid,
        limits(2, [1, 1]),
        spectrum(&[1000.0, 2000.0], &[800.0, 0.0]),
        SpectralIllumination::Uniform(&[1e300, f64::MAX]),
    )
    .unwrap();
    close(radiance[0], (1e300_f64.ln() - 800.0).exp());
    assert_eq!(radiance[1], f64::MAX);
    let invalid = render_density_spectral_xyz_box(
        &grid,
        limits(2, [1, 1]),
        spectrum(&[1000.0], &[1.0]),
        SpectralIllumination::Uniform(&[1.0]),
        &[1.0],
    );
    assert_eq!(
        invalid,
        Err(DvrError::InvalidObserverWavelength { lane: 0 })
    );
}

#[test]
fn spatial_refinement_preserves_the_same_light_and_density_fields() {
    let coarse = VolumeGrid::new([2, 1, 1], &[2.0, 0.0], [0.0; 3], [1.0; 3]);
    let fine = VolumeGrid::new(
        [4, 2, 1],
        &[2.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0],
        [1e12, -1e12, 0.0],
        [0.5, 0.5, 1.0],
    );
    let calibration = spectrum(&[550.0], &[0.5]);
    let first = render_density_spectral_radiance_box(
        &coarse,
        limits(2, [3, 2]),
        calibration,
        SpectralIllumination::Columns(&[7.0, 2.0]),
    )
    .unwrap();
    let second = render_density_spectral_radiance_box(
        &fine,
        limits(8, [3, 2]),
        calibration,
        SpectralIllumination::Columns(&[7.0, 7.0, 2.0, 2.0, 7.0, 7.0, 2.0, 2.0]),
    )
    .unwrap();
    for (actual, expected) in first.into_iter().zip(second) {
        close(actual, expected);
    }
}

#[test]
fn admission_rejects_bad_lighting_even_when_opaque_or_zero_weighted() {
    let grid = VolumeGrid::new([2, 1, 1], &[1000.0; 2], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[550.0], &[1.0]);
    for illumination in [
        SpectralIllumination::Uniform(&[]),
        SpectralIllumination::Columns(&[1.0]),
    ] {
        assert_eq!(
            render_density_spectral_radiance_box(
                &grid,
                limits(2, [1, 1]),
                calibration,
                illumination
            ),
            Err(DvrError::InvalidBackgroundSize)
        );
    }
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            render_density_spectral_xyz_box(
                &grid,
                limits(2, [1, 1]),
                calibration,
                SpectralIllumination::Columns(&[0.0, invalid]),
                &[0.0]
            ),
            Err(DvrError::InvalidBackgroundRadiance)
        );
        assert_eq!(
            render_density_spectral_xyz_box(
                &grid,
                limits(2, [1, 1]),
                calibration,
                SpectralIllumination::Uniform(&[1.0]),
                &[invalid]
            ),
            Err(DvrError::InvalidSpectralWeights)
        );
    }
    let mut refused = limits(2, [1, 1]);
    refused.max_spectral_evaluations = 1;
    assert_eq!(
        render_density_spectral_radiance_box(
            &grid,
            refused,
            calibration,
            SpectralIllumination::Columns(&[f64::NAN; 2])
        ),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 2,
            limit: 1
        })
    );
    let invalid_map = IndexedDensityAbsorptionSpectrum {
        wavelengths_nm: &[550.0],
        mass_absorption_m2_kg: &[1.0],
        cell_materials: &[0, 1],
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
    };
    assert_eq!(
        render_density_spectral_indexed_radiance_box(
            &grid,
            IndexedSpectralAbsorptionImageSettings {
                spectral: limits(2, [1, 1]),
                max_materials: 1
            },
            invalid_map,
            SpectralIllumination::Uniform(&[0.0])
        ),
        Err(DvrError::AbsorptionMaterialIndexOutOfBounds {
            cell: 1,
            material: 1,
            materials: 1
        })
    );
}
