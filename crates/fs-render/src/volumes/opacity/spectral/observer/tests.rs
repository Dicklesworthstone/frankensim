use super::*;
use crate::{
    spectral::{cie_x, cie_y, cie_z},
    volumes::{
        AbsorptionImageSettings, DensityAbsorptionSpectrum, SpectralAbsorptionImageSettings,
        VolumeGrid, render_density_spectral_box,
    },
};
use fs_math::det;

fn image(
    data: &[f64],
    width: usize,
    wavelengths: &[f64],
    absorption: &[f64],
) -> SpectralOpacityImage {
    render_density_spectral_box(
        &VolumeGrid::new([data.len(), 1, 1], data, [0.0; 3], [1.0; 3]),
        SpectralAbsorptionImageSettings {
            image: AbsorptionImageSettings {
                resolution: [width, 1],
                max_pixels: width,
                max_grid_cells: data.len(),
                max_cell_visits: 100,
            },
            max_spectral_evaluations: 100,
        },
        DensityAbsorptionSpectrum {
            wavelengths_nm: wavelengths,
            mass_absorption_m2_kg: absorption,
            density_unit_kg_m3: 1.0,
            length_unit_m: 1.0,
        },
    )
    .unwrap()
}

fn close(actual: [f64; 3], expected: [f64; 3]) {
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= expected.abs() * 3e-12,
            "{actual:e} != {expected:e}"
        );
    }
}

#[test]
fn cold_colored_absorber_projects_patterned_illumination_with_explicit_nm_weights() {
    let wavelengths = [450.0, 550.0, 650.0];
    let absorption = [0.25, 1.0, 4.0];
    let density = [0.5, 2.0];
    let weights = [25.0, 100.0, 50.0];
    let backgrounds = [[2.0, 5.0, 11.0], [13.0, 7.0, 3.0]];
    let opacity = image(&density, 2, &wavelengths, &absorption);
    let observed = opacity
        .approximate_transmitted_xyz_over_image(&weights, backgrounds.as_flattened())
        .unwrap();
    for pixel in 0..2 {
        let mut expected = [0.0; 3];
        for (channel, cmf) in [cie_x, cie_y, cie_z].into_iter().enumerate() {
            for lane in 0..3 {
                expected[channel] += backgrounds[pixel][lane]
                    * det::exp(-density[pixel] * absorption[lane])
                    * weights[lane]
                    * cmf(wavelengths[lane])
                    / y_integral();
            }
        }
        close(observed[pixel], expected);
    }
    let flat = opacity
        .approximate_transmitted_xyz(&weights, &backgrounds[0])
        .unwrap();
    let repeated = [backgrounds[0], backgrounds[0]];
    assert_eq!(
        flat,
        opacity
            .approximate_transmitted_xyz_over_image(&weights, repeated.as_flattened())
            .unwrap()
    );
    let twice = opacity
        .approximate_transmitted_xyz(&weights.map(|weight| weight * 2.0), &backgrounds[0])
        .unwrap();
    for pixel in 0..2 {
        close(twice[pixel], flat[pixel].map(|channel| channel * 2.0));
    }
}

#[test]
fn visible_projection_preserves_footprint_coverage_and_opaque_bright_backgrounds() {
    let opacity = image(&[1000.0, 0.0, 0.0, 0.0], 1, &[550.0], &[1.0]);
    let observed = opacity
        .approximate_transmitted_xyz(&[50.0], &[2.0])
        .unwrap();
    close(
        observed[0],
        [cie_x, cie_y, cie_z].map(|cmf| 0.75 * 2.0 * 50.0 * cmf(550.0) / y_integral()),
    );
    let opaque = image(&[1000.0], 1, &[550.0], &[1.0]);
    assert_eq!(opaque.samples()[0].transmittance(), 0.0);
    let observed = opaque
        .approximate_transmitted_xyz(&[50.0], &[1e300])
        .unwrap();
    let survived = det::exp(det::ln(1e300) - 1000.0);
    close(
        observed[0],
        [cie_x, cie_y, cie_z].map(|cmf| survived * 50.0 * cmf(550.0) / y_integral()),
    );
    assert!(observed[0].iter().all(|value| *value > 0.0));
}

#[test]
fn projection_rejects_ambiguous_layout_invalid_spectra_and_nonvisible_nodes() {
    let opacity = image(&[1.0, 2.0], 2, &[450.0, 650.0], &[1.0, 2.0]);
    assert_eq!(
        opacity.approximate_transmitted_xyz(&[1.0, 1.0], &[1.0]),
        Err(DvrError::InvalidBackgroundSize)
    );
    assert_eq!(
        opacity.approximate_transmitted_xyz_over_image(&[1.0, 1.0], &[1.0; 3]),
        Err(DvrError::InvalidBackgroundSize)
    );
    assert_eq!(
        opacity.approximate_transmitted_xyz(&[1.0], &[1.0; 2]),
        Err(DvrError::InvalidSpectralWeights)
    );
    for invalid in [-1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            opacity.approximate_transmitted_xyz(&[invalid, 1.0], &[1.0; 2]),
            Err(DvrError::InvalidSpectralWeights)
        );
        assert_eq!(
            opacity.approximate_transmitted_xyz(&[0.0, 0.0], &[1.0, invalid]),
            Err(DvrError::InvalidBackgroundRadiance)
        );
    }
    let infrared = image(&[1.0], 1, &[900.0], &[1.0]);
    assert!(infrared.samples()[0].transmittance() > 0.0);
    assert_eq!(
        infrared.approximate_transmitted_xyz(&[0.0], &[1.0]),
        Err(DvrError::InvalidObserverWavelength { lane: 0 })
    );
}
