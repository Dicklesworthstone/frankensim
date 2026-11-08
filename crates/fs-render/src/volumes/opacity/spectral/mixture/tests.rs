use super::*;
use crate::volumes::{
    AbsorptionImageSettings, IndexedDensityAbsorptionSpectrum,
    IndexedSpectralAbsorptionImageSettings, VolumeGrid, render_density_spectral_indexed_box,
};

fn settings() -> AbsorptionMixtureSettings {
    AbsorptionMixtureSettings {
        max_components: 4,
        max_mixtures: 4,
        max_spectral_evaluations: 64,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn preparation_keeps_mass_basis_and_original_point_spectra() {
    let wavelengths = [450.0, 650.0];
    let components = [2.0, 4.0, 10.0, 20.0];
    // A pure component, total-density mass fractions, carrier-relative
    // loadings (deliberately not normalized), and explicitly no absorber.
    let recipes = [1.0, 0.0, 0.25, 0.75, 2.0, 0.5, 0.0, 0.0];
    let palette =
        prepare_density_absorption_mixtures(&wavelengths, &components, &recipes, settings())
            .unwrap();
    assert_eq!(palette, [2.0, 4.0, 8.0, 16.0, 9.0, 18.0, 0.0, 0.0]);
}

#[test]
fn homogeneous_mass_mixing_and_unresolved_material_coverage_remain_distinct() {
    let wavelengths = [550.0];
    let components = [0.0, 2.0 * core::f64::consts::LN_2];
    let homogeneous =
        prepare_density_absorption_mixtures(&wavelengths, &components, &[0.5, 0.5], settings())
            .unwrap();
    let render = |density: &[f64],
                  material: &[u32],
                  coefficients: &[f64],
                  density_unit_kg_m3: f64,
                  depth_m: f64| {
        render_density_spectral_indexed_box(
            &VolumeGrid::new(
                [density.len(), 1, 1],
                density,
                [0.0; 3],
                [1.0, 1.0, depth_m],
            ),
            IndexedSpectralAbsorptionImageSettings {
                spectral: super::super::SpectralAbsorptionImageSettings {
                    image: AbsorptionImageSettings {
                        resolution: [1, 1],
                        max_pixels: 1,
                        max_grid_cells: density.len(),
                        max_cell_visits: density.len() as u64,
                    },
                    max_spectral_evaluations: density.len() as u64,
                },
                max_materials: 2,
            },
            IndexedDensityAbsorptionSpectrum {
                wavelengths_nm: &wavelengths,
                mass_absorption_m2_kg: coefficients,
                cell_materials: material,
                density_unit_kg_m3,
                length_unit_m: 1.0,
            },
        )
        .unwrap()
        .samples()[0]
    };
    let mixed = render(&[1.0], &[0], &homogeneous, 1.0, 1.0);
    let covered = render(&[1.0, 1.0], &[0, 1], &components, 1.0, 1.0);
    close(mixed.transmittance(), 0.5);
    close(covered.transmittance(), 0.625);
    close(mixed.opacity(), 0.5);
    // Physical density and path units retain the same areal reference mass.
    let rescaled = render(&[0.5], &[0], &homogeneous, 1000.0, 0.002);
    close(rescaled.transmittance(), mixed.transmittance());
}

#[test]
fn preparation_admits_full_work_and_validates_zero_weight_entries() {
    let wavelengths = [450.0, 650.0];
    let mut limits = settings();
    limits.max_spectral_evaluations = 7;
    // Two mixtures × two components × two wavelengths: refuse before scanning
    // the later invalid values, including recipes that would be transparent.
    assert_eq!(
        prepare_density_absorption_mixtures(
            &wavelengths,
            &[1.0, 2.0, f64::NAN, 4.0],
            &[1.0, 0.0, 0.0, 0.0],
            limits,
        ),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 8,
            limit: 7
        }),
    );
    limits.max_spectral_evaluations = 8;
    assert_eq!(
        prepare_density_absorption_mixtures(
            &wavelengths,
            &[1.0, 2.0, f64::NAN, 4.0],
            &[1.0, 0.0, 0.0, 0.0],
            limits,
        ),
        Err(DvrError::InvalidAbsorptionSpectrum),
    );
    for weight in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            prepare_density_absorption_mixtures(
                &wavelengths,
                &[1.0, 2.0, 3.0, 4.0],
                &[0.0, weight],
                settings(),
            ),
            Err(DvrError::InvalidAbsorptionMixtureWeights),
        );
    }
    assert_eq!(
        prepare_density_absorption_mixtures(
            &wavelengths,
            &[1.0, 2.0, 3.0, 4.0],
            &[1.0],
            settings()
        ),
        Err(DvrError::InvalidAbsorptionMixtureWeights),
    );
    assert_eq!(
        prepare_density_absorption_mixtures(
            &wavelengths,
            &[1.0, 2.0, 3.0, 4.0],
            &[1.0, 0.0],
            AbsorptionMixtureSettings {
                max_components: 1,
                ..settings()
            }
        ),
        Err(DvrError::AbsorptionComponentBudgetExceeded {
            requested: 2,
            limit: 1
        }),
    );
    assert_eq!(
        prepare_density_absorption_mixtures(
            &wavelengths,
            &[1.0, 2.0],
            &[1.0, 0.5],
            AbsorptionMixtureSettings {
                max_mixtures: 1,
                ..settings()
            }
        ),
        Err(DvrError::AbsorptionMaterialBudgetExceeded {
            requested: 2,
            limit: 1
        }),
    );
}

#[test]
fn preparation_retains_representable_tiny_absorption_and_refuses_lost_positive_terms() {
    let wavelengths = [550.0];
    let tiny =
        prepare_density_absorption_mixtures(&wavelengths, &[f64::MIN_POSITIVE], &[0.5], settings())
            .unwrap();
    assert_eq!(tiny, [f64::MIN_POSITIVE * 0.5]);
    let image = render_density_spectral_indexed_box(
        &VolumeGrid::new([1, 1, 1], &[1e300], [0.0; 3], [1.0; 3]),
        IndexedSpectralAbsorptionImageSettings {
            spectral: super::super::SpectralAbsorptionImageSettings {
                image: AbsorptionImageSettings {
                    resolution: [1, 1],
                    max_pixels: 1,
                    max_grid_cells: 1,
                    max_cell_visits: 1,
                },
                max_spectral_evaluations: 1,
            },
            max_materials: 1,
        },
        IndexedDensityAbsorptionSpectrum {
            wavelengths_nm: &wavelengths,
            mass_absorption_m2_kg: &tiny,
            cell_materials: &[0],
            density_unit_kg_m3: 1.0,
            length_unit_m: 1.0,
        },
    )
    .unwrap();
    let depth = f64::MIN_POSITIVE * 0.5 * 1e300;
    close(image.samples()[0].optical_depth(), depth);
    assert!(image.samples()[0].opacity() > 0.0);
    for (coefficients, weights) in [
        (vec![f64::from_bits(1)], vec![0.5]),
        (vec![f64::MAX], vec![2.0]),
        (vec![f64::MAX, f64::MAX], vec![0.75, 0.75]),
    ] {
        assert_eq!(
            prepare_density_absorption_mixtures(&wavelengths, &coefficients, &weights, settings()),
            Err(DvrError::UnrepresentableAbsorptionMixture),
        );
    }
}
