use super::*;
use crate::volumes::{
    AbsorptionImageSettings, render_density_spectral_box, render_density_spectral_transmittance_box,
};

fn limits(cells: usize, resolution: [usize; 2]) -> IndexedSpectralAbsorptionImageSettings {
    IndexedSpectralAbsorptionImageSettings {
        spectral: SpectralAbsorptionImageSettings {
            image: AbsorptionImageSettings {
                resolution,
                max_pixels: 128,
                max_grid_cells: cells,
                max_cell_visits: 10_000,
            },
            max_spectral_evaluations: 10_000,
        },
        max_materials: 16,
    }
}

fn spectrum<'a>(
    wavelengths: &'a [f64],
    coefficients: &'a [f64],
    materials: &'a [u32],
) -> IndexedDensityAbsorptionSpectrum<'a> {
    IndexedDensityAbsorptionSpectrum {
        wavelengths_nm: wavelengths,
        mass_absorption_m2_kg: coefficients,
        cell_materials: materials,
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
fn g1_indexed_layers_integrate_each_material_with_its_original_density() {
    let grid = VolumeGrid::new([1, 1, 2], &[2.0, 3.0], [0.0; 3], [1.0, 1.0, 0.4]);
    let optics = IndexedDensityAbsorptionSpectrum {
        density_unit_kg_m3: 2.0,
        length_unit_m: 0.5,
        ..spectrum(&[450.0, 650.0], &[0.2, 1.0, 1.7, 0.1], &[0, 1])
    };
    let image = render_density_spectral_indexed_box(&grid, limits(2, [1, 1]), optics).unwrap();
    assert_eq!(image.wavelengths_nm(), &[450.0, 650.0]);
    // Independent two-layer Beer--Lambert law: (2*kappa_0 + 3*kappa_1)*0.4.
    for (sample, tau) in image.samples().iter().zip([2.2_f64, 0.92]) {
        close(sample.optical_depth(), tau);
        close(sample.transmittance(), (-tau).exp());
        close(sample.opacity(), -(-tau).exp_m1());
    }
    // Separately averaging density and material coefficients gives tau=1.9.
    assert!((image.samples()[0].transmittance() - (-1.9_f64).exp()).abs() > 0.03);

    let reversed = VolumeGrid::new([1, 1, 2], &[3.0, 2.0], [0.0; 3], grid.cell);
    let reverse_optics = IndexedDensityAbsorptionSpectrum {
        cell_materials: &[1, 0],
        ..optics
    };
    let same =
        render_density_spectral_indexed_box(&reversed, limits(2, [1, 1]), reverse_optics).unwrap();
    for (sample, reference) in same.samples().iter().zip(image.samples()) {
        close(sample.optical_depth(), reference.optical_depth());
    }
    let changed =
        render_density_spectral_indexed_box(&grid, limits(2, [1, 1]), reverse_optics).unwrap();
    close(changed.samples()[0].optical_depth(), 1.6);
    close(changed.samples()[1].optical_depth(), 1.28);
}

#[test]
fn g1_indexed_material_coverage_is_filtered_after_spectral_attenuation() {
    let density = [1.0; 8];
    let materials = [0, 1, 1, 1, 0, 1, 1, 1];
    let grid = VolumeGrid::new([4, 2, 1], &density, [0.0; 3], [1.0; 3]);
    let optics = spectrum(&[450.0, 650.0], &[1000.0, 0.0, 0.0, 1000.0], &materials);
    let image = render_density_spectral_indexed_box(&grid, limits(8, [1, 1]), optics).unwrap();
    // The two materials absorb opposite wavelengths. Arithmetic averaging of
    // their coefficients would incorrectly make both wavelengths opaque.
    for (sample, expected) in image.samples().iter().zip([0.75_f64, 0.25]) {
        close(sample.transmittance(), expected);
        close(sample.opacity(), 1.0 - expected);
        close(sample.optical_depth(), -expected.ln());
    }
    let fine_materials: Vec<_> = (0..4)
        .flat_map(|y| (0..8).map(move |x| materials[x / 2 + 4 * (y / 2)]))
        .collect();
    let fine = VolumeGrid::new([8, 4, 1], &[1.0; 32], [0.0; 3], [0.5, 0.5, 1.0]);
    let refined = render_density_spectral_indexed_box(
        &fine,
        limits(32, [1, 1]),
        IndexedDensityAbsorptionSpectrum {
            cell_materials: &fine_materials,
            ..optics
        },
    )
    .unwrap();
    for (sample, reference) in refined.samples().iter().zip(image.samples()) {
        close(sample.transmittance(), reference.transmittance());
        close(sample.opacity(), reference.opacity());
    }
}

#[test]
fn g3_uniform_palette_matches_uniform_transport_and_preserves_si_rescaling() {
    let density = [0.5, 1.5, 1.0, 2.0];
    let grid = VolumeGrid::new([2, 1, 2], &density, [0.0; 3], [1.0, 1.0, 0.25]);
    let indexed = spectrum(&[450.0, 650.0], &[0.1, 0.7], &[0; 4]);
    let uniform = DensityAbsorptionSpectrum {
        wavelengths_nm: indexed.wavelengths_nm,
        mass_absorption_m2_kg: indexed.mass_absorption_m2_kg,
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
    };
    let request = limits(4, [3, 2]);
    let image = render_density_spectral_indexed_box(&grid, request, indexed).unwrap();
    let reference = render_density_spectral_box(&grid, request.spectral, uniform).unwrap();
    for (sample, reference) in image.samples().iter().zip(reference.samples()) {
        close(sample.transmittance(), reference.transmittance());
        close(sample.opacity(), reference.opacity());
        close(sample.optical_depth(), reference.optical_depth());
    }
    let direct =
        render_density_spectral_indexed_transmittance_box(&grid, request, indexed).unwrap();
    let uniform_direct =
        render_density_spectral_transmittance_box(&grid, request.spectral, uniform).unwrap();
    assert_eq!(direct.len(), 12);
    for ((value, sample), uniform) in direct.iter().zip(image.samples()).zip(uniform_direct) {
        assert_eq!(*value, sample.transmittance());
        close(*value, uniform);
    }

    let rescaled_density = density.map(|value| 2.0 * value);
    let rescaled_grid = VolumeGrid::new([2, 1, 2], &rescaled_density, [0.0; 3], [0.5, 0.5, 0.125]);
    let rescaled = render_density_spectral_indexed_box(
        &rescaled_grid,
        request,
        IndexedDensityAbsorptionSpectrum {
            density_unit_kg_m3: 0.5,
            length_unit_m: 2.0,
            ..indexed
        },
    )
    .unwrap();
    assert_eq!(rescaled, image);
}

#[test]
fn g1_indexed_range_stable_products_preserve_thin_and_opaque_lanes() {
    let grid = VolumeGrid::new([1, 1, 2], &[f64::MAX; 2], [0.0; 3], [1.0; 3]);
    let finite = render_density_spectral_indexed_box(
        &grid,
        limits(2, [1, 1]),
        spectrum(&[500.0], &[1e-308, 2e-308], &[0, 1]),
    )
    .unwrap();
    let tau = f64::MAX * 1e-308 + f64::MAX * 2e-308;
    close(finite.samples()[0].optical_depth(), tau);
    close(finite.samples()[0].transmittance(), (-tau).exp());

    let density = 1e-320;
    let thin_density = [density, 1e-297];
    let thin_grid = VolumeGrid::new([1, 1, 2], &thin_density, [0.0; 3], [1.0, 1.0, 1e300]);
    let image = render_density_spectral_indexed_box(
        &thin_grid,
        limits(2, [1, 1]),
        spectrum(&[450.0, 650.0], &[1.0, 0.0, 0.0, 0.8], &[0, 1]),
    )
    .unwrap();
    let thin_depth = density * 1e300;
    assert_eq!(image.samples()[0].transmittance(), 1.0);
    close(image.samples()[0].opacity(), thin_depth);
    close(image.samples()[0].optical_depth(), thin_depth);
    let opaque_depth = (thin_density[1] * 0.8) * 1e300;
    assert_eq!(image.samples()[1].transmittance(), 0.0);
    close(image.samples()[1].optical_depth(), opaque_depth);
    close(
        image.samples()[1].attenuate(1e300).unwrap(),
        (1e300_f64.ln() - opaque_depth).exp(),
    );
}

#[test]
fn g0_indexed_admission_bounds_every_material_and_all_original_cell_indices() {
    let grid = VolumeGrid::new([1, 1, 2], &[0.0, 1.0], [0.0; 3], [1.0; 3]);
    let optics = spectrum(&[450.0, 650.0], &[0.4, 0.1, 0.2, 0.7], &[0, 1]);
    let request = limits(2, [1, 1]);
    assert_eq!(
        render_density_spectral_indexed_box(
            &grid,
            IndexedSpectralAbsorptionImageSettings {
                max_materials: 1,
                ..request
            },
            optics,
        ),
        Err(DvrError::AbsorptionMaterialBudgetExceeded {
            requested: 2,
            limit: 1
        }),
    );
    assert_eq!(
        render_density_spectral_indexed_box(
            &grid,
            request,
            IndexedDensityAbsorptionSpectrum {
                cell_materials: &[0],
                ..optics
            },
        ),
        Err(DvrError::InvalidAbsorptionMaterialMap),
    );
    // An invalid material in a zero-density cell cannot silently disappear.
    assert_eq!(
        render_density_spectral_indexed_box(
            &grid,
            request,
            IndexedDensityAbsorptionSpectrum {
                cell_materials: &[u32::MAX, 0],
                ..optics
            },
        ),
        Err(DvrError::AbsorptionMaterialIndexOutOfBounds {
            cell: 0,
            material: u32::MAX,
            materials: 2,
        }),
    );
    for coefficients in [&[0.4, 0.1, 0.2][..], &[0.4, 0.1, f64::NAN, 0.7][..]] {
        assert_eq!(
            render_density_spectral_indexed_transmittance_box(
                &grid,
                request,
                IndexedDensityAbsorptionSpectrum {
                    mass_absorption_m2_kg: coefficients,
                    cell_materials: &[0, 0],
                    ..optics
                },
            ),
            Err(DvrError::InvalidAbsorptionSpectrum),
        );
    }
    assert_eq!(
        render_density_spectral_indexed_box(&grid, request, spectrum(&[], &[], &[0, 0]),),
        Err(DvrError::InvalidAbsorptionSpectrum),
    );
}

#[test]
fn g0_indexed_spectral_work_counts_layers_before_scanning_cells_or_allocating() {
    let grid = VolumeGrid::new([3, 2, 2], &[1.0; 12], [0.0; 3], [1.0; 3]);
    let optics = spectrum(&[450.0, 650.0], &[0.1, 0.7], &[0; 12]);
    // 4 horizontal overlaps * 4 vertical overlaps * 2 layers * 2 wavelengths.
    let exact = IndexedSpectralAbsorptionImageSettings {
        spectral: SpectralAbsorptionImageSettings {
            image: AbsorptionImageSettings {
                max_cell_visits: 32,
                ..limits(12, [2, 3]).spectral.image
            },
            max_spectral_evaluations: 64,
        },
        max_materials: 1,
    };
    assert_eq!(
        render_density_spectral_indexed_box(&grid, exact, optics)
            .unwrap()
            .samples()
            .len(),
        12,
    );
    let invalid = VolumeGrid::new([3, 2, 2], &[f64::NAN; 12], [0.0; 3], [1.0; 3]);
    let invalid_optics = IndexedDensityAbsorptionSpectrum {
        cell_materials: &[u32::MAX; 12],
        ..optics
    };
    let insufficient = IndexedSpectralAbsorptionImageSettings {
        spectral: SpectralAbsorptionImageSettings {
            max_spectral_evaluations: 63,
            ..exact.spectral
        },
        ..exact
    };
    let expected = Err(DvrError::SpectralEvaluationBudgetExceeded {
        requested: 64,
        limit: 63,
    });
    assert_eq!(
        render_density_spectral_indexed_box(&invalid, insufficient, invalid_optics),
        expected,
    );
    assert_eq!(
        render_density_spectral_indexed_transmittance_box(&invalid, insufficient, invalid_optics,),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 64,
            limit: 63
        }),
    );
}
