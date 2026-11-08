use super::*;
use crate::volumes::TransferPoint;

fn settings(resolution: [usize; 2], cells: usize) -> AbsorptionImageSettings {
    AbsorptionImageSettings {
        resolution,
        max_pixels: resolution[0] * resolution[1],
        max_grid_cells: cells,
        max_cell_visits: 100_000,
    }
}

fn calibration() -> DensityAbsorption {
    DensityAbsorption {
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
        mass_absorption_m2_kg: 1.0,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs() + f64::MIN_POSITIVE,
        "{actual:e} != {expected:e}"
    );
}

/// G1: unresolved coverage averages surviving light, not optical depth/density.
#[test]
fn clear_and_opaque_columns_retain_their_exact_coverage() {
    let values = [0.0, 1000.0];
    let grid = VolumeGrid::new([2, 1, 1], &values, [0.0; 3], [1.0; 3]);
    let pixel = render_density_transport_box(&grid, settings([1, 1], 2), calibration()).unwrap()[0];
    assert_eq!(pixel.transmittance(), 0.5);
    assert_eq!(pixel.opacity(), 0.5);
    close(pixel.optical_depth(), 2.0_f64.ln());
    assert_eq!(pixel.composite([10.0, 4.0, 2.0]).unwrap(), [5.0, 2.0, 1.0]);
    assert_eq!(
        render_density_transmittance_box(&grid, settings([1, 1], 2), calibration()).unwrap(),
        vec![0.5]
    );
}

/// G1: a three-cell field and two pixels share the middle cell equally.
#[test]
fn nondivisible_pixels_resolve_exact_two_dimensional_areas() {
    let values = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let grid = VolumeGrid::new([3, 2, 1], &values, [-8.0, 14.0, 5.0], [0.7, 8.0, 0.5]);
    let image = render_density_transport_box(&grid, settings([2, 3], 6), calibration()).unwrap();
    let transmission = values.map(|rho| (-0.5 * rho).exp());
    let rows = [
        [
            (2.0 * transmission[0] + transmission[1]) / 3.0,
            (transmission[1] + 2.0 * transmission[2]) / 3.0,
        ],
        [
            (2.0 * transmission[3] + transmission[4]) / 3.0,
            (transmission[4] + 2.0 * transmission[5]) / 3.0,
        ],
    ];
    for x in 0..2 {
        close(image[x].transmittance(), rows[0][x]);
        close(image[x + 2].transmittance(), 0.5 * (rows[0][x] + rows[1][x]));
        close(image[x + 4].transmittance(), rows[1][x]);
    }
}

/// G1: integrate each original layered column before filtering its light.
#[test]
fn emission_and_attenuation_keep_the_original_layer_correlations() {
    let points = [
        TransferPoint { scalar: 0.0, extinction: 0.3, source_rgb: [10.0, 0.0, 4.0] },
        TransferPoint { scalar: 1.0, extinction: 1.2, source_rgb: [0.0, 7.0, 2.0] },
    ];
    let transfer = TransferFunction::new(points.to_vec()).unwrap();
    let values = [0.0, 1.0, 1.0, 0.0];
    let grid = VolumeGrid::new([2, 1, 2], &values, [0.0; 3], [1.0; 3]);
    let pixel = render_transfer_transport_box(&grid, &transfer, settings([1, 1], 4)).unwrap()[0];
    let (a, b) = ((-0.3_f64).exp(), (-1.2_f64).exp());
    close(pixel.transmittance(), a * b);
    for channel in 0..3 {
        let red_front = points[0].source_rgb[channel] * (1.0 - a)
            + a * points[1].source_rgb[channel] * (1.0 - b);
        let green_front = points[1].source_rgb[channel] * (1.0 - b)
            + b * points[0].source_rgb[channel] * (1.0 - a);
        close(pixel.emission_rgb()[channel], 0.5 * (red_front + green_front));
    }
}

/// G3: neither lateral cell refinement nor SI unit choice changes the image.
#[test]
fn refinement_translation_and_units_preserve_physical_pixels() {
    let values = [0.0, 0.4, 3.0];
    let refined = [0.0, 0.0, 400.0, 400.0, 3000.0, 3000.0];
    let grid = VolumeGrid::new([3, 1, 1], &values, [0.0; 3], [2.0, 1.0, 0.8]);
    let fine = VolumeGrid::new([6, 1, 1], &refined, [100.0, -50.0, 4.0], [1.0, 1.0, 0.8]);
    let a = render_density_transport_box(&grid, settings([2, 1], 3), calibration()).unwrap();
    let b = render_density_transport_box(
        &fine, settings([2, 1], 6),
        DensityAbsorption { density_unit_kg_m3: 0.001, ..calibration() },
    ).unwrap();
    for (a, b) in a.into_iter().zip(b) {
        close(a.transmittance(), b.transmittance());
        close(a.opacity(), b.opacity());
        close(a.optical_depth(), b.optical_depth());
    }
    let moved = VolumeGrid::new([3, 1, 1], &values, [1e15; 3], [2.0, 1.0, 0.8]);
    assert_eq!(
        render_density_transport_box(&grid, settings([7, 1], 3), calibration()).unwrap(),
        render_density_transport_box(&moved, settings([7, 1], 3), calibration()).unwrap()
    );
}

#[test]
fn thin_opacity_and_bright_backgrounds_survive_footprint_averaging() {
    let thin_values = [0.0, 1e-20];
    let thin = VolumeGrid::new([2, 1, 1], &thin_values, [0.0; 3], [1.0; 3]);
    let pixel = render_density_transport_box(&thin, settings([1, 1], 2), calibration()).unwrap()[0];
    assert_eq!(pixel.transmittance(), 1.0);
    close(pixel.opacity(), 5e-21);
    close(pixel.optical_depth(), 5e-21);
    let opaque_values = [800.0, 801.0];
    let opaque = VolumeGrid::new([2, 1, 1], &opaque_values, [0.0; 3], [1.0; 3]);
    let pixel = render_density_transport_box(&opaque, settings([1, 1], 2), calibration()).unwrap()[0];
    assert_eq!(pixel.transmittance(), 0.0);
    assert_eq!(pixel.opacity(), 1.0);
    close(pixel.optical_depth(), 800.0 - ((1.0 + (-1.0_f64).exp()) / 2.0).ln());
    let expected = 0.5 * ((1e300_f64.ln() - 800.0).exp() + (1e300_f64.ln() - 801.0).exp());
    for actual in pixel.composite([1e300; 3]).unwrap() {
        close(actual, expected);
    }
}

#[test]
fn homogeneous_furnace_and_subnormal_sources_are_preserved() {
    let values = [0.0; 15];
    let grid = VolumeGrid::new([5, 3, 1], &values, [0.0; 3], [1.0; 3]);
    for source in [[2.0, 5.0, 11.0], [f64::from_bits(1); 3], [f64::MAX; 3]] {
        let transfer = TransferFunction::new(vec![TransferPoint {
            scalar: 0.0, extinction: 1000.0, source_rgb: source,
        }]).unwrap();
        let pixels = render_transfer_transport_box(&grid, &transfer, settings([7, 2], 15)).unwrap();
        for pixel in pixels {
            assert_eq!(pixel.emission_rgb(), source);
            assert_eq!(pixel.composite(source).unwrap(), source);
        }
    }
}

#[test]
fn exact_footprint_work_is_admitted_before_any_output() {
    let values = [0.0; 60];
    let grid = VolumeGrid::new([5, 3, 4], &values, [0.0; 3], [1.0; 3]);
    // (5+7-1)*(3+2-1)*4 = 176 column-cell visits.
    let limits = AbsorptionImageSettings {
        max_cell_visits: 176, ..settings([7, 2], 60)
    };
    assert_eq!(
        render_density_transmittance_box(&grid, limits, calibration()).unwrap(),
        vec![1.0; 14]
    );
    assert_eq!(
        render_density_transmittance_box(
            &grid, AbsorptionImageSettings { max_cell_visits: 175, ..limits }, calibration(),
        ),
        Err(DvrError::CellVisitBudgetExceeded { requested: 176, limit: 175 })
    );
    let invalid = VolumeGrid { data: &[f64::NAN; 60], ..grid };
    assert_eq!(
        render_density_transmittance_box(
            &invalid, AbsorptionImageSettings { max_cell_visits: 0, ..limits }, calibration(),
        ),
        Err(DvrError::CellVisitBudgetExceeded { requested: 176, limit: 0 })
    );
}

#[test]
fn integer_footprints_cover_every_pixel_without_gaps_or_double_counting() {
    for cells in 1..32 {
        for pixels in 1..32 {
            let mut visits = 0;
            for pixel in 0..pixels {
                let segments = overlaps(pixel, pixels, cells).collect::<Vec<_>>();
                assert!(segments.iter().all(|&(cell, weight)| cell < cells && weight > 0));
                assert_eq!(segments.iter().map(|&(_, weight)| weight).sum::<usize>(), cells);
                visits += segments.len() as u128;
            }
            assert_eq!(visits, overlap_count(cells, pixels));
        }
    }
    assert_eq!(
        overlaps(usize::MAX - 1, usize::MAX, usize::MAX).collect::<Vec<_>>(),
        vec![(usize::MAX - 1, usize::MAX)]
    );
}
