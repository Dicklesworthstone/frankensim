use super::*;
use crate::volumes::TransferPoint;

fn settings(cells: usize) -> AbsorptionImageSettings {
    AbsorptionImageSettings {
        resolution: [1, 1],
        max_pixels: 1,
        max_grid_cells: cells,
        max_cell_visits: cells as u64,
    }
}

fn close(actual: f64, expected: f64) {
    let tolerance = 2e-12 * expected.abs() + f64::MIN_POSITIVE;
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn constant(extinction: f64, source_rgb: [f64; 3]) -> TransferFunction {
    TransferFunction::new(vec![TransferPoint {
        scalar: 0.0,
        extinction,
        source_rgb,
    }])
    .unwrap()
}

fn calibration(kappa: f64) -> DensityAbsorption {
    DensityAbsorption {
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
        mass_absorption_m2_kg: kappa,
    }
}

#[test]
fn two_layer_sources_follow_camera_order_and_background_attenuation() {
    let transfer = TransferFunction::new(vec![
        TransferPoint { scalar: 0.0, extinction: 0.7, source_rgb: [10.0, 0.0, 3.0] },
        TransferPoint { scalar: 1.0, extinction: 1.1, source_rgb: [0.0, 4.0, 1.0] },
    ])
    .unwrap();
    // The camera faces -z: index 1 is the foreground.
    let values = [1.0, 0.0];
    let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0; 3]);
    let pixel = render_transfer_transport(&grid, &transfer, settings(2)).unwrap()[0];
    let front = -(-0.7_f64).exp_m1();
    let back = (-0.7_f64).exp() * -(-1.1_f64).exp_m1();
    let expected = [10.0 * front, 4.0 * back, 3.0 * front + back];
    for (actual, expected) in pixel.emission_rgb().into_iter().zip(expected) {
        close(actual, expected);
    }
    close(pixel.transmittance(), (-1.8_f64).exp());
    close(pixel.opacity(), -(-1.8_f64).exp_m1());
    let background = [2.0, 3.0, 5.0];
    let composite = pixel.composite(background).unwrap();
    for channel in 0..3 {
        close(composite[channel], expected[channel] + background[channel] * (-1.8_f64).exp());
    }
}

#[test]
fn source_equal_to_background_preserves_radiance() {
    let source = [2.0, 4.0, 8.0];
    let transfer = TransferFunction::new(vec![
        TransferPoint { scalar: 0.0, extinction: 0.01, source_rgb: source },
        TransferPoint { scalar: 1.0, extinction: 4.0, source_rgb: source },
    ])
    .unwrap();
    let values = [0.0, 0.1, 0.5, 0.7, 1.0, 0.3, 0.1];
    let grid = VolumeGrid::new([1, 1, 7], &values, [0.0; 3], [0.2; 3]);
    let pixel = render_transfer_transport(&grid, &transfer, settings(7)).unwrap()[0];
    for (actual, expected) in pixel.composite(source).unwrap().into_iter().zip(source) {
        close(actual, expected);
    }
}

#[test]
fn thin_luminous_cells_preserve_opacity_below_transmittance_precision() {
    let values = [0.0];
    let grid = VolumeGrid::new([1; 3], &values, [0.0; 3], [1.0; 3]);
    let pixel = render_transfer_transport(&grid, &constant(1e-20, [1e12; 3]), settings(1))
        .unwrap()[0];
    assert_eq!(pixel.transmittance(), 1.0);
    close(pixel.opacity(), 1e-20);
    for actual in pixel.emission_rgb() {
        close(actual, 1e-8);
    }
}

#[test]
fn vacuum_and_opaque_limits_remain_finite() {
    let values = [0.0];
    let grid = VolumeGrid::new([1; 3], &values, [0.0; 3], [2.0; 3]);
    let clear = render_transfer_transport(&grid, &constant(0.0, [4.0; 3]), settings(1))
        .unwrap()[0];
    assert_eq!(clear.emission_rgb(), [0.0; 3]);
    assert_eq!(clear.opacity(), 0.0);
    assert_eq!(clear.composite([2.0; 3]).unwrap(), [2.0; 3]);
    let opaque = render_transfer_transport(
        &grid, &constant(f64::MAX, [4.0; 3]), settings(1),
    )
    .unwrap()[0];
    assert_eq!(opaque.emission_rgb(), [4.0; 3]);
    assert_eq!(opaque.transmittance(), 0.0);
    assert_eq!(opaque.opacity(), 1.0);
}

#[test]
fn subdividing_physical_depth_preserves_transfer() {
    let transfer = constant(0.8, [1.0, 3.0, 5.0]);
    let coarse_values = [0.0];
    let fine_values = vec![0.0; 1024];
    let coarse = VolumeGrid::new([1; 3], &coarse_values, [0.0; 3], [1.0, 1.0, 2.0]);
    let fine = VolumeGrid::new([1, 1, 1024], &fine_values, [0.0; 3], [1.0, 1.0, 2.0 / 1024.0]);
    let a = render_transfer_transport(&coarse, &transfer, settings(1)).unwrap()[0];
    let b = render_transfer_transport(&fine, &transfer, settings(1024)).unwrap()[0];
    close(a.transmittance(), b.transmittance());
    close(a.opacity(), b.opacity());
    for (actual, expected) in a.emission_rgb().into_iter().zip(b.emission_rgb()) {
        close(actual, expected);
    }
}

#[test]
fn translation_keeps_integer_pixel_column_ownership() {
    let values = [0.0, 1.0, 0.0];
    let transfer = TransferFunction::new(vec![
        TransferPoint { scalar: 0.0, extinction: 0.1, source_rgb: [1.0, 0.0, 0.0] },
        TransferPoint { scalar: 1.0, extinction: 2.0, source_rgb: [0.0, 1.0, 0.0] },
    ])
    .unwrap();
    let a = VolumeGrid::new([3, 1, 1], &values, [0.0; 3], [1.0; 3]);
    let b = VolumeGrid::new([3, 1, 1], &values, [1e15; 3], [1.0; 3]);
    let settings = AbsorptionImageSettings {
        resolution: [7, 1], max_pixels: 7, max_grid_cells: 3, max_cell_visits: 7,
    };
    assert_eq!(
        render_transfer_transport(&a, &transfer, settings).unwrap(),
        render_transfer_transport(&b, &transfer, settings).unwrap()
    );
}

#[test]
fn mass_density_calibration_obeys_si_optical_depth() {
    let values = [0.25, 1.0];
    let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0, 1.0, 2.0]);
    let units = DensityAbsorption {
        density_unit_kg_m3: 1000.0, length_unit_m: 0.01, mass_absorption_m2_kg: 0.2,
    };
    let image = render_density_transmittance(&grid, settings(2), units).unwrap();
    close(image[0], (-5.0_f64).exp());
    let twice = render_density_transmittance(
        &grid, settings(2), DensityAbsorption { length_unit_m: 0.02, ..units },
    )
    .unwrap();
    close(twice[0], image[0] * image[0]);
    assert_eq!(
        render_density_transmittance(
            &grid, settings(2), DensityAbsorption { mass_absorption_m2_kg: 0.0, ..units },
        )
        .unwrap(),
        vec![1.0]
    );
}

#[test]
fn rescaling_density_units_and_mesh_keeps_physical_opacity() {
    let coarse_values = [2.0];
    let fine_values = [2000.0; 8];
    let coarse = VolumeGrid::new([1; 3], &coarse_values, [0.0; 3], [1.0, 1.0, 4.0]);
    let fine = VolumeGrid::new([1, 1, 8], &fine_values, [0.0; 3], [1.0, 1.0, 0.5]);
    let a = render_density_transmittance(&coarse, settings(1), calibration(0.1)).unwrap()[0];
    let b = render_density_transmittance(
        &fine, settings(8), DensityAbsorption { density_unit_kg_m3: 0.001, ..calibration(0.1) },
    )
    .unwrap()[0];
    close(a, b);
}

#[test]
fn opacity_product_does_not_lose_finite_scale_compensation() {
    let values = [1e200];
    let grid = VolumeGrid::new([1; 3], &values, [0.0; 3], [1.0, 1.0, 1e200]);
    let units = DensityAbsorption {
        density_unit_kg_m3: 1e-200, length_unit_m: 1e-200, mass_absorption_m2_kg: 1.0,
    };
    close(render_density_transmittance(&grid, settings(1), units).unwrap()[0], (-1.0_f64).exp());
}

#[test]
fn bright_sources_and_backgrounds_survive_subnormal_foreground_transmission() {
    for depth in [744.0, 800.0] {
        let values = [1.0, 0.0];
        let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0; 3]);
        let transfer = TransferFunction::new(vec![
            TransferPoint { scalar: 0.0, extinction: depth, source_rgb: [0.0; 3] },
            TransferPoint { scalar: 1.0, extinction: 1.0, source_rgb: [1e300; 3] },
        ])
        .unwrap();
        let pixel = render_transfer_transport(&grid, &transfer, settings(2)).unwrap()[0];
        let background_expected = (300.0 * 10.0_f64.ln() - depth).exp();
        let source_expected = background_expected * -(-1.0_f64).exp_m1();
        for actual in pixel.emission_rgb() {
            close(actual, source_expected);
        }
        let foreground_values = [0.0];
        let foreground = VolumeGrid::new([1; 3], &foreground_values, [0.0; 3], [1.0; 3]);
        let pixel = render_transfer_transport(
            &foreground, &constant(depth, [0.0; 3]), settings(1),
        )
        .unwrap()[0];
        assert!(pixel.transmittance() < f64::MIN_POSITIVE);
        assert_eq!(pixel.optical_depth(), depth);
        for actual in pixel.composite([1e300; 3]).unwrap() {
            close(actual, background_expected);
        }
    }
}

#[test]
fn budgets_and_invalid_optics_refuse_without_an_image() {
    let values = [1.0, 2.0];
    let grid = VolumeGrid::new([1, 1, 2], &values, [0.0; 3], [1.0; 3]);
    let transfer = constant(1.0, [1.0; 3]);
    assert!(matches!(
        render_transfer_transport(
            &grid, &transfer, AbsorptionImageSettings { max_cell_visits: 1, ..settings(2) },
        ),
        Err(DvrError::CellVisitBudgetExceeded { .. })
    ));
    assert!(matches!(
        render_transfer_transport(
            &grid, &transfer, AbsorptionImageSettings { max_grid_cells: 1, ..settings(2) },
        ),
        Err(DvrError::GridCellBudgetExceeded { .. })
    ));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            render_density_transmittance(&grid, settings(2), calibration(invalid)),
            Err(DvrError::InvalidDensityAbsorption)
        );
    }
    let pixel = render_transfer_transport(&grid, &transfer, settings(2)).unwrap()[0];
    assert_eq!(
        pixel.composite([f64::NAN, 0.0, 0.0]),
        Err(DvrError::InvalidBackgroundRadiance)
    );
}
