use super::*;

fn settings(cells: usize) -> SpectralBeamDepositionSettings {
    SpectralBeamDepositionSettings {
        max_grid_cells: cells,
        max_spectral_evaluations: 16 * cells as u64,
        max_output_values: 2 * cells,
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
        (actual - expected).abs() <= expected.abs().max(1e-300) * 4e-12,
        "{actual:e} != {expected:e}"
    );
}

fn balance(result: &SpectralBeamDeposition) {
    for lane in 0..result.wavelengths_nm().len() {
        close(
            result.incident_power_w()[lane],
            result.absorbed_power_w()[lane] + result.transmitted_power_w()[lane],
        );
    }
    close(
        result.total_absorbed_power_w(),
        result.cell_absorbed_power_w().iter().sum(),
    );
    close(
        result.total_transmitted_power_w(),
        result.column_transmitted_power_w().iter().sum(),
    );
    close(
        result.total_incident_power_w(),
        result.total_absorbed_power_w() + result.total_transmitted_power_w(),
    );
}

#[test]
fn g1_spectral_slabs_retain_si_power_cell_order_and_explicit_quadrature() {
    let grid = VolumeGrid::new([1, 1, 2], &[2.0, 4.0], [0.0; 3], [2.0, 3.0, 0.4]);
    let calibration = DensityAbsorptionSpectrum {
        density_unit_kg_m3: 5.0,
        length_unit_m: 0.5,
        ..spectrum(&[500.0, 10000.0], &[0.25, 0.5])
    };
    for direction in [BeamPropagation::PositiveZ, BeamPropagation::NegativeZ] {
        let result = deposit_density_spectral_beam(
            &grid,
            settings(2),
            calibration,
            SpectralBeamIrradiance::Uniform(&[4.0, 10.0]),
            &[2.0, 3.0],
            direction,
        )
        .unwrap();
        assert_eq!(result.dimensions(), [1, 1, 2]);
        assert_eq!(result.wavelengths_nm(), &[500.0, 10000.0]);
        assert_eq!(result.quadrature_weights_nm(), &[2.0, 3.0]);
        assert_eq!(result.propagation(), direction);
        // Column area is (2*0.5)*(3*0.5)=1.5 m².
        assert_eq!(result.incident_power_w(), &[12.0, 45.0]);
        let mut expected_cells = [0.0; 2];
        for (lane, power) in [12.0, 45.0].into_iter().enumerate() {
            let depth = [[0.5_f64, 1.0], [1.0, 2.0]][lane];
            let front = usize::from(direction == BeamPropagation::NegativeZ);
            let back = 1 - front;
            expected_cells[front] += power * -(-depth[front]).exp_m1();
            expected_cells[back] += power * (-depth[front]).exp() * -(-depth[back]).exp_m1();
            close(
                result.transmitted_power_w()[lane],
                power * (-depth[0] - depth[1]).exp(),
            );
        }
        for (actual, expected) in result.cell_absorbed_power_w().iter().zip(expected_cells) {
            close(*actual, expected);
        }
        balance(&result);
        let totals = [
            result.total_incident_power_w(),
            result.total_absorbed_power_w(),
            result.total_transmitted_power_w(),
        ];
        let flat = result.into_power_buffer();
        assert_eq!(flat.len(), 6);
        assert_eq!(&flat[..3], &totals);
    }
}

#[test]
fn g1_indexed_layers_and_spatial_irradiance_match_independent_beer_power_integrals() {
    let density = [1.0, 2.0, 3.0, 4.0];
    let indices = [0, 1, 1, 0];
    let coefficients = [0.25, 1.0, 0.5, 0.0];
    let source = [2.0, 10.0, 20.0, 4.0];
    let weights = [2.0, 5.0];
    let grid = VolumeGrid::new([2, 1, 2], &density, [0.0; 3], [1.0; 3]);
    let result = deposit_density_spectral_indexed_beam(
        &grid,
        IndexedSpectralBeamDepositionSettings {
            deposition: settings(4),
            max_materials: 2,
        },
        IndexedDensityAbsorptionSpectrum {
            wavelengths_nm: &[450.0, 900.0],
            mass_absorption_m2_kg: &coefficients,
            cell_materials: &indices,
            density_unit_kg_m3: 1.0,
            length_unit_m: 1.0,
        },
        SpectralBeamIrradiance::Columns(&source),
        &weights,
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    let mut expected_cells = [0.0; 4];
    for column in 0..2 {
        let mut transmitted = 0.0;
        for lane in 0..2 {
            let power = source[column * 2 + lane] * weights[lane];
            let a = density[column] * coefficients[indices[column] as usize * 2 + lane];
            let b = density[column + 2] * coefficients[indices[column + 2] as usize * 2 + lane];
            expected_cells[column] += power * -(-a).exp_m1();
            expected_cells[column + 2] += power * (-a).exp() * -(-b).exp_m1();
            transmitted += power * (-a - b).exp();
        }
        close(result.column_transmitted_power_w()[column], transmitted);
    }
    for (actual, expected) in result.cell_absorbed_power_w().iter().zip(expected_cells) {
        close(*actual, expected);
    }
    balance(&result);
}

#[test]
fn g3_fixed_irradiance_preserves_total_power_under_spatial_refinement_and_unit_changes() {
    let coarse = VolumeGrid::new([1; 3], &[5.0], [0.0; 3], [1.0; 3]);
    let fine = VolumeGrid::new([2, 2, 4], &[5.0; 16], [0.0; 3], [0.5, 0.5, 0.25]);
    let calibration = spectrum(&[550.0], &[0.1]);
    let first = deposit_density_spectral_beam(
        &coarse,
        settings(1),
        calibration,
        SpectralBeamIrradiance::Uniform(&[10.0]),
        &[2.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    let second = deposit_density_spectral_beam(
        &fine,
        settings(16),
        calibration,
        SpectralBeamIrradiance::Columns(&[10.0; 4]),
        &[2.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    assert_eq!(first.total_incident_power_w(), 20.0);
    assert_eq!(second.total_incident_power_w(), 20.0);
    close(
        first.total_absorbed_power_w(),
        second.total_absorbed_power_w(),
    );
    close(
        first.total_transmitted_power_w(),
        second.total_transmitted_power_w(),
    );
    balance(&second);
    let changed_units = VolumeGrid::new([1; 3], &[1.0], [1e12, -1e12, 0.0], [2.0; 3]);
    let rescaled = deposit_density_spectral_beam(
        &changed_units,
        settings(1),
        DensityAbsorptionSpectrum {
            density_unit_kg_m3: 5.0,
            length_unit_m: 0.5,
            ..calibration
        },
        SpectralBeamIrradiance::Uniform(&[10.0]),
        &[2.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    assert_eq!(
        first.cell_absorbed_power_w(),
        rescaled.cell_absorbed_power_w()
    );
    assert_eq!(first.incident_power_w(), rescaled.incident_power_w());
}

#[test]
fn g1_cumulative_depth_resolves_many_individually_sub_ulp_thin_layers() {
    let density = vec![1.0; 10_000];
    let grid = VolumeGrid::new([1, 1, density.len()], &density, [0.0; 3], [1.0; 3]);
    let result = deposit_density_spectral_beam(
        &grid,
        settings(density.len()),
        spectrum(&[550.0], &[1e-18]),
        SpectralBeamIrradiance::Uniform(&[1e10]),
        &[1.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    assert_eq!((-1e-18_f64).exp(), 1.0);
    close(
        result.total_absorbed_power_w(),
        1e10 * -(-1e-14_f64).exp_m1(),
    );
    close(
        result.total_transmitted_power_w(),
        1e10 * (-1e-14_f64).exp(),
    );
    assert!(result.total_transmitted_power_w() < 1e10);
    assert!(
        result
            .cell_absorbed_power_w()
            .iter()
            .all(|power| *power > 0.0)
    );
    balance(&result);
}

#[test]
fn g1_physical_power_preserves_thin_absorption_bright_opaque_survival_and_small_area() {
    let thin = VolumeGrid::new([1; 3], &[1e-200], [0.0; 3], [1.0; 3]);
    let result = deposit_density_spectral_beam(
        &thin,
        settings(1),
        spectrum(&[550.0], &[1e-200]),
        SpectralBeamIrradiance::Uniform(&[1e300]),
        &[1.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    close(result.total_absorbed_power_w(), 1e-100);
    let opaque = VolumeGrid::new([1; 3], &[800.0], [0.0; 3], [1.0; 3]);
    let result = deposit_density_spectral_beam(
        &opaque,
        settings(1),
        spectrum(&[550.0], &[1.0]),
        SpectralBeamIrradiance::Uniform(&[1e300]),
        &[1.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    close(
        result.total_transmitted_power_w(),
        (1e300_f64.ln() - 800.0).exp(),
    );
    assert!(result.total_transmitted_power_w() > 0.0);
    balance(&result);
    let tiny_area = VolumeGrid::new([1; 3], &[1.0], [0.0; 3], [1e-200, 1e-200, 1.0]);
    let result = deposit_density_spectral_beam(
        &tiny_area,
        settings(1),
        spectrum(&[550.0], &[core::f64::consts::LN_2]),
        SpectralBeamIrradiance::Uniform(&[1e300]),
        &[1e100],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    close(result.total_incident_power_w(), 1.0);
    close(result.total_absorbed_power_w(), 0.5);
    close(result.total_transmitted_power_w(), 0.5);
}

#[test]
fn g1_partial_coverage_absorbs_only_the_actual_correlated_incident_power() {
    let grid = VolumeGrid::new([2, 1, 1], &[0.0, 2.0], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[550.0], &[core::f64::consts::LN_2]);
    let result = deposit_density_spectral_beam(
        &grid,
        settings(2),
        calibration,
        SpectralBeamIrradiance::Columns(&[2.0, 20.0]),
        &[1.0],
        BeamPropagation::PositiveZ,
    )
    .unwrap();
    close(result.total_incident_power_w(), 22.0);
    close(result.total_absorbed_power_w(), 15.0);
    close(result.total_transmitted_power_w(), 7.0);
    assert_eq!(result.cell_absorbed_power_w()[0], 0.0);
    // Zero quadrature contributes no power but retains validated source nodes.
    let zero = deposit_density_spectral_beam(
        &grid,
        settings(2),
        calibration,
        SpectralBeamIrradiance::Columns(&[2.0, 20.0]),
        &[0.0],
        BeamPropagation::NegativeZ,
    )
    .unwrap();
    assert_eq!(zero.into_power_buffer(), vec![0.0; 7]);
}

#[test]
fn g0_power_sources_work_limits_and_original_material_map_are_admitted_atomically() {
    let grid = VolumeGrid::new([2, 1, 1], &[0.0, 1.0], [0.0; 3], [1.0; 3]);
    let calibration = spectrum(&[550.0], &[0.0]);
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            deposit_density_spectral_beam(
                &grid,
                settings(2),
                calibration,
                SpectralBeamIrradiance::Columns(&[0.0, invalid]),
                &[0.0],
                BeamPropagation::PositiveZ
            ),
            Err(DvrError::InvalidBeamIrradiance)
        );
        assert_eq!(
            deposit_density_spectral_beam(
                &grid,
                settings(2),
                calibration,
                SpectralBeamIrradiance::Uniform(&[0.0]),
                &[invalid],
                BeamPropagation::PositiveZ
            ),
            Err(DvrError::InvalidSpectralWeights)
        );
    }
    assert_eq!(
        deposit_density_spectral_beam(
            &grid,
            settings(2),
            calibration,
            SpectralBeamIrradiance::Columns(&[1.0]),
            &[1.0],
            BeamPropagation::PositiveZ
        ),
        Err(DvrError::InvalidBeamIrradianceSize)
    );
    let mut limit = settings(2);
    limit.max_output_values = 3;
    assert_eq!(
        deposit_density_spectral_beam(
            &grid,
            limit,
            calibration,
            SpectralBeamIrradiance::Columns(&[f64::NAN; 2]),
            &[1.0],
            BeamPropagation::PositiveZ
        ),
        Err(DvrError::BeamOutputBudgetExceeded {
            requested: 4,
            limit: 3
        })
    );
    limit = settings(2);
    limit.max_spectral_evaluations = 1;
    assert_eq!(
        deposit_density_spectral_beam(
            &grid,
            limit,
            calibration,
            SpectralBeamIrradiance::Columns(&[f64::NAN; 2]),
            &[1.0],
            BeamPropagation::PositiveZ
        ),
        Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: 2,
            limit: 1
        })
    );
    assert_eq!(
        deposit_density_spectral_indexed_beam(
            &grid,
            IndexedSpectralBeamDepositionSettings {
                deposition: settings(2),
                max_materials: 1
            },
            IndexedDensityAbsorptionSpectrum {
                wavelengths_nm: &[550.0],
                mass_absorption_m2_kg: &[0.0],
                cell_materials: &[0, 1],
                density_unit_kg_m3: 1.0,
                length_unit_m: 1.0,
            },
            SpectralBeamIrradiance::Uniform(&[0.0]),
            &[0.0],
            BeamPropagation::PositiveZ
        ),
        Err(DvrError::AbsorptionMaterialIndexOutOfBounds {
            cell: 1,
            material: 1,
            materials: 1
        })
    );
    let one = VolumeGrid::new([1; 3], &[0.0], [0.0; 3], [1.0; 3]);
    assert_eq!(
        deposit_density_spectral_beam(
            &one,
            settings(1),
            calibration,
            SpectralBeamIrradiance::Uniform(&[f64::MAX]),
            &[2.0],
            BeamPropagation::PositiveZ
        ),
        Err(DvrError::UnrepresentableBeamPower)
    );
}
