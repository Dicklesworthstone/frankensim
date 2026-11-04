//! Original public field/material/cell-solve integration, with authored laws.
use super::*;
use fs_couple::layered_wave_heating::spatial_interference::{
    SpatialInterferenceCellQuadrature, SpatialInterferenceHeating,
    SpatialInterferenceHeatingControls, SpatialInterferenceHeatingLimits,
};
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposure, AngularSpectrumExposureEnsemble,
    AngularSpectrumExposureEnsembleLimits, AngularSpectrumExposureLayerLimits,
    AngularSpectrumExposureLimits, AngularSpectrumGrid, AngularSpectrumHeatDepthSample,
    AngularSpectrumIntegratedHeatLimits, AngularSpectrumLimits, PlanarAngularSpectrum,
};
use fs_math::dd::Dd;

fn limits() -> SpatialInterferenceHeatingLimits {
    SpatialInterferenceHeatingLimits {
        maximum_cells: 4,
        maximum_preparation_workspace_bytes: 64 << 20,
        layer: AngularSpectrumExposureLayerLimits {
            maximum_frequency_groups: 16,
            maximum_frequency_layer_mode_pairs: 1024,
            maximum_point_line_pairs: 1024,
            maximum_point_reservoir_frequency_pairs: 8192,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        cell_heat: AngularSpectrumIntegratedHeatLimits {
            maximum_depth_samples: 33,
            maximum_point_line_pairs: 1 << 16,
            maximum_fft_work: 0,
            maximum_spatial_point_frequency_pairs: 1 << 20,
            maximum_reservoir_frequency_pairs: 1 << 16,
            maximum_generated_terms: 1024,
            maximum_workspace_bytes: 16 << 20,
        },
    }
}
fn depth(width: f64) -> Vec<AngularSpectrumHeatDepthSample> {
    (0..=16)
        .map(|i| AngularSpectrumHeatDepthSample {
            depth_m: if i == 16 {
                width
            } else {
                width * i as f64 / 16.0
            },
            weight_m: width / 48.0
                * if i == 0 || i == 16 {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                },
        })
        .collect()
}

/// G1/G3/G4: oblique cross-frequency patch interference becomes native spatial
/// cell energy. The thermal audit retains authored heat and internal conduction;
/// an independent time quadrature checks physical source integration.
#[test]
fn finite_patch_beats_reach_native_depth_cells_once_and_keep_the_accepted_clock() {
    let exterior = state(&card("spatial-exterior", 1.0, 0.0), false, 100_000.0);
    let film_card = card("spatial-film", 2.7, 0.6);
    let model = debye(&film_card, 100_000.0);
    let width = 0.6 / k();
    let area = 0.2;
    let native = body(&[&film_card], &[width], &[2], area);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &[DielectricResponse::Debye(&model)],
        || false,
    )
    .unwrap();
    let films: Vec<_> = (0..2)
        .map(|i| owner.optical_layers().layer(i).unwrap())
        .collect();
    let other_stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA + 50.0,
        &films,
    )
    .unwrap();
    let normal = [[C64::ZERO, C64::from_re(100.0)]; 4];
    let oblique: Vec<_> = (0..4)
        .map(|i| {
            let phase = -(i as f64) * core::f64::consts::FRAC_PI_2;
            [C64::ZERO, C64::new(60.0 * phase.cos(), 60.0 * phase.sin())]
        })
        .collect();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [1.0, 1.0],
    };
    let angular = AngularSpectrumLimits {
        maximum_grid_points: 16,
        maximum_workspace_bytes: 8 << 20,
        maximum_fft_work: 1 << 20,
    };
    let beams = [
        PlanarAngularSpectrum::prepare(owner.optical_layers(), grid, &normal, angular, || false)
            .unwrap(),
        PlanarAngularSpectrum::prepare(&other_stack, grid, &oblique, angular, || false).unwrap(),
    ];
    let factors = [AngularSpectrumExposure::try_new(
        &beams[0],
        &beams,
        AngularSpectrumExposureLimits {
            maximum_lines: 4,
            maximum_point_line_pairs: 1024,
            maximum_point_frequency_pairs: 1024,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &beams[0],
        None,
        &factors,
        AngularSpectrumExposureEnsembleLimits {
            maximum_factors: 2,
            maximum_point_line_pairs: 1024,
            maximum_point_frequency_pairs: 1024,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    let samples = depth(width / 2.0);
    let rules = [SpatialInterferenceCellQuadrature {
        depth_samples: &samples,
    }; 2];
    let footprint = AngularSpectrumAperture::try_new([0.15, -0.2], [0.5, 0.4]).unwrap();
    let heating = SpatialInterferenceHeating::try_new(
        &owner,
        &ensemble,
        &rules,
        footprint,
        "laboratory clock",
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(heating.cell_spectra().len(), 2);
    assert!(core::ptr::eq(heating.original_owner(), &owner));
    assert!(core::ptr::eq(heating.original_ensemble(), &ensemble));
    assert_ne!(
        heating.cell_spectra()[0].evaluate(0.0).unwrap(),
        heating.cell_spectra()[1].evaluate(0.0).unwrap()
    );

    // A translated patch of the SAME area sees a different cross-mode beat.
    let shifted = SpatialInterferenceHeating::try_new(
        &owner,
        &ensemble,
        &rules,
        AngularSpectrumAperture::try_new([1.15, -0.2], [0.5, 0.4]).unwrap(),
        "laboratory clock",
        limits(),
        || false,
    )
    .unwrap();
    assert!(
        (heating.cell_spectra()[0].evaluate(0.017).unwrap()
            - shifted.cell_spectra()[0].evaluate(0.017).unwrap())
        .abs()
            > 1e-4
    );
    // An area mismatch is a geometry error, never a request to spread heat.
    assert!(
        SpatialInterferenceHeating::try_new(
            &owner,
            &ensemble,
            &rules,
            AngularSpectrumAperture::try_new([0.0; 2], [1.0, 0.4]).unwrap(),
            "laboratory clock",
            limits(),
            || false
        )
        .is_err()
    );
    assert!(
        SpatialInterferenceHeating::try_new(
            &owner,
            &ensemble,
            &rules,
            footprint,
            "laboratory clock",
            SpatialInterferenceHeatingLimits {
                maximum_cells: 1,
                ..limits()
            },
            || false
        )
        .is_err()
    );

    let authored = [0.5, 0.0];
    let accepted = heating
        .initial_state(initial(&bound, &authored), Dd::from_f64(0.017), || false)
        .unwrap();
    let before = accepted.thermal().cells().clone();
    let temporal = SpatialInterferenceHeatingControls {
        maximum_thermal_step_s: 0.01,
        maximum_source_ripple_j: 1.0,
    };
    let dt = 0.002;
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let step = heating
        .advance_to(
            &accepted,
            dt,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    for (i, spectrum) in heating.cell_spectra().iter().enumerate() {
        let energy = (0..=256)
            .map(|j| {
                let weight = if j == 0 || j == 256 {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * spectrum.evaluate(0.017 + dt * j as f64 / 256.0).unwrap()
            })
            .sum::<f64>()
            * dt
            / (3.0 * 256.0);
        near(step.cell_energy_j[i].1, energy, 2e-13);
        assert_eq!(step.cell_energy_j[i].0, owner.cell_ids()[i]);
        assert_eq!(
            step.state.thermal().cells()[&owner.cell_ids()[i]].boundary,
            accepted.thermal().cells()[&owner.cell_ids()[i]].boundary
        );
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.absorbed_energy_j.to_f64() + authored.iter().sum::<f64>() * dt,
        1e-8,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-8,
    );
    near(
        (step.state.source_time_s() - accepted.source_time_s()).to_f64(),
        dt,
        1e-18,
    );
    assert_eq!(step.state.thermal().time_s(), dt);
    let next = heating
        .advance_to(
            &step.state,
            2.0 * dt,
            temporal,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    near(
        (next.state.source_time_s() - accepted.source_time_s()).to_f64(),
        2.0 * dt,
        1e-18,
    );
    assert_ne!(next.cell_energy_j, step.cell_energy_j);
    assert!(
        heating
            .advance_to(
                &accepted,
                dt,
                SpatialInterferenceHeatingControls {
                    maximum_source_ripple_j: 0.0,
                    ..temporal
                },
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || false
            )
            .is_err()
    );
    assert!(
        heating
            .advance_to(
                &accepted,
                dt,
                temporal,
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || true
            )
            .is_err()
    );
    assert_eq!(accepted.thermal().time_s(), 0.0);
    assert_eq!(accepted.thermal().cells(), &before);
}
