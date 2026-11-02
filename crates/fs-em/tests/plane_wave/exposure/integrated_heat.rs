//! G1/G3/G4 physically integrated loads through the original film owner.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureHeatChannel as Channel, AngularSpectrumHeatAreaQuadrature as Areas,
    AngularSpectrumHeatDepthSample as DepthSample, AngularSpectrumHeatVolumeQuadrature as Volume,
    AngularSpectrumIntegratedHeatLimits,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn integrated_limits() -> AngularSpectrumIntegratedHeatLimits {
    AngularSpectrumIntegratedHeatLimits {
        maximum_depth_samples: 128,
        maximum_point_line_pairs: 1 << 20,
        maximum_fft_work: 1 << 24,
        maximum_spatial_point_frequency_pairs: 1 << 20,
        maximum_reservoir_frequency_pairs: 1 << 16,
        maximum_generated_terms: 8192,
        maximum_workspace_bytes: 64 << 20,
    }
}
fn simpson(thickness: f64, intervals: usize) -> Vec<DepthSample> {
    assert_eq!(intervals % 2, 0);
    (0..=intervals)
        .map(|i| DepthSample {
            depth_m: if i == intervals {
                thickness
            } else {
                thickness * i as f64 / intervals as f64
            },
            weight_m: thickness / (3.0 * intervals as f64)
                * if i == 0 || i == intervals {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                },
        })
        .collect()
}

#[test]
fn g1_physical_mask_retains_beats_that_cancel_over_the_original_periodic_cell() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let stacks = [2.0, 2.3].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let oblique: Vec<_> = (0..4)
        .map(|x| {
            let phase = -(x as f64) * core::f64::consts::FRAC_PI_2;
            let z = C64::new(phase.cos(), phase.sin());
            [C64::new(-0.3, 0.6) * z, C64::new(0.1, -0.3) * z]
        })
        .collect();
    let lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        prepare(&stacks[1], &oblique),
    ];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let layer = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
    let samples = simpson(0.31, 8);
    let complete = Volume {
        depth_span_m: [0.0, 0.31],
        depth_samples: &samples,
        transverse: Areas::PeriodicCell,
    };
    let areas = [24.0, 12.0, 0.0, 24.0];
    let partial = Volume {
        transverse: Areas::CellAreasM2(&areas),
        ..complete
    };
    let full = layer
        .integrated_heat_spectrum(complete, Channel::Total, integrated_limits(), || false)
        .unwrap();
    let masked = layer
        .integrated_heat_spectrum(partial, Channel::Total, integrated_limits(), || false)
        .unwrap();
    assert!(core::ptr::eq(masked.layer(), &layer));
    assert_eq!(masked.grid(), grid());
    assert_eq!(masked.channel(), Channel::Total);
    assert_eq!(masked.quadrature().depth_samples.as_ptr(), samples.as_ptr());
    near(full.represented_area_m2_dd().to_f64(), 96.0, 2e-14);
    near(masked.represented_area_m2_dd().to_f64(), 60.0, 2e-14);
    near(
        masked.represented_volume_m3_dd().to_f64(),
        60.0 * 0.31,
        2e-13,
    );

    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let mut full_reference = 0.0;
        let mut masked_reference = 0.0;
        for sample in &samples {
            let field = layer.field(sample.depth_m, window, || false).unwrap();
            full_reference += sample.weight_m * field.absorbed_power_per_depth_w_m();
            masked_reference += sample.weight_m
                * field
                    .dissipation()
                    .iter()
                    .zip(areas)
                    .map(|(q, area)| q.total_heat_w_m3 * area)
                    .sum::<f64>();
        }
        near(full.observe_window(window).unwrap(), full_reference, 3e-11);
        near(
            masked.observe_window(window).unwrap(),
            masked_reference,
            3e-11,
        );
    }
    let full_start = full.evaluate(0.0).unwrap();
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for time in [0.0, 1.0, 2.0, 3.0, 4.0] {
        // Different transverse modes are orthogonal on the full periodic
        // cell. Their cross-frequency beat survives on the physical mask.
        near(full.evaluate(time).unwrap(), full_start, 3e-11);
        let value = masked.evaluate(time).unwrap();
        minimum = minimum.min(value);
        maximum = maximum.max(value);
    }
    assert!(maximum - minimum > 1e-3);
    let zero = [0.0; 4];
    let dark = layer
        .integrated_heat_spectrum(
            Volume {
                transverse: Areas::CellAreasM2(&zero),
                ..complete
            },
            Channel::Total,
            integrated_limits(),
            || false,
        )
        .unwrap();
    assert_eq!(dark.represented_volume_m3_dd().to_f64(), 0.0);
    assert_eq!(dark.evaluate(0.7).unwrap(), 0.0);
}

#[test]
fn g3_depth_refinement_matches_original_absorbed_watts_and_owner_survives_refusal() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let thickness = 1.2;
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let lines = [constant(&stack, [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let cycle = HarmonicWindow::new(0.0, core::f64::consts::PI).unwrap();
    let original = lines[0]
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 4096,
                maximum_workspace_bytes: 16 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap();
    let expected = original.absorbed_power_w();
    let mut errors = Vec::new();
    for intervals in [8, 32, 64] {
        let samples = simpson(thickness, intervals);
        let quadrature = Volume {
            depth_span_m: [0.0, thickness],
            depth_samples: &samples,
            transverse: Areas::PeriodicCell,
        };
        let heat = layer
            .integrated_heat_spectrum(quadrature, Channel::Total, integrated_limits(), || false)
            .unwrap();
        errors.push((heat.observe_window(cycle).unwrap() - expected).abs());
    }
    assert!(errors[1] < 0.03 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.15 * errors[1], "{errors:?}");
    assert!(errors[2] < 1e-5 * (1.0 + expected.abs()), "{errors:?}");

    let samples = simpson(thickness, 8);
    let quadrature = Volume {
        depth_span_m: [0.0, thickness],
        depth_samples: &samples,
        transverse: Areas::PeriodicCell,
    };
    let heat = layer
        .integrated_heat_spectrum(quadrature, Channel::Total, integrated_limits(), || false)
        .unwrap();
    let value = heat.evaluate(0.17).unwrap();
    let exact = AngularSpectrumIntegratedHeatLimits {
        maximum_depth_samples: samples.len(),
        maximum_point_line_pairs: layer.point_line_pairs() * samples.len(),
        maximum_fft_work: layer.field_fft_work() * samples.len() as u64,
        maximum_spatial_point_frequency_pairs: 4 * samples.len(),
        maximum_reservoir_frequency_pairs: 2,
        maximum_generated_terms: 2,
        maximum_workspace_bytes: heat.preparation_workspace_bytes(),
    };
    let replay = layer
        .integrated_heat_spectrum(quadrature, Channel::Total, exact, || false)
        .unwrap();
    assert_eq!(replay.evaluate(0.17).unwrap(), value);
    for limits in [
        AngularSpectrumIntegratedHeatLimits {
            maximum_depth_samples: samples.len() - 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_point_line_pairs: exact.maximum_point_line_pairs - 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_fft_work: exact.maximum_fft_work - 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_spatial_point_frequency_pairs: 4 * samples.len() - 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_reservoir_frequency_pairs: 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_generated_terms: 1,
            ..exact
        },
        AngularSpectrumIntegratedHeatLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            ..exact
        },
    ] {
        assert!(
            layer
                .integrated_heat_spectrum(quadrature, Channel::Total, limits, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        layer.integrated_heat_spectrum(quadrature, Channel::Total, exact, || {
            polls += 1;
            polls == 15
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        heat.evaluate_cancellable(0.17, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        heat.observe_window_cancellable(cycle, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(heat.evaluate(0.17).unwrap(), value);
    assert_eq!(
        layer
            .integrated_heat_spectrum(quadrature, Channel::Total, exact, || false)
            .unwrap()
            .evaluate(0.17)
            .unwrap(),
        value
    );
}
