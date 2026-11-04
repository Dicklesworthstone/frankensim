//! Continuous finite-patch heating through the original film/material owner.
use super::*;
use fs_em::plane_wave::angular_spectrum::AngularSpectrumAperture;

#[test]
fn g1_continuous_rectangle_heat_partitions_full_cell_and_keeps_local_sum_tones() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let stacks = [2.0, 2.3].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::DebyeWithOhmic {
                    dielectric: &dielectric,
                    conductor: &conductor,
                },
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let samples: Vec<_> = (0..4)
        .map(|x| {
            let sign = if x % 2 == 0 { 1.0 } else { -1.0 };
            [C64::ZERO, C64::new(0.8 * sign, 0.3 * sign)]
        })
        .collect();
    let lines = [prepare(&stacks[0], &samples), prepare(&stacks[1], &samples)];
    let banks = [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let layer = banks[0].prepare_layer(0, heat_limits(), || false).unwrap();
    let proper = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &banks,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let proper_layer = proper.prepare_layer(0, heat_limits(), || false).unwrap();
    let depth = simpson(0.31, 8);
    let base = Volume {
        depth_span_m: [0.0, 0.31],
        depth_samples: &depth,
        transverse: Areas::ContinuousPeriodicCell,
    };
    let zero_fft = AngularSpectrumIntegratedHeatLimits {
        maximum_fft_work: 0,
        ..integrated_limits()
    };
    let rectangle = |x, width| Volume {
        transverse: Areas::ContinuousRectangle(
            AngularSpectrumAperture::try_new([x, -0.2], [width, 3.0]).unwrap(),
        ),
        ..base
    };
    let full = layer
        .integrated_heat_spectrum(base, Channel::Total, zero_fft, || false)
        .unwrap();
    let whole = layer
        .integrated_heat_spectrum(rectangle(0.4, 32.0), Channel::Total, zero_fft, || false)
        .unwrap();
    let left = layer
        .integrated_heat_spectrum(rectangle(0.4, 7.0), Channel::Total, zero_fft, || false)
        .unwrap();
    let right = layer
        .integrated_heat_spectrum(rectangle(7.4, 25.0), Channel::Total, zero_fft, || false)
        .unwrap();
    let proper_left = proper_layer
        .integrated_heat_spectrum(rectangle(0.4, 7.0), Channel::Total, zero_fft, || false)
        .unwrap();
    let ohmic = layer
        .integrated_heat_spectrum(rectangle(0.4, 7.0), Channel::Ohmic, zero_fft, || false)
        .unwrap();
    let poles = layer
        .integrated_heat_spectrum(rectangle(0.4, 7.0), Channel::Polarization, zero_fft, || {
            false
        })
        .unwrap();
    assert!(core::ptr::eq(left.layer(), &layer));
    near(left.represented_area_m2_dd().to_f64(), 21.0, 1e-14);
    near(left.represented_volume_m3_dd().to_f64(), 21.0 * 0.31, 1e-14);
    let mut local_sum_effect = 0.0_f64;
    for time in [0.0, 0.17, 0.73, 1.3] {
        let expected = full.evaluate(time).unwrap();
        near(whole.evaluate(time).unwrap(), expected, 3e-11);
        near(
            left.evaluate(time).unwrap() + right.evaluate(time).unwrap(),
            expected,
            3e-11,
        );
        near(
            left.evaluate(time).unwrap(),
            ohmic.evaluate(time).unwrap() + poles.evaluate(time).unwrap(),
            3e-11,
        );
        local_sum_effect = local_sum_effect
            .max((left.evaluate(time).unwrap() - proper_left.evaluate(time).unwrap()).abs());
    }
    assert!(
        local_sum_effect > 1e-3,
        "finite patch lost physical real-field sum tones"
    );
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        near(
            whole.observe_window(window).unwrap(),
            full.observe_window(window).unwrap(),
            3e-11,
        );
        near(
            left.observe_window(window).unwrap() + right.observe_window(window).unwrap(),
            full.observe_window(window).unwrap(),
            3e-11,
        );
    }
}

#[test]
fn g4_continuous_rectangle_admits_mode_pair_work_and_empty_banks_without_cache() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response: DielectricResponse::DebyeWithOhmic {
                dielectric: &dielectric,
                conductor: &conductor,
            },
            thickness_m: 0.31,
        }],
    )
    .unwrap();
    let lines = [constant(&stack, [C64::new(0.8, 0.2), C64::ZERO])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let depth = simpson(0.31, 2);
    let volume = Volume {
        depth_span_m: [0.0, 0.31],
        depth_samples: &depth,
        transverse: Areas::ContinuousRectangle(
            AngularSpectrumAperture::try_new([0.3, -0.4], [2.0, 0.7]).unwrap(),
        ),
    };
    let exact = AngularSpectrumIntegratedHeatLimits {
        maximum_fft_work: 0,
        maximum_spatial_point_frequency_pairs: depth.len() * (4 * 4 + 1),
        ..integrated_limits()
    };
    let heat = layer
        .integrated_heat_spectrum(volume, Channel::Total, exact, || false)
        .unwrap();
    assert!(heat.evaluate(0.17).unwrap() > 0.0);
    assert!(
        layer
            .integrated_heat_spectrum(
                volume,
                Channel::Total,
                AngularSpectrumIntegratedHeatLimits {
                    maximum_spatial_point_frequency_pairs: exact
                        .maximum_spatial_point_frequency_pairs
                        - 1,
                    ..exact
                },
                || false
            )
            .is_err()
    );
    assert!(
        layer
            .integrated_heat_spectrum(
                volume,
                Channel::Total,
                AngularSpectrumIntegratedHeatLimits {
                    maximum_workspace_bytes: heat.preparation_workspace_bytes() - 1,
                    ..exact
                },
                || false
            )
            .is_err()
    );
    assert!(matches!(
        layer.integrated_heat_spectrum(volume, Channel::Total, exact, || true),
        Err(AngularSpectrumError::Cancelled)
    ));

    let empty = AngularSpectrumExposure::try_new(&lines[0], &[], limits(), || false).unwrap();
    let dark = empty.prepare_layer(0, heat_limits(), || false).unwrap();
    let no_work = AngularSpectrumIntegratedHeatLimits {
        maximum_point_line_pairs: 0,
        maximum_fft_work: 0,
        maximum_spatial_point_frequency_pairs: 0,
        maximum_reservoir_frequency_pairs: 0,
        maximum_generated_terms: 0,
        ..integrated_limits()
    };
    let cell = dark
        .integrated_heat_spectrum(
            Volume {
                transverse: Areas::ContinuousPeriodicCell,
                ..volume
            },
            Channel::Total,
            no_work,
            || false,
        )
        .unwrap();
    // Equal workspace to darkness on the full cell proves no rectangular
    // W²+H² cache is prepared when no modal-frequency products can exist.
    let patch = dark
        .integrated_heat_spectrum(
            volume,
            Channel::Total,
            AngularSpectrumIntegratedHeatLimits {
                maximum_workspace_bytes: cell.preparation_workspace_bytes(),
                ..no_work
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        patch.preparation_workspace_bytes(),
        cell.preparation_workspace_bytes()
    );
    assert_eq!(patch.evaluate(0.17).unwrap(), 0.0);
}
