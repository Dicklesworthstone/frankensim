//! Carrier-aware irreversible heat through the actual dispersive film owner.
use super::*;
use fs_em::plane_wave::angular_spectrum::AngularSpectrumAperture;

fn rotation(angle: f64) -> C64 {
    let (sine, cosine) = angle.sin_cos();
    C64::new(cosine, sine)
}
fn sinc(value: f64) -> f64 {
    if value == 0.0 {
        1.0
    } else {
        value.sin() / value
    }
}

#[test]
fn g1_carrier_heat_integrates_the_full_cell_and_finite_patches_with_physical_sum_tones() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let omega = 2.0;
    let carrier = [0.23, -0.17];
    let stack = PlanarDielectricLayers::try_new(
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
    .unwrap();
    let lines = [PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid(),
        carrier,
        &[[C64::ZERO, C64::new(0.8, 0.3)]; 4],
        AngularSpectrumLimits {
            maximum_grid_points: 4,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()];
    let banks = [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let layer = banks[0].prepare_layer(0, heat_limits(), || false).unwrap();
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &banks,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let proper_layer = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
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
    let depth = simpson(0.31, 4);
    let electric: Vec<_> = depth
        .iter()
        .map(|sample| {
            original
                .field(sample.depth_m, || false)
                .unwrap()
                .field()
                .electric_v_per_m()[0]
        })
        .collect();
    let volume = Volume {
        depth_span_m: [0.0, 0.31],
        depth_samples: &depth,
        transverse: Areas::ContinuousPeriodicCell,
    };
    let exact = AngularSpectrumIntegratedHeatLimits {
        maximum_fft_work: 0,
        maximum_spatial_point_frequency_pairs: depth.len() * (4 * 4 + 1),
        ..integrated_limits()
    };
    let full = layer
        .integrated_heat_spectrum(volume, Channel::Total, exact, || false)
        .unwrap();
    let whole_area = AngularSpectrumAperture::try_new([0.0; 2], [32.0, 3.0]).unwrap();
    let left_area = AngularSpectrumAperture::try_new([0.0; 2], [7.0, 3.0]).unwrap();
    let right_area = AngularSpectrumAperture::try_new([7.0, 0.0], [25.0, 3.0]).unwrap();
    let patch = |aperture| {
        layer
            .integrated_heat_spectrum(
                Volume {
                    transverse: Areas::ContinuousRectangle(aperture),
                    ..volume
                },
                Channel::Total,
                exact,
                || false,
            )
            .unwrap()
    };
    let whole = patch(whole_area);
    let left = patch(left_area);
    let right = patch(right_area);
    let proper = proper_layer
        .integrated_heat_spectrum(
            volume,
            Channel::Total,
            AngularSpectrumIntegratedHeatLimits {
                maximum_spatial_point_frequency_pairs: depth.len() * 4,
                ..exact
            },
            || false,
        )
        .unwrap();

    // Independent analytic area integral of Re(E*exp(i*omega*t-i*k.r))².
    // E comes from the original internal Maxwell field at x=y=0. The Debye
    // current uses J=i*omega*delta*E/(1+i*omega*tau), with heat tau/delta*J².
    let reference = |aperture: AngularSpectrumAperture, temporal_sum: C64, proper: bool| {
        let origin = aperture.origin_m();
        let extent = aperture.extent_m();
        let phase = -2.0
            * (carrier[0] * (origin[0] + 0.5 * extent[0])
                + carrier[1] * (origin[1] + 0.5 * extent[1]));
        let spatial_sum =
            rotation(phase).scale(sinc(carrier[0] * extent[0]) * sinc(carrier[1] * extent[1]));
        let mean_square = |values: [C64; 3]| {
            let mut difference = 0.0;
            let mut sum = C64::ZERO;
            for value in values {
                difference += value.norm_sq();
                sum = sum + value * value;
            }
            0.5 * difference
                + if proper {
                    0.0
                } else {
                    0.5 * (sum * spatial_sum * temporal_sum).re
                }
        };
        let current = C64::new(0.0, omega * 0.7) / C64::new(1.0, omega * 0.4);
        depth
            .iter()
            .zip(&electric)
            .map(|(sample, &e)| {
                sample.weight_m
                    * (0.3 * mean_square(e) + (0.4 / 0.7) * mean_square(e.map(|v| v * current)))
            })
            .sum::<f64>()
            * aperture.area_m2()
    };
    let mut real_sum_effect = 0.0_f64;
    for time in [0.0, 0.17, 0.73, 1.3] {
        let temporal = rotation(2.0 * omega * time);
        near(
            full.evaluate(time).unwrap(),
            reference(whole_area, temporal, false),
            3e-12,
        );
        near(
            whole.evaluate(time).unwrap(),
            full.evaluate(time).unwrap(),
            3e-12,
        );
        near(
            left.evaluate(time).unwrap(),
            reference(left_area, temporal, false),
            3e-12,
        );
        near(
            left.evaluate(time).unwrap() + right.evaluate(time).unwrap(),
            full.evaluate(time).unwrap(),
            3e-12,
        );
        near(
            proper.evaluate(time).unwrap(),
            reference(whole_area, temporal, true),
            3e-12,
        );
        real_sum_effect = real_sum_effect
            .max((full.evaluate(time).unwrap() - proper.evaluate(time).unwrap()).abs());
    }
    assert!(real_sum_effect > 1e-3);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let scale = match profile {
            HarmonicWindowProfile::Uniform => sinc(omega * window.duration_s()),
            HarmonicWindowProfile::Triangular => sinc(0.5 * omega * window.duration_s()).powi(2),
        };
        let temporal =
            rotation(2.0 * omega * (window.start_s() + 0.5 * window.duration_s())).scale(scale);
        near(
            full.observe_window(window).unwrap(),
            reference(whole_area, temporal, false),
            3e-12,
        );
    }
    // Nonzero real carriers require all spatial mode pairs even on a full cell.
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
                || false,
            )
            .is_err()
    );
}
