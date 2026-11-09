use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

const DEPTH: [AngularSpectrumHeatDepthSample; 1] = [AngularSpectrumHeatDepthSample {
    depth_m: 0.5,
    weight_m: 1.0,
}];

fn quad(area: AngularSpectrumHeatAreaQuadrature<'_>) -> AngularSpectrumHeatVolumeQuadrature<'_> {
    AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, 1.0],
        depth_samples: &DEPTH,
        transverse: area,
    }
}
fn grid(width: usize) -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width,
        height: 1,
        spacing_m: [1.0; 2],
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual / expected - 1.0).abs() < relative,
        "actual={actual:.17e}, expected={expected:.17e}, relative={relative:e}"
    );
}
fn point_mean(
    frequencies: &[f64],
    amplitudes: &[[Dd; 2]],
    proper: bool,
    duration: f64,
    epoch: Dd,
    conductivity: f64,
) -> ScaledDd {
    let quadrature = quad(AngularSpectrumHeatAreaQuadrature::PeriodicCell);
    let volume = admit_volume(grid(1), 1.0, quadrature, &mut || false).unwrap();
    positive_bank(
        frequencies,
        proper,
        HarmonicWindow::new(0.0, duration).unwrap(),
        epoch,
        quadrature,
        volume,
        None,
        1,
        |_| ScaledDd::from_f64(conductivity),
        |_, _| [ScaledDd::from_f64(1.0), ScaledDd::ZERO],
        |group, _, _| {
            let mut field = [[ScaledDd::ZERO; 2]; 3];
            field[0] = amplitudes[group].map(ScaledDd::from_dd);
            Ok(vec![field])
        },
        &mut || false,
    )
    .unwrap()
}

#[test]
fn proper_volume_heat_keeps_fourth_order_temporal_dark_leakage() {
    let carrier = 2.0_f64.powi(49);
    let omega = [carrier, carrier + 1.0, carrier + 2.0];
    let amplitudes = [1.0, -2.0, 1.0].map(|a| [Dd::from_f64(a), Dd::ZERO]);
    for duration in [1e-6, 1e-9, 1e-12] {
        let heat = point_mean(
            &omega,
            &amplitudes,
            true,
            duration,
            Dd::from_f64(-0.5 * duration),
            1.0,
        );
        close(
            heat.physical().unwrap().to_f64(),
            duration.powi(4) / 160.0,
            2e-12,
        );
    }
}

#[test]
fn real_odd_periodic_cell_keeps_both_spatial_quadratures_and_opposite_modes() {
    let frequencies = [1.0, 2.0, 3.0];
    let amplitudes = [1.0, -2.0, 1.0];
    let duration = 1e-9;
    let quadrature = quad(AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell);
    let volume = admit_volume(grid(3), 1.0, quadrature, &mut || false).unwrap();
    for imaginary in [false, true] {
        let value = positive_bank(
            &frequencies,
            false,
            HarmonicWindow::new(0.0, duration).unwrap(),
            Dd::from_f64(-duration / 2.0),
            quadrature,
            volume,
            None,
            1,
            |_| ScaledDd::from_f64(1.0),
            |_, _| [ScaledDd::from_f64(1.0), ScaledDd::ZERO],
            |group, _, _| {
                let mut result = vec![[[ScaledDd::ZERO; 2]; 3]; 3];
                let a = ScaledDd::from_f64(amplitudes[group]);
                if imaginary {
                    result[1][0][1] = a;
                    result[2][0][1] = -a;
                } else {
                    result[1][0][0] = a;
                    result[2][0][0] = a;
                }
                Ok(result)
            },
            &mut || false,
        )
        .unwrap()
        .physical()
        .unwrap()
        .to_f64();
        // 2 cos(kx) or 2 sin(kx), squared spatial mean=2; area=3.
        close(value, 3.0 * duration.powi(4) / 40.0, 5e-12);
    }
}

fn simpson_weight(index: usize, count: usize) -> f64 {
    if index == 0 || index == count {
        1.0
    } else if index % 2 == 0 {
        2.0
    } else {
        4.0
    }
}

#[test]
fn phased_off_center_real_rectangle_matches_independent_time_and_space_integral() {
    let omega = [1.0, 2.0, 3.0];
    let k = [0.0, -core::f64::consts::PI];
    let carrier = [1.7, 0.8];
    let aperture = AngularSpectrumAperture::try_new([0.37, -0.21], [0.19, 0.4]).unwrap();
    let quadrature = quad(AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(
        aperture,
    ));
    let volume = admit_volume(grid(2), 1.0, quadrature, &mut || false).unwrap();
    let rectangle = rectangle::PreparedRectangle::prepare(
        grid(2),
        aperture,
        carrier,
        |i| k[i],
        |_| 0.0,
        &mut || false,
    )
    .unwrap();
    let coefficients = [
        [C64::new(0.8, -0.4), C64::new(-0.7, 0.9)],
        [C64::new(-1.1, 0.2), C64::new(0.3, -0.6)],
        [C64::new(0.5, 0.7), C64::new(-0.2, 0.8)],
    ];
    // A complex, frequency-dependent reservoir transfer is retained before
    // projection, as are the original source and translated-aperture phases.
    let response = omega.map(|w| C64::new(1.0, -0.3 * w).scale(1.0 / (1.0 + 0.09 * w * w)));
    let start = 0.73;
    let duration = 0.04;
    let actual = positive_bank(
        &omega,
        false,
        HarmonicWindow::new(0.0, duration).unwrap(),
        Dd::from_f64(start),
        quadrature,
        volume,
        Some(&rectangle),
        1,
        |_| ScaledDd::from_f64(0.7),
        |group, _| parts(response[group]),
        |group, _, cancelled| {
            let raw = coefficients[group]
                .map(|a| {
                    let mut field = [[Dd::ZERO; 2]; 6];
                    field[0] = [Dd::from_f64(a.re), Dd::from_f64(a.im)];
                    field
                })
                .to_vec();
            rectangle.centered_electric(raw, cancelled)
        },
        &mut || false,
    )
    .unwrap()
    .physical()
    .unwrap()
    .to_f64();
    let (nt, nx, ny) = (80, 80, 40);
    let origin = aperture.origin_m();
    let extent = aperture.extent_m();
    let mut oracle = 0.0;
    for it in 0..=nt {
        let time = start + duration * it as f64 / nt as f64;
        for ix in 0..=nx {
            let x = origin[0] + extent[0] * ix as f64 / nx as f64;
            for iy in 0..=ny {
                let y = origin[1] + extent[1] * iy as f64 / ny as f64;
                let mut electric = 0.0;
                for group in 0..3 {
                    for mode in 0..2 {
                        let phase =
                            omega[group] * time - (k[mode] + carrier[0]) * x - carrier[1] * y;
                        let amplitude = coefficients[group][mode] * response[group];
                        electric += amplitude.re * phase.cos() - amplitude.im * phase.sin();
                    }
                }
                oracle += simpson_weight(it, nt)
                    * simpson_weight(ix, nx)
                    * simpson_weight(iy, ny)
                    * electric
                    * electric;
            }
        }
    }
    oracle *= 0.7 * aperture.area_m2() / (27.0 * (nt * nx * ny) as f64);
    close(actual, oracle, 3e-10);
}

#[test]
fn retained_dd_epoch_keeps_a_sub_ulp_accepted_step_in_the_source_phase() {
    let frequencies = [1.0, 2.0, 3.0];
    let epoch = 2.0_f64.powi(40);
    let offset = 1e-5;
    let duration = 1e-9;
    let weights = [1.0, -2.0, 1.0];
    let amplitudes = core::array::from_fn::<_, 3, _>(|group| {
        phase_rad_s_dd_from_clock_dd(
            Dd::from_f64(frequencies[group] - frequencies[0]),
            Dd::from_f64(epoch),
        )
        .unwrap()
        .map(|v| v * Dd::from_f64(weights[group]))
    });
    let heat = point_mean(
        &frequencies,
        &amplitudes,
        true,
        duration,
        Dd::from_pair(epoch, offset),
        1.0,
    );
    // Stable integral of (offset+t)^4/2; the next sin expansion is <2e-11.
    let expected = 0.5
        * (offset.powi(4)
            + 2.0 * offset.powi(3) * duration
            + 2.0 * offset.powi(2) * duration.powi(2)
            + offset * duration.powi(3)
            + duration.powi(4) / 5.0);
    close(heat.physical().unwrap().to_f64(), expected, 1e-8);
    let rounded = point_mean(
        &frequencies,
        &amplitudes,
        true,
        duration,
        Dd::from_f64(epoch),
        1.0,
    );
    assert!(rounded.physical().unwrap().to_f64() < expected * 1e-10);
}

#[test]
fn physical_scales_meet_before_square_and_interval_range_is_restored() {
    let frequencies = [1.0, 2.0, 3.0];
    let amplitudes = |value| {
        [
            [Dd::from_f64(value), Dd::ZERO],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
        ]
    };
    let large = point_mean(
        &frequencies,
        &amplitudes(1.5e154),
        true,
        1e-9,
        Dd::ZERO,
        1.0,
    );
    close(large.physical().unwrap().to_f64(), 1.125e308, 4e-15);
    let tiny = point_mean(
        &frequencies,
        &amplitudes(1e-200),
        true,
        1e-9,
        Dd::ZERO,
        1e200,
    );
    close(tiny.physical().unwrap().to_f64(), 0.5e-200, 4e-15);
    let watts = point_mean(
        &frequencies,
        &amplitudes(1.5e200),
        true,
        1e-100,
        Dd::ZERO,
        1.0,
    );
    assert!(watts.physical().is_err());
    close(
        (watts * ScaledDd::from_f64(1e-100))
            .physical()
            .unwrap()
            .to_f64(),
        1.125e300,
        4e-15,
    );
}

fn generous() -> AngularSpectrumHeatWindowLimits {
    AngularSpectrumHeatWindowLimits {
        maximum_depth_samples: usize::MAX,
        maximum_point_line_pairs: usize::MAX,
        maximum_fft_work: u64::MAX,
        maximum_temporal_projection_terms: usize::MAX,
        maximum_spatial_products: usize::MAX,
        maximum_fallback_terms: usize::MAX,
        maximum_workspace_bytes: usize::MAX,
    }
}

#[test]
fn explicit_work_admission_bounds_each_new_cost_before_field_reconstruction() {
    let mut work = Work::default();
    work.add_bank(3, Some(6), false, true, true, 4, 2, 2)
        .unwrap();
    work.add_bank(2, None, true, true, true, 4, 2, 2).unwrap();
    assert_eq!(work.temporal, 1728);
    assert_eq!(work.spatial, 2112);
    assert_eq!(work.fallback, 3);
    let exact = AngularSpectrumHeatWindowLimits {
        maximum_point_line_pairs: 14,
        maximum_fft_work: 0,
        maximum_temporal_projection_terms: work.temporal,
        maximum_spatial_products: work.spatial,
        maximum_fallback_terms: work.fallback,
        maximum_workspace_bytes: 300 + work.bytes,
        ..generous()
    };
    work.admit(exact, 7, 99, 100, 200, 2, true).unwrap();
    for which in 0..5 {
        let mut cap = exact;
        match which {
            0 => cap.maximum_point_line_pairs -= 1,
            1 => cap.maximum_temporal_projection_terms -= 1,
            2 => cap.maximum_spatial_products -= 1,
            3 => cap.maximum_fallback_terms -= 1,
            _ => cap.maximum_workspace_bytes -= 1,
        }
        assert!(work.admit(cap, 7, 99, 100, 200, 2, true).is_err());
    }
    assert!(work.admit(exact, 7, 1, 100, 200, 2, false).is_err());
    let mut impossible = Work::default();
    assert!(
        impossible
            .add_bank(usize::MAX, None, false, false, false, 4, 2, 2)
            .is_err()
    );
    let frequencies = [1.0, 2.0, 3.0];
    let quadrature = quad(AngularSpectrumHeatAreaQuadrature::PeriodicCell);
    let volume = admit_volume(grid(1), 1.0, quadrature, &mut || false).unwrap();
    let mut calls = 0;
    let result = positive_bank(
        &frequencies,
        false,
        HarmonicWindow::new(0.0, 1e-9).unwrap(),
        Dd::ZERO,
        quadrature,
        volume,
        None,
        1,
        |_| ScaledDd::from_f64(1.0),
        |_, _| [ScaledDd::from_f64(1.0), ScaledDd::ZERO],
        |_, _, _| {
            calls += 1;
            unreachable!()
        },
        &mut || true,
    );
    assert!(matches!(result, Err(AngularSpectrumError::Cancelled)));
    assert_eq!(calls, 0);
}

#[test]
fn sampled_partial_cells_keep_authored_depth_weights_and_reuse_original_planes() {
    let frequencies = [1.0, 2.0, 3.0];
    let depths = [
        AngularSpectrumHeatDepthSample {
            depth_m: 0.2,
            weight_m: 0.4,
        },
        AngularSpectrumHeatDepthSample {
            depth_m: 0.8,
            weight_m: 0.6,
        },
    ];
    let areas = [1.0, 0.25];
    let quadrature = AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, 1.0],
        depth_samples: &depths,
        transverse: AngularSpectrumHeatAreaQuadrature::CellAreasM2(&areas),
    };
    let volume = admit_volume(grid(2), 1.0, quadrature, &mut || false).unwrap();
    let mut calls = 0;
    let actual = positive_bank(
        &frequencies,
        true,
        HarmonicWindow::new(0.0, 0.02).unwrap(),
        Dd::from_f64(0.4),
        quadrature,
        volume,
        None,
        2,
        |reservoir| ScaledDd::from_f64(if reservoir == 0 { 3.0 } else { 1.0 }),
        |_, reservoir| {
            if reservoir == 0 {
                [ScaledDd::from_f64(1.0), ScaledDd::ZERO]
            } else {
                [ScaledDd::ZERO, ScaledDd::from_f64(2.0)]
            }
        },
        |group, depth, _| {
            calls += 1;
            Ok((0..2)
                .map(|point| {
                    let mut field = [[ScaledDd::ZERO; 2]; 3];
                    if group == 0 {
                        field[0][0] = ScaledDd::from_f64((point + 1) as f64 * (1.0 + depth));
                    }
                    field
                })
                .collect())
        },
        &mut || false,
    )
    .unwrap()
    .physical()
    .unwrap()
    .to_f64();
    assert_eq!(calls, frequencies.len() * depths.len());
    // Half population * (sigma + |reservoir transfer|²) * masked area * depth.
    close(
        actual,
        0.5 * 7.0 * (1.0 + 0.25 * 4.0) * (0.4 * 1.2 * 1.2 + 0.6 * 1.8 * 1.8),
        1e-14,
    );
}

#[test]
fn unsupported_profiles_keep_original_spectrum_and_dd_epoch_law() {
    let window =
        HarmonicWindow::with_profile(0.0, 0.07, HarmonicWindowProfile::Triangular).unwrap();
    assert!(
        temporal_modes(window, &[1.0, 2.0, 3.0], true)
            .unwrap()
            .is_none()
    );
    let epoch = Dd::from_pair(2.0_f64.powi(40), 1e-5);
    let terms = vec![
        AngularHarmonicSpectrumTerm {
            angular_frequency_rad_s: Dd::ZERO,
            peak_amplitude: [Dd::from_f64(2.0), Dd::ZERO],
        },
        AngularHarmonicSpectrumTerm {
            angular_frequency_rad_s: Dd::from_f64(0.2),
            peak_amplitude: [Dd::from_f64(0.3), Dd::from_f64(-0.1)],
        },
    ];
    let actual = fallback_mean(terms, window, epoch, &mut || false)
        .unwrap()
        .physical()
        .unwrap()
        .to_f64();
    let phase = phase_rad_s_dd_from_clock_dd(Dd::from_f64(0.2), epoch).unwrap();
    let amplitude = [
        Dd::from_f64(0.3) * phase[0] + Dd::from_f64(0.1) * phase[1],
        Dd::from_f64(0.3) * phase[1] - Dd::from_f64(0.1) * phase[0],
    ];
    let expected = Dd::from_f64(2.0)
        + window
            .product_kernel(0.2, 0.0)
            .unwrap()
            .mean_real_product_parts_dd(amplitude, [Dd::ONE, Dd::ZERO])
            .unwrap();
    close(actual, expected.to_f64(), 3e-15);
}
