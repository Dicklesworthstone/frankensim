use super::*;

fn values(reading: AngularSpectrumExposureApertureProjectionReading) -> [f64; 3] {
    [
        reading.electric_mean_squared_v2_m2,
        reading.magnetic_mean_squared_a2_m2,
        reading.electric_magnetic_mean_product_w_m2,
    ]
}

fn fields() -> Vec<Vec<[C64; 6]>> {
    (0..2)
        .map(|f| {
            (0..2)
                .map(|m| {
                    core::array::from_fn(|c| {
                        C64::new(
                            (1 + f + c) as f64 / (7.0 * (m + 1) as f64),
                            (c as f64 - 2.0) * (m as f64 - 1.5) / (11.0 * (f + 1) as f64),
                        )
                    })
                })
                .collect()
        })
        .collect()
}

#[test]
fn chromatic_retarder_changes_circular_selection_at_each_original_frequency() {
    let wave = [
        C64::ONE,
        C64::new(0.0, 1.0),
        C64::ZERO,
        C64::new(0.0, -1.0),
        C64::ONE,
        C64::ZERO,
    ];
    let fields = parts(&[vec![wave], vec![wave]]);
    let r = core::f64::consts::FRAC_1_SQRT_2;
    // A 45-degree linear analyzer after an ideal differential delay:
    // quarter-wave at omega=1, half-wave at omega=2. a=J†d, not Jd.
    let quarter = [C64::from_re(r), C64::new(0.0, r), C64::ZERO];
    let half = [C64::from_re(r), C64::from_re(-r), C64::ZERO];
    let samples = [
        (
            1.0,
            AngularSpectrumExposureApertureJonesProjection::try_new(quarter, quarter).unwrap(),
        ),
        (
            2.0,
            AngularSpectrumExposureApertureJonesProjection::try_new(half, half).unwrap(),
        ),
    ];
    let spectrum =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
    let area = AngularSpectrumAperture::try_new([0.2, -0.7], [0.4, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.17, core::f64::consts::TAU).unwrap();
    for proper in [false, true] {
        let (actual, constant) = if proper {
            (
                observe::<ProperKernel, _>(
                    &[0.0],
                    &[0.0],
                    &[1.0, 2.0],
                    &fields,
                    area,
                    window,
                    Projected::spectral(spectrum),
                ),
                observe::<ProperKernel, _>(
                    &[0.0],
                    &[0.0],
                    &[1.0, 2.0],
                    &fields,
                    area,
                    window,
                    Projected::new(samples[0].1),
                ),
            )
        } else {
            (
                observe::<RealKernel, _>(
                    &[0.0],
                    &[0.0],
                    &[1.0, 2.0],
                    &fields,
                    area,
                    window,
                    Projected::spectral(spectrum),
                ),
                observe::<RealKernel, _>(
                    &[0.0],
                    &[0.0],
                    &[1.0, 2.0],
                    &fields,
                    area,
                    window,
                    Projected::new(samples[0].1),
                ),
            )
        };
        // |a1†E|²/2=1, |a2†E|²/2=1/2; the common period removes
        // inter-line beats. A wavelength-independent quarter-wave model gives 2.
        close(actual.electric_mean_squared_v2_m2, 1.5, 4e-15);
        close(actual.magnetic_mean_squared_a2_m2, 1.5, 4e-15);
        close(actual.electric_magnetic_mean_product_w_m2, 0.0, 4e-15);
        close(constant.electric_mean_squared_v2_m2, 2.0, 4e-15);
    }
}

#[test]
fn spectral_common_delay_moves_subcycle_shutters_and_beats_on_the_original_clock() {
    let omega = [1.3_f64, 2.1];
    let kx = [0.2, -0.7];
    let ky = [0.3];
    let fields = parts(&fields());
    let base = AngularSpectrumExposureApertureJonesProjection::try_new(
        [C64::new(0.3, -0.7), C64::new(-0.4, 0.2), C64::new(0.1, 0.3)],
        [C64::new(-0.2, 0.1), C64::new(0.6, -0.4), C64::new(0.5, 0.2)],
    )
    .unwrap();
    let delay = 0.43;
    let samples = omega.map(|frequency| {
        let (s, c) = (frequency * delay).sin_cos();
        let phase = C64::new(c, s);
        (
            frequency,
            AngularSpectrumExposureApertureJonesProjection::try_new(
                base.electric_weights().map(|z| z * phase),
                base.magnetic_weights().map(|z| z * phase),
            )
            .unwrap(),
        )
    });
    let spectrum =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.71, 0.19, profile).unwrap();
        let shifted =
            HarmonicWindow::with_profile(window.start_s() - delay, window.duration_s(), profile)
                .unwrap();
        for proper in [false, true] {
            let (actual, expected) = if proper {
                (
                    observe::<ProperKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::spectral(spectrum),
                    ),
                    observe::<ProperKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields,
                        area,
                        shifted,
                        Projected::new(base),
                    ),
                )
            } else {
                (
                    observe::<RealKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::spectral(spectrum),
                    ),
                    observe::<RealKernel, _>(
                        &kx,
                        &ky,
                        &omega,
                        &fields,
                        area,
                        shifted,
                        Projected::new(base),
                    ),
                )
            };
            for (actual, expected) in values(actual).into_iter().zip(values(expected)) {
                close(actual, expected, 8e-14 * expected.abs().max(1.0));
            }
        }
    }
}

#[test]
fn dispersive_elliptical_channels_match_independent_space_time_phase_quadrature() {
    let omega = [1.3, 2.1];
    let kx = [0.2, -0.7];
    let ky = [0.3];
    let fields = fields();
    let samples = [
        (
            omega[0],
            AngularSpectrumExposureApertureJonesProjection::try_new(
                [C64::new(0.3, -0.7), C64::new(-0.4, 0.2), C64::new(0.1, 0.3)],
                [C64::new(-0.2, 0.1), C64::new(0.6, -0.4), C64::new(0.5, 0.2)],
            )
            .unwrap(),
        ),
        (
            omega[1],
            AngularSpectrumExposureApertureJonesProjection::try_new(
                [C64::new(-0.1, 0.6), C64::new(0.2, -0.4), C64::new(0.3, 0.7)],
                [
                    C64::new(0.8, -0.3),
                    C64::new(-0.6, 0.5),
                    C64::new(0.4, -0.2),
                ],
            )
            .unwrap(),
        ),
    ];
    let spectrum =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
    let projected: Vec<Vec<[C64; 6]>> = fields
        .iter()
        .zip(samples)
        .map(|(line, (_, analyzer))| {
            line.iter()
                .map(|v| {
                    let mut result = [C64::ZERO; 6];
                    for c in 0..3 {
                        result[0] = result[0] + analyzer.electric_weights()[c].conj() * v[c];
                        result[3] = result[3] + analyzer.magnetic_weights()[c].conj() * v[c + 3];
                    }
                    result
                })
                .collect()
        })
        .collect();
    let fields_dd = parts(&fields);
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for proper in [false, true] {
            let actual = if proper {
                observe::<ProperKernel, _>(
                    &kx,
                    &ky,
                    &omega,
                    &fields_dd,
                    area,
                    window,
                    Projected::spectral(spectrum),
                )
            } else {
                observe::<RealKernel, _>(
                    &kx,
                    &ky,
                    &omega,
                    &fields_dd,
                    area,
                    window,
                    Projected::spectral(spectrum),
                )
            };
            let expected = super::super::super::tests::reference_matrix(
                &kx, &ky, &omega, &projected, area, window, proper,
            );
            for (actual, expected) in
                values(actual)
                    .into_iter()
                    .zip([expected[0][0], expected[3][3], expected[0][3]])
            {
                close(actual, expected, 3e-7 * expected.abs().max(1.0));
            }
        }
    }
}

#[test]
fn constant_spectral_response_matches_original_real_and_jones_paths_bit_for_bit() {
    let omega = [1.3, 2.1];
    let fields = parts(&fields());
    let area = AngularSpectrumAperture::try_new([-0.19, 0.23], [0.7, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.17, 0.43).unwrap();
    let real =
        AngularSpectrumExposureApertureProjection::try_new([0.6, -0.8, 0.0], [0.3, 0.4, -0.5])
            .unwrap();
    for base in [
        real.into(),
        AngularSpectrumExposureApertureJonesProjection::try_new(
            [C64::new(0.3, -0.7), C64::new(-0.4, 0.2), C64::new(0.1, 0.3)],
            [C64::new(-0.2, 0.1), C64::new(0.6, -0.4), C64::new(0.5, 0.2)],
        )
        .unwrap(),
    ] {
        let samples = omega.map(|frequency| (frequency, base));
        let spectrum =
            AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
        for proper in [false, true] {
            let (actual, expected) = if proper {
                (
                    observe::<ProperKernel, _>(
                        &[0.2, -0.7],
                        &[0.3],
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::spectral(spectrum),
                    ),
                    observe::<ProperKernel, _>(
                        &[0.2, -0.7],
                        &[0.3],
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::new(base),
                    ),
                )
            } else {
                (
                    observe::<RealKernel, _>(
                        &[0.2, -0.7],
                        &[0.3],
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::spectral(spectrum),
                    ),
                    observe::<RealKernel, _>(
                        &[0.2, -0.7],
                        &[0.3],
                        &omega,
                        &fields,
                        area,
                        window,
                        Projected::new(base),
                    ),
                )
            };
            assert_eq!(
                values(actual).map(f64::to_bits),
                values(expected).map(f64::to_bits)
            );
        }
    }
}

#[test]
fn spectral_support_is_exact_and_admission_is_cancellable() {
    let projection = AngularSpectrumExposureApertureJonesProjection::try_new(
        [C64::ONE, C64::ZERO, C64::ZERO],
        [C64::ZERO; 3],
    )
    .unwrap();
    for frequencies in [
        vec![],
        vec![0.0],
        vec![-1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![1.0, 1.0],
        vec![2.0, 1.0],
    ] {
        let samples: Vec<_> = frequencies.into_iter().map(|f| (f, projection)).collect();
        assert!(AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).is_err());
    }
    let samples = [(1.0, projection), (2.0, projection)];
    let spectrum =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
    assert_eq!(spectrum.samples(), &samples);
    assert_eq!(spectrum.projection_at(1.0).unwrap(), projection);
    for missing in [0.5, 1.0_f64.next_up(), 1.5, 2.0_f64.next_up(), f64::NAN] {
        assert!(spectrum.projection_at(missing).is_err());
    }
    let projected = Projected::spectral(spectrum);
    projected
        .admit_frequencies(&[1.0, 2.0], &mut || false)
        .unwrap();
    assert!(
        projected
            .admit_frequencies(&[1.0, 1.5], &mut || false)
            .is_err()
    );
    assert!(matches!(
        projected.admit_frequencies(&[1.0], &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let many: Vec<_> = (1..=65).map(|f| (f as f64, projection)).collect();
    let mut polls = 0;
    assert!(matches!(
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&many, || {
            polls += 1;
            polls == 3
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
}
