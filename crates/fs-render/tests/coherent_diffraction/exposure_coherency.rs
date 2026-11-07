//! Full Stokes over the original chromatic pupil, source banks and shutter.
use super::*;
use fs_render::camera::PolarizedDiffractionCoherency;
use fs_render::polarization::Stokes;

fn polar(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}

fn scale(value: C64, factor: f64) -> C64 {
    C64::new(value.re * factor, value.im * factor)
}

fn matrix(value: PolarizedDiffractionCoherency) -> [f64; 4] {
    let [x, y] = value.axis_powers();
    let cross = value.cross();
    [x, y, cross.re, cross.im]
}

fn matrix_close(actual: PolarizedDiffractionCoherency, expected: [f64; 4], tolerance: f64) {
    let magnitude = (expected[0] + expected[1]).max(1e-280);
    for (a, b) in matrix(actual).into_iter().zip(expected) {
        assert!(
            (a - b).abs() <= tolerance * magnitude,
            "{actual:?} != {expected:?}"
        );
    }
    let s = actual.stokes();
    for (a, b) in [s.i, s.q, s.u, s.v].into_iter().zip([
        expected[0] + expected[1],
        expected[0] - expected[1],
        2.0 * expected[2],
        2.0 * expected[3],
    ]) {
        assert!((a - b).abs() <= 2.0 * tolerance * magnitude);
    }
    assert!(Stokes::try_new(s.i, s.q, s.u, s.v).is_some());
}

#[test]
fn g1_shutter_stokes_matches_complex_beats_and_original_analyzers() {
    let f0 = 299_792_458.0 / LAMBDA;
    let distance = 0.037;
    let at = [0.0302, -0.0001];
    let a = [CoherentJonesPupilWave {
        tangent_direction: [0.03, 0.0],
        amplitude_sqrt_w: [C64::new(0.8, -0.2), C64::ZERO],
    }];
    let b = [CoherentJonesPupilWave {
        tangent_direction: [0.0299, 0.00015],
        amplitude_sqrt_w: [C64::ZERO, C64::new(-0.3, 0.9)],
    }];
    for camera in [focused(), defocused(2.4)] {
        for delta in [1e9, 0.13 * f0] {
            let f1 = f0 + delta;
            let beat = f1 - f0;
            let field = |wave: CoherentJonesPupilWave, f: f64, axis: usize| {
                let offset = std::array::from_fn(|i| at[i] - wave.tangent_direction[i]);
                let wavelength = 299_792_458.0 / f * 1e9;
                let h = if camera.defocused_diffraction().is_some() {
                    camera.defocused_diffraction_amplitude(offset, wavelength)
                } else {
                    camera.in_focus_diffraction_amplitude(offset, wavelength)
                }
                .unwrap();
                h * wave.amplitude_sqrt_w[axis]
            };
            let x = field(a[0], f0, 0);
            let y = field(b[0], f1, 1);
            for (start, duration) in [(-0.13 / beat, 0.2 / beat), (0.0, 1.0 / beat)] {
                let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
                let phase = core::f64::consts::TAU
                    * beat
                    * (distance * (1.0 + 0.5 * (at[0] * at[0] + at[1] * at[1])) / 299_792_458.0
                        - start
                        - 0.5 * duration);
                let cross = scale(
                    x * y.conj() * polar(-phase),
                    sinc_reference(core::f64::consts::PI * beat * duration),
                );
                let expected = [x.norm_sq(), y.norm_sq(), cross.re, cross.im];
                for tones in [
                    [jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])],
                    [jones_tone(f1, &b, &[]), jones_tone(f0, &a, &[])],
                ] {
                    let query = camera
                        .prepare_coherent_jones_diffraction_exposure(
                            &tones,
                            distance,
                            shutter,
                            EXPOSURE_LIMITS,
                            || false,
                        )
                        .unwrap();
                    let measured = query.observe_coherency(at, || false).unwrap();
                    matrix_close(measured.coherency_density_w, expected, 3e-10);
                    assert_eq!(
                        measured.fluctuation_coherency_density_w.stokes(),
                        Stokes::unpolarized(0.0)
                    );
                    let s = measured.coherency_density_w.stokes();
                    close(
                        query
                            .observe(at, JONES_TOTAL, || false)
                            .unwrap()
                            .power_density_w,
                        s.i,
                        3e-13,
                    );
                    for (mode, expected) in [
                        ([C64::ONE, C64::ONE], 0.5 * (s.i + s.u)),
                        ([C64::ONE, C64::new(0.0, 1.0)], 0.5 * (s.i - s.v)),
                    ] {
                        let detector = CoherentJonesDetector::analyzer(mode).unwrap();
                        close(
                            query
                                .observe(at, detector, || false)
                                .unwrap()
                                .power_density_w,
                            expected,
                            3e-12,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g3_temporal_depolarization_retains_covariance_identity_and_frame_rotation() {
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = f0 + 1e9;
    let beat = f1 - f0;
    let a = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let b = [jones_wave([0.0; 2], [C64::ZERO, C64::new(0.0, f0 / f1)])];
    let angle = 0.37_f64;
    let rotate = |wave: CoherentJonesPupilWave| {
        let [x, y] = wave.amplitude_sqrt_w;
        CoherentJonesPupilWave {
            amplitude_sqrt_w: [
                scale(x, angle.cos()) + scale(y, angle.sin()),
                scale(y, angle.cos()) - scale(x, angle.sin()),
            ],
            ..wave
        }
    };
    let ar = a.map(rotate);
    let br = b.map(rotate);
    for width in [0.001, 0.3, 1.0] {
        let shutter = CoherentDiffractionShutter::try_new(-0.17 / beat, width / beat).unwrap();
        let query = |tones: &[CoherentJonesDiffractionTone<'_>]| {
            camera
                .prepare_coherent_jones_diffraction_exposure(
                    tones,
                    0.03,
                    shutter,
                    EXPOSURE_LIMITS,
                    || false,
                )
                .unwrap()
                .observe_coherency([0.0; 2], || false)
                .unwrap()
        };
        let expected = query(&[jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])]);
        let s = expected.coherency_density_w.stokes();
        assert!((s.degree() - sinc_reference(core::f64::consts::PI * width).abs()).abs() < 3e-12);
        let shared_a = [JonesPupilCoherenceFactor { waves: &a }];
        let shared_b = [JonesPupilCoherenceFactor { waves: &b }];
        let shared = query(&[
            jones_tone(f0, &[], &shared_a),
            jones_tone(f1, &[], &shared_b),
        ]);
        assert_eq!(
            shared.mean_coherency_density_w.stokes(),
            Stokes::unpolarized(0.0)
        );
        matrix_close(
            shared.fluctuation_coherency_density_w,
            matrix(expected.coherency_density_w),
            3e-13,
        );
        let independent_a = [
            JonesPupilCoherenceFactor { waves: &a },
            JonesPupilCoherenceFactor { waves: &[] },
        ];
        let independent_b = [
            JonesPupilCoherenceFactor { waves: &[] },
            JonesPupilCoherenceFactor { waves: &b },
        ];
        let independent = query(&[
            jones_tone(f0, &[], &independent_a),
            jones_tone(f1, &[], &independent_b),
        ]);
        assert_eq!(independent.coherency_density_w.cross(), C64::ZERO);
        close(independent.coherency_density_w.stokes().i, s.i, 3e-13);
        let rotated = query(&[jones_tone(f0, &ar, &[]), jones_tone(f1, &br, &[])])
            .coherency_density_w
            .stokes();
        let expected = s.rotate_axes(angle);
        for (actual, expected) in [rotated.i, rotated.q, rotated.u, rotated.v]
            .into_iter()
            .zip([expected.i, expected.q, expected.u, expected.v])
        {
            assert!((actual - expected).abs() < 3e-12 * s.i);
        }
    }
}

#[test]
fn g1_finite_pixel_stokes_matches_separable_chromatic_pupil_and_energy() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = 1.07 * f0;
    let beat = f1 - f0;
    let distance = 0.037;
    let centre = [0.03, 0.01];
    let offsets = [
        [centre[0] - 2.2 * SCALE, centre[1] + 0.4 * SCALE],
        [centre[0] + 1.4 * SCALE, centre[1] - 0.7 * SCALE],
    ];
    let amplitude = C64::new(0.3, 0.7);
    let a = [CoherentJonesPupilWave {
        tangent_direction: offsets[0],
        amplitude_sqrt_w: [C64::ONE, C64::ZERO],
    }];
    let b = [CoherentJonesPupilWave {
        tangent_direction: offsets[1],
        amplitude_sqrt_w: [C64::ZERO, amplitude],
    }];
    let tones = [jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])];
    let shutter = CoherentDiffractionShutter::try_new(-0.13 / beat, 0.37 / beat).unwrap();
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            distance,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let footprint = CoherentDiffractionPixel::try_new(
        [centre[0] - 5.1 * SCALE, centre[1] - 3.3 * SCALE],
        [centre[0] + 4.3 * SCALE, centre[1] + 2.8 * SCALE],
    )
    .unwrap();
    let measured = query
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let half_side = RADIUS * core::f64::consts::FRAC_1_SQRT_2;
    let wavelengths = [299_792_458.0 / f0, 299_792_458.0 / f1];
    let delta_k = core::f64::consts::TAU * beat / 299_792_458.0;
    let shape = |tone: usize, axis: usize, q: f64| {
        sinc_reference(
            core::f64::consts::TAU * half_side * (q - offsets[tone][axis]) / wavelengths[tone],
        )
    };
    let low = footprint.minimum_tangent();
    let high = footprint.maximum_tangent();
    let diagonal = std::array::from_fn::<_, 2, _>(|tone| {
        4.0 * half_side * half_side / wavelengths[tone].powi(2)
            * (if tone == 0 { 1.0 } else { amplitude.norm_sq() })
            * (0..2)
                .map(|axis| {
                    integrate_reference(low[axis], high[axis], |q| shape(tone, axis, q).powi(2))
                })
                .product::<f64>()
    });
    let mut cross = amplitude.conj();
    for axis in 0..2 {
        let integral = |imaginary: bool| {
            integrate_reference(low[axis], high[axis], |q| {
                let phase = -0.5 * delta_k * distance * (q * q - centre[axis] * centre[axis]);
                shape(0, axis, q)
                    * shape(1, axis, q)
                    * (if imaginary { phase.sin() } else { phase.cos() })
            })
        };
        cross = cross * C64::new(integral(false), integral(true));
    }
    let phase = -delta_k * distance * (1.0 + 0.5 * (centre[0] * centre[0] + centre[1] * centre[1]))
        + core::f64::consts::TAU * beat * (shutter.start_s() + 0.5 * shutter.duration_s());
    cross = scale(
        cross * polar(phase),
        4.0 * half_side * half_side / (wavelengths[0] * wavelengths[1])
            * sinc_reference(core::f64::consts::PI * beat * shutter.duration_s()),
    );
    let expected = [diagonal[0], diagonal[1], cross.re, cross.im];
    matrix_close(measured.average_power.coherency_w, expected, 3e-8);
    matrix_close(
        measured.coherency_j,
        expected.map(|x| x * shutter.duration_s()),
        3e-8,
    );
    let old = query
        .integrate_pixel(footprint, JONES_TOTAL, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    close(
        measured.average_power.coherency_w.stokes().i,
        old.average_power.power_w,
        3e-13,
    );
    assert_eq!(
        measured.average_power.pupil_evaluations,
        old.average_power.pupil_evaluations
    );
    assert_eq!(
        measured.average_power.pupil_evaluations,
        2 * measured.average_power.quadrature_samples
    );
    assert_eq!(
        measured.pair_contractions,
        14 * measured.average_power.quadrature_samples
    );
    let refined = query
        .integrate_pixel_coherency(
            footprint,
            CoherentDiffractionExposurePixelLimits {
                spatial: CoherentDiffractionPixelLimits {
                    maximum_phase_span_rad: 0.5,
                    ..PIXEL_LIMITS
                },
                ..EXPOSURE_PIXEL_LIMITS
            },
            || false,
        )
        .unwrap();
    matrix_close(
        refined.average_power.coherency_w,
        matrix(measured.average_power.coherency_w),
        2e-11,
    );
}

#[test]
fn g1_monochromatic_and_weak_coherencies_survive_the_original_shutter() {
    let camera = focused();
    let frequency = 299_792_458.0 / LAMBDA;
    let weak = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1e-120)])];
    let tones = [jones_tone(frequency, &weak, &[])];
    let shutter = CoherentDiffractionShutter::try_new(1e9, 1e-10).unwrap();
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let monochromatic = camera
        .prepare_polarized_coherent_diffraction(
            LAMBDA_NM,
            CoherentJonesPupilIllumination {
                mean: &weak,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let measured = query.observe_coherency([0.0; 2], || false).unwrap();
    matrix_close(
        measured.coherency_density_w,
        matrix(
            monochromatic
                .observe([0.0; 2], || false)
                .unwrap()
                .coherency_density_w,
        ),
        3e-12,
    );
    let footprint = pixel([-0.4; 2], [0.6; 2]);
    let measured_pixel = query
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    matrix_close(
        measured_pixel.average_power.coherency_w,
        matrix(
            monochromatic
                .integrate_pixel(footprint, PIXEL_LIMITS, || false)
                .unwrap()
                .coherency_w,
        ),
        3e-12,
    );
    let second = [jones_wave(
        [0.0; 2],
        [C64::new(0.2, 0.3), C64::new(-0.3e-120, 0.2e-120)],
    )];
    let spectral = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[
                jones_tone(frequency, &weak, &[]),
                jones_tone(frequency + 1e9, &second, &[]),
            ],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let spectral_coherency = spectral
        .observe_coherency([0.0; 2], || false)
        .unwrap()
        .coherency_density_w;
    let vertical = CoherentJonesDetector::analyzer([C64::ZERO, C64::ONE]).unwrap();
    let vertical_power = spectral
        .observe([0.0; 2], vertical, || false)
        .unwrap()
        .power_density_w;
    assert!((spectral_coherency.axis_powers()[1] / vertical_power - 1.0).abs() < 2e-12);
    for result in [
        measured.coherency_density_w,
        measured_pixel.average_power.coherency_w,
        measured_pixel.coherency_j,
        spectral_coherency,
    ] {
        let axes = result.axis_powers();
        assert!(axes[1] > 0.0);
        assert!((axes[1] / (axes[0] * 1e-240) - 1.0).abs() < 2e-12);
        assert_eq!(result.stokes().i, result.stokes().q);
        assert!(result.stokes().v < 0.0);
    }
    let strong = [jones_wave(
        [0.0; 2],
        [C64::new(1e16, -1e16), C64::new(-1e16, -1e16)],
    )];
    let opposite = [CoherentJonesPupilWave {
        amplitude_sqrt_w: strong[0].amplitude_sqrt_w.map(|a| C64::new(-a.re, -a.im)),
        ..strong[0]
    }];
    let small = [jones_wave(
        [0.0; 2],
        [C64::new(0.3, -0.2), C64::new(-0.4, 0.1)],
    )];
    let observe = |tones: &[CoherentJonesDiffractionTone<'_>]| {
        camera
            .prepare_coherent_jones_diffraction_exposure(
                tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap()
            .observe_coherency([0.0; 2], || false)
            .unwrap()
            .coherency_density_w
    };
    let expected = observe(&[jones_tone(frequency, &small, &[])]);
    let actual = observe(&[
        jones_tone(frequency, &strong, &[]),
        jones_tone(frequency, &small, &[]),
        jones_tone(frequency, &opposite, &[]),
    ]);
    matrix_close(actual, matrix(expected), 3e-12);
    assert_eq!(
        observe(&[
            jones_tone(frequency, &strong, &[]),
            jones_tone(frequency, &opposite, &[])
        ])
        .stokes(),
        Stokes::unpolarized(0.0)
    );
}

#[test]
fn g1_full_stokes_converges_to_independent_time_quadrature() {
    let camera = defocused(1.7);
    let f0 = 299_792_458.0 / LAMBDA;
    let frequencies = [f0, f0 + 1.3e9, f0 - 0.7e9];
    let waves = [
        [jones_wave(
            [-0.3, 0.2],
            [C64::new(0.8, -0.2), C64::new(-0.1, 0.4)],
        )],
        [jones_wave(
            [0.6, -0.1],
            [C64::new(-0.2, 0.7), C64::new(0.6, 0.1)],
        )],
        [jones_wave(
            [-0.1, -0.3],
            [C64::new(0.3, 0.4), C64::new(-0.5, 0.2)],
        )],
    ];
    let tones = std::array::from_fn::<_, 3, _>(|i| jones_tone(frequencies[i], &waves[i], &[]));
    let shutter = CoherentDiffractionShutter::try_new(-0.13e-9, 0.41e-9).unwrap();
    let distance = 0.037;
    let at = [0.8 * SCALE, -0.6 * SCALE];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            distance,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let expected = matrix(
        query
            .observe_coherency(at, || false)
            .unwrap()
            .coherency_density_w,
    );
    let fields = std::array::from_fn::<_, 3, _>(|i| {
        let h = camera
            .defocused_diffraction_amplitude(
                std::array::from_fn(|axis| at[axis] - waves[i][0].tangent_direction[axis]),
                299_792_458.0 / frequencies[i] * 1e9,
            )
            .unwrap();
        waves[i][0].amplitude_sqrt_w.map(|a| h * a)
    });
    let delay = distance * (1.0 + 0.5 * (at[0] * at[0] + at[1] * at[1])) / 299_792_458.0;
    let mut previous = f64::INFINITY;
    for samples in [32, 64, 128] {
        let mut integral = [0.0; 4];
        for k in 0..samples {
            let time = shutter.start_s() + (k as f64 + 0.5) * shutter.duration_s() / samples as f64;
            let mut field = [C64::ZERO; 2];
            for i in 0..3 {
                let rotation =
                    polar(core::f64::consts::TAU * (frequencies[i] - f0) * (delay - time));
                for axis in 0..2 {
                    field[axis] = field[axis] + fields[i][axis] * rotation;
                }
            }
            let cross = field[0] * field[1].conj();
            for (sum, value) in integral.iter_mut().zip([
                field[0].norm_sq(),
                field[1].norm_sq(),
                cross.re,
                cross.im,
            ]) {
                *sum += value / samples as f64;
            }
        }
        let error = integral
            .into_iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
            / (expected[0] + expected[1]);
        assert!(error < previous * 0.45);
        previous = error;
    }
    assert!(previous < 1e-4);
}

#[test]
fn g0_g4_coherency_preflights_original_budgets_and_reuses_after_cancellation() {
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let waves = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 0.3)])];
    let tones = [
        jones_tone(f0, &waves, &[]),
        jones_tone(f0 + 1e9, &waves, &[]),
    ];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.2e-9).unwrap();
    let original_budget = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: 6,
                ..EXPOSURE_LIMITS
            },
            || false,
        )
        .unwrap();
    assert!(
        original_budget
            .observe([0.0; 2], JONES_TOTAL, || false)
            .is_ok()
    );
    let mut polls = 0;
    assert!(matches!(
        original_budget.observe_coherency([0.0; 2], || {
            polls += 1;
            false
        }),
        Err(CoherentDiffractionError::Limit {
            quantity: "Jones coherency pair contractions",
            requested: 14,
            maximum: 6
        })
    ));
    assert_eq!(polls, 1);
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: 14,
                ..EXPOSURE_LIMITS
            },
            || false,
        )
        .unwrap();
    let point = query.observe_coherency([0.0; 2], || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        query.observe_coherency([0.0; 2], || {
            polls += 1;
            polls == 6
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(query.observe_coherency([0.0; 2], || false).unwrap(), point);
    let footprint = pixel([-0.4; 2], [0.6; 2]);
    let baseline = query
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    assert_eq!(
        baseline.pair_contractions,
        14 * baseline.average_power.quadrature_samples
    );
    let mut polls = 0;
    assert!(matches!(
        query.integrate_pixel_coherency(
            footprint,
            CoherentDiffractionExposurePixelLimits {
                maximum_pair_contractions: baseline.pair_contractions - 1,
                ..EXPOSURE_PIXEL_LIMITS
            },
            || {
                polls += 1;
                false
            }
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel exposure pair contractions",
            ..
        })
    ));
    assert_eq!(polls, 1);
    let mut polls = 0;
    assert_eq!(
        query.integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || {
            polls += 1;
            polls == 31
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query
            .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
            .unwrap(),
        baseline
    );
    for a in [1e200, 1e-200] {
        let range = [jones_wave([0.0; 2], [C64::new(a, 0.0), C64::ZERO])];
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &[jones_tone(f0, &range, &[])],
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        assert_eq!(
            query.observe_coherency([0.0; 2], || false),
            Err(CoherentDiffractionError::Unrepresentable)
        );
    }
    let tiny = [jones_wave(
        [0.0; 2],
        [C64::new(f64::from_bits(1), 0.0), C64::ZERO],
    )];
    let underflow = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[jones_tone(f0, &tiny, &[])],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let at = [0.075, 0.0];
    let h = camera
        .in_focus_diffraction_amplitude(at, LAMBDA_NM)
        .unwrap();
    assert!(h.re != 0.0 || h.im != 0.0);
    assert!(h.re.abs() < 0.5 && h.im.abs() < 0.5);
    assert_eq!(
        underflow.observe_coherency(at, || false),
        Err(CoherentDiffractionError::Unrepresentable)
    );
    assert_eq!(
        underflow.observe(at, JONES_TOTAL, || false),
        Err(CoherentDiffractionError::Unrepresentable)
    );
    let empty = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[],
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_tones: 0,
                maximum_pair_contractions: 0,
                sources: CoherentDiffractionLimits {
                    maximum_waves: 0,
                    maximum_factors: 0,
                },
            },
            || false,
        )
        .unwrap();
    let dark = empty
        .integrate_pixel_coherency(
            footprint,
            CoherentDiffractionExposurePixelLimits {
                maximum_pair_contractions: 0,
                spatial: CoherentDiffractionPixelLimits {
                    maximum_samples: 0,
                    maximum_pupil_evaluations: 0,
                    ..PIXEL_LIMITS
                },
            },
            || false,
        )
        .unwrap();
    assert_eq!(dark.coherency_j.stokes(), Stokes::unpolarized(0.0));
    assert_eq!(dark.pair_contractions, 0);
}
