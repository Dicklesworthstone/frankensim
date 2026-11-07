#![cfg(feature = "chart-backends")]
//! G1 pupil-integral references and G3 coherence/phase invariants.

use fs_geom::{Point3, Vec3};
use fs_math::c64::C64;
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, CoherentDiffractionError,
    CoherentDiffractionExposureLimits, CoherentDiffractionExposurePixelLimits,
    CoherentDiffractionLimits, CoherentDiffractionPixel, CoherentDiffractionPixelLimits,
    CoherentDiffractionPixelPower, CoherentDiffractionShutter, CoherentDiffractionTone,
    CoherentJonesDetector, CoherentJonesDiffractionTone, CoherentJonesPupilIllumination,
    CoherentJonesPupilWave, CoherentPupilIllumination, CoherentPupilWave, DefocusedDiffraction,
    FocusedDiffraction, JonesPupilCoherenceFactor, PhysicalCamera, PupilCoherenceFactor,
};

const JONES_TOTAL: CoherentJonesDetector = CoherentJonesDetector::POLARIZATION_INSENSITIVE;

fn jones_wave(q: [f64; 2], amplitude_sqrt_w: [C64; 2]) -> CoherentJonesPupilWave {
    CoherentJonesPupilWave {
        tangent_direction: q.map(|v| v * SCALE),
        amplitude_sqrt_w,
    }
}

fn jones_tone<'a>(
    frequency_hz: f64,
    mean: &'a [CoherentJonesPupilWave],
    factors: &'a [JonesPupilCoherenceFactor<'a>],
) -> CoherentJonesDiffractionTone<'a> {
    CoherentJonesDiffractionTone {
        frequency_hz,
        illumination: CoherentJonesPupilIllumination { mean, factors },
    }
}

#[test]
fn g1_jones_detector_orthogonality_malus_and_complex_analyzers() {
    let camera = focused();
    let frequency = 299_792_458.0 / LAMBDA;
    let shutter = CoherentDiffractionShutter::try_new(3.0, 0.02).unwrap();
    let one = core::f64::consts::PI * RADIUS * RADIUS / (LAMBDA * LAMBDA);
    let orthogonal = [
        jones_wave([0.0; 2], [C64::ONE, C64::ZERO]),
        jones_wave([0.0; 2], [C64::ZERO, C64::new(-1.0, 0.0)]),
    ];
    let tones = [jones_tone(frequency, &orthogonal, &[])];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let plus = CoherentJonesDetector::analyzer([C64::ONE, C64::ONE]).unwrap();
    let minus = CoherentJonesDetector::analyzer([C64::ONE, C64::new(-1.0, 0.0)]).unwrap();
    close(
        query
            .observe([0.0; 2], JONES_TOTAL, || false)
            .unwrap()
            .power_density_w,
        2.0 * one,
        3e-12,
    );
    assert_eq!(
        query
            .observe([0.0; 2], plus, || false)
            .unwrap()
            .power_density_w,
        0.0
    );
    close(
        query
            .observe([0.0; 2], minus, || false)
            .unwrap()
            .power_density_w,
        2.0 * one,
        3e-12,
    );

    let horizontal = [orthogonal[0]];
    let tones = [jones_tone(frequency, &horizontal, &[])];
    let horizontal = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    for angle in [0.0_f64, 0.23, 0.71, core::f64::consts::FRAC_PI_2] {
        let analyzer = CoherentJonesDetector::analyzer([
            C64::new(angle.cos(), 0.0),
            C64::new(angle.sin(), 0.0),
        ])
        .unwrap();
        close(
            horizontal
                .observe([0.0; 2], analyzer, || false)
                .unwrap()
                .power_density_w,
            one * angle.cos().powi(2),
            3e-12,
        );
    }
    let circular = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1.0)])];
    let tones = [jones_tone(frequency, &circular, &[])];
    let circular = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let matched = CoherentJonesDetector::analyzer([C64::ONE, C64::new(0.0, 1.0)]).unwrap();
    let crossed = CoherentJonesDetector::analyzer([C64::ONE, C64::new(0.0, -1.0)]).unwrap();
    close(
        circular
            .observe([0.0; 2], matched, || false)
            .unwrap()
            .power_density_w,
        2.0 * one,
        3e-12,
    );
    assert_eq!(
        circular
            .observe([0.0; 2], crossed, || false)
            .unwrap()
            .power_density_w,
        0.0
    );
    // Analyzer magnitude/global phase is not gain or an extra source phase.
    let same_mode =
        CoherentJonesDetector::analyzer([C64::new(0.0, 1e308), C64::new(-1e308, 0.0)]).unwrap();
    assert_eq!(
        circular.observe([0.0; 2], same_mode, || false).unwrap(),
        circular.observe([0.0; 2], matched, || false).unwrap()
    );
}

#[test]
fn g1_jones_analyzer_reveals_only_authored_cross_frequency_coherence() {
    let camera = defocused(2.4);
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
    for delta in [1e9, 0.13 * f0] {
        let f1 = f0 + delta;
        let beat = f1 - f0;
        let field = |wave: CoherentJonesPupilWave, frequency: f64, component: usize| {
            camera
                .defocused_diffraction_amplitude(
                    std::array::from_fn(|axis| at[axis] - wave.tangent_direction[axis]),
                    299_792_458.0 / frequency * 1e9,
                )
                .unwrap()
                * wave.amplitude_sqrt_w[component]
        };
        let x = field(a[0], f0, 0);
        let y = field(b[0], f1, 1);
        let diagonal = x.norm_sq() + y.norm_sq();
        for (start, duration) in [(0.0, 0.2 / beat), (-0.13 / beat, 1.0 / beat)] {
            let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
            let tones = [jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])];
            let query = camera
                .prepare_coherent_jones_diffraction_exposure(
                    &tones,
                    distance,
                    shutter,
                    EXPOSURE_LIMITS,
                    || false,
                )
                .unwrap();
            close(
                query
                    .observe(at, JONES_TOTAL, || false)
                    .unwrap()
                    .power_density_w,
                diagonal,
                3e-12,
            );
            let phase = core::f64::consts::TAU
                * beat
                * (distance * (1.0 + 0.5 * (at[0] * at[0] + at[1] * at[1])) / 299_792_458.0
                    - start
                    - 0.5 * duration);
            for analyzer_phase in [0.0_f64, 0.67] {
                let mode_y = C64::new(analyzer_phase.cos(), analyzer_phase.sin());
                let analyzer = CoherentJonesDetector::analyzer([C64::ONE, mode_y]).unwrap();
                let cross = y * x.conj() * mode_y.conj() * C64::new(phase.cos(), phase.sin());
                let expected = 0.5 * diagonal
                    + cross.re * sinc_reference(core::f64::consts::PI * beat * duration);
                close(
                    query
                        .observe(at, analyzer, || false)
                        .unwrap()
                        .power_density_w,
                    expected,
                    3e-10,
                );
                let af = [JonesPupilCoherenceFactor { waves: &a }];
                let bf = [JonesPupilCoherenceFactor { waves: &b }];
                let shared = [jones_tone(f0, &[], &af), jones_tone(f1, &[], &bf)];
                let shared = camera
                    .prepare_coherent_jones_diffraction_exposure(
                        &shared,
                        distance,
                        shutter,
                        EXPOSURE_LIMITS,
                        || false,
                    )
                    .unwrap()
                    .observe(at, analyzer, || false)
                    .unwrap();
                assert_eq!(shared.mean_power_density_w, 0.0);
                close(shared.fluctuation_power_density_w, expected, 3e-10);
                let af = [
                    JonesPupilCoherenceFactor { waves: &a },
                    JonesPupilCoherenceFactor { waves: &[] },
                ];
                let bf = [
                    JonesPupilCoherenceFactor { waves: &[] },
                    JonesPupilCoherenceFactor { waves: &b },
                ];
                let independent = [jones_tone(f0, &[], &af), jones_tone(f1, &[], &bf)];
                let independent = camera
                    .prepare_coherent_jones_diffraction_exposure(
                        &independent,
                        distance,
                        shutter,
                        EXPOSURE_LIMITS,
                        || false,
                    )
                    .unwrap()
                    .observe(at, analyzer, || false)
                    .unwrap();
                close(independent.power_density_w, 0.5 * diagonal, 3e-12);
            }
        }
    }
}

#[test]
fn g1_jones_pixel_matches_square_pupil_analyzer_integrals_and_complementary_power() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let waves = [
        jones_wave([-2.2, 0.0], [C64::ONE, C64::ZERO]),
        jones_wave([2.2, 0.0], [C64::ZERO, C64::new(-1.0, 0.0)]),
    ];
    let tones = [jones_tone(299_792_458.0 / LAMBDA, &waves, &[])];
    let shutter = CoherentDiffractionShutter::try_new(-0.07, 0.04).unwrap();
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let footprint = pixel([-7.1, -4.3], [8.3, 2.8]);
    let plus = CoherentJonesDetector::analyzer([C64::ONE, C64::ONE]).unwrap();
    let minus = CoherentJonesDetector::analyzer([C64::ONE, C64::new(-1.0, 0.0)]).unwrap();
    let total = query
        .integrate_pixel(footprint, JONES_TOTAL, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let selected = query
        .integrate_pixel(footprint, plus, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let complement = query
        .integrate_pixel(footprint, minus, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    // Independent analytic square-pupil transform and separable Simpson rule.
    let sinc = |q: f64| sinc_reference(q * core::f64::consts::FRAC_1_SQRT_2);
    let y = integrate_reference(-4.3, 2.8, |q| sinc(q).powi(2));
    let norm = y / (2.0 * core::f64::consts::PI.powi(2));
    let expected_total =
        norm * integrate_reference(-7.1, 8.3, |q| sinc(q + 2.2).powi(2) + sinc(q - 2.2).powi(2));
    let expected_selected =
        norm * integrate_reference(-7.1, 8.3, |q| 0.5 * (sinc(q + 2.2) - sinc(q - 2.2)).powi(2));
    close(total.average_power.power_w, expected_total, 3e-10);
    close(selected.average_power.power_w, expected_selected, 3e-10);
    close(
        selected.energy_j,
        shutter.duration_s() * expected_selected,
        3e-10,
    );
    close(
        selected.average_power.power_w + complement.average_power.power_w,
        total.average_power.power_w,
        3e-13,
    );
    assert_eq!(
        total.average_power.pupil_evaluations,
        2 * total.average_power.quadrature_samples
    );
    assert_eq!(
        selected.average_power.pupil_evaluations,
        total.average_power.pupil_evaluations
    );
    assert_eq!(total.pair_contractions, 2 * selected.pair_contractions);
    let refined = query
        .integrate_pixel(
            footprint,
            plus,
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
    close(refined.energy_j, selected.energy_j, 2e-11);
}

#[test]
fn g3_jones_projection_retains_weak_fields_and_rotates_with_its_basis() {
    let camera = defocused(2.4);
    let frequency = 299_792_458.0 / LAMBDA;
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.1).unwrap();
    let a = [jones_wave([0.2, -0.1], [C64::new(1e16, -2e16), C64::ZERO])];
    let small = [jones_wave([0.2, -0.1], [C64::new(0.5, 0.75), C64::ZERO])];
    let b = [jones_wave([0.2, -0.1], [C64::ZERO, C64::new(-1e16, 2e16)])];
    let tones = [
        jones_tone(frequency, &a, &[]),
        jones_tone(frequency, &small, &[]),
        jones_tone(frequency, &b, &[]),
    ];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    assert_eq!(query.frequency_count(), 1);
    let analyzer = CoherentJonesDetector::analyzer([C64::ONE, C64::ONE]).unwrap();
    let at = [0.5 * SCALE, -0.3 * SCALE];
    let h = camera
        .defocused_diffraction_amplitude([0.3 * SCALE, -0.2 * SCALE], LAMBDA_NM)
        .unwrap();
    let measured = query
        .observe(at, analyzer, || false)
        .unwrap()
        .power_density_w;
    let expected = 0.5 * h.norm_sq() * C64::new(0.5, 0.75).norm_sq();
    assert!(measured > 0.0);
    close(measured, expected, 2e-12);

    let angle = 0.37_f64;
    let cosine = C64::new(angle.cos(), 0.0);
    let sine = C64::new(angle.sin(), 0.0);
    let rotate = |z: [C64; 2]| [z[0] * cosine + z[1] * sine, z[1] * cosine - z[0] * sine];
    let source = [jones_wave(
        [0.2, -0.1],
        [C64::new(0.8, 0.2), C64::new(-0.3, 0.5)],
    )];
    let rotated = [CoherentJonesPupilWave {
        amplitude_sqrt_w: rotate(source[0].amplitude_sqrt_w),
        ..source[0]
    }];
    let mode = [C64::new(0.4, -0.2), C64::new(0.3, 0.7)];
    let observe = |waves: &[CoherentJonesPupilWave], mode| {
        camera
            .prepare_coherent_jones_diffraction_exposure(
                &[jones_tone(frequency, waves, &[])],
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap()
            .observe(at, CoherentJonesDetector::analyzer(mode).unwrap(), || false)
            .unwrap()
            .power_density_w
    };
    close(
        observe(&rotated, rotate(mode)),
        observe(&source, mode),
        3e-13,
    );
}

#[test]
fn g0_g4_jones_admission_exact_work_and_cancelled_reuse() {
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.2e-9).unwrap();
    for mode in [
        [C64::ZERO; 2],
        [C64::new(f64::NAN, 0.0), C64::ONE],
        [C64::ONE, C64::new(0.0, f64::INFINITY)],
    ] {
        assert_eq!(
            CoherentJonesDetector::analyzer(mode),
            Err(CoherentDiffractionError::InvalidPolarizationAnalyzer)
        );
    }
    let waves = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1.0)])];
    let tones = [
        jones_tone(f0, &waves, &[]),
        jones_tone(f0 + 1e9, &waves, &[]),
    ];
    assert!(matches!(
        camera.prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: 5,
                ..EXPOSURE_LIMITS
            },
            || false,
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "Jones exposure pair contractions",
            requested: 6,
            maximum: 5
        })
    ));
    let equal = [tones[0], tones[0]];
    assert!(
        camera
            .prepare_coherent_jones_diffraction_exposure(
                &equal,
                0.03,
                shutter,
                CoherentDiffractionExposureLimits {
                    maximum_pair_contractions: 2,
                    ..EXPOSURE_LIMITS
                },
                || false,
            )
            .is_ok()
    );
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let baseline = query.observe([0.0; 2], JONES_TOTAL, || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        query
            .integrate_pixel(
                pixel([-0.1; 2], [0.1; 2]),
                JONES_TOTAL,
                CoherentDiffractionExposurePixelLimits {
                    maximum_pair_contractions: 1,
                    ..EXPOSURE_PIXEL_LIMITS
                },
                || {
                    polls += 1;
                    false
                }
            )
            .unwrap_err(),
        CoherentDiffractionError::Limit {
            quantity: "pixel exposure pair contractions",
            requested: 96,
            maximum: 1,
        }
    );
    assert_eq!(
        polls, 1,
        "pair work must be refused before pupil evaluation"
    );
    let mut polls = 0;
    assert_eq!(
        query.integrate_pixel(
            pixel([-1.0; 2], [1.0; 2]),
            JONES_TOTAL,
            EXPOSURE_PIXEL_LIMITS,
            || {
                polls += 1;
                polls == 11
            }
        ),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
        baseline
    );
    let factors = [JonesPupilCoherenceFactor { waves: &waves }];
    let inconsistent = [tones[0], jones_tone(f0 + 1e9, &[], &factors)];
    assert!(matches!(
        camera.prepare_coherent_jones_diffraction_exposure(
            &inconsistent,
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        ),
        Err(CoherentDiffractionError::InconsistentCoherenceFactors)
    ));
    let invalid = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, f64::NAN)])];
    assert!(matches!(
        camera.prepare_coherent_jones_diffraction_exposure(
            &[jones_tone(f0, &invalid, &[])],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        ),
        Err(CoherentDiffractionError::NonFiniteAmplitude)
    ));
    let empty = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[],
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_tones: 0,
                sources: CoherentDiffractionLimits {
                    maximum_waves: 0,
                    maximum_factors: 0,
                },
                maximum_pair_contractions: 0,
            },
            || false,
        )
        .unwrap()
        .integrate_pixel(
            pixel([-1.0; 2], [1.0; 2]),
            JONES_TOTAL,
            CoherentDiffractionExposurePixelLimits {
                spatial: CoherentDiffractionPixelLimits {
                    maximum_samples: 0,
                    maximum_pupil_evaluations: 0,
                    ..PIXEL_LIMITS
                },
                maximum_pair_contractions: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(empty.energy_j, 0.0);
    assert_eq!(empty.pair_contractions, 0);
}

const EXPOSURE_LIMITS: CoherentDiffractionExposureLimits = CoherentDiffractionExposureLimits {
    maximum_tones: 16,
    sources: LIMITS,
    maximum_pair_contractions: 1000,
};

const EXPOSURE_PIXEL_LIMITS: CoherentDiffractionExposurePixelLimits =
    CoherentDiffractionExposurePixelLimits {
        spatial: PIXEL_LIMITS,
        maximum_pair_contractions: 1_000_000,
    };

fn tone<'a>(
    frequency_hz: f64,
    mean: &'a [CoherentPupilWave],
    factors: &'a [PupilCoherenceFactor<'a>],
) -> CoherentDiffractionTone<'a> {
    CoherentDiffractionTone {
        frequency_hz,
        illumination: CoherentPupilIllumination { mean, factors },
    }
}

fn sinc_reference(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

#[test]
fn g1_spectral_shutter_matches_beat_visibility_and_absolute_sensor_phase() {
    let f0 = 299_792_458.0 / LAMBDA;
    let distance = 0.037;
    let waves_a = [CoherentPupilWave {
        tangent_direction: [0.03, 0.0],
        amplitude_sqrt_w: C64::new(0.8, -0.2),
    }];
    let waves_b = [CoherentPupilWave {
        tangent_direction: [0.0299, 0.00015],
        amplitude_sqrt_w: C64::new(-0.3, 0.9),
    }];
    let at = [0.0302, -0.0001];
    for camera in [focused(), defocused(2.4)] {
        for delta in [1.0e9, 0.13 * f0] {
            let f1 = f0 + delta;
            let beat = f1 - f0;
            let tones = [tone(f0, &waves_a, &[]), tone(f1, &waves_b, &[])];
            let amplitude = |waves: &[CoherentPupilWave], frequency: f64| {
                let offset =
                    std::array::from_fn(|axis| at[axis] - waves[0].tangent_direction[axis]);
                let wavelength = 299_792_458.0 / frequency * 1e9;
                let pupil = if camera.defocused_diffraction().is_some() {
                    camera.defocused_diffraction_amplitude(offset, wavelength)
                } else {
                    camera.in_focus_diffraction_amplitude(offset, wavelength)
                }
                .unwrap();
                pupil * waves[0].amplitude_sqrt_w
            };
            let a = amplitude(&waves_a, f0);
            let b = amplitude(&waves_b, f1);
            for (start, duration) in [
                (0.0, 0.25 / beat),
                (0.13 / beat, 1.0 / beat),
                (-0.2 / beat, 3.7 / beat),
            ] {
                let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
                let query = camera
                    .prepare_coherent_diffraction_exposure(
                        &tones,
                        distance,
                        shutter,
                        EXPOSURE_LIMITS,
                        || false,
                    )
                    .unwrap();
                let phase = core::f64::consts::TAU
                    * beat
                    * (distance * (1.0 + 0.5 * (at[0] * at[0] + at[1] * at[1])) / 299_792_458.0
                        - start
                        - duration * 0.5);
                let rotation = C64::new(phase.cos(), phase.sin());
                let cross = b * a.conj() * rotation;
                let expected = a.norm_sq()
                    + b.norm_sq()
                    + 2.0 * cross.re * sinc_reference(core::f64::consts::PI * beat * duration);
                let measured = query.observe(at, || false).unwrap();
                // The widely separated lines exercise independent wavelength and
                // defocus transfers; their independent f64 reference phase has
                // a ~1e-11-radian argument-rounding floor in this fixture.
                close(measured.power_density_w, expected, 3e-10);
                assert_eq!(measured.fluctuation_power_density_w, 0.0);
            }
        }
    }
}

#[test]
fn g3_spectral_coherence_equal_lines_and_incoherent_banks_keep_their_normalization() {
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = f0 + 1e9;
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.2e-9).unwrap();
    let positive = [wave([0.0; 2], C64::ONE)];
    let negative = [wave([0.0; 2], C64::new(-1.0, 0.0))];
    let zero: [CoherentPupilWave; 0] = [];
    let coherent = [tone(f0, &positive, &[]), tone(f1, &negative, &[])];
    let power = camera
        .prepare_coherent_diffraction_exposure(&coherent, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap()
        .power_density_w;
    let factor_a = [PupilCoherenceFactor { waves: &positive }];
    let factor_b = [PupilCoherenceFactor { waves: &negative }];
    let same_factor = [tone(f0, &[], &factor_a), tone(f1, &[], &factor_b)];
    let measured = camera
        .prepare_coherent_diffraction_exposure(&same_factor, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap();
    close(measured.fluctuation_power_density_w, power, 2e-13);
    assert_eq!(measured.mean_power_density_w, 0.0);
    let factor_a = [
        PupilCoherenceFactor { waves: &positive },
        PupilCoherenceFactor { waves: &zero },
    ];
    let factor_b = [
        PupilCoherenceFactor { waves: &zero },
        PupilCoherenceFactor { waves: &negative },
    ];
    let separate = [tone(f0, &[], &factor_a), tone(f1, &[], &factor_b)];
    let incoherent = camera
        .prepare_coherent_diffraction_exposure(&separate, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap()
        .power_density_w;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let expected = area * ((f0 / 299_792_458.0).powi(2) + (f1 / 299_792_458.0).powi(2));
    close(incoherent, expected, 3e-12);
    assert!(power < 0.5 * incoherent);

    // A one-Hz separation beside an optical carrier, observed for 0.1 ns.
    // The beat energy is far below one binary64 step of either isolated power;
    // converting each window kernel to f64 before contraction would erase it.
    let tiny_duration = 1e-10;
    let tiny_shutter = CoherentDiffractionShutter::try_new(
        0.03 / 299_792_458.0 - 0.5 * tiny_duration,
        tiny_duration,
    )
    .unwrap();
    let equal_opposite = [wave([0.0; 2], C64::new(-f0 / (f0 + 1.0), 0.0))];
    let tiny_tones = [
        tone(f0, &positive, &[]),
        tone(f0 + 1.0, &equal_opposite, &[]),
    ];
    let tiny = camera
        .prepare_coherent_diffraction_exposure(
            &tiny_tones,
            0.03,
            tiny_shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap()
        .power_density_w;
    let leading =
        area * (f0 / 299_792_458.0).powi(2) * (core::f64::consts::TAU * tiny_duration).powi(2)
            / 12.0;
    assert!(tiny > 0.0);
    close(tiny, leading, 2e-10);

    // Equal-frequency fields must coalesce before any diagonal/cross moments.
    let strong = [wave([0.0; 2], C64::new(1e16, -2e16))];
    let weak = [wave([0.0; 2], C64::new(0.5, 0.75))];
    let opposite = [wave([0.0; 2], C64::new(-1e16, 2e16))];
    let entries = [
        tone(f0, &strong, &[]),
        tone(f0, &weak, &[]),
        tone(f0, &opposite, &[]),
    ];
    let query = camera
        .prepare_coherent_diffraction_exposure(&entries, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    assert_eq!(query.frequency_count(), 1);
    close(
        query.observe([0.0; 2], || false).unwrap().power_density_w,
        intensity(&camera, &weak, &[], [0.0; 2]),
        3e-12,
    );
    let footprint = pixel([-2.0; 2], [2.0; 2]);
    let collected = query
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    close(
        collected.average_power.power_w,
        pixel_power(&camera, &weak, &[], footprint).power_w,
        3e-12,
    );
    close(
        collected.energy_j,
        collected.average_power.power_w * shutter.duration_s(),
        2e-15,
    );
}

#[test]
fn g1_spectral_pixel_integrates_chromatic_square_pupils_and_sensor_chirp() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let centre = [0.03, 0.0];
    let minimum = [centre[0] - SCALE, -SCALE];
    let maximum = [centre[0] + SCALE, SCALE];
    let footprint = CoherentDiffractionPixel::try_new(minimum, maximum).unwrap();
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = f0 + 1e9;
    let beat = f1 - f0;
    let distance = 1e6;
    let duration = 0.2 / beat;
    let midpoint = distance * (1.0 + 0.5 * centre[0] * centre[0]) / 299_792_458.0;
    let shutter = CoherentDiffractionShutter::try_new(midpoint - duration * 0.5, duration).unwrap();
    let waves = [CoherentPupilWave {
        tangent_direction: centre,
        amplitude_sqrt_w: C64::ONE,
    }];
    let tones = [tone(f0, &waves, &[]), tone(f1, &waves, &[])];
    let query = camera
        .prepare_coherent_diffraction_exposure(&tones, distance, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let measured = query
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    assert!(
        measured.average_power.subdivisions[0] > 100,
        "sensor chirp must influence the spatial work bound"
    );

    // Independent separable square-pupil Fourier integral, evaluated with
    // Simpson in each axis. It uses no production transfer or shutter kernel.
    let a = RADIUS * core::f64::consts::FRAC_1_SQRT_2;
    let lambda = [299_792_458.0 / f0, 299_792_458.0 / f1];
    let delta_k = core::f64::consts::TAU * beat / 299_792_458.0;
    let shape = |line: usize, axis: usize, coordinate: f64| {
        sinc_reference(core::f64::consts::TAU * a * (coordinate - centre[axis]) / lambda[line])
    };
    let mut diagonal = 0.0;
    for line in 0..2 {
        let mut norm = (2.0 * a / lambda[line]).powi(2);
        for axis in 0..2 {
            norm *= integrate_reference(minimum[axis], maximum[axis], |q| {
                shape(line, axis, q).powi(2)
            });
        }
        diagonal += norm;
    }
    let mut cross = C64::ONE;
    for axis in 0..2 {
        let real = integrate_reference(minimum[axis], maximum[axis], |q| {
            shape(0, axis, q)
                * shape(1, axis, q)
                * (0.5 * delta_k * distance * (q * q - centre[axis] * centre[axis])).cos()
        });
        let imaginary = integrate_reference(minimum[axis], maximum[axis], |q| {
            shape(0, axis, q)
                * shape(1, axis, q)
                * (0.5 * delta_k * distance * (q * q - centre[axis] * centre[axis])).sin()
        });
        cross = cross * C64::new(real, imaginary);
    }
    let reference_phase =
        core::f64::consts::TAU * beat * (midpoint - shutter.start_s() - 0.5 * duration);
    cross = cross * C64::new(reference_phase.cos(), reference_phase.sin());
    let expected = diagonal
        + 8.0 * a * a / (lambda[0] * lambda[1])
            * cross.re
            * sinc_reference(core::f64::consts::PI * beat * duration);
    close(measured.average_power.power_w, expected, 3e-8);
    close(measured.energy_j, expected * duration, 3e-8);
    let refined = query
        .integrate_pixel(
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
    close(
        refined.average_power.power_w,
        measured.average_power.power_w,
        2e-11,
    );
}

#[test]
fn g0_g4_spectral_exposure_refuses_invalid_models_and_preflights_work() {
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let waves = [wave([0.0; 2], C64::ONE)];
    let tones = [tone(f0, &waves, &[]), tone(f0 + 1e9, &waves, &[])];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.2e-9).unwrap();
    let query = camera
        .prepare_coherent_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let footprint = pixel([-1.0; 2], [1.0; 2]);
    let baseline = query
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        query.integrate_pixel(
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
        query.integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || {
            polls += 1;
            polls == 31
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query
            .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
            .unwrap(),
        baseline
    );
    for frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            camera.prepare_coherent_diffraction_exposure(
                &[tone(frequency, &waves, &[])],
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false
            ),
            Err(CoherentDiffractionError::InvalidFrequency)
        ));
    }
    assert!(matches!(
        camera.prepare_coherent_diffraction_exposure(
            &[tone(1.0, &waves, &[])],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::Camera(_))
    ));
    let factors = [PupilCoherenceFactor { waves: &waves }];
    assert!(matches!(
        camera.prepare_coherent_diffraction_exposure(
            &[tones[0], tone(f0, &waves, &factors)],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::InconsistentCoherenceFactors)
    ));
    assert!(matches!(
        camera.prepare_coherent_diffraction_exposure(
            &tones,
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                sources: CoherentDiffractionLimits {
                    maximum_waves: 1,
                    ..LIMITS
                },
                ..EXPOSURE_LIMITS
            },
            || false
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pupil waves",
            ..
        })
    ));
    assert!(CoherentDiffractionShutter::try_new(0.0, 0.0).is_err());
    assert!(CoherentDiffractionShutter::try_new(f64::INFINITY, 1.0).is_err());
    let physical = PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        2.0,
        Aperture::try_circular(RADIUS).unwrap(),
    )
    .unwrap()
    .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
    .unwrap();
    assert!(matches!(
        physical.prepare_coherent_diffraction_exposure(
            &tones,
            0.05,
            shutter,
            EXPOSURE_LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::InvalidSensorPlane)
    ));
    physical
        .prepare_coherent_diffraction_exposure(
            &tones,
            0.05 / (1.0 - 0.05 / 2.0),
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let empty = camera
        .prepare_coherent_diffraction_exposure(&[], 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    assert_eq!(
        empty
            .integrate_pixel(
                footprint,
                CoherentDiffractionExposurePixelLimits {
                    spatial: CoherentDiffractionPixelLimits {
                        maximum_samples: 0,
                        maximum_pupil_evaluations: 0,
                        ..PIXEL_LIMITS
                    },
                    maximum_pair_contractions: 0,
                },
                || false
            )
            .unwrap()
            .energy_j,
        0.0
    );
}

#[path = "coherent_diffraction/polarization.rs"]
mod polarization;
#[path = "coherent_diffraction/exposure_coherency.rs"]
mod exposure_coherency;

const RADIUS: f64 = 0.001;
const LAMBDA_NM: f64 = 550.0;
const LAMBDA: f64 = LAMBDA_NM * 1e-9;
const SCALE: f64 = LAMBDA / (core::f64::consts::TAU * RADIUS);
const LIMITS: CoherentDiffractionLimits = CoherentDiffractionLimits {
    maximum_waves: 32,
    maximum_factors: 8,
};

fn base(aperture: Aperture) -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_vertical_fov(0.5).unwrap(),
        2.0,
        aperture,
    )
    .unwrap()
}

fn focused() -> PhysicalCamera {
    base(Aperture::try_circular(RADIUS).unwrap())
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
        .unwrap()
}

fn defocused(phase: f64) -> PhysicalCamera {
    let depth = 1.0 / (0.5 + phase * LAMBDA / (core::f64::consts::PI * RADIUS * RADIUS));
    base(Aperture::try_circular(RADIUS).unwrap())
        .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(depth), 0.0).unwrap())
        .unwrap()
}

fn wave(q: [f64; 2], a: C64) -> CoherentPupilWave {
    CoherentPupilWave {
        tangent_direction: q.map(|x| x * SCALE),
        amplitude_sqrt_w: a,
    }
}

fn intensity(
    camera: &PhysicalCamera,
    mean: &[CoherentPupilWave],
    factors: &[PupilCoherenceFactor<'_>],
    at: [f64; 2],
) -> f64 {
    camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination { mean, factors },
            LIMITS,
            || false,
        )
        .unwrap()
        .observe(at, || false)
        .unwrap()
        .power_density_w
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-20),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn g1_same_pupil_fields_add_before_squaring_and_covariance_columns_add_power() {
    let camera = focused();
    let unit = wave([0.0, 0.0], C64::ONE);
    let one = core::f64::consts::PI * RADIUS * RADIUS / (LAMBDA * LAMBDA);
    close(intensity(&camera, &[unit], &[], [0.0; 2]), one, 2e-12);
    close(
        intensity(&camera, &[unit, unit], &[], [0.0; 2]),
        4.0 * one,
        2e-12,
    );
    let opposite = wave([0.0; 2], C64::new(-1.0, 0.0));
    assert_eq!(intensity(&camera, &[unit, opposite], &[], [0.0; 2]), 0.0);
    let quadrature = wave([0.0; 2], C64::new(0.0, 1.0));
    close(
        intensity(&camera, &[unit, quadrature], &[], [0.0; 2]),
        2.0 * one,
        2e-12,
    );
    let a = [unit];
    let b = [opposite];
    let factors = [
        PupilCoherenceFactor { waves: &a },
        PupilCoherenceFactor { waves: &b },
    ];
    close(
        intensity(&camera, &[], &factors, [0.0; 2]),
        2.0 * one,
        2e-12,
    );
    close(
        intensity(&camera, &[unit], &factors, [0.0; 2]),
        3.0 * one,
        2e-12,
    );
}

#[test]
fn g3_airy_side_lobe_sign_and_weak_residual_survive_field_summation() {
    let camera = focused();
    // At q=5 the circular transfer is in its negative first side lobe.
    // Using sqrt(intensity) would turn the cross term constructive.
    let waves = [wave([0.0; 2], C64::ONE), wave([5.0, 0.0], C64::ONE)];
    let on_axis = camera
        .in_focus_diffraction_amplitude([0.0; 2], LAMBDA_NM)
        .unwrap();
    let lobe = camera
        .in_focus_diffraction_amplitude([-5.0 * SCALE, 0.0], LAMBDA_NM)
        .unwrap();
    assert!(lobe.re < 0.0);
    let coherent = intensity(&camera, &waves, &[], [0.0; 2]);
    close(coherent, (on_axis + lobe).norm_sq(), 3e-12);
    assert!(coherent < on_axis.norm_sq() + lobe.norm_sq());

    let strong = wave([0.2, -0.1], C64::new(1e16, -2e16));
    let small = wave([0.2, -0.1], C64::new(0.5, 0.75));
    let cancel = wave([0.2, -0.1], C64::new(-1e16, 2e16));
    let at = [0.5 * SCALE, -0.3 * SCALE];
    close(
        intensity(&camera, &[strong, small, cancel], &[], at),
        intensity(&camera, &[small], &[], at),
        2e-12,
    );
}

#[test]
fn g3_common_phase_and_unitary_covariance_factor_changes_preserve_the_image() {
    let camera = defocused(2.4);
    let a = wave([-0.7, 0.2], C64::new(0.8, -0.3));
    let b = wave([1.1, -0.4], C64::new(-0.1, 0.6));
    let common_phase = C64::new(0.6, 0.8);
    let rotated = [a, b].map(|mut v| {
        v.amplitude_sqrt_w = v.amplitude_sqrt_w * common_phase;
        v
    });
    let aa = [a];
    let bb = [b];
    let diagonal = [
        PupilCoherenceFactor { waves: &aa },
        PupilCoherenceFactor { waves: &bb },
    ];
    let scale = core::f64::consts::FRAC_1_SQRT_2;
    let plus = [a, b].map(|mut v| {
        v.amplitude_sqrt_w = v.amplitude_sqrt_w.scale(scale);
        v
    });
    let minus = [
        plus[0],
        CoherentPupilWave {
            amplitude_sqrt_w: plus[1].amplitude_sqrt_w.scale(-1.0),
            ..plus[1]
        },
    ];
    let mixed = [
        PupilCoherenceFactor { waves: &plus },
        PupilCoherenceFactor { waves: &minus },
    ];
    for q in [[0.0, 0.0], [0.3, -0.6], [2.1, 1.0], [5.0, -2.0]] {
        let at = q.map(|x| x * SCALE);
        close(
            intensity(&camera, &rotated, &[], at),
            intensity(&camera, &[a, b], &[], at),
            3e-12,
        );
        close(
            intensity(&camera, &[], &mixed, at),
            intensity(&camera, &[], &diagonal, at),
            3e-12,
        );
    }
}

// Independent midpoint area integral of the actual authored pupil field.
// t=r² is uniform in disk area. No production PSF, Bessel function, or
// amplitude/phase kernel is used by this oracle.
fn disk_field(waves: &[CoherentPupilWave], at: [f64; 2], phase: f64) -> C64 {
    let radial = 192;
    let azimuth = 256;
    let mut sum = C64::ZERO;
    for i in 0..radial {
        let t = (i as f64 + 0.5) / radial as f64;
        for j in 0..azimuth {
            let angle = core::f64::consts::TAU * (j as f64 + 0.5) / azimuth as f64;
            let u = [t.sqrt() * angle.cos(), t.sqrt() * angle.sin()];
            for w in waves {
                let argument = phase * t
                    - ((at[0] - w.tangent_direction[0]) * u[0]
                        + (at[1] - w.tangent_direction[1]) * u[1])
                        / SCALE;
                sum = sum + w.amplitude_sqrt_w * C64::new(argument.cos(), argument.sin());
            }
        }
    }
    sum.scale(core::f64::consts::PI.sqrt() * RADIUS / LAMBDA / (radial * azimuth) as f64)
}

#[test]
fn g1_tilted_wave_interference_matches_the_original_complex_pupil_integral() {
    let waves = [
        wave([-1.3, 0.2], C64::new(0.8, -0.3)),
        wave([0.9, -0.4], C64::new(-0.4, 0.7)),
        wave([0.1, 0.8], C64::new(0.2, 0.4)),
    ];
    let at = [0.2 * SCALE, -0.3 * SCALE];
    for phase in [0.0, 2.1, -2.1] {
        let camera = if phase == 0.0 {
            focused()
        } else {
            defocused(phase)
        };
        let prepared = camera
            .prepare_coherent_diffraction(
                LAMBDA_NM,
                CoherentPupilIllumination {
                    mean: &waves,
                    factors: &[],
                },
                LIMITS,
                || false,
            )
            .unwrap();
        let actual = prepared.observe(at, || false).unwrap();
        let reference = disk_field(&waves, at, phase);
        assert!(
            (actual.mean_field - reference).abs() < 3e-5 * reference.abs(),
            "phase {phase}: actual {:?}, pupil integral {reference:?}",
            actual.mean_field
        );
        close(actual.power_density_w, reference.norm_sq(), 6e-5);
    }
}

#[test]
fn g0_g4_original_admission_work_limits_and_cancellation_remain_effective() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE), wave([1.0, 0.0], C64::ONE)];
    let input = CoherentPupilIllumination {
        mean: &waves,
        factors: &[],
    };
    let limited = CoherentDiffractionLimits {
        maximum_waves: 1,
        ..LIMITS
    };
    assert!(matches!(
        camera.prepare_coherent_diffraction(LAMBDA_NM, input, limited, || false),
        Err(CoherentDiffractionError::Limit {
            quantity: "pupil waves",
            requested: 2,
            maximum: 1
        })
    ));
    let bad = [CoherentPupilWave {
        amplitude_sqrt_w: C64::new(f64::NAN, 0.0),
        ..waves[0]
    }];
    assert!(matches!(
        camera.prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &bad,
                factors: &[]
            },
            LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::NonFiniteAmplitude)
    ));
    assert!(matches!(
        camera.prepare_coherent_diffraction(0.0, input, LIMITS, || false),
        Err(CoherentDiffractionError::Camera(
            CameraError::InvalidDiffraction
        ))
    ));
    assert!(matches!(
        camera.prepare_coherent_diffraction(LAMBDA_NM, input, LIMITS, || true),
        Err(CoherentDiffractionError::Cancelled)
    ));
    let prepared = camera
        .prepare_coherent_diffraction(LAMBDA_NM, input, LIMITS, || false)
        .unwrap();
    assert_eq!(prepared.wave_count(), 2);
    assert_eq!(prepared.factor_count(), 0);
    assert_eq!(prepared.wavelength_nm(), LAMBDA_NM);
    let baseline = prepared.observe([0.0; 2], || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        prepared.observe([0.0; 2], || {
            polls += 1;
            polls == 3
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(prepared.observe([0.0; 2], || false).unwrap(), baseline);
    assert_eq!(
        prepared.observe([0.11, 0.0], || false),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
    assert_eq!(
        prepared.observe([f64::NAN, 0.0], || false),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
}

#[test]
fn g0_empty_banks_are_dark_and_finite_overflow_is_refused() {
    let camera = focused();
    let empty = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &[],
                factors: &[],
            },
            CoherentDiffractionLimits {
                maximum_waves: 0,
                maximum_factors: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        empty.observe([0.0; 2], || false).unwrap().power_density_w,
        0.0
    );
    let waves = [wave([0.0; 2], C64::new(1e200, 0.0))];
    let large = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    assert_eq!(
        large.observe([0.0; 2], || false),
        Err(CoherentDiffractionError::Unrepresentable)
    );
}

const PIXEL_LIMITS: CoherentDiffractionPixelLimits = CoherentDiffractionPixelLimits {
    maximum_phase_span_rad: 1.0,
    maximum_samples: 100_000,
    maximum_pupil_evaluations: 500_000,
};

fn pixel(low_q: [f64; 2], high_q: [f64; 2]) -> CoherentDiffractionPixel {
    CoherentDiffractionPixel::try_new(low_q.map(|v| v * SCALE), high_q.map(|v| v * SCALE)).unwrap()
}

fn pixel_power(
    camera: &PhysicalCamera,
    mean: &[CoherentPupilWave],
    factors: &[PupilCoherenceFactor<'_>],
    footprint: CoherentDiffractionPixel,
) -> CoherentDiffractionPixelPower {
    camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination { mean, factors },
            LIMITS,
            || false,
        )
        .unwrap()
        .integrate_pixel(footprint, PIXEL_LIMITS, || false)
        .unwrap()
}

// Independent 1-D composite Simpson rule, used only for the separable analytic
// square-pupil sinc transfer. It shares no production pupil or pixel kernel.
fn integrate_reference(low: f64, high: f64, integrand: impl Fn(f64) -> f64) -> f64 {
    let n = 8192;
    let step = (high - low) / n as f64;
    let mut sum = integrand(low) + integrand(high);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * integrand(low + step * i as f64);
    }
    sum * step / 3.0
}

#[test]
fn g1_finite_pixel_resolves_a_dark_centre_and_matches_square_pupil_power() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let waves = [
        wave([-2.2, 0.0], C64::ONE),
        wave([2.2, 0.0], C64::new(-1.0, 0.0)),
    ];
    let footprint = pixel([-7.1, -4.3], [8.3, 2.8]);
    let prepared = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let result = prepared
        .integrate_pixel(footprint, PIXEL_LIMITS, || false)
        .unwrap();
    let sinc = |q: f64| {
        let phase = q * core::f64::consts::FRAC_1_SQRT_2;
        if phase == 0.0 {
            1.0
        } else {
            phase.sin() / phase
        }
    };
    let x = integrate_reference(-7.1, 8.3, |q| (sinc(q + 2.2) - sinc(q - 2.2)).powi(2));
    let y = integrate_reference(-4.3, 2.8, |q| sinc(q).powi(2));
    // sqrt(area)/lambda for the square, times the tangent-to-q Jacobian.
    let expected = x * y / (2.0 * core::f64::consts::PI.powi(2));
    close(result.power_w, expected, 3e-10);
    assert!(result.power_w > 1.0);
    let centre = prepared.observe([0.0; 2], || false).unwrap();
    assert!(centre.power_density_w * footprint.area_tangent() < 1e-20 * result.power_w);
    assert_eq!(result.pupil_evaluations, 2 * result.quadrature_samples);
    assert_eq!(result.fluctuation_power_w, 0.0);
    let refined = prepared
        .integrate_pixel(
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_phase_span_rad: 0.5,
                ..PIXEL_LIMITS
            },
            || false,
        )
        .unwrap();
    close(refined.power_w, expected, 5e-12);
    assert!(refined.quadrature_samples > result.quadrature_samples);
}

#[test]
fn g3_pixel_power_preserves_coherent_nulls_and_independent_factor_powers() {
    let footprint = pixel([-3.5, -2.5], [4.5, 3.5]);
    let unit = [wave([0.2, -0.1], C64::ONE)];
    let opposite = [wave([0.2, -0.1], C64::new(-1.0, 0.0))];
    let factors = [
        PupilCoherenceFactor { waves: &unit },
        PupilCoherenceFactor { waves: &opposite },
    ];
    for camera in [focused(), defocused(2.4)] {
        let one = pixel_power(&camera, &unit, &[], footprint).power_w;
        let constructive = pixel_power(&camera, &[unit[0], unit[0]], &[], footprint);
        close(constructive.power_w, 4.0 * one, 2e-13);
        let dark = pixel_power(&camera, &[unit[0], opposite[0]], &[], footprint);
        assert_eq!(dark.power_w, 0.0);
        let mixed = pixel_power(&camera, &unit, &factors, footprint);
        close(mixed.mean_power_w, one, 2e-13);
        close(mixed.fluctuation_power_w, 2.0 * one, 2e-13);
        close(mixed.power_w, 3.0 * one, 2e-13);
    }
}

#[test]
fn g3_adjacent_pixels_add_and_weak_interference_survives_area_integration() {
    let camera = base(Aperture::try_annular(RADIUS, 0.4 * RADIUS).unwrap())
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
        .unwrap();
    let waves = [
        wave([-0.7, 0.2], C64::new(0.8, -0.3)),
        wave([1.1, -0.4], C64::new(-0.1, 0.6)),
    ];
    let whole = pixel_power(&camera, &waves, &[], pixel([-4.0, -2.0], [4.0, 2.0]));
    let left = pixel_power(&camera, &waves, &[], pixel([-4.0, -2.0], [0.0, 2.0]));
    let right = pixel_power(&camera, &waves, &[], pixel([0.0, -2.0], [4.0, 2.0]));
    close(whole.power_w, left.power_w + right.power_w, 2e-10);
    let strong = wave([0.2, -0.1], C64::new(1e16, -2e16));
    let weak = wave([0.2, -0.1], C64::new(0.5, 0.75));
    let cancel = wave([0.2, -0.1], C64::new(-1e16, 2e16));
    let footprint = pixel([-2.0; 2], [2.0; 2]);
    close(
        pixel_power(&camera, &[strong, weak, cancel], &[], footprint).power_w,
        pixel_power(&camera, &[weak], &[], footprint).power_w,
        3e-12,
    );
}

#[test]
fn g0_g4_pixel_admission_work_preflight_and_cancelled_reuse() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let prepared = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let footprint = pixel([-1.0; 2], [1.0; 2]);
    assert_eq!(footprint.minimum_tangent(), [-SCALE; 2]);
    assert_eq!(footprint.maximum_tangent(), [SCALE; 2]);
    let baseline = prepared
        .integrate_pixel(footprint, PIXEL_LIMITS, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        prepared.integrate_pixel(
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_samples: 0,
                ..PIXEL_LIMITS
            },
            || {
                polls += 1;
                false
            },
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel quadrature samples",
            ..
        })
    ));
    assert_eq!(polls, 1, "work refusal precedes every field sample");
    assert!(matches!(
        prepared.integrate_pixel(
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_pupil_evaluations: baseline.pupil_evaluations - 1,
                ..PIXEL_LIMITS
            },
            || false,
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel pupil evaluations",
            ..
        })
    ));
    let mut polls = 0;
    assert_eq!(
        prepared.integrate_pixel(footprint, PIXEL_LIMITS, || {
            polls += 1;
            polls == 31
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        prepared
            .integrate_pixel(footprint, PIXEL_LIMITS, || false)
            .unwrap(),
        baseline
    );
    for phase in [0.0, -1.0, 1.01, f64::NAN] {
        assert_eq!(
            prepared.integrate_pixel(
                footprint,
                CoherentDiffractionPixelLimits {
                    maximum_phase_span_rad: phase,
                    ..PIXEL_LIMITS
                },
                || false,
            ),
            Err(CoherentDiffractionError::InvalidPixelResolution)
        );
    }
    assert_eq!(
        CoherentDiffractionPixel::try_new([0.0; 2], [0.0, 0.1]),
        Err(CoherentDiffractionError::InvalidPixel)
    );
    assert_eq!(
        CoherentDiffractionPixel::try_new([f64::NAN, 0.0], [0.01; 2]),
        Err(CoherentDiffractionError::InvalidPixel)
    );
    assert_eq!(
        CoherentDiffractionPixel::try_new([0.0; 2], [0.08; 2]),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
    let empty = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &[],
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let dark = empty
        .integrate_pixel(
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_samples: 0,
                maximum_pupil_evaluations: 0,
                ..PIXEL_LIMITS
            },
            || false,
        )
        .unwrap();
    assert_eq!(dark.power_w, 0.0);
    assert_eq!(dark.quadrature_samples, 0);
}

#[path = "coherent_diffraction/pupil.rs"]
mod pupil;
