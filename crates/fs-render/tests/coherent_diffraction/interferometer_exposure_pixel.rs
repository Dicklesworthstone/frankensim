use super::*;
use fs_math::dd::Dd;
use fs_render::camera::{CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits};

const C: f64 = 299_792_458.0;
const F0: f64 = 562_949_953_421_312.0;

fn receiver(
    coordinates: [[f64; 2]; 2],
    distance: f64,
    fractions: [f64; 2],
    path: f64,
    phase: f64,
) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new(coordinates, distance, fractions, path, phase)
        .unwrap()
}

fn wavelength_nm(frequency: f64) -> f64 {
    (Dd::from_f64(C) / Dd::from_f64(frequency) * Dd::from_f64(1e9)).to_f64()
}

#[test]
fn g1_spectral_relay_pixel_matches_monochromatic_power_and_retains_weak_ports() {
    let camera = defocused(1.2);
    let mean = [wave([0.2, -0.1], C64::new(0.8, 0.4))];
    let noise = [wave([-0.6, 0.3], C64::new(-0.2, 0.5))];
    let factors = [PupilCoherenceFactor { waves: &noise }];
    let tones = [tone(F0, &mean, &factors)];
    let shutter = CoherentDiffractionShutter::try_new(-3.0, 0.125).unwrap();
    let prepared = camera
        .prepare_coherent_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let monochromatic = camera
        .prepare_coherent_diffraction(wavelength_nm(F0), tones[0].illumination, LIMITS, || false)
        .unwrap();
    let detector = receiver(
        [[-0.4 * SCALE, 0.1 * SCALE], [0.7 * SCALE, -0.2 * SCALE]],
        0.03,
        [0.35, 0.55],
        -0.17 * C / F0,
        0.31,
    );
    let footprint = pixel([-1.3, -0.7], [1.1, 0.9]);
    let actual = prepared
        .integrate_interferometer_pixel(detector, footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let expected = monochromatic
        .integrate_interferometer_pixel(detector, footprint, PIXEL_LIMITS, || false)
        .unwrap();
    for arm in 0..2 {
        close(
            actual.average_power.arm_power_w[arm],
            expected.arm_power_w[arm],
            3e-12,
        );
        close(
            actual.average_power.output_power_w[arm],
            expected.output_power_w[arm],
            3e-12,
        );
        close(
            actual.arm_energy_j[arm],
            shutter.duration_s() * expected.arm_power_w[arm],
            3e-12,
        );
        close(
            actual.output_energy_j[arm],
            shutter.duration_s() * expected.output_power_w[arm],
            3e-12,
        );
    }
    close(
        actual.average_power.mutual_coherence_w.re,
        expected.mutual_coherence_w.re,
        3e-12,
    );
    close(
        actual.average_power.mutual_coherence_w.im,
        expected.mutual_coherence_w.im,
        3e-12,
    );
    close(
        actual.average_power.fringe_visibility,
        expected.fringe_visibility,
        3e-12,
    );

    let other = [wave([0.5, -0.2], C64::new(0.3, -0.1))];
    let tones = [tone(F0, &mean, &[]), tone(F0 + 1.0, &other, &[])];
    let prepared = camera
        .prepare_coherent_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let native = prepared
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    for phase in [0.0_f64, 1e-20] {
        let actual = prepared
            .integrate_interferometer_pixel(
                receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase),
                footprint,
                EXPOSURE_PIXEL_LIMITS,
                || false,
            )
            .unwrap();
        close(
            actual.average_power.output_power_w[0],
            native.average_power.power_w,
            3e-12,
        );
        if phase == 0.0 {
            assert_eq!(actual.average_power.output_power_w[1], 0.0);
        } else {
            let expected = native.average_power.power_w * (0.5 * phase).sin().powi(2);
            assert!(actual.average_power.output_power_w[1] > 0.0);
            close(
                actual.average_power.output_power_w[1] / expected,
                1.0,
                3e-12,
            );
        }
    }
}

#[test]
fn g1_space_and_shutter_average_matches_separable_square_pupil_reference() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let frequencies = [F0, F0 + 4_294_967_296.0];
    let delta = frequencies[1] - frequencies[0];
    let amplitudes = [C64::ONE, C64::new(0.7, 0.2)];
    let waves = amplitudes.map(|amplitude| [wave([0.0; 2], amplitude)]);
    let tones = [
        tone(frequencies[0], &waves[0], &[]),
        tone(frequencies[1], &waves[1], &[]),
    ];
    let centers_q = [[-0.7, 0.2], [0.7, -0.2]];
    let centers = centers_q.map(|q| q.map(|value| value * SCALE));
    let distance = 30.0;
    let fractions = [0.4_f64, 0.6];
    let path = 0.13 * C / F0;
    let phase = -0.23;
    let detector = receiver(centers, distance, fractions, path, phase);
    let low = [-2.5, -1.1];
    let high = [2.5, 1.3];
    let footprint = pixel(low, high);
    let start = distance / C - 0.0625 / delta;
    let mut observed = [0.0; 2];
    for (case, duration) in [0.125 / delta, 1.0 / delta].into_iter().enumerate() {
        let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
        let prepared = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                distance,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let actual = prepared
            .integrate_interferometer_pixel(detector, footprint, EXPOSURE_PIXEL_LIMITS, || false)
            .unwrap();

        // Independent analytic square-pupil transform, analytic uniform-shutter
        // mean, and composite Simpson area integration. Production pupil,
        // temporal kernels and spatial Gauss nodes are not used by this oracle.
        let moment = |left: usize, right: usize| {
            let mut total = C64::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    let fi = frequencies[i];
                    let fj = frequencies[j];
                    let beat = fi - fj;
                    let cycles = beat * duration;
                    if cycles != 0.0 && cycles.fract() == 0.0 {
                        continue;
                    }
                    let temporal = sinc_reference(core::f64::consts::PI * cycles);
                    let piston = core::f64::consts::TAU
                        * (beat * (distance / C - start - 0.5 * duration)
                            + (fi * left as f64 - fj * right as f64) * path / C)
                        + (left as f64 - right as f64) * phase;
                    let scale = 2.0 * RADIUS * RADIUS * fi * fj / (C * C)
                        * SCALE
                        * SCALE
                        * (fractions[left] * fractions[right]).sqrt();
                    let mut term = amplitudes[i]
                        * amplitudes[j].conj()
                        * C64::new(piston.cos(), piston.sin()).scale(scale * temporal);
                    for axis in 0..2 {
                        let integrand = |q: f64| {
                            let x = centers[left][axis] + q * SCALE;
                            let y = centers[right][axis] + q * SCALE;
                            let shape = sinc_reference(
                                core::f64::consts::SQRT_2 * core::f64::consts::PI * RADIUS * fi / C
                                    * x,
                            ) * sinc_reference(
                                core::f64::consts::SQRT_2 * core::f64::consts::PI * RADIUS * fj / C
                                    * y,
                            );
                            let chirp =
                                core::f64::consts::PI * distance / C * (fi * x * x - fj * y * y);
                            C64::new(chirp.cos(), chirp.sin()).scale(shape)
                        };
                        term = term
                            * C64::new(
                                integrate_reference(low[axis], high[axis], |q| integrand(q).re),
                                integrate_reference(low[axis], high[axis], |q| integrand(q).im),
                            );
                    }
                    total = total + term;
                }
            }
            total
        };
        let arms = [moment(0, 0).re, moment(1, 1).re];
        let mutual = moment(0, 1);
        let outputs = [
            0.5 * (arms[0] + arms[1]) + mutual.re,
            0.5 * (arms[0] + arms[1]) - mutual.re,
        ];
        for arm in 0..2 {
            close(actual.average_power.arm_power_w[arm], arms[arm], 2e-9);
            close(actual.average_power.output_power_w[arm], outputs[arm], 2e-9);
            close(actual.output_energy_j[arm], outputs[arm] * duration, 2e-9);
        }
        close(actual.average_power.mutual_coherence_w.re, mutual.re, 2e-8);
        close(actual.average_power.mutual_coherence_w.im, mutual.im, 2e-8);
        close(
            actual.average_power.fringe_visibility,
            2.0 * mutual.abs() / (arms[0] + arms[1]),
            2e-8,
        );
        assert!(actual.average_power.fringe_visibility < 0.2);
        assert_eq!(
            actual.pair_contractions,
            18 * actual.average_power.quadrature_samples
        );
        observed[case] = actual.average_power.output_power_w[0];
        if case == 0 {
            let refined = prepared
                .integrate_interferometer_pixel(
                    detector,
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
            assert!(
                refined.average_power.quadrature_samples > actual.average_power.quadrature_samples
            );
            close(refined.average_power.output_power_w[0], outputs[0], 2e-10);
            close(refined.average_power.mutual_coherence_w.re, mutual.re, 2e-9);
        }
    }
    assert!((observed[0] - observed[1]).abs() > 0.2 * observed[1]);
}

#[test]
fn g3_spectral_pixels_add_independent_banks_and_partitions_and_preserve_long_phase() {
    let camera = focused();
    let a0 = [wave([-0.3, 0.1], C64::new(1.0, -0.1))];
    let a1 = [wave([0.2, -0.4], C64::new(0.4, 0.3))];
    let b0 = [wave([0.7, -0.1], C64::new(-0.2, 0.4))];
    let b1 = [wave([-0.4, 0.6], C64::new(0.3, -0.2))];
    let factors0 = [
        PupilCoherenceFactor { waves: &a0 },
        PupilCoherenceFactor { waves: &b0 },
    ];
    let factors1 = [
        PupilCoherenceFactor { waves: &a1 },
        PupilCoherenceFactor { waves: &b1 },
    ];
    let mixed = [tone(F0, &[], &factors0), tone(F0 + 1.0, &[], &factors1)];
    let separate = [
        [tone(F0, &a0, &[]), tone(F0 + 1.0, &a1, &[])],
        [tone(F0, &b0, &[]), tone(F0 + 1.0, &b1, &[])],
    ];
    let footprint = pixel([-1.3, -0.8], [1.2, 0.9]);
    let centers = [[-0.4 * SCALE, 0.2 * SCALE], [0.6 * SCALE, -0.1 * SCALE]];
    let collect = |tones: &[CoherentDiffractionTone<'_>], region, epoch: f64, path: f64| {
        camera
            .prepare_coherent_diffraction_exposure(
                tones,
                0.03,
                CoherentDiffractionShutter::try_new(epoch, 0.125).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap()
            .integrate_interferometer_pixel(
                receiver(centers, 0.03, [0.4, 0.6], path, 0.2),
                region,
                EXPOSURE_PIXEL_LIMITS,
                || false,
            )
            .unwrap()
    };
    let values = |result: fs_render::camera::CoherentInterferometerPixelExposure| {
        let p = result.average_power;
        [
            p.arm_power_w[0],
            p.arm_power_w[1],
            p.output_power_w[0],
            p.output_power_w[1],
            p.mutual_coherence_w.re,
            p.mutual_coherence_w.im,
        ]
    };
    let expected = values(collect(&mixed, footprint, 0.0, 0.25 * C));
    for parts in [
        [
            collect(&separate[0], footprint, 0.0, 0.25 * C),
            collect(&separate[1], footprint, 0.0, 0.25 * C),
        ],
        [
            collect(&mixed, pixel([-1.3, -0.8], [0.0, 0.9]), 0.0, 0.25 * C),
            collect(&mixed, pixel([0.0, -0.8], [1.2, 0.9]), 0.0, 0.25 * C),
        ],
    ] {
        let [a, b] = parts.map(values);
        for i in 0..6 {
            close(expected[i], a[i] + b[i], 2e-10);
        }
    }
    for sign in [-1.0_f64, 1.0] {
        let shifted = values(collect(
            &mixed,
            footprint,
            sign * 1_099_511_627_776.0,
            (sign * 1_048_576.0 + 0.25) * C,
        ));
        for i in 0..6 {
            close(shifted[i], expected[i], 2e-12);
        }
    }
}

#[test]
fn g1_sub_ulp_relay_offsets_retain_close_tone_sensor_chirp() {
    let camera = focused();
    let center = [0.0625, 0.0];
    let wave = [CoherentPupilWave {
        tangent_direction: center,
        amplitude_sqrt_w: C64::ONE,
    }];
    let tones = [tone(F0, &wave, &[]), tone(F0 + 1.0, &wave, &[])];
    let distance = C * 2.0_f64.powi(62);
    let half_width = 2.0_f64.powi(-60);
    assert_eq!(center[0] + half_width, center[0]);
    assert_eq!(center[0] - half_width, center[0]);
    let shutter = CoherentDiffractionShutter::try_new(-0.125, 0.25).unwrap();
    let prepared = camera
        .prepare_coherent_diffraction_exposure(&tones, distance, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let footprint = CoherentDiffractionPixel::try_new([-half_width; 2], [half_width; 2]).unwrap();
    let actual = prepared
        .integrate_interferometer_pixel(
            receiver([center; 2], distance, [0.5; 2], 0.0, 0.0),
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    // d/c=2^62 and center=2^-4 put the centre beat phase at an exact integer
    // cycle. Across the x half-width the one-Hz beat sweeps pi/2. The square
    // pupil coordinate changes are sub-ulp and its envelope is constant here;
    // the residual quadratic phase is below 1e-16 radians.
    let field =
        [F0, F0 + 1.0].map(|frequency| core::f64::consts::PI.sqrt() * RADIUS * frequency / C);
    let spatial = sinc_reference(core::f64::consts::FRAC_PI_2);
    let temporal = sinc_reference(core::f64::consts::FRAC_PI_4);
    let expected_density =
        field[0] * field[0] + field[1] * field[1] + 2.0 * field[0] * field[1] * spatial * temporal;
    close(
        actual.average_power.output_power_w[0] / (footprint.area_tangent() * expected_density),
        1.0,
        5e-11,
    );
    assert_eq!(actual.average_power.output_power_w[1], 0.0);
    assert_eq!(actual.average_power.subdivisions, [4, 1]);
    let point = prepared
        .observe_interferometer(
            receiver([center; 2], distance, [0.5; 2], 0.0, 0.0),
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: 4,
                maximum_pair_contractions: 18,
            },
            || false,
        )
        .unwrap();
    assert!(
        actual.average_power.output_power_w[0]
            < 0.9 * footprint.area_tangent() * point.output_power_density_w[0]
    );
}

#[test]
fn g0_g4_spectral_relay_pixel_preflights_all_work_and_reuses_after_cancellation() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let tones = [tone(F0, &waves, &[]), tone(F0 + 1.0, &waves, &[])];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.125).unwrap();
    let prepared = camera
        .prepare_coherent_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let footprint = pixel([-0.7; 2], [0.8; 2]);
    let centers = [[0.0; 2], [0.7 * SCALE, 0.2 * SCALE]];
    let detector = receiver(centers, 0.03, [0.5; 2], 0.0, 0.4);
    let expected = prepared
        .integrate_interferometer_pixel(detector, footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let exact = CoherentDiffractionExposurePixelLimits {
        spatial: CoherentDiffractionPixelLimits {
            maximum_samples: expected.average_power.quadrature_samples,
            maximum_pupil_evaluations: expected.average_power.pupil_evaluations,
            ..PIXEL_LIMITS
        },
        maximum_pair_contractions: expected.pair_contractions,
    };
    assert_eq!(
        prepared
            .integrate_interferometer_pixel(detector, footprint, exact, || false)
            .unwrap(),
        expected
    );
    for limits in [
        CoherentDiffractionExposurePixelLimits {
            spatial: CoherentDiffractionPixelLimits {
                maximum_samples: exact.spatial.maximum_samples - 1,
                ..exact.spatial
            },
            ..exact
        },
        CoherentDiffractionExposurePixelLimits {
            spatial: CoherentDiffractionPixelLimits {
                maximum_pupil_evaluations: exact.spatial.maximum_pupil_evaluations - 1,
                ..exact.spatial
            },
            ..exact
        },
        CoherentDiffractionExposurePixelLimits {
            maximum_pair_contractions: exact.maximum_pair_contractions - 1,
            ..exact
        },
    ] {
        let mut polls = 0;
        assert!(matches!(
            prepared.integrate_interferometer_pixel(detector, footprint, limits, || {
                polls += 1;
                false
            }),
            Err(CoherentDiffractionError::Limit { .. })
        ));
        assert_eq!(polls, 1);
    }
    let mut polls = 0;
    assert_eq!(
        prepared.integrate_interferometer_pixel(detector, footprint, exact, || {
            polls += 1;
            polls == 7
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        prepared
            .integrate_interferometer_pixel(detector, footprint, exact, || false)
            .unwrap(),
        expected
    );
    let one = prepared
        .integrate_interferometer_pixel(
            receiver(centers, 0.03, [1.0, 0.0], 0.0, 0.4),
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    let native = prepared
        .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    close(
        one.average_power.arm_power_w[0],
        native.average_power.power_w,
        3e-12,
    );
    close(
        one.average_power.output_power_w[0],
        0.5 * native.average_power.power_w,
        3e-12,
    );
    assert_eq!(
        one.average_power.output_power_w[0],
        one.average_power.output_power_w[1]
    );
    assert_eq!(
        one.average_power.pupil_evaluations,
        2 * one.average_power.quadrature_samples
    );
    assert_eq!(
        one.pair_contractions,
        3 * one.average_power.quadrature_samples
    );
    assert_eq!(one.average_power.mutual_coherence_w, C64::ZERO);
    let zero_limits = CoherentDiffractionExposurePixelLimits {
        spatial: CoherentDiffractionPixelLimits {
            maximum_samples: 0,
            maximum_pupil_evaluations: 0,
            ..PIXEL_LIMITS
        },
        maximum_pair_contractions: 0,
    };
    let closed = prepared
        .integrate_interferometer_pixel(
            receiver(centers, 0.03, [0.0; 2], 0.0, 0.0),
            footprint,
            zero_limits,
            || false,
        )
        .unwrap();
    assert_eq!(closed.average_power.subdivisions, [0; 2]);
    assert_eq!(closed.average_power.output_power_w, [0.0; 2]);
    assert_eq!(closed.average_power.pupil_evaluations, 0);
    assert_eq!(closed.pair_contractions, 0);
    let empty = camera
        .prepare_coherent_diffraction_exposure(&[], 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap()
        .integrate_interferometer_pixel(detector, footprint, zero_limits, || false)
        .unwrap();
    assert_eq!(empty.average_power.quadrature_samples, 0);
    assert_eq!(empty.output_energy_j, [0.0; 2]);
    let outside = receiver([[0.099, 0.0], [0.0; 2]], 0.03, [0.0, 1.0], 0.0, 0.0);
    let offsets = CoherentDiffractionPixel::try_new([0.0; 2], [0.002, 0.001]).unwrap();
    assert_eq!(
        prepared.integrate_interferometer_pixel(outside, offsets, EXPOSURE_PIXEL_LIMITS, || false),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
    assert_eq!(
        prepared.integrate_interferometer_pixel(
            receiver(centers, 0.04, [0.0; 2], 0.0, 0.0),
            footprint,
            zero_limits,
            || false
        ),
        Err(CoherentDiffractionError::InvalidSensorPlane)
    );
}
