//! Independent physical references for exact-frequency passive interferometer arms.
use super::*;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits, CoherentJonesArm,
    CoherentJonesInterferometerResponse, CoherentJonesPupilLayer, CoherentJonesPupilLimits,
    CoherentJonesPupilResponse, PolarizedCoherentInterferometerExposureSample,
    PolarizedDiffractionCoherency,
};

const C: f64 = 299_792_458.0;
const F0: f64 = 562_949_953_421_312.0;
const DELTA: f64 = 4_294_967_296.0;
const RECEIVER_LIMITS: CoherentInterferometerExposureLimits =
    CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 1000,
        maximum_pair_contractions: 10_000,
    };
const TWO_FREQUENCY_CACHE: usize = 2 * core::mem::size_of::<[CoherentJonesArm; 2]>();

fn receiver(phase: f64) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase).unwrap()
}

fn layer(power: [f64; 2], phase: [f64; 2], angle: f64) -> CoherentJonesPupilLayer {
    CoherentJonesPupilLayer::try_new(power, phase, angle, [0.0; 2]).unwrap()
}

fn arm(layers: &[CoherentJonesPupilLayer]) -> CoherentJonesArm {
    CoherentJonesArm::try_new(layers).unwrap()
}

fn response(frequency_hz: f64, arms: [CoherentJonesArm; 2]) -> CoherentJonesInterferometerResponse {
    CoherentJonesInterferometerResponse { frequency_hz, arms }
}

fn matrix(value: PolarizedDiffractionCoherency) -> [f64; 4] {
    let [x, y] = value.axis_powers();
    [x, y, value.cross().re, value.cross().im]
}

fn matrix_close(actual: PolarizedDiffractionCoherency, expected: [f64; 4], tolerance: f64) {
    let scale = (expected[0] + expected[1]).max(1e-280);
    for (a, e) in matrix(actual).into_iter().zip(expected) {
        assert!(
            (a - e).abs() <= tolerance * scale,
            "{actual:?} != {expected:?}"
        );
    }
}

fn physical_close(
    actual: PolarizedCoherentInterferometerExposureSample,
    expected: PolarizedCoherentInterferometerExposureSample,
) {
    for port in 0..2 {
        close(
            actual.arm_power_density_w[port],
            expected.arm_power_density_w[port],
            4e-12,
        );
        close(
            actual.arm_energy_density_j[port],
            expected.arm_energy_density_j[port],
            4e-12,
        );
        matrix_close(
            actual.output_coherency_density_w[port],
            matrix(expected.output_coherency_density_w[port]),
            4e-12,
        );
        matrix_close(
            actual.output_coherency_energy_density_j[port],
            matrix(expected.output_coherency_energy_density_j[port]),
            4e-12,
        );
    }
    let scale = expected.arm_power_density_w.iter().sum::<f64>().max(1e-280);
    assert!(
        (actual.mutual_coherence_density_w.re - expected.mutual_coherence_density_w.re).abs()
            < 4e-12 * scale
    );
    assert!(
        (actual.mutual_coherence_density_w.im - expected.mutual_coherence_density_w.im).abs()
            < 4e-12 * scale
    );
    close(actual.fringe_visibility, expected.fringe_visibility, 4e-12);
}

#[test]
fn g3_achromatic_arm_cache_matches_existing_point_and_pixel_with_patterned_pupils() {
    let camera = defocused(0.8);
    let a = [jones_wave(
        [0.1, -0.2],
        [C64::new(0.8, 0.2), C64::new(0.1, 0.5)],
    )];
    let b = [jones_wave(
        [-0.2, 0.1],
        [C64::new(0.3, -0.1), C64::new(-0.2, 0.4)],
    )];
    let tones = [jones_tone(F0, &a, &[]), jones_tone(F0 + DELTA, &b, &[])];
    let shutter = CoherentDiffractionShutter::try_new(-0.2 / DELTA, 0.3 / DELTA).unwrap();
    let mut source = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let layers =
        [CoherentJonesPupilLayer::try_new([0.7, 0.2], [0.1, 0.8], 0.3, [800.0, -350.0]).unwrap()];
    source
        .set_pupil_responses(
            &[
                CoherentJonesPupilResponse {
                    frequency_hz: F0,
                    layers: &layers,
                },
                CoherentJonesPupilResponse {
                    frequency_hz: F0 + DELTA,
                    layers: &layers,
                },
            ],
            CoherentJonesPupilLimits {
                maximum_layers: 2,
                maximum_modes: 12,
                maximum_pupil_evaluations: 6,
            },
            || false,
        )
        .unwrap();
    let arms = [
        arm(&[
            layer([0.8, 0.3], [0.2, -0.4], 0.3),
            layer([1.0; 2], [0.0, 0.7], -0.2),
        ]),
        arm(&[layer([0.4, 0.9], [0.1, 0.6], -0.4)]),
    ];
    let detector = CoherentDiffractionInterferometer::try_new(
        [[-0.3 * SCALE, 0.2 * SCALE], [0.5 * SCALE, -0.1 * SCALE]],
        0.03,
        [0.4, 0.6],
        0.12 * C / F0,
        -0.3,
    )
    .unwrap();
    let responses = [response(F0 + DELTA, arms), response(F0, arms)];
    let prepared = source
        .prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE, || false)
        .unwrap();
    assert_eq!(prepared.frequency_count(), 2);
    assert_eq!(prepared.cached_arm_bytes(), TWO_FREQUENCY_CACHE);
    assert_eq!(
        prepared.observe(RECEIVER_LIMITS, || false).unwrap(),
        source
            .observe_interferometer(detector, arms, RECEIVER_LIMITS, || false)
            .unwrap()
    );
    let footprint = pixel([-0.8, -0.6], [0.7, 0.9]);
    assert_eq!(
        prepared
            .integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || false)
            .unwrap(),
        source
            .integrate_interferometer_pixel(
                detector,
                arms,
                footprint,
                EXPOSURE_PIXEL_LIMITS,
                || false
            )
            .unwrap()
    );
}

#[derive(Clone, Copy)]
struct PrincipalLayer {
    power: [f64; 2],
    phase: [f64; 2],
    angle: f64,
}

impl PrincipalLayer {
    fn native(self) -> CoherentJonesArm {
        arm(&[layer(self.power, self.phase, self.angle)])
    }

    // Independent principal-axis projection and transmission, using ordinary
    // complex arithmetic and elementary trigonometry, not native arm matrices.
    fn apply_reference(self, field: [C64; 2]) -> [C64; 2] {
        let (s, c) = self.angle.sin_cos();
        let mut principal = [
            field[0].scale(c) + field[1].scale(s),
            field[1].scale(c) - field[0].scale(s),
        ];
        for axis in 0..2 {
            principal[axis] = principal[axis]
                * C64::new(self.phase[axis].cos(), self.phase[axis].sin())
                    .scale(self.power[axis].sqrt());
        }
        [
            principal[0].scale(c) - principal[1].scale(s),
            principal[0].scale(s) + principal[1].scale(c),
        ]
    }
}

// Exact two-line analytic shutter average. The one complex beat value is
// <exp(+i 2*pi*delta*(d/c-t))>; diagonal terms have unit temporal weight.
fn covariance(left: [C64; 2], right: [C64; 2], beat: C64) -> C64 {
    left[0] * right[0].conj()
        + left[1] * right[1].conj()
        + left[1] * right[0].conj() * beat
        + left[0] * right[1].conj() * beat.conj()
}

fn reference_coherency(fields: [[C64; 2]; 2], beat: C64) -> [f64; 4] {
    let x = fields.map(|field| field[0]);
    let y = fields.map(|field| field[1]);
    let cross = covariance(x, y, beat);
    [
        covariance(x, x, beat).re,
        covariance(y, y, beat).re,
        cross.re,
        cross.im,
    ]
}

#[test]
fn g1_chromatic_transmission_retardance_and_axes_match_independent_beat_reference() {
    let camera = focused();
    let frequencies = [F0, F0 + DELTA];
    let amplitudes = [
        [C64::ONE, C64::new(0.25, 0.4)],
        [C64::new(0.45, -0.2), C64::new(0.7, 0.1)],
    ];
    let waves = amplitudes.map(|amplitude| [jones_wave([0.0; 2], amplitude)]);
    let tones = [
        jones_tone(frequencies[0], &waves[0], &[]),
        jones_tone(frequencies[1], &waves[1], &[]),
    ];
    let optics = [
        [
            PrincipalLayer {
                power: [0.7, 0.2],
                phase: [0.2, -0.5],
                angle: 0.27,
            },
            PrincipalLayer {
                power: [0.3, 0.8],
                phase: [-0.1, 0.6],
                angle: -0.38,
            },
        ],
        [
            PrincipalLayer {
                power: [0.1, 0.9],
                phase: [0.8, 0.1],
                angle: -0.44,
            },
            PrincipalLayer {
                power: [0.8, 0.15],
                phase: [0.4, -0.9],
                angle: 0.63,
            },
        ],
    ];
    let responses = [
        response(frequencies[1], optics[1].map(PrincipalLayer::native)),
        response(frequencies[0], optics[0].map(PrincipalLayer::native)),
    ];
    let fractions = [0.35_f64, 0.55];
    let path = 0.19 * C / F0;
    let phase = 0.31;
    let distance = 0.03;
    let detector =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], distance, fractions, path, phase)
            .unwrap();
    for width in [0.15_f64, 1.0] {
        let duration = width / DELTA;
        let start = -0.13 / DELTA;
        let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
        let source = camera
            .prepare_coherent_jones_diffraction_exposure(
                &tones,
                distance,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let prepared = source
            .prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE, || false)
            .unwrap();
        let actual = prepared.observe(RECEIVER_LIMITS, || false).unwrap();
        let sinc = if width == 1.0 {
            0.0
        } else {
            sinc_reference(core::f64::consts::PI * width)
        };
        let angle = core::f64::consts::TAU * DELTA * (distance / C - start - 0.5 * duration);
        let beat = C64::new(angle.cos(), angle.sin()).scale(sinc);
        let fields: [[[C64; 2]; 2]; 2] = std::array::from_fn(|arm| {
            std::array::from_fn(|frequency| {
                let peak = core::f64::consts::PI.sqrt() * RADIUS * frequencies[frequency] / C;
                let input =
                    amplitudes[frequency].map(|value| value.scale(peak * fractions[arm].sqrt()));
                let shifted = if arm == 0 {
                    0.0
                } else {
                    core::f64::consts::TAU * frequencies[frequency] * path / C + phase
                };
                optics[frequency][arm]
                    .apply_reference(input)
                    .map(|value| value * C64::new(shifted.cos(), shifted.sin()))
            })
        });
        let arms = fields.map(|field| {
            let coherency = reference_coherency(field, beat);
            coherency[0] + coherency[1]
        });
        let mutual = covariance(
            fields[0].map(|field| field[0]),
            fields[1].map(|field| field[0]),
            beat,
        ) + covariance(
            fields[0].map(|field| field[1]),
            fields[1].map(|field| field[1]),
            beat,
        );
        for (port, sign) in [1.0, -1.0].into_iter().enumerate() {
            let output = std::array::from_fn(|frequency| {
                std::array::from_fn(|axis| {
                    (fields[0][frequency][axis] + fields[1][frequency][axis].scale(sign))
                        .scale(core::f64::consts::FRAC_1_SQRT_2)
                })
            });
            let expected = reference_coherency(output, beat);
            matrix_close(actual.output_coherency_density_w[port], expected, 4e-12);
            matrix_close(
                actual.output_coherency_energy_density_j[port],
                expected.map(|value| value * duration),
                4e-12,
            );
            close(actual.arm_power_density_w[port], arms[port], 4e-12);
            close(
                actual.arm_energy_density_j[port],
                arms[port] * duration,
                4e-12,
            );
        }
        close(actual.mutual_coherence_density_w.re, mutual.re, 4e-12);
        close(actual.mutual_coherence_density_w.im, mutual.im, 4e-12);
        close(
            actual.fringe_visibility,
            2.0 * mutual.abs() / (arms[0] + arms[1]),
            4e-12,
        );
        close(
            actual
                .output_coherency_density_w
                .iter()
                .map(|value| value.stokes().i)
                .sum(),
            arms.iter().sum(),
            4e-12,
        );
        assert_eq!(
            (actual.pupil_evaluations, actual.pair_contractions),
            (4, 52)
        );
        if width == 1.0 {
            // Over a complete beat period, passive spectral loss lowers the
            // collected power without the finite-window redistribution caveat.
            let clear = frequencies
                .into_iter()
                .zip(amplitudes)
                .map(|(frequency, amplitude)| {
                    let peak = core::f64::consts::PI.sqrt() * RADIUS * frequency / C;
                    peak * peak * (amplitude[0].norm_sq() + amplitude[1].norm_sq())
                })
                .sum::<f64>()
                * (fractions[0] + fractions[1]);
            assert!(arms.iter().sum::<f64>() < 0.8 * clear);
        }
    }
}

#[test]
fn g1_common_spectral_arm_phase_changes_coherent_beats_but_not_independent_lines() {
    let camera = focused();
    let first = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let second = [jones_wave(
        [0.0; 2],
        [C64::from_re(F0 / (F0 + DELTA)), C64::ZERO],
    )];
    let coherent = [
        jones_tone(F0, &first, &[]),
        jones_tone(F0 + DELTA, &second, &[]),
    ];
    let factors0 = [
        JonesPupilCoherenceFactor { waves: &first },
        JonesPupilCoherenceFactor { waves: &[] },
    ];
    let factors1 = [
        JonesPupilCoherenceFactor { waves: &[] },
        JonesPupilCoherenceFactor { waves: &second },
    ];
    let independent = [
        jones_tone(F0, &[], &factors0),
        jones_tone(F0 + DELTA, &[], &factors1),
    ];
    let duration = 0.25 / DELTA;
    let shutter = CoherentDiffractionShutter::try_new(0.03 / C - 0.5 * duration, duration).unwrap();
    let clear = CoherentJonesArm::IDENTITY;
    let inverted = arm(&[layer([1.0; 2], [core::f64::consts::PI; 2], 0.0)]);
    let peak = core::f64::consts::PI * RADIUS * RADIUS * (F0 / C).powi(2);
    let temporal = sinc_reference(core::f64::consts::FRAC_PI_4);
    for (tones, correlated) in [(&coherent, true), (&independent, false)] {
        let source = camera
            .prepare_coherent_jones_diffraction_exposure(
                tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        for (high_frequency_arm, sign) in [(clear, 1.0), (inverted, -1.0)] {
            let responses = [
                response(F0, [clear; 2]),
                response(F0 + DELTA, [high_frequency_arm; 2]),
            ];
            let prepared = source
                .prepare_interferometer(receiver(0.0), &responses, TWO_FREQUENCY_CACHE, || false)
                .unwrap();
            let actual = prepared.observe(RECEIVER_LIMITS, || false).unwrap();
            let expected = 2.0 * peak * (1.0 + if correlated { sign * temporal } else { 0.0 });
            matrix_close(
                actual.output_coherency_density_w[0],
                [expected, 0.0, 0.0, 0.0],
                4e-12,
            );
            assert_eq!(actual.output_coherency_density_w[1].stokes().i, 0.0);
            close(actual.arm_power_density_w[0], 0.5 * expected, 4e-12);
            close(actual.arm_power_density_w[1], 0.5 * expected, 4e-12);
            close(actual.fringe_visibility, 1.0, 4e-12);
            let collected = prepared
                .integrate_pixel(pixel([-0.2; 2], [0.2; 2]), EXPOSURE_PIXEL_LIMITS, || false)
                .unwrap();
            assert_eq!(collected.output_coherency_w[1].stokes().i, 0.0);
            close(
                collected.output_coherency_w[0].stokes().i,
                collected.arm_power_w.iter().sum(),
                4e-12,
            );
        }
    }
}

#[test]
fn g3_response_order_and_equal_tone_grouping_preserve_one_physical_arm_spectrum() {
    let camera = focused();
    let positive = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 0.4)])];
    let negative = [jones_wave(
        [0.0; 2],
        [C64::new(-1.0, 0.0), C64::new(0.0, -0.4)],
    )];
    let remaining = [jones_wave(
        [0.1, -0.2],
        [C64::new(0.3, -0.2), C64::new(0.8, 0.1)],
    )];
    let tones = [
        jones_tone(F0, &positive, &[]),
        jones_tone(F0 + DELTA, &remaining, &[]),
        jones_tone(F0, &negative, &[]),
    ];
    let reordered = [tones[1], tones[2], tones[0]];
    let surviving = [tones[1]];
    let shutter = CoherentDiffractionShutter::try_new(-0.3 / DELTA, 0.2 / DELTA).unwrap();
    let arms0 = [
        arm(&[layer([0.2, 0.7], [0.9, -0.5], 0.4)]),
        CoherentJonesArm::IDENTITY,
    ];
    let arms1 = [
        arm(&[layer([0.8, 0.1], [-0.2, 0.7], -0.3)]),
        arm(&[layer([0.5, 0.9], [0.4, -0.6], 0.2)]),
    ];
    let responses = [response(F0, arms0), response(F0 + DELTA, arms1)];
    let reversed_responses = [response(F0 + DELTA, arms1), response(F0, arms0)];
    let detector = receiver(-0.4);
    let source = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let prepared = source
        .prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE, || false)
        .unwrap();
    assert_eq!(prepared.frequency_count(), 2);
    assert_eq!(prepared.cached_arm_bytes(), TWO_FREQUENCY_CACHE);
    let expected = prepared.observe(RECEIVER_LIMITS, || false).unwrap();
    let reversed = source
        .prepare_interferometer(detector, &reversed_responses, TWO_FREQUENCY_CACHE, || false)
        .unwrap();
    assert_eq!(
        reversed.observe(RECEIVER_LIMITS, || false).unwrap(),
        expected
    );
    for source_tones in [&reordered[..], &surviving[..]] {
        let source = camera
            .prepare_coherent_jones_diffraction_exposure(
                source_tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let entries = if source_tones.len() == 1 {
            &responses[1..]
        } else {
            &responses[..]
        };
        let prepared = source
            .prepare_interferometer(detector, entries, TWO_FREQUENCY_CACHE, || false)
            .unwrap();
        physical_close(
            prepared.observe(RECEIVER_LIMITS, || false).unwrap(),
            expected,
        );
    }
}

#[test]
fn g0_g4_exact_arm_cache_admission_and_cancelled_reuse() {
    let camera = focused();
    let waves = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.2, 0.3)])];
    let tones = [
        jones_tone(F0, &waves, &[]),
        jones_tone(F0 + DELTA, &waves, &[]),
    ];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.125 / DELTA).unwrap();
    let source = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let arms = [CoherentJonesArm::IDENTITY; 2];
    let responses = [response(F0, arms), response(F0 + DELTA, arms)];
    let detector = receiver(0.3);
    let mut polls = 0;
    assert!(matches!(
        source.prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE - 1, || { polls += 1; false }),
        Err(CoherentDiffractionError::Limit { quantity: "Jones interferometer cached arm bytes", requested, maximum })
            if requested == TWO_FREQUENCY_CACHE && maximum == TWO_FREQUENCY_CACHE - 1
    ));
    assert_eq!(polls, 1);
    for invalid in [
        vec![],
        vec![response(F0, arms)],
        vec![response(F0, arms), response(F0, arms)],
        vec![response(F0, arms), response(F0 + 1.0, arms)],
        vec![
            response(F0, arms),
            response(F0 + DELTA, arms),
            response(F0 + 1.0, arms),
        ],
    ] {
        assert!(matches!(
            source.prepare_interferometer(detector, &invalid, usize::MAX, || false),
            Err(CoherentDiffractionError::InconsistentJonesArmSpectrum)
        ));
    }
    for frequency in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let invalid = [response(F0, arms), response(frequency, arms)];
        assert!(matches!(
            source.prepare_interferometer(detector, &invalid, usize::MAX, || false),
            Err(CoherentDiffractionError::InvalidFrequency)
        ));
    }
    let wrong_sensor =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.04, [0.5; 2], 0.0, 0.0)
            .unwrap();
    assert!(matches!(
        source.prepare_interferometer(wrong_sensor, &responses, TWO_FREQUENCY_CACHE, || false),
        Err(CoherentDiffractionError::InvalidSensorPlane)
    ));
    let mut polls = 0;
    assert!(matches!(
        source.prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE, || {
            polls += 1;
            polls == 2
        }),
        Err(CoherentDiffractionError::Cancelled)
    ));
    let prepared = source
        .prepare_interferometer(detector, &responses, TWO_FREQUENCY_CACHE, || false)
        .unwrap();
    let expected = prepared.observe(RECEIVER_LIMITS, || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        prepared.observe(RECEIVER_LIMITS, || {
            polls += 1;
            polls == 5
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    let footprint = pixel([-0.3; 2], [0.4; 2]);
    let mut polls = 0;
    assert_eq!(
        prepared.integrate_pixel(footprint, EXPOSURE_PIXEL_LIMITS, || {
            polls += 1;
            polls == 7
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        prepared.observe(RECEIVER_LIMITS, || false).unwrap(),
        expected
    );
    assert_eq!(prepared.cached_arm_bytes(), TWO_FREQUENCY_CACHE);
    let empty_source = camera
        .prepare_coherent_jones_diffraction_exposure(&[], 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let empty = empty_source
        .prepare_interferometer(detector, &[], 0, || false)
        .unwrap();
    assert_eq!(empty.frequency_count(), 0);
    assert_eq!(empty.cached_arm_bytes(), 0);
    let zero = empty
        .observe(
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: 0,
                maximum_pair_contractions: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(zero.arm_power_density_w, [0.0; 2]);
    assert_eq!((zero.pupil_evaluations, zero.pair_contractions), (0, 0));
    let zero_pixel = empty
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
            || false,
        )
        .unwrap();
    assert_eq!(zero_pixel.quadrature_samples, 0);
    assert_eq!(zero_pixel.output_coherency_w[0].stokes().i, 0.0);
}
