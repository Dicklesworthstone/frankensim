use super::*;
use fs_math::dd::Dd;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits,
    CoherentInterferometerExposureSample,
};

const C: f64 = 299_792_458.0;
const F0: f64 = 562_949_953_421_312.0; // 2^49 Hz; visible vacuum wavelength.
const RECEIVER_LIMITS: CoherentInterferometerExposureLimits =
    CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 128,
        maximum_pair_contractions: 1000,
    };

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

fn center_field(camera: &PhysicalCamera, frequency: f64, waves: &[CoherentPupilWave]) -> C64 {
    camera
        .prepare_coherent_diffraction(
            reference_value(Dd::from_f64(C) / Dd::from_f64(frequency) * Dd::from_f64(1e9)),
            CoherentPupilIllumination {
                mean: waves,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap()
        .observe([0.0; 2], || false)
        .unwrap()
        .mean_field
}

// Independent two-line analytic-field integral after removing one common
// carrier. Only the low difference frequency appears in this reference.
// Retain the reference field algebra in DD: a half-cycle OPD makes Re(W01)
// the small difference of bright powers, where a binary64 oracle loses the
// residual. This reference never calls the production temporal kernel.
type ReferenceField = [Dd; 2];

fn reference_value(value: Dd) -> f64 {
    value.hi + value.lo
}

fn reference_field(value: C64) -> ReferenceField {
    [Dd::from_f64(value.re), Dd::from_f64(value.im)]
}

fn reference_mul(a: ReferenceField, b: ReferenceField) -> ReferenceField {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn reference_split(value: C64) -> ReferenceField {
    reference_field(value).map(|part| part * Dd::from_f64(0.5).sqrt())
}

fn covariance(a: [ReferenceField; 2], b: [ReferenceField; 2], beat: C64) -> ReferenceField {
    let mut result = [Dd::ZERO; 2];
    for (i, j, window) in [
        (0, 0, C64::ONE),
        (1, 1, C64::ONE),
        (1, 0, beat),
        (0, 1, beat.conj()),
    ] {
        let term = reference_mul(
            reference_mul(a[i], [b[j][0], -b[j][1]]),
            reference_field(window),
        );
        for axis in 0..2 {
            result[axis] = result[axis] + term[axis];
        }
    }
    result
}

fn beat_mean(delta: f64, start: f64, duration: f64, propagation_time: f64) -> C64 {
    let cycles = delta * duration;
    // The exact integral over any complete beat period is zero. Evaluating
    // binary64 sin(pi*cycles) would create a false ~epsilon bright-power term.
    if cycles != 0.0 && cycles.fract() == 0.0 {
        return C64::ZERO;
    }
    let half = core::f64::consts::PI * cycles;
    let phase = -core::f64::consts::TAU * delta * (start + 0.5 * duration - propagation_time);
    C64::new(phase.cos(), phase.sin()).scale(sinc_reference(half))
}

fn check_reference(
    sample: CoherentInterferometerExposureSample,
    a: [ReferenceField; 2],
    b: [ReferenceField; 2],
    beat: C64,
    duration: f64,
) {
    let arms = [covariance(a, a, beat)[0], covariance(b, b, beat)[0]];
    let cross = covariance(a, b, beat);
    let split = Dd::from_f64(0.5).sqrt();
    let plus =
        std::array::from_fn(|i| std::array::from_fn(|axis| (a[i][axis] + b[i][axis]) * split));
    let minus =
        std::array::from_fn(|i| std::array::from_fn(|axis| (a[i][axis] - b[i][axis]) * split));
    let ports = [
        covariance(plus, plus, beat)[0],
        covariance(minus, minus, beat)[0],
    ];
    for i in 0..2 {
        close(
            sample.arm_power_density_w[i],
            reference_value(arms[i]),
            3e-12,
        );
        close(
            sample.output_power_density_w[i],
            reference_value(ports[i]),
            3e-12,
        );
        close(
            sample.arm_energy_density_j[i],
            reference_value(arms[i] * Dd::from_f64(duration)),
            3e-12,
        );
        close(
            sample.output_energy_density_j[i],
            reference_value(ports[i] * Dd::from_f64(duration)),
            3e-12,
        );
    }
    // The native DD trigonometric contraction can leave epsilon^2-scaled
    // residuals at an exact integral null. Keep the relative comparison and
    // bound only this absolute roundoff; never use a binary64 bright-power floor.
    let floor = 64.0 * f64::EPSILON * f64::EPSILON * reference_value(arms[0] + arms[1]);
    for (actual, expected) in [
        (
            sample.mutual_coherence_density_w.re,
            reference_value(cross[0]),
        ),
        (
            sample.mutual_coherence_density_w.im,
            reference_value(cross[1]),
        ),
    ] {
        assert!(
            (actual - expected).abs() <= 3e-12 * expected.abs() + floor,
            "mutual coherence actual {actual:.16e}, expected {expected:.16e}, floor {floor:.3e}",
        );
    }
    close(
        sample.fringe_visibility,
        reference_value(
            (cross[0] * cross[0] + cross[1] * cross[1]).sqrt() * Dd::from_f64(2.0)
                / (arms[0] + arms[1]),
        ),
        3e-12,
    );
    close(
        sample.output_power_density_w.iter().sum(),
        sample.arm_power_density_w.iter().sum(),
        3e-12,
    );
}

#[test]
fn g1_monochromatic_exposure_matches_the_original_complex_two_arm_receiver() {
    let camera = defocused(1.4);
    let mean = [wave([0.1, -0.2], C64::new(0.8, -0.3))];
    let random = [wave([-0.4, 0.1], C64::new(-0.2, 0.5))];
    let factors = [PupilCoherenceFactor { waves: &random }];
    let tones = [tone(F0, &mean, &factors)];
    let scalar = camera
        .prepare_coherent_diffraction(
            C / F0 * 1e9,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &factors,
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let coordinates = [[0.1 * SCALE, 0.0], [0.8 * SCALE, 0.2 * SCALE]];
    for duration in [1e-17, 0.04] {
        let query = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                0.03,
                CoherentDiffractionShutter::try_new(2.0_f64.powi(40), duration).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        for path in [0.0, C / (4.0 * F0)] {
            let receiver = receiver(coordinates, 0.03, [0.25, 0.75], path, 0.27);
            let original = scalar
                .observe_interferometer(receiver, 4, || false)
                .unwrap();
            let sample = query
                .observe_interferometer(receiver, RECEIVER_LIMITS, || false)
                .unwrap();
            for i in 0..2 {
                close(
                    sample.arm_power_density_w[i],
                    original.arm_power_density_w[i],
                    3e-12,
                );
                close(
                    sample.output_power_density_w[i],
                    original.output_power_density_w[i],
                    3e-12,
                );
                close(
                    sample.output_energy_density_j[i],
                    original.output_power_density_w[i] * duration,
                    3e-12,
                );
            }
            close(
                sample.mutual_coherence_density_w.re,
                original.mutual_coherence_density_w.re,
                3e-12,
            );
            close(
                sample.mutual_coherence_density_w.im,
                original.mutual_coherence_density_w.im,
                3e-12,
            );
            close(sample.fringe_visibility, original.fringe_visibility, 3e-12);
            assert_eq!(sample.pupil_evaluations, 4);
            assert_eq!(sample.pair_contractions, 10);
        }
    }
}

#[test]
fn g1_two_tone_beat_washout_and_complex_mutual_coherence_match_analytic_fields() {
    let camera = focused();
    let wave_a = [wave([0.0; 2], C64::ONE)];
    let wave_b = [wave([0.0; 2], C64::ONE)];
    let delta = 2.0_f64.powi(30);
    let tones = [tone(F0, &wave_a, &[]), tone(F0 + delta, &wave_b, &[])];
    let fields = [
        center_field(&camera, F0, &wave_a),
        center_field(&camera, F0 + delta, &wave_b),
    ];
    let a = fields.map(reference_split);
    let distance = C * 2.0_f64.powi(-33); // delta*d/c = 1/8 exactly.
    let path = C * 2.0_f64.powi(-31); // Carrier: whole cycles; second line: half a cycle.
    let start = 0.25 / delta;
    let mut visibilities = [0.0; 2];
    for (index, duration) in [0.125 / delta, 4.0 / delta].into_iter().enumerate() {
        let query = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                distance,
                CoherentDiffractionShutter::try_new(start, duration).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let beat = beat_mean(delta, start, duration, 0.125 / delta);
        for phase in [0.0_f64, 0.37] {
            let rotation = C64::new(phase.cos(), phase.sin());
            let b = [
                reference_mul(a[0], reference_field(rotation)),
                reference_mul(a[1], reference_field(rotation.scale(-1.0))),
            ];
            let sample = query
                .observe_interferometer(
                    receiver([[0.0; 2]; 2], distance, [0.5; 2], path, phase),
                    RECEIVER_LIMITS,
                    || false,
                )
                .unwrap();
            check_reference(sample, a, b, beat, duration);
            visibilities[index] = sample.fringe_visibility;
            assert_eq!(sample.pair_contractions, 18);
        }
    }
    assert!(visibilities[0] > 0.8);
    assert!(visibilities[1] < 1e-5);
}

#[test]
fn g3_shared_spectral_factors_retain_beats_but_independent_sources_do_not() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let shared = [PupilCoherenceFactor { waves: &waves }];
    let independent_a = [
        PupilCoherenceFactor { waves: &waves },
        PupilCoherenceFactor { waves: &[] },
    ];
    let independent_b = [
        PupilCoherenceFactor { waves: &[] },
        PupilCoherenceFactor { waves: &waves },
    ];
    let deterministic = [tone(F0, &waves, &[]), tone(F0 + 4.0, &waves, &[])];
    let correlated = [tone(F0, &[], &shared), tone(F0 + 4.0, &[], &shared)];
    let independent = [
        tone(F0, &[], &independent_a),
        tone(F0 + 4.0, &[], &independent_b),
    ];
    let mut powers = [0.0; 3];
    for (index, tones) in [&deterministic, &correlated, &independent]
        .into_iter()
        .enumerate()
    {
        let query = camera
            .prepare_coherent_diffraction_exposure(
                tones,
                0.03,
                CoherentDiffractionShutter::try_new(0.0, 0.0625).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let sample = query
            .observe_interferometer(
                receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, 0.0),
                RECEIVER_LIMITS,
                || false,
            )
            .unwrap();
        powers[index] = sample.output_power_density_w[0];
        assert_eq!(sample.output_power_density_w[1], 0.0);
        close(sample.fringe_visibility, 1.0, 1e-14);
    }
    assert_eq!(powers[0], powers[1]);
    let independent_power = center_field(&camera, F0, &waves).norm_sq()
        + center_field(&camera, F0 + 4.0, &waves).norm_sq();
    close(powers[2], independent_power, 3e-12);
    assert!(powers[0] > 1.5 * powers[2]);
}

#[test]
fn g3_equal_frequency_coalescing_and_multitone_dark_ports_keep_weak_fields() {
    let camera = focused();
    let strong = [wave([0.0; 2], C64::new(1e16, 0.0))];
    let weak = [wave([0.0; 2], C64::new(0.25, -0.125))];
    let negative = [wave([0.0; 2], C64::new(-1e16, 0.0))];
    let second = [wave([0.0; 2], C64::new(0.4, 0.2))];
    let expanded = [
        tone(F0, &strong, &[]),
        tone(F0, &weak, &[]),
        tone(F0, &negative, &[]),
        tone(F0 + 4.0, &second, &[]),
    ];
    let reduced = [tone(F0, &weak, &[]), tone(F0 + 4.0, &second, &[])];
    let shutter = CoherentDiffractionShutter::try_new(2.0_f64.powi(40), 0.05).unwrap();
    let expanded = camera
        .prepare_coherent_diffraction_exposure(&expanded, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    let reduced = camera
        .prepare_coherent_diffraction_exposure(&reduced, 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    assert_eq!(expanded.frequency_count(), 2);
    let source_power = reduced.observe([0.0; 2], || false).unwrap().power_density_w;
    for phase in [0.0, 1e-20] {
        let receiver = receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase);
        let a = expanded
            .observe_interferometer(receiver, RECEIVER_LIMITS, || false)
            .unwrap();
        let b = reduced
            .observe_interferometer(receiver, RECEIVER_LIMITS, || false)
            .unwrap();
        for i in 0..2 {
            close(
                a.output_power_density_w[i],
                b.output_power_density_w[i],
                3e-12,
            );
        }
        assert_eq!(a.pair_contractions, b.pair_contractions);
        if phase == 0.0 {
            assert_eq!(a.output_power_density_w[1], 0.0);
            close(a.output_power_density_w[0], source_power, 3e-12);
        } else {
            assert!(a.output_power_density_w[1] > 0.0);
            close(
                a.output_power_density_w[1] / (0.25 * source_power * phase * phase),
                1.0,
                3e-12,
            );
        }
    }
}

#[test]
fn g1_long_optical_paths_keep_authored_close_frequency_beats_on_the_source_clock() {
    let camera = focused();
    let first = [wave([0.0; 2], C64::ONE)];
    let second = [wave([0.0; 2], C64::new(0.4, 0.3))];
    let tones = [tone(F0, &first, &[]), tone(F0 + 1.0, &second, &[])];
    let a = [
        center_field(&camera, F0, &first),
        center_field(&camera, F0 + 1.0, &second),
    ]
    .map(reference_split);
    let distance = C / 16.0;
    let duration = 0.125;
    for clock_sign in [-1.0, 1.0] {
        let start = clock_sign * 2.0_f64.powi(40) + 0.125;
        let query = camera
            .prepare_coherent_diffraction_exposure(
                &tones,
                distance,
                CoherentDiffractionShutter::try_new(start, duration).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        // Both the carrier and all whole-second propagation delays contribute
        // exact complete cycles. The second tone alone retains an eighth turn.
        let beat = beat_mean(1.0, start.fract(), duration, 1.0 / 16.0);
        for path_sign in [-1.0, 1.0] {
            let path = path_sign * C * (2.0_f64.powi(20) + 0.125);
            let phase = path_sign * core::f64::consts::FRAC_PI_4;
            let b = [
                a[0],
                reference_mul(a[1], reference_field(C64::new(phase.cos(), phase.sin()))),
            ];
            let sample = query
                .observe_interferometer(
                    receiver([[0.0; 2]; 2], distance, [0.5; 2], path, 0.0),
                    RECEIVER_LIMITS,
                    || false,
                )
                .unwrap();
            check_reference(sample, a, b, beat, duration);
        }
        let weak_phase = 1e-20;
        for path_sign in [-1.0, 1.0] {
            let sample = query
                .observe_interferometer(
                    receiver(
                        [[0.0; 2]; 2],
                        distance,
                        [0.5; 2],
                        path_sign * C * 2.0_f64.powi(40),
                        weak_phase,
                    ),
                    RECEIVER_LIMITS,
                    || false,
                )
                .unwrap();
            let source = query.observe([0.0; 2], || false).unwrap().power_density_w;
            assert!(sample.output_power_density_w[1] > 0.0);
            close(
                sample.output_power_density_w[1] / (source * 0.25 * weak_phase * weak_phase),
                1.0,
                3e-12,
            );
        }
    }
}

#[test]
fn g0_g4_spectral_receiver_preflights_work_and_reuses_after_cancellation() {
    let camera = focused();
    let waves = [wave([0.0; 2], C64::ONE)];
    let tones = [tone(F0, &waves, &[]), tone(F0 + 4.0, &waves, &[])];
    let query = camera
        .prepare_coherent_diffraction_exposure(
            &tones,
            0.03,
            CoherentDiffractionShutter::try_new(0.0, 0.1).unwrap(),
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let both = receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, 0.0);
    for limits in [
        CoherentInterferometerExposureLimits {
            maximum_pupil_evaluations: 3,
            maximum_pair_contractions: 18,
        },
        CoherentInterferometerExposureLimits {
            maximum_pupil_evaluations: 4,
            maximum_pair_contractions: 17,
        },
    ] {
        let mut polls = 0;
        assert!(matches!(
            query.observe_interferometer(both, limits, || {
                polls += 1;
                false
            }),
            Err(CoherentDiffractionError::Limit { .. })
        ));
        assert_eq!(polls, 1);
    }
    let exact = CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 4,
        maximum_pair_contractions: 18,
    };
    let original = query.observe_interferometer(both, exact, || false).unwrap();
    assert_eq!(original.pupil_evaluations, 4);
    assert_eq!(original.pair_contractions, 18);
    let mut polls = 0;
    assert_eq!(
        query.observe_interferometer(both, exact, || {
            polls += 1;
            polls == 3
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe_interferometer(both, exact, || true),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe_interferometer(both, exact, || false).unwrap(),
        original
    );
    assert_eq!(
        query.observe_interferometer(
            receiver([[0.0; 2]; 2], 0.031, [0.5; 2], 0.0, 0.0),
            exact,
            || false
        ),
        Err(CoherentDiffractionError::InvalidSensorPlane),
    );
    let one = query
        .observe_interferometer(
            receiver([[0.0; 2]; 2], 0.03, [1.0, 0.0], 0.0, 0.0),
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: 2,
                maximum_pair_contractions: 3,
            },
            || false,
        )
        .unwrap();
    assert_eq!(one.pupil_evaluations, 2);
    assert_eq!(one.pair_contractions, 3);
    assert_eq!(one.output_power_density_w[0], one.output_power_density_w[1]);
    assert_eq!(one.mutual_coherence_density_w, C64::ZERO);
    assert_eq!(one.fringe_visibility, 0.0);
    let closed = query
        .observe_interferometer(
            receiver([[0.0; 2]; 2], 0.03, [0.0; 2], f64::MAX, 0.0),
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: 0,
                maximum_pair_contractions: 0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(closed.output_power_density_w, [0.0; 2]);
    assert_eq!(closed.output_energy_density_j, [0.0; 2]);
    assert_eq!(closed.pupil_evaluations, 0);
    assert_eq!(closed.pair_contractions, 0);
}
