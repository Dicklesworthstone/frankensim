//! Original-source polarization, coherence and shutter laws at both ports.
use super::*;
use fs_math::dd::Dd;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits, CoherentJonesArm,
    CoherentJonesPupilLayer, CoherentJonesPupilLimits, CoherentJonesPupilResponse,
    PolarizedCoherentInterferometerExposureSample, PolarizedDiffractionCoherency,
};
use fs_render::polarization::Mueller;

const C: f64 = 299_792_458.0;
const F0: f64 = 562_949_953_421_312.0;
const RECEIVER_LIMITS: CoherentInterferometerExposureLimits =
    CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 1000,
        maximum_pair_contractions: 10_000,
    };

fn wavelength(frequency: f64) -> f64 {
    let value = Dd::from_f64(C) / Dd::from_f64(frequency) * Dd::from_f64(1e9);
    value.hi + value.lo
}

fn receiver(phase: f64) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase).unwrap()
}

fn layer(power: [f64; 2], phase: [f64; 2], angle: f64) -> CoherentJonesPupilLayer {
    CoherentJonesPupilLayer::try_new(power, phase, angle, [0.0; 2]).unwrap()
}

fn arm(layers: &[CoherentJonesPupilLayer]) -> CoherentJonesArm {
    CoherentJonesArm::try_new(layers).unwrap()
}

fn matrix(coherency: PolarizedDiffractionCoherency) -> [f64; 4] {
    let [x, y] = coherency.axis_powers();
    [x, y, coherency.cross().re, coherency.cross().im]
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

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual > 0.0 && expected > 0.0, "{actual} != {expected}");
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn closure(sample: PolarizedCoherentInterferometerExposureSample) {
    let inputs = sample.arm_power_density_w;
    let outputs = sample.output_coherency_density_w.map(|c| c.stokes().i);
    let scale = inputs.into_iter().chain(outputs).fold(0.0_f64, f64::max);
    if scale > 0.0 {
        let residual = inputs.iter().map(|p| p / scale).sum::<f64>()
            - outputs.iter().map(|p| p / scale).sum::<f64>();
        assert!(residual.abs() < 4e-12, "{sample:?}");
    }
}

#[test]
fn g1_monochromatic_receiver_matches_full_native_polarized_ports() {
    let camera = defocused(1.4);
    let mean = [jones_wave(
        [0.1, -0.2],
        [C64::new(0.8, -0.3), C64::new(0.2, 0.5)],
    )];
    let random = [jones_wave(
        [-0.4, 0.1],
        [C64::new(-0.2, 0.5), C64::new(0.7, -0.1)],
    )];
    let factors = [JonesPupilCoherenceFactor { waves: &random }];
    let tones = [jones_tone(F0, &mean, &factors)];
    let mono = camera
        .prepare_polarized_coherent_diffraction(
            wavelength(F0),
            CoherentJonesPupilIllumination {
                mean: &mean,
                factors: &factors,
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let arms = [
        arm(&[layer([0.7, 0.3], [0.2, -0.3], 0.23)]),
        arm(&[layer([0.6, 1.0], [-0.1, 0.7], -0.35)]),
    ];
    let coordinates = [[0.1 * SCALE, 0.0], [0.8 * SCALE, 0.2 * SCALE]];
    for duration in [1e-17, 0.07] {
        let shutter = CoherentDiffractionShutter::try_new(2.0_f64.powi(40), duration).unwrap();
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        for path in [0.0, C / (4.0 * F0)] {
            let receiver = CoherentDiffractionInterferometer::try_new(
                coordinates,
                0.03,
                [0.25, 0.75],
                path,
                0.27,
            )
            .unwrap();
            let expected = mono
                .observe_interferometer(receiver, arms, 4, || false)
                .unwrap();
            let sample = query
                .observe_interferometer(receiver, arms, RECEIVER_LIMITS, || false)
                .unwrap();
            for port in 0..2 {
                close(
                    sample.arm_power_density_w[port],
                    expected.arm_power_density_w[port],
                    4e-12,
                );
                matrix_close(
                    sample.output_coherency_density_w[port],
                    matrix(expected.output_coherency_density_w[port]),
                    4e-12,
                );
                matrix_close(
                    sample.output_coherency_energy_density_j[port],
                    matrix(expected.output_coherency_density_w[port]).map(|v| v * duration),
                    4e-12,
                );
                close(
                    sample.arm_energy_density_j[port],
                    expected.arm_power_density_w[port] * duration,
                    4e-12,
                );
            }
            close(
                sample.mutual_coherence_density_w.re,
                expected.mutual_coherence_density_w.re,
                4e-12,
            );
            close(
                sample.mutual_coherence_density_w.im,
                expected.mutual_coherence_density_w.im,
                4e-12,
            );
            close(sample.fringe_visibility, expected.fringe_visibility, 4e-12);
            assert_eq!(
                (sample.pupil_evaluations, sample.pair_contractions),
                (4, 28)
            );
            closure(sample);
        }
    }
}

#[test]
fn g1_polarization_beats_and_analyzers_match_independent_shutter_integrals() {
    let camera = focused();
    let delta = 4_294_967_296.0;
    let f1 = F0 + delta;
    let x = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let y = [jones_wave([0.0; 2], [C64::ZERO, C64::new(0.0, F0 / f1)])];
    let tones = [jones_tone(F0, &x, &[]), jones_tone(f1, &y, &[])];
    let select_x = layer([1.0, 0.0], [0.0; 2], 0.0);
    let select_y = layer([0.0, 1.0], [0.0; 2], 0.0);
    let analyzer = layer([1.0, 0.0], [0.0; 2], core::f64::consts::FRAC_PI_4);
    let hx = camera
        .in_focus_diffraction_amplitude([0.0; 2], wavelength(F0))
        .unwrap();
    let hy = camera
        .in_focus_diffraction_amplitude([0.0; 2], wavelength(f1))
        .unwrap()
        * y[0].amplitude_sqrt_w[1];
    let phase = 0.37;
    let path = C / (4.0 * F0);
    let receiver =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.5; 2], path, phase)
            .unwrap();
    for width in [0.05, 0.5, 1.0, 2.0] {
        let duration = width / delta;
        let start = -0.17 / delta;
        let shutter = CoherentDiffractionShutter::try_new(start, duration).unwrap();
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let sample = query
            .observe_interferometer(
                receiver,
                [arm(&[select_x]), arm(&[select_y])],
                RECEIVER_LIMITS,
                || false,
            )
            .unwrap();
        // Independent analytic integral of exp(+i*2*pi*delta*t). The common
        // carrier is removed algebraically; no production temporal kernel or
        // receiver contractions enter this reference.
        let sinc = if width == 1.0 || width == 2.0 {
            0.0
        } else {
            sinc_reference(core::f64::consts::PI * width)
        };
        let angle = core::f64::consts::TAU * delta * (start + 0.5 * duration - 0.03 / C)
            - phase
            - core::f64::consts::TAU * f1 * path / C;
        let cross = (hx * hy.conj()).scale(0.25 * sinc) * C64::new(angle.cos(), angle.sin());
        for port in 0..2 {
            let sign = if port == 0 { 1.0 } else { -1.0 };
            let expected = [
                0.25 * hx.norm_sq(),
                0.25 * hy.norm_sq(),
                sign * cross.re,
                sign * cross.im,
            ];
            matrix_close(sample.output_coherency_density_w[port], expected, 3e-12);
            matrix_close(
                sample.output_coherency_energy_density_j[port],
                expected.map(|v| v * duration),
                3e-12,
            );
            assert!(
                (sample.output_coherency_density_w[port].stokes().degree() - sinc.abs()).abs()
                    < 4e-12
            );
        }
        assert_eq!(sample.mutual_coherence_density_w, C64::ZERO);
        assert_eq!(sample.fringe_visibility, 0.0);
        assert_eq!(
            (sample.pupil_evaluations, sample.pair_contractions),
            (4, 52)
        );
        closure(sample);
        let analysed = query
            .observe_interferometer(
                receiver,
                [arm(&[select_x, analyzer]), arm(&[select_y, analyzer])],
                RECEIVER_LIMITS,
                || false,
            )
            .unwrap();
        let roundoff = 64.0 * f64::EPSILON * f64::EPSILON;
        assert!((analysed.fringe_visibility - sinc.abs()).abs() <= 5e-12 * sinc.abs() + roundoff);
        for port in 0..2 {
            let expected = Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4)
                .apply(sample.output_coherency_density_w[port].stokes());
            let actual = analysed.output_coherency_density_w[port].stokes();
            for (a, e) in [actual.i, actual.q, actual.u, actual.v]
                .into_iter()
                .zip([expected.i, expected.q, expected.u, expected.v])
            {
                assert!((a - e).abs() < 4e-12 * expected.i);
            }
        }
        closure(analysed);
    }
}

#[test]
fn g3_shared_columns_interfere_but_analyzers_do_not_correlate_independent_sources() {
    let camera = focused();
    let x = [jones_wave([0.0; 2], [C64::ONE, C64::ZERO])];
    let y = [jones_wave([0.0; 2], [C64::ZERO, C64::new(0.3, 0.7)])];
    let shared_x = [JonesPupilCoherenceFactor { waves: &x }];
    let shared_y = [JonesPupilCoherenceFactor { waves: &y }];
    let independent_x = [
        JonesPupilCoherenceFactor { waves: &x },
        JonesPupilCoherenceFactor { waves: &[] },
    ];
    let independent_y = [
        JonesPupilCoherenceFactor { waves: &[] },
        JonesPupilCoherenceFactor { waves: &y },
    ];
    let analyzer = layer([1.0, 0.0], [0.0; 2], core::f64::consts::FRAC_PI_4);
    let arms = [
        arm(&[layer([1.0, 0.0], [0.0; 2], 0.0), analyzer]),
        arm(&[layer([0.0, 1.0], [0.0; 2], 0.0), analyzer]),
    ];
    let shutter = CoherentDiffractionShutter::try_new(0.02, 0.01).unwrap();
    let observe = |tones: &[CoherentJonesDiffractionTone<'_>], phase| {
        camera
            .prepare_coherent_jones_diffraction_exposure(
                tones,
                0.03,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap()
            .observe_interferometer(receiver(phase), arms, RECEIVER_LIMITS, || false)
            .unwrap()
    };
    for phase in [0.0, 0.9, core::f64::consts::PI] {
        let mean = observe(
            &[jones_tone(F0, &x, &[]), jones_tone(F0 + 8.0, &y, &[])],
            phase,
        );
        let shared = observe(
            &[
                jones_tone(F0, &[], &shared_x),
                jones_tone(F0 + 8.0, &[], &shared_y),
            ],
            phase,
        );
        let independent = observe(
            &[
                jones_tone(F0, &[], &independent_x),
                jones_tone(F0 + 8.0, &[], &independent_y),
            ],
            phase,
        );
        for port in 0..2 {
            matrix_close(
                shared.output_coherency_density_w[port],
                matrix(mean.output_coherency_density_w[port]),
                3e-13,
            );
            let expected = 0.5 * independent.arm_power_density_w.iter().sum::<f64>();
            close(
                independent.output_coherency_density_w[port].stokes().i,
                expected,
                3e-13,
            );
        }
        assert_eq!(independent.mutual_coherence_density_w, C64::ZERO);
        assert_eq!(independent.fringe_visibility, 0.0);
        assert!(mean.fringe_visibility > 0.9);
        closure(independent);
    }
}

#[test]
fn g3_equal_lines_cancel_before_products_and_dark_ports_keep_weak_axis_power() {
    let camera = focused();
    let a = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1e-100)])];
    let negative = [jones_wave(
        [0.0; 2],
        [C64::new(-1.0, 0.0), C64::new(0.0, -1e-100)],
    )];
    let cancel_tones = [jones_tone(F0, &a, &[]), jones_tone(F0, &negative, &[])];
    let shutter = CoherentDiffractionShutter::try_new(2.0_f64.powi(40), 0.0625).unwrap();
    let cancelled = camera
        .prepare_coherent_jones_diffraction_exposure(
            &cancel_tones,
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    assert_eq!(cancelled.frequency_count(), 1);
    let sample = cancelled
        .observe_interferometer(
            receiver(0.7),
            [CoherentJonesArm::IDENTITY; 2],
            RECEIVER_LIMITS,
            || false,
        )
        .unwrap();
    assert_eq!(sample.arm_power_density_w, [0.0; 2]);
    assert!(
        sample
            .output_coherency_density_w
            .iter()
            .all(|c| c.stokes().i == 0.0)
    );
    assert_eq!(sample.pair_contractions, 14);

    let tones = [jones_tone(F0, &a, &[]), jones_tone(F0 + 16.0, &a, &[])];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let baseline = query
        .observe_coherency([0.0; 2], || false)
        .unwrap()
        .coherency_density_w;
    let phase = 1e-20;
    let sample = query
        .observe_interferometer(
            receiver(phase),
            [CoherentJonesArm::IDENTITY; 2],
            RECEIVER_LIMITS,
            || false,
        )
        .unwrap();
    for axis in 0..2 {
        relative(
            sample.output_coherency_density_w[0].axis_powers()[axis],
            baseline.axis_powers()[axis],
            3e-12,
        );
        relative(
            sample.output_coherency_density_w[1].axis_powers()[axis],
            baseline.axis_powers()[axis] * (0.25 * phase * phase),
            3e-11,
        );
    }
    let weak = sample.output_coherency_density_w[1];
    assert_eq!(weak.stokes().i, weak.stokes().q);
    assert!(weak.axis_powers()[1] > 0.0);
    closure(sample);
}

#[test]
fn g1_patterned_spectral_pupil_retains_native_modes_and_port_coherencies() {
    let camera = defocused(0.8);
    let a = [jones_wave(
        [0.1, -0.2],
        [C64::new(0.8, 0.2), C64::new(0.1, 0.5)],
    )];
    let b = [jones_wave(
        [-0.2, 0.1],
        [C64::new(0.3, -0.1), C64::new(-0.2, 0.4)],
    )];
    let tones = [jones_tone(F0, &a, &[]), jones_tone(1.03 * F0, &b, &[])];
    let mut query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            CoherentDiffractionShutter::try_new(0.0, 1e-13).unwrap(),
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let layers =
        [CoherentJonesPupilLayer::try_new([0.7, 0.2], [0.1, 0.8], 0.3, [800.0, -350.0]).unwrap()];
    let responses = [
        CoherentJonesPupilResponse {
            frequency_hz: F0,
            layers: &layers,
        },
        CoherentJonesPupilResponse {
            frequency_hz: 1.03 * F0,
            layers: &layers,
        },
    ];
    query
        .set_pupil_responses(
            &responses,
            CoherentJonesPupilLimits {
                maximum_layers: 2,
                maximum_modes: 12,
                maximum_pupil_evaluations: 6,
            },
            || false,
        )
        .unwrap();
    let at = [0.7 * SCALE, -0.3 * SCALE];
    let expected = query
        .observe_coherency(at, || false)
        .unwrap()
        .coherency_density_w;
    let receiver =
        CoherentDiffractionInterferometer::try_new([at; 2], 0.03, [0.5; 2], 0.0, 0.0).unwrap();
    let sample = query
        .observe_interferometer(
            receiver,
            [CoherentJonesArm::IDENTITY; 2],
            RECEIVER_LIMITS,
            || false,
        )
        .unwrap();
    matrix_close(
        sample.output_coherency_density_w[0],
        matrix(expected),
        5e-12,
    );
    assert_eq!(sample.output_coherency_density_w[1].stokes().i, 0.0);
    assert_eq!(sample.pupil_evaluations, 12);
    assert_eq!(sample.pair_contractions, 52);
    closure(sample);
}

#[test]
fn g0_g4_exact_active_work_limits_empty_receiver_and_cancelled_reuse() {
    let camera = focused();
    let a = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 0.5)])];
    let factors = [JonesPupilCoherenceFactor { waves: &a }];
    let tones = [
        jones_tone(F0, &a, &factors),
        jones_tone(F0 + 8.0, &a, &factors),
    ];
    let shutter = CoherentDiffractionShutter::try_new(-0.1, 0.02).unwrap();
    // A prepared ordinary Jones detector uses 12 real contractions here;
    // the explicit receiver limit admits the more expensive two-port query.
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            shutter,
            CoherentDiffractionExposureLimits {
                maximum_pair_contractions: 12,
                ..EXPOSURE_LIMITS
            },
            || false,
        )
        .unwrap();
    let exact = CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 8,
        maximum_pair_contractions: 104,
    };
    let r = receiver(0.4);
    let arms = [CoherentJonesArm::IDENTITY; 2];
    let expected = query
        .observe_interferometer(r, arms, exact, || false)
        .unwrap();
    assert_eq!(
        (expected.pupil_evaluations, expected.pair_contractions),
        (8, 104)
    );
    for bad in [
        CoherentInterferometerExposureLimits {
            maximum_pupil_evaluations: 7,
            ..exact
        },
        CoherentInterferometerExposureLimits {
            maximum_pair_contractions: 103,
            ..exact
        },
    ] {
        let mut polls = 0;
        assert!(matches!(
            query.observe_interferometer(r, arms, bad, || {
                polls += 1;
                false
            }),
            Err(CoherentDiffractionError::Limit { .. })
        ));
        assert_eq!(polls, 1, "work must refuse before pupil evaluation");
    }
    let single =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.0, 0.7], 0.1, 0.4)
            .unwrap();
    let single = query
        .observe_interferometer(
            single,
            arms,
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: 4,
                maximum_pair_contractions: 28,
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        (single.pupil_evaluations, single.pair_contractions),
        (4, 28)
    );
    assert_eq!(
        single.output_coherency_density_w[0],
        single.output_coherency_density_w[1]
    );
    assert_eq!(single.fringe_visibility, 0.0);
    closure(single);
    let closed =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.03, [0.0; 2], 0.0, 0.0)
            .unwrap();
    let zero = CoherentInterferometerExposureLimits {
        maximum_pupil_evaluations: 0,
        maximum_pair_contractions: 0,
    };
    let sample = query
        .observe_interferometer(closed, arms, zero, || false)
        .unwrap();
    assert_eq!((sample.pupil_evaluations, sample.pair_contractions), (0, 0));
    assert_eq!(sample.arm_power_density_w, [0.0; 2]);
    let empty = camera
        .prepare_coherent_jones_diffraction_exposure(&[], 0.03, shutter, EXPOSURE_LIMITS, || false)
        .unwrap();
    assert_eq!(
        empty
            .observe_interferometer(r, arms, zero, || false)
            .unwrap()
            .pupil_evaluations,
        0
    );
    let wrong =
        CoherentDiffractionInterferometer::try_new([[0.0; 2]; 2], 0.031, [0.5; 2], 0.0, 0.0)
            .unwrap();
    assert!(matches!(
        query.observe_interferometer(wrong, arms, exact, || false),
        Err(CoherentDiffractionError::InvalidSensorPlane)
    ));
    for stop in [1, 15, 40] {
        let mut polls = 0;
        assert!(matches!(
            query.observe_interferometer(r, arms, exact, || {
                polls += 1;
                polls >= stop
            }),
            Err(CoherentDiffractionError::Cancelled)
        ));
    }
    assert_eq!(
        query
            .observe_interferometer(r, arms, exact, || false)
            .unwrap(),
        expected
    );
}

#[test]
fn g0_balanced_port_normalization_keeps_large_finite_power_representable() {
    let camera = focused();
    let h = camera
        .in_focus_diffraction_amplitude([0.0; 2], wavelength(F0))
        .unwrap();
    let amplitude = (1.2e308_f64 / h.norm_sq()).sqrt();
    let a = [jones_wave([0.0; 2], [C64::new(amplitude, 0.0), C64::ZERO])];
    let tones = [jones_tone(F0, &a, &[])];
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &tones,
            0.03,
            CoherentDiffractionShutter::try_new(0.0, 0.1).unwrap(),
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let sample = query
        .observe_interferometer(
            receiver(0.0),
            [CoherentJonesArm::IDENTITY; 2],
            RECEIVER_LIMITS,
            || false,
        )
        .unwrap();
    relative(
        sample.output_coherency_density_w[0].axis_powers()[0],
        1.2e308,
        3e-13,
    );
    assert_eq!(sample.output_coherency_density_w[1].stokes().i, 0.0);
    assert_eq!(sample.fringe_visibility, 1.0);
    closure(sample);
}
