use super::*;
use fs_math::dd::Dd;
use fs_render::camera::{
    CoherentDiffractionInterferometer, CoherentInterferometerExposureLimits, CoherentJonesArm,
    CoherentJonesPupilLayer, CoherentJonesPupilLimits, CoherentJonesPupilResponse,
    PolarizedCoherentInterferometerPixelExposure, PolarizedDiffractionCoherency,
};
use fs_render::polarization::{Mueller, Stokes};

const C: f64 = 299_792_458.0;
const POINT_LIMITS: CoherentInterferometerExposureLimits = CoherentInterferometerExposureLimits {
    maximum_pupil_evaluations: 1024,
    maximum_pair_contractions: 100_000,
};

fn receiver(
    centers: [[f64; 2]; 2],
    distance: f64,
    powers: [f64; 2],
    phase: f64,
) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new(centers, distance, powers, 0.0, phase).unwrap()
}

fn layer(power: [f64; 2], phase: [f64; 2], angle: f64) -> CoherentJonesPupilLayer {
    CoherentJonesPupilLayer::try_new(power, phase, angle, [0.0; 2]).unwrap()
}

fn stokes_values(s: Stokes) -> [f64; 4] {
    [s.i, s.q, s.u, s.v]
}

fn compare_coherency(
    actual: PolarizedDiffractionCoherency,
    expected: PolarizedDiffractionCoherency,
) {
    let scale = expected.stokes().i.max(1e-280);
    let a = [
        actual.axis_powers()[0],
        actual.axis_powers()[1],
        actual.cross().re,
        actual.cross().im,
    ];
    let e = [
        expected.axis_powers()[0],
        expected.axis_powers()[1],
        expected.cross().re,
        expected.cross().im,
    ];
    for (actual, expected) in a.into_iter().zip(e) {
        assert!(
            (actual - expected).abs() <= 3e-11 * scale,
            "{actual:e} != {expected:e}"
        );
    }
}

fn closure(sample: PolarizedCoherentInterferometerPixelExposure, duration: f64) {
    let incoming = sample.arm_power_w.iter().sum::<f64>();
    let outgoing = sample
        .output_coherency_w
        .iter()
        .map(|c| c.stokes().i)
        .sum::<f64>();
    close(outgoing, incoming, 3e-12);
    for port in 0..2 {
        close(
            sample.arm_energy_j[port],
            sample.arm_power_w[port] * duration,
            3e-14,
        );
        let energy = stokes_values(sample.output_coherency_j[port].stokes());
        let power = stokes_values(sample.output_coherency_w[port].stokes());
        for (e, p) in energy.into_iter().zip(power) {
            assert!((e - p * duration).abs() <= 3e-14 * incoming * duration);
        }
    }
}

#[test]
fn g1_joint_receiver_matches_native_chromatic_patterned_pupil_pixel() {
    let camera = defocused(1.2);
    let f0 = C / LAMBDA;
    let f1 = f0 + 1e10;
    let mean0 = [jones_wave(
        [0.2, -0.1],
        [C64::new(0.7, 0.2), C64::new(-0.1, 0.3)],
    )];
    let mean1 = [jones_wave(
        [-0.3, 0.4],
        [C64::new(-0.2, 0.8), C64::new(0.6, -0.1)],
    )];
    let noise0 = [jones_wave(
        [-0.5, 0.2],
        [C64::new(0.2, 0.1), C64::new(0.1, -0.2)],
    )];
    let noise1 = [jones_wave(
        [0.1, -0.4],
        [C64::new(-0.1, 0.2), C64::new(0.3, 0.2)],
    )];
    let factors0 = [JonesPupilCoherenceFactor { waves: &noise0 }];
    let factors1 = [JonesPupilCoherenceFactor { waves: &noise1 }];
    let tones = [
        jones_tone(f0, &mean0, &factors0),
        jones_tone(f1, &mean1, &factors1),
    ];
    let shutter = CoherentDiffractionShutter::try_new(-0.17 / (f1 - f0), 0.23 / (f1 - f0)).unwrap();
    let mut query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let patterned = [CoherentJonesPupilLayer::try_new(
        [0.9, 0.6],
        [0.2, -0.1],
        0.31,
        [0.2 / RADIUS, -0.3 / RADIUS],
    )
    .unwrap()];
    query
        .set_pupil_responses(
            &[
                CoherentJonesPupilResponse {
                    frequency_hz: f0,
                    layers: &patterned,
                },
                CoherentJonesPupilResponse {
                    frequency_hz: f1,
                    layers: &patterned,
                },
            ],
            CoherentJonesPupilLimits {
                maximum_layers: 4,
                maximum_modes: 16,
                maximum_pupil_evaluations: 128,
            },
            || false,
        )
        .unwrap();
    let footprint = pixel([-1.1, -0.8], [1.3, 0.9]);
    let native = query
        .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let result = query
        .integrate_interferometer_pixel(
            receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0),
            [CoherentJonesArm::IDENTITY; 2],
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    compare_coherency(
        result.output_coherency_w[0],
        native.average_power.coherency_w,
    );
    compare_coherency(result.output_coherency_j[0], native.coherency_j);
    assert_eq!(result.output_coherency_w[1].stokes().i, 0.0);
    close(result.fringe_visibility, 1.0, 3e-14);
    assert_eq!(
        result.pupil_evaluations,
        2 * native.average_power.pupil_evaluations
    );
    closure(result, shutter.duration_s());
}

#[test]
fn g1_spatial_and_temporal_polarization_beats_match_independent_square_integrals() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let frequencies = [(1_u64 << 49) as f64, ((1_u64 << 49) + (1_u64 << 40)) as f64];
    let beat = frequencies[1] - frequencies[0];
    let distance = 50.0;
    let centers = [[-0.7 * SCALE, 0.0], [0.7 * SCALE, 0.0]];
    let left = [CoherentJonesPupilWave {
        tangent_direction: centers[0],
        amplitude_sqrt_w: [C64::ONE, C64::ZERO],
    }];
    let right = [CoherentJonesPupilWave {
        tangent_direction: centers[1],
        amplitude_sqrt_w: [C64::ZERO, C64::ONE],
    }];
    let tones = [
        jones_tone(frequencies[0], &left, &[]),
        jones_tone(frequencies[1], &right, &[]),
    ];
    let arms = [
        CoherentJonesArm::try_new(&[layer([1.0, 0.0], [0.0; 2], 0.0)]).unwrap(),
        CoherentJonesArm::try_new(&[layer([0.0, 1.0], [0.0; 2], 0.0)]).unwrap(),
    ];
    let low = [-2.8, -1.1];
    let high = [3.1, 1.3];
    let footprint = pixel(low, high);
    let phase = 0.37;
    let limits = CoherentDiffractionExposurePixelLimits {
        maximum_pair_contractions: 10_000_000,
        ..EXPOSURE_PIXEL_LIMITS
    };
    for periods in [0.25, 1.0] {
        let shutter = CoherentDiffractionShutter::try_new(distance / C, periods / beat).unwrap();
        let query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &tones,
                distance,
                shutter,
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        let result = query
            .integrate_interferometer_pixel(
                receiver(centers, distance, [0.5; 2], phase),
                arms,
                footprint,
                limits,
                || false,
            )
            .unwrap();
        // Independent analytic square-pupil transform and analytic rectangular
        // shutter, with composite Simpson spatial integration. Neither the
        // production pupil, interferometer nor temporal kernel is this oracle.
        let shape = |f: f64, q: f64| {
            sinc_reference(
                core::f64::consts::TAU * RADIUS * f / C
                    * SCALE
                    * q
                    * core::f64::consts::FRAC_1_SQRT_2,
            )
        };
        let mut axes = frequencies.map(|f| 0.25 * 2.0 * (RADIUS * f / C).powi(2) * SCALE * SCALE);
        for axis in 0..2 {
            for component in 0..2 {
                axes[component] *= integrate_reference(low[axis], high[axis], |q| {
                    shape(frequencies[component], q).powi(2)
                });
            }
        }
        let mut cross = C64::from_re(
            0.5 * RADIUS * RADIUS * frequencies[0] / C * frequencies[1] / C * SCALE * SCALE,
        );
        for axis in 0..2 {
            let field_phase = |q: f64| {
                let u = q * SCALE;
                core::f64::consts::PI * distance / C
                    * (-beat * u * u
                        + 2.0
                            * (frequencies[0] * centers[0][axis]
                                - frequencies[1] * centers[1][axis])
                            * u)
            };
            cross = cross
                * C64::new(
                    integrate_reference(low[axis], high[axis], |q| {
                        shape(frequencies[0], q) * shape(frequencies[1], q) * field_phase(q).cos()
                    }),
                    integrate_reference(low[axis], high[axis], |q| {
                        shape(frequencies[0], q) * shape(frequencies[1], q) * field_phase(q).sin()
                    }),
                );
        }
        let time = (Dd::from_f64(shutter.start_s()) - Dd::from_f64(distance) / Dd::from_f64(C)
            + Dd::from_f64(shutter.duration_s()) * Dd::from_f64(0.5))
            * Dd::from_f64(beat);
        let center_phase = core::f64::consts::TAU * time.to_f64()
            - core::f64::consts::PI * distance / C * beat * centers[1][0].powi(2)
            - phase;
        let temporal = if periods == 1.0 {
            0.0
        } else {
            sinc_reference(core::f64::consts::PI * periods)
        };
        cross = cross * C64::new(center_phase.cos(), center_phase.sin()) * C64::from_re(temporal);
        for port in 0..2 {
            let actual = result.output_coherency_w[port];
            let sign = if port == 0 { 1.0 } else { -1.0 };
            for (a, e) in actual.axis_powers().into_iter().zip(axes) {
                close(a, e, 3e-10);
            }
            let scale = axes.iter().sum::<f64>();
            assert!((actual.cross().re - sign * cross.re).abs() < 2e-9 * scale);
            assert!((actual.cross().im - sign * cross.im).abs() < 2e-9 * scale);
            assert!(actual.cross().abs() < 0.2 * (axes[0] * axes[1]).sqrt());
        }
        assert_eq!(result.mutual_coherence_w, C64::ZERO);
        assert_eq!(result.fringe_visibility, 0.0);
        assert!(result.subdivisions[0] > 2 * result.subdivisions[1]);
        closure(result, shutter.duration_s());
    }
}

#[test]
fn g3_joint_pixel_partition_and_passive_analyzer_preserve_full_stokes() {
    let camera = focused();
    let f0 = C / LAMBDA;
    let f1 = f0 + 2e11;
    let a = [jones_wave([-0.2, 0.3], [C64::ONE, C64::new(0.1, 0.7)])];
    let b = [jones_wave(
        [0.4, -0.1],
        [C64::new(-0.3, 0.6), C64::new(0.8, -0.2)],
    )];
    let tones = [jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])];
    let shutter = CoherentDiffractionShutter::try_new(-0.1 / (f1 - f0), 0.3 / (f1 - f0)).unwrap();
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let detector = receiver(
        [[-0.6 * SCALE, 0.2 * SCALE], [0.7 * SCALE, -0.1 * SCALE]],
        0.03,
        [0.35, 0.55],
        0.6,
    );
    let footprint = pixel([-1.2, -0.8], [1.4, 1.1]);
    let plain = query
        .integrate_interferometer_pixel(
            detector,
            [CoherentJonesArm::IDENTITY; 2],
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    let angle = 0.41;
    let retardance = 0.73;
    let transformed = CoherentJonesArm::try_new(&[
        layer([1.0; 2], [0.0, retardance], 0.0),
        layer([1.0, 0.0], [0.0; 2], angle),
    ])
    .unwrap();
    let result = query
        .integrate_interferometer_pixel(
            detector,
            [transformed; 2],
            footprint,
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    for port in 0..2 {
        let expected = Mueller::linear_polarizer(angle).apply(
            Mueller::linear_retarder(0.0, retardance)
                .apply(plain.output_coherency_w[port].stokes()),
        );
        for (actual, expected) in stokes_values(result.output_coherency_w[port].stokes())
            .into_iter()
            .zip(stokes_values(expected))
        {
            assert!((actual - expected).abs() < 4e-12 * plain.output_coherency_w[port].stokes().i);
        }
    }
    assert!(result.arm_power_w.iter().sum::<f64>() < plain.arm_power_w.iter().sum::<f64>());
    let left = query
        .integrate_interferometer_pixel(
            detector,
            [transformed; 2],
            pixel([-1.2, -0.8], [0.13, 1.1]),
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    let right = query
        .integrate_interferometer_pixel(
            detector,
            [transformed; 2],
            pixel([0.13, -0.8], [1.4, 1.1]),
            EXPOSURE_PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    for port in 0..2 {
        let a = stokes_values(left.output_coherency_w[port].stokes());
        let b = stokes_values(right.output_coherency_w[port].stokes());
        let whole = stokes_values(result.output_coherency_w[port].stokes());
        for i in 0..4 {
            assert!((a[i] + b[i] - whole[i]).abs() < 3e-11 * whole[0]);
        }
        close(
            left.arm_power_w[port] + right.arm_power_w[port],
            result.arm_power_w[port],
            3e-11,
        );
    }
    let mutual = left.mutual_coherence_w + right.mutual_coherence_w;
    assert!(
        (mutual - result.mutual_coherence_w).abs() < 3e-11 * plain.arm_power_w.iter().sum::<f64>()
    );
    closure(result, shutter.duration_s());
}

#[test]
fn g3_joint_pixel_retains_weak_dark_channels_and_preflights_all_work() {
    let camera = focused();
    let source = [jones_wave([0.1, -0.2], [C64::ONE, C64::from_re(1e-100)])];
    let tones = [jones_tone(C / LAMBDA, &source, &[])];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.02).unwrap();
    let query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let phase = 1e-20;
    let detector = receiver([[0.0; 2]; 2], 0.03, [0.5; 2], phase);
    let footprint = pixel([-0.8, -0.6], [0.9, 0.7]);
    let arms = [CoherentJonesArm::IDENTITY; 2];
    let result = query
        .integrate_interferometer_pixel(detector, arms, footprint, EXPOSURE_PIXEL_LIMITS, || false)
        .unwrap();
    let bright = result.output_coherency_w[0].axis_powers();
    let dark = result.output_coherency_w[1].axis_powers();
    for axis in 0..2 {
        close(
            dark[axis] / bright[axis] / (0.25 * phase * phase),
            1.0,
            3e-12,
        );
    }
    close(bright[1] / bright[0] / 1e-200, 1.0, 3e-12);
    let point = query
        .observe_interferometer(detector, arms, POINT_LIMITS, || false)
        .unwrap();
    assert_eq!(
        result.pupil_evaluations,
        point.pupil_evaluations * result.quadrature_samples
    );
    assert_eq!(
        result.pair_contractions,
        point.pair_contractions * result.quadrature_samples
    );
    let exact = CoherentDiffractionExposurePixelLimits {
        spatial: CoherentDiffractionPixelLimits {
            maximum_samples: result.quadrature_samples,
            maximum_pupil_evaluations: result.pupil_evaluations,
            ..PIXEL_LIMITS
        },
        maximum_pair_contractions: result.pair_contractions,
    };
    assert_eq!(
        query
            .integrate_interferometer_pixel(detector, arms, footprint, exact, || false)
            .unwrap(),
        result
    );
    for which in 0..3 {
        let mut limits = exact;
        match which {
            0 => limits.spatial.maximum_samples -= 1,
            1 => limits.spatial.maximum_pupil_evaluations -= 1,
            _ => limits.maximum_pair_contractions -= 1,
        }
        assert!(matches!(
            query.integrate_interferometer_pixel(detector, arms, footprint, limits, || false),
            Err(CoherentDiffractionError::Limit { .. })
        ));
    }
    let mut polls = 0;
    assert_eq!(
        query.integrate_interferometer_pixel(detector, arms, footprint, exact, || {
            polls += 1;
            polls > 23
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query
            .integrate_interferometer_pixel(detector, arms, footprint, exact, || false)
            .unwrap(),
        result
    );
    let closed = query
        .integrate_interferometer_pixel(
            receiver([[0.0; 2]; 2], 0.03, [0.0; 2], 0.0),
            arms,
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
    assert_eq!(closed.output_coherency_w[0].stokes().i, 0.0);
    assert_eq!(
        (
            closed.quadrature_samples,
            closed.pupil_evaluations,
            closed.pair_contractions
        ),
        (0, 0, 0)
    );
    assert_eq!(closed.subdivisions, [0; 2]);
    closure(result, shutter.duration_s());
}
