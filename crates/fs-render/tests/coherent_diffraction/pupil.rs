//! Physical-pupil references for passive patterned Jones layers.
use super::*;
use fs_render::camera::{
    CoherentJonesPupilLayer, CoherentJonesPupilLimits, CoherentJonesPupilResponse,
};

const PUPIL_LIMITS: CoherentJonesPupilLimits = CoherentJonesPupilLimits {
    maximum_layers: 16,
    maximum_modes: 128,
    maximum_pupil_evaluations: 5000,
};
type Layer = ([f64; 2], [f64; 2], f64, [f64; 2]);

fn layer(spec: Layer) -> CoherentJonesPupilLayer {
    CoherentJonesPupilLayer::try_new(spec.0, spec.1, spec.2, spec.3).unwrap()
}

fn response(
    frequency_hz: f64,
    layers: &[CoherentJonesPupilLayer],
) -> CoherentJonesPupilResponse<'_> {
    CoherentJonesPupilResponse {
        frequency_hz,
        layers,
    }
}

// Direct physical principal-axis rotation, attenuation and phase. No Fourier
// modes, original pupil kernel or production matrix algebra enter this reference.
fn local_layer(field: [C64; 2], spec: Layer, at_m: [f64; 2]) -> [C64; 2] {
    let angle = spec.2 + spec.3[0] * at_m[0] + spec.3[1] * at_m[1];
    let (s, c) = angle.sin_cos();
    let local = [
        field[0].scale(c) + field[1].scale(s),
        field[1].scale(c) - field[0].scale(s),
    ];
    let value: [C64; 2] = std::array::from_fn(|axis| {
        local[axis] * C64::new(spec.1[axis].cos(), spec.1[axis].sin()).scale(spec.0[axis].sqrt())
    });
    [
        value[0].scale(c) - value[1].scale(s),
        value[0].scale(s) + value[1].scale(c),
    ]
}

fn project(field: [C64; 2], mode: [C64; 2]) -> C64 {
    (mode[0].conj() * field[0] + mode[1].conj() * field[1])
        .scale(1.0 / (mode[0].norm_sq() + mode[1].norm_sq()).sqrt())
}

#[test]
fn g1_uniform_jones_pupil_preserves_layer_order_and_weak_dichroic_transmission() {
    let camera = defocused(2.4);
    let f = 299_792_458.0 / LAMBDA;
    let waves = [jones_wave(
        [0.2, -0.1],
        [C64::new(0.6, 0.2), C64::new(-0.3, 0.7)],
    )];
    let tones = [jones_tone(f, &waves, &[])];
    let shutter = CoherentDiffractionShutter::try_new(0.0, 0.02).unwrap();
    let mut query = camera
        .prepare_coherent_jones_diffraction_exposure(&tones, 0.03, shutter, EXPOSURE_LIMITS, || {
            false
        })
        .unwrap();
    let at = [0.7 * SCALE, -0.4 * SCALE];
    let h = camera
        .defocused_diffraction_amplitude([0.5 * SCALE, -0.3 * SCALE], LAMBDA_NM)
        .unwrap();
    let specs = [
        ([1.0, 0.0], [0.0, 0.0], 0.0, [0.0; 2]),
        (
            [1.0, 1.0],
            [0.0, core::f64::consts::FRAC_PI_2],
            core::f64::consts::FRAC_PI_4,
            [0.0; 2],
        ),
    ];
    // Avoid an ideal dark analyzer here: an f64 direct-rotation reference
    // cannot resolve the same ~1e-32 residual as the production DD path.
    let mode = [C64::ONE, C64::new(0.0, 0.6)];
    let detector = CoherentJonesDetector::analyzer(mode).unwrap();
    let mut measured = [0.0; 2];
    for (case, order) in [[specs[0], specs[1]], [specs[1], specs[0]]]
        .into_iter()
        .enumerate()
    {
        let layers = order.map(layer);
        query
            .set_pupil_responses(&[response(f, &layers)], PUPIL_LIMITS, || false)
            .unwrap();
        let transformed = order
            .into_iter()
            .fold(waves[0].amplitude_sqrt_w, |field, spec| {
                local_layer(field, spec, [0.0; 2])
            });
        let expected = h.norm_sq() * project(transformed, mode).norm_sq();
        measured[case] = query
            .observe(at, detector, || false)
            .unwrap()
            .power_density_w;
        close(measured[case], expected, 3e-12);
    }
    assert!((measured[0] - measured[1]).abs() > 0.1 * measured[0].max(measured[1]));

    let vertical = [jones_wave([0.0; 2], [C64::ZERO, C64::ONE])];
    let mut weak = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[jones_tone(f, &vertical, &[])],
            0.03,
            shutter,
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let clear = weak
        .observe(at, JONES_TOTAL, || false)
        .unwrap()
        .power_density_w;
    let attenuator = [layer(([1.0, 1e-200], [0.0; 2], 0.0, [0.0; 2]))];
    weak.set_pupil_responses(&[response(f, &attenuator)], PUPIL_LIMITS, || false)
        .unwrap();
    let transmitted = weak
        .observe(at, JONES_TOTAL, || false)
        .unwrap()
        .power_density_w;
    assert!(transmitted > 0.0);
    assert!((transmitted / clear / 1e-200 - 1.0).abs() < 3e-12);
    weak.set_pupil_responses(&[response(f, &[])], PUPIL_LIMITS, || false)
        .unwrap();
    close(
        weak.observe(at, JONES_TOTAL, || false)
            .unwrap()
            .power_density_w,
        clear,
        3e-13,
    );
}

fn square_integral(
    waves: &[CoherentJonesPupilWave],
    at: [f64; 2],
    phase: f64,
    layers: &[Layer],
    n: usize,
) -> [C64; 2] {
    let a = RADIUS * core::f64::consts::FRAC_1_SQRT_2;
    let step = 2.0 * a / n as f64;
    let k = core::f64::consts::TAU / LAMBDA;
    let weight = |i| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let mut sum = [C64::ZERO; 2];
    for iy in 0..=n {
        for ix in 0..=n {
            let u = [-a + ix as f64 * step, -a + iy as f64 * step];
            let mut incident = [C64::ZERO; 2];
            for wave in waves {
                let phase =
                    k * (wave.tangent_direction[0] * u[0] + wave.tangent_direction[1] * u[1]);
                let rotation = C64::new(phase.cos(), phase.sin());
                for (value, amplitude) in incident.iter_mut().zip(wave.amplitude_sqrt_w) {
                    *value = *value + amplitude * rotation;
                }
            }
            let transmitted = layers
                .iter()
                .fold(incident, |field, &spec| local_layer(field, spec, u));
            let phase = phase * (u[0] * u[0] + u[1] * u[1]) / (RADIUS * RADIUS)
                - k * (at[0] * u[0] + at[1] * u[1]);
            let propagation = C64::new(phase.cos(), phase.sin()).scale(weight(ix) * weight(iy));
            for axis in 0..2 {
                sum[axis] = sum[axis] + transmitted[axis] * propagation;
            }
        }
    }
    // Input field is a/sqrt(area); Fresnel amplitude supplies 1/lambda.
    sum.map(|value| value.scale(step * step / (9.0 * LAMBDA * 2.0 * a)))
}

#[test]
fn g1_patterned_jones_pupil_converges_to_independent_physical_aperture_integral() {
    let f = 299_792_458.0 / LAMBDA;
    let waves = [
        jones_wave([-1.1, 0.2], [C64::new(0.8, -0.1), C64::new(0.3, 0.2)]),
        jones_wave([0.6, -0.7], [C64::new(-0.2, 0.4), C64::new(0.6, -0.1)]),
    ];
    let specs = [
        ([0.9, 0.35], [0.1, 1.3], 0.2, [1.8 / RADIUS, -0.7 / RADIUS]),
        (
            [0.65, 1.0],
            [-0.4, 0.8],
            -0.3,
            [-0.4 / RADIUS, 1.1 / RADIUS],
        ),
    ];
    let layers = specs.map(layer);
    let at = [0.7 * SCALE, -0.4 * SCALE];
    let mode = [C64::new(0.4, 0.1), C64::new(0.2, -0.8)];
    let detector = CoherentJonesDetector::analyzer(mode).unwrap();
    for phase in [0.0, 2.1, -2.1] {
        let depth = 1.0 / (0.5 + phase * LAMBDA / (core::f64::consts::PI * RADIUS * RADIUS));
        let camera =
            base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
                .with_defocused_diffraction(
                    DefocusedDiffraction::try_new(Some(depth), 0.0).unwrap(),
                )
                .unwrap();
        let mut query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &[jones_tone(f, &waves, &[])],
                0.03,
                CoherentDiffractionShutter::try_new(0.0, 0.01).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        query
            .set_pupil_responses(&[response(f, &layers)], PUPIL_LIMITS, || false)
            .unwrap();
        let actual = query
            .observe(at, detector, || false)
            .unwrap()
            .power_density_w;
        let estimates = [32, 64, 128]
            .map(|n| project(square_integral(&waves, at, phase, &specs, n), mode).norm_sq());
        let errors = estimates.map(|value| (value - actual).abs());
        assert!(
            errors[1] < 0.2 * errors[0] && errors[2] < 0.2 * errors[1],
            "phase {phase}: physical integral convergence {errors:?}"
        );
        close(actual, estimates[2], 3e-7);
        let field = square_integral(&waves, at, phase, &specs, 128);
        close(
            query
                .observe(at, JONES_TOTAL, || false)
                .unwrap()
                .power_density_w,
            field[0].norm_sq() + field[1].norm_sq(),
            3e-7,
        );
        // The concurrent full-coherency owner must see exactly the same
        // patterned fields and statistics as both detector readouts.
        let coherency = query
            .observe_coherency(at, || false)
            .unwrap()
            .coherency_density_w;
        let [x, y] = coherency.axis_powers();
        close(
            x + y,
            query
                .observe(at, JONES_TOTAL, || false)
                .unwrap()
                .power_density_w,
            3e-13,
        );
        let analyzed = (x * mode[0].norm_sq()
            + y * mode[1].norm_sq()
            + 2.0 * (mode[0].conj() * coherency.cross() * mode[1]).re)
            / (mode[0].norm_sq() + mode[1].norm_sq());
        close(analyzed, actual, 3e-13);
        close(x, field[0].norm_sq(), 3e-7);
        close(y, field[1].norm_sq(), 3e-7);
        let cross = field[0] * field[1].conj();
        close(coherency.cross().re, cross.re, 3e-7);
        close(coherency.cross().im, cross.im, 3e-7);
    }
}

#[test]
fn g1_chromatic_jones_pupil_retains_common_phase_in_finite_exposure_beats() {
    let camera = defocused(2.4);
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = 1.13 * f0;
    let beat = f1 - f0;
    let distance = 0.037;
    let a = [jones_wave(
        [-0.7, 0.1],
        [C64::new(0.8, -0.2), C64::new(0.3, 0.4)],
    )];
    let b = [jones_wave(
        [0.6, -0.2],
        [C64::new(-0.2, 0.7), C64::new(0.6, -0.3)],
    )];
    let specs = [
        ([0.85, 0.4], [0.7, -0.2], 0.37, [0.0; 2]),
        ([1.0, 0.25], [-0.3, 1.1], -0.21, [0.0; 2]),
    ];
    let first = [layer(specs[0])];
    let second = [layer(specs[1])];
    let mode = [C64::new(0.8, 0.1), C64::new(-0.2, 0.5)];
    let detector = CoherentJonesDetector::analyzer(mode).unwrap();
    let at = [0.2 * SCALE, -0.3 * SCALE];
    let field = |wave: CoherentJonesPupilWave, frequency: f64, spec| {
        let h = camera
            .defocused_diffraction_amplitude(
                std::array::from_fn(|axis| at[axis] - wave.tangent_direction[axis]),
                299_792_458.0 / frequency * 1e9,
            )
            .unwrap();
        project(local_layer(wave.amplitude_sqrt_w, spec, [0.0; 2]), mode) * h
    };
    let x = field(a[0], f0, specs[0]);
    let y = field(b[0], f1, specs[1]);
    for periods in [0.2, 1.0, 3.7] {
        let start = 0.13 / beat;
        let duration = periods / beat;
        let mut query = camera
            .prepare_coherent_jones_diffraction_exposure(
                &[jones_tone(f0, &a, &[]), jones_tone(f1, &b, &[])],
                distance,
                CoherentDiffractionShutter::try_new(start, duration).unwrap(),
                EXPOSURE_LIMITS,
                || false,
            )
            .unwrap();
        // Frequency labels, not caller ordering, select the material response.
        query
            .set_pupil_responses(
                &[response(f1, &second), response(f0, &first)],
                PUPIL_LIMITS,
                || false,
            )
            .unwrap();
        let phase = core::f64::consts::TAU
            * beat
            * (distance * (1.0 + 0.5 * (at[0] * at[0] + at[1] * at[1])) / 299_792_458.0
                - start
                - 0.5 * duration);
        let expected = x.norm_sq()
            + y.norm_sq()
            + 2.0
                * (y * x.conj() * C64::new(phase.cos(), phase.sin())).re
                * sinc_reference(core::f64::consts::PI * beat * duration);
        close(
            query
                .observe(at, detector, || false)
                .unwrap()
                .power_density_w,
            expected,
            3e-10,
        );
    }
}

#[test]
fn g3_polarization_grating_steers_chromatically_and_preserves_pixel_power() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let period = 20.0 * LAMBDA;
    let screen = [layer((
        [1.0; 2],
        [0.0, core::f64::consts::PI],
        0.21,
        [core::f64::consts::PI / period, 0.0],
    ))];
    let i = integrate_reference(-5.0, 5.0, |q| {
        sinc_reference(q * core::f64::consts::FRAC_1_SQRT_2).powi(2)
    });
    let expected_pixel = i * i / (2.0 * core::f64::consts::PI.powi(2));
    for frequency in [299_792_458.0 / LAMBDA, 1.1 * 299_792_458.0 / LAMBDA] {
        let wavelength = 299_792_458.0 / frequency;
        let scale = wavelength / (core::f64::consts::TAU * RADIUS);
        for spin in [-1.0, 1.0] {
            let waves = [jones_wave(
                [0.0; 2],
                [
                    C64::new(core::f64::consts::FRAC_1_SQRT_2, 0.0),
                    C64::new(0.0, spin * core::f64::consts::FRAC_1_SQRT_2),
                ],
            )];
            let mut query = camera
                .prepare_coherent_jones_diffraction_exposure(
                    &[jones_tone(frequency, &waves, &[])],
                    0.03,
                    CoherentDiffractionShutter::try_new(0.0, 0.03).unwrap(),
                    EXPOSURE_LIMITS,
                    || false,
                )
                .unwrap();
            query
                .set_pupil_responses(&[response(frequency, &screen)], PUPIL_LIMITS, || false)
                .unwrap();
            let converted =
                CoherentJonesDetector::analyzer([C64::ONE, C64::new(0.0, -spin)]).unwrap();
            let shift = spin * wavelength / period;
            let peak = query
                .observe([shift, 0.0], converted, || false)
                .unwrap()
                .power_density_w;
            close(
                peak,
                2.0 * RADIUS * RADIUS / (wavelength * wavelength),
                3e-12,
            );
            assert!(
                peak > 1e4
                    * query
                        .observe([-shift, 0.0], converted, || false)
                        .unwrap()
                        .power_density_w
            );
            let footprint = CoherentDiffractionPixel::try_new(
                [shift - 5.0 * scale, -5.0 * scale],
                [shift + 5.0 * scale, 5.0 * scale],
            )
            .unwrap();
            let energy = query
                .integrate_pixel(footprint, converted, EXPOSURE_PIXEL_LIMITS, || false)
                .unwrap();
            close(energy.average_power.power_w, expected_pixel, 3e-10);
            close(energy.energy_j, 0.03 * expected_pixel, 3e-10);
            assert!(energy.average_power.power_w < 1.0);
            assert_eq!(
                energy.average_power.pupil_evaluations,
                3 * energy.average_power.quadrature_samples
            );
            let full = query
                .integrate_pixel_coherency(footprint, EXPOSURE_PIXEL_LIMITS, || false)
                .unwrap();
            let powers = full.average_power.coherency_w.axis_powers();
            close(powers[0] + powers[1], energy.average_power.power_w, 3e-12);
            let joules = full.coherency_j.axis_powers();
            close(joules[0] + joules[1], energy.energy_j, 3e-12);
            assert_eq!(
                full.average_power.pupil_evaluations,
                3 * full.average_power.quadrature_samples
            );
        }
    }
}

#[test]
fn g0_g4_jones_pupil_budgets_and_failed_replacement_preserve_original_optics() {
    for spec in [
        ([-0.1, 1.0], [0.0; 2], 0.0, [0.0; 2]),
        ([1.1, 1.0], [0.0; 2], 0.0, [0.0; 2]),
        ([1.0; 2], [0.0, f64::NAN], 0.0, [0.0; 2]),
        ([1.0; 2], [0.0; 2], f64::INFINITY, [0.0; 2]),
        ([1.0; 2], [0.0; 2], 0.0, [f64::NAN, 0.0]),
    ] {
        assert_eq!(
            CoherentJonesPupilLayer::try_new(spec.0, spec.1, spec.2, spec.3),
            Err(CoherentDiffractionError::InvalidJonesPupil)
        );
    }
    let camera = focused();
    let f0 = 299_792_458.0 / LAMBDA;
    let f1 = f0 + 1e9;
    let waves = [jones_wave([0.0; 2], [C64::ONE, C64::new(0.0, 1.0)])];
    let mut query = camera
        .prepare_coherent_jones_diffraction_exposure(
            &[jones_tone(f0, &waves, &[]), jones_tone(f1, &waves, &[])],
            0.03,
            CoherentDiffractionShutter::try_new(0.0, 0.2e-9).unwrap(),
            EXPOSURE_LIMITS,
            || false,
        )
        .unwrap();
    let baseline = query.observe([0.0; 2], JONES_TOTAL, || false).unwrap();
    let screen = [layer(([0.8, 0.4], [0.1, 1.2], 0.2, [1.0 / RADIUS, 0.0]))];
    let responses = [response(f0, &screen), response(f1, &screen)];
    for bad in [
        vec![response(f0, &screen)],
        vec![response(f0, &screen), response(f0, &screen)],
        vec![response(f0, &screen), response(f1 + 1.0, &screen)],
    ] {
        assert_eq!(
            query.set_pupil_responses(&bad, PUPIL_LIMITS, || false),
            Err(CoherentDiffractionError::InconsistentJonesPupilSpectrum)
        );
        assert_eq!(
            query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
            baseline
        );
    }
    for limits in [
        CoherentJonesPupilLimits {
            maximum_layers: 1,
            ..PUPIL_LIMITS
        },
        CoherentJonesPupilLimits {
            maximum_modes: 4,
            ..PUPIL_LIMITS
        },
        CoherentJonesPupilLimits {
            maximum_pupil_evaluations: 5,
            ..PUPIL_LIMITS
        },
    ] {
        assert!(matches!(
            query.set_pupil_responses(&responses, limits, || false),
            Err(CoherentDiffractionError::Limit { .. })
        ));
        assert_eq!(
            query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
            baseline
        );
    }
    let excessive = [layer(([1.0; 2], [0.0, 1.0], 0.0, [1e8 / RADIUS, 0.0]))];
    assert_eq!(
        query.set_pupil_responses(
            &[response(f0, &excessive), response(f1, &screen)],
            PUPIL_LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
    let mut polls = 0;
    assert_eq!(
        query.set_pupil_responses(&responses, PUPIL_LIMITS, || {
            polls += 1;
            polls == 14
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
        baseline
    );
    query
        .set_pupil_responses(&responses, PUPIL_LIMITS, || false)
        .unwrap();
    let masked = query.observe([0.0; 2], JONES_TOTAL, || false).unwrap();
    let underflow = [layer(([1e-200; 2], [0.0; 2], 0.0, [0.0; 2])); 4];
    assert_eq!(
        query.set_pupil_responses(
            &[response(f0, &underflow), response(f1, &underflow)],
            PUPIL_LIMITS,
            || false
        ),
        Err(CoherentDiffractionError::Unrepresentable)
    );
    assert_eq!(
        query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
        masked
    );
    let mut polls = 0;
    assert!(matches!(
        query.integrate_pixel(
            pixel([-0.1; 2], [0.1; 2]),
            JONES_TOTAL,
            CoherentDiffractionExposurePixelLimits {
                spatial: CoherentDiffractionPixelLimits {
                    maximum_pupil_evaluations: 0,
                    ..PIXEL_LIMITS
                },
                ..EXPOSURE_PIXEL_LIMITS
            },
            || {
                polls += 1;
                false
            }
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel pupil evaluations",
            ..
        })
    ));
    assert_eq!(polls, 1);
    let mut polls = 0;
    assert_eq!(
        query.observe([0.0; 2], JONES_TOTAL, || {
            polls += 1;
            polls == 11
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query.observe([0.0; 2], JONES_TOTAL, || false).unwrap(),
        masked
    );

    // Exactly crossed principal channels block all light, without a numerical
    // clipping threshold or wasted zero-field pupil evaluations.
    let crossed = [
        layer(([1.0, 0.0], [0.0; 2], 0.0, [0.0; 2])),
        layer(([0.0, 1.0], [0.0; 2], 0.0, [0.0; 2])),
    ];
    query
        .set_pupil_responses(
            &[response(f0, &crossed), response(f1, &crossed)],
            CoherentJonesPupilLimits {
                maximum_modes: 1,
                maximum_pupil_evaluations: 0,
                ..PUPIL_LIMITS
            },
            || false,
        )
        .unwrap();
    let dark = query
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
    assert_eq!(dark.energy_j, 0.0);
    assert_eq!(dark.average_power.pupil_evaluations, 0);
}
