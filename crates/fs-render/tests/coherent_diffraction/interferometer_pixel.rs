use super::*;
use fs_render::camera::{CoherentDiffractionInterferometer, CoherentInterferometerPixelPower};

fn receiver(
    coordinates: [[f64; 2]; 2],
    distance: f64,
    powers: [f64; 2],
    path: f64,
    phase: f64,
) -> CoherentDiffractionInterferometer {
    CoherentDiffractionInterferometer::try_new(coordinates, distance, powers, path, phase).unwrap()
}

#[test]
fn g1_coincident_detector_integrates_native_power_and_retains_weak_dark_port() {
    let camera = focused();
    let mean = [wave([0.2, -0.1], C64::new(1.0, 0.3))];
    let noise = [wave([-0.4, 0.6], C64::new(0.4, -0.2))];
    let factors = [PupilCoherenceFactor { waves: &noise }];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &factors,
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let footprint = pixel([-1.3, -0.9], [1.1, 0.7]);
    let native = query
        .integrate_pixel(footprint, PIXEL_LIMITS, || false)
        .unwrap();
    for phase in [0.0_f64, 0.7, 1e-20] {
        let result = query
            .integrate_interferometer_pixel(
                receiver([[0.0; 2]; 2], 0.03, [0.5; 2], 0.0, phase),
                footprint,
                PIXEL_LIMITS,
                || false,
            )
            .unwrap();
        close(result.arm_power_w[0], 0.5 * native.power_w, 3e-12);
        close(result.arm_power_w[1], 0.5 * native.power_w, 3e-12);
        close(result.output_power_w.iter().sum(), native.power_w, 3e-12);
        close(
            result.output_power_w[0],
            native.power_w * (0.5 * phase).cos().powi(2),
            3e-12,
        );
        close(result.fringe_visibility, 1.0, 1e-14);
        assert_eq!(result.pupil_evaluations, 4 * result.quadrature_samples);
        if phase == 0.0 {
            assert_eq!(result.output_power_w[1], 0.0);
            assert_eq!(result.mutual_coherence_w.im, 0.0);
        } else {
            assert!(result.output_power_w[1] > 0.0);
            // Ratio comparison remains sensitive even when dark power is far
            // below the absolute tolerances appropriate for a bright detector.
            close(
                result.output_power_w[1] / (native.power_w * (0.5 * phase).sin().powi(2)),
                1.0,
                3e-12,
            );
            assert!(result.mutual_coherence_w.im < 0.0);
        }
    }
}

#[test]
fn g1_spatial_fringe_average_matches_independent_square_pupil_integrals() {
    let camera =
        base(Aperture::try_regular_polygon(RADIUS, 4, core::f64::consts::FRAC_PI_4).unwrap())
            .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
            .unwrap();
    let mean = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let centers_q = [[-0.8, 0.2], [0.8, -0.3]];
    let centers = centers_q.map(|center| center.map(|q| q * SCALE));
    let low = [-2.0, -1.2];
    let high = [2.3, 1.5];
    let footprint = pixel(low, high);
    let distance = 50.0;
    let powers = [0.35_f64, 0.55];
    let path = 0.11 * LAMBDA;
    let phase = 0.23;
    let detector = receiver(centers, distance, powers, path, phase);
    let measured = query
        .integrate_interferometer_pixel(detector, footprint, PIXEL_LIMITS, || false)
        .unwrap();

    // Analytic Fourier transform of an axis-aligned square pupil, integrated
    // with the independent composite Simpson rule. No production field,
    // interferometer or pixel quadrature participates in this reference.
    let shape = |arm: usize, axis: usize, q: f64| {
        sinc_reference((q + centers_q[arm][axis]) * core::f64::consts::FRAC_1_SQRT_2)
    };
    let norm = 1.0 / (2.0 * core::f64::consts::PI.powi(2));
    let mut arms = powers.map(|power| norm * power);
    let mut cross = C64::from_re(norm * (powers[0] * powers[1]).sqrt());
    for axis in 0..2 {
        for (arm, power) in arms.iter_mut().enumerate() {
            *power *= integrate_reference(low[axis], high[axis], |q| shape(arm, axis, q).powi(2));
        }
        let slope = core::f64::consts::TAU / LAMBDA
            * distance
            * (centers[1][axis] - centers[0][axis])
            * SCALE;
        let real = integrate_reference(low[axis], high[axis], |q| {
            shape(0, axis, q) * shape(1, axis, q) * (slope * q).cos()
        });
        let imaginary = integrate_reference(low[axis], high[axis], |q| {
            -shape(0, axis, q) * shape(1, axis, q) * (slope * q).sin()
        });
        cross = cross * C64::new(real, imaginary);
    }
    let radius = centers.map(|c| c[0] * c[0] + c[1] * c[1]);
    let center_phase =
        core::f64::consts::TAU / LAMBDA * (0.5 * distance * (radius[1] - radius[0]) + path) + phase;
    cross = cross * C64::new(center_phase.cos(), -center_phase.sin());
    let outputs = [
        0.5 * (arms[0] + arms[1]) + cross.re,
        0.5 * (arms[0] + arms[1]) - cross.re,
    ];
    for arm in 0..2 {
        close(measured.arm_power_w[arm], arms[arm], 3e-10);
        close(measured.output_power_w[arm], outputs[arm], 3e-10);
    }
    close(measured.mutual_coherence_w.re, cross.re, 2e-8);
    close(measured.mutual_coherence_w.im, cross.im, 2e-8);
    close(
        measured.fringe_visibility,
        2.0 * cross.abs() / (arms[0] + arms[1]),
        2e-8,
    );
    assert!(measured.fringe_visibility < 0.1);
    assert!(
        query
            .observe_interferometer(detector, 2, || false)
            .unwrap()
            .fringe_visibility
            > 0.9
    );
    let pupil_only_x = (2.0 * (high[0] - low[0])).ceil() as usize;
    assert!(
        measured.subdivisions[0] > 2 * pupil_only_x,
        "relative sensor chirp must set spatial work"
    );
    assert_eq!(measured.pupil_evaluations, 2 * measured.quadrature_samples);

    let refined = query
        .integrate_interferometer_pixel(
            detector,
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_phase_span_rad: 0.5,
                ..PIXEL_LIMITS
            },
            || false,
        )
        .unwrap();
    assert!(refined.quadrature_samples > measured.quadrature_samples);
    close(refined.output_power_w[0], outputs[0], 1e-11);
    close(refined.output_power_w[1], outputs[1], 1e-11);
    close(refined.mutual_coherence_w.re, cross.re, 5e-10);
    close(refined.mutual_coherence_w.im, cross.im, 5e-10);
}

#[test]
fn g3_factor_powers_and_mutual_coherence_add_across_pixel_partitions() {
    let camera = focused();
    let left = [wave([-1.1, 0.2], C64::new(0.8, -0.3))];
    let right = [wave([0.9, -0.1], C64::new(-0.2, 0.7))];
    let factors = [
        PupilCoherenceFactor { waves: &left },
        PupilCoherenceFactor { waves: &right },
    ];
    let detector = receiver(
        [[-SCALE, 0.0], [SCALE, 0.3 * SCALE]],
        0.03,
        [0.4, 0.6],
        0.0,
        0.4,
    );
    let collect = |mean: &[CoherentPupilWave], factors: &[PupilCoherenceFactor<'_>], footprint| {
        camera
            .prepare_coherent_diffraction(
                LAMBDA_NM,
                CoherentPupilIllumination { mean, factors },
                LIMITS,
                || false,
            )
            .unwrap()
            .integrate_interferometer_pixel(detector, footprint, PIXEL_LIMITS, || false)
            .unwrap()
    };
    let footprint = pixel([-2.0, -1.0], [2.0, 1.0]);
    let mixed = collect(&[], &factors, footprint);
    let independent = [
        collect(&left, &[], footprint),
        collect(&right, &[], footprint),
    ];
    let partition = [
        collect(&[], &factors, pixel([-2.0, -1.0], [0.0, 1.0])),
        collect(&[], &factors, pixel([0.0, -1.0], [2.0, 1.0])),
    ];
    let additive = |parts: [CoherentInterferometerPixelPower; 2]| {
        for i in 0..2 {
            close(
                mixed.arm_power_w[i],
                parts[0].arm_power_w[i] + parts[1].arm_power_w[i],
                2e-10,
            );
            close(
                mixed.output_power_w[i],
                parts[0].output_power_w[i] + parts[1].output_power_w[i],
                2e-10,
            );
        }
        close(
            mixed.mutual_coherence_w.re,
            parts[0].mutual_coherence_w.re + parts[1].mutual_coherence_w.re,
            2e-10,
        );
        close(
            mixed.mutual_coherence_w.im,
            parts[0].mutual_coherence_w.im + parts[1].mutual_coherence_w.im,
            2e-10,
        );
    };
    additive(independent);
    additive(partition);
    close(
        mixed.output_power_w.iter().sum(),
        mixed.arm_power_w.iter().sum(),
        1e-14,
    );
    assert!(mixed.fringe_visibility < 1.0);
}

#[test]
fn g0_g4_complete_footprints_active_arm_work_and_cancelled_reuse() {
    let camera = focused();
    let mean = [wave([0.0; 2], C64::ONE)];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let footprint = pixel([-1.0; 2], [1.0; 2]);
    let centers = [[0.0; 2], [2.0 * SCALE, 0.0]];
    let balanced = receiver(centers, 0.03, [0.5; 2], 0.0, 0.0);
    let expected = query
        .integrate_interferometer_pixel(balanced, footprint, PIXEL_LIMITS, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        query.integrate_interferometer_pixel(
            balanced,
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_samples: expected.quadrature_samples - 1,
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
    assert_eq!(polls, 1);
    assert!(matches!(
        query.integrate_interferometer_pixel(
            balanced,
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_pupil_evaluations: expected.pupil_evaluations - 1,
                ..PIXEL_LIMITS
            },
            || {
                polls += 1;
                false
            },
        ),
        Err(CoherentDiffractionError::Limit {
            quantity: "pixel pupil evaluations",
            ..
        })
    ));
    assert_eq!(polls, 2);
    polls = 0;
    assert_eq!(
        query.integrate_interferometer_pixel(balanced, footprint, PIXEL_LIMITS, || {
            polls += 1;
            polls == 7
        }),
        Err(CoherentDiffractionError::Cancelled)
    );
    assert_eq!(
        query
            .integrate_interferometer_pixel(balanced, footprint, PIXEL_LIMITS, || false)
            .unwrap(),
        expected
    );

    let blocked = query
        .integrate_interferometer_pixel(
            receiver(centers, 0.03, [1.0, 0.0], 0.0, 0.7),
            footprint,
            PIXEL_LIMITS,
            || false,
        )
        .unwrap();
    let native = query
        .integrate_pixel(footprint, PIXEL_LIMITS, || false)
        .unwrap();
    assert_eq!(blocked.subdivisions, native.subdivisions);
    assert_eq!(blocked.pupil_evaluations, blocked.quadrature_samples);
    assert_eq!(blocked.mutual_coherence_w, C64::ZERO);
    assert_eq!(blocked.fringe_visibility, 0.0);
    close(blocked.arm_power_w[0], native.power_w, 1e-14);
    close(blocked.output_power_w[0], 0.5 * native.power_w, 1e-14);
    assert_eq!(blocked.output_power_w[0], blocked.output_power_w[1]);
    let dark = query
        .integrate_interferometer_pixel(
            receiver(centers, 0.03, [0.0; 2], 0.0, 0.0),
            footprint,
            CoherentDiffractionPixelLimits {
                maximum_samples: 0,
                maximum_pupil_evaluations: 0,
                ..PIXEL_LIMITS
            },
            || false,
        )
        .unwrap();
    assert_eq!(dark.subdivisions, [0; 2]);
    assert_eq!(dark.pupil_evaluations, 0);
    assert_eq!(dark.quadrature_samples, 0);
    assert_eq!(dark.output_power_w, [0.0; 2]);
    // Even an inactive declared arm must admit its complete translated footprint.
    let outside = receiver([[0.099, 0.0], [0.0; 2]], 0.03, [0.0, 1.0], 0.0, 0.0);
    let offset = CoherentDiffractionPixel::try_new([0.0; 2], [0.002, 0.001]).unwrap();
    assert_eq!(
        query.integrate_interferometer_pixel(outside, offset, PIXEL_LIMITS, || false),
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    );
}

#[test]
fn g3_finite_ports_survive_large_unnormalized_combiner_power() {
    let camera = focused();
    let one = core::f64::consts::PI * RADIUS * RADIUS / (LAMBDA * LAMBDA);
    // Each arm and quadrature port is representable; the two-arm sum is not.
    let mean = [wave(
        [0.0; 2],
        C64::from_re((1.5 * (f64::MAX / one)).sqrt()),
    )];
    let query = camera
        .prepare_coherent_diffraction(
            LAMBDA_NM,
            CoherentPupilIllumination {
                mean: &mean,
                factors: &[],
            },
            LIMITS,
            || false,
        )
        .unwrap();
    let detector = receiver(
        [[0.0; 2]; 2],
        0.03,
        [0.5; 2],
        0.0,
        core::f64::consts::FRAC_PI_2,
    );
    let point = query.observe_interferometer(detector, 2, || false).unwrap();
    assert!(
        point
            .arm_power_density_w
            .iter()
            .all(|p| p.is_finite() && *p > 0.7 * f64::MAX)
    );
    assert!(
        point
            .output_power_density_w
            .iter()
            .all(|p| p.is_finite() && *p > 0.7 * f64::MAX)
    );
    close(point.fringe_visibility, 1.0, 1e-14);
    let result = query
        .integrate_interferometer_pixel(detector, pixel([-0.1; 2], [0.1; 2]), PIXEL_LIMITS, || {
            false
        })
        .unwrap();
    assert!(
        result
            .output_power_w
            .iter()
            .all(|p| p.is_finite() && *p > 0.0)
    );
    close(
        result.output_power_w.iter().sum(),
        result.arm_power_w.iter().sum(),
        2e-14,
    );
    close(result.fringe_visibility, 1.0, 1e-14);
}
