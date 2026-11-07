#![cfg(feature = "chart-backends")]
//! G1 pupil-integral references and G3 coherence/phase invariants.

use fs_geom::{Point3, Vec3};
use fs_math::c64::C64;
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, CoherentDiffractionError, CoherentDiffractionLimits,
    CoherentDiffractionPixel, CoherentDiffractionPixelLimits, CoherentDiffractionPixelPower,
    CoherentPupilIllumination, CoherentPupilWave, DefocusedDiffraction, FocusedDiffraction,
    PhysicalCamera, PupilCoherenceFactor,
};

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
