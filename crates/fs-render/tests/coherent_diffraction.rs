#![cfg(feature = "chart-backends")]
//! G1 pupil-integral references and G3 coherence/phase invariants.

use fs_geom::{Point3, Vec3};
use fs_math::c64::C64;
use fs_render::camera::{
    Aperture, CameraError, CameraProjection, CoherentDiffractionError, CoherentDiffractionLimits,
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
