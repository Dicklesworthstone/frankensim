//! Full Jones/coherency behavior on the original physical pupil and pixel rule.
use super::*;
use fs_render::camera::{
    JonesPupilCoherenceFactor, CoherentJonesPupilIllumination, CoherentJonesPupilWave,
    PolarizedCoherentDiffractionSample,
};
use fs_render::polarization::{Mueller, Stokes};

fn jones(q: [f64; 2], amplitude_sqrt_w: [C64; 2]) -> CoherentJonesPupilWave {
    CoherentJonesPupilWave {
        tangent_direction: q.map(|x| x * SCALE),
        amplitude_sqrt_w,
    }
}

fn observe(
    camera: &PhysicalCamera,
    mean: &[CoherentJonesPupilWave],
    factors: &[JonesPupilCoherenceFactor<'_>],
    at: [f64; 2],
) -> PolarizedCoherentDiffractionSample {
    camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean, factors }, LIMITS, || false,
    ).unwrap().observe(at, || false).unwrap()
}

fn components(s: Stokes) -> [f64; 4] {
    [s.i, s.q, s.u, s.v]
}

fn stokes_close(actual: Stokes, expected: Stokes, tolerance: f64) {
    for (a, b) in components(actual).into_iter().zip(components(expected)) {
        assert!((a - b).abs() <= tolerance * expected.i.max(1e-280),
            "{actual:?} != {expected:?}");
    }
    assert!(Stokes::try_new(actual.i, actual.q, actual.u, actual.v).is_some());
}

#[test]
fn g1_one_axis_jones_preserves_scalar_point_and_pixel_results_for_original_pupils() {
    let annular = base(Aperture::try_annular(RADIUS, 0.4 * RADIUS).unwrap())
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap()).unwrap();
    let polygon = base(Aperture::try_regular_polygon(RADIUS, 4, 0.3).unwrap())
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap()).unwrap();
    let scalar = [
        wave([-0.3, 0.2], C64::new(0.7, -0.2)),
        wave([0.6, -0.1], C64::new(-0.1, 0.4)),
    ];
    let vector = scalar.map(|w| CoherentJonesPupilWave {
        tangent_direction: w.tangent_direction,
        amplitude_sqrt_w: [w.amplitude_sqrt_w, C64::ZERO],
    });
    for camera in [focused(), defocused(2.4), annular, polygon] {
        let original = camera.prepare_coherent_diffraction(
            LAMBDA_NM, CoherentPupilIllumination { mean: &scalar, factors: &[] },
            LIMITS, || false,
        ).unwrap();
        let prepared = camera.prepare_polarized_coherent_diffraction(
            LAMBDA_NM, CoherentJonesPupilIllumination { mean: &vector, factors: &[] },
            LIMITS, || false,
        ).unwrap();
        for at in [[0.0, 0.0], [0.7 * SCALE, -0.4 * SCALE]] {
            let expected = original.observe(at, || false).unwrap();
            let actual = prepared.observe(at, || false).unwrap();
            assert_eq!(actual.mean_field[0], expected.mean_field);
            assert_eq!(actual.mean_field[1], C64::ZERO);
            assert_eq!(actual.coherency_density_w.axis_powers(), [expected.power_density_w, 0.0]);
            assert_eq!(actual.coherency_density_w.stokes(),
                Stokes { i: expected.power_density_w, q: expected.power_density_w, u: 0.0, v: 0.0 });
        }
        let footprint = pixel([-0.7, -0.4], [0.8, 0.6]);
        let expected = original.integrate_pixel(footprint, PIXEL_LIMITS, || false).unwrap();
        let actual = prepared.integrate_pixel(footprint, PIXEL_LIMITS, || false).unwrap();
        assert_eq!(actual.coherency_w.axis_powers(), [expected.power_w, 0.0]);
        assert_eq!(actual.subdivisions, expected.subdivisions);
        assert_eq!(actual.quadrature_samples, expected.quadrature_samples);
        assert_eq!(actual.pupil_evaluations, expected.pupil_evaluations);
    }
}

#[test]
fn g1_orthogonal_fields_keep_total_intensity_and_interfere_after_explicit_projection() {
    let camera = focused();
    let peak = intensity(&camera, &[wave([0.0; 2], C64::ONE)], &[], [0.0; 2]);
    let x = jones([0.0; 2], [C64::ONE, C64::ZERO]);
    for phase in [C64::ONE, C64::new(0.0, 1.0), C64::new(-1.0, 0.0), C64::new(0.0, -1.0)] {
        let y = jones([0.0; 2], [C64::ZERO, phase]);
        let coherent = observe(&camera, &[x, y], &[], [0.0; 2]).coherency_density_w.stokes();
        stokes_close(coherent, Stokes { i: 2.0 * peak, q: 0.0,
            u: 2.0 * peak * phase.re, v: -2.0 * peak * phase.im }, 2e-13);
        let analysed = Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4).apply(coherent).i;
        close(analysed, peak * (1.0 + phase.re), 2e-12);
        let a = [x]; let b = [y];
        let independent = [
            JonesPupilCoherenceFactor { waves: &a },
            JonesPupilCoherenceFactor { waves: &b },
        ];
        let mixed = observe(&camera, &[], &independent, [0.0; 2]).coherency_density_w.stokes();
        assert_eq!(mixed, Stokes::unpolarized(2.0 * peak));
        close(Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4).apply(mixed).i, peak, 2e-13);
    }
    let opposite = jones([0.0; 2], [C64::new(-1.0, 0.0), C64::ZERO]);
    assert_eq!(observe(&camera, &[x, opposite], &[], [0.0; 2]).coherency_density_w.stokes(),
        Stokes::unpolarized(0.0));
}

fn combine(a: C64, b: C64, x: f64, y: f64) -> C64 {
    C64::new(x * a.re + y * b.re, x * a.im + y * b.im)
}

fn factor_mix(a: CoherentJonesPupilWave, b: CoherentJonesPupilWave, x: f64, y: f64) -> CoherentJonesPupilWave {
    assert_eq!(a.tangent_direction, b.tangent_direction);
    CoherentJonesPupilWave { tangent_direction: a.tangent_direction,
        amplitude_sqrt_w: std::array::from_fn(|i| combine(a.amplitude_sqrt_w[i], b.amplitude_sqrt_w[i], x, y)) }
}

fn rotate(w: CoherentJonesPupilWave) -> CoherentJonesPupilWave {
    let [x, y] = w.amplitude_sqrt_w;
    CoherentJonesPupilWave { amplitude_sqrt_w: [combine(x, y, 0.6, 0.8), combine(x, y, -0.8, 0.6)], ..w }
}

#[test]
fn g3_covariance_factorization_and_common_polarization_rotation_preserve_physical_results() {
    let camera = defocused(1.7);
    let mean = [jones([0.1, -0.2], [C64::new(0.2, 0.1), C64::new(-0.1, 0.3)])];
    let a = [
        jones([-0.5, 0.2], [C64::new(0.7, -0.2), C64::new(0.1, 0.4)]),
        jones([0.4, -0.3], [C64::new(-0.3, 0.1), C64::new(0.5, -0.2)]),
    ];
    let b = [
        jones([-0.5, 0.2], [C64::new(-0.4, 0.3), C64::new(0.2, 0.1)]),
        jones([0.4, -0.3], [C64::new(0.6, -0.1), C64::new(-0.3, 0.4)]),
    ];
    let mixed_a = std::array::from_fn::<_, 2, _>(|i| factor_mix(a[i], b[i], 0.6, 0.8));
    let mixed_b = std::array::from_fn::<_, 2, _>(|i| factor_mix(a[i], b[i], -0.8, 0.6));
    let original = [JonesPupilCoherenceFactor { waves: &a }, JonesPupilCoherenceFactor { waves: &b }];
    let mixed = [JonesPupilCoherenceFactor { waves: &mixed_a }, JonesPupilCoherenceFactor { waves: &mixed_b }];
    let rotated_a = a.map(rotate); let rotated_b = b.map(rotate); let rotated_mean = mean.map(rotate);
    let rotated = [JonesPupilCoherenceFactor { waves: &rotated_a }, JonesPupilCoherenceFactor { waves: &rotated_b }];
    for at in [[0.0; 2], [0.8 * SCALE, -0.6 * SCALE]] {
        let expected = observe(&camera, &mean, &original, at);
        let actual = observe(&camera, &mean, &mixed, at);
        stokes_close(actual.coherency_density_w.stokes(), expected.coherency_density_w.stokes(), 3e-12);
        stokes_close(actual.fluctuation_coherency_density_w.stokes(),
            expected.fluctuation_coherency_density_w.stokes(), 3e-12);
        assert_eq!(actual.mean_coherency_density_w, expected.mean_coherency_density_w);
        let actual = observe(&camera, &rotated_mean, &rotated, at).coherency_density_w.stokes();
        stokes_close(actual, expected.coherency_density_w.stokes().rotate_axes(0.8_f64.atan2(0.6)), 3e-12);
    }
    let footprint = pixel([-0.4, -0.3], [0.7, 0.6]);
    let integrated = |mean, factors| camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean, factors }, LIMITS, || false,
    ).unwrap().integrate_pixel(footprint, PIXEL_LIMITS, || false).unwrap();
    let expected = integrated(&mean, &original);
    let actual = integrated(&mean, &mixed);
    stokes_close(actual.coherency_w.stokes(), expected.coherency_w.stokes(), 3e-12);
    stokes_close(actual.fluctuation_coherency_w.stokes(),
        expected.fluctuation_coherency_w.stokes(), 3e-12);
    assert_eq!(actual.mean_coherency_w, expected.mean_coherency_w);
    let actual = integrated(&rotated_mean, &rotated);
    stokes_close(actual.coherency_w.stokes(),
        expected.coherency_w.stokes().rotate_axes(0.8_f64.atan2(0.6)), 3e-12);
}

fn point_reference(camera: &PhysicalCamera, waves: &[CoherentJonesPupilWave], at: [f64; 2]) -> [f64; 4] {
    let mut field = [C64::ZERO; 2];
    for wave in waves {
        let h = camera.in_focus_diffraction_amplitude(
            std::array::from_fn(|i| at[i] - wave.tangent_direction[i]), LAMBDA_NM,
        ).unwrap();
        for i in 0..2 { field[i] = field[i] + h * wave.amplitude_sqrt_w[i]; }
    }
    let power = field.map(|e| e.re * e.re + e.im * e.im);
    let cross_re = field[0].re * field[1].re + field[0].im * field[1].im;
    let cross_im = field[0].im * field[1].re - field[0].re * field[1].im;
    [power[0] + power[1], power[0] - power[1], 2.0 * cross_re, 2.0 * cross_im]
}

#[test]
fn g1_positive_finite_pixel_integrates_coherency_and_resolves_spatial_depolarization() {
    let camera = focused();
    let waves = [
        jones([-1.4, 0.0], [C64::ONE, C64::ZERO]),
        jones([1.4, 0.0], [C64::ZERO, C64::new(0.0, 1.0)]),
    ];
    let prepared = camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &waves, factors: &[] }, LIMITS, || false,
    ).unwrap();
    let footprint = pixel([-3.0, -2.0], [3.0, 2.0]);
    let result = prepared.integrate_pixel(
        footprint, CoherentDiffractionPixelLimits { maximum_phase_span_rad: 0.25, ..PIXEL_LIMITS }, || false,
    ).unwrap();
    let expected = result.coherency_w.stokes();
    assert!(expected.i > 0.0 && expected.degree() < 0.999);
    assert_eq!(result.pupil_evaluations, 2 * result.quadrature_samples);
    let [x, y] = result.coherency_w.axis_powers();
    let cross = result.coherency_w.cross();
    assert!(cross.re.hypot(cross.im) <= x.sqrt() * y.sqrt() * (1.0 + 2e-13));
    let mut previous = f64::INFINITY;
    for n in [32, 64, 128] {
        let mut sum = [0.0; 4];
        for j in 0..n { for i in 0..n {
            let at = [(-3.0 + 6.0 * (i as f64 + 0.5) / n as f64) * SCALE,
                (-2.0 + 4.0 * (j as f64 + 0.5) / n as f64) * SCALE];
            let value = point_reference(&camera, &waves, at);
            for k in 0..4 { sum[k] += value[k]; }
        }}
        let reference = sum.map(|s| s * footprint.area_tangent() / (n * n) as f64);
        let error = reference.into_iter().zip(components(expected))
            .map(|(a, b)| (a - b).abs() / expected.i).fold(0.0, f64::max);
        assert!(error < previous * 0.45);
        previous = error;
    }
    assert!(previous < 1e-4);
}

#[test]
fn g1_dd_jones_cancellation_retains_weak_fields_and_individual_axis_powers() {
    let camera = focused();
    let single = [jones([0.0; 2], [C64::new(0.3, -0.2), C64::new(-0.4, 0.1)])];
    let strong = jones([0.0; 2], [C64::new(1e16, -1e16), C64::new(-1e16, -1e16)]);
    let opposite = CoherentJonesPupilWave { amplitude_sqrt_w: strong.amplitude_sqrt_w.map(|a| C64::new(-a.re, -a.im)), ..strong };
    let expected = observe(&camera, &single, &[], [0.0; 2]);
    let actual = observe(&camera, &[strong, single[0], opposite], &[], [0.0; 2]);
    stokes_close(actual.coherency_density_w.stokes(), expected.coherency_density_w.stokes(), 3e-12);
    let weak = observe(&camera, &[jones([0.0; 2], [C64::ONE, C64::new(0.0, 1e-120)])], &[], [0.0; 2]);
    let s = weak.coherency_density_w.stokes();
    let axes = weak.coherency_density_w.axis_powers();
    assert_eq!(s.i, s.q);
    assert!(axes[1] > 0.0);
    assert!((axes[1] / (axes[0] * 1e-240) - 1.0).abs() < 2e-12);
    assert!(s.v < 0.0);
    let weak_waves = [jones([0.0; 2], [C64::ONE, C64::new(0.0, 1e-120)])];
    let prepared = camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &weak_waves, factors: &[] }, LIMITS, || false,
    ).unwrap();
    let measured = prepared.integrate_pixel(pixel([-0.4; 2], [0.6; 2]), PIXEL_LIMITS, || false).unwrap();
    let axes = measured.coherency_w.axis_powers();
    assert!(axes[1] > 0.0);
    assert!((axes[1] / (axes[0] * 1e-240) - 1.0).abs() < 2e-12);
}

#[test]
fn g0_jones_admission_exact_work_limits_and_cancellation_remain_transactional() {
    let camera = focused();
    let waves = [jones([0.0; 2], [C64::ONE, C64::ZERO]), jones([0.3, -0.2], [C64::ZERO, C64::ONE])];
    let illumination = CoherentJonesPupilIllumination { mean: &waves, factors: &[] };
    assert!(matches!(camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, illumination, LIMITS, || true), Err(CoherentDiffractionError::Cancelled)));
    assert!(matches!(camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, illumination, CoherentDiffractionLimits { maximum_waves: 1, ..LIMITS }, || false),
        Err(CoherentDiffractionError::Limit { quantity: "pupil waves", .. })));
    let factors = [JonesPupilCoherenceFactor { waves: &waves }];
    assert!(matches!(camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &[], factors: &factors },
        CoherentDiffractionLimits { maximum_factors: 0, ..LIMITS }, || false),
        Err(CoherentDiffractionError::Limit { quantity: "source factors", .. })));
    let bad = [jones([0.0; 2], [C64::ONE, C64::new(f64::NAN, 0.0)])];
    assert!(matches!(camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &bad, factors: &[] }, LIMITS, || false),
        Err(CoherentDiffractionError::NonFiniteAmplitude)));
    let prepared = camera.prepare_polarized_coherent_diffraction(LAMBDA_NM, illumination, LIMITS, || false).unwrap();
    assert_eq!(prepared.wavelength_nm(), LAMBDA_NM);
    assert_eq!(prepared.wave_count(), 2);
    assert_eq!(prepared.factor_count(), 0);
    let point = prepared.observe([0.0; 2], || false).unwrap();
    let mut polls = 0;
    assert_eq!(prepared.observe([0.0; 2], || { polls += 1; polls == 3 }),
        Err(CoherentDiffractionError::Cancelled));
    assert_eq!(prepared.observe([0.0; 2], || false).unwrap(), point);
    assert_eq!(prepared.observe([0.11, 0.0], || false),
        Err(CoherentDiffractionError::OutsideParaxialDomain));
    let footprint = pixel([-0.4; 2], [0.6; 2]);
    let result = prepared.integrate_pixel(footprint, PIXEL_LIMITS, || false).unwrap();
    assert_eq!(result.pupil_evaluations, 2 * result.quadrature_samples);
    assert!(matches!(prepared.integrate_pixel(footprint,
        CoherentDiffractionPixelLimits { maximum_pupil_evaluations: result.pupil_evaluations - 1, ..PIXEL_LIMITS },
        || false), Err(CoherentDiffractionError::Limit { quantity: "pixel pupil evaluations", .. })));
    let mut polls = 0;
    assert_eq!(prepared.integrate_pixel(footprint, PIXEL_LIMITS, || { polls += 1; polls == 4 }),
        Err(CoherentDiffractionError::Cancelled));
    assert_eq!(prepared.integrate_pixel(footprint, PIXEL_LIMITS, || false).unwrap(), result);
    for a in [1e200, 1e-200] {
        let range = [jones([0.0; 2], [C64::new(a, 0.0), C64::ZERO])];
        let prepared = camera.prepare_polarized_coherent_diffraction(
            LAMBDA_NM, CoherentJonesPupilIllumination { mean: &range, factors: &[] }, LIMITS, || false,
        ).unwrap();
        assert_eq!(prepared.observe([0.0; 2], || false), Err(CoherentDiffractionError::Unrepresentable));
    }
    // The represented pupil transfer is nonzero but below 1/2: multiplying
    // it by the least positive authored field rounds the whole term to zero.
    let at = [0.075, 0.0];
    let h = camera.in_focus_diffraction_amplitude(at, LAMBDA_NM).unwrap();
    assert!(h.re != 0.0 || h.im != 0.0);
    assert!(h.re.abs() < 0.5 && h.im.abs() < 0.5);
    let tiny = [jones([0.0; 2], [C64::new(f64::from_bits(1), 0.0), C64::ZERO])];
    let underflow = camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &tiny, factors: &[] }, LIMITS, || false,
    ).unwrap();
    assert_eq!(underflow.observe(at, || false), Err(CoherentDiffractionError::Unrepresentable));
    let empty = camera.prepare_polarized_coherent_diffraction(
        LAMBDA_NM, CoherentJonesPupilIllumination { mean: &[], factors: &[] },
        CoherentDiffractionLimits { maximum_waves: 0, maximum_factors: 0 }, || false,
    ).unwrap();
    let dark = empty.integrate_pixel(footprint,
        CoherentDiffractionPixelLimits { maximum_samples: 0, maximum_pupil_evaluations: 0, ..PIXEL_LIMITS }, || false).unwrap();
    assert_eq!(dark.coherency_w.stokes(), Stokes::unpolarized(0.0));
    assert_eq!(dark.quadrature_samples, 0);
}
