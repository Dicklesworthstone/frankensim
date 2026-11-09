//! G1/G2/G3 controls of the exact field, not a geometrical shadow model.
use super::*;

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:?}, expected {expected:?}, error {:e}, tolerance {tolerance:e}",
        (actual - expected).abs(),
    );
}

fn model(k: f64, direction: [f64; 3], boundary: IdealPlaneBoundary) -> HalfPlaneDiffraction {
    let medium = FluidAcousticState::try_new(1.3, 340.0).unwrap();
    HalfPlaneDiffraction::new(
        medium,
        k * medium.sound_speed_m_s() / core::f64::consts::TAU,
        direction,
        boundary,
    )
    .unwrap()
}

fn pressure(model: &HalfPlaneDiffraction, point: [f64; 3]) -> C64 {
    model
        .evaluate_harmonic_line(point, C64::ONE, || false)
        .unwrap()
        .pressure_pa
}

#[test]
fn fresnel_transition_matches_independent_landmarks() {
    // Independently evaluated using scipy.special.fresnel for F(t) and
    // scipy.special.wofz for exp(-i*t²)F(-t). No production recurrence or
    // subtraction of asymptotic C/S values generates the tail oracle.
    for (t, re, im) in [
        (0.0, 0.5, 0.0),
        (1e-6, 5.00000398942280433e-1, -3.98942280401299709e-7),
        (0.1, 5.40026808909984113e-1, -3.97608492894355667e-2),
        (0.5, 7.14776582125083770e-1, -1.81679513294920825e-1),
        (1.0, 9.84632105972107929e-1, -2.37073818320497065e-1),
        (1.499999, 1.16919476025600022, -4.82514523798185935e-2),
    ] {
        let expected = C64::new(re, im);
        close(
            transition_series(t, &mut || false).unwrap(),
            expected,
            4e-16,
        );
        close(
            transition_series(-t, &mut || false).unwrap(),
            C64::ONE - expected,
            4e-16,
        );
    }
    for (t, re, im) in [
        (1.5, 1.43826448361685788e-1, 1.01336017308885862e-1),
        (2.0, 1.07023941538385070e-1, 8.56229479358879275e-2),
        (4.0, 5.12630554492778875e-2, 4.81894852510921928e-2),
        (8.0, 2.51239550619161527e-2, 2.47347183136743427e-2),
        (20.0, 9.98597691564088313e-3, 9.96104360744638588e-3),
        (100.0, 1.99481112261318998e-3, 1.99461165148046883e-3),
    ] {
        close(
            scaled_tail(t, &mut || false).unwrap(),
            C64::new(re, im),
            2e-15,
        );
    }
    let t = 1.5;
    close(
        transition_series(t, &mut || false).unwrap(),
        C64::ONE - cis(Dd::from_f64(t * t)) * scaled_tail(t, &mut || false).unwrap(),
        2e-15,
    );
}

#[test]
fn full_field_matches_independent_sommerfeld_reference() {
    // A3 of https://pmc.ncbi.nlm.nih.gov/articles/PMC4844525/ evaluated with
    // scipy.special.fresnel, k=2.3/m and arrival azimuth 1.1 radians.
    let direction = [-1.1_f64.cos(), -1.1_f64.sin(), 0.0];
    let samples = [
        (
            [1.0, 0.3, 0.0],
            [
                C64::new(1.1171887056460457, -1.4362736962217908),
                C64::new(-1.0200577517201594, -0.5735649986662512),
            ],
        ),
        (
            [-0.5, 0.2, 0.0],
            [
                C64::new(1.0335748382196528, 0.2120443021143807),
                C64::new(1.040490788839605, -0.4374338685520189),
            ],
        ),
        (
            [0.4, -0.2, 0.0],
            [
                C64::new(-0.0017640426844132764, 0.4245871806640792),
                C64::new(0.01639347505470328, 0.044834064642683774),
            ],
        ),
        (
            [0.0, -2.0, 0.0],
            [
                C64::new(0.13959782890513206, -0.3551434396735196),
                C64::new(0.029671806504970068, -0.20035309162847456),
            ],
        ),
        (
            [1e-8, 1e-8, 0.0],
            [
                C64::new(1.0001602871482642, -0.00016029758084159466),
                C64::new(4.0705960238307526e-5, -4.0726462071967145e-5),
            ],
        ),
    ];
    for (i, boundary) in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ]
    .into_iter()
    .enumerate()
    {
        let wave = model(2.3, direction, boundary);
        for (point, expected) in samples {
            close(pressure(&wave, point), expected[i], 4e-15);
        }
    }
}

#[test]
fn both_sheet_faces_enforce_the_physical_boundary() {
    for direction in [[-0.48, -0.64, 0.6], [-0.48, 0.64, 0.6]] {
        for boundary in [
            IdealPlaneBoundary::Rigid,
            IdealPlaneBoundary::PressureRelease,
        ] {
            let wave = model(7.0, direction, boundary);
            for x in [1e-8, 0.02, 0.3, 5.0, 1000.0] {
                for y in [0.0, -0.0] {
                    let field = wave
                        .evaluate([x, y, 0.12], C64::new(0.6, -0.8), || false)
                        .unwrap();
                    match boundary {
                        IdealPlaneBoundary::Rigid => {
                            assert_eq!(field.particle_velocity_m_s[1], C64::ZERO);
                            assert_eq!(field.active_intensity_w_m2[1], 0.0);
                        }
                        IdealPlaneBoundary::PressureRelease => {
                            assert_eq!(field.pressure_pa, C64::ZERO);
                            assert_eq!(field.particle_velocity_m_s[0], C64::ZERO);
                            assert_eq!(field.particle_velocity_m_s[2], C64::ZERO);
                        }
                    }
                }
            }
            // There is no screen at x<0: no artificial branch discontinuity.
            close(
                pressure(&wave, [-0.4, 0.0, 0.1]),
                pressure(&wave, [-0.4, -0.0, 0.1]),
                0.0,
            );
            close(
                pressure(&wave, [-0.4, 1e-10, 0.1]),
                pressure(&wave, [-0.4, -1e-10, 0.1]),
                3e-9,
            );
        }
    }
}

#[test]
fn analytic_velocity_matches_independent_pressure_derivatives() {
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let wave = model(11.0, [-0.48, -0.64, 0.6], boundary);
        let rho_omega = wave.acoustic.density_kg_m3() * core::f64::consts::TAU * wave.frequency_hz;
        for point in [
            [0.3, 0.2, 0.1],
            [0.4, -0.2, -0.3],
            [-0.15, -0.25, 0.2],
            [1e-4, 3e-4, 0.01],
        ] {
            let field = wave.evaluate(point, C64::ONE, || false).unwrap();
            let h = 2e-5 * det::hypot(point[0], point[1]).min(1.0);
            for axis in 0..3 {
                let shifted = |shift: f64| {
                    let mut p = point;
                    p[axis] += shift;
                    pressure(&wave, p)
                };
                let grad = (shifted(-2.0 * h) - shifted(-h).scale(8.0) + shifted(h).scale(8.0)
                    - shifted(2.0 * h))
                .scale(1.0 / (12.0 * h));
                let expected = C64::new(grad.im, -grad.re).scale(1.0 / rho_omega);
                close(
                    field.particle_velocity_m_s[axis],
                    expected,
                    3e-8 * expected.abs().max(1e-3),
                );
            }
        }
    }
}

#[test]
fn full_wave_field_satisfies_helmholtz_without_paraxial_assumptions() {
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let wave = model(8.0, [-0.48, -0.64, 0.6], boundary);
        for point in [[0.3, 0.2, 0.1], [0.4, -0.2, -0.3], [-0.15, -0.25, 0.2]] {
            let p = pressure(&wave, point);
            let mut laplacian = C64::ZERO;
            let h = 1e-4;
            for axis in 0..3 {
                let mut lo = point;
                let mut hi = point;
                lo[axis] -= h;
                hi[axis] += h;
                laplacian = laplacian
                    + (pressure(&wave, lo) + pressure(&wave, hi) - p.scale(2.0))
                        .scale(1.0 / (h * h));
            }
            close(laplacian + p.scale(wave.k * wave.k), C64::ZERO, 2e-5);
        }
    }
}

#[test]
fn shadow_and_transition_fields_have_the_correct_far_limits() {
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let wave = model(1.0, [0.0, -1.0, 0.0], boundary);
        let mut previous = f64::INFINITY;
        for radius in [100.0, 10000.0, 1000000.0] {
            // On the geometrical shadow boundary, the incident transition is
            // exactly 1/2; the separate reflected edge field tends to zero.
            let transition = pressure(&wave, [0.0, -radius, 0.0]) * cis(Dd::from_f64(-radius));
            let error = (transition - C64::from_re(0.5)).abs();
            assert!(error < previous * 0.11);
            previous = error;
            let shadow = pressure(&wave, [0.8 * radius, -0.6 * radius, 0.0]);
            assert!(shadow.abs() > 0.0 && shadow.abs() < 1.0 / det::sqrt(radius));
            let lit = pressure(&wave, [-0.8 * radius, -0.6 * radius, 0.0]);
            close(
                lit,
                cis(Dd::from_f64(0.6 * radius)),
                1.0 / det::sqrt(radius),
            );
        }
        assert!(previous < 3e-4);
    }
}

#[test]
fn edge_phase_material_response_and_existing_observation_line_agree() {
    let wave = model(7.0, [-0.48, -0.64, 0.6], IdealPlaneBoundary::Rigid);
    let point = [0.1, -0.2, 0.3];
    let amplitude = C64::new(0.4, -0.7);
    let field = wave.evaluate(point, amplitude, || false).unwrap();
    let line = wave
        .evaluate_harmonic_line(point, amplitude, || false)
        .unwrap();
    assert_eq!(line.frequency_hz, wave.frequency_hz);
    assert_eq!(line.pressure_pa, field.pressure_pa);
    assert_eq!(line.particle_velocity_m_s, field.particle_velocity_m_s);
    let delta_z = 0.23;
    let translated = wave
        .evaluate([point[0], point[1], point[2] + delta_z], amplitude, || {
            false
        })
        .unwrap();
    let phase = cis(Dd::from_f64(wave.k * wave.direction[2] * delta_z));
    close(translated.pressure_pa, phase * field.pressure_pa, 3e-15);
    for axis in 0..3 {
        close(
            translated.particle_velocity_m_s[axis],
            phase * field.particle_velocity_m_s[axis],
            3e-17,
        );
        assert!(
            (translated.active_intensity_w_m2[axis] - field.active_intensity_w_m2[axis]).abs()
                < 3e-17
        );
    }
    let dense = HalfPlaneDiffraction::new(
        FluidAcousticState::try_new(2.6, 340.0).unwrap(),
        wave.frequency_hz,
        wave.direction,
        wave.boundary,
    )
    .unwrap();
    let dense_field = dense.evaluate(point, amplitude, || false).unwrap();
    assert_eq!(dense_field.pressure_pa, field.pressure_pa);
    for axis in 0..3 {
        close(
            dense_field.particle_velocity_m_s[axis].scale(2.0),
            field.particle_velocity_m_s[axis],
            0.0,
        );
    }
}

#[test]
fn edge_domain_and_cancellation_refuse_without_partial_fields() {
    let medium = FluidAcousticState::try_new(1.3, 340.0).unwrap();
    for direction in [
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, -2.0, 0.0],
        [f64::NAN, -1.0, 0.0],
    ] {
        assert!(
            HalfPlaneDiffraction::new(medium, 1000.0, direction, IdealPlaneBoundary::Rigid)
                .is_err()
        );
    }
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            HalfPlaneDiffraction::new(
                medium,
                frequency,
                [0.0, -1.0, 0.0],
                IdealPlaneBoundary::Rigid
            )
            .is_err()
        );
    }
    let wave = model(2.0, [0.0, -1.0, 0.0], IdealPlaneBoundary::Rigid);
    assert!(wave.evaluate([0.0, 0.0, 0.2], C64::ONE, || false).is_err());
    assert!(
        wave.evaluate([1.0, f64::NAN, 0.0], C64::ONE, || false)
            .is_err()
    );
    assert!(
        wave.evaluate([1.0, 0.1, 0.0], C64::new(f64::NAN, 0.0), || false)
            .is_err()
    );
    for after in [0, 1, 3] {
        let mut polls = 0;
        let result = wave.evaluate([1.0, 0.1, 0.0], C64::ONE, || {
            polls += 1;
            polls > after
        });
        assert!(matches!(result, Err(AeroacError::Cancelled { .. })));
    }
}
