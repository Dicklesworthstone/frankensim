//! Polarization correlations through original grouped and two-face sources.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureApertureJonesProjection, AngularSpectrumExposureApertureJonesSpectrum,
    AngularSpectrumExposureApertureProjection,
};

#[test]
fn g1_dispersive_analyzers_follow_original_frequency_banks_on_both_faces() {
    let medium = material(1.0, 1.0);
    let stacks = [1.0, 2.0].map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let circular = [C64::ONE, C64::new(0.0, 1.0)];
    let lines = [
        constant(&stacks[0], circular),
        constant(&stacks[1], circular),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    assert_eq!(bank.angular_frequencies_rad_s(), &[1.0, 2.0]);
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let illumination: Vec<_> = lines
        .iter()
        .map(|line| {
            line.prepare_illumination(
                &[[C64::ZERO; 2]; 4],
                AngularSpectrumIlluminationLimits {
                    maximum_layer_mode_pairs: 4096,
                    maximum_workspace_bytes: 8 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
        })
        .collect();
    let two_face =
        AngularSpectrumExposure::try_new_illumination(&lines[0], &illumination, limits(), || false)
            .unwrap();
    let two_face_factors =
        [
            AngularSpectrumExposure::try_new_illumination(
                &lines[0],
                &illumination,
                limits(),
                || false,
            )
            .unwrap(),
        ];
    let two_face_ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &two_face_factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let r = core::f64::consts::FRAC_1_SQRT_2;
    let quarter = [C64::from_re(r), C64::new(0.0, r), C64::ZERO];
    let half = [C64::from_re(r), C64::from_re(-r), C64::ZERO];
    let samples = [
        (
            bank.angular_frequencies_rad_s()[0],
            AngularSpectrumExposureApertureJonesProjection::try_new(quarter, quarter).unwrap(),
        ),
        (
            bank.angular_frequencies_rad_s()[1],
            AngularSpectrumExposureApertureJonesProjection::try_new(half, half).unwrap(),
        ),
    ];
    let spectrum =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples, || false).unwrap();
    let area = AngularSpectrumAperture::try_new([0.2, -0.7], [0.4, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.17, core::f64::consts::TAU).unwrap();
    for (side, face) in [
        (Side::Incident, Face::TotalLeft),
        (Side::Transmitted, Face::OutgoingRight),
    ] {
        for actual in [
            bank.measure_aperture_spectral_projection(
                side,
                0.19,
                area,
                window,
                spectrum,
                aperture_limits(),
                || false,
            )
            .unwrap(),
            ensemble
                .measure_aperture_spectral_projection(
                    side,
                    0.19,
                    area,
                    window,
                    spectrum,
                    aperture_limits(),
                    || false,
                )
                .unwrap(),
            two_face
                .measure_illumination_aperture_spectral_projection(
                    face,
                    0.19,
                    area,
                    window,
                    spectrum,
                    aperture_limits(),
                    || false,
                )
                .unwrap(),
            two_face_ensemble
                .measure_illumination_aperture_spectral_projection(
                    face,
                    0.19,
                    area,
                    window,
                    spectrum,
                    aperture_limits(),
                    || false,
                )
                .unwrap(),
        ] {
            // One unit from the quarter-wave circular selection plus one half
            // from the second frequency's half-wave response; no beat survives.
            near(actual.electric_mean_squared_v2_m2, 1.5, 4e-14);
            near(actual.magnetic_mean_squared_a2_m2, 1.5, 4e-14);
            near(actual.electric_magnetic_mean_product_w_m2, 0.0, 4e-14);
        }
    }
    let incomplete =
        AngularSpectrumExposureApertureJonesSpectrum::try_new(&samples[..1], || false).unwrap();
    assert!(
        bank.measure_aperture_spectral_projection(
            Side::Incident,
            0.0,
            area,
            window,
            incomplete,
            aperture_limits(),
            || false,
        )
        .is_err()
    );
    assert!(
        two_face_ensemble
            .measure_illumination_aperture_spectral_projection(
                Face::TotalLeft,
                0.0,
                area,
                window,
                incomplete,
                aperture_limits(),
                || false,
            )
            .is_err()
    );
}

#[test]
fn g1_circular_analyzers_distinguish_native_maxwell_handedness() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::new(0.0, 1.0)])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors,
        &[AngularSpectrumFactorStatistics::ProperComplex],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let area = AngularSpectrumAperture::try_new([0.2, -0.7], [0.4, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.17, core::f64::consts::TAU).unwrap();
    let r = core::f64::consts::FRAC_1_SQRT_2;
    for handedness in [-1.0, 1.0] {
        let weights = [C64::from_re(r), C64::new(0.0, handedness * r), C64::ZERO];
        let analyzer =
            AngularSpectrumExposureApertureJonesProjection::try_new(weights, weights).unwrap();
        for side in [Side::Incident, Side::Transmitted] {
            for actual in [
                bank.measure_aperture_projection(
                    side,
                    0.19,
                    area,
                    window,
                    analyzer,
                    aperture_limits(),
                    || false,
                )
                .unwrap(),
                ensemble
                    .measure_aperture_projection(
                        side,
                        0.19,
                        area,
                        window,
                        analyzer,
                        aperture_limits(),
                        || false,
                    )
                    .unwrap(),
            ] {
                let expected = if handedness == 1.0 { 1.0 } else { 0.0 };
                near(actual.electric_mean_squared_v2_m2, expected, 3e-14);
                near(actual.magnetic_mean_squared_a2_m2, expected, 3e-14);
                near(actual.electric_magnetic_mean_product_w_m2, 0.0, 3e-14);
            }
        }
    }
}

#[test]
fn g1_grouped_original_sources_retain_crossed_analyzer_residuals() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let tiny = 1e-20;
    let lines = [
        constant(&stack, [C64::ONE, C64::ONE]),
        constant(&stack, [C64::from_re(tiny), C64::ZERO]),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let area = AngularSpectrumAperture::try_new([0.2, -0.7], [0.4, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.17, core::f64::consts::TAU).unwrap();
    let analyzer =
        AngularSpectrumExposureApertureProjection::try_new([1.0, -1.0, 0.0], [-1.0, -1.0, 0.0])
            .unwrap();
    for side in [Side::Incident, Side::Transmitted] {
        let direct = bank
            .measure_aperture_projection(
                side,
                0.0,
                area,
                window,
                analyzer,
                aperture_limits(),
                || false,
            )
            .unwrap();
        let expected = 0.5 * tiny * tiny;
        near(direct.electric_mean_squared_v2_m2 / expected, 1.0, 3e-13);
        near(direct.magnetic_mean_squared_a2_m2 / expected, 1.0, 3e-13);
        near(
            direct.electric_magnetic_mean_product_w_m2 / expected,
            -1.0,
            3e-13,
        );
        let moments = bank
            .measure_aperture_moments(side, 0.0, area, window, aperture_limits(), || false)
            .unwrap();
        near(moments.electric_second_moment_v2_m2()[0][1], 0.5, 2e-14);
        near(moments.magnetic_second_moment_a2_m2()[0][1], -0.5, 2e-14);
        near(
            moments.electric_magnetic_second_moment_w_m2()[0][0],
            -0.5,
            2e-14,
        );
        let original = bank
            .measure_aperture(side, 0.0, area, window, aperture_limits(), || false)
            .unwrap();
        for (a, b) in values(moments.reading().mean)
            .into_iter()
            .zip(values(original.mean))
        {
            near(a, b, 3e-14);
        }
        near(
            moments.reading().normal_power_w,
            original.normal_power_w,
            3e-14,
        );
    }
}

#[test]
fn g1_two_face_mixed_source_correlations_match_real_phase_ensemble_quadrature() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let left = constant(&stack, [C64::new(0.4, -0.3), C64::new(0.2, 0.7)]);
    let lines = [left
        .prepare_illumination(
            &[[C64::new(-0.2, 0.6), C64::new(0.3, 0.1)]; 4],
            AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 4096,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap()];
    let mean =
        AngularSpectrumExposure::try_new_illumination(&left, &lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new_illumination(&left, &lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new_illumination(&left, &lines, limits(), || false).unwrap(),
    ];
    let laws = [
        AngularSpectrumFactorStatistics::Real,
        AngularSpectrumFactorStatistics::ProperComplex,
    ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left,
        Some(&mean),
        &factors,
        &laws,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let area = AngularSpectrumAperture::try_new([-0.2, 0.7], [0.83, 0.41]).unwrap();
    let window =
        HarmonicWindow::with_profile(0.11, 0.37, HarmonicWindowProfile::Triangular).unwrap();
    let analyzer =
        AngularSpectrumExposureApertureProjection::try_new([0.6, -0.8, 0.0], [0.8, 0.6, 0.0])
            .unwrap();
    for side in [
        Face::IncomingLeft,
        Face::OutgoingLeft,
        Face::TotalLeft,
        Face::IncomingRight,
        Face::OutgoingRight,
        Face::TotalRight,
    ] {
        let native = lines[0].field(side, 0.17, || false).unwrap();
        let e = native.electric_v_per_m()[0];
        let h = native.magnetic_a_per_m()[0];
        let amplitudes = [e[0], e[1], e[2], h[0], h[1], h[2]];
        let expected = phase_reference(amplitudes, window);
        let moments = ensemble
            .measure_illumination_aperture_moments(
                side,
                0.17,
                area,
                window,
                aperture_limits(),
                || false,
            )
            .unwrap();
        for (matrix, row, column) in [
            (moments.electric_second_moment_v2_m2(), 0, 0),
            (moments.magnetic_second_moment_a2_m2(), 3, 3),
            (moments.electric_magnetic_second_moment_w_m2(), 0, 3),
        ] {
            for i in 0..3 {
                for j in 0..3 {
                    near(matrix[i][j], expected[row + i][column + j], 4e-12);
                }
            }
        }
        let direct = ensemble
            .measure_illumination_aperture_projection(
                side,
                0.17,
                area,
                window,
                analyzer,
                aperture_limits(),
                || false,
            )
            .unwrap();
        let contracted = moments.project(analyzer).unwrap();
        near(
            direct.electric_mean_squared_v2_m2,
            contracted.electric_mean_squared_v2_m2,
            3e-14,
        );
        near(
            direct.magnetic_mean_squared_a2_m2,
            contracted.magnetic_mean_squared_a2_m2,
            3e-14,
        );
        near(
            direct.electric_magnetic_mean_product_w_m2,
            contracted.electric_magnetic_mean_product_w_m2,
            3e-14,
        );
        let a = [C64::new(0.3, -0.7), C64::new(-0.4, 0.2), C64::new(0.1, 0.3)];
        let b = [C64::new(-0.2, 0.1), C64::new(0.6, -0.4), C64::new(0.5, 0.2)];
        let complex = AngularSpectrumExposureApertureJonesProjection::try_new(a, b).unwrap();
        let projected_e = (0..3).fold(C64::ZERO, |sum, i| sum + a[i].conj() * e[i]);
        let projected_h = (0..3).fold(C64::ZERO, |sum, i| sum + b[i].conj() * h[i]);
        let expected = phase_reference(
            [
                projected_e,
                C64::ZERO,
                C64::ZERO,
                projected_h,
                C64::ZERO,
                C64::ZERO,
            ],
            window,
        );
        let actual = ensemble
            .measure_illumination_aperture_projection(
                side,
                0.17,
                area,
                window,
                complex,
                aperture_limits(),
                || false,
            )
            .unwrap();
        near(actual.electric_mean_squared_v2_m2, expected[0][0], 4e-12);
        near(actual.magnetic_mean_squared_a2_m2, expected[3][3], 4e-12);
        near(
            actual.electric_magnetic_mean_product_w_m2,
            expected[0][3],
            4e-12,
        );
    }
}

// The source is mean + real z1 + proper z2, with all three carrying the same
// original phasor. Enumerate their actual joint second moments, then integrate
// real signals directly; no production covariance or harmonic kernel is used.
fn phase_reference(amplitudes: [C64; 6], window: HarmonicWindow) -> [[f64; 6]; 6] {
    let steps = 512;
    let mut sums = [[0.0; 6]; 6];
    for real in [-1.0, 1.0] {
        for proper in [
            C64::ONE,
            C64::from_re(-1.0),
            C64::new(0.0, 1.0),
            C64::new(0.0, -1.0),
        ] {
            let coefficient = C64::from_re(1.0 + real) + proper;
            for i in 0..=steps {
                let u = i as f64 / steps as f64;
                let t = window.start_s() + window.duration_s() * u;
                let (s, c) = (2.0 * t).sin_cos();
                let values = amplitudes.map(|a| (a * coefficient * C64::new(c, s)).re);
                let weight = if i == 0 || i == steps {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let weight = weight * 2.0 * (1.0 - (2.0 * u - 1.0).abs());
                for a in 0..6 {
                    for b in 0..6 {
                        sums[a][b] += weight * values[a] * values[b];
                    }
                }
            }
        }
    }
    sums.map(|row| row.map(|v| v / (8.0 * 3.0 * steps as f64)))
}
