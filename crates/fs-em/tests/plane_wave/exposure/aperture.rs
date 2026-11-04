//! Finite-area shutters through the original Maxwell and covariance owners.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposureApertureLimits, AngularSpectrumExposurePoint,
    AngularSpectrumFactorStatistics, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn aperture_limits() -> AngularSpectrumExposureApertureLimits {
    AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 4096,
        maximum_workspace_bytes: 8 << 20,
    }
}

fn values(point: AngularSpectrumExposurePoint) -> [f64; 9] {
    let e = point.electric_mean_squared_v2_m2;
    let h = point.magnetic_mean_squared_a2_m2;
    let p = point.poynting_w_m2;
    [e[0], e[1], e[2], h[0], h[1], h[2], p[0], p[1], p[2]]
}

#[test]
fn g1_joint_shutter_matches_original_cell_power_and_detector_partitions() {
    let medium = material(1.0, 1.0);
    let film = material(2.25, 1.0);
    // Both frequencies have propagating incident modes; the low-index exit
    // also carries evanescent modes with complex, quadrature-shifted E/H.
    let exit = material(0.01, 1.0);
    let stacks = [1.7, 2.3].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&exit),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let lines: Vec<_> = stacks
        .iter()
        .enumerate()
        .map(|(i, stack)| {
            let phase = C64::new(0.4 + 0.2 * i as f64, -0.3);
            prepare(
                stack,
                &[
                    [C64::ONE, phase],
                    [phase, C64::ZERO],
                    [C64::new(-0.2, 0.4), C64::from_re(-0.6)],
                    [C64::ZERO, phase],
                ],
            )
        })
        .collect();
    let mean = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        Some(&mean),
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let width = grid().width as f64 * grid().spacing_m[0];
    let height = grid().height as f64 * grid().spacing_m[1];
    let full = AngularSpectrumAperture::try_new([0.0; 2], [width, height]).unwrap();
    let left = AngularSpectrumAperture::try_new([0.0; 2], [11.0, height]).unwrap();
    let right = AngularSpectrumAperture::try_new([11.0, 0.0], [width - 11.0, height]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.51, profile).unwrap();
        for side in [
            Side::Incident,
            Side::Reflected,
            Side::TotalIncidentSide,
            Side::Transmitted,
        ] {
            let original = mean.field(side, 0.19, window, || false).unwrap();
            let measured = mean
                .measure_aperture(side, 0.19, full, window, aperture_limits(), || false)
                .unwrap();
            near(measured.normal_power_w, original.normal_power_w(), 5e-12);
            let original = ensemble.field(side, 0.19, window, || false).unwrap();
            let measured = ensemble
                .measure_aperture(side, 0.19, full, window, aperture_limits(), || false)
                .unwrap();
            near(measured.normal_power_w, original.normal_power_w(), 1e-11);
            let a = ensemble
                .measure_aperture(side, 0.19, left, window, aperture_limits(), || false)
                .unwrap();
            let b = ensemble
                .measure_aperture(side, 0.19, right, window, aperture_limits(), || false)
                .unwrap();
            near(
                a.normal_power_w + b.normal_power_w,
                measured.normal_power_w,
                1e-11,
            );
            for ((a, b), total) in values(a.mean)
                .into_iter()
                .zip(values(b.mean))
                .zip(values(measured.mean))
            {
                near(
                    (a * left.area_m2() + b * right.area_m2()) / full.area_m2(),
                    total,
                    3e-13,
                );
            }
        }
    }
}

#[test]
fn g1_proper_two_frequency_two_angle_source_has_joint_detector_visibility() {
    let medium = material(1.0, 1.0);
    let stacks = [2.0, 3.0].map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let kx = core::f64::consts::TAU / (grid().width as f64 * grid().spacing_m[0]);
    let samples: Vec<_> = (0..grid().width)
        .map(|i| {
            let angle = -kx * grid().spacing_m[0] * i as f64;
            let (s, c) = angle.sin_cos();
            [C64::ZERO, C64::new(c, s)]
        })
        .collect();
    let lines = [
        constant(&stacks[0], [C64::ZERO, C64::ONE]),
        prepare(&stacks[1], &samples),
    ];
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
    let admittance = (9.0 - kx * kx).sqrt() / 3.0;
    let sinc = |x: f64| x.sin() / x;
    for (width, duration) in [(5.7, 1.4), (32.0, 1.4), (5.7, core::f64::consts::TAU)] {
        let aperture = AngularSpectrumAperture::try_new([4.1, -0.2], [width, 0.9]).unwrap();
        let window = HarmonicWindow::new(0.37, duration).unwrap();
        let fringe = sinc(kx * width / 2.0)
            * sinc(duration / 2.0)
            * ((0.37 + duration / 2.0) - kx * (4.1 + width / 2.0)).cos();
        let reading = ensemble
            .measure_aperture(
                Side::Incident,
                0.0,
                aperture,
                window,
                aperture_limits(),
                || false,
            )
            .unwrap();
        near(
            reading.mean.electric_mean_squared_v2_m2[1],
            1.0 + fringe,
            2e-13,
        );
        near(
            reading.mean.magnetic_mean_squared_a2_m2[0],
            0.5 * (1.0 + admittance * admittance) + admittance * fringe,
            2e-13,
        );
        let normal = 0.5 * (1.0 + admittance) * (1.0 + fringe);
        near(reading.mean.poynting_w_m2[2], normal, 2e-13);
        near(
            reading.mean.poynting_w_m2[0],
            (kx / 6.0) * (1.0 + fringe),
            2e-13,
        );
        near(reading.normal_power_w, normal * aperture.area_m2(), 2e-12);
    }
}

#[test]
fn g1_two_face_mean_and_mixed_factors_match_original_uniform_source_moments() {
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
    let statistics = [
        AngularSpectrumFactorStatistics::Real,
        AngularSpectrumFactorStatistics::ProperComplex,
    ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left,
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let area = AngularSpectrumAperture::try_new([-0.2, 0.7], [0.83, 0.41]).unwrap();
    let window =
        HarmonicWindow::with_profile(0.11, 0.37, HarmonicWindowProfile::Triangular).unwrap();
    for side in [
        Face::IncomingLeft,
        Face::OutgoingLeft,
        Face::TotalLeft,
        Face::IncomingRight,
        Face::OutgoingRight,
        Face::TotalRight,
    ] {
        let original = mean
            .illumination_field(side, 0.17, window, || false)
            .unwrap();
        let observed = mean
            .measure_illumination_aperture(side, 0.17, area, window, aperture_limits(), || false)
            .unwrap();
        for (a, b) in values(observed.mean)
            .into_iter()
            .zip(values(original.points()[0]))
        {
            near(a, b, 3e-14);
        }
        let original = ensemble
            .illumination_field(side, 0.17, window, || false)
            .unwrap();
        let observed = ensemble
            .measure_illumination_aperture(side, 0.17, area, window, aperture_limits(), || false)
            .unwrap();
        for (a, b) in values(observed.mean)
            .into_iter()
            .zip(values(original.points()[0]))
        {
            near(a, b, 7e-14);
        }
        near(
            observed.normal_power_w,
            original.points()[0].poynting_w_m2[2] * area.area_m2(),
            3e-14,
        );
    }
    assert!(
        ensemble
            .measure_aperture(
                Side::Incident,
                0.0,
                area,
                window,
                AngularSpectrumExposureApertureLimits {
                    maximum_mode_frequency_pairs: 47,
                    ..aperture_limits()
                },
                || false
            )
            .is_err()
    );
    assert!(matches!(
        mean.measure_illumination_aperture(
            Face::TotalRight,
            0.0,
            area,
            window,
            aperture_limits(),
            || true
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let dark = AngularSpectrumExposure::try_new(&left, &[], limits(), || false).unwrap();
    let darkness = dark
        .measure_aperture(
            Side::Incident,
            0.0,
            area,
            window,
            AngularSpectrumExposureApertureLimits {
                maximum_mode_frequency_pairs: 0,
                ..aperture_limits()
            },
            || false,
        )
        .unwrap();
    assert_eq!(values(darkness.mean), [0.0; 9]);
    assert_eq!(darkness.normal_power_w, 0.0);
}
