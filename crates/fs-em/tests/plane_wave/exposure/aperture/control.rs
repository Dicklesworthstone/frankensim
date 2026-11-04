//! Original Maxwell/material preparation and finite-footprint source controls.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumApertureControlLimits, PreparedAngularSpectrumApertureEnergyControl,
};

fn control_limits() -> AngularSpectrumApertureControlLimits {
    AngularSpectrumApertureControlLimits {
        maximum_controls: 4,
        maximum_point_line_pairs: 4096,
        maximum_source_field_entries: 4096,
        maximum_kernel_pairs: 4096,
        maximum_workspace_bytes: 16 << 20,
    }
}

#[test]
fn original_film_fields_and_source_perturbations_match_the_prepared_energy_objective() {
    let medium = material(1.0, 1.0);
    let film = material(2.25, 1.0);
    let exit = material(0.01, 1.0);
    let stacks = [2.0, 3.0].map(|omega| {
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
    let fixed_samples = [
        [C64::ONE, C64::new(0.4, -0.3)],
        [C64::new(-0.1, 0.7), C64::ZERO],
        [C64::new(-0.2, 0.4), C64::from_re(-0.6)],
        [C64::ZERO, C64::new(0.6, 0.2)],
    ];
    let source_samples = [
        [C64::new(0.2, 0.3), C64::new(-0.2, 0.1)],
        [C64::new(0.3, -0.1), C64::ONE],
        [C64::new(-0.3, -0.4), C64::ZERO],
        [C64::new(0.7, 0.4), C64::new(-0.2, 0.4)],
    ];
    let fixed_lines: Vec<_> = stacks.iter().map(|s| prepare(s, &fixed_samples)).collect();
    let source_lines: Vec<_> = stacks.iter().map(|s| prepare(s, &source_samples)).collect();
    let fixed = AngularSpectrumExposure::try_new(&fixed_lines[0], &fixed_lines, limits(), || false)
        .unwrap();
    let controls =
        [
            AngularSpectrumExposure::try_new(&fixed_lines[0], &source_lines, limits(), || false)
                .unwrap(),
        ];
    let aperture = AngularSpectrumAperture::try_new([0.3, -0.7], [0.83, 0.41]).unwrap();
    let window =
        HarmonicWindow::with_profile(0.19, 0.47, HarmonicWindowProfile::Triangular).unwrap();
    let coefficient = C64::new(0.3, -0.7);
    for side in [Side::TotalIncidentSide, Side::Transmitted] {
        let actual = |a: C64| {
            let samples: Vec<_> = fixed_samples
                .iter()
                .zip(&source_samples)
                .map(|(fixed, source)| [fixed[0] + a * source[0], fixed[1] + a * source[1]])
                .collect();
            let lines: Vec<_> = stacks.iter().map(|s| prepare(s, &samples)).collect();
            let bank =
                AngularSpectrumExposure::try_new(&fixed_lines[0], &lines, limits(), || false)
                    .unwrap();
            bank.measure_aperture_stress(side, 0.17, aperture, window, aperture_limits(), || false)
                .unwrap()
                .mean_field_energy_density_j_m3()
        };
        let mut prepared = PreparedAngularSpectrumApertureEnergyControl::new(
            &fixed,
            &controls,
            side,
            0.17,
            aperture,
            window,
            control_limits(),
            || false,
        )
        .unwrap();
        let result = prepared.value_gradient(&[coefficient], || false).unwrap();
        near(
            result.mean_field_energy_density_j_m3,
            actual(coefficient),
            5e-12,
        );
        let gradient = result.source_gradient[0];
        let h = 1e-4;
        near(
            gradient.re,
            (actual(coefficient + C64::from_re(h)) - actual(coefficient - C64::from_re(h)))
                / (2.0 * h),
            2e-9,
        );
        near(
            gradient.im,
            (actual(coefficient + C64::new(0.0, h)) - actual(coefficient - C64::new(0.0, h)))
                / (2.0 * h),
            2e-9,
        );
        assert_eq!(prepared.frequencies_rad_s(), &[2.0, 3.0]);
        assert_eq!(prepared.source_field_entries(), 16);
        assert_eq!(prepared.kernel_pairs(), 56);
    }
    for limits in [
        AngularSpectrumApertureControlLimits {
            maximum_controls: 0,
            ..control_limits()
        },
        AngularSpectrumApertureControlLimits {
            maximum_point_line_pairs: 0,
            ..control_limits()
        },
        AngularSpectrumApertureControlLimits {
            maximum_source_field_entries: 0,
            ..control_limits()
        },
        AngularSpectrumApertureControlLimits {
            maximum_kernel_pairs: 0,
            ..control_limits()
        },
        AngularSpectrumApertureControlLimits {
            maximum_workspace_bytes: 0,
            ..control_limits()
        },
    ] {
        assert!(
            PreparedAngularSpectrumApertureEnergyControl::new(
                &fixed,
                &controls,
                Side::Incident,
                0.0,
                aperture,
                window,
                limits,
                || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedAngularSpectrumApertureEnergyControl::new(
            &fixed,
            &controls,
            Side::Incident,
            0.0,
            aperture,
            window,
            control_limits(),
            || true,
        )
        .is_err()
    );
}

#[test]
fn two_face_standing_fields_control_together_and_material_identity_is_retained() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let left = constant(&stack, [C64::ONE, C64::ZERO]);
    let illumination = [left
        .prepare_illumination(
            &[[C64::from_re(-1.0), C64::ZERO]; 4],
            AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 4096,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap()];
    let fixed =
        AngularSpectrumExposure::try_new_illumination(&left, &illumination, limits(), || false)
            .unwrap();
    let controls =
        [
            AngularSpectrumExposure::try_new_illumination(&left, &illumination, limits(), || false)
                .unwrap(),
        ];
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let mut prepared = PreparedAngularSpectrumApertureEnergyControl::new_illumination(
        &fixed,
        &controls,
        Face::TotalLeft,
        0.0,
        aperture,
        window,
        control_limits(),
        || false,
    )
    .unwrap();
    let result = prepared.value_gradient(&[C64::ZERO], || false).unwrap();
    near(result.mean_field_energy_density_j_m3, 1.0, 3e-13);
    near(result.source_gradient[0].re, 2.0, 3e-13);
    let result = prepared
        .value_gradient(&[C64::from_re(-1.0)], || false)
        .unwrap();
    assert_eq!(result.mean_field_energy_density_j_m3, 0.0);
    assert_eq!(result.source_gradient, &[C64::ZERO]);
    let foreign_material = material(1.0, 1.0);
    let foreign_stack = PlanarDielectricLayers::try_new(
        &foreign_material,
        DielectricResponse::Lossless(&foreign_material),
        1.0,
        &[],
    )
    .unwrap();
    let foreign_lines = [constant(&foreign_stack, [C64::ONE, C64::ZERO])];
    let foreign =
        [
            AngularSpectrumExposure::try_new(&foreign_lines[0], &foreign_lines, limits(), || false)
                .unwrap(),
        ];
    assert!(
        PreparedAngularSpectrumApertureEnergyControl::new_illumination(
            &fixed,
            &foreign,
            Face::TotalLeft,
            0.0,
            aperture,
            window,
            control_limits(),
            || false,
        )
        .is_err()
    );
}
