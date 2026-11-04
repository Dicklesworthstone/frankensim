//! Public integration through original Maxwell fields and material admission.
use super::*;
use fs_em::plane_wave::angular_spectrum::AngularSpectrumExposureStressLimits;

#[test]
fn finite_stress_patches_partition_the_original_continuous_cell_tensor() {
    let medium = material(1.0, 1.0);
    let film = material(2.25, 1.0);
    let stacks = [1.7, 2.3].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
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
        .map(|stack| {
            prepare(
                stack,
                &[
                    [C64::ONE, C64::new(0.4, -0.3)],
                    [C64::new(-0.1, 0.7), C64::ZERO],
                    [C64::new(-0.2, 0.4), C64::from_re(-0.6)],
                    [C64::ZERO, C64::new(0.6, 0.2)],
                ],
            )
        })
        .collect();
    let mean = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        Some(&mean),
        &factors,
        &[AngularSpectrumFactorStatistics::ProperComplex],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let width = grid().width as f64 * grid().spacing_m[0];
    let height = grid().height as f64 * grid().spacing_m[1];
    let whole = AngularSpectrumAperture::try_new([0.0; 2], [width, height]).unwrap();
    let left = AngularSpectrumAperture::try_new([0.0; 2], [0.37 * width, height]).unwrap();
    let right =
        AngularSpectrumAperture::try_new([0.37 * width, 0.0], [0.63 * width, height]).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.51, profile).unwrap();
        for side in [Side::TotalIncidentSide, Side::Transmitted] {
            let native = ensemble
                .stress_field(
                    side,
                    0.19,
                    window,
                    AngularSpectrumExposureStressLimits {
                        maximum_workspace_bytes: 8 << 20,
                    },
                    || false,
                )
                .unwrap();
            let result = ensemble
                .measure_aperture_stress(side, 0.19, whole, window, aperture_limits(), || false)
                .unwrap();
            let a = ensemble
                .measure_aperture_stress(side, 0.19, left, window, aperture_limits(), || false)
                .unwrap();
            let b = ensemble
                .measure_aperture_stress(side, 0.19, right, window, aperture_limits(), || false)
                .unwrap();
            let total = result.integrated_stress_n().unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    near(total[i][j], native.integrated_stress_n()[i][j], 3e-11);
                    near(
                        total[i][j],
                        a.integrated_stress_n().unwrap()[i][j]
                            + b.integrated_stress_n().unwrap()[i][j],
                        3e-11,
                    );
                }
            }
            near(
                result.mean_field_energy_density_j_m3() * whole.area_m2(),
                a.mean_field_energy_density_j_m3() * left.area_m2()
                    + b.mean_field_energy_density_j_m3() * right.area_m2(),
                3e-11,
            );
            let original_flux = ensemble
                .measure_aperture(side, 0.19, whole, window, aperture_limits(), || false)
                .unwrap();
            near(
                result.normal_power_w().unwrap(),
                original_flux.normal_power_w,
                3e-11,
            );
        }
    }
}

#[test]
fn original_two_face_standing_node_has_stress_without_energy_transport() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let left = constant(&stack, [C64::ONE, C64::ZERO]);
    let sources = [left
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
    let bank =
        AngularSpectrumExposure::try_new_illumination(&left, &sources, limits(), || false).unwrap();
    let area = AngularSpectrumAperture::try_new([0.2, -0.7], [0.4, 0.9]).unwrap();
    let window = HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap();
    let result = bank
        .measure_illumination_aperture_stress(
            Face::TotalLeft,
            0.0,
            area,
            window,
            aperture_limits(),
            || false,
        )
        .unwrap();
    near(result.mean_field_energy_density_j_m3(), 1.0, 3e-13);
    near(result.normal_power_w().unwrap(), 0.0, 3e-13);
    near(result.mean_maxwell_stress_pa()[2][2], -1.0, 3e-13);
    assert!(
        bank.measure_illumination_aperture_stress(
            Face::TotalLeft,
            0.0,
            area,
            window,
            aperture_limits(),
            || true
        )
        .is_err()
    );
    assert!(
        bank.measure_illumination_aperture_stress(
            Face::TotalLeft,
            0.0,
            area,
            window,
            AngularSpectrumExposureApertureLimits {
                maximum_mode_frequency_pairs: 0,
                ..aperture_limits()
            },
            || false
        )
        .is_err()
    );
    let replay = bank
        .measure_illumination_aperture_stress(
            Face::TotalLeft,
            0.0,
            area,
            window,
            aperture_limits(),
            || false,
        )
        .unwrap();
    assert_eq!(result, replay);
}
