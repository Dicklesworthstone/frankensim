//! G1 finite-shutter Maxwell stress of an arbitrarily directed plane wave.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposureApertureLimits,
    AngularSpectrumExposureStressLimits, AngularSpectrumFactorStatistics,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

#[test]
fn g1_carrier_stress_retains_real_spatial_sums_and_proper_cancellation() {
    let medium = material(1.0, 1.0);
    let source =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let cell = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [5.0, 3.0],
    };
    let carrier = [0.3_f64, 0.2];
    let kt = (carrier[0] * carrier[0] + carrier[1] * carrier[1]).sqrt();
    let q = (1.0 - kt * kt).sqrt();
    let direction = [carrier[0], carrier[1], q];
    let te = [
        C64::from_re(-carrier[1] / kt),
        C64::from_re(carrier[0] / kt),
    ];
    let lines = [PlanarAngularSpectrum::prepare_with_carrier(
        &source,
        cell,
        carrier,
        &[te; 2],
        AngularSpectrumLimits {
            maximum_grid_points: 2,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let bank = &factors[0];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors,
        &[AngularSpectrumFactorStatistics::ProperComplex],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let distance = 0.27;
    let area = 30.0;
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.37, profile).unwrap();
        let time_weight = match profile {
            HarmonicWindowProfile::Uniform => sinc(window.duration_s()),
            HarmonicWindowProfile::Triangular => sinc(window.duration_s() * 0.5).powi(2),
        };
        let time_phase = window.start_s() + 0.5 * window.duration_s() - q * distance;
        let spatial_weight = sinc(carrier[0] * 10.0) * sinc(carrier[1] * 3.0);
        let center_phase = time_phase - carrier[0] * 5.0 - carrier[1] * 1.5;
        // For unit TE E and epsilon=mu=1, the real Maxwell stress is
        // -k_hat outer k_hat * cos²(t-k.r). The rectangle and shutter
        // average of its double-frequency cosine separates into sinc factors.
        let mean_square = 0.5 * (1.0 + (2.0 * center_phase).cos() * time_weight * spatial_weight);
        let exact = AngularSpectrumExposureStressLimits {
            maximum_workspace_bytes: bank.stress_workspace_bytes().unwrap(),
        };
        let stress = bank
            .stress_field(Side::Transmitted, distance, window, exact, || false)
            .unwrap();
        let proper = ensemble
            .stress_field(
                Side::Transmitted,
                distance,
                window,
                AngularSpectrumExposureStressLimits {
                    maximum_workspace_bytes: ensemble.stress_workspace_bytes().unwrap(),
                },
                || false,
            )
            .unwrap();
        let patch = bank
            .measure_aperture_stress(
                Side::Transmitted,
                distance,
                AngularSpectrumAperture::try_new([0.0; 2], [10.0, 3.0]).unwrap(),
                window,
                AngularSpectrumExposureApertureLimits {
                    maximum_mode_frequency_pairs: 4,
                    maximum_workspace_bytes: 2 << 20,
                },
                || false,
            )
            .unwrap();
        let integrated = stress.integrated_stress_n();
        let integrated_proper = proper.integrated_stress_n();
        let integrated_patch = patch.integrated_stress_n().unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let tensor = -direction[i] * direction[j];
                near(integrated[i][j], area * tensor * mean_square, 3e-12);
                near(integrated_proper[i][j], area * tensor * 0.5, 3e-12);
                near(integrated[i][j], integrated_patch[i][j], 3e-12);
                for (point, actual) in stress.points().iter().enumerate() {
                    let phase = time_phase - carrier[0] * point as f64 * cell.spacing_m[0];
                    let local = 0.5 * (1.0 + (2.0 * phase).cos() * time_weight);
                    near(actual.maxwell_stress_pa[i][j], tensor * local, 2e-13);
                }
            }
        }
        assert!(
            bank.stress_field(
                Side::Transmitted,
                distance,
                window,
                AngularSpectrumExposureStressLimits {
                    maximum_workspace_bytes: exact.maximum_workspace_bytes - 1
                },
                || false
            )
            .is_err()
        );
        assert!(matches!(
            bank.stress_field(Side::Transmitted, distance, window, exact, || true),
            Err(AngularSpectrumError::Cancelled)
        ));
        assert!(
            bank.boundary_stress_resultant_n([0.0; 2], window, exact, || false)
                .is_err()
        );
    }
}
