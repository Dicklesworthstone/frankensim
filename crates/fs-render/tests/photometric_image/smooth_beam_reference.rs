//! Independent single-scattering references for a smooth, sub-degree beam.
//! Intensity contours guide quadrature; discrepancies remain heuristic and
//! neither fixture calibration nor omitted multiple scattering is certified.

use super::*;

const EXPONENT: f64 = 20_000.0;
const PEAK_CD: f64 = 10_000.0;
const LENGTH_M: f64 = 10.0;

// Source at (-1,0,5), optical axis tilted by beta from +X toward +Z.
// theta=atan(t-5) gives dt/r²=dtheta and the relative intensity is
// cos(theta-beta)^n. This oracle uses only standard elementary functions,
// declared SI coefficients, and Simpson integration in an independent variable.
fn reference_cd_m2(beta: f64, cutoff: f64, panels: u32) -> f64 {
    let lo = (-5.0_f64).atan().max(beta - cutoff);
    let hi = 5.0_f64.atan().min(beta + cutoff);
    let step = (hi - lo) / f64::from(panels);
    let value = |theta: f64| {
        (EXPONENT * (theta - beta).cos().ln() - 0.15 * (5.0 + theta.tan() + 1.0 / theta.cos()))
            .exp()
    };
    let mut sum = value(lo) + value(hi);
    for index in 1..panels {
        sum += if index % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(index) * step);
    }
    PEAK_CD * 0.1 / (4.0 * PI) * sum * step / 3.0
}

fn fixture(scale: f64, beta: f64, cutoff: f64, rotated: bool) -> (PhotometricSource, Ray) {
    let angle = PI / 2.0 - beta + if rotated { PI / 2.0 } else { 0.0 };
    let (sine, cosine) = (0.5 * angle).sin_cos();
    let (position, ray) = if rotated {
        (
            [2.0 + 5.0 * scale, -1.0, 3.0 + scale],
            Ray {
                origin: Point3::new(2.0, -1.0, 3.0),
                dir: Vec3::new(1.0, 0.0, 0.0),
            },
        )
    } else {
        ([-scale, 0.0, 5.0 * scale], forward_ray())
    };
    // Independent integral of I_axis * cos(theta)^n over the retained cone.
    // Do not derive this flux from the production peak or normalization helper.
    let flux = 2.0 * PI * PEAK_CD * scale * scale * (1.0 - cutoff.cos().powf(EXPONENT + 1.0))
        / (EXPONENT + 1.0);
    let source = PhotometricSource::try_new(
        RigidTransform::try_new([0.0, sine, 0.0, cosine], position).unwrap(),
        flux,
        1.0,
        AngularDistribution::cosine_power(EXPONENT, cutoff).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap();
    let axis = Vec3::new(angle.sin(), 0.0, angle.cos());
    close(
        source.intensity_cd(axis).unwrap(),
        PEAK_CD * scale * scale,
        3e-12,
    );
    (source, ray)
}

fn check_result(result: &LuminanceEstimate, expected: f64) {
    close(result.luminance_cd_m2, expected, 5e-5);
    assert_eq!(result.surface_cd_m2, 0.0);
    assert_eq!(result.numerical.kind, fs_evidence::NumericalKind::NoClaim);
    assert!(result.total_error_bound_cd_m2.is_none());
    assert!(
        result
            .single_scatter
            .omitted_multiple_scattering_bound_cd_m2
            .is_none()
    );
    assert!(result.single_scatter.tolerance_met);
    assert!(result.single_scatter.panels <= ScatterSettings::default().max_panels);
    assert!(
        result.single_scatter.evaluations < 4 * u64::from(ScatterSettings::default().max_panels)
    );
}

#[test]
fn g1_smooth_power_beam_recovers_light_missed_by_distance_only_refinement() {
    let expected = reference_cd_m2(0.0, PI / 2.0, 65_536);
    close(reference_cd_m2(0.0, PI / 2.0, 32_768), expected, 2e-11);
    close(expected, 0.573_446_921_506_254_4, 2e-11);
    with_cx(|cx, _| {
        let (source, ray) = fixture(1.0, 0.0, PI / 2.0, false);
        let sources = [source];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        // The old distance intervals [0,1,3,4,5,6,7,9,10] put four/eight
        // midpoint samples in each interval. They returned ~7.38e-67/~9.43e-17
        // cd/m² and accepted the default discrepancy. The pure contour-helper
        // regression executes that old sample pattern against this oracle.
        let result = scene
            .luminance_along_ray(cx, &ray, LENGTH_M, ScatterSettings::default())
            .unwrap();
        check_result(&result, expected);
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .unwrap()
                > 0.0
        );
        let replay = scene
            .luminance_along_ray(cx, &ray, LENGTH_M, ScatterSettings::default())
            .unwrap();
        assert_eq!(
            result.luminance_cd_m2.to_bits(),
            replay.luminance_cd_m2.to_bits()
        );
        assert_eq!(result.single_scatter.panels, replay.single_scatter.panels);
    });
}

#[test]
fn g3_smooth_beam_retains_tilt_cutoff_flux_and_si_rigid_covariance() {
    let (beta, cutoff) = (0.35, 0.03);
    let expected = reference_cd_m2(beta, cutoff, 32_768);
    close(reference_cd_m2(beta, cutoff, 16_384), expected, 2e-11);
    with_cx(|cx, _| {
        let mut baseline = None;
        for (scale, rotated) in [(1.0, false), (0.125, true), (8.0, true)] {
            let (source, ray) = fixture(scale, beta, cutoff, rotated);
            let sources = [source];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(0.05 / scale, 0.1 / scale, PhaseFunction::Isotropic)
                    .unwrap(),
                TransportOptions {
                    ray_epsilon_m: 1e-6 * scale,
                    max_distance_m: 20.0 * scale,
                    ..TransportOptions::default()
                },
                cx,
            )
            .unwrap();
            let result = scene
                .luminance_along_ray(cx, &ray, LENGTH_M * scale, ScatterSettings::default())
                .unwrap();
            check_result(&result, expected);
            close(
                result.luminance_cd_m2,
                *baseline.get_or_insert(result.luminance_cd_m2),
                5e-10,
            );
        }
    });
}

#[test]
fn g4_smooth_contours_share_the_cap_and_inactive_sources_add_no_work() {
    with_cx(|cx, gate| {
        let (source, ray) = fixture(1.0, 0.0, PI / 2.0, false);
        let inactive = [source.filtered([0.0; SPECTRAL_BINS]).unwrap()];
        let sources = [source.clone(), source];
        let medium = HomogeneousMedium::try_new(0.05, 0.1, PhaseFunction::Isotropic).unwrap();
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            medium.clone(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let settings = ScatterSettings {
            initial_panels: 1,
            max_panels: 38,
            ..ScatterSettings::default()
        };
        // Seven distance cuts plus thirty power-contour crossings. Identical
        // fixtures share those knots. One panel per interval supplies no
        // refinement comparison and cannot qualify the reported discrepancy.
        let result = scene
            .luminance_along_ray(cx, &ray, LENGTH_M, settings)
            .unwrap();
        assert_eq!(
            (
                result.single_scatter.panels,
                result.single_scatter.evaluations
            ),
            (38, 38)
        );
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .is_none()
        );
        assert!(!result.single_scatter.tolerance_met);
        assert!(matches!(
            scene.luminance_along_ray(
                cx,
                &ray,
                LENGTH_M,
                ScatterSettings {
                    max_panels: 37,
                    ..settings
                },
            ),
            Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: 37 })
        ));
        let complete = scene
            .luminance_along_ray(cx, &ray, LENGTH_M, ScatterSettings::default())
            .unwrap();
        check_result(&complete, 2.0 * reference_cd_m2(0.0, PI / 2.0, 32_768));
        let inactive_scene =
            PhotometricScene::try_new(&[], &inactive, medium, TransportOptions::default(), cx)
                .unwrap();
        let dark = inactive_scene
            .luminance_along_ray(
                cx,
                &ray,
                LENGTH_M,
                ScatterSettings {
                    max_panels: 1,
                    ..settings
                },
            )
            .unwrap();
        assert_eq!(dark.luminance_cd_m2, 0.0);
        assert_eq!(dark.single_scatter.panels, 1);
        gate.request();
        assert!(matches!(
            scene.luminance_along_ray(cx, &ray, LENGTH_M, settings),
            Err(PhotometryError::Cancelled)
        ));
    });
}
