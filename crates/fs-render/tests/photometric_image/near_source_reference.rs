//! Independent airlight reference for a nearby isotropic point fixture.
//! Distance partitions resolve inverse-square variation without claiming a
//! rigorous quadrature enclosure or measured-fixture calibration.

use super::*;

const IMPACT_M: f64 = 0.01;
const SOURCE_Z_M: f64 = 5.110918;
const LENGTH_M: f64 = 10.0;

fn reference_cd_m2(panels: u32) -> f64 {
    // theta=atan((t-c)/b), r=b/cos(theta), dt/r²=dtheta/b. Both propagation
    // legs attenuate. This oracle uses no production phase/source/observer,
    // distance partition, Beer kernel, or integration coordinate mapping.
    let lo = (-SOURCE_Z_M / IMPACT_M).atan();
    let hi = ((LENGTH_M - SOURCE_Z_M) / IMPACT_M).atan();
    let step = (hi - lo) / f64::from(panels);
    let value = |theta: f64| {
        let t = SOURCE_Z_M + IMPACT_M * theta.tan();
        let r = IMPACT_M / theta.cos();
        (-0.01 * (t + r)).exp()
    };
    let mut sum = value(lo) + value(hi);
    for index in 1..panels {
        sum += if index % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(index) * step);
    }
    100.0 * 0.01 / (4.0 * PI * IMPACT_M) * sum * step / 3.0
}

fn check_unknown_bounds(value: &LuminanceEstimate) {
    assert_eq!(value.surface_cd_m2, 0.0);
    assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
    assert!(value.total_error_bound_cd_m2.is_none());
    assert!(
        value
            .single_scatter
            .omitted_multiple_scattering_bound_cd_m2
            .is_none()
    );
}

#[test]
fn g1_nearby_isotropic_fixture_resolves_an_aliased_inverse_square_peak() {
    let expected = reference_cd_m2(32_768);
    close(reference_cd_m2(16_384), expected, 3e-10);
    with_cx(|cx, _| {
        let medium = HomogeneousMedium::try_new(0.0, 0.01, PhaseFunction::Isotropic).unwrap();
        let settings = ScatterSettings::default();
        // The bare callback has no fixture geometry. Its 32/64-panel estimates
        // still agree within the requested discrepancy while missing half the
        // light; this regression preserves that limitation explicitly.
        let aliased = fs_render::photometry::medium::single_scattered_luminance(
            &medium,
            &forward_ray(),
            LENGTH_M,
            settings,
            cx,
            |point, _| {
                let r = IMPACT_M.hypot(point.z - SOURCE_Z_M);
                Ok::<_, ()>(100.0 * (-0.01 * r).exp() / (4.0 * PI * r * r))
            },
        )
        .unwrap();
        assert!(aliased.tolerance_met);
        assert_eq!(aliased.panels, 64);
        assert!((aliased.luminance_cd_m2 / expected - 1.0).abs() > 0.45);
        let sources = [source(
            [IMPACT_M, 0.0, SOURCE_Z_M],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let scene =
            PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                .unwrap();
        let result = scene
            .luminance_along_ray(cx, &forward_ray(), LENGTH_M, settings)
            .unwrap();
        close(result.luminance_cd_m2, expected, 5e-5);
        assert!(result.single_scatter.tolerance_met);
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .unwrap()
                > 0.0
        );
        assert!(result.single_scatter.panels <= settings.max_panels);
        assert!(result.single_scatter.evaluations < 4 * u64::from(settings.max_panels));
        check_unknown_bounds(&result);
    });
}

#[test]
fn g3_near_source_airlight_preserves_si_scale_and_rigid_geometry() {
    let expected = reference_cd_m2(32_768);
    with_cx(|cx, _| {
        let mut baseline = None;
        for scale in [1.0, 0.125, 8.0] {
            // Rotate +Z to +X, translate the entire scene, and scale all
            // lengths. I/r², sigma_t*r and sigma_s*dt remain invariant.
            let ray = Ray {
                origin: Point3::new(2.0, -1.0, 3.0),
                dir: Vec3::new(1.0, 0.0, 0.0),
            };
            let sources = [source(
                [2.0 + SOURCE_Z_M * scale, -1.0, 3.0 - IMPACT_M * scale],
                100.0 * scale * scale,
                SourceSpectrum::equal_energy(),
            )];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(0.0, 0.01 / scale, PhaseFunction::Isotropic).unwrap(),
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
            close(result.luminance_cd_m2, expected, 5e-5);
            close(
                result.luminance_cd_m2,
                *baseline.get_or_insert(result.luminance_cd_m2),
                3e-12,
            );
            assert!(result.single_scatter.tolerance_met);
            check_unknown_bounds(&result);
        }
    });
}

#[test]
fn g4_near_source_partitions_share_budgets_and_preserve_singularity_refusal() {
    with_cx(|cx, gate| {
        let sources = [
            source(
                [IMPACT_M, 0.0, SOURCE_Z_M],
                100.0,
                SourceSpectrum::equal_energy(),
            ),
            source(
                [-IMPACT_M, 0.0, LENGTH_M - SOURCE_Z_M],
                100.0,
                SourceSpectrum::equal_energy(),
            ),
        ];
        let medium = HomogeneousMedium::try_new(0.0, 0.01, PhaseFunction::Isotropic).unwrap();
        let scene = PhotometricScene::try_new(
            &[],
            &sources[..1],
            medium.clone(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let settings = ScatterSettings {
            initial_panels: 1,
            max_panels: 20,
            ..ScatterSettings::default()
        };
        // Closest approach plus nine distance scales on either side gives
        // exactly nineteen knots. Twenty panels permit no refinement comparison.
        let result = scene
            .luminance_along_ray(cx, &forward_ray(), LENGTH_M, settings)
            .unwrap();
        assert_eq!(
            (
                result.single_scatter.panels,
                result.single_scatter.evaluations
            ),
            (20, 20)
        );
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .is_none()
        );
        assert!(!result.single_scatter.tolerance_met);
        check_unknown_bounds(&result);
        assert!(matches!(
            scene.luminance_along_ray(
                cx,
                &forward_ray(),
                LENGTH_M,
                ScatterSettings {
                    max_panels: 19,
                    ..settings
                },
            ),
            Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: 19 })
        ));
        let combined = PhotometricScene::try_new(
            &[],
            &sources,
            medium.clone(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            combined.luminance_along_ray(cx, &forward_ray(), LENGTH_M, settings),
            Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: 20 })
        ));
        let singular = [source(
            [0.0, 0.0, SOURCE_Z_M],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let singular_scene =
            PhotometricScene::try_new(&[], &singular, medium, TransportOptions::default(), cx)
                .unwrap();
        assert!(matches!(
            singular_scene.luminance_along_ray(
                cx,
                &forward_ray(),
                LENGTH_M,
                ScatterSettings::default(),
            ),
            Err(PhotometryError::SourceTooClose { source: 0 })
        ));
        gate.request();
        assert!(matches!(
            scene.luminance_along_ray(cx, &forward_ray(), LENGTH_M, settings),
            Err(PhotometryError::Cancelled)
        ));
    });
}
