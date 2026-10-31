//! A resolved near-camera scattering layer inside a long, optically thick ray.
//! The reference is first-order airlight, not a multiple-scattering solution.

use super::*;

fn reference_cd_m2(panels: u32) -> (f64, f64) {
    // Source at z=-a, camera at z=0 looking along +Z: r=a+t. Both
    // propagation legs attenuate, so
    // L = I*sigma_s/(4*pi) * integral exp[-sigma_t*(a+2t)]/(a+t)^2 dt.
    // For a=.01 m, sigma_t=100/m, substitute u=t/a. Integrate only u=0..40
    // with independent Simpson weights; the complete ray ends at u=10000.
    // No production source, phase, attenuation or quadrature kernel is used.
    assert!(panels > 0 && panels % 2 == 0);
    let step = 40.0 / f64::from(panels);
    let integrand = |u: f64| (-1.0 - 2.0 * u).exp() / (1.0 + u).powi(2);
    let mut sum = integrand(0.0) + integrand(40.0);
    for panel in 1..panels {
        sum += if panel % 2 == 0 { 2.0 } else { 4.0 } * integrand(f64::from(panel) * step);
    }
    let prefactor = 100.0 * 80.0 / (4.0 * PI * 0.01);
    // For u>=40, (1+u)^-2 <=41^-2. Extending the remaining exponential
    // integral to infinity bounds the omitted positive tail analytically.
    // This does not certify Simpson arithmetic or the production integrator.
    let tail_upper = prefactor * (-81.0_f64).exp() / (2.0 * 41.0_f64.powi(2));
    (prefactor * sum * step / 3.0, tail_upper)
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
fn g1_thick_haze_resolves_near_camera_airlight_and_reports_qualified_refinement() {
    let (expected, omitted_tail) = reference_cd_m2(16_384);
    close(reference_cd_m2(8192).0, expected, 2.0e-9);
    assert!(omitted_tail > 0.0 && omitted_tail < 2.0e-34);
    with_cx(|cx, _| {
        let sources = [source(
            [0.0, 0.0, -0.01],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(20.0, 80.0, PhaseFunction::Isotropic).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let settings = ScatterSettings::default();
        let result = scene
            .luminance_along_ray(cx, &forward_ray(), 100.0, settings)
            .unwrap();
        // About 6495.34889475 cd/m² lives in the first few centimetres.
        // Uniform 32/64-panel distance grids formerly sampled only negligible
        // light and incorrectly accepted their tiny mutual discrepancy.
        close(result.luminance_cd_m2, expected, 5.0e-5);
        assert!(result.single_scatter.tolerance_met);
        assert!(result.single_scatter.panels > settings.initial_panels);
        assert!(result.single_scatter.panels <= settings.max_panels);
        assert!(result.single_scatter.evaluations < 4 * u64::from(settings.max_panels));
        check_unknown_bounds(&result);

        // Distance intervals refine independently; their transmitted absolute
        // discrepancies are summed. Half the aggregate count is not a previous
        // grid. The generic single-interval regression checks exact differences.
        let discrepancy = result
            .single_scatter
            .estimated_quadrature_error_cd_m2
            .unwrap();
        assert!(discrepancy > 0.0);
        assert!(
            discrepancy <= settings.absolute_tolerance_cd_m2
                || discrepancy / result.single_scatter.luminance_cd_m2
                    <= settings.relative_tolerance
        );

        for max_panels in [1, 2] {
            assert!(matches!(
                scene.luminance_along_ray(
                    cx,
                    &forward_ray(),
                    100.0,
                    ScatterSettings {
                        initial_panels: 1,
                        max_panels,
                        ..settings
                    },
                ),
                Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: rejected })
                    if rejected == max_panels
            ));
        }
    });
}

#[test]
fn g3_thick_haze_preserves_si_length_extinction_and_intensity_scaling() {
    let expected = reference_cd_m2(16_384).0;
    with_cx(|cx, _| {
        let mut baseline = None;
        for scale in [1.0, 1.0 / 1024.0, 1024.0] {
            let sources = [source(
                [0.0, 0.0, -0.01 * scale],
                100.0 * scale * scale,
                SourceSpectrum::equal_energy(),
            )];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(20.0 / scale, 80.0 / scale, PhaseFunction::Isotropic)
                    .unwrap(),
                TransportOptions {
                    ray_epsilon_m: 1.0e-6 * scale,
                    max_distance_m: 256.0 * scale,
                    ..TransportOptions::default()
                },
                cx,
            )
            .unwrap();
            let value = scene
                .luminance_along_ray(
                    cx,
                    &forward_ray(),
                    100.0 * scale,
                    ScatterSettings::default(),
                )
                .unwrap();
            // I/r², sigma_t*distance and sigma_s*dt all remain unchanged.
            close(value.luminance_cd_m2, expected, 5.0e-5);
            close(
                value.luminance_cd_m2,
                *baseline.get_or_insert(value.luminance_cd_m2),
                2.0e-12,
            );
            assert!(value.single_scatter.tolerance_met);
            check_unknown_bounds(&value);
        }
    });
}
