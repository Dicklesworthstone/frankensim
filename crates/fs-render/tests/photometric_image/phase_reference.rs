//! Independent single-scattering references for narrow tabulated phase lobes.
//! A phase density is normalized over solid angle; its narrow support is not
//! a fixture beam, and must be resolved even for an isotropic point fixture.

use super::*;

const HALF_WIDTH: f64 = 0.001;
const CROSSING_M: f64 = 1.0 / 16.0;
const CANDELA: f64 = 100.0;

fn settings() -> ScatterSettings {
    ScatterSettings {
        initial_panels: 2,
        max_panels: 4096,
        absolute_tolerance_cd_m2: 0.0,
        relative_tolerance: 1.0e-8,
    }
}

fn phase(center: f64) -> PhaseFunction {
    PhaseFunction::tabulated(&[
        (-1.0, 0.0),
        (center - HALF_WIDTH, 0.0),
        (center, 1.0),
        (center + HALF_WIDTH, 0.0),
        (1.0, 0.0),
    ])
    .unwrap()
}

fn source_z(center: f64) -> f64 {
    CROSSING_M + center / (1.0 - center * center).sqrt()
}

fn declared_phase_per_sr(cosine: f64, center: f64) -> f64 {
    // The triangle has area HALF_WIDTH in dmu. Its solid-angle measure
    // is 2*pi*dmu, independently of the production table normalization.
    (1.0 - (cosine - center).abs() / HALF_WIDTH).max(0.0) / (2.0 * PI * HALF_WIDTH)
}

fn reference_cd_m2(center: f64) -> f64 {
    // Camera (0,0,s), source (1,0,c): r²=1+(s-c)² and mu=(c-s)/r.
    // theta=atan(s-c) gives ds/r²=dtheta and mu=-sin(theta).
    // Integrate in angle space with Simpson weights, independently of the
    // production distance-space midpoint/Beer panels, phase evaluator and
    // boundary solver. Split only at the declared phase endpoints and peak.
    const PANELS: u32 = 4096;
    let c = source_z(center);
    let knots = [
        -(center + HALF_WIDTH).asin(),
        -center.asin(),
        -(center - HALF_WIDTH).asin(),
    ];
    let value = |theta: f64| {
        let s = c + theta.tan();
        let r = 1.0 / theta.cos();
        declared_phase_per_sr(-theta.sin(), center) * (-0.15 * (s + r)).exp()
    };
    let integral: f64 = knots
        .windows(2)
        .map(|pair| {
            let step = (pair[1] - pair[0]) / f64::from(PANELS);
            let mut sum = value(pair[0]) + value(pair[1]);
            for i in 1..PANELS {
                sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(pair[0] + f64::from(i) * step);
            }
            sum * step / 3.0
        })
        .sum();
    CANDELA * 0.1 * integral
}

fn check_unknown_bounds(result: &LuminanceEstimate) {
    assert_eq!(result.surface_cd_m2, 0.0);
    assert_eq!(result.numerical.kind, fs_evidence::NumericalKind::NoClaim);
    assert!(result.total_error_bound_cd_m2.is_none());
    assert!(
        result
            .single_scatter
            .omitted_multiple_scattering_bound_cd_m2
            .is_none()
    );
}

#[test]
fn g1_narrow_forward_and_backward_phase_lobes_match_independent_airlight() {
    with_cx(|cx, _| {
        for center in [0.8, -0.8] {
            let c = source_z(center);
            // These are precisely the original unpartitioned 2/4-panel
            // nodes. Both refinement levels miss the complete phase lobe.
            for panels in [2, 4] {
                for i in 0..panels {
                    let s = (f64::from(i) + 0.5) / f64::from(panels);
                    let r = (1.0 + (s - c).powi(2)).sqrt();
                    assert_eq!(declared_phase_per_sr((c - s) / r, center), 0.0);
                }
            }
            let sources = [source(
                [1.0, 0.0, c],
                CANDELA,
                SourceSpectrum::equal_energy(),
            )];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(0.05, 0.1, phase(center)).unwrap(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let result = scene
                .luminance_along_ray(cx, &forward_ray(), 1.0, settings())
                .unwrap();
            let expected = reference_cd_m2(center);
            assert!(expected > 2.0);
            close(result.luminance_cd_m2, expected, 2.0e-7);
            assert!(result.single_scatter.tolerance_met);
            assert!(result.single_scatter.panels <= settings().max_panels);
            assert!(result.single_scatter.evaluations < 4 * u64::from(settings().max_panels));
            check_unknown_bounds(&result);
        }
    });
}

#[test]
fn g3_phase_lobe_partition_preserves_rigid_placement_and_si_scaling() {
    with_cx(|cx, _| {
        let center = 0.8;
        let c = source_z(center);
        let expected = reference_cd_m2(center);
        for scale in [1.0 / 1024.0, 1.0, 1024.0] {
            // Explicit +90-degree Y rotation and translation. The reference
            // uses neither the production rigid transform nor partition roots.
            let ray = Ray {
                origin: Point3::new(2.0 * scale, -3.0 * scale, 4.0 * scale),
                dir: Vec3::new(1.0, 0.0, 0.0),
            };
            let sources = [source(
                [(2.0 + c) * scale, -3.0 * scale, 3.0 * scale],
                CANDELA * scale * scale,
                SourceSpectrum::equal_energy(),
            )];
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::try_new(0.05 / scale, 0.1 / scale, phase(center)).unwrap(),
                TransportOptions {
                    ray_epsilon_m: 1.0e-6 * scale,
                    max_distance_m: 32.0 * scale,
                    ..TransportOptions::default()
                },
                cx,
            )
            .unwrap();
            let result = scene
                .luminance_along_ray(cx, &ray, scale, settings())
                .unwrap();
            // I scales as length²; sigma scales as 1/length. Optical depth,
            // incident illuminance, and the complete airlight remain fixed.
            close(result.luminance_cd_m2, expected, 2.0e-7);
            assert!(result.single_scatter.tolerance_met);
            check_unknown_bounds(&result);
        }
    });
}

#[test]
fn g4_phase_knots_share_the_global_budget_and_keep_missing_refinement_unknown() {
    with_cx(|cx, _| {
        let sources = vec![
            source(
                [1.0, 0.0, source_z(0.8)],
                CANDELA,
                SourceSpectrum::equal_energy(),
            );
            4
        ];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::try_new(0.05, 0.1, phase(0.8)).unwrap(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            scene.luminance_along_ray(
                cx,
                &forward_ray(),
                1.0,
                ScatterSettings {
                    max_panels: 3,
                    ..settings()
                },
            ),
            Err(PhotometryError::ScatteringPartitionBudgetExceeded { max_panels: 3 })
        ));
        // Three distinct phase knots imply four intervals, shared by all
        // four coincident fixtures. One panel per interval cannot supply a
        // refinement discrepancy, even though the positive lobe is resolved.
        // The closest source distance exceeds the 1 m ray length; distance
        // partitioning therefore adds no mandatory intervals to this test.
        let result = scene
            .luminance_along_ray(
                cx,
                &forward_ray(),
                1.0,
                ScatterSettings {
                    initial_panels: 1,
                    max_panels: 4,
                    ..settings()
                },
            )
            .unwrap();
        assert!(result.luminance_cd_m2 > 0.0);
        assert_eq!(result.single_scatter.panels, 4);
        assert_eq!(result.single_scatter.evaluations, 4);
        assert!(
            result
                .single_scatter
                .estimated_quadrature_error_cd_m2
                .is_none()
        );
        assert!(!result.single_scatter.tolerance_met);
        check_unknown_bounds(&result);
    });
}
