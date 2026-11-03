//! Native analytic emitting geometry through a real refracting interface.

use super::*;

#[test]
fn g1_native_sphere_source_refracts_with_area_density_and_complementary_camera_weight() {
    with_cx(|cx| {
        let sphere = GeometryInstance::try_new(
            9, hash_domain("native-sphere-area-camera", b"original-chart"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, -2.0), radius: 0.5,
            }), RigidTransform::identity(),
        ).unwrap();
        let geometry = InstanceScene::try_new(vec![
            panel(1, 0.0, 1.0, 10.0), panel(2, 0.0, -0.5, 10.0), sphere,
        ]).unwrap();
        let air = glass(1.0, 0.1);
        let solid = glass(1.5, 0.2);
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&solid)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &reflectances,
            vec![
                diffuse(1, 0),
                SurfaceBinding { object_id: 2, interaction: Interaction::Dielectric {
                    outside: 0, inside: 1, surface: DielectricSurface::SMOOTH,
                }},
                opaque(9),
            ], 100_000,
        ).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let selected = source.sample_area_with_context(cx, [0.0, 0.0]).unwrap();
        assert_eq!(selected.triangle_index(), None);
        assert_eq!(selected.position_m(), Point3::new(0.0, 0.0, -1.5));
        assert_eq!(selected.normal(), Vec3::new(0.0, 0.0, 1.0));
        near(source.total_area_m2(), core::f64::consts::PI, 2e-14);
        near(selected.pdf_area_per_m2(), 1.0 / core::f64::consts::PI, 2e-14);
        let continuation = |event| match event {
            0 => [0.0; 3],
            1 => [0.0, 0.0, 0.999],
            _ => panic!("only the diffuse receiver and smooth interface scatter"),
        };
        let baseline = camera.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
            settings().surfaces, continuation).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        for represented in [true, false] {
            let direction = if represented { 1.0 } else { -1.0 };
            let seeds = [seed(Vec3::new(0.0, 0.0, direction))];
            let bindings = [AreaLightBinding {
                source: &source, receiver_object_id: 1, source_side: BoundarySide::Positive,
                initial_medium: 1, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let mut native_draws = 0;
            let mut area_draws = 0;
            let result = owner.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace,
                |event| { native_draws += 1; continuation(event) },
                |event, binding| {
                    assert_eq!((event, binding), (0, 0));
                    area_draws += 1;
                    [0.0, 0.0]
                }).unwrap();
            // Independent paraxial Snell map about this exact axial ray:
            // one metre in n=1.5 followed by 1.5 metres in n=1. Each source
            // angle therefore spans 1 + 1.5*1.5 metres at the receiver.
            let jacobian = (1.0_f64 + 1.5 * 1.5).powi(2);
            let transmission = 1.0 - ((1.5_f64 - 1.0) / (1.5 + 1.0)).powi(2);
            let eta2 = (1.0_f64 / 1.5).powi(2);
            let light_pdf = eta2 * jacobian / core::f64::consts::PI;
            let continuation_pdf = transmission / core::f64::consts::PI;
            let denominator = light_pdf * light_pdf + continuation_pdf * continuation_pdf;
            let light_weight = light_pdf * light_pdf / denominator;
            let hit_weight = continuation_pdf * continuation_pdf / denominator;
            // Source intensity is Le*A, Fresnel power is T, and Beer depth
            // includes the original camera prefix once. The sampled camera
            // interface cancels T/qT, then local source radiance uses eta².
            let raw_area = RHO * LE * transmission / jacobian * (-0.45_f64).exp();
            let raw_hit = RHO * LE * eta2 * (-0.45_f64).exp();
            let expected_area = if represented { raw_area * light_weight } else { 0.0 };
            let expected_hit = raw_hit * if represented { hit_weight } else { 1.0 };
            near(baseline.radiance_w_m2_sr_nm(), raw_hit, 4e-6);
            near(result.area_radiance_w_m2_sr_nm(), expected_area, 8e-6);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + expected_hit, 8e-6);
            assert_eq!(result.radiance_sample().path(), baseline.path());
            assert_eq!((native_draws, area_draws, result.area_samples()), (2, 1, 1));
            assert_eq!(result.images(), u32::from(represented));
            let repeated = owner.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace, continuation, |_, _| [0.0, 0.0]).unwrap();
            assert_eq!(result.radiance_sample(), repeated.radiance_sample());
            assert_eq!(result.area_radiance_w_m2_sr_nm(), repeated.area_radiance_w_m2_sr_nm());
        }
    });
}
