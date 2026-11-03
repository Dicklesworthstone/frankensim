//! Geometric support belongs to the deterministic light proposal, not radiance.

use super::*;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSceneError;
use crate::polarization::dielectric::scene::curved::camera::area_lights::AreaEmitterError;

#[test]
fn g1_closed_emitter_excluded_seeds_keep_native_hits_and_standalone_strictness() {
    with_cx(|cx| {
        // Four original faces, with the outward front face first. Its native
        // camera hit is visible, even when the authored light seed leaves the
        // opposite source half space and therefore has no proposal support.
        let lamp = GeometryInstance::try_new(9,
            hash_domain("closed-area-support", b"tetrahedron"),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-0.3, -0.3, -1.0], [0.3, -0.3, -1.0],
                    [0.0, 0.6, -1.0], [0.0, 0.0, -2.0]],
                vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]],
            )), RigidTransform::identity()).unwrap();
        let geometry = InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), lamp]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(cx, &geometry, &media,
            &[], &[], &reflectances, vec![diffuse(1, 0), opaque(9)], 100_000).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 4).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, -1.0))];
        let bindings = [AreaLightBinding {
            source: &source, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let samples = [(4.0 / 9.0) * AREA / source.total_area_m2(), 0.5];
        let selected = source.sample_area(samples).unwrap();
        assert_eq!(selected.triangle_index(), 0);
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
            settings().surfaces, |_| [0.0; 3]).unwrap();
        let result = owner.sample_ray_with_samples(cx, camera_ray(), 550.0, 0, settings(),
            &mut workspace, |_| [0.0; 3], |event, binding| {
                assert_eq!((event, binding), (0, 0)); samples
            }).unwrap();
        assert_eq!(result.area_radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(result.images(), 0);
        near(result.radiance_sample().radiance_w_m2_sr_nm(),
            RHO * LE * (-0.6_f64).exp(), 3e-6);
        assert_eq!(result.radiance_sample().path(), baseline.path());
        assert_eq!(result.radiance_sample().radiance_w_m2_sr_nm(),
            baseline.radiance_w_m2_sr_nm());
        assert!(result.trace_evaluations() >= 3);
        assert!(result.field_evaluations() > 0 && result.query_work() > 0);

        // The standalone irradiance integral has no complementary continuation.
        // It keeps the explicit unresolved-seed refusal for this same sample.
        let mut images = BundleImageWorkspace::try_new(1, 1).unwrap();
        let receiver = BundleConnectionTarget::try_diffuse_receiver(1,
            Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
        assert!(matches!(source.incident_irradiance(cx, selected, BoundarySide::Positive,
            0, Vec3::new(0.0, 1.0, 0.0), 550.0, receiver, settings().connection,
            &seeds, image_options(1), &mut images),
            Err(AreaEmitterError::UnresolvedSeed { seed: 0, error:
                BundleConnectionError::Scene(CurvedBundleSceneError::SourceSideExcluded { .. }),
            })));
    });
}

#[test]
fn g1_occluded_area_draw_preserves_the_visible_emitter_hit_partition() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 0.0, 1.0, 1.0),
            panel(5, 0.075, 0.0, 0.03),
            panel(9, 0.0, -1.0, 0.3),
        ]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(cx, &geometry, &media,
            &[], &[], &reflectances, vec![diffuse(1, 0), opaque(5), opaque(9)], 100_000).unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        // The first sample is (0.15,0,-1). This seed hits the actual blocker at
        // (0.075,0,0). At the camera's different source hit (0,0,-1), the same
        // authored seed and its Newton refinements remain clear of the blocker.
        let seeds = [seed(Vec3::new(-0.15, 0.0, 2.0))];
        let bindings = [AreaLightBinding {
            source: &source, receiver_object_id: 1, source_side: BoundarySide::Positive,
            initial_medium: 0, source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds, image_options: image_options(1), position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = camera.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
            settings().surfaces, |_| [0.0; 3]).unwrap();
        let p_light = 4.0 / AREA;
        let p_bsdf = 1.0 / core::f64::consts::PI;
        let hit_weight = p_bsdf * p_bsdf / (p_light * p_light + p_bsdf * p_bsdf);
        let expected_hit = RHO * LE * (-0.6_f64).exp() * hit_weight;
        let mut hit_terms = [0.0; 2];
        for (index, samples) in [[121.0 / 144.0, 4.0 / 11.0], [4.0 / 9.0, 0.5]]
            .into_iter().enumerate()
        {
            let result = owner.sample_ray_with_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace, |_| [0.0; 3], |_, _| samples).unwrap();
            if index == 0 {
                assert_eq!(result.area_radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(result.images(), 0);
            } else {
                assert!(result.area_radiance_w_m2_sr_nm() > 0.0);
                assert_eq!(result.images(), 1);
            }
            hit_terms[index] = result.radiance_sample().radiance_w_m2_sr_nm()
                - result.area_radiance_w_m2_sr_nm();
            near(hit_terms[index], expected_hit, 4e-6);
            assert_eq!(result.radiance_sample().path(), baseline.path());
            assert_eq!(result.area_samples(), 1);
        }
        near(hit_terms[0], hit_terms[1], 4e-6);
    });
}

