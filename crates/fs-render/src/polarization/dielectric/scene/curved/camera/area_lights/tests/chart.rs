use super::*;

fn sphere_source(placement: RigidTransform) -> GeometryInstance {
    GeometryInstance::try_new(2, hash_domain("native-area-emitter-test", b"sphere-source"),
        SharedGeometry::chart(SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 0.5 }),
        placement).unwrap()
}

#[test]
fn g1_chart_source_preserves_native_area_placement_and_exact_preparation_owner() {
    with_cx(|cx| {
        let half = 0.5_f64.sqrt();
        let placement = RigidTransform::try_new([0.0, half, 0.0, half], [1.0, 2.0, 3.0]).unwrap();
        let geometry = InstanceScene::try_new(vec![sphere_source(placement)]).unwrap();
        let material = glass(1.3, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }], 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(2,
            super::super::super::RadianceSpectrum::try_constant(4.0).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let limits = AreaEmitterLimits { max_triangles: 0, max_preparation_work: 1, max_lookup_work: 1 };
        let source = CurvedAreaEmitter::try_new_with_limits(cx, &scene, &emitters[0], limits).unwrap();
        let duplicate = CurvedAreaEmitter::from_camera(cx, &camera, 2, 0).unwrap();
        assert!(source.support_overlaps(&duplicate));
        assert!(core::ptr::eq(source.source_instance(), &geometry.instances()[0]));
        assert!(core::ptr::eq(source.source_emitter(), &emitters[0]));
        near(source.total_area_m2(), core::f64::consts::PI, 1e-14);
        assert!(matches!(source.sample_area([0.25, 0.5]), Err(AreaEmitterError::ContextRequired)));
        assert!(matches!(source.sample_area_with_context(cx, [1.0, 0.5]),
            Err(AreaEmitterError::InvalidSample)));
        assert!(matches!(CurvedAreaEmitter::try_new_faces(cx, &scene, &emitters[0], &[0], limits),
            Err(AreaEmitterError::UnsupportedFaceSubset)));
        assert!(matches!(CurvedAreaEmitter::try_new_with_limits(cx, &scene, &emitters[0],
            AreaEmitterLimits { max_preparation_work: 0, ..limits }),
            Err(AreaEmitterError::SourceWorkLimit { limit: 0 })));

        let center = placement.transform_point(Point3::new(0.0, 0.0, 0.0));
        for draw in [[0.0, 0.3], [0.25, 0.125], [0.75, 0.6]] {
            let sample = source.sample_area_with_context(cx, draw).unwrap();
            assert_eq!(sample.object_id(), 2);
            assert_eq!(sample.triangle_index(), None);
            near(sample.pdf_area_per_m2(), 1.0 / core::f64::consts::PI, 1e-14);
            let radius = sample.position_m().delta_from(center);
            near(norm(radius), 0.5, 1e-14);
            near(radius.dot(sample.normal()), 0.5, 1e-14);
            let seeds = [BundleImageSeed { initial_direction: sample.normal(),
                policy: BundleBranchPolicy::Sequence(&[]) }];
            assert!(matches!(duplicate.launch(sample, BoundarySide::Positive, 0,
                Vec3::new(0.0, 1.0, 0.0), 550.0, &seeds), Err(AreaEmitterError::SampleOwnerMismatch)));
        }
        // The sphere's native north pole rotates into +X; its geometric
        // normal and physical area follow the original proper-rigid owner.
        let pole = source.sample_area_with_context(cx, [0.0, 0.0]).unwrap();
        near(pole.position_m().x, 1.5, 1e-14);
        near(pole.position_m().y, 2.0, 1e-14);
        near(pole.position_m().z, 3.0, 1e-14);
        near(pole.normal().x, 1.0, 1e-14);
        assert!(pole.normal().y.abs() < 1e-14 && pole.normal().z.abs() < 1e-14);
    });
}

#[test]
fn g1_native_spherical_emitter_irradiance_uses_area_density_and_original_departure() {
    with_cx(|cx| {
        let point = Point3::new(0.0, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![sphere_source(RigidTransform::identity()),
            detector(point)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }, detector_binding()],
            100_000).unwrap();
        let emitter = SurfaceEmitter::new(2,
            super::super::super::RadianceSpectrum::try_constant(4.0).unwrap(), EmissionSides::Front);
        let source = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 0).unwrap();
        let sample = source.sample_area_with_context(cx, [0.0, 0.0]).unwrap();
        assert_eq!(sample.triangle_index(), None);
        assert_eq!(sample.position_m(), Point3::new(0.0, 0.0, 0.5));
        let target = BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
        let seeds = [seed(0.0)];
        let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
        let result = source.incident_irradiance(cx, sample, BoundarySide::Positive, 0,
            Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(), &seeds,
            BundleImageOptions { max_images: 1, ..BundleImageOptions::default() }, &mut workspace).unwrap();
        // Le=4, source cosine=1, native sphere area=pi, distance=2 m,
        // absorption=0.2/m. Uniform host n=1.3 adds no extra source-leg factor.
        near(result.spectral_irradiance_w_m2_nm(), core::f64::consts::PI * (-0.4_f64).exp(), 2e-6);
        assert_eq!(result.image_count(), 1);
    });
}
