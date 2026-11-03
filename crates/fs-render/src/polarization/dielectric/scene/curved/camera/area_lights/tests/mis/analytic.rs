//! Actual chart contacts retain their native area measure under rigid placement.

use super::*;
use crate::dielectric::DielectricSurface;
use crate::polarization::dielectric::scene::curved::camera::{CurvedCameraScene, RadianceSpectrum};
use fs_geom::fixtures::SphereChart;

#[test]
fn g1_placed_sphere_emitter_hits_recover_native_density_and_reciprocal_mis() {
    with_cx(|cx| {
        let radius = 0.75;
        let area = 4.0 * core::f64::consts::PI * radius * radius;
        let angle = 0.37_f64;
        let placement = RigidTransform::try_new([0.0, angle.sin(), 0.0, angle.cos()],
            [10.0, -3.0, 2.0]).unwrap();
        let geometry = InstanceScene::try_new(vec![GeometryInstance::try_new(2,
            hash_domain("chart-area-mis", b"original-sphere"),
            SharedGeometry::chart(SphereChart { center: Point3::new(0.25, -0.3, 0.4), radius }),
            placement).unwrap()]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }],
            100_000).unwrap();
        let emitters = [SurfaceEmitter::new(2, RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 2, 0).unwrap();
        let other_preparation = CurvedAreaEmitter::from_camera(cx, &camera, 2, 0).unwrap();
        let no_lookup = CurvedAreaEmitter::try_new_with_limits(cx, &scene, &emitters[0],
            AreaEmitterLimits { max_lookup_work: 0, ..AreaEmitterLimits::default() }).unwrap();
        near(source.total_area_m2(), area, 1e-14);
        assert!(matches!(source.source_instance().geometry(), SharedGeometry::Chart(_)));
        assert!(source.support_overlaps(&other_preparation));

        for draws in [[0.0, 0.0], [0.25, 0.125], [0.5, 0.25], [0.9, 0.7]] {
            let sampled = source.sample_area_with_context(cx, draws).unwrap();
            assert_eq!(sampled.triangle_index(), None);
            let normal = sampled.normal();
            let ray = Ray { origin: sampled.position_m().offset(normal.scale(1.25)),
                dir: normal.scale(-1.0) };
            let hit = camera.sample_ray_with_samples(cx, ray, 550.0, 0, options().scene.trace,
                |_| panic!("opaque original sphere consumes no surface draw")).unwrap();
            let contact = hit.path().last_contact().unwrap();
            let position = hit.path().beam().ray().origin;
            let recovered = source.sample_at_contact_if_supported(cx, contact, position).unwrap().unwrap();
            assert_eq!(recovered.triangle_index(), None);
            assert_eq!(recovered.position_m(), position);
            near(recovered.pdf_area_per_m2(), 1.0 / area, 1e-14);
            near(recovered.normal().dot(normal), 1.0, 1e-12);
            assert!(matches!(no_lookup.sample_at_contact_if_supported(cx, contact, position),
                Err(AreaEmitterError::SourceWorkLimit { limit: 0 })));

            let reference = if normal.x.abs() < 0.9 { Vec3::new(1.0, 0.0, 0.0) }
                else { Vec3::new(0.0, 1.0, 0.0) };
            let target = BundleConnectionTarget::try_volume_receiver(0,
                position.offset(normal.scale(2.0)), normal.scale(-1.0), 2e-7).unwrap();
            let connected = scene.connect_bundle_from_surface_with_policy(cx,
                BundleConnectionLaunch {
                    ray: Ray { origin: position, dir: normal },
                    transverse_axis: cross(normal, reference), wavelength_nm: 550.0,
                    normalized_source_stokes: Stokes::unpolarized(1.0), initial_medium: 0,
                }, BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(), target,
                options(), BundleBranchPolicy::Sequence(&[])).unwrap();
            let weights = source.volume_mis_weights(recovered, &connected, target,
                options().projection, 0.2).unwrap();
            near(weights.light_log_pdf_solid_angle().exp(), 4.0 / area, 5e-6);
            near(weights.light_weight() + weights.continuation_weight(), 1.0, 1e-15);
            assert!(matches!(other_preparation.volume_mis_weights(recovered, &connected, target,
                options().projection, 0.2), Err(AreaEmitterError::SampleOwnerMismatch)));
        }
    });
}

#[test]
fn g1_emitting_dielectric_sphere_keeps_each_original_side_and_its_radiance_calibration() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![GeometryInstance::try_new(2,
            hash_domain("chart-area-mis", b"emitting-glass-sphere"),
            SharedGeometry::chart(SphereChart { center: Point3::new(0.0, 0.0, 0.0), radius: 0.5 }),
            RigidTransform::identity()).unwrap()]).unwrap();
        let air = glass(1.0, 0.0);
        let solid = glass(1.5, 0.2);
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&solid)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 2, interaction: Interaction::Dielectric {
                outside: 0, inside: 1, surface: DielectricSurface::SMOOTH,
            }}], 100_000).unwrap();
        let emitter = SurfaceEmitter::new(2, RadianceSpectrum::try_constant(4.0).unwrap(), EmissionSides::Both);
        let source = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 0).unwrap();
        let sample = source.sample_area_with_context(cx, [0.0, 0.0]).unwrap();
        let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
        for (side, medium, z, direction, optical_depth) in [
            (BoundarySide::Positive, 0, 0.9, 1.0, 0.0_f64),
            (BoundarySide::Negative, 1, 0.1, -1.0, 0.08),
        ] {
            let target = BundleConnectionTarget::try_volume_receiver(medium,
                Point3::new(0.0, 0.0, z), Vec3::new(0.0, 0.0, -direction), 2e-7).unwrap();
            let seeds = [BundleImageSeed { initial_direction: Vec3::new(0.0, 0.0, direction),
                policy: BundleBranchPolicy::Sequence(&[]) }];
            let result = source.incident_irradiance(cx, sample, side, medium,
                Vec3::new(0.0, 1.0, 0.0), 550.0, target, options(), &seeds,
                BundleImageOptions { max_images: 1, ..BundleImageOptions::default() },
                &mut workspace).unwrap();
            // Calibrated Le=4 is emitted into the selected adjacent medium.
            // Area=pi, distance=.4 and cosine=1 give 25*pi before Beer loss.
            // The source itself adds neither a Fresnel interaction, a side
            // selection factor nor an index conversion within the same medium.
            near(result.spectral_irradiance_w_m2_nm(),
                25.0 * core::f64::consts::PI * (-optical_depth).exp(), 5e-6);
            assert_eq!(result.image_count(), 1);
        }
    });
}
