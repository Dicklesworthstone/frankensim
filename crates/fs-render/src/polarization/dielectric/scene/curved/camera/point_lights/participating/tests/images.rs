use super::*;
use super::super::images::{
    CurvedRayleighImageLightCameraScene, PointLightImageGroup, PointLightImageReceiver,
};
use crate::dielectric::DielectricSurface;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
use fs_geom::fixtures::SphereChart;

#[test]
fn g1_ball_lens_camera_sums_distinct_images_once_and_accepts_exact_trace_budget() {
    with_cx(|cx| {
        // These three finite seeds resolve the three regular ball-lens images
        // independently checked by the native image-discovery regressions.
        let target = Point3::new(0.03, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                1,
                hash_domain("image-point-camera-test", b"unit-sphere"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
            GeometryInstance::try_new(
                9,
                hash_domain("image-point-camera-test", b"diffuse-receiver"),
                SharedGeometry::mesh(TriMesh::new(
                    vec![
                        [target.x - 0.5, -0.5, target.z],
                        [target.x + 0.5, -0.5, target.z],
                        [target.x, 0.5, target.z],
                    ],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let glass = |index, absorption| {
            DielectricGlass::new(
                CauchyIor::try_constant(index).unwrap(),
                BeerLambertAbsorption::try_constant(absorption).unwrap(),
                GlassProvenance::Custom,
            )
        };
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let reflectance = [SpectralReflectance::constant(REFLECTANCE).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectance,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 1,
                        reflectance: 0,
                    },
                },
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let volumes = [
            Some(CurvedRayleighMedium::try_new(&spectrum, 0.3, &phase).unwrap()),
            None,
        ];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -3.0]).unwrap(),
            100.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
        )
        .unwrap();
        let seed = |angle: f64| BundleImageSeed {
            initial_direction: Vec3::new(angle.sin(), 0.0, angle.cos()),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        };
        let seeds = [seed(-0.18), seed(0.065), seed(0.14), seed(-0.18)];
        let ray = Ray {
            origin: Point3::new(target.x, 0.0, 2.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut limits = settings();
        limits.point_lights.surfaces.max_distance_m = 8.0;
        limits.point_lights.surfaces.event.surface_tolerance_m = 1e-8;
        limits.point_lights.connection.scene.trace = limits.point_lights.surfaces;
        limits.point_lights.max_connections = 3;
        limits.point_lights.max_trace_evaluations = 512;
        limits.max_source_free_flight_draws = 0;
        limits.max_source_candidates = 0;

        // The diffuse continuation is oblique enough to miss the lens. Source
        // replay must not consume or alter these same camera variates.
        let singles = [seeds[0], seeds[1], seeds[2]].map(|seed| {
            let branches = [PointLightBranch {
                source: &source,
                receiver_object_id: 9,
                initial_medium: 1,
                initial_direction: seed.initial_direction,
                source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
                policy: seed.policy,
                position_tolerance_m: 2e-7,
            }];
            let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
            let single =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
            single
                .sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    1,
                    limits,
                    |_| [0.8; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source scattering must consume no tracking variates"),
                )
                .unwrap()
        });
        for single in &singles {
            assert!(single.fixture_radiance_w_m2_sr_nm() > 0.0);
            assert_eq!(single.connections(), 1);
            assert_eq!(single.volume_events(), 0);
            assert_eq!(single.free_flight_draws(), 0);
            assert_eq!(single.radiance_sample().path(), singles[0].radiance_sample().path());
        }
        let expected = singles
            .iter()
            .map(|sample| sample.fixture_radiance_w_m2_sr_nm())
            .sum::<f64>();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let owner =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let groups = [PointLightImageGroup {
            source: &source,
            receiver: PointLightImageReceiver::DiffuseObject(9),
            initial_medium: 1,
            source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 3,
                min_separation_rad: 0.01,
                duplicate_tolerance_rad: 1e-5,
                max_trace_evaluations: 512,
            },
            position_tolerance_m: 2e-7,
        }];
        let images = CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &groups).unwrap();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(3, 4).unwrap();
        let sample = |limits, workspace: &mut CurvedTransmittanceImageWorkspace| {
            images
                .sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    1,
                    limits,
                    workspace,
                    |_| [0.8; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source scattering must consume no tracking variates"),
                )
                .unwrap()
        };
        let combined = sample(limits, &mut workspace);
        assert_eq!(combined.connections(), 3);
        assert_eq!(combined.source_free_flight_draws(), 0);
        assert_eq!(combined.source_candidates(), 0);
        assert_eq!(combined.radiance_sample().path(), singles[0].radiance_sample().path());
        near(combined.fixture_radiance_w_m2_sr_nm(), expected);
        near(combined.radiance_sample().radiance_w_m2_sr_nm(), expected);
        // Each reference includes one uniform survival replay. The duplicate
        // contributes its geometric solve, but neither a replay nor power.
        let expected_traces = singles.iter().map(|s| s.trace_evaluations()).sum::<u32>()
            + singles[0].trace_evaluations() - 1;
        assert_eq!(combined.trace_evaluations(), expected_traces);

        // Cached contractions require no further traces after the native
        // search and its three replays exactly exhaust the shared budget.
        limits.point_lights.max_trace_evaluations = combined.trace_evaluations();
        let exact = sample(limits, &mut workspace);
        assert_eq!(exact, combined);
    });
}
