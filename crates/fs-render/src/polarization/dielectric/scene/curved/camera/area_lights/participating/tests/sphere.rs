//! Analytic emitting charts in the original refracting molecular camera.

use super::*;
use fs_geom::fixtures::SphereChart;

#[test]
fn g1_native_sphere_source_crosses_two_scattering_media_before_molecular_vertex() {
    with_cx(|cx| {
        let interface = GeometryInstance::try_new(
            1, hash_domain("molecular-sphere-source-interface", b"native"),
            SharedGeometry::mesh(TriMesh::new(
                vec![[1.0, -10.0, -10.0], [1.0, 10.0, -10.0], [1.0, 0.0, 20.0]],
                vec![[0, 1, 2]],
            )), RigidTransform::identity(),
        ).unwrap();
        let sphere = GeometryInstance::try_new(
            9, hash_domain("molecular-sphere-area-source", b"original-chart"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(2.0, 0.0, 0.5), radius: 0.5,
            }), RigidTransform::identity(),
        ).unwrap();
        let geometry = InstanceScene::try_new(vec![interface, sphere]).unwrap();
        let air = DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(), GlassProvenance::Custom);
        let solid = DielectricGlass::new(CauchyIor::try_constant(1.5).unwrap(),
            BeerLambertAbsorption::try_constant(0.3).unwrap(), GlassProvenance::Custom);
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&solid)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![
                SurfaceBinding { object_id: 1, interaction: Interaction::Dielectric {
                    outside: 1, inside: 0, surface: DielectricSurface::SMOOTH,
                }},
                SurfaceBinding { object_id: 9, interaction: Interaction::Opaque },
            ], 100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(LE).unwrap(), EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let selected = source.sample_area_with_context(cx, [0.5, 0.5]).unwrap();
        assert_eq!(selected.triangle_index(), None);
        near(source.total_area_m2(), core::f64::consts::PI);
        near(selected.pdf_area_per_m2(), 1.0 / core::f64::consts::PI);
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let air_scattering = MassScatteringSpectrum::constant(0.7).unwrap();
        let solid_scattering = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [
            Some(CurvedRayleighMedium::try_new(&air_scattering, 1.0, &phase).unwrap()),
            Some(CurvedRayleighMedium::try_new(&solid_scattering, 1.0, &phase).unwrap()),
        ];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let flight = |draw| {
            if draw == 0 { [1.0 - (-0.35_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0] }
            else { [1.0 - 1e-12, 0.5, 0.0, 0.0] }
        };
        let continuation = |event| {
            assert_eq!(event, 0, "only the smooth interface scatters at a surface");
            [0.0, 0.0, 0.999]
        };
        let baseline = molecular.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
            CurvedRayleighTraceSettings { surfaces: settings().area_lights.surfaces, max_volume_events: 8 },
            continuation, flight).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        for represented in [true, false] {
            let direction = if represented { -1.0 } else { 1.0 };
            let seeds = seeds(Vec3::new(direction, 0.0, 0.0));
            let bindings = [MolecularAreaLightBinding { initial_medium: 1, ..binding(&source, &seeds) }];
            let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
                cx, &molecular, &surface, &bindings).unwrap();
            let mut surface_draws = 0;
            let mut area_draws = 0;
            let result = owner.sample_ray_with_tracking_samples(cx, camera_ray(), 550.0, 0,
                settings(), &mut workspace,
                |event| { surface_draws += 1; continuation(event) }, flight,
                |_, _| panic!("no surface receiver"), |event, index| {
                    assert_eq!((event, index), (0, 0));
                    area_draws += 1;
                    [0.5, 0.5]
                }, |_, _, _| panic!("both source-medium survivals are analytic"),
            ).unwrap();
            let phase_pdf = 3.0 / (16.0 * core::f64::consts::PI);
            let transmission = 1.0 - ((1.5_f64 - 1.0) / (1.5 + 1.0)).powi(2);
            let eta2 = (1.0_f64 / 1.5).powi(2);
            // Half a metre in source glass, then one metre in air. The native
            // ray-transverse angular footprint follows differentiated Snell:
            // J = (0.5 + 1.5*1)^2. Receiving points have no surface cosine.
            let jacobian = (0.5_f64 + 1.5).powi(2);
            let light_pdf = eta2 * jacobian / core::f64::consts::PI;
            let (light_weight, hit_weight) = partition(light_pdf, transmission * phase_pdf);
            // Camera-prefix Beer depth .1, source Beer depth .35, and
            // independent source scattering depth .7*1 + .4*.5 = .9.
            // Survival is physical attenuation, absent from the partition.
            let raw_area = LE * core::f64::consts::PI * transmission / jacobian
                * phase_pdf * (-1.35_f64).exp();
            let raw_hit = LE * eta2 * (-0.45_f64).exp();
            let expected_area = if represented { raw_area * light_weight } else { 0.0 };
            let expected_hit = raw_hit * if represented { hit_weight } else { 1.0 };
            near(baseline.radiance_sample().radiance_w_m2_sr_nm(), raw_hit);
            near(result.area_radiance_w_m2_sr_nm(), expected_area);
            near(result.radiance_sample().radiance_w_m2_sr_nm(), expected_area + expected_hit);
            assert_eq!(result.radiance_sample().path(), baseline.radiance_sample().path());
            assert_eq!(result.free_flight_draws(), baseline.free_flight_draws());
            assert_eq!(result.volume_events(), 1);
            assert_eq!((surface_draws, area_draws, result.area_samples()), (1, 1, 1));
            assert_eq!(result.images(), u32::from(represented));
            assert_eq!((result.source_free_flight_draws(), result.source_candidates()), (0, 0));
        }
    });
}
