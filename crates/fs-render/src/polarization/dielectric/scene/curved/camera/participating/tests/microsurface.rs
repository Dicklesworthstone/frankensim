use super::*;
use crate::dielectric::{DielectricEvent, DielectricSurface};
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithDielectricCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, IsothermalVolumeEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::microsurface::SmithDielectricRayleighScene;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::polarization::dielectric::transport::microsurface::dielectric::SmithDielectricError;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;

fn wide_plane(id: u64, z: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("org.frankensim.test.smith-volume", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-100.0, -100.0, z], [100.0, -100.0, z], [0.0, 100.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn glass_bindings() -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface: rough_glass(),
            },
        },
        SurfaceBinding {
            object_id: 8,
            interaction: Interaction::Opaque,
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Opaque,
        },
    ]
}

fn rough_glass() -> DielectricSurface {
    DielectricSurface::try_anisotropic_rough(0.9, 0.6, Vec3::new(1.0, 0.0, 0.0)).unwrap()
}

fn glass_ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 1.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    }
}

fn limits() -> CurvedRayleighTraceSettings {
    let mut limits = settings();
    limits.surfaces.max_distance_m = 64.0;
    limits
}

fn facets(path: u32, facet: u32) -> [f64; 4] {
    fs_rand::philox::philox4x32_10([path, facet, 31, 0], [0x6173_676c, 0x6469_656c])
        .map(|v| (f64::from(v) + 0.5) / 4_294_967_296.0)
}

#[test]
fn g1_zero_scattering_retains_complete_glass_camera_bits_and_all_stream_addresses() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            wide_plane(1, 0.0),
            wide_plane(8, -2.0),
            wide_plane(9, 3.0),
        ])
        .unwrap();
        let air = glass(1.0, 0.2);
        let material = glass(1.5, 0.3);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, glass_bindings(), 100_000)
                .unwrap();
        let sources = [
            SurfaceEmitter::new(
                8,
                RadianceSpectrum::try_constant(7.0).unwrap(),
                EmissionSides::Both,
            ),
            SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(4.0).unwrap(),
                EmissionSides::Both,
            ),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let laws = [None, None];
        let transport = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let microscopic = SmithWalkSettings::try_new(512).unwrap();
        let full = SmithDielectricCameraScene::new(&camera, microscopic);
        let volume = SmithDielectricRayleighCameraScene::new(&participating, microscopic);
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        let mut completed = 0;
        for path in 0..32 {
            let plain = full.sample_ray_analyzed(
                cx,
                glass_ray(),
                550.0,
                0,
                analyzer,
                limits().surfaces,
                83,
                path,
            );
            let scattering =
                volume.sample_ray_analyzed(cx, glass_ray(), 550.0, 0, analyzer, limits(), 83, path);
            assert_eq!(plain, scattering.clone().map(|s| *s.radiance_sample()));
            if let Ok(sample) = scattering {
                completed += 1;
                assert_eq!((sample.volume_events(), sample.free_flight_draws()), (0, 0));
            }
        }
        assert!(completed >= 24);
        let plain = full.sample_ray_with_samples(
            cx,
            glass_ray(),
            550.0,
            0,
            limits().surfaces,
            |_| panic!("glass does not consume ordinary draws"),
            |_, _| panic!("glass does not consume opaque draws"),
            |_, facet| facets(8, facet),
        );
        let scattering = volume.sample_ray_with_tracking_samples(
            cx,
            glass_ray(),
            550.0,
            0,
            limits(),
            |_| panic!("glass does not consume ordinary draws"),
            |_| panic!("zero scattering does not consume flights"),
            |_, _| panic!("glass does not consume opaque draws"),
            |_, facet| facets(8, facet),
        );
        assert_eq!(plain, scattering.map(|s| *s.radiance_sample()));
    });
}

#[test]
fn g1_physical_collision_and_multifacet_glass_preserve_medium_beer_lte_and_limits() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            wide_plane(1, 0.0),
            wide_plane(8, -2.0),
            wide_plane(9, 3.0),
        ])
        .unwrap();
        let air = glass(1.0, 0.2);
        let material = glass(1.5, 0.3);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, glass_bindings(), 100_000)
                .unwrap();
        let sources = [
            SurfaceEmitter::new(
                8,
                RadianceSpectrum::try_constant(7.0).unwrap(),
                EmissionSides::Both,
            ),
            SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(4.0).unwrap(),
                EmissionSides::Both,
            ),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let hot_camera = CurvedCameraScene::try_new(cx, &native, &sources, None)
            .unwrap()
            .with_isothermal_volume_emitters(
                cx,
                &[
                    IsothermalVolumeEmitter::try_new(0, 2300.0, ContentHash([118; 32])).unwrap(),
                    IsothermalVolumeEmitter::try_new(1, 2100.0, ContentHash([119; 32])).unwrap(),
                ],
            )
            .unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.5).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [
            Some(CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap()),
            None,
        ];
        let transport = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let participating_hot =
            CurvedRayleighCameraScene::try_new(cx, &hot_camera, &transport).unwrap();
        let microscopic = SmithWalkSettings::try_new(512).unwrap();
        let full = SmithDielectricRayleighCameraScene::new(&participating, microscopic);
        let hot = SmithDielectricRayleighCameraScene::new(&participating_hot, microscopic);
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        let launch = analyzer.launch(&camera, glass_ray(), 550.0, 0).unwrap();
        let collision = phase
            .sample_polarized(launch.advance(0.25).unwrap(), [0.9, 0.2])
            .unwrap()
            .beam;
        assert!(collision.ray().dir.z < -0.25);
        let middle = -collision.ray().origin.z / collision.ray().dir.z;
        let arrival = collision.advance(middle).unwrap();
        let mut choices = [None, None];
        for path in 0..16_384 {
            let sample = arrival
                .sample_smith_dielectric(
                    cx,
                    Vec3::new(0.0, 0.0, 1.0),
                    material,
                    rough_glass(),
                    microscopic,
                    |facet| facets(path, facet),
                )
                .unwrap();
            if sample.facet_events() >= 2 && sample.beam().ray().dir.z.abs() > 0.4 {
                choices[usize::from(sample.event() == DielectricEvent::Transmission)]
                    .get_or_insert((path, sample));
            }
            if choices.iter().all(Option::is_some) {
                break;
            }
        }
        let flight = |ordinal| {
            if ordinal == 0 {
                [1.0 - (-0.125_f64).exp(), 0.9, 0.2, 0.0]
            } else {
                [1.0 - (-10.0_f64).exp(), 0.5, 0.0, 0.0]
            }
        };
        for (path, expected) in choices.map(|p| p.expect("both full glass branches")) {
            let transmitted = expected.event() == DielectricEvent::Transmission;
            let (index, absorption, light, distance, temperature) = if transmitted {
                (1.5, 0.3, 7.0, 2.0, 2100.0)
            } else {
                (1.0, 0.2, 4.0, 3.0, 2300.0)
            };
            let tail = distance / expected.beam().ray().dir.z.abs();
            let endpoint = expected.beam().advance(tail).unwrap();
            let mut addresses = Vec::new();
            let actual = full
                .sample_ray_analyzed_with_tracking_samples(
                    cx,
                    glass_ray(),
                    550.0,
                    0,
                    analyzer,
                    limits(),
                    |_| panic!("no ordinary response"),
                    flight,
                    |_, _| panic!("no opaque response"),
                    |surface, facet| {
                        addresses.push((surface, facet));
                        facets(path, facet)
                    },
                )
                .unwrap();
            assert!(addresses.iter().all(|(surface, _)| *surface == 0));
            assert_eq!(addresses.len(), expected.facet_events() as usize + 1);
            assert_eq!(actual.volume_events(), 1);
            assert_eq!(actual.free_flight_draws(), if transmitted { 2 } else { 3 });
            assert_eq!(
                actual.radiance_sample().path().medium_index(),
                usize::from(transmitted)
            );
            close(
                actual.radiance_sample().radiance_w_m2_sr_nm(),
                light * endpoint.power().i / (index * index),
                3e-6,
            );
            close(
                actual.radiance_sample().path().traveled_m(),
                0.25 + middle + tail,
                3e-6,
            );
            close(
                actual.radiance_sample().path().optical_path_length_m(),
                0.25 + middle + index * tail,
                5e-6,
            );
            let hot_actual = hot
                .sample_ray_analyzed_with_tracking_samples(
                    cx,
                    glass_ray(),
                    550.0,
                    0,
                    analyzer,
                    limits(),
                    |_| panic!("no ordinary response"),
                    flight,
                    |_, _| panic!("no opaque response"),
                    |_, facet| facets(path, facet),
                )
                .unwrap();
            let b0 = planck_radiance_w_m2_sr_nm(2300.0, 550.0).unwrap();
            let b1 = planck_radiance_w_m2_sr_nm(temperature, 550.0).unwrap();
            let thermal = b0
                * (launch.power().i * -(-0.05_f64).exp_m1()
                    + collision.power().i * -(-0.2 * middle).exp_m1())
                + b1 * expected.beam().power().i * -(-absorption * tail).exp_m1();
            close(
                hot_actual.radiance_sample().radiance_w_m2_sr_nm()
                    - actual.radiance_sample().radiance_w_m2_sr_nm(),
                thermal,
                2e-7,
            );
            let forward = SmithDielectricRayleighScene::new(&transport, microscopic)
                .trace_with_tracking_samples(
                    cx,
                    launch,
                    0,
                    limits(),
                    |_| panic!("no ordinary response"),
                    flight,
                    |_, _| panic!("no opaque response"),
                    |_, facet| facets(path, facet),
                )
                .unwrap();
            assert_eq!(forward.path(), actual.radiance_sample().path());
            let bounded = SmithDielectricRayleighCameraScene::new(
                &participating,
                SmithWalkSettings::try_new(1).unwrap(),
            );
            assert!(matches!(
                bounded.sample_ray_analyzed_with_tracking_samples(
                    cx,
                    glass_ray(),
                    550.0,
                    0,
                    analyzer,
                    limits(),
                    |_| panic!("no ordinary response"),
                    flight,
                    |_, _| panic!("no opaque response"),
                    |_, facet| facets(path, facet),
                ),
                Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                    SceneError::SmithDielectric(SmithDielectricError::FacetLimit { .. })
                )))
            ));
        }
        let mut capped = limits();
        capped.max_volume_events = 1;
        assert_eq!(
            full.sample_ray_with_tracking_samples(
                cx,
                glass_ray(),
                550.0,
                0,
                capped,
                |_| panic!("volume cap first"),
                flight,
                |_, _| panic!("volume cap first"),
                |_, _| panic!("volume cap first"),
            ),
            Err(CameraRadianceError::Unresolved(
                PathTermination::VolumeLimit
            ))
        );
        gate.request();
        assert!(
            full.sample_ray_with_tracking_samples(
                cx,
                glass_ray(),
                550.0,
                0,
                limits(),
                |_| panic!("cancelled"),
                |_| panic!("cancelled"),
                |_, _| panic!("cancelled"),
                |_, _| panic!("cancelled"),
            )
            .is_err()
        );
    });
}
