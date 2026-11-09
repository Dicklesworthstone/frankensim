use super::*;
use crate::dielectric::{DielectricEvent, DielectricSurface};
use crate::polarization::dielectric::transport::microsurface::dielectric::SmithDielectricError;

fn draws(path: u32, facet: u32) -> [f64; 4] {
    let block = fs_rand::philox::philox4x32_10([path, facet, 31, 0], [0x6173_676c, 0x6469_656c]);
    core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
}
fn glass_bindings(surface: DielectricSurface) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface,
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

#[test]
fn full_glass_camera_preserves_final_medium_beer_etendue_and_geometric_distance() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0), plane(8, -2.0), plane(9, 2.0)]).unwrap();
        let air = glass(1.0, 0.2);
        let material = glass(1.5, 0.3);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let surface =
            DielectricSurface::try_anisotropic_rough(0.9, 0.6, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, glass_bindings(surface), 100_000)
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
        let limits = SmithWalkSettings::try_new(512).unwrap();
        let smith = SmithDielectricCameraScene::new(&camera, limits);
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        let incoming = analyzer.launch(&camera, ray(), 550.0, 0).unwrap();
        let mut paths = [None, None];
        for path in 0..16_384 {
            let sample = incoming
                .sample_smith_dielectric(
                    cx,
                    Vec3::new(0.0, 0.0, 1.0),
                    material,
                    surface,
                    limits,
                    |facet| draws(path, facet),
                )
                .unwrap();
            let transmitted = sample.event() == DielectricEvent::Transmission;
            if sample.facet_events() >= 2
                && sample.beam().ray().dir.z.abs() > 0.4
                && (transmitted || sample.transmissions() >= 2)
            {
                let index = usize::from(transmitted);
                paths[index].get_or_insert((path, sample));
            }
            if paths.iter().all(Option::is_some) {
                break;
            }
        }
        for (path, expected) in
            paths.map(|p| p.expect("multiple-facet reflection and transmission"))
        {
            let mut called = Vec::new();
            let actual = smith
                .sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| panic!("full rough glass borrowed an ordinary surface variate"),
                    |_, _| panic!("glass borrowed an opaque variate"),
                    |surface, facet| {
                        called.push((surface, facet));
                        draws(path, facet)
                    },
                )
                .unwrap();
            assert_eq!(called.len(), expected.facet_events() as usize + 1);
            assert!(called.iter().all(|(surface, _)| *surface == 0));
            let transmitted = expected.event() == DielectricEvent::Transmission;
            let (index, absorption, emission, object) = if transmitted {
                (1.5, 0.3, 7.0, 8)
            } else {
                (1.0, 0.2, 4.0, 9)
            };
            let tail = 2.0 / expected.beam().ray().dir.z.abs();
            let transmission = (-0.2_f64 - absorption * tail).exp();
            let power = expected.beam().power().i * transmission;
            near(
                actual.radiance_w_m2_sr_nm(),
                emission * power / (index * index),
                5e-6,
            );
            near(actual.path().beam().power().i, power, 4e-6);
            near(actual.path().traveled_m(), 1.0 + tail, 4e-6);
            near(
                actual.path().optical_path_length_m(),
                1.0 + index * tail,
                6e-6,
            );
            assert_eq!(actual.path().medium_index(), usize::from(transmitted));
            assert_eq!(
                actual.endpoint(),
                RadianceEndpoint::Surface { object_id: object }
            );
            assert_eq!(actual.path().surface_events(), 2);
            let bounded =
                SmithDielectricCameraScene::new(&camera, SmithWalkSettings::try_new(1).unwrap());
            assert_eq!(
                bounded.sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| panic!("ordinary draw"),
                    |_, _| panic!("opaque draw"),
                    |_, facet| draws(path, facet)
                ),
                Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                    SceneError::SmithDielectric(SmithDielectricError::FacetLimit {
                        facet_events: 1
                    })
                )))
            );
        }
    });
}

#[test]
fn full_glass_opt_in_preserves_smooth_camera_and_four_wavelength_raster() {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0), plane(8, -2.0), plane(9, 2.0)]).unwrap();
        let air = glass(1.0, 0.0);
        let material = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            glass_bindings(DielectricSurface::SMOOTH),
            100_000,
        )
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
        let smith =
            SmithDielectricCameraScene::new(&camera, SmithWalkSettings::try_new(8).unwrap());
        for id in 0..16 {
            assert_eq!(
                smith
                    .sample_ray(cx, ray(), 550.0, 0, settings(), 61, id)
                    .unwrap(),
                camera
                    .sample_ray(cx, ray(), 550.0, 0, settings(), 61, id)
                    .unwrap()
            );
        }
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.1).unwrap(),
            1.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let analyzer = CameraAnalyzer::try_circular(Vec3::new(1.0, 0.0, 0.0), true).unwrap();
        assert_eq!(
            raster
                .sample_pixel_smith_dielectric(cx, &smith, 0, None, settings(), 37, 4, 9)
                .unwrap(),
            raster
                .sample_pixel(cx, &camera, 0, settings(), 37, 4, 9)
                .unwrap()
        );
        assert_eq!(
            raster
                .sample_pixel_smith_dielectric(cx, &smith, 0, Some(analyzer), settings(), 37, 4, 9)
                .unwrap(),
            raster
                .sample_pixel_analyzed(cx, &camera, 0, analyzer, settings(), 37, 4, 9)
                .unwrap()
        );
        gate.request();
        assert_eq!(
            smith.sample_ray(cx, ray(), 550.0, 0, settings(), 61, 9),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        );
    });
}
