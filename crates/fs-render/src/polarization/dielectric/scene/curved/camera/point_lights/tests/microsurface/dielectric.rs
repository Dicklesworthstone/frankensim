//! Native-scene composition tests for full glass and calibrated point sources.
//! The numerical phase/height law has independent tests in the core evaluator.

use super::*;
use crate::dielectric::DielectricSurface;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};

#[test]
fn g1_full_glass_point_connections_keep_both_source_media_beer_and_single_eta_conversion() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, 0.0, 20.0),
            plane(8, -2.0, 0.0, 20.0),
            plane(9, 2.0, 0.0, 20.0),
        ])
        .unwrap();
        let materials = [glass(1.0, 0.12), glass(1.5, 0.24)];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let surface =
            DielectricSurface::try_anisotropic_rough(0.7, 0.9, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface,
                    },
                },
                opaque(8),
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let limits = SmithWalkSettings::try_new(4).unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.37).unwrap();
        for camera_medium in [0, 1] {
            let sign = if camera_medium == 0 { 1.0 } else { -1.0 };
            let ray = Ray {
                origin: Point3::new(0.0, 0.0, 0.5 * sign),
                dir: Vec3::new(0.0, 0.0, -sign),
            };
            let launch = analyzer.launch(&camera, ray, 550.0, camera_medium).unwrap();
            for source_medium in [0, 1] {
                let source_sign = if source_medium == 0 { 1.0 } else { -1.0 };
                let lamp = source(Point3::new(0.0, 0.0, source_sign));
                let direction = Vec3::new(0.0, 0.0, -source_sign);
                let mut selected = branch(&lamp, 1, direction);
                selected.initial_medium = source_medium;
                let branches = [selected];
                let point = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
                let full = SmithDielectricPointLightCameraScene::new(&point, limits);
                let draw = [0.5, 0.0, 0.0, 0.99];
                let incoming = PolarizedBeam::try_new(
                    Point3::new(0.0, 0.0, 0.0),
                    launch.frame(),
                    550.0,
                    materials[camera_medium],
                    launch.power(),
                )
                .unwrap();
                let factor = incoming
                    .evaluate_smith_dielectric_connection(
                        cx,
                        Vec3::new(0.0, 0.0, 1.0),
                        materials[1 - camera_medium],
                        surface,
                        PolarizationFrame::try_new(direction, selected.source_transverse_axis)
                            .unwrap(),
                        Stokes::unpolarized(1.0),
                        limits,
                        |_| draw,
                    )
                    .unwrap()
                    .factor_per_sr();
                let mut continuation_ids = Vec::new();
                let mut evaluation_ids = Vec::new();
                let sample = full
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray,
                        550.0,
                        camera_medium,
                        analyzer,
                        settings(),
                        |_| panic!("full glass consumed an ordinary draw"),
                        |_, _| panic!("full glass consumed an opaque continuation"),
                        |event, facet| {
                            continuation_ids.push((event, facet));
                            draw
                        },
                        |_, _, _| panic!("full glass consumed an opaque evaluation"),
                        |event, binding, facet| {
                            evaluation_ids.push((event, binding, facet));
                            draw
                        },
                    )
                    .unwrap();
                let absorption: [f64; 2] = [0.12, 0.24];
                let expected = launch.power().i
                    * factor
                    * intensity(&lamp)
                    * (-0.5 * absorption[camera_medium] - absorption[source_medium]).exp();
                assert!(expected > 0.0);
                near(sample.fixture_radiance_w_m2_sr_nm(), expected);
                near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
                assert_eq!(sample.connections(), 1);
                assert_eq!(
                    sample.radiance_sample().path().medium_index(),
                    1 - camera_medium
                );
                assert_eq!(continuation_ids, [(0, 0), (0, 1)]);
                assert_eq!(evaluation_ids, [(0, 0, 0), (0, 0, 1)]);
            }
        }
    });
}

#[test]
fn g0_full_glass_point_refusal_never_returns_partial_direct_light() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, 0.0, 20.0),
            plane(8, -2.0, 0.0, 20.0),
            plane(9, 2.0, 0.0, 20.0),
        ])
        .unwrap();
        let materials = [glass(1.0, 0.0), glass(1.5, 0.0)];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::try_rough(1.0).unwrap(),
                    },
                },
                opaque(8),
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(Point3::new(0.0, 0.0, -1.0));
        let mut selected = branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0));
        selected.initial_medium = 1;
        let branches = [selected];
        let point = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let full = SmithDielectricPointLightCameraScene::new(
            &point,
            SmithWalkSettings::try_new(16).unwrap(),
        );
        let mut continued = 0;
        let result = full.sample_ray_with_samples(
            cx,
            Ray {
                origin: Point3::new(0.0, 0.0, 0.5),
                dir: Vec3::new(0.0, 0.0, -1.0),
            },
            550.0,
            0,
            settings(),
            |_| [0.5; 3],
            |_, _| [0.5; 3],
            |_, _| {
                continued += 1;
                [0.5, 0.0, 0.0, 0.99]
            },
            |_, _, _| panic!("opaque evaluation"),
            |_, _, _| {
                gate.request();
                [0.5, 0.0, 0.0, 0.99]
            },
        );
        assert!(result.is_err());
        assert_eq!(
            continued, 0,
            "a refused source estimate must not advance continuation"
        );
    });
}
