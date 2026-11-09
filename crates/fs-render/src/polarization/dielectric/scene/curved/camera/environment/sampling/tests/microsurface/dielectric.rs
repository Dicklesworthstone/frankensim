//! Full native scenes for the shared glass/environment composition.
//! These exercise the production geometry, glass walk, source and radiometry.

use super::*;
use crate::conductor::{CauchyCoatingLayer, CoherentPeriodicCoating};
use crate::dielectric::DielectricCoating;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithDielectricCameraScene;

fn glass_draw(transmitted: bool) -> [f64; 4] {
    [0.5, 0.0, 0.0, if transmitted { 0.999 } else { 0.001 }]
}

#[test]
fn full_glass_environment_connections_partition_both_sides_and_keep_native_camera_power() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let glass = material(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&glass),
        ];
        let surface = DielectricSurface::try_rough(0.7).unwrap();
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 1,
                    surface,
                },
            }],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, Some(source)).unwrap();
        let environment = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
        let limit = SmithWalkSettings::try_new(8).unwrap();
        let owner = SmithDielectricEnvironmentCameraScene::new(&environment, limit);
        for transmitted in [false, true] {
            let draw = [0.5, if transmitted { 0.25 } else { 0.75 }];
            let mut continuation = Vec::new();
            let mut connection = Vec::new();
            let observed = owner
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("ordinary glass draw"),
                    |event| {
                        assert_eq!(event, 0);
                        draw
                    },
                    |_, _| panic!("opaque continuation"),
                    |_, _| panic!("opaque connection"),
                    |event, facet| {
                        continuation.push((event, facet));
                        glass_draw(transmitted)
                    },
                    |event, facet| {
                        connection.push((event, facet));
                        glass_draw(transmitted)
                    },
                )
                .unwrap();
            let baseline = SmithDielectricCameraScene::new(&camera, limit)
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings().surfaces,
                    |_| panic!("ordinary glass draw"),
                    |_, _| panic!("opaque draw"),
                    |_, _| glass_draw(transmitted),
                )
                .unwrap();
            assert_eq!(observed.radiance_sample().path(), baseline.path());
            assert_eq!(baseline.path().medium_index(), usize::from(transmitted));
            assert_eq!(continuation, vec![(0, 0), (0, 1)]);
            assert!(
                connection
                    .iter()
                    .enumerate()
                    .all(|(i, &(e, f))| e == 0 && f == i as u32)
            );
            assert!(!connection.is_empty());
            let selected = source.sample_direction(550.0, draw).unwrap().unwrap();
            let incoming = camera
                .launch(camera_ray(), 550.0, 0)
                .unwrap()
                .advance(1.0)
                .unwrap();
            let full = incoming
                .evaluate_smith_dielectric_connection(
                    cx,
                    Z.scale(-1.0),
                    glass,
                    surface,
                    selected.frame,
                    selected.stokes,
                    limit,
                    |_| glass_draw(transmitted),
                )
                .unwrap();
            let guide = incoming
                .evaluate_dielectric_connection(
                    Z.scale(-1.0),
                    glass,
                    surface,
                    selected.frame,
                    Stokes::unpolarized(1.0),
                )
                .unwrap()
                .pdf;
            let light_weight = power_weight(selected.pdf, guide).unwrap();
            let camera_weight = power_weight(guide, selected.pdf).unwrap();
            near(light_weight + camera_weight, 1.0, 2e-13);
            near(observed.environment_endpoint_weight(), camera_weight, 4e-6);
            near(
                observed.environment_radiance_w_m2_sr_nm(),
                full.factor_per_sr() * selected.radiance * light_weight / selected.pdf,
                4e-6,
            );
            near(
                observed.radiance_sample().radiance_w_m2_sr_nm(),
                5.0 + observed.environment_radiance_w_m2_sr_nm()
                    + (baseline.radiance_w_m2_sr_nm() - 5.0) * camera_weight,
                4e-6,
            );
        }
    });
}

#[test]
fn equal_index_bare_and_coated_rough_glass_keep_the_complete_path_only_measure() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let material = material(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&material),
            SceneMedium::Homogeneous(&material),
        ];
        let provenance =
            ConductorSource::try_new(ContentHash([149; 32]), ConductorDataStatus::Representative)
                .unwrap();
        let layers =
            [
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.0, 100e-9)
                    .unwrap(),
            ];
        let coating = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(&layers, 1, provenance).unwrap(),
        )
        .unwrap();
        for coated in [false, true] {
            let surface = DielectricSurface::try_rough(0.7).unwrap();
            let surface = if coated {
                surface.with_coating(coating)
            } else {
                surface
            };
            let native = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface,
                    },
                }],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let environment = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
            let limit = SmithWalkSettings::try_new(8).unwrap();
            for transmitted in [false, true] {
                let mut facets = 0;
                let observed = SmithDielectricEnvironmentCameraScene::new(&environment, limit)
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        550.0,
                        0,
                        settings(),
                        |_| panic!("ordinary glass draw"),
                        |_| panic!("mixed/null law must remain path-only"),
                        |_, _| panic!("opaque continuation"),
                        |_, _| panic!("opaque connection"),
                        |_, _| {
                            facets += 1;
                            glass_draw(transmitted)
                        },
                        |_, _| panic!("mixed/null fixed-direction connection"),
                    )
                    .unwrap();
                let baseline = SmithDielectricCameraScene::new(&camera, limit)
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        550.0,
                        0,
                        settings().surfaces,
                        |_| panic!("ordinary"),
                        |_, _| panic!("opaque"),
                        |_, _| glass_draw(transmitted),
                    )
                    .unwrap();
                assert_eq!(observed.radiance_sample(), &baseline);
                assert_eq!(observed.environment_endpoint_weight(), 1.0);
                assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(observed.environment_samples(), 0);
                assert_eq!(observed.visibility_queries(), 0);
                assert_eq!(facets, if coated { 2 } else { 0 });
            }
        }
    });
}

#[test]
fn unsupported_infinity_is_checked_on_the_sampled_glass_side_before_facet_work() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let surface = DielectricSurface::try_rough(0.7).unwrap();
        for (outside_absorption, inside_absorption, source_transmitted, exit_transmitted) in [
            (0.1, 0.1, true, true),
            (0.0, 0.1, true, false),
            (0.1, 0.0, false, true),
        ] {
            let air = material(1.0, outside_absorption);
            let glass = material(1.5, inside_absorption);
            let media = [
                SceneMedium::Homogeneous(&air),
                SceneMedium::Homogeneous(&glass),
            ];
            let native = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface,
                    },
                }],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let environment = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
            let limit = SmithWalkSettings::try_new(8).unwrap();
            let mut environment_draws = 0;
            let result = SmithDielectricEnvironmentCameraScene::new(&environment, limit)
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    |_| panic!("ordinary glass draw"),
                    |_| {
                        environment_draws += 1;
                        [0.5, if source_transmitted { 0.25 } else { 0.75 }]
                    },
                    |_, _| panic!("opaque continuation"),
                    |_, _| panic!("opaque connection"),
                    |_, _| glass_draw(exit_transmitted),
                    |_, _| panic!("unsupported selected medium consumed a connection facet"),
                )
                .unwrap();
            assert_eq!(result.environment_radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(result.visibility_queries(), 0);
            let neither = outside_absorption > 0.0 && inside_absorption > 0.0;
            assert_eq!(environment_draws, u32::from(!neither));
            assert_eq!(result.environment_samples(), environment_draws);
            if neither {
                assert_eq!(result.environment_endpoint_weight(), 1.0);
            } else {
                // The other side's NEE is still represented in expectation.
                // A failed side in this draw must not reset its camera partner.
                assert!(result.environment_endpoint_weight() > 0.0);
                assert!(result.environment_endpoint_weight() < 1.0);
            }
        }
    });
}
