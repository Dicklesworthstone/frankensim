//! Original scene/source composition regressions for complete Smith lighting.
use super::*;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithCameraScene;
use crate::polarization::dielectric::transport::microsurface::{SmithBeamError, SmithWalkSettings};

fn metal() -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.3, 3.0).unwrap()
        }),
        ConductorSource::try_new(ContentHash([147; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}
fn conductor() -> Interaction {
    Interaction::RoughConductor {
        medium: 0,
        optics: 0,
        surface: ConductorSurface::try_rough(1.0).unwrap(),
    }
}

#[test]
fn smith_environment_combines_complementary_native_paths_and_keeps_emission_unweighted() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let metals = [metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: conductor(),
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
        let limit = SmithWalkSettings::try_new(2).unwrap();
        let owner = SmithEnvironmentCameraScene::new(&environment, limit);
        let mut continuation_draws = 0;
        let mut evaluation_draws = 0;
        let observed = owner
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                |_| panic!("Smith conductor must not consume an ordinary surface draw"),
                |event| {
                    assert_eq!(event, 0);
                    [0.5, 0.75]
                },
                |event, facet| {
                    assert_eq!(event, 0);
                    assert_eq!(facet, continuation_draws);
                    continuation_draws += 1;
                    [0.5, 0.0, 0.0]
                },
                |event, facet| {
                    assert_eq!(event, 0);
                    assert_eq!(facet, evaluation_draws);
                    evaluation_draws += 1;
                    [0.5, 0.0, 0.0]
                },
            )
            .unwrap();
        let baseline = SmithCameraScene::new(&camera, limit)
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings().surfaces,
                |_| panic!("ordinary stream"),
                |_, _| [0.5, 0.0, 0.0],
            )
            .unwrap();
        let reflection = ((0.3_f64 - 1.0).powi(2) + 9.0) / ((0.3_f64 + 1.0).powi(2) + 9.0);
        let le = source.evaluate(Z.scale(-1.0), 550.0).unwrap();
        // At this normal direction alpha=1 gives q1=p_environment=1/(4*pi).
        near(observed.environment_endpoint_weight(), 0.5, 2e-13);
        near(
            observed.environment_radiance_w_m2_sr_nm(),
            0.5 * reflection * le,
            4e-6,
        );
        near(
            observed.radiance_sample().radiance_w_m2_sr_nm(),
            5.0 + reflection * le,
            4e-6,
        );
        assert_eq!(observed.radiance_sample().path(), baseline.path());
        assert_eq!((continuation_draws, evaluation_draws), (2, 2));
        assert_eq!(
            (
                observed.environment_samples(),
                observed.visibility_queries()
            ),
            (1, 1)
        );
        let analyzed = owner
            .sample_ray_analyzed_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                CameraAnalyzer::try_linear(X, 0.0).unwrap(),
                settings(),
                |_| [0.0; 3],
                |_| [0.5, 0.75],
                |_, _| [0.5, 0.0, 0.0],
                |_, _| [0.5, 0.0, 0.0],
            )
            .unwrap();
        near(
            analyzed.radiance_sample().radiance_w_m2_sr_nm(),
            0.5 * observed.radiance_sample().radiance_w_m2_sr_nm(),
            4e-6,
        );
    });
}

#[test]
fn smith_environment_native_blockers_and_delta_chains_do_not_lose_source_paths() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let air = material(1.0, 0.0);
        let glass = material(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&glass),
        ];
        let metals = [metal()];
        let limit = SmithWalkSettings::try_new(8).unwrap();
        let continuing = |_: u32, facet: u32| {
            if facet == 0 {
                [0.5, 0.2, 0.0]
            } else {
                [0.000_001, 0.0, 0.0]
            }
        };
        let phi = (-0.6_f64).atan2(0.8).rem_euclid(core::f64::consts::TAU);
        for mirror in [false, true] {
            let geometry =
                InstanceScene::try_new(vec![patch(1, 0.8, 0.4, 0.15), patch(9, 0.0, 1.0, 0.3)])
                    .unwrap();
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: if mirror {
                            Interaction::Dielectric {
                                outside: 1,
                                inside: 0,
                                surface: DielectricSurface::SMOOTH,
                            }
                        } else {
                            Interaction::Opaque
                        },
                    },
                    SurfaceBinding {
                        object_id: 9,
                        interaction: conductor(),
                    },
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let environment = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
            let owner = SmithEnvironmentCameraScene::new(&environment, limit);
            let actual = owner
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    |_| [0.0; 3],
                    |_| [0.5, phi / core::f64::consts::TAU],
                    continuing,
                    |_, _| [0.5, 0.0, 0.0],
                )
                .unwrap();
            let baseline = SmithCameraScene::new(&camera, limit)
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings().surfaces,
                    |_| [0.0; 3],
                    continuing,
                )
                .unwrap();
            assert_eq!(actual.environment_radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(actual.environment_endpoint_weight(), 1.0);
            assert_eq!(actual.radiance_sample(), &baseline);
            assert_eq!(actual.radiance_sample().path().surface_events(), 2);
        }
    });
}

#[test]
fn smith_environment_refuses_the_complete_sample_when_light_evaluation_exceeds_its_limit() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let metals = [metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: conductor(),
            }],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
        let environment = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
        let owner =
            SmithEnvironmentCameraScene::new(&environment, SmithWalkSettings::try_new(1).unwrap());
        let error = owner
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                |_| panic!("ordinary draw"),
                |_| [0.5, 0.75],
                |_, _| panic!("no continuation after failed connection"),
                |_, facet| {
                    if facet == 0 {
                        [0.5, 0.8, 0.125]
                    } else {
                        [0.5, 0.1, 0.375]
                    }
                },
            )
            .unwrap_err();
        assert!(matches!(
            error,
            EnvironmentCameraError::Receiver(AreaLightCameraError::Camera(
                CameraRadianceError::Scene(CurvedSceneError::Scene(
                    crate::polarization::dielectric::scene::SceneError::Smith(
                        SmithBeamError::FacetLimit {
                            reflected_facets: 1
                        }
                    )
                ))
            ))
        ));
    });
}
#[test]
fn smith_environment_unsupported_absorbing_leg_keeps_continuation_without_light_draws() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let absorbing = material(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&absorbing)];
        let metals = [metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: conductor(),
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
        let limits = SmithWalkSettings::try_new(1).unwrap();
        let owner = SmithEnvironmentCameraScene::new(&environment, limits);
        let mut facets = 0;
        let observed = owner
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                |_| panic!("ordinary draw"),
                |_| panic!("unsupported source leg must not sample an environment direction"),
                |_, _| {
                    facets += 1;
                    [0.5, 0.0, 0.0]
                },
                |_, _| panic!("unsupported source leg must not sample connection facets"),
            )
            .unwrap();
        let baseline = SmithCameraScene::new(&camera, limits)
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings().surfaces,
                |_| panic!("ordinary draw"),
                |_, _| [0.5, 0.0, 0.0],
            )
            .unwrap();
        assert_eq!(facets, 2);
        assert_eq!(observed.radiance_sample(), &baseline);
        assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(observed.environment_endpoint_weight(), 1.0);
        assert_eq!(
            (
                observed.environment_samples(),
                observed.visibility_queries()
            ),
            (0, 0)
        );
        assert!(observed.radiance_sample().radiance_w_m2_sr_nm() > 0.0);
    });
}

mod dielectric;
