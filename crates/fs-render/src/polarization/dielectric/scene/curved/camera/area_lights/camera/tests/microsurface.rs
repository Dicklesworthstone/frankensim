//! Direct area support and native endpoint partition for complete Smith metals.
use super::super::microsurface::SmithAreaLightCameraScene;
use super::*;
use crate::polarization::dielectric::scene::SceneError;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithCameraScene;
use crate::polarization::dielectric::transport::microsurface::{SmithBeamError, SmithWalkSettings};

fn flat_draw(_: u32, _: u32) -> [f64; 3] {
    [0.5, 0.0, 0.0]
}

#[test]
fn g1_smith_area_uses_full_response_and_owns_only_direct_emitter_hits() {
    with_cx(|cx| {
        let x = 0.2;
        let ray = Ray {
            origin: Point3::new(x, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let geometry =
            InstanceScene::try_new(vec![panel(1, x, 1.0, 1.0), panel(9, x, -1.0, 0.3)]).unwrap();
        let material = glass(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [metal()];
        let detail = SurfaceDetail::try_new(0.2, 3.0, BrushingPattern::Radial, 0.0, None).unwrap();
        let laws = [
            (
                Interaction::RoughConductor {
                    optics: 0,
                    medium: 0,
                    surface: ConductorSurface::try_rough(0.2).unwrap(),
                },
                (0.2, 0.2),
            ),
            (
                Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail,
                },
                detail.alpha_pair(),
            ),
        ];
        for (interaction, (ax, ay)) in laws {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &conductors,
                &[],
                &[],
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction,
                    },
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let emitters = [
                SurfaceEmitter::new(
                    1,
                    RadianceSpectrum::try_constant(1.25).unwrap(),
                    EmissionSides::Both,
                ),
                emitter(9),
            ];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            // Deliberately wrong authored direction: this opt-in owns its
            // unique point-to-contact direct proposal, not the old image seed.
            let seeds = [seed(Vec3::new(1.0, 0.0, 1.0))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side: BoundarySide::Positive,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let limits = SmithWalkSettings::try_new(1).unwrap();
            let smith = SmithAreaLightCameraScene::try_new(cx, &owner, limits).unwrap();
            let baseline = SmithCameraScene::new(&camera, limits)
                .sample_ray_with_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    settings().surfaces,
                    |_| panic!("single-facet draw used"),
                    flat_draw,
                )
                .unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
            let mut continuation = Vec::new();
            let mut connections = Vec::new();
            let actual = smith
                .sample_ray_with_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("single-facet draw used"),
                    |event, facet| {
                        continuation.push((event, facet));
                        flat_draw(event, facet)
                    },
                    |event, binding| {
                        assert_eq!((event, binding), (0, 0));
                        [4.0 / 9.0, 0.5]
                    },
                    |event, binding, channel, facet| {
                        connections.push((event, binding, channel, facet));
                        flat_draw(event, facet)
                    },
                )
                .unwrap();
            assert_eq!(continuation, [(0, 0), (0, 1)]);
            assert_eq!(connections, [(0, 0, 0, 0), (0, 0, 0, 1)]);
            assert_eq!(actual.radiance_sample().path(), baseline.path());
            let fresnel = ((1.0_f64 - 0.2).powi(2) + 9.0) / ((1.0_f64 + 0.2).powi(2) + 9.0);
            let direct =
                fresnel * LE * AREA / (16.0 * core::f64::consts::PI * ax * ay) * (-0.3_f64).exp();
            let primary = 1.25 * (-0.1_f64).exp();
            near(actual.area_radiance_w_m2_sr_nm(), direct, 4e-6);
            near(
                actual.radiance_sample().radiance_w_m2_sr_nm(),
                primary + direct,
                4e-6,
            );
            near(
                baseline.radiance_w_m2_sr_nm(),
                primary + fresnel * LE * (-0.3_f64).exp(),
                4e-6,
            );
            assert_eq!((actual.area_samples(), actual.images()), (1, 1));
            assert!(
                actual.trace_evaluations() >= 1
                    && actual.field_evaluations() > 0
                    && actual.query_work() > 0
            );
            let repeat = smith
                .sample_ray_with_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("ordinary draw"),
                    flat_draw,
                    |_, _| [4.0 / 9.0, 0.5],
                    |event, _, _, facet| flat_draw(event, facet),
                )
                .unwrap();
            assert_eq!(actual, repeat);
            // A primary source hit has no competing receiver proposal.
            let primary = smith
                .sample_ray_with_samples(
                    cx,
                    Ray {
                        dir: Vec3::new(0.0, 0.0, -1.0),
                        ..ray
                    },
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("ordinary draw"),
                    |_, _| panic!("facet"),
                    |_, _| panic!("area"),
                    |_, _, _, _| panic!("connection"),
                )
                .unwrap();
            near(
                primary.radiance_sample().radiance_w_m2_sr_nm(),
                LE * (-0.1_f64).exp(),
                4e-6,
            );
            assert_eq!(primary.area_samples(), 0);
        }
    });
}

#[test]
fn g1_smith_area_unselected_original_faces_and_sides_keep_full_emission() {
    with_cx(|cx| {
        let mut vertices = Vec::new();
        for x in [3.0, 0.0] {
            vertices.extend([[x - 0.3, -0.3, -1.0], [x + 0.3, -0.3, -1.0], [x, 0.6, -1.0]]);
        }
        let lamp = GeometryInstance::try_new(
            9,
            hash_domain("smith-two-source-faces", b"native"),
            SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2], [3, 4, 5]])),
            RigidTransform::identity(),
        )
        .unwrap();
        let geometry = InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), lamp]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &optics,
            &[],
            &[],
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::RoughConductor {
                        optics: 0,
                        medium: 0,
                        surface: ConductorSurface::try_rough(0.2).unwrap(),
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(LE).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let limits = SmithWalkSettings::try_new(1).unwrap();
        let baseline = SmithCameraScene::new(&camera, limits)
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings().surfaces,
                |_| [0.0; 3],
                flat_draw,
            )
            .unwrap();
        for (face, side) in [(0, BoundarySide::Positive), (1, BoundarySide::Negative)] {
            let source = CurvedAreaEmitter::from_camera_faces(
                cx,
                &camera,
                9,
                &[face],
                AreaEmitterLimits::default(),
            )
            .unwrap();
            let bindings = [AreaLightBinding {
                source: &source,
                receiver_object_id: 1,
                source_side: side,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let smith = SmithAreaLightCameraScene::try_new(cx, &owner, limits).unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
            let result = smith
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| [0.0; 3],
                    flat_draw,
                    |_, _| [4.0 / 9.0, 0.5],
                    |event, _, _, facet| flat_draw(event, facet),
                )
                .unwrap();
            assert_eq!(result.radiance_sample().path(), baseline.path());
            near(
                result.radiance_sample().radiance_w_m2_sr_nm() - result.area_radiance_w_m2_sr_nm(),
                baseline.radiance_w_m2_sr_nm(),
                5e-6,
            );
            if side == BoundarySide::Negative {
                assert_eq!(result.area_radiance_w_m2_sr_nm(), 0.0);
            }
        }
    });
}

#[test]
fn g1_direct_bundle_visibility_stops_at_mirrors_and_dielectrics_with_native_work() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 0.0, 1.0, 1.0),
            panel(5, 0.0, 0.0, 1.0),
            panel(9, 0.0, -1.0, 0.3),
        ])
        .unwrap();
        let air = glass(1.0, 0.0);
        let glass = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&glass),
        ];
        let optics = [metal()];
        for blocker in [
            Interaction::SmoothConductor {
                optics: 0,
                medium: 0,
            },
            Interaction::Dielectric {
                outside: 1,
                inside: 0,
                surface: DielectricSurface::SMOOTH,
            },
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &optics,
                &[],
                &[],
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::RoughConductor {
                            optics: 0,
                            medium: 0,
                            surface: ConductorSurface::try_rough(0.2).unwrap(),
                        },
                    },
                    SurfaceBinding {
                        object_id: 5,
                        interaction: blocker,
                    },
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let emitters = [emitter(9)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let sample = light.sample_area([4.0 / 9.0, 0.5]).unwrap();
            let seeds = [BundleImageSeed {
                initial_direction: Vec3::new(0.0, 0.0, 1.0),
                policy: BundleBranchPolicy::Sequence(&[]),
            }];
            let (launch, source, _) = light
                .launch(
                    sample,
                    BoundarySide::Positive,
                    0,
                    Vec3::new(0.0, 1.0, 0.0),
                    550.0,
                    &seeds,
                )
                .unwrap();
            let target = BundleConnectionTarget::try_direct_surface_receiver(
                1,
                Point3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -1.0),
                2e-7,
            )
            .unwrap();
            let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
            let catalog = scene
                .discover_bundle_images_from_surface_with_workspace(
                    cx,
                    launch,
                    source,
                    target,
                    settings().connection,
                    &seeds,
                    image_options(1),
                    &mut workspace,
                )
                .unwrap();
            catalog.validate_proposal_support().unwrap();
            assert!(catalog.images().is_empty());
            assert_eq!(catalog.trace_evaluations(), 1);
            assert!(catalog.field_evaluations() > 0 && catalog.query_work() > 0);
            assert!(matches!(catalog.seed_results()[0].outcome(),
                crate::polarization::dielectric::scene::curved::bundle::connection::images::BundleImageSeedOutcome::Miss(
                    BundleConnectionError::ReceiverNotReached(PathTermination::Opaque))));
        }
    });
}

#[test]
fn g1_smith_area_preserves_unrepresented_transmission_chains() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 0.0, 1.0, 1.0),
            panel(2, 0.0, -0.2, 1.0),
            panel(3, 0.0, -0.4, 1.0),
            panel(9, 0.0, -1.0, 0.3),
        ])
        .unwrap();
        let air = glass(1.0, 0.1);
        let slab = glass(1.5, 0.1);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&slab),
        ];
        let optics = [metal()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &optics,
            &[],
            &[],
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::RoughConductor {
                        optics: 0,
                        medium: 0,
                        surface: ConductorSurface::try_rough(0.2).unwrap(),
                    },
                },
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 3,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let limits = SmithWalkSettings::try_new(1).unwrap();
        let smith = SmithAreaLightCameraScene::try_new(cx, &owner, limits).unwrap();
        let baseline = SmithCameraScene::new(&camera, limits)
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings().surfaces,
                |_| [0.5; 3],
                flat_draw,
            )
            .unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
        let result = smith
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                &mut workspace,
                |_| [0.5; 3],
                flat_draw,
                |_, _| [4.0 / 9.0, 0.5],
                |_, _, _, _| panic!("blocked source must not evaluate a Smith response"),
            )
            .unwrap();
        assert!(baseline.radiance_w_m2_sr_nm() > 0.0);
        assert_eq!(result.radiance_sample(), &baseline);
        assert_eq!(result.area_radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(
            (
                result.area_samples(),
                result.images(),
                result.trace_evaluations()
            ),
            (1, 0, 1)
        );
        assert!(result.field_evaluations() > 0 && result.query_work() > 0);
    });
}

#[test]
fn g4_smith_area_rejects_incompatible_receivers_and_refuses_facet_exhaustion() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        for interaction in [
            diffuse(1, 0).interaction,
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::try_rough(1.0).unwrap(),
            },
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &optics,
                &[],
                &reflectances,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction,
                    },
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let emitters = [emitter(9)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side: BoundarySide::Positive,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let smith = SmithAreaLightCameraScene::try_new(
                cx,
                &owner,
                SmithWalkSettings::try_new(1).unwrap(),
            );
            if matches!(interaction, Interaction::DepolarizingDiffuse { .. }) {
                assert!(matches!(
                    smith,
                    Err(AreaLightCameraError::InvalidBinding { index: 0 })
                ));
                continue;
            }
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
            let failed = smith.unwrap().sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                flat_draw,
                |_, _| [4.0 / 9.0, 0.5],
                |_, _, _, facet| {
                    if facet == 0 {
                        [0.5, 0.8, 0.125]
                    } else {
                        [0.5, 0.1, 0.375]
                    }
                },
            );
            assert!(matches!(
                failed,
                Err(AreaLightCameraError::Camera(CameraRadianceError::Scene(
                    CurvedSceneError::Scene(SceneError::Smith(SmithBeamError::FacetLimit {
                        reflected_facets: 1
                    }))
                )))
            ));
        }
    });
}
