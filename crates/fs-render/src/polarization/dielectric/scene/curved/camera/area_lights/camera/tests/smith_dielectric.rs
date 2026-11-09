//! Native finite-source integration for the complete two-sided glass law.
use super::super::microsurface::SmithDielectricAreaLightCameraScene;
use super::*;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithDielectricCameraScene;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;

fn full_glass_binding(surface: DielectricSurface) -> SurfaceBinding {
    // The panel's physical normal is +z; the camera below it occupies inside 0.
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::Dielectric {
            outside: 1,
            inside: 0,
            surface,
        },
    }
}

#[test]
fn g1_full_glass_area_matches_native_response_on_both_sides_and_owns_direct_hit() {
    with_cx(|cx| {
        let air = glass(1.0, 0.1);
        let material = glass(1.5, 0.2);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let surface =
            DielectricSurface::try_anisotropic_rough(0.35, 0.65, Vec3::new(1.0, 0.0, 0.0)).unwrap();
        let limits = SmithWalkSettings::try_new(1).unwrap();
        for transmitted in [false, true] {
            let source_z = if transmitted { 3.0 } else { -1.0 };
            let source_side = if transmitted {
                BoundarySide::Negative
            } else {
                BoundarySide::Positive
            };
            let source_medium = usize::from(transmitted);
            let branch = if transmitted { 0.999 } else { 0.001 };
            let draw = |_: u32| [0.5, 0.0, 0.0, branch];
            let geometry =
                InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, source_z, 0.3)])
                    .unwrap();
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![full_glass_binding(surface), opaque(9)],
                100_000,
            )
            .unwrap();
            let emitters = [
                SurfaceEmitter::new(
                    1,
                    RadianceSpectrum::try_constant(1.25).unwrap(),
                    EmissionSides::Both,
                ),
                SurfaceEmitter::new(
                    9,
                    RadianceSpectrum::try_constant(LE).unwrap(),
                    EmissionSides::Both,
                ),
            ];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0 - source_z))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side,
                initial_medium: source_medium,
                source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let full = SmithDielectricAreaLightCameraScene::try_new(cx, &owner, limits).unwrap();
            let baseline = SmithDielectricCameraScene::new(&camera, limits)
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings().surfaces,
                    |_| panic!("ordinary draw"),
                    |_, _| panic!("opaque draw"),
                    |_, f| draw(f),
                )
                .unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
            let mut continuation = Vec::new();
            let mut connections = Vec::new();
            let actual = full
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("ordinary draw"),
                    |_, _| panic!("opaque draw"),
                    |event, facet| {
                        continuation.push((event, facet));
                        draw(facet)
                    },
                    |event, binding| {
                        assert_eq!((event, binding), (0, 0));
                        [4.0 / 9.0, 0.5]
                    },
                    |_, _, _, _| panic!("opaque connection"),
                    |event, binding, channel, facet| {
                        connections.push((event, binding, channel, facet));
                        draw(facet)
                    },
                )
                .unwrap();
            assert_eq!(actual.radiance_sample().path(), baseline.path());
            assert_eq!(
                actual.radiance_sample().path().medium_index(),
                source_medium
            );
            assert_eq!(continuation, [(0, 0), (0, 1)]);
            assert!(
                connections
                    .iter()
                    .all(|&(e, b, c, _)| (e, b, c) == (0, 0, 0))
            );
            assert!(!connections.is_empty());
            let source_frame = PolarizationFrame::try_new(
                Vec3::new(0.0, 0.0, if transmitted { -1.0 } else { 1.0 }),
                Vec3::new(1.0, 0.0, 0.0),
            )
            .unwrap();
            // Isolate integration: independent native kernel evaluation is
            // multiplied by the analytic flat-source geometry and two Beer legs.
            let response = camera
                .launch(camera_ray(), 550.0, 0)
                .unwrap()
                .evaluate_smith_dielectric_connection(
                    cx,
                    Vec3::new(0.0, 0.0, 1.0),
                    material,
                    surface,
                    source_frame,
                    Stokes::unpolarized(1.0),
                    limits,
                    draw,
                )
                .unwrap();
            let absorption = if transmitted { 0.2 } else { 0.1 };
            let expected =
                response.factor_per_sr() * LE * AREA / 4.0 * (-0.1_f64 - 2.0 * absorption).exp();
            near(actual.area_radiance_w_m2_sr_nm(), expected, 5e-6);
            near(
                actual.radiance_sample().radiance_w_m2_sr_nm(),
                1.25 * (-0.1_f64).exp() + expected,
                5e-6,
            );
            assert_eq!((actual.area_samples(), actual.images()), (1, 1));
            // A primary source hit is never part of this receiver's partition.
            let primary = full
                .sample_ray_with_samples(
                    cx,
                    Ray {
                        origin: Point3::new(0.0, 0.0, source_z - 0.2),
                        dir: Vec3::new(0.0, 0.0, 1.0),
                    },
                    550.0,
                    source_medium,
                    settings(),
                    &mut workspace,
                    |_| panic!("ordinary"),
                    |_, _| panic!("opaque"),
                    |_, _| panic!("glass"),
                    |_, _| panic!("area"),
                    |_, _, _, _| panic!("opaque connection"),
                    |_, _, _, _| panic!("glass connection"),
                )
                .unwrap();
            near(
                primary.radiance_sample().radiance_w_m2_sr_nm(),
                LE * (-0.2 * absorption).exp(),
                5e-6,
            );
            assert_eq!(primary.area_samples(), 0);
        }
    });
}

#[test]
fn g0_equal_index_coated_and_bare_glass_retain_complete_continuation_measure() {
    use crate::conductor::{CauchyCoatingLayer, CoherentPeriodicCoating};
    use crate::dielectric::DielectricCoating;
    with_cx(|cx| {
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([57; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let cell =
            [
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.0, 100e-9)
                    .unwrap(),
            ];
        let bare = DielectricSurface::try_rough(0.8).unwrap();
        let coated = bare.with_coating(
            DielectricCoating::try_new(CoherentPeriodicCoating::try_new(&cell, 1, source).unwrap())
                .unwrap(),
        );
        let a = glass(1.5, 0.1);
        let b = glass(1.5, 0.2);
        let media = [SceneMedium::Homogeneous(&a), SceneMedium::Homogeneous(&b)];
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, 3.0, 0.3)]).unwrap();
        for surface in [bare, coated] {
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![full_glass_binding(surface), opaque(9)],
                100_000,
            )
            .unwrap();
            let emitters = [SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(LE).unwrap(),
                EmissionSides::Both,
            )];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, -1.0))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side: BoundarySide::Negative,
                initial_medium: 1,
                source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let limits = SmithWalkSettings::try_new(8).unwrap();
            let full = SmithDielectricAreaLightCameraScene::try_new(cx, &owner, limits).unwrap();
            let draw = |_: u32, _: u32| {
                assert!(surface.coating().is_some(), "bare null drew a facet");
                [0.5, 0.0, 0.0, 0.999]
            };
            let baseline = SmithDielectricCameraScene::new(&camera, limits)
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings().surfaces,
                    |_| panic!("ordinary"),
                    |_, _| panic!("opaque"),
                    draw,
                )
                .unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
            let actual = full
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("ordinary"),
                    |_, _| panic!("opaque"),
                    draw,
                    |_, _| panic!("mixed measure sampled an area source"),
                    |_, _, _, _| panic!("opaque connection"),
                    |_, _, _, _| panic!("mixed measure made a continuous connection"),
                )
                .unwrap();
            assert_eq!(actual.radiance_sample(), &baseline);
            assert_eq!(actual.area_samples(), 0);
            assert_eq!(actual.area_radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(actual.radiance_sample().path().medium_index(), 1);
        }
    });
}

#[test]
fn g0_glass_area_refuses_incomplete_connection_without_returning_partial_radiance() {
    with_cx(|cx| {
        let air = glass(1.0, 0.0);
        let material = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let surface = DielectricSurface::try_rough(0.6).unwrap();
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, 3.0, 0.3)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![full_glass_binding(surface), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(LE).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, -1.0))];
        let bindings = [AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Negative,
            initial_medium: 1,
            source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let full = SmithDielectricAreaLightCameraScene::try_new(
            cx,
            &owner,
            SmithWalkSettings::try_new(1).unwrap(),
        )
        .unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 1).unwrap();
        let result = full.sample_ray_with_samples(
            cx,
            camera_ray(),
            550.0,
            0,
            settings(),
            &mut workspace,
            |_| panic!("ordinary"),
            |_, _| panic!("opaque"),
            |_, _| panic!("continuation before refusal"),
            |_, _| [4.0 / 9.0, 0.5],
            |_, _, _, _| panic!("opaque connection"),
            |_, _, _, _| [0.0, 0.0, 0.0, 0.999],
        );
        assert!(matches!(
            result,
            Err(AreaLightCameraError::Camera(CameraRadianceError::Scene(
                CurvedSceneError::Scene(
                    crate::polarization::dielectric::scene::SceneError::SmithDielectric(_)
                )
            )))
        ));
    });
}
