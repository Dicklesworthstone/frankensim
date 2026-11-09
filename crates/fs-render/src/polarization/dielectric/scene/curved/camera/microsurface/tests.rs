use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::charts::TriMesh;
use crate::conductor::microsurface::{SmithEvent, SmithGgx};
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorOptics,
    ConductorSource, ConductorSurface,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::raster::CurvedCameraRaster;
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceEndpoint, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, SceneError, SurfaceBinding,
};
use crate::polarization::dielectric::transport::microsurface::SmithBeamError;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn metal() -> ConductorOptics {
    let source =
        ConductorSource::try_new(ContentHash([97; 32]), ConductorDataStatus::Representative)
            .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap()
}
fn plane(id: u64, z: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("org.frankensim.test.smith-camera", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-100.0, -100.0, z], [100.0, -100.0, z], [0.0, 100.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn bindings(surface: ConductorSurface) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface,
            },
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Opaque,
        },
    ]
}
fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 12.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 1.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    }
}
fn sources() -> [SurfaceEmitter; 2] {
    [
        SurfaceEmitter::new(
            1,
            RadianceSpectrum::try_constant(1.25).unwrap(),
            EmissionSides::Both,
        ),
        SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Both,
        ),
    ]
}
const DRAWS: [[f64; 3]; 3] = [[0.5, 0.8, 0.125], [0.5, 0.1, 0.375], [0.5, 0.2, 0.3]];
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn g1_two_facet_camera_retains_mueller_beer_and_emission_before_scattering() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.0), plane(9, 2.0)]).unwrap();
        let sources = sources();
        let bare = metal();
        let coated = bare.with_coating(
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        for index in [1.0, 1.33] {
            let material = glass(index, 0.2);
            let media = [SceneMedium::Homogeneous(&material)];
            for optics in [[bare], [coated]] {
                let native = CurvedPolarizedScene::try_new_with_conductors(
                    cx,
                    &geometry,
                    &media,
                    &optics,
                    bindings(ConductorSurface::try_rough(1.0).unwrap()),
                    100_000,
                )
                .unwrap();
                let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
                let smith = SmithCameraScene::new(&camera, SmithWalkSettings::try_new(2).unwrap());
                for wavelength in [450.0, 650.0] {
                    let mut draws = Vec::new();
                    let actual = smith
                        .sample_ray_analyzed_with_samples(
                            cx,
                            ray(),
                            wavelength,
                            0,
                            analyzer,
                            settings(),
                            |_| panic!("Smith metal borrowed a single-facet variate"),
                            |event, facet| {
                                draws.push((event, facet));
                                DRAWS[facet as usize]
                            },
                        )
                        .unwrap();
                    assert_eq!(draws, [(0, 0), (0, 1), (0, 2)]);
                    assert_eq!(
                        actual.endpoint(),
                        RadianceEndpoint::Surface { object_id: 9 }
                    );
                    assert_eq!(actual.path().surface_events(), 2);
                    // Independent traversal of the two actual facet normals;
                    // the camera driver and Smith beam adapter are not the oracle.
                    let mut expected = analyzer.launch(&camera, ray(), wavelength, 0).unwrap();
                    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
                    let mut depth = 0.0;
                    for draw in DRAWS {
                        let step = model.sample(expected.ray().dir, depth, draw).unwrap();
                        match step.event {
                            SmithEvent::Reflection {
                                depth: next, facet, ..
                            } => {
                                expected =
                                    expected.reflect_conductor(facet, &optics[0]).unwrap().beam;
                                depth = next;
                            }
                            SmithEvent::Escape => break,
                        }
                    }
                    let distance = 1.0 + 2.0 / expected.ray().dir.z;
                    let expected_radiance = 1.25 * analyzer.row().i * (-0.2_f64).exp()
                        + 4.0 * expected.power().i * (-0.2 * distance).exp();
                    near(actual.radiance_w_m2_sr_nm(), expected_radiance, 3e-6);
                    near(actual.path().traveled_m(), distance, 3e-6);
                    near(
                        actual.path().optical_path_length_m(),
                        index * distance,
                        4e-6,
                    );
                    assert_eq!(actual.launch_phase_index(), index);
                }
            }
        }
    });
}

#[test]
fn g4_facet_budget_and_cancellation_refuse_whole_camera_samples() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.0), plane(9, 2.0)]).unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &optics,
            bindings(ConductorSurface::try_rough(1.0).unwrap()),
            100_000,
        )
        .unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let bounded = SmithCameraScene::new(&camera, SmithWalkSettings::try_new(1).unwrap());
        assert_eq!(
            bounded.sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("ordinary draw"),
                |_, facet| DRAWS[facet as usize],
            ),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Smith(SmithBeamError::FacetLimit {
                    reflected_facets: 1
                },)
            )))
        );
        let smith = SmithCameraScene::new(&camera, SmithWalkSettings::try_new(8).unwrap());
        let mut calls = 0;
        let cancelled = smith.sample_ray_with_samples(
            cx,
            ray(),
            550.0,
            0,
            settings(),
            |_| panic!("ordinary draw"),
            |_, facet| {
                calls += 1;
                gate.request();
                DRAWS[facet as usize]
            },
        );
        assert_eq!(calls, 1);
        assert_eq!(
            cancelled,
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        );
    });
}

#[test]
fn g3_g5_smooth_limit_matches_native_camera_and_spectral_raster_exactly() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.0), plane(9, 2.0)]).unwrap();
        let material = glass(1.33, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &optics,
            bindings(ConductorSurface::SMOOTH),
            100_000,
        )
        .unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let smith = SmithCameraScene::new(&camera, SmithWalkSettings::try_new(1).unwrap());
        let expected = camera
            .sample_ray(cx, ray(), 550.0, 0, settings(), 61, 9)
            .unwrap();
        let actual = smith
            .sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("smooth ordinary draw"),
                |_, _| panic!("smooth microscopic draw"),
            )
            .unwrap();
        assert_eq!(actual, expected);
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
        let expected = raster
            .sample_pixel(cx, &camera, 0, settings(), 37, 4, 9)
            .unwrap();
        let actual = raster
            .sample_pixel_smith(cx, &smith, 0, settings(), 37, 4, 9)
            .unwrap();
        assert_eq!(actual, expected);
        let analyzer = CameraAnalyzer::try_circular(Vec3::new(1.0, 0.0, 0.0), true).unwrap();
        assert_eq!(
            raster
                .sample_pixel_smith_analyzed(cx, &smith, 0, analyzer, settings(), 37, 4, 9)
                .unwrap(),
            raster
                .sample_pixel_analyzed(cx, &camera, 0, analyzer, settings(), 37, 4, 9)
                .unwrap(),
        );
        assert_eq!(
            actual,
            raster
                .sample_pixel_smith(cx, &smith, 0, settings(), 37, 4, 9)
                .unwrap()
        );
        for lane in actual.spectral_samples() {
            assert_eq!(lane.path().termination(), PathTermination::Opaque);
        }
    });
}

#[test]
fn g1_brushed_camera_preserves_the_complete_anisotropic_radiance_path() {
    use crate::surface_detail::{BrushingPattern, SurfaceDetail};
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.0), plane(9, 2.0)]).unwrap();
        let material = glass(1.33, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let sources = sources();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        let ray = Ray {
            origin: Point3::new(1.0, 0.0, 1.0),
            ..ray()
        };
        let limits = SmithWalkSettings::try_new(64).unwrap();
        let draws = |event| match event {
            0 => [0.5, 0.95, 0.125],
            1 => [0.5, 0.1, 0.375],
            _ => [0.01, 0.05, 0.25],
        };
        for tilt in [0.0, 0.2] {
            let detail =
                SurfaceDetail::try_new(0.25, 4.0, BrushingPattern::Radial, tilt, None).unwrap();
            let mut bound = bindings(ConductorSurface::try_rough(0.25).unwrap());
            bound[0].interaction = Interaction::BrushedConductor {
                optics: 0,
                medium: 0,
                detail,
            };
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx, &geometry, &media, &optics, bound, 100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
            let smith = SmithCameraScene::new(&camera, limits);
            let incoming = analyzer.launch(&camera, ray, 550.0, 0).unwrap();
            let expected = incoming
                .sample_instanced_smith_brushed_conductor_for_transport(
                    cx,
                    &optics[0],
                    &detail,
                    RigidTransform::identity(),
                    Point3::new(1.0, 0.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    limits,
                    draws,
                    true,
                )
                .unwrap()
                .unwrap();
            assert!(expected.facet_events() > 1);
            let distance = 1.0 + 2.0 / expected.beam().ray().dir.z;
            let mut called = Vec::new();
            let actual = smith
                .sample_ray_analyzed_with_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| panic!("brushed Smith path used a single-facet draw"),
                    |surface, facet| {
                        called.push((surface, facet));
                        draws(facet)
                    },
                )
                .unwrap();
            assert_eq!(called.len(), expected.facet_events() as usize + 1);
            assert!(called.iter().all(|(surface, _)| *surface == 0));
            assert_eq!(
                actual.endpoint(),
                RadianceEndpoint::Surface { object_id: 9 }
            );
            assert_eq!(actual.path().surface_events(), 2);
            near(actual.path().traveled_m(), distance, 3e-6);
            near(actual.path().optical_path_length_m(), 1.33 * distance, 4e-6);
            let expected_radiance = 1.25 * analyzer.row().i * (-0.2_f64).exp()
                + 4.0 * expected.beam().power().i * (-0.2 * distance).exp();
            near(actual.radiance_w_m2_sr_nm(), expected_radiance, 3e-6);
        }
    });
}

mod dielectric;
