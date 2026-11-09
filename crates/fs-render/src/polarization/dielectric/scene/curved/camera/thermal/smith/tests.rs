use super::*;
use crate::charts::TriMesh;
use crate::conductor::{ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::SurfaceBinding;
use crate::polarization::dielectric::scene::curved::camera::microsurface::SmithDielectricCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::surface_detail::{BrushingPattern, SurfaceDetail};
use crate::thermal_surface::CoatingThermalEmission;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
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
        );
        f(&cx, &gate);
    });
}

fn metal(eta: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, 3.0).unwrap()
        }),
        ConductorSource::try_new(ContentHash([103; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}
fn thermal(optics: ConductorOptics) -> SmithThermalEmitter {
    let source =
        CoatingThermalEmission::try_opaque(&[], &[], optics, 2400.0, ContentHash([104; 32]))
            .unwrap();
    SmithThermalEmitter::try_new(
        1,
        BoundarySide::Positive,
        SmithCoatingThermalEmission::try_new(source, optics).unwrap(),
    )
    .unwrap()
}
fn plane(id: u64, z: f64, transform: RigidTransform) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("smith-thermal-camera", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-100.0, -100.0, z], [100.0, -100.0, z], [0.0, 100.0, z]],
            vec![[0, 1, 2]],
        )),
        transform,
    )
    .unwrap()
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.3, 0.2, 1.0),
        dir: Vec3::new(0.0, 0.0, -1.0),
    }
}
fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 4,
        max_distance_m: 64.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}
fn draw(path: u32, facet: u32) -> [f64; 3] {
    let b = philox4x32_10([path, facet, 0x7468_6572, 0], [37, 97]);
    core::array::from_fn(|i| (f64::from(b[i]) + 0.5) / 4_294_967_296.0)
}

#[test]
fn g1_native_thermal_faces_preserve_independent_complete_continuation_and_physical_prefix() {
    with_cx(|cx, _| {
        let host = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&host)];
        let optics = [metal(0.3)];
        let emitter = [thermal(optics[0])];
        let transform =
            RigidTransform::try_new([0.0, 0.0, 0.3_f64.sin(), 0.3_f64.cos()], [0.0; 3]).unwrap();
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, transform),
            plane(9, 3.0, RigidTransform::identity()),
        ])
        .unwrap();
        let detail =
            SurfaceDetail::try_new(0.4, 2.0, BrushingPattern::Circular, 0.0, None).unwrap();
        let rough = ConductorSurface::try_rough(0.8).unwrap();
        let microsurface = SmithWalkSettings::try_new(512).unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
        for interaction in [
            Interaction::SmoothConductor {
                optics: 0,
                medium: 0,
            },
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: rough,
            },
            Interaction::BrushedConductor {
                optics: 0,
                medium: 0,
                detail,
            },
        ] {
            let bindings = vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction,
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
            ];
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx, &geometry, &media, &optics, bindings, 100_000,
            )
            .unwrap();
            let sources = [SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(4.0).unwrap(),
                EmissionSides::Both,
            )];
            let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
            let ordinary = SmithDielectricCameraScene::new(&camera, microsurface);
            let thermal =
                SmithThermalCameraScene::try_new(cx, &camera, &emitter, microsurface).unwrap();
            let empty = SmithThermalCameraScene::try_new(cx, &camera, &[], microsurface).unwrap();
            let (path_id, plain) = (0..1024)
                .find_map(|id| {
                    ordinary
                        .sample_ray_analyzed_with_samples(
                            cx,
                            ray(),
                            550.0,
                            0,
                            analyzer,
                            settings(),
                            |_| panic!("ordinary draw"),
                            |_, facet| draw(id, facet),
                            |_, _| panic!("glass draw"),
                        )
                        .ok()
                        .map(|sample| (id, sample))
                })
                .unwrap();
            let prefix = analyzer
                .launch(&camera, ray(), 550.0, 0)
                .unwrap()
                .advance(1.0)
                .unwrap();
            let unit = normalized_camera(prefix, Stokes::unpolarized(1.0)).unwrap();
            for thermal_id in [31, 47] {
                let expected = match interaction {
                    Interaction::SmoothConductor { .. } => emitter[0].source.sample(
                        cx,
                        unit,
                        Vec3::new(0.0, 0.0, 1.0),
                        ConductorSurface::SMOOTH,
                        microsurface,
                        |_| panic!("smooth source draws nothing"),
                    ),
                    Interaction::RoughConductor { .. } => emitter[0].source.sample(
                        cx,
                        unit,
                        Vec3::new(0.0, 0.0, 1.0),
                        rough,
                        microsurface,
                        |facet| draw(thermal_id, facet),
                    ),
                    Interaction::BrushedConductor { .. } => {
                        let local = transform
                            .inverse()
                            .unwrap()
                            .transform_point(Point3::new(0.3, 0.2, 0.0));
                        let frame = unit
                            .smith_instanced_brushed_frame(
                                &detail,
                                transform,
                                local,
                                Vec3::new(0.0, 0.0, 1.0),
                            )
                            .unwrap()
                            .unwrap();
                        let (ax, ay) = detail.alpha_pair_for_mark(frame.mark_covered);
                        emitter[0].source.sample_anisotropic(
                            cx,
                            unit,
                            frame.shading_normal,
                            frame.tangent,
                            ax,
                            ay,
                            microsurface,
                            |facet| draw(thermal_id, facet),
                        )
                    }
                    _ => unreachable!(),
                }
                .unwrap();
                let mut addresses = Vec::new();
                let actual = thermal
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("ordinary draw"),
                        |surface, facet| {
                            assert_eq!(surface, 0);
                            draw(path_id, facet)
                        },
                        |_, _| panic!("glass draw"),
                        |surface, facet| {
                            addresses.push((surface, facet));
                            draw(thermal_id, facet)
                        },
                    )
                    .unwrap();
                assert_eq!(actual.path(), plain.path());
                assert!(addresses.iter().all(|(surface, _)| *surface == 0));
                if matches!(interaction, Interaction::SmoothConductor { .. }) {
                    assert!(addresses.is_empty());
                }
                let expected = expected.radiance_w_m2_sr_nm() * prefix.power().i;
                let observed = actual.radiance_w_m2_sr_nm() - plain.radiance_w_m2_sr_nm();
                assert!(
                    (observed - expected).abs() <= 4e-7 * expected.max(1.0),
                    "{observed} != {expected}"
                );
            }
            assert_eq!(
                empty
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("ordinary draw"),
                        |_, facet| draw(path_id, facet),
                        |_, _| panic!("glass draw"),
                        |_, _| panic!("empty thermal table must not draw"),
                    )
                    .unwrap(),
                plain
            );
        }
    });
}

#[test]
fn g1_thermal_binding_refuses_wrong_optics_tilt_overlap_and_nonvacuum_without_partial_output() {
    with_cx(|cx, gate| {
        let optics = [metal(0.3)];
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, RigidTransform::identity())]).unwrap();
        let rough = ConductorSurface::try_rough(0.8).unwrap();
        let limits = SmithWalkSettings::try_new(128).unwrap();
        let good = thermal(optics[0]);
        let wrong = thermal(metal(0.4));
        for (index, tilt) in [(1.0, 0.0), (1.0, 0.01), (1.3, 0.0)] {
            let host = DielectricGlass::new(
                CauchyIor::try_constant(index).unwrap(),
                BeerLambertAbsorption::CLEAR,
                GlassProvenance::Custom,
            );
            let media = [SceneMedium::Homogeneous(&host)];
            let interaction = if tilt == 0.0 {
                Interaction::RoughConductor {
                    optics: 0,
                    medium: 0,
                    surface: rough,
                }
            } else {
                Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail: SurfaceDetail::try_new(0.4, 2.0, BrushingPattern::Circular, tilt, None)
                        .unwrap(),
                }
            };
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &optics,
                vec![SurfaceBinding {
                    object_id: 1,
                    interaction,
                }],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
            assert!(
                SmithThermalCameraScene::try_new(cx, &camera, &[wrong.clone()], limits).is_err()
            );
            assert!(
                SmithThermalCameraScene::try_new(
                    cx,
                    &camera,
                    &[good.clone(), good.clone()],
                    limits
                )
                .is_err()
            );
            let source = [good.clone()];
            let prepared = SmithThermalCameraScene::try_new(cx, &camera, &source, limits);
            if tilt > 0.0 {
                assert!(prepared.is_err());
                continue;
            }
            let prepared = prepared.unwrap();
            if index != 1.0 {
                assert!(matches!(
                    prepared.sample_ray_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        settings(),
                        |_| panic!("invalid host"),
                        |_, _| panic!("thermal refusal precedes continuation"),
                        |_, _| panic!("invalid host"),
                        |_, _| panic!("invalid host precedes thermal draws"),
                    ),
                    Err(SmithThermalCameraError::Thermal {
                        object_id: 1,
                        error: SmithThermalError::VacuumRequired
                    })
                ));
            }
            let overlap = [SurfaceEmitter::new(
                1,
                RadianceSpectrum::try_constant(1.0).unwrap(),
                EmissionSides::Both,
            )];
            let occupied = CurvedCameraScene::try_new(cx, &native, &overlap, None).unwrap();
            assert!(SmithThermalCameraScene::try_new(cx, &occupied, &source, limits).is_err());
            if index != 1.0 {
                gate.request();
                assert!(matches!(
                    prepared.sample_ray_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        settings(),
                        |_| panic!("cancelled"),
                        |_, _| panic!("cancelled"),
                        |_, _| panic!("cancelled"),
                        |_, _| panic!("cancelled"),
                    ),
                    Err(SmithThermalCameraError::Camera(CameraRadianceError::Scene(
                        CurvedSceneError::Scene(SceneError::Cancelled)
                    )))
                ));
            }
        }
    });
}
