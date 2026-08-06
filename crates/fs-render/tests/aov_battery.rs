//! G0/G3/G5 battery for aligned cinematic AOV accumulation and EXR artifacts.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, RunError, RunId, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::aov::{
    CinematicAovConfig, CinematicAovError, CinematicAovLimits, CinematicAovProfile,
    CinematicAovProvenance, validity,
};
use fs_render::camera::{AnimatedCamera, Aperture, CutSide, PhysicalCamera};
use fs_render::charts::TriMesh;
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::instances::{RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    AdaptiveDecision, AdaptiveSamplingConfig, Camera, CinematicAovExecutionError, DirectStrategy,
    Film, Material, Primitive, RectLight, RenderExecutionConfig, RenderExecutionError,
    RenderWorkerPool, Sampler, Scene, Settings, Shape, render_cinematic,
    render_cinematic_adaptive_with_aovs, render_cinematic_range_with_aovs,
    render_cinematic_with_aovs, render_cinematic_with_aovs_execution, trace_cinematic_pixel_sample,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate(&gate, f)
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xa0_51,
                kernel_id: 6,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn hash(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.cinematic-aov", label.as_bytes())
}

fn quad() -> TriMesh {
    TriMesh::new(
        vec![
            [-0.75, -0.75, 0.0],
            [0.75, -0.75, 0.0],
            [0.75, 0.75, 0.0],
            [-0.75, 0.75, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn backface_quad() -> TriMesh {
    TriMesh::new(
        vec![
            [-0.75, -0.75, 0.0],
            [0.75, -0.75, 0.0],
            [0.75, 0.75, 0.0],
            [-0.75, 0.75, 0.0],
        ],
        vec![[0, 2, 1], [0, 3, 2]],
    )
}

fn left_pixel_quad() -> TriMesh {
    // At z=0 with the 2x1 test camera, the left pixel footprint spans
    // x in [-1, 0) and y in [-0.5, 0.5]. Keep the right edge just below zero
    // so every sub-pixel sample in pixel 0 hits and every sample in pixel 1
    // misses, including a hypothetical zero jitter value.
    TriMesh::new(
        vec![
            [-1.1, -0.75, 0.0],
            [-1.0e-9, -0.75, 0.0],
            [-1.0e-9, 0.75, 0.0],
            [-1.1, 0.75, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn scene(animated: bool) -> Scene {
    let geometry = SharedGeometry::mesh(quad());
    let shape = if animated {
        let start = TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap();
        let end = TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.5, 0.0, 0.0]).unwrap(),
            [1.0, 0.0, 0.0],
        )
        .unwrap();
        Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                101,
                hash("geometry"),
                geometry,
                RigidTransformTrajectory::try_new(vec![start, end]).unwrap(),
            )
            .unwrap(),
        )
    } else {
        Shape::Mesh(quad())
    };
    let reflectance = lift_rgb([0.7, 0.4, 0.2]);
    let emission = (reflectance, 2.0);
    Scene {
        primitives: vec![
            Primitive {
                shape,
                material: Material::Lambertian { reflectance },
                emission: Some(emission),
            },
            Primitive {
                shape: Shape::Mesh(TriMesh::new(
                    vec![
                        [9.0, -0.5, 0.0],
                        [10.0, -0.5, 0.0],
                        [10.0, 0.5, 0.0],
                        [9.0, 0.5, 0.0],
                    ],
                    vec![[0, 1, 2], [0, 2, 3]],
                )),
                material: Material::Lambertian { reflectance },
                emission: Some(emission),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(9.0, -0.5, 0.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 1,
            emission,
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.25,
        },
    }
}

fn camera_at(scene: &Scene, eye: Point3) -> AnimatedCamera {
    let physical = PhysicalCamera::try_legacy_compatible(
        eye,
        scene.camera.forward,
        scene.camera.up,
        scene.camera.half_tan,
        2.0,
        Aperture::try_circular(0.0).unwrap(),
    )
    .unwrap();
    AnimatedCamera::try_static(701, 0.0, 1.0, physical).unwrap()
}

fn camera(scene: &Scene) -> AnimatedCamera {
    camera_at(scene, scene.camera.eye)
}

fn settings(spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0xa0_51_0001,
    }
}

fn execution(
    tile_width: u32,
    tile_height: u32,
    workers: usize,
    memory_limit_bytes: u64,
    run: u64,
) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(
        tile_width,
        tile_height,
        workers,
        memory_limit_bytes,
        RunId(run),
    )
    .expect("valid cinematic AOV execution policy")
}

struct CancellingChart {
    gate: Arc<CancelGate>,
    evaluations: Arc<AtomicUsize>,
    cancel_at: usize,
}

impl Chart for CancellingChart {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        let evaluation = self.evaluations.fetch_add(1, Ordering::SeqCst) + 1;
        if evaluation == self.cancel_at {
            self.gate.request();
        }
        SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.5,
        }
        .eval(point, cx)
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "cinematic-aov-cancelling-chart"
    }
}

struct PanickingChart;

impl Chart for PanickingChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        panic!("declared cinematic AOV tile panic")
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn trace_value_enclosure(
        &self,
        _point: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::exact(0.0)
    }

    fn name(&self) -> &'static str {
        "cinematic-aov-panicking-chart"
    }
}

fn shutter() -> ShutterInterval {
    ShutterInterval::resolve(
        0.5,
        0.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 1_024 },
        ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
    )
    .unwrap()
}

fn config_with_times(
    profile: CinematicAovProfile,
    previous_frame_time_s: f64,
    frame_time_s: f64,
    next_frame_time_s: f64,
) -> CinematicAovConfig {
    CinematicAovConfig::new(
        profile,
        CinematicAovProvenance::try_new(
            12,
            frame_time_s,
            previous_frame_time_s,
            next_frame_time_s,
            hash("trajectory"),
            hash("scene"),
            hash("composition"),
        )
        .unwrap(),
        CinematicAovLimits::default(),
    )
}

fn config(profile: CinematicAovProfile) -> CinematicAovConfig {
    config_with_times(profile, 0.0, 0.5, 1.0)
}

fn config_with_limits(
    profile: CinematicAovProfile,
    limits: CinematicAovLimits,
) -> CinematicAovConfig {
    CinematicAovConfig::new(
        profile,
        CinematicAovProvenance::try_new(
            12,
            0.5,
            0.0,
            1.0,
            hash("trajectory"),
            hash("scene"),
            hash("composition"),
        )
        .unwrap(),
        limits,
    )
}

fn assert_film_bits_eq(left: &Film, right: &Film) {
    assert_eq!((left.width, left.height), (right.width, right.height));
    assert_eq!(left.spp_done, right.spp_done);
    assert_eq!(left.time_mode, right.time_mode);
    for (left, right) in left.xyz.iter().zip(&right.xyz) {
        assert_eq!(left.map(f64::to_bits), right.map(f64::to_bits));
    }
}

fn channel<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a fs_img::Channel {
    decoded
        .channels
        .iter()
        .find(|channel| channel.name == name)
        .unwrap_or_else(|| panic!("missing EXR channel {name}"))
}

fn attribute<'a>(decoded: &'a fs_img::DecodedExr, name: &str) -> &'a fs_img::ExrAttribute {
    decoded
        .attributes
        .iter()
        .find(|attribute| attribute.name == name)
        .unwrap_or_else(|| panic!("missing EXR attribute {name}"))
}

fn closed_slab_mesh(z_front: f64, z_back: f64) -> TriMesh {
    let half = 4.0;
    TriMesh::new(
        vec![
            [-half, -half, z_back],
            [half, -half, z_back],
            [half, half, z_back],
            [-half, half, z_back],
            [-half, -half, z_front],
            [half, -half, z_front],
            [half, half, z_front],
            [-half, half, z_front],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    )
}

fn dielectric_scene() -> Scene {
    let white = lift_rgb([1.0, 1.0, 1.0]);
    let emission = (white, 3.0);
    let glass = Material::Dielectric {
        glass: DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    };
    let emitter = TriMesh::new(
        vec![
            [-2.0, -2.0, -1.0],
            [2.0, -2.0, -1.0],
            [2.0, 2.0, -1.0],
            [-2.0, 2.0, -1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    );
    Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Mesh(closed_slab_mesh(0.0, -0.2)),
                material: glass,
                emission: None,
            },
            Primitive {
                shape: Shape::Mesh(emitter),
                material: Material::Lambertian { reflectance: white },
                emission: Some(emission),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(-2.0, -2.0, -1.0),
            edge_u: Vec3::new(4.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 4.0, 0.0),
            prim: 1,
            emission,
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.1,
        },
    }
}

#[test]
fn g5_aov_capture_preserves_cinematic_beauty_bits() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(3);
    let (legacy, with_aovs) = with_cx(|cx| {
        (
            render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter()).unwrap(),
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter(),
                config(CinematicAovProfile::FinalDiagnostic),
            )
            .unwrap(),
        )
    });
    assert_film_bits_eq(&legacy, with_aovs.beauty());
    assert_eq!(with_aovs.sample_count(0), Some(3));
    assert_eq!(with_aovs.primary_count(0), Some(3));
}

#[test]
fn g5_parallel_aovs_and_exr_are_bit_exact_across_workers_tiles_and_schedules() {
    let scene = scene(true);
    let camera = camera(&scene);
    for (profile, sampler) in [
        (CinematicAovProfile::DailyCore, Sampler::Iid),
        (CinematicAovProfile::FinalDiagnostic, Sampler::OwenSobol),
    ] {
        let mut settings = settings(3);
        settings.width = 9;
        settings.height = 7;
        settings.sampler = sampler;
        let aov_config = config(profile);
        let serial = with_cx(|cx| {
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter(),
                aov_config,
            )
            .expect("serial cinematic AOV oracle")
        });
        let serial_exr = serial.to_exr().expect("serial AOV EXR");
        for workers in [1, 2, 8] {
            let mut policy = execution(4, 3, workers, 64 << 20, 0xa051_1000 + workers as u64);
            if workers == 8 {
                policy = policy
                    .with_quantum_weights(vec![1, 3, 2, 5, 4, 7, 6, 8])
                    .expect("valid deliberately skewed worker schedule");
            }
            let parallel = with_cx(|cx| {
                render_cinematic_with_aovs_execution(
                    &scene,
                    &camera,
                    CutSide::After,
                    cx,
                    &settings,
                    shutter(),
                    aov_config,
                    &policy,
                )
                .expect("tile-parallel cinematic AOV render")
            });
            assert_eq!(
                parallel.film, serial,
                "complete AOV state drifted: profile={profile:?} workers={workers}"
            );
            assert_eq!(
                parallel.film.to_exr().unwrap(),
                serial_exr,
                "AOV EXR bytes drifted: profile={profile:?} workers={workers}"
            );
            assert_eq!(parallel.report.requested_workers, workers);
            assert_eq!(parallel.report.executor.declared_run, policy.run_id());
            assert_eq!(parallel.report.memory.used_bytes, 0);
        }
    }
}

#[test]
fn g5_parked_crew_reuses_workers_for_complete_aov_frames() {
    let scene = scene(true);
    let camera = camera(&scene);
    let mut settings = settings(2);
    settings.width = 8;
    settings.height = 6;
    let daily = config(CinematicAovProfile::DailyCore);
    let diagnostic = config(CinematicAovProfile::FinalDiagnostic);
    let first_policy = execution(3, 2, 4, 64 << 20, 0xa051_2001);
    let second_policy = execution(5, 4, 4, 64 << 20, 0xa051_2002);
    let first_serial = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            daily,
        )
        .unwrap()
    });
    let second_serial = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            diagnostic,
        )
        .unwrap()
    });
    with_cx(|cx| {
        let pool = RenderWorkerPool::new(&first_policy, cx.mode(), 0xa051_2ced);
        pool.with_parked_crew_local(|parked| {
            let first = parked
                .render_cinematic_with_aovs(
                    &scene,
                    &camera,
                    CutSide::After,
                    cx,
                    &settings,
                    shutter(),
                    daily,
                    &first_policy,
                )
                .unwrap();
            let second = parked
                .render_cinematic_with_aovs(
                    &scene,
                    &camera,
                    CutSide::After,
                    cx,
                    &settings,
                    shutter(),
                    diagnostic,
                    &second_policy,
                )
                .unwrap();
            assert_eq!(first.film, first_serial);
            assert_eq!(second.film, second_serial);
            assert_eq!(first.report.executor.declared_run, first_policy.run_id());
            assert_eq!(second.report.executor.declared_run, second_policy.run_id());
        });
    });
}

#[test]
fn g4_parallel_aov_refusal_cancellation_and_panic_publish_no_film() {
    let healthy_scene = scene(true);
    let healthy_camera = camera(&healthy_scene);
    let mut render_settings = settings(4);
    render_settings.width = 8;
    render_settings.height = 6;
    let aov_config = config(CinematicAovProfile::DailyCore);
    let one_byte = execution(3, 2, 4, 1, 0xa051_3001);
    let refused = with_cx(|cx| {
        render_cinematic_with_aovs_execution(
            &healthy_scene,
            &healthy_camera,
            CutSide::After,
            cx,
            &render_settings,
            shutter(),
            aov_config,
            &one_byte,
        )
    });
    assert!(matches!(
        refused,
        Err(CinematicAovExecutionError::Execution(
            RenderExecutionError::Memory(_)
        ))
    ));

    let bad_times = config_with_times(CinematicAovProfile::DailyCore, 0.0, 0.75, 1.0);
    let normal_policy = execution(3, 2, 4, 64 << 20, 0xa051_3002);
    let bad_config = with_cx(|cx| {
        render_cinematic_with_aovs_execution(
            &healthy_scene,
            &healthy_camera,
            CutSide::After,
            cx,
            &render_settings,
            shutter(),
            bad_times,
            &normal_policy,
        )
    });
    assert_eq!(
        bad_config,
        Err(CinematicAovExecutionError::Aov(
            CinematicAovError::ReferenceTimesDoNotCoverShutter
        ))
    );

    let gate = Arc::new(CancelGate::new());
    let evaluations = Arc::new(AtomicUsize::new(0));
    let mut cancelling_scene = scene(false);
    cancelling_scene.primitives[0].shape = Shape::Chart(Box::new(CancellingChart {
        gate: Arc::clone(&gate),
        evaluations: Arc::clone(&evaluations),
        cancel_at: 4,
    }));
    let cancelling_camera = camera(&cancelling_scene);
    let cancelled = with_gate(&gate, |cx| {
        render_cinematic_with_aovs_execution(
            &cancelling_scene,
            &cancelling_camera,
            CutSide::After,
            cx,
            &render_settings,
            shutter(),
            aov_config,
            &normal_policy,
        )
    });
    assert!(evaluations.load(Ordering::SeqCst) >= 4);
    assert!(matches!(
        cancelled,
        Err(CinematicAovExecutionError::Execution(
            RenderExecutionError::Tracer(_)
        ))
    ));

    let mut panicking_scene = scene(false);
    panicking_scene.primitives[0].shape = Shape::Chart(Box::new(PanickingChart));
    let panicking_camera = camera(&panicking_scene);
    let panicked = with_cx(|cx| {
        render_cinematic_with_aovs_execution(
            &panicking_scene,
            &panicking_camera,
            CutSide::After,
            cx,
            &render_settings,
            shutter(),
            aov_config,
            &normal_policy,
        )
    });
    assert!(matches!(
        panicked,
        Err(CinematicAovExecutionError::Execution(
            RenderExecutionError::Executor(RunError::TilePanicked { .. })
        ))
    ));
}

#[test]
fn g0_owned_denoise_guides_are_exact_exr_planes_with_background_zeros() {
    let foreground_scene = scene(true);
    let foreground_camera = camera(&foreground_scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &foreground_scene,
            &foreground_camera,
            CutSide::After,
            cx,
            &settings(3),
            shutter(),
            config(CinematicAovProfile::DailyCore),
        )
        .unwrap()
    });
    let guides = film.denoise_guides().unwrap();
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert_eq!((guides.width(), guides.height()), (1, 1));
    for (guide, exr_name) in [
        (guides.motion_prev_x(), "motion.prev.X"),
        (guides.motion_prev_y(), "motion.prev.Y"),
        (guides.axial_depth_m(), "depth.Z"),
        (guides.normal_x(), "normal.X"),
        (guides.normal_y(), "normal.Y"),
        (guides.normal_z(), "normal.Z"),
        (guides.primary_coverage(), "primary.coverage"),
        (guides.variance_luminance(), "variance.Y"),
    ] {
        assert_eq!(
            guide,
            channel(&decoded, exr_name).data.as_slice(),
            "{exr_name}"
        );
    }

    let mut background = scene(false);
    background.primitives.remove(0);
    background.lights[0].prim = 0;
    let background_camera = camera(&background);
    let background_film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &background,
            &background_camera,
            CutSide::After,
            cx,
            &settings(2),
            shutter(),
            config(CinematicAovProfile::DailyCore),
        )
        .unwrap()
    });
    let background_guides = background_film.denoise_guides().unwrap();
    for plane in [
        background_guides.motion_prev_x(),
        background_guides.motion_prev_y(),
        background_guides.axial_depth_m(),
        background_guides.normal_x(),
        background_guides.normal_y(),
        background_guides.normal_z(),
        background_guides.primary_coverage(),
        background_guides.variance_luminance(),
    ] {
        assert!(plane.iter().all(|value| value.to_bits() == 0));
    }

    let beauty_only = with_cx(|cx| {
        render_cinematic_with_aovs(
            &foreground_scene,
            &foreground_camera,
            CutSide::After,
            cx,
            &settings(1),
            shutter(),
            config(CinematicAovProfile::BeautyOnly),
        )
        .unwrap()
    });
    assert_eq!(
        beauty_only.denoise_guides(),
        Err(CinematicAovError::DenoiseGuidesUnavailable)
    );
}

#[test]
fn g0_final_exr_round_trip_has_exact_channels_palettes_and_units() {
    let scene = scene(true);
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(2),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let first = film.to_exr().unwrap();
    let second = film.to_exr().unwrap();
    assert_eq!(first, second, "EXR encoding must be byte deterministic");
    let decoded = fs_img::read_exr(&first).unwrap();
    assert_eq!(decoded.channels.len(), 30);
    assert!(
        decoded
            .channels
            .windows(2)
            .all(|pair| pair[0].name < pair[1].name)
    );
    assert!(decoded.channels.iter().all(|plane| {
        plane.ty == fs_img::PixelType::Float && plane.data.iter().all(|value| value.is_finite())
    }));
    assert!((channel(&decoded, "depth.Z").data[0] - 2.0).abs() < 1.0e-6);
    assert_eq!(channel(&decoded, "id.object").data, [1.0]);
    assert_eq!(channel(&decoded, "id.material").data, [1.0]);
    assert_eq!(channel(&decoded, "samples").data, [2.0]);
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_eq!(
        mask & (validity::PRIMARY
            | validity::ALBEDO
            | validity::PREVIOUS_MOTION
            | validity::OBJECT_ID
            | validity::MATERIAL_ID
            | validity::CONTRIBUTION_SPLIT),
        validity::PRIMARY
            | validity::ALBEDO
            | validity::PREVIOUS_MOTION
            | validity::OBJECT_ID
            | validity::MATERIAL_ID
            | validity::CONTRIBUTION_SPLIT
    );
    assert!(channel(&decoded, "motion.prev.X").data[0] < 0.0);
    assert_eq!(
        attribute(&decoded, "frankensim.aov.authority").value,
        b"raw-estimate"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.render.shutter").value,
        b"convention=front-loaded;distribution=stratified-counter-v1;strata=1024"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.aov.objectPalette").value,
        b"0=unavailable;1=101"
    );
    let material_palette =
        std::str::from_utf8(&attribute(&decoded, "frankensim.aov.materialPalette").value).unwrap();
    assert!(material_palette.starts_with("0=unavailable;1="));
}

#[test]
fn g0_every_profile_has_its_exact_frozen_channel_schema() {
    let scene = scene(true);
    let camera = camera(&scene);
    let cases: &[(CinematicAovProfile, &[&str])] = &[
        (CinematicAovProfile::BeautyOnly, &["B", "G", "R"]),
        (
            CinematicAovProfile::DailyCore,
            &[
                "B",
                "G",
                "R",
                "albedo.B",
                "albedo.G",
                "albedo.R",
                "depth.Z",
                "motion.prev.X",
                "motion.prev.Y",
                "normal.X",
                "normal.Y",
                "normal.Z",
                "primary.coverage",
                "variance.Y",
            ],
        ),
        (
            CinematicAovProfile::FinalDiagnostic,
            &[
                "B",
                "G",
                "R",
                "albedo.B",
                "albedo.G",
                "albedo.R",
                "depth.Z",
                "diagnostic.validity",
                "direct.B",
                "direct.G",
                "direct.R",
                "emission.B",
                "emission.G",
                "emission.R",
                "id.material",
                "id.object",
                "indirect.B",
                "indirect.G",
                "indirect.R",
                "motion.prev.X",
                "motion.prev.Y",
                "normal.X",
                "normal.Y",
                "normal.Z",
                "normal_geom.X",
                "normal_geom.Y",
                "normal_geom.Z",
                "primary.coverage",
                "samples",
                "variance.Y",
            ],
        ),
    ];

    for &(profile, expected_names) in cases {
        let film = with_cx(|cx| {
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings(1),
                shutter(),
                config(profile),
            )
            .unwrap()
        });
        let first = film.to_exr().unwrap();
        assert_eq!(
            first,
            film.to_exr().unwrap(),
            "{profile:?} EXR bytes drifted"
        );
        let decoded = fs_img::read_exr(&first).unwrap();
        let actual_names = decoded
            .channels
            .iter()
            .map(|channel| channel.name.as_str())
            .collect::<Vec<_>>();
        assert_eq!(actual_names, expected_names, "{profile:?} schema drifted");
        assert!(
            decoded
                .channels
                .iter()
                .all(|channel| channel.ty == fs_img::PixelType::Float)
        );
    }
}

#[test]
fn g0_true_background_exports_finite_zero_surface_aovs_and_zero_ids() {
    let mut background = scene(false);
    // Retain the admitted rectangular emitter at x=9 while removing the
    // camera-facing primitive. The primary ray therefore observes a genuine
    // black background without violating the tracer's finite-emitter
    // admission contract.
    background.primitives.remove(0);
    background.lights[0].prim = 0;
    let camera = camera(&background);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &background,
            &camera,
            CutSide::After,
            cx,
            &settings(3),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    assert_eq!(film.sample_count(0), Some(3));
    assert_eq!(film.primary_count(0), Some(0));

    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    for plane in &decoded.channels {
        assert!(plane.data.iter().all(|value| value.is_finite()));
        match plane.name.as_str() {
            "samples" => assert_eq!(plane.data, [3.0]),
            "diagnostic.validity" => {
                assert_eq!(plane.data, [validity::CONTRIBUTION_SPLIT as f32]);
            }
            _ => assert_eq!(plane.data, [0.0], "background {} was nonzero", plane.name),
        }
    }
    assert_eq!(channel(&decoded, "primary.coverage").data, [0.0]);
    assert_eq!(channel(&decoded, "id.object").data, [0.0]);
    assert_eq!(channel(&decoded, "id.material").data, [0.0]);
    assert_eq!(
        attribute(&decoded, "frankensim.aov.objectPalette").value,
        b"0=unavailable"
    );
    assert!(
        attribute(&decoded, "frankensim.aov.materialPalette")
            .value
            .starts_with(b"0=unavailable;1="),
        "the off-camera emitter remains in the exact scene palette"
    );
}

#[test]
fn g0_foreground_background_edge_keeps_every_aov_on_its_beauty_pixel() {
    let mut scene = scene(true);
    let stationary = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(0.0, RigidTransform::identity(), [0.0, 0.0, 0.0]).unwrap(),
        TransformKeyframe::try_new(1.0, RigidTransform::identity(), [0.0, 0.0, 0.0]).unwrap(),
    ])
    .unwrap();
    scene.primitives[0].shape = Shape::AnimatedInstance(
        AnimatedGeometryInstance::try_new(
            101,
            hash("left-pixel-geometry"),
            SharedGeometry::mesh(left_pixel_quad()),
            stationary,
        )
        .unwrap(),
    );
    let camera = camera(&scene);
    let mut edge_settings = settings(4);
    edge_settings.width = 2;
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &edge_settings,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    assert_eq!(film.beauty().width, 2);
    assert_eq!(film.beauty().height, 1);
    assert_eq!(film.primary_count(0), Some(4));
    assert_eq!(film.primary_count(1), Some(0));

    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert!(decoded.channels.iter().all(|plane| {
        plane.data.len() == 2 && plane.data.iter().all(|sample| sample.is_finite())
    }));
    assert_eq!(channel(&decoded, "samples").data, [4.0, 4.0]);
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0, 0.0]);
    assert_eq!(channel(&decoded, "id.object").data, [1.0, 0.0]);
    assert_eq!(channel(&decoded, "id.material").data, [1.0, 0.0]);
    assert!(channel(&decoded, "depth.Z").data[0] > 0.0);
    assert_eq!(channel(&decoded, "depth.Z").data[1], 0.0);
    for name in [
        "albedo.R",
        "albedo.G",
        "albedo.B",
        "normal.Z",
        "normal_geom.Z",
    ] {
        let plane = &channel(&decoded, name).data;
        assert_ne!(plane[0], 0.0, "foreground {name} must be populated");
        assert_eq!(plane[1], 0.0, "background {name} leaked across the edge");
    }
    for name in ["R", "G", "B", "emission.R", "emission.G", "emission.B"] {
        assert_eq!(
            channel(&decoded, name).data[1],
            0.0,
            "background {name} leaked across the edge"
        );
    }
    assert!(
        ["R", "G", "B"]
            .into_iter()
            .any(|name| channel(&decoded, name).data[0] != 0.0),
        "foreground beauty must be present at pixel 0"
    );
    let foreground_validity = channel(&decoded, "diagnostic.validity").data[0] as u32;
    let background_validity = channel(&decoded, "diagnostic.validity").data[1] as u32;
    assert_ne!(foreground_validity & validity::PRIMARY, 0);
    assert_ne!(foreground_validity & validity::OBJECT_ID, 0);
    assert_ne!(foreground_validity & validity::MATERIAL_ID, 0);
    assert_eq!(background_validity, validity::CONTRIBUTION_SPLIT);
}

#[test]
fn g0_camera_visible_emission_split_reconstructs_each_beauty_sample_and_exr() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(4);
    with_cx(|cx| {
        for sample_index in 0..settings.spp {
            let sample = trace_cinematic_pixel_sample(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter(),
                0,
                sample_index,
            )
            .unwrap();
            assert!(sample.primary.is_some());
            assert!(sample.emission_xyz.iter().any(|value| *value > 0.0));
            assert_eq!(sample.direct_xyz, [0.0; 3]);
            assert_eq!(sample.indirect_xyz, [0.0; 3]);
            for component in 0..3 {
                let reconstructed = sample.direct_xyz[component]
                    + sample.indirect_xyz[component]
                    + sample.emission_xyz[component];
                let tolerance = 16.0
                    * f64::EPSILON
                    * sample.xyz[component]
                        .abs()
                        .max(reconstructed.abs())
                        .max(1.0);
                assert!(
                    (sample.xyz[component] - reconstructed).abs() <= tolerance,
                    "sample {sample_index} XYZ component {component} failed reconstruction"
                );
            }
        }
    });

    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    for suffix in ["R", "G", "B"] {
        let beauty = channel(&decoded, suffix).data[0];
        let reconstructed = channel(&decoded, &format!("direct.{suffix}")).data[0]
            + channel(&decoded, &format!("indirect.{suffix}")).data[0]
            + channel(&decoded, &format!("emission.{suffix}")).data[0];
        let tolerance = 4.0 * f32::EPSILON * beauty.abs().max(reconstructed.abs()).max(1.0);
        assert!((beauty - reconstructed).abs() <= tolerance);
        assert_eq!(channel(&decoded, &format!("direct.{suffix}")).data, [0.0]);
        assert_eq!(channel(&decoded, &format!("indirect.{suffix}")).data, [0.0]);
    }
}

#[test]
fn g3_adaptive_aovs_preserve_the_exact_uniform_prefix_and_export_terminal_counts() {
    let scene = scene(true);
    let camera = camera(&scene);
    let maximum = settings(4);
    let policy = AdaptiveSamplingConfig::try_new(2, 1, 1.0e30, 0.0, 0.0).unwrap();
    let adaptive = with_cx(|cx| {
        render_cinematic_adaptive_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &maximum,
            policy,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    assert_eq!(adaptive.beauty().sample_counts(), [2]);
    assert_eq!(
        adaptive.beauty().decisions(),
        [AdaptiveDecision::ErrorThreshold]
    );

    let prefix_settings = settings(2);
    let uniform_prefix = with_cx(|cx| {
        render_cinematic(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &prefix_settings,
            shutter(),
        )
        .unwrap()
    });
    for (adaptive, uniform) in adaptive.beauty().xyz_sums()[0]
        .iter()
        .zip(uniform_prefix.xyz[0])
    {
        assert_eq!(adaptive.to_bits(), uniform.to_bits());
    }

    let first = adaptive.to_exr().unwrap();
    assert_eq!(first, adaptive.to_exr().unwrap());
    let decoded = fs_img::read_exr(&first).unwrap();
    assert_eq!(channel(&decoded, "samples").data, [2.0]);
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    assert_eq!(
        attribute(&decoded, "frankensim.render.sampleMode").value,
        b"adaptive"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.render.spp").value,
        b"per-pixel-channel"
    );
    assert_eq!(
        attribute(&decoded, "frankensim.render.sppCeiling").value,
        b"4"
    );
    let adaptive_policy =
        std::str::from_utf8(&attribute(&decoded, "frankensim.render.adaptive").value).unwrap();
    assert!(adaptive_policy.starts_with("version=1;minimum=2;batch=1;"));

    let daily = with_cx(|cx| {
        render_cinematic_adaptive_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &maximum,
            policy,
            shutter(),
            config(CinematicAovProfile::DailyCore),
        )
        .unwrap()
    });
    let daily = fs_img::read_exr(&daily.to_exr().unwrap()).unwrap();
    assert!(
        daily
            .channels
            .iter()
            .all(|channel| channel.name != "samples")
    );
    assert_eq!(
        attribute(&daily, "frankensim.render.spp").value,
        b"per-pixel-unexported-by-profile"
    );
}

#[test]
fn g0_shading_normal_matches_the_face_forwarded_beauty_frame() {
    let mut scene = scene(false);
    scene.primitives[0].shape = Shape::Mesh(backface_quad());
    scene.primitives[0].emission = None;
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(1),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert!(channel(&decoded, "normal.Z").data[0] > 0.999);
    assert!(channel(&decoded, "normal_geom.Z").data[0] < -0.999);
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_eq!(mask & validity::AUTHORED_SHADING_NORMAL, 0);
}

#[test]
fn g0_dielectric_first_hit_retains_surface_diagnostics_without_forging_albedo() {
    let scene = dielectric_scene();
    let camera = camera(&scene);
    let mut settings = settings(2);
    settings.max_depth = 3;
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert_eq!(channel(&decoded, "primary.coverage").data, [1.0]);
    assert!((channel(&decoded, "depth.Z").data[0] - 1.0).abs() <= 1.0e-6);
    for name in ["albedo.R", "albedo.G", "albedo.B"] {
        assert_eq!(channel(&decoded, name).data, [0.0]);
    }
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_ne!(mask & validity::PRIMARY, 0);
    assert_ne!(mask & validity::MATERIAL_ID, 0);
    assert_eq!(mask & validity::ALBEDO, 0);
    assert!(channel(&decoded, "id.material").data[0] > 0.0);

    for prefix in ["normal", "normal_geom"] {
        let squared_norm = ["X", "Y", "Z"]
            .into_iter()
            .map(|suffix| {
                let value = channel(&decoded, &format!("{prefix}.{suffix}")).data[0];
                value * value
            })
            .sum::<f32>();
        assert!((squared_norm - 1.0).abs() <= 8.0 * f32::EPSILON);
    }
    assert!(
        ["indirect.R", "indirect.G", "indirect.B"]
            .into_iter()
            .any(|name| channel(&decoded, name).data[0] > 0.0),
        "the transmitted emitter should be classified as a multi-bounce contribution"
    );
    for suffix in ["R", "G", "B"] {
        let beauty = channel(&decoded, suffix).data[0];
        let reconstructed = channel(&decoded, &format!("direct.{suffix}")).data[0]
            + channel(&decoded, &format!("indirect.{suffix}")).data[0]
            + channel(&decoded, &format!("emission.{suffix}")).data[0];
        let tolerance = 4.0 * f32::EPSILON * beauty.abs().max(reconstructed.abs()).max(1.0);
        assert!((beauty - reconstructed).abs() <= tolerance);
    }
}

#[test]
fn g3_progressive_ranges_equal_straight_through_and_mismatch_rolls_back() {
    let mut scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(4);
    let straight = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let mut progressive = fs_render::aov::CinematicAovFilm::try_new(
        1,
        1,
        config(CinematicAovProfile::FinalDiagnostic),
    )
    .unwrap();
    with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            0,
            2,
            shutter(),
        )
        .unwrap();
    });

    let before = progressive.to_exr().unwrap();
    let changed_camera = camera_at(&scene, Point3::new(0.125, 0.0, 2.0));
    let camera_error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &changed_camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            2,
            3,
            shutter(),
        )
    });
    assert_eq!(
        camera_error,
        Err(CinematicAovError::ProgressiveBindingMismatch)
    );
    assert_eq!(progressive.to_exr().unwrap(), before);

    scene.primitives[0].material = Material::Lambertian {
        reflectance: lift_rgb([0.1, 0.8, 0.3]),
    };
    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            2,
            3,
            shutter(),
        )
    });
    assert_eq!(error, Err(CinematicAovError::ProgressiveBindingMismatch));
    assert_eq!(progressive.to_exr().unwrap(), before);

    scene.primitives[0].material = Material::Lambertian {
        reflectance: lift_rgb([0.7, 0.4, 0.2]),
    };
    with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            2,
            4,
            shutter(),
        )
        .unwrap();
    });
    assert_film_bits_eq(straight.beauty(), progressive.beauty());
    assert_eq!(straight.to_exr().unwrap(), progressive.to_exr().unwrap());
}

#[test]
fn g4_reference_time_refusal_precedes_any_aov_publication() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(2);
    let bad_config = config_with_times(CinematicAovProfile::FinalDiagnostic, 0.0, 0.75, 1.0);
    let mut film = fs_render::aov::CinematicAovFilm::try_new(1, 1, bad_config).unwrap();
    let error = with_cx(|cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut film,
            0,
            2,
            shutter(),
        )
    });
    assert_eq!(
        error,
        Err(CinematicAovError::ReferenceTimesDoNotCoverShutter)
    );
    assert_eq!(film.beauty().spp_done, 0);
    assert_eq!(film.sample_count(0), Some(0));
    assert_eq!(film.primary_count(0), Some(0));
    assert!(film.palette().is_none());
    assert_eq!(film.to_exr(), Err(CinematicAovError::UnboundFilm));
}

#[test]
fn g4_pre_cancelled_render_leaves_every_public_plane_uncommitted() {
    let scene = scene(true);
    let camera = camera(&scene);
    let settings = settings(2);
    let mut film = fs_render::aov::CinematicAovFilm::try_new(
        1,
        1,
        config(CinematicAovProfile::FinalDiagnostic),
    )
    .unwrap();
    let gate = CancelGate::new();
    gate.request();
    let error = with_gate(&gate, |cx| {
        render_cinematic_range_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut film,
            0,
            2,
            shutter(),
        )
    });
    assert!(matches!(error, Err(CinematicAovError::Tracer(_))));
    assert_eq!(film.beauty().spp_done, 0);
    assert_eq!(film.sample_count(0), Some(0));
    assert_eq!(film.primary_count(0), Some(0));
    assert!(film.palette().is_none());
    assert_eq!(film.to_exr(), Err(CinematicAovError::UnboundFilm));
}

#[test]
fn g0_direct_mesh_object_identity_and_motion_are_unavailable_not_forged() {
    let scene = scene(false);
    let camera = camera(&scene);
    let film = with_cx(|cx| {
        render_cinematic_with_aovs(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings(1),
            shutter(),
            config(CinematicAovProfile::FinalDiagnostic),
        )
        .unwrap()
    });
    let decoded = fs_img::read_exr(&film.to_exr().unwrap()).unwrap();
    assert_eq!(channel(&decoded, "id.object").data, [0.0]);
    assert_eq!(channel(&decoded, "motion.prev.X").data, [0.0]);
    let mask = channel(&decoded, "diagnostic.validity").data[0] as u32;
    assert_eq!(mask & validity::OBJECT_ID, 0);
    assert_eq!(mask & validity::PREVIOUS_MOTION, 0);
    assert_ne!(mask & validity::PRIMARY, 0);
    assert_ne!(mask & validity::MATERIAL_ID, 0);
}

#[test]
fn g0_each_export_buffer_category_has_an_independent_fail_closed_ceiling() {
    let scene = scene(true);
    let camera = camera(&scene);
    let defaults = CinematicAovLimits::default();
    let limits = |planes, metadata, scratch, output| {
        CinematicAovLimits::try_new(
            defaults.max_pixels(),
            defaults.max_retained_bytes(),
            planes,
            metadata,
            scratch,
            output,
            defaults.max_palette_entries(),
        )
        .unwrap()
    };
    let render = |limits| {
        with_cx(|cx| {
            render_cinematic_with_aovs(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings(1),
                shutter(),
                config_with_limits(CinematicAovProfile::FinalDiagnostic, limits),
            )
            .unwrap()
        })
    };

    assert!(matches!(
        render(limits(
            1,
            defaults.max_export_metadata_bytes(),
            defaults.max_exr_encoder_scratch_bytes(),
            defaults.max_encoded_exr_bytes(),
        ))
        .to_exr(),
        Err(CinematicAovError::ExportMemoryLimit { .. })
    ));
    assert!(matches!(
        render(limits(
            defaults.max_export_plane_bytes(),
            1,
            defaults.max_exr_encoder_scratch_bytes(),
            defaults.max_encoded_exr_bytes(),
        ))
        .to_exr(),
        Err(CinematicAovError::ExportMetadataMemoryLimit { .. })
    ));
    assert!(matches!(
        render(limits(
            defaults.max_export_plane_bytes(),
            defaults.max_export_metadata_bytes(),
            1,
            defaults.max_encoded_exr_bytes(),
        ))
        .to_exr(),
        Err(CinematicAovError::ExrEncoderScratchLimit { .. })
    ));
    assert!(matches!(
        render(limits(
            defaults.max_export_plane_bytes(),
            defaults.max_export_metadata_bytes(),
            defaults.max_exr_encoder_scratch_bytes(),
            1,
        ))
        .to_exr(),
        Err(CinematicAovError::EncodedExrMemoryLimit { .. })
    ));
}
